/** @file gsC1BasisFunction.cpp

    @brief Searches for C1 Basisfunctions on a multipatch domain

    This file is part of the G+Smo library and based on the file
    gsGcontFunctions.cpp from K. Birner and A. Mantazafaris.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#include <gismo.h>
#include <gismo_dev.h> // Damit auf den unsupported Ordner zugegriffen werden kann!

using namespace gismo;

// +++++++++++++++++++++++++++++++++ HELPER-FUNCTION +++++++++++++++++++++++++++++++++++++

void alpha_S(const gsMultiPatch<> mp,
                  const index_t side,
                  const index_t degree,
                  gsGeometry<>::Ptr & alpha); // return alpha

void beta(const gsMultiPatch<> mp,
          const index_t degree,
                  gsGeometry<>::Ptr & beta); // return beta

void beta_L_R(const gsGeometry<>::Ptr alpha_L,
              const gsGeometry<>::Ptr alpha_R,
              const gsGeometry<>::Ptr beta_S,
              const index_t deg_beta_L,
              const index_t deg_beta_R,
              gsGeometry<>::uPtr & beta_L,
              gsGeometry<>::uPtr & beta_R); // return beta_L and beta_R

void gamma(const gsMatrix<> uv,
           real_t & gamma)
{
    gsFunctionExpr<> f("1000/(8167 + 60*v - 407*v*v + 516*v*v*v + 407*v*v*v*v)",1);
    real_t v = uv.at(1);
    gamma = 10000/(8167 + 60*v - 407*v*v + 516*v*v*v + 407*v*v*v*v);
    gamma = 1;
}

void g_S(const gsMultiPatch<> mp,
         const gsMultiPatch<> mp_tilde,
         const gsBSpline<> b_L,
         const gsBSpline<> b_R,
         const gsBSpline<> a_L,
         const gsBSpline<> a_R,
         const index_t side,
         gsTensorBSpline<2>  & g,
         gsTensorBSpline<2>  & g1);



// ++++++++++++++++++++++++++++++++++++++++ MAIN +++++++++++++++++++++++++++++++++++++++++
int main(int argc, char *argv[])
{
    bool plot = true;

    // ========= INPUTS =========
    std::string inputFile("inputData/GCFunction.xml"); // Geometry, two 2D-Patches
    gsFileData<> fileData(inputFile); // Store the inputFile in fileData

    gsOptionList optionList; // Option List (Degree, Refinements, etc...)
    fileData.getId(0,optionList);

    gsInfo << "Input parameters are: \n" << optionList << "\n";
    // ==========================

    // ======== GEOMETRY ========
    gsMultiPatch<> multiPatch;
    fileData.getId(110,multiPatch); // Multipatch id=100 or 110 or 120
    multiPatch.computeTopology();

    // 0 == only one interface exists
    gsInfo << "Got Multipatch:\n" << multiPatch << "\n";
    gsInfo << "Got Interface:\n" << multiPatch.bInterface(0) << "\n\n";
    // ==========================


    // ======== BASIS ========
    // Refinement
    index_t numRefine = optionList.getInt("numRefine"); // number of knot insertion
    // number of the max. polynom degree
    index_t polynomDegree = optionList.getInt("polynomDegree");

    // Get the degree of the input // Only for the first patch TODO
    index_t degreeInput = multiPatch.patch(0).basis().maxDegree();

    // Regularity within patch // Only regularity < polynomDegree !
    index_t regularity = optionList.getInt("regularity");
    // Refinements
    multiPatch.degreeElevate(polynomDegree-degreeInput);
    multiPatch.uniformRefine(numRefine,polynomDegree - regularity);

    gsInfo << "Basis functions (of the first patch): \n" << multiPatch.patch(0).basis() << "\n\n";
    // =======================

    gsInfo << "Computing alpha_S and beta... \n\n";
    // Compute alpha_S and beta
    const index_t left = 1; // Patch 111
    const index_t right = 0; // Patch 112

    const index_t degree_alpha_left = optionList.getInt("deg_alpha_left");
    const index_t degree_alpha_right = optionList.getInt("deg_alpha_right");
    const index_t degree_beta = optionList.getInt("deg_beta");

    gsGeometry<>::Ptr alpha_left, alpha_right, beta_S;

    alpha_S(multiPatch,left,degree_alpha_left,alpha_left); // alpha_left
    alpha_S(multiPatch,right,degree_alpha_right,alpha_right); // alpha_right

    beta(multiPatch,degree_beta,beta_S); // beta

    gsInfo << "alpha: \n" << alpha_left->coefs();

    // =======================
    // COMPUTE BETA_L and BETA_R
    gsGeometry<>::uPtr beta_L, beta_R;

    const index_t deg_beta_L = optionList.getInt("deg_beta_left");
    const index_t deg_beta_R = optionList.getInt("deg_beta_right");

    beta_L_R(alpha_left,alpha_right,beta_S,deg_beta_L,deg_beta_R,beta_L,beta_R); // beta_L, beta_R

    gsMatrix<> greville = beta_S->basis().anchors();
    gsMatrix<> randomPoints(1,10);

    randomPoints.setRandom();
    randomPoints = randomPoints.array().abs();

    gsInfo << "Random Points: " << randomPoints << "\n";


    gsInfo << "conditiontest: \n" << alpha_left->eval(randomPoints).transpose().asDiagonal() *  beta_R->eval(randomPoints).transpose()
        - alpha_right->eval(randomPoints).transpose().asDiagonal() * beta_L->eval(randomPoints).transpose()
        - beta_S->eval(randomPoints).transpose() << "\n";

    if (plot)
    {
        // Patch
        gsWriteParaview(multiPatch,"geometry",1000,true);

        // Functionen
        gsWriteParaview(beta_L.operator*(), "beta_L");
        gsWriteParaview(beta_R.operator*(), "beta_R");
        gsWriteParaview(beta_S.operator*(), "beta");

        gsWriteParaview(alpha_right.operator*(), "alpha_R");
    }

    gsMultiPatch<> multiPatch_tilde = multiPatch;

//    gsTensorBSpline<2, real_t> * splineLeftPtr = dynamic_cast<gsTensorBSpline<2, real_t> *> (&multiPatch_tilde.patch(left));
//    gsTensorBSpline<2, real_t> * splineRightPtr = dynamic_cast<gsTensorBSpline<2, real_t> *> (&multiPatch_tilde.patch(right));

    //splineLeftPtr->insertKnot(0.5,1,1);
    //splineLeftPtr->insertKnot(0.7,1,2);

    //splineRightPtr->insertKnot(0.5,1,1);
    //splineRightPtr->insertKnot(0.7,1,1);

    //gsInfo << "iFace: " << multiPatch.basis(left) << "\n";

    // splineLeftPtr->continuityIncrease(1,1); // (int i, int dir)  // TODO: where is this member? @weinmueller
    // splineRightPtr->continuityIncrease(1,1); // (int i, int dir) // TODO: where is this member? @weinmueller

    gsInfo << "iFace: " << multiPatch.basis(left) << "\n";
    gsInfo << "iFace: " << multiPatch.basis(right) << "\n";

    gsWriteParaview(multiPatch,"multiPatch",1000,true);

    gsInfo << "Basis: \n" << multiPatch.patch(left).basis().component(1).function(2).eval(greville) << "\n";

    gsTensorBSpline<2> g_left, g_right,g1_left, g1_right;
    gsBSpline<> b_L,b_R,a_L,a_R, b;

    b_L = dynamic_cast<gsBSpline<> &> (*beta_L);
    b_R = dynamic_cast<gsBSpline<> &> (*beta_R);
    a_L = dynamic_cast<gsBSpline<> &> (*alpha_left);
    a_R = dynamic_cast<gsBSpline<> &> (*alpha_right);
    b = dynamic_cast<gsBSpline<> &> (*beta_S);

    g_S(multiPatch,multiPatch_tilde,b_L,b_R,a_L,a_R,left,g_left,g1_left);
    g_S(multiPatch,multiPatch_tilde,b_L,b_R,a_L,a_R,right,g_right,g1_right);

    gsMatrix<> greville2d(2,randomPoints.size());

    greville2d.setZero();
    greville2d.block(1,0,1,randomPoints.size()) = randomPoints;

    greville2d(0,0) = 0;
    greville2d(1,0) = 0;

    gsInfo << "Greville2d: " << greville2d << "\n";
    gsInfo << "Ableitung: " << g_left.deriv(greville2d) << "\n";
    gsInfo << "Ableitung: " << g_right.deriv(greville2d) << "\n";

    gsInfo << "G1 condition: " << alpha_right->eval(randomPoints).cwiseProduct(g_left.deriv(greville2d).topRows(1))
        - alpha_left->eval(randomPoints).cwiseProduct(g_right.deriv(greville2d).topRows(1))
        + beta_S->eval(randomPoints).cwiseProduct(g_right.deriv(greville2d).bottomRows(1)) << "\n";

    gsFileData<> xml;
    xml << g_right;
    xml << a_R;
    xml << b_L;
    xml << b_R;
    xml << b;
    xml.save("functions");

    gsInfo << "coefs b_L: " << b_L.coefs() << "\n";
    gsInfo << "coefs b_R: " << b_R.coefs() << "\n";
    gsInfo << "coefs b: " << b.coefs() << "\n";

} // END OF MAIN

void g_S(const gsMultiPatch<> mp,
         const gsMultiPatch<> mp_tilde,
         const gsBSpline<> b_L,
         const gsBSpline<> b_R,
         const gsBSpline<> a_L,
         const gsBSpline<> a_R,
         const index_t side,
         gsTensorBSpline<2>  & g,
         gsTensorBSpline<2>  & g1)
{
    const index_t n_tilde = mp_tilde.patch(side).basis().component(1).size();
    gsInfo << "n_tilde: " << n_tilde << "\n";

    const index_t d_alpha = math::max(a_L.degree(0),a_R.degree(0));
    gsInfo << "d_alpha: " << d_alpha << "\n";

    // Basis for g_1^S
    gsKnotVector<> kv (0,1,2,mp.patch(side).degree(0) + 1 - d_alpha,1); // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> bsp_g1_bar(kv); // (u,v)
    const index_t n_bar = bsp_g1_bar.size();
    gsInfo << "n_bar: " << n_bar << "\n";
    gsInfo << "bsp_g1_bar: " << bsp_g1_bar << "\n";

    // Basis for g_0^S
    //gsKnotVector<> kv (0,1,2,degree + 1,2); // first,last,interior,mult_ends,mult_interior,degree
    //gsTensorBSplineBasis<2,real_t> bsp_g(kv,kv); // (u,v)
    gsTensorBSplineBasis<2,real_t> bsp_g0 = dynamic_cast<gsTensorBSplineBasis<2, real_t> & >(mp.patch(side).basis());
    gsTensorBSplineBasis<2,real_t> bsp_g1 = dynamic_cast<gsTensorBSplineBasis<2, real_t> & >(mp.patch(side).basis());



    gsInfo << "Basis Tilde: " << mp_tilde << "\n";

    gsGeometry<>::Ptr g0_temp, g1_temp;

    gsMatrix<> greville = mp.patch(side).basis().anchors(); // (u,v)
    gsMatrix<> ev_L, ev_R, ev1_L, ev1_R, ev;

    //gsInfo << "Greville: \n" << greville << "\n";

    gsBasis<> & bsp = mp.patch(side).basis().component(1);
    gsBasis<> & bsp_tilde = mp_tilde.patch(side).basis().component(1);
    gsInfo << "Basis for bsp_g0: " << bsp_tilde << "\n";

    //gsInfo << "bsp: " << bsp << "\n";
    //gsInfo << "bsp_tilde: " << bsp_tilde << "\n";

    //gsWriteParaview(bsp,"bsp",1000);
    const std::string baseName("square" + util::to_string(side));
    gsParaviewCollection collection(baseName);
    std::string fileName;

    gsMultiPatch<> BasisG1, BasisG1_right;
    for (index_t i = 2; i < 3; i++)
    {
        // Equation
        if (side == 1) // left
        {

            ev_L = bsp_tilde.evalSingle(i, greville.bottomRows(1)) // N_tilde_i
                            .cwiseProduct(
                                bsp.evalSingle(0, greville.topRows(1))
                                + bsp.evalSingle(1, greville.topRows(1)))
                + b_L.eval(greville.bottomRows(1))
                     .cwiseProduct(bsp_tilde.derivSingle(i, greville.bottomRows(1)).cwiseProduct(
                         bsp.evalSingle(1, greville.topRows(1)) / 9));

            gsInfo << "ev_L \n" << ev_L << "\n\n";

            gsFileData<> xml;

            xml <<  b_L.eval(greville.bottomRows(1));

            xml.save("beta_L_test");

            gsInfo << "ev_L: \n" << b_L.eval(greville.bottomRows(1)) << "\n\n";

            g0_temp = bsp_g0.interpolateData(ev_L, greville);

            g = dynamic_cast<gsTensorBSpline<2> &> (*g0_temp);

            gsWriteParaview(g, "basis_gL", 1000);

            gsField<> test(mp.patch(side), g);

            gsWriteParaview(test, "test_L", 1000);

            BasisG1.addPatch(g);
            fileName = baseName + util::to_string(i);
            gsWriteParaview<>(test, fileName, 10000, true);
            collection.addTimestep(fileName,i,"0.vts");

            gsInfo << fileName << "\n";
        }
        else if (side == 0) // right
        {
            ev_R = bsp_tilde.evalSingle(i, greville.bottomRows(1)) // N_tilde_i
                            .cwiseProduct(
                                bsp.evalSingle(0, greville.topRows(1)) + bsp.evalSingle(1, greville.topRows(1)))
                + b_R.eval(greville.bottomRows(1))
                     .cwiseProduct(bsp_tilde.derivSingle(i, greville.bottomRows(1)).cwiseProduct(
                         bsp.evalSingle(1, greville.topRows(1)) / 9));

            gsFileData<> xml;

            xml <<  b_R.eval(greville.bottomRows(1));

            xml.save("beta_R_test");

            g0_temp = bsp_g0.interpolateData(ev_R, greville);

            g = dynamic_cast<gsTensorBSpline<2> &> (*g0_temp);

            gsWriteParaview(g, "basis_gR", 1000);

            gsField<> test2(mp.patch(side), g);

            gsWriteParaview(test2, "test_R", 1000);

            BasisG1_right.addPatch(g);
            fileName = baseName + util::to_string(i);
            gsWriteParaview<>(test2, fileName, 10000, true);
            collection.addTimestep(fileName,i,"0.vts");
        }
    }

    collection.save();

    const std::string baseName2("FieldG1_2" + util::to_string(side));
    gsParaviewCollection collection2(baseName2);
    std::string fileName2;
    for (index_t j = 0; j < n_bar  ; j++)
    {
        if (side == 1) // left
        {
            ev1_L = a_L.eval(greville.bottomRows(1)).cwiseProduct(
                bsp_g1_bar.evalSingle(j,greville.bottomRows(1)).cwiseProduct(
                    bsp.evalSingle(1, greville.topRows(1))/9));

            //gsInfo << "alpha: " << a_L.eval(greville.bottomRows(1)) << "\n";

            g1_temp = bsp_g1.interpolateData(ev1_L, greville);
            g1 = dynamic_cast<gsTensorBSpline<2> &> (*g1_temp);

            gsField<> test3(mp.patch(side), g1);
            gsWriteParaview(test3, "g1_L", 1000);

            fileName2 = baseName2 + "g1_L" + util::to_string(j);
            gsWriteParaview<>(test3, fileName2, 1000, true);
            collection2.addTimestep(fileName2,j,"0.vts");

            BasisG1.addPatch(g1);
        }
        else if (side == 0) // right
        {
            ev1_R = a_R.eval(greville.bottomRows(1)).cwiseProduct(
                bsp_g1_bar.evalSingle(j,greville.bottomRows(1)).cwiseProduct(
                    bsp.evalSingle(1, greville.topRows(1))/9));

            g1_temp = bsp_g1.interpolateData(ev1_R, greville);
            g1 = dynamic_cast<gsTensorBSpline<2> &> (*g1_temp);

            gsField<> test3(mp.patch(side), g1);
            gsWriteParaview(test3, "g1_R", 1000);

            fileName2 = baseName2 + "g1_R" + util::to_string(j);
            gsWriteParaview<>(test3, fileName2, 1000, true);
            collection2.addTimestep(fileName2,j,"0.vts");

            // Export XML
            gsFileData<> xml;
            xml << g1;
            xml.save("BasisG1");

            BasisG1_right.addPatch(g1);
        }
    }

    for(index_t i = 2; i < 8; i++) // TODO
    {
        for(index_t j = 0; j < 8; j++)
        {
            if (side == 1)
            {
                ev = bsp_g0.component(0).evalSingle(i,greville.topRows(1)).cwiseProduct(
                    bsp_g0.component(1).evalSingle(j,greville.bottomRows(1)));

                g0_temp = bsp_g0.interpolateData(ev, greville);
                g = dynamic_cast<gsTensorBSpline<2> &> (*g0_temp);

                BasisG1.addPatch(g);

                gsField<> test3(mp.patch(side), g);
                gsWriteParaview(test3, "blabla", 1000);
            }

            if (side == 0)
            {
                ev = bsp_g0.component(0).evalSingle(i,greville.topRows(1)).cwiseProduct(
                    bsp_g0.component(1).evalSingle(j,greville.bottomRows(1)));

                g0_temp = bsp_g0.interpolateData(ev, greville);
                g = dynamic_cast<gsTensorBSpline<2> &> (*g0_temp);

                BasisG1_right.addPatch(g);
            }
        }
    }

    collection2.save();
    if (side == 1)
    {
        gsWriteParaview(BasisG1,"BasisG1",1000); // reihenfolge wichtig



        //g0_temp->continuityIncrease(1,0);
        //g0_temp->continuityIncrease(1,1);

        //gsInfo << "g: \n" << g.coefs() << "\n";

        gsField<> test(mp,BasisG1);

        gsWriteParaview(test,"test",1000);

        gsMultiBasis<> doppelt(BasisG1);


        const index_t n = BasisG1.patches().size();
        const index_t m = BasisG1.patch(1).coefs().size(); // Alle coefs gleich lange!!
        gsInfo << "n: " << n;

        gsMatrix<> matrixCoeffs(m,n);

        for(index_t i = 0; i < n; i++)
            matrixCoeffs.block(0,i,m,1) = BasisG1.patch(i).coefs();

        gsSparseMatrix<> sparseCoeffs = matrixCoeffs.sparseView();

        gsFileData<> xml;
        xml << doppelt.basis(0);
        xml << sparseCoeffs;
        xml.save("MultibasisG1_left");

        gsInfo << "Gespeichert \n";
    }

    if (side == 0)
    {
        gsMultiBasis<> doppelt(BasisG1_right);

        const index_t n = BasisG1_right.patches().size();
        const index_t m = BasisG1_right.patch(1).coefs().size(); // Alle coefs gleich lange!!
        gsInfo << "n: " << n;

        gsMatrix<> matrixCoeffs(m,n);

        for(index_t i = 0; i < n; i++)
            matrixCoeffs.block(0,i,m,1) = BasisG1_right.patch(i).coefs();

        gsSparseMatrix<> sparseCoeffs = matrixCoeffs.sparseView();

        gsFileData<> xml;
        xml << doppelt.basis(0);
        xml << sparseCoeffs;
        xml.save("MultibasisG1_right");

        gsInfo << "Gespeichert \n";
    }


} // END of g_S


void alpha_S(const gsMultiPatch<> mp,
                  const index_t side,
                  const index_t degree,
                  gsGeometry<>::Ptr & alpha) // return alpha
{
    const boundaryInterface iFace = mp.bInterface(0); // assume only one interface
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));

    // Basis for alphaS
    gsKnotVector<> kv (0,1,0,degree+1); // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> bsp(kv);

    bsp.uniformRefine(0);

    gsInfo << "The basis for alphaS is: " << bsp.detail() << "\n";
    //gsInfo << "The dim for alphaS is: " << bsp.size() << "\n";

    // Compute the minimum points for the basis and therefore
    // determine the alpha. TODO MORE FLEIXBLE
    gsMatrix<> uv1, uv2, ev1, ev2;
    gsMatrix<> greville = bsp.anchors();
    //gsInfo << "greville: " << greville <<  "\n";


    uv1.setZero(2,greville.cols());
    uv1.bottomRows(1) = greville;

    //gsInfo << "uv1: " << uv1 << "\n";

    ifaceMap.eval_into(uv1,uv2); // The correspond point on other patch u2

    //gsInfo << "---------------------------------------------\n";
    //gsInfo <<"Pair: ("<< uv1.transpose() <<"), ("<<uv2.transpose()<<")\n";

    if(side == 1) // left
    {
        const gsGeometry<> & P1 = mp.patch(iFace.first().patch); // patch right == 111

        for (index_t i = 0; i < uv1.cols(); i++)
        {
            P1.jacobian_into(uv1.col(i), ev1);
            real_t gamma1;
            gamma(uv1.col(i), gamma1);
            uv1(0, i) = gamma1 * ev1.determinant(); // erste spalte: alphaL an stelle zweiter spalte
        }

        //gsInfo << "uv1: " << uv1 << "\n";

        alpha = bsp.interpolateData(uv1.topRows(1),uv1.bottomRows(1));
    }
    else if(side == 0) // right
    {
        const gsGeometry<> & P2 = mp.patch(iFace.second().patch); // patch right == 111

        for (index_t i = 0; i < uv2.cols(); i++)
        {
            P2.jacobian_into(uv2.col(i), ev2);
            real_t gamma1;
            gamma(uv2.col(i), gamma1);
            uv2(0, i) = gamma1 * ev2.determinant(); // erste spalte: alphaL an stelle zweiter spalte
        }

        //gsInfo << "uv2: " << uv2 << "\n";

        alpha = bsp.interpolateData(uv2.topRows(1),uv2.bottomRows(1));
    }
} // END OF alpha_S

void beta(const gsMultiPatch<> mp,
                 const index_t degree,
                 gsGeometry<>::Ptr & beta) // return beta
{
    const index_t d = mp.parDim();
    const boundaryInterface iFace = mp.bInterface(0); // assume only one interface
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));
    const gsGeometry<> & P1 = mp.patch(iFace.first().patch); // patch right
    const gsGeometry<> & P2 = mp.patch(iFace.second().patch); // patch left

    gsMatrix<> uv1, uv2, ev1, ev2, D0(d,d);

    // Basis for beta
    gsKnotVector<> kv (0,1,0,degree + 1); // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> bsp(kv);

    bsp.uniformRefine(0);

    gsInfo << "The basis for beta is: " << bsp.detail() << "\n";
    //gsInfo << "The dim for beta is: " << bsp.size() << "\n";

    // Compute the minimum points for the basis and therefore
    // determine the alpha. TODO MORE FLEIXBLE
    gsMatrix<> greville = bsp.anchors();
    //gsInfo << "greville: " << greville <<  "\n";

    uv1.setZero(2,greville.cols());
    uv1.bottomRows(1) = greville;

    ifaceMap.eval_into(uv1,uv2); // The correspond point on other patch u2

    //gsInfo << "---------------------------------------------\n";
    //gsInfo <<"Pair: ("<< uv1.transpose() <<"), ("<<uv2.transpose()<<")\n";

    for(index_t i = 0; i < uv1.cols(); i++)
    {
        P1.jacobian_into(uv1.col(i),ev1);
        P2.jacobian_into(uv2.col(i),ev2);

        D0.col(0) = ev1.col(0); // (DuFL, *)
        D0.col(1) = ev2.col(0); // (*,DvFR)

        real_t gamma1;
        gamma(uv1.col(i), gamma1); // ?

        uv1(0,i) = gamma1 * D0.determinant();
    }

    //gsInfo << "uv1: " << uv1 << "\n";

    beta = bsp.interpolateData(uv1.topRows(1),uv1.bottomRows(1));
} // END of beta


void beta_L_R(const gsGeometry<>::Ptr alpha_L,
              const gsGeometry<>::Ptr alpha_R,
              const gsGeometry<>::Ptr beta_S,
              const index_t deg_beta_L,
              const index_t deg_beta_R,
              gsGeometry<>::uPtr & beta_L,
              gsGeometry<>::uPtr & beta_R)
{
    gsInfo << "Computing beta_left and beta_right... \n\n";
    // Basis for beta_left and beta_right
    gsKnotVector<> kvBeta_left (0,1,0,deg_beta_L+1); // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> bspBeta_left(kvBeta_left);
    gsKnotVector<> kvBeta_right (0,1,0,deg_beta_R+1); // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> bspBeta_right(kvBeta_right);

    gsInfo << "basis for beta_left: " << bspBeta_left.detail() << "\n";
    gsInfo << "basis for beta_right: " << bspBeta_right.detail() << "\n";

    index_t m = beta_S->basis().anchors().size();
    index_t nL = bspBeta_left.anchors().size();
    index_t nR = bspBeta_right.anchors().size();

    gsMatrix<> greville = beta_S->basis().anchors();

    //gsInfo << "Greville points: " << greville << "\n";

    gsMatrix<> Alpha_L = alpha_L->eval(greville).transpose().asDiagonal(); // m x m
    gsMatrix<> Alpha_R = alpha_R->eval(greville).transpose().asDiagonal(); // m x m

    gsMatrix<> Beta = beta_S->eval(greville).transpose(); // m x 1

    //gsInfo << "ALPHA_LEFT: " << Alpha_L  << "\n";
    //gsInfo << "ALPHA_RIGHT: " << Alpha_R  << "\n";
    //gsInfo << "Beta: " << Beta << "\n";

    // Basis function of beta_right
    gsMatrix<> N_L = bspBeta_left.eval(greville).transpose(); // m x n
    gsMatrix<> N_R = bspBeta_right.eval(greville).transpose(); // m x n
    //gsInfo << "N_Left: " << N_L << "\n";
    //gsInfo << "N_Right: " << N_R << "\n";

    // Right side of Ax = B lambda
    gsMatrix<> B_L = alpha_L->eval(greville).replicate(nL,1).cwiseProduct(N_L.transpose());
    gsMatrix<> B_R = alpha_R->eval(greville).replicate(nR,1).cwiseProduct(N_R.transpose());

    //gsInfo << "B_L: " << B_L << "\n";
    //gsInfo << "B_R: " << B_R << "\n";

    // Compute the Assemble matrix:
    gsExprAssembler<> Ass_L(1,1), Ass_R(1,1);
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> unitPatch;
    unitPatch = gsNurbsCreator<>::BSplineSquareGrid(1,1,1);

    gsMultiBasis<> unitBasis(unitPatch);
    unitBasis.uniformRefine(9);
    Ass_L.setIntegrationElements(unitBasis);
    Ass_R.setIntegrationElements(unitBasis);

    typedef gsExprAssembler<>::space       space;

    space u_L = Ass_L.getSpace(bspBeta_left); // Basis function
    space u_R = Ass_R.getSpace(bspBeta_right);

    Ass_L.initSystem();
    Ass_L.assemble(u_L * u_L.tr());

    Ass_R.initSystem();
    Ass_R.assemble(u_R * u_R.tr());

    gsMatrix<> A_L = Ass_L.matrix().toDense();
    gsMatrix<> A_R = Ass_R.matrix().toDense();

    //gsInfo << "A_R: " << A_R << "\n";
    //gsInfo << "A_L: " << A_L << "\n";

    gsSparseSolver<>::BiCGSTABILUT solver;

    //gsMatrix<> Matrix = (Alpha_L*N_R*(A_R.inverse())*B_R - Alpha_R*N_L*(A_L.inverse())*B_L);

    //gsInfo << "Matrix: \n" << Matrix << "\n";

    //gsMatrix<> lambda = Matrix.inverse() * Beta;

    //gsInfo << "lambda: " << lambda << "\n";

    //gsMatrix<> x = A_L.inverse() * B_L * lambda;

    //gsInfo << "x: " << x << "\n";

    // TEST INVERSE:

    gsMatrix<> gross(nL+nR+m,nL+nR+m), grossrhs(nL+nR+m,1);
    gross.setZero();
    gross.block(0,0,nL,nL) = 2*A_L;
    gross.block(0,nL+nR,nL,m) = -B_L;
    gross.block(nL,nL,nR,nR) = 2*A_R;
    gross.block(nL,nL+nR,nR,m) = B_R;
    gross.block(nL+nR,0,m,nL) = -Alpha_R*N_L;
    gross.block(nL+nR,nL,m,nR) = Alpha_L*N_R;

    grossrhs.setZero();
    grossrhs.block(nL+nR,0,m,1) = Beta;

    //gsInfo << "Gross: \n" << gross << "\n";

    //gsInfo << "grossrhs: \n" << grossrhs << "\n";

    // Computes a factorization (LU,QR) and solves Ax=b for the unknown x using
    // this factorization
    //x= A.partialPivLu().solve(b);
    //x= A.fullPivLu().solve(b);
    //x= A.colPivHouseholderQr().solve(b);

    gsVector<> xlambda;
    xlambda = gross.fullPivLu().solve(grossrhs);

    //gsInfo <<  "xlambda: \n" << xlambda << "\n";

    gsMatrix<> point(1,1);
    point << 0.5;

    //gsInfo << "conditiontest: \n" << alpha_L->eval(greville).transpose().asDiagonal() *  bspBeta_right.eval(greville).transpose() * xlambda.block(nL,0,nR,1)
    //    - alpha_R->eval(greville).transpose().asDiagonal() * bspBeta_left.eval(greville).transpose() * xlambda.block(0,0,nL,1)
    //    - beta_S->eval(greville).transpose() << "\n";

    beta_L = bspBeta_left.makeGeometry(xlambda.block(0,0,nL,1));
    beta_R = bspBeta_right.makeGeometry(xlambda.block(nL,0,nR,1));

    //gsInfo << "conditiontest: \n" << alpha_L->eval(greville).transpose().asDiagonal() *  beta_right->eval(greville).transpose()
    //    - alpha_R->eval(greville).transpose().asDiagonal() * beta_left->eval(greville).transpose()
    //    - beta_S->eval(greville).transpose() << "\n";

} // END OF beta_L_R
