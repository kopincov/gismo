/** @file gsGCFunctions.cpp

    @brief Searches for G1 functions on a multipatch domain

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

// Notation for functions:
std::size_t computeDimFormula(const gsMultiPatch<> & mp, bool randomize = false);

void G1Matrix(const gsMultiPatch<> & mp,
              const index_t pdeg,
              gsMatrix<> & mat,
              gsDofMapper & map);


void localBasis_incr(const gsMultiPatch<> & mp,
                     const index_t pdeg,
                     gsDofMapper & map,
                     gsMatrix<> & result, bool bdr = false);

void order_pln(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind,
    const index_t st, const index_t e);

inline gsMatrix<> columns(const gsMatrix<> & mat,
                          const gsVector<index_t> & ind)
{
    gsMatrix<> sm;
    mat.submatrixCols(ind, sm);
    return sm;
}


/*
void basisPlnDegree3(const gsMultiPatch<> & mp,
                     const index_t pdeg,
                     gsSparseMatrix<> & result);
*/

// ++++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++++
int main(int argc, char *argv[])
{
    // ========= INPUTS =========
    std::string inputFile("inputData/GCFunction.xml"); // Geometry, two 2D-Patches
    gsFileData<> fileData(inputFile); // Store the inputFile in fileData

    gsOptionList optionList; // Option List (Degree, Refinements, etc...)
    fileData.getId(0,optionList);

    gsInfo << "Input parameters are: \n" << optionList << "\n";
    // ==========================


    // ======== GEOMETRY ========
    gsMultiPatch<> multiPatch;
    fileData.getId(100,multiPatch); // Multipatch id=100
    multiPatch.computeTopology();

    const boundaryInterface & iFace = multiPatch.bInterface(0); // 0 == only one interface exists

    gsInfo << "Got Multipatch:\n" << multiPatch << "\n";
    gsInfo << "Got Interface:\n" << iFace << "\n\n";
    // ==========================


    // ======== BASIS ========
    // Refinement
    index_t numRefine = optionList.getInt("numRefine"); // number of knot insertion
    index_t polynomDegree = optionList.getInt("polynomDegree"); // number of the max. polynom degree

    // Get the degree of the input
    index_t degreeInput = multiPatch.patch(0).basis().maxDegree(); // Only for the first patch TODO

    // Regularity within patch
    index_t regularity = optionList.getInt("regularity"); // Only regularity < polynomDegree !

    // Refinements
    multiPatch.degreeElevate(polynomDegree-degreeInput);
    multiPatch.uniformRefine(numRefine,polynomDegree - regularity);

    // Create Multibasis
    gsMultiBasis<> multiBasis;
    multiBasis = gsMultiBasis<>(multiPatch);

    gsInfo << "Basis functions (of the first patch): \n" << multiBasis[0] << "\n\n";
    // =======================

    /* ==================================================
     *              Compute Dimension Formula
     *
     *
     *
     * ==================================================
     */
    index_t dimensionFormula = optionList.getInt("dimensionFormula"); // True or false

    if (dimensionFormula == 1) // 1 == True
    {
        computeDimFormula(multiPatch, false);
    }


    /* ==================================================
    *                   get Basis of G1
    *
    *
    *
    * ==================================================
    */
    gsInfo << "Working basis: \n" << multiBasis[0] << "\n";
    index_t kdim, odim;
    gsDofMapper map;
    gsMatrix<> Ker;

    index_t getbasis = optionList.getInt("getbasis"); // True or false

    if (getbasis == 1) // 1 == True
    {
        localBasis_incr(multiPatch,polynomDegree,map,Ker,false); // Last: Compute boundary functions
    }






} // +++++++++++++++++++++++++++++++ END OF MAIN +++++++++++++++++++++++++++++++


// TODO computeDimFormula()
std::size_t computeDimFormula(const gsMultiPatch<> & multiPatch, bool randomize)
{
    gsInfo << "Dimensions formula is now computing... \n\n";

    const boundaryInterface & iFace = multiPatch.bInterface(0); // 0 == only one interface exists

    index_t pdeg = multiPatch.basis(0).maxDegree(); // 0 == only first basis

    // Some notation
    gsMatrix<index_t> kDim(1,1), rValues(1,1), fDim(3,3), val(2,9);
    gsMultiPatch<> tmp;

    const index_t d = multiPatch.parDim(); // Dimension of the parameter space

    gsMatrix<> G;
    gsDofMapper map;
    G1Matrix(multiPatch,pdeg,G,map);

    //gsInfo << "G: " << G << "\n";

    kDim(0,0) = G.fullPivLu().dimensionOfKernel();
    // rank values
    //const index_t onInt = tmp.patch(tmp.iBegin()->first().patch).basis()
    //    .boundary(iFace.first().side()).size();
    rValues(0,0) = G.cols()  - kDim(0,0);
    gsInfo << "kDim: " << kDim << "\n";
    gsInfo << "rValues: " << rValues << "\n";

/*    for (index_t p = 0; p!= 3; p++) // pdeg ... pdeg+p --> p = 0,1,2
    {
        for (index_t k = 0; k!=3; k++) // k = 0,1,2
        {
            gsInfo << "p = " << pdeg + p << ", k = " << k << "\n";

            tmp = multiPatch;
            tmp.degreeElevate(p);
            tmp.uniformRefine(k,pdeg+p-1); // multiplicity = deg-1

            gsMatrix<> G;

            G1Matrix(tmp,pdeg,G,map);


        }
    } */

} // End of computeDimFormula


void G1Matrix(const gsMultiPatch<> & mp,
              const index_t pdeg,
              gsMatrix<> & mat,
              gsDofMapper & map)
{
    gsInfo << "G1 Matrix is computing... \n\n";

    const index_t d = mp.parDim();
    const boundaryInterface iFace = mp.bInterface(0); // assume only one interface
    const gsGeometry<> & P1 = mp.patch(iFace.first().patch); // patch 1
    const gsGeometry<> & P2 = mp.patch(iFace.second().patch); // patch 2
    const gsBasis<> & B1 = mp.basis(0); // basis of patch 1
    const gsBasis<> & B2 = mp.basis(1); // basis of patch 2
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));

    gsInfo << "Daten sind wie gefolgt gegeben: \n \n" <<
    "iFace: " << iFace << "\n" <<
    "P1: " << P1 << "\n" <<
    "P2: " << P2 << "\n" <<
    "B1: " << B1 << "\n" <<
    "B2: " << B2 << "\n\n";

    gsMatrix<> pt1, gr, u2, ev1, ev2, D0(d,d);
    gsMatrix<index_t> act1, act2;

    gsVector<index_t> sz(2);
    sz[0] = B1.size();
    sz[1] = B2.size();
    map = gsDofMapper(sz);
    // glue interface
    B1.matchWith(iFace, B2, act1, act2);
    map.matchDofs(iFace.first().patch, act1, iFace.second().patch, act2);
    // mark dofs
    map.markBoundary(iFace.first().patch, act1);//interface
    act1 = B1.boundaryOffset(iFace.first() .side(), 1);
    map.markBoundary(iFace.first().patch, act1); //first adj. face
    act2 = B2.boundaryOffset(iFace.second().side(), 1);
    map.markBoundary(iFace.second().patch, act2);//second adj. face
    map.finalize();
    //map.print(); // elim = second kind; free = first kind

    // the determinant lives in this basis
    const index_t l1 = iFace.first().direction(); // "orthogonal" direction p1
    const index_t l2 = iFace.second().direction(); // "orthogonal" direction p2
    gsBasis<>::uPtr detB = B1.boundaryBasis(iFace.first().side()); // Basis of boundary
    //detB->reduceContinuity(1); // Greville abscissae
    //detB->degreeElevate(d*pdeg-1);
    detB->anchors_into(gr); // gr = 0 1/k ... k-1/k 1
    //gsInfo << "gr: " << gr << "\n";
    pt1.resize(d, gr.cols());
    pt1.topRows(l1) = gr.topRows(l1); // maybe 3D
    pt1.row(l1).setConstant( iFace.first().parameter() ? 1 : 0); // 0 == false // Set zero, x-coordinates
    pt1.bottomRows(d-l1-1) = gr.bottomRows(d-l1-1); // y coordinates
    gsInfo << "Points: \n" << pt1 << "\n\n"; // Points in parameter space at the interface


    const index_t N = gr.cols();
    const index_t M = map.boundarySize();
    mat.setZero(N, M);

    gsInfo << "N: " << N << " M: " << M << "\n\n";

    gsMatrix<> u1, tmp;// temporary variable
    gr.resize(d+1,1); //gr[0]... gr[d]  --> values of alpha, beta, gamma, delta

    for (index_t i = 0; i != N; i++) // Each point on the interface
    {
        u1 = pt1.col(i); // One point on the interface u1
        ifaceMap.eval_into(u1,u2); // The correspond point on other patch u2
        P1.jacobian_into(u1,ev1);
        P2.jacobian_into(u2,ev2);
        gsInfo << "---------------------------------------------\n";
        gsInfo <<"Pair: ("<< u1.transpose() <<"), ("<<u2.transpose()<<")\n";

        gsInfo << "Jacobian DF(R): \n" << ev1 << "\n"
        << "Jacobian DF(L): \n" << ev2 << "\n\n";

        // FL = (FL1, FL2); FR = (FR1, FR2)
        // FL : [0,1] x [0,1] --> Omega^(i)
        //
        // | duFL1 duFR1 dvFL1 |
        // | duFL2 duFR2 dvFL2 | = 0
        // | dugL  dugR  dvgR  |
        //
        //
        // alpha^R * (duFL1, duFL2)^T - alpha^L * (duFR1, duFR2)^T + beta * (dvFL1, dvFL2)^T = 0
        //
        // with
        //
        // alpha^R = | DuFR DvFR |; alpha^L = | DuFL DvFR |; beta = | DuFL DuFR |
        //
        if (0) // Kathis rechnung 0 == False
        {
            // First compute alpha^R and beta
            int sgn = 1;
            D0.col(d-1) = ev2.col(l2); // (DuFL)
            //D0.col(0) = ev1.col(0); // (DuFR, *)
            for (index_t k = 0; k != d; ++k) // k = 0,1
            {
                for (index_t t = 0; t != d - 1; ++t) // t = 0,1
                    D0.col(t) = ev1.col(t + (index_t)(t >= k));
                //D0.col(1-k) = ev2.col(1-k);
                gr.at(k) = sgn * D0.determinant();
                sgn *= -1;
                //gsInfo << "D0: " << D0 << "\n"; // first alpha^L, then beta
            }
            // Second: compute alpha^R
            gr.at(d)  = sgn * ev1.determinant(); // | DuFR DvFR |
            //D0.col(0) = ev2.col(0);
            //D0.col(1) = ev1.col(0);
            //gr.at(d)  = sgn * D0.determinant();
            //gsInfo << "gr: " << gr << "\n\n"; // g(0) = alpha^R, g(1) = alpha^L, g(2) = beta

            B1.deriv_into(u1, ev1); // Evaluates the first partial derivatives of the nonzero basis function.
            B1.active_into(u1, act1);
            map.localToGlobal(act1, iFace.first().patch, act1); // Global position of the contol points
            B2.deriv_into(u2, ev2);
            B2.active_into(u2, act2);
            map.localToGlobal(act2, iFace.second().patch, act2);
            //gsInfo << "ev1: " << ev1 << "\n";
            //gsInfo << "act: " << act1 << "\n";

            for (index_t j = 0; j != act1.rows(); ++j) // == act2.rows() // Jede Basis function
            {
                const index_t jj1 = act1.at(j); // Globaler index
                if (map.is_boundary_index(jj1)) // Patch 1
                {
                    const index_t bjj = map.global_to_bindex(jj1); // Localer index

                    mat(i, bjj) +=  gr.col(0).topRows(d).dot( ev1.col(0).segment(d*j,d) );

                    // mat(i,bjj) += alpha^L *  + alpha^R *
                    // mat(i,bjj) += gr.row(0).dot(ev1.col(0).segment(d*j,1)); // col(0) wg matrix
                    //gsInfo << "bjj: " << bjj << "\n";
                    //gsInfo << "ev1: " << ev1 << "\n";
                    //gsInfo << "act1: " << ev1.col(0).segment(d*j,d) << "\n";
                }
                const index_t jj2 = act2.at(j);
                if (map.is_boundary_index(jj2)) // Patch 2
                {
                    const index_t bjj = map.global_to_bindex(jj2);

                    //gsInfo << "bjj2: " << bjj << "\n";
                    // mat(i,bjj) += beta *
                    mat(i,bjj) += ev2.at(d*j + l2) * gr.at(d);
                }
            }
        } // END kathis rechnung

        if (1)
        {
            // First compute alpha^R and beta
            int sgn = 1;
            D0.col(0) = ev1.col(0); // (DuFR, *)
            for (index_t k = 0; k != d; ++k) // k = 0,1
            {
                D0.col(1-k) = ev2.col(1-k);
                gr.at(k) = sgn * D0.determinant();
                sgn *= -1;
                gsInfo << "D0: " << D0 << "\n"; // first alpha^L, then beta
            }
            D0.col(0) = ev2.col(0);
            D0.col(1) = ev1.col(0);
            gr.at(d)  = sgn * D0.determinant();
            gsInfo << "gr: " << gr << "\n\n"; // g(0) = alpha^R, g(1) = alpha^L, g(2) = beta

            B1.deriv_into(u1, ev1); // Evaluates the first partial derivatives of the nonzero basis function.
            B1.active_into(u1, act1);
            map.localToGlobal(act1, iFace.first().patch, act1); // Global position of the contol points
            B2.deriv_into(u2, ev2);
            B2.active_into(u2, act2);
            map.localToGlobal(act2, iFace.second().patch, act2);
            //gsInfo << "ev1: " << ev1 << "\n";
            //gsInfo << "act: " << act1 << "\n";

            for (index_t j = 0; j != act1.rows(); ++j) // == act2.rows() // Jede Basis function
            {
                const index_t jj1 = act1.at(j); // Globaler index
                if (map.is_boundary_index(jj1)) // Patch 1
                {
                    const index_t bjj = map.global_to_bindex(jj1); // Localer index

                    mat(i, bjj) += gr.at(0) * ev1.at(d*j); // alpha^R * d_u F^L
                    mat(i, bjj) += gr.at(2) * ev1.at(d*j + 1); // beta * d_v F^L
                }
                const index_t jj2 = act2.at(j);
                if (map.is_boundary_index(jj2)) // Patch 2
                {
                    const index_t bjj = map.global_to_bindex(jj2);

                    mat(i, bjj) += gr.at(1) * ev2.at(d*j); // alpha^L * d_u F^R
                }
            }
        }

    }

    gsInfo << "G1 matrix: \n" << mat  << "\n";

} // End of G1Matrix


void localBasis_incr(const gsMultiPatch<> & mp,
                     const index_t pdeg,
                     gsDofMapper & map,
                     gsMatrix<> & result, bool bdr)
{
    gsMatrix<> M;
    G1Matrix(mp, pdeg, M, map); // matrix and mapper

    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mp.patch(0).basis();
    const gsBasis<> & B2 = mp.patch(1).basis();

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;
    // Returns the indices of the basis functions that are nonzero at the domain boundary.
    r1 = B1.boundaryOffset(iFace.first() .side(), 1); // left
    r2 = B1.boundaryOffset(iFace.first() .side(), 0); // interface
    r3 = B2.boundaryOffset(iFace.second().side(), 1); // right

    order_pln(iFace, r1, r2, r3, map, B1, ind, 0, B1.component(0).size());

    gsInfo << "ind: " << ind << "\n";
    gsMatrix<> Mperm, Ker;

    Mperm = columns(M,ind);

    Ker = Mperm.fullPivLu().kernel();

    result.setZero(3*r1.size(), Ker.cols() );

    for (index_t m=0; m!= Ker.rows(); ++m)
        result.row(ind[m]) = Ker.row(m);

}

void order_pln(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind,
    const index_t st, const index_t e)
{
    const index_t sz = B1.component(0).size();
    index_t m = 0;
    ind.resize(3*(e-st));
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz,1);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz,1);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz,1);

    for(index_t i=st; i < e; ++i) // e = 3, st = 0, i = 0,1,2
    {
        ind(m++) = map.bindex(rr1(i,0), iFace.first() .patch);
        ind(m++) = map.bindex(rr2(i,0), iFace.first() .patch);
        ind(m++) = map.bindex(rr3(i,0), iFace.second().patch);

        gsInfo << "bindex: " << map.bindex(rr1(i,0), iFace.first() .patch) << "\n";
        gsInfo << "bindex2: " << map.bindex(rr2(i,0), iFace.first() .patch) << "\n";
        gsInfo << "bindex3: " << map.bindex(rr3(i,0), iFace.second() .patch) << "\n";
    }
}

/*
void basisPlnDegree3(const gsMultiPatch<> & mp,
                     const index_t pdeg,
                     gsSparseMatrix<> & result)
{
    gsMatrix<> M;
    gsDofMapper map;
    G1Matrix(mp, pdeg, M, map);// matrix and mapper

    //const index_t d = mp.parDim();//=2
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mp.patch(0).basis();
    const gsBasis<> & B2 = mp.patch(1).basis();
    const index_t sP1 = B1.size();
    const index_t sP2 = B2.size();
    const index_t sz = B1.component(0).size();

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;
    // Returns the indices of the basis functions that are nonzero at the domain boundary.
    r1 = B1.boundaryOffset(iFace.first() .side(), 1); // left
    r2 = B1.boundaryOffset(iFace.first() .side(), 0); // interface
    r3 = B2.boundaryOffset(iFace.second().side(), 1); // right

    //gsInfo << "r1: " << r1 << "\n"
    //<< "r2: " << r2 << "\n"
    //<< "r3: " << r3 << "\n";

    gsVector<index_t> p(2);
    p(0) = iFace.first() .patch;
    p(1) = iFace.second().patch;

    gsInfo << "p: " << p << "\n";

    const index_t nKnots = (B1.component(1).size() - B1.component(0).degree(0)-1) /
        (B1.component(0).degree(0)-1); // = 0 (bei p=2,3)

    gsInfo << "nKnots A: " << B1.component(0).size() << "\n"
    << "nKnots B: " << B1.component(1).degree(0) << "\n";


    gsMatrix<> subm, subk;
    index_t c = 0;
    const index_t dgr=B1.component(0).degree(0);

    gsInfo << "nKnots: " << nKnots << "\n"
    << "dgr: " << dgr << "\n";

    //result = gsSparseMatrix<>(sP1+sP2, 7+2*nKnots );
    result = gsSparseMatrix<>(sP1+sP2, 7+2*nKnots +
        2*(dgr+1+nKnots*(dgr-1)-2)*(dgr+1+nKnots*(dgr-1))); // 13 bei p=2, 23 bei p=3

    std::vector<std::pair<index_t,index_t> > pi;

    gsInfo << "result: " <<  7+2*nKnots +
        2*(dgr+1+nKnots*(dgr-1)-2)*(dgr+1+nKnots*(dgr-1)) << "\n";

    // MISSING

} // End of basisPlnDegree3
*/
