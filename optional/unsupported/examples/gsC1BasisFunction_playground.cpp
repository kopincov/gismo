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


void alpha(const gsMultiPatch<> multiPatch,
           const index_t side, // 0 == right, 1 == left
           const gsMatrix<> uv, // Point uv (parameter space), where alpha is computes
           real_t & alpha); // return alpha (left for side=0, right for side=1)

void beta(const gsMultiPatch<> multiPatch,
          const gsMatrix<> uv, // Point uv (parameter space), where beta is computes
          real_t & beta); // return beta

void beta_S(const gsMultiPatch<> mp,
            const index_t side,
            const gsMatrix<> uv1, // Point uv (parameter space), where beta_S is computes
            real_t & beta_S); // return beta_S

void gamma(const gsMatrix<> uv,
           real_t & gamma)
{
    gsFunctionExpr<> f("1000/(8167 + 60*v - 407*v*v + 516*v*v*v + 407*v*v*v*v)",1);
    real_t v = uv.at(1);
    gamma = 10000/(8167 + 60*v - 407*v*v + 516*v*v*v + 407*v*v*v*v);
}

// Funktioniert nur mit poynomgrad 1
void alpha_function(const gsMultiPatch<> mp,
           gsMatrix<> & a); // return function alpha

void beta_function(const gsMultiPatch<> mp,
                   gsMatrix<> & b);

// new way
void alpha_spline(const gsMultiPatch<> mp,
                  const index_t side,
                  gsGeometry<>::uPtr & alpha); // return alpha

void beta_spline(const gsMultiPatch<> mp,
                 gsGeometry<>::uPtr & beta); // return beta

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
    fileData.getId(110,multiPatch); // Multipatch id=100 or 110
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


    // ========= Compute alpha, beta ========
    real_t alpha_L, alpha_R, beta1, gamma1, beta_L; //, beta_R;
    gsMatrix<> uv, result; // First component should be zero
    const index_t left = 0; // Patch 111
    const index_t right = 1; // Patch 112
    uv.setZero(2,1);
    uv.at(0) = 0;
    uv.at(1) = 0.5;
    real_t v = uv.at(1);

    alpha(multiPatch,left,uv,alpha_L);
    alpha(multiPatch,right,uv,alpha_R);
    beta(multiPatch,uv,beta1);
    gamma(uv,gamma1);

    gsInfo << "---------------------------------------------\n";
    gsInfo << "alpha left: " << gamma1*alpha_L << " paper: " << -3*(9+v)*0.5 << "\n"
        << "alpha right: " << gamma1*alpha_R << " paper: " << -3*(-7+v)*0.5 << "\n"
        << "beta: " << gamma1*beta1 << " paper: " << (15-32*v+v*v)/12 << "\n"
        << "gamma: " << gamma1 << " paper: " << 10000/(8167 + 60*v - 407*v*v + 516*v*v*v + 407*v*v*v*v) << "\n\n";


    beta_S(multiPatch,left,uv,beta_L);

    gsMatrix<> alpha2(2,2);
    alpha_function(multiPatch,alpha2);

    gsFunctionExpr<> alpha_function_L(util::to_string(alpha2(0,0))+" + x*" + util::to_string(alpha2(0,1)),1);
    gsFunctionExpr<> alpha_function_R(util::to_string(alpha2(1,0))+" + x*" + util::to_string(alpha2(1,1)),1);

    gsInfo << "Funktion alpha_L: " << alpha_function_L << "\n"
        << "Funktion alpha_R: " << alpha_function_R << "\n";


    gsMatrix<> beta2(1,3), uv2(1,1);

    uv2.at(0) = 0.5;

    beta_function(multiPatch,beta2);


    gsFunctionExpr<> beta_func(util::to_string(beta2(0,0))+" + x*"
        + util::to_string(beta2(0,1)) + " + x*x*" + util::to_string(beta2(0,2)),1);

    beta_func.eval_into(uv2,result);

    gsInfo << "Function beta: " << result << "\n";

    // +++++++++++++++++++++++++++++++++ new way ++++++++++++++++++++++
    gsInfo << "\n\n+++++++++++++++++++++++++++++++++ \n \n";

    // Basis for alphaS
    gsKnotVector<> kv (0,1,0,5); // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> bsp(kv);

    bsp.uniformRefine(0);

    gsInfo << "The basis for alphaS is: " << bsp.detail() << "\n";
    gsInfo << "The dim for alphaS is: " << bsp.size() << "\n";

    gsMatrix<> vals(1,3),pts(1,3), pts2(1,3), ergebnis, ergebnis2, ergebnis3, ergebnis4;

//    vals << -13.5,-13.75,-14,-14.25,-14.5,-14.75;
//    pts << 0,1/6,2/6,3/6,4/6,5/6;

    // Problem: man erwartet so viele punkte wie polynomgrad
    vals << -13.5,-14.25,-14.75; // Werte der determinante
    pts << 0,0.5,0.83333333333; // Feste Punkte

    gsMatrix<> greville = bsp.anchors();
    gsInfo << "greville: " << greville <<  "\n";

    gsMatrix<> evaluate  = alpha_function_L.eval( greville );
    gsInfo << "evaluate:" << evaluate <<  "\n";

    pts2 << 0, 1, 0.1;

    gsGeometry<>::uPtr test = (bsp.interpolateData(evaluate,greville));

    gsInfo << "interpolation: " << test->coefs() << "\n";

    test->eval_into(pts2,ergebnis);

    gsInfo << "interpolation: " << ergebnis << "\n";

    gsGeometry<>::uPtr alpha_spline_left, alpha_spline_right, beta_spline_S;

    alpha_spline(multiPatch,left,alpha_spline_left);
    alpha_spline(multiPatch,right,alpha_spline_right);

    gsInfo << "alpha_spline_left: " << alpha_spline_left->coefs() << "\n"
        << "alpha_spline_right: " << alpha_spline_right->coefs() << "\n";

    alpha_spline_left->eval_into(pts2,ergebnis2);
    alpha_spline_right->eval_into(pts2,ergebnis3);

    gsInfo << "interpolation: " << ergebnis2 << "\n";
    gsInfo << "interpolation: " << ergebnis3 << "\n";

    beta_spline(multiPatch,beta_spline_S);

    gsInfo << "beta_spline_S: " << beta_spline_S->coefs() << "\n";
    beta_spline_S->eval_into(pts2,ergebnis4);

    gsInfo << "interpolation: " << ergebnis4 << "\n";

    // Test of integral
    // Basis for beta
    gsKnotVector<> kvBeta (0,1,0,4); // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> bspBeta(kvBeta);

    gsInfo << "basis: " << bspBeta.detail() << "\n";

    gsExprEvaluator<> ev;
    gsExprAssembler<> A(1,1);

    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> unitPatch;
    unitPatch = gsNurbsCreator<>::BSplineSquareGrid(1, 1,1);

    gsMultiBasis<> unitBasis(unitPatch);
    unitBasis.uniformRefine(9);
    ev.setIntegrationElements(unitBasis);

    A.setIntegrationElements(unitBasis);

    typedef gsExprEvaluator<real_t>::variable variable;
    typedef gsExprAssembler<>::variable    variableAss;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    typedef gsExprAssembler<>::geometryMap geometryMap;
    variableAss fbetaR = A.getCoeff(alpha_spline_right->function(0));

    space u = A.getSpace(bspBeta);

//    geometryMap G = A.getMap(unitPatch);

    real_t lambda = -0.00108428;
    gsSparseSolver<>::CGDiagonal solver;
    gsMatrix<> solVectorA, solVectorB, rhsAnew, rhsBnew;

    gsInfo << "u: " << u << "\n";

    A.initSystem();
    A.assemble(u * u.tr(), lambda * u * fbetaR.val());

    solver.compute( A.matrix() );

    gsMatrix<> matrixA = A.matrix().toDense();
    gsMatrix<> rhsA = A.rhs();

    gsInfo << "inverse: " << matrixA.inverse() << "\n";

    gsMatrix<> mat1 = bspBeta.eval(bspBeta.anchors());
    gsMatrix<> mat2 = (alpha_spline_right->function(0).eval(bspBeta.anchors())).replicate<4,1>();


    rhsAnew = mat1 * mat2 ;

    gsInfo << "rhsAnew: " << alpha_spline_left->eval(bspBeta.anchors()) << "\n";


    // A.rhs_set(rhsAnew); // TODO: where is this member? @weinmueller

    solVectorA = solver.solve(A.rhs());

    variableAss fbetaL2 = A.getCoeff(alpha_spline_left->function(0));

    A.initSystem();
    A.assemble(u * u.tr(), - lambda * u * fbetaL2);

    solver.compute( A.matrix() );

    rhsBnew << lambda*13.5,0,0,0,
                1,1,1,1;

    // A.rhs_set(rhsBnew); // TODO: where is this member? @weinmueller

    solVectorB = solver.solve(A.rhs());

    gsMatrix<> lam, points;

    bspBeta.eval_into(bspBeta.anchors(),lam);

    gsMatrix<> Nx = solVectorA.transpose()*lam;
    gsMatrix<> My = solVectorB.transpose()*lam;

    //gsMatrix<> Nenner = alpha_function_L.function(0).eval(points);

    gsInfo << "Nx: " << Nx << "\n"
           << "My: " << My << "\n";



//    variable fbeta = ev.getVariable(alpha_spline_left->function(0));
//    variable fbetaL = ev.getVariable(beta_spline_S->function(0));
//    variable fbeta_allg = ev.getVariable(bspBeta.basis(0));
    //ev.integral(fbeta);
    //ev.integral(fbeta_allg);

    gsInfo << "lam: " << matrixA << "\n";
    gsInfo<< "Result old: " << rhsA <<"\n";


} // END OF MAIN

void alpha_spline(const gsMultiPatch<> mp,
           const index_t side,
           gsGeometry<>::uPtr & alpha) // return alpha
{
    const boundaryInterface iFace = mp.bInterface(0); // assume only one interface
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));

    gsInfo << "Daten sind wie gefolgt gegeben: \n \n" <<
        "iFace: " << iFace << "\n";

    // Basis for alphaS
    gsKnotVector<> kv (0,1,0,3); // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> bsp(kv);

    bsp.uniformRefine(0);

    gsInfo << "The basis for alphaS is: " << bsp.detail() << "\n";
    gsInfo << "The dim for alphaS is: " << bsp.size() << "\n";

    // Compute the minimum points for the basis and therefore
    // determine the alpha. TODO MORE FLEIXBLE
    gsMatrix<> uv1, uv2, ev1, ev2;
    gsMatrix<> greville = bsp.anchors();
    gsInfo << "greville: " << greville <<  "\n";


    uv1.setZero(2,greville.cols());
    uv1.bottomRows(1) = greville;

    gsInfo << "uv1: " << uv1 << "\n";

    ifaceMap.eval_into(uv1,uv2); // The correspond point on other patch u2

    gsInfo << "---------------------------------------------\n";
    gsInfo <<"Pair: ("<< uv1.transpose() <<"), ("<<uv2.transpose()<<")\n";

    if(side == 0) // left
    {
        const gsGeometry<> & P1 = mp.patch(iFace.first().patch); // patch right == 111

        for (index_t i = 0; i < uv1.cols(); i++)
        {
            P1.jacobian_into(uv1.col(i), ev1);
            real_t gamma1;
            gamma(uv1.col(i), gamma1);
            uv1(0, i) = gamma1 * ev1.determinant(); // erste spalte: alphaL an stelle zweiter spalte
        }

        gsInfo << "uv1: " << uv1 << "\n";

        alpha = bsp.interpolateData(uv1.topRows(1),uv1.bottomRows(1));
    }
    else if(side == 1) // right
    {
        const gsGeometry<> & P2 = mp.patch(iFace.second().patch); // patch right == 111

        for (index_t i = 0; i < uv2.cols(); i++)
        {
            P2.jacobian_into(uv2.col(i), ev2);
            real_t gamma1;
            gamma(uv2.col(i), gamma1);
            uv2(0, i) = gamma1 * ev2.determinant(); // erste spalte: alphaL an stelle zweiter spalte
        }

        gsInfo << "uv2: " << uv2 << "\n";

        alpha = bsp.interpolateData(uv2.topRows(1),uv2.bottomRows(1));
    }
}

void beta_spline(const gsMultiPatch<> mp,
        gsGeometry<>::uPtr & beta) // return beta
{
    const index_t d = mp.parDim();
    const boundaryInterface iFace = mp.bInterface(0); // assume only one interface
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));
    const gsGeometry<> & P1 = mp.patch(iFace.first().patch); // patch right
    const gsGeometry<> & P2 = mp.patch(iFace.second().patch); // patch left

    gsMatrix<> uv1, uv2, ev1, ev2, D0(d,d);

    // Basis for beta
    gsKnotVector<> kv (0,1,0,3); // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> bsp(kv);

    bsp.uniformRefine(0);

    gsInfo << "The basis for beta is: " << bsp.detail() << "\n";
    gsInfo << "The dim for beta is: " << bsp.size() << "\n";

    // Compute the minimum points for the basis and therefore
    // determine the alpha. TODO MORE FLEIXBLE
    gsMatrix<> greville = bsp.anchors();
    gsInfo << "greville: " << greville <<  "\n";

    uv1.setZero(2,greville.cols());
    uv1.bottomRows(1) = greville;

    ifaceMap.eval_into(uv1,uv2); // The correspond point on other patch u2

    gsInfo << "---------------------------------------------\n";
    gsInfo <<"Pair: ("<< uv1.transpose() <<"), ("<<uv2.transpose()<<")\n";

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

    gsInfo << "uv1: " << uv1 << "\n";

    beta = bsp.interpolateData(uv1.topRows(1),uv1.bottomRows(1));
}

void alpha(const gsMultiPatch<> mp,
           const index_t side,
           const gsMatrix<> uv1, // Point uv (parameter space), where alpha is computes
           real_t & alpha) // return alpha
{
    const boundaryInterface iFace = mp.bInterface(0); // assume only one interface
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));

    gsInfo << "Daten sind wie gefolgt gegeben: \n \n" <<
           "iFace: " << iFace << "\n";

    gsMatrix<> uv2, ev1, ev2;

    ifaceMap.eval_into(uv1,uv2); // The correspond point on other patch u2

    gsInfo << "---------------------------------------------\n";
    gsInfo <<"Pair: ("<< uv1.transpose() <<"), ("<<uv2.transpose()<<")\n";

    if (side == 0) // left
    {
        const gsGeometry<> & P1 = mp.patch(iFace.first().patch); // patch right == 111

        P1.jacobian_into(uv1,ev1);
        alpha = ev1.determinant();
    }
    else if (side == 1) // right
    {
        const gsGeometry<> & P2 = mp.patch(iFace.second().patch); // patch left == 112

        P2.jacobian_into(uv2,ev2);
        alpha = ev2.determinant();
    }

}

void beta(const gsMultiPatch<> mp,
          const gsMatrix<> uv1, // Point uv (parameter space), where beta is computes
          real_t & beta) // return beta
{
    const index_t d = mp.parDim();
    const boundaryInterface iFace = mp.bInterface(0); // assume only one interface
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));
    const gsGeometry<> & P1 = mp.patch(iFace.first().patch); // patch right
    const gsGeometry<> & P2 = mp.patch(iFace.second().patch); // patch left

    gsMatrix<> uv2, ev1, ev2, D0(d,d);

    ifaceMap.eval_into(uv1,uv2); // The correspond point on other patch u2

    P1.jacobian_into(uv1,ev1);
    P2.jacobian_into(uv2,ev2);

    gsInfo << "Daten sind wie gefolgt gegeben: \n \n" <<
           "iFace: " << iFace << "\n";

    gsInfo << "---------------------------------------------\n";
    gsInfo <<"Pair: ("<< uv1.transpose() <<"), ("<<uv2.transpose()<<")\n";

    D0.col(0) = ev1.col(0); // (DuFL, *)
    D0.col(1) = ev2.col(0); // (*,DvFR)
    beta = D0.determinant();
}

void beta_S(const gsMultiPatch<> mp,
           const index_t side,
           const gsMatrix<> uv1, // Point uv (parameter space), where beta_S is computes
           real_t & beta_S) // return beta_S
{
    gsMatrix<> points(2,1);
    points << 3,0;

    gsMultiBasis<> b(mp);
    gsFunctionExpr<> f("3*x*x + 3",2);
    gsExprEvaluator<> ev;

    gsFunctionExpr<> beta("(15-32*x+x*x)/12",1);
    gsFunctionExpr<> alphaL("-3*(9+x)/2",1);
    gsFunctionExpr<> firstTerm("(15-32*x+x*x)/12*1/(-3*(9+x)/2)",2);

    gsFunctionExpr<> newton("x+y","x-2*y",2);


    gsFunctionExpr<> f3("2*x",1);
    gsFunctionExpr<> g1(f3.expression() + " + 1",1);

    gsMatrix<> ergenbnis;

    g1.eval_into(points.transpose(),ergenbnis);


    gsInfo << "Composition: " << ergenbnis  << "\n";

    gsGeometry<>::uPtr f2 = gsNurbsCreator<>::NurbsQuarterAnnulus();
    gsVector<> x = gsVector<>::vec(0.5, 0.5),
        y = gsVector<>::vec(0.1, 1.9);

    gsVector<> arg(2), value(2);

    arg << 0.1, 0.1; // Startwert(=0) macht manchmal probleme
    value << 1,0; // rhs

    int iter = f2->newtonRaphson(y, x, true);

    newton.newtonRaphson(value,arg,false);

    gsInfo << "arg: " << arg << "\n";

    gsMatrix<> testets;
    f.eval_into(points,testets);

    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> unitPatch;
    unitPatch = gsNurbsCreator<>::BSplineSquareGrid(1, 1,1);

    b.uniformRefine(5);
    gsMultiBasis<> unitBasis(unitPatch);
    unitBasis.uniformRefine(5);
    ev.setIntegrationElements(unitBasis);


    typedef gsExprEvaluator<real_t>::variable    variable;
    typedef gsExprEvaluator<real_t>::geometryMap geometryMap;
    geometryMap G = ev.getMap(mp);
    gsInfo << "test Map \n" << G << "\n" ;
    variable v = ev.getVariable(firstTerm);
    gsInfo << "test \n" ;
    ev.min(v); // Compute min of finite elements (stored in setIntegrationElements)

    gsInfo << "ev1: " << ev.allValues() << "\n";

    //gsInfo << "minimum: " << minimum << "\n";
}

// Funktioniert nur mit poynomgrad 1
void alpha_function(const gsMultiPatch<> mp,
           gsMatrix<> & a) // return function alpha
{
    const boundaryInterface iFace = mp.bInterface(0); // assume only one interface
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));

    gsInfo << "Daten sind wie gefolgt gegeben: \n \n" <<
           "iFace: " << iFace << "\n";

    gsMatrix<> uv2, ev1, ev2;
    gsMatrix<> uv1(2,2);

    uv1 << 0, 0.5, 0, 0.1; // x y x y x y ...

    gsInfo << "uv2: " << uv1 << "\n";

    ifaceMap.eval_into(uv1.transpose(),uv2); // The correspond point on other patch u2

    gsInfo << "uv2: " << uv2 << "\n";

    const gsGeometry<> & P1 = mp.patch(iFace.first().patch); // patch right == 111

    for (index_t i = 0; i < uv1.rows(); i++)
    {
        P1.jacobian_into(uv1.row(i).transpose(), ev1);
        real_t gamma1;
        gamma(uv1.row(i).transpose(), gamma1);
        uv1(i, 0) = gamma1 * ev1.determinant(); // erste spalte: alphaL an stelle zweiter spalte
    }
    real_t b_L = (uv1(0,0)-uv1(1,0))/(uv1(0,1)-uv1(1,1));
    real_t a_L = uv1(0,0) - b_L*uv1(0,1);

    a(0,0) = a_L;
    a(0,1) = b_L;

    const gsGeometry<> & P2 = mp.patch(iFace.second().patch); // patch right == 112

    for (index_t i = 0; i < uv2.rows(); i++)
    {
        P2.jacobian_into(uv2.col(i), ev2);
        real_t gamma1;
        gamma(uv2.col(i), gamma1);
        uv2(0, i) = gamma1 * ev2.determinant(); // erste spalte: alphaL an stelle zweiter spalte
    }
    real_t b_R = (uv2(0,0)-uv2(0,1))/(uv2(1,0)-uv2(1,1));
    real_t a_R = uv2(0,0) - b_R*uv2(1,0);

    a(1,0) = a_R;
    a(1,1) = b_R;
}

// Newton Interpolationspolynom Grad 2
void beta_function(const gsMultiPatch<> mp,
                   gsMatrix<> & b)
{
    const index_t d = mp.parDim();
    const boundaryInterface iFace = mp.bInterface(0); // assume only one interface
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));
    const gsGeometry<> & P1 = mp.patch(iFace.first().patch); // patch right
    const gsGeometry<> & P2 = mp.patch(iFace.second().patch); // patch left

    gsMatrix<> uv1(2,3), uv2, ev1, ev2, D0(d,d);

    uv1 << 0,  0,  0, // u
           0, 0.5, 1; // v

    ifaceMap.eval_into(uv1,uv2); // The correspond point on other patch u2

    gsInfo << "uv1: " << uv1 << "\n"
        << "uv2: " << uv2 << "\n";

    for (index_t i = 0; i < uv1.cols(); i++)
    {
        P1.jacobian_into(uv1.col(i),ev1);
        P2.jacobian_into(uv2.col(i),ev2);

        gsInfo << "ev1: " << ev1 << "\n ev2: " << ev2 << "\n";

        D0.col(0) = ev1.col(0); // (DuFL, *)
        D0.col(1) = ev2.col(0); // (*,DvFR)

        real_t gamma1;
        gamma(uv2.col(i), gamma1);

        uv1(0,i) = gamma1 * D0.determinant();
    }


    gsInfo << "uv1: " << uv1 << "\n";

    // Polynominterpolation Grad 2
    b.at(0) = uv1(0,0); // a0
    b.at(1) = (uv1(0,1) - uv1(0,0))/(uv1(1,1) - uv1(1,0));  // a1
    b.at(2) = (uv1(0,2) - uv1(0,0) - b.at(1)*(uv1(1,2) - uv1(1,0)))/
        ((uv1(1,2) - uv1(1,0))*(uv1(1,2) - uv1(1,1)));

    // ausmultipliziert:
    b.at(0) += b.at(1) * (- uv1(1,0)) + b.at(2) * uv1(1,0) * uv1(1,1);
    b.at(1) += b.at(2) * (- uv1(1,0)) - b.at(2) * uv1(1,1);

}










