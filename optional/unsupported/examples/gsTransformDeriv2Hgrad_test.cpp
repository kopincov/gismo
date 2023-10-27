/** @file gsTransformDeriv2Hgrad_test.cpp

    @brief Unit testing the gsTransformDeriv2Hgrad function in gsGeometryEvaluator.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger, A. Bressan
*/

#include <gismo.h>
#include <iostream>
#include <gsCore/gsGeometryEvaluator.hpp>

using namespace gismo;


// helper function that gives back a geometry. it is a square
// parameterized by [s,t]->(s^2,t^2)
gsTensorBSpline<2,real_t> getGeom()
{
    gsKnotVector<real_t> knots(0,1,0,3,0);
    gsMatrix<real_t> coefs(9,2);
    coefs << 0, 0,  0, 0,  1, 0,
             0, 0,  0, 0,  1, 0,
             0, 1,  0, 1,  1, 1;
    gsTensorBSplineBasis<2,real_t> basis(knots,knots);
    return gsTensorBSpline<2,real_t>(basis,coefs);
}

// helper function that gives back a geometry. it is a square
// parameterized by [s,t]->(s^2,t,s^2)
gsTensorBSpline<2,real_t> getGeomSurf()
{
    gsKnotVector<real_t> knots(0,1,0,3,0);
    gsMatrix<real_t> coefs(9,3);
    coefs << 0, 0,   0,    0, 0,   0,    1, 0,   1,
             0, 0.5, 0,    0, 0.5, 0,    1, 0.5, 1,
             0, 1,   0,    0, 1,   0,    1, 1,   1;
    gsTensorBSplineBasis<2,real_t> basis(knots,knots);
    return gsTensorBSpline<2,real_t>(basis,coefs);
}

struct TestData{
    typename gsGeometry<real_t>::uPtr geom;

    gsFunctionExpr<real_t>* df_hat;
    gsFunctionExpr<real_t>* ddf_hat;
    gsFunctionExpr<real_t>* ddf;

    ~TestData()
    {
        delete df_hat;
        delete ddf_hat;
        delete ddf;
    }

    // function that tests the transform with given geometry and parametric function derivatives.
    bool testTransformLaplaceHgrad(real_t tol)
    {
        bool passed=true;
        typename gsGeometryEvaluator<real_t>::uPtr eval(getEvaluator(NEED_MEASURE|NEED_VALUE|NEED_GRAD_TRANSFORM|NEED_2ND_DER, *geom));
        gsMatrix<real_t> result;
        gsVector<real_t,2> start,end;
        start << 0.1,0.1;
        end << 0.9,0.9;
        gsVector<unsigned,2> np;
        np << 10,10;
        gsMatrix<real_t> u=gsPointGrid<real_t>(start,end,np);
        eval->evaluateAt(u);
        gsMatrix<real_t> dfeval = df_hat->eval(u);
        gsMatrix<real_t> ddfeval = ddf_hat->eval(u);
        for(int j = 0;j<u.cols();++j)
        {
            eval->transformDeriv2Hgrad(j,dfeval,ddfeval,result);
            passed = passed&&gsAllCloseAbsolute(result,(ddf->eval(eval->value(j))).transpose(),tol);
        }

        return passed;
    }
};


int main(int argc, char *argv[])
{
    bool passed = true;
    std::vector<TestData*> testdata;
    real_t tol = 2000*math::limits::epsilon();

    // First test:
    gsTensorBSpline<2,real_t>square2D=getGeom();
    //gsFunctionExpr<real_t> f_hat1("exp(x^2)+exp(y^2)",2);
    TestData test1;
    test1.geom= square2D.clone();
    test1.df_hat=new gsFunctionExpr<real_t>("2*x*exp(x^2)","2*y*exp(y^2)",2);
    test1.ddf_hat=new gsFunctionExpr<real_t>("2*exp(x^2)+4*x^2*exp(x^2)","2*exp(y^2)+4*y^2*exp(y^2)","0",2);
    test1.ddf=new gsFunctionExpr<real_t>("exp(x)","exp(y)","0",2);
    testdata.push_back(&test1);

    // Second test:
    //gsFunctionExpr<real_t> f_hat2("y*exp(x^2)",2);
    TestData test2;
    test2.geom= square2D.clone();
    test2.df_hat=new gsFunctionExpr<real_t>("2*x*y*exp(x^2)","exp(x^2)",2);
    test2.ddf_hat=new gsFunctionExpr<real_t>("(2+4*x^2)*y*exp(x^2)","0","2*x*exp(x^2)",2);
    test2.ddf=new gsFunctionExpr<real_t>("sqrt(y)*exp(x)","-exp(x)/(4*y*sqrt(y))","exp(x)/(2*sqrt(y))",2);
    testdata.push_back(&test2);

    // Third test:
    gsTensorBSpline<2,real_t>square3D=getGeomSurf();
    //gsFunctionExpr<real_t> f_hat3("4*x^4",2);
    TestData test3;
    test3.geom= square3D.clone();
    test3.df_hat=new gsFunctionExpr<real_t>("16*x^3","0",2);
    test3.ddf_hat=new gsFunctionExpr<real_t>("48*x^2","0","0","0",2);
    std::vector<std::string> exprVec;
    exprVec.push_back("2");
    exprVec.push_back("0");
    exprVec.push_back("2");
    exprVec.push_back("0");
    exprVec.push_back("2");
    exprVec.push_back("0");
    test3.ddf=new gsFunctionExpr<real_t>(exprVec,3);
    testdata.push_back(&test3);

    for(unsigned i = 0;i<testdata.size();++i)
    {
        bool passed_i = testdata[i]->testTransformLaplaceHgrad(tol);
        passed = passed && passed_i;
        gsInfo << "Test " << i << " has ";
        if(!passed_i)
            gsInfo << "not ";
        gsInfo << "passed.\n";
    }
    return !passed;
}
