/** @file gsExprIntegral_test.cpp

    @brief Testing integral computation using the expression evaluator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

#include <gsAssembler/gsExprEvaluator.h>

using namespace gismo;

int main(int argc, char *argv[]) 
{
    //std::string fn("surfaces/sphere1.xml"); // todo: test
    std::string fn("planar/two_squares.xml");
    //std::string fn("planar/quarter_annulus_2p.xml");

    gsCmdLine cmd("Testing expression evaluator.");
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)", fn);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;
    gsReadFile<>(fn, mp);
    mp.computeTopology();
    gsMultiBasis<> b(mp);
    b.uniformRefine(1);
    //b.degreeElevate();
    //b.basis(0).component(0).uniformRefine();

    gsFunctionExpr<> a_("x+y",2);
    gsFunctionExpr<> b_("x+y",2);

    // Initiate the expression evaluator
    gsExprEvaluator<real_t> ev;

    // Set the parameter mesh as the integration mesh
    ev.setIntegrationElements(b);

    // Define integrant variables
    typedef gsExprEvaluator<real_t>::element     element;
    typedef gsExprEvaluator<real_t>::geometryMap geometryMap;
    typedef gsExprEvaluator<real_t>::variable    variable;
    element     e = ev.getElement();
    geometryMap G = ev.getMap(mp);
    variable    u = ev.getVariable(a_, G);
    variable    v = ev.getVariable(b_);
        
    //------------- Evaluation on a point grid
    
    // Construct a tensor-product point grid
    const gsMatrix<> param = mp.patch(0).parameterRange();
    gsGridIterator<real_t,CUBE> grid(param, 12);
    gsInfo<< "* Jacobian determinant values on tensor-product grid:\n";
    ev.eval( jac(G), grid );
    gsInfo<< "  Result: "<< ev.allValues().transpose() <<"\n";

    gsInfo<< "* Derivative values on tensor-product grid:\n";
    ev.eval( grad(u) * jac(G), grid );
    gsInfo<< "  Result: "<< ev.allValues().transpose() <<"\n";

    //------------- Maximum and minimum value

    gsInfo<< "* The maximum value of [ meas(G) ]:\n";
    ev.max( meas(G) );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    gsInfo<< "* The minimum value of [ meas(G) ]:\n";
    ev.min( meas(G) );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    //------------- Geometric quantities
    
    gsInfo<< "* The area of the domain [ meas(G) ]:\n";
    ev.integral( meas(G) );
    gsInfo<< "  Result: "<< ev.value() <<"\n";
    
    gsInfo<< "* The area of the domain, assuming codim=0 [ jac(G).det() ] (!) :\n";
    ev.integral( jac(G).det() );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    gsInfo<< "* The boundary area (eg. perimeter) of the domain [ nv(G).norm() ]:\n";
    // Note: meas(G) is different than nv(G).norm()
    ev.integralBdr( nv(G).norm() );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    gsInfo<< "* The area of the parameter domain [ 1 ] :\n";
    ev.integral( 1 );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    //------------- L2 Norm

    gsInfo<< "* The squared L2 norm [ u.sqr() * meas(G) ]:\n";
    ev.integral( u.sqr() * meas(G) );
    
    gsInfo<< "  Result: "<< ev.value() <<"\n";    
    ev.calcSqrt();
    gsInfo<< "  sqrt  : "<< ev.value() <<"\n";

    gsInfo<< "* The squared L2 distance [ (u-v).sqNorm() * meas(G) ]:\n";
    ev.integral( (u-v).sqNorm() * meas(G) );
    gsInfo<< "  Result: "<< ev.value() <<"\n";
    
    gsInfo<< "* The squared L2 norm ||u-v+2*v*u|| [ (u-v+2*v*u).sqNorm() * meas(G) ]:\n";
    ev.integral( (u-v+2*v*u).sqNorm() * meas(G) );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    gsInfo<< "* The squared L2 distance on the boundary [ u.sqNorm() * nv(G).norm() ]:\n";
    ev.integralBdr( u.sqNorm() * nv(G).norm() );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    // gsInfo<< "* The cubed L3 distance [  ]:\n";
    // ev.integral( (u-v).abs().pow(3)*(u-v).abs().pow(3) * meas(G) );
    // gsInfo<< "  Result: "<< ev.value() <<"\n";

    //------------- H1 Norm

    gsInfo<< "* The squared H1 seminorm [ igrad(u,G).sqNorm() * meas(G) ]:\n";
    ev.integral( igrad(u,G).sqNorm() * meas(G) );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    gsInfo<< "* The squared H1 norm [ (igrad(u,G).sqNorm()+u.sqNorm() ) * meas(G) ]:\n";
    ev.integral( (igrad(u,G).sqNorm() + u.sqNorm() ) * meas(G) );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    gsInfo<< "* The squared H1 norm per element :\n";
    ev.integralElWise( (igrad(u,G).sqNorm() + u.sqNorm() ) * meas(G) );
    gsInfo<< "  Result (elwise sum): "<< ev.allValues().sum() <<"\n";
    gsInfo<< "  Result (global)    : "<< ev.value() <<"\n";
    ev.calcSqrt();
    gsInfo<< "  sqrt (elwise): "<< ev.allValues().transpose() <<"\n";
    gsInfo<< "  sqrt (global): "<< ev.value() <<"\n";
    
    //------------- H2 Norm

    gsInfo<< "* The squared H2 seminorm [ ihess(u,G).sqNorm() * meas(G) ]:\n";
    ev.integral( ihess(u,G).sqNorm() * meas(G) );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    gsInfo<< "* The squared H2 norm [ ihess(u,G).sqNorm()+igrad(u,G).sqNorm()+u.sqNorm() ) * meas(G) ]:\n";
    ev.integral( ( ihess(u,G).sqNorm() + igrad(u,G).sqNorm() + u.sqNorm() ) * meas(G) );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    //------------- Error estimates

    gsInfo<< "* Poisson residual [ (ilapl(u,G) + v).sqr() * meas(G) ]:\n";
    // u is the trial solution
    // v is the rhs
    ev.integral((ilapl(u,G) + v).sqNorm() * meas(G) );
    gsInfo<< "  Result: "<< ev.value() <<"\n";

    gsInfo<< "* Poisson residual estimator [ e.diam().sqr()*(ilapl(u,G)+v).sqr()*meas(G) ]:\n";
    ev.integralElWise( e.diam().sqr() * (ilapl(u,G) + v).sqr() * meas(G) );
    gsInfo<< "  Result: "<< ev.allValues().transpose() <<"\n";
    gsInfo<< "  Global: "<< ev.value() <<"\n"; //== ev.allValues().sum()
    // todo: add boundary contributions and interface contributions

    return EXIT_SUCCESS;
}
