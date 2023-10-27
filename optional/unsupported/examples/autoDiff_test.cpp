/** @file autoDiff_test.cpp

    @brief Testing automatic differentiation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>
#include <gismo.h>
#include <gsAutoDiff.h> // external

using namespace gismo;

int main()
{
    // The type of our autodiff scalar
    // Second argument is number of indep. variables
    //typedef ad::DScalar2<real_t, 2> DScalar;
    typedef ad::DScalar2<real_t, -1> DScalar;

    // Note: For first derivative only use
    // typedef ad::DScalar1<real_t, -1> DScalar;

    // The type of storage of gradient/hessian
    //typedef DScalar::Gradient_t Gradient_t;
    //typedef DScalar::Gradient_t Hessian_t;

    /* Do a multidimensional newton iteration for the function
       f(x0, x1) = -log(1-x0-x1) - log(x0) - log(x1)

       starting at x0 = .85 and x1 = .05
    */
    gsVector<real_t,2> x = gsVector<>::vec(0.85, 0.05);
    
    //Define a vector of auto-diff scalars
    typedef gsVector<DScalar> DVector; 
    // Note: obscure error internal error: assertion failed at: "shared/cfe/edgcpfe/lower_name.c", line 10541
    // on the intel compiler for static vector dimension
    DVector u(2);

    gsInfo <<"Newton iteration with starting point "<< x.transpose() <<"\n";
    //gsInfo <<"Initial auto-diff vector "<< u.transpose() <<"\n";

    for (int i=0; i<10; ++i) 
    {
        /* Print the current iterate -- it should end
           up converging to [1/3, 1/3] */
        gsInfo << "Iteration " << i << "  (x=" << x.transpose() << ")" << "\n";
        
		/* Get a reference to the independent variables. These are basically
           large data structures, which contain the value, gradient and
           hessian of the two functions
           f(x0,x1) = x0   (=> gradient=[1; 0], hessian = [0 0; 0 0])
           f(x0,x1) = x1   (=> gradient=[0; 1], hessian = [0 0; 0 0])
           evaluated at the current x0,x1 values.
        */

        u[0].setVariable(0, 2, x[0]); 
        u[1].setVariable(1, 2, x[1]);
        //equiv: 
        //DScalar::Initialize(x, u);// initialize vector u

        /* Compute f(x1, x1). This additionally propagates first and
           second derivative information through builtin operators and
           math library calls.
        */
        DScalar Fx = -log(1-u[0]-u[1]) - log(u[0]) - log(u[1]);
        // equiv: DScalar Fx = -log(1-x0-x1) - log(x0) - log(x1);

        /* Print the local derivative information */
        gsInfo << "Current Fx   " << Fx << "\n";

        /* Compute the Cholesky decomposition of f's Hessian
           and use it to take a full Newton step */
        x -= Fx.getHessian().llt().solve(Fx.getGradient());
    }

    return ((x.array() - 1.0/3.0).abs() <1e-8).all() ? 0 : 1 ;
}
