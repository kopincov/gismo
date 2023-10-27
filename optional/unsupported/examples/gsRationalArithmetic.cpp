/** @file gsRationalArithmetic.cpp

    @brief Tutorial for using GMP within G+Smo.
    Focus on rational numbers, examples as convertion from other types, basic mathematical expressions,
    matrizes and solving linear systems.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): K. Birner A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{

#ifdef GISMO_WITH_GMP
  
    // increase the printing accuracy of the output
    gsInfo.precision(20);

    
    double tmp(0);
    tmp = 0.5;
    gsInfo <<tmp <<"\n";

    // floating point with multiprecision
    mpf_class ff(tmp);
    gsInfo <<ff <<"\n";

    // integer number
    mpq_class zz("1234567891234567812345678912345678");
    gsInfo <<zz <<"\n";

    // real_t --> is actually mpq_class here (GISMO_COEFF_TYPE)
    // because we compiled with mpq_class
    // make rational from double
    real_t a(tmp);
    gsInfo <<a <<"\n";

    // make rational from numerator + denominator
    real_t b(16, 10000);
    gsInfo <<b <<"\n";
    // we can canonicalize the number
    b.canonicalize();
    gsInfo <<b <<"\n";

    //basic mathematical functions work intuitively
    real_t p = a*b;
    real_t n = a+b;
    real_t o = a/b;
    gsInfo << "a*b= " << p <<"\n";
    gsInfo << "a+b= " << n <<"\n";
    gsInfo << "a/b= " << o <<"\n";

    //canonicalizing is ALWAYS important!
    p.canonicalize();
    n.canonicalize();
    o.canonicalize();
    gsInfo << "Simplified a*b= " << p <<"\n";
    gsInfo << "Simplified a+b= " << n <<"\n";
    gsInfo << "Simplified a/b= " << o <<"\n";

    //Comparing two variables
    gsInfo << "Is a < b? " << ((a<b) ? "Yes" : "No") << "\n";

    // Get a double approximation of the rational
    gsInfo <<b.get_d() <<"\n";

    // get the numerator
    gsInfo <<b.get_num() <<"\n";

    // get the denominator
    gsInfo <<b.get_den() <<"\n";

    // Produce random rational numbers
    gmp_randclass X(gmp_randinit_default);
    X.seed(time(0));

    // fill in a (dense) matrix with random numbers
    index_t s = 5;
    gsMatrix<> m(s,s);
    for (index_t i = 0; i!=m.size(); ++i)
    {
        m.at(i) = mpq_class(X.get_z_bits(20),X.get_z_bits(20));
        m.at(i).canonicalize(); // !important 
    }

    // print out the matrix
    gsDebugVar( m );

    // make a vector full of ones
    gsVector<> bb = gsVector<>::Ones(s);
    gsDebugVar( bb.transpose() );

    // Solve the system m*x=b with LU
    gsMatrix<> x = m.fullPivLu().solve(bb);
    gsDebugVar(x.transpose());

    // Check that the solution is correct (must be == b)
    gsDebugVar( (m*x).transpose() );

    // We make a sparse matrix
    gsSparseMatrix<> sm = m.sparseView(); //(epsilon,reference);

    // solve again with LU
    gsSparseSolver<>::LU slv;
    x = slv.compute(sm).solve(bb);
    // check the result
    gsDebugVar( (sm*x).transpose() );

#else
    gsInfo <<"Enable the GISMO_WITH_GMP option to execute this test.\n";
#endif
    return 0;
}
