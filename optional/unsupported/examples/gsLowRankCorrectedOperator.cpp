/** @file gsLowRankCorrectedOperator.cpp

    @brief Provides test examples for gsSolver/gsLowRankCorrectedOp.h,
    which implements the Sherman-Morrison-Woodbury formula

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, S. Takacs
*/

#include <gismo.h>
#include <gsSolver/gsLowRankCorrectedOp.h>
#include <gsSolver/gsKronecker.h>
#include <gsSolver/gsTimedOp.h>

#include <iostream>


using namespace gismo;

void setOnes( gsSparseMatrix<>& U, index_t m, index_t n ) {
    for( index_t i=0; i<m; i++ )
        for( index_t j=0; j<n; j++ )
            U(i,j)=1.;
}


bool testInvertingLowRankCorrectedOperator( index_t d )
{
    gsInfo << "Test 1 with d = " << d << " ..." << "\n";

    gsSparseMatrix<real_t> U(d,1), V(d,1), Q(1,1);
    gsSparseMatrix<real_t> A(d,d);
    gsMatrix<real_t> x(d,1), result(d,1), x1(d,1);

    gsInfo << "    Setup A, Q and U ... " << std::flush;
    setOnes( U, d, 1 );
    setOnes( V, d, 1 );
    Q.setIdentity();

    for( index_t i=0; i<d; i++ )
        A(i,i) = 2.;

    gsInfo << "done.\n    Setup Ainv ... " << std::flush;
    gsLinearOperator<>::Ptr Ainv = gsTimedOp<>::make( "Cholesky(A)", makeSparseCholeskySolver( A ) ); 

    gsInfo << "done.\n    Setup SMW ... " << std::flush;
    gsLinearOperator<>::Ptr SMW = gsLowRankCorrectedOp<>::make( Ainv, Q, U, V );

    gsInfo << "done.\n    Setup x ... " << std::flush;
    x.setRandom( d, 1 );

    gsInfo << "done.\n" << std::flush;
    
    if( d < 10. ) gsInfo << "    x0 : " << x.transpose() << std::endl;

    gsInfo << "    Apply SMW ... " << std::flush;
    SMW->apply( x, result );

    gsInfo << "done.\n" << std::flush;
    if( d < 10. ) gsInfo << "    res: " << result.transpose() << std::endl;

    gsInfo << "    Compute x1 ... " << std::flush;
    x1.noalias() = A * result - U * U.transpose() * result; // no Q here as its only 1
    gsInfo << "done.\n" << std::flush;

    if( d < 10. ) gsInfo << "    x1 : " << x1.transpose() << std::endl;

    real_t err = (x-x1).norm();

    gsInfo << "    Error: " << err << std::endl;

    bool p = err <= 1.e-7;

    if( p )
        gsInfo << "...test passed." << std::endl;
    else
        gsInfo << "...test failed." << std::endl;

    return p;
}

bool testInvertingLowRankCorrectedOperatorTwo()
{
    gsInfo << "Test 2 with d = 9 ..." << "\n";

    gsSparseMatrix<real_t> U(9,2), Q(2,2);
    gsSparseMatrix<real_t> A(9,9);
    gsMatrix<real_t> x(9,1), result(9,1), x1(9,1);

    gsInfo << "    Setup A, Q and U ... " << std::flush;
    setOnes( U, 9, 2 );
    Q.setIdentity();

    for( index_t i=0; i<5; i++ )
    {
        U(i,1) = 0.;
    }
    
    for( index_t i=0; i<9; i++ )
        A(i,i) = 2.;

    gsInfo << "done.\n    Setup Ainv ... " << std::flush;
    gsLinearOperator<>::Ptr Ainv = gsTimedOp<>::make( "Cholesky(A)", makeSparseCholeskySolver( A ) ); 

    gsInfo << "done.\n    Setup SMW ... " << std::flush;
    gsLinearOperator<>::Ptr SMW = gsLowRankCorrectedOp<>::make( Ainv, Q, U, U );

    gsInfo << "done.\n    Setup x ... " << std::flush;
    x.setRandom( 9, 1 );

    gsInfo << "done.\n" << std::flush;
    
    gsInfo << "    x0 : " << x.transpose() << std::endl;

    gsInfo << "    Apply SMW ... " << std::flush;
    SMW->apply( x, result );

    gsInfo << "done.\n" << std::flush;
    gsInfo << "    res: " << result.transpose() << std::endl;

    gsInfo << "    Compute x1 ... " << std::flush;
    x1.noalias() = A * result - U * U.transpose() * result; // no Q here as its only 1
    gsInfo << "done.\n" << std::flush;

    gsInfo << "    x1 : " << x1.transpose() << std::endl;

    real_t err = (x-x1).norm();

    gsInfo << "    Error: " << err << std::endl;

    bool p = err <= 1.e-7;

    if( p )
        gsInfo << "...test passed." << std::endl;
    else
        gsInfo << "...test failed." << std::endl;

    return p;
}

bool testInvertingLowRankCorrectedOperatorThree( index_t d )
{
    gsInfo << "Test 3 with d = " << d << " ..." << "\n";

    gsSparseMatrix<real_t> U(d,1), Q(1,1);
    gsSparseMatrix<real_t> A(d,d);
    gsMatrix<real_t> x(d,1), result(d,1), x1(d,1);

    gsInfo << "    Setup A, Q and U ... " << std::flush;
    setOnes( U, d, 1 );
    Q.setIdentity();

    for( index_t i=0; i<d; i++ )
        A(i,i) = 2.;

    for( index_t i=0; i<d-1; i++ )
    {
        A(i,i+1) = -1.;
        A(i+1,i) = -1.;
    }

    gsInfo << "done.\n    Setup Ainv ... " << std::flush;
    gsLinearOperator<>::Ptr Ainv = gsTimedOp<>::make( "Cholesky(A)", makeSparseCholeskySolver( A ) ); 

    gsInfo << "done.\n    Setup SMW ... " << std::flush;
    gsLinearOperator<>::Ptr SMW = gsLowRankCorrectedOp<>::make( Ainv, Q, U, U );

    gsInfo << "done.\n    Setup x ... " << std::flush;
    x.setRandom( d, 1 );

    gsInfo << "done.\n" << std::flush;
    
    if( d < 10. ) gsInfo << "    x0 : " << x.transpose() << std::endl;

    gsInfo << "    Apply SMW ... " << std::flush;
    SMW->apply( x, result );

    gsInfo << "done.\n" << std::flush;
    if( d < 10. ) gsInfo << "    res: " << result.transpose() << std::endl;

    gsInfo << "    Compute x1 ... " << std::flush;
    x1.noalias() = A * result - U * U.transpose() * result;// no Q here as its only 1
    gsInfo << "done.\n" << std::flush;

    if( d < 10. ) gsInfo << "    x1 : " << x1.transpose() << std::endl;

    real_t err = (x-x1).norm();

    gsInfo << "    Error: " << err << std::endl;

    bool p = err <= 1.e-7;

    if( p )
        gsInfo << "...test passed." << std::endl;
    else
        gsInfo << "...test failed." << std::endl;

    return p;
}


bool testInvertingLowRankCorrectedOperatorFour( index_t d )
{
    gsInfo << "Test 4 with d = " << d << " ..." << "\n";

    gsSparseMatrix<real_t> U(d*d,1), Q(1,1);
    gsSparseMatrix<real_t> A(d,d);
    gsMatrix<real_t> x(d*d,1), result(d*d,1), x1(d*d,1);

    gsInfo << "    Setup A, Q and U ... " << std::flush;
    setOnes( U, d*d, 1 );
    Q.setIdentity();

    for( index_t i=0; i<d; i++ )
        A(i,i) = 2.;

    for( index_t i=0; i<d-1; i++ )
    {
        A(i,i+1) = -1.;
        A(i+1,i) = -1.;
    }

    gsInfo << "done.\n    Setup Ainv ... " << std::flush;
    gsLinearOperator<>::Ptr Ainv1 = gsTimedOp<>::make( "Cholesky(A1)", makeSparseCholeskySolver( A ) ); 
    gsLinearOperator<>::Ptr Ainv2 = gsTimedOp<>::make( "Cholesky(A2)", makeSparseCholeskySolver( A ) ); 

    gsInfo << "done.\n    Setup SMW ... " << std::flush;
    gsLinearOperator<>::Ptr SMW = gsLowRankCorrectedOp<>::make( gsTimedOp<>::make( "Cholesky( A(x)A )", gsKroneckerOp<>::make( Ainv1, Ainv2 ) ), Q, U, U );

    gsInfo << "done.\n    Setup x ... " << std::flush;
    x.setRandom( d*d, 1 );

    gsInfo << "done.\n" << std::flush;
    
    if( d*d < 10. ) gsInfo << "    x0 : " << x.transpose() << std::endl;

    gsInfo << "    Apply SMW ... " << std::flush;
    SMW->apply( x, result );

    gsInfo << "done.\n" << std::flush;
    if( d*d < 10. ) gsInfo << "    res: " << result.transpose() << std::endl;

    gsSparseMatrix<real_t> Akron = A.kron(A);
    
    gsInfo << "    Compute x1 ... " << std::flush;
    x1.noalias() = Akron * result - U * U.transpose() * result; // no Q here as its only 1
    gsInfo << "done.\n" << std::flush;

    if( d*d < 10. ) gsInfo << "    x1 : " << x1.transpose() << std::endl;

    real_t err = (x-x1).norm();

    gsInfo << "    Error: " << err << std::endl;

    bool p = err <= 1.e-4;

    if( p )
        gsInfo << "...test passed." << std::endl;
    else
        gsInfo << "...test failed." << std::endl;

    return p;
}

int main(int argc, char *argv[])
{
    gsCmdLine cmd("Tests gsLowRankCorrectedOp.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    srand( (unsigned)time( NULL ) );

    bool ok = true;

    ok = testInvertingLowRankCorrectedOperator( 9 ) && ok;
    //ok = testInvertingLowRankCorrectedOperator( 1000 ) && ok;
    ok = testInvertingLowRankCorrectedOperatorTwo() && ok;
    ok = testInvertingLowRankCorrectedOperatorThree( 9 ) && ok;
    //ok = testInvertingLowRankCorrectedOperatorThree( 1000 ) && ok;
    ok = testInvertingLowRankCorrectedOperatorFour( 9 ) && ok;
    //ok = testInvertingLowRankCorrectedOperatorFour( 100 ) && ok;

    return (ok ? 0 : 1);
}


