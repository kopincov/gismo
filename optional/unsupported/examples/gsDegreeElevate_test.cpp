
#include <iostream>

#include <gismo.h>

using namespace gismo;



void degTest(gsMatrix<> coefs, int deg, int elevationAmt)
{
    static int testNum = 0;
    testNum++;
    gsInfo << "Test " << testNum << "..." << std::flush;

    // construct the spline
    int numPts = coefs.rows();
    gsKnotVector<> kv(0.0, 1.0, numPts - deg - 1, deg + 1);
    gsBSpline<> bs(kv, give(coefs));
    GISMO_ASSERT(bs.basis().size() == numPts, "Unexpected B-Spline basis");

    // construct test points and values
    gsMatrix<> evalPts = gsPointGrid<real_t>(0.0, 1.0, 100);
    gsMatrix<> originalValues = bs.eval(evalPts);

    bs.degreeElevate(elevationAmt);

    // check the new points
    gsMatrix<> newValues = bs.eval(evalPts);
    gsMatrix<> resid =  newValues - originalValues;
    GISMO_ASSERT(bs.basis().degree() == deg + elevationAmt, "Degree was not elevated correctly.");
    GISMO_ASSERT(resid.norm() < 0.0001, "Degree elevated curve did not match original curve.");

    gsInfo << " OK.\n";
}

int main()
{
    // coeffecient matrix for a B-spline with 3 control points in 2 dimensions
    gsMatrix<> coefs1(3, 2);
    coefs1 << 0, 0, 1, 0, 1, 1;

    // coefficient matrix for a B-spline with 4 control points in 3 dimensions
    gsMatrix<> coefs2(4, 3);
    coefs2 << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;

    degTest(coefs1, 1, 1);
    degTest(coefs1, 1, 2);
    degTest(coefs1, 2, 1);
    degTest(coefs1, 2, 2);
    degTest(coefs2, 1, 1);
    degTest(coefs2, 1, 2);
    degTest(coefs2, 2, 1);
    degTest(coefs2, 2, 2);
    degTest(coefs2, 3, 1);
    degTest(coefs2, 3, 2);

}



