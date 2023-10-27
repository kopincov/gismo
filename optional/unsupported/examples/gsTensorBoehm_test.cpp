// Functions test gsTensorBoehm (knot insertion for arbitrary dimension).
// I assume that gsBoehm (knot insertion with one dimension) is correct and
// test gsTensorBoehm function with a help of the gsBoehm function.
//
// I test only 2D and 3D case.
//
// Author: Jaka Speh

#include <iostream>
#include <vector>

#include <gismo.h>


using namespace gismo;
using gismo::bspline::getIndex;
using gismo::bspline::getLastIndex;


// TESTING FUNCTION
// function for printing std::vector in gsInfo
template <typename TT>
void _print_std_vector(std::vector<TT> vec)
{
    for (typename std::vector<TT>::iterator it = vec.begin(); it != vec.end(); ++it)
    {
        gsInfo << *it << " ";
    }

    gsInfo << "\n";

}



// TESTING FUNCTION
// Builds strides from vector of dimension of the coefficients. Inverse
// function of the dimensionsOfCoeff
//
// \param dim_of_coeff - vector of the dimensions og the coefficients
//
// \returns vector of the strides
//
// \example
// buildStrides([2, 4, 6])
// => [1, 2, 8]
gsVector<unsigned> buildStrides(const std::vector<unsigned>& dim_of_coeff)
{

    unsigned length = dim_of_coeff.size();
    gsVector<unsigned> strides(length);

    strides[0] = 1;
    for (unsigned i = 1; i < length; ++i)
    {
        strides[i] = strides[i - 1] * dim_of_coeff[i - 1];
    }

    return strides;
}



/*
// TESTING FUNCTION
// builds coefficients matrix which has
// \param x coefficients in x direction
// \param y coefficients in y direction
// \param z coefficients in z direction
//
// Example:
// _buildTestCoeffs(2, 4, 6)
//
// returns matrix (output below has two matrix rows on each line):
//
// (0   0   0),  (1   0   0)  \
// (0   1   0),  (1   1   0)   \  z == 0
// (0   2   0),  (1   2   0)   /
// (0   3   0),  (1   3   0)  /
//
// (0   0   1),  (1   0   1) \
// (0   1   1),  (1   1   1)  \  z == 1
// (0   2   1),  (1   2   1)  /
// (0   3   1),  (1   3   1) /
//
// and rows for z = 2, 3, 4, 5.
*/
gsMatrix<> _buildTestCoeffs(int X, int Y, int Z)
{
    gsMatrix<> coef(X * Y * Z, 3);

    for (int z = 0; z < Z; z++)
        for (int y = 0; y < Y; y++)
            for (int x = 0; x < X; x++)
                coef.row(z * X * Y + y * X + x) << x * 1.0, y * 1.0, z * 1.0;

    return coef;
}



// TESTING FUNCTION
// builds coefficients matrix which has
// \param x coefficients in x direction
// \param y coefficients in y direction
//
// Example:
// _buildTestCoeffs2D(2, 4)
//
// returns matrix (output below)
//
// (0, 0), (1, 0)
// (0, 1), (1, 1)
// (0, 2), (1, 2)
// (0, 3), (1, 3)
gsMatrix<> _buildTestCoeffs2D(int X, int Y)
{
    gsMatrix<> coef(X * Y, 2);

    for (int y = 0; y < Y; y++)
        for (int x = 0; x < X; x++)
                coef.row(y * X + x) << x * 1.0, y * 1.0;

    return coef;
}



// TESTING FUNCTION
// tests gsTensorBoehm function for 3D and 2D case
//
// \param kv1 - knot vector in x direction
// \param kv2 - knot vector in y direction
// \param kv3 - knot vector in z direction
// \param direction - direction in which we insert a knot
// \param val - value of the knot, we will insert
// \param mult - how many times will we insert a knot
// \param d - dimension (should work for 2 or 3)
// \param output - if we print anything on gsInfo
//
// \return number of failed tests
//
// It it not the prettiest function on the world, but it does its job.
// If you are testing for dimension 2 (d = 2), then you must provide a
// dummy knot vector 3 (kv3) for third dimension (that's because I
// added test for 2D case later). Sorry :-).
unsigned _testTensorBoehm(gsKnotVector<> kv1,
                          gsKnotVector<> kv2,
                          gsKnotVector<> kv3,
                          int direction,
                          real_t val,
                          int mult,
                          int d = 3,
                          bool output = true)
{
    GISMO_ASSERT(direction < d,
                 "Can not insert into non existing direction.");

    // size of points in certain direction
    std::vector<unsigned> dim_of_coeff(d);
    std::vector<gsKnotVector<> > knots(d);

    knots[0] = kv1;
    knots[1] = kv2;

    if (2 < d)
        knots[2] = kv3;


    // building dimension of the coefficient (number of points) and the strides
    dim_of_coeff[0] = kv1.size() - kv1.degree() - 1;
    dim_of_coeff[1] = kv2.size() - kv2.degree() - 1;

    if (2 < d)
        dim_of_coeff[2] = kv3.size() - kv3.degree() - 1;

    std::vector<unsigned> dim_of_coeff2(dim_of_coeff);
    dim_of_coeff2[direction] += mult;


    gsVector<unsigned> strides = buildStrides(dim_of_coeff);
    gsVector<unsigned> strides2 = buildStrides(dim_of_coeff2);


    if (output)
    {
        gsInfo << "\n\n=============================================\n"
                  << "FUNCTION: _testTensorBoehm:\n"
                  << "Knot vectors:" << "\n"
                  << kv1.detail() << "\n"
                  << kv2.detail() << "\n"
                  << kv3.detail() << "\n"
                  << "direction: " << direction
                  << "\nknot value: " << val
                  << "\nknot multiplicity: " << mult
                  << "\ndimension: " << d << "\n";
    }

    // control points
    gsMatrix<> coeff;
    if (d == 3)
    {
        coeff = _buildTestCoeffs(dim_of_coeff[0],
                                 dim_of_coeff[1],
                                 dim_of_coeff[2]);
    }
    else if (d == 2)
    {
        coeff = _buildTestCoeffs2D(dim_of_coeff[0], dim_of_coeff[1]);
    }
    else
    {
        gsInfo<<"Function only works for d = 2 or d = 3.\n";
        return 1;
    }

    gsMatrix<> coeff_copy(coeff);

    // number of failed insertions
    unsigned failed = 0;

    if (output)
    {
        gsInfo << "Dimension of coefficients"
                  << "(number of points in directions): ";
        _print_std_vector<unsigned>(dim_of_coeff);
        gsInfo << "Strides: ";
        gsInfo << strides.transpose();
        gsInfo << "\nCalling gsTensorBoehm function.\n" << "\n";
    }

    gsTensorBoehm(knots[direction], coeff, val, direction, strides, mult, false);

    // number of control points in given direction
    unsigned t = dim_of_coeff[direction];

    // indices
    gsVector<int> index2(d);
    for (int i = 0; i < d; i++)
        index2[i] = 0;

    gsVector<int> index(index2);

    // necessary for computation of the indices
    gsVector<int> first(index2);
    gsVector<int> last(d);
    gsVector<int> last2(d);
    getLastIndex(strides2, coeff.rows(), last2);
    last2[direction] = 0;
    getLastIndex(strides, coeff_copy.rows(), last);
    last[direction] = 0;

    do
    {

        if (output)
        {
            gsInfo << "-----------------------------" << "\n";
            gsInfo << "-----------------------------" << "\n";
            gsInfo << "-----------------------------" << "\n";
        }

        int ind = getIndex(strides, index);
        gsMatrix<> coeff1(t, d);

        // copy proper coefficient to do 1D Boehm algorithm
        for (unsigned i = 0; i < t; ++i)
        {
            coeff1.row(i) = coeff_copy.row(ind + i * strides[direction]);

        }

        gsBoehm<real_t, gsKnotVector<>, gsMatrix<> >(
                    knots[direction],
                    coeff1,
                    val,
                    mult,
                    false);


        // collect proper coefficients for comparison with gsBoehm
        unsigned ind2 = getIndex(strides2, index2);
        nextCubePoint(index2, first, last2);

        gsMatrix<> temporary_mat(dim_of_coeff2[direction], d);


        for (unsigned i = 0; i < dim_of_coeff2[direction]; ++i)
        {
            temporary_mat.row(i) = coeff.row(ind2 + i * strides2[direction]);
        }


        // looking at the norm of difference betwen two coefficients matrices
        gsMatrix<> difference = coeff1 - temporary_mat;
        bool tmp_bol = difference.norm() < 0.00000000001;


        if (output)
        {
            gsInfo << "Coefficients gsBoehm:\n"
                      << coeff1;

            if (tmp_bol)
            {
                gsInfo << "\nSUCCESS, norm (difference "
                          << "between gsBoehm and geTensorBoehm): "
                          << difference.norm() << "\n";
            }
            else
            {
                gsInfo << "\n----------------------------" << "\n";
                gsInfo << "Coefficients gsTensorBoehm:\n";

                gsInfo << "\nAre equal: " << tmp_bol
                          << "\nDifference norm: " << difference.norm()
                          << "\n" << difference << "\n";
            }
        }

        if (!tmp_bol)
            failed++;

    } while(nextCubePoint(index, first, last));

    if (output)
    {
        gsInfo << "\nFailed tries: " << failed << "\n" << "\n";
    }

    return failed;
}



// TESTING FUNCTION
// function tests nextPointIndex
//
// \param direction - direction we will omit in nextPointIndex function
//
void _testNextCubePoint(int direction)
{

    std::vector<int> position (3, 0);

    std::vector<int> start (3, 0);

    std::vector<int> end (3, 0);

    end[0] = 1;
    end[1] = 3;
    end[2] = 5;

    end[direction] = 0;

    do
    {
        _print_std_vector<int>(position);

    } while (nextCubePoint<std::vector<int> > (position, start, end));

}



// TESTING FUNCTION
// Tests tensor Boehm algorithm in several cases for 2D and 3D.
//
// \param output - if we want any output
//
// \return number of failed tests
unsigned mainTestTensorBoehm(bool output)
{
    unsigned failed = 0;

    gsKnotVector<> kv1(0, 1, 0, 2);
    gsKnotVector<> kv2(0, 1, 1, 3);
    gsKnotVector<> kv3(0, 1, 2, 4);
    gsKnotVector<> kv4(0, 1, 1, 4);



    // 3D tests

    // insert a knot, with multiplicity 0, 1-time
    failed += _testTensorBoehm(kv1, kv2, kv3, 0, 0.5, 1, 3, output);

    // insert a knot, with multiplicity 0, 2-times
    failed += _testTensorBoehm(kv1, kv2, kv3, 1, 0.2, 2, 3, output);

    // insert a knot, with multiplicity 1, 1-time
    failed += _testTensorBoehm(kv1, kv2, kv3, 1, 0.5, 1, 3, output);

    // insert a knot, with multiplicity 0, 3 times
    failed += _testTensorBoehm(kv1, kv2, kv3, 2, 0.5, 3, 3, output);

    // insert a knot, with multiplicity 0, 3 times
    failed += _testTensorBoehm(kv1, kv2, kv4, 2, 0.2, 3, 3, output);

    // insert a knot, with multiplicity 1, 2-times
    failed += _testTensorBoehm(kv1, kv2, kv4, 2, 0.5, 2, 3, output);


    // 2D tests

    // insert a knot, with multiplicity 0, 1-time
    failed += _testTensorBoehm(kv1, kv2, kv3, 0, 0.5, 1, 2, output);

    // insert a knot, with multiplicity 0, 2-times
    failed += _testTensorBoehm(kv1, kv2, kv3, 1, 0.2, 2, 2, output);

    // insert a knot, with multiplicity 1, 1-time
    failed += _testTensorBoehm(kv1, kv2, kv3, 1, 0.5, 1, 2, output);

    // insert a knot, with multiplicity 0, 3 times
    failed += _testTensorBoehm(kv1, kv3, kv2, 1, 0.5, 3, 2, output);

    // insert a knot, with multiplicity 0, 3 times
    failed += _testTensorBoehm(kv1, kv4, kv3, 1, 0.2, 3, 2, output);

    // insert a knot, with multiplicity 1, 2-times
    failed += _testTensorBoehm(kv1, kv4, kv2, 1, 0.5, 2, 2, output);


    return failed;
}



// TESTING FUNCTION
// Tests cube points - starting indices where we perfom Boehm algorithm. It
// prints out all the relavant indices
//
void mainTestCubePoints()
{
    // if you want long output, set here to true
    // testing NextCubePoint

    for (unsigned i = 0; i < 3; ++i) {
        gsInfo << "\n\n"
                  << "===================================================\n"
                  << "Testing if indices are in right order for direction "
                  << i << ": " << "\n";

        _testNextCubePoint(i);
    }
}



// TESTING FUNCTION
// tests gsTensorBoehmRefine function for 3D and 2D case
//
// \param kv1 - knot vector in x direction
// \param kv2 - knot vector in y direction
// \param kv3 - knot vector in z direction
// \param direction - direction in which we will refine knots
// \param insert_knots - vector of knots we will insert
// \param d - dimension (should work for 2 and 3, and not for other dimensions
// \param output - if we print anything on gsInfo
// \param change_knots - if we want to update knots
//
// \return number of failed tests

// It is pretty ugly function. If you are testing for dimension 2 (d = 2) you
// must provide a dummy knot vector 3 (kv3) for third dimension (that' because
// I copied from _testTensorBoehm function). Sorry.
// This function (and _testTensorBoehm) could be written better, but because
// it is a testing function I did not do it.
unsigned _testTensorBoehmRefine(gsKnotVector<> kv1,
                                gsKnotVector<> kv2,
                                gsKnotVector<> kv3,
                                int direction,
                                std::vector<real_t> insert_knots,
                                int d = 3,
                                bool output = true,
                                bool change_knots = false)
{
    GISMO_ASSERT(direction < d,
                 "Can not insert into non existing direction.");

    // number of points in certain direction
    std::vector<unsigned> dim_of_coeff(d);
    std::vector<gsKnotVector<> > knots(d);

    knots[0] = kv1;
    knots[1] = kv2;

    if (2 < d)
        knots[2] = kv3;

    // building dimension of the coefficient (number of points) and the strides
    dim_of_coeff[0] = kv1.size() - kv1.degree() - 1;
    dim_of_coeff[1] = kv2.size() - kv2.degree() - 1;

    if (2 < d)
        dim_of_coeff[2] = kv3.size() - kv3.degree() - 1;

    std::vector<unsigned> dim_of_coeff2(dim_of_coeff);
    dim_of_coeff2[direction] += insert_knots.size();

    gsVector<unsigned> strides = buildStrides(dim_of_coeff);
    gsVector<unsigned> strides2 = buildStrides(dim_of_coeff2);

    if (output)
    {
        gsInfo << "\n\n=============================================\n"
                  << "FUNCTION: _testTensorBoehmRefine:\n"
                  << "Knot vectors:" << "\n"
                  << kv1.detail() << "\n"
                  << kv2.detail() << "\n"
                  << kv3.detail() << "\n"
                  << "direction: " << direction
                  << "\nknots: ";

        _print_std_vector<>(insert_knots);

        gsInfo << "number of knots: " << insert_knots.size()
                  << "\ndimension: " << d << "\n";
    }

    // control points
    gsMatrix<> coefs;
    if (d == 3)
    {
        coefs = _buildTestCoeffs(dim_of_coeff[0],
                                 dim_of_coeff[1],
                                 dim_of_coeff[2]);
    }
    else if (d == 2)
    {
        coefs = _buildTestCoeffs2D(dim_of_coeff[0], dim_of_coeff[1]);
    }
    else
    {
        gsInfo<<"Function only works for d = 2 or d = 3.\n";
        return 1;
    }

    gsMatrix<> coefs_copy(coefs);

    // number of failed insertions
    unsigned failed = 0;

    if (output)
    {
        gsInfo << "Dimension of coefficients"
                  << "(number of points in directions): ";
        _print_std_vector<unsigned>(dim_of_coeff);
        gsInfo << "Strides: ";
        gsInfo << strides.transpose();
        gsInfo << "\n\nCalling gsTensorBoehmRefine function.\n" << "\n";
    }

    std::vector<real_t>::const_iterator itBegin = insert_knots.begin();
    std::vector<real_t>::const_iterator itEnd = insert_knots.end();
    // copy knots in case we will change it in the gsTensorBoehmRefine function
    gsKnotVector<> copy_knots = knots[direction];
    gsTensorBoehmRefine(knots[direction], coefs, direction, strides,
                        itBegin, itEnd, change_knots);

    if (output && change_knots)
    {
        gsInfo << "You set variable change_knots on true.  "
                  << "Here are changed knots: \n"
                  << knots[direction].detail() << "\n";
    }


    // number of control points in given direction
    unsigned t = dim_of_coeff[direction];

    // indices
    gsVector<int> index2(d);
    for (int i = 0; i < d; i++)
        index2[i] = 0;

    gsVector<int> index(index2);

    // necessary for computation of the indices
    gsVector<int> first(index2);
    gsVector<int> last(d);
    gsVector<int> last2(d);
    getLastIndex(strides2, coefs.rows(), last2);
    getLastIndex(strides, coefs_copy.rows(), last);
    last2[direction] = 0;
    last[direction] = 0;

    do
    {
        if (output)
        {
            gsInfo << "\n\n\n-----------------------------" << "\n";
            gsInfo << "-----------------------------" << "\n";
            gsInfo << "-----------------------------" << "\n";
        }

        int ind = getIndex(strides, index);
        gsMatrix<> coefs1(t, d);

        // copy proper coefficients to do 1D Boehm refine algorithm
        for (unsigned i = 0; i < t; ++i)
        {
            coefs1.row(i) = coefs_copy.row(ind + i * strides[direction]);

        }

        // 1D algorithm
        gsBoehmRefine(
                    copy_knots,
                    coefs1,
                    copy_knots.degree(),
                    insert_knots.begin(),
                    insert_knots.end(),
                    false);

        if (output)
        {
            gsInfo << "gsBoehmRefine: \n"
                      << coefs1 << "\n";
        }

        // collect proper coefficients for comparison with gsBoehmRefine (1D)
        int ind2 = getIndex(strides2, index2);
        nextCubePoint(index2, first, last2);

        gsMatrix<> temporary_mat(dim_of_coeff2[direction], d);

        for (unsigned i = 0; i < dim_of_coeff2[direction]; ++i)
            temporary_mat.row(i) = coefs.row(ind2 + i * strides2[direction]);


        // looking at the norm of difference betwen two coefficients matrices
        gsMatrix<> difference = coefs1 - temporary_mat;
        bool tmp_bol = difference.norm() < 0.00000000001;

        if (output)
        {
            if (tmp_bol)
            {
                gsInfo << "\nSUCCESS, norm (difference "
                          << "between gsBoehmRefine and geTensorBoehmRefine): "
                          << difference.norm() << "\n";
            }
            else
            {
                gsInfo << "\n----------------------------"
                          << "Coefficients gsTensorBoehmRefine:\n"
                          << temporary_mat
                          << "\nAre equal: " << tmp_bol
                          << "\nDifference norm: " << difference.norm()
                          << "\n" << difference << "\n";
            }
        }

        if (!tmp_bol)
            failed++;

    } while(nextCubePoint(index, first, last));

    if (output)
    {
        gsInfo << "\nFailed tries: " << failed << "\n" << "\n";
    }

    return failed;
}



// TESTING FUNCTION
// Tests tensor Boehm refine algorithm in several cases for 2 and 3D.
//
// \param output - if we want any output
//
// \return number of failed tests
unsigned mainTestTensorBoehmRefine(bool output)
{

    unsigned failed = 0; // number of failed tests

    // knot vectors
    gsKnotVector<> kv1(0, 1, 0, 2);
    gsKnotVector<> kv2(0, 1, 1, 3);
    gsKnotVector<> kv3(0, 1, 2, 4);
    gsKnotVector<> kv4(0, 1, 1, 4);

    // knots we insert
    std::vector<real_t> vec1(1, 0.5);
    std::vector<real_t> vec2(2);
    std::vector<real_t> vec3(10);

    for (int i = 0; i < 10; i++)
        vec3[i] = i * 1.0 / 10 + 0.001;



    // -------------------------------------------------------------------------
    //  3D tests
    // -------------------------------------------------------------------------


    // inserting one knot

    // here we get multiple knots
    failed += _testTensorBoehmRefine(kv1, kv2, kv3, 1, vec1, 3, output);

    failed += _testTensorBoehmRefine(kv1, kv2, kv3, 0, vec1, 3, output);
    failed += _testTensorBoehmRefine(kv1, kv2, kv3, 2, vec1, 3, output);

    vec1[0]= 0.1;
    failed += _testTensorBoehmRefine(kv1, kv2, kv3, 1, vec1, 3, output);


    // inserting two knots
    vec2[0] = 0.3;
    vec2[1] = 0.6;

    failed += _testTensorBoehmRefine(kv1, kv2, kv3, 0, vec2, 3, output);
    failed += _testTensorBoehmRefine(kv1, kv2, kv3, 1, vec2, 3, output);
    failed += _testTensorBoehmRefine(kv1, kv2, kv3, 2, vec2, 3, output);

    // insert two equal knots
    vec2[0] = 0.5;
    vec2[1] = 0.5;
    failed += _testTensorBoehmRefine(kv1, kv2, kv3, 2, vec2, 3, output);


    // inserting 10 knots
    failed += _testTensorBoehmRefine(kv1, kv2, kv3, 0, vec3, 3, output);
    failed += _testTensorBoehmRefine(kv1, kv2, kv3, 1, vec3, 3, output);
    failed += _testTensorBoehmRefine(kv1, kv2, kv4, 2, vec3, 3, output);



    // -------------------------------------------------------------------------
    // 2D tests
    // -------------------------------------------------------------------------


    // inserting one knot

    failed += _testTensorBoehmRefine(kv4, kv2, kv3, 0, vec1, 2, output);
    failed += _testTensorBoehmRefine(kv4, kv2, kv3, 1, vec1, 2, output);


    // inserting two knots
    vec2[0] = 0.5;
    vec2[1] = 0.7;
    failed += _testTensorBoehmRefine(kv4, kv2, kv3, 0, vec2, 2, output);

    vec2[0] = 0.5;
    vec2[1] = 0.7;
    failed += _testTensorBoehmRefine(kv4, kv2, kv3, 1, vec2, 2, output);


    // inserting 10 knots
    failed += _testTensorBoehmRefine(kv4, kv3, kv3, 0, vec3, 2, output);
    failed += _testTensorBoehmRefine(kv4, kv3, kv3, 1, vec3, 2, output);

    // if you want to see if the knots are updated properly
    // output must be set on true, to see the changes
    // failed += _testTensorBoehmRefine(kv4, kv3, kv3, 1, vec3, 2, true, true);


    return failed;
}



int main() {

    //mainTestCubePoints();



    unsigned failed1 = mainTestTensorBoehm(false);

    gsInfo << "\nNumber of failed tests for Tensor Boehm algorithm: "
              << failed1 << "\n" << "\n";



    unsigned failed2 = mainTestTensorBoehmRefine(false);

    gsInfo << "Number of failed tests for Tensor Boehm Refine algoritm: "
              << failed2 << "\n\n"
              << "Number of all failed tests: "
              << failed1 + failed2 << "\n"
              << "\n";



    return failed1 + failed2;
}

