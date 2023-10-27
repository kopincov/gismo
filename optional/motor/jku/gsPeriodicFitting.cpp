#include <iostream>
#include <string>
#include <math.h>

#include <gismo.h>

using namespace gismo;


void getPoints(const int n,
               gsMatrix<>& pars,
               gsMatrix<>& pts,
               const real_t a,
               const real_t b)
{
    pars.resize(1, n);
    pts.resize(2, n);
    for (int i = 0; i != n; i++)
    {
        real_t phi = i * 1.0 / (n - 1);

        pars(0, i) = phi;
        pts(0, i) = a * math::cos(phi * 2 * EIGEN_PI);
        pts(1, i) = b * math::sin(phi * 2 * EIGEN_PI);
    }
}

std::vector<real_t> getPeriodicKnots(const int degree,
                                     const int numKnots)
{
    real_t h = 1.0 / numKnots;

    std::vector<real_t> knots;

    for (int d = degree; 0 < d; d--)
    {
        knots.push_back(-1 * d * h);
    }

    for (int i = 0; i != numKnots + 1; i++)
    {
        knots.push_back(i * h);
    }

    for (int d = 1; d != degree + 1; d++)
    {
        knots.push_back(1 + d * h);
    }

    return knots;
}

gsKnotVector<> getPeriodicKnotVector(const int degree,
                                     const int numKnots)
{
    std::vector<real_t> knots = getPeriodicKnots(degree, numKnots);
    return gsKnotVector<>(knots, degree);
}

void print(const real_t& el)
{
    std::cout << el << " ";
}

int main(int argc, char *argv[])
{    

    int n = 20;
    int degree = 2;
    int numKnots = 15;
    real_t a = 2;
    real_t b = 1;
    
    gsMatrix<> pars;
    gsMatrix<> pts;

    getPoints(n, pars, pts, a, b);

    gsWriteParaviewPoints(pts, "fittingPoints");

    gsKnotVector<> kv = getPeriodicKnotVector(degree, numKnots);
    
    gsCurveFitting<> fitting(pars.transpose(), pts.transpose(), kv, true);
    fitting.compute_periodic();

    gsBSpline<> curve = fitting.curve();

    gsWriteParaview(curve, "fittingResult");
    
    real_t error;
    fitting.computeApproxError(error);

    std::cout << "Error: " << error << std::endl;
    
    return 0;
}
