
#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>


using namespace gismo;

int main()
{
    bool passed = true;

    gsGeometry<>::uPtr f = gsNurbsCreator<>::NurbsQuarterAnnulus();
    gsVector<> x = gsVector<>::vec(0.5, 0.5),
               y = gsVector<>::vec(0.1, 1.9);

    // For any function
    gsInfo << "-- Invert point on function.\n";
    int iter = f->newtonRaphson(y, x, true);

    gsMatrix<> fx = f->eval(x);

    gsInfo << "Result:     " << x .transpose() << "\n";
    gsInfo << "Value:      " << fx.transpose() << "\n";
    real_t res = (fx - y).norm();
    gsInfo << "Res.norm:   " <<  res << "\n";
    gsInfo << "Iterations: " << iter << "\n";

    passed = passed && (iter >= 0) && (res <= 1e-5);

    // For a gsGeometry
    gsInfo << "-- Invert point on geometry.\n";
    gsMatrix<> points(2,2);
    points.col(0) = gsVector<>::vec(0.2, 1.7);
    points.col(1) = gsVector<>::vec(2.0, 0.0);

    gsMatrix<> params;
    f->invertPoints(points, params);
    fx = f->eval(params);
    gsInfo << "Result:     " << params.asRowVector() << "\n";
    gsInfo << "Value:      " << fx    .asRowVector() << "\n";
    res = (fx - points).norm();
    gsInfo << "Res.norm:   " << res << "\n";
    passed = passed && (res <= 1e-5);

    gsInfo << "-- Invert point on a 3D curve.\n";   
    f = gsReadFile<>("curves3d/bspline3d_curve_02.xml");    
    points = gsVector<>::LinSpaced(5, 0.3, 0.6);
    points.transposeInPlace();
    gsInfo << "Original parameters: " << points << "\n";
    fx = f->eval(points);
    f->invertPoints(fx, params);
    gsInfo << "Inversed parameters: " << params << "\n";
    gsInfo << "Error:   " << (points - params).cwiseAbs() << "\n";
    res = (params - points).norm();
    passed = passed && (res <= 1e-5);

    return passed ? 0 : 1;
}
