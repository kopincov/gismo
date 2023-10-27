/* Tests tensor-product B-splines that are periodic in one of the coordinate directions.
 *
 * Author: Dominik Mokris, dominik.mokris@jku.at
 *
 *
*/

#include <iostream>
#include <string>
#include <vector>

#include <gismo.h>


using namespace gismo;

int callParaviewOnExit(const gsGeometry<>& geometry,
                        std::string const & filename,
                        unsigned numEvaluationPoints = 100)
{
    gsWriteParaview( geometry, filename, numEvaluationPoints );
    char cmd1[100];
    strcpy(cmd1,("paraview " + filename + ".vts\0 &").c_str() );
    return system(cmd1);
}

int main(int argc, char *argv[])
{
    bool plot = false; // If set to true, paraview file is generated and launched on exit

    int test_result = 0; // Number of failed tests in this file.

    gsKnotVector<> knotVector0(0,1,1,3);//start,end,interior knots, start/end multiplicites of knots
    gsKnotVector<> knotVector1(0,1,1,3);
    gsTensorBSplineBasis<2,real_t> tensorBasis (new gsBSplineBasis<>(knotVector0), new gsBSplineBasis<>(knotVector1));
    gsMatrix<> uValues(2,6);
    uValues << 0,0.1,1,0.5,0.6,1, 0,0,0,0.5,0.6,1;
    gsMatrix<index_t> activeResults;
    gsMatrix<> valueResults;
    tensorBasis.active_into(uValues,activeResults);
    tensorBasis.eval_into(uValues,valueResults);

    gsMatrix<> coefs(16,3);
    coefs <<    1,1,0, 1,2,1, 1,3,0, 1,4,1,
                2,1,1, 2,2,0, 2,3,1, 2,4,0,
                3,1,0, 3,2,1, 3,3,0, 3,4,1,
                4,1,1, 4,2,0, 4,3,1, 4,4,0;
    gsTensorBSpline<2, real_t> surface(tensorBasis, coefs);

    surface.setPeriodic(0);
    gsInfo << "coefs after periodification:\n" << surface.coefs() << "\n";

    if( plot )
        gsInfo << "paraview returned: " << callParaviewOnExit(surface, "PeriodicSurface", 100) << "\n";

    if( test_result > 0 )
    {
        gsWarn << test_result << " failed tests. Please, try and correct the code.\n";
        return 1;
    }
    else
        return 0;
}
