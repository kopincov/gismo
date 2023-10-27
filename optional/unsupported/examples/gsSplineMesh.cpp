
#include <gismo.h>
#include <gismo_dev.h>

#include <iostream>



using namespace gismo;

int main(int argc, char *argv[])
{
    std::string fn = "saddle.xml";
    index_t numPoints = 1000;
    
    gsCmdLine cmd("Produce points on a parameterized geometry.");
    cmd.addPlainString("filename", "File containing the input geometry", fn);
    cmd.addInt("s","samples", 
               "Number of samples (approximately)", numPoints);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Read in the input file
    gsGeometry<>::uPtr geo = gsReadFile<>(fn);

    gsMatrix<> ab = geo->parameterRange() ;
    gsVector<> a = ab.col(0);
    gsVector<> b = ab.col(1);
    gsVector<unsigned> numPointsCoordWise = uniformSampleCount(a,b, numPoints );

    // Parameter values
    gsMatrix<> uv = gsPointGrid(a,b,numPointsCoordWise) ;
    // Corresponding points on the geometry
    gsMatrix<> xyz = geo->eval  ( uv ) ;

    gsInfo<<"Writing the point-data to XML file \"spline_points.xml\"..\n";
    gsFileData<> fd;

    // Write parameter values
    fd<< uv;  // id 0

    // Write corresponding points on the geometry
    fd<< xyz; // id 1

    fd.dump("spline_points");

    return 0;
}
