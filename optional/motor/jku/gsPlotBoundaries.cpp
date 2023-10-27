#include <fstream>
#include <sstream>

#include <gismo.h>
#include <gismo_dev.h>

//#include <gsModeling/gsMultiFitting.h>
//#include <gsModeling/gsAdaptiveMultiFitting.h>

#include <gsIO/gsIOUtils.h>

using namespace gismo;

gsMultiPatch<> extractBoundaries(gsMultiPatch<>& mp)
{
    mp.computeTopology(1e-4);
    gsMultiPatch<> result;

    for (gsMultiPatch<>::const_biterator it = mp.bBegin();
         it != mp.bEnd(); ++it)
    {
        const gsGeometry<> & g = mp.patch( it->patch );
        result.addPatch( g.boundary(it->side()) );
    }

    result.computeTopology(1e-4);

    return result;
}

void addLinesToMesh(gsMesh<>& mesh,
                    gsGeometry<>& geom,
                    const gsMatrix<real_t>& samples)
{
    typedef gsMesh<real_t>::VertexHandle Vertex;

    gsMatrix<> result;
    geom.eval_into(samples, result);

    gsVector<real_t> vertex(result.rows());
    vertex = result.col(0);
    mesh.addVertex(vertex);
    for (int col = 1; col != result.cols(); ++col)
    {
        Vertex first = mesh.vertices().back();
        vertex = result.col(col);
        mesh.addVertex(vertex);
        Vertex second = mesh.vertices().back();
        mesh.addEdge(first, second);
    }
}


void addLines(gsMultiPatch<>& mp,
              gsMesh<>& mesh,
              const int resolution)
{
    for (unsigned i = 0; i != mp.nPatches(); i++)
    {
        gsGeometry<>& geom = mp.patch(i);
        const gsMatrix<real_t>& support  = geom.support();
        const gsVector<real_t>& supStart = support.col(0);
        const gsVector<real_t>& supEnd = support.col(1);

        const gsMatrix<real_t>& samples = gsPointGrid(supStart(0), supEnd(0), resolution);

        addLinesToMesh(mesh, geom, samples);
    }
}


void write(gsMultiPatch<>& mp,
           const std::string output,
           const int resolution)
{
    gsMesh<> mesh;

    addLines(mp, mesh, resolution);
    
    gsWriteParaview(mesh, output);
}


int main(int argc, char* argv[])
{
    std::string multipatchFile("");
    std::string output("");
    int resolution = 20;

    gsCmdLine cmd("Plotting boundaries");

    cmd.addString("m", "m", "Multipatch file", multipatchFile);
    cmd.addString("o", "out", "Output file", output);
    cmd.addInt("r", "r", "Resolution", resolution);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "---------------------------------------------------------\n\n"
           << "Input Arguments: \n\n"
           << "Multipatch file: " << multipatchFile << "\n\n"
           << "Output: " << output << "\n\n"
           << "Resolution: " << resolution << "\n\n"
           << "---------------------------------------------------------\n"
           << std::endl;

    gsFileData<> fd(multipatchFile);
    gsMultiPatch<>* mp = fd.getAnyFirst< gsMultiPatch<> >().release();

    gsMultiPatch<> boundaries = extractBoundaries(*mp);

    write(boundaries, output, resolution);

    return 0;
}

        
    
    
