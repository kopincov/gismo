// Author Jaka Speh, Angelos Mantzaflaris
//
// Plots the elements (mesh) of a geometry or of a basis.
//
// Usage:
// $ gsPlotMesh inputXmlFile [-o outputParaviewFile=mesh] [-p] [-H] [-s numberOfSamples=10]
//
// if the input is geometry, physical mesh is the output
// if the input is geometry and -p switch is on, the computational mesh (basis
// mesh) is the output
// if the input is geometry you can set number of samples per element side
// for viewing (you don't get a line, but a linear interpolant)
// if the input is THBSpline and -H switch is on, mesh for each level is written
// in its own file
//
// Examples:
//
// $ gsPlotMesh basis.xml
//
// $ gsPlotMesh geometry.xml
//
// $ gsPlotMesh geometry.xml -p
//


#include <iostream>
#include <string>

#include <gismo.h>

#include <gsIO/gsIOUtils.h>


using namespace gismo;

int main(int argc, char* argv[])
{

    std::string inName;
    std::string outName("mesh");
    bool param = false;
    index_t samples = 10;
    bool hierarchical = false;
    
    gsCmdLine cmd("Reading gsSolid");
    cmd.addPlainString("filename", "Input file", inName);
    cmd.addString("o", "outFile", "Output filename", outName);
    cmd.addSwitch("parametric", "If the input is geometry this switch ", param);
    cmd.addSwitch("hierarchical", "If the input is THB-Spline this switch writes "
                  "mesh for each leven in its own file.", hierarchical);
    cmd.addInt("s", "samples",
               "Number of samples per element side to use for viewing", samples);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (inName.empty())
    {
        gsInfo << "Please provide an input filename or use -h for help.\n";
        return 0;
    }

    gsFileData<>  filedata(inName);
    gsMesh<>      mesh;

    if (hierarchical && filedata.has<gsGeometry<> >())
    {
        gsGeometry<>::uPtr geom = filedata.getAnyFirst<gsGeometry<> >();

        std::vector<gsMesh<> > meshes;
        bool success = makeHierarchicalMesh<> (geom->basis(),
                                                     meshes, samples);

        if (success)
        {

            if (!param)
            {
                for (unsigned index = 0; index < meshes.size(); index++)
                {
                    geom->evaluateMesh(meshes[index]);
                }
            }

            gsWriteParaview(meshes, outName);

            return 0;
        }
    }


    if (filedata.has<gsBasis<> >())
    {
        gsBasis<>::uPtr basis = filedata.getAnyFirst<gsBasis<> >();

        makeMesh<>(*basis, mesh, samples);
    }
    else if (filedata.has<gsGeometry<> >())
    {
        gsGeometry<>::uPtr geom = filedata.getAnyFirst<gsGeometry<> >();

        makeMesh<>( geom->basis(), mesh, samples);

        if (!param)
        {
            geom->evaluateMesh(mesh);
        }
    }

    gsWriteParaview(mesh, outName);
    gsInfo << "Output file: " << outName << " is written." << "\n";

    return 0;
}




/// Returns new mesh. Vertices of the new mesh are
///
/// { geom(v) | v vertex of input mesh }
// template <typename TT>
// void evaluateMesh(gsMesh<TT>& mesh, const gsGeometry<TT>& geom)
// {
//     for (size_t index = 0; index < mesh.vertex.size(); index++)
//     {
//         gsMatrix<TT> result;
//         gsMatrix<TT> u = mesh.vertex[index]->coords.topRows(geom.parDim());
//         geom.eval_into(u, result);

//         int rows = result.rows();
//         result.conservativeResize(3, 1);
//         if (rows < 3)
//           result.bottomRows(3 - rows).setZero();

//         for (int i = 0; i < 3; i++)
//         {
//             mesh.vertex[index]->coords[i] = result(i, 0);
//         }
//     }
// }
