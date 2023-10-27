#include <fstream>
#include <sstream>

#include <gismo.h>
#include <gismo_dev.h>

//#include <gsModeling/gsMultiFitting.h>
//#include <gsModeling/gsAdaptiveMultiFitting.h>

#include <gsIO/gsIOUtils.h>

using namespace gismo;


real_t convertToDouble(const std::string& word)
{
    std::istringstream convert(word);
    real_t number = 0;

    if ( !(convert >> number) )
    {
        std::cout << "Error converting parameter: " << word
                  << " to double..." << std::endl;
    }

    return number;
}

void readCoordinates(const std::string& coordinatesFile,
                     std::vector< real_t >& coordsX,
                     std::vector< real_t >& coordsY)
{
    std::ifstream file;
    file.open(coordinatesFile.c_str());

    bool x = true;
    std::string word;

    if (file.is_open())
    {
        while (file >> word)
        {
            if (word == "x:")
            {
                x = true;
            }
            else if (word == "y:")
            {
                x = false;
            }
            else
            {
                const real_t param = convertToDouble(word);

                if (x)
                {
                    coordsX.push_back(param);
                }
                else
                {
                    coordsY.push_back(param);
                }
            }
        }
    }
    else
    {
        std::cout << "Error at opening file." << std::endl;
    }


}

void addLinesToMesh(gsMesh<real_t>& mesh,
                    gsGeometry<real_t>& geom,
                    const gsMatrix<real_t>& samples,
                    int dir,
                    const std::vector<real_t>& xParams)
{
    typedef gsMesh<real_t>::VertexHandle Vertex;

    gsMatrix<real_t> ones(1, samples.cols());
    ones.setOnes();

    int dir2 = (dir == 0) ? 1 : 0;

    for (std::size_t index = 0; index != xParams.size(); index++)
    {
        const real_t constant = xParams[index];
        gsMatrix<real_t> params(2, samples.cols());

        params.row(dir) = constant * ones;
        params.row(dir2) = samples;

        gsMatrix<real_t> result;
        geom.eval_into(params, result);

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
}


void addParameterLines(gsMultiPatch<>& mp,
                       gsMesh<>& mesh,
                       const std::vector< real_t >& coordsX,
                       const std::vector< real_t >& coordsY,
                       const int resolution)
{
    for (unsigned i = 0; i != mp.nPatches(); i++)
    {

        gsGeometry<>& geom = mp.patch(i);
        const gsMatrix<real_t>& support  = geom.support();
        const gsVector<real_t>& supStart = support.col(0);
        const gsVector<real_t>& supEnd = support.col(1);

        const gsMatrix<real_t>& samples = gsPointGrid(supStart(1), supEnd(1), resolution);

        addLinesToMesh(mesh, geom, samples, 0, coordsX);
        addLinesToMesh(mesh, geom, samples, 1, coordsY);
    }
}

void writeParameterLines(gsMultiPatch<>& mp,
                         const std::string& coordinatesFile,
                         const std::string& output,
                         const int resolution)
{

    std::vector< real_t > coordsX, coordsY;
    readCoordinates(coordinatesFile, coordsX, coordsY);

    gsMesh<real_t> mesh;

    addParameterLines(mp, mesh, coordsX, coordsY, resolution);

    std::cout << "Writing parameter lines to: " << output << std::endl;
    gsWriteParaview(mesh, output);
}

void writeBoundaries(gsMultiPatch<>& mp,
                     const std::string& output,
                     const int resolution)
{
    std::vector< real_t > coords;
    coords.push_back(0.0);
    coords.push_back(1.0);

    gsMesh<real_t> mesh;
    addParameterLines(mp, mesh, coords, coords, resolution);

    std::string out = output + "_Boundary";
    std::cout << "Writing boundary to: " << output << std::endl;
    gsWriteParaview(mesh, out);
}

void writeMeshes(gsMultiPatch<>& mp,
                 const std::string& output)
{
    for (unsigned i = 0; i != mp.nPatches(); i++)
    {
        gsGeometry<>& geom = mp.patch(i);

        gsMesh<> mesh;
        makeMesh<> (geom.basis(), mesh, 0);

        const std::string out = output + util::to_string(i);
        std::cout << "Writing mesh: " << out << std::endl;
        gsWriteParaview(mesh, out);
    }
}

int main(int argc, char* argv[])
{
    std::string multipatchFile("");
    std::string coordinatesFile("");
    std::string output("");
    int resolution = 20;

    gsCmdLine cmd("Fitting");

    cmd.addString("m", "m", "Multi patch file", multipatchFile);
    cmd.addString("o", "out", "Output file", output);
    cmd.addString("c", "c", "Coordianates file", coordinatesFile);
    cmd.addInt("r", "r", "Resolution", resolution);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "---------------------------------------------------------\n\n"
           << "Input Arguments: \n\n"
           << "Multipatch file: " << multipatchFile << "\n\n"
           << "Coordinates file: " << coordinatesFile << "\n\n"
           << "Output: " << output << "\n\n"
           << "Resolution: " << resolution << "\n\n"
           << "---------------------------------------------------------\n"
           << std::endl;

    gsFileData<> fd(multipatchFile);
    gsMultiPatch<>* mp = fd.getAnyFirst< gsMultiPatch<> >().release();

    writeParameterLines(*mp, coordinatesFile, output, resolution);
    writeBoundaries(*mp, output, resolution);
    writeMeshes(*mp, output);

    return 0;
}
