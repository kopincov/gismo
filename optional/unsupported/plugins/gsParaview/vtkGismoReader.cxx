
#include "vtkGismoReader.h"

#include "vtkIndent.h"

#include "vtkObjectFactory.h"
#include "vtkStructuredGrid.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <gsIO/gsIOUtils.h>


vtkStandardNewMacro(vtkGismoReader);

//----------------------------------------------------------------------------
vtkGismoReader::vtkGismoReader() 
    : FileName(NULL), 
      Resolution(1000), 
      PlotType(1), 
      ControlNet(false), 
      m_patches(NULL), 
      m_isOpen(false)
{
    // for debugging
    // gsDebug << "vtkGismoReader::Constructor" << std::endl;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkGismoReader::~vtkGismoReader()
{
    // for debugging
    // gsDebug << "vtkGismoReader::Destructor\n" << std::endl;
    this->SetFileName(0);
    if (m_patches) 
        delete m_patches;
}

//----------------------------------------------------------------------------
void vtkGismoReader::PrintSelf(ostream& os, vtkIndent indent)
{
    gsDebug << "vtkGismoReader::PrintSelf" << std::endl;
    this->Superclass::PrintSelf(os, indent);
    os << indent << "FileName: "
       << (this->FileName ? this->FileName : "(NULL)") << endl;

    os << indent << "G+Smo " << m_patches << endl;
}

//-----------------------------------------------------------------------------
int vtkGismoReader::CanReadFile(const char *fname)
{
    // for debugging
    // gsDebug << "vtkGismoReader::CanReadFile" << std::endl;
    //gsDebugVar(fname);
    
    using namespace gismo;

    std::string file(fname);
    gsFileData<> fileData(file);
    if (fileData.has< gsGeometry<> >() || fileData.has< gsMultiPatch<> >())
    {
        return 1;
    }
    else 
    {
        return 0;
    }
}

//----------------------------------------------------------------------------
vtkStructuredGrid* vtkGismoReader::GetOutput()
{
    // for debugging
    // gsDebug << "vtkGismoReader::GetOutput" << std::endl;
    return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkStructuredGrid* vtkGismoReader::GetOutput(int idx)
{
    // for debugging
    // gsDebug << "vtkGismoReader::GetOutput(int idx)" << std::endl;
    return vtkStructuredGrid::SafeDownCast( this->GetOutputDataObject(idx) );
}


//----------------------------------------------------------------------------
const char* vtkGismoReader::GetDataSetName()
{
    // for debugging
    // gsDebug << "vtkGismoReader::GetDataSetName" << std::endl;
    return "GismoData";
}

//-----------------------------------------------------------------------------

void vtkGismoReader::EvaluateGeometry(const gismo::gsGeometry<>& geom,
                                      const int numSamples,
                                      gismo::gsVector<unsigned>& numPoints,
                                      gismo::gsMatrix<>& points,
                                      gismo::gsMatrix<>& normals) const 
{
    using namespace gismo;
    
    const unsigned geoDim = geom.geoDim();
    const unsigned parDim = geom.parDim();
    
    gsMatrix<> support = geom.parameterRange();
    gsVector<> suppStart = support.col(0);
    gsVector<> suppEnd = support.col(1);
    
    numPoints = uniformSampleCount(suppStart, suppEnd, numSamples); 
    gsMatrix<> params = gsPointGrid(suppStart, suppEnd, numPoints);

    geom.eval_into(params, points);
    
    if (parDim == 2)
    {
        gsMapData<> fd;
        fd.flags = NEED_NORMAL;
        fd.points = params;
        geom.computeMap(fd);

        gsVector<> n = fd.normal(0);
    
        normals.resize(n.rows(), params.cols());

        for (index_t i = 0; i != params.cols(); i++)
        {
            n = fd.normal(i);
            n /= n.norm();
            normals.col(i) = n;
        }
    }
    //
    // The following code fails because there is no implementation of 
    // outerNormal_into member function for geometries. 
    // TO DO: Uncomment the code when outerNormal_into is implemented.
    // 
    // else if (parDim == 3)
    // {
    //     gsVector<> zero(3);
    //     zero.setZero();
        
    //     gsMatrix<> outerNormal;
    //     gsVector<> normal_tmp;
        
    //     for (index_t i = 0; i != params.cols(); i++)
    //     {
    //         const gsMatrix<>& param = params.col(i);

    //         std::cout << "i: " << i << std::endl;

    //         if (isOnBoundary(param, support))
    //         {
    //             std::cout << "inside if" << std::endl;
    //             geom.outerNormal_into(param, outerNormal);
    //             normal_tmp = outerNormal.col(0);
    //             normal_tmp /= normal_tmp.norm();
    //             normals.col(i) = normal_tmp;
    //         }
    //         else
    //         {
    //             std::cout << "inside else" << std::endl;
    //             normals.col(i) = zero;
    //         }
    //     }
    // }
    

    if (3 - parDim > 0)
    {
        numPoints.conservativeResize(3);
        numPoints.bottomRows(3 - parDim).setOnes();
    }


    if (3 - geoDim > 0)
    {
        points.conservativeResize(3, points.cols());
        points.bottomRows(3 - geoDim).setZero();
    }
}                                      
                                      
bool vtkGismoReader::isOnBoundary(const gismo::gsMatrix<>& param,
                                  const gismo::gsMatrix<>& support)
{
    // assumption: param is one column matrix
    
    using namespace gismo;
    
    // std::cout << "param: \n"
    //           << param << "\n"
    //           << "support: \n"
    //           << support << std::endl;
    
    for (index_t row = 0; row != param.rows(); row++)
    {
        if (param(row, 0) == support(row, 0) || param(row, 0) == support(row, 1))
        {
            return true;
        }
    }
    return false;
}


//-----------------------------------------------------------------------------
int vtkGismoReader::RequestData(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector)
{
    // for debugging
    // gsDebug << "vtkGismoReader::RequestData" << std::endl;
    // gsDebugVar(PlotType);
    // gsDebugVar(ControlNet);
    
    // open file and read multipatch only if it is neccessary
    if (!m_isOpen)
    {
        m_patches = new gismo::gsMultiPatch<>();
        gismo::gsReadFile<>(FileName, *m_patches);
        m_isOpen = true;
    }

    // get the info object
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
  
    // get the ouptut
    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // data where we store the output
    vtkMultiBlockDataSet* vtkData = vtkMultiBlockDataSet::New();

    bool plotGeometry = false;
    bool plotMesh = false;
    
    if (PlotType == PlotGeometry || PlotType == PlotGeometryAndMesh)
    {
        plotGeometry = true;
    }
    
    if (PlotType == PlotMesh || PlotType == PlotGeometryAndMesh)
    {
        plotMesh = true;
    }

    for (size_t i = 0; i != m_patches->nPatches(); i++)
    {
        if (plotGeometry)
        {
            getGeometryData(i, vtkData);
        }

        if (plotMesh)
        {
            getMeshData(i, vtkData);
        }
        
        if (ControlNet)
        {
            getControlNetData(i, vtkData);
        }

    }
    
    output->ShallowCopy(vtkData);
    vtkData->Delete();
    return 1;


}

//-----------------------------------------------------------------------------
void vtkGismoReader::getControlNetData(const size_t patch,
                                       vtkMultiBlockDataSet* vtkData) const
{
    gismo::gsMesh<> mesh;
    m_patches->patch(patch).controlNet(mesh);
    getMeshData(mesh, vtkData);
}




//-----------------------------------------------------------------------------
void vtkGismoReader::getMeshData(const size_t patch,
                                 vtkMultiBlockDataSet* vtkData) const
{
    gismo::gsMesh<> mesh;
    gismo::makeMesh<>(m_patches->patch(patch).basis(), mesh, 10);
    m_patches->patch(patch).evaluateMesh(mesh);
    
    getMeshData(mesh, vtkData);
}


//-----------------------------------------------------------------------------
void vtkGismoReader::getMeshData(const gismo::gsMesh<>& mesh,
                                 vtkMultiBlockDataSet* vtkData) const
{
    // points
    vtkSmartPointer< vtkPoints > vtkPts = vtkSmartPointer<vtkPoints>::New();    
    for (size_t i = 0; i != mesh.numVertices(); i++)
    {
        vtkPts->InsertPoint(i, 
                            mesh.vertex(i).x(), 
                            mesh.vertex(i).y(),
                            mesh.vertex(i).z());
    }
    
    // lines
    vtkSmartPointer< vtkCellArray > vtkLines = vtkSmartPointer< vtkCellArray >::New();
    for (size_t i = 0; i != mesh.numEdges(); i++)
    {
        vtkLines->InsertNextCell(2);
        vtkLines->InsertCellPoint(mesh.edges()[i].source->getId());
        vtkLines->InsertCellPoint(mesh.edges()[i].target->getId());
    }
    
    // cells
    vtkSmartPointer< vtkCellArray > vtkCells = vtkSmartPointer< vtkCellArray >::New();
    for (size_t i = 0; i != mesh.numFaces(); i++)
    {
        vtkCells->InsertNextCell(mesh.faces()[i]->vertices.size());
        
        for (size_t j = 0; j != mesh.faces()[i]->vertices.size(); j++)
        {
            vtkCells->InsertCellPoint(mesh.faces()[i]->vertices[j]->getId());
        }
    }

    vtkSmartPointer< vtkPolyData > vtkMesh = vtkSmartPointer< vtkPolyData >::New();
    vtkMesh->SetPoints(vtkPts);
    vtkMesh->SetLines(vtkLines);               
    vtkMesh->SetPolys(vtkCells);
    
    vtkData->SetBlock(vtkData->GetNumberOfBlocks(), vtkMesh);
}


//-----------------------------------------------------------------------------
void vtkGismoReader::getGeometryData(const size_t patch,
                                     vtkMultiBlockDataSet* vtkData) const
{
    using namespace gismo;
    gsVector<unsigned> numPts; // number of points in each direction
    gsMatrix<> pts;
    gsMatrix<> normals;
    
    EvaluateGeometry(m_patches->patch(patch), Resolution, numPts, pts, normals);

    // make points
    vtkSmartPointer< vtkPoints > vtkPts = vtkSmartPointer<vtkPoints>::New();    
        
    for (index_t col = 0; col != pts.cols(); col++)
    {
        vtkPts->InsertPoint(col, pts(0, col), pts(1, col), pts(2, col));
    }

    // put the points into grid
    vtkSmartPointer< vtkStructuredGrid > vtkGrid = vtkSmartPointer< vtkStructuredGrid >::New();

    vtkGrid->SetDimensions(numPts[0], numPts[1], numPts[2]);
    vtkGrid->SetPoints(vtkPts);  

    if (normals.size() != 0) // if there are normals add normals
    {
        vtkDoubleArray* vtkNormals = vtkDoubleArray::New();
        vtkNormals->SetNumberOfComponents(normals.rows());
        vtkNormals->SetNumberOfTuples(normals.cols());
        vtkNormals->SetName("Normals");
    
        for (index_t col = 0; col != normals.cols(); col++)
        {
            for (index_t row = 0; row != normals.rows(); row++)
            {
                vtkNormals->SetComponent(col, row, normals(row, col));
            }
        }
    


        vtkGrid->GetPointData()->SetNormals(vtkNormals);
        vtkNormals->Delete();
    }

    // std::cout << "\n\nInfo about points: \n\n" << std::endl;
    // vtkPts->PrintSelf(std::cout, vtkIndent(0));

    // std::cout << "\n\nInfo about structured grid: \n\n" << std::endl;
    // vtkGrid->PrintSelf(std::cout, vtkIndent(0));

    vtkData->SetBlock(vtkData->GetNumberOfBlocks(), vtkGrid);

}



//-----------------------------------------------------------------------------
/*
int vtkGismoReader::RequestInformation(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector)
{
    gsDebug << "vtkGismoReader::RequestInformation \n";
    return 1;
}
*/
