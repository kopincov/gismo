
#pragma once

#include <vtkMultiBlockDataSetAlgorithm.h>

#include <gismo.h>
#include <gismo_dev.h>

class vtkStructuredGrid;


/**
   @brief Paraview Reader 

   Related readings:
   - http://www.paraview.org/Wiki/Writing_ParaView_Readers
   - http://www.paraview.org/Wiki/Plugin_HowTo
   - http://www.kitware.com/media/html/WritingAParaViewReaderPlug-in.html
   - http://www.vtk.org/Wiki/ParaView/Examples/Plugins/Reader
   - http://www.paraview.org/Wiki/ParaView/Users_Guide/Loading_Plugins
   - http://www.paraview.org/Wiki/VTK/Tutorials/Composite_Datasets

*/

enum GismoPlotType
{
    PlotGeometry = 1,
    PlotMesh = 2,
    PlotGeometryAndMesh = 3
};

class VTK_EXPORT vtkGismoReader : public vtkMultiBlockDataSetAlgorithm
{
public:
    static vtkGismoReader *New();
    vtkTypeMacro(vtkGismoReader,vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    //Description:
    //Choose file to read
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    //Description:
    //Check file for suitability to this reader
    int CanReadFile(const char* fname);

    //Description:
    //Set the sampling parameter
    //void SetResolution(int);
    vtkSetClampMacro(Resolution, int, 1, VTK_INT_MAX);
    vtkGetMacro(Resolution, int);

    // Description:
    // Specify the type of plot.
    vtkSetClampMacro(PlotType, int, PlotGeometry, PlotGeometryAndMesh);
    vtkGetMacro(PlotType, int);

    // Description:
    // If true also the control net is plotted. 
    vtkSetMacro(ControlNet, bool);
    vtkBooleanMacro(ControlNet, bool);

    const char* GetDataSetName();
    vtkStructuredGrid *GetOutput();
    vtkStructuredGrid *GetOutput(int idx);

protected:

    vtkGismoReader();
    ~vtkGismoReader();

    int RequestData(vtkInformation *,
                    vtkInformationVector **,
                    vtkInformationVector *);

    // int RequestInformation(vtkInformation *,
    //                       vtkInformationVector **,
    //                       vtkInformationVector *);

protected:

    char * FileName; 
    int Resolution;
    int PlotType;
    bool ControlNet;
    
    gismo::gsMultiPatch<> * m_patches;
    bool m_isOpen;

private:
    vtkGismoReader(const vtkGismoReader&);  // Not implemented.
    void operator=(const vtkGismoReader&);  // Not implemented.

    // \brief Evaluates the geometry on uniform parameters
    void EvaluateGeometry(const gismo::gsGeometry<>& geom,
                          const int numSamples,
                          gismo::gsVector<unsigned>& numPoints,
                          gismo::gsMatrix<>& points,
                          gismo::gsMatrix<>& normals) const;

    static bool isOnBoundary(const gismo::gsMatrix<>& param,
                             const gismo::gsMatrix<>& support);

    void getGeometryData(const size_t patch, vtkMultiBlockDataSet* vtkData) const;

    void getMeshData(const size_t patch, vtkMultiBlockDataSet* vtkData) const;

    void getMeshData(const gismo::gsMesh<>& mesh, vtkMultiBlockDataSet* vtkData) const;
    
    void getControlNetData(const size_t patch, vtkMultiBlockDataSet* vtkData) const;


};

