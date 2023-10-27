#include <gismo.h>

#include <iostream>
#include <gsCore/gsGeometrySlice.h>

using namespace gismo;


// test: gsTensorBSpline.h :   slice; constructCoefsForSlice
//       gsGeometry.h :        getIsoParametricSlice
//       gsGeometrySlice.h
//       gsHTensorBasis.h:     getBoxesAlongSlice
//       gsTHBSplineBasis.h:   boundaryBasis; basisSlice
//       gsTHBSpline:          slice
//       gsHBSplineBasis.h:    boundaryBasis; basisSlice
//       gsHBSpline:           slice

// struct to store the data for the tests
struct TestData{
    gsGeometry<real_t>* sourceGeom;
    gsFunction<real_t>* targetGeom;
    int dir_fixed;
    real_t par;
    int targetDomainDim;
    std::string testDescription;

    ~TestData()
    {
        delete sourceGeom;
        delete targetGeom;
    }

    bool checkSlice(real_t eps)
    {
        if(targetDomainDim!=targetGeom->domainDim())
            return false;
        gsMatrix<real_t> result1,result2;
        gsMatrix<real_t> testMatrixSource,testMatrixTarget;
        if(targetDomainDim==0)
        {
            testMatrixSource.resize(1,1);
            testMatrixSource(0,0)=par;
            testMatrixTarget=testMatrixSource;
        }
        else
        {
            gsVector<real_t> start,end;
            start.setZero(targetDomainDim);
            end.setOnes(targetDomainDim);
            gsVector<unsigned> np;
            np.setConstant(targetDomainDim,1,10);
            testMatrixTarget=gsPointGrid<real_t>(start,end,np);
            testMatrixSource.resize(testMatrixTarget.rows()+1,testMatrixTarget.cols());
            testMatrixSource.topRows(dir_fixed)=testMatrixTarget.topRows(dir_fixed);
            testMatrixSource.row(dir_fixed)=gsVector<real_t>::Constant(testMatrixTarget.cols(),par);
            testMatrixSource.bottomRows(targetDomainDim-dir_fixed)=testMatrixTarget.bottomRows(targetDomainDim-dir_fixed);
        }
        sourceGeom->eval_into(testMatrixSource,result1);
        targetGeom->eval_into(testMatrixTarget,result2);
        return gsAllCloseAbsolute(result1,result2,eps);
    }
};

// create a new testData object for the geometry slice
TestData* createGeometrySliceTest(int dir_fixed,real_t par,
                         const gsGeometry<real_t>* sourceGeom)
{
    TestData* testData=new TestData();
    testData->sourceGeom= sourceGeom->clone().release();
    gsGeometrySlice<real_t> slice = testData->sourceGeom->getIsoParametricSlice(dir_fixed,par);
    testData->targetGeom= slice.clone().release();
    testData->dir_fixed=dir_fixed;
    testData->par=par;
    testData->targetDomainDim=sourceGeom->domainDim()-1;
    return testData;
}

// create a new testData object for the tensor bspline slice
template<int d>
TestData* createTensorBSplineSliceTest(int dir_fixed,real_t par,
                         const gsTensorBSpline<d,real_t>* sourceGeom)
{
    TestData* testData=new TestData();
    typename gsTensorBSpline<d,real_t>::uPtr clone = sourceGeom->clone();
    typename gsTensorBSpline<d,real_t>::BoundaryGeometryType target;
    clone->slice(dir_fixed,par,target);
    testData->sourceGeom=clone.release();
    testData->targetGeom= target.clone().release();
    testData->dir_fixed=dir_fixed;
    testData->par=par;
    testData->targetDomainDim=d-1;
    return testData;
}

// create a new testData object for the thbspline slice
template<int d>
TestData* createTHBSplineSliceTest(int dir_fixed,real_t par,
                         const gsTHBSpline<d,real_t>* sourceGeom)
{
    TestData* testData=new TestData();
    typename gsTHBSpline<d,real_t>::uPtr clone = sourceGeom->clone();
    typename gsTHBSpline<d,real_t>::BoundaryGeometryType target;
    clone->slice(dir_fixed,par,target);
    testData->sourceGeom=clone.release();
    testData->targetGeom= target.clone().release();
    testData->dir_fixed=dir_fixed;
    testData->par=par;
    testData->targetDomainDim=d-1;
    return testData;
}

// create a new testData object for the hbspline slice
template<int d>
TestData* createHBSplineSliceTest(int dir_fixed,real_t par,
                         const gsHBSpline<d,real_t>* sourceGeom)
{
    TestData* testData=new TestData();
    typename gsHBSpline<d,real_t>::uPtr clone = sourceGeom->clone();
    testData->sourceGeom=clone.release();
    typename gsHBSpline<d,real_t>::BoundaryGeometryType target;
    clone->slice(dir_fixed,par,target);
    testData->targetGeom = target.clone().release();
    testData->dir_fixed=dir_fixed;
    testData->par=par;
    testData->targetDomainDim=d-1;
    return testData;
}

// add test examples to the list of testData
void addTestExamplesToList(int dir,real_t par,const std::vector<gsGeometry<real_t>* >& geometryList,
                           const std::vector<std::string>& paths,
                        std::vector<TestData* >& testData)
{
    for(unsigned i = 0;i<geometryList.size();++i)
    {
        int dim=0;
        gsTensorBSpline<2,real_t>* tensorBSplineGeometry2D = dynamic_cast<gsTensorBSpline<2,real_t>*>(geometryList[i]);
        if(NULL!=tensorBSplineGeometry2D && 2>=dir)
        {
            dim=2;
            TestData* data = createTensorBSplineSliceTest<2>(dir,par,tensorBSplineGeometry2D);
            data->testDescription="Path: "+paths[i]+", TensorBSplineSlice in 2D test";
            testData.push_back(data);
        }
        gsTensorBSpline<3,real_t>* tensorBSplineGeometry3D = dynamic_cast<gsTensorBSpline<3,real_t>*>(geometryList[i]);
        if(NULL!=tensorBSplineGeometry3D && 3>=dir)
        {
            dim=3;
            TestData* data = createTensorBSplineSliceTest<3>(dir,par,tensorBSplineGeometry3D);
            data->testDescription="Path: "+paths[i]+", TensorBSplineSlice in 3D test";
            testData.push_back(data);
        }
        gsTHBSpline<2,real_t>* tHBSplineGeometry2D = dynamic_cast<gsTHBSpline<2,real_t>*>(geometryList[i]);
        if(NULL!=tHBSplineGeometry2D && 2>=dir)
        {
            dim=2;
            TestData* data = createTHBSplineSliceTest<2>(dir,par,tHBSplineGeometry2D);
            data->testDescription="Path: "+paths[i]+", THBSplineSlice in 2D test";
            testData.push_back(data);
        }
        gsTHBSpline<3,real_t>* tHBSplineGeometry3D = dynamic_cast<gsTHBSpline<3,real_t>*>(geometryList[i]);
        if(NULL!=tHBSplineGeometry3D && 3>=dir)
        {
            dim=3;
            TestData* data = createTHBSplineSliceTest<3>(dir,par,tHBSplineGeometry3D);
            data->testDescription="Path: "+paths[i]+", THBSplineSlice in 3D test";
            testData.push_back(data);
        }
        gsHBSpline<2,real_t>* hBSplineGeometry2D = dynamic_cast<gsHBSpline<2,real_t>*>(geometryList[i]);
        if(NULL!=hBSplineGeometry2D && 2>=dir)
        {
            dim=2;
            TestData* data = createHBSplineSliceTest<2>(dir,par,hBSplineGeometry2D);
            data->testDescription="Path: "+paths[i]+", HBSplineSlice in 2D test";
            testData.push_back(data);
        }
        gsHBSpline<3,real_t>* hBSplineGeometry3D = dynamic_cast<gsHBSpline<3,real_t>*>(geometryList[i]);
        if(NULL!=hBSplineGeometry3D && 3>=dir)
        {
            dim=3;
            TestData* data = createHBSplineSliceTest<3>(dir,par,hBSplineGeometry3D);
            data->testDescription="Path: "+paths[i]+", HBSplineSlice in 3D test";
            testData.push_back(data);
        }
        if(dim>=dir)
        {
            TestData* data = createGeometrySliceTest(dir,par,geometryList[i]);
            data->testDescription="Path: "+paths[i]+", GeometrySlice test";
            testData.push_back(data);
        }
    }
}

// read in the geometrys from the path-list
void addGeometrysToList(const std::vector<std::string>& paths,std::vector<gsGeometry<real_t>* >&geometryList)
{
    for(unsigned i = 0;i<paths.size();++i)
    {
        gsInfo << paths[i] << "\n" ;
        gsGeometry<>::uPtr geom = gsReadFile<>(paths[i]);
        if(geom)
            geometryList.push_back(geom.release());
    }
}

int main()
{
    std::vector<gsGeometry<real_t>* > geometryList;
    std::vector<TestData*> testData;
    std::vector<std::string> paths;

    //surfaces
    paths.push_back("surfaces/shroud_kaplan.xml"); // tensorbspline
     // hb
    paths.push_back("planar/square_thb.xml");      // thb
    //volumes
    paths.push_back("volumes/cube.xml");           // tensorbspline
    paths.push_back("volumes/GshapedVolume.xml");  // tensorbspline
    // hb
    paths.push_back("volumes/waveTHB_refined.xml");// thb

    addGeometrysToList(paths,geometryList);
    gsInfo << geometryList.size() << " Geometrys read in.\n";

    int dir=0;
    real_t par=0.5;
    gsInfo << "Creating appropriate test examples for dir=" << dir << " and par=" << par <<"... \n";
    addTestExamplesToList(dir,par,geometryList,paths,testData);
    dir=1;
    par=0.12;
    gsInfo << "Creating appropriate test examples for dir=" << dir << " and par=" << par <<"... \n";
    addTestExamplesToList(dir,par,geometryList,paths,testData);

    bool passed = true;
    real_t eps=std::numeric_limits<real_t>::epsilon()*1e5;
    gsInfo << "Run all the tests for the geometries: \n";
    for(unsigned i = 0;i<testData.size();++i)
    {
        bool passed_i = testData[i]->checkSlice(eps);
        if(passed_i)
            gsInfo<<"Testdata "<<i<<" passed.\n";
        else
            gsInfo<<"Testdata "<<i<<": "<<testData[i]->testDescription<<" failed.\n";
        passed = passed && passed_i;
    }

    gsInfo << "Freeing the memories... \n";
    freeAll(geometryList);
    freeAll(testData);

    if(passed)
        gsInfo << "All test passed!\n";
    return passed ? 0 : 1;
}
