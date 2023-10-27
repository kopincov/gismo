
#include "otherUtils.h"
#include <gsUtils/gsPointGrid.h>
#include <gsIO/gsParaviewCollection.h>
#include <gsCore/gsFuncData.h>
#include <gsCore/gsGeometry.h>

namespace gismo {

template <typename T,typename TT>
void expandTo3 (TT &M, T val)
{
    const index_t oldR=M.rows();
    if (oldR<3)
    {
        M.conservativeResize(3,M.cols());
        M.bottomRows(3-oldR).setConstant(val);
    }

}

void writeVTSdata(std::ofstream &file, const std::string &name, const gsMatrix<real_t> &data)
{
    file<<"<DataArray type=\"Float32\" ";
    if (name.length()>0)
        file<<"Name=\""<<name<<"\" format=\"ascii\" ";
    file<<" NumberOfComponents=\""<< data.rows() <<"\">\n";
    for ( index_t j=0; j<data.cols(); ++j)
        for ( index_t i=0; i<data.rows(); ++i)
            file<< data(i,j)<<" ";
    file<<"</DataArray>\n";
}

void writeVTSpreamble(std::ofstream &file, const gsMatrix<> &pts, const gsVector<unsigned> &np)
{
    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
    file <<"<StructuredGrid WholeExtent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
    file <<"<Piece Extent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
    file <<"<Points>\n";

    writeVTSdata(file, "", pts);

    file <<"</Points>\n";
    file <<"<PointData Scalars=\"pressure\" Vectors=\"velocity\">\n";
}

void writeVTSpostamble(std::ofstream &file)
{
    file <<"</PointData>\n";
    file <<"</Piece>\n";
    file <<"</StructuredGrid>\n";
    file <<"</VTKFile>\n";
}


void printPatchFunctions(const gsRemappedBasis &space, const std::string &filename, size_t numPoints, size_t patch)
{
    const size_t size=space.size();

    gsVector<real_t> lower=space.getSelector().patch(patch).getBoundingBox().col(0);
    gsVector<real_t> upper=space.getSelector().patch(patch).getBoundingBox().col(1);

    gsVector<unsigned> np = uniformSampleCount(lower,upper, numPoints );
    gsMatrix<> pts=gsPointGrid(lower,upper,np);
    gsMatrix<> ptsW=pts;
    gsMatrix<> sum;
    expandTo3(ptsW, static_cast<real_t>(0));
    expandTo3(np,static_cast<unsigned>(1));
    gsFuncData<> val(NEED_VALUE);
    val.patchId=patch;
    gsMatrix<> coef;
    coef.setZero(size,1);
    sum.setZero (1,pts.cols());

    std::ofstream file(filename.c_str());


    const int numDigits=static_cast<int>(math::ceil(math::log10(static_cast<double>(size))));
    writeVTSpreamble(file,ptsW,np);

    for (size_t b=0; b<size;++b)
    {
        coef(b)=1;
        gsRemappedBasis *phi=space.makeReduced(coef,false);
        phi->compute(pts,val);
        delete phi;
        sum+=val.values[0];
        std::stringstream field;
        field<<"basis_"<<std::setfill('0')<<std::setw(numDigits)<<b;
        writeVTSdata(file, field.str().c_str(), val.values[0]);
        coef(b)=0;
    }
    writeVTSdata(file, "sum", sum);

    writeVTSpostamble(file);
    file.close();
}


void printAllFunctions(const gsRemappedBasis &space, const std::string &filename, size_t numPoints)
{
    size_t numPatches=space.getSelector().patchNum();
    if (numPatches==1)
    {
        printPatchFunctions(space,filename+".vts",numPoints,0);
    }
    else
    {
        gsParaviewCollection collection(filename);
        for (size_t p=0;p<numPatches;++p)
        {
            std::stringstream patchFile;
            patchFile<<filename<<"_patch_"<<p<<".vts";
            printPatchFunctions(space,patchFile.str(),numPoints,p);
            collection.addPart(patchFile.str());
       }
       collection.save();
    }
}



} // gismo
