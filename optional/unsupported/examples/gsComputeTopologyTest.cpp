/** @file gsComputeTopologyTest.cpp

    Test the topology reconstruction algorithm.
    It creates a two patch geometry starting from an interface.
    Then it calls computeTopology and checks that the computed interface
    equals the original one and that the number of boundaries is correct.

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):   A. Bressan
    Created on:  2014-11-26
*/


#include  <gismo.h>

#include  "interfaceTestUtils.h"


using namespace gismo;

bool temp_pass;
bool passed=true;
#define TEST(a)\
    temp_pass = (a); passed = (passed && temp_pass);\
    gsInfo << (temp_pass? "TEST OK\n":"TEST FAIL\n");

void printMapping (const gsMatrix<> &orig, const gsMatrix<> &mapped)
{
    gsInfo<<"F:";
    for(int i=0;i< orig.rows()-1;++i)
    {
        if (orig(i,0)==orig(i,1))
            gsInfo<<"{"<<orig(i,0)<<"}x";
        else
            gsInfo<<"["<<orig(i,0)<<","<<orig(i,1)<<"]x";
    }
    if (orig(orig.rows()-1,0)==orig(orig.rows()-1,1))
        gsInfo<<"{"<<orig(orig.rows()-1,0)<<"}->";
    else
        gsInfo<<"["<<orig(orig.rows()-1,0)<<","<<orig(orig.rows()-1,1)<<"]->";
    for(int i=0;i< mapped.rows()-1;++i)
    {
        if (mapped(i,0)==mapped(i,1))
            gsInfo<<"{"<<mapped(i,0)<<"}x";
        else
            gsInfo<<"["<<mapped(i,0)<<","<<mapped(i,1)<<"]x";
    }
    if (mapped(mapped.rows()-1,0)==mapped(mapped.rows()-1,1))
        gsInfo<<"{"<<mapped(mapped.rows()-1,0)<<"}  ";
    else
        gsInfo<<"["<<mapped(mapped.rows()-1,0)<<","<<mapped(mapped.rows()-1,1)<<"]  ";
}



gsBasis<> & getDimLinearBasis (int dim)
{
    gsKnotVector<> kv(0,1,0,2,0);
    switch(dim)
    {
    case 1:
        return *(new gsBSplineBasis<>(kv) );
    case 2:
        return *(new gsTensorBSplineBasis<2>(kv,kv) );
    case 3:
        return *(new gsTensorBSplineBasis<3>(kv,kv,kv) );
    default:
        GISMO_ERROR("we only support dimension 1,2,3");
    }
}

bool checkInterfaces(int deg ,const boundaryInterface &bi1, const boundaryInterface &bi2)
{
    const int dim = bi1.dirMap(bi1.first()).rows();
    boundaryInterface other;
    if ( bi1.first() == bi2.first() )
        other=bi2;
    else if (bi1.first() == bi2.second())
        other = bi2.getInverse();
    else return false;
    for (int i=0; i<dim; ++i)
    {
        bool ok = true;
        ok = ok && (bi1.dirMap(bi1.first())(i) == other.dirMap(other.first())(i));
        ok = ok && (bi1.dirOrientation(bi1.first())(i) == other.dirOrientation(other.first())(i));
        if (!ok && deg<=i) // skip differences on degenerate axis as the solution is not unique
            return false;
    }
    return bi1.second() == other.second();
}


bool test (int dim,int deg)
{
    bool pass=true;

    gsBasis<>     &basis    = getDimLinearBasis(dim);

    gsMatrix<> boxG=basis.support();
    gsMatrix<> box1=boxG;
    for (index_t i=0; i<deg; ++i) box1(i,1)=0;
    gsMatrix<> box2;

    gsMatrix<> corG=getCorners(boxG);
    gsMatrix<> cor1=getCorners(box1);
    gsMatrix<> cor2=corG;

    // the coefficients will be changed in the inner loop
    gsGeometry<>::uPtr geo1 = basis.makeGeometry(cor1.transpose());
    gsGeometry<>::uPtr geo2 = basis.makeGeometry(cor2.transpose());

    boundaryInterface current;
    firstInterface(current,dim);
    index_t count=0;
    do
    {
        ++count;
        gsMultiPatch<> mpr;
        mpr.addPatch(geo1->clone());
        mpr.addPatch(geo2->clone());

        // compute coefficients for patch1 to avoid that two opposite faces coincide in the degenerate case
        cor1=getCorners(box1);
        std::vector<boxCorner> notM1;
        current.first().opposite().getContainedCorners(dim,notM1);
        for (index_t i=0; i< math::min(static_cast<int>(notM1.size()), dim);++i )
            cor1(deg,notM1[i]-1)+=static_cast<real_t>(i+1)/40;
        mpr.patch(0).coefs()=cor1.transpose();

        // compute coefficients for the second patch
        gsAffineFunction<real_t>  interfaceMap=mpr.getMapForInterface(current.getInverse());
        interfaceMap.eval_into(corG,cor2);
        cor2.topRows(deg).setZero();
        // avoid that two opposite faces of box2 coincide
        std::vector<boxCorner> notM2;
        current.second().opposite().getContainedCorners(dim,notM2);
        for (index_t i=0; i< math::min(static_cast<int>(notM2.size()), dim);++i )
            cor2(deg, notM2[i]-1)-=static_cast<real_t>(i+1)/40;
        mpr.patch(1).coefs()=cor2.transpose();

        mpr.computeTopology(10*math::limits::epsilon());
        std::vector<boundaryInterface> inter=mpr.interfaces();
        std::vector<patchSide>         bound=mpr.boundaries();
        if (inter.size()!=1                                                  // only two faces coincide
                || !checkInterfaces(deg,current,inter[0])                        // the faces and orientation are correct
                || static_cast<size_t>(2*dim*2) != 2*inter.size()+bound.size())  // all other faces are boundaries
        {
            box1=(mpr.patch(0).eval(corG)).transpose();
            box2=(mpr.patch(1).eval(corG)).transpose();
            gsInfo
                <<"test failed on"
                <<"\nbox1:\n"<<box1
                <<"\nbox2:\n"<<box2
                <<"\nnumber of interfaces: "<<inter.size();
            if (inter.size()==1)
            {
                boundaryInterface computed = inter[0].first()==current.first()?inter[0] : inter[0].getInverse();
                gsInfo
                    <<"\nmap: "<<current.first().side()<<"->"<<current.second().side()
                    <<" vs "<<computed.first().side()<<"->"<<computed.second().side()
                    <<"\npermutation:\n"<<current.dirMap(current.first()).transpose()
                    <<" vs "<<computed.dirMap(current.first()).transpose()
                    <<"\norientation:\n"<<current.dirOrientation(current.first()).transpose()
                    <<" vs "<< computed.dirOrientation(current.first()).transpose() <<"\n";
            }
            pass=false;
            break;
        }
    } while(nextInterface(current,dim));
    gsInfo<<" "<<count<<" ";
    delete &basis;
    return pass;
}

int main()
{
    gsInfo<<"Test interfaces between segments\n"<<std::flush;
    TEST(test(1,0));
    gsInfo<<std::flush<<"\nTest interfaces between rectangles\n"<<std::flush;
    TEST(test(2,0));
    gsInfo<<std::flush<<"degenerate 1"<<std::flush; TEST(test(2,1));
    gsInfo<<std::flush<<"\nTest interfaces between cuboids\n"<<std::flush;
    TEST(test(3,0));
    gsInfo<<std::flush<<"degenerate 1"<<std::flush; TEST(test(3,1));
    gsInfo<<std::flush<<"degenerate 2"<<std::flush; TEST(test(3,2));
    return !passed;
}

