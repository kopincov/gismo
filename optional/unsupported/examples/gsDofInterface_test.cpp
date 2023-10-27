/**  gsInterfaceGlue.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):   A. Bressan, A. Manzaflaris
    Created on:  2015-02-02

    Test the functions that are identified by the gsDofMapper actually
    agree on the interface.
*/


#include <gismo.h>

#include "interfaceTestUtils.h"

using namespace gismo;
bool temp_pass;
bool passed=true;

#define TEST(a)\
    temp_pass = (a); passed = (passed && temp_pass);\
    gsInfo << (temp_pass? "TEST OK\n":"TEST FAIL\n");
bool test (int dim);

int main ()
{
    gsInfo<<"Test interfaces between segments\n"<<std::flush;
    TEST(test(1));
    gsInfo<<std::flush<<"\nTest interfaces between rectangles\n"<<std::flush;
    TEST(test(2));
    gsInfo<<std::flush<<"\nTest interfaces between cuboids\n"<<std::flush;
    TEST(test(3));
    return !passed;
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

bool test (int dim)
{
    bool pass=true;
    gsBasis<>     &basis    = getDimLinearBasis(dim);

    gsMatrix<>     box=basis.support();
    gsMatrix<>     cor=getCorners(box);
    gsMatrix<> cor2, _min, _max, cp1, cp2;
    gsGeometry<>::uPtr geo = basis.makeGeometry(cor.transpose());

    gsMatrix<index_t> index1, index2;
    gsDofMapper dofM;
    gsVector<index_t> np;
    np.setConstant(dim,3*dim);

    boundaryInterface current;
    firstInterface(current,dim);
    index_t count=0;
    do
    {
        ++count;
        gsMultiPatch<> mpr;
        mpr.clear();
        mpr.addPatch(geo->clone());
        mpr.addPatch(geo->clone());
        // change coefficients for the second patch
        gsAffineFunction<real_t>  interfaceMapInv=mpr.getMapForInterface(current.getInverse());
        gsAffineFunction<real_t>  interfaceMap=mpr.getMapForInterface(current);

        interfaceMapInv.eval_into(cor,cor2);
        mpr.patch(1).coefs()=cor2.transpose();
        mpr.computeTopology(10*math::limits::epsilon());

        // construct the multibasis
        gsMultiBasis<> mb(mpr);
        for (int i=0; i<dim; ++i)
        {
            mb[0].component(i).uniformRefine(i);
            mb[1].component(current.dirMap()(i)).uniformRefine(i);
        }

        if (count == 0)
            gsInfo<< "Dofs per patch: "<< mb.size(0) <<"\n";

        // set up evaluation points
        _min=box.col(0);
        _max=box.col(1);
        if (current.first().parameter() )
            _min(current.first().direction())=_max(current.first().direction());
        else
            _max(current.first().direction())=_min(current.first().direction());
        gsGridIterator<real_t, CUBE> gridPoint(_min, _max, np );
        
        index1 = mb[0].boundary(current.first().side());
        index2 = mb[1].boundary(current.second().side());

        mb.getMapper(true,dofM);

        for (index_t i=0; i<index1.size(); ++i)
        {
            const index_t gIndex1=dofM.index(index1(i,0),0);
            for (index_t j=0; j<index2.size();++j)
            {
                const index_t gIndex2=dofM.index(index2(j,0),1);
                if (gIndex1!=gIndex2)
                    continue;
                for(gridPoint.reset(); gridPoint; ++gridPoint )
                {
                    const gsMatrix<> & p = *gridPoint;
                    mb[0].evalSingle_into(index1(i,0), p, cp1);
                    mb[1].evalSingle_into(index2(j,0),interfaceMap.eval(p), cp2);
                    if ( !gsAllCloseAbsolute(cp1,cp2,100*std::numeric_limits<real_t>::epsilon()))
                    {
                        gsDebugVar( current );
                        gsDebugVar( index1(i,0) );
                        gsDebugVar( index2(j,0) );
                        gsDebugVar( cp1.transpose() );
                        gsDebugVar( cp2.transpose() );
                        GISMO_ERROR("gluing in dof mapper is not correct");
                    }
                 }
            }
        }

    } while(nextInterface(current,dim));
    gsInfo<<" "<<count<<" ";
    delete &basis; // (!) scary 
    return pass;
}
