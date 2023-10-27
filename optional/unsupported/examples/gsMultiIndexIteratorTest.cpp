/*
* gsMultiIndexIteratorTest.cpp created on 10.07.2014
*
* Author: Andrea Bressan
*
* This file is part of the G+SMO library
*/

#include <gismo.h>

#include <gsUtils/gsCompositionIterator.h>
#include <gsUtils/gsSimplexIterator.h>


using std::flush;
using namespace gismo;

///TODO: check the setToFlat
///TODO: check the setToMulti
///TODO: make the test work also for small cubes

bool mytemp_var;
#define TEST(a)\
    mytemp_var = (a);\
    passed = (passed && mytemp_var);\
    gsInfo << (mytemp_var? "OK":"FAIL")<<flush;\


template <int dim>
bool test_tensor_grid  ( index_t low_s, index_t upp_s);
template <int dim>
bool test_compositions ( index_t sum);

int main()
{
    bool passed = true;
    TEST(test_tensor_grid<1> (100, 108));
    TEST(test_tensor_grid<2> (100, 107));
    TEST(test_tensor_grid<3> (100, 106));
    TEST(test_tensor_grid<4> (100, 105));
    TEST(test_tensor_grid<5> (100, 105));

    for(int i=1;i<8;++i)
    {
        TEST(test_compositions<1>(i));
        TEST(test_compositions<2>(i));
        TEST(test_compositions<3>(i));
        TEST(test_compositions<4>(i));
        TEST(test_compositions<5>(i));
    }

    gsInfo<<"\nTOTAL:"; TEST(passed); gsInfo<<"\n";
    return ( passed ? 0 : 1 );
}

template <int dim>
bool test_tensor_grid( index_t low_s, index_t upp_s)
{
    gsVector<index_t,dim>  low, upp, str;
    for (index_t i=0; i<dim; ++i)
    {
        low(i)=low_s+i;
        upp(i)=upp_s+2*i;
    }
    
    gsGridIterator<index_t,CUBE,dim>     cube(low,upp);
    gsGridIterator<index_t,CUBE,-1 >     cube_d(low,upp);        
    gsGridIterator<index_t,BDR,dim>      boundary(low,upp);
    gsGridIterator<index_t,BDR,-1>       boundary_d(low,upp);
    gsGridIterator<index_t,VERTEX,dim>   vertices(low,upp);
    gsGridIterator<index_t,VERTEX,-1>    vertices_d(low,upp);

    bool passed      = true; // all ok?
    bool countOk     = true; // are the points the right number?
    bool fixedDimOk  = true; // dimension specialized iterator and generic one agree?
    bool convMTFOk   = true; // conversion multi to flat
    long int  count = 0;
    
    gsInfo<<"\n\nTesting tensor iterator of dimension "<<dim;
    gsInfo<<"\nlowCorner:("<<low.transpose()<<")  uppCorner:("<<upp.transpose()<<")\n";

    gsInfo<< "  Test Grid Points:";
    str = cube.strides();
    for (cube_d.reset(), cube.reset(); cube && cube_d; ++cube, ++cube_d )
    {
        fixedDimOk = fixedDimOk && *cube == *cube_d;
        convMTFOk  = convMTFOk  &&  str.dot(*cube - cube.lower()) == count;
        count++;
    }
    countOk = countOk && count==(upp-low).prod()
        && count==cube.numPoints()
        && ( (bool)cube == (bool)cube_d);
    gsInfo<< "\n    PointCount: "; TEST(countOk);
    gsInfo<< "; IndexConversion: "; TEST( convMTFOk );
    gsInfo <<"; Specialization: "; TEST(fixedDimOk); gsInfo<<"\n";

    countOk = fixedDimOk = true;
    count = 0;
    gsInfo<< "  Test Vertices Points:";
    for (vertices.reset(), vertices_d.reset(); vertices && vertices_d; ++vertices, ++vertices_d)
    {
        fixedDimOk = fixedDimOk && *vertices == *vertices_d;
        count++;
    }
    countOk = countOk && count==(1<<dim)
        && count==vertices.numPoints()
        && ( (bool)vertices == (bool)vertices_d );
    gsInfo<< "\n    PointCount: "; TEST(countOk);
    gsInfo <<"; Specialization: "; TEST(fixedDimOk); gsInfo<<"\n";

    countOk = fixedDimOk = true;
    count = 0;
    gsInfo<< "  Test Boundary Points:";
    for (boundary.reset(), boundary_d.reset(); boundary && boundary_d; ++boundary, ++boundary_d)
    {
        fixedDimOk= fixedDimOk && *boundary == *boundary_d;
        count++;
    }
    countOk = countOk && count==(upp-low).prod()-((upp-low).array()-2).prod()
        && count==boundary.numPoints() && ((bool)boundary == (bool)boundary_d);
    gsInfo<< "\n    PointCount: "; TEST(countOk);
    gsInfo <<"; Specialization: "; TEST(fixedDimOk); gsInfo<<"\n";

    return passed;
}

template <int dim>
bool test_compositions ( index_t sum)
{
    gsInfo<<"\n\nTesting composition iterator of dimension "<<dim<<" for sum "<< sum;
    gsInfo<<"\nand simplex iterator of size "<<sum<<" in Z^"<<dim<<"\n";

    gsCompositionIterator<index_t,dim+1> compositions(sum);
    gsCompositionIterator<index_t>       COMPOSITIONS(sum,dim+1);
    gsSimplexIterator<index_t,dim>       simplex(sum);
    gsSimplexIterator<index_t>           SIMPLEX(sum,dim);

    bool passed = true;

    bool fixedDimOk = true;
    bool simOkMult  = true;
    long int  count = 0;

    for (compositions.reset(), simplex.reset(), COMPOSITIONS.reset(), SIMPLEX.reset();
         compositions && simplex && COMPOSITIONS && SIMPLEX;
         ++compositions, ++simplex, ++COMPOSITIONS, ++SIMPLEX )
    {
        fixedDimOk = fixedDimOk && *compositions==*COMPOSITIONS;
        fixedDimOk = fixedDimOk && *simplex==*SIMPLEX;

        simOkMult  = simOkMult  && compositions->block(1,0,dim,1)==*simplex;
        simOkMult  = simOkMult  && COMPOSITIONS->block(1,0,dim,1)==*SIMPLEX;

        count++;
    }
    fixedDimOk = fixedDimOk && (bool)SIMPLEX==(bool)simplex && (bool)compositions==(bool)COMPOSITIONS;
    gsInfo<< "\n    Compositions  :";  gsInfo<<" "; TEST( simOkMult );
    TEST( count== compositions.numPoints() ); gsInfo<<" ";
    gsInfo<< "\n    Simplex:      :";
    TEST( (bool)compositions==(bool)simplex );
    gsInfo<<" "; TEST( count== simplex.numPoints() );
    gsInfo<< "\n    Specialization: "; TEST(fixedDimOk); gsInfo<<"\n";

    return passed;
}
