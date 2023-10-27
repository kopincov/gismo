/*
* gsPointGridTest.cpp created on 15.07.2014
*
* Author: Andrea Bressan
*
* This file is part of the G+SMO library
*/

#include <gismo.h>

#include <gsUtils/gsSimplexIterator.h>

#include <iostream>

using namespace gismo;

bool passed;
#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL")<<"\n";\


int main (void)
{

gsInfo<<"\n\nTEST ON TRIANGLE GRID\n\n";

gsVector<real_t,2> v0,v1,v2;
gsMatrix<real_t,2,3> V;

v0(0)=0; v0(1)=0;
v1(0)=10; v1(1)=10;
v2(0)=20; v2(1)=0;

V.col(0)=v0;
V.col(1)=v1;
V.col(2)=v2;
gsSimplexIterator<real_t>     mySimplexPoint(V,5);

index_t i = 0;
for (; mySimplexPoint; ++mySimplexPoint )
{
    gsInfo<<"point "<<i++<<": "<<mySimplexPoint->transpose()
          << (mySimplexPoint.isBoundary() ? " on boundary\n" : "\n");
}
GISMO_ENSURE(i == mySimplexPoint.numPoints() ,"Something went wrong");
    
gsInfo<<"\n\nTEST ON TENSOR GRID\n\n";

gsVector<index_t,2> numP;
numP.setConstant(2,5);
gsGridIterator<real_t,CUBE>     myTensorPoint(v0,v1,numP);

i = 0;
for (; myTensorPoint; ++myTensorPoint)
{
    gsInfo<<"point "<<i++<<": "<<myTensorPoint->transpose()<<"\n";
}
GISMO_ENSURE(i == myTensorPoint.numPoints() ,"Something went wrong");

gsInfo<<"\n\nTEST ON TENSOR GRID BOUNDARY\n\n";

gsGridIterator<real_t,BDR,-1>  myB(v0,v1,numP);

i = 0;
for (; myB; ++myB)
{
    gsInfo<<"point "<<i++<<": "<<myB->transpose()<<"\n";
    GISMO_ASSERT(myB.isBoundary(),"Something went wrong");
}
GISMO_ENSURE(i == myB.numPoints() ,"Something went wrong");

gsInfo<<"\n\nTEST ON TENSOR GRID 3D\n\n";

gsVector<index_t,3> numP3d;
numP3d.setConstant(3,5);
gsVector<real_t,3>  low,upp;
low.setConstant(3,10);
upp.setConstant(3,20);

gsGridIterator<real_t,CUBE>     myTP3d(low,upp,numP3d);

i = 0;
for (myTP3d.reset(); myTP3d; ++myTP3d)
{
    gsInfo<<"point "<<i++<<": "<<myTP3d->transpose()<<"\n";
}
GISMO_ENSURE(i == myTP3d.numPoints() ,"Something went wrong");

gsInfo<<"\n\nTEST ON TENSOR GRID BOUNDARY 3D\n\n";

gsGridIterator<real_t,BDR>  myB3d(low,upp,numP3d);

i = 0;
for (; myB3d; ++myB3d)
{
    gsInfo<<"point "<<i++<<": "<<myB3d->transpose()<<"\n";
    GISMO_ASSERT(myB3d.isBoundary(),"Something went wrong");
}
GISMO_ENSURE(i == myB3d.numPoints() ,"Something went wrong");

gsInfo<<"\n\nVertices of 3D cube\n\n";

gsGridIterator<real_t,VERTEX>  myV3d(low,upp,numP3d);

i = 0;
for (; myV3d; ++myV3d)
{
    gsInfo<<"point "<<i++<<": "<<myV3d->transpose()<<"\n";
    GISMO_ASSERT(myV3d.isBoundary(),"Something went wrong");
}
GISMO_ENSURE(i == myV3d.numPoints() ,"Something went wrong");


gsInfo<<"\n\nGrid of 3D points given coordinate-wise\n\n";

std::vector<gsMatrix<int> > cwise(3, gsVector<int>::LinSpaced(2,1,2));
cwise[1] *= 2;
cwise[2] *= 3;
i = 0;
gsGridIterator<int,CWISE>  myCwise(cwise);
for (; myCwise; ++myCwise)
{
    gsInfo<<"point "<<i++<<": "<<myCwise->transpose()<<"\n";
}
GISMO_ENSURE(i == myCwise.numPoints() ,"Something went wrong");

return 0;

}
