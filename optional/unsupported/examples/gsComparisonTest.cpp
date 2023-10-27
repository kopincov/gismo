/*
* gsComparisonTest.cpp created on 11.09.2014
*
* Author: Andrea Bressan
*
* This file is part of the G+SMO library
*/

#include <gismo.h>

#include <float.h>

using namespace gismo;

int main ()
{
real_t eps=2*math::limits::epsilon();
real_t a=1, b=1+.5, c=1+eps;
gsMatrix<> A(3,3);
gsMatrix<> B(3,3);
A<<0,1,2,3,4,5,6,7,8;
B<<0,1,2,3,4+7*FLT_EPSILON,5,6,7,8;

int error = 0;

if (! (gsClose(a,b,c)==true )) error |= 0x1;
if (! (gsClose(a,c,b)==true )) error |= 0x2;
if (! (gsClose(b,a,c)==true )) error |= 0x4;

if (! (gsAllCloseRelativeToMax(A,B,eps)!=true )) error |= 0x8;
//if (! (gsAllCloseRelativeToMax(A,B,eps*0.875)!=false )) error |= 0x10;


return error;
}
