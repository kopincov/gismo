#include <iostream>

#include <gismo.h>

#include <gsUtils/gsDerivatives.hpp>



using namespace gismo;


bool passed=true;

#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL\n");\

int main ()
{
    bool tpass=true;
    gsInfo<<"Test number of derivatives:\n";
    tpass = tpass && computedDerivativeSize<1,1,1> () == computedDerivativeSize(1,1,1);
    tpass = tpass && computedDerivativeSize<1,1,3> () == computedDerivativeSize(1,1,3);
    tpass = tpass && computedDerivativeSize<1,1,5> () == computedDerivativeSize(1,1,5);
    tpass = tpass && computedDerivativeSize<2,1,0> () == computedDerivativeSize(2,1,0);
    tpass = tpass && computedDerivativeSize<2,1,1> () == computedDerivativeSize(2,1,1);
    tpass = tpass && computedDerivativeSize<2,1,2> () == computedDerivativeSize(2,1,2);
    tpass = tpass && computedDerivativeSize<2,1,3> () == computedDerivativeSize(2,1,3);
    tpass = tpass && computedDerivativeSize<2,1,4> () == computedDerivativeSize(2,1,4);
    tpass = tpass && computedDerivativeSize<3,1,0> () == computedDerivativeSize(3,1,0);
    tpass = tpass && computedDerivativeSize<3,1,1> () == computedDerivativeSize(3,1,1);
    tpass = tpass && computedDerivativeSize<3,1,2> () == computedDerivativeSize(3,1,2);
    tpass = tpass && computedDerivativeSize<3,1,3> () == computedDerivativeSize(3,1,3);
    tpass = tpass && computedDerivativeSize<3,1,4> () == computedDerivativeSize(3,1,4);
    tpass = tpass && computedDerivativeSize<3,1,5> () == computedDerivativeSize(3,1,5);
    TEST(tpass);


    gsInfo<<"Domain dimension = 2, order = 2"<<"\n";
    gsInfo<<"derivativeToIndex<2,1>(0,0)=Dxx "
       <<derivativeToIndex<2,1>(0,0)<<"\n";
    gsInfo<<"derivativeToIndex<2,1>(0,1)=Dxy "
       <<derivativeToIndex<2,1>(0,1)<<"\n";
    gsInfo<<"derivativeToIndex<2,1>(1,1)=Dyy "
       <<derivativeToIndex<2,1>(1,1)<<"\n";

    gsInfo<<"Domain dimension = 2, order = 3"<<"\n";
    gsInfo<<"derivativeToIndex<2,1>(0,0,0)=Dxxx "
       <<derivativeToIndex<2,1>(0,0,0)<<"\n";
    gsInfo<<"derivativeToIndex<2,1>(0,0,1)=Dxxy "
       <<derivativeToIndex<2,1>(0,0,1)<<"\n";
    gsInfo<<"derivativeToIndex<2,1>(0,1,1)=Dxyy "
       <<derivativeToIndex<2,1>(0,1,1)<<"\n";
    gsInfo<<"derivativeToIndex<2,1>(1,1,1)=Dyyy "
       <<derivativeToIndex<2,1>(1,1,1)<<"\n";

    gsInfo<<"Domain dimension = 3, order = 2"<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(0,0)=Dxx "
       <<derivativeToIndex<3,1>(0,0)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(0,1)=Dxy "
       <<derivativeToIndex<3,1>(0,1)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(0,2)=Dxz "
       <<derivativeToIndex<3,1>(0,2)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(1,1)=Dyy "
       <<derivativeToIndex<3,1>(1,1)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(1,2)=Dyz "
       <<derivativeToIndex<3,1>(1,2)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(2,2)=Dzz "
       <<derivativeToIndex<3,1>(2,2)<<"\n";

    gsInfo<<"Domain dimension = 3, order = 3"<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(0,0,0)=Dxxx "
       <<derivativeToIndex<3,1>(0,0,0)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(0,0,1)=Dxxy "
       <<derivativeToIndex<3,1>(0,0,1)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(0,0,2)=Dxxz "
       <<derivativeToIndex<3,1>(0,0,2)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(0,1,1)=Dxyy "
       <<derivativeToIndex<3,1>(0,1,1)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(0,1,2)=Dxyz "
       <<derivativeToIndex<3,1>(0,1,2)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(0,2,2)=Dxzz "
       <<derivativeToIndex<3,1>(0,2,2)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(1,1,1)=Dyyy "
       <<derivativeToIndex<3,1>(1,1,1)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(1,1,2)=Dyyz "
       <<derivativeToIndex<3,1>(1,1,2)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(1,2,2)=Dyzz "
       <<derivativeToIndex<3,1>(1,2,2)<<"\n";
    gsInfo<<"derivativeToIndex<3,1>(2,2,2)=Dzzz "
       <<derivativeToIndex<3,1>(2,2,2)<<"\n";


    gsInfo<<"Domain dimension = 4, order = 2"<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,0)=Dxx "
       <<derivativeToIndex<4,1>(0,0)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,1)=Dxy "
       <<derivativeToIndex<4,1>(0,1)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,2)=Dxz "
       <<derivativeToIndex<4,1>(0,2)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,3)=Dxt "
       <<derivativeToIndex<4,1>(0,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(1,1)=Dyy "
       <<derivativeToIndex<4,1>(1,1)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(1,2)=Dyz "
       <<derivativeToIndex<4,1>(1,2)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(1,3)=Dyt "
       <<derivativeToIndex<4,1>(1,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(2,2)=Dzz "
       <<derivativeToIndex<4,1>(2,2)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(2,3)=Dzt "
       <<derivativeToIndex<4,1>(2,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(3,3)=Dtt "
       <<derivativeToIndex<4,1>(3,3)<<"\n";

    gsInfo<<"Domain dimension = 4, order = 3"<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,0,0)=Dxxx "
       <<derivativeToIndex<4,1>(0,0,0)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,0,1)=Dxxy "
       <<derivativeToIndex<4,1>(0,0,1)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,0,2)=Dxxz "
       <<derivativeToIndex<4,1>(0,0,2)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,0,3)=Dxxt "
       <<derivativeToIndex<4,1>(0,0,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,1,1)=Dxyy "
       <<derivativeToIndex<4,1>(0,1,1)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,1,2)=Dxyz "
       <<derivativeToIndex<4,1>(0,1,2)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,1,3)=Dxyt "
       <<derivativeToIndex<4,1>(0,1,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,2,2)=Dxzz "
       <<derivativeToIndex<4,1>(0,2,2)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,2,3)=Dxzt "
       <<derivativeToIndex<4,1>(0,2,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(0,3,3)=Dxtt "
       <<derivativeToIndex<4,1>(0,3,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(1,1,1)=Dyyy "
       <<derivativeToIndex<4,1>(1,1,1)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(1,1,2)=Dyyz "
       <<derivativeToIndex<4,1>(1,1,2)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(1,1,3)=Dyyt "
       <<derivativeToIndex<4,1>(1,1,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(1,2,2)=Dyzz "
       <<derivativeToIndex<4,1>(1,2,2)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(1,2,3)=Dyzt "
       <<derivativeToIndex<4,1>(1,2,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(1,3,3)=Dytt "
       <<derivativeToIndex<4,1>(1,3,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(2,3,2)=Dzzz "
       <<derivativeToIndex<4,1>(2,2,2)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(2,2,3)=Dzzt "
       <<derivativeToIndex<4,1>(2,2,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(2,3,3)=Dztt "
       <<derivativeToIndex<4,1>(2,3,3)<<"\n";
    gsInfo<<"derivativeToIndex<4,1>(3,3,3)=Dttt "
       <<derivativeToIndex<4,1>(3,3,3)<<"\n";

    return (passed ? 0: 1);
}
