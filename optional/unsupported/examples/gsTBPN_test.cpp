
#include <gismo.h>

#include <gsRemappedBasis/gsHLR.h>
#include <gsRemappedBasis/gsTBPN.h>
#include <gsRemappedBasis/gsTHB.h>

#include <gsRemappedBasis/otherUtils.h>

using namespace gismo;


struct Test
{
    std::vector<gsBoxList::basisPtr> bases;
    gsBoxList             boxes;

    Test()
        : boxes(2)
        {}
};

typedef Test (*testFun) ();

Test test1();
Test test2();
Test test3();
Test test4();

testFun tests[] = {test1,test2,test3,test4};
int numTests=sizeof(tests)/sizeof(testFun);

int main (int argn, char** args)
{
    bool passed=true;
    for (int t=0; t<numTests;++t)
    {
        std::cout<<"Test "<<t<<std::endl;
        std::stringstream base;
        base<<"test_"<<t<<"_";
        std::string domName=base.str()+"dom";
        std::string funName=base.str()+"fun";

        Test     curData = tests[t]();
        gsTBPN<2> basis(curData.bases,curData.boxes);
        basis.exportSelectorToTex(domName);

        /* // Note: this requires a proper latex installation
        std::string pdfCommand("pdflatex ");
        pdfCommand+=domName+".tex  >/dev/null";
        int pdfErr=system(pdfCommand.c_str());
        GISMO_UNUSED(pdfErr);
        */

        passed = passed && basis.test();

        printAllFunctions(basis,funName,1000);
    }

    return !passed;
}



Test test1()
{
    Test result;

    gsMatrix<> box(2,2);

    // example 1
    const int deg1=2;

    gsKnotVector<> knots1x(0,1,5,deg1+1,1);
    gsKnotVector<> knots1xr=knots1x;
    knots1xr.uniformRefine();

    result.bases.push_back( memory::make_shared(new gsTensorBSplineBasis<2>(knots1x,knots1x )));
    result.bases.push_back( memory::make_shared(new gsTensorBSplineBasis<2>(knots1xr,knots1x)));
    result.bases.push_back( memory::make_shared(new gsTensorBSplineBasis<2>(knots1x,knots1xr)));
    result.bases.push_back( memory::make_shared(new gsTensorBSplineBasis<2>(knots1xr,knots1xr)));

    box<<knots1x[0],knots1x[10],
         knots1x[0],knots1x[5];
    result.boxes.append(box,0);

    box<<knots1x[5],knots1x[10],
         knots1x[0],knots1x[5];
    result.boxes.append(box,1);

    box<<knots1x[0],knots1x[5],
         knots1x[5],knots1x[10];
    result.boxes.append(box,2);

    box<<knots1x[5],knots1x[10],
         knots1x[5],knots1x[10];
    result.boxes.append(box,3);

    return result;
}


Test test2()
{
    Test result;

    gsMatrix<> box(2,2);


    const int deg2=2;

    gsKnotVector<> knots2x(0,1.5,5,deg2+1,1);
    gsKnotVector<> knots2xr=knots2x;
    knots2xr.uniformRefine();
    gsKnotVector<> knots2y(0,1,3,deg2+1,1);
    gsKnotVector<> knots2yr=knots2y;
    knots2yr.uniformRefine();

    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots2x,knots2y)));
    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots2xr,knots2y)));
    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots2x,knots2yr)));

    box<<knots2x[4],knots2x[6],
         knots2y[0],knots2y[8];
    result.boxes.append(box,0);

    box<<knots2x[0],knots2x[4],
         knots2y[0],knots2y[8];
    result.boxes.append(box,1);

    box<<knots2x[6],knots2x[10],
         knots2y[0],knots2y[8];
    result.boxes.append(box,2);

    return result;
}


Test test3()
{
    Test result;

    gsMatrix<> box(2,2);

    const int deg3=2;

    gsKnotVector<> knots3x(0,1,3,deg3+1,1);
    gsKnotVector<> knots3xr=knots3x;
    knots3xr.uniformRefine();
    gsKnotVector<> knots3xrr=knots3xr;
    knots3xrr.uniformRefine();

    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots3x,knots3x)));
    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots3xr,knots3xr)));
    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots3xrr,knots3xrr)));

    box<<knots3x[0],knots3x[8],
         knots3x[0],knots3x[8];
    result.boxes.append(box,0);

    box<<knots3x[3],knots3x[8],
         knots3x[3],knots3x[8];
    result.boxes.append(box,1);

    box<<knots3xr[4],knots3xr[7],
         knots3xr[4],knots3xr[7];
    result.boxes.append(box,2);

    return result;
}

Test test4()
{
    Test result;

    gsMatrix<> box(2,2);

    const int deg4=2;

    gsKnotVector<> knots4x(0,1.4,6,deg4+1,1);
    gsKnotVector<> knots4xr=knots4x;
    knots4xr.uniformRefine();
    gsKnotVector<> knots4y(0,1,4,deg4+1,1);
    gsKnotVector<> knots4yr=knots4y;
    knots4yr.uniformRefine();
    gsKnotVector<> knots4yrr=knots4yr;
    knots4yrr.uniformRefine();

    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots4x,knots4y)));
    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots4xr,knots4y)));
    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots4x,knots4yr)));
    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots4xr,knots4yr)));
    result.bases.push_back(memory::make_shared(new gsTensorBSplineBasis<2>(knots4xr,knots4yrr)));


    box<<knots4x[0],knots4x[5],
         knots4y[4],knots4y[9];
    result.boxes.append(box,1);


    box<<knots4x[5],knots4x[11],
         knots4y[0],knots4y[4];
    result.boxes.append(box,2);

    box<<knots4x[5],knots4x[8],
         knots4y[4],knots4y[6];
    result.boxes.append(box,3);


    box<<knots4xr[11],knots4xr[14],
         knots4y[4],knots4y[5];
    result.boxes.append(box,4);
    return result;
}
