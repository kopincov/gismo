/** @file gsNorm_test.cpp

    @brief testing the different norms inheriting from gsNorm (gsNormL2,gsNormL2Boundary,gsSeminormH1,gsSeminormH2).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>
#include <iostream>
#include <gsAssembler/gsNormL2.h>
#include <gsAssembler/gsSeminormH1.h>
#include <gsAssembler/gsSeminormH2.h>
#include <gsAssembler/gsNormL2Boundary.h>

using namespace gismo;

// This test uses the identity as geometrical mapping, which relies on the correct implementation
// of the transformation methods. These are tested in a different test.

real_t tolerance=1000*std::numeric_limits<real_t>::epsilon();

// helper function that gives back a geometry. it is the identity
void getGeomId(gsTensorBSpline<2,real_t>& square)
{
    gsKnotVector<real_t> knots(0,1,0,2,0);
    gsMatrix<real_t> coefs(4,2);
    coefs << 0, 0,  1, 0,
             0, 1,  1, 1;
    gsTensorBSplineBasis<2,real_t> basis(knots,knots);
    square=gsTensorBSpline<2,real_t>(basis,coefs);
}

// helper function that gives back a constant function
void getFunctionConstant(gsTensorBSpline<2,real_t>& square)
{
    gsKnotVector<real_t> knots(0,1,0,2,0);
    gsMatrix<real_t> coefs(4,1);
    coefs << 1, 1,  1, 1;
    gsTensorBSplineBasis<2,real_t> basis(knots,knots);
    square=gsTensorBSpline<2,real_t>(basis,coefs);
}

// helper function that gives back test function
void getFunctionTest2(gsTensorBSpline<2,real_t>& square)
{
    gsKnotVector<real_t> knots(0,1,0,4,0);
    gsMatrix<real_t> coefs(16,1);
    coefs << 1, 0, 0, 0,
             0, 2, 0, 0,
             0, 0, 4, 0,
             0, 2, 0, 0;
    gsTensorBSplineBasis<2,real_t> basis(knots,knots);
    square=gsTensorBSpline<2,real_t>(basis,coefs);
}

// function that tests the h2 seminorm
bool testSeminormH2(real_t exactValue,const gsField<real_t>* field,const gsFunction<real_t>* f)
{
    bool passed = true;
    gsSeminormH2<real_t>* h2Seminorm;
    if(f==NULL)
        h2Seminorm=new gsSeminormH2<real_t>(*field); // if no f is given, we use the other constructor
    else
        h2Seminorm=new gsSeminormH2<real_t>(*field,*f);
    if( !gsClose<real_t>(h2Seminorm->compute(),exactValue, tolerance ) )
        passed = false;
    delete h2Seminorm;
    return passed;
}

// function that tests the h1 seminorm
bool testSeminormH1(real_t exactValue,const gsField<real_t>* field,
                    const gsFunction<real_t>* f)
{
    bool passed = true;
    gsSeminormH1<real_t>* h1Seminorm;
    if(f==NULL)
        h1Seminorm=new gsSeminormH1<real_t>(*field); // if no f is given, we use the other constructor
    else
        h1Seminorm=new gsSeminormH1<real_t>(*field,*f);
    if( !gsClose<real_t>(h1Seminorm->compute(),exactValue, tolerance) )
        passed = false;
    delete h1Seminorm;
    return passed;
}

// function that tests the l2 norm
bool testNormL2(real_t exactValue,const gsField<real_t>* field,const gsFunction<real_t>* f)
{
    bool passed = true;
    gsNormL2<real_t>* L2Norm;
    if(f==NULL)
        L2Norm=new gsNormL2<real_t>(*field); // if no f is given, we use the other constructor
    else
        L2Norm=new gsNormL2<real_t>(*field,*f);
    if( !gsClose<real_t>(L2Norm->compute(),exactValue, tolerance) )
        passed = false;
    delete L2Norm;
    return passed;
}

// function that tests the l2 boundary norm
bool testNormL2Boundary(real_t exactValue,const gsField<real_t>* field,const gsFunction<real_t>* f)
{
    bool passed = true;
    real_t normValue = 0;
    GISMO_ASSERT(field->nPatches() == 1, "Test for boundary norm assumes single patch");

    //Get contribution from every side
    if(f==NULL)
        // if no f is given, we use the other constructor
        normValue = gsNormL2Boundary<real_t>(*field).compute();
    else
        normValue = gsNormL2Boundary<real_t>(*field,*f).compute();

    if( !gsClose<real_t>(normValue,exactValue, tolerance) )
        passed = false;
    return passed;
}

int main(int argc, char *argv[])
{
    bool passed = true;
    std::vector<gsField<real_t>*>fields;
    std::vector<gsFunctionWithDerivatives<real_t>*> fs;
    std::vector<real_t>exactValuesH2,exactValuesH1,exactValuesL2, exactValuesL2Boundary;

    // First test: norm of identity geometry mapping + constant function
    gsTensorBSpline<2,real_t>geomId1,funcId1;
    getGeomId(geomId1);
    gsMultiPatch<real_t> mpGeomId1(geomId1);
    getFunctionConstant(funcId1);
    gsMultiPatch<real_t> mpFuncId1(funcId1);
    gsField<real_t>fieldId1(mpGeomId1,mpFuncId1);
    fields.push_back(&fieldId1);
    fs.push_back(NULL);
    exactValuesL2Boundary.push_back(2.0);
    exactValuesL2.push_back(1.0);
    exactValuesH1.push_back(0.0);
    exactValuesH2.push_back(0.0);

    // Second test: norm of identity geometry mapping + bezier function
    gsTensorBSpline<2,real_t>geomId2,funcTest2;
    getGeomId(geomId2);
    gsMultiPatch<real_t> mpGeomId2(geomId2);
    getFunctionTest2(funcTest2);
    gsMultiPatch<real_t> mpFuncTest2(funcTest2);
    gsField<real_t>fieldId2(mpGeomId2,mpFuncTest2);
    fields.push_back(&fieldId2);
    fs.push_back(NULL);
    exactValuesL2Boundary.push_back(pow(22.0/35.0,0.5));
    exactValuesL2.push_back(1.0/35.0*pow(991.0/2.0,1.0/2.0));
    exactValuesH1.push_back(pow(69.0/14.0,1.0/2.0));
    exactValuesH2.push_back(12.0/5.0*pow(186.0/7.0,1.0/2.0));

    // Third test: distance of identity geometry mapping + bezier function to
    // given function: x^2+4*y^2+3*x*y
    gsTensorBSpline<2,real_t>geomId3,funcTest3;
    getGeomId(geomId3);
    gsMultiPatch<real_t> mpGeomId3(geomId3);
    getFunctionTest2(funcTest3);
    gsMultiPatch<real_t> mpFuncTest3(funcTest3);
    gsField<real_t>fieldId3(mpGeomId3,mpFuncTest3);
    fields.push_back(&fieldId3);
    gsFunctionExpr<real_t> vf3("x^2-4*y^2-3*x*y",2);
    gsFunctionExpr<real_t> df3("2*x-3*y","-8*y-3*x",2);
    gsFunctionExpr<real_t> ddf3("2","-8","-3",2);
    gsFunctionWithDerivatives<real_t> f3(vf3,df3,ddf3);
    
    fs.push_back(&f3);
    exactValuesL2Boundary.push_back(2.0*std::pow(1147.0/105.0,0.5));
    exactValuesL2.push_back(1.0/210.0*std::pow(373627.0,1.0/2.0));
    exactValuesH1.push_back(1.0/2.0*std::pow(18667.0/105.0,1.0/2.0));
    exactValuesH2.push_back(1.0/5.0*std::pow(41309.0/7.0,1.0/2.0));

    for(unsigned i = 0;i<fields.size();++i)
    {
        bool passed_i0 = testNormL2Boundary(exactValuesL2Boundary[i],fields[i],fs[i]);
        bool passed_i1 = testNormL2(exactValuesL2[i],fields[i],fs[i]);
        bool passed_i2 = testSeminormH1(exactValuesH1[i],fields[i],fs[i]);
        bool passed_i3 = testSeminormH2(exactValuesH2[i],fields[i],fs[i]);
        passed = passed && passed_i0 && passed_i1 && passed_i2 && passed_i3;
        gsInfo << "Test " << i << " has ";
        if(!(passed_i0 && passed_i1 && passed_i2 && passed_i3))
            gsInfo << "not ";
        gsInfo << "passed.\n";
    }
    return !passed;
}
