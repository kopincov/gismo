/*
* gsDivPerservingTransformation.cpp created on 18.08.2014
*
* Author: Andrea Bressan, Jarle Sogn
*
* This file is part of the G+SMO library
*/


#include <iostream>
#include <gismo.h>
#include <gismo_dev.h>

#include <gsCore/gsVectorValuedFunctionSet.h>
#include <gsCore/gsTransformedFuncSet.h>



using std::flush;
using std::sqrt;
using std::ostringstream;
using std::string;
using namespace gismo;

#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL\n")<<std::flush;

bool passed=true;

real_t GISMO_EPS = 10*std::numeric_limits<real_t>::epsilon();

void run_test (gsGeometry<> &geo);

int main ()
{
    gsGeometry<>::uPtr geo;
    gsTensorBSpline<2,real_t>::uPtr geo2D;

    gsInfo << "Staring 2D test: \n";

    gsInfo << "\n\nTesting big square (scaled geometry):\n";
    geo   =  gsNurbsCreator<>::BSplineSquare(2.0,0,0);
    run_test(*geo);

    gsInfo << "\n\nTesting translated unit square (translated geometry):\n";
    geo   =  gsNurbsCreator<>::BSplineSquare(1.0, 0.5, 1.5);
    run_test(*geo);

    gsInfo << "\n\nTesting rotated rectange (2D geometry):\n";
    geo   =   gsNurbsCreator<>::BSplineRectangle(0.0, 0.0, 2.0, 1.0, 90);
    run_test(*geo);

    gsInfo << "\n\nTesting Fat Quarter Annulus (2D geometry):\n";
    geo   =   gsNurbsCreator<>::BSplineFatQuarterAnnulus();
    run_test(*geo);

    gsInfo << "2D test comlete: ";  TEST(passed);


    gsInfo << "\nStaring 3D test: \n";

    gsInfo << "\n\nTesting  big cube (scaled geometry):\n";
    geo   =  gsNurbsCreator<>::BSplineCube(2.0,0,0,0);
    run_test(*geo);

    gsInfo << "\n\nTesting translated unit cube (translated geometry):\n";
    geo   =  gsNurbsCreator<>::BSplineCube(1.0, 0.5, 1.5,-1.0);
    run_test(*geo);


    gsInfo << "\n\nTesting lifted Quarter Annulus (3D geometry):\n";
    geo2D =  gsNurbsCreator<>::BSplineFatQuarterAnnulus();
    geo   =  gsNurbsCreator<>::lift3D(*geo2D, 1.0);
    run_test(*geo);

    //This test is for testing non-constant 3-direction
    gsInfo << "\n\nTesting lifted and rotated Quarter Annulus (3D geometry):\n";
    gsVector<> xAxis(3); xAxis(0) = 1; xAxis(1) = 0; xAxis(2) = 0;
    geo2D =  gsNurbsCreator<>::BSplineFatQuarterAnnulus();
    geo   =  gsNurbsCreator<>::lift3D(*geo2D, 1.0);
    geo->rotate(1.5707963267948966, xAxis);
    run_test(*geo);

    gsInfo << "3D test comlete: ";  TEST(passed);

    gsInfo << "\nStaring surface test: \n";

    gsInfo << "\n\nTesting on surface:\n";
    geo   =  gsReadFile<>("saddle.xml");
    gsWarn<<"TEST SUPPRESSED AS CODE IS NOT READY";
//    run_test(*geo);

    gsInfo << "surface test comlete: ";  TEST(passed);

    return passed==1?0:1;
}


void computeDivSillyWay (short_t domDim,const gsMatrix<real_t> &derivs, gsMatrix<real_t> &out)
{
    int step=domDim+1;
    int block_size=domDim*domDim;
    int blocks=derivs.rows()/block_size;

    out.setZero(blocks,derivs.cols());

    for (int bl=0;bl<blocks;++bl)
        for (int c=0;c<domDim;++c)
            out.row(bl)+=derivs.row(bl*block_size+c*step);
}


void run_test (gsGeometry<> &geo)
{
    gsBasis<>::domainIter domIt = geo.basis().makeDomainIterator();
    gsVector<index_t> numNodes;
    numNodes.resize(geo.parDim());
    numNodes.setConstant(3);

    //domIt->computeQuadratureRule(numNodes);
    gsGaussRule<> quRule(numNodes);
    gsMatrix<> quNodes;
    gsVector<> quWeights;
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

    std::vector<gsBasis<> *> space_vec(geo.parDim(), &geo.basis());

    gsVBasis<real_t> my_no_trans(space_vec);
    gsDivConformingTFS my_div_conf(my_no_trans);
    gsMapData<> geo_vals(NEED_VALUE | NEED_GRAD_TRANSFORM | NEED_2ND_DER | NEED_MEASURE);
    gsFuncData<> bas_vals_div(NEED_VALUE | NEED_DERIV | NEED_DIV);

    bool d_ok = true; // divergence
    gsMatrix<> div_Error;
    gsMatrix<> div_silly;

    real_t l2_error = 0;

    for (; domIt->good(); domIt->next())
    {
        geo_vals.points = quNodes;
        geo.computeMap(geo_vals);
        my_div_conf.compute(geo_vals, bas_vals_div);
        computeDivSillyWay(geo.domainDim(), bas_vals_div.values[1], div_silly);

        div_Error = bas_vals_div.divs - div_silly;

        real_t l2_error_loc = div_Error.squaredNorm();
        l2_error += l2_error_loc;
        if (l2_error_loc > GISMO_EPS)
        {
            gsDebug << "div_Piola:\n " << div_silly << "\n";
            gsDebug << "div_No:\n " << bas_vals_div.divs << "\n";
            gsDebug << "div_Piola - div_No:\n " << div_Error << "\n";
            gsDebug << "quNodes:\n " << quNodes << "\n\n\n";
        }
        d_ok = d_ok && l2_error_loc < GISMO_EPS;

    }
    gsInfo << "l_2 error is: " << math::sqrt(l2_error) << "\n";
    gsInfo << "DIVS ";
    TEST(d_ok);
}
