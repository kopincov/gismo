/*
* gsBasisEvalTest2.cpp created on 05.08.2014
*
* Author: Andrea Bressan
*
* This file is part of the G+SMO library
*/


#include <gismo.h>
#include <gismo_dev.h>
#include <gsCore/gsBasisEvaluator.h>
using namespace gismo;


#include <iostream>

;
using std::flush;
using std::ostringstream;
using std::string;


#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL\n");

bool passed=true;

void run_test (gsBasis<>::domainIter &domIt, gsGeometry<> &geo, const gsMatrix<> &quNodes, gsFuncData<double> &bdata );
bool check_gradients (gsMatrix<> a, gsMatrix<> b);

int main ()
{
    gsVector<index_t> numNodes;
    gsGeometry<>::uPtr geo;
    gsBasis<>::domainIter domIt;
    gsQuadRule<> quRule;
    gsMatrix<> quNodes;
    gsVector<> quWeights;
    gsFuncData<> bdata(NEED_DERIV | NEED_ACTIVE | SAME_ELEMENT);

    geo = gsNurbsCreator<>::BSplineFatCircle();// BSpline Curve
    geo->uniformRefine(2);
    domIt = geo->basis().makeDomainIterator();
    numNodes.resize(1);
    numNodes.setConstant(1);
    quRule = gsGaussRule<>(numNodes);
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    run_test(domIt, *geo, quNodes, bdata);

    geo = gsNurbsCreator<>::BSplineFatQuarterAnnulus();// Planar TensorBSpline patch
    geo->uniformRefine(2);
    domIt = geo->basis().makeDomainIterator();
    numNodes.resize(2);
    numNodes.setConstant(1);
    quRule = gsGaussRule<>(numNodes);
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    run_test(domIt, *geo, quNodes, bdata);
    domIt.reset();

    geo = gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(2));// Planar TensorBSpline patch
    geo->uniformRefine(2);
    domIt = geo->basis().makeDomainIterator();
    numNodes.resize(2);
    numNodes.setConstant(1);
    quRule = gsGaussRule<>(numNodes);
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    run_test(domIt, *geo, quNodes, bdata);
    domIt.reset();


    geo = gsNurbsCreator<>::BSplineCube(static_cast<short_t>(2)); // 3D TensorBSpline patch
    geo->uniformRefine(2);
    domIt = geo->basis().makeDomainIterator();
    numNodes.resize(3);
    numNodes.setConstant(1);
    quRule = gsGaussRule<>(numNodes);
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    run_test(domIt, *geo, quNodes, bdata);
    domIt.reset();

    geo = gsReadFile<>("thbs_01.xml");  // 2D hier. truncated patch
    domIt = geo->basis().makeDomainIterator();
    numNodes.resize(2);
    numNodes.setConstant(1);
    quRule = gsGaussRule<>(numNodes);
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    run_test(domIt, *geo, quNodes, bdata);
    domIt.reset();

    return passed == 1 ? 0 : 1;
}


void run_test (gsBasis<>::domainIter &domIt, gsGeometry<> &geo, const gsMatrix<> &quNodes, gsFuncData<double> &bdata )
{
    gsBasis<> &space = geo.basis();
    gsBasisEvaluator<> &my_eval = *makeBasisEvaluator(space, NEED_VALUE | NEED_GRAD, &geo);

    gsGeometryEvaluator<>::uPtr geoEval(getEvaluator(NEED_VALUE | NEED_GRAD_TRANSFORM, geo));

    gsInfo << "\n*********** Running for " << space
           << "on " << geo << "\n";

    bool v_ok = true; // valued
    bool d_ok = true; // derivatives
    bool t_ok = true; // transformed gradients
    bool s_ok = true;
    GISMO_UNUSED(s_ok); // seecond derivatives

    gsMatrix<> grads;

    for (; domIt->good(); domIt->next())
    {
        geo.basis().compute(quNodes, bdata);
        geoEval->evaluateAt(quNodes);
        my_eval.evaluateAt(quNodes);

        v_ok = v_ok && gsAllCloseRelativeToMax(my_eval.values(), bdata.values[0], 1e-13);
        d_ok = v_ok && gsAllCloseRelativeToMax(my_eval.derivs(), bdata.values[1], 1e-13);

        geoEval->transformGradients(0, bdata.values[1], grads);
        my_eval.evaluateAt(quNodes, *geoEval);
        v_ok = v_ok && gsAllCloseRelativeToMax(my_eval.values(), bdata.values[0], 1e-13);
        t_ok = t_ok && check_gradients(my_eval.derivs(), grads);
        // s_ok= s_ok && my_eval.getDerivs (2) == domIt->basisDerivs (2);
    }
    gsInfo << "VALUES "; TEST(v_ok);
    gsInfo << "DERIVS "; TEST(d_ok);
    gsInfo << "DERIV2 " << "NOT_RUN" << "\n"; // TEST(s_ok);
    gsInfo << "TGRADS "; TEST(t_ok);
    delete &my_eval;
}

bool check_gradients (gsMatrix<> derivs, gsMatrix<> grads)
{
    unsigned row_num=grads.rows();
    unsigned active_num=grads.cols();
    bool ok=true;
    for (unsigned i=0; i<active_num && ok;++i)
    {
        ok = ok && gsAllCloseRelativeToMax(grads.col(i), derivs.block(i*row_num,0,row_num,1),1e-13);
    }
    return ok;
}


// the following test functions are templated in order to be able to transform between the
// domain iterator format and the standard format for transformed derivatives

template <short_t ParDim, short_t TarDim> void run_scalar_test(char * file_basename);


template <short_t ParDim, short_t TarDim>
void run_test_scalar ()
{
    run_scalar_test<ParDim,TarDim>("basis_eval_identity_test_");
    run_scalar_test<ParDim,TarDim>("basis_eval_affine_test_");
    run_scalar_test<ParDim,TarDim>("basis_eval_complex_test_");
}


template <short_t ParDim, short_t TarDim>
void run_scalar_test(char * file_basename)
{
    gsGeometry<> *geo = NULL;
    gsVector<index_t> numNodes;
    numNodes.setConstant(ParDim, 1);

    std::string fileName = file_basename;
    fileName += util::to_string(ParDim);
    fileName += "_";
    fileName += util::to_string(TarDim);
    fileName += ".xml";
    gsInfo << "Running test for scalar basis, data file: " << fileName << "\n";


    if (gsFileManager::fileExists(fileName))
    {
        gsGaussRule<> quRule(numNodes);
        gsMatrix<> quNodes;
        gsVector<> quWeights;
        gsFuncData<> bdata(NEED_DERIV | NEED_ACTIVE | SAME_ELEMENT);

        gsFileData<> dataFile(fileName);

        geo = dataFile.template getFirst<gsGeometry<> >().release();
        geo->uniformRefine(2);
        const gsBasis<> *basis = &geo->basis();

        gsBasis<>::domainIter domIt = basis->makeDomainIterator();
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
        gsGeometryEvaluator<>::uPtr geoEval(getEvaluator(0, *geo));
        gsBasisEvaluator<> *evalNoT = makeBasisEvaluator(*basis, static_cast<unsigned>(NEED_VALUE | NEED_GRAD));
        gsBasisEvaluator<> *evalIGA = makeBasisEvaluator(*basis, static_cast<unsigned>(NEED_VALUE | NEED_GRAD), geo, INVERSE_COMPOSITION);

        bool all_ok = true;
        bool a_ok = true; // active
        bool v_ok = true; // valued
        bool d_ok = true; // derivatives

        // TODO add test for transformed gradients
        // TODO add test for second derivatives
        // TODO add test for laplacian
        // TODO add test with try for things that should fail

        // TODO work with more than a single quadrature point

        for (; domIt->good() && all_ok; domIt->next())
        {
            basis->compute(quNodes, bdata);
            geoEval->evaluateAt(quNodes);
            evalNoT->evaluateAt(quNodes);
            evalIGA->evaluateAt(quNodes);

            a_ok = evalNoT->actives() == evalIGA->actives() && evalNoT->actives() == bdata.actives;//note: empty
            v_ok = evalNoT->values() == evalIGA->values() && evalNoT->values() == bdata.values[0];
            d_ok = evalNoT->derivs() == evalIGA->derivs() && evalNoT->derivs() == bdata.values[1];
            all_ok = a_ok && v_ok && d_ok;
        }
        gsInfo << "ACTIVES "; TEST(a_ok);
        gsInfo << "VALUES  "; TEST(v_ok);
        gsInfo << "DERIVS  "; TEST(d_ok);

        delete evalNoT;
        delete evalIGA;
    }
    else
    {
        gsInfo << "File " << fileName << " not found. Test Skipped." << "\n";
    }
}
