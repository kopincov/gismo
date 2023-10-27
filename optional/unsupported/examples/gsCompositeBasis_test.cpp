/** @file gsMappedBasis_test.h

    @brief File testing the gsMappedBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsMSplines/gsMappedBasis.h>
#include <gsSmoothPatches/gsCompositeBSplineBasis.h>
#include <gsSmoothPatches/gsCompositeIncrSmoothnessBasis.h>
#include <gsSmoothPatches/gsCompositeIncrSmoothnessGeom.h>
#include <gsSmoothPatches/gsCompositeUtils.h>

#include <iostream>
#include <cmath>


using std::flush;
using std::vector;
using std::pair;
using std::make_pair;
using namespace gismo;

//========================================== TIME TEST ========================================//

double test_TT_compositeBasisTime(std::string path, gsMatrix<real_t>& testMatrix)
{
    gsStopwatch clock;

    gsMultiPatch<> mp;
    gsReadFile<>(path, mp);

    gsCompositeBSplineBasis<2,real_t> base_new(mp);

    gsMatrix<> result;

    gsInfo<<"Running time test for "<< testMatrix.cols() <<" points("<<base_new.nPatches()<<" patches).\n";

    clock.restart();
    for(size_t i=0;i<base_new.nPatches();i++)
    {
        mp.patch(i).basis().eval_into(testMatrix,result);
    }
    const double ctime = clock.stop();
    gsInfo << "Bsplines time: "<<ctime<<"\n";

    clock.restart();
    for(size_t i=0;i<base_new.nPatches();i++)
    {
        base_new.eval_into(i,testMatrix,result);

        //base_new.active_into(i,u,result);
    }
    const double ctime2 = clock.stop();
    gsInfo << "Composite basis time: "<<ctime2<<"\n";

    return ctime2;
}

//========================================== UNIT TESTS WITH AUTOMATED CHECKS ========================================//

bool test_UTWAC_EvalAllDersSingle_into(std::string path, gsMatrix<real_t>& testMatrix, real_t eps)
{
    gsInfo << "----------- evalAllDersSingle_into test -----------\n\n";
    bool passed = true;
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }

    gsCompositeBSplineBasis<2,real_t> base_new(mp);
    gsMatrix<> result,result2;
    for(size_t i=0;i<base_new.nPatches();i++)
    {
        for(index_t j=0;j<base_new.size();j++)
        {
            base_new.evalSingle_into(i,j,testMatrix,result);
            base_new.evalAllDersSingle_into(i,j,testMatrix,0,result2);
            passed = gsAllCloseAbsolute(result,result2,eps) && passed;
        }
    }
    if(!passed)
        gsInfo << "evalAllDersSingle_into was not equal to evalSingle_into\n\n" ;
    else
        gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool test_UTWAC_EvalAllDers_into(std::string path, gsMatrix<real_t>& testMatrix, real_t eps)
{
    gsInfo << "----------- evalAllDers_into test -----------\n\n";
    bool passed = true;
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeBSplineBasis<2,real_t> base_new(mp);
    gsMatrix<> result;
    std::vector<gsMatrix<> > result2(1);
    for(size_t i=0;i<base_new.nPatches();i++)
    {
        base_new.eval_into(i,testMatrix,result);
        base_new.evalAllDers_into(i,testMatrix,0,result2);
        passed = gsAllCloseAbsolute(result,result2[0],eps) && passed;
    }
    if(!passed)
        gsInfo << "evalAllDers_into was not equal to eval_into\n\n" ;
    else
        gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool test_UTWAC_GlobalLocalEvaluation(std::string path, gsMatrix<real_t>& testMatrix,real_t eps)
{
    gsInfo<< "\n------------ Local vs Global ------------\n\n";
    bool passed = true;
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp,-1);
    gsMatrix<> resultLocal,resultGlobal;
    gsMultiPatch<> multipatch = geom.exportToPatches();
    gsInfo << "Geom testing...\n";
    geom.setBasis(0);
    int i = 0;
    do
    {
        gsInfo << "patch" << i <<":\n";
        geom.eval_into(testMatrix,resultGlobal);
        multipatch.patch(i).eval_into(testMatrix,resultLocal);
        passed = gsAllCloseAbsolute(resultGlobal,resultLocal,eps) && passed;
        gsInfo << "...eval passed\n";
        geom.deriv_into(testMatrix,resultGlobal);
        multipatch.patch(i).deriv_into(testMatrix,resultLocal);
        passed = gsAllCloseAbsolute(resultGlobal,resultLocal,eps) && passed;
        gsInfo << "...deriv passed\n";
        geom.deriv2_into(testMatrix,resultGlobal);
        multipatch.patch(i).deriv2_into(testMatrix,resultLocal);
        passed = gsAllCloseAbsolute(resultGlobal,resultLocal,eps) && passed;
        gsInfo << "...deriv2 passed\n";
        i++;
    }while(geom.next());
    if(!passed)
        gsInfo << "Global was not equal to Local\n\n" ;
    else
        gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return passed;
}

bool test_UTWAC_refinement(std::string path, gsMatrix<real_t>& testMatrix,real_t eps)
{
    gsInfo<< "\n------------ Refinement Test ------------\n\n";
    bool passed = true;
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp);
    gsMatrix<> result,resultRef;
    geom.setBasis(0);
    geom.setBasis(0);
    gsMatrix<> boxes(2,2);
    boxes(0,0) = 0.1;
    boxes(0,1) = 0.9;
    boxes(1,0) = 0.1;
    boxes(1,1) = 0.9;
    gsCompositeIncrSmoothnessGeom<2,real_t> geomRef(geom);
    geomRef.refine(0,boxes);
    geomRef.setBasis(0);
    geom.setBasis(0);
    do
    {
        geom.eval_into(testMatrix,result);
        geomRef.eval_into(testMatrix,resultRef);
        passed = gsAllCloseAbsolute(result,resultRef,eps) && passed;
        geomRef.next();
    }while(geom.next());
    if(!passed)
        gsInfo << "Refined was not equal to Local\n\n" ;
    else
        gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return passed;
}

bool test_UTWAC_uniformRefinement(std::string path, gsMatrix<real_t>& testMatrix,real_t eps)
{
    gsInfo<< "\n------------ uniform Refinement Test ------------\n\n";
    bool passed = true;
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp);
    gsMatrix<> result,resultRef;
    gsCompositeIncrSmoothnessGeom<2,real_t> geomRef(geom);
    geomRef.uniformRefine(2);
    geomRef.setBasis(0);
    geom.setBasis(0);
    do
    {
        geom.eval_into(testMatrix,result);
        geomRef.eval_into(testMatrix,resultRef);
        passed = gsAllCloseAbsolute(result,resultRef,eps) && passed;
        geomRef.next();
    }while(geom.next());
    if(!passed)
        gsInfo << "Refined was not equal to Local\n\n" ;
    else
        gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return passed;
}

//========================================== UNIT TESTS THAT ONLY PRINT RESULTS ========================================//

bool test_UTPR_compositeBasis(std::string path, gsMatrix<real_t>& testMatrix)
{
    gsInfo << "----------- gsMappedBasis -----------\n\n";
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeBSplineBasis<2,real_t> base_new(mp);
    patchCorner pc_special = patchCorner(0,boundary::northeast);
    gsInfo << "is special: " << base_new.isSpecialVertex(pc_special) << "\n" ;
    gsInfo << "gsMatrix<T> eval(int patch, int global_BF, const gsMatrix<TT> & u) const" << "\n" ;
    gsMatrix<index_t> result;
    for(size_t i=0;i<base_new.nPatches();i++)
    {
        base_new.active_into(i,testMatrix,result);
        gsInfo << "patch "<< i << "\n" <<result<<"\n" ;
    }
    gsMatrix<> result2;
    for(index_t i=0;i<base_new.size();i++)
    {
        base_new.evalSingle_into(0,i,testMatrix,result2);
        gsInfo << "globalbasisfunc "<< i << "\n" <<result2<<"\n" ;
    }
    gsInfo << "gsMatrix<> eval(int patch, const gsMatrix<TT> & u) const" << "\n" ;
    base_new.eval_into(0,testMatrix,result2);
    gsInfo<<result2<<"\n" ;
    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool test_UTPR_singleBasis(std::string path, gsMatrix<real_t>& testMatrix)
{
    gsInfo << "----------- gsSingleBasis test -----------\n\n";
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeBSplineBasis<2,real_t> base_new(mp);
    gsMappedSingleBasis<2,real_t> b = base_new.getMappedSingleBasis(0);
    gsMatrix<> result;
    int i = 0;
    do
    {
        b.eval_into(testMatrix,result);
        gsInfo << "patch "<< i++ << "\n" <<result<<"\n" ;
    }while(b.next());
    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool test_UTPR_compositeGeom(std::string path, gsMatrix<real_t>& testMatrix, bool plot,unsigned nsamples,int & plotNr)
{
    gsInfo << "----------- gsMappedSpline test -----------\n\n";
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp);
    geom.setBasis(0);
    gsMatrix<> result2;
    int i = 0;
    do
    {
        geom.eval_into(testMatrix,result2);
        gsInfo << "patch "<< i++ << "\n" <<result2<<"\n" ;
    }while(geom.next());
    std::stringstream ss;
    ss<<"compositeSurf"<<plotNr++;
    if(plot)
    {
        gsWriteParaview( geom.exportToPatches(), ss.str(), nsamples);
        gsInfo << "writing done: " <<ss.str()<< "\n" ;
    }
    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool test_UTPR_multiPatchConstructor(std::string path, bool plot, unsigned nsamples, int & plotNr)
{
    gsInfo<< "\n------------ multiPatchConstructor Test ------------\n\n";
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp,-1,4);
    std::stringstream ss;
    ss<<"compositeSurf"<<plotNr++;
    geom.setBasis(0);
    if(plot)
    {
        gsStopwatch time;
        gsWriteParaview( geom.exportToPatches(), ss.str(), nsamples, true, true);
        //gsWriteParaview( geom, ss.str(), nsamples, true, true);
        const double timeNeeded = time.stop();
        gsInfo << "time: " << timeNeeded << "\n" ;
        ss<<"multipatch";
        gsWriteParaview( mp , ss.str(), nsamples);
        const double timeNeeded2 = time.stop();
        gsInfo << "time2: " << timeNeeded2-timeNeeded << "\n" ;
        gsInfo << "writing done: " <<ss.str()<< "\n" ;
    }
    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool test_UTPR_differentWeights(std::string path, bool plot, unsigned nsamples, int & plotNr)
{
    gsInfo<< "\n------------ Interface-weights Test ------------\n\n";
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp,2);
    gsCompositeIncrSmoothnessBasis<2,real_t> & basis = dynamic_cast<gsCompositeIncrSmoothnessBasis<2,real_t>&>(geom.getCompBasis());
    patchSide ps(0,boundary::south);
    basis.setWeight(ps,2.0);
    ps=patchSide(0,boundary::east);
    basis.setWeight(ps,5.0);
    ps=patchSide(1,boundary::east);
    basis.setWeight(ps,3.0);
    basis.updateTopol();
    std::stringstream ss;
    ss<<"compositeSurf"<<plotNr++;
    geom.setBasis(0);
    if(plot)
    {
        gsStopwatch time;
        gsWriteParaview( geom.exportToPatches(), ss.str(), nsamples);
        const double timeNeeded = time.stop();
        gsInfo << "time: " << timeNeeded << "\n" ;
        ss<<"multipatch";
        gsWriteParaview( mp , ss.str(), nsamples);
        const double timeNeeded2 = time.stop();
        gsInfo << "time2: " << timeNeeded2-timeNeeded << "\n" ;
        gsInfo << "writing done: " <<ss.str()<< "\n" ;
    }
    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

//========================================== EXAMPLES OF USE CASES ========================================//

void do_L_shape(std::string path, bool plot, unsigned nsamples, int & plotNr)
{
    plot=true;
    gsInfo<< "\n------------ L-Shape Test ------------\n\n";
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return;
    }
    mp.patch(0).uniformRefine(1);
    mp.patch(1).uniformRefine(1);
    patchCorner pc1(0,boundary::southwest);
    patchCorner pc2(0,boundary::southeast);
    std::vector<patchCorner> cornerList;
    cornerList.push_back(pc1);
    cornerList.push_back(pc2);
    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp,cornerList,-1);
    geom.uniformRefine();
    geom.uniformRefine();
    geom.smoothEverything();
    gsMultiPatch<> mp2 = geom.exportToPatches();
    gsCompositeIncrSmoothnessGeom<2,real_t> geom2(mp2,cornerList,-1);
    std::stringstream ss;
    ss<<"compositeSurf"<<plotNr++;
    geom.setBasis(0);
    if(plot)
    {
        gsStopwatch time;
        gsWriteParaview( geom2.exportToPatches(), ss.str(), nsamples);
        const double timeNeeded = time.stop();
        gsInfo << "time: " << timeNeeded << "\n" ;
        ss<<"multipatch";
        gsWriteParaview( mp2 , ss.str(), nsamples,true);
        const double timeNeeded2 = time.stop();
        gsInfo << "time2: " << timeNeeded2-timeNeeded << "\n" ;
        gsInfo << "writing done: " <<ss.str()<< "\n" ;
    }
}

bool do_MTU_FIT(std::string path, bool plot, unsigned nsamples, int & plotNr)
{
    gsInfo<< "\n------------ MTUFIT Test ------------\n\n";
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }

    gsTensorBSpline<2> * patch;
    std::vector<real_t> knots;
    int deg;
    real_t new_knot;

    patch= dynamic_cast<gsTensorBSpline<2> *>(&(mp.patch(0)));
    new_knot = 1.320;
    knots.insert(knots.end(),patch->basis().component(0).knots().size(),new_knot);
    deg = patch->basis().component(0).knots().degree();
    patch->basis().component(0).knots()=gsKnotVector<>(give(knots), deg);

    patch = dynamic_cast<gsTensorBSpline<2> *>(&(mp.patch(1)));
    new_knot = .70;
    knots.insert(knots.end(),patch->basis().component(0).knots().size(), new_knot);
    deg = patch->basis().component(0).knots().degree();
    patch->basis().component(0).knots()=gsKnotVector<>(knots, deg);

    patch = dynamic_cast<gsTensorBSpline<2> *>(&(mp.patch(2)));
    new_knot = 1.476;
    knots.insert(knots.end(),patch->basis().component(0).knots().size(), new_knot);
    deg = patch->basis().component(0).knots().degree();
    patch->basis().component(0).knots()=gsKnotVector<>(give(knots), deg);

    patch = dynamic_cast<gsTensorBSpline<2> *>(&(mp.patch(3)));
    new_knot = .013;
    knots.insert(knots.end(),patch->basis().component(0).knots().size(), new_knot);
    deg = patch->basis().component(0).knots().degree();
    patch->basis().component(0).knots()=gsKnotVector<>(give(knots), deg);

    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp,3);
    //geom.uniformRefine(2);
    std::stringstream ss;
    ss<<"compositeSurf"<<plotNr++;
    geom.setBasis(0);
    if(plot)
    {
        gsStopwatch time;
        gsWriteParaview( geom.exportToPatches(), ss.str(), nsamples*100);
        const double timeNeeded = time.stop();
        gsInfo << "time: " << timeNeeded << "\n" ;
        ss<<"multipatch";
        gsWriteParaview( mp , ss.str(), nsamples*100);
        const double timeNeeded2 = time.stop();
        gsInfo << "time2: " << timeNeeded2-timeNeeded << "\n" ;
        gsInfo << "writing done: " <<ss.str()<< "\n" ;
    }
    gsFileData<> fd;
    gsMultiPatch<> mp2 = geom.exportToPatches() ;
    fd<< mp2;
    fd.dump("makeMultipatch_output");

    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool smoothenJakaExample()
{
    gsInfo<< "\n------------ smoothenJaka Test ------------\n\n";
    std::string path = "Jaka/2Patch.xml";
    unsigned nsamples = 100;
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }

    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp);
    std::stringstream ss;
    ss<<"JakaSurf";
    gsWriteParaview( geom.exportToPatches(), ss.str(), nsamples*100);
    gsFileData<> fd;
    gsMultiPatch<> mp2 = geom.exportToPatches() ;
    fd<< mp2;
    fd.dump("JakaSurf");

    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool periodicJakaExample()
{
//    gsKnotVector<> kv1(0, 1, 1, 2, 1);
//    gsKnotVector<> kv2(0, 1, 1, 2, 1);

    gsKnotVector<> kv1(0, 1, 1, 3, 1);
    gsKnotVector<> kv2(0, 1, 1, 3, 1);

    gsTensorBSplineBasis<2> tensorBasis(kv1, kv2);

    gsBoxTopology topology(2, 1);
    topology.addInterface(0, boundary::left, 0, boundary::right);
    topology.addInterface(0, boundary::up, 0, boundary::down);
    topology.addAutoBoundaries();
    std::vector< gsTensorBSplineBasis<2,real_t>* > vectorOfBasis;
    vectorOfBasis.push_back(&tensorBasis);

    gsCompositeBSplineBasis<2, real_t> periodicBasis(vectorOfBasis, topology);

    const gsBasis<>& basis = periodicBasis.getBase(0);

    gsInfo << "basis: " << basis << "\n";
    const gsSparseMatrix<> local2Global = periodicBasis.getMapper().asMatrix();

    gsInfo << "local2Global: \n"
              << local2Global.toDense() << "\n"
              << "sizeof: " << local2Global.rows() << " x " << local2Global.cols() << "\n";

    return true;
}

bool printSingleBasisFunc(std::string path,std::vector<int> basisFuncs,bool smooth, bool plot, unsigned nsamples, int &plotNr)
{
    gsInfo<< "\n------------ printSingleBasisFuncs Test ------------\n\n";
    if(!plot)
        return true;
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp,smooth?-1:0,4);
    gsMatrix<> coefs = geom.coefs();
    coefs.conservativeResize(coefs.rows(),3);
    coefs.col(2).setZero();
    for(size_t i = 0;i<basisFuncs.size();i++)
        if(basisFuncs[i]<coefs.rows())
            coefs(basisFuncs[i],2)=1;
    geom.coefs()=coefs;

    std::stringstream ss;
    ss<<"compositeSurf"<<plotNr++;
    geom.setBasis(0);
    gsStopwatch time;
    gsWriteParaview( geom.exportToPatches(), ss.str(), nsamples, true);
    const double timeNeeded = time.stop();
    gsInfo << "time: " << timeNeeded << "\n" ;
    gsInfo << "writing done: " <<ss.str()<< "\n" ;
    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

//========================================== MAIN FUNCTION ========================================//

int main(int argc, char *argv[])
{
    bool passed = true;

    bool plot       = false;
    bool timeTest   = false;
    bool utwac      = true;
    bool utpr       = false;
    index_t nsamples    = 1000;

    std::string fn ("planar/two_squares.xml");

    gsCmdLine cmd("Composite basis tests.");
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)", fn);
    cmd.addSwitch("t", "time", "Test evaluation time", timeTest);
    cmd.addSwitch("w", "utwac", "unit tests with automated checks", utwac);
    cmd.addSwitch("r", "utpr", "unit tests that print results",utpr);
    cmd.addInt("s","samples","Number of samples to use for viewing",nsamples);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    int plotNr = 0;
    real_t eps = math::sqrt(math::limits::epsilon())/100;

    gsVector<real_t,2> start,end;
    start << 0,0;
    end << 1,1;
    gsVector<unsigned,2> np;
    np << 10,10;
    gsMatrix<real_t> testMatrix=gsPointGrid<real_t>(start,end,np);

    std::vector<std::string> path_vec;
    path_vec.push_back("planar/four_squares.xml");
    //path_vec.push_back("planar/multipatch_triangle4.xml");
    //path_vec.push_back("planar/four_squares_as_one.xml");
    //path_vec.push_back("planar/multipatch_tunnel.xml");
//    path_vec.push_back("surfaces/multipatch_triangle2d.xml");
//    path_vec.push_back("planar/multipatch_2patch.xml");
//    path_vec.push_back("planar/two_squares.xml");
//    path_vec.push_back("surfaces/multipatch_2patch.xml");
    path_vec.push_back("planar/multipatch_triangle.xml");
//    path_vec.push_back("surfaces/multipatch_triangle.xml");
    //path_vec.push_back("planar/multipatch_triangle5_deg3.xml");
//    path_vec.push_back("surfaces/three_quads.xml");
//    path_vec.push_back("surfaces/multipatch_AirPassage.xml");
    //path_vec.push_back("planar/four_squares4.xml");

    if ( timeTest )
    {
        test_TT_compositeBasisTime(path_vec[0], testMatrix);
        return 0;
    }

    if( utwac )
        for(std::vector<std::string>::const_iterator it=path_vec.begin();it!=path_vec.end();++it)
        {
            passed = test_UTWAC_EvalAllDersSingle_into(*it,testMatrix,eps) && passed;
            passed = test_UTWAC_EvalAllDers_into(*it,testMatrix,eps) && passed;
            passed = test_UTWAC_GlobalLocalEvaluation(*it,testMatrix,eps) && passed;
            passed = test_UTWAC_refinement(*it,testMatrix,eps) && passed;
            passed = test_UTWAC_uniformRefinement(*it,testMatrix,eps) && passed;
        }
    if( utpr )
        for(std::vector<std::string>::const_iterator it=path_vec.begin();it!=path_vec.end();++it)
        {
            passed = test_UTPR_compositeBasis(*it,testMatrix) && passed;
            passed = test_UTPR_singleBasis(*it,testMatrix) && passed;
            passed = test_UTPR_compositeGeom(*it,testMatrix,plot,nsamples,plotNr) && passed;
            passed = test_UTPR_multiPatchConstructor(*it,plot,nsamples,plotNr) && passed;
            passed = test_UTPR_differentWeights(*it,plot,nsamples,plotNr) && passed;
        }
    //    std::vector<int> basisFuncs;
    //    basisFuncs.push_back(86);
    //    basisFuncs.push_back(98);
    //    basisFuncs.push_back(30);
    //    basisFuncs.push_back(59);
    //    basisFuncs.push_back(110);
    //    basisFuncs.push_back(120);
    //    passed = printSingleBasisFunc("planar/multipatch_triangle.xml",basisFuncs,true,true,nsamples,plotNr);
    //    basisFuncs.resize(0);
        //    basisFuncs.push_back(105);
        //    basisFuncs.push_back(19);
        //    basisFuncs.push_back(65);
        //    basisFuncs.push_back(32);
        //    basisFuncs.push_back(110);
        //    basisFuncs.push_back(120);
        //    passed = printSingleBasisFunc("planar/multipatch_triangle.xml",basisFuncs,false,true,nsamples,plotNr);

    //    passed = do_MTU_FIT("MTU/crescendo2.xml",plot,nsamples,plotNr) && passed;
    //    do_L_shape("planar/lshape2d_2patches2.xml",plot,nsamples,plotNr);
    //    smoothenJakaExample();
    //    periodicJakaExample();
    return passed ? 0 : 1;
}
