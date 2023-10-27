/** @file gsCompositeH_test.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsMSplines/gsMappedBasis.h>
#include <gsSmoothPatches/gsCompositeHBasis.h>
#include <gsSmoothPatches/gsCompositeIncrSmoothnessBasis.h>
#include <gsSmoothPatches/gsCompositeIncrSmoothnessBasis.hpp>
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
    try
    {
        gsReadFile<>(path, mp);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeHBasis<2,real_t> base_new(mp);

    gsMatrix<> result;

    gsInfo<<"Running time test for "<< testMatrix.cols() <<" points("<<base_new.nPatches()<<" patches).\n";

    clock.restart();
    for(size_t i=0;i<base_new.nPatches();i++)
    {
        mp.patch(i).basis().eval_into(testMatrix,result);
    }
    const double ctime = clock.stop();
    gsInfo << "Hsplines time: "<<ctime<<"\n";

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
    boxes(0,0) = 0.6;
    boxes(0,1) = 0.9;
    boxes(1,0) = 0.6;
    boxes(1,1) = 0.9;
    gsCompositeIncrSmoothnessGeom<2,real_t> geomRef(geom);
    geomRef.refine(0,boxes);
    geomRef.refine(0,boxes);
    geomRef.setBasis(0);
    geom.setBasis(0);
    unsigned nr=0;
    do
    {
        geom.eval_into(testMatrix,result);
        geomRef.eval_into(testMatrix,resultRef);
        passed = gsAllCloseAbsolute(result,resultRef,eps) && passed;
        geomRef.next();
        gsInfo << " >======= patch:" << nr++ << "\n";
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
    geomRef.uniformRefine();
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

bool test_UTPR_multiPatchConstructor(std::string path, bool plot, unsigned nsamples, int & plotNr, int minEVDistance)
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
    gsCompositeIncrSmoothnessGeom<2,real_t> geom(mp,minEVDistance);
    std::stringstream ss;
    ss<<"compositeSurf"<<plotNr++;
    geom.setBasis(0);
    if(plot)
    {
        gsStopwatch time;
        gsWriteParaview( geom.exportToPatches(), ss.str(), nsamples, true);
        const double timeNeeded = time.stop();
        gsInfo << "time: " << timeNeeded << "\n" ;
        ss<<"multipatch";
        gsWriteParaview( mp , ss.str(), nsamples,true);
        const double timeNeeded2 = time.stop();
        gsInfo << "time2: " << timeNeeded2-timeNeeded << "\n" ;
        gsInfo << "writing done: " <<ss.str()<< "\n" ;
    }
    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool test_UTPR_compositeHBasis(std::string path)
{
    gsInfo << "----------- gsCompositeHBasis -----------\n\n";
    gsHTensorBasis<2> * hbs = gsFileData<>( path ).getFirst< gsTHBSplineBasis<2> >().release();
    if(hbs==NULL)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsBoxTopology topol(2,1);
    topol.addAutoBoundaries();
    std::vector<gsHTensorBasis<2,real_t>*> bases;
    bases.push_back(hbs);
    gsCompositeHBasis<2,real_t> basis(bases,topol);
    gsInfo << "thb: " << hbs->size() << " composite: "  << basis.size();
    delete hbs;
    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool test_UTPR_compositeHBasis2(std::string path)
{
    gsInfo << "----------- gsCompositeHBasis2 -----------\n\n";
    gsHTensorBasis<2> * hbs = gsFileData<>( path ).getFirst< gsTHBSplineBasis<2> >().release();
    if(hbs==NULL)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    //hbs->uniformRefine();
    std::vector<gsHTensorBasis<2> *> bases;
    bases.push_back(hbs->clone().release());
    bases.push_back(hbs->clone().release());
    std::vector<boundaryInterface> interfaces;
    boundaryInterface binterface = boundaryInterface( patchSide(0,boundary::east),patchSide(1,boundary::west),
                                                      static_cast<short_t>(2));
    interfaces.push_back(binterface);
    std::vector<patchSide> boundaries;
    patchSide ps0_north(0,boundary::north);
    patchSide ps0_east(0,boundary::west);
    patchSide ps0_south(0,boundary::south);
    patchSide ps1_north(1,boundary::north);
    patchSide ps1_east(1,boundary::east);
    patchSide ps1_south(1,boundary::south);
    boundaries.push_back(ps0_north);
    boundaries.push_back(ps0_south);
    boundaries.push_back(ps0_east);
    boundaries.push_back(ps1_north);
    boundaries.push_back(ps1_south);
    boundaries.push_back(ps1_east);
    gsBoxTopology topol(2,2,boundaries,interfaces);
    gsCompositeHBasis<2,real_t> basis(bases,topol);
    basis.uniformRefine();
    //basis._knotsMatchNeighbours(0,patches);
    //basis.print(gsInfo);
    //basis.uniformRefine();
    gsInfo << "thb: " << hbs->size() << " composite: "  << basis.size();
    freeAll(bases);
    delete hbs;
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
    bool utwac      = false; // note: =true fails!
    bool utpr       = false;
    index_t nsamples    = 1000;

    gsMultiPatch<> mp;

    std::string fn("planar/two_squares.xml");

    gsCmdLine cmd("Composite basis tests.");
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)", fn);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("time", "Test evaluation time", timeTest);
    cmd.addInt("s","samples", "Number of samples to use for viewing", nsamples);
    cmd.addSwitch("utwac", "unit tests with automated checks", utwac);
    cmd.addSwitch("utpr", "unit tests that print results", utpr);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsReadFile<>(fn, mp);

    int plotNr = 0;
    real_t eps = math::sqrt(math::limits::epsilon())/100;

    gsVector<real_t,2> start,end;
    start << 0,0;
    end << 1,1;
    gsVector<unsigned,2> np;
    np << 10,10;
    gsMatrix<real_t> testMatrix=gsPointGrid<real_t>(start,end,np);

    std::vector<std::string> path_vec;
    //path_vec.push_back("planar/thbs_multipatch_01.xml");
    //path_vec.push_back("surfaces/thbs_multipatch_02.xml");
    path_vec.push_back("planar/multipatch_tunnel_thb.xml");

    if ( timeTest )
    {
        test_TT_compositeBasisTime("planar/thbs_multipatch_01.xml", testMatrix);
        return 0;
    }

    if( utwac )
        for(std::vector<std::string>::const_iterator it=path_vec.begin();it!=path_vec.end();++it)
        {
            passed = test_UTWAC_refinement(*it,testMatrix,eps) && passed;
            passed = test_UTWAC_uniformRefinement(*it,testMatrix,eps) && passed;
        }

    utpr=true;
    if( utpr )
        for(std::vector<std::string>::const_iterator it=path_vec.begin();it!=path_vec.end();++it)
        {
            passed = test_UTPR_compositeHBasis("basis_thbs_01.xml") && passed;
            passed = test_UTPR_compositeHBasis2("basis_thbs_01.xml") && passed;
            passed = test_UTPR_multiPatchConstructor(*it,plot,nsamples,plotNr,-1);
        }

    return passed ? 0 : 1;
}
