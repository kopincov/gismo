/** @file gsCompositeMeasure_test.h

    @brief File testing the gsCompositeMeasure class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/
#include <gismo.h>

#include <gsSmoothPatches/gsCompositeAssemblerUtils.h>
#include <gsSmoothPatches/gsCompositeUtils.h>

//#include "gsMSplines/gsMappedBasis.h"
//#include "gsSmoothPatches/gsCompositeBSplineBasis.h"
//#include "gsMSplines/gsMappedSpline.h"
#include "gsSmoothPatches/gsCompositeIncrSmoothnessGeom.h"

//#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsRecipeAssemblerAdaptive/gsAdaptiveSolver.h>
//#include <gsRecipeAssemblerAdaptive/gsCompositeBasisSpaceRefiners.h>
//#include <gsRecipeAssemblerAdaptive/gsTensorBasisSpaceRefiners.h>

//#include <gsMSplines/gsWeightMapper.h>


;
using std::flush;
using std::vector;
using std::pair;
using std::make_pair;
using std::ios;
using namespace gismo;

//========================================== Child Classes for the adaptive Solver ========================================//

class gsFittingErrorEstimator : public gsErrorEstimator
{
    gsQualityMeasureWeights        &m_weights;
    gsMappedSpline<2,real_t>      &m_compGeom;

public:
    gsFittingErrorEstimator(gsQualityMeasureWeights &weights,gsMappedSpline<2,real_t>& compGeom) :
        gsErrorEstimator(), m_weights(weights),m_compGeom(compGeom)
    {
        size_t fullSize = 0;
        for(size_t i = 0;i<m_weights.m_points.size();++i)
            fullSize+=m_weights.m_pars[i].cols();
        m_local.resize(1,fullSize);
        m_total.resize(1,1);
    }

    virtual void computeEstimate(
        const std::vector<gsPhysicalSpace*>  & /*space*/,
        const std::vector<gsMatrix<real_t> > & /*coefs*/)
    {
        calculateLocalErrors();
        calculateTotalErrors();
        gsInfo << "total error: " << m_total << "\n";
    }

private:
    void calculateLocalErrors()
    {
        gsMatrix<real_t> res;
        index_t pos = 0;
        for(size_t i = 0;i<m_weights.m_points.size();++i)
        {
            m_compGeom.setBasis(i);
            m_compGeom.eval_into(m_weights.m_pars[i],res);
            gsMatrix<real_t> tmp = (res-m_weights.m_points[i]);
            for(index_t point = 0;point<m_weights.m_points[i].cols();++point)
            {
                m_local(0,pos)=tmp.col(point).lpNorm<1>();
                pos++;
            }
        }
    }

    void calculateTotalErrors()
    {
        m_total=m_local.rowwise().sum();
    }
};

class gsHTensorBasisRefinerAdaptiveOptimizer : public gsHTensorBasisRefiner
{
protected:
    gsMultiPatch<>                          &m_domain;
    gsCompositeIncrSmoothnessGeom<2,real_t> &m_compGeom;
    gsQualityMeasureWeights                 &m_weights;
public:
    gsHTensorBasisRefinerAdaptiveOptimizer(
            gsMultiPatch<real_t> &domain,
            gsCompositeIncrSmoothnessGeom<2,real_t>  &compGeom,
            gsQualityMeasureWeights &weights
            )
        : m_domain(domain),m_compGeom(compGeom),m_weights(weights)
    {
        updateDomain(m_compGeom.coefs());
    }

    virtual std::vector<gsPhysicalSpace*> getSpaces()const
    {
        gsPhysicalSpaceScalar physicalSpace = gsPhysicalSpaceScalar(m_compGeom.getCompBasis().getBases(),m_domain,
                                                                    NO_TRANSFORMATION,m_compGeom.getCompBasis().getMapper());
        std::vector<gsPhysicalSpaceScalar*> phySpaces;
        phySpaces.push_back(&physicalSpace);
        phySpaces.push_back(&physicalSpace);
        gsPhysicalSpace* physicalSpaceVector = new gsPhysicalSpaceVector(phySpaces);
        std::vector<gsPhysicalSpace* > vectorSpaces;
        vectorSpaces.push_back(physicalSpaceVector);
        return vectorSpaces;
    }

    virtual void updateSpaces (const gsMatrix<real_t>& marked)
    {
        std::vector<std::vector<index_t> > boxes;
        std::vector<index_t> pboxes;
        boxes.push_back(pboxes);
        index_t localp = 0, patch = 0;
        gsVector<unsigned> ext(2);
        ext.setOnes();
        for(index_t p = 0;p<marked.cols();++p)
        {
            std::vector<index_t> box;
            if(marked(0,p)>0)
            {
                gsVector<real_t> point = m_weights.m_pars[patch].col(localp);
                const gsHTensorBasis<2, real_t>* hTensorBasis =
                        dynamic_cast< const gsHTensorBasis<2,real_t>* > (&m_compGeom.getCompBasis().getBase(patch));
                getBoxPerPoint<2>(*hTensorBasis,point,box);
                addExtension(*hTensorBasis,ext,box);
                liftLevel(*hTensorBasis,1,box);
                addBox(boxes[patch],box);
            }
            localp++;
            if(m_weights.m_pars[patch].cols()<=localp)
            {
                localp=0;
                patch++;
                boxes.push_back(pboxes);
            }
        }
        boxes.pop_back();
        m_compGeom.refineElements(boxes);
        m_domain=m_compGeom.exportToPatches();
    }

    void addBox(std::vector<index_t>& boxes,std::vector<index_t>& box)
    {
        boxes.insert(boxes.end(), box.begin(), box.end());
    }

    void updateDomain(const gsMatrix<real_t> &newCoefs)
    {
        m_compGeom.coefs()=newCoefs;
        m_domain=m_compGeom.exportToPatches();
    }

    gsMatrix<real_t> & getCurCoefs()
    {
        return m_compGeom.coefs();
    }

    gsMultiPatch<> & getDomain()
    {
        return m_domain;
    }

};

class gsAdaptiveSolverOptimizeGeometry : public gsAdaptiveSolver
{
private:
    gsMultiPatch<real_t>                    m_domain;
    gsHTensorBasisRefinerAdaptiveOptimizer *m_spaceRefinerAO;
    gsRecipeAssemblerQMOpt2D               *m_QMOptAssembler;
    gsRecipeAssemblerQMOpt2D               *m_areaCalculator;
    gsQualityMeasureWeights                *m_areaWeights;
    gsQualityMeasureWeights                &m_weights;
    bool                                    m_fixBoundaries;
    bool                                    m_isFittingProblem;

    gsMatrix<real_t>                        m_curCoefs;
    real_t                                  m_curVal;
    size_t                                  m_dampingSteps;

    std::vector<real_t>                     m_times;
    std::vector<size_t>                     m_dofs;


public:
    gsAdaptiveSolverOptimizeGeometry(gsCompositeIncrSmoothnessGeom<2,real_t> * compGeom,gsQualityMeasureWeights& weights,bool fixBoundaries)
        : gsAdaptiveSolver(),m_areaCalculator(NULL),m_weights(weights),m_fixBoundaries(fixBoundaries)
    {
        real_t maxError = 0.0001;
        m_domain=compGeom->exportToPatches();
        gsMappedBasis<2,real_t>& compBasis = compGeom->getCompBasis();
        std::vector<gsBasis<real_t>*> bases;
        for(size_t i =0;i<compBasis.nPatches();++i)
            bases.push_back(&compBasis.getBase(i));
        m_isFittingProblem = m_weights.m_points.size()!=0 && m_weights.m_points[0].size()!=0;
        gsEigenBiCGSTABILUT<real_t> * eigensolver = new gsEigenBiCGSTABILUT<real_t>;

        m_QMOptAssembler = new gsRecipeAssemblerQMOpt2D(m_domain,bases,compBasis.getMapper(),weights,true,true,true,false,false,m_isFittingProblem);
        m_assembler      = m_QMOptAssembler;
        m_solver         = eigensolver;
        m_solverEli      = m_solver;

        m_spaceRefinerAO = new gsHTensorBasisRefinerAdaptiveOptimizer(m_domain,*compGeom,weights);
        m_refiner        = m_spaceRefinerAO;
        if(m_isFittingProblem)
        {
            m_estimator      = new gsFittingErrorEstimator(weights,*compGeom);
            m_marker         = new gsMarkerAbsoluteThreshold(maxError);
        }

        m_criteria       = new gsStopCriteriaIterationOnly(5);

        eigensolver->setTolerance(1e-10);
        eigensolver->setMaxIterations(10000);

        m_dampingSteps = 0;
        m_curVal=std::numeric_limits<real_t>::max();
        m_curCoefs=m_spaceRefinerAO->getCurCoefs();
        updateGeometry();
    }

    ~gsAdaptiveSolverOptimizeGeometry()
    {
        delete m_assembler;
        delete m_solver;
        delete m_estimator;
        delete m_refiner;
        delete m_marker;
        delete m_criteria;
    }

    virtual void initStep()
    {
        m_assembler->reset();
        if(m_marker)
            m_marker->reset();
        if(m_estimator)
            m_estimator->reset();
        m_curCoefs=m_spaceRefinerAO->getCurCoefs();
    }

    void adaptiveSolve ()
    {
        m_criteria->reset();
        gsMatrix<real_t> vals;
        this->initStep();
        this->assemble();
        if(!m_QMOptAssembler->getQualityMeasureVals(vals))
            return;
        size_t steps=0;
        real_t dampingFactor=1;
        if(0<= vals(0,0) && vals(0,0)<m_curVal)
            m_curVal=vals(0,0);
        while(true)
        {
            if(!this->solve())
            {
                gsWarn<<"Exiting...\n";
                m_internalErrorFlag=true;
                return;
            }
            updateGeometry(true,dampingFactor);
            if(m_isFittingProblem)
            {
                this->errorEstimate();
                this->mark();
                this->refine();
            }
            this->initStep();
            this->assemble();
            if(!m_QMOptAssembler->getQualityMeasureVals(vals))
                return;
            if(steps<m_dampingSteps && vals(0,0)>m_curVal)
            {
                gsInfo << "step " << steps+1 << " damping. \n";
                steps++;
                dampingFactor/=2;
            }
            else
            {
                steps=0;
                dampingFactor=1;
                m_curCoefs=m_spaceRefinerAO->getCurCoefs();
                m_curVal=vals(0,0);
                gsInfo << "curVal " << m_curVal << " . \n";
                if(m_criteria->stop())
                    break;
                if(steps>=m_dampingSteps&&m_dampingSteps!=0)
                {
                    // damping failed
                    gsInfo << "\nDamping Failed: Aborting...\n";
                    break;
                }
            }
        }
    }

    void updateGeometry(bool newCoefs=false,real_t dampingFactor=1)
    {
        if(newCoefs)
        {
            this->m_solCoefs[0].resizeLike(m_curCoefs);
            this->m_spaceRefinerAO->updateDomain(this->m_solCoefs[0]*dampingFactor+m_curCoefs);
        }
        std::vector<gsPhysicalSpace* > spaces = m_spaceRefinerAO->getSpaces();
        if(m_weights.m_weightArea)
        {
            m_areaWeights = new gsQualityMeasureWeights (spaces.size());
            m_areaCalculator=new gsRecipeAssemblerQMOpt2D(m_spaceRefinerAO->getDomain(),spaces,*m_areaWeights,
                                                          false,false,false,true,false,false);
        }
    }

    void assemble()
    {
        gsStopwatch time;
        if(m_weights.m_weightArea)
            calculateAreas();
        if(m_fixBoundaries)
            m_QMOptAssembler->fixBoundaries();
        gsAdaptiveSolver::assemble();
        m_times.push_back(time.stop());
    }

    real_t getCurVal()
    {
        return m_curVal;
    }

    void calculateAreas()
    {
        m_areaCalculator->assemble();
        m_areaCalculator->getAreas(m_weights.m_areas);
    }

    void logProgress ()
    {
        m_dofs.push_back(m_assembler->getSpace()[0]->getMapper()->getNrOfTargets());
    }

    gsVector<real_t> getTimes() const
    {
        gsAsConstVector<real_t> times(m_times);
        return times;
    }

    gsVector<size_t> getDofs() const
    {
        gsAsConstVector<size_t> dofs(m_dofs);
        return dofs;
    }

};

void printMPToParaview(const gsMultiPatch<real_t>& mp,std::string name,int nr)
{
    for(size_t i = 0;i<mp.nPatches();++i)
    {
        gsMatrix<real_t> pts,u;
        u.resize(2,15);
        for(unsigned j = 0;j<10;++j)
        {
            u(0,j)=j/9.;
            u(1,j)=0;
            if(j<5)
            {
                u(1,j+10)=j/4.;
                u(0,j+10)=0;
            }
        }
        mp.patch(i).eval_into(u,pts);
        gsWriteParaviewPoints(pts,name+util::to_string(i)+"pts");
    }
    gsWriteParaview(mp,name,nr,true);
}

void checkMatrix(gsMatrix<real_t>& mat,unsigned patch)
{
    for(int i =0;i<mat.rows();++i)
    {
        for(int j=0;j<mat.cols();++j)
        {
            if(mat(i,j)<0 || mat(i,j)>1)
                gsInfo << "wrong parameter " << mat(i,j) <<" in patch " <<patch<<" at ("<<i<<","<<j<<")\n";
        }
    }
}

void checkPhysMatrix(gsMatrix<real_t>& mat,unsigned patch)
{
    for(int i =0;i<mat.rows();++i)
    {
        for(int j=0;j<mat.cols();++j)
        {
            if(mat(i,j)<-10 || mat(i,j)>10)
                gsInfo << "wrong point " << mat(i,j) <<" in patch " <<patch<<" at ("<<i<<","<<j<<")\n";
        }
    }
}

bool test_measureSolvingHatPatch(unsigned patch)
{
    std::vector<gsGeometry<real_t> *> patches;
    const size_t degree=2;
    gsKnotVector<real_t>kv(0,1,0,degree+1,1,degree);
    gsMatrix<real_t> coefs;
    const size_t coefsInDir=kv.size()-degree-1;
    coefs.setZero(coefsInDir*coefsInDir,2);
    patches.push_back(new gsTensorBSpline<2,real_t>(kv,kv,coefs));
    patches[0]->uniformRefine(8);

    gsMultiPatch<real_t> target(patches);
    target.addAutoBoundaries();
    gsQualityMeasureWeights weights(1);
    weights.m_length=0.000001;
    weights.m_approxPoints=1;
    gsFileData<real_t> dataInsidePatch( "Antonella/HAT_Final/SamplingHATpatches.xml" );
    gsFileData<real_t> dataBoundary( "Antonella/HAT_Final/BoundarySamplingHAT.xml" );
    gsMatrix<real_t>::uPtr parsInside, parsBoundary;
    gsMatrix<real_t>::uPtr pointsInside, pointsBoundary;
    const int id_par=2*patch;
    const int id_point=2*patch+1;
    parsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_par);
    pointsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_point);
    parsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_par);
    pointsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_point);
    if(parsInside && parsBoundary)
    {
        const index_t totalCols = parsInside->cols()+parsBoundary->cols();
        const index_t totalRows = parsInside->rows();
        weights.m_pars[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_points[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_pars[0].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
        weights.m_pars[0].rightCols(parsBoundary->cols())=parsBoundary->rightCols(parsBoundary->cols());
        weights.m_points[0].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
        weights.m_points[0].rightCols(pointsBoundary->cols())=pointsBoundary->rightCols(pointsBoundary->cols());
    }
    else if(parsInside)
    {
        weights.m_pars[0]=*parsInside;
        weights.m_points[0]=*pointsInside;
    }
    else if(parsBoundary)
    {
        weights.m_pars[0]=*parsBoundary;
        weights.m_points[0]=*pointsBoundary;
    }
    gsInfo << "\n";

    // real_t time;
    // unsigned dof;
    gsCompositeIncrSmoothnessGeom<2,real_t> compGeom(target,-1);

    gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,false);
    adaptiveSolver.adaptiveSolve();

    gsWriteParaviewPoints(weights.m_points[0],"pointsOfPatch"+util::to_string(patch));
    gsMultiPatch<real_t> mp2 = compGeom.exportToPatches();
    printMPToParaview(mp2,"firstOptSecond"+util::to_string(patch),10000);
    return true;
}

bool test_measureSolvingHat()
{
    std::string path="Antonella/HAT_Final/TemplateHat_patches2.xml";
    gsMultiPatch<> topolGeom;
    try
    {
        gsReadFile<>(path, topolGeom);
    }
    catch (std::runtime_error&)
    {
        return false;
    }

    if(topolGeom.empty())
    {
        gsInfo << "No gsMultiPatch in file." << "\n" ;
        return false;
    }

    std::vector<gsGeometry<real_t> *> patches;
    const size_t degree=2;
    for(size_t i = 0;i<topolGeom.nPatches();++i)
    {
        gsKnotVector<real_t>kv(0,1,0,degree+1,1,degree);
        gsMatrix<real_t> coefs;
        const size_t coefsInDir=kv.size()-degree-1;
        coefs.setZero(coefsInDir*coefsInDir,2);
        //patches.push_back(new gsTensorBSpline<2,real_t>(kv,kv,coefs));
        gsTensorBSpline<2,real_t> * patch = new gsTensorBSpline<2,real_t>(kv,kv,coefs);
        patch->uniformRefine(1);
        patches.push_back(new gsTHBSpline<2,real_t>(*patch));
    }
    gsMultiPatch<real_t> target(patches,topolGeom.boundaries(),topolGeom.interfaces());
    gsQualityMeasureWeights weights(target.nPatches());
    weights.m_uniformity=0.0001;
    weights.m_approxPoints=1;
    gsFileData<real_t> dataInsidePatch( "Antonella/HAT_Final/HATPuntiDENTRO.xml" );
    gsFileData<real_t> dataBoundaryHoles( "Antonella/HAT_Final/CappelloBordiBuchi.xml" );
    gsFileData<real_t> dataBoundary( "Antonella/HAT_Final/CappelloBordi.xml" );
    gsMatrix<real_t>::uPtr parsInside,  parsBoundary,  parsBoundaryHoles;
    gsMatrix<real_t>::uPtr pointsInside,pointsBoundary,pointsBoundaryHoles;
    for(size_t i = 0;i<target.nPatches();++i)
    {
        const int id_par=2*i;
        const int id_point=2*i+1;
        parsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_par);
        pointsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_point);
        parsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_par);
        pointsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_point);
        parsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_par);
        pointsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_point);
        if(parsInside && parsBoundaryHoles)
        {
            const index_t totalCols = parsInside->cols()+parsBoundaryHoles->cols();
            const index_t totalRows = parsInside->rows();
            weights.m_pars[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_points[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_pars[i].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
            weights.m_pars[i].rightCols(parsBoundaryHoles->cols())=parsBoundaryHoles->rightCols(parsBoundaryHoles->cols());
            weights.m_points[i].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
            weights.m_points[i].rightCols(pointsBoundaryHoles->cols())=pointsBoundaryHoles->rightCols(pointsBoundaryHoles->cols());
        }
        else if(parsInside && parsBoundary)
        {
            const index_t totalCols = parsInside->cols()+parsBoundary->cols();
            const index_t totalRows = parsInside->rows();
            weights.m_pars[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_points[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_pars[i].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
            weights.m_pars[i].rightCols(parsBoundary->cols())=parsBoundary->rightCols(parsBoundary->cols());
            weights.m_points[i].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
            weights.m_points[i].rightCols(pointsBoundary->cols())=pointsBoundary->rightCols(pointsBoundary->cols());
        }
        else if(parsInside)
        {
            weights.m_pars[i]=*parsInside;
            weights.m_points[i]=*pointsInside;
        }
        else if(parsBoundary)
        {
            weights.m_pars[i]=*parsBoundary;
            weights.m_points[i]=*pointsBoundary;
        }
        else if(parsBoundaryHoles)
        {
            weights.m_pars[i]=*parsBoundaryHoles;
            weights.m_points[i]=*pointsBoundaryHoles;
        }
        gsInfo << "\n";
    }
    // real_t time;
    // unsigned dof;
    gsCompositeIncrSmoothnessGeom<2,real_t> compGeom(target,-1);

    gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,false);
    adaptiveSolver.adaptiveSolve();

    for(size_t i =0;i<target.nPatches();++i)
    {
        if(weights.m_pars[i].rows()!=0 && weights.m_pars[i].cols()!=0)
            gsWriteParaviewPoints(weights.m_points[i],"pointsOfPatch"+util::to_string(i));
    }
    gsMultiPatch<real_t> mp2 = compGeom.exportToPatches();
    printMPToParaview(mp2,"firstOptSecond",10000);
    return true;
}

bool test_measureSolvingQuestionMarkPatch(unsigned patch)
{
    gsInfo << "solve for patch: " << patch << "\n";
    std::vector<gsGeometry<real_t> *> patches;
    const size_t degree=2;
    gsKnotVector<real_t>kv(0,1,0,degree+1,1,degree);
    gsMatrix<real_t> coefs;
    const size_t coefsInDir=kv.size()-degree-1;
    coefs.setZero(coefsInDir*coefsInDir,2);
    patches.push_back(new gsTensorBSpline<2,real_t>(kv,kv,coefs));
    patches[0]->uniformRefine(8);

    gsMultiPatch<real_t> target(patches);
    target.addAutoBoundaries();
    gsQualityMeasureWeights weights(1);
    //weights.m_length=0.000001;
    weights.m_uniformity=0.00001;
    weights.m_approxPoints=1;
    gsFileData<real_t> dataInsidePatch( "Antonella/QuestionMark/InnerPoints.xml" );
    gsFileData<real_t> dataBoundaryHoles( "Antonella/QuestionMark/OnlyHOLESpoints.xml" );
    gsFileData<real_t> dataBoundary( "Antonella/QuestionMark/Boundary.xml" );
    gsMatrix<real_t>::uPtr parsInside,  parsBoundary,  parsBoundaryHoles;
    gsMatrix<real_t>::uPtr pointsInside,pointsBoundary,pointsBoundaryHoles;
    const int id_par=2*patch;
    const int id_point=2*patch+1;
    parsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_par);
    pointsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_point);
    parsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_par);
    pointsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_point);
    parsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_par);
    pointsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_point);
    if(parsInside && parsBoundaryHoles)
    {
        const index_t totalCols = parsInside->cols()+parsBoundaryHoles->cols();
        const index_t totalRows = parsInside->rows();
        weights.m_pars[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_points[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_pars[0].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
        weights.m_pars[0].rightCols(parsBoundaryHoles->cols())=parsBoundaryHoles->rightCols(parsBoundaryHoles->cols());
        weights.m_points[0].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
        weights.m_points[0].rightCols(pointsBoundaryHoles->cols())=pointsBoundaryHoles->rightCols(pointsBoundaryHoles->cols());
    }
    else if(parsInside && parsBoundary)
    {
        const index_t totalCols = parsInside->cols()+parsBoundary->cols();
        const index_t totalRows = parsInside->rows();
        weights.m_pars[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_points[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_pars[0].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
        weights.m_pars[0].rightCols(parsBoundary->cols())=parsBoundary->rightCols(parsBoundary->cols());
        weights.m_points[0].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
        weights.m_points[0].rightCols(pointsBoundary->cols())=pointsBoundary->rightCols(pointsBoundary->cols());
    }
    else if(parsInside)
    {
        weights.m_pars[0]=*parsInside;
        weights.m_points[0]=*pointsInside;
    }
    else if(parsBoundary)
    {
        weights.m_pars[0]=*parsBoundary;
        weights.m_points[0]=*pointsBoundary;
    }
    else if(parsBoundaryHoles)
    {
        weights.m_pars[0]=*parsBoundaryHoles;
        weights.m_points[0]=*pointsBoundaryHoles;
    }
    gsInfo << "\n";

    // real_t time;
    // unsigned dof;
    gsCompositeIncrSmoothnessGeom<2,real_t> compGeom(target,-1);

    gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,false);
    adaptiveSolver.adaptiveSolve();

    gsWriteParaviewPoints(weights.m_points[0],"pointsOfPatch"+util::to_string(patch));
    gsMultiPatch<real_t> mp2 = compGeom.exportToPatches();
    printMPToParaview(mp2,"firstOptSecond"+util::to_string(patch),10000);
    return true;
}

bool test_measureSolvingQuestionMark()
{
    std::string path="Antonella/QuestionMark/2question_mark_template.xml";
    gsMultiPatch<> topolGeom;
    try
    {
        gsReadFile<>(path, topolGeom);
    }
    catch (std::runtime_error&)
    {
        return false;
    }

    if(topolGeom.empty())
    {
        gsInfo << "No gsMultiPatch in file." << "\n" ;
        return false;
    }

    std::vector<gsGeometry<real_t> *> patches;
    const size_t degree=2;
    for(size_t i = 0;i<topolGeom.nPatches();++i)
    {
        gsKnotVector<real_t>kv(0,1,0,degree+1,1,degree);
        gsMatrix<real_t> coefs;
        const size_t coefsInDir=kv.size()-degree-1;
        coefs.setZero(coefsInDir*coefsInDir,2);
        patches.push_back(new gsTensorBSpline<2,real_t>(kv,kv,coefs));
        patches[i]->uniformRefine(8);
    }
    gsMultiPatch<real_t> target(patches,topolGeom.boundaries(),topolGeom.interfaces());
    gsQualityMeasureWeights weights(target.nPatches());
    weights.m_uniformity=0.0001;
    weights.m_approxPoints=1;
    gsFileData<real_t> dataInsidePatch( "Antonella/QuestionMark/InnerPoints.xml" );
    gsFileData<real_t> dataBoundaryHoles( "Antonella/QuestionMark/OnlyHOLESpoints.xml" );
    gsFileData<real_t> dataBoundary( "Antonella/QuestionMark/Boundary.xml" );
    gsMatrix<real_t>::uPtr parsInside,  parsBoundary,  parsBoundaryHoles;
    gsMatrix<real_t>::uPtr pointsInside,pointsBoundary,pointsBoundaryHoles;
    for(size_t i = 0;i<target.nPatches();++i)
    {
        const int id_par=2*i;
        const int id_point=2*i+1;
        parsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_par);
        pointsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_point);
        parsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_par);
        pointsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_point);
        parsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_par);
        pointsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_point);
        if(parsInside && parsBoundaryHoles)
        {
            const index_t totalCols = parsInside->cols()+parsBoundaryHoles->cols();
            const index_t totalRows = parsInside->rows();
            weights.m_pars[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_points[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_pars[i].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
            weights.m_pars[i].rightCols(parsBoundaryHoles->cols())=parsBoundaryHoles->rightCols(parsBoundaryHoles->cols());
            weights.m_points[i].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
            weights.m_points[i].rightCols(pointsBoundaryHoles->cols())=pointsBoundaryHoles->rightCols(pointsBoundaryHoles->cols());
        }
        else if(parsInside && parsBoundary)
        {
            const index_t totalCols = parsInside->cols()+parsBoundary->cols();
            const index_t totalRows = parsInside->rows();
            weights.m_pars[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_points[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_pars[i].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
            weights.m_pars[i].rightCols(parsBoundary->cols())=parsBoundary->rightCols(parsBoundary->cols());
            weights.m_points[i].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
            weights.m_points[i].rightCols(pointsBoundary->cols())=pointsBoundary->rightCols(pointsBoundary->cols());
        }
        else if(parsInside)
        {
            weights.m_pars[i]=*parsInside;
            weights.m_points[i]=*pointsInside;
        }
        else if(parsBoundary)
        {
            weights.m_pars[i]=*parsBoundary;
            weights.m_points[i]=*pointsBoundary;
        }
        else if(parsBoundaryHoles)
        {
            weights.m_pars[i]=*parsBoundaryHoles;
            weights.m_points[i]=*pointsBoundaryHoles;
        }
        gsInfo << "\n";
    }
    // real_t time;
    // unsigned dof;
    gsCompositeIncrSmoothnessGeom<2,real_t> compGeom(target,0);

    gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,false);
    adaptiveSolver.adaptiveSolve();

    printMPToParaview(compGeom.exportToPatches(),"firstOptSecond"+util::to_string(0),10000);

    for(size_t i =0;i<target.nPatches();++i)
    {
        if(weights.m_pars[i].rows()!=0 && weights.m_pars[i].cols()!=0)
            gsWriteParaviewPoints(weights.m_points[i],"pointsOfPatch"+util::to_string(i));
    }

    return true;
}

bool test_measureSolvingPlatePatch(unsigned patch)
{
    std::cout << "solve for patch: " << patch << std::endl;
    std::vector<gsGeometry<real_t> *> patches;
    const size_t degree=2;
    gsKnotVector<real_t>kv(0,1,0,degree+1,1,degree);
    gsMatrix<real_t> coefs;
    const size_t coefsInDir=kv.size()-degree-1;
    coefs.setZero(coefsInDir*coefsInDir,2);
    gsTensorBSpline<2,real_t> * patchGeom = new gsTensorBSpline<2,real_t>(kv,kv,coefs);
    patchGeom->uniformRefine(1);
    patches.push_back(new gsTHBSpline<2,real_t>(*patchGeom));

    gsMultiPatch<real_t> target(patches);
    target.addAutoBoundaries();
    gsQualityMeasureWeights weights(1);
    //weights.m_length=0.000001;
    weights.m_uniformity=0.00001;
    weights.m_approxPoints=1;
    gsFileData<real_t> dataInsidePatch( "Antonella/Plate/PatchSampled.xml" );
    gsFileData<real_t> dataBoundaryHoles( "Antonella/Plate/HolesB_M_Plate_patch.xml" );
    gsFileData<real_t> dataBoundary( "Antonella/Plate/Boundary.xml" );
    gsMatrix<real_t>::uPtr parsInside,  parsBoundary,  parsBoundaryHoles;
    gsMatrix<real_t>::uPtr pointsInside,pointsBoundary,pointsBoundaryHoles;
    const int id_par=2*patch;
    const int id_point=2*patch+1;
    parsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_par);
    pointsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_point);
    parsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_par);
    pointsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_point);
    parsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_par);
    pointsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_point);
    if(parsInside && parsBoundaryHoles)
    {
        const index_t totalCols = parsInside->cols()+parsBoundaryHoles->cols();
        const index_t totalRows = parsInside->rows();
        weights.m_pars[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_points[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_pars[0].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
        weights.m_pars[0].rightCols(parsBoundaryHoles->cols())=parsBoundaryHoles->rightCols(parsBoundaryHoles->cols());
        weights.m_points[0].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
        weights.m_points[0].rightCols(pointsBoundaryHoles->cols())=pointsBoundaryHoles->rightCols(pointsBoundaryHoles->cols());
    }
    else if(parsInside && parsBoundary)
    {
        const index_t totalCols = parsInside->cols()+parsBoundary->cols();
        const index_t totalRows = parsInside->rows();
        weights.m_pars[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_points[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_pars[0].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
        weights.m_pars[0].rightCols(parsBoundary->cols())=parsBoundary->rightCols(parsBoundary->cols());
        weights.m_points[0].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
        weights.m_points[0].rightCols(pointsBoundary->cols())=pointsBoundary->rightCols(pointsBoundary->cols());
    }
    else if(parsInside)
    {
        weights.m_pars[0]=*parsInside;
        weights.m_points[0]=*pointsInside;
    }
    else if(parsBoundary)
    {
        weights.m_pars[0]=*parsBoundary;
        weights.m_points[0]=*pointsBoundary;
    }
    else if(parsBoundaryHoles)
    {
        weights.m_pars[0]=*parsBoundaryHoles;
        weights.m_points[0]=*pointsBoundaryHoles;
    }
    std::cout << std::endl;

    // real_t time;
    // unsigned dof;
    gsCompositeIncrSmoothnessGeom<2,real_t> compGeom(target,-1);

    gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,false);
    adaptiveSolver.adaptiveSolve();

    gsWriteParaviewPoints(weights.m_points[0],"pointsOfPatch"+util::to_string(patch));
    gsMultiPatch<real_t> mp2 = compGeom.exportToPatches();
    printMPToParaview(mp2,"firstOptSecond"+util::to_string(patch),10000);
    return true;
}

bool test_measureSolvingPlate()
{
    std::string path="Antonella/Plate/2Template_M_Plate_patches.xml";
    gsMultiPatch<> topolGeom;
    try
    {
        gsReadFile<>(path, topolGeom);
    }
    catch (std::runtime_error&)
    {
        return false;
    }

    if(topolGeom.empty())
    {
        gsInfo << "No gsMultiPatch in file." << "\n" ;
        return false;
    }

    std::vector<gsGeometry<real_t> *> patches;
    const size_t degree=2;
    for(size_t i = 0;i<topolGeom.nPatches();++i)
    {
        gsKnotVector<real_t>kv(0,1,0,degree+1,1,degree);
        gsMatrix<real_t> coefs;
        const size_t coefsInDir=kv.size()-degree-1;
        coefs.setZero(coefsInDir*coefsInDir,2);
        gsTensorBSpline<2,real_t> * patchGeom = new gsTensorBSpline<2,real_t>(kv,kv,coefs);
        patchGeom->uniformRefine(1);
        patches.push_back(new gsTHBSpline<2,real_t>(*patchGeom));
    }
    gsMultiPatch<real_t> target(patches,topolGeom.boundaries(),topolGeom.interfaces());
    gsQualityMeasureWeights weights(target.nPatches());
    weights.m_uniformity=0.0001;
    weights.m_approxPoints=1;
    gsFileData<real_t> dataInsidePatch( "Antonella/Plate/PatchSampled.xml" );
    gsFileData<real_t> dataBoundaryHoles( "Antonella/Plate/Boundary.xml" );
    //gsFileData<real_t> dataBoundaryHoles( "Antonella/Plate/HolesB_M_Plate_patch.xml" );
    gsFileData<real_t> dataBoundary( "Antonella/Plate/Boundary.xml" );
    gsMatrix<real_t>::uPtr parsInside,  parsBoundary,  parsBoundaryHoles;
    gsMatrix<real_t>::uPtr pointsInside,pointsBoundary,pointsBoundaryHoles;
    for(size_t i = 0;i<target.nPatches();++i)
    {
        const int id_par=2*i;
        const int id_point=2*i+1;
        parsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_par);
        pointsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_point);
        parsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_par);
        pointsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_point);
        parsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_par);
        pointsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_point);
        if(parsInside && parsBoundaryHoles)
        {
            const index_t totalCols = parsInside->cols()+parsBoundaryHoles->cols();
            const index_t totalRows = parsInside->rows();
            weights.m_pars[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_points[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_pars[i].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
            weights.m_pars[i].rightCols(parsBoundaryHoles->cols())=parsBoundaryHoles->rightCols(parsBoundaryHoles->cols());
            weights.m_points[i].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
            weights.m_points[i].rightCols(pointsBoundaryHoles->cols())=pointsBoundaryHoles->rightCols(pointsBoundaryHoles->cols());
        }
        else if(parsInside && parsBoundary)
        {
            const index_t totalCols = parsInside->cols()+parsBoundary->cols();
            const index_t totalRows = parsInside->rows();
            weights.m_pars[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_points[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_pars[i].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
            weights.m_pars[i].rightCols(parsBoundary->cols())=parsBoundary->rightCols(parsBoundary->cols());
            weights.m_points[i].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
            weights.m_points[i].rightCols(pointsBoundary->cols())=pointsBoundary->rightCols(pointsBoundary->cols());
        }
        else if(parsInside)
        {
            weights.m_pars[i]=*parsInside;
            weights.m_points[i]=*pointsInside;
        }
        else if(parsBoundary)
        {
            weights.m_pars[i]=*parsBoundary;
            weights.m_points[i]=*pointsBoundary;
        }
        else if(parsBoundaryHoles)
        {
            weights.m_pars[i]=*parsBoundaryHoles;
            weights.m_points[i]=*pointsBoundaryHoles;
        }
        std::cout << std::endl;
    }
    // real_t time;
    // unsigned dof;
    gsCompositeIncrSmoothnessGeom<2,real_t> compGeom(target,0);

    gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,false);
    adaptiveSolver.adaptiveSolve();

    printMPToParaview(compGeom.exportToPatches(),"firstOptSecond"+util::to_string(0),10000);

    for(size_t i =0;i<target.nPatches();++i)
    {
        if(weights.m_pars[i].rows()!=0 && weights.m_pars[i].cols()!=0)
            gsWriteParaviewPoints(weights.m_points[i],"pointsOfPatch"+util::to_string(i));
    }

    return true;
}

bool test_measureSolvingMountingPlatePatch(unsigned patch)
{
    std::cout << "solve for patch: " << patch << std::endl;
    std::vector<gsGeometry<real_t> *> patches;
    const size_t degree=2;
    gsKnotVector<real_t>kv(0,1,0,degree+1,1,degree);
    gsMatrix<real_t> coefs;
    const size_t coefsInDir=kv.size()-degree-1;
    coefs.setZero(coefsInDir*coefsInDir,2);
    gsTensorBSpline<2,real_t> * patchGeom = new gsTensorBSpline<2,real_t>(kv,kv,coefs);
    patchGeom->uniformRefine(1);
    patches.push_back(new gsTHBSpline<2,real_t>(*patchGeom));

    gsMultiPatch<real_t> target(patches);
    target.addAutoBoundaries();
    gsQualityMeasureWeights weights(1);
    //weights.m_length=0.000001;
    weights.m_uniformity=0.00001;
    weights.m_approxPoints=1;
    gsFileData<real_t> dataInsidePatch( "Antonella/MountingPlate/PatchSampled.xml" );
    gsFileData<real_t> dataBoundaryHoles( "Antonella/MountingPlate/BoundaryHoles.xml" );
    gsFileData<real_t> dataBoundary( "Antonella/MountingPlate/Boundary.xml" );
    gsMatrix<real_t>::uPtr parsInside,  parsBoundary,  parsBoundaryHoles;
    gsMatrix<real_t>::uPtr pointsInside,pointsBoundary,pointsBoundaryHoles;
    const int id_par=2*patch;
    const int id_point=2*patch+1;
    parsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_par);
    pointsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_point);
    parsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_par);
    pointsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_point);
    parsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_par);
    pointsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_point);
    if(parsInside && parsBoundaryHoles)
    {
        const index_t totalCols = parsInside->cols()+parsBoundaryHoles->cols();
        const index_t totalRows = parsInside->rows();
        weights.m_pars[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_points[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_pars[0].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
        weights.m_pars[0].rightCols(parsBoundaryHoles->cols())=parsBoundaryHoles->rightCols(parsBoundaryHoles->cols());
        weights.m_points[0].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
        weights.m_points[0].rightCols(pointsBoundaryHoles->cols())=pointsBoundaryHoles->rightCols(pointsBoundaryHoles->cols());
    }
    else if(parsInside && parsBoundary)
    {
        const index_t totalCols = parsInside->cols()+parsBoundary->cols();
        const index_t totalRows = parsInside->rows();
        weights.m_pars[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_points[0]=gsMatrix<real_t>(totalRows,totalCols);
        weights.m_pars[0].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
        weights.m_pars[0].rightCols(parsBoundary->cols())=parsBoundary->rightCols(parsBoundary->cols());
        weights.m_points[0].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
        weights.m_points[0].rightCols(pointsBoundary->cols())=pointsBoundary->rightCols(pointsBoundary->cols());
    }
    else if(parsInside)
    {
        weights.m_pars[0]=*parsInside;
        weights.m_points[0]=*pointsInside;
    }
    else if(parsBoundary)
    {
        weights.m_pars[0]=*parsBoundary;
        weights.m_points[0]=*pointsBoundary;
    }
    else if(parsBoundaryHoles)
    {
        weights.m_pars[0]=*parsBoundaryHoles;
        weights.m_points[0]=*pointsBoundaryHoles;
    }
    std::cout << std::endl;

    // real_t time;
    // unsigned dof;
    gsCompositeIncrSmoothnessGeom<2,real_t> compGeom(target,-1);

    gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,false);
    adaptiveSolver.adaptiveSolve();

    gsWriteParaviewPoints(weights.m_points[0],"pointsOfPatch"+util::to_string(patch));
    gsMultiPatch<real_t> mp2 = compGeom.exportToPatches();
    printMPToParaview(mp2,"firstOptSecond"+util::to_string(patch),10000);
    return true;
}

bool test_measureSolvingMountingPlate()
{
    std::string path="Antonella/MountingPlate/2Template_M_Plate_patches.xml";
    gsMultiPatch<> topolGeom;
    try
    {
        gsReadFile<>(path, topolGeom);
    }
    catch (std::runtime_error&)
    {
        return false;
    }

    if(topolGeom.empty())
    {
        gsInfo << "No gsMultiPatch in file." << "\n" ;
        return false;
    }
    std::vector<gsGeometry<real_t> *> patches;
    const size_t degree=2;
    for(size_t i = 0;i<topolGeom.nPatches();++i)
    {
        gsKnotVector<real_t>kv(0,1,0,degree+1,1,degree);
        gsMatrix<real_t> coefs;
        const size_t coefsInDir=kv.size()-degree-1;
        coefs.setZero(coefsInDir*coefsInDir,2);
        gsTensorBSpline<2,real_t> * patchGeom = new gsTensorBSpline<2,real_t>(kv,kv,coefs);
        patchGeom->uniformRefine(1);
        patches.push_back(new gsTHBSpline<2,real_t>(*patchGeom));
    }
    gsMultiPatch<real_t> target(patches,topolGeom.boundaries(),topolGeom.interfaces());
    gsQualityMeasureWeights weights(target.nPatches());
    weights.m_uniformity=0.001;
    weights.m_approxPoints=1;
    gsFileData<real_t> dataInsidePatch( "Antonella/MountingPlate/PatchSampled.xml" );
    gsFileData<real_t> dataBoundaryHoles( "Antonella/MountingPlate/BoundaryHoles.xml" );
    gsFileData<real_t> dataBoundary( "Antonella/MountingPlate/Boundary.xml" );
    gsMatrix<real_t>::uPtr parsInside,  parsBoundary,  parsBoundaryHoles;
    gsMatrix<real_t>::uPtr pointsInside,pointsBoundary,pointsBoundaryHoles;
    for(size_t i = 0;i<target.nPatches();++i)
    {
        const int id_par=2*i;
        const int id_point=2*i+1;
        parsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_par);
        pointsInside = dataInsidePatch.getId<gsMatrix<real_t> >(id_point);
        parsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_par);
        pointsBoundary = dataBoundary.getId<gsMatrix<real_t> >(id_point);
        parsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_par);
        pointsBoundaryHoles = dataBoundaryHoles.getId<gsMatrix<real_t> >(id_point);
        if(parsInside && parsBoundaryHoles)
        {
            const index_t totalCols = parsInside->cols()+parsBoundaryHoles->cols();
            const index_t totalRows = parsInside->rows();
            weights.m_pars[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_points[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_pars[i].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
            weights.m_pars[i].rightCols(parsBoundaryHoles->cols())=parsBoundaryHoles->rightCols(parsBoundaryHoles->cols());
            weights.m_points[i].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
            weights.m_points[i].rightCols(pointsBoundaryHoles->cols())=pointsBoundaryHoles->rightCols(pointsBoundaryHoles->cols());
        }
        else if(parsInside && parsBoundary)
        {
            const index_t totalCols = parsInside->cols()+parsBoundary->cols();
            const index_t totalRows = parsInside->rows();
            weights.m_pars[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_points[i]=gsMatrix<real_t>(totalRows,totalCols);
            weights.m_pars[i].leftCols(parsInside->cols())=parsInside->leftCols(parsInside->cols());
            weights.m_pars[i].rightCols(parsBoundary->cols())=parsBoundary->rightCols(parsBoundary->cols());
            weights.m_points[i].leftCols(pointsInside->cols())=pointsInside->leftCols(pointsInside->cols());
            weights.m_points[i].rightCols(pointsBoundary->cols())=pointsBoundary->rightCols(pointsBoundary->cols());
        }
        else if(parsInside)
        {
            weights.m_pars[i]=*parsInside;
            weights.m_points[i]=*pointsInside;
        }
        else if(parsBoundary)
        {
            weights.m_pars[i]=*parsBoundary;
            weights.m_points[i]=*pointsBoundary;
        }
        else if(parsBoundaryHoles)
        {
            weights.m_pars[i]=*parsBoundaryHoles;
            weights.m_points[i]=*pointsBoundaryHoles;
        }
        std::cout << std::endl;
    }
    // real_t time;
    // unsigned dof;
    gsCompositeIncrSmoothnessGeom<2,real_t> compGeom(target,0);

    gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,false);
    adaptiveSolver.adaptiveSolve();

    printMPToParaview(compGeom.exportToPatches(),"firstOptSecond"+util::to_string(0),10000);

    for(size_t i =0;i<target.nPatches();++i)
    {
        if(weights.m_pars[i].rows()!=0 && weights.m_pars[i].cols()!=0)
            gsWriteParaviewPoints(weights.m_points[i],"pointsOfPatch"+util::to_string(i));
    }

    return true;
}



bool test_measureSolvingExample1()
{
    //std::string path2="planar/two_squares_differentSize.xml";
    //std::string path2="Mario/square_deg4.xml";
    std::string path2="TopologyOptimization/MultiPatch_0.xml";
    gsMultiPatch<> distortedGeom;   // TODO: remove? never used!
    try
    {
        gsReadFile<>(path2, distortedGeom);
    }
    catch (std::runtime_error&)
    {
        return false;
    }

    if(distortedGeom.empty())
    {
        gsInfo << "No gsMultiPatch in file." << "\n" ;
        return false;
    }
    // set up problem
    gsInfo<< "\n------------ measureSolving ------------\n\n";
    std::string path=path2;
    //std::string path="Mario/square_deg4.xml";
    //std::string path="planar/four_squares2.xml";
    //std::string path="surfaces/multipatch_triangle2d.xml";
    //std::string path="planar/multipatch_tunnel_thb.xml";
    //std::string path="surfaces/multipatch_AirPassage.xml";
    //std::string path="Antonella/mp1010.xml";
    //std::string path="planar/two_squares_differentSize_distorted.xml";
    //std::string path="planar/two_squares_differentSize.xml";
    //std::string path="planar/puzzleEasy.xml";
    gsMultiPatch<> mp;
    try
    {
        gsReadFile<>(path2, mp);
    }
    catch (std::runtime_error&)
    {
        return false;
    }

    if(mp.empty())
    {
        gsInfo << "No gsMultiPatch in file." << "\n" ;
        return false;
    }
//    mp->patch(0).degreeElevate(1,1);
//    mp->patch(1).degreeElevate(1,-1);
//    for(size_t j=0;j<mp->nPatches();j++)
//        mp->patch(j).uniformRefine();
    for(size_t j=0;j<mp.nPatches();j++)
        mp.patch(j).embed(2);
//    mp->patch(0).uniformRefine(2);
//    mp->patch(1).uniformRefine(2);
    gsCompositeIncrSmoothnessGeom<2,real_t> compGeom2(mp,0);
//    std::vector<real_t> knots_u;
//    std::vector<real_t> knots_v;
//    for(unsigned i = 1;i<9;++i)
//        knots_v.push_back((static_cast<real_t>(i))/9);
//    static_cast<gsCompositeBSplineBasis<2,real_t>&>(compGeom2.getCompBasis()).refine_withCoefs(compGeom2.coefs(),0,knots_u,knots_v,true);
    gsMultiPatch<> mp3 = compGeom2.exportToPatches();
//    mp3->degreeElevate(1);
//    gsFileData<real_t> fd;
//    fd<<*mp3;
//    fd.dump("puzzleEasy.xml");
    gsCompositeIncrSmoothnessGeom<2,real_t> compGeom(mp3,0);

    // real_t time;
    // unsigned dof;

    gsMultiPatch<real_t> mp2 = compGeom.exportToPatches();
    gsWriteParaview(mp2,"firstOptStart");
    gsQualityMeasureWeights weights(mp.nPatches());
    weights.m_length=1;
    weights.m_orthogonality=0; //0.1;
    weights.m_uniformity=0;
    weights.m_skewness=0;
    weights.m_eccentricity=0;
    weights.m_area=0; //1;
    weights.m_selfIntersect=0;
    weights.m_areaInverse=0;
    weights.m_epsilon=1e-10;
    weights.m_weightArea=false;
    weights.m_approxPoints=0;
    gsVector<real_t,2> start,end;
    start << 0.1,0.1;
    end << 0.9,0.9;
    gsVector<unsigned,2> np;
    np << 10,10;
//    weights.m_pars[0] = gsPointGrid<real_t>(start,end,np);
//    weights.m_pars[1] = gsPointGrid<real_t>(start,end,np);
//    weights.m_points[0] = distortedGeom->patch(0).eval(weights.m_pars[0]);
//    weights.m_points[1] = distortedGeom->patch(1).eval(weights.m_pars[1]);
//    for(size_t i = 0;i<compGeom.getCompBasis().nPatches();++i)
//    {
//        gsVector<real_t,2> start,end;
//        start << 0.2,0.2;
//        end << 0.9,0.9;
//        gsVector<unsigned,2> np;
//        np << 10,10;
//        if(i!=6)
//            weights.m_pars[i] = gsPointGrid<real_t>(start,end,np);
//        else
//        {
//            gsVector<real_t,2> start2,end2;
//            start2 << 0.1,0.1;
//            end2 << 0.8,0.8;
//            weights.m_pars[i] = gsPointGrid<real_t>(start2,end2,np);
//        }
//        gsMatrix<real_t> * input;
//        std::stringstream ss;
//        ss<<"Antonella/PatchPoints/Patch"<<i<<"_Points.xml";
//        gsInfo << ss.str() << "\n";
//        input = gsReadFile<>(ss.str());
//        GISMO_ASSERT(input,"matrix not found in file");
//        weights.m_points[i]=*input;
//        delete input;
//    }

//    gsFitting<real_t> fitting(weights.m_pars[0],weights.m_points[0], mp->basis(0));
//    fitting.compute();
//    gsWriteParaview(*fitting.result(),"resultFittingClass",1000);

    gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,true);
    adaptiveSolver.adaptiveSolve();

    gsMultiPatch<real_t> mp4 = compGeom.exportToPatches();
    gsWriteParaview(mp4,"optResult");
    // Print out the L2 errors and element sizes
    /*
    gsInfo << std::scientific;
    gsInfo << "\ndofs: "<<dof;
    gsInfo << "\ntimes: "<<time;
    gsInfo << "\n";
    gsInfo << std::fixed;
    */

    //gsWriteParaviewPoints(weights.m_points[0],"patch0Points");
    //gsWriteParaview(compGeom.exportToPatches()->patch(0),"patch0",1000);
    //gsWriteParaviewPoints(weights.m_points[1],"patch2Points");
    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

gsMultiPatch<> splitGeometry(const gsGeometry<> *geom1,real_t split_parV)
{
    gsTensorBSpline<2,real_t>::uPtr cloneTBSpline = memory::convert_ptr<gsTensorBSpline<2,real_t> >(geom1->clone());
    if(!cloneTBSpline)
    {
        gsInfo<<"no gsTensorBSpline"<<"\n";
        GISMO_ERROR("no gsTensorBSpline");
    }
    std::vector<std::vector<real_t> > refineKnots;
    std::vector<real_t> knots;
    refineKnots.push_back(knots);// no knots inserted in u-dir
    unsigned mult = cloneTBSpline->basis().knots(1).multiplicity(split_parV);
    for(unsigned i = 0;i<cloneTBSpline->basis().degree(1)-mult;++i)
        knots.push_back(split_parV);
    refineKnots.push_back(knots);// knots inserted in v-dir
    cloneTBSpline->basis().refine_withCoefs(cloneTBSpline->coefs(),refineKnots);

    std::vector<real_t> kvec1;
    std::vector<real_t> kvec2;
    kvec2.push_back(split_parV);
    for(size_t i = 0;i<cloneTBSpline->basis().knots(1).size();++i)
    {
        const real_t curKnot = cloneTBSpline->basis().knots(1)[i];
        if(curKnot<split_parV)
            kvec1.push_back(curKnot);
        else if(curKnot>split_parV)
            kvec2.push_back(curKnot);
        else
        {
            kvec1.push_back(curKnot);
            kvec2.push_back(curKnot);
        }
    }
    kvec1.push_back(kvec1.back());

    gsKnotVector<> knots1(give(kvec1), cloneTBSpline->basis().degree(0));
    gsKnotVector<> knots2(give(kvec2), cloneTBSpline->basis().degree(0));
    knots1.transform(0,1);
    knots2.transform(0,1);
    gsTensorBSplineBasis<2,real_t> tbasis1(cloneTBSpline->basis().knots(0),knots1);
    gsTensorBSplineBasis<2,real_t> tbasis2(cloneTBSpline->basis().knots(0),knots2);

    gsMatrix<real_t> coefs=cloneTBSpline->coefs();
    gsMatrix<real_t> coefs1=coefs.bottomRows(tbasis1.size());
    gsMatrix<real_t> coefs2=coefs.topRows(tbasis2.size());

    std::vector<gsGeometry<>* >patches;
    patches.push_back(new gsTensorBSpline<2,real_t>(tbasis1,coefs1));
    patches.push_back(new gsTensorBSpline<2,real_t>(tbasis2,coefs2));
    gsMultiPatch<> mp(patches);
    mp.computeTopology();
    return mp;
}

bool test_SplitInvariance()
{
    gsInfo<< "\n------------ measureSolving :: test_SplitInvariance ------------\n\n";

    //std::string pathMP="planar/two_squares_differentSize.xml";
    std::string pathSP="planar/one_square_bent.xml";

    gsMultiPatch<> mpSP;
    try
    {
        gsReadFile<>(pathSP, mpSP);
    }
    catch (std::runtime_error&)
    {
        return false;
    }

    if(mpSP.empty())
    {
        gsInfo << "No gsMultiPatch in file." << "\n" ;
        return false;
    }

    gsMultiPatch<> mpMP = splitGeometry(&mpSP.patch(0),0.7);

    for(size_t j=0;j<mpSP.nPatches();j++)
    {
        mpSP.patch(j).embed(2);
        mpSP.patch(j).degreeElevate();
        mpSP.patch(j).uniformRefine(1);
    }
    for(size_t j=0;j<mpMP.nPatches();j++)
    {
        mpMP.patch(j).embed(2);
        mpMP.patch(j).degreeElevate();
        mpMP.patch(j).uniformRefine(1);
    }

    gsInfo << std::flush;

    for(unsigned i = 1;i<=8;++i)
    {
        gsCompositeIncrSmoothnessGeom<2,real_t> geomSP(mpSP,0);
        gsCompositeIncrSmoothnessGeom<2,real_t> geomMP(mpMP,0);

        gsQualityMeasureWeights weightsMP(geomMP.getCompBasis().nPatches());
        weightsMP.m_length=0;
        weightsMP.m_orthogonality=0;
        weightsMP.m_uniformity=0;
        weightsMP.m_skewness=0;
        weightsMP.m_eccentricity=0;
        weightsMP.m_area=0;
        weightsMP.m_selfIntersect=0;
        weightsMP.m_areaInverse=0;
        weightsMP.m_epsilon=1e-10;
        weightsMP.m_weightArea=true;

        switch(i)
        {
        case 1:
            gsInfo<<"length: ";
            weightsMP.m_length=1;
            break;
        case 2:
            gsInfo<<"ortho: ";
            weightsMP.m_orthogonality=1;
            break;
        case 3:
            gsInfo<<"uniformity: ";
            weightsMP.m_uniformity=1;
            break;
        case 4:
            gsInfo<<"skewness: ";
            weightsMP.m_skewness=1;
            break;
        case 5:
            gsInfo<<"eccentricity: ";
            weightsMP.m_eccentricity=1;
            break;
        case 6:
            gsInfo<<"area: ";
            weightsMP.m_area=1;
            break;
        case 7:
            gsInfo<<"selfIntersect: ";
            weightsMP.m_selfIntersect=1;
            break;
        case 8:
            gsInfo<<"areaInverse: ";
            weightsMP.m_areaInverse=1;
            break;
        default:;
        }

        gsQualityMeasureWeights weightsSP(1);
        weightsSP.m_length        = weightsMP.m_length;
        weightsSP.m_orthogonality = weightsMP.m_orthogonality;
        weightsSP.m_uniformity    = weightsMP.m_uniformity;
        weightsSP.m_skewness      = weightsMP.m_skewness;
        weightsSP.m_eccentricity  = weightsMP.m_eccentricity;
        weightsSP.m_area          = weightsMP.m_area;
        weightsSP.m_selfIntersect = weightsMP.m_selfIntersect;
        weightsSP.m_areaInverse   = weightsMP.m_areaInverse;
        weightsSP.m_epsilon       = weightsMP.m_epsilon;
        weightsSP.m_weightArea    = weightsMP.m_weightArea;

        gsMatrix<real_t> *qualityMeasureVals1=new gsMatrix<real_t>();
        gsMultiPatch<real_t> domain1 = geomSP.exportToPatches();
        gsMappedBasis<2,real_t> * compBasis1 = &geomSP.getCompBasis();
        gsMatrix<real_t> *qualityMeasureVals2=new gsMatrix<real_t>();
        gsMultiPatch<real_t> domain2 = geomMP.exportToPatches();
        gsMappedBasis<2,real_t> * compBasis2 = &geomMP.getCompBasis();

        if(weightsMP.m_weightArea)
        {
            calculateAreas(compBasis1,domain1,*qualityMeasureVals1);
            weightsSP.m_areas=*qualityMeasureVals1;
            //gsInfo << "\n" << "areas: " << qualityMeasureVals1->transpose() << "\n";

            calculateAreas(compBasis2,domain2,*qualityMeasureVals2);
            weightsMP.m_areas=*qualityMeasureVals2;
            //gsInfo << "\n" << "areas: " << qualityMeasureVals2->transpose() << "\n";
        }

        calculateMeasures(compBasis1,domain1,*qualityMeasureVals1,weightsSP);
        calculateMeasures(compBasis2,domain2,*qualityMeasureVals2,weightsMP);
        //    gsInfo<<"qualityMeasureValsSP: "<<*qualityMeasureVals1<<"\n";
        //    gsInfo<<"qualityMeasureValsMP: "<<*qualityMeasureVals2<<"\n";

        gsInfo<<"Difference: "<<(*qualityMeasureVals2-*qualityMeasureVals1)(0,0)<<"\n";
    }

    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

bool test_ScalingInvariance()
{
    gsInfo<< "\n------------ measureSolving :: test_ScalingInvariance ------------\n\n";

    //std::string pathMP="planar/two_squares_differentSize.xml";
    std::string pathSP="planar/multipatch_tunnel.xml";
    //std::string pathSP="planar/one_square_bent.xml";
    real_t scaling=10;

    gsMultiPatch<> mpMP_single;
    try
    {
        gsReadFile<>(pathSP, mpMP_single);
    }
    catch (std::runtime_error&)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsMultiPatch<> mpMP = mpMP_single;//splitGeometry(mpMP_single.patch(0),0.2);

    for(size_t j=0;j<mpMP.nPatches();j++)
    {
        mpMP.patch(j).embed(2);
        mpMP.patch(j).degreeElevate();
        mpMP.patch(j).uniformRefine(1);
    }

    gsMultiPatch<> mpSP = gsMultiPatch<>(mpMP);
    for(size_t j=0;j<mpSP.nPatches();j++)
    {
        mpSP.patch(j).coefs()*=scaling;
    }

    gsInfo << std::flush;

    for(unsigned i = 9;i<=9;++i)
    {
        gsCompositeIncrSmoothnessGeom<2,real_t> geomSP(mpMP,0);
        gsCompositeIncrSmoothnessGeom<2,real_t> geomMP(mpSP,0);

        gsQualityMeasureWeights weightsMP(geomMP.getCompBasis().nPatches());
        weightsMP.m_epsilon=1e-10;
        weightsMP.m_weightArea=true;

        switch(i)
        {
        case 1:
            gsInfo<<"length: \n";
            weightsMP.m_length=1;
            break;
        case 2:
            gsInfo<<"ortho: \n";
            weightsMP.m_orthogonality=1;
            break;
        case 3:
            gsInfo<<"uniformity: \n";
            weightsMP.m_uniformity=1;
            break;
        case 4:
            gsInfo<<"skewness: \n";
            weightsMP.m_skewness=1;
            break;
        case 5:
            gsInfo<<"eccentricity: \n";
            weightsMP.m_eccentricity=1;
            break;
        case 6:
            gsInfo<<"area: \n";
            weightsMP.m_area=1;
            break;
        case 7:
            gsInfo<<"selfIntersect: \n";
            weightsMP.m_selfIntersect=1;
            break;
        case 8:
            gsInfo<<"areaInverse: \n";
            weightsMP.m_areaInverse=1;
            break;
        case 9:
            std::cout<<"mixed functionals: \n";
            weightsMP.m_length=1;
//            weightsMP.m_uniformity=1;
            weightsMP.m_orthogonality=1;
//            weightsMP.m_area=1;
            weightsMP.m_eccentricity=1;
            break;
        default:;
        }

        gsQualityMeasureWeights weightsSP(geomMP.getCompBasis().nPatches());
        weightsSP.m_length        = weightsMP.m_length;
        weightsSP.m_orthogonality = weightsMP.m_orthogonality;
        weightsSP.m_uniformity    = weightsMP.m_uniformity;
        weightsSP.m_skewness      = weightsMP.m_skewness;
        weightsSP.m_eccentricity  = weightsMP.m_eccentricity;
        weightsSP.m_area          = weightsMP.m_area;
        weightsSP.m_selfIntersect = weightsMP.m_selfIntersect;
        weightsSP.m_areaInverse   = weightsMP.m_areaInverse;
        weightsSP.m_epsilon       = weightsMP.m_epsilon;
        weightsSP.m_weightArea    = weightsMP.m_weightArea;

        gsMatrix<real_t> *qualityMeasureVals1=new gsMatrix<real_t>();
        gsMultiPatch<real_t> domain1 = geomSP.exportToPatches();
        gsMappedBasis<2,real_t> * compBasis1 = &geomSP.getCompBasis();
        gsMatrix<real_t> *qualityMeasureVals2=new gsMatrix<real_t>();
        gsMultiPatch<real_t> domain2 = geomMP.exportToPatches();
        gsMappedBasis<2,real_t> * compBasis2 = &geomMP.getCompBasis();

        if(weightsMP.m_weightArea)
        {
            calculateAreas(compBasis1,domain1,*qualityMeasureVals1);
            weightsSP.m_areas=*qualityMeasureVals1;
            gsInfo << "areas: " << qualityMeasureVals1->transpose() << "\n";

            calculateAreas(compBasis2,domain2,*qualityMeasureVals2);
            weightsMP.m_areas=*qualityMeasureVals2;
            gsInfo << "areas: " << qualityMeasureVals2->transpose() << "\n";
        }

        // unsigned dof;
        // real_t time;

        gsAdaptiveSolverOptimizeGeometry adaptiveSolver1(&geomMP,weightsMP,true);
        adaptiveSolver1.adaptiveSolve();
        gsAdaptiveSolverOptimizeGeometry adaptiveSolver2(&geomSP,weightsSP,true);
        adaptiveSolver2.adaptiveSolve();

        real_t diff=0;
        for(size_t j=0;j<mpSP.nPatches();j++)
        {
            diff+=(mpSP.patch(j).coefs()*(1/scaling)-mpMP.patch(j).coefs()).norm();
        }
        gsInfo<<"Difference: "<<diff<<"\n";
    }

    gsInfo << "test passed\n\n" ;
    gsInfo << flush;
    return true;
}

void printTemplate()
{
    std::string path="Antonella/QuestionMark/2question_mark_template.xml";
    gsMultiPatch<> * topolGeom = (  (gsMultiPatch<>::uPtr)gsReadFile<>(path)  ).release();
    if(topolGeom==NULL)
    {
        std::cout << "File could not be read." << "\n" ;
        return;
    }

    printMPToParaview(*topolGeom,"template",1000);
}

int main(int argc, char *argv[])
{
//    bool passed = true;

    bool plot       = false;
    bool timeTest   = false;
    index_t nsamples    = 1000;

    // Read input
    gsCmdLine cmd("Composite basis tests.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot    );
    cmd.addSwitch("time", "Test evaluation time"          , timeTest);
    cmd.addInt("s","samples", "Number of samples to use for viewing", nsamples);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsWarn << "Test DISABLED.\n";

//    int plotNr = 0;
//    real_t eps = math::sqrt(math::numeric_limits<real_t>::epsilon())/100;

//    int i = 0;
    //test_measureSolvingExample1();
    //test_SplitInvariance();
    //test_ScalingInvariance();
    //test_measureSolvingHat();
    //for(unsigned i = 0;i<=15;++i)
    //    test_measureSolvingHatPatch(i);
    //for(unsigned i = 1;i<2;++i)
      //  test_measureSolvingPlatePatch(i);
    //test_measureSolvingPlate();
    //for(unsigned i = 0;i<26;++i)
      //      test_measureSolvingMountingPlatePatch(i);
    //test_measureSolvingMountingPlate();
    //test_measureSolvingHatPatch(16);
    //test_measureSolvingQuestionMark();
    //test_measureSolvingQuestionMarkPatch(16);
    //printTemplate();
    //for(unsigned i = 0;i<=32;++i)
    //    test_measureSolvingQuestionMarkPatch(i);


    //std::cout << plot << timeTest << nsamples << std::endl;
    return 0;
}
