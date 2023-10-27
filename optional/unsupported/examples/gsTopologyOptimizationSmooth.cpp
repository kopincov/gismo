#include <gismo.h>
#include <gismo_dev.h>

#include <gsUtils/gsTopologyGraph.h>
#include <algorithm>
#include <math.h>

#include <gsSmoothPatches/gsCompositeAssemblerUtils.h>
#include <gsSmoothPatches/gsCompositeUtils.h>

#include "gsSmoothPatches/gsCompositeBSplineBasis.h"
#include "gsSmoothPatches/gsCompositeIncrSmoothnessGeom.h"

#include "gsMSplines/gsMappedSpline.h"
#include "gsMSplines/gsMappedBasis.h"

#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsRecipeAssemblerAdaptive/gsAdaptiveSolver.h>
#include <gsRecipeAssemblerAdaptive/gsCompositeBasisSpaceRefiners.h>
#include <gsRecipeAssemblerAdaptive/gsTensorBasisSpaceRefiners.h>

#include <gsMSplines/gsWeightMapper.h>

using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::pair;
using std::make_pair;
using std::ios;
using namespace gismo;
using namespace std;

typedef Graph::PsCouple PsCouple;
typedef Graph::Couple Couple;

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
        const std::vector<gsPhysicalSpace*>  &space,
        const std::vector<gsMatrix<real_t> > &coefs
        )
    {
        calculateLocalErrors();
        calculateTotalErrors();
        std::cout << "total error: " << m_total << std::endl;
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
    gsAdaptiveSolverOptimizeGeometry(gsCompositeIncrSmoothnessGeom<2,real_t> * compGeom,gsQualityMeasureWeights& weights,bool fixBoundaries,unsigned steps)
        : gsAdaptiveSolver(),m_areaCalculator(NULL),m_weights(weights),m_fixBoundaries(fixBoundaries)
    {
        real_t maxError = 0.0001;
        m_domain=compGeom->exportToPatches();
        gsMappedBasis<2,real_t>& compBasis = compGeom->getCompBasis();
        std::vector<gsBasis<real_t>*> bases;
        for(size_t i =0;i<compBasis.nPatches();++i)
            bases.push_back(&compBasis.getBase(i));
        m_isFittingProblem = m_weights.m_points.size()!=0 && m_weights.m_points[0].size()!=0;
        gsSparseSolver<real_t>::BiCGSTABILUT *eigensolver = new gsSparseSolver<real_t>::BiCGSTABILUT();
//        gsEigenBiCGSTABILUT<real_t>* eigensolver = new gsEigenBiCGSTABILUT<real_t>();

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

        m_criteria       = new gsStopCriteriaIterationOnly(steps);

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
                std::cout << "step " << steps+1 << " damping. \n";
                steps++;
                dampingFactor/=2;
            }
            else
            {
                steps=0;
                dampingFactor=1;
                m_curCoefs=m_spaceRefinerAO->getCurCoefs();
                m_curVal=vals(0,0);
                //std::cout << "curVal " << m_curVal << " . \n";
                if(m_criteria->stop())
                    break;
                if(steps>=m_dampingSteps&&m_dampingSteps!=0)
                {
                    // damping failed
                    std::cout << "\nDamping Failed: Aborting...\n";
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

void printEdges(std::vector<Couple> edges)
{
    std::cout << "edges: ";
    for(unsigned i = 0;i<edges.size();++i)
        std::cout << edges[i].first << " " << edges[i].second << "\n";
    std::cout << std::endl;
}

void printVector(std::vector<unsigned> vec)
{
    std::cout << "vector: ";
    for(unsigned i = 0;i<vec.size();++i)
        std::cout << vec[i] <<" ";
    std::cout << std::endl;
}

bool isCorner(int patch1, int patch2, Graph g,std::vector<unsigned>& sides)
{
    size_t totalNr=0;
    if(sides.size()==1)
    {
        return false;
    }
    for(size_t i = 0;i<sides.size();++i)
    {
        for(size_t j = 0;j<sides[i]-1;++j)
        {
            std::vector<size_t> conn1 = g.vertexConnections(totalNr);
            std::vector<size_t> conn2 = g.vertexConnections(totalNr+1);
            GISMO_ASSERT(conn1.size()==1 && conn2.size()==1,
                         "Every boundary vertex can be only connected to one inner vertex");
            if( (conn1[0]==static_cast<unsigned>(patch1) && conn2[0]==static_cast<unsigned>(patch2)) ||
                    (conn1[0]==static_cast<unsigned>(patch2) && conn2[0]==static_cast<unsigned>(patch1)) )
                return false;
            totalNr++;
        }
        totalNr++;
    }
    return true;
}

bool calculateSmoothStartingSol(gsMultiPatch<real_t>& mp,gsCompositeIncrSmoothnessGeom<2,real_t>& compGeom,Graph g,std::vector<unsigned>& sides,int smoothDeg)
{
    // make list of all boundary corners and give it
    std::vector<patchCorner> corners;
    std::vector<patchSide> bounds = mp.boundaries();
    for(size_t k = 0;k<bounds.size();++k)
    {
        std::vector<boxCorner> boxCorners;
        bounds[k].side().getContainedCorners(2,boxCorners);
        for(size_t j=0;j<boxCorners.size();++j)
        {
            patchCorner c(bounds[k].patch,boxCorners[j]);
            corners.push_back(c);
        }
    }
    std::vector<patchCorner> actualCorners;
    std::vector<patchCorner> cornerList;
    for(size_t k=0;k<corners.size();++k)
    {
        cornerList.clear();
        patchCorner c=corners[k];
        mp.getCornerList(c,cornerList);
        if(cornerList.size()>2)
            actualCorners.push_back(c);
        else if(cornerList.size()==2)
        {
            int patch1=cornerList[0].patch;
            int patch2=cornerList[1].patch;
            if(isCorner(patch1,patch2,g,sides))
                    actualCorners.push_back(c);
        }
    }
    //mp.uniformRefine();
    compGeom = gsCompositeIncrSmoothnessGeom<2,real_t>(mp,actualCorners,smoothDeg);

    gsQualityMeasureWeights weights(mp.nPatches());
    weights.m_length=1;
    weights.m_orthogonality=0;
    weights.m_uniformity=1;
    weights.m_skewness=0;
    weights.m_eccentricity=0;
    weights.m_area=0;
    weights.m_selfIntersect=0;
    weights.m_areaInverse=0;
    weights.m_epsilon=1e-10;
    weights.m_weightArea=false;
    weights.m_approxPoints=0;

    gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,true,1);
    adaptiveSolver.adaptiveSolve();
    return true;
}

bool jacobianTest(gsMultiPatch<real_t>& mp)
{
    bool testPassed=true;
    gsVector<real_t,2> start,end;
    start << 0,0;
    end << 1,1;
    gsVector<unsigned,2> np;
    np << 100,100;
    gsMatrix<real_t> u = gsPointGrid<real_t>(start,end,np);
    gsMatrix<real_t> res;
    real_t det;
    for(size_t i = 0;i<mp.nPatches();++i)
    {
        mp.patch(i).jacobian_into(u,res);
        for(int j = 0;j<res.cols()/2;++j)
        {
            det = res(0,2*j)*res(1,2*j+1)-res(1,2*j)*res(0,2*j+1);
            if(det<=0)
            {
                testPassed=false;
                //std::cout << "Jacobian test failed for patch "<<i<<".\n";
                break;
            }
        }
    }
    return testPassed;
}

bool valid(Graph g,std::vector<unsigned>& sides)
{
    size_t totalNr=0;
    for(size_t i = 0;i<sides.size();++i)
    {
        for(size_t j = 0;j<sides[i]-1;++j)
        {
            std::vector<size_t> conn1 = g.vertexConnections(totalNr);
            std::vector<size_t> conn2 = g.vertexConnections(totalNr+1);
            GISMO_ASSERT(conn1.size()==1 && conn2.size()==1,
                         "Every boundary vertex can be only connected to one inner vertex");
            if(conn1[0]==conn2[0])
                return false;
            totalNr++;
        }
        totalNr++;
    }
    if(sides.size()==1)
    {
        std::vector<size_t> conn1 = g.vertexConnections(sides[0]-1);
        std::vector<size_t> conn2 = g.vertexConnections(0);
        GISMO_ASSERT(conn1.size()==1 && conn2.size()==1,
                     "Every boundary vertex can be only connected to one inner vertex");
        if(conn1[0]==conn2[0])
            return false;
    }
    return true;
}

std::pair<std::pair<int,int>,real_t> findBestGeom(std::vector<gsMatrix<real_t> >& boundaryMatrices,std::vector<unsigned>& sides,std::vector<unsigned>& smoothSides,
                                          unsigned minFaces,unsigned maxFaces,std::string name,
                                          gsQualityMeasureWeights* givenWeights = NULL,
                                          int smoothDeg=0,int refine=0, bool print = false)
{
    size_t boundarySize = boundaryMatrices.size();
    std::string file = "TopologyMatrizes/Top"+util::to_string(boundarySize)+".xml";

    std::vector<Graph> graphs;
    Graph::readGraphListFromFile(gsFileManager::find(file),graphs);
    std::vector<std::vector<Graph> > graphPartitions;
    Graph::getPartitionsOfGraphVector(graphs,graphPartitions,minFaces+boundarySize,maxFaces+boundarySize);

    std::vector<std::vector<gsMultiPatch<real_t> > >mps;
    std::vector<gsMultiPatch<real_t> >mps_temp;
    std::vector<std::vector<bool> >valids;
    std::vector<bool> valids_temp;
    for(size_t i = 0;i<graphPartitions.size();++i)
    {
        mps_temp.clear();
        ConstructMultiPatchFromGraph::graphsToMPs(graphPartitions[i],boundaryMatrices,mps_temp);
        mps.push_back(mps_temp);
        valids_temp.clear();
        for(size_t j = 0;j<graphPartitions[i].size();++j)
        {
            valids_temp.push_back(valid(graphPartitions[i][j],sides));
        }
        valids.push_back(valids_temp);
    }
    int curmin = -1;
    int curpatch = -1;
    real_t curval = 100000000;
    unsigned toposConsidered=0;
    for(size_t i = 0;i<mps.size();++i)
    {
        for(size_t j = 0;j<mps[i].size();++j)
        {
            //std::cout << " i: " <<i << " j: " <<j<<std::endl;
            if(6!=mps[i][j].nPatches()||j!=24)
                continue;
            if(!valids[i][j])
                continue;
            toposConsidered++;
            gsCompositeIncrSmoothnessGeom<2,real_t> compGeom;
            if(refine>0)
                mps[i][j].uniformRefine(refine);
            calculateSmoothStartingSol(mps[i][j],compGeom,graphPartitions[i][j],smoothSides,smoothDeg);
            if(j==156&&i==6&&false)
            {
                std::string filename=name+"_6_"+util::to_string(j)+"_start3";
                gsWriteParaview(compGeom.exportToPatches(),filename,1000,false,false);
            }
            gsQualityMeasureWeights weights(mps[i][j].nPatches());
            if(givenWeights!=NULL)
            {
                weights.m_length=givenWeights->m_length;
                weights.m_orthogonality=givenWeights->m_orthogonality;
                weights.m_uniformity=givenWeights->m_uniformity;
                weights.m_skewness=givenWeights->m_skewness;
                weights.m_eccentricity=givenWeights->m_eccentricity;
                weights.m_area=givenWeights->m_area;
                weights.m_selfIntersect=givenWeights->m_selfIntersect;
                weights.m_areaInverse=givenWeights->m_areaInverse;
                weights.m_epsilon=givenWeights->m_epsilon;
                weights.m_weightArea=givenWeights->m_weightArea;
                weights.m_approxPoints=givenWeights->m_approxPoints;
            }
            else
            {
                weights.m_length=0;
                weights.m_orthogonality=0.5; //0.1;
                weights.m_uniformity=0.1;
                weights.m_skewness=0;
                weights.m_eccentricity=0;
                weights.m_area=2; //1;
                weights.m_selfIntersect=0;
                weights.m_areaInverse=0;
                weights.m_epsilon=1e-10;
                weights.m_weightArea=false;
                weights.m_approxPoints=0;
            }

            gsAdaptiveSolverOptimizeGeometry adaptiveSolver(&compGeom,weights,true,25);
            adaptiveSolver.adaptiveSolve();
            real_t val = adaptiveSolver.getCurVal();
            gsMultiPatch<> mp = compGeom.exportToPatches();
            bool regular = jacobianTest(mp);
            if(print)
            {
                std::string filename=name+"_"+util::to_string(mp.nPatches())+"_"+util::to_string(j);
                if(!regular)
                    filename=filename+"_nr";
                if(true||regular)
                {
                    gsWriteParaview(mp,filename,500,false,false);
                    std::cout << filename << " val: " << val << std::endl;
                }
//                if(mp->nPatches()==9&&j==7353)
//                {
//                    gsFileData<> fd;
//                    fd<< mp ;
//                    fd.dump(filename);
//                }
            }

            if( val<curval && regular )
            {
                curval=val;
                curmin=j;
                curpatch=mp.nPatches();
            }
        }
    }
    std::cout << "Topologies considered: "  << toposConsidered << std::endl;
    return std::make_pair(std::make_pair(curpatch,curmin),curval);
}

void scaleToSquare(std::vector<gsMatrix<real_t> >& boundaryMatrices,real_t size=1)
{
    real_t xmin=1000000000,ymin=1000000000,xmax=-1000000000,ymax=-1000000000;
    real_t xval,yval;
    for(size_t i = 0;i<boundaryMatrices.size();++i)
    {
        for(index_t j = 0;j<boundaryMatrices[i].rows();++j)
        {
            xval=boundaryMatrices[i](j,0);
            yval=boundaryMatrices[i](j,1);
            if(xval<xmin)
                xmin=xval;
            if(xval>xmax)
                xmax=xval;
            if(yval<ymin)
                ymin=yval;
            if(yval>ymax)
                ymax=yval;
        }
    }
    real_t xUnitFactor=xmax-xmin;
    real_t yUnitFactor=ymax-ymin;
    real_t factor = xUnitFactor>yUnitFactor ? xUnitFactor : yUnitFactor;
    factor/=size;
    std::cout << "Factor: " << factor << std::endl;
    for(size_t i = 0;i<boundaryMatrices.size();++i)
    {
        for(index_t j = 0;j<boundaryMatrices[i].rows();++j)
        {
            xval=boundaryMatrices[i](j,0);
            yval=boundaryMatrices[i](j,1);
            xval=(xval-xmin)/factor;
            yval=(yval-ymin)/factor;
            boundaryMatrices[i](j,0)=xval;
            boundaryMatrices[i](j,1)=yval;
        }
    }
}

void comparisonTripplet()
{
    std::vector<string> descs;
    std::vector<int> mins_nr;
    std::vector<int> mins_patch;
    std::vector<real_t> mins_val;
    int smoothDeg=1;
    for(unsigned cases = 1;cases<=3;++cases)
    {
        unsigned maxFaces=6;
        std::vector<unsigned> sides;
        std::vector<unsigned> smoothSides;
        string desc;
        std::vector<gsMatrix<real_t> > boundaryMatrices;
        switch(cases)
        {
        case 1: // hexagon example
        {
            desc="Hexagon";
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            smoothSides=sides;
            gsMatrix<real_t> coefs(5,2);
            coefs<< 0, 4, 1, 3, 2, 2, 3, 1, 4, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 4, 0, 5, 0, 6, 0, 7, 0, 8, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 8, 0, 9, 1, 10, 2, 11, 3, 12, 4;
            boundaryMatrices.push_back(coefs);
            coefs<< 12, 4, 11, 5, 10, 6, 9, 7, 8, 8;
            boundaryMatrices.push_back(coefs);
            coefs<< 8, 8, 7, 8, 6, 8, 5, 8, 4, 8;
            boundaryMatrices.push_back(coefs);
            coefs<< 4, 8, 3, 7, 2, 6, 1, 5, 0, 4;
            boundaryMatrices.push_back(coefs);
            break;
        }
        case 2: // tunnel example
        {
            desc="Tunnel";
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            smoothSides=sides;
            gsMatrix<real_t> coefs(5,2);
            coefs<< 0, 0, 1, 0, 2, 0, 3, 0, 4, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 4, 0, 4, 2, 6, 3, 8, 2, 8, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 8, 0, 9, 0, 10, 0, 11, 0, 12, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 12, 0, 12, 1, 12, 2, 12, 3, 12, 4;
            boundaryMatrices.push_back(coefs);
            coefs<< 12, 4, 9, 4, 6, 4, 3, 4, 0, 4;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 4, 0, 3, 0, 2, 0, 1, 0, 0;
            boundaryMatrices.push_back(coefs);
            break;
        }
        case 3: // shoe example
        {
            desc="Yacht";
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            smoothSides=sides;
            gsMatrix<real_t> coefs(5,2);
            coefs<< 0,0,2,0,4,0,6,0,8,0;
            boundaryMatrices.push_back(coefs);
            coefs<< 8,0,9,0.5,10,1,11,2,12,4;
            boundaryMatrices.push_back(coefs);
            coefs<< 12,4,10,3.5,9,3.5,7,4.5,6,5;
            boundaryMatrices.push_back(coefs);
            coefs<< 6,5,5,4,4,3.5,3,4,2,4.5;
            boundaryMatrices.push_back(coefs);
            coefs<< 2,4.5,1.5,4.5,1,4.5,0.5,4.5,0,4.5;
            boundaryMatrices.push_back(coefs);
            coefs<< 0,4.5,0,3.5,0,2,0,1,0,0;
            boundaryMatrices.push_back(coefs);
            break;
        }
        default:
            return;
        }
        gsQualityMeasureWeights weights(1);
        weights.m_length=0;//0.1;
        weights.m_orthogonality=0.1;//0.5; //0.1;
        weights.m_uniformity=0.1;
        weights.m_skewness=0;
        weights.m_eccentricity=0;
        weights.m_area=2; //1;
        weights.m_selfIntersect=0;
        weights.m_areaInverse=0;
        weights.m_epsilon=1e-10;
        weights.m_weightArea=false;
        weights.m_approxPoints=0;
        //scaleToUnitSquare(boundaryMatrices);
        std::pair<std::pair<int,int>,real_t> res = findBestGeom(boundaryMatrices,sides,smoothSides,1, maxFaces,desc,&weights,smoothDeg,0,true);
        descs.push_back(desc);
        mins_nr.push_back(res.first.second);
        mins_patch.push_back(res.first.first);
        mins_val.push_back(res.second);

    }
    for(unsigned cases = 0;cases<descs.size();++cases)
    {
        std::cout << descs[cases] << ": Patch " <<mins_patch[cases]<<", nr: "<< mins_nr[cases] << " with val " << mins_val[cases] <<" \n";
    }
}

void comparisonSlightChangeOfGeometry()
{
    std::vector<string> descs;
    std::vector<int> mins_nr;
    std::vector<int> mins_patch;
    std::vector<real_t> mins_val;
    int smoothDeg=1;
    for(unsigned cases = 1;cases<=2;++cases)
    {
        unsigned maxFaces=4;
        std::vector<unsigned> sides;
        std::vector<unsigned> smoothSides;
        string desc;
        std::vector<gsMatrix<real_t> > boundaryMatrices;
        switch(cases)
        {
        case 1: // plate with hole example
        {
            desc="PlateWithHole";
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(2);
            sides.push_back(1);
            smoothSides=sides;
            gsMatrix<real_t> coefs(5,2);
            coefs<< 0, 8, 2, 8, 4, 7.5, 5.5, 6.5, 6, 6;
            boundaryMatrices.push_back(coefs);
            coefs<< 6, 6, 6.5, 5.5, 7.5, 4, 8, 2, 8, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 8, 0, 10, 0, 12, 0, 14, 0, 16, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 16, 0, 16, 3, 16, 5, 16, 8, 16, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 16, 10, 16, 13, 16, 14, 16, 15, 16, 16;
            boundaryMatrices.push_back(coefs);
            coefs<< 16, 16, 15, 16, 14, 16, 13, 16, 10, 16;
            boundaryMatrices.push_back(coefs);
            coefs<< 10, 16, 8, 16, 5, 16, 3, 16, 0, 16;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 16, 0, 14, 0, 12, 0, 10, 0, 8;
            boundaryMatrices.push_back(coefs);
            break;
        }
        case 2: // plate with hole skewed example
        {
            desc="PlateWithHoleSkewed";
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(2);
            sides.push_back(1);
            smoothSides=sides;
            gsMatrix<real_t> coefs(5,2);
            coefs<< 0, 8, 2, 8, 4, 7.5, 5.5, 6.5, 6, 6;
            boundaryMatrices.push_back(coefs);
            coefs<< 6, 6, 6.5, 5.5, 7.5, 4, 8, 2, 8, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 8, 0, 10, 0, 12, 0, 14, 0, 16, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 16, 0, 16, 1, 16, 2, 16, 3, 16, 4;
            boundaryMatrices.push_back(coefs);
            coefs<< 16, 4, 16, 5, 16, 7, 16, 8, 16, 9;
            boundaryMatrices.push_back(coefs);
            coefs<< 16, 9, 15, 9, 14, 9, 13, 9, 10, 9;
            boundaryMatrices.push_back(coefs);
            coefs<< 10, 9, 8, 9, 5, 9, 3, 9, 0, 9;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 9, 0, 8.75, 0, 8.5, 0, 8.25, 0, 8;
            boundaryMatrices.push_back(coefs);
           break;
        }
        default:
            return;
        }
        gsQualityMeasureWeights weights(1);
        weights.m_length=0;//0.1;
        weights.m_orthogonality=5;//0.5; //0.1;
        weights.m_uniformity=1;
        weights.m_skewness=0;
        weights.m_eccentricity=0;
        weights.m_area=0.1; //1;
        weights.m_selfIntersect=0;
        weights.m_areaInverse=0;
        weights.m_epsilon=1e-10;
        weights.m_weightArea=false;
        weights.m_approxPoints=0;
        std::pair<std::pair<int,int>,real_t> res = findBestGeom(boundaryMatrices,sides,smoothSides,1, maxFaces,desc,&weights,smoothDeg,0,true);
        descs.push_back(desc);
        mins_nr.push_back(res.first.second);
        mins_patch.push_back(res.first.first);
        mins_val.push_back(res.second);

    }
    for(unsigned cases = 0;cases<descs.size();++cases)
    {
        std::cout << descs[cases] << ": Patch " <<mins_patch[cases]<<", nr: "<< mins_nr[cases] << " with val " << mins_val[cases] <<" \n";
    }
}

void examples()
{
    std::vector<string> descs;
    std::vector<int> mins_nr;
    std::vector<int> mins_patch;
    std::vector<real_t> mins_val,times;
    int smoothDeg=1;
    for(unsigned cases = 5;cases<=5;++cases)
    {
        unsigned maxFaces=0;
        std::vector<unsigned> sides;
        std::vector<unsigned> smoothSides;
        string desc;
        std::vector<gsMatrix<real_t> > boundaryMatrices;
        gsQualityMeasureWeights weights(1);
        switch(cases)
        {
        case 1: // blade example
        {
            desc="Blade";
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            smoothSides=sides;
            gsMatrix<real_t> coefs(5,2);
            coefs<< 0,0,1,-1,2,-2,3,-3,4,-4;
            boundaryMatrices.push_back(coefs);
            coefs<< 4,-4,4,-5,4,-6,4,-7,4,-8;
            boundaryMatrices.push_back(coefs);
            coefs<< 4,-8,6,-6,8,-4,10,-2,12,0;
            boundaryMatrices.push_back(coefs);
            coefs<< 12,0,10,2,8,4,6,6,4,8;
            boundaryMatrices.push_back(coefs);
            coefs<< 4,8,4,7,4,6,4,5,4,4;
            boundaryMatrices.push_back(coefs);
            coefs<< 4,4,3,3,2,2,1,1,0,0;
            boundaryMatrices.push_back(coefs);

            //weights.m_orthogonality=1; //0.1;
            //weights.m_uniformity=0.1;
            //weights.m_area=1.5; //1;

            weights.m_orthogonality=0.5; //0.1;
            weights.m_uniformity=0.1;
            weights.m_area=2; //1;

            maxFaces=7;
            break;
        }
        case 2: // hammer example
        {
            desc="Hammer";
            sides.push_back(2);
            sides.push_back(3);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(3);
            smoothSides.push_back(2);
            smoothSides.push_back(2);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(2);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(2);
            gsMatrix<real_t> coefs(5,2);
            coefs<< 2,0,2.5,0,3,0,3.5,0,4,0;
            boundaryMatrices.push_back(coefs);
            coefs<< 4,0,4.5,0,5,0,5.5,0,6,0;
            boundaryMatrices.push_back(coefs);
            coefs<< 6,0,6,1,6,2,6,3,6,4;
            boundaryMatrices.push_back(coefs);
            coefs<< 6,4,6,5,6,6,6,7,6,8;
            boundaryMatrices.push_back(coefs);
            coefs<< 6,8,7,8,8,8,9,8,10,8;
            boundaryMatrices.push_back(coefs);
            coefs<< 10,8,9.5,9,9,10,8.5,11,8,12;
            boundaryMatrices.push_back(coefs);
            coefs<< 8,12,7,12,6,12,5,12,4,12;
            boundaryMatrices.push_back(coefs);
            coefs<< 4,12,3,12,2,12,1,12,0,12;
            boundaryMatrices.push_back(coefs);
            coefs<< 0,12,0,11,0,10,0,9,0,8;
            boundaryMatrices.push_back(coefs);
            coefs<< 0,8,0.5,8,1,8,1.5,8,2,8;
            boundaryMatrices.push_back(coefs);
            coefs<< 2,8,2,7,2,6,2,5,2,4;
            boundaryMatrices.push_back(coefs);
            coefs<< 2,4,2,3,2,2,2,1,2,0;
            boundaryMatrices.push_back(coefs);

            weights.m_orthogonality=0.5; //0.1;
            weights.m_uniformity=0.1;
            weights.m_area=2; //1;

            maxFaces=7;
            break;
        }
        case 3: // ice cone example
        {
            desc="IceCone";
            sides.push_back(6);
            smoothSides.push_back(1);
            smoothSides.push_back(3);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            gsMatrix<real_t> coefs(5,2);
            coefs<< 0.8,10.2,0.7,10.3,-0.2,9.7,-0.5,8.5,0,8;
            boundaryMatrices.push_back(coefs);
            coefs<< 0,8,0.25,7,0.5,6,0.75,5.75,0.8,5.5;
            boundaryMatrices.push_back(coefs);
            coefs<< 0.8,5.5,0.85,5.25,2,5.25,3.15,5.25,3.2,5.5;
            boundaryMatrices.push_back(coefs);
            coefs<< 3.2,5.5,3.25,5.75,3.5,6,3.75,7,4,8;
            boundaryMatrices.push_back(coefs);
            coefs<< 4,8,4.5,8.5,4.2,9.7,3.3,10.3,3.2,10.2;
            boundaryMatrices.push_back(coefs);
            coefs<< 3.2,10.2,3,10.5,2,10.75,1,10.5,0.8,10.2;
            boundaryMatrices.push_back(coefs);


            //weights.m_orthogonality=0.5; //0.1;
            //weights.m_uniformity=0.1;
            //weights.m_eccentricity=1;
            weights.m_area=5; //1;
            weights.m_length=1;//0.1;

            maxFaces=8;
            break;
        }
        case 4: // car example
        {
            desc="Car";
            sides.push_back(5);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            smoothSides.push_back(1);
            gsMatrix<real_t> coefs(5,2);
            coefs<< 0,0,0.5,0,1,0,1.5,0,2,0;
            boundaryMatrices.push_back(coefs);
            coefs<< 2,0,2,0.5,3,0.9,4,0.5,4,0;
            boundaryMatrices.push_back(coefs);
            coefs<< 4,0,5,0,6,0,7,0,8,0;
            boundaryMatrices.push_back(coefs);
            coefs<< 8,0,8,0.5,9,0.9,10,0.5,10,0;
            boundaryMatrices.push_back(coefs);
            coefs<< 10,0,10.5,0,11,0,11.5,0,12,0;
            boundaryMatrices.push_back(coefs);
            coefs<< 12,0,12,0.5,12,1,12,1.5,12,2;
            boundaryMatrices.push_back(coefs);
            coefs<< 12,2,11.5,2.1,10.5,2.2,9.5,2.3,9,2.4;
            boundaryMatrices.push_back(coefs);
            coefs<< 9,2.4,8.5,2.8,8,3.2,7.5,3.6,7,4;
            boundaryMatrices.push_back(coefs);
            coefs<< 7,4,6.5,4,4.5,4,3.5,4,3,4;
            boundaryMatrices.push_back(coefs);
            coefs<< 3,4,2.75,3.6,2.5,3.2,2.25,2.8,2,2.4;
            boundaryMatrices.push_back(coefs);
            coefs<< 2,2.4,1.5,2.3,1,2.2,0.5,2.1,0,2;
            boundaryMatrices.push_back(coefs);
            coefs<< 0,2,0,1.5,0,1,0,0.5,0,0;
            boundaryMatrices.push_back(coefs);


            weights.m_orthogonality=0.5; //0.1;
            weights.m_uniformity=0.1;
            //weights.m_eccentricity=1;
            weights.m_area=2; //1;
            //weights.m_length=1;//0.1;

            maxFaces=8;
            break;
        }
        case 5: // puzzle Piece example inward
        {
            desc="Puzzle Piece IW";
            sides.push_back(3);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(3);
            smoothSides=sides;
            gsMatrix<real_t> coefs(6,2);
            coefs<<0,12,0,11,0,10,0,8,0,7,0,6;
            boundaryMatrices.push_back(coefs);
            coefs<<0,6,0,5,0,4,0,3,0,2,0.5,1.5;
            boundaryMatrices.push_back(coefs);
            coefs<<0.5,1.5,1,1,3,2.5,5,2.5,6,1,6,0;
            boundaryMatrices.push_back(coefs);
            coefs<<6,0,7,0,8,0,10,0,11,0,12,0;
            boundaryMatrices.push_back(coefs);
            coefs<<12,0,12,1,12,2,12,4,12,5,12,6;
            boundaryMatrices.push_back(coefs);
            coefs<<12,6,11,6,9.5,7,9.5,9,11,11,10.5,11.5;
            boundaryMatrices.push_back(coefs);
            coefs<<10.5,11.5,10,12,9,12,8,12,7,12,6,12;
            boundaryMatrices.push_back(coefs);
            coefs<<6,12,5,12,4,12,2,12,1,12,0,12;
            boundaryMatrices.push_back(coefs);

            //weights.m_orthogonality=0.5; //0.1;
            weights.m_uniformity=0.1;
            //weights.m_eccentricity=1;
            weights.m_area=2; //1;
            //weights.m_length=1;//0.1;

            maxFaces=6;
            break;
        }
        case 6: // puzzle Piece example outward
        {
            desc="Puzzle Piece OW";
            sides.push_back(3);
            sides.push_back(3);
            sides.push_back(3);
            sides.push_back(3);
            smoothSides=sides;
            gsMatrix<real_t> coefs(6,2);
            coefs<<0,12,0,11,0,10,0,8,0,7,0,6;
            boundaryMatrices.push_back(coefs);
            coefs<<0,6,0,5,0,4,0,3,0,2,-0.5,1.5;
            boundaryMatrices.push_back(coefs);
            coefs<<-0.5,1.5,-1,1,-3,2.5,-5,2.5,-6,1,-6,0;
            boundaryMatrices.push_back(coefs);
            coefs<<-6,0,-4,0,-2,0,0,0,2,0,4,0;
            boundaryMatrices.push_back(coefs);
            coefs<<4,0,5,0,6,0,7,0,8,0,9,0;
            boundaryMatrices.push_back(coefs);
            coefs<<9,0,9.5,0,10,0,11,0,11.5,0,12,0;
            boundaryMatrices.push_back(coefs);
            coefs<<12,0,12,0.5,12,1,12,2,12,2.5,12,3;
            boundaryMatrices.push_back(coefs);
            coefs<<12,3,12,4,12,5,12,6,12,7,12,8;
            boundaryMatrices.push_back(coefs);
            coefs<<12,8,12,10,12,12,12,14,12,16,12,18;
            boundaryMatrices.push_back(coefs);
            coefs<<12,18,11,18,9.5,17,9.5,15,11,13,10.5,12.5;
            boundaryMatrices.push_back(coefs);
            coefs<<10.5,12.5,10,12,9,12,8,12,7,12,6,12;
            boundaryMatrices.push_back(coefs);
            coefs<<6,12,5,12,4,12,2,12,1,12,0,12;
            boundaryMatrices.push_back(coefs);

            weights.m_uniformity=0.01;
            weights.m_area=15; //1;

            maxFaces=9;
            break;
        }
        default:
            return;
        }
        scaleToSquare(boundaryMatrices,12);
        gsStopwatch time;
        std::pair<std::pair<int,int>,real_t> res = findBestGeom(boundaryMatrices,sides,smoothSides,maxFaces, maxFaces,desc,&weights,smoothDeg,0,true);
        times.push_back(time.stop());
        descs.push_back(desc);
        mins_nr.push_back(res.first.second);
        mins_patch.push_back(res.first.first);
        mins_val.push_back(res.second);
    }
    for(unsigned cases = 0;cases<descs.size();++cases)
    {
        std::cout << descs[cases] << ": Patch " <<mins_patch[cases]<<", nr: "<< mins_nr[cases] << " with val " << mins_val[cases] << " in " << times[cases] <<"sec. \n";
    }
}

void twoBlades()
{
    std::vector<string> descs;
    std::vector<int> mins_nr;
    std::vector<int> mins_patch;
    std::vector<real_t> mins_val,times;
    int smoothDeg=1;
    for(unsigned cases = 1;cases<=1;++cases)
    {
        unsigned maxFaces=9;
        std::vector<unsigned> sides;
        std::vector<unsigned> smoothSides;
        string desc;
        std::vector<gsMatrix<real_t> > boundaryMatrices;
        switch(cases)
        {
        case 1: // Blade1 example
        {
            desc="Blade1";
            sides.push_back(1);
            sides.push_back(3);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(3);
            sides.push_back(1);
            sides.push_back(1);
            smoothSides=sides;
            gsMatrix<real_t> coefs(4,2);
            coefs<< 0, 0, 1, 1, 2, 2, 3, 3;
            boundaryMatrices.push_back(coefs);
            coefs<< 3, 3, 2.8, 3.2, 5, 6.1, 6, 7.2;
            boundaryMatrices.push_back(coefs);
            coefs<< 6, 7.2, 7, 8.3, 8.5, 9.5, 10, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 10, 10, 11.5, 10.5, 13, 10.2, 13, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 10, 14, 10, 16, 10, 17, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 17, 10, 17, 12, 17, 16, 17, 18;
            boundaryMatrices.push_back(coefs);
            coefs<< 17, 18, 16, 18, 14, 18, 13, 18;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 18, 13, 17.8, 12.5, 17, 11.7, 16.2;
            boundaryMatrices.push_back(coefs);
            coefs<< 11.7, 16.2, 11.324167, 15.824169, 8.889348, 14.82165, 8.2, 14.7;
            boundaryMatrices.push_back(coefs);
            coefs<< 8.2, 14.7, 6.5, 14.4, 3.2, 10.8, 3, 11;
            boundaryMatrices.push_back(coefs);
            coefs<< 3, 11, 2, 10, 1, 9, 0, 8;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 8, 0, 6, 0, 2, 0, 0;
            boundaryMatrices.push_back(coefs);
            break;
        }
        case 2: // Blade2 example
        {
            desc="Blade2";
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(1);
            smoothSides=sides;
            gsMatrix<real_t> coefs(4,2);
            coefs<< 0, 0, 1.5, 1.5, 3, 3, 4.3, 4.3;
            boundaryMatrices.push_back(coefs);
            coefs<< 4.3, 4.3, 4.2, 4.4, 6.2, 7.8, 8, 9.9;
            boundaryMatrices.push_back(coefs);
            coefs<< 8, 9.9, 9.8, 12, 11.9, 14.1, 12, 14;
            boundaryMatrices.push_back(coefs);
            coefs<< 12, 14, 12.25, 14.25, 12.75, 14.75, 13, 15;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 15, 13, 18, 13, 20, 13, 23;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 23, 12.75, 22.75, 12.25, 22.25, 12, 22;
            boundaryMatrices.push_back(coefs);
            coefs<< 12, 22, 12.1, 21.9, 9.6, 19.3, 8, 17.2;
            boundaryMatrices.push_back(coefs);
            coefs<< 8, 17.2, 6.4, 15.1, 4.4, 12.2, 4.3, 12.3;
            boundaryMatrices.push_back(coefs);
            coefs<< 4.3, 12.3, 3, 11, 1.5, 9.5, 0, 8;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 8, 0, 5, 0, 3, 0, 0;
            boundaryMatrices.push_back(coefs);
            break;
        }
        case 3: // Blade1a example
        {
            desc="Blade1a";
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(1);
            smoothSides=sides;
            gsMatrix<real_t> coefs(5,2);
            coefs<< 0, 0, 0.5, 0.5,1.5,1.5, 2.5, 2.5, 3, 3;
            boundaryMatrices.push_back(coefs);
            coefs<< 3, 3, 2.8, 3.2, 5, 6.1,7, 8.3,7.75,8.9;
            boundaryMatrices.push_back(coefs);
            coefs<< 7.75,8.9, 8.5,9.5,11.5, 10.5, 13, 10.2, 13, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 10, 14, 10,15,10, 16, 10, 17, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 17, 10, 17, 12,17,14, 17, 16, 17, 18;
            boundaryMatrices.push_back(coefs);
            coefs<< 17, 18, 16, 18,15,18, 14, 18, 13, 18;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 18, 13, 17.8, 12.5, 17, 11.32, 15.6,9,14.9;
            boundaryMatrices.push_back(coefs);
            coefs<< 9,14.9, 7.84, 14.59, 6.5, 14.1, 3.2, 10.8, 3, 11;
            boundaryMatrices.push_back(coefs);
            coefs<< 3, 11, 2.5, 10.5, 1.5, 9.5,0.5,8.5, 0, 8;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 8, 0, 6,0,4, 0, 2, 0, 0;
            boundaryMatrices.push_back(coefs);
            break;
        }
        case 4: // Blade1b example
        {
            desc="Blade1b";
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(2);
            smoothSides=sides;
            gsMatrix<real_t> coefs(5,2);
            coefs<< 0, 0, 0.5, 0.5,1.5,1.5, 2.5, 2.5, 3, 3;
            boundaryMatrices.push_back(coefs);
            coefs<< 3, 3, 2.8, 3.2, 5, 6.1,7, 8.3,7.75,8.9;
            boundaryMatrices.push_back(coefs);
            coefs<< 7.75,8.9, 8.5,9.5,11.5, 10.5, 13, 10.2, 13, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 10, 14, 10,15,10, 16, 10, 17, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 17, 10, 17, 12,17,13, 17, 14, 17, 14;
            boundaryMatrices.push_back(coefs);
            coefs<< 17, 14, 17, 15,17,16, 17, 17, 17, 18;
            boundaryMatrices.push_back(coefs);
            coefs<< 17, 18, 16, 18,15,18, 14, 18, 13, 18;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 18, 13, 17.8, 12.5, 17, 11.32, 15.6,9,14.9;
            boundaryMatrices.push_back(coefs);
            coefs<< 9,14.9, 7.84, 14.59, 6.5, 14.1, 3.2, 10.8, 3, 11;
            boundaryMatrices.push_back(coefs);
            coefs<< 3, 11, 2.5, 10.5, 1.5, 9.5,0.5,8.5, 0, 8;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 8, 0, 7,0,6, 0, 5, 0, 4;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 4, 0, 3,0,2, 0, 1, 0, 0;
            boundaryMatrices.push_back(coefs);
            break;
        }
        case 5: // Blade2a example
        {
            desc="Blade2a";
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(2);
            sides.push_back(1);
            sides.push_back(2);
            smoothSides=sides;
            gsMatrix<real_t> coefs(4,2);
            coefs<< 0, 0, 1.5, 1.5, 3, 3, 4.3, 4.3;
            boundaryMatrices.push_back(coefs);
            coefs<< 4.3, 4.3, 4.2, 4.4, 6.2, 7.8, 8, 9.9;
            boundaryMatrices.push_back(coefs);
            coefs<< 8, 9.9, 9.8, 12, 11.9, 14.1, 12, 14;
            boundaryMatrices.push_back(coefs);
            coefs<< 12, 14, 12.25, 14.25, 12.75, 14.75, 13, 15;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 15, 13, 16, 13, 18, 13, 19;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 19, 13, 20, 13, 22, 13, 23;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 23, 12.75, 22.75, 12.25, 22.25, 12, 22;
            boundaryMatrices.push_back(coefs);
            coefs<< 12, 22, 12.1, 21.9, 9.6, 19.3, 8, 17.2;
            boundaryMatrices.push_back(coefs);
            coefs<< 8, 17.2, 6.4, 15.1, 4.4, 12.2, 4.3, 12.3;
            boundaryMatrices.push_back(coefs);
            coefs<< 4.3, 12.3, 3, 11, 1.5, 9.5, 0, 8;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 8, 0, 7, 0, 5, 0, 4;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 4, 0, 3, 0, 1, 0, 0;
            boundaryMatrices.push_back(coefs);
            break;
        }
        default:
            return;
        }
        gsQualityMeasureWeights weights(1);
        weights.m_length=0;
        weights.m_orthogonality=1;//1;
        weights.m_uniformity=1;//1;
        weights.m_skewness=1;//1;
        weights.m_eccentricity=0;
        weights.m_area=5; //1;
        weights.m_selfIntersect=0;
        weights.m_areaInverse=0;
        weights.m_epsilon=1e-10;
        weights.m_weightArea=false;
        weights.m_approxPoints=0;
        scaleToSquare(boundaryMatrices,12);
        gsStopwatch time;
        std::pair<std::pair<int,int>,real_t> res = findBestGeom(boundaryMatrices,sides,smoothSides,1,maxFaces,desc,&weights,smoothDeg,0,true);
        times.push_back(time.stop());
        descs.push_back(desc);
        mins_nr.push_back(res.first.second);
        mins_patch.push_back(res.first.first);
        mins_val.push_back(res.second);

    }
    for(unsigned cases = 0;cases<descs.size();++cases)
    {
        std::cout << descs[cases] << ": Patch " <<mins_patch[cases]<<", nr: "<< mins_nr[cases] << " with val " << mins_val[cases] << " in " << times[cases] <<"sec. \n";
    }
}
void newExamples()
{
    std::vector<string> descs;
    std::vector<int> mins_nr;
    std::vector<int> mins_patch;
    std::vector<real_t> mins_val,times;
    int smoothDeg=1;
    int refine=0;
    for(unsigned cases = 1;cases<=1;++cases)
    {
        unsigned maxFaces;
        std::vector<unsigned> sides;
        std::vector<unsigned> smoothSides;
        string desc;
        std::vector<gsMatrix<real_t> > boundaryMatrices;
        gsQualityMeasureWeights weights(1);
        switch(cases)
        {
        case 1: // tunnel example
        {
            desc="Tunnel-new";
            sides.push_back(2);
            sides.push_back(2);
            sides.push_back(2);
            sides.push_back(2);
            sides.push_back(2);
            smoothSides=sides;
            refine=1;
            maxFaces=7;
            gsMatrix<real_t> coefs(4,2);
            coefs<< 0, 0, 1, 0, 3, 0, 4, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 4, 0, 5, 0, 7, 0, 8, 0;
            boundaryMatrices.push_back(coefs);
            coefs<< 8, 0, 8, 1, 8.5, 3, 9, 4;
            boundaryMatrices.push_back(coefs);
            coefs<< 9, 4, 9.5, 5, 11, 6, 12, 6;
            boundaryMatrices.push_back(coefs);
            coefs<< 12, 6, 12, 7, 12, 8, 12, 9;
            boundaryMatrices.push_back(coefs);
            coefs<< 12, 9, 12, 10, 12,11, 12, 12;
            boundaryMatrices.push_back(coefs);
            coefs<< 12, 12, 10, 11, 6, 10, 4, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 4, 10, 3, 10, 1, 10, 0, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 10, 0, 9, 0, 7, 0, 6;
            boundaryMatrices.push_back(coefs);
            coefs<< 0, 6, 0, 4.5, 0, 1.5, 0, 0;
            boundaryMatrices.push_back(coefs);
            weights.m_orthogonality=1;//0.5; //0.1;
            weights.m_uniformity=1;
            weights.m_area=1; //1;
            break;
        }
        case 2: // Blade1 example
        {
            desc="Blade1-new";
            sides.push_back(6);
            smoothSides=sides;
            refine=0;
            maxFaces=9;
            gsMatrix<real_t> coefs(4,2);
            coefs<< 3, 3, 3.2, 2.8, 6.5, 6.4, 8.2, 6.7;
            boundaryMatrices.push_back(coefs);
            coefs<< 8.2, 6.7, 8.889348, 6.82165, 11.324167, 7.824169, 11.7, 8.2;
            boundaryMatrices.push_back(coefs);
            coefs<< 11.7, 8.2, 12.5, 9, 13, 9.8, 13, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 13, 10, 13, 10.2, 11.5, 10.5, 10, 10;
            boundaryMatrices.push_back(coefs);
            coefs<< 10, 10,8.5, 9.5, 7, 8.3, 6, 7.2;
            boundaryMatrices.push_back(coefs);
            coefs<< 6, 7.2, 5, 6.1, 2.8, 3.2, 3, 3;
            boundaryMatrices.push_back(coefs);
            weights.m_length=0;//0.1;
            weights.m_orthogonality=0.1;//0.5; //0.1;
            weights.m_uniformity=0.1;
            weights.m_area=2; //1;
            break;
        }
        default:
            return;
        }
        scaleToSquare(boundaryMatrices,12);
        gsStopwatch time;
        std::pair<std::pair<int,int>,real_t> res = findBestGeom(boundaryMatrices,sides,smoothSides,1,maxFaces,desc,&weights,smoothDeg,refine,true);
        times.push_back(time.stop());
        descs.push_back(desc);
        mins_nr.push_back(res.first.second);
        mins_patch.push_back(res.first.first);
        mins_val.push_back(res.second);

    }
    for(unsigned cases = 0;cases<descs.size();++cases)
    {
        std::cout << descs[cases] << ": Patch " <<mins_patch[cases]<<", nr: "<< mins_nr[cases] << " with val " << mins_val[cases] << " in " << times[cases] <<"sec. \n";
    }
}

int main()
{

//   comparisonTripplet();
//    comparisonSlightChangeOfGeometry();
//   examples();
//    twoBlades();
//   newExamples();
    return 0;
}
