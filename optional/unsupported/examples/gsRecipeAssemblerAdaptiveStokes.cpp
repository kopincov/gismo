/** @file

    @brief AdaptiveSolver based on RecipeAssembler

    Author(s): A. Bressan
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsRecipeAssembler/gsRecipeAssemblerStokes.h>
#include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>

#include <gsRecipeAssemblerAdaptive/gsAdaptiveSolver.h>
#include <gsRecipeAssemblerAdaptive/gsCompositeBasisSpaceRefiners.h>
#include <gsRecipeAssemblerAdaptive/gsTensorBasisSpaceRefiners.h>

#include <gsNurbs/gsNurbsCreator.h>
#include <gsUtils/gsExportMatrix.h>

#include <gsCore/gsTransformedFuncSet.h>

#include <algorithm>

using namespace gismo;


void addAllDirichletBoundaries(const gsMultiPatch<>& mp, const gsFunction<real_t> * g, gsBoundaryConditions<real_t> & bcInfo)
{
    std::vector<patchSide> boundaries = mp.boundaries();
    for(unsigned i = 0;i<boundaries.size();++i)
        bcInfo.addCondition(boundaries[i].patch, boundaries[i].side(),  condition_type::dirichlet, g->clone().release());
}


struct OutputOption
{
    OutputOption (std::string name, index_t points)
        : fileName(name), numPoints(points)
        {}
    std::string fileName;
    index_t     numPoints;
};

class gsAdaptiveSolverWithParaviewPrint : public gsAdaptiveSolver
{
protected:
    const gsStokesPde<real_t>  *m_pde;
    const gsMultiPatch<>       &m_domain;
    OutputOption                m_opt;
public:
    gsAdaptiveSolverWithParaviewPrint(
            gsRecipeAssemblerStokes &assembler,
            gsSparseSolver<real_t>                &solverSys,
            gsSparseSolver<real_t>                &solverEli,
            gsErrorEstimator        &estimator,
            gsMarker                &marker,
            gsSpaceRefiner          &refiner,
            gsStopCriteria          &crit,
            OutputOption             options,
            gsStokesPde<real_t>     *pde=NULL)
        : gsAdaptiveSolver (assembler,solverSys, solverEli,estimator,marker,refiner,crit), m_pde(pde), m_domain(assembler.pde().domain()), m_opt(options)
    {}

    virtual void logProgress()
    {
        gsMatrix<real_t> localRatio;
        if (m_pde->solution(0)!=NULL)
            localRatio = m_estimator->getLocalErrorEstimate().row(0).array()/m_estimator->getLocalErrorEstimate().row(1).array();
        else localRatio.setConstant(1,m_estimator->getLocalErrorEstimate().cols(),1);


        size_t dofsP, dofsV;
        {
            std::vector<gsPhysicalSpace*> vec;
            vec = m_refiner->getSpaces();
            dofsP = vec[pressure]->getMapper()->getNrOfTargets();
            dofsV = vec[velocity]->getMapper()->getNrOfTargets();
        }

        gsInfo << "\niteration " << m_criteria->getCurStepNumber()<<":"
                  << "\n\t number of dofs: " << m_assembler->getSysSize()<<"="<<dofsV<<"+"<<dofsP
                  << "\n\t estimError: " << math::sqrt(m_estimator->getTotalErrorEstimate()(0,0));
        if (m_pde->solution(0))
        {
            gsInfo << "\n\t totalError: " << math::sqrt(m_estimator->getTotalErrorEstimate()(1,0))<< "\n\t ratio     : " << math::sqrt(m_estimator->getTotalErrorEstimate()(0,0)/m_estimator->getTotalErrorEstimate()(1,0))
                      << "\n\t max local ratio: " << localRatio.maxCoeff()
                      << "\n\t min local ratio: " << localRatio.minCoeff();
        }
        gsInfo << std::endl;


        gsMatrix<real_t> points;
        gsMatrix<real_t> phPoints;
        gsMatrix<real_t> vValues;
        gsMatrix<real_t> vDiv;
        gsMatrix<real_t> vLap;
        gsMatrix<real_t> pValues;
        gsMatrix<real_t> pDeriv;



        const index_t numPatches=m_domain.nPatches();
        const index_t domDim=m_domain.dim();
        const index_t tarDim=m_domain.geoDim();

        // paraview

        std::stringstream fileNameBase;
        fileNameBase<<m_opt.fileName<<m_criteria->getCurStepNumber();

        gsInfo<<"\n"<<fileNameBase.str()<<"\n";
        gsParaviewCollection collectionData(fileNameBase.str()+"_data");
        gsParaviewCollection collectionMesh(fileNameBase.str()+"_mesh");
        for (index_t p=0;p<numPatches;++p)
        {
            gsVector<real_t> lower=m_domain[p].support().col(0);
            gsVector<real_t> upper=m_domain[p].support().col(1);
            gsVector<unsigned> np = uniformSampleCount(lower,upper, m_opt.numPoints );
            points = gsPointGrid(lower,upper,np);

            evaluateSolutionAtPoints(p,points,phPoints,vValues,vDiv,vLap,pValues,pDeriv);
            //
            if ( 3 - domDim > 0 )
            {
                np.conservativeResize(3);
                np.bottomRows(3-domDim).setOnes();
            }

            if ( 3 - tarDim > 0 )
            {
                phPoints.conservativeResize(3,phPoints.cols() );
                phPoints.bottomRows(3-tarDim).setZero();
            }

            const index_t sizeDiff=3-vValues.rows();
            vValues.conservativeResize(3,vValues.cols() );
            vValues.bottomRows( sizeDiff ).setZero();
            //

            std::stringstream fileNameData;
            fileNameData<<fileNameBase.str()<<"_data"<<"_P"<<p<<".vts";
            std::stringstream fileNameMesh;
            fileNameMesh<<fileNameBase.str()<<"_mesh"<<"_P"<<p;

            std::ofstream file(fileNameData.str().c_str());
            if ( ! file.is_open() )
                gsInfo<<"Problem opening "<<fileNameData.str()<<"\n";

            file << std::fixed; // no exponents
            file << std::setprecision (5);

            file <<"<?xml version=\"1.0\"?>\n";
            file <<"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
            file <<"<StructuredGrid WholeExtent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
            file <<"<Piece Extent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
            file <<"<PointData Scalars=\"pressure\" Vectors=\"velocity\">\n";

            writeDataArray(file, "v",         vValues);
            writeDataArray(file, "p",         pValues);
            writeDataArray(file, "pGrad",     pDeriv);
            writeDataArray(file, "vDiv",      vDiv);
            writeDataArray(file, "vLap",      vLap);
            writeDataArray(file, "force",    -vLap-pDeriv);
            writeDataArray(file, "R1",  m_pde->force()->eval(points) + vLap + pDeriv);


            file <<"</PointData>\n";
            file <<"<Points>\n";

            writeDataArray(file, "",  phPoints);

            file <<"</Points>\n";
            file <<"</Piece>\n";
            file <<"</StructuredGrid>\n";
            file <<"</VTKFile>\n";

            file.close();
            collectionData.addPart(fileNameData.str());
            gsInfo<<"\tadded "<<fileNameData.str()<<"\n";

            gsMesh<real_t> msh;
            makeMesh<real_t>(*static_cast<gsStokesRefiner*>(m_refiner)->pressureBasis(p), msh,8);
            m_domain.patch(p).evaluateMesh(msh);
            gsWriteParaview(msh, fileNameMesh.str(), false);
            collectionMesh.addPart(fileNameMesh.str(),".vtp");
        }
        collectionData.save();
        collectionMesh.save();
    }

    void writeDataArray(std::ofstream &file, const std::string &name, const gsMatrix<real_t> &data)
    {
        file<<"<DataArray type=\"Float32\" ";
        if (name.length()>0)
            file<<"Name=\""<<name<<"\" format=\"ascii\" ";
        file<<" NumberOfComponents=\""<< data.rows() <<"\">\n";
        for ( index_t j=0; j<data.cols(); ++j)
            for ( index_t i=0; i<data.rows(); ++i)
                file<< data(i,j)<<" ";
        file<<"</DataArray>\n";
    }

    void evaluateSolutionAtPoints(
            index_t patch, const gsMatrix<real_t> &points,
            gsMatrix<real_t> &phPoints,
            gsMatrix<real_t> &vValues,
            gsMatrix<real_t> &vDiv,
            gsMatrix<real_t> &vLap,
            gsMatrix<real_t> &pValues,
            gsMatrix<real_t> &pDeriv
            )
    {
        const index_t tarDim=m_domain.geoDim();
        const index_t numPts=points.cols();

        gsPhysicalSpace::spacePtr vSpace = m_refiner->getSpaces()[0]->getPatchSpace(patch);
        gsPhysicalSpace::spacePtr pSpace = m_refiner->getSpaces()[1]->getPatchSpace(patch);
        const gsGeometry<>       &gParam = m_domain[patch];

        gsFuncData<>       vEval(NEED_ACTIVE|NEED_VALUE|NEED_DIV|NEED_LAPLACIAN);
        gsFuncData<>       pEval(NEED_ACTIVE|NEED_VALUE|NEED_GRAD);
        gsMapData<real_t>  gEval(
            NEED_VALUE
            |static_cast<gsTransformedFuncSet<real_t>*>(vSpace.get())->getGeoFlags(vEval.flags)
            |static_cast<gsTransformedFuncSet<real_t>*>(pSpace.get())->getGeoFlags(pEval.flags)
            );

        gsWeightMapper<real_t> *tmpMap=m_refiner->getSpaces()[0]->getMapper();
        const gsMatrix<real_t> &vCoefs=tmpMap->asMatrix()*m_solCoefs[0];
        delete tmpMap; tmpMap=m_refiner->getSpaces()[1]->getMapper();
        const gsMatrix<real_t> &pCoefs=tmpMap->asMatrix()*m_solCoefs[1];
        delete tmpMap;

        phPoints.resize (tarDim, numPts);
        vValues.setZero (tarDim, numPts);
        vDiv.setZero    (1,      numPts);
        vLap.setZero    (tarDim, numPts);

        pValues.setZero (1,      numPts);
        pDeriv.setZero  (tarDim, numPts);

        gEval.points=points;
        gParam.computeMap(gEval);
        vSpace->compute(gEval.points ,vEval);
        pSpace->compute(gEval.points ,pEval);

        phPoints=gEval.values[0];


        for (index_t p=0; p<numPts; ++p)
        {

            gsMatrix<index_t>::Column vAct  = vEval.actives.col(p);
            gsMatrix<real_t>::Column   vBVal = vEval.values[0].col(p);
            gsMatrix<real_t>::Column   vBLap = vEval.laplacians.col(p);
            const index_t vActNum=vAct.rows();
            for (index_t a=0; a<vActNum; ++a)
            {
                vValues.col(p)+=vCoefs(vAct(a),0)*vBVal.block (a*tarDim,0,tarDim,1);
                vLap.col(p)   +=vCoefs(vAct(a),0)*vBLap.block (a*tarDim,0,tarDim,1);
                vDiv(0,p)     +=vCoefs(vAct(a),0)*vEval.divs(a,p);
            }
            gsMatrix<index_t>::Column pAct  = pEval.actives.col(p);
            gsMatrix<real_t>::Column   pBVal = pEval.values[0].col(p);
            gsMatrix<real_t>::Column   pBDer = pEval.values[1].col(p);
            const index_t pActNum=pAct.rows();
            for (index_t a=0; a<pActNum; ++a)
            {
                pValues(0,p)  +=pCoefs(pAct(a),0)*pBVal(a,0);
                pDeriv.col(p) +=pCoefs(pAct(a),0)*pBDer.block(a*tarDim,0,tarDim,1);
            }
        }
    }
};

class gsStokesResidualErrorEstimator : public gsErrorEstimator , public gsRecipeAssembler
{
protected:
    const gsStokesPde<real_t>   &m_pde;
    gsMatrix<real_t>             m_vCoef;
    gsMatrix<real_t>             m_pCoef;
    std::vector<real_t>          m_cellResidual;
    std::vector<real_t>          m_cellError;
    gsMatrix<real_t>             m_perActive;
public:
    gsStokesResidualErrorEstimator(const gsStokesPde<real_t> &pde)
        : gsRecipeAssembler(pde.domain()), m_pde(pde)
    {}

    virtual void computeEstimate(
            const std::vector<gsPhysicalSpace*>  &space,
            const std::vector<gsMatrix<real_t> > &coefs
            )
    {
        setSpace(space);
        gsWeightMapper<real_t> *tmp=m_space[0]->getMapper();
        m_vCoef=tmp->asMatrix()*coefs[0];
        delete tmp;
        tmp=m_space[1]->getMapper();
        m_pCoef=tmp->asMatrix()*coefs[1];
        delete tmp;
        assemble ();
    }

    virtual void init()
    {
        gsRecipeAssembler::m_shiftsSource.resize(3);
        gsRecipeAssembler::m_shiftsSource[0]=0;
        gsRecipeAssembler::m_shiftsSource[1]=0;
        gsRecipeAssembler::m_shiftsSource[2]=0;

        m_cellResidual.clear();
        m_cellError.clear();
        m_total.setZero(2,1);
    }

    virtual void assemble ()
    {
        this->init                   ();
        for (int patch=0; patch<m_domain.size(); ++patch)
        {
            SpaceList             spList = this->getPatchSpaces(patch);
            gsIntegrationRule intRule  = this->getPatchIntegration(patch);
            gsRecipe<real_t>   recipe  = this->getPatchRecipe(patch, *intRule.subdomains);
            recipe.assemble(intRule.rule,*intRule.subdomains,m_domain[patch],spList);
        }
        this->postProcess            ();

        m_local.resize( (m_pde.solution(0)!=NULL ? 2 : 1),m_cellResidual.size());
        m_local.row(0)=gsAsConstMatrix<real_t>(m_cellResidual,1,m_cellResidual.size());
        if (m_pde.solution(0) != NULL)
            m_local.row(1)=gsAsConstMatrix<real_t>(m_cellError,1,m_cellError.size());
        m_cellResidual.clear();
        m_cellError.clear();
    }


    virtual gsRecipe<real_t> getPatchRecipe (index_t patch)
    {
        GISMO_NO_IMPLEMENTATION
    }

    virtual gsRecipe<real_t> getPatchRecipe (index_t patch, const gsDomainIterator<> &elem)
    {
        gsRecipe<real_t>           result;
        gsRecipeIngredient<real_t> ingr;

        {
            ingr.setOperator(new StokesForceResidual(m_vCoef, m_pCoef, m_pde.rhs(), elem));
            ingr.setTestSpace(1);
            ingr.setUnknownSpace(0); // ugly hack to use the pressure too
            ingr.setRule(new gsNormRule<>(*m_map, 0, true, false, m_total, m_perActive, &m_cellResidual));
        }
        result.add(ingr);

        if (m_pde.solution(0) != NULL)
        {
            ingr.setOperator(new gsDistL2(m_vCoef, m_pde.solution(0)));
            ingr.setTestSpace(1);
            ingr.setUnknownSpace(0); // ugly hack to use the pressure too
            ingr.setRule(new gsNormRule<>(*m_map, 1, true, false, m_total, m_perActive, &m_cellError));
            result.add(ingr);
        }

        return result;
    }

    virtual gsRecipe<real_t> getBoundaryRecipe (patchSide s)
    {
        gsRecipe<real_t>           result;
        return result;
    }
public:
    class StokesForceResidual : public gsDistanceOperator
    {
        const gsMatrix<real_t>   &m_vCoefs;
        const gsMatrix<real_t>   &m_pCoefs;
        const gsDomainIterator<> &m_elem;
    public:
        StokesForceResidual (const gsMatrix<real_t> &vCoefs, const gsMatrix<real_t> &pCoefs, const gsFunction<real_t> *func, const gsDomainIterator<> &elem)
            : gsDistanceOperator(vCoefs,func), m_vCoefs(vCoefs), m_pCoefs(pCoefs), m_elem(elem)
        {}

        void  pointEval (
                const gsPointFuncData<real_t>  &testSpace,
                const gsPointFuncData<real_t>  &unknownSpace,
                const gsPointMapData <real_t>  &geoEval,
                gsRecipeAccumulator<real_t>     result
                ) const
        {
            gsMatrix<real_t> laplaceU =computeLocValues2(unknownSpace.laplacian(),unknownSpace.actives(),m_vCoefs);
            gsMatrix<real_t> divU     =computeLocValues2(unknownSpace.div(),unknownSpace.actives(),m_vCoefs);
            gsMatrix<real_t> gradP    =computeLocValues2(testSpace.deriv(),testSpace.actives().array(),m_pCoefs);
            gsMatrix<real_t> valF;

            if (m_func)
                m_func->eval_into(geoEval.value(),valF);
            else
                valF.setZero(laplaceU.rows(),laplaceU.cols());

            real_t diam=(m_elem.lowerCorner()-m_elem.upperCorner()).squaredNorm();

            result += (valF+laplaceU+gradP).colwise().squaredNorm()*diam  // force residual
                        + divU.colwise().squaredNorm()                      // divergence residual
                        ;
        }

        virtual unsigned  testSpaceNeeds() const {return NEED_GRAD;}
        virtual unsigned  unknownSpaceNeeds() const {return NEED_LAPLACIAN|NEED_DIV|NEED_VALUE;}
        virtual unsigned  geometryNeeds()      const {return NEED_VALUE ;}
    };
};


gsMultiBasis<>* getStartBasis(gsMultiPatch<> &domain)
{
    gsMultiBasis<> *startBasis=new gsMultiBasis<>(domain);
    for (size_t p=0; p<startBasis->nBases();++p)
    {
        const index_t deg=startBasis->patchBases()[p]->maxDegree();
        for (index_t comp=0; comp<domain.dim(); ++comp)
        {
            const index_t dirDeg=startBasis->patchBases()[p]->degree(comp);
            startBasis->patchBases()[p]->degreeElevate(deg-dirDeg,comp);
        }
    }
    if (startBasis->degree()<2)
        startBasis->degreeElevate();
    return startBasis;
}



int main (int argc, char**argv)
{
    // read arguments

    index_t problem  = 3;
    index_t domainId = 0;
    index_t itMax    = 1;
    index_t outputPoints=1000;
    std::string outputName("adaptiveStokes_I");
    bool useUniform=false;
    bool plot = false;

    // Read input
    gsCmdLine cmd("Testing adaptive structure.");
    cmd.addInt    ("p", "problem",   "code of the example problems", problem);
    cmd.addInt    ("d", "domain" ,   "code of the example domain",   domainId);
    cmd.addInt    ("i", "iteration", "maximum number of iterations", itMax);
    cmd.addString ("o", "output",    "output filename",              outputName);
    cmd.addInt    ("P", "Points",    "numper of samples per patch",  outputPoints);
    cmd.addSwitch ("u", "uniform",   "use uniform refinement",       useUniform);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // configure console output

    gsInfo<<std::scientific;
    gsInfo<<std::setprecision(2);

    // setup problem

    gsMultiPatch<> *domains[]=
    {
        ((gsMultiPatch<>::uPtr)gsReadFile<>("planar/thbs_multipatch_01.xml")).release(),
        new gsMultiPatch<real_t>(gsNurbsCreator<real_t>::BSplineSquareGrid(1, 1, 1)),
        new gsMultiPatch<real_t>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus()),
        ((gsMultiPatch<>::uPtr)gsReadFile<>("planar/thbPuzzlePiece.xml")).release()
    };
    index_t numDom=sizeof(domains)/sizeof(gsMultiPatch<> *);
    GISMO_ASSERT(domainId<numDom,"wrong domain id");
    domains[domainId]->computeTopology();


    gsMultiPatch<>       *domain=domains[domainId];
    gsMultiBasis<>       *basis=getStartBasis(*domain);

    gsPde<>              *pde=NULL;
    gsRecipeAssembler    *assembler=NULL;
    gsSparseSolver<real_t>             *sysSolver=NULL;
    gsErrorEstimator     *estimator=NULL;
    gsMarker             *marker=NULL;
    gsStokesRefiner      *refiner=NULL;
    gsAdaptiveSolver     *solver=NULL;

    gsSparseSolver<real_t>::LU eliSolver;

    gsFunctionExpr<> f, u, p;

    switch(problem)
    {
    case 3:
    {
        gsBoundaryConditions<>  bcInfo;
        f = gsFunctionExpr<>("+2*cos(x+y)+2*sin(x-y)","-2*cos(x+y)+2*sin(x-y)",2);
        u = gsFunctionExpr<>("cos(x+y)+sin(x-y)","-1-cos(x+y)+sin(x-y)",2);
        p = gsFunctionExpr<>("0",2);

        std::vector<gsFunction<real_t>*> svec(2);
        svec[0] = const_cast<gsFunctionExpr<real_t>*>(&u);
        svec[1] = const_cast<gsFunctionExpr<real_t>*>(&p);

        pde       = new gsStokesPde<real_t>(*domain, bcInfo, &f, NULL);
        addAllDirichletBoundaries(pde->domain(),&u,pde->boundaryConditions());
        assembler = new gsRecipeAssemblerStokes(*static_cast<gsStokesPde<real_t>*>(pde));
        sysSolver = new gsSparseSolver<real_t>::QR;
        estimator = new gsErrorEstimatorPerCellExact(pde->domain(),svec);
        refiner   = new gsStokesRefiner(pde->domain(),basis->patchBases(),gsStokesRefiner::TH,0,false,3);
    }
        break;

    case 4:
    {
        f = gsFunctionExpr<>(
                    "0",
                    "0",2);
        gsFunctionExpr<> upperBoundary(
                    ".5 -abs(x-.5)",
                    "0",2);
        GISMO_ASSERT(domain->nPatches()==1,"ERROR");
        gsBoundaryConditions<>  bcInfo;
        bcInfo.addCondition(patchSide(0,1),condition_type::dirichlet,f.clone().release());
        bcInfo.addCondition(patchSide(0,2),condition_type::dirichlet,f.clone().release());
        bcInfo.addCondition(patchSide(0,3),condition_type::dirichlet,f.clone().release());
        bcInfo.addCondition(patchSide(0,4),condition_type::dirichlet,upperBoundary.clone().release());

        pde       = new gsStokesPde<real_t>( *domain, bcInfo, &f);
        assembler = new gsRecipeAssemblerStokes(*static_cast<gsStokesPde<real_t>*>(pde));
        sysSolver = new gsSparseSolver<real_t>::QR;
        estimator = new gsStokesResidualErrorEstimator(*static_cast<gsStokesPde<real_t>*>(pde));
        refiner   = new gsStokesRefiner(pde->domain(),basis->patchBases(),gsStokesRefiner::TH,0,false,3);
    }
        break;


    case 5: //Note: Buggy due to local vars (eg. force)
    {
        // construct T geometry
        gsWarn<<"domain option ignored";

        std::vector<gsGeometry<>*> patches(4);
        gsKnotVector<> knots(0,1,0,2);
        gsTensorBSplineBasis<2> geoBasis(knots,knots);
        gsMatrix<> coefs(4,2);
        coefs<<0,0,1,0,0,1,1,1;

        patches[0]=geoBasis.makeGeometry(coefs).release();
        coefs.col(1).array()+=1;
        coefs.col(0).array()+=-1;
        patches[1]=geoBasis.makeGeometry(coefs).release();
        coefs.col(0).array()+=1;
        patches[2]=geoBasis.makeGeometry(coefs).release();
        coefs.col(0).array()+=1;
        patches[3]=geoBasis.makeGeometry(coefs).release();
        coefs.col(0).array()+=1;

        domain=new gsMultiPatch<> (patches);
        //       freeAll(patches);
        domain->computeTopology();

        // update compute basis
        delete basis;
        basis=getStartBasis(*domain);

        // set up boundary conditions
        gsFunctionExpr<> zero(
                    "0",
                    "0",2);
        gsFunctionExpr<> flow(
                    "(y-1)*(y-2)",
                    "0",2);

        gsFunctionExpr<> force(
                    "1",
                    "0",2);
        gsBoundaryConditions<>  bcInfo;
        bcInfo.addCondition(patchSide(0,1),condition_type::dirichlet,zero.clone().release());
        bcInfo.addCondition(patchSide(0,2),condition_type::dirichlet,zero.clone().release());
        bcInfo.addCondition(patchSide(0,3),condition_type::dirichlet,zero.clone().release());

        bcInfo.addCondition(patchSide(1,1),condition_type::neumann,zero.clone().release());
        bcInfo.addCondition(patchSide(1,3),condition_type::dirichlet,zero.clone().release());
        bcInfo.addCondition(patchSide(1,4),condition_type::dirichlet,zero.clone().release());

        bcInfo.addCondition(patchSide(2,4),condition_type::dirichlet,zero.clone().release());

        bcInfo.addCondition(patchSide(3,2),condition_type::neumann,zero.clone().release());
        bcInfo.addCondition(patchSide(3,3),condition_type::dirichlet,zero.clone().release());
        bcInfo.addCondition(patchSide(3,4),condition_type::dirichlet,zero.clone().release());


        pde       = new gsStokesPde<real_t>( *domain, bcInfo, &force);
        assembler = new gsRecipeAssemblerStokes(*static_cast<gsStokesPde<real_t>*>(pde));
        sysSolver = new gsSparseSolver<real_t>::QR;
        estimator = new gsStokesResidualErrorEstimator(*static_cast<gsStokesPde<real_t>*>(pde));
        marker    = new gsMarkerFraction(0.2);
        refiner   = new gsStokesRefiner(pde->domain(),basis->patchBases(),gsStokesRefiner::TH,0,false,3);
        refiner->setC0(patchCorner(0,boundary::northwest));
        refiner->setC0(patchCorner(3,boundary::southwest));
        refiner->update();
    }
        break;
    case 6:
    {   // u is not H1
        gsBoundaryConditions<>  bcInfo;
        f = gsFunctionExpr<>("11","12",2);
        u = gsFunctionExpr<>("-3*(x-y)^2","-3*(x-y)^2",2);
        p = gsFunctionExpr<>("x",2);

        pde       = new gsStokesPde<real_t>(*domain, bcInfo, &f, NULL);
        addAllDirichletBoundaries(pde->domain(),&u,pde->boundaryConditions());
        assembler = new gsRecipeAssemblerStokes(*static_cast<gsStokesPde<real_t>*>(pde));
        sysSolver = new gsSparseSolver<real_t>::QR;
        estimator = new gsStokesResidualErrorEstimator(*static_cast<gsStokesPde<real_t>*>(pde));
        refiner   = new gsStokesRefiner(pde->domain(),basis->patchBases(),gsStokesRefiner::TH,0,false,3);
    }
        break;
    default:
        GISMO_ERROR("unknown choice");
    }
    gsInfo<<"\n\nSOLVING PROBLEM "<<problem<<"\n\n";


    // honor options

    refiner->setUniform(useUniform);

    // solve

    if (!marker)
        marker = new gsMarkerRelativeThreshold (0.1);
    gsStopCriteriaIterationAndTotalError criteria(estimator,itMax,1e-6);

    //if (plot)
    {
        OutputOption outOpt(outputName, outputPoints);
        solver = new gsAdaptiveSolverWithParaviewPrint(*static_cast<gsRecipeAssemblerStokes*>(assembler),*sysSolver,eliSolver,*estimator,*marker,*refiner,criteria, outOpt, static_cast<gsStokesPde<real_t> *>(pde));
        solver->adaptiveSolve();
    }

    // free memory

    delete basis;
    delete pde;
    delete assembler;
    delete estimator;
    delete marker;
    delete refiner;
    delete sysSolver;
    delete solver;

    gsMultiPatch<> **tmp=std::find(domains,domains+numDom,domain);
    if (tmp==domains+numDom)
        delete domain;
    freeAll(domains,domains+numDom);

    return 0;
}


/*
    // Poisson
    case 1:
    {
        gsBoundaryConditions<>  bcInfo;
        gsFunctionExpr<> u("x^2",2);
        gsFunctionExpr<> f("-2",2);

        pde       = new gsPoissonPde<real_t>( *domain, bcInfo, &f, &u);
        addAllDirichletBoundaries(pde->domain(),pde->solution(0),pde->boundaryConditions());
        assembler = new gsRecipeAssemblerPoisson(*static_cast<gsPoissonPde<real_t>*>(pde));
        sysSolver = new gsSparseSolver<real_t>::LU;
        estimator = new gsErrorEstimatorPerCellExact(pde->domain(),pde->solutions());
        refiner   = new gsCompositeHBasisRefiner(pde->domain(), basis->patchBases() );
    }
        break;
    case 2:
    {
        gsBoundaryConditions<>  bcInfo;
        gsFunctionExpr<> u(
                    "sin(pi*x*1)*sin(pi*y*2)+pi/10",
                    "sin(pi*x*3)*sin(pi*y*4)-pi/10",2);
        gsFunctionExpr<> f(
                    "((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",
                    "((pi*3)^2 + (pi*4)^2)*sin(pi*x*3)*sin(pi*y*4)",2);

        pde       = new gsPoissonPde<real_t>( *domain, bcInfo, &f, &u);
        addAllDirichletBoundaries(pde->domain(),pde->solution(0),pde->boundaryConditions());
        assembler = new gsRecipeAssemblerPoisson(*static_cast<gsPoissonPde<real_t>*>(pde));
        sysSolver = new gsSparseSolver<real_t>::CGDiagonal;
        estimator = new gsErrorEstimatorPerCellExact(pde->domain(),pde->solutions());
        refiner   = new gsMultiBasisRefiner(pde->domain(), *basis);
    }
        break;
*/
