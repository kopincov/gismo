/** @file gsCompositeBasis_AssemblerExamples.h

    @brief File testing the gsMappedBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsSmoothPatches/gsCompositeAssemblerUtils.h>
#include <gsSmoothPatches/gsCompositeUtils.h>

#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMappedSpline.h>

#include <gsSmoothPatches/gsCompositeBSplineBasis.h>
#include <gsSmoothPatches/gsCompositeIncrSmoothnessBasis.h>
#include <gsSmoothPatches/gsCompositeIncrSmoothnessGeom.h>

#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsRecipeAssemblerAdaptive/gsAdaptiveSolver.h>
#include <gsRecipeAssemblerAdaptive/gsCompositeBasisSpaceRefiners.h>
#include <gsRecipeAssemblerAdaptive/gsTensorBasisSpaceRefiners.h>

#include <gsSolver/gsSolverUtils.h>

#include <gsUtils/gsExportMatrix.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>



using std::flush;
using std::vector;
using std::pair;
using std::make_pair;
using std::ios;
using std::cout;
using namespace gismo;

//========================================== UTILS ========================================//

void printMatrixToFile(gsSparseMatrix<real_t>& matrix,std::string filename)
{
    std::ofstream myfile(filename.data(),std::ios::out | std::ios::trunc);
    if(myfile.is_open())
    {
        myfile << std::fixed << std::setprecision(16);
        for(index_t i = 0;i<matrix.rows();++i)
        {
            for(index_t j = 0;j<matrix.cols();++j)
                myfile << matrix(i,j) << " ";
            myfile << "\n";
        }
        myfile.close();
    }
    std::cout << "matrix " << filename << " written" <<std::endl;
}

void printVector(const std::vector<real_t>& vector)
{
    for(unsigned i = 0;i<vector.size();++i)
        gsInfo << vector[i] << " ";
}

void printSummaryOfNormsAndTimes(std::vector<std::string>& normNames,gsMatrix<real_t>& norms,gsVector<real_t>& times,gsVector<size_t>& dofs,gsVector<real_t>& conditions)
{
    std::vector<real_t> norm,ratios;
    std::cout << "norms:\n";
    for(index_t j = 0;j<norms.rows();++j)
    {
        for(index_t i = 0;i<norms.cols();++i)
            norm.push_back(norms(j,i));
        getConvergenceRatios(norm,ratios);
        cout << std::scientific;
        std::cout << normNames[j] << ": ";
        printVector(norm);
        std::cout << std::endl;
        cout.unsetf(ios::fixed | ios::scientific);
        std::cout << normNames[j] << "-ratios: ";
        printVector(ratios);
        std::cout << std::endl;
        norm.clear();
        ratios.clear();
    }
    std::cout << "times:\n" << times.transpose() << std::endl;
    std::cout << "dofs:\n" << dofs.transpose() << std::endl;
    std::cout << "conditionnumbers:\n" << conditions.transpose() << std::endl;
    cout.unsetf(ios::fixed | ios::scientific);
    cout << "\n\n";
    return;
}

//========================================== Child Classes for the adaptive Solver ========================================//

class gsCompositeIncrSmoothnessBasisRefinerUniformRefine : public gsSpaceRefiner
{
protected:
    const gsMultiPatch<>     &m_domain;
    gsCompositeIncrSmoothnessBasis<2,real_t>& m_basis;
public:
    gsCompositeIncrSmoothnessBasisRefinerUniformRefine(
            const gsMultiPatch<>             &domain,
            gsCompositeIncrSmoothnessBasis<2,real_t>& basis
            )
        : m_domain(domain),m_basis(basis)
    {

    }

    virtual std::vector<gsPhysicalSpace*> getSpaces()const
    {
        std::vector<gsPhysicalSpace*> result;
        result.push_back(new gsPhysicalSpaceScalar (m_basis.getBases(),m_domain,INVERSE_COMPOSITION,m_basis.getMapper()));
        return result;
    }

    virtual void updateSpaces (const gsMatrix<real_t>& ) // markedCells)
    {
        m_basis.uniformRefine();
    }
};

class gsCompositeBasisMapFromFileRefiner : public gsSpaceRefiner
{
protected:
    std::string                          m_pathStart;
    std::string                          m_pathEnd;
    std::string                          m_mapPathStart;
    std::string                          m_mapPathEnd;
    mutable gsMappedBasis<2,real_t>      m_curBasis;
    mutable unsigned                     m_curIter;
    gsMultiPatch<real_t>                &m_mp;
    gsBoundaryConditions<real_t>        &m_bcInfo;
    gsFunction<real_t>                  &m_exactSolution;
public:
    gsCompositeBasisMapFromFileRefiner(
            std::string                   pathStart,
            std::string                   pathEnd,
            std::string                   mapPathStart,
            std::string                   mapPathEnd,
            gsMultiPatch<real_t>         &mp,
            gsBoundaryConditions<real_t> &bcInfo,
            gsFunction<real_t>           &exactSolution
            )
        : m_pathStart(pathStart), m_pathEnd(pathEnd), m_mapPathStart(mapPathStart),
          m_mapPathEnd(mapPathEnd), m_curIter(0), m_mp(mp), m_bcInfo(bcInfo), m_exactSolution(exactSolution)
    {
        readInCurBasisAndDomain();
    }

    virtual std::vector<gsPhysicalSpace*> getSpaces()const
    {
        std::vector<gsPhysicalSpace*> result;
        result.push_back(new gsPhysicalSpaceScalar (m_curBasis.getBases(),m_mp,INVERSE_COMPOSITION,m_curBasis.getMapper()));
        return result;
    }

    virtual void updateSpaces (const gsMatrix<real_t>& ) // markedCells)
    {
        m_curIter++;
        readInCurBasisAndDomain();
    }

    gsMultiPatch<> & getDomain()
    {
        return m_mp;
    }

private:
    void readInCurBasisAndDomain() const
    {
        std::string path=m_pathStart+util::to_string(m_curIter)+m_pathEnd;
        std::string mapPath=m_pathStart+util::to_string(m_curIter)+m_mapPathEnd;
        gsMultiPatch<>::uPtr mp = gsReadFile<>(path);
        if(!mp)
        {
            gsInfo << "File could not be read." << "\n" ;
            GISMO_ERROR("File not there.");
        }
        mp->computeTopology();
        m_mp=*mp;
        m_curBasis.~gsMappedBasis<2,real_t>();
        new (&m_curBasis) gsMappedBasis<2,real_t>(m_mp,mapPath);
        m_bcInfo.clear();
        addAllDirichletBoundaries(m_mp,m_exactSolution,m_bcInfo);
    }
};

gsMatrix<real_t> calculateNorm2OfFunc(gsFunction<real_t>* func,gsMultiPatch<real_t>& domain,gsMatrix<real_t>&normMat)
{
    std::vector<gsFunction<real_t>*> functions;
    functions.push_back(func);
    std::vector<gsPhysicalSpace*> spaces;
    gsMultiBasis<real_t> basis(domain);
    gsPhysicalSpaceScalar space(basis,domain,INVERSE_COMPOSITION);
    spaces.push_back(&space);
    std::vector<gsMatrix<real_t> >coefs;
    gsMatrix<real_t> zero;
    gsWeightMapper<real_t>* mapper = space.getMapper();
    zero.setZero(mapper->getNrOfTargets(),1);
    coefs.push_back(zero);
    gsRecipeDistance distComp(domain,spaces, coefs,functions);
    distComp.setAllSeminorms(normMat);
    distComp.assemble();
    gsMatrix<real_t> seminorms = distComp.getAllDistances();
    for(index_t norm_i=1;norm_i<seminorms.rows();++norm_i)
    {
        seminorms.row(norm_i)+=seminorms.row(norm_i-1); // make norms out of seminorms
    }
    delete mapper;
    return seminorms;
}

class gsAdaptiveSolverLogData : public gsAdaptiveSolver
{
private:
    std::vector<gsMatrix<real_t> > m_absoluteNorms;
    gsMatrix<real_t> m_normOfExactSolution;
    std::vector<real_t> m_times;
    std::vector<size_t> m_dofs;
    std::vector<real_t> m_conditions;
    std::string m_description;
    bool m_conditionNeeded;


public:
    gsAdaptiveSolverLogData(gsRecipeAssembler &assembler, gsSparseSolver<real_t> &solverSys, gsSparseSolver<real_t> *solverEli,gsErrorEstimator &estimator,gsMarker *marker,
                                 gsSpaceRefiner &refiner, gsStopCriteria &crit,gsMatrix<real_t>& normExactSol,std::string description,bool getConditions=false)
        : gsAdaptiveSolver(),m_normOfExactSolution(normExactSol),m_description(description),m_conditionNeeded(getConditions)
    {
        m_assembler = &assembler;
        m_solver    = &solverSys;
        m_solverEli =  solverEli ? solverEli : &solverSys;
        m_estimator = &estimator;
        m_marker    =  marker;
        m_refiner   = &refiner;
        m_criteria  = &crit;
    }

    void assemble()
    {
        gsStopwatch time;
        gsAdaptiveSolver::assemble();
        m_times.push_back(time.stop());

        if(m_conditionNeeded)
        {
            gsSparseMatrix<real_t> sysMat = m_assembler->getSystemMatrix();
            gsVector<real_t> diagVec = sysMat.diagonal();
            for(int i = 0;i<diagVec.rows();++i)
            {
                diagVec(i)=1/sqrt(diagVec(i));
            }
            sysMat = diagVec.asDiagonal()*sysMat*diagVec.asDiagonal();
            std::stringstream ss;
            ss << m_description << "/mat" << m_criteria->getCurStepNumber();
            exportMatrixToASCII(ss.str(),sysMat);
            //printMatrixToFile(sysMat,ss.str());
            real_t cond = 0.0;//gsSolverUtils<real_t>::conditionNumber(sysMat);
            m_conditions.push_back(cond);
        }
    }

    virtual void mark()
    {
        if(m_marker)
            gsAdaptiveSolver::mark();
    }

    virtual void refine()
    {
        if(m_marker)
            gsAdaptiveSolver::refine();
        else
        {
            gsMatrix<real_t> cells;
            m_assembler->reset();
            m_refiner->updateSpaces(cells);
        }
    }

    void logProgress ()
    {
        m_absoluteNorms.push_back(m_estimator->getTotalErrorEstimate());
        m_dofs.push_back(m_assembler->getSpace()[0]->getMapper()->getNrOfTargets());
    }

    gsMatrix<real_t> postprocessNorms() const
    {
        gsMatrix<real_t> summary(m_absoluteNorms[0].rows(),m_absoluteNorms.size());
        for(size_t step = 0;step<m_absoluteNorms.size();++step)
        {
            gsMatrix<real_t> stepNorm = m_absoluteNorms[step];
            for(index_t norm_i=1;norm_i<stepNorm.rows();++norm_i)
            {
                stepNorm.row(norm_i)+=stepNorm.row(norm_i-1);
            }
            for(index_t norm_i=0;norm_i<stepNorm.rows();++norm_i)
            {
                summary(norm_i,step)=math::sqrt(stepNorm(norm_i,0))/math::sqrt(m_normOfExactSolution(norm_i,0));
            }
        }
        return summary;
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

    gsVector<real_t> getConditions() const
    {
        if(m_conditionNeeded)
        {
            gsAsConstVector<real_t> conds(m_conditions);
            return conds;
        }
        else
        {
            gsVector<real_t> conds(0);
            return conds;
        }


    }
};

//========================================== TESTS ========================================//

bool test_poissonSolvingExample(unsigned switch_var, int & plotNr)
{
    // set up problem
    gsInfo<< "\n------------ poissonSolving ------------\n\n";
    std::string path;
    gsMultiPatch<> * mp = NULL;
    gsFunctionExpr<real_t> g,dg,f;
    int iterations=3;
    switch(switch_var)
    {
    case 1:
        gsInfo << "triangle domain, degree 2\n";
        path="surfaces/multipatch_triangle2d.xml";
        g=gsFunctionExpr<real_t>("sin(x)*cos(y)",2);
        dg=gsFunctionExpr<real_t>("cos(x)*cos(y)","-sin(x)*sin(y)",2);
        f=gsFunctionExpr<real_t>("2*sin(x)*cos(y)",2);
        iterations = 3;
        break;
    case 2:
        gsInfo << "triangle domain, degree 2\n";
        path="surfaces/multipatch_triangle2d.xml";
        g=gsFunctionExpr<real_t>("1/sqrt(100000)*(-(x/12)-y)*y*(14-4*x-y)*(3-x/3-y)*(3+x/8-y)*(-(13/2)-(9*x)/4-y)",2);
        dg=gsFunctionExpr<real_t>("(y*(180*(x)^(4)+4*(x)^(3)*(518+277*y)-15*(x)^(2)*(1694+y*(-2371+569*y))-24*(-3+y)*(1092+y*(-392+y*(-1281+314*y)))-4*x*(354+y*(54419+y*(-36003+6220*y)))))/(115200.*sqrt(10))",
                                   "(36*(x)^(5)+(x)^(4)*(518+554*y)-5*(x)^(3)*(1694+y*(-4742+1707*y))-576*(-3+y)*y*(1092+y*(-593+3*y*(-33+4*y)))-2*(x)^(2)*(354+y*(108838+y*(-108009+24880*y)))-24*x*(-3276+y*(4536+y*(10353+2*y*(-4446+785*y)))))/(115200.*sqrt(10))",2);
        f=gsFunctionExpr<real_t>(" (943488+277*x^4+x^3*(11855-8175*y)-1654404*y+146906*y^2+227526*y^3-29720*y^4-2*x^2*(54419-109563*y+36489*y^2)-21*x*(2592+13042*y-16937*y^2+3995*y^3))/(57600*sqrt(10))*(-1)",2);
        iterations = 3;
        break;
//    case 3:
//        gsInfo << "four patch domain, degree 2\n";
//        path="planar/four_squares2.xml";
//        g=gsFunctionExpr<real_t>("sin(x)*cos(y)",2);
//        dg=gsFunctionExpr<real_t>("cos(x)*cos(y)","-sin(x)*sin(y)",2);
//        f=gsFunctionExpr<real_t>("2*sin(x)*cos(y)",2);
//        iterations = 3;
//        break;
//    case 4:
//        gsInfo << "four patch domain, degree 2\n";
//        path="planar/four_squares2.xml";
//        g=gsFunctionExpr<real_t>("1/sqrt(100000)*(-(x/12)-y)*y*(14-4*x-y)*(3-x/3-y)*(3+x/8-y)*(-(13/2)-(9*x)/4-y)",2);
//        dg=gsFunctionExpr<real_t>("(y*(180*(x)^(4)+4*(x)^(3)*(518+277*y)-15*(x)^(2)*(1694+y*(-2371+569*y))-24*(-3+y)*(1092+y*(-392+y*(-1281+314*y)))-4*x*(354+y*(54419+y*(-36003+6220*y)))))/(115200.*sqrt(10))",
//                                   "(36*(x)^(5)+(x)^(4)*(518+554*y)-5*(x)^(3)*(1694+y*(-4742+1707*y))-576*(-3+y)*y*(1092+y*(-593+3*y*(-33+4*y)))-2*(x)^(2)*(354+y*(108838+y*(-108009+24880*y)))-24*x*(-3276+y*(4536+y*(10353+2*y*(-4446+785*y)))))/(115200.*sqrt(10))",2);
//        f=gsFunctionExpr<real_t>(" (943488+277*x^4+x^3*(11855-8175*y)-1654404*y+146906*y^2+227526*y^3-29720*y^4-2*x^2*(54419-109563*y+36489*y^2)-21*x*(2592+13042*y-16937*y^2+3995*y^3))/(57600*sqrt(10))*(-1)",2);
//        iterations = 3;
//        break;
    case 5:
    default:
        gsInfo << "four patch domain, degree 1\n";
        path="planar/four_squares.xml";
        g=gsFunctionExpr<real_t>("x^2+y^2",2);
        dg=gsFunctionExpr<real_t>("2x",
                                   "2y",2);
        f=gsFunctionExpr<real_t>("-4",2);
        iterations = 3;
        break;
    }

    mp = (  (gsMultiPatch<>::uPtr)gsReadFile<>(path)  ).release();
    if(mp==NULL)
    {
        std::cout << "File could not be read." << "\n" ;
        return false;
    }
    gsCompositeIncrSmoothnessBasis<2,real_t> * compBasis=getCompBasisFromMultiPatch<2,real_t>(*mp);
    if(compBasis==NULL)
        return false;

    gsFunctionWithDerivatives<real_t> g_withDeriv(g,dg);
    gsBoundaryConditions<real_t> bcInfo;
    addAllDirichletBoundaries(*mp,g_withDeriv,bcInfo);
    gsPoissonPde<real_t>poissonPde(*mp,bcInfo,f,&g_withDeriv);//gets everything from the refiner

    gsCompositeIncrSmoothnessBasisRefinerUniformRefine *refiner =new gsCompositeIncrSmoothnessBasisRefinerUniformRefine(*mp,*compBasis);
    gsRecipeAssemblerPoisson *assembler =new gsRecipeAssemblerPoisson(poissonPde);
    assembler->setDirichletStrategy(dirichlet::elimination);
    assembler->setSpace(refiner->getSpaces());
    assembler->setZeroAverage(false);

    std::vector<gsFunction<real_t>*> functions;
    functions.push_back(&g_withDeriv);
    gsErrorEstimatorPerCellExact    estimator(*mp,functions);
    gsMatrix<real_t> normMat(2,2);
    normMat(0,0)=0;
    normMat(0,1)=2;
    normMat(1,0)=1;
    normMat(1,1)=2;
    estimator.setAllSeminorms(normMat);
    gsMatrix<real_t> norm2OfExactSol = calculateNorm2OfFunc(&g_withDeriv,*mp,normMat);

    gsStopCriteriaIterationOnly     stopCriteria(iterations+1);
    gsEigenCGDiagonal<real_t>               eigCG;
    eigCG.setMaxIterations(100000000);
    eigCG.setTolerance(1e-10);

    gsAdaptiveSolverLogData solver(*assembler,eigCG,NULL,estimator,NULL,*refiner,stopCriteria,norm2OfExactSol,"PoissonSolving");
    solver.adaptiveSolve();
    gsMatrix<real_t> norms = solver.postprocessNorms();
    gsVector<real_t> times = solver.getTimes();
    gsVector<size_t> dofs = solver.getDofs();
    gsVector<real_t> conds = solver.getConditions();
    std::vector<std::string> normNames;
    normNames.push_back("l2");
    normNames.push_back("h1");
    printSummaryOfNormsAndTimes(normNames,norms,times,dofs,conds);

    delete refiner;
    delete assembler;
    delete mp;
    delete compBasis;

    plotNr++;
    return true;
}



//========================================== G1-EXAMPLES ========================================//

bool test_fittingExampleG1(int switch_var, int & plotNr)
{
    std::cout << "=======================================================\n";
    std::stringstream ss;
    ss << "Fitting on a ";
    //int switch_var = 3;
    gsFunctionExpr<real_t> f; // function to fit
    std::string pathStart,pathEnd,mapPathEnd;
    int iterLimit=-1;
    switch(switch_var)
    {
    case 1:
        ss << "2-patch domain, degree 4\n";
        pathStart = "Mario/Florian_TwoPatches_deg4/multipatch_bilinear";
        pathEnd = ".xml";
        mapPathEnd = "_map.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 2:
        ss << "3-patch domain, degree 3\n";
        pathStart = "Mario/Florian_ThreePatches_deg3/multipatch_three";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 3:
        ss << "3-patch domain, degree 4\n";
        pathStart = "Mario/Florian_ThreePatches_deg4/multipatch_three";
        pathEnd = "non.xml";
        mapPathEnd = "_map.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 4:
        ss << "4-patch domain, degree 3\n";
        pathStart = "Mario/Florian_FourPatches_deg3/multipatch_four";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 5:
        ss << "4-patch domain, degree 4\n";
        pathStart = "Mario/Florian_FourPatches_deg4/multipatch_four";
        pathEnd = "non.xml";
        mapPathEnd = "_map.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 6:
        ss << "5-patch domain, degree 3\n";
        pathStart = "Mario/Florian_FivePatches_deg3/multipatch_five";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 7:
        ss << "5-patch domain, degree 4\n";
        pathStart = "Mario/Florian_FivePatches_deg4/multipatch_five";
        pathEnd = "non.xml";
        mapPathEnd = "_map.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 8:
        ss << "domain with hole, degree 4\n";
        pathStart = "Mario/Florian_HoleDomain/multipatch_hole";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        f=gsFunctionExpr<real_t>("((-5+x)*(5+x)*(-5+y)*(5+y)*(-2.25+x^2+y^2))/1000",2);
        iterLimit=1;
        break;
    case 9:
        ss << "general domain, degree 4\n";
        pathStart = "Mario/Florian_GeneralDomain/multipatch_general";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        f=gsFunctionExpr<real_t>("-((-5+x)*(5+x)*(12+5*x-3*y)*(-4+y)*(-3+2*y)*(7+2*y)*(33+4*(x)^(2)+28*y+4*(y)^(2)))/240000.",2);
        iterLimit=2;
        break;
    case 10:
        ss << "general domain new, degree 4\n";
        pathStart = "Mario/Florian_GeneralDomainNew2/multipatch_generalnew";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        f=gsFunctionExpr<real_t>("-((-5+x)*(5+x)*(12+5*x-3*y)*(-4+y)*(-3+2*y)*(7+2*y)*(33+4*(x)^(2)+28*y+4*(y)^(2)))/240000.",2);
        iterLimit=2;
        break;
    case 11:
        ss << "3-patch domain, degree 3, v2\n";
        pathStart = "Mario/Florian_ThreePatches_deg3_v2/multipatch_three";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 12:
        ss << "4-patch domain, degree 3, v2\n";
        pathStart = "Mario/Florian_FourPatches_deg3_v2/multipatch_four";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 13:
        ss << "5-patch domain, degree 3, v2\n";
        pathStart = "Mario/Florian_FivePatches_deg3_v2/multipatch_five";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 14:
        ss << "3-patch domain, degree 4, v2\n";
        pathStart = "Mario/Florian_ThreePatches_deg4_v2/multipatch_three";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 15:
        ss << "4-patch domain, degree 4, v2\n";
        pathStart = "Mario/Florian_FourPatches_deg4_v2/multipatch_four";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 16:
        ss << "5-patch domain, degree 4, v2\n";
        pathStart = "Mario/Florian_FivePatches_deg4_v2/multipatch_five";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)",2);
        iterLimit=5;
        break;
    case 17:
        ss << "general domain, degree 4, v2: Fitting\n";
        pathStart = "Mario/Florian_GeneralDomain_v2/multipatch_general";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        f=gsFunctionExpr<real_t>("((-5 + x)*(5 + x)*(-3.5 - y)*(1.5 - y)*(4 - y)*(4 + (5*x)/3. - y)*(-4 + x^2 + (3.5 + y)^2))/5000.",2);
        iterLimit=3;
        break;
    case 18:
        ss << "general domain, degree 4, v2: Fitting2\n";
        pathStart = "Mario/Florian_GeneralDomain_v2/multipatch_general";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        f=gsFunctionExpr<real_t>("2*Cos(2*x)*Sin(2*y)",2);
        iterLimit=3;
        break;
    case 19:
        ss << "circle domain, degree 4, v2\n";
        pathStart = "Mario/Florian_CircleDomain_v2/multipatch_hole";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        f=gsFunctionExpr<real_t>("((-4-y)*(5-y)*(23-(9*x)/2.-y)*(23+(9*x)/2.-y)(-2+(x)^(2)+(y)^(2)))/20000.",2);
        iterLimit=2;
        break;
    default:
        GISMO_ERROR("No such case.");
        break;
    }
    std::cout << ss.str() << std::flush;

    gsFunctionExpr<real_t> f_copy=f;
    gsMultiPatch<real_t> emptyMP;
    gsBoundaryConditions<real_t> emptyBC;
    gsCompositeBasisMapFromFileRefiner *refiner =new gsCompositeBasisMapFromFileRefiner(pathStart,pathEnd,pathStart,mapPathEnd,emptyMP,emptyBC,f_copy);
    gsRecipeAssembler *assembler =new gsRecipeAssemblerFitting(refiner->getDomain(),f_copy,refiner->getSpaces());

    std::vector<gsFunction<real_t>*> functions;
    functions.push_back(&f_copy);
    gsErrorEstimatorPerCellExact    estimator(refiner->getDomain(),functions);
    gsMatrix<real_t> normMat(1,2);
    normMat(0,0)=0;
    normMat(0,1)=2;
    estimator.setAllSeminorms(normMat);
    gsMatrix<real_t> norm2OfExactSol = calculateNorm2OfFunc(&f_copy,refiner->getDomain(),normMat);

    gsStopCriteriaIterationOnly     stopCriteria(iterLimit+1);
    gsEigenCGIdentity<real_t>               eigCG;
    eigCG.setMaxIterations(100000000);
    eigCG.setTolerance(1e-10);
    gsEigenCGIdentity<real_t>               eigCGEli;

    std::stringstream ss2;
    ss2<<"fitting"<<switch_var;
    gsAdaptiveSolverLogData solver(*assembler,eigCG,&eigCGEli,estimator,NULL,*refiner,stopCriteria,norm2OfExactSol,ss2.str(),true);
    solver.adaptiveSolve();
    if(!solver.succeed())
    {
        gsWarn << "Adaptive Solver failed. Exiting...\n";
        return false;
    }
    gsMatrix<real_t> norms = solver.postprocessNorms();
    gsVector<real_t> times = solver.getTimes();
    gsVector<size_t> dofs = solver.getDofs();
    gsVector<real_t> conds = solver.getConditions();
    std::vector<std::string> normNames;
    normNames.push_back("l2");
    printSummaryOfNormsAndTimes(normNames,norms,times,dofs,conds);

    delete refiner;
    delete assembler;

    plotNr++;
    return true;
}

bool test_poissonSolvingExampleG1(int switch_var, int & plotNr)
{
    std::cout << "=======================================================\n";
    std::stringstream ss;
    ss << "Poisson solving on a ";
    //switch_var = 1;
    std::string pathStart,pathEnd,mapPathEnd;
    gsFunctionExpr<real_t> g; // exact solution
    gsFunctionExpr<real_t> dg; // derivative of exact solution
    gsFunctionExpr<real_t> f; // right hand side
    int iterLimit=-1;
    switch(switch_var)
    {
    case 1:
        ss << "2-patch domain, degree 4\n";
        pathStart = "Mario/Florian_TwoPatches_deg4/multipatch_bilinear";
        pathEnd = ".xml";
        mapPathEnd = "_map.xml";
        g=gsFunctionExpr<real_t>("1/sqrt(100000)*(-(x/12)-y)*y*(14-4*x-y)*(3-x/3-y)*(3+x/8-y)*(-(13/2)-(9*x)/4-y)",2);
        dg=gsFunctionExpr<real_t>("(y*(180*(x)^(4)+4*(x)^(3)*(518+277*y)-15*(x)^(2)*(1694+y*(-2371+569*y))-24*(-3+y)*(1092+y*(-392+y*(-1281+314*y)))-4*x*(354+y*(54419+y*(-36003+6220*y)))))/(115200.*Sqrt(10))",
                                   "(36*(x)^(5)+(x)^(4)*(518+554*y)-5*(x)^(3)*(1694+y*(-4742+1707*y))-576*(-3+y)*y*(1092+y*(-593+3*y*(-33+4*y)))-2*(x)^(2)*(354+y*(108838+y*(-108009+24880*y)))-24*x*(-3276+y*(4536+y*(10353+2*y*(-4446+785*y)))))/(115200.*Sqrt(10))",2);
        f=gsFunctionExpr<real_t>(" (943488+277*x^4+x^3*(11855-8175*y)-1654404*y+146906*y^2+227526*y^3-29720*y^4-2*x^2*(54419-109563*y+36489*y^2)-21*x*(2592+13042*y-16937*y^2+3995*y^3))/(57600*sqrt(10))*(-1)",2);
        iterLimit=5;
        break;
    case 2:
        ss << "3-patch domain, degree 3\n";
        pathStart = "Mario/Florian_ThreePatches_deg3/multipatch_three";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        g=gsFunctionExpr<real_t>("1/500000*(-(1147/4)+77*x-y)*(-(1147/308)+(39*x/77)-y)*(-(601/155)-(16*x)/31-y)*(-(477/10)-14x-y)*((5507/1440)+(37*x/72)-y)*(31/8-x/2-y)",2);
        dg=gsFunctionExpr<real_t>("(441389437313011-5272968224115*y+8*(-298666368000*(x)^(5)-28000*(x)^(4)*(-4133821+672325*y)+480*(x)^(3)*(52101274671+5*y*(-48329457+640401050*y))+(y)^(2)*(-8076277817723+20*y*(58621313017+200*(59499360-10828163*y)*y))+12*(x)^(2)*(-565843921171+5*y*(-14511590237+20*y*(-46516301+64021780*y)))-2*x*(229800161083395+y*(2082084269689+10*y*(-1675159861883+20*y*(207397103+3706248990*y))))))/2.749824e15",
                                   "(31*(65076108692243-567148797931020*y)+4*(-(x*(5272968224115+8*x*(2082084269689+20*x*(14511590237+10*x*(144988371+18825100*x)))))+16*x*(-8076277817723+10*x*(1675159861883+680*x*(-2736253+56505975*x)))*y+12*(-15296479392949+40*x*(58621313017+10*x*(-207397103+128043560*x)))*(y)^(2)-480*(-624025745233+400*x*(-39666240+123541633*x))*(y)^(3)-8000*(-1149261923+216563260*x)*(y)^(4)+32997888000*(y)^(5)))/1.0999296e16",2);
        f=gsFunctionExpr<real_t>(" -1/109992960000000*(-322888230451989 - 231952185600*x^4 - 16017154149832*y + 46664906041472*y^2 + 1444508432256*y^3 - 467800293120*y^4 + 256*x^3*(532218639 + 97939840*y) - 32*x^2*(-1085290327639 + 5968486140*y + 42841100160*y^2) - 16*x*(594656194871 - 105863334672*y - 54886603152*y^2 + 14252015360*y^3))",2);
        iterLimit=5;
        break;
    case 3:
        ss << "3-patch domain, degree 4\n";
        pathStart = "Mario/Florian_ThreePatches_deg4/multipatch_three";
        pathEnd = "non.xml";
        mapPathEnd = "_map.xml";
        g=gsFunctionExpr<real_t>("1/500000*(-(1147/4)+77*x-y)*(-(1147/308)+(39*x/77)-y)*(-(601/155)-(16*x)/31-y)*(-(477/10)-14x-y)*((5507/1440)+(37*x/72)-y)*(31/8-x/2-y)",2);
        dg=gsFunctionExpr<real_t>("(441389437313011-5272968224115*y+8*(-298666368000*(x)^(5)-28000*(x)^(4)*(-4133821+672325*y)+480*(x)^(3)*(52101274671+5*y*(-48329457+640401050*y))+(y)^(2)*(-8076277817723+20*y*(58621313017+200*(59499360-10828163*y)*y))+12*(x)^(2)*(-565843921171+5*y*(-14511590237+20*y*(-46516301+64021780*y)))-2*x*(229800161083395+y*(2082084269689+10*y*(-1675159861883+20*y*(207397103+3706248990*y))))))/2.749824e15",
                                   "(31*(65076108692243-567148797931020*y)+4*(-(x*(5272968224115+8*x*(2082084269689+20*x*(14511590237+10*x*(144988371+18825100*x)))))+16*x*(-8076277817723+10*x*(1675159861883+680*x*(-2736253+56505975*x)))*y+12*(-15296479392949+40*x*(58621313017+10*x*(-207397103+128043560*x)))*(y)^(2)-480*(-624025745233+400*x*(-39666240+123541633*x))*(y)^(3)-8000*(-1149261923+216563260*x)*(y)^(4)+32997888000*(y)^(5)))/1.0999296e16",2);
        f=gsFunctionExpr<real_t>(" -1/109992960000000*(-322888230451989 - 231952185600*x^4 - 16017154149832*y + 46664906041472*y^2 + 1444508432256*y^3 - 467800293120*y^4 + 256*x^3*(532218639 + 97939840*y) - 32*x^2*(-1085290327639 + 5968486140*y + 42841100160*y^2) - 16*x*(594656194871 - 105863334672*y - 54886603152*y^2 + 14252015360*y^3))",2);
        iterLimit=5;
        break;
    case 4:
        ss << "4-patch domain, degree 3\n";
        pathStart = "Mario/Florian_FourPatches_deg3/multipatch_four";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        g=gsFunctionExpr<real_t>("((-40.95 - (19*x)/2. - y)*(30.88 - (37*x)/5. - y)*(4.134285714285714 - (6*x)/35. - y)*(-4.225 - x/12. - y)*(4.076190476190476 + (5*x)/42. - y)*(-4.2875 + x/8. - y)*(-50.6 + 12*x - y)*(52.93333333333333 + (37*x)/3. - y))/1.e8",2);
        dg=gsFunctionExpr<real_t>("(-183266112727112768808+93639600000000*(x)^(7)+388500000*(x)^(6)*(5638646+318315*y)-150000*(x)^(5)*(1137334053991+5*y*(38770857663+14530646650*y))+12500*(x)^(4)*(-152716930485143+25*y*(-156401552027+349513306994*y+7224284480*(y)^(2)))+1000*(x)^(3)*(70984719289425039+5*y*(-118089672499773+50*y*(-29650643938913+293423504802*y+882217595180*(y)^(2))))+75*(x)^(2)*(677492733286137563+5*y*(11600640434709637+10*y*(-981143799014747+10*y*(-12665036490069+50*y*(20291359787+6678463216*y)))))-5*y*(4654806746679084698+5*y*(-595081862571168173+5*y*(-21886279724378373+10*y*(156538493064151+160*y*(377450354536+45*y*(156065933+4364500*y))))))-20*x*(61409432176835309771+5*y*(-231054258707148997+5*y*(-273537297379319971+5*y*(686445838949453+10*y*(155942373681973+100*y*(-5142492283+817279476*y)))))))/5.292e20",
                                   "(-24*(5557475750222879756+52643682537326664905*y)+5*(3533296500000*(x)^(7)-25000*(x)^(6)*(38770857663+29061293300*y)+12500*(x)^(5)*(-156401552027+4*y*(174756653497+5418213360*y))+250*(x)^(4)*(-118089672499773+100*y*(-29650643938913+y*(440135257203+1764435190360*y)))+25*(x)^(3)*(11600640434709637+20*y*(-981143799014747+5*y*(-37995109470207+200*y*(20291359787+8348079020*y))))+30*(y)^(2)*(199526503227539871+50*y*(9930626942414458+5*y*(-17989001088697+24*y*(-141086396887+122500*y*(2701+384*y)))))-10*(x)^(2)*(-231054258707148997+5*y*(-547074594758639942+5*y*(2059337516848359+40*y*(155942373681973+25*y*(-25712461415+4903676856*y)))))-x*(4654806746679084698+5*y*(-1190163725142336346+5*y*(-65658839173135119+40*y*(156538493064151+400*y*(188725177268+9*y*(468197799+15275750*y))))))))/5.292e20",2);
        f=gsFunctionExpr<real_t>("(124581851221627307657 + 148859223125000*(x)^(6) - 4148168841948843050*y - 18010387744699264525*(y)^(2) + 220723238033909125*(y)^(3) + 226672406402041250*(y)^(4) - 732147110375000*(y)^(5) + 87342334500000*(y)^(6) - 12500*(x)^(5)*(227331388801 + 13804395780*y) - 6250*(x)^(4)*(-36474648262859 - 282855215484*y + 4857386171580*(y)^(2)) - 125*(x)^(3)*(-4035482408717607 - 458151870715570*y + 235630732858000*(y)^(2) + 37004458320000*(y)^(3)) - 125*(x)^(2)*(139889122623174041 - 2767875551846997*y - 18251735602592280*(y)^(2) + 139451974270600*(y)^(3) + 252406086414000*(y)^(4)) - 5*x*(1313780031214790431 + 251151901193160075*y - 97066558885728675*(y)^(2) - 24596791548991750*(y)^(3) + 339547972912500*(y)^(4) + 233944560600000*(y)^(5)))/2.646e19",2);
        iterLimit=5;
        break;
    case 5:
        ss << "4-patch domain, degree 4\n";
        pathStart = "Mario/Florian_FourPatches_deg4/multipatch_four";
        pathEnd = "non.xml";
        mapPathEnd = "_map.xml";
        g=gsFunctionExpr<real_t>("((-40.95 - (19*x)/2. - y)*(30.88 - (37*x)/5. - y)*(4.134285714285714 - (6*x)/35. - y)*(-4.225 - x/12. - y)*(4.076190476190476 + (5*x)/42. - y)*(-4.2875 + x/8. - y)*(-50.6 + 12*x - y)*(52.93333333333333 + (37*x)/3. - y))/1.e8",2);
        dg=gsFunctionExpr<real_t>("(-183266112727112768808+93639600000000*(x)^(7)+388500000*(x)^(6)*(5638646+318315*y)-150000*(x)^(5)*(1137334053991+5*y*(38770857663+14530646650*y))+12500*(x)^(4)*(-152716930485143+25*y*(-156401552027+349513306994*y+7224284480*(y)^(2)))+1000*(x)^(3)*(70984719289425039+5*y*(-118089672499773+50*y*(-29650643938913+293423504802*y+882217595180*(y)^(2))))+75*(x)^(2)*(677492733286137563+5*y*(11600640434709637+10*y*(-981143799014747+10*y*(-12665036490069+50*y*(20291359787+6678463216*y)))))-5*y*(4654806746679084698+5*y*(-595081862571168173+5*y*(-21886279724378373+10*y*(156538493064151+160*y*(377450354536+45*y*(156065933+4364500*y))))))-20*x*(61409432176835309771+5*y*(-231054258707148997+5*y*(-273537297379319971+5*y*(686445838949453+10*y*(155942373681973+100*y*(-5142492283+817279476*y)))))))/5.292e20",
                                   "(-24*(5557475750222879756+52643682537326664905*y)+5*(3533296500000*(x)^(7)-25000*(x)^(6)*(38770857663+29061293300*y)+12500*(x)^(5)*(-156401552027+4*y*(174756653497+5418213360*y))+250*(x)^(4)*(-118089672499773+100*y*(-29650643938913+y*(440135257203+1764435190360*y)))+25*(x)^(3)*(11600640434709637+20*y*(-981143799014747+5*y*(-37995109470207+200*y*(20291359787+8348079020*y))))+30*(y)^(2)*(199526503227539871+50*y*(9930626942414458+5*y*(-17989001088697+24*y*(-141086396887+122500*y*(2701+384*y)))))-10*(x)^(2)*(-231054258707148997+5*y*(-547074594758639942+5*y*(2059337516848359+40*y*(155942373681973+25*y*(-25712461415+4903676856*y)))))-x*(4654806746679084698+5*y*(-1190163725142336346+5*y*(-65658839173135119+40*y*(156538493064151+400*y*(188725177268+9*y*(468197799+15275750*y))))))))/5.292e20",2);
        f=gsFunctionExpr<real_t>("(124581851221627307657 + 148859223125000*(x)^(6) - 4148168841948843050*y - 18010387744699264525*(y)^(2) + 220723238033909125*(y)^(3) + 226672406402041250*(y)^(4) - 732147110375000*(y)^(5) + 87342334500000*(y)^(6) - 12500*(x)^(5)*(227331388801 + 13804395780*y) - 6250*(x)^(4)*(-36474648262859 - 282855215484*y + 4857386171580*(y)^(2)) - 125*(x)^(3)*(-4035482408717607 - 458151870715570*y + 235630732858000*(y)^(2) + 37004458320000*(y)^(3)) - 125*(x)^(2)*(139889122623174041 - 2767875551846997*y - 18251735602592280*(y)^(2) + 139451974270600*(y)^(3) + 252406086414000*(y)^(4)) - 5*x*(1313780031214790431 + 251151901193160075*y - 97066558885728675*(y)^(2) - 24596791548991750*(y)^(3) + 339547972912500*(y)^(4) + 233944560600000*(y)^(5)))/2.646e19",2);
        iterLimit=5;
        break;
    case 6:
        ss << "5-patch domain, degree 3\n";
        pathStart = "Mario/Florian_FivePatches_deg3/multipatch_five";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        g=gsFunctionExpr<real_t>("-1/100000000*(256/5-32*x/3-y)*(883/290+x/29-y)*(873/230-10*x/23-y)*(89/20+x/2-y)*(61/10+x-y)*(-617/90-13*x/9-y)*(-179/10-23*x/5-y)*(-839/190-2*x/19-y)*(-837/170+9*x/17-y)*(-108/5+9*x/2-y)",2);
        dg=gsFunctionExpr<real_t>("(5184*(-24+5*x)*(-60759448049004644133+10*x*(9224871288649613571+100*x*(37604179460108437+10*x*(-85805440864572+5*x*(-23585288349969+2*x*(-344361452181+40*x*(2225326239+88527705*x+747500*(x)^(2))))))))-135*(17023674641706872782653+20*x*(-1695509837688518024354+5*x*(-41960889498677262899+100*x*(284438575105026496+5*x*(2833672964522191+12*x*(-96263698034234+5*x*(-270006365587+80*x*(969070628+45871135*x))))))))*y-250*(3209323891457105431287+4*x*(-1063905721177541644621+20*x*(-10863701126609174361+5*x*(311839709514378744+25*x*(4409023124019543+4*x*(-39053380365891+50*x*(-47676560757+30484236416*x)))))))*(y)^(2)+500*(414137231579450182133+100*x*(-8503573356019086696+5*x*(-169890059042411955+4*x*(19655998498441548+25*x*(66279852581471+60*x*(-104658896724+22349641013*x))))))*(y)^(3)+25000*(669962652693164829+4*x*(-252981723406319757+20*x*(-3480355330614159+5*x*(-149651754464292+25*x*(1579945392551+246871418796*x)))))*(y)^(4)-50000*(36497337774982877+20*x*(-13208323166449558+15*x*(-26473933807981+40*x*(345529412146+198697687975*x))))*(y)^(5)-500000*(-1618566060310711+60*x*(-11342941558297+40*x*(-39209104117+12314778830*x)))*(y)^(6)+15000000*(-8007076526467+20*x*(-437203258412+139024822065*x))*(y)^(7)+50000000*(-475397417077+52105762340*x)*(y)^(8)-621668112000000000*(y)^(9))/5.816907e24",
                                   "(-9*(-24+5*x)*(-16077059186642225498133+4*x*(-1369336665190955685819+5*x*(48913670472575889613+20*x*(946609999700553672+5*x*(-5001444039303987+100*x*(-28130131110147+x*(156037145859+40*x*(2620771809+91742270*x))))))))-20*(145441291720013932895181+5*x*(3209323891457105431287+2*x*(-1063905721177541644621+40*x*(-3621233708869724787+5*x*(77959927378594686+5*x*(4409023124019543+10*x*(-13017793455297+100*x*(-6810937251+3810529552*x))))))))*y+300*(4238909567658291807147+x*(414137231579450182133+100*x*(-4251786678009543348+5*x*(-56630019680803985+4*x*(4913999624610387+5*x*(66279852581471+300*x*(-17443149454+3192805859*x)))))))*(y)^(2)+4000*(75855763148894346513+5*x*(669962652693164829+2*x*(-252981723406319757+40*x*(-1160118443538053+5*x*(-37412938616073+5*x*(1579945392551+205726182330*x))))))*(y)^(3)-50000*(1644683832982118973+x*(36497337774982877+20*x*(-6604161583224779+5*x*(-26473933807981+60*x*(172764706073+79479075190*x)))))*(y)^(4)-600000*(18992823414703953+x*(-1618566060310711+10*x*(-34028824674891+40*x*(-78418208234+18472168245*x))))*(y)^(5)+21000000*(64762860701939+x*(-8007076526467+20*x*(-218601629206+46341607355*x)))*(y)^(6)+80000000*(1985626097571+x*(-475397417077+26052881170*x))*(y)^(7)-2700000000*(-500083579+414445408*x)*(y)^(8)-116338140000000000*(y)^(9))/1.1633814e24",2);
        f=gsFunctionExpr<real_t>("1/290845350000000000000000*(-1379811496395894844297041 -155651117600000000*x^8+864730263236693704359840*y+280762575505560121770050*y^2-103493125039153665390000*y^3-15509526178059563535000*y^4+2700446270433556400000*y^5+295002065997385500000*y^6-3857597549580000000*y^7-131496409150000000*y^8+30000000*x^7*(128772210277+28075757565*y)+500000*x^6*(143683134330177-37677356351640*y+1674088194200*y^2)-150000*x^5*(5469773657590877-881503681940180*y-363665639267200*y^2+91909227341000*y^3)+60000*x^4*(-122502180638535341+57058996209605910*y-2114069679695175*y^2-1517941634890000*y^3+164796886687500*y^4)+1000*x^3*(32138185716817468134-8072709978165256725*y-3596653694255435100*y^2+463768931947260000*y^3+102522732097750000*y^4+5351490076650000*y^5)+50*x^2*(4308273807913771404301-2427512232578220313200*y-245340946898105477400*y^2+191051227159820224000*y^3+614771067304890000*y^4-3376333469858400000*y^5+28614863488000000*y^6)+5*x*(-61796118531102956360211+23753557112026366446720*y+6355368408723164231400*y^2-1214423672961888545000*y^3+103570695822040290000*y^4-42502401974347800000*y^5-12370109179348000000*y^6+386547891750000000*y^7))*(-1)",2);
        iterLimit=5;
        break;
    case 7:
        ss << "5-patch domain, degree 4\n";
        pathStart = "Mario/Florian_FivePatches_deg4/multipatch_five";
        pathEnd = "non.xml";
        mapPathEnd = "_map.xml";
        g=gsFunctionExpr<real_t>("-1/100000000*(256/5-32*x/3-y)*(883/290+x/29-y)*(873/230-10*x/23-y)*(89/20+x/2-y)*(61/10+x-y)*(-617/90-13*x/9-y)*(-179/10-23*x/5-y)*(-839/190-2*x/19-y)*(-837/170+9*x/17-y)*(-108/5+9*x/2-y)",2);
        dg=gsFunctionExpr<real_t>("(5184*(-24+5*x)*(-60759448049004644133+10*x*(9224871288649613571+100*x*(37604179460108437+10*x*(-85805440864572+5*x*(-23585288349969+2*x*(-344361452181+40*x*(2225326239+88527705*x+747500*(x)^(2))))))))-135*(17023674641706872782653+20*x*(-1695509837688518024354+5*x*(-41960889498677262899+100*x*(284438575105026496+5*x*(2833672964522191+12*x*(-96263698034234+5*x*(-270006365587+80*x*(969070628+45871135*x))))))))*y-250*(3209323891457105431287+4*x*(-1063905721177541644621+20*x*(-10863701126609174361+5*x*(311839709514378744+25*x*(4409023124019543+4*x*(-39053380365891+50*x*(-47676560757+30484236416*x)))))))*(y)^(2)+500*(414137231579450182133+100*x*(-8503573356019086696+5*x*(-169890059042411955+4*x*(19655998498441548+25*x*(66279852581471+60*x*(-104658896724+22349641013*x))))))*(y)^(3)+25000*(669962652693164829+4*x*(-252981723406319757+20*x*(-3480355330614159+5*x*(-149651754464292+25*x*(1579945392551+246871418796*x)))))*(y)^(4)-50000*(36497337774982877+20*x*(-13208323166449558+15*x*(-26473933807981+40*x*(345529412146+198697687975*x))))*(y)^(5)-500000*(-1618566060310711+60*x*(-11342941558297+40*x*(-39209104117+12314778830*x)))*(y)^(6)+15000000*(-8007076526467+20*x*(-437203258412+139024822065*x))*(y)^(7)+50000000*(-475397417077+52105762340*x)*(y)^(8)-621668112000000000*(y)^(9))/5.816907e24",
                                   "(-9*(-24+5*x)*(-16077059186642225498133+4*x*(-1369336665190955685819+5*x*(48913670472575889613+20*x*(946609999700553672+5*x*(-5001444039303987+100*x*(-28130131110147+x*(156037145859+40*x*(2620771809+91742270*x))))))))-20*(145441291720013932895181+5*x*(3209323891457105431287+2*x*(-1063905721177541644621+40*x*(-3621233708869724787+5*x*(77959927378594686+5*x*(4409023124019543+10*x*(-13017793455297+100*x*(-6810937251+3810529552*x))))))))*y+300*(4238909567658291807147+x*(414137231579450182133+100*x*(-4251786678009543348+5*x*(-56630019680803985+4*x*(4913999624610387+5*x*(66279852581471+300*x*(-17443149454+3192805859*x)))))))*(y)^(2)+4000*(75855763148894346513+5*x*(669962652693164829+2*x*(-252981723406319757+40*x*(-1160118443538053+5*x*(-37412938616073+5*x*(1579945392551+205726182330*x))))))*(y)^(3)-50000*(1644683832982118973+x*(36497337774982877+20*x*(-6604161583224779+5*x*(-26473933807981+60*x*(172764706073+79479075190*x)))))*(y)^(4)-600000*(18992823414703953+x*(-1618566060310711+10*x*(-34028824674891+40*x*(-78418208234+18472168245*x))))*(y)^(5)+21000000*(64762860701939+x*(-8007076526467+20*x*(-218601629206+46341607355*x)))*(y)^(6)+80000000*(1985626097571+x*(-475397417077+26052881170*x))*(y)^(7)-2700000000*(-500083579+414445408*x)*(y)^(8)-116338140000000000*(y)^(9))/1.1633814e24",2);
        f=gsFunctionExpr<real_t>("1/290845350000000000000000*(-1379811496395894844297041 -155651117600000000*x^8+864730263236693704359840*y+280762575505560121770050*y^2-103493125039153665390000*y^3-15509526178059563535000*y^4+2700446270433556400000*y^5+295002065997385500000*y^6-3857597549580000000*y^7-131496409150000000*y^8+30000000*x^7*(128772210277+28075757565*y)+500000*x^6*(143683134330177-37677356351640*y+1674088194200*y^2)-150000*x^5*(5469773657590877-881503681940180*y-363665639267200*y^2+91909227341000*y^3)+60000*x^4*(-122502180638535341+57058996209605910*y-2114069679695175*y^2-1517941634890000*y^3+164796886687500*y^4)+1000*x^3*(32138185716817468134-8072709978165256725*y-3596653694255435100*y^2+463768931947260000*y^3+102522732097750000*y^4+5351490076650000*y^5)+50*x^2*(4308273807913771404301-2427512232578220313200*y-245340946898105477400*y^2+191051227159820224000*y^3+614771067304890000*y^4-3376333469858400000*y^5+28614863488000000*y^6)+5*x*(-61796118531102956360211+23753557112026366446720*y+6355368408723164231400*y^2-1214423672961888545000*y^3+103570695822040290000*y^4-42502401974347800000*y^5-12370109179348000000*y^6+386547891750000000*y^7))*(-1)",2);
        iterLimit=5;
        break;
    case 8:
        ss << "4-patch domain, degree 3, modified rhs\n";
        pathStart = "Mario/Florian_FourPatches_deg3/multipatch_four";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        g=gsFunctionExpr<real_t>("((856 + 25*x - 210*y)*(-343 + 10*x - 80*y)*(794 + 185*x - 15*y)*(-253 + 60*x - 5*y)*(819 + 190*x + 20*y)*(-772 + 185*x + 25*y)*(507 + 10*x + 120*y)*(-1447 + 60*x + 350*y))/2.646e16",2);
        dg=gsFunctionExpr<real_t>("(-183266112727112768808+93639600000000*(x)^(7)+388500000*(x)^(6)*(5638646+318315*y)-150000*(x)^(5)*(1137334053991+5*y*(38770857663+14530646650*y))+12500*(x)^(4)*(-152716930485143+25*y*(-156401552027+349513306994*y+7224284480*(y)^(2)))+1000*(x)^(3)*(70984719289425039+5*y*(-118089672499773+50*y*(-29650643938913+293423504802*y+882217595180*(y)^(2))))+75*(x)^(2)*(677492733286137563+5*y*(11600640434709637+10*y*(-981143799014747+10*y*(-12665036490069+50*y*(20291359787+6678463216*y)))))-5*y*(4654806746679084698+5*y*(-595081862571168173+5*y*(-21886279724378373+10*y*(156538493064151+160*y*(377450354536+45*y*(156065933+4364500*y))))))-20*x*(61409432176835309771+5*y*(-231054258707148997+5*y*(-273537297379319971+5*y*(686445838949453+10*y*(155942373681973+100*y*(-5142492283+817279476*y)))))))/5.292e15",
                                   "(-24*(5557475750222879756+52643682537326664905*y)+5*(3533296500000*(x)^(7)-25000*(x)^(6)*(38770857663+29061293300*y)+12500*(x)^(5)*(-156401552027+4*y*(174756653497+5418213360*y))+250*(x)^(4)*(-118089672499773+100*y*(-29650643938913+y*(440135257203+1764435190360*y)))+25*(x)^(3)*(11600640434709637+20*y*(-981143799014747+5*y*(-37995109470207+200*y*(20291359787+8348079020*y))))+30*(y)^(2)*(199526503227539871+50*y*(9930626942414458+5*y*(-17989001088697+24*y*(-141086396887+122500*y*(2701+384*y)))))-10*(x)^(2)*(-231054258707148997+5*y*(-547074594758639942+5*y*(2059337516848359+40*y*(155942373681973+25*y*(-25712461415+4903676856*y)))))-x*(4654806746679084698+5*y*(-1190163725142336346+5*y*(-65658839173135119+40*y*(156538493064151+400*y*(188725177268+9*y*(468197799+15275750*y))))))))/5.292e15",2);
        f=gsFunctionExpr<real_t>("(124581851221627307657-4148168841948843050*y+5*(x*(-1313780031214790431+25*x*(-139889122623174041+x*(4035482408717607+50*x*(36474648262859-454662777602*x+23817475700*(x)^(2)))))-25*x*(10046076047726403+x*(-2767875551846997+10*x*(-45815187071557+60*x*(-23571267957+2300732630*x))))*y-5*(720415509787970581+5*x*(-3882662355429147+40*x*(-456293390064807+25*x*(235630732858+242869308579*x))))*(y)^(2)-25*(-1765785904271273+10*x*(-98387166195967+20*x*(697259871353+185022291600*x)))*(y)^(3)-250*(-181337925121633+150*x*(9054612611+168270724276*x))*(y)^(4)-25000*(5857176883+9357782424*x)*(y)^(5)+17468466900000*(y)^(6)))/2.646e14",2);
        iterLimit=5;
        break;
    case 9:
        ss << "4-patch domain, degree 4, modified rhs\n";
        pathStart = "Mario/Florian_FourPatches_deg4/multipatch_four";
        pathEnd = "non.xml";
        mapPathEnd = "_map.xml";
        g=gsFunctionExpr<real_t>("((856 + 25*x - 210*y)*(-343 + 10*x - 80*y)*(794 + 185*x - 15*y)*(-253 + 60*x - 5*y)*(819 + 190*x + 20*y)*(-772 + 185*x + 25*y)*(507 + 10*x + 120*y)*(-1447 + 60*x + 350*y))/2.646e16",2);
        dg=gsFunctionExpr<real_t>("(-183266112727112768808+93639600000000*(x)^(7)+388500000*(x)^(6)*(5638646+318315*y)-150000*(x)^(5)*(1137334053991+5*y*(38770857663+14530646650*y))+12500*(x)^(4)*(-152716930485143+25*y*(-156401552027+349513306994*y+7224284480*(y)^(2)))+1000*(x)^(3)*(70984719289425039+5*y*(-118089672499773+50*y*(-29650643938913+293423504802*y+882217595180*(y)^(2))))+75*(x)^(2)*(677492733286137563+5*y*(11600640434709637+10*y*(-981143799014747+10*y*(-12665036490069+50*y*(20291359787+6678463216*y)))))-5*y*(4654806746679084698+5*y*(-595081862571168173+5*y*(-21886279724378373+10*y*(156538493064151+160*y*(377450354536+45*y*(156065933+4364500*y))))))-20*x*(61409432176835309771+5*y*(-231054258707148997+5*y*(-273537297379319971+5*y*(686445838949453+10*y*(155942373681973+100*y*(-5142492283+817279476*y)))))))/5.292e15",
                                   "(-24*(5557475750222879756+52643682537326664905*y)+5*(3533296500000*(x)^(7)-25000*(x)^(6)*(38770857663+29061293300*y)+12500*(x)^(5)*(-156401552027+4*y*(174756653497+5418213360*y))+250*(x)^(4)*(-118089672499773+100*y*(-29650643938913+y*(440135257203+1764435190360*y)))+25*(x)^(3)*(11600640434709637+20*y*(-981143799014747+5*y*(-37995109470207+200*y*(20291359787+8348079020*y))))+30*(y)^(2)*(199526503227539871+50*y*(9930626942414458+5*y*(-17989001088697+24*y*(-141086396887+122500*y*(2701+384*y)))))-10*(x)^(2)*(-231054258707148997+5*y*(-547074594758639942+5*y*(2059337516848359+40*y*(155942373681973+25*y*(-25712461415+4903676856*y)))))-x*(4654806746679084698+5*y*(-1190163725142336346+5*y*(-65658839173135119+40*y*(156538493064151+400*y*(188725177268+9*y*(468197799+15275750*y))))))))/5.292e15",2);
        f=gsFunctionExpr<real_t>("(124581851221627307657-4148168841948843050*y+5*(x*(-1313780031214790431+25*x*(-139889122623174041+x*(4035482408717607+50*x*(36474648262859-454662777602*x+23817475700*(x)^(2)))))-25*x*(10046076047726403+x*(-2767875551846997+10*x*(-45815187071557+60*x*(-23571267957+2300732630*x))))*y-5*(720415509787970581+5*x*(-3882662355429147+40*x*(-456293390064807+25*x*(235630732858+242869308579*x))))*(y)^(2)-25*(-1765785904271273+10*x*(-98387166195967+20*x*(697259871353+185022291600*x)))*(y)^(3)-250*(-181337925121633+150*x*(9054612611+168270724276*x))*(y)^(4)-25000*(5857176883+9357782424*x)*(y)^(5)+17468466900000*(y)^(6)))/2.646e14",2);
        iterLimit=5;
        break;
    case 10:
        ss << "domain with hole, degree 4\n";
        pathStart = "Mario/Florian_HoleDomain/multipatch_hole";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        g=gsFunctionExpr<real_t>("((-5+x)*(5+x)*(-5+y)*(5+y)*(-2.25+x^2+y^2))/1000",2);
        dg=gsFunctionExpr<real_t>("(x*(-25+y^2)*(-109+8*x^2+4*y^2))/2000",
                                   "((-25+x^2)*y*(-109+4*x^2+8*y^2))/2000",2);
        f=gsFunctionExpr<real_t>("(-5450-4*x^4+809*y^2-4*y^4+x^2*(809-48*y^2))/2000",2);
        iterLimit=1;
        break;
    case 11:
        ss << "domain with hole, degree 4\n";
        pathStart = "Mario/Florian_HoleDomain3/multipatch_hole";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        g=gsFunctionExpr<real_t>("((-5+x)*(5+x)*(-5+y)*(5+y)*(-2.25+x^2+y^2))/1000",2);
        dg=gsFunctionExpr<real_t>("(x*(-25+y^2)*(-109+8*x^2+4*y^2))/2000",
                                   "((-25+x^2)*y*(-109+4*x^2+8*y^2))/2000",2);
        f=gsFunctionExpr<real_t>("(-5450-4*x^4+809*y^2-4*y^4+x^2*(809-48*y^2))/2000",2);
        iterLimit=1;
        break;
    case 12:
        ss << "domain with hole, degree 4\n";
        pathStart = "Mario/Florian_HoleDomain2/multipatch_hole";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        g=gsFunctionExpr<real_t>("((-5+x)*(5+x)*(-5+y)*(5+y)*(-2.25+x^2+y^2))/1000",2);
        dg=gsFunctionExpr<real_t>("(x*(-25+y^2)*(-109+8*x^2+4*y^2))/2000",
                                   "((-25+x^2)*y*(-109+4*x^2+8*y^2))/2000",2);
        f=gsFunctionExpr<real_t>("(-5450-4*x^4+809*y^2-4*y^4+x^2*(809-48*y^2))/2000",2);
        iterLimit=2;
        break;
    case 13:
        ss << "general domain, degree 4\n";
        pathStart = "Mario/Florian_GeneralDomain/multipatch_general";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        g=gsFunctionExpr<real_t>("-((-5+x)*(5+x)*(12+5*x-3*y)*(-4+y)*(-3+2*y)*(7+2*y)*(33+4*(x)^(2)+28*y+4*(y)^(2)))/240000.",2);
        dg=gsFunctionExpr<real_t>("((-4+y)*(-3+2*y)*(7+2*y)*(-100*(x)^(4)+48*(x)^(3)*(-4+y)+125*(3+2*y)*(11+2*y)-15*(x)^(2)*(-67+4*y*(7+y))+6*x*(-4+y)*(-67+4*y*(7+y))))/240000.",
                                   "-((-5+x)*(5+x)*(-24*(x)^(2)*(-4+y)*(-37-4*y+8*(y)^(2))+20*(x)^(3)*(-53+4*y*(-4+3*y))+5*x*(603+8*y*(-353+2*y*(-57+5*y*(4+y))))-6*(-4+y)*(-45+2*y*(-787+4*y*(-40+y*(29+6*y))))))/240000.",2);
        f=gsFunctionExpr<real_t>("(401289+80*(x)^(5)*(-2+3*y)+36*(x)^(4)*(7-8*(-3+y)*y)+20*(x)^(3)*(687+2*y*(-529+20*y*(1+2*y)))-3*(x)^(2)*(287+4*y*(3429+2*y*(-831-52*y+42*(y)^(2))))+3*y*(14540+y*(-101951+4*y*(449-4*y*(-439+y+(y)^(2)))))+5*x*(18416+y*(40509+4*y*(-3459+4*y*(-382+3*y*(5+y))))))/120000.",2);
        iterLimit=2;
        break;
    case 14:
        ss << "general domain new, degree 4\n";
        pathStart = "Mario/Florian_GeneralDomainNew2/multipatch_generalnew";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        g=gsFunctionExpr<real_t>("-((-5+x)*(5+x)*(12+5*x-3*y)*(-4+y)*(-3+2*y)*(7+2*y)*(33+4*(x)^(2)+28*y+4*(y)^(2)))/240000.",2);
        dg=gsFunctionExpr<real_t>("((-4+y)*(-3+2*y)*(7+2*y)*(-100*(x)^(4)+48*(x)^(3)*(-4+y)+125*(3+2*y)*(11+2*y)-15*(x)^(2)*(-67+4*y*(7+y))+6*x*(-4+y)*(-67+4*y*(7+y))))/240000.",
                                   "-((-5+x)*(5+x)*(-24*(x)^(2)*(-4+y)*(-37-4*y+8*(y)^(2))+20*(x)^(3)*(-53+4*y*(-4+3*y))+5*x*(603+8*y*(-353+2*y*(-57+5*y*(4+y))))-6*(-4+y)*(-45+2*y*(-787+4*y*(-40+y*(29+6*y))))))/240000.",2);
        f=gsFunctionExpr<real_t>("(401289+80*(x)^(5)*(-2+3*y)+36*(x)^(4)*(7-8*(-3+y)*y)+20*(x)^(3)*(687+2*y*(-529+20*y*(1+2*y)))-3*(x)^(2)*(287+4*y*(3429+2*y*(-831-52*y+42*(y)^(2))))+3*y*(14540+y*(-101951+4*y*(449-4*y*(-439+y+(y)^(2)))))+5*x*(18416+y*(40509+4*y*(-3459+4*y*(-382+3*y*(5+y))))))/120000.",2);
        iterLimit=2;
        break;
    case 15:
        ss << "3-patch domain, degree 3, v2\n";
        pathStart = "Mario/Florian_ThreePatches_deg3_v2/multipatch_three";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        g=gsFunctionExpr<real_t>("1/500000*(-(1147/4)+77*x-y)*(-(1147/308)+(39*x/77)-y)*(-(601/155)-(16*x)/31-y)*(-(477/10)-14x-y)*((5507/1440)+(37*x/72)-y)*(31/8-x/2-y)",2);
        dg=gsFunctionExpr<real_t>("(441389437313011-5272968224115*y+8*(-298666368000*(x)^(5)-28000*(x)^(4)*(-4133821+672325*y)+480*(x)^(3)*(52101274671+5*y*(-48329457+640401050*y))+(y)^(2)*(-8076277817723+20*y*(58621313017+200*(59499360-10828163*y)*y))+12*(x)^(2)*(-565843921171+5*y*(-14511590237+20*y*(-46516301+64021780*y)))-2*x*(229800161083395+y*(2082084269689+10*y*(-1675159861883+20*y*(207397103+3706248990*y))))))/2.749824e15",
                                   "(31*(65076108692243-567148797931020*y)+4*(-(x*(5272968224115+8*x*(2082084269689+20*x*(14511590237+10*x*(144988371+18825100*x)))))+16*x*(-8076277817723+10*x*(1675159861883+680*x*(-2736253+56505975*x)))*y+12*(-15296479392949+40*x*(58621313017+10*x*(-207397103+128043560*x)))*(y)^(2)-480*(-624025745233+400*x*(-39666240+123541633*x))*(y)^(3)-8000*(-1149261923+216563260*x)*(y)^(4)+32997888000*(y)^(5)))/1.0999296e16",2);
        f=gsFunctionExpr<real_t>(" -1/109992960000000*(-322888230451989 - 231952185600*x^4 - 16017154149832*y + 46664906041472*y^2 + 1444508432256*y^3 - 467800293120*y^4 + 256*x^3*(532218639 + 97939840*y) - 32*x^2*(-1085290327639 + 5968486140*y + 42841100160*y^2) - 16*x*(594656194871 - 105863334672*y - 54886603152*y^2 + 14252015360*y^3))",2);
        iterLimit=5;
        break;
    case 16:
        ss << "4-patch domain, degree 3, v2\n";
        pathStart = "Mario/Florian_FourPatches_deg3_v2/multipatch_four";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        g=gsFunctionExpr<real_t>("((-40.95 - (19*x)/2. - y)*(30.88 - (37*x)/5. - y)*(4.134285714285714 - (6*x)/35. - y)*(-4.225 - x/12. - y)*(4.076190476190476 + (5*x)/42. - y)*(-4.2875 + x/8. - y)*(-50.6 + 12*x - y)*(52.93333333333333 + (37*x)/3. - y))/1.e8",2);
        dg=gsFunctionExpr<real_t>("(-183266112727112768808+93639600000000*(x)^(7)+388500000*(x)^(6)*(5638646+318315*y)-150000*(x)^(5)*(1137334053991+5*y*(38770857663+14530646650*y))+12500*(x)^(4)*(-152716930485143+25*y*(-156401552027+349513306994*y+7224284480*(y)^(2)))+1000*(x)^(3)*(70984719289425039+5*y*(-118089672499773+50*y*(-29650643938913+293423504802*y+882217595180*(y)^(2))))+75*(x)^(2)*(677492733286137563+5*y*(11600640434709637+10*y*(-981143799014747+10*y*(-12665036490069+50*y*(20291359787+6678463216*y)))))-5*y*(4654806746679084698+5*y*(-595081862571168173+5*y*(-21886279724378373+10*y*(156538493064151+160*y*(377450354536+45*y*(156065933+4364500*y))))))-20*x*(61409432176835309771+5*y*(-231054258707148997+5*y*(-273537297379319971+5*y*(686445838949453+10*y*(155942373681973+100*y*(-5142492283+817279476*y)))))))/5.292e20",
                                   "(-24*(5557475750222879756+52643682537326664905*y)+5*(3533296500000*(x)^(7)-25000*(x)^(6)*(38770857663+29061293300*y)+12500*(x)^(5)*(-156401552027+4*y*(174756653497+5418213360*y))+250*(x)^(4)*(-118089672499773+100*y*(-29650643938913+y*(440135257203+1764435190360*y)))+25*(x)^(3)*(11600640434709637+20*y*(-981143799014747+5*y*(-37995109470207+200*y*(20291359787+8348079020*y))))+30*(y)^(2)*(199526503227539871+50*y*(9930626942414458+5*y*(-17989001088697+24*y*(-141086396887+122500*y*(2701+384*y)))))-10*(x)^(2)*(-231054258707148997+5*y*(-547074594758639942+5*y*(2059337516848359+40*y*(155942373681973+25*y*(-25712461415+4903676856*y)))))-x*(4654806746679084698+5*y*(-1190163725142336346+5*y*(-65658839173135119+40*y*(156538493064151+400*y*(188725177268+9*y*(468197799+15275750*y))))))))/5.292e20",2);
        f=gsFunctionExpr<real_t>("(124581851221627307657 + 148859223125000*(x)^(6) - 4148168841948843050*y - 18010387744699264525*(y)^(2) + 220723238033909125*(y)^(3) + 226672406402041250*(y)^(4) - 732147110375000*(y)^(5) + 87342334500000*(y)^(6) - 12500*(x)^(5)*(227331388801 + 13804395780*y) - 6250*(x)^(4)*(-36474648262859 - 282855215484*y + 4857386171580*(y)^(2)) - 125*(x)^(3)*(-4035482408717607 - 458151870715570*y + 235630732858000*(y)^(2) + 37004458320000*(y)^(3)) - 125*(x)^(2)*(139889122623174041 - 2767875551846997*y - 18251735602592280*(y)^(2) + 139451974270600*(y)^(3) + 252406086414000*(y)^(4)) - 5*x*(1313780031214790431 + 251151901193160075*y - 97066558885728675*(y)^(2) - 24596791548991750*(y)^(3) + 339547972912500*(y)^(4) + 233944560600000*(y)^(5)))/2.646e19",2);
        iterLimit=5;
        break;
    case 17:
        ss << "5-patch domain, degree 3, v2\n";
        pathStart = "Mario/Florian_FivePatches_deg3_v2/multipatch_five";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        g=gsFunctionExpr<real_t>("-1/100000000*(256/5-32*x/3-y)*(883/290+x/29-y)*(873/230-10*x/23-y)*(89/20+x/2-y)*(61/10+x-y)*(-617/90-13*x/9-y)*(-179/10-23*x/5-y)*(-839/190-2*x/19-y)*(-837/170+9*x/17-y)*(-108/5+9*x/2-y)",2);
        dg=gsFunctionExpr<real_t>("(5184*(-24+5*x)*(-60759448049004644133+10*x*(9224871288649613571+100*x*(37604179460108437+10*x*(-85805440864572+5*x*(-23585288349969+2*x*(-344361452181+40*x*(2225326239+88527705*x+747500*(x)^(2))))))))-135*(17023674641706872782653+20*x*(-1695509837688518024354+5*x*(-41960889498677262899+100*x*(284438575105026496+5*x*(2833672964522191+12*x*(-96263698034234+5*x*(-270006365587+80*x*(969070628+45871135*x))))))))*y-250*(3209323891457105431287+4*x*(-1063905721177541644621+20*x*(-10863701126609174361+5*x*(311839709514378744+25*x*(4409023124019543+4*x*(-39053380365891+50*x*(-47676560757+30484236416*x)))))))*(y)^(2)+500*(414137231579450182133+100*x*(-8503573356019086696+5*x*(-169890059042411955+4*x*(19655998498441548+25*x*(66279852581471+60*x*(-104658896724+22349641013*x))))))*(y)^(3)+25000*(669962652693164829+4*x*(-252981723406319757+20*x*(-3480355330614159+5*x*(-149651754464292+25*x*(1579945392551+246871418796*x)))))*(y)^(4)-50000*(36497337774982877+20*x*(-13208323166449558+15*x*(-26473933807981+40*x*(345529412146+198697687975*x))))*(y)^(5)-500000*(-1618566060310711+60*x*(-11342941558297+40*x*(-39209104117+12314778830*x)))*(y)^(6)+15000000*(-8007076526467+20*x*(-437203258412+139024822065*x))*(y)^(7)+50000000*(-475397417077+52105762340*x)*(y)^(8)-621668112000000000*(y)^(9))/5.816907e24",
                                   "(-9*(-24+5*x)*(-16077059186642225498133+4*x*(-1369336665190955685819+5*x*(48913670472575889613+20*x*(946609999700553672+5*x*(-5001444039303987+100*x*(-28130131110147+x*(156037145859+40*x*(2620771809+91742270*x))))))))-20*(145441291720013932895181+5*x*(3209323891457105431287+2*x*(-1063905721177541644621+40*x*(-3621233708869724787+5*x*(77959927378594686+5*x*(4409023124019543+10*x*(-13017793455297+100*x*(-6810937251+3810529552*x))))))))*y+300*(4238909567658291807147+x*(414137231579450182133+100*x*(-4251786678009543348+5*x*(-56630019680803985+4*x*(4913999624610387+5*x*(66279852581471+300*x*(-17443149454+3192805859*x)))))))*(y)^(2)+4000*(75855763148894346513+5*x*(669962652693164829+2*x*(-252981723406319757+40*x*(-1160118443538053+5*x*(-37412938616073+5*x*(1579945392551+205726182330*x))))))*(y)^(3)-50000*(1644683832982118973+x*(36497337774982877+20*x*(-6604161583224779+5*x*(-26473933807981+60*x*(172764706073+79479075190*x)))))*(y)^(4)-600000*(18992823414703953+x*(-1618566060310711+10*x*(-34028824674891+40*x*(-78418208234+18472168245*x))))*(y)^(5)+21000000*(64762860701939+x*(-8007076526467+20*x*(-218601629206+46341607355*x)))*(y)^(6)+80000000*(1985626097571+x*(-475397417077+26052881170*x))*(y)^(7)-2700000000*(-500083579+414445408*x)*(y)^(8)-116338140000000000*(y)^(9))/1.1633814e24",2);
        f=gsFunctionExpr<real_t>("1/290845350000000000000000*(-1379811496395894844297041 -155651117600000000*x^8+864730263236693704359840*y+280762575505560121770050*y^2-103493125039153665390000*y^3-15509526178059563535000*y^4+2700446270433556400000*y^5+295002065997385500000*y^6-3857597549580000000*y^7-131496409150000000*y^8+30000000*x^7*(128772210277+28075757565*y)+500000*x^6*(143683134330177-37677356351640*y+1674088194200*y^2)-150000*x^5*(5469773657590877-881503681940180*y-363665639267200*y^2+91909227341000*y^3)+60000*x^4*(-122502180638535341+57058996209605910*y-2114069679695175*y^2-1517941634890000*y^3+164796886687500*y^4)+1000*x^3*(32138185716817468134-8072709978165256725*y-3596653694255435100*y^2+463768931947260000*y^3+102522732097750000*y^4+5351490076650000*y^5)+50*x^2*(4308273807913771404301-2427512232578220313200*y-245340946898105477400*y^2+191051227159820224000*y^3+614771067304890000*y^4-3376333469858400000*y^5+28614863488000000*y^6)+5*x*(-61796118531102956360211+23753557112026366446720*y+6355368408723164231400*y^2-1214423672961888545000*y^3+103570695822040290000*y^4-42502401974347800000*y^5-12370109179348000000*y^6+386547891750000000*y^7))*(-1)",2);
        iterLimit=5;
        break;
    case 18:
        ss << "3-patch domain, degree 4, v2\n";
        pathStart = "Mario/Florian_ThreePatches_deg4_v2/multipatch_three";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        g=gsFunctionExpr<real_t>("1/500000*(-(1147/4)+77*x-y)*(-(1147/308)+(39*x/77)-y)*(-(601/155)-(16*x)/31-y)*(-(477/10)-14x-y)*((5507/1440)+(37*x/72)-y)*(31/8-x/2-y)",2);
        dg=gsFunctionExpr<real_t>("(441389437313011-5272968224115*y+8*(-298666368000*(x)^(5)-28000*(x)^(4)*(-4133821+672325*y)+480*(x)^(3)*(52101274671+5*y*(-48329457+640401050*y))+(y)^(2)*(-8076277817723+20*y*(58621313017+200*(59499360-10828163*y)*y))+12*(x)^(2)*(-565843921171+5*y*(-14511590237+20*y*(-46516301+64021780*y)))-2*x*(229800161083395+y*(2082084269689+10*y*(-1675159861883+20*y*(207397103+3706248990*y))))))/2.749824e15",
                                   "(31*(65076108692243-567148797931020*y)+4*(-(x*(5272968224115+8*x*(2082084269689+20*x*(14511590237+10*x*(144988371+18825100*x)))))+16*x*(-8076277817723+10*x*(1675159861883+680*x*(-2736253+56505975*x)))*y+12*(-15296479392949+40*x*(58621313017+10*x*(-207397103+128043560*x)))*(y)^(2)-480*(-624025745233+400*x*(-39666240+123541633*x))*(y)^(3)-8000*(-1149261923+216563260*x)*(y)^(4)+32997888000*(y)^(5)))/1.0999296e16",2);
        f=gsFunctionExpr<real_t>(" -1/109992960000000*(-322888230451989 - 231952185600*x^4 - 16017154149832*y + 46664906041472*y^2 + 1444508432256*y^3 - 467800293120*y^4 + 256*x^3*(532218639 + 97939840*y) - 32*x^2*(-1085290327639 + 5968486140*y + 42841100160*y^2) - 16*x*(594656194871 - 105863334672*y - 54886603152*y^2 + 14252015360*y^3))",2);
        iterLimit=5;
        break;
    case 19:
        ss << "4-patch domain, degree 4, v2\n";
        pathStart = "Mario/Florian_FourPatches_deg4_v2/multipatch_four";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        g=gsFunctionExpr<real_t>("((-40.95 - (19*x)/2. - y)*(30.88 - (37*x)/5. - y)*(4.134285714285714 - (6*x)/35. - y)*(-4.225 - x/12. - y)*(4.076190476190476 + (5*x)/42. - y)*(-4.2875 + x/8. - y)*(-50.6 + 12*x - y)*(52.93333333333333 + (37*x)/3. - y))/1.e8",2);
        dg=gsFunctionExpr<real_t>("(-183266112727112768808+93639600000000*(x)^(7)+388500000*(x)^(6)*(5638646+318315*y)-150000*(x)^(5)*(1137334053991+5*y*(38770857663+14530646650*y))+12500*(x)^(4)*(-152716930485143+25*y*(-156401552027+349513306994*y+7224284480*(y)^(2)))+1000*(x)^(3)*(70984719289425039+5*y*(-118089672499773+50*y*(-29650643938913+293423504802*y+882217595180*(y)^(2))))+75*(x)^(2)*(677492733286137563+5*y*(11600640434709637+10*y*(-981143799014747+10*y*(-12665036490069+50*y*(20291359787+6678463216*y)))))-5*y*(4654806746679084698+5*y*(-595081862571168173+5*y*(-21886279724378373+10*y*(156538493064151+160*y*(377450354536+45*y*(156065933+4364500*y))))))-20*x*(61409432176835309771+5*y*(-231054258707148997+5*y*(-273537297379319971+5*y*(686445838949453+10*y*(155942373681973+100*y*(-5142492283+817279476*y)))))))/5.292e20",
                                   "(-24*(5557475750222879756+52643682537326664905*y)+5*(3533296500000*(x)^(7)-25000*(x)^(6)*(38770857663+29061293300*y)+12500*(x)^(5)*(-156401552027+4*y*(174756653497+5418213360*y))+250*(x)^(4)*(-118089672499773+100*y*(-29650643938913+y*(440135257203+1764435190360*y)))+25*(x)^(3)*(11600640434709637+20*y*(-981143799014747+5*y*(-37995109470207+200*y*(20291359787+8348079020*y))))+30*(y)^(2)*(199526503227539871+50*y*(9930626942414458+5*y*(-17989001088697+24*y*(-141086396887+122500*y*(2701+384*y)))))-10*(x)^(2)*(-231054258707148997+5*y*(-547074594758639942+5*y*(2059337516848359+40*y*(155942373681973+25*y*(-25712461415+4903676856*y)))))-x*(4654806746679084698+5*y*(-1190163725142336346+5*y*(-65658839173135119+40*y*(156538493064151+400*y*(188725177268+9*y*(468197799+15275750*y))))))))/5.292e20",2);
        f=gsFunctionExpr<real_t>("(124581851221627307657 + 148859223125000*(x)^(6) - 4148168841948843050*y - 18010387744699264525*(y)^(2) + 220723238033909125*(y)^(3) + 226672406402041250*(y)^(4) - 732147110375000*(y)^(5) + 87342334500000*(y)^(6) - 12500*(x)^(5)*(227331388801 + 13804395780*y) - 6250*(x)^(4)*(-36474648262859 - 282855215484*y + 4857386171580*(y)^(2)) - 125*(x)^(3)*(-4035482408717607 - 458151870715570*y + 235630732858000*(y)^(2) + 37004458320000*(y)^(3)) - 125*(x)^(2)*(139889122623174041 - 2767875551846997*y - 18251735602592280*(y)^(2) + 139451974270600*(y)^(3) + 252406086414000*(y)^(4)) - 5*x*(1313780031214790431 + 251151901193160075*y - 97066558885728675*(y)^(2) - 24596791548991750*(y)^(3) + 339547972912500*(y)^(4) + 233944560600000*(y)^(5)))/2.646e19",2);
        iterLimit=5;
        break;
    case 20:
        ss << "5-patch domain, degree 4, v2\n";
        pathStart = "Mario/Florian_FivePatches_deg4_v2/multipatch_five";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        g=gsFunctionExpr<real_t>("-1/100000000*(256/5-32*x/3-y)*(883/290+x/29-y)*(873/230-10*x/23-y)*(89/20+x/2-y)*(61/10+x-y)*(-617/90-13*x/9-y)*(-179/10-23*x/5-y)*(-839/190-2*x/19-y)*(-837/170+9*x/17-y)*(-108/5+9*x/2-y)",2);
        dg=gsFunctionExpr<real_t>("(5184*(-24+5*x)*(-60759448049004644133+10*x*(9224871288649613571+100*x*(37604179460108437+10*x*(-85805440864572+5*x*(-23585288349969+2*x*(-344361452181+40*x*(2225326239+88527705*x+747500*(x)^(2))))))))-135*(17023674641706872782653+20*x*(-1695509837688518024354+5*x*(-41960889498677262899+100*x*(284438575105026496+5*x*(2833672964522191+12*x*(-96263698034234+5*x*(-270006365587+80*x*(969070628+45871135*x))))))))*y-250*(3209323891457105431287+4*x*(-1063905721177541644621+20*x*(-10863701126609174361+5*x*(311839709514378744+25*x*(4409023124019543+4*x*(-39053380365891+50*x*(-47676560757+30484236416*x)))))))*(y)^(2)+500*(414137231579450182133+100*x*(-8503573356019086696+5*x*(-169890059042411955+4*x*(19655998498441548+25*x*(66279852581471+60*x*(-104658896724+22349641013*x))))))*(y)^(3)+25000*(669962652693164829+4*x*(-252981723406319757+20*x*(-3480355330614159+5*x*(-149651754464292+25*x*(1579945392551+246871418796*x)))))*(y)^(4)-50000*(36497337774982877+20*x*(-13208323166449558+15*x*(-26473933807981+40*x*(345529412146+198697687975*x))))*(y)^(5)-500000*(-1618566060310711+60*x*(-11342941558297+40*x*(-39209104117+12314778830*x)))*(y)^(6)+15000000*(-8007076526467+20*x*(-437203258412+139024822065*x))*(y)^(7)+50000000*(-475397417077+52105762340*x)*(y)^(8)-621668112000000000*(y)^(9))/5.816907e24",
                                   "(-9*(-24+5*x)*(-16077059186642225498133+4*x*(-1369336665190955685819+5*x*(48913670472575889613+20*x*(946609999700553672+5*x*(-5001444039303987+100*x*(-28130131110147+x*(156037145859+40*x*(2620771809+91742270*x))))))))-20*(145441291720013932895181+5*x*(3209323891457105431287+2*x*(-1063905721177541644621+40*x*(-3621233708869724787+5*x*(77959927378594686+5*x*(4409023124019543+10*x*(-13017793455297+100*x*(-6810937251+3810529552*x))))))))*y+300*(4238909567658291807147+x*(414137231579450182133+100*x*(-4251786678009543348+5*x*(-56630019680803985+4*x*(4913999624610387+5*x*(66279852581471+300*x*(-17443149454+3192805859*x)))))))*(y)^(2)+4000*(75855763148894346513+5*x*(669962652693164829+2*x*(-252981723406319757+40*x*(-1160118443538053+5*x*(-37412938616073+5*x*(1579945392551+205726182330*x))))))*(y)^(3)-50000*(1644683832982118973+x*(36497337774982877+20*x*(-6604161583224779+5*x*(-26473933807981+60*x*(172764706073+79479075190*x)))))*(y)^(4)-600000*(18992823414703953+x*(-1618566060310711+10*x*(-34028824674891+40*x*(-78418208234+18472168245*x))))*(y)^(5)+21000000*(64762860701939+x*(-8007076526467+20*x*(-218601629206+46341607355*x)))*(y)^(6)+80000000*(1985626097571+x*(-475397417077+26052881170*x))*(y)^(7)-2700000000*(-500083579+414445408*x)*(y)^(8)-116338140000000000*(y)^(9))/1.1633814e24",2);
        f=gsFunctionExpr<real_t>("1/290845350000000000000000*(-1379811496395894844297041 -155651117600000000*x^8+864730263236693704359840*y+280762575505560121770050*y^2-103493125039153665390000*y^3-15509526178059563535000*y^4+2700446270433556400000*y^5+295002065997385500000*y^6-3857597549580000000*y^7-131496409150000000*y^8+30000000*x^7*(128772210277+28075757565*y)+500000*x^6*(143683134330177-37677356351640*y+1674088194200*y^2)-150000*x^5*(5469773657590877-881503681940180*y-363665639267200*y^2+91909227341000*y^3)+60000*x^4*(-122502180638535341+57058996209605910*y-2114069679695175*y^2-1517941634890000*y^3+164796886687500*y^4)+1000*x^3*(32138185716817468134-8072709978165256725*y-3596653694255435100*y^2+463768931947260000*y^3+102522732097750000*y^4+5351490076650000*y^5)+50*x^2*(4308273807913771404301-2427512232578220313200*y-245340946898105477400*y^2+191051227159820224000*y^3+614771067304890000*y^4-3376333469858400000*y^5+28614863488000000*y^6)+5*x*(-61796118531102956360211+23753557112026366446720*y+6355368408723164231400*y^2-1214423672961888545000*y^3+103570695822040290000*y^4-42502401974347800000*y^5-12370109179348000000*y^6+386547891750000000*y^7))*(-1)",2);
        iterLimit=5;
        break;
    case 21:
        ss << "general domain, degree 4, v2\n";
        pathStart = "Mario/Florian_GeneralDomain_v2/multipatch_general";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        g=gsFunctionExpr<real_t>("((-5+x)*(5+x)*(-3.5-y)*(1.5-y)*(4-y)*(4+(5*x)/3.-y)(-4+(x)^(2)+(3.5+y)^(2)))/5000.",2);
        dg=gsFunctionExpr<real_t>("((-5+x)*x*(5+x)*(-3.5-y)*(1.5-y)*(4-y)*(4+(5*x)/3.-y))/2500.+((-5+x)*(5+x)*(-3.5-y)*(1.5-y)*(4-y)(-4+(x)^(2)+(3.5+y)^(2)))/3000.+((-5+x)*(-3.5-y)*(1.5-y)*(4-y)*(4+(5*x)/3.-y)(-4+(x)^(2)+(3.5+y)^(2)))/5000.+((5+x)*(-3.5-y)*(1.5-y)*(4-y)*(4+(5*x)/3.-y)(-4+(x)^(2)+(3.5+y)^(2)))/5000.",
                                   "((-5+x)*(5+x)*(-3.5-y)*(1.5-y)*(4-y)*(4+(5*x)/3.-y)(3.5+y))/2500.-((-5+x)*(5+x)*(-3.5-y)*(1.5-y)*(4-y)(-4+(x)^(2)+(3.5+y)^(2)))/5000.-((-5+x)*(5+x)*(-3.5-y)*(1.5-y)*(4+(5*x)/3.-y)(-4+(x)^(2)+(3.5+y)^(2)))/5000.-((-5+x)*(5+x)*(-3.5-y)*(4-y)*(4+(5*x)/3.-y)(-4+(x)^(2)+(3.5+y)^(2)))/5000.-((-5+x)*(5+x)*(1.5-y)*(4-y)*(4+(5*x)/3.-y)(-4+(x)^(2)+(3.5+y)^(2)))/5000.",2);
        f=gsFunctionExpr<real_t>("(401289+43620*y-305853*(y)^(2)+5388*(y)^(3)+21072*(y)^(4)-48*(y)^(5)-48*(y)^(6)+80*(x)^(5)*(-2+3*y)-36*(x)^(4)*(-7-24*y+8*(y)^(2))+20*(x)^(3)*(687-1058*y+40*(y)^(2)+80*(y)^(3))-3*(x)^(2)*(287+13716*y-6648*(y)^(2)-416*(y)^(3)+336*(y)^(4))+5*x*(18416+40509*y-13836*(y)^(2)-6112*(y)^(3)+240*(y)^(4)+48*(y)^(5)))/120000.",2);
        iterLimit=3;
        break;
    case 22:
        ss << "circle domain, degree 4, v2\n";
        pathStart = "Mario/Florian_CircleDomain_v2/multipatch_hole";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        g=gsFunctionExpr<real_t>("((-4-y)*(5-y)*(23-(9*x)/2.-y)*(23+(9*x)/2.-y)(-2+(x)^(2)+(y)^(2)))/20000.",2);
        dg=gsFunctionExpr<real_t>("(x*(-4-y)*(5-y)*(23-(9*x)/2.-y)*(23+(9*x)/2.-y))/10000.+(9*(-4-y)*(5-y)*(23-(9*x)/2.-y)*(-2+(x)^(2)+(y)^(2)))/40000.-(9*(-4-y)*(5-y)*(23+(9*x)/2.-y)(-2+(x)^(2)+(y)^(2)))/40000.",
                                   "((-4-y)*(5-y)*(23-(9*x)/2.-y)*(23+(9*x)/2.-y)*y)/10000.-((-4-y)*(5-y)*(23-(9*x)/2.-y)*(-2+(x)^(2)+(y)^(2)))/20000.-((-4-y)*(5-y)*(23+(9*x)/2.-y)(-2+(x)^(2)+(y)^(2)))/20000.-((-4-y)*(23-(9*x)/2.-y)*(23+(9*x)/2.-y)(-2+(x)^(2)+(y)^(2)))/20000.-((5-y)*(23-(9*x)/2.-y)*(23+(9*x)/2.-y)(-2+(x)^(2)+(y)^(2)))/20000.",2);
        f=gsFunctionExpr<real_t>("(92320+81*(x)^(4)-7222*y-17274*(y)^(2)+1987*(y)^(3)+17*(y)^(4)-(x)^(2)*(13722+165*y-948*(y)^(2)))/40000.",2);
        iterLimit=2;
        break;
    default:
        GISMO_ERROR("No such case.");
        break;
    }
    std::cout << ss.str() << std::flush;

    gsFunctionWithDerivatives<real_t> g_withDeriv(g,dg);
    gsMultiPatch<real_t> emptyMP;
    gsBoundaryConditions<real_t> emptyBC;
    gsPoissonPde<real_t>poissonPde(emptyMP,emptyBC,f,&g_withDeriv);//gets everything from the refiner

    gsCompositeBasisMapFromFileRefiner *refiner =new gsCompositeBasisMapFromFileRefiner(pathStart,pathEnd,pathStart,mapPathEnd,poissonPde.domain(),poissonPde.boundaryConditions(),g_withDeriv);
    gsRecipeAssemblerPoisson *assembler =new gsRecipeAssemblerPoisson(poissonPde);
    assembler->setDirichletStrategy(dirichlet::elimination);
    assembler->setSpace(refiner->getSpaces());
    assembler->setZeroAverage(false);

    std::vector<gsFunction<real_t>*> functions;
    functions.push_back(&g_withDeriv);
    gsErrorEstimatorPerCellExact    estimator(refiner->getDomain(),functions);
    gsMatrix<real_t> normMat(2,2);
    normMat(0,0)=0;
    normMat(0,1)=2;
    normMat(1,0)=1;
    normMat(1,1)=2;
    estimator.setAllSeminorms(normMat);
    gsMatrix<real_t> normOfExactSol = calculateNorm2OfFunc(&g_withDeriv,refiner->getDomain(),normMat);

    gsStopCriteriaIterationOnly     stopCriteria(iterLimit+1);
    gsEigenCGDiagonal<real_t>               eigCG;
    eigCG.setMaxIterations(100000000);
    eigCG.setTolerance(1e-10);

    std::stringstream ss2;
    ss2<<"poisson"<<switch_var;
    gsAdaptiveSolverLogData solver(*assembler,eigCG,NULL,estimator,NULL,*refiner,stopCriteria,normOfExactSol,ss2.str(),true);
    solver.adaptiveSolve();
    if(!solver.succeed())
    {
        gsWarn << "Adaptive Solver failed. Exiting...\n";
        return false;
    }
    gsMatrix<real_t> norms = solver.postprocessNorms();
    gsVector<real_t> times = solver.getTimes();
    gsVector<size_t> dofs = solver.getDofs();
    gsVector<real_t> conds = solver.getConditions();
    std::vector<std::string> normNames;
    normNames.push_back("l2");
    normNames.push_back("h1");
    printSummaryOfNormsAndTimes(normNames,norms,times,dofs,conds);

    delete refiner;
    delete assembler;

    plotNr++;
    return true;
}

bool test_biharmonicExampleG1(int switch_var,bool plot, int & plotNr)
{
    std::cout << "=======================================================\n";
    std::stringstream ss;
    ss << "Solving the biharmonic equation on a ";
    //int switch_var = 7;
    std::string pathStart,pathEnd,mapPathEnd;
    gsFunctionExpr<real_t> g; // exact solution
    gsFunctionExpr<real_t> dg; // derivative of exact solution
    gsFunctionExpr<real_t> ddg; // 2nd derivative of exact solution
    gsFunctionExpr<real_t> f; // right hand side
    int startIter=0,iterLimit=-1;
    switch(switch_var)
    {
    case 1:
        ss << "2-patch domain, degree 4\n";
        pathStart = "Mario/Florian_TwoPatches_deg4/multipatch_bilinear";
        pathEnd = ".xml";
        mapPathEnd = "_map.xml";
        startIter=1;
        g=gsFunctionExpr<real_t>("((14 - 4*x - y)^(2)*(-6.5 - (9*x)/4. - y)^2*(3 - x/3. - y)^(2)*(-x/12. - y)^(2)*(3 + x/8. - y)^(2)*(y)^(2))/100000",2); // 1
        dg=gsFunctionExpr<real_t>("((24+x-8*y)*(y)^(2)*(-14+4*x+y)*(-9+x+3*y)*(26+9*x+4*y)*(x+12*y)*(180*(x)^(4)+4*(x)^(3)*(518+277*y)-15*(x)^(2)*(1694+y*(-2371+569*y))-24*(-3+y)*(1092+y*(-392+y*(-1281+314*y)))-4*x*(354+y*(54419+y*(-36003+6220*y)))))/6.63552e10",
                                   "((24+x-8*y)*y*(-14+4*x+y)*(-9+x+3*y)*(26+9*x+4*y)*(x+12*y)*(36*(x)^(5)+(x)^(4)*(518+554*y)-5*(x)^(3)*(1694+y*(-4742+1707*y))-576*(-3+y)*y*(1092+y*(-593+3*y*(-33+4*y)))-2*(x)^(2)*(354+y*(108838+y*(-108009+24880*y)))-24*x*(-3276+y*(4536+y*(10353+2*y*(-4446+785*y))))))/6.63552e10",2);
        ddg=gsFunctionExpr<real_t>("((y)^(2)*(2*(-14+4*x+y)^(2)*(-9+x+3*y)^(2)*(26+9*x+4*y)^(2)*(288+3*(x)^(2)+8*(48-11*y)*y+12*x*(6+y))+8*(24+x-8*y)*(-14+4*x+y)*(-9+x+3*y)*(26+9*x+4*y)*(x+12*y)*(x+2*(6+y))*(-166+108*(x)^(2)+y*(-321+79*y)+x*(-692+266*y))+(24+x-8*y)^(2)*(x+12*y)^(2)*(24*x*(43843+x*(26941+30*x*(-346+27*x)))+24*x*(14849+x*(-28787+3990*x))*y+(-54647+6*x*(-74563+23377*x))*(y)^(2)+6*(-15423+10939*x)*(y)^(3)+9433*(y)^(4)+4*(-559859+386703*y))))/6.63552e10",
                                    "((-14+4*x+y)^(2)*(-9+x+3*y)^(2)*(26+9*x+4*y)^(2)*((x)^(2)*(24+x)^(2)+24*x*(24+x)*(72+x)*y-96*(-5184+x*(144+11*x))*(y)^(2)-7680*(72+x)*(y)^(3)+138240*(y)^(4))+4*(24+x-8*y)*y*(-14+4*x+y)*(-9+x+3*y)*(26+9*x+4*y)*(x+12*y)*((x)^(2)-288*(-2+y)*y+8*x*(3+y))*(133*(x)^(2)+x*(-321+158*y)+6*(-137+6*(-7+y)*y))+(24+x-8*y)^(2)*(y)^(2)*(x+12*y)^(2)*(23377*(x)^(4)+2*(x)^(3)*(-74563+32817*y)+(x)^(2)*(-54647+6*y*(-46269+9433*y))+12*x*(90597+2*y*(-6621+y*(-6903+790*y)))+36*(-4163+6*y*(3969+y*(-107+10*(-14+y)*y)))))/6.63552e10",
                                    "(y*(12960*(x)^(9)+324*(x)^(8)*(1036+831*y)-16*(x)^(7)*(170758+y*(-855399+128111*y))-7*(x)^(6)*(8825896+y*(370404+y*(-17609364+6179525*y)))+6*(x)^(5)*(76668340+y*(-476836806+y*(394147922+3*y*(-30508435+659673*y))))+20*(x)^(4)*(23361996+y*(705759870+y*(-1975975898+y*(1593588240+y*(-511347813+58209998*y)))))+13824*(-3+y)^(2)*y*(1788696+y*(-1852032+y*(-3233342+y*(2626554+y*(-195741-87792*y+6908*(y)^(2))))))+16*(x)^(3)*(-332847324+y*(1469182986+y*(5680225514+y*(-11547130575+y*(7347287523+14*y*(-141646149+14071088*y))))))+3456*x*(-3+y)*(-1192464+y*(2465736+y*(83636424+y*(-103877202+y*(30357811+2*y*(2206399+3*y*(-465327+41210*y)))))))+144*(x)^(2)*(-2319408+y*(-1031880276+y*(2291524536+y*(-573019985+y*(-1143280404+y*(840737765+144*y*(-1463239+126335*y)))))))))/6.63552e10",2);
        f=gsFunctionExpr<real_t>("(-326013*(x)^(8) - 6*(x)^(7)*(-4626117 + 5820533*y) + (x)^(6)*(581659435 - 453347142*y + 8979117*(y)^(2)) + 12*(x)^(5)*(-1003433267 + 3185880066*y - 2463115080*(y)^(2) + 539471541*(y)^(3)) + 3*(x)^(4)*(11743792728 - 97145412660*y + 152834063170*(y)^(2) - 81072512300*(y)^(3) + 13830173945*(y)^(4)) + 6*(x)^(3)*(27654041072 - 18094827680*y - 176867619620*(y)^(2) + 244268428620*(y)^(3) - 105387539885*(y)^(4) + 14510328525*(y)^(5)) + 3*(x)^(2)*(-217317993696 + 1389553239456*y - 1546898048964*(y)^(2) - 78453020700*(y)^(3) + 683577414495*(y)^(4) - 266711936886*(y)^(5) + 28986769847*(y)^(6)) + 12*x*(-47375834688 - 136733927904*y + 794323156956*(y)^(2) - 742422243980*(y)^(3) + 108560774575*(y)^(4) + 101767828104*(y)^(5) - 37311029779*(y)^(6) + 3262055060*(y)^(7)) + 12*(223057545984 - 1302675768576*y + 1333416480804*(y)^(2) + 74258561244*(y)^(3) - 443818288763*(y)^(4) + 91940509914*(y)^(5) + 22882187233*(y)^(6) - 7593579540*(y)^(7) + 533926888*(y)^(8)))/1.65888e10",2); // 1
        iterLimit=5;
        break;
    case 2:
        ss << "3-patch domain, degree 3\n";
        pathStart = "Mario/Florian_ThreePatches_deg3/multipatch_three";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        startIter=1;
        g=gsFunctionExpr<real_t>("((-47.7 - 14*x - y)^(2)*(-3.8774193548387097 - (16*x)/31. - y)^2*(3.875 - x/2. - y)^(2)*(-3.7240259740259742 + (39*x)/77. - y)^2*(3.8243055555555556 + (37*x)/72. - y)^2*(-286.75 + 77*x - y)^(2))/2.5e11",2); // 2
        dg=gsFunctionExpr<real_t>("((5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-441389437313011+5272968224115*y+8*(298666368000*(x)^(5)+28000*(x)^(4)*(-4133821+672325*y)-480*(x)^(3)*(52101274671+5*y*(-48329457+640401050*y))+12*(x)^(2)*(565843921171+5*y*(14511590237+20*(46516301-64021780*y)*y))+(y)^(2)*(8076277817723+20*y*(-58621313017+200*y*(-59499360+10828163*y)))+2*x*(229800161083395+y*(2082084269689+10*y*(-1675159861883+20*y*(207397103+3706248990*y)))))))/1.5123064061952e31",
                                   "((5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(31*(-65076108692243+567148797931020*y)+4*(x*(5272968224115+8*x*(2082084269689+20*x*(14511590237+10*x*(144988371+18825100*x))))-16*x*(-8076277817723+10*x*(1675159861883+680*x*(-2736253+56505975*x)))*y+12*(15296479392949-40*x*(58621313017+10*x*(-207397103+128043560*x)))*(y)^(2)+480*(-624025745233+400*x*(-39666240+123541633*x))*(y)^(3)+8000*(-1149261923+216563260*x)*(y)^(4)-32997888000*(y)^(5))))/6.0492256247808e31",2);
        ddg=gsFunctionExpr<real_t>("(2*(5507+740*x-1440*y)^(2)*(1147-308*x+4*y)^(2)*(1147-156*x+308*y)^(2)*(139686329327-200681825360*y+25*(16*x*(-362657382+x*(-101620989+28000*x*(355+84*x)))+8*x*(-373082571+20*x*(4913781+1257200*x))*y+(-653061619+240*x*(5025693+1529365*x))*(y)^(2)+110*(4996423+2317980*x)*(y)^(3)+48905650*(y)^(4)))+8*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-489329699+26666640*(x)^(2)+8*y*(15792077+4325450*y)-8*x*(8077099+8740640*y))*(-330574+75085*y+25*(672*(x)^(2)+931*(y)^(2)+4*x*(355+449*y)))+(-31+4*x+8*y)^(2)*(477+140*x+10*y)^(2)*(601+80*x+155*y)^(2)*(-3441*(-35572636019847+73096675692352*y)+16*(74073925926000*(x)^(4)-44444400*(x)^(3)*(8077099+8740640*y)+4*(y)^(2)*(-172228484745029+20*y*(13798425242429+959704939205*y))+6*(x)^(2)*(-478459593600089+80*y*(3519648832631+1435589920070*y))+3*x*(3988696521770833-8*y*(-206187448112451+10*y*(10414257900957+1896520058240*y))))))/3.780766015488e30",
                                    "(-32*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-184879+179600*(x)^(2)+60*y*(19717+620*y)+40*x*(15017+9310*y))*(-9535011+17481280*(x)^(2)+4*y*(31783721+166320*y)-4*x*(15792077+8650900*y))+1573679923200*(-31+4*x+8*y)^(2)*(477+140*x+10*y)^(2)*(601+80*x+155*y)^(2)*(((11119+113140*x-221760*y)*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y))/(2.45887488e10)+2*(-3.7240259740259742+(39*x)/77.-y)^2*(3.8243055555555556+(37*x)/72.-y)^2+2*(286.75-77*x+y)^(2)*((-350064735599+200*x*(15438523+96003266*x))/(1.22943744e10)-((11119+113140*x)*y)/18480.+6*(y)^(2)))+(5507+740*x-1440*y)^(2)*(1147-308*x+4*y)^(2)*(1147-156*x+308*y)^(2)*(-10479303116099-1317338496540*y+40*(4*x*(-41626408888+5*x*(-653061619+40*x*(10051386+1529365*x)))+330*x*(130950371+20*x*(4996423+772660*x))*y+15*(3452990809+40*x*(58793797+9781130*x))*(y)^(2)+62000*(59151+18620*x)*(y)^(3)+57660000*(y)^(4))))/6.0492256247808e31",
                                    "(263879111196180480000000*(x)^(10)-200704000000*(x)^(9)*(-8649747410767+243871656618475*y)-1102267*(-249908537057869493753722+7303900366010386651975527*y)-38707200000*(x)^(8)*(290585708203906+3*y*(-197029885967433+80994558247885*y))+8192000*(x)^(7)*(-9481688358976200867+5*y*(78477855207928463297+50*y*(-12621679664055075+221243488766886308*y)))+2867200*(x)^(6)*(-280977579189855473287+5*y*(-22493549631768585111+5*y*(-2747278652088440511+20*y*(-91855653718880188+46208541740778675*y))))-92160*(x)^(5)*(75483626813060289160351+25*y*(-20046466182105866760727+100*y*(11594160804861041637+20*y*(-887558443306536565+4*y*(1720896663954649+71090383325622549*y)))))-7680*(x)^(4)*(-4964239415059682722723998+5*y*(2296897469493749744928169+25*y*(4891483627964772371247+40*y*(5266351819578766236+5*y*(-13129448820140085149+80*y*(-26146153193635401+10391360459106683*y))))))+1280*(x)^(3)*(319898962438204378868630125+y*(-3820925583992735261392046391+10*y*(-4783100735284830070720533+20*y*(2105413765947766789130671+40*y*(-281075782689312175923+5*y*(-527029110724154050557+8*y*(1102027571944353971+2738805988130334580*y)))))))+(y)^(2)*(617841001190024474584165568361-80*y*(-21074739527441629699914866180+y*(2607549751484824916038727667+16*y*(91673462096616111761339829+125*y*(-107496658047501834366731+80*y*(-225835710267402201780+y*(31492159069935608973+80*y*(12239879471796649+40941370928304*y))))))))+288*(x)^(2)*(-920093625250723826423869547+10*y*(892691419041710198843985595+y*(-47417212378914652851450213+20*y*(-6547575133692972562670908+5*y*(131944367478582199929979+32*y*(2394158422480936839957+5*y*(-79873367397145988827+200*y*(-56248512753087052+12046162278210363*y))))))))+4*x*(-1197623698331623912863154084953+y*(12321028396612468035229975389749+40*y*(16186340135808633803985384165+4*y*(-13097715815752721726226123093+20*y*(-24778019926500117047201121+y*(43068787960711552962518799+160*y*(8538946273367625114111+10*y*(-695886202065453168197+5*y*(-4066084107601811649+10710407262566600*y))))))))))/7.561532030976e30",2);
        f=gsFunctionExpr<real_t>("(775*(10070636445024446987294496713 + 1580805041435296880732831616*y) + 32*(839212170045134580480000*(x)^(8) + 3072000*(x)^(7)*(-186153527560215013 + 23630822384505615*y) - 12800*(x)^(6)*(7397845103057974754479 + 600*y*(-1978168936115700 + 521427067157777179*y)) - 23040*(x)^(5)*(-2591139538241579069079 + 5*y*(-19502249110340706949 + 40*y*(-532005332078551427 + 231162700790668990*y))) + 240*(x)^(4)*(20629697508869638693714899 + 80*y*(-2832250234210447845948 + 5*y*(361036991177109628031 + 80*y*(851520842693099156 + 1397012256516649535*y)))) + 480*(x)^(3)*(-5011292996454793182398419 + y*(172125474399395900332421 + 80*y*(6788570295304368834777 + 10*y*(-121146247699783170709 + 40*y*(-983617023935690419 + 258889960348724901*y))))) + (y)^(2)*(-113362049203585367125403749833 + 80*y*(-133919062100907890758222638 + y*(119712510000704219793039297 + 16*y*(506422009899324159806439 + 5*y*(-20050351692839232072401 + 120*y*(-4976219400498134306 + 701059526547943945*y)))))) + 6*(x)^(2)*(-13034991506843954692628435037 + 40*y*(-12911412490822047011077716 + y*(71942117064863810449017711 + 800*y*(3016805471122684680351 + y*(-7306648318334971703039 + 8*y*(-22005304021302476127 + 13952322443747600185*y)))))) + 6*x*(3583617805503366573901766935 + y*(-539997863237214241774925283 + 40*y*(-32876604623363072658989403 + 2*y*(3430167819779391335035679 + 240*y*(5681012924106078824787 + y*(-1148410005247339716587 + 200*y*(-521520979887370115 + 80705670151250134*y)))))))))/1.5123064061952e29",2); // 2
        iterLimit=5;
        break;
    case 3:
        ss << "3-patch domain, degree 4\n";
        pathStart = "Mario/Florian_ThreePatches_deg4/multipatch_three";
        pathEnd = "non.xml";
        mapPathEnd = "_map.xml";
        g=gsFunctionExpr<real_t>("((-47.7 - 14*x - y)^(2)*(-3.8774193548387097 - (16*x)/31. - y)^2*(3.875 - x/2. - y)^(2)*(-3.7240259740259742 + (39*x)/77. - y)^2*(3.8243055555555556 + (37*x)/72. - y)^2*(-286.75 + 77*x - y)^(2))/2.5e11",2); // 2
        dg=gsFunctionExpr<real_t>("((5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-441389437313011+5272968224115*y+8*(298666368000*(x)^(5)+28000*(x)^(4)*(-4133821+672325*y)-480*(x)^(3)*(52101274671+5*y*(-48329457+640401050*y))+12*(x)^(2)*(565843921171+5*y*(14511590237+20*(46516301-64021780*y)*y))+(y)^(2)*(8076277817723+20*y*(-58621313017+200*y*(-59499360+10828163*y)))+2*x*(229800161083395+y*(2082084269689+10*y*(-1675159861883+20*y*(207397103+3706248990*y)))))))/1.5123064061952e31",
                                   "((5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(31*(-65076108692243+567148797931020*y)+4*(x*(5272968224115+8*x*(2082084269689+20*x*(14511590237+10*x*(144988371+18825100*x))))-16*x*(-8076277817723+10*x*(1675159861883+680*x*(-2736253+56505975*x)))*y+12*(15296479392949-40*x*(58621313017+10*x*(-207397103+128043560*x)))*(y)^(2)+480*(-624025745233+400*x*(-39666240+123541633*x))*(y)^(3)+8000*(-1149261923+216563260*x)*(y)^(4)-32997888000*(y)^(5))))/6.0492256247808e31",2);
        ddg=gsFunctionExpr<real_t>("(2*(5507+740*x-1440*y)^(2)*(1147-308*x+4*y)^(2)*(1147-156*x+308*y)^(2)*(139686329327-200681825360*y+25*(16*x*(-362657382+x*(-101620989+28000*x*(355+84*x)))+8*x*(-373082571+20*x*(4913781+1257200*x))*y+(-653061619+240*x*(5025693+1529365*x))*(y)^(2)+110*(4996423+2317980*x)*(y)^(3)+48905650*(y)^(4)))+8*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-489329699+26666640*(x)^(2)+8*y*(15792077+4325450*y)-8*x*(8077099+8740640*y))*(-330574+75085*y+25*(672*(x)^(2)+931*(y)^(2)+4*x*(355+449*y)))+(-31+4*x+8*y)^(2)*(477+140*x+10*y)^(2)*(601+80*x+155*y)^(2)*(-3441*(-35572636019847+73096675692352*y)+16*(74073925926000*(x)^(4)-44444400*(x)^(3)*(8077099+8740640*y)+4*(y)^(2)*(-172228484745029+20*y*(13798425242429+959704939205*y))+6*(x)^(2)*(-478459593600089+80*y*(3519648832631+1435589920070*y))+3*x*(3988696521770833-8*y*(-206187448112451+10*y*(10414257900957+1896520058240*y))))))/3.780766015488e30",
                                    "(-32*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-184879+179600*(x)^(2)+60*y*(19717+620*y)+40*x*(15017+9310*y))*(-9535011+17481280*(x)^(2)+4*y*(31783721+166320*y)-4*x*(15792077+8650900*y))+1573679923200*(-31+4*x+8*y)^(2)*(477+140*x+10*y)^(2)*(601+80*x+155*y)^(2)*(((11119+113140*x-221760*y)*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y))/(2.45887488e10)+2*(-3.7240259740259742+(39*x)/77.-y)^2*(3.8243055555555556+(37*x)/72.-y)^2+2*(286.75-77*x+y)^(2)*((-350064735599+200*x*(15438523+96003266*x))/(1.22943744e10)-((11119+113140*x)*y)/18480.+6*(y)^(2)))+(5507+740*x-1440*y)^(2)*(1147-308*x+4*y)^(2)*(1147-156*x+308*y)^(2)*(-10479303116099-1317338496540*y+40*(4*x*(-41626408888+5*x*(-653061619+40*x*(10051386+1529365*x)))+330*x*(130950371+20*x*(4996423+772660*x))*y+15*(3452990809+40*x*(58793797+9781130*x))*(y)^(2)+62000*(59151+18620*x)*(y)^(3)+57660000*(y)^(4))))/6.0492256247808e31",
                                    "(263879111196180480000000*(x)^(10)-200704000000*(x)^(9)*(-8649747410767+243871656618475*y)-1102267*(-249908537057869493753722+7303900366010386651975527*y)-38707200000*(x)^(8)*(290585708203906+3*y*(-197029885967433+80994558247885*y))+8192000*(x)^(7)*(-9481688358976200867+5*y*(78477855207928463297+50*y*(-12621679664055075+221243488766886308*y)))+2867200*(x)^(6)*(-280977579189855473287+5*y*(-22493549631768585111+5*y*(-2747278652088440511+20*y*(-91855653718880188+46208541740778675*y))))-92160*(x)^(5)*(75483626813060289160351+25*y*(-20046466182105866760727+100*y*(11594160804861041637+20*y*(-887558443306536565+4*y*(1720896663954649+71090383325622549*y)))))-7680*(x)^(4)*(-4964239415059682722723998+5*y*(2296897469493749744928169+25*y*(4891483627964772371247+40*y*(5266351819578766236+5*y*(-13129448820140085149+80*y*(-26146153193635401+10391360459106683*y))))))+1280*(x)^(3)*(319898962438204378868630125+y*(-3820925583992735261392046391+10*y*(-4783100735284830070720533+20*y*(2105413765947766789130671+40*y*(-281075782689312175923+5*y*(-527029110724154050557+8*y*(1102027571944353971+2738805988130334580*y)))))))+(y)^(2)*(617841001190024474584165568361-80*y*(-21074739527441629699914866180+y*(2607549751484824916038727667+16*y*(91673462096616111761339829+125*y*(-107496658047501834366731+80*y*(-225835710267402201780+y*(31492159069935608973+80*y*(12239879471796649+40941370928304*y))))))))+288*(x)^(2)*(-920093625250723826423869547+10*y*(892691419041710198843985595+y*(-47417212378914652851450213+20*y*(-6547575133692972562670908+5*y*(131944367478582199929979+32*y*(2394158422480936839957+5*y*(-79873367397145988827+200*y*(-56248512753087052+12046162278210363*y))))))))+4*x*(-1197623698331623912863154084953+y*(12321028396612468035229975389749+40*y*(16186340135808633803985384165+4*y*(-13097715815752721726226123093+20*y*(-24778019926500117047201121+y*(43068787960711552962518799+160*y*(8538946273367625114111+10*y*(-695886202065453168197+5*y*(-4066084107601811649+10710407262566600*y))))))))))/7.561532030976e30",2);
        f=gsFunctionExpr<real_t>("(775*(10070636445024446987294496713 + 1580805041435296880732831616*y) + 32*(839212170045134580480000*(x)^(8) + 3072000*(x)^(7)*(-186153527560215013 + 23630822384505615*y) - 12800*(x)^(6)*(7397845103057974754479 + 600*y*(-1978168936115700 + 521427067157777179*y)) - 23040*(x)^(5)*(-2591139538241579069079 + 5*y*(-19502249110340706949 + 40*y*(-532005332078551427 + 231162700790668990*y))) + 240*(x)^(4)*(20629697508869638693714899 + 80*y*(-2832250234210447845948 + 5*y*(361036991177109628031 + 80*y*(851520842693099156 + 1397012256516649535*y)))) + 480*(x)^(3)*(-5011292996454793182398419 + y*(172125474399395900332421 + 80*y*(6788570295304368834777 + 10*y*(-121146247699783170709 + 40*y*(-983617023935690419 + 258889960348724901*y))))) + (y)^(2)*(-113362049203585367125403749833 + 80*y*(-133919062100907890758222638 + y*(119712510000704219793039297 + 16*y*(506422009899324159806439 + 5*y*(-20050351692839232072401 + 120*y*(-4976219400498134306 + 701059526547943945*y)))))) + 6*(x)^(2)*(-13034991506843954692628435037 + 40*y*(-12911412490822047011077716 + y*(71942117064863810449017711 + 800*y*(3016805471122684680351 + y*(-7306648318334971703039 + 8*y*(-22005304021302476127 + 13952322443747600185*y)))))) + 6*x*(3583617805503366573901766935 + y*(-539997863237214241774925283 + 40*y*(-32876604623363072658989403 + 2*y*(3430167819779391335035679 + 240*y*(5681012924106078824787 + y*(-1148410005247339716587 + 200*y*(-521520979887370115 + 80705670151250134*y)))))))))/1.5123064061952e29",2); // 2
        iterLimit=5;
        break;
    case 4:
        ss << "4-patch domain, degree 3\n";
        pathStart = "Mario/Florian_FourPatches_deg3/multipatch_four";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        startIter=1;
        g=gsFunctionExpr<real_t>("((856 + 25*x - 210*y)^(2)*(794 + 185*x - 15*y)^(2)*(253 - 60*x + 5*y)^(2)*(819 + 190*x + 20*y)^(2)*(-772 + 185*x + 25*y)^(2)*(343 - 10*x + 80*y)^(2)*(507 + 10*x + 120*y)^(2)*(-1447 + 60*x + 350*y)^(2))/7.001316e42",2);
        dg=gsFunctionExpr<real_t>("((856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-183266112727112768808+93639600000000*(x)^(7)+388500000*(x)^(6)*(5638646+318315*y)-150000*(x)^(5)*(1137334053991+5*y*(38770857663+14530646650*y))+12500*(x)^(4)*(-152716930485143+25*y*(-156401552027+349513306994*y+7224284480*(y)^(2)))+1000*(x)^(3)*(70984719289425039+5*y*(-118089672499773+50*y*(-29650643938913+293423504802*y+882217595180*(y)^(2))))+75*(x)^(2)*(677492733286137563+5*y*(11600640434709637+10*y*(-981143799014747+10*y*(-12665036490069+50*y*(20291359787+6678463216*y)))))-5*y*(4654806746679084698+5*y*(-595081862571168173+5*y*(-21886279724378373+10*y*(156538493064151+160*y*(377450354536+45*y*(156065933+4364500*y))))))-20*x*(61409432176835309771+5*y*(-231054258707148997+5*y*(-273537297379319971+5*y*(686445838949453+10*y*(155942373681973+100*y*(-5142492283+817279476*y))))))))/7.001316e41",
                                   "((856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-24*(5557475750222879756+52643682537326664905*y)+5*(3533296500000*(x)^(7)-25000*(x)^(6)*(38770857663+29061293300*y)+12500*(x)^(5)*(-156401552027+4*y*(174756653497+5418213360*y))+250*(x)^(4)*(-118089672499773+100*y*(-29650643938913+y*(440135257203+1764435190360*y)))+25*(x)^(3)*(11600640434709637+20*y*(-981143799014747+5*y*(-37995109470207+200*y*(20291359787+8348079020*y))))+30*(y)^(2)*(199526503227539871+50*y*(9930626942414458+5*y*(-17989001088697+24*y*(-141086396887+122500*y*(2701+384*y)))))-10*(x)^(2)*(-231054258707148997+5*y*(-547074594758639942+5*y*(2059337516848359+40*y*(155942373681973+25*y*(-25712461415+4903676856*y)))))-x*(4654806746679084698+5*y*(-1190163725142336346+5*y*(-65658839173135119+40*y*(156538493064151+400*y*(188725177268+9*y*(468197799+15275750*y)))))))))/7.001316e41",2);
        ddg=gsFunctionExpr<real_t>("((819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(4*(-3+100*x-820*y)*(167+4440*x-365*y)*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)+(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(-5872151+15000*(x)^(2)-300*x*(3+820*y)+80*y*(949+12605*y))+(856+25*x-210*y)^(2)*(343-10*x+80*y)^(2)*(-178355327+29570400*(x)^(2)-13320*x*(-167+365*y)+5*y*(-55462+39965*y)))+(856+25*x-210*y)^(2)*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(343-10*x+80*y)^(2)*(40*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(319+24*x+214*y)*(-1447+60*x+350*y)*(967+14060*x+1690*y)+(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(-1777002527+296525400*(x)^(2)+42180*x*(967+1690*y)+60*y*(290448+71035*y))+4*(819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(-6259523+21600*(x)^(2)+1800*x*(319+214*y)+20*y*(172951+82445*y)))+20*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-9685978+444000*(x)^(3)+(x)^(2)*(5070-5515950*y)+y*(54496707+25*(13225-49548*y)*y)+2*x*(-132371273+5*y*(287721+1551850*y)))*(-2726354163+16872000*(x)^(3)+300*(x)^(2)*(1127087+762350*y)+10*y*(-257314061+40*y*(148557+180125*y))+40*x*(-260893019+5*y*(646909+3134030*y))))/1.4002632e41",
                                    "((819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(200*(-71+82*x-672*y)*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(7+73*x-6*y)*(-253+60*x-5*y)+(856+25*x-210*y)^(2)*(343-10*x+80*y)^(2)*(-1204067+199825*(x)^(2)+x*(30560-32850*y)+450*y*(-7+3*y))+4*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(-98526263+252100*(x)^(2)+50400*y*(71+336*y)-20*x*(14807+206640*y)))+(856+25*x-210*y)^(2)*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(343-10*x+80*y)^(2)*(8*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(1007+1690*x+200*y)*(-1447+60*x+350*y)*(381+1070*x+8400*y)+(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(-24276671+4262100*(x)^(2)+600*y*(1007+100*y)+60*x*(59951+16900*y))+4*(819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(1648900*(x)^(2)+180*x*(78963+149800*y)+27*(-22818637+2800*y*(127+1400*y))))-20*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(x*(-54496707+5*x*(-287721+367730*x))-250*x*(2645+62074*x)*y+450*(713+8258*x)*(y)^(2)-201600*(y)^(3)+4*(6617497+67949189*y))*(-3*(406850873+3587051620*y)+10*(7623500*(x)^(3)+7500*(y)^(2)*(1707+224*y)+10*(x)^(2)*(646909+6268060*y)+x*(-257314061+120*y*(99038+180125*y)))))/1.4002632e41",
                                    "(8736*(-5889591170879306018956271073243080922+22925776618634104814916613389531613295*y)+25*(62035588301512500000000000*(x)^(14)-4312350000000000*(x)^(13)*(3324258357687+2740557975655*y)-1803750000000*(x)^(12)*(307513783018520481+110*y*(1043969817180299+105368781544875*y))+18750000000*(x)^(11)*(1616209809586420291851+y*(1238021454387970462733+2*y*(183253666887968429619+58954243271100176660*y)))+137500000*(x)^(10)*(3860051208614600588030816+5*y*(-305396154884365519380817+50*y*(-2492215107508311154383+20*y*(-30895915857524455283+324570185189587235*y))))-6250000*(x)^(9)*(2469981477768200761696850107+5*y*(-170548179992646678377359407+10*y*(-15860356643192086745823063+20*y*(-245996004761701292488629+100*y*(462499575852494684963+191973545558897389557*y)))))+5625000*(x)^(8)*(-7933898947401529305690216499+10*y*(7139097275043354971244075328+y*(236968812099024590167053471+20*y*(-38254237482625477235750804+25*y*(-28832614495174597470411+8*y*(5612513556492807328773+96858236340211903900*y))))))+1250000*(x)^(7)*(-90008106283614443392212407123+y*(-8432585990927600005430502072037+2*y*(73725322169426088714622191129+20*y*(33267867471633597699708373449+50*y*(-11859166675181326539831076+y*(-37751660152923343565540931+10*y*(49339595481767764811939+78022255196391533949800*y)))))))+43750*(x)^(6)*(84576571841663318667353827852358+5*y*(-119270021041172613069243491657009+2*y*(-5267006730739703066455523697141+20*y*(365096653088222551419951669109+50*y*(720304994113427524676916996+y*(-262145884199479055265345723+10*y*(-3766966184007510154199879+40*y*(4858565047001895189461+1654278215577084353790*y))))))))-9375*(x)^(5)*(-3944665406579716882953309808152639+y*(-59631788330605795969226108804524601+2*y*(1175294921656112003132297904715299+20*y*(247260307856921203334807082003489+50*y*(-133087026003182851788371621127+y*(-280966662222325225695067978251+10*y*(462831431276262238433359389+1000*y*(543109951341832710961281+5*y*(-445873082654359563939+60986370928455606104*y)))))))))-1250*(x)^(4)*(86667944396611326125152191968589291+5*y*(-92711306107668717189746378411325493+2*y*(-6137440113826370431243828948588794+5*y*(1221053113991921726225936571149394+25*y*(7199094000253875689933523616543+4*y*(-487098523544301830556854507817+5*y*(-19400717898846330756546159071+20*y*(48841879472535873837618516+5*y*(3460946741435033281024677+200*y*(102652295928016054801+11312713918519394172*y))))))))))+20*(x)^(3)*(-43716237256277211654117871080508759761+5*y*(-100504857453883898313383545165878863533+10*y*(494455112403789248890687184221654578+5*y*(339712344190683828714051811107583207+25*y*(-435352475067945137636850632278106+y*(-775588456340703311348419990354869+10*y*(1551638818925981834611256221577+40*y*(37254802836739713592089212011+60*y*(-5055376602138171854610588+25*y*(5433891171018067058939+4*y*(-30756759099767012807+2078936708753422440*y)))))))))))-3*(y)^(2)*(-517187094791941512207875623483647488772+5*y*(77217754052102660111443949529475785072+25*y*(633448828506561189543976944087298763+5*y*(-27988221239870592477014513717078484+y*(-7102809297961140075487637866868499+20*y*(17876842976762970788426265473984+5*y*(1368773374720284161388360799551+200*y*(64129244261551769175819087+32*y*(-194672732009669889214129+90*y*(-68395720077821388488+13125*y*(-68241190290897+784*y*(2293809571+56115000*y))))))))))))+3*(x)^(2)*(357008814843850566704563621317729509638+5*y*(-316328078054897574117811933975431492486+5*y*(-10541312754448979091768159561196028679+5*y*(1754523182153765065935875096330768944+25*y*(12679802486482816224180486557374431+10*y*(-300432851590565251710084441841767+2*y*(-34770495496633093130101989295657+40*y*(43348784030460636922030771701+5*y*(3181431815859714378054765579+200*y*(131321929928275706385031+18*y*(273910185345199868311+240*y*(31713612288757402+1497839540845375*y))))))))))))+x*(5975401493509874489718402548998234307088-25*y*(-2425969083067941644958089230894270786564+y*(134823314409260395756442128075036963626+y*(415920671685559246029616085390606549878+25*y*(-600479315071868234907855631100928349+y*(-961533574482894251649317866257559143+10*y*(2260764901252985584405566308733511+200*y*(9504643318524626632512636282979+y*(-116341862738567034959898413679+800*y*(-9882119516256426019546181+6*y*(-2056855346844450996833+180*y*(-6902072649829414053+1750*y*(-12582105727243+1373281449960*y)))))))))))))))/7.001316e40",2);
        f=gsFunctionExpr<real_t>("(200483015939927424805749401753922095656048-20792963040042998776819931131812603569760*y+5*(85796184373560849187500000000*(x)^(12)-7500000000*(x)^(11)*(188763309224212970923+5082269320562618300*y)-1875000000*(x)^(10)*(-69162007436240087786159+8*y*(2773202004406440477929+3894193987465008424545*y))+25000000*(x)^(9)*(-149909623908054866557383907+400*y*(-63214168016399689867944+5*y*(14015207218248956806026+598454433048046209625*y)))+7500000*(x)^(8)*(19473177179241234640028364387+50*y*(31862975136666387602469487+10*y*(-16355455723645337460846387+10*y*(-112982613228639206194+109202032280810305350935*y))))+750000*(x)^(7)*(1572537911999558169581729836189+20*y*(12594309484008206716874150172+5*y*(-4165815250858270258696456113+100*y*(-4733092593717614720471873+40*y*(40903113642087276773771+5340069121376618277447*y)))))+12500*(x)^(6)*(-2271988464656030168577251066789687+20*y*(2944743138851483077922538850602+5*y*(9910947021964690206354762839607+700*y*(-292940854633575377776655719+40*y*(-28505650653407207205195069+279668611635311480169144*y+503637821082003617841650*(y)^(2))))))+7500*(x)^(5)*(-11404202153620510282266286367953621+40*y*(-41430614895086613420012978375721+10*y*(7300076807925748182243512523507+5*y*(207798978778270540024515420577+5*y*(-15704605651331332650308349669+20*y*(-131192453670640178453814723+100*y*(90447654906181480116769+18546532481001706820208*y)))))))+375*(x)^(4)*(3679660093973470770621986393970644679+100*y*(-2109005546190532388780547854257611+5*y*(-2849431322847162531737397627418393+20*y*(4557430842446960617713592202719+5*y*(2392150695428852486560506087695+4*y*(-10301361341987745035225899043+50*y*(-276806349418032697590524561+100*y*(18342483418404099256747+12322864260311585068187*y))))))))+100*(x)^(3)*(15056257746185269082445666631937465423+125*y*(20980001691361741808527279021965765+2*y*(-15841241175542331931316536052702589+10*y*(-285162413610154714654165106097818+15*y*(6409883038339352313893789206357+4*y*(333508619420475236878379216642+25*y*(-1023407103086551980758209087+20*y*(-12877828294431805893821999+5*y*(22595898331228390795931+8841866686884050087760*y)))))))))+5*(y)^(2)*(-3551531163988261907495742883855209660492+5*y*(34381021327838633385934483932486896924+5*y*(11829420487260040622572021242890193237+20*y*(-16049295832071177139434109298296377+5*y*(-2419016927062141523648868713984609+300*y*(122010981956094749795023301659+y*(44005037843774091904242393211+80*y*(-3472649528000044723949072+15*y*(7892045651005192672523+40*y*(-2397054208507204837+142550086209543870*y))))))))))-30*(x)^(2)*(565149859396440035224265578183981607208+5*y*(-8437287002057048996772266807335117994+5*y*(-8716574128893158922537895479104568879+50*y*(6918879623452658688632669271673658+5*y*(2958282868548894096977145766528735+4*y*(-16560015850331253540195841546284+5*y*(-3557510113118293272483686100191+60*y*(652163400886916781761946087+125*y*(2624182264153707638925741+8*y*(-999457787323601761019+135161035690724068988*y))))))))))-12*x*(530925557520583884799110854193758960732+5*y*(24376442401877963241185423320934974674+5*y*(-5999125722417512746535516369169733891+25*y*(-54851286351950879916482325313472751+5*y*(3098631074456770371333310355311893+8*y*(98897719114710296309142369775414+5*y*(-1337911001230730072008322933889+10*y*(-40430425105286506983844478437+250*y*(1301566201057324263051717+2*y*(363125623376994231239461+4*y*(377139568951883178101+44110620057322458600*y)))))))))))))/1.4002632e39",2);
        iterLimit=5;
        break;
    case 5:
        ss << "4-patch domain, degree 4\n";
        pathStart = "Mario/Florian_FourPatches_deg4/multipatch_four";
        pathEnd = "non.xml";
        mapPathEnd = "_map.xml";
        g=gsFunctionExpr<real_t>("((856 + 25*x - 210*y)^(2)*(794 + 185*x - 15*y)^(2)*(253 - 60*x + 5*y)^(2)*(819 + 190*x + 20*y)^(2)*(-772 + 185*x + 25*y)^(2)*(343 - 10*x + 80*y)^(2)*(507 + 10*x + 120*y)^(2)*(-1447 + 60*x + 350*y)^(2))/7.001316e42",2);
        dg=gsFunctionExpr<real_t>("((856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-183266112727112768808+93639600000000*(x)^(7)+388500000*(x)^(6)*(5638646+318315*y)-150000*(x)^(5)*(1137334053991+5*y*(38770857663+14530646650*y))+12500*(x)^(4)*(-152716930485143+25*y*(-156401552027+349513306994*y+7224284480*(y)^(2)))+1000*(x)^(3)*(70984719289425039+5*y*(-118089672499773+50*y*(-29650643938913+293423504802*y+882217595180*(y)^(2))))+75*(x)^(2)*(677492733286137563+5*y*(11600640434709637+10*y*(-981143799014747+10*y*(-12665036490069+50*y*(20291359787+6678463216*y)))))-5*y*(4654806746679084698+5*y*(-595081862571168173+5*y*(-21886279724378373+10*y*(156538493064151+160*y*(377450354536+45*y*(156065933+4364500*y))))))-20*x*(61409432176835309771+5*y*(-231054258707148997+5*y*(-273537297379319971+5*y*(686445838949453+10*y*(155942373681973+100*y*(-5142492283+817279476*y))))))))/7.001316e41",
                                   "((856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-24*(5557475750222879756+52643682537326664905*y)+5*(3533296500000*(x)^(7)-25000*(x)^(6)*(38770857663+29061293300*y)+12500*(x)^(5)*(-156401552027+4*y*(174756653497+5418213360*y))+250*(x)^(4)*(-118089672499773+100*y*(-29650643938913+y*(440135257203+1764435190360*y)))+25*(x)^(3)*(11600640434709637+20*y*(-981143799014747+5*y*(-37995109470207+200*y*(20291359787+8348079020*y))))+30*(y)^(2)*(199526503227539871+50*y*(9930626942414458+5*y*(-17989001088697+24*y*(-141086396887+122500*y*(2701+384*y)))))-10*(x)^(2)*(-231054258707148997+5*y*(-547074594758639942+5*y*(2059337516848359+40*y*(155942373681973+25*y*(-25712461415+4903676856*y)))))-x*(4654806746679084698+5*y*(-1190163725142336346+5*y*(-65658839173135119+40*y*(156538493064151+400*y*(188725177268+9*y*(468197799+15275750*y)))))))))/7.001316e41",2);
        ddg=gsFunctionExpr<real_t>("((819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(4*(-3+100*x-820*y)*(167+4440*x-365*y)*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)+(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(-5872151+15000*(x)^(2)-300*x*(3+820*y)+80*y*(949+12605*y))+(856+25*x-210*y)^(2)*(343-10*x+80*y)^(2)*(-178355327+29570400*(x)^(2)-13320*x*(-167+365*y)+5*y*(-55462+39965*y)))+(856+25*x-210*y)^(2)*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(343-10*x+80*y)^(2)*(40*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(319+24*x+214*y)*(-1447+60*x+350*y)*(967+14060*x+1690*y)+(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(-1777002527+296525400*(x)^(2)+42180*x*(967+1690*y)+60*y*(290448+71035*y))+4*(819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(-6259523+21600*(x)^(2)+1800*x*(319+214*y)+20*y*(172951+82445*y)))+20*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-9685978+444000*(x)^(3)+(x)^(2)*(5070-5515950*y)+y*(54496707+25*(13225-49548*y)*y)+2*x*(-132371273+5*y*(287721+1551850*y)))*(-2726354163+16872000*(x)^(3)+300*(x)^(2)*(1127087+762350*y)+10*y*(-257314061+40*y*(148557+180125*y))+40*x*(-260893019+5*y*(646909+3134030*y))))/1.4002632e41",
                                    "((819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(200*(-71+82*x-672*y)*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(7+73*x-6*y)*(-253+60*x-5*y)+(856+25*x-210*y)^(2)*(343-10*x+80*y)^(2)*(-1204067+199825*(x)^(2)+x*(30560-32850*y)+450*y*(-7+3*y))+4*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(-98526263+252100*(x)^(2)+50400*y*(71+336*y)-20*x*(14807+206640*y)))+(856+25*x-210*y)^(2)*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(343-10*x+80*y)^(2)*(8*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(1007+1690*x+200*y)*(-1447+60*x+350*y)*(381+1070*x+8400*y)+(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(-24276671+4262100*(x)^(2)+600*y*(1007+100*y)+60*x*(59951+16900*y))+4*(819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(1648900*(x)^(2)+180*x*(78963+149800*y)+27*(-22818637+2800*y*(127+1400*y))))-20*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(x*(-54496707+5*x*(-287721+367730*x))-250*x*(2645+62074*x)*y+450*(713+8258*x)*(y)^(2)-201600*(y)^(3)+4*(6617497+67949189*y))*(-3*(406850873+3587051620*y)+10*(7623500*(x)^(3)+7500*(y)^(2)*(1707+224*y)+10*(x)^(2)*(646909+6268060*y)+x*(-257314061+120*y*(99038+180125*y)))))/1.4002632e41",
                                    "(8736*(-5889591170879306018956271073243080922+22925776618634104814916613389531613295*y)+25*(62035588301512500000000000*(x)^(14)-4312350000000000*(x)^(13)*(3324258357687+2740557975655*y)-1803750000000*(x)^(12)*(307513783018520481+110*y*(1043969817180299+105368781544875*y))+18750000000*(x)^(11)*(1616209809586420291851+y*(1238021454387970462733+2*y*(183253666887968429619+58954243271100176660*y)))+137500000*(x)^(10)*(3860051208614600588030816+5*y*(-305396154884365519380817+50*y*(-2492215107508311154383+20*y*(-30895915857524455283+324570185189587235*y))))-6250000*(x)^(9)*(2469981477768200761696850107+5*y*(-170548179992646678377359407+10*y*(-15860356643192086745823063+20*y*(-245996004761701292488629+100*y*(462499575852494684963+191973545558897389557*y)))))+5625000*(x)^(8)*(-7933898947401529305690216499+10*y*(7139097275043354971244075328+y*(236968812099024590167053471+20*y*(-38254237482625477235750804+25*y*(-28832614495174597470411+8*y*(5612513556492807328773+96858236340211903900*y))))))+1250000*(x)^(7)*(-90008106283614443392212407123+y*(-8432585990927600005430502072037+2*y*(73725322169426088714622191129+20*y*(33267867471633597699708373449+50*y*(-11859166675181326539831076+y*(-37751660152923343565540931+10*y*(49339595481767764811939+78022255196391533949800*y)))))))+43750*(x)^(6)*(84576571841663318667353827852358+5*y*(-119270021041172613069243491657009+2*y*(-5267006730739703066455523697141+20*y*(365096653088222551419951669109+50*y*(720304994113427524676916996+y*(-262145884199479055265345723+10*y*(-3766966184007510154199879+40*y*(4858565047001895189461+1654278215577084353790*y))))))))-9375*(x)^(5)*(-3944665406579716882953309808152639+y*(-59631788330605795969226108804524601+2*y*(1175294921656112003132297904715299+20*y*(247260307856921203334807082003489+50*y*(-133087026003182851788371621127+y*(-280966662222325225695067978251+10*y*(462831431276262238433359389+1000*y*(543109951341832710961281+5*y*(-445873082654359563939+60986370928455606104*y)))))))))-1250*(x)^(4)*(86667944396611326125152191968589291+5*y*(-92711306107668717189746378411325493+2*y*(-6137440113826370431243828948588794+5*y*(1221053113991921726225936571149394+25*y*(7199094000253875689933523616543+4*y*(-487098523544301830556854507817+5*y*(-19400717898846330756546159071+20*y*(48841879472535873837618516+5*y*(3460946741435033281024677+200*y*(102652295928016054801+11312713918519394172*y))))))))))+20*(x)^(3)*(-43716237256277211654117871080508759761+5*y*(-100504857453883898313383545165878863533+10*y*(494455112403789248890687184221654578+5*y*(339712344190683828714051811107583207+25*y*(-435352475067945137636850632278106+y*(-775588456340703311348419990354869+10*y*(1551638818925981834611256221577+40*y*(37254802836739713592089212011+60*y*(-5055376602138171854610588+25*y*(5433891171018067058939+4*y*(-30756759099767012807+2078936708753422440*y)))))))))))-3*(y)^(2)*(-517187094791941512207875623483647488772+5*y*(77217754052102660111443949529475785072+25*y*(633448828506561189543976944087298763+5*y*(-27988221239870592477014513717078484+y*(-7102809297961140075487637866868499+20*y*(17876842976762970788426265473984+5*y*(1368773374720284161388360799551+200*y*(64129244261551769175819087+32*y*(-194672732009669889214129+90*y*(-68395720077821388488+13125*y*(-68241190290897+784*y*(2293809571+56115000*y))))))))))))+3*(x)^(2)*(357008814843850566704563621317729509638+5*y*(-316328078054897574117811933975431492486+5*y*(-10541312754448979091768159561196028679+5*y*(1754523182153765065935875096330768944+25*y*(12679802486482816224180486557374431+10*y*(-300432851590565251710084441841767+2*y*(-34770495496633093130101989295657+40*y*(43348784030460636922030771701+5*y*(3181431815859714378054765579+200*y*(131321929928275706385031+18*y*(273910185345199868311+240*y*(31713612288757402+1497839540845375*y))))))))))))+x*(5975401493509874489718402548998234307088-25*y*(-2425969083067941644958089230894270786564+y*(134823314409260395756442128075036963626+y*(415920671685559246029616085390606549878+25*y*(-600479315071868234907855631100928349+y*(-961533574482894251649317866257559143+10*y*(2260764901252985584405566308733511+200*y*(9504643318524626632512636282979+y*(-116341862738567034959898413679+800*y*(-9882119516256426019546181+6*y*(-2056855346844450996833+180*y*(-6902072649829414053+1750*y*(-12582105727243+1373281449960*y)))))))))))))))/7.001316e40",2);
        f=gsFunctionExpr<real_t>("(200483015939927424805749401753922095656048-20792963040042998776819931131812603569760*y+5*(85796184373560849187500000000*(x)^(12)-7500000000*(x)^(11)*(188763309224212970923+5082269320562618300*y)-1875000000*(x)^(10)*(-69162007436240087786159+8*y*(2773202004406440477929+3894193987465008424545*y))+25000000*(x)^(9)*(-149909623908054866557383907+400*y*(-63214168016399689867944+5*y*(14015207218248956806026+598454433048046209625*y)))+7500000*(x)^(8)*(19473177179241234640028364387+50*y*(31862975136666387602469487+10*y*(-16355455723645337460846387+10*y*(-112982613228639206194+109202032280810305350935*y))))+750000*(x)^(7)*(1572537911999558169581729836189+20*y*(12594309484008206716874150172+5*y*(-4165815250858270258696456113+100*y*(-4733092593717614720471873+40*y*(40903113642087276773771+5340069121376618277447*y)))))+12500*(x)^(6)*(-2271988464656030168577251066789687+20*y*(2944743138851483077922538850602+5*y*(9910947021964690206354762839607+700*y*(-292940854633575377776655719+40*y*(-28505650653407207205195069+279668611635311480169144*y+503637821082003617841650*(y)^(2))))))+7500*(x)^(5)*(-11404202153620510282266286367953621+40*y*(-41430614895086613420012978375721+10*y*(7300076807925748182243512523507+5*y*(207798978778270540024515420577+5*y*(-15704605651331332650308349669+20*y*(-131192453670640178453814723+100*y*(90447654906181480116769+18546532481001706820208*y)))))))+375*(x)^(4)*(3679660093973470770621986393970644679+100*y*(-2109005546190532388780547854257611+5*y*(-2849431322847162531737397627418393+20*y*(4557430842446960617713592202719+5*y*(2392150695428852486560506087695+4*y*(-10301361341987745035225899043+50*y*(-276806349418032697590524561+100*y*(18342483418404099256747+12322864260311585068187*y))))))))+100*(x)^(3)*(15056257746185269082445666631937465423+125*y*(20980001691361741808527279021965765+2*y*(-15841241175542331931316536052702589+10*y*(-285162413610154714654165106097818+15*y*(6409883038339352313893789206357+4*y*(333508619420475236878379216642+25*y*(-1023407103086551980758209087+20*y*(-12877828294431805893821999+5*y*(22595898331228390795931+8841866686884050087760*y)))))))))+5*(y)^(2)*(-3551531163988261907495742883855209660492+5*y*(34381021327838633385934483932486896924+5*y*(11829420487260040622572021242890193237+20*y*(-16049295832071177139434109298296377+5*y*(-2419016927062141523648868713984609+300*y*(122010981956094749795023301659+y*(44005037843774091904242393211+80*y*(-3472649528000044723949072+15*y*(7892045651005192672523+40*y*(-2397054208507204837+142550086209543870*y))))))))))-30*(x)^(2)*(565149859396440035224265578183981607208+5*y*(-8437287002057048996772266807335117994+5*y*(-8716574128893158922537895479104568879+50*y*(6918879623452658688632669271673658+5*y*(2958282868548894096977145766528735+4*y*(-16560015850331253540195841546284+5*y*(-3557510113118293272483686100191+60*y*(652163400886916781761946087+125*y*(2624182264153707638925741+8*y*(-999457787323601761019+135161035690724068988*y))))))))))-12*x*(530925557520583884799110854193758960732+5*y*(24376442401877963241185423320934974674+5*y*(-5999125722417512746535516369169733891+25*y*(-54851286351950879916482325313472751+5*y*(3098631074456770371333310355311893+8*y*(98897719114710296309142369775414+5*y*(-1337911001230730072008322933889+10*y*(-40430425105286506983844478437+250*y*(1301566201057324263051717+2*y*(363125623376994231239461+4*y*(377139568951883178101+44110620057322458600*y)))))))))))))/1.4002632e39",2);
        iterLimit=5;
        break;
    case 6:
        ss << "5-patch domain, degree 3\n";
        pathStart = "Mario/Florian_FivePatches_deg3/multipatch_five";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        startIter=1;
        g=gsFunctionExpr<real_t>("((51.2 - (32*x)/3. - y)^2*(-17.9 - (23*x)/5. - y)^2*(-6.855555555555555 - (13*x)/9. - y)^2*(3.7956521739130435 - (10*x)/23. - y)^2*(-4.41578947368421 - (2*x)/19. - y)^2*(3.0448275862068965 + x/29. - y)^(2)*(4.45 + x/2. - y)^(2)*(-4.923529411764706 + (9*x)/17. - y)^2*(6.1 + x - y)^(2)*(-21.6 + (9*x)/2. - y)^2)/1.e16",2); // 3
        dg=gsFunctionExpr<real_t>("((883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(5184*(-24+5*x)*(-60759448049004644133+10*x*(9224871288649613571+100*x*(37604179460108437+10*x*(-85805440864572+5*x*(-23585288349969+2*x*(-344361452181+40*x*(2225326239+88527705*x+747500*(x)^(2))))))))-135*(17023674641706872782653+20*x*(-1695509837688518024354+5*x*(-41960889498677262899+100*x*(284438575105026496+5*x*(2833672964522191+12*x*(-96263698034234+5*x*(-270006365587+80*x*(969070628+45871135*x))))))))*y-250*(3209323891457105431287+4*x*(-1063905721177541644621+20*x*(-10863701126609174361+5*x*(311839709514378744+25*x*(4409023124019543+4*x*(-39053380365891+50*x*(-47676560757+30484236416*x)))))))*(y)^(2)+500*(414137231579450182133+100*x*(-8503573356019086696+5*x*(-169890059042411955+4*x*(19655998498441548+25*x*(66279852581471+60*x*(-104658896724+22349641013*x))))))*(y)^(3)+25000*(669962652693164829+4*x*(-252981723406319757+20*x*(-3480355330614159+5*x*(-149651754464292+25*x*(1579945392551+246871418796*x)))))*(y)^(4)-50000*(36497337774982877+20*x*(-13208323166449558+15*x*(-26473933807981+40*x*(345529412146+198697687975*x))))*(y)^(5)-500000*(-1618566060310711+60*x*(-11342941558297+40*x*(-39209104117+12314778830*x)))*(y)^(6)+15000000*(-8007076526467+20*x*(-437203258412+139024822065*x))*(y)^(7)+50000000*(-475397417077+52105762340*x)*(y)^(8)-621668112000000000*(y)^(9)))/1.69182035233245e49",
                                   "-((883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(9*(-24+5*x)*(-16077059186642225498133+4*x*(-1369336665190955685819+5*x*(48913670472575889613+20*x*(946609999700553672+5*x*(-5001444039303987+100*x*(-28130131110147+x*(156037145859+40*x*(2620771809+91742270*x))))))))+20*(145441291720013932895181+5*x*(3209323891457105431287+2*x*(-1063905721177541644621+40*x*(-3621233708869724787+5*x*(77959927378594686+5*x*(4409023124019543+10*x*(-13017793455297+100*x*(-6810937251+3810529552*x))))))))*y-300*(4238909567658291807147+x*(414137231579450182133+100*x*(-4251786678009543348+5*x*(-56630019680803985+4*x*(4913999624610387+5*x*(66279852581471+300*x*(-17443149454+3192805859*x)))))))*(y)^(2)-4000*(75855763148894346513+5*x*(669962652693164829+2*x*(-252981723406319757+40*x*(-1160118443538053+5*x*(-37412938616073+5*x*(1579945392551+205726182330*x))))))*(y)^(3)+50000*(1644683832982118973+x*(36497337774982877+20*x*(-6604161583224779+5*x*(-26473933807981+60*x*(172764706073+79479075190*x)))))*(y)^(4)+600000*(18992823414703953+x*(-1618566060310711+10*x*(-34028824674891+40*x*(-78418208234+18472168245*x))))*(y)^(5)-21000000*(64762860701939+x*(-8007076526467+20*x*(-218601629206+46341607355*x)))*(y)^(6)-80000000*(1985626097571+x*(-475397417077+26052881170*x))*(y)^(7)+2700000000*(-500083579+414445408*x)*(y)^(8)+116338140000000000*(y)^(9)))/3.3836407046649e48",2);
        ddg=gsFunctionExpr<real_t>("(80*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(9*(-543914379+1936272910*y)+100*(162*x*(-1003582+5*x*(-993+3568*x+50*(x)^(2)))-18*x*(-3185437+643605*x+61400*(x)^(2))*y+(-59423569+20*x*(601672+190725*x))*(y)^(2)-50*(25691+84020*x)*(y)^(3)+1294400*(y)^(4)))*(2619201987768-935774142455*y+50*(16*x*(-2700305626+x*(-469551171+80*x*(1108018+37375*x)))+x*(-11430373379+20*x*(132519969+24500680*x))*y+4*(-1888260208+5*x*(102280583+48691290*x))*(y)^(2)+10*(2825021+47307530*x)*(y)^(3)+50886600*(y)^(4)))+25*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*((883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(837-90*x+170*y)^(2)*(81*(-5687+60*x*(13+10*x))+180*(167-330*x)*y+15700*(y)^(2))+8*(883+10*x-290*y)*(-837+90*x-170*y)*(117+180*x-110*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(9*(-11809+60*x*(293+5*x))-200*(1485+296*x)*y+104900*(y)^(2))+4*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(27*(-38127151491+30252883840*y)+40*(243*x*(-8844351+5*x*(749023+50*x*(586+5*x)))-360*x*(7615641+10*x*(217413+3700*x))*y+10*(-75702173+60*x*(4329504+266245*x))*(y)^(2)-300*(7743339+1596890*x)*(y)^(3)+421028250*(y)^(4))))+4*(883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(837-90*x+170*y)^(2)*(20*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(-3344+7360*x+1145*y)*(-5422463+240*x*(8226+325*x)+2406380*y+649600*x*y+780500*(y)^(2))+(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*(256*(-1932479+1380*x*(-209+230*x))+320*(-81373+79005*x)*y+1863025*(y)^(2))+25*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(-59816709411407-65166891181640*y+40*(6*x*(-163190328849+20*x*(865547273+6500*x*(16452+325*x)))+20*x*(3646882989+40*x*(143693421+5278000*x))*y+5*(-325804009+360*x*(135707027+12173420*x))*(y)^(2)+300*(536897933+131866100*x)*(y)^(3)+21616698250*(y)^(4)))))/1.69182035233245e49",
                                    "(-40*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(8106539481+9*x*(-193627291+10*x*(-3185437+10*x*(42907+3070*x)))+171788588*y-20*x*(-59423569+20*x*(300836+63575*x))*y+60*(-19223711+25*x*(25691+42010*x))*(y)^(2)-800*(-127477+64720*x)*(y)^(3)+9860000*(y)^(4))*(1804179436833-1545589725720*y+4*(x*(-187154828491+5*x*(-11430373379+40*x*(44173323+6125170*x)))+40*x*(-3776520416+5*x*(102280583+32460860*x))*y+300*(-347459898+x*(2825021+23653765*x))*(y)^(2)+6000*(-507835+339244*x)*(y)^(3)+147487500*(y)^(4)))+4*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*((883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(837-90*x+170*y)^(2)*(-2327+3925*(x)^(2)+300*y*(31+2*y)-20*x*(794+165*y))+400*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(-31+11*x-4*y)*(-115123+1480*(x)^(2)+x*(14850-10490*y)+3*y*(-8451+4930*y))+2*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(983008657193+3696921147420*y+20*(26624500*(x)^(4)+600*(x)^(3)*(1443168-798445*y)+15*(y)^(2)*(-3897676511+98600*y*(-8451+2465*y))+10*(x)^(2)*(-75702173+90*y*(-7743339+2806855*y))+x*(-70039005912+10*y*(2650634259+10*(167638491-51715700*y)*y)))))+(883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(837-90*x+170*y)^(2)*(8*(179+46*x+10*y)*(-768+160*x+15*y)*(-999+458*x+60*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(-4919953+324800*(x)^(2)+20*y*(294019+58995*y)+20*x*(120319+78050*y))+(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*(298084*(x)^(2)+60*x*(-16589+1374*y)+9*(-72407+60*y*(-333+10*y)))+4*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(-241540077242153-193437613834440*y+40*(3652026000*(x)^(4)+600*(x)^(3)*(135707027+65933050*y)+15*(y)^(2)*(47746822063+65550*y*(588038+58995*y))+5*(x)^(2)*(-325804009+60*y*(1610693799+432333965*y))+x*(-3152734415003+30*y*(5512681989+10*y*(3241247549+511617750*y)))))))/6.7672814093298e47",
                                    "(-11664*(24-5*x)^(2)*(3795508851241428081609835296757521036403467+8*x*(-617924889466712236003969418712490006166196+5*x*(-117799961390115888244607713916466201335763+10*x*(-1153754925035254244930228232162930630276+5*x*(190783077681420820239961284719819672741+400*x*(94136542817661227100453777189555504+5*x*(-2335776705045331973430661489837887+5*x*(-139662379522414132011851610140952+x*(2989945885086361921720201255521+40*x*(55027619254438667298341299236+5*x*(193110516368558089453986589+10*x*(-7645478313504331939059012+5*x*(-65763204980290218177741+160*x*(19183304824653527736+5*x*(350631176011193151+40*x*(197697136981902+1302969589675*x))))))))))))))))-3645*(-24+5*x)*(-133217222313490691110038761522914324054767573+4*x*(30707381976125011223697233902455332398970319+20*x*(1191165278871544872251481597207897295736307+25*x*(6059617791835945705996637645576510043311+5*x*(-587254050168334401901216944803729861137+80*x*(-2857926878484637578930592938424051364+5*x*(-1924196692092154206887727370076932+x*(24731098175816427171948650852354142+5*x*(298808202379124496086373676196557+20*x*(-4061264643115522705627614150193+20*x*(-9322530128165611705629893263+5*x*(254749554884573818720518079+5*x*(-2504987654656513967213119+40*x*(-3735831775608908438463+80*x*(11411323124051659378+15*x*(23233173342897719+23271419454370*x))))))))))))))))*y+675*(25695529994036350691195886821931982550755068915+2*x*(-21222561349774070294281120924390666690935287463+5*x*(-1470192296359704686899216335439906402872612009+20*x*(32828438928021100531095137055388422092801442+25*x*(375771878930388547769390107169649955359629+8*x*(-9297559430780970006311141091744818659836+5*x*(-518824950630865729493051643271960121031+40*x*(1239337523806790172465481593355752474+5*x*(70298941463740785448635833708320698+25*x*(-165266538000479350228250648119671+x*(-34434605624479501963715560448477+20*x*(127192476002036628020827168206+5*x*(1247946149098604414584565817+200*x*(-1222897611049032857462391+10*x*(1161178534000944312057+320*x*(1266531504681501663+41525418232972513*x))))))))))))))))*(y)^(2)+1000*(454908891261705405293290728846770879955211962+x*(1192377384387916217783146947465730421600899239+10*x*(4323767282325281414116112562531325444590531+20*x*(-3432367826491467794474142636672332254331339+200*x*(-829274239270937985012716379428979999543+x*(1866368217671403190852503064272450591753+5*x*(14660694531649725548851755369383919711+20*x*(-814066796197809791333767300226627538+25*x*(-845706731082957996585568457594364+5*x*(136341428904494516644467134448151+2*x*(-523686404450488294038346109853+100*x*(-6401149892629696416060739017+20*x*(17276974088396173674640689+x*(893137738722155892236061+100*x*(-308595463037841925575+35323694640674483104*x)))))))))))))))*(y)^(3)-12500*(203547927545787429973647834403689761163930297+4*x*(-88700405875258812633570461844349224850692594+5*x*(-5088884454361057926285129444449369229933843+20*x*(102037029733931730830901970205868315022536+25*x*(1139056577223180365690095879519796153237+20*x*(-5997254477044199421896067464297433513+4*x*(-675585137903823909090519779730570562+25*x*(-395141985092314217490706644985944+x*(588273024309087438974071601399403+100*x*(304204911714318416453141128308+x*(-62276343587707471314832158463+20*x*(-20909933786942601669098856+5*x*(2866213319543382630659383+60*x*(512820645506886052273+550407119691908322750*x))))))))))))))*(y)^(4)-150000*(-820912876560040055454432453797226841034706+5*x*(592450044345189421336013296999840871075323+80*x*(1311506726167945985871175023744316137447+x*(-1803079355393752025371711836650070599656+25*x*(-7267115612631854965848892570016054753+x*(4807123355570829949579856201052443019+80*x*(6152635135507867630226220089665361+5*x*(-309672464598060704832373938469532+5*x*(-7725830568381436150542669393672+25*x*(30980980740115059043000786097+16*x*(172350682542795431883899746+25*x*(-796881288663087402354876+5*x*(7391466654891332268515+236379089094056247129*x)))))))))))))*(y)^(5)+1750000*(80773870313750739897812091863957297184759+2*x*(-92643044110063114738520207042823764986917+25*x*(-719117177967701738403704851801130970015+4*x*(77634810668446222272798278038044225162+25*x*(766259340515986553579744402873911759+16*x*(-1214932437321533156667985363337886+125*x*(-15564052439637071236626757180113+2*x*(-844742087491991094270169700424+x*(80097624826001896502857819221+50*x*(405967284015213374071442325+x*(-10201816907430873985157249+60*x*(12140065785830258340462+8268842088053710100275*x))))))))))))*(y)^(6)-20000000*(539426160540983772054245755685428914630+x*(-1790941405357047567376315431779649452431+10*x*(-18091672375436750457831888150781079637+20*x*(1591626137629073876880994869869520241+400*x*(263805077291529704155031159940002+x*(-168796794272635732890337453811646+5*x*(-4533024644089582279726522153397+50*x*(6604069465120570812393238904+x*(2427734883696585912749793852+25*x*(5319031552190841163885487-970751473560998415598774*x+54273906819203417628840*(x)^(2)))))))))))*(y)^(7)-225000000*(8279588480765903113524666213672961341+4*x*(-14802325646602365620237113770635840653+40*x*(-24629286857699112254909947458500148+5*x*(4750519337020804929328849712882428+25*x*(34961665971044349794100633258143+4*x*(-97995562129368317132554465683+10*x*(-34798545644587604875845990133+15*x*(-276055521378564648930042248+5*x*(-3817133550493210740612501+100*x*(13829622398451154550069+2124211143820855593988*x))))))))))*(y)^(8)+2500000000*(138552077559380492720100083552259186+x*(-518129521994851713755497848220513949+100*x*(211060585084324952228717340164817+4*x*(181997514524528155615390486535897+25*x*(-191528781724603748936085265405+x*(-218178188813585401598670114969+80*x*(-316192049271703378831743693+10*x*(3242759559694185931373896+5*x*(43148931129595775808171+31046168331245185047725*x)))))))))*(y)^(9)+27500000000*(-4170899306537626970232096726120063+2*x*(-5594389308271109900804197339653193+25*x*(22937196753854554456622278215477+4*x*(2782364666669841709901242202438+25*x*(12655247052698824086276011299+8*x*(-362306615452479131928662592+5*x*(-21680317993242996375333143+180*x*(799650858117221408548+2385269629366586971095*x))))))))*(y)^(10)-300000000000*(32197381414060723715679384929442+x*(-74083180728190222571633627817181+10*x*(1787978184548109081564181844403+20*x*(35839990317274869045377525197+200*x*(-38454844495213594178907107+x*(-3614482475652951149560803+5*x*(-92854244493612083309023+50260732806860363055480*x)))))))*(y)^(11)-3250000000000*(-1386001602410594435486708400885+4*x*(-255050973614285344755027951184+5*x*(19773915409078463607981525441+20*x*(16772096711582432292766752+25*x*(-218069407031470220897891+12*x*(-18393983188864256510007+2078881834840572039740*x))))))*(y)^(12)+35000000000000*(9527124293161804308265821318+x*(-5150850242757708061721573957+40*x*(48233241498574572082579983+10*x*(-54691291943440779941534+25*x*(-15874909021313338081109+2545236255807717734007*x)))))*(y)^(13)+1125000000000000*(-37864715645896890114243841+2*x*(-4200264675491427861837169+5*x*(1072372985577804280726881+20*x*(-11343858139877762710678+726455486755366925275*x))))*(y)^(14)-4000000000000000*(1240395513980207205085266+x*(-341925624484588249279597+250*x*(-92745212797983746355+48508078875627116716*x)))*(y)^(15)-127500000000000000*(967519220929799627413-521836616853940347052*x+64740269323974620280*(x)^(2))*(y)^(16)+8100000000000000000*(203631613331080255+52205345473267412*x)*(y)^(17)+68707526255022096000000000000000000*(y)^(18))/1.69182035233245e48",2);
        f=gsFunctionExpr<real_t>("(243*(50497052605010560134063815435000252508733486047-199014893481648118666510199942247945744424463900*y)+20*(3*x*(212407118668070050400837154153816241340492165898+5*x*(-18417110037621892049763610258305244835011844133+20*x*(-680471534018737869090296316337848469115026083+80*x*(1175459346529485922139611944122104414210906+5*x*(172112960799268614739693592859142901011962+5*x*(-1716048929701797561985523204366578140681+25*x*(-61540159884695710531482095952281403942+x*(844705254259547966483851405161850047+100*x*(13311970415181511628653825075630350+x*(21645699412432041655030946513971+20*x*(-6840738135683690390772677187123+10*x*(-5518070219015182905821289939+100*x*(27808905458626020861295272+x*(366953170538786299641243+40*x*(-348306547122340607949+11916893274269801150*x)))))))))))))))-30*x*(33334846578718803712836745187308061439617260371+20*x*(-1235062268868930927860816357677686404097882117+100*x*(-4162363803420690873199447358185652490590699+20*x*(59354126418490270907082165869470140161808+x*(17547035779165263316468961476545420914483+5*x*(-494316158274430219657562901683749179639+10*x*(-13028169299228303533675304608094508296+25*x*(40844288345705546058646593265867851+10*x*(877548408929790707935605862947371+4*x*(-11672414521425367319178577561041+20*x*(-78110130302860872116052661922+125*x*(30690393923047724427250491+2*x*(795065881948686441543877+30*x*(-967170707027605271513+16218672343410725572*x))))))))))))))*y+180*(546784633884459797216264789039179116365722161+20*x*(-1261343022606605000579609325372221472989979+5*x*(-7910031470487984440969074309689579727435089+100*x*(249456105019200303851154822522557683342+5*x*(2658778978542364937119463883012840594164+x*(11837726480363995088877802606471616078+5*x*(-30236551928764998335683280960199343251+400*x*(-577270751437851340335378275133493+50*x*(34382298422256974073635525703362+x*(-157782044335771471639308494798+x*(-281100749404461107701887202019+200*x*(73673101806336391876342283+x*(2883637905851808887742981+25*x*(-18997479527234621613806+528402799876198962725*x))))))))))))))*(y)^(2)+100*(9017290723668017476849352046018796366923229347+250*x*(9492778262890475574690917992561811624648949+2*x*(-3672016492380678559820549609323939052626833+10*x*(-98747911651863524864358293730517276615679+10*x*(1883360543722071961245940708117682150673+2*x*(348889663344272393905935290004078564951+100*x*(-83055242212311260548420447421567781+10*x*(-9829995116001099856347546784538811+5*x*(-42441595642999529627470808931927+10*x*(1989955793891110454441595706607+30*x*(711566141441201205351729+350*x*(-731935742391105870558027+51939955730717659255126*x+739754669154418789140*(x)^(2)))))))))))))*(y)^(3)-3000*(21673824097841548582257881194967159573305707+100*x*(67307062718311535701271568890582266771212+5*x*(-27087802405204651505569485170902870948459+10*x*(-259091505734641385884563506842733972788+5*x*(57715317851430521102905951386608284001+40*x*(115217684517907150380276089395585652+5*x*(-10349458004456435064473344099802011+20*x*(-47464361636532608838961727858717+25*x*(291261872814751752245258403545+8*x*(2661511493441287143316932561+x*(-189918650537643461574473239+25*x*(-61934189881380958360838+30952138707475426534555*x))))))))))))*(y)^(4)+150000*(-755560363725838359148560498360055616719431+10*x*(-9398389697373280484814341607493095354511+8*x*(1415094952708344193780829970553508062152+5*x*(48606220060322909274751192041254372533+240*x*(-27043157640897968667989883659190327+5*x*(-2782368102220777091017371419980194+5*x*(-43768869279446291152036106405198+5*x*(2205506363056528264195298327218+5*x*(61168109026728820251894542829+70*x*(-16747992317543951268646941+40*x*(-96319517268418945191139+271324349632049063355*x)))))))))))*(y)^(5)+200000*(29104182198910131515974604531722603398751+5*x*(1597458383515281765931390171701584744226+5*x*(-665921640312003440225668625532789833643+400*x*(-73846786692495739528068979647246829+5*x*(24170024772592512246368424419143773+20*x*(96536989159256977948711411909977+10*x*(-2280584745509951477934817897207+2*x*(-223890845213618597042804274732+25*x*(66422983416824852434601577+5*x*(13567676776997201615030462+627502089657674564902863*x))))))))))*(y)^(6)-3000000*(-2183141472399162879850734636193543514253+10*x*(5203847932316006680924487269794620689+60*x*(425753650834482140037016531334563507+50*x*(398272255709993952277980984515093+10*x*(-11159307986682461224941277919733+2*x*(-2740157694693315875799124637527+10*x*(-24112942836407663765239735167+10*x*(208477385435870197417335587+25*x*(2707721917645645200751029+117992421195564066894710*x)))))))))*(y)^(7)+60000000*(-1858439803665643870582093506184009744+25*x*(-29764471284054843178141569434243102+5*x*(15236519769154931499868486299686379+20*x*(-52362908417205495508877028176011+10*x*(-3641558368086216508169976721159+12*x*(-4920755046711964831324968936+125*x*(26616532254093248644321669+2*x*(2536617248288587614014584+69553362429410636411685*x))))))))*(y)^(8)-100000000*(1814384361101580973443885929252954683+50*x*(-6787952704150918425822996742550169+20*x*(-155913858053357011016072814192261+20*x*(857024383932366856343012817511+180*x*(411909157530444376154960557+x*(146029709080024790265934293+125*x*(-71501768135188090199793+497222670838439389888*x)))))))*(y)^(9)-6000000000*(612080638310556845374925633967631+10*x*(-27057853148622611920984009081132+5*x*(5488885379474337745883239715903+40*x*(-29852035998534663201955877811+35*x*(-69102380756204146241290989+5*x*(2444066777750734913761658+126567460808479706709525*x))))))*(y)^(10)+30000000000*(66035767323965342858727737494787+10*x*(-2515527471708063637361944555387+30*x*(-7813199029223015950393815429+50*x*(90692135802509159620384297+30*x*(-191718947417303873969225+9917257402117231729682*x)))))*(y)^(11)+100000000000*(1068028978938837866922026856161+40*x*(-11388528454603134464531889237+35*x*(72287855979385413512289024+25*x*(-625337140047818638548434+42753880513474669329135*x))))*(y)^(12)+3000000000000*(-286768377678985797203887497+50*x*(9139797036319291657602643+120*x*(-19308984937116675326462+438275889262669832745*x)))*(y)^(13)+60000000000000*(29731229844053438771703+25*x*(9231283041521873999582+1049228045018022708195*x))*(y)^(14)+300000000000000*(11704072636745431716719+1160590434881456227170*x)*(y)^(15)+79296548109497025481500000000000000*(y)^(16)))/1.69182035233245e47",2); // 3
        iterLimit=5;
        break;
    case 7:
        ss << "5-patch domain, degree 4\n";
        pathStart = "Mario/Florian_FivePatches_deg4/multipatch_five";
        pathEnd = "non.xml";
        mapPathEnd = "_map.xml";
        g=gsFunctionExpr<real_t>("((51.2 - (32*x)/3. - y)^2*(-17.9 - (23*x)/5. - y)^2*(-6.855555555555555 - (13*x)/9. - y)^2*(3.7956521739130435 - (10*x)/23. - y)^2*(-4.41578947368421 - (2*x)/19. - y)^2*(3.0448275862068965 + x/29. - y)^(2)*(4.45 + x/2. - y)^(2)*(-4.923529411764706 + (9*x)/17. - y)^2*(6.1 + x - y)^(2)*(-21.6 + (9*x)/2. - y)^2)/1.e16",2); // 3
        dg=gsFunctionExpr<real_t>("((883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(5184*(-24+5*x)*(-60759448049004644133+10*x*(9224871288649613571+100*x*(37604179460108437+10*x*(-85805440864572+5*x*(-23585288349969+2*x*(-344361452181+40*x*(2225326239+88527705*x+747500*(x)^(2))))))))-135*(17023674641706872782653+20*x*(-1695509837688518024354+5*x*(-41960889498677262899+100*x*(284438575105026496+5*x*(2833672964522191+12*x*(-96263698034234+5*x*(-270006365587+80*x*(969070628+45871135*x))))))))*y-250*(3209323891457105431287+4*x*(-1063905721177541644621+20*x*(-10863701126609174361+5*x*(311839709514378744+25*x*(4409023124019543+4*x*(-39053380365891+50*x*(-47676560757+30484236416*x)))))))*(y)^(2)+500*(414137231579450182133+100*x*(-8503573356019086696+5*x*(-169890059042411955+4*x*(19655998498441548+25*x*(66279852581471+60*x*(-104658896724+22349641013*x))))))*(y)^(3)+25000*(669962652693164829+4*x*(-252981723406319757+20*x*(-3480355330614159+5*x*(-149651754464292+25*x*(1579945392551+246871418796*x)))))*(y)^(4)-50000*(36497337774982877+20*x*(-13208323166449558+15*x*(-26473933807981+40*x*(345529412146+198697687975*x))))*(y)^(5)-500000*(-1618566060310711+60*x*(-11342941558297+40*x*(-39209104117+12314778830*x)))*(y)^(6)+15000000*(-8007076526467+20*x*(-437203258412+139024822065*x))*(y)^(7)+50000000*(-475397417077+52105762340*x)*(y)^(8)-621668112000000000*(y)^(9)))/1.69182035233245e49",
                                   "-((883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(9*(-24+5*x)*(-16077059186642225498133+4*x*(-1369336665190955685819+5*x*(48913670472575889613+20*x*(946609999700553672+5*x*(-5001444039303987+100*x*(-28130131110147+x*(156037145859+40*x*(2620771809+91742270*x))))))))+20*(145441291720013932895181+5*x*(3209323891457105431287+2*x*(-1063905721177541644621+40*x*(-3621233708869724787+5*x*(77959927378594686+5*x*(4409023124019543+10*x*(-13017793455297+100*x*(-6810937251+3810529552*x))))))))*y-300*(4238909567658291807147+x*(414137231579450182133+100*x*(-4251786678009543348+5*x*(-56630019680803985+4*x*(4913999624610387+5*x*(66279852581471+300*x*(-17443149454+3192805859*x)))))))*(y)^(2)-4000*(75855763148894346513+5*x*(669962652693164829+2*x*(-252981723406319757+40*x*(-1160118443538053+5*x*(-37412938616073+5*x*(1579945392551+205726182330*x))))))*(y)^(3)+50000*(1644683832982118973+x*(36497337774982877+20*x*(-6604161583224779+5*x*(-26473933807981+60*x*(172764706073+79479075190*x)))))*(y)^(4)+600000*(18992823414703953+x*(-1618566060310711+10*x*(-34028824674891+40*x*(-78418208234+18472168245*x))))*(y)^(5)-21000000*(64762860701939+x*(-8007076526467+20*x*(-218601629206+46341607355*x)))*(y)^(6)-80000000*(1985626097571+x*(-475397417077+26052881170*x))*(y)^(7)+2700000000*(-500083579+414445408*x)*(y)^(8)+116338140000000000*(y)^(9)))/3.3836407046649e48",2);
        ddg=gsFunctionExpr<real_t>("(80*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(9*(-543914379+1936272910*y)+100*(162*x*(-1003582+5*x*(-993+3568*x+50*(x)^(2)))-18*x*(-3185437+643605*x+61400*(x)^(2))*y+(-59423569+20*x*(601672+190725*x))*(y)^(2)-50*(25691+84020*x)*(y)^(3)+1294400*(y)^(4)))*(2619201987768-935774142455*y+50*(16*x*(-2700305626+x*(-469551171+80*x*(1108018+37375*x)))+x*(-11430373379+20*x*(132519969+24500680*x))*y+4*(-1888260208+5*x*(102280583+48691290*x))*(y)^(2)+10*(2825021+47307530*x)*(y)^(3)+50886600*(y)^(4)))+25*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*((883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(837-90*x+170*y)^(2)*(81*(-5687+60*x*(13+10*x))+180*(167-330*x)*y+15700*(y)^(2))+8*(883+10*x-290*y)*(-837+90*x-170*y)*(117+180*x-110*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(9*(-11809+60*x*(293+5*x))-200*(1485+296*x)*y+104900*(y)^(2))+4*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(27*(-38127151491+30252883840*y)+40*(243*x*(-8844351+5*x*(749023+50*x*(586+5*x)))-360*x*(7615641+10*x*(217413+3700*x))*y+10*(-75702173+60*x*(4329504+266245*x))*(y)^(2)-300*(7743339+1596890*x)*(y)^(3)+421028250*(y)^(4))))+4*(883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(837-90*x+170*y)^(2)*(20*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(-3344+7360*x+1145*y)*(-5422463+240*x*(8226+325*x)+2406380*y+649600*x*y+780500*(y)^(2))+(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*(256*(-1932479+1380*x*(-209+230*x))+320*(-81373+79005*x)*y+1863025*(y)^(2))+25*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(-59816709411407-65166891181640*y+40*(6*x*(-163190328849+20*x*(865547273+6500*x*(16452+325*x)))+20*x*(3646882989+40*x*(143693421+5278000*x))*y+5*(-325804009+360*x*(135707027+12173420*x))*(y)^(2)+300*(536897933+131866100*x)*(y)^(3)+21616698250*(y)^(4)))))/1.69182035233245e49",
                                    "(-40*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(8106539481+9*x*(-193627291+10*x*(-3185437+10*x*(42907+3070*x)))+171788588*y-20*x*(-59423569+20*x*(300836+63575*x))*y+60*(-19223711+25*x*(25691+42010*x))*(y)^(2)-800*(-127477+64720*x)*(y)^(3)+9860000*(y)^(4))*(1804179436833-1545589725720*y+4*(x*(-187154828491+5*x*(-11430373379+40*x*(44173323+6125170*x)))+40*x*(-3776520416+5*x*(102280583+32460860*x))*y+300*(-347459898+x*(2825021+23653765*x))*(y)^(2)+6000*(-507835+339244*x)*(y)^(3)+147487500*(y)^(4)))+4*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*((883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(837-90*x+170*y)^(2)*(-2327+3925*(x)^(2)+300*y*(31+2*y)-20*x*(794+165*y))+400*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(-31+11*x-4*y)*(-115123+1480*(x)^(2)+x*(14850-10490*y)+3*y*(-8451+4930*y))+2*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(983008657193+3696921147420*y+20*(26624500*(x)^(4)+600*(x)^(3)*(1443168-798445*y)+15*(y)^(2)*(-3897676511+98600*y*(-8451+2465*y))+10*(x)^(2)*(-75702173+90*y*(-7743339+2806855*y))+x*(-70039005912+10*y*(2650634259+10*(167638491-51715700*y)*y)))))+(883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(837-90*x+170*y)^(2)*(8*(179+46*x+10*y)*(-768+160*x+15*y)*(-999+458*x+60*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(-4919953+324800*(x)^(2)+20*y*(294019+58995*y)+20*x*(120319+78050*y))+(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*(298084*(x)^(2)+60*x*(-16589+1374*y)+9*(-72407+60*y*(-333+10*y)))+4*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(-241540077242153-193437613834440*y+40*(3652026000*(x)^(4)+600*(x)^(3)*(135707027+65933050*y)+15*(y)^(2)*(47746822063+65550*y*(588038+58995*y))+5*(x)^(2)*(-325804009+60*y*(1610693799+432333965*y))+x*(-3152734415003+30*y*(5512681989+10*y*(3241247549+511617750*y)))))))/6.7672814093298e47",
                                    "(-11664*(24-5*x)^(2)*(3795508851241428081609835296757521036403467+8*x*(-617924889466712236003969418712490006166196+5*x*(-117799961390115888244607713916466201335763+10*x*(-1153754925035254244930228232162930630276+5*x*(190783077681420820239961284719819672741+400*x*(94136542817661227100453777189555504+5*x*(-2335776705045331973430661489837887+5*x*(-139662379522414132011851610140952+x*(2989945885086361921720201255521+40*x*(55027619254438667298341299236+5*x*(193110516368558089453986589+10*x*(-7645478313504331939059012+5*x*(-65763204980290218177741+160*x*(19183304824653527736+5*x*(350631176011193151+40*x*(197697136981902+1302969589675*x))))))))))))))))-3645*(-24+5*x)*(-133217222313490691110038761522914324054767573+4*x*(30707381976125011223697233902455332398970319+20*x*(1191165278871544872251481597207897295736307+25*x*(6059617791835945705996637645576510043311+5*x*(-587254050168334401901216944803729861137+80*x*(-2857926878484637578930592938424051364+5*x*(-1924196692092154206887727370076932+x*(24731098175816427171948650852354142+5*x*(298808202379124496086373676196557+20*x*(-4061264643115522705627614150193+20*x*(-9322530128165611705629893263+5*x*(254749554884573818720518079+5*x*(-2504987654656513967213119+40*x*(-3735831775608908438463+80*x*(11411323124051659378+15*x*(23233173342897719+23271419454370*x))))))))))))))))*y+675*(25695529994036350691195886821931982550755068915+2*x*(-21222561349774070294281120924390666690935287463+5*x*(-1470192296359704686899216335439906402872612009+20*x*(32828438928021100531095137055388422092801442+25*x*(375771878930388547769390107169649955359629+8*x*(-9297559430780970006311141091744818659836+5*x*(-518824950630865729493051643271960121031+40*x*(1239337523806790172465481593355752474+5*x*(70298941463740785448635833708320698+25*x*(-165266538000479350228250648119671+x*(-34434605624479501963715560448477+20*x*(127192476002036628020827168206+5*x*(1247946149098604414584565817+200*x*(-1222897611049032857462391+10*x*(1161178534000944312057+320*x*(1266531504681501663+41525418232972513*x))))))))))))))))*(y)^(2)+1000*(454908891261705405293290728846770879955211962+x*(1192377384387916217783146947465730421600899239+10*x*(4323767282325281414116112562531325444590531+20*x*(-3432367826491467794474142636672332254331339+200*x*(-829274239270937985012716379428979999543+x*(1866368217671403190852503064272450591753+5*x*(14660694531649725548851755369383919711+20*x*(-814066796197809791333767300226627538+25*x*(-845706731082957996585568457594364+5*x*(136341428904494516644467134448151+2*x*(-523686404450488294038346109853+100*x*(-6401149892629696416060739017+20*x*(17276974088396173674640689+x*(893137738722155892236061+100*x*(-308595463037841925575+35323694640674483104*x)))))))))))))))*(y)^(3)-12500*(203547927545787429973647834403689761163930297+4*x*(-88700405875258812633570461844349224850692594+5*x*(-5088884454361057926285129444449369229933843+20*x*(102037029733931730830901970205868315022536+25*x*(1139056577223180365690095879519796153237+20*x*(-5997254477044199421896067464297433513+4*x*(-675585137903823909090519779730570562+25*x*(-395141985092314217490706644985944+x*(588273024309087438974071601399403+100*x*(304204911714318416453141128308+x*(-62276343587707471314832158463+20*x*(-20909933786942601669098856+5*x*(2866213319543382630659383+60*x*(512820645506886052273+550407119691908322750*x))))))))))))))*(y)^(4)-150000*(-820912876560040055454432453797226841034706+5*x*(592450044345189421336013296999840871075323+80*x*(1311506726167945985871175023744316137447+x*(-1803079355393752025371711836650070599656+25*x*(-7267115612631854965848892570016054753+x*(4807123355570829949579856201052443019+80*x*(6152635135507867630226220089665361+5*x*(-309672464598060704832373938469532+5*x*(-7725830568381436150542669393672+25*x*(30980980740115059043000786097+16*x*(172350682542795431883899746+25*x*(-796881288663087402354876+5*x*(7391466654891332268515+236379089094056247129*x)))))))))))))*(y)^(5)+1750000*(80773870313750739897812091863957297184759+2*x*(-92643044110063114738520207042823764986917+25*x*(-719117177967701738403704851801130970015+4*x*(77634810668446222272798278038044225162+25*x*(766259340515986553579744402873911759+16*x*(-1214932437321533156667985363337886+125*x*(-15564052439637071236626757180113+2*x*(-844742087491991094270169700424+x*(80097624826001896502857819221+50*x*(405967284015213374071442325+x*(-10201816907430873985157249+60*x*(12140065785830258340462+8268842088053710100275*x))))))))))))*(y)^(6)-20000000*(539426160540983772054245755685428914630+x*(-1790941405357047567376315431779649452431+10*x*(-18091672375436750457831888150781079637+20*x*(1591626137629073876880994869869520241+400*x*(263805077291529704155031159940002+x*(-168796794272635732890337453811646+5*x*(-4533024644089582279726522153397+50*x*(6604069465120570812393238904+x*(2427734883696585912749793852+25*x*(5319031552190841163885487-970751473560998415598774*x+54273906819203417628840*(x)^(2)))))))))))*(y)^(7)-225000000*(8279588480765903113524666213672961341+4*x*(-14802325646602365620237113770635840653+40*x*(-24629286857699112254909947458500148+5*x*(4750519337020804929328849712882428+25*x*(34961665971044349794100633258143+4*x*(-97995562129368317132554465683+10*x*(-34798545644587604875845990133+15*x*(-276055521378564648930042248+5*x*(-3817133550493210740612501+100*x*(13829622398451154550069+2124211143820855593988*x))))))))))*(y)^(8)+2500000000*(138552077559380492720100083552259186+x*(-518129521994851713755497848220513949+100*x*(211060585084324952228717340164817+4*x*(181997514524528155615390486535897+25*x*(-191528781724603748936085265405+x*(-218178188813585401598670114969+80*x*(-316192049271703378831743693+10*x*(3242759559694185931373896+5*x*(43148931129595775808171+31046168331245185047725*x)))))))))*(y)^(9)+27500000000*(-4170899306537626970232096726120063+2*x*(-5594389308271109900804197339653193+25*x*(22937196753854554456622278215477+4*x*(2782364666669841709901242202438+25*x*(12655247052698824086276011299+8*x*(-362306615452479131928662592+5*x*(-21680317993242996375333143+180*x*(799650858117221408548+2385269629366586971095*x))))))))*(y)^(10)-300000000000*(32197381414060723715679384929442+x*(-74083180728190222571633627817181+10*x*(1787978184548109081564181844403+20*x*(35839990317274869045377525197+200*x*(-38454844495213594178907107+x*(-3614482475652951149560803+5*x*(-92854244493612083309023+50260732806860363055480*x)))))))*(y)^(11)-3250000000000*(-1386001602410594435486708400885+4*x*(-255050973614285344755027951184+5*x*(19773915409078463607981525441+20*x*(16772096711582432292766752+25*x*(-218069407031470220897891+12*x*(-18393983188864256510007+2078881834840572039740*x))))))*(y)^(12)+35000000000000*(9527124293161804308265821318+x*(-5150850242757708061721573957+40*x*(48233241498574572082579983+10*x*(-54691291943440779941534+25*x*(-15874909021313338081109+2545236255807717734007*x)))))*(y)^(13)+1125000000000000*(-37864715645896890114243841+2*x*(-4200264675491427861837169+5*x*(1072372985577804280726881+20*x*(-11343858139877762710678+726455486755366925275*x))))*(y)^(14)-4000000000000000*(1240395513980207205085266+x*(-341925624484588249279597+250*x*(-92745212797983746355+48508078875627116716*x)))*(y)^(15)-127500000000000000*(967519220929799627413-521836616853940347052*x+64740269323974620280*(x)^(2))*(y)^(16)+8100000000000000000*(203631613331080255+52205345473267412*x)*(y)^(17)+68707526255022096000000000000000000*(y)^(18))/1.69182035233245e48",2);
        f=gsFunctionExpr<real_t>("(243*(50497052605010560134063815435000252508733486047-199014893481648118666510199942247945744424463900*y)+20*(3*x*(212407118668070050400837154153816241340492165898+5*x*(-18417110037621892049763610258305244835011844133+20*x*(-680471534018737869090296316337848469115026083+80*x*(1175459346529485922139611944122104414210906+5*x*(172112960799268614739693592859142901011962+5*x*(-1716048929701797561985523204366578140681+25*x*(-61540159884695710531482095952281403942+x*(844705254259547966483851405161850047+100*x*(13311970415181511628653825075630350+x*(21645699412432041655030946513971+20*x*(-6840738135683690390772677187123+10*x*(-5518070219015182905821289939+100*x*(27808905458626020861295272+x*(366953170538786299641243+40*x*(-348306547122340607949+11916893274269801150*x)))))))))))))))-30*x*(33334846578718803712836745187308061439617260371+20*x*(-1235062268868930927860816357677686404097882117+100*x*(-4162363803420690873199447358185652490590699+20*x*(59354126418490270907082165869470140161808+x*(17547035779165263316468961476545420914483+5*x*(-494316158274430219657562901683749179639+10*x*(-13028169299228303533675304608094508296+25*x*(40844288345705546058646593265867851+10*x*(877548408929790707935605862947371+4*x*(-11672414521425367319178577561041+20*x*(-78110130302860872116052661922+125*x*(30690393923047724427250491+2*x*(795065881948686441543877+30*x*(-967170707027605271513+16218672343410725572*x))))))))))))))*y+180*(546784633884459797216264789039179116365722161+20*x*(-1261343022606605000579609325372221472989979+5*x*(-7910031470487984440969074309689579727435089+100*x*(249456105019200303851154822522557683342+5*x*(2658778978542364937119463883012840594164+x*(11837726480363995088877802606471616078+5*x*(-30236551928764998335683280960199343251+400*x*(-577270751437851340335378275133493+50*x*(34382298422256974073635525703362+x*(-157782044335771471639308494798+x*(-281100749404461107701887202019+200*x*(73673101806336391876342283+x*(2883637905851808887742981+25*x*(-18997479527234621613806+528402799876198962725*x))))))))))))))*(y)^(2)+100*(9017290723668017476849352046018796366923229347+250*x*(9492778262890475574690917992561811624648949+2*x*(-3672016492380678559820549609323939052626833+10*x*(-98747911651863524864358293730517276615679+10*x*(1883360543722071961245940708117682150673+2*x*(348889663344272393905935290004078564951+100*x*(-83055242212311260548420447421567781+10*x*(-9829995116001099856347546784538811+5*x*(-42441595642999529627470808931927+10*x*(1989955793891110454441595706607+30*x*(711566141441201205351729+350*x*(-731935742391105870558027+51939955730717659255126*x+739754669154418789140*(x)^(2)))))))))))))*(y)^(3)-3000*(21673824097841548582257881194967159573305707+100*x*(67307062718311535701271568890582266771212+5*x*(-27087802405204651505569485170902870948459+10*x*(-259091505734641385884563506842733972788+5*x*(57715317851430521102905951386608284001+40*x*(115217684517907150380276089395585652+5*x*(-10349458004456435064473344099802011+20*x*(-47464361636532608838961727858717+25*x*(291261872814751752245258403545+8*x*(2661511493441287143316932561+x*(-189918650537643461574473239+25*x*(-61934189881380958360838+30952138707475426534555*x))))))))))))*(y)^(4)+150000*(-755560363725838359148560498360055616719431+10*x*(-9398389697373280484814341607493095354511+8*x*(1415094952708344193780829970553508062152+5*x*(48606220060322909274751192041254372533+240*x*(-27043157640897968667989883659190327+5*x*(-2782368102220777091017371419980194+5*x*(-43768869279446291152036106405198+5*x*(2205506363056528264195298327218+5*x*(61168109026728820251894542829+70*x*(-16747992317543951268646941+40*x*(-96319517268418945191139+271324349632049063355*x)))))))))))*(y)^(5)+200000*(29104182198910131515974604531722603398751+5*x*(1597458383515281765931390171701584744226+5*x*(-665921640312003440225668625532789833643+400*x*(-73846786692495739528068979647246829+5*x*(24170024772592512246368424419143773+20*x*(96536989159256977948711411909977+10*x*(-2280584745509951477934817897207+2*x*(-223890845213618597042804274732+25*x*(66422983416824852434601577+5*x*(13567676776997201615030462+627502089657674564902863*x))))))))))*(y)^(6)-3000000*(-2183141472399162879850734636193543514253+10*x*(5203847932316006680924487269794620689+60*x*(425753650834482140037016531334563507+50*x*(398272255709993952277980984515093+10*x*(-11159307986682461224941277919733+2*x*(-2740157694693315875799124637527+10*x*(-24112942836407663765239735167+10*x*(208477385435870197417335587+25*x*(2707721917645645200751029+117992421195564066894710*x)))))))))*(y)^(7)+60000000*(-1858439803665643870582093506184009744+25*x*(-29764471284054843178141569434243102+5*x*(15236519769154931499868486299686379+20*x*(-52362908417205495508877028176011+10*x*(-3641558368086216508169976721159+12*x*(-4920755046711964831324968936+125*x*(26616532254093248644321669+2*x*(2536617248288587614014584+69553362429410636411685*x))))))))*(y)^(8)-100000000*(1814384361101580973443885929252954683+50*x*(-6787952704150918425822996742550169+20*x*(-155913858053357011016072814192261+20*x*(857024383932366856343012817511+180*x*(411909157530444376154960557+x*(146029709080024790265934293+125*x*(-71501768135188090199793+497222670838439389888*x)))))))*(y)^(9)-6000000000*(612080638310556845374925633967631+10*x*(-27057853148622611920984009081132+5*x*(5488885379474337745883239715903+40*x*(-29852035998534663201955877811+35*x*(-69102380756204146241290989+5*x*(2444066777750734913761658+126567460808479706709525*x))))))*(y)^(10)+30000000000*(66035767323965342858727737494787+10*x*(-2515527471708063637361944555387+30*x*(-7813199029223015950393815429+50*x*(90692135802509159620384297+30*x*(-191718947417303873969225+9917257402117231729682*x)))))*(y)^(11)+100000000000*(1068028978938837866922026856161+40*x*(-11388528454603134464531889237+35*x*(72287855979385413512289024+25*x*(-625337140047818638548434+42753880513474669329135*x))))*(y)^(12)+3000000000000*(-286768377678985797203887497+50*x*(9139797036319291657602643+120*x*(-19308984937116675326462+438275889262669832745*x)))*(y)^(13)+60000000000000*(29731229844053438771703+25*x*(9231283041521873999582+1049228045018022708195*x))*(y)^(14)+300000000000000*(11704072636745431716719+1160590434881456227170*x)*(y)^(15)+79296548109497025481500000000000000*(y)^(16)))/1.69182035233245e47",2); // 3
        iterLimit=5;
        break;
    case 8:
        ss << "3-patch domain, degree 3, v2\n";
        pathStart = "Mario/Florian_ThreePatches_deg3_v2/multipatch_three";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        startIter=1;
        g=gsFunctionExpr<real_t>("((-47.7 - 14*x - y)^(2)*(-3.8774193548387097 - (16*x)/31. - y)^2*(3.875 - x/2. - y)^(2)*(-3.7240259740259742 + (39*x)/77. - y)^2*(3.8243055555555556 + (37*x)/72. - y)^2*(-286.75 + 77*x - y)^(2))/2.5e11",2); // 2
        dg=gsFunctionExpr<real_t>("((5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-441389437313011+5272968224115*y+8*(298666368000*(x)^(5)+28000*(x)^(4)*(-4133821+672325*y)-480*(x)^(3)*(52101274671+5*y*(-48329457+640401050*y))+12*(x)^(2)*(565843921171+5*y*(14511590237+20*(46516301-64021780*y)*y))+(y)^(2)*(8076277817723+20*y*(-58621313017+200*y*(-59499360+10828163*y)))+2*x*(229800161083395+y*(2082084269689+10*y*(-1675159861883+20*y*(207397103+3706248990*y)))))))/1.5123064061952e31",
                                   "((5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(31*(-65076108692243+567148797931020*y)+4*(x*(5272968224115+8*x*(2082084269689+20*x*(14511590237+10*x*(144988371+18825100*x))))-16*x*(-8076277817723+10*x*(1675159861883+680*x*(-2736253+56505975*x)))*y+12*(15296479392949-40*x*(58621313017+10*x*(-207397103+128043560*x)))*(y)^(2)+480*(-624025745233+400*x*(-39666240+123541633*x))*(y)^(3)+8000*(-1149261923+216563260*x)*(y)^(4)-32997888000*(y)^(5))))/6.0492256247808e31",2);
        ddg=gsFunctionExpr<real_t>("(2*(5507+740*x-1440*y)^(2)*(1147-308*x+4*y)^(2)*(1147-156*x+308*y)^(2)*(139686329327-200681825360*y+25*(16*x*(-362657382+x*(-101620989+28000*x*(355+84*x)))+8*x*(-373082571+20*x*(4913781+1257200*x))*y+(-653061619+240*x*(5025693+1529365*x))*(y)^(2)+110*(4996423+2317980*x)*(y)^(3)+48905650*(y)^(4)))+8*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-489329699+26666640*(x)^(2)+8*y*(15792077+4325450*y)-8*x*(8077099+8740640*y))*(-330574+75085*y+25*(672*(x)^(2)+931*(y)^(2)+4*x*(355+449*y)))+(-31+4*x+8*y)^(2)*(477+140*x+10*y)^(2)*(601+80*x+155*y)^(2)*(-3441*(-35572636019847+73096675692352*y)+16*(74073925926000*(x)^(4)-44444400*(x)^(3)*(8077099+8740640*y)+4*(y)^(2)*(-172228484745029+20*y*(13798425242429+959704939205*y))+6*(x)^(2)*(-478459593600089+80*y*(3519648832631+1435589920070*y))+3*x*(3988696521770833-8*y*(-206187448112451+10*y*(10414257900957+1896520058240*y))))))/3.780766015488e30",
                                    "(-32*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-184879+179600*(x)^(2)+60*y*(19717+620*y)+40*x*(15017+9310*y))*(-9535011+17481280*(x)^(2)+4*y*(31783721+166320*y)-4*x*(15792077+8650900*y))+1573679923200*(-31+4*x+8*y)^(2)*(477+140*x+10*y)^(2)*(601+80*x+155*y)^(2)*(((11119+113140*x-221760*y)*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y))/(2.45887488e10)+2*(-3.7240259740259742+(39*x)/77.-y)^2*(3.8243055555555556+(37*x)/72.-y)^2+2*(286.75-77*x+y)^(2)*((-350064735599+200*x*(15438523+96003266*x))/(1.22943744e10)-((11119+113140*x)*y)/18480.+6*(y)^(2)))+(5507+740*x-1440*y)^(2)*(1147-308*x+4*y)^(2)*(1147-156*x+308*y)^(2)*(-10479303116099-1317338496540*y+40*(4*x*(-41626408888+5*x*(-653061619+40*x*(10051386+1529365*x)))+330*x*(130950371+20*x*(4996423+772660*x))*y+15*(3452990809+40*x*(58793797+9781130*x))*(y)^(2)+62000*(59151+18620*x)*(y)^(3)+57660000*(y)^(4))))/6.0492256247808e31",
                                    "(263879111196180480000000*(x)^(10)-200704000000*(x)^(9)*(-8649747410767+243871656618475*y)-1102267*(-249908537057869493753722+7303900366010386651975527*y)-38707200000*(x)^(8)*(290585708203906+3*y*(-197029885967433+80994558247885*y))+8192000*(x)^(7)*(-9481688358976200867+5*y*(78477855207928463297+50*y*(-12621679664055075+221243488766886308*y)))+2867200*(x)^(6)*(-280977579189855473287+5*y*(-22493549631768585111+5*y*(-2747278652088440511+20*y*(-91855653718880188+46208541740778675*y))))-92160*(x)^(5)*(75483626813060289160351+25*y*(-20046466182105866760727+100*y*(11594160804861041637+20*y*(-887558443306536565+4*y*(1720896663954649+71090383325622549*y)))))-7680*(x)^(4)*(-4964239415059682722723998+5*y*(2296897469493749744928169+25*y*(4891483627964772371247+40*y*(5266351819578766236+5*y*(-13129448820140085149+80*y*(-26146153193635401+10391360459106683*y))))))+1280*(x)^(3)*(319898962438204378868630125+y*(-3820925583992735261392046391+10*y*(-4783100735284830070720533+20*y*(2105413765947766789130671+40*y*(-281075782689312175923+5*y*(-527029110724154050557+8*y*(1102027571944353971+2738805988130334580*y)))))))+(y)^(2)*(617841001190024474584165568361-80*y*(-21074739527441629699914866180+y*(2607549751484824916038727667+16*y*(91673462096616111761339829+125*y*(-107496658047501834366731+80*y*(-225835710267402201780+y*(31492159069935608973+80*y*(12239879471796649+40941370928304*y))))))))+288*(x)^(2)*(-920093625250723826423869547+10*y*(892691419041710198843985595+y*(-47417212378914652851450213+20*y*(-6547575133692972562670908+5*y*(131944367478582199929979+32*y*(2394158422480936839957+5*y*(-79873367397145988827+200*y*(-56248512753087052+12046162278210363*y))))))))+4*x*(-1197623698331623912863154084953+y*(12321028396612468035229975389749+40*y*(16186340135808633803985384165+4*y*(-13097715815752721726226123093+20*y*(-24778019926500117047201121+y*(43068787960711552962518799+160*y*(8538946273367625114111+10*y*(-695886202065453168197+5*y*(-4066084107601811649+10710407262566600*y))))))))))/7.561532030976e30",2);
        f=gsFunctionExpr<real_t>("(775*(10070636445024446987294496713 + 1580805041435296880732831616*y) + 32*(839212170045134580480000*(x)^(8) + 3072000*(x)^(7)*(-186153527560215013 + 23630822384505615*y) - 12800*(x)^(6)*(7397845103057974754479 + 600*y*(-1978168936115700 + 521427067157777179*y)) - 23040*(x)^(5)*(-2591139538241579069079 + 5*y*(-19502249110340706949 + 40*y*(-532005332078551427 + 231162700790668990*y))) + 240*(x)^(4)*(20629697508869638693714899 + 80*y*(-2832250234210447845948 + 5*y*(361036991177109628031 + 80*y*(851520842693099156 + 1397012256516649535*y)))) + 480*(x)^(3)*(-5011292996454793182398419 + y*(172125474399395900332421 + 80*y*(6788570295304368834777 + 10*y*(-121146247699783170709 + 40*y*(-983617023935690419 + 258889960348724901*y))))) + (y)^(2)*(-113362049203585367125403749833 + 80*y*(-133919062100907890758222638 + y*(119712510000704219793039297 + 16*y*(506422009899324159806439 + 5*y*(-20050351692839232072401 + 120*y*(-4976219400498134306 + 701059526547943945*y)))))) + 6*(x)^(2)*(-13034991506843954692628435037 + 40*y*(-12911412490822047011077716 + y*(71942117064863810449017711 + 800*y*(3016805471122684680351 + y*(-7306648318334971703039 + 8*y*(-22005304021302476127 + 13952322443747600185*y)))))) + 6*x*(3583617805503366573901766935 + y*(-539997863237214241774925283 + 40*y*(-32876604623363072658989403 + 2*y*(3430167819779391335035679 + 240*y*(5681012924106078824787 + y*(-1148410005247339716587 + 200*y*(-521520979887370115 + 80705670151250134*y)))))))))/1.5123064061952e29",2); // 2
        iterLimit=5;
        break;
    case 9:
        ss << "4-patch domain, degree 3, v2\n";
        pathStart = "Mario/Florian_FourPatches_deg3_v2/multipatch_four";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        startIter=1;
        g=gsFunctionExpr<real_t>("((856 + 25*x - 210*y)^(2)*(794 + 185*x - 15*y)^(2)*(253 - 60*x + 5*y)^(2)*(819 + 190*x + 20*y)^(2)*(-772 + 185*x + 25*y)^(2)*(343 - 10*x + 80*y)^(2)*(507 + 10*x + 120*y)^(2)*(-1447 + 60*x + 350*y)^(2))/7.001316e42",2);
        dg=gsFunctionExpr<real_t>("((856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-183266112727112768808+93639600000000*(x)^(7)+388500000*(x)^(6)*(5638646+318315*y)-150000*(x)^(5)*(1137334053991+5*y*(38770857663+14530646650*y))+12500*(x)^(4)*(-152716930485143+25*y*(-156401552027+349513306994*y+7224284480*(y)^(2)))+1000*(x)^(3)*(70984719289425039+5*y*(-118089672499773+50*y*(-29650643938913+293423504802*y+882217595180*(y)^(2))))+75*(x)^(2)*(677492733286137563+5*y*(11600640434709637+10*y*(-981143799014747+10*y*(-12665036490069+50*y*(20291359787+6678463216*y)))))-5*y*(4654806746679084698+5*y*(-595081862571168173+5*y*(-21886279724378373+10*y*(156538493064151+160*y*(377450354536+45*y*(156065933+4364500*y))))))-20*x*(61409432176835309771+5*y*(-231054258707148997+5*y*(-273537297379319971+5*y*(686445838949453+10*y*(155942373681973+100*y*(-5142492283+817279476*y))))))))/7.001316e41",
                                   "((856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-24*(5557475750222879756+52643682537326664905*y)+5*(3533296500000*(x)^(7)-25000*(x)^(6)*(38770857663+29061293300*y)+12500*(x)^(5)*(-156401552027+4*y*(174756653497+5418213360*y))+250*(x)^(4)*(-118089672499773+100*y*(-29650643938913+y*(440135257203+1764435190360*y)))+25*(x)^(3)*(11600640434709637+20*y*(-981143799014747+5*y*(-37995109470207+200*y*(20291359787+8348079020*y))))+30*(y)^(2)*(199526503227539871+50*y*(9930626942414458+5*y*(-17989001088697+24*y*(-141086396887+122500*y*(2701+384*y)))))-10*(x)^(2)*(-231054258707148997+5*y*(-547074594758639942+5*y*(2059337516848359+40*y*(155942373681973+25*y*(-25712461415+4903676856*y)))))-x*(4654806746679084698+5*y*(-1190163725142336346+5*y*(-65658839173135119+40*y*(156538493064151+400*y*(188725177268+9*y*(468197799+15275750*y)))))))))/7.001316e41",2);
        ddg=gsFunctionExpr<real_t>("((819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(4*(-3+100*x-820*y)*(167+4440*x-365*y)*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)+(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(-5872151+15000*(x)^(2)-300*x*(3+820*y)+80*y*(949+12605*y))+(856+25*x-210*y)^(2)*(343-10*x+80*y)^(2)*(-178355327+29570400*(x)^(2)-13320*x*(-167+365*y)+5*y*(-55462+39965*y)))+(856+25*x-210*y)^(2)*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(343-10*x+80*y)^(2)*(40*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(319+24*x+214*y)*(-1447+60*x+350*y)*(967+14060*x+1690*y)+(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(-1777002527+296525400*(x)^(2)+42180*x*(967+1690*y)+60*y*(290448+71035*y))+4*(819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(-6259523+21600*(x)^(2)+1800*x*(319+214*y)+20*y*(172951+82445*y)))+20*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-9685978+444000*(x)^(3)+(x)^(2)*(5070-5515950*y)+y*(54496707+25*(13225-49548*y)*y)+2*x*(-132371273+5*y*(287721+1551850*y)))*(-2726354163+16872000*(x)^(3)+300*(x)^(2)*(1127087+762350*y)+10*y*(-257314061+40*y*(148557+180125*y))+40*x*(-260893019+5*y*(646909+3134030*y))))/1.4002632e41",
                                    "((819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(200*(-71+82*x-672*y)*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(7+73*x-6*y)*(-253+60*x-5*y)+(856+25*x-210*y)^(2)*(343-10*x+80*y)^(2)*(-1204067+199825*(x)^(2)+x*(30560-32850*y)+450*y*(-7+3*y))+4*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(-98526263+252100*(x)^(2)+50400*y*(71+336*y)-20*x*(14807+206640*y)))+(856+25*x-210*y)^(2)*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(343-10*x+80*y)^(2)*(8*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(1007+1690*x+200*y)*(-1447+60*x+350*y)*(381+1070*x+8400*y)+(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(-24276671+4262100*(x)^(2)+600*y*(1007+100*y)+60*x*(59951+16900*y))+4*(819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(1648900*(x)^(2)+180*x*(78963+149800*y)+27*(-22818637+2800*y*(127+1400*y))))-20*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(x*(-54496707+5*x*(-287721+367730*x))-250*x*(2645+62074*x)*y+450*(713+8258*x)*(y)^(2)-201600*(y)^(3)+4*(6617497+67949189*y))*(-3*(406850873+3587051620*y)+10*(7623500*(x)^(3)+7500*(y)^(2)*(1707+224*y)+10*(x)^(2)*(646909+6268060*y)+x*(-257314061+120*y*(99038+180125*y)))))/1.4002632e41",
                                    "(8736*(-5889591170879306018956271073243080922+22925776618634104814916613389531613295*y)+25*(62035588301512500000000000*(x)^(14)-4312350000000000*(x)^(13)*(3324258357687+2740557975655*y)-1803750000000*(x)^(12)*(307513783018520481+110*y*(1043969817180299+105368781544875*y))+18750000000*(x)^(11)*(1616209809586420291851+y*(1238021454387970462733+2*y*(183253666887968429619+58954243271100176660*y)))+137500000*(x)^(10)*(3860051208614600588030816+5*y*(-305396154884365519380817+50*y*(-2492215107508311154383+20*y*(-30895915857524455283+324570185189587235*y))))-6250000*(x)^(9)*(2469981477768200761696850107+5*y*(-170548179992646678377359407+10*y*(-15860356643192086745823063+20*y*(-245996004761701292488629+100*y*(462499575852494684963+191973545558897389557*y)))))+5625000*(x)^(8)*(-7933898947401529305690216499+10*y*(7139097275043354971244075328+y*(236968812099024590167053471+20*y*(-38254237482625477235750804+25*y*(-28832614495174597470411+8*y*(5612513556492807328773+96858236340211903900*y))))))+1250000*(x)^(7)*(-90008106283614443392212407123+y*(-8432585990927600005430502072037+2*y*(73725322169426088714622191129+20*y*(33267867471633597699708373449+50*y*(-11859166675181326539831076+y*(-37751660152923343565540931+10*y*(49339595481767764811939+78022255196391533949800*y)))))))+43750*(x)^(6)*(84576571841663318667353827852358+5*y*(-119270021041172613069243491657009+2*y*(-5267006730739703066455523697141+20*y*(365096653088222551419951669109+50*y*(720304994113427524676916996+y*(-262145884199479055265345723+10*y*(-3766966184007510154199879+40*y*(4858565047001895189461+1654278215577084353790*y))))))))-9375*(x)^(5)*(-3944665406579716882953309808152639+y*(-59631788330605795969226108804524601+2*y*(1175294921656112003132297904715299+20*y*(247260307856921203334807082003489+50*y*(-133087026003182851788371621127+y*(-280966662222325225695067978251+10*y*(462831431276262238433359389+1000*y*(543109951341832710961281+5*y*(-445873082654359563939+60986370928455606104*y)))))))))-1250*(x)^(4)*(86667944396611326125152191968589291+5*y*(-92711306107668717189746378411325493+2*y*(-6137440113826370431243828948588794+5*y*(1221053113991921726225936571149394+25*y*(7199094000253875689933523616543+4*y*(-487098523544301830556854507817+5*y*(-19400717898846330756546159071+20*y*(48841879472535873837618516+5*y*(3460946741435033281024677+200*y*(102652295928016054801+11312713918519394172*y))))))))))+20*(x)^(3)*(-43716237256277211654117871080508759761+5*y*(-100504857453883898313383545165878863533+10*y*(494455112403789248890687184221654578+5*y*(339712344190683828714051811107583207+25*y*(-435352475067945137636850632278106+y*(-775588456340703311348419990354869+10*y*(1551638818925981834611256221577+40*y*(37254802836739713592089212011+60*y*(-5055376602138171854610588+25*y*(5433891171018067058939+4*y*(-30756759099767012807+2078936708753422440*y)))))))))))-3*(y)^(2)*(-517187094791941512207875623483647488772+5*y*(77217754052102660111443949529475785072+25*y*(633448828506561189543976944087298763+5*y*(-27988221239870592477014513717078484+y*(-7102809297961140075487637866868499+20*y*(17876842976762970788426265473984+5*y*(1368773374720284161388360799551+200*y*(64129244261551769175819087+32*y*(-194672732009669889214129+90*y*(-68395720077821388488+13125*y*(-68241190290897+784*y*(2293809571+56115000*y))))))))))))+3*(x)^(2)*(357008814843850566704563621317729509638+5*y*(-316328078054897574117811933975431492486+5*y*(-10541312754448979091768159561196028679+5*y*(1754523182153765065935875096330768944+25*y*(12679802486482816224180486557374431+10*y*(-300432851590565251710084441841767+2*y*(-34770495496633093130101989295657+40*y*(43348784030460636922030771701+5*y*(3181431815859714378054765579+200*y*(131321929928275706385031+18*y*(273910185345199868311+240*y*(31713612288757402+1497839540845375*y))))))))))))+x*(5975401493509874489718402548998234307088-25*y*(-2425969083067941644958089230894270786564+y*(134823314409260395756442128075036963626+y*(415920671685559246029616085390606549878+25*y*(-600479315071868234907855631100928349+y*(-961533574482894251649317866257559143+10*y*(2260764901252985584405566308733511+200*y*(9504643318524626632512636282979+y*(-116341862738567034959898413679+800*y*(-9882119516256426019546181+6*y*(-2056855346844450996833+180*y*(-6902072649829414053+1750*y*(-12582105727243+1373281449960*y)))))))))))))))/7.001316e40",2);
        f=gsFunctionExpr<real_t>("(200483015939927424805749401753922095656048-20792963040042998776819931131812603569760*y+5*(85796184373560849187500000000*(x)^(12)-7500000000*(x)^(11)*(188763309224212970923+5082269320562618300*y)-1875000000*(x)^(10)*(-69162007436240087786159+8*y*(2773202004406440477929+3894193987465008424545*y))+25000000*(x)^(9)*(-149909623908054866557383907+400*y*(-63214168016399689867944+5*y*(14015207218248956806026+598454433048046209625*y)))+7500000*(x)^(8)*(19473177179241234640028364387+50*y*(31862975136666387602469487+10*y*(-16355455723645337460846387+10*y*(-112982613228639206194+109202032280810305350935*y))))+750000*(x)^(7)*(1572537911999558169581729836189+20*y*(12594309484008206716874150172+5*y*(-4165815250858270258696456113+100*y*(-4733092593717614720471873+40*y*(40903113642087276773771+5340069121376618277447*y)))))+12500*(x)^(6)*(-2271988464656030168577251066789687+20*y*(2944743138851483077922538850602+5*y*(9910947021964690206354762839607+700*y*(-292940854633575377776655719+40*y*(-28505650653407207205195069+279668611635311480169144*y+503637821082003617841650*(y)^(2))))))+7500*(x)^(5)*(-11404202153620510282266286367953621+40*y*(-41430614895086613420012978375721+10*y*(7300076807925748182243512523507+5*y*(207798978778270540024515420577+5*y*(-15704605651331332650308349669+20*y*(-131192453670640178453814723+100*y*(90447654906181480116769+18546532481001706820208*y)))))))+375*(x)^(4)*(3679660093973470770621986393970644679+100*y*(-2109005546190532388780547854257611+5*y*(-2849431322847162531737397627418393+20*y*(4557430842446960617713592202719+5*y*(2392150695428852486560506087695+4*y*(-10301361341987745035225899043+50*y*(-276806349418032697590524561+100*y*(18342483418404099256747+12322864260311585068187*y))))))))+100*(x)^(3)*(15056257746185269082445666631937465423+125*y*(20980001691361741808527279021965765+2*y*(-15841241175542331931316536052702589+10*y*(-285162413610154714654165106097818+15*y*(6409883038339352313893789206357+4*y*(333508619420475236878379216642+25*y*(-1023407103086551980758209087+20*y*(-12877828294431805893821999+5*y*(22595898331228390795931+8841866686884050087760*y)))))))))+5*(y)^(2)*(-3551531163988261907495742883855209660492+5*y*(34381021327838633385934483932486896924+5*y*(11829420487260040622572021242890193237+20*y*(-16049295832071177139434109298296377+5*y*(-2419016927062141523648868713984609+300*y*(122010981956094749795023301659+y*(44005037843774091904242393211+80*y*(-3472649528000044723949072+15*y*(7892045651005192672523+40*y*(-2397054208507204837+142550086209543870*y))))))))))-30*(x)^(2)*(565149859396440035224265578183981607208+5*y*(-8437287002057048996772266807335117994+5*y*(-8716574128893158922537895479104568879+50*y*(6918879623452658688632669271673658+5*y*(2958282868548894096977145766528735+4*y*(-16560015850331253540195841546284+5*y*(-3557510113118293272483686100191+60*y*(652163400886916781761946087+125*y*(2624182264153707638925741+8*y*(-999457787323601761019+135161035690724068988*y))))))))))-12*x*(530925557520583884799110854193758960732+5*y*(24376442401877963241185423320934974674+5*y*(-5999125722417512746535516369169733891+25*y*(-54851286351950879916482325313472751+5*y*(3098631074456770371333310355311893+8*y*(98897719114710296309142369775414+5*y*(-1337911001230730072008322933889+10*y*(-40430425105286506983844478437+250*y*(1301566201057324263051717+2*y*(363125623376994231239461+4*y*(377139568951883178101+44110620057322458600*y)))))))))))))/1.4002632e39",2);
        iterLimit=5;
        break;
    case 10:
        ss << "5-patch domain, degree 3, v2\n";
        pathStart = "Mario/Florian_FivePatches_deg3_v2/multipatch_five";
        pathEnd = "_degree3.xml";
        mapPathEnd = "_map_degree3.xml";
        startIter=1;
        g=gsFunctionExpr<real_t>("((51.2 - (32*x)/3. - y)^2*(-17.9 - (23*x)/5. - y)^2*(-6.855555555555555 - (13*x)/9. - y)^2*(3.7956521739130435 - (10*x)/23. - y)^2*(-4.41578947368421 - (2*x)/19. - y)^2*(3.0448275862068965 + x/29. - y)^(2)*(4.45 + x/2. - y)^(2)*(-4.923529411764706 + (9*x)/17. - y)^2*(6.1 + x - y)^(2)*(-21.6 + (9*x)/2. - y)^2)/1.e16",2); // 3
        dg=gsFunctionExpr<real_t>("((883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(5184*(-24+5*x)*(-60759448049004644133+10*x*(9224871288649613571+100*x*(37604179460108437+10*x*(-85805440864572+5*x*(-23585288349969+2*x*(-344361452181+40*x*(2225326239+88527705*x+747500*(x)^(2))))))))-135*(17023674641706872782653+20*x*(-1695509837688518024354+5*x*(-41960889498677262899+100*x*(284438575105026496+5*x*(2833672964522191+12*x*(-96263698034234+5*x*(-270006365587+80*x*(969070628+45871135*x))))))))*y-250*(3209323891457105431287+4*x*(-1063905721177541644621+20*x*(-10863701126609174361+5*x*(311839709514378744+25*x*(4409023124019543+4*x*(-39053380365891+50*x*(-47676560757+30484236416*x)))))))*(y)^(2)+500*(414137231579450182133+100*x*(-8503573356019086696+5*x*(-169890059042411955+4*x*(19655998498441548+25*x*(66279852581471+60*x*(-104658896724+22349641013*x))))))*(y)^(3)+25000*(669962652693164829+4*x*(-252981723406319757+20*x*(-3480355330614159+5*x*(-149651754464292+25*x*(1579945392551+246871418796*x)))))*(y)^(4)-50000*(36497337774982877+20*x*(-13208323166449558+15*x*(-26473933807981+40*x*(345529412146+198697687975*x))))*(y)^(5)-500000*(-1618566060310711+60*x*(-11342941558297+40*x*(-39209104117+12314778830*x)))*(y)^(6)+15000000*(-8007076526467+20*x*(-437203258412+139024822065*x))*(y)^(7)+50000000*(-475397417077+52105762340*x)*(y)^(8)-621668112000000000*(y)^(9)))/1.69182035233245e49",
                                   "-((883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(9*(-24+5*x)*(-16077059186642225498133+4*x*(-1369336665190955685819+5*x*(48913670472575889613+20*x*(946609999700553672+5*x*(-5001444039303987+100*x*(-28130131110147+x*(156037145859+40*x*(2620771809+91742270*x))))))))+20*(145441291720013932895181+5*x*(3209323891457105431287+2*x*(-1063905721177541644621+40*x*(-3621233708869724787+5*x*(77959927378594686+5*x*(4409023124019543+10*x*(-13017793455297+100*x*(-6810937251+3810529552*x))))))))*y-300*(4238909567658291807147+x*(414137231579450182133+100*x*(-4251786678009543348+5*x*(-56630019680803985+4*x*(4913999624610387+5*x*(66279852581471+300*x*(-17443149454+3192805859*x)))))))*(y)^(2)-4000*(75855763148894346513+5*x*(669962652693164829+2*x*(-252981723406319757+40*x*(-1160118443538053+5*x*(-37412938616073+5*x*(1579945392551+205726182330*x))))))*(y)^(3)+50000*(1644683832982118973+x*(36497337774982877+20*x*(-6604161583224779+5*x*(-26473933807981+60*x*(172764706073+79479075190*x)))))*(y)^(4)+600000*(18992823414703953+x*(-1618566060310711+10*x*(-34028824674891+40*x*(-78418208234+18472168245*x))))*(y)^(5)-21000000*(64762860701939+x*(-8007076526467+20*x*(-218601629206+46341607355*x)))*(y)^(6)-80000000*(1985626097571+x*(-475397417077+26052881170*x))*(y)^(7)+2700000000*(-500083579+414445408*x)*(y)^(8)+116338140000000000*(y)^(9)))/3.3836407046649e48",2);
        ddg=gsFunctionExpr<real_t>("(80*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(9*(-543914379+1936272910*y)+100*(162*x*(-1003582+5*x*(-993+3568*x+50*(x)^(2)))-18*x*(-3185437+643605*x+61400*(x)^(2))*y+(-59423569+20*x*(601672+190725*x))*(y)^(2)-50*(25691+84020*x)*(y)^(3)+1294400*(y)^(4)))*(2619201987768-935774142455*y+50*(16*x*(-2700305626+x*(-469551171+80*x*(1108018+37375*x)))+x*(-11430373379+20*x*(132519969+24500680*x))*y+4*(-1888260208+5*x*(102280583+48691290*x))*(y)^(2)+10*(2825021+47307530*x)*(y)^(3)+50886600*(y)^(4)))+25*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*((883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(837-90*x+170*y)^(2)*(81*(-5687+60*x*(13+10*x))+180*(167-330*x)*y+15700*(y)^(2))+8*(883+10*x-290*y)*(-837+90*x-170*y)*(117+180*x-110*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(9*(-11809+60*x*(293+5*x))-200*(1485+296*x)*y+104900*(y)^(2))+4*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(27*(-38127151491+30252883840*y)+40*(243*x*(-8844351+5*x*(749023+50*x*(586+5*x)))-360*x*(7615641+10*x*(217413+3700*x))*y+10*(-75702173+60*x*(4329504+266245*x))*(y)^(2)-300*(7743339+1596890*x)*(y)^(3)+421028250*(y)^(4))))+4*(883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(837-90*x+170*y)^(2)*(20*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(-3344+7360*x+1145*y)*(-5422463+240*x*(8226+325*x)+2406380*y+649600*x*y+780500*(y)^(2))+(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*(256*(-1932479+1380*x*(-209+230*x))+320*(-81373+79005*x)*y+1863025*(y)^(2))+25*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(-59816709411407-65166891181640*y+40*(6*x*(-163190328849+20*x*(865547273+6500*x*(16452+325*x)))+20*x*(3646882989+40*x*(143693421+5278000*x))*y+5*(-325804009+360*x*(135707027+12173420*x))*(y)^(2)+300*(536897933+131866100*x)*(y)^(3)+21616698250*(y)^(4)))))/1.69182035233245e49",
                                    "(-40*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(8106539481+9*x*(-193627291+10*x*(-3185437+10*x*(42907+3070*x)))+171788588*y-20*x*(-59423569+20*x*(300836+63575*x))*y+60*(-19223711+25*x*(25691+42010*x))*(y)^(2)-800*(-127477+64720*x)*(y)^(3)+9860000*(y)^(4))*(1804179436833-1545589725720*y+4*(x*(-187154828491+5*x*(-11430373379+40*x*(44173323+6125170*x)))+40*x*(-3776520416+5*x*(102280583+32460860*x))*y+300*(-347459898+x*(2825021+23653765*x))*(y)^(2)+6000*(-507835+339244*x)*(y)^(3)+147487500*(y)^(4)))+4*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*((883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(837-90*x+170*y)^(2)*(-2327+3925*(x)^(2)+300*y*(31+2*y)-20*x*(794+165*y))+400*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(-31+11*x-4*y)*(-115123+1480*(x)^(2)+x*(14850-10490*y)+3*y*(-8451+4930*y))+2*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(983008657193+3696921147420*y+20*(26624500*(x)^(4)+600*(x)^(3)*(1443168-798445*y)+15*(y)^(2)*(-3897676511+98600*y*(-8451+2465*y))+10*(x)^(2)*(-75702173+90*y*(-7743339+2806855*y))+x*(-70039005912+10*y*(2650634259+10*(167638491-51715700*y)*y)))))+(883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(837-90*x+170*y)^(2)*(8*(179+46*x+10*y)*(-768+160*x+15*y)*(-999+458*x+60*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(-4919953+324800*(x)^(2)+20*y*(294019+58995*y)+20*x*(120319+78050*y))+(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*(298084*(x)^(2)+60*x*(-16589+1374*y)+9*(-72407+60*y*(-333+10*y)))+4*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(-241540077242153-193437613834440*y+40*(3652026000*(x)^(4)+600*(x)^(3)*(135707027+65933050*y)+15*(y)^(2)*(47746822063+65550*y*(588038+58995*y))+5*(x)^(2)*(-325804009+60*y*(1610693799+432333965*y))+x*(-3152734415003+30*y*(5512681989+10*y*(3241247549+511617750*y)))))))/6.7672814093298e47",
                                    "(-11664*(24-5*x)^(2)*(3795508851241428081609835296757521036403467+8*x*(-617924889466712236003969418712490006166196+5*x*(-117799961390115888244607713916466201335763+10*x*(-1153754925035254244930228232162930630276+5*x*(190783077681420820239961284719819672741+400*x*(94136542817661227100453777189555504+5*x*(-2335776705045331973430661489837887+5*x*(-139662379522414132011851610140952+x*(2989945885086361921720201255521+40*x*(55027619254438667298341299236+5*x*(193110516368558089453986589+10*x*(-7645478313504331939059012+5*x*(-65763204980290218177741+160*x*(19183304824653527736+5*x*(350631176011193151+40*x*(197697136981902+1302969589675*x))))))))))))))))-3645*(-24+5*x)*(-133217222313490691110038761522914324054767573+4*x*(30707381976125011223697233902455332398970319+20*x*(1191165278871544872251481597207897295736307+25*x*(6059617791835945705996637645576510043311+5*x*(-587254050168334401901216944803729861137+80*x*(-2857926878484637578930592938424051364+5*x*(-1924196692092154206887727370076932+x*(24731098175816427171948650852354142+5*x*(298808202379124496086373676196557+20*x*(-4061264643115522705627614150193+20*x*(-9322530128165611705629893263+5*x*(254749554884573818720518079+5*x*(-2504987654656513967213119+40*x*(-3735831775608908438463+80*x*(11411323124051659378+15*x*(23233173342897719+23271419454370*x))))))))))))))))*y+675*(25695529994036350691195886821931982550755068915+2*x*(-21222561349774070294281120924390666690935287463+5*x*(-1470192296359704686899216335439906402872612009+20*x*(32828438928021100531095137055388422092801442+25*x*(375771878930388547769390107169649955359629+8*x*(-9297559430780970006311141091744818659836+5*x*(-518824950630865729493051643271960121031+40*x*(1239337523806790172465481593355752474+5*x*(70298941463740785448635833708320698+25*x*(-165266538000479350228250648119671+x*(-34434605624479501963715560448477+20*x*(127192476002036628020827168206+5*x*(1247946149098604414584565817+200*x*(-1222897611049032857462391+10*x*(1161178534000944312057+320*x*(1266531504681501663+41525418232972513*x))))))))))))))))*(y)^(2)+1000*(454908891261705405293290728846770879955211962+x*(1192377384387916217783146947465730421600899239+10*x*(4323767282325281414116112562531325444590531+20*x*(-3432367826491467794474142636672332254331339+200*x*(-829274239270937985012716379428979999543+x*(1866368217671403190852503064272450591753+5*x*(14660694531649725548851755369383919711+20*x*(-814066796197809791333767300226627538+25*x*(-845706731082957996585568457594364+5*x*(136341428904494516644467134448151+2*x*(-523686404450488294038346109853+100*x*(-6401149892629696416060739017+20*x*(17276974088396173674640689+x*(893137738722155892236061+100*x*(-308595463037841925575+35323694640674483104*x)))))))))))))))*(y)^(3)-12500*(203547927545787429973647834403689761163930297+4*x*(-88700405875258812633570461844349224850692594+5*x*(-5088884454361057926285129444449369229933843+20*x*(102037029733931730830901970205868315022536+25*x*(1139056577223180365690095879519796153237+20*x*(-5997254477044199421896067464297433513+4*x*(-675585137903823909090519779730570562+25*x*(-395141985092314217490706644985944+x*(588273024309087438974071601399403+100*x*(304204911714318416453141128308+x*(-62276343587707471314832158463+20*x*(-20909933786942601669098856+5*x*(2866213319543382630659383+60*x*(512820645506886052273+550407119691908322750*x))))))))))))))*(y)^(4)-150000*(-820912876560040055454432453797226841034706+5*x*(592450044345189421336013296999840871075323+80*x*(1311506726167945985871175023744316137447+x*(-1803079355393752025371711836650070599656+25*x*(-7267115612631854965848892570016054753+x*(4807123355570829949579856201052443019+80*x*(6152635135507867630226220089665361+5*x*(-309672464598060704832373938469532+5*x*(-7725830568381436150542669393672+25*x*(30980980740115059043000786097+16*x*(172350682542795431883899746+25*x*(-796881288663087402354876+5*x*(7391466654891332268515+236379089094056247129*x)))))))))))))*(y)^(5)+1750000*(80773870313750739897812091863957297184759+2*x*(-92643044110063114738520207042823764986917+25*x*(-719117177967701738403704851801130970015+4*x*(77634810668446222272798278038044225162+25*x*(766259340515986553579744402873911759+16*x*(-1214932437321533156667985363337886+125*x*(-15564052439637071236626757180113+2*x*(-844742087491991094270169700424+x*(80097624826001896502857819221+50*x*(405967284015213374071442325+x*(-10201816907430873985157249+60*x*(12140065785830258340462+8268842088053710100275*x))))))))))))*(y)^(6)-20000000*(539426160540983772054245755685428914630+x*(-1790941405357047567376315431779649452431+10*x*(-18091672375436750457831888150781079637+20*x*(1591626137629073876880994869869520241+400*x*(263805077291529704155031159940002+x*(-168796794272635732890337453811646+5*x*(-4533024644089582279726522153397+50*x*(6604069465120570812393238904+x*(2427734883696585912749793852+25*x*(5319031552190841163885487-970751473560998415598774*x+54273906819203417628840*(x)^(2)))))))))))*(y)^(7)-225000000*(8279588480765903113524666213672961341+4*x*(-14802325646602365620237113770635840653+40*x*(-24629286857699112254909947458500148+5*x*(4750519337020804929328849712882428+25*x*(34961665971044349794100633258143+4*x*(-97995562129368317132554465683+10*x*(-34798545644587604875845990133+15*x*(-276055521378564648930042248+5*x*(-3817133550493210740612501+100*x*(13829622398451154550069+2124211143820855593988*x))))))))))*(y)^(8)+2500000000*(138552077559380492720100083552259186+x*(-518129521994851713755497848220513949+100*x*(211060585084324952228717340164817+4*x*(181997514524528155615390486535897+25*x*(-191528781724603748936085265405+x*(-218178188813585401598670114969+80*x*(-316192049271703378831743693+10*x*(3242759559694185931373896+5*x*(43148931129595775808171+31046168331245185047725*x)))))))))*(y)^(9)+27500000000*(-4170899306537626970232096726120063+2*x*(-5594389308271109900804197339653193+25*x*(22937196753854554456622278215477+4*x*(2782364666669841709901242202438+25*x*(12655247052698824086276011299+8*x*(-362306615452479131928662592+5*x*(-21680317993242996375333143+180*x*(799650858117221408548+2385269629366586971095*x))))))))*(y)^(10)-300000000000*(32197381414060723715679384929442+x*(-74083180728190222571633627817181+10*x*(1787978184548109081564181844403+20*x*(35839990317274869045377525197+200*x*(-38454844495213594178907107+x*(-3614482475652951149560803+5*x*(-92854244493612083309023+50260732806860363055480*x)))))))*(y)^(11)-3250000000000*(-1386001602410594435486708400885+4*x*(-255050973614285344755027951184+5*x*(19773915409078463607981525441+20*x*(16772096711582432292766752+25*x*(-218069407031470220897891+12*x*(-18393983188864256510007+2078881834840572039740*x))))))*(y)^(12)+35000000000000*(9527124293161804308265821318+x*(-5150850242757708061721573957+40*x*(48233241498574572082579983+10*x*(-54691291943440779941534+25*x*(-15874909021313338081109+2545236255807717734007*x)))))*(y)^(13)+1125000000000000*(-37864715645896890114243841+2*x*(-4200264675491427861837169+5*x*(1072372985577804280726881+20*x*(-11343858139877762710678+726455486755366925275*x))))*(y)^(14)-4000000000000000*(1240395513980207205085266+x*(-341925624484588249279597+250*x*(-92745212797983746355+48508078875627116716*x)))*(y)^(15)-127500000000000000*(967519220929799627413-521836616853940347052*x+64740269323974620280*(x)^(2))*(y)^(16)+8100000000000000000*(203631613331080255+52205345473267412*x)*(y)^(17)+68707526255022096000000000000000000*(y)^(18))/1.69182035233245e48",2);
        f=gsFunctionExpr<real_t>("(243*(50497052605010560134063815435000252508733486047-199014893481648118666510199942247945744424463900*y)+20*(3*x*(212407118668070050400837154153816241340492165898+5*x*(-18417110037621892049763610258305244835011844133+20*x*(-680471534018737869090296316337848469115026083+80*x*(1175459346529485922139611944122104414210906+5*x*(172112960799268614739693592859142901011962+5*x*(-1716048929701797561985523204366578140681+25*x*(-61540159884695710531482095952281403942+x*(844705254259547966483851405161850047+100*x*(13311970415181511628653825075630350+x*(21645699412432041655030946513971+20*x*(-6840738135683690390772677187123+10*x*(-5518070219015182905821289939+100*x*(27808905458626020861295272+x*(366953170538786299641243+40*x*(-348306547122340607949+11916893274269801150*x)))))))))))))))-30*x*(33334846578718803712836745187308061439617260371+20*x*(-1235062268868930927860816357677686404097882117+100*x*(-4162363803420690873199447358185652490590699+20*x*(59354126418490270907082165869470140161808+x*(17547035779165263316468961476545420914483+5*x*(-494316158274430219657562901683749179639+10*x*(-13028169299228303533675304608094508296+25*x*(40844288345705546058646593265867851+10*x*(877548408929790707935605862947371+4*x*(-11672414521425367319178577561041+20*x*(-78110130302860872116052661922+125*x*(30690393923047724427250491+2*x*(795065881948686441543877+30*x*(-967170707027605271513+16218672343410725572*x))))))))))))))*y+180*(546784633884459797216264789039179116365722161+20*x*(-1261343022606605000579609325372221472989979+5*x*(-7910031470487984440969074309689579727435089+100*x*(249456105019200303851154822522557683342+5*x*(2658778978542364937119463883012840594164+x*(11837726480363995088877802606471616078+5*x*(-30236551928764998335683280960199343251+400*x*(-577270751437851340335378275133493+50*x*(34382298422256974073635525703362+x*(-157782044335771471639308494798+x*(-281100749404461107701887202019+200*x*(73673101806336391876342283+x*(2883637905851808887742981+25*x*(-18997479527234621613806+528402799876198962725*x))))))))))))))*(y)^(2)+100*(9017290723668017476849352046018796366923229347+250*x*(9492778262890475574690917992561811624648949+2*x*(-3672016492380678559820549609323939052626833+10*x*(-98747911651863524864358293730517276615679+10*x*(1883360543722071961245940708117682150673+2*x*(348889663344272393905935290004078564951+100*x*(-83055242212311260548420447421567781+10*x*(-9829995116001099856347546784538811+5*x*(-42441595642999529627470808931927+10*x*(1989955793891110454441595706607+30*x*(711566141441201205351729+350*x*(-731935742391105870558027+51939955730717659255126*x+739754669154418789140*(x)^(2)))))))))))))*(y)^(3)-3000*(21673824097841548582257881194967159573305707+100*x*(67307062718311535701271568890582266771212+5*x*(-27087802405204651505569485170902870948459+10*x*(-259091505734641385884563506842733972788+5*x*(57715317851430521102905951386608284001+40*x*(115217684517907150380276089395585652+5*x*(-10349458004456435064473344099802011+20*x*(-47464361636532608838961727858717+25*x*(291261872814751752245258403545+8*x*(2661511493441287143316932561+x*(-189918650537643461574473239+25*x*(-61934189881380958360838+30952138707475426534555*x))))))))))))*(y)^(4)+150000*(-755560363725838359148560498360055616719431+10*x*(-9398389697373280484814341607493095354511+8*x*(1415094952708344193780829970553508062152+5*x*(48606220060322909274751192041254372533+240*x*(-27043157640897968667989883659190327+5*x*(-2782368102220777091017371419980194+5*x*(-43768869279446291152036106405198+5*x*(2205506363056528264195298327218+5*x*(61168109026728820251894542829+70*x*(-16747992317543951268646941+40*x*(-96319517268418945191139+271324349632049063355*x)))))))))))*(y)^(5)+200000*(29104182198910131515974604531722603398751+5*x*(1597458383515281765931390171701584744226+5*x*(-665921640312003440225668625532789833643+400*x*(-73846786692495739528068979647246829+5*x*(24170024772592512246368424419143773+20*x*(96536989159256977948711411909977+10*x*(-2280584745509951477934817897207+2*x*(-223890845213618597042804274732+25*x*(66422983416824852434601577+5*x*(13567676776997201615030462+627502089657674564902863*x))))))))))*(y)^(6)-3000000*(-2183141472399162879850734636193543514253+10*x*(5203847932316006680924487269794620689+60*x*(425753650834482140037016531334563507+50*x*(398272255709993952277980984515093+10*x*(-11159307986682461224941277919733+2*x*(-2740157694693315875799124637527+10*x*(-24112942836407663765239735167+10*x*(208477385435870197417335587+25*x*(2707721917645645200751029+117992421195564066894710*x)))))))))*(y)^(7)+60000000*(-1858439803665643870582093506184009744+25*x*(-29764471284054843178141569434243102+5*x*(15236519769154931499868486299686379+20*x*(-52362908417205495508877028176011+10*x*(-3641558368086216508169976721159+12*x*(-4920755046711964831324968936+125*x*(26616532254093248644321669+2*x*(2536617248288587614014584+69553362429410636411685*x))))))))*(y)^(8)-100000000*(1814384361101580973443885929252954683+50*x*(-6787952704150918425822996742550169+20*x*(-155913858053357011016072814192261+20*x*(857024383932366856343012817511+180*x*(411909157530444376154960557+x*(146029709080024790265934293+125*x*(-71501768135188090199793+497222670838439389888*x)))))))*(y)^(9)-6000000000*(612080638310556845374925633967631+10*x*(-27057853148622611920984009081132+5*x*(5488885379474337745883239715903+40*x*(-29852035998534663201955877811+35*x*(-69102380756204146241290989+5*x*(2444066777750734913761658+126567460808479706709525*x))))))*(y)^(10)+30000000000*(66035767323965342858727737494787+10*x*(-2515527471708063637361944555387+30*x*(-7813199029223015950393815429+50*x*(90692135802509159620384297+30*x*(-191718947417303873969225+9917257402117231729682*x)))))*(y)^(11)+100000000000*(1068028978938837866922026856161+40*x*(-11388528454603134464531889237+35*x*(72287855979385413512289024+25*x*(-625337140047818638548434+42753880513474669329135*x))))*(y)^(12)+3000000000000*(-286768377678985797203887497+50*x*(9139797036319291657602643+120*x*(-19308984937116675326462+438275889262669832745*x)))*(y)^(13)+60000000000000*(29731229844053438771703+25*x*(9231283041521873999582+1049228045018022708195*x))*(y)^(14)+300000000000000*(11704072636745431716719+1160590434881456227170*x)*(y)^(15)+79296548109497025481500000000000000*(y)^(16)))/1.69182035233245e47",2); // 3
        iterLimit=5;
        break;
    case 11:
        ss << "3-patch domain, degree 4, v2\n";
        pathStart = "Mario/Florian_ThreePatches_deg4_v2/multipatch_three";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        startIter=1;
        g=gsFunctionExpr<real_t>("((-47.7 - 14*x - y)^(2)*(-3.8774193548387097 - (16*x)/31. - y)^2*(3.875 - x/2. - y)^(2)*(-3.7240259740259742 + (39*x)/77. - y)^2*(3.8243055555555556 + (37*x)/72. - y)^2*(-286.75 + 77*x - y)^(2))/2.5e11",2); // 2
        dg=gsFunctionExpr<real_t>("((5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-441389437313011+5272968224115*y+8*(298666368000*(x)^(5)+28000*(x)^(4)*(-4133821+672325*y)-480*(x)^(3)*(52101274671+5*y*(-48329457+640401050*y))+12*(x)^(2)*(565843921171+5*y*(14511590237+20*(46516301-64021780*y)*y))+(y)^(2)*(8076277817723+20*y*(-58621313017+200*y*(-59499360+10828163*y)))+2*x*(229800161083395+y*(2082084269689+10*y*(-1675159861883+20*y*(207397103+3706248990*y)))))))/1.5123064061952e31",
                                   "((5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(31*(-65076108692243+567148797931020*y)+4*(x*(5272968224115+8*x*(2082084269689+20*x*(14511590237+10*x*(144988371+18825100*x))))-16*x*(-8076277817723+10*x*(1675159861883+680*x*(-2736253+56505975*x)))*y+12*(15296479392949-40*x*(58621313017+10*x*(-207397103+128043560*x)))*(y)^(2)+480*(-624025745233+400*x*(-39666240+123541633*x))*(y)^(3)+8000*(-1149261923+216563260*x)*(y)^(4)-32997888000*(y)^(5))))/6.0492256247808e31",2);
        ddg=gsFunctionExpr<real_t>("(2*(5507+740*x-1440*y)^(2)*(1147-308*x+4*y)^(2)*(1147-156*x+308*y)^(2)*(139686329327-200681825360*y+25*(16*x*(-362657382+x*(-101620989+28000*x*(355+84*x)))+8*x*(-373082571+20*x*(4913781+1257200*x))*y+(-653061619+240*x*(5025693+1529365*x))*(y)^(2)+110*(4996423+2317980*x)*(y)^(3)+48905650*(y)^(4)))+8*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-489329699+26666640*(x)^(2)+8*y*(15792077+4325450*y)-8*x*(8077099+8740640*y))*(-330574+75085*y+25*(672*(x)^(2)+931*(y)^(2)+4*x*(355+449*y)))+(-31+4*x+8*y)^(2)*(477+140*x+10*y)^(2)*(601+80*x+155*y)^(2)*(-3441*(-35572636019847+73096675692352*y)+16*(74073925926000*(x)^(4)-44444400*(x)^(3)*(8077099+8740640*y)+4*(y)^(2)*(-172228484745029+20*y*(13798425242429+959704939205*y))+6*(x)^(2)*(-478459593600089+80*y*(3519648832631+1435589920070*y))+3*x*(3988696521770833-8*y*(-206187448112451+10*y*(10414257900957+1896520058240*y))))))/3.780766015488e30",
                                    "(-32*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y)*(-31+4*x+8*y)*(477+140*x+10*y)*(601+80*x+155*y)*(-184879+179600*(x)^(2)+60*y*(19717+620*y)+40*x*(15017+9310*y))*(-9535011+17481280*(x)^(2)+4*y*(31783721+166320*y)-4*x*(15792077+8650900*y))+1573679923200*(-31+4*x+8*y)^(2)*(477+140*x+10*y)^(2)*(601+80*x+155*y)^(2)*(((11119+113140*x-221760*y)*(5507+740*x-1440*y)*(-1147+156*x-308*y)*(-1147+308*x-4*y))/(2.45887488e10)+2*(-3.7240259740259742+(39*x)/77.-y)^2*(3.8243055555555556+(37*x)/72.-y)^2+2*(286.75-77*x+y)^(2)*((-350064735599+200*x*(15438523+96003266*x))/(1.22943744e10)-((11119+113140*x)*y)/18480.+6*(y)^(2)))+(5507+740*x-1440*y)^(2)*(1147-308*x+4*y)^(2)*(1147-156*x+308*y)^(2)*(-10479303116099-1317338496540*y+40*(4*x*(-41626408888+5*x*(-653061619+40*x*(10051386+1529365*x)))+330*x*(130950371+20*x*(4996423+772660*x))*y+15*(3452990809+40*x*(58793797+9781130*x))*(y)^(2)+62000*(59151+18620*x)*(y)^(3)+57660000*(y)^(4))))/6.0492256247808e31",
                                    "(263879111196180480000000*(x)^(10)-200704000000*(x)^(9)*(-8649747410767+243871656618475*y)-1102267*(-249908537057869493753722+7303900366010386651975527*y)-38707200000*(x)^(8)*(290585708203906+3*y*(-197029885967433+80994558247885*y))+8192000*(x)^(7)*(-9481688358976200867+5*y*(78477855207928463297+50*y*(-12621679664055075+221243488766886308*y)))+2867200*(x)^(6)*(-280977579189855473287+5*y*(-22493549631768585111+5*y*(-2747278652088440511+20*y*(-91855653718880188+46208541740778675*y))))-92160*(x)^(5)*(75483626813060289160351+25*y*(-20046466182105866760727+100*y*(11594160804861041637+20*y*(-887558443306536565+4*y*(1720896663954649+71090383325622549*y)))))-7680*(x)^(4)*(-4964239415059682722723998+5*y*(2296897469493749744928169+25*y*(4891483627964772371247+40*y*(5266351819578766236+5*y*(-13129448820140085149+80*y*(-26146153193635401+10391360459106683*y))))))+1280*(x)^(3)*(319898962438204378868630125+y*(-3820925583992735261392046391+10*y*(-4783100735284830070720533+20*y*(2105413765947766789130671+40*y*(-281075782689312175923+5*y*(-527029110724154050557+8*y*(1102027571944353971+2738805988130334580*y)))))))+(y)^(2)*(617841001190024474584165568361-80*y*(-21074739527441629699914866180+y*(2607549751484824916038727667+16*y*(91673462096616111761339829+125*y*(-107496658047501834366731+80*y*(-225835710267402201780+y*(31492159069935608973+80*y*(12239879471796649+40941370928304*y))))))))+288*(x)^(2)*(-920093625250723826423869547+10*y*(892691419041710198843985595+y*(-47417212378914652851450213+20*y*(-6547575133692972562670908+5*y*(131944367478582199929979+32*y*(2394158422480936839957+5*y*(-79873367397145988827+200*y*(-56248512753087052+12046162278210363*y))))))))+4*x*(-1197623698331623912863154084953+y*(12321028396612468035229975389749+40*y*(16186340135808633803985384165+4*y*(-13097715815752721726226123093+20*y*(-24778019926500117047201121+y*(43068787960711552962518799+160*y*(8538946273367625114111+10*y*(-695886202065453168197+5*y*(-4066084107601811649+10710407262566600*y))))))))))/7.561532030976e30",2);
        f=gsFunctionExpr<real_t>("(775*(10070636445024446987294496713 + 1580805041435296880732831616*y) + 32*(839212170045134580480000*(x)^(8) + 3072000*(x)^(7)*(-186153527560215013 + 23630822384505615*y) - 12800*(x)^(6)*(7397845103057974754479 + 600*y*(-1978168936115700 + 521427067157777179*y)) - 23040*(x)^(5)*(-2591139538241579069079 + 5*y*(-19502249110340706949 + 40*y*(-532005332078551427 + 231162700790668990*y))) + 240*(x)^(4)*(20629697508869638693714899 + 80*y*(-2832250234210447845948 + 5*y*(361036991177109628031 + 80*y*(851520842693099156 + 1397012256516649535*y)))) + 480*(x)^(3)*(-5011292996454793182398419 + y*(172125474399395900332421 + 80*y*(6788570295304368834777 + 10*y*(-121146247699783170709 + 40*y*(-983617023935690419 + 258889960348724901*y))))) + (y)^(2)*(-113362049203585367125403749833 + 80*y*(-133919062100907890758222638 + y*(119712510000704219793039297 + 16*y*(506422009899324159806439 + 5*y*(-20050351692839232072401 + 120*y*(-4976219400498134306 + 701059526547943945*y)))))) + 6*(x)^(2)*(-13034991506843954692628435037 + 40*y*(-12911412490822047011077716 + y*(71942117064863810449017711 + 800*y*(3016805471122684680351 + y*(-7306648318334971703039 + 8*y*(-22005304021302476127 + 13952322443747600185*y)))))) + 6*x*(3583617805503366573901766935 + y*(-539997863237214241774925283 + 40*y*(-32876604623363072658989403 + 2*y*(3430167819779391335035679 + 240*y*(5681012924106078824787 + y*(-1148410005247339716587 + 200*y*(-521520979887370115 + 80705670151250134*y)))))))))/1.5123064061952e29",2); // 2
        iterLimit=5;
        break;
    case 12:
        ss << "4-patch domain, degree 4, v2\n";
        pathStart = "Mario/Florian_FourPatches_deg4_v2/multipatch_four";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        startIter=1;
        g=gsFunctionExpr<real_t>("((856 + 25*x - 210*y)^(2)*(794 + 185*x - 15*y)^(2)*(253 - 60*x + 5*y)^(2)*(819 + 190*x + 20*y)^(2)*(-772 + 185*x + 25*y)^(2)*(343 - 10*x + 80*y)^(2)*(507 + 10*x + 120*y)^(2)*(-1447 + 60*x + 350*y)^(2))/7.001316e42",2);
        dg=gsFunctionExpr<real_t>("((856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-183266112727112768808+93639600000000*(x)^(7)+388500000*(x)^(6)*(5638646+318315*y)-150000*(x)^(5)*(1137334053991+5*y*(38770857663+14530646650*y))+12500*(x)^(4)*(-152716930485143+25*y*(-156401552027+349513306994*y+7224284480*(y)^(2)))+1000*(x)^(3)*(70984719289425039+5*y*(-118089672499773+50*y*(-29650643938913+293423504802*y+882217595180*(y)^(2))))+75*(x)^(2)*(677492733286137563+5*y*(11600640434709637+10*y*(-981143799014747+10*y*(-12665036490069+50*y*(20291359787+6678463216*y)))))-5*y*(4654806746679084698+5*y*(-595081862571168173+5*y*(-21886279724378373+10*y*(156538493064151+160*y*(377450354536+45*y*(156065933+4364500*y))))))-20*x*(61409432176835309771+5*y*(-231054258707148997+5*y*(-273537297379319971+5*y*(686445838949453+10*y*(155942373681973+100*y*(-5142492283+817279476*y))))))))/7.001316e41",
                                   "((856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-24*(5557475750222879756+52643682537326664905*y)+5*(3533296500000*(x)^(7)-25000*(x)^(6)*(38770857663+29061293300*y)+12500*(x)^(5)*(-156401552027+4*y*(174756653497+5418213360*y))+250*(x)^(4)*(-118089672499773+100*y*(-29650643938913+y*(440135257203+1764435190360*y)))+25*(x)^(3)*(11600640434709637+20*y*(-981143799014747+5*y*(-37995109470207+200*y*(20291359787+8348079020*y))))+30*(y)^(2)*(199526503227539871+50*y*(9930626942414458+5*y*(-17989001088697+24*y*(-141086396887+122500*y*(2701+384*y)))))-10*(x)^(2)*(-231054258707148997+5*y*(-547074594758639942+5*y*(2059337516848359+40*y*(155942373681973+25*y*(-25712461415+4903676856*y)))))-x*(4654806746679084698+5*y*(-1190163725142336346+5*y*(-65658839173135119+40*y*(156538493064151+400*y*(188725177268+9*y*(468197799+15275750*y)))))))))/7.001316e41",2);
        ddg=gsFunctionExpr<real_t>("((819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(4*(-3+100*x-820*y)*(167+4440*x-365*y)*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)+(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(-5872151+15000*(x)^(2)-300*x*(3+820*y)+80*y*(949+12605*y))+(856+25*x-210*y)^(2)*(343-10*x+80*y)^(2)*(-178355327+29570400*(x)^(2)-13320*x*(-167+365*y)+5*y*(-55462+39965*y)))+(856+25*x-210*y)^(2)*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(343-10*x+80*y)^(2)*(40*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(319+24*x+214*y)*(-1447+60*x+350*y)*(967+14060*x+1690*y)+(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(-1777002527+296525400*(x)^(2)+42180*x*(967+1690*y)+60*y*(290448+71035*y))+4*(819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(-6259523+21600*(x)^(2)+1800*x*(319+214*y)+20*y*(172951+82445*y)))+20*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(-9685978+444000*(x)^(3)+(x)^(2)*(5070-5515950*y)+y*(54496707+25*(13225-49548*y)*y)+2*x*(-132371273+5*y*(287721+1551850*y)))*(-2726354163+16872000*(x)^(3)+300*(x)^(2)*(1127087+762350*y)+10*y*(-257314061+40*y*(148557+180125*y))+40*x*(-260893019+5*y*(646909+3134030*y))))/1.4002632e41",
                                    "((819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(200*(-71+82*x-672*y)*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(7+73*x-6*y)*(-253+60*x-5*y)+(856+25*x-210*y)^(2)*(343-10*x+80*y)^(2)*(-1204067+199825*(x)^(2)+x*(30560-32850*y)+450*y*(-7+3*y))+4*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(-98526263+252100*(x)^(2)+50400*y*(71+336*y)-20*x*(14807+206640*y)))+(856+25*x-210*y)^(2)*(794+185*x-15*y)^(2)*(253-60*x+5*y)^(2)*(343-10*x+80*y)^(2)*(8*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(1007+1690*x+200*y)*(-1447+60*x+350*y)*(381+1070*x+8400*y)+(507+10*x+120*y)^(2)*(-1447+60*x+350*y)^(2)*(-24276671+4262100*(x)^(2)+600*y*(1007+100*y)+60*x*(59951+16900*y))+4*(819+190*x+20*y)^(2)*(-772+185*x+25*y)^(2)*(1648900*(x)^(2)+180*x*(78963+149800*y)+27*(-22818637+2800*y*(127+1400*y))))-20*(856+25*x-210*y)*(-343+10*x-80*y)*(794+185*x-15*y)*(-253+60*x-5*y)*(819+190*x+20*y)*(-772+185*x+25*y)*(507+10*x+120*y)*(-1447+60*x+350*y)*(x*(-54496707+5*x*(-287721+367730*x))-250*x*(2645+62074*x)*y+450*(713+8258*x)*(y)^(2)-201600*(y)^(3)+4*(6617497+67949189*y))*(-3*(406850873+3587051620*y)+10*(7623500*(x)^(3)+7500*(y)^(2)*(1707+224*y)+10*(x)^(2)*(646909+6268060*y)+x*(-257314061+120*y*(99038+180125*y)))))/1.4002632e41",
                                    "(8736*(-5889591170879306018956271073243080922+22925776618634104814916613389531613295*y)+25*(62035588301512500000000000*(x)^(14)-4312350000000000*(x)^(13)*(3324258357687+2740557975655*y)-1803750000000*(x)^(12)*(307513783018520481+110*y*(1043969817180299+105368781544875*y))+18750000000*(x)^(11)*(1616209809586420291851+y*(1238021454387970462733+2*y*(183253666887968429619+58954243271100176660*y)))+137500000*(x)^(10)*(3860051208614600588030816+5*y*(-305396154884365519380817+50*y*(-2492215107508311154383+20*y*(-30895915857524455283+324570185189587235*y))))-6250000*(x)^(9)*(2469981477768200761696850107+5*y*(-170548179992646678377359407+10*y*(-15860356643192086745823063+20*y*(-245996004761701292488629+100*y*(462499575852494684963+191973545558897389557*y)))))+5625000*(x)^(8)*(-7933898947401529305690216499+10*y*(7139097275043354971244075328+y*(236968812099024590167053471+20*y*(-38254237482625477235750804+25*y*(-28832614495174597470411+8*y*(5612513556492807328773+96858236340211903900*y))))))+1250000*(x)^(7)*(-90008106283614443392212407123+y*(-8432585990927600005430502072037+2*y*(73725322169426088714622191129+20*y*(33267867471633597699708373449+50*y*(-11859166675181326539831076+y*(-37751660152923343565540931+10*y*(49339595481767764811939+78022255196391533949800*y)))))))+43750*(x)^(6)*(84576571841663318667353827852358+5*y*(-119270021041172613069243491657009+2*y*(-5267006730739703066455523697141+20*y*(365096653088222551419951669109+50*y*(720304994113427524676916996+y*(-262145884199479055265345723+10*y*(-3766966184007510154199879+40*y*(4858565047001895189461+1654278215577084353790*y))))))))-9375*(x)^(5)*(-3944665406579716882953309808152639+y*(-59631788330605795969226108804524601+2*y*(1175294921656112003132297904715299+20*y*(247260307856921203334807082003489+50*y*(-133087026003182851788371621127+y*(-280966662222325225695067978251+10*y*(462831431276262238433359389+1000*y*(543109951341832710961281+5*y*(-445873082654359563939+60986370928455606104*y)))))))))-1250*(x)^(4)*(86667944396611326125152191968589291+5*y*(-92711306107668717189746378411325493+2*y*(-6137440113826370431243828948588794+5*y*(1221053113991921726225936571149394+25*y*(7199094000253875689933523616543+4*y*(-487098523544301830556854507817+5*y*(-19400717898846330756546159071+20*y*(48841879472535873837618516+5*y*(3460946741435033281024677+200*y*(102652295928016054801+11312713918519394172*y))))))))))+20*(x)^(3)*(-43716237256277211654117871080508759761+5*y*(-100504857453883898313383545165878863533+10*y*(494455112403789248890687184221654578+5*y*(339712344190683828714051811107583207+25*y*(-435352475067945137636850632278106+y*(-775588456340703311348419990354869+10*y*(1551638818925981834611256221577+40*y*(37254802836739713592089212011+60*y*(-5055376602138171854610588+25*y*(5433891171018067058939+4*y*(-30756759099767012807+2078936708753422440*y)))))))))))-3*(y)^(2)*(-517187094791941512207875623483647488772+5*y*(77217754052102660111443949529475785072+25*y*(633448828506561189543976944087298763+5*y*(-27988221239870592477014513717078484+y*(-7102809297961140075487637866868499+20*y*(17876842976762970788426265473984+5*y*(1368773374720284161388360799551+200*y*(64129244261551769175819087+32*y*(-194672732009669889214129+90*y*(-68395720077821388488+13125*y*(-68241190290897+784*y*(2293809571+56115000*y))))))))))))+3*(x)^(2)*(357008814843850566704563621317729509638+5*y*(-316328078054897574117811933975431492486+5*y*(-10541312754448979091768159561196028679+5*y*(1754523182153765065935875096330768944+25*y*(12679802486482816224180486557374431+10*y*(-300432851590565251710084441841767+2*y*(-34770495496633093130101989295657+40*y*(43348784030460636922030771701+5*y*(3181431815859714378054765579+200*y*(131321929928275706385031+18*y*(273910185345199868311+240*y*(31713612288757402+1497839540845375*y))))))))))))+x*(5975401493509874489718402548998234307088-25*y*(-2425969083067941644958089230894270786564+y*(134823314409260395756442128075036963626+y*(415920671685559246029616085390606549878+25*y*(-600479315071868234907855631100928349+y*(-961533574482894251649317866257559143+10*y*(2260764901252985584405566308733511+200*y*(9504643318524626632512636282979+y*(-116341862738567034959898413679+800*y*(-9882119516256426019546181+6*y*(-2056855346844450996833+180*y*(-6902072649829414053+1750*y*(-12582105727243+1373281449960*y)))))))))))))))/7.001316e40",2);
        f=gsFunctionExpr<real_t>("(200483015939927424805749401753922095656048-20792963040042998776819931131812603569760*y+5*(85796184373560849187500000000*(x)^(12)-7500000000*(x)^(11)*(188763309224212970923+5082269320562618300*y)-1875000000*(x)^(10)*(-69162007436240087786159+8*y*(2773202004406440477929+3894193987465008424545*y))+25000000*(x)^(9)*(-149909623908054866557383907+400*y*(-63214168016399689867944+5*y*(14015207218248956806026+598454433048046209625*y)))+7500000*(x)^(8)*(19473177179241234640028364387+50*y*(31862975136666387602469487+10*y*(-16355455723645337460846387+10*y*(-112982613228639206194+109202032280810305350935*y))))+750000*(x)^(7)*(1572537911999558169581729836189+20*y*(12594309484008206716874150172+5*y*(-4165815250858270258696456113+100*y*(-4733092593717614720471873+40*y*(40903113642087276773771+5340069121376618277447*y)))))+12500*(x)^(6)*(-2271988464656030168577251066789687+20*y*(2944743138851483077922538850602+5*y*(9910947021964690206354762839607+700*y*(-292940854633575377776655719+40*y*(-28505650653407207205195069+279668611635311480169144*y+503637821082003617841650*(y)^(2))))))+7500*(x)^(5)*(-11404202153620510282266286367953621+40*y*(-41430614895086613420012978375721+10*y*(7300076807925748182243512523507+5*y*(207798978778270540024515420577+5*y*(-15704605651331332650308349669+20*y*(-131192453670640178453814723+100*y*(90447654906181480116769+18546532481001706820208*y)))))))+375*(x)^(4)*(3679660093973470770621986393970644679+100*y*(-2109005546190532388780547854257611+5*y*(-2849431322847162531737397627418393+20*y*(4557430842446960617713592202719+5*y*(2392150695428852486560506087695+4*y*(-10301361341987745035225899043+50*y*(-276806349418032697590524561+100*y*(18342483418404099256747+12322864260311585068187*y))))))))+100*(x)^(3)*(15056257746185269082445666631937465423+125*y*(20980001691361741808527279021965765+2*y*(-15841241175542331931316536052702589+10*y*(-285162413610154714654165106097818+15*y*(6409883038339352313893789206357+4*y*(333508619420475236878379216642+25*y*(-1023407103086551980758209087+20*y*(-12877828294431805893821999+5*y*(22595898331228390795931+8841866686884050087760*y)))))))))+5*(y)^(2)*(-3551531163988261907495742883855209660492+5*y*(34381021327838633385934483932486896924+5*y*(11829420487260040622572021242890193237+20*y*(-16049295832071177139434109298296377+5*y*(-2419016927062141523648868713984609+300*y*(122010981956094749795023301659+y*(44005037843774091904242393211+80*y*(-3472649528000044723949072+15*y*(7892045651005192672523+40*y*(-2397054208507204837+142550086209543870*y))))))))))-30*(x)^(2)*(565149859396440035224265578183981607208+5*y*(-8437287002057048996772266807335117994+5*y*(-8716574128893158922537895479104568879+50*y*(6918879623452658688632669271673658+5*y*(2958282868548894096977145766528735+4*y*(-16560015850331253540195841546284+5*y*(-3557510113118293272483686100191+60*y*(652163400886916781761946087+125*y*(2624182264153707638925741+8*y*(-999457787323601761019+135161035690724068988*y))))))))))-12*x*(530925557520583884799110854193758960732+5*y*(24376442401877963241185423320934974674+5*y*(-5999125722417512746535516369169733891+25*y*(-54851286351950879916482325313472751+5*y*(3098631074456770371333310355311893+8*y*(98897719114710296309142369775414+5*y*(-1337911001230730072008322933889+10*y*(-40430425105286506983844478437+250*y*(1301566201057324263051717+2*y*(363125623376994231239461+4*y*(377139568951883178101+44110620057322458600*y)))))))))))))/1.4002632e39",2);
        iterLimit=5;
        break;
    case 13:
        ss << "5-patch domain, degree 4, v2\n";
        pathStart = "Mario/Florian_FivePatches_deg4_v2/multipatch_five";
        pathEnd = "_degree4.xml";
        mapPathEnd = "_map_degree4.xml";
        startIter=1;
        g=gsFunctionExpr<real_t>("((51.2 - (32*x)/3. - y)^2*(-17.9 - (23*x)/5. - y)^2*(-6.855555555555555 - (13*x)/9. - y)^2*(3.7956521739130435 - (10*x)/23. - y)^2*(-4.41578947368421 - (2*x)/19. - y)^2*(3.0448275862068965 + x/29. - y)^(2)*(4.45 + x/2. - y)^(2)*(-4.923529411764706 + (9*x)/17. - y)^2*(6.1 + x - y)^(2)*(-21.6 + (9*x)/2. - y)^2)/1.e16",2); // 3
        dg=gsFunctionExpr<real_t>("((883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(5184*(-24+5*x)*(-60759448049004644133+10*x*(9224871288649613571+100*x*(37604179460108437+10*x*(-85805440864572+5*x*(-23585288349969+2*x*(-344361452181+40*x*(2225326239+88527705*x+747500*(x)^(2))))))))-135*(17023674641706872782653+20*x*(-1695509837688518024354+5*x*(-41960889498677262899+100*x*(284438575105026496+5*x*(2833672964522191+12*x*(-96263698034234+5*x*(-270006365587+80*x*(969070628+45871135*x))))))))*y-250*(3209323891457105431287+4*x*(-1063905721177541644621+20*x*(-10863701126609174361+5*x*(311839709514378744+25*x*(4409023124019543+4*x*(-39053380365891+50*x*(-47676560757+30484236416*x)))))))*(y)^(2)+500*(414137231579450182133+100*x*(-8503573356019086696+5*x*(-169890059042411955+4*x*(19655998498441548+25*x*(66279852581471+60*x*(-104658896724+22349641013*x))))))*(y)^(3)+25000*(669962652693164829+4*x*(-252981723406319757+20*x*(-3480355330614159+5*x*(-149651754464292+25*x*(1579945392551+246871418796*x)))))*(y)^(4)-50000*(36497337774982877+20*x*(-13208323166449558+15*x*(-26473933807981+40*x*(345529412146+198697687975*x))))*(y)^(5)-500000*(-1618566060310711+60*x*(-11342941558297+40*x*(-39209104117+12314778830*x)))*(y)^(6)+15000000*(-8007076526467+20*x*(-437203258412+139024822065*x))*(y)^(7)+50000000*(-475397417077+52105762340*x)*(y)^(8)-621668112000000000*(y)^(9)))/1.69182035233245e49",
                                   "-((883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(9*(-24+5*x)*(-16077059186642225498133+4*x*(-1369336665190955685819+5*x*(48913670472575889613+20*x*(946609999700553672+5*x*(-5001444039303987+100*x*(-28130131110147+x*(156037145859+40*x*(2620771809+91742270*x))))))))+20*(145441291720013932895181+5*x*(3209323891457105431287+2*x*(-1063905721177541644621+40*x*(-3621233708869724787+5*x*(77959927378594686+5*x*(4409023124019543+10*x*(-13017793455297+100*x*(-6810937251+3810529552*x))))))))*y-300*(4238909567658291807147+x*(414137231579450182133+100*x*(-4251786678009543348+5*x*(-56630019680803985+4*x*(4913999624610387+5*x*(66279852581471+300*x*(-17443149454+3192805859*x)))))))*(y)^(2)-4000*(75855763148894346513+5*x*(669962652693164829+2*x*(-252981723406319757+40*x*(-1160118443538053+5*x*(-37412938616073+5*x*(1579945392551+205726182330*x))))))*(y)^(3)+50000*(1644683832982118973+x*(36497337774982877+20*x*(-6604161583224779+5*x*(-26473933807981+60*x*(172764706073+79479075190*x)))))*(y)^(4)+600000*(18992823414703953+x*(-1618566060310711+10*x*(-34028824674891+40*x*(-78418208234+18472168245*x))))*(y)^(5)-21000000*(64762860701939+x*(-8007076526467+20*x*(-218601629206+46341607355*x)))*(y)^(6)-80000000*(1985626097571+x*(-475397417077+26052881170*x))*(y)^(7)+2700000000*(-500083579+414445408*x)*(y)^(8)+116338140000000000*(y)^(9)))/3.3836407046649e48",2);
        ddg=gsFunctionExpr<real_t>("(80*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(9*(-543914379+1936272910*y)+100*(162*x*(-1003582+5*x*(-993+3568*x+50*(x)^(2)))-18*x*(-3185437+643605*x+61400*(x)^(2))*y+(-59423569+20*x*(601672+190725*x))*(y)^(2)-50*(25691+84020*x)*(y)^(3)+1294400*(y)^(4)))*(2619201987768-935774142455*y+50*(16*x*(-2700305626+x*(-469551171+80*x*(1108018+37375*x)))+x*(-11430373379+20*x*(132519969+24500680*x))*y+4*(-1888260208+5*x*(102280583+48691290*x))*(y)^(2)+10*(2825021+47307530*x)*(y)^(3)+50886600*(y)^(4)))+25*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*((883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(837-90*x+170*y)^(2)*(81*(-5687+60*x*(13+10*x))+180*(167-330*x)*y+15700*(y)^(2))+8*(883+10*x-290*y)*(-837+90*x-170*y)*(117+180*x-110*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(9*(-11809+60*x*(293+5*x))-200*(1485+296*x)*y+104900*(y)^(2))+4*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(27*(-38127151491+30252883840*y)+40*(243*x*(-8844351+5*x*(749023+50*x*(586+5*x)))-360*x*(7615641+10*x*(217413+3700*x))*y+10*(-75702173+60*x*(4329504+266245*x))*(y)^(2)-300*(7743339+1596890*x)*(y)^(3)+421028250*(y)^(4))))+4*(883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(837-90*x+170*y)^(2)*(20*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(-3344+7360*x+1145*y)*(-5422463+240*x*(8226+325*x)+2406380*y+649600*x*y+780500*(y)^(2))+(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*(256*(-1932479+1380*x*(-209+230*x))+320*(-81373+79005*x)*y+1863025*(y)^(2))+25*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(-59816709411407-65166891181640*y+40*(6*x*(-163190328849+20*x*(865547273+6500*x*(16452+325*x)))+20*x*(3646882989+40*x*(143693421+5278000*x))*y+5*(-325804009+360*x*(135707027+12173420*x))*(y)^(2)+300*(536897933+131866100*x)*(y)^(3)+21616698250*(y)^(4)))))/1.69182035233245e49",
                                    "(-40*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(179+46*x+10*y)*(-768+160*x+15*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(8106539481+9*x*(-193627291+10*x*(-3185437+10*x*(42907+3070*x)))+171788588*y-20*x*(-59423569+20*x*(300836+63575*x))*y+60*(-19223711+25*x*(25691+42010*x))*(y)^(2)-800*(-127477+64720*x)*(y)^(3)+9860000*(y)^(4))*(1804179436833-1545589725720*y+4*(x*(-187154828491+5*x*(-11430373379+40*x*(44173323+6125170*x)))+40*x*(-3776520416+5*x*(102280583+32460860*x))*y+300*(-347459898+x*(2825021+23653765*x))*(y)^(2)+6000*(-507835+339244*x)*(y)^(3)+147487500*(y)^(4)))+4*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*((883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(837-90*x+170*y)^(2)*(-2327+3925*(x)^(2)+300*y*(31+2*y)-20*x*(794+165*y))+400*(883+10*x-290*y)*(-837+90*x-170*y)*(89+10*x-20*y)*(61+10*x-10*y)*(-216+45*x-10*y)*(-31+11*x-4*y)*(-115123+1480*(x)^(2)+x*(14850-10490*y)+3*y*(-8451+4930*y))+2*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(983008657193+3696921147420*y+20*(26624500*(x)^(4)+600*(x)^(3)*(1443168-798445*y)+15*(y)^(2)*(-3897676511+98600*y*(-8451+2465*y))+10*(x)^(2)*(-75702173+90*y*(-7743339+2806855*y))+x*(-70039005912+10*y*(2650634259+10*(167638491-51715700*y)*y)))))+(883+10*x-290*y)^(2)*(89+10*x-20*y)^(2)*(61+10*x-10*y)^(2)*(216-45*x+10*y)^(2)*(837-90*x+170*y)^(2)*(8*(179+46*x+10*y)*(-768+160*x+15*y)*(-999+458*x+60*y)*(617+130*x+90*y)*(839+20*x+190*y)*(-873+100*x+230*y)*(-4919953+324800*(x)^(2)+20*y*(294019+58995*y)+20*x*(120319+78050*y))+(617+130*x+90*y)^(2)*(839+20*x+190*y)^(2)*(-873+100*x+230*y)^(2)*(298084*(x)^(2)+60*x*(-16589+1374*y)+9*(-72407+60*y*(-333+10*y)))+4*(179+46*x+10*y)^(2)*(-768+160*x+15*y)^(2)*(-241540077242153-193437613834440*y+40*(3652026000*(x)^(4)+600*(x)^(3)*(135707027+65933050*y)+15*(y)^(2)*(47746822063+65550*y*(588038+58995*y))+5*(x)^(2)*(-325804009+60*y*(1610693799+432333965*y))+x*(-3152734415003+30*y*(5512681989+10*y*(3241247549+511617750*y)))))))/6.7672814093298e47",
                                    "(-11664*(24-5*x)^(2)*(3795508851241428081609835296757521036403467+8*x*(-617924889466712236003969418712490006166196+5*x*(-117799961390115888244607713916466201335763+10*x*(-1153754925035254244930228232162930630276+5*x*(190783077681420820239961284719819672741+400*x*(94136542817661227100453777189555504+5*x*(-2335776705045331973430661489837887+5*x*(-139662379522414132011851610140952+x*(2989945885086361921720201255521+40*x*(55027619254438667298341299236+5*x*(193110516368558089453986589+10*x*(-7645478313504331939059012+5*x*(-65763204980290218177741+160*x*(19183304824653527736+5*x*(350631176011193151+40*x*(197697136981902+1302969589675*x))))))))))))))))-3645*(-24+5*x)*(-133217222313490691110038761522914324054767573+4*x*(30707381976125011223697233902455332398970319+20*x*(1191165278871544872251481597207897295736307+25*x*(6059617791835945705996637645576510043311+5*x*(-587254050168334401901216944803729861137+80*x*(-2857926878484637578930592938424051364+5*x*(-1924196692092154206887727370076932+x*(24731098175816427171948650852354142+5*x*(298808202379124496086373676196557+20*x*(-4061264643115522705627614150193+20*x*(-9322530128165611705629893263+5*x*(254749554884573818720518079+5*x*(-2504987654656513967213119+40*x*(-3735831775608908438463+80*x*(11411323124051659378+15*x*(23233173342897719+23271419454370*x))))))))))))))))*y+675*(25695529994036350691195886821931982550755068915+2*x*(-21222561349774070294281120924390666690935287463+5*x*(-1470192296359704686899216335439906402872612009+20*x*(32828438928021100531095137055388422092801442+25*x*(375771878930388547769390107169649955359629+8*x*(-9297559430780970006311141091744818659836+5*x*(-518824950630865729493051643271960121031+40*x*(1239337523806790172465481593355752474+5*x*(70298941463740785448635833708320698+25*x*(-165266538000479350228250648119671+x*(-34434605624479501963715560448477+20*x*(127192476002036628020827168206+5*x*(1247946149098604414584565817+200*x*(-1222897611049032857462391+10*x*(1161178534000944312057+320*x*(1266531504681501663+41525418232972513*x))))))))))))))))*(y)^(2)+1000*(454908891261705405293290728846770879955211962+x*(1192377384387916217783146947465730421600899239+10*x*(4323767282325281414116112562531325444590531+20*x*(-3432367826491467794474142636672332254331339+200*x*(-829274239270937985012716379428979999543+x*(1866368217671403190852503064272450591753+5*x*(14660694531649725548851755369383919711+20*x*(-814066796197809791333767300226627538+25*x*(-845706731082957996585568457594364+5*x*(136341428904494516644467134448151+2*x*(-523686404450488294038346109853+100*x*(-6401149892629696416060739017+20*x*(17276974088396173674640689+x*(893137738722155892236061+100*x*(-308595463037841925575+35323694640674483104*x)))))))))))))))*(y)^(3)-12500*(203547927545787429973647834403689761163930297+4*x*(-88700405875258812633570461844349224850692594+5*x*(-5088884454361057926285129444449369229933843+20*x*(102037029733931730830901970205868315022536+25*x*(1139056577223180365690095879519796153237+20*x*(-5997254477044199421896067464297433513+4*x*(-675585137903823909090519779730570562+25*x*(-395141985092314217490706644985944+x*(588273024309087438974071601399403+100*x*(304204911714318416453141128308+x*(-62276343587707471314832158463+20*x*(-20909933786942601669098856+5*x*(2866213319543382630659383+60*x*(512820645506886052273+550407119691908322750*x))))))))))))))*(y)^(4)-150000*(-820912876560040055454432453797226841034706+5*x*(592450044345189421336013296999840871075323+80*x*(1311506726167945985871175023744316137447+x*(-1803079355393752025371711836650070599656+25*x*(-7267115612631854965848892570016054753+x*(4807123355570829949579856201052443019+80*x*(6152635135507867630226220089665361+5*x*(-309672464598060704832373938469532+5*x*(-7725830568381436150542669393672+25*x*(30980980740115059043000786097+16*x*(172350682542795431883899746+25*x*(-796881288663087402354876+5*x*(7391466654891332268515+236379089094056247129*x)))))))))))))*(y)^(5)+1750000*(80773870313750739897812091863957297184759+2*x*(-92643044110063114738520207042823764986917+25*x*(-719117177967701738403704851801130970015+4*x*(77634810668446222272798278038044225162+25*x*(766259340515986553579744402873911759+16*x*(-1214932437321533156667985363337886+125*x*(-15564052439637071236626757180113+2*x*(-844742087491991094270169700424+x*(80097624826001896502857819221+50*x*(405967284015213374071442325+x*(-10201816907430873985157249+60*x*(12140065785830258340462+8268842088053710100275*x))))))))))))*(y)^(6)-20000000*(539426160540983772054245755685428914630+x*(-1790941405357047567376315431779649452431+10*x*(-18091672375436750457831888150781079637+20*x*(1591626137629073876880994869869520241+400*x*(263805077291529704155031159940002+x*(-168796794272635732890337453811646+5*x*(-4533024644089582279726522153397+50*x*(6604069465120570812393238904+x*(2427734883696585912749793852+25*x*(5319031552190841163885487-970751473560998415598774*x+54273906819203417628840*(x)^(2)))))))))))*(y)^(7)-225000000*(8279588480765903113524666213672961341+4*x*(-14802325646602365620237113770635840653+40*x*(-24629286857699112254909947458500148+5*x*(4750519337020804929328849712882428+25*x*(34961665971044349794100633258143+4*x*(-97995562129368317132554465683+10*x*(-34798545644587604875845990133+15*x*(-276055521378564648930042248+5*x*(-3817133550493210740612501+100*x*(13829622398451154550069+2124211143820855593988*x))))))))))*(y)^(8)+2500000000*(138552077559380492720100083552259186+x*(-518129521994851713755497848220513949+100*x*(211060585084324952228717340164817+4*x*(181997514524528155615390486535897+25*x*(-191528781724603748936085265405+x*(-218178188813585401598670114969+80*x*(-316192049271703378831743693+10*x*(3242759559694185931373896+5*x*(43148931129595775808171+31046168331245185047725*x)))))))))*(y)^(9)+27500000000*(-4170899306537626970232096726120063+2*x*(-5594389308271109900804197339653193+25*x*(22937196753854554456622278215477+4*x*(2782364666669841709901242202438+25*x*(12655247052698824086276011299+8*x*(-362306615452479131928662592+5*x*(-21680317993242996375333143+180*x*(799650858117221408548+2385269629366586971095*x))))))))*(y)^(10)-300000000000*(32197381414060723715679384929442+x*(-74083180728190222571633627817181+10*x*(1787978184548109081564181844403+20*x*(35839990317274869045377525197+200*x*(-38454844495213594178907107+x*(-3614482475652951149560803+5*x*(-92854244493612083309023+50260732806860363055480*x)))))))*(y)^(11)-3250000000000*(-1386001602410594435486708400885+4*x*(-255050973614285344755027951184+5*x*(19773915409078463607981525441+20*x*(16772096711582432292766752+25*x*(-218069407031470220897891+12*x*(-18393983188864256510007+2078881834840572039740*x))))))*(y)^(12)+35000000000000*(9527124293161804308265821318+x*(-5150850242757708061721573957+40*x*(48233241498574572082579983+10*x*(-54691291943440779941534+25*x*(-15874909021313338081109+2545236255807717734007*x)))))*(y)^(13)+1125000000000000*(-37864715645896890114243841+2*x*(-4200264675491427861837169+5*x*(1072372985577804280726881+20*x*(-11343858139877762710678+726455486755366925275*x))))*(y)^(14)-4000000000000000*(1240395513980207205085266+x*(-341925624484588249279597+250*x*(-92745212797983746355+48508078875627116716*x)))*(y)^(15)-127500000000000000*(967519220929799627413-521836616853940347052*x+64740269323974620280*(x)^(2))*(y)^(16)+8100000000000000000*(203631613331080255+52205345473267412*x)*(y)^(17)+68707526255022096000000000000000000*(y)^(18))/1.69182035233245e48",2);
        f=gsFunctionExpr<real_t>("(243*(50497052605010560134063815435000252508733486047-199014893481648118666510199942247945744424463900*y)+20*(3*x*(212407118668070050400837154153816241340492165898+5*x*(-18417110037621892049763610258305244835011844133+20*x*(-680471534018737869090296316337848469115026083+80*x*(1175459346529485922139611944122104414210906+5*x*(172112960799268614739693592859142901011962+5*x*(-1716048929701797561985523204366578140681+25*x*(-61540159884695710531482095952281403942+x*(844705254259547966483851405161850047+100*x*(13311970415181511628653825075630350+x*(21645699412432041655030946513971+20*x*(-6840738135683690390772677187123+10*x*(-5518070219015182905821289939+100*x*(27808905458626020861295272+x*(366953170538786299641243+40*x*(-348306547122340607949+11916893274269801150*x)))))))))))))))-30*x*(33334846578718803712836745187308061439617260371+20*x*(-1235062268868930927860816357677686404097882117+100*x*(-4162363803420690873199447358185652490590699+20*x*(59354126418490270907082165869470140161808+x*(17547035779165263316468961476545420914483+5*x*(-494316158274430219657562901683749179639+10*x*(-13028169299228303533675304608094508296+25*x*(40844288345705546058646593265867851+10*x*(877548408929790707935605862947371+4*x*(-11672414521425367319178577561041+20*x*(-78110130302860872116052661922+125*x*(30690393923047724427250491+2*x*(795065881948686441543877+30*x*(-967170707027605271513+16218672343410725572*x))))))))))))))*y+180*(546784633884459797216264789039179116365722161+20*x*(-1261343022606605000579609325372221472989979+5*x*(-7910031470487984440969074309689579727435089+100*x*(249456105019200303851154822522557683342+5*x*(2658778978542364937119463883012840594164+x*(11837726480363995088877802606471616078+5*x*(-30236551928764998335683280960199343251+400*x*(-577270751437851340335378275133493+50*x*(34382298422256974073635525703362+x*(-157782044335771471639308494798+x*(-281100749404461107701887202019+200*x*(73673101806336391876342283+x*(2883637905851808887742981+25*x*(-18997479527234621613806+528402799876198962725*x))))))))))))))*(y)^(2)+100*(9017290723668017476849352046018796366923229347+250*x*(9492778262890475574690917992561811624648949+2*x*(-3672016492380678559820549609323939052626833+10*x*(-98747911651863524864358293730517276615679+10*x*(1883360543722071961245940708117682150673+2*x*(348889663344272393905935290004078564951+100*x*(-83055242212311260548420447421567781+10*x*(-9829995116001099856347546784538811+5*x*(-42441595642999529627470808931927+10*x*(1989955793891110454441595706607+30*x*(711566141441201205351729+350*x*(-731935742391105870558027+51939955730717659255126*x+739754669154418789140*(x)^(2)))))))))))))*(y)^(3)-3000*(21673824097841548582257881194967159573305707+100*x*(67307062718311535701271568890582266771212+5*x*(-27087802405204651505569485170902870948459+10*x*(-259091505734641385884563506842733972788+5*x*(57715317851430521102905951386608284001+40*x*(115217684517907150380276089395585652+5*x*(-10349458004456435064473344099802011+20*x*(-47464361636532608838961727858717+25*x*(291261872814751752245258403545+8*x*(2661511493441287143316932561+x*(-189918650537643461574473239+25*x*(-61934189881380958360838+30952138707475426534555*x))))))))))))*(y)^(4)+150000*(-755560363725838359148560498360055616719431+10*x*(-9398389697373280484814341607493095354511+8*x*(1415094952708344193780829970553508062152+5*x*(48606220060322909274751192041254372533+240*x*(-27043157640897968667989883659190327+5*x*(-2782368102220777091017371419980194+5*x*(-43768869279446291152036106405198+5*x*(2205506363056528264195298327218+5*x*(61168109026728820251894542829+70*x*(-16747992317543951268646941+40*x*(-96319517268418945191139+271324349632049063355*x)))))))))))*(y)^(5)+200000*(29104182198910131515974604531722603398751+5*x*(1597458383515281765931390171701584744226+5*x*(-665921640312003440225668625532789833643+400*x*(-73846786692495739528068979647246829+5*x*(24170024772592512246368424419143773+20*x*(96536989159256977948711411909977+10*x*(-2280584745509951477934817897207+2*x*(-223890845213618597042804274732+25*x*(66422983416824852434601577+5*x*(13567676776997201615030462+627502089657674564902863*x))))))))))*(y)^(6)-3000000*(-2183141472399162879850734636193543514253+10*x*(5203847932316006680924487269794620689+60*x*(425753650834482140037016531334563507+50*x*(398272255709993952277980984515093+10*x*(-11159307986682461224941277919733+2*x*(-2740157694693315875799124637527+10*x*(-24112942836407663765239735167+10*x*(208477385435870197417335587+25*x*(2707721917645645200751029+117992421195564066894710*x)))))))))*(y)^(7)+60000000*(-1858439803665643870582093506184009744+25*x*(-29764471284054843178141569434243102+5*x*(15236519769154931499868486299686379+20*x*(-52362908417205495508877028176011+10*x*(-3641558368086216508169976721159+12*x*(-4920755046711964831324968936+125*x*(26616532254093248644321669+2*x*(2536617248288587614014584+69553362429410636411685*x))))))))*(y)^(8)-100000000*(1814384361101580973443885929252954683+50*x*(-6787952704150918425822996742550169+20*x*(-155913858053357011016072814192261+20*x*(857024383932366856343012817511+180*x*(411909157530444376154960557+x*(146029709080024790265934293+125*x*(-71501768135188090199793+497222670838439389888*x)))))))*(y)^(9)-6000000000*(612080638310556845374925633967631+10*x*(-27057853148622611920984009081132+5*x*(5488885379474337745883239715903+40*x*(-29852035998534663201955877811+35*x*(-69102380756204146241290989+5*x*(2444066777750734913761658+126567460808479706709525*x))))))*(y)^(10)+30000000000*(66035767323965342858727737494787+10*x*(-2515527471708063637361944555387+30*x*(-7813199029223015950393815429+50*x*(90692135802509159620384297+30*x*(-191718947417303873969225+9917257402117231729682*x)))))*(y)^(11)+100000000000*(1068028978938837866922026856161+40*x*(-11388528454603134464531889237+35*x*(72287855979385413512289024+25*x*(-625337140047818638548434+42753880513474669329135*x))))*(y)^(12)+3000000000000*(-286768377678985797203887497+50*x*(9139797036319291657602643+120*x*(-19308984937116675326462+438275889262669832745*x)))*(y)^(13)+60000000000000*(29731229844053438771703+25*x*(9231283041521873999582+1049228045018022708195*x))*(y)^(14)+300000000000000*(11704072636745431716719+1160590434881456227170*x)*(y)^(15)+79296548109497025481500000000000000*(y)^(16)))/1.69182035233245e47",2); // 3
        iterLimit=5;
        break;
    default:
        GISMO_ERROR("No such case.");
        break;
    }
    std::cout << ss.str() << std::flush;

    size_t steps = iterLimit+1-startIter;
    gsMatrix<real_t> norms(3,steps);
    gsVector<real_t> times(steps),conditionNrs(steps);
    gsVector<size_t> dofs(steps);

    for(int i = startIter;i<=iterLimit;++i)
    {
        std::string path=pathStart+util::to_string(i)+pathEnd;
        std::string mapPath=pathStart+util::to_string(i)+mapPathEnd;
        gsMultiPatch<> * mp = (  (gsMultiPatch<>::uPtr)gsReadFile<>(path)  ).release();
        if(mp==NULL)
        {
            std::cout << "File could not be read." << "\n" ;
            return false;
        }
        mp->computeTopology();
        gsMappedBasis<2,real_t> compBasis(*mp,mapPath);

        // Define Boundary conditions
        gsBoundaryConditions<real_t> bcInfo;
        gsFunctionExpr<real_t> g_copy=g;
        gsFunctionExpr<real_t> dg_copy=dg;
        gsFunctionExpr<real_t> ddg_copy=ddg;
        gsFunctionExpr<real_t> f_copy=f;
        addAllDirichletBoundaries(*mp,g_copy,bcInfo);

        real_t h,l2error,h1error,h2error,time,conditionNr;
        size_t dof;

        std::stringstream ss2;
        ss2<<"biharmonic"<<switch_var;
        solveBiharmonicCollectResults_oneStep(&compBasis,bcInfo,mp,g_copy,dg_copy,ddg_copy,f_copy,i,
                                                             h,l2error,h1error,h2error,dof,conditionNr,time,plotNr,ss2.str(),plot);

        norms(0,i-startIter)=l2error;
        norms(1,i-startIter)=h1error;
        norms(2,i-startIter)=h2error;
        times(i-startIter)=time;
        dofs(i-startIter)=dof;
        conditionNrs(i-startIter)=conditionNr;

        //gsInfo << "test passed\n\n" ;
        gsInfo << flush;
        delete mp;
    }

    std::vector<std::string> normNames;
    normNames.push_back("l2");
    normNames.push_back("h1");
    normNames.push_back("h2");
    printSummaryOfNormsAndTimes(normNames,norms,times,dofs,conditionNrs);

    plotNr++;
    return true;
}

bool test_fittingExampleG1Volume(int switch_var, int & plotNr)
{
    std::cout << "=======================================================\n";
    std::cout << "Fitting on a ";
    //int switch_var = 3;
    gsFunctionExpr<real_t> f; // function to fit
    std::string pathStart,pathEnd,mapPathEnd;
    int iterLimit=-1;
    switch(switch_var)
    {
    case 1:
        std::cout << "2-patch domain - volumes, degree 3\n";
        pathStart = "Katharina/2Patches/multipatch";
        pathEnd = ".xml";
        mapPathEnd = "_map.xml";
        f=gsFunctionExpr<real_t>("2*cos(2*x)*sin(2*y)*cos(2*z)",3);
        iterLimit=0;
        break;
    default:
        GISMO_ERROR("No such case.");
        break;
    }

    gsFunctionExpr<real_t> f_copy=f;
    gsMultiPatch<real_t> emptyMP;
    gsBoundaryConditions<real_t> emptyBC;
    gsCompositeBasisMapFromFileRefiner *refiner =new gsCompositeBasisMapFromFileRefiner(pathStart,pathEnd,pathStart,mapPathEnd,emptyMP,emptyBC,f_copy);
    gsMultiPatch<real_t> solution1 = refiner->getDomain();
    std::cout << "dimension: " << solution1.dim() << std::endl;
    gsWriteParaview(solution1,"KatharinaSolution0",27);
    gsRecipeAssembler *assembler =new gsRecipeAssemblerFitting(refiner->getDomain(),f_copy,refiner->getSpaces());

    std::vector<gsFunction<real_t>*> functions;
    functions.push_back(&f_copy);
    gsErrorEstimatorPerCellExact    estimator(refiner->getDomain(),functions);
    gsMatrix<real_t> normMat(1,2);
    normMat(0,0)=0;
    normMat(0,1)=2;
    estimator.setAllSeminorms(normMat);
    gsMatrix<real_t> norm2OfExactSol = calculateNorm2OfFunc(&f_copy,refiner->getDomain(),normMat);

    gsStopCriteriaIterationOnly     stopCriteria(iterLimit+1);
    gsEigenCGIdentity<real_t>       eigCG;
    eigCG.setMaxIterations(100000000);
    eigCG.setTolerance(1e-10);
    gsEigenCGIdentity<real_t>       eigCGEli;

    gsAdaptiveSolverLogData solver(*assembler,eigCG,&eigCGEli,estimator,NULL,*refiner,stopCriteria,norm2OfExactSol,"VolumeFitting");
    solver.adaptiveSolve();
    if(!solver.succeed())
    {
        gsWarn << "Adaptive Solver failed. Exiting...\n";
        return false;
    }
    gsMatrix<real_t> norms = solver.postprocessNorms();
    gsVector<real_t> times = solver.getTimes();
    gsVector<size_t> dofs = solver.getDofs();
    gsVector<real_t> conds = solver.getConditions();
    std::vector<std::string> normNames;
    normNames.push_back("l2");
    printSummaryOfNormsAndTimes(normNames,norms,times,dofs,conds);
    gsMultiPatch<real_t> solution = refiner->getDomain();
    gsWriteParaview(solution,"KatharinaSolution1",2000);


    delete refiner;
    delete assembler;

    plotNr++;
    return true;
}

//========================================== MAIN ========================================//

int main(int argc, char *argv[])
{
    bool passed = true;

    bool plot = false;
    bool g1geom = false;
    bool allTests = false;

    gsCmdLine cmd("Composite basis tests.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("g1geom", "Tests the g1-geometries", g1geom);
    cmd.addSwitch("testAll", "Run all tests other than g1-geometries", allTests);
	try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    int plotNr = 0;
    //plot=true;


    passed = test_poissonSolvingExample(5, plotNr) && passed;

    if(allTests)
    {
        unsigned poissonTest_limit=4;
        for(unsigned i = 1;i<=poissonTest_limit;++i)
            passed = test_poissonSolvingExample(i, plotNr) && passed;
    }

    if(g1geom)
    {
        unsigned fitting_limit=19;
        for(unsigned i = 1;i<=fitting_limit;++i)
            passed = test_fittingExampleG1(i, plotNr) && passed;

        unsigned poisson_limit=22;
        for(unsigned i = 21;i<=poisson_limit;++i)
            passed = test_poissonSolvingExampleG1(i, plotNr) && passed;

        unsigned biharmonic_limit=13;
        for(unsigned i = 8;i<=biharmonic_limit;++i)
            passed = test_biharmonicExampleG1(i,plot,plotNr) && passed;
    }

    //test_fittingExampleG1(10, plotNr);
    //test_poissonSolvingExampleG1(14, plotNr);

    //test_fittingExampleG1Volume(1, plotNr);

    return passed ? 0 : 1;
}
