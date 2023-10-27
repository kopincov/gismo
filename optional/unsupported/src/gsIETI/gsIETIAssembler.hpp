/**  gsIETIAssembler.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on:  2014-12-03
*/

#include <gsIETI/gsIETIAssembler.h>

#include <gsCore/gsField.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsAssembler/gsQuadrature.h>

#include <gsMultiGrid/gsMultiGridAdapter.h>
#include <gsSolver/gsTwoLevel.h>
#include <gsSolver/gsKronecker.h>
#include <gsSolver/gsCompositePrecOp.h>

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsAssembler/gsGenericAssembler.h>

#include <gsSolver/gsPowerIteration.h>

#include <gsCore/gsGeometryEvaluator.h>

namespace gismo {


template<class T>
gsIETIAssembler<T>::gsIETIAssembler(gsAssembler<T> &assembler)

    :m_assembler(&assembler), m_patches(assembler.patches()), m_bConditions(assembler.pde().boundaryConditions()),m_timings(m_patches.nPatches())
{
#ifdef _OPENMP
    Eigen::initParallel();
#endif
    m_useGivenMats = false;
    m_isInit=false;
    setOptions(defaultOptions());

    size_t nBasis = m_assembler->numMultiBasis();
    //fixme: needs to be reworked if using several bases.
    for(size_t c=0;c< m_assembler->system().dofMappers().size();c++)
    {
        if(c<nBasis)
            m_basis.push_back(m_assembler->multiBasis(c));
        else
            m_basis.push_back(m_assembler->multiBasis(nBasis-1));
    }

    //Partition the Boundary conditions
    m_bc_loc.resize(m_patches.nPatches());
    for(size_t np=0; np<m_patches.nPatches();np++)
        m_bConditions.getConditionsForPatch(np,m_bc_loc[np]);
}

template<class T>
gsOptionList gsIETIAssembler<T>::defaultOptions()
{
    gsOptionList opt;

    opt.addSwitch("ExtraTiming", "enable extraTimings", false);
    opt.addSwitch("NoMinimumEnergy","choose energy minimizing prim subspaces", false);
    opt.addSwitch("EnforceRescaling","enforces rescaling of local matrices",  false);
    opt.addSwitch("NonlinearMode","enables the usage for nonlinear PDEs",  false);
    opt.addSwitch("CalcRhsNorm","calculates the norm of the original rhs",  false);
    opt.addSwitch("SaddlePoint","use saddle point for inexact", false);
    opt.addSwitch("NoAssemblyOfRHS","IETI does not assemble the rhs, used for preconditioner", false);
    opt.addSwitch("NonSymmetric","Use IETI for non symmetric matrices", false);

    opt.addString("Strategy", "Choosen strategy","C");
    opt.addString("SolverKii", "chosen solver for Kii","D");
    opt.addString("SolverKC", "chosen solver for KC","D");
    opt.addString("SolverKrr", "chosen solver for Krr","D");
    opt.addString("Scaling", "chosen scaling","coeff");

    opt.addInt("nRhs", "specify the number of rhs to be solved", 1);


    opt.addInt   ("MG_KC.Cycles",             "Number of multi-grid cycles", 1);
    opt.addInt   ("MG_KC.Presmooth",          "Number of pre-smoothing steps", 1);
    opt.addInt   ("MG_KC.Postsmooth",         "Number of post-smoothing steps", 1);
    opt.addReal  ("MG_KC.Damping",            "Damping factor for the smoother (handed over to smoother)", 1.);
    opt.addReal  ("MG_KC.OuterDamping",       "Damping factor for the smoother (globally)", 1.);
    opt.addString("MG_KC.Smoother",            "Smoothing method",  "GaussSeidel");
    opt.addSwitch("MG_KC.UseP1Projection",    "enable the P1 projection of coarser grids, then uses TwoLevel", false);
    opt.addSwitch("MG_KC.BasedOnRankOneCorrection", "Setup the mass smoother (if chosen) based on a rank one approximation of the matrices", false );
    opt.addInt   ("MG_KC.TL.Presmooth",          "Number of pre-smoothing steps", 1);
    opt.addInt   ("MG_KC.TL.Postsmooth",         "Number of post-smoothing steps", 1);
    opt.addReal  ("MG_KC.TL.CoarseDamping",      "Damping factor for the two level smoother (handed over to smoother)", 1.);
    opt.addReal  ("MG_KC.TL.Damping",        "Damping factor for the smoother (handed over to smoother)", 1.);
    opt.addReal  ("MG_KC.TL.OuterDamping",       "Damping factor for the smoother (globally)", 1.);
    opt.addString("MG_KC.TL.Smoother",       "Smoothing method", "SubspaceCorrectedMassSmoother");
    opt.addSwitch("MG_KC.TL.BasedOnRankOneCorrection", "Setup the mass smoother (if chosen) based on a rank one approximation of the matrices", false );

    opt.addInt   ("MG_Kii.Cycles",             "Number of multi-grid cycles", 1);
    opt.addInt   ("MG_Kii.Presmooth",          "Number of pre-smoothing steps", 1);
    opt.addInt   ("MG_Kii.Postsmooth",         "Number of post-smoothing steps", 1);
    opt.addReal  ("MG_Kii.Damping",            "Damping factor for the smoother (handed over to smoother)", 1.);
    opt.addReal  ("MG_Kii.OuterDamping",       "Damping factor for the smoother (globally)", 1.);
    opt.addString("MG_Kii.Smoother",             "Smoothing method",  "GaussSeidel");
    opt.addSwitch("MG_Kii.UseP1Projection","enable the P1 projection of coarser grids, then uses TwoLevel", false);
    opt.addSwitch("MG_Kii.BasedOnRankOneCorrection", "Setup the mass smoother (if chosen) based on a rank one approximation of the matrices", false );
    opt.addInt   ("MG_Kii.TL.Presmooth",          "Number of pre-smoothing steps", 1);
    opt.addInt   ("MG_Kii.TL.Postsmooth",         "Number of post-smoothing steps", 1);
    opt.addReal  ("MG_Kii.TL.CoarseDamping",      "Damping factor for the two level smoother (handed over to smoother)", 1.);
    opt.addReal  ("MG_Kii.TL.Damping",         "Damping factor for the smoother (handed over to smoother)", 1.);
    opt.addReal  ("MG_Kii.TL.OuterDamping",       "Damping factor for the smoother (globally)", 1.);
    opt.addString("MG_Kii.TL.Smoother",       "Smoothing method", "SubspaceCorrectedMassSmoother");
    opt.addSwitch("MG_Kii.TL.BasedOnRankOneCorrection", "Setup the mass smoother (if chosen) based on a rank one approximation of the matrices", false );

    opt.addReal ("Regularization", " regularization parameter for local neumann problems (inexact problems) ", 1.e-2);

    opt.addReal ("MG_Kii.RHS.Tolerance",            "Stopping criterion for cg", 1.e-10);
    opt.addInt  ("MG_Kii.RHS.MaxIterations",        "Stopping criterion for cg", 200);
    opt.addReal ("MG_Kii.Prec.Tolerance",           "Stopping criterion for cg", 1.e-2);
    opt.addInt  ("MG_Kii.Prec.MaxIterations",       "Stopping criterion for cg", 2);
    opt.addReal ("MG_KC.Solve.Tolerance",           "Stopping criterion for cg", 1.e-10);
    opt.addInt  ("MG_KC.Solve.MaxIterations",       "Stopping criterion for cg", 200);
    opt.addReal ("MG_KC.Basis.Tolerance",           "Stopping criterion for cg", 1.e-12);
    opt.addInt  ("MG_KC.Basis.MaxIterations",       "Stopping criterion for cg", 200);
    opt.addReal ("MG_KC.Basis.PreconditionScaling1","K-Scaling for SchöberlZuLehner", 1.24);
    opt.addReal ("MG_KC.Basis.PreconditionScaling2","S-Scaling for SchöberlZuLehner", 0.99/opt.getReal("MG_KC.Basis.PreconditionScaling1"));


    //if saddlePoint == true:
    //MG_KC.Solve.Tolerance = 1.e-2
    //MG_KC.Solve.MaxIterations = 5

    return opt;
}

template<class T>
void gsIETIAssembler<T>::setOptions(const gsOptionList &opt)
{
    if(m_isInit)
        printWarn("You have already initialized IETI, you should be carefull, when you edit the options!");

    m_IETIoptions.opt.update(opt,gsOptionList::addIfUnknown); // A copy is made here

    //   m_IETIoptions.opt.update(m_assembler->options(),gsOptionList::addIfUnknown); //extract the options from the assembler also in the IETI namespace
    m_IETIoptions.KiiSolver = IETILocalSolver::chooseSolver(m_IETIoptions.opt.getString("SolverKii"));
    m_IETIoptions.KrrSolver = IETILocalSolver::chooseSolver(m_IETIoptions.opt.getString("SolverKrr"));
    m_IETIoptions.KCSolver =  IETILocalSolver::chooseSolver(m_IETIoptions.opt.getString("SolverKC"));
    m_IETIoptions.strat.strat = primalDofMethod::chooseMethod(m_IETIoptions.opt.getString("Strategy"));
    m_IETIoptions.scal = IETIPrecondScaling::chooseMethod(m_IETIoptions.opt.getString("Scaling"));

    info.nRhs = m_IETIoptions.opt.getInt("nRhs");

#if defined(GISMO_WITH_PARDISO)
 //   print("Using Pardiso!");
    //m_options.needRescaling = false;
#elif defined(_OPENMP)
    //  print("Using SparseLU!");
#elif defined (GISMO_WITH_SUPERLU)
    //  print("Using SuperLU!");
    //m_options.needRescaling = false;
#endif

}

template<class T>
void gsIETIAssembler<T>::init()
{

    //check the primal dof method
    checkPrimalDofStrategy(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy"),m_patches.parDim());

    //initialize the patchInterfaceSides and the local dimensions of the basis
    init_patchInterfaceSides();

    //Initialize the mappers
    init_mappers();

    //initialize data for jump operator
    init_jumpOperatorData();

    //initialize the boundary Dofs and all the bookkeeping included (interior vs. boundary/interface Dofs)
    init_boundaryDofs();

    //initialize the primal Dofs and all the bookkeeping included
    init_primalDofs();

    //populate the info structure
    setInfo();

    //initialize the remaining Dofs() and all the bookkeeping included
    if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
        init_remainingDofs(); //needs to be done after setInfo.

    init_InitialGuess();

    //since #267 the columns of the rhs of m_system is required
    // TODO:: Is this the right place for this??
    m_assembler->system().rhs().conservativeResize(Eigen::NoChange, info.nRhs);
    
    m_LU_Spp = NULL;
    m_isInit=true;
}

template<class T>
void gsIETIAssembler<T>::init_mappers()
{
    size_t cDim = m_basis.size();

    //Create different mappers;
    m_locDofsMapper.resize(m_patches.nPatches());
    for(size_t np = 0; np <m_patches.nPatches();np++)
        m_locDofsMapper[np].reserve(cDim);
    m_primalDofMapper.resize(cDim);
    m_stdMapper.resize(cDim);

    if ( m_assembler->options().getInt("DirichletStrategy") == dirichlet::elimination)
    {
        for(size_t c =0; c<cDim;c++)
        {
            //Create localDofMappers (where the boundary is eliminated)
            for(size_t np = 0; np<m_patches.nPatches();np++)
            {
                gsMultiBasis<T> mBs(m_basis[c].basis(np));//transform to MultiBasis

                gsDofMapper locdofMapper(mBs, m_bc_loc[np]);
                locdofMapper.finalize();
                m_locDofsMapper[np].push_back(locdofMapper);
            }
            m_basis[c].getMapper(true,m_bConditions, m_stdMapper[c]);
          //  m_stdMapper[c] = m_assembler->system().colMapper(c);

            //The primalDofMapper still needs to be finalized -> init_PrimalDofs
            m_primalDofMapper[c].init(m_basis[c], m_bConditions);
        }

    }
    else
    {
        for(size_t c =0; c<cDim;c++)
        {
            //Create localDofMappers (where the boundary is eliminated)
            for(size_t np = 0; np<m_patches.nPatches();np++)
            {
                gsMultiBasis<T> mBs(m_basis[c].basis(np));//transform to MultiBasis
                gsDofMapper locdofMapper(mBs);
                locdofMapper.finalize();
                m_locDofsMapper[np].push_back(locdofMapper);
            }
            m_basis[c].getMapper(true, m_stdMapper[c]);
           // m_stdMapper[c] = m_assembler->system().colMapper(c);

            //The primalDofMapper still needs to be finalized -> init_PrimalDofs
            m_primalDofMapper[c].init(m_basis[c]);
        }

    }
    unsigned offset = 0, offset2 = 0, offsetb = 0,offsetb2=0;
    for(size_t np = 0;np<m_patches.nPatches();np++)
    {
        offset = offsetb = 0;
        for(size_t c =0; c<cDim;c++)
        {
            m_locDofsMapper[np][c].setShift(offset);
            offset+=m_locDofsMapper[np][c].freeSize();

            m_locDofsMapper[np][c].setBoundaryShift(offsetb);
            offsetb+=m_locDofsMapper[np][c].boundarySize();
        }
    }
    for(size_t c =0; c<cDim;c++)
    {
        m_stdMapper[c].setShift(offset2);
        offset2+=m_stdMapper[c].freeSize();

        m_stdMapper[c].setBoundaryShift(offsetb2);
        offsetb2+=m_stdMapper[c].boundarySize();
    }


    if(m_patches.parDim() ==3 && m_IETIoptions.strat.contains(primalDofMethod::edges))
    {
        gsVector<index_t> entr(m_patches.nPatches());
        for(size_t np=0;np<m_patches.nPatches();np++)
            entr[np]=12;
        //Each patch in 3D consists of 12 edges
        m_edgeMapper = gsDofMapper(entr);

        // Set the entries of Edge
        InstantiateEdge();
    }
}


template<class T>
void gsIETIAssembler<T>::setInfo()
{
    info.numberPatches = m_patches.nPatches();
    info.dim = m_patches.parDim(); //Maybe dim()
    info.cDim = m_basis.size();

    //this is now directly set in setOptions, use the option "nRhs".
    //info.nRhs = m_assembler->system().rhsCols(); // before const 1

    info.dofTotal=0;
    info.origSystemSize =0;
    for(size_t c=0;c<info.cDim;c++)
    {
        info.dofTotal += m_primalDofMapper[c].freeSize();
        info.origSystemSize +=m_stdMapper[c].freeSize();
    }
    info.dofTotalB=0;
    info.dofTotalI=0;

    info.dofsB.resize(m_patches.nPatches(),0);
    info.dofsI.resize(m_patches.nPatches(),0);
    info.dofsP.resize(m_patches.nPatches(),0);
    info.dofs.resize(m_patches.nPatches(),0);


    int diffMethods=3;
    info.dofTotalPtype.resize(diffMethods,0);
    info.dofsPtype.resize(m_patches.nPatches());

    //set size on patch np for dofs, boundaryDofs, interiorDofs
    for (size_t np=0; np < m_patches.nPatches(); ++np )
    {
        info.dofsPtype[np].resize(diffMethods);

        for(size_t c=0;c<info.cDim;c++)
            info.dofs[np]+=(m_locDofsMapper[np][c].freeSize());

        info.dofsB[np]=(m_boundDofs[np].size());
        info.dofsI[np]=(info.dofs[np]-info.dofsB[np]);
        info.dofTotalB+= info.dofsB[np];
        info.dofTotalI+= info.dofsI[np];
    }


    //set size for primal dofs
    info.dofTotalP=0;
    int method;
    if(m_IETIoptions.strat.contains( primalDofMethod::vertices))
    {
        method = primalDofMethod::getMethodIndex(primalDofMethod::vertices);

        for (size_t np=0; np < m_patches.nPatches(); ++np )
            info.dofsPtype[np][method] = m_primalVdofs[np].size();
        for(size_t c=0;c<info.cDim;c++)
            info.dofTotalPtype[method]+=m_primalDofMapper[c].coupledSize();

        info.dofTotalP +=  info.dofTotalPtype[method];
    }
    if((info.dim == 2 && m_IETIoptions.strat.contains( primalDofMethod::edges ) ) ||  (info.dim ==3 && m_IETIoptions.strat.contains(primalDofMethod::faces) ))
    {
        if(info.dim == 2)
            method =primalDofMethod::getMethodIndex(primalDofMethod::edges);
        else
            method =primalDofMethod::getMethodIndex(primalDofMethod::faces);

        info.dofTotalPtype[method]=0;
        for (size_t np=0; np < m_patches.nPatches(); ++np )
        {
            info.dofsPtype[np][method] = info.cDim*m_patchAveragesSides[np].size();
            info.dofTotalPtype[method]+=info.dofsPtype[np][method] ;
        }
        info.dofTotalPtype[method]-= info.cDim*m_patches.nInterfaces();
        info.dofTotalP +=   info.dofTotalPtype[method];
    }
    if(info.dim ==3 && ( m_IETIoptions.strat.contains(primalDofMethod::edges)))
    {
        method = primalDofMethod::getMethodIndex(primalDofMethod::edges);
        for (size_t np=0; np < m_patches.nPatches(); ++np )
            info.dofsPtype[np][method] = info.cDim*m_primalEdges[np].size();
        info.dofTotalPtype[method]= info.cDim*m_edgeMapper.coupledSize();
        info.dofTotalP +=   info.dofTotalPtype[method];
    }


    //set size for the remaining dofs and the sum of primal dofs on patch np
    info.dofsR.resize(m_patches.nPatches());
    for (size_t np=0; np < m_patches.nPatches(); ++np )
    {
        unsigned sum =0;
        for(int i = 0; i<diffMethods;i++)
            sum+= info.dofsPtype[np][i];
        info.dofsP[np]=sum;

        //remaining are all dofs minus primal vertex dofs (only for noMinimalenergy)
        info.dofsR[np]=info.dofs[np] - info.dofsPtype[np][primalDofMethod::getMethodIndex(primalDofMethod::vertices)];
    }

    //for minimal energy, we only consider boundary dofs
    if(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
        info.dofs = info.dofsB;

    //number of lagrange multiplier
    info.lagrangeMult = m_lagrangeTable.size();

    //number of global Interface dofs
    info.dofInterface= m_globalConnectionTable.size();
}

template<class T>
void gsIETIAssembler<T>::init_jumpOperatorData()
{
    size_t cDim = m_basis.size();
    std::set<index_t> freeElimDof;
    for ( typename gsMultiPatch<T>::const_iiterator it =
          m_patches.iBegin(); it != m_patches.iEnd(); ++it )
    {
        boundaryInterface iFace=(*it);

        patchSide p1 = iFace.first();
        patchSide p2 = iFace.second();

        gsMatrix<index_t> b1, b2;

        for(size_t c=0;c<cDim;c++)
        {
            m_basis[c].basis(p1.patch).matchWith(iFace,
                                                 m_basis[c].basis(p2.patch), b1, b2);

            gsMatrix<index_t> glob1;
            m_stdMapper[c].localToGlobal(b1,p1.patch,glob1);


            for(int i=0;i<glob1.rows();i++)
            {
                GISMO_ASSERT(m_stdMapper[c].index(b1(i,0),p1.patch)==m_stdMapper[c].index(b2(i,0),p2.patch), "Corner points are not coupled correctly!");
                bool isfree1 = m_locDofsMapper[p1.patch][c].is_free(b1(i,0));
                bool isfree2 = m_locDofsMapper[p2.patch][c].is_free(b2(i,0));

                //if both are free just add them
                if(isfree1 && isfree2)
                {
                    m_globalConnectionTable[glob1(i,0)].insert(patchDof(p1.patch,compCalc(p1.patch,b1(i,0),c)) );
                    m_globalConnectionTable[glob1(i,0)].insert(patchDof(p2.patch,compCalc(p2.patch,b2(i,0),c)) );
                }
                else if(isfree1 || isfree2)
                {
                    //if one of the two are not free check if this global dof already has a dirichlet dof
                    if(freeElimDof.find(glob1(i,0))==freeElimDof.end())
                    {
                        //if not then add  it to the free-eliminated dofs and add both dof to the map
                        freeElimDof.insert(glob1(i,0));
                        m_globalConnectionTable[glob1(i,0)].insert(patchDof(p1.patch,compCalc(p1.patch,b1(i,0),c)) );
                        m_globalConnectionTable[glob1(i,0)].insert(patchDof(p2.patch,compCalc(p2.patch,b2(i,0),c)) );
                    }
                    else
                    {
                        //if a dirichlet dof was already there, then only add the free dof.
                        if(isfree1)
                            m_globalConnectionTable[glob1(i,0)].insert(patchDof(p1.patch,compCalc(p1.patch,b1(i,0),c)) );
                        if(isfree2)
                            m_globalConnectionTable[glob1(i,0)].insert(patchDof(p2.patch,compCalc(p2.patch,b2(i,0),c)) );
                    }
                }

            }//for i

        }
    }

    m_lagrangeTable.reserve(m_globalConnectionTable.size());

    bool flag=false;

    bool isfreeA, isfreeB;
    for(std::map<index_t,std::set<patchDof> >::iterator it=m_globalConnectionTable.begin(); it!=m_globalConnectionTable.end(); ++it)
    {
        flag = false;
        if(freeElimDof.find((*it).first)!=freeElimDof.end())
            flag=true;
        for(std::set<patchDof>::iterator itA=(*it).second.begin(); itA!=(*it).second.end(); ++itA)
        {
            if(itA==(*it).second.end())
                break;

            //This can be improved
            if(flag)
                isfreeA = m_locDofsMapper[(*itA).first][0].is_free((*itA).second);
            else
                isfreeA = true;

            for(std::set<patchDof>::iterator itB=itA; itB!=(*it).second.end(); ++itB)
            {
                if(itB==itA)
                    continue;

                //This can be improved
                if(flag)
                    isfreeB = m_locDofsMapper[(*itB).first][0].is_free((*itB).second);
                else
                    isfreeB = true;

                if(isfreeA && isfreeB)
                    m_lagrangeTable.push_back(std::pair<patchDof,patchDof>(*itA,*itB));
                else if(isfreeA || isfreeB)
                {
                    m_lagrangeTable.push_back(std::pair<patchDof,patchDof>(*itA,*itB));
                    m_freeElimLagr.push_back(m_lagrangeTable.size()-1);

                }
            }
        }
    }

}

template<class T>
void gsIETIAssembler<T>::init_patchInterfaceSides()
{
    m_patchISides.resize(m_patches.nPatches());

    std::vector<boundaryInterface> ifaces = m_patches.interfaces();
    for(std::vector<boundaryInterface>::iterator it = ifaces.begin(); it!=ifaces.end();it++)
    {
        patchSide side = (*it).first();
        patchSide side2 = (*it).second();

        m_patchISides[side.patch].push_back(side);
        m_patchISides[side2.patch].push_back(side2);

    }
    m_patchAveragesSides= m_patchISides;
    /*
    int d = m_patches.parDim();
    for(size_t np=0; np<m_patches.nPatches();++np)
    {
        const typename gsBoundaryConditions<T>::bcContainer& bcs = m_bc_loc[np].dirichletSides();
        for(boxSide side=boxSide::getFirst(d);side!=boxSide::getEnd(d);++side)
        {
            bool contains = false;
            for(typename gsBoundaryConditions<T>::bcContainer::const_iterator bc = bcs.begin();bc!=bcs.end();bc++)
            {
                if(bc->side()==side)
                {
                    contains = true;
                    break;
                }
            }
            if(!contains && std::find(m_patchISides[np].begin(), m_patchISides[np].end(), side) == m_patchISides[np].end() )
                m_patchAveragesSides[np].push_back(patchSide(np,side));
        }

    }
    //*/

}
template<class T>
void gsIETIAssembler<T>::init_primalDofs()
{

    int dim = m_patches.parDim();
    size_t cDim = m_basis.size();
    m_primalVdofs.resize(m_patches.nPatches());
    m_primalEdges.resize(m_patches.nPatches());

    std::vector<boundaryInterface> ifaces = m_patches.interfaces();
    std::vector<boxCorner> corners, corners2;
    gsMatrix<index_t> b1,b2;
    gsVector<int>  bSize(dim-1);
    for(int i=0; i<dim-1;i++)
        bSize[i]=2;

    for(std::vector<boundaryInterface>::iterator it = ifaces.begin(); it!=ifaces.end();it++)
    {
        boundaryInterface iFace = (*it);
        patchSide side = iFace.first();
        patchSide side2 = iFace.second();

        if(m_IETIoptions.strat.contains(primalDofMethod::vertices))
        {

            corners.clear();
            corners2.clear();
            side.getContainedCorners(dim,corners);
            side2.getContainedCorners(dim,corners2);

            b1.setZero(corners.size(),1);
            b2.setZero(corners2.size(),1);

            for(size_t c=0;c<cDim;c++)
            {
                for(size_t i=0;i<corners.size();i++)
                {
                    b1(i,0) =m_basis[c][side.patch].functionAtCorner(corners[i]);
                    b2(i,0) =m_basis[c][side2.patch].functionAtCorner(corners2[i]);
                }
                iFace.reorderCorners(b1);

                for(size_t i=0;i<corners.size();i++)
                {
                    GISMO_ASSERT(m_stdMapper[c].index(b1(i,0),side.patch)==m_stdMapper[c].index(b2(i,0),side2.patch), "Corner points are not coupled correctly!");

                    if(m_locDofsMapper[side.patch][c].is_free(b1(i,0) )&& m_locDofsMapper[side2.patch][c].is_free(b2(i,0) ) )
                    {
                        m_primalVdofs[side.patch].push_sorted_unique(compCalc(side.patch,b1(i,0),c));
                        m_primalVdofs[side2.patch].push_sorted_unique(compCalc(side2.patch,b2(i,0),c));

                        m_primalDofMapper[c].matchDof(side.patch,b1(i,0),side2.patch,b2(i,0));
                        //m_primalDofMapper[c].matchDof(side.patch,compCalc(side.patch,b1(i,0),c),side2.patch,compCalc(side2.patch,b2(i,0),c));
                    }
                }


            }

        }
        if( dim == 3&& m_IETIoptions.strat.contains(primalDofMethod::edges))
        {
            std::vector<Edge> e1, e2;
            e1.reserve(4);
            e2.reserve(4);

            getEdgesFromBoundary(side,e1);
            getEdgesFromBoundary(side2,e2);

            gsMatrix<index_t> num(4,1);
            num << 0,1,2,3;
            iFace.reorderCorners(num); // get permutation of corners


            for(size_t i = 0; i<e1.size();i++)
            {
                if(!e1[num(i,0)].isEliminated && !e2[i].isEliminated && e1[num(i,0)].number >2 && e2[i].number>2)
                {
                    //for(size_t c=0;c<cDim;c++)
                    //{
                    m_primalEdges[side.patch].push_sorted_unique(e1[num(i,0)]);
                    m_primalEdges[side2.patch].push_sorted_unique(e2[i]);
                    //}
                    m_edgeMapper.matchDof(side.patch,e1[num(i,0)].edgeNumber,side2.patch,e2[i].edgeNumber);
                }


            }

        }

    }

    unsigned off = 0;
    std::vector<index_t> pShifts;

    for(size_t c=0;c<cDim;c++)
    {
        pShifts.push_back(off);
        m_primalDofMapper[c].finalize();
        m_primalDofMapper[c].setShift(off);
        off+=m_primalDofMapper[c].freeSize();
    }

    m_edgeMapper.finalize();


    //populate the mappings for the Primal dofs (loc->glob)
    m_pDofsLoc2Glob.resize(m_patches.nPatches());

    size_t c = 0;
    if(m_IETIoptions.strat.contains(primalDofMethod::vertices))
    {
        for(size_t np=0; np<m_patches.nPatches();np++)
        {
            for(std::vector<index_t>::iterator it=m_primalVdofs[np].begin();  it!=m_primalVdofs[np].end();it++)
            {
                c=getComp(np,*it);

                // the number of standard dofs (not coupled and not eliminated
                int offset =(m_primalDofMapper[c].freeSize()-m_primalDofMapper[c].coupledSize());
                int idx= m_primalDofMapper[c].index(compCalcBack(np,*it),np)-pShifts[c]-offset; //so the primal dofs have get a index starting from zero.

                //add the other components.
                for(size_t cc = 0;cc<c;cc++)
                    idx+= m_primalDofMapper[cc].coupledSize();

                m_pDofsLoc2Glob[np].push_back(idx);
            }
        }


    }

    // m_pDofsLoc2Glob[np][offset : end] also corresponds to patchISides[np][1:end]
    int offset = 0;
    for(size_t cSizeT1=0;cSizeT1<cDim;cSizeT1++)
        offset+= m_primalDofMapper[cSizeT1].coupledSize();

    if((dim == 2 && m_IETIoptions.strat.contains(primalDofMethod::edges)) || (dim == 3 && m_IETIoptions.strat.contains(primalDofMethod::faces) ) )
    {
        for(std::vector<boundaryInterface>::iterator it = ifaces.begin(); it!=ifaces.end();it++)
        {
            patchSide side = (*it).first();
            patchSide side2 = (*it).second();

            for(size_t cSizeT2 = 0;cSizeT2<cDim;cSizeT2++)
            {
                m_pDofsLoc2Glob[side.patch].push_back(offset);
                m_pDofsLoc2Glob[side2.patch].push_back(offset);
                offset++;
            }
        }
        for(size_t np=0; np<m_patches.nPatches();np++)
        {
            for(size_t i=0; i<m_patchAveragesSides[np].size() - m_patchISides[np].size();i++)
                for(size_t cSizeT2 = 0;cSizeT2<cDim;cSizeT2++)
                {
                    m_pDofsLoc2Glob[np].push_back(offset);
                    offset++;
                }
        }

    }


    if(dim ==3 && m_IETIoptions.strat.contains(primalDofMethod::edges))
    {
        for(size_t np=0; np<m_patches.nPatches();np++)
        {

            //delete the duplicates in m_primalVdofs
            for(std::vector<Edge>::iterator it=m_primalEdges[np].begin();  it!=m_primalEdges[np].end();it++)
            {
                // if(*it>=m_primalEdges[np].size())
                //continue;

                //the number of non matched edges (to start from 0)
                int offset2 =(m_edgeMapper.freeSize()-m_edgeMapper.coupledSize());
                // the index starts from the number of Vertex and Face dofs (==offset).
                int idx= m_edgeMapper.index((*it).edgeNumber,np)-offset2+offset;

                for(size_t cSizeT3=0; cSizeT3<cDim;cSizeT3++)
                    m_pDofsLoc2Glob[np].push_back(idx+cSizeT3*m_edgeMapper.coupledSize());
            }
        }
    }

}

template<class T>
void gsIETIAssembler<T>::init_remainingDofs()
{
    size_t cDim = m_basis.size();
    std::vector<index_t>::iterator it;
    unsigned iR,sum,idx;

    m_glob2PrimalRemainingIndex.resize(m_patches.nPatches());
    m_globIsPrimalIndex.resize(m_patches.nPatches());

    for(size_t np=0; np<m_patches.nPatches();np++)
    {
        iR=sum=0;

        for(size_t c=0;c<cDim;c++)
            sum+=m_locDofsMapper[np][c].freeSize();

        m_globIsPrimalIndex[np].resize(sum);
        m_glob2PrimalRemainingIndex[np].resize(sum);

        for(size_t c=0; c<cDim;c++)
            for(int i = 0; i<m_locDofsMapper[np][c].size();i++)
            {
                idx = m_locDofsMapper[np][c].index(i);
                if(m_locDofsMapper[np][c].is_free_index(idx))
                {
                    it = std::find(m_primalVdofs[np].begin(), m_primalVdofs[np].end(), compCalc(np,i,c));
                    if(it!=m_primalVdofs[np].end())
                    {
                        m_globIsPrimalIndex[np][idx]=true;
                        m_glob2PrimalRemainingIndex[np][idx]=std::distance(m_primalVdofs[np].begin(), it);
                    }
                    else
                    {
                        m_globIsPrimalIndex[np][idx]=false;
                        m_glob2PrimalRemainingIndex[np][idx]=iR;
                        iR++;
                    }
                }
            }
    }


    m_remDofs.resize(info.numberPatches);
    m_remDofIsAlsoBoundDof.resize(info.numberPatches);

    for(size_t np=0; np<info.numberPatches;np++)
    {
        m_remDofs[np].reserve(info.dofsR[np]);
        m_remDofIsAlsoBoundDof[np].reserve(info.dofsR[np]);

        //loop over all dofs, included the dirichlet ones
        for(size_t c = 0; c< info.cDim;c++)
        {
            for(index_t d = 0; d<info.dofs[np]+m_locDofsMapper[np][c].boundarySize();d++)
            {
                int idx2 = m_primalDofMapper[c].index(d,np);
                if(!m_primalDofMapper[c].is_coupled_index(idx2) && !m_primalDofMapper[c].is_boundary_index(idx2))
                {
                    m_remDofs[np].push_back(compCalc(np,d,c));

                    //also test if d is a interface dof, i.e. is in m_boundDofs[np]:
                    if(std::find(m_boundDofs[np].begin(), m_boundDofs[np].end(), compCalc(np,d,c))!=m_boundDofs[np].end())
                        m_remDofIsAlsoBoundDof[np].push_back(true);
                    else
                        m_remDofIsAlsoBoundDof[np].push_back(false);

                }
            }
        }
    }
}

template<class T>
void gsIETIAssembler<T>::init_boundaryDofs()
{
    unsigned cDim = m_basis.size();

    m_boundDofs.resize(m_patches.nPatches());


    for ( typename gsMultiPatch<T>::iiterator it =
          m_patches.iBegin(); it != m_patches.iEnd(); ++it )
    {
        patchSide p1 = (*it).first();
        patchSide p2 = (*it).second();

        for(size_t c=0; c<cDim;c++)
        {
            gsMatrix<index_t> b1 = m_basis[c].basis(p1.patch).boundary(p1);
            gsMatrix<index_t> b2 = m_basis[c].basis(p2.patch).boundary(p2);

            for(int i=0;i<b1.rows();i++)
            {
                if(m_locDofsMapper[p1.patch][c].is_free((b1)(i,0)))
                    m_boundDofs[p1.patch].push_sorted_unique(compCalc(p1.patch,(b1)(i,0),c));
            }
            for(int i=0;i<b2.rows();i++)
            {
                if(m_locDofsMapper[p2.patch][c].is_free((b2)(i,0)))
                    m_boundDofs[p2.patch].push_sorted_unique(compCalc(p2.patch,(b2)(i,0),c));
            }
        }
    }
    //*/
    /*
    for(size_t np=0; np<m_patches.nPatches();++np)
    {
        for(size_t i=0; i<m_patchAveragesSides[np].size();++i)
        {
            patchSide p1 = m_patchAveragesSides[np][i];
            for(size_t c=0; c<cDim;c++)
            {
                gsMatrix<index_t> b1 = m_basis[c].basis(p1.patch).boundary(p1);
                for(int i=0;i<b1.rows();i++)
                {
                    if(m_locDofsMapper[p1.patch][c].is_free((b1)(i,0)))
                        m_boundDofs[p1.patch].push_sorted_unique(compCalc(p1.patch,(b1)(i,0),c));
                }
            }
        }
    }
//*/
    std::vector<index_t>::iterator it;
    unsigned iI,iR,sum,idx;

    m_glob2BoundInteriorIndex.resize(m_patches.nPatches());
    m_globIsBoundIndex.resize(m_patches.nPatches());

    for(size_t np=0; np<m_patches.nPatches();np++)
    {
        iI =iR=sum=0;

        for(size_t c=0;c<cDim;c++)
            sum+=m_locDofsMapper[np][c].freeSize();

        m_globIsBoundIndex[np].resize(sum);
        m_glob2BoundInteriorIndex[np].resize(sum);

        for(size_t c=0; c<cDim;c++)
        {
            for(int i = 0; i<m_locDofsMapper[np][c].size();i++)
            {
                idx = m_locDofsMapper[np][c].index(i);
                if(m_locDofsMapper[np][c].is_free_index(idx))
                {
                    it = std::find(m_boundDofs[np].begin(), m_boundDofs[np].end(), compCalc(np,i,c));

                    if(it!=m_boundDofs[np].end())
                    {
                        m_globIsBoundIndex[np][idx]=true;
                        m_glob2BoundInteriorIndex[np].at(idx)=std::distance(m_boundDofs[np].begin(), it);
                    }
                    else
                    {
                        m_globIsBoundIndex[np][idx]=false;
                        m_glob2BoundInteriorIndex[np].at(idx)=iI;
                        iI++;
                    }
                }
            }
        }
    }
}

template<class T>
void gsIETIAssembler<T>::init_InitialGuess()
{
    if(m_IETIoptions.KCSolver == IETILocalSolver::Direct)
        return;

    gsVector<index_t> sizes( m_patches.parDim());
    gsMatrix<index_t> slice_l, slice_r ;
    m_primalPerm.resize(info.numberPatches);
    m_primalCorners.resize(info.numberPatches);
    gsVector<index_t> perm;
    for(size_t np = 0; np<info.numberPatches;++np)
    {
        m_primalPerm[np].resize(m_patches.parDim());

        switch(m_patches.parDim()){
        case 1:
        {
            const gsTensorBasis<1,T>& TB=  static_cast<const gsTensorBasis<1,T>& >( m_basis[0][np]);
            sizes[0]= cast<int,index_t>(TB.size(0));
            break;
        }
        case 2:
        {
            const gsTensorBasis<2,T>& TB=  static_cast<const gsTensorBasis<2,T>& >( m_basis[0][np]);
            TB.size_cwise(sizes);
            break;
        }
        case 3:
        {
            const gsTensorBasis<3,T>& TB=  static_cast<const gsTensorBasis<3,T>& >( m_basis[0][np]);
            TB.size_cwise(sizes);
            break;
        }
        };
        gsVector<bool> good = gsVector<bool>::Constant(m_patches.parDim(),true);
        for(typename gsBoundaryConditions<T>::const_iterator dIt=m_bc_loc[np].dirichletBegin(); dIt!=m_bc_loc[np].dirichletEnd(); ++dIt)
            good(dIt->side().direction())=false;

        for(index_t d=0; d<m_patches.parDim();++d)
        {
            if(good[d])
                perm.setLinSpaced(info.dofsB[np]+info.dofsI[np]+info.dofsP[np],0,info.dofsB[np]+info.dofsI[np]+info.dofsP[np]);
            else
            {
                perm.setConstant(1,0);
                m_primalPerm[np][d] = Eigen::PermutationMatrix<Dynamic,Dynamic,index_t>(perm);
                continue;
            }

            for(int i=0; i<(int)sizes(d)/2;++i)
            {

                switch(m_patches.parDim()){
                case 1:
                {
                    slice_l.resize(1,1);
                    (slice_l)(1,1) = i;
                    slice_r.resize(1,1);
                    (slice_r)(1,1) = (sizes(d)-i-1);
                    break;
                }
                case 2:
                {
                    const gsTensorBasis<2,T>& TB=  static_cast<const gsTensorBasis<2,T>& >( m_basis[0][np]);
                    slice_l = TB.coefSlice(d,i);
                    slice_r = TB.coefSlice(d,sizes(d)-i-1);
                    break;
                }
                case 3:
                {
                    const gsTensorBasis<3,T>& TB=  static_cast<const gsTensorBasis<3,T>& >( m_basis[0][np]);
                    slice_l = TB.coefSlice(d,i);
                    slice_r = TB.coefSlice(d,sizes(d)-i-1);
                    break;
                }
                };

                for(index_t j=0; j<slice_l.size();++j)
                {
                    int a = m_locDofsMapper[np][0].index((slice_l)(j));
                    int b = m_locDofsMapper[np][0].index((slice_r)(j));

                    bool af = m_locDofsMapper[np][0].is_free_index(a);
                    bool bf = m_locDofsMapper[np][0].is_free_index(b);

                    if(af && bf)
                    {
                        perm( a ) = cast<int,index_t>(b);
                        perm( b) = cast<int,index_t>(a);
                    }
                    else if(!af && !bf)
                        continue;
                    else if(!af || !bf)
                        GISMO_ERROR("Something went wrong with permutation init for primal initial guess");

                }
                //TODO: switch the primal variables! -> do this in application
                m_primalPerm[np][d] = Eigen::PermutationMatrix<Dynamic,Dynamic,index_t>(perm);

            }

        }

        if(m_IETIoptions.strat.contains(primalDofMethod::vertices))
        {
            m_primalCorners[np].reserve(m_primalVdofs[np].size());
            gsSortedVector<index_t>::iterator it;
            for(index_t i=boxCorner::getFirst(m_patches.parDim()); i<boxCorner::getEnd(m_patches.parDim());++i)
            {
                it = m_primalVdofs[np].find_it_or_fail(m_basis[0][np].functionAtCorner(boxCorner(i)));
                if( it!=m_primalVdofs[np].end())
                    m_primalCorners[np].push_back(boxCorner(i));

            }
        }

    }

}

template<class T>
void gsIETIAssembler<T>::assembleInit()
{
    GISMO_ASSERT(m_isInit, "IETI method not initialized, call .init() before .assemble()");

    m_scalingsKC.resize(info.numberPatches);
    m_scalingsKii.resize(info.numberPatches);
    m_scalingsKrr.resize(info.numberPatches);


    m_Kbb.resize(info.numberPatches);
    m_Kib.resize(info.numberPatches);
    m_Kbi.resize(info.numberPatches);
    m_LU_Kii.resize(info.numberPatches);
    m_rhs_d.resize(info.numberPatches);
    m_rhs_p.setZero(info.dofTotalP,info.nRhs);

    m_dirDofs.resize(info.numberPatches);

    m_permMat.resize(info.numberPatches);
    m_permInvMat.resize(info.numberPatches);
    if(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
    {
        m_LU_KC.resize(info.numberPatches);
        m_C.resize(info.numberPatches);

        m_rhs_i.resize(info.numberPatches);
        m_Phi.resize(info.numberPatches);
        m_PhiA.resize(info.numberPatches);

    }
    else
    {
        m_Krp.resize(info.numberPatches);
        m_Kpp.resize(info.numberPatches);
        m_LU_Krr.resize(info.numberPatches);
    }

    for(size_t c=0; c<info.cDim;++c)
        if(!m_useGivenMats) m_assembler->computeDirichletDofs(c);

    if(m_IETIoptions.opt.getSwitch("NonlinearMode"))
        m_assembler->homogenizeFixedDofs(-1);

    m_maxEigKii.resize(info.numberPatches, -1);
    switch (m_IETIoptions.KiiSolver) {
    case IETILocalSolver::Multigrid:
        m_inexact_Kii_MG.resize(info.numberPatches);
        m_Kii.resize(info.numberPatches);
        m_KiiSolve.resize(info.numberPatches);
        break;
    case IETILocalSolver::FastDiagonalization:
        m_inexact_Kii_BS.resize(info.numberPatches);
        m_Kii.resize(info.numberPatches);
        m_KiiSolve.resize(info.numberPatches);
        break;
    default:
        break;
    }


    m_maxEigKC.resize(info.numberPatches, -1);
    switch (m_IETIoptions.KCSolver) {
    case IETILocalSolver::Multigrid:
        m_inexact_KC.resize(info.numberPatches);
        m_Kii.resize(info.numberPatches);
        m_K.resize(info.numberPatches);
        m_S.resize(info.numberPatches);
        m_KC.resize(info.numberPatches);
        m_exC.resize(info.numberPatches);
        //m_Kreg.resize(info.numberPatches);
        m_KCBasis.resize(info.numberPatches);
        m_KCSolve.resize(info.numberPatches);
        // m_KM.resize(info.numberPatches);
        //  m_projK.resize(info.numberPatches);
        break;
    case IETILocalSolver::FastDiagonalization:
        m_inexact_KC.resize(info.numberPatches);
        m_KC.resize(info.numberPatches);
        m_S.resize(info.numberPatches);
        m_Kii.resize(info.numberPatches);
        m_K.resize(info.numberPatches);
        m_exC.resize(info.numberPatches);
        m_KCBasis.resize(info.numberPatches);
        m_KCSolve.resize(info.numberPatches);
        // m_Kreg.resize(info.numberPatches);
        //      m_KM.resize(info.numberPatches);
        break;
    default:
        break;
    }

    if(m_IETIoptions.opt.getSwitch("SaddlePoint"))
        m_Kii.resize(info.numberPatches);
}


template<class T>
void gsIETIAssembler<T>::reserveSpace(std::vector<gsMatrix<T> >& dirDofs, std::vector<gsMatrix<T> >& rhs_loc)
{
    size_t sumD, sumR;

    for(size_t np=0; np<info.numberPatches;np++)
    {
        sumD=sumR=0;
        for(size_t c=0; c<info.cDim;c++)
        {
            sumD+= m_locDofsMapper[np][c].boundarySize();
            sumR+= m_locDofsMapper[np][c].freeSize();
        }


        dirDofs[np].setZero(sumD,info.nRhs);
        rhs_loc[np].setZero(sumR,info.nRhs);

    }
}

template <typename T>
void gsIETIAssembler<T>::giveAssembledMatrices(std::vector<gsSparseMatrix<T> > &matrices, const gsMatrix<T>& rhs)
{
    m_useGivenMats = true;
    m_tempMatrix.resize(matrices.size());
    m_tempRhs.resize(matrices.size());

    m_tempMatrix.swap(matrices);
    for(size_t j =0; j<matrices.size();++j)
        extractPatch(j,rhs,m_tempRhs[j]);
}

template<typename T>
void gsIETIAssembler<T>::setNewRhs(const gsMatrix<T>& rhs)
{
    m_rhs_p.setZero(info.dofTotalP,info.nRhs);

    std::vector<gsMatrix<T> > rhs_loc(m_patches.nPatches());
    for(size_t j =0; j<m_patches.nPatches();++j)
    {
        extractPatch(j,rhs,rhs_loc[j]);
    //    gsInfo<<"Rhs on Patch: "<<j<<"\n "<<rhs_loc[j].transpose()<<"\n"<<std::flush;
        assembleRhs(rhs_loc, j);
    }
    assembleRhsFreeElim();
}

template <typename T>
void gsIETIAssembler<T>::extractPatch(size_t np, const gsMatrix<T>& rhs, gsMatrix<T>& rhsLocal) const
{
    const gsDofMapper& locMapper = m_locDofsMapper[np][0];
    const gsDofMapper& globMapper = m_assembler->system().colMapper(0);
    rhsLocal.setZero(locMapper.freeSize(),rhs.cols());

    for(size_t i=0; i<(size_t)locMapper.size();++i)
    {
        unsigned glIdxAss = globMapper.index(i,np);
        unsigned glIdx =  m_stdMapper[0].index(i,np);
        unsigned idx = locMapper.index(i);
        if(!globMapper.is_free_index(glIdx)) continue;

        if(m_globIsBoundIndex[np][idx])
        {
            const std::set<patchDof>& set = (m_globalConnectionTable.find(glIdx))->second;
           // for(std::set<patchDof>::const_iterator it = set.begin();it!=set.end();++it)
            rhsLocal.row(idx) = rhs.row(glIdxAss)/set.size();
        }
        else
            rhsLocal.row(idx) = rhs.row(glIdxAss);
    }
}

template<class T>
void gsIETIAssembler<T>::assembleLocal(gsAssembler<T>* A, gsSparseMatrix<T>& matrix, gsMatrix<T>& rhs, gsMatrix<T>& dirDofs, size_t np, const gsMultiPatch<T> &curSol)
{
    gsStopwatch time;
    // -- Assembling
    if(!m_useGivenMats)
    {
        gsBasisRefs<T> mBsRefs(m_basis,np);

        gsPde<T>* pde = m_assembler->pde().restrictToPatch(np);
        A->initialize(*pde,mBsRefs,m_assembler->options());
        A->system().rhs().conservativeResize(Eigen::NoChange, rhs.cols());

        if(m_IETIoptions.opt.getSwitch("NonlinearMode"))
        {
            for(size_t c =0; c<info.cDim;c++)
                A->computeDirichletDofs(c);
            A->homogenizeFixedDofs(-1);
            A->assemble(curSol.patch(np));

        }
        else
            A->assemble();

        const gsSparseSystem<T> & system = A->system();
        if(system.symmetry())
            matrix = A->matrix().template selfadjointView<Lower>();
        else
            matrix = A->matrix();

        rhs.topRows(A->numDofs()) = system.rhs();
        // gsDebugVar(matrix.toDense());
        //  gsDebugVar(rhs.transpose());

        index_t i = 0;
        for(size_t c=0; c<info.cDim;c++)
        {

            if(A->fixedDofs(c).rows() != 0)
                dirDofs.block(i, 0, A->fixedDofs(c).rows(), A->fixedDofs(c).cols()) = A->fixedDofs(c);

            i+=A->fixedDofs(c).rows();
        }
        delete pde;
    }
    else
    {
        matrix.swap(m_tempMatrix[np]);
        rhs.topRows(m_tempRhs[np].rows()) = m_tempRhs[np];
        dirDofs.setZero();
    }



#pragma omp critical
    m_timings.totalAssemble[np]+= time.stop();
}


template<class T>
void gsIETIAssembler<T>::assemble(const gsMultiPatch<T> &curSol)
{
    assembleInit();
    bool assembleRHS = !m_IETIoptions.opt.getSwitch("NoAssemblyOfRHS");
    std::vector<typename gsSparseMatrix<T>::Ptr > matrices(info.numberPatches);
    std::vector<typename gsSparseMatrix<T>::Ptr > matrices2(info.numberPatches);
    std::vector<gsMatrix<T> > rhs_loc(info.numberPatches);
    std::vector<gsMatrix<T> > spp_loc(info.numberPatches);

    gsStopwatch time, time2;
    Eigen::setNbThreads(1);

    reserveSpace(m_dirDofs,rhs_loc);

    if(!m_IETIoptions.opt.getSwitch("ExtraTiming")) //The Fast One
    {
#pragma omp parallel
        {
            gsSparseMatrix<T> matrix, matrix2;

            gsAssembler<T>* A = m_assembler->create();
#pragma omp for  schedule(static,1) nowait
            for (size_t np=0; np < info.numberPatches; ++np )
            {
                assembleLocal(A, matrix, rhs_loc[np],m_dirDofs[np],np,curSol);

                // -- Reordering
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                {
                    matrix2 = matrix;
                    makeReorderingPrimalRem(matrix2,np);
                }
                makeReorderingBoundInt(matrix,np);

                //For Preconditioning and processing the rhs and solution ICO minimum energy
                assembleKiiKibKbb(matrix, np);

                // The main matrices for applying the system matrix
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                    assembleKrrKrpKpp(matrix2,np);
                else
                {
                    assembleC(np);
                    assembleLUofKC(matrix,np);
                }


                //assemble local Spp
                assembleSppLoc(spp_loc[np], np);

                //assemble local rhs
                if(assembleRHS) assembleRhs(rhs_loc, np);

            } //end-for

            //free memory
            delete A;

        }//end-parallel

    }//end if
    else // The slow one
    {
        for (size_t np=0; np < info.numberPatches; ++np )
        {
            matrices[np] = gsSparseMatrix<T>().moveToPtr();
            matrices2[np] = gsSparseMatrix<T>().moveToPtr();
        }
#pragma omp parallel
        {
            gsAssembler<T>* A = m_assembler->create();

#pragma omp for  schedule(static, 1)
            for (size_t np=0; np < info.numberPatches; ++np )
                assembleLocal(A, *matrices[np], rhs_loc[np],m_dirDofs[np], np,curSol);

            delete A;


            printTime(time,"Time for assemling all patch local matrices: ");

            // -- Reordering
#pragma omp for  schedule(static, 1)
            for (size_t np=0; np < info.numberPatches; ++np )
            {
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                {
                    *matrices2[np] = *matrices[np];
                    makeReorderingPrimalRem(*matrices2[np],np);
                }
                makeReorderingBoundInt(*matrices[np],np);
            }

            printTime(time,"Time for reordering all patch local matrices: ");

            //------------------------------------------------------------------------------/
            //For Preconditioning and processing the rhs and solution ICO minimum energy


#pragma omp for  schedule(static, 1)
            for (size_t np=0; np < info.numberPatches; ++np )
                assembleKiiKibKbb(*matrices[np], np);

            printTime(time,"Time for calculating LU of Kii : ");


            // The main matrices for applying the system matrix
            if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
            {
#pragma omp for schedule(static, 1)
                for (size_t np=0; np < info.numberPatches; ++np )
                    assembleKrrKrpKpp(*matrices2[np],np);

                printTime(time,"Time for computing LU factorization of Krr : ");
            }
            else
            {
#pragma omp for  schedule(static, 1)
                for (size_t np=0; np < info.numberPatches; ++np )
                    assembleC(np);

                printTime(time,"Time for computing C : ");

#pragma omp for  schedule(static, 1)
                for (size_t np=0; np < info.numberPatches; ++np )
                    assembleLUofKC(*matrices[np],np);

                printTime(time,"Time for computing LU factorization of KC : ");
            }


#pragma omp for  schedule(static, 1)
            for (size_t np=0; np < info.numberPatches; ++np )
                assembleSppLoc(spp_loc[np], np);

            printTime(time,"Time for assemling local Spp: ");

#pragma omp for  schedule(static, 1)
            for (size_t np=0; np < info.numberPatches; ++np )
                if(assembleRHS) assembleRhs(rhs_loc, np);

            printTime(time,"Time for assemling rhs: ");

        }
    }
    printTime(time2,"Time for doing parallel assembling part: ");

    //------------Do Serial stuff--------------
    Eigen::setNbThreads(0);

    //1. Rhs
    if(assembleRHS) assembleRhsFreeElim();

    //2. Spp
    assembleSpp(spp_loc);

    if(m_IETIoptions.opt.getSwitch("CalcRhsNorm"))
        calcRhsNorm(rhs_loc);

    printTime(time2,"Time for doing serial assembling part: ");
}

template<class T>
void gsIETIAssembler<T>::constructKC(gsSparseMatrix<T> & KC, const gsSparseMatrix<T> &C, IETILocalSolver::solvers type, size_t np)
{
    switch (type) {
    case IETILocalSolver::Direct:
    {
        gsVector<int> nonZerosPerCol(info.dofsB[np]+info.dofsI[np]+ info.dofsP[np]);
        KC.conservativeResize(info.dofsB[np]+info.dofsI[np]+ info.dofsP[np],info.dofsB[np]+info.dofsI[np]+ info.dofsP[np]);
        for(index_t i=0; i<KC.outerSize();i++)
            if(i<info.dofsB[np])
                nonZerosPerCol(i) = cast<double,int>((KC.innerVector(i).nonZeros()+C.innerVector(i).nonZeros())*1.333);
            else if(i>=info.dofsB[np]+info.dofsI[np])
                nonZerosPerCol(i) = cast<double,int>(info.dofsB[np]*1.333);
            else
                nonZerosPerCol(i) = cast<double,int>(KC.innerVector(i).nonZeros()*1.333);

        KC.reserve(nonZerosPerCol);

        for (int k=0; k<C.outerSize(); ++k)
            for (typename gsSparseMatrix<T>::InnerIterator it(C,k); it; ++it)
            {
                KC.insert(info.dofsB[np]+info.dofsI[np]+it.row(),it.col()) = it.value();
                KC.insert(it.col(),info.dofsB[np]+info.dofsI[np]+it.row()) = it.value();
            }
        break;
    }
    case IETILocalSolver::FastDiagonalization:
    case IETILocalSolver::Multigrid:
    {
        /*
        for(index_t i=0; i<KC.outerSize();i++)
            if((unsigned)i<info.dofsB[np]+info.dofsI[np])
            {
                int nC=0;
                if(m_globIsBoundIndex[np][i])
                    nC=C.col(m_glob2BoundInteriorIndex[np][i]).nonZeros();
                nonZerosPerCol(i) = cast<double,int>((KC.innerVector(i).nonZeros()+nC)*1.333);
            }
            else
                nonZerosPerCol(i) = cast<double,int>(info.dofsB[np]*1.333);

        KC.reserve(nonZerosPerCol);

        for (index_t k=0; k<C.outerSize(); ++k)
            for (typename gsSparseMatrix<T>::InnerIterator it(C,k); it; ++it)
            {
                int c=getComp(np,it.col());
                index_t col = m_locDofsMapper[np][c].index(m_boundDofs[np][it.col()]);
                KC.insert(info.dofsB[np]+info.dofsI[np]+it.row(),col) = it.value();
                KC.insert(col,info.dofsB[np]+info.dofsI[np]+it.row()) = it.value();
            }
            */
        m_exC[np].resize(m_C[np].rows(),info.dofsB[np]+info.dofsI[np]);
        gsVector<int> nonZerosPerCol(info.dofsB[np]+info.dofsI[np]);
        for(index_t i=0; i<info.dofsB[np]+info.dofsI[np];i++)
        {
            int nC=0;
            if(m_globIsBoundIndex[np][i])
                nC=m_C[np].col(m_glob2BoundInteriorIndex[np][i]).nonZeros();
            nonZerosPerCol(i) = cast<double,int>(nC*1.333);
        }

        m_exC[np].reserve(nonZerosPerCol);
        for (index_t k=0; k<m_C[np].outerSize(); ++k)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_C[np],k); it; ++it)
            {
                int c=getComp(np,it.col());
                index_t col = m_locDofsMapper[np][c].index(m_boundDofs[np][it.col()]);
                m_exC[np].insert(it.row(),col) = it.value();
            }

    }
        break;
    }

    //cd ..  gsInfo<<"KC:\n"<<KC.toDense()<<"\n\n";
    //gsSparseMatrix<T> Ct = C.transpose();
    // KC.bottomLeftCorner(info.dofsP[np],info.dofsB[np])=C;
    // KC.topRightCorner(info.dofsB[np],info.dofsP[np])=Ct;

    KC.makeCompressed();
}

template<class T>
void gsIETIAssembler<T>::assembleKiiKibKbb(gsSparseMatrix<T> &matrix, size_t np)
{
    gsVector<index_t> splitting(3);
    splitting[0]= info.dofsB[np],splitting[1]= info.dofsI[np];splitting[2]=matrix.rows()-info.dofsI[np]-info.dofsB[np];

    typename gsSparseMatrix<T>::BlockView view = matrix.blockView(splitting, splitting);

    m_Kbb[np] = gsSparseMatrix<T>(view(0,0)).moveToPtr();
    m_Kib[np] = gsSparseMatrix<T>(view(1,0)).moveToPtr();
    m_Kbi[np] = gsSparseMatrix<T>(view(0,1)).moveToPtr();
    m_LU_Kii[np]=NULL;

    //  gsInfo<<"patch: "<<np<<"\n Kii: \n"<<view(1,1).toDense()<<"\n\n Kbb: \n"<<view(0,0)<<"\n\n Kib: \n"<<view(1,0)<<"\n\n Kbi: \n"<<view(0,1)<<"\n\n";

    if(m_IETIoptions.opt.getSwitch("EnforceRescaling") || m_IETIoptions.opt.getSwitch("SaddlePoint") ||
            m_IETIoptions.KiiSolver == IETILocalSolver::FastDiagonalization   || m_IETIoptions.KiiSolver == IETILocalSolver::Multigrid )
        m_Kii[np]=memory::make_shared(new gsSparseMatrix<T>(view(1,1))); // std::make_shared<gsSparseMatrix<T> >(view(1,1));//memory::make_shared(new gsSparseMatrix<T>(view(1,1)));


    if(info.dofsI[np] ==0)
        return;

    if(m_IETIoptions.opt.getSwitch("EnforceRescaling") )
    {
        gsSparseMatrix<T> Kii = view(1,1);
        m_scalingsKii[np].computeRef(Kii);
    }

    switch (m_IETIoptions.KiiSolver)
    {
    case IETILocalSolver::Direct:
        m_LU_Kii[np]=new sparseSPDfact(view(1,1));
        if(!m_LU_Kii[np]->succeed())
	{
            gsWarn<<"Sparse LLT not succeeded on patch "<<np<<"\n!";
	}
        break;
    case IETILocalSolver::FastDiagonalization:
    {
        gsBoundaryConditions<T> bc = m_bc_loc[np];
        gsConstantFunction<> zero(0.0, info.dim);
        for(size_t c=0; c<info.cDim;++c)
            for(std::vector<patchSide>::iterator it = m_patchISides[np].begin();it!=m_patchISides[np].end();it++)
                if(it->patch== static_cast<index_t>(np))
                    bc.addCondition(0,it->side(),condition_type::dirichlet,zero,c);
	if(m_IETIoptions.opt.getSwitch("MG_Kii.BasedOnRankOneCorrection"))
	   m_inexact_Kii_BS[np] = gsPatchPreconditionersCreator2<T>::fastDiagonalizationWithGeoOp(m_patches.patch(np),m_basis.front()[np],bc,0);
	else
	   m_inexact_Kii_BS[np] = gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(m_basis.front()[np],bc,gsAssembler<T>::defaultOptions(),0);
	
	if(m_maxEigKii[np]==-1)
	{
        m_maxEigKii[np] = powerIteration<T>(makeMatrixOp(m_Kii[np]),m_inexact_Kii_BS[np]);
	}

        m_KiiSolve[np]=memory::make_unique(new gsConjugateGradient<>(*m_Kii[np], m_inexact_Kii_BS[np]));
        m_KiiSolve[np]->setOptions(m_IETIoptions.opt.getGroup("MG_Kii.RHS"));

        break;
    }
    case IETILocalSolver::Multigrid:
    {
        gsBoundaryConditions<T> bc = m_bc_loc[np];
        gsConstantFunction<> zero(0.0, info.dim);
        for(size_t c=0; c<info.cDim;++c)
            //for(std::vector<patchSide>::iterator it = m_patchISides[np].begin();it!=m_patchISides[np].end();it++)
            for(std::vector<patchSide>::iterator it = m_patchAveragesSides[np].begin();it!=m_patchAveragesSides[np].end();it++)
                if(it->patch== static_cast<index_t>(np))
                    bc.addCondition(0,it->side(),condition_type::dirichlet,zero,c);

        if(!m_IETIoptions.opt.getSwitch("MG_Kii.UseP1Projection"))
        {
            gsMultiPatch<T>mp(m_patches[np]);
            m_inexact_Kii_MG[np]=
                    getMultiGrid<T>(gsMultiBasis<>(m_basis[0][np]),
                    m_Kii[np],bc,m_IETIoptions.opt.getGroup("MG_Kii"), &mp);
        }
        else
        {
            gsBasis<T>& P1Basis = (gsBasis<T>&)*m_basis[0][np].clone().release();
            P1Basis.degreeDecrease(P1Basis.maxDegree()-1);

            memory::unique_ptr< gsPde<T> > pde(m_assembler->pde().restrictToPatch(np));
            gsAssembler<T>* A = m_assembler->create();

            pde->boundaryConditions() = bc;
            A->initialize(*pde,P1Basis,m_assembler->options());
            A->assemble();

            gsMultiPatch<T>mp(m_patches[np]);
            gsMultiGridOp<>::Ptr mg=getMultiGrid<T>(gsMultiBasis<>(P1Basis),
                    give(A->system().matrix()),bc,m_IETIoptions.opt.getGroup("MG_Kii"), &mp);


            gsSparseMatrix<real_t> Mrect;
            assembleParameterMassForTensorProductSpace<>( m_basis[0][np], P1Basis, bc, Mrect);


            std::vector< gsSparseMatrix<real_t> > Morigs(m_patches[np].geoDim());
            for ( index_t i=0; i<m_patches[np].geoDim(); ++i )
            {
                assembleParameterMass(m_basis[0][np].component(i), Morigs[i]);
                handleDirichletConditions(Morigs[i],bc,1+2*i,2+2*i);
            }

            gsSparseMatrix<> B_tilde_full;
            if (false)
            {
                std::vector< gsSparseMatrix<> > B_tilde;
                std::vector< gsSparseMatrix<> > B_l2compl;
                std::vector< gsSparseMatrix<> > B_tilde2;
                constructTildeSpaceBasis( m_basis[0][np], bc, B_tilde, B_l2compl );

                for ( index_t i=0; i<m_patches[np].geoDim(); ++i )
                    Morigs[i] = B_tilde[i].transpose() * Morigs[i] * B_tilde[i];

                B_tilde2.resize(m_patches[np].geoDim());
                for ( index_t i=0; i<m_patches[np].geoDim(); ++i )
                    B_tilde2[i] = B_tilde[m_patches[np].geoDim()-1-i];

                B_tilde_full = kroneckerProduct( B_tilde2 );
                Mrect = B_tilde_full.transpose() * Mrect;
            }

            std::vector< gsLinearOperator<>::Ptr> massInvs(m_patches[np].geoDim());
            for ( index_t i=0; i<m_patches[np].geoDim(); ++i )
                massInvs[m_patches[np].geoDim()-1-i] = makeSparseCholeskySolver(Morigs[i]);
            gsKroneckerOp<>::Ptr massInv = gsKroneckerOp<>::make( massInvs );

            // the call of assembleParameterMassInverseForTensorProductSpace has to be splitted in order to also allow the Tilde space projection
            // gsKroneckerOp<>::Ptr massInv;
            //  assembleParameterMassInverseForTensorProductSpace(m_basis[0][np],bc,massInv);

            gsLinearOperator<>::Ptr massSmoother;
            if(m_IETIoptions.opt.getSwitch("MG_Kii.TL.BasedOnRankOneCorrection"))
                massSmoother = makeSubspaceCorrectedMassSmootherOperator(m_patches[np],m_basis[0][np], m_IETIoptions.opt.getReal("MG_Kii.TL.Damping"), bc);
            else
                massSmoother = makeSubspaceCorrectedMassSmootherOperator(m_basis[0][np], m_IETIoptions.opt.getReal("MG_Kii.TL.Damping"), bc);

            gsPreconditionerOp<>::Ptr ps = gsCompositePrecOp<>::make(
                        gsPreconditionerFromOp<>::make( makeMatrixOp(  m_Kii[np] ),massSmoother,m_IETIoptions.opt.getReal("MG_Kii.TL.OuterDamping") ),
                        makeGaussSeidelOp( m_Kii[np] )
                        );
            gsTwoLevel::Ptr tw;
            if(false)
                tw = gsTwoLevel::make(
                            makeMatrixOp(  m_Kii[np]),
                            ps,
                            // makeSubspaceCorrectedMassSmootherOperator(m_basis[0][np], m_options.Kii_optionsMG.damping, bc),
                            give(Mrect),
                            massInv,
                            give(B_tilde_full),
                            mg);
            else
                tw = gsTwoLevel::make(
                            makeMatrixOp(  m_Kii[np]),
                            ps,
                            //makeSubspaceCorrectedMassSmootherOperator(m_basis[0][np], m_options.Kii_optionsMG.damping, bc),
                            give(Mrect),
                            massInv,
                            mg);
            tw->setOptions(m_IETIoptions.opt.getGroup("MG_Kii.TL"));
            m_inexact_Kii_MG[np] = tw;
            delete A;
            delete & P1Basis;

        }

        m_KiiSolve[np]=memory::make_unique(new gsConjugateGradient<>(*m_Kii[np], m_inexact_Kii_MG[np]));
        m_KiiSolve[np]->setOptions(m_IETIoptions.opt.getGroup("MG_Kii.RHS"));

        break;
    }
    }


    if(false && (np==1 || np==2))
    {
        gsInfo<<"Mat Kii:\n"<<m_Kii[np]->toDense()<<"\n\n";
        gsMatrix<T> mat;
        m_inexact_Kii_MG[np]->toMatrix(mat);
        gsInfo<<"MG: Kii: \n"<<mat<<"\n\n";
    }
}

template<class T>
void gsIETIAssembler<T>::getSchur(gsMatrix<T> &schur, size_t np)
{
    /*
    gsMatrix<T> d(m_Kbb[np].rows(),1);
    d.array()= m_Kbb[np].diagonal().array();
    gsSparseMatrix<T> d_inv(m_Kbb[np].rows(),m_Kbb[np].cols());
    d_inv.setIdentity();
    for(index_t i =0;i<m_Kbb[np].rows();i++)
        d_inv(i,i)/=(d(i));

    schur = (m_C[np]*d_inv*m_C[np].transpose().toDense());
*/
    /*
    gsMatrix<T> d(m_K[np].rows(),1);
    d.array()= m_K[np].diagonal().array();
    gsSparseMatrix<T> d_inv(m_K[np].rows(),m_K[np].cols());
    d_inv.setIdentity();
    for(size_t i =0;i<m_K[np].rows();i++)
        d_inv(i,i)/=d(i);

    schur = (m_exC[np]*d_inv*m_exC[np].transpose().toDense());
*/
    /*
    sparseSPDfact Kbb_LU(m_Kbb[np]);
    schur = (m_C[np]*Kbb_LU.solve(m_C[np].transpose().toDense()));
*/

    /*
    gsMatrix<T> d(m_Kii[np].rows(),1);
    d.array()= m_Kii[np].diagonal().array();
    gsSparseMatrix<T> d_inv(m_Kii[np].rows(),m_Kii[np].cols());
    d_inv.setIdentity();
    for(size_t i =0;i<m_Kii[np].rows();i++)
        d_inv(i,i)/=d(i);
    sparseLUfact Kbb_LU(m_Kbb[np]-m_Kib[np].transpose()*d_inv*m_Kib[np]);
    schur = (m_C[np]*Kbb_LU.solve(m_C[np].transpose().toDense()));
*/

    /*
    sparseSPDfact kiiinv(m_Kii[np]);
           // gsGenericAssembler<T> genAss(m_patches.patch(np), m_basis.front()[np],gsAssembler<>::defaultOptions(),&m_bc_loc[np]);
    Eigen::PartialPivLU<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > Sb(m_Kbb[np].toDense()-m_Kib[np].transpose().toDense()*kiiinv.solve(m_Kib[np].toDense()));
   schur = m_C[np]*Sb.solve(m_C[np].transpose().toDense());
*/
    /*
    gsGenericAssembler<T> genAss(m_patches.patch(np), m_basis.front()[np],gsAssembler<>::defaultOptions(),&m_bc_loc[np]);
    sparseSPDfact kinv(m_K[np]+1.e-5*genAss.assembleMass());
    schur = m_exC[np]*kinv.solve(m_exC[np].transpose().toDense());
    */


    gsMatrix<T> col;
    schur.setZero(info.dofsP[np],info.dofsP[np]);
    gsPreconditionerOp<T>* mg =  dynamic_cast< gsPreconditionerOp<T>*>((m_inexact_KC[np]->getOperator(0,0).get()));
    if(mg !=NULL)
    {
      /*
        for(size_t i=0; i<info.dofsP[np];++i)
        {
            col.setZero(m_exC[np].cols(),1);
            for( int p=0; p<1;++p )
                mg->step(m_exC[np].row(i).transpose().toDense(), col);
            schur.col(i) = m_exC[np]*col;
        }
        */
              for(index_t i=0; i<info.dofsP[np];++i)
        {
        col.setZero(m_exC[np].cols(),1);
        m_inexact_KC[np]->getOperator(0,0)->apply(m_exC[np].row(i).transpose().toDense(), col);
        schur.col(i) = m_exC[np]*col;
	}
        
    }
    else
    {
        /*
        gsSparseMatrix<T> d_inv(m_Kbb[np].rows(),m_Kbb[np].cols());
        d_inv.setIdentity();
        for(index_t i =0;i<m_Kbb[np].rows();i++)
            d_inv(i,i)/=m_Kbb[np](i,i);

        schur = (m_C[np]*d_inv*m_C[np].transpose().toDense());
        */
        for(index_t i=0; i<info.dofsP[np];++i)
        {
            col.setZero(m_exC[np].cols(),1);
            m_inexact_KC[np]->getOperator(0,0)->apply(m_exC[np].row(i).transpose().toDense(), col);
            schur.col(i) = m_exC[np]*col;
        }

    }
  //  gsInfo<<"Finished schur\n";

}

template<class T>
void gsIETIAssembler<T>::setupMGforKC(size_t np)
{
    typename gsPreconditionerOp<T>::Ptr K_Precond;
    gsBoundaryConditions<T> bc = m_bc_loc[np];


    m_exC[np].resize(m_C[np].rows(),info.dofsB[np]+info.dofsI[np]);
    gsVector<int> nonZerosPerCol(info.dofsB[np]+info.dofsI[np]);
    for(index_t i=0; i<info.dofsB[np]+info.dofsI[np];i++)
    {
        int nC=0;
        if(m_globIsBoundIndex[np][i])
            nC=m_C[np].col(m_glob2BoundInteriorIndex[np][i]).nonZeros();
        nonZerosPerCol(i) = cast<double,int>(nC*1.333);
    }

    m_exC[np].reserve(nonZerosPerCol);
    for (index_t k=0; k<m_C[np].outerSize(); ++k)
        for (typename gsSparseMatrix<T>::InnerIterator it(m_C[np],k); it; ++it)
        {
            int c=getComp(np,it.col());
            index_t col = m_locDofsMapper[np][c].index(m_boundDofs[np][it.col()]);
            m_exC[np].insert(it.row(),col) = it.value();
        }

    if(!m_IETIoptions.opt.getSwitch("MG_KC.UseP1Projection"))
    {
        gsMultiPatch<T>mp (m_patches[np]);
        if(m_bc_loc[np].dirichletSides().size()==0)
        {
            gsStopwatch time;
            gsSparseMatrix<T> M;
            assembleParameterMassForTensorProductSpace(m_basis[0][np],m_bc_loc[np],M);
            m_timings.assemblingMass[np]+=time.stop();
            T H = (m_patches.patch(np).coefAtCorner(boxCorner::getFirst(m_patches.parDim())) - m_patches.patch(np).coefAtCorner(boxCorner::getLast(m_patches.parDim()))).norm();
	  //  T H = 1;
            K_Precond= getMultiGrid<T>(gsMultiBasis<>(m_basis[0][np]),
                    *m_K[np]+math::pow(H,m_patches.parDim())*m_IETIoptions.opt.getReal("Regularization")*M,bc,m_IETIoptions.opt.getGroup("MG_KC"), &mp);
        }
        else
            K_Precond= getMultiGrid<T>(gsMultiBasis<>(m_basis[0][np]),
                    m_K[np],bc,m_IETIoptions.opt.getGroup("MG_KC"), &mp);

    }
    else
    {

        //TODO: use coarsening or set it outside IETI algorithm
        gsBasis<T>& P1Basis = (gsBasis<T>&)*m_basis[0][np].clone().release();
        P1Basis.degreeDecrease(P1Basis.maxDegree()-1);

        memory::unique_ptr< gsPde<T> > pde(m_assembler->pde().restrictToPatch(np));
        gsAssembler<T>* A = m_assembler->create();
        pde->boundaryConditions() = bc;
        A->initialize(*pde,P1Basis,m_assembler->options());
        A->assemble();


        gsMultiGridOp<>::Ptr mg;
        gsMultiPatch<T>mp(m_patches[np]);
        if(m_bc_loc[np].dirichletSides().size()==0)
        {
            // gsInfo<<"Mass: \n"<<genAss.assembleMass().toDense()<<"\n \n";
            gsStopwatch time;
            gsSparseMatrix<T> M;

            assembleParameterMassForTensorProductSpace(P1Basis,m_bc_loc[np],M);
            m_timings.assemblingMass[np]+=time.stop();
            T H = (m_patches.patch(np).coefAtCorner(boxCorner::getFirst(m_patches.parDim())) - m_patches.patch(np).coefAtCorner(boxCorner::getLast(m_patches.parDim()))).norm();
	   // T H = 1;
            mg= getMultiGrid<T>(gsMultiBasis<>(P1Basis),
                    A->matrix()+math::pow(H,m_patches.parDim())*m_IETIoptions.opt.getReal("Regularization")*M,bc,m_IETIoptions.opt.getGroup("MG_KC"), &mp);
        }
        else
            mg= getMultiGrid<T>(gsMultiBasis<>(P1Basis),
                    give(A->system().matrix()),bc,m_IETIoptions.opt.getGroup("MG_KC"), &mp);

        gsKroneckerOp<>::Ptr massInv;
        assembleParameterMassInverseForTensorProductSpace(m_basis[0][np],bc,massInv);

        gsSparseMatrix<real_t> Mrect;
        assembleParameterMassForTensorProductSpace<>( m_basis[0][np], P1Basis, bc, Mrect);

        gsPreconditionerOp<>::Ptr ps;

        gsSparseMatrix<T> M;
        gsLinearOperator<>::Ptr massSmoother;
        if(m_IETIoptions.opt.getSwitch("MG_KC.TL.BasedOnRankOneCorrection"))
            massSmoother = makeSubspaceCorrectedMassSmootherOperator(m_patches[np],m_basis[0][np], m_IETIoptions.opt.getReal("MG_KC.TL.Damping"), bc);
        else
            massSmoother = makeSubspaceCorrectedMassSmootherOperator(m_basis[0][np], m_IETIoptions.opt.getReal("MG_KC.TL.Damping"), bc);


        if(m_bc_loc[np].dirichletSides().size()==0)
        {
            T H = (m_patches.patch(np).coefAtCorner(boxCorner::getFirst(m_patches.parDim())) - m_patches.patch(np).coefAtCorner(boxCorner::getLast(m_patches.parDim()))).norm();
	   // T H = 1;
            assembleParameterMassForTensorProductSpace(m_basis[0][np],m_bc_loc[np],M);
            typename gsSparseMatrix<T>::Ptr Kreg= memory::make_shared(new gsSparseMatrix<T>(*m_K[np]+ math::pow(H,m_patches.parDim())*m_IETIoptions.opt.getReal("Regularization")*M));
            ps= gsCompositePrecOp<>::make(
                        gsPreconditionerFromOp<>::make( makeMatrixOp(  Kreg ), massSmoother,m_IETIoptions.opt.getReal("MG_KC.TL.OuterDamping")),
                        makeGaussSeidelOp( Kreg )
                        );
            K_Precond= gsTwoLevel::make(
                        makeMatrixOp(  Kreg ),
                        ps,
                        give(Mrect),
                        massInv,
                        mg);
        }
        else
        {
            ps= gsCompositePrecOp<>::make(
                        gsPreconditionerFromOp<>::make( makeMatrixOp(  m_K[np] ), massSmoother,m_IETIoptions.opt.getReal("MG_KC.TL.OuterDamping")),
                        makeGaussSeidelOp( m_K[np] )
                        );
            K_Precond= gsTwoLevel::make(
                        makeMatrixOp(  m_K[np]),
                        ps,
                        give(Mrect),
                        massInv,
                        mg);
        }
        K_Precond->setOptions(m_IETIoptions.opt.getGroup("MG_KC.TL"));
        // tw->setCoarseDamping(m_options.KC_optionsMG.coarseDamping);

        delete A;
        delete &P1Basis;
    }

    m_inexact_KC[np] = gsBlockOp<>::make(2,2);
    m_inexact_KC[np]->addOperator(0,0,K_Precond);

    getSchur(m_S[np],np);
    m_inexact_KC[np]->addOperator(1,1,makePartialPivLUSolver(m_S[np]));
}

template<class T>
void gsIETIAssembler<T>::assembleLUofKC(gsSparseMatrix<T> &matrix, size_t np)
{
    const gsSparseMatrix<T>& C = getC(np);

    if(matrix.cols() ==0)
        return;

    m_LU_KC[np]=NULL;
    switch (m_IETIoptions.KCSolver)
    {
    case IETILocalSolver::Direct:
    {

        gsSparseMatrix<T>& KC = matrix;
        constructKC(KC,C,IETILocalSolver::Direct,np);
        //   gsInfo<<"KC: \n"<<KC.toDense()<<"\n\n";
        if(m_IETIoptions.opt.getSwitch("EnforceRescaling"))
            m_scalingsKC[np].computeRef(KC);
        m_LU_KC[np]=new sparseLUfact(KC);
        if(!m_LU_KC[np]->succeed())
            gsWarn<<"Sparse LU not succeeded on patch "<<np<<"\n!"<<std::flush;
        break;
    }
    case IETILocalSolver::FastDiagonalization:
    {
        gsBoundaryConditions<T> bc = m_bc_loc[np];

        constructKC(*m_K[np],C,IETILocalSolver::FastDiagonalization,np);

        m_inexact_KC[np] = gsBlockOp<>::make(2,2);
        T reg = m_bc_loc[np].dirichletSides().size()==0 ?m_IETIoptions.opt.getReal("Regularization") : 0;
	
	if(m_IETIoptions.opt.getSwitch("MG_KC.BasedOnRankOneCorrection"))
	    m_inexact_KC[np]->addOperator(0,0,gsPatchPreconditionersCreator2<T>::fastDiagonalizationWithGeoOp(m_patches.patch(np),m_basis.front()[np],bc,reg));
	else
	    m_inexact_KC[np]->addOperator(0,0,gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(m_basis.front()[np],bc,gsAssembler<T>::defaultOptions(),reg));

        getSchur(m_S[np],np);
        m_inexact_KC[np]->addOperator(1,1,makePartialPivLUSolver(m_S[np]));
	if(m_maxEigKC[np]==-1)
	{
        m_maxEigKC[np] = powerIteration<T>(makeMatrixOp(m_K[np]),m_inexact_KC[np]->getOperator(0,0));
	}

        m_KC[np] = gsBlockOp<>::make(2,2);
        m_KC[np]->addOperator(0,0,makeMatrixOp(*m_K[np]));
        m_KC[np]->addOperator(1,0,makeMatrixOp(m_exC[np]));
        m_KC[np]->addOperator(0,1,makeMatrixOp(m_exC[np].transpose()));

        m_KCBasis[np]=gsBramblePasciakCG<gsBPCG_Types::SchoeberlZulehner>::make(m_KC[np],m_inexact_KC[np]->getOperator(0,0),m_inexact_KC[np]->getOperator(1,1));
        m_KCBasis[np]->setOptions(m_IETIoptions.opt.getGroup("MG_KC.Basis"));


        m_KCSolve[np]=memory::make_shared(new gsConjugateGradient<>(*m_K[np],m_inexact_KC[np]->getOperator(0,0)));
        m_KCSolve[np]->setOptions(m_IETIoptions.opt.getGroup("MG_KC.Solve"));
        break;
    }
    case IETILocalSolver::Multigrid:
    {
        setupMGforKC(np);

        m_KC[np] = gsBlockOp<>::make(2,2);
        m_KC[np]->addOperator(0,0,makeMatrixOp(*m_K[np]));
        m_KC[np]->addOperator(1,0,makeMatrixOp(m_exC[np]));
        m_KC[np]->addOperator(0,1,makeMatrixOp(m_exC[np].transpose()));

        m_KCBasis[np]=gsBramblePasciakCG<gsBPCG_Types::SchoeberlZulehner>::make(m_KC[np],m_inexact_KC[np]->getOperator(0,0),m_inexact_KC[np]->getOperator(1,1));
        m_KCBasis[np]->setOptions(m_IETIoptions.opt.getGroup("MG_KC.Basis"));


        m_KCSolve[np]=memory::make_shared(new gsConjugateGradient<>(*m_K[np],m_inexact_KC[np]->getOperator(0,0)));
        m_KCSolve[np]->setOptions(m_IETIoptions.opt.getGroup("MG_KC.Solve"));
        break;

    }

    }
}


template<class T>
void gsIETIAssembler<T>::assembleSppLoc(gsMatrix<T>& spp_loc, size_t np)
{
    if(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy") )
    {
        gsMatrix<T> rhs, solution;
        index_t rowPhi=info.dofsB[np];
        if(m_IETIoptions.opt.getSwitch("SaddlePoint"))
            rowPhi= info.dofsB[np]+info.dofsI[np];

        if(info.dofsP[np]==0)
        {
            m_PhiA[np].setZero(rowPhi,info.dofsP[np]);
            m_Phi[np].setZero(rowPhi,info.dofsP[np]);
            return;
        }
        rhs.setZero(info.dofsB[np]+info.dofsI[np]+info.dofsP[np], info.dofsP[np]);
        rhs.bottomRows(info.dofsP[np]).setIdentity(info.dofsP[np],info.dofsP[np]);

        gsStopwatch time;
        solveKC<true,true>(np,rhs,solution);
        m_timings.KC_basis(np)+=time.stop();

        //TODO: Needs cleanup
        
        if(true|| m_IETIoptions.opt.getSwitch("SaddlePoint"))
        {
            m_PhiA[np] = solution.topRows(info.dofsB[np]+info.dofsI[np]); // have it in [normal] order
        }

        m_Phi[np] = solution.topRows(rowPhi);

        storeSpp(-solution.bottomRows(info.dofsP[np]),spp_loc);

        //---------------------------------------------------------
        /*
        // debug the BPCG constants:
        gsMatrix<T> out(info.dofsI[np]+info.dofsB[np],info.dofsB[np]+info.dofsI[np]);
        gsMatrix<T>col;
        for(int i=0; i<info.dofsB[np]+info.dofsI[np];++i)
        {
            m_inexact_KC[np]->getOperator(0,0)->apply(m_K[np].col(i),col);
            out.col(i)=col;
        }
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,Dynamic,Dynamic>  > eig(out);
        T max = eig.eigenvalues().maxCoeff();
        gsMatrix<T> PT = gsMatrix<T>::Identity(info.dofsI[np]+info.dofsB[np],info.dofsB[np]+info.dofsI[np])-m_exC[np].transpose()*m_PhiA[np].transpose();
        gsMatrix<T> P = gsMatrix<T>::Identity(info.dofsI[np]+info.dofsB[np],info.dofsB[np]+info.dofsI[np])-m_PhiA[np]*m_exC[np];
       // gsInfo<<std::fixed<<std::setprecision(14)<<std::endl;
        //gsInfo<<"out:\n"<<out<<"\n"<<std::flush;
        gsMatrix<T> Kd = (PT*m_K[np]*P);
     //   gsInfo<<Kd<<std::endl;
        for(int i=0; i<info.dofsB[np]+info.dofsI[np];++i)
        {
            m_inexact_KC[np]->getOperator(0,0)->apply(Kd.col(i),col);
            out.col(i)=col;
        }
       // gsInfo<<PT<<"\n\n"<<P<<"\n";
       // gsInfo<<"out:\n"<<out<<"\n"<<std::flush;
        Eigen::EigenSolver<Eigen::Matrix<real_t,Dynamic,Dynamic>  > eig2(out);
        gsMatrix<T> mins= eig2.eigenvalues().array().abs2().sqrt();
        T min=mins.maxCoeff();
        for(index_t i=0; i<mins.cols();++i)
            if(mins(i)<min && mins(i)>1.e-12)
                min = mins(i);

        gsInfo<<"np: "<<np<<"   Min: "<< min<<"  Max: "<<max<<" num ker: "<<(mins.array()<1.e-12).count()<<" alpha: "<<min/1.23<<"\n";
*/
        //---------------------------------------------------------

    }
    else
    {
        if(info.dofsP[np]==0)
            return;
        gsMatrix<T> temp;
        gsMatrix<T> KrpI = (*m_Krp[np])*gsMatrix<T>::Identity(info.dofsP[np],info.dofsP[np]);
        solveKrr(np,KrpI, temp);

        temp.transposeInPlace();
        gsMatrix<T> temp_2 = temp* (*m_Krp[np]);
        temp_2.transposeInPlace();
        storeSpp(m_Kpp[np] - temp_2,spp_loc);
    }

}

template <class T>
void gsIETIAssembler<T>::getPrimalInitialGuess(const gsMatrix<T> &allSol, gsMatrix<T> & sol, index_t c, int np) const
{
    sol.setZero();
    //  sol.topRows(info.dofsB[np]+info.dofsI[np])=m_exC[np].row(c).transpose();

    return;
    /*
    if( c==0 || c==info.dofsPtype[np][0] || c==info.dofsPtype[np][1]+info.dofsPtype[np][0] )
    {
        //sol.setRandom();
        sol.setZero();
        return;
    }
    if(c<info.dofsPtype[np][0])
    {
        boxCorner have(m_primalCorners[np][c-1]);
        boxCorner need(m_primalCorners[np][c]);

        gsVector<bool> have_b, need_b;
        have.parameters_into(m_patches[np].parDim(),have_b);
        need.parameters_into(m_patches[np].parDim(),need_b);
        sol = allSol.col(c-1);
        for(int d=0; d<m_patches[np].parDim();++d)
            if(have_b(d) ^ need_b(d) )
                sol=m_primalPerm[np][d]*sol;
    }
    else if(info.dim == 2 ||
            ( info.dim == 3 && c >= info.dofsPtype[np][1]+info.dofsPtype[np][0]))
    {
        index_t numb = info.dim ==2 ? info.dofsPtype[np][0] : info.dofsPtype[np][1]+info.dofsPtype[np][0];
        patchSide need =  m_patchAveragesSides[np][c-numb];
        bool found=false;
        for(index_t j=0; j<c-numb;++j)
        {
            const patchSide & have =  m_patchAveragesSides[np][j];
            if(have.direction() == need.direction())
            {
                sol=m_primalPerm[np][have.direction()]*allSol.col(j+numb);
                found=true;
                break;
            }

        }
        if(!found)
        {
            // sol.setRandom();
            sol.setZero();
            return;
        }

    }
    else
    {
        Edge need = m_primalEdges[np][c-info.dofsPtype[np][0]];

        bool found=false;
        for(int j=0; j<c-info.dofsPtype[np][0];++j)
        {
            const Edge & have = m_primalEdges[np][j];
            bool flag=true;
            for(int i=0; i<2;++i)
                if(have.directions[i] != need.directions[i] )
                    flag = false;
            if(flag)
            {
                sol= allSol.col(j+info.dofsPtype[np][0]);
                for(int i=0; i<2;++i)
                    if(have.parameters[i] ^ need.parameters[i])
                        sol=m_primalPerm[np][have.directions[i]]*sol;
                found = true;
                break;
            }
            else
                continue;
        }
        if(!found)
        {
            // sol.setRandom();
            sol.setZero();
            return;
        }

    }
*/
    gsSparseMatrix<T> CC = m_C[np]*(m_C[np].transpose());
    // gsInfo<<"CCT:\n"<<CC.toDense()<<"\n";
    //  gsMatrix<T> applyK;
    //  m_KC[np]->getOperator(0,0)->apply(sol.topRows(info.dofsB[np]+info.dofsI[np]), applyK);
    // gsInfo<<"rhs2: \n"<<-m_exC[np]*applyK<<"\n";
    //   gsInfo<<"rhs: \n"<<-m_exC[np]*(m_K[np]*sol.topRows(info.dofsB[np]+info.dofsI[np]))<<"\n";
    sol.bottomRows(info.dofsP[np])= CC.toDense().llt().solve(-m_exC[np]*(*m_K[np]*sol.topRows(info.dofsB[np]+info.dofsI[np])));
    // gsInfo<<"guess for Lagrange: \n"<<sol.bottomRows(info.dofsP[np])<<"\n";
}

template<class T>
void gsIETIAssembler<T>::storeSpp(const gsMatrix<T>& in, gsMatrix<T>& out)
{
    out = in;
}

template<class T>
void gsIETIAssembler<T>::assembleSpp(std::vector<gsMatrix<T> >& spp_loc)
{
    m_LU_Spp = NULL;
    if(info.dofTotalP==0)
        return;

    gsSparseMatrix<T> Spp(info.dofTotalP,info.dofTotalP);

    unsigned max = *max_element(info.dofsP.begin(),info.dofsP.end());
    unsigned min = *min_element(info.dofsP.begin(),info.dofsP.end());

    gsInfo<<"pDofs: min: "<<min<<", max: "<<max<<"\n";

    Spp.reservePerColumn(math::ipow(3,info.dim)*max);

    for(size_t np=0; np<info.numberPatches;np++)
    {
        for (index_t k=0; k<spp_loc[np].outerSize(); ++k)
            for (typename gsMatrix<T>::InnerIterator it(spp_loc[np],k); it; ++it)
                Spp.coeffRef(m_pDofsLoc2Glob[np][it.row()],m_pDofsLoc2Glob[np][it.col()])+= it.value();
    }

    //gsInfo<<"Spp: \n"<<Spp.toDense()<<std::endl<<std::endl;
    Spp.makeCompressed();
    if(m_IETIoptions.opt.getSwitch("SaddlePoint"))
        m_Spp = Spp;

    // Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real_t,Dynamic,Dynamic>  > eig(Spp.toDense());
    // gsInfo<<"Eigenvalues of Spp:\n"<<eig.eigenvalues()<<std::endl<<std::endl;

    m_LU_Spp=new sparseSPDfact(Spp);
    if(!m_LU_Spp->succeed())
        gsWarn<<"Sparse LLT of Spp not succeeded\n!";

}

template<class T>
void gsIETIAssembler<T>::assembleKrrKrpKpp(gsSparseMatrix<T> &matrix, size_t np)
{

    gsVector<index_t> splitting(2);
    splitting[0] = info.dofsP[np],splitting[1]= info.dofsR[np];
    typename gsSparseMatrix<T>::BlockView view = matrix.blockView(splitting, splitting);
    m_Kpp[np] = view(0,0);
    m_Krp[np] = gsSparseMatrix<T>(view(1,0)).moveToPtr();

    gsSparseMatrix<T> Krr = view(1,1);

    if(Krr.cols()==0)
        return;

    if(m_IETIoptions.opt.getSwitch("EnforceRescaling") )
        m_scalingsKrr[np].computeRef(Krr);

    switch (m_IETIoptions.KrrSolver)
    {
    case IETILocalSolver::Direct:
        m_LU_Krr[np]=new sparseSPDfact(Krr);
        break;
    default:
        GISMO_ERROR("Inexact solvers for Krr not implemented");
    }


}

template<class T>
void gsIETIAssembler<T>::assembleC(size_t np)
{
    gsMatrix<T> val(info.dofsB[np],1);

    m_C[np].resize(info.dofsP[np],info.dofsB[np]);
    int estimated_num=0;

    if(m_IETIoptions.strat.contains(primalDofMethod::vertices))
        estimated_num+=info.dofsP[np];
    if((info.dim == 2 && m_IETIoptions.strat.contains(primalDofMethod::edges)) || (info.dim == 3 && m_IETIoptions.strat.contains(primalDofMethod::faces)))
        estimated_num+= info.cDim*m_patchAveragesSides[np].size()* cast<T,int>( math::floor(math::sqrt(T(info.dofTotal))) );
    if(info.dim ==3 && m_IETIoptions.strat.contains(primalDofMethod::edges))
        estimated_num+= info.dofsPtype[np][primalDofMethod::getMethodIndex(primalDofMethod::edges)]*m_primalEdges[np].front().number;

    std::vector<Trip> tripletList;
    tripletList.reserve(estimated_num);

    int offset=0;
    if(m_IETIoptions.strat.contains(primalDofMethod::vertices))
    {
        for(size_t i=0; i<m_primalVdofs[np].size();i++)
        {
            int c = getComp(np,m_primalVdofs[np][i]);
            int idx = m_glob2BoundInteriorIndex[np][m_locDofsMapper[np][c].index(compCalcBack(np,m_primalVdofs[np][i]))];
            tripletList.push_back(Trip(i,idx,1));
        }
        offset+=m_primalVdofs[np].size();
    }


    if((info.dim == 2 && m_IETIoptions.strat.contains(primalDofMethod::edges)) || (info.dim == 3 && m_IETIoptions.strat.contains(primalDofMethod::faces)))
    {
        int method;
        info.dim == 2 ? method = primalDofMethod::getMethodIndex(primalDofMethod::edges) : method = primalDofMethod::getMethodIndex(primalDofMethod::faces);
        for(size_t i = 0; i< m_patchAveragesSides[np].size();i++)
        {
            patchSide side = m_patchAveragesSides[np][i];

            for(size_t c=0;c<info.cDim;c++)
            {
                val.setZero();
                getInterfaceAverageValue(side, val,c);

                for(int j = 0; j<val.rows();j++)
                    //  if(val(j,0)!=0 && m_stdMapper[c].is_free(m_boundDofs[np][j],np))
                    if(val(j,0)!=0 && m_locDofsMapper[np][c].is_free(m_boundDofs[np][j]))
                        tripletList.push_back(Trip(offset+info.cDim*i+c,j,val(j,0)));
            }
        }
        offset+=info.dofsPtype[np][method];
    }
    if(info.dim ==3 && m_IETIoptions.strat.contains(primalDofMethod::edges))
    {
        int method = m_IETIoptions.strat.contains(primalDofMethod::edges);

        for(size_t e=0;e<m_primalEdges[np].size();e++)
        {
            Edge edge = m_primalEdges[np][e];

            for(size_t c=0;c<info.cDim;c++)
            {
                typename gsBasis<T>::uPtr bbasis = m_basis[c].basis(np).boundaryBasis(edge.side3D);
                gsMatrix<index_t> bI = m_basis[c].basis(np).boundary(edge.side3D);
                gsMatrix<index_t> bbI = bbasis->boundary(edge.side2D);

                gsMatrix<T> valMatrix(bbI.rows(),1);
                valMatrix.setZero();

                getFaceEdgeAverageValue(edge,valMatrix,c,true);
                for(int j = 0;j<bbI.rows();j++)
                {
                    unsigned entr = (bI)((bbI)(j,0),0);
                    if(valMatrix(j,0)!=0 && m_locDofsMapper[np][c].is_free(entr))
                    {
                        int idx = m_glob2BoundInteriorIndex[np][m_locDofsMapper[np][c].index(entr)];
                        tripletList.push_back(Trip(offset+info.cDim*e+c,idx, valMatrix(j,0)));
                    }
                }
            }
        }
        offset+=info.dofsPtype[np][method];

    }
    m_C[np].setFromTriplets(tripletList.begin(), tripletList.end());
    //if(np==0)
    //    gsInfo<<m_C[np].toDense()<<"\n\n";
    GISMO_ASSERT(m_C[np].rows()<=m_C[np].cols(), "More primal varibles than interface dofs, cannot be linear independet. Refine your domain, or choose less primal dofs (different strategy)");
    //gsInfo<<"Patch: "<<np<<std::endl;
    //gsInfo<<m_C[np].toDense()<<std::endl<<std::endl;

}

template<class T>
void gsIETIAssembler<T>::assembleRhs(const std::vector<gsMatrix<T> >& rhs_loc, size_t np)
{

    if(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy") )
    {
        gsMatrix<T> rhs_b(info.dofsB[np], info.nRhs);
        m_rhs_d[np].setZero(info.dofsB[np],info.nRhs);
        m_rhs_i[np].setZero(info.dofsI[np],info.nRhs);

        unsigned idx;

        for(int k=0;k<rhs_loc[np].rows();k++)
        {
            idx = m_glob2BoundInteriorIndex[np][k];
            if(m_globIsBoundIndex[np][k] ==false)
                m_rhs_i[np].row(idx)=rhs_loc[np].row(k);
            else
                rhs_b.row(idx)=rhs_loc[np].row(k);
        }
        if(!m_IETIoptions.opt.getSwitch("SaddlePoint"))
        {
            gsMatrix<T> temp;
            solveKii<true>(np,m_rhs_i[np], temp);
           // temp.transposeInPlace();
          //  temp = temp*m_Kib[np];
          //  temp.transposeInPlace();
            temp = *m_Kbi[np]*temp;
            rhs_b-=temp; //fb - Kbi*Kii^-1*fi
        }
        else
        {
            rhs_b.conservativeResize(info.dofsB[np]+info.dofsI[np],rhs_b.cols());
            rhs_b.bottomRows(info.dofsI[np])= m_rhs_i[np];
            //   gsDebugVar(rhs_b);
        }

        embeddingTrans(rhs_b,np,m_rhs_p,m_rhs_d[np]);
        if(m_IETIoptions.opt.getSwitch("SaddlePoint"))
        {
            m_rhs_i[np]= m_rhs_d[np].bottomRows(info.dofsI[np]);
            m_rhs_d[np].conservativeResize(info.dofsB[np],m_rhs_d[np].cols());
        }
    }
    else
    {
        m_rhs_d[np].setZero(info.dofsR[np],info.nRhs);

        gsMatrix<T> rhs_p;
        rhs_p.setZero(info.dofsP[np],info.nRhs);

        for(size_t i=0;i<m_remDofs[np].size();i++)
        {
            int c = getComp(np,m_remDofs[np][i]);
            m_rhs_d[np].row(i) = rhs_loc[np].row(m_locDofsMapper[np][c].index(compCalcBack(np,m_remDofs[np][i])));
        }

        for(size_t k=0;k<m_primalVdofs[np].size();k++)
        {
            int c = getComp(np,m_primalVdofs[np][k]);
            rhs_p.row(k)= rhs_loc[np].row( m_locDofsMapper[np][c].index(compCalcBack(np,m_primalVdofs[np][k])));
        }

        assemblePrimal(rhs_p,np,m_rhs_p);
    }

}


template<class T>
void gsIETIAssembler<T>::assembleRhsFreeElim()
{
    //assemble the rhs for the coupled eliminated and free dofs (as lagrange mult.)
    if(!nCoupledElimDofs())
    {
        int sign;
        m_rhs_dir.setZero(m_lagrangeTable.size(),info.nRhs);
        for(std::vector<index_t>::iterator it = m_freeElimLagr.begin(); it!=m_freeElimLagr.end();it++)
        {
            patchDof p1= m_lagrangeTable[*it].first;
            patchDof p2 = m_lagrangeTable[*it].second;
            p1.first > p2.first ? sign = 1 : sign = -1;

            int c = getComp(p2.first,p2.second); // p1 and p2 must have the same component
            if(m_locDofsMapper[p1.first][c].is_free(compCalcBack(p1.first,p1.second)))
                m_rhs_dir.row(*it)=sign*m_dirDofs[p2.first].row(m_locDofsMapper[p2.first][c].bindex(compCalcBack(p2.first,p2.second))); //<<< continue HERE.
            else if(m_locDofsMapper[p2.first][c].is_free(compCalcBack(p2.first,p2.second)))
                m_rhs_dir.row(*it)=-sign*m_dirDofs[p1.first].row(m_locDofsMapper[p1.first][c].bindex(compCalcBack(p1.first,p1.second)));
            else
            {
                GISMO_ERROR("BOTH DOFS ARE DIRICHLET DOFS; THIS IS NOT ALLOWED");
            }
        }
    }
}


template<class T>
void gsIETIAssembler<T>:: getInterfaceAverageValue(const patchSide& side, gsMatrix<T>& average,int c)
{
    int patch = side.patch;
    //The boundary basis functions
    average.setZero(info.dofsB[patch],1);
    const gsBasis<T> & basis = m_basis[c].basis(patch);

    //the basis entries of the boundary of the boundary
    //const gsMatrix<index_t> bbI = basis.boundary( );
    //const gsMatrix<index_t> bI = basis.boundary( );

    typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(side);
    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_VALUE|NEED_JACOBIAN, m_patches[patch]));

    // Temporaries
    gsMatrix<T> quNodes, basisVals;
    gsVector<T> quWeights;

    gsMatrix<index_t> actives;
    index_t numActive;

    //get quadrature degree
    gsVector<index_t> numQuNodes( m_basis[c].dim());
    for( int i=0; i < m_basis[c].dim(); i++)
        numQuNodes[i] = 2*(m_basis[c].basis(patch).degree(i)+1);
    numQuNodes[side.direction()] = 1;

    gsGaussRule<T> QuRule(numQuNodes);

    // T measure = 0;
    gsVector<T> normal;

    // Create the iterator along the given part boundary.
    for(; bdryIter->good(); bdryIter->next() )
    {
        QuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                      quNodes, quWeights);

        // Compute the active basis functions
        basis.active_into(quNodes.col(0) , actives);
        numActive = actives.rows();

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt( quNodes);

        //Evaluate the basis
        basis.eval_into(quNodes, basisVals);

        // Do the actual assembly:
        for( index_t i=0; i < numActive; ++i)
        {
            int glIdx = m_locDofsMapper[patch][c].index(actives(i));
            if(m_stdMapper[c].is_free(actives(i),patch) && m_globIsBoundIndex[patch][glIdx] )
            {
                const int ii = m_glob2BoundInteriorIndex[patch][m_locDofsMapper[patch][c].index(actives(i))];
                for(index_t k=0; k<quWeights.rows();k++)
                {
                    geoEval->outerNormal(k,side,normal);
                    const T weight_k = quWeights[k]* normal.norm();

                    average(ii,0) += weight_k*basisVals(i,k);
                }
            }
        }

        /*
        //calculate the measure of the face
        for(index_t k=0; k<quWeights.rows();k++)
        {
            geoEval->outerNormal(k,side,normal);
            measure += quWeights[k] * normal.norm();
        }
        */

    }



    //divide by the measure to get average;
    T sum = average.sum();
    GISMO_ASSERT(sum!=0, "Zero row in contraint matrix C" );
    average /=sum;
}

template<class T>
void gsIETIAssembler<T>:: getFaceEdgeAverageValue(const Edge& edge, gsMatrix<T>& average, int c, bool eliminateCorners)
{
    int patch = edge.side3D.patch;

    //The boundary basis functions
    const typename gsBasis<T>::uPtr basis = m_basis[c].basis(patch).boundaryBasis(edge.side3D);
    const typename gsBasis<T>::uPtr bbasis = basis->boundaryBasis(edge.side2D);

    //the basis entries of the boundary of the boundary
    const gsMatrix<index_t> bbbI = bbasis->allBoundary( );

    typename gsBasis<T>::domainIter bdryIter = m_basis[c].basis(patch).makeDomainIterator(edge.side3D);
    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_JACOBIAN, m_patches[patch]));

    // Temporaries
    gsMatrix<T> quNodes, quNodesd_,basisVals;
    gsVector<T> quWeights;

    gsMatrix<index_t> actives;
    index_t numActive;

    //get quadrature degree
    gsVector<index_t> numQuNodes(m_basis[c].dim());
    for( int i=0; i < m_basis[c].dim(); i++)
        numQuNodes[i] = 2*(m_basis[c].basis(patch).degree(i)+2);
    numQuNodes[edge.directions[0]] = 1;
    numQuNodes[edge.directions[1]] = 1;

    gsGaussRule<T> QuRule(numQuNodes);
    int thirdDir=0;
    for(int i=0; i<3;i++)
        if(edge.directions[0]!=i && edge.directions[1]!=i)
        {
            thirdDir = i;
            break;
        }
    gsMatrix<T> jacobian;

    T measure =0;
    // Create the iterator along the given part boundary.
    for(; bdryIter->good(); bdryIter->next() )
    {
        gsMatrix<T> lC = bdryIter->lowerCorner();
        gsMatrix<T> uC = bdryIter->upperCorner();

        //Extract the edge cubes
        if(edge.side2D.parameter() == false)
        {
            if(math::min(lC(edge.directions[1]),uC(edge.directions[1]))!=0)
                continue;
            else
                uC(edge.directions[1])=0;
        }
        else
        {
            if(math::max(lC(edge.directions[1]),uC(edge.directions[1]))!=1)
                continue;
            else
                lC(edge.directions[1])=1;
        }


        QuRule.mapTo( lC, uC, quNodes, quWeights);


        //reduce the dimension of the quadrature points to fit the dimension of the boundary
        int numRows = quNodes.rows()-1;
        int numCols = quNodes.cols();
        quNodesd_ = quNodes;

        int start =0;
        if(edge.directions[0]<edge.directions[1])
            start =1;
        for(int i=start;i<start+2;i++)
        {
            int delRow = edge.directions[i%2];
            if(delRow<numRows)
                quNodesd_.block(delRow,0,numRows-delRow,numCols) = quNodesd_.block(delRow+1,0,numRows-delRow,numCols);
            quNodesd_.conservativeResize(numRows,numCols);

            numRows--;
        }


        // Compute the active basis functions
        bbasis->active_into(quNodesd_.col(0), actives);
        numActive = actives.rows();

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt( quNodes );

        //Evaluate the basis
        bbasis->eval_into(quNodesd_, basisVals);

        // Do the actual assembly:
        for( index_t i=0; i < numActive; ++i)
        {
            const int ii = actives(i);

            for(index_t k=0; k<quWeights.rows();k++)
            {
                jacobian= geoEval->jacobian(k);
                const T weight_k = quWeights[k]* jacobian.col(thirdDir).norm();

                average(ii,0) += weight_k*basisVals(i,k);
            }
        }

        //calculate the measure of the edge
        for(index_t k=0; k<quWeights.rows();k++)
        {
            jacobian= geoEval->jacobian(k);
            measure += quWeights[k] * jacobian.col(thirdDir).norm();
        }

    }

    //ignore the boundary basis functions
    if(eliminateCorners)
        for(int i=0;i<bbbI.rows();i++)
            average.row((bbbI)(i,0)).setZero();

    //divide by the measure to get average;
    average /=measure;
}

template<class T>
void gsIETIAssembler<T>::checkPrimalDofStrategy(bool isMinimalEnergy,int dim)
{
    if(isMinimalEnergy==false && m_IETIoptions.strat.strat!= primalDofMethod::A)
    {
        printWarn("Your primal dof strategy is not supported by the method, choose isMinimalEnergy as true to use your strategy. Strategy is changed to A.");
        m_IETIoptions.strat.strat=primalDofMethod::A;
    }
    if(dim == 3)
    {
        if(m_IETIoptions.strat.strat== primalDofMethod::A)
        {
            printWarn("Only vertex primal dofs in 3D lead to not optimal iteration number. In order to speed up the algorithm, use a different stretegy, e.g. B or C");
        }
    }
    if(dim ==2)
    {
        if(m_IETIoptions.strat.contains(primalDofMethod::faces))
        {
            printWarn("Your primal dof strategy is not supported in this dimension. For a 2D object, only vertices and edges are available! Strategy is changed to C");
            m_IETIoptions.strat.strat = primalDofMethod::C;
        }
    }
}

template<class T>
void gsIETIAssembler<T>::processSolution(const gsMatrix<T>& uP, const std::vector<gsMatrix<T> >& u2, gsMatrix<T>& solVec) const
{
    unsigned totSize = 0;
    for(size_t c = 0;c<info.cDim;c++)
        totSize+= m_stdMapper[c].freeSize();

    solVec.resize(totSize,info.nRhs);

    std::vector<gsDofMapper > stdMapper(info.cDim);
    for(size_t c=0; c<info.cDim;c++)
        stdMapper[c]= m_assembler->system().colMapper(c);//m_stdMapper[c]; //FIXME: m_assembler->dofMapper(c);


    if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
    {
#pragma omp parallel
        {
            gsMatrix<T> xPi;
#pragma omp for schedule(static, 1) nowait
            for(size_t np=0;np<info.numberPatches;np++)
            {
                distributePrimal(uP,np,xPi);

                for(size_t c=0;c<info.cDim;c++)
                {
                    size_t sz  = m_basis[c][np].size();
                    for (size_t i = 0; i < sz; ++i)
                    {
                        int idx = m_locDofsMapper[np][c].index(i);
                        int globI = stdMapper[c].index(i,np);

                        if(stdMapper[c].is_free(i,np))
                        {
                            if(m_globIsPrimalIndex[np][idx])
                                solVec.row(globI) = xPi.row(m_glob2PrimalRemainingIndex[np][idx]);
                            else
                                solVec.row(globI) = u2[np].row(m_glob2PrimalRemainingIndex[np][idx]);
                        }
                    }
                }
            }

        }
    }
    else
    {
#pragma omp parallel
        {
            gsMatrix<T> xI, temp, w;

#pragma omp for schedule(static, 1) nowait
            for(size_t np=0; np<info.numberPatches;np++)
            {
                embedding(uP,u2[np],np,w);

                temp.noalias() = m_rhs_i[np]-*m_Kib[np]*w;
                solveKii<true>(np,temp, xI) ;

                for(size_t c=0;c<info.cDim;c++)
                {
                    size_t sz  =  m_basis[c][np].size();
                    for (size_t i = 0;i < sz ; ++i)
                    {
                        int idx = m_locDofsMapper[np][c].index(i);
                        int globI = stdMapper[c].index(i,np);
                        if(stdMapper[c].is_free(i,np))
                        {
                            if(m_globIsBoundIndex[np][idx])
                                solVec.row(globI) = w.row(m_glob2BoundInteriorIndex[np][idx]);
                            else
                                solVec.row(globI)  = xI.row(m_glob2BoundInteriorIndex[np][idx]);
                        }
                    }
                }
            }
        }
    }


}

template<class T>
void gsIETIAssembler<T>::processSolution(typename gsMatrix<T>::BlockView& view, gsMatrix<T>& solVec) const
{
    unsigned totSize = 0;
    for(size_t c = 0;c<info.cDim;c++)
        totSize+= m_stdMapper[c].freeSize();

    solVec.resize(totSize,info.nRhs);

    std::vector<gsDofMapper > stdMapper(info.cDim);
    for(size_t c=0; c<info.cDim;c++)
        stdMapper[c]= m_assembler->system().colMapper(c);//m_stdMapper[c]; //FIXME: m_assembler->dofMapper(c);


    if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
    {
#pragma omp parallel
        {
            gsMatrix<T> xPi;
#pragma omp for schedule(static, 1) nowait
            for(size_t np=0;np<info.numberPatches;np++)
            {
                distributePrimal(view(info.numberPatches),np,xPi);

                for(size_t c=0;c<info.cDim;c++)
                {
                    size_t sz  = m_basis[c][np].size();
                    for (size_t i = 0; i < sz; ++i)
                    {
                        int idx = m_locDofsMapper[np][c].index(i);
                        int globI = stdMapper[c].index(i,np);

                        if(stdMapper[c].is_free(i,np))
                        {
                            if(m_globIsPrimalIndex[np][idx])
                                solVec.row(globI) = xPi.row(m_glob2PrimalRemainingIndex[np][idx]);
                            else
                                solVec.row(globI) = view(np).row(m_glob2PrimalRemainingIndex[np][idx]);
                        }
                    }
                }
            }

        }
    }
    else
    {
#pragma omp parallel
        {
            gsMatrix<T> w;

#pragma omp for schedule(static, 1) nowait
            for(size_t np=0; np<info.numberPatches;np++)
            {
                embedding(view(info.numberPatches),view(np),np,w);

                for(size_t c=0;c<info.cDim;c++)
                {
                    size_t sz  =  m_basis[c][np].size();
                    for (size_t i = 0;i < sz ; ++i)
                    {
                        int idx = m_locDofsMapper[np][c].index(i);
                        int globI = stdMapper[c].index(i,np);
                        if(stdMapper[c].is_free(i,np))
                        {
                            if(m_globIsBoundIndex[np][idx])
                                solVec.row(globI) = w.row(m_glob2BoundInteriorIndex[np][idx]);
                            else
                                solVec.row(globI) = w.row(info.dofsB[np]+m_glob2BoundInteriorIndex[np][idx]);
                        }
                    }
                }
            }
        }
    }

}


//w has size info.dofsB[np]
template<class T>
void gsIETIAssembler<T>::embedding(gsMatrix<T> const & xP,std::vector<gsMatrix<T> >const & x2,std::vector<gsMatrix<T> >& w) const
{
    for(size_t np=0; np<info.numberPatches;np++)
        embedding(xP,x2[np],np,w[np]);
}

template<class T>
void gsIETIAssembler<T>::embedding(gsMatrix<T> const & xP, gsMatrix<T> const & x2, size_t np, gsMatrix<T>& w) const
{
    gsStopwatch time;
    int nRhs = x2.cols();
    gsMatrix<T> xPi;

    distributePrimal(xP,np,xPi);
    if(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
    {
        w=x2;
        if(info.dofsP[np]!=0)
            w.noalias()+=m_Phi[np]*xPi;
    }
    else
    {
        unsigned globI;
        w.setZero(info.dofsB[np],nRhs);
        for (index_t i=0;i<info.dofsB[np]; ++i)
        {
            globI = m_locDofsMapper[np][getComp(np,m_boundDofs[np][i])].index(compCalcBack(np,m_boundDofs[np][i]));

            if(m_globIsPrimalIndex[np][globI])
                w.row(i) = xPi.row(m_glob2PrimalRemainingIndex[np][globI]);
            else
                w.row(i) = x2.row(m_glob2PrimalRemainingIndex[np][globI]);


        }

    }
    m_timings.embedding[np]+=time.stop();
}

//w has size info.dofsB[np]
template<class T>
void gsIETIAssembler<T>::embeddingTrans(const std::vector<gsMatrix<T> >& w,gsMatrix<T>&  uP, std::vector<gsMatrix<T> > &  u2) const
{
    int nRhs = w.front().cols();
    uP.setZero(info.dofTotalP,nRhs);

    for(size_t np=0; np<info.numberPatches;np++)
        embeddingTrans(w[np],np,uP,u2[np]);
}

//w has size info.dofsB[np]
template<class T>
void gsIETIAssembler<T>::embeddingTrans(const gsMatrix<T> & w,size_t np, gsMatrix<T>&  uP, gsMatrix<T> &  u2) const
{
    gsStopwatch time;
    gsMatrix<T> xPi;
    int nRhs = w.cols();
    if(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
    {
        u2=w;
        gsMatrix<T> temp_trans = (m_Phi[np].topRows(u2.rows())).transpose()*u2;
        u2.topRows(info.dofsB[np]) -= (m_C[np].transpose())*temp_trans;

        assemblePrimal(temp_trans,np,uP);

    }
    else
    {
        xPi.setZero(info.dofsP[np],nRhs);

        u2.setZero(info.dofsR[np],nRhs);

        for(size_t k=0;k<m_primalVdofs[np].size();k++)
            xPi.row(k) += w.row( m_glob2BoundInteriorIndex[np][m_locDofsMapper[np][getComp(np,m_primalVdofs[np][k])].index(compCalcBack(np,m_primalVdofs[np][k]))]);

        assemblePrimal(xPi,np,uP);

        for(size_t i=0;i<m_remDofs[np].size();i++)
            if(m_remDofIsAlsoBoundDof[np][i])
                u2.row(i) = w.row(m_glob2BoundInteriorIndex[np][m_locDofsMapper[np][getComp(np,m_remDofs[np][i])].index(compCalcBack(np,m_remDofs[np][i]))]);

    }
    m_timings.embeddingT[np]+=time.stop();
}

template< typename T>
void gsIETIAssembler<T>::applyDP(unsigned np,const gsMatrix<T>& xP , gsMatrix<T>& uD) const
{
    gsMatrix<T> wi, xPi;
    distributePrimal(xP,np,xPi);
    gsMatrix<T> temp = m_Phi[np]*xPi;
    gsMatrix<T> ev = *m_Kib[np]*temp;
    solveKii<true>(np,ev,wi);
    temp = *m_Kbb[np]*temp-*m_Kbi[np]*wi;
    gsMatrix<T> u1 = m_C[np].transpose()*m_Phi[np].transpose()*temp;
    uD = temp - u1;
}

template< typename T>
void gsIETIAssembler<T>::applyDP(const gsMatrix<T>& xP , std::vector<gsMatrix<T> >& uD) const
{
    for(size_t np=0; np<info.numberPatches;np++)
    {
        applyDP(np,xP,uD[np]);
    }
}

template< typename T>
void gsIETIAssembler<T>::applyPD(unsigned np,const gsMatrix<T>& uD, gsMatrix<T>& xP) const
{
    gsMatrix<T> wi;

    gsMatrix<T> ev = *m_Kib[np]*uD;
    solveKii<true>(np,ev,wi);
    ev = m_Phi[np].transpose()*(*m_Kbb[np]*uD-*m_Kbi[np]*wi);
    assemblePrimal(ev,np,xP);
}

template< typename T>
void gsIETIAssembler<T>::applyPD(const std::vector<gsMatrix<T> >& uD, gsMatrix<T>& xP) const
{
    xP.setZero();
    for(size_t np=0; np<info.numberPatches;np++)
    {
        applyPD(np,uD[np],xP);
    }
}

template< typename T>
template<bool exact, bool basis>
void gsIETIAssembler<T>::solveKC(unsigned np, gsMatrix<T>& rhs, gsMatrix<T>& sol) const
{
    GISMO_ASSERT(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy"), "Not available for this algorithm, choose NoMinimumEnergy as false");
    sol.setZero(rhs.rows(),rhs.cols());

    if(rhs.rows()==0)
        return;

    if(m_IETIoptions.opt.getSwitch("EnforceRescaling"))
        for(int col=0;col<rhs.cols();col++)
            rhs.col(col) = m_scalingsKC[np].LeftScaling().cwiseProduct(rhs.col(col));
    gsStopwatch time;
    bool isBI=true;
    bool needsProjection = false;

    switch (m_IETIoptions.KCSolver)
    {
    case IETILocalSolver::Direct:
        sol=m_LU_KC[np]->solve(rhs);
        break;
    case IETILocalSolver::FastDiagonalization:
    {
        rhs.topRows(info.dofsB[np]+info.dofsI[np])=m_permInvMat[np]*rhs.topRows(info.dofsB[np]+info.dofsI[np]);
        isBI=false;
        if(basis)
        {
            for(int r=0; r<rhs.cols();r++)
            {
                const gsMatrix<T>& rhsi=rhs.col(r);
                gsMatrix<T> soli=sol.col(r);
                getPrimalInitialGuess(sol,soli,r,np);

                m_KCBasis[np]->solve(rhsi,soli);

                //minres.solve(rhsi,soli);
#pragma omp critical
                {
                    m_timings.averageItNumKCBasis[np]+=m_KCBasis[np]->iterations();
                    m_timings.numKCItBasis[np]++;
                }
                gsInfo<<"#rhs: "<<r<<", KC - BPCG-BS iterations: "<<m_KCBasis[np]->iterations()<<" - BS residual: "<<m_KCBasis[np]->error()<<"\n";
                sol.col(r) = soli;

            }
        }
        else
        {
            needsProjection=true;
            for(int r=0; r<rhs.cols();r++)
            {
                const gsMatrix<T>& rhsi=rhs.col(r);
                gsMatrix<T> soli(info.dofsB[np]+info.dofsI[np],1);
                soli.setZero();

                if(exact)
                {
                    m_KCSolve[np]->solve(rhsi.topRows(info.dofsB[np]+info.dofsI[np]),soli);
#pragma omp critical
                    {
                        m_timings.averageItNumKCSpace[np]+=m_KCSolve[np]->iterations();
                        m_timings.numKCItSpace[np]++;
                    }
                    gsInfo<<"#rhs: "<<r<<", KC - CG-BS iterations: "<<m_KCSolve[np]->iterations()<<" - BS residual: "<<m_KCSolve[np]->error()<<"\n";
                }
                else
                {
                    gsMatrix<T> rhs_i=rhsi.topRows(info.dofsB[np]+info.dofsI[np]);
                    int max_iter = m_IETIoptions.opt.getInt("MG_KC.Solve.MaxIterations");
                    gsPreconditionerOp<T>* stOp =  dynamic_cast<gsPreconditionerOp<T>*>((m_inexact_KC[np]->getOperator(0,0).get()));
                    if(stOp!=NULL)
                    {
                        for( int i=0; i<max_iter;++i )
                            stOp->step(rhs_i, soli);
                    }
                    else
                    {
		      /*
                        gsConjugateGradient<> cg(*m_K[np],m_inexact_KC[np]->getOperator(0,0));
                        cg.setMaxIterations(max_iter);
                        cg.solve(rhs_i,soli);
*/
                        
                        m_inexact_KC[np]->getOperator(0,0)->apply(rhs_i,soli);
                        gsMatrix<T> soli_old;
                        gsMatrix<T> temp;
                        for( int i=1; i<max_iter;++i )
                        {
                            soli_old = soli;
                            temp = rhs_i - *m_K[np]*soli_old;
                            m_inexact_KC[np]->getOperator(0,0)->apply(temp,soli_old);
                            soli += (real_t)1/m_maxEigKC[np]*soli_old;
                        }
                        
			//gsInfo<<"KC - CG-BS iterations: "<<cg.iterations()<<" - BS residual: "<<cg.error()<<"\n";
 
                    }
                }
                sol.block(0,r,info.dofsB[np]+info.dofsI[np],1) = soli;
            }
        }

        break;
    }
    case IETILocalSolver::Multigrid:
    {
        //[B,I] -> [normal]
        rhs.topRows(info.dofsB[np]+info.dofsI[np])=m_permInvMat[np]*rhs.topRows(info.dofsB[np]+info.dofsI[np]);
        isBI=false;
        /*
        gsMatrix<T> dummy(info.dofTotalP,rhs.cols());
        gsMatrix<T> out(info.dofsP[np]+info.dofsB[np]+info.dofsI[np],info.dofsP[np]+info.dofsB[np]+info.dofsI[np]);
        for(int i=0; i<info.dofsP[np]+info.dofsB[np]+info.dofsI[np];++i)
        {
            gsMatrix<T> col;
           m_inexact_KC[np]->apply(gsMatrix<T>::Identity(info.dofsP[np]+info.dofsB[np]+info.dofsI[np],info.dofsP[np]+info.dofsB[np]+info.dofsI[np]).col(i),col);
            out.col(i) = col;
        }
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real_t,Dynamic,Dynamic>  > eig(out.block(0,0,info.dofsB[np]+info.dofsI[np],info.dofsB[np]+info.dofsI[np]));
        gsInfo<<out<<"\n\n";
        gsInfo<<"Eigenvalues of MG_S:\n"<<eig.eigenvalues()<<std::endl<<std::endl;
        */

      
        if(exact&&!basis && m_bc_loc[np].dirichletSides().size()==0)
        {
            orthProjToDualofDualsubspace proj(m_PhiA[np],m_exC[np],info,np);
            gsMatrix<T> temp;
            proj.apply(rhs.topRows(info.dofsB[np]+info.dofsI[np]),temp);
            rhs.topRows(info.dofsB[np]+info.dofsI[np])=temp;
        }
        
        if(basis)
        {
            /*
            gsMinimalResidual<> minres(m_KC[np],m_inexact_KC[np]);
            minres.setInexactResidual(false);
//*/
            /*
            gsGMRes<> minres(m_KC[np],nonSymBlockPreconditioner::make(m_inexact_KC[np]->getOperator(0,0),m_inexact_KC[np]->getOperator(1,1),makeMatrixOp(m_exC[np])));
//*/
            // gsVector<T,2> scaling;
            //  gsBramblePasciakCG<gsBPCG_Types::SchoeberlZulehner> minres(m_KC[np],m_inexact_KC[np]->getOperator(0,0),m_inexact_KC[np]->getOperator(1,1));
            //*/


            for(index_t r=0; r<rhs.cols();r++)
            {
                const gsMatrix<T>& rhsi=rhs.col(r);
                gsMatrix<T> soli=sol.col(r);
                getPrimalInitialGuess(sol,soli,r,np);

                gsMatrix<T> hist;
                m_KCBasis[np]->solveDetailed(rhsi,soli,hist);
                //m_KCBasis[np]->solve(rhsi,soli);

                //minres.solve(rhsi,soli);
#pragma omp critical
                {
                    m_timings.averageItNumKCBasis[np]+=m_KCBasis[np]->iterations();
                    m_timings.numKCItBasis[np]++;
                }
                gsInfo<<"#rhs: "<<r<<", KC - BPCG-MG iterations: "<<m_KCBasis[np]->iterations()<<" - MG residual: "<<m_KCBasis[np]->error()<<"\n";
		if(m_KCBasis[np]->iterations() > 100)
		  gsInfo<<hist<<"\n\n";
                sol.col(r) = soli;

            }
        }
        else
        {
            needsProjection=true;
            for(index_t r=0; r<rhs.cols();r++)
            {
                const gsMatrix<T>& rhsi=rhs.col(r);
                gsMatrix<T> soli(info.dofsB[np]+info.dofsI[np],1); //=sol.block(0,r,info.dofsB[np]+info.dofsI[np],1);
                soli.setZero();
                // soli.conservativeResize(info.dofsB[np]+info.dofsI[np],soli.cols());

                if(exact)
                {
                    // soli.setZero();
                    //soli.setRandom();

                    /*
                    std::vector<gsLinearOperator<>::Ptr > ops, ops2;

                    if(m_bc_loc[np].dirichletSides().size()!=0)
                        ops.push_back(makeMatrixOp(m_K[np])) ;
                    else
                    {
                        //     ops.push_back(permutation::make(m_permMat[np]));
                          //   ops.push_back(orthProjToDualsubspace::make(m_PhiA[np],m_exC[np],info,np));
                        //     ops.push_back(permutation::make(m_permInvMat[np]));
                        ops.push_back(makeMatrixOp(m_K[np])) ;
                        //     ops.push_back(permutation::make(m_permMat[np]));
                        //     ops.push_back(orthProjToDualofDualsubspace::make(m_PhiA[np],m_C[np],info,np));
                        //      ops.push_back(permutation::make(m_permInvMat[np]));
                    }
                    gsProductOp<>::Ptr conOp1 = gsProductOp<>::make(ops);

                    if(m_bc_loc[np].dirichletSides().size()!=0)
                        ops2.push_back(m_inexact_KC[np]->getOperator(0,0)) ;
                    else
                    {
                        //ops2.push_back(permutation::make(m_permMat[np]));
                        //ops2.push_back(orthProjToDualofDualsubspace::make(m_Phi[np],m_C[np],info,np));
                        //ops2.push_back(permutation::make(m_permInvMat[np]));
                        ops2.push_back(m_inexact_KC[np]->getOperator(0,0)) ;
                        // ops2.push_back(permutation::make(m_permMat[np]));
                       // ops2.push_back(orthProjToDualsubspace::make(m_PhiA[np],m_exC[np],info,np)); //not working... idk why
                        //   ops2.push_back(permutation::make(m_permInvMat[np]));
                    }
                    gsProductOp<>::Ptr conOp2 = gsProductOp<>::make(ops2);
*/
                    /*
                    gsMatrix<T> out(info.dofsB[np]+info.dofsI[np],info.dofsB[np]+info.dofsI[np]);
                    for(int i=0; i<info.dofsB[np]+info.dofsI[np];++i)
                    {
                        gsMatrix<T> col;
                        conOp2->apply(gsMatrix<T>::Identity(info.dofsB[np]+info.dofsI[np],info.dofsB[np]+info.dofsI[np]).col(i),col);
                        out.col(i) = col;
                    }
*/
                    //  gsInfo<<out<<"\n\n";

                    //   gsMatrix<T> mat;
                    //   conOp1.toMatrix(mat);
                    //   gsInfo<<mat<<"\n\n";

                    // unsigned f = 50;

                    /*
                        orthProjToDualofDualsubspace proj(m_Phi[np],m_C[np],info,np);
                         proj.apply(rhsi.topRows(info.dofsB[np]+info.dofsI[np]),rs);
                      gsInfo<<"Norm: "<<(rs-rhsi.topRows(info.dofsB[np]+info.dofsI[np])).norm()<<"\n";

                        cg.solve(rhsi.topRows(info.dofsB[np]+info.dofsI[np]),soli,conOp2);
                     cg.solve(rhsi.topRows(info.dofsB[np]+info.dofsI[np]),soli,*gsIdentityOp::make(m_inexact_KC[np]->getOperator(0,0)->rows()));
                    gsInfo<<"rhs: "<<rhsi.topRows(info.dofsB[np]+info.dofsI[np]).transpose()<<"\n";
                    cg.solve(rhsi.topRows(info.dofsB[np]+info.dofsI[np]),soli,*(m_inexact_KC[np]->getOperator(0,0)));

                       cg.setPreconditioner(conOp2);
                      gsMatrix<T> hist;

                    gsInfo<<std::setprecision(16)<<std::fixed;
                    gsInfo<<rhsi.topRows(info.dofsB[np]+info.dofsI[np])<<"\n\n";
                    gsInfo<<m_PhiA[np]<<"\n\n";
                    gsInfo<<m_exC[np].toDense()<<"\n\n";
                    gsInfo<<m_K[np].toDense()<<"\n\n";

                     cg.solveDetailed(rhsi.topRows(info.dofsB[np]+info.dofsI[np]),soli,hist);
                     gsInfo<<"#rhs: "<<r<<", KC - CG-MG iterations: "<<cg.iterations()<<" - MG residual: "<<cg.error()<<"\n"<<hist<<"\n";
                    unsigned k=0;
                    */
		    real_t res;
		    gsMatrix<T> hist;
		    
		   // if(m_bc_loc[np].dirichletSides().size()==0)
		    //  m_KCSolve[np]->setMaxIterations(12);
		    
                    m_KCSolve[np]->solveDetailed(rhsi.topRows(info.dofsB[np]+info.dofsI[np]),soli,hist);
		   
		    res = m_KCSolve[np]->error();
		    index_t iter = m_KCSolve[np]->iterations();
#pragma omp critical
                    {
                        m_timings.averageItNumKCSpace[np]+=m_KCSolve[np]->iterations();
                        m_timings.numKCItSpace[np]++;
                    }
                   // int count = 0;
                    /*
		    real_t rhs_norm;
                    f(m_KCSolve[np]->iterations() >= 11)
		      gsInfo<<hist<<"\n\n";
                    gsMatrix<T> resVec;
                    if(res > m_KCSolve[np]->tolerance())
		    {
		      soli = m_permMat[np]*soli;
		      projectToDualSubspaceA(soli,np);
		      soli = m_permInvMat[np]*soli;
		      resVec = rhsi.topRows(info.dofsB[np]+info.dofsI[np]) - (*m_K[np])*soli;
		      //resVec = m_permMat[np]*resVec;
		     // projectToDualOfDualSubpace(resVec,np);
		     // resVec = m_permInvMat[np]*resVec;
		       rhs_norm = m_KCSolve[np]->rhsNorm();
		      gsInfo<<"Old residual is: "<<resVec.norm()/rhs_norm<<"\n";
		      
		    }
		    
                    while(res > m_KCSolve[np]->tolerance() && iter <200)
                    {
		      gsMatrix<T> x_delta;
		      x_delta.setZero(info.dofsB[np]+info.dofsI[np],1);
		      m_KCSolve[np]->solveDetailed(resVec,x_delta,hist);
		      x_delta = m_permMat[np]*x_delta;
		      projectToDualSubspaceA(x_delta,np);
		      x_delta  = m_permInvMat[np]*x_delta;
		      soli = soli+ x_delta;
		      resVec = rhsi.topRows(info.dofsB[np]+info.dofsI[np]) - (*m_K[np])*soli;
		       //resVec = m_permMat[np]*resVec;
		     // projectToDualOfDualSubpace(resVec,np);
		     // resVec = m_permInvMat[np]*resVec;
		      res = resVec.norm()/rhs_norm;
		      gsInfo<<"New residual is: "<<res<<"\n";
		      count ++;
		      iter += m_KCSolve[np]->iterations();
		      #pragma omp critical
                    {
                        m_timings.averageItNumKCSpace[np]+=m_KCSolve[np]->iterations();
                        m_timings.numKCItSpace[np]++;
                    }
                    }
                   */ 
                    /*
                    for(k=0; k<(unsigned)max_iter/f;k++)
                    {
                      //  if((k!=0 && cg.error()<tol*1.e2))

                      //  else
                       //     cg.setPreconditioner(m_inexact_KC[np]->getOperator(0,0));
                        gsMatrix<T> hist;
                        //cg.solveDetailed(rhsi.topRows(info.dofsB[np]+info.dofsI[np]),soli,hist);

                        //gsInfo<<"error in "<<k<<" : "<<cg.error()<<"\n"<<hist<<"\n";


                        gsMatrix<T> temp, K2;
                        conOp1->apply(soli,temp);
                        gsInfo<<"Ku:\n"<<temp<<"\n\n";
                        conOp1->toMatrix(K2);
                        gsInfo<<"K2:\n"<<K2<<"\n\n";
                        gsInfo<<"u:\n"<<soli<<"\n\n";
                        gsInfo<<"real residual in "<<k<<" :  "<<(temp-rhsi.topRows(info.dofsB[np]+info.dofsI[np])).norm()<<"\n"<<std::flush;

                        // gsInfo<<"error in "<<k<<" : "<<cg.error()<<"\n";
                        if(cast<int,index_t>(cg.iterations())<f)
                            break;

                        //for giving correct number of it. in the output.
                        if(k+1 == (unsigned)max_iter/f)
                            break;

                        projectToDualSubspace(soli,np);
                    }
*/
		    
                    //k*f+
                    gsInfo<<"#rhs: "<<r<<", KC - CG-MG iterations: "<<iter<<" - MG residual: "<<res<<"\n";

                }
                else
                {
                    gsMatrix<T> rhs_i=rhsi.topRows(info.dofsB[np]+info.dofsI[np]);
                    int max_iter = m_IETIoptions.opt.getInt("MG_KC.Solve.MaxIterations");
                    gsPreconditionerOp<T>* stOp =  dynamic_cast<gsPreconditionerOp<T>*>((m_inexact_KC[np]->getOperator(0,0).get()));
                    if(stOp!=NULL)
                    {
                        for( int i=0; i<max_iter;++i )
                            stOp->step(rhs_i, soli);
                    }
                    else
                    {
                        m_inexact_KC[np]->getOperator(0,0)->apply(rhs_i,soli);
                        gsMatrix<T> soli_old,temp;
                        for( int i=1; i<max_iter;++i )
                        {
                            soli_old = soli;
                            temp =rhs_i - *m_K[np]*soli_old;
                            m_inexact_KC[np]->getOperator(0,0)->apply(temp,soli_old);
                            soli += 0.9*soli_old;
                        }
                    }
                    // real_t norm = (rhsi.topRows(info.dofsB[np]+info.dofsI[np]) - mg.matrix() * soli).norm();
                    // gsInfo<<"#rhs: "<<r<<", KC - MG iterations: "<<max_iter<<"\n";
                    // gsInfo<<" - MG residual: "<<norm<<"\n";
                    // gsInfo<<"\n";
                }
                //  soli.conservativeResizeLike(rhsi);
                //sol.col(r) = soli;
                sol.block(0,r,info.dofsB[np]+info.dofsI[np],1) = soli;
            }
        }


        break;
    }
    }
    if(!isBI)//[normal] -> [B,I]
        sol.topRows(info.dofsB[np]+info.dofsI[np])=m_permMat[np]*sol.topRows(info.dofsB[np]+info.dofsI[np]);

    if(needsProjection) //Projection onto V_\Delta\Delta
        projectToDualSubspace(sol,np);

    if(m_IETIoptions.opt.getSwitch("EnforceRescaling"))
        for(int col=0;col<rhs.cols();col++)
        {
            sol.col(col) =  m_scalingsKC[np].RightScaling().cwiseProduct(sol.col(col));
            rhs.col(col) = m_scalingsKC[np].LeftScaling().cwiseQuotient(rhs.col(col));
        }

#pragma omp critical
    if(!exact)
        m_timings.KC_subspace[np]+=time.stop();
}

template< typename T>
template<bool exact>
void gsIETIAssembler<T>::solveKii(unsigned np, gsMatrix<T>& rhs, gsMatrix<T>& sol) const
{
    sol.setZero(rhs.rows(),rhs.cols());
    if(rhs.rows()==0)
        return;
    gsStopwatch time;
    if(m_IETIoptions.opt.getSwitch("EnforceRescaling"))
        for(int col=0;col<rhs.cols();col++)
            rhs.col(col) = m_scalingsKii[np].LeftScaling().cwiseProduct(rhs.col(col));

    switch (m_IETIoptions.KiiSolver)
    {
    case IETILocalSolver::Direct:
        sol=m_LU_Kii[np]->solve(rhs);
        break;
    case IETILocalSolver::FastDiagonalization:
    {
        if(exact)
        {
	    if(false && (np==1 || np==2 || np==20))
            {
                gsMatrix<T> error, eigs;
                m_KiiSolve[np]->setCalcEigenvalues(true);
                m_KiiSolve[np]->solveDetailed(rhs,sol,error);
                gsInfo<<"ErrorHistory:\n"<<error<<"\n\n";
                m_KiiSolve[np]->getEigenvalues(eigs);
                gsInfo<<"Eigenvalues:\n"<<eigs.maxCoeff()<< " - "<<eigs.minCoeff()<<"\n\n";
                m_KiiSolve[np]->setCalcEigenvalues(false);
            }
	  
	    
	  
            m_KiiSolve[np]->solve(rhs,sol);
            gsInfo<<"Kii - BS iterations: "<<m_KiiSolve[np]->iterations()<<" - BS residual: "<<m_KiiSolve[np]->error()<<"\n";
#pragma omp critical
            {
                m_timings.averageItNumKii[np]+=m_KiiSolve[np]->iterations();
                m_timings.numKiiIt[np]++;
            }
        }
        else
        {
            gsPreconditionerOp<T>* stOp =  dynamic_cast<gsPreconditionerOp<T>*>((m_inexact_Kii_BS[np].get()));
            int maxIter = m_IETIoptions.opt.getInt("MG_Kii.Prec.MaxIterations");

	/*    
            gsConjugateGradient<> cg(*m_Kii[np],m_inexact_Kii_BS[np]);
            cg.setMaxIterations(maxIter);
            cg.solve(rhs,sol);
         //   gsInfo<<"Kii - BS iterations: "<<cg.iterations()<<" - BS residual: "<<cg.error()<<"\n";
          */  

	    
            if(stOp!=NULL)
                for( int i=0; i<maxIter;++i )
                    stOp->step(rhs, sol);
            else
            {
                m_inexact_Kii_BS[np]->apply(rhs,sol);
                gsMatrix<T> sol_old;
                gsMatrix<T> temp;
                for( int i=1; i<maxIter;++i )
                {
                    sol_old = sol;
                    temp =rhs- *m_Kii[np]*sol_old;
                    m_inexact_Kii_BS[np]->apply(temp,sol_old);
                    sol += (real_t)1/m_maxEigKii[np]*sol_old;
                }

		
           //     m_inexact_Kii_BS[np]->apply(sol,temp);
           //     gsInfo<<"Kii - BS iterations: "<<maxIter<<"  ";
           //     gsInfo<<" - BS residual: "<<(rhs - temp).norm()<<"\n";
           //     gsInfo<<"\n";
		
	    }
	  
           
        }
        break;
    }
    case IETILocalSolver::Multigrid:
        if(exact)
        {
            if(false && (np==1 || np==2))
            {
                gsMatrix<T> error, eigs;
                m_KiiSolve[np]->setCalcEigenvalues(true);
                m_KiiSolve[np]->solveDetailed(rhs,sol,error);
                gsInfo<<"ErrorHistory:\n"<<error<<"\n\n";
                m_KiiSolve[np]->getEigenvalues(eigs);
                gsInfo<<"Eigenvalues:\n"<<eigs.maxCoeff()<< " - "<<eigs.minCoeff()<<"\n\n";
                m_KiiSolve[np]->setCalcEigenvalues(false);
            }
            else
                m_KiiSolve[np]->solve(rhs,sol);
#pragma omp critical
            {
                m_timings.averageItNumKii[np]+=m_KiiSolve[np]->iterations();
                m_timings.numKiiIt[np]++;
            }

            gsInfo<<"Kii - MG iterations: "<<m_KiiSolve[np]->iterations()<<" - MG residual: "<<m_KiiSolve[np]->error()<<"\n";
        }
        else
        {
            gsPreconditionerOp<T>* stOp =  dynamic_cast<gsPreconditionerOp<T>*>((m_inexact_Kii_MG[np].get()));
            int maxIter = m_IETIoptions.opt.getInt("MG_Kii.Prec.MaxIterations");
            if(stOp!=NULL)
                for( int i=0; i<maxIter;++i )
                    stOp->step(rhs, sol);
            else
                GISMO_ERROR("You chose some strange preconditioner, it should be a stepable operator");


                        //real_t norm = (rhs - m_inexact_Kii_MG[np]->matrix() * sol).norm();
                        //   gsInfo<<"Kii - MG iterations: "<<max_iter<<"\n";
                        //gsInfo<<" - MG residual: "<<norm<<"\n";
                        //   gsInfo<<"\n";
        }


        break;
    }

    if(m_IETIoptions.opt.getSwitch("EnforceRescaling"))
        for(int col=0;col<rhs.cols();col++)
        {
            sol.col(col) =  m_scalingsKii[np].RightScaling().cwiseProduct(sol.col(col));
            rhs.col(col) = m_scalingsKii[np].LeftScaling().cwiseQuotient(rhs.col(col));
        }
#pragma omp critical
    {
        if(exact)
            m_timings.Kii_schur[np]+=time.stop();
        else
            m_timings.Kii_precond[np]+=time.stop();
    }
}

template<typename T>
void gsIETIAssembler<T>::solveKrr(unsigned np, gsMatrix<T>& rhs, gsMatrix<T>& sol) const
{
    GISMO_ASSERT(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"), "Not available for this algorithm, enable NoMinimumEnergy");

    sol.setZero(rhs.rows(),rhs.cols());
    if(rhs.rows()==0)
        return;

    if(m_IETIoptions.opt.getSwitch("EnforceRescaling"))
        for(int col=0;col<rhs.cols();col++)
            rhs.col(col) = m_scalingsKrr[np].LeftScaling().cwiseProduct(rhs.col(col));

    switch (m_IETIoptions.KrrSolver)
    {
    case IETILocalSolver::Direct:
        sol=m_LU_Krr[np]->solve(rhs);
        break;
    default:
        GISMO_ERROR("Inexact solvers not implemented for Krr");
    }

    if(m_IETIoptions.opt.getSwitch("EnforceRescaling"))
        for(int col=0;col<rhs.cols();col++)
        {
            sol.col(col) =  m_scalingsKrr[np].RightScaling().cwiseProduct(sol.col(col));
            rhs.col(col) = m_scalingsKrr[np].LeftScaling().cwiseQuotient(rhs.col(col));
        }
}

template<typename T>
void gsIETIAssembler<T>::solveSpp(const gsMatrix<T>& rhs, gsMatrix<T>& sol) const
{
    sol.setZero(rhs.rows(),rhs.cols());
    if(rhs.rows()==0)
        return;

    sol=m_LU_Spp->solve(rhs);
}

} // namespace gismo
