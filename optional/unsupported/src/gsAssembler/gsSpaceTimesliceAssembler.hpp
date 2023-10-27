/**  gsSpaceTimesliceAssembler.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    Created on:  2017-06-06
*/

#pragma once
#include <gsIO/gsOptionList.h>
#include <gsAssembler/gsSpaceTimesliceAssembler.h>
#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>
#include <gsPde/gsSpaceTimePoissonPde.h>
#include <gsPde/gsPoissonHeterogeneousPde.h>
#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsSolver/gsKronecker.h>
#include <gsAssembler/gsVisitorSpaceTimeSlice.h>
#include <gsUtils/gsStopwatch.h>
#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsAssembler/gsGenericAssembler.h>


namespace gismo {

///Extract the part of the geometry \a inp which points in \a desiredDirection
template<typename T>
void gsSpaceTimesliceAssembler<T>::reduceSpline(gsGeometry<T> & inp, int desiredDirection, gsBSpline<T>& result)
{
    int d = inp.parDim();
    switch(d)
    {
    case 2:
    {
        typename gsTensorBSpline<2, T>::uPtr tb = memory::convert_ptr<gsTensorBSpline<2, T> >(inp.clone());
        tb->swapDirections(0,desiredDirection);
        internal::Slicer<2,T>::reduceSpline(*tb, result);
        break;
    }
    case 3:
    {
        typename gsTensorBSpline<3,T>::uPtr tb = memory::convert_ptr<gsTensorBSpline<3,T> >(inp.clone());
        tb->swapDirections(0,desiredDirection);
        internal::Slicer<3,T>::reduceSpline(*tb,result);
        break;
    }
    case 4:
    {
        typename gsTensorBSpline<4,T>::uPtr tb = memory::convert_ptr<gsTensorBSpline<4,T> >(inp.clone());
        tb->swapDirections(0,desiredDirection);
        internal::Slicer<4,T>::reduceSpline(*tb,result);
        break;
    }
    }
    //embedding does not work here
    result.setCoefs(result.coefs().col(d-1));
}



template<typename T>
void gsSpaceTimesliceAssembler<T>::refresh()
{
    // Check for coherency
    GISMO_ASSERT(this->check(), "Incoherent data in Assembler");

    GISMO_ASSERT(1==m_bases.size(), "Expecting a single discrete space "
                                    "for standard scalar Galerkin");

    m_TimeSliceForPatch.clear(); m_tops.clear();m_Spatial_Domain.clear();
    m_Spatial_Bases.clear();m_Spatial_boundary_conditions.clear();m_Temporal_Domain.clear();
    m_Temporal_Bases.clear();m_dGInterfaces.clear();m_patchIndicesForTimeSlice.clear();
    m_hasAssembledTP = false;
    std::vector<patchSide> lowerBounds;
    if(!m_timeIsLastDirection)
    {
        const typename gsBoundaryConditions<T>::bcContainer init = m_pde_ptr->bc().container("Initial");
        lowerBounds.resize(init.size());
        for(size_t i=0; i<init.size();++i)
            lowerBounds[i] = patchSide(init[i].patch(),init[i].side());
        std::sort(lowerBounds.begin(),lowerBounds.end());
    }
    else
    {
        real_t tVal = m_pde_ptr->domain().patch(0).coefAtCorner(boxCorner(0))(m_pde_ptr->domain().patch(0).parDim()-1);
        for(size_t np=0; np<m_pde_ptr->domain().nPatches();++np)
        {
            if(m_pde_ptr->domain().patch(np).coefAtCorner(boxCorner(0))(m_pde_ptr->domain().patch(0).parDim()-1) == tVal)
                lowerBounds.push_back(patchSide(np,boxSide(m_pde_ptr->domain().patch(0).parDim()-1,0)));
        }
    }
    m_TimeSliceForPatch.resize(m_pde_ptr->domain().nPatches());
    m_tops.reserve(lowerBounds.size());
    for(size_t i=0; i< lowerBounds.size();++i)
    {
        patchSide& bc_ = lowerBounds[i];
        switch(m_pde_ptr->domain().parDim())
        {
        case 2:
        {
            gsTensorBSpline<2,T>& geo = dynamic_cast<gsTensorBSpline<2,T>& > (m_pde_ptr->domain().patch(bc_.patch));
            typename gsTensorBSpline<2,T>::BoundaryGeometryType tb_new;
            geo.slice(1,0,tb_new);
            tb_new.embed(1);
            m_Spatial_Domain.addPatch(tb_new);
            break;
        }
        case 3:
        {
            gsTensorBSpline<3,T>& geo = dynamic_cast<gsTensorBSpline<3,T>& > (m_pde_ptr->domain().patch(bc_.patch));
            typename gsTensorBSpline<3,T>::BoundaryGeometryType tb_new;
            geo.slice(2,0,tb_new);
            tb_new.embed(2);
            m_Spatial_Domain.addPatch(tb_new);
            break;
        }
        case 4:
        {
            gsTensorBSpline<4,T>& geo = dynamic_cast<gsTensorBSpline<4,T>& > (m_pde_ptr->domain().patch(bc_.patch));
            typename gsTensorBSpline<4,T>::BoundaryGeometryType tb_new;
            geo.slice(3,0,tb_new);
            tb_new.embed(3);
            m_Spatial_Domain.addPatch(tb_new);
            break;
        }
        };

        //Assume the same boundary conditions for all timesteps
        gsBoundaryConditions<T> bc_spatial_;
        m_pde_ptr->bc().getConditionsForPatch(bc_.patch,bc_spatial_);
        for(typename gsBoundaryConditions<T>::const_iterator it= bc_spatial_.begin("Neumann"); it!=bc_spatial_.end("Neumann");++it)
            m_Spatial_boundary_conditions.add(i,it->side(),it->ctype(),it->function(),it->unknown(),it->unkComponent(),it->parametric());
        for(typename gsBoundaryConditions<T>::const_iterator it= bc_spatial_.begin("Dirichlet"); it!=bc_spatial_.end("Dirichlet");++it)
            m_Spatial_boundary_conditions.add(i,it->side(),it->ctype(),it->function(),it->unknown(),it->unkComponent(),it->parametric());

        patchSide side;
        patchSide result = bc_;
        int slice = 0;
        do
        {
            if(i==0)
            {
                //create the right size of the Basis vector (size == nSlices)
                m_Spatial_Bases.push_back(gsMultiBasis<T>());
                m_patchIndicesForTimeSlice.push_back(std::vector<index_t>());

                //Make the temporal basis and domain
                gsBSpline<T> TimeDomain;
                reduceSpline(m_pde_ptr->domain().patch(result.patch),result.direction(),TimeDomain);

                m_Temporal_Domain.addPatch(TimeDomain.clone());
                m_Temporal_Bases.addBasis(m_bases.front().basis(result.patch).component(result.direction()).clone());
            }

            //Make the spatial multibasis for each time slice
            m_Spatial_Bases[slice].addBasis(m_bases.front().basis(result.patch).boundaryBasis(result.side()).release());
            m_patchIndicesForTimeSlice[slice].push_back(result.patch);
            m_TimeSliceForPatch[result.patch] =slice;

            //save the old side
            side=result;
            slice ++;
        }
        while(m_pde_ptr->domain().getNeighbour(patchSide(side.patch,side.opposite()),result)); // move on in time
        m_tops.push_back(patchSide(side.patch,side.opposite()));
    }

    m_Spatial_Domain.computeTopology();
    m_Temporal_Domain.computeTopology();
    gsMultiBasis<T> b;
    for(size_t slice = 0; slice < m_Temporal_Domain.nPatches();++slice)
        m_Spatial_Bases[slice].setTopology(m_Spatial_Domain.topology());
    m_Temporal_Bases.setTopology(m_Temporal_Domain.topology());

    for(typename gsMultiPatch<T>::iiterator it = m_pde_ptr->domain().iBegin(); it!=m_pde_ptr->domain().iEnd();++it)
    {
        if(m_TimeSliceForPatch[it->first().patch]!=m_TimeSliceForPatch[it->second().patch])
            m_dGInterfaces.push_back(*it);
    }


    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper mapper;
    m_bases.front().getMapper(
                (dirichlet::strategy)(m_options.getInt("DirichletStrategy")),
                (iFace::none),
                this->pde().bc(), mapper, 0,false);


    for(typename gsMultiPatch<T>::iiterator it= m_pde_ptr->domain().iBegin();it!=m_pde_ptr->domain().iEnd();++it)
    {
        boundaryInterface bI = *it;
        if(m_TimeSliceForPatch[bI.first().patch] == m_TimeSliceForPatch[bI.second().patch])
            m_bases.front().matchInterface(bI, mapper);
    }
    mapper.finalize();
    if ( 0 == mapper.freeSize() ) // Are there any interior dofs ?
        gsWarn << " No internal DOFs, zero sized system.\n";


    //permuation for creating a tensor product representation for each time-slice
    //and create also the nunmber of global elements in each slice
    gsVector<index_t> perm = gsVector<index_t>::Constant(mapper.freeSize(),-1);
    m_sizes.setZero(m_Temporal_Bases.nBases());
    int iter = 0;
    for(size_t slice=0; slice<m_patchIndicesForTimeSlice.size();++slice)
    {
        gsDofMapper spatialMapper;
        m_Spatial_Bases[slice].getMapper(
                    (dirichlet::strategy)(m_options.getInt("DirichletStrategy")),
                    (iFace::strategy)m_options.getInt("InterfaceStrategy"),
                    m_Spatial_boundary_conditions, spatialMapper, 0,true);



        for(int t = 0; t < m_Temporal_Bases.basis(slice).size();++t)
        {
            //   int iterC = 0;
            for(size_t s = 0; s< m_patchIndicesForTimeSlice[slice].size();++s)
            {
                int size = m_Spatial_Bases[slice].basis(s).size();
                for(int i=t*size; i<(t+1)*size;++i)
                {
                    index_t mapped  = mapper.index(i,m_patchIndicesForTimeSlice[slice][s]);
                    index_t spatialMapped = spatialMapper.index(i%size,s);
                    if(mapper.is_free_index(mapped) && perm(mapped)==-1)
                    {
                        perm(mapped) = spatialMapped + iter;
                        m_sizes[slice]++;

                    }
                }
            }
            iter+= spatialMapper.freeSize();
        }
    }

    mapper.markCoupledAsTagged();
    mapper.permuteFreeDofs(perm);

    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(mapper);//1,1

    //   m_permToTensor = Permutation(perm);

    //get a rough estimate for the physical meshsize in time (tau) and space (h)
    m_h.resize(m_patchIndicesForTimeSlice.size()), m_tau.resize(m_patchIndicesForTimeSlice.size());
    for(size_t slice=0; slice<m_patchIndicesForTimeSlice.size();++slice)
    {
        T len = (m_Temporal_Domain.patch(slice).coefAtCorner(boxCorner::getFirst(1)) - m_Temporal_Domain.patch(slice).coefAtCorner(boxCorner::getLast(1))).norm();
        m_tau[slice] = m_Temporal_Bases.basis(slice).getMaxCellLength()*len;

        len = (m_Spatial_Domain.patch(0).coefAtCorner(boxCorner::getFirst(1)) - m_Spatial_Domain.patch(0).coefAtCorner(boxCorner::getLast(1))).norm();
        m_h[slice] = m_Spatial_Bases[slice].basis(0).getMaxCellLength()*len;
    }

    m_patchesForAssembling.resize(m_Temporal_Domain.nPatches());
    for(size_t slice=0; slice< m_Temporal_Domain.nPatches();++slice)
        m_patchesForAssembling[slice].resize(m_Spatial_Domain.nPatches(),true);

}


template<typename T>
void gsSpaceTimesliceAssembler<T>::assembleTPMats(bool assemblePatchwise)
{
    gsStopwatch time;
    m_hasAssembledTP = true;

    gsSpaceTimePoissonPde<T>* STpde= dynamic_cast<gsSpaceTimePoissonPde<T>* >(m_pde_ptr.get());
    gsPoissonHeterogeneousPde<T> spatialPde(m_Spatial_Domain,m_Spatial_boundary_conditions,*STpde->rhs_x(),*STpde->getAlpha());
    gsOptionList optTime = gsAssembler<T>::defaultOptions();
    optTime.setInt("DirichletStrategy" , dirichlet::none );
    optTime.setInt("InterfaceStrategy", iFace::none  );
    gsOptionList opt = gsAssembler<T>::defaultOptions();
    opt.setInt("DirichletStrategy" , dirichlet::none );
    opt.setInt("InterfaceStrategy", iFace::glue  );
    int nSlice = m_Temporal_Domain.nPatches();
    int sliceCount = 0;
    //m_TPData.M1t.resize(nSlice); m_TPData.K1t.resize(nSlice);
    //m_TPData.Mht.resize(nSlice); m_TPData.Kht.resize(nSlice);
    m_TPData.Mx.resize(nSlice); m_TPData.Mt.resize(nSlice); m_TPData.Kt.resize(nSlice);
    m_TPData.Kx.resize(nSlice); m_TPData.MWx.resize(nSlice); m_TPData.Nt.resize(nSlice);
    if(assemblePatchwise)
    {
        m_TPData.PatchKx.resize(nSlice);m_TPData.PatchMx.resize(nSlice);m_TPData.PatchMWx.resize(nSlice);
    }
    gsMatrix<T> rhs;
    rhs.setZero(this->system().colMapper(0).size(),1);

    typedef typename gsExprAssembler<T>::variable variable;
    typedef typename gsExprAssembler<T>::space    space   ;
    typedef typename gsExprAssembler<T>::geometryMap geoMap;
    gsBoundaryConditions<T> bc;
    std::vector<typename gsMatrix<T>::uPtr> fxs(m_Temporal_Bases.nBases());
    //  std::vector<gsMatrix<T> > rhs(m_Temporal_Bases.nBases());



    for(size_t n = 0; n< m_Temporal_Bases.nBases(); ++n)
    {
        gsGenericAssembler<T> genSpatial(m_Spatial_Domain,m_Spatial_Bases[n],opt,&m_Spatial_boundary_conditions);
        gsPoissonHeterogeneousAssembler<T> spatialAssembler(spatialPde,m_Spatial_Bases[n],dirichlet::none,iFace::glue);
        spatialAssembler.options().setInt("DirichletValues",dirichlet::homogeneous); //just something, that the boudary condition is not evaluated.
        gsSparseSystem<T> sysMWX;

        gsGenericAssembler<T> genTemporal(m_Temporal_Domain.patch(n),m_Temporal_Bases.basis(n),optTime);

        if(assemblePatchwise)
        {
            m_TPData.PatchKx[n].resize(m_Spatial_Domain.nPatches());m_TPData.PatchMx[n].resize(m_Spatial_Domain.nPatches());m_TPData.PatchMWx[n].resize(m_Spatial_Domain.nPatches());
        }
        if(false && std::count(m_patchesForAssembling[n].begin(), m_patchesForAssembling[n].end(), true)==0)
        {
            sliceCount+=genTemporal.numDofs()*spatialAssembler.numDofs();
            continue;
        }
        real_t h = m_tau[n];
        //Calculate the temporal matrices
        time.restart();

        gsSparseMatrix<T> Kt = genTemporal.assembleStiffness();
        gsSparseMatrix<T> Mt = genTemporal.assembleMass();
        timings(1)+=time.stop();

        // calculate the point evaluation at t_n-1 -> 0
        gsSparseMatrix<T> Ptnn(m_Temporal_Bases.basis(n).size(),m_Temporal_Bases.basis(n).size());
        Ptnn.coeffRef(0,0) = 1;

        // calculate the mixed point evaluation at t_n-1 for Basis[n] and Basis[n-1]
        //calculate the mixed spatial integeral on the timeslap interface for Basis[n] and Basis[n-1]
        if(n>0)
        {
            time.restart();
            gsSparseMatrix<T> Ptn1n(m_Temporal_Bases.basis(n).size(),m_Temporal_Bases.basis(n-1).size());
            Ptn1n.coeffRef(0,m_Temporal_Bases.basis(n-1).size()-1) = 1;
            m_TPData.Nt[n] = Ptn1n.moveToPtr();

            gsExprAssembler<T> massAssemblerI(2,2);
            geoMap G1 = massAssemblerI.getMap(m_Spatial_Domain);

            space un = massAssemblerI.getSpace(m_Spatial_Bases[n], 1, 0); // args: splines, bc, dim, id
            un.addBc( bc.get("Dirichlet",0) );
            space un1 = massAssemblerI.getSpace(m_Spatial_Bases[n-1], 1, 1); // args: splines, bc, dim, id
            un1.addBc( bc.get("Dirichlet",1) );

            m_Spatial_Bases[n].totalSize()>m_Spatial_Bases[n-1].totalSize() ?
                        massAssemblerI.setIntegrationElements(m_Spatial_Bases[n]) :massAssemblerI.setIntegrationElements(m_Spatial_Bases[n-1]);
            massAssemblerI.setOptions(opt);
            massAssemblerI.initSystem();
            if(!assemblePatchwise)
            {
                massAssemblerI.assemble(un*un1.tr()*meas(G1));
                m_TPData.MWx[n] = gsSparseMatrix<T>(massAssemblerI.matrixBlockView()(0,1)).moveToPtr();
            }
            timings(6) = time.stop();

            std::vector<gsDofMapper> mappers(2);
            m_Spatial_Bases[n-1].getMapper(true,mappers[0],true);
            m_Spatial_Bases[n].getMapper(true,mappers[1],true);
            gsVector<index_t> dims = gsVector<index_t>::Constant(1,1);
            sysMWX = gsSparseSystem<T>(mappers,dims);
            sysMWX.reserve(m_Spatial_Bases[n],opt,1);
        }

        gsExprAssembler<T> dtuvAssembler(1,1);

        time.restart();
        gsMultiBasis<T> basis(m_Temporal_Bases.basis(n));
        geoMap G = dtuvAssembler.getMap(m_Temporal_Domain.patch(n));
        space u = dtuvAssembler.getSpace(basis, 1, 0); // args: splines, bc, dim, id
        u.addBc( bc.get("Dirichlet",0) );
        variable ft = dtuvAssembler.getCoeff(*STpde->rhs_t(),G); // force
        dtuvAssembler.setOptions(optTime);

        dtuvAssembler.setIntegrationElements(basis); // subdivided mesh
        dtuvAssembler.initSystem();
        dtuvAssembler.assemble(u * igrad(u,G).tr()*meas(G),(u+m_theta*h*igrad(u,G))*ft*meas(G)); //

        timings(2)+=time.stop();

        gsSparseMatrix<T> Tt = dtuvAssembler.matrix(); //this is really Dtu*v
        const gsMatrix<T>& f_t = dtuvAssembler.rhs();
        //m_TPData.M1t[n] = Mt.moveToPtr();
        //m_TPData.Mht[n] = (gsSparseMatrix<T>(m_theta*h*Tt.transpose())).moveToPtr();
        //m_TPData.K1t[n] = (gsSparseMatrix<T>(Tt + Ptnn)).moveToPtr();
        //m_TPData.Kht[n] = (gsSparseMatrix<T>(m_theta*h*Kt)).moveToPtr();

//      m_TPData.Kt[n] = (gsSparseMatrix<T>(*m_TPData.K1t[n]+*m_TPData.Kht[n])).moveToPtr();
//      m_TPData.Mt[n] = (gsSparseMatrix<T>(*m_TPData.M1t[n]+*m_TPData.Mht[n])).moveToPtr();

        m_TPData.Kt[n] = (gsSparseMatrix<T>(Tt + Ptnn+m_theta*h*Kt)).moveToPtr();
        m_TPData.Mt[n] = (gsSparseMatrix<T>(Mt+m_theta*h*gsSparseMatrix<T>(Tt.transpose()))).moveToPtr();


        //Spatial Matrices
        if(!assemblePatchwise)
        {
            time.restart();
            spatialAssembler.assemble();
            timings(3)+=time.stop();

            time.restart();
            genSpatial.assembleMass();
            timings(4)+=time.stop();
        }



        if(assemblePatchwise)
        {
            spatialAssembler.system().rhs().setZero(spatialAssembler.system().colMapper(0).freeSize(),1);
            genSpatial.system().rhs().setZero(spatialAssembler.system().colMapper(0).freeSize(),1);
            spatialAssembler.computeDirichletDofs(0);

            spatialAssembler.system().reserve(m_Spatial_Bases[n],opt,1);
            genSpatial.system().reserve(m_Spatial_Bases[n],opt,1);

            if(n>0)sysMWX.rhs().setZero(sysMWX.colMapper(0).freeSize(),1);
            for(size_t np=0; np<m_Spatial_Domain.nPatches();++np)
            {
                if(m_patchesForAssembling[n][np]==false)
                    continue;

                gsPoissonHeterogeneousPde<T>* patchPde = dynamic_cast<gsPoissonHeterogeneousPde<T>*>(spatialPde.restrictToPatch(np));
                gsPoissonHeterogeneousAssembler<T> patchSpatialAssembler(*patchPde,m_Spatial_Bases[n].basis(np),dirichlet::none,iFace::glue);
                patchSpatialAssembler.options().setInt("DirichletValues",dirichlet::homogeneous);

                gsBoundaryConditions<T> patchBC;
                spatialPde.bc().getConditionsForPatch(np,patchBC);
                gsGenericAssembler<T> patchGenSpatial(m_Spatial_Domain.patch(np),m_Spatial_Bases[n].basis(np),opt,&patchBC);
                patchGenSpatial.options().setInt("DirichletValues",dirichlet::homogeneous);
                patchGenSpatial.options().setInt("DirichletStrategy",dirichlet::none);
                if(n>0)
                {
                    time.restart();
                    gsExprAssembler<T> patchMassAssemblerI(2,2);
                    gsMultiBasis<T> bn = gsMultiBasis<T>(m_Spatial_Bases[n][np]);
                    gsMultiBasis<T> bn1 = gsMultiBasis<T>(m_Spatial_Bases[n][np]);
                    gsMultiPatch<T> p = gsMultiPatch<T>(m_Spatial_Domain.patch(np));

                    geoMap G1 = patchMassAssemblerI.getMap(p);
                    /*
                    variable un = patchMassAssemblerI.setSpace(bn,patchBC , 1, 0); // args: splines, bc, dim, id
                    gsBoundaryConditions<T> patchBC2 = patchBC;
                    for(typename gsBoundaryConditions<T>::iterator dIt = patchBC2.dirichletBegin();dIt!=patchBC2.dirichletEnd();++dIt)
                        dIt->m_unknown  =1;

                    space un1 = patchMassAssemblerI.getSpace(bn1,1, 1); // args: splines, bc, dim, id
                    un1.addBc( patchBC.get("Dirichlet",1) );
                    */
                    space un = patchMassAssemblerI.getSpace(bn , 1, 0); // args: splines, bc, dim, id
                    space un1 = patchMassAssemblerI.getSpace(bn1 , 1, 1); // args: splines, bc, dim, id

                    patchMassAssemblerI.setOptions(opt);
                    patchMassAssemblerI.options().setInt("DirichletValues",dirichlet::homogeneous);
                    patchMassAssemblerI.options().setInt("DirichletStrategy",dirichlet::none);
                    m_Spatial_Bases[n][np].size()>m_Spatial_Bases[n-1][np].size() ?
                                patchMassAssemblerI.setIntegrationElements(bn) :patchMassAssemblerI.setIntegrationElements(bn1);

                    patchMassAssemblerI.initSystem();
                    patchMassAssemblerI.assemble(un*un1.tr()*meas(G1));
                    m_TPData.PatchMWx[n][np] = gsSparseMatrix<T>(patchMassAssemblerI.matrixBlockView()(0,1)).moveToPtr();

                    gsMatrix<index_t> actives1(m_TPData.PatchMWx[n][np]->rows(),1);
                    gsMatrix<index_t> actives2(m_TPData.PatchMWx[n][np]->cols(),1);
                    size_t iter=0;
                    for(index_t i =0; i<m_Spatial_Bases[n-1][np].size();++i)
                        if(sysMWX.colMapper(0).is_free(i,np))
                            actives1(iter++,0)=i;
                    iter =0;
                    for(index_t i =0; i<m_Spatial_Bases[n][np].size();++i)
                        if(sysMWX.rowMapper(0).is_free(i,np))
                            actives2(iter++,0)=i;

                    sysMWX.mapRowIndices(actives1,np,actives1,0);
                    sysMWX.mapColIndices(actives2,np,actives2,0);
                    sysMWX.pushSparse(*m_TPData.PatchMWx[n][np],gsMatrix<T>::Zero(m_TPData.PatchMWx[n][np]->rows(),1),actives1,actives2,spatialAssembler.fixedDofs());
                    timings(6)+=time.stop();
                }

                time.restart();
                patchSpatialAssembler.assemble();
                timings(3)+=time.stop();
                patchGenSpatial.assembleMass();
                timings(4)+=time.stop();

                m_TPData.PatchKx[n][np] = patchSpatialAssembler.system().matrix().moveToPtr();
                m_TPData.PatchMx[n][np] = patchGenSpatial.system().matrix().moveToPtr();

                gsMatrix<index_t> actives(patchSpatialAssembler.system().colMapper(0).freeSize(),1);
                size_t iter=0;
                for(size_t i =0; i<patchSpatialAssembler.system().colMapper(0).mapSize();++i)
                    if(patchSpatialAssembler.system().colMapper(0).is_free(i))
                        actives(iter++,0)=i;

                spatialAssembler.system().mapColIndices(actives,np,actives,0);

                time.restart();
                spatialAssembler.system().pushSparse(*m_TPData.PatchKx[n][np],patchSpatialAssembler.rhs(),actives,actives,spatialAssembler.fixedDofs());
                timings(3)+=time.stop();
                time.restart();
                genSpatial.system().pushSparse(*m_TPData.PatchMx[n][np],genSpatial.rhs(),actives,actives,genSpatial.fixedDofs());
                timings(4)+=time.stop();


                delete patchPde;

            }
        }



        m_TPData.Kx[n] = spatialAssembler.system().matrix().moveToPtr();
        m_TPData.Mx[n] = genSpatial.system().matrix().moveToPtr();
        if(assemblePatchwise)
            m_TPData.MWx[n] = sysMWX.matrix().moveToPtr();

        fxs[n] = spatialAssembler.system().rhs().moveToPtr();
        //Build the right hand side;

        time.restart();
        gsMatrix<T> temp = f_t.kron(*fxs[n]);


        if(n==0)
        {
            spatialAssembler.system().rhs().setZero( m_TPData.Kx[n]->rows(),1);

            for(size_t np=0; np<m_Spatial_Domain.nPatches();++np)
            {
                if(m_patchesForAssembling[n][np]==false)
                    continue;

                //TODO: this assumes that the time is always up!
                gsMultiPatch<T> mp(m_Spatial_Domain[np]);
                mp.embed(m_Spatial_Domain.parDim()+1);
                gsVector<T> transl;
                transl.setZero(m_Spatial_Domain.parDim()+1);
                transl(m_Spatial_Domain.parDim()) = m_Temporal_Domain.patch(0).coef(0,0);
                mp.patch(0).translate(transl);

                /*
                    gsBoundaryConditions<T> BC;
                    const typename gsBoundaryConditions<T>::bcContainer initBC = m_pde_ptr->bc().container("Initial");
                    for(size_t ii = 0; ii<initBC.size();++ii)
                    {
                        boundary_condition<T> bc_ = initBC[ii];
                        BC.add(ii,bc_.side(),bc_.ctype(),bc_.function(),bc_.unknown(),bc_.unkComponent(),bc_.parametric());
                    }
    */
                gsExprAssembler<T> ass(1,1);
                gsMultiBasis<T> bb(m_Spatial_Bases[n].basis(np));
                geoMap GS = ass.getMap(mp);

                space uS = ass.getSpace(bb , 1, 0); // args: splines, bc, dim, id
                variable ftS = ass.getCoeff(*m_pde_ptr->bc().begin("Initial")->function(),GS); // force
                gsOptionList optft = gsAssembler<T>::defaultOptions();
                optft.setInt("DirichletValues" , dirichlet::l2Projection );
                //optft.setInt("DirichletStrategy" , dirichlet::homogeneous );
                optft.setInt("DirichletStrategy" , dirichlet::none );
                optft.setInt("InterfaceStrategy", iFace::glue  );
                ass.setOptions(optft);

                ass.setIntegrationElements(bb); // subdivided mesh
                ass.initSystem();

                ass.assemble( uS * ftS*meas(GS));

                gsMatrix<index_t> actives(ass.rhs().rows(),1);
                size_t iter=0;
                for(index_t i =0; i<m_Spatial_Bases[n].basis(np).size();++i)
                    if(spatialAssembler.system().colMapper(0).is_free(i,np))
                        actives(iter++,0)=i;

                spatialAssembler.system().mapColIndices(actives,np,actives,0);
                spatialAssembler.system().pushToRhs(ass.rhs(),actives,0);
            }
            const gsMatrix<T> & init = spatialAssembler.system().rhs();
            gsMatrix<T> unit;
            unit.setZero(m_Temporal_Bases.size(n),1);
            unit(0,0) = 1;

            gsMatrix<T> result = unit.kron(init);
            temp+=result; //Add the inital condition
        }
        //gsInfo<<"K on Patch: "<<np<<" in slice "<<n<<"\n"<<Ks[n]->toDense()<<"\n";
        //gsInfo<<"M on Patch: "<<np<<" in slice "<<n<<"\n"<<Mxs[n]->toDense()<<"\n";
        //gsInfo<<"fx on Patch: "<<np<<" in slice "<<n<<"\n"<<*fxs[n]<<"\n";


        rhs.block(sliceCount,0,temp.rows(),1) = temp;
        sliceCount+=temp.rows();
        timings(5)+=time.stop();



        //gsInfo<<"KMt on Patch: "<<np<<" in slice "<<n<<"\n"<<KMt.toDense()<<"\n";
        //gsInfo<<"MxTt on Patch: "<<np<<" in slice "<<n<<"\n"<<MxTt.toDense()<<"\n";
        //gsInfo<<"KTt on Patch: "<<np<<" in slice "<<n<<"\n"<<KTt.toDense()<<"\n";
        //gsInfo<<"F on Patch: "<<np<<" in slice "<<n<<"\n"<<F<<"\n";

    }
    time.restart();
    //When eliminating the Boundary conditions, one has to be carefull when they depend on t.
    //Here we do the ellimination of the dofs manually.
    gsDofMapper mapperG = this->system().colMapper(0);
    std::vector<gsDofMapper> mapperElims(m_Temporal_Bases.nBases());
    std::vector<gsDofMapper> mappers(m_Temporal_Bases.nBases());
    m_Spatial_Bases[0].getMapper(dirichlet::elimination,iFace::glue,m_Spatial_boundary_conditions,mapperElims[0],0);
    sliceCount = 0;

    std::vector<std::vector< std::vector<std::pair<index_t,index_t> > > >preImages(m_Temporal_Bases.nBases());

    for(size_t n = 0; n< m_Temporal_Bases.nBases(); ++n)
    {
        gsPoissonHeterogeneousAssembler<T> spatialAssembler(spatialPde,m_Spatial_Bases[n],dirichlet::none,iFace::glue);
        spatialAssembler.options().setInt("DirichletValues",dirichlet::homogeneous); //just something, that the boudary condition is not evaluated.

        if(false && std::count(m_patchesForAssembling[n].begin(), m_patchesForAssembling[n].end(), true)==0)
        {
            sliceCount+=m_Temporal_Bases.basis(n).size()*spatialAssembler.numDofs();
            continue;
        }

        mappers[n] = spatialAssembler.system().colMapper(0);
        preImages[n].resize(mappers[n].size());

        std::vector<bool> elimCols;
        elimCols.resize(mapperElims[n].size()*m_Temporal_Bases[n].size(),false);

        for(size_t np = 0; np<m_patchIndicesForTimeSlice[n].size();++np)
            for(int i=0; i<m_Spatial_Bases[n].basis(np).size();++i)
                mappers[n].preImage(mappers[n].index(i,np),preImages[n][mappers[n].index(i,np)]); //here is np, because we use the spatial mapper

        gsDofMapper mapper1;
        if(n<m_Temporal_Bases.nBases()-1)
        {
            gsPoissonHeterogeneousAssembler<T> spatialAssembler1(spatialPde,m_Spatial_Bases[n+1],dirichlet::none,iFace::glue);
            spatialAssembler1.options().setInt("DirichletValues",dirichlet::homogeneous); //just something, that the boudary condition is not evaluated.
            m_Spatial_Bases[n+1].getMapper(dirichlet::elimination,iFace::glue,m_Spatial_boundary_conditions,mapperElims[n+1],0);
            mapper1 = spatialAssembler.system().colMapper(0);
        }



        time.restart();
        for(size_t np = 0; np<m_patchIndicesForTimeSlice[n].size();++np)
        {
            int size = m_Spatial_Bases[n].basis(np).size();
            for(int i=0; i<size;++i)
            {

                index_t gE = mapperElims[n].index(i,np);
                if(mapperElims[n].is_boundary_index(gE))
                {
                    index_t g = mappers[n].index(i,np);
                    if (elimCols[g] ==true)
                        continue;

                    elimCols[g]=true;

                    gsSparseMatrix<T> Acol, Bcol;
                    gsSparseMatrix<T> temp1;
                    Acol  = m_TPData.Kt[n]->kron(gsSparseMatrix<T>(m_TPData.Mx[n]->col(g)));
                    temp1 = m_TPData.Mt[n]->kron(gsSparseMatrix<T>(m_TPData.Kx[n]->col(g)));
                    Acol+=temp1;


                    for (int k=0; k<Acol.outerSize(); ++k)
                        for (typename gsSparseMatrix<T>::InnerIterator it(Acol,k); it; ++it)
                        {
                            int idxGB = mapperG.bindex(m_Spatial_Bases[n].basis(np).size()*it.col()+i,m_patchIndicesForTimeSlice[n][np]);

                            rhs(sliceCount+it.row(),0)-= m_ddof.front()(idxGB)*it.value(); //f[n]
                        }


                    if(n<m_Temporal_Bases.nBases()-1)
                    {
                        Bcol = (m_TPData.Nt[n+1])->kron(gsSparseMatrix<T>(m_TPData.MWx[n+1]->col(g)));


                        for (int k=0; k<Bcol.outerSize(); ++k)
                            for (typename gsSparseMatrix<T>::InnerIterator it(Bcol,k); it; ++it)
                            {
                                int idxGB = mapperG.bindex(m_Spatial_Bases[n].basis(np).size()*it.col()+i,m_patchIndicesForTimeSlice[n][np]);
                                //corresponds to n timeslap but entry f[n+1]. The block is given by -B!!!!
                                rhs(sliceCount+Acol.rows()+it.row(),0)-= m_ddof.front()(idxGB,0)*(-it.value());
                            }
                    }

                }
            }

        }
        timings(7)+= time.stop();
        //Resize the spatial matrices

        //Reorder the eliminated Dofs last
        time.restart();
        gsVector<index_t> perm (mapperElims[n].size());
        std::vector<std::pair<index_t, index_t> > result;
        int iterF, iterB;
        iterF = iterB=0;
        for(int i=0; i<perm.rows();++i)
        {
            //mappers[n].preImage(i,result);
            result = preImages[n][i];

            if(mapperElims[n].is_free(result.front().second,result.front().first))
                perm(i) = iterF++;
            else
                perm(i) = mapperElims[n].freeSize()+iterB++;
        }


        Permutation permut(perm);
        gsSparseMatrix<T> temp1;
        temp1=   m_TPData.Mx[n]->twistedBy(permut);
        temp1.uncompress();
        temp1.conservativeResize(mapperElims[n].freeSize(),mapperElims[n].freeSize());
        temp1.makeCompressed();
        m_TPData.Mx[n] = temp1.moveToPtr();
        temp1 =  m_TPData.Kx[n]->twistedBy(permut);
        temp1.uncompress();
        temp1.conservativeResize(mapperElims[n].freeSize(),mapperElims[n].freeSize()) ;
        m_TPData.Kx[n] =temp1.moveToPtr();
        temp1.makeCompressed();

        if(n<m_Temporal_Bases.nBases()-1)
        {
            temp1 = ((*m_TPData.MWx[n+1])*permut.inverse());
            temp1.uncompress();
            temp1.conservativeResize(m_TPData.MWx[n+1]->rows(),mapperElims[n].freeSize());
            temp1.makeCompressed();
            m_TPData.MWx[n+1] = temp1.moveToPtr();
        }
        if(n>0)
        {
            temp1 = permut*(*m_TPData.MWx[n]);
            temp1.uncompress();
            temp1.conservativeResize(mapperElims[n].freeSize(),m_TPData.MWx[n]->cols());
            temp1.makeCompressed();
            m_TPData.MWx[n] = temp1.moveToPtr(); //(gsSparseMatrix<T>(permut))*(*MWx[n+1]);
            //  MWx[n]->applyOnTheLeft(permut);

        }




        if(assemblePatchwise)
        {
            for(size_t np = 0; np< m_Spatial_Domain.nPatches();++np)
            {
                if(m_patchesForAssembling[n][np]==false)
                    continue;

                gsPoissonHeterogeneousPde<T>* patchPde = dynamic_cast<gsPoissonHeterogeneousPde<T>*>(spatialPde.restrictToPatch(np));
                gsMultiBasis<T>mb(m_Spatial_Bases[0][np]);
                gsDofMapper mapper = mb.getMapper(dirichlet::elimination,iFace::none,patchPde->boundaryConditions(),0);

                gsVector<index_t> permP(m_Spatial_Bases[n][np].size());
                for(index_t i=0;i<m_Spatial_Bases[n][np].size();++i)
                    permP(i) = mapper.mapIndex(i);


                temp1=   m_TPData.PatchMx[n][np]->twistedBy(permP.asPermutation());
                temp1.uncompress();
                temp1.conservativeResize(mapper.freeSize(),mapper.freeSize());
                temp1.makeCompressed();
                m_TPData.PatchMx[n][np] = temp1.moveToPtr();
                temp1 =  m_TPData.PatchKx[n][np]->twistedBy(permP.asPermutation());
                temp1.uncompress();
                temp1.conservativeResize(mapper.freeSize(),mapper.freeSize()) ;
                m_TPData.PatchKx[n][np] =temp1.moveToPtr();
                temp1.makeCompressed();

                if(n<m_Temporal_Bases.nBases()-1)
                {
                    temp1 = ((*m_TPData.PatchMWx[n+1][np])*permP.asPermutation().inverse());
                    temp1.uncompress();
                    temp1.conservativeResize(m_TPData.PatchMWx[n+1][np]->rows(),mapper.freeSize());
                    temp1.makeCompressed();
                    m_TPData.PatchMWx[n+1][np] = temp1.moveToPtr();
                }
                if(n>0)
                {
                    temp1 = permP.asPermutation()*(*m_TPData.PatchMWx[n][np]);
                    temp1.uncompress();
                    temp1.conservativeResize(mapper.freeSize(),m_TPData.PatchMWx[n][np]->cols());
                    temp1.makeCompressed();
                    m_TPData.PatchMWx[n][np] = temp1.moveToPtr();
                }
                delete patchPde;
            }
        }

        sliceCount+=m_Temporal_Bases.basis(n).size()*mappers[n].size();

        timings(8)+= time.stop();
    }

    time.restart();
    gsMatrix<index_t> permFull(this->system().colMapper(0).size(),1);
    int iterF,iterB, currSlice,iSlice;
    iterF = iterB= currSlice = iSlice = 0;
    std::vector<std::pair<index_t, index_t> > result;
    for(int i=0; i<permFull.rows();++i)
    {
        int sIdx = iSlice % mappers[currSlice].size();
        int tIdx = (int)iSlice /mappers[currSlice].size();
        result = preImages[currSlice][sIdx];

        if(this->system().colMapper(0).is_free(result.front().second+tIdx*m_Spatial_Bases[currSlice].basis(result.front().first).size(),m_patchIndicesForTimeSlice[currSlice][result.front().first]))
            permFull(i) = iterF++;
        else
            permFull(i) = this->system().colMapper(0).freeSize()+iterB++;

        if((++iSlice)==mappers[currSlice].size()*m_Temporal_Bases[currSlice].size())
        {
            iSlice = 0;
            currSlice++;
        }
    }

    Permutation permutFull(permFull);
    rhs = permutFull*rhs; // This is the full vector!!

    rhs.conservativeResize(this->system().colMapper(0).freeSize(),1);
    m_TPData.rhs.resize(m_Temporal_Bases.nBases());

    gsVector<index_t> col(1); col<<1;
    typename gsMatrix<T>::BlockView view = rhs.blockView(m_sizes,col);
    for(size_t n = 0; n< m_Temporal_Bases.nBases(); ++n){
        m_TPData.rhs[n] = view(n,0);

        //   gsInfo<<"K on Slice: "<<n<<" in slice "<<n<<"\n"<<m_TPData.Kx[n]->toDense()<<"\n";
        //   gsInfo<<"M on Slice: "<<n<<" in slice "<<n<<"\n"<<m_TPData.Mx[n]->toDense()<<"\n";
        //   gsInfo<<"fx on Slice: "<<n<<" in slice "<<n<<"\n"<<m_TPData.rhs[n].transpose()<<"\n";
    }
    timings(9) = time.stop();
}


template<typename T>
void gsSpaceTimesliceAssembler<T>::assembleClassical(bool assemblePatchwise)
{
    gsStopwatch time;
    m_hasAssembledTP = false;
    gsAssemblerOptions assOpt;
    //typename gsBasis<T>::uPtr bBasis = m_bases[0][0].boundaryBasis(boundary::west);
    const index_t nz = 2*assOpt.numColNz(m_bases[0][0]); //+assOpt.numColNz(*bBasis);
    m_system.setZero();
    m_system.reserve(nz,m_pde_ptr->numRhs());

    gsSpaceTimePoissonPde<T>* STpde= dynamic_cast<gsSpaceTimePoissonPde<T>* >(m_pde_ptr.get());
    gsPoissonHeterogeneousPde<T> spatialPde(m_Spatial_Domain,m_Spatial_boundary_conditions,*STpde->rhs_x(),*STpde->getAlpha());

    std::vector<typename gsSparseMatrix<T>::uPtr> Ks(m_Temporal_Bases.nBases()), Mxs(m_Temporal_Bases.nBases());
    std::vector<typename gsMatrix<T>::uPtr> fxs(m_Temporal_Bases.nBases());


    gsOptionList opt = gsAssembler<T>::defaultOptions();
    opt.setInt("DirichletStrategy" , dirichlet::none );
    opt.setInt("InterfaceStrategy", iFace::none  );

    for(size_t n = 0; n< m_Temporal_Bases.nBases(); ++n)
    {
        real_t h = m_tau[n];
        //Calculate the temporal matrices
        time.restart();
        gsGenericAssembler<T> genTemporal(m_Temporal_Domain.patch(n),m_Temporal_Bases.basis(n),opt);
        gsSparseMatrix<T> KtpMt = genTemporal.assembleStiffness();
        const gsSparseMatrix<T>& Mt = genTemporal.assembleMass();
        timings(1)+=time.stop();

        KtpMt*=m_theta*h;

        gsExprAssembler<T> dtuvAssembler(1,1);
        gsBoundaryConditions<T> bc;


        typedef typename gsExprAssembler<T>::variable variable;
        typedef typename gsExprAssembler<T>::space    space;
        typedef typename gsExprAssembler<T>::geometryMap geoMap;

        time.restart();
        gsMultiBasis<T> basis(m_Temporal_Bases.basis(n));
        geoMap G = dtuvAssembler.getMap(m_Temporal_Domain.patch(n));
        space u = dtuvAssembler.getSpace(basis, 1, 0); // args: splines, bc, dim, id
        u.addBc( bc.get("Dirichlet",0) );
        variable ft = dtuvAssembler.getCoeff(*STpde->rhs_t(),G); // force
        dtuvAssembler.setOptions(opt);

        dtuvAssembler.setIntegrationElements(basis); // subdivided mesh
        dtuvAssembler.initSystem();
        dtuvAssembler.assemble(u * igrad(u,G).tr()*meas(G),(u+m_theta*h*igrad(u,G))*ft*meas(G)); //
        timings(2)+=time.stop();

        const gsSparseMatrix<T>& Tt = dtuvAssembler.matrix().transpose();
        KtpMt+=Tt;
        const gsMatrix<T>& f_t = dtuvAssembler.rhs();

        //gsInfo<<"T on slice "<<n<<"\n"<<Tt.toDense()<<"\n";
        //gsInfo<<"M on slice "<<n<<"\n"<<Mt.toDense()<<"\n";
        //gsInfo<<"f_t on slice "<<n<<"\n"<<f_t<<"\n";

        //calculate the spatial matrices
        for(size_t np = 0; np<m_Spatial_Bases[n].nBases();++np)
        {
            const gsPoissonHeterogeneousPde<T>& ppde = *dynamic_cast<gsPoissonHeterogeneousPde<T>*>((spatialPde.restrictToPatch(np)));
            //TODO: implement the reusing of already computed spatial matrices!
            gsPoissonHeterogeneousAssembler<T> spatialAssembler(ppde,m_Spatial_Bases[n].basis(np),dirichlet::none,iFace::none);
            spatialAssembler.options().setInt("DirichletValues",dirichlet::homogeneous); //just something, that the boudary condition is not evaluated.

            time.restart();
            spatialAssembler.assemble();
            timings(3)+=time.stop();


            gsBoundaryConditions<T> bc_patch;
            m_Spatial_boundary_conditions.getConditionsForPatch(np,bc_patch);
            gsGenericAssembler<T> genSpatial(m_Spatial_Domain.patch(np),m_Spatial_Bases[n].basis(np),opt,&bc_patch);
            time.restart();
            genSpatial.assembleMass();
            timings(4)+=time.stop();


            Ks[n] = spatialAssembler.system().matrix().moveToPtr();
            fxs[n] = spatialAssembler.system().rhs().moveToPtr();
            Mxs[n] = genSpatial.system().matrix().moveToPtr();

            //gsInfo<<"K on Patch: "<<np<<" in slice "<<n<<"\n"<<Ks[n]->toDense()<<"\n";
            //gsInfo<<"M on Patch: "<<np<<" in slice "<<n<<"\n"<<Mxs[n]->toDense()<<"\n";
            //gsInfo<<"fx on Patch: "<<np<<" in slice "<<n<<"\n"<<*fxs[n]<<"\n";

            //Build the Space-Time matrices for each timeslice
            time.restart();
            gsSparseMatrix<T> KMt,MxTt,KTt, result;
            KMt=Mt.kron(*Ks[n]);
            MxTt=KtpMt.kron(*Mxs[n]);
            KTt=gsSparseMatrix<T>(Tt.transpose()).kron(*Ks[n]);

            //Build the right hand side;
            gsMatrix<T> F;
            F=f_t.kron(*fxs[n]);
            timings(5)+=time.stop();

            //gsInfo<<"KMt on Patch: "<<np<<" in slice "<<n<<"\n"<<KMt.toDense()<<"\n";
            //gsInfo<<"MxTt on Patch: "<<np<<" in slice "<<n<<"\n"<<MxTt.toDense()<<"\n";
            //gsInfo<<"KTt on Patch: "<<np<<" in slice "<<n<<"\n"<<KTt.toDense()<<"\n";
            //gsInfo<<"F on Patch: "<<np<<" in slice "<<n<<"\n"<<F<<"\n";

            result = KMt + MxTt +m_theta*h*KTt;


            index_t sz = m_bases.front().basis(m_patchIndicesForTimeSlice[n][np]).size();
            gsMatrix<index_t> actives;
            m_system.mapColIndices(gsVector<index_t>::LinSpaced(sz,0,(sz-1)),m_patchIndicesForTimeSlice[n][np],actives);

            //gsInfo<<"Matrix on Patch: "<<np<<" in slice "<<n<<"\n"<<result.toDense()<<"\n";
            //gsInfo<<"Rhs on Patch: "<<np<<" in slice "<<n<<"\n"<<F<<"\n";
            time.restart();
            //m_system.push(result,F,actives,m_ddof.front());
            m_system.pushSparse(result.transpose(),F,actives,actives,m_ddof.front());
            timings(6)+=time.stop();

            //
            //gsInfo<<"System after Patch: "<<np<<" in slice "<<n<<"\n"<<m_system.matrix().toDense().transpose()<<"\n"<<std::flush;
            //gsInfo<<"Rhs after Patch: "<<np<<" in slice "<<n<<"\n"<<m_system.rhs()<<"\n"<<std::flush;

            delete &ppde;
        }
    }

    time.restart();
    gsAssembler<T>::template push<gsVisitorSpaceTimeInitial<T> >(m_pde_ptr->bc().container("Initial"));
    timings(7)+=time.stop();

    time.restart();
    gsVisitorSpaceTimeInterface<T> visitor(*m_pde_ptr);
    for ( std::vector<boundaryInterface>::iterator it=m_dGInterfaces.begin(); it!=m_dGInterfaces.end();++it)
    {
        const boundaryInterface & iFace = //recover master element
                ( m_bases[0][it->first() .patch].numElements(it->first() .side() ) <
                m_bases[0][it->second().patch].numElements(it->second().side() ) ?
                    it->getInverse() : *it );

        this->apply(visitor, iFace);
    }
    timings(8)+=time.stop();

    Base::finalize();

}


template<typename T>
void gsSpaceTimesliceAssembler<T>::printTimings()
{
    if(m_hasAssembledTP)
    {
        gsInfo<<"Time required to compute DDof values: "<<timings(0)<<"\n";
        gsInfo<<"Time required to assemble Kt and Mt: "<<timings(1)<<"\n";
        gsInfo<<"Time required to assemble Time-T and Rhs: "<<timings(2)<<"\n";
        gsInfo<<"Time required to assemble Spatial-K and Rhs: "<<timings(3)<<"\n";
        gsInfo<<"Time required to assemble Spatial-Mass: "<<timings(4)<<"\n";
        gsInfo<<"Time required to assemble initial condition: "<<timings(5)<<"\n";
        gsInfo<<"Time required to assemble interface: "<<timings(6)<<"\n";
        gsInfo<<"Time required to update rhs with BC: "<<timings(7)<<"\n";
        gsInfo<<"Time required to reorder/shrink spatial matrices: "<<timings(8)<<"\n";
        gsInfo<<"Time required to reorder/shrink rhs: "<<timings(9)<<"\n";
    }
    else
    {
        gsInfo<<"Time required to compute DDof values: "<<timings(0)<<"\n";
        gsInfo<<"Time required to assemble Time-Masses: "<<timings(1)<<"\n";
        gsInfo<<"Time required to assemble Time-T and Rhs: "<<timings(2)<<"\n";
        gsInfo<<"Time required to assemble Spatial-K and Rhs: "<<timings(3)<<"\n";
        gsInfo<<"Time required to assemble Spatial-Mass: "<<timings(4)<<"\n";
        gsInfo<<"Time required to build Konecker: "<<timings(5)<<"\n";
        gsInfo<<"Time required to push: "<<timings(6)<<"\n";
        gsInfo<<"Time required to assemble initial: "<<timings(7)<<"\n";
        gsInfo<<"Time required to assemble dG: "<<timings(8)<<"\n";
    }
}


template<typename T>
gsMatrix<T> gsSpaceTimesliceAssembler<T>::assembleTemporalMassForTimeSlice(size_t n)
{
    gsOptionList opt = gsAssembler<T>::defaultOptions();
    gsBoundaryConditions<T> bc;
    gsGenericAssembler<T> genAss(m_Temporal_Domain.patch(n),m_Temporal_Bases[n],opt,&bc);
    genAss.assembleMass();
    return genAss.system().matrix().toDense();
}


/// computes the restrictions from two timeslices (a,b)(b,c) to one timeslice (a,c)
/// u_coarse[n/2] = transfer1*u_fine[n] + transfer2*u_fine[n+1]
template<typename T>
void gsSpaceTimesliceAssembler<T>::computeTransferForTwoSlices(size_t n, gsMatrix<T>& transfer1, gsMatrix<T>& transfer2)
{
    bool singleTimeSlap = m_Temporal_Bases.nBases()==1;
    int n1 = m_Temporal_Bases.basis(n).size();
    int n2 = !singleTimeSlap? m_Temporal_Bases.basis(n+1).size() : n1;
    GISMO_ASSERT(n1%n2 == 0 || n2%n1 == 0, "The two basis do not have nested spaces.");

    //Take the finer of the two basis as coarse basis
    gsBasis<T>& finerBasis = n1>=n2?m_Temporal_Bases.basis(n) :m_Temporal_Bases.basis(n+1);
    gsMultiBasis<T> finerBasisMP(finerBasis);


    typename gsBSplineBasis<T>::uPtr basisC_1 = memory::convert_ptr<gsBSplineBasis<T> >(finerBasis.clone());
    typename gsBSplineBasis<T>::KnotVectorType & knot1 = basisC_1->knots(0);
    knot1.affineTransformTo(0,2);
    gsMultiBasis<T> basisC1(*basisC_1);

    typename gsBSplineBasis<T>::uPtr basisC_2 = memory::convert_ptr<gsBSplineBasis<T> >(finerBasis.clone());
    typename gsBSplineBasis<T>::KnotVectorType & knot2 = basisC_2->knots(0);
    knot2.affineTransformTo(-1,1);
    gsMultiBasis<T> basisC2(*basisC_2);

    gsBoundaryConditions<T> bc;
    typedef typename gsExprAssembler<T>::variable variable;
    typedef typename gsExprAssembler<T>::space    space;
    typedef typename gsExprAssembler<T>::geometryMap geoMap;


    gsExprAssembler<T> massAssembler1(2,2);

    geoMap G1 = massAssembler1.getMap(m_Temporal_Domain[n]);
    gsMultiBasis<T> basis1(m_Temporal_Bases.basis(n));

    space uF1 = massAssembler1.getSpace(basis1, 1, 0); // args: splines, dim, id
    uF1.addBc( bc.get("Dirichlet",0) );
    space uC1 = massAssembler1.getSpace(basisC1, 1, 1); // args: splines, dim, id
    uC1.addBc( bc.get("Dirichlet",1) );

    massAssembler1.setIntegrationElements(finerBasisMP); // subdivided mesh
    massAssembler1.initSystem();
    massAssembler1.assemble(uF1 * uC1.tr()*meas(G1));


    gsExprAssembler<T> massAssembler2(2,2);

    geoMap G2 = massAssembler2.getMap(!singleTimeSlap ? m_Temporal_Domain[n+1] :m_Temporal_Domain[n] );
    gsMultiBasis<T> basis2(!singleTimeSlap ? m_Temporal_Bases.basis(n+1):m_Temporal_Bases.basis(n) );

    space uF2 = massAssembler2.getSpace(basis2, 1, 0); // args: splines, dim, id
    uF2.addBc( bc.get("Dirichlet",0) );
    space uC2 = massAssembler2.getSpace(basisC2, 1, 1); // args: splines, dim, id
    uC2.addBc( bc.get("Dirichlet",1) );

    massAssembler2.setIntegrationElements(finerBasisMP); // subdivided mesh
    massAssembler2.initSystem();
    massAssembler2.assemble(uF2 * uC2.tr()*meas(G2));

    gsMatrix<T> tr1 = massAssembler1.matrix().block(0,basis1[0].size(),basisC1[0].size(),basisC1[0].size()).toDense();
    gsMatrix<T> tr2 = massAssembler2.matrix().block(0,basis2[0].size(),basisC2[0].size(),basisC2[0].size()).toDense();

    gsMatrix<T> M1 = assembleTemporalMassForTimeSlice(n);
    gsMatrix<T> M2 = assembleTemporalMassForTimeSlice(!singleTimeSlap ? n+1 : n);

    transfer1 = (M1.inverse()*tr1).transpose();
    transfer2 = (M2.inverse()*tr2).transpose();

    // gsInfo<<"\n"<<transfer1<<"\n\n";
    //  gsInfo<<"\n"<<transfer2<<"\n\n"<<std::flush;

}

template<typename T>
typename gsSpaceTimesliceAssembler<T>::TPData gsSpaceTimesliceAssembler<T>::getTensorProductData(std::vector<index_t>& slices)
{
    GISMO_ASSERT(m_hasAssembledTP, "you have not assembled the TP matrices.");
    TPData result;
    result.rhs.resize(slices.size());
    //  result.K1t.resize(slices.size()); //always give the whole time stuff
    //  result.Kht.resize(slices.size()); //always give the whole time stuff
    //  result.Kt.resize(slices.size());  //always give the whole time stuff
    result.Kx.resize(slices.size());
    if(m_TPData.PatchKx.size() != 0)result.PatchKx.resize(slices.size());
    //  result.M1t.resize(slices.size()); //always give the whole time stuff
    //  result.Mht.resize(slices.size()); //always give the whole time stuff
    //  result.Mt.resize(slices.size());  //always give the whole time stuff
    result.Mx.resize(slices.size());
    if(m_TPData.PatchMx.size() != 0)result.PatchMx.resize(slices.size());
    if(m_TPData.PatchMWx.size() != 0)result.PatchMWx.resize(slices.size());
    //  result.Nt.resize(slices.size());  //always give the whole time stuff
    result.MWx.resize(slices.size());

    //result.K1t = m_TPData.K1t;
    //result.Kht = m_TPData.Kht;
    result.Kt = m_TPData.Kt;
    // result.M1t = m_TPData.M1t;
    //result.Mht = m_TPData.Kht;
    result.Mt = m_TPData.Mt;
    result.Nt = m_TPData.Nt;
    for(size_t i=0; i<slices.size();++i)
    {
        result.rhs[i]=(m_TPData.rhs[slices[i]]);
        //  result.K1t[i]=(m_TPData.K1t[slices[i]]);
        //  result.Kht[i]=(m_TPData.Kht[slices[i]]);
        //   result.Kt[i]=(m_TPData.Kt[slices[i]]);
        result.Kx[i]=(m_TPData.Kx[slices[i]]);
        //    result.M1t[i]=(m_TPData.M1t[slices[i]]);
        //    result.Mht[i]=(m_TPData.Mht[slices[i]]);
        //    result.Mt[i]=(m_TPData.Mt[slices[i]]);
        result.Mx[i]=(m_TPData.Mx[slices[i]]);

        if(m_TPData.PatchKx.size() != 0)result.PatchKx[i]=(m_TPData.PatchKx[slices[i]]);
        if(m_TPData.PatchMx.size() != 0)result.PatchMx[i]=(m_TPData.PatchMx[slices[i]]);
        if(slices[i]>0 && m_TPData.PatchMWx.size() != 0)result.PatchMWx[i]=(m_TPData.PatchMWx[slices[i]]);
        if(slices[i]>0) result.MWx[i]=(m_TPData.MWx[slices[i]]);
    }
    return result;
}


} // namespace gismo
