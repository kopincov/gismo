/** @file gsSpaceTimesliceAssembler.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer, M. Neumüller
    Created on: 2017-06-06
*/


#pragma once

#include <gsIO/gsOptionList.h>
#include <gsPde/gsSpaceTimePoissonPde.h>
#include <gsAssembler/gsAssembler.h>
#include <gsSolver/gsKronecker.h>
#include <gsAssembler/gsVisitorSpaceTimeSlice.h>
#include <gsUtils/gsStopwatch.h>
#include <gsNurbs/gsTensorBSplineBasis.h>



namespace gismo {


namespace internal{

template<int d, typename T> struct Slicer {
    static void reduceSpline(gsTensorBSpline<d,T> & inp, gsBSpline<T>& result)
    {
        typename gsTensorBSpline<d,T>::BoundaryGeometryType tb_new;
        inp.slice(d-1,0,tb_new);
        Slicer<d-1,T>::reduceSpline(tb_new,result);
    }
};

template<typename T>
struct Slicer<2,T>
{
    static void reduceSpline(gsTensorBSpline<2,T> & inp, gsBSpline<T>& result)
    {
        inp.slice(1,0,result);
    }
};

}

template<class T>
class gsSpaceTimesliceAssembler : public gsAssembler<T>
{
    typedef gsAssembler<T> Base;

    template<int d, typename > friend struct internal::Slicer;
    /// Default Constructor
    gsSpaceTimesliceAssembler(real_t theta): m_theta(theta) {}
public:

    struct TPData
    {
        std::vector<typename gsSparseMatrix<T>::Ptr > Mx, Mt, Kt, Kx, MWx, Nt;
        std::vector<gsMatrix<T> > rhs;
        std::vector<std::vector<typename gsSparseMatrix<T>::Ptr> > PatchMx, PatchKx, PatchMWx;
    };


    typedef Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic, index_t> Permutation;

    gsSpaceTimesliceAssembler( const gsSpaceTimePoissonPde<T> & pde,
                               const gsMultiBasis<T>          & bases,
                               real_t            theta,
                               dirichlet::strategy           dirStrategy,
                               iFace::strategy               intStrategy = iFace::glue,
                               bool timeIsLastDirection = false
            ): m_timeIsLastDirection(timeIsLastDirection), m_theta(theta)
    {

        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        Base::initialize(pde, bases, m_options);
        this->system().rhs().conservativeResize(Eigen::NoChange, 1); // adaption to commit #267

    }

    virtual gsAssembler<T>* clone() const
    {
        return new gsSpaceTimesliceAssembler<T>(*this);
    }

    virtual gsAssembler<T>* create() const
    {
        return new gsSpaceTimesliceAssembler<T>(m_theta);
    }


    ///Extract the part of the geometry \a inp which points in \a desiredDirection
    void reduceSpline(gsGeometry<T> & inp, int desiredDirection, gsBSpline<T>& result);



    virtual void refresh();
    virtual void assemble() {assemble(false,false);}
    virtual void assemble(bool assembleTP, bool assemblePatchwise=false)
    {
        GISMO_ASSERT(m_system.initialized(),
                     "Sparse system is not initialized, call initialize() or refresh()");

        timings.setZero(10);

        gsStopwatch time;
        Base::computeDirichletDofs(0);
        timings(0)+=time.stop();

        if(assembleTP)
            assembleTPMats(assemblePatchwise);
        else
            assembleClassical(assemblePatchwise);
    }

    virtual void assembleTPMats(bool assemblePatchwise);
    virtual void assembleClassical(bool assemblePatchwise);

    void printTimings();

    gsMatrix<T> assembleTemporalMassForTimeSlice(size_t n);

    /// computes the restrictions from two timeslices (a,b)(b,c) to one timeslice (a,c)
    /// u_coarse[n/2] = transfer1*u_fine[n] + transfer2*u_fine[n+1]
    void computeTransferForTwoSlices(size_t n, gsMatrix<T>& transfer1, gsMatrix<T>& transfer2);

    real_t getTheta() {return m_theta;}
    real_t getTau(int slice) {return m_tau[slice];}
    real_t getH(int slice) {return m_h[slice];}

    std::vector<boundaryInterface>& getdGInterfaces() {return m_dGInterfaces;}
    std::vector<patchSide>&  getInitialBoundary() {return m_initials;}
    std::vector<patchSide>&  getTopBoundary() {return m_tops;}
    const gsVector<index_t>& getBlockViewSizes() {return m_sizes;}
    size_t getNTimeSlices() {return m_Temporal_Bases.nBases();}
    size_t getNSpacePatches() {return m_Spatial_Bases.front().nBases();}

    const gsMultiBasis<T>& getTimeBases() const {return m_Temporal_Bases;}
    const gsMultiBasis<T>& getSpaceBases(int slice) const {return m_Spatial_Bases[slice];}
    typename gsMultiPatch<T>::Ptr getSlice(int slice) const
    {
        // std::pair<gsMultiBasis<T>,gsSpaceTimePoissonPde<T>::Ptr > res;
        typename gsMultiPatch<T>::Ptr mp = memory::make_shared(new gsMultiPatch<T>());
        for(size_t i=0; i<m_patchIndicesForTimeSlice[slice].size();++i)
        {
            mp->addPatch(m_pde_ptr->domain().patch(m_patchIndicesForTimeSlice[slice][i]));

            //   res.first.addBasis(m_bases[m_patchIndicesForTimeSlice[slice][i]]);
        }
        //   gsSpaceTimePoissonPde<T>* STpde= dynamic_cast<gsSpaceTimePoissonPde<T>* >(m_pde_ptr.get());
        mp->computeTopology();
        //   res.first.setTopology(mp.topology());
        //  res.second=memory::make_unique(
        //              new gsSpaceTimePoissonPde<T>(mp,m_Spatial_boundary_conditions,STpde->rhs(),STpde->rhs_x(),STpde->rhs_t(),*STpde->getAlpha()));
        return mp;
    }

    memory::unique_ptr<gsPoissonHeterogeneousPde<T> > getSpatialPde() const
    {
        const gsSpaceTimePoissonPde<T>& stpde= *dynamic_cast<const gsSpaceTimePoissonPde<T>* >(m_pde_ptr.get());
        return memory::make_unique(
                    new gsPoissonHeterogeneousPde<T>(m_Spatial_Domain,m_Spatial_boundary_conditions,*stpde.rhs_x(),*stpde.getAlpha())
                    );
    }

    const gsBoundaryConditions<T>& getSpaceBC() const {return m_Spatial_boundary_conditions;}

    TPData& getTensorProductData() {GISMO_ASSERT(m_hasAssembledTP, "you have not assembled the TP matrices."); return m_TPData;}
    TPData getTensorProductData(std::vector<index_t>& slices);

    // select which patch should be assembled (patches[slice][np] = true/false, slice = 1,..., nSlices, np = 1,...,nSpatialPatches
    void setPatchesForAssembling(std::vector<std::vector<bool> >& patches) {m_patchesForAssembling = patches;}

protected:
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

    bool m_timeIsLastDirection;
    real_t m_theta;
    std::vector<real_t> m_h, m_tau;

    gsMultiPatch<T> m_Spatial_Domain;
    gsMultiPatch<T> m_Temporal_Domain;
    std::vector<gsMultiBasis<T> > m_Spatial_Bases;
    gsMultiBasis<T> m_Temporal_Bases;
    gsBoundaryConditions<T> m_Spatial_boundary_conditions;

    std::vector<std::vector<index_t> >m_patchIndicesForTimeSlice;
    std::vector<index_t>m_TimeSliceForPatch;

    std::vector<boundaryInterface> m_dGInterfaces;

    std::vector<patchSide> m_initials, m_tops;

    gsVector<real_t> timings;
    gsVector<index_t> m_sizes;

    bool m_hasAssembledTP;
    std::vector<std::vector<bool> > m_patchesForAssembling;

    TPData m_TPData;

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSpaceTimesliceAssembler.hpp)
#endif
