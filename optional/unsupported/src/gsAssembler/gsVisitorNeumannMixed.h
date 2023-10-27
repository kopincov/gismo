
#pragma once


#include <gsCore/gsVectorValuedFunctionSet.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsCore/gsTransformPure.h>
#include <gsCore/gsDebug.h>

namespace gismo
{

/**
Class for the stokes based on gsAssemblerBase2
**/
GISMO_DEPRECATED
template <class T>
        class gsVisitorNeumannMixed : public gsVisitorNeumann<T>
{
public:
    typedef gsVisitorNeumann<T> Base;
public:
    GISMO_DEPRECATED
    gsVisitorNeumannMixed(const gsFunction<T> & neudata, boundary::side s) :
        Base(neudata,s)
    { }
    void localToGlobal(const gsDofMapper     & mapper,
                       const gsMatrix<T>     & eliminatedDofs,
                       const index_t           patchIndex,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        const index_t usz = mapper.freeSize();
        // Local Dofs to global dofs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();
        // Push element contribution to the global load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            // convert local dof index to global dof index
            const unsigned jj = actives(j);
            if (mapper.is_free_index(jj))
                for (index_t s = 0; s!=localRhs.cols(); ++s)
                    rhsMatrix(jj + s*usz,0) += localRhs(j,s);
        }
    }
private:
    // Neumann function
    using Base::neudata_ptr;
    using Base::side;

    // Basis values
    using Base::basisData;
    using Base::actives;

    // Normal and Neumann values
    using Base::unormal;
    using Base::neuData;

    // Local matrix and rhs
    using Base::localMat;
    using Base::localRhs;
};

/** @brief
    Implementation of a Neumann BC for the Stokes Equation. Covers
    also the Divconforming case. Note: it is required, that each
    boundary condition stores only a scalar function. The algorithm
    will crash if not so.

    However, you are allowed to loop only once around the boundary.
    The visitor will look by itself for the remaining components.

    It adds the following term to the linear term.
    \f[ \nabla u \cdot \mathbf{n} = g_N  \f]
*/
template <class T>
class gsVisitorStokesNeumann
{
public:

    gsVisitorStokesNeumann(const gsPde<T> & pde, const boundary_condition<T> & s)
        : side(s.side())
    {

        //Here we assume that all components of v on this side are Neumann dofs
        neuData_ptr.resize(pde.unknownDim()(0));

        //We look for the remaining components of the Neumann BC for this unknown
        //and store them
        for ( typename gsBoundaryConditions<T>::const_iterator
              it = pde.bc().neumannBegin();
              it != pde.bc().neumannEnd(); ++it )
        {
            if((*it).ps == s.ps)
                neuData_ptr[(*it).unknown()]= (*it).function().get();
        }
        for(size_t i=0; i<neuData_ptr.size();++i )
            GISMO_ASSERT(neuData_ptr[i]->targetDim()==1, "The Algorithm cannot handle vector valued BC, each component should be and individual BC ");
    }


    /// Initialize the Visitor
    void initialize(const gsBasisRefs<T> & basis,
                    const index_t patchIndex,
                    const gsAssemblerOptions & options,
                    gsQuadRule<T>    & rule,
                    unsigned         & evFlagsMap )
    {

        //Store the options
        m_opt = options;
        d = m_opt.info.unknownDims(0);

        // Setup Quadrature (harmless slicing occurs)
        rule = gsGaussRule<T>(basis, options.quA, options.quB, side.direction() );

        // Set Basis evaluation flag
        basisDataV.addFlags( SAME_ELEMENT| NEED_ACTIVE | NEED_VALUE);

        // Set Geometry evaluation flags
        evFlagsMap |= NEED_VALUE | NEED_MEASURE|NEED_OUTER_NORMAL;

        //Choose if V is divconforming or not and add corresponding flags
        if(m_opt.info.getTransform(0)==transform::Hdiv)
            evFlagsMap |= gsTransformPureDivConforming<T>::getGeoFlags(basisDataV.flags);
        else
            evFlagsMap |= gsTransformPureGradConforming<T>::getGeoFlags(basisDataV.flags);

        neuData.addFlags(SAME_ELEMENT|NEED_VALUE);

        //We build our vector basis
        std::vector<gsBasis<T>* > v_bases(d);
        //The const cast is required in order to extract certein components of basisRefs and
        //store them in a std::vector to call the gsBasis. It would be nice to have a direct
        //pipeline from gsBasisRefs to gsVBasis, where the used components are specified.
        for(int comp_v = 0; comp_v < d; ++comp_v) // all components for unknown 0
            v_bases[comp_v]=const_cast<gsBasis<T>*>(&basis[m_opt.info.getBasis(0,comp_v)]);

        //Set shifts to zero, there are stored in gsSparseSystem
        std::vector<index_t> shifts(d,0);
        vBasis = memory::make_unique(new gsVBasis<T>(v_bases,shifts));
    }

    /// Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
                         const gsMapData<T> & mapData)
    {
        //dependent on the transformation, eval the basis and transform
        if(m_opt.info.getTransform(0)==transform::Hdiv)
            gsTransformPureDivConforming<T>::transform(mapData,*vBasis,basisDataV);
        else
            gsTransformPureGradConforming<T>::transform(mapData,*vBasis,basisDataV);

        //The vector valued Neumann values
        neuVals.resize(d,mapData.points.cols());

        // Evaluate the Neumann data
        for(index_t i=0; i<d;++i)
        {
            neuData_ptr[i]->compute(mapData.values[0], neuData);
            neuVals.row(i)=neuData.values[0];
        }

        // Initialize local matrix/rhs
        localRhs_u.setZero(basisDataV.actives.rows(), 1);

    }

    /// assemble the local contributions
    inline void assemble(gsDomainIterator<T>    & element,
                         const gsMapData<T>     & mapData,
                         const gsVector<T>      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const gsMatrix<T> & vVals = basisDataV.eval(k); // d x n_actV
            const gsMatrix<T> & neuVal = neuVals.col(k); // d x 1

            // Compute the outer normal vector on the side
            const T weight =  quWeights[k] * mapData.outNormal(k).norm();

            // Sum up quadrature point evaluations
            localRhs_u.noalias() += weight * (vVals.transpose()*neuVal);
        }
    }

    /// push the local contributions to the sparse system
    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        //Make vectors for block structure
        gsVector<index_t> dBlocks(d),oneCol(1);
        for(int i =0; i<d;++i)
            dBlocks(i) = getActiveBlockLength(i);
        oneCol<<1;

        //deduce block structure
        typename gsMatrix<index_t>::BlockView actives_v_Block = basisDataV.actives.blockView(dBlocks, oneCol);
        typename gsMatrix<T>::BlockView rhs_Block = localRhs_u.blockView(dBlocks,oneCol);

        // Map patch-local DoFs to global DoFs, a blockview cannot be overridden (idk why)
        std::vector<gsMatrix<index_t> > actives_v(d);
        for(int i=0; i!=d;++i)
            system.mapColIndices(actives_v_Block(i,0), patchIndex,actives_v[i],i);

        // Add contributions to the system matrix and right-hand side
        for(index_t i=0; i!=d;++i)
            system.pushToRhs( rhs_Block(i,0),actives_v[i],i);

    }
protected:
    /// This auxiliary function calculates the blocklength of the actives for each component
    /// Might be different for each component.
    index_t getActiveBlockLength(index_t component)
    {
        index_t i=1, start=0;
        component--;
        for(i=1; i<basisDataV.actives.rows();++i)
            if(basisDataV.actives(i,0)<=basisDataV.actives(i-1,0))
            {
                component--;
                if(component == -1)
                    start = i;
                else if(component== -2)
                    return i-start;
            }
        //End of actives reached (i is already increased)
        return i-start;
    }

private:
    gsAssemblerOptions m_opt;

    // Velocity vector dimension
    index_t d;

    // Side
    boxSide side;

    // Neumann function
    std::vector<gsFunction<T> *> neuData_ptr;

    //velocity vector basis
    memory::unique_ptr<gsVBasis<T> > vBasis;

private:

    //Basis and function Data
    gsFuncData<T> basisDataV;
    gsFuncData<T> neuData;
    gsMatrix<T> neuVals;


protected:
    // Local values of the right hand side
    gsMatrix<T> localRhs_u;
};





} // namespace gismo
