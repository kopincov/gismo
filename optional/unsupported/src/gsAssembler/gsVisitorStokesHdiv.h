
#include<gsPde/gsStokesPde.h>
#include <gsCore/gsVectorValuedFunctionSet.h>

#pragma once



/** @brief Visitor for the Stokes problem using Div conforming elements

   ( A11   A12  A13   B1^T )      (rhs_u1)
   ( A21   A22  A23   B2^T )  and (rhs_u2)
   ( A31   A32  A33   B3^T )      (rhs_u3)
   ( B1  B2  B3   0        )      (rhs_p )

   Note that the rhs u must be a vector valued function.
*/
namespace gismo
{

template <class T>
class gsVisitorStokesHdiv
{
public:

    /**
     * @brief gsVisitorStokesHdiv Constructor for the stokes visitor
     * @param pde the Stokes PDE
     */
    gsVisitorStokesHdiv(const gsPde<T> & pde)
    {
        const gsStokesPde<T>* spde = static_cast<const gsStokesPde<T>* >(&pde);
        rhs_ptr =spde->rhs() ;
        m_viscosity =spde->viscocity();

    }

    /// Initialize the Visitor
    void initialize(const gsBasisRefs<T> & basis,
                    const index_t patchIndex,
                    const gsAssemblerOptions & options,
                    gsQuadRule<T>    & rule,
                    unsigned         & evFlagsMap
                    )
    {
        //Store the options
        m_info = options.info;
        d = m_info.unknownDims(0); // number of velocity unknowns

        // Setup Quadrature
        rule = gsGaussRule<T>(basis[0], options.quA, options.quB);// harmless slicing occurs here


        // Set Basis evaluation flag
        unsigned evFlagsFunc = SAME_ELEMENT| NEED_ACTIVE | NEED_VALUE;
        basisDataV.addFlags(evFlagsFunc|NEED_DERIV|NEED_DIV);
        basisDataP.addFlags(evFlagsFunc);

        // Set Geometry evaluation flags
        evFlagsMap |= NEED_VALUE | NEED_MEASURE;
        evFlagsMap |= gsTransformPureDivConforming<T>::getGeoFlags(basisDataV.flags);
        evFlagsMap |= gsTransformPureGradConforming<T>::getGeoFlags(basisDataP.flags);

        // Set Function evaluation flag
        funcData.addFlags(SAME_ELEMENT|NEED_VALUE);

        //We build our vector basis
        std::vector<gsBasis<T>* > v_bases(d);
        //The const cast is required in order to extract certein components of basisRefs and
        //store them in a std::vector to call the gsBasis. It would be nice to have a direct
        //pipeline from gsBasisRefs to gsVBasis, where the used components are specified.
        for(int comp_v = 0; comp_v < d; ++comp_v) // all components for unknown 0
            v_bases[comp_v]=const_cast<gsBasis<T>*>(&basis[m_info.getBasis(0,comp_v)]);

        //Set shifts to zero, there are stored in gsSparseSystem
        std::vector<index_t> shifts(d,0);
        vBasis = memory::make_unique(new gsVBasis<T>(v_bases,shifts));
    }


    /// Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
                         const gsMapData<T> & mapData)
    {

        //extract the presure basis, unknown 1, its zero component (the only one).
        const gsBasis<T>& pBasis = basisRefs[m_info.getBasis(1,0)];

        //Do Piola (we use the vBasis form the initialize
        gsTransformPureDivConforming<T>::transform(mapData,*vBasis,basisDataV);
        //Do transformation with inverse geometrical map
        gsTransformPureGradConforming<T>::transform(mapData,pBasis,basisDataP);

        unsigned numActV = basisDataV.actives.rows();
        unsigned numActP = basisDataP.actives.rows();

        // Evaluate right-hand side at the geometry points
        rhs_ptr->compute(mapData.values[0], funcData );

        // Initialize local matrix/rhs
        localMatA .setZero(numActV, numActV);
        localRhs_u.setZero(numActV, 1);
        localMatB.setZero(numActP, numActV);
    }

    /// assemble the local contributions
    inline void assemble(gsDomainIterator<T>    & element,
                         const gsMapData<T>     & mapData,
                         gsVector<T> const      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const gsMatrix<T> & bGradsV= basisDataV.deriv(k); //d*d x n_actV
            const gsMatrix<T> & bDivV= basisDataV.div(k);     // 1 x n_actV
            const gsMatrix<T> & bValsV = basisDataV.eval(k);  // d x n_actV
            const gsMatrix<T> & bValsP = basisDataP.eval(k);  // 1 x n_actP
            const gsMatrix<T> & rhsVals = funcData.eval(k);   // d x 1

            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * mapData.measure(k);

            // Right-hand side
            localRhs_u +=  weight * (bValsV.transpose() *  rhsVals);

            // Local block A
            localMatA.noalias()  += weight * m_viscosity * (bGradsV.transpose() * bGradsV);

            // Local blocks B_i
            localMatB.noalias() += weight * ( bValsP.transpose() * bDivV );
        }
    }

    /// push the local contributions to the sparse system
    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        //Make vectors for block structure
        gsVector<index_t> dBlocks(d),oneBlock(1),oneCol(1);
        for(int i =0; i<d;++i)
            dBlocks(i) = getActiveBlockLength(i);
        oneBlock<<basisDataP.actives.rows();
        oneCol<<1;

        //deduce block structure
        typename gsMatrix<index_t>::BlockView actives_v_Block = basisDataV.actives.blockView(dBlocks, oneCol);
        typename gsMatrix<T>::BlockView A_Block = localMatA.blockView(dBlocks,dBlocks);
        typename gsMatrix<T>::BlockView B_Block = localMatB.blockView(oneBlock,dBlocks);
        typename gsMatrix<T>::BlockView rhs_Block = localRhs_u.blockView(dBlocks,oneCol);

        // Map patch-local DoFs to global DoFs,  a blockview cannot be overridden (idk why)
        std::vector<gsMatrix<index_t> > actives_v(d);
        for(int i=0; i!=d;++i)
            system.mapColIndices(actives_v_Block(i,0), patchIndex,actives_v[i],i);
        system.mapColIndices(basisDataP.actives, patchIndex, basisDataP.actives,d);

        // Add contributions to the system matrix and right-hand side
        for(index_t i=0; i!=d;++i)
        {
            for(index_t j=0; j<d;++j)
                system.pushToMatrix(A_Block(i,j), actives_v[i],actives_v[j], eliminatedDofs[j], i, j);

            system.pushToMatrix(B_Block(0,i), basisDataP.actives, actives_v[i], eliminatedDofs[i], d, i);
            system.pushToMatrix(B_Block(0,i).transpose(),actives_v[i],basisDataP.actives, eliminatedDofs[d], i, d);

            system.pushToRhs( rhs_Block(i,0),actives_v[i],i);
        }
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



protected:
    gsAssemblerSetup m_info;

    // Velocity vector dimension
    index_t d;

    // viscosity
    T m_viscosity;

    // Right hand side
    const gsFunction<T> * rhs_ptr;

    //velocity vector basis
    memory::unique_ptr<gsVBasis<T> > vBasis;

protected:

    //Basis and function Data
    gsFuncData<T> basisDataV;
    gsFuncData<T> basisDataP;
    gsFuncData<T> funcData;

protected:
    // Local matrices
    gsMatrix<T> localMatA;
    gsMatrix<T> localMatB;
    gsMatrix<T> localRhs_u; //, localRhs_p;

    // Local values of the right hand side
    gsMatrix<T> rhsVals;
};


} // namespace gismo

