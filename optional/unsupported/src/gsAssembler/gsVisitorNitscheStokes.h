

#include <gsPde/gsStokesPde.h>
#include <gsCore/gsVectorValuedFunctionSet.h>
#include <gsCore/gsTransformPure.h>

#pragma once


namespace gismo
{


/** \brief Visitor for the weak imposition of the dirichlet boundary condition.
  This Visitor also covers the use DivConforming Transformation.
  Note: it is required, that each

  boundary condition stores only a scalar function. The algorithm
  will crash if not so.

  However, you are allowed to loop only once around the boundary.
  The visitor will look by itself for the remaining components.
 */
template <class T>
class gsVisitorNitscheStokes
{
public:

    /**
     * @brief gsVisitorNitscheStokes Constructor for the visitor
     * @param pde the Stokes PDE
     * @param s the considered BC
     */
    gsVisitorNitscheStokes(const gsPde<T> & pde, const boundary_condition<T> & s)
        : side(s.side())
    {
        const gsStokesPde<T>* spde = static_cast<const gsStokesPde<T>*>(&pde);

        m_viscosity =spde->viscocity();

        //Here we assume that all components of v on this side are Dirichlet dofs
        dirdata_ptr.resize(pde.unknownDim()(0));

        //We look for the remaining components of the Neumann BC for this unknown
        //and store them
        for ( typename gsBoundaryConditions<T>::const_iterator
              it = pde.bc().dirichletBegin();
              it != pde.bc().dirichletEnd(); ++it )
        {
            if((*it).ps == s.ps)
                dirdata_ptr[(*it).unknown()]= (*it).function().get();
        }
        for(size_t i=0; i<dirdata_ptr.size();++i )
            GISMO_ASSERT(dirdata_ptr[i]->targetDim()==1, "The Algorithm cannot handle vector valued BC, each component should be and individual BC ");
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

        // Compute penalty parameter
        const int deg = basis[m_opt.info.getBasis(0,1)].degree(0);
        penalty = T(5.0)*(deg + 1); //(deg + basis.dim()) * (deg + 1) * T(2.0);

        // Set Geometry evaluation flags
        evFlagsMap |= NEED_VALUE | NEED_MEASURE |NEED_GRAD_TRANSFORM|NEED_OUTER_NORMAL;

        // Set Basis evaluation flag
        unsigned evFlagsFunc = SAME_ELEMENT| NEED_ACTIVE | NEED_VALUE;

        basisDataV.addFlags(evFlagsFunc|NEED_DERIV);
        basisDataP.addFlags(evFlagsFunc);

         //Choose if V is divconforming or not and add corresponding flags
        if(m_opt.info.getTransform(0)==transform::Hdiv)
            evFlagsMap |= gsTransformPureDivConforming<T>::getGeoFlags(basisDataV.flags);
        else
            evFlagsMap |= gsTransformPureGradConforming<T>::getGeoFlags(basisDataV.flags);

        evFlagsMap |= gsTransformPureGradConforming<T>::getGeoFlags(basisDataP.flags);

        dirData.addFlags(SAME_ELEMENT|NEED_VALUE);

        //We build our vector basis
        std::vector<gsBasis<T>* > v_bases(d);
        //The const cast is required in order to extract certein components of basisRefs and
        //store them in a std::vector to call the gsBasis. It would be nice to have a direct
        //pipeline from gsBasisRefs to gsVBasis, where the used components are specified.
        for(int comp_v = 0; comp_v < d; ++comp_v) // all components for unknown 0
            v_bases[comp_v]=const_cast<gsBasis<T>*>(&basis[m_opt.info.getBasis(0,comp_v)]);

        //Set shifts to zero, there are stored in gsSparseSystem
        std::vector<index_t> shifts(d,0);
        vBasis = memory::make_unique(new gsVBasis<T>(v_bases, shifts));
    }

    /// Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
                         const gsMapData<T> & mapData)
    {
        // set the preasure basis
        const gsBasis<T>& pBasis = basisRefs[m_opt.info.getBasis(1,0)]; //unknown 1, its zero component (the only one)

        //dependent on the transformation, eval the basis and transform
        if(m_opt.info.getTransform(0)==transform::Hdiv)
            gsTransformPureDivConforming<T>::transform(mapData,*vBasis,basisDataV);
        else
            gsTransformPureGradConforming<T>::transform(mapData,*vBasis,basisDataV);

        //The GradcConforming for preasure
        gsTransformPureGradConforming<T>::transform(mapData,pBasis,basisDataP);

        unsigned numActV = basisDataV.actives.rows();
        unsigned numActP = basisDataP.actives.rows();

        //The vector valued Neumann values
        m_dirVals.resize(d,mapData.points.cols());

        // Evaluate the Dirichlet data
        for(index_t i=0; i<d;++i)
        {
            dirdata_ptr[i]->compute(mapData.values[0], dirData);
            m_dirVals.row(i)=dirData.values[0];
        }

        // Initialize local matrix/rhs
        localMatA .setZero(numActV, numActV);
        localRhs_u.setZero(numActV, 1);
        localRhs_p.setZero(numActP, 1);
        localMatB.setZero(numActP, numActV);


        unormalLarge.setZero(d*d,d);
    }

    /// assemble the local contributions
    inline void assemble(gsDomainIterator<T>    & element,
                         const gsMapData<T>     & mapData,
                         const gsVector<T>      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const gsMatrix<T> & vGrads= basisDataV.deriv(k); //dim*dim x n_actV
            const gsMatrix<T> & vVals = basisDataV.eval(k); // 1 x n_actV
            const gsMatrix<T> & pVals = basisDataP.eval(k); // 1 x n_actP
            const gsMatrix<T> & dirVals = m_dirVals.col(k); //dim x 1.

            // Compute the outer normal vector on the side
            unormal = mapData.outNormal(k); //dim x 1
            const T weight =  quWeights[k] * unormal.norm();

            //We store the normal, because the result would be const and cannot be normalized!
            unormal.normalize();
            //Make it larger, be able to be multiplied with vGrads
            if(d==2)
                unormalLarge.block(0,0,d,1) = unormalLarge.block(d,1,d,1)= unormal;
            else if(d==3)
                unormalLarge.block(0,0,d,1) = unormalLarge.block(d,1,d,1)= unormalLarge.block(2*d,2,d,1)= unormal;

            // Get penalty parameter
            const T mu = penalty / element.getPerpendicularCellSize();

            // Sum up quadrature point evaluations
            localRhs_u.noalias() -= weight * m_viscosity* ((vGrads.transpose() * unormalLarge - mu * vVals.transpose() )
                                                           * dirVals );

            localRhs_p.noalias() -=weight  * pVals.transpose()* unormal.transpose() * dirVals;

            localMatA.noalias() -= weight * m_viscosity* ( vVals.transpose() * (unormalLarge.transpose() * vGrads)
                                                           +  (vVals.transpose()* (unormalLarge.transpose() * vGrads)).transpose()
                                                           -  mu * vVals.transpose() * vVals );

            localMatB.noalias() -= weight* pVals.transpose()*( unormal.transpose() * vVals );
        }
    }

    /// push the local contributions
    inline void localToGlobal(const index_t                     patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>               & system)
    {
        //Make vectors for block structure
        gsVector<index_t> dBlocks(d),oneCol(1);
        for(int i =0; i<d;++i)
            dBlocks(i) = getActiveBlockLength(i);
        oneCol<<1;

        //deduce block structure
        typename gsMatrix<index_t>::BlockView actives_v_Block = basisDataV.actives.blockView(dBlocks, oneCol);
        typename gsMatrix<T>::BlockView A_Block = localMatA.blockView(dBlocks,dBlocks);
        typename gsMatrix<T>::BlockView rhs_Block = localRhs_u.blockView(dBlocks,oneCol);


        // Map patch-local DoFs to global DoFs, a blockview cannot be overridden (idk why)
        std::vector<gsMatrix<index_t> > actives_v(d);
        for(int i=0; i!=d;++i)
            system.mapColIndices(actives_v_Block(i,0), patchIndex,actives_v[i],i);

        // Add contributions to the system matrix and right-hand side
        for(index_t i=0; i!=d;++i)
        {
            system.pushToRhs( rhs_Block(i,0),actives_v[i],i);
            for(index_t j=0; j<d;++j)
                system.pushToMatrix(A_Block(i,j), actives_v[i],actives_v[j], eliminatedDofs[j], i, j);
        }



        //If you use nitsche for the incorporating also the normal, do additionally the following
        //or you use just nitsche for the classical stokes
        if (!(m_opt.dirStrategy == dirichlet::eliminatNormal))
        {
            gsMatrix<index_t>& actives_p = basisDataP.actives;

            //auxiliary vector for blockstructures
            gsVector<index_t> oneBlock(1);
            oneBlock<<actives_p.rows();
            typename gsMatrix<T>::BlockView B_Block = localMatB.blockView(oneBlock,dBlocks);

            //map the preasure actives
            system.mapColIndices(actives_p, patchIndex, actives_p,d);

            //Push to matrix and RHS
            for(index_t i=0; i!=d;++i)
            {
                system.pushToMatrix(B_Block(0,i), actives_p, actives_v[i], eliminatedDofs[i], d, i);
                system.pushToMatrix(B_Block(0,i).transpose(),actives_v[i],actives_p, eliminatedDofs[d], i, d);
            }
            system.pushToRhs( localRhs_p,actives_p,d);
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
private:

    gsAssemblerOptions m_opt;

    // Velocity vector dimension
    index_t d;

    // Side
    boxSide side;

    // Dirichlet function
    std::vector<gsFunction<T> *> dirdata_ptr;

    //velocity vector basis
    memory::unique_ptr<gsVBasis<T> > vBasis;

    // Penalty constant
    T penalty;

    //viscosity
    T m_viscosity;
private:

    //Basis and function Data
    gsFuncData<T> basisDataV;
    gsFuncData<T> basisDataP;

    gsFuncData<T> dirData;
    gsMatrix<T> m_dirVals;

    // Normal
    gsVector<T> unormal;

    // ( n , 0 )
    // ( 0,  n )
    gsMatrix<T> unormalLarge;


protected:
    // Local matrices
    gsMatrix<T> localMatA;
    gsMatrix<T> localMatB;
    gsMatrix<T> localRhs_u, localRhs_p;


};

} // namespace gismo
