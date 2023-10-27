/** @file gsVisitorPlateMixedPsiq.h

    @brief Boundary element visitor for the contribution of Psiq boundaries

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
*/


namespace gismo
{


template <class T>
class gsVisitorPlateMixedPsiq
{
public:

    gsVisitorPlateMixedPsiq(boxSide side)
        : side(side)
    { }


    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule,
                    unsigned         & evFlags )
    {
        // Setup Quadrature (harmless slicing occurs)
        rule = gsGaussRule<T>(basis, options.getReal("quA"), options.getInt("quB"), side.direction());

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_OUTER_NORMAL | NEED_GRAD_TRANSFORM;;
    }


    // Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const   & basis, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval,
                         // todo: add element here for efficiency
                         gsMatrix<T>            & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis[0].active_into(quNodes.col(0) , actives);
        numActive = actives.rows();

        // Evaluate basis functions on element
        basis[0].evalAllDers_into( quNodes, 1, basisData);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Initialize local matrix/rhs
        localMat.setZero(numActive, 2);
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];

        for (index_t k = 0; k < quWeights.rows(); ++k)
        {

            // Compute the outer normal vector on the side
            geoEval.outerNormal(k, side, unormal);

            // Multiply quadrature weight by the measure of normal
            const T weight = quWeights[k] * unormal.norm();

            unormal.normalize();
            gsVector<T,2> tangent;
            tangent(0) = -unormal(1);
            tangent(1) = unormal(0);

            // run over boundary counter-clockwise (as direction of tangent)
            //const T   sgn = sideOrientation(side) * geoEval.orientation();

            for (index_t i = 0; i < numActive; i++)
            {         
                localMat.row(i) += -weight * bVals(i,k)*unormal.transpose();
            }
        }
    }


    inline void localToGlobal(const index_t patchIndex,
                              gsSparseSystem<T> & system,
                              gsMatrix<T> & psiq)
    {
        size_t r = 0;
        system.mapColIndices(actives, patchIndex, actives_p, 0);

        const gsDofMapper & rowMap = system.rowMapper(r);

        for (index_t i = 0; i != numActive; ++i)
        {
            const int ii =  actives_p(i);
            if ( rowMap.is_free_index(actives_p(i)) )
            {
                psiq.row(ii) += localMat.row(i);
            }
        }

    }


    inline gsMatrix<T> getLocalMat()
    {
        return localMat;
    }

protected:


    // Neumann function
    const gsFunction<T> * pData_ptr;
    boxSide side;
    double penalty;

    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<index_t> actives;
    gsMatrix<index_t> actives_p;
    gsMatrix<T>		   physGrad;
    index_t            numActive;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> pVals;

    // Local matrix and rhs
    gsMatrix<T> localMat;


};


} // namespace gismo
