/** @file gsVisitorDg2.h

    @brief A DG interface visitor for the Poisson problem .

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer, A. Mantzflaris, S. Moore
*/

#pragma once

#include <gsPde/gsPoissonHeterogeneousPde.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{
/** @brief
    A new Implementation of a interface condition for the
    discontinuous Galerkin Assembler.

    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can also be imposed weakly (i.e Nitsche ).

    This visitor is also capable of handling inhomogeneous diffusion
    coefficients, however they have to be constant on each patch. It
    would need further development in order to make this right.
    (problem with assembling the penalty terms)

    Furthermore, one can choose between different dg version (symmetrized,
    non-symmetrized, anti-symmetric) by setting the symmetrizedDG value to
    1, 0, -1.

*/


//TODO: Clean up, check mappings
template <class T, int symmetrizedDG = 1>
class gsVisitorDg2
{
public:

    gsVisitorDg2(const gsPde<T> & pde, boxSide s1, boxSide s2)
        : m_domain(&pde.domain()),side1(s1), side2(s2), isHomogeneous(false)
    {
        const gsPoissonHeterogeneousPde<T>* ppde = static_cast<const gsPoissonHeterogeneousPde<T> *>(&pde);
        if(ppde == NULL)
            isHomogeneous = true;
        else
            m_alpha = ppde->getAlpha();

    }


    /** \brief Visitor for adding the interface conditions for the
 * interior penalty symmetrizeded of the poisson problem.
 *
 * This visitor adds the following term to the left-hand side (bilinear form).
 * \f[ - \{\nabla u\} \cdot \mathbf{n} [ v ]
 *     - \{\nabla v\} \cdot \mathbf{n} [ u ]
 *     + \alpha [ u ][  v ] \f].
 * Where \f[v\f] is the test function and \f[ u \f] is trial function.
 */
    //keep compatibility with other interface....
    gsVisitorDg2(T _penalty, boxSide s1, boxSide s2) :
        penalty(_penalty), side1(s1), side2(s2), d(0), isHomogeneous(true){ }


    gsVisitorDg2(T _penalty,const boundaryInterface & bI,const gsPiecewiseFunction<T>* alpha) :
        m_alpha(alpha), penalty(_penalty), side1(bI.first()), side2(bI.second()), d(0),isHomogeneous(false)
    {
        //gsInfo<< "Using dG with \n symmetrized "<<symmetrized<<"\n homogeneous diffusion coeff "<<isHomogeneous<<std::endl;
    }

    gsVisitorDg2(T _penalty,const boundaryInterface & bI) :
        penalty(_penalty), side1(bI.first()), side2(bI.second()), d(0), isHomogeneous(true)
    {
        //gsInfo<< "Using dG with \n symmetrized "<<symmetrized<<"\n homogeneous diffusion coeff "<<isHomogeneous<<std::endl;
    }

    void initialize(const gsBasis<T> & basis1,
                    const gsBasis<T> & ,
                    gsQuadRule<T> & rule,
                    unsigned & evFlags
                    )
    {

        //increase the penalty parameter for unsymmetrized version
        if(symmetrizedDG==0)
            penalty*=3;

        d = basis1.dim();
        const int dir = side1.direction();
        gsVector<index_t> numQuadNodes ( d );
        for (int i = 0; i < basis1.dim(); ++i)
            numQuadNodes[i] = basis1.degree(i) + 1;
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM;

        // Compute penalty parameter
        const int deg = basis1.maxDegree();
        penalty = (deg + basis1.dim()) * (deg + 1) * T(2.0);

    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T> const       & B1, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval1,
                         gsBasis<T> const       & B2, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval2,
                         gsMatrix<T>            & quNodes1,
                         gsMatrix<T>            & quNodes2)
    {
        // Compute the active basis functions
        B1.active_into(quNodes1.col(0), actives1);
        B2.active_into(quNodes2.col(0), actives2);
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        // Evaluate basis functions and their first derivatives
        B1.evalAllDers_into( quNodes1, 1, basisData1);
        B2.evalAllDers_into( quNodes2, 1, basisData2);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval1.evaluateAt(quNodes1);
        geoEval2.evaluateAt(quNodes2);

        if(!isHomogeneous)
        {
            m_alpha->piece(geoEval1.id()).eval_into(geoEval1.values(), alphaVals1);
            m_alpha->piece(geoEval2.id()).eval_into(geoEval2.values(), alphaVals2);
        }
        else
        {
            alphaVals1 = Eigen::Matrix<T, Dynamic,Dynamic>::Constant(1,quNodes1.cols(),1);
            alphaVals2 = Eigen::Matrix<T, Dynamic,Dynamic>::Constant(1,quNodes2.cols(),1);
        }

        // Initialize local matrices
        B11.setZero(numActive1, numActive1); B12.setZero(numActive1, numActive2);
        E11.setZero(numActive1, numActive1); E12.setZero(numActive1, numActive2);
        B22.setZero(numActive2, numActive2); B21.setZero(numActive2, numActive1);
        E22.setZero(numActive2, numActive2); E21.setZero(numActive2, numActive1);

        //TODO: This should move to initialize, but no information about the two patch indices is there.
        if(m_domain==NULL)
            m_patchDiameter1 = m_patchDiameter2 = 1;
        else
        {
             m_patchDiameter1 = (m_domain->patch(geoEval1.id()).coefAtCorner(boxCorner::getFirst(1)) - m_domain->patch(geoEval1.id()).coefAtCorner(boxCorner::getLast(1))).norm();
             m_patchDiameter2 = (m_domain->patch(geoEval2.id()).coefAtCorner(boxCorner::getFirst(1)) - m_domain->patch(geoEval2.id()).coefAtCorner(boxCorner::getLast(1))).norm();
        }
    }

    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element1,
                         gsDomainIterator<T>    & element2,
                         gsGeometryEvaluator<T> & geoEval1,
                         gsGeometryEvaluator<T> & geoEval2,
                         gsVector<T>            & quWeights)
    {
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector from patch1
            geoEval1.outerNormal(k, side1, unormal);

            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T weight = quWeights[k] * unormal.norm();

            // Compute the outer unit normal vector from patch1 in place
            unormal.normalize();

            // Take blocks of values and derivatives of basis functions
            //const typename gsMatrix<T>::Column val1 = basisData1.col(k);//bug
            const typename gsMatrix<T>::Block val1 = basisData1[0].block(0,k,numActive1,1);
            gsMatrix<T> & grads1 = basisData1[1];// all grads
            //const typename gsMatrix<T>::Column val2 = basisData2.col(k);//bug
            const typename gsMatrix<T>::Block val2 = basisData2[0].block(0,k,numActive2,1);
            gsMatrix<T> & grads2 = basisData2[1];// all grads

            // Transform the basis gradients
            geoEval1.transformGradients(k, grads1, phGrad1);
            geoEval2.transformGradients(k, grads2, phGrad2);

            // Compute element matrices
            const T c1     = weight * T(0.5);

            N1.noalias()   = alphaVals1(0,k) * unormal.transpose() * phGrad1;
            N2.noalias()   = alphaVals2(0,k) * unormal.transpose() * phGrad2;

            const T h1     = m_patchDiameter1*element1.getCellSize();
            const T h2     = m_patchDiameter2*element2.getCellSize();
            const T h12    = 2*(1./h1 + 1./h2);

            B11.noalias() += c1 * ( val1 * N1 );
            B12.noalias() += c1 * ( val1 * N2 );
            B22.noalias() -= c1 * ( val2 * N2 );
            B21.noalias() -= c1 * ( val2 * N1 );

            //diffusion coefficients for the penalty terms are included in the localToGlobal method.
            const T c2     = weight*penalty*(h12) ;

            E11.noalias() += c2 * ( val1 * val1.transpose() );
            E12.noalias() += c2 * ( val1 * val2.transpose() );
            E22.noalias() += c2 * ( val2 * val2.transpose() );
            E21.noalias() += c2 * ( val2 * val1.transpose() );

            /*
             //Right implemented, this could bring an improvement with gaps.
             // Remaining stuff from gaps, keeping for the while.
            if(true)
            {
                E21.noalias() += c2 * (-h1*h1*h1( val2 * N1 ) );
                E12.noalias() += c2 * (h1*h1*h1*( val1 * N2 ) );
            }
            */

        }
    }

    //map local to global
    void localToGlobal(const gsDofMapper  & mapper,
                       const gsMatrix<T>     & eliminatedDofs,
                       const index_t patch1,
                       const index_t patch2,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        /*
         *Remaining stuff from gaps, keeping for the while
        int sgn = 1;
        if(!meanValue && isGap)
            sgn = -1;
        */

        // Local DoFs to global DoFs
        mapper.localToGlobal(actives1, patch1, actives1);
        mapper.localToGlobal(actives2, patch2, actives2);
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        // TODO: This only works for constant alpha on each patch, one would need additional
        // Exx in order to perform a correct implementation. Here the average is used.
        //T alpha = (2*alphaVals1(0,0)*alphaVals2(0,0))/(alphaVals1(0,0)+alphaVals2(0,0));
        T alpha = (alphaVals1(0,0)+alphaVals2(0,0))/2;

        // Push element contributions 1-2 to the global matrix and load vector
        for (index_t j=0; j!=numActive1; ++j)
        {
            const index_t jj1 = actives1(j); // N1_j
            if ( mapper.is_free_index(jj1) )
            {
                for (index_t i=0; i!=numActive1; ++i)
                {
                    const index_t  ii1 = actives1(i); // N1_i

                    if ( mapper.is_free_index(ii1) )
                        //  if ( jj1 <= ii1 )
                        sysMatrix( ii1, jj1 ) -=  B11(i,j) + symmetrizedDG*B11(j,i) - alpha*E11(i,j);
                    else
                    {
                        rhsMatrix.row(jj1).noalias() += (B11(i,j) + symmetrizedDG*B11(j,i) - alpha*E11(i,j)) *
                                eliminatedDofs.row( mapper.global_to_bindex(ii1) );
                    }
                }
                for (index_t i=0; i!=numActive2; ++i)
                {
                    const index_t  ii2 = actives2(i); // N2_i

                    if ( mapper.is_free_index(ii2) )
                        //  if ( jj1 <= ii2 )
                        sysMatrix( ii2, jj1)  -=  B21(i,j) +symmetrizedDG*B12(j,i) + alpha*E21(i,j);
                    //sysMatrix( ii2, jj1)  -=  sgn*B21(i,j) +sgn*symmetrized*B12(j,i) + alpha*E21(i,j);

                    else
                    {
                        rhsMatrix.row(jj1).noalias() += (B21(i,j) +symmetrizedDG*B12(j,i) + alpha*E21(i,j)) *
                                eliminatedDofs.row( mapper.global_to_bindex(ii2) );
                    }
                }
            }
        }
        // Push element contributions 2-1 to the global matrix and load vector
        for (index_t j=0; j!=numActive2; ++j)
        {
            const index_t jj2 = actives2(j); // N2_j
            if ( mapper.is_free_index(jj2) )
            {
                for (index_t i=0; i!=numActive2; ++i)
                {
                    const index_t  ii2 = actives2(i); // N2_i

                    if ( mapper.is_free_index(ii2) )
                        //  if ( jj2 <= ii2 )
                        sysMatrix( ii2, jj2 ) -=  B22(i,j) + symmetrizedDG*B22(j,i) - alpha*E22(i,j);
                    //sysMatrix( ii2, jj2 ) +=  sgn*B22(i,j) + sgn*symmetrized*B22(j,i) - alpha*E22(i,j);

                    else
                    {
                        rhsMatrix.row(jj2).noalias() += (B22(i,j) + symmetrizedDG*B22(j,i) - alpha*E22(i,j)) *
                                eliminatedDofs.row( mapper.global_to_bindex(ii2) );
                        //rhsMatrix.row(jj2).noalias() += (sgn*B22(i,j) + sgn*symmetrized*B22(j,i) - alpha*E22(i,j)) *
                        //       eliminatedDofs.row( mapper.global_to_bindex(ii2) );
                    }
                }

                for (index_t i=0; i!=numActive1; ++i)
                {
                    const index_t  ii1 = actives1(i); // N1_i

                    if ( mapper.is_free_index(ii1) )
                        // if ( jj2 <= ii1 )
                        sysMatrix( ii1, jj2)  -=  B12(i,j) +symmetrizedDG*B21(j,i) + alpha*E12(i,j);
                    //sysMatrix( ii1, jj2)  -=  sgn*B12(i,j) +sgn*symmetrized*B21(j,i) + alpha*E12(i,j);
                    else
                    {
                        //rhsMatrix.row(jj2).noalias() += (sgn*B12(i,j) +sgn*symmetrized*B21(j,i) + alpha*E12(i,j)) *
                        //        eliminatedDofs.row( mapper.global_to_bindex(ii1) );
                        rhsMatrix.row(jj2).noalias() += (B12(i,j) +symmetrizedDG*B21(j,i) + alpha*E12(i,j)) *
                                eliminatedDofs.row( mapper.global_to_bindex(ii1) );
                    }
                }
            }
        }
    }

    inline void localToGlobal(const index_t patchIndex1,
                              const index_t patchIndex2,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives1, patchIndex1, actives1);
        system.mapColIndices(actives2, patchIndex2, actives2);

        // Add contributions to the system matrix and right-hand side

        // TODO: This only works for constant alpha on each patch, one would need additional
        // Exx in order to perform a correct implementation. Here the average is used.
        //T alpha = (2*alphaVals1(0,0)*alphaVals2(0,0))/(alphaVals1(0,0)+alphaVals2(0,0));
        T alpha = (alphaVals1(0,0)+alphaVals2(0,0))/2;

        gsMatrix<T> localRhs1,localRhs2;
        localRhs1.setZero(actives1.rows(),1);
        localRhs2.setZero(actives2.rows(),1);

        system.push(-B11 - symmetrizedDG*B11.transpose() + alpha*E11, localRhs1,actives1,actives1,eliminatedDofs.front(),0,0);
        system.push(-B21 - symmetrizedDG*B12.transpose() - alpha*E21, localRhs2,actives2,actives1,eliminatedDofs.front(),0,0);
        system.push(-B12 - symmetrizedDG*B21.transpose() - alpha*E12, localRhs1,actives1,actives2,eliminatedDofs.front(),0,0);
        system.push(-B22 - symmetrizedDG*B22.transpose() + alpha*E22, localRhs2,actives2,actives2,eliminatedDofs.front(),0,0);


    }

    //map local to global for the IETIdG algorithm
    void localToGlobalIETI(const gsDofMapper  & mapper,
                           const gsMatrix<T>     & eliminatedDofs,
                           const gsMatrix<index_t>& activeExtra,
                           gsSparseMatrix<T>     & sysMatrix,
                           gsMatrix<T>           & rhsMatrix)
    {
        //Always use 1, they actives1 is swaped with actives2 for the second call.
        mapper.localToGlobal(actives1, 0, actives1);

        const index_t numActive = B11.rows(); //here we cannot use the rows of active1, because it is modified outside.

        index_t numActiveExtra = activeExtra.rows();

        // TODO: This only works for constant alpha on each patch (anyway, IETI can only work with this in a nice way),
        // one would need additional Exx in order to perform a correct implementation.
        T alpha = (alphaVals1(0,0)/2);

        // Push element contributions 1-2 to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const index_t jj1 = actives1(j); // N1_j
            if ( mapper.is_free_index(jj1) )
            {
                for (index_t i=0; i!=numActive; ++i)
                {
                    const index_t  ii1 = actives1(i); // N1_i
                    if ( mapper.is_free_index(ii1) )
                        //if ( jj1 <= ii1 )
                        sysMatrix( ii1, jj1 ) -=  B11(i,j) + B11(j,i) - alpha*E11(i,j);
                    else
                    {
                        rhsMatrix.row(jj1).noalias() += (B11(i,j) + B11(j,i) - alpha*E11(i,j)) *
                                eliminatedDofs.row( mapper.global_to_bindex(ii1) );
                    }
                }

                for (index_t i=0; i!=numActiveExtra; ++i)
                {
                    const index_t  ii2 = actives1(numActive+activeExtra(i)); // N2_i
                    if ( mapper.is_free_index(ii2) )
                        //if ( jj1 <= ii2 )
                        sysMatrix( ii2, jj1)  -=  B21(activeExtra(i),j) + alpha*E21(activeExtra(i),j);
                    else
                    {
                        rhsMatrix.row(jj1).noalias() += (B21(activeExtra(i),j) + alpha*E21(activeExtra(i),j)) *
                                eliminatedDofs.row( mapper.global_to_bindex(ii2));
                    }
                }
            }
        }
        for (index_t j=0; j!=numActiveExtra; ++j)
        {
            const index_t  jj1 = actives1(numActive+activeExtra(j)); // N2_i

            if ( mapper.is_free_index(jj1) )
            {
                for (index_t i=0; i!=numActive; ++i)
                {
                    const index_t  ii1 = actives1(i); // N1_i
                    if ( mapper.is_free_index(ii1) )
                        //if ( jj1 <= ii1 )
                        sysMatrix( ii1, jj1 ) -=  B21(activeExtra(j),i) + alpha*E21(activeExtra(j),i);
                    else
                    {
                        rhsMatrix.row(jj1).noalias() += (B21(activeExtra(j),i) + alpha*E21(activeExtra(j),i)) *
                                eliminatedDofs.row( mapper.global_to_bindex(ii1) );
                    }
                }

                for (index_t i=0; i!=numActiveExtra; ++i)
                {
                    const index_t  ii2 = actives1(numActive+activeExtra(i)); // N2_i
                    if ( mapper.is_free_index(ii2) )
                        //if ( jj1 <= ii2 )
                        sysMatrix( ii2, jj1)  -= - alpha*E22(activeExtra(i),activeExtra(j));
                    else
                    {
                        rhsMatrix.row(jj1).noalias() += (-alpha*E22(activeExtra(i),activeExtra(j))) *
                                eliminatedDofs.row( mapper.global_to_bindex(ii2));
                    }
                }
            }
        }
    }

    // assemble on element
    inline void assembleSingleSide(gsDomainIterator<T>    & element1,
                         gsDomainIterator<T>    & element2,
                         gsGeometryEvaluator<T> & geoEval1,
                         gsGeometryEvaluator<T> & geoEval2,
                         gsVector<T>            & quWeights)
    {
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector from patch1
            geoEval1.outerNormal(k, side1, unormal);

            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T weight = quWeights[k] * unormal.norm();

            // Compute the outer unit normal vector from patch1 in place
            unormal.normalize();

            // Take blocks of values and derivatives of basis functions
            //const typename gsMatrix<T>::Column val1 = basisData1.col(k);//bug
            const typename gsMatrix<T>::Block val1 = basisData1[0].block(0,k,numActive1,1);
            gsMatrix<T> & grads1 = basisData1[1];// all grads
            //const typename gsMatrix<T>::Column val2 = basisData2.col(k);//bug
            const typename gsMatrix<T>::Block val2 = basisData2[0].block(0,k,numActive2,1);

            // Transform the basis gradients
            geoEval1.transformGradients(k, grads1, phGrad1);

            // Compute element matrices
            const T c1     = weight * T(0.5);

            N1.noalias()   = alphaVals1(0,k) * unormal.transpose() * phGrad1;

            const T h1     = m_patchDiameter1*element1.getCellSize();
            const T h2     = m_patchDiameter2*element2.getCellSize();
            const T h12    = 2*(1./h1 + 1./h2);

            B11.noalias() += c1 * ( val1 * N1 );
            B21.noalias() -= c1 * ( val2 * N1 );

            //diffusion coefficients for the penalty terms are included in the localToGlobal method.
            const T c2     = weight*penalty*(h12) ;

            E11.noalias() += c2 * ( val1 * val1.transpose() );
            E22.noalias() += c2 * ( val2 * val2.transpose() );
            E21.noalias() += c2 * ( val2 * val1.transpose() );

            /*
             //Right implemented, this could bring an improvement with gaps.
             // Remaining stuff from gaps, keeping for the while.
            if(true)
            {
                E21.noalias() += c2 * (-h1*h1*h1( val2 * N1 ) );
                E12.noalias() += c2 * (h1*h1*h1*( val1 * N2 ) );
            }
            */

        }
    }


    /**
     * @brief getActives return a reference to the pointer to the set of active basis functions.
     *     (the same effect, as they would be public). Used in the IETIdG Algorithm.
     * @param activesParam1
     * @param activesParam2
     */
    void getActives(gsMatrix<index_t>*& activesParam1, gsMatrix<index_t>*& activesParam2)
    {
        activesParam1 = &(this->actives1);
        activesParam2 = &(this->actives2);
    }

    /**
     * @brief revert swaps the content of B11 with B22, ... etc., used in the IETIdg Algorithm
     * in order to use one localToGlobal function, which is then called twice.
     */
    void revert(bool extra =false)
    {
        B11.swap(B22);
        B21.swap(B12);
        E11.swap(E22);
        E21.swap(E12);
        actives1.swap(actives2);
        alphaVals1.swap(alphaVals2);
        if(extra)
        {
            std::swap(side1,side2);
            basisData1.swap(basisData2);
        }
    }

private:
    //Diffusion coefficient
    const gsPiecewiseFunction<T>* m_alpha;

    const gsMultiPatch<T>* m_domain;

    // Penalty constant
    T penalty;

    // Side
    boxSide side1, side2;

    // dimension of the problem
    unsigned d;

    // versions of dG: 1 .. symmetrizied, 0 .. unsymmetrized, -1 .. antisymmetrized(?)

    // if we use homogeneous diffusion coefficients (==1).
    bool isHomogeneous;
private:

    T m_patchDiameter1,m_patchDiameter2;

    // Basis values etc
    gsMatrix<T>        alphaVals1, alphaVals2;
    std::vector<gsMatrix<T> > basisData1, basisData2;
    gsMatrix<T>        phGrad1   , phGrad2;
    gsMatrix<index_t> actives1  , actives2;

    // Outer normal
    gsVector<T> unormal;

    // Auxiliary element matrices
    gsMatrix<T> B11, B12, E11, E12, N1,
    B22, B21, E22, E21, N2;
};


} // namespace gismo
