/** @file gsMasoStress.hpp

    @brief Provides implementation for the Masonry problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):  Y. Xia, A. Mantzaflaris
*/
// #include <gsSelfSuppSurf/gsVisitorNonLinLoad.h>
// #include <vector>
// using std::vector;

namespace gismo
{
/*
template<class T> 
void gsMasoStress::CalStress(const   gsMultiPatch<T> & m_curSolution, 
                             vector<gsMatrix<T> > & stressPatches)
{
    
         
    gsVisitorCDR2<T> visitCDR( *m_rhsFun, *m_coeff_b, *m_coeff_c, m_flagStabilization);
    for (index_t patchIndex = 0; patchIndex< m_curSolution->nPatches(); ++patchIndex)
    {  
        const gsBasisRefs<T> bases(m_bases, patchIndex);
        const gsDofMappers mappers(m_dofMappers);
        gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
        gsVector<T> quWeights; // Temp variable for mapped weights
        gsQuadRule<T> QuRule; // Reference Quadrature rule
        unsigned evFlags(0);

        // Initialize 
        visitCDR.initialize(bases, QuRule, evFlags);
       
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_patches[patchIndex]));

        // Initialize domain element iterator 
        typename gsBasis<T>::domainIter domIt = bases[0].makeDomainIterator(side);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );
            // Perform required evaluations on the quadrature nodes
            visitCDR.evalResultHessian(bases, *geoEval, quNodes);


        }
        


    }
   
}
 */
// template<class T> 
// void gsMasoStress<T>::CalGeoHessian()
// {
// 
// 
// 
// }
/*
template<class T> 
void gsMasoStress<T>::applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<index_t> acts;

    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
        //m_patches
        if ( m_pLoads[i].parametric )
        {
            m_bases.front().basis(m_pLoads[i].patch).active_into( m_pLoads[i].point, acts );
            m_bases.front().basis(m_pLoads[i].patch).eval_into  ( m_pLoads[i].point, bVals);
        }
        else
        {
            gsWarn<< "Point loads parametric for now.\n";
        }

        // translate patch-local indices to global dof indices
        m_dofMappers.front().localToGlobal(acts, m_pLoads[i].patch, acts);

        for (index_t k=0; k < acts.rows(); ++k)
        {
            m_rhs(acts(k,0), 0) += bVals(k,0) * m_pLoads[i].value[2];
        }
    }

}*/

}//namespace gismo