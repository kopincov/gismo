/** @file gsVisitorNonLinShell.h

    @brief Element visitor for nonlinear elasticity on thin shell.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, Y. Xia
*/

#pragma once

#include <gsAssembler/gsVisitorCDR.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

/** 
    @brief The visitor computes element non-linear load
    The Pde of the problem is : 
    \[ - \left[ {{k_1}\frac{{{\partial ^2}z}}{{\partial {x^2}}} 
    + 2{k_2}\frac{{{\partial ^2}z}}{{\partial x\partial y}} 
    + {k_3}\frac{{{\partial ^2}z}}{{\partial {y^2}}}} \right] 
    = {k_4}\sqrt {1 + {z_{,x}}^2 + {z_{,y}}^2} \]

    The Newton iteration formula of IGA solver is:
    \[\left( {{z^{\left( {n + 1} \right)}} 
    - {z^{\left( n \right)}}} \right)f'\left( {{z^{\left( n \right)}}} \right) 
    =  - f\left( {{z^{\left( n \right)}}} \right)\]

    in which
    \[f\left( z \right) = {S_{ij}} \cdot {z_j} - 
    {k_4}\int\limits_{{\Omega _e}} {{N_i} \cdot R{\rm{d}}\Omega } \]

    \[f'\left( z \right) = {S_{ij}} - {k_4}\int\limits_{{\Omega _e}} 
    {{R^{ - 1}} \cdot {N_i} \cdot \left( {\frac{{\partial z}}{{\partial x}}
    \cdot \frac{{\partial {N_j}}}{{\partial x}} + 
    \frac{{\partial z}}{{\partial y}} \cdot \frac{{\partial {N_j}}}
    {{\partial y}}} \right){\rm{d}}\Omega } \]
    
    and

    \[{S_{ij}} = {k_1}\int\limits_{{\Omega _e}} {\frac{{\partial {N_i}}}
    {{\partial x}}\frac{{\partial {N_j}}}{{\partial x}}{\rm{d}}\Omega }  
    + 2{k_2}\int\limits_{{\Omega _e}} {\frac{{\partial {N_i}}}{{\partial y}}\frac{{\partial {N_j}}}
    {{\partial x}}{\rm{d}}\Omega }  + {k_3}\int\limits_{{\Omega _e}} {\frac{{\partial {N_i}}}
    {{\partial y}}\frac{{\partial {N_j}}}{{\partial y}}{\rm{d}}\Omega } \]

    \[R = \sqrt {1 + \left[ {{{\left( {\frac{{\partial {N_j}}}
    {{\partial x}}{z_j}} \right)}^2} + {{\left( {\frac{{\partial {N_j}}}
    {{\partial y}}{z_j}} \right)}^2}} \right]} \]

    coef_b and coef_a should always set as zero.
    \tparam T coefficient type
    
    \ingroup SelfSuppSurf



*/
template <class T>
class gsVisitorNonLinLoad : public gsVisitorCDR<T> //to do: rename to VisitorNewton..
{
public:
    typedef gsVisitorCDR<T> Base;
public:

    /// Constructor with thickness, material parameters and the surface force as inputs.
    gsVisitorNonLinLoad(const gsMultiPatch<T> & curSolution, 
        const gsFunction<T> &coef_A = gsConstantFunction<T> (1, 0, 0, 1, 2),
        const gsFunction<T> &loa = gsConstantFunction<T>(1,2),
        const gsFunction<T> &coef_b = gsConstantFunction<T>(0,0,2),
        const gsFunction<T> &coef_c = gsConstantFunction<T>(0,2)) 
    : Base(loa, coef_A, coef_b, coef_c), sol_ptr(&curSolution)
    { 
    }
    
    /// Evaluate on element.
    inline void evaluate(const gsBasis<T>    & basis,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T>   & quNodes)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);
        numActive = actives.rows();
        
        // Evaluate basis functions on element
        basis.evalAllDers_into( md.points, 1, basisData);
        
        // Compute image of Gauss nodes under geometry mapping as well
        // as Jacobians
        geo.computeMap(md);

        // Evaluate current solution gradients
        sol_ptr->patch(geo.id()).jacobian_into(md.points, solGrads);

        // get the gradients to columns
        solGrads.transposeInPlace();
        solGrads.resize(md.points.rows(), md.points.cols() );

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive );
        localRhs.setZero(numActive, 1         );
    }

    /// Assembles the local matrix and right hand side
    inline void assemble(gsDomainIterator<T> & /*element*/,
                         const gsVector<T>   & quWeights)
    {
        const gsMatrix<T> & bVals  = basisData[0];
        const gsMatrix<T> & bGrads = basisData[1];

        const gsMatrix<T> & quCoords = md.values[0];
        
        gsMatrix<T> coeVal, rigVal;

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);
            
            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads, physBasisGrad);

            // Compute current solution physical gradients at k as a
            // Dim x NumActive matrix
            transformGradients(md, k, solGrads, physSolGrad);

            // Non-linear function 
            const T lv = - math::sqrt(1.0 + physSolGrad.squaredNorm());

            //coefficients on numerical integration points           
            coeff_A_ptr->eval_into(quCoords.col(k), coeVal);
            coeVal.resize(2, 2);

            rhs_ptr->eval_into(quCoords.col(k), rigVal);

            //gsDebug << "coef " << coeVal << coeVal.rows() << "\t" << coeVal.cols()  << "\n";
            //gsDebug << "right " << rigVal << "\t" << rigVal(0,0) << "\n";

            const gsMatrix<T> & rv =  physBasisGrad.transpose() * coeVal * physBasisGrad;

            localMat.noalias() += weight * (rv + rigVal(0,0) * ( 1.0 / lv ) 
                                 * bVals.col(k) * physSolGrad.transpose() * physBasisGrad );   

            localRhs.noalias() -= weight * (physBasisGrad.transpose() * coeVal * physSolGrad
                                  + rigVal(0,0) * lv * bVals.col(k));
        }
    }
    
    /// Local to global with homogeneous boundary
    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> > & /*eliminatedDofs*/,
                              gsSparseSystem<T>     & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives, 0, 0);
    }


protected:

    const gsMultiPatch<T> * sol_ptr;
    const gsFunction<T> *coeff_ptr;
    const gsFunction<T> *load_ptr;

protected:

    // Basis values
    using Base::basisData;
    using Base::actives;
    using Base::physBasisGrad;
    using Base::numActive;
    using Base::rhsVals;
    using Base::coeff_A_ptr;
    using Base::coeff_b_ptr;
    using Base::coeff_c_ptr;
    using Base::rhs_ptr;

    gsMatrix<T> physSolGrad, solGrads;

protected:
    // Local matrices
    using Base::localMat;
    using Base::localRhs;
    using Base::md;
};


template <class T>
class gsVisitorMasonryRhs : public gsVisitorNonLinLoad<T>
{
public:
    typedef gsVisitorNonLinLoad<T> Base;
public:

    /// Constructor with thickness, material parameters and the surface force as inputs.
    gsVisitorMasonryRhs(
        const gsMultiPatch<T> & curSolution, 
        const gsFunction<T>   & coef_A,
        const gsFunction<T>           & rhs) 
    : Base(curSolution,rhs,rhs)
    { 
        rhs_ptr = &rhs;
        coeff_A_ptr = & coef_A;
        sol_ptr = &curSolution;
    }
    
    /// Assembles the local matrix and right hand side
    inline void assemble(gsDomainIterator<T> & /*element*/,
                         const gsVector<T>   & quWeights)
    {
        const gsMatrix<T> & bVals  = basisData[0];
        const gsMatrix<T> & bGrads = basisData[1];
        const gsMatrix<T> & quCoords = md.values[0];
        
        gsMatrix<T> rigVal, coeVal;

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);
            
            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads, physBasisGrad);

            // Compute current solution physical gradients at k as a
            // Dim x NumActive matrix
            transformGradients(md, k, solGrads, physSolGrad);

            // Non-linear function 
            const T lv = math::sqrt(1.0 + physSolGrad.squaredNorm());

            //coefficients on numerical integration points           
            coeff_A_ptr->eval_into(quCoords.col(k), coeVal);
            coeVal.resize(2, 2);

            rhs_ptr->eval_into(quCoords.col(k), rigVal);

            localRhs.noalias() += ( weight * rigVal.value() * lv ) * bVals.col(k) ;
            localMat.noalias() += weight * (physBasisGrad.transpose() * coeVal * physBasisGrad);
        }

        //gsDebugVar(sum);
    }
     
    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        gsVisitorCDR<T>::localToGlobal(patchIndex, eliminatedDofs, system);
    }


protected:

     using Base::sol_ptr;

protected:

    // Basis values
    using Base::basisData;
    using Base::actives;
    using Base::physBasisGrad;
    using Base::numActive;
    using Base::rhsVals;
    using Base::coeff_A_ptr;
    using Base::coeff_b_ptr;
    using Base::coeff_c_ptr;
    using Base::rhs_ptr;

    using Base::solGrads;
    using Base::physSolGrad;

protected:
    // Local matrices
    using Base::localMat;
    using Base::localRhs;
    using Base::md;
};


} // namespace gismo

