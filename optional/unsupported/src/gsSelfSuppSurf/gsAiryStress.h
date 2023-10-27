/** @file gsMasoStress.h

    @brief Provides assembler for calculation of the stress in Masonry problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):  Y. Xia, A. Mantzaflaris.
*/
#pragma once
#include <gsSelfSuppSurf/gsMasonry.h>
#include <gsAssembler/gsVisitorCDR.h>
#include <gsAssembler/gsCDRAssembler.h>
#include <gsSelfSuppSurf/gsVisitorNonLinLoad.h>
#include <gsPde/gsPointLoads.h>
#include <vector>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

/** @brief
    Simple class to evaluate the diffusion coefficient coming from an Airy stress functional.
    
    \ingroup function
*/
template <typename T>
class gsAiryStressCoeff : public gsFunction<T>
{
public:
    gsAiryStressCoeff() : m_phi(NULL) { }

    /// Construct by a stress functional \a phi
    explicit gsAiryStressCoeff(const gsFunction<T> & phi)
    : m_phi(&phi)
    { }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT( m_phi!=NULL, "No function given.");
        result.resize(4,u.cols() );
        gsMatrix<T> tmp;
        
        for (index_t i = 0; i!=u.cols(); ++i )
        {
            m_phi->deriv2_into(u.col(i), tmp);

            result(0,i) =  tmp(1,0);
            result(1,i) = 
            result(2,i) = -tmp(2,0);
            result(3,i) =  tmp(0,0);
        }
    }
    
    void setFunction(const gsFunction<T> & phi)
    {
        m_phi = &phi;
    }

    virtual short_t domainDim () const {return 2;}

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "Airy stress represented by: "<< *m_phi <<"\n";
        return os; 
    }


    private:
        const gsFunction<T> * m_phi;

};


    /** @brief
    Implementation of an (multiple right-hand side) masonry assembler.    
    \ingroup Assembler
*/
    
template<class T>
class gsVisitorAiryStress : public gsVisitorCDR<T>
{
    typedef gsVisitorCDR<T> base;

public:
    gsVisitorAiryStress(const gsFunction<T> & rhs,
        const gsFunction<T> & coeff_a,
        const gsFunction<T> & coeff_b,
        const gsFunction<T> & coeff_c,
        const gsMatrix<T>   & coor_z,
        stabilizerCDR::method flagStabilization = stabilizerCDR::none):
    base(rhs, coeff_a, coeff_b, coeff_c, flagStabilization )
    {
        m_zCoors = coor_z;
    }


    void initialize(const gsBasis<T>    & basis,
                    const index_t         /*patchIndex*/,
                    const gsOptionList  & options,
                          gsQuadRule<T> & rule)
    {
        // Setup Quadrature
        rule = gsGaussRule<T>(basis, options);// harmless slicing occurs here
        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER | NEED_JACOBIAN ;
    }

    inline void evaluate(const gsBasis<T>    & basis,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T>   & quNodes)
    {
        md.points = quNodes;

        //call evaluate function of gsVisitorCDR
        base::evaluate(basis, geo, md.points);
        //Evaluate the first derivatives of the bases
        basis.deriv_into(md.points, m_BasisGrad);
        //Evaluate the second derivatives of the bases
        basis.deriv2_into(md.points, m_basisHessian);
        gsDebug << "md.points" << "\t" << md.points << "\t";
        gsDebug << "gard" << "\n";
        gsDebug << m_BasisGrad.rows() << m_BasisGrad.cols() << "\n";
        gsDebug << "hessian" << "\n";
        gsDebug << m_basisHessian.rows() << m_basisHessian.cols()  << "\n";
        gsDebug << m_basisHessian << "\n";
        //calculate the physical Hessian matrix
        unsigned geoDim = 2;//there should be a better way
        unsigned nrOfActives = actives.rows();
        m_phyHessian.setZero((geoDim*(geoDim+1)/2), nrOfActives * md.points.cols());
        gsMatrix<T> pointPhyHessian;
        gsMatrix<T> pointPhyGrad;
        m_phyGrad.setZero(geoDim, nrOfActives * md.points.cols());
        for (index_t i=0; i<md.points.cols(); ++i)
        { 
            transformGradients(md, i, m_BasisGrad, pointPhyGrad);
            transformDeriv2Hgrad(md, i,m_BasisGrad, m_basisHessian, pointPhyHessian);
            gsDebugVar(pointPhyHessian);
            pointPhyHessian.transposeInPlace();
            m_phyHessian.block(0,i*nrOfActives,(geoDim*(geoDim+1)/2),nrOfActives) = pointPhyHessian;
            m_phyGrad.block(0,i*nrOfActives,geoDim,nrOfActives) = pointPhyGrad;
        }
        gsDebug << m_phyHessian.rows() << m_phyHessian.cols() << "\n";
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T>   & quWeights)
    {
        gsMatrix<T> & basisVals   = basisData[0];
        gsMatrix<T> &  basisGrads = basisData[1];

        const gsMatrix<index_t> & actDofEle = actives; // are you sure?
        gsDebug << "actDofEle " << actDofEle;
        const unsigned geoD = element.dim();

        m_zCoorEle.setZero(actDofEle.rows(),1);    
        index_t nrOfActives = actDofEle.size();
        for (index_t i=0; i<nrOfActives; ++i)
        {
            m_zCoorEle(i) = m_zCoors(actDofEle.at(i));
        }
        gsDebugVar(m_zCoorEle);
        //Get the z coordinates of the control points
        
        coeff_A_vals.setZero(4, quWeights.rows());
        for (index_t i=0; i<quWeights.rows(); ++i)
        {
            coeff_A_vals.block(0,i,1,1) = m_phyHessian.block(1, i * nrOfActives, 1, nrOfActives) * m_zCoorEle ;
            gsDebugVar(m_phyHessian.block(1, i * nrOfActives, 1, nrOfActives) );
            coeff_A_vals.block(1,i,1,1) = -1 * m_phyHessian.block(2, i * nrOfActives, 1, nrOfActives) * m_zCoorEle ;
            coeff_A_vals.block(2,i,1,1) = coeff_A_vals.block(1, i, 1, 1);
            coeff_A_vals.block(3,i,1,1) = m_phyHessian.block(0, i * nrOfActives, 1, nrOfActives) * m_zCoorEle ;
        }
        //Calculate the physical gradients
        gsDebugVar(coeff_A_vals);
        m_Z_phyGrad.setZero(2, quWeights.rows());
        for (index_t i=0; i<quWeights.rows(); ++i)
        {   
            m_Z_phyGrad.col(i) = m_phyGrad.block(0,i*nrOfActives, 2, nrOfActives)  * m_zCoorEle;
        }
        gsDebugVar(m_Z_phyGrad);
        //modify the values of loads
        gsDebug << "right hand side values" << rhsVals << "\n";
        for (index_t i=0; i<nrOfActives; ++i)
        {
            //uncomment this two lines if not in the linear test
            rhsVals(0, i) *= math::sqrt(1+ m_Z_phyGrad(0, i) * m_Z_phyGrad(0, i) 
                 + m_Z_phyGrad(1, i) * m_Z_phyGrad(1, i));
        }
        gsDebug << "right hand side values" << rhsVals << "\n";
        gsDebugVar(coeff_A_vals);        

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);

            // Compute physical gradients at k as a Dim x numActive matrix
            transformGradients(md, k, basisGrads, physBasisGrad);
            gsMatrix<T> tmp_A = coeff_A_vals.col(k);
            tmp_A.resize(geoD,geoD);            
            gsMatrix<T> b_basisGrads = coeff_b_vals.col(k).transpose() * physBasisGrad;
            
            localRhs.noalias() += weight * ( basisVals.col(k) * rhsVals.col(k).transpose() ) ;

            // ( N x d ) * ( d x d ) * ( d x N ) = N x N
            localMat.noalias() += weight * (physBasisGrad.transpose() * ( tmp_A * physBasisGrad) );
            // ( N x 1 ) * ( 1 x N) = N x N
            localMat.noalias() += weight * (basisVals.col(k) * b_basisGrads);
            // ( scalar ) * ( N x 1 ) * ( 1 x N ) = N x N
            localMat.noalias() += weight * coeff_c_vals(0,k) * (basisVals.col(k) * basisVals.col(k).transpose());                    
        }
    }

/*
    inline void evalResultHessian(gsBasis<T> const       & basis, 
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        //get solution gradients
        base::evaluate(basis, geoEval, quNodes);

        for (index_t k = 0; k < quNodes.rows(); ++k) // loop over quadrature nodes
        {           
            // Compute current solution physical Hessian at k as a
            geoEval.transformGradients(k, m_SolGrads, m_PhysSolGrad);        ///??? need modification   
        }
    }
  

    void calPhysicalHessian(gsBasis<T> const       & basis, 
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        //call evaluate function of gsVisitorCDR
        base::evaluate(basis, geoEval, quNodes);

        //Evaluate the first derivatives of the bases
        //geoEval.deriv_into(quNodes, m_BasisGrad);
        basis.deriv_into(quNodes, m_BasisGrad);
        //Evaluate the second derivatives of the bases
        basis.deriv2_into(quNodes, m_basisHessian);

        //calculate the physical Hessian matrix
        unsigned geoDim = 3;
        //unsigned nrOfActives = basis.numActive(); actives.rows();
        unsigned nrOfActives = actives.rows();

        m_phyHessian.setZero((geoDim*(geoDim+1)/2), nrOfActives * quNodes.cols());
        gsMatrix<T> pointPhyHessian;
        gsMatrix<T> pointPhyGrad;
        m_phyGrad.setZero(geoDim, nrOfActives * quNodes.cols());
        for (index_t i=0; i<quNodes.cols(); ++i)
        { 
            geoEval.transformGradients(i, m_BasisGrad, pointPhyGrad);
            geoEval.transformDeriv2Hgrad(i,m_BasisGrad, m_basisHessian, pointPhyHessian);
            gsDebug << (pointPhyHessian) << "\t pointPhyHessian " << "\n";
            pointPhyHessian.transposeInPlace();
            m_phyHessian.block(0,i*nrOfActives,(geoDim*(geoDim+1)/2),nrOfActives) = pointPhyHessian;
            m_phyGrad.block(0,i*nrOfActives,geoDim,nrOfActives) = pointPhyGrad;
        }
    }
  */

protected:
    gsMatrix<T> m_BasisGrad;
    gsMatrix<T> m_basisHessian;
    gsMatrix<T> m_phyHessian;
    gsMatrix<T> m_phyGrad;
    gsMatrix<T> m_Z_phyGrad;
    gsMatrix<T> m_zCoors;
    gsMatrix<T> m_zPhyHessian;
    gsMatrix<T> m_zCoorEle;
    
    using base::coeff_A_vals;
    using base::actives;
    using base::rhsVals;
    using base::basisData;
    using base::physBasisGrad;

    using base::numActive;

    using base::coeff_b_vals;
    using base::coeff_c_vals;
    using base::rhs_ptr;
    using base::localMat;
    using base::localRhs;
    using base::md;
};


template <class T>
class gsAiryStress: public gsCDRAssembler<T>
{
public:
    typedef gsCDRAssembler<T> Base;

    public:

    //gsAiryStress() { }
                  
    /** @brief
        Constructor of the object.
    */
    gsAiryStress( gsMultiPatch<T> const         & patches,
                    gsBoundaryConditions<T> const & bconditions,
                    const gsFunction<T>           & rhs,
                    const std::vector<gsMatrix<> > & initialZCoor,                    
                    const gsPointLoads<T>         & pLoads  = gsPointLoads<T>()         
                    )
                    : Base( patches,
                    gsMultiBasis<T>(patches),
                    bconditions,
                    rhs,
                    m_artifiFun,
                    coeffB,
                    coeffC)
        {
            gsVector<T> tem_vec9(4);
            tem_vec9 << 0,0,0,0;
            rhsFun = &rhs;
            m_artifiFun.setValue(tem_vec9, 2);
            tem_vec9.setZero(2);
            coeffB.setValue(tem_vec9,2);
            coeffC.setValue(0,2);

            m_pLoads = pLoads;
            m_zCoor = initialZCoor;
        }

public:


//     void CalStress( const   gsMultiPatch<T> & m_curSolution, 
//         std::vector<gsMatrix<T> > & stressPatches )
//     {    }

      /// Main assembly routine
    void assemble()
    {
        if (0 == this->numDofs() ) // Are there any interior dofs ?
        {
            gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
            return;
        }

        m_system.reserve(m_bases[0], m_options, this->pde().numRhs());
        this->computeDirichletDofs();

        // Pre-allocate non-zero elements for each column of the
        // sparse matrix
        //m_system.matrix().setZero();

        //this->constructSolution(m_zCoor, m_zCoorPatch);
//         m_zCoorPatch = m_patches;
//         m_zCoorPatch.patch(0).embed_Nd(3);
//         m_zCoorPatch.patch(0).coefs().rightCols(1)= m_zCoor[0].leftCols(1);
//         gsWriteParaview(m_zCoorPatch, "temp_paraview", 1000);

        gsMatrix<T> zPatchCoor;
        // Assemble volume stiffness and load vector integrals
        for (size_t np=0; np < m_pde_ptr->domain().nPatches(); ++np )
        {
            //Assemble stiffness matrix and rhs for the local patch
            // with index np and add to m_matrix and m_rhs
            //zPatchCoor = m_zCoorPatch.patch(np).coefs();
            
            gsVisitorAiryStress<T> visitor( *rhsFun, m_artifiFun, coeffB, coeffC, m_zCoor[np]);
            this->apply(visitor, np);
        }
        
        // Enforce Neumann boundary conditions
        Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );

        const int dirStr = m_options.getInt("DirichletStrategy");
        // If requested, enforce Dirichlet boundary conditions by Nitsche's method
        if ( dirStr == dirichlet::nitsche )
            Base::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());
        // If requested, enforce Dirichlet boundary conditions by diagonal penalization
        else if ( dirStr == dirichlet::penalize )
            Base::penalizeDirichletDofs();

        // Assembly is done, compress the matrix
        Base::finalize();
    }

    void applyLoads()
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
            m_system.colMapper(0).localToGlobal(acts, m_pLoads[i].patch, acts);

            gsMatrix<T> & m_rhs = m_system.rhs();
            for (index_t k=0; k < acts.rows(); ++k)
            {
                m_rhs(acts(k,0), 0) += bVals(k,0) * m_pLoads[i].value[2];
            }
        }
    }

protected:
    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;


public:
    gsConstantFunction<T> coeff_b;  
    gsConstantFunction<T> coeff_c; 
    gsConstantFunction<T> m_artifiFun;  

public:
    const gsFunction<T> *coeff_ptr;
    const gsFunction<T> * rhsFun;
    gsConstantFunction<T> coeffB;
    gsConstantFunction<T> coeffC;
    gsPointLoads<T>  m_pLoads;   
    std::vector<gsMatrix<> > m_zCoor;
    gsMultiPatch<T> m_zCoorPatch;

};


}//namespace gismo

#include <gsSelfSuppSurf/gsAiryStress.hpp>
