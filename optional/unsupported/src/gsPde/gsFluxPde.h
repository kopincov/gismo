/** @file gsFluxPde.h

    @brief Describes a gernic Flux PDE.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/


#pragma once

#include <gsPde/gsPdeWithCoeff.h>
#include <gsCore/gsPiecewiseFunction.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

/** @brief
    A PDE of the form

    \f$-\text{div}(F(\nabla\mathbf{u}))=\mathbf{f} \f$,

    with a given nonlinear Flux F(\nabla u).

    This class describes a this PDE, with an arbitrary right-hand side
    function.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T> class gsFlux;

template<class T>
class gsFluxPde : public gsPde<T>, public gsPdeWithCoeff<T>
{

public:

    gsFluxPde( ) { }


    /// Constructor
    gsFluxPde(const gsMultiPatch<T>         &domain,
                 const gsBoundaryConditions<T> &bc,
                 const gsFlux<T>            *flux,
                 const gsPiecewiseFunction<T>  &rhs)
    : gsPde<T>(domain,bc), m_rhs(rhs)
    {
        m_flux.resize(domain.nPatches());
        for(size_t np=0; np<domain.nPatches();np++)
        {
            m_flux[np]= flux->clone().release();
            m_flux[np]->setPatch(np);
        }
        m_unknownDim.setOnes(1);
    }

    ~gsFluxPde()
    {
        freeAll(m_flux);
    }

    /**
     * @brief gives the number of rhs functions of the PDEs
     */
    virtual int numRhs() const
    {
        return m_rhs.piece(0).targetDim();
    }

    void setCurSolution(const gsMultiPatch<T>& curSol) const
    {
        for(size_t np=0; np<m_domain.nPatches();np++)
            m_flux[np]->setCurSolution(curSol);
    }

    const gsFunction<T> *    rhs()      const { return &m_rhs.piece(0); }

    const std::vector<gsFlux<T>* > & getFlux() const {return m_flux;}
    
    virtual int numUnknowns() const     {return 1;}

    virtual bool isSymmetric () const { gsWarn<<"Function is gsPde::isSymmetric should not be used!!"; return true;}

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
       //TODO:
        return os;
    }


//this may be not threadsave .. depends on gsFlux
    virtual gsFluxPde<T>* restrictToPatch(unsigned np) const
    {
        gsBoundaryConditions<T> bc;
        m_boundary_conditions.getConditionsForPatch(np,bc);
        return new gsFluxPde<T>(m_domain.patch(np),bc,m_flux.front(), m_rhs);
    }
//change to gsMatrix<T>
    virtual T getCoeffForIETI(unsigned np) const {
        return m_flux[np]->getCoeffForIETI();
    }
protected:
    using gsPde<T>::m_unknownDim;
    using gsPde<T>::m_domain;
    using gsPde<T>::m_boundary_conditions;

    std::vector<gsFlux<T>* > m_flux;

    gsPiecewiseFunction<T> m_rhs;
}; // class gsFluxPde


/**
 * Interface for a generic nonlinear Flux on a single patch.
 * The instance must implement the function:
 * -isLinear()
 * -eval_into()
 * -deriv_into()
 *
 *  isLinear() should return if the considered flux is linear (for optimization)
 *  eval_into() evaluated the flux at given quadrature points
 *  deriv_into() evaluated the linearization of the flux in given quadrature points
 *
 *  Before the latter two are called, setCurSolution() and setGeometryEvaluator()
 *  must be called, since the flux needs the current solution (nonlinear dependence) and
 *  the Geometry Evaluator for the physical quadrature points and the transformation
 *  of the gradients
 *
 */
template<class T>
class gsFlux : public gsFunctionSet<T>
{
public:
    /// Shared pointer for gsFlux
    typedef memory::shared_ptr< gsFlux > Ptr;

    /// Unique pointer for gsFlux
    typedef memory::unique_ptr< gsFlux > uPtr;

    index_t size() const {return 1;}

    short_t domainDim() const {return m_targetDim;}

    short_t targetDim() const {return m_targetDim;}

public: // --- Do before calling eval and deriv ----//

    //is called at the construction of the FluxPDE
    void setPatch(int patch) {m_patch= patch;}

    void setCurSolution(const gsMultiPatch<T> & curSolution) const
    {m_curSolution = &curSolution.patch(m_patch);}

    void setMapData(const gsMapData<T> & md) { m_md = md;}

public:
    //Access the current geometry for a patch
    gsGeometry<T>& getCurSolution() {return *m_curSolution;}


public:// --- Abstract Methods ----//
    GISMO_UPTR_FUNCTION_PURE(gsFlux, clone)

    //TODO: this might be removed.--- returns true if Flux is linear
    virtual bool isLinear() const = 0;

    virtual void eval_into(const gsMatrix<T> &u, gsMatrix<T> &result) = 0;

    virtual void deriv_into(const gsMatrix<T> &u, gsMatrix<T> &result) = 0;

//    virtual void eval_into(const gsMapData<T> &md, gsMatrix<T> &result) = 0;
//
//    virtual void deriv_into(const gsMapData<T> &md, gsMatrix<T> &result) = 0;

public: // --- optional non implement Methods ----//

    //implement if you want to use coefficient scaling in the IETI method
    T getCoeffForIETI() const {GISMO_NO_IMPLEMENTATION}

protected:
     mutable const gsGeometry<T> * m_curSolution;

     int m_targetDim;

     int m_patch;

     gsMapData<T>  m_md;
};


/**
 * A simple example for a Flux:
 * \f$ F(\nabla u):= \alpha(x)\nabla u \f$. This leads just to the Poisson Equation.
 *  So the nonlinear Newton solver must converge in 1 step.
 */
template<class T>
class gsPoissonFlux : public gsFlux<T>
{
public:
    /// Shared pointer for gsPoissonFlux
    typedef memory::shared_ptr< gsPoissonFlux > Ptr;

    /// Unique pointer for gsPoissonFlux
    typedef memory::unique_ptr< gsPoissonFlux > uPtr;

    gsPoissonFlux(short_t targetDim, const gsPiecewiseFunction<T>& alpha)
        :m_alpha(alpha)
    {
        m_targetDim =targetDim;
    }

    GISMO_CLONE_FUNCTION(gsPoissonFlux)

    void deriv_into(const gsMatrix<T> &u, gsMatrix<T> &result)
    {
        int d = m_targetDim;
        gsMatrix<T> alphaVals;
        m_alpha.piece(m_patch).eval_into(m_md.values[0],alphaVals);
        result.setZero(d*d,u.cols());

        for(int c = 0;c!= d;++c)
        {
            for(int j=0;j!=d;++j)
                if(j==c)
                    result.row(c*d+j)=alphaVals;
        }
    }

    void eval_into(const gsMatrix<T> &u, gsMatrix<T> &result)
    {
        gsMatrix<T> alphaVals, grads, transGrad;
        result.setZero(u.rows(), u.cols());

        m_curSolution->deriv_into(u,grads);
        for(int k=0; k!=u.cols();++k)
        {
            transformGradients(m_md, k,grads,transGrad);
            result.col(k) = transGrad.col(0);
        }

        m_alpha.piece(m_patch).eval_into(m_md.values[0],alphaVals);

        for(index_t p = 0; p!=u.cols(); ++p )
            result.col(p)= result.col(p)*T(alphaVals(0,p));
    }

//    void deriv_into(const gsMapData<T> &md, gsMatrix<T> &result)
//    {
//        int d = m_targetDim;
//        gsMatrix<T> alphaVals;
//        m_alpha.piece(m_patch).eval_into(md.values[0], alphaVals);
//        result.setZero(d * d, md.points.cols());
//
//        for (int c = 0; c != d; ++c)
//        {
//            for (int j = 0; j != d; ++j)
//                if (j == c)
//                    result.row(c * d + j) = alphaVals;
//        }
//    }
//
//    void eval_into(const gsMapData<T> &md, gsMatrix<T> &result)
//    {
//        gsMatrix<T> alphaVals, grads, transGrad;
//        result.setZero(md.points.rows(), md.points.cols());
//
//        m_curSolution->deriv_into(md.points, grads);
//        for (int k = 0; k != md.points.cols(); ++k)
//        {
//            transformGradients(md, k, grads, transGrad);
//            result.col(k) = transGrad.col(0);
//        }
//
//        m_alpha.piece(m_patch).eval_into(md.values[0], alphaVals);
//
//        for (index_t p = 0; p != md.points.cols(); ++p)
//            result.col(p) = result.col(p) * T(alphaVals(0, p));
//    }

    bool isLinear() const {return true;}

    T getCoeffForIETI(unsigned np) const
    {
        gsMatrix<T> result;
        //m_targetdim should be replaced here by m_pardim
        m_alpha.piece(np).eval_into(gsMatrix<T>::Zero(m_targetDim,1),result);
        return result(0,0);
    }

protected:
    gsPiecewiseFunction<T> m_alpha;

    using gsFlux<T>::m_patch;
    using gsFlux<T>::m_targetDim;
    using gsFlux<T>::m_curSolution;
    using gsFlux<T>::m_md;
};

template<class T>
class gspLaplaceFlux : public gsFlux<T>
{
public:
    /// Shared pointer for gspLaplaceFlux
    typedef memory::shared_ptr< gspLaplaceFlux > Ptr;

    /// Unique pointer for gspLaplaceFlux
    typedef memory::unique_ptr< gspLaplaceFlux > uPtr;

    gspLaplaceFlux(short_t targetDim, T eps, T p)
        :m_eps(eps), m_p(p)
    {
        m_targetDim =targetDim;
    }

    GISMO_CLONE_FUNCTION(gspLaplaceFlux)

    void deriv_into(const gsMatrix<T> &u, gsMatrix<T> &result)
    {
        gsMatrix<T> grads, transGrad;
        int d = m_targetDim;
        result.setZero(d*d,u.cols());

        m_curSolution->deriv_into(u,grads);
        for(int k=0; k!=u.cols();++k)
        {
            transformGradients(m_md, k,grads,transGrad);
            result(0,k)= (m_p-1)*transGrad(0,0)*transGrad(0,0) + transGrad(1,0)*transGrad(1,0)+ m_eps*m_eps;
            result(1,k)= (m_p -2)*transGrad(0,0)*transGrad(1,0);
            result(2,k)= result(1,k);
            result(3,k)= (m_p-1)*transGrad(1,0)*transGrad(1,0) + transGrad(0,0)*transGrad(0,0)+ m_eps*m_eps;

            result.col(k)*=math::pow((m_eps*m_eps + transGrad.col(0).squaredNorm()),(m_p-4)/2);
        }
    }

    void eval_into(const gsMatrix<T> &u, gsMatrix<T> &result)
    {
        gsMatrix<T> grads, transGrad;
        result.setZero(u.rows(), u.cols());

        m_curSolution->deriv_into(u,grads);
        for(int k=0; k!=u.cols();++k)
        {
            transformGradients(m_md, k,grads,transGrad);
            T val = math::pow((m_eps*m_eps+ transGrad.col(0).squaredNorm()),(m_p-2)/2);
            result.col(k) = transGrad.col(0)*val;
        }
    }

//    void deriv_into(const gsMapData<T> &md, gsMatrix<T> &result)
//    {
//        gsMatrix<T> grads, transGrad;
//        int d = m_targetDim;
//        result.setZero(d * d, md.points.cols());
//
//        m_curSolution->deriv_into(md.points, grads);
//        for (int k = 0; k != md.points.cols(); ++k)
//        {
//            transformGradients(md, k, grads, transGrad);
//            result(0, k) =
//                (m_p - 1) * transGrad(0, 0) * transGrad(0, 0) + transGrad(1, 0) * transGrad(1, 0) + m_eps * m_eps;
//            result(1, k) = (m_p - 2) * transGrad(0, 0) * transGrad(1, 0);
//            result(2, k) = result(1, k);
//            result(3, k) =
//                (m_p - 1) * transGrad(1, 0) * transGrad(1, 0) + transGrad(0, 0) * transGrad(0, 0) + m_eps * m_eps;
//
//            result.col(k) *= math::pow((m_eps * m_eps + transGrad.col(0).squaredNorm()), (m_p - 4) / 2);
//        }
//    }
//
//    void eval_into(const gsMapData<T> &md, gsMatrix<T> &result)
//    {
//        gsMatrix<T> grads, transGrad;
//        result.setZero(md.points.rows(), md.points.cols());
//
//        m_curSolution->deriv_into(md.points, grads);
//        for (int k = 0; k != md.points.cols(); ++k)
//        {
//            transformGradients(md, k, grads, transGrad);
//            T val = math::pow((m_eps * m_eps + transGrad.col(0).squaredNorm()), (m_p - 2) / 2);
//            result.col(k) = transGrad.col(0) * val;
//        }
//    }

    bool isLinear() const {return false;}

    T getCoeffForIETI(unsigned np) const
    {
        return math::pow((m_eps*m_eps),(m_p-2)/2);
    }

protected:
    T m_eps;
    T m_p;

    using gsFlux<T>::m_patch;
    using gsFlux<T>::m_targetDim;
    using gsFlux<T>::m_curSolution;
    using gsFlux<T>::m_md;
};

} // namespace gismo
