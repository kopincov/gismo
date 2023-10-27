
#pragma once

#include <gsPde/gsPde.h>

namespace gismo
{


template <class T> class gsFunction;

/** @brief
    A stationary linear elasticity PDE.

    This class describes a stationary linear elasticity PDE, with an arbitrary right-hand side
    function.

    In the 2D case, the plain strain problem is considered.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsLinearElasticityPde : public gsPde<T>
{
protected:
    gsLinearElasticityPde( ) { }
    using gsPde<T>::m_domain;
    using gsPde<T>::m_unknownDim;

public:
    gsLinearElasticityPde(
            const gsMultiPatch<T>          &domain,
            const gsBoundaryConditions<T>  &bc,
            gsFunction<T>       *force,
            const T             YoungsModulus = 10000,
            const T             PoissonsRatio = 0.3,
            gsFunction<T>       *source = NULL
        )
        :
            gsPde<T>(domain,bc),
            m_YoungsModulus(YoungsModulus),
            m_PoissonsRatio( PoissonsRatio)
    {
        m_force  = force  ? force->clone().release()  : NULL;
        m_source = source ? source->clone().release() : NULL;

        m_unknownDim.setConstant(1,m_domain.dim());

        m_lambda = m_YoungsModulus * m_PoissonsRatio / ( (1.+m_PoissonsRatio)*(1.-2.*m_PoissonsRatio)) ;
        m_mu     = m_YoungsModulus / (2.*(1.+m_PoissonsRatio)) ;
    }

    ~gsLinearElasticityPde( )
    {
        delete m_force;
        delete m_source;
    }

    const gsFunction<T>* rhs() const
    { return m_force; }
    const gsFunction<T>* force() const
    { return m_force; }
    const gsFunction<T>* source() const
    { return m_source; }

    T lambda() const   { return m_lambda; }
    T mu() const       { return m_mu; }


    // / Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os<<"Linear elasticity equation:\n"
          <<"-div sigma(u) = f,\n"
          <<" u=gD "
          <<"with:\n";
        if ( m_force )
        os<<"Force  function f= "<< *m_force <<".\n";
        return os;
    }
    /// Consistency check
    bool check()
    {
        return true;
    }
protected:
    const gsFunction<T> * m_force;
    const gsFunction<T> * m_source;

public:
    T m_YoungsModulus;
    T m_PoissonsRatio;
    T m_lambda;
    T m_mu;
}; // class gsLinearElasticityPde

} // namespace gismo
