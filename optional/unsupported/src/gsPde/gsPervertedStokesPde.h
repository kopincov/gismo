
#pragma once

#include <gsPde/gsPde.h>

namespace gismo
{


template <class T> class gsFunction;

/** @brief
    A perverted Stokes PDE.

    This class describes a perverted Stokes PDE, with an arbitrary right-hand side
    function and optionally a known solution.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsPervertedStokesPde : public gsPde<T>
{
protected:
    gsPervertedStokesPde( ) { }
    using gsPde<T>::m_domain;
    using gsPde<T>::m_unknownDim;

public:
    gsPervertedStokesPde(
        const gsMultiPatch<T>          &domain,
        const gsBoundaryConditions<T>  &bc,
        gsFunction<T>       *force,
        gsFunction<T>       *source = NULL,
        const T                    viscosity = 1
        )
        :
            gsPde<T>(domain,bc),    m_viscosity(viscosity)

    {
        m_force  = force ? force->clone().release() : NULL;
        m_source = source ? source->clone().release() : NULL;

        m_unknownDim.resize(4);
        m_unknownDim[0] = m_domain.dim();
        m_unknownDim[1] = 1;
        m_unknownDim[2] = m_domain.dim();
        m_unknownDim[3] = 1;
    }

    ~gsPervertedStokesPde( )
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

    T viscocity() const                 { return m_viscosity; }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os<<"Perverted Stokes's :\n"
          <<"ask Jarle Sogn!\n"
          <<"with:\n";
        if ( m_force )
        os<<"Force  function f= "<< *m_force <<".\n";
        if ( m_source )
        os<<"Source function g= "<< *m_source <<".\n";
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

    T m_viscosity;
}; // class gsPervertedStokesPde

} // namespace gismo
