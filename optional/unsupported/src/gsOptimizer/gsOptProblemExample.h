/** @file gsOptProblemExample.h

    @brief Example setting up an optimization problem

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsIpopt/gsOptProblem.h>

namespace gismo
{


/** 
 * @brief 
 * Simple optimization example, to demonstrate the definition of an
 * optimization problem using the base class gsOptProblem.
 *
 *  This class implements the following NLP.
 *
 * min_x f(x) = -(x1-2)^2      (objetive function)
 *  s.t.
 *       0 = x0^2 + x1 - 1     (constraint)
 *       -1 <= x0 <= 1         (variable bounds)
 *
 */

template <typename T>
class gsOptProblemExample : public gsOptProblem<T>
{
public:

    gsOptProblemExample()
    {
        m_numDesignVars  = 2;
        m_numConstraints = 1;
        m_numConJacNonZero = 2;

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);
        
        // x0 has a lower bound of -1 and an upper bound of 1
        m_desLowerBounds[0] = -1.0;
        m_desUpperBounds[0] =  1.0;
        
        // x1 has no upper or lower bound, so we set them to
        // a large negative and a large positive number.
        // The value that is interpretted as -/+infinity can be
        // set in the options, but it defaults to -/+1e19
        m_desLowerBounds[1] = -1.0e19;
        m_desUpperBounds[1] =  1.0e19;

        m_conLowerBounds.resize(m_numConstraints);
        m_conUpperBounds.resize(m_numConstraints);

        // we have one equality constraint, so we set the bounds on
        // this constraint to be equal (and zero).
        m_conLowerBounds[0] = 
        m_conUpperBounds[0] = 0;
        
        // we initialize x in bounds, in the upper right quadrant
        m_curDesign.resize(2,1);
        m_curDesign(0,0) = 0.5;
        m_curDesign(1,0) = 1.5;        
        
        // 
        m_conJacRows.resize(m_numConJacNonZero);
        m_conJacCols.resize(m_numConJacNonZero);

        // element at 0,0: grad_{x0} g_{1}(x)
        m_conJacRows[0] = 0;
        m_conJacCols[0] = 0;
        // element at 0,1: grad_{x0} g_{1}(x)
        m_conJacRows[1] = 0;
        m_conJacCols[1] = 1;
    }


public:

    T evalObj( const gsAsConstVector<T> & u ) const
    {
        // return the value of the objective function
        const T x1 = u(1,0);
        return -(x1 - 2.0) * (x1 - 2.0);
    }
    
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const  
    {
        result.resize(m_numDesignVars,1);

        // grad_{x0} f(x): x0 is not in the objective
        result(0,0)  = 0.0;

        // grad_{x1} f(x):
        const T x1 = u(1,0);
        result(1,0)  =   -2.0*(x1 - 2.0);
    }
    
    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        // return the value of the constraints: g(x)
        result.resize(m_numConstraints,1);

        const T x0 = u(0,0);
        const T x1 = u(1,0);
        result[0]  = -(x0*x0 + x1 - 1.0);
    }
    
    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        result.resize(m_numDesignVars,1);

        const T x0 = u(0,0);

        // element at 0,0: grad_{x0} g_{1}(x)
        result[0] = -2.0 * x0;
        
        // element at 0,1: grad_{x1} g_{1}(x)
        result[1] = -1.0;
    }

private:

    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;
};

} // end namespace gismo
