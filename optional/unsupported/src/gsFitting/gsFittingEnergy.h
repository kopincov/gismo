/** @file gsFittingEnergy.h

    @brief Contains the class gsFittingBase that permits to perform linear fitting.
    TODO use gsAssembler and gsPde?

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsMatrix/gsMatrix.h>
#include <gsFitting/gsFittingParam.h>

#include <gsFitting/gsFittingIdConstr.h>
#include <gsFitting/gsFittingQuadrature.h>
#include <gsFitting/gsLeastSquares.h>
#include <gsFitting/gsLeastSquaresLMult.h>

namespace gismo
{

template<class GenBasis, class T=real_t>
class gsFittingQuadrEnergy
{

    template <unsigned d>
    struct ContainerQuadr
    {
        typedef std::vector< gsFittingQuadrature<d, GenBasis, 2, T> > type;
        typedef std::vector< gsFittingQuadrature<d, GenBasis,
                                                 2, T> >* Ptype;
    };

public:
    gsFittingQuadrEnergy(bool deform_min)
    : m_1d(), m_2d(), m_3d()
    {   m_deform_minim = deform_min;  }

    /// If true, we solve the system independently on each dimension
    bool m_split_dim;

    /// If true, all the energies are quadratic
    bool m_isLinear;

    /// Quadrature domains can have dimension up to d.

    /// Quadrature of dimension 1. Generally used in the 2D
    /// case for fitting some boundary
    typename ContainerQuadr<1>::type m_1d;

    /// Quadrature of dimension 2. Generally used:
    ///  - in the 3D case for fitting some boundary
    ///  - in the 2D case for smoothing energies
    typename ContainerQuadr<2>::type m_2d;

    /// Quadrature of dimension 3 (can only be used in the 3D case).
    /// Generally used for smoothing energies or L2 projection.
    typename ContainerQuadr<3>::type m_3d;

    /// True if we minimize the deformation
    /// (only used in the case where there is an energy
    ///  depending on the metric)
    bool m_deform_minim;

public:

    gsFittingQuadrEnergy();

    /// Returns the maximum of the errors
    T maxError();
    /// Returns the minimum of the errors
    T minError();
    /// Returns the sum of the errors
    T totalError();

    /// Return the maximum of the smoothing energies
    T maxEnergySmoothing();
    /// Return the minimum of the smoothing energies
    T minEnergySmoothing();
    /// Return the sum of the smoothing energies
    T totalEnergySmoothing();

    /// Check if the integrands are linear and if they are
    /// used for smoothing of for fitting data. Put the results
    /// in the appropriate attributes.
    void set_properties_integrands();

    /// Return true if the energy has some smoothing component
    bool has_smoothing();

    /// Add the integrand integr to the set of integrands
    /// (contained in the set of quadrature).
    /// TODO: if we consider multiple quadrature of the same dimension,
    ///       the domain should be checked (find the quadrature that
    ///       has the whole domain)
    template<unsigned d>
    void addGlobalIntegrand(gsFittingIntegrand<d, T>* integr,
                            index_t dim_im, GenBasis& basis);

    void assemble(gsFittingSystem<GenBasis, T>& system);

    void set_current_map(gsFunctionSet<T>* current);

    /// Sets the quadrature of dimension 1.
    /// Generally used in the 2D case for fitting some boundary
    template<unsigned d> void addQuadrature
    (gsFittingQuadrature <d,GenBasis,2,T>& quadr)
    {   m_1d.push_back(quadr);  }

    /// Sets the quadrature of dimension 2
    /// Generally used:
    ///  - in the 3D case for fitting some boundary
    ///  - in the 2D case for smoothing energies
    void addQuadrature
    (gsFittingQuadrature <2,GenBasis,2,T>& quadr)
    {  m_2d.push_back(quadr);  }

    /// Sets the quadrature of dimension 3 (can only be used
    /// in the 3D case). Generally used for smoothing energies.
    void addQuadrature
    (gsFittingQuadrature <3,GenBasis,2,T>& quadr)
    {  m_3d.push_back(quadr);   }


    void copy_non_smoothing(gsFittingQuadrEnergy<GenBasis, T>
                            & new_ener);
    template<unsigned d> static void
    copy_non_smoothing(typename ContainerQuadr<d>::type& origin,
                       typename ContainerQuadr<d>::type& copy);

    /**
       Return the list of quadrature of a given dimension
       Any other way to do this without pointor?
     */

private:
    inline void _getQuadratures(typename ContainerQuadr<1>::Ptype& res)
    {   res = &m_1d;   }
    inline void _getQuadratures(typename ContainerQuadr<2>::Ptype& res)
    {   res = &m_2d;   }
    inline void _getQuadratures(typename ContainerQuadr<3>::Ptype& res)
    {   res = &m_3d;   }

public:
    template<unsigned d> typename
    ContainerQuadr<d>::type&
    getQuadratures()
    {
        typename ContainerQuadr<d>::Ptype Pres;
        _getQuadratures(Pres);
        return *Pres;
    }

    index_t sizeQuadrature(index_t d)
    {
        GISMO_ASSERT(d < 4 && d > 0, "Wrong index");
        if(d == 1)
            return m_1d.size();
        else if(d == 2)
            return m_2d.size();
        else if(d == 3)
            return m_3d.size();
        else
            return -1;
    }

    gsFittingQuadratureGen<GenBasis, 2, T>&
    getQuadrature(index_t d, index_t ind)
    {
        GISMO_ASSERT(ind >= 0  && ind < sizeQuadrature(d),
                     "Wrong index");
        if(d == 1)
            return m_1d[ind];
        else if(d == 2)
            return m_2d[ind];
        else
            return m_3d[ind];
    }


    /// Multiply the coefficient of every integrand of type
    /// smoothing by ratio.
    void decrease_coeff_smoothing(T ratio);

    /// Multiply the coefficient of every integrand
    /// of quadratic energies of type smoothing by ratio.
    void decrease_linear_coeff_smoothing(T ratio);

    /// Adds a Tikhonov regularization term.
    /// Currently used for the projection onto the constraints
    template<unsigned d> void
    add_regularization_term(T coeff, index_t dim, GenBasis& basis);

    /// Compute all the energies (errors + smoothing energies)
    void computeErrors(index_t type=1);

    /// The identity mapping has been modified (or simply set).
    /// Actualize this modification on the energies using it,
    /// in the case where this identity is used.
    /// Most probably, these energies are metric related and we
    /// minimize the metric associated with the deformation
    void actualize_identity_mapping(gsFunctionSet<T>* ident);

};  /// gsFittingQuadrEnergy



/**
   @brief
*/
template<class GenBasis, class T=real_t>
class gsFittingEnergy
{
protected:

    /// Performs the least squares methods on some points.
    /// NULL is the user did not give any.
    gsLeastSquaresLMult<GenBasis, T> m_least_squares;

    gsFittingQuadrEnergy<GenBasis, T> m_quadrEner;

    /// The dimension of the image
    /// (dimension of the space where the target domain is embedded)
    index_t m_dim;

    /// Should we print messages during the process
    bool m_print_messages;

public:
    gsFittingEnergy(index_t dim, GenBasis& basis,
                    gsFittingParam<T>& param,
                    const gsPointContainer<T>& pts_LS
                    = gsPointContainer<T>())
    : m_least_squares(basis, pts_LS, param.output,
                      param.export_points, param.coeff_lagrangeM),
      m_quadrEner(param.deform_min)
    {
        m_dim = dim;
        m_print_messages = param.print_messages;
    }

    /// Returns the maximum of the errors
    T maxError()
    {
        T maxError = 0.;
        if(! m_least_squares.isEmpty())
            maxError = m_least_squares.maxError();
        return std::max(m_quadrEner.maxError(), maxError);
    }

    /// Returns the minimum of the errors
    T minError()
    {
        T minError = m_quadrEner.minError();
        if(m_least_squares.isEmpty())
            return minError;
        else
            return std::min(minError,
                            m_least_squares.minError() );
    }
    /// Returns the sum of the errors
    T totalError()
    {
        T totalError = 0.;
        if(! m_least_squares.isEmpty())
            totalError = m_least_squares.totalError();
        return m_quadrEner.totalError() +
            totalError;
    }

    void decreaseLMult(T val)
    {
        if(! m_least_squares.isEmpty())
            m_least_squares.decreaseLMult(val);
    }

    /// Return true if we can split the dimensions
    /// during the computation
    bool canSplitDimension()
    {  return m_quadrEner.m_split_dim;  }

    /// Sets the least squares method
    void resetPointsLS(std::vector<gsMatrix<T> >& pts,
                       std::vector<gsMatrix<T> >& param)
    {   m_least_squares.resetPointsLS(pts, param);  }
    void resetPointsLS(gsMatrix<T>& pts, gsMatrix<T>& param)
    {   m_least_squares.resetPointsLS(pts, param);  }

    void resetPointsLS(gsPointContainer<T>& pts)
    {   m_least_squares.resetPointsLS(pts);  }

    void addPointsLS(std::vector< gsMatrix<T> >& pts,
                     std::vector< gsMatrix<T> >& param)
    {   m_least_squares.addPointsLS(pts, param);  }
    void addPointsLS(gsMatrix<T>& pts, gsMatrix<T>& param,
                     index_t ind)
    {   m_least_squares.addPointsLS(pts, param, ind);  }

    void addPointsLS(gsPointContainer<T>& pts)
    {   m_least_squares.addPointsLS(pts);  }

    /// Compute all the energies (errors + smoothing energies)
    void computeErrors(index_t type=1)
    {
        if(! m_least_squares.isEmpty())
            m_least_squares.computeErrors(type);
        quadrEner().computeErrors(type);
    }

    /// Return true if the energy has some smoothing component
    bool has_smoothing()
    {  return m_quadrEner.has_smoothing();  }

    /// Return true if the energy is quadratic
    bool isLinear(){  return m_quadrEner.m_isLinear;  }

    void assemble(gsFittingSystem<GenBasis, T>& system)
    {
        if(! m_least_squares.isEmpty())
            m_least_squares.assemble(system);
        m_quadrEner.assemble(system);
    }

    void set_current_map(gsFunctionSet<T>* current)
    {
        if(! m_least_squares.isEmpty())
            m_least_squares.set_current_map(current);
        m_quadrEner.set_current_map(current);
    }

    /// Returns the energy obtained by quadrature
    inline gsFittingQuadrEnergy<GenBasis, T>& quadrEner()
    {   return m_quadrEner;  }

    /// Returns the dimension of the image
    index_t dim_im(){    return m_dim;  }

    void copy_non_smoothing(gsFittingEnergy<GenBasis, T>& new_ener)
    {
        if(! m_least_squares.isEmpty())
            new_ener.m_least_squares.resetPointsLS
                (m_least_squares);
        m_quadrEner.copy_non_smoothing(new_ener.m_quadrEner);
    }

    void print_error(bool compute);

    void actualize_basis();

    gsLeastSquares<GenBasis, T>& leastSquares()
    {  return m_least_squares;  }

    void adaptMultipliers(gsFunctionSet<T>* sol,
                          gsFunctionSet<T>* proj
                          /*, gsSparseMatrix<T>& mat*/ )
    {
        if(! m_least_squares.isEmpty())
            m_least_squares.adaptMultipliers(sol, proj);
    }

}; /// gsFittingEnergy

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFittingEnergy.hpp)
#endif
