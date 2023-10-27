/** @file gsFittingBase.h

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

#include <gsCore/gsLinearAlgebra.h>
#include <gsMatrix/gsMatrix.h>

#include <gsFitting/gsFittingParam.h>

#include <gsFitting/gsFittingIdConstr.h>
#include <gsFitting/gsFittingEnergy.h>
#include <gsFitting/gsFittingUtilsGen.h>

namespace gismo
{


/**
   @brief
   Class for performing the minimization of a given energy that can contain fitting elements
   and smoothing elements. The fitting component of the energy can be of two types:
    - a gsLeastSquares element that contain a certain number of points (and associated parameters) that must be fitted
    - a gsFittingIntegrand element (contained in a gsFittingQuadrature element) for the case where a L2 distance to a given Geometry is minimized.
    The smoothing component is given by one or multiple gsFittingIntegrand (contained in a gsFittingQuadrature element, usually of dimension d).
    Smoothing energies given in gsFittingIntegrand (inside a gsFittingQuadrature) can be either linear or nonlinear.
    Only linear systems are solved, although nonlinear energies can be given.
    These nonlinear energies are linearized near some current position m_current.

    The basis can be:
    - a gsBasis, for single patch optimization;
    - a gsMultiBasis for multipatch optimization.
   GenBasis is the type of basis used (only gsBasis<T> and gsMultiBasis<T> are used)
   This class is also used for optimization performed iteratively (see gsFittingIter)
**/
template<short_t d, class GenBasis, class T=real_t>
class gsFittingBase
{

protected:
    /// Initialize all the different attributes of the class
    void init_gsFittingBase(bool remove_result, unsigned dim_im,
                            gsFittingParam<T>& param);
    /// Set the current map (used for the linearization of the
    /// nonlinear energies)
    void set_current_map(gsFunctionSet<T>* current);

public:

    /// Constructors

    /// Set only the least squares method. The smoothing
    /// (if one is used) has to be added manually
    /// before the computation.
    /// If no energy is given when this constructor is called,
    /// these must be added manually before calling computeLinear
    explicit gsFittingBase(unsigned dim_im, GenBasis & basis,
                           bool remove_result,
                           gsFittingParam<T>& param,
                           const gsPointContainer<T>& pts_LS);
    explicit gsFittingBase(unsigned dim_im, GenBasis & basis,
                           bool remove_result,
                           gsFittingParam<T>& param);


    /// Used in the case where the template mapping
    /// (initial mapping) is used
    /// and not given: construct this mapping
    void constructId();

    /// Destructor
    virtual ~gsFittingBase()
    {
        if ( m_remove_result && m_result_linear != NULL)
            delete m_result_linear;
        /*   if(m_bc != NULL)
        {
               /// We free all the geometries associated
            /// to the boundaries
            for ( typename gsBoundaryConditions<T>::const_iterator
                      iter = m_bc->dirichletBegin();
                  iter != m_bc->dirichletEnd(); ++iter )
            {
                gsGeometry<T>* geom = static_cast<gsGeometry<T>*>
                    (iter->function().get());
                delete geom;
                }

            delete m_bc;
        }  */
    }
    //  gsBoundaryConditions<T>* get_bc() {  return *m_bc;  }

public:

    /// Computes the minimization. The result is set in m_result_linear.
    void computeLinear();
    virtual void compute();

    /// If the basis has been modified, most likely in the case
    /// of hierarchical fitting, actualise the data
    /// depending on this basis. In the current implementation,
    /// we only resample points of the least squares method,
    /// if needed, take into account these additional
    /// degrees of freedom
    void actualize_basis();

    /// The identity mapping has been modified (or simply set).
    /// Actualize this modification on the energies using it,
    /// in the case where this identity is used.
    /// Most probably, these energies are metric related and we
    /// minimize the metric associated with the deformation
    void actualize_identity_mapping()
    {   quadrEner().actualize_identity_mapping(m_identity);   }

    /// Adds to the matrix A_mat terms for the continuous
    /// component of the energy, i.e., the quadrature
    void assembleContinuous();

public:
    /// Returns the computed approximation
    virtual gsFunctionSet<T>* result() const
    { return m_result_linear; }

    /// Returns the basis of the approximation
    GenBasis & basis() {  return m_basis;  }

    void computeErrors(int type=1)
    {
        GISMO_ASSERT(m_current != NULL,
                     "Error: there is no solution");
        m_energy.computeErrors(type);
    }

    /// Returns the template (identity mapping)
    gsFunctionSet<T>* get_template(){ return m_identity; }

    /// Sets the template (identity map) as the current map.
    void template_as_current_map()
    {  this->set_current_map(m_identity); }

    /// Print the energies
    void print_error(bool compute = true)
    {   m_energy.print_error(compute);  }

    /// Reset the result to NULL
    void reset()
    {  m_result_linear = NULL;  }

    /// Return true if the energy has some smoothing component
    bool has_smoothing()
    {   return quadrEner().has_smoothing();   }


    /// Multiply the coefficient of every integrand of type
    /// smoothing by ratio.
    void decrease_coeff_smoothing(T ratio)
    {  quadrEner().decrease_coeff_smoothing(ratio);  }

    /// Multiply the coefficient of every integrand
    /// of quadratic energies of type smoothing by ratio.
    void decrease_linear_coeff_smoothing(T ratio)
    {  quadrEner().
            decrease_linear_coeff_smoothing(ratio);  }

    void decreaseLMult(T coeff_reduc)
    {  m_energy.decreaseLMult(coeff_reduc);  }

    /// Adds a Tikhonov regularization term.
    /// Currently used for the projection onto the constraints
    void add_regularization_term(T coeff)
    {   quadrEner().template
            add_regularization_term<d>(coeff, m_dim, m_basis);   }

    gsLeastSquares<GenBasis, T>& leastSquares()
    {  return m_energy.leastSquares();  }

    /*  void add_fixed_border(std::vector<gsMatrix<T> >& pts,
                          std::vector<gsMatrix<T> >& param,
                          std::vector<int>& ind);

    void add_fixed_border(gsMultiPatch<T>& border,
                          gsMultiPatch<T>& temp_border,
                          std::vector<int>& ind);*/

    /// Returns the dimension of the image
    int dim_im(){    return m_dim;  }

    gsFittingEnergy<GenBasis, T>& energy()
    {   return m_energy;    }

    inline gsFittingQuadrEnergy<GenBasis, T>& quadrEner()
    {   return m_energy.quadrEner();    }

    inline T maxEnergySmoothing()
    {  return quadrEner().maxEnergySmoothing();   }

    inline bool isLinear()
    { return m_energy.isLinear();  }

    inline bool canSplitDimension()
    { return m_energy.canSplitDimension(); }

    void set_properties_integrands();

    void resetPointsLS(std::vector<gsMatrix<T> >& pts,
                       std::vector<gsMatrix<T> >& param)
    {  m_energy.resetPointsLS(pts, param);  }

protected:

    gsFittingEnergy<GenBasis, T> m_energy;

    gsFittingSystem<GenBasis, T> m_system;

    /*  /// The Dirichlet boundary condition, if any
        gsBoundaryConditions<T>* m_bc;  */

    /// The path of the output
    std::string m_output;

    /// The dimension of the image
    /// (dimension of the space where the target domain is embedded)
    unsigned m_dim;

    /// Pointer keeping the basis generating the space on which
    /// the minimization is performed
    GenBasis& m_basis;

    /// Constructs a template mapping in case
    /// it is not given when calling compute()
    gsFittingIdConstr<d, GenBasis, T> m_identity_constr;

    /// The geometry resulting from the computation.
    /// It can be:
    /// - of type gsGeometry if we consider a single patch
    /// - of type gsMultiPatch if we consider a multipatch
    gsFunctionSet<T>* m_result_linear;

    /// The current mapping. We linearise the nonlinear
    /// energies near this mapping for the minimization
    gsFunctionSet<T>* m_current;

    /// True if we free the result when we destruct the object
    bool m_remove_result;

    /// Template domain: defines the identity mapping in case
    ///   the optimization is operated on the deformation
    gsFunctionSet<T>* m_identity;

    /// Is the template given at the construction
    bool m_identity_given;

    /// Should we print messages during the process
    bool m_print_messages;

    /// Is the system linearized (at each step, we compute
    /// a displacement) or not?
    bool m_is_displacement;

    /// True if we minimize the deformation
    /// (only used in the case where there is an energy
    ///  depending on the metric)
    bool m_deform_minim;


}; /// class gsFittingBase


}// namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFittingBase.hpp)
#endif
