/** @file gsFittingLinSyst.h

    @brief Contains the class gsFittingLinSyst that permits to perform linear fitting

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsFitting/gsFittingParam.h>

#include <gsFitting/gsFittingIdConstr.h>
#include <gsFitting/gsFittingQuadrature.h>

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
class _gsFittingLinSyst
{

protected:
    /// Initialize all the different attributes of the class
    void init_gsFittingLinSyst(GenBasis &basis, bool remove_result,
                           unsigned dim_im, gsFittingParam<T>& param);
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
    explicit _gsFittingLinSyst(unsigned dim_im, GenBasis & basis,
                               bool remove_result,
                               gsFittingParam<T>& param,
                               gsLeastSquares<GenBasis, T>*
                               least_squares = NULL);

    /// Sets the least squares method
    void set_least_squares(gsLeastSquares<GenBasis, T>*
                           least_squares)
    {   m_least_squares = least_squares;  }

    /// Sets the quadrature of dimension 1.
    /// Generally used in the 2D case for fitting some boundary
    void setQuadrature
    (gsFittingQuadrature <1,GenBasis,2,T>& quadr)
    {
        m_quadr_1d = new gsFittingQuadrature <1,GenBasis,2,T>(quadr);
    }

    void decreaseLMult(T val)
    {
        if(m_least_squares != NULL)
            m_least_squares->decreaseLMult(val);
    }

    /// Sets the quadrature of dimension 2
    /// Generally used:
    ///  - in the 3D case for fitting some boundary
    ///  - in the 2D case for smoothing energies
    void setQuadrature
    (gsFittingQuadrature <2,GenBasis,2,T>& quadr)
    {
        m_quadr_2d = new gsFittingQuadrature <2,GenBasis,2,T>(quadr);
    }

    /// Sets the quadrature of dimension 3 (can only be used
    /// in the 3D case). Generally used for smoothing energies.
    void setQuadrature
    (gsFittingQuadrature <3,GenBasis,2,T>& quadr)
    {
        m_quadr_3d = new gsFittingQuadrature <3,GenBasis,2,T>(quadr);
    }

    gsFittingQuadrature<1, GenBasis, 2, T>*
    quadr_1d() {   return m_quadr_1d;   }
    gsFittingQuadrature<2, GenBasis, 2, T>*
    quadr_2d() {   return m_quadr_2d;   }
    gsFittingQuadrature<3, GenBasis, 2, T>*
    quadr_3d() {   return m_quadr_3d;   }

    /// Add the integrand integr to the set of integrands
    /// (contained in the set of quadrature).
    /// TODO: if we consider multiple quadrature of the same dimension,
    ///       the domain should be checked (find the quadrature that
    ///       has the whole domain)
    void addGlobalIntegrand(gsFittingIntegrand<d,T>* integr);

    /// Used in the case where the template mapping (initial mapping) is used
    /// and not given: construct this mapping
    void construct_template_mapping();


    /// Destructor
    virtual ~_gsFittingLinSyst()
    {
        if ( m_remove_result && m_result_linear != NULL)
            delete m_result_linear;
    }

public:

    /// Computes the minimization. The result is set in m_result_linear.
    void computeLinear();
    void compute();

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
    void actualize_identity_mapping();

    /// Return true if the energy has some smoothing component
    bool has_smoothing();

    /// Adds to the matrix A_mat terms for the continuous
    /// component of the energy, i.e., the quadrature
    void assembleContinuous(gsSparseEntries<T> & A_mat,
                            gsMatrix<T>& m_B);

    /// Assembles system for the least square fit.
    void assembleSystem(gsSparseEntries<T>& A_mat,
                        gsMatrix<T>& B);

public:
    /// Returns the computed approximation
    virtual gsFunctionSet<T>* result() const
    { return m_result_linear; }

    /// Returns the basis of the approximation
    GenBasis & getBasis() {  return *m_basis;  }

    /// Returns the template (identity mapping)
    gsFunctionSet<T>* get_template(){ return m_template; }

    /// Sets the template (identity map) as the current map.
    void template_as_current_map()
    {  this->set_current_map(m_template); }

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

    /// Compute all the energies (errors + smoothing energies)
    void computeErrors(int type=1);

    int dim_im(){    return m_dim;  }

    /// Print the energies
    void print_error(bool compute = true);

    /// Reset the result to NULL
    void reset()
    {  m_result_linear = NULL;  }

    /// Return true if we can split the dimensions
    /// during the computation
    bool canSplitDimension(){  return m_split_dim;  }

    /// Return true if the energy is quadratic
    bool isLinear(){  return m_isLinear;  }

    /// In the case where the dimensions cannot be considered
    /// separatly, the global index has to be modified.
    /// Take this into account if needed.
    int add_dimension_to_index(int global_index);

    /// Check if the integrands are linear and if they are
    /// used for smoothing of for fitting data. Put the results
    /// in the appropriate attributes.
    void set_properties_integrands();

    /// Multiply the coefficient of every integrand of type
    /// smoothing by ratio.
    void decrease_coeff_smoothing(T ratio);

    /// Multiply the coefficient of every integrand
    /// of quadratic energies of type smoothing by ratio.
    void decrease_linear_coeff_smoothing(T ratio);

    /// Returns the least squares pointor.
    gsLeastSquares<GenBasis, T>* leastSquares()
    { return m_least_squares;  }

    /// Adds a Tikhonov regularization term.
    /// Currently used for the projection onto the constraints
    void add_regularization_term(T coeff);

    /// Construct the geometry of the result from the coefficients
    /// obtained by solving the linear system.
    void makeGeometry(gsMatrix<T>& coefs);

    /// In the case of multipatch, actualize the mapper if
    /// the basis have been modified
    void actualize_multipatch(bool repair)
    {  m_local_global.actualize_mapper(repair);  }

    /// Return the global index from the index in some patch.
    inline int localToGlobal(int ind, int patch, bool split_dim)
    {  return m_local_global.
            localToGlobal(ind, patch, split_dim);  }

    /// Returns the number of degrees of freedom
    inline int sizeDof(){ return m_local_global.sizeDof(); }

protected:
    /// Calls the mapper when necessary (multipatch)
    /// to obtain a global index from a local index
    gsLocalGlobal<GenBasis, T> m_local_global;

    /// The path of the output
    std::string m_output;


    /// The dimension of the image
    /// (dimension of the space where the target domain is embedded)
    unsigned m_dim;

    /// If true, we solve the system independently on each dimension
    bool m_split_dim;

    /// If true, all the energies are quadratic
    bool m_isLinear;

    /// Pointer keeping the basis generating the space on which
    /// the minimization is performed
    GenBasis * m_basis;

    /// Performs the least squares methods on some points.
    /// NULL is the user did not give any.
    gsLeastSquares<GenBasis, T>* m_least_squares;

    /// Constructs a template mapping in case
    /// it is not given when calling compute()
    gsFittingIdConstr<d, GenBasis, T> m_template_constr;

    /// <---------- The quadrature part of the energy ----------->

    /// Quadrature domains can have dimension up to d.

    /// Quadrature of dimension 1. Generally used in the 2D
    /// case for fitting some boundary
    gsFittingQuadrature<1, GenBasis, 2, T>* m_quadr_1d;

    /// Quadrature of dimension 2. Generally used:
    ///  - in the 3D case for fitting some boundary
    ///  - in the 2D case for smoothing energies
    gsFittingQuadrature<2, GenBasis, 2, T>* m_quadr_2d;

    /// Quadrature of dimension 3 (can only be used in the 3D case).
    /// Generally used for smoothing energies or L2 projection.
    gsFittingQuadrature<3, GenBasis, 2, T>* m_quadr_3d;

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
    gsFunctionSet<T>* m_template;

    /// Is the template given at the construction
    bool m_template_given;

    /// True if we minimize the deformation
    /// (only used in the case where there is an energy
    ///  depending on the metric)
    bool m_deform_minim;

    /// Should we print messages during the process
    bool m_print_messages;

    /// Is the system linearized (at each step, we compute
    /// a displacement) or not?
    bool m_is_displacement;

    /// The matrix associated with the linear system.
    /// Used in compute linear (and possibly for updating
    /// the Lagrange multipliers)
    gsSparseMatrix<T> m_A_mat;

    /// The maximal smoothing energy value
    T m_energy_max;

}; /// class _gsFittingLinSyst



template<short_t d, class GenBasis, class T = real_t>
class gsFittingLinSyst : public _gsFittingLinSyst<d, GenBasis, T>
{
public:
    typedef _gsFittingLinSyst<d, GenBasis, T> Base;

    /// Constructors
    explicit gsFittingLinSyst(unsigned dim_im, GenBasis & basis,
                          bool remove_result,
                          gsFittingParam<T>& param,
                          gsLeastSquares<GenBasis, T>*
                              least_squares = NULL);

public:

    /// The number of different bases, i.e., the number of patches.
    inline int nBases();


    using Base::m_least_squares;
    using Base::m_basis;
    using Base::m_result_linear;
};   /// class gsFittingLinSyst


template<short_t d, class T>
class gsFittingLinSyst<d, gsBasis<T>, T> :
        public _gsFittingLinSyst<d, gsBasis<T>, T>
{
public:
    typedef _gsFittingLinSyst<d, gsBasis<T>, T>  Base;


    /// Constructors
    explicit gsFittingLinSyst(unsigned dim_im, gsBasis<T> & basis,
                          bool remove_result,
                          gsFittingParam<T>& param,
                          gsLeastSquares<gsBasis<T>, T>*
                              least_squares = NULL);


public:
    /// The number of different bases, i.e., the number of patches.
    inline int nBases(){ return 1; }

protected:
    using Base::m_basis;
    using Base::m_least_squares;
    using Base::m_template;
    using Base::m_result_linear;
    using Base::m_dim;
    using Base::m_current;

    using Base::m_quadr_1d;
    using Base::m_quadr_2d;
    using Base::m_quadr_3d;
};  /// class gsFittingLinSyst


template<short_t d, class T>
class gsFittingLinSyst<d, gsMultiBasis<T>, T>
    : public _gsFittingLinSyst<d, gsMultiBasis<T>, T>
{
public:
    typedef _gsFittingLinSyst<d, gsMultiBasis<T>, T> Base;

    /// Constructors
    explicit gsFittingLinSyst(unsigned dim_im,
                              gsMultiBasis<T> & basis,
                              bool remove_result,
                              gsFittingParam<T>& param,
                              gsLeastSquares<gsMultiBasis<T>, T>*
                              least_squares = NULL);


    /// The number of different bases, i.e., the number of patches.
    inline int nBases(){  return m_basis->nBases();  }

protected:
    using Base::m_result_linear;
    using Base::m_basis;
    using Base::m_template;
    using Base::m_least_squares;
    using Base::m_dim;
    using Base::m_current;

    using Base::m_quadr_1d;
    using Base::m_quadr_2d;
    using Base::m_quadr_3d;

/*private:
  gsDofMapper m_mapper; */
};  /// class gsFittingLinSyst

}// namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFittingLinSyst.hpp)
#endif
