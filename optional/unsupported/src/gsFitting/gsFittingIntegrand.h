/** @file gsFittingIntegrand.h

    @brief Provides the virtual gsFittingIntegrand class and some linear integrands classes inheriting gsFittingIntegrand.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsGeometry.h>

#include <gsFitting/gsFittingIterator.h>

namespace gismo
{


/**
   @brief Virtual class representing the integrand of an energy
   that is a fonctionnal depending on some mapping u
   (to be constructed).
   A simple example may be the square norm of the gradient of u.
   The domain is traveled in gsFittingQuadrature and then each
   integrand is called to fill the matrix and the RHS. The integrand
   can also be called to compute the value of the energy for a given
   function u (a solution for example).
   For each elementary cell of the domain, the points of the
   quadrature and the values (image + derivatives) of the active basis
   functions at those points must be set calling setLocalValues.
   Once these values are set, the virtual function actualize()
   is called.
   Then, the elements to be added to the matrix and to the RHS
   can be computed by calling the virtual functions compute_matrix
   and compute_RHS.

   In the case there is at least one nonlinear component in the energy,
   the system must be linearized near some position (say u_0), which means that
   the solution is a "small displacement" of this current position.
   In that case m_have_current_mapping is true and the elements to be
   added to the matrix and to the RHS depend on this current position
   (even in the linear case). In this case, the values and
   derivatives of the current position have to be set using
   setLocalValuesCurr. Once these values are set, the virtual
   function actualizeCurr is called. The elements to be added to the
   matrix and to the RHS are then obtained in a same fashion
   (nothing has to be changed for that part).

   In a same manner, to compute the energy, we call at each
   elementary cell the function setLocalValuesError to set the
   points used for the quadrature and values of the current solution
   at these points. Then, the virtual function actualize_energy is
   called. The function compute_energy can finally be called to
   compute the energy at these points.

   Note that the quadrature can be computed on a set different
   from the domain on which is defined the basis.
   Then, to the quadrature points can be associated parameters called
   points_param_space. This is for example used for quadrature on
   spaces of lower dimension (boundary for example). If the mapping
   is the identity (the quadrature is performed on the same space)
   then points_param_space are simply set to NULL and not used.

   o: the maximal order needed for the computation of the
      derivatives (o between 0 and 2)
   d: the dimension of the parameter space
   (equals to the dimension of the domain of the basis functions
   in the case there is no parameter space)
   \ingroup
**/
template<short_t d, class T=real_t>
class gsFittingIntegrand
{

public:
    /// Constructor: initialize the properties of the integrand
    ///  (linear or not, smoothing or fitting type, can we split
    ///   the dimensions in the computation if we have this energy)
    gsFittingIntegrand(bool is_smoothing, bool split_dim,
                       bool isLinear);

    /// Destructor
    virtual ~gsFittingIntegrand(){}

    /// Compute the value to be added to the matrix
    /// for the point ind_pt and the basis (base_i, base_j),
    /// respectively primal and dual.
    /// Note that ind_pt is the position of the point in the
    /// set of points given when calling setLocalValues. Also,
    /// base_i and base_j are the positions of the bases in the set
    /// of active basis given in a same manner.
    /// In the case where the dimensions are not split,
    /// the value also depends on the dimension considered for the
    /// primal element (base_i) and the dual element (base_j)
    /// respectively given by d1 and d2.
    /// If this energy can split dimension, the value
    /// 0 should be returned when d1 != d2.
    /// Note that if the dimensions are split, we suppose that
    /// no energy depend on the dimension and this
    /// function is only called once with d1 = d2 = 0
    virtual T compute_matrix(int base_i, int base_j, int ind_pt,
                             int d1, int d2)
    { GISMO_NO_IMPLEMENTATION }

    /// Compute the value to be added to the RHS at the position
    /// for the point ind_pt and the dual basis base_j. See
    /// compute_matrix for more details on these integers.
    /// The computation is set in res.
    virtual void compute_RHS(int base_j, int ind_pt,
                             gsMatrix<T>& res)
    { GISMO_NO_IMPLEMENTATION }

    /// Returns the value of the energy at the point ind_pt
    /// (see compute_matrix for more details on this value).
    virtual T compute_energy(int ind_pt){ GISMO_NO_IMPLEMENTATION }


    /********* Informations about the integrand   ********/
    /// Returns true if, according to this energy, we can split
    /// the dimensions.
    bool canSplitDimension(){   return m_split_dim;   }

    /// Returns true if this energy is linear.
    bool isLinear() {   return m_isLinear;   }

    /// Returns true if this energy is of type smoothing.
    bool is_smoothing(){  return m_is_smoothing;  }

    /// If the energy is of smoothing type, multiply the coefficient by ratio
    virtual void decrease_coeff_smoothing(T ratio){  }
    virtual void decrease_linear_coeff_smoothing(T ratio){  }


/******** functions called when entering in a cell:
        * set all the values that will be used      ***********/

    /// Sets the following data used for the computation of the
    /// RHS and of the matrix:
    /// - points of the quadrature
    /// - the points associated to them in the parameter space
    /// - the values for each basis function active in this cell
    /// - the number of basis functions actives in this cell
    /// - the dimension of the geometry (that can be greater
    ///   or equal to the dimension d of the domain of the basis)
    /// - the first and second derivatives of the basis functions
    /// - the current patch (unused in the single patch case)
    /// - have_current_mapping which is true in the case where
    ///   the current system is obtained by linearization
    void setLocalValues(gsMatrix<T>* points,
                        gsMatrix<T>* points_param_space,
                        gsMatrix<T>* image,
                        unsigned nActives,
                        unsigned dim_im,
                        gsMatrix<T>* der1 = NULL,
                        gsMatrix<T>* der2 = NULL,
                        unsigned patch = 0,
                        bool have_current_mapping = false);


    /// Sets the following data used for the computation of the energy
    /// - points of the quadrature
    /// - the points associated to them in the parameter space
    /// - the values of the current solution at these points
    /// - the first and second derivatives of the current solution
    /// - the current patch (unused in the single patch case)
    void setLocalValuesError(unsigned m_dim, gsMatrix<T>* points,
                             gsMatrix<T>* points_param_space,
                             gsMatrix<T>* image,
                             gsMatrix<T>* der1 = NULL,
                             gsMatrix<T>* der2 = NULL,
                             unsigned patch = 0);

    /// In the case where the system is linearized,
    /// sets the value of the current position (solution)
    /// as well as it's first and second derivatives at
    /// the quadrature points
    void setLocalValuesCurr(gsMatrix<T>* image,
                            gsMatrix<T>* der1,
                            gsMatrix<T>* der2);

    /// In the case where the energy depends on the metric,
    /// the optimization can be performed on the deformation, i.e.,
    /// the difference between the mapping and some identity mapping.
    /// This identity mapping can be set after the construction of
    /// the object by implementing this virtual function.
    /// This function is called once the
    /// identity mapping is computed.
    virtual void
    set_identity_mapping(gsFunctionSet<T>* identity_map){  }

    /// Must be implemented. These actualize... functions may be
    /// used to compute some values relative to the energy
    /// before calling compute_matrix, compute_RHS
    /// and compute_energy.

    /// Called in the end setLocalValues
    /// (where we set the data related to the basis functions)
    virtual void actualize(){ GISMO_NO_IMPLEMENTATION }

    /// Called in the end setLocalValuesCurr
    /// (where we set the data related to the current solution)
    virtual void actualizeCurr(){  }

    /// Called in the end setLocalValuesError
    /// (where we set the data related to the current solution
    ///  before computing the energy)
    virtual void actualize_energy(){ GISMO_NO_IMPLEMENTATION }

    /// Copies the current integrand
    virtual gsFittingIntegrand<d, T>* copy()
    { GISMO_NO_IMPLEMENTATION }

protected:
    /// We denote in the following:
    ///  - Np the number of points used for the quadrature

    /// The points (and their position in the parameter space, if any)
    /// in the current elementary cell.
    gsMatrix<T>* m_points;  /// SIZE: m_dim X Np
    gsMatrix<T>* m_points_param_space;   /// SIZE: d X Np


    /// The images and derivatives of the active basis functions
    /// at the positions m_points
    gsMatrix<T>* m_image; /// SIZE: (m_nActives, Np)
    gsMatrix<T>* m_der1;  /// SIZE: (m_nActives * m_stride[0], Np)
    gsMatrix<T>* m_der2;  /// SIZE: (m_nActives * m_stride[1], Np)


    /// The images and derivatives of the current solution
    /// at the positions m_points
    gsMatrix<T>* m_image_curr;  /// SIZE: (m_dim, Np)
    gsMatrix<T>* m_der1_curr; /// SIZE: (m_dim * m_stride[0], Np)
    gsMatrix<T>* m_der2_curr; /// SIZE: (m_dim * m_stride[1], Np)

    /// The patch of the current elementary cell
    /// (only used in case of multipatches)
    unsigned m_patch;

    /// The number of basis that are active
    unsigned m_nActives;

    /// True if the current mapping is used
    bool m_have_current_mapping;

    /// Stride for the first and the second derivatives.
    /// Equals to {d, d*(d-1)/2}
    unsigned m_stride[2];

    /// True if this energy does not prevent from splitting the
    /// dimensions in the computation
    bool m_split_dim;

    /// True if this integrand is of type smoothing.
    /// False if it is associated with a fitting
    bool m_is_smoothing;

    /// True if the energy is quadratic
    bool m_isLinear;

    /// The dimension of the image
    unsigned m_dim;

}; /// gsFittingIntegrand


/**
   @brief A class containing the "identity". Used by integrands
   having energies depending on the deformation
*/
template<unsigned o = 1, class T=real_t>
class gsFittingId
{
public:
    /// Constructor
    gsFittingId(gsFunctionSet<T>* map = NULL)
    {   m_map = map;  }

    /// Sets the reference to the identity mapping
    void set_identity_mapping(gsFunctionSet<T>* map)
    { m_map = map; }

    /// Actualize: compute the value and derivatives of the identity
    /// at the quadrature points.
    void actualize(int patch, gsMatrix<T>& points);

    /// Return false if the identity is used.
    inline bool is_null(){  return m_map == NULL; }

    /// Return the values of the identity at the points of
    /// the current elementary cell.
    inline gsMatrix<T>& image_id()
    {  return m_image_id;   }

    /// Return the derivatives of the identity at the points of
    /// the current elementary cell.
    inline gsMatrix<T>& der1_id()
    {  return m_der1_id;   }

protected:
    /// A reference to the identity mapping
    gsFunctionSet<T>* m_map;

    /// The values and the derivatives of the identity mapping
    /// at the points of the current elementary cell.
    gsMatrix<T> m_image_id;
    gsMatrix<T> m_der1_id;

}; /// gsFittingID

/*
  @brief: Class used for representing the integrand of the energy
  given by the squared norm of the first and second derivatives.
  This class is also used for Tikhonov regularization.
  In that case, m_isLinear must be set to false.
  In the linear case, the energy depends on u_0 + h, with u_0 the current
  solution and h the displacement searched.
  In the Tikhonov regularization case, the energy only depends on h, the displacement searched.
*/
template<short_t d, unsigned o = 2, class T=real_t>
class gsFittingIntegrandLin : public gsFittingIntegrand<d, T>
{
public:
    /// Constructor.
    /// isLinear is true when we minimize this energy.
    /// isLinear is false when this class is used
    /// for Tikhonov regularization.
    gsFittingIntegrandLin(T coeff_linear, T coeff_grad,
                          T coeff_hess, bool isLinear) :
    gsFittingIntegrand<d, T>::
    gsFittingIntegrand(true, true, isLinear),
    m_identity_map(NULL)
    {
        m_coeff_linear = coeff_linear;
        m_coeff_grad = coeff_grad * m_coeff_linear;
        m_coeff_hess = coeff_hess * m_coeff_linear;
    }

    /// Resets the global linear coefficient
    void setCoeffGlobal(T coeff_global)
    {
        m_coeff_grad *= coeff_global / m_coeff_linear;
        m_coeff_hess *= coeff_global / m_coeff_linear;
        m_coeff_linear = coeff_global;
    }

    void decrease_linear_coeff_smoothing(T ratio)
    {
        decrease_coeff_smoothing(ratio);
    }

    void decrease_coeff_smoothing(T ratio)
    {
        m_coeff_grad *= ratio;
        m_coeff_hess *= ratio;
        m_coeff_linear *= ratio;
    }

private:
    typedef gsFittingIntegrand<d, T> Base;

    using Base::m_stride;
    using Base::m_image;
    using Base::m_points;
    using Base::m_patch;
    using Base::m_der1;
    using Base::m_der2;
    using Base::m_der1_curr;
    using Base::m_der2_curr;
    using Base::m_have_current_mapping;
    using Base::m_dim;
    using Base::m_isLinear;

    /// Importance of the first and second derivatives in the energy
    T m_coeff_grad;
    T m_coeff_hess;
    T m_coeff_linear;

    /// In the case where we minimize on the gradient
    /// minus the gradient of the identity
    gsFittingId<1, T> m_identity_map;

    /// Computes the value to be added to the matrix
    /// for the gradient energy
    T computeGrad(int base_i, int base_j, int ind_pt);

    /// Computes the square norm of the gradient
    /// of the current solution
    T energyGrad(int ind_pt);

    /// Computes the value to be added to the matrix
    /// for the Hessian energy
    T computeHess(int base_i, int base_j, int ind_pt);

    /// Computes the square norm of the Hessian
    /// of the current solution
    T energyHess(int ind_pt);

public:
    T compute_matrix(int base_i, int base_j, int ind_pt,
                     int d1, int d2);

    /// Computes the value to be added to the RHS
    /// for the gradient energy
    void RHS_grad(int base_i, int ind_pt, gsMatrix<T>& res);

    /// Computes the value to be added to the RHS
    /// for the Hessian energy
    void RHS_hess(int base_i, int ind_pt, gsMatrix<T>& res);

    inline void compute_RHS(int base_i, int ind_pt,
                            gsMatrix<T>& res)
    {
        res.setZero(m_dim, 1);
        if(m_have_current_mapping && m_isLinear)
        {
            RHS_grad(base_i, ind_pt, res);
            RHS_hess(base_i, ind_pt, res);
        } else if(! m_identity_map.is_null())
            RHS_grad_id(base_i, ind_pt, res);
    }

    /// Computes the RHS for the gradient energy in the
    /// case where there is no current mapping
    /// (this energy is given by the Id mapping).
    void RHS_grad_id(int base_i, int ind_pt,
                     gsMatrix<T>& res);


    inline T compute_energy(int ind_pt)
    {
        if(m_isLinear)
           return (m_coeff_grad*energyGrad(ind_pt)
                   + m_coeff_hess*energyHess(ind_pt))
           / m_coeff_linear;
        else
            return 0.;

    }

    /// We compute the values and derivatives of the identity
    void actualize()
    {  m_identity_map.actualize(m_patch, *m_points);  }

    void actualize_energy()
    { actualize(); }

    void set_identity_mapping(gsFunctionSet<T>* map)
    {
        if(m_isLinear)
            m_identity_map.set_identity_mapping(map);
    }

    gsFittingIntegrandLin<d, o, T>* copy();

}; /// gsFittingIntegrandLin


/// The integrand associated with the L^2 distance minimization
template<short_t d, class GenGeom, class T=real_t>
class gsFittingIntegrandL2Dist : public gsFittingIntegrand<d, T>
{
public:
    typedef gsFittingIntegrand<d, T> Base;

    gsFittingIntegrandL2Dist(int dimAmbient, GenGeom *f,
                             bool isLinear = true,
                             T coeff = 1.) :
    gsFittingIntegrand<d, T>::gsFittingIntegrand(false, true,
                                                 isLinear)
    {
        m_f = f;
        m_dimAmbient = dimAmbient;
        m_coeff = coeff;
    }
    void decrease_linear_coeff_smoothing(T ratio)
    {
        decrease_coeff_smoothing(ratio);
    }

    void decrease_coeff_smoothing(T ratio)
    {
        m_coeff *= ratio;
    }

private:
    using Base::m_image;
    using Base::m_image_curr;
    using Base::m_patch;
    using Base::m_points;
    using Base::m_points_param_space;
    using Base::m_have_current_mapping;
    using Base::m_dim;

    /// The target function
    GenGeom *m_f;

public:
    T compute_matrix(int base_i, int base_j, int ind_pt,
                     int d1, int d2);
    void compute_RHS(int base_i, int ind_pt, gsMatrix<T>& res);

    T compute_energy(int ind_pt);

    /// Computes the values and derivatives of the function m_f
    /// that must be fitted
    void actualize();
    void actualize_energy(){  actualize();  }

    gsFittingIntegrandL2Dist<d, GenGeom, T>* copy();

private:
    gsMatrix<T> m_image_f;

    int m_dimAmbient;

    T m_coeff;

}; /// gsFittingIntegrandL2Dist


/// The integrand associated with the L^2 norm minimization
/// (only used to compute the L2 norm, the value of the matrix
///  and of the RHS is not implemented)
template<short_t d, class T=real_t>
class gsFittingIntegrandL2Norm : public gsFittingIntegrand<d, T>
{
public:
    typedef gsFittingIntegrand<d, T> Base;

    gsFittingIntegrandL2Norm() :
    gsFittingIntegrand<d, T>::gsFittingIntegrand(false, true, true)
    {   }

private:
    using Base::m_image;
    using Base::m_image_curr;
    using Base::m_patch;
    using Base::m_points;
    using Base::m_points_param_space;
    using Base::m_have_current_mapping;
    using Base::m_dim;

public:
    T compute_matrix(int base_i, int base_j, int ind_pt,
                     int d1, int d2){ return 0.; }
    void compute_RHS(int base_i, int ind_pt, gsMatrix<T>& res){  }

    T compute_energy(int ind_pt);

    void actualize(){  }
    void actualize_energy(){   }

    gsFittingIntegrandL2Norm<d, T>* copy()
    {  return new gsFittingIntegrandL2Norm<d, T>();   }

}; /// gsFittingIntegrandL2Norm

}// namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFittingIntegrand.hpp)
#endif

//////////////////////////////////////////////////
//////////////////////////////////////////////////
