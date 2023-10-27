/** @file gsFittingIntegrandNL.h

    @brief Contains the integrands of nonlinear energies

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsFitting/gsFittingIterator.h>

namespace gismo
{


/**
   Virtual class inheriting gsFittingIntegrand. This class
   represents the integrand of nonlinear energies
*/
template<short_t d, class T=real_t>
class gsFittingIntegrandNL : public gsFittingIntegrand<d, T>
{

public:
    gsFittingIntegrandNL(T coeff_NL) :
    gsFittingIntegrand<d, T>::
    gsFittingIntegrand(true, false, false)
    {  m_coeff_NL = coeff_NL;}
    virtual ~gsFittingIntegrandNL(){}

    virtual gsFittingIntegrand<d, T>* copy()
    { GISMO_NO_IMPLEMENTATION }

    /// Returns the signs of the permutation
    /// that transforms {1,2,3} to {i,j,k}
    inline static int sgn(int i, int j, int k)
    {
        if(i+1 == j || j+1 == k)
            return 1;
        else return -1;
    }
protected:
    /// The importance of the NL energies
    T m_coeff_NL;

}; /// gsFittingIntegrandNL


/**
   Class representing the integrand of the energy defined
   by the square norm of the metric defined by the
   mapping (or the deformation)
*/
template<short_t d, class T=real_t>
class gsFittingIntegrandOrth : public gsFittingIntegrandNL<d, T>
{

    typedef gsFittingIntegrand<d, T> Base;
    typedef gsFittingIntegrandNL<d, T> BaseNL;
    using BaseNL::m_coeff_NL;

public:
    /// Constructor
    gsFittingIntegrandOrth(T coeff_NL, T coeff_orth) :
    BaseNL::gsFittingIntegrandNL(coeff_NL)
    {  m_coeff_orth = coeff_orth * m_coeff_NL;  }

    /// Destructor
    virtual ~gsFittingIntegrandOrth(){}

    gsFittingIntegrand<d, T>* copy()
    { return new gsFittingIntegrandOrth
            (m_coeff_NL, m_coeff_orth/m_coeff_NL);  }

    /// Computes some values (metric of the current
    /// solution or deformation)
    /// used for computing the matrix and the RHS
    void actualizeCurr();

    /// Computes some values (metric of the current
    /// solution or deformation)
    /// used for computing the matrix and the RHS
    void set_identity_mapping(gsFunctionSet<T>* map)
    {  m_identity_map.set_identity_mapping(map);  }

    /// Computes the image and derivatives of the identity
    void actualize()
    {  m_identity_map.actualize(m_patch, *m_points);  }

    /// Computes some values (metric of the current
    /// solution or deformation) used for computing the energy
    void actualize_energy()
    {
        m_identity_map.actualize(m_patch, *m_points);
        computeMetric();
    }

    T compute_matrix(int base_i, int base_j, int ind_pt,
                     int d1, int d2);
    void compute_RHS(int base_i, int ind_pt,
                     gsMatrix<T>& res);

    /// Returns the index associated with the basis and the
    /// derivative in the matrix m_mixed_metric
    /// der_ij: the index of the expected second derivative
    inline unsigned getPosMixed(int ind_basis, int der_ij)
    {  return ind_basis*m_dim*(d-1)*d/2 + der_ij*m_dim;  }

    /// Returns the index associated with the basis and the
    /// derivative in the matrix m_mixed_metric_diag
    inline unsigned getPosMixedDiag(int ind_basis, int der)
    {  return ind_basis*m_dim*m_stride[0] + der*m_dim;  }

    /// Computes the metric of the current solution
    void computeMetric();


    T compute_energy(int ind_pt);

    void decrease_coeff_smoothing(T ratio)
    {
        m_coeff_NL *= ratio;
        m_coeff_orth *= ratio;
    }

protected:
    using Base::m_stride;
    using Base::m_der1;
    using Base::m_der1_curr;
    using Base::m_nActives;
    using Base::m_points;
    using Base::m_dim;
    using Base::m_patch;

    /// The metric of the current solution (or deformation)
    /// at each point of the current elementary cell.
    /// Diagonal of the metric
    gsMatrix<T> m_metric_diag; /// SIZE: (d, nPoints)
    /// Other elements of the metric
    gsMatrix<T> m_metric; /// SIZE: (d*(d-1)/2, nPoints)

    /// Some values used to compute the derivatives of the energy
    gsMatrix<T> m_mixed_metric; /// SIZE: (d*(d-1)/2 * nBasis * d, nPoints)
    gsMatrix<T> m_mixed_metric_diag;  /// SIZE: (d* nBasis * d, nPoints)

    /// In the case where we minimize the gradient
    /// minus the gradient of the identity
    gsFittingId<1, T> m_identity_map;

    /// The proportion of this energy in the global nonlinear energy
    T m_coeff_orth;

}; /// gsFittingIntegrandOrth


/**
   Class representing the integrand of the Winslow energy
   in 2D.
*/
template<class T=real_t>
class gsFittingIntegrandW : public gsFittingIntegrandNL<2, T>
{

    typedef gsFittingIntegrand<2, T> Base;
    typedef gsFittingIntegrandNL<2, T> BaseNL;

public:
    /// Constructor
    gsFittingIntegrandW(T coeff_NL, T coeff_winslow) :
    BaseNL::gsFittingIntegrandNL(coeff_NL)
    {  m_coeff_winslow = coeff_winslow * m_coeff_NL;  }

    /// Destructor
    virtual ~gsFittingIntegrandW(){}

    gsFittingIntegrand<2, T>* copy()
    { return new gsFittingIntegrandW
            (m_coeff_NL, m_coeff_winslow/m_coeff_NL);  }

    void actualizeCurr();

    /// Computes the normals and their norms at each
    /// point of the current elementary cell
    void actualize_energy(){  setNormals(); }
    void actualize(){  }

    /// Computes the normals and their norms at each
    /// point of the current elementary cell
    void setNormals();


    T compute_matrix(int base_i, int base_j, int ind_pt,
                     int d1, int d2);
    void compute_RHS(int base_i, int ind_pt,
                     gsMatrix<T>& res);
    T compute_energy(int ind_pt);


    void decrease_coeff_smoothing(T ratio)
    {
        m_coeff_NL *= ratio;
        m_coeff_winslow *= ratio;
    }

private:
    inline void add_cross_prod(T a, int pos,
                               const gsVector<T>& vect,
                               gsMatrix<T> &res)
    {
        unsigned j = 0;
        for(int i = 0;i < 3;i++)
        {
            for(int k = 0;k < 3;k++)
            {
                if(k != i && k != pos)
                {
                    j = k;
                    break;
                }
            }
            if(i != pos)
                res(i, 0) += a*vect(j,0) * BaseNL::sgn(pos, j, i);
        }
    }

    inline T dot_prod_normal_mixedCross(int base, int dim,
                                        int ind_pt)
    {
        T res = 0.;
        for(unsigned dim2 = 0;dim2 < 3;dim2++)
        {
            res += m_normal(dim2, ind_pt)
                * m_mixed_pv(m_dim*3*base + dim*3
                             + dim2, ind_pt);
        }
        return res;
    }

    inline T dot_prod_normal_crossBasis(int base_i, int dim_i,
                                        int base_j, int dim_j,
                                        int ind_pt)
    {
        T val = 0.;
        int k = 0;
        if(dim_i == dim_j)
            return 0.;

        for(int dim = 0;dim < 3;dim++)
        {
            if(dim != dim_i && dim != dim_j)
            {
                k = dim;
                break;
            }
        }
        val = (*m_der1)(base_i * m_stride[0], ind_pt)
            * (*m_der1)(base_j * m_stride[0] + 1, ind_pt)
            - (*m_der1)(base_j * m_stride[0], ind_pt)
            * (*m_der1)(base_i * m_stride[0] + 1, ind_pt);
        return m_normal(k, ind_pt)
            * BaseNL::sgn(dim_i, dim_j, k) * val;
    }

    inline T dot_prod_mixedCross(int base_i, int dim_i,
                                 int base_j, int dim_j,
                                 int ind_pt)
    {
        T res = 0.;
        for(unsigned dim2 = 0;dim2 < 3;dim2++)
        {
            res += m_mixed_pv(m_dim*3*base_i + dim_i*3
                             + dim2, ind_pt)
                * m_mixed_pv(m_dim*3*base_j + dim_j*3
                             + dim2, ind_pt);
        }
        return res;
    }

protected:
    using Base::m_stride;
    using Base::m_der1;
    using Base::m_der1_curr;
    using Base::m_nActives;
    using Base::m_points;
    using Base::m_dim;
    using Base::m_patch;
    using BaseNL::m_coeff_NL;

    gsMatrix<T> m_mixed_pv;
    gsMatrix<T> m_mixed_grad;

    /// 1 over the norm of the normal
    /// (Dx cross Dy) at each point
    gsMatrix<T> m_inverse_area;

    /// the "normal" at each point
    /// (Dx cross Dy)
    gsMatrix<T> m_normal;

    /// the norm of the gradient at each point
    gsMatrix<T> m_squareNorm;

    /// the gradient at each point
    std::vector<gsVector3d<T> > m_Dx;
    std::vector<gsVector3d<T> > m_Dy;

    /// The proportion of Winslow energy in the
    /// global nonlinear energy
    T m_coeff_winslow;

}; /// gsFittingIntegrandW



/**
   Class representing the integrand of the Winslow energy
   in 3D.
*/
template<class T=real_t>
class gsFittingIntegrandW3D : public gsFittingIntegrandNL<3, T>
{

    typedef gsFittingIntegrand<3, T> Base;
    typedef gsFittingIntegrandNL<3, T> BaseNL;

public:
    /// Constructor
    gsFittingIntegrandW3D(T coeff_NL, T coeff_winslow3D) :
    BaseNL::gsFittingIntegrandNL(coeff_NL)
    {  m_coeff_winslow3D = coeff_winslow3D * m_coeff_NL;  }

    /// Destructor
    virtual ~gsFittingIntegrandW3D(){}

    gsFittingIntegrand<3, T>* copy()
    { return new gsFittingIntegrandW3D
            (m_coeff_NL, m_coeff_winslow3D/m_coeff_NL);  }

    void actualizeCurr();
    void actualize_energy()
    {
        computeMetric();
        computeValues();
    }
    void actualize(){  }

    T compute_matrix(int base_i, int base_j, int ind_pt,
                     int d1, int d2);
    void compute_RHS(int base_i, int ind_pt,
                     gsMatrix<T>& res);
    T compute_energy(int ind_pt)
    {
        return ( m_e1(ind_pt, 0) - m_e2(ind_pt, 0) )
            * m_inv_deter(ind_pt, 0);
    }


    void decrease_coeff_smoothing(T ratio)
    {
        m_coeff_NL *= ratio;
        m_coeff_winslow3D *= ratio;
    }

    void computeValues();

    inline void setIndDiagMetric(int der, int* ind)
    {
        if(der == 0){
            ind[0] = 3;
            ind[1] = 5;
        }
        else if(der == 1){
            ind[0] = 0;
            ind[1] = 5;
        }
        else{
            ind[0] = 0;
            ind[1] = 3;
        }
    }

    inline int complementary(int ind1, int ind2)
    {
        if(ind1 != 0 && ind2 != 0)
            return 0;
        else if(ind1 != 1 && ind2 != 1)
            return 1;
        else
            return 2;
    }

    inline T valMetric(int ind1, int ind2, int pt)
    {
        if(ind1 == 0)
            return m_metric(ind2, pt);
        else if(ind1 == 1)
            return m_metric(2 + ind2, pt);
        else
            return m_metric(5, pt);
    }

    void computeFirstDer();
    void computeSecondDer();
    void computeDerDeter();
    void computeMetric();

protected:
    using Base::m_stride;
    using Base::m_der1;
    using Base::m_der1_curr;
    using Base::m_nActives;
    using Base::m_points;
    using Base::m_dim;
    using Base::m_patch;
    using BaseNL::m_coeff_NL;


    /// the metric at each point
    gsMatrix<T> m_metric; /// SIZE: (6, nPoints)

    /// The determinant of the gradient at each point
    gsMatrix<T> m_deter; /// SIZE: (nPoints, 1)
    gsMatrix<T> m_inv_deter; /// SIZE: (nPoints, 1)

    /// The derivative of the gradient at each point
    gsMatrix<T> m_Ddeter; /// SIZE: (3 * nBasis, nPoints)

    /// The second derivative of the gradient at each point
    gsMatrix<T> m_DDdeter; /// SIZE: (3 * nBasis * nBasis, nPoints)

    /// The energy takes the form $\int (e_1(u) - e_2(u)) / det|\nabla(u)|$

    /// Contains the value of $e_1$ at each point
    gsMatrix<T> m_e1; /// SIZE: (nPoints, 1)
    /// Contains the value of $e_2$ at each point
    gsMatrix<T> m_e2; /// SIZE: (nPoints, 1)

    /// Contains the derivative of $e_1$ at each point
    /// for each basis element
    gsMatrix<T> m_De1; /// SIZE: (3 * nBasis, nPoints)
    /// Contains the derivative of $e_2$ at each point
    /// for each basis element
    gsMatrix<T> m_De2; /// SIZE: (3 * nBasis, nPoints)

    /// Contains the second derivative of $e_1$ at each point
    /// for each couple of basis element
    /// (primal and dual basis elements)
    gsMatrix<T> m_DDe1; /// SIZE: (9 * nBasis * nBasis, nPoints)
    gsMatrix<T> m_DDe2; /// SIZE: (9 * nBasis * nBasis, nPoints)

    /// The proportion of the 3D Winslow energy in the
    /// global nonlinear energy
    T m_coeff_winslow3D;

}; /// gsFittingIntegrandW3D


}// namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFittingIntegrand.hpp)
#endif

//////////////////////////////////////////////////
//////////////////////////////////////////////////
