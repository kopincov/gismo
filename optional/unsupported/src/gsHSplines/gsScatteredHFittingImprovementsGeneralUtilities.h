/** @file gsScatteredHFittingImprovementsGeneralUtilities.h

    @brief Improvements of scattered fitting.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#pragma once

//#include <gsHSplines/gsScatteredHFittingRefinement.h>
//#include <gsHSplines/gsScatteredHFittingCoefficients.h>


namespace gismo
{
    

/** \brief biggest value of an std::vector
 */
template <class T>
T maximum_value(const std::vector<T>& vector)
{
  if(vector.size() == 0)
    return 0;

  T maximum = vector[0];
  
    for(int j = 1; j < vector.size(); j++)
    {
      if(maximum < vector[j])
      maximum = vector[j];
    }
    
  return maximum;
}


/** \brief biggest value of an gsVector
 */
template <class T>
T gs_maximum_value(const gsVector<T>& vector)
{
  if(vector.size() == 0)
    return 0;

  T maximum = vector(0);
  
    for(int j = 1; j < vector.size(); j++)
    {
      if(maximum < vector(j))
      maximum = vector(j);
    }
    
  return maximum;
}


/** \brief smallest value of an std::vector
 */
template <class T>
T minimum_value(const std::vector<T>& vector)
{
    if(vector.size() == 0)
        return 0;
    
    T minimum = vector[0];
    
    for(int j = 1; j < vector.size(); j++)
    {
        if( minimum > vector[j])
            minimum = vector[j];
    }
    
    return minimum;
}


/** \brief smallest value of a gsVector
 */
template <class T>
T gs_minimum_value(const gsVector<T>& vector)
{
  if(vector.size() == 0)
    return 0;

  T minimum = vector(0);
  
    for(int j = 1; j < vector.size(); j++)
    {
      if(minimum > vector(j))
      minimum = vector(j);
    }
    
  return minimum;
}


/** \brief biggest value of an std::set
 */
template <class T>
T set_maximum_value(const std::set<T>& set)
{
    T maximum;
    if(!set.empty())
        maximum = *set.rbegin();
    else
    {
        std::cout << "empty set." << std::endl;
        maximum = 0;
    }
    
    return maximum;
}


/** \brief smallest value of an std::set
 */
template <class T>
T set_minimum_value(const std::set<T>& set)
{
    T minimum;
    if(!set.empty())
        minimum = *set.begin();
    else
    {
        std::cout << "empty set." << std::endl;
        minimum = 0;
    }
    
    return minimum;
}


/** \brief average value of an std::vector
 */
template <class T>
T arithmetic_average_stdVector(const std::vector<T>& vector)
{
    T average;
    T counting = vector[0];
    for( int i = 1; i < vector.size(); i++)
        counting = counting + vector[i];
    
    average = counting/vector.size();
    return average;
}


/** \brief Insetion of elements at the beginning of an std::vector
 */
template <class T>
void push_front(std::vector<T>& vector, T element)
{
    std::vector<T> result(vector.size()+1);
    for(size_t i=0; i<vector.size(); ++i)
        result[i+1]=vector[i];
    result[0] = element;
    vector=result;
}



/** \brief elements of std::vector \em vector are stored in to std::set \em set
 */
template<class T>
std::set<T> vector_to_set(const std::vector<T>& vector)
{
    std::set<T> set;
    
    for(unsigned k = 0; k < vector.size(); k++)
        set.insert(vector[k]);
    
    return set;
}

/** \brief elements of std::set \em set are stored in std::vector \em vector
*/
template <class T>
std::vector<T> set_to_vector(const std::set<T>& set) // Values of a set stored in a std::vector
{
    typename std::set<T>::const_iterator it;
    std::vector<T> vector;
    
    for (it = set.begin(); it != set.end(); ++it)
      vector.push_back(*it);
    
    return vector;
}

/** \brief elemets of gsMatrix \em matrix are stored in std::set \em set
 */
template <class T>
std::set<T> matrix_to_set(const gsMatrix<T>& matrix)
{
    std::set<T> set;
    for(int i = 0; i < matrix.rows(); i++)
    {
        for(int j = 0; j < matrix.cols(); j++)
            set.insert(matrix(i,j));
    }
    
    return set;
            
}

/** \brief elemets of gsMatrix \em matrix are stored in std::vector \em vector
 */
template <class T>
std::vector<T> matrix_to_vector(const gsMatrix<T>& matrix)
{
    std::vector<T> vector;
    for(int i = 0; i < matrix.rows(); i++)
    {
        for(int j = 0; j < matrix.cols(); j++)
            vector.push_back(matrix(i,j));
    }
    
    return vector;
}


// called from get_coefficients_print_LS_surface
void data_least_square_approx(const gsTHBSplineBasis<2>& basis,
                              const gsMatrix<>& parameters, 
                              const int basis_fun_index,
                              int& degreeLS,
                              std::vector<index_t>& indicesLS,
                              const int cp_magnitude,
                              int& enlarging_nr,
                              bool circle
                             )
{
    enlarging_nr = 0;
    gsMatrix<> support = basis.support(basis_fun_index);
    gsMatrix<> center = get_center(support);
    real_t radius = (support.col(0) - support.col(1)).norm() / 2;

    real_t u_diff = std::abs(support(0,0) - support(0,1)) / 2;
    real_t v_diff = std::abs(support(1,0) - support(1,1)) / 2;

    const int degree = basis.maxDegree();
    
    const int K = maximum_enlarging_factor(basis, basis_fun_index);
    
    for (int k = 1; k <= K; k++)
    {
        // std::cout << "  enlarging k = " << k << std::endl;

	if(circle)
	    indicesLS = get_points_inside_circle(parameters, center, radius * k);
	else
	    indicesLS = get_points_inside_rectangle(parameters, center, u_diff * k, v_diff * k);
        
        const int num_LS_points = static_cast<int>(indicesLS.size());
        
        if (num_LS_points <= cp_magnitude) // if ( num_LS_points == 0)
        {
            continue;
        }
        else
        {
            
	        int d = degree;
            while ( num_LS_points < (d + 1) * (d + 2) / 2)
            {
                d -= 1;
            }
            degreeLS = d;
            enlarging_nr = k;
            break;
        }
    }
}



gsTensorBSplineBasis<2, real_t> tensorbasis_constructor(real_t u_min, real_t u_max, real_t v_min, real_t v_max, int interior, int degree, int multEnd)
{
    gsKnotVector<> kx(u_min, u_max, interior, multEnd); // knots vector, x direction
    gsKnotVector<> ky(v_min, v_max, interior, multEnd); // knots vector y direction
    gsTensorBSplineBasis<2, real_t> tens_basis(kx, ky); // tensor-product basis

    return tens_basis;
}


/** elements of a std::vector \em vector are stored in a gsMatrix \em matrix 
 */
template<class T>
gsMatrix<T> copy_in_gs_matrix(std::vector<T>& vector)
{
    gsMatrix<T> matrix(vector.size(), 1);
    for(size_t j=0; j<vector.size(); j++)
    matrix(j,0) = vector[j];
    return matrix;
}


// avoiding errors due to rounding
bool almostEqual(real_t left, real_t right)
{
    //return (math::abs(left-right) < 1e-12);
    return (fabs(left-right) < 1e-12);
}

// avoiding errors due to rounding
template <class T>
bool almostEqual(const gsMatrix<T>& left, const gsMatrix<T>& right)
{
    for(size_t i=0; i<left.rows(); i++)
        for(size_t j=0; j<left.cols(); j++)
            if( !almostEqual(left(i,j), right(i,j)) )
                return false;
            
    return true;
}

// avoiding errors due to rounding
bool betweenEqual(real_t value, real_t left, real_t right)
{
    if( (left < value) && (value < right) )
        return true;
    else
        return(almostEqual(left, value) || almostEqual(value, right));
}



    
template <class T> // Existing function of Eigen.
T scalr_product(std::vector<T> vector1,
                std::vector<T> vector2
                )
{
    T sp = 0;
    for(unsigned i = 0; i < vector1.size(); i++)
        sp = sp + vector1[i] * vector2[i];
    
    return sp;
    
}

/** \brief returns true if the point of coordinates (\em x, \em y) is inside the support \em supp.
 */
bool is_inside_the_support_2D(const gsMatrix<>& supp,
                              const real_t x,
                              const real_t y
                             )
{
    return (
            ((supp(0,0) <= x)||almostEqual(supp(0,0), x))&&
            ((x <= supp(0,1))||almostEqual(x, supp(0,1)))&&
            ((supp(1,0) <= y)||almostEqual(supp(1,0), y))&&
            ((y <= supp(1,1))||almostEqual(supp(1,1), y))
           );
}


/** \brief indices of parametrized points \em parameters inside the support \em support
 */
std::vector<index_t> parameters_inside_support_idx(const gsMatrix<>& support,
                                               const gsMatrix<>& parameters
                                              )
{
    std::vector<index_t> indices;
    for(int p = 0; p < parameters.cols(); p++)
    {
        if(is_inside_the_support_2D(support, parameters(0,p), parameters(1,p)))
           indices.push_back(p);
    }
    
    return indices;
}


}
