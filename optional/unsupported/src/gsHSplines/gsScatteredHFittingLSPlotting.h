/** @file gsScatteredHFittingLSPlotting.h

    @brief Improvements of scattered fitting.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#pragma once

#include <gsHSplines/gsScatteredHFittingRefinement.h>
#include <gsHSplines/gsScatteredHFittingCoefficients.h>
#include <gsHSplines/gsScatteredHFittingImprovementsGeneralUtilities.h>
#include <gsHSplines/gsScatteredHFittingImprovements.h>
#include <gsIO/gsWriteParaview.h>

namespace gismo
{

/** \brief regular grid made by points belonging to the whole working domain, e.g. [0, 1]x[0, 1].
 */
void generate_points_whole_domain(const gsTHBSplineBasis<2>& basis, 
                                        const int basis_fun_index,
                                        gsMatrix<>& points,
                                        const int level
                                       ) // For fitting LS surface
{
   //const int level = basis.levelOf(basis_fun_index);
   const gsTensorBSplineBasis<2>& ten_basis = *basis.getBases()[level];
   const int degree = ten_basis.maxDegree();

   const gsMatrix<> support = basis.support(basis_fun_index);   
   
   const real_t factor = 1.0 / (2 * (degree + 1));
   const gsMatrix<> diff = factor * (support.col(1) - support.col(0));
   
   //const int n = 14; to fix in accordance to the considered function.
   //gsInfo << "Basis max degree: " << degree << std::endl;
   const int n = ten_basis.component(0).size();
   
   gsVector<unsigned> num_points(2);
   num_points << n, n;
   
   gsVector<> start(2);
   start << 0, 0;
   start(0) += diff(0, 0);
   start(1) += diff(1, 0);

   gsVector<> end(2);
   end << 1, 1;
   end(0) -= diff(0, 0);
   end(1) -= diff(1, 0);

   points = gsPointGrid(start, end, num_points);
}


/** \brief collocation matrix for the local fitting of \em basis_fun_index
 * The functions considered are the ones active on a regular point grid, computed either on the whole domain (whole_domain = true),
 * or on the suport of basis function given in input - the index of the basis_function is given in input.
 * the functions actually used for bulding the collocation matrix are store in a std::set.
 */
gsMatrix<> assemble_system_coefsLS_to_tensorbasis_fun(const gsTHBSplineBasis<2>& basis, 
                                                      const int basis_fun_index,
                                                      gsMatrix<>& gridPoints,
                                                      gsMatrix<index_t>& indices,
                                                      int& basis_fun_tensor_index,
                                                      std::set<index_t>& employed_functions,
                                                      const int level,
                                                      bool whole_domain
                                                     )
{
   //const int level = basis.levelOf(basis_fun_index);
   const gsTensorBSplineBasis<2>& ten_basis = *basis.getBases()[level];
   const unsigned tensor_index = basis.flatTensorIndexOf(basis_fun_index);
   
   gsMatrix<> tmpPoints;
   
   if(whole_domain)
       generate_points_whole_domain(basis, basis_fun_index, tmpPoints, level);
   else
       tmpPoints = generate_points(basis, basis_fun_index);
   
   ten_basis.active_into(tmpPoints, indices);
   
   std::set<index_t> set_indices = matrix_to_set(indices);
   std::map<index_t, index_t> basis_to_cols = compute_ordering(set_indices);

    gsMatrix<real_t> domain(2,2);
    domain << 0, 1, 0, 1;
    
    if(whole_domain)
        gridPoints = computeGridPoints(ten_basis, domain, set_indices);
    else
        gridPoints = computeGridPoints(ten_basis, basis.support(basis_fun_index), 
                                  set_indices);
  
   gsMatrix<> eval;
   ten_basis.eval_into(gridPoints, eval);
   ten_basis.active_into(gridPoints, indices);

   
   gsMatrix<> A(gridPoints.cols(), set_indices.size());
   A.setZero();
   
   
   for (int point_index = 0; point_index != eval.cols(); point_index++)
   {
       for (int index = 0; index != eval.rows(); index++)
       {
           const real_t basis_value = eval(index, point_index);
           const unsigned basis_index = indices(index, point_index);
           if (basis_to_cols.find(basis_index) != basis_to_cols.end())
           {
               const int col = basis_to_cols[basis_index];
               A(point_index, col) = basis_value;
               employed_functions.insert(basis_index);
           }
       }
   }


   
   basis_fun_tensor_index = basis_to_cols[tensor_index];
   return A;
}


/** \brief computation of coefficient for the function basis \em k .
 *         the data for the local fitting of \em k are printed in a paraview file.
 */
gsMatrix<> get_coefficient_pol(const gsTHBSplineBasis<2>& basis, 
                               const gsMatrix<>& parameters, 
                               const gsMatrix<>& points, 
                               const int k, 
                               int degreeLS, 
                               const std::vector<index_t>& indicesLS,
                               const real_t sigma,
                               int iteration,
                               const index_t cp_num,
                               const int it_max
                              )
{
    typedef Eigen::Matrix<real_t, Dynamic, Dynamic> matrix_t;
    
    
    
    while (true)
    {
        gsMatrix<> A;
        gsMatrix<> B; 

        assemble_regular_LS_system(basis, k, indicesLS, parameters, points, 
                                   degreeLS, A, B);
        

        Eigen::JacobiSVD<matrix_t> svdA(A.transpose(), 
                                        Eigen::ComputeThinU | Eigen::ComputeThinV);
        gsMatrix<> coefficientsLS = svdA.solve(B);
        
        // std::cout << k << " -> " << std::endl;

        // ------------------------------------------------------------
        gsMatrix<> gridPoints;
        gsMatrix<> gridPoints_ls;
        gsMatrix<index_t> indices;
        gsMatrix<index_t> indices_ls;
        std::set<index_t> employed_functions;
        std::set<index_t> employed_functions_ls;
        int row_index = -1;
        int row_index_ls = -1;
        
        const int level=basis.levelOf(k);
        
        gsMatrix<> S = assemble_system_coefsLS_to_tensorbasis_fun(basis, k, gridPoints, indices, row_index, employed_functions, level, false);
        gsMatrix<> W = assemble_system_coefsLS_to_tensorbasis_fun(basis, k, gridPoints_ls, indices_ls, row_index_ls, employed_functions_ls, 0, true);
    
        const gsMatrix<> support = basis.support(k);
        gsMatrix<> P = polynomial_basis_evaluation(degreeLS, gridPoints, support);
        
        gsMatrix<> support_w(2,2);
        support_w << 0, 1, 0, 1;
        gsMatrix<> P_ls = polynomial_basis_evaluation(degreeLS, gridPoints_ls, support);
        
        gsMatrix<> rhs = P * coefficientsLS;
        gsMatrix<> rhs_ls = P_ls * coefficientsLS;
        
    
        // S*coefficients = P*coefficientsLS
        // W*coefficients_ls = P_ls*coefficientsLS
        
        gsMatrix<> coefficients = S.fullPivHouseholderQr().solve(rhs);
        gsMatrix<> coefficients_ls = W.fullPivHouseholderQr().solve(rhs_ls);
        
        
       
        gsMatrix<> coefficient = coefficients.row(row_index);
        
        
        gsMatrix<> AAT_1 = (A * A.transpose()).inverse();
        gsMatrix<> G = S.inverse() * P * AAT_1 * A;
        
        
        Eigen::JacobiSVD<matrix_t> svdG(G);
        real_t largest_singular_value = svdG.singularValues()(0);

        if (sigma < largest_singular_value)
        {
            degreeLS--;
            if (degreeLS == -1)
            {
		GISMO_ERROR("Choose a larger value of sigma!");
                return coefficient;
            }
        }
        else
        {
            if( k == cp_num && iteration ==  it_max - 1 )
            {
                gsMatrix<> controlPoint = coefficient.transpose();
                gsWriteParaviewPoints(controlPoint, "CP_" + internal::to_string(k));
                
                gsMatrix<> restrPar(2, indicesLS.size());
                gsMatrix<> restrPoint(3, indicesLS.size());
                for(size_t i = 0; i < indicesLS.size(); i++)
                {
                    restrPar.col(i) = parameters.col(indicesLS[i]);
                    restrPoint.col(i) = points.col(indicesLS[i]);
                }
                gsWriteParaviewPoints(restrPar, "parameters_fun" + internal::to_string(k));
                gsWriteParaviewPoints(restrPoint, "points_fun" + internal::to_string(k));
                
                int level = basis.levelOf(k);
                gsTensorBSplineBasis<2> tensor_basis = *basis.getBases()[level];
                gsMatrix<> coef_pol(tensor_basis.size(), 3);
                coef_pol.setZero();
                std::vector<index_t> vector_employed_functions = set_to_vector(employed_functions);
                for( unsigned i = 0 ; i < vector_employed_functions.size(); i++)
                    coef_pol.row(vector_employed_functions[i]) = coefficients.row(i);
            
                gsTensorBSpline<2> tens_bspline(tensor_basis, coef_pol);
                gsMatrix<real_t> domain = basis.support(k);
                gsWriteParaview<real_t>(tens_bspline, domain, "iteration" + std::to_string(iteration) + "restr_to" + std::to_string(cp_num), 10000);
                
                gsTensorBSplineBasis<2> coarse_basis = *basis.getBases()[0];
                gsMatrix<> coef_pol_ls(coarse_basis.size(), 3);
                coef_pol_ls.setZero();
                std::vector<index_t> vector_employed_functions_ls = set_to_vector(employed_functions_ls);
                for( unsigned i = 0 ; i < vector_employed_functions_ls.size(); i++)
                    coef_pol_ls.row(vector_employed_functions_ls[i]) = coefficients_ls.row(i);
                
                gsTensorBSpline<2> tens_bspline_ls(coarse_basis, coef_pol_ls);
                gsWriteParaview(tens_bspline_ls, "iteration" + std::to_string(iteration) + "_LSsurfaceOf" + std::to_string(cp_num), 10000, false, true);
            }
            return coefficient;
        }

    } 
}




/** \brief  computation of coefficients for THBSplineBasis \em basis.
*/
gsMatrix<> get_coefficients_print_LS_surface(const gsTHBSplineBasis<2>& basis,
                                             const gsMatrix<>& parameters,
                                             const gsMatrix<>& points,
                                             const real_t sigma,
                                             int iteration,
                                             int cp_magnitude,
                                             const index_t cp_num,
                                             const int it_max,
                                             bool circle = true
                                            )
{

    gsMatrix<> coefficients(basis.size(), points.rows());

    for (index_t k = 0; k != basis.size(); k++)
    {
        // LS - least squares
        int degreeLS = 0;
        std::vector<index_t> indicesLS;
        int enlarging_nr = 0;
        
        data_least_square_approx(basis, parameters, k, degreeLS, indicesLS, cp_magnitude, enlarging_nr,
                                 circle);

        gsMatrix<> coefficient = get_coefficient_pol(basis, parameters, points, k, 
                                                 degreeLS, indicesLS, sigma, iteration, cp_num, it_max);

        coefficients.row(k) = coefficient.row(0);
    }
    
    return coefficients;
}




} // namespace gismo
