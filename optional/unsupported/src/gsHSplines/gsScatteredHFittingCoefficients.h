/** @file gsScatteredHFittingCoefficients.h

    @brief Scattered fitting using hierarchical splines -- computation of coefficients.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Bracco, C. Giannelli, J. Speh, D. Mokris
*/

#pragma once

#include <algorithm>
#include <cmath>

#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSpline.h>

namespace gismo {

class TransformationToUnitSquare
{
public:
    TransformationToUnitSquare(const gsMatrix<>& square)
        : dx(square(0, 1) - square(0, 0)),
          dy(square(1, 1) - square(1, 0)),
          origin_x(square(0, 0)),
          origin_y(square(1, 0))
    {
        
    }

    void transform(const real_t x,
                   const real_t y,
                   real_t& out_x,
                   real_t& out_y)
    {
        out_x = (x - origin_x) / dx;
        out_y = (y - origin_y) / dy;
    }
    
private:
    real_t dx;
    real_t dy;
    real_t origin_x;
    real_t origin_y;
};


gsMatrix<> get_center(const gsMatrix<>& support)
{
    gsMatrix<> c(2, 1);
    
    c(0, 0) = (support(0, 0) + support(0, 1)) / 2;
    c(1, 0) = (support(1, 0) + support(1, 1)) / 2;
    
    return c;
}


std::vector<index_t> get_points_inside_rectangle(const gsMatrix<>& parameters,
					     const gsMatrix<>& center,
					     const real_t u_diff,
					     const real_t v_diff)
{
    std::vector<index_t> indices;
    
    for (int i = 0; i != parameters.cols(); i++)
    {
        if( (std::fabs(parameters(0,i) - center(0,0)) <= u_diff) &&
            (std::fabs(parameters(1,i) - center(1,0)) <= v_diff))
        {
            indices.push_back(i);
        }
    }
    
    return indices;
}
					     

std::vector<index_t> get_points_inside_circle(const gsMatrix<>& parameters, 
                                          const gsMatrix<>& center,
					  const real_t radius)
{
    std::vector<index_t> indices;
    
    for (int i = 0; i != parameters.cols(); i++)
    {
        const gsMatrix<> point = parameters.col(i);

	if ((point - center).squaredNorm() <= radius * radius)
        {
            indices.push_back(i);
	}
    }
    
    return indices;
}

real_t maximum_interval(const std::vector<real_t>& breaks)
{
    // assumption uniform knots, maximum interval is the 
    // same as second knot

  // TODO: What if breaks[0] is not zero?
    
    return breaks[1];
}

int maximum_enlarging_factor(const gsTHBSplineBasis<2>& basis, 
                             const index_t basis_fun_index)
{
    // TO DO, compute this for each level in advance
    
    // real_t diagonal_of_element = 0;
    
    int level = basis.levelOf(basis_fun_index);
    const gsTensorBSplineBasis<2>& ten_basis = *basis.getBases()[level];
    const int degree = ten_basis.maxDegree();
    const gsMatrix<> support = basis.support(basis_fun_index);
    const real_t radius = (support.col(0) - support.col(1)).norm() / 2;
    
    const gsKnotVector<>& kv1 = ten_basis.knots(0);
    const gsKnotVector<>& kv2 = ten_basis.knots(1);

    const std::vector<real_t> breaks1 = kv1.breaks();
    const std::vector<real_t> breaks2 = kv2.breaks();

    const real_t cell_side_x = maximum_interval(breaks1);
    const real_t cell_side_y = maximum_interval(breaks2);

    const real_t cell_side_x_square = cell_side_x * cell_side_x;
    const real_t cell_side_y_square = cell_side_y * cell_side_y;
    
    const real_t cell_diag = math::sqrt(cell_side_x_square + cell_side_y_square);
    
    // delta^{ell-1}
    // look at paper Cesare Bracco: adaptive scattered data fitting...
    const real_t delta_ell_1 = cell_diag * (degree + 1);
    
    const int K_ell = cast<real_t,index_t>(math::ceil(2 * delta_ell_1 / radius) + 1);
    
    return K_ell;
}

void get_data_for_least_square_approx(const gsTHBSplineBasis<2>& basis,
                                      const gsMatrix<>& parameters, 
                                      const int basis_fun_index,
                                      int& degreeLS, 
                                      std::vector<index_t>& indicesLS,
				      bool circle)
{
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
        
        const int num_LS_points = static_cast<index_t>(indicesLS.size());
        
        if (num_LS_points == 0)
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
            break;
        }
    }
}


void assemble_LS_system(const gsTHBSplineBasis<2>& basis,
                        const int basis_fun_index,
                        const std::vector<index_t>& indices, 
                        const gsMatrix<>& parameters, 
                        const gsMatrix<>& points, 
                        const int degree, 
                        gsMatrix<>& A, 
                        gsMatrix<>& B)
{
    const int dimension = (degree + 1) * (degree + 2) / 2;
    A.resize(dimension, indices.size());
    B.resize(indices.size(), points.rows());
    
    const gsMatrix<> support = basis.support(basis_fun_index);
    TransformationToUnitSquare transformation(support);
    
    // const real_t dx = support(0, 1) - support(0, 0);
    // const real_t dy = support(1, 1) - support(1, 0);
    // const real_t origin_x = support(0, 0);
    // const real_t origin_y = support(1, 0);
    
    for (size_t i = 0; i != indices.size(); i++)
    {
        const int index = indices[i];
        real_t xx, yy;
        const real_t x = parameters(0, index);
        const real_t y = parameters(1, index);
        transformation.transform(x, y, xx, yy);
        // const real_t xx = (x - origin_x) / dx;
        // const real_t yy = (y - origin_y) / dy;
        
        int j = 0;
        for (int exp_x = 0; exp_x <= degree; exp_x++)
        {
            for (int exp_y = 0; exp_y <= degree - exp_x; exp_y++)
            {
                A(j, i) = math::pow(xx, exp_x) * math::pow(yy, exp_y);
                j++;
            }
        }
        B.row(i) = points.col(index).transpose();
    }
}

int compute_rank(const gsMatrix<>& A)
{
    // Eigen::JacobiSVD< Eigen::Matrix<real_t, Dynamic, Dynamic> > svd(A);
    Eigen::FullPivLU< Eigen::Matrix<real_t, Dynamic, Dynamic> > lu_decomp(A);
    return lu_decomp.rank();
}

void assemble_regular_LS_system(const gsTHBSplineBasis<2>& basis,
                                const int basis_fun_index,
                                const std::vector<index_t>& indices, 
                                const gsMatrix<>& parameters, 
                                const gsMatrix<>& points, 
                                int& deg, 
                                gsMatrix<>& A, 
                                gsMatrix<>& B)
{
    const int k = basis_fun_index;
    
    while (true)
    {
        assemble_LS_system(basis, k, indices, parameters, points, deg, A, B);
            
        if (A.rows() == compute_rank(A))
        {
            break;
        }
        else
        {
            deg--;
        }
    }
}


gsMatrix<> generate_points(const gsTHBSplineBasis<2>& basis, 
                           const int basis_fun_index)
{
   const int level = basis.levelOf(basis_fun_index);
   const gsTensorBSplineBasis<2>& ten_basis = *basis.getBases()[level];
   const int degree = ten_basis.maxDegree();

   const gsMatrix<> support = basis.support(basis_fun_index);   
   
   const real_t factor = 1.0 / (2 * (degree + 1));
   const gsMatrix<> diff = factor * (support.col(1) - support.col(0));
   
   const int n = degree + 1;
   
   gsVector<unsigned> num_points(2);
   num_points << n, n;
   
   gsVector<> start = support.col(0);
   start(0) += diff(0, 0);
   start(1) += diff(1, 0);

   gsVector<> end = support.col(1);
   end(0) -= diff(0, 0);
   end(1) -= diff(1, 0);
   
   gsMatrix<> points = gsPointGrid(start, end, num_points);
   
   // std::cout << "support: \n" << support << "\n"
   //           << "points:\n" << points << "\n"
   //           << "size: " << points.rows() << " x " << points.cols() << std::endl;

   return points;
}


std::set<index_t> matrix_to_set(const gsMatrix<index_t>& matrix)
{
    std::set<index_t> set;
    for (int row = 0; row != matrix.rows(); row++)
    {
        for (int col = 0; col != matrix.cols(); col++)
        {
            set.insert(matrix(row, col));
        }
    }
    
    return set;
}


std::map<index_t, index_t> compute_ordering(const std::set<index_t>& indices)
{
    int row_index = 0;
    std::set<index_t>::iterator it;
    std::map<index_t, index_t> ordering;
    for (it = indices.begin(); it != indices.end(); ++it)
    {
        ordering[*it] = row_index;
        row_index++;
    }
    
    return ordering;
}

int number_different_basis_function(const gsTensorBSplineBasis<2>& ten_basis,
                                    const std::set<index_t>& set_indices, 
                                    const int direction)
{
    std::set<index_t>::iterator it;
    std::set<index_t> indices_in_direction;
    
    for (it = set_indices.begin(); it != set_indices.end(); ++it)
    {
        gsVector<index_t> tensor_index = ten_basis.tensorIndex(*it);
        indices_in_direction.insert(tensor_index(direction));
    }
    
    return static_cast<index_t>(indices_in_direction.size());
}

gsMatrix<> computeGridPoints(const gsTensorBSplineBasis<2>& ten_basis, 
                             const gsMatrix<>& support, 
                             const std::set<index_t>&set_indices)
{
   
    
   int num_basis_fun_x = number_different_basis_function(ten_basis, set_indices, 0);
   int num_basis_fun_y = number_different_basis_function(ten_basis, set_indices, 1);
   
   gsVector<unsigned> num_points(2);
   num_points << num_basis_fun_x, num_basis_fun_y;
   
   gsVector<> start = support.col(0);
   gsVector<> end = support.col(1);
   
   return gsPointGrid(start, end, num_points);
}


gsMatrix<> assemble_system_coefsLS_to_basis_fun(const gsTHBSplineBasis<2>& basis, 
                                                const int basis_fun_index, 
                                                gsMatrix<>& gridPoints,
                                                int& basis_fun_tensor_index)
{
   const int level = basis.levelOf(basis_fun_index);
   const gsTensorBSplineBasis<2>& ten_basis = *basis.getBases()[level];
   const unsigned tensor_index = basis.flatTensorIndexOf(basis_fun_index);
   
   gsMatrix<> tmpPoints = generate_points(basis, basis_fun_index);

   gsMatrix<index_t> indices;
   ten_basis.active_into(tmpPoints, indices);
   
   std::set<index_t> set_indices = matrix_to_set(indices);
   std::map<index_t, index_t> basis_to_cols = compute_ordering(set_indices);

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
           }
       }
   }

   // std::for_each(set_indices.begin(), set_indices.end(), print);
   // std::cout << std::endl;
   // std::cout << "tensor_index:" << tensor_index << "\n"
   //           << "col: " << basis_to_cols[tensor_index] << std::endl;
   
   basis_fun_tensor_index = basis_to_cols[tensor_index];
   return A;
}


gsMatrix<> polynomial_basis_evaluation(const int degreeLS, 
                                       const gsMatrix<>& gridPoints, 
                                       const gsMatrix<>& support)
{
    const int dimension = (degreeLS + 1) * (degreeLS + 2) / 2;
    gsMatrix<> P(gridPoints.cols(), dimension);
    TransformationToUnitSquare transformation(support);
    
    for (int point_index = 0; point_index != gridPoints.cols(); point_index++)
    {
        int j = 0;
        real_t xx, yy;
        const real_t x = gridPoints(0, point_index);
        const real_t y = gridPoints(1, point_index);
        transformation.transform(x, y, xx, yy);
        
        for (int exp_x = 0; exp_x <= degreeLS; exp_x++)
        {
            for (int exp_y = 0; exp_y <= degreeLS - exp_x; exp_y++)
            {
                P(point_index, j) = math::pow(xx, exp_x) * math::pow(yy, exp_y);
                j++;
            }
        }
    }
    
    return P;
}


gsMatrix<> get_coefficient(const gsTHBSplineBasis<2>& basis, 
                           const gsMatrix<>& parameters, 
                           const gsMatrix<>& points, 
                           const int k, 
                           int degreeLS, 
                           const std::vector<index_t>& indicesLS,
                           const real_t sigma)
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
        int row_index = -1;
  
    
        gsMatrix<> S = assemble_system_coefsLS_to_basis_fun(basis, k, gridPoints, 
                                                            row_index);
        const gsMatrix<> support = basis.support(k);    
        gsMatrix<> P = polynomial_basis_evaluation(degreeLS, gridPoints, support);
        
        gsMatrix<> rhs = P * coefficientsLS;
        
        gsMatrix<> coefficients = S.fullPivHouseholderQr().solve(rhs);
        gsMatrix<> coefficient = coefficients.row(row_index);
        
        // gsInfo << std::setprecision(16) << coefficient << "\n";
        
        // ------------------------------------------------------------
        
        //std::cout << "size S: " << S.rows() << " x " << S.cols() << "\n"
        //          << "size P: " << P.rows() << " x " << P.cols() << "\n"
        //          << "size A: " << A.rows() << " x " << A.cols() << "\n";
        gsMatrix<> AAT_1 = (A * A.transpose()).inverse();
        gsMatrix<> G = S.inverse() * P * AAT_1 * A;
        
        // std::cout << "size G: " << G.rows() << " x " << G.cols() << "\n";
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
            //std::cout << "degree: " << degreeLS << std::endl;
            return coefficient;
        }
    }
}


// 2018-06-12: New parameter circle; if true, enlargment has the form
// of a circle and of a rectangle otherwise.
gsMatrix<> get_coefficients(const gsTHBSplineBasis<2>& basis, 
                            const gsMatrix<>& parameters, 
                            const gsMatrix<>& points, 
                            const real_t sigma,
			                bool circle = true)
{
    gsMatrix<> coefficients(basis.size(), points.rows());

    for (index_t k = 0; k != basis.size(); k++)
    {
        // LS - least squares
        int degreeLS = 0;
        std::vector<index_t> indicesLS;

        get_data_for_least_square_approx(basis, parameters, k, degreeLS, indicesLS, circle);

        gsMatrix<> coefficient = get_coefficient(basis, parameters, points, k, 
                                                 degreeLS, indicesLS, sigma);

        coefficients.row(k) = coefficient.row(0);
    }


    return coefficients;
}

}
