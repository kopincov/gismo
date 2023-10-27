/** @file gsScatteredHFittingRefinement.h

    @brief Scattered fitting using hierarchical splines -- refinement.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Bracco, C. Giannelli, J. Speh, D. Mokris
*/

#pragma once

#include <vector>
#include <set>

namespace gismo {

gsMatrix<> compute_errors(const gsTHBSpline<2>& QI, 
                          const gsMatrix<>& parameters, 
                          const gsMatrix<>& points)
{
    gsMatrix<> eval;
    QI.eval_into(parameters, eval);
    
    gsMatrix<> errors(1, eval.cols());
    
    for (int col = 0; col != eval.cols(); col++)
    {
        errors(0, col) = (eval.col(col) - points.col(col)).norm();
    }
    return errors;
}


real_t percentage_of_points_below_threshold(const gsMatrix<>& errors,
                                            const real_t threshold)
{
    int num_pts_below_threshold = 0;
    for (int i = 0; i != errors.cols(); i++)
    {
        if (errors(0, i) < threshold)
        {
            num_pts_below_threshold++;
        }
    }
    
    return 1.0 * num_pts_below_threshold / errors.cols();
}


bool is_point_inside_support(const gsMatrix<>& parameter, 
                             const gsMatrix<>& support)
{
    const real_t x = parameter(0, 0);
    const real_t y = parameter(1, 0);
    
    return support(0, 0) <= x && x < support(0, 1) &&
        support(1, 0) <= y && y < support(1, 1);
}

bool is_point_inside_support(const real_t x,
                             const real_t y,
                             const gsMatrix<>& support)
{
    return support(0, 0) <= x && x < support(0, 1) &&
        support(1, 0) <= y && y < support(1, 1);
}

index_t num_pts_inside_support(const gsMatrix<>& parameters, const gsMatrix<>& support)
{
    int num_pts_inside_support = 0;
    
    for (int col = 0; col != parameters.cols(); col++)
    {
        if (is_point_inside_support(parameters(0,col), parameters(1,col), support))
        {
            num_pts_inside_support++;
        }
    }
    return num_pts_inside_support;    
}


bool there_are_at_least_n_pts_in_support(const gsMatrix<>& parameters, const gsMatrix<>& support, int n)
{
    int num_pts_inside_support = 0;

    for (int col = 0; col != parameters.cols(); col++)
    {
        //gsMatrix<> parameter = parameters.col(col);
        //if (is_point_inside_support(parameter, support))
        if (is_point_inside_support(parameters(0,col), parameters(1,col), support))
        {
            num_pts_inside_support++;
            if(num_pts_inside_support>=n)
                return true;
        }
    }
    return false;
}


// The version in the paper.
bool support_contains_at_least_n_loc_points(const gsTHBSplineBasis<2>& basis, 
					    const int basis_fun_index, 
					    const index_t n_loc, 
					    const gsMatrix<>& parameters,
					    const gsMatrix<>& errors,
					    const real_t threshold)
{
    const gsMatrix<> support = basis.support(basis_fun_index);
    return (n_loc <= num_pts_inside_support(parameters,support));
}

// The first improvement.
bool each_support_quadrant_contains_at_least_n_loc_points(const gsTHBSplineBasis<2>& basis, 
							  const int basis_fun_index, 
							  const index_t n_loc, 
							  const gsMatrix<>& parameters,
							  const gsMatrix<>& errors,
							  const real_t threshold)
{
    const gsMatrix<> support = basis.support(basis_fun_index);
    
    const real_t x_min = support(0,0);
    const real_t x_max = support(0,1);
    const real_t x_mid = (x_min + x_max)/2;

    const real_t y_min = support(1,0);
    const real_t y_max = support(1,1);
    const real_t y_mid = (y_min + y_max)/2;

    gsMatrix<> top_left_supp(2,2);
    top_left_supp << x_min, x_mid, y_mid, y_max;
    const index_t top_left = num_pts_inside_support(parameters,top_left_supp);

    gsMatrix<> top_right_supp(2,2);
    top_right_supp << x_mid, x_max, y_mid, y_max;
    const index_t top_right = num_pts_inside_support(parameters,top_right_supp);

    gsMatrix<> bot_left_supp(2,2);
    bot_left_supp << x_min, x_mid, y_min, y_mid;
    const index_t bot_left = num_pts_inside_support(parameters,bot_left_supp);

    gsMatrix<> bot_right_supp(2,2);
    bot_right_supp << x_mid, x_max, y_min, y_mid;
    const index_t bot_right = num_pts_inside_support(parameters,bot_right_supp);

    const index_t controlSum = top_left+top_right+bot_left+bot_right;
    if(controlSum !=num_pts_inside_support(parameters,support))
      std::cerr << "Bug, " << controlSum << "!=" << num_pts_inside_support(parameters,support) << "!" << std::endl;

    return(
	   top_left >= (n_loc / 4) &&
	   top_right >= (n_loc / 4) &&
	   bot_left >= (n_loc / 4) &&
	   bot_right >= (n_loc / 4));
}

// The second improvement.
bool each_child_support_contains_at_least_n_loc_points(const gsTHBSplineBasis<2>& basis, 
						       const int basis_fun_index, 
						       const index_t n_loc, 
						       const gsMatrix<>& parameters,
						       const gsMatrix<>& errors,
						       const real_t threshold)
{
    if(basis.degree(0) != 3 || basis.degree(1) != 3)
    {
	std::cerr << "This strategy is implemented only for degree 3 so far." << std::endl;
	return false;
    }

    // We gather the knots of the next (finer) level inside the support.
    // TODO: Is there a more elegant method?
    const gsMatrix<> support = basis.support(basis_fun_index);
    // std::cout << "Support:" << std::endl << support << std::endl;
    // std::cin.get();
    
    const real_t u_min = support(0,0);
    const real_t u_max = support(0,1);
    const real_t v_min = support(1,0);
    const real_t v_max = support(1,1);
 
    std::vector<real_t> u_knots, v_knots;
    size_t deg = basis.degree(0);
    size_t ord = 2 * (deg + 1);
    for(size_t i = 0; i <= ord; ++i)
    {
        u_knots.push_back(((ord-i) * u_min + i * u_max)/ord);
        v_knots.push_back(((ord-i) * v_min + i * v_max)/ord);
    }

    // Now we make sure that support of each finer function contains at least n_loc points.
    for(size_t i = 0; i <= deg + 1; ++i)
    {
        for(size_t j = 0; j <= deg + 1; ++j)
        {
            gsMatrix<> smallSupp(2,2);
            smallSupp << u_knots[i], u_knots[i+deg+1], v_knots[j], v_knots[j+deg+1];
            if( num_pts_inside_support(parameters, smallSupp) <= n_loc )
                return false;
        }
    }
    return true;
}

// Dispatch function based on the strategy
bool should_be_refined(const gsTHBSplineBasis<2>& basis, 
		       const int basis_fun_index, 
		       const index_t n_loc, 
		       const gsMatrix<>& parameters,
		       const gsMatrix<>& errors,
		       const real_t threshold,
		       int strategy)
{
    switch(strategy)
    {
    case 0:
      return support_contains_at_least_n_loc_points(basis,basis_fun_index,n_loc,parameters, errors, threshold);
    case 1:
      return each_support_quadrant_contains_at_least_n_loc_points(basis,basis_fun_index,n_loc,parameters, errors, threshold);
    case 2:
      return each_child_support_contains_at_least_n_loc_points(basis,basis_fun_index,n_loc,parameters, errors, threshold);
    default:
      std::cerr << "Unknown option to should_be_refined." << std::endl;
      return false;
    }
}


std::vector<index_t> get_marked_functions(gsTHBSplineBasis<2>& basis, 
                                           const gsMatrix<>& parameters, 
                                           const gsMatrix<>& errors, 
                                           const real_t threshold,
                                           const index_t n_loc,
					   int strategy)
{
    std::set<index_t> pre_marked_functions;
    
    for (int point_index = 0; point_index != parameters.cols(); point_index++)
    {
        const real_t error = errors(0, point_index);
        if (threshold < error)
        {
            gsMatrix<> parameter = parameters.col(point_index);
            gsMatrix<index_t> active_functions;
            basis.active_into(parameter, active_functions);
            for (int row = 0; row != active_functions.rows(); row++)
            {
                pre_marked_functions.insert(active_functions(row, 0));
            }
        }
    }

    std::set<index_t>::iterator it;
    std::vector<index_t> marked_functions;
    for (it = pre_marked_functions.begin(); it != pre_marked_functions.end(); ++it)
    {
        if (should_be_refined(basis, *it, n_loc, parameters, errors, threshold, strategy))
            marked_functions.push_back(*it);
    }

    return marked_functions;
}

void refine_marked_functions(gsTHBSplineBasis<2>& basis, 
                             const std::vector<index_t>& marked_functions)
{
    std::vector<index_t> boxes;
    for (size_t fun = 0; fun != marked_functions.size(); fun++)
    {
        const unsigned flat_index = basis.flatTensorIndexOf(marked_functions[fun]);
        const int level = basis.levelOf(marked_functions[fun]);
        const gsTensorBSplineBasis<2>& ten_basis = *basis.getBases()[level];
        //std::cout << "global index" << marked_functions[fun] << "level" << level << std::endl;
        
        gsMatrix<index_t, 2, 2> support_indices;
        ten_basis.elementSupport_into(flat_index, support_indices);
        
        boxes.push_back(level + 1);
        boxes.push_back(support_indices(0, 0) * 2);
        boxes.push_back(support_indices(1, 0) * 2);
        boxes.push_back(support_indices(0, 1) * 2);
        boxes.push_back(support_indices(1, 1) * 2);        
    }
    
    basis.refineElements(boxes);
}

    // void print(unsigned x)
    // {
    //     std::cout << " " << x;
    // }

void scatteredHFittingRefine(gsTHBSplineBasis<2>& basis, 
                             const gsMatrix<>& parameters, 
                             const gsMatrix<>& errors, 
                             const real_t threshold,
                             const index_t n_loc,
			                 int strategy)
{
    std::cout << "Strategy: " << strategy << std::endl;

    // std::cout << "errors: \n" << std::setprecision(16) 
    //           << errors.transpose() << std::endl;
    
    // basis.uniformRefine();
    std::vector<index_t> marked_functions = get_marked_functions(basis, parameters, 
                                                                  errors, threshold, 
                                                                  n_loc, strategy);
    
    refine_marked_functions(basis, marked_functions);
    
    // std::cout << "marked functions: ";
    // std::for_each(marked_functions.begin(), marked_functions.end(), print);
    // std::cout << std::endl;
}

}

