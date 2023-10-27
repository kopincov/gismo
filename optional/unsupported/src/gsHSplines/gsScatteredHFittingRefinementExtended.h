/** @file gsScatteredHFittingRefinementExtended.h

    @brief Improvements of Scattered fitting with hierarchical splines  -- refinement

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore, D. Mokris
*/

#pragma once

#include <gismo.h>
#include <gsSmoothPatches/gsCompositeHBasis.h>
//#include <gsHSplines/gsScatteredHFittingImprovementsStiffness.h>
//#include <gsHSplines/gsScatteredHFittingImprovementsGeneralUtilities.h>
#include <gsHSplines/gsScatteredHFittingRefinement.h>

namespace gismo{

/*gsMatrix<> enlarged_support_with_knot_span(const gsMatrix<>& support,
                                           const gsKnotVector<>& kv0,
                                           const gsKnotVector<>& kv1,
                                           const int k)
{
    gsMatrix<> new_support(2,2);
    
    gsKnotVector<real_t>::iterator it_u_left  = kv0.iFind(support(0,0));
    gsKnotVector<real_t>::iterator it_u_right = kv0.iFind(support(0,1));
    gsKnotVector<real_t>::iterator it_v_left  = kv1.iFind(support(1,0));
    gsKnotVector<real_t>::iterator it_v_right = kv1.iFind(support(1,1));

    real_t Ax = support(0,0); 
    real_t Bx = support(0,1);
    real_t Ay = support(1,0);
    real_t By = support(1,1);
    
    if(!almostEqual(Ax,0))
        Ax = *( it_u_left - k + 1 );
    if(!almostEqual(Bx,1))
        Bx = *( it_u_right + k -1 );
    if(!almostEqual(Ay,0))
        Ay = *( it_v_left - k + 1 );
    if(!almostEqual(By,1))
        By = *( it_v_right + k - 1);
    
    new_support << Ax, Bx, Ay, By;
    
    return new_support;
}
*/


/** \brief refinement strategy: if there are at least n_loc points inside the enlargment of \em basis_fun_index, the function is marked as to be refined.
 */
bool enlarged_support_contains_at_least_n_loc_points(const gsTHBSplineBasis<2>& basis,
                                                     const int basis_fun_index,
                                                     std::vector<index_t>& enlargement,
                                                     const index_t n_loc)
{
    if(basis_fun_index >= basis.size())
        gsInfo << "Wrong index." << std::endl;
    
    int num_points = enlargement[basis_fun_index]; 
    
    //return (n_loc <= num_points); 
    
    if(n_loc > num_points)
    {
        gsInfo << "Function: " << basis_fun_index << ". " ;
        gsInfo << "There are: " << num_points << "points inside the enlarged support." << std::endl;
        return false;
    }
    
    return true;
}



/** \brief refinement strategies
 */
bool should_be_refined_with_enlargment(const gsTHBSplineBasis<2>& basis,
                                       const int basis_fun_index,
                                       const int circle,
                                       std::vector<index_t>& enlargement,
                                       const index_t n_loc, 
                                       const gsMatrix<>& parameters,
                                       const gsMatrix<>& errors,
                                       const real_t threshold,
                                       int strategy)
{
    switch(strategy)
    {
    case 0:
      return support_contains_at_least_n_loc_points(basis, basis_fun_index, n_loc, parameters, errors, threshold);
    case 1:
      return each_support_quadrant_contains_at_least_n_loc_points(basis, basis_fun_index, n_loc, parameters, errors, threshold);
    case 2:
      return each_child_support_contains_at_least_n_loc_points(basis, basis_fun_index, n_loc, parameters, errors, threshold);
    case 3:
        return enlarged_support_contains_at_least_n_loc_points(basis, basis_fun_index, enlargement, n_loc);
    default:
      std::cerr << "Unknown option to should_be_refined." << std::endl;
      return false;
    }
}



/** \brief basis functions marked as to be refined on the interface of a periodic surface.
 * TODO: functions on the interface are added only if they are active on the corners of a support of a basis function marked as to be refine on the interface: is it enough? 
 *          better to chek the points inside the support of the basis functions maked as to be refined on the interface.
*/
std::set<index_t> active_functions_on_the_interface(const gsTHBSplineBasis<2>&basis,
                                                     const std::vector<index_t>& marked_functions,
                                                     const int dir
                                                    )
{
    std::set<index_t> interface;
    gsMatrix<index_t> periodicMatrix;
    gsMatrix<real_t> tmpSup;
    
    for(unsigned i = 0; i < marked_functions.size(); i++)
    {
        tmpSup = basis.support(marked_functions[i]);
        gsMatrix<real_t> point(2,2);
        
        if(almostEqual(0, tmpSup(dir,0)))
        {
            interface.insert(marked_functions[i]);
            point << tmpSup(0,0), tmpSup(1,0), tmpSup(0,dir), tmpSup(1,1-dir);
        }
        else if(almostEqual(1, tmpSup(dir,1)))
        {
            interface.insert(marked_functions[i]);
            point << tmpSup(0,1), tmpSup(1,1), tmpSup(0,1-dir), tmpSup(1-dir,dir);
        }
        else
            continue;
        
        periodicMatrix = basis.active(point);
        for(int row = 0; row < periodicMatrix.rows(); row ++)
            for(int col = 0; col < periodicMatrix.cols(); col ++)
                interface.insert(periodicMatrix(row, col));
    }
    
    return interface;
}



/** \brief If a basis function is maked as to be refined on the interface, also its counterpart basis function on the other side of the interface is marked as to be refined.
 */
void functions_on_the_interface_with_periodicity(const gsTHBSplineBasis<2>& basis, const std::vector<index_t>& marked_functions, std::set<index_t>& interface)
{
    std::vector<index_t> marked_functions_interface = set_to_vector(interface);
    for(unsigned s = 0; s<marked_functions_interface.size(); s++)
    {
        unsigned basis_fun_index = marked_functions_interface[s];
        const int level = basis.levelOf(basis_fun_index);
        gsTensorBSplineBasis<2> tensBasis = *basis.getBases()[level];
        int tensIdx = basis.flatTensorIndexOf(basis_fun_index, level);
        gsVector<index_t, 2> compTensIdx = tensBasis.tensorIndex(tensIdx);
        const unsigned addIdx = tensBasis.index(tensBasis.component(0).size() - 1 - compTensIdx[0], compTensIdx[1]);
        int addfun = basis.flatTensorIndexToHierachicalIndex(addIdx, level);
                
        interface.insert(addfun);
    }
}

/** \brief functions marked as to be refined on the interface
 */
std::vector<index_t> get_marked_functions_on_the_interface(const gsTHBSplineBasis<2>& basis,
                                                            const std::vector<index_t>& marked_functions,
                                                            const int dir
                                                           )
{
    std::vector<index_t> marked_functions_interface;
    std::set<index_t> interface = active_functions_on_the_interface(basis, marked_functions, dir);
    functions_on_the_interface_with_periodicity(basis, marked_functions, interface);
    
    /*gsWriteParaview(basis, "check_BASIS");
    gsInfo << "Marked functions on the interface: " << std::endl;
    for(std::set<index_t>::iterator it = interface.begin(); it != interface.end(); ++ it)
        gsInfo << *it << std::endl;*/
        
    for(std::vector<index_t>::const_iterator it = marked_functions.begin(); it != marked_functions.end(); ++it)
        interface.insert(*it);
    
    std::vector<index_t> result = set_to_vector(interface);
    return result;
    
}


/** \brief functions marked as to be refined for periodic surfaces
 */
std::vector<index_t> get_marked_functions_periodic(gsTHBSplineBasis<2>& basis, 
                                                    const gsMatrix<>& parameters, 
                                                    const gsMatrix<>& errors, 
                                                    const real_t threshold,
                                                    const index_t n_loc,
                                                    const int dir,
                                                    int strategy,
                                                    const int circle,
                                                    std::vector<index_t>& enlargement)
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
        if (should_be_refined_with_enlargment(basis, *it, circle, enlargement, n_loc, parameters, errors, threshold, strategy))
            marked_functions.push_back(*it);
    }
    
    std::vector<index_t> marked_functions_interface = get_marked_functions_on_the_interface(basis, marked_functions, dir);

    return marked_functions_interface;
}


/** \brief Hierarchical refinement for periodic surfaces
 */
void scatteredHFittingRefine_periodic(gsTHBSplineBasis<2>& basis, 
                                      const gsMatrix<>& parameters, 
                                      const gsMatrix<>& errors, 
                                      const real_t threshold,
                                      const index_t n_loc,
                                      const int dir,
                                      int strategy,
                                      const int circle,
                                      std::vector<index_t>& enlargement)

{
    
    std::cout << "Strategy: " << strategy << std::endl;

    // basis.uniformRefine();
    std::vector<index_t> marked_functions = get_marked_functions_periodic(basis, parameters, 
                                                                           errors, threshold, 
                                                                           n_loc, dir, strategy, circle, enlargement);
    
    
    refine_marked_functions(basis, marked_functions);
    
    // std::cout << "marked functions: ";
    // std::for_each(marked_functions.begin(), marked_functions.end(), print);
    // std::cout << std::endl;
}



/** \brief functions marked as to be refined
 */
std::vector<index_t> get_marked_functions(gsTHBSplineBasis<2>& basis, 
                                           const gsMatrix<>& parameters, 
                                           const gsMatrix<>& errors, 
                                           const real_t threshold,
                                           const index_t n_loc,
                                           int strategy,
                                           const int circle,
                                           std::vector<index_t>& enlargement)
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
    
    //gsInfo << "There are " << pre_marked_functions.size() << " pre-marked functions." << std::endl;

    std::set<index_t>::iterator it;
    std::vector<index_t> marked_functions;
    for (it = pre_marked_functions.begin(); it != pre_marked_functions.end(); ++it)
    {
        if (should_be_refined_with_enlargment(basis, *it, circle, enlargement, n_loc, parameters, errors, threshold, strategy))
            marked_functions.push_back(*it);
    }
    
    //gsInfo << "There are " << marked_functions.size() << " marked functions." << std::endl;

    return marked_functions;
}

/** \brief Hierrchical refinement
 */
void scatteredHFittingRefine(gsTHBSplineBasis<2>& basis, 
                             const gsMatrix<>& parameters, 
                             const gsMatrix<>& errors, 
                             const real_t threshold,
                             const index_t n_loc,
                             int strategy,
                             const int circle,
                             std::vector<index_t>& enlargement)

{
    
    std::cout << "Strategy: " << strategy << std::endl;

    // std::cout << "errors: \n" << std::setprecision(16) 
    //           << errors.transpose() << std::endl;
    
    
    
    
    
    // basis.uniformRefine();
    std::vector<index_t> marked_functions = get_marked_functions(basis, parameters, 
                                                                  errors, threshold, 
                                                                  n_loc, strategy, circle, enlargement);
    
    
    refine_marked_functions(basis, marked_functions);
    
    // std::cout << "marked functions: ";
    // std::for_each(marked_functions.begin(), marked_functions.end(), print);
    // std::cout << std::endl;
}

}//namespace gismo
