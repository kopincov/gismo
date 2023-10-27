/** @file gsScatteredHFittingCoefficientsStiffness.h

    @brief Evaluation of Coefficients with Stiffness

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S.Imperatore
*/

#pragma once

#include <gismo.h> 

#include <gsHSplines/gsScatteredHFittingCoefficients.h>
#include <gsHSplines/gsScatteredHFittingImprovementsGeneralUtilities.h> // matrix_to_set

// Uncomment for plotting with Parasolid.
// #include <gsParasolid/gsWriteParasolid.h>

namespace gismo{

    
 /** \brief growing of the number of basis functions for the restriction consistently with \em degree
 *         In order to build a gsTensorBSplineBasis of degree \em degree, \em degree+1 basis functions are needed.
 *         If the size of \em chosenIndices is too little, the missing basis functions are added.
 */
void check_consistent_basis(std::set<index_t>& chosenIndices,
                            const size_t degree,
                            const index_t basis_component_size)
{
    
    // TODO: Consider growing in a balanced way.
    while( chosenIndices.size() < degree + 1 )
    {
        index_t tmpMax = *std::max_element(chosenIndices.begin(), chosenIndices.end());
        if( tmpMax < basis_component_size - 1 )
            chosenIndices.insert(tmpMax + 1);
        else
        {
            index_t tmpMin = *std::min_element(chosenIndices.begin(), chosenIndices.end());
            if( tmpMin > 0 )
                chosenIndices.insert(tmpMin - 1);
            else
                gsInfo << "Problem with the basis restriction: we should not get here." << std::endl;
        }
    }
}
    
/** \brief construction of gsTensorBSplineBasis restiction of gsTensorBSplineBasis \em tensBasis : check size.
 *         It can be that the active functions on the resticted data for the local fitting are not enough to build a gsTensorBSplineBasis accordingly to \em degree.
 *         If that is the case, the missing basis functions are added.       
 * 
 * \param[in]  tensBasis  const gsTensorBSpline
 * \param[in]  indicesR   const indexContainer with indices of the basis functions of \em tensBasis involved in the local fitting
 * \param[in]  degree     of the local fitting
 * \param[in]  chosen     global flat tensor index of the basis function involved in the local fitting
 * \param[out] tensChosen global tensor index corresponding to \em chosen 
 * 
 */
gsVector<index_t,2> selection_uvIndex(const gsTensorBSplineBasis<2>& tensBasis,
                                       const std::set<index_t>& indicesR,
                                       unsigned degree,
                                       index_t chosen,
                                       std::set<index_t>& selection_u,
                                       std::set<index_t>& selection_v)
{
    gsVector<index_t, 2> tmpTensIndex;
    gsVector<index_t, 2> tensChosen;
    
    for(typename std::set<index_t>::const_iterator it=indicesR.begin(); it!=indicesR.end(); ++it)
    {
        tmpTensIndex = tensBasis.tensorIndex(*it);
        selection_u.insert(tmpTensIndex(0));
        selection_v.insert(tmpTensIndex(1));
        
        if(*it == chosen)
            tensChosen << tmpTensIndex[0], tmpTensIndex[1];
    }
    
    if( (selection_u.size() < degree + 1) || (selection_v.size() < degree +1) )
    {
        check_consistent_basis(selection_u, degree, tensBasis.component(0).size());
        check_consistent_basis(selection_v, degree, tensBasis.component(1).size());
    }
    
    return tensChosen;
}



/** \brief construction of gsTensorBSplineBasis restiction of gsTensorBSplineBasis \em tensBasis : check shape.
 *         It can be that the union of the supports of the active function on the restricted data for the fitting does not have a rectangular shape.
 *         In order to deal with the G+smo implementation of gsTensorBSplineBasis, the missing basis functions are added.
 */
std::set<index_t> check_rectangle(const gsTensorBSplineBasis<2>& tensBasis,
                                   std::set<index_t>& indices)
{
    gsVector<index_t, 2> tmpIdx;
    std::set<index_t> indices_u;
    std::set<index_t> indices_v;
    std::set<index_t> check_u;
    std::set<index_t> check_v;
    std::set<index_t> check;
    
    for(std::set<index_t>::iterator it=indices.begin(); it!=indices.end(); ++it)
    {
        tmpIdx = tensBasis.tensorIndex(*it);
        indices_u.insert(tmpIdx(0));
        indices_v.insert(tmpIdx(1));
    }
    
    for(index_t u = *indices_u.begin(); u <= *(--indices_u.end()); ++u)
        check_u.insert(u);
    for(index_t v = *indices_v.begin(); v <= *(--indices_v.end()); ++v)
        check_v.insert(v);
    
    for(std::set<index_t>::iterator it = check_u.begin(); it != check_u.end(); ++it)
        for(std::set<index_t>::iterator jt = check_v.begin(); jt != check_v.end(); ++jt)
            check.insert(tensBasis.index(*it,*jt));
        
    return check;
}

/** \brief gsMatrix \em restrPar contains the parameters for the local fitting
 *         gsMatrix \em restrPt contains the points for the local fitting
 */
void restrictionOfData(const std::set<index_t>& indices,
                       const gsMatrix<>& parameters,
                       const gsMatrix<>& points,
                       gsMatrix<>& restrPar,
                       gsMatrix<>& restrPt)
{
    restrPar.resize(2, indices.size());
    restrPt.resize(3, indices.size());
    int k = 0;
    for(std::set<index_t>::const_iterator it=indices.begin(); it!=indices.end(); ++it)
    {
        restrPar.col(k) = parameters.col(*it);
        restrPt.col(k) = points.col(*it);
        k++;
    }
}


/** \brief tensor indices of the basis functions of \em tensBasis, active on the points of the local fitting.
 *         std::set \em selection_u tensor indices in \a u direction
 *         std::set \em selection_v tensor indices in \a v direction
 */
void restrictionOfBasis(const gsTensorBSplineBasis<2>& tensBasis,
                        const int tensIdx,
                        const unsigned fitting_deg,
                        const gsMatrix<>& restrPar,
                        std::set<index_t>& selection_u,
                        std::set<index_t>& selection_v,
                        gsVector<index_t,2>& selected)
{
    std::set<index_t> setR;
    
    gsMatrix<index_t> restrActive;
    tensBasis.active_into(restrPar, restrActive);
    
    std::set<index_t> set_active_fun = matrix_to_set(restrActive);
    set_active_fun.insert(tensIdx);
    setR = check_rectangle(tensBasis, set_active_fun);
    selected = selection_uvIndex(tensBasis, setR, fitting_deg, tensIdx, selection_u, selection_v);
}


/** \brief std::set of data for the local fitting.
 *         The fitting area is a circle of center \em center and radius \em radius
 */
std::set<index_t> get_data_inside_circle(const gsMatrix<>& parameters,
                                       const gsMatrix<>& center,
                                       const real_t radius)
{
    std::set<index_t> indices;
    
    for (int i = 0; i != parameters.cols(); i++)
    {
        const gsMatrix<> point = parameters.col(i);
        if ((point - center).squaredNorm() <= radius * radius)
        {
            indices.insert(i);
        }
    }
    
    return indices;
}

/** \brief std::set of data for the local fitting.
 *         The fitting area is a rectangle of center \em center and sides \em u_diff and \em v_diff
 */
std::set<index_t> get_data_inside_rectangle(const gsMatrix<>& parameters,
                                             const gsMatrix<>& center,
                                             const real_t u_diff,
                                             const real_t v_diff)
{
    std::set<index_t> indices;
    
    for (int i = 0; i != parameters.cols(); i++)
    {
        if( (std::fabs(parameters(0,i) - center(0,0)) <= u_diff) &&
            (std::fabs(parameters(1,i) - center(1,0)) <= v_diff))
        {
            indices.insert(i);
        }
    }
    
    return indices;
}

/** \brief std::set of data for the local fitting.
 *         The fitting area is defined enlarging \em support with respect to gsKnotVector \em kv0 and gsKnotvector \em kv1
 */
std::set<index_t> get_data_inside_knot_span(const gsMatrix<>& parameters,
                                             const gsMatrix<>& support,
                                             const gsMatrix<>& center,
                                             const int k,
                                             const gsKnotVector<>& kv0,
                                             const gsKnotVector<>& kv1,
                                             real_t& u_diff,
                                             real_t& v_diff)
{
    std::set<index_t> indices;
    
    if( k == 1 )
    {
        u_diff = std::fabs(support(0,0) - support(0,1)) / 2;
        v_diff = std::fabs(support(1,0) - support(1,1)) / 2;
        
        indices = get_data_inside_rectangle(parameters, center, u_diff, v_diff);
    }
    
    else
    {
        // Assuming uniform knot vectors:
        u_diff += kv0.maxIntervalLength();
        v_diff += kv1.maxIntervalLength();
        
        if(u_diff == 0 || v_diff == 0)
        {
            gsInfo << "Warning: error due to rounding";
        }
        
        indices = get_data_inside_rectangle(parameters, center, u_diff, v_diff);
    }
    
    return indices;
}

/** \brief knots inside the interval [\em left, \em right]
 */
void knots_within_iterval(gsKnotVector<real_t>& knots,
                          const real_t left,
                          const real_t right)
{
    if( left == right )
        gsInfo << "Warning: No existing interval: " << "[" << left << ", " << right << "]" << std::endl;
    
    std::vector<real_t> repeated_knots;
    for(typename gsKnotVector<real_t>::const_iterator it=knots.begin(); it!=knots.end(); ++it)
        if( (*it > left  || almostEqual(*it,left)  ) &&
            (*it < right || almostEqual(*it,right) ) )
            repeated_knots.push_back(*it);
        
    knots = gsKnotVector<real_t>(repeated_knots);
}


/** \brief knot spans in \a u and \a v direction bounding the area of the local fitting.
 */
void knots_perimeter(const gsTensorBSplineBasis<2, real_t>& tensBasis,
                     const gsMatrix<>& center,
                     const real_t half_b,
                     const real_t half_h,
                     gsKnotVector<real_t>& uKnots,
                     gsKnotVector<real_t>& vKnots)
{
    uKnots = tensBasis.component(0).knots();
    vKnots = tensBasis.component(1).knots();
    
    knots_within_iterval(uKnots, center(0,0) - half_b, center(0,0) + half_b);
    knots_within_iterval(vKnots, center(1,0) - half_h, center(1,0) + half_h);
}


/** \brief number of knot spans of the parametric domain, with respect to the area of the local fitting, containing at least one point.
 */
int parameters_inside_pieces(const gsMatrix<>& parameters,
                             const std::set<index_t>& indices,
                             gsKnotVector<real_t>& uKnots,
                             gsKnotVector<real_t>& vKnots)
{
    int pieces = 0;
    
    if( uKnots.size() == 0 || vKnots.size() == 0)
    {
        std::cerr << "Warning: error due to rounding." << std::endl;
        return pieces;
    }
        
    for(size_t u = 0; u < uKnots.size(); u++)
    {
        for(size_t v = 0; v < vKnots.size(); v++)
        {
            for(std::set<index_t>::const_iterator it = indices.begin(); it != indices.end(); ++it)
            {
                real_t x = parameters(0,*it);
                real_t y = parameters(1,*it);
                
                if((uKnots[u]<=x)&&(x<=uKnots[u+1])&&(vKnots[v]<=y)&&(y<=vKnots[v+1]))
                {
                    pieces++;
                    break;
                }
            }
        }
    }
    return pieces;
}


/** \brief density of points with respect to the area of the local fitting
 * (computed as the number of knot spans in the parametric domain coinaining at least one point over the tot number of knot span.)
 * 
 */
real_t parameters_distribution(const gsTensorBSplineBasis<2, real_t>& tensBasis,
                               const gsMatrix<>& parameters,
                               const std::set<index_t>& indices,
                               const gsMatrix<>& center,
                               const real_t half_b,
                               const real_t half_h,
                               const int basis_fun_index)
{
    gsKnotVector<real_t> uPieces, vPieces;
    knots_perimeter(tensBasis, center, half_b, half_h, uPieces, vPieces);

    real_t totPieces = (uPieces.size()-1) * (vPieces.size()-1);
    real_t numPieces = parameters_inside_pieces(parameters, indices, uPieces, vPieces);
    
    return numPieces/totPieces;
}


/** \brief std::set \em indices contains the indices necessary for the local fitting of the control point corresponding to the basis function \em basis_fun_index
 */
std::set<index_t> indices_for_fitting(const gsTHBSplineBasis<2, real_t>& basis,
                                       const int basis_fun_index,
                                       const gsMatrix<>& parameters,
                                       const int cp_magnitude,
                                       const real_t hole,
                                       const int circle,
                                       const unsigned degree,
                                       unsigned fitting_deg,
                                       std::vector<index_t>& enlargement) 
{
    std::set<index_t> indices;
    
    const gsMatrix<> support = basis.support(basis_fun_index);
    const gsMatrix<> center = get_center(support);
    real_t radius = (support.col(0) - support.col(1)).norm() / 2;
    real_t u_diff = std::fabs(support(0,0) - support(0,1)) / 2;
    real_t v_diff = std::fabs(support(1,0) - support(1,1)) / 2;
    const int level = basis.levelOf(basis_fun_index);
    const gsKnotVector<>& kv0 = basis.getBases()[level]->knots(0);
    const gsKnotVector<>& kv1 = basis.getBases()[level]->knots(1);
    const gsTensorBSplineBasis<2>& tbasis = *basis.getBases()[level];
    
    int k = 1;
    int num_points = 0;
    real_t distribPar = 0;
    while( ( num_points <= cp_magnitude ) || ( distribPar < hole ) )
    {
        if(circle == 0)
        {
            indices = get_data_inside_circle(parameters, center, radius * k);
            distribPar = 1;
        }
        else if(circle == 1)
        {
            indices = get_data_inside_rectangle(parameters, center, u_diff * k, v_diff * k );
            if( distribPar < hole )
            {
                if( u_diff == 0 || v_diff == 0 )
                {
                    gsInfo << "Warning: error due to rounding" << std::endl;
                    break;
                }
                else
                    distribPar = parameters_distribution(tbasis, parameters, indices, center, u_diff * k, v_diff * k, basis_fun_index);
            }
        }
        else
        {
            indices = get_data_inside_knot_span(parameters, support, center, k, kv0, kv1, u_diff, v_diff);
            if( distribPar < hole )
            {
                    distribPar = parameters_distribution(tbasis, parameters, indices, center, u_diff, v_diff, basis_fun_index);
            }
        }
        
        num_points = indices.size();
        k++;
    }
    
    enlargement.push_back(num_points);

    int d = degree;
    while( ( num_points < (d + 1)*(d + 2)/2 )&&( 0 < d ) )
        d -= 1;
    fitting_deg = d;

    return indices;
}

/** \brief (index) position of \em elem inside the container \em container 
*/
template <class T, class container>
size_t position(T elem,
                  const container& cont)
{
    return std::distance(cont.begin(), std::find(cont.begin(), cont.end(), elem));
}

/** \brief global tenrorIndex for the restricted basis.
 * 
 * \param[in]  selected const gsVector containing the tensorIndex of the basis function (selected[0], selected[1]) with respect to the complete gsTensorBSplineBasis the basis function belongs to
 * \param[in]  selection_u const std::set functions in the \a u direction of the restricted gsTensorBSplineBasis
 * \param[out] index global tensorIndex of the basis function indexed as (selected[0], selected[1]) with respect to the restricted gsTensorBSplineBasis 
 */
unsigned restrictionOfTensIdx(const gsVector<index_t,2>& selected,
                              const std::set<index_t>& selection_u,
                              const std::set<index_t>& selection_v)
{
    index_t selected_0 = selected(0);
    index_t selected_1 = selected(1);
    
    unsigned pos_u = position(selected_0, selection_u);
    unsigned pos_v = position(selected_1, selection_v);
    
    return pos_v * selection_u.size() + pos_u;
}

/** \brief returns gsKnotVector cointaining the knots of \em bigKV indexed between \em minFunIndex and \em maxFunIndex
 * \param[in]  bigKV gsKnotVector of knots
 * \param[in]  minFunIndex index of the first knot of the restriction
 * \param[in]  maxFunIndex index of the last knot of the restriction
 * \param[out] result gsKnotVector with restricted knots between \em minFunIndex and \em maxFunIndex
 */
gsKnotVector<real_t> restrictionOfKnots(const gsKnotVector<real_t>& bigKV,
                                        int minFunIndex,
                                        int maxFunIndex)
{
    gsKnotVector<real_t> result(bigKV.degree(),
                                bigKV.begin() + minFunIndex,
                                bigKV.end() - (bigKV.size() - maxFunIndex - 1) + bigKV.degree() + 1);
    return result;
}


/** \brief returns the restriction of \em bigTensorBasis to its functions with tensor indices \em indicesR_u in \a u direction and \em indicesR_v in \a v direction.
 * \param[in]  bigTensorBasis const gsTensorBSplineBasis from which the restricted basis is computed
 * \param[in]  indicesR_u const std::set of functions for the restriction in \a u direction
 * \param[in]  indicesR_v const std::set of functions for the restriction in \a v direction
 * \param[out] result gsTensorBSplineBasis with functions tensor-indexed by \em indicesR_u and \em indicesR_v
 */
gsTensorBSplineBasis<2, real_t> restrictionOfTensorBasis(const gsTensorBSplineBasis<2, real_t>& bigTensorBasis,
                                                         const std::set<index_t>& indicesR_u,
                                                         const std::set<index_t>& indicesR_v)
{

    gsTensorBSplineBasis<2, real_t> result(restrictionOfKnots(bigTensorBasis.component(0).knots(), *indicesR_u.begin(), *(--indicesR_u.end()) ),
                                           restrictionOfKnots(bigTensorBasis.component(1).knots(), *indicesR_v.begin(), *(--indicesR_v.end()) ));
    
    return result;
}


/** \brief computation of coefficient for the function basis \em basis_fun_index with smoothing term.
 */
gsMatrix<> get_coefficient_not_scaled_stiff(const gsTHBSplineBasis<2>& basis,
                                            const gsMatrix<>& parameters,
                                            const gsMatrix<>& points,
                                            int iteration,
                                            const int basis_fun_index,
                                            const int cp_magnitude,
                                            const real_t hole,
                                            const real_t lambda,
                                            const int circle,
                                            std::vector<index_t>& enlargement)
{
    gsMatrix<> restrPar;
    gsMatrix<> restrPt;
    gsVector<index_t,2> selected;
    gsTensorBSplineBasis<2, real_t> restrB;
    int restrTensorIndex;
    std::set<index_t> selection_u;
    std::set<index_t> selection_v;
    
    int level = basis.levelOf(basis_fun_index); 
    const short_t degree = basis.maxDegree();
    int tensIdx = basis.flatTensorIndexOf(basis_fun_index, level);
    
    short_t fitting_deg = degree;
    std::set<index_t> indices = indices_for_fitting(basis, basis_fun_index, parameters, cp_magnitude, hole, circle, degree, fitting_deg, enlargement);
   
    restrictionOfData(indices, parameters, points, restrPar, restrPt);

    gsTensorBSplineBasis<2>& tensBasis = *basis.getBases()[level];
    const int reduction = degree - fitting_deg;
    tensBasis.degreeReduce(reduction);
    
    restrictionOfBasis(tensBasis, tensIdx, fitting_deg, restrPar, selection_u, selection_v, selected);
    
    restrTensorIndex = restrictionOfTensIdx(selected, selection_u, selection_v);
    restrB = restrictionOfTensorBasis(tensBasis, selection_u, selection_v);
    
    gsFitting<real_t> fitting_fun(restrPar, restrPt, restrB);
    fitting_fun.compute(lambda);
    
    //const gsTensorBSpline<2>& resulting_geometry = *(static_cast<gsTensorBSpline<2,real_t>*>(fitting_fun.result()));
    //gsMatrix<> coefficients = resulting_geometry.coefs();
    //gsMatrix<> coefficient = coefficients.row(restrTensorIndex).transpose();
    /*extensions::gsWritePK_SHEET(resulting_geometry, "patch.xmt_txt");
    gsInfo << "Continue?" << std::endl;
    std::cin.get();*/
    
    /*if((basis_fun_index == 0) && (iteration == 0))
    {
        std::string filename="iteration" + std::to_string(iteration);
        std::string surfName= filename + "_surfaceCP" + std::to_string(basis_fun_index);
        std::string ptsName = filename + "_points" + std::to_string(basis_fun_index);
        gsWriteParaview(resulting_geometry, surfName, 1000, false, true);
        gsWriteParaviewPoints(restrPt, ptsName);
        gsWriteParaviewPoints(coefficient, filename + "CP" + std::to_string(basis_fun_index));
    }*/
    
    return fitting_fun.result()->coef(restrTensorIndex);

    //return coefficient;
}


/** \brief  computation of coefficients for THBSplineBasis \em basis.
*/
gsMatrix<> get_coefficients_not_scaled_stiff(const gsTHBSplineBasis<2>& basis, 
                                             const gsMatrix<>& parameters, 
                                             const gsMatrix<>& points,
                                             int iteration,
                                             const int cp_magnitude,
                                             const real_t hole,
                                             const real_t lambda,
                                             const int circle,
                                             std::vector<index_t>& enlargement)
{
    
    gsMatrix<> coefficients(basis.size(), points.rows());

    for (index_t k = 0; k != basis.size(); k++)
    {
        
        gsMatrix<> coefficient = get_coefficient_not_scaled_stiff(basis, parameters, points, iteration, k, cp_magnitude, hole, lambda, circle, enlargement);
        //gsInfo << "Function basis: " << k << ", coefficient: " << coefficient << std::endl;

        coefficients.row(k) = coefficient.row(0);
        /*gsMatrix<> printCoefficient = coefficient.transpose();
        gsWriteParaviewPoints(printCoefficient, "iteration" + std::to_string(iteration) + "_control point" + std::to_string(k));*/
    }
    
    return coefficients;
}

} //namespace gismo
