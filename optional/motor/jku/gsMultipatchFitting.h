/** @file gsMultipatchFitting.h

    @brief ...

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

// This code is in experimental stage, it should be refactored, cleaned
// improved readability...

// Ordinary people should not look at it :-).




#pragma once

#include <gsCore/gsGeometry.h>
#include <gsIO/gsIOUtils.h>
//#include "gsFittingUtils.h"

namespace gismo
{

// ================================================================================
// fitting 
    
inline
void assembleSystem(const gsMultiBasis<>& basis,
                    const std::vector< gsMatrix<> >& params,
                    const std::vector< gsMatrix<> >& points,
                    const gsDofMapper& mapper,
                    gsSparseMatrix<>& A,
                    gsMatrix<>& B)
{
    for (std::size_t patch = 0; patch != basis.nBases(); patch++)
    {
        const gsMatrix<>& par = params[patch];
        const gsMatrix<>& pts = points[patch];

        const index_t numPts = pts.cols();


        for (index_t k = 0; k != numPts; k++)
        {
            const gsMatrix<>& currPar = par.col(k);

            gsMatrix<>values;
            basis[patch].eval_into(currPar, values);

            gsMatrix<index_t> actives;
            basis[patch].active_into(currPar, actives);
            const int numActive = actives.rows();

            for (index_t i = 0; i != numActive; i++)
            {
                index_t I = mapper.index(static_cast<index_t>(actives(i)),
                                         static_cast<index_t>(patch));
                B.row(I) += values(i, 0) * pts.col(k).transpose();

                for (index_t j = 0; j != numActive; j++)
                {
                    index_t J  = mapper.index(static_cast<index_t>(actives(j)),
                                              static_cast<index_t>(patch));
                    A(I, J) += values(i, 0) * values(j, 0);
                }
            }
        }
    }
}


inline
void applySmoothing(const gsMultiBasis<>& basis,
                    const gsDofMapper& mapper,
                    const real_t lambda,
                    gsSparseMatrix<>& A)
{
    for (std::size_t patch = 0; patch != basis.nBases(); patch++)
    {
        const gsBasis<>& b = basis[patch];
        const int dim = b.dim();
        const int stride = dim * (dim + 1) / 2;

        gsVector<int> numNodes(dim);

        for (int i = 0; i != dim; i++)
        {
            numNodes[i] = b.degree(i);
        }

        gsGaussRule<> quRule(numNodes);
        gsMatrix<> quNodes, der2, localA;
        gsVector<> quWeights;
        gsMatrix<index_t> actives;

        gsBasis<>::domainIter domIt = b.makeDomainIterator();

        for (; domIt->good(); domIt->next())
        {
            quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
            b.deriv2_into(quNodes, der2);

            b.active_into(domIt->center, actives);
            const index_t numActive = actives.rows();
            localA.setZero(numActive, numActive);

            for (index_t k = 0; k != quWeights.rows(); ++k)
            {
                const real_t weight = quWeights[k] * lambda;

                for (index_t i = 0; i != numActive; i++)
                {
                    for (index_t j = 0; j != numActive; j++)
                    {
                        real_t localAij = 0;

                        for (index_t s = 0; s != stride; s++)
                        {
                            if (s < dim)
                            {
                                localAij += der2(i * stride + s, k) *
                                            der2(j * stride + s, k);
                            }
                            else
                            {
                                localAij += 2 * der2(i * stride + s, k) *
                                                der2(j * stride + s, k);
                            }
                        }

                        localA(i, j) += weight * localAij;
                    }
                }
            }

            for (index_t i = 0;i != numActive; ++i)
            {
                const index_t I = mapper.index(static_cast<index_t>(actives(i)),
                                               static_cast<index_t>(patch));

                for (index_t j = 0; j != numActive; ++j)
                {
                    const index_t J = mapper.index(static_cast<index_t>(actives(j)),
                                                   static_cast<index_t>(patch));
                    A(I, J) += localA(i, j);
                }
            }
        }
    }
}

// Boundary constraints
inline
void applyConstraints(const gsMultiBasis<>& basis,
                      const gsDofMapper& mapper,
                      const std::vector< gsMatrix<> >& params,
                      const std::vector< gsMatrix<> >& points,
                      const real_t penalty,
                      gsSparseMatrix<>& A,
                      gsMatrix<>& B)
{

    const std::size_t numPatches = params.size();//basis.nBases();
    for (std::size_t patch = 0; patch != numPatches; patch++)
    {
        const gsMatrix<>& par = params[patch];
        const gsMatrix<>& pts = points[patch];

        const index_t numPts = pts.cols();

        for (index_t k = 0; k != numPts; k++)
        {
            const gsMatrix<>& currPar = par.col(k);

            gsMatrix<>values;
            basis[patch].eval_into(currPar, values);

            gsMatrix<index_t> actives;
            basis[patch].active_into(currPar, actives);
            const int numActive = actives.rows();

            for (index_t i = 0; i != numActive; i++)
            {
                index_t I = mapper.index(static_cast<index_t>(actives(i)),
                                         static_cast<index_t>(patch));
                B.row(I) += penalty * values(i, 0) * pts.col(k).transpose();

                for (index_t j = 0; j != numActive; j++)
                {
                    index_t J  = mapper.index(static_cast<index_t>(actives(j)),
                                              static_cast<index_t>(patch));
                    A(I, J) += penalty * values(i, 0) * values(j, 0);
                }
            }
        }
    }
}

inline
gsMultiPatch<> constructSolution(const gsMultiBasis<>& basis,
                                 const gsMatrix<>& globalCoefs)
{
    gsDofMapper mapper;
    basis.getMapper(true, mapper);

    gsMultiPatch<> mp;

    for (std::size_t patch = 0; patch != basis.nBases(); patch++)
    {

        const int numBasisFun = basis.size(patch);
        gsMatrix<> coefs(numBasisFun, globalCoefs.cols());
        for (int i = 0; i != numBasisFun; i++)
        {
            const int globalI = mapper.index(i, static_cast<int>(patch));
            coefs.row(i) = globalCoefs.row(globalI);
        }
	
        gsGeometry<>::uPtr t;

        try
        {
            const gsTHBSplineBasis< 2 >& b =
                    dynamic_cast< const gsTHBSplineBasis< 2 >& >(basis[patch]);
            t = gsGeometry<>::uPtr(new gsTHBSpline< 2 >(b, coefs));
        }
        catch (const std::bad_cast& e) {
            try
            {
                const gsTensorBSplineBasis< 2 >& b =
                        dynamic_cast< const gsTensorBSplineBasis< 2 >& >(basis[patch]);
                t = gsGeometry<>::uPtr(new gsTensorBSpline< 2 >(b, coefs));
            }
            catch (const std::bad_cast& e) {
                try
                {
                    const gsTHBSplineBasis< 3 >& b =
                            dynamic_cast< const gsTHBSplineBasis< 3 >& >(basis[patch]);
                    t = gsGeometry<>::uPtr(new gsTHBSpline< 3 >(b, coefs));
                }
                catch(const std::bad_cast& e) {
                    const gsTensorBSplineBasis< 3 >& b =
                            dynamic_cast< const gsTensorBSplineBasis< 3 >& >(basis[patch]);
                    t = gsGeometry<>::uPtr(new gsTensorBSpline< 3 >(b, coefs));
                }
            }
        }
	
        mp.addPatch(give(t));
    }

    mp.computeTopology();

    return mp;
}

inline
void fittingMatrix(const gsMultiBasis<>& basis,
                   const std::vector< gsMatrix<> >& params,
                   const std::vector< gsMatrix<> >& points,
                   const real_t lambda,
                   gsSparseMatrix<>& A,
                   gsMatrix<>& B)
{
    gsDofMapper mapper;
    basis.getMapper(true, mapper);
    index_t dofs = mapper.size();
    const int dimension = points[0].rows();

    A.resize(dofs, dofs);
    A.setZero();
    B.setZero(dofs, dimension);

    assembleSystem(basis, params, points, mapper, A, B);

    if (0 < lambda)
    {
        applySmoothing(basis, mapper, lambda, A);
    }

}

// Constrained boundary
inline
void fittingMatrix(const gsMultiBasis<>& basis,
                   const std::vector< gsMatrix<> >& params,
                   const std::vector< gsMatrix<> >& points,
                   const real_t lambda,
                   const std::vector< gsMatrix<> >& constrParams,
                   const std::vector< gsMatrix<> >& constrPoints,
                   const real_t penalty,
                   gsSparseMatrix<>& A,
                   gsMatrix<>& B)
{
    gsDofMapper mapper;
    basis.getMapper(true, mapper);
    index_t dofs = mapper.size();
    const int dimension = points[0].rows();

    A.resize(dofs, dofs);
    A.setZero();
    B.setZero(dofs, dimension);

    assembleSystem(basis, params, points, mapper, A, B);

    if (0 < lambda)
    {
        applySmoothing(basis, mapper, lambda, A);
    }

    if (0 < penalty) // +++
    {
        applyConstraints(basis, mapper, constrParams, constrPoints, penalty, A, B);
    }

}

// ================================================================================
// refinement of multi basis
template <unsigned d>
bool isCellAlreadyInserted(const gsVector<index_t, d>& a_cell,
                           const std::vector<index_t>& cells)
{
    for (std::size_t i = 0; i != cells.size(); i+= a_cell.rows())
    {
        int commonEntries = 0;
        for (index_t col = 0; col != a_cell.rows(); col++)
        {
            if (cells[i + col] == a_cell[col])
            {
                commonEntries++;
            }
        }

        if (commonEntries == a_cell.rows())
        {
            return true;
        }
    }

    return false;
}

void append(std::vector<index_t>& boxes,
            const gsVector<index_t>& box)
{
    for (index_t row = 0; row != box.rows(); row++)
    {
        boxes.push_back(box[row]);
    }
}

template<unsigned d>
void appendBox(std::vector<index_t>& boxes,
               std::vector<index_t>& cells,
               gsTHBSplineBasis< d >& basis,
               const gsVector<>& parameter,
               const std::vector<index_t>& ext)
{
    const int maxLvl = basis.maxLevel();
    const gsTensorBSplineBasis< d >& tensorBasis =
                *(basis.getBases()[maxLvl]);

    gsVector<index_t, d> a_cell;
    for (unsigned dim = 0; dim != d; dim++)
    {
        const gsKnotVector<>& kv = tensorBasis.component(dim).knots();
        a_cell(dim) = kv.uFind(parameter(dim)).uIndex();
    }

    if (!isCellAlreadyInserted<d>(a_cell, cells))
    {
        append(cells, a_cell);

        // get level of a cell
        gsVector<index_t, d> a_cell_up = a_cell + gsVector<index_t, d>::Ones();
        const int cell_lvl = basis.tree().query3(a_cell, a_cell_up, maxLvl) + 1;

        gsVector<index_t> box(2 * d + 1);
        box(0) = cell_lvl;
        for (unsigned dim = 0; dim != d; dim++)
        {
            unsigned lowIndex = 0;
            if (cell_lvl < maxLvl)
            {
                const unsigned shift = maxLvl - cell_lvl;
                lowIndex = (a_cell(dim) >> shift);
            }
            else
            {
                const unsigned shift = cell_lvl - maxLvl;
                lowIndex = (a_cell(dim) << shift);
            }

            // appply extensions
            const unsigned numBreaks = basis.numBreaks(cell_lvl, dim) - 1;

            const unsigned uext = static_cast<unsigned>(ext[dim]);
            const unsigned low = (lowIndex > uext) ? (lowIndex - uext): 0;
            const unsigned upp = (lowIndex + uext + 1 < numBreaks) ?
                        (lowIndex + uext + 1): numBreaks;

            box[1 + dim] = low;
            box[1 + 2 + dim] = upp;
        }
        append(boxes, box);
    }

}

template<unsigned d>
std::vector<index_t> getBoxes(gsTHBSplineBasis< d >& basis,
                               const gsMatrix<>& params,
                               const gsVector<>& errors,
                               const real_t threshold,
                               const std::vector<index_t>& ext)
{
    std::vector<index_t> cells;
    std::vector<index_t> boxes;

    for (index_t index = 0; index != errors.rows(); index++)
    {
        if (threshold <= errors(index))
        {
            appendBox<d>(boxes, cells, basis, params.col(index), ext);
        }
    }

    return boxes;
}

template<unsigned d>
void refinePatch(gsBasis<>& basis_in,
                 const gsMatrix<>& params,
                 const gsVector<>& errors,
                 const real_t threshold,
                 const std::vector<int>& ext)
{
    gsTHBSplineBasis< d >& basis = dynamic_cast< gsTHBSplineBasis< d >& >(basis_in);

    std::vector<index_t> boxes = getBoxes<d>(basis, params, errors, threshold, ext);
    basis.refineElements(boxes);

    std::cout << "  Inserted " << boxes.size() / 5 << " boxes. " << std::endl;
}

template<unsigned d>
void refine(gsMultiBasis<>& basis,
            const std::vector< gsMatrix<> >& bParams,
            const std::vector< gsVector<> >& errors,
            const real_t threshold,
            const std::vector<int>& ext)
{
    for (std::size_t patch = 0; patch != basis.nBases(); patch++)
    {
        refinePatch<d>(basis[patch], bParams[patch], errors[patch], threshold, ext);
    }
}

// ================================================================================
// fixing interfaces    
    
void repairInterfaces(gsMultiBasis<>& basis)
{
    const gsBoxTopology& topology = basis.topology();
    std::vector< boundaryInterface > interfaces = topology.interfaces();
    basis.repairInterfaces(interfaces);
}

// ================================================================================
// errors
    
void computePatchErrors(const gsGeometry<>& geom,
                        const gsMatrix<>& bParams,
                        const gsMatrix<>& bPoints,
                        gsVector<>& errors)
{
    gsMatrix<> eval;
    geom.eval_into(bParams, eval);

    errors.resize(bParams.cols());
    gsVector<> p;

    for (index_t col = 0; col != eval.cols(); col++)
    {
        p = eval.col(col) - bPoints.col(col);
        errors(col) = p.norm();
    }
}


void computeErrors(const gsMultiPatch<>& solution,
                   const std::vector< gsMatrix<> >& bParams,
                   const std::vector< gsMatrix<> >& bPoints,
                   std::vector< gsVector<> >& errors)
{
    for (std::size_t i = 0; i != bParams.size(); i++)
    {
        gsVector<> patchErrors;
        computePatchErrors(solution[i], bParams[i], bPoints[i], patchErrors);
        errors.push_back(patchErrors);
    }
}

void print_tmp(const real_t error,
               const int belowThreshold,
               const int all)
{
    std::cout << "  Fitting error: " << error << "\n"
              << "  Below threshold: " << (belowThreshold * 1.0) / all << std::endl;
}
    
real_t analyseErrors(const std::vector< gsVector<> >& errors,
                     const real_t threshold,
                     const int verbose)
{
    int belowThreshold = 0;
    int all = 0;
    real_t maxError = 0;

    for (std::size_t i = 0; i != errors.size(); i++)
    {
        int belowThresholdI = 0;
        const gsVector<>& errorsI = errors[i];
        int allI = errorsI.rows();
        real_t maxI = 0;

        for (index_t row = 0; row != allI; row++)
        {
            const real_t err = errorsI(row);

            if (err < threshold)
            {
                belowThresholdI++;
            }

            if (maxI < err)
            {
                maxI = err;
            }
        }

        if (verbose == 2)
        {
            print_tmp(maxI, belowThresholdI, allI);
        }

        belowThreshold += belowThresholdI;
        all += allI;

        if (maxError < maxI)
        {
            maxError = maxI;
        }
    }

    if (1 <= verbose)
    {
        print_tmp(maxError, belowThreshold, all);
    }

    return maxError;
}

// ================================================================================
// saving

inline
void saveDataSolution(const gsMultiPatch<>& solution,
                      const std::string& output)
{
    gsFileData<> fd;
    fd << solution;

    const std::string out = output + "_Multipatch";
    std::cout << "Saving multipatch to: " << out << std::endl;
    fd.dump(out);

    std::cout << std::endl;
}

inline
void saveDataMeshAndGeometry(const gsMultiPatch<>& solution,
                             const std::string& output)
{
    for (std::size_t patch = 0; patch != solution.nPatches(); patch++)
    {
        const gsGeometry<> & geom = solution.patch(patch);
        gsMesh<> mesh;

        makeMesh<>(geom.basis(), mesh, 3);
        geom.evaluateMesh(mesh);

        std::string out = output + "_Mesh_" + util::to_string(patch);
        std::cout << "Saving mesh to: " << out << std::endl;
        gsWriteParaview(mesh, out);
    }
    std::cout << std::endl;

    std::string out = output + "_Geometry_";
    std::cout << "Saving geometry to: " << out << "\n" << std::endl;
    gsWriteParaview(solution, out);


}

/*

// Weight function - w : [0, 1]^1 -> [1, c]
// w(u, v; c, p) = -(c - 1)*(2^(4*p)) (u^p)*(1 - u)^p*(v^p)*(1 - v)^p + c
real_t weight_function(const real_t u, const real_t v, const real_t c, const int pow)
{
    return -(c-1)*(math::pow(2, 4*pow))*math::pow(u, pow)*math::pow(1-u, pow)*math::pow(v, pow)*math::pow(1-v, pow)+c;
}

gsMatrix<> weights(const gsMatrix<>& params,
                   const real_t c, const int pow)
{
    const index_t numPoints = params.cols();
    gsMatrix<> weights(1, numPoints);
    weights.setZero();

    for (int i = 0; i != numPoints; i++)
    {
        weights(0, i) = weight_function(params(0, i), params(1, i), c, pow);
    }

    return weights;
}

void getWeights(std::vector< gsMatrix<> >& weights,
                const std::string& filename,
                const real_t c, const int pow)
{
    const int num_of_patches = 3;//number_of_patches(filename);

    gsFileData<> fd;

    for (int i = 0; i != num_of_patches; i++)
    {
        const int id_par = 2 * i;

        gsMatrix<> par;
        fd.getId< gsMatrix<> >(id_par, par);

        gsMatrix<> w = weights(par, c, pow);
        weights.push_back(w);
    }
}
*/
}

