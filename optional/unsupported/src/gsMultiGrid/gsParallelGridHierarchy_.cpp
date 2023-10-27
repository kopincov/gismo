/** @file gsParallelGridHierarchy_.cpp

    @brief Coarsening algorithms for knot vectors and bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_MPI
#include <gsMultiGrid/gsParallelGridHierarchy.hpp>
namespace gismo
{

CLASS_TEMPLATE_INST gsParallelGridHierarchy<real_t>;

TEMPLATE_INST
gsCoarseSolverAdapter<real_t>::Ptr constructCoarseSolver<real_t>(
        const std::vector<gsMatrixOp<gsSparseMatrix<real_t> >::Ptr >& mat_coarse,
        const gsParallelGlobalLocalHandler::Ptr& handler_coarse,
        const gsSortedVector<size_t>& myPatches,
        const std::vector<gsDofMapper>& patchLocalMappers,
        const gsDofMapper& globalMapper);

TEMPLATE_INST
gsCoarseSolverAdapter<real_t>::Ptr constructCoarseSolver<real_t>(
        const gsMatrixOp<gsSparseMatrix<real_t> >::Ptr& mat_coarse,
        const gsParallelGlobalLocalHandler::Ptr& handler_coarse);


} // namespace gismo

#endif
