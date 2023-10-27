/** @file gsParallelMultiPatchSmoothers.cpp

    @brief Provides Multigrid smoothers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/


#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_MPI
#include <gsMultiGrid/gsParallelMultiPatchPreconditioners.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsParallelAdditivePreconditionerOp<real_t>;


TEMPLATE_INST
std::vector< std::vector< gsSparseMatrix<real_t> > > getPiecewiseTransfers(
        const gsMultiPatch<real_t>& mp,
        const gsMultiBasis<real_t>& mb,
        const gsBoundaryConditions<real_t>& bc,
        const gsOptionList& opt,
        const gsSortedVector<size_t>& myPatches,
        const gsParallelGlobalLocalHandler& A_handler,
        const gsMpiComm comm,
        gsSortedVector<index_t>& interfaceDofs
        );

TEMPLATE_INST
std::pair< std::vector< gsSparseMatrix<real_t> >, std::vector< typename gsLinearOperator<real_t>::Ptr > > setupPiecewisePreconditioner(
        const gsParallelOperator<real_t>& A,
        const gsParallelGlobalLocalHandler& A_handler,
        const std::vector< typename gsLinearOperator<real_t>::Ptr >& ops,
        const gsMultiPatch<real_t>& mp,
        const gsMultiBasis<real_t>& mb,
        const gsBoundaryConditions<real_t>& bc,
        const gsOptionList& opt,
        const gsSortedVector<size_t>& myPatches,
        const gsMpiComm comm
        );


} // namespace gismo

#endif
