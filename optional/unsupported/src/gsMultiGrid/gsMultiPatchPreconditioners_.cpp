#include<gsMultiGrid/gsMultiPatchPreconditioners.h>
#include<gsMultiGrid/gsMultiPatchPreconditioners.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsAdditivePreconditionerOp<real_t>;
CLASS_TEMPLATE_INST gsMultiplicativePreconditionerOp<real_t>;


TEMPLATE_INST
std::vector< gsSparseMatrix<real_t> > getPatchwiseTransfers(
                            const gsMultiBasis<real_t>& mb,
                            const gsBoundaryConditions<real_t>& bc,
                            const gsOptionList& opt
                        );

TEMPLATE_INST
std::vector< std::vector< gsSparseMatrix<real_t> > > getPiecewiseTransfers(
                            const gsMultiBasis<real_t>& mb,
                            const gsBoundaryConditions<real_t>& bc,
                            const gsOptionList& opt
                        );

TEMPLATE_INST
std::pair< std::vector< gsSparseMatrix<real_t> >, std::vector< gsLinearOperator<real_t>::Ptr > > setupPiecewisePreconditioner(
                            gsSparseMatrix<real_t> A,
                            const std::vector< gsLinearOperator<real_t>::Ptr >& ops,
                            const gsMultiBasis<real_t>& mb,
                            const gsBoundaryConditions<real_t>& bc,
                            const gsOptionList& opt
                        );

TEMPLATE_INST
std::vector< gsLinearOperator<real_t>::Ptr> getLocalExactSolvers(
                    const gsSparseMatrix<real_t>& A,
                    const std::vector< gsSparseMatrix<real_t> >& transfers
                    );

TEMPLATE_INST
std::vector< std::vector< std::pair< gsBasis<real_t>::Ptr, gsSparseMatrix<real_t> > > > constructPieces(
                                    const gsMultiBasis<real_t>& mb,
                                    const std::vector< gsVector<index_t> >& locals,
                                    index_t totalNumberDof,
                                    bool combineVertices
                                );

TEMPLATE_INST
std::vector< std::vector< std::pair< gsBasis<real_t>::Ptr, gsSparseMatrix<real_t> > > > constructPieces(
                                    const gsMultiBasis<real_t>& mb,
                                    const gsDofMapper& dm,
                                    bool combineVertices
                                );

TEMPLATE_INST
std::vector< std::vector< std::pair< gsBasis<real_t>::Ptr, gsSparseMatrix<real_t> > > > constructPieces(
                                    const gsMultiBasis<real_t>& mb,
                                    const gsBoundaryConditions<real_t>& bc,
                                    const gsOptionList& opt,
                                    bool combineVertices
                                );

} // namespace gismo
