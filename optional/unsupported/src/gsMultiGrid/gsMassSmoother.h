/** @file gsMassSmoother.h

    @brief Provides Multigrid smoothers based on mass matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, C. Hofreither
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsMultiPatch.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsSolver/gsLinearOperator.h>
#include <gsMultiGrid/gsMultiPatchPreconditioners.h>
#include <gsSolver/gsPatchPreconditionersCreator.h>

namespace gismo {

/// @brief Implementation of a weighted mass smoother for tensor-B-splines
///
/// \ingroup Solver
GISMO_EXPORT
gsLinearOperator<>::Ptr makeMassSmootherOperator(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc);


/// @brief Implementation of arXiv:1512.07091 for one and two dimensional tensor-B-splines
///
/// \ingroup Solver
GISMO_EXPORT
gsLinearOperator<>::Ptr makeBoundaryCorrectedMassSmootherOperator(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc);

/// @brief Implementation of G+S-Report 45/2015 for tensor-B-splines
/// The smoother is adjusted to the problem \f$-\Delta u + \alpha u = f\f$.
/// The parameter damping is used for the inverse inequality: \f$ K_0 \rightarrow 1/{damping} M_0 \f$
/// The parameter truncate truncates the boundary contributions (using the exponential decay of the inverse of the mass matrix)
///
/// \ingroup Solver
inline gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherOperator(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc, real_t alpha = 1.)
{
    return gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(basis,bc,gsAssembler<>::defaultOptions(),damping,alpha);
}

/// @brief Implementation of G+S-Report 45/2015 for tensor-B-splines
/// The smoother is adjusted to the problem \f$-\Delta u + \alpha u = f\f$.
/// The parameter damping is used for the inverse inequality: \f$ K_0 \rightarrow 1/{damping} M_0 \f$
/// The parameter truncate truncates the boundary contributions (using the exponential decay of the inverse of the mass matrix)
///
/// \ingroup Solver
inline std::vector< gsLinearOperator<>::Ptr > makeSubspaceCorrectedMassSmootherOperators( const gsMultiBasis<>& mb, real_t damping, const gsBoundaryConditions<>& bc, real_t alpha = 1.) {
    const index_t nBases = mb.nBases();
    std::vector< gsLinearOperator<>::Ptr > result(nBases);
    for (index_t j=0; j<nBases; ++j)
    {
        gsBoundaryConditions<> localbc;
        bc.getConditionsForPatch(j, localbc);
        result[j] = makeSubspaceCorrectedMassSmootherOperator(mb[j], damping, localbc, alpha);
    }
    return result;
}

/// @brief Implementation of G+S-Report 45/2015 for tensor-B-splines
/// The smoother is adjusted to the problem \f$-\Delta u + \alpha u = f\f$.
/// The parameter damping is used for the inverse inequality: \f$ K_0 \rightarrow 1/{damping} M_0 \f$
/// The parameter truncate truncates the boundary contributions (using the exponential decay of the inverse of the mass matrix)
///
/// \ingroup Solver
inline std::vector< gsLinearOperator<>::Ptr > makeSubspaceCorrectedMassSmootherOperatorsDirichlet( const gsMultiBasis<>& mb, real_t damping, real_t alpha = 1. ) {
    const index_t nBases = mb.nBases();
    std::vector< gsLinearOperator<>::Ptr > result(nBases);
    for (index_t j=0; j<nBases; ++j)
    {
        index_t d = mb[j].dim();

        gsBoundaryConditions<> localbc;
        for( index_t ps=0; ps < 2*d; ++ps )
            localbc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

        result[j] = makeSubspaceCorrectedMassSmootherOperator(mb[j], damping, localbc, alpha); //TODO: CHECK IF THIS NEEDS IMPROVEMENT (DIRBC)
    }
    return result;
}


/// @brief Implementation of G+S-Report 45/2015 for tensor-B-splines
/// The smoother is adjusted to the problem \f$-\Delta u + \alpha u = f\f$.
/// The parameter damping is used for the inverse inequality: \f$ K_0 \rightarrow 1/{damping} M_0 \f$
/// The parameter truncate truncates the boundary contributions (using the exponential decay of the inverse of the mass matrix)
///
/// This variant makes use of a 1D approximation of the domain
///
/// \ingroup Solver
GISMO_EXPORT
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherOperator(const gsGeometry<>& domain, const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc, real_t alpha = 1., bool truncate = false);

/// @brief Implementation of G+S-Report 45/2015 for tensor-B-splines
/// The smoother is adjusted to the problem \f$-\Delta u + \alpha u = f\f$.
/// The parameter damping is used for the inverse inequality: \f$ K_0 \rightarrow 1/{damping} M_0 \f$
/// The parameter truncate truncates the boundary contributions (using the exponential decay of the inverse of the mass matrix)
///
/// This variant makes use of a 1D approximation of the domain
///
/// \ingroup Solver
inline std::vector< gsLinearOperator<>::Ptr > makeSubspaceCorrectedMassSmootherOperators( const gsMultiPatch<>& domain, const gsMultiBasis<>& mb, real_t damping, const gsBoundaryConditions<>& bc, real_t alpha = 1., bool truncate = false ) {
    const index_t nBases = mb.nBases();
    std::vector< gsLinearOperator<>::Ptr > result(nBases);
    for (index_t j=0; j<nBases; ++j)
    {
        gsBoundaryConditions<> localbc;
        bc.getConditionsForPatch(j, localbc);
        result[j] = makeSubspaceCorrectedMassSmootherOperator(domain.patch(j), mb[j], damping, localbc, alpha, truncate);
    }
    return result;
}

/// @brief The smoother is adjusted to the problem \f$\Delta\Delta u = f\f$.
///
/// \ingroup Solver
GISMO_EXPORT
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherOperatorBiharmonic(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc);

/// @brief The smoother is adjusted to the problem \f$\Delta\Delta u = f\f$.
///
/// This variant makes use of a 1D approximation of the domain
///
/// \ingroup Solver
GISMO_EXPORT
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherOperatorBiharmonic(const gsGeometry<>& domain, const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc);

/// @brief The smoother is adjusted to the problem \f$\Delta\Delta u = f\f$. including all terms
///
/// \ingroup Solver
GISMO_EXPORT
gsLinearOperator<>::Ptr makeSubspaceCorrectedMassSmootherOperatorFullBiharmonic(const gsBasis<>& basis, real_t damping, const gsBoundaryConditions<real_t>& bc);

/// @brief Compute a basis for S-tilde and one for its orthogonal complement
/// If \em odd is true  then the  odd-derivatives at the boundary are zero (standard).
/// If \em odd is false then the even-derivatives at the boundary are zero.
///
/// \ingroup Solver
inline void constructTildeSpaceBasis(const gsBasis<>& basis, const gsBoundaryConditions<real_t>& bc, std::vector< gsSparseMatrix<> >& B_tilde, std::vector< gsSparseMatrix<> >& B_l2compl, const bool odd = true)
{
    gsOptionList opt = gsAssembler<>::defaultOptions();
    if(!odd)
        opt.addSwitch("UseVanishingEvenDerivatives","UseVanishingEvenDerivatives",true);

    std::pair< std::vector< gsSparseMatrix<> >, std::vector< gsSparseMatrix<> > >
        result = gsPatchPreconditionersCreator<>::getTildeSpaceBasisTransformation(basis,bc,opt);

    B_tilde = give( result.first );
    B_l2compl = give( result.second );
}

} // end of namespace gismo
