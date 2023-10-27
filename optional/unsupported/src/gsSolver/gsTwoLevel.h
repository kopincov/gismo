/** @file gsTwoLevel.h

    @brief A two-level preconditioner/solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#pragma once

#include <gsCore/gsField.h>
#include <gsCore/gsDofMapper.h>
#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{


class gsTwoLevel : public gsPreconditionerOp<>
{

public:

    /// Shared pointer for gsTwoLevel
    typedef memory::shared_ptr<gsTwoLevel> Ptr;

    /// Unique pointer for gsTwoLevel
    typedef memory::unique_ptr<gsTwoLevel> uPtr;

    /// Base class
    typedef gsPreconditionerOp<> Base;

    /// Pointer to gsLinearOperator
    typedef gsLinearOperator<>::Ptr OpPtr;

    /// Matrix type
    typedef gsSparseMatrix<> SpMatrix;

    gsTwoLevel( SpMatrix mat,
                OpPtr smoother,
                SpMatrix Mrect,
                OpPtr massInv,
                OpPtr op
       )
       : m_mat(give(mat)), m_mat_op(makeMatrixOp(m_mat)),
       m_smoother(gsPreconditionerFromOp<>::make(makeMatrixOp(m_mat),smoother)), m_pre_smooth(1), m_post_smooth(1),
       m_cycles(1), m_Mrect(give(Mrect)), m_massInv(massInv), m_tildeEmbed(false), m_B_tilde_full(),
       m_op(op), m_coarseDamping(1)
    {}

    static uPtr make( SpMatrix mat, OpPtr smoother, SpMatrix Mrect, OpPtr massInv, OpPtr op )
        { return memory::make_unique( new gsTwoLevel( mat, smoother, Mrect, massInv, op ) ); }

    gsTwoLevel( SpMatrix mat,
                gsPreconditionerOp<>::Ptr smoother,
                SpMatrix Mrect,
                OpPtr massInv,
                OpPtr op
       )
       : m_mat(give(mat)), m_mat_op(makeMatrixOp(m_mat)),
       m_smoother(smoother), m_pre_smooth(1), m_post_smooth(1),
       m_cycles(1), m_Mrect(give(Mrect)), m_massInv(massInv), m_tildeEmbed(false), m_B_tilde_full(),
       m_op(op), m_coarseDamping(1)
    {}

    static uPtr make( SpMatrix mat, gsPreconditionerOp<>::Ptr smoother, SpMatrix Mrect, OpPtr massInv, OpPtr op )
        { return memory::make_unique( new gsTwoLevel( mat, smoother, Mrect, massInv, op ) ); }

    gsTwoLevel( SpMatrix mat,
                OpPtr smoother,
                SpMatrix Mrect,
                OpPtr massInv,
                SpMatrix B_tilde_full,
                OpPtr op
       )
       : m_mat(give(mat)), m_mat_op(makeMatrixOp(m_mat)),
       m_smoother(gsPreconditionerFromOp<>::make(makeMatrixOp(m_mat),smoother)), m_pre_smooth(1), m_post_smooth(1),
       m_cycles(1), m_Mrect(give(Mrect)), m_massInv(massInv), m_tildeEmbed(true), m_B_tilde_full(give(B_tilde_full)),
       m_op(op), m_coarseDamping(1)
    {}

    static uPtr make( SpMatrix mat, OpPtr smoother, SpMatrix Mrect, OpPtr massInv, SpMatrix B_tilde_full, OpPtr op )
        { return memory::make_unique( new gsTwoLevel( mat, smoother, Mrect, massInv, B_tilde_full, op ) ); }

    gsTwoLevel( SpMatrix mat,
                gsPreconditionerOp<>::Ptr smoother,
                SpMatrix Mrect,
                OpPtr massInv,
                SpMatrix B_tilde_full,
                OpPtr op
       )
       : m_mat(give(mat)), m_mat_op(makeMatrixOp(m_mat)),
       m_smoother(smoother), m_pre_smooth(1), m_post_smooth(1),
       m_cycles(1), m_Mrect(give(Mrect)), m_massInv(massInv), m_tildeEmbed(true), m_B_tilde_full(give(B_tilde_full)),
       m_op(op), m_coarseDamping(1)
    {}

    static uPtr make( SpMatrix mat, gsPreconditionerOp<>::Ptr smoother, SpMatrix Mrect, OpPtr massInv, SpMatrix B_tilde_full, OpPtr op )
        { return memory::make_unique( new gsTwoLevel( mat, smoother, Mrect, massInv, B_tilde_full, op ) ); }

    gsTwoLevel( OpPtr mat_op,
                OpPtr smoother,
                SpMatrix Mrect,
                OpPtr massInv,
                OpPtr op
       )
       : m_mat_op(mat_op),
       m_smoother(gsPreconditionerFromOp<>::make(mat_op,smoother)), m_pre_smooth(1), m_post_smooth(1),
       m_cycles(1), m_Mrect(give(Mrect)), m_massInv(massInv), m_tildeEmbed(false), m_B_tilde_full(),
       m_op(op), m_coarseDamping(1)
    {}

    static uPtr make( OpPtr mat_op, OpPtr smoother, SpMatrix Mrect, OpPtr massInv, OpPtr op )
        { return memory::make_unique( new gsTwoLevel( mat_op, smoother, Mrect, massInv, op ) ); }

    gsTwoLevel( OpPtr mat_op,
                gsPreconditionerOp<>::Ptr smoother,
                SpMatrix Mrect,
                OpPtr massInv,
                OpPtr op
       )
       : m_mat_op(mat_op),
       m_smoother(smoother), m_pre_smooth(1), m_post_smooth(1),
       m_cycles(1), m_Mrect(give(Mrect)), m_massInv(massInv), m_tildeEmbed(false), m_B_tilde_full(),
       m_op(op), m_coarseDamping(1)
    {}

    static uPtr make( OpPtr mat_op, gsPreconditionerOp<>::Ptr smoother, SpMatrix Mrect, OpPtr massInv, OpPtr op )
        { return memory::make_unique( new gsTwoLevel( mat_op, smoother, Mrect, massInv, op ) ); }

    gsTwoLevel( OpPtr mat_op,
                OpPtr smoother,
                SpMatrix Mrect,
                OpPtr massInv,
                SpMatrix B_tilde_full,
                OpPtr op
       )
       : m_mat_op(mat_op),
       m_smoother(gsPreconditionerFromOp<>::make(mat_op,smoother)), m_pre_smooth(1), m_post_smooth(1),
       m_cycles(1), m_Mrect(give(Mrect)), m_massInv(massInv), m_tildeEmbed(true), m_B_tilde_full(give(B_tilde_full)),
       m_op(op), m_coarseDamping(1)
    {}

    static uPtr make( OpPtr mat_op, OpPtr smoother, SpMatrix Mrect, OpPtr massInv, SpMatrix B_tilde_full, OpPtr op )
        { return memory::make_unique( new gsTwoLevel( mat_op, smoother, Mrect, massInv, B_tilde_full, op ) ); }

    gsTwoLevel( OpPtr mat_op,
                gsPreconditionerOp<>::Ptr smoother,
                SpMatrix Mrect,
                OpPtr massInv,
                SpMatrix B_tilde_full,
                OpPtr op
       )
       : m_mat_op(mat_op),
       m_smoother(smoother), m_pre_smooth(1), m_post_smooth(1),
       m_cycles(1), m_Mrect(give(Mrect)), m_massInv(massInv), m_tildeEmbed(true), m_B_tilde_full(give(B_tilde_full)),
       m_op(op), m_coarseDamping(1)
    {}

    static uPtr make( OpPtr mat_op, gsPreconditionerOp<>::Ptr smoother, SpMatrix Mrect, OpPtr massInv, SpMatrix B_tilde_full, OpPtr op )
        { return memory::make_unique( new gsTwoLevel( mat_op, smoother, Mrect, massInv, B_tilde_full, op ) ); }

public:

    void step(const gsMatrix<real_t>& rhs, gsMatrix<real_t>& x) const
    {
        gsMatrix<real_t> update, res, tmp, tmp2;
        // PRESMOOTH
        for ( index_t i=0; i<m_pre_smooth; ++i )
        {
            m_smoother->step(rhs, x);
        }

        // RAW-GRID CORRECTION
        for ( index_t i=0; i<m_cycles; ++i )
        {
            m_mat_op->apply( x, res );
            res -= rhs;
            if (m_tildeEmbed) res = m_B_tilde_full.transpose() * res;
            m_massInv->apply(res,tmp);
            tmp2.noalias() = - m_Mrect.transpose() * tmp;

            m_op->apply(tmp2,tmp);

            tmp2.noalias() = m_Mrect * tmp;
            m_massInv->apply(tmp2,update);
            if (m_tildeEmbed) update = m_B_tilde_full * update;
            x+= m_coarseDamping * update;
        }

        // POSTSMOOTH
        for ( index_t i=0; i<m_post_smooth; ++i )
        {
            m_smoother->stepT(rhs, x);
        }
    }

    void setPreSmooth( index_t i ) { m_pre_smooth = i; }
    void setPostSmooth( index_t i ) { m_post_smooth = i; }
    void setCycles( index_t i ) { m_cycles = i; }
    void setCoarseDamping( real_t tau ) { m_coarseDamping = tau; }

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addInt   ("NumPreSmooth"     , "Number of pre-smoothing steps", 1                                  );
        opt.addInt   ("NumPostSmooth"    , "Number of post-smoothing steps", 1                                 );
        opt.addInt   ("NumCycles"        , "Number of cycles (usually 1 for V-cycle or 2 for W-cycle)", 1      );
        opt.addReal  ("CoarseDamping"    , "Damping parameter for the Two-Level method", 1e-10   );
        return opt;
    }

    /// @brief Set the options based on a gsOptionList
    virtual void setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_pre_smooth       = opt.askInt   ("NumPreSmooth"     , m_pre_smooth      );
        m_post_smooth      = opt.askInt   ("NumPostSmooth"    , m_post_smooth     );
        m_cycles           = opt.askInt   ("NumCycles"        , m_cycles          );
        m_coarseDamping    = opt.askReal  ("CoarseDamping"    , m_coarseDamping   );
    }

    index_t rows() const { return m_mat_op->rows(); }
    index_t cols() const { return m_mat_op->cols(); }
    OpPtr underlyingOp() const { return m_mat_op; }

protected:

    // Matrix to be solved for
    SpMatrix m_mat; //can be stored
    OpPtr m_mat_op;

    // Smoother
    gsPreconditionerOp<>::Ptr m_smoother;
    index_t m_pre_smooth;
    index_t m_post_smooth;

    // Correction
    index_t m_cycles;
    SpMatrix m_Mrect;
    OpPtr m_massInv;
    bool m_tildeEmbed;
    SpMatrix m_B_tilde_full;
    OpPtr m_op;
    real_t m_coarseDamping;

}; // class gsTwoLevel

}  // namespace gismo
