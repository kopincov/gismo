/** @file gsL2GMapper.h

    General local to global mapper for PDE assembling

    Here there is an interface and some implementations that represent
    the act of adding the element-wise contribution of a matrix to
    the global object.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, F. Buchegger
**/

#pragma once

#include <gsMSplines/gsWeightMapper.h>
#include <gsRecipeAssembler/gsLocalToGlobal.h>

namespace gismo {

template <typename Writer>
class gsL2GMapper
        : public gsLocalToGlobalMapper<real_t>
{
protected:
    typedef typename WriterDestType<Writer>::Destination Wr;
    typedef typename WriterDestType<Writer>::Argument    ArgT;

    const gsWeightMapper<real_t> &m_tu;
    const gsWeightMapper<real_t> &m_tt;

    Wr m_writer;
public:
    gsL2GMapper(const gsWeightMapper<real_t> & tu, const gsWeightMapper<real_t> & tt, ArgT writer)
        : m_tu(tu), m_tt(tt), m_writer(writer)
    {}

    virtual void store(
            const gsMatrix<index_t> &activeTest,
            const gsMatrix<index_t> &activeUnknown,
            const gsMatrix<real_t>   &locMat
            )
    {
        for (int r=0; r < locMat.rows (); ++r)
            for (int c=0; c < locMat.cols (); ++c)
                for (typename gsWeightMapper<real_t>::Iterator row_it = m_tt.fastSourceToTarget(activeTest(r,0));row_it;++row_it)
                    for (typename gsWeightMapper<real_t>::Iterator col_it = m_tu.fastSourceToTarget(activeUnknown(c,0));col_it;++col_it)
                        m_writer.add(row_it.index(),col_it.index(),row_it.weight()*col_it.weight()*locMat(r,c));
    }
};

template <typename Writer>
class gsL2GMapperRhs
        : public gsLocalToGlobalMapper<real_t>
{
protected:
    typedef typename WriterDestType<Writer>::Destination Wr;
    typedef typename WriterDestType<Writer>::Argument    ArgT;

    const gsWeightMapper<real_t> &m_tt;

    Wr m_writer;
public:
    gsL2GMapperRhs(const gsWeightMapper<real_t> & tt, ArgT writer)
        : m_tt(tt), m_writer(writer)
    {}

    virtual void store(
            const gsMatrix<index_t> &activeTest,
            const gsMatrix<index_t> &activeUnknown,
            const gsMatrix<real_t>   &locMat
            )
    {
        for (int r=0; r < locMat.rows (); ++r)
            for (int c=0; c < locMat.cols (); ++c)
                for (typename gsWeightMapper<real_t>::Iterator row_it = m_tt.fastSourceToTarget(activeTest(r,0));row_it;++row_it)
                    m_writer.add(row_it.index(),c,row_it.weight()*locMat(r,c));
    }
};

/**
    \brief Simply add the matrix, no index conversion.

    \tparam T          type of the scalar coefficients for the local matrix
    \tparam MatrixT    type of the destination matrix
*/
template <typename MatrixT=gsSparseMatrix<> >
class gsL2GPlain:
        public gsLocalToGlobalMapper<typename MatrixT::Scalar>
{
private:
    MatrixT     &m_dest;
    index_t      m_shift_u;
    index_t      m_shift_v;
public:
    gsL2GPlain (MatrixT &storage,index_t shift_u=0,index_t shift_v=0)
        : m_dest (storage),m_shift_u(shift_u),m_shift_v(shift_v)
    {}
    virtual void store(
            const gsMatrix<index_t> &activeTest,
            const gsMatrix<index_t> &activeUnknown,
            const gsMatrix<typename MatrixT::Scalar>        &locMat
            )
    {
        m_dest.block(m_shift_u,m_shift_v,locMat.rows(),locMat.cols())+=locMat;
    }
};


}
