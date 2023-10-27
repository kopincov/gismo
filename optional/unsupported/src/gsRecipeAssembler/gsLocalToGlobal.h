﻿/** @file gsLocalToGlobal.h

    @brief General local to global mapper for PDE assembling

    Here there is an interface and some implementations that represent
    the act of adding the element-wise contribution of a matrix to
    the global object.

    In order to reduce code duplication the base implementation
    gsL2GBase is templated over 3 policy classes: the conversion of
    local active indexes to global active indexes for the unknown
    and the test space and the write method.

    This allows us to generate many variations avoiding code duplications.
    The most common required methods correspond to explicit instantiations
    of gsL2GBase.


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsDofMapper.h>

namespace gismo {



//////////////////// Writer policies

/**
    \brief Writers are composable objects that
    
**/

template <typename MatrixT>
class gsBaseWriter
{
public:
    enum { Flags=MatrixT::Flags };
    typedef typename MatrixT::Scalar Scalar;
protected:
    MatrixT &m_dest;
public:
    gsBaseWriter(MatrixT &dest) : m_dest(dest) {}
    void add(index_t r, index_t c, Scalar v) {m_dest(r,c)+=v;}
};

template <typename T>
class gsBaseWriter<gsSparseEntries<T> >
{
public:
    enum { Flags=Eigen::ColMajor };
    typedef T Scalar;
protected:
    gsSparseEntries<T> &m_dest;
public:
    gsBaseWriter(gsSparseEntries<T> &dest) : m_dest(dest) {}
    void add(index_t r, index_t c, Scalar v) {m_dest.add(r,c,v);}
};

template < class T,  int _Options, typename I >
class gsBaseWriter<typename Eigen::Block<Eigen::SparseMatrix<T,_Options,I> > >
{
public:
    enum { Flags=Eigen::ColMajor };
    typedef T Scalar;
protected:
    typename Eigen::Block<Eigen::SparseMatrix<T,_Options,I> > m_dest;
public:
    gsBaseWriter(typename Eigen::Block<Eigen::SparseMatrix<T,_Options,I> > dest) : m_dest(dest) {}
    void add(index_t r, index_t c, Scalar v) {m_dest.coeffRef(r,c)+=v;}
};

template <typename Dest>
struct WriterDestType
{
    typedef Dest Destination;
    typedef Dest Argument;
};
template <class T, int _Rows, int _Cols, int _Options>
struct WriterDestType<gsMatrix<T,_Rows,_Cols,_Options> >
{
    typedef gsBaseWriter<gsMatrix<T,_Rows,_Cols,_Options> > Destination;
    typedef gsMatrix<T,_Rows,_Cols,_Options>& Argument;

};
template <class T, int _Options, typename _Index>
struct WriterDestType<gsSparseMatrix<T,_Options,_Index> >
{
    typedef gsBaseWriter<gsSparseMatrix<T,_Options,_Index> > Destination;
    typedef gsSparseMatrix<T>& Argument;
};


template < class T, int _Rows, int _Cols, int _Options >
struct WriterDestType<typename Eigen::Block<Eigen::Matrix<T,_Rows,_Cols,_Options> > >
{
    typedef gsBaseWriter<typename Eigen::Block<Eigen::Matrix<T,_Rows,_Cols,_Options> > >  Destination;
    typedef typename Eigen::Block<Eigen::Matrix<T,_Rows,_Cols,_Options> >                 Argument;
};

template < class T,  int _Options, typename I >
struct WriterDestType<typename Eigen::Block<Eigen::SparseMatrix<T,_Options,I> > >
{
    typedef gsBaseWriter<typename Eigen::Block<Eigen::SparseMatrix<T,_Options,I> > >  Destination;
    typedef typename Eigen::Block<Eigen::SparseMatrix<T,_Options,I> >           Argument;
};



/**
 \brief Do nothing
**/
template <typename ScalarT>
class gsNullWriter
{
public:
    enum {Flags = 0};
    typedef ScalarT Scalar;
public:
    static gsNullWriter<ScalarT> unique_instance;
    void add(index_t /*r*/,
             index_t /*c*/,
             ScalarT /*value*/)
    {}
};

template <typename ScalarT>
gsNullWriter<ScalarT>  gsNullWriter<ScalarT>::unique_instance;

/**
 \brief By using this writer the upper triangular part of
 the destination matrix  is unaffected.
 
**/
template <typename Writer>
class gsSymmetricWriter
{
protected:
    typedef typename WriterDestType<Writer>::Destination Dest;
    typedef typename WriterDestType<Writer>::Argument    ArgT;

    Dest m_dest;
public:
    enum {Flags = Dest::Flags};
    typedef typename Dest::Scalar Scalar;
public:
    gsSymmetricWriter (ArgT dest) : m_dest(dest) {}

    void add(index_t r, index_t c, Scalar value)
    {
        if (r>=c)
            m_dest.add(r,c,value);
    }
};

/**
 \brief This writer is intended to add Lagrange multipliers to a system matrix:
 it writes both B and B^t.
 
**/
template <typename Writer1, typename Writer2=Writer1>
class gsMultiplierWriter
{
protected:
    typedef typename WriterDestType<Writer1>::Destination Dest1;
    typedef typename WriterDestType<Writer1>::Argument    ArgT1;
    typedef typename WriterDestType<Writer2>::Destination Dest2;
    typedef typename WriterDestType<Writer2>::Argument    ArgT2;

    Dest1 m_dest1;
    Dest2 m_dest2;
public:
    enum {Flags= Dest1::Flags};
    typedef typename Dest1::Scalar Scalar;
public:
    gsMultiplierWriter (ArgT1 dest1, ArgT2 dest2) : m_dest1(dest1), m_dest2(dest2) {}
    gsMultiplierWriter (ArgT1 dest) : m_dest1(dest), m_dest2(dest) {}

    void add(index_t r, index_t c, Scalar v)
    {
        m_dest1.add(r,c,v);
        m_dest2.add(c,r,v);
    }
};


template <typename Writer>
class gsCoeffWriter
{
protected:
    typedef typename WriterDestType<Writer>::Destination Dest;
    typedef typename WriterDestType<Writer>::Argument    ArgT;
public:
    enum {Flags = Dest::Flags};
    typedef typename Dest::Scalar Scalar;
protected:
    Dest   m_dest;
    Scalar m_coef;
public:
    gsCoeffWriter (ArgT dest, Scalar coef) : m_dest(dest), m_coef(coef) {}

    void add(index_t r, index_t c, Scalar value)
    {
        m_dest.add(r,c,value*m_coef);
    }
};


template <typename Writer>
class gsBoundaryWriter
{
protected:
    typedef typename WriterDestType<Writer>::Destination Dest;
    typedef typename WriterDestType<Writer>::Argument    ArgT;
public:
    enum {Flags = Dest::Flags};
    typedef typename Dest::Scalar Scalar;
protected:
    Dest       m_dest;
    index_t    m_tshift;
    index_t    m_ushift;
public:
    gsBoundaryWriter (ArgT dest, index_t tshift,index_t  ushift)
        : m_dest(dest),m_tshift(tshift) ,m_ushift(ushift)
    {}

    void add(index_t r, index_t c, Scalar value)
    {
        if (r>=m_tshift && c >=m_ushift)
            m_dest.add(r-m_tshift,c-m_ushift,value);
    }
};


template <typename Writer1, typename Writer2>
class gsMultipleWriter
{
protected:
    typedef typename WriterDestType<Writer1>::Destination Dest1;
    typedef typename WriterDestType<Writer1>::Argument    ArgT1;
    typedef typename WriterDestType<Writer2>::Destination Dest2;
    typedef typename WriterDestType<Writer2>::Argument    ArgT2;

    typedef typename Dest1::Scalar Scalar;

    Dest1   m_dest1;
    Dest2   m_dest2;
public:
    enum {Flags = Dest1::Flags};
public:
    gsMultipleWriter (ArgT1 dest1,ArgT2 dest2) : m_dest1(dest1), m_dest2(dest2) {}

    void add(index_t r, index_t c, Scalar value)
    {
        m_dest1.add(r,c,value);
        m_dest2.add(r,c,value);
    }
};

/**
 \brief This writer shift the position written to by a fixed amount.
 
**/
template <typename Writer>
class gsShiftWriter
{
protected:
    typedef typename WriterDestType<Writer>::Destination Dest;
    typedef typename WriterDestType<Writer>::Argument    ArgT;

    Dest m_dest;
    index_t m_shift_r, m_shift_c;
public:
    enum {Flags= Dest::Flags};
    typedef typename Dest::Scalar Scalar;
public:
    gsShiftWriter(ArgT dest, index_t shift_r, index_t shift_c)
        : m_dest(dest), m_shift_r(shift_r), m_shift_c(shift_c)
    {}
    void add(index_t r, index_t c, Scalar v)
    {
        m_dest.add(r+m_shift_r,c+m_shift_c,v);
    }
};

/**
 \brief This writer intercept writes to degrees of freedom
 greater then a fixed amount and perform then to a second destination.

 This is the mechanism used to write the matrix that produces the rhs_modification
 when some degrees of freedom are fixed.
**/
template <typename Writer1, typename Writer2>
class gsMatAndRhsModWriter
{
protected:
    typedef typename WriterDestType<Writer1>::Destination Dest1;
    typedef typename WriterDestType<Writer1>::Argument    ArgT1;
    typedef typename WriterDestType<Writer2>::Destination Dest2;
    typedef typename WriterDestType<Writer2>::Argument    ArgT2;
    Dest1 m_dest1;
    Dest2 m_dest2;
    index_t m_rlimit;
    index_t m_climit;
public:
    enum {Flags= Dest1::Flags};
    typedef typename Dest1::Scalar Scalar;
public:
    gsMatAndRhsModWriter(index_t limit, ArgT1 dest1, ArgT2 dest2)
        : m_dest1(dest1), m_dest2(dest2), m_rlimit(limit),m_climit(limit)
    {}
    gsMatAndRhsModWriter(index_t r_limit, index_t c_limit, ArgT1 dest1, ArgT2 dest2)
        : m_dest1(dest1), m_dest2(dest2), m_rlimit(r_limit),m_climit(c_limit)
    {}

    void add(index_t r, index_t c, Scalar v)
    {
        if (r<m_rlimit)
        {
            if (c<m_climit)
                m_dest1.add(r,c,v);
            else
                m_dest2.add(r,c-m_climit,v);
        }
    }
};


/**
 \brief This writer intercept writes to degrees of freedom
 greater then a fixed amount and perform then to a second destination.

 This is the mechanism used to write the matrix that produces the rhs_modification
 when some degrees of freedom are fixed.
**/
template <typename Writer>
class gsAccumulatorWriter
{
protected:
    typedef typename WriterDestType<Writer>::Destination Dest;
    typedef typename WriterDestType<Writer>::Argument    ArgT;
    Dest m_dest;
    index_t m_r;
    index_t m_c;
public:
    enum {Flags= Dest::Flags};
    typedef typename Dest::Scalar Scalar;
public:
    gsAccumulatorWriter(index_t r, index_t c, ArgT dest)
        : m_dest(dest), m_r(r), m_c(c)
    {}
    void add(index_t /*r*/,
             index_t /*c*/,
             Scalar v)
    {
        m_dest.add(m_r,m_c,v);
    }
};


/**
    \brief Interface to add a local contribution to the global matrix

    gsLocalToGlobalMapper objects describe how to add the element
    contributions to the global problem matrix.

    Different techniques are coded in different sub-classes.
    The sub-classes come usually in pairs:
        -one for the main matrix
        -one for the right-hand-side
    The name of sub-classes for the right-hand-side ends in Rhs and
    they adopt a one to one correspondence of the columns.

    \tparam T type of the scalar coefficients for the local matrix

    \ingroup Assembler
**/

template < typename T = real_t >
class gsLocalToGlobalMapper{
public:
    virtual ~gsLocalToGlobalMapper()
    {}
    /**
       \brief Adds contribution from a local matrix to a destination matrix

       \param[in] activeTest    input indices for test space (rows)
       \param[in] activeUnknown input indices for solution space (columns)
       \param[in] locMat        local matrix to store to the global
     */
    virtual void store(
            const gsMatrix<index_t> &activeTest,
            const gsMatrix<index_t> &activeUnknown,
            const gsMatrix<T>        &locMat
            )
    {
        GISMO_UNUSED(activeTest);
        GISMO_UNUSED(activeUnknown);
        GISMO_UNUSED(locMat);
    }
};


/**
 *  use active index directly
 */
template < typename Writer>
class gsL2GMapperActive : public gsLocalToGlobalMapper<typename Writer::Scalar>
{
protected:
    typedef typename Writer::Scalar                      Scalar;
    typedef typename WriterDestType<Writer>::Destination Dest;
    typedef typename WriterDestType<Writer>::Argument    ArgT;

    Dest m_wr;
public:
    gsL2GMapperActive(ArgT dest)
        : m_wr(dest)
    {
    }
    /**
       \brief Adds contribution from a local matrix to a destination matrix

       \param[in] activeTest    input indices for test space (rows)
       \param[in] activeUnknown input indices for solution space (columns)
       \param[in] locMat        local matrix to store to the global
     */
    virtual void store(
            const gsMatrix<index_t> &activeTest,
            const gsMatrix<index_t> &activeUnknown,
            const gsMatrix<Scalar>   &locMat
            )
    {
        if ( (Writer::Flags & Eigen::RowMajor) == Eigen::RowMajor )
        {
            for (index_t r=0; r < locMat.rows (); ++r)
            {
                index_t rr=activeTest(r);
                for (index_t c=0; c < locMat.cols (); ++c)
                {
                    index_t cc=activeUnknown(c);
                    m_wr.add(rr,cc,locMat(r,c));
                }
            }
        }
        else if ( (Writer::Flags & Eigen::ColMajor) == Eigen::ColMajor )
        {
            for (index_t c=0; c < locMat.cols (); ++c)
            {
                index_t cc=activeUnknown(c);
                for (index_t r=0; r < locMat.rows (); ++r)
                {
                    index_t rr=activeTest(r);
                    m_wr.add(rr,cc,locMat(r,c));
                }
            }
        }
    }
};

template <typename MatrixT>
memory::shared_ptr<gsL2GMapperActive<gsBaseWriter<MatrixT> > >  getWriter(MatrixT &mat)
{
    return memory::shared_ptr<gsL2GMapperActive<gsBaseWriter<MatrixT> > >(new gsL2GMapperActive<gsBaseWriter<MatrixT> > (gsBaseWriter<MatrixT>(mat)));
}

template <typename MatrixT>
memory::shared_ptr<gsL2GMapperActive<gsCoeffWriter<MatrixT> > > getWriterWithCoef(MatrixT &mat, typename MatrixT::Scalar coef)
{
    return memory::shared_ptr<gsL2GMapperActive<gsCoeffWriter<MatrixT> > >(new gsL2GMapperActive<gsCoeffWriter<MatrixT> > (gsCoeffWriter<MatrixT>(mat,coef)));
}

template <typename MatrixT>
memory::shared_ptr<gsL2GMapperActive<gsCoeffWriter<gsMultiplierWriter<MatrixT,MatrixT> > > > getSymmetricWriterWithCoef(MatrixT &mat, typename MatrixT::Scalar coef)
{
    return memory::shared_ptr<gsL2GMapperActive<gsCoeffWriter<gsMultiplierWriter<MatrixT,MatrixT> > > >
    (new gsL2GMapperActive<gsCoeffWriter<gsMultiplierWriter<MatrixT,MatrixT> > >(gsCoeffWriter<gsMultiplierWriter<MatrixT,MatrixT> >(gsMultiplierWriter<MatrixT,MatrixT>(mat,mat),coef)));
}


template <typename MatrixT>
memory::shared_ptr<gsL2GMapperActive<gsMultiplierWriter<MatrixT,MatrixT> > > getMultiplierWriter(MatrixT &mat,MatrixT &mat2)
{
    return memory::shared_ptr<gsL2GMapperActive<gsMultiplierWriter<MatrixT,MatrixT> > >
    (new gsL2GMapperActive<gsMultiplierWriter<MatrixT,MatrixT> >(gsMultiplierWriter<MatrixT,MatrixT> (mat,mat2)));
}


template <typename MatrixT>
memory::shared_ptr<gsL2GMapperActive<gsCoeffWriter<gsMultiplierWriter<gsNullWriter<typename MatrixT::Scalar>,MatrixT> > > > getTransposeWriterWithCoef(MatrixT &mat, typename MatrixT::Scalar coef)
{
    return memory::shared_ptr<gsL2GMapperActive<gsCoeffWriter<gsMultiplierWriter<gsNullWriter<typename MatrixT::Scalar>,MatrixT> > > >
    (new gsL2GMapperActive<gsCoeffWriter<gsMultiplierWriter<gsNullWriter<typename MatrixT::Scalar>,MatrixT> > >(gsCoeffWriter<gsMultiplierWriter<gsNullWriter<typename MatrixT::Scalar>,MatrixT> >(gsMultiplierWriter<gsNullWriter<typename MatrixT::Scalar>,MatrixT>(gsNullWriter<typename MatrixT::Scalar>::unique_instance,mat),coef)));
}

template <typename MatrixT, typename MatrixT2>
memory::shared_ptr<gsL2GMapperActive<gsMatAndRhsModWriter<MatrixT,MatrixT2> > >  getSplitWriter(MatrixT &mat,MatrixT2 &rmat, index_t limit)
{
    return memory::shared_ptr<gsL2GMapperActive<gsMatAndRhsModWriter<MatrixT,MatrixT2> > >(new gsL2GMapperActive<gsMatAndRhsModWriter<MatrixT,MatrixT2> > (gsMatAndRhsModWriter<MatrixT,MatrixT2> (limit,limit,mat,rmat)));
}
template <typename MatrixT>
memory::shared_ptr<gsL2GMapperActive<gsMatAndRhsModWriter<MatrixT,gsNullWriter<typename MatrixT::Scalar> > > >  getSplitWriterRhs(MatrixT &mat, index_t limit)
{
    return memory::shared_ptr<gsL2GMapperActive<gsMatAndRhsModWriter<MatrixT,gsNullWriter<typename MatrixT::Scalar> > > >(new gsL2GMapperActive<gsMatAndRhsModWriter<MatrixT,gsNullWriter<typename MatrixT::Scalar> > > (gsMatAndRhsModWriter<MatrixT,gsNullWriter<typename MatrixT::Scalar> > (limit,INT_MAX,mat,gsNullWriter<typename MatrixT::Scalar>::unique_instance)));
}


/**
 *  write according to the position in a vector
 */
template < typename Writer>
class gsL2GMapperFilter : public gsLocalToGlobalMapper<typename Writer::Scalar>
{
protected:
    typedef typename Writer::Scalar                      Scalar;
    typedef typename WriterDestType<Writer>::Destination Dest;
    typedef typename WriterDestType<Writer>::Argument    ArgT;

    Dest m_wr;
    const std::vector<index_t> &m_rows;
    const std::vector<index_t> &m_cols;
public:
    gsL2GMapperFilter(ArgT dest,const std::vector<index_t> &whatR,const std::vector<index_t> &whatC)
        : m_wr(dest),m_rows(whatR),m_cols(whatC)
    {
    }

    virtual void store(
            const gsMatrix<index_t> &activeTest,
            const gsMatrix<index_t> &activeUnknown,
            const gsMatrix<Scalar>   &locMat
            )
    {
        std::vector<index_t>::const_iterator Rbeg=m_rows.begin();
        std::vector<index_t>::const_iterator Rit=Rbeg;
        std::vector<index_t>::const_iterator Rend=m_rows.end();

        std::vector<index_t>::const_iterator Cbeg=m_cols.begin();
        std::vector<index_t>::const_iterator Cit=Cbeg;
        std::vector<index_t>::const_iterator Cend=m_cols.end();

        for (index_t r=0; r< activeTest.rows();++r)
        {
            Rit=std::lower_bound(Rit,Rend,static_cast<index_t>(activeTest(r)));
            if (*Rit!=static_cast<index_t>(activeTest(r)))
                continue;
            Cit=Cbeg;
            for (index_t c=0; c< activeUnknown.rows();++c)
            {
                Cit=std::lower_bound(Cit,Cend,
                                     static_cast<index_t>(activeUnknown(c)));
                if (*Cit!=static_cast<index_t>(activeUnknown(c)))
                    continue;
                m_wr.add(Rit-Rbeg,Cit-Cbeg,locMat(r,c));
            }
        }

    }
};

template < typename Writer>
class gsL2GMapperFilterRow : public gsLocalToGlobalMapper<typename Writer::Scalar>
{
protected:
    typedef typename Writer::Scalar                      Scalar;
    typedef typename WriterDestType<Writer>::Destination Dest;
    typedef typename WriterDestType<Writer>::Argument    ArgT;

    Dest m_wr;
    const std::vector<index_t> &m_rows;
public:
    gsL2GMapperFilterRow(ArgT dest,const std::vector<index_t> &what)
        : m_wr(dest),m_rows(what)
    {
    }

    virtual void store(
            const gsMatrix<index_t> &activeTest,
            const gsMatrix<index_t> &activeUnknown,
            const gsMatrix<Scalar>   &locMat
            )
    {
        std::vector<index_t>::const_iterator Rbeg=m_rows.begin();
        std::vector<index_t>::const_iterator Rit=Rbeg;
        std::vector<index_t>::const_iterator Rend=m_rows.end();


        for (index_t r=0; r< activeTest.rows();++r)
        {
            Rit=std::lower_bound(Rit,Rend,static_cast<index_t>(activeTest(r)) );
            if (*Rit!=static_cast<index_t>(activeTest(r)))
                continue;
            for (index_t c=0; c< activeUnknown.rows();++c)
                m_wr.add(Rit-Rbeg,activeUnknown(c),locMat(r,c));
        }

    }
};


template <typename MatrixT>
memory::shared_ptr<gsLocalToGlobalMapper<typename MatrixT::Scalar> > getFilteredWriter(MatrixT &mat, const std::vector<index_t> &what)    { return memory::shared_ptr<gsLocalToGlobalMapper<typename MatrixT::Scalar> >(new gsL2GMapperFilter<gsBaseWriter<MatrixT> > (mat,what,what));}

template <typename MatrixT>
memory::shared_ptr<gsLocalToGlobalMapper<typename MatrixT::Scalar> > getFilteredRowWriter(MatrixT &mat, const std::vector<index_t> &what) { return memory::shared_ptr<gsLocalToGlobalMapper<typename MatrixT::Scalar> >(new gsL2GMapperFilterRow<gsBaseWriter<MatrixT> > (mat,what));}




/**
 *  use dof mapper
 */
template < typename Writer>
class gsL2GMapperDofMapper : public gsLocalToGlobalMapper<typename Writer::Scalar>
{
protected:
    typedef typename Writer::Scalar                      Scalar;
    typedef typename WriterDestType<Writer>::Destination Dest;
    typedef typename WriterDestType<Writer>::Argument    ArgT;

    Dest m_wr;
    const gsDofMapper &m_dofT;
    const gsDofMapper &m_dofU;
    const int m_patch;
public:
    gsL2GMapperDofMapper(ArgT dest, const gsDofMapper &dofMap,int patch=0)
        : m_wr(dest), m_dofT(dofMap), m_dofU(dofMap), m_patch(patch)
    {}
    gsL2GMapperDofMapper(ArgT dest, const gsDofMapper &dofMapT,const gsDofMapper &dofMapU
                         ,int patch=0)
        : m_wr(dest), m_dofT(dofMapT), m_dofU(dofMapU), m_patch(patch)
    {}



    /**
       \brief Adds contribution from a local matrix to a destination matrix

       \param[in] activeTest    input indices for test space (rows)
       \param[in] activeUnknown input indices for solution space (columns)
       \param[in] locMat        local matrix to store to the global
     */
    virtual void store(
            const gsMatrix<index_t> &activeTest,
            const gsMatrix<index_t> &activeUnknown,
            const gsMatrix<Scalar>   &locMat
            )
    {
        if ( (Writer::Flags & Eigen::RowMajor) == Eigen::RowMajor )
        {
            for (index_t r=0; r < locMat.rows (); ++r)
            {
                index_t rr=m_dofT.index(activeTest(r),m_patch);
                for (index_t c=0; c < locMat.cols (); ++c)
                {
                    index_t cc=m_dofU.index(activeUnknown(c),m_patch);
                    m_wr.add(rr,cc,locMat(r,c));
                }
            }
        }
        else if ( (Writer::Flags & Eigen::ColMajor) == Eigen::ColMajor )
        {
            for (index_t c=0; c < locMat.cols (); ++c)
            {
                index_t cc=m_dofU.index(activeUnknown(c),m_patch);
                for (index_t r=0; r < locMat.rows (); ++r)
                {
                    index_t rr=m_dofT.index(activeTest(r),m_patch);
                    m_wr.add(rr,cc,locMat(r,c));
                }
            }
        }
    }
};





//////////////////////////////////////////////////////////////////////////////
///       OLD  COMPLEX BUT WORKING CODE
//////////////////////////////////////////////////////////////////////////////


/**
    \brief this templated class contains the code that loops over the
    local matrix entries and add them to the global matrix.

    The actual behavior of the loop is determined by 3 policy classes:
        - one for the transformation of the row indexes,
        - one for the transformation of the column indexes,
        - one for the actual addition

    The classes implementing the transformation of the indexes must
    contain a method

    unsigned operator() (unsigned \em pos, const gsMatrix<unsigned> &\em active)

    pos is the row in the matrix, active contains the local \em active function indexes
    so that the row \em pos of the local matrix correspond to the function active(\em pos).
    See gsIdentityMapping, gsActiveMapping and gsDOFMappedMapping.


    The classes implementing the writing must contain a method

    void add(index_t \em r, index_t \em c, MatrixT::Scalar \em v)

    where \em r and \em c are the indexes of the coefficient to write to
    and \em v is the value to add to the global matrix.
    Writer classes are templated on the destination so that they can be chained.

    See gsDW, gsSimmetricWriter, gsSimmetricWriter and gsMultiplierWriter.

    \ingroup Assembler
**/

template <typename TransT, typename TransU, typename Writer>
class gsL2GBase
        : public gsLocalToGlobalMapper<typename Writer::Scalar>
{
protected:
    typedef typename Writer::Scalar T;
private:
    TransU m_tu;
    TransT m_tt;
    Writer m_wr;
public:
    gsL2GBase( TransT tt, TransU tu,  Writer wr)
        : m_tu(tu), m_tt(tt), m_wr(wr)
    {}
    gsL2GBase( Writer wr)
        : m_wr(wr)
    {}

    virtual void store(
            const gsMatrix<index_t> &activeTest,
            const gsMatrix<index_t> &activeUnknown,
            const gsMatrix<T>        &locMat
            )
    {
        if ( (Writer::Flags & Eigen::RowMajor) == Eigen::RowMajor )
        {
            for (index_t r=0; r < locMat.rows (); ++r)
            {
                index_t rr=m_tt(r,activeTest);
                for (index_t c=0; c < locMat.cols (); ++c)
                {
                    index_t cc=m_tu(c,activeUnknown);
                    m_wr.add(rr,cc,locMat(r,c));
                }
            }
        }
        else if ( (Writer::Flags & Eigen::ColMajor) == Eigen::ColMajor )
        {
            for (index_t c=0; c < locMat.cols (); ++c)
            {
                index_t cc=m_tu(c,activeUnknown);
                for (index_t r=0; r < locMat.rows (); ++r)
                {
                    index_t rr=m_tt(r,activeTest);
                    m_wr.add(rr,cc,locMat(r,c));
                }
            }
        }

    }
};


/**
   \brief Maps the coefficient according to the reported position in the local matrix.
   It is used for multiple rhs vectors.
**/
class gsIdentityMapping
{
public:
    unsigned operator() (unsigned v, const gsMatrix<index_t> &) {return v;}
};

/**
   \brief Maps the coefficient according to the index of the active function.
   This is what is usually taught in FEM courses.
**/
class gsActiveMapping
{
public:
    unsigned operator() (unsigned v, const gsMatrix<index_t> &active) {return active(v,0);}
};

/**
   \brief Maps the coefficient according to a gsDofMapper.
   As above it uses the index of the active function, but
   perform an additional steps (done by the dof mapper)
   to identify common degrees of freedom on patch interfaces
   and to move fixed degrees of freedom to the last indexes.
   This allow to eliminate Dirichlet degrees of freedom.
**/
template <typename T>
class gsDOFMappedMapping
{
private:
    const gsDofMapper &m_dof_mapper;
    const unsigned        m_patch_id;
public:
    gsDOFMappedMapping(const gsDofMapper &dof_mapper, unsigned patch_id)
        : m_dof_mapper(dof_mapper), m_patch_id( patch_id)
    {}
    unsigned operator() (unsigned v, const gsMatrix<index_t> &active) {return m_dof_mapper.index(active(v,0), m_patch_id);}
};




/**
        \brief Uses the active indexes as the indexes of the global matrix.

        Adds the coefficient (i,j) of the local matrix to the coefficient
        (r,s) of the global matrix where
            -s is the active function i of the test space on the element
            -t is the active function j of the unknown space on the element

        \tparam T          type of the scalar coefficients for the local matrix
        \tparam MatrixT    type of the destination matrix
**/
template <typename MatrixT=gsSparseMatrix<> >
class gsL2GActives :  public gsL2GBase<gsActiveMapping,gsActiveMapping, typename WriterDestType<MatrixT>::Destination >
{
protected:
    typedef typename WriterDestType<MatrixT>::Destination Dest;
    typedef gsL2GBase<gsActiveMapping,gsActiveMapping, Dest> BaseClass;
public:
    gsL2GActives(MatrixT &dest) : BaseClass(dest) {}
};

/**
        \brief Uses the active indexes as the indexes of the global matrix.

        Adds the coefficient (i,j) of the local matrix to the coefficient
        (r,j) of the global matrix where
            -s is the active function i of the test space on the element

        \tparam T          type of the scalar coefficients for the local matrix
        \tparam MatrixT    type of the destination matrix
**/
template <typename MatrixT=gsMatrix<> >
class gsL2GActivesRhs :  public gsL2GBase<gsActiveMapping,gsIdentityMapping, typename WriterDestType<MatrixT>::Destination >
{
protected:
    typedef typename WriterDestType<MatrixT>::Destination Dest;
    typedef gsL2GBase<gsActiveMapping,gsIdentityMapping, Dest> BaseClass;
public:
    gsL2GActivesRhs(MatrixT &dest) : BaseClass(dest) {}
};


/**
        \brief Adds a local matrix and its transpose as a Lagrange Multiplier

        The coefficient (i,j) of the local matrix is added to the coefficients
        (vS+s,cS+j) and (cS+j,vS+s) of the global matrix where:
            - s is the active function i of the test space on the element
            - vS is the constructor argument varShift
            - cS is the constructor argument conShift
        If the local matrix is A the matrix added to the global looks like
        0  0  0  0
        0  0  0  A
        0  0  0  0
        0  At 0  0
        where the varShift is the starting row of A in the above matrix and
        con Shift the starting column.

        \tparam T          type of the scalar coefficients for the local matrix
        \tparam MatrixT    type of the destination matrix
**/
template <typename MatrixT=gsSparseMatrix<> >
class gsL2GActiveMultiplier :  public gsL2GBase<
        gsActiveMapping,
        gsActiveMapping,
        gsMultiplierWriter<gsShiftWriter<MatrixT>,gsShiftWriter<MatrixT> >
        >
{
public:
    typedef gsShiftWriter<MatrixT> SW;
    typedef gsMultiplierWriter<SW,SW> MW;
    typedef gsL2GBase<gsActiveMapping,gsActiveMapping,MW> BaseClass;
    gsL2GActiveMultiplier(
            MatrixT &dest,
            unsigned start_r,
            unsigned start_c)
        :   BaseClass(MW(SW(dest, start_r,start_c),SW(dest, start_c,start_r))) {}
    gsL2GActiveMultiplier(
            MatrixT &dest,
            unsigned start_r,
            unsigned start_c,
            unsigned start_rt,
            unsigned start_ct
            )
        : BaseClass(MW(SW(dest, start_r,start_c),SW(dest, start_rt,start_ct))) {}

};


/**
        \brief Uses two gsDofMapper to get the indexes of the global matrix

        This L2G supports elimination of Dirichlet degrees of freedom and gluing
        of the common degrees of freedom on patch interfaces.
        Because of this it needs two destination matrices:
            -a global matrix M for the operator
            -a global matrix RHS_MOD that stores the linear mapping from the
             Dirichlet values to the RHS

        Given a position (i,j) in the local matrix the transformed indexes s, and t
        are obtained by the tSpaceDMap and uSpaceDMap respectively.

        Then the coefficient (i,j) can be either be:
            - ignored if s is a fixed degree of freedom
            - multiplied by the provided coefficient for the Dirichlet degrees of freedom
              and added to the s row of RHS, if s if a not fixed and t is fixed
            - added to M in position (s,t), if both s and t are not fixed degrees of freedom

        \tparam T          type of the scalar coefficients for the local matrix
        \tparam MatrixT    type of the destination matrix
    */
template <typename MatrixT1=gsSparseMatrix<> , typename MatrixT2=gsSparseMatrix<> >
class gsL2GMapped :  public gsL2GBase<
        gsDOFMappedMapping<typename MatrixT1::Scalar>,
        gsDOFMappedMapping<typename MatrixT1::Scalar>,
        gsMatAndRhsModWriter<MatrixT1 ,MatrixT2 >
        >
{
protected:
    typedef typename MatrixT1::Scalar T;
    typedef gsDOFMappedMapping<typename MatrixT1::Scalar> DOF;
    typedef gsMatAndRhsModWriter<MatrixT1 ,MatrixT2 > MW;
    typedef gsL2GBase<DOF,DOF,MW> BaseClass;
public:
    /**
            @brief Construct the object using the following parameters
            @param[in] mat_dest   is a reference to global problem matrix
            @param[in] rhs_dest   is a reference to the rhs matrix
            @param[in] tSpaceDMap is the gsDofMapper for the test basis,
                                 it corresponds to rows
            @param[in] uSpaceDMap is the gsDofMapper for the unknown basis,
                                 it corresponds to columns
            @param[in] patchID    the patch argument to pass to the gsDofMappers

    **/
    gsL2GMapped(    MatrixT1 &mat_dest,
                    MatrixT2 &rhs_dest,
                    const gsDofMapper   &tSpaceDMap,
                    const gsDofMapper   &uSpaceDMap,
                    int               patchID
                    ) : BaseClass(
                            gsDOFMappedMapping<T>(tSpaceDMap,patchID),
                            gsDOFMappedMapping<T>(uSpaceDMap,patchID),
                            gsMatAndRhsModWriter<MatrixT1,MatrixT2>(tSpaceDMap.freeSize(),uSpaceDMap.freeSize(),mat_dest,rhs_dest)
                            ) {}
};


/**
        \brief Uses a gsDofMapper to get the row index of the global matrix

        The coefficient (i,j) of the local matrix is either
            -ignored
            -added to the coefficient (s,j) of the global destination matrix
        depending on whether the index s obtained by the provided gsDofMapper
        corresponds to a Dirichlet degree of freedom or to not.

        \tparam T          type of the scalar coefficients for the local matrix
        \tparam MatrixT    type of the destination matrix
**/
template <typename MatrixT=gsMatrix<> >
class gsL2GMappedRhs :  public gsL2GBase<
        gsDOFMappedMapping<typename MatrixT::Scalar >,
        gsIdentityMapping,
        gsMatAndRhsModWriter<MatrixT,gsNullWriter<typename MatrixT::Scalar> >
        >
{
protected:
    typedef typename MatrixT::Scalar T;
    typedef gsMatAndRhsModWriter<MatrixT ,gsNullWriter<T> > MARHS;
    typedef gsL2GBase< gsDOFMappedMapping<T>, gsIdentityMapping, MARHS> BaseClass;
public:
    /**
            @brief Construct the object using the following parameters
            @param[in] mat_dest is a reference to the global matrix
            @param[in] tSpaceDMap is the gsDofMapper for the test basis,
                                  it corresponds to rows
            @param[in] patchID    the patch argument to pass to the gsDofMappers
**/
    gsL2GMappedRhs( MatrixT &mat_dest,
                    const gsDofMapper   &tSpaceDMap,
                    int                     patchID
                    ) : BaseClass(
                            gsDOFMappedMapping<T>(tSpaceDMap, patchID),
                            gsIdentityMapping(),
                            gsMatAndRhsModWriter<MatrixT,gsNullWriter<typename MatrixT::Scalar>  >(tSpaceDMap.freeSize(), INT_MAX, mat_dest,gsNullWriter<typename MatrixT::Scalar>::unique_instance)
                            ) {}
};


template <typename MatrixT1=gsSparseMatrix<> , typename MatrixT2=gsSparseMatrix<> >
class gsL2GMappedMultiplier :  public gsL2GBase<
        gsDOFMappedMapping<typename MatrixT1::Scalar >,
        gsDOFMappedMapping<typename MatrixT1::Scalar >,
        gsMultiplierWriter<
        gsMatAndRhsModWriter<gsShiftWriter<MatrixT1>,gsShiftWriter<MatrixT2> >,
        gsMatAndRhsModWriter<gsShiftWriter<MatrixT1>,gsShiftWriter<MatrixT2> > >
        >
{
protected:
    typedef typename MatrixT1::Scalar T;
    typedef gsShiftWriter<MatrixT1> SM;
    typedef gsShiftWriter<MatrixT2> SR;
    typedef gsDOFMappedMapping<T> DOF;
    typedef gsMatAndRhsModWriter<SM,SR> SpW;
    typedef gsMultiplierWriter<SpW,SpW> MW;
    typedef gsL2GBase<DOF,DOF,MW> BaseClass;
public:

    gsL2GMappedMultiplier(    MatrixT1 &mat_dest,
                              MatrixT2 &rhs_dest,
                              const gsDofMapper   &tSpaceDMap,
                              const gsDofMapper   &uSpaceDMap,
                              int               patchID,
                              unsigned          start_r,
                              unsigned          start_c
                              ) : BaseClass(
                                      gsDOFMappedMapping<T>(tSpaceDMap,patchID),
                                      gsDOFMappedMapping<T>(uSpaceDMap,patchID),
                                      MW(
                                          SpW(tSpaceDMap.freeSize(),uSpaceDMap.freeSize(),
                                              SM(mat_dest,start_r,start_c),
                                              SR(rhs_dest,start_r,0)
                                              ),
                                          SpW(
                                              uSpaceDMap.freeSize(),tSpaceDMap.freeSize(),
                                              SM(mat_dest,start_c,start_r),
                                              SR(rhs_dest,start_c,0)
                                              )
                                          )
                                      ) {}
};



} // namespace gismo
