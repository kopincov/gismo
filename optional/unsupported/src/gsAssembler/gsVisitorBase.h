/** @file gsVisitorBase.h

    @brief Base class for all visitors, tentative implementations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, S. Kleiss, A. Mantzaflaris
*/


#include <gsCore/gsFuncData.h>
#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsDofMapper.h>
#include <gsMSplines/gsWeightMapper.h>

#include <gsCore/gsDomainIterator.h>
//#include <gsAssembler/gsLocalToGlobal.h>

namespace gismo {


template <typename T>
class gsLocalMatrices
{
    index_t                   m_patch;
    std::vector<gsMatrix<index_t> > m_active;
    std::vector<gsMatrix<T> > m_values;
};


template <typename T>
class gsVisitorBase
{
protected:
    std::vector<gsFuncData<T> > m_data;
    gsMapData<T>                m_geo;
    gsLocalMatrices<T>          m_result;
public:
    void evaluate (
            const std::vector<gsFunctionSet<T>*> &basis,
            const gsFunctionSet<T>               &geometry,
            const gsMatrix<T>                    &points,
            const gsVector<T>                    &weights
            )
    {
        geometry.compute(points,m_geo);
        for (size_t b=0;b<m_data.size();++b)
            basis[b]->compute(m_geo,m_data[b]);
        if (m_geo.measures.rows()+m_geo.measures.cols())
            m_geo.measures.array()*=weights.transpose().array();
    }

    virtual const gsLocalMatrices<T> & assemble (const gsDomainIterator<T> *element) = 0;
};


// TODO visitor interface
// TODO visitor boundary



// just a try to implement the loop
template <typename T, typename Visitor, typename Result>
class BBBassembler
{
    std::vector<gsFunctionSet<T>*> m_discretization;
    gsFunction<T>                  m_geo;

    void integrate (gsDomainIterator<T> &domIt, gsQuadRule<T> &quad, Visitor &vis, Result &res)
    {
        gsMatrix<T> qNodes;
        gsVector<T> qWeights;
        for (;domIt.good(); domIt.next() )
        {
            quad.mapTo(domIt.upperCorner(),domIt.lowerCorner(),qNodes,qWeights);
            vis.evaluate(m_discretization,m_geo,qNodes,qWeights);
            res+=vis.assemble(domIt);
        }
    }

};





template <typename T, typename Writer>
void l2gDofMap (int patch, const gsDofMapper &dof, const gsMatrix<T> &locMat, const gsMatrix<index_t> &active, WriterDestType<Writer>::Argument wrA)
{
    const WriterDestType<Writer>::Destination m_wr(wrA);
    for (index_t c=0; c < locMat.cols (); ++c)
    {
        index_t cc=dof.index(active(c),patch);
        for (index_t r=0; r < locMat.rows (); ++r)
        {
            index_t rr=dof.index(active(r),patch);
            m_wr.add(rr,cc,locMat(r,c));
        }
    }
}



// store the result in one big matrix using a dof mapper.
// works only for scalar problems at the moment
template <typename T>
class gsOneBigMatrixDofMapper
{
public:
    typedef gsSparseMatrix<T>                   SMatT;
    typedef gsMatrix<T>                         FMatT;
    typedef gsMatAndRhsModWriter<SMatT,SMatT>            SysW;
    typedef gsMatAndRhsModWriter<FMatT,gsNullWriter<T> > RhsW;
public:
    gsOneBigMatrixDofMapper(const gsMultiBasis<T> &basis, index_t rhsNum=1)
        : m_map(basis), m_sysWriter(m_map.freeSize(),m_sys,m_mod), m_rhsWriter(m_map.freeSize(),m_rhs)
    {
        initMatrices(basis);
    }

    void clear()
    {
        m_sys.setZero();
        m_mod.setZero();
        m_rhs.setZero(m_rhs.rows(), m_rhs.cols());
    }

    void initMatrices(const gsMultiBasis<T> &basis)
    {
        const index_t freeS=m_map.freeSize();
        const index_t elimS=m_map.boundarySize();
        // set function sizes
        m_sys.resize(freeS, freeS);
        m_mod.resize(freeS, elimS);
        m_rhs.resize(freeS, m_rhsNum);
        // TODO reserve
    }

    void operator+= (const gsLocalMatrices<T> &loc)
    {
        GISMO_ASSERT(loc.m_active.size()==1,"wrong number of spaces");
        GISMO_ASSERT(loc.m_values.size()==2,"expected 2 local matrices: system matrix and rhs matrix");

        l2gDofMap(loc.patch,m_map,loc.m_values[0],loc.m_actives[0],m_sysWriter);
        l2gDofMap(loc.patch,m_map,loc.m_values[1],loc.m_actives[0],m_rhsWriter);
    }

          gsSparseMatrix<T> & getMatrix()       {return m_sys;}
    const gsSparseMatrix<T> & getMatrix() const {return m_sys;}
          gsSparseMatrix<T> & getRhs()          {return m_rhs;}
    const gsSparseMatrix<T> & getRhs()    const {return m_rhs;}
          gsSparseMatrix<T> & getMod()          {return m_mod;}
    const gsSparseMatrix<T> & getMod()    const {return m_mod;}
protected:
    index_t            m_rhsNum; // number of rhs columns
    gsDofMapper        m_map;    // mapper from patch local basis function Ids to global functions Ids
    gsSparseMatrix<T>  m_sys;    // out sys matrix
    gsMatrix<T>        m_rhs;    // out rhs matrix
    gsSparseMatrix<T>  m_mod;    // out rhsmod matrix: the linear map between eliminated dofs and rhs

    SysW  m_sysWriter;
    RhsW  m_rhsWriter;
};


} // namespace gismo
