/** @file gsRemappedBasis.hpp

    @brief Implementation of gsRemappedBasis.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
**/

#include <gsRemappedBasis/gsRemappedBasis.h>
#include <gsRemappedBasis/gsVectorUtils.h>
#include <gsMapUtils/gsMapperUtils.h>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsMapUtils/gsMapFactoryMultiBasis.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

// Constructors

gsRemappedBasis::gsRemappedBasis( const gsWeightMapper<real_t>              &repr,
                                  const gsSelector                          &sele,
                                  const std::vector<gsFunctionSet<real_t>* >&basis ):
    m_repr ( repr ),
    m_sele ( sele ),
    m_isBasis(true)
{
    for (size_t b=0;b<basis.size();++b)
        m_basis.push_back(gsFunctionSet<real_t>::Ptr(basis[b]->clone()));

    init();
}

gsRemappedBasis::gsRemappedBasis( const gsWeightMapper<real_t>  &repr,
                                  const gsSelector              &sele,
                                  const std::vector<basisPtr>   &basis ):
    m_repr ( repr ),
    m_sele ( sele ),
    m_basis( basis),
    m_isBasis(true)
{
    init();
}

gsRemappedBasis::gsRemappedBasis(const gsRemappedBasis &other)
    : m_repr(other.m_repr),
      m_sele(other.m_sele),
      m_basis(other.m_basis),
      m_shift(other.m_shift),
      m_info(other.m_info),
      m_isBasis(other.m_isBasis)
{}

// protected
gsRemappedBasis::gsRemappedBasis()
    : m_isBasis(true)
{}

void gsRemappedBasis::init()
{
    GISMO_ASSERT (m_basis.size(), "cannot construct a gsRemappedBasis without a local basis\n");
    m_info = m_basis[0]->dimensions();

    checkDimAndInitShifts();
    m_repr.optimize();
}

// gsFunctionSet interface

index_t gsRemappedBasis::size() const
{
    return m_repr.getNrOfTargets();
}

short_t gsRemappedBasis::domainDim () const
{
    return m_info.first;
}

short_t gsRemappedBasis::targetDim () const
{
    return m_info.second;
}

void gsRemappedBasis::compute(const gsMapData<real_t> &geo,    gsFuncData<real_t> & result) const
{
    compute(geo.points,result);
}

void gsRemappedBasis::eval_into(const gsMatrix<real_t> &points, gsMatrix<real_t> &result) const
{
    gsFuncData<real_t> output(NEED_VALUE);
    compute(points,output);
    std::swap(result,output.values[0]);
}

void gsRemappedBasis::compute(const gsMatrix<real_t>   &points,
                              gsFuncData<real_t> &result) const
{
    const unsigned   flags=result.flags;
    const patchIdT   patch=result.patchId;

    result.dim = m_info;

    if (flags & SAME_ELEMENT || points.cols()==1 )
    {    // all points are in the same element

        result.values.resize(result.maxDeriv()+1);
        const size_t basisId=m_sele.getBasisAt(points.col(0),patch);

        gsFuncData<real_t> locBasisEval(result.flags | NEED_ACTIVE);
        m_basis[basisId]->compute(points,locBasisEval);

        gsWeightMapper<real_t>::IndexContainer target;
        gsMatrix<real_t>   Rm;
        getActiveAndMatrix(locBasisEval.actives,basisId,target,Rm);

        if (flags & NEED_ACTIVE )
            result.actives=gsAsConstVector<index_t>(target);
        if (flags & NEED_VALUE )
            remapValues(Rm,locBasisEval.values[0], result.values[0]);
        if (flags & NEED_DERIV )
            remapValues(Rm,locBasisEval.values[1], result.values[1]);
        if (flags & NEED_DIV)
            remapValues(Rm,locBasisEval.divs, result.divs);
        if (flags & NEED_CURL)
            remapValues(Rm,locBasisEval.curls, result.curls);
        if (flags & NEED_DERIV2)
            remapValues(Rm,locBasisEval.values[2], result.values[2]);
        if (flags & NEED_LAPLACIAN)
            remapValues(Rm,locBasisEval.laplacians, result.laplacians);
    }
    else
    {   // points can be in different bezier elements,
        // compute in each point and then join the results
        result.values.resize(result.maxDeriv()+1);
        const index_t numPoints=points.cols();
        // this implementation is faster because it operates on each point once,
        // but it requires a copy of the data and thus it uses double the memory
        // probably ok as for 100 000 points we talk of about 50 MB for values and first derivatives 2D

        // if memory is a problem the memory requirement can be halved by computing only the
        // the actives in a first loop over the points, then resizing the destination matrices
        // and lastly re-loop over the points to compute the values and copy them to the
        // destination matrices

        index_t maxActive=0;
        index_t maxDiv=0;
        index_t maxCurl=0;
        index_t maxLap=0;
        std::vector<index_t> maxValues(result.maxDeriv()+1,0);

        // this other version makes a copy of all data and thus causes
        // higher memory usage
        std::vector<gsFuncData<real_t> > resultAtPoint(numPoints);

        gsWeightMapper<real_t>::IndexContainer source;
        gsWeightMapper<real_t>::IndexContainer target;
        gsMatrix<real_t>   Rm;

        for (index_t p=0; p<numPoints; ++p)
        {
            resultAtPoint[p].values.resize(result.maxDeriv()+1);

            const size_t basisId=m_sele.getBasisAt(points.col(p),patch);
            gsFuncData<real_t> locBasisEval(result.flags | NEED_ACTIVE | SAME_ELEMENT);
            m_basis[basisId]->compute(points.col(p),locBasisEval);

            gsWeightMapper<real_t>::IndexContainer target2;
            gsMatrix<real_t>   Rm2;
            getActiveAndMatrix(locBasisEval.actives,basisId,target2,Rm2);

            maxActive=math::max<index_t>(maxActive,target2.size() );

            if (flags & NEED_ACTIVE)
            {
                resultAtPoint[p].actives=gsAsConstVector<index_t>(target2);
            }
            if (flags & NEED_VALUE )
            {
                remapValues(Rm2,locBasisEval.values[0], resultAtPoint[p].values[0]);
                maxValues[0]=math::max(maxValues[0],resultAtPoint[p].values[0].rows());
            }
            if (flags & NEED_DERIV )
            {
                remapValues(Rm2,locBasisEval.values[1], resultAtPoint[p].values[1]);
                maxValues[1]=math::max(maxValues[1],resultAtPoint[p].values[1].rows());
            }
            if (flags & NEED_DIV)
            {
                remapValues(Rm2,locBasisEval.divs, resultAtPoint[p].divs);
                maxDiv=math::max(maxDiv,resultAtPoint[p].divs.rows());
            }
            if (flags & NEED_CURL)
            {
                remapValues(Rm2,locBasisEval.curls, resultAtPoint[p].curls);
                maxCurl=math::max(maxCurl,resultAtPoint[p].curls.rows());
            }
            if (flags & NEED_DERIV2)
            {
                remapValues(Rm2,locBasisEval.values[2], resultAtPoint[p].values[2]);
                maxValues[2]=math::max(maxValues[2],resultAtPoint[p].values[2].rows());
            }
            if (flags & NEED_LAPLACIAN)
            {
                remapValues(Rm2,locBasisEval.laplacians, resultAtPoint[p].laplacians);
                maxLap=math::max(maxLap,resultAtPoint[p].laplacians.rows());
            }
        }
        // resize destination matrices
        if ( flags & NEED_ACTIVE )
            result.actives.resize(maxActive,numPoints);
        if (flags & NEED_VALUE)
            result.values[0].resize(maxValues[0],numPoints);
        if (flags & NEED_DERIV)
            result.values[1].resize(maxValues[1],numPoints);
        if (flags & NEED_DIV)
            result.divs.resize(maxDiv,numPoints);
        if (flags & NEED_CURL)
            result.curls.resize(maxCurl,numPoints);
        if (flags & NEED_DERIV2)
            result.values[2].resize(maxValues[2],numPoints);
        if (flags & NEED_LAPLACIAN)
            result.laplacians.resize(maxLap,numPoints);

        // copy the values
        for (index_t p=0; p<numPoints && flags&NEED_ACTIVE; ++p)
            copyCompleteCol(p, result.actives, resultAtPoint[p].actives );
        for (index_t p=0; p<numPoints && flags&NEED_VALUE; ++p)
            copyCompleteCol(p, result.values[0], resultAtPoint[p].values[0]);
        for (index_t p=0; p<numPoints && flags&NEED_DERIV; ++p)
            copyCompleteCol(p, result.values[1], resultAtPoint[p].values[1]);
        for (index_t p=0; p<numPoints && flags&NEED_DERIV2; ++p)
            copyCompleteCol(p, result.values[2], resultAtPoint[p].values[2]);
        for (index_t p=0; p<numPoints && flags&NEED_DIV; ++p)
            copyCompleteCol(p, result.divs, resultAtPoint[p].divs);
        for (index_t p=0; p<numPoints && flags&NEED_CURL; ++p)
            copyCompleteCol(p, result.curls, resultAtPoint[p].curls);
        for (index_t p=0; p<numPoints && flags&NEED_LAPLACIAN; ++p)
            copyCompleteCol(p, result.laplacians, resultAtPoint[p].laplacians);
        return;
    }
}

void gsRemappedBasis::getActiveAndMatrix(gsMatrix<index_t> &localAct, patchIdT basisId, gsWeightMapper<real_t>::IndexContainer &target, gsMatrix<real_t> &Rm) const
{
    localAct.array()+=m_shift[basisId];
    const index_t* dat= reinterpret_cast<index_t*>(localAct.data());
    const gsWeightMapper<real_t>::IndexContainer source(dat, dat+localAct.rows());

    if (m_isBasis)
    {
        m_repr.fastSourceToTarget(source,target);
        m_repr.getLocalMap(source,target,Rm);
    }
    else
    {
        const index_t numFunc = m_repr.getNrOfTargets()/m_info.second;
        target.resize(numFunc);
        for (index_t i=0;i<numFunc;++i)
            target[i]=i;
        m_repr.getLocalMap(source,Rm);
    }
}

template <typename T>
void gsRemappedBasis::copyCompleteCol (index_t c, gsMatrix<T> & dest, const gsMatrix<T> &data) const
{
    const index_t diff=dest.rows()-data.rows();
    dest.col(c).topRows(data.rows())=data.col(0);
    if (diff)
        dest.col(c).bottomRows(diff).setZero();
}


void gsRemappedBasis::remapValues (const gsMatrix<real_t>  &matrix,
                                   const gsMatrix<real_t> &input,
                                   gsMatrix<real_t> &output
                                   ) const
{
    const index_t numDelAct=matrix.rows();
    const index_t numOwnAct=matrix.cols();
    const index_t numRows=input.rows()/numDelAct;

    output.resize(numOwnAct*numRows,input.cols());
    for (index_t c=0;c<input.cols();++c)
        gsAsMatrix<real_t>(output.data()+c*output.rows(),numRows,numOwnAct)= gsAsConstMatrix<real_t>(input.data()+c*input.rows(),numRows,numDelAct)*matrix;
}

// Utils

const gsSelector& gsRemappedBasis::getSelector () const
{
    return m_sele;
}

void gsRemappedBasis::exportSelectorToTex(std::string filename) const
{
    m_sele.exportToTex(filename);
}

const gsWeightMapper<real_t>& gsRemappedBasis::getMapper () const
{
    return m_repr;
}

const std::vector<gsRemappedBasis::basisPtr> &gsRemappedBasis::getBases() const
{
    return m_basis;
}


void gsRemappedBasis::checkDimAndInitShifts()
{
    m_shift.resize(m_basis.size()+1);
    m_shift[0]=0;
    for (size_t b=0; b<m_basis.size();++b)
    {
        m_shift[b+1]=m_shift[b]+m_basis[b]->size();
    }
}

void gsRemappedBasis::reduce (const gsDofMapper &coeff)
{
    index_t maxRow=coeff.mapSize();
    index_t maxCol=0;

    gsSparseEntries<real_t> e;

    for(index_t i=0;i<maxRow;++i)
    {
        const index_t c=coeff.index(i);
        maxCol=math::max(maxCol,c);
        e.add(i,c,1);
    }
    gsSparseMatrix<real_t> mat(maxRow,maxCol);
    mat.setFromTriplets(e.begin(),e.end());
    reduce(mat);
}

std::vector<gsFunctionSet<real_t>::Ptr> cloneToShared (const gsMultiBasis<real_t> &input)
{
    std::vector<gsFunctionSet<real_t>::Ptr> result;
    for (size_t p=0;p<input.nBases();++p)
    {
        result.push_back(gsFunctionSet<real_t>::Ptr(input[p].clone()));
    }
    return result;
}

gsRemappedBasis* gsRemappedBasis::makeMultiPatch(const gsMultiBasis<real_t> &input, const gsWeightMapper<real_t> *gluing )
{
    return makeMultiPatch(cloneToShared(input),gluing);
}

gsRemappedBasis* gsRemappedBasis::makeMultiPatch(const gsMultiBasis<real_t> &input)
{
    gsWeightMapper<real_t> *mapper = gsMapFactoryMultiBasis(input).makeMapper();
    gsRemappedBasis *result = gsRemappedBasis::makeMultiPatch(cloneToShared(input), mapper);
    delete mapper;
    return result;
}

gsRemappedBasis* gsRemappedBasis::makeMultiPatch(const std::vector<basisPtr> input, const gsWeightMapper<real_t> *gluing )
{
    if (!input.size())
        return NULL;

    gsRemappedBasis* result = new gsRemappedBasis();
    result->m_info = input[0]->dimensions();

    // temporary data
    gsRemappedBasis* tmp;

    std::vector<gsWeightMapper<real_t>*> mapper;
    std::vector<bool>      toFree;

    gsMatrix<real_t> defaultBox(input[0]->domainDim(),2);
    for (index_t r=0;r<defaultBox.rows();++r)
        defaultBox.row(r)<<std::numeric_limits<real_t>::min(),std::numeric_limits<real_t>::max();

    index_t totPatch=0; // number of patches

    // build vector of bases,
    for (size_t b=0; b<input.size();++b)
    {
        tmp=dynamic_cast<gsRemappedBasis*>(input[b].get());
        if ( tmp == NULL )
        {
            gsBasis<real_t> *oldBasis = dynamic_cast<gsBasis<real_t>*>(input[b].get());
            if (oldBasis)
                result->m_sele.patch(totPatch).initConstant(static_cast<basisIdT>(result->m_basis.size()), oldBasis->support());
            else
                result->m_sele.patch(totPatch).initConstant(static_cast<basisIdT>(result->m_basis.size()), defaultBox);
            result->m_basis.push_back(input[b]);
            ++totPatch;
            // make idMatrix
            gsSparseMatrix<real_t> id(input[b]->size(),input[b]->size());
            for (index_t i=0;i<input[b]->size();++i)
                id.insert(i,i)=1;
            mapper.push_back(new gsWeightMapper<real_t>(id));
            toFree.push_back(true);
        }
        else
        {
            basisIdT basisIdShift=result->m_basis.size();
            result->m_basis.insert(result->m_basis.end(),tmp->m_basis.begin(),tmp->m_basis.end());
            for ( size_t p=0; p<tmp->m_sele.patchNum();++p)
            {
                result->m_sele.patch(totPatch)=shiftIds(tmp->m_sele.patch(p), basisIdShift);
                ++totPatch;
            }
            mapper.push_back(&tmp->m_repr);
            toFree.push_back(false);
        }
    }

    // check that all bases have the same domain dim
    try {
        result->checkDimAndInitShifts();
    }
    catch (...)
    {
        delete result;
        return NULL;
    }

    // construct the representation matrix
    // first step do the standard coomposition
    std::vector<index_t> colShifts;
    combineMappers(mapper,true,colShifts,result->m_repr);
    result->m_repr.optimize();

    // second step multiply by gluing
    if (gluing!=NULL)
    {
        gluing->optimize();
        result->m_repr*=gluing->asMatrix();
    }

    for (size_t m=0;m<mapper.size();++m)
    {
        if (toFree[m])
            delete mapper[m];
    }
    return result;
}



} // namespace gismo

