/** @file gsGeometryEvaluator.h

    @brief Provides implementation of BasisEvaluator class.
    
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
*/

#pragma once

#include <gsCore/gsBasisEvaluator.h>
#include <gsCore/gsGeoTransform.hpp>
#include <gsCore/gsGeometry.h>
#include <gsCore/gsDomainIterator.h>



namespace gismo {


template <typename T, short_t ParDim, short_t TarDim, int AmbDim>
gsBasisEvaluator<T> * makeBasisEvaluatorWithDims (const std::vector<gsBasis<T> *> &basis, const std::vector<index_t>& shifts, unsigned flags, ValueTransformationType geoTrans )
{
    const gsBasis<T> * basisPtr[TarDim];
    for (int i=0; i<TarDim;++i)
        basisPtr[i]=basis[i];
    switch (geoTrans)
    {
    case NO_TRANSFORMATION:
        return new gsGenericBasisEvaluator<T,ParDim,TarDim,gsGeoNoTransform<T,ParDim,TarDim,AmbDim> >(basisPtr,shifts,flags);
    case INVERSE_COMPOSITION:
        return new gsGenericBasisEvaluator<T,ParDim,TarDim,gsGeoGradPreservingTransform<T,ParDim,TarDim,AmbDim> >(basisPtr,shifts,flags);
    case DIV_CONFORMING:
        return new gsGenericBasisEvaluator<T,ParDim,TarDim,gsGeoDivPreservingTransform<T,ParDim,TarDim,AmbDim> >(basisPtr,shifts,flags);
    case CURL_CONFORMING:
        return new gsGenericBasisEvaluator<T,ParDim,TarDim,gsGeoCurlPreservingTransform<T,ParDim,TarDim,AmbDim> >(basisPtr,shifts,flags);
    default:
        gsWarn<<"I do not know this way to transform to the physical domain!\nUse constants from gsCore/gsBasisEvaluator.h.\n";
        return NULL;
    }
}


template <typename T, short_t ParDim, short_t TarDim>
gsBasisEvaluator<T> * makeBasisEvaluatorWithParAndTarDim (const std::vector<gsBasis<T> *> &basis, const std::vector<index_t>& shifts, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans )
{

    if (geo == NULL)
    {
        return makeBasisEvaluatorWithDims<T,ParDim,TarDim,ParDim>( basis, shifts, flags,  geoTrans);
    }
    else if( geo->coefDim() < ParDim)
    {
        gsWarn<<"Impossible to make basis evaluator with a parametrization from R^n to R^m with m<n.\n";
        return NULL;
    }
    else switch (geo->coefDim())
    {
    case 1:
        return makeBasisEvaluatorWithDims<T,ParDim,TarDim,1>( basis, shifts, flags, geoTrans);
    case 2:
        return makeBasisEvaluatorWithDims<T,ParDim,TarDim,2>( basis, shifts, flags, geoTrans);
    case 3:
        return makeBasisEvaluatorWithDims<T,ParDim,TarDim,3>( basis, shifts, flags, geoTrans);
    case 4:
        return makeBasisEvaluatorWithDims<T,ParDim,TarDim,4>( basis, shifts, flags, geoTrans);
    default:
        gsWarn<<"Cannot construct basis evaluator for basis on R^n n>=5\nuse the appropriate instantiation of gsGenericBasisEvaluator!\n";
        return NULL;
    }
}


template <typename T, short_t TarDim>
gsBasisEvaluator<T> * makeBasisEvaluatorWithTarDim (const std::vector<gsBasis<T> *> &basis, const std::vector<index_t>& shifts, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans )
{
    switch (basis[0]->dim())
    {
    case 1:
        return makeBasisEvaluatorWithParAndTarDim<T,1,TarDim>( basis, shifts, flags, geo,  geoTrans);
    case 2:
        return makeBasisEvaluatorWithParAndTarDim<T,2,TarDim>( basis, shifts, flags, geo,  geoTrans);
    case 3:
        return makeBasisEvaluatorWithParAndTarDim<T,3,TarDim>( basis, shifts, flags, geo,  geoTrans);
    case 4:
        return makeBasisEvaluatorWithParAndTarDim<T,4,TarDim>( basis, shifts, flags, geo,  geoTrans);
    default:
        gsWarn<<"Cannot construct basis evaluator for basis on R^n n>=5\nuse the appropriate instantiation of gsGenericBasisEvaluator!\n";
        return NULL;
    }
}


template <typename T>
gsBasisEvaluator<T> * makeBasisEvaluator (const std::vector<gsBasis<T> *> &basis, const std::vector<index_t>& shifts, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans )
{
    switch (basis.size())
    {
    case 0:
        gsWarn<<"Cannot make evaluator without basis!\n";
        return NULL;
    case 1:
        return makeBasisEvaluatorWithTarDim<T,1>(basis,shifts,flags,geo,geoTrans);
    case 2:
        return makeBasisEvaluatorWithTarDim<T,2>(basis,shifts,flags,geo,geoTrans);
    case 3:
        return makeBasisEvaluatorWithTarDim<T,3>(basis,shifts,flags,geo,geoTrans);
    case 4:
        return makeBasisEvaluatorWithTarDim<T,4>(basis,shifts,flags,geo,geoTrans);
    default:
        gsWarn<<"Cannot make evaluator with target dimension >=4.\nUse the appropriate instantiation of gsGenericEvaluator.\n";
        return NULL;
    }
}

template <typename T>
gsBasisEvaluator<T> * makeBasisEvaluator (const std::vector<gsBasis<T> *> &basis, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans )
{
    index_t accu=0;
    std::vector<index_t> shifts;
    for (size_t b=0; b<basis.size();++b)
    {
        shifts.push_back(accu);
        accu+=basis[b]->size();
    }
    return makeBasisEvaluator(basis,shifts,flags,geo,geoTrans);
}


template <typename T>
gsBasisEvaluator<T> * makeBasisEvaluator ( const gsBasis<T> &basis, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans )
{
    return makeBasisEvaluator(basis,0,flags,geo,geoTrans);
}

template <typename T>
gsBasisEvaluator<T> * makeBasisEvaluator ( const gsBasis<T> &basis, index_t  shift, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans)
{
    std::vector<gsBasis<T>* >  basis_Array(1,const_cast<gsBasis<T>*>(&basis));
    std::vector<index_t>             shifts(1,shift);
    return makeBasisEvaluator(basis_Array,shifts,flags,geo,geoTrans);
}



template <typename T, short_t ParDim, short_t TarDim, typename geometryTransform >
gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::gsGenericBasisEvaluator( const gsBasis<T> *(&basis)[TarDim], unsigned flags)
    : gsBasisEvaluator<T>(flags)
{
    m_active_shift[0]=0;
    for (int i =1;i<TarDim;++i)
        m_active_shift[i]=m_active_shift[i-1]+basis[i-1]->size();
    init(basis,flags);
}

template <typename T, short_t ParDim, short_t TarDim, typename geometryTransform >
gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::gsGenericBasisEvaluator( const gsBasis<T> *(&basis)[TarDim], const std::vector<index_t>& shifts, unsigned flags)
    : gsBasisEvaluator<T>(flags)
{
    GISMO_ASSERT(TarDim==shifts.size(),"Target dimension does not fit shifts size.");
    for (int i =0;i<TarDim;++i)
        m_active_shift[i]=shifts[i];
    init(basis,flags);
}

template <typename T, short_t ParDim, short_t TarDim, typename geometryTransform >
void gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::init( const gsBasis<T> *(&basis)[TarDim], unsigned flags)
{
    for (int i=0;i<TarDim;++i)
        m_basis[i]=basis[i];
    setFlags(flags);
    m_parDim=ParDim;
    m_tarDim=TarDim;
    m_spaceDim=m_active_shift[TarDim-1]+m_basis[TarDim-1]->size();
}

template <typename T, short_t ParDim, short_t TarDim, typename geometryTransform >
void gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::setFlags   ( unsigned newFlags)
{
    m_flags=geometryTransform::addAuxiliaryFlags(newFlags);
    m_geo_flags=geometryTransform::getGeometryFlags(m_flags);

    if (m_flags & NEED_2ND_DER)
        m_max_deriv=2;
    else if (m_flags & (NEED_GRAD | NEED_JACOBIAN) )
        m_max_deriv=1;
    else if (m_flags & NEED_VALUE)
        m_max_deriv=0;
    else
        m_max_deriv=-1;
}


template <typename T, short_t ParDim, short_t TarDim, typename geometryTransform >
void gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::addActiveShift (unsigned shift)
{
    for (int i =TarDim-1;i>=0;--i)
    {
        m_active_shift[i]+=shift;
    }
}

template <typename T, short_t ParDim, short_t TarDim, typename geometryTransform >
gsDomainIterator<T>* gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::getDomainIterator(const boxSide s) const
{
    return m_basis[0]->makeDomainIterator(s).release();
}

template <typename T, short_t ParDim, short_t TarDim, typename geometryTransform >
gsVector<unsigned>  gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::getDegree() const
{
    gsVector<unsigned> result;
    result.setZero(ParDim);
    for (index_t i=0;i<TarDim;++i)
    {
        for (index_t d=0;d<ParDim;++d)
        {
            result(d)=std::max<unsigned>(result(d),m_basis[i]->degree(d));
        }
    }
    return result;
}


template <typename T, short_t ParDim, short_t TarDim, typename geometryTransform >
void gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::evaluateAt ( const gsMatrix<T> &points)
{
    int                 active_num[TarDim];
    gsMatrix<index_t>   active[TarDim];
    unsigned            tot_active(0);

    std::vector<gsMatrix<T> > tmp;
    index_t c;

    for (short_t i =0;i<TarDim;++i)
    {
        m_basis[i]->active_into(points.col(0),active[i]);
        active_num[i]=active[i].rows();
        tot_active+=active_num[i];
        m_basis[i]->evalAllDers_into(points, m_max_deriv, tmp);//basis values

        // copy back to one matrix
        c = 0;
        for (int k =0; k<=m_max_deriv; ++k)
            c += tmp[k].rows();
        m_basis_vals[i].resize(c, points.cols());

        c = 0;
        for (int k =0; k<=m_max_deriv; ++k)
        {
            m_basis_vals[i].middleRows(c,tmp[k].rows() ).swap( tmp[k] );
            c+= tmp[k].rows();
        }

    }
    m_actives.resize(tot_active,1);

    int size;
    int start = tot_active;
    for (int i =TarDim-1;i>=0;--i)
    {
        size  = active[i].rows();
        start -=size;
        m_actives.block(start,0,size,1)=(active[i].array()+m_active_shift[i]);
    }

    if (this->m_flags & NEED_VALUE)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeValues(this, NULL, m_basis_vals, active_num, m_values);
    if (this->m_flags & NEED_GRAD)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeGrads(this, NULL, m_basis_vals, active_num, m_derivs);
    if (this->m_flags & NEED_JACOBIAN)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeJacobians(this, NULL, m_basis_vals, active_num, m_jacobians);
    if (this->m_flags & NEED_DIV)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeDivs(this,  NULL, m_basis_vals, active_num,m_divs);
    if (this->m_flags & NEED_CURL)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeCurls(this, NULL, m_basis_vals, active_num, m_curls);
    if (this->m_flags & NEED_2ND_DER)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeSecDers(this, NULL, m_basis_vals, active_num, m_2ndDers);
    if (this->m_flags & NEED_LAPLACIAN)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeLaplacians(this, NULL, m_basis_vals, active_num, m_laps);
}


template <typename T, short_t ParDim, short_t TarDim, typename geometryTransform >
void gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::evaluateAt ( const gsMatrix<T> &points, const gsGeometryEvaluator<T> &geoEval)
{
    int                 active_num[TarDim];
    gsMatrix<index_t>   active[TarDim];
    unsigned            tot_active(0);

    for (int i =0;i<TarDim;++i)
    {
        m_basis[i]->active_into(points.col(0),active[i]);
        active_num[i]=active[i].rows();
        tot_active+=active_num[i];
    }
    m_actives.resize(tot_active,1);

    // stack together actives together with offsets
    int size;
    int start = tot_active;
    for (int i =TarDim-1;i>=0;--i)
    {
        size  = active[i].rows();
        start -=size;
        m_actives.block(start,0,size,1)=(active[i].array()+m_active_shift[i]);
    }

    std::vector<gsMatrix<T> > tmp;
    index_t c;
    for (int i =0;i<TarDim && m_max_deriv>=0;++i)
    {
        m_basis[i]->evalAllDers_into(points, m_max_deriv, tmp);

        // copy back to one matrix
        c = 0;
        for (int k =0; k<=m_max_deriv; ++k)
            c += tmp[k].rows();
        m_basis_vals[i].resize(c, points.cols());

        c = 0;
        for (int k =0; k<=m_max_deriv; ++k)
        {
            m_basis_vals[i].middleRows(c,tmp[k].rows() ).swap( tmp[k] );
            c+= tmp[k].rows();
        }
    }


    if (this->m_flags & NEED_VALUE)
        geometryTransform::computeValues(this, &geoEval, m_basis_vals, active_num, m_values);
    if (this->m_flags & NEED_GRAD)
        geometryTransform::computeGrads(this, &geoEval, m_basis_vals, active_num, m_derivs);
    if (this->m_flags & NEED_JACOBIAN)
        geometryTransform::computeJacobians(this, &geoEval, m_basis_vals, active_num, m_jacobians);
    if (this->m_flags & NEED_DIV)
        geometryTransform::computeDivs(this,  &geoEval, m_basis_vals, active_num, m_divs);
    if (this->m_flags & NEED_CURL)
        geometryTransform::computeCurls(this, &geoEval, m_basis_vals, active_num,m_curls);
    if (this->m_flags & NEED_2ND_DER)
        geometryTransform::computeSecDers(this, &geoEval, m_basis_vals, active_num,m_2ndDers);
    if (this->m_flags & NEED_LAPLACIAN)
        geometryTransform::computeLaplacians(this, &geoEval, m_basis_vals, active_num,m_laps);
}


template <typename T, short_t ParDim, short_t TarDim, typename geometryTransform >
std::vector<index_t> gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::getBoundary(boxSide bs,index_t offset) const
{
    gsMatrix<index_t> piece;
    std::vector<index_t> r;
    for (int i=0; i<TarDim; ++i)
    {
        piece=m_basis[i]->boundaryOffset(bs,offset);
        const index_t size=piece.rows();
        for(index_t j = 0;j<size;++j)
            r.push_back(piece.at(j)+m_active_shift[i]);
    }
    return r;
}


} // namespace gismo


