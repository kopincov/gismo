/**
 * @file gsPhysicalEvaluator.h
 *
 * \brief The physical evaluator is an object the wraps a parametric evaluator and
 *        and automatically computes the transformed values on the physical domain.
 *
 *        At the moment only 2 physical evaluators are provided:
 *          1) grad conforming for which the space on the physical domain is defined by
 *             \f[ \phi(x) = \hat\phi(\G^{-1}(x) ) \f], where \f[\phi\f] is the
 *             mapped function and  \f[G\f] is the parametrization of the physical domain.
 *             This is the most common mapping in IGA.
 *          2) div conforming for which the (vector valued) space on the physical domain is defined by
 *             \f[\phi(x) =  detJ^{-1}(x) J(\G^{-1}(x))\hat\phi(\G^{-1}(x) ) \f], where \f[\phi\f] is the
 *             mapped function, \f[G\f] is the parametrization of the physical domain and \f[J\f] the
 *             Jacobian of the parametrization.
 *
 *   This file is part of the G+Smo library.
 *
 *   This Source Code Form is subject to the terms of the Mozilla Public
 *   License, v. 2.0. If a copy of the MPL was not distributed with this
 *   file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *   Author(s): A. Bressan
**/

#pragma once

#include <gsCore/gsTransformPure.h>
#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <typename T>
class gsTransformedFuncSet : public gsFunctionSet<T>
{
public:
    typedef memory::shared_ptr<gsFunctionSet<T> > spacePtr;
public:
    gsTransformedFuncSet()
        :activeShift(0)
        {}
    index_t activeShift;
    virtual unsigned                 getGeoFlags (unsigned flags) const=0;
};





template <typename T, typename transform=gsTransformPureGradConforming<T> >
class gsTransformedFuncSetImpl : public gsTransformedFuncSet<T>
{
protected:
    // TODO move to shared pointers
    const gsFunctionSet<T>  &m_funValue;
public:

    gsTransformedFuncSetImpl (const gsFunctionSet<T> &func)
        : m_funValue(func)
    {
        gsTransformedFuncSet<T>::activeShift=0;
    }

public:
    virtual unsigned                 getGeoFlags       (unsigned flags) const {return transform::getGeoFlags(flags);}
    virtual const gsFunctionSet<T>&  getBaseFunctionSet () const {return m_funValue;}


    virtual index_t size () const {return m_funValue.size();}

    virtual short_t domainDim () const {return m_funValue.domainDim();}
    virtual short_t targetDim () const {return m_funValue.targetDim();}

    virtual void compute        (const gsMatrix <T> &points, gsFuncData<T> &result) const
    {
        GISMO_ERROR("transformed function sets need gsMapData input");
    }
    virtual void compute        (const gsMapData<T> &geoData, gsFuncData<T> &result) const
    {
        transform::transform(geoData,m_funValue,result);
        result.actives.array()+=gsTransformedFuncSet<T>::activeShift;
    }
    
    std::ostream &print(std::ostream &os) const { os << "gsTransformedFuncSetImpl\n"; return os; }
};


typedef gsTransformedFuncSetImpl<real_t, gsTransformPureGradConforming<real_t> > gsGradConformingTFS;
typedef gsTransformedFuncSetImpl<real_t, gsTransformPureDivConforming<real_t> >  gsDivConformingTFS;
typedef gsTransformedFuncSetImpl<real_t, gsTransformPureRestriction<real_t>  >   gsRestrictTFS;
typedef gsTransformedFuncSetImpl<real_t, gsTransformPureIdentity<real_t> >       gsIdentityTFS;


template <typename T=real_t>
struct gsFunctionSetDiff :public gsTransformedFuncSet<T>
{
    const gsFunctionSet<T>&m_a;
    const gsFunctionSet<T>&m_b;
public:
    gsFunctionSetDiff (   const gsFunctionSet<T>&a,  const gsFunctionSet<T>&b)
        : m_a(a), m_b(b)
    {}
    index_t     size() const {return m_a.size();}
    short_t   targetDim() const {return m_a.targetDim();}
    short_t   domainDim() const {return m_a.domainDim();}

    void difference (  gsFuncData<T> &r, const gsFuncData<T> &t) const
    {
        for (int d=0; d<r.maxDeriv()+1; ++d)
            r.values[d]-=t.values[d];
        if (t.flags & NEED_DIV)
            r.divs-=t.divs;
        if (t.flags & NEED_CURL)
            r.curls-=t.curls;
        if (t.flags & NEED_LAPLACIAN)
            r.laplacians-=t.laplacians;
    }

    unsigned getGeoFlags (unsigned flags) const
    {
        unsigned result=0;
        const gsTransformedFuncSet<T> *at=dynamic_cast<const gsTransformedFuncSet<T> *>(&m_a);
        if (at!=NULL)
            result = at->getGeoFlags(flags);
        at = dynamic_cast<const gsTransformedFuncSet<T> *>(&m_b);
        if (at!=NULL)
            result |= at->getGeoFlags(flags);
        return result;
    }

    void compute (const gsMapData<T> &p, gsFuncData<T> &out) const
    {
        gsFuncData<T> tmp(out.flags);
        tmp.patchId=out.patchId;
        m_a.compute(p,out);
        m_b.compute(p,tmp);
        difference(out,tmp);
    }
    void compute (const gsMatrix<T> &p, gsFuncData<T> &out) const
    {
        gsFuncData<T> tmp(out.flags);
        tmp.patchId=out.patchId;
        m_a.compute(p,out);
        m_b.compute(p,tmp);
        difference(out,tmp);
    }

    std::ostream &print(std::ostream &os) const { os << "gsFunctionSetDiff\n"; return os; }
};


} // namespace
