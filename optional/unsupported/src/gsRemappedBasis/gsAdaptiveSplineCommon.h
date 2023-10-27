/** @file gsAdaptiveSplineCommon.h

    @brief
    Common code shared by the implementations of THB, HLR, DHB, and NTHB
    as remapped basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
*/


#pragma once

#include <gsRemappedBasis/gsTensorBasesUtils.h>
#include <gsRemappedBasis/gsDomainMapApply.hpp>
#include <gsRemappedBasis/gsRemappedBasis.h>
#include <gsRemappedBasis/gsDomainMapIterator.h>

namespace gismo
{

struct gsBoxRefinableInterface
{
    virtual void refine (const gsBoxList & boxes ) =0;

    virtual void refineWithTransfer (const gsBoxList & boxes,  gsSparseMatrix<real_t> &transfer ) =0;
    virtual void refineWithCoefs    (const gsBoxList & boxes,  gsMatrix<real_t> &coefs )
    {
        gsSparseMatrix<real_t> transfer;
        refineWithTransfer(boxes,transfer);
        coefs=transfer*coefs;
    }
};


/** @brief A class implementing the common features of the spline
    bases derived from the gsRemappedBasis:
    - gsTHB,
    - gsHLR,
    - gsTBPN,
    - gsDecoupledBasis.
 */
template <short_t DIM>
class gsAdaptiveSplineCommon : public gsRemappedBasis
{
public:
    using  gsRemappedBasis::m_basis;
    using  gsRemappedBasis::m_info;
    using  gsRemappedBasis::m_sele;
    using  gsRemappedBasis::m_repr;
    using  gsRemappedBasis::m_shift;

    typedef typename gsDomainMap::actionBase<DIM,2>::BoxT BoxT;
    typedef typename gsBSplineTraits<DIM,real_t>::Basis   tensorBasisT;
    typedef gsKnotVector<real_t>                          knotVectorT;
    
    typedef typename knotVectorT::const_iterator      knotIter;
    typedef typename std::reverse_iterator<knotIter>  backKnotIter;

    typedef std::vector<index_t>                       indexSetT;
    // types for multi-indices
    typedef typename gsTensorIndexDomain<DIM>::DomainT indexDomT;
    typedef typename gsTensorIndexDomain<DIM>::MIndexT indexVecT;
public:
    gsBasis<real_t>::domainIter makeDomainIterator(size_t p=0) const
    {
        return gsBasis<real_t>::domainIter(new gsDomainMapIterator<DIM>(m_sele.patch(p),m_basis));
    }


    /** Increases the size of the boxes so that they are aligned
     * with knots. It aligns the boxes according to their level
     * shifted by the optional parameter shift
     */
    void sanitizeBoxList (const std::vector<basisPtr>& basisArg, gsBoxList &boxes, int shift) const
    {
        for (size_t b=0; b<boxes.size();++b)
            sanitizeBox(basisArg,boxes.box(b),boxes.basisId(b)+shift);
    }

    /// See void sanitizeBoxList(std::vector<basisPtr>&, gsBoxList, int) const.
    void sanitizeBoxList (gsBoxList &boxes, int shift) const
    {
        sanitizeBoxList(m_basis,boxes,shift);
    }

    /// See sanitizeBoxList.
    template <typename boxT>
    void sanitizeBox (boxT box, int lvl) const
    {
        sanitizeBox<boxT>(m_basis,box,lvl);
    }

    /// See sanitizeBoxList.
    template <typename boxT>
    void sanitizeBox (const std::vector<basisPtr>& basis, boxT box, int lvl) const
    {
        typedef typename gsRemTypes<DIM>::knotVectorT::const_uiterator const_uiterator;
        typedef typename gsRemTypes<DIM>::tensorBasisT                 tensorBasisT;

        lvl = lvl>0? lvl: 0;
        tensorBasisT &_basis=*asBasis<DIM>(basis[lvl]);

        for (short_t d=0; d<m_info.first; ++d)
        {
            const_uiterator kL=_basis.knots(d).ubegin()+_basis.knots(d).uFind(box(d,0)).uIndex();
            const_uiterator kU=_basis.knots(d).ubegin()+_basis.knots(d).uFind(box(d,1)).uIndex();
            const_uiterator kM=_basis.knots(d).uend()-1;
            box(d,0)=*kL;
            if ( box(d,1)<= *kU )
                box(d,1)=*kU;
            else if (kU!=kM)
                box(d,1)=*(kU+1);
            else GISMO_ERROR("???");
        }
    }


    /** Structure for collecting the minimal and maximal level in a
        specified domain with the selector. For an
        example of usage, see gsHLR::representStride(...).
     */
    struct getMinAndMaxLevel : public gsDomainMap::actionBase<DIM,2>
    {
        typedef typename gsDomainMap::actionBase<DIM,2>::BoxT BoxT;

        basisIdT                                     fromLvl;
        basisIdT                                     minLvl;
        basisIdT                                     maxLvl;

        getMinAndMaxLevel(basisIdT from = 0)
            : fromLvl(from), minLvl(std::numeric_limits<basisIdT>::max()) , maxLvl(0)
        {}

        bool enter (const gsDomainMap*        map,
                    const BoxT              & /*box*/,
                    const gsDomainMap::NodeId node )
        {
            if ((*map)[node].isFork() )
                return true;
            const basisIdT curLvl=(*map)[node].data.space;

            minLvl=math::min(minLvl,curLvl);
            maxLvl=math::max(maxLvl,curLvl);
            if ( curLvl < fromLvl)
                return false;
            return true;
        }
    };


    /** Structure for collecting the influencing functions in a
     * specified area (i.e., those with a support intersecting the
     * area). For an example of usage, see
     * gsHLR::makeRepresentation(). */
    struct collectInfluencing : public gsDomainMap::actionBase<DIM,2>
    {
        typedef typename gsDomainMap::actionBase<DIM,2>::BoxT BoxT;

        std::vector<indexSetT>                        data;
        indexSetT                                     tmp;

        basisIdT                                     fromLvl;
        basisIdT                                     minLvl;
        basisIdT                                     maxLvl;

        const std::vector<basisPtr>                 *bases;
        std::vector<gsMatrix<index_t,DIM,2>, 
        typename gsMatrix<index_t,DIM,2>::aalloc>    limits;

        collectInfluencing(const std::vector<basisPtr> &b)
            :  fromLvl(0), bases(&b)
        {initLimit(); init();}
        collectInfluencing(const std::vector<basisPtr> &b, basisIdT  mLvl, const std::vector<indexDomT > &l)
            :  fromLvl(mLvl), bases(&b), limits(l)
        {init();}

        void initLimit()
        {
            limits.resize(bases->size());
            for (basisIdT b=0; b<bases->size(); ++b)
            {
                for (short_t dir=0;dir<DIM;++dir)
                {
                    const knotVectorT &knots = asBasis<DIM>((*bases)[b])->knots(dir);
                    const index_t      deg   = asBasis<DIM>((*bases)[b])->degree(dir);
                    limits[b](dir,0) = 0;
                    limits[b](dir,1) = knots.size() -deg-1;
                }
            }
        }

        void init()
        {
            data.resize(bases->size());
            for (size_t i=0; i<bases->size();++i)
                data[i].resize(0);
            minLvl=std::numeric_limits<basisIdT>::max();
            maxLvl=std::numeric_limits<basisIdT>::min();
        }

        bool enter (const gsDomainMap *map, const BoxT& box, const gsDomainMap::NodeId node )
        {
            if ((*map)[node].isFork())
                return true;

            const basisIdT curLvl = (*map)[node].data.space;

            minLvl=math::min(minLvl,curLvl);
            maxLvl=math::max(maxLvl,curLvl);
            if ( curLvl < fromLvl)
                return true;

            const tensorBasisT  &tensorBasis = *asBasis<DIM>((*bases)[curLvl]);
            gsMatrix<index_t, DIM,2> dom, subdom;

            for( short_t dir = 0; dir < DIM; ++dir )
            {
                const index_t deg=tensorBasis.degree(dir);
                const knotVectorT &knots=tensorBasis.knots(dir);
                dom.row(dir)<< 0,  knots.size() - deg - 1;
                const knotIter beg=knots.begin();
                subdom(dir,0)=std::upper_bound(beg+limits[curLvl](dir,0),beg+limits[curLvl](dir,1),box(dir,0)) - beg;
                subdom(dir,0)=std::max(limits[curLvl](dir,0), subdom(dir,0)- deg- 1);
                subdom(dir,1)=std::lower_bound(beg+limits[curLvl](dir,0),beg+limits[curLvl](dir,1),box(dir,1)) - beg;
            }
            gsTensorIndexSubDomain<DIM> MSD(gsTensorIndexDomain<DIM>(dom),subdom);
            tmp.resize(MSD.size()+data[curLvl].size());
            indexSetT::const_iterator last=std::set_union(MSD.flat_begin(),MSD.flat_end(),data[curLvl].begin(),data[curLvl].end(),tmp.begin());
            tmp.resize(last-tmp.begin());
            std::swap(tmp,data[curLvl]);
            return true;
        }
    };


    /** Collects the tensor indices of a domain for the given level
     * @param lvl.
     */
    gsTensorIndexDomain<DIM> getTensorIndexDomainForLevel(basisIdT lvl)
    {
        const tensorBasisT  &tensorBasis = *asBasis<DIM>(m_basis[lvl]);
        gsMatrix<index_t, DIM,2> dom;
        for( short_t dir = 0; dir < DIM; ++dir )
        {
            const index_t deg=tensorBasis.degree(dir);
            const knotVectorT &knots=tensorBasis.knots(dir);
            dom.row(dir) << 0, knots.size() - deg - 1;
        }
        return gsTensorIndexDomain<DIM>(dom);
    }


    /** A wrapper for a subset of a knot vector.
     */
    struct knotRange
    {
        typedef typename gsKnotVector<real_t>::const_iterator iter;
        typedef typename gsKnotVector<real_t>::const_reverse_iterator riter;

        iter beg;
        iter end;

        size_t            size() const  {return end-beg;}
        knotRange&        operator++()  {++beg; return *this;}
        bool              empty() const {return end==beg;}
        const real_t&     front() const {return *beg;}
        const real_t&     back()  const {return *(end-1);}
    };

    /**
       @brief The tensorCoefs struct
       Represent a tensor structure using a fixed memory layout determined
       by memSize
    **/
    struct tensorCoefs
    {
        typedef gsVector<index_t,DIM> Vector;
    private:
        void initStride()
        {
            stride(0)=1;
            for (short_t dir=1; dir<DIM; ++dir)
                stride(dir)=stride(dir-1)*memSize(dir-1);
        }

    public:
        Vector memSize;
        Vector stride;

        Vector curSize;
        std::vector<real_t>   data;

    public:
        real_t& operator[] (index_t i)       {return data[i];}
        const real_t& operator[] (index_t i) const {return data[i];}
        real_t& operator[] (const Vector &i)       {return data[i.dot(stride)];}
        const real_t& operator[] (const Vector &i) const {return data[i.dot(stride)];}

        void initOne (const Vector &size, real_t start=1)
        {
            memSize=size;
            curSize=Vector::Constant(1);
            initStride();
            data.resize(memSize.prod());
            if (data.size()>0)
                data[0]=start;
        }
    };

    /**
     * @brief updateCoefs Changes the coefficients so that they
     * represent the whole function on a finer knot vector.
     * @param coefs initial coefficients, will be rewritten. Warning,
     * the vector must be big enough to hold the representation in
     * the finest space.
     * @param stride tells after how many coefs to start a new line.
     * @param tmpSize
     * @param tmpStorage
     */
    void static updateCoefs(
            const knotRange            (&oldK)[DIM],
            const knotRange            (&newK)[DIM],
            tensorCoefs           &coefs,
            std::vector<real_t>         &tmpStorage)
    {
        index_t sizeSave;
        index_t deg;
        for (short_t dir=DIM-1; dir>0;--dir)
        {
            sizeSave=coefs.curSize(dir);

            deg = oldK[dir].size()-sizeSave-1;

            coefs.curSize(dir)=1;
            gsTensorIndexDomain<DIM> lines(coefs.curSize);
            coefs.curSize(dir)=newK[dir].size()-deg-1;
            typename gsTensorIndexDomain<DIM>::const_iterator it=lines.begin(), end=lines.end();

            for (;it!=end;++it)
            {
                index_t base=it->dot(coefs.stride);
                for (index_t i=0; i<sizeSave; ++i)
                    tmpStorage[i]=coefs[base+i*coefs.stride(dir)];
                updateCoefs1D(deg, oldK[dir], newK[dir], &tmpStorage[0]);
                for (index_t i=0; i<coefs.curSize(dir); ++i)
                    coefs[base+i*coefs.stride(dir)]=tmpStorage[i];
            }
        }
        sizeSave=coefs.curSize(0);
        deg = oldK[0].size()-sizeSave-1;

        coefs.curSize(0)=1;
        gsTensorIndexDomain<DIM> lines(coefs.curSize);
        coefs.curSize(0)=newK[0].size()-deg-1;
        typename gsTensorIndexDomain<DIM>::const_iterator it=lines.begin(), end=lines.end();

        for (;it!=end;++it)
            updateCoefs1D(deg, oldK[0], newK[0], &coefs[*it]); // continuous memory
    }


    /**
     * @brief updateCoefs1D inserts knots to bring coefficients
     *  relative to oldK to coefficients relative to newK
     * @param deg degree,
     * @param oldK old knots,
     * @param newK new knots,
     * @param coefs coefficients.
     */
    static void updateCoefs1D(index_t deg , const knotRange oldK, const knotRange newK, real_t *coefs)
    {
        typename knotRange::iter nK = newK.beg; //newKnot
        typename knotRange::iter oK = oldK.beg; //oldKnot

        while ( nK != newK.end && oK !=oldK.end )
        {
            if ( *oK > *nK )
            {
                insertKnot(deg, oldK, newK, oK, nK, coefs);
                ++nK;
            }
            else if ( *oK == *nK )
            {
                ++oK;
                ++nK;
            }
            else if ( *oK  < *nK )
                ++oK; // should never happen
        }
    }

    /** @brief Boehm algorithm in place.

        The memory can be allocated in order to increase coefs so
    better reserve in the first place. This implementation assumes
    that the knots from \a newKnot are being inserted in order into \a
    oldK; so that the current knot vector when inserting \a nK is \a
    [*newK.beg .... *(nK-1),*oK .... *(oldK.end-1)] */
    static inline void insertKnot(index_t deg , knotRange oldK, knotRange newK, typename knotRange::iter oK, typename knotRange::iter nK, real_t * cBeg)
    {
        const index_t mu = nK-newK.beg-1; // the new knot is between current knots mu and mu+1
        const index_t sM = std::min<index_t>(oldK.end-oK-1, deg-1); // max shift from oK
        const index_t rS = deg-1-sM;

        real_t * cEnd = cBeg+(mu+1)+(oldK.end-oK)-deg;
        const index_t cSiz = cEnd-cBeg;
        *(cBeg+cSiz-1)=0; // set the new coef to 0

        // shift coefficients from mu+1 to the end by 1
        real_t * movBeg = cBeg + std::min<ptrdiff_t>(mu+1, cSiz);
        real_t * movEnd = cEnd;
        if ( movEnd >movBeg)
#       ifdef _MSC_VER
            // avoid C4996 warning
            std::copy_backward(movBeg-1, movEnd-1,
                               stdext::unchecked_array_iterator<real_t*>(movEnd) );
#       else
            std::copy_backward(movBeg-1, movEnd-1, movEnd );
#       endif
        
        // coefficients i for i in [mu-deg+1, mu] are linear combinations of coefs i and i-1
        // so we loop backward from min(mu, coefs.size()-1) to max(1,mu-deg+1)
        // but if mu-deg+1 is <= 0 then we take special care of the first result
        // coefficients
        real_t * fstC = cBeg+std::max<ptrdiff_t>(0,mu-deg);
        real_t * lstC = movBeg-1;
        typename knotRange::iter gK(oK+sM);   // knots >  *nK
        typename knotRange::iter sK(nK-1-rS); // knots <= *nK
        while ( fstC != lstC )
        {
            *lstC = ( (*nK-*sK)*(*lstC) + (*gK-*nK)*( *(lstC-1) ) ) / (*gK-*sK);
            --lstC;
            --gK;
            --sK;
        }
        // is mu-deg < 0 then we should compute coef[0] with one term only
        if (mu<deg)
            *lstC *= (*nK-*sK) / (*gK-*sK);
    }

    inline void getTensorFunctionInfo(const index_t id, const basisIdT lvl, BoxT &support, indexDomT &mult, knotRange knot[DIM])
    {
        tensorBasisT               &_basis = *asBasis<DIM>(m_basis[lvl]);
        gsVector<index_t, DIM>     idM = _basis.tensorIndex(id).template cast<index_t>();
        getTensorFunctionInfo(_basis,idM,support,mult,knot);
    }

    inline void getTensorFunctionInfo(const basisPtr &basis, const index_t id, BoxT &support, indexDomT &mult, knotRange knot[DIM])
    {
        tensorBasisT            &_basis = *asBasis<DIM>(basis);
        gsVector<index_t, DIM>     idM = _basis.tensorIndex(id).template cast<index_t>();
        getTensorFunctionInfo(_basis, idM, support, mult, knot);
    }

    static inline void getTensorFunctionInfo(const tensorBasisT &_basis, const indexVecT &idM, BoxT &support, indexDomT &mult, knotRange knot[DIM])
    {
        for (short_t dir=0;dir<DIM;++dir)
        {
            getFuncInfo1D(_basis.component(dir),idM(dir),mult.row(dir),knot[dir]);
            support(dir,0)=knot[dir].front();
            support(dir,1)=knot[dir].back();

            // get function support and knots
            const int deg  =  _basis.degree(dir);
            knotIter  beg  =  _basis.knots(dir).begin()+idM(dir);
            knotIter  end  =  beg+deg+2;
            support(dir,0) = *beg;
            support(dir,1) = *(end-1);
            mult(dir,0)    = std::find_if(beg,end,std::bind1st(std::not_equal_to<real_t>(),support(dir,0)))-beg;
            mult(dir,1)    = std::find_if(backKnotIter(end),backKnotIter(beg),std::bind1st(std::not_equal_to<real_t>(),support(dir,1)))-backKnotIter(end);
            knot[dir].beg  = beg;
            knot[dir].end  = end;
        }
    }

    static inline void getFuncInfo1D(
            const typename gsRemTypes<DIM>::tensorBasisT::Basis_t    & _basis,
            const index_t                 &id,
            typename indexDomT::Row        mult,
            knotRange                     &knots
            )
    {
        const int deg  = _basis.degree();
        knots.beg      = _basis.knots ().begin()+id;
        knots.end      = knots.beg+deg+2;
        mult(0,0)      = std::find_if(knots.beg,knots.end,std::bind1st(std::not_equal_to<real_t>(),knots.front()))-knots.beg;
        mult(0,1)      = std::find_if(backKnotIter(knots.end),backKnotIter(knots.beg),std::bind1st(std::not_equal_to<real_t>(),knots.back()))-backKnotIter(knots.end);
    }


    inline void getTensorFunctionInfo(const knotRange k[DIM], BoxT &support, indexDomT &mult)
    {
        for (short_t dir=0;dir<DIM;++dir)
        {
            support(dir,0) = k[dir].front();
            support(dir,1) = k[dir].back();
            mult(dir,0)    = std::find_if(k[dir].beg,k[dir].end,std::bind1st(std::not_equal_to<real_t>(),support(dir,0)))-k[dir].beg;
            mult(dir,1)    = std::find_if(backKnotIter(k[dir].end),backKnotIter(k[dir].beg),std::bind1st(std::not_equal_to<real_t>(),support(dir,1)))-backKnotIter(k[dir].end);
        }
    }

    inline void getTensorLevelLocalView(const basisIdT lvl, const BoxT &area, const indexDomT &mult, knotRange knot[DIM], gsTensorIndexSubDomain<DIM>  &subDom)
    {
        tensorBasisT & _basis=*asBasis<DIM>(m_basis[lvl]);
        getTensorLevelLocalView(_basis,area,mult,knot,subDom);
    }
    static inline void getTensorLevelLocalView(basisPtr basis, const BoxT &area, const indexDomT &mult, knotRange knot[DIM], gsTensorIndexSubDomain<DIM>  &subDom)
    {
        tensorBasisT & _basis=*asBasis<DIM>(basis);
        getTensorLevelLocalView(_basis,area,mult,knot,subDom);
    }

    static inline void getTensorLevelLocalView(const tensorBasisT & _basis, const BoxT &area, const indexDomT &mult, knotRange knot[DIM], gsTensorIndexSubDomain<DIM>  &subDom)
    {
        indexDomT limit;
        indexVecT size;
        for (short_t dir=0;dir<DIM;++dir)
        {
            const knotVectorT &allKnots = _basis.knots(dir);
            const index_t      deg   = _basis.degree(dir);
            knotIter beg  = allKnots.begin();
            knotIter end  = allKnots.end();
            knot[dir].beg = std::upper_bound(beg,end,area(dir,0))-mult(dir,0);
            knot[dir].end = std::lower_bound(knot[dir].beg,end,area(dir,1))+mult(dir,1);

            limit(dir,0) = knot[dir].beg-beg;
            limit(dir,1) = knot[dir].end-beg-deg-1;
            size(dir)    = allKnots.size()-deg-1;
        }
        subDom=gsTensorIndexSubDomain<DIM>(gsTensorIndexDomain<DIM>(size),limit);
    }

    // Thomas Broke Something 42-times.
    /** Class for representing a tensor-product B-spline basis on a
     * subdomain. Corresponds to \cite bm2016, Section 4.1, paragraph
     * Implementation 2. It is used in gsDecoupledBasis.
     */
    class gsTBS42 : public gsFunctionSet<real_t>
    {
    public:
        std::ostream &print(std::ostream &os) const { return os;}

        gsTBS42( const typename gsRemTypes<DIM>::tensorBasisT &locBasis,
                 gsTensorIndexSubDomain<DIM> &subdomain,
                 basisIdT origLvl,
                 gsDomainMap leaf)
            : m_locBasis(locBasis), m_subdomain(subdomain), m_origLvl(origLvl), m_leaf(leaf)
        {}

    public: // Getters
        short_t domainDim() const {return DIM;}
        short_t targetDim() const {return 1;}
        index_t size() const {return m_subdomain.size();}
        basisIdT origLvl() const {return m_origLvl;}
        gsDomainMap leafDomain() const {return m_leaf;}
        const gsTensorIndexSubDomain<DIM> &subdomain() const {return m_subdomain;}

    public:
        /** Calls compute of the \a m_locBasis and writes the active
         * functions if required.
         */
        void compute(const gsMatrix<real_t> &pts, gsFuncData<real_t> &result) const
        {
            m_locBasis.compute(pts,result);
            if(result.flags & NEED_ACTIVE)
                transformActives(result.actives);
        }
        
    private:
        /** Transforms the indices of the active functions from the
         * globals (w.r.t. the whole basis) to the locals (w.r.t. the
         * basis on the subdomain). */
        void transformActives(gsMatrix<index_t>& actives) const
        {
            // TODO: check for correctness, changed from unsigned to index_t
            index_t *b = actives.data();
            index_t *e = b+actives.size();
            for( index_t *p=b; p < e; ++p)
                *p = m_subdomain.toPos(m_subdomain.toMulti(*p));
        }
    protected: // members
        const typename gsRemTypes<DIM>::tensorBasisT &m_locBasis;
        gsTensorIndexSubDomain<DIM> m_subdomain;
        basisIdT                    m_origLvl;
        gsDomainMap                 m_leaf;
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    gsTBS42& asTBS42(basisPtr ptr)
    {
        return *static_cast<gsTBS42*>(ptr.get());
    }

    /** Class for splitting the basis into leaves and thus making one
     * per leaf. Used in combination with the selector.
     */
    class makeBasisPerLeaf : public gsDomainMap::actionBase<DIM,2>
    {
    public:
        typedef typename gsDomainMap::actionBase<DIM,2>::BoxT BoxT;
    protected:
        const std::vector<basisPtr> &m_defBasis;
        std::vector<basisPtr>       &m_newBasis;
    public:
        makeBasisPerLeaf(const std::vector<basisPtr> &defBasis,
                         std::vector<basisPtr> &newBasis)
            : m_defBasis(defBasis), m_newBasis(newBasis)
        {}

        bool enter(gsDomainMap *map,
                   const BoxT &box,
                   gsDomainMap::NodeId node)
        {
            if( (*map)[node].isFork() )
                return true;
            const typename gsRemTypes<DIM>::tensorBasisT &locBasis=
                    *asBasis<DIM>(m_defBasis[(*map)[node].data.space]);
            indexDomT mult;
            for (short_t d = 0; d < DIM; ++ d)
                mult.row(d).setConstant(locBasis.degree(d)+1);
            knotRange knot[DIM];
            gsTensorIndexSubDomain<DIM> subDom;
            getTensorLevelLocalView(locBasis, box, mult, knot, subDom);
            m_newBasis.push_back(basisPtr(new gsTBS42(locBasis,subDom,(*map)[node].data.space,  gsDomainMap(locBasis.support(),box) )));
            (*map)[node].data.space=m_newBasis.size()-1;
            return true;
        }
    };
};

} // namespace gismo
