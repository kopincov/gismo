/** @file gsTBPN.h

    @brief Implementation of Truncated B-splines with Partially Nested
    hierarchy using gsRemappedBasis.
    These splines are a generalization of THB to the case in which the
    set of levels is not totally oredered. For instance, some levels can
    be refined only in one direction some in another.
    The hierarchy of the levels is described by a directed graph.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsRemappedBasis/gsAdaptiveSplineCommon.h>

namespace gismo {


/** Graph storing the ordering of the bases. If two bases are
 * connected with an edge, one of them is the refinement of the
 * other. */
template <int DIM>
class gsHierarchyGraph
{
    typedef typename gsRemTypes<DIM>::basisIdT      basisIdT;
    typedef typename gsRemTypes<DIM>::basisPtr      basisPtr;
public:
    gsMatrix<index_t> incidence;
    gsMatrix<index_t> connection;

    std::vector<basisIdT> coarsest;
    std::vector<basisIdT> fineset;

    gsHierarchyGraph()
    {}

    gsHierarchyGraph(const std::vector<basisPtr> &basis)
    {
        const basisIdT numBases=static_cast<basisIdT>(basis.size());
        incidence.setZero(numBases,numBases);

        for (basisIdT i=0; i<numBases-1;++i)
        {
            for (basisIdT j=0; j<numBases;++j)
            {
                if ( j==i)
                    continue;
                if ( isRefinement(basis[i], basis[j]) )
                    incidence(j,i)=1;
            }
        }
        connection=incidence;
        //        simplifyIncidence();
        //        computeConnections(incidence,connection);
    }

    void simplifyIncidence()
    {
        for (index_t c=0; c<incidence.cols();++c)
            removeIndirectPaths(c);
    }

    void removeIndirectPaths(index_t testColId)
    {
        gsMatrix<index_t>::Column col=incidence.col(testColId);

        for (index_t r=0; r<col.rows();++r)
        {
            if ( col(r)!=0 )
            {
                if (incidence(testColId,r)==1)
                    gsWarn<<"Spaces "<<r << " and "<<testColId<<" are the same\n"<<std::flush;
                else for (index_t c=0;c<incidence.cols();++c)
                {
                    if (c==testColId)
                        continue;
                    if ( incidence(r,c) && col(c) )
                    {
                        col(r)=0;
                        break;
                    }
                }
            }
        }
    }

    static void computeConnections (const gsMatrix<index_t> &inc, gsMatrix<index_t> &path)
    {
        // compute incidence to the numLvl power efficiently
        gsMatrix<index_t> tmp;

        path=inc;
        index_t log;
        index_t rem = inc.rows()-1;
        while (rem)
        {
            log = math::floor(log2(rem));
            rem = inc.rows()-(1<<log);

            tmp=inc;
            for (index_t exp=0; exp<log;++exp)
                tmp*=tmp;
            path*=tmp;
        }
    }

    static bool isRefinement(basisPtr coarser, basisPtr finer)
    {
        typename gsRemTypes<DIM>::tensorBasisT &coarserB=*asBasis<DIM>(coarser);
        typename gsRemTypes<DIM>::tensorBasisT &finerB=*asBasis<DIM>(finer);

        for (int i=0; i<coarserB.domainDim();++i)
        {
            int degDiff=finerB.degree(i)-coarserB.degree(i);
            if (degDiff<0)
                return false;
            if (!parKnots(degDiff,coarserB.knots(i), finerB.knots(i)))
                return false;
        }
        return true;
    }

    static bool parKnots(int degDiff, const gsKnotVector<real_t> &coarse, const gsKnotVector<real_t> &finer)
    {
        // check that all knots in par are in chl with multiplicity at least m + degDiff

        gsKnotVector<real_t>::const_uiterator parIt=coarse.ubegin();
        gsKnotVector<real_t>::const_uiterator parEn=coarse.uend();

        for (;parIt!=parEn;++parIt)
        {
            if (parIt.multiplicity() + degDiff > finer.multiplicity(*parIt))
                return false;
        }
        return true;
    }

};

/** Implementation of Truncated B-splines for partially nested
 * refinement. They were defined in \cite zore2016 and the
 * implementation is described in \cite bm2016. */
template <int DIM>
class gsTBPN : public gsAdaptiveSplineCommon<DIM>
{
    typedef typename gsAdaptiveSplineCommon<DIM>::tensorBasisT  tensorBasisT;
    typedef typename gsAdaptiveSplineCommon<DIM>::knotVectorT   knotVectorT;
    typedef typename knotVectorT::const_iterator                knotIter;
    typedef typename knotVectorT::reverse_iterator              backKnotIter;
    typedef typename gsAdaptiveSplineCommon<DIM>::basisPtr      basisPtr;
    typedef typename gsAdaptiveSplineCommon<DIM>::directionT    directionT;
    typedef typename gsAdaptiveSplineCommon<DIM>::basisIdT      basisIdT;

    using  gsAdaptiveSplineCommon<DIM>::m_basis;
    using  gsAdaptiveSplineCommon<DIM>::m_info;
    using  gsAdaptiveSplineCommon<DIM>::m_sele;
    using  gsAdaptiveSplineCommon<DIM>::m_repr;
    using  gsAdaptiveSplineCommon<DIM>::m_shift;

    typedef typename gsAdaptiveSplineCommon<DIM>::collectInfluencing collectInfluencing;
    typedef typename gsAdaptiveSplineCommon<DIM>::getMinAndMaxLevel  getMinAndMaxLevel;
    typedef typename gsAdaptiveSplineCommon<DIM>::indexSetT           indexSetT;
    typedef typename gsAdaptiveSplineCommon<DIM>::indexDomT          indexDomT;
    typedef typename gsAdaptiveSplineCommon<DIM>::indexVecT          indexVecT;
    typedef typename gsAdaptiveSplineCommon<DIM>::knotRange          knotRange;
    typedef typename gsAdaptiveSplineCommon<DIM>::tensorCoefs        tensorCoefs;

    using gsAdaptiveSplineCommon<DIM>::sanitizeBoxList;
    using gsAdaptiveSplineCommon<DIM>::getTensorFunctionInfo;
    using gsAdaptiveSplineCommon<DIM>::getTensorLevelLocalView;
    using gsAdaptiveSplineCommon<DIM>::updateCoefs;


    using gsRemappedBasis::checkDimAndInitShifts;
    using gsRemappedBasis::size;
public:
    gsTBPN (std::vector<basisPtr> basis, const gsBoxList &boxes)
    {
        m_info.first=DIM;
        m_info.second=1;

        m_basis=basis;

        gsBoxList mBoxList=boxes;
        sanitizeBoxList(mBoxList,-1);
        checkDimAndInitShifts();
        m_graph=gsHierarchyGraph<DIM>(m_basis);

        m_sele.initFromBoxesMax(mBoxList,asBasis<DIM>(m_basis[0])->support());
        makeRepresentation();

    }

    bool test()
    {
        bool passed=true;
        collectInfluencing influencing(m_basis);
        m_sele.patch(0).apply(influencing);
        gsMatrix<> rows=gsMatrix<>(m_repr.asMatrix()).rowwise().sum();
        for (size_t lvl=0; lvl<m_basis.size();++lvl)
        {
            typename indexSetT::const_iterator it = influencing.data[lvl].begin();
            typename indexSetT::const_iterator end = influencing.data[lvl].end();
            for(;it!=end;++it)
            {
                const real_t val=rows(*it+m_shift[lvl]);
                if ( math::abs(val-1)>0.000001)
                {
                    std::cout<<"lvl "<<lvl<< " fun "<<*it<<" sum to "<<val<< std::endl;
                    passed=false;
                }
            }
        }
        return passed;
    }


private:
    void makeRepresentation()
    {
        const int numLvl=static_cast<int>(m_basis.size());
        // per each level L lists all functions active on Omega_L
        collectInfluencing influencing(m_basis);
        m_sele.patch(0).apply(influencing);

        // Triplets in the representation matrix
        gsSparseEntries<real_t> repr;
        index_t count=0;

        for (int lvl=numLvl-1; lvl>=0; --lvl)
        {
            typename indexSetT::const_iterator it  = influencing.data[lvl].begin();
            typename indexSetT::const_iterator end = influencing.data[lvl].end();

            for (;it!=end;++it)
                analyzeFunction(static_cast<basisIdT>(lvl),*it, count,repr);
        }
        m_repr.asMatrix().resize(m_shift[m_basis.size()],count);
        m_repr.asMatrix().setFrom(repr);
    }

    void analyzeFunction(const basisIdT &lvl, const index_t &id, index_t &count, gsSparseEntries<real_t> &repr)
    {
        std::vector<indexSetT> slaves;
        getSlaveFunctions(lvl, id, slaves);
        for (index_t i=0; i<m_graph.connection.cols();++i)
        {
            if (static_cast<basisIdT>(i)!=lvl                              // Skip the level of the function we are analyzing.
                    && m_graph.connection(lvl,i)!=0 // Level i is coarser than level lvl
                    && slaves[i].size()>0)          // There are slave functions of level i
                return;                             // Function id from lvl is not selected.
        }
        repr.add(m_shift[lvl]+id, count, 1);
        std::vector<basisIdT> masterLvl(1,lvl);
        representFunction(lvl,id, slaves, count,repr,1, masterLvl, lvl);
        ++count;
    }



    void getSlaveFunctions(const basisIdT &lvl, const index_t &id, std::vector<indexSetT> & slaves)
    {
        gsMatrix<real_t,DIM,2> support = asBasis<DIM>(m_basis[lvl])->support(id); // function support as gsDomainMap
        gsDomainMap locLvl = getAreaOfLvl(support,lvl);   // representation of area at Lvl

        slaveGetter getter(m_basis, m_sele.patch(0), support);
        locLvl.apply(support,getter);
        std::swap(slaves,getter.data);
    }

    gsDomainMap getAreaOfLvl (const gsMatrix<real_t,DIM,2> &area, const basisIdT &lvl) const
    {
        gsDomainMap   result(m_sele.patch(0).getBoundingBox());
        result[result.root()].data.space=0;
        levelSelector select(result,lvl);
        m_sele.patch(0).apply(area, select);
        return result;
    }

    void representFunction(const basisIdT &lvl, const index_t &id, const std::vector<indexSetT> &slaves, index_t &count, gsSparseEntries<real_t> &repr, real_t fact,  std::vector<basisIdT> &masterLvl, const basisIdT preLvl)
    {
        // function description
        knotRange               oldK[DIM];
        gsMatrix<real_t,DIM,2>  support;
        indexDomT               mult;

        getTensorFunctionInfo   (id,lvl,support,mult,oldK);
        gsTensorIndexSubDomain<DIM> indexSubDom;
        tensorCoefs          coefs;
        std::vector<real_t>  tmpStorage;

        for (basisIdT dstLvl=0; dstLvl<static_cast<basisIdT>(m_graph.connection.rows());++dstLvl)
        {
            if (dstLvl==lvl) continue;
            if (m_graph.connection(dstLvl,lvl)!=0 && slaves[dstLvl].size()>0)
            {
                knotRange newK[DIM];
                getTensorLevelLocalView (dstLvl,support,mult,newK,indexSubDom);
                coefs.initOne           (indexSubDom.sizeM(),fact);
                tmpStorage.resize       (indexSubDom.sizeM().bottomRows(DIM-1).maxCoeff());
                updateCoefs             ( oldK, newK, coefs, tmpStorage);
                saveCoefs               (dstLvl, coefs, indexSubDom,slaves, count,repr, masterLvl, preLvl);
            }
        }
    }


    void saveCoefs(const basisIdT dstLvl,const tensorCoefs &coef, const gsTensorIndexSubDomain<DIM> &idSubDom, const std::vector<indexSetT> &slaves, index_t &count, gsSparseEntries<real_t> &repr,  std::vector<basisIdT> &masterLvl, const basisIdT preLevel)
    {

        typename indexSetT::const_iterator it=slaves[dstLvl].begin();
        typename indexSetT::const_iterator end=slaves[dstLvl].end();

        for (;it!=end;++it)
        {
            gsVector<index_t, DIM> tenId=idSubDom.toMulti(*it);
            if ( idSubDom.contains(tenId) )
            {
                real_t c=coef[tenId-idSubDom.lower()];
                if (c!=0)
                {
                    gsSparseEntries<real_t>::reverse_iterator pos = alreadyWritten(dstLvl, *it, count, repr, masterLvl.size());
                    if (pos!=repr.rbegin()+masterLvl.size() )
                    {
                        size_t n_pos=masterLvl.size()-(pos-repr.rbegin())-1;
                        if ( isSlave(masterLvl[n_pos],preLevel ) )
                        {
                            masterLvl[n_pos]=preLevel;
                            *pos=gsSparseEntries<real_t>::Triplet(pos->row(),pos->col(),c);
                        }
                        else if (masterLvl[n_pos]==preLevel)
                            *pos=gsSparseEntries<real_t>::Triplet(pos->row(),pos->col(),c+pos->value());
                        else
                        {
                        // GISMO_ASSERT( isSlave(preLevel,masterLvl[n_pos] ),"asdfa"  ); // this can happen in case of example 1
                        // gsWarn<< "Same function written through different paths\n";
                        }
                    }
                    else
                    {
                        masterLvl.push_back(preLevel );
                        repr.add(m_shift[dstLvl]+*it,count,c);
                        std::vector<indexSetT> slaves2;
                        getSlaveFunctions(dstLvl, *it, slaves2);
                        representFunction(dstLvl, *it, slaves2, count,repr,c, masterLvl, dstLvl);
                    }
                }
            }
        }
    }

    inline bool isSlave(basisIdT slave, basisIdT master)
    {
        return 0 != m_graph.connection(slave, master);
    }


    gsSparseEntries<real_t>::reverse_iterator alreadyWritten(const basisIdT lvl, const index_t id, const index_t &count, gsSparseEntries<real_t> &repr, size_t maxBack)
    {
        gsSparseEntries<real_t>::reverse_iterator it=repr.rbegin();
        gsSparseEntries<real_t>::reverse_iterator end=repr.rend();
        for (;it != end && it->col()==count; ++it)
            if (it->row()==m_shift[lvl]+id)
                return it;
        return repr.rbegin()+maxBack;
    }


    struct levelSelector: public gsDomainMap::actionBase<DIM,2>
    {
        typedef gsMatrix<real_t,DIM,2> BoxT;
        gsDomainMap                   &target;
        const basisIdT                &lvl;
        levelSelector ( gsDomainMap  &t, const basisIdT  &l)
            : target(t), lvl(l)
        {}
        bool enter  (const gsDomainMap *map, const gsMatrix<real_t,DIM,2> &bbox, gsDomainMap::NodeId &node )
        {
            if ((*map)[node].isFork())
                return true;
            if ((*map)[node].data.space!=lvl)
                return true;

            gsDomainMap::setBasisMax<DIM,2> setter(1);
            target.apply(bbox,setter);
            return true;
        }
    };


    struct slaveGetter : public gsDomainMap::actionBase<DIM,2>
    {
        typedef gsMatrix<real_t,DIM,2> BoxT;
        const std::vector<basisPtr>           &bases;
        gsDomainMap                   &target;
        gsMatrix<real_t,DIM,2>        &boundingBox;

        std::vector<indexSetT>  data;
        indexSetT               tmp;



        struct Collector: public gsDomainMap::actionBase<DIM,2>
        {
            typedef gsMatrix<real_t,DIM,2> BoxT;

            Collector (const std::vector<basisPtr> &bas, std::vector<indexSetT> &dat, boxSide s, indexSetT &t)
                :  bases(bas), data(dat), side(s), tmp(t)
            {}

            const std::vector<basisPtr> &bases;
            std::vector<indexSetT>      &data;
            boxSide side;
            indexSetT                   &tmp;

            bool enter (const gsDomainMap *map, const BoxT& box, const gsDomainMap::NodeId node )
            {
                if ((*map)[node].isFork())
                    return true;

                const basisIdT curLvl = (*map)[node].data.space;
                const tensorBasisT  &tensorBasis = *asBasis<DIM>(bases[curLvl]);

                gsMatrix<index_t, DIM,2> dom, subdom;
                for( int dir = 0; dir < DIM; ++dir )
                {
                    const index_t deg=tensorBasis.degree(dir);
                    const knotVectorT &knots=tensorBasis.knots(dir);
                    dom.row(dir)<< 0,  knots.size() - deg - 1;
                    const knotIter beg=knots.begin();
                    const knotIter end=knots.end();
                    real_t posB=box(dir,0);
                    real_t posE=box(dir,1);

                    if (side.direction()==dir )
                    {
                        if (side.parameter())
                            posE=posB;
                        else
                            posB=posE;
                    }
                    subdom(dir,0)=std::upper_bound(beg,end,posB) - beg;
                    subdom(dir,0)=std::max(static_cast<index_t>(0), subdom(dir,0)- deg - 1);
                    subdom(dir,1)=std::lower_bound(beg,end,posE) - beg;
                    subdom(dir,1)=std::max(subdom(dir,1), subdom(dir,0));
                }
                gsTensorIndexSubDomain<DIM> MSD(gsTensorIndexDomain<DIM>(dom),subdom);
                tmp.resize(MSD.size()+data[curLvl].size());
                typename indexSetT::const_iterator last=std::set_union(MSD.flat_begin(),MSD.flat_end(),data[curLvl].begin(),data[curLvl].end(),tmp.begin());
                tmp.resize(last-tmp.begin());
                std::swap(tmp,data[curLvl]);
                return true;
            }
        };

        slaveGetter(const std::vector<basisPtr> &bas, gsDomainMap &dom, gsMatrix<real_t,DIM,2>  &b)
            : bases(bas),target(dom), boundingBox(b)
        { init();}
        void init()
        {
            data.resize(bases.size());
            for (size_t i=0; i<bases.size();++i)
                data[i].resize(0);
        }
        bool enter  (const gsDomainMap *map, const gsMatrix<real_t,DIM,2> &bbox, gsDomainMap::NodeId &node )
        {
            if ((*map)[node].isFork())
                return true;
            if ((*map)[node].data.space==0)
                return true;

            gsMatrix<real_t,DIM,2> tBox=bbox;
            for (boxSide s=boxSide::getFirst(DIM); s<boxSide::getEnd(DIM);++s)
            {
                if ( tBox(s.direction(),s.parameter())==boundingBox(s.direction(),s.parameter()) )
                    continue;
                else
                {
                    const real_t a=tBox(s.direction(),  s.parameter() );
                    tBox(s.direction(), !s.parameter() ) = a;
                    tBox(s.direction(),  s.parameter() ) = s.parameter() ? math::nextafter(a,std::numeric_limits<real_t>::max()): math::nextafter(a,std::numeric_limits<real_t>::min());
                    Collector collect(bases, data, s,tmp);
                    target.apply(tBox, collect);
                    tBox.row(s.direction())=bbox.row(s.direction());
                }
            }
            return true;
        }
    };



protected:
    gsHierarchyGraph<DIM> m_graph;
    std::vector<index_t>  m_FromId;
    std::vector<basisIdT> m_FromLvl;

};

} // namespace
