/** @file gsHLR.h

    @brief Implementation of HLR-splines using gsRemappedBasis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsRemappedBasis/gsAdaptiveSplineCommon.h>

namespace gismo {

#define HLR_DEBUG 1

/** Implementation of the hierarchical LR-splines. HLR-splines are
 * defined in \cite bressan2015 and this implementation is described
 * in \cite bm2016, Section 4.4. */
template <int DIM>
class gsHLR : public gsAdaptiveSplineCommon<DIM>
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

public:
    using gsRemappedBasis::size;
    gsHLR (std::vector<basisPtr> basis, const gsBoxList &boxes)
        : m_activeT(m_basis)
    {
        m_info.first = DIM;
        m_info.second= 1;

        m_fullLvl=basis;
        m_basis=addIntermediateLevels(m_fullLvl);
        this->checkDimAndInitShifts();

        for (directionT dir =0; dir<DIM;++dir)
            m_deg(dir)=asBasis<DIM>(m_basis[0])->degree(dir);

        initDomains(boxes);
        makeRepresentation();
    }

    const gsDomainMap&      getDomain() const {return m_sele.patch(0);}

    int  degree(directionT dir=0) {return m_deg(dir);}

    void exportRequestToTex(std::string filename) const
    {
        m_reqDom.exportToTex(filename);
    }
    void exportDefinitionToTex(std::string filename) const
    {
        m_defDom.exportToTex(filename);
    }



protected:
    struct Contest
    {
        basisIdT   lvl;
        directionT dir;
    };

    std::vector<basisPtr> addIntermediateLevels(std::vector<basisPtr> basis )
    {
        std::vector<gsKnotVector<real_t> > knots(DIM);
        std::vector<basisPtr> result;

        result.push_back(basis[0]);
        for (directionT dir=0;dir<DIM;++dir)
            knots[dir]=asBasis<DIM>(basis[0])->knots(dir);
        directionT dir;
        for (size_t b=1;b<basis.size();++b)
        {
            for (dir=0;dir<DIM-1;++dir)
            {
                knots[dir]=asBasis<DIM>(basis[b])->knots(dir);
                result.push_back(basisPtr(new tensorBasisT(knots)));
            }
            knots[DIM-1]=asBasis<DIM>(basis[b])->knots(dir);
            result.push_back(basis[b]);
        }
        return result;
    }

        // we have three domains:
        //
        // the request domain: areas the user set to specific levels (where each level is refined in all directions)
        // the definition domain:areas obtained by the recursively applying the extension operator (here levels are sublevels)
        // the representation domain: areas used in the evaluation (here we use again levels of refinement in all directions)

    void initDomains(const gsBoxList &request)
    {
        gsMatrix<real_t> omega=asBasis<DIM>(m_basis[0])->support();

        gsBoxList mBoxList=request;
        this->sanitizeBoxList(mBoxList,0);
        m_reqDom.initFromBoxesMax(mBoxList,omega);

        initDefDomain ();
        initRepDomain ();
    }

    void initDefDomain ()
    {
        gsBoxList mBoxList=m_reqDom.asBoxList();
        gsMatrix<real_t,DIM,2> tmp;
        for (size_t b=mBoxList.size(); b>0; --b)
        {
            mBoxList.basisId(b-1)*=DIM;
            tmp=mBoxList.box(b-1);
            for (basisIdT lvl=mBoxList.basisId(b-1); lvl>0;--lvl)
            {
                extendBox(tmp,lvl-1, (lvl-1)%DIM);
                mBoxList.append(tmp,lvl-1);
            }
        }
        m_defDom.initFromBoxesMax(mBoxList,m_reqDom.patch(0).getBoundingBox());
    }

    void initRepDomain ()
    {
        gsBoxList mBoxList=m_defDom.asBoxList();
        gsMatrix<real_t,DIM,2> tmp;
        for (size_t b=mBoxList.size(); b>0; --b)
        {
            tmp=mBoxList.box(b-1);
            for (basisIdT lvl=mBoxList.basisId(b-1); lvl>0;--lvl)
            {
                extendBox(tmp,lvl-1, (lvl-1)%DIM);
                mBoxList.append(tmp,lvl);
            }
        }
        m_sele.initFromBoxesMax(mBoxList,m_defDom.patch(0).getBoundingBox());
    }


    void extendBox (gsMatrix<real_t,DIM,2> &box, basisIdT lvl, directionT dir)
    {
        const tensorBasisT &basis=*asBasis<DIM>(m_basis[lvl!=0? lvl-1:0]);

        knotIter kL = basis.knots(dir).begin();
        kL+=basis.knots(dir).uFind(box(dir,0)).lastAppearance()-m_deg(dir);
        knotIter kU = basis.knots(dir).begin();
        kU+=basis.knots(dir).uFind(box(dir,1)).firstAppearance()+m_deg(dir);
        knotIter kM = basis.knots(dir).end();
        box(dir,0)=*kL;
        if ( box(dir,1)<= *kU )
            box(dir,1)=*kU;
        else if (kU<kM)
            box(dir,1)=*(kU);
        else GISMO_ERROR("???");
    }


    struct project1D : public gsDomainMap::actionBase<DIM,2>
    {
        typedef gsMatrix<real_t,DIM,2> BoxT;

        gsDomainMap &target;
        directionT  dir;

        project1D(gsDomainMap &other, directionT d)
            : target(other), dir(d)
        {}

        bool enter  (const gsDomainMap *map, const gsMatrix<real_t,DIM,2> &bbox, gsDomainMap::NodeId &node )
        {
            if ((*map)[node].isFork())
                return true;
            typename gsDomainMap::setBasisMin<1,2>::BoxT boxC=bbox.row(dir);
            gsDomainMap::setBasisMin<1,2> action((*map)[node].data.space);
            target.apply(boxC, action);
            return true;
        }
    };

    void getLineSegments (const gsMatrix<real_t,DIM,2> & support, const basisIdT &lvl, gsBoxList &result)
    {
        const directionT dir = lvl%DIM;
        gsDomainMap seg(support.row(dir));
        gsMatrix<real_t,DIM,2> dom=support;
        dom.row(dir)=m_defDom.patch(0).getBoundingBox().row(dir);
        project1D action(seg,dir);
        m_defDom.patch(0).apply(dom,action);
        result=seg.asBoxList();
    }

    void collectFunctionsFromLine (const indexVecT &pos, const basisIdT &lvl, index_t &count, gsSparseEntries<real_t> &repr)
    {
        if (lvl!=m_basis.size()-1)
            collectsFunctionLine (pos, lvl, count, repr);
        else
            collectsFunctionLineFinest (pos, lvl, count, repr);
    }

    enum {
        ST_E = 0, // no knots yet
        ST_M = 3, // both coarse and fine knots
        ST_C = 1, // only coarse
        ST_F = 2  // only fine
    };

    void collectsFunctionLine (const indexVecT &pos, const basisIdT &lvl, index_t &count, gsSparseEntries<real_t> &repr)
    {
        const directionT        dir=lvl%DIM;
        knotRange               funK[DIM];

        // this function description
        gsMatrix<real_t,DIM,2>   supp;
        indexDomT                mult;
        this->getTensorFunctionInfo   (*asBasis<DIM>(m_basis[lvl]),pos,  supp, mult, funK);

        supp.row(dir)=m_defDom.patch(0).getBoundingBox().row(dir);
        gsBoxList segments(1);
        getLineSegments(supp,lvl,segments);

        knotRange coarse;
        knotRange fine;
        coarse.beg=asBasis<DIM>(m_basis[lvl])->knots(dir).begin();
        coarse.end=asBasis<DIM>(m_basis[lvl])->knots(dir).end();
        fine.beg=asBasis<DIM>(m_basis[lvl+1])->knots(dir).begin();
        fine.end=asBasis<DIM>(m_basis[lvl+1])->knots(dir).end();

        std::vector<real_t> actVal;
        basisIdT  curLvl=0;
        basisIdT  preLvl=0;
        index_t   numF=0;
        index_t   numC=0;

        int  state = ST_E;
        for(size_t s=0; s <segments.size(); ++s)
        {
            curLvl=segments.basisId(s);
            const real_t &sBeg=segments.box(s)(0,0);
            const real_t &sEnd=segments.box(s)(0,1);

            // handling interface
            const basisIdT begLvl=math::max(curLvl,preLvl);
            if (begLvl==lvl)
            {
                numC  += addSegmentLim(sBeg,coarse,actVal);
                numF   = 0;
                state |= ST_C;
            }
            else if (begLvl>lvl)
            {
                numF += addSegmentLim(sBeg,fine,actVal);
                numC  = 0;
                state |= ST_F;
            }
            representIfGreater(pos,lvl,state,numF,numC, funK,actVal,count,repr );

            // handling interior
            if (curLvl==lvl)
            {
                numC += addSegmentInt(sBeg,sEnd,coarse,actVal);
                numF  = 0;
                state |= ST_C;
            }
            else if (curLvl>lvl)
            {
                numF += addSegmentInt(sBeg,sEnd,fine,actVal);
                numC  = 0;
                state |= ST_F;
            }
            representIfGreater(pos,lvl,state,numF,numC, funK,actVal,count,repr );

            if (curLvl<lvl)
            {
                representStride (lvl, funK, actVal, count, repr);
                actVal.resize(0);
                state=ST_E;
                numC=0;
                numF=0;
            }

            preLvl=curLvl;
        }
        if (curLvl==lvl)
        {
            addSegmentLim(supp(dir,1),coarse, actVal );
            state |= ST_C;
        }
        else if (curLvl==lvl+1)
        {
            numF+=addSegmentLim(supp(dir,1),fine, actVal );
            representIfGreater(pos,lvl,state,numF,numC, funK,actVal,count,repr );
            state |= ST_F;
        }
        representStride (lvl, funK, actVal, count, repr);
    }


    void collectsFunctionLineFinest (const indexVecT &pos, const basisIdT &lvl, index_t &count, gsSparseEntries<real_t> &repr)
    {
        const directionT        dir=lvl%DIM;
        knotRange               funK[DIM];

        gsMatrix<real_t,DIM,2>   supp;
        indexDomT                mult;
        this->getTensorFunctionInfo   (*asBasis<DIM>(m_basis[lvl]),pos, supp, mult, funK);
        supp.row(dir)=m_defDom.patch(0).getBoundingBox().row(dir);

        gsBoxList segments(1);
        getLineSegments(supp,lvl,segments);

        knotRange coarse;
        coarse.beg=asBasis<DIM>(m_basis[lvl])->knots(dir).begin();
        coarse.end=asBasis<DIM>(m_basis[lvl])->knots(dir).end();

        std::vector<real_t> actVal;
        basisIdT  curLvl=0;
        basisIdT  preLvl=0;

        int  state = ST_E;
        for(size_t s=0; s <segments.size(); ++s)
        {
            curLvl=segments.basisId(s);
            const real_t &sBeg=segments.box(s)(0,0);
            const real_t &sEnd=segments.box(s)(0,1);

            // handling interface
            const basisIdT begLvl=math::max(curLvl,preLvl);
            if (begLvl==lvl)
            {
                addSegmentLim(sBeg,coarse,actVal);
                state |= ST_C;
            }
            // handling interior
            if (curLvl==lvl)
            {
                addSegmentInt(sBeg,sEnd,coarse,actVal);
                state |= ST_C;
            }

            if (curLvl<lvl)
            {
                representStride (lvl, funK, actVal, count, repr);
                actVal.resize(0);
                state=ST_E;
            }

            preLvl=curLvl;
        }
        if (curLvl==lvl)
        {
            addSegmentLim(supp(dir,1),coarse, actVal );
            state |= ST_C;
        }
        representStride (lvl, funK, actVal, count, repr);
    }



    index_t addSegmentLim(const real_t &pos, const knotRange &orig, std::vector<real_t> &knots)
    {
        typename knotRange::iter itt = orig.beg;
        typename knotRange::iter end = orig.end;
        itt=std::lower_bound(itt,end,pos);
        end=std::upper_bound(itt,end,pos);
        std::copy(itt,end,std::back_inserter(knots));
        return end-itt;
    }

    index_t addSegmentInt(const real_t &posB, const real_t &posE, const knotRange &orig, std::vector<real_t> &knots)
    {
        typename knotRange::iter itt = orig.beg;
        typename knotRange::iter end = orig.end;
        itt=std::upper_bound(itt,end,posB);
        end=std::lower_bound(itt,end,posE);
        std::copy(itt,end,std::back_inserter(knots));
        return end-itt;
    }

    void representIfGreater(
            const indexVecT & /*pos*/,
            const basisIdT &lvl,
            int &st,
            index_t &numF,
            index_t &numC,
            const knotRange (&funK)[DIM],
            std::vector<real_t> &knots,
            index_t &count,
            gsSparseEntries<real_t> &repr)
    {
        const directionT    dir = lvl%DIM;
        const int           deg = m_deg(dir);

        if (numF>deg+1)
        {
            std::vector<real_t> tmp(knots.end()-deg-1, knots.end());
            if (st!=ST_F)
            {
                knots.resize(knots.size()-numF+deg+1);
                representStride (lvl,funK, knots,count,repr);
            }
            numF=deg+1;
            numC=0;
            knots=tmp;
            st=ST_F;
        }


    }

    void representStride (const basisIdT &lvl, const knotRange (&funK)[DIM], const std::vector<real_t> &knots,index_t &count, gsSparseEntries<real_t> &repr)
    {
        const directionT    dir = lvl%DIM;
        const int           deg = m_deg(dir);
        indexDomT           mult;

        if (knots.size()<static_cast<size_t>(deg+2))
            return;

        knotRange mFunK[DIM];
        std::copy(funK,funK+DIM,mFunK);

        for (size_t b=0; b<knots.size()-deg-1;++b)
        {
            mFunK[dir].beg = knots.begin()+b;
            mFunK[dir].end = mFunK[dir].beg+deg+2;

            gsMatrix<real_t,DIM,2> support;
            this->getTensorFunctionInfo(mFunK,support,mult);

            getMinAndMaxLevel act(lvl);
            m_sele.patch(0).apply(support,act);

            if (act.minLvl>= lvl)
            {
                representFunction(mFunK,act.minLvl,act.maxLvl,count,repr);
                ++count;
            }
        }
    }


    void makeRepresentation()
    {
        const basisIdT numLvl=static_cast<basisIdT>(m_basis.size());

        gsTensorIndexDomain<DIM>    idxDom;
        gsTensorIndexSubDomain<DIM> idxSub;

        indexDomT                   lim;

        m_activeT=collectInfluencing(m_basis);
        m_sele.patch(0).apply(m_activeT);

        gsSparseEntries<real_t> repr;
        index_t count=0;
        for (basisIdT lvl=0; lvl<numLvl; ++lvl)
        {
            const int dir = lvl%DIM;
            idxDom = this->getTensorIndexDomainForLevel(lvl);
            lim=static_cast<indexDomT>(idxDom);
            lim(dir,1)=1;
            idxSub=gsTensorIndexSubDomain<DIM>(idxDom,lim);

            typename gsTensorIndexSubDomain<DIM>::const_iterator lineI=idxSub.begin();
            typename gsTensorIndexSubDomain<DIM>::const_iterator lineE=idxSub.end();

            for (;lineI!=lineE;++lineI)
                collectFunctionsFromLine(*lineI,lvl, count, repr);
        }
        m_repr.asMatrix().resize(m_shift[m_basis.size()],count);
        m_repr.asMatrix().setFrom(repr);
    }


    void representFunction(
            // function description
            const knotRange (&funK)[DIM],
            // saving options
            const basisIdT minLvl,
            const basisIdT maxLvl,
            // representation
            const index_t count,
            gsSparseEntries<real_t> &repr)
    {
        knotRange newK[DIM];
        knotRange oldK[DIM];
        std::copy(funK,funK+DIM,oldK);

        gsMatrix<real_t,DIM,2>  support; // the support of the function is the cartesian products of the intervals described by the rows
        indexDomT               mult;    // mult(dir,i) is the multiplicity of support(dir,i) in the knot vector of the function

        gsTensorIndexSubDomain<DIM> indexSubDom;

        this->getTensorFunctionInfo(oldK,support,mult);
        gsAdaptiveSplineCommon<DIM>::getTensorLevelLocalView (maxLvl,support,mult,newK,indexSubDom);

        // pre-allocate big data structures
        tensorCoefs coefs;
        coefs.initOne(indexSubDom.sizeM());
        std::vector<real_t>  tmpStorage(indexSubDom.sizeM().bottomRows(DIM-1).maxCoeff());

        // compute the representation in the next levels
        for (basisIdT dstLvl=minLvl; dstLvl<=maxLvl; ++dstLvl)
        {
            gsAdaptiveSplineCommon<DIM>::getTensorLevelLocalView(dstLvl,support,mult,newK,indexSubDom);
            gsAdaptiveSplineCommon<DIM>::updateCoefs     ( oldK, newK, coefs, tmpStorage);
            saveCoefs       (m_activeT.data[dstLvl],m_shift[dstLvl],count,indexSubDom,coefs,repr);
            std::copy(newK,newK+DIM,oldK); // updates old knots
        }
    }

    void static saveCoefs (
            const indexSetT                    &influencing,
            const index_t                      rowShift,
            const index_t                      column,
            const gsTensorIndexSubDomain<DIM> &idSubDom,
            tensorCoefs                       &coefs,
            gsSparseEntries<real_t>           &repr
            )
    {
        typename indexSetT::const_iterator inflIt = std::lower_bound(influencing.begin(),influencing.end(),idSubDom.minFlat());
        typename indexSetT::const_iterator inflEn = std::upper_bound(inflIt,influencing.end(), idSubDom.endFlat());

        while (inflIt!=inflEn )
        {
            gsVector<index_t, DIM>  tenId = idSubDom.toMulti(*inflIt);
            if (idSubDom.contains(tenId) )
            {
                real_t c=coefs[tenId-idSubDom.lower()];
                if (c!=0) repr.add(rowShift+*inflIt,column,c);
            }
            ++inflIt;
        }
    }


protected:
    gsSelector             m_defDom;
    gsSelector             m_reqDom;
    gsVector<int,DIM>      m_deg;
    std::vector<basisPtr>  m_fullLvl;

    std::vector<index_t>   m_fromLvl;

    collectInfluencing     m_activeT;
};





} // namespace gismo

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlUtils.h>
#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo {

namespace internal
{

template <int d>
gsHLR<d> *getHLRFromXML(gsXmlNode * node)
{
    // TODO add attribute to differentiate between dyadic and provided
    // bases list

    // Initialize the HBSplineBasis
    std::istringstream str;
    gsXmlNode * tmp;

    // Insert all boxes
    unsigned c;
    typename gsHLR<d>::basisIdT maxLevel=0;
    std::vector<unsigned int> all_boxes;
    for (tmp = node->first_node("box");
         tmp; tmp = tmp->next_sibling("box"))
    {
        typename gsHLR<d>::basisIdT level(atoi( tmp->first_attribute("level")->value() ));
        maxLevel=math::max(level,maxLevel);
        all_boxes.push_back(level);
        str.clear();
        str.str( tmp->value() );
        for( unsigned i = 0; i < 2*d; i++)
        {
            gsGetInt(str, c);
            all_boxes.push_back(c);
        }
    }

    // Construct the bases
    tmp = node->first_node("Basis");
    GISMO_ASSERT( tmp , "Expected to find a basis node.");
    std::vector<typename gsHLR<d>::basisPtr> bases;
    gsTensorBSplineBasis<d,real_t> *last=gsXml<gsTensorBSplineBasis<d,real_t> >::get(tmp);
    bases.push_back( basisPtr(last) );
    for (typename gsHLR<d>::basisIdT lvl=1; lvl <= maxLevel;++lvl)
    {
        last = last->clone();
        last->uniformRefine();
        bases.push_back(basisPtr(last));
    }

    gsBoxList boxes(bases,bases[0]->domainDim(),all_boxes);
    return new gsHLR<d>(bases,boxes);
}


template <int DIM>
gsXmlNode *putHLRFromXML(const gsHLR<DIM> & obj,
                         gsXmlTree & data)
{
    // Add a new node (without data)
    gsXmlNode* tp_node = internal::makeNode("Basis" , data);

    std::string typeName;
    typeName+="HLRSplineBasis";
    typeName+= DIM>1 ? to_string(DIM) :"";
    tp_node->append_attribute( makeAttribute("type", typeName, data) );

    // Write the component bases
    std::vector<typename gsHLR<DIM>::tensorBasisT*> bases = obj.getTensorBases();

    gsXmlNode * tmp = putTensorBasisToXml(*bases[0], data);
    tp_node->append_node(tmp);

    //Output boxes
    gsBoxList boxes = obj.getSelector().asBoxList();

    std::vector<index_t> boxEntries;
    gsMatrix<index_t> boxInd;
    for (size_t b=0; b< boxes.size(); ++b)
    {
        if (boxes.basisId(b)==0)
            continue;
        gsBoxList::toRefineElementFormat(static_cast<const gsFunctionSet<real_t>*>(bases[boxes.basisId(b)]), boxes.basisId(b), boxes.box(b),boxEntries);
        boxInd=gsAsConstMatrix<index_t>(boxEntries.data()+1,1,2*DIM);
        tmp = putMatrixToXml( boxInd , data, "box" );
        tmp->append_attribute( makeAttribute("level", to_string(boxEntries[0]), data ) );
        tp_node->append_node(tmp);
        boxEntries.clear();
    }
    // All set, return the basis
    return tp_node;
}

/// Get a Truncated Hierarchical B-spline basis from XML data
template<int d>
class gsXml< gsHLR<d> >
{
    typedef gsHLR<d> Object;
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(Object)
    static std::string tag  () { return "Basis"; }
    static std::string type () { return "HLRSplineBasis"+ (d>1 ? to_string(d):""); }

    static gsHLR<d> * get (gsXmlNode * node)
    {
        GISMO_ASSERT(
                    ( !strcmp( node->name(),"Basis") ) && ( !strcmp(node->first_attribute("type")->value(), internal::gsXml<gsHLR<d> >::type().c_str() ) ),
                    "Something is wrong with the XML data: There should be a node with a "<<
                    internal::gsXml<Object>::type().c_str()<<" Basis.");
        return getHLRFromXML<d>(node);
    }

    static gsXmlNode * put (const gsHLR<d> & obj,
                            gsXmlTree & data )
    {
        return putHLRFromXML<d>(obj,data);
    }
};


} // internal




} // namespace gismo
