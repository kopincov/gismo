/** @file gsTHB.h

    @brief New implementation of THB-splines using gsRemappedBasis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsRemappedBasis/gsAdaptiveSplineCommon.h>

namespace gismo {

    /** Implementation of THB-splines as described in \cite bm2016,
     * Section 4.1, Implementation 1. */

// TODO reduce memory usage in the case of many levels by
//      decoupling the bases used in the evaluation from the
//      the basis used in the definition.
//      The problem is that all functions of all levels get a row in the
//      representation matrix so that for the common case of and L levels
//      dyadic refinement the matrix as at least 2^(L*DIM)*(1-1/2^(DIM))^-1 (exact for degree 0).
//      Every function costs at least a pointer to the correspondong row
//      in the sparse matrix data structure in 2d 15 levels is over the
//      maximum we can do with 32 index_t; an empty representation matrix
//      would need 2^3 * 2*30 * 4/3 = 2^35/3 bytes = 11GB
//      A possible solution is to use more, but smaller spaces. even oner
//      per leaf of the domain map but containing only the influencing functions
//      so that all the empty rows are not used.

template <int DIM>
class gsTHB : public gsAdaptiveSplineCommon<DIM>, public gsBoxRefinableInterface
{
public:
    typedef memory::shared_ptr<gsTHB> Ptr;
    typedef memory::unique_ptr<gsTHB> uPtr;

    typedef typename gsAdaptiveSplineCommon<DIM>::tensorBasisT  tensorBasisT;
    typedef typename gsAdaptiveSplineCommon<DIM>::knotVectorT   knotVectorT;
    typedef typename knotVectorT::const_iterator                knotIter;
    typedef typename knotVectorT::reverse_iterator              backKnotIter;
    typedef typename gsAdaptiveSplineCommon<DIM>::basisPtr      basisPtr;
    typedef typename gsAdaptiveSplineCommon<DIM>::directionT    directionT;
    typedef typename gsAdaptiveSplineCommon<DIM>::basisIdT      basisIdT;

private:
    using  gsAdaptiveSplineCommon<DIM>::m_basis;
    using  gsAdaptiveSplineCommon<DIM>::m_info;
    using  gsAdaptiveSplineCommon<DIM>::m_sele;
    using  gsAdaptiveSplineCommon<DIM>::m_repr;
    using  gsAdaptiveSplineCommon<DIM>::m_shift;

    typedef typename gsAdaptiveSplineCommon<DIM>::collectInfluencing collectInfluencing;
    typedef typename gsAdaptiveSplineCommon<DIM>::getMinAndMaxLevel  getMinAndMaxLevel;
    typedef typename gsAdaptiveSplineCommon<DIM>::indexSetT          indexSetT;
    typedef typename gsAdaptiveSplineCommon<DIM>::indexDomT          indexDomT;
    typedef typename gsAdaptiveSplineCommon<DIM>::indexVecT          indexVecT;
    typedef typename gsAdaptiveSplineCommon<DIM>::knotRange          knotRange;
    typedef typename gsAdaptiveSplineCommon<DIM>::tensorCoefs        tensorCoefs;

public:
    using gsRemappedBasis::size;

public:
    gsTHB (std::vector<basisPtr> _basis, const gsBoxList &boxes, bool truncate=true)
        : m_truncate(truncate)
    {
        m_info.first  = DIM;
        m_info.second = 1;

        m_basis=_basis;
        m_selected.resize(m_basis.size());

        gsBoxList mBoxList=boxes;
        this->sanitizeBoxList(mBoxList,-1);
        this->checkDimAndInitShifts();

        m_sele.initFromBoxesMax(mBoxList,asBasis<DIM>(m_basis[0])->support());
        makeRepresentation();
    }

    const gsDomainMap&      getDomain() const {return m_sele.patch(0);}

    std::vector<typename gsRemTypes<DIM>::tensorBasisPtr> getTensorBases() const
    {
        std::vector<typename gsRemTypes<DIM>::tensorBasisPtr>  result;
        for (size_t b=0; b<m_basis.size();++b)
        {
            result.push_back(asBasis<DIM>(m_basis[b]));
        }
        return result;
    }

    /// Returns the degree of the basis with basisId 0 in the direction dir.
    short_t degree(directionT dir=0) {return asBasis<DIM>(m_basis[0])->degree(dir);}

    /**
       @brief levelAndIndex
        given a function index \a funcId returns the level \a level and the in-level-index \a inLevelIndex
        of the originating function.
       @param[in]  funcId
       @param[out] level
       @param[out] inLevelIndex
     */
    void levelAndIndexOf (index_t funcId, basisIdT &level, index_t &inLevelIndex) const
    {
        GISMO_ASSERT (0<=funcId && funcId <size() ,"Function id is out of bounds: you asked for "<<funcId<<" but the space has size "<< size() << ".");
        level=0;
        inLevelIndex=funcId;
        for ( ; m_selected[level].size()>static_cast<size_t>(funcId); ++level)
            inLevelIndex-=m_selected[level].size();
    }

    /**
       @brief functionIndexOf
        Given a level \a _level and an in-level-function index \a inLevelIndex,
        it returns the corresponding function index
        or -1 if the specified function is not selected.
       @param _level
       @param inLevelIndex
     */
    index_t functionIndexOf (basisIdT _level, index_t inLevelIndex) const
    {
        const basisIdT numLvl=static_cast<basisIdT>(m_basis.size());
        if (_level>=numLvl || _level<0 )
            return -1; // level does not exists
        std::vector<index_t>::const_iterator beg=m_selected[_level].begin();
        std::vector<index_t>::const_iterator end=m_selected[_level].end();
        std::vector<index_t>::const_iterator it =beg;
        it = std::lower_bound(beg,end,inLevelIndex);
        if (it != end && *it!=inLevelIndex)
            return -1; // not found
        index_t result=it-beg;
        for (basisIdT l=_level;l>0; --l)
            result+=m_selected[l-1].size();
        return result;
    }

    /// Returns a pointer to the basis of basisId \a lvl.
    const tensorBasisT* getBasisOfLevel (basisIdT lvl) const
    {
        if (0<= lvl && lvl<m_basis.size())
            return asBasis<DIM>(m_basis[lvl]);
        return NULL;
    }
public:

    /// Refines the selector \a m_sele based on \a boxes (in the format of \a gsBoxList).
    virtual void refine (const gsBoxList &boxes)
    {
        gsBoxList mBoxList=boxes;
        this->sanitizeBoxList(mBoxList,-1);
        m_sele.addBoxesMax(mBoxList);
        m_selected.resize(0);
        m_selected.resize(m_basis.size());
        makeRepresentation();
    }

    virtual void refineWithTransfer (const gsBoxList &boxes, gsSparseMatrix<real_t> &transfer)
    {
        std::vector<indexSetT> oldSelected(m_basis.size());
        for (size_t i=0; i< m_basis.size(); ++i)
            std::swap(oldSelected[i],m_selected[i]);

        gsBoxList mBoxList=boxes;
        this->sanitizeBoxList(mBoxList,-1);
        m_sele.addBoxesMax(mBoxList);
        makeRepresentation();

        m_sele.exportToTex("mshAfterRef");
        int err=system ("pdflatex mshAfterRef.tex >>/dev/null");
        GISMO_UNUSED(err);

        makeTransfer(oldSelected,transfer);
    }

    virtual void refineElements(std::vector<index_t> const & elem)
    {
        gsBoxList boxes(m_basis,DIM,elem);
        refine(boxes);
    }

    virtual void refineElementsWithTransfer(std::vector<index_t> const & elem, gsSparseMatrix<real_t> &tran)
    {
        gsBoxList boxes(m_basis,DIM,elem);
        refineWithTransfer(boxes, tran);
    }

protected:
    void makeTransfer (const std::vector<indexSetT> &oldSel, gsSparseMatrix<real_t> &transfer)
    {
        const basisIdT numLvl=static_cast<basisIdT>(m_basis.size());
        gsSparseEntries<real_t> tranT;
        index_t count=0;
        index_t trPos=0;

        getMinAndMaxLevel desc;
        for (basisIdT lvl=0; lvl<numLvl; ++lvl)
        {
            typename indexSetT::const_iterator it  = oldSel[lvl].begin();
            typename indexSetT::const_iterator end = oldSel[lvl].end();

            for (;it!=end;++it)
            {
                typename indexSetT::const_iterator pos, posE;
                pos = std::lower_bound(m_selected[lvl].begin(), m_selected[lvl].end(),*it);
                posE = m_selected[lvl].end();
                if ( pos!=posE  &&  *pos==*it )
                    tranT.add(pos-m_selected[lvl].begin()+m_lvlShifts[lvl],count,1);
                trPos=tranT.size();
                desc.maxLvl=lvl;
                desc.minLvl=lvl;
                desc.fromLvl=lvl;
                m_sele.patch(0).apply(asBasis<DIM>(m_basis[lvl])->support(*it),desc);
                representFunction( *it, lvl, desc.maxLvl, m_selected, oldSel, count, tranT);
                ++count;
                changeRowIndices(lvl,trPos,tranT);
            }
        }
        transfer.resize(m_repr.getNrOfTargets(), count);
        transfer.setFrom(tranT);
    }

    void changeRowIndices (const basisIdT &lvl, const index_t &trPos, gsSparseEntries<real_t> &tranT)
    {
        basisIdT curLvl=lvl;
        index_t id;
        for (size_t i=trPos; i<tranT.size(); ++i)
        {
            while (tranT[i].row()>m_shift[curLvl+1])
                ++curLvl;
            id=tranT[i].row()-m_shift[curLvl];
            tranT[i]=gsSparseEntries<real_t>::Triplet(functionIndexOf(curLvl,id), tranT[i].col(), tranT[i].value());
        }
    }


    void makeRepresentation()
    {
        const int numLvl=static_cast<int>(m_basis.size());
        // per each level L lists all functions active on Omega_L
        collectInfluencing influencing(m_basis);
        m_sele.patch(0).apply(influencing);

        // Triplets in the representation matrix
        gsSparseEntries<real_t> repr;

        index_t count=0;//Number of selected functions encountered in the current level so far.
        m_lvlShifts.resize(numLvl+1);

        // we build the space from the finest space
        // and then we invert the indices
        // to be compatible with the stable implementation
        // this is index HELL
        std::vector<index_t> lvlSize  (numLvl);
        std::vector<size_t>  tripletBeg;

        getMinAndMaxLevel desc;
        for (int lvl=numLvl-1; lvl>=0; --lvl)
        {
            count=0;
            tripletBeg.push_back(repr.size());
            typename indexSetT::const_iterator it  = influencing.data[lvl].begin();
            typename indexSetT::const_iterator end = influencing.data[lvl].end();

            for (;it!=end;++it)
            {
                desc.maxLvl = lvl;
                desc.minLvl = lvl;
                desc.fromLvl= lvl;
                //Write into desc the min and max level encoutered inside the support of the function with index *it:
                m_sele.patch(0).apply(asBasis<DIM>(m_basis[lvl])->support(*it),desc);

                if (desc.minLvl==static_cast<basisIdT>(lvl) )// the function with index *it from level lvl is selected because its support does not intersect any coarser level
                {
                    m_selected[lvl].push_back(*it);
                    repr.add(m_shift[lvl]+*it, count, 1);// It is represented in level lvl by itself with coeff 1.
                    if (desc.maxLvl> static_cast<basisIdT>(lvl) )
                        representFunction( *it, lvl, desc.maxLvl, influencing.data, m_selected , count, repr);
                    ++count;
                }
            }
            lvlSize[lvl]=count;
        }
        tripletBeg.push_back(repr.size());
        std::reverse(tripletBeg.begin(),tripletBeg.end());


        typedef gsSparseEntries<real_t>::Triplet Triplet;
        m_lvlShifts.resize(numLvl+1);
        m_lvlShifts[0]=0;
        for (basisIdT lvl=0;lvl<static_cast<basisIdT>(m_basis.size());++lvl)
        {
            m_lvlShifts[lvl+1]=m_lvlShifts[lvl]+lvlSize[lvl];
            for (size_t t=tripletBeg[lvl+1];t<tripletBeg[lvl];++t)
                repr[t]=Triplet(repr[t].row(), repr[t].col()+m_lvlShifts[lvl], repr[t].value());
        }
        m_repr.asMatrix().resize(m_shift[m_basis.size()],m_lvlShifts[numLvl]);
        m_repr.asMatrix().setFrom(repr);
    }

protected:
    std::vector<index_t>    m_lvlShifts;
    std::vector<indexSetT>  m_selected;
    bool                    m_truncate;
private:
    void representFunction(
            // function description
            const index_t  idF,
            const basisIdT lvl,
            // saving options
            const basisIdT maxLvl,
            const std::vector<indexSetT> & infl,
            const std::vector<indexSetT> & sele,
            // representation
            const index_t pos,
            gsSparseEntries<real_t> &repr)
    {
        knotRange oldK[DIM];
        knotRange newK[DIM];

        gsMatrix<real_t,DIM,2>  support; // the support of the function is the cartesian products of the intervals described by the rows
        indexDomT               mult;    // mult(dir,i) is the multiplicity of support(dir,i) in the knot vector of the function

        gsTensorIndexSubDomain<DIM> indexSubDom;

        this->getTensorFunctionInfo   (idF,lvl,support,mult,oldK);
        this->getTensorLevelLocalView (maxLvl,support,mult,newK,indexSubDom);

        // pre-allocate big data structures
        tensorCoefs coefs;
        coefs.initOne(indexSubDom.sizeM());
        std::vector<real_t> tmpStorage(indexSubDom.sizeM().bottomRows(DIM-1).maxCoeff());

        // compute the representation in the next levels using Boehm refinement and truncation
        for (basisIdT dstLvl=lvl+1; dstLvl<=maxLvl; ++dstLvl)
        {
            this->getTensorLevelLocalView(dstLvl,support,mult,newK,indexSubDom);
            this->updateCoefs     ( oldK, newK, coefs, tmpStorage);
            if (m_truncate)
                truncateAndSave (infl[dstLvl],sele[dstLvl],m_shift[dstLvl],pos,indexSubDom,coefs,repr);
            else
                saveOnly(infl[dstLvl],sele[dstLvl],m_shift[dstLvl],pos,indexSubDom,coefs,repr);
            std::copy(newK,newK+DIM,oldK); // updates old knots
        }
    }

    static void truncateAndSave (
            const indexSetT                    &influencing,
            const indexSetT                    &selected,
            const index_t                      rowShift,
            const index_t                      column,
            const gsTensorIndexSubDomain<DIM> &idSubDom,
            tensorCoefs                       &coefs,
            gsSparseEntries<real_t>           &repr
            )
    {
        // this function should be called with all data from the same level, that we assume to be L

        // influencing contains the indices of the functions of level L active on Omega_L
        // selected contains the indices of the functions of level L selected by the Kraft procedure
        // rowShift and column specify where to write the coefficients to save
        // idDomain converts between multi-indices and indices associated to functions of level L
        // idSubDom represents the box of indices that must be considered

        // coefs    are the full coefficients, those of selected functions are set to 0 by this function
        // repr     is a collection of triplets, this function adds one entry for each influencing non-selected function
        //          whose coefficient in coefs is non zero

        // truncate the coefficients of selected and save those of influencing

        typename indexSetT::const_iterator selIt = std::lower_bound(selected.begin(), selected.end(),idSubDom.minFlat());
        typename indexSetT::const_iterator selEn = std::upper_bound(selIt,selected.end(),idSubDom.endFlat());

        typename indexSetT::const_iterator inflIt = std::lower_bound(influencing.begin(),influencing.end(),idSubDom.minFlat());
        typename indexSetT::const_iterator inflEn = std::upper_bound(inflIt,influencing.end(), idSubDom.endFlat());

        while (inflIt!=inflEn && selIt!=selEn)
        {
            if (*inflIt<*selIt)
            {  // write influencing coefficient
                gsVector<index_t, DIM>  tenId = idSubDom.toMulti(*inflIt);
                if (idSubDom.contains(tenId) )
                {
                    real_t c=coefs[tenId-idSubDom.lower()];
                    if (c!=0) repr.add(rowShift+*inflIt,column,c);
                }
                ++inflIt;
            }
            else
            {   // truncate coefs from selected
                gsVector<index_t, DIM>  tenId = idSubDom.toMulti(*selIt);
                if (idSubDom.contains(tenId) )
                    coefs[tenId-idSubDom.lower()]=0;
                if(*inflIt==*selIt)
                    ++inflIt;
                ++selIt;
            }
        }
        while (inflIt!=inflEn) // complete saving
        {
            gsVector<index_t, DIM>  tenId = idSubDom.toMulti(*inflIt);
            if (idSubDom.contains(tenId) )
            {
                real_t c=coefs[tenId-idSubDom.lower()];
                if (c!=0) repr.add(rowShift+*inflIt,column,c);
            }
            ++inflIt;
        }
        while (selIt!=selEn ) // complete truncating
        {
            // truncate coefs
            gsVector<index_t, DIM>  tenId = idSubDom.toMulti(*selIt);
            if (idSubDom.contains(tenId) )
                coefs[tenId-idSubDom.lower()]=0;
            ++selIt;
        }
    }

    static void saveOnly (
            const indexSetT                   & influencing,
            const indexSetT                   & /*selected*/,
            const index_t                       rowShift,
            const index_t                       column,
            const gsTensorIndexSubDomain<DIM> & idSubDom,
            tensorCoefs                       & coefs,
            gsSparseEntries<real_t>           & repr
            )
    {
        typename indexSetT::const_iterator inflIt = std::lower_bound(influencing.begin(),influencing.end(),idSubDom.minFlat());
        typename indexSetT::const_iterator inflEn = std::upper_bound(inflIt,influencing.end(), idSubDom.endFlat());

        while (inflIt!=inflEn)
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




};

}

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlUtils.h>
#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo {

namespace internal
{

template <int d>
gsTHB<d> *getTHBFromXML(gsXmlNode * node)
{
    // TODO add attribute to differentiate between dyadic and provided
    // bases list

    // Initialize the HBSplineBasis
    std::istringstream str;
    gsXmlNode * tmp;

    // Insert all boxes
    unsigned c;
    typename gsTHB<d>::basisIdT maxLevel=0;
    std::vector<index_t> all_boxes;
    for (tmp = node->first_node("box");
         tmp; tmp = tmp->next_sibling("box"))
    {
        typename gsTHB<d>::basisIdT level(atoi( tmp->first_attribute("level")->value() ));
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
    std::vector<typename gsTHB<d>::basisPtr> bases;
    typename gsTensorBSplineBasis<d,real_t>::uPtr current(gsXml<gsTensorBSplineBasis<d,real_t> >::get(tmp));
    bases.push_back( gsBoxList::basisPtr(current->clone().release()) );
    for (typename gsTHB<d>::basisIdT lvl=1; lvl <= maxLevel;++lvl)
    {
        current->uniformRefine();
        bases.push_back(gsBoxList::basisPtr(current->clone().release()));
    }

    gsBoxList boxes(bases,bases[0]->domainDim(),all_boxes);
    return new gsTHB<d>(bases,boxes);
}


template <int DIM>
gsXmlNode *putTHBFromXML(const gsTHB<DIM> & obj,
                         gsXmlTree & data)
{
    // Add a new node (without data)
    gsXmlNode* tp_node = internal::makeNode("Basis" , data);

    std::string typeName;
    typeName+="THBSplineBasis";
    typeName+= DIM>1 ? to_string(DIM) :"";
    tp_node->append_attribute( makeAttribute("type", typeName, data) );

    // Write the component bases
    gsXmlNode * tmp = putTensorBasisToXml(*obj.getBasisOfLevel(0), data);
    tp_node->append_node(tmp);

    //Output boxes
    gsBoxList boxes = obj.getSelector().asBoxList();

    std::vector<index_t> boxEntries;
    gsMatrix<index_t> boxInd;
    for (size_t b=0; b< boxes.size(); ++b)
    {
        if (boxes.basisId(b)==0)
            continue;
        gsBoxList::toRefineElementFormat(obj.getBasisOfLevel(boxes.basisId(b)), boxes.basisId(b), boxes.box(b),boxEntries);
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
class gsXml< gsTHB<d> >
{
    typedef gsTHB<d> Object;
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(Object)
    static std::string tag  () { return "Basis"; }
    static std::string type () { return "THBSplineBasis"+ (d>1 ? to_string(d):""); }

    static gsTHB<d> * get (gsXmlNode * node)
    {
        GISMO_ASSERT(
                    ( !strcmp( node->name(),"Basis") ) && ( !strcmp(node->first_attribute("type")->value(), internal::gsXml<gsTHB<d> >::type().c_str() ) ),
                    "Something is wrong with the XML data: There should be a node with a "<<
                    internal::gsXml<Object>::type().c_str()<<" Basis.");
        return getTHBFromXML<d>(node);
    }

    static gsXmlNode * put (const gsTHB<d> & obj,
                            gsXmlTree & data )
    {
        return putTHBFromXML<d>(obj,data);
    }
};


} // internal




} // namespace gismo
