/** @file gsDecoupledBasis.h

    @brief Implementation of iteratively decoupled basis.
    This is a modification of TDHB (cf. Mokris, Juettler:
    TDHB-splines (...), CAGD, 2014), where decoupling is
    repeated in every level instead of truncation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris, A. Bressan
**/

#pragma once

#include <gsRemappedBasis/gsAdaptiveSplineCommon.h>

namespace gismo
{

/** Helper class for holding the decoupled functions (i.e., if a
 * B-spline (or a decoupled function) is decoupled, these will be the
 * results). */    
class decoupledFunction
{
public:
    decoupledFunction()
    {}
    decoupledFunction( index_t origLvl,
                       index_t origFun,
                       gsDomainMap support,
                       gsSparseVector<real_t> coeffs,
                       bool isSelected )
        : m_origLvl(origLvl),
          m_origFun(origFun),
          m_support(support),
          m_coeffs(coeffs),
          m_isSelected(isSelected)
    {
        // Nothing else to do?
    }

    // basisIdT for the levels
    operator index_t()
    {
        return m_origFun;
    }

    bool operator<(index_t i) const {return m_origFun < i;}

public: // data
    index_t                 m_origLvl;
    index_t                 m_origFun;
    gsDomainMap             m_support;
    gsSparseVector<real_t>  m_coeffs;
    bool                    m_isSelected;
};

/** Implementation of (repeatedly) decoupled hierarchical splines as
 * described in \cite bm2016, Section 4.3.*/
template <int DIM>
class gsDecoupledBasis : public gsAdaptiveSplineCommon<DIM>
{
    // Lovely typedefs and usings.
    using  gsAdaptiveSplineCommon<DIM>::m_basis;
    using  gsAdaptiveSplineCommon<DIM>::m_info;
    using  gsAdaptiveSplineCommon<DIM>::m_sele;
    using  gsAdaptiveSplineCommon<DIM>::m_repr;
    using  gsAdaptiveSplineCommon<DIM>::m_shift;

    typedef typename gsAdaptiveSplineCommon<DIM>::collectInfluencing collectInfluencing;
    typedef typename gsAdaptiveSplineCommon<DIM>::indexSetT          indexSetT;
    typedef typename gsAdaptiveSplineCommon<DIM>::indexDomT          indexDomT;
    typedef typename gsAdaptiveSplineCommon<DIM>::knotRange          knotRange;
    typedef typename gsAdaptiveSplineCommon<DIM>::tensorCoefs        tensorCoefs;
    typedef typename gsAdaptiveSplineCommon<DIM>::makeBasisPerLeaf   makeBasisPerLeaf;
    typedef typename gsAdaptiveSplineCommon<DIM>::gsTBS42            gsTBS42;
    typedef typename gsAdaptiveSplineCommon<DIM>::basisPtr basisPtr;

    typedef typename gsAdaptiveSplineCommon<DIM>::getMinAndMaxLevel getMinAndMaxLevel;

public: // constructors

    gsDecoupledBasis( std::vector< basisPtr > basis,
                      const gsBoxList &boxes )
    {
        m_info.first = DIM;
        m_info.second= 1;

        m_defBasis=basis;

        GISMO_ASSERT(boxes.maxId()<basis.size(), "You pushed a box with non-existent basisId.");

        gsBoxList mBoxList=boxes;
        this->sanitizeBoxList(m_defBasis,mBoxList,-1);
        //  checkDimAndInitShifts();

        m_domain.initFromBoxesMax(mBoxList,asBasis<DIM>(m_defBasis[0])->support());
        makeRepresentation();
    }

    void printSupport(const decoupledFunction& fun,
                      real_t& x,
                      real_t& y,
                      gsBoxList::basisIdT& lvl)
    {
        real_t left, right, bottom, top;
        gsBoxList support = fun.m_support.asBoxList();
        gsAsMatrix<real_t> currBox = support.box(0);
        left = 1;
        right = 0;
        bottom = 1;
        top = 0;
        //gsInfo << left << ", " << right << "\n" << bottom << ", " << top << std::endl ;
        //gsInfo << currBox << std::endl << std::endl;

        for(size_t p = 0; p < support.size(); ++p)
        {
            currBox = support.box(p);
            if(support.basisId(p)==1)
            {
                left = std::min<real_t>(left,currBox(0,0));
                right = std::max<real_t>(right,currBox(0,1));
                bottom = std::min<real_t>(bottom,currBox(1,0));
                top = std::max<real_t>(top,currBox(1,1));
            }
        }

        x = (left+right)/2;
        y = (bottom+top)/2;
        lvl = fun.m_origLvl;
    }

    void printComprehensive()
    {
        real_t x = 0;
        real_t y = 0;
        gsBoxList::basisIdT lvl = 17;
        std::ofstream fout;
        fout.open("decoupled.tex");
        if( fout.is_open() )
        {

            fout<<"\\documentclass{standalone}\n"
                  "\\usepackage{tikz}\n"
                  "\\begin{document}\n"
                  "\\begin{tikzpicture}[scale=20]\n";
            for(std::vector<decoupledFunction>::iterator it = m_comprehensiveVec.begin();
                it != m_comprehensiveVec.end();
                ++it)
                if( it->m_isSelected)
                {
                    printSupport(*it,x,y,lvl);
                    fout
                            << "\\node at (" << x << "," << y << ") {"
                            << lvl << "};\n";
                }
            fout << "\\end{tikzpicture}\n\n"
                    "\\end{document}\n";
            fout.close();
            int pdfErr=system("pdflatex decoupled.tex >/dev/null");
            pdfErr=system("okular decoupled.pdf >/dev/null");
            gsInfo << pdfErr;
        }
        else
            gsWarn << "Error opening file, output not saved!";

    }

    std::vector<size_t> levelSizes()
    {
        std::vector<size_t> result(m_defBasis.size(),0);
        for(typename std::vector<decoupledFunction>::iterator it = m_comprehensiveVec.begin();
            it != m_comprehensiveVec.end();
            ++it)
            if(it->m_isSelected)
                ++result[it->m_origLvl];
        return result;
    }

private:

    /// Creates the representation matrix.
    void makeRepresentation()
    {
        makeComprehensiveVec();      // compute the list of decoupled functions
        makeDomainAndBasis();
        makeRepresentationMatrix();  // computes a representation suitable for evaluation
    }

    void makeDomainAndBasis()
    {
        m_sele.patch(0)=m_domain.patch(0);
        makeBasisPerLeaf leafMaker(m_defBasis,m_basis);
        m_sele.patch(0).apply(leafMaker);
        // After setting m_basis we can finally init the shifts:
        gsRemappedBasis::checkDimAndInitShifts();
    }

    class LeavesIdCollector : public gsDomainMap::actionBase<DIM,2>
    {
    public:
        LeavesIdCollector(std::set<gsDomainMap::NodeId> &collectedLeaves)
            : m_collectedLeaves(collectedLeaves)
        {}
        std::set<gsDomainMap::NodeId> &m_collectedLeaves;
        typedef typename gsDomainMap::actionBase<DIM,2>::BoxT BoxT;
        bool enter(const gsDomainMap*  map,
                   const BoxT        & /*box*/,
                   gsDomainMap::NodeId node)
        {
            if( (*map)[node].isFork() )
                return true;
            m_collectedLeaves.insert(node);
            return true;
        }
    };


    // TODO: enter and enterLeaf
    class Gardener : public gsDomainMap::actionBase<DIM,2>
    {
    public:
        typedef typename gsDomainMap::actionBase<DIM,2>::BoxT BoxT;

        std::set<gsDomainMap::NodeId> collectedLeaves;
        gsDomainMap &m_seleP;
        Gardener(gsDomainMap &seleP )
            :m_seleP(seleP)
        {}
        bool enter(const gsDomainMap *map,
                   const BoxT &box,
                   gsDomainMap::NodeId node)
        {
            if( (*map)[node].isFork() || (*map)[node].data.space==0 )
                return true;

            LeavesIdCollector Michael(collectedLeaves);
            m_seleP.apply(box,Michael);
            return true;
        }
    };

    gsBoxList leaveSetToBoxList( const std::set<gsDomainMap::NodeId> &leaves )
    {
        gsBoxList result(DIM);
        for( typename std::set<gsDomainMap::NodeId>::iterator it = leaves.begin();
             it != leaves.end();
             ++it )
            result.append(m_sele.patch(0).getBoundingBox(*it),
                          m_sele.patch(0)[*it].data.space);
        return result;
    }

    void computeReprCoefs(size_t bId, const index_t &count, gsSparseEntries<real_t> &entries)
    {
        // GISMO_ASSERT( season == autumn, "No leaves to be found.\n");
        Gardener Thomas(m_sele.patch(0)); // Thomas collects leaves
        m_comprehensiveVec[bId].m_support.apply(Thomas);

        gsBoxList leavesInSupport = leaveSetToBoxList(Thomas.collectedLeaves);
        std::set<gsBoxList::basisIdT> leafIndices;
        for(size_t b = 0; b < leavesInSupport.size(); ++b)
            leafIndices.insert(leavesInSupport.basisId(b));

        gsBoxList::basisIdT currLvl =  m_comprehensiveVec[bId].m_origLvl;

        gsSparseVector<real_t> localCoefs(m_levelEnd[currLvl]-m_levelStart[currLvl]);
        localCoefs(bId-m_levelStart[currLvl])=1;


        representOnLeaves(currLvl,localCoefs,leafIndices,count,entries);
    }

    void representOnLeaves(const gsBoxList::basisIdT lvl, const gsSparseVector<> &coefs, 
                           std::set<gsBoxList::basisIdT> &leaves, const index_t &count, 
                           gsSparseEntries<real_t> &entries)
    {
        bool shouldContinue=false;
        for(typename std::set<gsBoxList::basisIdT>::iterator jt=leaves.begin(); jt != leaves.end();++jt)
        {
            gsTBS42     & _basis   = this->asTBS42(m_basis[*jt]);
            gsDomainMap   _leafSup = _basis.leafDomain();
            if(_basis.origLvl()==lvl)
            {
                gsTensorIndexSubDomain<DIM> dom=_basis.subdomain();
                gsSparseVector<>::InnerIterator it(coefs,0);
                for (;it;++it)
                {
                    const index_t origFun = m_comprehensiveVec[m_levelStart[lvl]+it.row()].m_origFun;
                    const index_t relFunId = dom.toPos(origFun);
                    gsDomainMap intersection=polytopeIntersect(m_comprehensiveVec[m_levelStart[lvl]+it.row()].m_support, _leafSup );
                    if(dom.contains(origFun) && !isEmpty(intersection))
                        entries.add(m_shift[*jt]+relFunId,
                                count,
                                it.value());
                }
            }
            else if( _basis.origLvl()>lvl)
                shouldContinue=true;
        }

        if (shouldContinue)
        {
            gsSparseVector<> newCoefs(m_levelEnd[lvl+1]-m_levelStart[lvl+1]);
            gsSparseVector<>::InnerIterator it(coefs,0);
            for (;it;++it)
            {
                newCoefs+=it.value()*m_comprehensiveVec[m_levelStart[lvl]+it.row()].m_coeffs;
            }

            GISMO_ASSERT(lvl+1 < static_cast<gsBoxList::basisIdT>(m_defBasis.size()),"What is the level of the remaining leaves?\n");

            representOnLeaves(lvl+1, newCoefs, leaves,count,entries);
        }
    }

    void makeRepresentationMatrix()
    {
        gsSparseEntries<real_t> entries;
        index_t           count=0;//number of decoupled functions written to the matrix so far.
        gsBoxList::basisIdT numLvl = static_cast<gsBoxList::basisIdT>(m_defBasis.size());
        for(gsBoxList::basisIdT lvl = 0; lvl < numLvl; ++lvl)
            for(size_t bId = m_levelStart[lvl]; bId < m_levelEnd[lvl]; ++bId)
                if(isSelected(bId))
                {
                    computeReprCoefs(bId,count,entries);
                    ++count;
                }
        m_repr.asMatrix().resize(m_shift.back(),count);
        m_repr.asMatrix().setFromTriplets(entries.begin(), entries.end());
    }

    void makeComprehensiveVec()
    {
        m_levelEnd.resize(m_defBasis.size()+1);
        m_levelStart.resize(m_defBasis.size()+1);
        m_levelEnd.back()=0;
        m_levelStart.back()=0;// Mark that n+1st level is empty.

        const gsBoxList::basisIdT numLvl=static_cast<gsBoxList::basisIdT>(m_defBasis.size());
        typename gsRemTypes<DIM>::tensorBasisT* currBasis;

        for( int lvl = numLvl-1; lvl >= 0; --lvl )
        {
            m_levelStart[lvl]=m_comprehensiveVec.size();
            currBasis = asBasis<DIM>(m_defBasis[lvl]);
            for( index_t i = 0; i < currBasis->size(); ++i )
                pushIntoComprehensiveVec(lvl, i);
            m_levelEnd[lvl]=m_comprehensiveVec.size();
        }
    }

    /// Find what it decouples into and push it to the vector.
    void pushIntoComprehensiveVec(gsBoxList::basisIdT lvl, index_t id)
    {
        const index_t maxLvl = m_defBasis.size()-1;
        if(static_cast<index_t>(lvl)==maxLvl)
            pushSimple(lvl, id);
        else
            pushDecoupled(lvl, id);
    }

    /// Represents a B-spline with id \a id from level \a lvl
    /// in terms of B-splines of level \a lvl + 1
    /// and provides the corresponding subdomain iterator.
    void representInFinerLvl( gsBoxList::basisIdT lvl, index_t id,
                              tensorCoefs& coefs,
                              gsTensorIndexSubDomain<DIM>& indexSubDom)
    {
        knotRange oldK[DIM];
        knotRange newK[DIM];
        gsMatrix<real_t,DIM,2>  support; // the support of the function is the cartesian product of the intervals described by the rows
        indexDomT               mult;    // mult(dir,i) is the multiplicity of support(dir,i) in the knot vector of the function

        this->getTensorFunctionInfo   (m_defBasis[lvl],id,support,mult,oldK);
        this->getTensorLevelLocalView (m_defBasis[lvl+1],support,mult,newK,indexSubDom);
        coefs.initOne(indexSubDom.sizeM());
        std::vector<real_t> tmpStorage(indexSubDom.sizeM().bottomRows(DIM-1).maxCoeff());
        this->updateCoefs(oldK, newK, coefs, tmpStorage);
    }


    /// Represents a B-spline with id \a id in terms of (decoupled) functions of level \a lvl + 1.
    gsSparseVector<real_t> finerRepre( gsBoxList::basisIdT lvl, index_t id )
    {
        tensorCoefs coefs;
        gsTensorIndexSubDomain<DIM> indexSubDom;
        representInFinerLvl( lvl, id, coefs, indexSubDom );

        typename std::vector< decoupledFunction >::iterator curr, lvlEnd;
        curr = m_comprehensiveVec.begin() + m_levelStart[lvl+1];
        lvlEnd = m_comprehensiveVec.begin() + m_levelEnd[lvl+1];

        gsSparseVector<real_t> functionsToGroup(m_levelEnd[lvl+1]-m_levelStart[lvl+1]);

        for(typename gsTensorIndexSubDomain<DIM>::flat_iterator it = indexSubDom.flat_begin();
            it != indexSubDom.flat_end();
            ++it )
        {
            real_t currCoef = coefs[it.multiIndex()-indexSubDom.lower()];
            if (curr->m_origFun < *it )
                ++curr;
            if (curr != lvlEnd && curr->m_origFun < *it )
                curr = std::lower_bound( curr, lvlEnd, *it );

            if( curr == lvlEnd )
                GISMO_ERROR( "No function found, weird.");
            if( *curr > *it ) // The function is not in the comprVec.
                GISMO_ERROR( "Function inside support not Refinement outside support, weird.");
            GISMO_ASSERT( curr->m_origFun == *it, "This should be guaranteed by the previous two checks." );

            do
            {
                functionsToGroup(curr-m_comprehensiveVec.begin()-m_levelStart[lvl+1])= currCoef;
                ++curr;
            }
            while (curr->m_origFun == *it);
        }

        return functionsToGroup;
    }

    bool isSelected(size_t bId)
    {
        return m_comprehensiveVec[bId].m_isSelected;
    }

    bool isSelected( gsBoxList::basisIdT lvl, gsDomainMap support )
    {
        gsBoxList supportBoxes(DIM);
        supportBoxes=support.asBoxList();
        getMinAndMaxLevel presentLvl(lvl);
        for (size_t b=0;b<supportBoxes.size();++b)
        {
            if(supportBoxes.basisId(b)!=0)
                m_domain.patch(0).apply(supportBoxes.box(b),presentLvl);
        }
        return presentLvl.minLvl == lvl;
    }

    struct getLeavesOfLevel : public gsDomainMap::actionBase<DIM,2>
    {
        typedef typename gsDomainMap::actionBase<DIM,2>::BoxT BoxT;
        gsBoxList &m_boxes;
        gsBoxList::basisIdT m_level;
        getLeavesOfLevel(gsBoxList &boxes,gsBoxList::basisIdT level)
            : m_boxes(boxes), m_level(level)
        {

        }

        bool enter(const gsDomainMap *map,
                   const BoxT &box,
                   gsDomainMap::NodeId node)
        {
            if( (*map)[node].isFork() )
                return true;
            if( (*map)[node].data.space <= m_level ) // <= instead of == anywhere else?
                m_boxes.append(box,m_level);
            return true;
        }
    };

    void pushDecoupled( gsBoxList::basisIdT lvl, index_t id )
    {
        const gsSparseVector<real_t> expansionInFinerDec = finerRepre( lvl, id );

        knotRange oldK[DIM];
        gsMatrix<real_t,DIM,2> wholeSupport;
        indexDomT mult;
        this->getTensorFunctionInfo(m_defBasis[lvl],id,wholeSupport,mult,oldK);

        gsBoxList leavesToProcess(DIM);
        getLeavesOfLevel searcher(leavesToProcess,lvl);
        m_domain.patch(0).apply(wholeSupport,searcher);

        if ( leavesToProcess.size()==0 )
            pushSimple(lvl,id);

        while( leavesToProcess.size()>0 )
        {
            m_comprehensiveVec.push_back(getFirstDecoupled(lvl,
                                                           id,
                                                           leavesToProcess,
                                                           expansionInFinerDec));
            m_comprehensiveVec.back().m_isSelected=isSelected(lvl,m_comprehensiveVec.back().m_support);
        }
    }

    static bool compareId(const decoupledFunction& left, const index_t &value)
    {
        return left.m_origFun < value;
    }

    decoupledFunction getFirstDecoupled(gsBoxList::basisIdT origLvl,
                                        index_t origId,
                                        gsBoxList& leavesToProcess,
                                        const gsSparseVector<real_t>& expansion )
    {
        gsDomainMap curLeaf;
        gsDomainMap intersection;
        size_t nextLeaf=0;

        decoupledFunction result;
        result.m_origLvl = origLvl;
        result.m_origFun = origId;
        // Set the support to be empty with the proper bounding box:
        result.m_support.initFromBoxesMax(gsBoxList(DIM),m_domain.patch(0).getBoundingBox());
        result.m_coeffs.resize(m_levelEnd[origLvl+1]-m_levelStart[origLvl+1]);
        while(nextLeaf<leavesToProcess.size())
        {
            // init iteration
            curLeaf=domainMapFromBox(leavesToProcess.box(nextLeaf));
            // grow decoupled function
            for(typename gsSparseVector<real_t>::InnerIterator it( expansion); it;++it)
            {
                const decoupledFunction& func = m_comprehensiveVec[it.row()+m_levelStart[origLvl+1]];
                intersection = polytopeIntersect(curLeaf,func.m_support);
                if (!isEmpty(intersection))
                {
                    result.m_support=polytopeUnion(result.m_support,func.m_support);
                    result.m_coeffs(it.row())=it.value();
                }
            }
            leavesToProcess.remove(nextLeaf);
            // decide next iteration or set nextLeaf to exit the loop
            nextLeaf=0;
            while (nextLeaf<leavesToProcess.size() )
            {
                intersection = polytopeIntersect(result.m_support, domainMapFromBox(leavesToProcess.box(nextLeaf)));
                if( !isEmpty(intersection) )
                    break;
                ++nextLeaf;
            }
        }
        return result;
    }

    bool isEmpty(const gsDomainMap& map) const
    {
        getMinAndMaxLevel presentLvl(0);
        map.apply(presentLvl);
        GISMO_ASSERT(presentLvl.maxLvl<=1,"I want a boolean map.");
        return presentLvl.maxLvl == 0;
    }

    // TODO: Consider using the  gsDomainMap constructor.
    inline gsDomainMap domainMapFromBox(const gsMatrix<real_t,DIM,2> &box)
    {
        gsBoxList boxes(DIM);
        boxes.append(box,1);
        gsDomainMap result;
        result.initFromBoxesMax(boxes,m_domain.patch(0).getBoundingBox());
        return result;
    }

    void pushSimple( gsBoxList::basisIdT lvl, index_t id )
    {
        gsMatrix<real_t> support = asBasis<DIM>(m_defBasis[lvl])->support(id);
        m_comprehensiveVec.push_back( decoupledFunction(lvl, id,
                                                        domainMapFromBox(support),
                                                        gsSparseVector<real_t>(comprLvlSize(lvl+1)),
                                                        isSelected(lvl,support) ));
    }

    inline index_t comprLvlSize(const gsBoxList::basisIdT lvl) const
    {
        return m_levelEnd[lvl]-m_levelStart[lvl];
    }

public:
    void debugWithTeX()
    {
        const int maxLvl = m_defBasis.size()-1;
        std::ofstream out;
        out.open("debugWithTex.tex");

        if(!out.is_open())
            GISMO_ERROR("I want to write!");
        out << std::setprecision(2);
        out << std::fixed;
        out<<"\\documentclass[a3paper,landscape]{article}\n"
             "\\usepackage[margin=2cm]{geometry}\n"
             "\\usepackage{tikz}\n"
             "\\begin{document}\n";

        out<<"\\section{Remapped basis info}\n"
             "Shifts: "<< gsAsMatrix<index_t>(m_shift)<< "\n\n";

        int selId=0;
        for(int lvl = 0; lvl<=maxLvl; ++lvl)
        {
            out<<"\\clearpage\\section*{Lvl "<<lvl<<"}\n";
            for(size_t f = m_levelStart[lvl]; f < m_levelEnd[lvl]; ++f)
            {
                const decoupledFunction& func = m_comprehensiveVec[f];
                out << "\\noindent ";
                if (func.m_isSelected)
                    out<< "S: "<<selId++;
                else
                    out<< "S: --";
                out << ", Id: "<< f-m_levelStart[lvl] << ", origId: " << func.m_origFun << ", origLvl: " << func.m_origLvl <<  (func.m_origLvl!=lvl ? "{\\Large WRONG}\n\n": "\n\n");
                out << "\\noindent{\\footnotesize";
                gsSparseVector<>::InnerIterator it(func.m_coeffs);
                if (it)
                {
                    out << "\\begin{tabular}{";
                    for(int r=0;r<func.m_coeffs.nonZeros();++r)
                        out << " c ";
                    out << "}\n";
                    out << it.row();
                    while (++it)
                    {
                        out  << " & " << it.row();
                    }
                    gsSparseVector<>::InnerIterator itit(func.m_coeffs);
                    out <<  "\\\\ \n" << itit.value();
                    while (++itit)
                    {
                        out  << " & " << itit.value();
                    }
                    out << "\\end{tabular}";
                }

                out << "}\\hfill\\begin{tikzpicture}[scale=1.4]\n";
                out << func.m_support.asBoxList();
                out << "\\end{tikzpicture}\n\n";
                out << "\\noindent\\rule{\\textwidth}{2pt}";
                out << "\n\n";
            }
        }
        out << "\n\\end{document}";

        m_domain.exportToTex("defDomain");
    }

    bool testReprMatrix()
    {
        std::vector<index_t> rowsThatSumUpToMoreThanOne;
        std::vector<index_t> selectedInvolved;

        const gsMatrix<real_t> rowSums = m_repr.asMatrix().toDense().rowwise().sum();
        for(index_t r =0; r < rowSums.rows(); ++r )
            if(math::abs(rowSums(r,0)-1)>1e-13)
                rowsThatSumUpToMoreThanOne.push_back(r);

        m_repr.sourceToTarget(rowsThatSumUpToMoreThanOne,selectedInvolved);

        if (rowsThatSumUpToMoreThanOne.size())
        {
            gsInfo << "rows: " << gsAsConstMatrix<index_t>(rowsThatSumUpToMoreThanOne) << std::endl;
            gsInfo << "sele: " << gsAsConstMatrix<index_t>(selectedInvolved) << std::endl;
            gsMatrix<real_t> map;
            m_repr.getLocalMap(rowsThatSumUpToMoreThanOne,selectedInvolved,map);
            gsInfo << "map:  " << map << std::endl;
        }
        return rowsThatSumUpToMoreThanOne.size() == 0;
    }

private: // data

    gsSelector                       m_domain;
    std::vector<basisPtr>            m_defBasis;
    std::vector< size_t >            m_levelStart;
    std::vector< size_t >            m_levelEnd;
    std::vector< decoupledFunction > m_comprehensiveVec;
};

} // namespace gismo
