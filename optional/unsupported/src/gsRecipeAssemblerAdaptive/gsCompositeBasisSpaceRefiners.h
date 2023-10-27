/** @file gsCompositeBasisSpaceRefiners.h

    @brief space refiners based on gsCompositeHBasis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/


#pragma once

#include <gsRecipeAssemblerAdaptive/gsSpaceRefiner.h>
#include <gsSmoothPatches/gsCompositeHBasis.h>
#include <gsSmoothPatches/gsCompositeAssemblerUtils.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsHBSplineBasis.h>

#include <gsUtils/gsExportMatrix.h>

namespace gismo {

class gsCompositeHBasisRefiner : public gsSpaceRefiner
{
protected:
    const gsMultiPatch<>     &m_domain;
    mutable gsCompositeHBasis<2,real_t> m_basis;
    int                       m_refExt;
    int                       m_smooth;
    bool                      m_truncate;
    bool                      m_debug;

public:
    gsCompositeHBasisRefiner(
            const gsMultiPatch<>           &domain,
            std::vector<gsBasis<>*>         bases,
            int                             refExt=-1,
            int                             smoothness=-1,
            bool                            truncate=false
            )
        : m_domain(domain),
        m_refExt(refExt),
        m_smooth(smoothness),
        m_truncate(truncate),
        m_debug(false)
    {
        initBasis(bases);
    }

    void initBasis(std::vector<gsBasis<>*> bases)
    {
        m_basis.~gsCompositeHBasis<2,real_t>();
        std::vector<gsHTensorBasis<2,real_t>*> temp;
        for (size_t i=0; i<bases.size(); ++i)
        {
            gsHTensorBasis<2,real_t>* cur=NULL;
            cur = dynamic_cast<gsHTensorBasis<2,real_t>*>(bases[i]);
            if (cur)
                temp.push_back(cur->clone().release());
            else
            {
                cur= m_truncate ? static_cast<gsHTensorBasis<2,real_t>*>(new gsTHBSplineBasis<2,real_t>(*bases[i]))
                    : static_cast<gsHTensorBasis<2,real_t>*>(new gsHBSplineBasis<2,real_t>(*bases[i]));
                if (cur)
                    temp.push_back(cur);
                else
                    GISMO_ERROR("no chance");
            }
        }
        new (&m_basis) gsCompositeHBasis<2,real_t>(temp,m_domain,m_smooth);
        if(m_refExt==-1)
            m_refExt=(m_basis.maxDegree()+1)/2;
        freeAll(temp);
    }

    void setDebug(bool flag)
    {
        m_debug=flag;
    }

    virtual std::vector<gsPhysicalSpace*> getSpaces()const
    {
        std::vector<gsPhysicalSpace*> result;
        result.push_back(new gsPhysicalSpaceScalar (m_basis.getBases(),m_domain,INVERSE_COMPOSITION,m_basis.getMapper()));
        return result;
    }

    virtual void updateSpaces (const gsMatrix<real_t>& marked)
    {
        gsMatrix<real_t> markedMax = marked.colwise().maxCoeff();

        // numMarked: Number of marked cells on current patch, also currently marked cell
        // poffset  : offset index for the first element on a patch
        // globalCount: counter for the current global element index
        index_t globalCount = 0;

        // refBoxes: contains marked boxes on a given patch
        std::vector<index_t> patchBoxes;
        for (size_t patch=0; patch < m_basis.nPatches(); ++patch )// for all patches
        {
            patchBoxes.clear();
            // for all elements in patch pn
            gsBasis<real_t>::domainIter domIt = m_basis.basis(patch).makeDomainIterator();
            for (; domIt->good(); domIt->next(),++globalCount)
                if( 0!= markedMax(0,globalCount) )// refine this element ?
                    addBox(patch,domIt,patchBoxes);
            // Refine all of the found refBoxes in this patchx
            m_basis.refineElements( patch, patchBoxes, false );
        }
        m_basis.repairPatches();
        m_basis.updateTopol();
    }

    virtual void addBox(size_t patch,const gsBasis<real_t>::domainIter &domIt, std::vector<index_t> &boxes) const
    {
        const size_t dim = m_basis.dim();
        const size_t begOffset = boxes.size()+1;
        const size_t endOffset = begOffset+dim;


        const gsVector<real_t> &center=domIt->centerPoint();
        unsigned lev=m_basis.basis(patch).getLevelAtPoint(center);

        boxes.push_back(lev+1);
        boxes.resize(boxes.size()+2*dim);
        for(size_t comp=0; comp<dim;++comp)
        {
            const long max = 2*(m_basis.basis(patch).getBases()[lev]->knots(comp).uSize()-1);
            const long idx = 2*m_basis.basis(patch).getBases()[lev]->knots(comp).uFind(center(comp)).uIndex();
            const long deg = m_refExt; //m_basis.basis(patch).getBases()  [lev]->degree(comp);

            unsigned beg = idx>deg ? idx-deg : 0;;
            unsigned end = idx+deg+1<max ? idx+deg+1 : max;

            boxes[begOffset+comp]=beg;
            boxes[endOffset+comp]=end;
        }
        printBox(boxes,begOffset-1);
    }

    void printBox( const std::vector<index_t> &desc, size_t offset) const
    {
        if (m_debug)
        {
            const size_t dim = m_domain.dim();
            const size_t begOff = offset+1;
            const size_t endOff = offset+dim+1;

            std::cout<<"-> l:"<<desc[offset]<<"; beg : ";
            for (size_t i=0; i<dim;++i)
                std::cout<<desc[begOff+i]<<", ";
            std::cout<<"; end : ";
            for (size_t i=0; i<dim;++i)
                std::cout<<desc[endOff+i]<<", ";
            std::cout<<"\n";
        }
    }
};



class gsStokesRefiner : public gsCompositeHBasisRefiner
{
public:
    enum
    {
        none=0,
        TH,
        SG,
        RT // mapper not implemented yet
    };
protected:
    int                                                m_method;
    mutable std::vector<gsCompositeHBasis<2,real_t>*>  m_basisV;
    int                                                m_enforceGradedMesh;
    bool                                               m_useUniform;
public:
    gsStokesRefiner(
            const gsMultiPatch<>           &domain,
            const std::vector<gsBasis<>*>  &bases,
            int                             method=TH,
            int                             smoothness=0,
            bool                            truncate=false,
            int                             levelSep=-1
            )
        :
          gsCompositeHBasisRefiner(domain, bases, 0,smoothness, truncate),
          m_method(method),
          m_enforceGradedMesh(levelSep),
          m_useUniform(false)
    {
        if(m_enforceGradedMesh<0)
            m_enforceGradedMesh=m_basis.maxDegree();
        makeVelocityBasis ();
    }

    void makeVelocityBasis () const
    {
        const int dim=m_domain.dim();
        freeAll(std::unique(m_basisV.begin(),m_basisV.end()),m_basisV.end());
        switch (m_method)
        {
        case TH:
        {
            gsCompositeHBasis<2,real_t>* basisV=m_basis.clone().release();
            basisV->degreeElevate();
            /// TODO check line multiplicity in all levels
            m_basisV=std::vector<gsCompositeHBasis<2,real_t>*>(dim,basisV);
        }
            break;
        case SG:
        {
            gsCompositeHBasis<2,real_t>* basisV=m_basis.clone().release();
            basisV->degreeIncrease(1,-1,false);
            basisV->uniformRefine();
            /// TODO check line multiplicity in all levels
            m_basisV=std::vector<gsCompositeHBasis<2,real_t>*>(dim,basisV);
        }
            break;
        case RT:
        {
            /// TODO remove gluing as it is wrong
            m_basisV.resize(dim);
            for (int comp=0; comp<dim;++comp)
            {
                gsCompositeHBasis<2,real_t>* basisV=m_basis.clone().release();
                basisV->degreeIncrease(1,comp);
                m_basisV[comp]=basisV;
            }
        }
            break;
        default:
            GISMO_ERROR("method not implemented");
        }
    }

    ~gsStokesRefiner()
    {
        freeAll(std::unique(m_basisV.begin(),m_basisV.end()),m_basisV.end());
    }

    virtual std::vector<gsPhysicalSpace*> getSpaces()const
    {
        const int dim=m_domain.dim();

        ValueTransformationType transform = (m_method==RT ? DIV_CONFORMING : INVERSE_COMPOSITION);

        std::vector<gsPhysicalSpace*> result;
        std::vector<gsPhysicalSpaceScalar*> velV;
        for (int comp=0;comp<dim;++comp)
        {
            velV.push_back(new gsPhysicalSpaceScalar (m_basisV[comp]->getBases(),m_domain,transform,m_basisV[comp]->getMapper()));
        }
        result.push_back(new gsPhysicalSpaceVector(velV));
        result.push_back(new gsPhysicalSpaceScalar (m_basis.getBases(),m_domain,INVERSE_COMPOSITION,m_basis.getMapper()));
        return result;
    }

    virtual void updateSpaces (const gsMatrix<real_t>& marked)
    {
        if (m_useUniform)
            m_basis.uniformRefine();
        else
            gsCompositeHBasisRefiner::updateSpaces(marked);
        makeVelocityBasis();
    }

    virtual void addBox(size_t patch,const gsBasis<real_t>::domainIter &domIt, std::vector<index_t> &boxes) const
    {
        const size_t dim = m_basis.dim();
        const size_t begOffset = boxes.size()+1;
        const size_t endOffset = begOffset+dim;


        const gsVector<real_t> &center=domIt->centerPoint();
        unsigned lev=m_basis.basis(patch).getLevelAtPoint(center);

        boxes.push_back(lev+1);
        boxes.resize(boxes.size()+2*dim);
        for(size_t comp=0; comp<dim;++comp)
        {
            const long max = 2*(m_basis.basis(patch).getBases()[lev]->knots(comp).uSize()-1);
            const long idx = 2*m_basis.basis(patch).getBases()[lev]->knots(comp).uFind(center(comp)).uIndex();
            const long deg = m_refExt; //m_basis.basis(patch).getBases()  [lev]->degree(comp);

            unsigned beg = idx>deg ? idx-deg : 0;;
            unsigned end = idx+deg+1<max ? idx+deg+1 : max;

            boxes[begOffset+comp]=beg;
            boxes[endOffset+comp]=end;
        }
        printBox(boxes,begOffset-1);
        appendBoxesForMeshGrading(patch,boxes);
    }

    void appendBoxesForMeshGrading(size_t patch, std::vector<index_t> &boxes) const
    {
        const size_t dim=m_domain.dim();
        const size_t blockSize=2*dim+1;

        const size_t originalStart=boxes.size()-blockSize;

        size_t begOffset=originalStart+1;
        size_t endOffset=begOffset+dim;

        std::vector<index_t> beg(boxes.data()+begOffset, boxes.data()+begOffset+dim);
        std::vector<index_t> end(boxes.data()+endOffset, boxes.data()+endOffset+dim);

        for (index_t lev=boxes[originalStart]-1;lev>0;--lev)
        {
            size_t curOffset=boxes.size();
            begOffset=curOffset+1;
            endOffset=begOffset+dim;

            boxes.resize(curOffset+blockSize);
            boxes[curOffset]=lev;

            for ( size_t comp=0; comp<dim; ++comp)
            {
                const index_t max = m_basis.basis(patch).getBases()[lev]->knots(comp).uSize()-1;

                beg[comp]= beg[comp]/2>= static_cast<index_t>(m_enforceGradedMesh) ? beg[comp]/2-m_enforceGradedMesh : 0;
                end[comp]= (end[comp]+1)/2+static_cast<index_t>(m_enforceGradedMesh)>max ? max: (end[comp]+1)/2+m_enforceGradedMesh ;

                boxes[begOffset+comp]=beg[comp];
                boxes[endOffset+comp]=end[comp];
            }
            printBox(boxes,curOffset);
        }
    }

    void setC0 (patchCorner pc) const {
        m_basis.setC0(pc);
    }

    void update () const
    {
        m_basis.updateTopol();
    }

    const gsBasis<real_t> * pressureBasis(size_t patch)
    {
        return &m_basis.basis(patch);
    }

    const gsBasis<real_t> * velocityBasis(size_t patch)
    {
        return &m_basisV[0]->basis(patch);
    }

    void setUniform (bool flag)
    {
        m_useUniform=flag;
    }

};



} // namespace gismo

