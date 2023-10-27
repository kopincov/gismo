/** @file gsSpaceRefiner.h

    @brief DESCRIPTION

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsRecipeAssemblerAdaptive/gsErrorEstimator.h>
#include <math.h>

namespace gismo {

class gsSpaceRefiner
{
public:
    gsSpaceRefiner()
    {}
    virtual ~gsSpaceRefiner()
    {}

    virtual std::vector<gsPhysicalSpace*> getSpaces () const = 0;
    virtual void updateSpaces (const gsMatrix<real_t>& marked) = 0;
};



class gsHTensorBasisRefiner : public gsSpaceRefiner
{
public:
    gsHTensorBasisRefiner()
    {}
    virtual ~gsHTensorBasisRefiner()
    {}

    virtual std::vector<gsPhysicalSpace*> getSpaces () const = 0;
    virtual void updateSpaces (const gsMatrix<real_t>& marked) = 0;
protected:
    template<unsigned d>
    static void getBoxPerPoint    (const gsHTensorBasis<d> &hTensorBasis, gsVector<real_t>& point, std::vector<index_t> &box)
    {
        const size_t dim = hTensorBasis.domainDim();
        const size_t begOffset = +1;
        const size_t endOffset = begOffset+dim;

        const unsigned lev=hTensorBasis.getLevelAtPoint(point);

        box.resize(2*dim+1);
        box[0]=(lev);
        for(size_t comp=0; comp<dim;++comp)
        {
            const long max = (hTensorBasis.getBases()[lev]->knots(comp).uSize()-1);
            const long idx = hTensorBasis.getBases()[lev]->knots(comp).uFind(point(comp,0)).uIndex();

            if (idx!=max)
            {
                box[begOffset+comp]=idx;
                box[endOffset+comp]=idx+1;
            }
            else if(idx!=0)
            {
                box[begOffset+comp]=idx-1;
                box[endOffset+comp]=idx;
            }
            else
                GISMO_ERROR("empty mesh");
        }
    }

    template< unsigned d>
    static void getBoxPerElement  (const gsHTensorBasis<d> &basis, index_t element, std::vector<unsigned> &box)
    {
        typename gsHTensorBasis<d>::domainIter iter=basis.makeDomainIterator();
        for (index_t eInd=0; eInd<element; ++eInd)
            ++iter;
        getBoxPerPoint(basis,iter.centerPoint(),box);
    }

    template< unsigned d>
    static void getBoxPerElement  (const gsHTensorBasis<d> &basis, const gsDomainIterator<real_t> &element, std::vector<unsigned> &box)
    {
        getBoxPerPoint(basis,element.centerPoint(),box);
    }

    template< unsigned d>
    static void getBoxPerFunction (const gsHTensorBasis<d> &hTensorBasis, index_t basisFuncIndex, std::vector<unsigned> &box)
    {
        const size_t dim = hTensorBasis.domainDim();
        const size_t begOffset = +1;
        const size_t endOffset = begOffset+dim;
        const unsigned lev=hTensorBasis.getLevelOf(basisFuncIndex);
        gsMatrix<real_t> supp = hTensorBasis.function(basisFuncIndex).support();

        box.resize(2*dim+1);
        box[0]=(lev);

        for(size_t comp=0; comp<dim;++comp)
        {
            const long max = (hTensorBasis.getBases()[lev]->knots(comp).uSize()-1);
            const long lowleft = hTensorBasis.getBases()[lev]->knots(comp).Uniquefindspan(supp(comp,0));
            const long upright = hTensorBasis.getBases()[lev]->knots(comp).Uniquefindspan(supp(comp,1));

            if (upright!=max)
            {
                box[begOffset+comp]=lowleft;
                box[endOffset+comp]=upright+1;
            }
            else if(lowleft!=0)
            {
                box[begOffset+comp]=lowleft-1;
                box[endOffset+comp]=upright;
            }
            else if(lowleft!=upright)
            {
                box[begOffset+comp]=lowleft;
                box[endOffset+comp]=upright;
            }
            else
                GISMO_ERROR("empty mesh");
        }
    }

    template< short_t d>
    static void liftLevel (const gsHTensorBasis<d> &hTensorBasis,size_t lift,std::vector<index_t> &boxes)
    {
        const size_t dim = hTensorBasis.domainDim();
        const size_t boxSize = 2*dim+1;
        const size_t shift = 1 << lift;
        for(size_t b = 0;b<boxes.size()/boxSize;++b)
        {
            const size_t begOffset = b*boxSize+1;
            const size_t endOffset = begOffset+dim;
            boxes[b*boxSize] += lift;
            for(size_t comp=0; comp<dim;++comp)
            {
                // here we could maybe get the new values from hTensorBasis.getBases()[lev+1]
                // to be able to cope with non-dyadic refinement
                boxes[begOffset+comp]*=shift;
                boxes[endOffset+comp]*=shift;
            }
        }
    }

    template< short_t d>
    static void addExtension (const gsHTensorBasis<d> &hTensorBasis,const gsVector<unsigned>& extensions,std::vector<index_t> &boxes)
    {
        const size_t dim = extensions.size();
        const size_t boxSize = 2*dim+1;
        for(size_t b = 0;b<boxes.size()/boxSize;++b)
        {
            const size_t begOffset = b*boxSize+1;
            const size_t endOffset = begOffset+dim;
            const size_t lev = boxes[b*boxSize];
            for(size_t comp=0; comp<dim;++comp)
            {
                const long max = (hTensorBasis.getBases()[lev]->knots(comp).uSize()-1);
                const int beg = boxes[begOffset+comp]-extensions(comp);
                const long end = boxes[endOffset+comp]+extensions(comp);
                boxes[begOffset+comp]=beg>=0   ? beg : 0;
                boxes[endOffset+comp]=end<=max ? end : max;
            }
        }
    }
};



}

