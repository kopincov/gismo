/** @file gsHFittingLvlConstrained.h

    @brief Adaptive fitting using hierarchical splines, where
    boxes can be specified, where different maxLvls of refinement are prescribed

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsModeling/gsFitting.h>
#include <gsHSplines/gsHTensorBasis.h>

namespace gismo {

/**
    \brief
    This class applies hierarchical fitting of parametrized point clouds.

    \tparam T coefficient type

    \ingroup HSplines
*/

template <short_t d, class T>
class gsHFittingLvlConstrained : public gsHFitting<d,T>
{
public:
    typedef typename gsBSplineTraits<d,T>::Basis tensorBasis;

public:
    /// Default constructor
    gsHFittingLvlConstrained();

    /**
        \brief
        Main constructor of the fitting class

        \param param_values a matrix containing the parameter values
        that parametrize the \a points

        \param points The points to be fitted

        \param basis  Hiearchical basis to use for fitting

        \param refin Percentage of errors to refine (if this strategy is chosen)

        \param extension Extension to apply to marked cells

        \param lambda Smoothing parameter
    */
    gsHFittingLvlConstrained(gsMatrix<T> const & param_values,
               gsMatrix<T> const & points,
               gsHTensorBasis<d,T> & basis,
               T refin, const std::vector<unsigned> & extension,
               T lambda,
               const std::vector<index_t>& constrainedBoxes)
    : gsHFitting<d,T>(param_values, points, basis,refin,extension,lambda)
    {
        std::vector<index_t> box;
        std::vector<index_t>::const_iterator first,last;
        short_t boxSize=2*d+1;
        for(unsigned i = 0;i<constrainedBoxes.size();i=i+boxSize)
        {
            int lvl = constrainedBoxes[i];
            first=constrainedBoxes.begin()+i;
            last=constrainedBoxes.begin()+(i+boxSize);
            box=std::vector<index_t>(first,last);
            constrainedBox constrainedBox;
            constrainedBox.m_box=box;
            constrainedBox.m_maxLvl=lvl;
            m_constrainedBoxes.push_back(constrainedBox);
        }
    }

protected:
    /// Appends a box around parameter to the boxes only if the box is not
    /// already in boxes
    virtual void appendBox(std::vector<index_t>& boxes,
                   std::vector<index_t>& cells,
                   const gsVector<T>& parameter);

private:
    bool isBoxViolatingConstrainedBoxes(const std::vector<index_t>& a_box)
    {
        for(unsigned i =0;i<m_constrainedBoxes.size();i++)
        {
            if(m_constrainedBoxes[i].m_maxLvl>= a_box[0])
            {
                continue;
            }
            if(boxesOverlapping(a_box,m_constrainedBoxes[i].m_box))
            {
                return true;
            }
        }
        return false;
    }

    bool boxesOverlapping(std::vector<index_t> box,std::vector<index_t> constrainedBox)
    {
        unsigned highestLvl = box[0]>constrainedBox[0] ? box[0] : constrainedBox[0];
        bringBoxToLvl(constrainedBox,highestLvl);
        bringBoxToLvl(box,highestLvl);

        bool intersect = true;
        for(short_t i = 0;i<d;++i)
            intersect = intersect && box[i+1+d]>constrainedBox[i+1] && box[i+1]<constrainedBox[i+1+d];
        return intersect;
    }

    bool smaller(const gsVector<unsigned,d>& vec1, const gsVector<unsigned,d>& vec2)
    {
        for(unsigned i = 0;i<d;++i)
            if(vec1(i)>=vec2(i))
                return false;
        return true;
    }

    bool smallerEqual(const gsVector<unsigned,d>& vec1, const gsVector<unsigned,d>& vec2)
    {
        for(unsigned i = 0;i<d;++i)
            if(vec1(i)>vec2(i))
                return false;
        return true;
    }

    void bringBoxToLvl(std::vector<index_t>& box,unsigned lvl)
    {
        gsTHBSplineBasis<d, T>* basis = static_cast< gsTHBSplineBasis<d,T>* > (this->m_basis);
        if(static_cast<index_t>(lvl)<box[0])
        {
            const tensorBasis & maxBasis = *(basis->getBases()[lvl]);

            T val_low,val_high;
            for(unsigned i = 0;i<m_constrainedBoxes.size();++i)
            {
                const unsigned tLvl = box[0];
                const tensorBasis & tBasis = *(basis->getBases()[tLvl]);

                if(tLvl>lvl)
                    continue;
                for (unsigned dim = 0; dim != d; dim++)
                {
                    const gsKnotVector<T> & tKv = tBasis.component(dim).knots();
                    val_low=*(tKv.ubegin()+(box[1+dim]));
                    val_high=*(tKv.ubegin()+(box[1+d+dim]));
                    const gsKnotVector<T> & maxKv = maxBasis.component(dim).knots();
                    box[1+dim] = maxKv.uFind(val_low).uIndex();
                    box[1+d+dim] = maxKv.uFind(val_high).uIndex();
                    if(val_high==*(tKv.uend()-1))
                        box[1+d+dim]++;
                }
                box[0]=lvl;
            }
        }
        else if(static_cast<index_t>(lvl)>box[0])
        {
            unsigned lvlDiff = lvl-box[0];
            for(unsigned i = 0;i<2*d;++i)
            {
                box[i+1]*=math::ipow(2,lvlDiff);
            }
        }

    }

private:
    struct constrainedBox
    {
        std::vector<index_t> m_box;
        int m_maxLvl;
    };

private:

    std::vector<constrainedBox> m_constrainedBoxes;

    using gsHFitting<d,T>::m_ext;
    using gsHFitting<d,T>::append;
    using gsHFitting<d,T>::isCellAlreadyInserted;
};

template <short_t d, class T>
void  gsHFittingLvlConstrained<d, T>::appendBox(std::vector<index_t>& boxes,
                                  std::vector<index_t>& cells,
                                  const gsVector<T>& parameter)
{
    gsTHBSplineBasis<d, T>* basis = static_cast< gsTHBSplineBasis<d,T>* > (this->m_basis);
    const index_t maxLvl = basis->maxLevel();
    const tensorBasis & tBasis = *(basis->getBases()[maxLvl]);

    // get a cell
    gsVector<index_t, d> a_cell;

    for (short_t dim = 0; dim != d; dim++)
    {
        const gsKnotVector<T> & kv = tBasis.component(dim).knots();
        a_cell(dim) = kv.uFind(parameter(dim)).uIndex();
    }

    // get level of a cell
    gsVector<index_t, d> a_cell_upp = a_cell + gsVector<index_t, d>::Ones();
    const index_t cell_lvl = basis->tree().query3(a_cell, a_cell_upp, maxLvl) + 1;

    if (isCellAlreadyInserted(a_cell, cells))
        return;

    append(cells, a_cell);

    // get the box
    std::vector<index_t> box(2 * d + 1);
    box[0] = cell_lvl;
    for (short_t dim = 0; dim != d; dim++)
    {
        const index_t numBreaks = basis->numBreaks(cell_lvl, dim) - 1 ;

        unsigned lowIndex = 0;
        if (cell_lvl < maxLvl)
        {
            const index_t shift = maxLvl - cell_lvl;
            lowIndex = (a_cell(dim) >> shift);
                //=(static_cast<util::make_unsigned<index_t>::type>(a_cell(dim)) >> shift);
        }
        else
        {
            const index_t shift = cell_lvl - maxLvl;
            lowIndex = (a_cell(dim) << shift);
        }

        // apply extensions
        index_t low = ( (lowIndex > m_ext[dim]) ? (lowIndex - m_ext[dim]) : 0 );
        index_t upp = ( (lowIndex + m_ext[dim] + 1 < numBreaks) ?
                             (lowIndex + m_ext[dim] + 1) : numBreaks );

        box[1 + dim    ] = low;
        box[1 + d + dim] = upp;
    }

    if (isBoxViolatingConstrainedBoxes(box))
        return;

    for(size_t i = 0;i<box.size();++i)
        boxes.push_back(box[i]);
}

}
