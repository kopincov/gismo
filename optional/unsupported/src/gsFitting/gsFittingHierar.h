/** @file gsFittingLocalHierarchical.h

    @brief Used for adaptive fitting with hierarchical splines

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsHSplines/gsHTensorBasis.h>
#include <gsFitting/gsFittingIterative.h>

namespace gismo {

/**
   \brief
   This class applies hierarchical refinement to a fitting.
   Called once the iterative fitting has converged.
*/
template <short_t d, class GenBasis, class T>
class gsFittingHierarchical :
        public gsFittingImproveAfterCV<d, GenBasis, T>
{
public:

    typedef typename gsBSplineTraits<d,T>::Basis tensorBasis;

public:
    /**
       \brief
       Main constructor
    */
    gsFittingHierarchical(GenBasis& basis,
                               T error_threshold, T extension,
                               int max_iter)
    : gsFittingImproveAfterCV<d, GenBasis, T>::gsFittingImproveAfterCV(),
    m_basis(basis)
    {
        T prop_refinement = 1.;
        initHFitting(prop_refinement, error_threshold,
                     extension, max_iter);
    }

    /// Checks where the basis should be refined and performs
    /// this refinement
    bool iterAdaptFitting(gsFittingIter<d, GenBasis, T>& fitting,
                          int iter);


public:

    /// Returns the refinement percentage
    T getRefPercentage() const
    {   return m_refinement;  }

    /// Returns the chosen cell extension
    const std::vector<unsigned> & get_extension() const 
    {   return m_extension;   }

    /// Sets the refinement percentage
    void setRefPercentage(double refPercent)
    {
        GISMO_ASSERT((refPercent >= 0) && (refPercent <= 1), "Invalid percentage" );
        m_refinement = refPercent;
    }

    /// Sets the cell extension
    void setExtension(T extension)
    {
        GISMO_ASSERT(extension > 0,
                      "Extension must be a positive number.");
        m_extension = std::vector<unsigned>
            (this->m_basis.dim(), extension);
    }

    /// Returns boxes which define refinement area.
    std::vector<std::vector<index_t> >
    getBoxes(gsLeastSquares<GenBasis, T>& leastSquares,
             const T threshold);

protected:
    /// Appends a box around parameter to the boxes only if the box is not
    /// already in boxes
    virtual void appendBox(std::vector<index_t>& boxes,
                           std::vector<index_t>& cells,
                           const gsVector<T>& parameter,
                           gsBasis<T>& basis);

    /// Identifies the threshold from where we should refine
    T setRefineThreshold(const std::vector<T>& errors);
    T setRefineThreshold(gsMatrix<T>& errors );

    /// Checks if a_cell is already inserted in container of cells
    static bool isCellAlreadyInserted(const gsVector<index_t, d>& a_cell,
                                      const std::vector<index_t>& cells);
    
    /// Appends a box to the end of boxes (This function also works for cells)
    static void append(std::vector<index_t>& boxes,
                       const gsVector<index_t>& box)
    {
        for (index_t col = 0; col != box.rows(); col++)
            boxes.push_back(box[col]);
    }

private:

    /// Called by the constructor. Sets the data of the object.
    void initHFitting(T refin, T error_threshold,
                      T extension, int max_iter)
    {
        GISMO_ASSERT((refin >= 0) && (refin <= 1),
                     "Refinement percentage must be between 0 and 1." );
        GISMO_ASSERT(extension > 0,
                      "Extension must be a positive number.");

        std::vector<unsigned> ext(d, extension);

        m_refinement    = refin;     //how many % to refine
        setExtension(extension);
        m_error_threshold = error_threshold;
        m_max_iter = max_iter;
        m_iter = 0;
    }

protected:

    /// How many % to refine in the first refinement - 0-1 interval
    T m_refinement;

    /// The maximal error tolerated. In case this threshold is
    /// not reached, the basis is refined
    T m_error_threshold;

    /// The maximal number of refinements
    int m_max_iter;

    /// The index of the iteration
    int m_iter;

    /// Size of the extension
    std::vector<unsigned> m_extension;

    /// The hierarchical basis (possibly a multipatch basis)
    GenBasis& m_basis;
};


template<short_t d, class GenBasis, class T>
bool gsFittingHierarchical<d, GenBasis, T>::
iterAdaptFitting(gsFittingIter<d, GenBasis, T>& fitting, int iter)
{
    // INVARIANT
    // look at iterativeRefine
    gsLeastSquares<GenBasis, T>& leastSquares
        = fitting.leastSquares();
    GISMO_ASSERT(! leastSquares.isEmpty(),
                 "Local refinement only implemented for the least square method (and no point have been given)");

    T max_error = leastSquares.maxError();

    if(m_iter >= m_max_iter)
        return true;
    if ( max_error > m_error_threshold)
    {
        // if err_treshold is -1, we refine the m_refinement percent of the whole domain
        /*    T threshold = (m_err_threshold >= 0) ? err_threshold
              : setRefineThreshold(pointErrors);*/

        unsigned nPatches = nBasesGen(fitting.basis());

        std::vector<std::vector<index_t> > boxes
            = getBoxes(leastSquares, m_error_threshold);
        for(unsigned patch = 0;patch < nPatches;patch++)
        {
            gsHTensorBasis<d, T>& basis =
                static_cast<gsHTensorBasis<d,T>&>
                (getBasisGen(this->m_basis, patch));
            /*      if(boxes.size()==0)
                    return false;   */

            basis.refineElements(boxes[patch]);
            if(fitting.is_displacement())
            {
                gsGeometry<T>& geom = getNthGeometry(m_basis, *fitting.m_current,
                                                     patch);
                geom.refineElements(boxes[patch]);
            }

            gsInfo << "inserted in patch " << patch << ": "
                   << boxes[patch].size() / (2 * d + 1)
               << " boxes.\n";
        }
        fitting.actualize_basis();

        m_iter++;
        fitting.exportNodes();
        return false;
    } else
    {
        gsInfo << "Tolerance for hierarchical refinement reached.\n";
        m_iter++;
        return true;
    }

}

template <short_t d, class GenBasis, class T>
std::vector<std::vector<index_t> >
gsFittingHierarchical<d, GenBasis, T>::
getBoxes(gsLeastSquares<GenBasis, T>& leastSquares,
         const T threshold)
{
    index_t nPatches = leastSquares.nPatches();
    // cells contains lower corners of elements marked for refinment from maxLevel 
    std::vector<std::vector<index_t> > cells(nPatches);
    
    // boxes contains elements marked for refinement from different levels,
    // format: { level lower-corners  upper-corners ... }
    std::vector<std::vector<index_t> >
        boxes(nPatches);

    gsPointIterator<T> it(leastSquares);

    for(;it.good();it.next())
    {
        index_t patch = it.patch();
        if (threshold <= leastSquares.error(it))
        {
            appendBox(boxes[patch], cells[patch],
                      it.currParam(), it.currentBasis(m_basis));
        }
    }
    
    return boxes;
}

/// Add a new box to boxes if parameter is not already included
/// boxes: all the boxes that will be added
/// cells: all the knot vectors to be refined
template <short_t d, class GenBasis, class T>
void gsFittingHierarchical<d, GenBasis, T>::
appendBox(std::vector<index_t>& boxes,
          std::vector<index_t>& cells,
          const gsVector<T>& parameter, gsBasis<T>& basis)
{
    gsTHBSplineBasis<d, T>& _basis =
        static_cast< gsTHBSplineBasis<d,T>& > (basis);
    const int maxLvl = _basis.maxLevel();
    const tensorBasis & tBasis = *(_basis.getBases()[maxLvl]);

    /// a_cell: the position of the given parameter in the knot vectors (d knot vectors)
    gsVector<index_t, d> a_cell;
    
    for (short_t dim = 0; dim != d; dim++)
    {
        const gsKnotVector<T> & kv = tBasis.component(dim).knots();
        a_cell(dim) = kv.uFind(parameter(dim)).uIndex();
    }
    
    if (!isCellAlreadyInserted(a_cell, cells))
    {
        append(cells, a_cell);
	
        // get level of a cell
        gsVector<index_t, d> a_cell_upp = a_cell + gsVector<index_t, d>::Ones();
        const int cell_lvl = _basis.tree().query3(a_cell, a_cell_upp, maxLvl) + 1;
	
        // get the box
        gsVector<index_t> box(2 * d + 1);
        box[0] = cell_lvl;
        for (short_t dim = 0; dim != d; dim++)
        {
            const unsigned numBreaks = _basis.numBreaks(cell_lvl, dim) - 1 ;
	    
            unsigned lowIndex = 0;
            if (cell_lvl < maxLvl)
            {
                const unsigned shift = maxLvl - cell_lvl;
                lowIndex = (a_cell(dim) >> shift);
            }
            else
            {
                const unsigned shift = cell_lvl - maxLvl;
                lowIndex = (a_cell(dim) << shift);
            }

            // apply extensions
            unsigned low = ( (lowIndex > m_extension[dim])
                             ? (lowIndex - m_extension[dim]) : 0 );
            unsigned upp = ( (lowIndex + m_extension[dim] + 1 < numBreaks)
                             ? (lowIndex + m_extension[dim] + 1) : numBreaks );
	    
            box[1 + dim    ] = low;
            box[1 + d + dim] = upp;
        }

        append(boxes, box);
    }
}

/// cells: vector of size N*d, with N the number of elements.
/// Checks that cells is not one of these elements
template <short_t d, class GenBasis, class T>
bool gsFittingHierarchical<d, GenBasis, T>::
isCellAlreadyInserted(const gsVector<index_t, d>& a_cell,
                      const std::vector<index_t>& cells)
{
    
    for (size_t i = 0; i != cells.size(); i += a_cell.rows())
    {
        int commonEntries = 0;
        for (index_t col = 0; col != a_cell.rows(); col++)
        {
            if (cells[i + col] == a_cell[col])
            {
                commonEntries++;
            }
        }
	
        if (commonEntries == a_cell.rows())
        {
            return true;
        }
    }
    
    return false;
}

template<short_t d, class GenBasis, class T>
T gsFittingHierarchical<d, GenBasis, T>::
setRefineThreshold(const std::vector<T>& errors )
{
    std::vector<T> errorsCopy = errors; 
    const size_t i = cast<T,size_t>(errorsCopy.size() * (1.0 - m_refinement));
    typename std::vector<T>::iterator pos = errorsCopy.begin() + i;
    std::nth_element(errorsCopy.begin(), pos, errorsCopy.end());
    return *pos;
}

template<short_t d, class GenBasis, class T> T
gsFittingHierarchical<d, GenBasis, T>::
setRefineThreshold(gsMatrix<T>& errors )
{
    unsigned s = errors.rows();
    std::vector<T> errorsCopy(s);
    for(unsigned j = 0;j < s;j++)
        errorsCopy.push_back(errors(j,0));

    const size_t i = cast<T,size_t>(errorsCopy.size() * (1.0 - m_refinement));
    typename std::vector<T>::iterator pos = errorsCopy.begin() + i;
    std::nth_element(errorsCopy.begin(), pos, errorsCopy.end());
    return *pos;
}


};// namespace gismo
