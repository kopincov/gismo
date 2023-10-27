/** @file gsRemappedBasis.h

    @brief Implementation of a basis that delegates evaluation
    to other basis locally.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
**/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsDofMapper.h>
//#include <gsRemappedBasis/gsFunctionComputeOnly.h>
#include <gsRemappedBasis/gsSelector.h>
#include <gsMSplines/gsWeightMapper.h>

namespace gismo {


/**
   @brief Implementation of the general spline basis described in \cite bm2016.

   Main idea
   =========

   The computational domain is partitioned into several disjoint \em
   subdomains. For each subdomain \f$ D_{\ell}\f$ a \em local \em
   basis is stored so that each \em global basis function can be on
   \f$ D_{\ell} \f$ expressed using the local basis.

   Example
   -------

   If the global basis are THB-splines with three levels, we can
   choose the subdomains \f$ D_0 = \Omega_0\f$, \f$D_1 = \Omega_1
   \setminus \Omega_0\f$ and \f$ D_2 = \Omega_2 \setminus \Omega_1
   \f$. For the \f$\ell\f$--th local basis we can take \em all the
   B-splines on the knot vector of level \f$\ell\f$.

   Ingredients
   -----------

   - Mapper \a m_repr is a wrapper around a matrix holding the
     coefficients expressing the \em global basis in terms of the \em
     local bases.
   - Selector \a m_sele finds in which subdomain a given point is.
   - Vector of local bases \a m_basis.
 */
//template<class T>
class GISMO_EXPORT gsRemappedBasis : public gsBasis<real_t>
{
public:
    /// Shared pointer for gsRemappedBasis
    typedef memory::shared_ptr< gsRemappedBasis > Ptr;

    /// Unique pointer for gsRemappedBasis
    typedef memory::unique_ptr< gsRemappedBasis > uPtr;

    typedef unsigned basisIdT;
    typedef short_t directionT;
    typedef size_t  patchIdT;
    typedef gsFunctionSet<real_t>::Ptr basisPtr;

public: // Constructors from components
    gsRemappedBasis( const gsWeightMapper<real_t>               &repr,
                     const gsSelector                           &sele,
                     const std::vector<gsFunctionSet<real_t>* > &basis );

    gsRemappedBasis( const gsWeightMapper<real_t>   &repr,
                     const gsSelector               &sele,
                     const std::vector<basisPtr>    &basis );

    gsRemappedBasis( const gsRemappedBasis &other);

public:

    virtual memory::unique_ptr<gsGeometry<real_t> > makeGeometry(gsMatrix<real_t> ) const // coefs ) const
    {
        return memory::unique_ptr<gsGeometry<real_t> >();
    }
    virtual std::ostream &print(std::ostream &os) const
    {
        os<<"No description";
        return os;
    }

    virtual index_t size()   const;

    virtual short_t domainDim () const;
    virtual short_t targetDim () const;
    virtual gsMatrix<real_t> support () const {return m_sele.patch(0).getBoundingBox();}

    virtual void compute(const gsMatrix<real_t>  &points, gsFuncData<real_t> & result) const;
    virtual void compute(const gsMapData<real_t> &geo,    gsFuncData<real_t> & result) const;

    virtual void eval_into(const gsMatrix<real_t> &points, gsMatrix<real_t> &result) const;

public: // gsRemappedBasis specific interface
    /// Obtain a new remapped function set whose basis-functions are the
    /// linear combination of the basis functions of the current speciefied by
    /// the columns of coefs
    template <typename MatrixT>
    void reduce (const MatrixT & coeff);
    void reduce (const gsDofMapper &coeff);

    /// as the function before but allows to choose if the returned object
    /// behaves as a basis or a function as the resulting target dimension.
    /// If the target dimension is -1 then the target dimension is preserved
    /// otherwise the target dimension of the resulting object is set to the provided value.
    template <typename MatrixT>
    gsRemappedBasis* makeReduced (const MatrixT & coeff, bool isBasis, short_t targetDim=-1) const;

    /// takes many single patch remapped bases and construct a multipatch remapped basis.
    static gsRemappedBasis* makeMultiPatch(const std::vector<basisPtr> input, const gsWeightMapper<real_t> *gluing =NULL);
    static gsRemappedBasis* makeMultiPatch(const gsMultiBasis<real_t> &input, const gsWeightMapper<real_t> *gluing );
    static gsRemappedBasis* makeMultiPatch(const gsMultiBasis<real_t> &input);

public: // Utility functions

    /// Exports the mesh to the file \a filename .tex, which you then compile with pdflatex and watch as a pdf.
    /// \param boundingBox matrix describing the outer boundary of the mesh.
    void exportSelectorToTex(std::string filename) const;

    /// retrieve the selector
    const gsSelector&             getSelector () const;
    /// retrieve the mapper
    const gsWeightMapper<real_t>& getMapper   () const;
    /// retrieve the local basis
    const std::vector<basisPtr>&  getBases    () const;

    GISMO_CLONE_FUNCTION(gsRemappedBasis)

protected: // Helpers for derived classes
    // Empty constructor is protected for use by derived classes.
    gsRemappedBasis();
    void checkDimAndInitShifts();

private: // Implementations
    void getActiveAndMatrix(gsMatrix<index_t> &localAct, patchIdT basisId, gsWeightMapper<real_t>::IndexContainer &globalAct, gsMatrix<real_t> &Rm) const;

    template <typename T>
    void copyCompleteCol (index_t col, gsMatrix<T> & dest, const gsMatrix<T> &data) const;

    void remapValues (
            const gsMatrix<real_t>  &matrix,
            const gsMatrix<real_t>  &input,
                  gsMatrix<real_t>  &output
            ) const;

    /// Interior common to the constructors.
    inline void init();

protected: // members
    gsWeightMapper<real_t>  m_repr; // matrix representing global functions as linear combinations
    gsSelector              m_sele; // choose which basis to use based on point
    std::vector<basisPtr>   m_basis;
    std::vector<index_t>    m_shift; // {Bi} set of local basis :  0, dim B1, dim B1+dim B2 ..., sum_i dim Bi

    std::pair<short_t, short_t> m_info;

    /** true if the output should use the basis format:
     *  i.e., suppress the zero functions
     *  false if also 0 functions must be reported */
    bool                    m_isBasis;
};


template <typename MatrixT>
void gsRemappedBasis::reduce (const MatrixT & coeff)
{
    m_repr*=coeff;
}

template <typename MatrixT>
gsRemappedBasis* gsRemappedBasis::makeReduced (const MatrixT & coeff, bool isBasis, short_t targetDim) const
{
    gsRemappedBasis * result = new gsRemappedBasis();
    result->m_repr.asMatrix().resize(m_repr.asMatrix().rows(),coeff.cols());
    result->m_repr.asMatrix() = m_repr.asMatrix() * coeff;
    result->m_repr.optimize();
    result->m_isBasis=isBasis;
    result->m_sele=m_sele;
    result->m_basis=m_basis;
    result->m_info=m_info;
    result->m_info.second = targetDim == -1? coeff.cols(): targetDim;
    result->m_shift=m_shift;
    return result;
}

} // namespace gismo
