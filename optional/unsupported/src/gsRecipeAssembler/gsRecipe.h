/** @file gsRecipe.h

    @brief Many PDEs have a block structure where different blocks correspond
    to common first or second order bilinear differencial operators.
    For instance in the Poisson equation (f$ A(u,w)=B(f,w) f$) there are two
    blocks: f$Af$ and f$Bf$ where f$Af$ is the integral of the scalar product
    of the gradients (hereafter called GradGrad operator) and f$Bf$ is the
    scalar product in f$L^2f$.

    In caso of the Stokes equation the blocks are 5:
    f{eqnarray}
        A(u,w) + B(w,p) &=& C(f,w)\n
        B^t(q,u)        &=& C(s,q)
    f{eqnarray}
    where f$Af$ is again the GradGrad operator, f$Bf$ is the integral of the
    divergence of the first argument multiplided by the second argument
    (hereafter called the Div operator) and f$Cf$ is the scalar product in f$L^2f$.

    A gsRecipe is a collection of bilinear operators associated with:
    -space ids;
    -writing rules.
    Each quaddruple (operator,testSpace,trialSpace,writingRule) is called
    ingredient.

    The purpose of the gsRecipe framework is to assemble system matrices
    for problems like the above without having to write a the assembling
    procedure, but out of pre-coded blocks.
    To assemble the Poisson problem system matrix it is enough to have a
    recipe f$Rf$ with one ingredient: {opGradGrad,0,0,writeToSparseMatrix}.
    Then call the assemble method with the following arguments

                   const gsQuadRule<T>     & quadRule,
                   gsDomainIterator<T>    & domain  ,
                   const gsGeometry<T>    & geometry,
                   std::vector<spacePtr>    spaces  ,
    where
        -quadRule specifies the per element quadrature,
        -domain specifies the element list,
        -geometry is the parametrization of the domain,
        -spaces contains a spacePtr to the discrete space.

    The gsRecipe does not deduce the elements from the spaces, so that it is
    possible to mix spaces defined on different grids, does not deduce the
    quadrature so that it is possible to over- and under-integrate according
    to personal taste. Lastly, and more interestingly, it does not map the
    spaces to the physical domain. The spaces must provide data on the physical
    domain directly. This is usually achieved by the use of gsTransformedFunctionSet.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
    Created on: 2014-10-28
*/

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBoundary.h>

#pragma once
namespace gismo {

// forward declarations
template <typename T=real_t> class gsRecipe;
template <typename T>        class gsBilinearOp;
template <typename T>        class gsFunctionSet;
template <typename T>        class gsLocalToGlobalMapper;
}

namespace gismo {

/**
    \brief group together an operator, a pair of indexes of discrete
    spaces and a gsLocalToGlobalMapper that specify the destination matrix.

        The model is to construct a matrix \em A whose coefficients are
       f$ a_{i,j} = \int_{\Omega} F(\psi_i, \phi_j) f$
    where the f$\psi f$ are a basis of the test space and the f$\phi f$
    are a basis of the unknown space.


    gsRecipeIngredient assumes to \em own the operator and the storage rule
    it is provided and will thus free them automatically.
**/

template<typename T = real_t>
class gsRecipeIngredient
{
public:
    typedef memory::shared_ptr<gsBilinearOp<T>           >   opHandle  ;
    typedef memory::shared_ptr< gsLocalToGlobalMapper<T> >   ruleHandle;

public:

    const opHandle   & getOp  () const { return m_op  ; }
    const ruleHandle & getRule() const { return m_rule; }

    unsigned tSpaceId() const { return m_tSpace; }
    unsigned uSpaceId() const { return m_uSpace; }
    
    void operator= (gsRecipeIngredient other);

private:
    opHandle                    m_op;
    ruleHandle                  m_rule;
    unsigned                    m_tSpace;
    unsigned                    m_uSpace;
public:
    gsRecipeIngredient();

    gsRecipeIngredient(opHandle op, ruleHandle rule,unsigned trial, unsigned test)
        : m_op(op), m_rule(ruleHandle(rule)), m_tSpace(test), m_uSpace(trial)
    {}
    gsRecipeIngredient(gsBilinearOp<T> *op,gsLocalToGlobalMapper<T> *rule, unsigned trial, unsigned test)
        : m_op(opHandle(op)), m_rule(ruleHandle(rule)), m_tSpace(test), m_uSpace(trial)
    {}

    /**
     * \brief setOperator changes the operator associated to
     *        this recipe entry.
     *
     * Since gsRecipeIngredient will own the given operator and
     * delete it automatically we recommend to construct the
     * proper operator on the spot like in:
     *
     * setOperator(new gsGradGradOp<>() );
     *
     * \param op a pointer to a gsBilinearOp object
     */
    void setOperator (gsBilinearOp<T> *op);

    /**
     * \brief setOperator changes the operator associated to
     *        this recipe entry.
     *
     * \param op a shared_ptr pointer to a gsBilinearOp object
     */
    void setOperator (const opHandle & op);

    /**
     * \brief setOperator changes the operator associated to
     *        this recipe entry.
     *
     * Since gsRecipeIngredient will own the given operator and
     * delete it automatically we recommend to construct the
     * proper operator on the spot like in:
     *
     * setRule(new gsL2GActives<>(sys_matrix) );
     *
     * \param rule a pointer to a gsLocalToGlobalMapper object
     */
    void setRule (gsLocalToGlobalMapper<T>   *rule);

    /**
     * \brief setOperator changes the operator associated to
     *        this recipe entry.
     *
     * \param rule a shared_ptr pointer to a gsLocalToGlobalMapper object
     */
    void setRule ( const ruleHandle &         rule);

    /**
     * \brief set the test space for this ingredient
     *        (the space varying with rows)
     * \param tSpace
     */
    void setTestSpace (unsigned tSpace);

    /**
     * \brief set the unknown space for this ingredient
     *        (the space varying with columns)
     * \param uSpace
     */
    void setUnknownSpace (unsigned uSpace);
};



template <typename T=real_t>
class gsElementList : public std::vector<std::pair<gsDomainIterator<T>*,size_t> >
{
public:
    void free()
    {
        for (size_t p=0;
                p<std::vector<std::pair<gsDomainIterator<T>*,size_t> >::size();
                ++p)
            delete std::vector<std::pair<gsDomainIterator<T>*,size_t> >::at(p).first;
    }
};


/**
    \brief A recipe describes a way to assemble a system matrix and/or right-hand side

    Internally, it is a container for a list of gsRecipeIngredients
**/
template <typename T>
class gsRecipe
{
public:
    typedef typename gsRecipeIngredient<T>::opHandle   opHandle  ;
    typedef typename gsRecipeIngredient<T>::ruleHandle ruleHandle;
    typedef memory::shared_ptr<gsFunctionSet<T> >      spacePtr;
    typedef std::vector<spacePtr>                      spaceList;

public:
    gsRecipe();
    gsRecipe(size_t size);

public:// to do: document

    gsRecipeIngredient<T>& operator[](size_t i);

    const gsRecipeIngredient<T>& operator[](size_t i) const;

    size_t size() const;

    void operator= (gsRecipe<T> other);

    void append (const gsRecipe & other);

    void add(const gsRecipeIngredient<T> & ingredient);
    void add(opHandle op, ruleHandle writer,unsigned trial, unsigned test);

    void resize (size_t new_size);

    /**
        \brief assemble the matrices according to the recipe for a single patch.

        \param quadRule is the quadrature rule used in numerical integration
        \param domain is a domain iterator that provides the integration
               elements: each of these is assumed to be fully contained
               in a Bezies element for the discrete spaces
        \param geometry is the parametrization of the domain, multipatch
               geometries must be handled manually outside of this function
        \param spaces is a list of shared pointers to function sets to which
               the bilinear operators in the recipe are applied

        This function asks each operator in the recipe to know what to
        compute on each node for each space and for the geometry.

        Note that this implementation can not compete with an hand tuned
        specific equation assembler as the interfaces can introduce
        some additional cost. Anyway it is fairly efficient.
    **/

    void assemble(const gsQuadRule<T>     & quadRule,
                   gsDomainIterator<T>    & domain  ,
                   const gsGeometry<T>    & geometry,
                   spaceList                spaces  ,
                   patchSide                side=patchSide(0,boundary::none)) const;

    /**
        \brief assemble the matrices according to the recipe on multipatch.

        Contrary to the previous function this function interface tries to
        remove the burden from the programmer. There are only 2 mandatory
        parameter:

        \param domain that is a multi-patch geometry

        \param spaces that contains the discrete spaces

        This function will try to deduce reasonable default values for the
        other parameters:

        \param sides is a list of parch sides on which the integration should
               be performed. By default it is assumed to be the list of all
               patch volumes.

        \param quadRule is the quadrature rule, if not specified, the
               number of quadrature nodes in each direction is set to be the
               maximum degree of the bases in spaces plus 1.
               A minimum of 2 quadrature nodes in each direction is used if
               there is no basis in spaces.

        \param mesh conveys the partitioning of each side into elements. It
               must be an appropriate domain iterator on the correspondin
               patchSide.
               If it is not provided the spaces are searched for the first
               physical space present and the mesh is taken from it.
               If this fails the assembler will crash.

        This function can be used to assemble boundary conditions, by passing
        the list of boundary faces.

        Note that this implementation can not compete with an hand tuned
        specific equation assembler as the interfaces can introduce
        some additional cost. Anyway it is fairly efficient.
    **/
    void assemble(const gsMultiPatch<T>    &domain,
                   spaceList                spaces,
                   gsVector<index_t>            numGaussPoints,
                   const gsElementList<T>   mesh) const;
private:
    std::vector< gsRecipeIngredient<T> > m_data;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include <gsRecipeAssembler/gsRecipe.hpp>
#endif
