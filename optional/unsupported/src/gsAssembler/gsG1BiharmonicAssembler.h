/** @file gsG1BiharmonicAssembler.h

    @brief Provides assembler for a homogenius Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, P. Weinmüller
*/


#pragma once

#include <gsAssembler/gsAssembler.h>

#include <gsPde/gsG1BiharmonicPde.h>

#include <gsAssembler/gsVisitorG1Biharmonic.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorNeumannBiharmonic.h>
//#include <gsAssembler/gsVisitorNitscheBiharmonic.h>

#include "gsIO/gsParaviewCollection.h"

namespace gismo
{

/** @brief
    Implementation of a homogeneous Biharmonic Assembler.

    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can only be enforced strongly (i.e Nitsche is
    not implemented).
*/
template <class T, class bhVisitor = gsVisitorBiharmonic<T> >
class gsG1BiharmonicAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
/** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] bconditions  is a gsBoundaryConditions object that holds boundary conditions on the form:
    \f[ \text{Dirichlet: } u = g \text{ on } \Gamma, \text{ and Neumann: } \nabla \Delta u \cdot \mathbf{n} = h \text{ on } \Gamma\f]
    \param[in] bconditions2 is a gsBoundaryConditions object that holds Neumann boundary conditions on the form:
    \f[\text{Neumann: } \nabla \Delta u \cdot \mathbf{n} = g\, \rightarrow \,(g,\nabla v \cdot \mathbf{n})_\Gamma, \f] where \f$ g \f$ is the Neumann data,
    \f$ v \f$ is the test function and \f$ \Gamma\f$ is the boundary side.
    \param[in] rhs is the right-hand side of the Biharmonic equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary in the \em bconditions object.
    \param[in] intStrategy option for the treatment of patch interfaces
*/
    gsG1BiharmonicAssembler(gsMultiPatch<T> const        & patches,
                           gsMultiBasis<T> const         & bases,
                           std::vector<gsMultiPatch<>>  & G1Basis,
                           gsBoundaryConditions<T> const & bconditions,
                           gsBoundaryConditions<T> const & bconditions2,
                           const gsFunction<T>           & rhs,
                           dirichlet::strategy           dirStrategy,
                           iFace::strategy               intStrategy = iFace::glue)
    : m_ppde(patches,bconditions,bconditions2,rhs,G1Basis)
    {
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        Base::initialize(m_ppde, bases, m_options);
    }

    void refresh();
    
    void assemble();

    void constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>& result, int unk = 0);

    void writeParaview(const gsField<T> & field,
                       std::string const & fn,
                       unsigned npts = 1000, bool mesh = false);

    std::vector< gsMultiPatch<>> G1Basis() {return * m_ppde.G1Basis();}

protected:

    // fixme: add constructor and remove this
    gsG1BiharmonicPde<T> m_ppde;

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

};

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::constructSolution(const gsMatrix<T>& solVector,
                                       gsMultiPatch<T>& result, int unk)
{
    const gsDofMapper & mapper = m_system.colMapper(unk); // unknown = 0
    result.clear();

    const index_t dim = ( 0!=solVector.cols() ? solVector.cols() :  m_ddof[unk].cols() );

    gsMatrix<T> coeffs;
    for (index_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
    {
        // Reconstruct solution coefficients on patch p
        const int sz  = m_bases[m_system.colBasis(unk)][p].size();
        coeffs.resize(sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if ( mapper.is_free(i, p) ) // DoF value is in the solVector
            {
                gsInfo << "mapper index: " << mapper.index(i, p) << "\n";
                coeffs.row(i) = solVector.row( mapper.index(i, p) );
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                //coeffs.row(i) = m_ddof[unk].row( mapper.bindex(i, p) ).head(dim);
                gsMatrix<T> d_0(1,1);
                d_0 << 0;
                coeffs.row(i) = d_0; // Dirichlet = 0
            }
        }

        result.addPatch( m_bases[m_system.colBasis(unk)][p].makeGeometry( give(coeffs) ) );

        // Basis G1:
        const index_t n = m_ppde.G1Basis()->at(p).nPatches();

        for (index_t i = 0; i < n; i++)
        {
            index_t m = m_ppde.G1Basis()->at(p).patch(i).basis().size();
            coeffs.resize(m, dim); // Coefficients pro Basis function

            coeffs = m_ppde.G1Basis()->at(p).patch(i).coefs();

            if ( mapper.is_free(i+sz, p) ) // DoF value is in the solVector
            {
                gsInfo << "mapper index 2: " << mapper.index(i+sz, p) << "\n";
                coeffs *=  solVector.row( mapper.index(i+sz, p) );
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                gsMatrix<T> d_0(1,1);
                //d_0 = m_ddof[unk].row( mapper.bindex(i, p) ).head(dim);
                d_0 << 0;
                coeffs *= d_0; // Dirichlet = 0
            }

            m_ppde.G1Basis()->at(p).patch(i).setCoefs(coeffs);
        }
    }
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::refresh() // Only for 2 Patches // TODO
{
    gsMatrix<index_t> act1, act2;

    const boundaryInterface & iFace = *m_ppde.patches().iBegin();// assume one single interface
    const gsAffineFunction<> ifaceMap(m_ppde.patches().getMapForInterface(iFace));


    // 1. Obtain a map from basis functions to matrix columns and rows
    gsVector<index_t> sz(m_ppde.patches().nPatches()); // n Patches
    for (index_t i = 0; i < sz.size(); i++)
        sz[i] = m_ppde.patches().patch(i).basis().size() + m_ppde.G1Basis()->at(i).nPatches();

    gsBasis<> & B1 = m_ppde.patches().patch(iFace.first().patch).basis();
    gsBasis<> & B2 = m_ppde.patches().patch(iFace.second().patch).basis();

    gsDofMapper map(sz);

    gsMatrix<index_t> boundary;
    boundary = B1.allBoundary();
    map.markBoundary(0,boundary); // Patch 0
    boundary = B2.allBoundary();
    map.markBoundary(1,boundary); // Patch 1

    for (index_t i = 0; i < m_ppde.patches().boundaries().size(); i++)
    {
        boundary = B1.boundaryOffset(m_ppde.patches().boundaries()[i],1);
        map.markBoundary(0,boundary);
        map.markBoundary(1,boundary);
    }

    // glue interface
    B1.matchWith(iFace, B2, act1, act2);
    // mark dofs
    map.markBoundary(iFace.first().patch, act1);//interface
    act1 = B1.boundaryOffset(iFace.first() .side(), 1);
    map.markBoundary(iFace.first().patch, act1); //first adj. face
    act2 = B2.boundaryOffset(iFace.second().side(), 1);
    map.markBoundary(iFace.second().patch, act2);//second adj. face

    index_t xx = 8; // numRef + 3

    gsVector<index_t > coupled_g1(xx+xx);
    coupled_g1.setLinSpaced(xx+xx,m_bases[0].basis(0).size(),m_bases[0].basis(0).size() + xx+xx-1);
    gsInfo << "coupled " << coupled_g1 << "\n";
    map.matchDofs(iFace.second().patch,coupled_g1,iFace.first().patch,coupled_g1);

    gsMatrix<index_t > boundary_g1(8,1);
    boundary_g1 << m_bases[0].basis(0).size(), m_bases[0].basis(0).size() +1 ,
        m_bases[0].basis(0).size() + xx-2, m_bases[0].basis(0).size() + xx-1,
        m_bases[0].basis(0).size() + xx, m_bases[0].basis(0).size() + xx+1,
        m_bases[0].basis(0).size() + xx+xx-2, m_bases[0].basis(0).size() + xx+xx-1;
    map.markBoundary(iFace.first().patch, boundary_g1);
    map.markBoundary(iFace.second().patch, boundary_g1);

    map.finalize();
    map.print(); // map without interface

    // =========================================
/*
    gsInfo << "act1: " << act1 << "\n";
    gsInfo << "act2: " << act2 << "\n";

    gsMatrix<unsigned> globals;
    map.localToGlobal(act1,iFace.first().patch,globals);
    gsInfo << "map: " << globals << "\n";

    //gsInfo << "map: " << map.asVector() << "\n";
*/

    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(map);//1,1

    gsInfo << "\nMatrix size: " << m_system.matrix().dim() << "\n";
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::assemble()
{
    GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
    m_system.reserve(nz, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    Base::computeDirichletDofs();

    // Assemble volume integrals
    Base::template push<bhVisitor >();
    
    // Neumann conditions of first kind
    //Base::template push<gsVisitorNeumann<T> >(
    //    m_ppde.bcFirstKind().neumannSides() );

    // Neumann conditions of second kind
    //Base::template push<gsVisitorNeumannBiharmonic<T> >(
    //   m_ppde.bcSecondKind().neumannSides() );

    if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
        gsWarn <<"DG option ignored.\n";

    /*
    // If requested, force Dirichlet boundary conditions by Nitsche's method
    this->template push<gsVisitorNitscheBiharmonic<T> >(
    m_ppde.bcSecondKind().dirichletSides() );
    */
    
    // Assembly is done, compress the matrix
    Base::finalize();
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::writeParaview(const gsField<T> & field,
                                                         std::string const & fn,
                                                         unsigned npts, bool mesh)
{
    const unsigned n = field.nPieces();
    gsParaviewCollection collection(fn);
    std::string fileName;

    for ( unsigned i=0; i < n; ++i ) // Patches
    {
        //const gsBasis<T> & dom = field.isParametrized() ?
        //                         field.igaFunction(i).basis() : field.patch(i).basis();

        fileName = fn + util::to_string(i);
        //writeSinglePatchField( field, i, fileName, npts );

        const gsFunction<T> & geometry = field.patch(i);
        const gsFunction<T> & parField = field.function(i);

        const int n = geometry.targetDim();
        const int d = geometry.domainDim();

        gsMatrix<T> ab = geometry.support();
        gsVector<T> a = ab.col(0);
        gsVector<T> b = ab.col(1);

        gsVector<unsigned> np = uniformSampleCount(a, b, npts);
        gsMatrix<T> pts = gsPointGrid(a, b, np);

        gsMatrix<T> eval_geo = geometry.eval(pts);//pts
        gsMatrix<T> eval_field = field.isParametric() ? parField.eval(pts) : parField.eval(eval_geo);

        // Hier g1 basis dazu addieren!!!
        //eval_field.setZero();

        gsWriteParaview(m_ppde.G1Basis()->at(0).patch(2),"testL");
        gsWriteParaview(m_ppde.G1Basis()->at(1).patch(2),"testR");

        const index_t m = m_ppde.G1Basis()->at(i).nPatches(); // number of g1 basis function
        for (index_t j = 0; j < m ; j++)
        {
            gsField<T> plotBasisG1(field.patch(i),m_ppde.G1Basis()->at(i).patch(j));

            //gsInfo << "function " << parField_g1 << "\n";
            if (i == 0)
                eval_field += plotBasisG1.value(pts);
            if (i == 1)
                eval_field += plotBasisG1.value(pts);
        }

        if ( 3 - d > 0 )
        {
            np.conservativeResize(3);
            np.bottomRows(3-d).setOnes();
        }
        else if (d > 3)
        {
            gsWarn<< "Cannot plot 4D data.\n";
            return;
        }

        if ( 3 - n > 0 )
        {
            eval_geo.conservativeResize(3,eval_geo.cols() );
            eval_geo.bottomRows(3-n).setZero();
        }
        else if (n > 3)
        {
            gsWarn<< "Data is more than 3 dimensions.\n";
        }

        if ( eval_field.rows() == 2)
        {
            eval_field.conservativeResize(3,eval_geo.cols() );
            eval_field.bottomRows(1).setZero(); // 3-field.dim()
        }

        gsWriteParaviewTPgrid(eval_geo, eval_field, np.template cast<index_t>(), fileName);



        collection.addPart(fileName, ".vts");
        if ( mesh )
        {
            fileName+= "_mesh";
            //writeSingleCompMesh(dom, field.patch(i), fileName);
            collection.addPart(fileName, ".vtp");
        }

    }
    collection.save();



}


} // namespace gismo



