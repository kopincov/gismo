
#pragma once

#include <iostream>
#include <gsAssembler/gsAssembler.h>
#include <gsPde/gsSurfacePoissonPde.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{    

template <class T> class gsSparseSystem0;

template<class T>
class gsAssembler0
{
public:

    /// Default empty constructor
    gsAssembler0() : m_geometry(NULL) { }
    
    /// Constructor using a geometry
    gsAssembler0(const gsGeometry<T> & geom)
    { 
        m_geometry = &geom;
    }
    
    virtual ~gsAssembler0() { } //destructor
    
public:

    /// Get the assembler's geometry.
    const gsGeometry<T> & geometry() const    { return *m_geometry; }

    /// Set the assembler's geometry.
    virtual void setGeometry(const gsGeometry<T> & geom)
    { 
        m_geometry = &geom;
    } 

    /** @brief
        Assemble a sparse linear system for the given PDE on the assembler's geometry.

        \param basis    the discretization basis
        \param mapper   the gsDofMapper which describes interior, interface, and boundary dofs
        \param ddof     values for eliminated Dirichlet dofs
        \param pde      the partial differential equation to be assembled
        \param patchIndex index of the patch currently being assembled
        \return the sparse linear system that was assembled
     */
    virtual gsSparseSystem0<T>
    assemble( const gsBasis<T>& basis, const gsDofMapper& mapper, 
              const gsMatrix<T> & ddof, const gsPde<T> & pde, int patchIndex=0) = 0;

    /** @brief
        Assemble the mass matrix for the given discretization basis on the assembler's geometry.

        \param basis    the discretization basis
        \param mapper   the gsDofMapper which describes interior, interface, and boundary dofs
        \param patchIndex index of the patch currently being assembled
        \return         the mass matrix as a sparse matrix
     */
    virtual gsSparseMatrix<T> *
    assembleMass( const gsBasis<T>& basis, const gsDofMapper& mapper, int patchIndex=0 )
    { GISMO_NO_IMPLEMENTATION }

    /// \brief Returns the stiffness matrix (for Laplace/Poisson problem) for a given basis.
    ///
    /// \warning If the problem is symmetric, only the lower part is filled!
    ///
    /// To get a full view use
    ///
    /// gsSparseMatrix<T> Afull;
    ///  Afull = memory::make_unique(A)->template selfadjointView<Lower>();
    ///
    /// or
    ///
    ///    gsSparseMatrix<> Afull;
    ///    Afull = A.selfadjointView<Lower>();
    ///
    ///
    ///
    /// Or to copy the matrix to the complete matrix use
    ///
    /// gsSparsematrix<T> Acomplete = A.selfadjointView<Lower>();
    ///
    ///
    virtual gsSparseMatrix<T> * stiffness( const gsBasis<T>& B )
    { GISMO_NO_IMPLEMENTATION }


    gsSparseMatrix<T> * stiffnessFull( const gsBasis<T>& B )
    {
        gsSparseMatrix<T> * L = new gsSparseMatrix<T>;
        *L = memory::make_unique( stiffness(B) )->template selfadjointView<Lower>();
        L->makeCompressed();
        return L;
    }

    /// \brief Returns the mass matrix for a given basis.
    ///
    /// \warning If the problem is symmetric, only the lower part is filled!
    virtual gsSparseMatrix<T> * massMatrix( const gsBasis<T>& B )
    { GISMO_NO_IMPLEMENTATION }

    gsSparseMatrix<T> * massMatrixFull( const gsBasis<T>& B )
    {
        gsSparseMatrix<T> * L = new gsSparseMatrix<T>;
        *L = memory::make_unique( massMatrix(B) )->template selfadjointView<Lower>();
        L->makeCompressed();
        return L;
    }

    virtual gsVector<T> * moments( const gsBasis<T>& B, gsFunction<T> const& f )
    { GISMO_NO_IMPLEMENTATION }

    /// Compute Boundary moments of basis B with respect to function f
    /// along side s of geometry
    virtual gsVector<T> * boundaryMoments( const gsBasis<T>& B,
                                           gsFunction<T> const& f,
                                           boxSide const& s )
    { GISMO_NO_IMPLEMENTATION }

    /// Add contribution of boundary condition \a bc to the system \a system
    /// \param B        discretization basis for the patch
    /// \param bc       the boundary condition to apply
    /// \param mapper   contains the global Dof map for the \a system
    /// \param system   the linear system where the contribution shall be added
    virtual void applyBoundary( const gsBasis<T> & B,
                        const boundary_condition<T> & bc,
                        const gsDofMapper& mapper,
                        gsSparseSystem0<T> & system )
    { GISMO_NO_IMPLEMENTATION }

    /// Add contribution of interface \a bi to the system \a system
    /// \param B1       discretization basis for the first patch
    /// \param B2       discretization basis for the second patch
    /// \param geo2     the adjacent patch to m_geometry
    /// \param bi       description of the interface which joins the two patches
    /// \param mapper   contains the global Dof map for the \a system
    /// \param system   the linear system where the contribution shall be added
    virtual void applyDG( const gsBasis<T> & B1, const gsBasis<T> & B2,
                          const gsGeometry<T> & geo2,
                          const boundaryInterface & bi,
                          const gsDofMapper& mapper,
                          gsSparseSystem0<T> & system )
    { GISMO_NO_IMPLEMENTATION }


    // generic helper functions for assembling routines


    /// Add contributions from local to global stiffness matrix, identity DOF mapping
    static void localToGlobal(const gsMatrix<T>& localStiffness, 
                              const gsMatrix<index_t>& localDofs,
                              gsSparseMatrix<T>& K, 
                              bool symmetric);

    /// Add contributions from local stiffness matrix/load vector to
    /// global stiffness matrix/load vector, eliminating Dirichlet BCs
    static void localToGlobal_withBC(const gsMatrix<T>& localStiffness, 
                                     const gsVector<T>& localRhs, 
                                     const gsMatrix<T>& dirbc, 
                                     const gsDofMapper& mapper, 
                                     const gsVector<index_t>& loc2glob,
                                     gsSparseMatrix<T>& K, gsVector<T>& f, 
                                     bool symmetric);

    /// Add contributions from local to global stiffness matrix, eliminating Dirichlet BCs
    static void localToGlobal_withBC(const gsMatrix<T>& localStiffness,
                                     const gsDofMapper& mapper,
                                     const gsVector<index_t>& loc2glob,
                                     gsSparseMatrix<T>& K,
                                     bool symmetric);

    /// Allocate a sparse matrix and reserve enough space in it for assembling using the given matrix
    static gsSparseMatrix<T> *initMatrix(int nDofs, const gsBasis<T>& basis, bool symmetric);

    /// Allocate sparse matrix and right-hand side and initialize them
    static gsSparseSystem0<T> initLinearSystem(int nDofs, const gsBasis<T>& basis, bool symmetric);

    T getMu(const gsBasis<T>& b)
    {
        // TODO: basis should provide a meshSize() method
        // BUG: Nitsche parameter should depend on size of local cell, not global mesh size
        const T h = math::pow( (T) b.size(), -1.0 / b.dim() );
        const T bdeg = (T)b.degree(0);
        return ( (bdeg+b.dim())* (bdeg+1) * 2.0 / h );
        //return ( 2.0 / (h * h) );
    }

// Data members
protected:

    const gsGeometry<T> * m_geometry;

}; // class gsAssembler0


/// @brief A linear system consisting of sparse matrix and right hand side.
///
/// This is a very simple class which only holds pointers to a sparse matrix
/// and a right-hand side vector. It does no processing on its own.
template <class T>
class gsSparseSystem0
{
public:
    /// Empty constructor.
    gsSparseSystem0()
    { m_matrix = 0; m_rhs = 0; }

    /// Constructor from given matrix and right-hand side
    gsSparseSystem0(gsSparseMatrix<T>* m, gsVector<T>* rhs)
    { m_matrix = m; m_rhs = rhs; }

    /// Return the matrix.
    gsSparseMatrix<T>* matrix() const  { return m_matrix; }

    /// Return the right-hand side.
    gsVector<T>*       rhs   () const  { return m_rhs; }

    /// Return the number of degrees of freedom in the system.
    index_t size()      { return this->matrix()->cols(); }

    /// Delete the matrix and right-hand side.
    void free()
    {
        delete m_matrix; m_matrix = 0;
        delete m_rhs;    m_rhs = 0;
    }

    /// Print the linear system to a stream.
    friend std::ostream &operator<<(std::ostream &os, const gsSparseSystem0<T> & system)
    { 
        os << "gsSparseSystem0:\n";
        if ( system.m_matrix )
            os << "* Left hand side:\n" << *system.m_matrix ;    
        if ( system.m_rhs )
            os << "Right hand side (transposed):\n"  << system.m_rhs->transpose() ;
        os<<"\n";
        return os;
    }

private:
    gsSparseMatrix<T> *m_matrix;
    gsVector<T> *m_rhs;
};


/** @brief
    Implementation of a surface assembler using Gauss quadrature.
*/
    
template<class T>
class gsGaussSurfaceAssembler : public gsAssembler0<T>
{
public:
    /// Default empty constructor
    gsGaussSurfaceAssembler() : gsAssembler0<T>() { }
    
    /// Construct by a provided geometry
    gsGaussSurfaceAssembler(const gsGeometry<T> & geom) : gsAssembler0<T>(geom)
    {
        this->setGeometry( geom );
    }
    
    ~gsGaussSurfaceAssembler()                 //destructor
    { }

    virtual void setGeometry(const gsGeometry<T> & geom)
    { 
      this->m_geometry = &geom;
      // Caution: integration points should be set wrt the degree of
      // the discretization basis, not the geometry map
    }

public:

    static gsVector<index_t> getNumIntNodesFor(const gsBasis<T>& b)
    {
        gsVector<index_t> numNodes( b.dim() );
        for (int i = 0; i < b.dim(); ++i)
            numNodes[i] = b.degree(i) + 1;
        return numNodes;
    }

    virtual gsSparseSystem0<T>
    assemble( const gsBasis<T>& basis, const gsDofMapper& mapper2, 
              const gsMatrix<T> & ddof, const gsPde<T> & pde, int patchIndex=0)
    {
        GISMO_NO_IMPLEMENTATION
    }

    virtual gsSparseSystem0<T>
    assemble( const gsBasis<T>& basis,
              const gsMatrix<T> & ddof, const gsPde<T> & pde, int patchIndex=0)
    {
        const gsSurfacePoissonPde<T> *surfacepoisson = 
            dynamic_cast<const gsSurfacePoissonPde<T>*>(&pde);
        if (surfacepoisson)
            return assembleSurfacePoisson(basis, mapper,ddof, *surfacepoisson, patchIndex);

        GISMO_ERROR("Unknown PDE type in assemble()");
    }

    /// Assembler for single patch Surface Poisson equation
    gsSparseSystem0<T> assembleSurfacePoisson( const gsBasis<T>& basis, 
                                              const gsDofMapper& mapper,
                                              const gsMatrix<T> & ddof, 
                                              const gsSurfacePoissonPde<T> & pde, 
                                              index_t patchIndex=0);

    void applyBoundary( const gsBasis<T>   & B,
                        const boundary_condition<T> & bc,
                        const gsDofMapper& mapper,
                        gsSparseSystem0<T> & system );
    
    void applyDG(const gsBasis<T> & B1, const gsBasis<T> & B2,
                  const gsGeometry<T> & geo2,
                  const boundaryInterface & bi,
                  const gsDofMapper& mapperParam,
                  gsSparseSystem0<T> & system );

private:

    /// Add contribution of Nitsche Dirichlet boundary to matrix K
    /// \param B is a boundary basis, \param f is the Dirichlet function
    void boundaryNitsche( const gsBasis<T>   & B,
                          const index_t patch,
                          const boxSide s,
                          const gsFunction<T> & f,
                          gsSparseSystem0<T> & system );

    /// Add contribution of Neumann boundary condition to he \a system
    /// \param B is a boundary basis, \param f is the Dirichlet function    
    void boundaryNeumann( const gsBasis<T>   & B,
                          const index_t patch,
                          const boxSide s,
                          const gsFunction<T> & f,

                          gsSparseSystem0<T> & system );

    static gsVector<index_t> getNumIntNodesForSide(const gsBasis<T>& b, int dir)
    {
        gsVector<index_t> numNodes ( b.dim() );
        for (int i = 0; i < dir; ++i)
            numNodes[i] = b.degree(i) + 1;
        numNodes[dir] = 1;
        for (int i = dir+1; i < b.dim(); ++i)
            numNodes[i] = b.degree(i) + 1;
        return numNodes;
    }

    static gsVector<index_t> getNumIntNodesForInterface(const gsBasis<T>& b1, const gsBasis<T>& b2,
                           const boundaryInterface & bi, bool left = true)
    {
        // assumes matching orientation
        gsVector<index_t> numNodes ( b1.dim() );
        const int dir = ( left ? bi.first().direction() : bi.second().direction() );
        for (int i = 0; i < dir; ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        numNodes[dir] = 1;
        for (int i = dir+1; i < b1.dim(); ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        return numNodes;
    }

    static gsVector<index_t> getNumIntNodesForCoupled(const gsBasis<T>& b1, const gsBasis<T>& b2)
    {
        gsVector<index_t> numNodes ( b1.dim() );
        for (int i = 0; i < b1.dim(); ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        return numNodes;
    }

public:

    gsDofMapper mapper;

}; // class gsGaussSurfaceAssembler

//////////////////////////////////////////////////
//////////////////////////////////////////////////


template<class T>
gsSparseSystem0<T>
gsGaussSurfaceAssembler<T>::assembleSurfacePoisson( const gsBasis<T>& basis, 
                                                    const gsDofMapper& mapper2, 
                                                    const gsMatrix<T> & ddof, 
                                                    const gsSurfacePoissonPde<T> & pde, 
                                                    index_t patchIndex)
{
    const short_t d = 2;

    gsVector<index_t> numNodes = getNumIntNodesFor( basis );
    gsSparseSystem0<T> sys = this->initLinearSystem( mapper2.freeSize(), basis, pde.isSymmetric() );
    
    gsMatrix<T> rhsVals;        // values of the right-hand side
    gsMatrix<T> trf_grads_k;    // transformed (physical) gradients of all basis functions at one quadrature node
    
    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell
    gsVector<T> localRhs;       // local load vector

    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM, this->geometry()));

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // Set the number of integration points for each element
    //domIt->computeQuadratureRule( numNodes );
    gsGaussRule<T> quRule(numNodes);
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;
    gsFuncData<T> bdata(NEED_DERIV|NEED_ACTIVE|SAME_ELEMENT);
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

    gsMatrix<index_t> activeFuncs;

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        basis.active_into(domIt->center, activeFuncs);
        mapper2.localToGlobal(activeFuncs, patchIndex, activeFuncs);
        const index_t numActive = activeFuncs.rows();
        basis.compute(quNodes, bdata);

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(quNodes);

        // evaluate right-hand side at the geometry points
        if (pde.rhs())
            pde.rhs()->eval_into( geoEval->values(), rhsVals );

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);
        localRhs.setZero(numActive);

        for (index_t k = 0; k < bdata.allValues().cols(); ++k)      // loop over quadrature nodes
        {
            const gsMatrix<T> & Jk = 
                geoEval->jacobians().block(0, k*d, d+1, d);
            gsMatrix<T> FirstFund = Jk.transpose()*Jk;
            const T weight = quWeights[k] * math::sqrt(FirstFund.determinant());
            FirstFund = FirstFund.inverse();

            trf_grads_k.resize( 2, numActive );
            for (index_t i = 0; i < numActive; ++i)
            {
                trf_grads_k.template block<2, 1>(0, i) =
                    bdata.values[1].template block<2, 1>(i * 2, k);
            }
            
            if (pde.rhs())
                localRhs += (weight * rhsVals(0,k)) * bdata.values[0].col(k);
            
            localStiffness.noalias() += weight * (trf_grads_k.transpose() *FirstFund* trf_grads_k
                                       //Plus a mass term
                                       + bdata.values[0].col(k) * bdata.values[0].col(k).transpose()
                                     );


        }  // end loop Gauss nodes
        
        // add contributions from local stiffness matrix to global
        // stiffness matrix and load vector
        this->localToGlobal_withBC(localStiffness, localRhs, ddof, mapper2, 
                                   activeFuncs, *sys.matrix(), *sys.rhs(), true);
    } // end loop over all domain elements

    sys.matrix()->makeCompressed();
    return sys;
}


template<class T> void 
gsGaussSurfaceAssembler<T>::applyDG( const gsBasis<T>        & B1, 
                                     const gsBasis<T>        & B2,
                                     const gsGeometry<T>     & geo2,
                                     const boundaryInterface & bi,
                                     const gsDofMapper       & mapperParam,
                                     gsSparseSystem0<T>       & system )
{
    const gsGeometry<T> & geo1 = *this->m_geometry;
    gsSparseMatrix<T> & lhs    = *system.matrix();
    const int d                = this->m_geometry->parDim() ;
    const T mu                 = this->getMu(B1);
    
    //gsDebug<<"Apply surf. DG on "<< bi <<"(mu="<<mu<<").\n";
    
    const index_t patch1       = bi.first().patch;
    const index_t patch2       = bi.second().patch;
    const boxSide side1        = bi.first().side();
    const boxSide side2        = bi.second().side();
    const T orient = sideOrientation(side1);
    //const T orient2 = sideOrientation(side2);
    //gsDebugVar( orient );
    //gsDebugVar( orient2 );

    const int dir1             = side1.direction();
    //const int dir2             = side2.direction();

    //GISMO_ASSERT( B1.component(!dir1).size() == B2.component(!dir2).size(), 
    //              "DG method not implemented yet for non matching interfaces");
    
    // Quadrature for boundary integrals
    gsVector<index_t> intNodes1 = getNumIntNodesForInterface( B1, B2, bi, true  );
    gsVector<index_t> intNodes2 = getNumIntNodesForInterface( B1, B2, bi, false );
    
    // Evaluators for the two patches
    typename gsGeometryEvaluator<T>::uPtr geoEval1(getEvaluator(NEED_VALUE| NEED_JACOBIAN| NEED_GRAD_TRANSFORM, geo1));
    typename gsGeometryEvaluator<T>::uPtr geoEval2(getEvaluator(NEED_VALUE| NEED_JACOBIAN| NEED_GRAD_TRANSFORM, geo2));
    
    // Temporaries
    gsMatrix<T> grads_k_1, grads_k_2;
    gsVector<T> unormal(d);
    
    gsMatrix<T> B11, B22, B12, B21, 
                E11, E22, E12, E21,
                N1 , N2;

    // iterator on grid cells on the "right"
    typename gsDomainIterator<T>::uPtr domIter2= B2.makeDomainIterator(side2);

    //int count = 0;
    // iterate over all boundary grid cells on the "left"

    gsGaussRule<T> quRule1(intNodes1);
    gsGaussRule<T> quRule2(intNodes2);
    gsMatrix<T> quNodes1, quNodes2;
    gsVector<T> quWeights1, quWeights2;
    gsFuncData<T> bdata1(NEED_DERIV|NEED_ACTIVE|SAME_ELEMENT);
    gsFuncData<T> bdata2(NEED_DERIV|NEED_ACTIVE|SAME_ELEMENT);


    for (typename gsDomainIterator<T>::uPtr domIter1 = B1.makeDomainIterator(side1); 
         domIter1->good(); domIter1->next())
    {
        // Compute the quadrature rule on both sides
        quRule1.mapTo(domIter1->lowerCorner(), domIter1->upperCorner(), quNodes1, quWeights1);
        quRule2.mapTo(domIter2->lowerCorner(), domIter2->upperCorner(), quNodes2, quWeights2);

        // Push forward the quad-points to the physical domain
        geoEval1->evaluateAt(quNodes1);
        geoEval2->evaluateAt(quNodes2);
        
        // Evaluate basis functions and their first derivatives
        // assuming numActive1=numActive2
        B1.compute(quNodes1, bdata1);
        mapperParam.localToGlobal( bdata1.actives, patch1, bdata1.actives);
        const gsMatrix<T> & ev1  = bdata1.values[0];
        const index_t numActive = bdata1.actives.rows();

        B2.compute(quNodes2, bdata2);
        mapperParam.localToGlobal( bdata2.actives, patch2, bdata2.actives);
        const gsMatrix<T> & ev2  = bdata2.values[0];
        
        B11.setZero(numActive, numActive); B22.setZero(numActive, numActive); 
        B12.setZero(numActive, numActive); B21.setZero(numActive, numActive);
        E11.setZero(numActive, numActive); E22.setZero(numActive, numActive); 
        E12.setZero(numActive, numActive); E21.setZero(numActive, numActive);

        // assuming domIter1->quNodes.cols() == domIter2->quNodes.cols()
        for (index_t k=0; k!= bdata1.allValues().cols(); ++k)
        {            
            // Compute first fund. form
            const gsMatrix<T> & jac1 = geoEval1->jacobians().block(0, k*d, d+1, d);
            const gsMatrix<T> & jac2 = geoEval2->jacobians().block(0, k*d, d+1, d);
            const gsMatrix<T> FirstFund1 = jac1.transpose()*jac1;
            const gsMatrix<T> FirstFund2 = jac2.transpose()*jac2;
            
            // Transform the basis gradients
            grads_k_1 = bdata1.values[1].col(k);
            grads_k_2 = bdata2.values[1].col(k);
            grads_k_1.resize( 2, numActive );
            grads_k_2.resize( 2, numActive );
            grads_k_1 = jac1 * FirstFund1.inverse() * grads_k_1;
            grads_k_2 = jac2 * FirstFund2.inverse() * grads_k_2;// param. assumed conforming    
            // *** Compute the co-normal vector 
            const gsVector<T,3> tangent  =  orient * jac1.template block<3, 1>(0,!dir1);
            // Check for zero tangent
            if ( tangent.squaredNorm() < 1e-10 ) 
            {
                gsDebug<< "Skip "<< geoEval1->values().col(k).transpose() <<"\n";
                continue;
            }

            // parametrization normal side 1
            gsVector<T> surfNormal;
            geoEval1->normal(k, surfNormal);
            //equiv: jac1.template block<3, 1>(0,0).cross(
            //       jac1.template block<3, 1>(0,1) ).normalized();
       
            unormal = tangent.cross( surfNormal.template head<3>() ).normalized() ;
            //gsDebugVar( unormal.transpose() );

            // Integral transformation and quadarature weight
            const T fff = quWeights1[k] * tangent.norm(); //unormal.norm();
            //gsDebugVar( tangent2.cross( u_surf_normal ).normalized().transpose() ) ;


            // Compute element matrices
            const gsMatrix<T> & val1 = ev1.col(k);
            const gsMatrix<T> & val2 = ev2.col(k);
            const T c1 = fff * T(0.5);
            N1.noalias()   = unormal.transpose() * grads_k_1;
            N2.noalias()   = unormal.transpose() * grads_k_2;
            B11.noalias() += c1 * ( val1 * N1 );
            B22.noalias() -= c1 * ( val2 * N2 );
            B12.noalias() += c1 * ( val1 * N2 );
            B21.noalias() -= c1 * ( val2 * N1 );
            const T c2 = fff * mu;
            E11.noalias() += c2 * ( val1 * val1.transpose() );
            E22.noalias() += c2 * ( val2 * val2.transpose() );
            E12.noalias() += c2 * ( val1 * val2.transpose() );
            E21.noalias() += c2 * ( val2 * val1.transpose() );
        }
        
        //gsDebugVar( N1 );
        
        // Push element contributions to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const index_t jj1 = bdata1.actives.at(j); // N1_j
            const index_t jj2 = bdata2.actives.at(j); // N2_j
            for (index_t i=0; i!=numActive; ++i)
            {
                // assuming symmetric problem
                const index_t  ii1 = bdata1.actives.at(i); // N1_i
                const index_t  ii2 = bdata2.actives.at(i); // N2_i
                
                if ( jj1 <= ii1 )
                    lhs( ii1, jj1 ) -=  B11(i,j) + B11(j,i) - E11(i,j);
                if ( jj2 <= ii2 ) 
                    lhs( ii2, jj2 ) -=  B22(i,j) + B22(j,i) - E22(i,j);
                if ( jj2 <= ii1 ) 
                    lhs( ii1, jj2)  -=  B12(i,j) + B21(j,i) + E12(i,j);
                if ( jj1 <= ii2 ) 
                    lhs( ii2, jj1)  -=  B21(i,j) + B12(j,i) + E21(i,j);
            }
        }
        
        domIter2->next();
    }
}

template<class T> void 
gsGaussSurfaceAssembler<T>::applyBoundary( const gsBasis<T>   & B,
                                           const boundary_condition<T> & bc,
                                           const gsDofMapper& mapper3,
                                           gsSparseSystem0<T> & system )
{    
    switch ( bc.type() )
    {
    case condition_type::dirichlet:
        boundaryNitsche(B, bc.patch(), bc.side(), *bc.function(), system);
        break;
    case condition_type::neumann:
        boundaryNeumann(B, bc.patch(), bc.side(), *bc.function(), system);
        break;
    default:
        gsWarn<<"Unknown boundary condition.\n";
    }
}


template<class T> void 
gsGaussSurfaceAssembler<T>::boundaryNitsche( const gsBasis<T> & B,
                                      const index_t patch,
                                      const boxSide s,
                                      const gsFunction<T> & f,
                                      gsSparseSystem0<T> & system )
{
    gsSparseMatrix<T> & lhs = *system.matrix();
    gsVector<T>       & rhs = *system.rhs();
    const int d   = this->m_geometry->parDim() ;

    const T mu = this->getMu(B);
    const int dir      = s.direction();

    //gsDebug<<"Nitsche boundary: side="<< s<<", patch="<<patch<<"(mu="<<mu<<").\n";

    // Quadrature for boundary integral: we fix coordinates for
    // direction = dir to be the fixed coordinate on the edge/face
    gsVector<index_t> bd_intNodes = getNumIntNodesForSide( B, s.direction() );
    
    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_VALUE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM, this->geometry()));
    const T orient = sideOrientation(s);

    // Temporaries
    gsMatrix<T> fev, grads_k;
    gsVector<T> unormal(d);

    // Local matrix and load vector
    gsMatrix<T> LM;
    gsVector<T> LB;

    typename gsDomainIterator<T>::uPtr domIter = B.makeDomainIterator(s);

    gsGaussRule<T> quRule(bd_intNodes);
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;
    gsFuncData<T> bdata(NEED_DERIV|NEED_ACTIVE|SAME_ELEMENT);

    // iterate over all boundary grid cells
    for (; domIter->good(); domIter->next())
    {
        // Compute the quadrature rule (nodes and weights)
        quRule.mapTo(domIter->lowerCorner(), domIter->upperCorner(), quNodes, quWeights);

        // Evaluate the geometry on the Gauss points
        geoEval->evaluateAt(quNodes);

        // Evaluate basis functions and their first derivatives
        B.compute(quNodes, bdata);
        mapper.localToGlobal(bdata.actives, patch, bdata.actives);
        const index_t numActive = bdata.actives.rows();
        const gsMatrix<T> & ev  = bdata.values[0];

        // Evaluate the Dirichlet data
        f.eval_into(geoEval->values(), fev);

        LM.setZero(numActive, numActive);
        LB.setZero(numActive);

        for (index_t k=0; k!= bdata.allValues().cols(); ++k) // For all quadrature points
        {
            //  *** Compute first fund. form
            const gsMatrix<T> & jac = geoEval->jacobians().block(0, k*d, d+1, d);
            const gsMatrix<T> FirstFund = jac.transpose()*jac;
                        
            //  *** Transform the basis gradients
            grads_k = bdata.values[1].col(k);
            grads_k.resize( 2, numActive );
            grads_k = jac * FirstFund.inverse()*grads_k;

            // *** Compute the co-normal vector 
            const gsVector<T,3> tangent  = orient * jac.template block<3, 1>(0,!dir);
            // Check for zero tangent
            if ( tangent.squaredNorm() < 1e-10 ) 
            {
                gsDebug<< "Skip "<< geoEval->values().col(k).transpose() <<"\n";
                continue;
            }

            gsVector<T> surfNormal;
            geoEval->normal(k, surfNormal);
            surfNormal.normalize();
            //equiv: jac1.template block<3, 1>(0,0).cross(
            //       jac1.template block<3, 1>(0,1) ).normalized();
            unormal = tangent.cross( surfNormal.template head<3>() ).normalized() ;

            // Integral transformation and quadarature weight
            const T fff = quWeights[k] * tangent.norm();
           
            // Sum up quadrature point evaluations
            LB.noalias() += fff * fev(0,k) * ( grads_k.transpose() * unormal - mu * ev.col(k) );
            LM.noalias() += fff * ( ev.col(k) * unormal.transpose() * grads_k
                                +  (ev.col(k) * unormal.transpose() * grads_k).transpose()
                                -  mu * ev.col(k) * ev.col(k).transpose() );
        }

        // Push element contribution to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const unsigned jj = bdata.actives.at(j);
            rhs[jj] -= LB[j];
            for (index_t i=0; i!=numActive; ++i)
            {
                const unsigned ii = bdata.actives(i);
                if ( jj <= ii ) // assuming symmetric problem
                    lhs( ii, jj ) -= LM(i,j);
            }
        }
    }
}


template<class T> void 
gsGaussSurfaceAssembler<T>::boundaryNeumann( const gsBasis<T> & B,
                                      const index_t patch,
                                      const boxSide s,
                                      const gsFunction<T> & f,
                                      gsSparseSystem0<T> & system )
{  
    //gsDebug<<"Neumann boundary: side="<< s<<", patch="<<patch<<"\n";
    const int d   = this->m_geometry->parDim() ;
    gsVector<T> & rhs = *system.rhs();

    // Quadrature for boundary integral: we fix coordinates for
    // direction = dir to be the fixed coordinate on the edge/face
    gsVector<index_t> bd_intNodes = getNumIntNodesForSide( B, s.direction() );

    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_VALUE | NEED_JACOBIAN, this->geometry()));

    // Temporaries
    gsMatrix<T> fev;
    gsVector<T> localRhs, unormal(d);

    gsGaussRule<T> quRule(bd_intNodes);
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;
    gsFuncData<T> bdata(NEED_VALUE|NEED_ACTIVE|SAME_ELEMENT);
    
    // iterate over all boundary grid cells
    for (typename gsDomainIterator<T>::uPtr domIter = B.makeDomainIterator(s); 
         domIter->good(); domIter->next())
    {
        // Compute the quadrature rule (nodes and weights)
        quRule.mapTo(domIter->lowerCorner(), domIter->upperCorner(), quNodes, quWeights);

        // Evaluate the geometry on the Gauss points
        geoEval->evaluateAt(quNodes);

        // Evaluate the basis functions
        mapper.localToGlobal(bdata.actives, patch, bdata.actives);
        const index_t numActive = bdata.actives.rows();
        B.compute(quNodes, bdata);

        // Evaluate the Neumann data
        f.eval_into(geoEval->values(), fev);

        localRhs.setZero(numActive);

        for (index_t k=0; k!= bdata.allValues().cols(); ++k) // For all quadrature points
        {
            // Compute the outer normal vector on the side
            geoEval->outerNormal(k, s, unormal);

            // Sum up quadrature evaluations
            const T fff = quWeights[k] * fev(0,k) *unormal.norm();
            localRhs.noalias() += fff * bdata.values[0].col(k);
        }
        
        // Push element contribution to the global load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            // convert local dof index to global dof index
            const unsigned jj = bdata.actives.at(j);
            rhs[jj] += localRhs[j];
        }
    }

}




template <class T>
void gsAssembler0<T>::localToGlobal(const gsMatrix<T>& localStiffness, 
                                   const gsMatrix<index_t>& localDofs,
                                   gsSparseMatrix<T>& K, 
                                   bool symmetric)
{
    const int numActive = localDofs.rows();

    for (index_t i = 0; i < numActive; ++i)
    {
        const int ii = localDofs(i,0);
        for (index_t j = 0; j < numActive; ++j)
        {
            const int jj = localDofs(j,0);
            // if matrix is symmetric, store only lower triangular part
            if (!symmetric || jj <= ii)
                K.coeffRef(ii, jj) += localStiffness(i, j);
        }
    }
}


template <class T>
void gsAssembler0<T>::localToGlobal_withBC(const gsMatrix<T>& localStiffness,
                                          const gsVector<T>& localRhs,
                                          const gsMatrix<T>& dirbc,
                                          const gsDofMapper& mapper,
                                          const gsVector<index_t>& loc2glob,
                                          gsSparseMatrix<T>& K,
                                          gsVector<T>& f,
                                          bool symmetric)
{
    const int numActive = loc2glob.size();

    for (index_t i=0; i < numActive; ++i)
    {
        const int ii = loc2glob[i];
        if ( mapper.is_free_index(ii) )
        {
            f[ii] += localRhs[i];

            for (index_t j=0; j < numActive; ++j)
            {
                const int jj = loc2glob[j];
                if ( mapper.is_free_index(jj) )
                {
                    // if matrix is symmetric, store only lower triangular part
                    if (!symmetric || jj <= ii)
                        K.coeffRef(ii, jj) += localStiffness(i, j);
                }
                else if ( mapper.is_boundary_index(jj) )        // Dirichlet boundary condition?
                {
                    f[ii] -= dirbc( mapper.global_to_bindex(jj) ) * localStiffness(i, j);
                }
            }
        }
    }
}

/// Add contributions from local to global matrix, eliminating Dirichlet BCs
template <class T>
void gsAssembler0<T>::localToGlobal_withBC(const gsMatrix<T>& localMatrix,
                                          const gsDofMapper& mapper,
                                          const gsVector<index_t>& loc2glob,
                                          gsSparseMatrix<T>& K,
                                          bool symmetric)
{
    const int numActive = loc2glob.size();

    for (index_t i=0; i < numActive; ++i)
    {
        const int ii = loc2glob[i];
        if ( !mapper.is_free_index(ii) )
            continue;

        for (index_t j=0; j < numActive; ++j)
        {
            const int jj = loc2glob[j];
            if ( !mapper.is_free_index(jj) )
                continue;

            // if matrix is symmetric, store only lower triangular part
            if (!symmetric || jj <= ii)
                K.coeffRef(ii, jj) += localMatrix(i, j);
        }
    }
}


template <class T>
gsSparseMatrix<T> * gsAssembler0<T>::initMatrix(int nDofs, const gsBasis<T>& basis, bool symmetric)
{
    gsSparseMatrix<T> * K = new gsSparseMatrix<T>(nDofs, nDofs);

    // estimate max nz per row (only valid for tensor product bases right now)
    unsigned nzRowsPerCol = 1;
    for (int i = 0; i < basis.dim(); ++i)
        nzRowsPerCol *= 2 * basis.degree(i) + 1;    // diagonal entry plus p interactions on each side

    // Unfortunately, this optimization doesn't work in general because Neumann/interface
    // dofs are not ordered the same way as interior dofs, thus they can have more nzs.
    //if (symmetric)
    //    nzRowsPerCol = (nzRowsPerCol + 1) / 2;
    if (nDofs > 0)
    K->reserve( gsVector<index_t>::Constant(nDofs, nzRowsPerCol) );
    return K;
}

template <class T>
gsSparseSystem0<T> gsAssembler0<T>::initLinearSystem(int nDofs, const gsBasis<T>& basis, bool symmetric)
{
    gsSparseMatrix<T> * K = initMatrix(nDofs, basis, symmetric);
    gsVector<T> * b = new gsVector<T>;
    b->setZero( nDofs );
    return gsSparseSystem0<T>( K, b );
}


} // namespace gismo

