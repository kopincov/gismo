/** @file gsLocalErrorEstimation.cpp

    @brief Simple testfile to estimate the error on an element.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;

template <class T, bool paramCoef = false>
class gsVisitorPoissonResidual
{
public:

    /** \brief Constructor for gsVisitorPoisson.
     */
    gsVisitorPoissonResidual(const gsPde<T> & pde,const gsField<T>& uh):m_uh(uh)
    {
        pde_ptr = static_cast<const gsPoissonPde<T>*>(&pde);
    }

    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule)
    {
        // Grab right-hand side for current patch
        rhs_ptr = &pde_ptr->rhs()->piece(patchIndex);

        // Setup Quadrature
        rule = gsGaussRule<T>(basis, options);// harmless slicing occurs here

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>    & basis,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T>   & quNodes)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);
        numActive = actives.rows();

        // Evaluate basis functions on element
        basis.evalAllDers_into(md.points, 1, basisData);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo.computeMap(md);

        // Evaluate right-hand side at the geometry points paramCoef
        // specifies whether the right hand side function should be
        // evaluated in parametric(true) or physical (false)
        rhs_ptr->eval_into((paramCoef ? md.points : md.values[0]), rhsVals);
        m_uh.function(0).deriv_into(md.values[0], gradUhVals);

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides
    }

    inline void assemble(gsDomainIterator<T>    & /*element*/,
                         gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];
        gsMatrix<T> & bGrads = basisData[1];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);

            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads, physGrad);

            localRhs.noalias() += weight * ( bVals.col(k) * rhsVals.col(k).transpose() - physGrad.transpose()*gradUhVals.col(k) ) ;
            localMat.noalias() += weight * (physGrad.transpose() * physGrad);
        }
    }

    /* old gsAssembler2.h */
//    void initialize(const gsBasis<T> & basis,
//                    const index_t patchIndex,
//                    const gsOptionList & options,
//                    gsQuadRule<T>    & rule,
//                    unsigned         & evFlags )
//    {
//        // Grab right-hand side for current patch
//        rhs_ptr = &pde_ptr->rhs()->piece(patchIndex);
//
//        // Setup Quadrature
//        rule = gsGaussRule<T>(basis, options);// harmless slicing occurs here
//
//        // Set Geometry evaluation flags
//        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
//    }
//
//    // Evaluate on element.
//    inline void evaluate(gsBasis<T> const       & basis,
//                         gsGeometryEvaluator<T> & geoEval,
//                         gsMatrix<T> const      & quNodes)
//    {
//        // Compute the active basis functions
//        // Assumes actives are the same for all quadrature points on the elements
//        basis.active_into(quNodes.col(0), actives);
//        numActive = actives.rows();
//
//        // Evaluate basis functions on element
//        basis.evalAllDers_into( quNodes, 1, basisData);
//
//        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
//        geoEval.evaluateAt(quNodes);// is this generic ??
//
//        // Evaluate right-hand side at the geometry points paramCoef
//        // specifies whether the right hand side function should be
//        // evaluated in parametric(true) or physical (false)
//        rhs_ptr->eval_into( (paramCoef ?  quNodes :  geoEval.values() ), rhsVals );
//        m_uh.function(0).deriv_into(geoEval.values(),gradUhVals);
//
//        // Initialize local matrix/rhs
//        localMat.setZero(numActive, numActive      );
//        localRhs.setZero(numActive, rhsVals.rows() );//multiple right-hand sides
//    }
//
//    inline void assemble(gsDomainIterator<T>    & element,
//                         gsGeometryEvaluator<T> & geoEval,
//                         gsVector<T> const      & quWeights)
//    {
//        gsMatrix<T> & bVals  = basisData[0];
//        gsMatrix<T> & bGrads = basisData[1];
//
//        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
//        {
//            // Multiply weight by the geometry measure
//            const T weight = quWeights[k] * geoEval.measure(k);
//
//            // Compute physical gradients at k as a Dim x NumActive matrix
//            geoEval.transformGradients(k, bGrads, physGrad);
//
//            localRhs.noalias() += weight * ( bVals.col(k) * rhsVals.col(k).transpose() - physGrad.transpose()*gradUhVals.col(k) ) ;
//            localMat.noalias() += weight * (physGrad.transpose() * physGrad);
//        }
//    }

    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives, eliminatedDofs.front(), 0, 0);
    }

protected:
    // Pointer to the pde data
    const gsPoissonPde<T> * pde_ptr;
    const gsField<T>& m_uh;

protected:
    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T>        physGrad;
    gsMatrix<index_t> actives;
    index_t numActive;

protected:
    // Right hand side ptr for current patch
    const gsFunction<T> * rhs_ptr;

    // Local values of the right hand side
    gsMatrix<T> rhsVals;
    gsMatrix<T> gradUhVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    gsMapData<T> md;
};


template<typename T>
class gsResidualPoissonAssembler : public gsPoissonAssembler<T>
{
    typedef gsPoissonAssembler<T> Base;
public:
    gsResidualPoissonAssembler( gsMultiPatch<T> const         & patches,
                                gsMultiBasis<T> const         & basis,
                                gsBoundaryConditions<T> const & bconditions,
                                const gsFunction<T>           & rhs,
                                const gsField<T>          & uh) : gsPoissonAssembler<T>(patches,basis,bconditions,rhs),m_uh(uh)
    {

    }

    virtual void assemble()
    {
        GISMO_ASSERT(m_system.initialized(),
                     "Sparse system is not initialized, call initialize() or refresh()");

        // Reserve sparse system
        m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        Base::computeDirichletDofs();

        // Clean the sparse system
        // m_system.setZero(); //<< this call leads to a quite significant performance degrade!

        // Assemble volume integrals
        gsVisitorPoissonResidual<T> vis(*m_pde_ptr,m_uh);
        Base::template push(vis);

        // Enforce Neumann boundary conditions
        //       Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );

        //     const int dirStr = m_options.getInt("DirichletStrategy");

        // Assembly is done, compress the matrix
        Base::finalize();
    }

private:
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
    const gsField<T>& m_uh;
};


template<typename T>
class gsElementErrorPlotter : public gsFunction<T>
{
public:
    gsElementErrorPlotter(const gsBasis<T>& mp, const std::vector<T>& errors ) : m_mp(mp),m_errors(errors)
    {

    }

    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& res) const
    {
        // Initialize domain element iterator -- using unknown 0
        res.setZero(1,u.cols());
        for(index_t i=0; i<u.cols();++i)
        {
            int iter =0;
            // Start iteration over elements

            typename gsBasis<T>::domainIter domIt = m_mp.makeDomainIterator();
            for (; domIt->good(); domIt->next() )
            {
                 bool flag = true;
                const gsVector<T>& low = domIt->lowerCorner();
                const gsVector<T>& upp = domIt->upperCorner();


                for(int d=0; d<domainDim();++d )
                {
                    if(low(d)> u(d,i) || u(d,i) > upp(d))
                    {
                        flag = false;
                        break;
                    }
                }
                if(flag)
                {
                     res(0,i) = m_errors[iter];
                     break;
                }
                iter++;
            }
        }
    }

    short_t domainDim() const { return m_mp.dim();}

private:
    const gsBasis<T>& m_mp;
    const std::vector<T>& m_errors;
};


int main(int argc, char *argv[])
{
    bool plot = false;

    gsCmdLine cmd("Calculates the elementwise residuum");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Right hand side and solution
    gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",2);
    gsFunctionExpr<> g("sin(pi*x*1)*sin(pi*y*2)",2);

    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<"\n\n";

    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches(*gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(1)));
    gsInfo << "The domain is a "<< patches <<"\n";

    gsBoundaryConditions<> bcInfo;
    for (gsMultiPatch<>::const_biterator bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
        bcInfo.addCondition( *bit, condition_type::dirichlet, &g );

    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );


    //Do two refinements and set degree to 3 (initial basis has p = 1)
    int numRefine  = 2;
    int numElevate = 1;

    for (int i = 0; i < numRefine; ++i)
        refine_bases.uniformRefine();

    int max_tmp = refine_bases.minCwiseDegree();
    max_tmp += numElevate;
    refine_bases.setDegree(max_tmp);
    refine_bases.basis(0).print(gsInfo);
    gsInfo<<"Degree of initial basis is set to "<<max_tmp<<"\n";

    gsPoissonAssembler<real_t> PoissonAssembler(patches,refine_bases,bcInfo,f,dirichlet::elimination);

    // Generate system matrix and load vector
    gsInfo<< "Assembling...\n";
    PoissonAssembler.assemble();
    gsInfo << "Assembled a system (matrix and load vector) with "<< PoissonAssembler.numDofs() << " dofs.\n";

    // Setup the direct Solver
    gsInfo << "Solving...\n";
    gsSparseSolver<>::SimplicialLDLT solver( PoissonAssembler.matrix() );
    gsMatrix<> u_vec = solver.solve( PoissonAssembler.rhs() );


    // Construct the solution as a scalar field
    gsMultiPatch<> mpsoluh;
    PoissonAssembler.constructSolution(u_vec, mpsoluh);
    gsField<> u_h( PoissonAssembler.patches(), mpsoluh);


    //Construct the finer basis
    gsMultiBasis<> fineBasis(refine_bases);
    fineBasis.degreeIncrease(1);
    fineBasis.reduceContinuity(fineBasis.degree()-2);

    fineBasis.basis(0).print(gsInfo);

    gsResidualPoissonAssembler<real_t> PoissonAssemblerError(patches,fineBasis,bcInfo,f,u_h);
    PoissonAssemblerError.options().setInt("DirichletValues",dirichlet::homogeneous);
    PoissonAssemblerError.assemble();
    // gsInfo<<PoissonAssemblerError.matrix().toDense()<<"\n";
    solver.compute(PoissonAssemblerError.matrix());
    gsMatrix<> e_vec = solver.solve( PoissonAssemblerError.rhs() );

    gsMultiPatch<> mpsoleh;
    PoissonAssemblerError.constructSolution(e_vec, mpsoleh);
    gsField<> e_h( PoissonAssemblerError.patches(), mpsoleh);

    gsSeminormH1<real_t> seminorm(e_h);
    seminorm.compute(true);
    std::vector<real_t> elementNorms = seminorm.elementNorms();


    gsElementErrorPlotter<real_t> err_eh(fineBasis.basis(0),elementNorms);


    gsTensorBSplineBasis<2,real_t>::uPtr b(memory::convert_ptr<gsTensorBSplineBasis<2,real_t> >(fineBasis.basis(0).clone()));
    // b.insertKnot(0.3,0,1); //xi,direction,multiplicity
    //  b.insertKnot(0.3,1,1);
    //... insert more

    b->print(gsInfo);

    gsResidualPoissonAssembler<real_t> PoissonAssemblerErrorNew(patches,*b,bcInfo,f,u_h);
    PoissonAssemblerErrorNew.options().setInt("DirichletValues",dirichlet::homogeneous);
    PoissonAssemblerErrorNew.assemble();
    solver.compute(PoissonAssemblerErrorNew.matrix());
    gsMatrix<> e_vec_new = solver.solve( PoissonAssemblerErrorNew.rhs() );

    gsMultiPatch<> mpsoleh_new;
    PoissonAssemblerErrorNew.constructSolution(e_vec_new, mpsoleh_new);
    gsField<> e_h_new( PoissonAssemblerErrorNew.patches(), mpsoleh_new);

    gsSeminormH1<real_t> seminorm_new(e_h_new);
    seminorm_new.compute(true);
    std::vector<real_t> elementNorms_new = seminorm_new.elementNorms();

    gsElementErrorPlotter<real_t> err_eh_new(*b,elementNorms_new);

    if (plot)
    {
        //! [Plot in Paraview]
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>(u_h, "solution_uh", 1000,true);
        gsWriteParaview<>(e_h, "error_eh", 1000,true);
        gsWriteParaview<>(e_h_new, "errorFine_eh", 1000,true);
        const gsField<> exact( PoissonAssembler.patches(), g, false );
        const gsField<> elemError_eh( PoissonAssembler.patches(), err_eh, false );
        const gsField<> elemError_eh_new( PoissonAssembler.patches(), err_eh_new, false );
        gsWriteParaview<>( exact, "solution_exact", 1000);
        gsWriteParaview<>( elemError_eh, "error_elem_eh", 1000);
        gsWriteParaview<>( elemError_eh_new, "errorFine_elem_eh", 1000);

        // Run paraview
        return system("paraview solution_uh.pvd &");
        //! [Plot in Paraview]
    }
    else
    {
        gsInfo<<"Quitting.. No output created, re-run with --plot to get a ParaView "
                "file containing Plotting image data.\n";
        return 0;
    }


}// end main
