/** @file gsOptElectricMotor.hpp

    @brief Provides implementation of an optimization problem.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): R. Schneckenleitner, A. Mantzaflaris
*/

#include <gsIO/gsFileManager.h>
#include <gsIO/gsWriteParaview.h>
#include <gsCore/gsField.h>
#include <gsModeling/gsCoonsPatch.h>
#include <gsModeling/gsSpringPatch.h>
#include <gsModeling/gsCrossApPatch.h>

// Dense jacobian too much for MUMPS solver
#define JAC_STRUCTURE

namespace gismo
{

void createInterfaceTopology(std::vector<std::pair<unsigned, unsigned > > & interfaces);
gsMultiPatch<> createInterfacePatch(const gsMultiPatch<> & mp);


template <typename T>
gsOptElectricMotor<T>::gsOptElectricMotor(gsMultiPatch<T> domain, gsMultiBasis<T> mbasis, gsPiecewiseFunction<T> magRel, gsPiecewiseFunction<T> magnetization, gsFunctionExpr<T> reference, gsPiecewiseFunction<T> penalization,
                                          gsOptionList opt, gsBoundaryConditions<> bcInfo, std::vector<size_t> desDomain, std::vector<size_t> desNeighbours, unsigned maxIter) :
                                                                                                               m_tmpPatches(domain),
                                                                                                               m_basis(give(mbasis)),
                                                                                                               //m_interfacepatch(createInterfacePatch(domain)), // createInterfacePatch(domain)
                                                                                                               m_reluctivity(magRel),
                                                                                                               m_reference(give(reference)),
                                                                                                               m_rhs(magnetization),
                                                                                                               m_penalization(penalization),
                                                                                                               m_opt(give(opt)),
                                                                                                               m_bc(give(bcInfo)),
                                                                                                               m_desDomain(give(desDomain)),
                                                                                                               m_desNeighbours(give(desNeighbours)),
                                                                                                               m_maxIterations(maxIter)

{
    // Set the default tolerance
    m_tolerance = T(1e-6);

    // Resize the vector to the number of patches
    m_colPoints.resize(domain.nPatches());
    m_corner.resize(domain.nPatches());
    m_patchOrientation.resize(domain.nPatches(), -1);
    //m_jacSolver.resize(domain.nPatches());

    // Define the basis for the determinant of the jacobian
    //m_jacBasis = *dynamic_cast<gsMultiBasis<T>* >(&m_basis.clone());
    m_jacBasis = *m_basis.clone();

    for(size_t np = 0; np < domain.nPatches(); np++) {
        m_corner[np].resize(1<<domain.dim());
        m_jacBasis.basis(np).reduceContinuity(1);

        for (short_t k = 0; k < m_tmpPatches.dim(); k++)
            m_jacBasis.basis(np).component(k).degreeElevate((m_tmpPatches.dim() - 1) * m_basis.degree(k) - 1);
    }

    //m_mapper = gsDofMapper(m_basis);
    ///===================================NEW=======================================================
    ///SIZE Reduction
    m_basis.getMapper(iFace::glue, m_mapper, false);
    m_numDesignVars = 0;

    for(index_t np = 0; np < static_cast<index_t>(domain.nPatches()); np++)
    {
        gsMatrix<T> tmpMat = domain[np].coefs();

        if( !(std::find(m_desDomain.begin(), m_desDomain.end(), np) != m_desDomain.end()) &&
            !(std::find(m_desNeighbours.begin(), m_desNeighbours.end(), np) != m_desNeighbours.end()) ) {
            for (index_t dof = 0; dof < tmpMat.rows(); dof++)
                m_mapper.eliminateDof(dof, np);
        }
        else
        {
            if((std::find(m_desDomain.begin(), m_desDomain.end(), np) != m_desDomain.end())) {
                gsMatrix<index_t> bMat;
                size_t start = 0;

                for (std::vector<size_t>::iterator p2 = m_desNeighbours.begin(); p2 != m_desNeighbours.end(); p2++) {
                    const boundaryInterface *bI = domain.findInterface(np, *p2);

                    if (bI != NULL) {
                        gsMatrix<index_t> singleBnd;
                        if(bI->first().patch == np){
                            singleBnd = m_basis.basis(np).boundary(bI->first()); //gsInfo << "Dirmap: " << bI->dirMap(bI->first()) << "\n";
                            }
                        else{
                            singleBnd = m_basis.basis(np).boundary(bI->second()); //gsInfo << "Dirmap: " << bI->dirMap(bI->second()) << "\n";
                            }

                        //gsInfo << "patch: " << np << "\n";
                        //gsInfo << "boundary: " << singleBnd << "\n";

                        bMat.conservativeResize(start + singleBnd.rows(), 1);
                        bMat.block(start, 0, singleBnd.rows(), singleBnd.cols()) = singleBnd;
                        m_numDesignVars += 2 * singleBnd.rows();
                        start += singleBnd.rows();
                        //gsInfo << "bMat: " << bMat << "\n";
                    }
                }

                for (index_t r = 0; r < tmpMat.rows(); r++)
                    if ((unsigned(r) != bMat.array()).all()) {
                            //gsInfo << "Eliminated: " << r << "\n";
                            m_mapper.eliminateDof(r, np);
                    }

                // Do not count the corner coefficients as design variables
                m_numDesignVars -= 2*m_corner[np].size(); // TODO: Remove this line when T-junctions are allowed
            }
            else
            {
                if ((std::find(m_desNeighbours.begin(), m_desNeighbours.end(), np) != m_desNeighbours.end()))
                {
                    gsMatrix<index_t> bMat;
                    size_t start = 0;

                    for (std::vector<size_t>::iterator p2 = m_desDomain.begin(); p2 != m_desDomain.end(); p2++)
                    {
                        const boundaryInterface *bI = domain.findInterface(np, *p2);

                        if (bI != NULL)
                        {
                            gsMatrix<index_t> singleBnd;
                            if(bI->first().patch == np){
                                singleBnd = m_basis.basis(np).boundary(bI->first()); //gsInfo << "Dirmap2: " << bI->dirMap(bI->first()) << "\n";
                                }
                            else{
                                singleBnd = m_basis.basis(np).boundary(bI->second()); //gsInfo << "Dirmap2: " << bI->dirMap(bI->second()) << "\n";
                                }



                            bMat.conservativeResize(start + singleBnd.rows(), 1);
                            bMat.block(start, 0, singleBnd.rows(), singleBnd.cols()) = singleBnd;
                            start += singleBnd.rows();
                        }
                    }

                    for (index_t r = 0; r < tmpMat.rows(); r++)
                        if ( (r != bMat.array()).all() ) {
                            m_mapper.eliminateDof(r, np);
                        }
                }
            }
        }
    }

    m_mapper.finalize();

    GISMO_ASSERT(2*m_mapper.freeSize() == m_numDesignVars, "Freesize does not correspond to no. of design vars, free dofs: " << m_mapper.freeSize() << " and design vars: " <<m_numDesignVars);

    // Once save the preImages
    m_curDesign.resize(2*m_mapper.freeSize(), 1);
    std::vector<std::pair<index_t, index_t> > preIm;

    for(index_t i = 0; i < m_mapper.freeSize(); i++)
    {
        m_mapper.preImage(i, preIm);
        m_preImages.push_back(preIm);

        GISMO_ASSERT(preIm.size() == 2, "Dof is a corner or not an interface");

        //gsInfo << "The patches " << preIm[0].first << " and the dof: " << preIm[0].second <<"\n";
        //gsInfo << "The boundary " << domain[preIm[0].first].basis().allBoundary() <<"\n";

        // no loop since the coefficients must be the same
        for(index_t d = 0; d < m_tmpPatches.dim(); d++)
            m_curDesign(i + d * m_mapper.freeSize(), 0) = domain[preIm[0].first].coefs()(preIm[0].second, d);


    }

    //gsInfo << "The dof is free? " << m_mapper.is_free(0, 60) <<"\n"; //part of the design domain

    ///===================================END NEW=======================================================

    // Initialize m_colPoints and specify the number of constraints
    // ... and save the corners
    for(size_t np = 0; np < domain.nPatches(); np++) {
        m_jacBasis.basis(np).anchors_into(m_colPoints[np]);

        for (index_t i = 0; i!= m_corner[np].size(); ++i)
            m_corner[np][i] = m_jacBasis.basis(np).functionAtCorner(i);
    }

    // Set the current design
    //m_curDesign = multipatchToVector(domain);

    /// START change
    /// Old version (Big version)
    //m_numDesignVars = m_curDesign.size();

    // Design bounds
    m_desLowerBounds.setConstant(m_numDesignVars, -1.0e19);
    m_desUpperBounds.setConstant(m_numDesignVars,  1.0e19);

    // number of constraints
    m_numConstraints = 0;
    //m_numConJacNonZero = 0; // Only needed here without contraints

    // constraints for the design domains
    for(std::vector<size_t>::iterator it = m_desDomain.begin(); it != m_desDomain.end(); it++)
    {
        gsVector<T> firstBound = gsVector<T>::Constant(m_jacBasis.basis(*it).size(), 1.0e-8);
        gsVector<T> secondBound = gsVector<T>::Constant(m_jacBasis.basis(*it).size(), 1.0e19);

        bool isPositive = computeInitialjacobianDetCoefs(*it);
        m_patchOrientation.push_back(isPositive);

        m_conLowerBounds.conservativeResize(m_numConstraints + m_jacBasis.basis(*it).size());
        m_conUpperBounds.conservativeResize(m_numConstraints + m_jacBasis.basis(*it).size());

        // Constraint bounds
        if(isPositive)
        {
            m_conLowerBounds.segment(m_numConstraints, m_jacBasis.basis(*it).size()) = firstBound;
            m_conUpperBounds.segment(m_numConstraints, m_jacBasis.basis(*it).size()) = secondBound;
        }
        else
        {
            m_conLowerBounds.segment(m_numConstraints, m_jacBasis.basis(*it).size()) = T(-1.0) * secondBound;
            m_conUpperBounds.segment(m_numConstraints, m_jacBasis.basis(*it).size()) = T(-1.0) * firstBound;
        }
        m_numConstraints += m_jacBasis.basis(*it).size();

        // store the factorizations of the collocation matrices
        // cheaper than to compute it every time
        //gsSparseMatrix<T> colMat;
        //m_jacBasis.basis(*it).collocationMatrix(m_colPoints[*it], colMat);
        //m_jacSolver[*it].compute(colMat);
    }

    // constraints for the neighbours of the design domains
    for(std::vector<size_t>::iterator it = m_desNeighbours.begin(); it != m_desNeighbours.end(); it++)
    {
        gsVector<T> firstBound = gsVector<T>::Constant(m_jacBasis.basis(*it).size(), 1.0e-8);
        gsVector<T> secondBound = gsVector<T>::Constant(m_jacBasis.basis(*it).size(), 1.0e19);

        bool isPositive = computeInitialjacobianDetCoefs(*it);
        m_patchOrientation.push_back(isPositive);

        m_conLowerBounds.conservativeResize(m_numConstraints + m_jacBasis.basis(*it).size());
        m_conUpperBounds.conservativeResize(m_numConstraints + m_jacBasis.basis(*it).size());

        // Constraint bounds
        if(isPositive)
        {
            m_conLowerBounds.segment(m_numConstraints, m_jacBasis.basis(*it).size()) = firstBound;
            m_conUpperBounds.segment(m_numConstraints, m_jacBasis.basis(*it).size()) = secondBound;
        }
        else
        {
            m_conLowerBounds.segment(m_numConstraints, m_jacBasis.basis(*it).size()) = T(-1.0) * secondBound;
            m_conUpperBounds.segment(m_numConstraints, m_jacBasis.basis(*it).size()) = T(-1.0) * firstBound;
        }
        m_numConstraints += m_jacBasis.basis(*it).size(); // = m_colPoints[*it]


        // store the factorizations of the collocation matrices
        // cheaper than to compute it every time
        //gsSparseMatrix<T> colMat;
        //m_jacBasis.basis(*it).collocationMatrix(m_colPoints[*it], colMat);
        //m_jacSolver[*it].compute(colMat);
    }

    //m_conLowerBounds.setConstant(m_numConstraints,  1.0e-8);
    //m_conUpperBounds.setConstant(m_numConstraints,  1.0e19);


#ifdef JAC_STRUCTURE
    computeJacStructure();
#else
    gsOptProblem<T>::computeJacStructure();
#endif

    numIterations = 0;
    finalObjective = 0;
}


template <typename T>
T gsOptElectricMotor<T>::evalObj(const gsAsConstVector<T> & u) const
{
    // START: Update the temporary multipatch domain with the input u
    updateTempPatches(u);
    // END: Update the temporary multipatch domain with the input u

    // START: Construct an interface domain with the current design
//    gsMultiPatch<> InterfaceMPCopy(m_tmpPatches);
//    createInterfacePatch(InterfaceMPCopy);
//    m_interfacepatch = InterfaceMPCopy;
    // END: Construct an interface domain with the current design

    // Compute the coefficients of the solution Vector
    gsMagnetostaticPde<T> magPDE(m_tmpPatches, m_bc, m_rhs, m_reluctivity);
    gsMagnetostaticAssembler<T> magAss(magPDE, m_basis, m_opt, m_rhs);
    gsMatrix<T> state = solveState(magAss);

    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(m_basis);
    geometryMap G = A.getMap(m_tmpPatches);

    typedef gsExprAssembler<real_t>::space space;

    space phi = A.getSpace(m_basis);
    phi.setInterfaceCont(0);
    phi.addBc(m_bc.get("Dirichlet"));

    A.setOptions(m_opt);

    A.initSystem();

    solution u_f = A.getSolution(phi, state);

    gsExprEvaluator<T> J(A);
    variable u_df = J.getVariable(m_reference, G);
    J.integralInterface( (((grad(u_f) * jac(G).inv()) * (tv(G)/tv(G).norm())).tr() - u_df).sqr() * nv(G).norm(), (createInterfacePatch(m_tmpPatches)).interfaces());

    return J.value();
}

template <typename T>
void gsOptElectricMotor<T>::gradObj_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const
{
    //gsInfo << m_tmpPatches[60].coefs() << "\n now after: \n";
    updateTempPatches(u);
    //gsInfo << m_tmpPatches[60].coefs() << "\n";
    const index_t dimLen = u.size() / m_tmpPatches.dim();

    gsMultiPatch<T> solAsMultipatch;
    gsMultiPatch<T> solution_uh, solution_ph, shapeGradient;
    gsVector<T> sg_vector, red_vector;
    //gsMatrix<T> red_vector;

    // Compute the solution u of the state equation
    gsMagnetostaticPde<T> magPDE(m_tmpPatches, m_bc, m_rhs, m_reluctivity);
    solution_uh = solveStateForMultipatch(magPDE, m_opt, m_rhs);

    // Compute the solution p of the adjoint equation
    gsMagnetostaticAdjointPde<T> magAdjPDE(m_tmpPatches, m_bc, solution_uh, m_reference, m_reluctivity);
    solution_ph = solveAdjoint(magAdjPDE, m_tmpPatches, m_opt);

    // Compute the shape gradient
    gsMagnetostaticShapeDerivPde<T> shapeDerivPDE(m_tmpPatches, m_bc, m_penalization, solution_uh, solution_ph, m_reluctivity, m_rhs, 2);
    shapeGradient = calcShapeDerivative(shapeDerivPDE, m_opt, m_bc);

    // Start: Make a vector out of the multipatch
    sg_vector = multipatchToVector(shapeGradient);

    /// For the smaller, new approach
    red_vector.resize(m_numDesignVars);

    for(size_t freeDof = 0; freeDof < m_preImages.size(); freeDof++) // m_preImages.size() = m_numDesignVars
    {
        std::vector<std::pair<index_t, index_t> > preImage = m_preImages[freeDof];

        index_t patch = preImage[0].first;
        index_t dof = preImage[0].second;

        //gsInfo << "patch: " << patch << "\n";

        //gsMatrix<T> &cf = shapeGradient.patch(patch).coefs();
        gsMatrix<T> cf = shapeGradient.patch(patch).coefs();
        //gsInfo << "the row: " << cf.row(dof) << "\n";

        for (short_t d = 0; d < m_tmpPatches.dim(); d++)
            red_vector(freeDof + d * dimLen) = cf(dof, d);
        //red_vector.row(freeDof) = cf.row(dof);

    }

    //gsAsMatrix<> domain(sg_vector.data(), m_numDesignVars/m_tmpPatches.dim(), m_tmpPatches.dim());
    gsAsMatrix<> domain(red_vector.data(), red_vector.size()/m_tmpPatches.dim(), m_tmpPatches.dim());

    /// END

    // Store the domain as a vector
    // Compute M as a proper scaling factor for the shape gradient. Avoid division by 0
    // M = sqrt( \nabla J(row, 0)^2 + \nabla J(row, 1)^2 ) (for 2D)

    T M = domain.rowwise().norm().maxCoeff();

    //gsInfo << "Gradient: " << domain.block(0,0,m_tmpPatches[0].coefs().rows(),2) << "\n";

    if(M == 0.)
        M = 1.;
    /// Older, bigger approach
    //sg_vector /= M;
    //result = sg_vector.asVector();
    /// Newer, smaller approach
    red_vector /= M;
    result = red_vector.asVector();

    //gsInfo<< "Exact   Gradient: "<< result.transpose() <<"\n";
    //gsOptProblem<T>::gradObj_into(u,result);
    //gsInfo<< "Approx. Gradient: "<< result.transpose() <<"\n";


    // End: Make a vector out of the domain
}

template <typename T>
void gsOptElectricMotor<T>::evalCon_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
{
    //gsInfo << "I am here \n";
    // get the coefficients
    jacobianDetCoefs(u, result);

    // set the corners either to 1.0 or to -1.0 depending on the orientation of the patches
    for(std::vector<size_t>::const_iterator it = m_desDomain.begin(); it != m_desDomain.end(); it++) {
        // Disactivate the corners, since they have constant coefficient
        if(m_patchOrientation[*it])
            for (index_t i = 0; i != m_corner[*it].size(); ++i)
                result[m_corner[*it][i]] = T(1.0);
        else
            for (index_t i = 0; i != m_corner[*it].size(); ++i)
                result[m_corner[*it][i]] = T(-1.0);
    }

    for(std::vector<size_t>::const_iterator it = m_desNeighbours.begin(); it != m_desNeighbours.end(); it++) {
        // Disactivate the corners, since they have constant coefficient
        if(m_patchOrientation[*it])
            for (index_t i = 0; i != m_corner[*it].size(); ++i)
                result[m_corner[*it][i]] = T(1.0);
        else
            for (index_t i = 0; i != m_corner[*it].size(); ++i)
                result[m_corner[*it][i]] = T(-1.0);
    }

    //gsAsMatrix<> tmpmat(result.data(), result.rows()/4, 4);
    //gsInfo << "Result: " << tmpmat << "\n";
}

template <typename T>
void gsOptElectricMotor<T>::jacobianDetCoefs(const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
{
    //gsInfo << "This is patch 0: " << m_tmpPatches[0].coefs() << "\n";
    updateTempPatches(u);
    //gsInfo << "This is patch 0: " << m_tmpPatches[0].coefs() << "\n";

    index_t r = 0;
    /// Compute the jacobian determinant coefficients on the design domain ...
    for(std::vector<size_t>::const_iterator it = m_desDomain.begin(); it != m_desDomain.end(); it++)
    {
        gsVector<T> b(m_colPoints[*it].cols());
        gsMatrix<T> tmpJac;

        for (index_t i = 0; i != m_colPoints[*it].cols(); ++i)
        {
            m_basis.basis(*it).jacobianFunc_into(m_colPoints[*it].col(i), m_tmpPatches[*it].coefs(), tmpJac);
            b[i] = tmpJac.determinant();
        }

        // Solve for the constraints
        // note: result is expected to be allocated at least
        // m_numConstraints numbers

        gsSparseMatrix<T> colMat;
        m_jacBasis.basis(*it).collocationMatrix(m_colPoints[*it], colMat);
        m_jacSolver.compute(colMat);

        result.segment(r, m_colPoints[*it].cols()) = m_jacSolver.solve(b);
        r += m_colPoints[*it].cols();
    }

    /// ... also for the neighbours --> for tests!!!
    for(std::vector<size_t>::const_iterator it = m_desNeighbours.begin(); it != m_desNeighbours.end(); it++)
    {
        gsVector<T> b(m_colPoints[*it].cols());
        gsMatrix<T> tmpJac;

        for (index_t i = 0; i != m_colPoints[*it].cols(); ++i)
        {
            m_basis.basis(*it).jacobianFunc_into(m_colPoints[*it].col(i), m_tmpPatches[*it].coefs(), tmpJac);
            b[i] = tmpJac.determinant();
        }

        // Solve for the constraints
        // note: result is expected to be allocated at least
        // m_numConstraints numbers

        gsSparseMatrix<T> colMat;
        m_jacBasis.basis(*it).collocationMatrix(m_colPoints[*it], colMat);
        m_jacSolver.compute(colMat);

        result.segment(r, m_colPoints[*it].cols()) = m_jacSolver.solve(b);
        r += m_colPoints[*it].cols();
    }
}

template <typename T>
void gsOptElectricMotor<T>::jacobCon_into(const gsAsConstVector<T> & u, gsAsVector<T> &result) const // only one design domain TODO: check code in case that all magnets are considered
{
    updateTempPatches(u);

    //TODO: Sparse implementation
#ifdef JAC_STRUCTURE
/*
    const index_t dim = m_tmpPatches.dim();
    gsMatrix<T> JacAdj(dim, dim);
    gsMatrix<T> gradB;

    const index_t sz = m_mapper.freeSize();

    // Jacobian is 4536 x 4536 matrix
    // Note: result.rows() == m_conJacRows.size() !!

    for(std::vector<size_t>::const_iterator it = m_desDomain.begin(); it != m_desDomain.end(); it++) {

        index_t shift = 0;

        // Space for right-and side and solution vector
        gsMatrix<T> currDV(m_colPoints[*it].cols(), dim), b(m_colPoints[*it].cols(), dim);
        gsInfo << "DV rows: " << currDV.rows() << "\n";

        for (index_t k = 0; k != m_basis.basis(*it).size(); ++k) // for all control points
        {
            // get index of the design control point
            const index_t kk = m_mapper.index(k, *it);

            if (m_mapper.is_free_index(kk)) // interior ?
            {
                for (index_t i = 0; i != m_colPoints[*it].cols(); ++i) // for all constraints on detJ
                {
                    // Skip corner constraints
                    if ((m_corner[*it].array() == i).any()) {
                        b.row(i).setZero();
                        continue; // for
                    }

                    if (m_basis.basis(*it).isActive(k, m_colPoints[*it].col(i))) {
                        // Get Jacobian
                        m_basis.basis(*it).jacobianFunc_into(m_colPoints[*it].col(i), m_tmpPatches[*it].coefs(), JacAdj);

                        // Adjugate (cofactor transpose)
                        JacAdj.adjugateInPlace();

                        // Gradient of B_k
                        m_basis.basis(*it).derivSingle_into(k, m_colPoints[*it].col(i), gradB);

                        for (index_t s = 0; s != dim; ++s) // for all components
                            b(i, s) = (JacAdj.col(s) * gradB.transpose()).trace();
                    }
                }

                // Solve for the kk-th sensitivities
                gsSparseMatrix<T> colMat;
                m_jacBasis.basis(*it).collocationMatrix(m_colPoints[*it], colMat);
                m_jacSolver.analyzePattern(colMat);
                m_jacSolver.factorize(colMat);
                currDV = m_jacSolver.solve(b);

                // Fill in sparse result
                for (index_t s = 0; s != dim; ++s) // for all components
                {
                    const index_t rr = s * sz + kk;

                    for (index_t m = m_marker[rr]; m != m_marker[rr + 1]; ++m)// for all entries
                    {
                        //gsInfo << "m: " << m << "\n";
                        //gsInfo << "dv: " << currDV.row(m+shift) << "\n";
                        result[m] = currDV(m_conJacRows[m], s);
                    }
                }
            }
            //gsInfo << "shift: " << shift << "\n";
        }
    }

#else
    */
    ///Implementation of the derivatives
    gsInfo << "evaluating the derivatives of the constraint! " << "\n";
    const index_t dim = m_tmpPatches.dim();
    index_t sz = m_mapper.freeSize();

    gsMatrix<T> JacAdj(dim, dim);
    gsMatrix<T> gradB;
    gsMatrix<index_t> activeB;

    index_t r = 0;

    gsAsMatrix<T> JC(result.data(), m_numConstraints, dim*sz);
    gsMatrix<T> b;

    for(std::vector<size_t>::const_iterator it = m_desDomain.begin(); it != m_desDomain.end(); it++) {
        // Initialize right-hand side to zero
        // Assumption: Every design domain has the same number of design variables
        b.setZero(m_colPoints[*it].cols(), dim * sz);

        // Evaluate right hand side b
        for (index_t i = 0; i != m_colPoints[*it].cols(); ++i) {
            // m_geoBasis->compute(active, deriv, jacob, m_tmpCoefs);
            // gsBasis<T>::linearComb(active, evals, m_tmpCoefs, result);
            // gsBasis<T>::jacobianFromGradients(active, grads, m_tmpCoefs, result);

            // Get Jacobian determinant and inverse
            m_basis.basis(*it).jacobianFunc_into(m_colPoints[*it].col(i), m_tmpPatches[*it].coefs(), JacAdj);
            //Adjugate (cofactor transpose)
            JacAdj.adjugateInPlace();

            // Active basis function at point i
            //const int numActive =
            m_basis.basis(*it).active_into(m_colPoints[*it].col(i), activeB);
            const int numActive = activeB.rows();

            // Basis Gradients matrix
            m_basis.basis(*it).deriv_into(m_colPoints[*it].col(i), gradB);

            for (index_t k = 0; k != numActive; ++k) // For all active         //numactive instead of m_basis.basis(*it).size()
            {
                if (m_mapper.is_free(activeB(k, 0), *it)) { //activeB(k, 0) instead of k
                    // get index of the design control point index
                    index_t kk = m_mapper.index(activeB(k, 0), *it); //activeB(k, 0) instead of k

                    for (index_t s = 0; s != dim; ++s) // for all components
                    {
                        //gsInfo << "Index: " << kk % 8 << "\n";
                        b(i, s * sz + kk) =
                                // JacDet*( JacInv.col(s) * gradB.template block<d,1>(k*d,0).transpose() ).trace();
                                (JacAdj.col(s) * gradB.block(k * dim, 0, dim, 1).transpose()).trace();
                    }
                }

            }
        }

        // Solve for the Constraint Jacobian (sensitivities)
        //gsAsMatrix<T> JC(result.data(), m_numConstraints, b.cols());
        // Note: result may contain extra data
        // ( ie. m_numConstraints>b.cols() ), so we fill the
        // principal sub-matrix
        gsSparseMatrix<T> colMat;
        m_jacBasis.basis(*it).collocationMatrix(m_colPoints[*it], colMat);
        m_jacSolver.compute(colMat);

        /// new start
        gsMatrix<T> jacSens = m_jacSolver.solve(b);

        // Disactivate the corners by setting their sensitivities to 0
        for (index_t i = 0; i != m_corner[*it].size(); ++i)
            jacSens.row(m_corner[*it][i]).setZero();

        //jacSens.transposeInPlace();

        //for (index_t row = 0; row < b.rows(); row++)
        //{
        //    for (index_t s = 0; s != dim; s++)
        //    {
        JC.block(r, 0, b.rows(), b.cols()) = jacSens;
        //    }

        r += b.rows();
        //}

        //JC.block(r, 0, b.rows(), b.cols()) = m_jacSolver.solve(b);

        //JC.block(0, 0, b.rows(), b.cols()) = lusolver.solve(b);
        // Deactivate the corners by setting their sensitivities to 0
        //for (index_t i = 0; i != m_corner[*it].size(); ++i)
            //JC(m_corner[*it][i], m_numDesignVars - 1) = 0.0;
        //    JC.block(r, 0, b.rows(), b.cols()).row(m_corner[*it][i]).setZero();

    }

    for(std::vector<size_t>::const_iterator it = m_desNeighbours.begin(); it != m_desNeighbours.end(); it++)
    {
        // Initialize right-hand side to zero
        //gsMatrix<T> b = gsMatrix<T>::Zero(m_colPoints[*it].cols(), dim * sz / T(m_desNeighbours.size()));
        b.setZero(m_colPoints[*it].cols(), dim * sz);

        // Evaluate right hand side b
        for (index_t i = 0; i != m_colPoints[*it].cols(); ++i)
        {
            // m_geoBasis->compute(active, deriv, jacob, m_tmpCoefs);
            // gsBasis<T>::linearComb(active, evals, m_tmpCoefs, result);
            // gsBasis<T>::jacobianFromGradients(active, grads, m_tmpCoefs, result);

            // Get Jacobian determinant and inverse
            m_basis.basis(*it).jacobianFunc_into(m_colPoints[*it].col(i), m_tmpPatches[*it].coefs(), JacAdj);
            //Adjugate (cofactor transpose)
            JacAdj.adjugateInPlace();

            // Active basis function at point i
            //const int numActive =
            m_basis.basis(*it).active_into(m_colPoints[*it].col(i), activeB);
            const int numActive = activeB.rows();

            // Basis Gradients matrix
            m_basis.basis(*it).deriv_into(m_colPoints[*it].col(i), gradB);

            for (index_t k = 0; k != numActive; ++k) // For all active
            {
                if(m_mapper.is_free(activeB(k, 0), *it))
                {
                    // get index of the design control point index
                    index_t kk = m_mapper.index(activeB(k, 0), *it); //activeB(k, 0);

                    for (index_t s = 0; s != dim; ++s) // for all components
                    {
                        //gsInfo << "Index: " << s * sz / m_desNeighbours.size() + kk % m_desNeighbours.size() << "\n";
                        b(i, s * sz + kk) =
                                // JacDet*( JacInv.col(s) * gradB.template block<d,1>(k*d,0).transpose() ).trace();
                                (JacAdj.col(s) * gradB.block(k * dim, 0, dim, 1).transpose()).trace();
                    }
                }

            }
        }

        // Solve for the Constraint Jacobian (sensitivities)
        //gsAsMatrix<T> JC(result.data(), m_numConstraints, b.cols());
        // Note: result may contain extra data
        // ( ie. m_numConstraints>b.cols() ), so we fill the
        // principal sub-matrix
        gsSparseMatrix<T> colMat;
        m_jacBasis.basis(*it).collocationMatrix(m_colPoints[*it], colMat);
        m_jacSolver.compute(colMat);

        /// new start
        gsMatrix<T> jacSens = m_jacSolver.solve(b);

        // Disactivate the corners by setting their sensitivities to 0
        for (index_t i = 0; i != m_corner[*it].size(); ++i)
            jacSens.row(m_corner[*it][i]).setZero();

        JC.block(r, 0, b.rows(), b.cols()) = jacSens;

            //r += sz / T(m_desNeighbours.size());
        r += b.rows();

    }

    GISMO_ASSERT(r == m_numConstraints, "Something went wrong with the dense implementation of the jacobian constraints");

    //check the order of the jacobian of the contraints
    //std::vector<std::pair<index_t, index_t> > preIm;
    //m_mapper.preImage(3, preIm);
    //for(int i = 0; i < preIm.size(); i++)
    //    gsInfo << "patch: " << preIm[i].first << "\n";


    // finite differences for test purposes
    // -> works!! derivative checker matches in all components
    /*
    real_t htol = 6e-7;
    gsMatrix<T> pert = u;
    gsVector<T> con(m_numConstraints);
    gsAsConstVector<T> pertVec(pert.data(), u.size());
    gsAsVector<T> resultCopy(con.data(), m_numConstraints);
    //gsMatrix<T> JM, jm;

    for(index_t i = 0; i < m_numDesignVars; i++)
    {
        pert(i) = u(i) + htol;
        evalCon_into(pertVec, resultCopy);
        const gsMatrix<T> JM = resultCopy;

        pert(i) = u(i) - htol;
        evalCon_into(pertVec, resultCopy);
        const gsMatrix<T> jm = resultCopy;

        for(index_t j = 0; j < m_numConstraints; j++)
        {
            result(i*m_numConstraints+j) = (JM(j) - jm(j)) / (2*htol);
        }

        pert(i) = u(i);
    }
    */


#endif
}

template <typename T>
void gsOptElectricMotor<T>::solve()
{
#ifdef GISMO_WITH_IPOPT
    gsOptProblem<T>::solve();
#else
    //*********
    // Optimization loop
    // calling evalObj_into, gradObj_into

    // Create an empty gsAsVector object
    gsVector<T> curDesignCopy(m_numDesignVars), tmpvector(m_numDesignVars), constrVector(m_numConstraints);
    gsAsVector<T> shapeGradAsVec(tmpvector.data(), m_numDesignVars);
    gsAsConstVector<T> copydomAsConstVec(curDesignCopy.data(), m_numDesignVars);
    gsAsVector<T> constrAsVec(constrVector.data(), m_numConstraints);
    /// !!Does not work! Bug?
    //gsAsConstVector<T> curDesignAsConstVector(m_curDesign.data(), m_numDesignVars);

    // Error tolerance
    T currVal, newVal; // --> member m_tolerance;
    unsigned acc_iter = 0;

    for(unsigned k = 0; k < m_maxIterations; k++)
    {

        gsInfo << "This is iteration no.: " << k << "\n";
        // Evaluate the objective
        //gsInfo << "New curr design: " << m_curDesign << "\n";
        gsAsConstVector<T> curDesignAsConstVector(m_curDesign.data(), m_numDesignVars);
        currVal = evalObj(curDesignAsConstVector);

        // Compute the shape gradient
        gradObj_into(curDesignAsConstVector, shapeGradAsVec);

        newVal = std::numeric_limits<T>::max();

        // Line search point count
        unsigned t = 0;

        // Initial stepsize
        const T h = 1.; // 0.05

        // Do line search
        while(currVal <= newVal &&t<50)
        {
            curDesignCopy.noalias() = m_curDesign - ( h/(1<<t) ) * shapeGradAsVec; // /M -> division is already done in evalGrad_into
            evalCon_into(copydomAsConstVec, constrAsVec);

            if((constrAsVec.array() <= m_conUpperBounds.array()).all() && (constrAsVec.array() >= m_conLowerBounds.array()).all())
            {
                newVal = evalObj(copydomAsConstVec);
            }
            else
            {
                t++;
                continue;
            }

            gsInfo << "New functional value: " << newVal << " and old one: " <<  currVal <<  "\n";

            if( newVal < currVal)
            {
                // Update the current
                m_curDesign = copydomAsConstVec;
                // curDesignAsConstVector = copydomAsConstVec; // equivalent ??

                finalObjective = newVal;

                break;
            }
            t++;
        }


        if (50==t) gsWarn << "Line search failed after " << t << "steps.\n";

        if( (currVal - finalObjective) <= m_tolerance)
        {
            acc_iter++;
            if(acc_iter == 4)
                break;
        }

        numIterations++;
    }
#endif
}

template <typename T>
gsMultiPatch<T> gsOptElectricMotor<T>::solveStateForMultipatch(gsMagnetostaticPde <T> magPDE,
                                                               gsOptionList opt,
                                                               gsPiecewiseFunction <T> f) const
{
    gsMagnetostaticAssembler<T> magAss(magPDE, m_basis, opt, f);

    gsMatrix<T> state = solveState(magAss);
    gsMultiPatch<> result;

    magAss.constructSolution(state, result);
    return result;
}

template <typename T>
gsMatrix<T> gsOptElectricMotor<T>::solveState(gsMagnetostaticAssembler <T> & magAss) const
{
    //for(size_t i = 0; i < 3; i++)
    //    multiBasis.uniformRefine();

    gsStopwatch time;
    magAss.assembleFull();
    gsSparseMatrix<> K = magAss.matrix();
    //gsInfo << "Eigenvalues of the full matrix: \n";
    //gsInfo << K.toDense().eigenvalues().real().minCoeff() << "\n";
    //gsInfo<<"Start calculating eigenvalues, #dofs: "<<K.cols()<<"\n";
    //Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<real_t> > eigenSolver(K.cols());
    //eigenSolver.compute(K);
    //gsInfo<<"Condition number of K:"<<eigenSolver.eigenvalues().maxCoeff()/eigenSolver.eigenvalues().minCoeff()<<"\n\n";

    gsSparseSolver<>::LU solver;
    solver.compute(K);
    // !!! magAss.rhs() gives only the 0th piece of the function !!!!!
    gsMatrix<> sol = solver.solve( magAss.getRhsFull() );
    //gsInfo << "Time needed to solve the equation on this domain with a direct solver: " << time << "\n";

    time.restart();
    gsIETIAssembler<real_t> EM_IETIAss(magAss);

    gsOptionList defaultOpt = gsIETIAssembler<real_t>::defaultOptions();
    //defaultOpt.setString("Scaling", "mult");

    EM_IETIAss.setOptions(defaultOpt);

    EM_IETIAss.init();

    EM_IETIAss.assemble();

    //Setup the linearoperator for the iterative solver
    gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(EM_IETIAss);
    solv->init();

    //Setup the preconditioner
    gsScaledDirichletPrecond<real_t>::Ptr precDir = gsScaledDirichletPrecond<real_t>::make(EM_IETIAss);
    //Setup the CG
    gsConjugateGradient<> PCG(solv,precDir);
    PCG.setMaxIterations(100);
    PCG.setCalcEigenvalues(true);

    //Intial guess
    gsMatrix<> solVector(EM_IETIAss.systemSize(),EM_IETIAss.numberRhs());
    solVector.setZero();

    //Solve the IETI system, we obtain the lagrange multipliers
    PCG.solve(solv->getRhs(),solVector);
    time.stop();
    //gsInfo << "Time needed to solve the equation on this domain with the IETI solver: " << time << "\n";

    //Information about the eigenvalues of the preconditioned system
    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);
    gsInfo<<"eigenvalues: "<<"\n"<<eigs.minCoeff()<<" -- "<<eigs.maxCoeff()<<"\n";
    //gsInfo<<"Number of iterations: "<<PCG.iterations()<<"\n";

    //Backsubstituting the lagrange multipliers in order to obtain the solution from Ku=f
    gsMatrix<> solution;
    gsMultiPatch<> result;

    solv->calculateSolution(solVector,solution);

    return solution;

}

template <typename T>
gsMultiPatch<T> gsOptElectricMotor<T>::solveAdjoint(gsMagnetostaticAdjointPde<T> magPDE, gsMultiPatch<T> mp, gsOptionList opt) const
{
    //for(size_t i = 0; i < 3; i++)
    //    multiBasis.uniformRefine();

    std::vector<std::pair<unsigned, unsigned> > ints;
    createInterfaceTopology(ints);
    gsMagnetostaticAdjointAssembler<T> magAss(magPDE, m_basis, opt, ints);

    gsStopwatch time;

    magAss.assemble();
    gsSparseMatrix<> K = magAss.matrix();
    gsSparseSolver<>::LU solver;
    solver.compute(K);
    // !!! magAss.rhs() gives only the 0th piece of the function !!!!!
    gsMatrix<> sol = solver.solve( magAss.getRhsFull() );
    //gsInfo << "Time needed to solve the equation on this domain with a direct solver: " << time << "\n";

    time.restart();
    gsIETIAssembler<real_t> EM_IETIAss(magAss);

    gsOptionList defaultOpt = gsIETIAssembler<real_t>::defaultOptions();

    EM_IETIAss.setOptions(defaultOpt);
    EM_IETIAss.init();

    EM_IETIAss.assemble();

    //Setup the linearoperator for the iterative solver
    gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(EM_IETIAss);
    solv->init();

    //Setup the preconditioner
    gsScaledDirichletPrecond<real_t>::Ptr precDir = gsScaledDirichletPrecond<real_t>::make(EM_IETIAss);
    //Setup the CG
    gsConjugateGradient<> PCG(solv,precDir);
    PCG.setMaxIterations(100);
    PCG.setCalcEigenvalues(true);

    //Intial guess
    gsMatrix<> solVector(EM_IETIAss.systemSize(),EM_IETIAss.numberRhs());
    solVector.setZero();

    //Solve the IETI system, we obtain the lagrange multipliers
    PCG.solve(solv->getRhs(),solVector);
    time.stop();
    //gsInfo << "Time needed to solve the equation on this domain with the IETI solver: " << time << "\n";

    //Information about the eigenvalues of the preconditioned system
    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);
    gsInfo<<"eigenvalues: "<<"\n"<<eigs.minCoeff()<<" -- "<<eigs.maxCoeff()<<"\n";
    //gsInfo<<"Number of iterations: "<<PCG.iterations()<<"\n";

    //Backsubstituting the lagrange multipliers in order to obtain the solution from Ku=f
    gsMatrix<> solution;
    gsMultiPatch<> result;

    solv->calculateSolution(solVector,solution);
    //gsInfo << "The difference between the solution techniques is: " << (solution - sol).norm() << "\n";

    magAss.constructSolution(solution, result);

    return result;
}

template <typename T>
gsMultiPatch<T> gsOptElectricMotor<T>::calcShapeDerivative(gsMagnetostaticShapeDerivPde<T> shapeD, gsOptionList opt, gsBoundaryConditions<> bcInfo) const
{
    //gsMultiBasis<> multiBasis(shapeD.domain());

    //for(size_t i = 0; i < 1; i++)
    //    multiBasis.uniformRefine();

    gsMagnetostaticShapeDerivAssembler_decoupled<T> magAss(shapeD, m_basis, opt);

    gsStopwatch time;
    magAss.assemble();
    gsSparseMatrix<> K = magAss.matrix();
    gsSparseSolver<>::LU solver;
    solver.compute(K);
    // !!! magAss.rhs() gives only the 0th piece of the function !!!!!
    gsMatrix<> sol = solver.solve(magAss.getRhsFull());
    //gsInfo << "Time needed to solve the equation on this domain with a direct solver: " << time << "\n";

    time.restart();
    gsIETIAssembler<real_t> EM_IETIAss(magAss);

    gsOptionList defaultOpt = gsIETIAssembler<real_t>::defaultOptions();
    defaultOpt.setInt("nRhs", 2);
    defaultOpt.setString("Scaling", "coeff");
    defaultOpt.setString("Strategy", "C");

    EM_IETIAss.setOptions(defaultOpt);
    EM_IETIAss.init();

    EM_IETIAss.assemble();

    //Setup the linearoperator for the iterative solver
    gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(EM_IETIAss);
    solv->init();

    //Setup the preconditioner
    gsScaledDirichletPrecond<real_t>::Ptr precDir = gsScaledDirichletPrecond<real_t>::make(EM_IETIAss);
    //Setup the CG
    gsConjugateGradient<> PCG(solv, precDir);
    PCG.setMaxIterations(100);
    PCG.setCalcEigenvalues(true);

    //Intial guess
    gsMatrix<> solVector(EM_IETIAss.systemSize(), EM_IETIAss.numberRhs());
    solVector.setZero();

    //Solve the IETI system, we obtain the lagrange multipliers
    for (index_t i = 0; i < EM_IETIAss.numberRhs(); i++) {
        gsMatrix<> tmpSolution;
        PCG.solve(solv->getRhs(i), tmpSolution);
        solVector.col(i) = tmpSolution;
    }
    time.stop();
    //gsInfo << "Time needed to solve the equation on this domain with the IETI solver: " << time << "\n";

    //Information about the eigenvalues of the preconditioned system
    gsMatrix<> eigs;
    PCG.getEigenvalues(eigs);
    gsInfo << "eigenvalues: " << "\n" << eigs.minCoeff() << " -- " << eigs.maxCoeff() << "\n";
    //gsInfo << "Number of iterations: " << PCG.iterations() << "\n";

    //Backsubstituting the lagrange multipliers in order to obtain the solution from Ku=f
    gsMatrix<> solution;
    gsMultiPatch<> result;

    solv->calculateSolution(solVector, solution);
    //gsInfo << "The difference between the solution techniques is: " << (solution - sol).norm() << "\n";

    magAss.constructSolution(solution, result, 0);

    return result;
}

template <typename T>
bool gsOptElectricMotor<T>::computeInitialjacobianDetCoefs(const size_t designPatch)
{
    // get the coefficients
    gsVector<T> b(m_colPoints[designPatch].cols());
    gsMatrix<T> tmpJac;

    for (index_t i = 0; i != m_colPoints[designPatch].cols(); ++i) {
        m_basis.basis(designPatch).jacobianFunc_into(m_colPoints[designPatch].col(i), m_tmpPatches[designPatch].coefs(), tmpJac);
        b[i] = tmpJac.determinant();
    }

    if( (b.array() > 0).all() )
        return true;
    else
        return false;
}

template <typename T>
void gsOptElectricMotor<T>::updateTempPatches(const gsAsConstVector<T> & u) const
{
    typedef typename gsGeometry<T>::uPtr geoPtr;
    int dimLen = u.size() / m_tmpPatches.dim(); // control coefficients per dimension

    for (size_t freeDof = 0; freeDof < m_preImages.size(); freeDof++) // m_preImages.size() = m_numDesVars
    {
        std::vector<std::pair<index_t, index_t> > preImage = m_preImages[freeDof];
        for (size_t nPre = 0; nPre < preImage.size(); nPre++) {
            index_t patch = preImage[nPre].first;
            index_t dof = preImage[nPre].second;

            gsMatrix<T> &cf = m_tmpPatches.patch(patch).coefs();

            for (short_t d = 0; d < m_tmpPatches.dim(); d++)
                cf(dof, d) = u(freeDof + d * dimLen);
        }
    }
    //gsInfo << "This is patch 60: " << m_tmpPatches[60].coefs() << "\n";
    patchFromBoundary(u);
    //gsInfo << "This is also patch 60: " << m_tmpPatches[60].coefs() << "\n";

    /*
    typedef typename gsGeometry<T>::uPtr geoPtr;
    int dimLen = u.size() / m_tmpPatches.dim(); // control coefficients per dimension

    gsInfo << "Length: " << u.size() << "\n";
    index_t c = 0;
    for(index_t np = 0; np < m_tmpPatches.nPatches(); np++)
    {
            gsMatrix<T> &cf = m_tmpPatches.patch(np).coefs();
            for (int d = 0; d < m_tmpPatches.dim(); d++)
                cf.col(d) = u.segment(c + d * dimLen, cf.rows());

            c += cf.rows();
    }
*/
}


template<typename T>
void gsOptElectricMotor<T>::patchFromBoundary(const gsAsConstVector<T> & u) const
{
    // Reorder the inner control points of the design domain
    for (std::vector<size_t>::const_iterator ddom = m_desDomain.begin(); ddom != m_desDomain.end(); ++ddom)
    {
        gsMultiPatch<T> desDom(m_tmpPatches[*ddom]);
        desDom.computeTopology();
        std::vector<patchSide> sides = desDom.boundaries();
        std::vector<gsGeometry<T>* > bnd;

        for(size_t i = sides.size()-1; i < -1UL; --i) // --> this direction preserves the right order of the coefficients
        {
            gsGeometry<T>* geo = (m_tmpPatches[*ddom].boundary(sides[i])).release();
            bnd.push_back(geo);
        }
        gsMultiPatch<T> boundarymp(bnd);
        boundarymp.closeGaps(1e-12);

        gsSpringPatch<T> newPatch(boundarymp);
        newPatch.compute();

        m_tmpPatches.patch(*ddom).coefs() = newPatch.result().coefs();
    }

    //Reorder the inner control points of the neighbouring elements
    for(std::vector<size_t>::const_iterator ndom = m_desNeighbours.begin(); ndom != m_desNeighbours.end(); ++ndom)
    {
        gsMultiPatch<T> neighbour(m_tmpPatches[*ndom]);
        neighbour.computeTopology();
        std::vector<patchSide> sides = neighbour.boundaries();
        std::vector<gsGeometry<T>* > bnd;

        for(size_t i = sides.size()-1; i < -1UL; --i) // --> this direction preserves the right order of the coefficients
        {
            gsGeometry<T>* geo = (m_tmpPatches[*ndom].boundary(sides[i])).release();
            bnd.push_back(geo);
        }
        gsMultiPatch<T> boundarymp(bnd);
        boundarymp.closeGaps(1e-12);

        gsSpringPatch<T> newPatch(boundarymp);
        newPatch.compute();

        m_tmpPatches.patch(*ndom).coefs() = newPatch.result().coefs();
    }

    m_tmpPatches.closeGaps();
    m_tmpPatches.computeTopology();

    //gsWriteParaview(m_tmpPatches, "Current");
}


template <typename T>
gsMatrix<T> gsOptElectricMotor<T>::multipatchToVector(const gsMultiPatch<T> &mp) const
{
    gsMatrix<> copyDomain;

    for(size_t np = 0; np < mp.nPatches(); np++)
    {
        gsMatrix<T> & cf = (mp.patch(np)).coefs();
        copyDomain.conservativeResize(copyDomain.rows() + cf.rows(), m_tmpPatches.dim());

        // Append the matrix
        for(index_t j = 0; j < cf.cols(); j++)
            copyDomain.col(j).tail(cf.rows()) << cf.col(j);


    }

    return ( copyDomain.asVector() );


    // m_numDesignVars not defined yet
/*
    gsMatrix<T> result(m_numDesignVars, 1);
    gsAsMatrix<T> rr(result.data(),m_numDesignVars/2,2);
    index_t c = 0;
    for(size_t np = 0; np < mp.nPatches(); np++)
    {
        gsMatrix<T> & cf = mp.patch(np).coefs();
        rr.middleRows(c,cf.rows() ) = cf;
        c+= cf.rows();
    }
    return result;
    */
}

    /*
template <typename T>
void gsOptElectricMotor<T>::computeJacStructure()
{
    //to try: typedef gsSortedVector<std::pair<index_t,index_t> > pairList;
    typedef std::set<std::pair<index_t,index_t> > pairList;
    gsMatrix<index_t> actRow, actCol;
    pairList pairs; // pair: (DesVar,Constr)
    int sz = m_mapper.freeSize();

    //gsInfo << "full mapper: " << m_mapper.size() << "\n";
    //gsInfo << "free mapper: " << sz << "\n";

    m_numConJacNonZero = 0;
    for(std::vector<size_t>::iterator it = m_desDomain.begin(); it != m_desDomain.end(); it++)
    {
        typename gsBasis<T>::domainIter domIt = m_basis.basis(*it).makeDomainIterator();

        for (; domIt->good(); domIt->next())
        {
            m_jacBasis.basis(*it).active_into(domIt->center, actRow); //Returns the indices of active (non-zero) basis functions at the center, as a list of indices, in actRow.
            m_basis.basis(*it).active_into(domIt->center, actCol); // In this case actRow and actCol are column vectors

            for (index_t j = 0; j != actCol.rows(); ++j) {
                const index_t jj = m_mapper.index(actCol(j), *it); // Returns the global dof index

                if (m_mapper.is_free_index(jj))
                    for (index_t i = 0; i != actRow.rows(); ++i)
                        for (index_t s = 0; s != m_tmpPatches.dim(); ++s) // for all components
                            pairs.insert(std::make_pair(s * sz + jj, actRow(i)));
            }
        }

        // sz now bacomes the size of the constraint jacobian related to detJ
        sz = pairs.size();
        m_conJacRows.reserve(sz);
        m_conJacCols.reserve(sz);

        // the marker encodes the columns: we don't really need
        // m_conJacCols anymore, due to sorting wrt design vars!
        m_marker.setZero(m_tmpPatches.dim() * sz + 1);

        for (pairList::const_iterator pit = pairs.begin(); pit != pairs.end(); ++pit) {
            if ( (m_corner[*it].array() == pit->second).any() )
                continue;

            m_conJacRows.push_back(pit->second);//constraint
            m_conJacCols.push_back(pit->first);//design var

            ++m_marker[pit->first + 1];
        }

        // Complete the marking
        for (index_t k = 1; k <= m_tmpPatches.dim() * sz; ++k)
            m_marker[k] += m_marker[k - 1];

        //m_numConJacNonZero = m_conJacRows.size();

    }
    m_numConJacNonZero += m_conJacRows.size();

}
    */

template <typename T>
void gsOptElectricMotor<T>::computeJacStructure()
{
    index_t desVars = 0;

    for(std::vector<size_t>::iterator it = m_desDomain.begin(); it != m_desDomain.end(); it++)
    {
        desVars += 2 * m_tmpPatches[*it].coefs().rows();
    }

    //m_conJacRows.resize(m_numConstraints*desVars);
    //m_conJacCols.resize(m_numConstraints*desVars);
    m_conJacRows.resize(m_numConstraints*m_numDesignVars);
    m_conJacCols.resize(m_numConstraints*m_numDesignVars);
    index_t c = 0;
    for ( index_t j = 0; j < m_numDesignVars; ++j)
        for ( index_t i = 0; i < m_numConstraints; ++i)
        {
            m_conJacRows[c  ] = i;
            m_conJacCols[c++] = j;
        }

    m_numConJacNonZero = m_conJacRows.size();

}


template <typename T>
void gsOptElectricMotor<T>::setTolerance(const T eps)
{
    m_tolerance = eps;
}

void createInterfaceTopology(std::vector<std::pair<unsigned, unsigned > > & interfaces)
{
    interfaces.clear();

    std::pair<unsigned, unsigned> interface_pair;

    for(int quad = 0; quad < 4; quad++) {
        interface_pair = std::make_pair(44+93*quad, 27+93*quad);
        interfaces.push_back(interface_pair);
        interface_pair = std::make_pair(43+93*quad, 28+93*quad);
        interfaces.push_back(interface_pair);
        interface_pair = std::make_pair(42+93*quad, 29+93*quad);
        interfaces.push_back(interface_pair);
        interface_pair = std::make_pair(41+93*quad, 30+93*quad);
        interfaces.push_back(interface_pair);
        interface_pair = std::make_pair(40+93*quad, 31+93*quad);
        interfaces.push_back(interface_pair);
        interface_pair = std::make_pair(32+93*quad, 39+93*quad);
        interfaces.push_back(interface_pair);
        interface_pair = std::make_pair(38+93*quad, 33+93*quad);
        interfaces.push_back(interface_pair);
        interface_pair = std::make_pair(37+93*quad, 34+93*quad);
        interfaces.push_back(interface_pair);
        interface_pair = std::make_pair(36+93*quad, 35+93*quad);
        interfaces.push_back(interface_pair);
    }

/*
    std::pair<unsigned, unsigned> interface_pair(44, 27);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(43, 28);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(42, 29);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(41, 30);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(40, 31);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(32, 39);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(38, 33);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(37, 34);
    interfaces.push_back(interface_pair);
    interface_pair = std::make_pair(36, 35);
    interfaces.push_back(interface_pair);
    */
}

gsMultiPatch<> createInterfacePatch(const gsMultiPatch<> & mp)
{
    gsMultiPatch<> InterfaceMP(mp);
    InterfaceMP.clearTopology();

    for(int quad = 0; quad < 4; quad++) {
        InterfaceMP.addInterface(&InterfaceMP.patch(44+93*quad), 2, &InterfaceMP.patch(27+93*quad), 3);
        InterfaceMP.addInterface(&InterfaceMP.patch(43+93*quad), 4, &InterfaceMP.patch(28+93*quad), 3);
        InterfaceMP.addInterface(&InterfaceMP.patch(42+93*quad), 4, &InterfaceMP.patch(29+93*quad), 4);
        InterfaceMP.addInterface(&InterfaceMP.patch(41+93*quad), 4, &InterfaceMP.patch(30+93*quad), 4);
        InterfaceMP.addInterface(&InterfaceMP.patch(40+93*quad), 4, &InterfaceMP.patch(31+93*quad), 3);
        InterfaceMP.addInterface(&InterfaceMP.patch(32+93*quad), 3, &InterfaceMP.patch(39+93*quad), 3);
        InterfaceMP.addInterface(&InterfaceMP.patch(38+93*quad), 4, &InterfaceMP.patch(33+93*quad), 3);
        InterfaceMP.addInterface(&InterfaceMP.patch(37+93*quad), 4, &InterfaceMP.patch(34+93*quad), 3);
        InterfaceMP.addInterface(&InterfaceMP.patch(36+93*quad), 4, &InterfaceMP.patch(35+93*quad), 1);
    }

    return InterfaceMP;
    /*
    InterfaceMP.addInterface(&InterfaceMP.patch(44), 2, &InterfaceMP.patch(27), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(43), 4, &InterfaceMP.patch(28), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(42), 4, &InterfaceMP.patch(29), 4);
    InterfaceMP.addInterface(&InterfaceMP.patch(41), 4, &InterfaceMP.patch(30), 4);
    InterfaceMP.addInterface(&InterfaceMP.patch(40), 4, &InterfaceMP.patch(31), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(32), 3, &InterfaceMP.patch(39), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(38), 4, &InterfaceMP.patch(33), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(37), 4, &InterfaceMP.patch(34), 3);
    InterfaceMP.addInterface(&InterfaceMP.patch(36), 4, &InterfaceMP.patch(35), 1);
     */
}



} // end namespace gismo
