/** @file gsExprAssembler_test.cpp

    @brief Diffuse interface method

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris and Ct Wu
*/

#include <gismo.h>

#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>

using namespace gismo;
using namespace gismo::expr;


// testing solution evaluation
void SolutionEvaluator(gsBasis<> &domainBasis, const gsMatrix<> &coef, const gsMatrix<> &x, gsMatrix<> &val);

// TYPES: 1 completely int. / 0 interface / -1 completeely outsid
int CalElmType(gsFunctionExpr<>& disFunc, gsMatrix<> box);

// calculates active Dofs for the elmTypes
void CalActiveDofs(gsBasis<> &Basis, gsMultiBasis<> &Mesh, std::vector<int> &elmType, std::vector<bool> &dofStatus);

//adaptive ref.
bool AdaptiveRefinement(gsMultiBasis<>& Basis, gsFunctionExpr<>& disFunc, real_t eps);


int main(int argc, char *argv[]) 
{
    // Command line arguments
    bool plot = false;
    index_t nRefine = 1; // number of uniform refine steps to intial knot vector.
    index_t degree = 2;
    real_t eps = 0.03;
    bool c0Basis = false;
    gsCmdLine cmd("Diffuse interface method.");
    cmd.addSwitch("plot", "Create a ParaView visualization", plot);
    cmd.addSwitch("c0Basis", "Use C0 basis functions", c0Basis);
    cmd.addInt("r", "nRefine", "# of inserted knots in uniform refinement", nRefine);
    cmd.addInt("d", "degree", "Degree of Spline space", degree);
    cmd.addReal("e", "epsilon", "Epsilon parameter", eps);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    if (-1==eps) eps = ( 0!=nRefine ?  0.01 / nRefine : 0.03 );
    
    // Parameters
    real_t nu = 1, g = 9.81; // K = 1, alpha = 1, p = 1
    
    // Body force
    gsFunctionExpr<> f("(-4.905 + (1 - 4.905)*y) * cos(x)",
                       "(4.905 - 2 - 4.905*y + (0.5 - 2.4525)*y^2)*sin(x)", 2);
    // Velocity 
    gsFunctionExpr<> u_S ("(-4.905 +(1-4.905)*y)*cos(x)",
                          "(-1-4.905*y+(0.5 - 2.4525)*y^2)*sin(x)", 2);
    // Pressure
    gsConstantFunction<> p_S (1., 2); 
    // Piezometric head
    gsFunctionExpr<> phi_D ("exp(y)*sin(x) + 1/9.81", 2); 
    // Phase field
    gsFunctionExpr<> c("0.5 * (1+tanh((y-0.) / v) )", 2);
    c.set_v(eps);
    
    gsFunctionExpr<> dc("0", " (0.5/v) * ( 1-tanh((y-0.) / v)^2 )", 2);
    dc.set_v(eps);
    
    // Domain indicator
    gsFunctionExpr<> dind("if( y>0, 1, 0)", 2);
    
    // Tangent vector
    gsConstantFunction<> tau(1., 0., 2);

    // boundary data
    gsFunctionExpr<> t_east("-1-2*(-4.905 + (1-4.905)*y)*sin(x)",
                            " (-1-4.905*y+(0.5 - 2.4525)*y^2 + 1-4.905)*cos(x)", 2);
    gsFunctionExpr<> t_west("1+2*(-4.905 + (1-4.905)*y)*sin(x)",
                            "-(-1-4.905*y+(.5 - 2.4525)*y^2 + 1-4.905)*cos(x)", 2);
    gsFunctionExpr<> h_east("-exp(y)*cos(x)", 2);
    gsFunctionExpr<> h_west(" exp(y)*cos(x)", 2);
    
    // Knots first direction
    gsKnotVector<> Kx(0,EIGEN_PI, 2, degree+1);
    // Knots second direction
    gsKnotVector<> Ky(-1, 1, 1, degree+1);
    
    //0227 CT test for c0 basis functions
    if (c0Basis) {
        if (nRefine != 0) {
            Kx.uniformRefine(nRefine);
            Ky.uniformRefine(nRefine);
        }
        if (degree > 1) {
            Kx.increaseMultiplicity();
            Ky.increaseMultiplicity();
        }
    }
    //
    
    // Spline space of degree 2
    gsTensorBSplineBasis<2> spline2(give(Kx), give(Ky));
    
    // Domain
    gsTensorBSpline<2> spdomain( give(spline2), spline2.anchors().transpose() );
    
    gsMultiPatch<> domain(spdomain);
    domain.computeTopology();
    gsMultiBasis<> sp2(domain);
    
    // refine uniformly
    if (nRefine != 0 && !c0Basis)
        sp2.uniformRefine(nRefine);
    
    // Spline space of degree 3
    gsMultiBasis<> sp3 = sp2;
    sp3.degreeElevate();

    gsInfo << "Epsilon : " << eps << "\n";
    gsInfo << "Dofs of degree " << degree+1 << " Spline basis: " << sp3.totalSize() << "\n";
    gsInfo << "Dofs of degree " << degree << " Spline basis: " << sp2.totalSize() << "\n";
    
//    std::vector<unsigned> uDofs;
//    std::vector<unsigned> pDofs;
//    ActiveDofs(sp3.basis(0), uDofs);
//    ActiveDofs(sp2.basis(0), pDofs);
//    
//    gsInfo << "Active dofs of degree " << degree+1 << " Spline basis: " << uDofs.size() << "\n";
//    gsInfo << "Active dofs of degree " << degree << " Spline basis: " << pDofs.size() << "\n";
    
    
    // HBSpline basis
    gsHBSplineBasis<2> hb(sp2.basis(0));
    
    // Generate adaptive refined mesh for integration
    gsMultiBasis<> subEl(hb);
    gsFunctionExpr<> disFunc("y", 2);
    
    gsInfo << "Generating adaptive refined mesh for intergation:\n";
//    AdaptiveRefinement2(nRefine, subEl, disFunc, eps);
    AdaptiveRefinement(subEl, disFunc, eps);
    
    //----------------------------------------------------------------
    //Identify element type
    gsBasis<>::domainIter domIt = subEl.basis(0).makeDomainIterator();
    std::vector<int> elmType;
    
    gsMatrix<> element(sp3.dim(), 2);
    
    for (; domIt->good(); domIt->next())
    {
        element.col(0) = domIt->lowerCorner();
        element.col(1) = domIt->upperCorner();
        
        elmType.push_back(CalElmType(disFunc, element));
    }
#if __cplusplus > 199711L 
    elmType.shrink_to_fit();
#else 
    // shrink to fit idiom 
    std::vector<int>(elmType).swap(elmType); 
#endif 
    
    //Cal active dofs
    std::vector<bool> uDofs;
    std::vector<bool> pDofs;
    
    CalActiveDofs(sp3.basis(0), subEl, elmType, uDofs);
    CalActiveDofs(sp2.basis(0), subEl, elmType, pDofs);
    
    int cnt = 0;
    for (size_t i = 0; i < uDofs.size(); ++i)
    {
        if (uDofs[i])
            cnt++;
    }
    gsInfo << "Active dofs of degree " << degree+1 << " Spline basis: " << cnt << "\n";

    cnt = 0;
    for (size_t i = 0; i < pDofs.size(); ++i)
    {
        if (pDofs[i])
            cnt++;
    }
    gsInfo << "Active dofs of degree " << degree << " Spline basis: " << cnt << "\n";
    //----------------------------------------------------------------
    
    // Boundary conditions
    gsBoundaryConditions<> bc;
    bc.add(0, boundary::north, "Dirichlet", &u_S   , 0, -1, true);
    bc.add(0, boundary::east , "Traction" , &t_east, 0, -1, true);
    bc.add(0, boundary::west , "Traction" , &t_west, 0, -1, true);    
    // bc.addCondition(0, boundary::east , condition_type::dirichlet, &u_S, 0, true);
    // bc.addCondition(0, boundary::west , condition_type::dirichlet, &u_S, 0, true);

    bc.add(0, boundary::south,"Dirichlet", &phi_D , 1, -1, true);
    bc.add(0, boundary::east, "Neumann"  , &h_east, 1, -1, true);
    bc.add(0, boundary::west, "Neumann"  , &h_west, 1, -1, true);
    // bc.addCondition(0, boundary::east , condition_type::dirichlet, &phi_D, 1, true);
    // bc.addCondition(0, boundary::west , condition_type::dirichlet, &phi_D, 1, true);

    //todo remove all dofs strictly inside the fictitious domain
    //the system still can be solved since the phase field solutions we used are no exactly 0 in the fictitious domain
    
    // Setup system assembler
    gsExprAssembler<real_t> assembler(3,3);
    typedef gsExprAssembler<real_t>::variable variable;
    typedef gsExprAssembler<real_t>::space    space;
    typedef gsExprAssembler<real_t>::solution solution;
    space    u_N = assembler.getSpace(sp3, 2, 0); // args: splines, dim, id
    u_N.addBc( bc.get("Dirichlet",0));
    u_N.setInterfaceCont(0);
    space    f_D = assembler.getSpace(sp3, 1, 1);
    f_D.addBc( bc.get("Dirichlet",1));
    f_D.setInterfaceCont(0);
    space    p_N = assembler.getSpace(sp2, 1, 2); // args: splines, dim, id
    variable f_N = assembler.getCoeff(f        ); // force
    variable ph  = assembler.getCoeff(c        ); // phase field
    variable dph = assembler.getCoeff(dc       ); // explicit phase field gradient
    variable t1  = assembler.getCoeff(tau      ); // tangent vector
    variable Ind = assembler.getCoeff(dind     ); // domain indicator
    
    assembler.setIntegrationElements(subEl); // subdivided mesh
    assembler.initSystem();
    gsInfo<<"Matrix blocks:\n"<< assembler.matrixBlockView() <<"\n";

    gsMatrix<> sol;
    /*
    sol.resizeLike(assembler.rhs());
    sol.setZero();
    solution uu = assembler.getSolution(u_N, sol);
    uu.insert(spdomain);
    */

    // Newmann BC data
    variable gg = assembler.getBdrFunction(); // binds to functions passed to BCs

    // (28a)
    assembler.assemble(
        (nu/2) * ph.val() * ( jac(u_N)+jac(u_N).tr() ) % ( jac(u_N)+jac(u_N).tr() ) ,
        u_N * f_N * ph.val() );
    assembler.assemble( - div(u_N) * p_N.tr() * ph.val() );
    //assembler.assembleLhs( - g * u_N * grad(ph).tr() * f_D.tr()  );
    assembler.assemble( - g * (u_N * dph) * f_D.tr() );
    //assembler.assembleLhs( (u_N*t1) * (u_N*t1).tr() * grad(ph).norm()  );
    assembler.assemble( (u_N*t1) * (u_N*t1).tr() * dph.norm()  );
    assembler.assemble(bc.get("Traction",0), u_N * gg * Ind.val() );
    // (28b)
    assembler.assemble( p_N * div(u_N).tr() * ph.val() );
    // (28c)
    assembler.assemble( grad(f_D) * grad(f_D).tr()  * ( 1 - ph.val() ) );
    //assembler.assembleLhs( f_D  * grad(ph) * u_N.tr() );
    assembler.assemble( f_D  * dph.tr() * u_N.tr() );
    assembler.assemble(bc.get("Neumann",1), -f_D * gg.val() * (1-Ind.val()) );
    
    // Solve the linear system
    //gsSparseSolver<>::BiCGSTABDiagonal solver;
    gsSparseSolver<>::QR solver;

    sol = solver.compute(assembler.matrix()).solve( assembler.rhs() );
    solution _uN = assembler.getSolution(u_N, sol);
    solution _fD = assembler.getSolution(f_D, sol);
    solution _pN = assembler.getSolution(p_N, sol);

    gsExprEvaluator<real_t> ev(assembler);
    variable _uS = ev.getVariable(u_S);

    gsGeometry<>::uPtr ub = _uN.extractPiece(0);
    variable _ub = ev.getVariable(*ub);
    
    real_t relL2Norm = math::sqrt( ev.integral( (_uN - _uS).sqNorm() * Ind.val() ) );
    real_t relH1Norm = math::sqrt( ev.integral( (fjac(_ub) - fjac(_uS)).sqNorm() * Ind.val() ) );
    
    relL2Norm = relL2Norm / math::sqrt( ev.integral( _uS.sqNorm() * Ind.val() ) );
    relH1Norm = relH1Norm / math::sqrt( ev.integral( fjac(_uS).sqNorm() * Ind.val() ) );
    
    gsInfo<< "* relative L2 error : " << relL2Norm << "\n";
    gsInfo<< "* relative H1 error : " << relH1Norm << "\n";
    
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsField<> velocity(domain, u_S);
        gsWriteParaview<>(velocity, "velocity", 1000);
        // gsField<> pressure(domain, p_S);
        // gsWriteParaview<>(pressure, "pressure", 1000);
        // gsField<> piezometric(domain, phi_D);
        // gsWriteParaview<>(piezometric, "piezometric", 1000);        
        gsField<> phaseField(domain, c);
        gsWriteParaview<>(phaseField, "phaseField", 1000, true);

        // plot computed solution
        ev.writeParaview( _uN * Ind.val() - (1-Ind.val())*grad(_fD).tr(), "velocity_sol");

        // plot error
        ev.writeParaview( ( (_uN - _uS) * Ind.val() ).norm() / _uS.norm() , "velocity_err");

        // Extract solutions
        // gsGeometry<>::uPtr uN = _uN.extract(0);
        // gsGeometry<>::uPtr fD = _fD.extract(0);
        
        gsWriteParaview<>(subEl, domain, "subEl_mesh", 100);
        
        gsWriteParaview(sp3.basis(0), "BasisFcn");

        // Run paraview
        return system("paraview velocity_sol.pvd &");
    }

    return EXIT_SUCCESS;
}

void SolutionEvaluator(gsBasis<> &domainBasis, const gsMatrix<> &coef,
                       const gsMatrix<> &x, gsMatrix<> &val)
{
    // x: dim by nPt matrix contains the coordinates of points to be evulated
    // coef: nDofs by nComp matrix contains the solution coefficients
    // val: nPt by nComp matrix contains the evulated results
    
    int nPt = x.cols();
    int nComp = coef.cols();
    int nAct;
    
    val.setZero(nPt, nComp);
    
    gsMatrix<index_t> activeBasisId; // nAct by nPt
    gsMatrix<real_t> activeBasisVal; // nAct by nPt
    
    domainBasis.active_into(x, activeBasisId);
    domainBasis.eval_into(x, activeBasisVal);
    
    nAct = activeBasisId.rows();
    
    for (int i = 0; i < nPt; ++i)
    {
        for (int j = 0; j < nAct; ++j)
        {
            val.row(i) += activeBasisVal(j,i)*coef.row(activeBasisId(j,i));
        }
    }
    
}

int CalElmType(gsFunctionExpr<>& disFunc, gsMatrix<> box)
{
    // _disFunc: distance function
    // _box: dim-by-2 matrix containing the lower & upper corners of a box
    int dim = box.rows(), elmType = -100;
    real_t eps = 1e-10;
    gsMatrix<> BOX;
    
    if (dim == 2)
    {
        BOX.resize(dim, 4);
        BOX.col(0) = box.col(0);
        BOX.col(1) = box.col(1);
        
        BOX(0,2) = box(0,1);
        BOX(1,2) = box(1,0);
        
        BOX(0,3) = box(0,0);
        BOX(1,3) = box(1,1);
    }
    else if (dim == 3)
    {
        BOX.resize(dim, 8);
        BOX.col(0) = box.col(0);
        BOX.col(1) = box.col(1);
        
        BOX(0,2) = box(0,0);
        BOX(1,2) = box(1,1);
        BOX(2,2) = box(2,2);
        
        BOX(0,3) = box(0,1);
        BOX(1,3) = box(1,2);
        BOX(2,3) = box(2,0);
        
        BOX(0,4) = box(0,2);
        BOX(1,4) = box(1,0);
        BOX(2,4) = box(2,1);
        
        BOX(0,5) = box(0,2);
        BOX(1,5) = box(1,1);
        BOX(2,5) = box(2,0);
        
        BOX(0,6) = box(0,1);
        BOX(1,6) = box(1,0);
        BOX(2,6) = box(2,2);
        
        BOX(0,7) = box(0,0);
        BOX(1,7) = box(1,2);
        BOX(2,7) = box(2,1);
    }
    
    gsMatrix<> dis = disFunc.eval(BOX);
    
    real_t temp = dis.cwiseSign().sum();
    
    if (dim == 2)
    {
        if (temp == 4)
            elmType =  1;
        else if (temp == -4)
            elmType =  -1;
        else
            elmType = 0;
    }
    else if (dim == 3)
    {
        if (temp == 8)
            elmType = 1;
        else if (temp == -8)
            elmType = -1;
        else
            elmType = 0;
    }
    
    //If any vertex is on the interface and the others are in the fictitious domain
    for (int i = 0; i < dis.cols(); ++i)
    {
        if ( math::abs(dis(i)) < eps )
        {
            elmType = 0;
//            if (dis.sum() < 0)
//            {
//                elmType = -1;
//            }
        }
    }
    //todo: modify the algorithm, may have exception
    
    return elmType;
}

void CalActiveDofs(gsBasis<> &Basis, gsMultiBasis<> &Mesh, std::vector<int> &elmType, std::vector<bool> &dofStatus)
{
    gsBasis<>::domainIter domIt = Mesh.basis(0).makeDomainIterator();
    gsMatrix<index_t> activeFuncs;
    
    dofStatus.assign(Basis.size(), false);
    
    int cnt = 0;
    
    for (; domIt->good(); domIt->next())
    {
        if (elmType[cnt] >= 0)
        {
            Basis.active_into(domIt->center, activeFuncs);
            
            for (int i = 0; i < activeFuncs.rows(); ++i)
            {
                if (dofStatus[activeFuncs(i,0)] == false)
                {
                    dofStatus[activeFuncs(i,0)] = true;
                }
            }
        }
        cnt++;
    }
}

bool AdaptiveRefinement(gsMultiBasis<>& Basis, gsFunctionExpr<>& disFunc, real_t eps)
{
    int dim = Basis.dim();
    int patchId = 0;
    
    gsMatrix<> refBoxes, element;
    std::vector<real_t> boxX, boxY;
    element.resize(2, 2);
    
    real_t minElmSize = Basis.basis(0).getMinCellLength();
    int level = 1;
    
    // refinement loop
    while (minElmSize > eps)
    {
        gsBasis<>::domainIter domIt = Basis.basis(patchId).makeDomainIterator();
        
        int numMarked = 0;
        for (; domIt->good(); domIt->next())
        {
            element.col(0) = domIt->lowerCorner();
            element.col(1) = domIt->upperCorner();
            
            // Identify the intersected elements, can be replaced by other rules
            if ( CalElmType(disFunc, element) == 0 )
            {
                boxX.push_back( domIt->lowerCorner()(0) );
                boxX.push_back( domIt->upperCorner()(0) );
                boxY.push_back( domIt->lowerCorner()(1) );
                boxY.push_back( domIt->upperCorner()(1) );
                
                numMarked++;
            }
        }
        
        // creat box
        refBoxes.resize(dim, int(boxX.size()));
        for (size_t i = 0; i < boxX.size(); ++i)
        {
            refBoxes(0,i) = boxX[i];
            refBoxes(1,i) = boxY[i];
        }
        
        // refine mesh
        Basis.refine(patchId,refBoxes);
        
        refBoxes.clear();
        boxX.clear();
        boxY.clear();
        
        minElmSize = Basis.basis(0).getMinCellLength();
        
        gsInfo << "level : " << level << " " << "numMarked " << numMarked << ". minElmSize: " << minElmSize << "\n";
        
        level ++;
    }
    // refinement loop
    
    return true;
}
