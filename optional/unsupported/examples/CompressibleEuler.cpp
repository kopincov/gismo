
// This file is not a part of the library!

// NOT IMPLEMENTED YET

// Implementation of SSP RK method for time dependent convdiff problems, AFC stabilization used

// Author: A.Jaeschke, parts copied from ConvDiff.cpp
#define _USE_MATH_DEFINES

#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>
#include <gsAssembler/gsCDRAssembler.h>
#include <gsAssembler/gsGenericAssembler.h>
#include <iomanip>


using namespace gismo;
gsSparseMatrix<real_t> AddABlock(int sizei, int sizej, int ibegin,
                                 int jbegin, gsSparseMatrix<real_t> A);
gsMatrix<real_t> A_jx(gsMatrix<real_t> u);
gsMatrix<real_t> A_jy(gsMatrix<real_t> u);
gsMatrix<real_t> F_x(gsMatrix<real_t> u);
gsMatrix<real_t> F_y(gsMatrix<real_t> u);
gsMatrix<real_t> UtoW (gsMatrix<real_t> n, gsMatrix<real_t> u);
gsMatrix<real_t> WtoU (gsMatrix<real_t> n, gsMatrix<real_t> w);
gsMatrix<real_t> nA (gsMatrix<real_t> n, gsMatrix<real_t> u, gsMatrix<real_t> u_inf);
gsMatrix<real_t> F_n_(gsMatrix<real_t> n, gsMatrix<real_t> u, gsMatrix<real_t> u_inf);
gsMatrix<real_t> F_n_wall (gsMatrix<real_t> n, gsMatrix<real_t> u);
gsMatrix<real_t> F_n_subin (gsMatrix<real_t> n, gsMatrix<real_t> u,
                            real_t rho, real_t p, real_t v_tangential);
gsMatrix<real_t> F_n_subout (gsMatrix<real_t> n, gsMatrix<real_t> u, real_t p_out);
gsMatrix<real_t> F_n_supin (gsMatrix<real_t> n, real_t rho, real_t v_x,
                            real_t v_y, real_t E);
gsMatrix<real_t> F_n_supout (gsMatrix<real_t> n, gsMatrix<real_t> u);
gsMatrix<real_t> F_boundary (gsMultiPatch<real_t> mp, gsMultiBasis<> bases, gsField<> & U1,
                             gsField<> & U2, gsField<> & U3, gsField<> & U4);
gsMatrix<real_t> ConstL2(gsMatrix <real_t> INrhs, gsSparseMatrix <real_t> M, gsSparseMatrix <real_t> Ml );
gsMatrix<real_t> outer_normal(gsMatrix<real_t> points, gsMultiPatch<real_t> mp, boxSide a);
gsMatrix<real_t> mapping(gsMatrix<real_t> points, gsMultiPatch<real_t> mp);

int flag = 0;

int main(int argc, char *argv[])
{
    // Input options
    index_t numElevate  = 0;
    index_t numHref     = 0;
    index_t basisDegree = 0;
    bool plot       = false;
    index_t ntsteps = 300;
    real_t tau = 0.001;
    // Flag for SUPG-Stabilization
    //bool Flag_Stabilization = 0;
    //Flag_Stabilization = 1;

    gsCmdLine cmd("Testing compressible Euler problem.");
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("p","degree",
               "Degree of the basis functions to use for solving (will elevate or reduce the input)",
               basisDegree);
    cmd.addInt("e","degreeElevation",
               "Number of degree elevation steps to perform on the Geometry's basis before solving",
               numElevate);
    cmd.addInt("n","nTimesteps", "Number of time steps", ntsteps);
    cmd.addReal("t","Timestep","Time step size",tau);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFunctionExpr<real_t>  f("0.0", 2);
    gsFunctionExpr<real_t>  g("0.0", 2);
    gsFunctionExpr<real_t>  coeff_c("0", 2);
    gsFunctionExpr<real_t>  coeff_cM("1.0", 2);
    gsFunctionExpr<real_t>  coeff_AS("0.0","0","0","0.0", 2);
    gsFunctionExpr<real_t>  coeff_bS("0","0", 2);
    gsFunctionExpr<real_t>  coeff_AK("0","0","0","0", 2);
    gsFunctionExpr<real_t>  coeff_bCx("1.0","0.0", 2);
    gsFunctionExpr<real_t>  coeff_bCy("0.0","1.0", 2); 
    //Test of normals
    gsMultiPatch<real_t> mp;
    //gsReadFile<>("planar/testofnormals.xml", mp);
    gsReadFile<>("planar/square.xml", mp);
    //B1 & B2
    //gsMultiPatch<real_t> mp( *gsNurbsCreator<>::BSplineRectangle(-10.0,-10.0,10.0,10.0) );
    //B3
    //gsMultiPatch<real_t> mp( *gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.0) );

    gsMultiBasis<> bases(mp);

    if (basisDegree)
        bases.setDegree(basisDegree);
    else if (numElevate)
        bases.degreeElevate(numElevate);

    gsSparseSolver<>::LU solver;
    for (int i=0; i<numHref ; i++)
        bases.uniformRefine();


    // Basic version, all parameters are NOT time dependent
    gsBoundaryConditions<> BCs;
    //gsCDRAssembler<real_t> galerkinS(mp,bases,BCs,f,coeff_AS,coeff_bS,coeff_c,
    //  dirichlet::none, iFace::glue, 0);
    //galerkinS.assemble();
    gsCDRAssembler<real_t> galerkinCx(mp,bases,BCs,f,coeff_AK,coeff_bCx,coeff_c,
                                      dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinCx.assemble();
    gsCDRAssembler<real_t> galerkinCy(mp,bases,BCs,f,coeff_AK,coeff_bCy,coeff_c,
                                      dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinCy.assemble();
    gsCDRAssembler<real_t> galerkinM(mp,bases,BCs,f,coeff_AK,coeff_bS,coeff_cM,
                                     dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinM.assemble();

    int basisSize = galerkinM.matrix().innerSize();
    gsSparseMatrix<real_t> M (4*basisSize, 4*basisSize);
    M = AddABlock(4*basisSize, 4*basisSize, 0,0,galerkinM.matrix()) + AddABlock(4*basisSize, 4*basisSize, basisSize,basisSize ,galerkinM.matrix()) +
        AddABlock(4*basisSize, 4*basisSize, 2*basisSize, 2* basisSize ,galerkinM.matrix()) + AddABlock(4*basisSize, 4*basisSize, 3*basisSize, 3* basisSize ,galerkinM.matrix());
    gsSparseMatrix<real_t> K (4*basisSize, 4*basisSize);

    gsSparseMatrix<real_t> Ml_ (galerkinM.matrix());

    ///*
    for (int e=0; e<Ml_.outerSize(); ++e)
        for (gsSparseMatrix<real_t>::InnerIterator it(Ml_,e); it; ++it)
            {
                int k = it.row();
                int l = it.col();
                if((Ml_.at(k,l) != 0)&& k!=l)
                    {
                        Ml_.coeffRef(k,k) = Ml_.at(k,k) + Ml_.at(k,l);
                        Ml_.coeffRef(k,l) = 0.0;
                    }
            }

    // Initial conditions
    //B1&2
    //gsFunctionExpr<real_t>  init1("pow(1.0-0.09046534*exp(1.0-(x*x+y*y)),2.5)", 2);
    //B3
    gsFunctionExpr<real_t>  init1("if(x<0.5,1.0,0.125)", 2);
    gsCDRAssembler<real_t> galerkinIN1(mp,bases,BCs,init1,coeff_AK,coeff_bS,coeff_cM,
                                       dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinIN1.assemble();
    gsMatrix<real_t> u_init1 = ConstL2(galerkinIN1.rhs(), galerkinM.matrix(), Ml_);

    // B1
    //gsFunctionExpr<real_t>  init2("pow(1.0-0.09046534*exp(1.0-(x*x+y*y)),2.5)*(0.0-0.79577472*y*exp(0.5*(1.0-(x*x+y*y))))", 2);
    //B2
    //gsFunctionExpr<real_t>  init2("pow(1.0-0.09046534*exp(1.0-(x*x+y*y)),2.5)*(5.0-0.79577472*y*exp(0.5*(1.0-(x*x+y*y))))", 2);
    //B3
    gsFunctionExpr<real_t>  init2("0.0", 2);
    gsCDRAssembler<real_t> galerkinIN2(mp,bases,BCs,init2,coeff_AK,coeff_bS,coeff_cM,
                                       dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinIN2.assemble();
    gsMatrix<real_t> u_init2 = ConstL2(galerkinIN2.rhs(), galerkinM.matrix(), Ml_);

    // B1&2
    //gsFunctionExpr<real_t>  init3("pow(1.0-0.09046534*exp(1.0-(x*x+y*y)),2.5)*(0.0+0.79577472*x*exp(0.5*(1.0-(x*x+y*y))))", 2);
    //B3
    gsFunctionExpr<real_t>  init3("0.0", 2);
    gsCDRAssembler<real_t> galerkinIN3(mp,bases,BCs,init3,coeff_AK,coeff_bS,coeff_cM,
                                       dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinIN3.assemble();
    gsMatrix<real_t> u_init3 = ConstL2(galerkinIN3.rhs(), galerkinM.matrix(), Ml_);
    //B1
    //gsFunctionExpr<real_t>  init4("pow(1.0-0.09046534*exp(1.0-(x*x+y*y)),3.5)/0.4+pow(1.0-0.09046534*exp(1.0-(x*x+y*y)),2.5)*(((0.0-0.79577472*y*exp(0.5*(1.0-(x*x+y*y))))^2.0+(0.0+0.79577472*x*exp(0.5*(1.0-(x*x+y*y))))^2.0)/2.0)", 2);
    //B2
    //gsFunctionExpr<real_t>  init4("pow(1.0-0.09046534*exp(1.0-(x*x+y*y)),3.5)/0.4+pow(1.0-0.09046534*exp(1.0-(x*x+y*y)),2.5)*(((5.0-0.79577472*y*exp(0.5*(1.0-(x*x+y*y))))^2.0+(0.0+0.79577472*x*exp(0.5*(1.0-(x*x+y*y))))^2.0)/2.0)", 2);
    // B3
    gsFunctionExpr<real_t>  init4("if(x<0.5,1.0/0.4,0.1/0.4)", 2);
    gsCDRAssembler<real_t> galerkinIN4(mp,bases,BCs,init4,coeff_AK,coeff_bS,coeff_cM,
                                       dirichlet::none, iFace::glue, stabilizerCDR::none);
    galerkinIN4.assemble();
    gsMatrix<real_t> u_init4 = ConstL2(galerkinIN4.rhs(), galerkinM.matrix(), Ml_);
    ///*
    // for sure someting is wrong with u passed t the boundary conditions
    if (plot)
        {           //May be an issue with u_init4
            gsField<> sol4 = galerkinM.constructSolution(u_init1);
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>( sol4, "aaaINIT", 10000);//1000000);
        }
    gsMatrix <real_t> u_init (4*basisSize, 1);
    for (int k=0; k<4*basisSize; k++)
        {
            if (k<basisSize)
                u_init.coeffRef(k) = u_init1.at(k);
            else if (k<2*basisSize)
                u_init.coeffRef(k) = u_init2.at(k-basisSize);
            else if (k<3*basisSize)
                u_init.coeffRef(k) = u_init3.at(k-2*basisSize);
            else
                u_init.coeffRef(k) = u_init4.at(k-3*basisSize);
        }
    //*/
    /*
      gsMatrix <real_t> u_init (4*basisSize, 1);
      for (int k=0;k<4*basisSize; k++)
      {
      u_init(k) = 1.0;
      if (k>= bases.totalSize() && k< 3*bases.totalSize())
      u_init(k) = 0.0;
      }//*/
    //Time stepping loop
    solver.compute(M);
    gsParaviewCollection collection("aaaTIMELrho");
    gsParaviewCollection collection2("aaaTIMELv1");
    gsParaviewCollection collection3("aaaTIMELv2");
    int orderFlag = 2;
    gsMatrix<real_t> sol (u_init);
    if (orderFlag == 2)
        {
            gsMatrix<real_t> rhs_it1 (galerkinCx.rhs());
            gsMatrix<real_t> rhs_it2 (galerkinCx.rhs());
            gsMatrix<real_t> u1 (u_init);
            //int bases1d = bases[0].component(0)->size();
            //int bases2d = bases[0].component(1)->size();
            for (int i=0; i< ntsteps; i++)
                {

                    //Building up K
                    std::vector<Eigen::Triplet<real_t> > tripletList;
                    tripletList.reserve(galerkinCx.matrix().nonZeros());
                    for (int e=0; e<galerkinCx.matrix().outerSize(); ++e)
                        for (gsSparseMatrix<real_t>::InnerIterator it(galerkinCx.matrix(),e); it; ++it)
                            {
                                //int l = it.row();
                                //int k = it.col();
                                int k = it.row();
                                int l = it.col();
                                //fix solution
                                gsMatrix<real_t> solat (4,1);
                                solat << sol.at(l),sol.at(l+basisSize),sol.at(l+2*basisSize),sol.at(l+3*basisSize);
                                gsMatrix<real_t> k_kl = galerkinCx.matrix().at(l,k)*A_jx(solat) +
                                    galerkinCy.matrix().at(l,k)*A_jy(solat);
                                for (int o=0;o<4;o++)
                                    for (int oo=0;oo<4;oo++)
                                        tripletList.push_back(Eigen::Triplet<real_t> (k+o*basisSize,l+oo*basisSize,k_kl.coeff(o,oo)));
                            }
                    K.setFromTriplets(tripletList.begin(), tripletList.end());
                    //BCs in function:
                    gsMatrix <real_t> sol1 (bases.totalSize(),1);
                    gsMatrix <real_t> sol2 (bases.totalSize(),1);
                    gsMatrix <real_t> sol3 (bases.totalSize(),1);
                    gsMatrix <real_t> sol4 (bases.totalSize(),1);
                    for (int j=0;j<int(bases.totalSize());j++)
                        {
                            sol1.coeffRef(j) = sol.at(j);
                            sol2.coeffRef(j) = sol.at(j+bases.totalSize());
                            sol3.coeffRef(j) = sol.at(j+2*bases.totalSize());
                            sol4.coeffRef(j) = sol.at(j+3*bases.totalSize());
                        }
                    gsField<> U1 = galerkinM.constructSolution(sol1);
                    gsField<> U2 = galerkinM.constructSolution(sol2);
                    gsField<> U3 = galerkinM.constructSolution(sol3);
                    gsField<> U4 = galerkinM.constructSolution(sol4);

                    gsMatrix<real_t> S = F_boundary (mp, bases, U1, U2, U3, U4);
                    ///*
                    rhs_it1 = M*sol + tau*K*sol -tau*S;
                    u1 = solver.solve(rhs_it1);
                    sol = u1;
                    //rhs_it2 = 0.5*M*sol + 0.5*M*u1 + 0.5*tau*K*u1 - 0.5*tau*S;
                    //sol = solver.solve(rhs_it2);//*/
                    //rhs_it1 = M*sol + tau*K*sol - tau*S;
                    //sol = solver.solve(rhs_it1);
                    // plotting in ParaView

                    ///////FIX this later on
                    gsInfo << i << "/" << ntsteps << "\n";
                    ///*
                    for (int j=0;j<int(bases.totalSize());j++)
                        {
                            sol1.coeffRef(j) = sol.at(j);
                            sol2.coeffRef(j) = sol.at(j+bases.totalSize())/sol.at(j);
                            sol3.coeffRef(j) = sol.at(j+2*bases.totalSize())/sol.at(j);
                            //sol4.coeffRef(j) = sol.at(j+3*bases.totalSize());
                        }
                    /*gsInfo << "Solution" <<
                      "\n" <<
                      sol2 << "\n";*/
                    if (plot)
                        {
                            gsField<> U11 = galerkinM.constructSolution(sol1);
                            gsField<> U21 = galerkinM.constructSolution(sol2);
                            gsField<> U31 = galerkinM.constructSolution(sol3);
                            //gsField<> U41 = galerkinM.constructSolution(sol4);
                            std::string fileName = "aaaTIMErho" + util::to_string(i);
                            gsWriteParaview<>( U11, fileName, 10000);
                            fileName = fileName + "0";
                            collection.addTimestep(fileName,i,".vts");
                            /*fileName = "aaaTIMEv1" + util::to_string(i);
                              gsWriteParaview<>( U21, fileName, 10000);
                              fileName = fileName + "0";
                              collection2.addTimestep(fileName,i,".vts");
                              fileName = "aaaTIMEv2" + util::to_string(i);
                              gsWriteParaview<>( U21, fileName, 10000);
                              fileName = fileName + "0";
                              collection3.addTimestep(fileName,i,".vts");*/
                        }
                }

        }



    else if (orderFlag == 3)
        {
            gsMatrix<real_t> u1 (u_init);
            gsMatrix<real_t> u2 (u_init);
            for (int i=0; i< ntsteps; i++)
                {
                    //u1 = sol - tau*M_inv*K*sol + tau*M_inv*rhs;
                    //u2 = 0.75*sol + 0.25*u1 - 0.25*tau*M_inv*K*u1 + 0.25*tau*M_inv*rhs;
                    //sol = sol/3.0 + 2.0/3.0*u2 - 2.0/3.0*tau*M_inv*K*u2 + 2.0/3.0*tau*M_inv*rhs;
                }
        }
    else
        {
            gsInfo << "NOT IMPLEMENTED ORDER"<< "\n";
            return -1;
        }


    //Postprocessing
    if (plot)
        {
            collection.save();
            /*collection2.save();
              collection3.save();*/


            gsMatrix <real_t> sol1 (bases.totalSize(),1);
            gsMatrix <real_t> sol2 (bases.totalSize(),1);
            gsMatrix <real_t> sol3 (bases.totalSize(),1);
            gsMatrix <real_t> sol4 (bases.totalSize(),1);
            for (int j=0;j<int(bases.totalSize());j++)
                {
                    sol1.coeffRef(j) = sol.at(j);
                    sol2.coeffRef(j) = sol.at(j+bases.totalSize());
                    sol3.coeffRef(j) = sol.at(j+2*bases.totalSize());
                    sol4.coeffRef(j) = sol.at(j+3*bases.totalSize());
                }
            /*gsInfo << "Solution" <<
              "\n" <<
              sol2 << "\n";*/
            gsField<> U1 = galerkinM.constructSolution(sol1);
            gsField<> U2 = galerkinM.constructSolution(sol2);
            gsField<> U3 = galerkinM.constructSolution(sol3);
            gsField<> U4 = galerkinM.constructSolution(sol4);
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>( U1, "aaa_SOL_rhov1", 10000);
        }



    return 0;
}

gsSparseMatrix<real_t> AddABlock(int sizei, int sizej, int ibegin,
                                 int jbegin, gsSparseMatrix<real_t> A){
    //only for ColMajor Sparse Matrix
    int i,j;
    std::vector<Eigen::Triplet<real_t> > tripletList;
    tripletList.reserve(A.nonZeros());
    for (int e=0; e<A.outerSize(); ++e)
        for (gsSparseMatrix<real_t>::InnerIterator it(A,e); it; ++it)
            {
                i = it.row() + ibegin;
                j = it.col() + jbegin;

                tripletList.push_back(Eigen::Triplet<real_t> (i,j,it.value()));
            }

    gsSparseMatrix<real_t> res(sizei,sizej);
    res.setFromTriplets(tripletList.begin(), tripletList.end());
    return res;
}

gsMatrix<real_t> A_jx(gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,4);
    real_t gamma = 1.4;
    real_t gamma1 = (gamma-1.0)/2.0;
    real_t gamma2 = (gamma-3.0)/2.0;
    real_t v1 = u.at(1)/u.at(0);
    real_t v2 = u.at(2)/u.at(0);
    real_t H = (gamma*u.at(3) + (1.0-gamma)*(u.at(1)*u.at(1)+u.at(2)*u.at(2))/
                (2.0*u.at(0)))/u.at(0);
    res << 0, 1.0, 0, 0,
        gamma2*v1*v1+gamma1*v2*v2, (3.0-gamma)*v1, (1.0-gamma)*v2, gamma-1.0,
        -v1*v2, v2, v1, 0,
        gamma1*(v1*v1*v1+v2*v2*v1)-H*v1, H-(gamma-1.0)*v1*v1, (1.0-gamma)*v1*v2, gamma*v1;
    return res;
}

gsMatrix<real_t> A_jy(gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,4);
    real_t gamma = 1.4;
    real_t gamma1 = (gamma-1.0)/2.0;
    real_t gamma2 = (gamma-3.0)/2.0;
    real_t v1 = u.at(1)/u.at(0);
    real_t v2 = u.at(2)/u.at(0);
    real_t H = (gamma*u.at(3) + (1.0-gamma)*(u.at(1)*u.at(1)+u.at(2)*u.at(2))/(2.0*u.at(0)))/
        u.at(0);
    res << 0, 0, 1, 0,
        -v1*v2, v2, v1, 0,
        gamma2*v2*v2+gamma1*v1*v1, (1.0-gamma)*v1, (3.0-gamma)*v2, gamma-1.0,
        gamma1*(v1*v1*v2+v2*v2*v2)-H*v2, (1.0-gamma)*v1*v2, H-(gamma-1.0)*v2*v2, gamma*v2;
    return res;
}

gsMatrix<real_t> F_x(gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    real_t p = (gamma-1.0)*(u.at(3) - (u.at(1)*u.at(1)+u.at(2)*u.at(2))/(2.0*u.at(0)));
    res << u.at(1),
        u.at(1)*u.at(1)/u.at(0) + p,
        u.at(1)*u.at(2)/u.at(0),
        u.at(3)*u.at(1)/u.at(0) + p*u.at(1)/u.at(0);
    return res;
}

gsMatrix<real_t> F_y(gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    real_t p = (gamma-1.0)*(u.at(3) - (u.at(1)*u.at(1)+u.at(2)*u.at(2))/(2.0*u.at(0)));
    res << u.at(2),
        u.at(1)*u.at(2)/u.at(0),
        u.at(2)*u.at(2)/u.at(0)+p,
        u.at(3)*u.at(2)/u.at(0) + p*u.at(2)/u.at(0);
    return res;
}

gsMatrix<real_t> UtoW (gsMatrix<real_t> n, gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    real_t p = (gamma-1.0)*(u.at(3) - (u.at(1)*u.at(1)+u.at(2)*u.at(2))/(2.0*u.at(0)));
    real_t vn = n.at(0)*u.at(1)/u.at(0) + n.at(1)*u.at(2)/u.at(0);
    real_t vxi = n.at(1)*u.at(1)/u.at(0) - n.at(0)*u.at(2)/u.at(0);
    real_t c = math::sqrt(gamma*p/u.at(0));
    res << vn - (2.0*c)/(gamma-1.0),
        p/math::pow(u.at(0),gamma),
        vxi,
        vn + (2.0*c)/(gamma-1.0);
    return res;
}

gsMatrix<real_t> WtoU (gsMatrix<real_t> n, gsMatrix<real_t> w){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    real_t vn = (w.at(3)+w.at(0))/2.0;
    real_t vxi = w.at(2);
    real_t c = (gamma-1.0)/4.0*(w.at(3)-w.at(0));
    real_t rho = math::pow(c*c/gamma/w.at(1),1.0/(gamma-1.0));
    //real_t rho = math::pow(c*c/gamma*exp(-w.at(1)/c_v),1.0/(gamma-1.0));
    real_t p = rho*c*c/gamma;
    res << rho,
        rho*(vn*n.at(0)+vxi*n.at(1)),
        rho*(vn*n.at(1)-vxi*n.at(0)),
        p/(gamma-1.0)+rho/2.0*(vn*vn+vxi*vxi);//rho/2.0*(vn*vn+vxi*vxi); //not sure
    return res;
}

gsMatrix<real_t> nA (gsMatrix<real_t> n, gsMatrix<real_t> u, gsMatrix<real_t> u_inf){
    // Assumes n is unit vector
    gsMatrix<real_t> res (4,4);
    real_t gamma = 1.4;
    real_t v_ix = u.at(1)/u.at(0);
    real_t v_iy = u.at(2)/u.at(0);
    real_t v_jx = u_inf.at(1)/u_inf.at(0);
    real_t v_jy = u_inf.at(2)/u_inf.at(0);
    real_t v_ijx = (math::sqrt(u.at(0))*v_ix + math::sqrt(u_inf.at(0))*v_jx)/
        (math::sqrt(u.at(0)) + math::sqrt(u_inf.at(0)));
    real_t v_ijy = (math::sqrt(u.at(0))*v_iy + math::sqrt(u_inf.at(0))*v_jy)/
        (math::sqrt(u.at(0)) + math::sqrt(u_inf.at(0)));
    real_t v_ijn = n.at(0)*v_ijx + n.at(1)*v_ijy;
    real_t H_i = (gamma*u.at(3) + (1.0-gamma)*(u.at(1)*u.at(1)+u.at(2)*u.at(2))/
                  (2.0*u.at(0)))/u.at(0);
    real_t H_j = (gamma*u_inf.at(3) + (1.0-gamma)*
                  (u_inf.at(1)*u_inf.at(1)+u_inf.at(2)*u_inf.at(2))/(2*u_inf.at(0)))/u_inf.at(0);
    real_t H_ij = (math::sqrt(u.at(0))*H_i + math::sqrt(u_inf.at(0))*H_j)/
        (math::sqrt(u.at(0)) + math::sqrt(u_inf.at(0)));
    real_t c_ij = math::sqrt((gamma-1.0)*(H_ij-(v_ijx*v_ijx+v_ijy*v_ijy)/2.0));
    real_t q = math::sqrt(v_ijx*v_ijx+v_ijy*v_ijy)/2.0;
    real_t b2 = (gamma-1.0)/(c_ij*c_ij);
    real_t b1 = b2*q;
    gsMatrix<real_t> lambda (4,4);
    lambda << math::abs(v_ijn - c_ij), 0, 0, 0,
        0, math::abs(v_ijn), 0, 0,
        0, 0, math::abs(v_ijn + c_ij), 0,
        0, 0, 0, math::abs(v_ijn);
    gsMatrix<real_t> R(4,4);
    R << 1.0, 1.0, 1.0, 0,
        v_ijx-c_ij*n.at(0), v_ijx, v_ijx+c_ij*n.at(0), n.at(1),
        v_ijy-c_ij*n.at(1), v_ijy, v_ijy+c_ij*n.at(1), -n.at(0),
        H_ij-c_ij*v_ijn, q, H_ij+c_ij*v_ijn, v_ijx*n.at(1)-v_ijy*n.at(0);
    gsMatrix<real_t> L(4,4);
    L << 0.5*(b1+v_ijn/c_ij), 0.5*(-b2*v_ijx-n.at(0)/c_ij), 0.5*(-b2*v_ijy - n.at(1)/c_ij), 0.5*b2,
        1.0-b1, b2*v_ijx, b2*v_ijy, -b2,
        0.5*(b1 - v_ijn/c_ij), 0.5*(-b2*v_ijx+n.at(0)/c_ij), 0.5*(-b2*v_ijy+ n.at(1)/c_ij), 0.5*b2,
        n.at(0)*v_ijn-n.at(1)*v_ijx, n.at(1), -n.at(0), 0;
    // Old version of the last row causing problems (v_ijy-v_ijn*n.at(1))/n.at(0), n.at(1), (n.at(1)*n.at(1)-1.0)/n.at(0), 0;
    res = R*lambda*L;
    return res;
}

gsMatrix<real_t> F_n_ (gsMatrix<real_t> n, gsMatrix<real_t> u, gsMatrix<real_t> u_inf){
    gsMatrix<real_t> res (4,1);
    res = n.at(0)*(F_x(u)+F_x(u_inf))/2.0+n.at(1)*(F_y(u)+F_y(u_inf))/2.0
        - 0.5*nA(n,u,u_inf)*(u_inf-u);
    return res;
}

gsMatrix<real_t> F_n_wall (gsMatrix<real_t> n, gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,1);
    gsMatrix<real_t> u_inf (4,1);
    real_t vx = u.at(1)/u.at(0);
    real_t vy = u.at(2)/u.at(0);
    u_inf << u.at(0),
        u.at(0)*(vx - 2*n.at(0)*(n.at(0)*vx+n.at(1)*vy)),
        u.at(0)*(vy - 2*n.at(1)*(n.at(0)*vx+n.at(1)*vy)),
        u.at(3);
    res = F_n_(n,u,u_inf);
    return res;
}

gsMatrix<real_t> F_n_subin (gsMatrix<real_t> n, gsMatrix<real_t> u,
                            real_t rho, real_t p, real_t v_tangential){
    real_t gamma = 1.4;
    gsMatrix<real_t> res (4,1);
    gsMatrix<real_t> w = UtoW(n,u);
    w.coeffRef(2) = v_tangential;
    //temporary fix
    //w.coeffRef(1) = 1.0;
    w.coeffRef(1) = p/math::pow(rho,gamma);//c_v*math::log(p/math::pow(rho,gamma));
    w.coeffRef(0) = w.at(3)-4.0/(gamma-1.0)*math::sqrt(gamma*p/rho);
    gsMatrix<real_t> u_inf = WtoU(n,w);
    res = F_n_(n,u,u_inf);
    return res;
}

gsMatrix<real_t> F_n_subout (gsMatrix<real_t> n, gsMatrix<real_t> u, real_t p_out){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    real_t p = (gamma-1.0)*(u.at(3) - (u.at(1)*u.at(1)+u.at(2)*u.at(2))/(2.0*u.at(0)));
    gsMatrix<real_t> w = UtoW(n,u);
    //real_t rho_out;
    //rho_out = u.at(0)*(math::pow(p_out/p,1/gamma));
    w.coeffRef(0) = w.at(3) - 4.0/(gamma-1)*math::sqrt(gamma*p_out/u.at(0)*math::pow(p/p_out,1/gamma)) ;
    gsMatrix<real_t> u_inf = WtoU(n,w);
    res = F_n_(n,u,u_inf);
    return res;
}

gsMatrix<real_t> F_n_supin (gsMatrix<real_t> n, real_t rho, real_t v_x,
                            real_t v_y, real_t p_in){
    gsMatrix<real_t> res (4,1);
    real_t gamma = 1.4;
    gsMatrix<real_t> u_inf (4,1);
    //real_t v_n = n.at(0)*v_x + n.at(1)*v_y;
    //real_t v_t = n.at(1)*v_x - n.at(0)*v_y;
    real_t rhoE = p_in/(gamma-1.0)+rho*(v_x*v_x+v_y*v_y)/2.0;
    u_inf << rho,
        rho*v_x,
        rho*v_y,
        rhoE;
    res = F_n_(n,u_inf,u_inf);
    return res;
}

gsMatrix<real_t> F_n_supout (gsMatrix<real_t> n, gsMatrix<real_t> u){
    gsMatrix<real_t> res (4,1);
    res = F_n_(n,u,u);
    return res;
}
///*
gsMatrix<real_t> F_boundary (gsMultiPatch<real_t> mp, gsMultiBasis<> bases, gsField<> & U1,
                             gsField<> & U2, gsField<> & U3, gsField<> & U4)
{
    gsMatrix <real_t> res (4*bases.totalSize(),1);
    for (int k = 0; k < 4*int(bases.totalSize()); k++)
        res(k) = 0.0;
    int p = bases.degree();
    gsMatrix <real_t> x (14,1);
    x << -0.9862838086968121,
        -0.9284348836635741,
        -0.8272013150697653,
        -0.6872929048116853,
        -0.5152486363581548,
        -0.3191123689278905,
        -0.1080549487073448,
        0.108054948707344,
        0.319112368927890,
        0.515248636358154,
        0.687292904811686,
        0.827201315069765,
        0.928434883663574,
        0.986283808696812 ;
    gsMatrix <real_t> w (14,1);
    w << 0.035119460331751,
        0.080158087159761,
        0.121518570687903,
        0.157203167158193,
        0.185538397477938,
        0.205198463721295,
        0.215263853463158,
        0.215263853463158,
        0.205198463721296,
        0.185538397477938,
        0.157203167158194,
        0.121518570687903,
        0.080158087159761,
        0.035119460331751;
    //west boundary - subsonic inlet
    {
        //real_t rho = 1.0;
        //real_t press = 10.0;
        //real_t v_tangential = 0.0;
        gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(1);
        for (int i=0; i< bd_basis->size(); i++)              //Disabled: Not repeating corners
            {
                gsMatrix <real_t> interval = bd_basis->function(i).support();
                real_t subint = (interval.at(1) - interval.at(0))/(p+1.0);
                gsMatrix <real_t> pres (4,1);
                pres << 0.0, 0.0, 0.0, 0.0;
                for (int k=0; k<p+1;k++)
                    {
                        gsMatrix <real_t> pres0 (4,1);
                        pres0 << 0.0, 0.0, 0.0, 0.0;
                        gsMatrix <real_t> quad (2,14);
                        for (int l=0; l<14;l++)
                            {
                                quad.coeffRef(1,l) = ((interval.at(0)+subint*k)*(1.0-x.at(l))+
                                                      (interval.at(0)+subint*(k+1.0))*(1.0+x.at(l)))/2.0; //Maybe an error here
                                quad.coeffRef(0,l) = 0.0;
                            }
                        gsMatrix <real_t> u (4,14);
                        gsMatrix <real_t> u_val;
                        u_val = U1.value(quad);
                        u.row(0) = u_val.row(0);
                        u_val = U2.value(quad);
                        u.row(1) = u_val.row(0);
                        u_val = U3.value(quad);
                        u.row(2) = u_val.row(0);
                        u_val = U4.value(quad);
                        u.row(3) = u_val.row(0);
                        gsMatrix <real_t> normal (2,1);
                        normal = outer_normal(quad, mp, 1);
                        gsMatrix <real_t> funct = bd_basis->function(i).eval(quad.row(1));
                        gsMatrix <real_t> F_n (4,14);
                        quad = mapping(quad, mp);
                        real_t len_of = 0.0;
                        for (int l=0; l<14; l++)
                            {
                                // B1 & B3 
                                F_n.col(l) = F_n_wall(normal.col(l),u.col(l));
                                // B2 
                                //F_n.col(l) = F_n_supin (normal, u.at(0,l), u.at(1,l)/u.at(0,l), u.at(2,l)/u.at(0,l),
                                //  (1.4-1.0)*(u.at(3,l) - (u.at(1,l)*u.at(1,l)+u.at(2,l)*u.at(2,l))/(2.0*u.at(0,l))));
                            }
                        for (int l=0; l<14; l++)
                            {
                                pres0 = pres0 + w.at(l)*funct.at(l)*F_n.col(l)*0.5;
                            }
                        for (int l=0; l<13; l++)
                            {
                                len_of = len_of + math::sqrt((quad.coeffRef(0,l+1)-quad.coeffRef(0,l))*(quad.coeffRef(0,l+1)-quad.coeffRef(0,l))
                                                             + (quad.coeffRef(1,l+1)-quad.coeffRef(1,l))*(quad.coeffRef(1,l+1)-quad.coeffRef(1,l)));
                            }
                        pres = pres + pres0 * len_of;
                    }
                for (int k=0; k<4; k++)
                    {
                        res.coeffRef(i*bases[0].component(1).size()+k*bases.totalSize()) += pres.at(k);
                    }
            }
    }
    //east boundary - subsonic outlet
    {
        //real_t press = 1.0;
        gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(2);
        for (int i=0; i< bd_basis->size(); i++)          //Disabled: Not repeating corners
            {
                gsMatrix <real_t> interval = bd_basis->function(i).support();
                real_t subint = (interval.at(1) - interval.at(0))/(p+1.0);
                gsMatrix <real_t> pres (4,1);
                pres << 0.0, 0.0, 0.0, 0.0;
                for (int k=0; k<p+1;k++)
                    {
                        gsMatrix <real_t> pres0 (4,1);
                        pres0 << 0.0, 0.0, 0.0, 0.0;
                        gsMatrix <real_t> quad (2,14);
                        for (int l=0; l<14;l++)
                            {
                                quad.coeffRef(1,l) = ((interval.at(0)+subint*k)*(1.0-x.at(l))+
                                                      (interval.at(0)+subint*(k+1))*(1.0+x.at(l)))/2.0; //Maybe an error here
                                quad.coeffRef(0,l) = 1.0;
                            }
                        gsMatrix <real_t> quad_mapped = mp.patch(0).eval(quad);
                        gsMatrix <real_t> u (4,14);
                        gsMatrix <real_t> u_val;
                        u_val = U1.value(quad);
                        u.row(0) = u_val.row(0);
                        u_val = U2.value(quad);
                        u.row(1) = u_val.row(0);
                        u_val = U3.value(quad);
                        u.row(2) = u_val.row(0);
                        u_val = U4.value(quad);
                        u.row(3) = u_val.row(0);
                        gsMatrix <real_t> normal (2,1);
                        normal = outer_normal(quad, mp, 2);
                        gsMatrix <real_t> funct = bd_basis->function(i).eval(quad.row(1)); // this was not working
                        //gsMatrix <real_t> funct = bases[0].evalSingle((i+1)*bases[0].component(0).size(),quad);
                        gsMatrix <real_t> F_n (4,14);
                        quad = mapping(quad, mp);
                        real_t len_of = 0.0;
                        for (int l=0; l<14; l++)
                            {
                                // B1 & B3
                                F_n.col(l) = F_n_wall(normal.col(l),u.col(l));
                                // B2
                                //F_n.col(l) = F_n_supout(normal, u.col(l));
                            }
                        for (int l=0; l<14; l++)
                            {
                                pres0 = pres0 + w.at(l)*funct.at(l)*F_n.col(l)*0.5;
                            }
                        for (int l=0; l<13; l++)
                            {
                                len_of = len_of + math::sqrt((quad.coeffRef(0,l+1)-quad.coeffRef(0,l))*(quad.coeffRef(0,l+1)-quad.coeffRef(0,l))
                                                             + (quad.coeffRef(1,l+1)-quad.coeffRef(1,l))*(quad.coeffRef(1,l+1)-quad.coeffRef(1,l)));
                            }
                        pres = pres + pres0 * len_of;

                    }
                for (int k=0; k<4; k++)
                    {
                        res.coeffRef((i+1)*bases[0].component(0).size()-1+k*bases.totalSize()) += pres.at(k);
                    }
            }
    }
    //south boundary - wall
    {
        gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(3);
        for (int i=0; i< bd_basis->size(); i++)
            {
                gsMatrix <real_t> interval = bd_basis->function(i).support();
                real_t subint = (interval.at(1) - interval.at(0))/(p+1.0);
                gsMatrix <real_t> pres (4,1);
                pres << 0.0, 0.0, 0.0, 0.0;
                for (int k=0; k<p+1;k++)
                    {
                        gsMatrix <real_t> pres0 (4,1);
                        pres0 << 0.0, 0.0, 0.0, 0.0;
                        gsMatrix <real_t> quad (2,14);
                        for (int l=0; l<14;l++)
                            {
                                quad.coeffRef(0,l) = ((interval.at(0)+subint*k)*(1.0-x.at(l))+
                                                      (interval.at(0)+subint*(k+1))*(1.0+x.at(l)))/2.0; //Maybe an error here
                                quad.coeffRef(1,l) = 0.0;
                            }
                        gsMatrix <real_t> u (4,14);
                        gsMatrix <real_t> u_val;
                        u_val = U1.value(quad);       // EVALUATION TAKES PLACE IN THE PARAMETRIC DOMAIN
                        u.row(0) = u_val.row(0);
                        u_val = U2.value(quad);
                        u.row(1) = u_val.row(0);
                        u_val = U3.value(quad);
                        u.row(2) = u_val.row(0);
                        u_val = U4.value(quad);
                        u.row(3) = u_val.row(0);
                        gsMatrix <real_t> normal (2,1);
                        normal = outer_normal(quad, mp, 3);
                        gsMatrix <real_t> funct = bd_basis->function(i).eval(quad.row(0));
                        gsMatrix <real_t> F_n (4,14);
                        quad = mapping(quad, mp);
                        real_t len_of = 0.0;
                        for (int l=0; l<14; l++)
                            {
                                //B1&2&3
                                F_n.col(l) = F_n_wall(normal.col(l),u.col(l));
                            }
                        for (int l=0; l<14; l++)
                            {
                                pres0 = pres0 + w.at(l)*funct.at(l)*F_n.col(l)*0.5;
                            }
                        for (int l=0; l<13; l++)
                            {
                                len_of = len_of + math::sqrt((quad.coeffRef(0,l+1)-quad.coeffRef(0,l))*(quad.coeffRef(0,l+1)-quad.coeffRef(0,l))
                                                             + (quad.coeffRef(1,l+1)-quad.coeffRef(1,l))*(quad.coeffRef(1,l+1)-quad.coeffRef(1,l)));
                            }
                        pres = pres + pres0 * len_of;
                    }
                for (int k=0; k<4; k++)
                    {
                        res.coeffRef(i+k*bases.totalSize()) += pres.at(k);
                    }
            }
    }

    //north boundary - wall
    {
        gsBasis<>::uPtr bd_basis = bases[0].boundaryBasis(4);
        for (int i=0; i< bd_basis->size(); i++)
            {
                gsMatrix <real_t> interval = bd_basis->function(i).support();
                real_t subint = (interval.at(1) - interval.at(0))/(p+1.0);
                gsMatrix <real_t> pres (4,1);
                pres << 0.0, 0.0, 0.0, 0.0;
                //pres = 0;
                for (int k=0; k<p+1;k++)
                    {
                        gsMatrix <real_t> pres0 (4,1);
                        pres0 << 0.0, 0.0, 0.0, 0.0;
                        gsMatrix <real_t> quad (2,14);
                        for (int l=0; l<14;l++)
                            {
                                quad.coeffRef(0,l) = ((interval.at(0)+subint*k)*(1.0-x.at(l))+
                                                      (interval.at(0)+subint*(k+1))*(1.0+x.at(l)))/2.0; //Maybe an error here
                                quad.coeffRef(1,l) = 1.0;
                            }
                        gsMatrix <real_t> u (4,14);
                        gsMatrix <real_t> u_val;
                        u_val = U1.value(quad);
                        u.row(0) = u_val.row(0);
                        u_val = U2.value(quad);
                        u.row(1) = u_val.row(0);
                        u_val = U3.value(quad);
                        u.row(2) = u_val.row(0);
                        u_val = U4.value(quad);
                        u.row(3) = u_val.row(0);

                        gsMatrix <real_t> normal (2,1);
                        normal = outer_normal(quad, mp, 4);
                        gsMatrix <real_t> funct = bd_basis->function(i).eval(quad.row(0));
                        gsMatrix <real_t> F_n (4,14);
                        quad = mapping(quad, mp);
                        real_t len_of = 0.0;
                        for (int l=0; l<14; l++)
                            {
                                //B1&2&3
                                F_n.col(l) = F_n_wall(normal.col(l),u.col(l));
                            }
                        for (int l=0; l<14; l++)
                            {
                                pres0 = pres0 + w.at(l)*funct.at(l)*F_n.col(l)*0.5;
                            }
                        for (int l=0; l<13; l++)
                            {
                                len_of = len_of + math::sqrt((quad.coeffRef(0,l+1)-quad.coeffRef(0,l))*(quad.coeffRef(0,l+1)-quad.coeffRef(0,l))
                                                             + (quad.coeffRef(1,l+1)-quad.coeffRef(1,l))*(quad.coeffRef(1,l+1)-quad.coeffRef(1,l)));
                            }
                        pres = pres + pres0 * len_of;
                    }
                for (int k=0; k<4; k++)
                    {
                        res.coeffRef(bases.totalSize()-bd_basis->size()+i+k*bases.totalSize()) += pres.at(k);
                    }
            }
    }
    return res;
}
//*/

// Initial condition
gsMatrix<real_t> ConstL2 (gsMatrix <real_t>       INrhs,
                          gsSparseMatrix <real_t> Mc,
                          gsSparseMatrix <real_t> Ml )
{
    gsSparseSolver<>::LU solver;
    gsMatrix<real_t> u_l (INrhs);
    gsMatrix<real_t> u_h (INrhs);
    solver.compute(Mc);
    u_h = solver.solve(INrhs);
    solver.compute(Ml);
    u_l = solver.solve(INrhs);
    gsMatrix<real_t> u_init (u_l);
    gsMatrix<real_t> Pp (INrhs);
    gsMatrix<real_t> Pm (INrhs);
    gsMatrix<real_t> Qp (INrhs);
    gsMatrix<real_t> Qm (INrhs);
    gsMatrix<real_t> ff (INrhs);
    for (int k = 0; k < Mc.outerSize(); k++)
        {
            Pp.coeffRef(k) = 0.0;
            Pm.coeffRef(k) = 0.0;
            Qp.coeffRef(k) = 0.0;
            Qm.coeffRef(k) = 0.0;
            ff.coeffRef(k) = 0.0;
        }
    gsSparseMatrix<real_t> f (Mc);
    for (int e=0; e<f.outerSize(); ++e)
        for (gsSparseMatrix<real_t>::InnerIterator it(f,e); it; ++it)
            {
                int k = it.row();
                int l = it.col();
                if ((u_h(k)-u_h(l))*(u_l(k)-u_l(l)) > 0)
                    f.coeffRef(k,l) = Mc.coeff(k,l)*(u_h.at(k)-u_h.at(l));
                else
                    f.coeffRef(k,l) = 0.0;
                if (l!= k)
                    {
                        Pp.coeffRef(k) += math::max<real_t>(0.0,f.coeff(k,l));
                        Pm.coeffRef(k) += math::min<real_t>(0.0,f.coeff(k,l));
                        Qp.coeffRef(k)  = math::max<real_t>(Qp.coeffRef(k),u_l(l)-u_l(k));
                        Qm.coeffRef(k)  = math::min<real_t>(Qm.coeffRef(k),u_l(l)-u_l(k));
                    }
            }
    gsMatrix<real_t> Rp (INrhs);
    gsMatrix<real_t> Rm (INrhs);
    for (int k = 0; k < f.outerSize(); k++)
        {
            if (Pp.at(k)==0)
                Rp.coeffRef(k) = 1.0;
            else
                Rp.coeffRef(k) = math::min<real_t>(1.0,Ml.coeff(k,k)*Qp.at(k)/(Pp.at(k)));
            if (Pm.at(k)==0)
                Rm.coeffRef(k) = 1.0;
            else
                Rm.coeffRef(k) = math::min<real_t>(1.0,Ml.coeff(k,k)*Qm.at(k)/(Pm.at(k)));
        }
    for (int e=0; e<f.outerSize(); ++e)
        for (gsSparseMatrix<real_t>::InnerIterator it(f,e); it; ++it)
            {
                int k = it.row();
                int l = it.col();
                if (f.coeff(k,l) > 0.0)
                    {
                        ff.coeffRef(k) += math::min(Rp.at(k),Rm.at(l))*f.coeff(k,l);
                    }
                else
                    {
                        ff.coeffRef(k) += math::min(Rm.at(k),Rp.at(l))*f.coeff(k,l);
                    }
            }
    for (int k=0; k<f.outerSize();k++)
        ff.coeffRef(k) = ff.at(k)/Ml.coeff(k,k);
    u_init = u_l + ff;
    return u_init;
}


gsMatrix<real_t> outer_normal(gsMatrix<real_t> points, gsMultiPatch<real_t> mp, boxSide a){
    gsMapData<real_t> data;
    data.addFlags(NEED_OUTER_NORMAL);
    data.side = a;//boundary::north;
    data.points = points;
    mp.patch(0).computeMap(data);
    data.outNormals.colwise().normalize();
    return data.outNormals;
}

gsMatrix<real_t> mapping(gsMatrix<real_t> points, gsMultiPatch<real_t> mp){
    gsMapData<real_t> data;
    data.addFlags(NEED_VALUE);
    data.points = points;
    mp.patch(0).computeMap(data);
    return data.values.at(0);
}
