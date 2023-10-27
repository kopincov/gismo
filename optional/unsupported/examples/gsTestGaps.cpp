
/**  gsTestGaps.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):     chofer
    Created on:  2014-12-03

    DESCRIBE THE PURPOSE OF THE EXAMPLE HERE

    WHAT IS HAPPENING

    EXPECTED RESULT
*/

#include <iostream>
#include <gismo.h>
#include <gismo_dev.h>

#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>


using namespace gismo;


//This function should make the geometry more IGA and less FEM ;)
void prepareGeometry(gsMultiPatch<> & patches, int testcase)
{
    if(testcase == 2||testcase ==0 || testcase == 9)
    {
        //FOR TWO SQUARES!!!
        gsMatrix<> & coef1 = patches.patch(0).coefs();
        for(int i = 0; i<coef1.rows();i++)
        {
            if(coef1(i,0)==-1 && coef1(i,1)==0)
                coef1(i,1) =-0.2;
            if(coef1(i,0)==-1 && coef1(i,1)==1)
                coef1(i,1) =1.2;
        }

        gsMatrix<> & coef0 = patches.patch(1).coefs();
        for(int i = 0; i<coef0.rows();i++)
        {
            if(coef0(i,0)==1 && coef0(i,1)==0)
                coef0(i,1) =0.2;
            if(coef0(i,0)==1 && coef0(i,1)==1)
                coef0(i,1) =0.8;
        }
    }

    else if(testcase == 5 || testcase == 6 ||testcase == 7||testcase == 8)
    {
        int dim = patches.parDim();
        for(size_t np = 0; np<patches.nPatches();np++)
        {
            if(dim == 2)
            {
                gsTensorBSpline<2>* p = dynamic_cast<gsTensorBSpline<2>* >(&patches.patch(np));
                p->insertKnot(0.5,0,2);
            }
            else if(dim==3)
            {
                gsTensorBSpline<3>* p = dynamic_cast<gsTensorBSpline<3>* >(&patches.patch(np));
                p->insertKnot(0.5,0,2);
                for(int d = 2; d<dim;d++)
                    p->insertKnot(0.5,d,1);
            }
            else if(dim ==4)
            {
                gsTensorBSpline<4>* p = dynamic_cast<gsTensorBSpline<4>* >(&patches.patch(np));
                p->insertKnot(0.5,0,2);
                for(int d = 2; d<dim;d++)
                    p->insertKnot(0.5,d,1);
            }
            else
                GISMO_ERROR("Dimension to high");
        }
    }

    if(testcase == 9)
    {
        short_t dim = patches.parDim();
        for(size_t np = 0; np<patches.nPatches();np++)
        {
            if(dim == 2)
            {
                gsTensorBSpline<2>* p = dynamic_cast<gsTensorBSpline<2>* >(&patches.patch(np));
                p->insertKnot(0.25,0,2);
                p->insertKnot(0.5,0,1);
                p->insertKnot(0.75,0,2);
            }
            else if(dim==3)
            {
                gsTensorBSpline<3>* p = dynamic_cast<gsTensorBSpline<3>* >(&patches.patch(np));
                p->insertKnot(0.25,0,2);
                p->insertKnot(0.5,0,1);
                p->insertKnot(0.75,0,2);
            }

        }

    }
    if(testcase == 8)
    {
        patches.patch(1).uniformRefine();
    }
}

//This function makes the gap with gapsize gap
void makeGapInGeometry(gsMultiPatch<> & patches, int  testcase, const gsMultiPatch<> & backup, real_t gap,int dim)
{
    if(testcase == 0)
    {
        int sgn;
        for(size_t np=0;np<2;np++)
        {
            np==0 ? sgn=-1 : sgn=1;
            gsMatrix<> & coef = patches.patch(np).coefs(); //get reference
            coef = backup.patch(np).coefs(); //override with backup (THIS IS A COPY!)
            for(int i = 0; i<coef.rows();i++)
            {
                if(coef(i,0)==0)
                {
                    bool interior = true;
                    for(int k=1; k<dim;k++)
                        if(!(math::abs(coef(i,k))>0.1 && math::abs(coef(i,k))<0.9))
                            interior=false;
                    if(interior)
                        coef(i,0) =sgn*gap;
                }
            }
        }

    }
    else if(testcase == 1|| testcase ==2)
    {
        gsMatrix<> & coef = patches.patch(1).coefs(); //get reference
        coef = backup.patch(1).coefs(); //override with backup (THIS IS A COPY!)
        for(int i = 0; i<coef.rows();i++)
        {
            if(coef(i,0)==0)
            {
                bool interior = true;
                for(int k=1; k<dim;k++)
                    if(!(math::abs(coef(i,k))>0.1 && math::abs(coef(i,k))<0.9))
                        interior=false;
                if(interior)
                    coef(i,0) =gap;
            }
        }
    }
    else if(testcase == 3)
    {
        gsMatrix<> & coef0 = patches.patch(0).coefs(); //get reference
        coef0 = backup.patch(0).coefs(); //override with backup (THIS IS A COPY!)
        int count =0;
        for(int i = 0; i<coef0.rows();i++)
            if(coef0(i,0)==0 && (coef0(i,1)==0))
            {
                if(count ==0)
                    coef0(i,0) =gap;
                else if(count ==1)
                { coef0(i,0) =gap;coef0(i,1) =gap;}
                else if(count ==2)
                    coef0(i,1) =gap;
                else if(count ==3)
                { coef0(i,0) =-gap;coef0(i,1) =gap;}
                else if(count ==4)
                    coef0(i,0) = -gap;
                count++;
            }

        gsMatrix<> & coef1 = patches.patch(1).coefs(); //get reference
        coef1 = backup.patch(1).coefs(); //override with backup (THIS IS A COPY!)
        count =0;
        for(int i = 0; i<coef1.rows();i++)
            if(coef1(i,0)==0 && (coef1(i,1)==0))
            {
                if(count ==0)
                    coef1(i,0) =gap;
                else if(count ==1)
                { coef1(i,0) =gap;coef1(i,1) =-gap;}
                else if(count ==2)
                    coef1(i,1) =-gap;
                else if(count ==3)
                { coef1(i,0) =-gap;coef1(i,1) =-gap;}
                else if(count ==4)
                    coef1(i,0) = -gap;
                count++;
            }
    }
    else if(testcase==4)
    {
        gsMatrix<> sgn(4,2);
        sgn<< 1,1,1,-1,-1,-1,-1,1;
        for(size_t np =0; np<patches.nPatches();np++)
        {
            gsMatrix<> & coef = patches.patch(np).coefs(); //get reference
            coef = backup.patch(np).coefs(); //override with backup (THIS IS A COPY!)
            for(index_t i = 0; i<coef.rows();i++)
                for(int k=0; k<dim;k++)
                    if((coef(i,k)==0 && (math::abs(coef(i,(k+1)%2))>0.1 && math::abs(coef(i,(k+1)%2))<0.9)))
                        coef(i,k) = sgn(np,k)*gap;

        }
    }
    else if(testcase == 5)
    {
        gsMatrix<> & coef = patches.patch(1).coefs(); //get reference
        coef = backup.patch(1).coefs(); //override with backup (THIS IS A COPY!)
        for(index_t i = 0; i<coef.rows();i++)
        {
            bool interior = true;
            for(int k=2; k<dim;k++)
                if(!(math::abs(coef(i,k))>0.01 && math::abs(coef(i,k))<5.99))
                    interior=false;
            if(interior && math::abs(coef(i,0)+coef(i,1))<1.e-6 && coef(i,0) !=0 && coef(i,0) !=-3.5)
            {
                coef(i,0)+= gap;
                // coef(i,1)+= gap;
            }
        }
    }
    else if(testcase == 6 || testcase  == 8)
    {
        gsMatrix<> sign(1,2) ;
        sign<<-1, 1;
        for(unsigned np=0;np<2;np++)
        {
            gsMatrix<> & coef = patches.patch(np).coefs(); //get reference
            coef = backup.patch(np).coefs(); //override with backup (THIS IS A COPY!)
            for(index_t i = 0; i<coef.rows();i++)
            {
                bool interior = true;
                for(int k=2; k<dim;k++)
                    if(!(math::abs(coef(i,k))>0.01 && math::abs(coef(i,k))<5.99))
                        interior=false;
                if(interior && math::abs(coef(i,0)+coef(i,1))<1.e-6 && coef(i,0) !=0 && coef(i,0) !=-3.5)
                {
                    coef(i,0)+= sign(0,np)*gap;
                    // coef(i,1)+= sign(0,np)*gap;
                }
            }
        }
    }
    else if(testcase == 7 )
    {
        gsMatrix<> sign(1,2) ;
        sign<<0, 1;
        for(unsigned np=0;np<2;np++)
        {
            gsMatrix<> & coef = patches.patch(np).coefs(); //get reference
            coef = backup.patch(np).coefs(); //override with backup (THIS IS A COPY!)
            for(int i = 0; i<coef.rows();i++)
            {
                bool interior = true;
                for(int k=2; k<dim;k++)
                    if(!(math::abs(coef(i,k))>0.01 && math::abs(coef(i,k))<5.99))
                        interior=false;

                if(interior && math::abs(coef(i,0)+coef(i,1))<1.e-6 && coef(i,0) !=0 && coef(i,0) !=-3.5)
                {
                    coef(i,0)+= 0.5+sign(0,np)*gap;
                    //coef(i,1)+= 0.5+sign(0,np)*gap; \\ no movement in y-direction!
                }
            }
        }
    }
    else if(testcase == 9)
    {
        gsMatrix<> sign(1,2) ;
        sign<<1, -1;
        for(size_t np=0;np<patches.nPatches();np++)
        {
            gsMatrix<> & coef = patches.patch(np).coefs(); //get reference
            coef = backup.patch(np).coefs(); //override with backup (THIS IS A COPY!)
            for(int i = 0; i<coef.rows();i++)
            {
                if(math::abs(coef(i,0))< 1.e-4)
                {
                    bool interior = true;
                    for(int k=1; k<dim;k++)
                        if(!(math::abs(coef(i,k))>0.01 && math::abs(coef(i,k))<0.99))
                            interior=false;
                    if(interior && coef(i,1)<0.49 )//gap
                        coef(i,0) = sign(0,np)*gap;
                    else if(interior && coef(i,1)>0.51 ) //overlap
                        coef(i,0) =-sign(0,np)*gap;
                }
            }
        }

    }
}

real_t getInitialMeshsize(int testcase)
{
    switch(testcase){
    case 0:
    case 1:
    case 2:
    case 9:
        return 0.5;
    case 3:
        return 2;
    case 4:
        return 1;
    case 5:
    case 6:
    case 7:
    case 8:
        return 6;
    }
    return 1;
}

int main (int argc, char** args)
{

#ifdef _OPENMP
    int t;
    Eigen::initParallel();
    const char * nProcs = getenv("OMP_NUM_THREADS");
    if(nProcs != NULL)
        sscanf( nProcs, "%d", &t);
    else
        t = 1; //one OpenMP thread
    Eigen::setNbThreads(t);
#endif

    bool test3D = false; //maxRefine =5;
    bool test4D = false; //maxRefine =4;
    bool plot = false;
    bool overlap = false;

    index_t testcase = 7;
    //initRef + maxRef should not be larger than 8
    index_t initRef = 0;
    index_t maxRefine = 2;
    index_t numElevate = 0;
    index_t iniScal = 2;
    index_t maxNGaps = 1;

    index_t solCase = 14;

    int result=0;
    int oSgn = 1;

    gsCmdLine cmd("Hi, I will test the convergence behaviour of IETI");
    cmd.addSwitch("d", "dim", "test3D", test3D);
    cmd.addSwitch("p", "plot", "Plot solution", plot);
    cmd.addSwitch("x", "dimfour", "test4D", test4D);
    cmd.addSwitch("o", "overlap", "overlap", overlap);
    cmd.addInt("r","refine","Number of refinements",maxRefine);
    cmd.addInt("c","case","Choosen testcase",testcase);
    cmd.addInt("e","elevate","maximal elevation",numElevate);
    cmd.addInt("i","initRefine","initial refine",initRef);
    cmd.addInt("g","gaps","max number of gaps",maxNGaps);
    cmd.addInt("s","solcase", "Choosen solutioncase",solCase);
    cmd.addInt("a", "initGap", "init scaling for gap ", iniScal);

    try { cmd.getValues(argc,args); } catch (int rv) { return rv; }


    iFace::strategy fStrat = iFace::dg;

    plot = test4D ? false : plot;
    test3D = test4D ? true : test3D;

    //gsVector<real_t> lambdas(3);
    //lambdas << 1,2,3;

    gsVector<real_t> lambdas(4+numElevate);
    for(int i=0; i<3+numElevate-1;i++)
        lambdas[i]=i+1;
    lambdas[3+numElevate-1]=3+numElevate-1+0.5;
    lambdas[3+numElevate]=3+numElevate;

    gsVector<real_t> gaps(8);
    gaps << 0, 0.1, 0.01,0.005, 0.001,0.0001,0.00005, 0.00001;


    real_t rho=1;
    int dim;
    test4D ? dim =4 : (test3D ? dim = 3 : dim=2);

    if(overlap)
        oSgn = -1;
    //////// Right-hand side and analytical solution ////////
    gsFunctionExpr<>* f,*g, * g_temp = 0;
    gsPiecewiseFunction<real_t>  alpha;
    switch(solCase)
    {
    default:
    case 0:
    {
        //Gap Paper for SIAM
        f = new gsFunctionExpr<>("((4*pi)^2 + (5*pi)^2)*sin(pi*x*5)*sin(pi*y*4)",dim);
        g = new gsFunctionExpr<>("sin(pi*x*5)*sin(pi*y*4)",dim);

        // Gap paper with 4D example
        /*
        if(dim==4)
        {
            f = new gsFunctionExpr<>("((pi/2)^2+(pi/8)^2+(pi/2)^2)*sin(pi*y/2)*sin(pi*x/8)*sin(pi*z/2)",dim);
            g = new gsFunctionExpr<>("sin(pi*y/2)*sin(pi*x/8)*sin(pi*z/2)",dim);
        }
        else
        {
            f = new gsFunctionExpr<>("((pi/2)^2+(pi/8)^2)*sin(pi*y/2)*sin(pi*x/8)",dim);
            g = new gsFunctionExpr<>("sin(pi*y/2)*sin(pi*x/8)",dim);
        }
*/
        gsFunctionExpr<> a("1",dim);
        alpha.addPiece(a);
        alpha.addPiece(a);
        break;
    }
    case 1:
    {
        f = new gsFunctionExpr<>("((pi*1)^2 + (pi*2)^2)*sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)",dim);
        g = new gsFunctionExpr<>("sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)+x+y",dim);
        // g->set_x_der("pi*cos(pi*(x+0.4))*sin(pi*(y+0.3)*2)+1",dim);
        // g->set_y_der("2*pi*sin(pi*(x+0.4))*cos(pi*(y+0.3)*2)+1",dim);

        gsFunctionExpr<> a("1",dim);
        alpha.addPiece(a);
        alpha.addPiece(a);
        break;
    }
    case 2:
    {
        GISMO_ASSERT(testcase!=4, "This rhs is only made for 2 patches, choose another one.");

        rho = 4* EIGEN_PI;

        std::ostringstream strs;
        strs << rho;
        std::string rho_str = strs.str();

        f = new gsFunctionExpr<>("if(x<0,-"+rho_str+"*exp(x),(4*pi)^2*sin(4*pi*x))",dim);
        g_temp = new gsFunctionExpr<>("if(x<0,"+rho_str+"*(exp(x)-1),sin(4*pi*x))",dim);
        g= new gsFunctionExpr<>("if(x<0,(exp(x)-1),sin(4*pi*x))",dim);
        // g->set_x_der("if(x<0,exp(x),4*pi*cos(4*pi*x))",dim);
        // g->set_y_der("0");

        gsFunctionExpr<> a("1",dim);
        alpha.addPiece(a);
        alpha.addPiece(a);
        break;
    }
    case 3:
    {
        GISMO_ASSERT(testcase!=4, "This rhs is only made for 2 patches, choose another one.");
        f = new gsFunctionExpr<>("if(x<0,4,-2)",dim);
        g = new gsFunctionExpr<>("if(x<0,2*(x-x^2),x^2+2*x)",dim);
        // g->set_x_der("if(x<0,2*(1-2*x),2*x+2)",dim);
        // g->set_y_der("0",dim);

        gsFunctionExpr<> a("1",dim);
        alpha.addPiece(a);
        alpha.addPiece(a);
        break;
    }
    case 4:
    {
        GISMO_ASSERT(testcase!=4, "This rhs is only made for 2 patches, choose another one.");
        f = new gsFunctionExpr<>("if(x<0,sin(x),20)",dim);
        g = new gsFunctionExpr<>("if(x<0,sin(x),10*x*(0.1-x))",dim);
        // g->set_x_der("if(x<0,cos(x),-20*(x-0.05))");
        // g->set_y_der("0");

        gsFunctionExpr<> a("1",dim);
        alpha.addPiece(a);
        alpha.addPiece(a);
        break;
    }
    case 5:
    {
        GISMO_ASSERT(testcase!=4, "This rhs is only made for 2 patches, choose another one.");

        rho = 4* EIGEN_PI;

        std::ostringstream strs;
        strs << rho;
        std::string rho_str = strs.str();

        f = new gsFunctionExpr<>("if(x<0,-"+rho_str+"*exp(x),(4*pi)^2*sin(4*pi*x))",dim);
        g = new gsFunctionExpr<>("if(x<0,(exp(x)-1),sin(4*pi*x))",dim);
        // g->set_x_der("if(x<0,exp(x),4*pi*cos(4*pi*x))",dim);
        // g->set_y_der("0",dim);

        gsFunctionExpr<> a1(rho_str,dim);
        gsFunctionExpr<> a2("1",dim);
        alpha.addPiece(a1);
        alpha.addPiece(a2);
        break;
    }
    case 6:
    {
        //Edge Singularity
        GISMO_ASSERT(testcase!=4, "This rhs is only made for 2 patches, choose another one.");
        real_t param_a=1.5;
        real_t param_reg =0.000001;
        std::stringstream strs;

        strs << param_a;
        std::string a_str = strs.str();
        strs.str(std::string());

        strs << std::fixed<<param_reg;
        std::string reg_str = strs.str();

        //f = new gsFunctionExpr<>("u:=x;v:=x;w:="+a_str+";if(x<0,-w*(w-1)*pow(v^2+"+reg_str+",(w/2.)-1),-pow(u^2+"+reg_str+",(w/2.)-1)*((w^2-w-u^2)*cos(u)-2*w*u*sin(u)))",dim);
        f = new gsFunctionExpr<>("w:="+a_str+";v:="+reg_str+";if(x<0,-w*(v+(w-1)*x^2)*pow(x^2+v,(w/2.)-2),-pow(x^2+v,(w/2.)-2)*(-(x^4+x^2*(2*v+w-w^2)+v*(v-w))*cos(x)-2*w*x*(v+x^2)*sin(x)))",dim);
        //g = new gsFunctionExpr<>("u:=x;v:=x;w:="+a_str+";if(x<0,pow(v^2+"+reg_str+",w/2),pow(u^2+"+reg_str+",w/2.)*cos(u))",dim);
        g = new gsFunctionExpr<>("w:="+a_str+";v:="+reg_str+";if(x<0,pow(x^2+v,w/2),pow(x^2+v,w/2.)*cos(x))",dim);
        // g->set_x_der("w:="+a_str+";v:="+reg_str+";if(x<0, w*x*pow(x^2+v,(w/2.)-1),w*x*pow(x^2+v,(w/2.)-1)*cos(x)-pow(x^2+v,w/2.)*sin(x))");
        // g->set_y_der("0");

        gsFunctionExpr<> a("1",dim);
        alpha.addPiece(a);
        alpha.addPiece(a);
        break;
    }
    case 7:
    {
        //Point Singularity!
        GISMO_ASSERT(testcase!=4, "This rhs is only made for 2 patches, choose another one.");
        real_t param_a=1;
        real_t param_reg =0;

        std::ostringstream strs;
        strs << param_a;
        std::string a_str = strs.str();
        strs.str(std::string());

        strs << param_reg;
        std::string reg_str = strs.str();
        f = new gsFunctionExpr<>("v:="+a_str+"; w:=pow(x-0.,2.)+pow(y-0.5,2.)+"+reg_str+"; -( v/2.*(v/2.-1.)*pow( w,v/2.-2.)*4.*(w-0.000)   +  2.*v*(pow(w,v/2.-1.)))",dim);
        g = new gsFunctionExpr<>("v:="+a_str+"; w:=pow(x-0.,2.)+pow(y-0.5,2.)+"+reg_str+"; pow( w,v/2.)",dim);
        // g->set_x_der("v:="+a_str+"; w:=pow(x-0.,2.)+pow(y-0.5,2.)+"+reg_str+"; v*pow( w,v/2.-1)*x",dim);
        // g->set_y_der("v:="+a_str+"; w:=pow(x-0.,2.)+pow(y-0.5,2.)+"+reg_str+"; v*pow( w,v/2.-1)*(y-0.5)",dim);

        gsFunctionExpr<> a("1",dim);
        alpha.addPiece(a);
        alpha.addPiece(a);
        break;
    }
    case 8:
    {
        //Simplier Edge Singularity
        GISMO_ASSERT(testcase!=4, "This rhs is only made for 2 patches, choose another one.");
        real_t param_a=2;
        real_t param_reg =1.e-6;
        std::stringstream strs;

        strs << param_a;
        std::string a_str = strs.str();
        strs.str(std::string());

        strs << param_reg;
        std::string reg_str = strs.str();

        f = new gsFunctionExpr<>("u:=x+"+reg_str+";v:=x-"+reg_str+";w:="+a_str+";if(x<0, -w*(w-1)*pow(v^2,(w/2)-1),-w*(w-1)*pow(v^2,(w/2)-1))",dim);
        g = new gsFunctionExpr<>("u:=x+"+reg_str+";v:=x-"+reg_str+";w:="+a_str+";if(x<0,pow(v^2,w/2),pow(u^2,w/2))",dim);
        // g->set_x_der("u:=x+"+reg_str+";v:=x-"+reg_str+";w:="+a_str+";if(x<0, w*v*pow(v^2,(w/2)-1),w*v*pow(v^2,(w/2)-1))",dim);
        // g->set_y_der("0");

        gsFunctionExpr<> a("1",dim);
        alpha.addPiece(a);
        alpha.addPiece(a);
        break;
    }
    case 9:
    {
        //Point Singularity 2!
        GISMO_ASSERT(testcase!=4, "This rhs is only made for 2 patches, choose another one.");
        real_t param_a=3;
        real_t param_reg =0;

        std::ostringstream strs;
        strs << param_a;
        std::string a_str = strs.str();
        strs.str(std::string());

        strs << param_reg;
        std::string reg_str = strs.str();
        f = new gsFunctionExpr<>("v:="+a_str+"; w:=pow(x-0.,2.)+pow(y-0.5,2.); pow(w,v/2. -1)*((-v^2+100*w)*cos(10*x)+20*v*x*sin(10*x))",dim);
        g = new gsFunctionExpr<>("v:="+a_str+"; w:=pow(x-0.,2.)+pow(y-0.5,2.); pow(w,v/2.)*cos(10*x)",dim);
        // g->set_x_der("v:="+a_str+"; w:=pow(x-0.,2.)+pow(y-0.5,2.); v*x*pow(w,(v/2.)-1)*cos(10*x)-10*pow(w,v/.2)*sin(10*x)",dim);
        // g->set_y_der("v:="+a_str+"; w:=pow(x-0.,2.)+pow(y-0.5,2.); v*pow( w,v/2.-1)*(y-0.5)*cos(10*x)",dim);

        gsFunctionExpr<> a("1",dim);
        alpha.addPiece(a);
        alpha.addPiece(a);
        break;
    }
    case 10:
    {
        rho = 4;

        std::ostringstream strs;
        strs << rho;
        std::string rho_str = strs.str();

        f = new gsFunctionExpr<>("v:="+rho_str+";if(x<0,2*pi^2*v*sin(pi*(x+y)),pi^2*(1+v^2)*sin(pi*(v*x+y)))",dim);
        g = new gsFunctionExpr<>("v:="+rho_str+"; if(x<0,sin(pi*(x+y)),sin(pi*(v*x+y)))",dim);
        // g->set_x_der("v:="+rho_str+";if(x<0,pi*cos(pi*(x+y)),pi*v*cos(pi*(v*x+y)))",dim);
        // g->set_y_der("v:="+rho_str+";if(x<0,pi*cos(pi*(x+y)),pi*cos(pi*(v*x+y)))",dim);

        gsFunctionExpr<> a1(rho_str,dim);
        gsFunctionExpr<> a2("1",dim);
        alpha.addPiece(a1);
        alpha.addPiece(a2);
        break;
    }
    case 11:
    {
        f = new gsFunctionExpr<>("v:=x+y;if(v<0,-2*3*3*exp(3*v),2*3*sin(v))",dim);
        g = new gsFunctionExpr<>("v:=x+y; if(v<0,(exp(3*v)-1),sin(v))",dim);
        gsFunctionExpr<> a1("1",dim);
        gsFunctionExpr<> a2("3",dim);
        alpha.addPiece(a1);
        alpha.addPiece(a2);
        break;
    }
    case 12:
    {
        f = new gsFunctionExpr<>("v:=x+y;if(v<0,2*2*exp(sin(v))*(sin(v)-cos(v)*cos(v)),2*4*sin(2*v))",dim);
        g = new gsFunctionExpr<>("v:=x+y; if(v<0,(exp(sin(v))-1),sin(2*v))",dim);
        gsFunctionExpr<> a1("2",dim);
        gsFunctionExpr<> a2("1",dim);
        alpha.addPiece(a1);
        alpha.addPiece(a2);
        break;
    }
    case 13:
    {
        f = new gsFunctionExpr<>("v:=x+y;if(v<0,2*pi*exp(sin(v))*(sin(v)-cos(v)*cos(v)),2*pi*pi*sin(pi*v))",dim);
        g = new gsFunctionExpr<>("v:=x+y; if(v<0,(exp(sin(v))-1),sin(pi*v))",dim);
        gsFunctionExpr<> a1("pi",dim);
        gsFunctionExpr<> a2("1",dim);
        alpha.addPiece(a1);
        alpha.addPiece(a2);
        break;
    }
    case 14:
    {
        f = new gsFunctionExpr<>("((pi/2)^2+(pi/8)^2)*sin(pi*y/2)*sin(pi*x/8)",dim);
        g = new gsFunctionExpr<>("sin(pi*y/2)*sin(pi*x/8)",dim);

        gsFunctionExpr<> a("1",dim);
        alpha.addPiece(a);
        alpha.addPiece(a);
        break;
    }
    case 15: //does not give right results
    {
        rho = 4* EIGEN_PI;

        std::ostringstream strs;
        strs << rho;
        std::string rho_str = strs.str();

        f = new gsFunctionExpr<>("if(x<0,-"+rho_str+"*exp(x),(4*pi)^2*sin(4*pi*x))",dim);
        g = new gsFunctionExpr<>("if(x<0,(exp(x)-1),sin(4*pi*x))",dim);
        gsFunctionExpr<> a1("if(x<=0,"+rho_str+",1)",dim);
        gsFunctionExpr<> a2("if(x<0,"+rho_str+",1)",dim);
        alpha.addPiece(a1);
        alpha.addPiece(a2);
        break;
    }
    }
    // Print out source function and solution
    gsInfo<<"Source function "<< *f << "\n";
    gsInfo<<"Exact solution "<< *g <<"\n" << "\n";

    /////////////////// Setup geometry ///////////////////
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;

    switch(testcase)
    {
    default:
    case 0: //2 squares
    case 1:
    case 2:
    {
        gsFileData<> fileData("gaps/2pSquare_lin_0.xml");
        if (!fileData.getFirst< gsMultiPatch<> >(patches))
            return 1;
        break;
    }
    case 3:
    {
        gsFileData<> fileData("gaps/Circle2p_0.xml");
        if (fileData.getFirst< gsMultiPatch<> >(patches))
        {
            //gsInfo<<"Degree Patch: "<<patches.basis(0).degree(0)<<" , "<<patches.basis(0).degree(1)<<"\n";
        }
        else
            return 1;
        break;
    }
    case 4:
    {
        gsFileData<> fileData("gaps/4pSquare_quadr_0.xml");
        if (!fileData.getFirst< gsMultiPatch<> >(patches))
            return 1;
        break;
    }
    case 5: //bumper
    case 6:
    case 7:
    case 8:
    {
        gsFileData<> fileData("planar/bumper.xml");
        if (!fileData.getFirst< gsMultiPatch<> >(patches))
            return 1;
        patches.computeTopology();
        break;
    }
    case 9:
    {
        gsFileData<> fileData("gaps/2pSquare_lin_0.xml");
        if (!fileData.getFirst< gsMultiPatch<> >(patches))
            return 1;
        break;
    }
    }

    if(test3D)
    {
        gsMultiPatch<real_t> patches3D;
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<2,real_t>* tb = static_cast<gsTensorBSpline<2,real_t>* >(&patches.patch(i));
            if(tb==NULL)
                GISMO_ERROR("not a 2D-domain");
            if(testcase>=5 && testcase <=8)
                        patches3D.addPatch(gsNurbsCreator<>::lift3D(*tb,6));
            else
                patches3D.addPatch(gsNurbsCreator<>::lift3D(*tb));
        }
        patches3D.computeTopology();

        patches = give(patches3D);
    }

    if(test4D)
    {
        gsMultiPatch<real_t> patches4D;
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<3,real_t>* tb = static_cast<gsTensorBSpline<3,real_t>* >(&patches.patch(i));
            if(tb==NULL)
                GISMO_ERROR("not a 3D-domain");
            if(testcase>=5 && testcase <=8)
                        patches4D.addPatch(gsNurbsCreator<>::lift4D(*tb,6));
            else
                patches4D.addPatch(gsNurbsCreator<>::lift4D(*tb));
        }
        patches4D.computeTopology();

        patches = give(patches4D);
    }

    if((testcase >=0 && testcase<=2) || testcase == 9)
    {
        patches.degreeElevate(1);
        patches.uniformRefine();
    }
    gsInfo<<patches.detail()<<std::endl;
    //Make some geometry more IgA
    prepareGeometry(patches, testcase);

    gsInfo<<"p0:\n"<< patches.patch(0).coefs()<<"\n";
    gsInfo<<"p1:\n"<< patches.patch(1).coefs()<<"\n";

    real_t h0 = getInitialMeshsize(testcase);
    //Backup for introducing gaps
    gsMultiPatch<> backup(patches);


    if(plot)
    {
        gsField<> gField(patches, *g,false);
        gsWriteParaview<>( gField, "gaps_ex_conf", 10000, true);
        result = system("paraview gaps_ex_conf.pvd &");
    }


    /////////////////// Setup boundary conditions ///////////////////
    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;

    for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
        if(solCase==2)
            bcInfo.addCondition(*it, condition_type::dirichlet, g_temp);
        else
            bcInfo.addCondition(*it, condition_type::dirichlet, g);


    gsInfo<<bcInfo<<"\n";
    ///////////////////////ASSEMBLE////////////////////////////////////

    gsMatrix<real_t> l2Error(maxRefine,maxNGaps);
    gsMatrix<real_t> h1error(maxRefine,maxNGaps);
    gsMatrix<real_t> dGerror(maxRefine,maxNGaps);
    gsMatrix<real_t> h(maxRefine,maxNGaps);

    gsMatrix<real_t> l2Factor(maxRefine,maxNGaps);
    gsMatrix<real_t> h1Factor(maxRefine,maxNGaps);
    gsMatrix<real_t> dgFactor(maxRefine,maxNGaps);


    gsField<> Sol, Sol_temp;
    gsMatrix<> sol;



    for(int gap = 0; gap<maxNGaps;gap++)
    {
        // Make the gap
        makeGapInGeometry(patches, testcase,backup,oSgn*gaps[gap],dim);

        //Exact solution
        gsField<> gField(patches, *g,false);

        //Make the basis for refinement
        //gsMultiBasis<> refine_bases(refine_bases_backup);

        gsMultiBasis<> refine_bases(patches);

        /////////////////// Basis setup ///////////////////
        gsMultiBasis<> refine_bases_backup;
        if(testcase == 3)
        {
            for(size_t np =0; np<patches.nPatches();np++)
            {
                gsTensorNurbsBasis<2> & NurbsBasis = dynamic_cast<gsTensorNurbsBasis<2> &>(patches.patch(np).basis());
                gsBasis<>::uPtr bptr = memory::convert_ptr<gsBasis<> >(NurbsBasis.source().clone());
                refine_bases_backup.addBasis(bptr.release());
            }
            GISMO_ERROR("DOES NOT WORK!, NO NURBS ALLOWED FOR THE MOMENT");
        }

        for(int i=0; i<initRef;i++)
            refine_bases.uniformRefine();

        // Elevate and p-refine the basis to order k + numElevate
        // where k is the highest degree in the bases
        if ( numElevate > -1 )
        {
            // Find maximum degree with respect to all the variables
            int max_tmp = refine_bases.maxDegree(0);
            for (int j = 1; j < patches.parDim(); ++j )
                if ( max_tmp < refine_bases.maxDegree(j) )
                    max_tmp = refine_bases.maxDegree(j);

            // Elevate all degrees uniformly
            // max_tmp += numElevate;

            refine_bases.setDegree(max_tmp);
            refine_bases.degreeElevate(numElevate);
        }
        gsInfo<<"Degree set to: "<<refine_bases.basis(0).degree(0)<<" , "<<refine_bases.basis(0).degree(1)<<"\n";

        //Start the loop over the refinements
        for(int ref = 0; ref<maxRefine;ref++)
        {
            gsInfo<<"Starting ref: "<<ref<<"\n";

            //gsPoissonAssembler2<real_t> assembler(patches,refine_bases,bcInfo,*f, dirichlet::nitsche,iFace::dg);
            //gsPoissonPde<real_t>ppde (*patches,bcInfo,*f);
            //gsPoissonAssembler<real_t> assembler(ppde,refine_bases,dirichlet::nitsche,fStrat);
            gsPoissonHeterogeneousPde<real_t>ppde(patches,bcInfo,*f,alpha);
            gsPoissonHeterogeneousAssembler<real_t>assembler (ppde,refine_bases,dirichlet::elimination,fStrat);

            gsInfo<<"dofs "<<assembler.numDofs()<<"\n";

            gsInfo<<"assembling\n"<<std::flush;
            assembler.assemble();
            sol.setZero(assembler.numDofs(),1);

#if defined(GISMO_WITH_PARDISO)
            gsSparseSolver<>::PardisoLU solver;
#elif defined(GISMO_WITH_SUPERLU)
            gsSparseSolver<>::SuperLU solver;
#else
            gsSparseSolver<>::LU solver;
#endif
            gsInfo<<"solving\n"<<std::flush;
            solver.compute(assembler.matrix());
            sol = solver.solve( assembler.rhs());
            if(solCase==2)
            {
                Sol_temp = assembler.constructSolution(sol);
                int p=0;

                gsDofMapper mapper = assembler.system().dofMappers()[0];
                const int sz  = refine_bases[p].size();
                for (index_t i = 0; i < sz; ++i)
                    if ( mapper.is_free(i, p) )
                        sol.row( mapper.index(i, p) )/=rho;

            }

            Sol = assembler.constructSolution(sol);

            l2Error(ref,gap) = Sol.distanceL2(*g);
            h1error(ref,gap) = Sol.distanceH1(*g);
            if(solCase==2)
                dGerror(ref,gap) = Sol_temp.distanceDG(*g_temp);
            else
                dGerror(ref,gap) = Sol.distanceDG(*g);

            h(ref,gap) = h0*1. / math::exp2(ref+initRef);

            refine_bases.uniformRefine();

            if (plot && ref == 2)
            {
                // Write approximate and exact solution to paraview files
                gsInfo<<"Plotting in Paraview...\n";
                gsWriteParaview<>(Sol, "gaps", 8000,true);
                gsWriteParaview<>(gField, "gaps_ex", 8000);

                // Run paraview
                result = system("paraview gaps.pvd &");
            }
        }


        l2Factor.row(0).setZero();
        h1Factor.row(0).setZero();
        dgFactor.row(0).setZero();
        for(int i=1;i<maxRefine;i++)
        {
            l2Factor(i,gap) = math::log(l2Error(i-1,gap)/l2Error(i,gap))/math::log(2.0);
            h1Factor(i,gap) = math::log(h1error(i-1,gap)/h1error(i,gap))/math::log(2.0);
            dgFactor(i,gap) = math::log(dGerror(i-1,gap)/dGerror(i,gap))/math::log(2.0);
        }
    }

    gsInfo<<"Gaps: "<<gaps.transpose()<<"\n"<<"\n";
    gsInfo<<"L2Error: "<<"\n"<<l2Error<<"\n"<<"\n";
    gsInfo<<"L2Factor: "<<"\n"<<l2Factor<<"\n"<<"\n";

    gsInfo<<"H1Error: "<<"\n"<<h1error<<"\n"<<"\n";
    gsInfo<<"H1Factor: "<<"\n"<<h1Factor<<"\n"<<"\n";
    gsInfo<<"DgError: "<<"\n"<<dGerror<<"\n"<<"\n";
    gsInfo<<"DGFactor: "<<"\n"<<dgFactor<<"\n"<<"\n";
    gsInfo<<"h: "<<"\n"<<h<<"\n"<<"\n";


    ///////////////////////ASSEMBLE////////////////////////////////////

    l2Error.setZero(maxRefine,lambdas.rows());
    h1error.setZero(maxRefine,lambdas.rows());
    dGerror.setZero(maxRefine,lambdas.rows());
    h.setZero(maxRefine,lambdas.rows());

    l2Factor.setZero(maxRefine,lambdas.rows());
    h1Factor.setZero(maxRefine,lambdas.rows());
    dgFactor.setZero(maxRefine,lambdas.rows());

    gsMatrix<> Gaps(maxRefine,lambdas.rows());
    maxNGaps = maxRefine;

    //gap = h/2


    for(int l = 0; l<lambdas.rows();l++)
    {
        real_t lambda = lambdas(l);
        gsMultiBasis<> refine_bases(patches);



        // Elevate and p-refine the basis to order k + numElevate
        // where k is the highest degree in the bases
        if ( numElevate > -1 )
        {
            // Find maximum degree with respect to all the variables
            int max_tmp = refine_bases.maxDegree(0);
            for (int j = 1; j < patches.parDim(); ++j )
                if ( max_tmp < refine_bases.maxDegree(j) )
                    max_tmp = refine_bases.maxDegree(j);

            // Elevate all degrees uniformly
            //max_tmp += numElevate;

            refine_bases.setDegree(max_tmp);
            refine_bases.degreeElevate(numElevate);
        }


        for(int i=0; i<initRef;i++)
            refine_bases.uniformRefine();

        gsInfo<<"Degree set to: "<<refine_bases.basis(0).degree(0)<<" , "<<refine_bases.basis(0).degree(1)<<"\n";

        for(int gap = 0; gap<maxNGaps;gap++)
        {
            Gaps(gap,l)= oSgn * h0 / math::pow((real_t)2.,lambda*(gap+initRef)+iniScal);
            // Make the gap
            makeGapInGeometry(patches, testcase,backup,Gaps(gap,l),dim);
            // gsInfo<<"p1:\n"<< patches.patch(1).coefs()<<"\n";
            //Exact solution
            gsField<> gField(patches, *g,false);

            gsInfo<<"Starting ref: "<<gap<<"\n";

            //gsPoissonAssembler2<real_t> assembler(patches,refine_bases,bcInfo,*f, dirichlet::nitsche,iFace::dg);
            //gsPoissonHeterogeneousAssembler<real_t> assembler(*patches,refine_bases,bcInfo,*f,alpha, dirichlet::elimination,fStrat);

            //gsPoissonPde<real_t>ppde (*patches,bcInfo,*f);
            //gsPoissonAssembler<real_t> assembler(ppde,refine_bases,dirichlet::nitsche,fStrat);

            gsPoissonHeterogeneousPde<real_t>ppde(patches,bcInfo,*f,alpha);
            gsPoissonHeterogeneousAssembler<real_t>assembler (ppde,refine_bases,dirichlet::elimination,fStrat);

            gsInfo<<"dofs "<<assembler.numDofs()<<"\n";
            gsInfo<<"assembling\n";
            assembler.assemble();

            sol.setZero(assembler.numDofs(),1);

#if defined(GISMO_WITH_PARDISO)
            gsSparseSolver<>::PardisoLU solver;
#elif defined(GISMO_WITH_SUPERLU)
            gsSparseSolver<>::SuperLU solver;
#else
            gsSparseSolver<>::LU solver;
#endif
            gsInfo<<"solving\n";
            solver.compute(assembler.matrix());
            sol = solver.solve( assembler.rhs());
            gsInfo<<"done\n";
            if(solCase==2)
            {
                Sol_temp = assembler.constructSolution(sol);
                int p=0;
                gsDofMapper mapper = assembler.system().dofMappers()[0];
                const int sz  = refine_bases[p].size();
                for (index_t i = 0; i < sz; ++i)
                    if ( mapper.is_free(i, p) )
                        sol.row( mapper.index(i, p) )/=rho;

            }

            Sol = assembler.constructSolution(sol);

            l2Error(gap,l) = Sol.distanceL2(*g);
            h1error(gap,l) = Sol.distanceH1(*g);
            if(solCase==2)
                dGerror(gap,l) = Sol_temp.distanceDG(*g_temp);
            else
                dGerror(gap,l) = Sol.distanceDG(*g);

            h(gap,l) = h0*1. / math::exp2(gap+initRef);

            refine_bases.uniformRefine();

            if (plot && gap == 2 && l==0)
            {

                // Write approximate and exact solution to paraview files
                gsInfo<<"Plotting in Paraview...\n";
                gsWriteParaview<>( Sol, "gaps", 8000, true);
                gsWriteParaview<>( gField, "gaps_ex", 8000);

                // Run paraview
                result = system("paraview gaps.pvd &");

            }
        }
    }
    l2Factor.row(0).setZero();
    h1Factor.row(0).setZero();
    dgFactor.row(0).setZero();
    for(int i=1;i<maxRefine;i++)
    {
        for(int l= 0; l<lambdas.rows();l++)
        {
            l2Factor(i,l) = math::log(l2Error(i-1,l)/l2Error(i,l))/math::log(2.0);
            h1Factor(i,l) = math::log(h1error(i-1,l)/h1error(i,l))/math::log(2.0);
            dgFactor(i,l) = math::log(dGerror(i-1,l)/dGerror(i,l))/math::log(2.0);
        }
    }

    gsInfo<<"h^lambda: "<<lambdas.transpose()<<"\n"<<"\n";
    gsInfo<<"Gaps: \n"<<Gaps.transpose()<<"\n"<<"\n";
    gsInfo<<"L2Error: "<<"\n"<<l2Error<<"\n"<<"\n";
    gsInfo<<"L2Factor: "<<"\n"<<l2Factor<<"\n"<<"\n";

    gsInfo<<"H1Error: "<<"\n"<<h1error<<"\n"<<"\n";
    gsInfo<<"H1Factor: "<<"\n"<<h1Factor<<"\n"<<"\n";
    gsInfo<<"DgError: "<<"\n"<<dGerror<<"\n"<<"\n";
    gsInfo<<"DGFactor: "<<"\n"<<dgFactor<<"\n"<<"\n";

    gsInfo<<"DGError for output: "<<"\n";
    for(int i=0; i< dGerror.rows();i++)
        gsInfo<<h(i,0)<<"   "<<dGerror.row(i)<<"\n";

    gsInfo<<"DGFactor for output: "<<"\n";
    for(int i=0; i< dgFactor.rows();i++)
        gsInfo<<h(i,0)<<"   "<<dgFactor.row(i)<<"\n";
    gsInfo<<"h: "<<"\n"<<h<<"\n"<<"\n";




    //return  1 for failures and 0 for success

    delete f;
    delete g;
    if(solCase==2)
        delete g_temp;

    return result;
}


