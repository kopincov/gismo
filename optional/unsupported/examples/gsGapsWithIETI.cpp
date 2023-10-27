

/**  gsGapsWithIETI.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):      C. Hofer
    Created on:  2015-01-20

    DESCRIBE THE PURPOSE OF THE EXAMPLE HERE

    WHAT IS HAPPENING

    EXPECTED RESULT
*/

#include <iomanip>
#include <sstream>
#include <iostream>
#include <gismo.h>
#include <gismo_dev.h>
#include <gsIETI/gsIETIUtils.h>
#include <gsIETI/gsIETIdGAssembler.h>
#include <gsIETI/gsIETISolver.h>
#include <gsIETI/gsIETIScaledDirichlet.h>
#include <gsSolver/gsConjugateGradient.h>
#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>

using namespace gismo;


/// A locally used geometry
gsMultiPatch<real_t> approximateQuarterAnnulus(int deg)
{
    gsGeometry<>::uPtr quann = gsNurbsCreator<>::NurbsQuarterAnnulus();

    gsKnotVector<> KV1(0,1, 0, deg+1);        // no interior knots in x direction
    gsKnotVector<> KV2(0,1, 0, deg+1);        // no interior knot in y direction

    gsTensorBSplineBasis<2> tbsp (give(KV1), give(KV2));
    gsMatrix<real_t> eval = quann->eval(tbsp.anchors());
    gsGeometry<>::uPtr approxGeom = tbsp.interpolateAtAnchors( eval );
    gsMultiPatch<real_t> mp(*approxGeom);

    //gsMultiPatch<real_t> res = mp.uniformSplit();
    return mp;
}

template<typename T> void printField(T t, const int& width, int prec = 3)
{
    gsInfo << std::left << std::setw(width) << std::setfill(' ') <<std::setprecision(prec)<< t;
}
void printElement(unsigned globDof, unsigned lagMult, unsigned locDofs, real_t Hh,const gsMatrix<index_t>& iter,const gsMatrix<real_t>& cond)
{
    printField(globDof,8);
    printField(lagMult,6);
    printField(locDofs,6);
    printField(cast<real_t,int>(Hh),6);
    for(index_t i=0; i<iter.cols();i++)
        printField(iter(0,i),12);

    for(index_t i=0; i<cond.cols();i++)
        printField(cond(0,i),12);

    gsInfo<<"\n";

}


void prepareInitialGeometry(gsMultiPatch<real_t> & patches,const int testcase, bool testhihj)
{
    switch(testcase)
    {
    case 1:
    case 2:
    case 3:
        patches.degreeElevate();
        patches.uniformRefine();
        break;
    case 4:
        if(patches.parDim()==3)
        {
            for(size_t np=0; np<patches.nPatches();np++)
            {
                gsSparseMatrix<real_t, RowMajor> transfer;
                gsSparseMatrix<real_t, RowMajor> transfers[3];
                gsMatrix<>& coefs = patches.patch(np).coefs();
                patches.patch(np).basis().component(0).uniformRefine_withTransfer(transfers[0],0,0);
                patches.patch(np).basis().component(1).uniformRefine_withTransfer(transfers[1],0,0);
                patches.patch(np).basis().component(2).uniformRefine_withTransfer(transfers[2]);
                tensorCombineTransferMatrices<3,real_t>(transfers,transfer);
                coefs = transfer * coefs;
            }
        }
        break;
    case 5:
        patches.degreeElevate();
        patches.uniformRefine(1);
        break;
    }


    if(!testhihj)
    {
        switch(testcase){
        //    case 2:
        //         for(size_t np = 0; np<patches.nPatches();np++)
        //             if(np%2==0)
        //                 patches.patch(np).uniformRefine();
        break;
        case 1:
        case 7:
            patches.patch(1).uniformRefine();
            patches.patch(2).uniformRefine();
            patches.patch(5).uniformRefine();
            patches.patch(6).uniformRefine();
            break;
        case 3:

            for(size_t np = 0; np<patches.nPatches();np++)
                if(np==0 || np == 2|| np == 4||np==6||np==8)
                    patches.patch(np).uniformRefine();

            break;
        case 4:
            for(size_t np = 0; np<patches.nPatches();np++)
                if(np%3 == 0 && patches.parDim()!=3)
                    patches.patch(np).uniformRefine();
            break;
        case 5:
        {
            gsMatrix<> & coef1 = patches.patch(0).coefs();
            for(index_t i = 0; i<coef1.rows();i++)
            {
                if(coef1(i,0)==-2 && coef1(i,1)==0)
                    coef1(i,1) =0.2;
                if(coef1(i,0)==-2 && coef1(i,1)==1)
                    coef1(i,1) =1.2;
                if(coef1(i,0)==-1 && coef1(i,1)==0)
                    coef1(i,1) =-0.2;
                if(coef1(i,0)==-1 && coef1(i,1)==1)
                    coef1(i,1) =0.8;
            }

            gsMatrix<> & coef2 = patches.patch(1).coefs();
            for(index_t i = 0; i<coef2.rows();i++)
            {
                if(coef2(i,0)==-1 && coef2(i,1)==0)
                    coef2(i,1) =-0.2;
                if(coef2(i,0)==-1 && coef2(i,1)==1)
                    coef2(i,1) =0.8;
            }

            gsMatrix<> & coef0 = patches.patch(2).coefs();
            for(index_t i = 0; i<coef0.rows();i++)
            {
                if(coef0(i,0)==1 && coef0(i,1)==0)
                    coef0(i,1) =0.2;
                if(coef0(i,0)==1 && coef0(i,1)==1)
                    coef0(i,1) =1.2;
            }
            gsMatrix<> & coef3 = patches.patch(3).coefs();
            for(index_t i = 0; i<coef3.rows();i++)
            {
                if(coef3(i,0)==1 && coef3(i,1)==0)
                    coef3(i,1) =0.2;
                if(coef3(i,0)==1 && coef3(i,1)==1)
                    coef3(i,1) =1.2;
                if(coef3(i,0)==2 && coef3(i,1)==0)
                    coef3(i,1) =-0.2;
                if(coef3(i,0)==2 && coef3(i,1)==1)
                    coef3(i,1) =0.8;
            }

            short_t dim = patches.parDim();
            for(size_t np = 2; np<patches.nPatches();np++)
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
                    p->insertKnot(0.25,0,1);
                    p->insertKnot(0.5,0,1);
                    p->insertKnot(0.75,0,1);
                }

            }
            if(dim==2)
            {
                patches.patch(1).uniformRefine();
                // }
                patches.patch(3).uniformRefine();
            }
            break;
        }
        case 6:
            patches.patch(1).uniformRefine();
            break;
        case 8:
        {
            patches.uniformRefine();
            break;
        }
        }


    }

}
void makeGapInGeometry(gsMultiPatch<> & patches, int  testcase, const gsMultiPatch<> & backup, real_t gap,int dim)
{
    for(size_t np=0;np<patches.nPatches();np++)
    {
        gsMatrix<> & coef = patches.patch(np).coefs(); //get reference
        coef = backup.patch(np).coefs(); //override with backup (THIS IS A COPY!)
    }

    switch(testcase){
    case 1:
        //case 2:
    case 3:
    {
        int k=0;
        gsVector<index_t> sgn(2);
        sgn<<1,-1;
        for(gsMultiPatch<>::iiterator iit = patches.iBegin(); iit!= patches.iEnd();++iit)
        {
            gsMatrix<index_t> bb1 = patches.basis(iit->first().patch).boundary(iit->first());
            gsMatrix<index_t> bb2 = patches.basis(iit->second().patch).boundary(iit->second());

            gsMatrix<> & coef1 = patches.patch(iit->first().patch).coefs(); //get reference
            gsMatrix<> & coef2 = patches.patch(iit->second().patch).coefs(); //get reference

            for(index_t i = 1; i<bb1.rows()-1;i++)
            {
                coef1((bb1)(i, 0), iit->first().direction()) +=
                    -sgn(k % 2) * gap * (
                        math::abs_diff(2 * i, bb1.rows() - 1) < 2 // math::abs(i - (bb1.rows() - 1.) / 2) < 0.75
                            && bb1.rows() > 3 ? 0.3 : 1);
            }
            for(index_t i = 1; i<bb2.rows()-1;i++)
            {
                coef2((bb2)(i, 0), iit->second().direction()) +=
                    sgn(k % 2) * gap * (
                        math::abs_diff(2 * i, bb2.rows() - 1) < 2 // math::abs(i - (bb2.rows() - 1.) / 2) < 0.75
                            && bb1.rows() > 3 ? 0.3 : 1); // bb1.rows not bb2.rows?
            }
            k++;
        }
        break;
    }
    case 5:
    {
        gsMatrix<> sign(1,4) ;
        sign<<1,1, -1,-1;
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
        break;
    }
    case 8:
    {
        gsVector<real_t> normal;
        unsigned evFlags = NEED_NORMAL;
        gsMatrix<real_t> midpoint(dim,1);

        bool firstTime = true;

        for(gsMultiPatch<>::iiterator iit = patches.iBegin(); iit!= patches.iEnd();++iit)
        {
            if((iit->first().patch % 2 !=0 &&  iit->second().patch % 2 == 0) || (iit->first().patch % 2 ==0  && iit->second().patch % 2 != 0) )
            {
                if(firstTime)
                {
                    patchSide s1=iit->first();
                    gsGeometryEvaluator<real_t>::uPtr geoEval(getEvaluator(evFlags, *(patches.patches()[s1.patch])));
                    midpoint.setConstant(dim,1,0.5);
                    midpoint(s1.direction(),0)=(int)s1.parameter();

                    geoEval->evaluateAt(midpoint);

                    geoEval->outerNormal(0, s1, normal);
                    normal.normalize();
                    firstTime=false;
                }



                gsMatrix<index_t> bb1 = patches.basis(iit->first().patch).boundary(iit->first());
                gsMatrix<index_t> bb2 = patches.basis(iit->second().patch).boundary(iit->second());

                gsBasis<>::uPtr b1 = patches.basis(iit->first().patch).boundaryBasis(iit->first());
                gsBasis<>::uPtr b2 = patches.basis(iit->second().patch).boundaryBasis(iit->second());

                gsMatrix<index_t> bbb1 = b1->allBoundary();
                gsMatrix<index_t> bbb2 = b2->allBoundary();

                gsMatrix<> & coef1 = patches.patch(iit->first().patch).coefs(); //get reference
                gsMatrix<> & coef2 = patches.patch(iit->second().patch).coefs(); //get reference


                for(index_t i = 0; i<bb1.rows();i++)
                {
                    bool interior=true;
                    for(index_t j=0; j<bbb1.rows();j++)
                        if((bbb1)(j,0)==i && math::abs(coef1((bb1)(i,0),0)-1.06066)>1.e-4 && math::abs(coef1((bb1)(i,0),1)-1.06066)>1.e-4 )
                        {
                            interior=false;
                            break;
                        }

                    if(interior && iit->first().patch % 2 == 1)
                        for(int d=0;d<dim;d++)
                            coef1((bb1)(i,0),d)  += gap*normal(d);

                }
                for(index_t i = 0; i<bb2.rows();i++)
                {
                    bool interior=true;
                    for(index_t j=0; j<bbb2.rows();j++)
                        if((bbb2)(j,0)==i && math::abs(coef2((bb2)(i,0),0)-1.06066)>1.e-4 && math::abs(coef2((bb2)(i,0),1)-1.06066)>1.e-4  )
                        {
                            interior=false;
                            break;
                        }

                    if(interior && iit->second().patch % 2 == 1)
                        for(int d=0;d<dim;d++)
                            coef2((bb2)(i,0),d)  += gap*normal(d);

                }

            }
        }
        break;
    }
    case 4:
    default:
    {
        int k=0;
        gsVector<index_t> sgn(2);
        sgn<<1,-1;
        gsVector<real_t> normal;
        unsigned evFlags = NEED_NORMAL;
        gsMatrix<real_t> midpoint(dim,1);


        for(gsMultiPatch<>::iiterator iit = patches.iBegin(); iit!= patches.iEnd();++iit)
        {
            patchSide s1=iit->first();
            gsGeometryEvaluator<real_t>::uPtr geoEval(getEvaluator(evFlags, *(patches.patches()[s1.patch])));
            midpoint.setConstant(dim,1,0.5);
            midpoint(s1.direction(),0)=(int)s1.parameter();

            geoEval->evaluateAt(midpoint);

            geoEval->outerNormal(0, s1, normal);
            normal.normalize();

            gsMatrix<index_t> bb1 = patches.basis(iit->first().patch).boundary(iit->first());
            gsMatrix<index_t> bb2 = patches.basis(iit->second().patch).boundary(iit->second());

            gsBasis<>::uPtr b1 = patches.basis(iit->first().patch).boundaryBasis(iit->first());
            gsBasis<>::uPtr b2 = patches.basis(iit->second().patch).boundaryBasis(iit->second());

            gsMatrix<index_t> bbb1 = b1->allBoundary();
            gsMatrix<index_t> bbb2 = b2->allBoundary();

            gsMatrix<> & coef1 = patches.patch(iit->first().patch).coefs(); //get reference
            gsMatrix<> & coef2 = patches.patch(iit->second().patch).coefs(); //get reference


            for(index_t i = 0; i<bb1.rows();i++)
            {
                bool interior=true;
                for(index_t j=0; j<bbb1.rows();j++)
                    if((bbb1)(j,0)==i)
                    {
                        interior=false;
                        break;
                    }

                if(interior)
                    for(int d=0;d<dim;d++)
                        coef1((bb1)(i,0),d)  += sgn(k%2)*gap*normal(d);

            }
            for(index_t i = 0; i<bb2.rows();i++)
            {
                bool interior=true;
                for(index_t j=0; j<bbb2.rows();j++)
                    if((bbb2)(j,0)==i)
                    {
                        interior=false;
                        break;
                    }

                if(interior)
                    for(int d=0;d<dim;d++)
                        coef2((bb2)(i,0),d)  += -sgn(k%2)*gap*normal(d);

            }

            k++;
        }
        break;
    }

    }/*end switch*/
    /*
    switch(testcase){
    case 1:


        for(size_t np=0;np<patches.nPatches();np++)
        {
            if(np==0||np==2)
            {
                gsVector<index_t> sgn(3);
                sgn << 1,0, -1;
                gsMatrix<> & coef = patches.patch(np).coefs(); //get reference
                coef = backup.patch(np).coefs(); //override with backup (THIS IS A COPY!)
                for(int i = 0; i<coef.rows();i++)
                {
                    if(coef(i,0)==1)
                    {
                        bool interior = true;
                        for(int k=1; k<dim;k++)
                            if(!(math::abs(coef(i,k))>0.1 && math::abs(coef(i,k))<0.9))
                                interior=false;
                        if(interior)
                            coef(i,0) +=-sgn(np)*gap;
                    }

                }
            }
            if(np==2||np==3)
            {
                gsVector<index_t> sgn(2);
                sgn << 1, -1;
                gsMatrix<> & coef = patches.patch(np).coefs(); //get reference
                coef = backup.patch(np).coefs(); //override with backup (THIS IS A COPY!)
                for(int i = 0; i<coef.rows();i++)
                {
                    if(coef(i,1)==1)
                    {
                        bool interior = true;
                        if(!(math::abs(coef(i,0))>1.1 && math::abs(coef(i,0))<1.9))
                                interior=false;
                        if(interior)
                            coef(i,1) +=sgn(np%2)*gap;
                    }
                }
            }
        }
        break;
    }
    */
}

real_t getInitialMeshsize(int testcase)
{
    switch(testcase){
    case 1:
    case 2:
    case 3:
        return 0.25;
    case 4:
        return 0.1;
    case 5:
        return 0.21;
    case 6:
        return 6;
    }
    return 1;
}


int main (int argc, char** args)
{
    index_t testcase =1;
    index_t max_ref = 2;
    index_t max_elev = 0;
    index_t max_inc = 0;
    index_t jump = 0;
    real_t dom_size =1;
    index_t nPrec =4;
    index_t maxNGaps = 1;
    bool plot = false;
    bool timings = false;
    bool test3D=false;
    bool testhihj = false;
    bool overlap = false;

    bool  comparism = false;

    index_t solCase = 0;
    index_t iniScal = 2;
    index_t initRef = 0;
    index_t oSgn = 1;

    std::string m = "B";
    std::string scaling = "coeff";

    gsCmdLine cmd("Hi, I will test the convergence behaviour of IETI");
    cmd.addSwitch("d", "dim", "test3D", test3D);
    // cmd.addSwitch("p","h-ref test or p-ref test",testDegree);
    cmd.addSwitch("p","plot", "Plot solution", plot);
    cmd.addSwitch("", "IETI.ExtraTiming", "enable extraTimings", timings);
    cmd.addSwitch("n", "nonconf", "Refine only certain patches, leads to a more and more nonconforming mesh",testhihj);
    cmd.addSwitch("o", "overlap", "overlap", overlap);
    cmd.addSwitch("comp","comparism",comparism);
    cmd.addInt("r","refine","Number of refinements",max_ref);
    cmd.addInt("c","case","Choosen testcase",testcase);
    cmd.addInt("e","elevate","maximal elevation",max_elev);
    cmd.addInt("i","increase","maximal elevation (keeping multiplicity)",max_inc);
    cmd.addInt("j","jumps","jumping coeff from 10^-j to 10^j",jump);
    cmd.addString("s", "IETI.Strategy", "Choosen strategy",m);
    cmd.addString("",  "IETI.Scaling", "chosen scaling",scaling);
    cmd.addInt("t","solcase", "Choosen solutioncase",solCase);
    cmd.addInt("f","initRefine","initial refine",initRef);
    cmd.addInt("g","gaps","max number of gaps",maxNGaps);

    cmd.addInt("a", "initGap", "init scaling for gap ", iniScal);

    try { cmd.getValues(argc,args); } catch (int rv) { return rv; }

    gsOptionList options = cmd.getOptionList();

    gsVector<real_t> lambdas(4+max_elev);
    for(int i=0; i<3+max_elev-1;i++)
        lambdas[i]=i+1;
    lambdas[3+max_elev-1]=3+max_elev-1+0.5;
    lambdas[3+max_elev]=3+max_elev;

    // gsVector<real_t> lambdas(1);
    //  lambdas[0]=0;

    gsVector<real_t> gaps(8);
    gaps << 0, 0.1, 0.01,0.005, 0.001,0.0001,0.00005, 0.00001;

    //test3D = true;
    int dim;
    test3D ? dim =3:dim =2;
    gsInfo<<"Dimension of Domain: "<<dim<<"\n";
    gsStopwatch time;

    std::stringstream ss;
    ss << jump;
    std::string jump_str = ss.str();
    gsInfo<<"Using jumping coeff: 10^"+jump_str+" - 10^-"+jump_str<<"\n";
    gsInfo<<"Using method "+m<<"\n";
    //--------------------------------


    /////////////////// Setup geometry ///////////////////
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;

    gsInfo<<"Preparing patch"<<"\n";

    switch(testcase){

    case 2:
    {
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,3, dom_size,-2);
        gsFileData<> fd;
        fd<< patches ;
        fd.dump("makeMultipatch_output");
    }
        break;
    case 1:
    case 7:
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,2, dom_size);
        break;
    case 3:
        patches = gsNurbsCreator<>::BSplineSquareGrid(3,3, dom_size);
        break;
    case 4:
    {
        gsFileData<> fileData("yeti_mp2.xml");
        if(!fileData.getFirst(patches))
            return 1;
        break;
    }
    case 5:
    {
        /*
        gsFileData<> fileData("gaps/2pSquare_lin_0.xml");
        if (fileData.has< gsMultiPatch<> >())
            patches = fileData.getFirst< gsMultiPatch<> >();
        else
            return 1;
            */
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,1, dom_size,-2);
        break;
    }
    case 6:
    {
        gsFileData<> fileData("planar/bumper.xml");
        if (!fileData.getFirst(patches))
            return 1;
        patches.computeTopology();
        break;
    }
    case 8:
    {
        patches = approximateQuarterAnnulus(2);
        patches = patches.uniformSplit();
        break;
    }
    default:
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,2, dom_size);
        break;
    }
    //////// Right-hand side and analytical solution ////////

    gsInfo<<"Preparing diffusion coefficient"<<"\n";
    gsPiecewiseFunction<real_t>  alpha;
    gsFunctionExpr<> a1("1.e-"+jump_str,dim);
    gsFunctionExpr<> a2("1.e+"+jump_str,dim);

    gsFunctionExpr<> a3("1",dim);
    gsFunctionExpr<> a4("3*pi/4",dim);

    switch(testcase){
    case 1:
        for(size_t np = 0; np<patches.nPatches();np++)
            if(np==0||np==3||np==4||np==7)
                alpha.addPiece(a1);
            else
                alpha.addPiece(a2);
        break;
    case 2:
    {
        real_t rho = 1.5* EIGEN_PI;
        real_t rho2= 2;
        if(dim==3)
        {
            rho=1*EIGEN_PI;
            rho2=1.5;
        }

        std::ostringstream strs;
        strs << rho;
        std::string rho_str = strs.str();
        strs.str("");
        strs.clear();
        strs<< rho2;
        std::string rho2_str = strs.str();
        gsFunctionExpr<> a5(rho_str,dim);
        gsFunctionExpr<> a6(rho2_str,dim);
        for(size_t np = 0; np<patches.nPatches();np++)
            if(np<math::floor(patches.nPatches()/2.))
                alpha.addPiece(a5);
            else
                alpha.addPiece(a6);
        break;
    }
    case 3:
        for(size_t np = 0; np<patches.nPatches();np++)
            if (np < 3)  //(np/3 == 0)
                alpha.addPiece(a3);
            else if (np < 6) //(np/3 == 1)
                alpha.addPiece(a4);
            else
                alpha.addPiece(a2);
        break;
    case 5:
    {

        real_t rho = 3* EIGEN_PI;
        real_t rho2= 3;
        if(dim==3)
        {
            rho=2*EIGEN_PI;
            rho2=2;
        }
        std::ostringstream strs;
        strs << rho;
        std::string rho_str = strs.str();
        strs.str("");
        strs.clear();
        strs<< rho2;
        std::string rho2_str = strs.str();
        gsFunctionExpr<> a5(rho_str,dim);
        gsFunctionExpr<> a6(rho2_str,dim);
        alpha.addPiece(a5);
        alpha.addPiece(a5);
        alpha.addPiece(a6);
        alpha.addPiece(a6);

    }
        break;
    case 6:
    {
        gsFunctionExpr<> a5("1",dim);
        gsFunctionExpr<> a6("pi/2",dim);
        alpha.addPiece(a5);
        alpha.addPiece(a6);
    }
        break;
    case 8:
    {
        for(int i=0; i<4;++i)
            alpha.addPiece(a1);
    }
    default:
        for(size_t np = 0; np<patches.nPatches();np++)
            if(np%2 == 0)
                alpha.addPiece(a1);
            else
                alpha.addPiece(a2);
        break;
    }

    // Define source function
    gsFunctionExpr<> f,g;

    switch(solCase)
    {
    case 0:
    {
        if(dim==2)
        {
            f= gsFunctionExpr<>("((pi*1)^2 + (pi*2)^2)*sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)",dim);
            g= gsFunctionExpr<>("sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)+x+y",dim);
        }
        else
        {
            f= gsFunctionExpr<>("((pi*1)^2 + (pi*2)^2+(0.5*pi)^2)*sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)*sin(0.5*pi*(z+0.6))",dim);
            g= gsFunctionExpr<>("sin(pi*(x+0.4))*sin(pi*(y+0.3)*2)*sin(0.5*pi*(z+0.6))",dim);
        }

        break;

    }
    case 1:
    {
        std::string r1 = a1.expression();
        std::string r2 = a2.expression();
        std::string r12 = r1+"*"+r2;

        real_t freq = 4* EIGEN_PI;

        std::ostringstream strs;
        strs << freq;
        std::string freq_s = strs.str();
        /*
        g= gsFunctionExpr<>("if(x<1 & y <1, "+r2+"*sin("+freq_s+"*(x-0))*sin("+freq_s+"*(y-0)),"+
                            "if(x <1 & y<2 , "+r1+"*sin("+freq_s+"*(x-0))*sin("+freq_s+"*(y-1)),"+
                            "if(x < 2 & y< 1, "+r1+"*sin("+freq_s+"*(x-1))*sin("+freq_s+"*(y-0)),"+
                            "if(x < 2 & y< 2, "+r2+"*sin("+freq_s+"*(x-1))*sin("+freq_s+"*(y-1)),"+
                            "if(x < 3 & y< 1, "+r2+"*sin("+freq_s+"*(x-2))*sin("+freq_s+"*(y-0)),"+
                            "if(x < 3 & y< 2, "+r1+"*sin("+freq_s+"*(x-2))*sin("+freq_s+"*(y-1)),"+
                            "if(x < 4 & y< 1, "+r1+"*sin("+freq_s+"*(x-3))*sin("+freq_s+"*(y-0)),"+
                            "if(x < 4 & y< 2, "+r2+"*sin("+freq_s+"*(x-3))*sin("+freq_s+"*(y-1)),0"+
                            "))))))))"
                            ,dim);
        f= gsFunctionExpr<>("if(x<1 & y<1, "+r12+"*2*"+freq_s+"^2*sin("+freq_s+"*(x-0))*sin("+freq_s+"*(y-0)),"+
                            "if(x <1 & y<2 , "+r12+"*2*"+freq_s+"^2*sin("+freq_s+"*(x-0))*sin("+freq_s+"*(y-1)),"+
                            "if(x < 2 & y< 1, "+r12+"*2*"+freq_s+"^2*sin("+freq_s+"*(x-1))*sin("+freq_s+"*(y-0)),"+
                            "if(x < 2 & y< 2, "+r12+"*2*"+freq_s+"^2*sin("+freq_s+"*(x-1))*sin("+freq_s+"*(y-1)),"+
                            "if(x < 3 & y< 1, "+r12+"*2*"+freq_s+"^2*sin("+freq_s+"*(x-2))*sin("+freq_s+"*(y-0)),"+
                            "if(x < 3 & y< 2, "+r12+"*2*"+freq_s+"^2*sin("+freq_s+"*(x-2))*sin("+freq_s+"*(y-1)),"+
                            "if(x < 4 & y< 1, "+r12+"*2*"+freq_s+"^2*sin("+freq_s+"*(x-3))*sin("+freq_s+"*(y-0)),"+
                            "if(x < 4 & y< 2, "+r12+"*2*"+freq_s+"^2*sin("+freq_s+"*(x-3))*sin("+freq_s+"*(y-1)),0"+
                            "))))))))"
                            ,dim);
        */
        g= gsFunctionExpr<>("v:=-sin(pi*x)*sin(pi*y);if(x < 1 & y< 1, "+r2+"*(sin("+freq_s+"*(x-0))*sin("+freq_s+"*(y-0))+v),"+
                            "if(x < 1 & y< 2, "+r1+"*(sin("+freq_s+"*(x-0))*sin("+freq_s+"*(y-1))+v),"+
                            "if(x < 2 & y< 1, "+r1+"*(sin("+freq_s+"*(x-1))*sin("+freq_s+"*(y-0))+v),"+
                            "if(x < 2 & y< 2, "+r2+"*(sin("+freq_s+"*(x-1))*sin("+freq_s+"*(y-1))+v),"+
                            "if(x < 3 & y< 1, "+r2+"*(sin("+freq_s+"*(x-2))*sin("+freq_s+"*(y-0))+v),"+
                            "if(x < 3 & y< 2, "+r1+"*(sin("+freq_s+"*(x-2))*sin("+freq_s+"*(y-1))+v),"+
                            "if(x < 4 & y< 1, "+r1+"*(sin("+freq_s+"*(x-3))*sin("+freq_s+"*(y-0))+v),"+
                            "if(x < 4 & y< 2, "+r2+"*(sin("+freq_s+"*(x-3))*sin("+freq_s+"*(y-1))+v),0"+
                            "))))))))",dim);
        f= gsFunctionExpr<>(r12+"*2*"+freq_s+"^2*sin("+freq_s+"*(x-0))*sin("+freq_s+"*(y-0))-"+
                            "2*pi^2*sin(pi*x)*sin(pi*y)",dim);
        break;
    }
    case 2:
    {
        GISMO_ASSERT(patches.nPatches()==4, "This rhs is only made for 2 patches, choose another one.");

        real_t rho = 3* EIGEN_PI;
        if(dim==3)
            rho=2*EIGEN_PI;
        std::ostringstream strs;
        strs << rho;
        std::string freq1 = strs.str();
        strs.str("");
        strs.clear();
        if(dim==2)
        {
            strs << 3;
            std::string freq2 = strs.str();
            f =  gsFunctionExpr<>("v:="+freq1+";u:="+freq2+";if(x<0,(1+u^2)*pi^2*v*sin(pi*(u*x+y)),pi^2*(1+v^2)*u*sin(pi*(v*x+y)))",dim);
            g =  gsFunctionExpr<>("v:="+freq1+";u:="+freq2+"; if(x<0,sin(pi*(u*x+y)),sin(pi*(v*x+y)))",dim);
        }
        else
        {
            strs << 2;
            std::string freq2 = strs.str();
            f =  gsFunctionExpr<>("v:="+freq1+";u:="+freq2+";if(x<0,(2+u^2)*pi^2*v*sin(pi*(u*x+y+z)),pi^2*(2+v^2)*u*sin(pi*(v*x+y+z)))",dim);
            g =  gsFunctionExpr<>("v:="+freq1+";u:="+freq2+"; if(x<0,sin(pi*(u*x+y+z)),sin(pi*(v*x+y+z)))",dim);
        }
        break;
    }
    case 3:
    {
        std::string r1 = a1.expression();
        std::string r2 = a2.expression();

        real_t freq_3 = 3;
        real_t gamma = 14*pow(10.,jump)-8.;

        std::ostringstream strs;
        strs << freq_3;
        std::string freq_s = strs.str();
        strs.str("");
        strs.clear();

        strs << gamma;
        std::string gamma_s = strs.str();
        g= gsFunctionExpr<>("v:=sin(7*pi/2*y);u:="+freq_s+"; if(x < 1, exp(-sin(u*pi*x))*v,"+
                            "if(x < 2, ((2*x^2-1)+(x-1)^2*(x-2)*"+gamma_s+")*v,"+
                            "7*exp(-cos((x*pi/2+pi/2)*u))*v"+
                            "))",dim);
        f= gsFunctionExpr<>("v:=sin(7*pi/2*y);u:="+freq_s+"; if(x < 1, -exp(-sin(u*pi*x))*v/4*pi^2*(-49+4*u^2*cos(u*pi*x)^2+4*u^2*sin(u*pi*x)),"+
                            "if(x < 2, u*pi*v/16*(-16+"+gamma_s+"*(32+49*pi^2*(x-2)*(x-1)^2-24*x)+49*pi^2*(2x^2-1)),"+
                            "-(8+"+gamma_s+")*pi^2*exp(-cos(u*pi*(x+1)/2))*v/8*(-49+u^2*cos(u*pi*(x+1)/2)+u^2 *sin(u*pi*(x+1)/2)^2 )"+
                            "))",dim);
        break;
    }
    case 4:
    {
        real_t rho = 1.5* EIGEN_PI;
        if(dim==3)
            rho=1*EIGEN_PI;
        std::ostringstream strs;
        strs << rho;
        std::string freq1 = strs.str();
        strs.str("");
        strs.clear();
        if(dim==2)
        {
            strs << 2;
            std::string freq2 = strs.str();
            f =  gsFunctionExpr<>("v:="+freq1+";u:="+freq2+";if(x<0,(1+u^2)*pi^2*v*sin(pi*(u*x+y)),pi^2*(1+v^2)*u*sin(pi*(v*x+y)))",dim);
            g =  gsFunctionExpr<>("v:="+freq1+";u:="+freq2+"; if(x<0,sin(pi*(u*x+y)),sin(pi*(v*x+y)))",dim);
        }
        else
        {
            strs << 1.5;
            std::string freq2 = strs.str();
            f =  gsFunctionExpr<>("v:="+freq1+";u:="+freq2+";if(x<0,(2+u^2)*pi^2*v*sin(pi*(u*x+y+z)),pi^2*(2+v^2)*u*sin(pi*(v*x+y+z)))",dim);
            g =  gsFunctionExpr<>("v:="+freq1+";u:="+freq2+"; if(x<0,sin(pi*(u*x+y+z)),sin(pi*(v*x+y+z)))",dim);
        }
        break;
    }
    case 5:
    {
        f =  gsFunctionExpr<>("v:=x+y;if(v>0,2*pi/2*exp(sin(v))*(sin(v)-cos(v)*cos(v)),2/4*pi*pi*sin(pi/2*v))",dim);
        g =  gsFunctionExpr<>("v:=x+y; if(v>0,(exp(sin(v))-1),sin(pi/2*v))",dim);
        if(dim==2)
        {
            f =  gsFunctionExpr<>("v:=x+y;if(v>0,2*pi/2*exp(sin(v))*(sin(v)-cos(v)*cos(v)),2/4*pi*pi*sin(pi/2*v))",dim);
            g =  gsFunctionExpr<>("v:=x+y; if(v>0,(exp(sin(v))-1),sin(pi/2*v))",dim);
        }
        break;
    }
    case 7:
    {
        f= gsFunctionExpr<>("0",dim);
        g= gsFunctionExpr<>("0",dim);
        break;
    }

    }
    //Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<"\n" << "\n";

    // for(size_t np = 0; np<patches.nPatches();np++)
    //    gsInfo<<"Alpha on "<<np<<" "<< alpha[np]<<"\n";



#ifdef _OPENMP
    // To get the number of threads from the environment
    int nThreads;
    const char * nProcs = getenv("OMP_NUM_THREADS");
    if(nProcs != NULL)
        sscanf( nProcs, "%d", &nThreads );
    else
        nThreads = 1; //one OpenMP thread
    nThreads = math::min(nThreads,(int)patches.nPatches());
    gsInfo<<"Using "<<nThreads<<" threads\n";
    std::ostringstream oss;
    oss << nThreads;

    setenv("OMP_NUM_THREADS", oss.str().c_str(),1);
#endif


    /////////////////// Setup boundary conditions ///////////////////
    // Define Boundary conditions
    gsInfo<<"Setting BC"<<"\n";
    gsBoundaryConditions<> bcInfo;

    gsFunctionExpr<> hEast, hSouth,hNorth, hWest;

    hWest = gsFunctionExpr<>("-pi*cos(pi*0.4)*sin(2*pi*(y+0.3))-1",dim);
    hSouth = gsFunctionExpr<>("-pi*2*cos(2*pi*0.3)*sin(pi*(x+0.4))-1",dim);

    switch(testcase){
    case 1:
        if(solCase==1)
        {
            for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
                bcInfo.addCondition(*it, condition_type::dirichlet, g);
        }
        else
        {
            hEast = gsFunctionExpr<>("1*pi*cos(pi*(4+0.4))*sin(2*pi*(y+0.3))+1",dim);
            hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(2+0.3))*sin(pi*(x+0.4))+1",dim);

            bcInfo.addCondition(6, boundary::east, condition_type::neumann, hEast);
            bcInfo.addCondition(7, boundary::east, condition_type::neumann, hEast);

            bcInfo.addCondition(1, boundary::north, condition_type::neumann, hNorth);
            bcInfo.addCondition(3, boundary::north, condition_type::neumann, hNorth);
            bcInfo.addCondition(5, boundary::north, condition_type::neumann, hNorth);
            bcInfo.addCondition(7, boundary::north, condition_type::neumann, hNorth);

            bcInfo.addCondition(0, boundary::south, condition_type::neumann, hSouth);
            bcInfo.addCondition(2, boundary::south, condition_type::neumann, hSouth);
            bcInfo.addCondition(4, boundary::south, condition_type::neumann, hSouth);
            bcInfo.addCondition(6, boundary::south, condition_type::neumann, hSouth);

            bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
            bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, g);
        }
        break;
    case 2:
        if(solCase==1||solCase==4)
        {
            for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
                // if((*it).patch==1) // otherwise the convergence rates are wrong, since we dont know the exact solution
                bcInfo.addCondition(*it, condition_type::dirichlet, g);
        }
        else
        {
            hEast = gsFunctionExpr<>("1*pi*cos(pi*(4+0.4))*sin(2*pi*(y+0.3))+1",dim);
            hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(2+0.3))*sin(pi*(x+0.4))+1",dim);

            bcInfo.addCondition(6, boundary::east, condition_type::neumann, hEast);
            bcInfo.addCondition(7, boundary::east, condition_type::neumann, hEast);

            bcInfo.addCondition(1, boundary::north, condition_type::neumann, hNorth);
            bcInfo.addCondition(3, boundary::north, condition_type::neumann, hNorth);
            bcInfo.addCondition(5, boundary::north, condition_type::neumann, hNorth);
            bcInfo.addCondition(7, boundary::north, condition_type::neumann, hNorth);

            bcInfo.addCondition(0, boundary::south, condition_type::neumann, hSouth);
            bcInfo.addCondition(2, boundary::south, condition_type::neumann, hSouth);
            bcInfo.addCondition(4, boundary::south, condition_type::neumann, hSouth);
            bcInfo.addCondition(6, boundary::south, condition_type::neumann, hSouth);

            bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, g);
            bcInfo.addCondition(1, boundary::west,  condition_type::neumann, hWest);
        }
        break;
    case 3:
        /*
        hEast = gsFunctionExpr<>("1*pi*cos(pi*(4+0.4))*sin(2*pi*(y+0.3))+1",dim);
        hNorth = gsFunctionExpr<>("pi*2*cos(2*pi*(4+0.3))*sin(pi*(x+0.4))+1",dim);

        for(int i=0;i<4;i++)
        {
            bcInfo.addCondition(i+12, boundary::east, condition_type::neumann, hEast);
            bcInfo.addCondition(i*4+3, boundary::north, condition_type::neumann, hNorth);
            bcInfo.addCondition(i*4, boundary::south, condition_type::neumann, hSouth);
            bcInfo.addCondition(i, boundary::west,  condition_type::dirichlet, g);
        }
        break;*/
    case 4:
    case 5:
    default:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
            // if((*it).patch==1) // otherwise the convergence rates are wrong, since we dont know the exact solution
            bcInfo.addCondition(*it, condition_type::dirichlet, g);
        break;

    }

    if(test3D)
    {
        gsMultiPatch<real_t> patches3D;
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<2,real_t>* tb = static_cast<gsTensorBSpline<2,real_t>* >(&patches.patch(i));
            patches3D.addPatch(gsNurbsCreator<>::lift3D(*tb));
        }
        patches3D.computeTopology();
        patches = give(patches3D);

        gsFunctionExpr<> h("2",dim);
        gsFunctionExpr<> hm("-2",dim);

        if(testcase==5)
            for(size_t np=0; np<patches.nPatches();++np)
            {
                bcInfo.addCondition(np, boundary::front,  condition_type::dirichlet, g);
                bcInfo.addCondition(np, boundary::back,  condition_type::dirichlet, g);
            }

        if(solCase==0)
            for(size_t np=0; np<patches.nPatches();++np)
            {
                //       bcInfo.addCondition(np, boundary::front,  condition_type::neumann, h);
                //      bcInfo.addCondition(np, boundary::back,  condition_type::neumann, hm);
                //    bcInfo.addCondition(np, boundary::front,  condition_type::dirichlet, g);
                //      bcInfo.addCondition(np, boundary::back,  condition_type::dirichlet, g);
            }

    }



    prepareInitialGeometry(patches, testcase, testhihj);

    gsFileData<> fd;
    fd<< patches ;
    fd.dump("makeMultipatch_output");
    real_t h0 = getInitialMeshsize(testcase);
    gsMultiPatch<> backup(patches);

    int result = 0;
    if(plot)
    {
        gsWriteParaview<>( patches, "IETI_patch", 1000, true);
        gsWriteParaview<>( patches, "IETI_patch2_", 1000);
        result = system("paraview IETI_patch.pvd &");
    }
    gsMatrix<> eigs, solVector, rhs, solution;
    gsField<> sol;

    ////////////////////// Refinement h and p //////////////////////
    // Refinement

    gsMatrix<index_t> iter(max_ref,nPrec);
    gsMatrix<real_t> Hh(max_ref,1);
    gsMatrix<real_t> hihj(max_ref,1);
    gsMatrix<index_t> globDof(max_ref,1);
    gsMatrix<index_t> lagMult(max_ref,1);
    gsMatrix<index_t>  locDof(max_ref,1);
    gsMatrix<real_t> cond(max_ref,nPrec);

    gsMatrix<real_t> l2Error = gsMatrix<real_t>::Zero(max_ref,maxNGaps);
    gsMatrix<real_t> h1error= gsMatrix<real_t>::Zero(max_ref,maxNGaps);
    gsMatrix<real_t> dGerror= gsMatrix<real_t>::Zero(max_ref,maxNGaps);
    gsMatrix<real_t> h= gsMatrix<real_t>::Zero(max_ref,maxNGaps);

    gsMatrix<real_t> l2Factor= gsMatrix<real_t>::Zero(max_ref,maxNGaps);
    gsMatrix<real_t> h1Factor= gsMatrix<real_t>::Zero(max_ref,maxNGaps);
    gsMatrix<real_t> dgFactor= gsMatrix<real_t>::Zero(max_ref,maxNGaps);

    gsMatrix<> Gaps(max_ref,lambdas.rows());
    gsMatrix<real_t> l2ErrorComp = gsMatrix<real_t>::Zero(max_ref,lambdas.rows());
    gsMatrix<real_t> h1errorComp= gsMatrix<real_t>::Zero(max_ref,lambdas.rows());
    gsMatrix<real_t> dGerrorComp= gsMatrix<real_t>::Zero(max_ref,lambdas.rows());

    gsMatrix<real_t> l2FactorComp= gsMatrix<real_t>::Zero(max_ref,lambdas.rows());
    gsMatrix<real_t> h1FactorComp= gsMatrix<real_t>::Zero(max_ref,lambdas.rows());
    gsMatrix<real_t> dgFactorComp= gsMatrix<real_t>::Zero(max_ref,lambdas.rows());

    cond.setZero();
    iter.setZero();
    l2Error.setZero();



    for(int gap = 0; gap<maxNGaps;gap++)
    {
        // Make the gap
        makeGapInGeometry(patches, testcase,backup,oSgn*gaps[gap],dim);

        //Exact solution
        gsField<> gField(patches, g,false);

        //Make the basis for refinement
        //gsMultiBasis<> refine_bases(refine_bases_backup);

        gsMultiBasis<> refine_bases(patches);

        /////////////////// Basis setup ///////////////////
        int max_tmp = refine_bases.maxDegree(0);
        for (int j = 1; j < patches.parDim(); ++j )
            if ( max_tmp < refine_bases.maxDegree(j) )
                max_tmp = refine_bases.maxDegree(j);

        // Elevate all degrees uniformly
        //max_tmp += max_elev;

        refine_bases.setDegree(max_tmp);
        refine_bases.degreeElevate(max_elev);

        for(int i=0; i<initRef;i++)
            refine_bases.uniformRefine();
        gsInfo<<"Degree set to: "<<refine_bases.basis(0).degree(0)<<" , "<<refine_bases.basis(0).degree(1)<<"\n";

        //Start the loop over the refinements
        for(int ref = 0; ref<max_ref;ref++)
        {
            gsInfo<<"Starting ref: "<<ref<<"\n";

            //gsPoissonAssembler<real_t> assembler(patches,refine_bases,bcInfo,*f, dirichlet::nitsche,iFace::dg);
            //gsPoissonPde<real_t>ppde (*patches,bcInfo,*f);
            //gsPoissonAssembler<real_t> assembler(ppde,refine_bases,dirichlet::nitsche,fStrat);

            gsPoissonHeterogeneousPde<real_t>ppde(patches,bcInfo,f,alpha);
            gsPoissonHeterogeneousAssembler<real_t>assembler (ppde,refine_bases,dirichlet::elimination,iFace::dg);


            time.restart();
            gsIETIdGAssembler<real_t> ass(assembler);
            ass.setOptions(options.getGroup("IETI"));
            ass.init();
            time.stop();
            gsInfo<<"Time for preparing the book-keeping: "<<time<<"\n";
            time.restart();
            ass.assemble();
            time.stop();
            gsInfo<<"Time for assembling everything: "<<time<<"\n";

            gsIETIInfo info = ass.getInfo();
            gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(ass);
            solv->init();
            solVector.setZero(ass.systemSize(),ass.numberRhs());
            rhs = solv->getRhs();

            gsScaledDirichletPrecond<real_t>::Ptr precDir =
                    gsScaledDirichletPrecond<real_t>::make(ass);
            gsConjugateGradient<> PCG(solv,precDir);
            PCG.setMaxIterations(200);
            PCG.setCalcEigenvalues(true);

            time.restart();
            PCG.solve(rhs,solVector);
            //gsInfo<<"sol:\n"<<solVector<<"\n";
            gsInfo<<"Time for solving PCG total: "<<time<<"\n";
            gsInfo<<"Time per PCG-Iteration: "<<time<<"/"<<PCG.iterations()<<"\n";
            solv->calculateSolution(solVector,solution);
            sol = assembler.constructSolution(solution);

            l2Error(ref,gap) = sol.distanceL2(g);
            h1error(ref,gap) = sol.distanceH1(g);
            dGerror(ref,gap) = sol.distanceDG(g);

            h(ref,gap) = h0*1. / math::exp2(ref+initRef);

            refine_bases.uniformRefine();
        }


        l2Factor.row(0).setZero();
        h1Factor.row(0).setZero();
        dgFactor.row(0).setZero();
        for(int i=1;i<max_ref;i++)
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


    l2Error.setZero(max_ref,lambdas.rows());
    h1error.setZero(max_ref,lambdas.rows());
    dGerror.setZero(max_ref,lambdas.rows());
    h.setZero(max_ref,lambdas.rows());

    l2Factor.setZero(max_ref,lambdas.rows());
    h1Factor.setZero(max_ref,lambdas.rows());
    dgFactor.setZero(max_ref,lambdas.rows());
    Gaps.setZero(max_ref,lambdas.rows());

    for(int l = 0; l<lambdas.rows();l++)
    {
        gsMultiBasis<> refine_bases(patches);
        real_t lambda = lambdas(l);
        // Elevate and p-refine the basis to order k + max_elev
        // where k is the highest degree in the bases
        // Find maximum degree with respect to all the variables
        int max_tmp = refine_bases.maxDegree(0);
        for (int j = 1; j < patches.parDim(); ++j )
            if ( max_tmp < refine_bases.maxDegree(j) )
                max_tmp = refine_bases.maxDegree(j);

        // Elevate all degrees uniformly
        //max_tmp += max_elev;

        refine_bases.setDegree(max_tmp);
        refine_bases.degreeElevate(max_elev);

        for(int i=0; i<initRef;i++)
            refine_bases.uniformRefine();

        gsInfo<<"Degree set to: "<<refine_bases.basis(0).degree(0)<<" , "<<refine_bases.basis(0).degree(1)<<"\n";


        for(int r = 0;r < max_ref;r++)
        {
            Gaps(r,l)= oSgn*h0*1./math::pow((real_t)2,lambda*(r+initRef)+iniScal);
            if(lambda == 0)
                Gaps(r,l) = 0;
            makeGapInGeometry(patches, testcase,backup,Gaps(r,l),dim);
            gsField<> gField(patches, g,false);

            gsInfo<<"Starting ref: "<<r<<"\n";

            gsPoissonHeterogeneousPde<real_t>ppde(patches,bcInfo,f,alpha);
            gsPoissonHeterogeneousAssembler<real_t>assembler (ppde,refine_bases,dirichlet::elimination,iFace::dg);


            time.restart();

            gsIETIdGAssembler<real_t> ass(assembler);
            ass.setOptions(options.getGroup("IETI"));
            ass.init();
            time.stop();
            gsInfo<<"Time for preparing the book-keeping: "<<time<<"\n";
            time.restart();
            ass.assemble();
            time.stop();
            gsInfo<<"Time for assembling everything: "<<time<<"\n";

            gsIETIInfo info = ass.getInfo();
            gsIETISolver<real_t>::Ptr solv = gsIETISolver<real_t>::make(ass);
            solv->init();
            const gsBasis<>& B= refine_bases.basis(0);
            if(!test3D)
                Hh(r,0) = math::sqrt( (real_t)(B.numElements()) ); //Hh domsize cancels
            else
                Hh(r,0) = math::pow( (real_t)(B.numElements()), 1.0/3);

            const gsBasis<>& B1= refine_bases.basis(1);

            locDof(r,0) = refine_bases.size(0);

            globDof(r,0)=info.dofTotal;
            lagMult(r,0)=info.lagrangeMult;

            if(!testhihj)
                refine_bases.uniformRefine();
            else
            {
                if(!test3D)
                    hihj(r,0) = math::sqrt( (real_t)(B.numElements())/(real_t)(B1.numElements()) ); //Hh domsize cancels
                else
                    hihj(r,0) = math::pow( (real_t)(B.numElements())/(real_t)(B1.numElements()), 1.0/3);

                for(size_t np = 0; np<patches.nPatches();np++)
                    if(np%3 == 0)
                        refine_bases.basis(np).uniformRefine();
            }

            for(int pre = 0;pre<nPrec;pre++)
            {
                if(pre == 0 && comparism)
                {
                    gsInfo<<"assembling\n"<<std::flush;
                    assembler.assemble();
                    solution.setZero(assembler.numDofs(),1);

#if defined(GISMO_WITH_PARDISO)
                    gsSparseSolver<>::PardisoLLT solver;
#elif defined(GISMO_WITH_SUPERLU)
                    gsSparseSolver<>::SuperLU solver;
#else
                    gsSparseSolver<>::LU solver;
#endif
                    gsInfo<<"solving\n"<<std::flush;
                    //  gsInfo<<"matrix\n"<<assembler.matrix().toDense()<<"\n";
                    // gsInfo<<"rhs\n"<<assembler.rhs()<<"\n";
                    solver.compute(assembler.matrix());
                    solution = solver.solve( assembler.rhs());

                    sol = assembler.constructSolution(solution);

                    l2ErrorComp(r,l) = sol.distanceL2(g);
                    h1errorComp(r,l) = sol.distanceH1(g);
                    dGerrorComp(r,l) = sol.distanceDG(g);
                    continue;
                }
                else if(pre == 0)
                    continue;

                gsInfo<<"starting preconditioner "<<pre<<"\n";
                if(solCase==1)
                    solVector.setRandom(ass.systemSize(),ass.numberRhs());
                else
                    solVector.setZero(ass.systemSize(),ass.numberRhs());
                rhs = solv->getRhs();

                gsScaledDirichletPrecond<real_t>::Ptr precDir =
                        gsScaledDirichletPrecond<real_t>::make(ass);
                gsConjugateGradient<> PCG(solv,precDir);
                PCG.setMaxIterations(200);
                PCG.setCalcEigenvalues(true);

                ass.getOptions().scal = static_cast<IETIPrecondScaling::strategy>(pre+1);

                time.restart();
                PCG.solve(rhs,solVector);
                //gsInfo<<"sol:\n"<<solVector<<"\n";
                gsInfo<<"Time for solving PCG total: "<<time<<"\n";
                gsInfo<<"Time per PCG-Iteration: "<<time<<"/"<<PCG.iterations()<<"\n";



                PCG.getEigenvalues(eigs);
                if(l==0)
                {
                    iter(r,pre)=PCG.iterations();
                    cond(r,pre)=eigs.maxCoeff()/eigs.minCoeff();
                }
                time.restart();
                solv->calculateSolution(solVector,solution);
                time.stop();
                gsInfo<<"Time for reconstructing solution: "<<time<<"\n";
                if(pre==nPrec-1)
                    sol = assembler.constructSolution(solution);
            }
            l2Error(r,l)= sol.distanceL2(g);
            h1error(r,l) = sol.distanceH1(g);
            dGerror(r,l) = sol.distanceDG(g);
            h(r,l) = h0*1. / math::exp2(r+initRef);

            if (plot && l==0)
            {
                // Write approximate and exact solution to paraview files
                gsInfo<<"Plotting in Paraview...\n";

                std::stringstream Str;
                Str << r;
                gsField<> diff =  gsFieldCreator<>::absError(sol,g);


                gsField<> alph( patches, alpha, false );
                gsWriteParaview<>( alph, "IETI_alpha", 1000);

                gsWriteParaview<>(patches, "gapsMesh"+Str.str()+"", 8000, true);
                gsWriteParaview<>(sol, "gaps"+Str.str(), 8000);
                gsWriteParaview<>( gField, "gaps_ex"+Str.str(), 8000);

                gsWriteParaview<>(diff, "diff"+Str.str(), 8000);

                // Run paraview
                std::string string("paraview gaps"+Str.str()+".pvd &");
                result = system(string.c_str());

            }
        }
    }
    l2Factor.row(0).setZero();
    h1Factor.row(0).setZero();
    dgFactor.row(0).setZero();
    for(int i=1;i<max_ref;i++)
    {
        for(index_t l= 0; l<lambdas.rows();l++)
        {
            l2Factor(i,l) = math::log(l2Error(i-1,l)/l2Error(i,l))/math::log(2.0);
            h1Factor(i,l) = math::log(h1error(i-1,l)/h1error(i,l))/math::log(2.0);
            dgFactor(i,l) = math::log(dGerror(i-1,l)/dGerror(i,l))/math::log(2.0);
        }
    }

    if(comparism)
    {
        l2FactorComp.row(0).setZero();
        h1FactorComp.row(0).setZero();
        dgFactorComp.row(0).setZero();
        for(index_t i=1;i<max_ref;i++)
        {
            for(index_t l= 0; l<lambdas.rows();l++)
            {
                l2FactorComp(i,l) = math::log(l2ErrorComp(i-1,l)/l2ErrorComp(i,l))/math::log(2.0);
                h1FactorComp(i,l) = math::log(h1errorComp(i-1,l)/h1errorComp(i,l))/math::log(2.0);
                dgFactorComp(i,l) = math::log(dGerrorComp(i-1,l)/dGerrorComp(i,l))/math::log(2.0);
            }
        }
        gsInfo<<"L2ErrorComp: "<<"\n"<<l2ErrorComp<<"\n"<<"\n";
        gsInfo<<"L2FactorComp: "<<"\n"<<l2FactorComp<<"\n"<<"\n";

        gsInfo<<"H1ErrorComp: "<<"\n"<<h1errorComp<<"\n"<<"\n";
        gsInfo<<"H1FactorComp: "<<"\n"<<h1FactorComp<<"\n"<<"\n";
        gsInfo<<"DgErrorComp: "<<"\n"<<dGerrorComp<<"\n"<<"\n";
        gsInfo<<"DGFactorComp: "<<"\n"<<dgFactorComp<<"\n"<<"\n";
    }
    gsInfo<<"Lambda: "<<lambdas.transpose()<<"\n\n";
    gsInfo<<"Gaps: "<<Gaps.transpose()<<"\n"<<"\n";

    gsInfo<<"L2Error: "<<"\n"<<l2Error<<"\n"<<"\n";
    gsInfo<<"L2Factor: "<<"\n"<<l2Factor<<"\n"<<"\n";

    gsInfo<<"H1Error: "<<"\n"<<h1error<<"\n"<<"\n";
    gsInfo<<"H1Factor: "<<"\n"<<h1Factor<<"\n"<<"\n";
    gsInfo<<"DgError: "<<"\n"<<dGerror<<"\n"<<"\n";
    gsInfo<<"DGFactor: "<<"\n"<<dgFactor<<"\n"<<"\n";

    gsInfo<<"DGError for output: "<<"\n";
    for(index_t i=0; i< dGerror.rows();i++)
        gsInfo<<h(i,0)<<"   "<<dGerror.row(i)<<"\n";

    gsInfo<<"DGFactor for output: "<<"\n";
    for(index_t i=0; i< dgFactor.rows();i++)
        gsInfo<<h(i,0)<<"   "<<dgFactor.row(i)<<"\n";
    gsInfo<<"h: "<<"\n"<<h<<"\n"<<"\n";

    gsInfo << "\n";
    printField("gDof",8);printField("lagM",6);printField("lDof",6); printField("H/h",6); printField("iF",12);printField("iPF-coeff",12);printField("iPF-stiff",12);printField("iPF-stiffM",12);printField("cF",12); printField("cPF-coeff",12);printField("cPF-stiff",12);printField("cPF-stiffM",12);
    gsInfo<<"\n";
    for(int i=0;i<max_ref;i++)
        printElement(globDof(i,0),lagMult(i,0),locDof(i,0),Hh(i,0),iter.row(i),cond.row(i));

    gsInfo<<"\n\n";
    if(testhihj)
        gsInfo <<hihj<<"\n";

    //	return  1 for failures and 0 for success
    return result;
}
