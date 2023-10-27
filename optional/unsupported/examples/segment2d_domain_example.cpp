/** @file segment2d_domain_example.cpp
    
    @brief Takes a planar multiply connected domain as input, a
    template as a reference domain and segment the planar domain
    accordingly to the segmentation of the template.  It needs to
    construct a map between those two objects and approximates its
    inverse via a tracing curve algorithm
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Falini, A. Mantzaflaris
*/


#include <gismo.h>
#include <gismo_dev.h>
#include <gsBem/gsBemLaplace.h>
//#include <gsBem/gsBemUtils.h>
#include <gsModeling/gsTraceCurve.hpp>

using namespace gismo;

template<class T>
std::pair<gsFunction<T>* , gsFunction<T>* >
mapto(gsPlanarDomain<T> & pdomain, gsTemplate<T> & my_template, gsMatrix<T> &Dirichlet,int numRefine);

void quadPointsMultConnected(gsPlanarDomain<real_t> *P, int&n, real_t &l, real_t&scalar, gsMatrix<real_t>&indice, real_t &areaWeight)
{

    int k=0;

    //1. creating mask of a shrinked curve
    gsBSpline<real_t>::uPtr b = memory::convert_ptr<gsBSpline<real_t> >( P->outer().singleCurve() );

    gsMatrix<real_t> mod = b->coefs();
    mod *=scalar;//0.85
    gsKnotVector<real_t> kv = b->basis().knots();

    gsBSpline<real_t> *interior = new gsBSpline<real_t>(kv, mod);


    gsCurveLoop<real_t> * inner = new gsCurveLoop<real_t>( interior );
    gsPlanarDomain<real_t> Idomain(inner);
    std::vector<gsMatrix<real_t> *> Giga;


    gsMatrix<real_t> tab(2,2);
    tab = Idomain.boundingBox();
    gsMatrix<real_t> storage;
    storage= uniformPointGrid<real_t>( tab.col(0), tab.col(1),n);
    indice.resize(2,storage.cols()); // matrix indice stores points of the grid which are
    // inside the planar domain
    indice.setZero();

    //2. inserting holes
    if(P->numHoles()!=0)
    {

        int dof = 0;

        std::vector<gsKnotVector<real_t> *> KV2;
        std::string sub_num;

        for(int v=1; v<=P->numHoles(); ++v)
        {
            dof = P->loop(v).singleCurve()->basis().size();
            gsMatrix<real_t> V(dof,2);

            V = P->loop(v).singleCurve()->coefs();
            //  gsInfo<<"\m matrix of cc for hole \n"<<V<<"\n";

            //*** The following is still NOT WORKING
            //  V *=1.26;//we can play around with this parameter, it depends how close to the original boundary we want to be.
            //   gsInfo<<"\n scaled V is \n"<<V<<"\n";


            Giga.push_back(&V);

            //gsInfo<<"\n matrix block \n "<<*Giga[v-1]<<"\n";

            gsBSpline<real_t>::uPtr Bb = memory::convert_ptr<gsBSpline<real_t> >(P->loop(v).singleCurve());
            KV2.push_back(&(Bb->basis().knots()));//size=to #holes
            //  gsInfo<<"\n knot vector "<<*KV2[v-1]<<"\n"<<"\n";

            gsBSpline<real_t>::uPtr holeI(new gsBSpline<real_t>(*KV2[v-1],*Giga[v-1]));
            std::ostringstream transform;
            transform << v;
            sub_num = transform.str();
            gsWriteParaview(*holeI,"buco"+sub_num,200);
            Idomain.insertHole(new gsCurveLoop<real_t>(holeI.release()));
        }

    }

    index_t j;
    for(j=0;j<storage.cols();j++)
    {
        if(  Idomain.inDomain(storage.col(j),1) && Idomain.inDomain( storage.col(j),0 ) ) //point inside the domain
        {
            indice.col(k).noalias() = storage.col(j);
            // gsInfo<<"indice (0,k) : " << indice(0,k) << "\n k is : "<< k<<"\n";
            k++;
        }
    }
    indice.conservativeResize(2,k);

    areaWeight = l*l;


    //   gsWriteParaview(Idomain,"Idomain",500);
}


template<class T>
void segment(gsPlanarDomain<T> & pdomain,gsMatrix<T> Dirichlet, gsMatrix<T> indice,
             int n_points, T tolerance,bool circle, gsMesh<T> &segmentation)
{
    // ***************************   Consider a template
    gsTemplate<T> * pmy_template;
    if(circle)
        pmy_template = new gsTemplate<T>(pdomain.numHoles() ); // Convex template with no holes

    else
    {
        bool sq=true;
        pmy_template= new gsTemplate<T>(sq,1);
    }

    gsTemplate<T>& my_template (*pmy_template);

    // ***************************   first step: solve Laplace problems


    // call : std::pair< gsFunction<T> *, gsFunction<T> *, > S = pdomain.mapto( template );
    // u.first, u.second are the components of the solution
    std::pair <gsFunction<T>* , gsFunction<T>*> u;

    u=mapto(pdomain,my_template,Dirichlet, 2); // 2 for amoeba_hole, 3 for austria_hole
    // 4 for Puzzle1 domain


    // indice points: points inside the planar domain at a suitable distance from
    // the outer boundary

    int k=0;
    k = indice.cols();
    gsMatrix<T> Image, Image1(1,k), Image2(1,k);
    u.first->eval_into(indice, Image1);
    u.second->eval_into(indice, Image2);
    //    gsInfo<<"\n x coordinates of points mapped : \n"<<Image1<<"\n";
    //    gsInfo<<"\n y coordinates of points mapped : \n"<<Image2<<"\n";
    Image.resize(2, Image1.cols() );
    Image.row(0) = Image1;
    Image.row(1) = Image2;

    gsMatrix<T> medium_point (1,1), middle, x, leftPart, rightPart;
    medium_point<< 0.5 ;
    // ***************************     second step: trace boundary curves, starting from a point
    // ***************************     in the middle, then applying the algorithm
    // ***************************     in both directions

    gsMatrix<T> param;

    for(int s=0;s!=(my_template.skeletonSize());++s)
    {

        middle=my_template.skeleton(s)->eval(medium_point);
        param = my_template.skeleton(s)->parameterRange();
        index_t j=0;
        T dist = (Image.col(j) - middle ).norm();
        for (index_t i=1; i!= Image.cols(); ++i)
        {
            if( ( Image.col(i) - middle ).norm() < dist )
                j = i;
        }//j gives you the index of the closest point to the mid point

        //   gsInfo<<"\n selected indice.col(j) "<<indice.col(j).transpose()<<"\n";
        gsTraceLine<T>(u, indice.col(j), middle, x );
        //  gsInfo<<"\n x from trace line is "<< x.transpose()<<"\n";

        gsInfo<<"\n TRACING Curve "<< s <<"\n\n"<<"\n";
        gsTraceCurvePart<T>(u,x,my_template.skeleton(s), param(0,1)/2, param(0,0), leftPart, n_points, tolerance );
        leftPart = leftPart.rowwise().reverse().eval();

        gsTraceCurvePart<T>(u,x,my_template.skeleton(s), param(0,1)/2, param(0,1), rightPart, n_points, tolerance );

        const index_t lastcol = rightPart.cols()-1;

        leftPart.conservativeResize( 2, leftPart.cols()+ lastcol );
        leftPart.rightCols( rightPart.cols()-1 ) = rightPart.block(0,1,2,lastcol);

        segmentation.addLine(leftPart);
//        gsInfo<<"\n left part \n "<< leftPart.transpose() <<"\n\n"
//                <<"\n";
    }

    // third step: construct Coon's patches


    //Clean Memory
    delete u.first;
    delete u.second;
    delete pmy_template;
}


template<class T>
std::pair<gsFunction<T>* , gsFunction<T>* >
mapto(gsPlanarDomain<T> & pdomain, gsTemplate<T> & my_template, gsMatrix<T> &Dirichlet,int numRefine)
{
    //  GISMO_ASSERT( my_template.domain().numLoops()== pdomain.numLoops(), "Template has different topology.");

    //******* basis vector of all basis functions of all curves in pd
    std::vector < gsBasis<T> * > basis;
    gsConstantFunction<T> f_1(0.0,1);
    gsConstantFunction<T> f_2 (0.0,1);
    std::vector<gsFunction<T>*> bc_0;
    std::vector<gsFunction<T>*> bc_1;
    gsMatrix<T> param(1,1);
    param<<1;
    typename gsBSpline<T>::uPtr tmpl_loop = memory::convert_ptr<gsBSpline<T> >(my_template.loop(0).singleCurve());

    for(index_t v=0; v<=pdomain.numHoles(); v++)
    {
        typename gsBSpline<T>::uPtr this_loop = memory::convert_ptr<gsBSpline<T> >( pdomain.loop(v).singleCurve() );

        basis.push_back(this_loop->basis().clone().release());
        if(v==0)
        {
            /* for the square template we need :
             *  gsKnotVector<T>sqKnots(0,1,3,2,1,1);
             *  in the place of
             *  tmpl_loop->basis().knots()*/
            bc_0.push_back( new gsBSpline<T>( tmpl_loop->basis().knots(),  tmpl_loop->coefs().col(0))  );
            bc_1.push_back( new gsBSpline<T>( tmpl_loop->basis().knots(),  tmpl_loop->coefs().col(1))  );
        }
        else if(v==1)
        {
            f_1.setValue(Dirichlet(v-1,0),1);
            f_2.setValue(Dirichlet(v,0),1);

            // **** RIGHT THINGS :
            bc_0.push_back(f_1.clone().release());//-0.334477 CII
            bc_1.push_back(f_2.clone().release());
        }
        else if(v==2)
        {
            f_1.setValue(Dirichlet(v,0),1);
            f_2.setValue(Dirichlet(v+1,0),1);
            bc_0.push_back(f_1.clone().release());
            bc_1.push_back(f_2.clone().release());
            /***************************************************************
             * please, do not delete the commented part
//             if(pdomain.numLoops()==1)
//             { */
            //             std::vector<gsFunction<T>*> bc_map(2);
            //             bc_map[0]=bc_0[0];
            //             bc_map[1]=bc_1[0];
            //            gsMatrix<T> av = averageValue(bc_map, tmpl_loop->basis().knots().breaks() );
            // //             gsMatrix<T> av (2,1);
            // //             av<< 0.5,
            // //                  0;
            //            bc_0.push_back( new gsConstantFunction<T> (av(0,0), 1 ) );
            //             bc_1.push_back( new gsConstantFunction<T> (av(1,0), 1 ) );
            //             }
            //             else
            //             {
            //                 bc_0.push_back( new gsConstantFunction<T> (*my_template.skeleton(0)->eval(param), 1 ) );
            //                 bc_1.push_back( new gsConstantFunction<T> (*my_template.skeleton(1)->eval(param), 1 ) );
            //             }
            //         }


            //             std::vector<gsFunction<T>*> bc_map(2);
            //             bc_map[0]=bc_0[0];
            //             bc_map[1]=bc_1[0];
            //             gsMatrix<T> av = averageValue(bc_map, tmpl_loop->basis().knots().breaks() );
            //             bc_0.push_back( new gsConstantFunction<T> (av(0,0), 1 ) );
            //             bc_1.push_back( new gsConstantFunction<T> (av(1,0), 1 )  );

            //***************************************/
        }

        if(v==3)//in case of 3 holes
        {
            f_1.setValue(Dirichlet(v+1,0),1);
            f_2.setValue(Dirichlet(v+2,0),1);

            bc_0.push_back(f_1.clone().release());
            bc_1.push_back(f_2.clone().release());

        }
    }


    for(index_t i=0; i<pdomain.numLoops() ; i++)
    {
        for (int k = 0; k < numRefine; ++k)
            basis[i]->uniformRefine();
    }

    int count = 0;
    for(index_t i=0; i< pdomain.numLoops() ; i++)
        count += basis[i]->size();

    // gsInfo<<"Number of DoFs: "<< count <<"\n";

    gsBemLaplace<T> Lsolver_1( &pdomain, bc_0, basis );

    gsBemLaplace<T> Lsolver_2( &pdomain, bc_1, basis );

    //////////////////********************
    gsBemSolution<T> * sol1= Lsolver_1.solve(true);
    gsBemSolution<T> * sol2= Lsolver_2.solve(true);
    ////////////////////*******************


    // Free the bases
    freeAll(basis);
    return std::make_pair(sol1, sol2);


}


int main(int argc, char *argv[])
{
    gsPlanarDomain<> * Pdomain = NULL;
    index_t n_points(10);
    real_t tolerance = 1e-4;
    index_t smpl(30);
    bool plot = false;
    std::vector<real_t> start_v;
    std::string fn = "planar/amoeba0_pdomain.xml";
    
    gsCmdLine cmd("Segmentation in quadrangular patches of a planar domain ");
    cmd.addMultiReal("i", "starting_point",
                     "point inside the computational domain from which looking for the starting point of the tracing curves",
                     start_v);
    cmd.addString("g","geometry","File containing Geometry (.axl, .txt)", fn);
    cmd.addSwitch("plot", "Plot result with ParaView ", plot);
    cmd.addInt("s","samples", "Number of samples", smpl);
    cmd.addInt("p","n_points", "Number of points traced per curve", n_points);
    cmd.addReal("t","tolerance", "Required accuracy", tolerance);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Pdomain =  gsReadFile<>( fn ) ;
    gsFileData<>  filedata(fn);
        
    if( filedata.has<gsPlanarDomain<> >() )
        Pdomain = filedata.getFirst< gsPlanarDomain<> >().release();
        
    else
    {
            
        if(filedata.has<gsCurve <> >() )
        {
            gsCurve<> * geo = filedata.getAnyFirst< gsCurve<> >().release();
            gsCurveLoop<> * cp = new gsCurveLoop<>( geo );
            Pdomain = new gsPlanarDomain<>( cp );
        }
            
            
        else
        {
            gsWarn<< "Did not find any planar domain or geometry in "<< fn<<", quitting.\n";
            return 1;
        }
    }
    
    if(start_v.empty())
    {
        start_v.push_back(1.5);
        start_v.push_back(2);
    }
    
    gsInfo << "------------------------------------------------------------"
        "\nInput Arguments: \n\n"
        "input geometry: " << fn <<
        "\n how many samples:  "<<smpl<<
        "\n number of points traced per curve: "<<2*n_points<<
        "\n tolerance: "<<tolerance<<"\n"
        "------------------------------------------------------------"
        "\n\n";


//    gsAsMatrix<> start (start_v,2,1);
    gsMatrix<> indice;
    real_t l = 0.12;
    real_t area;
    real_t scalar = 1.0; //gives you te original domain
    int n= 4000;
    quadPointsMultConnected(Pdomain,n,l,scalar,indice,area);


    bool c= true;
    gsMesh<> res;
    gsTemplate<real_t> Circle(2);

    gsMatrix<real_t> Dirichlet(4,1);// to be changed accordingly to #holes

    //HAT:
    Dirichlet<<    0.10,
        -0.55,
        0.25,
        0.55;
    
    
    segment<real_t>(*Pdomain,Dirichlet, indice, n_points,tolerance,c,res);
    
    //  gsInfo<<"\n res \n "<<res << "\n";
    
    int exitCommand = 0;
    if(plot)
    {

        gsWriteParaview<real_t>( *Pdomain, "Hat_pdomain", smpl) ;
        gsWriteParaview<real_t>(Circle.domain(),"CircleTemplate",smpl);
        gsWriteParaviewPoints<real_t>(indice,"sampledPdomain");
        gsWriteParaview<real_t>( res, "Segmentation60and90HAT");

        //run: paraview
        exitCommand = system("paraview Hat_pdomain.pvd &");
    }

    delete Pdomain;

    return exitCommand;
}
