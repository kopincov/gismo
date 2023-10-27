/** @file gsTriangularBezierBasis_test

    @brief Testing the TriangularBezier(Basis) classes.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): G. Kiss
*/

#include <gismo.h>

#include <iostream>

#include <gsBezier/gsTriangularBezierBasis.h>
#include <gsBezier/gsTriangularBezier.h>

using namespace gismo;

int main(int argc, char *argv[])
{

    gsTriangularBezierBasis<2,real_t> triang(2);


    gsFileData<> newdata;
    newdata << triang ;
    newdata.dump("triang_basis");

    //gsWriteParaview( triang, "trbasis", 10000);

    gsVector<> param;
    param.resize(2);
    param[0] = 0;
    param[1] = 0;
    gsVector<> res;
    res = triang.getBaricentricCoordinate(param);
    gsInfo<<"The barycentric coordinates are:\n "<<res[0]<<", "<<res[1]<<", "<<res[2]<<"\n"<<"\n";

    std::vector<gsVector<unsigned int> > compos;
    triang.getCompositions(compos);
    for(unsigned int i = 0; i < compos.size(); i++){
        for(int j = 0; j < compos[i].size(); j++){
            gsInfo<<compos[i][j]<<"  ";
        }
         gsInfo<<"\n";
    }

    gsMatrix<> para(2,2);
    para(0,0) = 0.1;
    para(1,0) = 0.1;
    para(0,1) = 0.4;
    para(1,1) = 0.4;
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    gsMatrix<> pts = uniformPointGrid(c0,c1, 10) ;
//    gsMatrix<> pts(2,3);
//        pts(0,0) =0.5;
//        pts(1,0) =0.125;

//        pts(0,1) =0.125;
//        pts(1,1) =0.5;


//        pts(0,2) =0.3;
//        pts(1,2) =0.3;
//    gsMatrix<> pts(2,9);// = uniformPointGrid(c0,c1, 2) ;
//    pts(0,0) =0.1;
//    pts(1,0) =0.1;
//    pts(0,1) =0.1;
//    pts(1,1) =0.2;
//    pts(0,2) =0.1;
//    pts(1,2) =0.3;
//    pts(0,3) =0.1;
//    pts(1,3) =0.4;
//    pts(0,4) =0.1;
//    pts(1,4) =0.5;
//    pts(0,5) =0.1;
//    pts(1,5) =0.6;
//    pts(0,6) =0.1;
//    pts(1,6) =0.7;
//    pts(0,7) =0.1;
//    pts(1,7) =0.8;
//    pts(0,8) =0.1;
//    pts(1,8) =0.9;
    gsMatrix<> baric(3, pts.cols());
    for(int i = 0; i < pts.cols();i++){
        baric.col(i) = triang.getBaricentricCoordinate(pts.col(i));
    }
    for(int i = 0; i < baric.cols();i++){
        for(int j = 0; j < baric.rows();j++){
            if(baric(j,i)<0){
                baric.removeCol(i);
                pts.removeCol(i);
                i--;
                break;
            }
        }
    }

    gsInfo<<"////////////EVALUATION//////////////////"<<"\n";
    gsMatrix<> result;
    gsMatrix<> sum (1, pts.cols());
    sum.setZero();
    for(int i = 0; i< triang.size(); i++){
        triang.evalSingle_into(i,pts,result);
        sum += result;
    }
    gsInfo<<"size is: "<<triang.size()<<"\n";
    gsInfo<<"sum is: "<<sum<<"\n";

    triang.eval_into(pts,result);
    gsInfo<<"evaluation result: "<<"\n"<<result<<"\n";
    gsInfo<<"column sum of results: "<<"\n"<< result.colwise().sum()<<"\n";


    gsInfo<<"Baricentric coordinates of evaluated points:\n"<<baric<<"\n"<<"\n";
    gsMatrix<> result_baric;

    gsInfo<<"Result of evaluation based on baric. coordinated:"<<"\n";
    for(int j = 0; j < triang.size();j++){
        triang.evalSingle_into_baric(j,baric,result_baric);
        gsInfo<<result_baric<<"\n";

    }


    gsInfo<<"////////////DERIVATIVE//////////////////"<<"\n";
    int function = 0;
    gsMatrix<> deriv1;
    triang.derivSingle_into(function,pts,deriv1);
    gsInfo<<"First derivative of function "<<function<<" is: "<<"\n"<<deriv1<<"\n"<<"\n";


    ///numerical first derivative
    gsMatrix<real_t, 2> u(2,2) ;
    real_t eps= 0.000000001;//epsilon distance from the evaluated point
    gsMatrix<> test_der1(1,2);
    gsMatrix<> der1(2,pts.cols());
    //for each evaluated point compute the finite difference
    for (int i = 0; i < pts.cols();i++){
        u(0,0) = pts(0,i)+eps;
        u(1,0) = pts(1,i);
        u(0,1) = pts(0,i)-eps;
        u(1,1) = pts(1,i);
        gsMatrix<> x ;

        triang.evalSingle_into(function,u,x);
        test_der1.col(0) = (x.col(0) - x.col(1))/(2*eps);

        u(0,0) = pts(0,i);
        u(1,0) = pts(1,i)+eps;
        u(0,1) = pts(0,i);
        u(1,1) = pts(1,i)-eps;
        gsMatrix<> y;
        triang.evalSingle_into(function,u,y);
        test_der1.col(1) = (y.col(0) - y.col(1))/(2*eps);


        der1(0,i) = test_der1(0,0);
        der1(1,i) = test_der1(0,1);
    }
    //gsInfo<<"Numerical evaluation of first derivative for function "<< function<<" is: "<<"\n"<<der1<<"\n";

    gsMatrix<> allderivs;
    triang.deriv_into(pts,allderivs);
    //gsInfo<<"First derivatives of all functions "<<"\n"<<allderivs<<"\n";


    gsInfo<<"//////////// 2ND DERIVATIVE //////////////////"<<"\n";
    gsMatrix<> deriv2;
    triang.deriv2Single_into(function,pts,deriv2);
    //gsInfo<<"Second derivative of function "<<function<<" is: "<<"\n"<<deriv2<<"\n"<<"\n";

    gsMatrix<> deriv2all;
    triang.deriv2_into(pts,deriv2all);

   // gsInfo<<"Second derivatives of all functions "<<"\n"<<deriv2all<<"\n";

    gsInfo<<"Support of functions "<<"\n"<<triang.support(1)<<"\n";

    gsInfo<<"//////////// Creating geometry //////////////////"<<"\n";

    int pos = 0;
    gsMatrix<>coeff(triang.size(),3);

    srand((unsigned)time(NULL));
    for(int i = 0 ; i <= triang.degree();i++){
        for(int j = 0 ; j <= triang.degree()-i;j++){
            coeff(pos,0) = (1.0/(triang.degree()+1))*j;
            coeff(pos,1) = (1.0/(triang.degree()+1))*i;
            coeff(pos,2) = (rand()%100)/100.0;
            pos++;
        }
    }
    gsInfo<<"coefs"<<"\n"<<coeff<<"\n";

    gsTriangularBezier<2,real_t> Tr_geom(triang,coeff);

    gsMatrix<> eval_geom = Tr_geom.eval(pts);

//    gsMatrix<> eval_geom2;
//    Tr_geom.eval_into_test(pts,eval_geom2);

//    for(int i = 0; i < eval_geom.cols();i++){
//        for(int j = 0; j < eval_geom.rows();j++){
//            if(abs(eval_geom(j,i)-eval_geom2(j,i))>0.000001){
//                gsInfo<<"evaluation problem"<<"\n";
//            }
//        }
//    }

    //subdividing Tr_geom into 3 patches
    gsVector<> split;
    split.resize(2);
    split[0] = 0.25;
    split[1] = 0.25;
    gsMultiPatch<> r;

    Tr_geom.subdivide(split,r);
    //plotting ste patches and the original geometry
    gsMesh<real_t> mesh;
    gsMesh<real_t> mesh1;
    gsMesh<real_t> mesh2;
    gsMesh<real_t> mesh3;
    Tr_geom.toMesh(mesh,500);
    r.patch(0).toMesh(mesh1,500);
    r.patch(1).toMesh(mesh2,500);
    r.patch(2).toMesh(mesh3,500);
    gsWriteParaview( mesh , "mesh_output" );
    gsWriteParaview( mesh1 , "patch1_output" );
    gsWriteParaview( mesh2 , "patch2_output" );
    gsWriteParaview( mesh3 , "patch3_output" );

    //"traditional" way of plotting- creates small bumps on pondaries of the mesh- mesh is rectangular
    gsFileData<> newdata2;
    newdata2 << Tr_geom ;
    newdata2.dump("triang_geometry");

    gsWriteParaview( Tr_geom, "trGeometry", 1000);

  return 0;
}

