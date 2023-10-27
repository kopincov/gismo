/** @file gsHTensorBasis.h

    @brief Provides implementation of HTensorBasis common operations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss
*/

#include <iostream>
#include <set>
#include <map>

#include <gismo.h>
#include <gismo_dev.h>

#include "common/random_refinement.h"


using namespace gismo;

int main()
{

    std::string bspl_filename;
    //bspl_filename= "face.xml"; //default example coarsening_test.xml coarstest
    bspl_filename= "coarstest.xml";
    //bspl_filename= "coarsening_test.xml";
    //bspl_filename= "coarsening_test_1.xml";

    memory::unique_ptr< gsGeometry<> > bspl;
    if ( ! bspl_filename.empty() )
    {
        gsFileData<>  data( bspl_filename );
        if ( data.has< gsTensorBSpline<2,real_t>  >() )
        {
          bspl = data.getFirst< gsTensorBSpline<2,real_t>  >();
        }else
        {
             gsInfo<<"wrong file." << "\n";
            return 1;
        }
    }
    int num_of_iter = 1;
    for(int ii = 0; ii < num_of_iter;ii++){
        gsInfo<<"/////////////////////////////////Face test HB///////////////////////////"<<"\n";
        gsInfo<<"Test number: "<<ii<<"\n";
//        gsHBSpline<2> * hbs = 0;
//        std::string thb_filename;

//        //thb_filename= "thbs_face_3levels.xml";
//        thb_filename= "hbs_01.xml";
//        gsFileData<>  data( thb_filename );
//        if ( data.has< gsHBSpline<2> >() )
//        {
//            hbs = data.getFirst< gsHBSpline<2> >();
//        }else{
//            gsInfo<<"wrong imput file"<<"\n";
//            return 0;
//        }

        gsHBSplineBasis<2> tt(bspl->basis());
        //gsHBSplineBasis<2> tt(hbs->basis());
        //std::vector<gsSortedVector<unsigned> > OX = tt.getXmatrix();

//       std::vector<unsigned>boxes = random_refinement(3,3,&tt);

//for face
//        std::vector<unsigned int> boxes;
//        boxes.push_back(1);
//        boxes.push_back(1);
//        boxes.push_back(12);
//        boxes.push_back(42);
//        boxes.push_back(28);

//        boxes.push_back(2);
//        boxes.push_back(16);
//        boxes.push_back(21);
//        boxes.push_back(25);
//        boxes.push_back(27);

//        boxes.push_back(3);
//        boxes.push_back(16);
//        boxes.push_back(22);
//        boxes.push_back(20);
//        boxes.push_back(27);

//        boxes.push_back(2);
//        boxes.push_back(6);
//        boxes.push_back(2);
//        boxes.push_back(16);
//        boxes.push_back(18);



        std::vector<index_t> boxes;
        boxes.push_back(1);
        boxes.push_back(1);
        boxes.push_back(5);
        boxes.push_back(12);
        boxes.push_back(17);

//        boxes.push_back(2);
//        boxes.push_back(7);
//        boxes.push_back(12);
//        boxes.push_back(8);
//        boxes.push_back(14);

        boxes.push_back(2);
        boxes.push_back(4);
        boxes.push_back(8);
        boxes.push_back(11);
        boxes.push_back(12);

        boxes.push_back(3);
        boxes.push_back(6);
        boxes.push_back(9);
        boxes.push_back(11);
        boxes.push_back(11);







//        boxes.push_back(2);
//        boxes.push_back(10);
//        boxes.push_back(20);
//        boxes.push_back(20);
//        boxes.push_back(31);

//        boxes.push_back(2);
//        boxes.push_back(5);
//        boxes.push_back(18);
//        boxes.push_back(15);
//        boxes.push_back(47);

//        boxes.push_back(3);
//        boxes.push_back(21);
//        boxes.push_back(28);
//        boxes.push_back(42);
//        boxes.push_back(31);

//        boxes.push_back(3);
//        boxes.push_back(11);
//        boxes.push_back(25);
//        boxes.push_back(23);
//        boxes.push_back(26);

//        boxes.push_back(3);
//        boxes.push_back(23);
//        boxes.push_back(9);
//        boxes.push_back(31);
//        boxes.push_back(35);

//        boxes.push_back(4);
//        boxes.push_back(39);
//        boxes.push_back(5);
//        boxes.push_back(74);
//        boxes.push_back(7);

//        boxes.push_back(5);
//        boxes.push_back(75);
//        boxes.push_back(253);
//        boxes.push_back(183);
//        boxes.push_back(308);

//        boxes.push_back(5);
//        boxes.push_back(14);
//        boxes.push_back(95);
//        boxes.push_back(73);
//        boxes.push_back(162);

//        boxes.push_back(5);
//        boxes.push_back(163);
//        boxes.push_back(27);
//        boxes.push_back(190);
//        boxes.push_back(61);
//        std::vector<unsigned int> boxes;
//        boxes.push_back(1);
//        boxes.push_back(2);
//        boxes.push_back(2);
//        boxes.push_back(3);
//        boxes.push_back(4);

//        boxes.push_back(1);
//        boxes.push_back(0);
//        boxes.push_back(18);
//        boxes.push_back(11);
//        boxes.push_back(24);

//        boxes.push_back(1);
//        boxes.push_back(7);
//        boxes.push_back(6);
//        boxes.push_back(10);
//        boxes.push_back(20);

//        boxes.push_back(3);
//        boxes.push_back(10);
//        boxes.push_back(6);
//        boxes.push_back(15);
//        boxes.push_back(92);


//        tt.refineElements(boxes);

        //tt.increaseMultiplicity(1,0,0.5);
        //tt.increaseMultiplicity(2,1,0.75);

//        gsFileData<> newdata;
//        newdata << tt ;
//        newdata.dump("HB_transfer_test");

        //compute transfer matrices
//        std::vector<gsMatrix<> > tr;
//        tt.transferbyLvl(tr);//transfer(TT);
//        gsInfo << "transfers ready"  << "\n";

//        gsMatrix<> coefs = bspl->coefs();
//        for(unsigned int i = 0; i < tr.size();i++)
//        {
//            coefs = tr[i]*coefs;

//        }
//        gsHBSpline<2> H (tt, coefs);
//        gsMatrix<> paramRange = H.parameterRange();
//        const real_t dist = computeMaximumDistance<real_t>(*bspl, H, paramRange.col(0), paramRange.col(1));
//        gsInfo << "Maximum distance: " << dist << "\n";
//        if (dist > 1e-12) {  gsInfo << "FAIL" << "\n"; break; }

//        gsInfo << "Direct transfer test for HB" << "\n";

//        //gsMatrix<> dir_transf = coarsening(OX,NX,transfer);
//        gsMatrix<> dir_transf;
//        tt.transfer(OX, dir_transf);
//        gsInfo << "transfers 2 ready"  << "\n";
//        gsMatrix<> coefs2 = bspl->coefs();
//        coefs2 = dir_transf*coefs2;
//        gsHBSpline<2> H2 (tt, coefs2);
//        gsMatrix<> paramRange = H2.parameterRange();
//        const real_t dist2 = computeMaximumDistance<real_t>(*bspl, H2, paramRange.col(0), paramRange.col(1));
//        gsInfo << "Maximum distance: " << dist2 << "\n";
//        if (dist2 > 1e-12) {  gsInfo << "FAIL" << "\n"; break;}

//        gsMatrix<> mult=  tr[1]*tr[0];
//        gsInfo<<"Differences: "<<"\n";
//        for(int i =0; i <  mult.rows();i++ ){
//            for(int j = 0; j < mult.cols();j++){
//                if(mult(i,j)!=dir_transf(i,j)){
//                    std:: gsInfo<<"mult = "<< mult(i,j)<<" direct= "<<dir_transf(i,j)<<" col= "<<j<<" row= "<<i<<"\n";
//                }
//            }
//        }
//        gsInfo<<mult<<"\n"<<"\n";
//        gsInfo<<dir_transf<<"\n";

        //transfer with coef test
//        gsInfo << "withCoef transfer test for HB" << "\n";
//        gsHBSplineBasis<2> tt_withCoef(bspl->basis());
//        gsMatrix<> coefs_withCoef = bspl->coefs();
//        tt_withCoef.refineElements_withCoefs(coefs_withCoef,boxes);
//        gsHBSpline<2> H_withCoef (tt_withCoef, coefs_withCoef);
//        const real_t dist_withCoef = computeMaximumDistance<real_t>(*bspl, H_withCoef, paramRange.col(0), paramRange.col(1));
//        gsInfo << "Maximum distance: " << dist_withCoef << "\n";
//        //gsInfo<<coefs_withCoef<<"\n";
//        //gsInfo<<tt_withCoef<<"\n";
//        if (dist_withCoef > 1e-12) {  gsInfo << "FAIL" << "\n"; break;}

//        gsInfo << "withCoef transfer2 test for HB" << "\n";
//        gsHBSplineBasis<2> tt_withCoef2(bspl->basis());
//        gsMatrix<> coefs_withCoef2 = bspl->coefs();
//        tt_withCoef2.refineElements_withCoefs2(coefs_withCoef2,boxes);
//        gsHBSpline<2> H_withCoef2 (tt_withCoef2, coefs_withCoef2);
//        const real_t dist_withCoef2 = computeMaximumDistance<real_t>(*bspl, H_withCoef2, paramRange.col(0), paramRange.col(1));
//        gsInfo << "Maximum distance: " << dist_withCoef2 << "\n";
//        //gsInfo<<coefs_withCoef2<<"\n";
//        //gsInfo<<tt_withCoef2<<"\n";
//        if (dist_withCoef2 > 1e-12) {  gsInfo << "FAIL" << "\n"; break;}
//        for (int ii = 0; ii < coefs_withCoef.rows(); ii++){
//            for (int jj = 0; jj < coefs_withCoef.cols(); jj++){
//                if(math::abs(coefs_withCoef(ii,jj)-coefs_withCoef2(ii,jj)) > 0.00001){
//                gsInfo<<coefs_withCoef<<"\n"<<"\n";
//                gsInfo<<coefs_withCoef2<<"\n"<<"\n";
//                break;
//                }
//            }
//        }




//        gsInfo<<"/////////////////////////////////Face test THB///////////////////////////"<<"\n";
//        gsTHBSpline<2> * hbs = 0;
//        std::string thb_filename;
//        //thb_filename= "thbs_face_3levels.xml";
//        thb_filename= "thbs_01.xml";
//        gsFileData<>  data( thb_filename );
//        if ( data.has< gsTHBSpline<2> >() )
//        {
//            hbs = data.getFirst< gsTHBSpline<2> >();
//        }else{
//            gsInfo<<"wrong imput file"<<"\n";
//            return 0;
//        }
        gsTHBSplineBasis<2> tt_thb(bspl->basis());
        //gsTHBSplineBasis<2> tt_thb = hbs->basis();

        //T_thb->refine(boxes);//refining the basis
        std::vector<gsSortedVector<index_t> > OX_thb = tt_thb.getXmatrix();

        tt_thb.refineElements(boxes);
        gsInfo<<"basis size "<<tt_thb.size()<<"\n";
        std::vector<gsSortedVector<index_t> > NX_thb = tt_thb.getXmatrix();
        //compute transfer matrices
        std::vector<gsSparseMatrix<> > tr_thb;
        tt_thb.transferbyLvl(tr_thb);//transfer_thb(T_thb);
        gsInfo << "transfers ready"  << "\n";

        gsMatrix<> coefs_thb = bspl->coefs();


        //gsMatrix<> coefs_thb = hbs->coefs();
        for(unsigned int i = 0; i < tr_thb.size();i++)
        {
            //gsInfo<<tr_thb[i].rows()<<" "<<tr_thb[i].cols()<<"    "<<bspl->basis().size()<<"\n";
            coefs_thb = tr_thb[i]*coefs_thb;

        }
        //gsInfo<<"b"<<"\n";
        gsTHBSpline<2> H_thb (tt_thb, coefs_thb);
        //gsInfo<<"c"<<"\n";
        gsMatrix<> paramRange_thb = H_thb.parameterRange();
        //gsInfo<<"d"<<"\n";
        const real_t dist_thb = computeMaximumDistance<real_t>(*bspl, H_thb, paramRange_thb.col(0), paramRange_thb.col(1));
        //gsInfo<<"e"<<"\n";
        gsInfo << "Maximum distance: " << dist_thb << "\n";
        if (dist_thb > 1e-12) {  gsInfo << "FAIL" << "\n"; break;}

        gsSparseMatrix<> dir_transf_thb;
        gsStopwatch time;

        tt_thb.transfer(OX_thb, dir_transf_thb);//= coarsening_thb(OX_thb,NX_thb,transfer_thb);
        const double transftime = time.stop();
        gsInfo  << "Time of computation in milliseconds: "<<"transf time "<<transftime <<"\n";
        gsInfo << "transfers 2 ready"  << "\n";
        gsMatrix<> coefs2_thb = bspl->coefs();

        //gsMatrix<> coefs2_thb = hbs->coefs();

        coefs2_thb = dir_transf_thb*coefs2_thb;
        gsTHBSpline<2> H2_thb (tt_thb, coefs2_thb);
                //gsMatrix<> paramRange_thb = H2_thb.parameterRange();
        const real_t dist2_thb = computeMaximumDistance<real_t>(*bspl, H2_thb, paramRange_thb.col(0), paramRange_thb.col(1));
        gsInfo << "Maximum distance: " << dist2_thb << "\n";
        //gsWriteParaview( H2_thb, "THB_test_thb", 500);
        //gsWriteParaview( tt_thb, "THB_test_basis", 500);
        //gsWriteParaview( *hbs, "THB_test_orig", 500);

//        gsFileData<> newdata;
//        newdata << tt_thb ;
//        newdata.dump("THB_transfer_test");
        if (dist2_thb > 1e-12) {  gsInfo << "FAIL" << "\n"; break;}

        gsInfo << "withCoef transfer test for THB" << "\n";
        gsTHBSplineBasis<2> tt_withCoef_thb(bspl->basis());
        gsMatrix<> coefs_withCoef_thb = bspl->coefs();
        //gsInfo<<"a"<<"\n";
        tt_withCoef_thb.refineElements_withCoefs(coefs_withCoef_thb,boxes);
        //gsInfo<<"b"<<"\n";
        gsTHBSpline<2> H_withCoef_thb (tt_withCoef_thb, coefs_withCoef_thb);
        //gsInfo<<"c"<<"\n";
        const real_t dist_withCoef_thb = computeMaximumDistance<real_t>(*bspl, H_withCoef_thb, paramRange_thb.col(0), paramRange_thb.col(1));
        gsInfo << "Maximum distance: " << dist_withCoef_thb << "\n";
        gsInfo<<coefs_withCoef_thb<<"\n";
        if (dist_withCoef_thb > 1e-12) {  gsInfo << "FAIL" << "\n"; break;}


        gsInfo << "withCoef transfer2 test for THB" << "\n";
        gsTHBSplineBasis<2> tt_withCoef_thb2(bspl->basis());
        gsMatrix<> coefs_withCoef_thb2 = bspl->coefs();
        //gsInfo<<"a"<<"\n";
        tt_withCoef_thb2.refineElements_withCoefs2(coefs_withCoef_thb2,boxes);
        //gsInfo<<"b"<<"\n";
        gsTHBSpline<2> H_withCoef_thb2 (tt_withCoef_thb2, coefs_withCoef_thb2);
        //gsInfo<<"c"<<"\n";
        const real_t dist_withCoef_thb2 = computeMaximumDistance<real_t>(*bspl, H_withCoef_thb2, paramRange_thb.col(0), paramRange_thb.col(1));
        gsInfo << "Maximum distance: " << dist_withCoef_thb2 << "\n";
        gsInfo<<coefs_withCoef_thb2<<"\n";

        for (int ii2 = 0; ii2 < coefs_withCoef_thb.rows(); ii2++){
            for (int jj = 0; jj < coefs_withCoef_thb.cols(); jj++){
                if(math::abs(coefs_withCoef_thb(ii2,jj)-coefs_withCoef_thb2(ii2,jj)) > 0.00001){
                gsInfo<<coefs_withCoef_thb(ii2,jj)<<"\n"<<"\n";
                gsInfo<<coefs_withCoef_thb2(ii2,jj)<<"\n"<<"\n";
                gsInfo<<"ii "<<ii2<<" jj "<<jj<<"\n";
                }
            }
        }
        if (dist_withCoef_thb2 > 1e-12) {  gsInfo << "FAIL" << "\n"; break;}



//        gsInfo << "UniformRefine transfer test for THB" << "\n";
//        gsTHBSplineBasis<2> tt_uniform(tt_withCoef_thb);
//        gsMatrix<> coefs_uniform = coefs_withCoef_thb;
//        tt_uniform.uniformRefine_withCoefs(coefs_uniform);
//        gsInfo<<"uniform size"<<tt_uniform.size()<<"\n";
////        gsWriteParaview( tt_withCoef_thb, "THB_test_basis_withCoef", 5000);
////        gsWriteParaview( tt_uniform, "THB_test_basis_uniform", 5000 );
//        gsTHBSpline<2> H_uniform (tt_uniform, coefs_uniform);
//        const real_t dist_uniform = computeMaximumDistance<real_t>(H_withCoef_thb, H_uniform, paramRange_thb.col(0), paramRange_thb.col(1));
//        gsInfo << "Maximum distance: " << dist_uniform  << "\n";
//        if (dist_uniform > 1e-12) {  gsInfo << "FAIL" << "\n"; break;}

//        gsMatrix<> mult=  tr_thb[1]*tr_thb[0];
//        gsInfo<<"mult size "<< mult.rows()<<" "<< mult.cols()<<"... dir_tran size: "<<dir_transf_thb.rows()<<" "<<dir_transf_thb.cols()<<"\n";

//        gsInfo<<"Differences: "<<"\n";
//        for(int i =0; i <  mult.rows();i++ ){
//            for(int j = 0; j < mult.cols();j++){
//                if(mult(i,j)!=dir_transf_thb(i,j)){
//                    std:: gsInfo<<"mult = "<< mult(i,j)<<" direct= "<<dir_transf_thb(i,j)<<"col= "<<j<<" row= "<<i<<"\n";
//                }
//            }
//        }
       // if (dist_thb > 1e-12) {  gsInfo << "FAIL" << "\n"; break;}
//        gsWriteParaview( H2_thb, "THB_test_thb", 500);
//        gsWriteParaview( tt_withCoef_thb, "THB_test_basis_withCoef", 1000);
//        gsWriteParaview( tt_uniform, "THB_test_basis_uniform", 1000 );
//        gsFileData<> newdata;
//        newdata << tt_withCoef_thb ;
//        newdata.dump("THB_transfer_test_withCoef");

//        gsFileData<> newdata1;
//        newdata1 << tt_uniform ;
//        newdata1.dump("THB_transfer_test_uniform");
//        gsWriteParaview( *hbs, "THB_test_orig", 500);
//        gsWriteParaview( *bspl, "THB_test_bspl", 5000);
    }



   return 0;
}
