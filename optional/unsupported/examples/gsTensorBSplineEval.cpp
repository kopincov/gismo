// We are testing evaluation of tensor BSplines. We compute point on Tensor
// BSpline (algorithm uses knot insertion) and compare with normal evaluation
// algorithm.
//
// Author: Jaka Speh

#include <iostream>
#include <gismo.h>

#include <gsNurbs/gsDeboor.hpp>



using namespace gismo;


void example01(unsigned type);
void example02(unsigned type);


int main()
{

    example01(0);

    example02(0);

//    example01(1);

//    example02(1);

    return 0;

}

template <typename KnotVectorType>
void compareEvaluationMultipleTimes3D(KnotVectorType& KV_1,
                                      KnotVectorType& KV_2,
                                      KnotVectorType& KV_3,
                                      gsMatrix<>& coefs)
{


    gsTensorBSpline<3> bs = gsTensorBSpline<3>( KV_1, KV_2, KV_3, coefs);

    gsMatrix<> para  = bs.support();
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);

    gsStopwatch time;

    unsigned tmp[] = {10, 100, 1000, 10000, 100000, 1000000};
    std::vector<unsigned> vec(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));
    std::vector<double> times1(sizeof(tmp) / sizeof(tmp[0]), 0);
    std::vector<double> times2(sizeof(tmp) / sizeof(tmp[0]), 0);

    // 50 times we repeat comparison and measure time
    for (unsigned ii = 0; ii < 50; ii++)
    {
    for (unsigned i = 0; i < vec.size(); ++i)
    {

        gsMatrix<> params = uniformPointGrid(c0, c1, vec[i]);

        gsTensorBSplineBasis<3, real_t> base = bs.basis();

        gsMatrix<> deboor_result;

        time.restart();
        gsTensorDeboor_v2<3, real_t, gsKnotVector<>, gsMatrix<> >
                (params, base, coefs, deboor_result);
        double t1 = time.stop();


        gsMatrix<> result;

        time.restart();
        bs.eval_into(params, result);
        double t2 = time.stop();


        times1[i] += t1 / 50;
        times2[i] += t2 / 50;
    }
    }

    for (unsigned i = 0; i < vec.size(); ++i)
    {
        gsInfo << "======================================="
                  << "Evaluation with " << vec[i] << " points: \n"
                  << "Old evaluation: " << times2[i] << " s \n"
                  << "New evaluation: " << times1[i] << " s \n"
                  //<< "Difference in accuracy: " << diff.norm() << "\n"
                  << "\n";
    }

//    gsInfo << "[";
//    for (unsigned i = 0; i < vec.size(); ++i)
//    {
//        gsInfo << times2[i];
//        gsInfo << ((i == vec.size() - 1) ? "" : ", ");
//    }
//    gsInfo << "]" << "\n";


//    gsInfo << "[";
//    for (unsigned i = 0; i < vec.size(); ++i)
//    {
//        gsInfo << times1[i];
//        gsInfo << ((i == vec.size() - 1) ? "" : ", ");
//    }
//    gsInfo << "]" << "\n";

}

template <typename KnotVectorType>
void compareEvaluationMultipleTimes2D(KnotVectorType& KV_1,
                                      KnotVectorType& KV_2,
                                      gsMatrix<>& coefs)
{
    //    gsInfo<<"Knot Vector 1"<< KV_1<< "\n";
    //    gsInfo<<"Knot Vector 2"<< KV_2<< "\n";

    gsTensorBSpline<2> bs = gsTensorBSpline<2>( KV_1, KV_2, coefs);

    gsMatrix<> para  = bs.support();
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);

    gsStopwatch time;

    unsigned tmp[] = {10, 100, 1000, 10000, 100000, 1000000};
    std::vector<unsigned> vec(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));
    std::vector<double> times1(sizeof(tmp) / sizeof(tmp[0]), 0);
    std::vector<double> times2(sizeof(tmp) / sizeof(tmp[0]), 0);

    // 50 times we repeat comparison and measure time
    for (unsigned ii = 0; ii < 50; ii++)
    {
    for (unsigned i = 0; i < vec.size(); ++i)
    {

        gsMatrix<> params = uniformPointGrid(c0, c1, vec[i]);

        gsTensorBSplineBasis<2, real_t> base = bs.basis();

        gsMatrix<> deboor_result;

        time.restart();
        gsTensorDeboor_v2<2, real_t, gsKnotVector<>, gsMatrix<> >
                (params, base, coefs, deboor_result);
        double t1 = time.stop();


        gsMatrix<> result;

        time.restart();
        bs.eval_into(params, result);
        double t2 = time.stop();

        times1[i] += t1 / 50;
        times2[i] += t2 / 50;

    }
    }

    for (unsigned i = 0; i < vec.size(); ++i)
    {
        gsInfo << "======================================="
                  << "Evaluation with " << vec[i] << " points: \n"
                  << "Old evaluation: " << times2[i] << " s \n"
                  << "New evaluation: " << times1[i] << " s \n"
                  //<< "Difference in accuracy: " << diff.norm() << "\n"
                  << "\n";
    }

//    gsInfo << "[";
//    for (unsigned i = 0; i < vec.size(); ++i)
//    {
//        gsInfo << times2[i];
//        gsInfo << ((i == vec.size() - 1) ? "" : ", ");
//    }
//    gsInfo << "]" << "\n";


//    gsInfo << "[";
//    for (unsigned i = 0; i < vec.size(); ++i)
//    {
//        gsInfo << times1[i];
//        gsInfo << ((i == vec.size() - 1) ? "" : ", ");
//    }
//    gsInfo << "]" << "\n";


}

// compare evaluation in 2D with knot vectors KV_1 and KV_2, coefficients
// coefs and number of (uniformly sampled) points n
template <typename KnotVectorType>
real_t compareTwoEvaluation2D(KnotVectorType& KV_1, KnotVectorType& KV_2,
                            gsMatrix<>& coefs, unsigned n)
{
    gsInfo<<"Knot Vector 1"<< KV_1<< "\n";
    gsInfo<<"Knot Vector 2"<< KV_2<< "\n";

    gsTensorBSpline<2> bs = gsTensorBSpline<2>( KV_1, KV_2, coefs);

    gsMatrix<> para  = bs.support();
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    gsMatrix<> params = uniformPointGrid(c0, c1, n);

    gsTensorBSplineBasis<2, real_t> base = bs.basis();

    gsMatrix<> deboor_result;
    gsTensorDeboor_v2<2, real_t, gsKnotVector<>, gsMatrix<> >
            (params, base, coefs, deboor_result);


    gsMatrix<> result;
    bs.eval_into(params, result);

    gsMatrix<> diff = deboor_result - result;

    return diff.norm();
}

// compare evaluation in 3D with knot vectors KV_1, KV_2 and KV_3, coefficients
// coefs and number of (uniformly sampled) points n
template <typename KnotVectorType>
real_t compareTwoEvaluation3D(KnotVectorType& KV_1, KnotVectorType& KV_2,
                              KnotVectorType& KV_3, gsMatrix<>& coefs,
                              unsigned n)
{
    gsInfo<<"Knot Vector 1"<< KV_1<< "\n";
    gsInfo<<"Knot Vector 2"<< KV_2<< "\n";
    gsInfo<<"Knot Vector 3"<< KV_3<< "\n";

    gsTensorBSpline<3> bs = gsTensorBSpline<3>( KV_1, KV_2, KV_3, coefs);

    gsMatrix<> para  = bs.support();
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    gsMatrix<> params = uniformPointGrid(c0, c1, n);

    gsTensorBSplineBasis<3, real_t> base = bs.basis();

    gsMatrix<> deboor_result;
    gsTensorDeboor_v2<3, real_t, gsKnotVector<>, gsMatrix<> >
            (params, base, coefs, deboor_result);


    gsMatrix<> result;
    bs.eval_into(params, result);

    gsMatrix<> diff = deboor_result - result;

    return diff.norm();
}

// type 0 - just accuracy
// type 1 - also time measurement
void example01(unsigned type)
{

    real_t tmp[] = {0, 0, 0, 0.2, 0.2, 0.2, 0.6, 1, 1, 1};
    std::vector<real_t> knots(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));

    gsKnotVector<> KV_1(knots, 2);
    //gsKnotVector<> KV_1(0,1,4,3);
    gsKnotVector<> KV_2(knots, 2);
    //gsKnotVector<> KV_2(0,1,4,3);


    gsMatrix<> coefs(49,3);
    coefs.setZero();
    coefs.row(0) << 0.6, 0.08, 0.88;
    coefs.row(1) << 1.0, 0.64, 0.56;
    coefs.row(2) << 0.44, 0.48, 0.72;
    coefs.row(3) << 1.84, 0.28, 1.76;
    coefs.row(4) << 1.88, 1.6, 0.28;
    coefs.row(5) << 1.64, 0.76, 1.16;
    coefs.row(6) << 0.92, 1.04, 1.28;
    coefs.row(7) << 0.08, 1.76, 0.64;
    coefs.row(8) << 0.04, 1.52, 0.96;
    coefs.row(9) << 1.8, 0.12, 0.56;
    coefs.row(10) << 1.32, 0.36, 1.92;
    coefs.row(11) << 1.36, 0.4, 1.44;
    coefs.row(12) << 0.6, 1.12, 0.96;
    coefs.row(13) << 1.28, 1.36, 1.32;
    coefs.row(14) << 1.2, 0.84, 1.64;
    coefs.row(15) << 0.76, 0.6, 1.96;
    coefs.row(16) << 0.64, 1.12, 0.32;
    coefs.row(17) << 0.2, 0.04, 0.96;
    coefs.row(18) << 0.04, 0.72, 1.24;
    coefs.row(19) << 1.48, 0.56, 1.64;
    coefs.row(20) << 0.0, 1.24, 1.48;
    coefs.row(21) << 0.16, 0.96, 0.04;
    coefs.row(22) << 1.04, 1.16, 0.04;
    coefs.row(23) << 0.68, 0.36, 0.92;
    coefs.row(24) << 1.64, 0.88, 1.08;
    coefs.row(25) << 1.24, 2.0, 0.48;
    coefs.row(26) << 0.64, 0.64, 0.4;
    coefs.row(27) << 0.96, 1.36, 1.68;
    coefs.row(28) << 0.52, 0.12, 0.6;
    coefs.row(29) << 0.4, 1.32, 1.2;
    coefs.row(30) << 0.96, 0.84, 1.44;
    coefs.row(31) << 1.96, 0.12, 0.56;
    coefs.row(32) << 0.68, 0.32, 0.32;
    coefs.row(33) << 0.4, 0.44, 1.28;
    coefs.row(34) << 0.64, 1.04, 1.12;
    coefs.row(35) << 1.2, 0.28, 1.56;
    coefs.row(36) << 1.64, 1.08, 1.16;
    coefs.row(37) << 0.44, 0.6, 0.88;
    coefs.row(38) << 0.04, 0.32, 1.0;
    coefs.row(39) << 0.2, 0.48, 0.6;
    coefs.row(40) << 1.96, 0.92, 0.84;
    coefs.row(41) << 0.16, 1.36, 1.32;
    coefs.row(42) << 1.72, 0.6, 1.32;
    coefs.row(43) << 1.0, 1.2, 0.4;
    coefs.row(44) << 0.72, 1.52, 1.56;
    coefs.row(45) << 0.8, 1.0, 1.96;
    coefs.row(46) << 1.64, 1.68, 0.08;
    coefs.row(47) << 0.0, 1.88, 1.68;
    coefs.row(48) << 0.72, 1.64, 0.2;



    if (type == 0)
    {
        int n = 10000;
        real_t error = compareTwoEvaluation2D<gsKnotVector<> >
                (KV_1, KV_2, coefs, n);

        gsInfo << "Number of points is " << n << "." << "\n";
        gsInfo << "Difference between algorithms: " << error << "\n";

    }
    else
    {
        compareEvaluationMultipleTimes2D<gsKnotVector<> >
                (KV_1, KV_2, coefs);
    }




//    for (int i = 0; i < params.cols(); i++)
//    {
//        gsMatrix<> deboor_result1;
//        gsTensorDeboor_v2<2, real_t, gsKnotVector<>, gsMatrix<> >
//                (params.col(i), base, coefs, deboor_result1);

////        gsInfo << "result: " << deboor_result1.rows() << " x " << deboor_result1.cols()
////                  << "\n" << deboor_result1 << "\n";

//        gsMatrix<> real_result1;
//        bs.eval_into(params.col(i), real_result1);

////        gsInfo << "real result: " << real_result1.rows() << " x " << real_result1.cols()
////                  << "\n" << real_result1 << "\n";


//        gsMatrix<> diff = deboor_result1 - real_result1;
//        gsInfo << "difference norm: p=" << params.col(i).transpose()
//                  << " diff: " << diff.norm() << "\n";
//    }

}


// type 0 - just accuracy
// type 1 - also time measurement
void example02(unsigned type)
{

    real_t tmp[] = {0, 0, 0, 0.2, 0.2, 0.6, 1, 1, 1};
    std::vector<real_t> knots(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));

//    gsKnotVector<> KV_1(2, knots);
//    //gsKnotVector<> KV_1(0,1,4,3);
//    gsKnotVector<> KV_2(2, knots);
//    gsKnotVector<> KV_3(0,1,3,4);

    gsKnotVector<> KV_1(0, 1, 1, 5);
    //gsKnotVector<> KV_1(0,1,4,3);
    gsKnotVector<> KV_2(0, 1, 1, 5);
    gsKnotVector<> KV_3(0, 1, 3, 5);

    gsMatrix<> coefs(6 * 6 * 8, 3);

    coefs.row(0) << 0.16, 1.52, 0.76;
    coefs.row(1) << 2.0, 0.6, 1.16;
    coefs.row(2) << 0.24, 0.72, 1.96;
    coefs.row(3) << 0.68, 0.2, 1.4;
    coefs.row(4) << 0.8, 1.76, 1.96;
    coefs.row(5) << 1.76, 0.12, 0.4;
    coefs.row(6) << 0.6, 0.32, 0.48;
    coefs.row(7) << 2.0, 0.76, 1.2;
    coefs.row(8) << 1.12, 0.48, 1.12;
    coefs.row(9) << 0.44, 1.56, 0.64;
    coefs.row(10) << 1.36, 1.0, 0.88;
    coefs.row(11) << 1.8, 0.84, 1.2;
    coefs.row(12) << 1.88, 1.08, 1.8;
    coefs.row(13) << 1.44, 1.8, 0.2;
    coefs.row(14) << 0.08, 0.76, 0.84;
    coefs.row(15) << 0.68, 0.24, 0.64;
    coefs.row(16) << 1.48, 0.52, 1.48;
    coefs.row(17) << 1.84, 1.72, 0.56;
    coefs.row(18) << 1.76, 1.56, 1.8;
    coefs.row(19) << 0.2, 0.72, 0.96;
    coefs.row(20) << 0.16, 1.64, 0.6;
    coefs.row(21) << 0.56, 1.68, 1.0;
    coefs.row(22) << 0.24, 0.8, 1.2;
    coefs.row(23) << 0.68, 1.36, 1.92;
    coefs.row(24) << 0.52, 1.44, 0.6;
    coefs.row(25) << 0.12, 0.28, 1.16;
    coefs.row(26) << 0.68, 1.72, 1.8;
    coefs.row(27) << 1.6, 1.52, 1.88;
    coefs.row(28) << 1.88, 1.96, 0.64;
    coefs.row(29) << 1.96, 0.48, 0.96;
    coefs.row(30) << 1.8, 1.6, 0.2;
    coefs.row(31) << 0.92, 1.88, 1.72;
    coefs.row(32) << 1.96, 1.0, 0.0;
    coefs.row(33) << 0.44, 0.28, 0.72;
    coefs.row(34) << 0.28, 0.36, 0.56;
    coefs.row(35) << 1.2, 1.76, 1.2;
    coefs.row(36) << 0.28, 1.92, 1.36;
    coefs.row(37) << 1.32, 1.76, 1.4;
    coefs.row(38) << 0.08, 0.96, 0.36;
    coefs.row(39) << 1.0, 1.48, 0.92;
    coefs.row(40) << 0.4, 0.2, 1.48;
    coefs.row(41) << 1.08, 0.12, 0.96;
    coefs.row(42) << 1.2, 2.0, 0.36;
    coefs.row(43) << 1.88, 1.8, 1.68;
    coefs.row(44) << 0.0, 1.68, 2.0;
    coefs.row(45) << 0.32, 0.52, 0.48;
    coefs.row(46) << 1.2, 1.88, 0.68;
    coefs.row(47) << 1.92, 1.52, 0.28;
    coefs.row(48) << 1.16, 0.96, 0.68;
    coefs.row(49) << 1.0, 0.24, 1.64;
    coefs.row(50) << 0.24, 0.0, 0.8;
    coefs.row(51) << 1.44, 1.04, 1.4;
    coefs.row(52) << 1.48, 1.24, 0.6;
    coefs.row(53) << 1.56, 0.12, 0.92;
    coefs.row(54) << 1.16, 1.36, 1.56;
    coefs.row(55) << 0.12, 0.8, 0.24;
    coefs.row(56) << 0.48, 0.4, 0.8;
    coefs.row(57) << 0.04, 0.16, 1.56;
    coefs.row(58) << 1.48, 1.84, 1.08;
    coefs.row(59) << 1.0, 0.16, 0.28;
    coefs.row(60) << 0.4, 1.12, 1.32;
    coefs.row(61) << 0.92, 0.4, 1.16;
    coefs.row(62) << 1.2, 0.56, 1.84;
    coefs.row(63) << 0.28, 1.76, 1.84;
    coefs.row(64) << 0.12, 0.0, 0.4;
    coefs.row(65) << 1.12, 0.76, 1.92;
    coefs.row(66) << 0.4, 1.04, 0.88;
    coefs.row(67) << 1.4, 1.16, 0.0;
    coefs.row(68) << 1.72, 1.2, 0.92;
    coefs.row(69) << 0.0, 1.36, 1.96;
    coefs.row(70) << 1.36, 0.52, 1.24;
    coefs.row(71) << 1.2, 1.28, 0.28;
    coefs.row(72) << 0.8, 1.64, 0.64;
    coefs.row(73) << 1.44, 1.24, 1.16;
    coefs.row(74) << 0.28, 0.4, 0.28;
    coefs.row(75) << 1.36, 1.48, 1.4;
    coefs.row(76) << 1.52, 0.36, 0.8;
    coefs.row(77) << 1.28, 1.56, 1.52;
    coefs.row(78) << 0.56, 1.16, 0.12;
    coefs.row(79) << 1.96, 1.52, 0.36;
    coefs.row(80) << 0.56, 1.72, 1.64;
    coefs.row(81) << 1.88, 1.6, 1.48;
    coefs.row(82) << 0.24, 0.32, 1.76;
    coefs.row(83) << 0.32, 0.88, 1.24;
    coefs.row(84) << 1.16, 1.52, 0.08;
    coefs.row(85) << 0.28, 1.8, 1.8;
    coefs.row(86) << 0.2, 0.04, 1.96;
    coefs.row(87) << 0.52, 0.04, 1.16;
    coefs.row(88) << 0.56, 1.96, 0.8;
    coefs.row(89) << 1.08, 1.72, 0.4;
    coefs.row(90) << 0.52, 1.2, 0.16;
    coefs.row(91) << 1.52, 1.64, 0.4;
    coefs.row(92) << 0.64, 1.6, 0.36;
    coefs.row(93) << 1.88, 0.68, 0.04;
    coefs.row(94) << 0.92, 1.76, 1.24;
    coefs.row(95) << 0.44, 1.56, 0.4;
    coefs.row(96) << 0.52, 1.32, 1.16;
    coefs.row(97) << 0.64, 1.2, 0.92;
    coefs.row(98) << 0.12, 0.04, 0.12;
    coefs.row(99) << 1.2, 0.8, 0.48;
    coefs.row(100) << 0.08, 0.28, 1.76;
    coefs.row(101) << 0.96, 0.44, 0.96;
    coefs.row(102) << 1.32, 1.48, 1.92;
    coefs.row(103) << 1.52, 1.76, 1.68;
    coefs.row(104) << 0.24, 1.32, 0.48;
    coefs.row(105) << 0.0, 0.08, 1.24;
    coefs.row(106) << 1.4, 0.24, 0.28;
    coefs.row(107) << 0.08, 0.92, 1.24;
    coefs.row(108) << 0.72, 1.72, 1.6;
    coefs.row(109) << 0.88, 0.32, 1.84;
    coefs.row(110) << 1.52, 1.36, 1.84;
    coefs.row(111) << 1.32, 0.2, 0.36;
    coefs.row(112) << 0.0, 0.08, 1.96;
    coefs.row(113) << 0.6, 0.72, 0.2;
    coefs.row(114) << 1.08, 1.72, 0.56;
    coefs.row(115) << 0.56, 0.64, 1.36;
    coefs.row(116) << 1.16, 1.8, 0.0;
    coefs.row(117) << 1.92, 1.16, 1.28;
    coefs.row(118) << 0.84, 0.64, 0.52;
    coefs.row(119) << 1.2, 0.52, 0.56;
    coefs.row(120) << 0.16, 0.32, 1.44;
    coefs.row(121) << 1.96, 1.6, 1.44;
    coefs.row(122) << 0.88, 1.52, 0.72;
    coefs.row(123) << 1.28, 0.64, 1.64;
    coefs.row(124) << 1.08, 1.96, 1.56;
    coefs.row(125) << 0.44, 1.12, 1.76;
    coefs.row(126) << 0.0, 0.52, 0.24;
    coefs.row(127) << 0.92, 1.32, 1.16;
    coefs.row(128) << 1.96, 1.28, 0.6;
    coefs.row(129) << 1.56, 0.92, 1.04;
    coefs.row(130) << 1.16, 1.32, 0.64;
    coefs.row(131) << 1.36, 1.12, 1.68;
    coefs.row(132) << 1.12, 1.96, 0.92;
    coefs.row(133) << 1.52, 1.76, 0.56;
    coefs.row(134) << 1.36, 1.92, 0.0;
    coefs.row(135) << 0.64, 0.56, 0.0;
    coefs.row(136) << 1.96, 0.96, 1.24;
    coefs.row(137) << 0.36, 1.68, 1.36;
    coefs.row(138) << 1.36, 1.72, 0.44;
    coefs.row(139) << 1.64, 1.72, 0.92;
    coefs.row(140) << 1.72, 1.92, 0.08;
    coefs.row(141) << 1.8, 0.44, 0.08;
    coefs.row(142) << 0.68, 0.68, 0.8;
    coefs.row(143) << 1.76, 1.96, 1.56;
    coefs.row(144) << 1.8, 1.76, 0.28;
    coefs.row(145) << 1.04, 1.48, 0.96;
    coefs.row(146) << 1.12, 1.48, 0.76;
    coefs.row(147) << 1.36, 0.48, 0.08;
    coefs.row(148) << 0.36, 0.48, 0.32;
    coefs.row(149) << 1.56, 0.48, 1.04;
    coefs.row(150) << 1.44, 1.92, 0.04;
    coefs.row(151) << 0.68, 0.88, 1.04;
    coefs.row(152) << 1.92, 1.68, 0.76;
    coefs.row(153) << 1.6, 0.8, 0.32;
    coefs.row(154) << 1.04, 0.16, 1.36;
    coefs.row(155) << 1.92, 1.92, 1.12;
    coefs.row(156) << 0.96, 1.56, 1.68;
    coefs.row(157) << 1.84, 1.04, 1.24;
    coefs.row(158) << 0.68, 0.56, 0.0;
    coefs.row(159) << 0.6, 1.56, 1.2;
    coefs.row(160) << 0.24, 0.36, 1.28;
    coefs.row(161) << 1.32, 1.72, 0.12;
    coefs.row(162) << 1.92, 1.04, 1.44;
    coefs.row(163) << 0.88, 0.84, 0.6;
    coefs.row(164) << 0.88, 1.88, 0.6;
    coefs.row(165) << 1.4, 1.32, 0.4;
    coefs.row(166) << 1.76, 0.56, 1.08;
    coefs.row(167) << 0.72, 0.28, 1.08;
    coefs.row(168) << 0.4, 0.96, 1.44;
    coefs.row(169) << 1.72, 0.4, 1.04;
    coefs.row(170) << 1.04, 1.84, 1.52;
    coefs.row(171) << 0.68, 1.28, 0.44;
    coefs.row(172) << 1.84, 0.84, 1.64;
    coefs.row(173) << 1.0, 1.32, 0.2;
    coefs.row(174) << 1.52, 1.0, 0.8;
    coefs.row(175) << 0.16, 0.4, 1.32;
    coefs.row(176) << 1.6, 0.88, 1.96;
    coefs.row(177) << 1.32, 1.68, 0.88;
    coefs.row(178) << 0.52, 0.36, 0.0;
    coefs.row(179) << 1.24, 1.92, 0.32;
    coefs.row(180) << 1.4, 0.68, 0.52;
    coefs.row(181) << 0.24, 1.36, 0.6;
    coefs.row(182) << 1.36, 1.96, 1.32;
    coefs.row(183) << 0.08, 1.64, 0.52;
    coefs.row(184) << 0.52, 0.12, 0.36;
    coefs.row(185) << 0.64, 0.4, 1.44;
    coefs.row(186) << 1.16, 1.12, 0.8;
    coefs.row(187) << 1.24, 1.68, 0.68;
    coefs.row(188) << 0.52, 1.32, 0.68;
    coefs.row(189) << 0.0, 0.28, 1.8;
    coefs.row(190) << 1.0, 1.12, 2.0;
    coefs.row(191) << 0.08, 0.52, 0.04;
    coefs.row(192) << 0.4, 0.32, 0.52;
    coefs.row(193) << 0.92, 1.56, 0.92;
    coefs.row(194) << 0.08, 1.28, 0.8;
    coefs.row(195) << 1.68, 1.92, 0.16;
    coefs.row(196) << 1.8, 1.4, 1.56;
    coefs.row(197) << 1.08, 1.0, 1.36;
    coefs.row(198) << 0.6, 0.0, 0.48;
    coefs.row(199) << 1.32, 1.64, 1.84;
    coefs.row(200) << 0.24, 1.44, 0.92;
    coefs.row(201) << 0.12, 0.28, 1.72;
    coefs.row(202) << 1.44, 1.48, 1.8;
    coefs.row(203) << 0.48, 0.92, 0.56;
    coefs.row(204) << 0.04, 0.76, 0.16;
    coefs.row(205) << 1.52, 0.8, 1.92;
    coefs.row(206) << 1.8, 1.92, 1.72;
    coefs.row(207) << 0.12, 0.48, 1.52;
    coefs.row(208) << 1.12, 1.92, 1.96;
    coefs.row(209) << 1.24, 1.32, 1.08;
    coefs.row(210) << 0.48, 1.76, 0.6;
    coefs.row(211) << 1.6, 1.0, 0.0;
    coefs.row(212) << 1.2, 1.48, 0.16;
    coefs.row(213) << 2.0, 0.96, 1.48;
    coefs.row(214) << 0.92, 1.24, 0.32;
    coefs.row(215) << 0.16, 1.56, 1.32;
    coefs.row(216) << 0.92, 1.8, 1.8;
    coefs.row(217) << 1.72, 0.96, 1.48;
    coefs.row(218) << 1.4, 1.52, 0.28;
    coefs.row(219) << 1.76, 0.12, 0.12;
    coefs.row(220) << 1.36, 1.8, 1.2;
    coefs.row(221) << 0.8, 0.36, 1.76;
    coefs.row(222) << 0.84, 0.96, 0.92;
    coefs.row(223) << 0.56, 1.6, 0.2;
    coefs.row(224) << 1.92, 0.36, 1.28;
    coefs.row(225) << 1.64, 0.36, 0.12;
    coefs.row(226) << 1.52, 1.8, 0.56;
    coefs.row(227) << 0.76, 0.48, 0.52;
    coefs.row(228) << 0.48, 1.56, 0.0;
    coefs.row(229) << 0.28, 0.0, 0.56;
    coefs.row(230) << 0.04, 0.4, 0.08;
    coefs.row(231) << 0.76, 1.36, 0.84;
    coefs.row(232) << 0.08, 0.48, 0.64;
    coefs.row(233) << 1.0, 0.0, 1.28;
    coefs.row(234) << 0.96, 1.48, 1.32;
    coefs.row(235) << 0.48, 0.52, 1.52;
    coefs.row(236) << 0.88, 1.04, 1.64;
    coefs.row(237) << 0.32, 1.92, 2.0;
    coefs.row(238) << 0.44, 0.48, 1.4;
    coefs.row(239) << 0.52, 1.08, 1.56;
    coefs.row(240) << 0.08, 0.04, 0.24;
    coefs.row(241) << 1.36, 1.24, 1.72;
    coefs.row(242) << 0.88, 0.88, 0.52;
    coefs.row(243) << 0.24, 1.64, 1.16;
    coefs.row(244) << 1.76, 0.2, 1.0;
    coefs.row(245) << 0.2, 0.24, 0.04;
    coefs.row(246) << 1.76, 1.16, 0.48;
    coefs.row(247) << 1.48, 1.6, 1.36;
    coefs.row(248) << 0.48, 1.52, 1.32;
    coefs.row(249) << 1.68, 2.0, 0.72;
    coefs.row(250) << 1.64, 1.64, 1.92;
    coefs.row(251) << 1.44, 0.68, 1.04;
    coefs.row(252) << 0.92, 1.48, 1.6;
    coefs.row(253) << 1.12, 1.68, 1.0;
    coefs.row(254) << 1.0, 1.72, 1.88;
    coefs.row(255) << 1.76, 0.76, 1.08;
    coefs.row(256) << 1.48, 0.88, 1.96;
    coefs.row(257) << 1.84, 0.96, 1.12;
    coefs.row(258) << 0.6, 1.04, 0.04;
    coefs.row(259) << 1.72, 1.04, 0.32;
    coefs.row(260) << 1.4, 0.2, 0.48;
    coefs.row(261) << 1.6, 1.24, 1.32;
    coefs.row(262) << 1.16, 0.56, 0.64;
    coefs.row(263) << 1.08, 1.04, 1.44;
    coefs.row(264) << 0.96, 1.76, 0.4;
    coefs.row(265) << 1.36, 0.76, 1.44;
    coefs.row(266) << 0.68, 1.52, 0.56;
    coefs.row(267) << 0.24, 1.52, 1.0;
    coefs.row(268) << 0.68, 1.32, 1.16;
    coefs.row(269) << 2.0, 0.88, 1.16;
    coefs.row(270) << 0.96, 1.48, 0.32;
    coefs.row(271) << 1.72, 0.24, 1.32;
    coefs.row(272) << 1.28, 0.24, 0.84;
    coefs.row(273) << 0.76, 0.88, 1.04;
    coefs.row(274) << 1.04, 0.32, 0.4;
    coefs.row(275) << 1.64, 0.32, 1.68;
    coefs.row(276) << 1.32, 1.72, 0.12;
    coefs.row(277) << 0.08, 1.72, 1.04;
    coefs.row(278) << 1.72, 1.6, 1.76;
    coefs.row(279) << 0.68, 1.72, 1.24;
    coefs.row(280) << 0.84, 0.44, 0.12;
    coefs.row(281) << 0.6, 1.0, 0.32;
    coefs.row(282) << 1.4, 2.0, 1.76;
    coefs.row(283) << 0.0, 0.52, 0.92;
    coefs.row(284) << 1.96, 1.32, 1.88;
    coefs.row(285) << 0.64, 0.28, 0.96;
    coefs.row(286) << 0.96, 1.92, 1.04;
    coefs.row(287) << 0.08, 1.0, 0.68;

    if (type == 0)
    {
        int n = 10000;
        real_t error = compareTwoEvaluation3D<gsKnotVector<> >
                (KV_1, KV_2, KV_3, coefs, n);

        gsInfo << "Number of points is " << n << "." << "\n";
        gsInfo << "Difference between algorithms: " << error << "\n";
    }
    else
    {
        compareEvaluationMultipleTimes3D(KV_1, KV_2, KV_3, coefs);
    }

}
