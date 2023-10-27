

#include <iostream>
#include <set>
#include <map>

#include <gismo.h>
#include <gismo_dev.h>

/// TODO this test is very sensitive to the choice of the numerical parameters
/// this make it difficult to have it working for different coefficient types
/// better to find a different test.

using namespace gismo;

int main(int argc, char *argv[])
{
    //unsigned np(1000);
    std::string filename = "thbs_16.xml";
    gsTHBSpline<2>::uPtr hbs;
    
    gsCmdLine cmd("Hi, give me a file (.xml) with some hierarchical splines.");
    cmd.addPlainString("filename", "File containing hierarchical splines (.xml).", filename);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // The file data
    gsFileData<>  data( filename );
    /*if ( data.has< gsHBSpline<2> >() )
    {
        gsInfo<<"The HB spline functions are not fully functional"<<"\n";
        return 0;
        //hbs = data.getFirst< gsHBSpline<2> >();
    }*/
    if ( data.has< gsTHBSpline<2> >() )
    {
        hbs = data.getFirst< gsTHBSpline<2> >();
    }
    gsInfo<< "  Got "<< *hbs << "\n";


    gsTHBSplineBasis<2>  & HB = hbs->basis();

    HB.printCharMatrix();

    gsInfo<<"\n"<<"Size of the basis "<<HB.size()<<"\n";
    gsMatrix<> para  = hbs->parameterRange();
    para(0,0)= 0.1;
    para(1,0)= 0.1;
    para(0,1)= 0.9;
    para(1,1)= 0.9;
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    gsMatrix<> pts = uniformPointGrid(c0,c1, 8) ;

    gsInfo<<"-------------------------------------------------------------------------------------------------------------"<<"\n";
    gsInfo<<"Comparing the derivatives of the THB surface with the \"numerical derivation\"(Finite Difference)"<<"\n";
    gsInfo<<"-------------------------------------------------------------------------------------------------------------"<<"\n";
    //test of the first derivative

   /*pts->resize(2,3);
    (*pts)(0,0) = 0.4;
    (*pts)(1,0) = 0.5;


    (*pts)(0,1) = 0.5;
    (*pts)(1,1) = 0.5;


    (*pts)(0,2) = 0.6;
    (*pts)(1,2) = 0.5;*/

    gsInfo<<"points: "<<"\n"<<pts<<"\n";

    // Compute first derivatives by a finite difference implementation and
    // by the gsTHBSplineBasis implementation and confront the results
    gsMatrix<> der1, hbs_der1;
    hbs->gsFunction<>::deriv_into(pts,der1); // finite differences
    hbs->deriv_into(pts,hbs_der1);           // implementation in gsTHBSplineBasis

    real_t deriv_max_diff=(der1-hbs_der1).lpNorm<Eigen::Infinity>();
    bool ok1 = deriv_max_diff<0.000001;
    if (ok1)
        gsInfo<<"The comparison between the finite difference method and the first derivative was successfull"<<"\n";
    else
    {
        for (index_t d=0;d<2;++d)
            for (index_t p=0;p<pts.cols();++p)
                if (math::abs(der1(d,p)-(hbs_der1)(d,p))>0.000001)
                    gsInfo<<"something wrong with the partial derivative "<< d<< " on point "<< p<<"\n";
    }

  ///////////////////////////////second derivatives test///////////////////////

  gsMatrix<> der2, hbs_der2;
  hbs->gsFunction<>::deriv2_into(pts,der2);
  hbs->deriv2_into(pts,hbs_der2);
  // gsInfo<<"Second derivative by finit differenc method "<<"\n"<<der2<<"\n";
  // gsInfo<<"Second derivative by THB spline derivative"  <<"\n"<<hbs_der2<<"\n";

  real_t deriv2_max_diff=(der2-hbs_der2).lpNorm<Eigen::Infinity>();
  bool ok2 = deriv2_max_diff<0.001;
  if (ok2)
       gsInfo<<"second derivatives are ok :)"<<"\n";
  else
  {
      for (index_t d=0;d<der2.rows();++d)
          for (index_t p=0;p<pts.cols();++p)
              if (math::abs( der2(d,p)-hbs_der2(d,p) )>0.001)
                  gsInfo<<"something wrong with the derivative in position "<< d << 
                      " on point "<< p <<",  "<< der2(d,p) <<" ~ "<< hbs_der2(d,p) <<"\n";
  }

  return (ok1 && ok2) ? 0 : 1;
}
















