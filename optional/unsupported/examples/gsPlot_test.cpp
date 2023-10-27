// example viewer

#include <iostream>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunctionExpr.h>

#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsNurbsCreator.h>

#include <gsHSplines/gsHBSplineBasis.h>

#include <gsIO/gsPlot.h>



using namespace gismo;

int main(int argc, char *argv[])
{  

    gsPlot<> plot ;

    // gsFunctionExpr<> g("(7*x+2*x^2-3*x^3)*cos(x)*sin(x)") ;
    // gsInfo<<" Plotting function g(x)="<< g <<"\n";
    // plot.function( g );
    
    // int p = 2;
    // gsKnotVector<> KV (0,1,p,p+1) ;
    // gsBSplineBasis<> bsp(KV,p);
    // gsInfo<<" Plotting "<< bsp <<"\n";
    // plot.basis( bsp );
    // gsTensorBasis<2,gsBSplineBasis<> > tbsp( bsp.clone() , bsp.clone() );
    // gsInfo<<" Plotting "<< tbsp <<"\n";
    // plot.basis( tbsp );
    
    // gsTensorNurbs<real_t,2> tgeo = * NurbsQuarterAnnulus<>();
    // gsInfo<<" Plotting "<< tgeo <<"\n";
    // plot.geometry( tgeo );


    // gsFunctionExpr<> f("x+y+z") ;
    // gsField<> fld  ( & tgeo , f , false ) ;
    //plot.field( fld );

    gsKnotVector<> KV (0, 1, 3,3, 1 ) ;
    gsBSplineBasis<> bsp( KV, 2 );
    gsTensorBasis<2,gsBSplineBasis<> > tbasis( new gsBSplineBasis<>(bsp) , new gsBSplineBasis<>(bsp) ) ;
    gsHBSplineBasis<2>  HB( tbasis , 3) ;
    gsVector<unsigned int> i1(2),i2(2);
    i1[0] = 4;
    i1[1] = 4;
    i2[0] = 12;
    i2[1] = 12;
//    HB.insert_box(i1,i2,1);
    i1[0] = 4;
    i1[1] = 4;
    i2[0] = 10;
    i2[1] = 10;
//    HB.insert_box(i1,i2,2);
    HB.set_xmatrix();
    gsInfo<<" Plotting "<< HB <<"\n";
    plot.basis( HB );
    

    return 0;
}
