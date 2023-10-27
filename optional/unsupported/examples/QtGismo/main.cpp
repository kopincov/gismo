
#include <QtGui/QApplication>

#include "gsMainWindow.h"

#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsBSpline.h>
#include "gsThbs/gsQuadTree.h"


using namespace gismo;


int main(int argc, char *argv[])
{


  // some evaluation points
  gsMatrix<real_t,1> u(1,3) ;

  u << 0.05, 0.12 , 0.75  ;
  
  // start with a knot vector
  gsKnotVector<> KV (0, 1, 9 ,3 ) ;

  // make B-Spline basis by knotvector and degree
  gsBSplineBasis<> * bsp = new gsBSplineBasis<>(KV, 2) ;
  std::cout<< "The "<< *bsp << std::endl;



    QApplication a(argc, argv);
    gsMainWindow w;
    w.show();



    return a.exec();
}


