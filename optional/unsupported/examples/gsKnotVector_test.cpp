/** @file gsKnotVector_test.cpp

    @brief File testing the knot vector classes.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

#include <iostream>

using namespace gismo;


int main()
{


  ///////////////////////////////////////////////////
  // Knot vector
  ///////////////////////////////////////////////////
  gsInfo<< "------------ Knot vectors tests ------------" << "\n";

  gsKnotVector<> KV (-1, 0, 3,3, 1 ) ;
  gsInfo<< "A knot vector:"<< KV << "\n";
  gsInfo<< "Number of spans:"<< KV.uSize() - 1 << "\n";
  std::vector<index_t> knmult = KV.multiplicities();
  gsInfo<< "Multiplicities:"<< gsAsVector<index_t>(knmult).transpose() << "\n";

  gsInfo<< "Full view:"<< KV.detail() << "\n";
  
  gsInfo<< "Matrix view:"<< KV.asMatrix() << "\n";

  gsInfo<< "Number of spans:"<< KV.uSize() - 1 << "\n";

  gsInfo<< "Copied knot vector:"<<  gsKnotVector<>(KV) << "\n";

  KV.addConstant(1);
  gsInfo<< "Shifted by 1 KV :"<< KV << "\n";

  gsInfo << "Breaks iterator (knot-vector):";
  for (gsKnotVector<>::const_uiterator it=KV.ubegin(); it!=KV.uend(); ++it)
    gsInfo << " " << *it;
  gsInfo << "\n";

  std::vector<real_t> breaks = KV.unique() ;
  gsKnotVector<> KV2( breaks, 2, 0 ) ;
  gsInfo<< "Knot vector with double knots: "<< KV2 << "\n";
  gsInfo<< "Number of spans:"<< KV.uSize() - 1 << "\n";

  gsInfo << "Breaks reverse iterator (knot-vector):";
  for (gsKnotVector<>::reverse_uiterator it=KV.urbegin(); it!=KV.urend(); ++it)
    gsInfo << " " << *it;
  gsInfo << "\n";

  gsInfo <<"\n";
  gsKnotVector<> uKV, gKV;
  gsInfo<< "Uniform knot-vector:\n";
  uKV.initUniform(11, 5, 1 );
  gsInfo<< uKV.detail() ;
  gsInfo<< "Graded knot-vector (factor 0.5):\n";
  gKV.initGraded(11, 4, 0.5, 1 );
  gsInfo<< gKV.detail() ;
  gsInfo <<"\n";

  // some evaluation points
  gsMatrix<real_t,1> u(1,5) ;
  u << 0, 0.4, 0.5 , 0.98, 1  ;

  gsBSplineBasis<> bsp(KV);
  gsInfo<< "The "<< bsp << "\n";
  gsInfo<< "Evaluate non-zero basis functions on "<< u <<": \n";
  gsInfo<<  bsp.eval( u )  <<"\n";
  gsMatrix<> grev = bsp.anchors() ;
  gsInfo<< "The Greville points: "<< "\n";
  gsInfo<<  grev << "\n";
  for ( index_t i = 0; i < bsp.size(); ++i)
      gsInfo<<  KV.greville(i) << " ";
  gsInfo<<  "\n";

  grev.transposeInPlace();

  gsBSpline<real_t> CB ( bsp , grev ) ;
  gsInfo<< "The "<< CB << "\n";

  gsTensorBSplineBasis<3, real_t>
      CTB(bsp.clone().release(), bsp.clone().release(), bsp.clone().release());
  gsInfo<< "A "<< CTB << "\n";

  gsBSpline<> B ( bsp , grev ) ;
  gsInfo<< "The "<< B << "\n";


  gsKnotVector<>::iterator knt, a;
  a= KV2.begin() + 4;
  for ( knt = KV2.begin(); knt != a ; ++knt)
      gsInfo<< *knt << ", ";  
  gsInfo<< "\n";  

  gsInfo<<" ------- "<< "\n";  
  index_t p = 2;

  gsMatrix<> cc2 = memory::make_unique( KV2.greville() )->transpose();


  gsMatrix<> cc = cc2.middleRows(2,p+1);

  gsInfo<< "coefs  "<< cc.transpose() << "\n";
  gsBoehmSingle ( KV2.begin()+2 , cc, p, (real_t)(0.3) );
  gsInfo<< "coefs  "<< cc.transpose() << "\n";

  gsInfo<<" ------- "<< "\n";  

  gsInfo<< "Gcoefs "<< cc2.transpose()  << "\n";
  gsInfo<< "KV2 "<< KV2 << "\n";
  gsBoehmSingle ( KV2, cc2 , (real_t)(0.8) );
  gsInfo<< "Gcoefs "<<  cc2.transpose()  << "\n";
  gsInfo<< "KV2 "<< KV2 << "\n";

  gsKnotVector<> KV3( breaks, 2, 0 ) ;
  gsMatrix<> cc3 = memory::make_unique( KV3.greville() )->transpose();

  gsInfo<<"\n"<< "KV3 "<< KV3 << "\n";
  gsInfo<< "coeff 3: "<< cc3.transpose()  << "\n";
  std::vector<real_t> knots;
  knots.push_back(0.8);  gsBoehmRefine(KV3,cc3,p, knots.begin(), knots.end());
  gsInfo<< "KV2 "<< KV2 << "\n";
  gsInfo<< "Gcoefs "<<  cc2.transpose()  << "\n";
  gsInfo<<"\n"<< "KV3 "<< KV3 << "\n";
  gsInfo<< "refined  "<< cc3.transpose()  << "\n";

  return 0;

}
