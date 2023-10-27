/** @file gsMatrixViewTest.cpp

    @brief ...

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gismo.h>
#include <gsSolver/gsLowRankCorrectedOp.h>
#include <gsSolver/gsKronecker.h>
#include <gsSolver/gsGradientMethod.h>

#include <iostream>
#include <gsUtils/gsStopwatch.h>


using namespace gismo;


struct MatrixView : public Eigen::Map< gsMatrix<real_t>::Base > {
    
    //MatrixView( real_t * data, int sz ) : Eigen::Map< gsVector<real_t>::Base >(data,sz) {}
    MatrixView( gsVector<real_t>& other ) : Eigen::Map< gsMatrix<real_t>::Base >(other.data(),other.rows(),1) {}
    MatrixView( gsVector<real_t>& other, int begin, int length ) : Eigen::Map< gsMatrix<real_t>::Base >(other.data()+begin,length,1) {}
    
    MatrixView( gsMatrix<real_t>& other ) : Eigen::Map< gsMatrix<real_t>::Base >(other.data(),other.rows(),other.cols()) {}
    MatrixView( gsMatrix<real_t>& other, int begin, int length ) : Eigen::Map< gsMatrix<real_t>::Base >(other.data()+begin,length,1) { GISMO_ENSURE (other.cols()==1,"This only works for one column. Think why!"); }
    
};


int main(int argc, char** argv)
{
    gsCmdLine cmd("This is just a test which is testing something.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    index_t jj=1;
    
    gsMatrix<real_t> v(10,jj);
    for (int i=0; i<10; ++i)
        for (int j=0; j<jj; ++j)
            v(i,j)=1+2*i+100*j;
    
    
    gsInfo << "v =" << v.transpose() << "\n";
    
    MatrixView w(v,2,3);
    
    gsInfo << "w =" << w.transpose() << "\n";
    
    w(0,0) = 99;
    
    gsInfo << "v =" << v.transpose() << "\n";
    gsInfo << "w =" << w.transpose() << "\n";
    
    MatrixView u(v,6,3);

    u = w;

    gsInfo << "u =" << u.transpose() << "\n";
    gsInfo << "v =" << v.transpose() << "\n";
    gsInfo << "w =" << w.transpose() << "\n";
    

    return 0;
}


