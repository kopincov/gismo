/** @file tutorial_make_triangular_topology.cpp

    @brief Constructs a multi patch, which represents a triangle.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

#include <gismo.h>

using namespace gismo;


int main(int argc, char* argv[])
{
    gsMultiPatch<> mp;


    gsVector<> bottom_left(2);   bottom_left  << -1,   -1;
    gsVector<> left(2);          left         << -0.5, 0.5;
    gsVector<> bottom(2);        bottom       <<  0,   -1;	
    gsVector<> mid(2);           mid          <<  0,    0;	
    gsVector<> right(2);         right        <<  0.5,  0.5;
    gsVector<> bottom_right(2);  bottom_right <<  1,   -1;
    gsVector<> top(2);           top          <<  0,    1;

    gsKnotVector<> kv1(0.0, 1.0, 0, 2);
    gsKnotVector<> kv2(0.0, 1.0, 0, 2);


    // ----------------------------------------------------------------------
    // patch 1
    gsMatrix<> coefs1(4, 2);
    coefs1 << 
	bottom_left(0), bottom_left(1),
	bottom(0), bottom(1),
	left(0), left(1),
	mid(0), mid(1);
    mp.addPatch(gsTensorBSpline<2>(kv1, kv2, coefs1));

    
    // ----------------------------------------------------------------------
    // patch 2
    gsMatrix<> coefs2(4, 2);
    coefs2 << 
	bottom(0), bottom(1),
	bottom_right(0), bottom_right(1),
	mid(0), mid(1),
	right(0), right(1);
    mp.addPatch(gsTensorBSpline<2>(kv1, kv2, coefs2));


    // ----------------------------------------------------------------------
    // patch 3    
    gsMatrix<> coefs3(4, 2);
    coefs3 << 
	left(0), left(1),
	mid(0), mid(1),
	top(0), top(1),
	right(0), right(1);
    mp.addPatch(gsTensorBSpline<2>(kv1, kv2, coefs3));
    
    
    // toopology
    mp.computeTopology();
    
    // saving 
    gsInfo << "Saving triangular topology.\n";
    gsFileData<> fd;
    fd << mp;
    fd.dump("triangular_topology");

    return 0;
}


