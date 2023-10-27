/**  gsInterfaceMapTest.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):   A. Bressan
    Created on:  2014-11-26

    Test the construction of the affine map between two patches.
    We enumerate all possible maps between two faces of 1,2 and
    3D boxes.
**/


#include  <gismo.h>
#include  <gsUtils/gsCombinatorics.h>

using namespace gismo;

void printMapping (const gsMatrix<> &orig, const gsMatrix<> &mapped)
{
    gsInfo<<"F:";
    for(int i=0;i< orig.rows()-1;++i)
    {
        if (orig(i,0)==orig(i,1))
            gsInfo<<"{"<<orig(i,0)<<"}x";
        else
            gsInfo<<"["<<orig(i,0)<<","<<orig(i,1)<<"]x";
    }
    if (orig(orig.rows()-1,0)==orig(orig.rows()-1,1))
        gsInfo<<"{"<<orig(orig.rows()-1,0)<<"}->";
    else
        gsInfo<<"["<<orig(orig.rows()-1,0)<<","<<orig(orig.rows()-1,1)<<"]->";
    for(int i=0;i< mapped.rows()-1;++i)
    {
        if (mapped(i,0)==mapped(i,1))
            gsInfo<<"{"<<mapped(i,0)<<"}x";
        else
            gsInfo<<"["<<mapped(i,0)<<","<<mapped(i,1)<<"]x";
    }
    if (mapped(mapped.rows()-1,0)==mapped(mapped.rows()-1,1))
        gsInfo<<"{"<<mapped(mapped.rows()-1,0)<<"}  ";
    else
        gsInfo<<"["<<mapped(mapped.rows()-1,0)<<","<<mapped(mapped.rows()-1,1)<<"]  ";
}

void setOrientation (index_t s, gsVector<bool> &o)
{
    for (index_t i=0; i<o.rows();++i)
    {
        o(i)= !(s& 1<<i);
    }
}

void test (int dim)
{
    gsMatrix<> box1(dim,2);
    gsMatrix<> box2(dim,2);
    gsMatrix<> origin(dim,2);
    gsMatrix<> target(dim,2);

    for (index_t i=0; i<2*dim; ++i)
    {
        box1(i/2,i%2)=i;
    }
    for (index_t i=0; i<2*dim; ++i)
    {
        box2(i/2,i%2)=2*dim+i;
    }

    gsVector<index_t> per(dim);
    firstPermutation(per);
    gsVector<bool>    o(dim);
    do
    {
        for ( index_t os=0; os< 1<<dim; os++)
        {
            setOrientation(os,o);
            gsInfo<<"P:["<<per.transpose()<<"], O:["<< o.transpose()<<"], ";
            gsAffineFunction<real_t> map(per,o, box1, box2 );
            index_t subdim=0;
            origin=box1;
            do {
                map.eval_into(origin,target);
                printMapping(origin,target);
                origin(subdim,1)= origin(subdim,0);
                gsInfo<<" ";
                ++subdim;
            } while (subdim<dim);
            gsInfo<<"\n";
        }
        gsInfo<<"\n";
    } while (  nextPermutation(per) );
}

int main()
{
    gsInfo<<
        "P: is the permutation of the axis\n"
        "O: is the 1 if the axis and detination have the same orientation\n"
        "F: is the map on the box, the map on lower dimensional faces is printed too\n\n"
        "THE RESULT MUST BE CHECKED MANUALLY\n\n";

    gsInfo<<"Test affine maps between segments\n\n";
    test(1);
    gsInfo<<"Test affine maps between rectangles\n\n";
    test(2);
    gsInfo<<"Test affine maps between cuboids\n\n";
    test(3);
//    gsInfo<<"Test affine maps between hyper-cuboids\n\n";
//    test(4);

    return 0;
}

