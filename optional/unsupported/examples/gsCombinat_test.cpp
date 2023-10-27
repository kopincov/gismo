/** @file gsCombinat_test.cpp

    @brief Test of utilities related to combinatorics.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, A. Mantzaflaris
*/

#include <gismo.h>

using std::flush;
using namespace gismo;

#define TEST(a)\
    passed = (passed && (a));\
    gsInfo << ((a)? "TEST OK\n":"TEST FAIL\n");\


int main()
{
    bool passed = true;
    unsigned r;


    for (int d = 1; d<4; d++)
    {


        int cstart=0;
        int cend  =4;
        int csub_start=2;
        int csub_end  =3;

        // we test the d-dimensional cube  [start, end]x...x[start, end]
        // in the subcube we test the subcube [sub_start, sub_end]x...x[sub_start, sub_end]

        unsigned clength;
        clength     = cend-cstart+1 >0?(cend-cstart+1):0 ;
        unsigned csub_length;
        csub_length = csub_end-csub_start+1>0 ? (csub_end-csub_start+1):0 ;
        unsigned points=1;
        unsigned sub_points=1;
        unsigned cvertexes=1;
        unsigned boundaries1=1;
        unsigned boundaries2=1;
        unsigned int_length1;
        int_length1 = clength-2>0 ?clength-2: 0;
        unsigned int_length2;
        int_length2 = clength-4>0 ?clength-4: 0;

        for (int i=0;i<d;i++)
        {
            points*=clength;
            sub_points*=csub_length;
            cvertexes*=2;
            boundaries1*=int_length1;
            boundaries2*=int_length2;

        }
        boundaries1=points-boundaries1;
        boundaries2=points-boundaries2;


        gsVector<unsigned> dim(d), beg(d), end(d);
        dim.setConstant(cend-cstart+1);
        beg.setConstant(csub_start);
        end.setConstant(csub_end);


        gsVector<unsigned> cube0(d), cube1(d), curr(d);
        cube0.setConstant(cstart);
        cube1.setConstant(cend);
/*
        gsTensorLatticeIterator<unsigned> my_cube(dim);
        gsTensorLatticeIterator<unsigned> sub_cube= my_cube.getSubTensorIterator (beg,end);



        gsInfo<<"----------------------------------------\n"
              "----------------------------------------\n" <<d<<
              "D cube test : ("<< cube0.transpose() <<"), ("<< cube1.transpose() <<")\n"<<
              "----------------------------------------\n"
              "----------------------------------------\n";

        gsInfo<<"----------------------------------------\n"
              "Lattice iterator version.\n"
              "----------------------------------------\n";

        gsInfo<< "Points: \n";
        r = 0;
        do
        {
            r++;
            gsInfo<< my_cube.flatIndex()<<":("<< my_cube.multiIndex().transpose() <<") "<<flush;
        } while ( my_cube.next() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST(r==points);

        gsInfo<< "Points in subcube: \n";
        r = 0;
        do
        {
            r++;
            gsInfo<< sub_cube.flatIndex()<<":("<< sub_cube.multiIndex().transpose() <<") "<<flush;
        } while ( sub_cube.next() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST(r==sub_points);
        gsInfo<< "\nThe same backward:\n";
        r = 0;
        do
        {
            r++;
            gsInfo<< sub_cube.flatIndex()<<":("<< sub_cube.multiIndex().transpose() <<") "<<flush;
        } while ( sub_cube.previous() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST(r==sub_points);

        gsInfo<< " Boundary points with iterator: \n";
        gsVector<unsigned> thickness;
        thickness.setZero(d);
        gsTensorBoundaryIterator<unsigned> boundary(my_cube,thickness);
        r = 0;
        do
        {
            r++;
            gsInfo<< boundary.flatIndex()<<":("<< boundary.multiIndex().transpose() <<") "<<flush;
        } while ( boundary.next() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST(r==boundaries1);
        gsInfo<< "\nSame backward:\n";
        r = 0;
        do
        {
            r++;
            gsInfo<< boundary.flatIndex()<<":("<< boundary.multiIndex().transpose() <<") "<<flush;
        } while ( boundary.previous() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST( r==boundaries1);


        gsInfo<< "Points with distance <=1 from the boundary: \n";
        thickness.setOnes(d);
        gsTensorBoundaryIterator<unsigned> thick_boundary(my_cube,thickness);
        r = 0;
        do
        {
            r++;
            gsInfo<< thick_boundary.flatIndex()<<":("<< thick_boundary.multiIndex().transpose() <<") "<<flush;
        } while ( thick_boundary.next() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST(r==boundaries2);
        gsInfo<< "\nSame backward:\n";
        r = 0;
        do
        {
            r++;
            gsInfo<< thick_boundary.flatIndex()<<":("<< thick_boundary.multiIndex().transpose() <<") "<<flush;
        } while ( thick_boundary.previous() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST( r==boundaries2);

        gsInfo<< " Verteces points with iterator of the subcube: \n";
        gsTensorVertexIterator<unsigned> vertexes(sub_cube);
        r = 0;
        do
        {
            r++;
            gsInfo<< vertexes.flatIndex()<<":("<< vertexes.multiIndex().transpose() <<") "<<flush;
        } while ( vertexes.next() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST(r==cvertexes);
        gsInfo<< "\nSame backward:\n";
        r = 0;
        do
        {
            r++;
            gsInfo<< vertexes.flatIndex()<<":("<< vertexes.multiIndex().transpose() <<") "<<flush;
        } while ( vertexes.previous() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST( r==cvertexes);
*/
        gsInfo<<"----------------------------------------\n"
              "Free function version.\n"
              "----------------------------------------\n";

        gsInfo<< " Lattice points: \n";
        curr = cube0;
        r = 0;
        do
        {
            r++;
            gsInfo<< "("<< curr.transpose() <<") ";
        } while ( nextCubePoint( curr, cube0, cube1 ) );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST(r==points);


        gsInfo<< " Vertices: \n";
        curr = cube0;
        r = 0;
        do
        {
            r++;
            gsInfo<< "("<< curr.transpose() <<") ";
        } while ( nextCubeVertex( curr, cube0, cube1 ) );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST( r==cvertexes);


        gsInfo<< " Cube boundary: \n";
        curr = cube0;
        r = 0;
        do
        {
            r++;
            gsInfo<< "("<< curr.transpose() <<") ";
        } while ( nextCubeBoundary( curr, cube0, cube1 ) );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST(r==boundaries1);


        gsVector<unsigned> offset;
        offset.setOnes(d);
        gsInfo << " Cube boundary points offseted by ("<<offset.transpose()<<"): \n";
        curr = cube0;
        r = 0;
        do
        {
            r++;
            gsInfo<< "("<< curr.transpose() <<") ";
        } while ( nextCubeBoundaryOffset( curr, cube0, cube1, offset ) );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST(r==boundaries2);


        gsInfo<< "Check iteration over the "<<d<<"-cube elements (faces):\n";
        for (int i = 0; i<=d; ++i)
        {
            firstCubeElement(curr, i);
            gsInfo<<"* The "<<numCubeElements(i,d);
            switch( i )
            {
            case 0:  gsInfo<<" vertices:\n"; break;
            case 1:  gsInfo<<" edges:\n"   ; break;
            case 2:  gsInfo<<" faces:\n"   ; break;
            default: gsInfo<<" faces of dim="<<i<<":\n "; break;
            }

            index_t ne = 0;
            do
            {
                ++ne;
                GISMO_ASSERT( dimCubeElement(curr) == i, "Error in Cube element." );
                gsInfo<<"("<<curr.transpose() <<") ";
            }
            while ( nextCubeElement(curr, i) );
            TEST( ne == numCubeElements(i,d) );
        }
        

        gsInfo<<" Cube isometries:\n";
        gsVector<index_t>  perm(d);
        gsVector<bool> flip(d), upp(d);
        upp.setOnes();
        gsVector<index_t>  isometry;
        r = 0;
        flip.setZero();
        do//for all binary sequences of length d
        {
            firstPermutation(perm);
            do //for all permutations of (0,..d-1)
            {
                cubeIsometry(flip, perm, isometry);
                gsInfo<<"("<< isometry.transpose()<<") ";
                ++r;
            }
            while( nextPermutation(perm) );               
        }
        while( nextCubeVertex(flip, upp) );
        gsInfo<<"\nListed "<< r <<" isometries. ";
        TEST( r == (1<<d)*factorial(d) );
        
    }// for int d

    for (int i=1; i<8; i++)
    {
/*        gsInfo<< "\nCombinations iterator.\n"<<"\n";
        gsSimplexLatticeIterator<unsigned> combinations(i,5);
        r = 0;
        do
        {
            r++;
            gsInfo<< combinations.flatIndex()<<":("<< combinations.multiIndex().transpose() <<") "<<flush;
        } while ( combinations.next() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST( r == numCompositions(5,i) );
        gsInfo<< "\nSame backward:\n"<<"\n";
        r = 0;
        do
        {
            r++;
            gsInfo<< combinations.flatIndex()<<":("<< combinations.multiIndex().transpose() <<") "<<flush;
        } while ( combinations.previous() );
        gsInfo<< "\nThey were "<< r <<".\n";
        TEST( r == numCompositions(5,i) );



        bool tpass=true;
        gsInfo<< "\n Check index conversion ";
        combinations.reset();
        do {
            tpass = (tpass && combinations.multiIndex(combinations.flatIndex())==combinations.multiIndex());
        } while ( combinations.next());
        TEST(tpass);

*/
        gsInfo<< " Compositions of "<< 5 <<" into "<< i <<" parts: \n";
        gsVector<unsigned> curr(i);
        firstComposition(5,i,curr);
        r = 0;
        do
        {
            r++;
            //gsInfo<< "("<< curr.transpose() <<") ";
        } while ( nextComposition(curr) );
        gsInfo<< "They were "<< r <<".\n";
        TEST( r == numCompositions(5,i) );

        gsVector<unsigned> mh(4);
        mh<< 4, 2, 1, 3;
        gsInfo<< " Multi-compositions of "<< mh.transpose() <<" into "<< i <<" parts: \n";
        gsMatrix<unsigned> cmh;
        firstMultiComposition(mh, i, cmh);
        r = 0;
        do
        {
            r++;
            //gsInfo<< r <<":\n"<< cmh <<"\n";
            //gsInfo <<"fact prod:\n"<< factorial(mh.sum())/cmh.unaryExpr(std::ptr_fun(factorial)).prod() <<"\n";

        } while ( nextMultiComposition(cmh) );
        gsInfo<< "They were "<< r <<".\n";
        TEST( r == numMultiCompositions(mh,i) );
    }

    gsInfo<< " Check binomials: ";
    gsVector<unsigned>  v;
    binomial_into(8,v);
    bool c=true;

    c= c &&( v[0] == binomial<unsigned>(8,0) ) && ( v[0]== (binomial<8,0>()) );
    c= c &&( v[1] == binomial<unsigned>(8,1) ) && ( v[1]== (binomial<8,1>()) );
    c= c &&( v[2] == binomial<unsigned>(8,2) ) && ( v[2]== (binomial<8,2>()) );
    c= c &&( v[3] == binomial<unsigned>(8,3) ) && ( v[3]== (binomial<8,3>()) );
    c= c &&( v[4] == binomial<unsigned>(8,4) ) && ( v[4]== (binomial<8,4>()) );
    c= c &&( v[5] == binomial<unsigned>(8,5) ) && ( v[5]== (binomial<8,5>()) );
    c= c &&( v[6] == binomial<unsigned>(8,6) ) && ( v[6]== (binomial<8,6>()) );
    c= c &&( v[7] == binomial<unsigned>(8,7) ) && ( v[7]== (binomial<8,7>()) );
    c= c &&( v[8] == binomial<unsigned>(8,8) ) && ( v[8]== (binomial<8,8>()) );

    TEST(c);

    return ( passed ? 0 : 1 );
}



