/** @brief Reads, refines and writes THB-splines to and from gismo xml files. Can be used to update files from the
    times, when the assumption that $Omega^k$ is a union of cells of level $k-1$ was not yet imposed. If you intend
    to do this, first ask Angelos and then, eventually, comment the line with adaptiveAlignedSplit(iBox, m_index_level) to
    adaptiveSplit(iBox) in gsHDomain.hpp, currently on lines 161 and 162.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Dominik Mokris, dominik.mokris@jku.at
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;

gsTHBSpline<2> * readTHBSplineFromFile( std::string filename )
{
  gsTHBSpline<2>::uPtr result;
  gsFileData< real_t > data( filename );
  result = data.getFirst< gsTHBSpline<2> >();
  return result.release();
}

void writeTHBSplineIntoFile( gsTHBSpline<2>& THBSpline, std::string filename )
{
  gsFileData<> newdata;
  newdata << THBSpline;
  newdata.dump( filename );
}

unsigned powFloor(unsigned num, unsigned pow)
{
    return num - (num % pow);
}

unsigned powCeil(unsigned num, unsigned pow )
{
    unsigned mod = num % pow;
    if( mod > 0 )
        return num - mod + pow;
    else
        return num;
    //return ( mod ? num - (num % pow) + 1 : num );
}

std::vector<index_t> getBoxesAndModifyThemToAvoidLShapes (const gsTHBSpline<2> & THBSpline )
{
    gsMatrix<index_t> b1;
    gsMatrix<index_t> b2;
    gsVector<index_t> level;
    gsHDomain<2> tree = THBSpline.basis().tree();
    tree.getBoxes( b1, b2, level );
    std::vector<index_t> boxes( 5 * level.size(), 0 );
    { /* the scope of i and j;
       I cannot declare two different types of variables inside for, so I have to declare them before.
       See http://stackoverflow.com/questions/2687392/is-it-possible-to-declare-two-variables-of-different-types-in-a-for-loop
       */
        index_t i;
        unsigned j;
        for( i=0, j=0 ; i < level.size(); i++, j+=5 )
        {
            //gsInfo << "before:" << level(i) <<","<< b1(i,0) <<","<< b1(i,1) <<","<< b2(i,0) <<","<< b2(i,1) << "\n";
            unsigned pow = math::ipow( 2, level(i) );
            unsigned lev = level(i);
            boxes[j] = lev;

            gsVector<index_t,2> glIndex;
            gsVector<index_t,2> locIndex;

            //gsInfo << "middle:"<<lev<<",";
            glIndex(0) = b1(i,0);
            glIndex(1) = b1(i,1);
            tree.global2localIndex(glIndex,lev,locIndex);
            //gsInfo << locIndex(0) <<","<< locIndex(1) <<",";
            boxes[j+1] = powFloor(locIndex(0),pow);
            boxes[j+2] = powFloor(locIndex(1),pow);

            glIndex(0) = b2(i,0);
            glIndex(1) = b2(i,1);
            tree.global2localIndex(glIndex,lev,locIndex);
            //gsInfo << locIndex(0) <<","<< locIndex(1) << "\n";
            boxes[j+3] = powCeil(locIndex(0),pow);
            boxes[j+4] = powCeil(locIndex(1),pow);

            //gsInfo << "after: " << boxes[j] <<","<< boxes[j+1] <<","<< boxes[j+2] <<","<< boxes[j+3] <<","<< boxes[j+4] << "\n";
        }
    } // end of the scope of i and j

    return boxes;
}

bool succeededInReadRefineAndWrite( std::string filenameIn, std::string filenameOut )
{
  gsTHBSpline<2> *pointer = readTHBSplineFromFile( filenameIn );
  if( pointer == 0 )
      return false;
  else
  {
      gsTHBSpline<2> THBSpline = *pointer;
      std::vector<index_t> boxes = getBoxesAndModifyThemToAvoidLShapes( THBSpline );
      THBSpline.refineElements( boxes );
      writeTHBSplineIntoFile( THBSpline, filenameOut );
      gsInfo << "Success, " << filenameOut << " done." << "\n";
      delete pointer;
      return true;
  }
}

bool MTUconversionsHaveSucceeded()
{
    return (
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/a02-9_e2_d3_t-6_e2-5.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/a02-9_e2_d3_t-6_e2-5_converted.xml" )) &&
                // Yields problem during reading:
                //( succeededInReadRefineAndWrite( "face.xml", "face_converted.xml" )) &&
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/fillet-irreg-9_e2_d3_t-6_e2-glob2.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/fillet-irreg-9_e2_d3_t-6_e2-glob2_converted.xml" )) &&
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/fillet-irreg-9_e2_d3_t-6_e2-glob2a.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/fillet-irreg-9_e2_d3_t-6_e2-glob2a_converted.xml" )) &&
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/fillet-irreg-9_e2_d3_t-6_e22.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/fillet-irreg-9_e2_d3_t-6_e22_converted.xml" )) &&
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/fillet-irreg-9_e2_d3_t-6_e22a.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/fillet-irreg-9_e2_d3_t-6_e22a_converted.xml" )) &&
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/gsThbs_MTU_00.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/gsThbs_MTU_00_converted.xml" )) &&
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/gsThbs_MTU_00_bump.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/gsThbs_MTU_00_bump_converted.xml" )) &&
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/gsThbs_MTU_blade_029.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/gsThbs_MTU_blade_029_converted.xml" )) &&
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/gsThbs_MTU_fillet_bump_01.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/gsThbs_MTU_fillet_bump_01_converted.xml" )) &&
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/gsThbs_MTU_fillet_bump_02.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/gsThbs_MTU_fillet_bump_02_converted.xml" )) &&
                // Yields problem during reading:
                //( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/gsThbs_adapt_face_10.xml",
                //                                 "../../Gabor_data_remade/conv_avoidL/gsThbs_adapt_face_10_converted.xml" )) &&
                ( succeededInReadRefineAndWrite( "../../Gabor_data_remade/gismo_files/thbs_face_3levels.xml",
                                                 "../../Gabor_data_remade/conv_avoidL/thbs_face_3levels_converted.xml" )));
}

int main()
{
    bool MTUConversion = false;
    if( MTUConversion ) // Uses files that are not part of the repository.
    {
        if(  MTUconversionsHaveSucceeded() )
            return 0;
        else
            return 1;
    }
    else // Default test.
    {
        std::string filenameIn = "thbs_00.xml";
        std::string filenameOut = "thbs_00_refined.xml";
        if( succeededInReadRefineAndWrite( filenameIn, filenameOut ) )
            return 0;
        else
            return -1;
    }
}
