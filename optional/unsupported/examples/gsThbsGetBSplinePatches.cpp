/** @brief Reads a THB surface from an xml file, converts it to a set of patches and compares values of several variables
 * to have a reasonable probability that everything went as it should. Should prevent the off-by one bug in level numbering
 * from reappearing and living in our source code again.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Dominik Mokris, dominik.mokris@jku.at
*/

#include <iostream>
#include <fstream>
#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;


bool answerWasCorrect( const gsMatrix<index_t>& yourAnswer, const gsMatrix<index_t>& correctAnswer )
{
    if( yourAnswer == correctAnswer )
        return true;
    else
    {
        gsInfo << "Wrong answer!" << "\n";
        gsInfo << "Your answer:" << "\n" << yourAnswer << "." << "\n" << "\n";
        gsInfo << "Correct answer:" << "\n" << correctAnswer << "." << "\n" << "\n";
        return false;
    }
}

bool facePatches_test( std::string filename )
{
    gsFileData<real_t>  data( filename );
    gsInfo<<"Running file: "<< filename  <<"\n"<<"\n";
    gsTHBSpline<2>::uPtr spline;
    spline = data.getFirst< gsTHBSpline<2> >();
    //gsInfo << *spline << "\n";
    gsTHBSplineBasis<2>  basis = spline->basis();

    gsMatrix<real_t> cp;
    gsVector<index_t> level;
    gsMatrix<index_t> lowerLefts, upperRights, nvertices;

    std::vector<std::vector<std::vector< std::vector<real_t> > > >* trim_curves =
      new std::vector<std::vector<std::vector< std::vector<real_t> > > >();

    basis.getBsplinePatches_trimming(spline->coefs(), cp, lowerLefts, upperRights, level, nvertices, *trim_curves); //automatic splitting

    // It would be good to somehow test also the trim curves and coefs and also cp to see that they haven't changed.

    delete trim_curves;

    gsMatrix<index_t> correct_lowerLefts(3,2);
    correct_lowerLefts <<0,8,0,0,6,24;
    gsMatrix<index_t> correct_upperRights(3,2);
    correct_upperRights << 16,12,24,48,18,42;
    gsVector<index_t> correct_level(3);
    correct_level << 1,2,3;
    gsMatrix<index_t> correct_nvertices(3,2);
    correct_nvertices << 7,4,15,27,15,21;

    return (    answerWasCorrect( lowerLefts, correct_lowerLefts ) &&
                answerWasCorrect( upperRights, correct_upperRights ) &&
                answerWasCorrect( level, correct_level ) &&
                answerWasCorrect( nvertices, correct_nvertices ));
}

/*void seeTheStuffInParaview( std::string filename )
{
    gsFileData<real_t> data(filename);
    gsTHBSpline<2>::uPtr THB;
    THB = data.getFirst< gsTHBSpline<2> >();
    gsWriteParaview( *THB, "thb_out",100,false, true);
    }*/

int main()
{
  std::string filename = "thbs_face_3levels.xml";
  //  seeTheStuffInParaview(filename);
  if( facePatches_test( filename ) )
  {
      gsInfo << "Test successful.\n"<<"Congratulations!\n";
      return 0;
  }
  else
  {
      gsInfo << "Test failed.\n"<<"Please, take the time and care to find and correct the mistakes.\n";
      return 1;
  }
}
