/** @file plot_div_preserving_trans.cpp

    @brief Example programs that reconstructs a 2D vector field using differnt tarsformations

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#include <iostream>
#include <gismo.h>
#include <gsCore/gsDivConSolution.h>



;
using namespace gismo;

int main(int argc, char *argv[])
{
    // Define Geometry
    gsMultiPatch<> * patches;
    patches = new gsMultiPatch<>( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    gsMultiPatch<> * patchesSquare;
    patchesSquare = new gsMultiPatch<>( *gsNurbsCreator<>::BSplineSquare(1.0,0,0) );


    // Define and refine discretization space by refining the basis of the geometry
    std::vector< gsMultiBasis<> >  refine_bases;
    refine_bases.push_back(gsMultiBasis<>( *patches ));//Basis for x direction of the vector field
    refine_bases.push_back(gsMultiBasis<>( *patches ));//Basis for y direction of the vector field

    int numRefine = 3;
    for (int i = 0; i < numRefine; ++i)
    {
        refine_bases[0].uniformRefine();
        refine_bases[1].uniformRefine();
    }

    refine_bases[0].degreeElevate(1,0);
    refine_bases[1].degreeElevate(1,0);



    gsInfo << "Degrees of basis for x direction: "<< refine_bases[0].degree(0,0)<< " and "<< refine_bases[0].degree(0,1)<<"\n";
    gsInfo << "Degrees of basis for y direction: "<< refine_bases[1].degree(0,0)<< " and "<< refine_bases[1].degree(0,1)<<"\n";

    int szX = refine_bases[0][0].size();
    int szY = refine_bases[1][0].size();

    //Coefficients for vector field in the parametric domain
    //In the parametric domain the field is (0,1) everywhere
    gsMatrix<> coeffs(szX+szY,1);
    for (index_t k=0; k< szX + szY; ++k)
    {
        if (k < szX)
            coeffs(k,0) = 0;
        else
            coeffs(k,0) = 1;
    }
    //Resizing the matrix since the different methods have differnt format
    gsMatrix<> coeffsIC = coeffs;
    coeffsIC.resize(szX,2);
    gsMatrix<> coeffsCP2 = coeffsIC;

    std::vector<gsBasis<> *> bases_vec;
    bases_vec.push_back(&refine_bases[0][0]);
    bases_vec.push_back(&refine_bases[1][0]);


    //Creating a vector field in the physical domain with
    //divergence preserving transformation (Piola).
    gsMultiPatch<> patchesTest = *patches;
    gsDivConSolution<real_t> funcPiola(coeffs,patchesTest[0],bases_vec);
    gsField<> fieldPiola( *patches, funcPiola);

    //Creating a vector field in the physical domain with
    //inverse composition (standard IGA)
    gsFunction<>::uPtr funcIGA= refine_bases[0][0].makeGeometry( give(coeffsIC) );
    gsField<> fieldIGA( *patches , *funcIGA) ;

    //Creating a vector field in the physical domain with
    //no transform (parametric space)
    gsFunction<>::uPtr funcNo= refine_bases[0][0].makeGeometry( give(coeffsCP2) );
    gsField<> fieldNoTrans( *patchesSquare , *funcNo) ;

    // Write solution to paraview files
    gsWriteParaview<>( fieldPiola, "transPiola", 500);
    gsWriteParaview<>( fieldIGA, "transIGA", 500);
    gsWriteParaview<>( fieldNoTrans, "transNo", 500);

    delete patches;
    delete patchesSquare;

    gsInfo << "Test is done: Exiting" << "\n";
    return 0;
}
