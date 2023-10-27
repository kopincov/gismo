/** @file gsOptParam

    @brief Parametrization by optimization

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

// #include <gsModeling/gsSpringPatch.h>
// #include <gsModeling/gsCoonsPatch.h>
#include <gsOptimizer/gsOptParamJacobian.h>
#include <gsOptimizer/gsOptParamWinslow.h>
#include <gsOptimizer/gsOptParamLiao.h>
//#include <gsOptimizer/gsOptParamAreaOrth.h>


using namespace gismo;

int main(int argc, char* argv[])
{
    // Reading the filename from command line argument
    std::string fn;
    index_t numKnots   = 0;
    index_t func = 0;

    gsCmdLine cmd("Constructs a spring patch given a domain boundary.");
    cmd.addPlainString("filename", "File containing boundary data", fn);
    cmd.addInt("r", "urefine", "initial uniform refinement steps" 
               "(num. of knots in each knot-span)", numKnots);
    
    cmd.addInt("f","func" ,"Functional to use for optimization\n" 
               " 0 MaxJac\n 1 Winslow\n 2 Liao \n 3 Area orthogonality", func);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    // Read in the input file
    gsGeometry<>::uPtr geo = gsReadFile<>(fn);

    // Uniform refine twice (ie. add numKnofts knots uniformly inside
    // every knot-span
    if (numKnots>0)
        geo->uniformRefine(numKnots);


    gsMatrix<> pts;
    gsVector<index_t> corner;
    int nPositive = 0;
    real_t val;
    corner.setZero(geo->parDim());
    do
    {
        pts = corner.cast<real_t>();
	gsDebugVar(pts);
	gsDebugVar(geo->jacobian(pts));
        val = geo->jacobian(pts).determinant();
        gsInfo<< "Corner ("<<corner.transpose() <<") :" << val <<"\n";

        if ( val >= -1e-12 )
            nPositive++;
    }
    while( nextCubeVertex(corner) );
    
    if ( nPositive == 0 )
    {
        gsWarn<< "Switching the orientation.\n"; 
        if(gsTensorBSpline<3> * r = dynamic_cast<gsTensorBSpline<3>*>(geo.get()))
	  r->swapDirections(0,1);
        if(gsTensorBSpline<2> * r = dynamic_cast<gsTensorBSpline<2>*>(geo.get()))
	  r->swapDirections(0,1);
        corner.setZero(geo->parDim());
        do
        {
            pts = corner.cast<real_t>();
            val = geo->jacobian(pts).determinant();
            gsInfo<< "Corner ("<<corner.transpose() <<") :" << val <<"\n";
        }
        while( nextCubeVertex(corner) );

    }
    else if ( nPositive < (1<< geo->parDim()) )
    {
        gsWarn<< "Value of the Jacobian at the corners have different signs, quitting.\n"; 
        return 0;
    }

    // Setup optimizer
    gsOptParam<real_t> * optimizer = NULL;
    std::string choice;
    switch (func)
    {
    case 0:
        optimizer = new gsOptParamJacobian<real_t>(*geo);
        choice = "jac";
        break;
    case 1:
        optimizer = new gsOptParamWinslow<real_t>(*geo);
        choice = "winslow";
        break;
    case 2:
        optimizer = new gsOptParamLiao<real_t>(*geo);
        choice = "liao";
        break;
    default:
        gsInfo<<"Bad option.\n";
        return 0;
        break;
    }

    // Get initial Jacobian
    gsVector<> tmp;
    optimizer->currentJacobianCoefs(tmp);

    // Keep initial design info
    real_t minC = tmp.minCoeff();
    real_t maxC = tmp.maxCoeff();

    // Print initial design info
    gsInfo << "\nInitial min/max Jac-coeff: "<< minC <<" / " <<maxC <<"\n\n";

    // Run optimizer
    optimizer->solve();

    // Print some details in the output
    optimizer->print(gsInfo);

    // Print initial design info again
    gsInfo << "\n\nInitial min/max Jac-coeff: "<< minC <<" / " <<maxC;

    // Print final design info
    optimizer->currentJacobianCoefs(tmp);
    gsInfo << "\nFinal min/max Jac-coeff: "
           << tmp.minCoeff() <<" / "
           << tmp.maxCoeff() <<"\n";
    gsInfo << "Number of iterations : " << optimizer->iterations() <<"\n";
    gsInfo << "Final objective value: " << optimizer->objective() <<"\n";

    // The final coefficients
    const gsMatrix<> & finalCoefs = optimizer->currentCoefs();

    if ( argc < 2 ) // Are we running the default example ?
    {
        const real_t dist = (finalCoefs - geo->basis().anchors().transpose() ).norm();
        gsInfo<<"Distance from the Greville points: "<< dist <<"\n";
        // Delete the memory 
        return ( dist < 1e-8 ? 0 : 1); // Check result
    }

    const gsBasis<> & Jbasis = optimizer->jacBasis();
    gsMatrix<> Jgraph = Jbasis.anchors().transpose();
    Jgraph.conservativeResize(Eigen::NoChange, 3);
    Jgraph.col(2) = tmp.array() / tmp.maxCoeff(); //scale down
    gsGeometry<>::uPtr finalJacobian = Jbasis.makeGeometry(Jgraph);
    std::string fname = std::string("optimized_jacobian_graph_") + choice;
    gsFileData<> outJacFile;
    outJacFile << *finalJacobian;
    outJacFile.dump(fname);
    gsInfo<<"Wrote final jacobian to "<<fname<<".xml .\n";

    // Write result to an xml file
    gsGeometry<>::uPtr geo2 = geo->basis().makeGeometry( finalCoefs );
    gsFileData<> outFile;
    outFile << *geo2;
    fname = std::string("optimized_geometry_") + choice;
    outFile.dump(fname);
    gsInfo<<"Wrote result to "<<fname<<".xml .\n";

 /*
    // Plot final value of jacobian
    gsMultiPatch<> geo2mp(geo2.clone());
    gsField<>::uPtr jacD = gsFieldCreator<>::jacDet(geo2mp);
    gsWriteParaview(*jacD, "optGeom_jacDet", 50000);
//*/

        // Free memory 
        delete optimizer;
        return 0;
}
