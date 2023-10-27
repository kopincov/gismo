/** gsHFittingLvlConstrained.cpp
 *
 * For testing of the class gsHFittingLvlConstrained.
 *
 * This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Author(s): F. Buchegger.
 */

#include <iostream>

#include <gismo.h>
#include <gsHSplines/gsHFittingLvlConstrained.h>


using namespace gismo;

bool checkResults( gsMatrix<> results )
{
    bool result = true;
    /*if( ( results(0,0) != 81 ) || (results(1,0) != 225 ) || (results(2,0) != 648) || (results(3,0) < 1710) )
    {
        // Suprisingly (for me), the number of the boxes differs based on the compiler.
        gsWarn << "Sizes of the bases seem to be incorrect.\n";
        result = false;
    }*/
    if( ( results(0,1) != 1 ) || (results(1,1) != 2 ) || (results(2,1) != 3) || (results(3,1) != 4) )
    {
        gsWarn << "Max inserted levels seem to be incorrect.\n";
        result = false;
    }
    if( ( results(0,2) > 6.0e-5 ) || (results(1,2) > 2.6e-5 ) || (results(2,2) > 2.6e-5 ) || (results(3,2) > 2.6e-5) )
    {
        gsWarn << "Minimum errors seem to be incorrect.\n";
        result = false;
    }
    if( ( results(0,3) > 3.2e-1 ) || (results(1,3) > 3.2e-1 ) || (results(2,3) > 1.3e-2) || (results(3,3) > 6.2e-3) )
    {
        gsWarn << "Maximum errors seem to be incorrect.\n";
        result = false;
    }
    if( ( results(0,4) > 1.0e-5 ) || (results(1,4) > 1.0e-5 ) || (results(2,4) > 1.0e-5) || (results(3,4) > 1.0e-5) )
    {
        gsWarn << "Values of the magical fourth item in results seem to be incorrect.\n";
        result = false;
    }
    if( ( results(0,5) > 1.3 ) || (results(1,5) > 3e-1 ) || (results(2,5) > 7e-4 ) || (results(3,5) > 2e-4) )
    {
        gsWarn << "Errors seem to be incorrect.\n";
        result = false;
    }
    return result;
}

bool checkBoxes( gsMatrix<unsigned> lowerLefts, gsMatrix<unsigned> upperRights, gsVector<unsigned> levels )
{
    gsMatrix<unsigned> correctLowerLefts(31,2);
    gsMatrix<unsigned> correctUpperRights(31,2);
    gsVector<unsigned> correctLevels(31);

    bool result = true;
    gsVector<unsigned> correctValues(31);
    correctValues << 42,42,38,38,38,36,36,36,32,32,32,32,20,20,26,26,20,20,16,16,16,16,10,10,10,4,4,4,8,4,0;
    correctLowerLefts.col(0) = correctValues;
    correctValues << 4,0,42,38,16,10,16,0,36,26,16,0,32,26,20,16,4,0,16,36,26,0,42,38,0,26,22,4,0,0,0;
    correctLowerLefts.col(1) = correctValues;
    correctValues << 48,48,42,42,42,42,38,42,36,36,36,36,32,32,32,32,32,32,26,20,20,20,16,16,16,10,10,10,10,8,4;
    correctUpperRights.col(0) = correctValues;
    correctValues << 48,4,48,42,38,16,48,10,48,36,26,16,48,32,26,20,16,4,26,48,36,16,48,42,38,48,26,22,4,4,48;
    correctUpperRights.col(1) = correctValues;
    correctLevels << 4,3,4,3,4,3,4,4,4,3,4,3,4,3,4,3,4,3,4,4,3,3,4,3,4,4,3,4,3,2,2;

    if( correctLowerLefts != lowerLefts )
    {
        gsWarn << "Lower left box corners seem to be incorrect.\n";
        gsInfo<<"Your answer:"<<"\n"<<lowerLefts<<"\n"<<"Correct answer:"<<"\n"<<correctLowerLefts<<"\n";
        result = false;
    }
    if( correctUpperRights != upperRights )
    {
        gsWarn << "Upper right box corners seem to be incorrect.\n";
        gsInfo<<"Your answer:"<<"\n"<<upperRights<<"\n"<<"Correct answer:"<<"\n"<<correctUpperRights<<"\n";
        result = false;
    }
    if( levels != correctLevels )
    {
        gsWarn << "Levels seem to be incorrect.\n";
        gsInfo<<"Your answer:"<<"\n"<<levels<<"\n"<<"Correct answer:"<<"\n"<<correctLevels<<"\n";
        result = false;
    }
    return result;
}


int main(int argc, char *argv[])
{
    bool plot = false;
    bool defaultExample = false; // Reminds whether we've used the default example.
    bool result = true; // Overall result of the test.
    std::string filename = "face.xml";
    gsGeometry<>::uPtr inputSurf;
    int np = 100;
    int nd = 100;
    int err_type = 1;
    int numIter = 4;
    real_t lambda = 0.0000001;
    real_t err_threshold = 0.0001;

    gsCmdLine cmd("Hi, give me a file (.xml) with some hierarchical splines.");
    cmd.addPlainString("filename", "File containing hierarchical splines (.xml).", filename);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // The file data
    gsFileData<>  data( filename );
    /*if ( data.has< gsHBSpline<2> >() )
    {
        gsInfo<<"The HB spline functions are not fully functional"<<"\n";
        return 0;
        //hbs = data.getFirst< gsHBSpline<2> >();
    }*/
    if ( data.has< gsTHBSpline<2,real_t> >() )
    {
        inputSurf = data.getFirst< gsTHBSpline<2,real_t> >();
    }
    else
    {
        if ( data.has< gsTensorBSpline<2,real_t>  >() )
        {
            inputSurf = data.getFirst< gsTensorBSpline<2,real_t>  >();
        }
        else
        {
            gsInfo<<"Wrong input file"<<"\n";
            return 0;
        }
    }

    gsInfo<< "  Got "<< *inputSurf << "\n";

    gsMatrix<> para  = inputSurf->parameterRange();
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    //the parameter values for the fitting
    gsMatrix<> inputParams = uniformPointGrid(c0,c1, nd);
    //gsInfo<< "Parameter values used for fitting: "<<"\n"<<inputParams<<"\n";
    //the evaluated values of the original surface
    gsMatrix<> inputValues = inputSurf->eval(inputParams);
    //gsInfo<<"Original surface points: "<<"\n"<< *inputValues<<"\n";

    // Commented, because potential error messages do not get displayed in CDash (threshold of 1024 bytes).
    //gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsInfo<<"Creating the initial basis for fitting"<<"\n";
    //gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    int deg_x = 3;
    int deg_y = 3;
    int intKnots = 2;
    gsKnotVector<> KnotVector_u (0, 1, intKnots, deg_x+1, 1);
    gsKnotVector<> KnotVector_v (0, 1, intKnots, deg_y+1, 1);
    gsInfo<<"Knot Vector"<<KnotVector_u<<"\n";
    gsInfo<<"Knot Vector"<<KnotVector_v<<"\n";

    gsTensorBSplineBasis<2> TensorBasis( KnotVector_u, KnotVector_v );

    int nLevels =5;
    //helps to create the refinement
    gsTHBSplineBasis<2>  THB (TensorBasis) ;

    gsInfo<<"Basis has degree "<< deg_x <<" and "<< deg_y <<". The number of levels is "<< nLevels<<"\n";
    gsInfo<<"The tree has "<< THB.tree().size() << " nodes.\n" << "\n";

    //gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    gsInfo<<"Creating the THB fitting object"<<"\n";
    //gsInfo<<"//////////////////////////////////////////////////////////"<<"\n";
    std::vector<unsigned> ext;
    ext.push_back(2);
    ext.push_back(2);

    std::vector<index_t> constrainedBoxes;
    constrainedBoxes.push_back(1);
    constrainedBoxes.push_back(0);
    constrainedBoxes.push_back(0);
    constrainedBoxes.push_back(4);
    constrainedBoxes.push_back(4);

    plot=true;

    gsHFittingLvlConstrained<2, real_t> ref(inputParams, inputValues , THB, 0.01, ext, lambda, constrainedBoxes);
    gsMatrix<> results(numIter,6);
    results.setZero();

    for(int i = 0; i < numIter; i++) // We control the iterations explicitly.
    {
        ref.iterativeRefine(1,0,err_threshold); // We should rather call nextIteration.
        gsGeometry<> * test;
        test = ref.result();
        gsFileData<> newdata;
        newdata << *test ;
        gsTHBSpline<2>  * resultSurf = static_cast< gsTHBSpline<2>  *> (test);


        //results(i,0) = hbs1->basis().size();
        std::vector<real_t> errors;
        ref.get_Error(errors, err_type);
        real_t error;
        ref.computeApproxError(error, 0);
        real_t min = 1000000;
        real_t max = -1000000;
        for(unsigned int j =0; j < errors.size();j++)
        {
            if(errors[j]>max)
            {
                max = errors[j];
            }
            if(errors[j]<min)
            {
                min = errors[j];
            }
        }
        results(i,0) = resultSurf->basis().size();
        results(i,1) = resultSurf->basis().maxLevel();
        results(i,2) = min;
        results(i,3) = max;
        results(i,5) = error;

        real_t num = 0;
        for(unsigned int j = 0; j < errors.size(); j++)
        {
            if(errors[j]< 0.00001){
                num++;
            }
        }
        results(i,4) = (num*100)/errors.size();
        gsInfo<<"iteration: "<<i<<"\n";

        if(plot)
        {
            std::string mfn("gsThbs_adapt_face_1");
            std::stringstream ss;
            ss<<mfn<<i;
            newdata.dump(ss.str());
            gsWriteParaview( *test , ss.str(), np);
            ss<<".xml";
            std::string mfn1("gsThbs_cnet_adapt_face_1");
            std::stringstream ss1;
            ss1<<mfn1<<i;
            gsWriteParaview( *resultSurf, ss1.str(), np, false, true);
            std::string mfn2("gsThbs_error_adapt_face_1");
            std::stringstream ss2;
            ss2<<mfn2<<i;
            gsMatrix<> hbs1_eval =resultSurf->eval(inputParams);
            plot_errors(inputValues, hbs1_eval, errors, ss2.str());
        }
    }
    //gsInfo<<results<<"\n";
    gsInfo.setf(std::ios::scientific);
    gsInfo.setf(std::ios::scientific);
    for(int i = 0; i < results.rows(); i++)
    {
        gsInfo<< i+1 <<" & " << cast<real_t,int>(results(i,0))<<" & "<< cast<real_t,int>(results(i,1))<<" & "<<results(i,2)<<" & "<<results(i,3)<<" & "<<results(i,4)<<" & "<<results(i,5)<<"\n";
    }
    gsInfo<<"\n";
    gsInfo<<"error type "<< err_type<<"\n";
    gsInfo<<"extension "<< ext[0]<<" "<<ext[1]<<"\n";
    gsInfo<<"Basis has degree "<< deg_x <<" and "<< deg_y <<". The number of levels is "<< nLevels<<"\n";
    gsInfo<<"The tree has "<< THB.tree().size() << " nodes.\n" << "\n";

    gsMatrix<index_t> ll;
    gsMatrix<index_t> ur;
    gsVector<index_t> lvl;
    THB.tree().getBoxes(ll,ur,lvl);
    if( defaultExample )
        result = checkResults(results); // && checkBoxes(ll,ur,lvl);

//    gsInfo << "BOXES:" << ll << ";\n" << ur << ";\n" << lvl << "\n";


    gsInfo<<"\n"<<"Finished"<<"\n";

    if( result )
        return 0;
    else
        return -1;
}
