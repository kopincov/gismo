/** @file gsMultipatchFitting.cpp

    @brief Constructs a planar or volumetric multi-patch from sets of
    boundary points.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Sajavicius, J. Speh
*/

#include <gismo.h>
#include <gismo_dev.h>

#include "gsMotorUtils.h"

#include "gsMultipatchFitting.h"
#include "gsParameterLines.h"

using namespace gismo;

gsMatrix<> solve(const gsSparseMatrix<>& A,
                 const gsMatrix<>& B);

gsMultiPatch<> fitting(const gsMultiBasis<>& basis,
                       const std::vector< gsMatrix<> >& params,
                       const std::vector< gsMatrix<> >& points,
                       const real_t lambda);

gsMultiPatch<> fitting(const gsMultiBasis<>& basis,
                       const std::vector< gsMatrix<> >& params,
                       const std::vector< gsMatrix<> >& points,
                       const real_t lambda,
                       const std::vector< gsMatrix<> >& constrParams,
                       const std::vector< gsMatrix<> >& constrPoints,
                       const real_t penalty);

// initializes multi basis with mp.nPatches() number of basis
// each basis is a set of bernstein polynomials of degree degree
template <unsigned d>
void makeMultiBasis(gsMultiBasis<>& mb,
                    const gsMultiPatch<>& mp,
                    const int degree)
{
    gsMultiBasis<>::BasisContainer bases;

    std::vector< gsKnotVector<> > knotVectors;
    gsKnotVector<> kv(0.0, 1.0, 11, degree + 1);

    for (unsigned i = 0; i < d; i++)
    {
        knotVectors.push_back(kv);
    }

    for (unsigned i = 0; i != mp.nPatches(); i++)
    {
        gsTensorBSplineBasis<d> tenBasis(knotVectors);
        bases.push_back(new gsTHBSplineBasis<d>(tenBasis));
    }

    std::vector< patchSide > boundaries = mp.boundaries();
    std::vector< boundaryInterface > interfaces = mp.interfaces();

    gsBoxTopology topology(d, mp.nPatches(), boundaries, interfaces);

    gsMultiBasis<> m(bases, topology);

    mb = m;
}

void getPoints(std::vector< gsMatrix<> >& params,
               std::vector< gsMatrix<> >& points,
               const int n,
               const std::string& filename);

int main(int argc, char* argv[])
{
    // assumption point file have inside parameters and points put
    // inside one after another
    std::string pointsFile(MOTOR_DATA_DIR "jku/butterfly-points.xml");
    std::string constrPointsFile("constraint_points_and_parameters.xml");
    std::string multipatchFile(MOTOR_DATA_DIR "jku/butterfly-template.xml");
    std::string output("mp_fitting");
    int degree = 2;
    real_t lambda = 1e-9;
    real_t penalty = 0;
    int numAdaptiveRefIter = 5;
    int numUniformRefIter = 0;
    real_t threshold = 1e-4;
    int extension = 1;

    gsCmdLine cmd("Constructs a planar or volumetric multi-patch from sets of boundary points");
    cmd.addString("o", "output", "Output file", output);
    cmd.addInt("r", "numRefIter", "Number of adaptive (local) refinement procedure iterations", numAdaptiveRefIter);
    cmd.addInt("", "numURefIter", "Number of uniform (global) refinement procedure iterations", numUniformRefIter);
    cmd.addReal("t", "threshold", "Error threshold", threshold);
    cmd.addInt("x", "extension", "Extension of the refinement", extension);
    cmd.addReal("c", "penalty", "Penalty parameter (for boundary constraints)", penalty);
    cmd.addReal("l", "lambda", "Smoothing parameter", lambda);
    cmd.addInt("d", "degree", "Degree", degree);
    cmd.addString("M", "multipatchFile", "File containing multi-patch (required for topology data)", multipatchFile);
    cmd.addString("C", "constrPointsFile", "File containing boundary constraint points and parameters", constrPointsFile);
    cmd.addString("P", "pointsFile", "File containing fitting points and parameters", pointsFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "---------------------------------------------------------\n\n"
           << "Input Arguments: \n\n"
           << "points file:            " << pointsFile << "\n"
           << "constraint points file: " << constrPointsFile << "\n"
           << "multipatch file:        " << multipatchFile << "\n"
           << "output:                 " << output << "\n"
           << "lambda:                 " << lambda << "\n"
           << "penalty:                " << penalty << "\n"
           << "numRefIter:             " << numAdaptiveRefIter << "\n"
           << "threshold:              " << threshold << "\n"
           << "extension:              " << extension << "\n"
           << "---------------------------------------------------------\n"
           << std::endl;
    gsFileData<> fd(multipatchFile);
    gsMultiPatch<>* mp = fd.getAnyFirst< gsMultiPatch<> >().release();

    std::vector< gsMatrix<> > params;
    std::vector< gsMatrix<> > points;

    getPoints(params, points, static_cast<int>(mp->nPatches()), pointsFile);
    gsWriteParaviewPoints(stdVectorToMatrix(points), "test_points");
    // Dimension
    unsigned d = points[0].rows();

    std::vector< gsMatrix<> > constrParams;
    std::vector< gsMatrix<> > constrPoints;

    if (0 < penalty)
    {
        gsFileData<> fdTmp(constrPointsFile);
        int numBoundaries = fdTmp.count< gsMatrix<> >() / 2;
        getPoints(constrParams, constrPoints, numBoundaries, constrPointsFile);
    }

    gsMultiBasis<> mb;

    if (d == 2)
    {
        makeMultiBasis<2>(mb, *mp, degree);
    }
    else if (d == 3)
    {
        makeMultiBasis<3>(mb, *mp, degree);
    }
    gsInfo << "mb:\n" << mb << std::endl;

    for (int iter = 1; iter != numUniformRefIter+1; iter++)
    {
        mb.uniformRefine();
    }

    std::vector<int> ext(3, extension);

    for (int iter = 1; iter != numAdaptiveRefIter+1; iter++)
    {
        std::cout << "Iteration: " << iter << " / " << numAdaptiveRefIter << std::endl;

        gsMultiPatch<> solution;
        if (0 < penalty)
        {
            solution = fitting(mb, params, points, lambda, constrParams, constrPoints, penalty);
        }
        else
        {
            solution = fitting(mb, params, points, lambda);
        }

        // errors
        std::vector< gsVector<> > errors;
        computeErrors(solution, params, points, errors);

        // analyse
        real_t maxError = analyseErrors(errors, threshold, 0);
        std::cout << "Max error: " << maxError << std::endl;

        // saving
        saveDataMeshAndGeometry(solution, output + "_" + util::to_string(iter));
        saveDataSolution(solution, output +  "_" + util::to_string(iter));
        writeParameterLines(solution, output +  "_PARAMETER_LINES_" + util::to_string(iter), 500, 7);

        if (maxError < threshold)
        {
            std::cout << "Max error is below threshold, aborting..." << std::endl;
            break;
        }

        if (iter != numAdaptiveRefIter)
        {
            if (d == 2)
            {
                refine<2>(mb, params, errors, threshold, ext);
            }
            else if (d == 3)
            {
                refine<3>(mb, params, errors, threshold, ext);
            }
            repairInterfaces(mb);
        }

    }

    return 0; 

}

gsMatrix<> solve(const gsSparseMatrix<>& A,
                 const gsMatrix<>& B)
{
    gsSparseSolver<>::BiCGSTABILUT solver(A);
    if ( solver.preconditioner().info() != Eigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
        return gsMatrix<>();
    }
    return solver.solve(B);
}

gsMultiPatch<> fitting(const gsMultiBasis<>& basis,
                       const std::vector< gsMatrix<> >& params,
                       const std::vector< gsMatrix<> >& points,
                       const real_t lambda)
{
    gsSparseMatrix<> A;
    gsMatrix<> B;

    fittingMatrix(basis, params, points, lambda, A, B);

    const gsMatrix<> X = solve(A, B);

    return constructSolution(basis, X);
}

gsMultiPatch<> fitting(const gsMultiBasis<>& basis,
                       const std::vector< gsMatrix<> >& params,
                       const std::vector< gsMatrix<> >& points,
                       const real_t lambda,
                       const std::vector< gsMatrix<> >& constrParams,
                       const std::vector< gsMatrix<> >& constrPoints,
                       const real_t penalty)
{
    gsSparseMatrix<> A;
    gsMatrix<> B;

    fittingMatrix(basis, params, points, lambda, constrParams, constrPoints, penalty, A, B);

    const gsMatrix<> X = solve(A, B);

    return constructSolution(basis, X);
}

void getPoints(std::vector< gsMatrix<> >& params,
               std::vector< gsMatrix<> >& points,
               const int n,
               const std::string& filename)
{
    gsFileData<> fd(filename);

    for (int i = 0; i != n; i++)
    {
        const int id_par = 2 * i;
        const int id_pts = 2 * i + 1;

        gsMatrix<>* par = fd.getId< gsMatrix<> >(id_par).release();
        gsMatrix<>* pts = fd.getId< gsMatrix<> >(id_pts).release();

        params.push_back(*par);
        points.push_back(*pts);
    }
}
