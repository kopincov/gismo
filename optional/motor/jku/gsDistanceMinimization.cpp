/** @file gsDistanceMinimization.cpp

    @brief ....

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh, S. Sajavicius
*/

#include <utility>

#include <gismo.h>
#include "gsMotorUtils.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

//void BasisToMatrix(const gsMatrix<>& basis, const gsMatrix<index_t> &actives, gsMatrix<>& result);

gsMatrix<> assembleMatrixTDM(const gsMatrix<>& basisValues,
                             const gsMatrix<index_t>& basisActives,
                             const gsMatrix<>& nrmls,
                             const gsMatrix<> &coefs);
gsMatrix<> assembleRHSTDM(const gsMatrix<>& basisValues,
                          const gsMatrix<index_t>& basisActives,
                          const gsMatrix<>& sampledPtsTargetCurve,
                          const gsMatrix<>& targetPoints,
                          const gsMatrix<>& nrmls,
                          const gsMatrix<>& coefs);
gsMatrix<> assembleMatrixTDMdirect(const gsMatrix<>& basisValues,
                                   const gsMatrix<index_t>& basisActives,
                                   const gsMatrix<>& nrmls,
                                   const gsMatrix<> &coefs);
gsMatrix<> assembleRHSTDMdirect(const gsMatrix<>& basisValues,
                                const gsMatrix<index_t>& basisActives,
                                const gsMatrix<>& sampledPtsTargetCurve,
                                const gsMatrix<>& targetPoints,
                                const gsMatrix<>& nrmls,
                                const gsMatrix<>& coefs);
gsMatrix<> assembleMatrixPDM(const gsMatrix<>& basisValues,
                             const gsMatrix<index_t>& basisActives,
                             const gsMatrix<>& sampledPtsTargetCurve,
                             const gsMatrix<>& targetPoints,
                             const gsMatrix<> &coefs);
gsMatrix<> assembleRHSPDM(const gsMatrix<>& basisValues,
                          const gsMatrix<index_t>& basisActives,
                          const gsMatrix<>& sampledPtsTargetCurve,
                          const gsMatrix<>& targetPoints,
                          const gsMatrix<>& coefs);
gsMatrix<> assembleMatrixPDMdirect(const gsMatrix<>& basisValues,
                             const gsMatrix<index_t>& basisActives,
                             const gsMatrix<> &coefs);
gsMatrix<> assembleRHSPDMdirect(const gsMatrix<>& basisValues,
                          const gsMatrix<index_t>& basisActives,
                          const gsMatrix<>& sampledPtsTargetCurve,
                          const gsMatrix<>& targetPoints,
                          const gsMatrix<>& coefs);
void applyTikhonovRegularization(real_t lambda,
                                 gsMatrix<>& A_mat);
void applySmoothing(real_t lambda,
                    gsMatrix<>& A_mat,
                    gsMatrix<>& rhs,
                    gsGeometry<>& map);
gsMatrix<> solveTest(const gsMatrix<>& A);
gsMatrix<> solveDisplacements(const gsMatrix<>& A,
                              const gsMatrix<>& b,
                              const bool debug = false,
                              const std::string debPrefix = "");
real_t computeFunctionalTDM(const gsMatrix<>& sampledPtsTargetCurve,
                            const gsMatrix<>& targetPoints,
                            const gsMatrix<>& nrmls);
real_t computeFunctionalTDM(const gsMatrix<>& sampledPtsTargetCurve,
                            const gsMatrix<>& templatePoints,
                            const gsGeometry<>& thbMap,
                            const gsMatrix<>& nrmls);
gsMatrix<> computeClosestPointsParameters(const gsGeometry<>& thbMap,
                                          const gsMatrix<>& sampledPtsTargetCurve,
                                          const gsMatrix<>& ptsTemplateCurve,
                                          const gsMatrix<>& sampleParameters);
gsMatrix<> normals(const gsMatrix<>& derivMapBasis,
                   const gsMatrix<index_t>& basisActives,
                   const gsMatrix<> &coefs,
                   const gsMatrix<>& derivTempBoundaryGeom);
void saveData(const gsGeometry<>& thbMap,
              const gsMultiPatch<>& tempBoundaryGeom,
              const std::string outputPrefix);
std::pair<real_t, real_t>
computePDMerror(const gsMatrix<>& sampledPtsTargetCurve,
                      const gsMatrix<>& templatePoints,
                      const gsGeometry<>& thbMap);
std::pair<real_t, real_t>
computeTDMerror(const gsMatrix<>& sampledPtsTargetCurve,
                const gsMatrix<>& templatePoints,
                const gsMatrix<>& derivTempBoundaryGeom,
                const gsGeometry<>& thbMap);

real_t computeMAXerror(const gsMatrix<>& sampledPtsTargetCurve,
                       const gsMatrix<>& templatePoints,
                       const gsGeometry<>& thbMap);

/**********************************************************************/

class gsDistanceMinimization
{
public:
    gsDistanceMinimization(gsGeometry<>& thbMap,
                           const gsMultiPatch<>& templateMultipatch,
                           const gsMultiPatch<>& targetMultipatch)
                           : m_thbMap(thbMap),
                             m_templateMultipatch(templateMultipatch),
                             m_targetMultipatch(targetMultipatch),
                             m_debug(false),
                             m_sampleTemplateParameters(),
                             //m_sampledPtsTemplateCurve(),
                             m_sampledPtsTargetCurve(),
                             m_sampledPtsTemplateMultipatch(), // NEW
                             m_sampledPtsTargetMultipatch(),   // NEW
                             //m_derivTempBoundaryGeom(), // Newly added
                             m_templatePoints(),
                             m_templatePointsMultipatch(),   // NEW
                             m_parametersMultipatch(),
                             m_nPatches(templateMultipatch.nPatches()),
                             m_indent("    "),
                             m_debStr(m_indent + "gsDistanceMinimization debug: ")
    {
    }

    // Sample numPoints points on the template geometry
    // and numPointsTarget points on the target geometry
    void computeTargetAndTemplatePoints(const int numPoints, const int numPointsTarget)
    {
        // Generating sample parameters
        m_sampleTemplateParameters = uniformParameters<2>(0.0, 1.0, numPoints);
        const gsMatrix<> sampleTemplateParameters = uniformParameters<2>(0.0, 1.0, numPointsTarget);

        for (std::size_t i = 0; i != m_nPatches ; i++)
        {
            // Sampling points on the template boundary curve
            const gsGeometry<>& templateBoundaryCurve = m_templateMultipatch.patch(i);
            gsMatrix<> ptsTemplateCurve = templateBoundaryCurve.eval(m_sampleTemplateParameters);
            m_sampledPtsTemplateMultipatch.push_back(ptsTemplateCurve);

            // Sampling points on the target boundary curve
            const gsGeometry<>& targetBoundaryCurve = m_targetMultipatch.patch(i);
            gsMatrix<> ptsTargetCurve = targetBoundaryCurve.eval(sampleTemplateParameters);
            m_sampledPtsTargetMultipatch.push_back(ptsTargetCurve);
        }
        m_sampledPtsTargetCurve = stdVectorToMatrix(m_sampledPtsTargetMultipatch);
    }

    // 
    void set_debug(const bool debug)
    {
        m_debug = debug;
    }

    //
    void computeParameters(const bool pcorrect)
    {
        m_parametersMultipatch.clear();
        for (std::size_t i = 0; i != m_nPatches; i++)
        {
            if (pcorrect)
            {
                const gsMatrix<> sampledPtsTemplateCurve = m_sampledPtsTemplateMultipatch.at(i);
                const gsMatrix<> sampledPtsTargetCurve   = m_sampledPtsTargetMultipatch.at(i);
                const gsMatrix<> parameters = computeClosestPointsParameters(m_thbMap, sampledPtsTargetCurve,
                                                                             sampledPtsTemplateCurve,
                                                                             m_sampleTemplateParameters);
                m_parametersMultipatch.push_back(parameters);
            }
            else
            {
                const int numPoints = m_sampledPtsTargetMultipatch.at(i).cols();
                const gsMatrix<> parameters = uniformParameters<2>(0.0, 1.0, numPoints);
                m_parametersMultipatch.push_back(parameters);
            }
        }
    }
  
    // Estimate accuracy of the initial solution (usually, the identity map)
    void estimateInitialSolution(const std::string& method)
    {
        std::vector< gsMatrix <> > targetPointsMultipatch = getTargetPointsMultipatch();
        m_templatePoints = stdVectorToMatrix(m_templatePointsMultipatch);

        computeDerivTempBoundaryGeom();

        typedef std::pair<real_t, real_t> pair_t;
        pair_t error;
        if (method== "tdm")
        {
            error = computeTDMerror(m_sampledPtsTargetCurve, m_templatePoints, m_derivTempBoundaryGeom, m_thbMap);
        }
        else if (method == "pdm")
        {
            error = computePDMerror(m_sampledPtsTargetCurve, m_templatePoints, m_thbMap);
        }
        gsInfo << m_indent << " Max error:    " << error.second << std::endl;
        gsInfo << m_indent << " L2 error:     " << error.first << std::endl;
        gsInfo << m_indent << " Eucl error:   " << computeMAXerror(m_sampledPtsTargetCurve, m_templatePoints, m_thbMap) << std::endl;

    }
    
    // Compute function with parameters given as an argument
    void compute(const int numIter,
                 const real_t lambdaTikhonov,
                 const real_t lambdaSmoothing,
                 const std::string& method,
                 const bool pcorrect,
                 const std::string& outputPrefix)
    {
        // Distance minimization (PDM or TDM) iterations
        for (int iter = 0; iter < numIter; iter++)
        {
            gsInfo << m_indent << method + " Iteration "
                   << iter + 1 << " / " << numIter << std::endl;

            // If method uses parameter correction, we recalculate parameters here
            if (pcorrect)
            {
                this->computeParameters(pcorrect);
            }
            std::vector< gsMatrix <> > targetPointsMultipatch = getTargetPointsMultipatch();
            m_templatePoints = stdVectorToMatrix(m_templatePointsMultipatch);
            gsMatrix<> targetPoints = stdVectorToMatrix(targetPointsMultipatch);
            if (m_debug)
            {
                std::string tempPtsName = makeDebugName(outputPrefix, "templatePoints", iter);
                std::string targPtsName = makeDebugName(outputPrefix, "targetPoints", iter);

                gsInfo << m_debStr << "saving: " << tempPtsName << " & "
                       << targPtsName << std::endl;

                gsWriteParaviewPoints(m_templatePoints, tempPtsName);
                gsWriteParaviewPoints(targetPoints, targPtsName);

            }

            // THB-spline map values on the template points
            gsMatrix<> basisValues;
            m_thbMap.basis().eval_into(m_templatePoints, basisValues);
            gsMatrix<index_t> basisActives;
            m_thbMap.basis().active_into(m_templatePoints, basisActives);

            computeDerivTempBoundaryGeom();

            // Derivatives of basis functions
            // tau_1_u(t_1) tau_1_u(t_2) ... tau_1_u(t_M)
            // tau_1_v(t_1) tau_1_v(t_2) ... tau_1_v(t_M)
            // tau_2_u(t_1) tau_2_u(t_2) ... tau_2_u(t_M)
            // tau_2_v(t_1) tau_2_v(t_2) ... tau_2_v(t_M)
            // ...          ...          ... ...
            // tau_N_u(t_1) tau_N_u(t_2) ... tau_N_u(t_M)
            // tau_N_v(t_1) tau_N_v(t_2) ... tau_N_v(t_M)
            gsMatrix<> derivMapBasis;
            m_thbMap.basis().deriv_into(m_templatePoints, derivMapBasis);

            gsMatrix<>& coefs = m_thbMap.coefs();
            const int M = derivMapBasis.cols();
            gsMatrix<> nrmls(2, M);
            nrmls.setZero();
            if (method == "tdm")
            {
                nrmls = normals(derivMapBasis, basisActives, coefs, m_derivTempBoundaryGeom);
            }

            int N = coefs.rows();
            gsMatrix<> A(2*N, 2*N);
            A.setZero();
            gsMatrix<> b(2*N, 1);
            b.setZero();
            if (method == "tdm")
            {
                A = assembleMatrixTDM(basisValues, basisActives, nrmls, coefs);
                b = assembleRHSTDM(basisValues, basisActives, m_sampledPtsTargetCurve, targetPoints, nrmls, coefs);
            }
            else if (method == "pdm")
            {
                A = assembleMatrixPDMdirect(basisValues, basisActives, coefs);
                b = assembleRHSPDMdirect(basisValues, basisActives, m_sampledPtsTargetCurve, targetPoints, coefs);
            }
            else
            {
                gsInfo << "Wrong value of --method option (should be pdm or tdm and it is case sensitive)" << std::endl;
            }

            if (lambdaTikhonov != 0)
            {
                applyTikhonovRegularization(lambdaTikhonov, A);
            }

            if (lambdaSmoothing > 0)
            {
                applySmoothing(lambdaSmoothing, A, b, m_thbMap);
            }

            gsMatrix<> delta = solveDisplacements(A, b);

            if (m_debug)
            {
                gsInfo << m_debStr << "Displacements max norm: "
                       << delta.lpNorm<Eigen::Infinity>() << std::endl;
            }

            coefs = coefs + delta;
            m_thbMap.setCoefs(coefs);

            typedef std::pair<real_t, real_t> pair_t;
            pair_t error;
            if (method == "tdm")
            {
                error = computeTDMerror(m_sampledPtsTargetCurve, m_templatePoints, m_derivTempBoundaryGeom, m_thbMap);
            }
            else if (method == "pdm")
            {
                error = computePDMerror(m_sampledPtsTargetCurve, m_templatePoints, m_thbMap);
            }
            gsInfo << m_indent << " Max error:    " << error.second << std::endl;
            gsInfo << m_indent << " L2 error:     " << error.first << std::endl;
            gsInfo << m_indent << " Eucl error:   " << computeMAXerror(m_sampledPtsTargetCurve, m_templatePoints, m_thbMap) << std::endl;
            
            saveIterationOutput(outputPrefix, iter);

            // Giving different names for files with (possibly) final results
            writeMap(m_thbMap, 0,"results_dm_final");

        }

    }

    //
    void
    computeDerivTempBoundaryGeom()
    {
        // Derivatives of template boundary curve (geometry) at parameter values
        // gamma_1'(t_1) gamma_1'(t_2) ... gamma_1'(t_M)
        // gamma_2'(t_1) gamma_2'(t_2) ... gamma_2'(t_M)
        std::vector< gsMatrix<> > derivTempBoundaryGeomMultipatch;
        for (size_t i = 0; i != m_nPatches; i++)
        {
            gsMatrix<> derivTempBoundaryGeomTmp;
            m_templateMultipatch.patch(i).deriv_into(m_parametersMultipatch.at(i), derivTempBoundaryGeomTmp);
            derivTempBoundaryGeomMultipatch.push_back(derivTempBoundaryGeomTmp);
        }
        m_derivTempBoundaryGeom = stdVectorToMatrix(derivTempBoundaryGeomMultipatch);
    }


    // 
    const gsMatrix<>& getTemplatePoints() const
    {
        return m_templatePoints;
    }

    //
    const gsMatrix<>& getSampledTargetPoints() const
    {
        return m_sampledPtsTargetCurve;
    }
    
    //
    const gsMatrix<>& getSampleParameters() const
    {
        return m_sampleTemplateParameters;
    }
    
    //
    const gsMatrix<>& getDerivTempBoundaryGeom() const
    {
        return m_derivTempBoundaryGeom;
    }

    
private:

    //
    void saveIterationOutput(const std::string& outputPrefix,
                             const int iter)
    {
        const std::string prefix = outputPrefix + "Iter" + util::to_string(iter);
        gsMultiPatch<> mappedMp = mapMultipatch(m_thbMap, m_templateMultipatch);

        writeMultipatch(mappedMp, prefix + "_points");
        writeMap(m_thbMap, prefix);
        writeKnots(m_thbMap, prefix);
    }

    //
    std::vector< gsMatrix <> > getTargetPointsMultipatch()
    {
        std::vector< gsMatrix <> > targetPointsMultipatch;
        m_templatePointsMultipatch.clear();
        //gsInfo << "m_thbMap.coefs():\n" << m_thbMap.coefs() << std::endl;
        for (std::size_t i = 0; i != m_nPatches; i++)
        {
            // Calculate points on the template boundary curve
            gsMatrix<> templatePoints;
            m_templateMultipatch.patch(i).eval_into(m_parametersMultipatch.at(i), templatePoints);
            m_templatePointsMultipatch.push_back(templatePoints);

            // Use THB-splines map to map these points (see previous step) on the target geometry
            gsMatrix<> targetPointsTmp;
            m_thbMap.eval_into(templatePoints, targetPointsTmp);
            targetPointsMultipatch.push_back(targetPointsTmp);

            //gsInfo << "ptsTemp / ptsTempParam (" << i << "):\n" << templatePoints << std::endl;
            //gsInfo << "ptsMappedTemp / ptsMappedTempParam (" << i << "):\n" << targetPointsTmp << std::endl;
        }

        return targetPointsMultipatch;
    }

    //
    std::string makeDebugName(const std::string& outputPrefix,
                              const std::string& name)
    {
        std::string newName = "debug_" + outputPrefix + name;
        return newName;
    }

    //
    std::string makeDebugName(const std::string& outputPrefix,
                              const std::string& name,
                              const int iteration)
    {
        return makeDebugName(outputPrefix + "Iter" + util::to_string(iteration), name);
    }
    
    // documentation for this variable
    gsGeometry<>& m_thbMap;

    const gsMultiPatch<>& m_templateMultipatch;

    const gsMultiPatch<>& m_targetMultipatch;

    bool m_debug;

    gsMatrix<> m_sampleTemplateParameters;
    
    //gsMatrix<> m_sampledPtsTemplateCurve;
    
    gsMatrix<> m_sampledPtsTargetCurve; // P_k

    std::vector< gsMatrix<> > m_sampledPtsTemplateMultipatch;

    std::vector< gsMatrix<> > m_sampledPtsTargetMultipatch;
    
    gsMatrix<> m_derivTempBoundaryGeom;

    gsMatrix<> m_templatePoints;

    std::vector< gsMatrix<> > m_templatePointsMultipatch;

    std::vector< gsMatrix<> > m_parametersMultipatch;

    std::size_t m_nPatches;
    
    const std::string m_indent;

    std::string m_debStr;
    
};

std::vector<real_t> computeErrors(const gsGeometry<>& thbMap,
                                  const gsMatrix<>& params,
                                  const gsMatrix<>& points)
{
    std::vector<real_t> pointErrors;
    pointErrors.reserve(params.cols());
        
    gsMatrix<> values;
    thbMap.eval_into(params, values);

    for (index_t i = 1; i < points.cols(); i++)
    {
        const real_t err = (points.col(i) - values.col(i)).norm();
        pointErrors.push_back(err);
    }

    return pointErrors;
}

void append(std::vector<index_t>& boxes,
            const gsVector<index_t>& box)
{
    for (index_t col = 0; col != box.rows(); col++)
    {
        boxes.push_back(box[col]);
    }
}

bool isCellAlreadyInserted(const gsVector<index_t, 2>& a_cell,
                           const std::vector<index_t>& cells)
{
    for (std::size_t i = 0; i != cells.size(); i += a_cell.rows())
    {
        int commonEntries = 0;
        for (index_t col = 0; col != a_cell.rows(); col++)
        {
            if (cells[i + col] == a_cell[col])
            {
                commonEntries++;
            }

        }

        if (commonEntries == a_cell.rows())
        {
            return true;
        }
    }

    return false;
}

void appendBox(const gsTHBSplineBasis<2>& basis,
               std::vector<index_t>& boxes,
               std::vector<index_t>& cells,
               const gsVector<>& param,
               const std::vector<index_t>& ext)
{

    const int maxLvl = basis.maxLevel();
    const gsTensorBSplineBasis<2>& tBasis = *(basis.getBases()[maxLvl]);

    gsVector<index_t, 2> a_cell;

    for (index_t dim = 0; dim != 2; dim++)
    {
        const gsKnotVector<>& kv = tBasis.component(dim).knots();
        a_cell(dim) = static_cast<index_t>(kv.uFind(param(dim)).uIndex());
    }

    if (!isCellAlreadyInserted(a_cell, cells))
    {
        append(cells, a_cell);

        // get level of a cell
        gsVector<index_t, 2> a_cell_upp = a_cell + gsVector<index_t, 2>::Ones();
        const int cell_lvl = basis.tree().query3(a_cell, a_cell_upp, maxLvl) + 1;

        // get the box
        gsVector<index_t> box(2 * 2 + 1);
        box[0] = cell_lvl;

        for (index_t dim = 0; dim != 2; dim++)
        {
            const index_t numBreaks = basis.numBreaks(cell_lvl, dim) - 1 ;

            index_t lowIndex = 0;
            if (cell_lvl < maxLvl)
            {
                const index_t shift = maxLvl - cell_lvl;
                lowIndex = (a_cell(dim) >> shift);
            }
            else
            {
                const index_t shift = cell_lvl - maxLvl;
                lowIndex = (a_cell(dim) << shift);
            }

            index_t low = ( (lowIndex > ext[dim]) ? (lowIndex - ext[dim]) : 0 );
            index_t upp = ( (lowIndex + ext[dim] + 1 < numBreaks) ?
                             (lowIndex + ext[dim] + 1) : numBreaks );

            box[1 + dim    ] = low;
            box[1 + 2 + dim] = upp;
        }

        append(boxes, box);
    }
}

std::vector<index_t> getBoxes(const gsTHBSplineBasis<2>& basis,
                               const gsMatrix<>& params,
                               const std::vector<real_t>& pointErrors,
                               const real_t threshold,
                               const std::vector<index_t>& ext)
{
    // look at gsHFitting for more info

    std::vector<index_t> cells;
    std::vector<index_t> boxes;

    for (std::size_t index = 0; index != pointErrors.size(); index++)
    {
        if (threshold <= pointErrors[index])
        {
            appendBox(basis, boxes, cells, params.col(index), ext);
        }
    }

    return boxes;
}


bool uniformRefine(gsTHBSpline<2>& thbMap)
{
    thbMap.basis().uniformRefine_withCoefs(thbMap.coefs());

    return true;
}

bool refine(gsTHBSpline<2>& thbMap,
            const gsMatrix<>& params,
            const gsMatrix<>& points,
            const int extension,
            const real_t threshold)
{
    
    std::vector<index_t> ext(2, extension);

    std::vector<real_t> pointErrors = computeErrors(thbMap, params, points);
    const real_t maxError = *std::max_element(pointErrors.begin(), pointErrors.end());

    if (threshold <= maxError)
    {
        std::vector<index_t> boxes = getBoxes(thbMap.basis(), params,
                                               pointErrors, threshold, ext);

        thbMap.basis().refineElements_withCoefs(thbMap.coefs(), boxes);

        return true;
    }
    else
    {
        return false;
    }
}

/*************** Main function ***************/

int main(int argc, char *argv[])
{

    // Options with default values
    std::string tempFile(MOTOR_DATA_DIR "jku/tdm_template_multipatch.xml");
    std::string geomFile(MOTOR_DATA_DIR "jku/tdm_target_multipatch.xml");
    std::string mapFile(MOTOR_DATA_DIR "jku/tdm_thb_map.xml");
    std::string method("tdm");

    std::string outputPrefix("results");
    real_t lambdaTikhonov = 1.0e-6;
    real_t lambdaSmoothing = 1.0e-12;
    int numPoints        = 1000;
    int numPointsTarget  = 100;
    int numIter          = 5;
    int numRefIter       = 5;
    bool adaptiveRef     = false;
    real_t threshold     = 1e-5;
    int extension = 1;

    bool pcorrect = false;
    bool debug    = false;
    
    gsCmdLine cmd("Tangent distance minimization (TDM)");
    cmd.addString("T", "templateFile", "Name of template geometry multipatch file (input)", tempFile);
    cmd.addString("G", "geometryFile", "Name of target geometry multipatch file (input)", geomFile);
    cmd.addString("M", "mapFile", "Name of initial map file (input)", mapFile);
    cmd.addString("A", "method", "Method of fitting (expected values: tdm, pdm)", method);
    cmd.addSwitch("p", "pcorrect", "Apply parameter correction (optimization)", pcorrect);
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);
    cmd.addReal("l", "lambdaT", "Regularization parameter", lambdaTikhonov);
    cmd.addReal("L", "lambdaS", "Smoothing parameter", lambdaSmoothing);
    cmd.addInt("n", "npts", "Number of sampling points (for sampling)", numPoints);
    cmd.addInt("m", "nptst", "Number of points fitting", numPointsTarget);
    cmd.addInt("i", "niter", "Number of iterations", numIter);
    cmd.addInt("r", "nriter", "Number of refinement iterations", numRefIter);
    cmd.addSwitch("a", "adaptiveRef", "Adaptive refinement", adaptiveRef);
    cmd.addReal("t", "threshold", "Threshold", threshold);
    cmd.addInt("x", "extension", "The extension of the refinement", extension);
    cmd.addSwitch("debug", "If debug, it outputs more data.", debug);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments: \n"
           << "Template file:     " << tempFile << "\n"
           << "Geometry file:     " << geomFile << "\n"
           << "Map file:          " << mapFile << "\n"
           << "Method:            " << method << "\n"
           << "Par. correction:   " << pcorrect << "\n"
           << "output:      Compute      " << outputPrefix << "\n"
           << "lambdaTikhonov:    " << lambdaTikhonov << "\n"
           << "lambdaSmoothing:   " << lambdaSmoothing << "\n"
           << "numPoints:         " << numPoints << "\n"       // Large number
           << "numPointsTarget:   " << numPointsTarget << "\n" // numPointsTarget << numPoints
           << "numIter:           " << numIter << "\n"
           << "numRefIter:        " << numRefIter << "\n"
           << "adaptiveRef:       " << adaptiveRef << "\n"
           << "threshold:         " << threshold << "\n"
           << "extension:         " << extension << "\n"
           << "Debug:             " << debug << "\n"
           << "--------------------------------------------------\n" << std::endl;

    const gsMultiPatch<>* templateM = readMultipatch(tempFile);
    const gsMultiPatch<>& templateMultipatch = *templateM;
    const gsMultiPatch<>* targetM   = readMultipatch(geomFile);
    const gsMultiPatch<>& targetMultipatch = *targetM;
    gsGeometry<>* mapPtr = readGeometry(mapFile);
    gsTHBSpline<2>& thbMap = static_cast<gsTHBSpline<2>&>(*mapPtr);
    
    const int templateNumPatches = templateMultipatch.nPatches();
    const int targetNumPatches   = targetMultipatch.nPatches();

    if ( templateNumPatches != targetNumPatches)
    {
        gsInfo << "Template geometry and target geometry multipatch files should have the same numbers of patches (curves).\n";
        return 1;
    }

    gsWriteParaview(targetMultipatch, outputPrefix + "_target_multipatch");
    saveData(thbMap, templateMultipatch, outputPrefix + "_initial");
    
    gsDistanceMinimization distanceMinimizer(thbMap, templateMultipatch, targetMultipatch);
    distanceMinimizer.set_debug(debug);
    distanceMinimizer.computeTargetAndTemplatePoints(numPoints, numPointsTarget);
    
    // We precalculate parameters(we need this for estimateInitialSolution())
    // If pcorrect == true, we recalculate parameters in each PDM/TDM iteration, if no - we keep them fixed
    distanceMinimizer.computeParameters(pcorrect);
    // Estimating the accuracy of the initial solution
    distanceMinimizer.estimateInitialSolution(method);

    // Computing in the initial spline space (before refinement)
    const std::string prefix = outputPrefix + "Ref" + util::to_string(0);
    gsInfo << "Initial spline space DOFs: " << thbMap.coefsSize() << std::endl;
    distanceMinimizer.compute(numIter, lambdaTikhonov, lambdaSmoothing, method, pcorrect, prefix);

    for (int r = 0; r != numRefIter; r++)
    {      
        // Refinement
        gsInfo << "Refinement: " << r + 1 << " / " << numRefIter << std::endl;       
        bool refined = false;
        if(adaptiveRef == true) // Adaptive refinement
        {
            refined = refine(thbMap, distanceMinimizer.getTemplatePoints(), distanceMinimizer.getSampledTargetPoints(),
                                  extension, threshold);
            gsInfo << "Adaptive refinement DOFs: " << thbMap.coefsSize() << std::endl;
        }
        else // Uniform refinement
        {
            refined = uniformRefine(thbMap);
            gsInfo << "Uniform refinement DOFs: " << thbMap.coefsSize() << std::endl;
        }
        if (!refined)
        {
            gsInfo << "Error is below the threshold (" << threshold << ")" << std::endl;
        } 
        gsInfo << std::endl;

        const std::string prefix = outputPrefix + "Ref" + util::to_string(r+1);
        distanceMinimizer.compute(numIter, lambdaTikhonov, lambdaSmoothing, method, pcorrect, prefix);

    }
    
    return 0;
}

/******************************/


/*************** Tangent distance minimization (TDM) **************/

gsMatrix<> assembleMatrixTDM(const gsMatrix<>& basisValues,
                             const gsMatrix<index_t>& basisActives,
                             const gsMatrix<>& nrmls,
                             const gsMatrix<>& coefs)
{
    const int M = basisValues.cols();
    const int N = coefs.rows();
    gsMatrix<> A(2*N, 2*N);
    A.setZero();

    //gsInfo << "nrmls (used in assembleMatrixTDM):\n" << nrmls << std::endl;

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {
        for (index_t idActive1 = 0; idActive1 < basisActives.rows(); idActive1++) // Works as l
        {
            const int idBasis1 = basisActives(idActive1, idPoint);

            for (index_t idActive2 = 0; idActive2 < basisActives.rows(); idActive2++) // Works as i
            {
                const int idBasis2 = basisActives(idActive2, idPoint);

                if (    ((idBasis1 > 0) || ((idBasis1 == 0) && (idActive1 == 0)))
                     && ((idBasis2 > 0) || ((idBasis2 == 0) && (idActive2 == 0))) )
                {
                    const real_t basisValue1 = basisValues(idActive1, idPoint);
                    const real_t basisValue2 = basisValues(idActive2, idPoint);
                    const real_t nrml0 = nrmls(0, idPoint);
                    const real_t nrml1 = nrmls(1, idPoint);

                    A(2*idBasis1,   2*idBasis2)   += basisValue2 * nrml0 * basisValue1 * nrml0; //
                    A(2*idBasis1,   2*idBasis2+1) += basisValue2 * nrml1 * basisValue1 * nrml0; //

                    A(2*idBasis1+1, 2*idBasis2)   += basisValue2 * nrml0 * basisValue1 * nrml1; //
                    A(2*idBasis1+1, 2*idBasis2+1) += basisValue2 * nrml1 * basisValue1 * nrml1; //
                }
            }
        }      
    }

    return A;
}

gsMatrix<> assembleRHSTDM(const gsMatrix<>& basisValues,
                          const gsMatrix<index_t>& basisActives,
                          const gsMatrix<>& sampledPtsTargetCurve,
                          const gsMatrix<>& targetPoints,
                          const gsMatrix<>& nrmls,
                          const gsMatrix<>& coefs)
{
    const int M = basisValues.cols();
    const int N = coefs.rows();
    gsMatrix<> rhs(2*N, 1);
    rhs.setZero();

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {

        const gsVector<> P = sampledPtsTargetCurve.col(idPoint);
        const gsVector<> targetPoint = targetPoints.col(idPoint);
        const gsVector<> nrml = nrmls.col(idPoint);
        const real_t tmp = nrml.dot(P - targetPoint);

        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
        {
            int idBasis = basisActives(idActive, idPoint);
            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) )
            {
                const real_t basisValue = basisValues(idActive, idPoint);
                rhs(2*idBasis,   0) += tmp*basisValue*nrml(0);
                rhs(2*idBasis+1, 0) += tmp*basisValue*nrml(1);
            }
        }
    }
    return rhs;
}

/*************** Tangent distance minimization (TDM) - direct approach (without Gauss--Newton) **************/

gsMatrix<> assembleMatrixTDMdirect(const gsMatrix<>& basisValues,
                                   const gsMatrix<index_t>& basisActives,
                                   const gsMatrix<>& nrmls,
                                   const gsMatrix<>& coefs)
{
    const int M = basisValues.cols();
    const int N = coefs.rows();
    gsMatrix<> A(2*N, 2*N);
    A.setZero();

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {
        for (index_t idActive1 = 0; idActive1 < basisActives.rows(); idActive1++) // Works as l
        {
            const int idBasis1 = basisActives(idActive1, idPoint);

            for (index_t idActive2 = 0; idActive2 < basisActives.rows(); idActive2++) // Works as i
            {
                const int idBasis2 = basisActives(idActive2, idPoint);

                if (    ((idBasis1 > 0) || ((idBasis1 == 0) && (idActive1 == 0)))
                     && ((idBasis2 > 0) || ((idBasis2 == 0) && (idActive2 == 0))) )
                {
                    const real_t basisValue1 = basisValues(idActive1, idPoint);
                    const real_t basisValue2 = basisValues(idActive2, idPoint);
                    const real_t nrml0 = nrmls(0, idPoint);
                    const real_t nrml1 = nrmls(1, idPoint);

                    A(2*idBasis1,   2*idBasis2)   += basisValue2 * nrml0 * basisValue1 * nrml0; //
                    A(2*idBasis1,   2*idBasis2+1) += basisValue2 * nrml1 * basisValue1 * nrml0; //

                    A(2*idBasis1+1, 2*idBasis2)   += basisValue2 * nrml0 * basisValue1 * nrml1; //
                    A(2*idBasis1+1, 2*idBasis2+1) += basisValue2 * nrml1 * basisValue1 * nrml1; //
                }
            }
        }
    }

    return A;
}

gsMatrix<> assembleRHSTDMdirect(const gsMatrix<>& basisValues,
                                const gsMatrix<index_t>& basisActives,
                                const gsMatrix<>& sampledPtsTargetCurve,
                                const gsMatrix<>& targetPoints,
                                const gsMatrix<>& nrmls,
                                const gsMatrix<>& coefs)
{
    const int M = basisValues.cols();
    const int N = coefs.rows();
    gsMatrix<> rhs(2*N, 1);
    rhs.setZero();

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {

        const gsVector<> P = sampledPtsTargetCurve.col(idPoint);
        const gsVector<> targetPoint = targetPoints.col(idPoint);
        const gsVector<> nrml = nrmls.col(idPoint);
        const real_t tmp = nrml.dot(P - targetPoint);

        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
        {
            int idBasis = basisActives(idActive, idPoint);
            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) )
            {
                const real_t basisValue = basisValues(idActive, idPoint);
                rhs(2*idBasis,   0) += tmp*basisValue*nrml(0);
                rhs(2*idBasis+1, 0) += tmp*basisValue*nrml(1);
            }
        }
    }
    return rhs;
}


/*************** Point distance minimization (PDM) ***************/
gsMatrix<> assembleMatrixPDM(const gsMatrix<>& basisValues,
                             const gsMatrix<index_t>& basisActives,
                             const gsMatrix<>& sampledPtsTargetCurve,
                             const gsMatrix<>& targetPoints,
                             const gsMatrix<>& coefs)
{
    const int M = basisValues.cols();
    const int N = coefs.rows();
    gsMatrix<> A(2*N, 2*N);
    A.setZero();

    gsMatrix<> nrmls = sampledPtsTargetCurve - targetPoints;
    normalize(nrmls);
    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {
        for (index_t idActive1 = 0; idActive1 < basisActives.rows(); idActive1++) // Works as l
        {
            const int idBasis1 = basisActives(idActive1, idPoint);

            for (index_t idActive2 = 0; idActive2 < basisActives.rows(); idActive2++) // Works as i
            {
                const int idBasis2 = basisActives(idActive2, idPoint);

                if (    ((idBasis1 > 0) || ((idBasis1 == 0) && (idActive1 == 0)))
                     && ((idBasis2 > 0) || ((idBasis2 == 0) && (idActive2 == 0))) )
                {
                    //gsInfo << "2*idBasis1: " << 2*idBasis1 << " 2*idBasis2: " << 2*idBasis2 << std::endl;
                    const real_t basisValue1 = basisValues(idActive1, idPoint);
                    const real_t basisValue2 = basisValues(idActive2, idPoint);
                    const real_t nrml0 = nrmls(0, idPoint);
                    const real_t nrml1 = nrmls(1, idPoint);

                    A(2*idBasis1,   2*idBasis2)   += basisValue2 * nrml0 * basisValue1 * nrml0; //
                    A(2*idBasis1,   2*idBasis2+1) += basisValue2 * nrml1 * basisValue1 * nrml0; //

                    A(2*idBasis1+1, 2*idBasis2)   += basisValue2 * nrml0 * basisValue1 * nrml1; //
                    A(2*idBasis1+1, 2*idBasis2+1) += basisValue2 * nrml1 * basisValue1 * nrml1; //
                }
            }
        }
    }

    return A;
}

gsMatrix<> assembleRHSPDM(const gsMatrix<>& basisValues,
                          const gsMatrix<index_t>& basisActives,
                          const gsMatrix<>& sampledPtsTargetCurve,
                          const gsMatrix<>& targetPoints,
                          const gsMatrix<>& coefs)
{
    const int M = basisValues.cols();
    const int N = coefs.rows();
    gsMatrix<> rhs(2*N, 1);
    rhs.setZero();

    gsMatrix<> nrmls = sampledPtsTargetCurve - targetPoints;
    gsMatrix<> vNorms = nrmls.colwise().norm();
    normalize(nrmls);

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {
        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
        {
            int idBasis = basisActives(idActive, idPoint);
            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) )
            {
                const real_t basisValue = basisValues(idActive, idPoint);
                rhs(2*idBasis,   0) += vNorms(0, idPoint) * basisValue * nrmls(0, idPoint);
                rhs(2*idBasis+1, 0) += vNorms(0, idPoint) * basisValue * nrmls(1, idPoint);
            }
        }
    }
    return rhs;
}

/*************** Point distance minimization (PDM) - direct approach (without Gauss--Newton) ***************/
gsMatrix<> assembleMatrixPDMdirect(const gsMatrix<>& basisValues,
                             const gsMatrix<index_t>& basisActives,
                             const gsMatrix<>& coefs)
{
    const int M = basisValues.cols();
    const int N = coefs.rows();
    gsMatrix<> A(2*N, 2*N);
    A.setZero();

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {
        for (index_t idActive1 = 0; idActive1 < basisActives.rows(); idActive1++) // Works as l
        {
            const int idBasis1 = basisActives(idActive1, idPoint);

            for (index_t idActive2 = 0; idActive2 < basisActives.rows(); idActive2++) // Works as i
            {
                const int idBasis2 = basisActives(idActive2, idPoint);

                if (    ((idBasis1 > 0) || ((idBasis1 == 0) && (idActive1 == 0)))
                     && ((idBasis2 > 0) || ((idBasis2 == 0) && (idActive2 == 0))) )
                {
                    const real_t basisValue1 = basisValues(idActive1, idPoint);
                    const real_t basisValue2 = basisValues(idActive2, idPoint);

                    A(2*idBasis1,   2*idBasis2)   += basisValue2 * basisValue1; //
                    A(2*idBasis1+1, 2*idBasis2+1) += basisValue2 * basisValue1; //
                }
            }
        }
    }

    return A;
}

gsMatrix<> assembleRHSPDMdirect(const gsMatrix<>& basisValues,
                                const gsMatrix<index_t>& basisActives,
                                const gsMatrix<>& sampledPtsTargetCurve,
                                const gsMatrix<>& targetPoints,
                                const gsMatrix<>& coefs)
{
    const int M = basisValues.cols();
    const int N = coefs.rows(); 
    gsMatrix<> rhs(2*N, 1);
    rhs.setZero();

    gsMatrix<> nrmls = sampledPtsTargetCurve - targetPoints;

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {
        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
        {
            int idBasis = basisActives(idActive, idPoint);

            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) )
            {
                const real_t basisValue = basisValues(idActive, idPoint);
                rhs(2*idBasis,   0) += basisValue * nrmls(0, idPoint);
                rhs(2*idBasis+1, 0) += basisValue * nrmls(1, idPoint);
            }
        }
    }

    return rhs;
}

/*************** Tikhonov regularization ***************/

void applyTikhonovRegularization(real_t lambda,
                                 gsMatrix<>& A_mat)
{
    for (index_t i = 0; i < A_mat.rows(); i++)
    {
        // Regularization terms are added to the diagonal of the matrix
        A_mat(i, i) += lambda;
    }
}

/*************** Smoothing ***************/

void applySmoothing(real_t lambda,
                    gsMatrix<>& A_mat,
                    gsMatrix<>& rhs,
                    gsGeometry<>& map)
{
    const gsBasis<>& basis = map.basis();
    const int dim = basis.dim();
    const int stride = dim*(dim+1)/2;

    gsVector<int> numNodes(dim);
    for ( int i = 0; i!= dim; ++i )
        numNodes[i] = basis.degree(i)+1;
    gsGaussRule<> QuRule( numNodes ); // Reference Quadrature rule
    gsMatrix<> quNodes, der2, mapDer2, localA, localRhs;
    gsVector<> quWeights;
    gsMatrix<index_t> actives;

    gsBasis<>::domainIter domIt = basis.makeDomainIterator();

    const real_t sqrt2 = math::sqrt((real_t)2);
    for (; domIt->good(); domIt->next() )
    {
        // Map the Quadrature rule to the element and compute basis derivatives
        QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );
        basis.deriv2_into(quNodes, der2);
        basis.active_into(domIt->center, actives);
        const index_t numActive = actives.rows();
        localA.setZero(numActive, numActive );
        localRhs.setZero(2*numActive, 1 );

        map.deriv2_into(quNodes, mapDer2);

        // perform the quadrature
        for (index_t k = 0; k < quWeights.rows(); ++k)
        {
            const real_t weight = quWeights[k] * lambda;

            for (index_t i=0; i!=numActive; ++i)
            {
                // Matrix
                for (index_t j=0; j!=numActive; ++j)
                {
                    real_t localAij = 0; // temporary variable

                    for (int s = 0; s < stride; s++)
                    {
                        // 1 -- pure derivatives, 2 -- mixed derivatives
                        const real_t factor = (s < dim) ? 1.0 : sqrt2;
                        localAij += factor * der2(i * stride + s, k) *
                                             der2(j * stride + s, k);
                    }

                    localA(i, j) += weight * localAij;

                }

                // Right-hand side
                //real_t localRhsi = 0; // temporary variable
                const int ii = actives(i, 0);
                for (int s = 0; s < stride; s++)
                {
                    // 1 -- pure derivatives, 2 -- mixed derivatives
                    const real_t factor = (s < dim) ? 1.0 : sqrt2;
                    rhs(2*ii, 0)   -= weight * factor * der2(i * stride + s, k) *
                                                        mapDer2(s, k);
                    rhs(2*ii+1, 0) -= weight * factor * der2(i * stride + s, k) *
                                                        mapDer2(stride + s, k);
                }

            }
        }

        for (index_t i=0; i!=numActive; ++i)
        {

            const int ii = actives(i,0);
            for (index_t j=0; j!=numActive; ++j)
            {
                const int jj = actives(j,0);
                A_mat( 2*ii,   2*jj )   += localA(i,j);
                A_mat( 2*ii+1, 2*jj+1 ) += localA(i,j);
            }
        }
    }
}

// Solves homogeneous linear system
// (function is used for testing)
gsMatrix<> solveTest(const gsMatrix<>& A)
{
    gsSparseMatrix<> sparseA = A.sparseView();
    sparseA.makeCompressed();

    gsMatrix<> b(A.cols(), 1);
    b.setZero();

    gsSparseSolver<>::LU solver(sparseA);

    gsMatrix<> res = solver.solve(b);

    return res;
}

/*************** Solving control point displacement (incremental updated) ***************/

gsMatrix<> solveDisplacements(const gsMatrix<>& A,
                              const gsMatrix<>& b,
                              const bool debug,
                              const std::string debPrefix)
{

    // (8.?) Configurating the linear solver
    gsSparseMatrix<> sparseA = A.sparseView();
    sparseA.makeCompressed();

    //----------------------EIGEN-ITERATIVE-SOLVERS----------------------//
    //gsSparseSolver<>::CGIdentity solver(sparseA);
    //gsSparseSolver<>::CGDiagonal solver(sparseA);
    //gsSparseSolver<>::BiCGSTABIdentity solver(sparseA);
    //gsSparseSolver<>::BiCGSTABDiagonal solver(sparseA);
    gsSparseSolver<>::BiCGSTABILUT solver(sparseA);

    if (solver.preconditioner().info() != Eigen::Success)
    {
        gsWarn << "The preconditioner failed. Aborting...\n";
        return gsMatrix<> ();
    }

    //----------------------EIGEN-DIRECT-SOLVERS----------------------//
    //gsSparseSolver<>::SimplicialLDLT solver(sparseA);
    //gsSparseSolver<>::QR solver(sparseA);
    //gsSparseSolver<>::LU solver(sparseA);
    //gsSparseSolver<>::fullPivLu solver(sparseA);

    gsMatrix<> tmpDelta = solver.solve(b);

    // (8.?) Solving the linear system, i.e. calculating displacements of the control points
    //gsMatrix<> tmpDelta = solver.solve(b);

    /*
    if (debug)
    {
        // Calculating and printing residuals
        //gsMatrix<> residuals = sparseA*tmpDelta-b;
        gsMatrix<> residuals = A*tmpDelta-b;
        gsInfo << debPrefix << "Max Lp norm residuals: "
               << residuals.lpNorm<Eigen::Infinity>() << std::endl;
    }
    */

    const int N = tmpDelta.rows()/2;
    gsMatrix<> delta (N, 2);
    delta.setZero();

    for (int i = 0; i < N; i++)
    {
        delta(i, 0) = tmpDelta(2*i,   0);
        delta(i, 1) = tmpDelta(2*i+1, 0);
    }

    return delta;
}


// wrong functional, fix it 
real_t computeFunctionalTDM(const gsMatrix<>& sampledPtsTargetCurve,
                            const gsMatrix<>& targetPoints,
                            const gsMatrix<>& nrmls)
{
    real_t F = 0.0;
    const int M = sampledPtsTargetCurve.cols();

    real_t error = 0.0;
    for (index_t idPoint = 0; idPoint < M; idPoint++)
    {

        const gsVector<> P = sampledPtsTargetCurve.col(idPoint);
        const gsVector<> targetPoint = targetPoints.col(idPoint);
        const gsVector<> nrml = nrmls.col(idPoint);
        const real_t tmp = (P - targetPoint).transpose() * nrml;

        F += tmp*tmp;

        real_t norm = (P - targetPoint).norm();
        if (norm > error)
        {
            error = norm;
        }

    }
    //gsInfo << "Max. Euclidian error: " << error << std::endl;
    gsInfo << "Max. Euclidian error: " << std::setprecision(0)
           << std::scientific << error << std::endl;

    return F;    
}

// wrong functional, fix it 
real_t computeFunctionalTDM(const gsMatrix<>& sampledPtsTargetCurve,
                         const gsMatrix<>& templatePoints,
                         const gsGeometry<>& thbMap,
                         const gsMatrix<>& nrmls)
{
    gsMatrix<> targetPoints;
    thbMap.eval_into(templatePoints, targetPoints);

    return computeFunctionalTDM(sampledPtsTargetCurve, targetPoints, nrmls);
}


std::pair<real_t, real_t>
computePDMerror(const gsMatrix<>& sampledPtsTargetCurve,
                const gsMatrix<>& templatePoints,
                const gsGeometry<>& thbMap)
{
    gsMatrix<> targetPoints;
    thbMap.eval_into(templatePoints, targetPoints);

    const int M = sampledPtsTargetCurve.cols();
    
    real_t errorMax = 0;
    real_t errorL2 = 0;

    for (index_t idPoint = 0; idPoint != M; idPoint++)
    {
        const gsVector<> P = sampledPtsTargetCurve.col(idPoint);
        const gsVector<> targetPoint = targetPoints.col(idPoint);
        const real_t err = (P - targetPoint).squaredNorm();
        errorL2 += err;
        
        if (err > errorMax * errorMax)
        {
            errorMax = math::sqrt(err);
        }
    }
    
    errorL2 = math::sqrt(errorL2 / M);
    
    return std::make_pair(errorL2, errorMax);
}

std::pair<real_t, real_t>
computeTDMerror(const gsMatrix<>& sampledPtsTargetCurve,
                const gsMatrix<>& templatePoints,
                const gsMatrix<>& derivTempBoundaryGeom,
                const gsGeometry<>& thbMap)
{
    gsMatrix<index_t> basisActives;
    thbMap.basis().active_into(templatePoints, basisActives);

    gsMatrix<> derivMapBasis;
    thbMap.basis().deriv_into(templatePoints, derivMapBasis);
    
    gsMatrix<> nrmls = normals(derivMapBasis, basisActives, thbMap.coefs(),
                               derivTempBoundaryGeom);
//    gsInfo << "nrmls:\n" << nrmls << std::endl;
    gsMatrix<> targetPoints;
    thbMap.eval_into(templatePoints, targetPoints);
    //gsInfo << "targetPoints (computeTDMerror):\n" << targetPoints << std::endl;

    real_t errorMax = 0;
    real_t errorL2 = 0;

    const int M = sampledPtsTargetCurve.cols();
    for (index_t idPoint = 0; idPoint != M; idPoint++)
    {
        const gsVector<> P = sampledPtsTargetCurve.col(idPoint);
        const gsVector<> targetPoint = targetPoints.col(idPoint);
        const gsVector<> N = nrmls.col(idPoint);

        //gsInfo << "diff: " << P - targetPoint << std::endl;
        const real_t err = math::abs( (P - targetPoint).dot(N) );
        //gsInfo << "err: " << err << std::endl;
        const real_t err2 = err * err;
        errorL2 += err2;

        if (err > errorMax)
        {
            errorMax = err;
        }
    }

    errorL2 = math::sqrt(errorL2*(1.0/M));
    return std::make_pair(errorL2, errorMax);
}

real_t computeMAXerror(const gsMatrix<>& sampledPtsTargetCurve,
                       const gsMatrix<>& templatePoints,
                       const gsGeometry<>& thbMap)
{
    gsMatrix<> targetPoints;
    thbMap.eval_into(templatePoints, targetPoints);

    const int M = sampledPtsTargetCurve.cols();

    real_t errorMax = 0;

    for (index_t idPoint = 0; idPoint != M; idPoint++)
    {
        const gsVector<> P = sampledPtsTargetCurve.col(idPoint);
        const gsVector<> targetPoint = targetPoints.col(idPoint);
        const real_t err = (P - targetPoint).squaredNorm();

        if (err > errorMax * errorMax)
        {
            errorMax = math::sqrt(err);
        }
    }

    return errorMax;
}


gsMatrix<> computeClosestPointsParameters(const gsGeometry<>& thbMap,
                                          const gsMatrix<>& sampledPtsTargetCurve,
                                          const gsMatrix<>& ptsTemplateCurve,
                                          const gsMatrix<>& sampleParameters)
{
    //gsInfo << "sampleParameters:\n" << sampleParameters << std::endl;
    // Use THB-splines map to map template boundary curve points on the target geometry
    gsMatrix<> mappedPtsTargetCurve = thbMap.eval(ptsTemplateCurve);

    // DEBUGGING_STEP
    gsWriteParaviewPoints(mappedPtsTargetCurve, "debug_mappedPtsTargetCurve");

    // Closest points computation
    const int M = sampledPtsTargetCurve.cols();
    // Parameters corresponding to the closest points:
    // [t_1 t_2 ... t_M]
    gsMatrix<> parameters(1, M);
    for (index_t i = 0; i != M; i++)
    {
        index_t ind = findClosestIndex(sampledPtsTargetCurve.col(i), mappedPtsTargetCurve);
        //gsInfo << "i: " << i << " ind: " << ind << std::endl;
        parameters(0, i) = sampleParameters(0, ind);
    }

    //gsInfo << "parameters:\n" << parameters << std::endl;
    return parameters;
}

gsMatrix<> normals(const gsMatrix<>& derivMapBasis,
                   const gsMatrix<index_t>& basisActives,
                   const gsMatrix<>& coefs,
                   const gsMatrix<>& derivTempBoundaryGeom)
{
    // Normals
    // [N_1 N_2 ... N_M] =
    // N(t_1)_0  N(t_2)_0 ... N(t_M)_0
    // N(t_1)_1  N(t_2)_1 ... N(t_M)_1
    const int M = derivMapBasis.cols();
    gsMatrix<> nrmls(2, M);
    nrmls.setZero();
    for (index_t idPoint = 0; idPoint < M; idPoint++)
    {
        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
        {
            int idBasis = basisActives(idActive, idPoint);
            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) ) //
            {
                const real_t tmpDeriv = (derivMapBasis(2*idActive,   idPoint)*derivTempBoundaryGeom(0, idPoint) +
                                         derivMapBasis(2*idActive+1, idPoint)*derivTempBoundaryGeom(1, idPoint));
                nrmls(0, idPoint) += -coefs(idBasis, 1)*tmpDeriv;
                nrmls(1, idPoint) +=  coefs(idBasis, 0)*tmpDeriv;
            }
        }
    }
    normalize(nrmls);
    //gsInfo << "nrmls (0):\n" << nrmls.row(0) << std::endl;
    //gsInfo << "nrmls (1):\n" << nrmls.row(1) << std::endl;
    return nrmls;
}

void saveData(const gsGeometry<>& thbMap,
              const gsMultiPatch<>& temp,
              const std::string outputPrefix)
{
    gsMultiPatch<> newMp = mapMultipatch(thbMap, temp);

    writeMultipatch(newMp, outputPrefix);
    writeKnots(thbMap, outputPrefix);
}

/* UNFINISHED
gsMatrix<> computeJacobian(const gsMatrix<>& basisValues,
                           const gsMatrix<index_t>& basisActives,
                           const gsMatrix<>& nrmls,
                           const gsMatrix<>& coefs,
                           const real_t lambda)
{
    const int M = basisValues.cols();
    const int N = coefs.rows();
    const sqrtLambda = math::sqrt(lambda);

    //gsMatrix<> J(M, 2*N);
    gsMatrix<> J(M+2*N, 2*N);
    J.setZero();

    for (index_t i = 0; i < M; i++) // rows: 1, 2, ..., M
    {
        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
        {
            int idBasis = basisActives(idActive, i);
            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) ) //
            {
                J(i, 2*idBasis)   = basisValues(idActive, i)*nrmls(0, i);
                J(i, 2*idBasis+1) = basisValues(idActive, i)*nrmls(1, i);
            }
        }
    }

    for (index_t i = 0; i < N; i++) // rows: M+1, M+2, ..., M+N
    {
        for (index_t j = 0; j < N; j++)
        {
            //const real_t coefNorm = (coefs.row(j)).norm();
            J(M+i, 2*j)   = sqrtLambda; // d q_l / d Delta d_(l, 1)
            J(M+i, 2*j+1) = 0.0; // d q_l / d Delta d_(l, 2)
        }
    }

    return J;
}
*/

/* UNFINISHED */
//gsMatrix<> computeQ(const gsMatrix<>& basisValues,
//                    const gsMatrix<index_t>& basisActives,
//                    const gsMatrix<>& targetPoints,
//                    const gsMatrix<>& nrmls,
//                    const gsMatrix<>& coefs,
//                    const real_t lambda)
//{

//    //gsInfo << "nrmls " << nrmls << std::endl;


//    const int M = targetPoints.cols();
//    const int N = coefs.rows();

//    gsMatrix<> q(M+2*N, 1);
//    //gsMatrix<> q(M, 1);
//    q.setZero();
//    // [q_1, q_2, ..., q_M]
//    for (index_t idPoint = 0; idPoint < M; idPoint++)
//    {
//        gsVector<> tmp(2);
//        tmp.setZero();
//        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
//        {
//            int idBasis = basisActives(idActive, idPoint);
//            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) ) //
//            {
//                tmp += (coefs.row(idBasis).transpose())*basisValues(idActive, idPoint);
//            }
//        }
//        //gsInfo << "tmp " << tmp << std::endl;
//        /*
//        gsInfo << "tmp  " << tmp.transpose() << "\n"
//               << "targetPoints.col(idPoint) " << targetPoints.col(idPoint).transpose() << "\n"
//               << "vec1 " << (tmp - targetPoints.col(idPoint)).transpose() << "\n"
//               << "vec2 " << nrmls.col(idPoint).transpose() << "\n"
//               << std::endl;
//        */
//        gsMatrix<> tmp2 = (tmp - targetPoints.col(idPoint)).transpose()*nrmls.col(idPoint);
//        //gsInfo << "tmp2 " << tmp2 << std::endl;
//        q(idPoint, 0) = tmp2(0, 0);
//    }
//    // [q_(M+1), q_(M+2), ..., q_(M+N)]
//    /*
//    for (index_t i = 0; i < N; i++)
//    {
//        q(M+i, 0) = 0.0;//l
//    }
//    */
//    //gsInfo << "q: " << q << std::endl;

//    return q;
//}

//void BasisToMatrix(const gsMatrix<>& basis, const gsMatrix<index_t>& actives, gsMatrix<>& result)
//{
//    for (index_t iPoint = 0; iPoint < actives.cols(); iPoint++) // For each point
//    {
//        for (index_t iBasis = 0; iBasis < actives.rows(); iBasis++) // For each basis function
//        {
//            index_t numActive = actives(iBasis, iPoint);
//            if ( (numActive > 0) || ((numActive == 0) && (iBasis == 0)) )
//            {
//                result(numActive, iPoint) = basis(iBasis, iPoint);
//            }
//        }
//    }
//}
