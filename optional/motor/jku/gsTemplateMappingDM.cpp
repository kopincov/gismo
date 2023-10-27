/** @file gsTemplateMappingDM.cpp

    @brief  Constructs a map which transforms (maps) the template shape
    to the target shape. The map is constructed using point distance
    minimization (PDM) or tangent distance minimization (TDM) procedures.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh, S. Sajavicius
*/

// ${GISMO_BIN}/gsTemplateMappingDMvol --pcorrect --adaptiveRef --templateFile wiggly/input/template_multipatch.xml --geometryFile wiggly/input/target_multipatch.xml --mapFile wiggly/input/thb_map_initial.xml --degree 3 --numKnots 11 --numIter 5 --numPtsT 200 --numPts 10000 --weightPDM 1 --weightTDM 0 --lambdaT 1e-2 --lambdaS 0 --numRefIter 0 --extension 1 --kmin -1.68 --kmax 1.68 -o results_dm

#include <utility>

#include <gismo.h>
#include "gsMotorUtils.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

template <unsigned dim>
gsSparseMatrix<> assembleMatrix(const gsMatrix<>& basisValues,
                          const gsMatrix<index_t>& basisActives,
                          const real_t weightPDM,
                          const real_t weightTDM,
                          const gsMatrix<>& nrmls,
                          const gsMatrix<>& coefs);

template <unsigned dim>
gsMatrix<> assembleRHS(const gsMatrix<>& basisValues,
                       const gsMatrix<index_t>& basisActives,
                       const real_t weightPDM,
                       const real_t weightTDM,
                       const gsMatrix<>& sampledPtsTargetCurve,
                       const gsMatrix<>& nrmls,
                       const gsMatrix<>& coefs);

template <unsigned dim>
gsMatrix<> assembleRHS(const gsMatrix<>& basisValues,
                       const gsMatrix<index_t>& basisActives,
                       const real_t weightPDM,
                       const real_t weightTDM,
                       const gsMatrix<>& sampledPtsTargetCurve,
                       const gsMatrix<>& targetPts,
                       const gsMatrix<>& nrmls,
                       const gsMatrix<>& coefs);

template <unsigned dim>
void assembleMatrixPDM(gsSparseMatrix<>& A,
                       const gsMatrix<>& basisValues,
                       const gsMatrix<index_t>& basisActives,
                       const real_t weightPDM);

template <unsigned dim>
void assembleRHSPDM(gsMatrix<>& rhs,
                    const gsMatrix<>& basisValues,
                    const gsMatrix<index_t>& basisActives,
                    const real_t weightPDM,
                    const gsMatrix<>& sampledPtsTargetCurve);

template <unsigned dim>
void assembleRHSPDMSolveUpdates(gsMatrix<>& rhs,
                    const gsMatrix<>& basisValues,
                    const gsMatrix<index_t>& basisActives,
                    const real_t weightPDM,
                    const gsMatrix<>& sampledPtsTargetCurve,
                    const gsMatrix<>& targetPts);

template <unsigned dim>
void assembleMatrixTDM(gsSparseMatrix<>& A,
                       const gsMatrix<>& basisValues,
                       const gsMatrix<index_t>& basisActives,
                       const real_t weightTDM,
                       const gsMatrix<>& nrmls);

template <unsigned dim>
void assembleRHSTDM(gsMatrix<>& rhs,
                    const gsMatrix<>& basisValues,
                    const gsMatrix<index_t>& basisActives,
                    const real_t weightTDM,
                    const gsMatrix<>& sampledPtsTargetCurve,
                    const gsMatrix<>& nrmls);

template <unsigned dim>
void assembleRHSTDMSolveUpdates(gsMatrix<>& rhs,
                    const gsMatrix<>& basisValues,
                    const gsMatrix<index_t>& basisActives,
                    const real_t weightTDM,
                    const gsMatrix<>& sampledPtsTargetCurve,
                    const gsMatrix<>& targetPts,
                    const gsMatrix<>& nrmls);

template <unsigned dim>
void applyTikhonovRegularization(gsSparseMatrix<>& A_mat,
                                 gsMatrix<>& rhs,
                                 const real_t lambda,
                                 const gsMatrix<>& coefs);

void applySmoothing(gsSparseMatrix<>& A_mat,
                    gsMatrix<>& rhs,
                    real_t lambda,
                    gsGeometry<>& map);

template <unsigned dim>
gsMatrix<> solveCoefs(const gsSparseMatrix<>& A,
                      const gsMatrix<>& b,
                      const bool debug = false,
                      const std::string debPrefix = "");

template <unsigned dim>
gsMatrix<> normals(const gsMatrix<>& derivMapBasis,
                   const gsMatrix<index_t>& basisActives,
                   const gsMatrix<> &coefs,
                   const gsMatrix<>& derivTempBoundaryGeom);

void saveData(const gsGeometry<>& thbMap,
              const gsMultiPatch<>& tempBoundaryGeom,
              const std::string outputPrefix);
real_t computeSquaredError(const gsMatrix<>& dataPoints,
                           const gsMatrix<>& params,
                           const gsGeometry<>& thbMap);
/**********************************************************************/
template<unsigned int dim>
class gsDistanceMinimization
{
public:
    gsDistanceMinimization(gsGeometry<>& thbMap,
                           gsTHBSplineBasis<dim>& thbBasis,
                           const gsMultiPatch<>& templateMultipatch,
                           const gsMultiPatch<>& targetMultipatch)
                           : m_thbMap(thbMap),
                             m_thbBasis(thbBasis),
                             m_templateMultipatch(templateMultipatch),
                             m_targetMultipatch(targetMultipatch),
                             m_nPatches( templateMultipatch.nPatches() ),
                             m_sampleTemplateParameters(),
                             m_sampledPtsTargetMultipatch(),
                             m_sampledPtsTargetCurve(),
                             m_sampledPtsTemplateMultipatch(),
                             m_parametersMultipatch(),
                             m_templatePointsMultipatch(),
                             m_templatePoints(),
                             m_numJacPts(),
                             m_numPlotPts(),
                             m_derivTempBoundaryGeomMultipatch(),
                             m_derivTempBoundaryGeom(),
                             m_debug(false),
                             m_indent("    "),
                             m_debStr(m_indent + "gsDistanceMinimization debug: ")
    {
    }

    /// Sample numPointsTemplate points on the template geometry
    /// and numPointsTarget points on the target geometry
    /// (assumption: numPointsTemplate >> numPointsTarget)
    void computeTargetAndTemplatePoints(const int numPointsTemplate, const int numPointsTarget)
    {
        // Generating sample parameters
        m_sampleTemplateParameters = uniformParameters<dim>(0.0, 1.0, numPointsTemplate);

        const gsMatrix<> sampleTargetParameters = uniformParameters<dim>(0.0, 1.0, numPointsTarget);

        for (std::size_t i = 0; i != m_nPatches ; i++)
        {
            // Sampling points on the template boundary curve
            const gsGeometry<>& templateBoundaryCurve = m_templateMultipatch.patch(i);
            gsMatrix<> ptsTemplateCurve = templateBoundaryCurve.eval(m_sampleTemplateParameters);
            m_sampledPtsTemplateMultipatch.push_back(ptsTemplateCurve);

            // Sampling points on the target boundary curve
            const gsGeometry<>& targetBoundaryCurve = m_targetMultipatch.patch(i);
            gsMatrix<> ptsTargetCurve = targetBoundaryCurve.eval(sampleTargetParameters);
            m_sampledPtsTargetMultipatch.push_back(ptsTargetCurve);
        }
        m_sampledPtsTargetCurve = stdVectorToMatrix(m_sampledPtsTargetMultipatch);

        //gsWriteParaviewPoints(m_sampledPtsTargetCurve, "debug_3d_m_sampledPtsTargetCurve");
        //gsWriteParaviewPoints(stdVectorToMatrix(m_sampledPtsTargetMultipatch), "debug_m_sampledPtsTargetMultipatch");
        //gsWriteParaviewPoints(stdVectorToMatrix(m_sampledPtsTemplateMultipatch), "debug_m_sampledPtsTemplateMultipatch");
    }

    /// Set debug mode
    void set_debug(const bool debug)
    {
        m_debug = debug;
    }

    /// Set map
    void set_map(const gsGeometry<>& thbMap)
    {
        m_thbMap = thbMap;
    }

    /// Set map
    void set_numJacPts(const unsigned numJacPts)
    {
        m_numJacPts = numJacPts;
    }

    /// ...
    void set_numPlotPts(const unsigned numPlotPts)
    {
        m_numPlotPts = numPlotPts;
    }

    /// Compute parameters
    void computeParameters(const bool pcorrect)
    {
        for (std::size_t iPatch= 0; iPatch != m_nPatches; iPatch++)
        {
            // If parameter correction is used, the initial closest points are computed
            if (pcorrect)
            {
                computeInitialClosestPointsParameters();
            }
            else
            {
                computeUniformPointsParameters();
            }
        }

        // If parameter correction is used, the initial closest points are improved
        if (pcorrect)
        {
            for (std::size_t iPatch = 0; iPatch != m_nPatches; iPatch++)
            {
                computeClosestPointsParameters();
            }
        }

    }

    /// Compute initial closest point parameters
    void computeInitialClosestPointsParameters()
    {
        m_parametersMultipatch.clear();
        for (std::size_t iPatch = 0; iPatch != m_nPatches; iPatch++)
        {
            const gsMatrix<> sampledPtsTargetCurve   = m_sampledPtsTargetMultipatch.at(iPatch);
            const gsMatrix<> ptsTemplateCurve = m_sampledPtsTemplateMultipatch.at(iPatch);
            // Use THB-splines map to map template boundary curve points on the target geometry
            gsMatrix<> mappedPtsTargetCurve = m_thbMap.eval(ptsTemplateCurve);
            // Closest points computation
            const int M = sampledPtsTargetCurve.cols();
            // Parameters corresponding to the closest points:
            // [t_1 t_2 ... t_M]
            gsMatrix<> parameters(dim-1, M);
            parameters.setZero();
            for (index_t i = 0; i != M; i++)
            {
                index_t ind = findClosestIndex(sampledPtsTargetCurve.col(i), mappedPtsTargetCurve);
                for (unsigned d = 0; d < dim-1; d++)
                {
                    parameters(d, i) = m_sampleTemplateParameters(d, ind);
                }
            }
            m_parametersMultipatch.push_back(parameters);
        }
    }

    /// Compute uniform (not optimized) point parameters
    void computeUniformPointsParameters()
    {
        m_parametersMultipatch.clear();
        for (std::size_t iPatch = 0; iPatch != m_nPatches; iPatch++)
        {
            const int numPoints = m_sampledPtsTargetMultipatch.at(iPatch).cols();
            const gsMatrix<> parameters = uniformParameters<dim>(0.0, 1.0, numPoints);
            m_parametersMultipatch.push_back(parameters);
        }
    }

    /// Compute optimal parameters (improves initial parameter values by Newton iteration)
    void computeClosestPointsParameters()
    {
        const unsigned numNewtonIt = 2;

        // Executing Newton iterations for optimal parameter computation
        for (unsigned iter=0; iter<numNewtonIt; iter++)
        {
            // Computes m_templatePoints and m_templatePointsMultipatch
            // corresponding to current optimal parameters (m_parametersMultipatch)
            computeTemplatePointsMultipatch();
            // Computes m_derivTempBoundaryGeom corresponding to current optimal parameters (m_parametersMultipatch)
            computeDerivTempBoundaryGeom();

            // Computes basisValues
            gsMatrix<> basisValues;
            basisValues.setZero();
            m_thbMap.basis().eval_into(m_templatePoints, basisValues);
            // Computes basisActives
            gsMatrix<index_t> basisActives;
            basisActives.setZero();
            m_thbMap.basis().active_into(m_templatePoints, basisActives);
            // Computes derMapBasis
            gsMatrix<> derivMapBasis;
            derivMapBasis.setZero();
            m_thbMap.basis().deriv_into(m_templatePoints, derivMapBasis);
            // Map coefficients
            gsMatrix<> coefs = m_thbMap.coefs();

            std::vector< gsMatrix<> > tmp_m_parametersMultipatch = m_parametersMultipatch;
            m_parametersMultipatch.clear();
            for (std::size_t iPatch = 0; iPatch != m_nPatches; iPatch++)
            {
                // Current optimal parameters corresponding to ith patch
                gsMatrix<> parameters = tmp_m_parametersMultipatch.at(iPatch);
                // Target points on the ith patch
                const gsMatrix<> sampledPtsTargetCurve = m_sampledPtsTargetMultipatch.at(iPatch);
                // Template points corresponding to the optimal parameters
                const gsMatrix<> sampledPtsTemplateCurve = m_templatePointsMultipatch.at(iPatch);

                const gsMatrix<> derivTempBoundary = m_derivTempBoundaryGeomMultipatch.at(iPatch);

                index_t M = sampledPtsTargetCurve.cols();
                unsigned numImprPts = 0;
                // Compute map derivatives
                gsMatrix<> derMap(dim*(dim-1), M);
                derMap.setZero();
                for (index_t idPoint = 0; idPoint < M; idPoint++)
                {
                    for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
                    {
                        int idBasis = basisActives(idActive, idPoint);
                        if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) ) //
                        {
                            real_t tmpDeriv = 0.0;

                            real_t tmpDeriv_u = 0.0;
                            real_t tmpDeriv_v = 0.0;

                            if (dim == 2)
                            {
                                tmpDeriv = (derivMapBasis(dim*idActive,   idPoint)*derivTempBoundary(0, idPoint) +
                                            derivMapBasis(dim*idActive+1, idPoint)*derivTempBoundary(1, idPoint));

                                for (unsigned d = 0; d < dim; d++)
                                {
                                    derMap(d, idPoint) += coefs(idBasis, d)*tmpDeriv;
                                    // This "if" is very questionable
                                    if ( math::isnan(derMap(d, idPoint)) )
                                    {
                                        derMap(d, idPoint) = 0.0;
                                    }
                                }

                            }
                            else if (dim ==3)
                            {
                                tmpDeriv_u = (derivMapBasis(dim*idActive,   idPoint)*derivTempBoundary(0, idPoint) +
                                              derivMapBasis(dim*idActive+1, idPoint)*derivTempBoundary(2, idPoint) +
                                              derivMapBasis(dim*idActive+2, idPoint)*derivTempBoundary(4, idPoint));
                                tmpDeriv_v = (derivMapBasis(dim*idActive,   idPoint)*derivTempBoundary(1, idPoint) +
                                              derivMapBasis(dim*idActive+1, idPoint)*derivTempBoundary(3, idPoint) +
                                              derivMapBasis(dim*idActive+2, idPoint)*derivTempBoundary(5, idPoint));
                                derMap(0, idPoint) += coefs(idBasis, 0)*tmpDeriv_u;
                                derMap(1, idPoint) += coefs(idBasis, 1)*tmpDeriv_u;
                                derMap(2, idPoint) += coefs(idBasis, 2)*tmpDeriv_u;
                                derMap(3, idPoint) += coefs(idBasis, 0)*tmpDeriv_v;
                                derMap(4, idPoint) += coefs(idBasis, 1)*tmpDeriv_v;
                                derMap(5, idPoint) += coefs(idBasis, 2)*tmpDeriv_v;
                            }
                        }
                    }

                    // Calculating update of the parameter
                    if (dim == 2)
                    {
                        const gsVector<> point = sampledPtsTargetCurve.col(idPoint);
                        const gsVector<> mapValue = m_thbMap.eval(sampledPtsTemplateCurve.col(idPoint));
                        const gsVector<> derMapValue = derMap.col(idPoint);
                        const real_t derMapNorm = derMapValue.norm();

                        // Store old version of parameter
                        const real_t parameter_old = parameters(0, idPoint);
                        // Compute an incremental update for the parameter
                        const real_t deltaPar = ((point-mapValue).dot(derMapValue))/(derMapNorm*derMapNorm);

                        parameters(0, idPoint) += deltaPar;
                        // Ensuring that all parameter values stay in [0, 1]
                        if (parameters(0, idPoint) < 0)
                        {
                            parameters(0, idPoint) = 0.0;
                        }
                        else if (parameters(0, idPoint) > 1)
                        {
                            parameters(0, idPoint) = 1.0;
                        }

                        // Checking, if the Newton iteration improves the accuracy

                        // The map value wrt. to the updated parameter
                        const gsVector<> mapValueNew = m_thbMap.eval((m_templateMultipatch.patch(iPatch)).eval(parameters.block(0, idPoint, 1, 1)));
                        const real_t errorNorm    = (point-mapValue).norm();
                        const real_t errorNormNew = (point-mapValueNew).norm();

                        if ( errorNorm < errorNormNew )
                        // If update does not improves the accuracy
                        {
                            parameters(0, idPoint) = parameter_old;
                        }
                        else
                        {
                            numImprPts++;
                        }
                    }
                    else if (dim == 3)
                    {
                        const gsVector<> point = sampledPtsTargetCurve.col(idPoint);
                        const gsVector<> mapValue = m_thbMap.eval(sampledPtsTemplateCurve.col(idPoint));
                        const gsVector<> derMapValue_u = derMap.block(0, idPoint, 3, 1);
                        const gsVector<> derMapValue_v = derMap.block(3, idPoint, 3, 1);
                        const real_t derMapNorm_u = derMapValue_u.norm();
                        const real_t derMapNorm_v = derMapValue_v.norm();

                        // Store old version of parameter
                        const gsMatrix<> parameter_old = parameters.block(0, idPoint, dim-1, 1);

                        // Solving linear system for the optimal parameter updates
                        gsMatrix<> mat(2, 2);
                        mat << derMapNorm_u*derMapNorm_u, derMapValue_u.dot(derMapValue_v),
                               derMapValue_u.dot(derMapValue_v), derMapNorm_v*derMapNorm_v;
                        gsSparseMatrix<> sparseMat = mat.sparseView();
                        sparseMat.makeCompressed();

                        gsMatrix<> rhs(2, 1);
                        rhs << (point-mapValue).dot(derMapValue_u),
                               (point-mapValue).dot(derMapValue_v);
                        gsSparseSolver<>::LU solver(sparseMat);
                        gsMatrix<> deltaPar = solver.solve(rhs);

                        parameters(0, idPoint) += deltaPar(0, 0);
                        parameters(1, idPoint) += deltaPar(1, 0);
                        // Ensuring that all parameter values stay in [0, 1]x[0, 1]
                        for (unsigned d = 0; d < dim-1; d++)
                        {
                            if (parameters(d, idPoint) < 0)
                            {
                                parameters(d, idPoint) = 0.0;
                            }
                            else if (parameters(d, idPoint) > 1)
                            {
                                parameters(d, idPoint) = 1.0;
                            }
                        }

                        // Checking, if the Newton iteration improves the accuracy

                        // The map value wrt. to the updated parameter
                        const gsVector<> mapValueNew = m_thbMap.eval((m_templateMultipatch.patch(iPatch)).eval(parameters.block(0, idPoint, dim-1, 1)));
                        const real_t errorNorm    = (point-mapValue).norm();
                        const real_t errorNormNew = (point-mapValueNew).norm();

                        if ( errorNorm < errorNormNew )
                        // If update does not improves the accuracy
                        {
                            parameters(0, idPoint) = parameter_old(0, 0);
                            parameters(1, idPoint) = parameter_old(1, 0);
                        }
                        else
                        {
                            numImprPts++;
                        }

                     }
                }

                m_parametersMultipatch.push_back(parameters);
            }
        }
    }

    /// Compute initial map
    void computeInitialMap(const gsTHBSplineBasis<dim>& thbBasis,
                           const real_t lambdaTikhonov,
                           const real_t lambdaSmoothing,
                           const int numIter,
                           const bool pcorrect,
                           const std::string& outputPrefix)
    {
        gsInfo << m_indent << " Iteration "
               << 0 << " / " << numIter << std::endl;

        // Fill m_parametersMultipatch with parameters uniformly
        // distributed on the unit interval
        this->computeParameters(false);

        // Re-calculate m_templatePoint:
        // the template points corresponding to the parameters in
        computeTemplatePointsMultipatch();

        // Calculate basisValues and basisActives for points from m_templatePoints
        gsMatrix<> basisValues;
        basisValues.setZero();
        thbBasis.eval_into(m_templatePoints, basisValues);
        gsMatrix<index_t> basisActives;
        basisActives.setZero();
        thbBasis.active_into(m_templatePoints, basisActives);

        // Initializing fitting matrix and rhs
        const int N = thbBasis.size();
        //gsMatrix<> A(dim*N, dim*N);
        gsSparseMatrix<> A(dim*N, dim*N);
        A.setZero();
        gsMatrix<> b(dim*N, 1);
        b.setZero();
        gsMatrix<> coefs(N, dim);
        coefs.setZero();

        // Assembling fitting matrix and rhs
        assembleMatrixPDM<dim>(A, basisValues, basisActives, 1.0);
        assembleRHSPDM<dim>(b, basisValues, basisActives, 1.0, m_sampledPtsTargetCurve); // m_sampledPtsTargetCurve = data points
        if (lambdaTikhonov != 0)
        {
            applyTikhonovRegularization<dim>(A, b, lambdaTikhonov, coefs); // coefs consists of zeros for this moment
        }
        if (lambdaSmoothing > 0)
        {
            applySmoothing(A, b, lambdaSmoothing, m_thbMap);
        }

        // Solving coefficients
        coefs = solveCoefs<dim>(A, b);

        if (m_debug)
        {
            gsInfo << m_debStr << "coefs max norm: "
                   << coefs.lpNorm<Eigen::Infinity>() << std::endl;
        }

        // Setting up initial map
        m_thbMap.setCoefs(coefs);

        // If method uses parameter correction, we re-calculate (optimize) parameters here
        if (pcorrect)
        {
            this->computeParameters(pcorrect);
        }

        // Updating m_templatePoints
        computeTemplatePointsMultipatch();

        real_t error = 0;
        // If method uses foot point computation step:
        if (pcorrect)
        {
            error = computeSquaredError(m_sampledPtsTargetCurve, m_templatePoints, m_thbMap);
        }
        // If method does not use foot point computation step:
        else
        {
            // Save m_templatePoints in temp. variable
            gsMatrix<> tmp_m_templatePoints = m_templatePoints;
            std::vector< gsMatrix<> > tmp_m_parametersMultipatch = m_parametersMultipatch;
            // We optimize parameters, but use them for measuring error only
            this->computeParameters(true);
            computeTemplatePointsMultipatch();
            error = computeSquaredError(m_sampledPtsTargetCurve, m_templatePoints, m_thbMap);
            m_parametersMultipatch = tmp_m_parametersMultipatch;
            m_templatePoints = tmp_m_templatePoints;
        }
        gsInfo << m_indent << " Sq. error:    " << error << std::endl;
        if (m_numJacPts > 0)
        {
            checkJacobianDeterminant(m_thbMap, m_numJacPts);
        }

        saveIterationOutput(outputPrefix, 0, m_numPlotPts);
        writeMap(m_thbMap, 0, outputPrefix + "_initial");
        saveData(m_thbMap, m_templateMultipatch, outputPrefix + "_initial");

        std::string targPtsName = makeDebugName("test_", "_m_templatePoints_", 0);
        std::vector< gsMatrix <> > targetPointsMultipatch = getTargetPointsMultipatch();

    }

    /// Upload initial map from the file
    //template <unsigned dim>
    void uploadInitialMap(const std::string& mapFile,
                          const int numIter,
                          const bool pcorrect,
                          const std::string& outputPrefix)
    {
        gsInfo << m_indent << " Iteration "
               << 0 << " / " << numIter << std::endl;

        // Fill m_parametersMultipatch with parameters uniformly
        // distributed on the unit interval
        this->computeParameters(false);

        // Re-calculate m_templatePoint:
        // the template points corresponding to the parameters in  m_parametersMultipatch
        computeTemplatePointsMultipatch();

        // Setting up initial map
        gsGeometry<>* mapPtr = readGeometry(mapFile);
        gsTHBSpline<dim>& thbMap = static_cast<gsTHBSpline<dim>&>(*mapPtr);
        set_map(thbMap);

        // If method uses parameter correction, we re-calculate (optimize) parameters here
        if (pcorrect)
        {
            this->computeParameters(pcorrect);
        }

        // Updating m_templatePoints
        computeTemplatePointsMultipatch();

        real_t error = 0;
        // If method uses foot point computation step:
        if (pcorrect)
        {
            error = computeSquaredError(m_sampledPtsTargetCurve, m_templatePoints, m_thbMap);
        }
        // If method does not use foot point computation step:
        else
        {
            // Save m_templatePoints in temp. variable
            gsMatrix<> tmp_m_templatePoints = m_templatePoints;
            std::vector< gsMatrix<> > tmp_m_parametersMultipatch = m_parametersMultipatch;
            // We optimize parameters, but use them for measuring error only
            this->computeParameters(true);
            computeTemplatePointsMultipatch();
            error = computeSquaredError(m_sampledPtsTargetCurve, m_templatePoints, m_thbMap);
            // We recover old (not optimized) parameters
            m_parametersMultipatch = tmp_m_parametersMultipatch;
            m_templatePoints = tmp_m_templatePoints;
        }
        gsInfo << m_indent << " Sq. error:    " << error << std::endl;
        if (m_numJacPts > 0)
        {
            checkJacobianDeterminant(m_thbMap, m_numJacPts);
        }

        saveIterationOutput(outputPrefix, 0, m_numPlotPts);
        writeMap(m_thbMap, 0, outputPrefix + "_initial");
        saveData(m_thbMap, m_templateMultipatch, outputPrefix + "_initial");

    }

    /// Executes distance minimization (PDM or TDM procedure)
    //template <unsigned dim>
    void compute(const int numIter,
                 const real_t weightPDM,
                 const real_t weightTDM,
                 const real_t lambdaTikhonov,
                 const real_t lambdaSmoothing,
                 const bool pcorrect,
                 const bool solveUpdates,
                 const std::string& outputPrefix)
    {
        // Distance minimization procedure iterations
        std::vector< gsMatrix <> > targetPointsMultipatch = getTargetPointsMultipatch();
        gsMatrix<> targetPoints = stdVectorToMatrix(targetPointsMultipatch);

        //gsWriteParaviewPoints(targetPoints, "debug_3d_targetPoints");
        //gsWriteParaviewPoints(m_templatePoints, "debug_3d_m_templatePoints");

        for (int iter = 1; iter <= numIter; iter++)
        {
            gsInfo << m_indent << " Iteration "
                   << iter << " / " << numIter << std::endl;

           // THB-spline map values on the template points
            gsMatrix<> basisValues;
            basisValues.setZero();
            m_thbMap.basis().eval_into(m_templatePoints, basisValues);
            gsMatrix<index_t> basisActives;
            basisActives.setZero();
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
            derivMapBasis.setZero();
            m_thbMap.basis().deriv_into(m_templatePoints, derivMapBasis);

            gsMatrix<>& coefs = m_thbMap.coefs();
            const int M = derivMapBasis.cols();
            gsMatrix<> nrmls(dim, M);
            nrmls.setZero();
            if (weightTDM > 0)
            {
                nrmls = normals<dim>(derivMapBasis, basisActives, coefs, m_derivTempBoundaryGeom);
            }
            int N = coefs.rows();
            gsSparseMatrix<> A(dim*N, dim*N);
            A.setZero();
            gsMatrix<> b(dim*N, 1);
            b.setZero();
            A = assembleMatrix<dim>(basisValues, basisActives, weightPDM, weightTDM, nrmls, coefs);
            if (solveUpdates)
            {
                b = assembleRHS<dim>(basisValues, basisActives, weightPDM, weightTDM, m_sampledPtsTargetCurve, targetPoints, nrmls, coefs);
                if (lambdaTikhonov != 0)
                {
                    applyTikhonovRegularization<dim>(A, b, lambdaTikhonov, coefs);
                }
                if (lambdaSmoothing > 0)
                {
                    applySmoothing(A, b, lambdaSmoothing, m_thbMap);
                }
                gsMatrix<> updCoefs = solveCoefs<dim>(A, b);
                coefs = coefs + updCoefs;
            }
            else
            {
                b = assembleRHS<dim>(basisValues, basisActives, weightPDM, weightTDM, m_sampledPtsTargetCurve, nrmls, coefs);
                if (lambdaTikhonov != 0)
                {
                    applyTikhonovRegularization<dim>(A, b, lambdaTikhonov, coefs);
                }
                if (lambdaSmoothing > 0)
                {
                    applySmoothing(A, b, lambdaSmoothing, m_thbMap);
                }
                coefs = solveCoefs<dim>(A, b);
            }
            if (m_debug)
            {
                gsInfo << m_debStr << "coefs max norm: "
                       << coefs.lpNorm<Eigen::Infinity>() << std::endl;
            }
            m_thbMap.setCoefs(coefs);

            // If method uses parameter correction, we recalculate parameters here
            if (pcorrect)
            {
                this->computeParameters(pcorrect);
            }
            computeTemplatePointsMultipatch();

            real_t error = 0;
            // If method uses foot point computation step:
            if (pcorrect)
            {
                error = computeSquaredError(m_sampledPtsTargetCurve, m_templatePoints, m_thbMap);
            }
            // If method does not use foot point computation step:
            else
            {
                // Save m_templatePoints in temp. variable
                gsMatrix<> tmp_m_templatePoints = m_templatePoints;
                std::vector< gsMatrix<> > tmp_m_parametersMultipatch = m_parametersMultipatch;
                // We optimize parameters, but use them for measuring error only
                this->computeParameters(true);
                computeTemplatePointsMultipatch();
                error = computeSquaredError(m_sampledPtsTargetCurve, m_templatePoints, m_thbMap);
                m_parametersMultipatch = tmp_m_parametersMultipatch;
                m_templatePoints = tmp_m_templatePoints;
            }
            gsInfo << m_indent << " Sq. error:    " << error << std::endl;
            if (m_numJacPts > 0)
            {
                checkJacobianDeterminant(m_thbMap, m_numJacPts);
            }
            saveIterationOutput(outputPrefix, iter, m_numPlotPts);

            std::string targPtsName = makeDebugName("test_", "_m_templatePoints_", iter+1);
            std::vector< gsMatrix <> > targetPointsMultipatch = getTargetPointsMultipatch();
            gsMatrix<> targetPoints = stdVectorToMatrix(targetPointsMultipatch);
            gsWriteParaviewPoints(targetPoints, targPtsName);

            // Giving different names for files with (possibly) final results
            writeMap(m_thbMap, 0,"results_dm_final");
        }
    }

    /// Computes derivatives of the template boundary curve
    void computeDerivTempBoundaryGeom()
    {
        // Derivatives of template boundary curve (geometry) at parameter values
        // gamma_1'(t_1) gamma_1'(t_2) ... gamma_1'(t_M)
        // gamma_2'(t_1) gamma_2'(t_2) ... gamma_2'(t_M)
        m_derivTempBoundaryGeomMultipatch.clear();
        for (size_t i = 0; i != m_nPatches; i++)
        {
            gsMatrix<> derivTempBoundaryGeomTmp;
            derivTempBoundaryGeomTmp.setZero();
            m_templateMultipatch.patch(i).deriv_into(m_parametersMultipatch.at(i), derivTempBoundaryGeomTmp);
            m_derivTempBoundaryGeomMultipatch.push_back(derivTempBoundaryGeomTmp);
        }
        m_derivTempBoundaryGeom = stdVectorToMatrix(m_derivTempBoundaryGeomMultipatch);
    }

    /// Returns m_templatePoints
    const gsMatrix<>& getTemplatePoints() const
    {
        return m_templatePoints;
    }

    /// Returns m_sampledPtsTargetCurve
    const gsMatrix<>& getSampledTargetPoints() const
    {
        return m_sampledPtsTargetCurve;
    }

private:

    /// Saves output of the iteration
    void saveIterationOutput(const std::string& outputPrefix,
                             const int iter,
                             const int numPts)
    {
        const std::string prefix = outputPrefix + "Iter" + util::to_string(iter);
        gsMultiPatch<> mappedMp = mapMultipatch(m_thbMap, m_templateMultipatch);

        writeMultipatch(mappedMp, prefix + "_points", numPts);
        writeMap(m_thbMap, prefix);
        writeKnots(m_thbMap, prefix);
    }

    /// Computes ...
    void computeTemplatePointsMultipatch()
    {
        m_templatePointsMultipatch.clear();
        for (std::size_t i = 0; i != m_nPatches; i++)
        {
            // Calculate points on the template boundary curve
            gsMatrix<> templatePoints;
            templatePoints.setZero();
            m_templateMultipatch.patch(i).eval_into(m_parametersMultipatch.at(i), templatePoints);
            m_templatePointsMultipatch.push_back(templatePoints);
        }
        m_templatePoints = stdVectorToMatrix(m_templatePointsMultipatch);
    }

    /// ...
    std::vector< gsMatrix <> > getTargetPointsMultipatch()
    {
        std::vector< gsMatrix <> > targetPointsMultipatch;
        for (std::size_t i = 0; i != m_nPatches; i++)
        {
            // Calculate points on the template boundary curve
            gsMatrix<> templatePoints;
            templatePoints.setZero();
            m_templateMultipatch.patch(i).eval_into(m_parametersMultipatch.at(i), templatePoints);
            // Use THB-splines map to map these points (see previous step) on the target geometry
            gsMatrix<> targetPointsTmp;
            targetPointsTmp.setZero();
            m_thbMap.eval_into(templatePoints, targetPointsTmp);
            targetPointsMultipatch.push_back(targetPointsTmp);
        }
        return targetPointsMultipatch;
    }

    /// ...
    std::string makeDebugName(const std::string& outputPrefix,
                              const std::string& name)
    {
        std::string newName = "debug_" + outputPrefix + name;
        return newName;
    }

    /// ...
    std::string makeDebugName(const std::string& outputPrefix,
                              const std::string& name,
                              const int iteration)
    {
        return makeDebugName(outputPrefix + "Iter" + util::to_string(iteration), name);
    }

    // THB-Spline map
    gsGeometry<>& m_thbMap;

    // THB-Spline basis used to define the THB-Spline map
    gsTHBSplineBasis<dim>& m_thbBasis;

    // Multipatch template geometry
    const gsMultiPatch<>& m_templateMultipatch;

    // Multipatch target geometry
    const gsMultiPatch<>& m_targetMultipatch;

    // Number of patches (in template or target geometry)
    std::size_t m_nPatches;

    // Parameters corresponding to the points sampled on each of the template boundaries
    // (number of points = numPointsTemplate)
    gsMatrix<> m_sampleTemplateParameters;

    // Points on the target boundaries stored as std::vector of several gsMatrix structures
    std::vector< gsMatrix<> > m_sampledPtsTargetMultipatch;

    // Points on the target boundaries stored in gsMatrix format
    // (the same as m_sampledPtsTargetMultipatch, just stored in one gsMatrix)
    gsMatrix<> m_sampledPtsTargetCurve; // P_k

    // Points on the template boundaries stored as std::vector of several gsMatrix structures
    std::vector< gsMatrix<> > m_sampledPtsTemplateMultipatch;

    // Parameters
    std::vector< gsMatrix<> > m_parametersMultipatch;

    // Points on the template boundaries corresponding to parameters from m_parametersMultipatch
    // and stored as std::vector of several gsMatrix structures
    std::vector< gsMatrix<> > m_templatePointsMultipatch;

    // Points on the template boundaries corresponding to parameters from m_parametersMultipatch
    // and in gsMatrix format
    // (the same as m_templatePointsMultipatch, just stored in one gsMatrix)
    gsMatrix<> m_templatePoints;

    // Number of points to check Jacobian
    unsigned m_numJacPts;

    // Number of points for Praview
    unsigned m_numPlotPts;

    //
    std::vector< gsMatrix<> > m_derivTempBoundaryGeomMultipatch;

    //
    gsMatrix<> m_derivTempBoundaryGeom;

    // Debug mode
    bool m_debug;

    // Auxiliary variable
    const std::string m_indent;

    // Auxiliary variable
    std::string m_debStr;

};

/// Computes errors
gsVector<> computeErrors(const gsGeometry<>& thbMap,
                         const gsMatrix<>& params,
                         const gsMatrix<>& points)
{
    gsVector<> pointErrors(params.cols());

    gsMatrix<> values;
    values.setZero();
    thbMap.eval_into(params, values);

    for (index_t i = 0; i < points.cols(); i++)
    {
        const real_t err = (points.col(i) - values.col(i)).norm();
        //const real_t err = (points.col(i) - values.col(i)).squaredNorm();
        pointErrors(i) = err;
    }

    return pointErrors;
}

/// ...
template <unsigned d>
bool isCellAlreadyInserted(const gsVector<index_t, d>& a_cell,
                           const std::vector<index_t>& cells)
{
    for (std::size_t i = 0; i != cells.size(); i+= a_cell.rows())
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

void append(std::vector<index_t>& boxes,
            const gsVector<index_t>& box)
{
    for (index_t row = 0; row != box.rows(); row++)
    {
        boxes.push_back(box[row]);
    }
}

template<unsigned d>
void appendBox(std::vector<index_t>& boxes,
               std::vector<index_t>& cells,
               gsTHBSplineBasis< d >& basis,
               const gsVector<>& parameter,
               const std::vector<index_t>& ext)
{
    const int maxLvl = basis.maxLevel();
    const gsTensorBSplineBasis< d >& tensorBasis =
                *(basis.getBases()[maxLvl]);

    gsVector<index_t, d> a_cell;
    for (unsigned dim = 0; dim != d; dim++)
    {
        const gsKnotVector<>& kv = tensorBasis.component(dim).knots();
        a_cell(dim) = static_cast< unsigned >(kv.uFind(parameter(dim)).uIndex());
    }

    if (!isCellAlreadyInserted<d>(a_cell, cells))
    {
        append(cells, a_cell);

        // get level of a cell
        gsVector<index_t, d> a_cell_up = a_cell + gsVector<index_t, d>::Ones();
        const int cell_lvl = basis.tree().query3(a_cell, a_cell_up, maxLvl) + 1;

        gsVector<index_t> box(2 * d + 1);
        box(0) = cell_lvl;

        for (unsigned dim = 0; dim != d; dim++)
        {
            const unsigned numBreaks = basis.numBreaks(cell_lvl, dim) - 1;

            unsigned lowIndex = 0;
            if (cell_lvl < maxLvl)
            {
                const unsigned shift = maxLvl - cell_lvl;
                lowIndex = (a_cell(dim) >> shift);
            }
            else
            {
                const unsigned shift = cell_lvl - maxLvl;
                lowIndex = (a_cell(dim) << shift);
            }

            const unsigned uext = static_cast<unsigned>(ext[dim]);
            const unsigned low = (lowIndex > uext) ? (lowIndex - uext): 0;
            const unsigned upp = (lowIndex + uext + 1 < numBreaks) ?
                                 (lowIndex + uext + 1): numBreaks;

            box[1 + dim] = low;
            box[1 + d + dim] = upp;
        }
        append(boxes, box);
    }
}

/// ...
template<unsigned d>
std::vector<index_t> getBoxes(gsTHBSplineBasis< d >& basis,
                               const gsMatrix<>& params,
                               const gsVector<>& errors,
                               const real_t threshold,
                               const std::vector<index_t>& ext)
{
    std::vector<index_t> cells;
    std::vector<index_t> boxes;

    for (index_t index = 0; index != errors.rows(); index++)
    {
        if (threshold <= errors(index))
        {
            appendBox<d>(boxes, cells, basis, params.col(index), ext);
        }
    }

    return boxes;
}

/// ...
template<unsigned d>
std::vector<index_t> getBoxesRelative(gsTHBSplineBasis< d >& basis,
                                       const gsMatrix<>& params,
                                       const gsVector<>& errors,
                                       const real_t percentage,
                                       const std::vector<index_t>& ext)
{
    std::vector<index_t> cells;
    std::vector<index_t> boxes;

    // Sorting errors
    gsVector<> errorsTmp = errors;

    int size = errorsTmp.size();
    gsMatrix<> ind_max(1, size);
    ind_max.setZero();
    real_t max_error = 0;
    std::size_t i_max = 0;
    for (index_t i = 0; i != size; i++)
    {
        // Finding index of max element in the rest of the array
        for (index_t ii = i; ii != size; ii++)
        {
            if (errorsTmp[ii] > max_error)
            {
                 i_max = ii;
            }
        }

        // Changing the values
        max_error = errorsTmp[i_max];
        errorsTmp[i_max] = errorsTmp[i];
        errorsTmp[i] = max_error;
        ind_max(0, i) = i_max;
    }

    index_t perc_num = percentage*size + 1;
    gsInfo << "perc_num: " << perc_num << std::endl;

    for (index_t index = 0; index != perc_num; index++)
    {
        appendBox<d>(boxes, cells, basis, params.col(ind_max(0, index)), ext);
    }

    return boxes;
}

/// ...
template<unsigned d>
bool uniformRefine(gsTHBSpline<d>& thbMap)
{
    thbMap.basis().uniformRefine_withCoefs(thbMap.coefs());

    return true;
}

/// ...
template<unsigned d>
bool refine(gsTHBSpline< d >& thbMap,
            const gsMatrix<>& params,
            const gsMatrix<>& points,
            const int extension,
            const real_t threshold,
            const real_t percentage)
{
    std::vector<index_t> ext(d, extension);

    const gsVector<> pointErrors = computeErrors(thbMap, params, points);
    const real_t maxError = *std::max_element(pointErrors.begin(), pointErrors.end());

    if (percentage > 0)
    {
        gsInfo << "Adaptive refinement using relative threshold" << std::endl;
        std::vector<index_t> boxes = getBoxesRelative<d>(thbMap.basis(), params, pointErrors, percentage, ext);
        thbMap.basis().refineElements_withCoefs(thbMap.coefs(), boxes);
        return true;
    }
    else
    {
        gsInfo << "Adaptive refinement using absolute threshold" << std::endl;
        if (threshold <= maxError)
        {
            std::vector<index_t> boxes = getBoxes<d>(thbMap.basis(), params, pointErrors, threshold, ext);
            thbMap.basis().refineElements_withCoefs(thbMap.coefs(), boxes);
            return true;
        }
        else
        {
            return false;
        }
    }
}

template<unsigned dim>
void runComputations(std::string tempFile, std::string geomFile, std::string mapFile,
                     int deg, int numInnerKnots, int numIter,
                     int numPointsTarget, int numPoints,
                     real_t weightPDM, real_t weightTDM,
                     bool pcorrect,
                     real_t lambdaTikhonov, real_t lambdaSmoothing,
                     bool adaptiveRef, int numRefIter, int extension,
                     real_t threshold, real_t percentage,
                     int numJacPts, int numPlotPts,
                     real_t knotMin, real_t knotMax,
                     bool solveUpdates, std::string outputPrefix, bool debug)
{

    const gsMultiPatch<>* templateM = readMultipatch(tempFile);
    const gsMultiPatch<>& templateMultipatch = *templateM;
    const gsMultiPatch<>* targetM   = readMultipatch(geomFile);
    const gsMultiPatch<>& targetMultipatch = *targetM;

    // Constructing initial basis
    // If initial map in the file is not given,
    // then knotMIn, knotMax, numInnerKnots and deg are expected to be given
    std::vector< gsKnotVector<> > knotVectors;
    gsKnotVector<> kv(knotMin, knotMax, numInnerKnots, deg+1);

    for (unsigned i = 0; i < dim; i++)
    {
        knotVectors.push_back(kv);
    }
    gsTensorBSplineBasis<dim> initBSplineBasis( knotVectors );
    gsTHBSplineBasis<dim> initTHBSplineBasis (initBSplineBasis);
    gsMatrix<> coefs(initTHBSplineBasis.size(), dim);
    coefs.setZero();
    gsTHBSpline<dim> thbMap(initTHBSplineBasis, coefs);

    const int templateNumPatches = templateMultipatch.nPatches();
    const int targetNumPatches   = targetMultipatch.nPatches();

    if ( templateNumPatches != targetNumPatches)
    {
        gsInfo << "Template geometry and target geometry multipatch files should have the same numbers of patches (curves).\n";
    }

    gsWriteParaview(templateMultipatch, outputPrefix + "_template_multipatch");
    gsWriteParaview(targetMultipatch, outputPrefix + "_target_multipatch");

    gsDistanceMinimization<dim> distanceMinimizer(thbMap, initTHBSplineBasis, templateMultipatch, targetMultipatch);

    distanceMinimizer.set_debug(debug);
    distanceMinimizer.set_numJacPts(numJacPts);
    distanceMinimizer.set_numPlotPts(numPlotPts);
    // Compute target points (# = numPointsTarget) and template points (# = numPoints)
    distanceMinimizer.computeTargetAndTemplatePoints(numPoints, numPointsTarget);

    if ( mapFile == "" )
    {
        // If map file is not given, the map is computed using PDM
        distanceMinimizer.computeInitialMap(initTHBSplineBasis, lambdaTikhonov, lambdaSmoothing, numIter, pcorrect, outputPrefix);
    }
    else
    {
        // Otherwise, the map is uploaded from the file
        distanceMinimizer.uploadInitialMap(mapFile, numIter, pcorrect, outputPrefix);
    }

    // Computing in the initial spline space (before refinement)
    const std::string prefix = outputPrefix + "Ref" + util::to_string(0);
    gsInfo << "Initial spline space DOFs: " << thbMap.coefsSize() << std::endl;
    distanceMinimizer.compute(numIter, weightPDM, weightTDM, lambdaTikhonov, lambdaSmoothing, pcorrect, solveUpdates, prefix);

    // Refinement
    for (int r = 0; r != numRefIter; r++)
    {
        gsInfo << "Refinement: " << r + 1 << " / " << numRefIter << std::endl;
        bool refined = false;
        if(adaptiveRef) // Adaptive refinement
        {
            refined = refine<dim>(thbMap, distanceMinimizer.getTemplatePoints(), distanceMinimizer.getSampledTargetPoints(),
                                  extension, threshold, percentage);
            gsInfo << "Adaptive refinement DOFs: " << thbMap.coefsSize() << std::endl;
        }
        else // Uniform refinement
        {
            refined = uniformRefine<dim>(thbMap);
            gsInfo << "Uniform refinement DOFs: " << thbMap.coefsSize() << std::endl;
        }
        if (!refined)
        {
            gsInfo << "Error is below the threshold (" << threshold << ")" << std::endl;
        }
        gsInfo << std::endl;

        const std::string prefix = outputPrefix + "Ref" + util::to_string(r+1);
        distanceMinimizer.compute(numIter, weightPDM, weightTDM, lambdaTikhonov, lambdaSmoothing, pcorrect, solveUpdates, prefix);
    }
}


int main(int argc, char *argv[])
{

    // Options with default values
    // Name of template geometry multipatch file (input)
    std::string tempFile(MOTOR_DATA_DIR "jku/tdm_template_multipatch.xml");

    // Name of target geometry multipatch file (input)
    std::string geomFile(MOTOR_DATA_DIR "jku/tdm_target_multipatch.xml");

    // Name of the initial map file (input)
    std::string mapFile(MOTOR_DATA_DIR "jku/tdm_thb_map.xml");

    // Prefix for all output files
    std::string outputPrefix("results");

    // Weight for PDM terms in the objective function
    real_t weightPDM = 0.0;

    // Weight for TDM terms in the objective function
    real_t weightTDM = 1.0;

   // Regularization parameter
    real_t lambdaTikhonov  = 1.0e-2;
    // Smoothing parameter
    real_t lambdaSmoothing = 0.0;

    // Number of sampling points (on the template geometry)
    int numPoints        = 1000;
    // Number of fitting points (on the target geometry)
    int numPointsTarget  = 100;
    // Degree of the map
    int deg              = 2;
    // Number of inner knots
    int numInnerKnots    = 7;
    // The value of the first knot (should be given if the initial map is not given in file only)
    real_t knotMin      = 0.0;
    // The value of the last knot (should be given if the initial map is not given in file only)
    real_t knotMax      = 1.0;
    // Number of fitting procedure iterations
    int numIter          = 5;
    // Number of refinement procedure iterations
    int numRefIter       = 5;
    bool adaptiveRef     = false;
    real_t threshold     = 1e-5;
    real_t percentage    = 0.0;
    int numJacPts        = 1000000;
    int numPlotPts       = 1000;
    int extension = 1;

    bool pcorrect = false;
    bool solveUpdates = false;
    bool debug    = false;

    gsCmdLine cmd("Constructs a map which transforms (maps) the template shape to the target shape. The map is constructed using point distance minimization (PDM) or tangent distance minimization (TDM) procedures.");
    cmd.addSwitch("debug", "If debug, it outputs more data.", debug);
    cmd.addString("o", "output", "Prefix for all output files", outputPrefix);
    cmd.addSwitch("", "updates", "Solve incremental updates of control points, instead of control points directly", solveUpdates);
    cmd.addReal("", "kmax", "The value of the last knot (should be given if the initial map is not given in file only)", knotMax);
    cmd.addReal("", "kmin", "The value of the first knot (should be given if the initial map is not given in file only)", knotMin);
    cmd.addInt("", "numPlotPts", "Number of points to Paraview", numPlotPts);
    cmd.addInt("", "numJacPts", "Number of points to check Jacobian", numJacPts);
    cmd.addReal("", "percentage", "Percentage (refinement)", percentage);
    cmd.addReal("t", "threshold", "Error threshold (refinement)", threshold);
    cmd.addInt("x", "extension", "Extension of the refinement", extension);
    cmd.addInt("r", "numRefIter", "Number of adaptive (local) refinement procedure iterations", numRefIter);
    cmd.addSwitch("a", "adaptiveRef", "Apply adaptive refinement", adaptiveRef);
    cmd.addReal("", "lambdaS", "Smoothing parameter", lambdaSmoothing);
    cmd.addReal("", "lambdaT", "Regularization parameter", lambdaTikhonov);
    cmd.addSwitch("p", "pcorrect", "Apply the foot point computation (parameter optimization)", pcorrect);
    cmd.addReal("", "weightTDM", "Weight for TDM terms in the objective function", weightTDM);
    cmd.addReal("", "weightPDM", "Weight for PDM terms in the objective function", weightPDM);
    cmd.addInt("n", "numPts", "Number of points sampled on the template shape (used for the foot point computation)", numPoints);
    cmd.addInt("m", "numPtsT", "Number of fitting points on the target shape", numPointsTarget);
    cmd.addInt("i", "numIter", "Number of fitting procedure iterations", numIter);
    cmd.addInt("k", "numKnots", "Number of interior knots", numInnerKnots);
    cmd.addInt("d", "degree", "Degree", deg);
    cmd.addString("M", "mapFile", "File containing initial map (optional)", mapFile);
    cmd.addString("G", "geometryFile", "File containing target shape multi-patch", geomFile);
    cmd.addString("T", "templateFile", "File containing template shape multi-patch", tempFile);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments: \n"
           << "Template file:     " << tempFile << "\n"
           << "Geometry file:     " << geomFile << "\n"
           << "Map file:          " << mapFile << "\n"
           << "Par. correction:   " << pcorrect << "\n"
           << "solveUpdates   :   " << solveUpdates << "\n"
           << "weightPDM:         " << weightPDM << "\n"
           << "weightTDM:         " << weightTDM << "\n"
           << "lambdaTikhonov:    " << lambdaTikhonov << "\n"
           << "lambdaSmoothing:   " << lambdaSmoothing << "\n"
           << "numPoints:         " << numPoints << "\n"       // Large number
           << "numPointsTarget:   " << numPointsTarget << "\n" // numPointsTarget << numPoints
           << "deg:               " << deg << "\n"
           << "numInnerKnots:     " << numInnerKnots << "\n"
           << "knotMin:           " << knotMin << "\n"
           << "knotMax:           " << knotMax << "\n"
           << "numIter:           " << numIter << "\n"
           << "numRefIter:        " << numRefIter << "\n"
           << "adaptiveRef:       " << adaptiveRef << "\n"
           << "threshold:         " << threshold << "\n"
           << "percentage:        " << percentage << "\n"
           << "extension:         " << extension << "\n"
           << "numJacPts:         " << numJacPts << "\n"
           << "numPlotPts:        " << numPlotPts << "\n"
           << "output:            " << outputPrefix << "\n"
           << "Debug:             " << debug << "\n"
           << "--------------------------------------------------\n" << std::endl;

    const gsMultiPatch<>* targetM   = readMultipatch(geomFile);
    const gsMultiPatch<>& targetMultipatch = *targetM;
    if (targetMultipatch[0].basis().dim() == 1)
    {
        runComputations<2>(tempFile, geomFile, mapFile,
                           deg, numInnerKnots, numIter,
                           numPointsTarget, numPoints,
                           weightPDM, weightTDM,
                           pcorrect,
                           lambdaTikhonov, lambdaSmoothing,
                           adaptiveRef, numRefIter, extension,
                           threshold, percentage,
                           numJacPts, numPlotPts,
                           knotMin, knotMax,
                           solveUpdates, outputPrefix, debug);
    }
    else if (targetMultipatch[0].basis().dim() == 2)
    {
        runComputations<3>(tempFile, geomFile, mapFile,
                           deg, numInnerKnots, numIter,
                           numPointsTarget, numPoints,
                           weightPDM, weightTDM,
                           pcorrect,
                           lambdaTikhonov, lambdaSmoothing,
                           adaptiveRef, numRefIter, extension,
                           threshold, percentage,
                           numJacPts, numPlotPts,
                           knotMin, knotMax,
                           solveUpdates, outputPrefix, debug);
    }

    

    return 0;
}

// Assembling fitting matrix
template <unsigned dim>
gsSparseMatrix<> assembleMatrix(const gsMatrix<>& basisValues,
                          const gsMatrix<index_t>& basisActives,
                          const real_t weightPDM,
                          const real_t weightTDM,
                          const gsMatrix<>& nrmls,
                          const gsMatrix<>& coefs)
{
    //const int M = basisValues.cols();
    const int N = coefs.rows();
    gsSparseMatrix<> A(dim*N, dim*N);
    A.setZero();

    if (weightPDM > 0)
    {
        assembleMatrixPDM<dim>(A, basisValues, basisActives, weightPDM);
    }
    if (weightTDM > 0)
    {
        assembleMatrixTDM<dim>(A, basisValues, basisActives, weightTDM, nrmls);
    }

    return A;
}

// Assembling rhs of the linear system
// (when control points are calculating)
template <unsigned dim>
gsMatrix<> assembleRHS(const gsMatrix<>& basisValues,
                       const gsMatrix<index_t>& basisActives,
                       const real_t weightPDM,
                       const real_t weightTDM,
                       const gsMatrix<>& sampledPtsTargetCurve,
                       const gsMatrix<>& nrmls,
                       const gsMatrix<>& coefs)
{
    //const int M = basisValues.cols();
    const int N = coefs.rows();
    gsMatrix<> rhs(dim*N, 1);
    rhs.setZero();

    if (weightPDM > 0)
    {
        assembleRHSPDM<dim>(rhs, basisValues, basisActives, weightPDM, sampledPtsTargetCurve);
    }
    if (weightTDM > 0)
    {
        assembleRHSTDM<dim>(rhs, basisValues, basisActives, weightTDM, sampledPtsTargetCurve, nrmls);
    }

    return rhs;
}

// Assembling rhs of the linear system
// (when control points updates are calculating)
template <unsigned dim>
gsMatrix<> assembleRHS(const gsMatrix<>& basisValues,
                       const gsMatrix<index_t>& basisActives,
                       const real_t weightPDM,
                       const real_t weightTDM,
                       const gsMatrix<>& sampledPtsTargetCurve,
                       const gsMatrix<>& targetPts,
                       const gsMatrix<>& nrmls,
                       const gsMatrix<>& coefs)
{
    //const int M = basisValues.cols();
    const int N = coefs.rows();
    gsMatrix<> rhs(dim*N, 1);
    rhs.setZero();

    if (weightPDM > 0)
    {
        assembleRHSPDMSolveUpdates<dim>(rhs, basisValues, basisActives, weightPDM, sampledPtsTargetCurve, targetPts);
    }
    if (weightTDM > 0)
    {
        assembleRHSTDMSolveUpdates<dim>(rhs, basisValues, basisActives, weightTDM, sampledPtsTargetCurve, targetPts, nrmls);
    }

    return rhs;
}


// Fitting matrix: point distance minimization (PDM)
template <unsigned dim>
void assembleMatrixPDM(gsSparseMatrix<> &A,
                       const gsMatrix<>& basisValues,
                       const gsMatrix<index_t>& basisActives,
                       const real_t weightPDM)
{
    const int M = basisValues.cols();
    //const int N = coefs.rows();

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

                    for (unsigned d = 0; d < dim; d++)
                    {
                       A(dim*idBasis1+d, dim*idBasis2+d) += 2 * weightPDM * basisValue2 * basisValue1; //
                    }
                }
            }
        }
    }


    
}

// rhs: point distance minimization (PDM)
template <unsigned dim>
void assembleRHSPDM(gsMatrix<>& rhs,
                    const gsMatrix<>& basisValues,
                    const gsMatrix<index_t>& basisActives,
                    const real_t weightPDM,
                    const gsMatrix<>& sampledPtsTargetCurve)
{
    const int M = basisValues.cols();
    //const int N = coefs.rows();

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {
        const gsVector<> P = sampledPtsTargetCurve.col(idPoint);
        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
        {
            int idBasis = basisActives(idActive, idPoint);
            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) )
            {
                const real_t basisValue = basisValues(idActive, idPoint);
                for (unsigned d = 0; d < dim; d++)
                {
                    rhs(dim*idBasis+d,   0) += 2 * weightPDM * basisValue * P(d, 0);
                }
            }
        }
    }
}

template <unsigned dim>
void assembleRHSPDMSolveUpdates(gsMatrix<>& rhs,
                    const gsMatrix<>& basisValues,
                    const gsMatrix<index_t>& basisActives,
                    const real_t weightPDM,
                    const gsMatrix<>& sampledPtsTargetCurve,
                    const gsMatrix<>& targetPts)
{
    const int M = basisValues.cols();

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {
        const gsVector<> P = sampledPtsTargetCurve.col(idPoint) - targetPts.col(idPoint);
        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
        {
            int idBasis = basisActives(idActive, idPoint);
            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) )
            {
                const real_t basisValue = basisValues(idActive, idPoint);
                for (unsigned d = 0; d < dim; d++)
                {
                    rhs(dim*idBasis+d, 0) += 2 * weightPDM * basisValue * P(d, 0) ;
                }
            }
        }
    }
}

// Fitting matrix: tangent distance minimization (TDM)
template <unsigned dim>
void assembleMatrixTDM(gsSparseMatrix<> &A,
                       const gsMatrix<>& basisValues,
                       const gsMatrix<index_t>& basisActives,
                       const real_t weightTDM,
                       const gsMatrix<>& nrmls)
{
    const int M = basisValues.cols();
    //const int N = coefs.rows();

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

                    if (dim == 2)
                    {
                        const real_t nrml0 = nrmls(0, idPoint);
                        const real_t nrml1 = nrmls(1, idPoint);

                        A(2*idBasis1,   2*idBasis2)   += weightTDM * basisValue2 * nrml0 * basisValue1 * nrml0; //
                        A(2*idBasis1,   2*idBasis2+1) += weightTDM * basisValue2 * nrml1 * basisValue1 * nrml0; //

                        A(2*idBasis1+1, 2*idBasis2)   += weightTDM * basisValue2 * nrml0 * basisValue1 * nrml1; //
                        A(2*idBasis1+1, 2*idBasis2+1) += weightTDM * basisValue2 * nrml1 * basisValue1 * nrml1; //
                    }
                    else if (dim == 3)
                    {
                        const real_t nrml0 = nrmls(0, idPoint);
                        const real_t nrml1 = nrmls(1, idPoint);
                        const real_t nrml2 = nrmls(2, idPoint);

                        A(3*idBasis1,   3*idBasis2)   += weightTDM * basisValue2 * nrml0 * basisValue1 * nrml0; //
                        A(3*idBasis1,   3*idBasis2+1) += weightTDM * basisValue2 * nrml1 * basisValue1 * nrml0; //
                        A(3*idBasis1,   3*idBasis2+2) += weightTDM * basisValue2 * nrml2 * basisValue1 * nrml0; //

                        A(3*idBasis1+1, 3*idBasis2)   += weightTDM * basisValue2 * nrml0 * basisValue1 * nrml1; //
                        A(3*idBasis1+1, 3*idBasis2+1) += weightTDM * basisValue2 * nrml1 * basisValue1 * nrml1; //
                        A(3*idBasis1+1, 3*idBasis2+2) += weightTDM * basisValue2 * nrml2 * basisValue1 * nrml1; //

                        A(3*idBasis1+2, 3*idBasis2)   += weightTDM * basisValue2 * nrml0 * basisValue1 * nrml2; //
                        A(3*idBasis1+2, 3*idBasis2+1) += weightTDM * basisValue2 * nrml1 * basisValue1 * nrml2; //
                        A(3*idBasis1+2, 3*idBasis2+2) += weightTDM * basisValue2 * nrml2 * basisValue1 * nrml2; //
                    }

                }
            }
        }
    }
}

// rhs: tangent distance minimization (TDM)
template <unsigned dim>
void assembleRHSTDM(gsMatrix<>& rhs,
                    const gsMatrix<>& basisValues,
                    const gsMatrix<index_t>& basisActives,
                    const real_t weightTDM,
                    const gsMatrix<>& sampledPtsTargetCurve,
                    const gsMatrix<>& nrmls)
{
    const int M = basisValues.cols();
    //const int N = coefs.rows();

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {
        const gsVector<> P = sampledPtsTargetCurve.col(idPoint);
        const gsVector<> nrml = nrmls.col(idPoint);
        const real_t tmp = P.dot(nrml);

        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
        {
            int idBasis = basisActives(idActive, idPoint);
            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) )
            {
                const real_t basisValue = basisValues(idActive, idPoint);
                for (unsigned d = 0; d < dim; d++)
                {
                    rhs(dim*idBasis+d, 0) += weightTDM * tmp * basisValue * nrml(d);
                }
            }
        }
    }
}

template <unsigned dim>
void assembleRHSTDMSolveUpdates(gsMatrix<>& rhs,
                    const gsMatrix<>& basisValues,
                    const gsMatrix<index_t>& basisActives,
                    const real_t weightTDM,
                    const gsMatrix<>& sampledPtsTargetCurve,
                    const gsMatrix<>& targetPts,
                    const gsMatrix<>& nrmls)
{
    const int M = basisValues.cols();
    //const int N = coefs.rows();

    for (index_t idPoint = 0; idPoint < M; idPoint++) // Works as k
    {
        const gsVector<> P = sampledPtsTargetCurve.col(idPoint) - targetPts.col(idPoint);
        const gsVector<> nrml = nrmls.col(idPoint);
        const real_t tmp = P.dot(nrml);

        for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
        {
            int idBasis = basisActives(idActive, idPoint);
            if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) )
            {
                const real_t basisValue = basisValues(idActive, idPoint);
                for (unsigned d = 0; d < dim; d++)
                {
                    rhs(dim*idBasis+d, 0) += weightTDM * tmp * basisValue * nrml(d);
                }
            }
        }
    }
}

/// Applies Tikhonov regularization
template <unsigned dim>
void applyTikhonovRegularization(gsSparseMatrix<> &A_mat,
                                 gsMatrix<>& rhs,
                                 const real_t lambda,
                                 const gsMatrix<>& coefs)
{

    //const int M = basisValues.cols();

    for (index_t i = 0; i < A_mat.rows(); i++)
    {
        // Regularization terms are added to the diagonal of the matrix
        A_mat(i, i) += 2*lambda;
    }

    for (index_t i = 0; i < coefs.rows(); i++)
    {
        // Regularization terms are added to the diagonal of the matrix
        for (unsigned d = 0; d < dim; d++)
        {
            rhs(dim*i+d, 0) += 2*lambda*coefs(i, d);
        }
    }
}

/// Applies smoothing
void applySmoothing(gsSparseMatrix<> &A_mat,
                    gsMatrix<>& rhs,
                    real_t lambda,
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

/// Solves coefficients (control points)
template <unsigned dim>
gsMatrix<> solveCoefs(const gsSparseMatrix<> &A,
                      const gsMatrix<>& b,
                      const bool debug,
                      const std::string debPrefix)
{

    // (8.?) Configurating the linear solver
    gsSparseMatrix<> sparseA = A;//.sparseView();
    //sparseA.makeCompressed();

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

    gsMatrix<> tmpSolution = solver.solve(b);

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

    const int N = tmpSolution.rows()/dim;
    gsMatrix<> solution (N, dim);
    solution.setZero();

    for (int i = 0; i < N; i++)
    {
        for (unsigned d = 0; d < dim; d++)
        {
            solution(i, d) = tmpSolution(dim*i+d, 0);
        }
    }

    return solution;
}
/*
// wrong functional, fix it
real_t computeFunctional(const gsMatrix<>& sampledPtsTargetCurve,
                         const gsMatrix<>& templatePoints,
                         const gsGeometry<>& thbMap,
                         const gsMatrix<>& nrmls)
{
    gsMatrix<> targetPoints;
    thbMap.eval_into(templatePoints, targetPoints);

    return computeFunctionalTDM(sampledPtsTargetCurve, targetPoints, nrmls);
}

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
    return F;
}
*/

/// Computes the squared \f$ l_2 \f$-error
real_t computeSquaredError(const gsMatrix<>& dataPoints,
                            const gsMatrix<>& params,
                            const gsGeometry<>& thbMap)
{
    gsMatrix<> targetPoints;
    targetPoints.setZero();
    thbMap.eval_into(params, targetPoints);

    const int M = dataPoints.cols();

    real_t error = 0;

    for (index_t idPoint = 0; idPoint != M; idPoint++)
    {
        const gsVector<> P = dataPoints.col(idPoint);
        const gsVector<> targetPoint = targetPoints.col(idPoint);
        const real_t err = (P - targetPoint).squaredNorm();
        error += err;
    }

    return error;
}

/// Computes normals
template<unsigned dim>
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
    gsMatrix<> nrmls(dim, M);
    nrmls.setZero();

    gsMatrix<> derivMap_u(dim, M);
    derivMap_u.setZero();

    gsMatrix<> derivMap_v(dim, M);
    derivMap_v.setZero();

    if (dim == 2)
    {
        for (index_t idPoint = 0; idPoint < M; idPoint++)
        {
            for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
            {
                int idBasis = basisActives(idActive, idPoint);
                if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) ) //
                {
                        const real_t tmpDeriv = (derivMapBasis(dim*idActive,   idPoint)*derivTempBoundaryGeom(0, idPoint) +
                                                 derivMapBasis(dim*idActive+1, idPoint)*derivTempBoundaryGeom(1, idPoint));
                        nrmls(0, idPoint) += -coefs(idBasis, 1)*tmpDeriv;
                        nrmls(1, idPoint) +=  coefs(idBasis, 0)*tmpDeriv;
                }
            }
        }
    }
    else if (dim == 3)
    {
        for (index_t idPoint = 0; idPoint < M; idPoint++)
        {
            for (index_t idActive = 0; idActive < basisActives.rows(); idActive++)
            {
                int idBasis = basisActives(idActive, idPoint);
                if ( (idBasis > 0) || ((idBasis == 0) && (idActive == 0)) ) //
                {
                    const real_t tmpDeriv_u = (derivMapBasis(dim*idActive,   idPoint)*derivTempBoundaryGeom(0, idPoint) +
                                               derivMapBasis(dim*idActive+1, idPoint)*derivTempBoundaryGeom(2, idPoint) +
                                               derivMapBasis(dim*idActive+2, idPoint)*derivTempBoundaryGeom(4, idPoint));
                    const real_t tmpDeriv_v = (derivMapBasis(dim*idActive,   idPoint)*derivTempBoundaryGeom(1, idPoint) +
                                               derivMapBasis(dim*idActive+1, idPoint)*derivTempBoundaryGeom(3, idPoint) +
                                               derivMapBasis(dim*idActive+2, idPoint)*derivTempBoundaryGeom(5, idPoint));

                    derivMap_u(0, idPoint) += coefs(idBasis, 0)*tmpDeriv_u;
                    derivMap_u(1, idPoint) += coefs(idBasis, 1)*tmpDeriv_u;
                    derivMap_u(2, idPoint) += coefs(idBasis, 2)*tmpDeriv_u;

                    derivMap_v(0, idPoint) += coefs(idBasis, 0)*tmpDeriv_v;
                    derivMap_v(1, idPoint) += coefs(idBasis, 1)*tmpDeriv_v;
                    derivMap_v(2, idPoint) += coefs(idBasis, 2)*tmpDeriv_v;
                 }
            }
            nrmls(0, idPoint) = derivMap_u(1, idPoint)*derivMap_v(2, idPoint) - derivMap_u(2, idPoint)*derivMap_v(1, idPoint);
            nrmls(1, idPoint) = derivMap_u(2, idPoint)*derivMap_v(0, idPoint) - derivMap_u(0, idPoint)*derivMap_v(2, idPoint);
            nrmls(2, idPoint) = derivMap_u(0, idPoint)*derivMap_v(1, idPoint) - derivMap_u(1, idPoint)*derivMap_v(0, idPoint);
        }
    }
    normalize(nrmls);

    return nrmls;
}

/// Saves data
void saveData(const gsGeometry<>& thbMap,
              const gsMultiPatch<>& temp,
              const std::string outputPrefix)
{
    gsMultiPatch<> newMp = mapMultipatch(thbMap, temp);

    writeMultipatch(newMp, outputPrefix);
    writeKnots(thbMap, outputPrefix);
}
