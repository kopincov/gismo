/** @file uwbINSAssemblerBase.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include "uwbINSBlockAssembler.h"
#include "uwbINSSolverParams.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerBase
{
public:
    uwbINSAssemblerBase(uwbINSSolverParams<T>& params) :
        m_blockAssembler(params)
    {
        m_bCavity = params.settings().get(constantsINS::cavity);
    }

    virtual ~uwbINSAssemblerBase()
    {
    }

protected:
    void initMembers() { GISMO_NO_IMPLEMENTATION }

    virtual void reinitMembers() { initMembers(); }

public:
    virtual void initialize() { GISMO_NO_IMPLEMENTATION }

    virtual void update(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_solution = solVector;
        m_blockAssembler.updateCurrentSolField(solVector);

        updateAssembly();

        fillMatrix();
        fillRhs();

        if (m_bCavity)
            applyDirichletConditions();
    }

    virtual void updatePicard(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_blockAssembler.updateCurrentSolField(solVector);

        updatePicardAssembly(solVector);

        fillMatrix();
        fillRhs();

        if (m_bCavity)
            applyDirichletConditions();
    }

    virtual void fillStokesSystem_into(gsSparseMatrix<T> & stokesMatrix, gsMatrix<T> & stokesRhs) const
    {
        m_blockAssembler.fillStokesSystem_into(stokesMatrix, stokesRhs);
    }

    virtual void constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>& result, gsVector<index_t> patchNumbers, std::vector<boxSide> sides, int unk, bool relative = false) const
    {
        m_blockAssembler.constructSolution(solVector, result, patchNumbers, sides, unk, relative);
    }

    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk, bool relative = false) const
    {
        return m_blockAssembler.constructSolution(solVector, unk, relative);
    }

    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk, gsVector<index_t> patchNumbers, std::vector<boxSide> sides, bool relative = false) const
    {
        return m_blockAssembler.constructSolution(solVector, unk, patchNumbers, sides, relative);
    }

    virtual gsField<T> constructSolutionCombined(const gsMatrix<T>& solVector, gsVector<size_t> relPatches) const
    {
        return m_blockAssembler.constructSolutionCombined(solVector, relPatches);
    }

    virtual T computeFlowRate(int patch, boxSide side, gsMatrix<T> solution) const
    {
        return m_blockAssembler.computeFlowRate(patch, side, solution);
    }

    virtual T computeDimensionlessWallDistance(gsMatrix<T> solution, gsVector<int> distancePatches, std::vector<boxSide> distanceSides, T viscosity, T reynoldsNumber, T uFreeStream, T maxYplus = 1.0, unsigned npts = 20, bool print = false, bool estimate = true) const
    {
        return m_blockAssembler.computeDimensionlessWallDistance(solution, distancePatches, distanceSides, viscosity, reynoldsNumber, uFreeStream, maxYplus, npts, print, estimate);
    }

    T computeDimensionlessWallDistance(gsMatrix<T> solution,  std::vector< gsMultiBasis<T> > basis, gsVector<int> distancePatches, std::vector<boxSide> distanceSides, T viscosity, T reynoldsNumber, T uFreeStream, T maxYplus = 1.0, unsigned npts = 20, bool print = false, bool estimate = true) const
    {
        return m_blockAssembler.computeDimensionlessWallDistance(solution, basis, distancePatches, distanceSides, viscosity, reynoldsNumber, uFreeStream, maxYplus, npts, print, estimate);
    }

    T computeAspectRatio(bool minAR = false) { return m_blockAssembler.computeAspectRatio(minAR); }

    void addPressureOutletCondition(int patch, boxSide side)
    {
        gsMatrix<unsigned> boundaryIndicesOutlet = m_blockAssembler.getBases().back().basis(patch).boundary(side);

        std::vector< gsMatrix< index_t > > boundaryDofsToEliminate;
        boundaryDofsToEliminate.resize(m_blockAssembler.getBases().back().nBases());
        for (size_t i = 0; i < boundaryDofsToEliminate.size(); ++i)
        {
            boundaryDofsToEliminate[i].setZero(0, 0);
            if ( i == static_cast<size_t>(patch) )
                boundaryDofsToEliminate[i] = boundaryIndicesOutlet;
        }

        m_blockAssembler.markDofsAsEliminatedZeros(boundaryDofsToEliminate, 1);

        reinitMembers();
    }

    void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk = 0)
    {
        m_blockAssembler.markDofsAsEliminatedZeros(boundaryDofs, unk);

        reinitMembers();
    }

    virtual const gsSparseMatrix<T> & matrix() const { GISMO_NO_IMPLEMENTATION }

    virtual const gsMatrix<T> & rhs() const { GISMO_NO_IMPLEMENTATION }

    virtual const gsSparseMatrix<T>& getVelocityMassMatrix()
    { 
        return m_blockAssembler.getBlockM();
    }

    virtual const gsSparseMatrix<T>& getPressureMassMatrix()
    {
        return m_blockAssembler.getBlockMp();
    }

    void preparePCDboundary(std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, int bcType = 0)
    {
        m_blockAssembler.preparePCDboundary(bndIn, bndOut, bndWall, bcType);
    }

    virtual void fillPCDblocks(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, int bcType, bool assembAp, bool assembFp, bool lumping)
    {
        Ap = m_blockAssembler.getPressurePoissonMatrix(assembAp, lumping);

        gsSparseMatrix<T> ApforFp;
        if (assembAp == assembFp)
            ApforFp = Ap;
        else
            ApforFp = m_blockAssembler.getPressurePoissonMatrix(assembFp, lumping);

        if(Fp.nonZeros()) // the matrix is not empty
            Fp += getViscosity() * ApforFp + m_blockAssembler.getPressureConvectionMatrix();
        else 
            Fp = getViscosity() * ApforFp + m_blockAssembler.getPressureConvectionMatrix();

        m_blockAssembler.applyPCDboundaryConditions(Ap, Fp, bcType);
    }

    void getPCDblocks(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, int bcType, bool assembAp, bool assembFp, bool lumping)
    {
        m_blockAssembler.preparePCDboundary(bndIn, bndOut, bndWall, bcType);
        this->fillPCDblocks(Ap, Fp, bcType, assembAp, assembFp, lumping);
    }

    virtual void evalElWiseForLocRef(const gsMatrix<T> & solVector, std::vector<gsVector<T> > & elWiseVals, bool outputInQuadPoints = false)
    { GISMO_NO_IMPLEMENTATION }

    void evalElWisePhysQuadPoints(std::vector<gsMatrix<T> > & elWiseQPts)
    {
        m_blockAssembler.evaluatePhysQuadraturePoints(elWiseQPts);
    }

protected:

    virtual void updateAssembly()
    {
        m_blockAssembler.assembleNonlinearPart(m_solution);
    }

    virtual void updatePicardAssembly(const gsMatrix<T> & solVector)
    { GISMO_NO_IMPLEMENTATION }

    virtual void fillBase() { GISMO_NO_IMPLEMENTATION }

    virtual void fillMatrix() { GISMO_NO_IMPLEMENTATION }

    virtual void fillRhs() { GISMO_NO_IMPLEMENTATION }

public:
    void plot2DVorticity(std::string const & fn, gsMatrix<T>& solution, unsigned npts)
    {
        gsParaviewCollection collection(fn);
        std::string fileName;

        gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);

        for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
        {
            fileName = fn + util::to_string(p);

            gsMatrix<T> geoVals, vorticityVals;
            gsVector<unsigned> np;

            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
            const int parDim = patch->domainDim();
            const int tarDim = patch->targetDim();

            gsMatrix<T> ab = patch->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);
            np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts = gsPointGrid(a, b, np);

            evalVorticity_singlePatch(p, pts, uSolField, vorticityVals);

            // Evaluate geometry at pts
            geoVals = patch->eval(pts);

            if (3 - parDim > 0)
            {
                np.conservativeResize(3);
                np.bottomRows(3 - parDim).setOnes();
            }
            if (3 - tarDim > 0)
            {
                geoVals.conservativeResize(3, geoVals.cols());
                geoVals.bottomRows(3 - tarDim).setZero();
            }

            gsWriteParaviewTPgrid(geoVals, vorticityVals, np.template cast<index_t>(), fileName);

            collection.addPart(fileName + ".vts");
        }
        collection.save();
    }

    void plotPressureCoefficient(std::string const & fn, gsMatrix<T>& solution,
                                   index_t referencePatchIndex, gsVector<T>& referencePoint, T freeStreamVelocity,
                                   T rho, unsigned npts)
    {
        gsParaviewCollection collection(fn);
        std::string fileName;

        //gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);
        gsField<T> pSolField = m_blockAssembler.constructSolution(solution, 1);

        gsMatrix<T> referenceSolPVal = pSolField.value(referencePoint, referencePatchIndex);

        for (unsigned int patchIndex = 0; patchIndex < m_blockAssembler.getPatches().nPatches(); ++patchIndex)
        {
            fileName = fn + util::to_string(patchIndex);

            gsMatrix<T> geoVals, pressCoeffVals;
            gsVector<unsigned> np;

            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(patchIndex);
            const int parDim = patch->domainDim();
            const int tarDim = patch->targetDim();

            gsMatrix<T> ab = patch->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);
            np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts = gsPointGrid(a, b, np);

            //----------
            //evalPressureCoeff_singlePatch(patchIndex, pts, uSolField, pSolField, pressCoeffVals);
            //gsMatrix<T> solUVals = uSolField.value(pts, patchIndex);
            gsMatrix<T> solPVals = pSolField.value(pts, patchIndex);

            index_t nEvalPoints = pts.cols();
            gsVector<T> pressureCoeffVals(nEvalPoints);
            // Evaluate pressure coefficient at pts
            for (int k = 0; k < nEvalPoints; k++)
            {
                T pCoeffK = (solPVals.coeff(k, 0) - referenceSolPVal(0, 0)) / rho * math::pow(freeStreamVelocity, 2);
                //T pCoeffK = (solPVals.coeff(k, 0) - referenceSolPVal(0, 0)) / rho * math::pow((solUVals.col(k)).norm(), 2);
                //T pCoeffK = solPVals.coeff(k, 0) / rho * math::pow((solUVals.col(k)).norm(), 2);
                pressureCoeffVals(k) = pCoeffK;
            }

            pressCoeffVals = pressureCoeffVals.transpose();

            //----------

            // Evaluate geometry at pts
            geoVals = patch->eval(pts);

            if (3 - parDim > 0)
            {
                np.conservativeResize(3);
                np.bottomRows(3 - parDim).setOnes();
            }
            if (3 - tarDim > 0)
            {
                geoVals.conservativeResize(3, geoVals.cols());
                geoVals.bottomRows(3 - tarDim).setZero();
            }

            gsWriteParaviewTPgrid(geoVals, pressCoeffVals, np.template cast<index_t>(), fileName);

            collection.addPart(fileName + ".vts");
        }
        collection.save();
    }

protected:
    void evalVorticity_singlePatch(index_t patchIndex, gsMatrix<T>& pts, const gsField<T>& uSolField, gsMatrix<T>& vortVals)
    {
        gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(patchIndex);
        gsBasisRefs<T> bases(m_blockAssembler.getBases(), patchIndex);

        const int tarDim = patch->targetDim();

        unsigned evFlags = NEED_GRAD_TRANSFORM;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, *patch));

        gsMatrix<index_t> actives;
        gsMatrix<T> parGrads, physGrad;
        bases.front().deriv_into(pts, parGrads);
        geoEval->evaluateAt(pts);

        gsMatrix<T> actCoeffsU;
        gsMatrix<T> uGrads;
        index_t nQuPoints = pts.cols();
        gsVector<T> vorticityVals(nQuPoints);
        // Evaluate vorticity at pts
        for (int k = 0; k < nQuPoints; k++)
        {
            // eval uGrads at all pts
            bases.front().active_into(pts.col(k), actives);
            int numActU = actives.rows();
            actCoeffsU.setZero(tarDim, numActU);

            for (int j = 0; j < numActU; j++)
                actCoeffsU.col(j) = uSolField.coefficientVector(patchIndex).row(actives(j)).transpose();

            geoEval->transformGradients(k, parGrads, physGrad);
            uGrads.noalias() = actCoeffsU * physGrad.transpose();

            T vorK = uGrads.coeff(1, 0) - uGrads.coeff(0, 1);
            vorticityVals(k) = vorK;
        }

        vortVals = vorticityVals.transpose();
    }
    
public:
    void setSolution(gsMatrix<T> solVector)
    {
        m_solution = solVector;
        m_blockAssembler.updateCurrentSolField(solVector, true);
    }

public:
    void plotShearStress(std::string const & fn, gsMatrix<T> & solution, unsigned npts)
    {
        gsParaviewCollection collectionM(fn);
        std::string fileNameM;

        gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);

        for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
        {
            fileNameM = fn + util::to_string(p);

            gsMatrix<T> geoVals, shearStressVals;
            gsVector<unsigned> np;

            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
            const int parDim = patch->domainDim();
            const int tarDim = patch->targetDim();

            gsMatrix<T> ab = patch->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);
            np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts = gsPointGrid(a, b, np);

            evalShearStress_singlePatch(p, pts, uSolField, shearStressVals);

            // Evaluate geometry at pts
            geoVals = patch->eval(pts);

            if (3 - parDim > 0)
            {
                np.conservativeResize(3);
                np.bottomRows(3 - parDim).setOnes();
            }
            if (3 - tarDim > 0)
            {
                geoVals.conservativeResize(3, geoVals.cols());
                geoVals.bottomRows(3 - tarDim).setZero();
            }

            gsWriteParaviewTPgrid(geoVals, shearStressVals, np.template cast<index_t>(), fileNameM);

            collectionM.addPart(fileNameM + ".vts");
        }
        collectionM.save();
    }

protected:
    void evalShearStress_singlePatch(index_t patchIndex, gsMatrix<T>& pts, gsField<T>& uSolField, gsMatrix<T>& shearStressVals)
    {
        gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(patchIndex);
        gsBasisRefs<T> bases(m_blockAssembler.getBases(), patchIndex);

        const int tarDim = patch->targetDim();

        unsigned evFlags = NEED_VALUE | NEED_GRAD_TRANSFORM;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, *patch));

        gsMatrix<index_t> activesU;

        gsMatrix<T> basisGradsU, physGradU;
        bases.front().deriv_into(pts, basisGradsU);

        geoEval->evaluateAt(pts);

        index_t nQuPoints = pts.cols();

        gsMatrix<T> actCoeffsU;
        gsMatrix<T> uGrads;
        gsMatrix<T> tauWvals;
        tauWvals.setZero(nQuPoints);
        // Evaluate shear stress at pts
        for (int k = 0; k < nQuPoints; k++)
        {
            // eval uGrads at all pts
            bases.front().active_into(pts.col(k), activesU);
            int numActU = activesU.rows();
            actCoeffsU.setZero(tarDim, numActU);

            for (int j = 0; j < numActU; j++)
                actCoeffsU.col(j) = uSolField.coefficientVector(patchIndex).row(activesU(j)).transpose();

            geoEval->transformGradients(k, basisGradsU, physGradU);
            uGrads.noalias() = actCoeffsU * physGradU.transpose();

            tauWvals(k) = uGrads(0, 1);//;uGrads.row(var);
        }

        shearStressVals = tauWvals.transpose();
    }

public:
    void saveCfAtWall(std::string const & fn, gsMatrix<T> & solution, int patchIndex, gsVector<T>& referencePoint, unsigned npts)
    {
        gsFileData<> fd;

        std::string fileNameM;

        gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);
        gsField<T> pSolField = m_blockAssembler.constructSolution(solution, 1);

        //gsField<T> uSolFieldParam = m_blockAssembler.constructSolutionParam(solution, 0);
        //gsField<T> pSolFieldParam = m_blockAssembler.constructSolutionParam(solution, 1);

        int referencePatch = 2;
        gsMatrix<T> referencePoint_par;
        const gsMultiPatch<T> patches(m_blockAssembler.getPatches());
        patches.patch(referencePatch).invertPoints(referencePoint, referencePoint_par);

        gsMatrix<T> referenceSolUVal = uSolField.value(referencePoint_par, referencePatch);
        gsMatrix<T> referenceSolPVal = pSolField.value(referencePoint_par, referencePatch);
        T refSolPVal = referenceSolPVal(0,0);

/*        gsInfo << "U reference value = \n" << uSolField.value(referencePoint, referencePatch) << "\n";
        gsInfo << "P reference value = " << pSolField.value(referencePoint, referencePatch) << "\n";
        gsInfo << "U reference value = \n" << uSolField.pvalue(referencePoint, referencePatch) << "\n";
        gsInfo << "P reference value = " << pSolField.pvalue(referencePoint, referencePatch) << "\n";
        gsInfo << "U reference value = \n" << uSolFieldParam.value(referencePoint, referencePatch) << "\n";
        gsInfo << "P reference value = " << pSolFieldParam.value(referencePoint, referencePatch) << "\n";
        gsInfo << "U reference value = \n" << uSolFieldParam.pvalue(referencePoint, referencePatch) << "\n";
        gsInfo << "P reference value = " << pSolFieldParam.pvalue(referencePoint, referencePatch) << "\n";
*/

        gsInfo << "\n\n====================================\n\n";
        gsInfo << "U reference value = \n" << referenceSolUVal << "\n";
        gsInfo << "P reference value = " << referenceSolPVal << "\n";
        //gsInfo << "referenceSolUVal(0,0) = " << referenceSolUVal(0,0) << "\n";
        gsInfo << "\n====================================\n\n";

        //m_blockAssembler.getPatches().patch(patchIndex).eval_into(referencePoint, referencePoint_phys);

        //gsInfo << "referencePoint = \n" << referencePoint << "\n";
        //gsInfo << "referencePoint_par = \n" << referencePoint_par << "\n";

        int referencePatchTrans = 0;
        gsVector<T>& refPointTrans = referencePoint;
        refPointTrans(0) = 40*0.0127;
        refPointTrans(1) = 0;
        gsMatrix<T> refPointTrans_par;
        patches.patch(referencePatchTrans).invertPoints(refPointTrans, refPointTrans_par);
        gsMatrix<T> referenceSolPValTrans = pSolField.value(refPointTrans_par, referencePatchTrans);
        T refSolPValTrans = referenceSolPValTrans(0,0);
        refSolPValTrans = (refSolPValTrans - refSolPVal) / (0.5 * referenceSolUVal(0,0) * referenceSolUVal(0,0));
        gsInfo << "refPointTrans = " << refPointTrans << "\n";
        gsInfo << "refPointTrans_par = " << refPointTrans_par << "\n";
        gsInfo << "referenceSolPValTrans = " << referenceSolPValTrans << "\n";
        gsInfo << "refSolPValTrans = " << refSolPValTrans << "\n";


        fileNameM = fn + ".csv";

        std::ofstream file;
        file.open(fileNameM);
        std::ofstream fileCp;
        fileCp.open("Cp.csv");
        //-----------------------------------------------------------------

        for(int i = 0; i < 1; i++)
        {
            int p = 2;

            gsMatrix<T> geoVals, shearStressVals, cpVals;
            gsVector<unsigned> np;

            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
            //const int parDim = patch->domainDim();
            //const int tarDim = patch->targetDim();

            boxSide side = boundary::south;

            typename gsGeometry<T>::uPtr geomBoundary = m_blockAssembler.getPatches().patch(p).boundary(side);

            gsMatrix<T> ab = geomBoundary->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);
            np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts_pom = gsPointGrid(a, b, np);
            gsMatrix<T> pts;
            pts.setZero(2,pts_pom.cols());
            pts.row(0) = pts_pom;

            evalShearStressAtWall_singlePatch(p, pts, uSolField, shearStressVals);

            // Evaluate geometry at pts
            geoVals = patch->eval(pts);

            geoVals = geoVals.transpose()/0.0127;
            shearStressVals = shearStressVals.transpose()*m_blockAssembler.getViscosity()/(0.5 * referenceSolUVal(0,0) * referenceSolUVal(0,0));

            cpVals = pSolField.value(pts, p);
            cpVals = cpVals.transpose();
            for(int i = 0; i < cpVals.rows(); i++)
                for(int j = 0; j < cpVals.cols(); j++)
                    cpVals(i, j) = ((cpVals(i, j) - refSolPVal) / (0.5 * referenceSolUVal(0,0) * referenceSolUVal(0,0))) - refSolPValTrans;


            gsMatrix<T> result;
            result.setZero(geoVals.rows(), 2);
            result.col(0) = geoVals.col(0);
            //result.col(1) = geoVals.col(1);
            result.col(1) = shearStressVals;

            file << "x\tCf\n";
            for (int i = 0; i < result.rows(); i++)
                file << result(i,0) << "\t" << result(i, 1) << "\n";

            fileCp << "x\tCp\n";
            for (int i = 0; i < result.rows(); i++)
                fileCp << result(i,0) << "\t" << cpVals(i, 0) << "\n";

        }
            //-----------------------------------------------------------------
        for(int i = 0; i < 1; i++)
        {
            int p = 0; //vertical

                gsMatrix<T> geoVals, shearStressVals, cpVals;
                gsVector<unsigned> np;

                gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
                //const int parDim = patch->domainDim();
                //const int tarDim = patch->targetDim();

                boxSide side = boundary::south;

                typename gsGeometry<T>::uPtr geomBoundary = m_blockAssembler.getPatches().patch(p).boundary(side);

                gsMatrix<T> ab = geomBoundary->support();
                gsVector<T> a = ab.col(0);
                gsVector<T> b = ab.col(1);
                np = uniformSampleCount(a, b, npts);
                gsMatrix<T> pts_pom = gsPointGrid(a, b, np);
                gsMatrix<T> pts;
                pts.setZero(2,pts_pom.cols());
                pts.row(1) = pts_pom;

                evalShearStressAtWall_singlePatch(p, pts, uSolField, shearStressVals, true);

                // Evaluate geometry at pts
                geoVals = patch->eval(pts);

                geoVals = geoVals.transpose()/0.0127;
                shearStressVals = shearStressVals.transpose()*m_blockAssembler.getViscosity()/(0.5 * referenceSolUVal(0,0) * referenceSolUVal(0,0));

                cpVals = pSolField.value(pts, p);
                cpVals = cpVals.transpose();
                for(int i = 0; i < cpVals.rows(); i++)
                    for(int j = 0; j < cpVals.cols(); j++)
                        cpVals(i, j) = ((cpVals(i, j) - refSolPVal) / (0.5 * referenceSolUVal(0,0) * referenceSolUVal(0,0))) - refSolPValTrans;

                gsMatrix<T> result;
                result.setZero(geoVals.rows(), 2);
                result.col(0) = geoVals.col(0);
                //result.col(1) = geoVals.col(1);
                result.col(1) = shearStressVals;

                //file << "x\tCf\n";
                for (int i = result.rows()-1; i >= 0; i--)
                    file << result(i,0) << "\t" << result(i, 1) << "\n";


                for (int i = result.rows()-1; i >= 0; i--)
                    fileCp << result(i,0) << "\t" << cpVals(i, 0) << "\n";

        }
            //-----------------------------------------------------------------
        for(int i = 0; i < 1; i++)
        {
                int p = 0;

                    gsMatrix<T> geoVals, shearStressVals, cpVals;
                    gsVector<unsigned> np;

                    gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
                    //const int parDim = patch->domainDim();
                    //const int tarDim = patch->targetDim();

                    boxSide side = boundary::south;

                    typename gsGeometry<T>::uPtr geomBoundary = m_blockAssembler.getPatches().patch(p).boundary(side);

                    gsMatrix<T> ab = geomBoundary->support();
                    gsVector<T> a = ab.col(0);
                    gsVector<T> b = ab.col(1);
                    np = uniformSampleCount(a, b, npts);
                    gsMatrix<T> pts_pom = gsPointGrid(a, b, np);
                    gsMatrix<T> pts;
                    pts.setZero(2,pts_pom.cols());
                    pts.row(0) = pts_pom;

                    evalShearStressAtWall_singlePatch(p, pts, uSolField, shearStressVals);

                    // Evaluate geometry at pts
                    geoVals = patch->eval(pts);

                    geoVals = geoVals.transpose()/0.0127;
                    shearStressVals = shearStressVals.transpose()*m_blockAssembler.getViscosity()/(0.5 * referenceSolUVal(0,0) * referenceSolUVal(0,0));

                    cpVals = pSolField.value(pts, p);
                    cpVals = cpVals.transpose();
                    /*gsInfo << "cpVals = " << cpVals << "\n";
                    gsInfo << "cpVals.rows() = " << cpVals.rows() << "\n";
                    gsInfo << "cpVals.cols() = " << cpVals.cols() << "\n";*/
                    for(int i = 0; i < cpVals.rows(); i++)
                        for(int j = 0; j < cpVals.cols(); j++)
                            cpVals(i, j) = ((cpVals(i, j) - refSolPVal) / (0.5 * referenceSolUVal(0,0) * referenceSolUVal(0,0))) - refSolPValTrans;
/*                    gsInfo << "pts = " << pts << "\n";
                    gsInfo << "cpVals = " << cpVals << "\n";
                    gsInfo << "cpVals.rows() = " << cpVals.rows() << "\n";
                    gsInfo << "cpVals.cols() = " << cpVals.cols() << "\n";
*/


                    gsMatrix<T> result;
                    result.setZero(geoVals.rows(), 2);
                    result.col(0) = geoVals.col(0);
                    //result.col(1) = geoVals.col(1);
                    result.col(1) = shearStressVals;

                    //file << "x\tCf\n";
                    for (int i = 0; i < result.rows(); i++)
                        file << result(i,0) << "\t" << result(i, 1) << "\n";

                    for (int i = 0; i < result.rows(); i++)
                        fileCp << result(i,0) << "\t" << cpVals(i, 0) << "\n";

        }


            file.close();
            fileCp.close();

            //fd << result;
            //fd.save(fileNameM + ".xml");
            gsInfo << "Cf csv saved\n";
    }

protected:
    void evalShearStressAtWall_singlePatch(index_t patchIndex, gsMatrix<T>& pts, gsField<T>& uSolField, gsMatrix<T>& shearStressVals, bool vertical = false)
    {
        gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(patchIndex);
        gsBasisRefs<T> bases(m_blockAssembler.getBases(), patchIndex);

        const int tarDim = patch->targetDim();

        unsigned evFlags = NEED_VALUE | NEED_GRAD_TRANSFORM;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, *patch));

        gsMatrix<index_t> activesU;

        gsMatrix<T> basisGradsU, physGradU;
        bases.front().deriv_into(pts, basisGradsU);

        geoEval->evaluateAt(pts);

        index_t nQuPoints = pts.cols();

        gsMatrix<T> actCoeffsU;
        gsMatrix<T> uGrads;
        gsMatrix<T> tauWvals;
        tauWvals.setZero(nQuPoints);
        // Evaluate shear stress at pts
        for (int k = 0; k < nQuPoints; k++)
        {
            // eval uGrads at all pts
            bases.front().active_into(pts.col(k), activesU);
            int numActU = activesU.rows();
            actCoeffsU.setZero(tarDim, numActU);

            for (int j = 0; j < numActU; j++)
                actCoeffsU.col(j) = uSolField.coefficientVector(patchIndex).row(activesU(j)).transpose();

            geoEval->transformGradients(k, basisGradsU, physGradU);
            uGrads.noalias() = actCoeffsU * physGradU.transpose();

            if (vertical)
                tauWvals(k) = uGrads(1, 0);
            else
                tauWvals(k) = uGrads(0, 1);//;uGrads.row(var);
        }

        shearStressVals = tauWvals.transpose();
    }

    /*public:
        void plotSkinFrictionCoeff(std::string const & fn, gsMatrix<T> & solution, unsigned npts)
        {
            gsParaviewCollection collectionM(fn);
            std::string fileNameM;

            gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);

            for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
            {
                fileNameM = fn + util::to_string(p);

                gsMatrix<T> geoVals, shearStressVals;
                gsVector<unsigned> np;

                gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
                const int parDim = patch->domainDim();
                const int tarDim = patch->targetDim();

                gsMatrix<T> ab = patch->support();
                gsVector<T> a = ab.col(0);
                gsVector<T> b = ab.col(1);
                np = uniformSampleCount(a, b, npts);
                gsMatrix<T> pts = gsPointGrid(a, b, np);

                evalShearStress_singlePatch(p, pts, uSolField, shearStressVals);

                // Evaluate geometry at pts
                geoVals = patch->eval(pts);

                if (3 - parDim > 0)
                {
                    np.conservativeResize(3);
                    np.bottomRows(3 - parDim).setOnes();
                }
                if (3 - tarDim > 0)
                {
                    geoVals.conservativeResize(3, geoVals.cols());
                    geoVals.bottomRows(3 - tarDim).setZero();
                }

                gsWriteParaviewTPgrid(geoVals, shearStressVals, np.template cast<index_t>(), fileNameM);

                collectionM.addPart(fileNameM + ".vts");
            }
            collectionM.save();
        }

    protected:
        void evalSkinFrictionCoeff_singlePatch(index_t patchIndex, gsMatrix<T>& pts, gsField<T>& uSolField, gsMatrix<T>& shearStressVals)
        {
            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(patchIndex);
            gsBasisRefs<T> bases(m_blockAssembler.getBases(), patchIndex);

            const T viscosity = m_blockAssembler.getViscosity();

            const int tarDim = patch->targetDim();

            unsigned evFlags = NEED_VALUE | NEED_GRAD_TRANSFORM;
            typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, *patch));

            gsMatrix<index_t> activesU;

            gsMatrix<T> basisGradsU, physGradU;
            bases.front().deriv_into(pts, basisGradsU);

            geoEval->evaluateAt(pts);

            index_t nQuPoints = pts.cols();

            gsMatrix<T> actCoeffsU;
            gsMatrix<T> uGrads;
            gsMatrix<T> tauWvals;
            tauWvals.setZero(nQuPoints);
            // Evaluate shear stress at pts
            for (int k = 0; k < nQuPoints; k++)
            {
                // eval uGrads at all pts
                bases.front().active_into(pts.col(k), activesU);
                int numActU = activesU.rows();
                actCoeffsU.setZero(tarDim, numActU);

                for (int j = 0; j < numActU; j++)
                    actCoeffsU.col(j) = uSolField.coefficientVector(patchIndex).row(activesU(j)).transpose();

                geoEval->transformGradients(k, basisGradsU, physGradU);
                uGrads.noalias() = actCoeffsU * physGradU.transpose();

                tauWvals(k) = viscosity * uGrads(0, 1) / (0.5 * );//;uGrads.row(var);
            }

            shearStressVals = tauWvals.transpose();
        }*/

    void applyDirichletConditions()
    {
        //gsInfo << "Apply Dirichlet for pressure for cavity problem\n";
        for (index_t patchI = 0; patchI < 1; patchI++)
        {
            index_t patchIndex = patchI; //1st,2nd,3th patch
            T dirValue = 0.0;

            int uDofs = getUdofs();

            const gsMultiBasis<T>& basis = m_blockAssembler.getBases().at(1);

            boxSide bSide;
            bSide = boundary::west;
            gsMatrix<index_t> boundary = basis.piece(patchIndex).boundary(bSide);

            gsSparseMatrix<T> mat(m_matrix.rows(), m_matrix.cols());
            mat.setIdentity();

            index_t k = 0;
            if (m_blockAssembler.getMappers().back().is_free(boundary.at(k), patchIndex)) // DoF value is in the solVector
            {
                index_t ii = m_blockAssembler.getMappers().back().index(boundary.at(k), patchIndex);
                mat.coeffRef(ii + 2 * uDofs, ii + 2 * uDofs) = 0.;
            }

            mat.prune(0,0);

            m_matrix = mat * m_matrix;

            if (m_blockAssembler.getMappers().back().is_free(boundary.at(k), patchIndex)) // DoF value is in the solVector
            {
                index_t ii = m_blockAssembler.getMappers().back().index(boundary.at(k), patchIndex);

                m_matrix.coeffRef(ii + 2 * uDofs, ii + 2 * uDofs) = 1.0;
                m_rhs(ii + 2 * uDofs, 0) = dirValue;
            }
        }
    }

public:
    void plotGradient(std::string const & fn, gsMatrix<T> & solution, unsigned npts)
    {
        std::string fnX = fn + "velocityUx";
        std::string fnY = fn + "velocityUy";
        std::string fnP = fn + "pressure";

        gsParaviewCollection collectionX(fnX);
        gsParaviewCollection collectionY(fnY);
        gsParaviewCollection collectionC(fnP);
        std::string fileNameX, fileNameY, fileNameP;

        gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);
        gsField<T> pSolField = m_blockAssembler.constructSolution(solution, 1);

        for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
        {
            fileNameX = fnX + util::to_string(p);
            fileNameY = fnY + util::to_string(p);
            fileNameP = fnP + util::to_string(p);

            gsMatrix<T> geoVals, gradientUxVals, gradientUyVals, gradientPVals;
            gsVector<unsigned> np;

            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
            const int parDim = patch->domainDim();
            const int tarDim = patch->targetDim();

            gsMatrix<T> ab = patch->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);
            np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts = gsPointGrid(a, b, np);

            evalGradient_singlePatch(p, pts, uSolField, pSolField, gradientUxVals, gradientUyVals, gradientPVals);

            // Evaluate geometry at pts
            geoVals = patch->eval(pts);

            if (3 - parDim > 0)
            {
                np.conservativeResize(3);
                np.bottomRows(3 - parDim).setOnes();
            }
            if (3 - tarDim > 0)
            {
                geoVals.conservativeResize(3, geoVals.cols());
                geoVals.bottomRows(3 - tarDim).setZero();
            }

            gsWriteParaviewTPgrid(geoVals, gradientUxVals, np.template cast<index_t>(), fileNameX);
            gsWriteParaviewTPgrid(geoVals, gradientUyVals, np.template cast<index_t>(), fileNameY);
            gsWriteParaviewTPgrid(geoVals, gradientPVals, np.template cast<index_t>(), fileNameP);

            collectionX.addPart(fileNameX + ".vts");
            collectionY.addPart(fileNameY + ".vts");
            collectionC.addPart(fileNameP + ".vts");
        }
        collectionX.save();
        collectionY.save();
        collectionC.save();
    }

    void evalGradient(gsMatrix<T> & solution, gsMatrix<T>& gradUxVals, gsMatrix<T>& gradUyVals, gsMatrix<T>& gradPVals, unsigned npts)
    {
        gsField<T> uSolField = m_blockAssembler.constructSolution(solution, 0);
        gsField<T> pSolField = m_blockAssembler.constructSolution(solution, 1);

        for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
        {
            gsMatrix<T> geoVals, gradientUxVals, gradientUyVals, gradientPVals;
            gsVector<unsigned> np;

            gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(p);
            const int parDim = patch->domainDim();
            const int tarDim = patch->targetDim();

            gsMatrix<T> ab = patch->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);
            np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts = gsPointGrid(a, b, np);

            evalGradient_singlePatch(p, pts, uSolField, pSolField, gradientUxVals, gradientUyVals, gradientPVals);

            // Evaluate geometry at pts
            geoVals = patch->eval(pts);

            if (3 - parDim > 0)
            {
                np.conservativeResize(3);
                np.bottomRows(3 - parDim).setOnes();
            }
            if (3 - tarDim > 0)
            {
                geoVals.conservativeResize(3, geoVals.cols());
                geoVals.bottomRows(3 - tarDim).setZero();
            }

            gradUxVals = gradientUxVals;
            gradUyVals = gradientUyVals;
            gradPVals = gradientPVals;
        }
    }

public:
    void evalGradient_singlePatch(index_t patchIndex, gsMatrix<T>& pts, gsField<T>& uSolField, gsField<T>& pSolField,
                                             gsMatrix<T>& gradUxVals, gsMatrix<T>& gradUyVals, gsMatrix<T>& gradPVals)
    {
        gsGeometry<T>* patch = &m_blockAssembler.getPatches().patch(patchIndex);
        gsBasisRefs<T> bases(m_blockAssembler.getBases(), patchIndex);

        const int tarDim = patch->targetDim();

        unsigned evFlags = NEED_GRAD_TRANSFORM;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, *patch));

        gsMatrix<index_t> activesU, activesP;

        gsMatrix<T> basisGradsU, physGradU;
        bases.front().deriv_into(pts, basisGradsU);
        gsMatrix<T> bHessianU;
        bases.front().deriv2_into(pts, bHessianU);

        gsMatrix<T> basisGradsP, physGradP;
        bases.back().deriv_into(pts, basisGradsP);

        geoEval->evaluateAt(pts);

        index_t nQuPoints = pts.cols();

        gsMatrix<T> actCoeffsU, actCoeffsP;
        gsMatrix<T> uGrads, pGrads;
        gsMatrix<T> gXVals;
        gXVals.setZero(nQuPoints, tarDim);
        gsMatrix<T> gYVals;
        gYVals.setZero(nQuPoints, tarDim);
        gsMatrix<T> gPvals;
        gPvals.setZero(nQuPoints, tarDim);
        // Evaluate gradient at pts
        for (int k = 0; k < nQuPoints; k++)
        {
            // eval uGrads at all pts
            bases.front().active_into(pts.col(k), activesU);
            int numActU = activesU.rows();
            actCoeffsU.setZero(tarDim, numActU);
            bases.back().active_into(pts.col(k), activesP);
            int numActP = activesP.rows();
            actCoeffsP.setZero(1, numActP);

            for (int j = 0; j < numActU; j++)
                actCoeffsU.col(j) = uSolField.coefficientVector(patchIndex).row(activesU(j)).transpose();
            for (int j = 0; j < numActP; j++)
                actCoeffsP.col(j) = pSolField.coefficientVector(patchIndex).row(activesP(j)).transpose();

            geoEval->transformGradients(k, basisGradsU, physGradU);
            uGrads.noalias() = actCoeffsU * physGradU.transpose();
            geoEval->transformGradients(k, basisGradsP, physGradP);
            pGrads.noalias() = actCoeffsP * physGradP.transpose();

            gXVals.row(k) = uGrads.row(0);
            gYVals.row(k) = uGrads.row(1);
            gPvals.row(k) = pGrads.row(0);
        }
        gradUxVals = gXVals.transpose();
        gradUyVals = gYVals.transpose();
        gradPVals = gPvals.transpose();
    }

public:
    //computes integral (p*n_p*n_s) over patch side, where p is pressure, n_p is unit outer normal of the patch side, n_s is the normal vector to the stream
    T computeLift(gsField<T> pSolField)
    {
        gsInfo << "\nCompute lift coefficient...\n";
        T rho = 1.;
        T dynamicP = 0.5 * rho * math::pow(m_freeStreamVel, 2);

        int tarDim = m_blockAssembler.getPatches().dim();
        T lift = 0.;
        for (int index = 0; index < m_liftPatches.rows(); index++)
        {
            const gsBasis<T> & basis = m_blockAssembler.getBases().at(1).basis(m_liftPatches[index]);

            gsVector<int> numQuadNodes(tarDim);
            const int dir = m_liftSides[index].direction();
            for (int i = 0; i < tarDim; ++i)
                numQuadNodes[i] = (2 * basis.degree(i) + 1);
            numQuadNodes[dir] = 1;

            // Setup Quadrature
            gsGaussRule<T> QuRule(numQuadNodes);// harmless slicing occurs here

            gsMatrix<T> quNodes; // Mapped nodes
            gsVector<T> quWeights; // Mapped weights
            unsigned evFlags = NEED_VALUE | NEED_OUTER_NORMAL;   // Evaluation flags for the Geometry map

            // Initialize geometry evaluator
             typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_blockAssembler.getPatches().patch(m_liftPatches[index])));

            typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(m_liftSides[index]);
            for (; domIt->good(); domIt->next())
            {
                // Compute the quadrature rule
                QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

                // Compute image of Gauss nodes under geometry mapping as well as Jacobians
                geoEval->evaluateAt(quNodes);

                // Evaluate solution on element nodes
                gsMatrix<T> solPVals = pSolField.value(quNodes, m_liftPatches[index]);

                for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
                {
                    // Compute the outer normal vector
                    gsVector<T> normal;
                    geoEval->outerNormal(k, m_liftSides[index], normal);
                    //lift force
                    lift += quWeights[k] * solPVals(k) * normal.dot(m_freeStream.transpose());//p*cos(eta)

                }
            }
        }
        //lift coefficient
        return (lift / (dynamicP * m_liftReferenceArea));
        gsInfo << "Lift coefficient computed.\n";
    }

    void setLiftParams(gsVector<int> patchNumbers, std::vector<boxSide> sides, T freeStreamVel, T referenceArea, gsVector<T> freeStream)
    {
        m_liftSides = sides;
        m_liftPatches = patchNumbers;
        m_freeStreamVel = freeStreamVel;
        m_liftReferenceArea = referenceArea;
        m_freeStream = freeStream;
    }

public:
    void plotResiduum(std::string const & fn, gsMatrix<T> & solution, gsMatrix<T>& solutionOld, gsMatrix<T>& TMsolution, unsigned npts)
    { GISMO_NO_IMPLEMENTATION }
    void plotResiduum(std::string const & fn, gsMatrix<T> & solution, gsMatrix<T>& solutionOld, unsigned npts)
    { GISMO_NO_IMPLEMENTATION }

    virtual gsMatrix<T> getSolution_full(const gsMatrix<T>& solVector) const
    { GISMO_NO_IMPLEMENTATION }

public:
    bool isInitialized() { return m_bInitialized; }
    const uwbINSBlockAssembler<T>& getBlockAssembler() const { return m_blockAssembler; }
    uwbINSBlockAssembler<T>& getBlockAssembler() { return m_blockAssembler; }
    virtual const gsMatrix<T>& getSolution() const { return m_solution; }
    virtual int numDofs() const { return m_blockAssembler.numDofs(); }
    virtual int getUdofs() const { return m_blockAssembler.getUdofs(); }
    virtual int getPdofs() const { return m_blockAssembler.getPdofs(); }
    virtual int getPshift() const { return m_blockAssembler.getPshift(); }
    int getTarDim() const { return m_blockAssembler.getTarDim(); }
    real_t getViscosity() const { return m_blockAssembler.getViscosity(); }
    bool isRotation() const { return m_blockAssembler.isRotation(); }

protected:
    uwbINSBlockAssembler<T>     m_blockAssembler;
    gsMatrix<T>                 m_solution;
    bool                        m_bInitialized;

    bool                        m_bCavity;

    gsSparseMatrix<T> m_matrix;
    gsMatrix<T> m_rhs;

    gsVector<int> m_liftPatches;
    std::vector<boxSide> m_liftSides;
    T m_freeStreamVel;
    T m_liftReferenceArea;
    gsVector<T> m_freeStream;

}; //uwbINSAssemblerBase

} //namespace gismo
