
#pragma once

#include <iomanip>
#include <ios>

#include <gsIO/gsLsdynaUtils.h>

namespace gismo
{

namespace lsdyna    
{

class gsLsdynaEigProblem
{
public:

    // --------------------------------------------------
    // constructors - destructors
    // --------------------------------------------------
    
    gsLsdynaEigProblem()
        :
        numEigValsExtr(200),
        flagConsisMassMat(1),
        flagImplExplAnalysis(1),
        dt0(1e-3),
        solverMethod(4),
        flagPrintSolver(2),
        terminationTime(5e-3),
        timeIntervalBinary(1e-3),
        timeIntervalAscii(1e-3)
    {
        
    }


    ~gsLsdynaEigProblem()
    {

    }

public:
    
    void write(std::ostream& out) const
    {
        setGlobalStreamState(out);
        writeProblemShell(out);
    }

    
    // --------------------------------------------------
    // setters
    // --------------------------------------------------
    
    void setNumberOfEigenvalsToExtract(const lsInt number)
    {
        numEigValsExtr = number;
    }

    void setConsistentMassMatFlag(const lsInt flag)
    {
        flagConsisMassMat = flag;
    }

    void setImplicitExplicitAnalysisFlag(const lsInt flag)
    {
        flagImplExplAnalysis = flag;
    }

    void setTimeStepForImplicitAnalysis(const lsFloat timeStep)
    {
        dt0 = timeStep;
    }

    void setSolverMethod(const lsInt method)
    {
        solverMethod = method;
    }

    void setPrintSolverFlag(const lsInt flag)
    {
        flagPrintSolver = flag;
    }

    void setTerminationTime(const lsFloat time)
    {
        terminationTime = time;
    }

    void setTimeIntervalForBinaryOutput(const lsFloat time)
    {
        timeIntervalBinary = time;
    }

    void setTimeIntervalForASCIIOutput(const lsFloat time)
    {
        timeIntervalAscii = time;
    }
    
private:

    void setGlobalStreamState(std::ostream& out) const 
    {
        out << std::uppercase; 
    }
    
    void writeProblemShell(std::ostream& out) const
    {
        write::CONTROL_IMPLICIT_EIGENVALUE(out, numEigValsExtr);

        write::CONTROL_IMPLICIT_CONSISTENT_MASS(out, flagConsisMassMat);

        write::CONTROL_IMPLICIT_GENERAL(out, flagImplExplAnalysis, dt0);

        write::CONTROL_IMPLICIT_SOLVER(out, solverMethod, flagPrintSolver);

        write::CONTROL_TERMINATION(out, terminationTime);

        write::DATABASE_BINARY_D3PLOT(out, timeIntervalBinary);

        write::DATABASE_GLSTAT(out, timeIntervalAscii);
    }

    // disable copy constructor and assignment operator
    gsLsdynaEigProblem(const gsLsdynaEigProblem& other);
    const gsLsdynaEigProblem operator=(const gsLsdynaEigProblem& other);

    // --------------------------------------------------
    // data memebers
    // --------------------------------------------------
    
private:
    // number of eigenvalues to extract
    lsInt numEigValsExtr;

    // consistent mass matrix flag
    lsInt flagConsisMassMat;

    // Implicit / Explicit analysis type flag
    lsInt flagImplExplAnalysis;

    // Initial time step size for implicit analysis
    lsFloat dt0;

    // Linear equation solver method
    lsInt solverMethod;

    // Linear solver print flag controls screen and message file output
    lsInt flagPrintSolver;

    // Termination time
    lsFloat terminationTime;

    // Time interval for binary output.
    lsFloat timeIntervalBinary;

    // Time interval for ascii files.
    lsFloat timeIntervalAscii;
};


} // namespace lsdyna

} // namespace gismo
