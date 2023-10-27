/** @file gsTimeParallelMultigrid.h

    @brief A time parallel multigrid class for space time formulations.
    Originally coded by Martin Neumueller and adapted to G+Smo by
    Christoph Hofer.


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer, M. Neumüller
    Created on: 2017-07-04
*/


#pragma once
#include<gsCore/gsConfig.h>
#include <ctime>

namespace gismo {

#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>
#include <gsMpi/gsMpiComm.h>

using std::vector;



//--------------------------------------------------------------------------------
//Funktionalität von STSlapOperator

//   void solveOnSTSla(Vector &u, const Vector &f, int timeLevel, int timeStep)
//       löst für den gegeben Zeitschritt timeStep und einem gegebenen Zeitlevel
//       das Gleichungsystem mit der rechtens f und speichert es in den Vector u.
//	     Dabei haben die Vectoren schon die richtige größe

//   void calcInitialVector(const Vector &u, Vector &f, int timeLevel, int timeStep)
//       berechnet für einen gegebenen Zeitschritt "timeStep" und ein gegebenes
//       Zeitlevel "timeLevel" die Rechte Seite für den nächsten Zeitschritt
//       "timeStep+1". Dabei ist der Vector u der Lösungsvektor vom Zeitschritt
//       "timeStep". Weiters wird auf den Vector f nur draufaddiert, er wird also
//       nicht überschrieben.

//   int getNDofs(void)
//       gibt die Anzahl der Freiheitsgrade für ein Zeitintervall zurück

//--------------------------------------------------------------------------------

template <class STSlapOperator> class STMultigridMPI
{
public:
    typedef  gsMatrix<real_t> Vector;



    struct MPI_Comm_Vars
    {
        int my_rank;
        int size;
        MPI_Comm Comm;

    };

    struct STSolVector
    {
        vector<vector< Vector> > u; // vector for solution
        vector<vector< Vector> > v; // vector for temporary computations
        vector<vector< Vector> > w; // vector for defects
        vector<vector< Vector> > x; // vector for the approximation on each space-time slap (for the smoother)
    };

    enum CoarseningStrategy
    {
        NO_SPACE_COARSENING=0,
        LOW_SPACE_COARSENING,
        MODERATE_SPACE_COARSENING,
        AGGRESSIVE_SPACE_COARSENING,
        AUTOMATIC_SPACE_COARSENING,
        STOCHASTIC_SPACE_COARSENING
    };


private:
    MPI_Comm_Vars MPICommVars_;

    int NTimeLevels_;
    int NSpaceLevels_;

    int timeSm1_,timeSm2_;
    int timeCycles_;
    double timeOmega_;
    int coarseLevel_;
    bool useHybridSmoother_;

    vector<int> spaceLevels_;
    vector<int> timeSteps_;

    vector<typename STSlapOperator::Ptr> A_;
    vector<int> spaceSizes_;

    vector<bool> activeTimeLevels_;
    vector<int> commForward_;
    vector<int> commBackward_;
    vector<int> MytimeSteps_;
    vector<vector<int> > MytimeStepsIdx_;

    MPI_Request send_request_;
    MPI_Status send_status_;
    MPI_Status receive_status_;
    MPI_Request receive_request_;
    bool bufferInitialized_;
    vector<double*> sendBuffer_;
    vector<double*> receiveBuffer_;
    vector<Vector> xold_;

    bool useTriDiag_;
    vector<Vector> xnew_;
    MPI_Request send_requestFwd_;
    MPI_Status send_statusFwd_;
    MPI_Status receive_statusFwd_;
    MPI_Request receive_requestFwd_;
    vector<double*> sendBufferFwd_;
    vector<double*> receiveBufferFwd_;

    void initActiveTimeLevels(void)
    {
        activeTimeLevels_.resize(NTimeLevels_);
        commForward_.resize(NTimeLevels_);
        commBackward_.resize(NTimeLevels_);
        MytimeSteps_.resize(NTimeLevels_);
        MytimeStepsIdx_.resize(NTimeLevels_);

        for(unsigned int i=0; i<activeTimeLevels_.size(); i++)
        {
            int timeSteps = pow(2,i);
            if(timeSteps>=MPICommVars_.size)
            {
                activeTimeLevels_[i] = true;
                commBackward_[i] = MPICommVars_.my_rank-1;
                if(MPICommVars_.my_rank==MPICommVars_.size-1) commForward_[i] = -1;
                else commForward_[i] = MPICommVars_.my_rank+1;
                int factor = timeSteps/MPICommVars_.size;
                MytimeSteps_[i] = factor;
                for(int j=0; j<MytimeSteps_[i];++j)
                    MytimeStepsIdx_[i].push_back(MPICommVars_.my_rank*MytimeSteps_[i]+j);
            }
            else
            {
                int factor = MPICommVars_.size/timeSteps;
                if(MPICommVars_.my_rank%factor==0)
                {
                    activeTimeLevels_[i] = true;
                    if(MPICommVars_.my_rank==MPICommVars_.size-factor) commForward_[i] = -1;
                    else commForward_[i] = MPICommVars_.my_rank+factor;
                    if(MPICommVars_.my_rank!=0) commBackward_[i] = MPICommVars_.my_rank-factor;
                    else commBackward_[i] = -1;
                    MytimeSteps_[i] = 1;

                    MytimeStepsIdx_[i].push_back(MPICommVars_.my_rank/factor*MytimeSteps_[i]);
                }
                else
                {
                    activeTimeLevels_[i] = false;
                    commForward_[i] = -1;
                    commBackward_[i] = -1;
                    MytimeSteps_[i] = 0;
                }
            }
        }
    }

public:
    STMultigridMPI(MPI_Comm_Vars MPICommVars, int timeLevels, int spaceLevels, gsOptionList opt)
    {
        MPICommVars_ = MPICommVars;
        NTimeLevels_ = timeLevels;
        NSpaceLevels_ = spaceLevels;

        gsOptionList def = defaultOptions();
        def.update(opt);
        timeSm1_ = def.getInt("NumPreSmooth");
        timeSm2_ = def.getInt("NumPostSmooth");
        timeCycles_ = def.getInt("NumCycles");
        timeOmega_ = def.getReal("Damping");
        coarseLevel_ = def.getInt("CoarseLevel");
        useHybridSmoother_ = def.getSwitch("Hybrid");
        useTriDiag_ = def.getSwitch("TriDiag");

        //init spaceLevels - the standard is no space coarsening
        spaceLevels_.resize(NTimeLevels_); for(int i=0; i<NTimeLevels_; i++) spaceLevels_[i] = NSpaceLevels_-1;

        //compute and store the number of timesteps in each level
        timeSteps_.resize(NTimeLevels_); for(int i=0; i<NTimeLevels_; i++) timeSteps_[i] = pow(2,i);

        A_.resize(NSpaceLevels_);
        spaceSizes_.resize(NSpaceLevels_); for(int i=0; i<NSpaceLevels_; i++) spaceSizes_[i] = -1;

        bufferInitialized_ = false;

        initActiveTimeLevels();
    }

    gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addReal("Tolerance", "tolerance for the space time MG", 1.e-8);
        opt.addInt("MaxIters", "maximal number of Iteration", 100);
        opt.addInt("NumPreSmooth", "number of smoothing steps", 2);
        opt.addInt("NumPostSmooth", "number of smoothing steps", 2);
        opt.addReal("Damping", "damping parameter omega", 0.5);
        opt.addInt("NumCycles", "number of cycles", 1);
        opt.addInt("CoarseLevel", "the coarsest level",0);
        opt.addSwitch("TriDiag", "handle TriDiag formulation", false);
        opt.addSwitch("Hybrid", "use hybrid smoother", false);

        return opt;
    }

    ~STMultigridMPI()
    {
        if(bufferInitialized_)
        {
            for(unsigned int i=0; i<receiveBuffer_.size(); i++)
            {
                delete [] receiveBuffer_[i];
                delete [] sendBuffer_[i];

                if(useTriDiag_)
                {
                    delete [] receiveBufferFwd_[i];
                    delete [] sendBufferFwd_[i];
                }
            }
        }
    }

    bool IamLastProcOnLevel(int level)
    {
        return ((activeTimeLevels_[level] && ((commForward_[level]<0 && commBackward_[level]>=0) || level==0)) || MPICommVars_.size==1);
    }

    void setCoarseLevel(int level) { coarseLevel_ = level; }
    void setSmoothingSteps(int sm1, int sm2) { timeSm1_ = sm1; timeSm2_ = sm2; }
    void setDampingParameter(double omega) { timeOmega_=omega; }
    void setCycleIndex(int cycles) { timeCycles_ = cycles; }
    void setHybridSmoother(bool setVar) { useHybridSmoother_ = setVar; }

    int getMyTimeSteps(int level) const { return MytimeSteps_[level]; }
    const std::vector<int> &  getMyTimeStepsIdx(int level) const {return MytimeStepsIdx_[level];}

    const vector<int> & getSpaceLevels(void) const { return spaceLevels_; }

    void setSTSlapOperator(typename STSlapOperator::Ptr A, int spaceLevel)
    {
        A_[spaceLevel] = A;
        spaceSizes_[spaceLevel] = A_[spaceLevel]->getNDofs();
    }

    typename STSlapOperator::Ptr getSTSlapOperator(int spaceLevel) const {return A_[spaceLevel]; }

    void setUpSpaceTimeHierarchy(int strategy)
    {
        if(strategy==NO_SPACE_COARSENING)
        {
            for(int i=0; i<NTimeLevels_; i++) spaceLevels_[i] = NSpaceLevels_-1;
        }
        else if(strategy==AUTOMATIC_SPACE_COARSENING)
        {
            int actSpaceLevel = NSpaceLevels_-1;
            for(int i=NTimeLevels_-1; i>=0; i--)
            {
                spaceLevels_[i] = actSpaceLevel;

                if(A_[actSpaceLevel]->checkIfSpaceCoarseningIsAllowed(i))
                {
                    if(actSpaceLevel>0) actSpaceLevel--;
                }
            }
        }
        else if(strategy==STOCHASTIC_SPACE_COARSENING)
        {
            int actSpaceLevel = 0;
            spaceLevels_[0] = actSpaceLevel;
            if(actSpaceLevel<NSpaceLevels_-1) actSpaceLevel++;

            for(int i=1; i<NTimeLevels_; i+=2)
            {
                spaceLevels_[i] = actSpaceLevel;
                if(i<NTimeLevels_-1) spaceLevels_[i+1] = actSpaceLevel;
                if(actSpaceLevel<NSpaceLevels_-1) actSpaceLevel++;
            }
        }
        else if(NSpaceLevels_ >= NTimeLevels_ || strategy==AGGRESSIVE_SPACE_COARSENING)
        {
            int actSpaceLevel = NSpaceLevels_-1;
            for(int i=NTimeLevels_-1; i>=0; i--)
            {
                spaceLevels_[i] = actSpaceLevel;
                if(actSpaceLevel>0) actSpaceLevel--;
            }
        }
        else if(strategy==MODERATE_SPACE_COARSENING)
        {
            int timeLevelBuffer = int(double(NTimeLevels_-NSpaceLevels_)/2.0);

            int diff = 0;
            if(diff >= timeLevelBuffer+1) diff = timeLevelBuffer;

            int actSpaceLevel = 0;
            for(int i=0; i<NTimeLevels_; i++)
            {
                spaceLevels_[i] = actSpaceLevel;
                if(actSpaceLevel<NSpaceLevels_-1 && i>=timeLevelBuffer+diff) actSpaceLevel++;
            }
        }
        else if(strategy==LOW_SPACE_COARSENING)
        {
            int actSpaceLevel = 0;
            for(int i=0; i<NTimeLevels_; i++)
            {
                spaceLevels_[i] = actSpaceLevel;
                if(actSpaceLevel<NSpaceLevels_-1) actSpaceLevel++;
            }
        }


    }

    void printSpaceTimeHierarchy(void) const
    {
        gsInfo << "\t";
        for(int TLevel=NTimeLevels_-1; TLevel>=0; TLevel--)  { gsInfo << TLevel; if(TLevel>0) gsInfo << "\t"; } gsInfo << std::endl;

        for(int SLevel=NSpaceLevels_-1; SLevel>=0; SLevel--)
        {
            gsInfo << SLevel << "\t";
            for(int TLevel=NTimeLevels_-1; TLevel>=0; TLevel--)
            {
                if(SLevel==spaceLevels_[TLevel])
                {
                    gsInfo << "X";
                }
                if(TLevel>0) gsInfo << "\t";
            }
            gsInfo << std::endl;
        }
    }

    void initSTSolVectors(STSolVector &vec, bool useRandomVector=false)
    {
        vec.u.resize(NTimeLevels_);
        vec.v.resize(NTimeLevels_);
        vec.w.resize(NTimeLevels_);
        vec.x.resize(NTimeLevels_);

        if(bufferInitialized_) for(unsigned int i=0; i<receiveBuffer_.size(); i++) { delete [] receiveBuffer_[i]; delete [] sendBuffer_[i]; }

        receiveBuffer_.resize(NTimeLevels_);
        sendBuffer_.resize(NTimeLevels_);
        xold_.resize(NTimeLevels_);

        if(useTriDiag_)
        {
            xnew_.resize(NTimeLevels_);
            sendBufferFwd_.resize(NTimeLevels_);
            receiveBufferFwd_.resize(NTimeLevels_);
        }

        if(useRandomVector) srand(time(NULL));

        for(int level=0; level<NTimeLevels_; level++)
        {
            int spaceLevel = spaceLevels_[level];

            vec.u[level].resize(MytimeSteps_[level]);
            vec.v[level].resize(MytimeSteps_[level]);
            vec.w[level].resize(MytimeSteps_[level]);
            vec.x[level].resize(MytimeSteps_[level]);

            for(int j=0; j<MytimeSteps_[level]; j++)
            {
                vec.u[level][j].setZero(spaceSizes_[spaceLevel],1);
                if(useRandomVector)
                {
                    for(int i=0; i<spaceSizes_[spaceLevel]; i++) vec.u[level][j](i,0) = ((double) rand()/(RAND_MAX));
                }
                vec.v[level][j].setZero(spaceSizes_[spaceLevel],1);
                vec.w[level][j].setZero(spaceSizes_[spaceLevel],1);
                vec.x[level][j].setZero(spaceSizes_[spaceLevel],1);
            }

            if(level<NTimeLevels_-1) if(activeTimeLevels_[level]==false && activeTimeLevels_[level+1]==true)
            {
                vec.u[level].resize(1);
                vec.u[level][0].setZero(spaceSizes_[spaceLevel],1);
                if(useRandomVector)
                {
                    for(int i=0; i<spaceSizes_[spaceLevel]; i++) vec.u[level][0](i,0) = ((double) rand()/(RAND_MAX));
                }
                vec.v[level].resize(1); vec.v[level][0].setZero(spaceSizes_[spaceLevel],1);
                vec.w[level].resize(1); vec.w[level][0].setZero(spaceSizes_[spaceLevel],1);
                vec.x[level].resize(1); vec.x[level][0].setZero(spaceSizes_[spaceLevel],1);
            }

            receiveBuffer_[level] = new double[spaceSizes_[spaceLevel]];
            sendBuffer_[level] = new double[spaceSizes_[spaceLevel]];
            xold_[level].setZero(spaceSizes_[spaceLevel],1);

            if(useTriDiag_)
            {
                xnew_[level].setZero(spaceSizes_[spaceLevel],1);
                receiveBufferFwd_[level] = new double[spaceSizes_[spaceLevel]];
                sendBufferFwd_[level] = new double[spaceSizes_[spaceLevel]];
            }
        }

        bufferInitialized_ = true;
    }

    void zeroSTSolVectors(STSolVector &vec, bool useRandomVector=false)
    {
        if(useRandomVector) srand(time(NULL));

        for(int level=0; level<NTimeLevels_; level++)
        {
            int spaceLevel = spaceLevels_[level];

            for(int j=0; j<MytimeSteps_[level]; j++)
            {
                if(useRandomVector)
                {
                    for(int i=0; i<spaceSizes_[spaceLevel]; i++) vec.u[level][j](i,0) = ((double) rand()/(RAND_MAX));
                }
                else vec.u[level][j].setZero();
                vec.v[level][j].setZero();
                vec.w[level][j].setZero();
                vec.x[level][j].setZero();
            }

            if(level<NTimeLevels_-1) if(activeTimeLevels_[level]==false && activeTimeLevels_[level+1]==true)
            {
                if(useRandomVector)
                {
                    for(int i=0; i<spaceSizes_[spaceLevel]; i++) vec.u[level][0][i] = ((double) rand()/(RAND_MAX));
                }
                else vec.u[level][0].setZero();
                vec.v[level][0].setZero();
                vec.w[level][0].setZero();
                vec.x[level][0].setZero();
            }

            xold_[level].setZero();

            if(useTriDiag_)
            {
                xnew_[level].setZero();
            }
        }
    }

    void STForwardSolve(STSolVector &vec, const vector<Vector> &fl, int timeLevel) const
    {
        //coarse grid solver, solves simply forward in time on the level=timeLevel
        int spaceLevel = spaceLevels_[timeLevel];
        int nTimeSteps = MytimeSteps_[timeLevel];
        vec.v[timeLevel][0] = fl[0];
        if(commBackward_[timeLevel]>=0) A_[spaceLevel]->calcInitialVector(xold_[timeLevel], vec.v[timeLevel][0], timeLevel, 0);
        for(int ti=0; ti<nTimeSteps; ti++)
        {
            A_[spaceLevel]->solveOnSTSlap(vec.u[timeLevel][ti], vec.v[timeLevel][ti], timeLevel, ti);
            if(ti<nTimeSteps-1)
            {
                vec.v[timeLevel][ti+1] = fl[ti+1];
                A_[spaceLevel]->calcInitialVector(vec.u[timeLevel][ti], vec.v[timeLevel][ti+1], timeLevel, ti+1);
            }
        }
    }

    int STForwardSolveWithStatistics(STSolVector &vec, const vector<Vector> &fl, int timeLevel, bool output = true) const
    {
        //coarse grid solver, solves simply forward in time on the level=timeLevel
        int spaceLevel = spaceLevels_[timeLevel];
        int nTimeSteps = MytimeSteps_[timeLevel];
        vec.v[timeLevel][0] = fl[0];
        if(commBackward_[timeLevel]>=0) A_[spaceLevel]->calcInitialVector(xold_[timeLevel], vec.v[timeLevel][0], timeLevel, 0);

        int iterations = 0;
        int maxIt = 0; int minIt = 0;
        for(int ti=0; ti<nTimeSteps; ti++)
        {
            int iter = A_[spaceLevel]->solveOnSTSlap(vec.u[timeLevel][ti], vec.v[timeLevel][ti], timeLevel, ti);
            iterations += iter;
            if(maxIt<iter) maxIt = iter;
            if(minIt>iter || minIt==0) minIt = iter;
            if(output) gsInfo << "average it: " << double(iterations)/double(ti+1) << "          \r"<< std::flush;

            if(ti<nTimeSteps-1)
            {
                vec.v[timeLevel][ti+1] = fl[ti+1];
                A_[spaceLevel]->calcInitialVector(vec.u[timeLevel][ti], vec.v[timeLevel][ti+1], timeLevel, ti);
            }
        }

        if(output) gsInfo << "average iterations: " << double(iterations)/double(nTimeSteps) << " (min: " << minIt << ", max: " << maxIt << ")" << std::endl;

        return maxIt;
    }

    int STForwardSolveBackBackWithStatistics(STSolVector &vec, const vector<Vector> &fl, int timeLevel, bool output = true) const
    {
        //coarse grid solver, solves simply forward in time on the level=timeLevel
        int spaceLevel = spaceLevels_[timeLevel];
        int nTimeSteps = MytimeSteps_[timeLevel];
        vec.v[timeLevel][0] = fl[0];
        //			if(commBackward_[timeLevel]>=0) A_[spaceLevel]->calcInitialVector(xold_[timeLevel], vec.v[timeLevel][0], timeLevel, 0);

        int iterations = 0;
        int maxIt = 0; int minIt = 0;
        for(int ti=0; ti<nTimeSteps; ti++)
        {
            int iter = A_[spaceLevel]->solveOnSTSlap(vec.u[timeLevel][ti], vec.v[timeLevel][ti], timeLevel, ti);
            iterations += iter;
            if(maxIt<iter) maxIt = iter;
            if(minIt>iter || minIt==0) minIt = iter;
            if(output) gsInfo << "average it: " << double(iterations)/double(ti+1) << "          \r" << std::flush;

            if(ti<nTimeSteps-1)
            {
                vec.v[timeLevel][ti+1] = fl[ti+1];
                A_[spaceLevel]->calcInitialVector(vec.u[timeLevel][ti], vec.v[timeLevel][ti+1], timeLevel, ti);

                if(ti>0)
                {
                    A_[spaceLevel]->calcPreInitialVector(vec.u[timeLevel][ti-1], vec.v[timeLevel][ti+1], timeLevel, ti);
                }
            }
        }

        if(output) gsInfo << "average iterations: " << double(iterations)/double(nTimeSteps) << " (min: " << minIt << ", max: " << maxIt << ")" << std::endl;

        return maxIt;
    }

    void computeResidual(const vector<Vector> &ul, const Vector &fl, Vector &residual, int timeLevel, int ti)
    {
        int spaceLevel = spaceLevels_[timeLevel];

        residual = fl;
        double sign = -1.0; A_[spaceLevel]->mult(ul[ti], residual, timeLevel, ti, sign);
        if(ti!=0) A_[spaceLevel]->calcInitialVector(ul[ti-1], residual, timeLevel, ti);
        if(commBackward_[timeLevel]>=0 && ti==0) A_[spaceLevel]->calcInitialVector(xold_[timeLevel], residual, timeLevel, ti);
    }
    void computeResidualTriDiag(const vector<Vector> &ul, const Vector &fl, Vector &residual, int timeLevel, int ti)
    {
        int spaceLevel = spaceLevels_[timeLevel];

        residual = fl;
        double sign = -1.0; A_[spaceLevel]->mult(ul[ti], residual, timeLevel, ti, sign);
        if(ti!=0) A_[spaceLevel]->calcInitialVector(ul[ti-1], residual, timeLevel, ti);
        if(commBackward_[timeLevel]>=0 && ti==0) A_[spaceLevel]->calcInitialVector(xold_[timeLevel], residual, timeLevel, ti);

        if(ti=MytimeSteps_[timeLevel]-1) A_[spaceLevel]->calcTerminalVector(ul[ti+1], residual, timeLevel, ti);
        if(commForward_[timeLevel]>=0 && ti==MytimeSteps_[timeLevel]-1) A_[spaceLevel]->calcTerminalVector(xnew_[timeLevel], residual, timeLevel, ti);
    }

    void SpaceTimeSmoother(STSolVector &vec, const vector<Vector> &fl, int timeLevel, int sm=1)
    {
        int spaceLevel = spaceLevels_[timeLevel];

        //jacobi - smoothing
        for(int k=0; k<sm; k++)
        {
            if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commBackward_[timeLevel], 926, MPICommVars_.Comm, &receive_request_);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=MytimeSteps_[timeLevel]-1; ti>=0; ti--)
            {
                //    gsInfo<<"\t \t rank: "<<gsMpi::worldRank()<<"start with timestep "<<ti<<"\n";
                //compute f - A*u and store it in "vec.v[timeLevel][ti]"
                computeResidual(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);
                //    gsInfo<<"\t \t rank: "<<gsMpi::worldRank()<<"Compute res done!\n";
                //compute D^{-1}(f - A*u) (D^{-1} is maybe an approximation of the exact inverse) and store it in "vec.x[timeLevel][ti]"
                vec.x[timeLevel][ti].setZero();
                A_[spaceLevel]->ApproximateSolve(vec.x[timeLevel][ti], vec.v[timeLevel][ti], timeLevel, ti);
                //   gsInfo<<"\t \t rank: "<<gsMpi::worldRank()<<"Approx Solve done!\n";

                //MPI communicate forward
                if(ti==MytimeSteps_[timeLevel]-1)
                {
                    if(commForward_[timeLevel]>=0)
                    {
                        //    gsInfo<<"\t \t rank: "<<gsMpi::worldRank()<<"Do forward comm: u: "<<vec.u[timeLevel][ti].rows()<<" - "<<vec.u[timeLevel][ti].cols()<<" and x: "<<vec.x[timeLevel][ti].rows()<<" - "<<vec.x[timeLevel][ti].cols()<<"\n";
                        if(k>0) MPI_Wait(&send_request_, &send_status_);
                        for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][ti](l,0) + timeOmega_ * vec.x[timeLevel][ti](l,0);
                        //     gsInfo<<"\t \t rank: "<<gsMpi::worldRank()<<"Do Send!\n";
                        MPI_Isend(sendBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 926, MPICommVars_.Comm, &send_request_);
                    }
                }
            }
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=0; ti<MytimeSteps_[timeLevel]; ti++)
            {
                //compute u = u + omega * D^{-1}(f - A*u)
                vec.u[timeLevel][ti] += timeOmega_ * vec.x[timeLevel][ti];
            }

            //MPI receive value from the past
            if(commBackward_[timeLevel]>=0)
            {
                //       gsInfo<<"\t \t rank: "<<gsMpi::worldRank()<<"wait for values from the past!\n";
                MPI_Wait(&receive_request_, &receive_status_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) xold_[timeLevel](l,0) = receiveBuffer_[timeLevel][l];
            }

            if(commForward_[timeLevel]>=0)
            {
                MPI_Wait(&send_request_, &send_status_);
            }
        }
    }

    void SpaceTimeHybridSmoother(STSolVector &vec, const vector<Vector> &fl, int timeLevel, int sm=1)
    {
        int spaceLevel = spaceLevels_[timeLevel];
        int nTimeSteps = MytimeSteps_[timeLevel];

        //Jacobi/Gauß-Seidel - smoothing (Jacobi over the processors and Gauß-Seidel for each processor)
        for(int k=0; k<sm; k++)
        {
            if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commBackward_[timeLevel], 0, MPICommVars_.Comm, &receive_request_);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=0; ti<nTimeSteps; ti++)
            {
                //compute f - A*u and store it in "vec.v[timeLevel][ti]"
                computeResidual(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);
            }
            //solve forward in time
            for(int ti=0; ti<nTimeSteps; ti++)
            {
                vec.x[timeLevel][ti].setZero();
                A_[spaceLevel]->ApproximateSolve(vec.x[timeLevel][ti], vec.v[timeLevel][ti], timeLevel, ti);

                if(ti<nTimeSteps-1)
                {
                    A_[spaceLevel]->calcInitialVector(vec.x[timeLevel][ti], vec.v[timeLevel][ti+1], timeLevel, ti);
                }
            }
            if(commForward_[timeLevel]>=0)
            {
                int ti = MytimeSteps_[timeLevel]-1;
                if(k>0) MPI_Wait(&send_request_, &send_status_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][ti](l,0) + timeOmega_ * vec.x[timeLevel][ti](l,0);
                MPI_Isend(sendBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 0, MPICommVars_.Comm, &send_request_);
            }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=0; ti<MytimeSteps_[timeLevel]; ti++)
            {
                //compute u = u + omega * L^{-1}(f - A*u)
                //					vec.u[timeLevel][ti] *= (1.0-timeOmega_);
                vec.u[timeLevel][ti] += timeOmega_ * vec.x[timeLevel][ti];
            }

            //MPI receive value from the past
            if(commBackward_[timeLevel]>=0)
            {
                MPI_Wait(&receive_request_, &receive_status_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) xold_[timeLevel](l,0) = receiveBuffer_[timeLevel][l];
            }
        }
    }

    void SpaceTimeHybridSmootherWB(STSolVector &vec, const vector<Vector> &fl, int timeLevel, int sm=1, int seqSteps=-1)
    {
        int spaceLevel = spaceLevels_[timeLevel];
        int nTimeSteps = MytimeSteps_[timeLevel];

        int myseqSteps = seqSteps;
        if(seqSteps<0) myseqSteps = nTimeSteps;
        if(seqSteps>nTimeSteps) myseqSteps = nTimeSteps;

        //Jacobi/Gauß-Seidel - smoothing (Jacobi over the processors and Gauß-Seidel for each processor)
        for(int k=0; k<sm; k++)
        {
            if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commBackward_[timeLevel], 0, MPICommVars_.Comm, &receive_request_);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int iterTi=0; iterTi<nTimeSteps; iterTi+=myseqSteps)
            {
                for(int ti=iterTi; (ti<iterTi+myseqSteps && ti<nTimeSteps); ti++)
                {
                    //compute f - A*u and store it in "vec.v[timeLevel][ti]"
                    computeResidual(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);
                }
                //solve forward in time
                for(int ti=iterTi; (ti<iterTi+myseqSteps && ti<nTimeSteps); ti++)
                {
                    vec.x[timeLevel][ti].setZero();
                    A_[spaceLevel]->ApproximateSolve(vec.x[timeLevel][ti], vec.v[timeLevel][ti], timeLevel, ti);

                    if(ti<nTimeSteps-1)
                    {
                        A_[spaceLevel]->calcInitialVector(vec.x[timeLevel][ti], vec.v[timeLevel][ti+1], timeLevel, ti);
                    }
                }
            }
            if(commForward_[timeLevel]>=0)
            {
                int ti = MytimeSteps_[timeLevel]-1;
                if(k>0) MPI_Wait(&send_request_, &send_status_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][ti](l,0) + timeOmega_ * vec.x[timeLevel][ti](l,0);
                MPI_Isend(sendBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 0, MPICommVars_.Comm, &send_request_);
            }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=0; ti<MytimeSteps_[timeLevel]; ti++)
            {
                //compute u = u + omega * L^{-1}(f - A*u)
                vec.u[timeLevel][ti] += timeOmega_ * vec.x[timeLevel][ti];
            }

            //MPI receive value from the past
            if(commBackward_[timeLevel]>=0)
            {
                MPI_Wait(&receive_request_, &receive_status_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) xold_[timeLevel](l,0) = receiveBuffer_[timeLevel][l];
            }
        }
    }


    void SpaceTimeSmootherTriDiag(STSolVector &vec, const vector<Vector> &fl, int timeLevel, int sm=1)
    {
        int spaceLevel = spaceLevels_[timeLevel];

        //jacobi - smoothing
        for(int k=0; k<sm; k++)
        {
            if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commBackward_[timeLevel], 0, MPICommVars_.Comm, &receive_request_);
            if(commForward_[timeLevel]>=0) MPI_Irecv(receiveBufferFwd_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 22, MPICommVars_.Comm, &receive_requestFwd_);

            //compute the first and the last vector first and then communicate the results
            int ti = 0;
            computeResidualTriDiag(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);
            vec.x[timeLevel][ti].setZero();	A_[spaceLevel]->ApproximateSolve(vec.x[timeLevel][ti], vec.v[timeLevel][ti], timeLevel, ti);
            if(commBackward_[timeLevel]>=0)
            {
                if(k>0) MPI_Wait(&send_requestFwd_, &send_statusFwd_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBufferFwd_[timeLevel][l] = vec.u[timeLevel][ti](l,0) + timeOmega_ * vec.x[timeLevel][ti](l,0);
                MPI_Isend(sendBufferFwd_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commBackward_[timeLevel], 22, MPICommVars_.Comm, &send_requestFwd_);
            }

            if(MytimeSteps_[timeLevel]>1)
            {
                ti = MytimeSteps_[timeLevel]-1;
                computeResidualTriDiag(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);
                vec.x[timeLevel][ti].setZero();	A_[spaceLevel]->ApproximateSolve(vec.x[timeLevel][ti], vec.v[timeLevel][ti], timeLevel, ti);
            }
            if(commForward_[timeLevel]>=0)
            {
                if(k>0) MPI_Wait(&send_request_, &send_status_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][ti](l,0) + timeOmega_ * vec.x[timeLevel][ti](l,0);
                MPI_Isend(sendBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 0, MPICommVars_.Comm, &send_request_);
            }


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=1; ti<MytimeSteps_[timeLevel]-1; ti++)
            {
                //compute f - A*u and store it in "vec.v[timeLevel][ti]"
                computeResidualTriDiag(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);

                //compute D^{-1}(f - A*u) (D^{-1} is maybe an approximation of the exact inverse) and store it in "vec.x[timeLevel][ti]"
                vec.x[timeLevel][ti].setZero();
                A_[spaceLevel]->ApproximateSolve(vec.x[timeLevel][ti], vec.v[timeLevel][ti], timeLevel, ti);
            }
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=0; ti<MytimeSteps_[timeLevel]; ti++)
            {
                //compute u = u + omega * D^{-1}(f - A*u)
                vec.u[timeLevel][ti] += timeOmega_ * vec.x[timeLevel][ti];
            }

            //MPI receive value from the past
            if(commBackward_[timeLevel]>=0)
            {
                MPI_Wait(&receive_request_, &receive_status_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) xold_[timeLevel](l,0) = receiveBuffer_[timeLevel][l];
            }

            //MPI receive value from the future
            if(commForward_[timeLevel]>=0)
            {
                MPI_Wait(&receive_requestFwd_, &receive_statusFwd_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) xnew_[timeLevel](l,0) = receiveBufferFwd_[timeLevel][l];
            }
        }
    }


    void SpaceTimeRestrictToAllLevels(STSolVector &vec, int timeLevel, double factor=1.0)
    {
        if(timeLevel>0)
        {
            int nTimeSteps = MytimeSteps_[timeLevel];
            int nTimeStepsCoarse = MytimeSteps_[timeLevel-1];
            int spaceLevel = spaceLevels_[timeLevel];
            int spaceLevelCoarse = spaceLevels_[timeLevel-1];

            //calc restriction
            if(nTimeSteps>1)
            {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
                for(int tiC=0; tiC<nTimeStepsCoarse; tiC++)
                {
                    if(spaceLevel!=spaceLevelCoarse)
                    {
                        A_[spaceLevelCoarse]->SpaceRestriction(vec.u[timeLevel][2*tiC], vec.v[timeLevel-1][tiC], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.v[timeLevel-1][tiC], vec.u[timeLevel-1][tiC], timeLevel-1, true);

                        A_[spaceLevelCoarse]->SpaceRestriction(vec.u[timeLevel][2*tiC+1], vec.v[timeLevel-1][tiC], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.v[timeLevel-1][tiC], vec.u[timeLevel-1][tiC], timeLevel-1, false);

                    }else
                    {
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.u[timeLevel][2*tiC], vec.u[timeLevel-1][tiC], timeLevel-1, true);
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.u[timeLevel][2*tiC+1], vec.u[timeLevel-1][tiC], timeLevel-1, false);
                    }
                    vec.u[timeLevel-1][tiC]*=factor;
                }
            }else
            {
                if(activeTimeLevels_[timeLevel-1])
                {
                    //compute restriction1 in vec.w[timeLevel-1][0]
                    //receive from future to past the coarser vector and add this one to my computed coarser vector
                    //add the coarser vectors together
                    MPI_Irecv(receiveBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commForward_[timeLevel], 1, MPICommVars_.Comm, &receive_request_);

                    if(spaceLevel!=spaceLevelCoarse)
                    {
                        A_[spaceLevelCoarse]->SpaceRestriction(vec.u[timeLevel][0], vec.v[timeLevel-1][0], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.v[timeLevel-1][0], vec.u[timeLevel-1][0], timeLevel-1, true);
                    }else
                    {
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.u[timeLevel][0], vec.u[timeLevel-1][0], timeLevel-1, true);
                    }

                    MPI_Wait(&receive_request_, &receive_status_);
                    for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) vec.u[timeLevel-1][0](l,0) += receiveBuffer_[timeLevel-1][l];

                    vec.u[timeLevel-1][0] *= factor;

                }else
                {
                    //compute coarser vector and send it to the past
                    //then receive from the past the update and store it in vec.u[timeLevel-1][0]
                    if(spaceLevel!=spaceLevelCoarse)
                    {
                        A_[spaceLevelCoarse]->SpaceRestriction(vec.u[timeLevel][0], vec.v[timeLevel-1][0], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.v[timeLevel-1][0], vec.u[timeLevel-1][0], timeLevel-1, true);
                    }else
                    {
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.u[timeLevel][0], vec.u[timeLevel-1][0], timeLevel-1, true);
                    }
                    for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) sendBuffer_[timeLevel-1][l] = vec.u[timeLevel-1][0](l,0);
                    MPI_Send(sendBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commBackward_[timeLevel], 1, MPICommVars_.Comm);
                }
            }

            //go down to the next coarser level
            if(activeTimeLevels_[timeLevel-1]) SpaceTimeRestrictToAllLevels(vec, timeLevel-1, factor);
        }
    }

    void SpaceTimeRestriction(int timeLevel, vector<vector<Vector> > &v, vector<vector<Vector> > &w)
    {
        int nTimeSteps = MytimeSteps_[timeLevel];
        int nTimeStepsCoarse = MytimeSteps_[timeLevel-1];
        int spaceLevel = spaceLevels_[timeLevel];
        int spaceLevelCoarse = spaceLevels_[timeLevel-1];

        //calc restriction
        if(nTimeSteps>1)
        {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int tiC=0; tiC<nTimeStepsCoarse; tiC++)
            {
                if(spaceLevel!=spaceLevelCoarse)
                {
                    A_[spaceLevelCoarse]->SpaceRestriction(v[timeLevel][2*tiC], v[timeLevel-1][tiC], spaceLevelCoarse, true);
                    A_[spaceLevelCoarse]->TimeRestriction1(v[timeLevel-1][tiC], w[timeLevel-1][tiC], timeLevel-1, true);

                    A_[spaceLevelCoarse]->SpaceRestriction(v[timeLevel][2*tiC+1], v[timeLevel-1][tiC], spaceLevelCoarse, true);
                    A_[spaceLevelCoarse]->TimeRestriction2(v[timeLevel-1][tiC], w[timeLevel-1][tiC], timeLevel-1, false);

                }else
                {
                    A_[spaceLevelCoarse]->TimeRestriction1(v[timeLevel][2*tiC], w[timeLevel-1][tiC], timeLevel-1, true);
                    A_[spaceLevelCoarse]->TimeRestriction2(v[timeLevel][2*tiC+1], w[timeLevel-1][tiC], timeLevel-1, false);
                }
                //					vec.u[timeLevel-1][tiC] = 0.0;
                //					xold_[timeLevel-1] = 0.0;
            }
        }else
        {
            if(activeTimeLevels_[timeLevel-1])
            {
                //compute restriction1 in vec.w[timeLevel-1][0]
                //receive from future to past the coarser vector and add this one to my computed coarser vector
                //add the coarser vectors together
                MPI_Irecv(receiveBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commForward_[timeLevel], 1, MPICommVars_.Comm, &receive_request_);

                if(spaceLevel!=spaceLevelCoarse)
                {
                    A_[spaceLevelCoarse]->SpaceRestriction(v[timeLevel][0], v[timeLevel-1][0], spaceLevelCoarse, true);
                    A_[spaceLevelCoarse]->TimeRestriction1(v[timeLevel-1][0], w[timeLevel-1][0], timeLevel-1, true);
                }else
                {
                    A_[spaceLevelCoarse]->TimeRestriction1(v[timeLevel][0], w[timeLevel-1][0], timeLevel-1, true);
                }
                //					vec.u[timeLevel-1][0] = 0.0;
                //					xold_[timeLevel-1] = 0.0;

                MPI_Wait(&receive_request_, &receive_status_);
                for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) w[timeLevel-1][0](l,0) += receiveBuffer_[timeLevel-1][l];
            }else if(commBackward_[timeLevel]>=0)
            {
                //compute coarser vector and send it to the past
                //then receive from the past the update and store it in vec.u[timeLevel-1][0]
                if(spaceLevel!=spaceLevelCoarse)
                {
                    A_[spaceLevelCoarse]->SpaceRestriction(v[timeLevel][0], v[timeLevel-1][0], spaceLevelCoarse, true);
                    A_[spaceLevelCoarse]->TimeRestriction2(v[timeLevel-1][0], w[timeLevel-1][0], timeLevel-1, true);
                }else
                {
                    A_[spaceLevelCoarse]->TimeRestriction2(v[timeLevel][0], w[timeLevel-1][0], timeLevel-1, true);
                }
                for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) sendBuffer_[timeLevel-1][l] = w[timeLevel-1][0](l,0);
                MPI_Send(sendBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commBackward_[timeLevel], 1, MPICommVars_.Comm);

                //					MPI_Recv(receiveBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commBackward_[timeLevel], 2, MPICommVars_.Comm, &receive_status_);
                //					for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) vec.u[timeLevel-1][0][l] = receiveBuffer_[timeLevel-1][l];
            }
        }
    }

    void SpaceTimeProlongation(int timeLevel, vector<vector<Vector> > &u, vector<vector<Vector> > &v, bool exchangeXoldBuffer = false)
    {
        int nTimeSteps = MytimeSteps_[timeLevel];
        int nTimeStepsCoarse = MytimeSteps_[timeLevel-1];
        int spaceLevel = spaceLevels_[timeLevel];
        int spaceLevelCoarse = spaceLevels_[timeLevel-1];

        if(exchangeXoldBuffer) if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commBackward_[timeLevel], 123, MPICommVars_.Comm, &receive_request_);


        //prolongation and correction
        if(nTimeSteps>1)
        {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int tiC=nTimeStepsCoarse-1; tiC>=0; tiC--)
            {
                A_[spaceLevelCoarse]->TimeProlongation1(u[timeLevel-1][tiC], v[timeLevel-1][tiC], timeLevel-1, true);
                if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(v[timeLevel-1][tiC], u[timeLevel][2*tiC], spaceLevelCoarse, false);
                else u[timeLevel][2*tiC] += v[timeLevel-1][tiC];

                A_[spaceLevelCoarse]->TimeProlongation2(u[timeLevel-1][tiC], v[timeLevel-1][tiC], timeLevel-1, true);
                if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(v[timeLevel-1][tiC], u[timeLevel][2*tiC+1], spaceLevelCoarse, false);
                else u[timeLevel][2*tiC+1] += v[timeLevel-1][tiC];

                if(exchangeXoldBuffer)
                {
                    if(commForward_[timeLevel]>=0 && tiC==nTimeStepsCoarse-1)
                    {
                        for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBuffer_[timeLevel][l] = u[timeLevel][MytimeSteps_[timeLevel]-1](l,0);
                        MPI_Isend(sendBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 123, MPICommVars_.Comm, &send_request_);
                    }
                }
            }
        }else
        {
            if(activeTimeLevels_[timeLevel-1])
            {
                //send update vector vec.u[timeLevel-1][0] to the future
                //compute prolongation1 and correction with the vector vec.u[timeLevel-1][0]
                for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) sendBuffer_[timeLevel-1][l] = u[timeLevel-1][0](l,0);
                MPI_Send(sendBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commForward_[timeLevel], 2, MPICommVars_.Comm);

                A_[spaceLevelCoarse]->TimeProlongation1(u[timeLevel-1][0], v[timeLevel-1][0], timeLevel-1, true);
                if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(v[timeLevel-1][0], u[timeLevel][0], spaceLevelCoarse, false);
                else u[timeLevel][0] += v[timeLevel-1][0];
            }else if(commBackward_[timeLevel]>=0)
            {
                //					cout << MPICommVars_.my_rank << " receive from " << commBackward_[timeLevel] << endl;

                MPI_Recv(receiveBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commBackward_[timeLevel], 2, MPICommVars_.Comm, &receive_status_);
                for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) u[timeLevel-1][0](l,0) = receiveBuffer_[timeLevel-1][l];

                //compute prolongation2 and correction with the vector vec.u[timeLevel-1][0]
                A_[spaceLevelCoarse]->TimeProlongation2(u[timeLevel-1][0], v[timeLevel-1][0], timeLevel-1, true);
                if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(v[timeLevel-1][0], u[timeLevel][0], spaceLevelCoarse, false);
                else u[timeLevel][0] += v[timeLevel-1][0];
            }

            if(exchangeXoldBuffer)
            {
                if(commForward_[timeLevel]>=0)
                {
                    for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBuffer_[timeLevel][l] = u[timeLevel][0](l,0);
                    MPI_Isend(sendBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 123, MPICommVars_.Comm, &send_request_);
                }
            }

        }

        //wait until the xold's are exchanged
        if(exchangeXoldBuffer)
        {
            if(commBackward_[timeLevel]>=0)
            {
                MPI_Wait(&receive_request_, &receive_status_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) xold_[timeLevel](l,0) = receiveBuffer_[timeLevel][l];
            }
            if(commForward_[timeLevel]>=0) MPI_Wait(&send_request_, &send_status_);
        }
    }



    void SpaceTimeMGMSchritt(STSolVector &vec, const vector<Vector> &fl, int timeLevel)
    {
        if(timeLevel==coarseLevel_)
        {
            STForwardSolve(vec, fl, timeLevel);
        }else
        {

            int nTimeSteps = MytimeSteps_[timeLevel];
            int nTimeStepsCoarse = MytimeSteps_[timeLevel-1];
            int spaceLevel = spaceLevels_[timeLevel];
            int spaceLevelCoarse = spaceLevels_[timeLevel-1];
            //  gsInfo<<"Rank "<<MPICommVars_.my_rank<<" start smoothing\n"<<std::flush;
            // MPI_Barrier(MPICommVars_.Comm);
            //pre smoothing
            if(useHybridSmoother_) SpaceTimeHybridSmootherWB(vec, fl, timeLevel, timeSm1_,8);
            else SpaceTimeSmoother(vec, fl, timeLevel, timeSm1_);
            //   gsInfo<<"Rank "<<MPICommVars_.my_rank<<" start residual\n"<<std::flush;
            //   MPI_Barrier(MPICommVars_.Comm);
            //compute defect
#ifdef USE_OPENMP
#pragma omp parallel for
#endif

            for(int ti=0; ti<nTimeSteps; ti++)
                computeResidual(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);

            // gsInfo<<"Rank "<<MPICommVars_.my_rank<<" calc restriction\n"<<std::flush;
            //  MPI_Barrier(MPICommVars_.Comm);
            //calc restriction
            if(nTimeSteps>1)
            {
                //gsInfo<<"Rank "<<MPICommVars_.my_rank<<" Case 1\n"<<std::flush;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
                for(int tiC=0; tiC<nTimeStepsCoarse; tiC++)
                {
                    if(spaceLevel!=spaceLevelCoarse)
                    {
                        A_[spaceLevelCoarse]->SpaceRestriction(vec.v[timeLevel][2*tiC], vec.v[timeLevel-1][tiC], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.v[timeLevel-1][tiC], vec.w[timeLevel-1][tiC], timeLevel-1, true);

                        A_[spaceLevelCoarse]->SpaceRestriction(vec.v[timeLevel][2*tiC+1], vec.v[timeLevel-1][tiC], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.v[timeLevel-1][tiC], vec.w[timeLevel-1][tiC], timeLevel-1, false);

                    }else
                    {
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.v[timeLevel][2*tiC], vec.w[timeLevel-1][tiC], timeLevel-1, true);
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.v[timeLevel][2*tiC+1], vec.w[timeLevel-1][tiC], timeLevel-1, false);
                    }
                    vec.u[timeLevel-1][tiC].setZero();
                    xold_[timeLevel-1].setZero();
                }
            }else
            {
                //  gsInfo<<"Rank "<<MPICommVars_.my_rank<<" Case 2\n"<<std::flush;
                if(activeTimeLevels_[timeLevel-1])
                {
                    //    gsInfo<<"Rank "<<MPICommVars_.my_rank<<" Case 2a\n"<<std::flush;
                    //compute restriction1 in vec.w[timeLevel-1][0]
                    //receive from future to past the coarser vector and add this one to my computed coarser vector
                    //add the coarser vectors together
                    MPI_Irecv(receiveBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commForward_[timeLevel], 1, MPICommVars_.Comm, &receive_request_);

                    if(spaceLevel!=spaceLevelCoarse)
                    {
                        A_[spaceLevelCoarse]->SpaceRestriction(vec.v[timeLevel][0], vec.v[timeLevel-1][0], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.v[timeLevel-1][0], vec.w[timeLevel-1][0], timeLevel-1, true);
                    }else
                    {
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.v[timeLevel][0], vec.w[timeLevel-1][0], timeLevel-1, true);
                    }
                    vec.u[timeLevel-1][0].setZero();
                    xold_[timeLevel-1].setZero();

                    MPI_Wait(&receive_request_, &receive_status_);
                    for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) vec.w[timeLevel-1][0](l,0) += receiveBuffer_[timeLevel-1][l];
                }else
                {
                    //   gsInfo<<"Rank "<<MPICommVars_.my_rank<<" Case 2b\n"<<std::flush;
                    //compute coarser vector and send it to the past
                    //then receive from the past the update and store it in vec.u[timeLevel-1][0]
                    if(spaceLevel!=spaceLevelCoarse)
                    {
                        A_[spaceLevelCoarse]->SpaceRestriction(vec.v[timeLevel][0], vec.v[timeLevel-1][0], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.v[timeLevel-1][0], vec.w[timeLevel-1][0], timeLevel-1, true);
                    }else
                    {
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.v[timeLevel][0], vec.w[timeLevel-1][0], timeLevel-1, true);
                    }
                    for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) sendBuffer_[timeLevel-1][l] = vec.w[timeLevel-1][0](l,0);
                    MPI_Send(sendBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commBackward_[timeLevel], 1, MPICommVars_.Comm);


                    MPI_Recv(receiveBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commBackward_[timeLevel], 2, MPICommVars_.Comm, &receive_status_);
                    for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) vec.u[timeLevel-1][0](l,0) = receiveBuffer_[timeLevel-1][l];
                }
            }
            // gsInfo<<"Rank "<<MPICommVars_.my_rank<<" start coarser\n"<<std::flush;
            // MPI_Barrier(MPICommVars_.Comm);
            //multigrid step on the coarser level
            if(activeTimeLevels_[timeLevel-1]) for(int k=0; k<timeCycles_; k++) SpaceTimeMGMSchritt(vec, vec.w[timeLevel-1], timeLevel-1);

            if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commBackward_[timeLevel], 123, MPICommVars_.Comm, &receive_request_);


            //gsInfo<<"Rank "<<MPICommVars_.my_rank<<" start prolong\n"<<std::flush;
            //MPI_Barrier(MPICommVars_.Comm);
            //prolongation and correction
            if(nTimeSteps>1)
            {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
                for(int tiC=nTimeStepsCoarse-1; tiC>=0; tiC--)
                {
                    A_[spaceLevelCoarse]->TimeProlongation1(vec.u[timeLevel-1][tiC], vec.v[timeLevel-1][tiC], timeLevel-1, true);
                    if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][tiC], vec.u[timeLevel][2*tiC], spaceLevelCoarse, false);
                    else vec.u[timeLevel][2*tiC] += vec.v[timeLevel-1][tiC];

                    A_[spaceLevelCoarse]->TimeProlongation2(vec.u[timeLevel-1][tiC], vec.v[timeLevel-1][tiC], timeLevel-1, true);
                    if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][tiC], vec.u[timeLevel][2*tiC+1], spaceLevelCoarse, false);
                    else vec.u[timeLevel][2*tiC+1] += vec.v[timeLevel-1][tiC];

                    if(commForward_[timeLevel]>=0 && tiC==nTimeStepsCoarse-1)
                    {
                        for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][MytimeSteps_[timeLevel]-1](l,0);
                        MPI_Isend(sendBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 123, MPICommVars_.Comm, &send_request_);
                    }
                }
            }else
            {
                if(activeTimeLevels_[timeLevel-1])
                {
                    //send update vector vec.u[timeLevel-1][0] to the future
                    //compute prolongation1 and correction with the vector vec.u[timeLevel-1][0]
                    for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) sendBuffer_[timeLevel-1][l] = vec.u[timeLevel-1][0](l,0);
                    MPI_Isend(sendBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commForward_[timeLevel], 2, MPICommVars_.Comm, &send_request_);

                    A_[spaceLevelCoarse]->TimeProlongation1(vec.u[timeLevel-1][0], vec.v[timeLevel-1][0], timeLevel-1, true);
                    if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][0], vec.u[timeLevel][0], spaceLevelCoarse, false);
                    else vec.u[timeLevel][0] += vec.v[timeLevel-1][0];
                }else
                {
                    //compute prolongation2 and correction with the vector vec.u[timeLevel-1][0]
                    A_[spaceLevelCoarse]->TimeProlongation2(vec.u[timeLevel-1][0], vec.v[timeLevel-1][0], timeLevel-1, true);
                    if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][0], vec.u[timeLevel][0], spaceLevelCoarse, false);
                    else vec.u[timeLevel][0] += vec.v[timeLevel-1][0];
                }

                if(commForward_[timeLevel]>=0)
                {
                    for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][0](l,0);
                    MPI_Isend(sendBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 123, MPICommVars_.Comm, &send_request_);
                }

            }

            //wait until the xold's are exchanged
            if(commBackward_[timeLevel]>=0)
            {
                MPI_Wait(&receive_request_, &receive_status_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) xold_[timeLevel](l,0) = receiveBuffer_[timeLevel][l];
            }
            if(commForward_[timeLevel]>=0) MPI_Wait(&send_request_, &send_status_);


            //      gsInfo<<"Rank "<<MPICommVars_.my_rank<<" start postsmooth\n"<<std::flush;
            //     MPI_Barrier(MPICommVars_.Comm);
            //post smoothing
            if(useHybridSmoother_) SpaceTimeHybridSmootherWB(vec, fl, timeLevel, timeSm2_,2);
            else SpaceTimeSmoother(vec, fl, timeLevel, timeSm2_);
        }
    }



    void SpaceTimeMGMSchrittTriDiag(STSolVector &vec, const vector<Vector> &fl, int timeLevel)
    {
        if(timeLevel==coarseLevel_)
        {
            if(coarseLevel_==0)
            {
                STForwardSolve(vec, fl, timeLevel);
            }else
            {
                int tempCoarseLevel = coarseLevel_;
                coarseLevel_ = 0;

                for(int k=0; k<20; k++) SpaceTimeMGMSchrittTriDiag(vec, fl, timeLevel);

                coarseLevel_ = tempCoarseLevel;
            }

        }else
        {

            int nTimeSteps = MytimeSteps_[timeLevel];
            int nTimeStepsCoarse = MytimeSteps_[timeLevel-1];
            int spaceLevel = spaceLevels_[timeLevel];
            int spaceLevelCoarse = spaceLevels_[timeLevel-1];

            //pre smoothing
            SpaceTimeSmootherTriDiag(vec, fl, timeLevel, timeSm1_);

            //compute defect
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=0; ti<nTimeSteps; ti++)
            {
                computeResidualTriDiag(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);
            }


            //calc restriction
            if(nTimeSteps>1)
            {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
                for(int tiC=0; tiC<nTimeStepsCoarse; tiC++)
                {
                    if(spaceLevel!=spaceLevelCoarse)
                    {
                        A_[spaceLevelCoarse]->SpaceRestriction(vec.v[timeLevel][2*tiC], vec.v[timeLevel-1][tiC], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.v[timeLevel-1][tiC], vec.w[timeLevel-1][tiC], timeLevel-1, true);

                        A_[spaceLevelCoarse]->SpaceRestriction(vec.v[timeLevel][2*tiC+1], vec.v[timeLevel-1][tiC], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.v[timeLevel-1][tiC], vec.w[timeLevel-1][tiC], timeLevel-1, false);

                    }else
                    {
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.v[timeLevel][2*tiC], vec.w[timeLevel-1][tiC], timeLevel-1, true);
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.v[timeLevel][2*tiC+1], vec.w[timeLevel-1][tiC], timeLevel-1, false);
                    }
                    vec.u[timeLevel-1][tiC].setZero();
                    xold_[timeLevel-1].setZero();
                    xnew_[timeLevel-1].setZero();
                }
            }else
            {
                if(activeTimeLevels_[timeLevel-1])
                {
                    //compute restriction1 in vec.w[timeLevel-1][0]
                    //receive from future to past the coarser vector and add this one to my computed coarser vector
                    //add the coarser vectors together
                    MPI_Irecv(receiveBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commForward_[timeLevel], 1, MPICommVars_.Comm, &receive_request_);

                    if(spaceLevel!=spaceLevelCoarse)
                    {
                        A_[spaceLevelCoarse]->SpaceRestriction(vec.v[timeLevel][0], vec.v[timeLevel-1][0], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.v[timeLevel-1][0], vec.w[timeLevel-1][0], timeLevel-1, true);
                    }else
                    {
                        A_[spaceLevelCoarse]->TimeRestriction1(vec.v[timeLevel][0], vec.w[timeLevel-1][0], timeLevel-1, true);
                    }
                    vec.u[timeLevel-1][0].setZero();
                    xold_[timeLevel-1].setZero();
                    xnew_[timeLevel-1].setZero();

                    MPI_Wait(&receive_request_, &receive_status_);
                    for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) vec.w[timeLevel-1][0](l,0) += receiveBuffer_[timeLevel-1][l];
                }else
                {
                    //compute coarser vector and send it to the past
                    //then receive from the past the update and store it in vec.u[timeLevel-1][0]
                    if(spaceLevel!=spaceLevelCoarse)
                    {
                        A_[spaceLevelCoarse]->SpaceRestriction(vec.v[timeLevel][0], vec.v[timeLevel-1][0], spaceLevelCoarse, true);
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.v[timeLevel-1][0], vec.w[timeLevel-1][0], timeLevel-1, true);
                    }else
                    {
                        A_[spaceLevelCoarse]->TimeRestriction2(vec.v[timeLevel][0], vec.w[timeLevel-1][0], timeLevel-1, true);
                    }
                    for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) sendBuffer_[timeLevel-1][l] = vec.w[timeLevel-1][0](l,0);
                    MPI_Send(sendBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commBackward_[timeLevel], 1, MPICommVars_.Comm);


                    MPI_Recv(receiveBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commBackward_[timeLevel], 2, MPICommVars_.Comm, &receive_status_);
                    for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) vec.u[timeLevel-1][0](l,0) = receiveBuffer_[timeLevel-1][l];
                }
            }

            //multigrid step on the coarser level
            if(activeTimeLevels_[timeLevel-1]) for(int k=0; k<timeCycles_; k++) SpaceTimeMGMSchrittTriDiag(vec, vec.w[timeLevel-1], timeLevel-1);

            if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commBackward_[timeLevel], 123, MPICommVars_.Comm, &receive_request_);
            if(commForward_[timeLevel]>=0) MPI_Irecv(receiveBufferFwd_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 456, MPICommVars_.Comm, &receive_requestFwd_);

            //prolongation and correcton
            if(nTimeSteps>1)
            {
                int tiC=nTimeStepsCoarse-1;
                A_[spaceLevelCoarse]->TimeProlongation1(vec.u[timeLevel-1][tiC], vec.v[timeLevel-1][tiC], timeLevel-1, true);
                if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][tiC], vec.u[timeLevel][2*tiC], spaceLevelCoarse, false);
                else vec.u[timeLevel][2*tiC] += vec.v[timeLevel-1][tiC];

                A_[spaceLevelCoarse]->TimeProlongation2(vec.u[timeLevel-1][tiC], vec.v[timeLevel-1][tiC], timeLevel-1, true);
                if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][tiC], vec.u[timeLevel][2*tiC+1], spaceLevelCoarse, false);
                else vec.u[timeLevel][2*tiC+1] += vec.v[timeLevel-1][tiC];

                if(commForward_[timeLevel]>=0)
                {
                    for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][MytimeSteps_[timeLevel]-1](l,0);
                    MPI_Isend(sendBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 123, MPICommVars_.Comm, &send_request_);
                }

                if(nTimeStepsCoarse>1)
                {
                    tiC=0;
                    A_[spaceLevelCoarse]->TimeProlongation1(vec.u[timeLevel-1][tiC], vec.v[timeLevel-1][tiC], timeLevel-1, true);
                    if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][tiC], vec.u[timeLevel][2*tiC], spaceLevelCoarse, false);
                    else vec.u[timeLevel][2*tiC] += vec.v[timeLevel-1][tiC];

                    A_[spaceLevelCoarse]->TimeProlongation2(vec.u[timeLevel-1][tiC], vec.v[timeLevel-1][tiC], timeLevel-1, true);
                    if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][tiC], vec.u[timeLevel][2*tiC+1], spaceLevelCoarse, false);
                    else vec.u[timeLevel][2*tiC+1] += vec.v[timeLevel-1][tiC];
                }

                if(commBackward_[timeLevel]>=0)
                {
                    for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBufferFwd_[timeLevel][l] = vec.u[timeLevel][0](l,0);
                    MPI_Isend(sendBufferFwd_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commBackward_[timeLevel], 456, MPICommVars_.Comm, &send_requestFwd_);
                }



#ifdef USE_OPENMP
#pragma omp parallel for
#endif
                for(int tiC=1; tiC<nTimeStepsCoarse-1; tiC++)
                {
                    A_[spaceLevelCoarse]->TimeProlongation1(vec.u[timeLevel-1][tiC], vec.v[timeLevel-1][tiC], timeLevel-1, true);
                    if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][tiC], vec.u[timeLevel][2*tiC], spaceLevelCoarse, false);
                    else vec.u[timeLevel][2*tiC] += vec.v[timeLevel-1][tiC];

                    A_[spaceLevelCoarse]->TimeProlongation2(vec.u[timeLevel-1][tiC], vec.v[timeLevel-1][tiC], timeLevel-1, true);
                    if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][tiC], vec.u[timeLevel][2*tiC+1], spaceLevelCoarse, false);
                    else vec.u[timeLevel][2*tiC+1] += vec.v[timeLevel-1][tiC];
                }
            }else
            {
                if(activeTimeLevels_[timeLevel-1])
                {
                    //send update vector vec.u[timeLevel-1][0] to the future
                    //compute prolongation1 and correction with the vector vec.u[timeLevel-1][0]
                    for(int l=0; l<spaceSizes_[spaceLevelCoarse]; l++) sendBuffer_[timeLevel-1][l] = vec.u[timeLevel-1][0](l,0);
                    MPI_Isend(sendBuffer_[timeLevel-1], spaceSizes_[spaceLevelCoarse], MPI_DOUBLE, commForward_[timeLevel], 2, MPICommVars_.Comm, &send_request_);

                    A_[spaceLevelCoarse]->TimeProlongation1(vec.u[timeLevel-1][0], vec.v[timeLevel-1][0], timeLevel-1, true);
                    if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][0], vec.u[timeLevel][0], spaceLevelCoarse, false);
                    else vec.u[timeLevel][0] += vec.v[timeLevel-1][0];
                }else
                {
                    //compute prolongation2 and correction with the vector vec.u[timeLevel-1][0]
                    A_[spaceLevelCoarse]->TimeProlongation2(vec.u[timeLevel-1][0], vec.v[timeLevel-1][0], timeLevel-1, true);
                    if(spaceLevel!=spaceLevelCoarse) A_[spaceLevelCoarse]->SpaceProlongation(vec.v[timeLevel-1][0], vec.u[timeLevel][0], spaceLevelCoarse, false);
                    else vec.u[timeLevel][0] += vec.v[timeLevel-1][0];
                }


                if(commForward_[timeLevel]>=0)
                {
                    for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][0](l,0);
                    MPI_Isend(sendBuffer_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commForward_[timeLevel], 123, MPICommVars_.Comm, &send_request_);
                }

                if(commBackward_[timeLevel]>=0)
                {
                    for(int l=0; l<spaceSizes_[spaceLevel]; l++) sendBufferFwd_[timeLevel][l] = vec.u[timeLevel][0](l,0);
                    MPI_Isend(sendBufferFwd_[timeLevel], spaceSizes_[spaceLevel], MPI_DOUBLE, commBackward_[timeLevel], 456, MPICommVars_.Comm, &send_requestFwd_);
                }
            }

            //wait until the xold's are exchanged
            if(commBackward_[timeLevel]>=0)
            {
                MPI_Wait(&receive_request_, &receive_status_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) xold_[timeLevel](l,0) = receiveBuffer_[timeLevel][l];
            }
            //wait until the xnew's are exchanged
            if(commForward_[timeLevel]>=0)
            {
                MPI_Wait(&receive_requestFwd_, &receive_statusFwd_);
                for(int l=0; l<spaceSizes_[spaceLevel]; l++) xnew_[timeLevel](l,0) = receiveBufferFwd_[timeLevel][l];
            }


            if(commForward_[timeLevel]>=0) MPI_Wait(&send_request_, &send_status_);
            if(commBackward_[timeLevel]>=0) MPI_Wait(&send_requestFwd_, &send_statusFwd_);



            //post smoothing
            SpaceTimeSmootherTriDiag(vec, fl, timeLevel, timeSm2_);
        }
    }

    double ComputeSpaceTimeNorm(const vector<Vector> &fl, int timeLevel)
    {
        //compute residuum
        double localError = 0.0;

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localError)
#endif
        for(int ti=0; ti<MytimeSteps_[timeLevel]; ti++)
        {
            localError += fl[ti].col(0).dot(fl[ti].col(0));
        }
        double error = -1.0;
        MPI_Allreduce(&localError, &error, 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);

        return sqrt(error);
    }

    double ComputeSpaceTimeResiduum(STSolVector &vec, const vector<Vector> &fl, int timeLevel)
    {
        //compute residuum
        double localError = 0.0;

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localError)
#endif
        for(int ti=0; ti<MytimeSteps_[timeLevel]; ti++)
        {
            computeResidual(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);
            localError += vec.v[timeLevel][ti].col(0).dot(vec.v[timeLevel][ti].col(0));
        }
        double error = -1.0;
        MPI_Allreduce(&localError, &error, 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);

        return sqrt(error);
    }

    double ComputeSpaceTimeResiduum(STSolVector &vec, const vector<Vector> &fl, int timeLevel, MPI_Comm SpaceComm)
    {
        //compute residuum
        double localError = 0.0;

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localError)
#endif
        Vector v_acc;
        int sRank;
        MPI_Comm_rank(SpaceComm, &sRank);
        for(int ti=0; ti<MytimeSteps_[timeLevel]; ti++)
        {
            computeResidual(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);
            A_[spaceLevels_[timeLevel]]->accumulate(vec.v[timeLevel][ti],v_acc,SpaceComm);
            localError += vec.v[timeLevel][ti].col(0).dot(v_acc.col(0));
        }
        double errorTloacal = -1.0;

        //sum up the error for each space-component
        MPI_Allreduce(&localError, &errorTloacal, 1, MPI_DOUBLE, MPI_SUM, SpaceComm);

        //sum up the error for each time-component
        double error = -1;
        MPI_Allreduce(&errorTloacal, &error, 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);

        return sqrt(math::abs(error));
    }

    double ComputeSpaceTimeResiduumTriDiag(STSolVector &vec, const vector<Vector> &fl, int timeLevel)
    {
        //compute residuum
        double localError = 0.0;

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localError)
#endif
        for(int ti=0; ti<MytimeSteps_[timeLevel]; ti++)
        {
            computeResidualTriDiag(vec.u[timeLevel], fl[ti], vec.v[timeLevel][ti], timeLevel, ti);
            localError += vec.v[timeLevel][ti].col(0).dot(vec.v[timeLevel][ti].col(0));
        }
        double error = -1.0;
        MPI_Allreduce(&localError, &error, 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);

        return sqrt(error);
    }

    double computeL2Difference(const vector<Vector> &u, const vector<Vector> &v, int timeLevel)
    {
        double localError = 0.0;

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localError)
#endif
        for(int ti=0; ti<MytimeSteps_[timeLevel]; ti++)
        {
            localError += (u[ti]-v[ti]).squaredNorm();
        }

        double error = -1.0;
        MPI_Allreduce(&localError, &error, 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);

        return sqrt(error);
    }


    int TimeMGM_MPI(STSolVector &vec, const vector<Vector> &fl, int maxit=100, double tol=1e-6, int useTimeLevel = -1, bool output=true, bool calcResiduum=true, bool zeroXResidualTest=false)
    {
        int timeLevel = NTimeLevels_-1;
        if(useTimeLevel>=0) timeLevel = useTimeLevel;

        int nTimeSteps = MytimeSteps_[timeLevel];

        double error0 = 1.0;
        if(zeroXResidualTest) error0 = ComputeSpaceTimeNorm(fl, timeLevel);
        else error0 = ComputeSpaceTimeResiduum(vec, fl, timeLevel);

        if(output) gsInfo << "\tInitial Error: " << error0 << std::endl;

        double rho = 1.0;
        double errorOld = error0;
        double error = 1.0;

        vector<double> rhos(maxit);
        double maxRho = 0.0;

        int i = 0;
        for(i=0; i<maxit; i++)
        {
            if(nTimeSteps>0) SpaceTimeMGMSchritt(vec, fl, timeLevel);
            // gsInfo<<"Rank "<<MPICommVars_.my_rank<<" finished step\n"<<std::flush;
            // MPI_Barrier(MPICommVars_.Comm);
            if(calcResiduum || output)
            {
                error = ComputeSpaceTimeResiduum(vec, fl, timeLevel);

                rho = error/errorOld;
                errorOld = error;

                if(MPICommVars_.my_rank==0) if(output) gsInfo << "\t\t\t" << i+1 << ": " << error/error0 << "  -->  " << rho  << " (" <<  pow(error/error0,1.0/double(i)) << ")" << std::endl;

                rhos[i] = rho;
                if(maxRho<rho && i>0) maxRho = rho;

                if(error/error0 < tol || error/error0>1e6) break;
            }
        }

        if(MPICommVars_.my_rank==0) if(output)
            gsInfo << "\t maximal contraction factor: " << maxRho << std::endl <<
                      "\t average contraction factor: " << pow(error/error0,1.0/double(i-1)) << std::endl;

        return i+1;
    }

    int TimeMGM_MPI(STSolVector &vec, const vector<Vector> &fl, MPI_Comm SpaceComm, int maxit=100, double tol=1e-6,int useTimeLevel=-1, bool output=true, bool calcResiduum=true)
    {
        int timeLevel = NTimeLevels_-1;
        if(useTimeLevel>=0) timeLevel = useTimeLevel;
        int nTimeSteps = MytimeSteps_[timeLevel];

        double error0 = 1.0;
        error0 = ComputeSpaceTimeResiduum(vec, fl, timeLevel,SpaceComm);
        if(output) gsInfo << "\tInitial Error: " << error0 << std::endl;

        double rho = 1.0;
        double errorOld = error0;
        double error = 1.0;

        vector<double> rhos(maxit);
        double maxRho = 0.0;



        int i = 0;
        for(i=0; i<maxit; i++)
        {
            if(nTimeSteps>0) SpaceTimeMGMSchritt(vec, fl, timeLevel);

            if(calcResiduum || output)
            {
                error = ComputeSpaceTimeResiduum(vec, fl, timeLevel, SpaceComm);
                rho = error/errorOld;
                errorOld = error;

                if(MPICommVars_.my_rank==0) if(output) gsInfo << "\t\t\t" << i+1 << ": " << error/error0 << "  -->  " << rho  << " (" <<  pow(error/error0,1.0/double(i)) << ")" << std::endl;

                rhos[i] = rho;
                if(maxRho<rho && i>0) maxRho = rho;

                if(error/error0 < tol || error/error0>1e6) break;
            }
        }

        if(MPICommVars_.my_rank==0) if(output)
            gsInfo << "\t maximal contraction factor: " << maxRho << std::endl <<
                      "\t average contraction factor: " << pow(error/error0,1.0/double(i-1)) << std::endl;

        return i+1;
    }



    unsigned int PGMRESSolve(STSolVector &vec, const vector<Vector> &fl,MPI_Comm SpaceComm,  int maxit=100, double tol=1e-6, bool output=true, bool calcResiduum=true, bool setZero = true)
    {
        if(maxit==0) return 0;

        int timeLevel = NTimeLevels_-1;

        int nTimeSteps = MytimeSteps_[timeLevel];
        int spaceLevel = spaceLevels_[timeLevel];
        int N = spaceSizes_[spaceLevel];

        //initial solution
        if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], N, MPI_DOUBLE, commBackward_[timeLevel], 789, MPICommVars_.Comm, &receive_request_);
        vector<Vector> x(nTimeSteps);
        Vector xOld(N,1); xOld.setZero();

        //communicate the the values between the processors
        if(commForward_[timeLevel]>=0)
        {
            if(setZero) for(int l=0; l<N; l++) sendBuffer_[timeLevel][l] = 0.0;
            else for(int l=0; l<N; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][nTimeSteps-1](l,0);
            MPI_Isend(sendBuffer_[timeLevel], N, MPI_DOUBLE, commForward_[timeLevel], 789, MPICommVars_.Comm, &send_request_);
        }

        if(setZero) for(int ti=0; ti<nTimeSteps; ti++) { x[ti].setZero(N,1); }
        else for(int ti=0; ti<nTimeSteps; ti++) x[ti] = vec.u[timeLevel][ti];

        //wait until the xold's are exchanged
        if(commBackward_[timeLevel]>=0)
        {
            MPI_Wait(&receive_request_, &receive_status_);
            for(int l=0; l<spaceSizes_[spaceLevel]; l++) xOld(l,0) = receiveBuffer_[timeLevel][l];
        }

        //			if(commForward_[timeLevel]>=0) MPI_Wait(&send_request_, &send_status_);


        //compute the residual
        vector<Vector> r(nTimeSteps);
        double rho0 = 1.0;
        {
            double localError = 0.0;

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localError)
#endif
            Vector r_acc;
            for(int ti=0; ti<nTimeSteps; ti++)
            {
                r[ti] = fl[ti];
                double sign = -1.0; A_[spaceLevel]->mult(x[ti], r[ti], timeLevel, ti, sign);
                if(ti!=0) A_[spaceLevel]->calcInitialVector(x[ti-1], r[ti], timeLevel, ti);
                if(commBackward_[timeLevel]>=0 && ti==0) A_[spaceLevel]->calcInitialVector(xOld, r[ti], timeLevel, ti);
                A_[spaceLevel]->accumulate(r[ti],r_acc,SpaceComm);
                localError += r_acc.col(0).dot(r[ti].col(0));
            }
            double errorTloacal = 0.0;
            MPI_Allreduce(&localError, &errorTloacal, 1, MPI_DOUBLE, MPI_SUM, SpaceComm);
            rho0 = 0.0;
            MPI_Allreduce(&errorTloacal, &rho0, 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);
            rho0 = sqrt(rho0);
        }

        if(rho0 < tol)  return 0;

        vector<vector<double> > H(maxit+1);
        vector<vector<Vector> > V(maxit+1);
        vector<vector<Vector> > Vp(maxit+1);
        Vector VpOld(N,1);
        vector<double> c(maxit+1);
        vector<double> s(maxit+1);
        vector<double> gamma(maxit+1);

        V[0].resize(nTimeSteps);
        for(int ti=0; ti<nTimeSteps; ti++) V[0][ti] = (1.0/rho0)*r[ti];
        gamma[0] = rho0;

        vector<Vector> w(nTimeSteps); for(int ti=0; ti<nTimeSteps; ti++) { w[ti].setZero(N,1);}
        double beta = 0.0;
        double tmp = 0.0;
        double relErrOld = 1.0;

        int j = 0;
        for(j=0; j < maxit; j++)
        {
            H[j].resize(j+2);
            Vp[j].resize(nTimeSteps);

            //apply preconditioner to V
            if(true)
            {
                for(int ti=0; ti<nTimeSteps; ti++)
                {
                    w[ti] = V[j][ti];
                    vec.u[timeLevel][ti].setZero();
                    xold_[timeLevel].setZero();

                }
                if(nTimeSteps>0) SpaceTimeMGMSchritt(vec, w, timeLevel);
                for(int ti=0; ti<nTimeSteps; ti++) Vp[j][ti] = vec.u[timeLevel][ti];

            }else
            {
                for(int ti=0; ti<nTimeSteps; ti++)
                {
                    w[ti] = V[j][ti];
                    vec.u[timeLevel][ti].setZero();
                    xold_[timeLevel].setZero();

                }
                SpaceTimeSmoother(vec, w, timeLevel, timeSm1_);

                for(int ti=0; ti<nTimeSteps; ti++) Vp[j][ti] = vec.u[timeLevel][ti];
            }
            if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], N, MPI_DOUBLE, commBackward_[timeLevel], 888, MPICommVars_.Comm, &receive_request_);

            //communicate the the values between the processors
            if(commForward_[timeLevel]>=0)
            {
                for(int l=0; l<N; l++) sendBuffer_[timeLevel][l] = Vp[j][nTimeSteps-1](l,0);
                MPI_Isend(sendBuffer_[timeLevel], N, MPI_DOUBLE, commForward_[timeLevel], 888, MPICommVars_.Comm, &send_request_);
            }


            //wait until the xold's are exchanged
            if(commBackward_[timeLevel]>=0)
            {
                MPI_Wait(&receive_request_, &receive_status_);
                for(int l=0; l<N; l++) VpOld(l,0) = receiveBuffer_[timeLevel][l];
            }

            //				if(commForward_[timeLevel]>=0) MPI_Wait(&send_request_, &send_status_);


            //compute w = A * Vp[j]
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=0; ti<nTimeSteps; ti++)
            {
                w[ti].setZero();
                double sign = -1.0; A_[spaceLevel]->mult(Vp[j][ti], w[ti], timeLevel, ti, sign);
                if(ti!=0) A_[spaceLevel]->calcInitialVector(Vp[j][ti-1], w[ti], timeLevel, ti);
                if(commBackward_[timeLevel]>=0 && ti==0) A_[spaceLevel]->calcInitialVector(VpOld, w[ti], timeLevel, ti);

                w[ti] *= -1.0;
            }


            Vector w_acc;
            for(int i=0; i<=j; i++)
            {
                //spord(V[i],w)
                double localH = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localH)
#endif
                for(int ti=0; ti<nTimeSteps; ti++)
                {
                    A_[spaceLevel]->accumulate(w[ti], w_acc, SpaceComm);
                    localH += V[i][ti].col(0).dot(w_acc.col(0));
                }
                double TloacalH = 0.0;
                MPI_Allreduce(&localH, &TloacalH, 1, MPI_DOUBLE, MPI_SUM, SpaceComm);
                H[j][i] = 0.0;
                MPI_Allreduce(&TloacalH, &H[j][i], 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
                for(int ti=0; ti<nTimeSteps; ti++) w[ti] += -H[j][i]*V[i][ti];
            }

            //sqrt(sprod(w,w))
            double localHH = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localHH)
#endif

            for(int ti=0; ti<nTimeSteps; ti++)
            {
                A_[spaceLevel]->accumulate(w[ti], w_acc, SpaceComm);
                localHH += w[ti].col(0).dot(w_acc.col(0));
            }
            double TloacalHH = 0.0;
            MPI_Allreduce(&localHH, &TloacalHH, 1, MPI_DOUBLE, MPI_SUM, SpaceComm);
            H[j][j+1] = 0.0;
            MPI_Allreduce(&TloacalHH, &H[j][j+1], 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);
            H[j][j+1] = sqrt(H[j][j+1]);


            for(int i=0; i<j; i++)
            {
                tmp = H[j][i];
                H[j][i] = c[i+1]*tmp + s[i+1]*H[j][i+1];
                H[j][i+1] = s[i+1]*tmp - c[i+1]*H[j][i+1];
            }
            beta = sqrt(H[j][j]*H[j][j] + H[j][j+1]*H[j][j+1]);

            s[j+1] = H[j][j+1]/beta;
            c[j+1] = H[j][j]/beta;
            H[j][j] = beta;
            gamma[j+1] = s[j+1]*gamma[j];
            gamma[j] = c[j+1]*gamma[j];

            if(MPICommVars_.my_rank==0) if(output)
            {
                double relErr = fabs(gamma[j+1])/rho0;
                gsInfo << "\t " << j+1 << ": " << relErr << "                          " << std::flush;
                if(int(log(relErrOld)/log(10)) > int(log(relErr)/log(10))) gsInfo << std::endl;
                relErrOld = relErr;
            }


            if(fabs(gamma[j+1]) < tol*rho0 || j == maxit-1) break;

            V[j+1].resize(nTimeSteps);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=0; ti<nTimeSteps; ti++)
            {
                w[ti] /= H[j][j+1];
                V[j+1][ti] = w[ti];
            }

        }

        /*
        for(int p=0; p<gsMpi::worldSize();++p)
        {
            if(p==gsMpi::worldRank())
            {

                gsInfo<<"rank: "<<gsMpi::worldRank()<<"- H: \n";
                for(int i=0; i<H.size();++i){
                    for(int j=0;j<H[i].size();++j)
                        gsInfo<<H[i][j]<<", ";
                    gsInfo<<"\n";
                }

            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        */


        if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], N, MPI_DOUBLE, commBackward_[timeLevel], 999, MPICommVars_.Comm, &receive_request_);

        vector<double> alpha(j+1);
        for(int i=j; i>=0; i--)
        {
            alpha[i] = gamma[i];
            for(int k=i+1; k<=j; k++) alpha[i]-= H[k][i]*alpha[k];
            alpha[i]/=H[i][i];
        }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for(int ti=0; ti<nTimeSteps; ti++)
        {
            for(int i=0;i<=j;i++) x[ti] +=alpha[i]* Vp[i][ti];
            vec.u[timeLevel][ti] = x[ti];
        }


        //communicate the the values between the processors
        if(commForward_[timeLevel]>=0)
        {
            for(int l=0; l<N; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][nTimeSteps-1](l,0);
            MPI_Isend(sendBuffer_[timeLevel], N, MPI_DOUBLE, commForward_[timeLevel], 999, MPICommVars_.Comm, &send_request_);
        }
        //wait until the xold's are exchanged
        if(commBackward_[timeLevel]>=0)
        {
            MPI_Wait(&receive_request_, &receive_status_);
            for(int l=0; l<N; l++) xold_[timeLevel](l,0) = receiveBuffer_[timeLevel][l];
        }
        if(commForward_[timeLevel]>=0) MPI_Wait(&send_request_, &send_status_);

        return j+1;
    }

    unsigned int PGMRESSolve(STSolVector &vec, const vector<Vector> &fl, int maxit=100, double tol=1e-6, bool output=true, bool calcResiduum=true, bool setZero = true)
    {
        if(maxit==0) return 0;

        int timeLevel = NTimeLevels_-1;

        int nTimeSteps = MytimeSteps_[timeLevel];
        int spaceLevel = spaceLevels_[timeLevel];
        int N = spaceSizes_[spaceLevel];

        //initial solution
        if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], N, MPI_DOUBLE, commBackward_[timeLevel], 789, MPICommVars_.Comm, &receive_request_);
        vector<Vector> x(nTimeSteps);
        Vector xOld(N,1); xOld.setZero();

        //communicate the the values between the processors
        if(commForward_[timeLevel]>=0)
        {
            if(setZero) for(int l=0; l<N; l++) sendBuffer_[timeLevel][l] = 0.0;
            else for(int l=0; l<N; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][nTimeSteps-1](l,0);
            MPI_Isend(sendBuffer_[timeLevel], N, MPI_DOUBLE, commForward_[timeLevel], 789, MPICommVars_.Comm, &send_request_);
        }

        if(setZero) for(int ti=0; ti<nTimeSteps; ti++) { x[ti].setZero(N,1); }
        else for(int ti=0; ti<nTimeSteps; ti++) x[ti] = vec.u[timeLevel][ti];

        //wait until the xold's are exchanged
        if(commBackward_[timeLevel]>=0)
        {
            MPI_Wait(&receive_request_, &receive_status_);
            for(int l=0; l<spaceSizes_[spaceLevel]; l++) xOld(l,0) = receiveBuffer_[timeLevel][l];
        }

        //			if(commForward_[timeLevel]>=0) MPI_Wait(&send_request_, &send_status_);


        //compute the residual
        vector<Vector> r(nTimeSteps);
        double rho0 = 1.0;
        {
            double localError = 0.0;

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localError)
#endif
            for(int ti=0; ti<nTimeSteps; ti++)
            {
                r[ti] = fl[ti];
                double sign = -1.0; A_[spaceLevel]->mult(x[ti], r[ti], timeLevel, ti, sign);
                if(ti!=0) A_[spaceLevel]->calcInitialVector(x[ti-1], r[ti], timeLevel, ti);
                if(commBackward_[timeLevel]>=0 && ti==0) A_[spaceLevel]->calcInitialVector(xOld, r[ti], timeLevel, ti);

                localError += r[ti].col(0).dot(r[ti].col(0));
            }
            rho0 = 0.0;
            MPI_Allreduce(&localError, &rho0, 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);
            rho0 = sqrt(rho0);
        }

        if(rho0 < tol)  return 0;

        vector<vector<double> > H(maxit+1);
        vector<vector<Vector> > V(maxit+1);
        vector<vector<Vector> > Vp(maxit+1);
        Vector VpOld(N,1);
        vector<double> c(maxit+1);
        vector<double> s(maxit+1);
        vector<double> gamma(maxit+1);

        V[0].resize(nTimeSteps);
        for(int ti=0; ti<nTimeSteps; ti++) V[0][ti] = (1.0/rho0)*r[ti];
        gamma[0] = rho0;

        vector<Vector> w(nTimeSteps); for(int ti=0; ti<nTimeSteps; ti++) { w[ti].setZero(N,1);}
        double beta = 0.0;
        double tmp = 0.0;
        double relErrOld = 1.0;

        int j = 0;
        for(j=0; j < maxit; j++)
        {
            H[j].resize(j+2);
            Vp[j].resize(nTimeSteps);

            //apply preconditioner to V
            if(true)
            {
                for(int ti=0; ti<nTimeSteps; ti++)
                {
                    w[ti] = V[j][ti];
                    vec.u[timeLevel][ti].setZero();
                    xold_[timeLevel].setZero();

                }
                if(nTimeSteps>0) SpaceTimeMGMSchritt(vec, w, timeLevel);
                for(int ti=0; ti<nTimeSteps; ti++) Vp[j][ti] = vec.u[timeLevel][ti];

            }else
            {
                for(int ti=0; ti<nTimeSteps; ti++)
                {
                    w[ti] = V[j][ti];
                    vec.u[timeLevel][ti].setZero();
                    xold_[timeLevel].setZero();

                }
                SpaceTimeSmoother(vec, w, timeLevel, timeSm1_);

                for(int ti=0; ti<nTimeSteps; ti++) Vp[j][ti] = vec.u[timeLevel][ti];
            }
            if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], N, MPI_DOUBLE, commBackward_[timeLevel], 888, MPICommVars_.Comm, &receive_request_);

            //communicate the the values between the processors
            if(commForward_[timeLevel]>=0)
            {
                for(int l=0; l<N; l++) sendBuffer_[timeLevel][l] = Vp[j][nTimeSteps-1](l,0);
                MPI_Isend(sendBuffer_[timeLevel], N, MPI_DOUBLE, commForward_[timeLevel], 888, MPICommVars_.Comm, &send_request_);
            }


            //wait until the xold's are exchanged
            if(commBackward_[timeLevel]>=0)
            {
                MPI_Wait(&receive_request_, &receive_status_);
                for(int l=0; l<N; l++) VpOld(l,0) = receiveBuffer_[timeLevel][l];
            }

            //				if(commForward_[timeLevel]>=0) MPI_Wait(&send_request_, &send_status_);


            //compute w = A * Vp[j]
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=0; ti<nTimeSteps; ti++)
            {
                w[ti].setZero();
                double sign = -1.0; A_[spaceLevel]->mult(Vp[j][ti], w[ti], timeLevel, ti, sign);
                if(ti!=0) A_[spaceLevel]->calcInitialVector(Vp[j][ti-1], w[ti], timeLevel, ti);
                if(commBackward_[timeLevel]>=0 && ti==0) A_[spaceLevel]->calcInitialVector(VpOld, w[ti], timeLevel, ti);

                w[ti] *= -1.0;
            }


            for(int i=0; i<=j; i++)
            {
                //spord(V[i],w)
                double localH = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localH)
#endif
                for(int ti=0; ti<nTimeSteps; ti++)
                {
                    localH += V[i][ti].col(0).dot(w[ti].col(0));
                }
                H[j][i] = 0.0;
                MPI_Allreduce(&localH, &H[j][i], 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
                for(int ti=0; ti<nTimeSteps; ti++) w[ti] += -H[j][i]*V[i][ti];
            }

            //sqrt(sprod(w,w))
            double localHH = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+: localHH)
#endif
            for(int ti=0; ti<nTimeSteps; ti++)
            {
                localHH += w[ti].col(0).dot(w[ti].col(0));
            }
            H[j][j+1] = 0.0;
            MPI_Allreduce(&localHH, &H[j][j+1], 1, MPI_DOUBLE, MPI_SUM, MPICommVars_.Comm);
            H[j][j+1] = sqrt(H[j][j+1]);


            for(int i=0; i<j; i++)
            {
                tmp = H[j][i];
                H[j][i] = c[i+1]*tmp + s[i+1]*H[j][i+1];
                H[j][i+1] = s[i+1]*tmp - c[i+1]*H[j][i+1];
            }
            beta = sqrt(H[j][j]*H[j][j] + H[j][j+1]*H[j][j+1]);

            s[j+1] = H[j][j+1]/beta;
            c[j+1] = H[j][j]/beta;
            H[j][j] = beta;
            gamma[j+1] = s[j+1]*gamma[j];
            gamma[j] = c[j+1]*gamma[j];

            if(MPICommVars_.my_rank==0) if(output)
            {
                double relErr = fabs(gamma[j+1])/rho0;
                gsInfo << "\t " << j+1 << ": " << relErr << "                          " << std::flush;
                if(int(log(relErrOld)/log(10)) > int(log(relErr)/log(10))) gsInfo << std::endl;
                relErrOld = relErr;
            }


            if(fabs(gamma[j+1]) < tol*rho0 || j == maxit-1) break;

            V[j+1].resize(nTimeSteps);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int ti=0; ti<nTimeSteps; ti++)
            {
                w[ti] /= H[j][j+1];
                V[j+1][ti] = w[ti];
            }

        }

        if(commBackward_[timeLevel]>=0) MPI_Irecv(receiveBuffer_[timeLevel], N, MPI_DOUBLE, commBackward_[timeLevel], 999, MPICommVars_.Comm, &receive_request_);

        vector<double> alpha(j+1);
        for(int i=j; i>=0; i--)
        {
            alpha[i] = gamma[i];
            for(int k=i+1; k<=j; k++) alpha[i]-= H[k][i]*alpha[k];
            alpha[i]/=H[i][i];
        }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for(int ti=0; ti<nTimeSteps; ti++)
        {
            for(int i=0;i<=j;i++) x[ti] +=alpha[i]* Vp[i][ti];
            vec.u[timeLevel][ti] = x[ti];
        }


        //communicate the the values between the processors
        if(commForward_[timeLevel]>=0)
        {
            for(int l=0; l<N; l++) sendBuffer_[timeLevel][l] = vec.u[timeLevel][nTimeSteps-1](l,0);
            MPI_Isend(sendBuffer_[timeLevel], N, MPI_DOUBLE, commForward_[timeLevel], 999, MPICommVars_.Comm, &send_request_);
        }
        //wait until the xold's are exchanged
        if(commBackward_[timeLevel]>=0)
        {
            MPI_Wait(&receive_request_, &receive_status_);
            for(int l=0; l<N; l++) xold_[timeLevel](l,0) = receiveBuffer_[timeLevel][l];
        }
        if(commForward_[timeLevel]>=0) MPI_Wait(&send_request_, &send_status_);

        return j+1;
    }

    int TimeMGMTriDiag_MPI(STSolVector &vec, const vector<Vector> &fl, int maxit=100, double tol=1e-6, bool output=true, bool calcResiduum=true)
    {
        if(!useTriDiag_) gsInfo << "WARNING: you should enable the flag useTriDiag_!!" << std::endl;

        int timeLevel = NTimeLevels_-1;

        double error0 = 1.0;

        double rho = 1.0;
        double errorOld = error0;
        double error = 1.0;

        vector<double> rhos(maxit);
        double maxRho = 0.0;

        double averRho = 1.0;

        int i = 0;
        for(i=0; i<maxit; i++)
        {
            SpaceTimeMGMSchrittTriDiag(vec, fl, timeLevel);

            if(calcResiduum || output)
            {
                error = ComputeSpaceTimeResiduumTriDiag(vec, fl, timeLevel);
                if(i==0)
                {
                    error0 = error;
                    errorOld = error;
                    if(error0<tol)
                    {
                        if(MPICommVars_.my_rank==0) if(output) gsInfo << "\t\t\t" << i+1 << ": " << error << "  -->  " << rho  << " (" <<  pow(error/error0,1.0/double(i)) << ")" << std::endl;
                        break;
                    }
                }

                rho = error/errorOld;
                errorOld = error;

                if(MPICommVars_.my_rank==0) if(output) gsInfo << "\t\t\t" << i+1 << ": " << error/error0 << "  -->  " << rho  << " (" <<  pow(error/error0,1.0/double(i)) << ")" << std::endl;

                if(error>0) averRho = pow(error/error0,1.0/double(i));

                rhos[i] = rho;
                if(maxRho<rho && i>0) maxRho = rho;

                if(error/error0 < tol || error/error0>1e6) break;
            }
        }

        if(MPICommVars_.my_rank==0) if(output)
            gsInfo << "\t maximal contraction factor: " << maxRho << std::endl <<
                      "\t average contraction factor: " << pow(error/error0,1.0/double(i)) << std::endl;



        //write something to a file
        /*
        if(MPICommVars_.my_rank==0 && false)
        {

            time_t rawtime;
            struct tm * timeinfo;
            char buffer[80];

            time (&rawtime);
            timeinfo = localtime(&rawtime);

            strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
            std::string strTime(buffer);

            double mue = A_[spaceLevels_[timeLevel]]->getCoerseningMeasure(timeLevel);
            double rho = A_[spaceLevels_[timeLevel]]->getDesignParameter();
            int rhoExp = int(round(log(rho)/log(10.0)));


            int sm = timeSm1_;
            int poly = A_[spaceLevels_[timeLevel]]->getOrder();

            std::string FEnd; std::stringstream out; out << "_rho" << rhoExp << "_p" << poly << "_s" << sm << "_Lx" << spaceLevels_[timeLevel] << "_Lt" << timeLevel; FEnd = out.str();

            char const* home = getenv("HOME");
            std::string homefolder = std::string(home)+std::string("/");
            std::string folder = homefolder+"Simulation/Measures/MGConergence/";

            std::string filenameMax = folder + "maxRho" + FEnd;
            std::string filenameAverage = folder + "averageRho" + FEnd;

            std::ofstream out_fileMax(filenameMax.c_str(), std::fstream::app);
            out_fileMax << strTime << "\t" << mue << "\t" << maxRho << std::endl;
            out_fileMax.close();

            std::ofstream out_fileAver(filenameAverage.c_str(), std::fstream::app);
            out_fileAver << strTime << "\t" << mue << "\t" << averRho << std::endl;
            out_fileAver.close();
        }
*/

        return i+1;
    }

    void buildGlobalSolution(const std::vector<gsMatrix<real_t> >& sol, const gsVector<index_t>& sizes, Vector& result)
    {
        int s=sizes.sum();
        result.setZero(s,1);
        gsVector<index_t> col(1);
        col<<1;
        gsMatrix<real_t>::BlockView view = result.blockView(sizes,col);

        //gsInfo<<"result rows: "<< s<<"\n";

        for(int i=0; i<MytimeSteps_[NTimeLevels_-1];++i)
        {
            //gsInfo<<"sol.u[NTimeLevels_-1][i]: "<< sol.u[NTimeLevels_-1][i].rows()<<" view("<<MPICommVars_.my_rank*MytimeSteps_[NTimeLevels_-1]+i<<",0)"<<view(MPICommVars_.my_rank*MytimeSteps_[NTimeLevels_-1]+i,0).rows()<<"\n";
            view(MPICommVars_.my_rank*MytimeSteps_[NTimeLevels_-1]+i,0) = sol[i];

        }

        gsMpiComm comm(MPICommVars_.Comm);
        comm.sum(result.data(),result.size(),0);

        //     MPI_Reduce(result.data(),result.data(),result.rows(),MPI_DOUBLE,MPI_SUM,0,MPICommVars_.Comm);

    }

};




#endif //GISMO_WITH_MPI


} // namespace gismo
