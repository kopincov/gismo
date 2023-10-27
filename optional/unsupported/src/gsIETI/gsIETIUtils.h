/** @file gsAssemblerUtils.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBoundary.h>
#include <gsCore/gsBoxTopology.h>
#include <gsIO/gsOptionList.h>
#include <iomanip>

namespace gismo
{


/// \brief Method for selecting primal degrees of freedom
struct primalDofMethod
{
    primalDofMethod() : strat(B) {}
    enum method{

        /// \brief  only vertices
        vertices = 1<<0,

        /// \brief  only edges
        edges = 1<<1,

        /// \brief  only faces
        faces = 1<<2

    };
    enum strategy{
        /// \brief strategy A: only vertices
        A = vertices,

        /// \brief strategy B: all methods
        B = vertices | edges | faces,

        /// \brief strategy C: vertices and edges
        C = vertices | edges,

        /// \brief strategy D: only edges
        D= edges,


        // for debugging purposes
        X= faces,

        Y = vertices |faces,

        Z = edges |faces
    };

    /// \brief checks if the method m is contained in the given strategy
    bool contains(method m) const
    {
        if((strat & m)== m)
            return true;
        else
            return false;
    }

    static int getMethodIndex(method m)
    {
        switch(m){
        case vertices:
            return 0;
            break;
        case edges:
            return 1;
            break;
        case faces:
            return 2;
            break;
        default:
            return -1;
        };
    }

    static primalDofMethod::strategy chooseMethod(std::string str)
    {
        if(str=="A")
            return primalDofMethod::A;
        else if(str == "B")
            return primalDofMethod::B;
        else if(str == "C")
            return primalDofMethod::C;
        else if(str == "D")
            return primalDofMethod::D;
        else if(str=="X")
            return primalDofMethod::X;
        else if(str == "Y")
            return primalDofMethod::Y;
        else if(str == "Z")
            return primalDofMethod::Z;

        GISMO_ERROR("Useage of a valid primal dof method, A,B,C,D,X,Y,Z");
    }

    strategy strat;
};

/// \brief  Method for choosing the scaling of the jump operator
struct IETIPrecondScaling
{
    enum strategy{

        /// \brief  no scaling
        none =0,

        /// \brief  scaling by the multiplicity of the dof
        multiplicity = 1,

        /// \brief  scaling by the diffusion coefficient
        coefficient = 2,

        /// \brief scaling by the diagonal entry of the stiffness matrix
        stiffness = 3,

        /// \brief scaling by ONE diagonal entry of the stiffness matrix, the same for all dofs on a patch
        stiffnessModified = 4
    };

    static strategy chooseMethod(std::string name)
    {
        if(name == "none" ) return none;
        if(name == "mult" ) return multiplicity;
        if(name == "coeff" ) return coefficient;
        if(name == "stiff" ) return stiffness;
        if(name == "stiffM") return stiffnessModified;
        GISMO_ERROR("No valid scaling method chosen.");
        return none;
    }
};

struct IETILocalSolver
{
    enum solvers
    {
        Direct = 0,
        FastDiagonalization,
        Multigrid,
    };

    static IETILocalSolver::solvers chooseSolver(std::string str)
    {
        if(str=="D")
            return IETILocalSolver::Direct;
        else if(str == "BS")
            return IETILocalSolver::FastDiagonalization;
        else if(str == "MG")
            return IETILocalSolver::Multigrid;
        GISMO_ERROR("Useage of a valid solver, D (direct), BS (FastDiagonalization), MG (Multigrid)");
    }

};


struct gsInexactOptions
{
    gsInexactOptions( int max_iter, real_t tol): max_iter(max_iter) , tol(tol) {}

    int max_iter;
    real_t tol;
};

struct gsDummyOptions: gsInexactOptions
{
    gsDummyOptions(): gsInexactOptions(-1,-1) {}

    int max_iter;
    real_t tol;
};
struct gsIETIBSOptions : gsInexactOptions
{
    gsIETIBSOptions( int max_iter, real_t tol): gsInexactOptions(max_iter,tol), regularization(1.e-8) {}

    real_t regularization;
};

namespace Smoother {
enum type {
    Richardson,
    Jacobi,
    GaussSeidel,
    MassRichardson,
    MassRichardsonBoundaryCorrection,
    MassRichardsonSubspaceCorrection,
    MassRichardsonSubspaceCorrectionTrunc,
    MassRichardsonSubspaceCorrectionGS,
    Generic
};
}

struct gsIETIMGOptions : gsInexactOptions
{
    gsIETIMGOptions ( int max_iter, real_t tol, int numLevels) : gsInexactOptions(max_iter,tol),  numLevels(numLevels), numRefine(1),
        numPreSmooth(1), numPostSmooth(1), cycles(1), damping(1),coarseDamping(1), smoother(Smoother::MassRichardsonSubspaceCorrection) {}
    int numLevels;
    int numRefine;
    int numPreSmooth;
    int numPostSmooth;
    int cycles;
    real_t damping;
    real_t coarseDamping;
    Smoother::type smoother;
};

struct gsIETITimings
{
    gsIETITimings() {}

    gsIETITimings(size_t nPatches): KC_basis(nPatches),KC_subspace(nPatches), Kii_precond(nPatches),
        Kii_schur(nPatches), projection(nPatches), embedding(nPatches), embeddingT(nPatches),
        assemblingMass(nPatches), totalAssemble(nPatches),
        averageItNumKCBasis(nPatches),averageItNumKCSpace(nPatches), averageItNumKii(nPatches), numKCItBasis(nPatches),  numKCItSpace(nPatches),
        numKiiIt(nPatches)
    {
        KC_basis.setZero();
        KC_subspace= Kii_precond= Kii_schur= projection=embedding= embeddingT= assemblingMass= totalAssemble=  averageItNumKCBasis= averageItNumKCSpace= averageItNumKii=
                numKCItBasis= numKCItSpace= numKiiIt = KC_basis;
    }

    gsVector<real_t> KC_basis;
    gsVector<real_t> KC_subspace;
    gsVector<real_t> Kii_precond;
    gsVector<real_t> Kii_schur;
    gsVector<real_t> projection;
    gsVector<real_t> embedding;
    gsVector<real_t> embeddingT;
    gsVector<real_t> assemblingMass;
    gsVector<real_t> totalAssemble;

    gsVector<real_t> averageItNumKCBasis;
    gsVector<real_t> averageItNumKCSpace;
    gsVector<real_t> averageItNumKii;
    gsVector<real_t> numKCItBasis;
    gsVector<real_t> numKCItSpace;
    gsVector<real_t> numKiiIt;



    void print(std::ostream& out)
    {
        std::ios  state(NULL);
        state.copyfmt(std::cout);

        out<<std::fixed<<std::setprecision(4);
        out<<std::setw(12)<<"K-Basis "<<KC_basis.transpose()<<"\n";
        out<<std::setw(12)<<"K-Space "<<KC_subspace.transpose()<<"\n";
        out<<std::setw(12)<<"Kii-prec "<<Kii_precond.transpose()<<"\n";
        out<<std::setw(12)<<"Kii-ex "<<Kii_schur.transpose()<<"\n";
        out<<std::setw(12)<<"proj "<<projection.transpose()<<"\n";
        out<<std::setw(12)<<"embedding "<<embedding.transpose()<<"\n";
        out<<std::setw(12)<<"embeddingT "<<embeddingT.transpose()<<"\n";
        out<<std::setw(12)<<"assMass "<<assemblingMass.transpose()<<"\n";
        out<<std::setw(12)<<"ass Time "<<totalAssemble.transpose()<<"\n";
        out<<"\n";
        out<<std::setprecision(2);
        out<<std::setw(12)<<"MGIt-KC-B "<<averageItNumKCBasis.cwiseQuotient(numKCItBasis).transpose()<<"\n";
        out<<std::setw(12)<<"MGIt-KC-S "<<averageItNumKCSpace.cwiseQuotient(numKCItSpace).transpose()<<"\n";
        out<<std::setw(12)<<"MGIt-Kii "<<averageItNumKii.cwiseQuotient(numKiiIt).transpose()<<"\n";

        out<<std::setprecision(4)<<"\n\n";
        out<<std::setw(12)<<"K-Basis "<<KC_basis.transpose().sum()<<"\n";
        out<<std::setw(12)<<"K-Space "<<KC_subspace.transpose().sum()<<"\n";
        out<<std::setw(12)<<"Kii-prec "<<Kii_precond.transpose().sum()<<"\n";
        out<<std::setw(12)<<"Kii-ex "<<Kii_schur.transpose().sum()<<"\n";
        out<<std::setw(12)<<"proj "<<projection.transpose().sum()<<"\n";
        out<<std::setw(12)<<"embedding "<<embedding.transpose().sum()<<"\n";
        out<<std::setw(12)<<"embeddingT "<<embeddingT.transpose().sum()<<"\n";
        out<<std::setw(12)<<"assMass "<<assemblingMass.transpose().sum()<<"\n";
        out<<std::setw(12)<<"ass Time "<<totalAssemble.transpose().sum()<<"\n";
        out<<"\n";
        out<<std::setprecision(2);
        out<<std::setw(12)<<"MGIt-KC-B "<<averageItNumKCBasis.sum()/numKCItBasis.sum()<<"\n";
        out<<std::setw(12)<<"MGIt-KC-S "<<averageItNumKCSpace.sum()/numKCItSpace.sum()<<"\n";
        out<<std::setw(12)<<"MGIt-Kii "<<averageItNumKii.sum()/numKiiIt.sum()<<"\n";
        out<<"\n";
        std::cout.copyfmt(state);
    }
};

struct gsIETIOptions
{
    gsIETIOptions():
        //isMinimumEnergy(true), extraTimings(false), needRescaling(false),nonLinearMode(false), getRhsNorm(false),
        strat(), scal(IETIPrecondScaling::coefficient),
       // solvingSaddlePoint(false),
        KiiSolver(IETILocalSolver::Direct),
        KCSolver(IETILocalSolver::Direct),
        KrrSolver(IETILocalSolver::Direct),
        nSppHolder(1)
        //Kii_optionsBS(3,1.e-2),
        //Krr_optionsBS(3,1.e-2),
        //KC_optionsBS(3,1.e-2),
        //Kii_optionsMG(3,1.e-2,1),
        //Krr_optionsMG(3,1.e-2,1),
        //KC_optionsMG(3,1.e-2,1)
    {}
/*
    const gsInexactOptions & getInexactKii() const
    {
        switch(KiiSolver)
        {
        case IETILocalSolver::Multigrid:
            return Kii_optionsMG;
        case IETILocalSolver::FastDiagonalization:
            return Kii_optionsBS;
        case IETILocalSolver::Direct:
            return Kii_optionsD;
        }
        return Kii_optionsD;
    }

    const gsInexactOptions & getInexactKC() const
    {
        switch(KCSolver)
        {
        case IETILocalSolver::Multigrid:
            return KC_optionsMG;
        case IETILocalSolver::FastDiagonalization:
            return KC_optionsBS;
        case IETILocalSolver::Direct:
            return KC_optionsD;
        }
        return KC_optionsD;
    }
    */
    // bool isMinimumEnergy;

    // bool extraTimings;
    // bool needRescaling;

    // bool nonLinearMode;
   //  bool getRhsNorm;

    gsOptionList opt; // is set
    primalDofMethod strat; // is set
    IETIPrecondScaling::strategy scal;
   // bool solvingSaddlePoint;  // is set

    IETILocalSolver::solvers KiiSolver;  // is set
    IETILocalSolver::solvers KCSolver;  // is set
    IETILocalSolver::solvers KrrSolver; // is set
    int nSppHolder;
    /*
    gsIETIBSOptions Kii_optionsBS;
    gsIETIBSOptions Krr_optionsBS;
    gsIETIBSOptions KC_optionsBS;

    gsIETIMGOptions Kii_optionsMG;
    gsIETIMGOptions Krr_optionsMG;
    gsIETIMGOptions KC_optionsMG;

    gsDummyOptions Kii_optionsD;
    gsDummyOptions Krr_optionsD;
    gsDummyOptions KC_optionsD;
    */
};


//------------------------------------------------------------------------------------------------------------------
/**
 * @brief The Edge class
 *
 * This class represents an Edge of a 3D cube in the parameter domain.
 * It handles information about the location and probable eliminaition due to dirichlet conditions
 */
class Edge
{
public:

    /// \brief  if the edge lies (partially) on a dirichlet boundary
    unsigned isEliminated;

    /// \brief  number of dof on this edge
    unsigned number;

    /// \brief  number of the edge on the cube
    int edgeNumber;

    /// \brief  to which patchSide (==face) the edge is associated
    patchSide side3D;

    /// \brief  to which patchSide (==face) the edge is associated too
    patchSide side3D_2;

    /// \brief  one which Side of the \p side3D the edge is associate
    boxSide side2D;

    /// \brief  directions of the parameter domain, which have to be set fixed to get the direction of the edge
    int directions[2];

    /// \brief  which of the directions have to be set to zero or one
    bool parameters[2];

    /// \brief the component this edge belongs to
    int component;

    Edge(bool isEliminated_, unsigned number_, int edgeNumber_, patchSide side3D_, boxSide side2D_)
        : isEliminated(isEliminated_),  number(number_), edgeNumber(edgeNumber_) ,side3D(side3D_), side2D(side2D_)
    {
        directions[0] = side3D.direction();
        int d;
        bool p;
        switch(directions[0])
        {
        case 0:
            d = side2D.direction()+1;
            break;
        case 1:
            d = side2D.direction()*2;
            break;
        case 2:
            d= side2D.direction();
            break;
        default:
            d=0;
            gsWarn<<"the dimension of the object seems to high!"<<std::endl;
            break;

        }

        p=side2D.parameter(); //this really works
        directions[1] = d;
        parameters[0] = side3D.parameter();
        parameters[1] = p;
        side3D_2=patchSide(side3D.patch,boxSide(d,p));
    }
    /**
     * @brief operator == checks for equality, based on the edge number
     * @param e the other edge
     * @return true if equal
     */
    bool operator==(const Edge& e) const {return edgeNumber == e.edgeNumber;}

    /**
     * @brief operator != checks for not equality, based on the edge number
     * @param e the other edge
     * @return true if not equal
     */
    bool operator!=(const Edge& e) const {return !(*this == e);}

    /**
     * @brief operator < checks for smaller, based on the edge number
     * @param e the other edge
     * @return  true if smaller
     */
    bool operator<(const Edge& e) const {return edgeNumber < e.edgeNumber;}

};




/**
 * @brief The gsIETIInfo struct stores information about the IETI method, e.g. dofs on the boundary, primal dofs...
 */
struct gsIETIInfo
{
    size_t numberPatches;
    int nRhs;

    int dim;
    size_t cDim;

    unsigned dofTotal;
    unsigned dofTotalP;
    unsigned dofTotalI;
    unsigned dofTotalB;

    unsigned origSystemSize;

    std::vector<int> dofTotalPtype;
    unsigned dofInterface;
    unsigned lagrangeMult;
    std::vector<index_t> dofs; //<- size of the local vectors dofsB or dofsB+dofsI
    std::vector<index_t> dofsB;
    std::vector<index_t> dofsI;
    std::vector<index_t> dofsP;
    std::vector<index_t> dofsR;
    std::vector<std::vector<int> > dofsPtype;

    //--------------------------------------------------------------
    void print() const{
        std::cout<< "number of Patches:  "<< numberPatches<<std::endl;
        std::cout<< "number of rhs:  "<< nRhs<<std::endl;
        std::cout<< "dimension of the object:  "<< dim<<std::endl;
        std::cout<< "Total number of free Dofs:  "<< dofTotal<<std::endl;
        std::cout<< "Total number of primal Dofs:  "<< dofTotalP<<std::endl;
        std::cout<< "global number of interface dofs: "<<dofInterface<<std::endl;
        std::cout<< "global number of lagrange multipiers: "<<lagrangeMult<<std::endl;
        std::cout<< "original system size: "<<origSystemSize<<std::endl;

        std::cout<< "Free Dofs in Patch i:  ";
        for(size_t i=0; i<dofs.size();i++)
            std::cout<<"i= "<<i<<":   " <<dofs[i]<<", ";
        std::cout<< std::endl;

        std::cout<< "Boundary Dofs in Patch i:  ";
        for(size_t i=0; i<dofsB.size();i++)
            std::cout<<"i= "<<i<<":   " <<dofsB[i]<<", ";
        std::cout<< std::endl;

        std::cout<< "Interior Dofs in Patch i:  ";
        for(size_t i=0; i<dofsI.size();i++)
            std::cout<<"i= "<<i<<":   " <<dofsI[i]<<", ";
        std::cout<< std::endl;

        std::cout<< "Primal Dofs in Patch i:  ";
        for(size_t i=0; i<dofsP.size();i++)
            std::cout<<"i= "<<i<<":   " <<dofsP[i]<<", ";
        std::cout<< std::endl;



    }
};
}

// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2012 Desire NUENTSA WAKAM <desire.nuentsa_wakam@inria.fr
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

/**
  * \ingroup IterativeSolvers_Module
  * \brief iterative scaling algorithm to equilibrate rows and column norms in matrices
  *
  * This class can be used as a preprocessing tool to accelerate the convergence of iterative methods
  *
  * This feature is  useful to limit the pivoting amount during LU/ILU factorization
  * The  scaling strategy as presented here preserves the symmetry of the problem
  * NOTE It is assumed that the matrix does not have empty row or column,
  *
  * Example with key steps
  * \code
  * VectorXd x(n), b(n);
  * SparseMatrix<double> A;
  * // fill A and b;
  * IterScaling<SparseMatrix<double> > scal;
  * // Compute the left and right scaling vectors. The matrix is equilibrated at output
  * scal.computeRef(A);
  * // Scale the right hand side
  * b = scal.LeftScaling().cwiseProduct(b);
  * // Now, solve the equilibrated linear system with any available solver
  *
  * // Scale back the computed solution
  * x = scal.RightScaling().cwiseProduct(x);
  * \endcode
  *
  * \tparam _MatrixType the type of the matrix. It should be a real square sparsematrix
  *
  * References : D. Ruiz and B. Ucar, A Symmetry Preserving Algorithm for Matrix Scaling, INRIA Research report RR-7552
  *
  * \sa \ref IncompleteLUT
  */
namespace Eigen {
using std::abs;

template<typename _MatrixType>
class IterScaling
{
public:
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
    typedef Matrix<Scalar, Dynamic,1> Vector;

public:
    IterScaling() { init(); }

    IterScaling(const MatrixType& matrix)
    {
        init();
        compute(matrix);
    }

    ~IterScaling() { }

    /**
     * Compute the left and right diagonal matrices to scale the input matrix @p mat
     *
     * FIXME This algorithm will be modified such that the diagonal elements are permuted on the diagonal.
     *
     * \sa LeftScaling() RightScaling()
     */
    void compute (MatrixType& mat)
    {
        int m = mat.rows();
        int n = mat.cols();

        eigen_assert((m == n) && "Please give a non - empty matrix");
        m_left.resize(m);
        m_right.resize(n);
        m_left.setOnes();
        m_right.setOnes();
        m_matrix = &mat;
        if(m==0)
        {
            m_isInitialized = true;
            return;
        }

        Vector Dr, Dc, DrRes, DcRes; // Temporary Left and right scaling vectors
        Dr.resize(m); Dc.resize(n);
        DrRes.resize(m); DcRes.resize(n);
        Scalar EpsRow = 1.0, EpsCol = 1.0;
        int its = 0;
        do
        { // Iterate until the infinite norm of each row and column is approximately 1
            // Get the maximum value in each row and column
            Dr.setZero(); Dc.setZero();
            for (int k=0; k<m_matrix->outerSize(); ++k)
            {
                for (typename MatrixType::InnerIterator it(*m_matrix, k); it; ++it)
                {
                    if ( Dr(it.row()) < abs(it.value()) )
                        Dr(it.row()) = abs(it.value());

                    if ( Dc(it.col()) < abs(it.value()) )
                        Dc(it.col()) = abs(it.value());
                }
            }
            for (int i = 0; i < m; ++i)
            {
                Dr(i) = gismo::math::sqrt(Dr(i));
                Dc(i) = gismo::math::sqrt(Dc(i));
            }
            // Save the scaling factors
            for (int i = 0; i < m; ++i)
            {
                m_left(i) /= Dr(i);
                m_right(i) /= Dc(i);
            }
            // Scale the rows and the columns of the matrix
            DrRes.setZero(); DcRes.setZero();
            for (int k=0; k<m_matrix->outerSize(); ++k)
            {
                for (typename MatrixType::InnerIterator it(*m_matrix, k); it; ++it)
                {
                    it.valueRef() = it.value()/( Dr(it.row()) * Dc(it.col()) );
                    // Accumulate the norms of the row and column vectors
                    if ( DrRes(it.row()) < abs(it.value()) )
                        DrRes(it.row()) = abs(it.value());

                    if ( DcRes(it.col()) < abs(it.value()) )
                        DcRes(it.col()) = abs(it.value());
                }
            }
            DrRes.array() = (1-DrRes.array()).abs();
            EpsRow = DrRes.maxCoeff();
            DcRes.array() = (1-DcRes.array()).abs();
            EpsCol = DcRes.maxCoeff();
            its++;

            // std::cout<<EpsRow<<std::endl<<std::endl;
        }while ( (EpsRow >m_tol || EpsCol > m_tol) && (its < m_maxits) );
        m_isInitialized = true;
    }
    /** Compute the left and right vectors to scale the vectors
     * the input matrix is scaled with the computed vectors at output
     *
     * \sa compute()
     */
    void computeRef (MatrixType& mat)
    {
        compute (mat);
        //mat = m_matrix;
    }
    /** Get the vector to scale the rows of the matrix
     */
    const Vector& LeftScaling() const
    {
        return m_left;
    }

    /** Get the vector to scale the columns of the matrix
     */
    const Vector& RightScaling() const
    {
        return m_right;
    }

    /** Set the tolerance for the convergence of the iterative scaling algorithm
     */
    void setTolerance(double tol)
    {
        m_tol = tol;
    }

protected:

    void init()
    {
        m_tol = 1e-5;
        m_maxits = 5;
        m_isInitialized = false;
    }

    MatrixType* m_matrix;
    mutable ComputationInfo m_info;
    bool m_isInitialized;
    Vector m_left; // Left scaling vector
    Vector m_right; // m_right scaling vector
    double m_tol;
    int m_maxits; // Maximum number of iterations allowed
};
}

