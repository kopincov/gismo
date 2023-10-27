
#pragma once

#include <iomanip>
#include <string>


namespace gismo
{

namespace lsdyna
{

    // typedefs
    
    typedef int lsInt;
    typedef double lsFloat;

    
namespace write
{

    // Purpose: Activate implicit eigenvalue analysis and define associated
    // input parameters (see also *CONTROL_IMPLICIT_GENERAL).
    //
    // Parameters:
    // Number of eigenvalues to extract. This must be specified. The other
    // parameters below are optional.
    void CONTROL_IMPLICIT_EIGENVALUE(std::ostream& out,
                                     const lsInt numEigVals)
    {
        out << "*CONTROL_IMPLICIT_EIGENVALUE\n";
        out << "       " << numEigVals << "\n";
    }
                                
    
    // Purpose: Use the consistent mass matrix in implicit dynamics and
    // eigenvalue solutions.
    //
    // Parameters:
    // Consistent mass matrix flag
    // EQ.0: Use the standard lumped mass formulation (DEFAULT)
    // EQ.1: Use the consistent mass matrix.
    void CONTROL_IMPLICIT_CONSISTENT_MASS(std::ostream& out,
                                         const lsInt flag)
    {
        out << "*CONTROL_IMPLICIT_CONSISTENT_MASS\n";
        out << flag << "\n";
    }

    // Purpose: Activate implicit analysis and define associated control
    // parameters. This keyword is required for all implicit analyses.
    //
    // Parameters:
    // imflag Implicit/Explicit analysis type flag
    //        EQ.0: explicit analysis
    //        EQ.1: implicit analysis
    //        ...
    // dto    Initial time step size for implicit analysis
    void CONTROL_IMPLICIT_GENERAL(std::ostream& out,
                                  const lsInt imflag,
                                  const lsFloat dt0)
    {
        out << "*CONTROL_IMPLICIT_GENERAL\n";
        out << "        " << imflag << " ";
        out << std::setprecision(3) << std::scientific;
        out << dt0 << "\n";
    }

    // Purpose: These optional cards apply to implicit calculations.
    //
    // Parameters:
    // lsolvr Linear equation solver method
    //        EQ.4: SMP parallel multi-frontal sparse solver (default).
    //        EQ.5: SMP parallel multi-frontal sparse solver, double precision
    //        ...
    // lprint Linear solver print flag controls screen and message file output
    //        EQ.0: no printing
    //        EQ.1: output summary statistics on memory, cpu requirements
    //        EQ.2: more statistics
    void CONTROL_IMPLICIT_SOLVER(std::ostream& out,
                                 const lsInt lsolvr,
                                 const lsInt lprint)
    {
        out << "*CONTROL_IMPLICIT_SOLVER\n"
            << lsolvr << "," << lprint << "\n";
    }

    // Purpose: Stop the job.
    //
    // Parameters:
    // endtim Termination time. Mandatory.
    // ...
    void CONTROL_TERMINATION(std::ostream& out,
                             const lsFloat endtim)
    {
        out << "*CONTROL_TERMINATION\n";
        out << std::setprecision(3) << std::scientific;
        out << endtim << "\n";
    }

    // Options for binary output.
    //
    // Parameters:
    // timeInterval This field defines the time interval between output states,
    //              DT, for all options except D3DUMP, RUNRSF, and D3DRLF.
    void DATABASE_BINARY_D3PLOT(std::ostream& out,
                                const lsFloat timeInterval)
    {
        out << "*DATABASE_BINARY_D3PLOT\n";
        out << std::setprecision(3) << std::scientific;
        out << timeInterval << "\n";
    }


    // Options for ASCII files include.
    //
    // Parameters:
    // dt Time interval between outputs. If dt is zero, no output is printed.
    void DATABASE_GLSTAT(std::ostream& out,
                         const lsFloat dt)
    {
        out << "*DATABASE_GLSTAT\n";
        out << std::setprecision(3) << std::scientific;
        out << dt << "\n";

    }


    // Purpose: The keyword, *KEYWORD, flags LS-DYNA that the input deck is a
    // keyword deck rather than the structured format, which has a strictly
    // defined format. This must be the first card in the input file.
    //
    // Parameters:
    // memory the memory size to be used in words
    void KEYWORD(std::ostream& out,
                 const lsInt memory)
    {
        out << "*KEYWORD " << memory << "M" << "\n";
    }

    // Purpose: Define job title.
    //
    // Parameters:
    // title Heading to appear on output and in output files.
    void TITLE(std::ostream& out,
               const std::string& title)
    {
        out << "*TITLE\n";
        out << title << "\n";
    }

    // Purpose: Define a curve [for example, load (ordinate value) versus time
    // (abscissa value)], often loosely referred to as a load curve.
    //
    // Parameters:
    // lcId     Load curve ID.
    // a1, a2   Abscissa values.
    // o1, o2   Ordinate (function) values.
    void DEFINE_CURVE(std::ostream& out,
                      const lsInt lcId,
                      const lsFloat a1,
                      const lsFloat a2,
                      const lsFloat o1,
                      const lsFloat o2)
    {
        std::streamsize prec = out.precision();
        out << std::setprecision(13);
        
        out << "*DEFINE_CURVE\n"
            << std::setw(10) << lcId << "\n"
            << std::setw(20) << a1 << std::setw(20) << o1 << "\n"
            << std::setw(20) << a2 << std::setw(20) << o2 << "\n";

        out << std::setprecision(prec);
    }

    // This is Material Type 1. This is an isotropic hypoelastic material and is
    // available for beam, shell, and solid elements in LS-DYNA.
    //
    // Parameters:
    // materialId   Material identification. A unique number or label not
    //              exceeding 8 characters must be specified.
    // others: self-explanatory
    void MAT_ELASTIC(std::ostream& out,
                     const lsInt materialId,
                     const lsFloat massDensity,
                     const lsFloat youngsModulus,
                     const lsFloat poissonsRatio)
    {
        std::streamsize prec = out.precision();
        out << std::setprecision(3);

        
        out << "*MAT_ELASTIC\n"
            << std::setw(10) << materialId <<
               std::setw(10) << massDensity <<
               std::setw(10) << youngsModulus <<
               std::setw(10) << poissonsRatio << "\n";

        out << std::setprecision(prec);
    }

    // Purpose: Define parts, i.e., combine material information, section
    // properties, hourglass type, thermal properties, and a flag for part
    // adaptivity.
    //
    // Parameters:
    // partId      Part identification. A unique number or label must be specified.
    // sectionId   Section identification defined in a *SECTION keyword.
    // materialId  Material identification defined in the *MAT section.
    // heading     Heading for the part
    void PART(std::ostream& out,
              const lsInt partId,
              const lsInt sectionId,
              const lsInt materialId,
              const std::string& heading = "")
    {
        out << "*PART\n"
            << heading << "\n"
            << std::setw(10) << partId
            << std::setw(10) << sectionId
            << std::setw(10) << materialId << "\n";
    }

    // Purpose: Define section properties for shell elements.
    //
    // Parameters:
    // sectionId        Section ID. SECID is referenced on the *PART card.
    // elementForm      Element formulation options...
    // shearCorrection  Shear correction factor which scales the transverse shear stress.
    // numIntPoints     Number of through thickness integration points.
    // shellThickness   Shell thickness at node n1.
    void SECTION_SHELL(std::ostream& out,
                       const lsInt sectionId,
                       const lsInt elementForm,
                       const lsFloat shearCorrection,
                       const lsInt numIntPoints,
                       const lsFloat shellThickness)
    {
        std::streamsize prec = out.precision();
        out << std::setprecision(3);
        
        out << "*SECTION_SHELL\n"
            << sectionId << ","
            << elementForm << ","
            << shearCorrection << ","
            << numIntPoints << "\n"
            << shellThickness << "\n";

        out << std::setprecision(prec);
    }

    // Purpose: Define a node and its coordinates in the global coordinate
    // system. Also, the boundary conditions in global directions can be
    // specified.
    //
    // Parameters:
    // points     columnwise coordinates of points
    // onBoundary if true then point is on boundary of the patch
    // offset     starting index for the nodes
    template <typename T>
    void NODES(std::ostream& out,
               const gsMatrix<T>& points,
               const gsVector<bool>& onBoundary,
               const index_t offset)
    {
        std::streamsize prec = out.precision();
        out << std::setprecision(9);

        out << "*NODES\n";

        for (index_t col = 0; col != points.cols(); col++)
        {
            out << std::setw(8) << (offset + col);
            
            for (index_t dim = 0; dim != points.rows(); dim++)
            {
                out << std::setw(16) << points(dim, col);
            }

            T translConst = (onBoundary(col)) ? 7.0 : 0.0;

            
            out << std::fixed << std::setprecision(0) << std::setw(7)
                << translConst << "." << std::setw(7) << 0. << ".\n";

            out << std::scientific << std::setprecision(9);
        }

        out << std::setprecision(prec);
    }


    //  Purpose: Define three, four, six, and eight node elements including 3D
    //  shells, membranes, 2D plane stress, plane strain, and axisymmetric
    //  solids.
    //
    // Parameters:
    // elementId Element ID. Chose a unique number with respect to other elements.
    // partId    Part ID, see *PART.
    // nodes     indices of the nodes
    // offset    starting index for the nodes
    void ELEMENT_SHELL(std::ostream& out,
                        const lsInt elementId,
                        const lsInt partId,
                        const gsVector<index_t>& nodes,
                        const index_t offset)
    {
        out << "*ELEMENT_SHELL\n"
            << elementId << "," << partId;

        for (index_t row = 0; row != nodes.rows(); row++)
        {
            out << "," << (offset + nodes(row));
        }

        out << "\n";

    }

    // Purpose: Define a general 3D shell formulation to be used in combination
    // with *ELEMENT_GENERALIZED_SHELL. The objective of this feature is to
    // allow the rapid prototyping of new shell element formulations by adding
    // them through the keyword input file.
    //
    // Parameters:
    //
    // elform Element Formulation ID referenced via *SECTION_SHELL to connect
    //        *ELEMENT_GENERALIZED_SHELL with the appropriate shell
    //        formulation. The chosen number needs to be greater or equal than
    //        1000.
    //
    // numIntPts Number of in-plane integration points.
    //
    // numNodes Number of nodes for this element formulation.
    //
    // imass  Option for lumping of mass matrix:
    //        EQ.0: row sum
    //        EQ.1: diagonal weighting.
    //
    // iform Shell formulation to be used
    //       0, 1, 2, 3
    void DEFINE_ELEMENT_GENERALIZED_SHELL_header(std::ostream& out,
                                                 const lsInt elform,
                                                 const lsInt numIntPts,
                                                 const lsInt numNodes,
                                                 const lsInt imass,
                                                 const lsInt iform)
    {
        out << "*DEFINE_ELEMENT_GENERALIZED_SHELL\n"
            << std::setw(10) << elform << std::setw(10) << numIntPts
            << std::setw(10) << numNodes << std::setw(10) << imass
            << std::setw(10) << iform << "\n";

    }
                                                 

    // Purpose: Look DEFINE_ELEMENT_GENERALIZED_SHELL_header.
    //
    // Parameters:
    // quWeight ... quadrature weights
    // values   ... values of all basis functions
    // derivs   ... derivative of all basis functions
    template <typename T>
    void DEFINE_ELEMENT_GENERALIZED_SHELL_values(std::ostream& out,
                                                 const gsVector<T>& quWeights,
                                                 const gsMatrix<T>& values,
                                                 const gsMatrix<T>& derivs)
    {
        std::streamsize prec = out.precision();
        
        for (index_t k = 0; k != quWeights.rows(); k++)
        {
            gsMatrix<T> der = derivs.col(k);
            der.resize(2, values.rows());

            out << "$ gauss point " << k + 1 << "\n";
            out << std::setprecision(13);
            out << std::setw(20) << quWeights(k) << "\n";

            for (index_t i = 0; i != values.rows(); i++)
            {
                out << std::setw(20) << values(i, k);
                for (index_t u = 0; u != 2; u++) // parametric directions
                {
                    out << std::setw(20) << der(u, i);
                }
                out << "\n";
            }
        }
        
        out << std::setprecision(prec);
    }
                                                 
                                          
    
} // end namespace write    


} // end namespace lsdyna
    
} // end namespace gismo














