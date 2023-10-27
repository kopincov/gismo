/** @file gsQualityMeasureTest.cpp

    @brief ....

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): ...
*/

#include <gismo.h>
#include "gsMotorIOUtils.h"
//#include "gsMotorUtils.h"
//#include "gsUtils/gsExportMatrix.h"

using namespace gismo;

int main(int argc, char *argv[])
{
    // This code should be revised (problems with cdash)
    return 0;

//    const real_t orthogonality = 1.0;
// //    const real_t skewness = 0.0;
//    const real_t eccentricity = 0.0;
//    const real_t uniformity = 0.0;
// //    const real_t length = 0.0;
// //    const real_t area = 0.0;
// //    const real_t intersection = 0.0;
// //    const real_t epsilon = 1e-7;
// //    const bool dumped = false;

//    // Options with default values
//    std::string inputFile(MOTOR_DATA_DIR "jku/target_input.xml");

//    gsCmdLine cmd("...");
//    cmd.addString("f", "input", "Name of input file", inputFile);

//    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

//    gsInfo << "Input arguments:\n"
//           << "Input file:        " << inputFile << "\n"
//           << "--------------------------------------------------\n" << std::endl;

//    gsGeometry<>* geom = readGeometry(inputFile);
//    gsGeometry<>& mGeom = *geom;

//    const int matrixSize = mGeom.basis().size();
//    const short_t parDim = mGeom.parDim();
// //    const short_t geoDim = mGeom.geoDim();

//    gsMatrix<> basisDer;
//    gsMatrix<> geomDer;
//    gsMatrix<> basis2Der;
//    gsMatrix<> geom2Der;

//    gsSparseMatrix<> A(matrixSize, matrixSize);
//    A.setZero();

//    gsMatrix<> b(matrixSize, 1);
//    b.setZero();

//    gsVector<int> numNodes(parDim);
//    for (int i = 0; i != parDim; ++i)
//    {
//        numNodes[i] = mGeom.basis().degree(i);
//    }

//    gsBasis<>::domainIter domIt = mGeom.basis().makeDomainIterator();
//    domIt->computeQuadratureRule( numNodes );
//    for (; domIt->good(); domIt->next())
//    {

//        gsMatrix<> geomDer;
//        gsMatrix<> geom2Der;

//        geom->jacobian_into(domIt->quNodes, geomDer);

//        if (0 < eccentricity || 0 < uniformity)
//        {
//            mGeom.basis().deriv2_into(domIt->quNodes, basis2Der);
//            mGeom.deriv2_into(domIt->quNodes, geom2Der);
//        }

// //        const index_t numActive = domIt->computeActiveFunctions().rows();

//        for (index_t quNode = 0; quNode < domIt->numQuNodes(); ++quNode)
//        {
//            if (0 < orthogonality)
//            {
//                gsMatrix<> Jort, q;
// //                getJort(Jort, numActive, basisDer, geomDer, quNode);
// //                getQort(q, geomDer, quNode);
// //                assemble(A, b, Jort, q, orthogonality, domIt, quNode);
//            }
// /*
//            if (0 < skewness)
//            {
//                gsMatrix<> Jskew, q;
// //                getJskew(Jskew, numActive, basisDer, geomDer, quNode);
// //                getQskew(q, geomDer, quNode);
// //                assemble(A, b, Jskew, q, skewness, domIt, quNode);
//            }

//            if (0 < eccentricity)
//            {
//                gsMatrix<> Jecc, q;
// //                getJecc(Jecc, numActive, basisDer, geomDer, basis2Der,
// //                        geom2Der, quNode);
// //                getQecc(q, geomDer, geom2Der, quNode);
// //                assemble(A, b, Jecc, q, eccentricity, domIt, quNode);
//            }

//            if (0 < intersection)
//            {
//                gsMatrix<> Jsint, q;
// //                getJsint(Jsint, numActive, basisDer, geomDer, epsilon, quNode);
// //                getQsint(q, geomDer, epsilon, quNode);
// //                assemble(A, b, Jsint, q, intersection, domIt, quNode);
//            }

//            if (0 < uniformity)
//            {
//                gsMatrix<> Juni, q;
// //                getJuni(Juni, numActive, basis2Der, quNode);
// //                getQuni(q, geom2Der, quNode);
// //                assemble(A, b, Juni, q, uniformity, domIt, quNode);
//            }

//            if (0 < length)
//            {
//                gsMatrix<> Jlen, q;
// //                getJlen(Jlen, numActive, basisDer, quNode);
// //                getQlen(q, geom2Der, quNode);
// //                assemble(A, b, Jlen, q, length, domIt, quNode);
//            }

//            if (0 < area)
//            {
//                gsMatrix<> Jarea, q;
// //                getJarea(Jarea, numActive, basisDer, geomDer, quNode);
// //                getQarea(q, geomDer, quNode);
// //                assemble(A, b, Jarea, q, area, domIt, quNode);
//            }
//*/
//        }

// //        A.makeCompressed();

//        gsSparseSolver<>::BiCGSTABILUT solver(A);
//        if (solver.preconditioner().info() != Eigen::Success)
//        {
//            gsWarn << "The preconditioner failed. Aborting...\n";
//            return -1;
//        }

//        gsMatrix<> x;
//        x = solver.solve(b);


//        return 0;
//    }


    /* **********************************************************
    template<class T>
    gsFitting<T>::  gsFitting(gsMatrix<T> const & param_values,
                              gsMatrix<T> const & points,
                              gsBasis<T>  & basis)
    {
        m_param_values = param_values;
        m_points = points;
        m_result = NULL;
        m_basis = &basis;
        m_points.transposeInPlace();
    }

    template <class T>
    void gsFitting<T>::assembleSystem(gsSparseMatrix<T>& A_mat,
                      gsMatrix<T>& m_B)
    {
        const int num_points = m_points.rows();

        //for computing the value of the basis function
        gsMatrix<T> value;
        gsMatrix<unsigned> actives;


        for(index_t k = 0; k < num_points; k++)
        {
            const gsMatrix<T>& curr_point = m_param_values.col(k);

            //computing the values of the basis functions at the current point
            m_basis->eval_into(curr_point, value);

            // which functions have been computed i.e. which are active
            m_basis->active_into(curr_point, actives);

        const index_t numActive = actives.rows();

            for (index_t i = 0; i != numActive; ++i)
            {
                const int ii = actives(i, 0);
                m_B.row(ii) += value(i, 0) * m_points.row(k);
                for (index_t j = 0; j != numActive; ++j)
                    A_mat(ii, actives(j, 0)) += value(i, 0) * value(j, 0);
            }
        }
    }
    ********************************************************** */

    //return 0;
}

