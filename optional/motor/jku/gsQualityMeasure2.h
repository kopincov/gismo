// File gsQualityMeasure.h part of G+SMo library
// Author: Jaka Speh
#pragma once

#include <gsCore/gsGeometry.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMath.h>

namespace gismo
{
/**
    \brief
    Class with wich you can improve the parameterization (geometry) by
    relocating (inner) control points. Control points are moved such that
    certain quality measures are minimized.

    We can optimize with respect to many functionals. For example, let
    \f$ \boldsymbol{g} \f$ be an input geometry

    \f[
    \boldsymbol{g} : (u, v)^T \in \boldsymbol{\Omega} \mapsto \mathbb R^m,
    \f]

    with (inner) control points \f$ \boldsymbol{P}^* \f$.

    Then the functionals are:
      - Orthogonality
        \f[
        \mathcal{Q}_o( \boldsymbol{P}^* )
        =
        \int_{ \boldsymbol{\Omega} }
        \left( \boldsymbol{g}_u \boldsymbol{g}_u \right)^2
        d \boldsymbol{\Omega}.
        \f]

      - Skewness
        \f[
        \mathcal{Q}_s( \boldsymbol{P}^* )
        =
        \int_{ \boldsymbol{\Omega} }
        \left(
        \frac{
        ( \boldsymbol{g}_u \boldsymbol{g}_v )^2
        }{
        \boldsymbol{g}_u \boldsymbol{g}_u \cdot
        \boldsymbol{g}_v \boldsymbol{g}_v
        }
        \right)^2
        d \boldsymbol{\Omega}.
        \f]

      - Eccentricity
        \f[
        \mathcal{Q}_e (\boldsymbol{P}^*)
        =
        \int_{\boldsymbol{\Omega} }
        \left(
        \frac{
        \boldsymbol{g}_u \boldsymbol{g}_{uu}
        }{
        \boldsymbol{g}_u \boldsymbol{g}_u
        }
        \right)^2
        +
        \left(
        \frac{
        \boldsymbol{g}_v \boldsymbol{g}_{vv}
        }{
        \boldsymbol{g}_v \boldsymbol{g}_v
        }
        \right)^2
        d \boldsymbol{\Omega}.
        \f]

      - Self intersection
        \f[
        \mathcal{Q}_i(\boldsymbol{P}^*)
        =
        \int_{\boldsymbol{\Omega}}
        \left(
        \frac{1}{(\boldsymbol{g}_u^2 + \varepsilon)^2} +
        \frac{1}{(\boldsymbol{g}_v^2 + \varepsilon)^2}
        \right)
        d \boldsymbol{\Omega}.
        \f]
        where \f$\varepsilon\f$ is a small tolearance (for example \f$10^{-7}\f$).

      - Uniformity
       \f[
        \mathcal{Q}_u(\boldsymbol{P}^*)
        =
        \int_{\boldsymbol{\Omega}}
        (\| \boldsymbol{g}_{uu} \|^2 + 2\| \boldsymbol{g}_{uv} \|^2 +
        \| \boldsymbol{g}_{vv} \|^2)
        d \boldsymbol{\Omega}.
       \f]

     - Length
       \f[
        \mathcal{Q}_{\ell}(\boldsymbol{P}^*)
        =
        \int_{\boldsymbol{\Omega}}
        (\| \boldsymbol{g}_{u} \|^2 + \| \boldsymbol{g}_{v} \|^2)
        d \boldsymbol{\Omega}.
       \f]

     - Area
       \f[
        \mathcal{Q}_{a}(\boldsymbol{P}^*)
        =
        \int_{\boldsymbol{\Omega}}
        |\boldsymbol{g}_u \times \boldsymbol{g}_v|^2
        d \boldsymbol{\Omega}.
       \f]

    The optimization is done with Gauss-Newton method. We can see the measures
    as
    \f[
    \mathcal{Q} (\boldsymbol{P}^*)
    =
    \int_{\boldsymbol{\Omega}}
    \sum_{i = 1}^m q_i(\boldsymbol{P}^*)^2
    d \boldsymbol{\Omega}
    \f]
    and solve it with Gauss-Newton method. Look at
    http://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm

    Let \f$\boldsymbol{q} = [q_1, q_2, \ldots, q_m]^T\f$. For the Gauss-Newton
    algorithm, we need  to compute the Jacobian matrix of $q$
    \f[
    J_{\boldsymbol{P}^*}
    =
    \left[
    \begin{array}{ccc}
    \frac{ \partial{q_1}(\boldsymbol{P}^*)}{\partial \boldsymbol{P}^*_0}
    &
    \cdots
    &
    \frac{ \partial{q_1}(\boldsymbol{P}^*)}{\partial \boldsymbol{P}^*_n}
    \\
    \vdots & & \vdots
    \\
    \frac{ \partial{q_m}(\boldsymbol{P}^*)}{\partial \boldsymbol{P}^*_0}
    &
    \cdots
    &
    \frac{ \partial{q_m}(\boldsymbol{P}^*)}{\partial \boldsymbol{P}^*_n}
    \end{array}
    \right]
    \f]
    in order to assemble the system and solve it for \f$\boldsymbol{P}^*\f$.

    Each functional has two methods
      - getJ(functional) with arguments:
           - number of active functions (\a numActive)
           - first and second derivatives of the basis and geometry, only where necessary
             (\a basisDer, \a geomDer, \a basis2Der, \geom2Der)
           - index of quadrature node (\a quNode)

        and output argument
           - J(functional) which returns Jacobian matrix of functional at quadrature
             node \a quNode
      - getQ(functional) with arguments
           - first and second derivatives of the geometry (\a geomDer, \geom2Der)
           - index of quadrature node (\a quNode)

        and output argument
           - \a q the functional \f$\boldsymbol{q}\f$, as defined above, evaluated at
             quadrature node \a quNode

 */

template <typename T>
class gsQualityMeasure2
{
public:

    /// \brief Quality measure constructor
    ///
    /// \param geom input geometry
    /// \param fixedBoundary if true optimization moves just inner control
    ///                      points holding boudnary control points fixed
    gsQualityMeasure2(gsGeometry<T>& geom,
                      const gsMatrix<T>& params,
                      const gsMatrix<T>& points,
                      const bool fixedBoundary = true)
        : mGeom(geom),
          mParams(params),
          mPoints(points),
          mFixedBoundary(fixedBoundary),
          mMatrixToBasis(),
          mBasisToMatrix()
    {

        if (mFixedBoundary)
        {
            setFixedMappings();
        }
        else
        {
            setIdentityMappings();
            //GISMO_ERROR("Version with not fixed boundary points is not "
            //            "implemented.");
        }

    }

    /// \brief Move (inner) control points such that linear combination of
    ///        quality measures is minimized.
    ///
    /// Input parameters weight individual quality measures.
    ///
    /// \param[in] orthogonality weight for quality measure: orthogonality
    /// \param[in] skewness weight for quality measure: skewness
    /// \param[in] eccentricity weight for quality measure: eccentricity
    /// \param[in] uniformity weight for quality measure: uniformity
    /// \param[in] length weight for quality measure: length
    /// \param[in] area weight for quality measure: area
    /// \param[in] intersection weight for quality measure self intersection
    /// \param[in] epsilon parameter used in quality measure self intersection
    /// \param[in] dumped true if we use dumped method
    void optimize(const T fitting,
                  const T orthogonality,
                  const T skewness,
                  const T eccentricity,
                  const T uniformity,
                  const T length,
                  const T area,
                  const T intersection,
                  const T epsilon = 1e-7,
                  const bool dumped = false)
    {
        const int matrixSize = mMatrixToBasis.size();
        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();

        gsMatrix<T> basisDer;
        gsMatrix<T> geomDer;
        gsMatrix<T> basis2Der;
        gsMatrix<T> geom2Der;

        gsSparseMatrix<T> A(matrixSize * geoDim, matrixSize * geoDim);
        A.setZero();

        gsMatrix<T> b(matrixSize * geoDim, 1);
        b.setZero();

        gsVector<int> numNodes(parDim);
        for (int i = 0; i != parDim; ++i)
        {
            numNodes[i] = mGeom.basis().degree(i);
        }

        typename gsBasis<T>::domainIter domIt = mGeom.basis().makeDomainIterator();
        gsGaussRule<T> quRule(numNodes);
        gsMatrix<T> quNodes;
        gsVector<T> quWeights;
        gsMatrix<index_t> actives;

        assembleFitting(A, fitting);
//        gsInfo << "A.nonZeros():\n" << A.nonZeros() << std::endl;
//        gsInfo << "A.nonZerosPerCol():\n" << A.nonZerosPerCol() << std::endl;

        for (; domIt->good(); domIt->next())
        {
            quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

            mGeom.basis().deriv_into(quNodes, basisDer);
            mGeom.jacobian_into(quNodes, geomDer);
            mGeom.basis().active_into(quNodes.col(0), actives);

            if (0 < eccentricity || 0 < uniformity)
            {
                mGeom.basis().deriv2_into(quNodes, basis2Der);
                mGeom.deriv2_into(quNodes, geom2Der);
            }

            const index_t numActive = actives.rows();

            for (index_t quNode = 0; quNode < quWeights.size(); ++quNode)
            {
                // SS
                if (0 < fitting)
                {
/*
                    gsMatrix<T> Jort, q;
                    getJort(Jort, numActive, basisDer, geomDer, quNode);
                    getQort(q, geomDer, quNode);
                    assemble(A, b, Jort, q, quWeights[quNode] * orthogonality, actives);
*/
                }

                if (0 < orthogonality)
                {
                    gsMatrix<T> Jort, q;
                    getJort(Jort, numActive, basisDer, geomDer, quNode);
                    getQort(q, geomDer, quNode);
                    assemble(A, b, Jort, q, quWeights[quNode] * orthogonality, actives);
                }

                if (0 < skewness)
                {
                    gsMatrix<T> Jskew, q;
                    getJskew(Jskew, numActive, basisDer, geomDer, quNode);
                    getQskew(q, geomDer, quNode);
                    assemble(A, b, Jskew, q, quWeights[quNode] * skewness, actives);
                }

                if (0 < eccentricity)
                {
                    gsMatrix<T> Jecc, q;
                    getJecc(Jecc, numActive, basisDer, geomDer, basis2Der,
                            geom2Der, quNode);
                    getQecc(q, geomDer, geom2Der, quNode);
                    assemble(A, b, Jecc, q, quWeights[quNode] * eccentricity, actives);
                }

                if (0 < intersection)
                {
                    gsMatrix<T> Jsint, q;
                    getJsint(Jsint, numActive, basisDer, geomDer, epsilon, quNode);
                    getQsint(q, geomDer, epsilon, quNode);
                    assemble(A, b, Jsint, q, quWeights[quNode] * intersection, actives);
                }

                if (0 < uniformity)
                {
                    gsMatrix<T> Juni, q;
                    getJuni(Juni, numActive, basis2Der, quNode);
                    getQuni(q, geom2Der, quNode);
                    assemble(A, b, Juni, q, quWeights[quNode] * uniformity, actives);
                }

                if (0 < length)
                {
                    gsMatrix<T> Jlen, q;
                    getJlen(Jlen, numActive, basisDer, quNode);
                    getQlen(q, geom2Der, quNode);
                    assemble(A, b, Jlen, q, quWeights[quNode] * length, actives);
                }

                if (0 < area)
                {
                    gsMatrix<T> Jarea, q;
                    getJarea(Jarea, numActive, basisDer, geomDer, quNode);
                    getQarea(q, geomDer, quNode);
                    assemble(A, b, Jarea, q, quWeights[quNode] * area, actives);
                }

            }
        }

        A.makeCompressed();


        gsSparseSolver<>::CGIdentity solver(A);
        //gsSparseSolver<>::CGDiagonal solver(A);
        //gsSparseSolver<>::BiCGSTABIdentity solver(A);
        //gsSparseSolver<>::BiCGSTABDiagonal solver(A);
        //gsSparseSolver<>::BiCGSTABILUT solver(A);

        if (solver.preconditioner().info() != Eigen::Success)
        {
            gsWarn << "The preconditioner failed. Aborting...\n";
            return;
        }

        //----------------------EIGEN-DIRECT-SOLVERS----------------------//
        //gsSparseSolver<>::SimplicialLDLT solver(A);
        //gsSparseSolver<>::QR solver(A);
        //gsSparseSolver<>::LU solver(A);
        //gsSparseSolver<>::fullPivLu solver(A);

        gsMatrix<T> x;
        x = solver.solve(b);

//        std::cout << "size x: " << x.rows() << " x " << x.cols() << std::endl;
//        std::cout << "A: \n" << A << "\n\n"
//                  << "b: \n" << b << "\n\n"
//                  << std::endl;





        gsMatrix<T>& coefs = mGeom.coefs();

        if (dumped) // Q: What is dumped?
        {
            gsMatrix<T> oldCoefs = mGeom.coefs(); // important copy old coefs !!

            T startFunc = functional(fitting, orthogonality, skewness, eccentricity,
                                     uniformity, length, area, intersection,
                                     epsilon);

            for (int i = 0; i != 10; i++)
            {
                const T factor = T(1.0) / T(1 << i);

                for (int row = 0; row != coefs.rows(); row++)
                {
                    const int index = mBasisToMatrix[row];
                    if (index == -1)
                    {
                        continue;
                    }
                    for (int dim = 0; dim != geoDim; ++dim)
                    {
                        coefs(row, dim) = oldCoefs(row, dim) -
                                factor * x(index * geoDim + dim, 0);
                    }
                }

                T newFunc = functional(fitting, orthogonality, skewness,
                                       eccentricity, uniformity, length,
                                       area, intersection, epsilon);

                if (newFunc < startFunc)
                {
                    break;
                }
            }
        }
        else
        {
            for (int row = 0; row != coefs.rows(); row++)
            {
                const int index = mBasisToMatrix[row];
                if (index == -1)
                {
                    continue;
                }
                for (int dim = 0; dim != geoDim; ++dim)
                {
                    coefs(row, dim) -= x(index * geoDim + dim, 0);
                }
            }
        }
    }


    /// \brief Computes value of linear combination of quality measures.
    ///
    /// Input parameters weight individual quality measures.
    ///
    /// For the desription of input parameters look at the method optimize()
    ///
    /// \return value of the combined quality measures
    T functional(const T fitting,
                 const T orthogonality,
                 const T skewness,
                 const T eccentricity,
                 const T uniformity,
                 const T length,
                 const T area,
                 const T intersection,
                 const T epsilon = 1e-7)
    {
        const short_t parDim = mGeom.parDim();
        //const short_t geoDim = mGeom.geoDim();

        gsVector<int> numNodes(parDim);
        for (int i = 0; i != parDim; ++i)
        {
            numNodes[i] = mGeom.basis().degree(i);
        }

        T functional = 0;

        typename gsBasis<T>::domainIter domIt = mGeom.basis().makeDomainIterator();

        gsGaussRule<T> quRule(numNodes);
        gsMatrix<T> quNodes;
        gsVector<T> quWeights;

        for (; domIt->good(); domIt->next())
        {
            quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

            gsMatrix<T> geomDer;
            gsMatrix<T> geom2Der;

            mGeom.jacobian_into(quNodes, geomDer);

            if (0 < eccentricity || 0 < uniformity)
            {
                mGeom.deriv2_into(quNodes, geom2Der);
            }

            for (index_t quNode = 0; quNode < quWeights.size(); ++quNode)
            {
                const T weight = quWeights[quNode];

                // SS
                if (0 < fitting)
                {
                    gsMatrix<T> q;
                    //getQort(q, geomDer, quNode);
                    functional += assembleFunctional(q, weight, fitting);
                }

                if (0 < orthogonality)
                {
                    gsMatrix<T> q;
                    getQort(q, geomDer, quNode);
                    functional += assembleFunctional(q, weight, orthogonality);
                }

                if (0 < skewness)
                {
                    gsMatrix<T> q;
                    getQskew(q, geomDer, quNode);
                    functional += assembleFunctional(q, weight, skewness);
                }

                if (0 < eccentricity)
                {
                    gsMatrix<T> q;
                    getQecc(q, geomDer, geom2Der, quNode);
                    functional += assembleFunctional(q, weight, eccentricity);
                }

                if (0 < intersection)
                {
                    gsMatrix<T> q;
                    getQsint(q, geomDer, epsilon, quNode);
                    functional += assembleFunctional(q, weight, intersection);
                }

                if (0 < uniformity)
                {
                    gsMatrix<T> q;
                    getQuni(q, geom2Der, quNode);
                    functional += assembleFunctional(q, weight, uniformity);
                }

                if (0 < area)
                {
                    gsMatrix<T> q;
                    getQarea(q, geomDer, quNode);
                    functional += assembleFunctional(q, weight, area);
                }

                if (0 < length)
                {
                    gsMatrix<T> q;
                    getQlen(q, geomDer, quNode);
                    functional += assembleFunctional(q, weight, length);
                }

            }
        }
        return functional;
    }


private:
    void getJort(gsMatrix<T>& Jort,
                 const index_t numActive,
                 const gsMatrix<T>& basisDer,
                 const gsMatrix<T>& geomDer,
                 const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();

        // ** why nit numQ = 1 since it is scalar ? -->(row)
        const int numQ = (parDim - 1) * parDim / 2;
        //const index_t numActive = domIt->computeActiveFunctions().rows();

        Jort.resize(numQ * geoDim, numActive);
        Jort.setZero();

        // computing Jort
        for (index_t i = 0; i != numActive; ++i)
        {
            int row = 0;
            for (index_t u = 0; u != parDim; ++u)
            {
                for (index_t v = u + 1; v != parDim; ++v)
                {
                    const int uBasInd = parDim * i + u;
                    const int vBasInd = parDim * i + v;

                    const int uGeomInd = quNode * parDim + u;
                    const int vGeomInd = quNode * parDim + v;

                    Jort.block(row * geoDim, i, geoDim, 1) =
                            basisDer(uBasInd, quNode) * geomDer.col(vGeomInd) +
                            basisDer(vBasInd, quNode) * geomDer.col(uGeomInd);
                    row++;
                }
            }
        }
    }

    /// Look at class gsQualityMeasure
    void getJskew(gsMatrix<T>& Jskew,
                  const index_t numActive,
                  const gsMatrix<T>& basisDer,
                  const gsMatrix<T>&  geomDer,
                  const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();

        const int numQ = (parDim - 1) * parDim / 2;

        Jskew.resize(numQ * geoDim, numActive);
        Jskew.setZero();

        // computing Jskew
        for (index_t i = 0; i != numActive; ++i)
        {
            int row = 0;
            for (index_t u = 0; u != parDim; ++u)
            {

                for (index_t v = u + 1; v != parDim; ++v)
                {
                    const int uBasInd = parDim * i + u;
                    const int vBasInd = parDim * i + v;

                    const int uGeomInd = quNode * parDim + u;
                    const int vGeomInd = quNode * parDim + v;

                    const gsVector<T>& Gu = geomDer.col(uGeomInd);
                    const gsVector<T>& Gv = geomDer.col(vGeomInd);

                    const T Nu = basisDer(uBasInd, quNode);
                    const T Nv = basisDer(vBasInd, quNode);

                    const T GuGv = ( Gu.transpose() * Gv ).value() ;

                    const T GuGu = ( Gu.transpose() * Gu ).value() ;

                    const T GvGv = ( Gv.transpose() * Gv ).value() ;

                    Jskew.block(row * geoDim, i, geoDim, 1) =
                            (2 * GuGv * (Nu * Gv + Nv * Gu) * GuGu * GvGv -
                            GuGv * GuGv * 2 * (Nu * GvGv * Gu + Nv * GuGu * Gv))
                            / (GuGu * GuGu * GvGv * GvGv);
                    row++;
                }
            }
        }
    }


    /// Look at class gsQualityMeasure
    void getJecc(gsMatrix<T>& Jecc,
                 const index_t numActive,
                 const gsMatrix<T>& basisDer,
                 const gsMatrix<T>& geomDer,
                 const gsMatrix<T>& basis2Der,
                 const gsMatrix<T>& geom2Der,
                 const index_t quNode)
    {

        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();
        const int stride = parDim * (parDim + 1) / 2;

        Jecc.resize(parDim * geoDim, numActive);
        Jecc.setZero();

        // computing Jecc
        for (index_t i = 0; i != numActive; ++i)
        {
            for (index_t u = 0; u != parDim; ++u)
            {
                const int basIndDer = parDim * i + u;
                // pure 2nd deriv of basis in u coordinate
                const int basInd2Der = stride * i + u;

                const int uGeomInd = quNode * parDim + u;

                const gsVector<T>& Gu = geomDer.col(uGeomInd);

                // pure 2nd deriv of geometry in u coordinate
                gsVector<T> Guu(geoDim);
                for (int dim = 0; dim != geoDim; ++dim)
                {
                    Guu(dim) = geom2Der(stride * dim + u, quNode);
                }

                const T Nu = basisDer(basIndDer, quNode);
                const T Nuu = basis2Der(basInd2Der, quNode);

                const T GuGu  = ( Gu.transpose() * Gu ).value() ;
                const T GuGuu = ( Gu.transpose() * Guu).value() ;

                Jecc.block(u * geoDim, i, geoDim, 1) =
                        (GuGu * (Nu * Guu + Nuu * Gu) -
                         GuGuu * 2 * Nu * Gu) / (GuGu * GuGu);

            }
        }
    }

    /// Look at class gsQualityMeasure
    void getJuni(gsMatrix<T>& Juni,
                 const index_t numActive,
                 const gsMatrix<T>& basis2Der,
                 const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();
        const int stride = parDim * (parDim + 1) / 2;

        const int numQ = stride * geoDim;
        // functionals are
        // (G_uu)_x, (G_uu)_y, (G_vv)_x, (G_vv)_y, sqrt2 (G_uv)_x, sqrt2 (G_uv)_y

        Juni.resize(numQ * geoDim, numActive);
        Juni.setZero();

        for (index_t i = 0; i != numActive; ++i)
        {
            for (int s = 0; s != stride; ++s)
            {
                const int basInd2Der = i * stride + s;
                T Nuu = basis2Der(basInd2Der, quNode);

                if (parDim <= s)
                {
                    Nuu *= M_SQRT2;
                }

                for (int dimQ = 0; dimQ != geoDim; dimQ++)
                {
                    const int qIndex = s * geoDim + dimQ;

                    Juni(qIndex * geoDim + dimQ, i) = Nuu;

                }
            }
        }
    }

    /// Look at class gsQualityMeasure
    void getJlen(gsMatrix<T>& Jlen,
                 const index_t numActive,
                 const gsMatrix<T>& basisDer,
                 const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();

        const int numQ = parDim * geoDim;
        // functionals are
        // (G_u)_x, (G_u)_y, (G_v)_x, (G_v)_y

        Jlen.resize(numQ * geoDim, numActive);
        Jlen.setZero();

        for (index_t i = 0; i != numActive; ++i)
        {
            for (int u = 0; u != parDim; ++u)
            {
                const int basIndDer = parDim * i + u;
                const T Nu = basisDer(basIndDer, quNode);

                for (int dimQ = 0; dimQ != geoDim; dimQ++)
                {
                    const int qIndex = u * geoDim + dimQ;
                    Jlen(qIndex * geoDim + dimQ, i) = Nu;
                }
            }
        }
    }

    /// Look at class gsQualityMeasure
    void getJarea(gsMatrix<T>& Jarea,
                  const index_t numActive,
                  const gsMatrix<T>& basisDer,
                  const gsMatrix<T>& geomDer,
                  const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();

        GISMO_ASSERT(parDim == 2, "Area functional is implemented only for "
                                  "parametric dimension 2\n");

        GISMO_ASSERT(geoDim == 2, "Area functional is implemented only for "
                                  "geometric dimension 2\n");

        const int numQ = 1;

        Jarea.resize(numQ * geoDim, numActive);
        Jarea.setZero();

        for (index_t i = 0; i != numActive; ++i)
        {
            const int uGeomInd = quNode * parDim + 0;
            const int vGeomInd = quNode * parDim + 1;

            const int uBasInd = i * parDim + 0;
            const int vBasInd = i * parDim + 1;

            const gsVector<T>& Gu = geomDer.col(uGeomInd);
            const gsVector<T>& Gv = geomDer.col(vGeomInd);

            const T Nu = basisDer(uBasInd, quNode);
            const T Nv = basisDer(vBasInd, quNode);

            const T area = Gu(0) * Gv(1) - Gu(1) * Gv(0);
            const T constant = 2 * area;

            Jarea(0, i) = constant * (Nu * Gv(1) - Gu(1) * Nv);
            Jarea(1, i) = constant * (Gu(0) * Nv - Nu * Gv(0));
        }
    }


    /// Look at class gsQualityMeasure
    void getJsint(gsMatrix<T>& Jsint,
                  const index_t numActive,
                  const gsMatrix<T>& basisDer,
                  const gsMatrix<T>&  geomDer,
                  const T epsilon,
                  const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();

        Jsint.resize(parDim * geoDim, numActive);
        Jsint.setZero();

        // computing Jecc
        for (index_t i = 0; i != numActive; ++i)
        {
            for (index_t u = 0; u != parDim; ++u)
            {
                const int uBasInd = parDim * i + u;

                const int uGeomInd = quNode * parDim + u;

                const gsVector<T>& Gu = geomDer.col(uGeomInd);
                const T Nu = basisDer(uBasInd, quNode);

                const T GuGu = ( Gu.transpose() * Gu).value() ;

                Jsint.block(u * geoDim, i, geoDim, 1) =
                        (-2 * Nu * Gu) /
                        ((GuGu * GuGu + epsilon) * (GuGu * GuGu + epsilon));
            }
        }
    }

    /// Look at class gsQualityMeasure
    void getQort(gsMatrix<T>& q,
                 const gsMatrix<T>& geomDer,
                 const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();
        const int numQ = (parDim - 1) * parDim / 2;

        q.resize(1, numQ);
        q.setZero();

        int col = 0;
        // computing q
        for (index_t u = 0; u != parDim; ++u)
        {
            for (index_t v = u + 1; v != parDim; ++v)
            {
                const int uGeomInd = quNode * parDim + u;
                const int vGeomInd = quNode * parDim + v;

                q(0, col) = geomDer.col(uGeomInd).transpose() * geomDer.col(vGeomInd);
                col++;
            }
        }
    }

    /// Look at class gsQualityMeasure
    void getQskew(gsMatrix<T>& q,
                  const gsMatrix<T>&  geomDer,
                  const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();
        const int numQ = (parDim - 1) * parDim / 2;

        q.resize(1, numQ);
        q.setZero();

        // computing q
        int col = 0;
        for (index_t u = 0; u != parDim; ++u)
        {
            for (index_t v = u + 1; v != parDim; ++v)
            {
                const int uGeomInd = quNode * parDim + u;
                const int vGeomInd = quNode * parDim + v;

                const gsVector<T>& Gu = geomDer.col(uGeomInd);
                const gsVector<T>& Gv = geomDer.col(vGeomInd);

                const T GuGv = ( Gu.transpose() * Gv ).value();
                const T GuGu = ( Gu.transpose() * Gu ).value();
                const T GvGv = ( Gv.transpose() * Gv ).value();

                q(0, col) = (GuGv * GuGv) / (GuGu * GvGv);
            }
        }
    }


    /// Look at class gsQualityMeasure
    void getQecc(gsMatrix<T>& q,
                 const gsMatrix<T>& geomDer,
                 const gsMatrix<T>& geom2Der,
                 const index_t quNode)
    {

        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();

        const int stride = parDim * (parDim + 1) / 2;

        q.resize(1, parDim);
        q.setZero();

        for (index_t u = 0; u != parDim; ++u)
        {
            const int uGeomInd = quNode * parDim + u;
            const gsVector<T>& Gu = geomDer.col(uGeomInd);

            // pure 2nd deriv of geometry in u coordinate
            gsVector<T> Guu(geoDim);
            for (int dim = 0; dim != geoDim; ++dim)
            {
                Guu(dim) = geom2Der(stride * dim + u, quNode);
            }

            const T GuGuu = ( Gu.transpose() * Guu ).value() ;
            const T GuGu  = ( Gu.transpose() * Gu  ).value() ;
            q(0, u) = GuGuu / GuGu;
        }
    }

    /// Look at class gsQualityMeasure
    void getQuni(gsMatrix<T>& q,
                 const gsMatrix<T>& geom2Der,
                 const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();
        const int stride = parDim * (parDim + 1) / 2;

        const int numQ = stride * geoDim;

        q.resize(1, numQ);
        q.setZero();

        for (int s = 0; s != stride; ++s)
        {
            gsVector<T> Gss(geoDim);
            for (int dim = 0; dim != geoDim; ++dim)
            {
                Gss(dim) = geom2Der(stride * dim + s, quNode);
                if (parDim <= s)
                {
                    Gss(dim) *= M_SQRT2;
                }
            }

            for (int dimQ = 0; dimQ != geoDim; dimQ++)
            {
                const int qIndex = s * geoDim + dimQ;

                q(0, qIndex) = Gss(dimQ);
            }
        }
    }

    /// Look at class gsQualityMeasure
    void getQlen(gsMatrix<T>& q,
                 const gsMatrix<T>& geomDer,
                 const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();

        const int numQ = parDim * geoDim;
        q.resize(1, numQ);
        q.setZero();

        for (int u = 0; u != parDim; ++u)
        {
            const int uGeomInd = quNode * parDim + u;
            const gsVector<T>& Gu = geomDer.col(uGeomInd);

            for (int dimQ = 0; dimQ != geoDim; dimQ++)
            {
                const int qIndex = u * geoDim + dimQ;
                q(0, qIndex) = Gu(dimQ);
            }
        }
    }

    /// Look at class gsQualityMeasure
    void getQarea(gsMatrix<T>& q,
                  const gsMatrix<T>& geomDer,
                  const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();

        GISMO_ASSERT(parDim == 2, "Area functional is implemented only for "
                                  "parametric dimension 2\n.");

        GISMO_ASSERT(mGeom.geoDim() == 2, "Area functional is implemented only for "
                                  "geometric dimension 2\n");

        const int numQ = 1;

        q.resize(1, numQ);
        q.setZero();

        const int uGeomInd = quNode * parDim + 0;
        const int vGeomInd = quNode * parDim + 1;

        const gsVector<T>& Gu = geomDer.col(uGeomInd);
        const gsVector<T>& Gv = geomDer.col(vGeomInd);

        q(0, 0) = Gu(0) * Gv(1) - Gu(1) * Gv(0);
    }

    /// Look at class gsQualityMeasure
    void getQsint(gsMatrix<T>& q,
                  const gsMatrix<T>&  geomDer,
                  const T epsilon,
                  const index_t quNode)
    {
        const short_t parDim = mGeom.parDim();

        q.resize(1, parDim);
        q.setZero();

        for (index_t u = 0; u != parDim; ++u)
        {
            const int uGeomInd = quNode * parDim + u;
            const gsVector<T>& Gu = geomDer.col(uGeomInd);

            const T GuGu = ( Gu.transpose() * Gu ).value();
            const T tmp = GuGu * GuGu + epsilon;

            q(0, u) = 1 / (tmp * tmp);
        }
    }



    /// \brief Sets mappings: matrix index to basis index and vise versa.
    ///
    /// Mappings don't take into account the boudary.
    /// The function populates mMatrixToBasis and mBasisToMatrix vectors
    ///
    /// mBasisToMatrix[i] = j  <=>  if j == -1 then i-th basis function is
    /// on boundary and doesn't have place in the matrix
    ///                             if j != -1 then i-th basis function
    /// corresponds to i-th row in a matrix (and i-th column)
    ///
    /// mMatrixToBasis[j] = i <=> j-th row (or column) of a matrix corresponds
    /// to i-th basis function
    void setFixedMappings()
    {
        gsMatrix<index_t> boundary = mGeom.basis().allBoundary();
        int size = mGeom.basis().size();
        int matrixSize = size - boundary.rows();

        mMatrixToBasis.resize(matrixSize);
        mMatrixToBasis.setZero();
        mBasisToMatrix.resize(size);
        mBasisToMatrix.setConstant(-1);

        int boundaryCounter = 0;
        int matrixCounter = 0;
        for (int basisFun = 0; basisFun != size; basisFun++)
        {
            if (boundaryCounter != boundary.rows() &&
                basisFun == static_cast<int>(boundary(boundaryCounter, 0)))
            {
                boundaryCounter++;
            }
            else
            {
                mMatrixToBasis(matrixCounter) = basisFun;
                mBasisToMatrix(basisFun) = matrixCounter;
                matrixCounter++;
            }
        }
    }

    void setIdentityMappings()
    {
        //gsMatrix<index_t> boundary = mGeom.basis().allBoundary();
        int size = mGeom.basis().size();
        //int matrixSize = size - boundary.rows();
        int matrixSize = size;

        mMatrixToBasis.resize(matrixSize);
        mMatrixToBasis.setZero();
        mBasisToMatrix.resize(size);
        mBasisToMatrix.setConstant(-1);

        //int boundaryCounter = 0;
        int matrixCounter = 0;
        for (int basisFun = 0; basisFun != size; basisFun++)
        {
            /*
            if (boundaryCounter != boundary.rows() &&
                basisFun == static_cast<int>(boundary(boundaryCounter, 0)))
            {
                boundaryCounter++;
            }
            else
            {
                mMatrixToBasis(matrixCounter) = basisFun;
                mBasisToMatrix(basisFun) = matrixCounter;
                matrixCounter++;
            }
            */
            mMatrixToBasis(matrixCounter) = basisFun;
            mBasisToMatrix(basisFun) = matrixCounter;
            matrixCounter++;
        }
    }

    /// \brief Assembles the functional \f$ \boldsymbol{q}\f$ into
    /// \f$ \sum_i q_i^2 * \mbox{weight} * \mbox{constant} \f$
    ///
    /// \param[in] q quality measure
    /// \param[in] weight the integration weight
    /// \param[in] constant how big is the contribution of this quality measure
    ///
    /// \return value of the quality measure
    T assembleFunctional(const gsMatrix<T>& q,
                         const T weight,
                         const T constant)
    {
        T value = 0;
        for (int col = 0; col != q.cols(); col++)
        {
            value += q(0, col) * q(0, col);
        }
        return weight * constant * value;
    }

    /***********************************************************/
    void assembleFitting(gsSparseMatrix<T>& A,
                         const T constant)
    {
//        const short_t parDim = mGeom.parDim();
        const short_t geoDim = mGeom.geoDim();

        const int numPoints = mPoints.cols();
        //const index_t numActive = domIt->computeActiveFunctions().rows();

        //for computing the value of the basis function
        gsMatrix<T> value;
        gsMatrix<index_t> actives;

        for(index_t k = 0; k < numPoints; k++)
        {
            const gsMatrix<T>& currentPoint = mParams.col(k);

            //computing the values of the basis functions at the current point
            mGeom.basis().eval_into(currentPoint, value);

            // which functions have been computed i.e. which are active
            mGeom.basis().active_into(currentPoint, actives);

            const index_t numActive = actives.rows();

            for (int dim = 0; dim != geoDim; dim++)
            {
                for (index_t i = 0; i != numActive; ++i)
                {
                    const int indexI = mBasisToMatrix[actives(i, 0)];
                    for (index_t j = 0; j != numActive; ++j)
                    {
                        const int indexJ = mBasisToMatrix[actives(j, 0)];
                        A(indexI*geoDim + dim, indexJ*geoDim + dim) += constant * value(i, 0) * value(j, 0);
                    }

                }
            }
/*
            for (index_t i = 0; i != numActive; ++i)
            {
                const int indexI = mBasisToMatrix[actives(i, 0)];
                for (index_t j = 0; j != numActive; ++j)
                {
                    const int indexJ = mBasisToMatrix[actives(j, 0)];

                    for (int dimu = 0; dimu != geoDim; dimu++)
                    {
                        for (int dimv = 0; dimv != geoDim; dimv++)
                        {
                            A(indexI * geoDim + dimu, indexJ * geoDim + dimv) += constant * value(i, 0) * value(j, 0);
                        }
                    }
                }
            }
*/



/*                const int indexI = mBasisToMatrix[actives(i, 0)];
                for (int dim = 0; dim != geoDim; dim++)
                {
                        A(indexI * geoDim + dim, indexJ * geoDim + dimv) +=
                            tmp * weight;
                }
*/
        }
    }
    /***********************************************************/


    /// \brief Assembles matrix \a A and right hand sie vector \a b from Jacobian
    /// matrix \A J and functional values \a q.
    ///
    /// We solve a system \fA \boldsymbol{P}^* = b$ to get new control points.
    ///
    /// \param[in,out] A matrix of a system
    /// \param[in, out] b right hand side of a system
    /// \param[in] J Jacobian matrix of a quality measure
    /// \param[in] q values of a quality measure
    /// \param[in] weight weight of this quality measure
    /// \param[in] actives non zero basis function
    void assemble(gsSparseMatrix<T>& A,
                  gsMatrix<T>& b,
                  const gsMatrix<T>& J,
                  const gsMatrix<T>& q,
                  const T weight,
                  const gsMatrix<index_t> & actives)
    {
        const short_t geoDim = mGeom.geoDim();

        for (index_t i = 0; i != actives.rows(); ++i)
        {
            const int indexI = mBasisToMatrix[actives(i, 0)];
            if (indexI == -1) // we are on boundary, skip
            {
                continue;
            }

            // b part
            for (int col = 0; col != q.cols(); col++)
            {
                for (int dim = 0; dim != geoDim; dim++)
                {
                    b(indexI * geoDim + dim) += weight * q(0, col) *
                            J(col * geoDim + dim, i);
                }
            }

            // A part
            for (index_t j = 0; j != actives.rows(); ++j)
            {
                const int indexJ = mBasisToMatrix[actives(j, 0)];
                if (indexJ == -1) // we are on boundary, skip
                {
                    continue;
                }

                // A part
                for (int col = 0; col != q.cols(); col++)
                {
                    for (int dimu = 0; dimu != geoDim; dimu++)
                    {
                        for (int dimv = 0; dimv != geoDim; dimv++)
                        {
                            const T tmp = J(col * geoDim + dimu, i) *
                                          J(col * geoDim + dimv, j);

                            A(indexI * geoDim + dimu, indexJ * geoDim + dimv) +=
                                tmp * weight;
                        }
                    }
                }
            }
        }
    }

    // disable copy constructor
    gsQualityMeasure2(const gsQualityMeasure2<T>& other);

    // disable assignment operator
    gsQualityMeasure2<T>& operator=(const gsQualityMeasure2<T>& other);

private:

    // data members

    /// \brief Input geometry, parameterization.
    gsGeometry<T>& mGeom;

    /// \brief The parameter values of the point cloud
    gsMatrix<T> mParams;

    /// \brief The points of the point cloud
    gsMatrix<T> mPoints;

    /// \brief Indicates if the boundary is fixed.
    ///
    /// If mFixedBoundary is true we move just inner control points when we
    /// optimize the geometry with quality measures.
    bool mFixedBoundary;

    /// \brief Maps index of a matrix (row / column) to index of a basis function.
    gsVector<int> mMatrixToBasis;

    /// \brief Maps index of a basis to index of a matrix (row / column).
    gsVector<int> mBasisToMatrix;

};


} // end namespace gismo
