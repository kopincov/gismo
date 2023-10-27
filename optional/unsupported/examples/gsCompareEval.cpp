
#include <iostream>
#include <gismo.h>



using namespace gismo;

enum
{
    NONE = 1U << 0,

    EVAL = 1U << 1,
    DERIV = 1U << 2,
    DERIV2 = 1U << 3,

    OLD = 1U << 5,
    NEW = 1U << 6,
    VERY_OLD = 1U << 7
};


// checks if A and B match (element wise) by absolute value,
// if they don't match by absolute value, we try to match it by relative value
// if A and B match, we don't write anything
void checkTolerance(const gsMatrix<>& A,
                    const gsMatrix<>& B,
                    const real_t tol)
{
    for (index_t col = 0; col != A.cols(); col++)
    {
        for (index_t row = 0; row != A.rows(); row++)
        {
            const real_t diff = math::abs(A(row, col) - B(row, col));

            if (tol < diff)
            {
                const real_t rel = diff / math::abs(A(row, col));

                if (tol < rel)
                {
                    gsInfo << "error: "
                              << std::setw(25) << std::setprecision(16)
                              << A(row, col)
                              << "|" << std::setw(25)
                              << std::setprecision(16)
                              << B(row, col)
                              << "| " << std::setprecision(6) << diff
                              << "\n";
                }
            }
        }
    }
}


template <unsigned d>
void compareEval(gsTHBSplineBasis<d, real_t>& basis,
                 const unsigned flags)
{
    const real_t tol = 1e-15;

//    gsInfo << "Tolerance is: " << tol << "\n"
//              << "If there is no output below everything is ok..."
//              << "\n";


    const short_t parDim = basis.dim();

    gsVector<int> numNodes(parDim);
    for (int i = 0; i != parDim; ++i)
    {
        numNodes[i] = basis.degree(i);
    }

    gsBasis<>::domainIter domIt = basis.makeDomainIterator();
    gsGaussRule<> quRule(numNodes);
    gsMatrix<> quNodes;
    gsVector<> quWeights;
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

    real_t maxEval = -1;
    real_t maxDeriv = -1;
    real_t maxDeriv2 = -1;

    for (; domIt->good(); domIt->next())
    {
        if (flags & EVAL)
        {
            gsMatrix<> eval1, eval2;//, eval3;

            basis.eval_into(quNodes, eval1);
            basis.fastEval_into(quNodes, eval2);
            //basis.slowEval_into(quNodes, eval3);

            gsMatrix<> diff = eval1 - eval2;
            const real_t norm = diff.norm();
            if (maxEval < norm)
            {
                maxEval = norm;
            }

            checkTolerance(eval1, eval2, tol);
            //checkTolerance(eval1, eval3, tol);
        }

        if (flags & DERIV)
        {
            gsMatrix<> deriv1, deriv2;
            basis.deriv_into(quNodes, deriv1);
            basis.fastDeriv_into(quNodes, deriv2);

            gsMatrix<> diff = deriv1 - deriv2;
            const real_t norm = diff.norm();
            if (maxDeriv < norm)
            {
                maxDeriv = norm;
            }

            checkTolerance(deriv1, deriv2, tol);
        }

        if (flags & DERIV2)
        {
            gsMatrix<> deriv1, deriv2;
            basis.deriv2_into(quNodes, deriv1);
            basis.fastDeriv2_into(quNodes, deriv2);

            gsMatrix<> diff = deriv1 - deriv2;
            const real_t norm = diff.norm();
            if (maxDeriv2 < norm)
            {
                maxDeriv2 = norm;
            }

            checkTolerance(deriv1, deriv2, tol);
        }
    }

//    gsInfo << "If there is some output above me, we have problems..."
//              << "\n";

//    gsInfo << "difference between implementations: \n"
//              << "max(norm( evalOld - evalNew )): " << maxEval << "\n"
//              << "max(norm( derivOld - derivNew )): " << maxDeriv << "\n"
//              << "max(norm( deriv2Old - deriv2New )): " << maxDeriv2 << "\n"
//              << "\n";

}

template <unsigned d>
double speedTest(gsTHBSplineBasis<d, real_t>& basis,
                 const unsigned flags)
{
    const short_t parDim = basis.dim();

    gsVector<int> numNodes(parDim);
    for (int i = 0; i != parDim; ++i)
    {
        numNodes[i] = basis.degree(i);
    }

    gsBasis<>::domainIter domIt = basis.makeDomainIterator();
    gsGaussRule<> quRule(numNodes);
    gsMatrix<> quNodes;
    gsVector<> quWeights;
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

    gsStopwatch time;
    const double start = time.stop();
    for (; domIt->good(); domIt->next())
    {
        if (flags & EVAL)
        {
            gsMatrix<> oldEval, newEval, veryOldEval;

            if (flags & OLD)
            {
                basis.eval_into(quNodes, oldEval);
            }
            else if (flags & NEW)
            {
                basis.fastEval_into(quNodes, newEval);
            }
//            else if (flags & VERY_OLD)
//            {
//                basis.slowEval_into(quNodes, veryOldEval);
//            }
        }

        if (flags & DERIV)
        {
            gsMatrix<> oldDeriv, newDeriv;

            if (flags & OLD)
            {
                basis.deriv_into(quNodes, oldDeriv);
            }
            else if (flags & NEW)
            {
                basis.fastDeriv_into(quNodes, newDeriv);
            }
        }

        if (flags & DERIV2)
        {
            gsMatrix<> oldDeriv2, newDeriv2;

            if (flags & OLD)
            {
                basis.deriv2_into(quNodes, oldDeriv2);
            }
            else if (flags & NEW)
            {
                basis.fastDeriv2_into(quNodes, newDeriv2);
            }
        }
    }
    const double end = time.stop();

    return end - start;
}


int main(int argc, char* argv[])
{
    std::string input("");
    bool eval = false;
    bool deriv = false;
    bool deriv2 = false;
    
    gsCmdLine cmd("Optimization");
    cmd.addPlainString("filename", "File containing the input mesh", input);
    cmd.addSwitch("eval", "Test evaluation.", eval);
    cmd.addSwitch("deriv", "Test first derivatives.", deriv);
    cmd.addSwitch("deriv2", "Test second derivatives.", deriv2);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "------------------------------------------------------------"
                 "\nInput Arguments: \n\n"
                 "Input: " << input << "\n\n"
                 "Eval flag: " << eval << "\n\n"
                 "Deriv flag: " << deriv << "\n\n"
                 "Deriv2 flag: " << deriv2 << "\n\n"
                 "------------------------------------------------------------"
                 "\n\n";

    unsigned flags = NONE;
    if (eval)
    {
        flags |= EVAL;
    }

    if (deriv)
    {
        flags |= DERIV;
    }

    if (deriv2)
    {
        flags |= DERIV2;
    }

    if (input != "")
    {
        double oldTime = 0;
        double newTime = 0;
        // double veryOldTime = 0;

        gsFileData<> data(input);

        if ( data.has< gsTHBSplineBasis<1, real_t> >() )
        {
            gsTHBSplineBasis<1, real_t>::uPtr thbs;
            thbs = data.getFirst< gsTHBSplineBasis<1, real_t> >();

            oldTime = speedTest<1>(*thbs, flags | OLD);
            newTime = speedTest<1>(*thbs, flags | NEW);
            //veryOldTime = speedTest<1>(*thbs, flags | VERY_OLD);
        }
        else if ( data.has< gsTHBSplineBasis<2, real_t> >() )
        {
            gsInfo << "2d" << "\n";

            gsTHBSplineBasis<2, real_t>::uPtr thbs;
            thbs = data.getFirst< gsTHBSplineBasis<2, real_t> >();

            oldTime = speedTest<2>(*thbs, flags | OLD);
            newTime = speedTest<2>(*thbs, flags | NEW);
            //veryOldTime = speedTest<2>(*thbs, flags | VERY_OLD);

        }
        else if ( data.has< gsTHBSplineBasis<3, real_t> >() )
        {
            gsInfo << "3d" << "\n";

            gsTHBSplineBasis<3, real_t>::uPtr thbs;
            thbs = data.getFirst< gsTHBSplineBasis<3, real_t> >();

            oldTime = speedTest<3>(*thbs, flags | OLD);
            newTime = speedTest<3>(*thbs, flags | NEW);
            //veryOldTime = speedTest<3>(*thbs, flags | VERY_OLD);
        }

        gsInfo << "Times: \n"
                  << "       old: " << oldTime << "\n"
                  << "       new: " << newTime << "\n"
                  // << "  very old: " << veryOldTime << "\n"
                  << "\n";

    }
    else
    {


        const char* inputs[] = {
            "basis_thbs_01.xml",
            "basis_thbs_05.xml",
//    #ifdef NDEBUG
//            "basis_thbs_12.xml",
//    #endif
            "basis_thbs_03.xml",
//            "basis_thbs_14.xml",
//            "basis_thbs_15.xml",
            "basis_thbs_08.xml",
            "basis_thbs_04.xml",
//          "basis_thbs_1D.xml",
//          "basis2d/thbs_basis_02.xml", // periodic
//          "basis2d/thbsBasisDeg3CenterRef.xml", // debug this...
            "basis_thbs_07.xml",
            "basis_thbs_02.xml",
//          "basis_thbs_11.xml", // to big...
//    #ifdef NDEBUG
//            "basis_thbs_09.xml",
//    #endif
            "basis_thbs.xml",
    #ifdef NDEBUG
//          "basis_thbs_10.xml",
            "basis3d/basis_thbs_03.xml",
    #endif
            "basis3d/basis_thbs_01.xml",
//            "basis3d/basis_thbs_02.xml", // periodic
//          "basisBug.xml",
            "basis_thbs_deriv2_test.xml",
            "basis_thbs_06.xml",
//            "jaka.xml",
            "basis1d/thbs1_basis_03.xml",
//          "basis1d/thbs1_basis_02.xml", // periodic
            "basis1d/thbs1_basis_01.xml",
            "basis1d/thbs1_basis_04.xml",
            "basis1d/thbs1_basis_05.xml",
//            "basis_thbs_13.xml"//,
            };

        gsInfo << "Now we will compare the old and new evaluation on "
                     "several files.\n"
                  << "If there is no additional output, then the evaluations "
                     "are the same.\n"
                  << "\n";

        for (unsigned i = 0; i < sizeof(inputs) / sizeof(inputs[0]); ++i)
        {
            std::string fileName(inputs[i]);

            gsInfo << "Compare evaluation: " << fileName << "\n";

            gsFileData<> data(fileName);
            if ( data.has< gsTHBSplineBasis<1, real_t> >() )
            {
                gsInfo << "1d" << "\n";

                gsTHBSplineBasis<1, real_t>::uPtr thbs;
                thbs = data.getFirst< gsTHBSplineBasis<1, real_t> >();

                compareEval<1>(*thbs, EVAL | DERIV | DERIV2);
            }
            else if ( data.has< gsTHBSplineBasis<2, real_t> >() )
            {
                gsInfo << "2d" << "\n";

                gsTHBSplineBasis<2, real_t>::uPtr thbs;
                thbs = data.getFirst< gsTHBSplineBasis<2, real_t> >();

                compareEval<2>(*thbs, EVAL | DERIV | DERIV2);
            }
            else if ( data.has< gsTHBSplineBasis<3, real_t> >() )
            {
                gsInfo << "3d" << "\n";

                gsTHBSplineBasis<3, real_t>::uPtr thbs;
                thbs = data.getFirst< gsTHBSplineBasis<3, real_t> >();

                compareEval<3>(*thbs, EVAL | DERIV | DERIV2);
            }
        }
    }

    return 0;
}



/* saving something for future reference, maybe I will need this

void slowEval_into(const gsMatrix<>& u,
                   gsMatrix<>& result)
{
    gsMatrix<index_t> indices;
    this->active_into(u, indices);

    result.setZero(indices.rows(), u.cols());

    for (int pt = 0; pt != u.cols(); pt++)
    {
        for (int ind = 0; ind != indices.rows(); ind++)
        {
            unsigned index = indices(ind, pt);
            if (ind != 0 && index == 0)
                break;


            unsigned lvl = static_cast<unsigned>(this->levelOf(index));
            unsigned tenIndex = this->flatTensorIndexOf(index, lvl);

            gsMatrix<unsigned, d, 2> elementInd(d, 2);
            this->m_bases[lvl]->elementSupport_into(tenIndex, elementInd);

            gsVector<unsigned, d> low = elementInd.col(0);
            gsVector<unsigned, d> high = elementInd.col(1);

            unsigned clevel = this->m_tree.query4(low, high, lvl);

            if (lvl != clevel)
            {
                this->m_tree.computeFinestIndex(low, lvl, low);
                this->m_tree.computeFinestIndex(high, lvl, high);


                const gsTensorBSplineBasis<d, real_t, gsCompactKnotVector<> >& base =
                    *this->m_bases[clevel];

                gsSparseVector<> pres(base.size());
                _representBasisFunction2(index, clevel, low, high, pres);

                gsMatrix<> res;
                gsTensorDeboor<d, real_t, gsCompactKnotVector<>, gsSparseVector<> >
                    (u.col(pt), base, pres, res);
                result(ind, pt) = res(0, 0);

            }
            else
            {

                gsMatrix<> res;
                this->m_bases[lvl]->evalSingle_into(tenIndex, u.col(pt), res);
                result(ind, pt) = res(0, 0);
            }

        }
    }
}


void _representBasisFunction2(
        const unsigned j,
        const unsigned pres_level,
        const gsVector<unsigned, d>& finest_low,
        const gsVector<unsigned, d>& finest_high,
        gsSparseVector<>& presentation)
{
    const unsigned cur_level = this->levelOf(j);

    // actual size of the coefficients
    gsVector<unsigned, d> act_size_of_coefs(d);
    act_size_of_coefs.fill(1);

    // number of new coefficients
    unsigned nmb_of_coefs = _updateSizeOfCoefs(cur_level, pres_level,
                                         finest_low, finest_high,
                                         act_size_of_coefs);


    gsVector<> one(1);
    one(0) = 1.0;
    gsMatrix<> coefs(nmb_of_coefs, 1);
    coefs.fill(0);
    coefs.row(0) = one;

    // vector of the numbers of the coefficients (in each dimension)
    // stored in coefs
    gsVector<unsigned, d> vec_nmb_of_coefs(d);
    vec_nmb_of_coefs.fill(1);

    unsigned tensor_index = this->flatTensorIndexOf(j, cur_level);

    // B-Spline vector tensor index
    gsVector<unsigned, d> bspl_vec_ti =
            this->m_bases[cur_level]->tensorIndex(tensor_index);


    // we need to separately save knot vectors because we will modify
    // them, when we proceed from one level on another
    std::vector<gsCompactKnotVector<> > vector_of_kv(d);

    // size of the coefficients that are affected in individual iteration
    gsVector<unsigned, d> cur_size_of_coefs(d);
    cur_size_of_coefs.fill(1);

    for (unsigned level = cur_level; level < pres_level; ++level)
    {
        _updateSizeOfCoefs(level, level + 1, finest_low,
                           finest_high, cur_size_of_coefs);

        // index of a support of the j-th basis function (l_low, l_high
        // on level, and l1_high, l1_low on level + 1)
        gsVector<unsigned, d> clow, chigh, fhigh, flow;

        this->m_tree.computeLevelIndex(finest_low, level, clow);
        this->m_tree.computeLevelIndex(finest_high, level, chigh);

        this->m_tree.computeLevelIndex(finest_low, level + 1, flow);
        this->m_tree.computeLevelIndex(finest_high, level + 1, fhigh);


        for (unsigned dim = 0; dim < d; ++dim)
        {

            std::vector<real_t> knots;
            gsCompactKnotVector<>& ckv =
                    this->m_bases[level]->component(dim).knots();
            gsCompactKnotVector<>& fkv =
                    this->m_bases[level + 1]->component(dim).knots();


            if (level == cur_level)
                vector_of_kv[dim] = ckv;

            this->_differenceBetweenKnotVectors(ckv, clow[dim], chigh[dim],
                                          fkv, flow[dim], fhigh[dim],
                                          knots);

            gsTensorBoehmRefineLocal<d,
                                    gsCompactKnotVector<>,
                                    gsMatrix<>,
                                    std::vector<real_t>::const_iterator>
                (vector_of_kv[dim], bspl_vec_ti[dim], coefs, vec_nmb_of_coefs,
                 act_size_of_coefs,
                 cur_size_of_coefs, dim, knots.begin(), knots.end(),
                 true);
        }

        _truncate(coefs, act_size_of_coefs, cur_size_of_coefs,
                  level + 1, bspl_vec_ti, cur_level, finest_low);
    }

    _saveNewBasisFunPresentation2(coefs, act_size_of_coefs,
                                 j, pres_level, finest_low, presentation);
}


void _saveNewBasisFunPresentation2(
        const gsMatrix<>& coefs,
        const gsVector<unsigned, d>& act_size_of_coefs,
        const unsigned j,
        const unsigned pres_level,
        const gsVector<unsigned, d>& finest_low,
        gsSparseVector<>& presentation)
{

    const unsigned level = this->levelOf(j);
    const unsigned tensor_index = this->flatTensorIndexOf(j, level);

    gsVector<unsigned, d> bspl_vec_ti =
            this->m_bases[level]->tensorIndex(tensor_index);


    // finer tensor index
    const unsigned f_ten_index = _basisFunIndexOnLevel(bspl_vec_ti, level,
                                  finest_low, pres_level);


    gsVector<unsigned, d> act_coefs_strides(d);
    bspline::buildCoeffsStrides<d>(act_size_of_coefs, act_coefs_strides);


    gsVector<unsigned, d> position(d);
    position.fill(0);


    gsVector<unsigned, d> first_point(position);
    gsVector<unsigned, d> last_point(d);
    bspline::getLastIndexLocal<d>(act_size_of_coefs, last_point);

    presentation.reserve(coefs.rows());

    do
    {
        // ten_index - (tensor) index of a bspline function with respect to
        //             the coef at position "position"
        // coef_index - (local) index of a ceofficient at position

        unsigned ten_index = f_ten_index;
        for (unsigned dim = 0; dim < d; dim++)
        {

            ten_index += position(dim) *
                    this->m_bases[pres_level]->stride(static_cast<int>(dim));
        }

        unsigned coef_index = bspline::getIndex<d>(act_coefs_strides, position);

        if (coefs(coef_index) != 0)
            presentation(ten_index) = coefs(coef_index);
            //this->m_presentation[j](ten_index) = coefs(coef_index);


    } while(nextCubePoint<gsVector<unsigned, d> > (position, first_point,
                                                   last_point));
}

*/
