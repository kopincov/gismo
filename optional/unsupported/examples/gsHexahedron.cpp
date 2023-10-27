// gsHexahedron.cpp
// Experiment: trying to parameterize the internal of gsVolumeBlock
//
// Author: Jaka Speh
//

// =============================================================================
//
// THIS IS WORK IN PROGRESS
//
// I am just saving my work
// =============================================================================





#include <iostream>
#include <vector>
#include <fstream>

#include <algorithm>

#include <gismo.h>


using namespace gismo;


/// Function returns pointers to all hexahedra in given solid
///
/// Hexahedron is volume with 6 faces and each face has 4 edges.
std::vector<gsVolumeBlock<>*> getHexahedra(gsSolid<>& solid)
{
    std::vector<gsVolumeBlock<>*> hexahedra;

    gsVolumeBlock<>* volume;

    for (int idVolume = 0; idVolume != solid.nVolumes(); idVolume++)
    {
        volume = solid.volume[idVolume];

        if (volume->isHexahedron())
        {
            hexahedra.push_back(volume);
        }
    }

    return hexahedra;
}


template < unsigned d>
class gsHierarchicalRefinement: public gsFitting<real_t>
{
public:
    gsHierarchicalRefinement():
        gsFitting<real_t> (),
        m_lambda(0),
        m_refinePercent(0),
        m_extension(0),
        m_errorType(0),
        m_aligned(false)
    { }

    gsHierarchicalRefinement(const gsMatrix<>& parameters,
                             const gsMatrix<>& points,
                             gsTHBSplineBasis<d, real_t>& basis,
                             std::vector<int>& extension,
                             const real_t refinePercent = 0.8,
                             const real_t lambda = 1e-9
                             ):
        gsFitting<real_t> (parameters, points, basis),
        m_lambda(lambda),
        m_refinePercent(refinePercent),
        m_extension(extension),
        m_errorType(0),
        m_aligned(true)
    {
        GISMO_ASSERT(m_extension.size() == d,
                     "Extension does not have correct size");
        GISMO_ASSERT(0 <= m_refinePercent && m_refinePercent <= 1,
                     "Refine percent must be between 0 and 1");

        this->compute(lambda);
    }


    /// Refines at areas with big error and does the fitting again.
    ///
    /// \param iterations number of times we refine and perform a fitting
    /// \param tolerance maximum error should be under the tolerance
    /// \param errorTreshold if negative, we refine m_refinePercent of errors
    void refine(const int iterations,
                const real_t tolerance,
                real_t errorTreshold,
                const bool output = false,
                std::ostream& out = gsInfo)
    {
        gsTHBSplineBasis<d, real_t>* basis =
                static_cast<gsTHBSplineBasis<d, real_t>*> (this->m_basis);

//        gsInfo << "basis: " << *basis << "\n"
//                  << "level: " << *(basis->m_bases[0]) << "\n";


        for (int iter = 0; iter < iterations; iter++)
        {
            std::vector<real_t> errors;
            this->get_Error(errors, m_errorType);
            real_t max = *std::max_element(errors.begin(), errors.end());

            if (output)
            {
                out << "max error: " << max << "\n"
                       "tolerance: " << tolerance << "\n";
            }

//            gsInfo << "Errors: \n";
//            for (size_t index = 0; index < errors.size(); index++)
//            {
//                gsInfo << errors[index] << "\n";
//            }




            if (tolerance < max)
            {
                if (errorTreshold < 0)
                {
                    errorTreshold = setErrorTreshold(errors);
                }

                std::vector<int> cells;
                selectCells(cells, errors, errorTreshold);

//                // delete
//                gsInfo << "Cells: \n";
//                for (size_t index = 0; index < cells.size(); index++)
//                {
//                    if (index % 3 == 0)
//                        gsInfo << "\n";

//                    gsInfo << cells[index] << " ";
//                }
//                gsInfo << "\n";
//                // delete


                gsVector<int> levels(cells.size() / d);
                getLevels(levels, cells);

//                gsInfo << "levels: \n" << levels << "\n";

                // boxes are defined as
                // [level,
                //  lowCorner_1, ..., lowCorner_d,
                //  uppCorner_1, ..., uppCorner_d, ... ]
                std::vector<unsigned> boxes (levels.size() * (2 * d + 1));

//                gsInfo << "boxes.size(): " << boxes.size() << "\n";

                getBoxes(boxes, cells, levels);

//                for (size_t index = 0; index < boxes.size(); index++)
//                {
//                    if (index % 7 == 0)
//                        gsInfo << "\n";

//                    gsInfo << boxes[index] << " ";
//                }


                basis->refineElements(boxes);

                if (output)
                {
                    out << "Boxes are inserted..." << "\n";
                }

//                gsInfo << "boxes are inserted..." << "\n";
//                gsInfo << "basis: " << *basis << "\n"
//                          << "level: " << *(basis->m_bases[1]) << "\n";

//                gsTHBSplineBasis<d, real_t>* myBasis =
//                        static_cast<gsTHBSplineBasis<d, real_t>*> (this->m_basis);

//                gsFileData<real_t> newdata;
//                newdata << *myBasis;

//                std::string xmlFile("debugBasis.xml");
//                newdata.dump(xmlFile);


                this->compute(m_lambda);

//                gsInfo << "we recomputed the fitting..." << "\n";

            }
            else
            {
                if (output)
                {
                    out << "Errors are under tolerance!\n"
                           "End fitting...\n";
                }
            }


        }
    }

    bool isAligned() const
    {
        return m_aligned;
    }

    void setAligned(const bool aligned)
    {
        m_aligned = aligned;
    }

    int getErrorType() const
    {
        return m_errorType;
    }

    void setErrorType(const int errorType)
    {
        m_errorType = errorType;
    }

    real_t getMaxError() const
    {
        std::vector<real_t> errors;
        this->get_Error(errors, m_errorType);
        return *std::max_element(errors.begin(), errors.end());
    }

    void isErrorOnEdges(real_t treshHold = 0.0001) const
    {

        std::vector<real_t> errors;
        this->get_Error(errors, m_errorType);

        int errorOnEdges = 0;
        real_t comError = 0.0;

        for (size_t col = 0; col < errors.size(); col++)
        {
            if (treshHold < errors[col])
            {
                unsigned isBorder = 0;
                for (int row = 0; row < 3; row++)
                {
                    if (this->m_param_values(row, col) == 0.0 ||
                            this->m_param_values(row, col) == 1.0)
                    {
                        isBorder++;
                    }
                }

                if (2 <= isBorder)
                {
                    comError += errors[col];
                    errorOnEdges++;
                }
            }
        }

        gsInfo << "errorOnEdges: " << errorOnEdges << "\n"
                     "total: " << errors.size() << "\n"
                     "percent: " << errorOnEdges * 1.0 / errors.size() << "\n"
                     "comError: " << comError << "\n";
    }

    gsMatrix<> getErrorPoints(real_t percent)
    {
        std::vector<real_t> errors;
        this->get_Error(errors, m_errorType);

        real_t tmp = m_refinePercent;
        m_refinePercent = percent;

        real_t treshold = setErrorTreshold(errors);

        m_refinePercent = tmp;

        int counter = 0;
        for (size_t er = 0; er < errors.size(); er++)
        {
            if (treshold < errors[er])
                counter++;
        }

        gsMatrix<> result;
        result.resize(this->m_points.cols(), counter);

        int index = 0;
        for (size_t er = 0; er < errors.size(); er++)
        {
            if (treshold < errors[er])
            {
                result.col(index) = this->m_points.row(er);
                index++;
            }

        }

        return result;
    }

    gsMatrix<> getGoodPoints(const real_t treshold = 1e-6) const
    {
        std::vector<real_t> errors;
        this->get_Error(errors, m_errorType);

        int counter = 0;
        for (size_t er = 0; er < errors.size(); er++)
        {
            if (errors[er] < treshold)
                counter++;
        }

        gsMatrix<> result;
        result.resize(this->m_points.cols(), counter);

        int index = 0;
        for (size_t er = 0; er < errors.size(); er++)
        {
            if (errors[er] < treshold)
            {
                result.col(index) = this->m_points.row(er);
                index++;
            }

        }

        return result;
    }

    gsMatrix<> getBadPoints(const real_t treshold = 1e-6) const
    {
        std::vector<real_t> errors;
        this->get_Error(errors, m_errorType);

        int counter = 0;
        for (size_t er = 0; er < errors.size(); er++)
        {
            if (treshold <= errors[er])
                counter++;
        }

        gsMatrix<> result;
        result.resize(this->m_points.cols(), counter);

        int index = 0;
        for (size_t er = 0; er < errors.size(); er++)
        {
            if (treshold <= errors[er])
            {
                result.col(index) = this->m_points.row(er);
                index++;
            }

        }

        return result;
    }

private:

    void getBoxes(std::vector<unsigned>& boxes,
                  const std::vector<int>& cells,
                  const gsVector<int>& levels)
    {
        gsTHBSplineBasis<d, real_t>* basis = static_cast<gsTHBSplineBasis<d, real_t>*>
                                                    (this->m_basis);

        const int maxLevel = basis->maxLevel();

//        gsInfo << "maxLevel: " << maxLevel << "\n";

        // for loop over levels of all cells
        for (int ell = 0; ell < levels.rows(); ell++)
        {

//            gsInfo << "=======================================\n"
//                         "box " << ell << "\n";

            const int lvl = levels[ell] + 1;

            const int offset = ell * (2 * d + 1);
            boxes[offset] = lvl;

            const gsTensorBSplineBasis<d,real_t>& tensorBasis =
                    *(basis->getBases()[lvl]);

            for (unsigned dim = 0; dim < d; dim++)
            {
                // find index

//                gsInfo << "level: " << lvl << "\n";

                const std::vector<real_t>& finestKnots =
                        basis->getBases()[maxLevel]->component(dim).knots().unique();

                // we add 1 because uniquefindspan returns previous knot if an
                // input is a knot in knot sequence
                const real_t knt = finestKnots[cells[ell + dim] + 1];

                const gsKnotVector<real_t>& ckv = tensorBasis.component(dim).knots();

                const int indx = ckv.uFind(knt).uIndex();

                // apply extensions

                int low = indx - m_extension[dim];
                int upp = indx + m_extension[dim] + 1;


                // check boundaries

                if (low < 0)
                {
                    low = 0;
                }

                // max knot index in level lvl
                const unsigned maxIndx = ckv.unique().size() - 1;

                if (maxIndx < static_cast<unsigned>(upp))
                {
                    upp = static_cast<int>(maxIndx);
                }


                if (m_aligned)
                {
                    const gsKnotVector<real_t>& coarseKV =
                            basis->getBases()[lvl - 1]->component(dim).knots();


//                    gsInfo << "cell: " << cells[ell + dim] << "\n"
//                              << "low: " << low << "    "
//                              << "upp: " << upp << "\n";

                    const real_t lowKnot = ckv.unique()[low];
                    const real_t uppKnot = ckv.unique()[upp];


                    // this works only for dyadic refinement

                    if (!coarseKV.has(lowKnot))
                    {
                        low--;
                        //low = coarseKV.Uniquefindspan(lowKnot) * 2;
                    }

                    if (!coarseKV.has(uppKnot))
                    {
                        //upp = (coarseKV.Uniquefindspan(uppKnot) + 1) * 2;
                        upp++;
                    }
                }

//                gsInfo << "cIndx: " << ell << "  dim: " << dim << "\n"
//                          << "  offset: " << offset << "  low: " << low
//                          << "  upp: " << upp << "\n" << "\n";

                boxes[offset + 1 + dim] = static_cast<unsigned>(low);
                boxes[offset + 1 + d + dim] = static_cast<unsigned>(upp);
            }

        } // for cells
    }


    void getLevels(gsVector<int>& levels,
                   const std::vector<int>& cells)
    {
        gsTHBSplineBasis<d, real_t>* basis = static_cast<gsTHBSplineBasis<d, real_t>*>
                                                    (this->m_basis);

        gsVector<unsigned, d> lower(d);
        gsVector<unsigned, d> upper(d);

        for (size_t cIndx = 0; cIndx < cells.size(); cIndx += d)
        {
            for (unsigned dim = 0; dim < d; dim++)
            {
                lower(dim) = static_cast<unsigned>(cells[cIndx + dim]);
                upper(dim) = lower(dim) + 1;
            }

            const int lvl = basis->tree().query3(lower,
                                                 upper,
                                                 basis->getBases().size() - 1);

            levels(cIndx / d) = lvl;
        }
    }


    void selectCells(std::vector<int>& cells,
                     const std::vector<real_t>& errors,
                     const real_t errorTreshold)
    {
        gsTHBSplineBasis<d, real_t>* basis = static_cast<gsTHBSplineBasis<d, real_t>*>
                                                       (this->m_basis);

        // temporarily save cell
        gsVector<int, d> tmp;
        const int maxLevel = basis->maxLevel();
        gsTensorBSplineBasis<d,real_t> tensorBasis =
                *(basis->getBases()[maxLevel]);


        for (size_t eIndx = 0; eIndx < errors.size(); eIndx++)
        {

            if (errorTreshold <= errors[eIndx])
            {

                // find the cell at smallest level

                for (unsigned dim = 0; dim < d; dim++)
                {
                    tmp[dim] = tensorBasis.component(dim).knots().Uniquefindspan(
                            this->m_param_values(dim, eIndx));
                }

                // check if the cell is already found

                bool insert = true;
                for (size_t cIndx = 0; cIndx < cells.size(); cIndx += d)
                {
                    unsigned matches = 0;
                    for (unsigned dim = 0; dim < d; dim++)
                    {
                        if (cells[cIndx + dim] == tmp[dim])
                        {
                            matches++;
                        }
                    }

                    if (matches == d)
                    {
                        insert = false;
                        break;
                    }
                }

                if (insert)
                {
                    for (unsigned dim = 0; dim < d; dim++)
                    {
                        cells.push_back(tmp[dim]);
                    }
                }

            } // if errors greater than errorTreshold

        } // for errors
    }


    real_t setErrorTreshold(std::vector<real_t> errors)
    {
        // we must make a copy of errors vector, because we must not change
        // the original errors

        std::sort(errors.begin(), errors.end());
        const size_t index =
                static_cast<size_t> (errors.size() * (1.0 - m_refinePercent));

        return errors[index];
    }


private:
    /// smoothing parameter
    real_t m_lambda;

    /// how many percent of errors to refine, it should be between [0, 1]
    real_t m_refinePercent;

    /// how many neighbouring cells we refine additionally (in every direction)
    std::vector<int> m_extension;

    /// error type
    int m_errorType;

    /// True if boxes in level ell are aligned with knots of level ell - 1
    bool m_aligned;
};




int main()
{
    gsInfo << "In development..." << "\n";

    return 0;

    /* // put // (a double-slash) in the beginning of the line to uncomment
    // parsing arguments

    std::string filename("solids/gsSolidSP2.xml");
    int n = 10;
    real_t lambda = 1e-9;
    int iterations = 1;
    real_t tolerance = 1e-6;
    int degree = 3;
    int internalKnots = 3;
    std::string output("hex");
    
    gsCmdLine cmd("Reading gsSolid");
    cmd.addString("f", "filename", "Input filename", filename);
    cmd.addInt("n", "numPts", "Number of sample points", n);
    cmd.addReal("l", "lambda", "Smoothing parameter", lambda);
    cmd.addInt("i", "iterations", "Number of iterations", iterations);
    cmd.addReal("t", "tolerance", "Maximum error should be under the tolerance",
                tolerance);
    cmd.addInt("d", "degree", "Degree of the B-Splines", degree);
    cmd.addInt("k", "numKnots", "Number of interior knots", internalKnots);
    cmd.addString("o", "output", "Output prefix", output);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "------------------------------------------------------------"
                 "\nInput Arguments: \n\n"
                 "Input filename: " << filename << "\n\n"
                 "Ouput prefix: " << output << "\n\n"
                 "Number of sample points: " << n << "\n\n"
                 "Lambda: " << lambda << "\n\n"
                 "Iterations: " << iterations << "\n\n"
                 "Tolerance: " << tolerance << "\n\n"
                 "Degree: " << degree << "\n\n"
                 "Number of internal knots: " << internalKnots << "\n"
                 "------------------------------------------------------------"
                 "\n\n";



    // reading file

    gsInfo << "Reading file..." << "\n";

    gsFileData<> data(filename);
    gsSolid<>* pSolid = NULL;

    if (data.has< gsSolid<> >())
    {
        pSolid = data.getFirst< gsSolid<> >();
    }

    // -------------------------------------------------------------------------
    // program

    std::vector< gsVolumeBlock<>* > hexahedra;
    hexahedra = getHexahedra(*pSolid);

    gsInfo << "Got " << hexahedra.size() << " hexahedra.\n" << "\n";

    for (size_t h = 0; h < hexahedra.size(); h++)
    {


        gsInfo << "  ===================================================\n"
                     "  Hexahedra " << h + 1 << " / " << hexahedra.size()
                  << "\n";

        std::string localPrefix = output + util::to_string(h);

        gsVolumeBlock<> volumeBlock = *(hexahedra[h]);

        std::string input = localPrefix + "_inputVolume";
        gsWriteParaview<>(volumeBlock, input);

        gsInfo << "  Getting boundary points...\n" << "\n";

        gsMatrix<> points;
        gsMatrix<> params;
        volumeBlock.getUniformBoundaryPoints(n, params, points);

        gsInfo << "  Got boundary points...\n" << "\n";

        // ---------------------------------------------------------------------
        // delete min, max, points

        gsInfo << "size: " << points.rows() << " x " << points.cols() << "\n";

        gsInfo << "min x: " << points.row(0).minCoeff() << "\n"
                  << "min y: " << points.row(1).minCoeff() << "\n"
                  << "min z: " << points.row(2).minCoeff() << "\n"
                  << "\n";

        gsInfo << "max x: " << points.row(0).maxCoeff() << "\n"
                  << "max y: " << points.row(1).maxCoeff() << "\n"
                  << "max z: " << points.row(2).maxCoeff() << "\n"
                  << "\n";

        // delete
        // ---------------------------------------------------------------------


//        for (int col = 0; col < points.cols(); col++)
//        {
//            gsInfo << "parameter: " << params.col(col).transpose() << "\n"
//                      << "points: " << points.col(col).transpose() << "\n"
//                      << "\n";
//        }


        std::string pts = localPrefix + "_points";
        gsWriteParaviewPoints<>(points, pts);


        gsKnotVector<> kv1(0, 1, internalKnots, degree + 1, 1);
        gsKnotVector<> kv2(0, 1, internalKnots, degree + 1, 1);
        gsKnotVector<> kv3(0, 1, internalKnots, degree + 1, 1);

        gsTensorBSplineBasis<3> tensorBasis(kv1, kv2, kv3);
        gsTHBSplineBasis<3>  thbs(tensorBasis);

        gsInfo << "  Build thbs object...\n";


        std::vector<int> ext;
        ext.push_back(2);
        ext.push_back(2);
        ext.push_back(2);

        gsInfo << "  Starting fitting...\n" << "\n";

        gsHierarchicalRefinement<3> ref(params,
                                                points,
                                                thbs,
                                                ext);

        real_t error = -1;
        ref.computeApproxError(error);

        gsInfo << "  Tensor B-Spline fitting complete. \n"
                  << "  Maximum error: " << ref.getMaxError() << "\n"
                     "  Aprox. error: " << error / points.cols() << "\n" << "\n";

        ref.isErrorOnEdges();

        gsTHBSpline<3, real_t>* tensorBSpline;
        tensorBSpline = static_cast<gsTHBSpline<3, real_t>*> (ref.result());

//        gsTHBSpline<3, real_t>* tensorToThb = dynamic_cast<gsTHBSpline<3, real_t>*>
//                (tensorBSpline);

        gsFileData<> tensorData;
        tensorData << *tensorBSpline;

        std::string xmlTensorFile = localPrefix + "_bspline.xml";
        tensorData.dump(xmlTensorFile);

        std::string paraTensorFile = localPrefix + "_bspline";

        gsWriteParaview(*tensorBSpline, paraTensorFile);

        // ---------------------------------------------------------------------
        // points, delete
        std::string paraErrorPoints = localPrefix + "_bsplineBadPoints";

        gsMatrix<> tensorBadPoints = ref.getBadPoints();
        gsMatrix<> tensorGoodPoints = ref.getGoodPoints();
        gsMatrix<>& tensorCoefs = tensorBSpline->coefs();

        gsWriteParaviewPoints(tensorBadPoints, paraErrorPoints);
        paraErrorPoints = localPrefix + "_bsplineGoodPoints";

        gsWriteParaviewPoints(tensorGoodPoints, paraErrorPoints);

        paraErrorPoints = localPrefix + "_bsplineControlPoints";

        tensorCoefs.transposeInPlace();
        gsWriteParaviewPoints(tensorCoefs, paraErrorPoints);

        // points, delete
        // ---------------------------------------------------------------------

        for (int it = 0; it < iterations; it++)
        {
            gsInfo << "    -----------------------------------------------\n"
                      << "    Refinement " << it + 1 << " / " << iterations
                      << "\n" << "\n";

            ref.refine(1, tolerance, tolerance);

            real_t errorApprox = -1;
            ref.computeApproxError(errorApprox);

            gsInfo << "    End of refinement...  \n"
                         "    Maximum error: " << ref.getMaxError() << "\n"
                         "    Aprox. error: " << errorApprox / points.cols() << "\n";

            ref.isErrorOnEdges();

            gsGeometry <>* tmpThbs;
            tmpThbs = ref.result();

            gsMatrix<>& coefs  = tmpThbs->coefs();


            // -----------------------------------------------------------------
            // delete min, max, coefs

            gsInfo << "size: " << coefs.rows() << " x " << coefs.cols() << "\n";

            gsInfo << "coefs min x: " << coefs.col(0).minCoeff() << "\n"
                      << "coefs min y: " << coefs.col(1).minCoeff() << "\n"
                      << "coefs min z: " << coefs.col(2).minCoeff() << "\n"
                      << "\n";

            gsInfo << "coefs max x: " << coefs.col(0).maxCoeff() << "\n"
                      << "coefs max y: " << coefs.col(1).maxCoeff() << "\n"
                      << "coefs max z: " << coefs.col(2).maxCoeff() << "\n"
                      << "\n";

            // -----------------------------------------------------------------






            gsFileData<> thbsData;
            thbsData << *tmpThbs;

            std::string xmlThbsFile = localPrefix + "_thbs" +
                                      util::to_string(it) + ".xml";

            thbsData.dump(xmlThbsFile);

            // basis
            xmlThbsFile = localPrefix + "_thbsBasis" +
                    util::to_string(it) + ".xml";

            gsBasis<>& b = tmpThbs->basis();

            gsFileData<> thbsDataBasis;
            thbsDataBasis << b;
            thbsDataBasis.dump(xmlThbsFile);



            std::string paraThbsFile = localPrefix + "_thbs" +
                                         util::to_string(it);
            gsWriteParaview(*tmpThbs, paraThbsFile);

            // -----------------------------------------------------------------
            // points, delete
            std::string paraErrorPoints2 = localPrefix + "_thbsBadPts" +
                                          util::to_string(it);
            gsMatrix<> badPoints = ref.getBadPoints();
            gsWriteParaviewPoints(badPoints, paraErrorPoints2);

            paraErrorPoints2 = localPrefix + "_thbsGoodPoints" +
                                                     util::to_string(it);

            gsMatrix<> goodPoints = ref.getGoodPoints();
            gsWriteParaviewPoints(goodPoints, paraErrorPoints2);


            paraErrorPoints2 = localPrefix + "_thbsControlPoints" +
                    util::to_string(it);

            coefs.transposeInPlace();
            gsWriteParaviewPoints(coefs, paraErrorPoints2);

            // points, delete
            // -----------------------------------------------------------------


            gsInfo << "    Number of elements: "
                      << tmpThbs->basis().numElements() << "\n"
                      << "    Maximum possible number of elements: "
                      << std::pow(
                             static_cast<real_t>((internalKnots + 1) *
                                                 (1 << (it + 1))),
                             3.0)
                      << "\n\n"
                      << "\n";

        }

        gsInfo << "\n" << "\n";

    }

    gsInfo << "End of program." << "\n";




//    std::ofstream filePoints("IMPORTANT_POINTS.txt");
//    std::ofstream fileParams("IMPORTANT_PARAMETERS.txt");

//    filePoints << points;
//    fileParams << params;



    delete pSolid;
    return 0;
    
    //*/
}




