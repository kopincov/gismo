#include <iostream>
#include <string>
#include <bitset>

#include <parasolid_kernel.h>

#include <gsNurbs/gsBSpline.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsHBSplineBasis.h>
#include <gsHSplines/gsHFitting.h>
#include <gsHSplines/gsHFittingLvlConstrained.h>

#include <gsIO/gsIOUtils.h>
#include <gsIO/gsWriteParaview.h>

#include <gsUtils/gsStopwatch.h>

#include <gsParasolid/gsReadParasolid.h>
#include <gsParasolid/gsWriteParasolid.h>

#include <gsMSplines/gsMappedBasis.h>
#include <gsSmoothPatches/gsCompositeBSplineBasis.h>
#include <gsSmoothPatches/gsCompositeHBasis.h>

using namespace gismo;


// ================================================================================
// delete from here

struct TransformationData
{
    gsGeometry<>* geometry;
    PK_CURVE_t evaluation_curve;
    bool reverse;
    double tmax;
    double tmin;
};

void transform_vector(PK_VECTOR_t position,
                      double* derivs,
                      gsGeometry<>* geom)
{
    gsMatrix<> input(3, 1);
    input(0, 0) = position.coord[0]; // everything is in meters now
    input(1, 0) = position.coord[1];
    input(2, 0) = position.coord[2];

    gsMatrix<> result;
    geom->eval_into(input, result);

    derivs[0] = position.coord[0] + result(0, 0); // everything is in meters now
    derivs[1] = position.coord[1] + result(1, 0);
    derivs[2] = position.coord[2] + result(2, 0);
}

PK_ERROR_code_t deformationEvaluator(PK_VECTOR_t position,
                                     PK_FACE_t face,
                                     PK_LOGICAL_t have_params,
                                     PK_UV_t params,
                                     PK_POINTER_t external_data,
                                     PK_VECTOR_t *const deformed_position)
{
    TransformationData* data = static_cast< TransformationData* >(external_data);
    transform_vector(position, deformed_position->coord, data->geometry);
    return 0;
}

PK_ERROR_code_t curveEvaluator(double parameter,
                               int n_deriv,
                               PK_HAND_t hand_direction,
                               PK_POINTER_t external_data,
                               double* derivs)
{
    TransformationData* data = static_cast< TransformationData* >(external_data);

    if (data->reverse)
    {
        double t = data->tmax + data->tmin - parameter;
        parameter = t;
    }

    PK_VECTOR_t position;
    PK_CURVE_eval(data->evaluation_curve, parameter, n_deriv, &position);
    transform_vector(position, derivs, data->geometry);

    return 0;
}

PK_CURVE_t transformEdge(PK_EDGE_t edge,
                         gsGeometry<>* transformation,
                         double tolerance)
{
    PK_BCURVE_create_fitted_o_t options;
    PK_BCURVE_create_fitted_o_m(options);

    PK_LOGICAL_t orientation;
    PK_CURVE_t evaluation_curve;

    {
        PK_VECTOR_t ends[2];
        PK_CLASS_t cls;
        PK_ERROR_t err = PK_EDGE_ask_geometry(edge, PK_LOGICAL_true, &evaluation_curve,
                                              &cls, ends, &(options.range), &orientation);
        if (err)
        {
            std::cout << "PK_ERROR: PK_EDGE_ask_geometry: " << err << std::endl;
        }

    }

    if (evaluation_curve == 0)
    {
        return 0;
    }

    TransformationData data;
    data.geometry = transformation;
    data.evaluation_curve = evaluation_curve;

    if (orientation == PK_LOGICAL_false)
    {
        data.reverse = true;
        data.tmin = options.range.value[0]; // are not set in every case
        data.tmax = options.range.value[1];
    }
    else
    {
        data.reverse = false;
    }

    PK_BCURVE_t curve;
    PK_BCURVE_fitted_fault_t fault;
    options.curve.type = PK_CURVE_general_user_c;
    options.curve.curve.user_curve.highest_deriv = 0;
    options.curve.curve.user_curve.dimension = 3;
    options.curve.curve.user_curve.eval_fn = curveEvaluator;
    options.curve.curve.user_curve.eval_data = static_cast<void*>(&data);
    options.curve.curve.user_curve.is_closed = PK_LOGICAL_false;
    options.curve.curve.user_curve.is_periodic = PK_LOGICAL_false;
    options.curve.curve.user_curve.n_discontinuities = 0;
    options.curve.curve.user_curve.discontinuities = 0;
    options.tolerance = tolerance;
    PK_ERROR_t err = PK_BCURVE_create_fitted(&options, &curve, &fault);
    if (err)
    {
        std::cout << "PK_ERROR: PK_BCURVE_create_fitted: " << err << std::endl;
    }

    if (fault.status != PK_BCURVE_fitted_success_c)
    {
        std::cout << "fitting status = " << fault.status << std::endl;
        std::cout << "curve fitting failed in edge transformation" << std::endl;
        // throw ("curve fitting failed in edge transformation");
    }

    return curve;
}


void transform_solid(PK_BODY_t body,
                     gsGeometry<>* transformation,
                     double tolerance)
{
    PK_FACE_t* faces;
    int nfaces;
    PK_BODY_ask_faces(body, &nfaces, &faces);

    PK_FACE_change_o_t options;
    PK_FACE_change_o_m(options);

    std::vector<PK_EDGE_t> innerEdges;
    std::vector<PK_EDGE_t> transformEdges;
    std::vector<PK_CURVE_t> transformedCurves;
    std::vector<double> tolerances;

    PK_FACE_t* edges;
    int nedges;
    PK_ERROR_t err = PK_BODY_ask_edges(body, &nedges, &edges);
    if (err)
    {
        std::cout << "Error PK_BODY_ask_edges: " << err << std::endl;
    }

    std::cout << "found " << nedges << " edges alltogether" << std::endl;


    for (int i = 0; i < nedges; i++)
    {
        int nf;
        PK_ERROR_t err = PK_EDGE_ask_faces(edges[i], &nf, 0);
        if (err)
        {
            std::cout << "Error PK_EDGE_ask_faces: " << err << std::endl;
        }

        if (nf > 1)
        {
            innerEdges.push_back(edges[i]);
        }
    }

    for (size_t i = 0; i < innerEdges.size(); i++)
    {
        PK_CURVE_t c = transformEdge(innerEdges[i], transformation, tolerance);
        if (c)
        {
            transformEdges.push_back(innerEdges[i]);
            transformedCurves.push_back(c);
            tolerances.push_back(tolerance);
        }
    }

    options.edge_data.n_edges = transformEdges.size();
    options.edge_data.edges = &(transformEdges[0]);
    options.edge_data.curves = &(transformedCurves[0]);
    options.edge_data.tolerances = &(tolerances[0]);
    options.edge_data.replace_use = PK_replace_use_yes_c;

    int mapping[nfaces];
    for (int i = 0; i < nfaces; i++)
    {
        mapping[i] = 0;
    }


    {
        TransformationData data;
        data.geometry = transformation;

        PK_TOPOL_track_r_t tracking;
        PK_TOPOL_local_r_t results;
        PK_FACE_change_t operation;

        PK_FACE_change_deform_o_t deformOptions;
        PK_FACE_change_deform_o_m(deformOptions);
        deformOptions.deform_uv = PK_deform_uv_face_box_c;
        // deformOptions.deform_uv =  PK_deform_uv_all_c;

        operation.op_type = PK_FACE_change_type_deform_c;
        operation.op_param.deform.eval_fn = deformationEvaluator;
        operation.op_param.deform.eval_data = static_cast<void*>(&data);
        operation.op_opts.deform = &deformOptions;
        options.merge_face = PK_LOGICAL_false;
        options.grow = PK_FACE_grow_no_c;
        PK_ERROR_t err = PK_FACE_change(nfaces, faces, mapping, 1, &operation, tolerance,
                                        &options, &tracking, &results);
        if (err)
        {
            std::cout << "PK_ERROR: PK_FACE_change: " << err << std::endl;
        }

        if (results.status != PK_local_status_ok_c)
        {
            std::cout << "failed: error = " << results.status << std::endl;
            std::cout << "PK_FACE_change failed" << std::endl;
            // throw "transformation failed";
        }
    }

}


// delete end
// ================================================================================

// Functions for expressing gismo surfaces as an pk_FSURF in PARASOLID.

#define FGEVSQ 2

// Index for surface positon P
inline
bool ev_P(int& n, int nu, int nv, int tri)
{
    n = 0;
    return true;
}

// Index for surface u derivative dP/du
inline
bool ev_dP_du(int& n, int nu, int nv, int tri)
{
    bool do_eval = false;
    if (nu > 0)
    {
        n = 3;
        do_eval = true;
    }
    return do_eval;
}

// Index surface 2nd u derivative d2P/du2
inline
bool ev_dP2_du2(int& n, int nu, int nv, int tri)
{
    bool do_eval = false;
    if (nu > 1)
    {
        n = 6;
        do_eval = true;
    }
    return do_eval;
}

// Index for surface v derivative dP/dv
inline
bool ev_dP_dv(int& n, int nu, int nv, int tri)
{
    int do_eval = false;
    if (nv > 0)
    {
        n = 3 * (nu + 1);
        do_eval = true;
    }
    return do_eval;
}

// Index surface 2nd v derivative d2P/dv2
inline
bool ev_dP2_dv2(int& n, int nu, int nv, int tri)
{
    bool do_eval = false;
    int triang = (tri == FGEVSQ) ? 0 : 1;;

    if (nv > 1)
    {
        n = 6 * nu + 3 * (2 - triang);
        do_eval = true;
    }
    return do_eval;
}

// Index for mixed derivative d2P/dudv
inline
int ev_dP2_dudv(int& n, int nu, int nv, int tri)
{
    bool do_eval = false;
    int triang = (tri == FGEVSQ) ? 0 : 1;

    if ((nu > 0 && nv > 0 && triang == 0) ||(nu > 1 && triang == 1 ))
    {
        n = 3 * (nu + 2);
        do_eval = true;
    }
    return do_eval;
}

// ================================================================================

void fill_matrix(gsMatrix<>& mat,
                 double* data,
                 int rows,
                 int cols)
{
    // this should be equivalent to
    // mat = gsAsMatrix<>(data,cols,rows).transpose();

    mat.resize(rows, cols);
    for (int r = 0; r != rows; r++)
    {
        for (int c = 0; c != cols; c++)
        {
            mat(r, c) = data[r * cols + c];
        }
    }
}

void fill_vector(gsVector<int>& vec,
                 int* data,
                 int length)
{
    vec = gsAsVector<int>(data, length);
}

void fill_stdvector(std::vector<unsigned>& vec,
                    int* data,
                    int length)
{
    vec.clear();
    for(int i = 0;i<length;++i)
    {
        vec.push_back(static_cast<unsigned>(data[i]));
    }
}

void fill_stdvector(std::vector<double>& vec,
                    double* data,
                    int length)
{
    vec.clear();
    for(int i = 0;i<length;++i)
    {
        vec.push_back(data[i]);
    }
}

double* make_c_array(const gsMatrix<>& matrix)
{
    // code below is equivalent to
    // double* result = new double[matrix.size()];
    // copyRange(matrix.data(), result, matrix.size());
    //
    // I prefer to use the code below, because it is independent
    // of row or column mode of Eigen.

    double* result = new double[matrix.rows() * matrix.cols()];
    for (int col = 0; col != matrix.cols(); col++)
    {
        for (int row = 0; row != matrix.rows(); row++)
        {
            result[col * matrix.rows() + row] = matrix(row, col);
        }
    }

    return result;
}

double* make_c_array(const std::vector<double>& vector)
{
    double* data = new double[vector.size()];

    for (size_t i = 0; i != vector.size(); i++)
    {
        data[i] = vector[i];
    }

    return data;
}

int* make_c_array_int(const gsVector<int>& vector)
{
    int* data = new int[vector.rows()];

    for (size_t i = 0; i != static_cast<unsigned>(vector.rows()); i++)
    {
        data[i] = vector(i);
    }

    return data;
}

std::vector<real_t> compute_errors(const gsMatrix<>& uv,
                                   const gsMatrix<>& xyz,
                                   const gsGeometry<>& geom)
{
    gsMatrix<> eval;
    geom.eval_into(uv, eval);

    std::vector<real_t> errors;
    gsVector<real_t, 3> dist;
    for (index_t col = 0; col != eval.cols(); col++)
    {
        dist = xyz.col(col) - eval.col(col);

        errors.push_back(dist.norm());
    }

    return errors;
}


real_t percent_point_below_threshold(const std::vector<real_t>& errors,
                                     const real_t threshold)
{
    int numErrors = errors.size();

    int belowThreshold = 0;
    for (size_t i = 0; i != errors.size(); i++)
    {
        if (errors[i] <= threshold)
        {
            belowThreshold++;
        }
    }

    return belowThreshold * 1.0 / numErrors;
}

gsBoxTopology get_periodic_topology()
{
    gsBoxTopology topology(2, 1);
    gsVector<index_t> dirMap(2);
    dirMap << 0, 1;
    gsVector<bool> dirOrientation(2);
    dirOrientation << true, true;

    patchSide ps1(0, boundary::left);
    patchSide ps2(0, boundary::right);

    boundaryInterface boundInter(ps1, ps2, dirMap, dirOrientation);
    topology.addInterface(boundInter);

    return topology;
}


template <typename HierSpline, unsigned d>
void refine_basis_function(void* pointer,
                           int basisFun)
{
    HierSpline* spline = static_cast< HierSpline* >(pointer);

    gsMatrix<unsigned, d, 2> support;
    spline->basis().elementSupport_into(basisFun, support);

    const unsigned level = spline->basis().levelOf(basisFun) + 1;

    std::vector<unsigned> box;
    box.push_back(level);
    for (unsigned corner = 0; corner != 2; corner++)
    {
        for (unsigned dim = 0; dim != d; dim++)
        {
            box.push_back(support(dim, corner) << 1); // assumes dyadic refinement
        }
    }

    spline->basis().refineElements_withCoefs(spline->coefs(), box);
}

template <unsigned d>
int get_index(void* pointer,
              const gsVector<int>& tensorIndex)
{
    typedef typename gsBSplineTraits<d, real_t>::Basis Basis;

    Basis* b = static_cast< Basis* >(pointer);

    gsVector<unsigned, d> vec;
    for (int row = 0; row != vec.rows(); row++)
    {
        vec(row) = static_cast<unsigned>(tensorIndex(row));
    }

    return static_cast<int>(b->index(vec));
}

template <unsigned d>
void get_tensor_index(void* pointer,
                      int index,
                      int** out_data,
                      int* out_length)
{
    typedef typename gsBSplineTraits<d, real_t>::Basis Basis;

    Basis* b = static_cast< Basis* >(pointer);

    gsVector<unsigned, d> tensorIndex = b->tensorIndex(index);

    int* vector = new int[d];
    for (int row = 0; row != tensorIndex.rows(); row++)
    {
        vector[row] = static_cast<int>(tensorIndex[row]);
    }

    for (int row = 0; row != static_cast<int>(d); row++)
    {
        std::cout << "  " << vector[row];
    }
    std::cout << std::endl;

    *out_data = vector;
    *out_length = tensorIndex.rows();
}

std::vector<unsigned> get_boxes(int* data,
                                int rows,
                                int cols)
{
    std::vector<unsigned> boxes;
    for (int r = 0; r != rows; r++)
    {
        for (int c = 0; c != cols; c++)
        {
            boxes.push_back(static_cast<unsigned>(data[r * cols + c]));
        }
    }

    std::cout << "Boxes: ";
    for (size_t i = 0; i != boxes.size(); i++)
    {
        std::cout << boxes[i] << " ";
    }
    std::cout << std::endl;

    return boxes;
}

template <unsigned d>
gsTensorBSpline<d>* bspline_fitting(gsTensorBSplineBasis<d>* basis,
                                    double* par,
                                    double* pts,
                                    int rowsPts,
                                    int cols,
                                    double lambda,
                                    int iterations,
                                    double threshold,
                                    double percentageThreshold,
                                    int errorType, //  0 == euclidian distance,
                                    //  1 == max norm
                                    double* maxError,
                                    double* percentileBelow)
{
    gsStopwatch clock;
    gsMatrix<> uv;
    fill_matrix(uv, par, static_cast<int>(d), cols);

    gsMatrix<> xyz;
    fill_matrix(xyz, pts, rowsPts, cols);

    gsFitting<> fitting(uv, xyz, *basis);

    gsTensorBSpline<d>* bspline = NULL;

    for (int i = 0; i != iterations; i++)
    {
        clock.restart();
        if (i != 0)
        {
            basis->uniformRefine();
        }

        fitting.compute(lambda);

        bspline = static_cast< gsTensorBSpline<d>* >(fitting.result());

        if (errorType == 0)
        {
            fitting.computeErrors();
        }
        else if (errorType == 1)
        {
            fitting.computeMaxNormErrors();
        }
        else
        {
            std::cout << "Unknown error type specifier." << std::endl;
            return bspline;
        }

        const std::vector<real_t>& errors = fitting.pointWiseErrors();

        *maxError = fitting.maxPointError();
        *percentileBelow = percent_point_below_threshold(errors, threshold);

        std::cout << "iteration " << i + 1 << " / " << iterations << "\n"
                  << "dofs: "     << bspline->coefs().rows() << "\n"
                  << "max error: " << *maxError << "\n"
                  << "percent of data below threshold: " << percentileBelow << "\n"
                  << "time for this step: " << clock.stop() << std::endl;

        if (*maxError < threshold)
        {
            std::cout << "Max error is below treshold: "
                      << *maxError << " < " << threshold << std::endl;
            break;
        }
        if (*percentileBelow > percentageThreshold)
        {
            std::cout << "Percentage of data below threshold is higher than the given limit: "
                      << *percentileBelow << " > " << percentageThreshold << std::endl;
            break;
        }
    }

    return new gsTensorBSpline<d>(*bspline);
}


template <unsigned d>
gsTHBSpline<d>* thb_fitting(gsTensorBSplineBasis<d>* basis,
                            double* par,
                            double* pts,
                            int rowsPts,
                            int cols,
                            double lambda,
                            int iterations,
                            double threshold,
                            double percentageThreshold,
                            double refThreshold,
                            double refPercent,
                            int extension,
                            int errorType, //  0 == euclidian distance,
                            //  1 == max norm
                            double* maxError,
                            double* percentileBelow)
{
    gsStopwatch clock;
    gsMatrix<> uv;
    fill_matrix(uv, par, static_cast<int>(d), cols);

    gsMatrix<> xyz;
    fill_matrix(xyz, pts, rowsPts, cols);

    gsTHBSplineBasis<d> THB (*basis);
    std::vector<unsigned> ext(d, extension);
    gsHFitting<d, real_t> fitting(uv, xyz, THB, refPercent, ext, lambda);

    gsTHBSpline<d>* thb = NULL;

    for (int i = 0; i != iterations; i++)
    {
        clock.restart();
        fitting.nextIteration(threshold, refThreshold);

        thb = static_cast< gsTHBSpline<d>* >(fitting.result());

        if (errorType == 0)
        {
            fitting.computeErrors();
        }
        else if (errorType == 1)
        {
            fitting.computeMaxNormErrors();
        }
        else
        {
            std::cout << "Unknown error type specifier." << std::endl;
            return thb;
        }

        const std::vector<real_t>& errors = fitting.pointWiseErrors();

        *maxError = fitting.maxPointError();
        *percentileBelow = percent_point_below_threshold(errors, threshold);

        std::cout << "iteration " << i + 1 << " / " << iterations << "\n"
                  << "dofs: "     << thb->coefs().rows() << "\n"
                  << "maxError: " << *maxError << "\n"
                  << "percent of data below treshold: " << percentileBelow << "\n"
                  << "time for this step: " << clock.stop() << std::endl;

        if (*maxError < threshold)
        {
            std::cout << "maxError is below treshold: "
                      << *maxError << " < " << threshold << std::endl;
            break;
        }
        if (*percentileBelow > percentageThreshold)
        {
            std::cout << "Percentage of data below threshold is higher than the given limit: "
                      << *percentileBelow << " > " << percentageThreshold << std::endl;
            break;
        }
    }

    return new gsTHBSpline<d>(*thb);
}

extern "C"
{

    PK_ERROR_code_t deformationEvaluatorC(PK_VECTOR_t position,
                                          PK_FACE_t face,
                                          PK_LOGICAL_t have_params,
                                          PK_UV_t params,
                                          PK_POINTER_t external_data,
                                          PK_VECTOR_t *const deformed_position)
    {
        return deformationEvaluator(position,face,have_params,params,external_data,deformed_position);
    }

    struct gs_Data
    {
        void* ptr;
        int dim;
    };


    gs_Data gs_make_uniform_bspline_basis(double start,
                                          double end,
                                          int numInteriorKnots,
                                          int degree,
                                          int interiorMultiplicity)
    {
        gs_Data basis;
        basis.ptr = new gsBSplineBasis<>(start, end, numInteriorKnots, degree, interiorMultiplicity);
        basis.dim = 1;

        return basis;
    }


    gs_Data gs_make_bspline_basis(double* data,
                                  int numKnots,
                                  int degree)
    {
        std::vector<double> knots;
        for (int i = 0; i != numKnots; i++)
        {
            knots.push_back(data[i]);
        }

        gsKnotVector<> kv(knots, degree);
        gs_Data basis = {.ptr = new gsBSplineBasis<>(kv), .dim = 1};
        return basis;
    }

    gs_Data gs_make_bspline_basis_2D(gs_Data basis1,
                                     gs_Data basis2)
    {
        gsBSplineBasis<>* b1 = static_cast< gsBSplineBasis<>* >(basis1.ptr);
        gsBSplineBasis<>* b2 = static_cast< gsBSplineBasis<>* >(basis2.ptr);

        gs_Data tensorBasis = {.ptr = new gsTensorBSplineBasis<2>(b1, b2),
                               .dim = 2};
        return tensorBasis;
    }

    gs_Data gs_make_bspline_basis_3D(gs_Data basis1,
                                     gs_Data basis2,
                                     gs_Data basis3)
    {
        gsBSplineBasis<>* b1 = static_cast< gsBSplineBasis<>* >(basis1.ptr);
        gsBSplineBasis<>* b2 = static_cast< gsBSplineBasis<>* >(basis2.ptr);
        gsBSplineBasis<>* b3 = static_cast< gsBSplineBasis<>* >(basis3.ptr);

        gs_Data tensorBasis = {.ptr = new gsTensorBSplineBasis<3> (b1, b2, b3),
                               .dim = 3};

        return tensorBasis;
    }


    gs_Data gs_load_bspline(char* filename)
    {
        //std::cout << "Loading file: " << filename << std::endl;

        gs_Data bspl;

        gsFileData<> data(filename);
        if (data.has< gsBSpline<> >())
        {
            bspl.ptr = data.getFirst< gsBSpline<> >();
            bspl.dim = 1;
        }
        else if (data.has< gsTensorBSpline<2> >())
        {
            bspl.ptr = data.getFirst< gsTensorBSpline<2> >();
            bspl.dim = 2;
        }
        else if (data.has< gsTensorBSpline<3> >())
        {
            bspl.ptr = data.getFirst< gsTensorBSpline<3> >();
            bspl.dim = 3;
        }
        else
        {
            bspl.ptr = NULL;
            bspl.dim = -1;
            std::cout << "File doesn't contain B-splines." << std::endl;
        }

        return bspl;
    }

    void gs_save_bspline(gs_Data spline,
                         char* filename)
    {
        gsFileData<> data;
        if (spline.dim == 1)
        {
            gsBSpline<>* geom = static_cast< gsBSpline<>* >(spline.ptr);
            data << *geom;
        }
        else if (spline.dim == 2)
        {
            gsTensorBSpline<2>* geom = static_cast< gsTensorBSpline<2>* >(spline.ptr);
            data << *geom;
        }
        else if (spline.dim == 3)
        {
            gsTensorBSpline<3>* geom = static_cast< gsTensorBSpline<3>* >(spline.ptr);
            data << *geom;
        }
        else
        {
            std::cout << "Dimension not supported." << std::endl;
        }

        data.dump(filename);
    }


    gs_Data gs_make_geometry(gs_Data basis,
                             double* data,
                             int num_coefs,
                             int dim)
    {
        gsMatrix<> coefs(num_coefs, dim);
        for (int row = 0; row != num_coefs; row++)
        {
            for (int col = 0; col != dim; col++)
            {
                coefs(row, col) = data[row * dim + col];
            }
        }
        gsBasis<>* b = static_cast< gsBasis<>* >(basis.ptr);

        gs_Data bspl = {.ptr = (b->makeGeometry(coefs)).release(),
                        .dim = basis.dim};

        return bspl;
    }

    gs_Data gs_get_bspline_basis(gs_Data bspl)
    {
        gs_Data basis;
        basis.dim = bspl.dim;

        if (bspl.dim == 1)
        {
            gsBSpline<>* geom = static_cast< gsBSpline<>* >(bspl.ptr);
            basis.ptr = new gsBSplineBasis<>(geom->basis());
        }
        else if (bspl.dim == 2)
        {
            gsTensorBSpline<2>* geom = static_cast< gsTensorBSpline<2>* >(bspl.ptr);
            basis.ptr = new gsTensorBSplineBasis<2>(geom->basis());
        }
        else if (bspl.dim == 3)
        {
            gsTensorBSpline<3>* geom = static_cast< gsTensorBSpline<3>* >(bspl.ptr);
            basis.ptr = new gsTensorBSplineBasis<3>(geom->basis());
        }
        else
        {
            basis.dim = -1;
            std::cout << "Dimension not supported." << std::endl;
        }

        return basis;
    }

    void gs_get_knot_vector(gs_Data basis,
                            int direction,
                            double** data,
                            int* numElements)
    {
        std::cout << "direction: " << direction << std::endl;

        if (!(0 <= direction && direction < basis.dim))
        {
            std::cout << "Direction is not between 0 and " << basis.dim << std::endl;
            return;
        }

        std::vector<double> knots;

        if (basis.dim == 1)
        {
            gsBSplineBasis<>* ptr = static_cast< gsBSplineBasis<>* >(basis.ptr);
            knots = ptr->knots().get();
        }
        else if (basis.dim == 2)
        {
            gsTensorBSplineBasis<2>* ptr = static_cast< gsTensorBSplineBasis<2>* >(basis.ptr);
            knots = ptr->knots(direction).get();
        }
        else if (basis.dim == 3)
        {
            gsTensorBSplineBasis<3>* ptr = static_cast< gsTensorBSplineBasis<3>* >(basis.ptr);
            knots = ptr->knots(direction).get();
        }
        else
        {
            std::cout << "Dimension not supported..." << std::endl;
            return;
        }

        *data = make_c_array(knots);
        *numElements = static_cast<int>(knots.size());
    }

    void gs_get_coefs(gs_Data spline,
                      const double** data,
                      int* rows,
                      int* cols)
    {
        gsGeometry<>* geom = static_cast< gsGeometry<>* >(spline.ptr);
        const gsMatrix<>& coefs = geom->coefs();
        *data = &coefs.data()[0];
        *rows = coefs.rows();
        *cols = coefs.cols();
    }

    int gs_get_degree(gs_Data basis,
                      int direction)
    {
        gsBasis<>* b = static_cast< gsBasis<>* >(basis.ptr);
        return b->degree(direction);
    }

    void gs_get_support(gs_Data basis,
                        const double** data)
    {
        gsBasis<>* b = static_cast< gsBasis<>* >(basis.ptr);
        const gsMatrix<>& support = b->support();
        *data = make_c_array(support);
    }

    void gs_get_basis_fun_support(gs_Data basis,
                                  int basisFunction,
                                  const double** data)
    {
        gsBasis<>* b = static_cast< gsBasis<>* >(basis.ptr);
        const gsMatrix<>& support = b->support(basisFunction);
        std::cout << "support: " << support << std::endl;
        *data = make_c_array(support);
    }

    int gs_get_size(gs_Data basis)
    {
        gsBasis<>* b = static_cast< gsBasis<>* >(basis.ptr);
        return b->size();
    }

    PK_BSURF_t gs_create_PK_BSURF(gs_Data bspline,
                                  int closedU,
                                  int closedV)
    {
        if (bspline.dim != 2)
        {
            std::cout << "Dimension " << bspline.dim << " is not equal to 2"
                      << std::endl;
        }

        PK_BSURF_t bsurf;
        gsTensorBSpline<2>* b = static_cast< gsTensorBSpline<2>* >(bspline.ptr);

        bool clU = false;
        if (closedU == 1)
        {
            clU = true;
        }

        bool clV = false;
        if (closedV == 1)
        {
            clV = true;
        }

        if (!extensions::createPK_BSURF(*b, bsurf, clU, clV))
        {
            std::cout << "pk::BSURF not created successfully" << std::endl;
        }

        return bsurf;
    }

    PK_BCURVE_t gs_create_PK_BCURVE(gs_Data bspline)
    {
        if (bspline.dim != 1)
        {
            std::cout << "Dimension " << bspline.dim << " is not equal to 1"
                      << std::endl;
        }

        PK_BCURVE_t bcurve;
        gsBSpline<>* b = static_cast< gsBSpline<>* >(bspline.ptr);
        extensions::createPK_BCURVE(*b, bcurve);
        return bcurve;
    }


    void gs_delete_bspline(gs_Data bspline)
    {
        void* ptr = bspline.ptr;
        const int dim = bspline.dim;
        if (dim == 1)
        {
            gsBSpline<>* bspl = static_cast< gsBSpline<>* >(ptr);
            delete bspl;
        }
        else if (dim == 2)
        {
            gsTensorBSpline<2>* bspl = static_cast< gsTensorBSpline<2>* >(ptr);
            delete bspl;
        }
        else if (dim == 3)
        {
            gsTensorBSpline<3>* bspl = static_cast< gsTensorBSpline<3>* >(ptr);
            delete bspl;
        }
        else
        {
            std::cout << "Data doesn't contain B-Splines" << std::endl;
            return;
        }
    }

    void gs_delete_bspline_basis(gs_Data basis)
    {
        void* ptr = basis.ptr;
        const int dim = basis.dim;
        if (dim == 1)
        {
            gsBSplineBasis<>* b = static_cast< gsBSplineBasis<>* >(ptr);
            delete b;
        }
        else if (dim == 2)
        {
            gsTensorBSplineBasis<2>* b = static_cast< gsTensorBSplineBasis<2>* >(ptr);
            delete b;
        }
        else if (dim == 3)
        {
            gsTensorBSplineBasis<3>* b = static_cast< gsTensorBSplineBasis<3>* >(ptr);
            delete b;
        }
        else
        {
            std::cout << "Data doesn't contain B-Spline-Basis" << std::endl;
            return;
        }
    }

    PK_ASSEMBLY_t gs_make_mesh(gs_Data data,
                               int sampling)
    {
        gsGeometry<>* geom = static_cast< gsGeometry<>* >(data.ptr);
        gsMesh<> mesh;

        makeMesh<>(geom->basis(), mesh, sampling);
        geom->evaluateMesh(mesh);

        PK_ASSEMBLY_t assembly;

        extensions::exportMesh(mesh, assembly);

        return assembly;
    }


    void gs_eval_geometry(gs_Data bspline,
                          double* data,
                          int rows,
                          int cols,
                          double** out_data,
                          int* out_rows,
                          int* out_cols,
                          int mode)
    {
        gsGeometry<>* geom = static_cast< gsGeometry<>* >(bspline.ptr);

        gsMatrix<> params;
        fill_matrix(params, data, rows, cols);

        gsMatrix<> result;
        if (mode == 0)
        {
            geom->eval_into(params, result);
        }
        else if (mode == 1)
        {
            geom->deriv_into(params, result);
        }
        else if (mode == 2)
        {
            geom->deriv2_into(params, result);
        }

        //  std::cout << "params: \n" << params << "\n"
        //        << "points: \n" << result << std::endl;

        *out_data = make_c_array(result);
        *out_rows = result.rows();
        *out_cols = result.cols();

    }

    void gs_delete_array(double** array)
    {
        delete[] *array;
    }

    // THB-splines below

    gs_Data gs_make_thb_spline_basis(gs_Data basis)
    {
        gs_Data thb;
        thb.dim = basis.dim;
        thb.ptr = NULL;

        if (basis.dim == 1)
        {
            gsBSplineBasis<>* b = static_cast< gsBSplineBasis<>* >(basis.ptr);
            // const gsKnotVector<>& kv = b->knots();
            // gsCompactKnotVector<> ckv(kv);
            // gsTHBSplineBasis<1>::tensorBasis basis(ckv);
            thb.ptr = new gsTHBSplineBasis<1>(*b);
        }
        else if (basis.dim == 2)
        {
            gsTensorBSplineBasis<2>* b = static_cast< gsTensorBSplineBasis<2>* >(basis.ptr);
            thb.ptr = new gsTHBSplineBasis<2>(*b);
        }
        else if (basis.dim == 3)
        {
            gsTensorBSplineBasis<3>* b = static_cast< gsTensorBSplineBasis<3>* >(basis.ptr);
            thb.ptr = new gsTHBSplineBasis<3>(*b);
        }
        else
        {
            std::cout << "Dimension not supported." << std::endl;
        }

        return thb;
    }


    gs_Data gs_load_thb_spline(char* filename)
    {
        //std::cout << "Loading file: " << filename << std::endl;

        gs_Data thb;

        gsFileData<> data(filename);
        if (data.has< gsTHBSpline<1> >())
        {
            thb.ptr = data.getFirst< gsTHBSpline<1> >();
            thb.dim = 1;
        }
        else if (data.has< gsTHBSpline<2> >())
        {
            thb.ptr = data.getFirst< gsTHBSpline<2> >();
            thb.dim = 2;
        }
        else if (data.has< gsTHBSpline<3> >())
        {
            thb.ptr = data.getFirst< gsTHBSpline<3> >();
            thb.dim = 3;
        }
        else
        {
            thb.ptr = NULL;
            thb.dim = -1;
            std::cout << "File doesn't contain THB-splines." << std::endl;
        }

        return thb;
    }

    void gs_save_thb_spline(gs_Data thb,
                            char* filename)
    {
        gsFileData<> data;
        if (thb.dim == 1)
        {
            gsTHBSpline<1>* geom = static_cast< gsTHBSpline<1>* >(thb.ptr);
            data << *geom;
        }
        else if (thb.dim == 2)
        {
            gsTHBSpline<2>* geom = static_cast< gsTHBSpline<2>* >(thb.ptr);
            data << *geom;
        }
        else if (thb.dim == 3)
        {
            gsTHBSpline<3>* geom = static_cast< gsTHBSpline<3>* >(thb.ptr);
            data << *geom;
        }
        else
        {
            std::cout << "Dimension not supported." << std::endl;
        }

        data.dump(filename);
    }


    gs_Data gs_get_thb_spline_basis(gs_Data thb)
    {
        gs_Data basis;
        basis.dim = thb.dim;

        if (thb.dim == 1)
        {
            gsTHBSpline<1>* geom = static_cast< gsTHBSpline<1>* >(thb.ptr);
            basis.ptr = new gsTHBSplineBasis<1>(geom->basis());
        }
        else if (thb.dim == 2)
        {
            gsTHBSpline<2>* geom = static_cast< gsTHBSpline<2>* >(thb.ptr);
            basis.ptr = new gsTHBSplineBasis<2>(geom->basis());
        }
        else if (thb.dim == 3)
        {
            gsTHBSpline<3>* geom = static_cast< gsTHBSpline<3>* >(thb.ptr);
            basis.ptr = new gsTHBSplineBasis<3>(geom->basis());
        }
        else
        {
            basis.dim = -1;
            std::cout << "Dimension not supported." << std::endl;
        }
        return basis;
    }

    int gs_get_max_hier_level(gs_Data basis)
    {
        if (basis.dim == 1)
        {
            gsHTensorBasis<1>* b = static_cast< gsHTensorBasis<1>* >(basis.ptr);
            return static_cast<int>(b->maxLevel());
        }
        else if (basis.dim == 2)
        {
            gsHTensorBasis<2>* b = static_cast< gsHTensorBasis<2>* >(basis.ptr);
            return static_cast<int>(b->maxLevel());
        }
        else if (basis.dim == 3)
        {
            gsHTensorBasis<3>* b = static_cast< gsHTensorBasis<3>* >(basis.ptr);
            return static_cast<int>(b->maxLevel());
        }
        else
        {
            std::cout << "Dimension not supported." << std::endl;
            return -1;
        }
    }

    int gs_hier_level_of_basis_function(gs_Data basis,
                                        int basis_fun)
    {
        std::cout << "dim: " << basis.dim << "\n"
                  << "basis_fun: " << basis_fun << std::endl;

        if (basis.dim == 1)
        {
            gsHTensorBasis<1>* b = static_cast< gsHTensorBasis<1>* >(basis.ptr);
            return b->levelOf(static_cast< unsigned >(basis_fun));
        }
        else if (basis.dim == 2)
        {
            gsHTensorBasis<2>* b = static_cast< gsHTensorBasis<2>* >(basis.ptr);
            return b->levelOf(static_cast< unsigned >(basis_fun));
        }
        else if (basis.dim == 3)
        {
            gsHTensorBasis<3>* b = static_cast< gsHTensorBasis<3>* >(basis.ptr);
            return b->levelOf(static_cast< unsigned >(basis_fun));
        }
        else
        {
            std::cout << "Dimension not supported." << std::endl;
            return -1;
        }
    }

    void gs_delete_thb_spline_basis(gs_Data basis)
    {
        void* ptr = basis.ptr;
        const int dim = basis.dim;
        if (dim == 1)
        {
            gsTHBSplineBasis<1>* b = static_cast< gsTHBSplineBasis<1>* >(ptr);
            delete b;
        }
        else if (dim == 2)
        {
            gsTHBSplineBasis<2>* b = static_cast< gsTHBSplineBasis<2>* >(ptr);
            delete b;
        }
        else if (dim == 3)
        {
            gsTHBSplineBasis<3>* b = static_cast< gsTHBSplineBasis<3>* >(ptr);
            delete b;
        }
        else
        {
            std::cout << "Data doesn't contain THB-Spline-Basis" << std::endl;
        }
    }

    void gs_delete_thb_spline(gs_Data thb)
    {
        void* ptr = thb.ptr;
        const int dim = thb.dim;
        if (dim == 1)
        {
            gsTHBSpline<1>* b = static_cast< gsTHBSpline<1>* >(ptr);
            delete b;
        }
        else if (dim == 2)
        {
            gsTHBSpline<2>* b = static_cast< gsTHBSpline<2>* >(ptr);
            delete b;
        }
        else if (dim == 3)
        {
            gsTHBSpline<3>* b = static_cast< gsTHBSpline<3>* >(ptr);
            delete b;
        }
        else
        {
            std::cout << "Data doesn't contain THB-Splines" << std::endl;
            return;
        }
    }

    gs_Data gs_get_bspline_basis_from_hier_basis(gs_Data thb,
                                                 int level)
    {
        const int dim = thb.dim;
        gs_Data b_basis;
        b_basis.dim = dim;
        b_basis.ptr = NULL;

        if (dim == 1)
        {
            gsHTensorBasis<1>* b = static_cast< gsHTensorBasis<1>* >(thb.ptr);
            //const gsTensorBSplineBasis<1, real_t, gsCompactKnotVector<> >* tbsb = b->getBases()[level];
            const gsHTensorBasis<1>::tensorBasis* tbsb = b->getBases()[level];
            //gsBSplineBasis<>* basis = new gsBSplineBasis<>(tbsb->knots(0));
            // b_basis.ptr = new gsTensorBSplineBasis<1>(basis);
            b_basis.ptr = new gsBSplineBasis<>(tbsb->knots(0));
        }
        else if (dim == 2)
        {
            gsHTensorBasis<2>* b = static_cast< gsHTensorBasis<2>* >(thb.ptr);
            const gsTensorBSplineBasis<2, real_t>* tbsb = b->getBases()[level];
            b_basis.ptr = new gsTensorBSplineBasis<2>(tbsb->knots(0), tbsb->knots(1));
        }
        else if (dim == 3)
        {
            gsHTensorBasis<3>* b = static_cast< gsHTensorBasis<3>* >(thb.ptr);
            const gsTensorBSplineBasis<3, real_t>* tbsb = b->getBases()[level];
            b_basis.ptr = new gsTensorBSplineBasis<3>(tbsb->knots(0), tbsb->knots(1), tbsb->knots(2));
        }
        else
        {
            std::cout << "Data doesn't contain THB-Spline-Basis" << std::endl;
        }

        return b_basis;
    }

    PK_ASSEMBLY_t gs_export_thb_to_bsplines(gs_Data thb,
                                            int* error,
                                            double* boxes,
                                            int boxes_length)
    {
        if (thb.dim != 2)
        {
            std::cout << "gs_export_THB_to_BSplines works only for surfaces." << std::endl;
        }

        gsTHBSpline<2>* surface = static_cast< gsTHBSpline<2>* >(thb.ptr);

        std::vector<double> boxes_vec;
        fill_stdvector(boxes_vec, boxes, boxes_length);
        PK_ASSEMBLY_t assembly;
        if (!extensions::exportTHBsurface(*surface,boxes_vec,assembly))
        {
            *error = -1;
            return 0;
        }

        *error = 1;
        return assembly;
    }

    void gs_load_points(char* filename,
                        double** params,
                        int* rowsParams,
                        double** points,
                        int* rowsPoints,
                        int* cols)
    {
        gsFileData<> data(filename);

        if (data.has< gsMatrix<> >())
        {
            gsMatrix<> par, pts;
            data.getId< gsMatrix<> >(0, par);
            data.getId< gsMatrix<> >(1, pts);

            *params = make_c_array(par);
            *points = make_c_array(pts);
            *rowsParams = par.rows();
            *rowsPoints = pts.rows();
            *cols = par.cols();
        }
        else
        {
            std::cout << "The file doesn't contain points." << std::endl;
            *params = NULL;
            *points = NULL;
            *rowsParams = -1;
            *rowsPoints = -1;
            *cols = -1;
        }
    }

    void gs_save_points(double* c_parameters,
                        int par_rows,
                        int par_cols,
                        double* c_points,
                        int pts_rows,
                        int pts_cols,
                        char* filename)
    {
        gsMatrix<> parameters;
        fill_matrix(parameters, c_parameters, par_rows, par_cols);

        gsMatrix<> points;
        fill_matrix(points, c_points, pts_rows, pts_cols);

        gsFileData<> fd;
        fd << parameters;
        fd << points;

        fd.dump(filename);

        std::cout << "Saved points into: " << filename << std::endl;
    }

    gs_Data gs_bspline_fitting(gs_Data b,
                               double* par,
                               int rowsPar,
                               double* pts,
                               int rowsPts,
                               int cols,
                               double lambda,
                               int iterations,
                               double threshold,
                               double percentageThreshold,
                               int errorType,
                               double* maxError,
                               double* percentileBelow)
    {
        gs_Data result;

        if (rowsPar == 2)
        {
            gsTensorBSplineBasis<2>* basis = static_cast< gsTensorBSplineBasis<2>* >(b.ptr);
            result.dim = 2;
            result.ptr = bspline_fitting<2>(basis, par, pts, rowsPts, cols, lambda,
                                            iterations, threshold,percentageThreshold,
                                            errorType, maxError,percentileBelow);
        }
        else if (rowsPar == 3)
        {
            gsTensorBSplineBasis<3>* basis = static_cast< gsTensorBSplineBasis<3>* >(b.ptr);
            result.dim = 3;
            result.ptr = bspline_fitting<3>(basis, par, pts, rowsPts, cols, lambda,
                                            iterations, threshold,percentageThreshold,
                                            errorType, maxError,percentileBelow);

        }
        else
        {
            std::cout << "Not supported dimension." << std::endl;
        }

        return result;
    }

    gs_Data gs_bspline_periodic_fitting(gs_Data b,
                                        double* par,
                                        int rowsPar,
                                        double* pts,
                                        int rowsPts,
                                        int cols,
                                        double lambda,
                                        int iterations,
                                        double threshold,
                                        double percentageThreshold,
                                        double* maxError,
                                        double* percentileBelow)
    {
        gsStopwatch clock;
        gs_Data result;
        result.dim = 2;
        result.ptr = NULL;

        if (rowsPar != 2 || b.dim != 2)
        {
            std::cout << "Periodic fitting is done only for surfaces\n"
                "Aborting..."
                      << std::endl;
            return result;
        }

        gsTensorBSplineBasis<2>* tensorBasis = static_cast< gsTensorBSplineBasis<2>* >(b.ptr);

        gsMatrix<> uv;
        fill_matrix(uv, par, rowsPar, cols);

        gsMatrix<> xyz;
        fill_matrix(xyz, pts, rowsPts, cols);

        std::vector< gsTensorBSplineBasis<2>* > vectorOfBasis;
        vectorOfBasis.push_back(tensorBasis);

        gsCompositeBSplineBasis<2, real_t> periodicBasis(vectorOfBasis, get_periodic_topology());
        gsTensorBSpline<2>* bspline = NULL;

        for (int i = 0; i != iterations; i++)
        {
            clock.restart();
            if (bspline == NULL)
            {
                delete bspline;
            }

            gsTensorBSplineBasis<2> basis = dynamic_cast< const gsTensorBSplineBasis<2>& >
                (periodicBasis.getBase(0));

            gsSparseMatrix<> global2local = periodicBasis.getMapper().asMatrix();

            gsFitting<> fitting(uv, xyz, basis);

            const int numBasis = basis.size();
            gsSparseMatrix<> A_local(numBasis, numBasis);
            gsMatrix<> B_local(numBasis, 3);
            A_local.setZero();
            B_local.setZero();

            fitting.assembleSystem(A_local, B_local);
            if (0 < lambda)
            {
                fitting.applySmoothing(lambda, A_local);
            }

            gsSparseMatrix<> A = global2local.transpose() * A_local * global2local;
            gsMatrix<> B = global2local.transpose() * B_local;

            A.makeCompressed();
            gsSparseSolver<real_t>::BiCGSTABILUT solver(A);
            if ( solver.preconditioner().info() != Eigen::Success )
            {
                gsWarn<<  "The preconditioner failed. Aborting.\n";
            }

            gsMatrix<> coeffs = global2local * solver.solve(B);

            bspline = new gsTensorBSpline<2>(basis, coeffs);

            const std::vector<double> errors = compute_errors(uv, xyz, *bspline);
            *maxError = *std::max_element(errors.begin(), errors.end());
            *percentileBelow = percent_point_below_threshold(errors, threshold);

            std::cout << "iteration " << i + 1 << " / " << iterations << "\n"
                      << "dofs: "     << global2local.cols() << "\n"
                      << "maxError: " << *maxError << "\n"
                      << "percent of data below treshold: " << percentileBelow << "\n"
                      << "time for this step: " << clock.stop() << std::endl;

            if (*maxError < threshold)
            {
                std::cout << "maxError is below treshold: "
                          << *maxError << " < " << threshold << std::endl;
                break;
            }

            if (*percentileBelow > percentageThreshold)
            {
                std::cout << "Percentage of data below threshold is higher than the given limit: "
                          << *percentileBelow << " > " << percentageThreshold << std::endl;
                break;
            }

            if (i != iterations - 1)
            {
                periodicBasis.uniformRefine();
            }
        }

        result.ptr = bspline;
        return result;
    }

    gs_Data gs_thb_spline_fitting(gs_Data b,
                                  double* par,
                                  int rowsPar,
                                  double* pts,
                                  int rowsPts,
                                  int cols,
                                  double lambda,
                                  int iterations,
                                  double threshold,
                                  double percentageThreshold,
                                  double refThreshold,
                                  double refPercent,
                                  int extension,
                                  int errorType,
                                  double* maxError,
                                  double* percentileBelow)
    {
        gs_Data result;

        if (rowsPar == 2)
        {
            gsTensorBSplineBasis<2>* basis = static_cast< gsTensorBSplineBasis<2>* >(b.ptr);
            result.dim = 2;
            result.ptr = thb_fitting<2>(basis, par, pts, rowsPts, cols, lambda,
                                        iterations, threshold, percentageThreshold, refThreshold,
                                        refPercent, extension, errorType, maxError,percentileBelow);
        }
        else if (rowsPar == 3)
        {
            gsTensorBSplineBasis<3>* basis = static_cast< gsTensorBSplineBasis<3>* >(b.ptr);
            result.dim = 3;
            result.ptr = thb_fitting<3>(basis, par, pts, rowsPts, cols, lambda,
                                        iterations, threshold, percentageThreshold, refThreshold,
                                        refPercent, extension, errorType, maxError,percentileBelow);
        }
        else
        {
            std::cout << "Not supported dimension." << std::endl;
        }

        return result;
    }

    gs_Data gs_thb_periodic_fitting(gs_Data b,
                                    double* par,
                                    int rowsPar,
                                    double* pts,
                                    int rowsPts,
                                    int cols,
                                    double lambda,
                                    int iterations,
                                    double threshold,
                                    double percentageThreshold,
                                    int extension,
                                    double* maxError,
                                    int* boxes,
                                    int boxes_length,
                                    double* percentileBelow)
    {
        std::cout << "start" << std::endl;
        gsStopwatch clock;

        gs_Data result;
        result.dim = 2;
        result.ptr = NULL;

        if (rowsPar != 2 || b.dim != 2)
        {
            std::cout << "Periodic fitting is done only for surfaces.\n"
                "Aborting..."
                      << std::endl;

            return result;
        }

        gsTensorBSplineBasis<2>* tensorBasis = static_cast< gsTensorBSplineBasis<2>* >(b.ptr);

        gsMatrix<> uv;
        fill_matrix(uv, par, rowsPar, cols);

        gsMatrix<> xyz;
        fill_matrix(xyz, pts, rowsPts, cols);

        gsTHBSplineBasis<2> thbBasis(*tensorBasis);
        std::vector<unsigned> ext(2, extension);

        // gsBoxTopology topology(2, 1);
        // gsVector<index_t> dirMap(2);
        // dirMap << 0, 1;
        // gsVector<bool> dirOrientation(2);
        // dirOrientation << true, true;

        // patchSide ps1(0, boundary::left);
        // patchSide ps2(0, boundary::right);

        // boundaryInterface boundInter(ps1, ps2, dirMap, dirOrientation);
        // topology.addInterface(boundInter);

        std::vector< gsHTensorBasis<2>* > vectorOfBasis;
        vectorOfBasis.push_back(&thbBasis);

        gsCompositeHBasis<2, real_t> periodicBasis(vectorOfBasis, get_periodic_topology());
        gsTHBSpline<2>* thb = NULL;

        std::vector<unsigned> boxes_vec;
        fill_stdvector(boxes_vec, boxes, boxes_length);

        for (int i = 0; i != iterations; i++)
        {
            clock.restart();
            if (thb == NULL)
            {
                delete thb;
            }

            gsTHBSplineBasis<2> basis = dynamic_cast< const gsTHBSplineBasis<2>& >
                (periodicBasis.getBase(0));

            gsSparseMatrix<> global2local = periodicBasis.getMapper().asMatrix();

            const int numBasis = basis.size();
            gsSparseMatrix<> A_local(numBasis, numBasis);
            gsMatrix<> B_local(numBasis, 3);
            A_local.setZero();
            B_local.setZero();

            gsHFittingLvlConstrained<2, real_t>fitting(uv,xyz,basis,1,ext,lambda,boxes_vec);

            fitting.assembleSystem(A_local, B_local);
            if (0 < lambda)
            {
                fitting.applySmoothing(lambda, A_local);
            }

            gsSparseMatrix<> A = global2local.transpose() * A_local * global2local;
            gsMatrix<> B = global2local.transpose() * B_local;

            A.makeCompressed();
            gsSparseSolver<>::BiCGSTABILUT solver(A);
            if ( solver.preconditioner().info() != Eigen::Success )
            {
                gsWarn<<  "The preconditioner failed. Aborting.\n";
            }

            gsMatrix<> coefs = global2local * solver.solve(B);

            thb = new gsTHBSpline<2>(basis, coefs);

            const std::vector<double> errors = compute_errors(uv, xyz, *thb);
            *maxError = *std::max_element(errors.begin(), errors.end());
            *percentileBelow = percent_point_below_threshold(errors, threshold);

            std::cout << "iteration " << i + 1 << " / " << iterations << "\n"
                      << "dofs: "     << global2local.cols() << "\n"
                      << "maxError: " << *maxError << "\n"
                      << "percent of data below treshold: " << *percentileBelow << "\n"
                      << "time for this step: " << clock.stop() << std::endl;

            if (*maxError < threshold)
            {
                std::cout << "maxError is below treshold: "
                          << *maxError << " < " << threshold << std::endl;
                break;
            }

            if (*percentileBelow > percentageThreshold)
            {
                std::cout << "Percentage of data below threshold is higher than the given limit: "
                          << *percentileBelow << " > " << percentageThreshold << std::endl;
                break;
            }

            if (i != iterations - 1)
            {
                std::vector<unsigned> boxes = fitting.getBoxes(errors, threshold);
                if ( boxes.size()>0 )
                {
                    periodicBasis.refineElements(0, boxes);
                    std::cout << "Inserted " << boxes.size() / 5 << " boxes." << std::endl;
                }
                else
                {
                    std::cout << "Nothing to insert - finishing fitting early..." << std::endl;
                    break;
                }
            }
        }

        result.ptr = thb;
        return result;
    }




    void gs_refine_thb_basis_function(gs_Data thb,
                                      int basisFun)
    {
        if (thb.dim == 2)
        {
            refine_basis_function<gsTHBSpline<2>, 2>(thb.ptr, basisFun);
        }
        else if (thb.dim == 3)
        {
            refine_basis_function<gsTHBSpline<3>, 3>(thb.ptr, basisFun);
        }
    }


    void gs_refine_thb_spline(gs_Data spline,
                              int* data,
                              int rows,
                              int cols)
    {
        std::vector<unsigned> boxes = get_boxes(data, rows, cols);

        if (spline.dim == 2)
        {
            gsTHBSpline<2>* thb = static_cast< gsTHBSpline<2>* >(spline.ptr);
            thb->basis().refineElements_withCoefs(thb->coefs(), boxes);
        }
        else if (spline.dim == 3)
        {
            gsTHBSpline<3>* thb = static_cast< gsTHBSpline<3>* >(spline.ptr);
            thb->basis().refineElements_withCoefs(thb->coefs(), boxes);
        }
    }

    gs_Data gs_convert_thb_to_bspline(gs_Data spline)
    {
        gs_Data result;
        if (spline.dim != 2)
        {
            std::cout << "gs_convert_thb_to_bspline not implemented for d != 2"
                      << std::endl;         // TODO
        }

        gsTHBSpline<2>* thb = static_cast< gsTHBSpline<2>* >(spline.ptr);
        gsTHBSpline<2> thb_copy(*thb);

        gsTensorBSpline<2>* bspline = new gsTensorBSpline<2>();
        thb_copy.convertToBSpline(*bspline);

        result.dim = 2;
        result.ptr = bspline;

        return result;
    }


    void* gs_extract_boundaries(gs_Data spline)
    {
        gsGeometry<>* geom = static_cast< gsGeometry<>* >(spline.ptr);

        gsMultiPatch<> mp(*geom);
        mp.computeTopology();

        gsMultiPatch<> result;

        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
        {
            const gsGeometry<>& g = mp.patch(it->patch);
            result.addPatch(g.boundary(it->side()));
        }

        return new gsMultiPatch<>(result);
    }


    int gs_number_of_patches(void* mp)
    {
        gsMultiPatch<>* multiPatch = static_cast< gsMultiPatch<>* >(mp);
        return static_cast<int>(multiPatch->nPatches());
    }


    gs_Data gs_get_thb_patch(void* mp, int index)
    {
        gsMultiPatch<>* multiPatch = static_cast< gsMultiPatch<>* >(mp);
        gsGeometry<>& geom = multiPatch->patch(index);

        // TO DO: make it dimension free
        gs_Data data;

        gsTHBSpline<2>& thb = dynamic_cast< gsTHBSpline<2>& >(geom);

        data.ptr = new gsTHBSpline<2>(thb);
        data.dim = 2;

        return data;
    }

    gs_Data gs_get_bspline_patch(void* mp, int index)
    {
        gsMultiPatch<>* multiPatch = static_cast< gsMultiPatch<>* >(mp);
        gsGeometry<>& geom = multiPatch->patch(index);

        // TO DO: make it dimension free
        gs_Data data;
        gsTensorBSpline<2>& bspl = dynamic_cast< gsTensorBSpline<2>& >(geom);

        data.ptr = new gsTensorBSpline<2>(bspl);
        data.dim = 2;
        return data;
    }


    void gs_delete_multipatch(void* mp)
    {
        gsMultiPatch<>* multiPatch = static_cast< gsMultiPatch<>* >(mp);
        delete multiPatch;
        mp = NULL;
    }

    int gs_index_of_bspline(gs_Data basis,
                            int* data,
                            int length)
    {
        gsVector<int> vec;
        fill_vector(vec, data, length);

        if (basis.dim == 2)
        {
            return get_index<2>(basis.ptr, vec);
        }
        else if (basis.dim == 3)
        {
            return get_index<3>(basis.ptr, vec);
        }
        else
        {
            std::cout << "Dimensions different than 2 and 3 are not supported" << std::endl;
            return -1;
        }
    }

    void gs_tensor_index_of_bspline(gs_Data basis,
                                    int index,
                                    int** out_data,
                                    int* out_length)
    {
        if (basis.dim == 2)
        {
            get_tensor_index<2>(basis.ptr, index, out_data, out_length);
        }
        else if (basis.dim == 3)
        {
            get_tensor_index<3>(basis.ptr, index, out_data, out_length);
        }
        else
        {
            std::cout << "Dimensions different than 2 and 3 are not supported" << std::endl;
        }
    }

    void gs_degree_elevate(gs_Data basis,
                           int amount,
                           int dir)
    {
        gsBasis<>* b = static_cast< gsBasis<>* >(basis.ptr);
        b->degreeElevate(amount, dir);
    }

    // --------------------------------------------------------------------------------
    // HB-splines

    gs_Data gs_make_hb_spline_basis(gs_Data basis)
    {
        gs_Data thb;
        thb.dim = basis.dim;
        thb.ptr = NULL;

        if (basis.dim == 1)
        {
            gsBSplineBasis<>* b = static_cast< gsBSplineBasis<>* >(basis.ptr);
            // const gsKnotVector<>& kv = b->knots();
            // gsCompactKnotVector<> ckv(kv);
            // gsTHBSplineBasis<1>::tensorBasis basis(ckv);
            thb.ptr = new gsHBSplineBasis<1>(*b);
        }
        else if (basis.dim == 2)
        {
            gsTensorBSplineBasis<2>* b = static_cast< gsTensorBSplineBasis<2>* >(basis.ptr);
            thb.ptr = new gsHBSplineBasis<2>(*b);
        }
        else if (basis.dim == 3)
        {
            gsTensorBSplineBasis<3>* b = static_cast< gsTensorBSplineBasis<3>* >(basis.ptr);
            thb.ptr = new gsHBSplineBasis<3>(*b);
        }
        else
        {
            std::cout << "Dimension not supported." << std::endl;
        }

        return thb;
    }

    gs_Data gs_load_hb_spline(char* filename)
    {
        //std::cout << "Loading file: " << filename << std::endl;

        gs_Data hb;

        gsFileData<> data(filename);
        if (data.has< gsHBSpline<1> >())
        {
            hb.ptr = data.getFirst< gsHBSpline<1> >();
            hb.dim = 1;
        }
        else if (data.has< gsHBSpline<2> >())
        {
            hb.ptr = data.getFirst< gsHBSpline<2> >();
            hb.dim = 2;
        }
        else if (data.has< gsHBSpline<3> >())
        {
            hb.ptr = data.getFirst< gsHBSpline<3> >();
            hb.dim = 3;
        }
        else
        {
            hb.ptr = NULL;
            hb.dim = -1;
            std::cout << "File doesn't contain HB-splines." << std::endl;
        }

        return hb;
    }

    void gs_save_hb_spline(gs_Data hb,
                           char* filename)
    {
        gsFileData<> data;
        if (hb.dim == 1)
        {
            gsHBSpline<1>* geom = static_cast< gsHBSpline<1>* >(hb.ptr);
            data << *geom;
        }
        else if (hb.dim == 2)
        {
            gsHBSpline<2>* geom = static_cast< gsHBSpline<2>* >(hb.ptr);
            data << *geom;
        }
        else if (hb.dim == 3)
        {
            gsHBSpline<3>* geom = static_cast< gsHBSpline<3>* >(hb.ptr);
            data << *geom;
        }
        else
        {
            std::cout << "Dimension not supported." << std::endl;
        }

        data.dump(filename);
    }

    gs_Data gs_get_hb_spline_basis(gs_Data hb)
    {
        gs_Data basis;
        basis.dim = hb.dim;

        if (hb.dim == 1)
        {
            gsHBSpline<1>* geom = static_cast< gsHBSpline<1>* >(hb.ptr);
            basis.ptr = new gsHBSplineBasis<1>(geom->basis());
        }
        else if (hb.dim == 2)
        {
            gsHBSpline<2>* geom = static_cast< gsHBSpline<2>* >(hb.ptr);
            basis.ptr = new gsHBSplineBasis<2>(geom->basis());
        }
        else if (hb.dim == 3)
        {
            gsHBSpline<3>* geom = static_cast< gsHBSpline<3>* >(hb.ptr);
            basis.ptr = new gsHBSplineBasis<3>(geom->basis());
        }
        else
        {
            basis.dim = -1;
            std::cout << "Dimension not supported." << std::endl;
        }
        return basis;
    }

    void gs_delete_hb_spline_basis(gs_Data basis)
    {
        void* ptr = basis.ptr;
        const int dim = basis.dim;
        if (dim == 1)
        {
            gsHBSplineBasis<1>* b = static_cast< gsHBSplineBasis<1>* >(ptr);
            delete b;
        }
        else if (dim == 2)
        {
            gsHBSplineBasis<2>* b = static_cast< gsHBSplineBasis<2>* >(ptr);
            delete b;
        }
        else if (dim == 3)
        {
            gsHBSplineBasis<3>* b = static_cast< gsHBSplineBasis<3>* >(ptr);
            delete b;
        }
        else
        {
            std::cout << "Data doesn't contain HB-Spline-Basis" << std::endl;
        }
    }

    void gs_delete_hb_spline(gs_Data thb)
    {
        void* ptr = thb.ptr;
        const int dim = thb.dim;
        if (dim == 1)
        {
            gsHBSpline<1>* b = static_cast< gsHBSpline<1>* >(ptr);
            delete b;
        }
        else if (dim == 2)
        {
            gsHBSpline<2>* b = static_cast< gsHBSpline<2>* >(ptr);
            delete b;
        }
        else if (dim == 3)
        {
            gsHBSpline<3>* b = static_cast< gsHBSpline<3>* >(ptr);
            delete b;
        }
        else
        {
            std::cout << "Data doesn't contain HB-Splines" << std::endl;
            return;
        }
    }

    void gs_refine_hb_basis_function(gs_Data hb,
                                     int basisFun)
    {
        if (hb.dim == 2)
        {
            refine_basis_function<gsHBSpline<2>, 2>(hb.ptr, basisFun);
        }
        else if (hb.dim == 3)
        {
            refine_basis_function<gsHBSpline<3>, 3>(hb.ptr, basisFun);
        }
    }

    void gs_refine_hb_spline(gs_Data spline,
                             int* data,
                             int rows,
                             int cols)
    {

        std::vector<unsigned> boxes = get_boxes(data, rows, cols);

        if (spline.dim == 2)
        {
            gsHBSpline<2>* thb = static_cast< gsHBSpline<2>* >(spline.ptr);
            thb->basis().refineElements_withCoefs(thb->coefs(), boxes);
        }
        else if (spline.dim == 3)
        {
            gsHBSpline<3>* thb = static_cast< gsHBSpline<3>* >(spline.ptr);
            thb->basis().refineElements_withCoefs(thb->coefs(), boxes);
        }
    }


    // --------------------------------------------------------------------------------
    // Functions for expressing THB spline as an pk_FSURF in PARASOLID.

    void gs_evaluate_surface_as_foreign(int low,
                                        int high,
                                        double u,
                                        double v,
                                        int nu,
                                        int nv,
                                        int triang,
                                        double* coord)
    {
        unsigned long long int_address = (unsigned long long) high << 32 | low;
        void* address = reinterpret_cast< void* >(int_address);
        gsGeometry<>* geometry = static_cast< gsGeometry<>* >(address);

        gsMatrix<real_t, 2, 1> uv;
        uv(0, 0) = u;
        uv(1, 0) = v;

        gsMatrix<> eval;
        gsMatrix<> der;
        gsMatrix<> der2;

        int n = 0;
        if (ev_P(n, nu, nv, triang))
        {
            geometry->eval_into(uv, eval);
        }

        if (ev_dP_du(n, nu, nv, triang) ||
            ev_dP_dv(n, nu, nv, triang))
        {
            geometry->deriv_into(uv, der);
        }

        if (ev_dP2_du2(n, nu, nv, triang) ||
            ev_dP2_dv2(n, nu, nv, triang) ||
            ev_dP2_dudv(n, nu, nv, triang))
        {
            geometry->deriv2_into(uv, der2);
        }

        n = 0;
        if (ev_P(n, nu, nv, triang)) // P
        {
            coord[n++] = eval(0, 0);
            coord[n++] = eval(1, 0);
            coord[n] = eval(2, 0);
        }

        if (ev_dP_du(n, nu, nv, triang)) // Pu
        {
            coord[n++] = der(0, 0);
            coord[n++] = der(2, 0);
            coord[n] = der(4, 0);
        }

        if (ev_dP2_du2(n, nu, nv, triang)) // Puu
        {
            coord[n++] = der2(0, 0);
            coord[n++] = der2(3, 0);
            coord[n] = der2(6, 0);
        }

        if (ev_dP_dv(n, nu, nv, triang)) // Pv
        {
            coord[n++] = der(1, 0);
            coord[n++] = der(3, 0);
            coord[n] = der(5, 0);
        }

        if (ev_dP2_dv2(n, nu, nv, triang)) // Pvv
        {
            coord[n++] = der2(1, 0);
            coord[n++] = der2(4, 0);
            coord[n] = der2(7, 0);
        }

        if (ev_dP2_dudv(n, nu, nv, triang)) // Puv
        {
            coord[n++] = der2(2, 0);
            coord[n++] = der2(5, 0);
            coord[n++] = der2(8, 0);
        }

    }

    PK_FSURF_t gs_create_foreign_geometry(gs_Data spline)
    {
        if (spline.dim != 2)
        {
            std::cout << "Error, given geometry must be surface." << std::endl;
            return 0;
        }

        if ((sizeof(spline.ptr) != 8) ||
            (sizeof(unsigned long long) != 8) ||
            (sizeof(unsigned int) != 4))
        {
            std::cout << "This is not the right platform." << std::endl;
            return 0;
        }

        PK_FSURF_sf_t sform;

        // --------------------------------------------------------------------------------
        static const char* foreignGeometryName = "pku:gismo";
        sform.key = const_cast< char* >(foreignGeometryName);
        sform.n_ints = 2;

        // --------------------------------------------------------------------------------
        // horrible but it works on MTU machines:
        // split an 64 bit pointer to two 32 bit ints
        unsigned long long data64 = reinterpret_cast<unsigned long long>(spline.ptr);
        unsigned int high = (unsigned int)((data64 & 0xFFFFFFFF00000000) >> 32);
        unsigned int low = (unsigned int)(data64 & 0xFFFFFFFF);

        int* ints = new int[2]; // yes here is a memory leak
        ints[0] = low;
        ints[1] = high;
        sform.ints = ints;
        // --------------------------------------------------------------------------------

        sform.n_doubles = 4;

        // --------------------------------------------------------------------------------
        gsGeometry<>* geom = static_cast< gsGeometry<>* >(spline.ptr);
        const gsMatrix<>& support = geom->support();

        double* reals = new double[4]; // yes
        reals[0] = support(0, 0);
        reals[1] = support(0, 1);
        reals[2] = support(1, 0);
        reals[3] = support(1, 1);
        sform.doubles = reals;
        // --------------------------------------------------------------------------------

        sform.space = 1;
        sform.transf = PK_ENTITY_null;


        PK_FSURF_t fsurf;
        PK_ERROR_code_t err = PK_FSURF_create(&sform, &fsurf);
        if (err)
        {
            std::cout << "Parasolid PK_FSURF_create: " << err << std::endl;
        }

        return fsurf;
    }


    void gs_export_to_vtk(gs_Data spline,
                          char* file,
                          int geom,
                          int knot_configuration,
                          int physical_knot_configuration,
                          int basis,
                          int control_points,
                          int control_polygon,
                          int geom_resolution,
                          int basis_resolution,
                          int cell_resolution)
    {
        gsGeometry<>* geometry = static_cast< gsGeometry<>* >(spline.ptr);

        std::string filename = file;

        if (geom)
        {
            std::string out = filename + "Geometry";
            std::cout << "Exporting: " << out << std::endl;
            gsWriteParaview(*geometry, out, geom_resolution);
        }

        if (knot_configuration)
        {
            gsMesh<> mesh;
            makeMesh(geometry->basis(), mesh, cell_resolution);
            std::string out = filename + "KnotConfiguration";
            std::cout << "Exporting: " << out << std::endl;
            gsWriteParaview(mesh, out);
        }

        if (physical_knot_configuration)
        {
            gsMesh<> mesh;
            makeMesh(geometry->basis(), mesh, cell_resolution);
            geometry->evaluateMesh(mesh);
            std::string out = filename + "PhysicalKnotConfiguration";
            std::cout << "Exporting: " << out << std::endl;
            gsWriteParaview(mesh, out);
        }

        if (basis)
        {
            std::string out = filename + "Basis";
            std::cout << "Exporting: " << out << std::endl;
            gsWriteParaview(geometry->basis(), out, basis_resolution);
        }

        if (control_points)
        {
            std::string out = filename + "ControlPoints";
            std::cout << "Exporting: " << out << std::endl;
            gsMatrix<> coefs = geometry->coefs().transpose();
            gsWriteParaviewPoints(coefs, out);
        }

        if (control_polygon)
        {
            std::string out = filename + "ControlPolygon";
            std::cout << "Exporting: " << out << std::endl;
            gsMesh<> controlNet;
            geometry->controlNet(controlNet);
            gsWriteParaview(controlNet, out);
        }
    }

    void gs_export_physical_knot_configuration_to_vtk(gs_Data spline,
                                                      char* file,
                                                      int cell_resolution)
    {
        gsGeometry<>* geometry = static_cast< gsGeometry<>* >(spline.ptr);

        gsMesh<> mesh;
        makeMesh(geometry->basis(), mesh, cell_resolution);
        geometry->evaluateMesh(mesh);
        std::string out = file;
        std::cout << "Exporting: " << out << std::endl;
        gsWriteParaview(mesh, out);
    }

    void gs_find_param_of_points(gs_Data spline,
                                 double* pts,
                                 int rowsPts,
                                 int cols,
                                 double** par,
                                 int** found,
                                 int smart_start,
                                 double damping_factor,
                                 int iterations,
                                 double accuracy)
    {
        if (cols != 3 || spline.dim != 3)
        {
            std::cout << "Finding parameters of points is done only for volumes\n"
                "Aborting..."
                      << std::endl;
            return;
        }
        gsGeometry<>* geometry = static_cast< gsGeometry<>* >(spline.ptr);

        gsMatrix<> xyz;
        fill_matrix(xyz, pts, rowsPts, cols);
        gsMatrix<> pars;
        fill_matrix(pars,*par,rowsPts,cols);
        gsVector<int> success;
        fill_vector(success,*found,rowsPts);

        gsVector<real_t> value(3),arg(3);
        int result;
        for(int i = 0;i<xyz.rows();++i)
        {
            if(success(i)==1)
                continue;
            value(0)=xyz(i,0);
            value(1)=xyz(i,1);
            value(2)=xyz(i,2);
            arg(0)=0.5;
            arg(1)=0.5;
            arg(2)=0.5;
            result = geometry->newtonRaphson(value,arg,true,accuracy,iterations);
            success(i)= result>=0 ? 1 : 0;
        }

        (*par) = make_c_array(pars);
        (*found) = make_c_array_int(success);
    }

    void gs_find_param_of_point(gs_Data spline,
                                double pt_x,
                                double pt_y,
                                double pt_z,
                                double start_x,
                                double start_y,
                                double start_z,
                                double* par_x,
                                double* par_y,
                                double* par_z,
                                int* found,
                                double damping_factor,
                                int iterations,
                                double accuracy)
    {
        // make sure that (cols == 3 && spline.dim == 3)
        gsGeometry<>* geometry = static_cast< gsGeometry<>* >(spline.ptr);
        gsVector<real_t> value(3),arg(3);
        value(0)=pt_x;
        value(1)=pt_y;
        value(2)=pt_z;
        arg(0)=start_x;
        arg(1)=start_y;
        arg(2)=start_z;
        *found = geometry->newtonRaphson(value,arg,true,accuracy,iterations);
        *par_x=arg(0);
        *par_y=arg(1);
        *par_z=arg(2);
    }


    // ================================================================================
    // delete below

    void gs_transform_solid(PK_BODY_t body,
                            gs_Data transformation,
                            double tolerance)
    {
        gsGeometry<>* transf = static_cast< gsGeometry<>* >(transformation.ptr);
        transform_solid(body, transf, tolerance);
    }

    // delete above
    // ================================================================================

}


