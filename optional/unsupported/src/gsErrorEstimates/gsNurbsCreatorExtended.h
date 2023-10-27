#pragma once

#include <gismo.h>
#include <gsCore/gsForwardDeclarations.h>

namespace gismo {

/**
   @brief Class gsNurbsCreatorExtended provides an extension to the exiting gsNurbsCreator class

   \ingroup Nurbs
*/

template<class T>
struct gsNurbsCreatorExtended : public gsNurbsCreator<T> {

    friend  struct gsNurbsCreator<T>;

    typedef memory::unique_ptr<gsGeometry<T> >  GeometryPtr;
    typedef typename gsTensorBSpline<2,T>::uPtr TensorBSpline2Ptr;
    typedef typename gsTensorNurbs<2,T>::uPtr   TensorNurbs2Ptr;
    typedef typename gsTensorNurbs<3,T>::uPtr   TensorNurbs3Ptr;

    static TensorNurbs2Ptr NurbsRectangleCurved(T low_x = 0, T low_y = 0, T upp_x = 1, T upp_y = 1) {
        gsKnotVector<T> KVx(low_x, upp_x, 0, 2);
        gsKnotVector<T> KVy(low_y, upp_y, 0, 3);

        gsMatrix<T> C(6, 2);

        C << low_x, low_y,
                upp_x, low_y,
                upp_x / 4, upp_y / 2,
                3 * upp_x / 4, upp_y / 2,
                low_x, upp_y,
                upp_x, upp_y;


        // Set weights
        gsMatrix<T> ww(6, 1);
        ww.setOnes();
        ww.at(2) = 0.707106781186548;
        ww.at(3) = 0.707106781186548;

        return TensorNurbs2Ptr(new gsTensorNurbs<2, T>(KVx, KVy, give(C), give(ww)));
    }

    static GeometryPtr BSplineRectangleCurved(int const & deg = 2) {

        GeometryPtr quann = gsNurbsCreatorExtended<T>::NurbsRectangleCurved();

        gsKnotVector<T> KV1(0, 1, 0, 2);
        gsKnotVector<T> KV2(0, 1, deg - 2, deg + 1);

        gsTensorBSplineBasis<2, T> tbsp(new gsBSplineBasis<T>(KV1), new gsBSplineBasis<T>(KV2));
        const gsMatrix<T> pts = tbsp.anchors();
        return tbsp.interpolateData(quann->eval(pts), pts);
    }

    static TensorNurbs3Ptr NurbsCubeCurved(T low_x = 0, T low_y = 0, T low_z = 0, T upp_x = 1, T upp_y = 1, T upp_z = 1) {
        gsKnotVector<T> KVx(low_x, upp_x, 0, 2);
        gsKnotVector<T> KVy(low_y, upp_y, 0, 2);
        gsKnotVector<T> KVz(low_z, upp_z, 0, 3);

        gsMatrix<T> C(12, 3);

        C << low_x, low_y, low_z,
                upp_x, low_y, low_z,
                low_x, upp_y, low_z,
                upp_x, upp_y, low_z,
                upp_x / 4, upp_y / 4, upp_z / 2,
                3 * upp_x / 4, upp_y / 4, upp_z / 2,
                upp_x / 4, 3 * upp_y / 4, upp_z / 2,
                3 * upp_x / 4, 3 * upp_y / 4, upp_z / 2,
                low_x, low_y, upp_z,
                upp_x, low_y, upp_z,
                low_x, upp_y, upp_z,
                upp_x, upp_y, upp_z;

        // Set weights
        gsMatrix<T> ww(12, 1);
        ww.setOnes();
        ww.at(4) = 0.707106781186548;
        ww.at(5) = 0.707106781186548;
        ww.at(6) = 0.707106781186548;
        ww.at(7) = 0.707106781186548;

        return TensorNurbs3Ptr(new gsTensorNurbs<3, T>(KVx, KVy, KVz, give(C), give(ww)));
    }

    static GeometryPtr BSplineCubeCurved(int const & deg = 2) {
        GeometryPtr quann = gsNurbsCreatorExtended<T>::NurbsCubeCurved();

        gsKnotVector<T> KV1(0, 1, 0, 2);
        gsKnotVector<T> KV2(0, 1, 0, 2);
        gsKnotVector<T> KV3(0, 1, deg - 2, deg + 1);

        gsTensorBSplineBasis<3, T> tbsp(new gsBSplineBasis<T>(KV1),
                                        new gsBSplineBasis<T>(KV2),
                                        new gsBSplineBasis<T>(KV3));
        const gsMatrix<T> pts = tbsp.anchors();
        return tbsp.interpolateData(quann->eval(pts), pts);
    }

    // Rectangle described by the identity mapping over the given parameter domain, using tensor product B-splines.
    static TensorBSpline2Ptr BSplineRectangleWithIncreasedMultiplicity( T low_x, T low_y, T upp_x, T upp_y)
    {
        // [low_x, low_y], 3 interior knots, multiplicity 3
        gsKnotVector<T> KVx (0, 1, 0, 3); // {0, 0, 0, 1, 1, 1}         => 6 - p(2) - 1 = 3
        gsKnotVector<T> KVy (0, 1, 0, 3); // {0, 0, 0, 0.5, 1, 1, 1}    =>

        //gsInfo << "KVy: \n";
        //KVy.print(gsInfo);
        //gsInfo << "\n\n";

        KVy.insert(0.5, 2);

        //gsInfo << "KVy: \n";
        //KVy.print(gsInfo);
        //gsInfo << "\n\n";

        // {0, 0, 0, 0.5, 0.5, 1, 1, 1}    =>    8 - p(2) - 1 = 5
        // we need 7 x 4 control points
        gsMatrix<T> C(15,2);
        C << low_x,      low_y,
                0.50*upp_x, low_y,
                upp_x,      low_y,
                low_x,      0.25*upp_y,
                0.50*upp_x, 0.25*upp_y,
                upp_x,      0.25*upp_y,
                low_x,      0.50*upp_y,
                0.5*upp_x,  0.50*upp_y,
                upp_x,      0.50*upp_y,
                low_x,      0.75*upp_y,
                0.50*upp_x, 0.75*upp_y,
                upp_x,      0.75*upp_y,
                low_x,      upp_y,
                0.50*upp_x, upp_y,
                upp_x,      upp_y;

        return TensorBSpline2Ptr(new gsTensorBSpline<2,T>(KVx, KVy, give(C)));
    }
};

}