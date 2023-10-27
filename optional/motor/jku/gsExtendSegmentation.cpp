/** @file gsExtendSegmentation

    @brief Extends the segmentation. 

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <gismo.h>
#include "gsMotorUtils.h"
#include "gsMotorIOUtils.h"

using namespace gismo;

gsMatrix<> getStartOrEnd(const bool start);
void getExtendedSegmentation(const gsGeometry<>& geom,
                             const gsGeometry<>& offset,
                             const gsMultiPatch<>& segmentation,
                             const real_t tolerance,
                             const bool start,
                             gsMultiPatch<>& extSegmentation);

int main(int argc, char *argv[])
{

    std::string inputCurve(MOTOR_DATA_DIR "jku/target_input.xml");
    std::string inputOffset(MOTOR_DATA_DIR "jku/target_input_offset.xml");
    std::string inputSegmentation(MOTOR_DATA_DIR "jku/target_segmentation.xml");
    std::string outputFile("target_extended_segmentation");
    real_t tolerance = 1e-9;
    bool start = false;

    gsCmdLine cmd("Extending the segmentation.");
    cmd.addString("c", "curveFile", "Name of the Curve's file", inputCurve);
    cmd.addString("f", "offsetFile", "Name of the oFfset's file", inputOffset);
    cmd.addString("s", "segmentationFile", "Name of the Segmentation's file", inputSegmentation);
    cmd.addString("o", "output", "Name of the output file", outputFile);    
    cmd.addReal("t", "tolerance", "The tolerance", tolerance);
    cmd.addSwitch("start", "Sample points from the beginning of the segmentation curves", start);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input arguments:\n"
           << "Input curve:         " << inputCurve << "\n"
           << "Input offset:        " << inputOffset << "\n"
           << "Input segmentation:  " << inputSegmentation << "\n"
           << "Output:              " << outputFile << "\n"
           << "Tolerance:           " << tolerance << "\n"
           << "Start:               " << start << "\n"
           << "--------------------------------------------------\n" << std::endl;

    gsFileData<> fdCurve(inputCurve);
    gsGeometry<>* geom = fdCurve.getAnyFirst< gsGeometry<> >().release();

    gsFileData<> fdOffset(inputOffset);
    gsGeometry<>* offset = fdOffset.getAnyFirst< gsGeometry<> >().release();

    gsFileData<> fdSegmentation(inputSegmentation);
    gsMultiPatch<>* segmentation = fdSegmentation.getAnyFirst< gsMultiPatch<> >().release();

    gsMultiPatch<> extSegmentation;

    getExtendedSegmentation(*geom, *offset, *segmentation, tolerance, start, extSegmentation);

    writeMultipatch(extSegmentation, outputFile);

    return 0;
}

    
gsMatrix<> getStartOrEnd(const bool start)
{
    gsMatrix<> mat(1, 1);

    if (start)
    {
        mat.setOnes();
    }
    else
    {
        mat.setZero();
    }

    return mat;
}

void getExtendedSegmentation(const gsGeometry<>& geom,
                             const gsGeometry<>& offset,
                             const gsMultiPatch<>& segmentation,
                             const real_t tolerance,
                             const bool start,
                             gsMultiPatch<>& extSegmentation)
{
    gsMatrix<> segParam = getStartOrEnd(start);
    
    for (index_t i = 0; i != segmentation.nPatches(); i++)
    {
        const gsGeometry<>& segLine = segmentation[i];
        gsMatrix<> point = segLine.eval(segParam);
        gsMatrix<> param = findParameter(offset, point, tolerance);

        gsMatrix<> points(2, 2);
        gsMatrix<> firstPoint = offset.eval(param);
        points.col(0) = firstPoint.col(0);
        gsMatrix<> secondPoint = geom.eval(param);
        points.col(1) = secondPoint;
        points.transposeInPlace();
        
        gsKnotVector<> kv(0, 1, 0, 2);
        extSegmentation.addPatch(gsBSpline<>(kv, points));
    }
}

