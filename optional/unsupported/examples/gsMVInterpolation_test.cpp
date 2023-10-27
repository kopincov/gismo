#define debug 1

#include <iostream>

#include <gismo.h>

#include <gsModeling/gsModelingUtils.hpp>
#include <gsSegment/gsMVInterpolation.h>
#include <gsSegment/gsVolumeSegment.h>

/*
#include <gsModeling/gsCurveLoop.h>
#include <gsSegment/gsMVInterpolation.h>
#include <gsNurbs/gsBSpline.h>
#include <gsModeling/gsModelingUtils.hpp>
#include <gsModeling/gsTraceCurve.hpp>
#include <gsSegment/gsVolumeSegment.h>
#include <gsIO/gsCmdLine.h>
*/

#include <fstream>


using namespace gismo;



int nGaussCases = 6;
int nGauss[] = {8, 16, 32, 64, 128, 128};
bool regPoly[] = {false, false, false, false, false, true};

gsCurveLoop<> *makeTestCurveLoop()
{
  gsInfo << "\nConstructing a curve loop from points in space\n";
  gsMatrix<> srcPts(6, 3);
  srcPts << 0, 0, 0,
      1, 0, 0,
      1, 0.5, 0,
      0.5, 0.5, 0,
      0.5, 1, 0,
      0, 1, 0;

  std::vector< gsVector3d<> *> points3d;
  for(index_t i = 0; i < srcPts.rows(); i++)
  {
    points3d.push_back(new gsVector3d<>(srcPts.row(i)));
  }
  std::vector<bool> isConvex;

  for (int i = 0; i < 6; i++) isConvex.push_back(i!=3);

  gsCurveLoop<> *curveLoop = new gsCurveLoop<>(points3d, isConvex, 0.1, NULL);

  for(index_t i = 0; i < srcPts.rows(); i++)
  {
    delete points3d[i];
  }
  return curveLoop;
}

gsTrimSurface<> *makeTestSurface(gsCurveLoop<> *curveLoop, gsMatrix<> &controlPoints)
{
  // set up a trimmed surface from the curve loop
  gsPlanarDomain<> * domain= new gsPlanarDomain<>(curveLoop);
  gsKnotVector<> KV1 = gsKnotVector<>(0, 1, 0, 2);//multiplicity at ends==1+1
  gsKnotVector<> KV2 = gsKnotVector<>(0, 1, 0, 2);//multiplicity at ends==1+1
  typedef gsTensorBSplineBasis<2> Basis; //dimension==2
  Basis *basis = new Basis(KV1, KV2);
  gsTensorBSpline<2,real_t>::Ptr masterSpline(new gsTensorBSpline<2,real_t>(*basis, give(controlPoints)));
  gsTrimSurface<> * trimSurface = new gsTrimSurface<>(masterSpline, domain);
  delete basis;
  return trimSurface;
}

gsMVInterpolation<> makeTestMVI(gsTrimSurface<> *trimSurface, const gsMatrix<> &corners, int nGauss2)
{
  gsMVInterpolation<> mvi(trimSurface, corners, nGauss2);

  return mvi;
}

void testCuttingCurve(int idxStart, int idxEnd, real_t w_reg, bool regular,
                      gsTrimSurface<> &trimSurface, gsMatrix<> &fitTimes, gsBSpline<> **imgCurve,
                      gsMatrix<> &fitPoints, gsBSpline<> &fitCurve)
{
  gsInfo << "Finding a cutting curve\n";
  real_t curveFraction = 0.9;
  int nFitPoints = 20;

  gsVolumeSegment<>::chooseCuttingCurve(trimSurface, idxStart, idxEnd,
                                           w_reg, 1.0, nFitPoints, 0.0001,
                                           false, regular, curveFraction, 5, imgCurve, fitCurve, fitPoints);

  fitTimes = gsPointGrid<real_t>((1 - curveFraction) / 2, (1 + curveFraction) / 2, nFitPoints);
}

void outputResults(std::string filename, gsTrimSurface<> &trimSurface,
                    gsMatrix<> &cornersNormal, gsMatrix<> &cornersReg, gsMatrix<> &gridPts,
                    std::vector<gsMatrix<> > &val, std::vector<int> &gridi,
                    std::vector<int> &gridj, int count, const gsMatrix<> fitTimes[], gsBSpline<> *imgCurve[],
                    const gsMatrix<> fitPoints[], const gsBSpline<> fitCurve[], const gsMatrix<> & evalCurveTimes)
{
  std::ofstream fout(filename.c_str());
  fout << std::fixed << std::setprecision(6);

  const gsCurveLoop<> & curveLoop = trimSurface.domain().loop(0);

  fout << "domCor={";
  for(int j = 0; j < curveLoop.size(); j++)
  {
    fout << "{" << curveLoop.curve(j).coefs()(0, 0) << "," << curveLoop.curve(j).coefs()(0, 1) << "},";
  }
  // repeat the first corner
  fout << "{" << curveLoop.curve(0).coefs()(0, 1) << "," << curveLoop.curve(0).coefs()(0, 1) << "}};\n";

  fout << "cornersNormal={";
  for(int j = 0; j < cornersNormal.rows(); j++)
  {
    fout << "{" << cornersNormal(j, 0) << "," << cornersNormal(j, 1) << "},";
  }
  fout << "{" << cornersNormal(0, 0) << "," << cornersNormal(0, 1) << "}};\n";

  fout << "cornersReg={";
  for(int j = 0; j < cornersReg.rows(); j++)
  {
    fout << "{" << cornersReg(j, 0) << "," << cornersReg(j, 1) << "},";
  }
  fout << "{" << cornersReg(0, 0) << "," << cornersReg(0, 1) << "}};\n";

  fout << "gridPts={\n";
  for(size_t j = 0; j < gridi.size(); j++)
  {
    fout << "{" << gridi[j] << "," << gridj[j] << "}";
    if(j < gridi.size() - 1) fout << ",\n";
  }
  fout << "};\n";

  fout << "dom={\n";
  for(int j = 0; j < gridPts.cols(); j++)
  {
    fout << "{" << gridPts(0, j) << "," << gridPts(1, j) << "}";
    if(j < gridPts.cols() - 1) fout << ",\n";
  }
  fout << "};\n";

  GISMO_ASSERT(val.size() == (size_t)nGaussCases, "Array size does not match expected number of cases");
  fout << "nGaussCases={\n";
  for(int idxV = 0; idxV < nGaussCases; idxV++)
  {
    fout << nGauss[idxV];
    if(idxV < nGaussCases - 1) fout << ",";
  }
  fout << "};\n\n";

  GISMO_ASSERT(val.size() == (size_t)nGaussCases, "Array size does not match expected number of cases");
  fout << "regPoly={\n";
  for(int idxV = 0; idxV < nGaussCases; idxV++)
  {
    fout << regPoly[idxV];
    if(idxV < nGaussCases - 1) fout << ",";
  }
  fout << "};\n\n";

  fout << "val={\n";
  for(size_t idxV = 0; idxV < val.size(); idxV++)
  {
    fout << "{";
    for(int j = 0; j < val[idxV].cols(); j++)
    {
      fout << "{" << val[idxV](0, j) << "," << val[idxV](1, j) << "}";
      if(j < val[idxV].cols() - 1) fout << ",\n";
    }
    fout << "}";
    if(idxV < val.size() - 1) fout << ",";
  }
  fout << "};\n\n";

  fout << "fitTimes={\n";
  for(index_t i = 0; i < count; i++)
  {
    GISMO_ASSERT(fitTimes[i].rows() == 1, "Unexpected matrix of times");
    fout << "{";
    for(index_t idxT = 0; idxT < fitTimes[i].cols(); idxT++)
    {
      fout << fitTimes[i](0, idxT);
      if(idxT < fitTimes[i].cols() - 1) fout << ",";
    }
    fout << "}";
    if(i < count - 1) fout << ",\n";
  }
  fout << "};\n\n";

  fout << "imgValues={\n";
  for(index_t i = 0; i < count; i++)
  {
    fout << "{";
    for(index_t idxT = 0; idxT < fitTimes[i].cols(); idxT++)
    {
      gsMatrix<> fitValue = imgCurve[i]->eval(fitTimes[i].block(0, idxT, 1, 1));
      GISMO_ASSERT(fitValue.rows() == 2 && fitValue.cols() == 1, "Unexpected result matrix dimensions");
      fout << "{" << fitValue(0, 0) << "," << fitValue(1, 0) << "}";
      if(idxT < fitTimes[i].cols() - 1) fout << ",";
    }
    fout << "}";
    if(i < count - 1) fout << ",\n";
  }
  fout << "};\n\n";

  fout << "fitPoints={\n";
  for(index_t i = 0; i < count; i++)
  {
    GISMO_ASSERT(fitPoints[i].rows() == 2, "Unexpected matrix of points");
    fout << "{";
    for(index_t idxT = 0; idxT < fitPoints[i].cols(); idxT++)
    {
      fout << "{" << fitPoints[i](0, idxT) << "," << fitPoints[i](1, idxT) << "}";
      if(idxT < fitPoints[i].cols() - 1) fout << ",\n";
    }
    fout << "}";
    if(i < count - 1) fout << ",\n";
  }
  fout << "};\n\n";

  fout << "ctrlPoints={\n";
  for(index_t i = 0; i < count; i++)
  {
    gsMatrix<> ctrlPts = fitCurve[i].coefs();
    GISMO_ASSERT(ctrlPts.cols() == 2, "Unexpected control point matrix");
    fout << "{";
    for(index_t j = 0; j < ctrlPts.rows(); j++)
    {
      fout << "{" << ctrlPts(j, 0) << "," << ctrlPts(j, 1) << "}";
      if(j < ctrlPts.rows() - 1) fout << ",\n";
    }
    fout << "}";
    if(i < count - 1) fout << ",\n";
  }
  fout << "};\n\n";

  fout << "resultValues={\n";
  for(index_t i = 0; i < count; i++)
  {
    fout << "{";
    for(index_t idxT = 0; idxT < evalCurveTimes.cols(); idxT++)
    {
      gsMatrix<> valMatrix = fitCurve[i].eval(evalCurveTimes.block(0, idxT, 1, 1));
      GISMO_ASSERT(valMatrix.rows() == 2 && valMatrix.cols() == 1, "Unexpected result matrix dimensions");
      fout << "{" << valMatrix(0, 0) << "," << valMatrix(1, 0) << "}";
      if(idxT < evalCurveTimes.cols() - 1) fout << ",\n";
    }
    fout << "}";
    if(i < count - 1) fout << ",\n";
  }
  fout << "};\n\n";
}

int main(int argc, char * argv[])
{
  index_t numGridSpaces = 10;
  gsCmdLine cmd("Test mean value interpolation.");
  cmd.addInt("s","samples", "Number of grid points per axis", numGridSpaces);
  
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  GISMO_ASSERT(sizeof(nGauss) / sizeof(int) == nGaussCases, "Array size does not match expected number of cases");

  gsCurveLoop<> *curveLoop = makeTestCurveLoop();
  gsMatrix<> pcp(4, 3);
  pcp << 0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      1, 1, 0;

  gsTrimSurface<> *trimSurface = makeTestSurface(curveLoop, pcp);
  gsMatrix<> cornersNormal, cornersReg;
  gsVolumeSegment<>::chooseConvexPolygon(*trimSurface, false, cornersNormal);
  gsVolumeSegment<>::chooseConvexPolygon(*trimSurface, true, cornersReg);

  gsInfo << "Finding mean value interpolant at a point\n";
  gsMatrix<> gridPts(2, (numGridSpaces + 1) * (numGridSpaces + 1));
  std::vector<int> gridi, gridj;
  index_t thisVecIdx = 0;
  for(int j = 0; j <= numGridSpaces; j++)
  {
    for(int i = 0; i <= numGridSpaces; i++)
    {
      real_t x = ((real_t)i) / numGridSpaces;
      real_t y = ((real_t)j) / numGridSpaces;
      if((x < 0.5 || y < 0.5) && x > 0 && y > 0 && x < 1 && y < 1) // test that this point is inside the L-shape
      {
        gridPts(0, thisVecIdx) = x;
        gridPts(1, thisVecIdx) = y;
        gridi.push_back(i);
        gridj.push_back(j);
        thisVecIdx++;
      }
    }
  }
  gridPts = gridPts.block(0, 0, 2, thisVecIdx).eval();
  gsInfo << "Generated " << gridPts.cols() << " points in the domain.\n";

  std::vector<gsMatrix<> > val;

  for(int idxGauss = 0; idxGauss < nGaussCases; idxGauss++)
  {
    gsInfo << "Performing MVI with " << nGauss[idxGauss] << " points per curve. (regular = " << regPoly[idxGauss] << ")\n" << std::flush;
    gsMVInterpolation<> mvi = makeTestMVI(trimSurface, regPoly[idxGauss]?cornersReg:cornersNormal, nGauss[idxGauss]);

    gsMatrix<> newVal;
    mvi.eval_into(gridPts, newVal);
    val.push_back(newVal);
  }

  const int curveCount = 7;
  gsMatrix<> fitTimes[curveCount];
  gsMatrix<> fitPoints[curveCount];
  gsBSpline<> fitCurve[curveCount];
  gsBSpline<> *imgCurve[curveCount];
  int i = 0;
  testCuttingCurve(1, 5, 0.0, false, *trimSurface, fitTimes[i], &imgCurve[i], fitPoints[i], fitCurve[i]); i++;
  testCuttingCurve(1, 5, 0.01, false, *trimSurface, fitTimes[i], &imgCurve[i], fitPoints[i], fitCurve[i]); i++;
  testCuttingCurve(1, 5, 0.1, false, *trimSurface, fitTimes[i], &imgCurve[i], fitPoints[i], fitCurve[i]); i++;
  testCuttingCurve(2, 4, 0.0, false, *trimSurface, fitTimes[i], &imgCurve[i], fitPoints[i], fitCurve[i]); i++;
  testCuttingCurve(2, 4, 0.01, false, *trimSurface, fitTimes[i], &imgCurve[i], fitPoints[i], fitCurve[i]); i++;
  testCuttingCurve(2, 4, 0.1, false, *trimSurface, fitTimes[i], &imgCurve[i], fitPoints[i], fitCurve[i]); i++;
  testCuttingCurve(0, 3, 0.0, false, *trimSurface, fitTimes[i], &imgCurve[i], fitPoints[i], fitCurve[i]); i++;
  GISMO_ASSERT(i == curveCount, "Incorrect number of cases");

  // Uncomment to output the result to txt file
  // gsMatrix<> evalCurveTimes = *gsPointGrid<real_t>(0, 1, 50);
  // outputResults("gsmvinb.txt", *trimSurface, cornersNormal, cornersReg, gridPts, val, gridi, gridj, curveCount, fitTimes, imgCurve, fitPoints, fitCurve, evalCurveTimes);

  for(int j = 0; j < curveCount; j++) delete imgCurve[j];
  delete trimSurface;
}
