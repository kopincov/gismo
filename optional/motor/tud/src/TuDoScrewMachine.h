#pragma once

#include <gismo.h>



using namespace gismo;
using namespace std;

class TuDoScrewMachine
{
public:
	TuDoScrewMachine();
	~TuDoScrewMachine();

	void				fitting_test();

	void				convertMatrixToSplineCurve(string filename, int numberOfLobesMale, int numberOfLobesFemale, double axialDistance, double rotationAngle, bool fit, int fitPoints);
	void				convertMatrixToSplineSurface(string filename, int numberOfLobesMale, int numberOfLobesFemale, double axialDistance, double length, double wrapAngleMale, double rotationAngle);
	

private:
	gsMatrix<>			closeContour(gsMatrix<> unclosedContour);
	gsMatrix<>			rotateContour(gsMatrix<> closedContour, double rotationAngle, double axialDistance);
	gsMatrix<>			calcCoefs(gsMatrix<> matrixClosed, double wrapAngle, double length, double axialDistance);
	void				buildSurface(gsMatrix<> matrixClosed, double wrapAngle, double length, double axialDistance, string fileName);
	void				buildCurve(gsMatrix<> matrixClosed, string fileName);
	void				fitCurve(gsMatrix<> matrixClosed, int fitPoints, string fileName);
};
