#include <TuDoScrewMachine.h>


TuDoScrewMachine::TuDoScrewMachine()
{
}

TuDoScrewMachine::~TuDoScrewMachine()
{
}

//this function loads a matrix from the xml-file and converts it to a spline curve
void TuDoScrewMachine::convertMatrixToSplineCurve(string filename, int numberOfLobesMale, int numberOfLobesFemale, double axialDistance, double rotationAngle, bool fit, int fitPoints)
{
	//int numberOfLobesMale = 4; //constant for this profile
	//int numberOfLobesFemale = 6; 
	//double axialDistance = 22.608000 + 33.912000;


	//load the matrices from file
	//gsFileData <> fdIn("SRM4+6.xml");
	gsFileData <> fdIn(filename);
	gsMatrix<> maleFromFile = *fdIn.getId<gsMatrix<> >(0);
	gsMatrix<> femaleFromFile = *fdIn.getId<gsMatrix<> >(1);
	gsMatrix<> casingFromFile = *fdIn.getId<gsMatrix<> >(2);

	//transpose the matrices
	gsMatrix<> maleTransposed = maleFromFile.transpose();
	gsMatrix<> femaleTransposed = femaleFromFile.transpose();
	gsMatrix<> CasingTransposed = casingFromFile.transpose();

	//close the matrices (casing is already closed)
	gsMatrix<> maleClosed = closeContour(maleTransposed);
	gsMatrix<> femaleClosed = closeContour(femaleTransposed);

	//rotate to the desired rotor position
	gsMatrix<> maleClosedRotated = rotateContour(maleClosed, rotationAngle, 0.0);
	gsMatrix<> femaleClosedRotated = rotateContour(femaleClosed, -rotationAngle*numberOfLobesMale / numberOfLobesFemale, axialDistance);

	if (fit == false)
	{
		//build the curves with all points initially given
		buildCurve(maleClosedRotated, "Male");
		buildCurve(femaleClosedRotated, "Female");
		buildCurve(CasingTransposed, "Casing");
	}
	else
	{
		//fit the curves to a certain number of points
		fitCurve(maleClosedRotated, fitPoints, "Male");
		fitCurve(femaleClosedRotated, fitPoints, "Female");
		fitCurve(CasingTransposed, fitPoints, "Casing");
	}
}

void TuDoScrewMachine::fitCurve(gsMatrix<> matrixClosed, int fitPoints, string fileName)
{
	int pointNumber = matrixClosed.rows();
	gsMatrix<> uMatrix(pointNumber, 1);
	gsKnotVector<> kvRotor(0, 1, pointNumber - 3, 3);
	for (int a = 0; a < pointNumber; a++)
	{
		uMatrix(a, 0) = kvRotor.at(a);
	}

	gsKnotVector<> kvFit(0, 1, fitPoints, 3);
	gsCurveFitting<> thatFit(uMatrix, matrixClosed, kvFit);
	thatFit.compute();
	const gsBSpline<> & curve = thatFit.curve();

	//computes the approximation error
	real_t error;
	thatFit.computeApproxError(error);
	cout << "The approximation error of the fitted curve is: " << error << "\n";

	gsWrite(curve, "FittedSpline_" + fileName);
	gsWriteParaview(curve, "FittedCurve_" + fileName, pointNumber + 1, false, false);
}

void	TuDoScrewMachine::buildCurve(gsMatrix<> matrixClosed, string fileName)
{
	//create the spline
	int pointNumber = matrixClosed.rows();

	gsKnotVector<> kvRotor(0, 1, pointNumber - 3, 3);
	gsBSpline<> rotorSpline(kvRotor, give(matrixClosed));
	gsWrite(rotorSpline, "CurveSpline_" + fileName);
	gsWriteParaview(rotorSpline, "Curve_" + fileName, pointNumber + 1, false, false);
}



//this function loads a matrix from the xml-file and converts it to a spline surface
void TuDoScrewMachine::convertMatrixToSplineSurface(string filename, int numberOfLobesMale, int numberOfLobesFemale, double axialDistance, double length, double wrapAngleMale, double rotationAngle)
{
	//int numberOfLobesMale = 4; //constant for this profile
	//int numberOfLobesFemale = 6; //constant for this profile
	//double axialDistance = 22.608000 + 33.912000; //constant for this profile. Defines the distance of the centers of the rotors (sum of the roll circles)

	double wrapAngleFemale = -numberOfLobesMale*1.0 / numberOfLobesFemale * wrapAngleMale;

	//load the matrices from file
	//gsFileData <> fdIn("SRM4+6.xml");
	gsFileData <> fdIn(filename);
	gsMatrix<> maleFromFile   = *fdIn.getId<gsMatrix<> >(0);
	gsMatrix<> femaleFromFile = *fdIn.getId<gsMatrix<> >(1);
	gsMatrix<> casingFromFile = *fdIn.getId<gsMatrix<> >(2);

	//transpose the matrices
	gsMatrix<> maleTransposed = maleFromFile.transpose();
	gsMatrix<> femaleTransposed = femaleFromFile.transpose();
	gsMatrix<> CasingTransposed = casingFromFile.transpose();

	//close the matrices (casing is already closed)
	gsMatrix<> maleClosed = closeContour(maleTransposed);
	gsMatrix<> femaleClosed = closeContour(femaleTransposed);

	//rotate to the desired rotor position
	gsMatrix<> maleClosedRotated = rotateContour(maleClosed, rotationAngle, 0.0);
	gsMatrix<> femaleClosedRotated = rotateContour(femaleClosed, -rotationAngle*numberOfLobesMale / numberOfLobesFemale, axialDistance);


	//build the surfaces
	buildSurface(maleClosedRotated, wrapAngleMale, length, 0.0, "surface Male");
	buildSurface(femaleClosedRotated, wrapAngleFemale, length, axialDistance, "surface Female");
	buildSurface(CasingTransposed, 0.0, length, 0.0, "surface Casing");
}

gsMatrix<>	TuDoScrewMachine::rotateContour(gsMatrix<> closedContour, double angle, double axialDistance)
{

	gsMatrix<> coefs(closedContour.rows(), 2);
	for (index_t row = 0; row < closedContour.rows(); row++)
	{
        coefs(row, 0) = (cos(angle / 180 * EIGEN_PI)*(closedContour(row, 0) - axialDistance) - sin(angle / 180 * EIGEN_PI)*closedContour(row, 1)) + axialDistance;
        coefs(row, 1) = sin(angle / 180 * EIGEN_PI)*(closedContour(row, 0) - axialDistance) + cos(angle / 180 * EIGEN_PI)*closedContour(row, 1);

	}
	//crossSections.push_back(matrix);
	return coefs;
}

//this function adds the first point of the matrix to its end
gsMatrix<>	TuDoScrewMachine::closeContour(gsMatrix<> unclosedContour)
{

	gsMatrix<> matrixClosed(unclosedContour.rows() + 1, unclosedContour.cols() + 1);
	for (int c = 0; c < unclosedContour.rows(); c++)
	{
		matrixClosed(c, 0) = unclosedContour(c, 0);
		matrixClosed(c, 1) = unclosedContour(c, 1);
		matrixClosed(c, 2) = 0;
	}
	matrixClosed(unclosedContour.rows(), 0) = unclosedContour(0, 0);
	matrixClosed(unclosedContour.rows(), 1) = unclosedContour(0, 1);
	matrixClosed(unclosedContour.rows(), 2) = 0;

	return matrixClosed;
}

//this function calculates the coefficients of the spline
gsMatrix<>	TuDoScrewMachine::calcCoefs(gsMatrix<> matrixClosed, double wrapAngle, double length, double axialDistance)
{
	// ---------------------- Casing (has no wrapping)
	if (wrapAngle == 0)
	{
		int intLength = (int)math::round(length);
		if (intLength < length)
			intLength++;
		gsMatrix<> coefs(matrixClosed.rows()*(intLength + 1), 3);
		for (int w = 0; w <= intLength; w++)
		{
			if (w >= length)
			{
				for (index_t row = 0; row < matrixClosed.rows(); row++)
				{
					coefs(matrixClosed.rows()*w + row, 0) = matrixClosed(row, 0);
					coefs(matrixClosed.rows()*w + row, 1) = matrixClosed(row, 1);
					coefs(matrixClosed.rows()*w + row, 2) = length;
				}
				break;
			}

			for (index_t row = 0; row < matrixClosed.rows(); row++)
			{
				coefs(matrixClosed.rows()*w + row, 0) = matrixClosed(row, 0);
				coefs(matrixClosed.rows()*w + row, 1) = matrixClosed(row, 1);
				coefs(matrixClosed.rows()*w + row, 2) = w;
			}
		}
		return coefs;
	}
	// ----------------------

	// ---------------------- male and female rotor
	double wrappingDirection = 1.0; //defines direction of wrapping: >0 for male, <0 for female
	if (wrapAngle < 0)
	{
		wrappingDirection = -1.0;
	}

	double absWrap = math::abs(wrapAngle);
	int absWrapInt = (int)math::round(absWrap);
	if (absWrapInt < absWrap)
		absWrapInt++;

	gsMatrix<> coefs(matrixClosed.rows()*(absWrapInt + 1), 3);

	//this for loop rotates the cross sections and sets a new z-coordinate
	//for the male rotor axialDistance is 0
	for (int w = 0; w <= absWrapInt; w++)
	{
		//this request is necessary for the case the rotor does not end with a full degree of wrapping
		if (w*1.0 / absWrap * length >= length)
		{
			for (index_t row = 0; row < matrixClosed.rows(); row++)
			{
                coefs(matrixClosed.rows()*w + row, 0) = (cos(absWrap*wrappingDirection / 180 * EIGEN_PI)*(matrixClosed(row, 0) - axialDistance) - sin(absWrap*wrappingDirection / 180 * EIGEN_PI)*matrixClosed(row, 1)) + axialDistance;
                coefs(matrixClosed.rows()*w + row, 1) = sin(absWrap*wrappingDirection / 180 * EIGEN_PI)*(matrixClosed(row, 0) - axialDistance) + cos(absWrap*wrappingDirection / 180 * EIGEN_PI)*matrixClosed(row, 1);
				coefs(matrixClosed.rows()*w + row, 2) = length;
			}
			break;
		}

		for (index_t row = 0; row < matrixClosed.rows(); row++)
		{
            coefs(matrixClosed.rows()*w + row, 0) = (cos(w*wrappingDirection / 180 * EIGEN_PI)*(matrixClosed(row, 0) - axialDistance) - sin(w*wrappingDirection / 180 * EIGEN_PI)*matrixClosed(row, 1)) + axialDistance;
            coefs(matrixClosed.rows()*w + row, 1) = sin(w*wrappingDirection / 180 * EIGEN_PI)*(matrixClosed(row, 0) - axialDistance) + cos(w*wrappingDirection / 180 * EIGEN_PI)*matrixClosed(row, 1);
			coefs(matrixClosed.rows()*w + row, 2) = w*1.0 / absWrap * length;
		}
		//crossSections.push_back(matrix);
	}
	// ----------------------
	return coefs;
}

//this function calculates the surfaces
//the procedure is adopted from a gismo-tutorial
void	TuDoScrewMachine::buildSurface(gsMatrix<> matrixClosed, double wrapAngle, double length, double axialDistance, string fileName)
{
	//Build surface of male rotor
	int n = matrixClosed.rows();
	int degree = 3;

	//construction of coefficients
	gsMatrix<> coefs = calcCoefs(matrixClosed, wrapAngle, length, axialDistance);
	int m = coefs.rows() / n;

	//Output as surface
	gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
	gsKnotVector<> kv2(0, 1, m - degree - 1, degree + 1);

	//construction of a basis
	gsTensorBSplineBasis<2, real_t > basis(kv1, kv2);

	//putting basis and coefficients toghether
	gsTensorBSpline<2, real_t >  surface(basis, coefs);

	//output
	gsWrite(surface, fileName);
	gsWriteParaview(surface, fileName, 200000);
}
