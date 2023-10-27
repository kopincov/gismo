/** @file gsConverterParXML.h

    @brief It allows to convert from .par to Gismo's .xml format.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include <gismo.h>

namespace gismo
{

/// Counts the number of lines in \a filename.
size_t countLines(const std::string& filename)
{ 
    std::ifstream fileIn;
    fileIn.open(filename);

    std::string line;
    size_t nl=0;
    while(std::getline(fileIn, line))
    {
        nl++;
    }
    fileIn.close();
    return nl;
}

/// Reads \a filename (a .par file) and saves its contents to \a parameters and \a points.
template <class T>
void fillTheMatrices(const std::string& filename,
                     gsMatrix<T>& parameters,
                     gsMatrix<T>& points)
{
    size_t nl=0;
    nl = countLines(filename);
    
    gsInfo << "The current file has " << nl << " lines." << std::endl;
  
    parameters.resize(2,nl);
    points.resize(3,nl);

    std::ifstream fileIn;
    fileIn.open(filename);
    real_t number;
    size_t i=0;
    std::string line;
    while(std::getline(fileIn, line))
    {
        std::istringstream myStream(line);
    
        myStream >> number;
        parameters(0,i)=number;
        myStream >> number;
        parameters(1,i)=number;

        myStream >> number;
        points(0,i)=number;
        myStream >> number;
        points(1,i)=number;
        myStream >> number;
        points(2,i)=number;

        i++;
    }
    fileIn.close();
}

/// Converts a parametrized point cloud from \a input (.par file) to \a output (.xml).
void parToXML(const std::string& input,
              const std::string& output
             )
{
    gsMatrix<> myMatrixPoints;
    gsMatrix<> myMatrixPar;
    fillTheMatrices(input, myMatrixPar, myMatrixPoints);
    
    gsInfo << "the parametrized matrix has " << myMatrixPar.rows() << " rows and " << myMatrixPar.cols() << " cols." << std::endl;
    
    gsInfo << "the poitns matrix has " << myMatrixPoints.rows() << " rows and " << myMatrixPoints.cols() << " cols." << std::endl;
    
    gsFileData<> out;
    
    out.add(myMatrixPar);
    out.add(myMatrixPoints);
  
    out.save(output);
}

/// Reads a parametrized point cloud from \a input (.xml) and saves them to \a parameters and \a points.
template <class T>
void readXML(const std::string& input,
             gsMatrix<T>& parameters,
             gsMatrix<T>& points)
{
    gsFileData<> fd(input);
    fd.getId< gsMatrix<> >(0, parameters);
    fd.getId< gsMatrix<> >(1, points);

    gsInfo << "Parameters has " << parameters.rows() << " rows and " << parameters.cols() << " cols." << std::endl;
    size_t nl=parameters.rows() + points.rows();
    gsInfo << "Value of nl: " << nl << std::endl;
}

/// Writes \a parameters and \a points into \a output (an .xml file).
template <class T>
void writeMatrixintoFile(const std::string& output,
                         gsMatrix<T>& parameters,
                         gsMatrix<T>& points)
{
    std::ofstream filestream;
    filestream.precision(16);
    filestream.open (output);
    for(index_t i=0; i < parameters.cols(); i++)
    {
        filestream << parameters(0,i)<< " " << parameters(1,i) << " " << points(0,i) << " " << points(1,i) << " " << points(2,i) << "\n";
    }
    filestream.close();
}

/// Converts \a input (.xml) into \a output (.par).
void XMLToPar(const std::string& input,
              const std::string& output)
{
    gsMatrix<> parameters, points;
    readXML(input, parameters, points);
    writeMatrixintoFile(output, parameters, points);
}

/// Maximum value in the column \a cl of the \a matrix.
template <class T>
T gs_maximum_value(gsMatrix<T>& matrix, index_t cl)
{
    if(matrix.size() == 0)
        return 0;

    if( cl >= matrix.cols())
    {
        gsInfo << "Invalid choice of column" << "\n";
        return 0;
    }
    
    T maximum = matrix(0,cl);
  
    for(int j = 1; j < matrix.rows(); j++)
    {
        if(maximum < matrix(j,cl))
            maximum = matrix(j,cl);
    }
    return maximum;
}

/// Minimum value in the column \a cl of \a matrix.
template <class T>
T gs_minimum_value(gsMatrix<T>& matrix,
                   index_t cl)
{
    if(matrix.size() == 0)
        return 0;
  
    if( cl >= matrix.cols())
    {
        gsInfo << "Invalid choice of column" << "\n";
        return 0;
    }
    
    T minimum = matrix(0,cl);

    for(index_t j = 1; j < matrix.rows(); j++)
    {
        if(minimum > matrix(j,cl))
            minimum = matrix(j,cl);
    }
    return minimum;
}

/// Scales \a cn so that \a n_min corresponds to 0 and \a n_max to 1.
real_t scalingPoints(real_t cn,
                     real_t n_min,
                     real_t n_max
                    )
{
  return (cn-n_min)/(n_max-n_min);
}


/// Reads the points from \a input (a .par file) and saves them to \a output (another .par file) scaled
/// so that the parameters are between 0 and 1.
void opt_scale(const std::string& input,
               const std::string& output)
{
    std::ifstream fileIn;
    fileIn.open(input);
    if(!fileIn.is_open())
    {
        gsInfo << "File not found." << "\n";
        return;
    }
    
    gsMatrix<real_t> parameters;
    gsMatrix<real_t> points;
    fillTheMatrices(input, parameters, points);
    
    // Because of gs_minimum_value and gs_maximum_value.
    parameters = parameters.transpose();
    points = points.transpose();
    
    gsInfo << "Size of parameters : " << parameters.rows() << " x " << parameters.cols() << "\n";
    gsInfo << "Size of points : " << points.rows() << " x " << points.cols() << "\n";

    
    real_t u_min = gs_minimum_value(parameters, 0);
    real_t u_max = gs_maximum_value(parameters, 0);
    real_t v_min = gs_minimum_value(parameters, 1);
    real_t v_max = gs_maximum_value(parameters, 1);
    
    gsInfo << "minimum u : " << u_min << "\n";
    gsInfo << "maximum u : " << u_max << "\n";
    gsInfo << "minimum v : " << v_min << "\n";
    gsInfo << "maximum v : " << v_max << "\n";

	
    bool scalePar = true;
    if(u_min == 0 && u_max == 1 && v_min == 0 && v_max == 1)
    scalePar = false;
    
    if(scalePar)
    {
        for( int s = 0; s < parameters.rows(); s++)
        {
            parameters(s,0) = scalingPoints( parameters(s,0), u_min, u_max );
            parameters(s,1) = scalingPoints( parameters(s,1), v_min, v_max );
        }
    }
    
    parameters = parameters.transpose();
    points = points.transpose();
    writeMatrixintoFile(output, parameters, points);
}

/// Rescale the data we given in \a input (.par), scales them so that the parametrization is between 0 and 1 and saves the result to \a output (.par).
void scaled_Parameters(const std::string& input,
                       const std::string& output)
{
    std::ifstream fileIn;
    fileIn.open(input);
    if(!fileIn.is_open())
    {
      gsInfo << "File not found." << "\n";
      return;
    }
    real_t u;
    real_t v;
    real_t x;
    real_t y;
    real_t z;
  
    real_t u_min;
    real_t u_max;
    real_t v_min;
    real_t v_max;
  
    bool scalePar = true;
  
    std::ofstream fileOut;
    fileOut.open(output);
    std::string line;
  
    std::vector<real_t> u_par;
    std::vector<real_t> v_par;

    // Read the whole file in to get u_min, u_max, v_min and v_max.
    while(std::getline(fileIn, line))
    {
       std::istringstream myStream(line);
       myStream.precision(16);
       myStream >> u;
       myStream >> v;
       myStream >> x;
       myStream >> y;
       myStream >> z;
     
       u_par.push_back(u);
       v_par.push_back(v);
    }
  
    fileIn.close();
  
    u_min = *std::min_element(u_par.begin(), u_par.end());
    u_max = *std::max_element(u_par.begin(), u_par.end());
  
    v_min = *std::min_element(v_par.begin(), v_par.end());
    v_max = *std::max_element(v_par.begin(), v_par.end());
  
  
  
    if(u_min ==0 && u_max ==1 && v_min ==0 && v_max == 1)
        scalePar = false;

    // Now re-read the file and scale the parameters right away.
    if(scalePar)
    {
        fileIn.open(input);
        std::ofstream fileOut;
        fileOut.open(output);
        std::string line;

        while(std::getline(fileIn, line))
        {
            std::istringstream myStream(line);
            myStream.precision(16);
            myStream >> u;
            myStream >> v;
            myStream >> x;
            myStream >> y;
            myStream >> z;
      
            u = scalingPoints(u, u_min, u_max);
            v = scalingPoints(v, v_min, v_max);
      
            fileOut.precision(16);
            fileOut << u << " " << v << " " << x << " " << y << " " << z << "\n";
        }
    }

    fileIn.close();
    fileOut.close();
}
    
/// Reads all points from \a input that have u between \a u_min and \a u_max and v between \a v_min and \a v_max and writes them to \a output (par).
void restrictedPar(const std::string& input,
                   const std::string& output,
                   real_t u_min,
                   real_t u_max,
                   real_t v_min,
                   real_t v_max,
                   bool scalePar = false,
                   real_t scaling = 1)
{
   
    if(u_min ==0 && u_max ==1 && v_min ==0 && v_max == 1){
      scalePar = false;
    }
  
    std::ifstream fileIn;
    fileIn.open(input);
    if(!fileIn.is_open())
    {
        gsInfo << "File not found." << "\n";
        return;
    }
    real_t u;
    real_t v;
    real_t x;
    real_t y;
    real_t z;
    int nl=0;
    std::ofstream fileOut;
    fileOut.open(output);
  
    std::string line;
    while(std::getline(fileIn, line))
    {
        std::istringstream myStream(line);
        myStream.precision(16);

        myStream >> u;
        if(u_min <= u && u <= u_max)
        {
            myStream >> v;
            if(v_min <= v && v <= v_max)
            {
                myStream >> x;
                myStream >> y;
                myStream >> z;

                if(scalePar)
                {
                    u = scalingPoints(u, u_min, u_max);
                    v = scalingPoints(v,v_min,v_max);
                }
                x *= scaling;
                y *= scaling;
                z *= scaling;

                fileOut << std::setprecision(16);
                fileOut << u << " " << v << " " << x << " " << y << " " << z << "\n";
                nl++;
            }
            else
            {
                myStream.clear();
            }
        }
    }
    fileIn.close();
    fileOut.close();
    gsInfo << "There are " << nl << " points.";
}

/** @brief Reads the parametrized data points from \a input (a .par file)
 * and writes to \a output (another .par file) those, whose parameter in
 * the direction \a dir (false = u, true = v) is not between \a left and \a right.
 */
void cutting_interval(const std::string& input,
                      const std::string& output,
                      real_t& left,
                      real_t& right,
                      bool dir
                     )
{
    
    std::ifstream fileIn;
    fileIn.open(input);
    if(!fileIn.is_open())
    {
        gsInfo << "File not found." << "\n";
        return;
    }
    
    real_t u;
    real_t v;
    real_t x;
    real_t y;
    real_t z;
    
    std::ofstream fileOut;
    fileOut.open(output);
    fileOut << std::setprecision(16);
    std::string line;
     
    while(std::getline(fileIn, line))
    {
        std::istringstream myStream(line);
        myStream.precision(16);

        myStream >> u;
        myStream >> v;
        myStream >> x;
        myStream >> y;
        myStream >> z;

        if(!dir)
        {
            if( (u < left)||(u > right) )
                fileOut << u << " " << v << " " << x << " " << y << " " << z << "\n";
        }
        else
        {
            if( (v < left)||(v > right) )
                fileOut << u << " " << v << " " << x << " " << y << " " << z << "\n";
        }
    }
    
    fileIn.close();
    fileOut.close();
}


/** @brief Reads the parametrized data points from \a input (a .par file)
 * and writes them into \a output (another .par file) with u and v swapped.
 */
void transpose_domain(const std::string& input,
                      const std::string& output)
{
    std::ifstream fileIn;
    fileIn.open(input);
    if(!fileIn.is_open())
    {
     gsInfo << "File not found." << "\n";
     return;
    }
    
    real_t u;
    real_t v;
    real_t x;
    real_t y;
    real_t z;
    
    std::ofstream fileOut;
    fileOut.open(output);
    fileOut << std::setprecision(16);
    std::string line;
     
    while(std::getline(fileIn, line))
    {
        std::istringstream myStream(line);
        myStream.precision(16);
        
        myStream >> u;
        myStream >> v;
        myStream >> x;
        myStream >> y;
        myStream >> z;
        
        fileOut << v << " " << u << " " << x << " " << y << " " << z << "\n";
    }
     
    fileIn.close();
    fileOut.close();
    
}

} // namespace gismo
