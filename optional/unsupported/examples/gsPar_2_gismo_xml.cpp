/** @file gsPar_2_gismo_xml.cpp

    @brief Converts a parametrised point cloud in .par format to a
    gismo xml file. Based on gsUV_XYZ_2_gismo_xml.cpp by C. Bracco,
    C. Giannelli and J. Speh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
*/

#include <fstream>
#include <vector>

#include <gismo.h>

using namespace gismo;


/**
   @brief Read the parameters and points from a par file.
   @param[in] filename Name of the input par file.
   @param[out] parameters Parameters read from @filename.
   @param[out] points Points read from @a filename.
   @return true iff successful.
*/
bool read_par(const std::string& filename, gsMatrix<>& parameters, gsMatrix<>& points)
{
    std::fstream stream;
    stream.open(filename.c_str());
    if(!stream)
    {
	gsWarn << "Failed to open the file " << filename << ", exiting.\n";
	return false;
    }	    
    
    real_t u, v, x, y, z;
    std::vector<real_t> vParameters;
    std::vector<real_t> vPoints;
    
    int nrLines = 0;
    while (stream >> u >> v >> x >> y >> z)
    {
	++nrLines;
	vParameters.push_back(u);
        vParameters.push_back(v);

	vPoints.push_back(x);
	vPoints.push_back(y);
	vPoints.push_back(z);
    }

    std::cout << "Read " << nrLines << " lines of input." << std::endl;

    if((3 * vParameters.size()) != (2 * vPoints.size()))
	gsWarn << "Read " << vParameters.size() / 2 << " parameters and\n"
	       << " and " << vPoints.size() / 3 << "points, \n"
	       << " i.e., the sizes are different and weird things might happen." << "\n\n";

    int num_cols = vPoints.size() / 3;
    parameters.resize(2, num_cols);
    points.resize(3, num_cols);
    // Warning: This way, the elements are not erased if the size of
    // the matrices is correct in the beginning. But as we fill all
    // the values afterwards anyway, this should not be a problem.

    for (int col = 0; col != num_cols; col++)
    {
        parameters(0, col) = vParameters[2 * col + 0];
        parameters(1, col) = vParameters[2 * col + 1];

        points(0, col) = vPoints[3 * col + 0];
        points(1, col) = vPoints[3 * col + 1];
        points(2, col) = vPoints[3 * col + 2];
    }
    stream.close();
    return true;
}



int main(int argc, char *argv[])
{
    std::string input("");
    std::string output("fitting_data");
    
    // ----------------------------------------------------------------------

    gsCmdLine cmd("Converting a parameterized point cloud from the .par format to G+Smo .xml.");

    cmd.addString("i", "input", "File in the .par format with the data.", input);
    cmd.addString("o", "output", "Output file", output);
    
    try { cmd.getValues(argc,argv); } catch( int rv ) { return rv; }
    
    // ----------------------------------------------------------------------
    
    gsInfo << "----------------------------------------\n"
           << "INPUT DATA:\n\n" 
           << "input:   " << input  << "\n\n"
           << "output:  " << output << "\n\n"
           << "----------------------------------------" << std::endl;

    // ----------------------------------------------------------------------
    
    if (input == "")
    {
        gsInfo << "No input. Exiting...\n";
        return 0;
    }
    
    gsMatrix<> uv;
    gsMatrix<> xyz;

    if( read_par(input, uv, xyz) )
    {
        gsFileData<> fd;
	fd << uv;
	fd << xyz;

	gsInfo << "Output: " << output << std::endl;
	fd.dump(output);

	return 0;
    }
    else
    {
	gsWarn << "Reading par file not successful, exiting." << std::endl;
	return -1;
    }
}

