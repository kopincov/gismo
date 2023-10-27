/** @file gsVTK.h

    @brief .vtk file converter

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S.Imperatore
*/

#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>

#include <gismo.h> // enough if only the stable (public) part of gismo is used.


namespace gismo
{

template <class T>
void heading_vtk(const gsMatrix<T>& matrix,
                 std::ofstream& fstream)
{
    
    fstream << "# vtk DataFile Version 3.0\n"
             << "Really cool data\n"
             << "ASCII\n"
             << "DATASET POLYDATA\n"
             << "POINTS " << matrix.cols() << " float\n";
}


template <class T>
void points_vtk(const gsMatrix<T>& matrix,
                std::ofstream& fstream
               )
{
    for( int j = 0; j < matrix.cols(); j++ )
    {
        for(index_t i = 0; i < matrix.rows(); i++ )
            fstream << matrix(i,j) << " ";
        
        fstream << "\n";
    }
}


template <class T>
void vertices_vtk(const gsMatrix<T>& matrix,
                  std::ofstream& fstream
                 )
{
    fstream << "VERTICES " << matrix.cols() << " " << 2*matrix.cols() << "\n";
    
    for(index_t i = 0; i < matrix.cols(); i++)
        fstream << 1 << " " << i << "\n";

}

template <class T>
void attributes_vtk(const std::vector<T>& info,
                    std::ofstream& fstream
                   )
{
    fstream << "POINT_DATA " << info.size() << "\n" << "SCALARS info double" << "\n" << "LOOKUP_TABLE default" << "\n";
    for( size_t s = 0; s < info.size(); s++ )
        fstream << info[s] << "\n";    
}

template <class T>
void vtk_converter(const gsMatrix<T>& matrix,
                   std::vector<T>& info,
                   const std::string& output
                  )
{
    std::ofstream fstream;
    fstream.precision(16);
    fstream.open(output);
    heading_vtk(matrix, fstream);
    points_vtk(matrix, fstream);
    vertices_vtk(matrix, fstream);
    attributes_vtk(info, fstream);
    
    fstream.close();
    
}

} // namespace gismo
