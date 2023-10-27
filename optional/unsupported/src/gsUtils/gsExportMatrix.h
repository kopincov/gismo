/** @file gsExportMatrix.h

    @brief export matrices to ASCII format

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <string>
#include <iostream>
#include <iomanip>
#include <limits>

#include <gsCore/gsDebug.h>
#include <gsCore/gsConfig.h>
#include <gsMatrix/gsMatrix.h>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsMatrix/gsVector.h>
#include <gsMatrix/gsSparseVector.h>


#pragma once

namespace gismo {

template<typename T>
void initStream(std::ostream &stream)
{
    stream.precision(REAL_DIG);
    stream<< std::scientific;
}

template <typename T>
void openFile(const std::string &filename, std::ofstream &file)
{
    file.open(filename.c_str());
    initStream<T>(file);
}

template <typename T>
void openFile(const std::string &filename, std::ifstream &file)
{
    file.open(filename.c_str());
}




template <typename MatrixT>
struct gsMatrixExportImplementation
{
    static void matrix_export(std::ostream &stream, const MatrixT &matrix)
    {
        GISMO_ERROR("Cannot export matrix of type "<< std::string(typeid(MatrixT).name()) <<",add your own template specialization, please.");
    }
};

template <class T, int _Rows, int _Cols, int _Options>
struct gsMatrixExportImplementation<Eigen::Matrix<T,_Rows,_Cols,_Options> >
{
    typedef Eigen::Matrix<T,_Rows,_Cols,_Options> MatrixT;
    static void matrix_export(std::ostream &stream, const MatrixT &matrix)
    {

        for (index_t r=0; r<matrix.rows(); ++r)
        {
            index_t c=0;
            for ( ; c<matrix.cols()-1; ++c)
            {
                stream<<matrix(r,c)<<"\t";
            }
            stream<<matrix(r,c)<<"\n";
        }
    }
};

template <typename T, int _Options, typename _Index>
struct gsMatrixExportImplementation<Eigen::SparseMatrix<T,_Options,_Index> >
{
    typedef Eigen::SparseMatrix<T,_Options,_Index> MatrixT;
    static void matrix_export(std::ostream &stream, const MatrixT &matrix)
    {
        // specify matrix size
        stream <<(matrix.rows())<<" "<<(matrix.cols())<<" 0\n";
        for (int k=0; k < matrix.outerSize(); ++k)
        {
            for (typename MatrixT::InnerIterator it(matrix,k); it; ++it)
            {
                if( it.value() != 0.0 )
                    stream <<(it.row()+1)<<" "<<(it.col()+1)<<" "<<it.value()<<"\n";
            }
        }
    }
};


/**
 * \brief write a matrix in ascii format to a file
 * the format is suitable to use for exchanging data with external programs
 * in particular the following matlab code is guaranteed to work:
 * -for full matrices
 *   ...
 *   M = dlmread('filename')
 *   ...
 * -for sparse matrices
 *   ...
 *   M = spconvert(load('filename'))
 *   ...
 */
template <typename MatrixT>
void exportMatrixToASCII (const std::string &filename, const MatrixT &matrix)
{
    typedef typename MatrixT::Scalar T;
    std::ofstream stream;
    openFile<T>(filename, stream);
    gsMatrixExportImplementation<typename MatrixT::Base>::matrix_export(stream,matrix);
    stream.close();
}

/**
 * \brief as the above, but writing the data to any ostream: a pipe for example.
 */
template <typename MatrixT>
void exportMatrixToASCII (std::ostream &stream, const MatrixT &matrix)
{
    typedef typename MatrixT::Scalar T;
    initStream<T>(stream); // should we remove it?
    gsMatrixExportImplementation<typename MatrixT::Base>::matrix_export(stream,matrix);
}

// template<typename Derived>
// void exportMatrixToASCII<Eigen::DenseBase<Derived> >(const std::string &filename, const MatrixT &matrix) { }

// template<typename Derived>
// void exportMatrixToASCII<Eigen::SparseBase<Derived> >(const std::string &filename, const MatrixT &matrix) { }

template <typename MatrixT>
struct gsMatrixImportImplementation
{
    static void matrix_import(std::istream &stream, const MatrixT &matrix)
    {
        GISMO_ERROR("Cannot import matrix of type "<< std::string(typeid(MatrixT).name()) <<", add your own template specialization, please.");
    }
};

template <class T, int _Rows, int _Cols, int _Options>
struct gsMatrixImportImplementation<Eigen::Matrix<T,_Rows,_Cols,_Options> >
{
    typedef Eigen::Matrix<T,_Rows,_Cols,_Options> MatrixT;
    typedef typename MatrixT::Scalar ST;
    typedef typename MatrixT::Index  IT;

    static void readRow(std::istream &stream,std::vector<T> &data, IT &numEntries)
    {
        ST val;
        numEntries=0;
        while (true) // read a column
        {
            if (stream.peek() == '\r' || stream.peek() == '\n')
            {
                stream.get();
                break;
            }
            if ( stream >> val )
            {
                data.push_back(val);
                ++numEntries;
            }
            if( !stream.good() )
                break;
        }
    }

    static void matrix_import(std::istream &stream, MatrixT &matrix)
    {
        std::vector<T> data;
        IT rows = 0;
        IT curCol = 0;
        IT cols = 0;

        while ( cols==0 && !stream.eof() ) // read up to the first non empty row
            readRow(stream,data,cols);
        if (cols == 0)
        {
            gsWarn << "Empty file.";
            matrix.resize(0,0);
            return;
        }
        ++rows;
        while (true) // read each remaining row
        {
            readRow(stream,data,curCol);
            if (cols==curCol)
                ++rows;
            else if(curCol>0)
            {
                gsWarn << "Your rows should be of the same size but, unfortunately, they aren't.";
                matrix.resize(0,0);
                return;
            }
            else
                break;
        }
        matrix=Eigen::Map<Eigen::Matrix<ST,-1,-1,Eigen::RowMajor> >(data.data(),rows,cols);
    }
};

template <typename T, int _Options, typename _Index>
struct gsMatrixImportImplementation<Eigen::SparseMatrix<T,_Options,_Index> >
{
    typedef Eigen::SparseMatrix<T,_Options,_Index> MatrixT;
    typedef Eigen::Triplet<T> tripletT;
    static void matrix_import(std::istream &stream, MatrixT &matrix)
    {
        std::vector<tripletT> triplets;
        //std::stringstream line;
        _Index row, col;
        T val;

        stream >> row;
        stream >> col;
        stream >> val; // the extra 0;
        if (!stream.good())
        {
            gsWarn<<"Empty stream gives empty matrix! Cannot tell you more details.\n";
            matrix.resize(0, 0);
            return;
        }

        matrix.resize(row, col);

        while( true )
        {
            stream >> row;
            stream >> col;
            stream >> val;
            if (!stream.good())
                break;
            //std::cout << row-1 << " " << col-1 << " " << val << std::endl;
            triplets.push_back(tripletT(row-1, col-1, val));
        }
        matrix.setFromTriplets(triplets.begin(), triplets.end());
        matrix.makeCompressed();
    }
};

template <typename MatrixT>
void importMatrixFromASCII (const std::string &filename, MatrixT &matrix)
{
    typedef typename MatrixT::Scalar T;
    std::ifstream stream;
    openFile<T>(filename, stream);
    gsMatrixImportImplementation<typename MatrixT::Base>::matrix_import(stream,matrix);
    stream.close();
}

template <typename MatrixT>
void importMatrixFromASCII (std::istream &stream, MatrixT &matrix)
{
    gsMatrixImportImplementation<typename MatrixT::Base>::matrix_import(stream,matrix);
}

} // gismo
