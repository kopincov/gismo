/** @file gsMatrixToFile.h

    @brief Allows to write sparse and dense matrices to files and to read them.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/

#include <gsMatrix/gsSparseMatrix.h>
#include <gsMatrix/gsMatrix.h>
#include <gsIO/gsFileData.h>
#include <gsIO/gsFileManager.h>

namespace gismo {

// https://scicomp.stackexchange.com/questions/21417/eigen-store-sparse-matrix-as-binary

// ti and td are just known numbers such that we can compare if integers and doubles are stored properly

template <typename T, int _Options, typename _Index>
bool writeToBinaryFile(const std::string& fn, Eigen::SparseMatrix<T, _Options, _Index>& m) {
    if(fn.empty())
    {
        gsWarn << "writeToBinaryFile: Filename is empty.\n";
        return false;
    }

    m.makeCompressed();

    std::fstream f;
    f.open(fn.c_str(), std::ios::binary | std::ios::out);

    if(!f.is_open())
    {
        gsWarn << "writeToBinaryFile: Cannot open file \"" << fn << "\".\n";
        return false;
    }
    T td;
    _Index ti, rows, cols, nonZeros, outerSize, innerSize;
    ti        = 1984;
    td        = (T)-3.14159;
    rows      = m.rows();
    cols      = m.cols();
    nonZeros  = m.nonZeros();
    outerSize = m.outerSize();
    innerSize = m.innerSize();

    f.write((const char*)&(ti),        sizeof(_Index));
    f.write((const char*)&(td),        sizeof(T     ));
    f.write((const char*)&(rows),      sizeof(_Index));
    f.write((const char*)&(cols),      sizeof(_Index));
    f.write((const char*)&(nonZeros),  sizeof(_Index));
    f.write((const char*)&(outerSize), sizeof(_Index));
    f.write((const char*)&(innerSize), sizeof(_Index));

    f.write((const char*)(m.valuePtr()),      sizeof(T)      * m.nonZeros() );
    f.write((const char*)(m.outerIndexPtr()), sizeof(_Index) * m.outerSize());
    f.write((const char*)(m.innerIndexPtr()), sizeof(_Index) * m.nonZeros() );

    f.close();
    return true;
}

template <typename T, int _Options, typename _Index>
bool readFromBinaryFile(const std::string& fn, Eigen::SparseMatrix<T, _Options, _Index>& m) {

    std::string fn2 = gsFileManager::find(fn);

    if(fn2.empty())
    {
        gsWarn << "readFromBinaryFile: Cannot find file \"" << fn << "\".\n";
        return false;
    }

    std::fstream f;
    f.open(fn2.c_str(), std::ios::binary | std::ios::in);

    if(!f.is_open())
    {
        gsWarn << "readFromBinaryFile: Cannot open file \"" << fn2 << "\".\n";
        return false;
    }
    T td;
    _Index ti, rows, cols, nnz, inSz, outSz;
    f.read((char*)&ti,    sizeof(_Index));
    f.read((char*)&td,    sizeof(T     ));
    f.read((char*)&rows,  sizeof(_Index));
    f.read((char*)&cols,  sizeof(_Index));
    f.read((char*)&nnz,   sizeof(_Index));
    f.read((char*)&inSz,  sizeof(_Index));
    f.read((char*)&outSz, sizeof(_Index));

    GISMO_ENSURE ( ti == 1984,        "That's not a sparse matrix." );
    GISMO_ENSURE ( td == (T)-3.14159, "That's not a sparse matrix." );

    m.resize(rows, cols);
    m.makeCompressed();
    m.resizeNonZeros(nnz);

    f.read((char*)(m.valuePtr()),      sizeof(T) * nnz       );
    f.read((char*)(m.outerIndexPtr()), sizeof(_Index) * outSz);
    f.read((char*)(m.innerIndexPtr()), sizeof(_Index) * nnz  );

    m.finalize();
    f.close();
    return true;
}

template <typename T, int _Rows, int _Cols>
bool writeToBinaryFile(const std::string& fn, Eigen::Matrix<T, _Rows, _Cols>& m) {
    typedef int _Index;
    if(fn.empty())
    {
        gsWarn << "writeToBinaryFile: Filename is empty.\n";
        return false;
    }

    std::fstream f;
    f.open(fn.c_str(), std::ios::binary | std::ios::out);

    if(!f.is_open())
    {
        gsWarn << "writeToBinaryFile: Cannot open file \"" << fn << "\".\n";
        return false;
    }
    T td;
    _Index ti, rows, cols;
    ti   = 2012;
    td   = (T)-2.71828;
    rows = m.rows();
    cols = m.cols();

    f.write((const char*)&(ti),   sizeof(_Index));
    f.write((const char*)&(td),   sizeof(T     ));
    f.write((const char*)&(rows), sizeof(_Index));
    f.write((const char*)&(cols), sizeof(_Index));

    f.write((const char*)(m.data()), sizeof(T) * rows * cols);

    f.close();
    return true;
}

template <typename T, int _Rows, int _Cols>
bool readFromBinaryFile(const std::string& fn, Eigen::Matrix<T, _Rows, _Cols>& m) {
    typedef int _Index;

    std::string fn2 = gsFileManager::find(fn);

    if(fn2.empty())
    {
        gsWarn << "readFromBinaryFile: Cannot find file \"" << fn << "\".\n";
        return false;
    }

    std::fstream f;
    f.open(fn2.c_str(), std::ios::binary | std::ios::in);

    if(!f.is_open())
    {
        gsWarn << "readFromBinaryFile: Cannot open file \"" << fn2 << "\".\n";
        return false;
    }
    T td;
    _Index ti, rows, cols;
    f.read((char*)&ti,   sizeof(_Index));
    f.read((char*)&td,   sizeof(T     ));
    f.read((char*)&rows, sizeof(_Index));
    f.read((char*)&cols, sizeof(_Index));

    GISMO_ENSURE( ti == 2012,        "That's not a dense matrix." );
    GISMO_ENSURE( td == (T)-2.71828, "That's not a dense matrix." );

    m.resize(rows, cols);

    f.read((char*)(m.data()), sizeof(T) * rows * cols);

    f.close();
    return true;
}

}
