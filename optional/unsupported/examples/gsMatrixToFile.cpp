/** @file gsMatrixToFile.cpp

    @brief Test for gsIO/gsMatrixToFile.h.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gismo.h>
#include <gsIO/gsMatrixToFile.h>
#include <ctime>

using namespace gismo;

int main( int argc, char** argv)
{
    std::string fn;
    std::string c = "rs";
    index_t rows = 5;
    index_t cols = 5;
    
    gsCmdLine cmd("Read/Write (sparse) matrices from/to file.");
    cmd.addString("f", "filename", "Filename",                                                                  fn  );
    cmd.addString("x", "command",  "Command to be exeuted: (rd,rs,wd,ws); r=read, w=write, s=sparse, d=dense.", c   );
    cmd.addInt   ("r", "rows",     "Number of rows (only for writing)",                                         rows);
    cmd.addInt   ("c", "cols",     "Number of columns (only for writing)",                                      cols);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (fn.empty())
    {
        gsInfo << "No filename provided. Invoke -h for help.\n";
        return 0;
    }

    if (rows < 1 || cols < 1)
    {
        gsInfo << "Rows and cols need to be > 0. Invoke -h for help.\n";
        return 1;
    }

    std::srand(std::time(0));
    
    if (c=="rd")
    {
        gsMatrix<> m;
        readFromBinaryFile(fn, m);
        gsInfo << "There is a " << m.rows() << "x" << m.cols() << " matrix:\n\n";
        if (m.rows() < 100 && m.cols() < 100)
            gsInfo << m << "\n\n";
        else
            gsInfo << "To large to display.\n\n";

    }
    else if (c=="rs")
    {
        gsSparseMatrix<> m;
        readFromBinaryFile(fn, m);
        gsInfo << "There is a " << m.rows() << "x" << m.cols() << " matrix:\n\n";
        if (m.rows() < 100 && m.cols() < 100)
            gsInfo << m.toDense() << "\n\n";
        else
            gsInfo << "To large to display.\n\n";
    }
    else if (c=="wd")
    {
        gsMatrix<> m(rows,cols);
        for (index_t i=0; i<rows; ++i)
            for (index_t j=0; j<cols; ++j)
                m(i,j) = rand() / ( RAND_MAX / 20.) - 10.;
        
        gsInfo << "I made up a " << m.rows() << "x" << m.cols() << " dense matrix:\n\n";
        if (m.rows() < 100 && m.cols() < 100)
            gsInfo << m << "\n\n";        
        else
            gsInfo << "To large to display.\n\n";
        
        writeToBinaryFile(fn, m);
        gsInfo << "Stored the matrix to the specified file.\n";
        
    }
    else if (c=="ws")
    {
        index_t number = rows * 4;
        gsSparseEntries<> se;
        for (index_t i=0; i<number; ++i)
                se.add(
                    rand() / ( RAND_MAX / rows),
                    rand() / ( RAND_MAX / cols),
                    rand() / ( RAND_MAX / 20.) - 10.
                );
        
        gsSparseMatrix<> m(rows,cols);
        m.setFromTriplets(se.begin(), se.end());
        m.makeCompressed();
        
        gsInfo << "I made up a " << m.rows() << "x" << m.cols() << " sparse matrix:\n\n";
        if (m.rows() < 100 && m.cols() < 100)
            gsInfo << m.toDense() << "\n\n";        
        else
            gsInfo << "To large to display.\n\n";
        
        writeToBinaryFile(fn, m);
        gsInfo << "Stored the matrix to the specified file.\n";
        
    }
    else
    {
        gsInfo << "Unknown command. Invoke -h for help.\n";
        return 1;
    }
    return 0;
    
}


