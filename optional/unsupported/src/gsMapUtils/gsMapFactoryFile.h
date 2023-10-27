/** @file gsMapFactoryFile.h

    @brief Constructs a gsWeightMapper from a gsSparseMatrix read from a File.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsMapUtils/gsMapFactory.h>

namespace gismo
{

class GISMO_EXPORT gsMapFactoryFile : public gsMapFactory
{
public:

    /// Constructor using a path to input file (full path is expected)
    gsMapFactoryFile(std::string pathToInputFile) :
        m_pathToInputFile(pathToInputFile){}

    virtual gsWeightMapper<real_t> * makeMapper() const;

private:
    const std::string m_pathToInputFile;
};

}
