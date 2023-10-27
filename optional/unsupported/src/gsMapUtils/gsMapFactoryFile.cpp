/** @file gsMapFactoryFile.cpp

    @brief Provides implementation of gsMapFactoryFile class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsMapUtils/gsMapFactoryFile.h>

#include <gsIO/gsFileData.h>
#include <gsIO/gsReadFile.h>

namespace gismo
{

gsWeightMapper<real_t> * gsMapFactoryFile::makeMapper() const
{
    gsWeightMapper<real_t> * mapper;
    gsSparseMatrix<real_t>::uPtr mapMatrix;
    mapMatrix = gsReadFile<>(m_pathToInputFile);
    if(!mapMatrix)
        std::cout << "file could not be read" << std::endl;
    mapper = new gsWeightMapper<real_t>(*mapMatrix);
    return mapper;
}

}
