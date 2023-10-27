/** @file gsMapFactoryFile.cpp

    @brief test gsMapFactoryFile class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
**/

#include "gismo_unittest.h"
#include <gsMSplines/gsWeightMapper.h>
#include <gsMapUtils/gsMapFactoryFile.h>

SUITE(gsMapper_Factory)
{
    TEST(gsMapFactory_File)
    {
        std::string path("matrices/identity_sparse.xml");
        std::string fullPath = gsFileManager::find(path);
        gsMapFactoryFile mapFactory(fullPath);
        gsWeightMapper<real_t>* mapper = mapFactory.makeMapper();
        gsSparseMatrix<real_t> matrix;
        gsReadFile<>(fullPath, matrix);
        CHECK(gsAllCloseAbsolute(matrix.toDense(),mapper->asMatrix().toDense(),std::numeric_limits<real_t>::epsilon()));
        delete mapper;
    }
}
