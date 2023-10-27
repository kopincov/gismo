/** @file gsFittingConstr.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include<gsFitting/gsFittingUtilsCheck.h>

namespace gismo
{

template<class T>
bool positiveDet(gsBasis<T>& basis,
                 gsFunctionSet<T>& func, T tol)
{
    gsGeometry<T>& map = getGeometry(basis, func);
    return positiveDet(map);
}

template<class T>
bool positiveDet(gsMultiBasis<T>& basis,
                 gsFunctionSet<T>& func, T tol)
{
    gsMultiPatch<T>& map = getGeometry(basis, func);
    return positiveDet(map);
}

template<class T>
bool positiveDet(gsMultiPatch<T>& func, T tol)
{
    unsigned nPatches = func.nPatches();
    for(unsigned i = 0;i < nPatches;i++)
    {
        if(! positiveDet(func.patch(i)) )
            return false;
    }
    return true;
}


template<class T>
bool positiveDet(gsGeometry<T>& func, T tol)
{
    int dim = func.parDim();
    GISMO_ASSERT(dim == 2 || dim == 3,
                 "Only implemented for dimension 2 and 3");
    int dim2 = dim * dim;
    gsMatrix<T> bb = func.basis().support();

    int nbpts = math::pow(2, dim);
    gsMatrix<T> del(dim, 1);
    gsMatrix<T> pts(dim, nbpts);
    gsMatrix<T> der(dim2, nbpts);
    T deter = 0.;
    T sgn = 1.;
    T eps = 1e-6;

    for(int i = 0;i < dim+1;i++)
        pts.col(i) = bb.col(0);
    pts.col(nbpts - 1) = bb.col(1);

    del = bb.col(1) - bb.col(0);
    for(int i = 0;i < dim;i++)
    {
        pts(i,i+1) += del(i,0) * (1 - eps);

        if(dim == 3)
        {
            pts.col(i+dim+1) = pts.col(i+1);
            int j = (i+1) % 3;
            pts(j, i+dim+1) += del(j,0) * (1 - eps);
        }
    }
    func.deriv_into(pts, der);
    for(int i = 0;i < nbpts;i++)
    {
        gsMatrix<T> mat(dim, dim);
        for(int j = 0;j < dim;j++)
        {
            for(int k = 0;k < dim;k++)
                mat(j, k) = der(dim*j+k, i);
        }
        deter = mat.determinant();
        if(i == 0)
        {
            if(deter > 0.)
                sgn = 1.;
            else
                sgn = -1.;
        }
        if(deter * sgn < tol)
            return false;
    }
    if(sgn < 0.)
        gsWarn << "The determinant is negative everywhere. This should not be a problem for optimization if the Winslow energy is not considered." << std::endl;
    return true;
}


} /// namespace gismo
