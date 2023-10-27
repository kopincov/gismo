#pragma once

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>

#include <gsSmoothPatches/gsCompositeIncrSmoothnessBasis.h>
#include <gsSmoothPatches/gsCompositeHBasis.h>
#include <gsSmoothPatches/gsCompositeBSplineBasis.h>

#include <gsIO/gsFileData.h>
#include <gsIO/gsReadFile.h>

namespace gismo
{

template<short_t d, class T>
gsCompositeIncrSmoothnessBasis<d,T> * getCompBasisFromMultiPatch(const gsMultiPatch<> & mp,int incrSmoothness = -1,int minEVDistance = -1 )
{
    gsCompositeIncrSmoothnessBasis<d,T> * compBasis=NULL;
    bool tensorBSpline= true;
    bool hTensor = true;
    std::vector<gsTensorBSplineBasis<d,T>* >tensorBases;
    for(size_t i = 0;i<mp.nPatches();++i)
    {
        tensorBases.push_back(dynamic_cast<gsTensorBSplineBasis<d,T> * >(& mp.basis(i)));
        tensorBSpline = tensorBSpline && tensorBases[i]!=NULL;
    }
    if(tensorBSpline)
        compBasis = (new gsCompositeBSplineBasis<d,T>(tensorBases,mp,incrSmoothness,minEVDistance));
    else
    {
        std::vector<gsHTensorBasis<d,T>* >hBases;
        for(size_t i = 0;i<mp.nPatches();++i)
        {
            hBases.push_back(dynamic_cast<gsHTensorBasis<d,T> * >(& mp.basis(i)));
            hTensor = hTensor && hBases[i]!=NULL;
        }
        if(hTensor)
            compBasis = (new gsCompositeHBasis<d,T>(hBases,mp,incrSmoothness,minEVDistance));
    }
    GISMO_ASSERT(tensorBSpline||hTensor,"No suitable basis for gsMappedBasis found.");
    return compBasis;
}

template<short_t d, class T>
gsCompositeIncrSmoothnessBasis<d,T> * getCompBasisFromMultiPatch_withCoefs(const gsMultiPatch<> & mp, std::vector<gsMatrix<T>* >&coefs,int incrSmoothness = -1,int minEVDistance = -1 )
{
    gsCompositeIncrSmoothnessBasis<d,T> * compBasis=NULL;
    bool tensorBSpline= true;
    bool hTensor = true;
    std::vector<gsTensorBSplineBasis<d,T>* >tensorBases;
    for(size_t i = 0;i<mp.nPatches();++i)
    {
        tensorBases.push_back(dynamic_cast<gsTensorBSplineBasis<d,T> * >(& mp.basis(i)));
        tensorBSpline = tensorBSpline && tensorBases[i]!=NULL;
    }
    if(tensorBSpline)
        compBasis = (new gsCompositeBSplineBasis<d,T>(tensorBases,mp,coefs,incrSmoothness,minEVDistance));
    else
    {
        std::vector<gsHTensorBasis<d,T>* >hBases;
        for(size_t i = 0;i<mp.nPatches();++i)
        {
            hBases.push_back(dynamic_cast<gsHTensorBasis<d,T> * >(& mp.basis(i)));
            hTensor = hTensor && hBases[i]!=NULL;
        }
        if(hTensor)
            compBasis = (new gsCompositeHBasis<d,T>(hBases,mp,coefs,incrSmoothness,minEVDistance));
    }
    GISMO_ASSERT(tensorBSpline||hTensor,"No suitable basis for gsMappedBasis found.");
    return compBasis;
}

} // end namespace gismo
