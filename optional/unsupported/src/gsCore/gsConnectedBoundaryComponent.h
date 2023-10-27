/** @file getConnectedBoundaryComponent.h

    @brief Provides test examples for multigrid algorithms

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/


namespace gismo {

std::vector<patchSide> getConnectedBoundaryComponent( const gsMultiPatch<> & mp, const patchSide& ps )
{
    const index_t np = mp.nPatches();
    const index_t dim = mp.dim();
    
    std::vector<patchSide> result;
    result.reserve(np * 2 * dim);
    result.push_back( ps );

    for ( size_t i = 0; i < result.size() /* can change */; ++i )
    {
        // get all neigbors form ps
        const patchSide & current = result[i];

        std::vector<patchCorner> corners;
        current.getContainedCorners( dim, corners );
        
        for ( size_t j = 0; j < corners.size(); ++j )
        {
            std::vector<patchCorner> pcout;
            mp.getCornerList( corners[j], pcout );                   // the same corners elsewhere
            for ( size_t k = 0; k < pcout.size(); ++k )
            {
                std::vector<patchSide> othersides;
                pcout[k].getContainingSides( dim, othersides );
                
                for ( size_t l = 0; l < othersides.size(); ++l )
                {
                    if ( std::find( result.begin(), result.end(), othersides[l] ) == result.end() )
                    {
                        if ( std::find( mp.bBegin(), mp.bEnd(), othersides[l] ) != mp.bEnd() )
                        {
                            result.push_back(othersides[l]);
                        }
                    }
                }
            }
        }
        
    }


    return result;
}

}
