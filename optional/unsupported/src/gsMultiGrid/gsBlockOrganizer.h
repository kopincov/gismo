/** @file gsBlockOrgaizer.h

@brief Organizes blocks for linear system of equations for block smoothers

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): J. Sogn
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsBasis.h>

namespace gismo
{

/** @brief
Utility class for organizing blocks for linear system of equations for block smoothers.

The block information is a std::vector<gsVector<index_t> > where
each gsVector<index_t> is a block from the indices in gsVector<index_t>.
 \ingroup Solver
 */
class gsBlockOrgaizer
{
public:

    ///Find the DoFs close to the boundary (with offset) and creats one big block.
    ///The remaining DoFs are put in blocks of size 1 (point blocks).
    static std::vector<gsVector<index_t> > boundaryBlock(const gsBasis<real_t> &basis,const gsDofMapper &mapper, index_t offset)
    {
        int sides;
        if (basis.dim() == 2)
            sides = 4;
        else if(basis.dim() == 3)
            sides = 6;
        else
            GISMO_ERROR("Not implementet for this dimention!");
        std::vector<gsVector<index_t> > blockInfo;

        index_t dofs = mapper.freeSize();

        //Find dof close to boundary and add them to list
        gsMatrix<index_t> boundaryIndex;
        std::vector<index_t> closeToBoundary;
        for (int  side = 1; side <= sides; ++side)
        {
            for (int i = 0; i<= offset; ++i)
            {
                boundaryIndex = basis.boundaryOffset(side,i);
                for (int k = 0; k< boundaryIndex.rows(); ++k)
                {
                    //Index corresponding to basis.
                    index_t basisIndex = boundaryIndex(k,0);
                    if (mapper.is_free(basisIndex))
                    {
                        index_t globalIndex = mapper.index(basisIndex);//Index inside system matrix
                        //Add dof to list if not alreay added.
                        if ( std::find(closeToBoundary.begin(), closeToBoundary.end(), globalIndex) == closeToBoundary.end())
                            closeToBoundary.push_back(globalIndex);
                    }
                }
            }
        }
        //Sort the list (not needed)
        //std::sort (closeToBoundary.begin(), closeToBoundary.end());

        //convert from std::vector, to eigen vector
        index_t blockSize = closeToBoundary.size();
        gsVector<index_t> closeToBoundaryVector(blockSize);
        for (index_t k = 0; k< blockSize; k++)
            closeToBoundaryVector(k) = closeToBoundary[k];

        blockInfo.push_back(closeToBoundaryVector);

        for (index_t k = 0; k< dofs; k++)
        {
            //If the DoF is in the main block, add it in a single block.
            if ( (std::find(closeToBoundary.begin(), closeToBoundary.end(), k) == closeToBoundary.end()) )
            {
                 gsVector<index_t> tempVector(1);
                 tempVector(0) = k;
                 blockInfo.push_back(tempVector);
            }
        }

        GISMO_ASSERT(checkConsistency(blockInfo,dofs),"Consistency test failed!");
        return blockInfo;
    }

    /// Creats blocks of size \a blockSize in chronological order.
    /// The last block has differnt size if the modulus is not zero.
    /// if blockSize is 1, we get normal point smoother.
    static std::vector<gsVector<index_t> > uniformBlocks(index_t blockSize, index_t size)
    {
        std::vector<gsVector<index_t> > blockInfo;
        const index_t numberOfBlocks = size/blockSize;
        const index_t modulus = size%blockSize;
        //const index_t sizeSquare = math::sqrt(size);
        for (index_t k = 0; k< numberOfBlocks; k++)
        {
            //index_t tmp = k/sizeSquare;
            gsVector<index_t> block(blockSize);
            for (index_t j = 0; j< blockSize; j++)
            {
                block(j) = blockSize*k + j;
                //block(j) = k+ tmp*sizeSquare + j*sizeSquare;//conecting dof in y direction: Used for testing
                //block(j) = blockSize*k-(k%2)+ j*2; //alternating DoF 02, 13, 46, 57, 8 10, 9 11: Used for testing
            }
            blockInfo.push_back(block);
        }
        //Fill the last block
        if (modulus > 0)
        {
            gsVector<index_t> block(modulus);
            for (index_t j = 0; j< modulus; j++)
                block(j) = blockSize*numberOfBlocks + j;

            blockInfo.push_back(block);
        }

        GISMO_ASSERT(checkConsistency(blockInfo,size),"Consistency test failed!");
        return blockInfo;
    }

    /// Checks if all DOFs are included and only appear once
    static bool checkConsistency(std::vector<gsVector<index_t> > & blockInfo, index_t size)
    {
        const index_t numberOfBlocks = blockInfo.size();
        gsVector<index_t> multiplicity;
        multiplicity.setZero(size);

        for (index_t k = 0; k< numberOfBlocks; k++)
            for (index_t j = 0; j < blockInfo[k].rows(); j++)
            {
                if (blockInfo[k](j) < size)
                    multiplicity(blockInfo[k](j)) += 1;
                else
                    return false;
            }
        return multiplicity.maxCoeff() == 1 && multiplicity.minCoeff() == 1;
    }

    private:

    // Do not allow to create objects
    gsBlockOrgaizer()  {}
    ~gsBlockOrgaizer() {}

};
} // namespace gismo
