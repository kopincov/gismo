/** @file gsRemappedBasis.cpp

    @brief Tests the ambitious gsRemappedBasis project.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
**/

#include "gismo_unittest.h"
#include <gsRemappedBasis/gsTensorBasesUtils.h>
#include <gsRemappedBasis/gsRemappedBasis.h>

SUITE(gsRemappedBasis)
{
    template <short_t dim>
    struct RemappedTensorProduct
    {
        typedef typename gsRemTypes<dim>::tensorBasisT basisT;
        basisT                    *tensor;
        gsRemappedBasis           *remapped;
        gsDomainIterator<real_t>  *domIt;

        RemappedTensorProduct()
        {
            // basis
            gsKnotVector<> k(0,1,3,4,2);
            std::vector<gsKnotVector<> > knots(dim,k);
            basisT temp(knots);
            tensor = new basisT(temp);

            // domain iterator
            domIt = tensor->makeDomainIterator().release();

            // mapper init
            gsWeightMapper<real_t> mapper(tensor->size(),tensor->size());
            mapper.asMatrix().setIdentity();

            // selector init
            gsSelector selector;
            gsBoxList          boxes(dim);
            gsMatrix<real_t>   bbox=tensor->support();
            boxes.append(bbox,0);
            selector.initFromBoxesMax(boxes,bbox);

            // vector of basis
            std::vector<gsFunctionSet<real_t>*> basis(dim,tensor);

            // remapped
            remapped = new gsRemappedBasis(mapper, selector ,basis);
        }

        ~RemappedTensorProduct()
        {
            delete tensor;
            delete remapped;
            delete domIt;
        }


        void testEvalSameElement()
        {
            gsMatrix<real_t> bbox(dim,2);
            gsMatrix<real_t> points;
            gsFuncData<real_t> tenData(NEED_ACTIVE | NEED_VALUE | NEED_DERIV | NEED_DERIV2 | SAME_ELEMENT);
            gsFuncData<real_t> remData(NEED_ACTIVE | NEED_VALUE | NEED_DERIV | NEED_DERIV2 | SAME_ELEMENT);

            gsDomainIterator<real_t> &elem = *domIt;
            for (elem.reset(); elem.good(); elem.next())
            {
                bbox.col(0) =  elem.lowerCorner();
                bbox.col(1).array() = 0.02*elem.lowerCorner().array()+0.98*elem.upperCorner().array(); // do not stay on the upper boundary

                points = gsPointGrid(bbox, 16);
                static_cast<gsFunctionSet<real_t>*>(tensor) ->compute(points, tenData);
                remapped ->compute(points, remData);

                CHECK( tenData.actives==remData.actives);
                CHECK( tenData.values[0]==remData.values[0]);
                CHECK( tenData.values[1]==remData.values[1]);
                CHECK( tenData.values[2]==remData.values[2]);
            }
        }

        void testEvalDifferentElements()
        {
            gsMatrix<real_t> bbox(dim,2);
            bbox.col(0).setZero();
            bbox.col(1).setConstant(1);

            gsMatrix<real_t> points = gsPointGrid(bbox, 300);
            gsFuncData<> tenData(NEED_ACTIVE | NEED_VALUE | NEED_DERIV | NEED_DERIV2);
            gsFuncData<> remData(NEED_ACTIVE | NEED_VALUE | NEED_DERIV | NEED_DERIV2);

            static_cast<gsFunctionSet<real_t>*>(tensor) ->compute(points, tenData);
            remapped ->compute(points, remData);

            CHECK( tenData.actives==remData.actives);
            CHECK( tenData.values[0]==remData.values[0]);
            CHECK( tenData.values[1]==remData.values[1]);
            CHECK( tenData.values[2]==remData.values[2]);
        }

    };


    typedef RemappedTensorProduct<1> RemappedTensorProduct1D;
    TEST_FIXTURE(RemappedTensorProduct1D, evalSameElement1D)
    { testEvalSameElement(); }
    TEST_FIXTURE(RemappedTensorProduct1D, evalManyElement1D)
    { testEvalDifferentElements(); }

    typedef RemappedTensorProduct<2> RemappedTensorProduct2D;
    TEST_FIXTURE(RemappedTensorProduct2D, evalSameElement2D)
    { testEvalSameElement(); }
    TEST_FIXTURE(RemappedTensorProduct2D, evalManyElement2D)
    { testEvalDifferentElements(); }

    typedef RemappedTensorProduct<3> RemappedTensorProduct3D;
    TEST_FIXTURE(RemappedTensorProduct3D, evalSameElement3D)
    { testEvalSameElement(); }
    TEST_FIXTURE(RemappedTensorProduct3D, evalManyElement3D)
    { testEvalDifferentElements(); }




    template <short_t dim>
    struct   Remapped2Basis
    {
        gsFunctionSet<real_t>     *tensor[2];
        gsRemappedBasis           *remapped;
        gsMatrix<real_t>           box1;
        gsMatrix<real_t>           box0;

        Remapped2Basis()
        {
            typedef gsTensorBSplineBasis<dim, real_t> basisT;
            gsKnotVector<> k(0,1,3,4,2);
            std::vector<gsKnotVector<> > knots(dim,k);
            tensor[0] = new basisT(knots);
            basisT *temp = new basisT(knots);
            temp ->degreeElevate();
            temp ->uniformRefine();
            tensor[1] = temp;

            // mapper
            const index_t totalSize=tensor[0]->size()+tensor[1]->size();
            gsWeightMapper<real_t> mapper(totalSize,totalSize);
            mapper.asMatrix().setIdentity();

            // selector
            gsSelector selector;
            gsBoxList          boxes(dim);
            gsMatrix<real_t>   bbox=static_cast<basisT*>(tensor[0])->support();
            box1.resizeLike(bbox);
            box0.resizeLike(bbox);

            box1.col(0)=(bbox.col(0)+bbox.col(1)).array()/2;
            box1.col(1)=bbox.col(1);
            box0.col(0)=bbox.col(0);
            box0.col(1)=0.99*box1.col(0).array(); // points in the top boundary are in level 1

            boxes.append(bbox,0);
            boxes.append(box1,1);
            selector.initFromBoxesMax(boxes,bbox);

            std::vector<gsFunctionSet<real_t>*> vec(2);
            vec[0]=tensor[0];
            vec[1]=tensor[1];
            remapped = new gsRemappedBasis(mapper, selector ,vec);
        }
        ~Remapped2Basis()
        {
            freeAll(tensor,tensor+2);
            delete remapped;
        }
    };

    typedef Remapped2Basis<2> Remapped2Basis2;
    TEST_FIXTURE( Remapped2Basis2 , basis2D)
    {
            gsFuncData<real_t> tenData(NEED_ACTIVE | NEED_VALUE | NEED_DERIV | NEED_DERIV2);
            gsFuncData<real_t> remData(NEED_ACTIVE | NEED_VALUE | NEED_DERIV | NEED_DERIV2);

            gsMatrix<real_t> points = gsPointGrid(box0, 300);
            tensor[0]->compute(points, tenData);
            remapped ->compute(points, remData);

            CHECK( tenData.actives==remData.actives);
            CHECK( tenData.values[0]==remData.values[0]);
            CHECK( tenData.values[1]==remData.values[1]);
            CHECK( tenData.values[2]==remData.values[2]);

            points = gsPointGrid(box1, 300);
            tensor[1]->compute(points, tenData);
            remapped ->compute(points, remData);

            CHECK( (tenData.actives.array()+tensor[0]->size()).matrix()==remData.actives );
            CHECK( tenData.values[0]==remData.values[0]);
            CHECK( tenData.values[1]==remData.values[1]);
            CHECK( tenData.values[2]==remData.values[2]);
    }

    // This was testing whether the export to tex gives plausible results.
    /*TEST_FIXTURE( Remapped2Basis2, exportMapToTex )
    {
        gsMatrix<real_t> bb;
        bb << 0, 1, 0, 1;
        remapped->exportMapToTex("mesh", bb);
    }*/

}


