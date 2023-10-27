/** @file gsRemappedTHB.cpp

    @brief Tests gsRemappedTHB, which is a part of the ambitious gsRemappedBasis project.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
**/

#include "gismo_unittest.h"
#include <gsRemappedBasis/gsTHB.h>
#include <gsUtils/gsExportMatrix.h>
#include <gsIO/gsIOUtils.h> // makeMesh

extern const char* data1;

SUITE(gsRemappedBasis)
{
    SUITE(gsRemappedTHB)
    {

        bool checkResultsValues(const gsFuncData<> &resultN, const gsFuncData<> &resultO)
        {
            const index_t derS=2;
            const index_t secS=3;
            const index_t maxDer=resultN.maxDeriv();

            bool ok=true;
            real_t tol = 100*std::numeric_limits<real_t>::epsilon();
            for (index_t p=0; p< resultO.values[0].cols()&& ok; ++p)
            {
                index_t oAI=0;
                index_t nAI=0;
                index_t aCol = (resultN.flags & SAME_ELEMENT) ? 0 : p;
                for (; oAI<resultO.actives.rows() && ok; ++oAI)
                {
                    if (resultN.actives(nAI,aCol) < resultO.actives(oAI,aCol))
                        ++nAI;
                    if (resultN.actives(nAI,aCol) > resultO.actives(oAI,aCol))
                    {
                        if (maxDer>=0) CHECK(resultO.values[0](oAI,p)==0);
                        if (maxDer>=1) CHECK(resultO.values[1].block(oAI*derS,p,derS,1).isZero());
                        if (maxDer>=2) CHECK(resultO.values[2].block(oAI*secS,p,secS,1).isZero());
                    }
                    else if (resultN.actives(nAI,aCol) == resultO.actives(oAI,aCol) )
                    {
                        real_t diffNormValue;
                        real_t diffNormDeriv;
                        real_t diffNormSec;

                        switch (maxDer)
                        {
                        default:
                        case 2:
                            diffNormSec   =  (resultO.values[2].block(oAI*secS,p,secS,1)-resultN.values[2].block(nAI*secS,p,secS,1)).norm()/math::max(resultO.values[2].block(oAI*secS,p,secS,1).norm(),(real_t)0.0001);
                            if ( diffNormSec >= tol )
                                ok = false;
                            CHECK(diffNormSec<=tol);
                        case 1:
                            diffNormDeriv = (resultO.values[1].block(oAI*derS,p,derS,1)-resultN.values[1].block(nAI*derS,p,derS,1)).norm()/math::max(resultO.values[1].block(oAI*derS,p,derS,1).norm(), (real_t)0.0001);
                            if ( diffNormDeriv >= tol )
                                ok = false;
                            CHECK(diffNormDeriv<=tol);
                        case 0:
                            diffNormValue = (resultO.values[0](oAI,p)-resultN.values[0](nAI,p)) / math::max(resultO.values[0].block(oAI,p,1,1).norm(), (real_t)0.0001);
                            diffNormValue = diffNormValue>0? diffNormValue : -diffNormValue;
                            if ( diffNormValue >= tol )
                                ok = false;
                            CHECK(diffNormValue<=tol);
                        case -1:
                            ;
                        }
                    }
                }
            }
            return ok;
        }


        void checkPerElement(gsFunctionSet<>* newB, gsFunctionSet<>* oldB, gsDomainIterator<> *it, int degree)
        {
            gsMatrix<> supp = oldB->support();

            gsVector<index_t> numQuad;
            numQuad.setConstant(newB->domainDim(),degree+1);
            gsGaussRule<> quad(numQuad);

            gsFuncData<> resultN(NEED_ACTIVE|NEED_VALUE|NEED_DERIV|NEED_DERIV2|SAME_ELEMENT);
            gsFuncData<> resultO(NEED_ACTIVE|NEED_VALUE|NEED_DERIV|NEED_DERIV2|SAME_ELEMENT);

            gsMatrix<> nodes;
            gsVector<> weights;
            for(;it->good();it->next())
            {
                quad.mapTo(it->lowerCorner(), it->upperCorner(),nodes, weights);
                newB->compute(nodes, resultN);
                oldB->compute(nodes, resultO);
                if (!checkResultsValues(resultN,resultO))
                    break;
            }
            delete it;
        }

        TEST(THB_comparison)
        {
            // typedefs
            typedef gsTensorBSplineBasis<2, real_t> BSbasisT;
            typedef gsHBSplineBasis<2, real_t>      HBbasisT;
            typedef gsTHBSplineBasis<2, real_t>     THbasisT;

            const int deg    = 2;
            const int maxLvl = 1;

            std::vector<real_t> knots(3);
            knots[0] = 0;
            knots[1] = 0.4;
            knots[2] = 1;
            gsKnotVector<real_t> KVx(knots, deg, 0); // knots, degree, regularity
            knots[1] = 0.6;
            gsKnotVector<real_t> KVy(knots, deg, 0);

            BSbasisT previousBasis(KVx, KVy);
            std::vector<gsBoxList::basisPtr> basisVector(1,gsBoxList::basisPtr(previousBasis.clone()));
            for (int l=0; l<maxLvl;++l)
            {
                previousBasis.uniformRefine();
                basisVector.push_back(gsBoxList::basisPtr(previousBasis.clone()));
            }

            // domain definition
            gsBoxList newBoxes(2);
            gsMatrix<real_t> box(2,2);
            box << 0.4, 1, 0.6, 1;
            newBoxes.append(box,1);

            std::vector<index_t> oldBoxes(5);
            oldBoxes[0]=1;
            oldBoxes[1]=deg+1;
            oldBoxes[2]=deg+1;
            oldBoxes[3]=deg+3;
            oldBoxes[4]=deg+3;

            // construct bases
            gsTHB<2>  newHB (basisVector, newBoxes,false);
            HBbasisT  oldHB (*dynamic_cast<gsBasis<real_t>*>(basisVector[0].get()));
            oldHB.refineElements(oldBoxes);

            gsTHB<2>  newTHB (basisVector, newBoxes,true);
            THbasisT  oldTHB (*dynamic_cast<gsBasis<real_t>*>(basisVector[0].get()));
            oldTHB.refineElements(oldBoxes);

            // define points
            const gsMatrix<real_t> points = gsPointGrid(previousBasis.support(),10);

            // storage for values
            gsFuncData<real_t> oldData(NEED_ACTIVE | NEED_VALUE | NEED_DERIV | NEED_DERIV2);
            gsFuncData<real_t> remData(NEED_ACTIVE | NEED_VALUE | NEED_DERIV | NEED_DERIV2);

            // compute HB and check
            oldHB.compute(points, oldData);
            newHB.compute(points, remData);

            CHECK( oldData.actives==remData.actives);
            for (int i=0; i<oldData.maxDeriv();++i)
                CHECK( gsAllCloseAbsolute(oldData.values[i],remData.values[i], 100*std::numeric_limits<real_t>::epsilon()) );

            // compute THB and check
            oldTHB.compute(points, oldData);
            newTHB.compute(points, remData);

            checkResultsValues(remData,oldData);
        }



        TEST( ReadAndCompare2d )
        {
            std::string files[] = {
                "basis2d/adapt_basis_0.xml",
                "basis2d/adapt_basis_1.xml",
                "basis2d/adapt_basis_2.xml"
            };
            size_t fileNum=sizeof(files)/sizeof(std::string);

            for(size_t fId=0; fId<fileNum;++fId)
            {
                gsFileData<> fileDataN(files[fId]);
                gsFileData<> fileDataO(files[fId]);
                gsTHB<2>::uPtr newB = fileDataN.getFirst< gsTHB<2> >();
                gsTHBSplineBasis<2>::uPtr oldB = fileDataO.getFirst< gsTHBSplineBasis<2> >();
                if (!oldB)
                {continue;}
                CHECK( (bool)newB );
                if (!newB)
                {continue;}
                CHECK(newB->size()==oldB->size());
//                {
//                    newB->exportSelectorToTex("newBmesh");
//                    if (system("pdflatex newBmesh.tex 2>/dev/null 1>&2 &"))
//                        std::cout<<"mesh not generated"<<std::endl;
//                    gsMesh<real_t> msh;
//                    makeMesh<real_t>(*oldB, msh,8);
//                    gsWriteParaview(msh, "oldBMesh", false);
//                }
                checkPerElement(newB.get(),oldB.get(),oldB->makeDomainIterator().release(),oldB->degree(0));
            }
        }


        TEST( ReadWriteAndCompare2d )
        {
            std::string files[] = {
                "basis2d/adapt_basis_0.xml",
                "basis2d/adapt_basis_1.xml",
                "basis2d/adapt_basis_2.xml"
            };
            size_t fileNum=sizeof(files)/sizeof(std::string);


            for(size_t fId=0; fId<fileNum;++fId)
            {
                gsFileData<> fileRead(files[fId]);
                gsFileData<> fileWrit;
                gsTHB<2>::uPtr readB = fileRead.getFirst< gsTHB<2> >();
                CHECK ( (bool)readB);
                fileWrit.add(*readB);
                fileWrit.dump("testBasisWrite");
                gsFileData<> fileRead2("testBasisWrite.xml");
                gsTHB<2>::uPtr writB = fileRead2.getFirst< gsTHB<2> >();
                CHECK ( (bool)writB);
                checkPerElement(readB.get(),writB.get(),writB->makeDomainIterator().release(), readB->degree());
            }
        }

        TEST( m_repr_Construction )
        {
            /* Computes a representation matrix and compares it with the known correct answer
             * for this particular mesh.
             */

            // Create boxes
            gsBoxList bList(2);
            gsMatrix<real_t> box(2,2);
            box << 0.000, 0.500, 0.000, 0.500;
            bList.append(box, 1);
            box << 0.500, 1.000, 0.000, 1.000;
            bList.append(box, 1);
            box << 0.250, 0.500, 0.000, 0.500;
            bList.append(box, 2);
            box << 0.500, 1.000, 0.000, 0.750;
            bList.append(box, 2);
            box << 0.500, 0.750, 0.375, 0.750;
            bList.append(box,3);

            /* The mesh should look like this now.
             * +---------------+-------+-------+
             * |               |       |       |
             * |               |       |       |
             * |               |       |       |
             * |               |   1   |   1   |
             * |               |       |       |
             * |               |       |       |
             * |               |       |       |
             * |       0       +-+-+-+-+---+---+
             * |               |3|3|3|3|   |   |
             * |               +-+-+-+-+ 2 | 2 |
             * |               |3|3|3|3|   |   |
             * |               +-+-+-+-+---+---+
             * |               |3|3|3|3|   |   |
             * |               +-+-+-+-+ 2 | 2 |
             * |               |3|3|3|3|   |   |
             * +-------+---+---+-+-+-+-+---+---+
             * |       |   |   |3|3|3|3|   |   |
             * |       | 2 | 2 +-+-+-+-+ 2 | 2 |
             * |       |   |   |3|3|3|3|   |   |
             * |   1   +---+---+-+-+-+-+---+---+
             * |       |   |   |   |   |   |   |
             * |       | 2 | 2 | 2 | 2 | 2 | 2 |
             * |       |   |   |   |   |   |   |
             * +-------+---+---+---+---+---+---+
             * |       |   |   |   |   |   |   |
             * |       | 2 | 2 | 2 | 2 | 2 | 2 |
             * |       |   |   |   |   |   |   |
             * |   1   +---+---+---+---+---+---+
             * |       |   |   |   |   |   |   |
             * |       | 2 | 2 | 2 | 2 | 2 | 2 |
             * |       |   |   |   |   |   |   |
             * +-------+---+---+---+---+---+---+
             */

            // Create gsTHB
            typedef gsKnotVector<real_t> gsCKV;
            typedef gsTensorBSplineBasis<2, real_t> basisT;

            std::vector<gsBoxList::basisPtr> basisVector;
            // Level 0
            gsCKV KV(0,1,1,1,4);
            basisVector.push_back(gsBoxList::basisPtr(new basisT(KV, KV)));
            // Level 1
            KV.uniformRefine();
            basisVector.push_back(gsBoxList::basisPtr(new basisT(KV, KV)));
            // Level 2
            KV.uniformRefine();
            basisVector.push_back(gsBoxList::basisPtr(new basisT(KV, KV)));
            // Level 3
            KV.uniformRefine();
            basisVector.push_back(gsBoxList::basisPtr(new basisT(KV, KV)));

            gsTHB<2> remappedBasis(basisVector, bList);


            // Read the saved correct answer.
            // Surely there is a more elegant way.
            gsSparseMatrix<real_t,Eigen::RowMajor,index_t> imported;
            std::stringstream stream;
            stream << data1;
            importMatrixFromASCII(stream, imported);

            // Compare your current representation matrix with the precomputed correct one.
            for( index_t i = 0; i < imported.rows(); ++ i )
                for( index_t j = 0; j < imported.cols(); ++j )
                    CHECK( imported(i,j) == remappedBasis.getMapper().asMatrix()(i,j));
        }
    }
}

const char* data1=
        "556 73 0\n"
        "21 1 1.000000000000000e+00\n"
        "26 2 1.000000000000000e+00\n"
        "33 3 1.000000000000000e+00\n"
        "73 4 1.000000000000000e+00\n"
        "74 5 1.000000000000000e+00\n"
        "77 6 1.000000000000000e+00\n"
        "78 7 1.000000000000000e+00\n"
        "82 8 1.000000000000000e+00\n"
        "83 9 1.000000000000000e+00\n"
        "84 10 1.000000000000000e+00\n"
        "85 11 1.000000000000000e+00\n"
        "88 12 1.000000000000000e+00\n"
        "89 13 1.000000000000000e+00\n"
        "93 14 1.000000000000000e+00\n"
        "94 15 1.000000000000000e+00\n"
        "95 16 1.000000000000000e+00\n"
        "96 17 1.000000000000000e+00\n"
        "99 18 1.000000000000000e+00\n"
        "100 19 1.000000000000000e+00\n"
        "104 20 1.000000000000000e+00\n"
        "105 21 1.000000000000000e+00\n"
        "106 22 1.000000000000000e+00\n"
        "107 23 1.000000000000000e+00\n"
        "110 24 1.000000000000000e+00\n"
        "111 25 1.000000000000000e+00\n"
        "117 26 1.000000000000000e+00\n"
        "118 27 1.000000000000000e+00\n"
        "128 28 1.000000000000000e+00\n"
        "129 29 1.000000000000000e+00\n"
        "139 30 1.000000000000000e+00\n"
        "140 31 1.000000000000000e+00\n"
        "150 32 1.000000000000000e+00\n"
        "151 33 1.000000000000000e+00\n"
        "161 34 1.000000000000000e+00\n"
        "162 35 1.000000000000000e+00\n"
        "172 36 1.000000000000000e+00\n"
        "173 37 1.000000000000000e+00\n"
        "321 38 1.000000000000000e+00\n"
        "322 39 1.000000000000000e+00\n"
        "323 40 1.000000000000000e+00\n"
        "324 41 1.000000000000000e+00\n"
        "340 42 1.000000000000000e+00\n"
        "341 43 1.000000000000000e+00\n"
        "342 44 1.000000000000000e+00\n"
        "343 45 1.000000000000000e+00\n"
        "359 46 1.000000000000000e+00\n"
        "360 47 1.000000000000000e+00\n"
        "361 48 1.000000000000000e+00\n"
        "362 49 1.000000000000000e+00\n"
        "378 50 1.000000000000000e+00\n"
        "379 51 1.000000000000000e+00\n"
        "380 52 1.000000000000000e+00\n"
        "381 53 1.000000000000000e+00\n"
        "397 54 1.000000000000000e+00\n"
        "398 55 1.000000000000000e+00\n"
        "399 56 1.000000000000000e+00\n"
        "400 57 1.000000000000000e+00\n"
        "416 58 1.000000000000000e+00\n"
        "417 59 1.000000000000000e+00\n"
        "418 60 1.000000000000000e+00\n"
        "419 61 1.000000000000000e+00\n"
        "435 62 1.000000000000000e+00\n"
        "436 63 1.000000000000000e+00\n"
        "437 64 1.000000000000000e+00\n"
        "438 65 1.000000000000000e+00\n"
        "454 66 1.000000000000000e+00\n"
        "455 67 1.000000000000000e+00\n"
        "456 68 1.000000000000000e+00\n"
        "457 69 1.000000000000000e+00\n"
        "473 70 1.000000000000000e+00\n"
        "474 71 1.000000000000000e+00\n"
        "475 72 1.000000000000000e+00\n"
        "476 73 1.000000000000000e+00\n";
