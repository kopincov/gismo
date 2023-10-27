/** @file gsKnotSelector.cpp

    @brief Tests gsKnotSelector, which is a part of the ambitious
    gsRemappedBasis project.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
**/

#include "gismo_unittest.h"
#include <gsRemappedBasis/gsSelector.h>


typedef gsSelector::NodeId NodeId;

// Andrea, sorry, Dominik has commented out your tikz drawings to shorten the output.

const gsBoxList::basisIdT invalidBasisId = gsBoxList::basisIdT(UINT_MAX);

template <typename Data>
std::ostream& operator<< (std::ostream &out, const std::vector<Data> &v)
{
    for(size_t iter=0; iter<v.size(); ++iter)
    {
        out<<v[iter]<<", ";
    }
    return out;
}

class data
{
public:
    gsDomainMap domainMap;

    NodeId root;

    NodeId left;
    NodeId right;

    NodeId lowerLeft;
    NodeId upperLeft;
    NodeId lowerRight;
    NodeId upperRight;

    data()
    {
        root = domainMap.root();
        domainMap[root].data.dir = 0;
        domainMap[root].data.par = 0.4;

        left = domainMap.append(root,true);
        domainMap[left].data.dir = 1;
        domainMap[left].data.par = 0.3;

        right = domainMap.append(root,false);
        domainMap[right].data.dir = 1;
        domainMap[right].data.par = 0.7;

        lowerLeft  = domainMap.append(left,true);
        upperLeft  = domainMap.append(left,false);
        lowerRight = domainMap.append(right,true);
        upperRight = domainMap.append(right,false);
        domainMap[lowerLeft].data.space = 2;
        domainMap[upperLeft].data.space = 1;
        domainMap[lowerRight].data.space = 3;
        domainMap[upperRight].data.space = 0;

        // we should init the domain bounding box but it is not us
    }
};


SUITE(gsSelector_gsDomainMap)
{

    TEST(gsBoxesToSelector)
    {
        gsMatrix<real_t> boundingBox(2,2);
        boundingBox<<0,1,0,1;

        gsBoxList bb(2);
        gsMatrix<real_t> dom(2,2);

        dom<<0,0.2,0,0.4;
        bb.append(dom,2);


        dom<<0,0.2,0.4,1;
        bb.append(dom,1);

        dom<<0.2,1,0,0.7;
        bb.append(dom,3);

        dom<<0.2,1,0.7,1;
        bb.append(dom,0);

        gsSelector sel;
        sel.initFromBoxes(bb,boundingBox);

        gsBoxList bl=sel.asBoxList();

        /*std::cout<<"\n\\documentclass{standalone}\\usepackage{tikz}\n\\begin{document}\\begin{tikzpicture}[scale=20]"
            <<bl
            <<"\\end{tikzpicture}\\end{document}\n"<<std::endl;*/
    }


    TEST(gsBoxesToSelector2)
    {
        gsMatrix<real_t> boundingBox(2,2);
        boundingBox<<0,1,0,1;

        gsBoxList bb(2);
        gsMatrix<real_t> dom(2,2);

        dom<<0.3,0.6,0.4,0.8;
        bb.append(dom,2);

        gsDomainMap domainMap;
        domainMap.initFromBoxesMax(bb,boundingBox);

//        std::cout<<"\n\n" << sel.domainMap <<"\n\n"<< std::endl;


        gsBoxList bl=domainMap.asBoxList();

        /*std::cout<<"\n\\documentclass{standalone}\\usepackage{tikz}\n\\begin{document}\\begin{tikzpicture}[scale=20]"
            <<bl
            <<"\\end{tikzpicture}\\end{document}\n"<<std::endl;*/
    }


    TEST(gsBoxesToSelector3)
    {
        gsMatrix<real_t> boundingBox(2,2);
        boundingBox<<0,1,0,1;

        gsBoxList bb(2);
        gsMatrix<real_t> dom(2,2);

        dom<<0.3,0.6,0.4,0.8;
        bb.append(dom,2);
        dom<<0.4,0.9,0.2,0.9;
        bb.append(dom,3);

        gsDomainMap domainMap;
        domainMap.initFromBoxesMax(bb,boundingBox);

//        std::cout<<"\n\n" << sel.domainMap <<"\n\n"<< std::endl;

        gsBoxList bl=domainMap.asBoxList();

        /*std::cout<<"\n\\documentclass{standalone}\\usepackage{tikz}\n\\begin{document}\\begin{tikzpicture}[scale=20]"
            <<bl
            <<"\\end{tikzpicture}\\end{document}\n"<<std::endl;*/
    }


}


SUITE(gsSelector)
{
    TEST (selector1)
    {
        size_t correctAnswers[] = {0, 0, 0, 0, 1, 1, 1, 1, 1, 1};
        
        typedef gsSelector::NodeId NodeId;
        gsDomainMap domainMap;
        NodeId root = domainMap.root();
        domainMap[root].data.dir = 0;
        domainMap[root].data.par = 0.5;
        NodeId left = domainMap.append(root,true);
        NodeId right = domainMap.append(root,false);
        domainMap[left].data.space = 0;
        domainMap[right].data.space = 1;

        gsMatrix<real_t> points(2,10);
        points<<
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
            1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1;

        std::vector<gsBoxList::basisIdT> results;
        domainMap.getBasisAt( points, results );

        CHECK( std::equal( results.begin(), results.end(), correctAnswers ));
    }

    TEST_FIXTURE (data, selector2)
    {
        size_t correctAnswers[] = {0, 1, 2, 3, 0, 0, 3, 3, 1, 3, 0};


        gsMatrix<real_t> evalPoints(2,11);
        evalPoints <<
            0.6, 0.3, 0.2, 0.5, 0.7, 0.4, 0.4, 0.4, 0.1, 0.4, 0.4,
            0.8, 0.4, 0.2, 0.5, 0.7, 0.9, 0.4, 0.1, 0.3, 0.3, 0.7;

        std::vector<gsBoxList::basisIdT> results;
        domainMap.getBasisAt( evalPoints, results );
        CHECK( std::equal( results.begin(), results.end(), correctAnswers ));
    }

    TEST_FIXTURE (data, cutBranchAway)
    {
        // Check children.
        CHECK(domainMap[root].first == left &&
              domainMap[root].second == right &&
              domainMap[left].first == lowerLeft &&
              domainMap[left].second == upperLeft &&
              domainMap[right].first == lowerRight &&
              domainMap[right].second == upperRight);

        // Check parents.
        CHECK(domainMap[lowerLeft].parent == left &&
              domainMap[upperLeft].parent == left &&
              domainMap[lowerRight].parent == right &&
              domainMap[upperRight].parent == right &&
              domainMap[left].parent == root &&
              domainMap[right].parent == root);

        // Erase one branch.
        domainMap.cutBranchAway( right, true );

        // Check children.
        CHECK(domainMap[root].second == upperRight &&
              domainMap[root].first == left &&
              domainMap[left].first == lowerLeft &&
              domainMap[left].second == upperLeft
              );

        // Check parents.
        CHECK(domainMap[upperRight].parent == root &&
              domainMap[lowerLeft].parent == left &&
              domainMap[upperLeft].parent == left &&
              domainMap[left].parent == root);

        // Check erased.
        CHECK(domainMap[lowerRight].parent == right);

        // Erase another branch.
        domainMap.cutBranchAway( left, false );

        // Check children.
        CHECK(domainMap[root].first == lowerLeft &&
              domainMap[root].second == upperRight);

        // Check parents.
        CHECK(domainMap[lowerLeft].parent == root &&
              domainMap[upperRight].parent == root);

        // Check erased.
        CHECK(domainMap[left].parent == lowerRight &&
              domainMap[upperLeft].parent == left &&
              domainMap[lowerRight].parent == right);
        // I (Dominik) am not completely sure whether the order of the
        // erased ones is machine independent.
    }

    gsBoxList::basisIdT getMax( std::vector< gsBoxList::basisIdT > vec )
    {
        if( vec.size() == 0 )
            return 0;
        return *std::max_element( vec.begin(), vec.end() );
    }

    gsBoxList::basisIdT getMin( std::vector< gsBoxList::basisIdT > vec )
    {
        if( vec.size() == 0 )
            return invalidBasisId;
        return *std::min_element( vec.begin(), vec.end() );
    }

    gsBoxList::basisIdT getSingle( std::vector< gsBoxList::basisIdT > vec )
    {
        CHECK( vec.size() == 1 );
        return vec[0];
    }

    void checkManyPoints( const gsMatrix<real_t> points,
                          const gsBoxList boxList,
                          const gsMatrix<real_t> boundingBox,
                          const std::vector<gsBoxList::basisIdT> minResult,
                          const std::vector<gsBoxList::basisIdT> maxResult )
    {
        gsSelector maxSelector, minSelector;
        maxSelector.initFromBoxesMax(boxList,boundingBox);
        minSelector.initFromBoxesMin(boxList,boundingBox);
        gsBoxList maxNewBoxes = maxSelector.asBoxList(  );
        gsBoxList minNewBoxes = minSelector.asBoxList(  );

        for( index_t p = 0; p < points.cols(); ++p )
        {
            std::vector< gsBoxList::basisIdT > boxes = boxList.getBasisAt( points.col(p) );
            CHECK( getMax( boxes ) == maxResult[p] );
            CHECK( getMin( boxes ) == minResult[p] );

            if(     (points.col(p).array() >= boundingBox.col(0).array()).all() &&
                    (points.col(p).array() <  boundingBox.col(1).array()).all() )
            {
                CHECK( maxSelector.getBasisAt( points.col(p) ) == maxResult[p] );
                CHECK( minSelector.getBasisAt( points.col(p) ) == minResult[p] );
                CHECK( getSingle( maxNewBoxes.getBasisAt( points.col(p) ) ) == maxResult[p] );
                CHECK( getSingle( minNewBoxes.getBasisAt( points.col(p) ) ) == minResult[p] );
            }
        }
    }

    void addPoint (gsMatrix<real_t> &m, const gsMatrix<real_t> &p)
    {
        m.conservativeResize(m.rows(), m.cols()+1);
        m.col( m.cols() - 1 ) = p;
    }

    TEST( checkPolicyMax )
    {
        gsBoxList boxList(2);
        gsMatrix<real_t> dom(2,2);

        dom<<0.0,1.0,0.4,0.6;
        boxList.append(dom,1);
        dom<<0.2,0.8,0.2,0.8;
        boxList.append(dom,2);
        dom<<0.4,0.6,0.0,1.0;
        boxList.append(dom,3);

        gsMatrix<real_t> boundingBox(2,2);
        boundingBox << 0.0,1.0,0.0,1.0;

        gsMatrix<real_t> point(2,1);
        gsMatrix<real_t> points(2,0);
        std::vector<gsBoxList::basisIdT> maxResult;
        std::vector<gsBoxList::basisIdT> minResult;

        point << 0.5, 1;
        addPoint(points,point);
        maxResult.push_back(0);
        minResult.push_back(invalidBasisId);

        point << 0.1, 0.9;
        addPoint(points,point);
        maxResult.push_back(0);
        minResult.push_back(invalidBasisId);

        point << 0.4, 0.9;
        addPoint(points,point);
        maxResult.push_back(3);
        minResult.push_back(3);

        point << 0.5, 0.9;
        addPoint(points,point);
        maxResult.push_back(3);
        minResult.push_back(3);


        point << 0.2, 0.7;
        addPoint(points,point);
        maxResult.push_back(2);
        minResult.push_back(2);
        point << 0.3, 0.7;
        addPoint(points,point);
        maxResult.push_back(2);
        minResult.push_back(2);
        point << 0.4, 0.7;
        addPoint(points,point);
        maxResult.push_back(3);
        minResult.push_back(2);
        point << 0.5, 0.7;
        addPoint(points,point);
        maxResult.push_back(3);
        minResult.push_back(2);
        point << 0.6, 0.7;
        addPoint(points,point);
        maxResult.push_back(2);
        minResult.push_back(2);
        point << 0.8, 0.7;
        addPoint(points,point);
        maxResult.push_back(0);
        minResult.push_back(invalidBasisId);

        point << 0.1, 0.6;
        addPoint(points,point);
        maxResult.push_back(0);
        minResult.push_back(invalidBasisId);

        point << 0.0, 0.5;
        addPoint(points,point);
        maxResult.push_back(1);
        minResult.push_back(1);
        point << 0.1, 0.5;
        addPoint(points,point);
        maxResult.push_back(1);
        minResult.push_back(1);
        point << 0.2, 0.5;
        addPoint(points,point);
        maxResult.push_back(2);
        minResult.push_back(1);
        point << 0.3, 0.5;
        addPoint(points,point);
        maxResult.push_back(2);
        minResult.push_back(1);
        point << 0.4, 0.5;
        addPoint(points,point);
        maxResult.push_back(3);
        minResult.push_back(1);
        point << 0.5, 0.5;
        addPoint(points,point);
        maxResult.push_back(3);
        minResult.push_back(1);
        point << 0.6, 0.5;
        addPoint(points,point);
        maxResult.push_back(2);
        minResult.push_back(1);
        point << 0.7, 0.5;
        addPoint(points,point);
        maxResult.push_back(2);
        minResult.push_back(1);
        point << 0.8, 0.5;
        addPoint(points,point);
        maxResult.push_back(1);
        minResult.push_back(1);
        point << 0.9, 0.5;
        addPoint(points,point);
        maxResult.push_back(1);
        minResult.push_back(1);

        point << 0.1, 0.4;
        addPoint(points,point);
        maxResult.push_back(1);
        minResult.push_back(1);
        point << 0.2, 0.4;
        addPoint(points,point);
        maxResult.push_back(2);
        minResult.push_back(1);
        point << 0.8, 0.4;
        addPoint(points,point);
        maxResult.push_back(1);
        minResult.push_back(1);

        point << 0.2, 0.3;
        addPoint(points,point);
        maxResult.push_back(2);
        minResult.push_back(2);
        point << 0.4, 0.3;
        addPoint(points,point);
        maxResult.push_back(3);
        minResult.push_back(2);
        point << 0.6, 0.3;
        addPoint(points,point);
        maxResult.push_back(2);
        minResult.push_back(2);

        point << 0.5, 0.0;
        addPoint(points,point);
        maxResult.push_back(3);
        minResult.push_back(3);

        checkManyPoints( points, boxList, boundingBox, minResult, maxResult );
    }

}
