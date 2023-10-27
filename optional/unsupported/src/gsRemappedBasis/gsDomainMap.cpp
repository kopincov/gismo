/** @file gsDomainMap.cpp

    @brief Implementation of gsDomainMap, which is used in gsSelector.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
*/

#include <gsRemappedBasis/gsDomainMap.h>
#include <gsRemappedBasis/gsDomainMapApply.hpp>

namespace gismo {


void gsDomainMap::getBasisAt (const gsMatrix<real_t> &points, std::vector<basisIdT> &out) const
{
    out.resize(points.cols());
    for (index_t p=0; p<points.cols();++p)
        out[p]=getBasisAt(points.col(p));
}




gsMatrix<real_t> gsDomainMap::getBoundingBox(NodeId node, const gsMatrix<real_t> &domain) const
{
    gsMatrix<real_t> result;
    result=domain;
    const_node_iterator_NFS pre(*this, node);
    const_node_iterator_NFS cur(*this, pre->parent);

    while(cur!=const_node_iterator_NFS(*this,no_node))
    {
        GISMO_ASSERT(cur->data.dir<domain.rows(), "dimension mismatch");
        if(cur->first==pre)
            result( cur->data.dir,1 )=math::min(result(cur->data.dir,1 ),cur->data.par);
        else // (cur->second==pre)
            result(cur->data.dir,0 )=math::max(result(cur->data.dir,0 ),cur->data.par);
        pre=cur;
        cur=const_node_iterator_NFS(*this,cur->parent);
    }
    return result;
}

gsBoxList gsDomainMap::asBoxList(const gsMatrix<> &boundingBox) const
{
    gsBoxList result(static_cast<directionT>(boundingBox.rows()));
    getBoxes<> getter(result);
    applyToSubtree(root(),boundingBox,getter);
    return result;
}


template <typename Action, typename Arg>
void applyToDomWithMatAndArg (gsDomainMap &dom, const typename selectBoxT<Action>::BoxT &mat, Arg arg )
{
    Action act(arg);
    dom.apply(mat,act);
}

template <template <int dim1,int dim2> class op, typename MatrixT, typename Arg>
void dispatchOpToDimWithMat ( gsDomainMap &dom, const MatrixT &mat, Arg arg, index_t domainDim)
{
    switch (domainDim)
    {
    case 1:
        applyToDomWithMatAndArg<op<1,2> >(dom,mat,arg);
        break;
    case 2:
        applyToDomWithMatAndArg<op<2,2> >(dom,mat,arg);
        break;
    case 3:
        applyToDomWithMatAndArg<op<3,2> >(dom,mat,arg);
        break;
    case 4:
        applyToDomWithMatAndArg<op<4,2> >(dom,mat,arg);
        break;
    default:
        applyToDomWithMatAndArg<op<-1,-1> >(dom,mat,arg);
        break;
    }
}



void gsDomainMap::initFromBoxesMax(
        const gsBoxList  &boxes,
        const gsMatrix<> &boundingBox
        )
{
    GISMO_ASSERT(boundingBox.rows()==boxes.domainDim(),"dimension mismatch");
    boxes.check();

    m_boundingBox=boundingBox;
    addBoxesMax(boxes);
}

void gsDomainMap::addBoxesMax(
        const gsBoxList  &boxes
        )
{
    const directionT domainDim=m_boundingBox.rows();
    GISMO_ASSERT(domainDim==boxes.domainDim(),"dimension mismatch");
    const basisIdT    maxLev=boxes.maxId();
    (*this)[m_root].data.space=0;

    for (basisIdT lvl=1; lvl<maxLev+1; ++lvl)
    {
        for (size_t b=0;b<boxes.size();++b)
        {
            if (boxes.basisId(b)!=lvl || empty(boxes.box(b)) )
                continue;
            dispatchOpToDimWithMat<setBasisMax>(*this,boxes.box(b),lvl,domainDim);
        }
    }
}

void gsDomainMap::initFromBoxes(
        const gsBoxList  &boxes,
        const gsMatrix<> &boundingBox
        )
{
    const directionT domainDim=boundingBox.rows();
    GISMO_ASSERT( domainDim==boxes.domainDim(),"dimension mismatch");
    boxes.check();

    m_boundingBox=boundingBox;
    (*this)[m_root].data.space=invalidBasisId;
    for (size_t b=0;b<boxes.size();++b)
    {
        dispatchOpToDimWithMat<setBasis>(*this,boxes.box(b),boxes.basisId(b),domainDim);
    }
}

void gsDomainMap::initConstant (
        basisIdT b,
        const gsMatrix<> &boundingBox
        )
{
    m_boundingBox=boundingBox;
    (*this)[m_root].data.space=b;
}

void gsDomainMap::initFromBoxesMin(
        const gsBoxList  &boxes,
        const gsMatrix<> &boundingBox
        )
{
    GISMO_ASSERT(boundingBox.rows()==boxes.domainDim(),"dimension mismatch");
    boxes.check();

    m_boundingBox=boundingBox;

    (*this)[m_root].data.space=invalidBasisId;
    addBoxesMin(boxes);
}

void gsDomainMap::addBoxesMin(
        const gsBoxList  &boxes
        )
{
    const directionT domainDim=m_boundingBox.rows();
    GISMO_ASSERT(domainDim==boxes.domainDim(),"dimension mismatch");
    const basisIdT    maxLev=boxes.maxId();
    for (basisIdT lvl=1; lvl<maxLev+1; ++lvl)
    {
        for (size_t b=0;b<boxes.size();++b)
        {
            if (boxes.basisId(b)!=lvl || empty(boxes.box(b)) )
                continue;
            dispatchOpToDimWithMat<setBasisMin>(*this,boxes.box(b),lvl,domainDim);
        }
    }
}




template <typename Action, typename Arg>
void applyToDomWithArg (gsDomainMap &dom, Arg arg )
{
    Action act(arg);
    dom.apply(act);
}

template <typename Action, typename Arg>
void applyToDomWithArg (const gsDomainMap &dom, Arg arg )
{
    Action act(arg);
    dom.apply(act);
}

template <template <int dim1,int dim2> class op, typename Arg>
void dispatchOpToDim ( gsDomainMap &dom, gsDomainMap::basisIdT arg, index_t domainDim)
{
    switch (domainDim)
    {
    case 1:
        applyToDomWithArg<op<1,2> >(dom,arg);
        break;
    case 2:
        applyToDomWithArg<op<2,2> >(dom,arg);
        break;
    case 3:
        applyToDomWithArg<op<3,2> >(dom,arg);
        break;
    case 4:
        applyToDomWithArg<op<4,2> >(dom,arg);
        break;
    default:
        applyToDomWithArg<op<-1,-1> >(dom,arg);
        break;
    }
}


template <template <int dim1,int dim2> class op>
gsDomainMap dispatchBoolOpToDim ( const gsDomainMap & left,
                                  const gsDomainMap & right )
{
    gsDomainMap::directionT domainDim = left.getBoundingBox().rows();
    GISMO_ASSERT(domainDim==right.getBoundingBox().rows(),"Dimension mismatch");
    gsDomainMap result=left;
    switch (domainDim)
    {
    case 1:
        applyToDomWithArg<applyToOther<op, 1, 2>, gsDomainMap &>(right,result);
        break;
    case 2:
        applyToDomWithArg<applyToOther<op, 2, 2>, gsDomainMap & >(right,result);
        break;
    case 3:
        applyToDomWithArg<applyToOther<op, 3, 2>, gsDomainMap & >(right,result);
        break;
    case 4:
        applyToDomWithArg<applyToOther<op, 4, 2>, gsDomainMap & >(right,result);
        break;
    default:
        applyToDomWithArg<applyToOther<op, -1, -1>, gsDomainMap & >(right,result);
        break;
    }
    return result;
}


GISMO_EXPORT gsDomainMap shiftIds ( const gsDomainMap & orig, int amount)
{
    gsDomainMap result=orig;
    dispatchOpToDim<gsDomainMap::shiftIds,int>(result,amount,orig.getBoundingBox().rows() );
    return result;
}


GISMO_EXPORT gsDomainMap polytopeUnion( const gsDomainMap & left,
                                        const gsDomainMap & right )
{
    return dispatchBoolOpToDim<gsDomainMap::setBasisMax>(left,right);
}


GISMO_EXPORT gsDomainMap polytopeIntersect( const gsDomainMap & left,
                                            const gsDomainMap & right )
{
    return dispatchBoolOpToDim<gsDomainMap::setBasisMin>(left,right);
}

GISMO_EXPORT gsDomainMap polytopeDifference( const gsDomainMap & left,
                                             const gsDomainMap & right )
{
    return dispatchBoolOpToDim<gsDomainMap::setDifference>(left,right);
}

} // namespace gismo
