/** @file gsHeVertex.h

    @brief Provides gsHeVertex class for a vertex of a gsHeMesh

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsUtils/gsHeMeshElement.h>
#include <gsUtils/gsHalfEdge.h>


namespace gismo {

template <class T>  class gsHalfEdge;

template <class T> 
class gsHeVertex  : public gsVertex<T>
{
private:
    typedef gsVertex<T> Base;
    typedef T scalar_t;
    typedef typename gsHeMeshElement<T>::gsHalfEdgeHandle gsHalfEdgeHandle;
private:
    int classification; //-1 = undefined
    gsHalfEdgeHandle edge;
    unsigned int valence;

public:
    gsHeVertex(T x, T y, T z = 0) : Base(x,y,z) {classification = -1;}

    virtual ~gsHeVertex(){ }

    inline int getClassification() { return classification;}

    void setClassification(int i)
    {
        if(i == -1 || i >= 0)
            classification = i;
        else
            std::cout << "Wrong classification";
    }

    inline gsHalfEdgeHandle getHalfEdge(){return edge;}
    inline void setHalfEdge(gsHalfEdgeHandle heh){edge = heh; classification = 0;}
    inline int getValence(){return valence;}
    inline void setValence(int i){valence = i;}

    bool equal(gsHeVertex<T> *rhs)
    {
        return (((this->x())==rhs->x())&&
                ((this->y())==rhs->y())&&
                ((this->z())==rhs->z()));
    }

    //using Base::operator==;

/*    bo4ol operator != (gsHeVertex<T> const & rhs)
    {
        return (((this->x())!=rhs.x())||
                ((this->y())!=rhs.y())||
                ((this->z())!=rhs.z()));
    }
*/

};

}// namespace gismo


