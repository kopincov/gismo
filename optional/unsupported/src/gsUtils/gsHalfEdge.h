/** @file gsHalfEdge.h

    @brief Provides gsHalfEdge class for an edge of a gsHeMesh

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen
*/

//http://www.cs.purdue.edu/homes/cmh/distribution/books/geo.html
// http://www.cs.utah.edu/~xchen/euler-doc/

#ifndef _HALFEDGE_H_
#define _HALFEDGE_H_

#include <iostream>

#include <gsUtils/gsHeMeshElement.h>
#include <gsUtils/gsHalfEdge.h>
#include <gsNurbs/gsKnotVector.h>


namespace gismo {


enum gsHEdSign{ PLUS=true, MINUS=false };

template <class T>
class gsHalfEdge : public gsHeMeshElement<T>
{
private:
    typedef gsHeMeshElement<T> MeshElement;
    typedef typename MeshElement::scalar_t scalar_t;
    typedef typename MeshElement::gsHeVertexHandle gsVertexHandle;
    typedef typename MeshElement::gsHalfEdgeHandle gsHalfEdgeHandle;
    typedef typename MeshElement::gsHalfFaceHandle gsFaceHandle;
    typedef gsVector3d<T> gsVector;
    typedef gsVector * gsVectorHandle;     

private:

    //gsEdgeHandle edge;#
    gsHalfEdgeHandle mate;

    gsVertexHandle source;
    gsHalfEdgeHandle next;

    gsFaceHandle face;// parent face

    // Additional properties
//    bool is_convex;
//    gsKnotVector<T>* kv; // Knot vector and degree of the trimming line (in parameter domain, not the trimming curve in the trimming surface)
//    gsMatrix<T>* cp; // Control points of  the trimming line



public:
    gsHalfEdge(const gsVertexHandle & v,
               const gsFaceHandle   & f,
               const int            & i,
               const bool           & // conv
               ) : MeshElement(i), source(v), face(f)
    {
        mate = 0;
        next = 0;
    }

    gsHalfEdge() { };
    virtual ~gsHalfEdge(){ };

    //getter methods{
    inline gsHalfEdgeHandle getMate(){return mate;}
    inline gsVertexHandle getSource(){return source;}
    inline gsHalfEdgeHandle getNext(){return next;}
    inline gsFaceHandle getFace(){return face;}
    inline void setMate(gsHalfEdgeHandle heh){mate = heh;}
    inline void setSource(gsVertexHandle vh){source = vh;}
    inline void setNext(gsHalfEdgeHandle heh){next = heh;}
    inline void setFace(gsFaceHandle fh){face = fh;}

/*    std::ostream &print(std::ostream &os) const
        {
            //os<<"gsHalfEdge from "<< *source ;
            return os;
        }; */

   /* bool operator < (gsHalfEdge const & rhs) const
    {
      return ( Xless<T>(this->source,rhs.source)

              ||
           ( (this->source->x() == rhs.source->x() &&
              this->source->y() == rhs.source->y() &&
              this->source->z() == rhs.source->z() ) &&
             Xless<T>(this->target,rhs.target) ));
    }
    */

    bool operator == (gsHalfEdge const & rhs) const
    {
        return (((this->source->x())==rhs.source->x())&&
                ((this->source->y())==rhs.source->y())&&
                ((this->source->z())==rhs.source->z())&&
                ((this->target->x())==rhs.target->x())&&
                ((this->target->y())==rhs.target->y())&&
                ((this->target->z())==rhs.target->z()));
    }
   /* bool operator> (gsHalfEdge const & rhs) const
    {
        return !(*this<rhs)&&!(*this==rhs);
    }
    */
    bool operator != (gsHalfEdge const & rhs) const
    {
        return !(*this==rhs);
    }

// template <class T>
// addHe(typename gsHeMeshElement<Vertex>::gsVertexHandle *v, gsHalfEdge<T> *h, gsHEdSign sign)
//  {
//     gsHalfEdge<T> *he;
//     if (h->edge == 0)
//         he = h;
//     else {
//         he = new gsHalfEdge<T>();
//         h->prev->next = he;
//         he->prev = h->prev;
//         h->prev = he;
//         he->next = h;
//     }

//     he->edge = e;
//     he->origin = v;
//     he->loop = h->loop;
//     if(sign == PLUS)
//         e->hed1 = he;
//     else 
//         e->hed2 = he;

//     return he;
// }

// template <class T>
// gsHalfEdge<T>* delHe(gsHalfEdge<T> *he)
//  {
//     if(he->edge == 0) {
//         delete he;
//     } else if (he->next == he) {
//         he->edge = 0;
//         return he;
//     }
//     else {
//         he->prev->next = he->next;
//         he->next->prev = he->prev;
//         delete he;
//         return he->prev;
//     }
// }

};//namespace gismo
};
#endif
