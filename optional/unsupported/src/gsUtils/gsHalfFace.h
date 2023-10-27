/** @file gsHalfFace.h

    @brief Provides gsHalfFace class for a face of a gsHeMesh

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen
*/

#pragma once

#include <list>
#include <gsUtils/gsHeMesh.h>
#include <gsUtils/gsHeMeshElement.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo 
{
    
template <class T>
class gsHalfFace : public gsHeMeshElement<T>
{
private:
    typedef gsHeMeshElement<T> MeshElement;
    typedef typename MeshElement::scalar_t scalar_t;
    typedef typename MeshElement::gsHeVertexHandle gsVertexHandle;
    typedef typename MeshElement::gsHalfEdgeHandle gsHalfEdgeHandle;
    typedef typename MeshElement::gsHalfFaceHandle gsHalfFaceHandle;
    typedef typename MeshElement::gsCellHandle gsCellHandle;
    typedef gsVector3d<T> gsVector;
    typedef gsVector * gsVectorHandle;
    
private:

    std::vector<gsVertexHandle> vertices;
    gsHalfFaceHandle mate;
    gsHalfFaceHandle next; //prev;
   // gsHalfEdgeHandle boundary;    // todo: make a *loop* member as done with gsSolid
    std::vector<gsHalfEdgeHandle> boundary;
    gsCellHandle parent;    
    int classification; //-1 = undefined; 0 = domain; 1 = extended; 2 = outer
    bool regular;

public:
    gsHalfFace(): MeshElement() {this->classification = -1;}

   // gsHalfFace(gsHalfEdgeHandle const & he, int const& i): MeshElement(i), boundary(he){ };
    
    gsHalfFace(std::vector<gsVertexHandle> const & verts, int const& i): MeshElement(i){this->vertices = verts; this->classification = -1; this->regular = false;}
    
    gsHalfFace(std::vector<gsHalfEdgeHandle> const &  hedges, int const& i)
        : MeshElement(i), boundary(hedges[0])
        {
            this->setBoundary( hedges );
            for ( typename std::vector<gsHalfEdgeHandle>::iterator
                  it = hedges.begin(); it!= hedges.end(); ++it)
            {
                gsHalfEdgeHandle edge = *it;
                this->vertices.push_back(edge->getSource());
            }
            this->classification = -1;
            this->regular = false;
        }

    virtual ~gsHalfFace();

    //getter-methods
    inline std::vector<gsVertexHandle> getVertices() {return vertices;}
    inline gsHalfFaceHandle getMate() {return mate;}
    inline gsHalfFaceHandle getNext() {return next;}
    inline std::vector<gsHalfEdgeHandle> getBoundary() {return boundary;}
    inline void setBoundary(std::vector<gsHalfEdgeHandle> const &  hedges){boundary = hedges;}
    inline gsCellHandle getParent() {return parent;}
    inline void setVertex(gsVertexHandle v) {vertices.push_back(v);}
    inline int getClassification() { return classification;}

    void setClassification(int i)
    {
        if(i == 0 || i == 1 || i == 2)
            classification = i;
        else
            std::cout << "Wrong classification";
    }
    inline bool getRegular() { return regular;}
    inline void setRegular(bool r) { regular = r;}

    /*
    // set values of the members: next, prev, face of the half-edges in the boundary of a face
    void setBoundary(std::vector<gsHalfEdgeHandle> const &  hedges)
    {

            typename std::vector<gsHalfEdgeHandle>::const_iterator p = hedges.begin();
            this->boundary = *p;
            (*p)->prev = hedges.back();
            (*p)->face = this;	    
            for ( typename std::vector<gsHalfEdgeHandle>::const_iterator 
                      it = hedges.begin()+1; it!= hedges.end(); ++it)
            {   
                (*it)->face = this;
                (*it)->prev = *p;
                (*p)->next  = *it;
                p++;
            }
            (*p)->next = hedges.front();

        gsHalfEdgeHandle ne = *hedges.begin();
        gsHalfEdgeHandle pr;
        ne->setFace(this);
        for ( typename std::vector<gsHalfEdgeHandle>::const_iterator
                            it = hedges.begin()+1; it!= hedges.end(); ++it)
        {
            ne->next = *it;
            this->boundary.push_back(ne);
            pr = ne;
            ne = *it;
            ne->face = this;
            ne->prev = pr;
        }
        ne->next = hedges.front();
        this->boundary.push_back(ne);
    }
*/
    std::ostream &print(std::ostream &os) const
    {
        os<<"gsHalfFace \n" ;
        return os;
    }
 
};

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////


template <typename T>
gsHalfFace<T>::~gsHalfFace()
{
    // delete boundary half-edges
}
    

}// namespace gismo
