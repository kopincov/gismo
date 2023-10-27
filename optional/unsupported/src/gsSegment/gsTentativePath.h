/** @file gsTentativePath.h

    @brief Provides declaration of gsTentativePath class. A candidate cutting
    loop.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): D.-M. Nguyen
*/

#pragma once

#include <gsUtils/gsGraph/gsGraphPath.h>

namespace gismo
{

  /** 
     @brief
 Class 
  */

  
template<class T>
class gsTentativePath: public gsGraphPath<T>
{

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Data members
protected:
  //todo: check if the line below is possible
  //gsDataStructure<T>* solid; //a data structure from which the edge graph is derived, candidates: gsMesh, gsHeMesh, gsSolid
  gsSolid<T>* solid; //a data structure from which the edge graph is derived, candidates: gsMesh, gsHeMesh, gsSolid
  unsigned tentativeVert; 
  gsMatrix<T> * adjMatrix; // adjacent matrix, only need upper half 
  /// the following are just NULL pointers in case of traditional graph
  gsVolSegOptions<T> * vsOption; // options for Volume Segmentation  
  
private:
  gsVector<T> RefPlaneCoefs; // reference plane coefficients ax+by+cz+d=0
  std::vector<unsigned> Path2tttVert; // path from *source* to *tentativeVert*: (v, ) tentativeVert,..., u_adjacent_to_source, source, target

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Constructors
public:

  /// Default empty constructor
  gsTentativePath(): gsGraphPath<T>() { };
  
  /// Construct Base
  gsTentativePath(int const & s, int const & t, gsVector<unsigned>::Ptr prev,gsMatrix<T>* adjMatrixI,gsVolSegOptions<T>* vsOptionI): gsGraphPath<T>(s,t,prev)
  {adjMatrix = adjMatrixI; vsOption = vsOptionI;std::vector<unsigned> Path2tttVert2;gsVector<T> RefPlaneCoefs2;}
  
  /// Sets this object's solid member.
  /// \param solidI a solid
  void setDataStructure(gsSolid<T>* solidI){solid=solidI;}
  
  /// Sets this object's tentativeVert member.
  /// \param tentativeVertI the index of a vertex.
  void setTentativeVertex(unsigned const & tentativeVertI){tentativeVert = tentativeVertI;}
  
  /// Compute the reference plane.
  void setRefPlaneCoefs()
  {   
    assert(Path2tttVert.empty()==false); // Make sure to setPath2tttVert() before executing this members
    unsigned u_adjacent_to_source = *(Path2tttVert.end()-3);
    gsVector3d<T> normal_to_RefPlane = ( solid->getVertexFromID(this->target)->coords-solid->getVertexFromID(this->source)->coords ).cross(
      solid->getVertexFromID(u_adjacent_to_source)->coords - solid->getVertexFromID(this->source)->coords ); 
    gsVector<T> abcd(4); // coefficients a, b, c, d
    abcd(0) = normal_to_RefPlane.x();
    abcd(1) = normal_to_RefPlane.y();
    abcd(2) = normal_to_RefPlane.z();
    abcd(3) = -( normal_to_RefPlane.x()*(solid->getVertexFromID(this->source)->coords).x() + 
      normal_to_RefPlane.y()*(solid->getVertexFromID(this->source)->coords).y() + 
      normal_to_RefPlane.z()*(solid->getVertexFromID(this->source)->coords).z()
    );
    RefPlaneCoefs = abcd;
    
    #ifdef debug
    #if (debug==3)
	std::cout<<"\n Normal to the reference plane: "<< normal_to_RefPlane;
	std::cout<<"\n Coefficients a, b, c, d: "<< abcd;
	std::cout<<std::endl;
    #endif
    #endif 
  };
  
  /// Compute Path2tttVert, the path as a list of vertices
  void setPath2tttVert(){std::vector<unsigned> path=this->computePath(this->previous,tentativeVert);path.push_back(this->target);Path2tttVert = path;};  

  /// Destructor
  ~gsTentativePath() { };

  
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// const members
public:  
  /// this is to specify whether a new edge [bestu,v], where bestu is the lastest found vertex, satisfying the constraints for the Dijkstra
  /// algorithm. That is whether this new edge is on the same face with the previous edges in the *path* from vert1 to bestu, 
  /// also whether it is on the same face with the edge [vert1 vert2] 
  bool isDiffFaces(unsigned const& v) const;   
  
  /// This is to calculate the weighted distance from a vertex v to the plane (vert2, path[0] = vert1, and path[1])
  T evalPlanarityCost(unsigned const & v) const
  {
    // the reference plane pertains to: vert2, vert1, and u_adjacent_to_vert1
    // define the cost as the total distance of the vertices in the set [path(2:end),v]    
    assert(RefPlaneCoefs.cols()>0); // check if RefPlaneCoefs has been initialized before  
    T x0 = (solid->getVertexFromID(v)->coords).x();
    T y0 = (solid->getVertexFromID(v)->coords).y();
    T z0 = (solid->getVertexFromID(v)->coords).z();
    T Distance = ( RefPlaneCoefs(0)*x0 + RefPlaneCoefs(1)*y0 + RefPlaneCoefs(2)*z0 + RefPlaneCoefs(3))/ RefPlaneCoefs.head(3).squaredNorm();
    T RelativeDistance;
    if (vsOption->PlanarityCostType==1)
      RelativeDistance = math::abs(Distance/( (solid->getBoundingBox())).getMaxSize() );
    else
    {
        RelativeDistance = 0;
        gsWarn <<"\n not yet\n";assert(0);
    };
    return RelativeDistance;
  };  


// Data members
private:


}; // class gsTentativePath


//////////////////////////////////////////////////
//////////////////////////////////////////////////

template <class T>
bool gsTentativePath<T>::isDiffFaces(unsigned const& v) const
{
  // check if the edge u-v is on the same face with any from *path*
  bool isDiffFaces=true;
  std::vector< gsSolidHalfFace<T>* > uvFace;
  std::vector< gsSolidHalfFace<T>* > pathFace;
  
  assert(Path2tttVert.empty()==false); // Make sure to setPath2tttVert() before executing this members
  
  if ( ( (*adjMatrix)(tentativeVert,v) == vsOption->costConvexEdge ) || ( (*adjMatrix)(tentativeVert,v) == vsOption->costNonConvexEdge ) )
    uvFace = solid->getVertexFromID(tentativeVert)->getFacesContaining2Vertices(solid->getVertexFromID(v),true);
  else uvFace = solid->getVertexFromID(tentativeVert)->getFacesContaining2Vertices(solid->getVertexFromID(v),false);
  
  for (unsigned i=0;i<=(Path2tttVert.size()-1)-1;i++)
  {    
    pathFace.clear();    
    if ( ( (*adjMatrix)(Path2tttVert[i],Path2tttVert[i+1]) == vsOption->costConvexEdge ) || ( (*adjMatrix)(Path2tttVert[i],Path2tttVert[i+1]) == vsOption->costNonConvexEdge ) )      
      pathFace = solid->getVertexFromID(Path2tttVert[i])->getFacesContaining2Vertices(solid->getVertexFromID(Path2tttVert[i+1]),true);
    else pathFace = solid->getVertexFromID(Path2tttVert[i])->getFacesContaining2Vertices(solid->getVertexFromID(Path2tttVert[i+1]),false);
    
    // now check if uvFace and pathFace intersects each other    
    for (unsigned i2=0;i2!=uvFace.size();i2++)
    {
      for (unsigned j=0;j!=pathFace.size();j++)
      {
	if (uvFace[i2]==pathFace[j]) { isDiffFaces=false; break; };
      };
    };
  };
  return isDiffFaces;
};


}; // namespace gismo
