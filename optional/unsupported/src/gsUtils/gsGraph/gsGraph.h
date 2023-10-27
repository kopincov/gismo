#pragma once

#include <gsUtils/gsGraph/gsGraphPath.h>
#include <gsUtils/gsGraph/gsGraphOptions.h>
#include <gsModeling/gsSolid.h>
#include <gsSegment/gsTentativePath.h>

namespace gismo
{

  /** 
      @brief Traditional Dijkstra algorithm and its extension to volume segmentation problem
  */

  
template<class T>
class gsGraph
{
public:
  typedef gsMatrix<T> gsMatrixT;
    
/// Data members
protected:
  gsMatrixT adjMatrix; // adjacent matrix, todo: only need upper half
  gsSolid<T> * solid;  
  
  gsGraphOptions<T> option;
  gsVolSegOptions<T> vsOption; // options for Volume Segmentation

public:
  /// Default empty constructor
  gsGraph()
  {
      solid = 0;
  }

  gsGraph(bool defaultExample);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Dijkstra algorithm
public:
  gsMatrix<T>& adjMat() {return adjMatrix;}
  gsGraphOptions<T>& graphOp() {return option;}
  gsVolSegOptions<T>& volSegOp() {return vsOption;}
  gsSolid<T>* getSolid() const {return solid;}

  /// traditional Dijkstra's algorithm
  virtual gsGraphPath<T>* Dijkstra(unsigned vert1,unsigned vert2)
  {return Dijkstra(vert1,vert2,0); };

  /// Dijkstra algorithm
  gsGraphPath<T>* Dijkstra(unsigned vert1,unsigned vert2, gsSolidHalfEdge<T>* he);
  
  /// Prints the object as a string.
  std::ostream &print(std::ostream &ost) const;   

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Local auxiliary members
private:
  /// augment Tentative Distance or not, based on the InfRep method to determine if *du + dist_uv < dv*
  bool augmentTentativeDistance(T du,T dist_uv,T dv);
  
  /// true if dist_v < dist_u based on the InfRep method
  bool isSmaller(T dist_v,T dist_u);
  
  /// Asign the value *value* a symmetric matrix *mat* values of entries given by rows *ind1* and columns *ind2* 
  static void SetMatrixEntries(gsMatrix<T>* mat,std::vector<unsigned> const & ind1,std::vector<unsigned> const & ind2,T value=-1);
  
  /// Asign the values *value* a symmetric matrix *mat* values of entries given by rows *ind1* and columns *ind2* 
  static void SetMatrixEntries(gsMatrix<T>* mat,std::vector<unsigned> const & ind1,std::vector<unsigned> const & ind2,std::vector<T> const & values);  

}; // class gsGraph


//=============================================================================
// SOURCE
//=============================================================================

// Define the functionality of the operator << acting on a class
template <class T>
std::ostream &operator<<(std::ostream &os, const gsGraph<T>& b)
{return b.print(os); };

template <class T>
std::ostream &gsGraph<T>::print(std::ostream &ost) const
    {
	ost << "gsGraph: \n";	  	  
	ost << "\tAdjacency matrix: \n" << this->adjMatrix << std::endl;
	return ost;
    };  

template <class T>
gsGraph<T>::gsGraph(bool defaultExample)
{
    if ( defaultExample==true )
    {
        // the example from http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm; source=0, target=4;
        gsMatrixT am(6,6);
        am.setOnes(); am = -am;
        // setting up the upper half first and then the lower one
        am(0,1)=7;am(0,2) = 9;am(0,5) = 14;
        am(1,0)=7;am(2,0) = 9;am(5,0) = 14;
        am(1,2)=10;am(1,3)=15;
        am(2,1)=10;am(3,1)=15;
        am(2,3)=11;am(2,5)=2;
        am(3,2)=11;am(5,2)=2;
        am(3,4)=6;
        am(4,3)=6;
        am(4,5)=9;
        am(5,4)=9;

        adjMatrix.swap( am );

        solid = 0;
    }
}; 

template <class T>
bool gsGraph<T>::augmentTentativeDistance(T du,T dist_uv,T dv)
{
  if (option.InfRep == -1 )
  {
    if (dist_uv==-1 || du==-1) return false;
    if (dv==-1) return true;
    return ( du + dist_uv < dv );    
  }
  else
  {
    gsWarn << "\n Not yet";
    exit(1);    
  };
};

template <class T>
bool gsGraph<T>::isSmaller(T dist_v,T dist_u)
{
  if (dist_v==-1) return false;
  if (dist_u==-1) return true;
  return ( dist_v < dist_u ); 
};

template <class T>
void gsGraph<T>::SetMatrixEntries(gsMatrix<T>* mat,std::vector<unsigned> const & ind1,std::vector<unsigned> const & ind2,T value)
{
  assert(ind1.size()==ind2.size());
  for (unsigned i=0;i!=ind1.size();i++)
  {
    (*mat)(ind1[i],ind2[i]) = value;
    (*mat)(ind2[i],ind1[i]) = value;
  };
};

template <class T>
void gsGraph<T>::SetMatrixEntries(gsMatrix<T>* mat,std::vector<unsigned> const & ind1,std::vector<unsigned> const & ind2,std::vector<T> const & values)
{
  assert(ind1.size()==ind2.size());
  assert(ind1.size()==values.size());
  for (unsigned i=0;i!=ind1.size();i++)
  {
    (*mat)(ind1[i],ind2[i]) = values[i];
    (*mat)(ind2[i],ind1[i]) = values[i];
  };  
};

template <class T>
gsGraphPath<T>* gsGraph<T>::Dijkstra(unsigned vert1,unsigned vert2, gsSolidHalfEdge<T> * he)
{
    T Infinity;  
    if (option.InfRep == -1 ) Infinity = -1; else { gsWarn << "\n Not yet"; exit(1); };
    
    unsigned nv = adjMatrix.rows(); // number of vertices
    gsVector<T> dist(nv); // distance vector
    gsVector<unsigned>::Ptr prev(new gsVector<unsigned>(nv)); // previous vector    
    gsVector<T>  plaCost(nv);// planarity costs
    
    // initialization for dist and prev
    for (unsigned i=0;i<=nv-1;i++)
    {
        if ((solid->vertex)[i]->hed->face->vol==he->face->vol) dist(i) = adjMatrix(vert1,i);
        else dist(i) = Infinity;
      (*prev)(i) = vert1;
      plaCost(i) = 0.;
    };    
    
    // For DiffFaces
    std::vector<unsigned> ind1;
    std::vector<unsigned> ind2;
    std::vector<T> store_ajdMat;    
    if ( vsOption.DiffFaces==true )
    {
      std::vector< gsSolidHeVertex<T>* > V12 = he->face->getVertices();   
      for (unsigned i=0; i!=V12.size();i++)
      {
	dist( V12[i]->getId() ) = Infinity;
      };
      V12.clear();
      V12 = he->mate->face->getVertices();   
      for (unsigned i=0; i!=V12.size();i++)
      {
	dist( V12[i]->getId() ) = Infinity;
      };      
     
      // assign edges in the same face with vert2, but not in the same face with vert1, and do not containing vert2 to
      // be inf, this is to make sure the near last edge can be extended to vert2   
      // - find faces incident to vert1 and faces incident to vert2
      std::vector< gsSolidHalfFace<T>* > vert2Faces = solid->getVertexFromID(vert2)->getHalfFaces();
      // he->face, he->mate->face
      std::vector<gsSolidHeVertex<T>*> facevert;
      vert2Faces.push_back(vert2Faces.at(0)); // a trick to make a cycle
      for (typename std::vector< gsSolidHalfFace<T>* >::const_iterator it=vert2Faces.begin();it!=vert2Faces.end();++it)
      {
	if ( (*it!=he->face) && (*it!=he->mate->face) )
	{
	  facevert = (*it)->getVertices();
	  for (typename std::vector<gsSolidHeVertex<T>*>::const_iterator iv1=facevert.begin();iv1<=facevert.end()-1-1;++iv1)
	  {
	    for (typename std::vector<gsSolidHeVertex<T>*>::const_iterator iv2=iv1+1;iv2<=facevert.end()-1;++iv2)
	    {
	      if ( ( *iv1!=solid->getVertexFromID(vert2) ) && ( *iv2!=solid->getVertexFromID(vert2) ) )
	      {		
		ind1.push_back( (*iv1)->getId() );
		ind2.push_back( (*iv2)->getId() );
		store_ajdMat.push_back( adjMatrix((*iv1)->getId(), (*iv2)->getId()) );
	      };
	    };
	  };
	  //*it
	}
      };
      SetMatrixEntries(&adjMatrix,ind1,ind2,Infinity);
    };  
    
    std::vector<unsigned> TentativeVert;
    for (unsigned i=0;i!=nv;i++){ if ( ((solid->vertex)[i]->hed->face->vol==he->face->vol) && (i!=vert1) ) TentativeVert.push_back(i); }


    //TentativeVert.erase( TentativeVert.begin() + vert1 ); // remove *vert1* from the the set of tentative vertices
    
    T planarCost_v(0);
    
    gsTentativePath<T> tttPath(vert1,vert2,prev,&adjMatrix,&vsOption);
    tttPath.setDataStructure(solid);
    
    /// search
    T planar_weight=vsOption.PlanarityCostWeight;
    unsigned  u;
    T du;
    unsigned iu;
    bool IsDiffFace;
    unsigned v;
    while (TentativeVert.empty()==false)
    {
        // find *u* in *TentativeVert* such that *d(u) = min{d(z): z in TentativeVert}*
        u = TentativeVert[0];
        du = dist(u);
        iu=0; // position of u in T
        for (unsigned i=1;i!=TentativeVert.size();i++)
        {
            v=TentativeVert[i];
            if ( isSmaller(dist(v),du) )
            {
                u = v; du = dist(v); iu=i;
            }
        }

#ifdef debug
#if (debug==10)
        gsWarn <<"\n The shortest path found is from " << vert1 <<" to the vertex: " << u;
        gsWarn <<"\n with cost: " << du << std::endl;
#endif
#endif
        // terminate the *while* loop if all remaining vertices are inaccessable from the source *vert1*
        if (option.InfRep==-1)
        {
            GISMO_ASSERT(du!=-1,"\n All remaining vertices are inaccessable from source");
        }
        else
        {
            gsWarn << "\n Pls modify the terminating command to make it suitable for the new *InfRep* value\n";
            exit(1);
        };

        // break the *while* loop if *u* is *vert2*
        if (u==vert2) break;

        // remove *u* from *TentativeVert*
        TentativeVert.erase(TentativeVert.begin()+iu);

        // code for ignoredPath here

        /// augment tentative distances
        IsDiffFace=true;

        tttPath.setPrevious(prev);
        tttPath.setTentativeVertex(u);
        tttPath.setPath2tttVert();
        tttPath.setRefPlaneCoefs(); // todo: can be improved by only computing this for each branch enamating from *vert1*

        for (unsigned i=0; i!=TentativeVert.size(); i++)
        {
            v=TentativeVert[i];
            if ( augmentTentativeDistance(du, adjMatrix(u,v), dist(v)) )
            {
                // note: path(path.size()-3) is the third last vertex in path, thus the u adjacent to vert1, as path(path.size()-2)=vert1, path(path.size()-1)=vert2
                //if ( (vsOption!=0)&&(vsOption->PlanarityCost==true)) planarCost_v = EvalPlanarityCost(v,path[path.size()-3],vert1,vert2,vsOption);
                if ( vsOption.PlanarityCost==true) planarCost_v = tttPath.evalPlanarityCost(v);
                // if dv > du + A(u,v) + + plaCost(u) + planarCost_v, now du must for sure be a finite number
                if ( augmentTentativeDistance(du + plaCost(u) + planar_weight * planarCost_v, adjMatrix(u,v), dist(v)) )
                {
                    if ( vsOption.DiffFaces==true ) IsDiffFace = tttPath.isDiffFaces(v);
                    if ( IsDiffFace )
                    {
                        plaCost(v) = plaCost(u) + planar_weight * planarCost_v;
                        (*prev)(v) = u;
                        dist(v) = du + adjMatrix(u,v) + plaCost(v);
                    };
                };
            };
        };

    };
    
    // For DiffFaces
    if ( vsOption.DiffFaces==true )
    {       
      // now return the original entries of the adjacent matrix 
      SetMatrixEntries(&adjMatrix,ind1,ind2,store_ajdMat);
    }; 
    return new gsGraphPath<T>(vert1,vert2,prev,dist(vert2));
}

} // namespace gismo
