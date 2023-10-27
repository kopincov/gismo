/** @file gsHeMesh.h

    @brief Provides gsHeMesh class - half-edge data structure

    This file is part of the G+Smo library. 

    Author(s): D. Kloimstein
*/

#pragma once

#include <set>

#include <gsUtils/gsMesh/gsBoundingBox.h>
#include <gsUtils/gsMesh/gsVertex.h>
#include <gsUtils/gsMesh/gsCell.h>
#include <gsUtils/gsHeVertex.h>
#include <gsUtils/gsHeMeshElement.h>
#include <gsUtils/gsHalfEdge.h>
#include <gsUtils/gsHalfFace.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsBoxSplines/gsBoxSplineBasis.h>

namespace gismo {

/** 
    A half-edge mesh
 */
template <class T >
class gsHeMesh : public gsHeMeshElement<T>
{
private:
    typedef gsHeMeshElement<T> MeshElement;
    typedef typename gsHeMeshElement<T>::scalar_t scalar_t;
    typedef typename gsHeMeshElement<T>::gsHeVertexHandle gsHalfVertexHandle;
    typedef typename gsHeMeshElement<T>::gsHalfEdgeHandle gsHalfEdgeHandle;
    typedef typename gsHeMeshElement<T>::gsHalfFaceHandle gsHalfFaceHandle;
    typedef typename gsHeMeshElement<T>::gsCellHandle gsCellHandle;
    typedef gsVector3d<T>* gsVectorHandle;
    typedef gsMatrix<T> gsMatrixT;
    
private:

    std::vector< gsVector<int>  > halfFace_indices;

    ///Number of vertices
    int numVertices;
    ///Number of halfedged
    int numHalfEdges;
    ///Number of halffaces
    int numHalfFaces;
    ///Number of cells
    int numCells;

    //std::map<int, gsVertexHandle> vertices;
    ///Set of vertices
    std::vector<gsHalfVertexHandle > vertex;
    ///Set of halfedges
    std::vector<gsHalfEdgeHandle >  edge;
    ///Set of halffaces
    std::vector<gsHalfFaceHandle >  face;
    ///Set of cells
    std::vector<gsCellHandle >      cell;
    ///Set of regular faces
    std::vector<gsHalfFaceHandle> regface;
    ///Set of halfedges, which defines the boundary of the mesh
    std::vector<gsHalfEdgeHandle > globalboundary;
    ///Holds the 12 values of the boxspline evaluation
    std::vector<T> phi;
    ///Holds the 12 values of the boxspline derivative evaluation with respect to the vector (1,0)
    std::vector<T> phi1;
    ///Holds the 12 values of the boxspline derivative evaluation with respect to the vector (0,1)
    std::vector<T> phi2;

public:
    /** \brief Constructs a new half-edge mesh
     *  \param i the index of the new half-edge mesh
     */

    gsHeMesh(int const & i = 0) : MeshElement(i)
    {
        numHalfFaces = 0;
        numVertices = 0;
        numCells = 0;
        numHalfEdges = 0;
    }
        
    /// \brief Destructor for the half-edge mesh
    ~gsHeMesh();

    //getter-setter methods
    ///Get the number of vertices
    inline int getnumVertices() {return vertex.size();}
    ///Get the number of halfedges
    inline int getnumHalfEdges() {return edge.size();}
    ///Get the number of halffaces
    inline int getnumHalfFaces() {return face.size();}
    ///Get the number of cells
    inline int getnumCells() {return cell.size();}
    ///Get the set of vertices
    inline std::vector<gsHalfVertexHandle >  getVertices() {return vertex;}
    ///Set vertices
    inline void setVertices(std::vector<gsHalfVertexHandle > ve) {vertex.swap(ve);}
    ///Get vertex at index i
    inline gsHalfVertexHandle getVertexAtIndex(int i) {return vertex[i];}
    ///Get vertex at id i
    inline gsHalfVertexHandle getVertexAtId(int i)
    {
        for ( typename std::vector<gsHalfVertexHandle>::iterator
                      it = vertex.begin(); it!= vertex.end(); ++it)
            {
                gsHalfVertexHandle v = *it;
                if(v->getId() == i)
                {
                    return v;
                    break;
                }
            }
        GISMO_ERROR("Vertex does not exist");
    }

    ///Get the set of halfedges
    inline std::vector<gsHalfEdgeHandle > getHalfEdges() {return edge;}
    ///Set halfedges
    inline void setEdges(std::vector<gsHalfEdgeHandle > ed) {edge.swap(ed);}
    ///Get halfedge at index i
    inline gsHalfEdgeHandle getHalfEdgeAtIndex(int i) {return edge[i];}
    ///Get the set of halffaces
    inline std::vector<gsHalfFaceHandle > getHalfFaces() {return face;}
    ///Get halfface at index i
    inline gsHalfFaceHandle getHalfFaceAtIndex(int i) {return face[i];}
    ///Set halffaces
    inline void setFaces(std::vector<gsHalfFaceHandle > fa) {face.swap(fa);}
    ///Get the set of cells
    inline std::vector<gsCellHandle > getCells() {return cell;}
    ///Get cell at index i
    inline gsCellHandle getCellAtIndex(int i) {return cell[i];}
    ///Set cells
    inline void setCells(std::vector<gsCellHandle > ce) {cell.swap(ce);}
    ///Get set of regular halffaces
    inline std::vector<gsHalfFaceHandle > getRegularHalfFaces() {return regface;}
    ///Get regular halfface at index i
    inline gsHalfFaceHandle getRegularHalfFaceAtIndex(int i) {return regface[i];}
    ///Set regular halffaces
    inline void setRegularFaces(std::vector<gsHalfFaceHandle > fa) {regface.swap(fa);}
    ///Get the set of halfedges, which defines the boundary of the mesh
    inline std::vector<gsHalfEdgeHandle > getGlobalBoundary() {return globalboundary;}
    ///Set halfedges, which defines the boundary of the mesh
    inline void setGlobalBoundary(std::vector<gsHalfEdgeHandle > gb) {globalboundary.swap(gb);}
    ///Get the boxspline evaluation values
    inline std::vector<T> getPhi() {return phi;}
    ///Get the boxspline derivative evaluation values with respect to the vector (1,0)
    inline std::vector<T> getPhi1() {return phi1;}
    ///Get the boxspline derivative evaluation values with respect to the vector (0,1)
    inline std::vector<T> getPhi2() {return phi2;}
    ///Set the boxspline evaluation values
    inline void setPhi(std::vector<T> p) {phi.swap(p);}
    ///Set the boxspline derivative evaluation values with respect to the vector (1,0)
    inline void setPhi1(std::vector<T> p1) {phi1.swap(p1);}
    ///Set the boxspline derivative evaluation values with respect to the vector (0,1)
    inline void setPhi2(std::vector<T> p2) {phi2.swap(p2);}
    ///Solve the partial differential equation
    void solveEquation(T x0=0, T x1=0, T x2=0, T x3=0, T y0=0, T y1=0, T y2=0, T y3=0, T x1y1=0, T x2y1=0, T x3y1=0, T x1y2=0, T x2y2=0, T x3y2=0, T x1y3=0, T x2y3=0, T x3y3=0);
    ///Calculate the LU-Decomposition from matrix A and solve the equation Ax=b
    gsMatrix<T> LUDecomposition(gsMatrix<T> A, gsMatrix<T> b);
    ///Interpolate at given points y
    gsMatrix<T> interpolate(std::vector<T> y, int option);

    ///Calculate the Distance of two vertices with the euclidean distance
    T VertexDistance(gsHalfVertexHandle v1, gsHalfVertexHandle v2);

    gsHalfFaceHandle getHalfFaceByBoundary(std::vector<gsHalfVertexHandle> boundary)
    {
        gsHalfFaceHandle emptyface;
        for ( typename std::vector<gsHalfFaceHandle>::iterator
                  it = face.begin(); it!= face.end(); ++it)
        {
            std::vector<gsHalfEdgeHandle> heh = (*it)->getBoundary();
            if((heh[0]->getSource()==boundary[0]||heh[0]->getSource()==boundary[1]||heh[0]->getSource()==boundary[2]) &&
               (heh[1]->getSource()==boundary[0]||heh[1]->getSource()==boundary[1]||heh[1]->getSource()==boundary[2]) &&
               (heh[2]->getSource()==boundary[0]||heh[2]->getSource()==boundary[1]||heh[2]->getSource()==boundary[2]))
            {
                return *it;
            }
        }
        return emptyface;
    }

    gsHalfFaceHandle getHalfFaceByBoundary(gsHalfVertexHandle v1, gsHalfVertexHandle v2, gsHalfVertexHandle v3)
    {
        gsHalfFaceHandle emptyface;
        for ( typename std::vector<gsHalfFaceHandle>::iterator
                  it = face.begin(); it!= face.end(); ++it)
        {
            std::vector<gsHalfEdgeHandle> heh = (*it)->getBoundary();
            if((heh[0]->getSource()==v1||heh[0]->getSource()==v2||heh[0]->getSource()==v3) &&
               (heh[1]->getSource()==v1||heh[1]->getSource()==v2||heh[1]->getSource()==v3) &&
               (heh[2]->getSource()==v1||heh[2]->getSource()==v2||heh[2]->getSource()==v3))
            {
                return *it;
            }
        }
        return emptyface;
    }
    /**
     * Initialize properties with methods setHeMate(), setFaceClassifications(),
     * setValence(), setRegularity(int i) and setBasisFactor() for the half-edge mesh
     */
    void initialize();
    ///Add a vertex to the half-edge mesh
    void addVertex(gsHalfVertexHandle v);

    ///Add vertex with the coordinates x,y,z to the half-edge mesh and return this vertex
    gsHalfVertexHandle addVertex(scalar_t const& x, scalar_t const& y, scalar_t const& z=0);
    ///Add a halfedge to the half-edge mesh
    void addHalfEdge(gsHalfEdgeHandle he);

    ///Add halfface with the vertices V as boundary to the half-edge mesh and return this halfface
    gsHalfFaceHandle addHalfFace(std::vector<gsHalfVertexHandle> V);
    
    // todo:  adds a half-face given indices of vertices
    //gsHalfFaceHandle addHalfFace(const std::vector<index_t> &  V);

    /*std::ostream &print(std::ostream &os) const
     {
         os<<"gsHeMesh with "<<numVertices<<" vertices and "<<numHalfFaces<<" half-faces.\n";
         for ( typename std::vector<gsHalfVertexHandle>::const_iterator
                   it = vertex.begin(); it!= vertex.end(); ++it)
             os<< **it ;
         return os;
     } */

    //------------------------------------------------------------------------------------------------------
    // Set up default values for the spline information of each trimming curve attached to each HE	
 //   void setDefaultTrimmingLine();
    
    //------------------------------------------------------------------------------------------------------
    ///Assigning mates for each halfedge
    void setHeMate();	    
    
    //------------------------------------------------------------------------------------------------------
    //generating full control points from four corners of a patch 
  //  gsMatrixT* simpleFaceCP(std::vector< gsVector3d<T>* > corner, gsKnotVector<T>* kv1, gsKnotVector<T>* kv2);

    ///Get the number of halfedges, which go to or go out of the vertex s
    int getValence(gsHalfVertexHandle s);
    ///Get the number of halfedges, which go to or go out of the vertex s for a limited area by indices
    int getValence(gsHalfVertexHandle s, std::list<int> indices);
    ///Set the number of halfedges, which go to or go out of the vertex s
    void setValence();

  /*  void getMate(gsHalfVertexHandle const & source, gsHalfVertexHandle const & target)
    {
        gsHalfEdgeHandle heh = findHalfEdge(source,target);
        heh = heh->mate;
        std::cout << heh->source->x();
        std::cout << heh->source->y()<<" -> ";
        std::cout << heh->target->x();
        std::cout << heh->target->y()<<"\n";
    } */

    ///Return the mass-matrix of the evaluated basisfunctions
    gsMatrix<T> getMassMatrix(std::vector<gsHalfFaceHandle> refmeshfaces, std::vector<gsHalfVertexHandle> refmeshvertices);
    ///Return the stiffness-matrix of the evaluated basisfunctions with derivative r
    gsMatrix<T> getStiffnessMatrix(std::vector<gsHalfFaceHandle> refmeshfaces, std::vector<gsHalfVertexHandle> refmeshvertices);
    ///Return the load vector of the evaluated basisfunctions with respect to the functionvalue
    gsMatrix<T> getLoadVector(std::vector<gsHalfFaceHandle> refmeshfaces, std::vector<gsHalfVertexHandle> refmeshvertices, T functionvalue);
    ///Return the 1-Ring neighborhood of the halfface fa
    std::list<int> get1Ring(gsHalfFaceHandle fa);
    ///Return the 1-Ring neighborhood of the vertex ve
    std::list<int> get1Ring(gsHalfVertexHandle ve);
    ///Return the 2-Ring neighborhood of the halfface fa
    std::list<int> get2Ring(gsHalfFaceHandle fa);
    ///Assign classifications 0=domain, 1=extended, 2=outer to all halffaces
    void setFaceClassifications();
    ///Return the set of halffaces with the classification as parameter
    std::vector<gsHalfFaceHandle> getFacesByClassification(int classification);
    /**
     * Set regularity for each halfface. A halfface is regular, if the valence of each
     * bounding vertex is equal to the parameter check
     */
    void setRegularity(int check);
    ///Set the number of basis functions
    void setBasisFactor();

    ///Set up the domain for calculations
    void buildReferenceConfiguration();

    ///Return type of vertex with type0 = inner vertex, type1=vertex with valence 3 and type2=vertex with valence 4
    int getVertexType(gsHalfVertexHandle v, std::list<int> indices);

    ///Transform the input into the referenceconfiguration
    std::vector<T> transformPoint(T bary1, T bary2, T bary3, gsHalfVertexHandle v1, gsHalfVertexHandle v2, gsHalfVertexHandle v3, gsHalfVertexHandle veval, std::vector<gsHalfFaceHandle> f, std::list<int> indices);

    ///Transform all points in the 1Ring neighborhood into the referenceconfiguration and evaluate them. Return the evaluation points.
    gsMatrix<T> evaluate1Ring(T bary1, T bary2, T bary3, gsHalfVertexHandle v1, gsHalfVertexHandle v2, gsHalfVertexHandle v3, gsHalfFaceHandle hf, std::vector<gsHalfFaceHandle> refmeshfaces);

    ///Transform all points in the 1Ring neighborhood into the referenceconfiguration and evaluate them. Return all basis functions of the evaluation.
    gsMatrix<T> basisfunctions(T bary1, T bary2, T bary3, gsHalfVertexHandle v1, gsHalfVertexHandle v2, gsHalfVertexHandle v3, gsHalfFaceHandle hf, std::vector<gsHalfFaceHandle> refmeshfaces);

    ///Store the boxspline evaluation values and their derivative values into phi, phi1 and phi2. Part of method initialize.
    void setBoxsplineTable();

    void getAllEdgesPrinted();

    void Testprint();

    //gsHalfEdgeHandle changeOrientation(gsHalfEdgeHandle heh);

    //std::vector<gsHalfEdgeHandle> changeOrientation(std::vector<gsHalfEdgeHandle> heh);

    //bool checkRightBoundaryOrientation();

};
 
} // namespace gismo

//#ifndef GISMO_BUILD_LIB
//#include GISMO_HPP_HEADER(gsHeMesh.hpp)
//#endif
