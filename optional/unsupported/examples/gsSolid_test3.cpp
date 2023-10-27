//#define NDEBUG
#define debug 1
#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>


using namespace gismo;


/// construct a 2d L-shaped face for testing gsSolid.
/**
 * \param sl the solid to add the face to.
 * \param vertices a std::vector of vertices making the face (anticlockwise from
 * outside by default) starting with the bottom-left corner of the L shape.
 * \param clockwise whether to reverse the order of the vertices in the domain
 * (to get the orientation right).
 */

void makeLFace(gsSolid<> &sl, std::vector<gsSolidHeVertex<>* > &vertices, bool clockwise)
{
  int curveDeg = 1;
  int patchDeg1 = 1, patchDeg2 = 1;

  // set up a 7x2 matrix with the (transposed) vectors of the L shape in the domain.
  // (last vector is equal to the first one)
  int numVertices = 6;
  gsMatrix<> DomCor = gsMatrix<>(numVertices + 1, 2);
  if(!clockwise)
  {
    DomCor << 0, 0, 2, 0, 2, 1, 1, 1, 1, 2, 0, 2, 0, 0;
  }
  else
  {
    DomCor << 0, 0, 0, 2, 1, 2, 1, 1, 2, 1, 2, 0, 0, 0;
  }

  // construct a loop of B-splines that are straight lines in the domain.
  // this code is nabbed from gsTrimSurface.h
  
  // Construct the Bezier knot vector for all trimming curves
  gsKnotVector<> kv(0, 1, 0, curveDeg+1);
  unsigned int ntcp = kv.size() - kv.degree() - 1;
  
  // Initialize a trimming curve loop
  gsCurveLoop<> * tloop = new gsCurveLoop<>();
  
  for (int ic = 0; ic < numVertices; ic++)
  {
    // Define a spline curve from the current vertex to the next one
    gsMatrix<> tcp(ntcp, 2);

    for (unsigned int i = 0; i <= ntcp - 1; i++)
    {
      for (unsigned int xi = 0; xi <= 1; xi++)
      {
        tcp(i, xi) = DomCor(ic, xi) + (real_t)(i) / ((real_t)(ntcp) - 1) * (DomCor(ic + 1, xi) - DomCor(ic, xi));
      }
    }

    gsBSpline<> * tcurve = new gsBSpline<>(0, 1, 0, curveDeg, give(tcp));
    tloop->insertCurve(tcurve);
  }
  // Then construct a planar domain with only an outer loop *tloop*
  gsPlanarDomain<> * domain1= new gsPlanarDomain<>(tloop);
    
  // Define the base NURBS surface
  gsKnotVector<> KV1 = gsKnotVector<>(0, 2, 0, patchDeg1+1);
  gsKnotVector<> KV2 = gsKnotVector<>(0, 2, 0, patchDeg2+1);
  
  // construct a B-Spline basis. this code nabbed from gsTensorBSpline
  gsBSplineBasis<> * Bu= new gsBSplineBasis<>(KV1);
  gsBSplineBasis<> * Bv= new gsBSplineBasis<>(KV2);
  typedef gsTensorBSplineBasis<2> Basis; //dimension==2
  Basis *tbasis = new Basis(Bu,Bv);
  
  // assuming the L-shape's vertices are all on one plane, we can easily
  // figure out the control points for a linear spline
  gsMatrix<> pcp(4, 3);
  for(int xi = 0; xi < 3; xi++) // loop over coordinates
  {
    pcp(0, xi) = vertices[0]->coords(xi);
    pcp(1, xi) = vertices[1]->coords(xi);
    pcp(2, xi) = vertices[5]->coords(xi);
    pcp(3, xi) = vertices[5]->coords(xi) + vertices[1]->coords(xi) - vertices[0]->coords(xi);
  }
  
  gsTensorBSpline<2>::Ptr tp1(new gsTensorBSpline<2,real_t>(*tbasis, give(pcp)));
  delete tbasis;
  
  gsTrimSurface<> * ts = new gsTrimSurface<>(tp1, domain1);
  sl.addFace(vertices, ts);
}

int main(int argc, char * argv[])
{
  ///---------------------------------------------------------------------
  /// Modeling an L shape using gsSolid

  gsSolid<> sl;

  // add coords to gsSolidHeVertex, not yet pointers to HEs
  sl.addHeVertex(0.,  0.,  0.);
  sl.addHeVertex(0.,  2.,  0.);
  sl.addHeVertex(1.,  0.,  0.);
  sl.addHeVertex(1.,  2.,  0.);
  sl.addHeVertex(0.,  1.,  1.);
  sl.addHeVertex(0.,  2.,  1.);
  sl.addHeVertex(1.,  1.,  1.);
  sl.addHeVertex(1.,  2.,  1.);
  sl.addHeVertex(0.,  0.,  2.);
  sl.addHeVertex(0.,  1.,  2.);
  sl.addHeVertex(1.,  0.,  2.);
  sl.addHeVertex(1.,  1.,  2.);

  // add faces with 4 vertices
  // bottom
  sl.addFace_4Vertices(sl.vertex[0], sl.vertex[1], sl.vertex[3], sl.vertex[2]);
  // face with vertices 2,3,7,6,11,10
  std::vector<gsSolidHeVertex<>* > frontFace;
  frontFace.push_back(sl.vertex[2]);
  frontFace.push_back(sl.vertex[3]);
  frontFace.push_back(sl.vertex[7]);
  frontFace.push_back(sl.vertex[6]);
  frontFace.push_back(sl.vertex[11]);
  frontFace.push_back(sl.vertex[10]);
  makeLFace(sl, frontFace, false);
  // face with vertices 0,8,9,4,5,1
  std::vector<gsSolidHeVertex<>* > backFace;
  backFace.push_back(sl.vertex[0]);
  backFace.push_back(sl.vertex[8]);
  backFace.push_back(sl.vertex[9]);
  backFace.push_back(sl.vertex[4]);
  backFace.push_back(sl.vertex[5]);
  backFace.push_back(sl.vertex[1]);
  makeLFace(sl, backFace, true);
  // sides with 4 edges
  sl.addFace_4Vertices(sl.vertex[0], sl.vertex[2], sl.vertex[10], sl.vertex[8]);
  sl.addFace_4Vertices(sl.vertex[1], sl.vertex[5], sl.vertex[7], sl.vertex[3]);
  sl.addFace_4Vertices(sl.vertex[4], sl.vertex[9], sl.vertex[11], sl.vertex[6]);
  // top faces
  sl.addFace_4Vertices(sl.vertex[4], sl.vertex[6], sl.vertex[7], sl.vertex[5]);
  sl.addFace_4Vertices(sl.vertex[8], sl.vertex[10], sl.vertex[11], sl.vertex[9]);


  // Now: only members 'mate' of halfedge is missing, we will set them up now
  sl.setHeMate();
  // Add volume
  sl.addVolume(sl.face);

  ///---------------------------------------------------------------------
  /// print out the gsSolid and export its data for visualization
  if (0)
  {
    gsInfo << "\n" << sl << "\n";

    // Output a paraview file and open it
    gsWriteParaview(sl , "SolidEdgeGraph");

    // Call paraview on exit
    char cmd[100];
    strcpy(cmd, "paraview SolidEdgeGraph.vtp\0");
    strcat(cmd, " &");
    //system(cmd);
  };

  ///----------------------------------------------------------------------------------------------------------------
  /// tests for gsGraph
  typedef gsGraph<real_t> gsGraph;

  typedef gsGraphPath<real_t> gsGraphPath;

  if (0)
  {
    gsGraph gr(true);
    gsInfo << gr;
    gsGraphPath * gp = gr.Dijkstra(0, 4);
    gsInfo << (*gp) << "\n";
  };

  ///----------------------------------------------------------------------------------------------------------------
  /// tests for gsVolumeSegment
  if (1)
  {
    typedef gsVolumeSegment<> gsVolumeSegment;

    //sl.edge[10]->is_convex = false;
    //sl.edge[10]->mate->is_convex = false;
    gsSolidHalfEdge<> * convEdge = sl.getVertexFromID(4)->getHalfEdge( sl.getVertexFromID(6) );
    convEdge->is_convex = false;
    convEdge->mate->is_convex = false;
    gsVolumeSegment vs(&sl);
    //vs.vsOption->update("DiffFaces",false);
    gsInfo << vs;
    //gsGraphPath<>* cutLoop = vs.FindCuttingLoop(sl.edge[0]);
    int v1 = 0;
    int v2 = 1;

// Use arguments instead.. tests take part on automated runs, so they should not have user input.
    gsInfo << "\n" << "Vertex 1: ";
//    cin >> v1;
    gsInfo << "\n" << "Vertex 2: ";
//    cin >> v2;

    gsGraphPath * cutLoop = vs.FindCuttingLoop(v1, v2);
    gsInfo << vs;
    gsInfo << *cutLoop << "\n";
    delete cutLoop;
  };

  if (1)
  {
    gsInfo << "\n" << sl << "\n";

    // Output a paraview file and open it
    //gsWriteParaview( sl , "solidEdgeGraph", 50, 2);
    gsWriteParaview(sl , "solidEdgeGraph", 50);

    // Call paraview on exit
    char cmd[100];
    strcpy(cmd, "paraview solidEdgeGraph.vtp\0");
    strcat(cmd, " &");
    //system(cmd);
  };



  return 0;
}





