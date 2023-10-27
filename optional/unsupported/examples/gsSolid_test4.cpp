//#define NDEBUG
#define debug 1
#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsModeling/gsModelingUtils.hpp>



using namespace gismo;


void makeLShapedSolid(gsSolidHalfFace<> *&frontFace, gsSolidHalfFace<> *&backFace, gsSolid<> &sl)
{
  // add coords to gsSolidHeVertex, not yet pointers to HEs
  real_t eps = 0.5;
  real_t y0,y1,y2,z0,z1,z2;
  y0 = 0.;
  y1 = 1.;
  y2 = 2.;
  z0 = 0.; 
  z1 = 1.;
  z2 = 2.;   
  sl.addHeVertex(0.,  y0,  z0);
  sl.addHeVertex(0.,  y2,  z0);
  sl.addHeVertex(1.,  y0,  z0);
  sl.addHeVertex(1.,  y2,  z0); // vertex 3
  sl.addHeVertex(0.,  y1,  z1); 
  sl.addHeVertex(0.,  y2-eps,  z1); // vertex 5
  sl.addHeVertex(1.,  y1,  z1); 
  sl.addHeVertex(1.,  y2-eps,  z1); // vertex 7
  sl.addHeVertex(0.,  y0,  z2); 
  sl.addHeVertex(0.,  y1,  z2-eps); // vertex 9
  sl.addHeVertex(1.,  y0,  z2);
  sl.addHeVertex(1.,  y1,  z2-eps); // vertex 11
   

  // add faces with 4 vertices
  // bottom
  sl.addFace_4Vertices(sl.vertex[0], sl.vertex[1], sl.vertex[3], sl.vertex[2]);
  // face with vertices 2,3,7,6,11,10
  std::vector<gsSolidHeVertex<>* > frontFaceVertices;
  frontFaceVertices.push_back(sl.vertex[2]);
  frontFaceVertices.push_back(sl.vertex[3]);
  frontFaceVertices.push_back(sl.vertex[7]);
  frontFaceVertices.push_back(sl.vertex[6]);
  frontFaceVertices.push_back(sl.vertex[11]);
  frontFaceVertices.push_back(sl.vertex[10]);
  frontFace = sl.addFace_PlanarPolygon(frontFaceVertices);
  // face with vertices 0,8,9,4,5,1
  std::vector<gsSolidHeVertex<>* > backFaceVertices;
  backFaceVertices.push_back(sl.vertex[0]);
  backFaceVertices.push_back(sl.vertex[8]);
  backFaceVertices.push_back(sl.vertex[9]);
  backFaceVertices.push_back(sl.vertex[4]);
  backFaceVertices.push_back(sl.vertex[5]);
  backFaceVertices.push_back(sl.vertex[1]);
  backFace = sl.addFace_PlanarPolygon(backFaceVertices);
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
  
  // gsInfo << "\nOriginal front face:\n" << *frontFace << "\n";
}

int main(int argc, char * argv[])
{
  ///---------------------------------------------------------------------
  /// Modeling an L shape using gsSolid
  gsSolidHalfFace<> *frontFace, *backFace;
  gsSolid<> sl;
  makeLShapedSolid(frontFace, backFace, sl);
  real_t weightCurveReg(0), weightSurfReg(0);
  index_t numSamples = 30;
  bool harmonic = false;
  
  gsCmdLine cmd("This test computes a cutting surface for an in-built shape.");
  cmd.addInt("s","samples", "Number of samples to use for viewing", numSamples);
  cmd.addReal("c","weight-curve-reg", 
              "Weighting for curve regularity (vs curve approximation)", weightCurveReg);
  cmd.addReal("u","weight-surf-reg",
              "Weighting for surface regularity (vs surface approximation)", weightSurfReg);
  cmd.addSwitch("harmonic",
                "Use harmonic functions instead of mean value interpolation", harmonic);
  
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  ///---------------------------------------------------------------------
  /// print out the gsSolid and export its data for visualization
  if (0)
  {
    gsInfo << "\n" << sl << "\n";

    // Output a paraview file and open it
    gsWriteParaview(sl , "SolidEdgeGraph");

    // Call paraview on exit
    char cmd2[100];
    strcpy(cmd2, "paraview SolidEdgeGraph.vtp\0");
    strcat(cmd2, " &");
    //system(cmd2);
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

  gsInfo << "\nOriginal solid:\n" << sl << "\n";

  gsInfo << "\noriginal front face domain defined as follows:";

  std::vector< gsCurve<> *> curves = frontFace->surf->domain().loop(0).curves();

  for (size_t i = 0; i < curves.size(); i++)
  {
    gsInfo << "\ncurve " << i << ":\n" << curves[i]->coefs();
  }

  // test the creation of an auxiliary edge and splitting an L-shaped face along it.

  // request a spline (in the domain) from frontface's vertex 0 to vertex 3
  // (the vertices with global id 2 and 6).
  gsBSpline<> * frontAuxE = new gsBSpline<>(frontFace->surf->cuttingCurve(0, 3));

  // split the face along the new spline
  gsSolidHalfFace<> * frontFace2 = sl.splitFace(frontFace, sl.vertex[2], sl.vertex[6], frontAuxE);

  gsInfo << "\n---Front face modified---\n";
  gsInfo << "\nNew front face:\n" << *frontFace2 << "\n";
  gsInfo << "\nModified front face:\n" << *frontFace << "\n\n";

  // output some info about the trimming loops
  gsInfo << "\nnew face domain defined as follows:";
  curves = frontFace2->surf->domain().loop(0).curves();

  for (size_t i = 0; i < curves.size(); i++)
  {
    gsInfo << "\ncurve " << i << ":\n" << curves[i]->coefs();
  }

  gsInfo << "\nmodified front face domain defined as follows:";
  curves = frontFace->surf->domain().loop(0).curves();

  for (size_t i = 0; i < curves.size(); i++)
  {
    gsInfo << "\ncurve " << i << ":\n" << curves[i]->coefs();
  }

  // do the same for the back face (split across vertices 0 and 4)
  gsBSpline<> * backAuxE = new gsBSpline<>(backFace->surf->cuttingCurve(0, 3));
  gsSolidHalfFace<> * backFace2 = sl.splitFace(backFace, sl.vertex[0], sl.vertex[4], backAuxE);

  gsInfo << "\n---Back face modified---\n";
  gsInfo << "\nNew back face:\n" << *backFace2 << "\n";
  gsInfo << "\nModified back face:\n" << *backFace << "\n\n";

  // output some info about the trimming loops
  gsInfo << "\nnew back face domain defined as follows:";
  curves = backFace2->surf->domain().loop(0).curves();

  for (size_t i = 0; i < curves.size(); i++)
  {
    gsInfo << "\ncurve " << i << ":\n" << curves[i]->coefs();
  }

  gsInfo << "\nmodified back face domain defined as follows:";
  curves = backFace->surf->domain().loop(0).curves();

  for (size_t i = 0; i < curves.size(); i++)
  {
    gsInfo << "\ncurve " << i << ":\n" << curves[i]->coefs();
  }

  gsInfo << "\nFinal solid:\n" << sl << "\n";

  
  //--------------------------------------------------------------------------------------------------
  // recreate the L-shape
  
  gsInfo << "\nResetting the solid.\n";
  
  gsSolid<> sl2;
  makeLShapedSolid(frontFace, backFace, sl2);
  // manually set nonconvex edges: edge 5 for Lshape
  sl2.edge[5]->is_convex = false;
  sl2.edge[5]->mate->is_convex = false;
  gsWriteParaview(sl2 , "solidEdgeGraph",50);


  // find a cutting loop
  gsVolumeSegment<> vs(&sl2);
  //gsCuttingLoop<> * cuttingLoop = vs.FindCuttingLoop(4, 6);
  gsCuttingLoop<> * cuttingLoop = vs.FindCuttingLoop(8, 10);
  gsInfo << "\n---------------------------\nCutting loop found: " << *cuttingLoop << "\n";

  // test automatic face splitting
  gsVolumeSegment<>::splitFaces(sl2, *cuttingLoop, weightCurveReg, 1.0, 20, 0.0001, harmonic, 0.5);
  
  gsInfo << "\nFinal solid:\n" << sl2 << "\n";
  
  // turn the cutting loop into a curve loop
  gsCurveLoop<> * curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl2, *cuttingLoop);
  gsInfo << "\n---Curve loop found---\n";
  curves = curveLoop->curves();

  // output some info about the curve loop
  gsInfo << "\nCurve loop info:\n";

  for (size_t i = 0; i < curves.size(); i++)
  {
    gsInfo << "\n==" << i << "==\n" << *curves[i];
  }

  // Call paraview on exit
  char cmd2[100];
  strcpy(cmd2, "paraview solidEdgeGraph.vtp\0");
  strcat(cmd2, " &");
  //system(cmd2);

  gsInfo << "\n" << "==================== gsCuttingLoop ====================" << "\n";

  gsInfo << " sample points in parameter domains: " << "\n";
  gsInfo << curveLoop->sample(2, 1) << "\n";

  gsInfo << " getVertices: Vertices of the cutting loop: " << "\n";
  gsInfo << cuttingLoop->getVertices() << "\n";
  gsInfo << " sample: sample points of the cutting loop: " << "\n";
  gsInfo << cuttingLoop->sample(2, 2) << "\n";
  gsInfo << " sampleNormal: sample Normals of the trimmed surface long the cutting loop: " << "\n";
  gsInfo << cuttingLoop->sampleNormal(3) << "\n";

  delete curveLoop;

  gsInfo << "\n" << "==================== cuttingSurface ====================" << "\n";
  const gsCurveLoop<>& cloop = sl2.face[0]->surf->boundaryLoop();
  
  int deg = 3;
  gsInfo << "\n" << "Base surface spline degree: " << deg << "\n";
  gsKnotVector<> kv(0,1,0,deg+1);
  gsInterpOption<> intopt(deg);
  intopt.wReg = weightSurfReg;
  gsTrimSurface<>* cface = vs.cuttingSurface(*cuttingLoop, cloop, kv,kv,intopt);
  gsInfo << "Control points of the resulting trimmed surface: \n" << convert2Zero<>(cface->getTP()->coefs()) << "\n";

  gsInfo << "\n" << "trimm curves: " << cface->boundaryLoop().curve(0) << "\n";
  gsInfo << "\n" << "trimm curves: " << cface->boundaryLoop().curve(0).coefs() << "\n";
  gsInfo << "\n" << "trimm curves: " << cface->boundaryLoop().curve(1).coefs() << "\n";
  gsInfo << "\n" << "trimm curves: " << cface->boundaryLoop().curve(2).coefs() << "\n";
  gsInfo << "\n" << "trimm curves: " << cface->boundaryLoop().curve(3).coefs() << "\n";

  gsInfo << "\n" << "==================== addFaceWithMate ====================" << "\n";
  std::vector<unsigned> graphPath = cuttingLoop->computePath();
  std::vector<gsSolidHeVertex<>*> cuttingLoopVerts;
  for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl2.vertex[graphPath[i]]);
  sl2.addFaceWithMate(cuttingLoopVerts, cface);
  gsInfo << sl2 << "\n\n";
  gsInfo << "\n" << "==================== begin list of edges with volume ====================" << "\n";
  for(std::vector<gsSolidHalfEdge<>*>::iterator iter = sl2.edge.begin(); iter != sl2.edge.end(); iter++)
  {
    if(!(*iter)->face->vol) gsInfo << "(vol NULL, face " << (*iter)->face->getId() << ") ";
    else gsInfo << "(vol " << (*iter)->face->vol->getId() << ", face " << (*iter)->face->getId() << ") ";
    gsInfo << (*iter)->source->getId() << " -- " <<
            (*iter)->next->source->getId() << "\n";
  }

  gsInfo << "\n" << "==================== end list of edges with volume ====================" << "\n";

  gsMesh<>::uPtr m = cface->toMesh(numSamples);
  gsInfo << "\n" << "Inherited gsMesh: " << *m << "\n";
  gsWriteParaview(*m , "cuttingSurface");

  gsTrimSurface<> *tSurf=vs.cuttingSurface(*cuttingLoop, cloop, kv,kv,intopt);
  gsInfo <<"\nResulting based surface: \n"<< *tSurf->getTP();
  gsInfo << "\nControl points of the resulting based surface: \n" << convert2Zero<>(tSurf->getTP()->coefs()) << "\n";

  // plot each (sub) volumes of the solid
  gsWriteParaview(sl2 , "solidEdgeGraph1", 50, 1);
  gsWriteParaview(sl2 , "solidEdgeGraph0", 50, 0);
  gsInfo << "\n test test test \n";

  delete cuttingLoop;
  delete tSurf;
  return 0;
}

