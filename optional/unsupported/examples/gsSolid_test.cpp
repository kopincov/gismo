//#define NDEBUG
#define debug 1
#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsModeling/gsModelingUtils.hpp>



using namespace gismo;


int main(int argc, char *argv[])
{
  std::string filename = "ffoSciencePark2.xml";
  index_t numSamples(30);
  index_t v1(0), v2(1);

  /// Load a gsSolid from a file and test volume segmentation on it.
  gsCmdLine cmd("This test opens a solid file and performs volume segmentation on it.");
  cmd.addPlainString("filename", "File containing solid data (.xml, .gsm)", filename);
  cmd.addInt("e","evertex", "The end vertex of the considered path", v2);
  cmd.addInt("b","bvertex", "The beginning vertex of the considered path", v1);
  cmd.addInt("s","samples", "Number of samples to use for viewing", numSamples);
  
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  gsFileData<> data( filename );
  if (!data.has< gsSolid<> >())
  {
    gsInfo << "Failed to open file\n";
    return 1;
  }
  //Linking error: data.getFirst< gsSolid<> >(sl);
  gsSolid<>::uPtr slh = data.getFirst< gsSolid<> >();
  gsSolid<> &sl = *slh;

  gsInfo<<"\n"<<"==================== gsSolid ===================="<<"\n";
  gsInfo<< "\n" << sl <<"\n";
  gsInfo<<"\n"<<" Trimming loop index of the fist halfedge: " << sl.edge[1]->trimLoopInd();

  gsInfo<<"\n"<<"==================== Face trimming loop information ===================="<<"\n";
  for(unsigned i = 0; i < sl.face.size(); i++)
  {
    gsInfo << i << ":\n";
    int n = sl.face[i]->surf->domain().outer().size();
    for(int j = 0; j < n; j++)
    {
      gsInfo << "(" << sl.face[i]->surf->vertexCoord(0, j).transpose() << "), ";
    }
    gsInfo << "\n";
  }
  gsInfo<<"\n"<<"==================== gsWriteParaview ===================="<<"\n";
  gsWriteParaview( sl , "solidEdgeGraph", 50);

  gsInfo<<"\n"<<"==================== gsGraph ===================="<<"\n";
  //       gsGraph<> gr(true);
  //       gsInfo << gr;
  //       gsGraphPath<>* gp = gr.Dijkstra(0,4);
  //       gsInfo << (*gp) << "\n";




  gsInfo<<"\n"<<"==================== gsVolumeSegment ===================="<<"\n";
  // manually set nonconvex edges
  //    sl.edge[5]->is_convex = false;
  //    sl.edge[5]->mate->is_convex = false;
  gsVolumeSegment<> vs(&sl);
  //vs.vsOption->update("DiffFaces",false);
  gsInfo << vs;
  //gsGraphPath<>* cutLoop = vs.FindCuttingLoop(sl.edge[0]);

  gsCuttingLoop<>* cutLoop = vs.FindCuttingLoop(v1,v2);
  gsInfo << vs;
  gsInfo << *cutLoop<< "\n";


  gsInfo << "\n" << "==================== cuttingSurface ====================" << "\n";
  //    for (int i=0;i<=7;i++)
  //    gsInfo << "Control points of the boundary trimmed surface "<<i<<" : \n" <<
  //            convert2Zero<>(sl.face[i]->surf->getTP()->coefs()) << "\n";
  gsInfo << "Control points of the boundary trimmed surface "<<3<<" : \n" <<
          convert2Zero<>(sl.face[3]->surf->getTP()->coefs()) << "\n";

  gsVolumeSegment<>::splitFaces(sl, *cutLoop, 0.0, 0.1, 20, 0.001, false, 0.5);
  gsWriteParaview( sl , "solidEdgeGraph_withCuttingEdges", 50);
  gsCurveLoop<> * curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cutLoop);
  //gsCurveLoop<> * curveLoop = sl.face[2]->surf->boundaryLoop();
  gsTrimSurface<>* cface;
  gsInfo<<"\n"<< "Resulting curve loop "<< *curveLoop;

  gsInfo << " getVertices: Vertices of the cutting loop: " << "\n";
  gsInfo << cutLoop->getVertices()<<"\n";
  gsInfo << " sample: sample points of the cutting loop: " << "\n";
  gsInfo << cutLoop->sample(3,2);

  int deg = 3;
  gsInfo << "\n" << "Base surface spline degree: " << deg << "\n";
  gsKnotVector<> kv(0,1,0,deg+1);
  gsInterpOption<> intopt(deg);
  cface = vs.cuttingSurface(*cutLoop, *curveLoop, kv,kv,intopt);
  gsInfo << "Control points of the resulting trimmed surface: \n" << convert2Zero<>(cface->getTP()->coefs()) << "\n";

  gsInfo << "\n" << "trimm curves: " << cface->boundaryLoop().curve(0) << "\n";
  for(unsigned i = 0; i < cface->boundaryLoop().curves().size(); i++)
  {
    gsInfo << "\n" << "trimm curves (" << i << ": " << cface->boundaryLoop().curve(i).coefs() << "\n";
  }


  gsMesh<>::uPtr m = cface->toMesh(numSamples);
  gsInfo << "\n" << "Inherited gsMesh: " << *m << "\n";
  gsWriteParaview(*m , "cuttingSurface");

  delete cface;
  delete cutLoop;
  delete curveLoop;

  gsInfo<<"\n";return 0;
}
