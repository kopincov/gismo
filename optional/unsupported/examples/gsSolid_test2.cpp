
// EXAMPLES
//#define NDEBUG
#define debug 10
#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;


int main(int argc, char *argv[])
{
  std::string filename =
      //"ffoTriDoublePyramid.xml";
      //"ffoTriBiPyramid.xml";
      //"ffoVaseHalfed.xml";
      //"ffoSciencePark2.xml";
      "ffoSciencePark2_lowTail.xml";

  index_t numSamples(50);
  index_t v1 = 0;
  index_t v2 = 1;

  /// Load a gsSolid from a file and test volume segmentation on it.
  gsCmdLine cmd("This test opens a solid file and performs volume segmentation on it.");
  cmd.addPlainString("filename", "File containing the input mesh", filename);
  cmd.addInt("s","samples", "Number of samples to use for viewing", numSamples);
  cmd.addInt("b","bvertex", "The beginning vertex of the considered path", v1);
  cmd.addInt("e","evertex", "The end vertex of the considered path", v2);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  filename = gsFileManager::find(filename);

  gsFileData<>  data( filename );
  if (!data.has< gsSolid<> >())
  {
    gsInfo << "Failed to open file\n";
    return 1;
  }
  std::string basename = gsFileManager::getFilename(filename);
  //Linking error: data.getFirst< gsSolid<> >(sl);
  gsSolid<>::uPtr slh = data.getFirst< gsSolid<> >();
  gsSolid<> &sl = *slh;

  gsInfo << "\n" << sl <<"\n";

  std::vector<int> ncEdgeV1; // ncEdgeV1[i] and ncEdgeV2[i] are the two vertices of one nonconvex edge
  std::vector<int> ncEdgeV2;
  unsigned edgePointn=200;
  int deg = 3; // spline degree of segmenting surfaces
  int vol;
  gsSolidHalfEdge<> * he;
  gsCurveLoop<> * curveLoop;
  std::vector< gsSolidHalfEdge<>* > nce;
  int nceNo;
  real_t linewidth(0.05);
  gsVolumeSegment<> vs;
  gsCuttingLoop<>* cuttingLoop;
  gsTrimSurface<>* cface;
  gsKnotVector<> kv(0,1,0,deg+1);
  gsInterpOption<> intopt(deg);
  gsMesh<>::uPtr m;
  std::vector<unsigned> graphPath;
  std::vector<gsSolidHeVertex<>*> cuttingLoopVerts;
  //------------------------------------------------------------------------------
  // Triangular Biparamid
  if (basename=="ffoTriBiPyramid.xml")
  {
    sl.checkStructure();
    gsWriteParaview<real_t>( sl , "tbpEG1", edgePointn, 0, linewidth/3);
    ncEdgeV1.clear();
    ncEdgeV1.push_back(2);ncEdgeV2.push_back(3);
    ncEdgeV1.push_back(1);ncEdgeV2.push_back(3);
    ncEdgeV1.push_back(0);ncEdgeV2.push_back(3);
    nce = sl.detectNonConvexEdges(ncEdgeV1,ncEdgeV2);

    sl.checkStructure();
    // ----------- step 1 --------------------
    nceNo = 0;
    vol = sl.vertex[0]->hed->face->vol->getId();
    he = 0;
    for (size_t i=0;i<sl.edge.size();i++)
    {
        if ( (sl.edge)[i]->face->vol->getId()==vol && (sl.edge)[i] == nce[nceNo]) {he = (sl.edge)[i];break;}
    }
    sl.checkStructure();
    assert(he!=0); // Halfedge not found
    sl.handleImpedingEdges(he);
    sl.checkStructure();
    gsVolumeSegment<> vst(&sl);
    vst.volSegOp().costAuxEdge = 2;
    vst.volSegOp().PlanarityCostWeight = 100;
    sl.checkStructure();
    cuttingLoop = vst.FindCuttingLoop(he);
    gsInfo<<"\n Segmenting loop found: "<< *cuttingLoop;
    sl.checkStructure();
    gsVolumeSegment<>::splitFaces(sl, *cuttingLoop, 0.1, 0.1, 20, 0.001, false, 0.5);
    //
    graphPath = cuttingLoop->computePath();
    gsWriteParaview( sl , "tbpSL1", edgePointn, 0, linewidth, gsVector3d<>(0,0,0), 0, 10, 5, graphPath);
    //
    curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cuttingLoop);
    cface = vst.cuttingSurface(*cuttingLoop, *curveLoop, kv,kv,intopt);
    delete curveLoop;
    m = cface->toMesh(numSamples);
    gsWriteParaview(*m , "tbpCS1");

    graphPath = cuttingLoop->computePath();
    cuttingLoopVerts.clear();
    for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl.vertex[graphPath[i]]);
    sl.addFaceWithMate(cuttingLoopVerts, cface);

    gsWriteParaview(sl , "tbpEG11", edgePointn, vol, linewidth);
    gsWriteParaview(sl , "tbpEG12", edgePointn, sl.numVolumes-1, linewidth);
    // -- cut --
    delete cuttingLoop;
  }
  //------------------------------------------------------------------------------
  // Triangular Double Paramid
  else if (basename=="ffoTriDoublePyramid.xml")
  {
    gsWriteParaview<real_t>( sl , "tdpEG1", edgePointn, 0, linewidth/3);
    ncEdgeV1.clear();
    ncEdgeV1.push_back(2);ncEdgeV2.push_back(3);
    ncEdgeV1.push_back(1);ncEdgeV2.push_back(3);
    ncEdgeV1.push_back(0);ncEdgeV2.push_back(3);
    nce = sl.detectNonConvexEdges(ncEdgeV1,ncEdgeV2);

    // ----------- step 1 --------------------
    nceNo = 0;
    vol = sl.vertex[0]->hed->face->vol->getId();
    he = 0;
    for (size_t i=0;i<sl.edge.size();i++)
    {
        if ( (sl.edge)[i]->face->vol->getId()==vol && (sl.edge)[i] == nce[nceNo]) {he = (sl.edge)[i];break;}
    }
    assert(he!=0); // Halfedge not found
    sl.handleImpedingEdges(he);
    gsVolumeSegment<> vst(&sl);
    vst.volSegOp().costAuxEdge = 2;
    vst.volSegOp().PlanarityCostWeight = 100;
    cuttingLoop = vst.FindCuttingLoop(he);
    gsInfo<<"\n Segmenting loop found: "<< *cuttingLoop;
    gsVolumeSegment<>::splitFaces(sl, *cuttingLoop, 0.1, 0.1, 20, 0.001, false, 0.5);

    //
    graphPath = cuttingLoop->computePath();
    gsWriteParaview( sl , "tdpEGsl1", edgePointn, 0, linewidth, gsVector3d<>(0,0,0), 0, 10, 5, graphPath);
    //
    curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cuttingLoop);
    cface = vst.cuttingSurface(*cuttingLoop, *curveLoop, kv,kv,intopt);
    delete curveLoop;
    m = cface->toMesh(numSamples);
    gsWriteParaview(*m , "tdpCS1");

    graphPath = cuttingLoop->computePath();
    cuttingLoopVerts.clear();
    for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl.vertex[graphPath[i]]);
    sl.addFaceWithMate(cuttingLoopVerts, cface);

    gsWriteParaview(sl , "tdpEG11", edgePointn, vol, linewidth);
    gsWriteParaview(sl , "tdpEG12", edgePointn, sl.numVolumes-1, linewidth);
    // -- cut --
    delete cuttingLoop;
  }
  //------------------------------------------------------------------------------
  // Halfed Vase
  else if (basename=="ffoVaseHalfed.xml")
  {
    gsWriteParaview( sl , "vhEG1", edgePointn, 0, linewidth);
    ncEdgeV1.clear();
    ncEdgeV1.push_back(5);ncEdgeV2.push_back(11);
    ncEdgeV1.push_back(3);ncEdgeV2.push_back(5);
    ncEdgeV1.push_back(4);ncEdgeV2.push_back(5);
    nce = sl.detectNonConvexEdges(ncEdgeV1,ncEdgeV2);
    //gsVolumeSegment<> vs0 = gsVolumeSegment<>(&sl);
    //vs=vs0;

    // ----------- step 1 --------------------
    nceNo = 0;
    vol = sl.vertex[5]->hed->face->vol->getId();
    he = 0;
    for (size_t i=0;i<sl.edge.size();i++)
    {
        if ( (sl.edge)[i]->face->vol->getId()==vol && (sl.edge)[i] == nce[nceNo]) {he = (sl.edge)[i];break;}
    }
    assert(he!=0); // Halfedge not found
    //gsVolumeSegment<> vs5 = gsVolumeSegment<>(&sl); //TODO: only do this for now, should not completely re-compute the adjacent matrix but update it
    //vs=vs5;
    // -- cut --
    //sl.impedingEdges(he);
    //sl.insertNewVertex(he);
    sl.handleImpedingEdges(he);
    sl.insertNewVertex(sl.vertex[0]->getHalfEdge(sl.vertex[1]));
    gsVolumeSegment<> nvs1 = gsVolumeSegment<>(&sl);
    vs=nvs1;
    //vs.volSegOp().update(100,"PlanarityCostWeight");
    //gsInfo<<"\n"<< "\n ------------ here :"<< vs.volSegOp().PlanarityCostWeight<<"\n";
    vs.volSegOp().costAuxEdge = 0;
    vs.volSegOp().PlanarityCostWeight = 10000;
    cuttingLoop = vs.FindCuttingLoop(he);
    gsInfo<<"\n Segmenting loop found: "<< *cuttingLoop;
    gsVolumeSegment<>::splitFaces(sl, *cuttingLoop, 0.01, 0.1, 20, 0.001, false, 0.9);
    curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cuttingLoop);
    cface = vs.cuttingSurface(*cuttingLoop, *curveLoop, kv,kv,intopt);
    delete curveLoop;
    m = cface->toMesh(numSamples);
    gsWriteParaview(*m , "vhCS1");

    graphPath = cuttingLoop->computePath();
    gsWriteParaview( sl , "vhEGsl1", edgePointn, 0, linewidth, gsVector3d<>(0,0,0), 0, 10, 5, graphPath);
    cuttingLoopVerts.clear();
    for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl.vertex[graphPath[i]]);
    sl.addFaceWithMate(cuttingLoopVerts, cface);

    gsWriteParaview(sl , "vhEG11", edgePointn, vol, linewidth);
    gsWriteParaview(sl , "vhEG12", edgePointn, sl.numVolumes-1, linewidth);
    // -- cut --
    delete cuttingLoop;
  }
  //------------------------------------------------------------------------------
  // Science Park 2
  else if (basename=="ffoSciencePark2.xml" || basename =="ffoSciencePark2_lowTail.xml")
  {
    ncEdgeV1.clear();
    ncEdgeV1.push_back(13);ncEdgeV2.push_back(22);
    ncEdgeV1.push_back(12);ncEdgeV2.push_back(29);
    ncEdgeV1.push_back(13);ncEdgeV2.push_back(17);
    ncEdgeV1.push_back(17);ncEdgeV2.push_back(19);
    ncEdgeV1.push_back(10);ncEdgeV2.push_back(21);
    ncEdgeV1.push_back(5);ncEdgeV2.push_back(14);
    ncEdgeV1.push_back(16);ncEdgeV2.push_back(17);
    //ncEdgeV1.push_back(17);ncEdgeV2.push_back(19);
    ncEdgeV1.push_back(18);ncEdgeV2.push_back(19); // edge 7

    nce = sl.detectNonConvexEdges(ncEdgeV1,ncEdgeV2);
    gsWriteParaview( sl , "spEG0", edgePointn, 0, (real_t)(0.1));
    gsVolumeSegment<> vs0 = gsVolumeSegment<>(&sl);
    vs=vs0;

    // ----------- step 0 --------------------
    v1 = 13; v2 = 22;
    cuttingLoop = vs.FindCuttingLoop(v1,v2);
    gsVolumeSegment<>::splitFaces(sl, *cuttingLoop, 0.1, 0.1, 20, 0.001, false, 0.5);
    curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cuttingLoop);

    intopt.update("wReg",.01);
    intopt.update("wdPoint",1);
    intopt.update("wdNormal",.01);
    cface = vs.cuttingSurface(*cuttingLoop, *curveLoop, kv,kv,intopt,true);
    delete curveLoop;
    m = cface->toMesh(numSamples);
    gsWriteParaview(*m , "spSS0");

    graphPath = cuttingLoop->computePath();
    for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl.vertex[graphPath[i]]);
    sl.addFaceWithMate(cuttingLoopVerts, cface);
    gsWriteParaview( sl , "spEGloop0", 3, 0, linewidth, gsVector3d<>(0,0,0), 0, 10, 7, graphPath);

    gsWriteParaview(sl , "spEG00", edgePointn, 0, (real_t)(0.1));
    //gsWriteParaview(sl , "spEG01", edgePointn, 1, 0.1, gsVector3d<>(-1.5,0,1.5));
    gsWriteParaview(sl , "spEG01", edgePointn, 1, (real_t)(0.1));
    delete cuttingLoop;


    // ----------- step 01 --------------------
    cuttingLoop = vs.FindCuttingLoop(nce[1]);
    gsVolumeSegment<>::splitFaces(sl, *cuttingLoop, 0.1, 0.1, 20, 0.001, false, 0.5);
    curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cuttingLoop);

//  	gsKnotVector<> kv(0,1,0,deg+1);
//  	gsInterpOption<> intopt(deg);
    cface = vs.cuttingSurface(*cuttingLoop, *curveLoop, kv,kv,intopt);
    delete curveLoop;
    m = cface->toMesh(numSamples);
    gsWriteParaview(*m , "spCS01");

    graphPath = cuttingLoop->computePath();
    cuttingLoopVerts.clear();
    for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl.vertex[graphPath[i]]);
    sl.addFaceWithMate(cuttingLoopVerts, cface);

    gsWriteParaview(sl , "spEG010", edgePointn, 1, (real_t)(0.1));
    gsWriteParaview(sl , "spEG011", edgePointn, 2, (real_t)(0.1));
    delete cuttingLoop;


    // ----------- step 010: volume 1 --------------------
    he = nce[0];
    vol = sl.vertex[24]->hed->face->vol->getId();
//  	for (size_t i=0;i<sl.edge.size();i++)
//  	{
//      	if ( (sl.edge)[i]->face->vol->getId()==vol && (sl.edge)[i]->source == he->source && (sl.edge)[i]->mate->source == he->mate->source)
//      	{
//          	he = (sl.edge)[i];
//          	break;
//      	}
//  	} // this only works when removing nonconvex properties of the new edges
    for (size_t i=0;i<sl.edge.size();i++)
    {
        if ( (sl.edge)[i]->face->vol->getId()==vol && (sl.edge)[i]->source->getId() == 24
            && (sl.edge)[i]->mate->source->getId() == 27)
        {
            he = (sl.edge)[i];
            break;
        }
    }

    gsVolumeSegment<> vs1 = gsVolumeSegment<>(&sl); //TODO: only do this for now, should not completely re-compute the adjacent matrix but update it
    vs=vs1;
    //return 0;
    cuttingLoop = vs.FindCuttingLoop(he);
    gsVolumeSegment<>::splitFaces(sl, *cuttingLoop, 0.1, 0.1, 20, 0.001, false, 0.5);
    curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cuttingLoop);

//  	gsKnotVector<> kv(0,1,0,deg+1);
//  	gsInterpOption<> intopt(deg);
    cface = vs.cuttingSurface(*cuttingLoop, *curveLoop, kv,kv,intopt);
    delete curveLoop;
    m = cface->toMesh(numSamples);
    gsWriteParaview(*m , "spCS010");

    graphPath = cuttingLoop->computePath();
    cuttingLoopVerts.clear();
    for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl.vertex[graphPath[i]]);
    sl.addFaceWithMate(cuttingLoopVerts, cface);

    gsWriteParaview(sl , "spEG0100", edgePointn, vol, (real_t)(0.1));
    gsWriteParaview(sl , "spEG0101", edgePointn, sl.numVolumes-1, (real_t)(0.1));
    delete cuttingLoop;


    // ----------- step 00: volume 0 --------------------

    //he = nce[6]; // v 16 and 17
    he = nce[5]; // v 5 and 14
    vol = sl.vertex[16]->hed->face->vol->getId();
//  	for (size_t i=0;i<sl.edge.size();i++)
//  	{
//      	if ( (sl.edge)[i]->face->vol->getId()==vol && (sl.edge)[i]->source == he->source && (sl.edge)[i]->mate->source == he->mate->source)
//      	{
//          	he = (sl.edge)[i];
//          	break;
//      	}
//  	} // this only works when removing nonconvex properties of the new edges
    for (size_t i=0;i<sl.edge.size();i++)
    {
        if ( (sl.edge)[i]->face->vol->getId()==vol && (sl.edge)[i]->source->getId() == 5
            && (sl.edge)[i]->mate->source->getId() == 14)
        {
            he = (sl.edge)[i];
            break;
        }
    }
    gsVolumeSegment<> vs2 = gsVolumeSegment<>(&sl); //TODO: only do this for now, should not completely re-compute the adjacent matrix but update it
    vs=vs2;
    //return 0;
    cuttingLoop = vs.FindCuttingLoop(he);
    gsInfo<<"\n Segmenting loop found: "<< *cuttingLoop;
    gsVolumeSegment<>::splitFaces(sl, *cuttingLoop, 0.1, 0.1, 20, 0.001, false, 0.5);
    curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cuttingLoop);
    //curveLoop = sl.face[0]->surf->boundaryLoop();
    cface = vs.cuttingSurface(*cuttingLoop, *curveLoop, kv,kv,intopt);
    delete curveLoop;
    m = cface->toMesh(numSamples);
    gsWriteParaview(*m , "spCS00");

    graphPath = cuttingLoop->computePath();
    cuttingLoopVerts.clear();
    for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl.vertex[graphPath[i]]);
    sl.addFaceWithMate(cuttingLoopVerts, cface);

    gsWriteParaview(sl , "spEG000", edgePointn, vol, (real_t)(0.1));
    gsWriteParaview(sl , "spEG001", edgePointn, sl.numVolumes-1, (real_t)(0.1));
    delete cuttingLoop;

    // ----------- step 000 --------------------

    he = nce[4]; // v 10 and 21
    vol = sl.vertex[20]->hed->face->vol->getId();
    for (size_t i=0;i<sl.edge.size();i++)
    {
        if ( (sl.edge)[i]->face->vol->getId()==vol && (sl.edge)[i]->source->getId() == 10
            && (sl.edge)[i]->mate->source->getId() == 21)
        {
            he = (sl.edge)[i];
            break;
        }
    }
    gsVolumeSegment<> vs3 = gsVolumeSegment<>(&sl); //TODO: only do this for now, should not completely re-compute the adjacent matrix but update it
    vs=vs3;
    //return 0;
    cuttingLoop = vs.FindCuttingLoop(he);
    gsInfo<<"\n Segmenting loop found: "<< *cuttingLoop;
    gsVolumeSegment<>::splitFaces(sl, *cuttingLoop, 0.1, 0.1, 20, 0.001, false, 0.5);
    curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cuttingLoop);
    //curveLoop = sl.face[0]->surf->boundaryLoop();
    cface = vs.cuttingSurface(*cuttingLoop, *curveLoop, kv,kv,intopt);
    delete curveLoop;
    m = cface->toMesh(numSamples);
    gsWriteParaview(*m , "spCS000");

    graphPath = cuttingLoop->computePath();
    cuttingLoopVerts.clear();
    for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl.vertex[graphPath[i]]);
    sl.addFaceWithMate(cuttingLoopVerts, cface);

    gsWriteParaview(sl , "spEG0000", edgePointn, vol, (real_t)(0.1));
    gsWriteParaview(sl , "spEG0001", edgePointn, sl.numVolumes-1, (real_t)(0.1));
    delete cuttingLoop;

    // ----------- step 0001 --------------------
    nceNo = 6;
    vol = sl.numVolumes-1;
    he = 0;
    for (size_t i=0;i<sl.nHalfEdges();i++)
    {
        if ( (sl.edge)[i]->face->vol->getId()==vol && (sl.edge)[i]->isEquiv(nce[nceNo])==true )
        {he = (sl.edge)[i];break;}
    }
    assert(he!=0); // Halfedge not found
    gsVolumeSegment<> vs4 = gsVolumeSegment<>(&sl); //TODO: only do this for now, should not completely re-compute the adjacent matrix but update it
    vs=vs4;
    cuttingLoop = vs.FindCuttingLoop(he);
    gsInfo<<"\n Segmenting loop found: "<< *cuttingLoop;
    gsVolumeSegment<>::splitFaces(sl, *cuttingLoop, 0.1, 0.1, 20, 0.001, false, 0.5);
    graphPath = cuttingLoop->computePath();
    gsWriteParaview( sl , "spSL0001", edgePointn, vol, linewidth, gsVector3d<>(0,0,0), 0, 10, 5, graphPath);

    curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cuttingLoop);
    cface = vs.cuttingSurface(*cuttingLoop, *curveLoop, kv,kv,intopt);
    delete curveLoop;
    m = cface->toMesh(numSamples);
    gsWriteParaview(*m , "spCS0001");

    sl.checkStructure();

    cuttingLoopVerts.clear();
    for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl.vertex[graphPath[i]]);
    sl.addFaceWithMate(cuttingLoopVerts, cface);
    sl.checkStructure();

    gsWriteParaview(sl , "spEG00010", edgePointn, vol, (real_t)(0.1));
    gsWriteParaview(sl , "spEG00011", edgePointn, sl.numVolumes-1, (real_t)(0.1));
    delete cuttingLoop;

    // ----------- step 00011 --------------------
    nceNo = 7;
    vol = sl.numVolumes-1;
    he = 0;
    for (size_t i=0;i<sl.nHalfEdges();i++)
    {
        if ( (sl.edge)[i]->face->vol->getId()==vol && (sl.edge)[i]->isEquiv(nce[nceNo])==true )
        {he = (sl.edge)[i];break;}
    }
    assert(he!=0); // Halfedge not found
    real_t tol=1e-7;
    gsSolid<>::gsSolidHeVertexHandle tenso,tenta,so,ta;
    so = (sl.vertex)[7];
    ta = (sl.vertex)[11];
    gsSolidHalfEdge<>* he1 = 0;
    for (size_t i=0;i<sl.nHalfEdges();i++)
    {
        tenso = (sl.edge)[i]->source;
        tenta = (sl.edge)[i]->target();
        if ( (sl.edge)[i]->face->vol->getId()==vol &&  so->isEquiv(tenso,tol) && ta->isEquiv(tenta,tol) )
        {he1 = (sl.edge)[i];break;}
    }
    assert(he1!=0); // Halfedge not found
    sl.checkStructure();
    sl.handleImpedingEdges(he);
    sl.insertNewVertex(he1);
    gsVolumeSegment<> vs5 = gsVolumeSegment<>(&sl); //TODO: only do this for now, should not completely re-compute the adjacent matrix but update it
    vs=vs5;
    cuttingLoop = vs.FindCuttingLoop(he);
    gsInfo<<"\n Segmenting loop found: "<< *cuttingLoop;
    gsVolumeSegment<>::splitFaces(sl, *cuttingLoop, 0.1, 0.1, 20, 0.001, false, 0.5);
    curveLoop = gsVolumeSegment<>::curveLoopFromCuttingLoop(sl, *cuttingLoop);
    cface = vs.cuttingSurface(*cuttingLoop, *curveLoop, kv,kv,intopt);
    delete curveLoop;
    m = cface->toMesh(numSamples);
    gsWriteParaview(*m , "spCS00011");

    graphPath = cuttingLoop->computePath();
    gsWriteParaview( sl , "spSL00011", edgePointn, vol, linewidth, gsVector3d<>(0,0,0), 0, 10, 5, graphPath);
    cuttingLoopVerts.clear();
    for(size_t i = 0; i < graphPath.size(); i++) cuttingLoopVerts.push_back(sl.vertex[graphPath[i]]);
    sl.addFaceWithMate(cuttingLoopVerts, cface);

    gsWriteParaview(sl , "spEG000110", edgePointn, vol, (real_t)(0.1));
    gsWriteParaview(sl , "spEG000111", edgePointn, sl.numVolumes-1, (real_t)(0.1));
    delete cuttingLoop;
  }
  else
  {
    gsInfo << "I am only working for the files "
        <<"ffoTriDoublePyramid.xml, "
        <<"ffoTriBiPyramid.xml, "
        <<"ffoVaseHalfed.xml, "
        <<"ffoSciencePark2.xml, "
        <<"ffoSciencePark2_lowTail.xml"
      <<"\n";
  }

  gsFileData<> newdata;
  newdata << sl;
  newdata.dump("dump_write");

  gsInfo<<"\n";return 0;
}
