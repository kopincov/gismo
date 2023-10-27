
#include <iostream>

#include <gismo.h>

#include <gsUtils/gsHeMesh.h>
#include <gsUtils/gsHeMesh.hpp>
#include <gsIO/gsWriteParaviewDev.h>
#include <gsBoxSplines/gsBoxSplineBasis.h>

using namespace gismo;


typedef gsMesh<> Mesh;

int main(int argc, char *argv[])
{
    bool plot = false; // If user gives --plot, paraview file is generated and launched on exit
    gsCmdLine cmd("Testing gsMesh.");    
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //HeMash-Testing
    typedef gsHeMesh<> HeMesh;
    typedef gsMeshElement<>::gsHeVertexHandle gsVertexHandle;
    //typedef gsMeshElement<>::gsHalfFaceHandle gsHalfFaceHandle;
    HeMesh hm;


    std::vector<gsVertexHandle> vh;

    HeMesh mesh;

    mesh.addVertex(1,0); //1 0
    mesh.addVertex(1,1); //2 1

    mesh.addVertex(2,0); //4 2
    mesh.addVertex(2,1); //5 3
    mesh.addVertex(3,1); //6 4

    mesh.addVertex(1,2); //8 5
    mesh.addVertex(2,2); //9 6
    mesh.addVertex(3,2); //10 7
    mesh.addVertex(4,2); //11 8

    mesh.addVertex(2,3); //13 9
    mesh.addVertex(3,3); //14 10
    mesh.addVertex(4,3); //15 11


    vh.push_back(mesh.getVertexAtId(0));vh.push_back(mesh.getVertexAtId(2));vh.push_back(mesh.getVertexAtId(3));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(3));vh.push_back(mesh.getVertexAtId(1));vh.push_back(mesh.getVertexAtId(0));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(2));vh.push_back(mesh.getVertexAtId(4));vh.push_back(mesh.getVertexAtId(3));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(1));vh.push_back(mesh.getVertexAtId(3));vh.push_back(mesh.getVertexAtId(6));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(6));vh.push_back(mesh.getVertexAtId(5));vh.push_back(mesh.getVertexAtId(1));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(3));vh.push_back(mesh.getVertexAtId(4));vh.push_back(mesh.getVertexAtId(7));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(7));vh.push_back(mesh.getVertexAtId(6));vh.push_back(mesh.getVertexAtId(3));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(4));vh.push_back(mesh.getVertexAtId(8));vh.push_back(mesh.getVertexAtId(7));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(5));vh.push_back(mesh.getVertexAtId(6));vh.push_back(mesh.getVertexAtId(9));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(6));vh.push_back(mesh.getVertexAtId(7));vh.push_back(mesh.getVertexAtId(10));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(10));vh.push_back(mesh.getVertexAtId(9));vh.push_back(mesh.getVertexAtId(6));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(7));vh.push_back(mesh.getVertexAtId(8));vh.push_back(mesh.getVertexAtId(11));
    mesh.addHalfFace(vh);vh.clear();

    vh.push_back(mesh.getVertexAtId(11));vh.push_back(mesh.getVertexAtId(10));vh.push_back(mesh.getVertexAtId(7));
    mesh.addHalfFace(vh);vh.clear();

    mesh.initialize();




    HeMesh mesh2; // equal to reference configuration
    mesh2.buildReferenceConfiguration();
    std::vector<double> y;
/*    y.push_back(0);
    y.push_back(1);
    y.push_back(2);
    y.push_back(3);
    y.push_back(4);
    y.push_back(5);
    y.push_back(6);
    y.push_back(7);
    y.push_back(8);
    y.push_back(9);
    y.push_back(10); */

/*    y.push_back(0);
    y.push_back(0.01);
    y.push_back(0.04);
    y.push_back(0.09);
    y.push_back(0.16);
    y.push_back(0.25);
    y.push_back(0.36);
    y.push_back(0.49);
    y.push_back(0.64);
    y.push_back(0.81);
    y.push_back(1); */

/*    y.push_back(0);
    y.push_back(0.001);
    y.push_back(0.008);
    y.push_back(0.027);
    y.push_back(0.064);
    y.push_back(0.125);
    y.push_back(0.216);
    y.push_back(0.343);
    y.push_back(0.512);
    y.push_back(0.729);
    y.push_back(1); */

/*    y.push_back(0);
    y.push_back(2*0.001);
    y.push_back(2*0.008);
    y.push_back(2*0.027);
    y.push_back(2*0.064);
    y.push_back(2*0.125);
    y.push_back(2*0.216);
    y.push_back(2*0.343);
    y.push_back(2*0.512);
    y.push_back(2*0.729);
    y.push_back(2*1); */

/*    y.push_back(0);
    y.push_back(0.015625);
    y.push_back(0.125);
    y.push_back(0.421875);
    y.push_back(1); */
//    mesh.interpolate(y,3);

    mesh.solveEquation(5,0,1);

    //Initialize refmeshfaces for 1Ring-evaluation tests
/*    HeMesh refmesh;
    refmesh.buildReferenceConfiguration();
    std::vector<gsHalfFaceHandle> refmeshfaces = refmesh.getHalfFaces();

    double inc1 = 0.0;
    double inc2 = 0.0;
    double inc3 = 0.0;
*/
    //Test 1: Input mesh is the reference configuration without the points (0,0),(0,1),(0,2),(1,3),(2,4),(3,4),(4,4) - evaluate the triangle (2,2),(3,2),(2,1)


/*

    std::vector<gsVertexHandle> vert = mesh.getVertices();
    for (std::vector<gsVertexHandle>::iterator
              it = vert.begin(); it!= vert.end(); ++it)
    {
        gsInfo << "Index: " << (*it)->getId() << " X: " << (*it)->x() << " Y: " << (*it)->y() << "\n \n";
    }


    gsInfo << "Triangle: X: " << mesh.getVertexAtId(6)->x() << " Y: " << mesh.getVertexAtId(6)->y() << "\n";
    gsInfo << "Triangle: X: " << mesh.getVertexAtId(3)->x() << " Y: " << mesh.getVertexAtId(3)->y() << "\n";
    gsInfo << "Triangle: X: " << mesh.getVertexAtId(7)->x() << " Y: " << mesh.getVertexAtId(7)->y() << "\n \n";
    gsInfo << "General triangle: " << "\n";
    mesh.evaluate1Ring(1,0.0,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)),refmeshfaces);
    mesh.evaluate1Ring(0,1.0,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)),refmeshfaces);
    mesh.evaluate1Ring(0,0.0,1.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)),refmeshfaces);

    //Inner triangle:
    gsInfo << "Graphics[{Red, Polygon[{{2, 2}, {2, 1}, {3, 2}}],";
    gsInfo << "Blue, Polygon[{";
    mesh.evaluate1Ring(0.5,0.25,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.5,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.25,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 1:
    gsInfo << "Gray, Polygon[{";
    mesh.evaluate1Ring(0.5,0.25,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.5,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.5,0.5,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 2:
    gsInfo << "Green, Polygon[{";
    mesh.evaluate1Ring(0.25,0.25,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.0,0.75,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.0,0.25,0.75,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 3:
    gsInfo << "Yellow, Polygon[{";
    mesh.evaluate1Ring(0.0,0.0,1,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.0,0.75,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.0,0.25,0.75,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 4:
    gsInfo << "Black, Polygon[{";
    mesh.evaluate1Ring(0.25,0.25,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.0,0.75,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.5,0.0,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 5:
    gsInfo << "Orange, Polygon[{";
    mesh.evaluate1Ring(0.25,0.25,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.0,0.25,0.75,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.0,0.5,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 6:
    gsInfo << "White, Polygon[{";
    mesh.evaluate1Ring(0.5,0.25,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.25,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.5,0.0,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 7:
    gsInfo << "Cyan, Polygon[{";
    mesh.evaluate1Ring(0.5,0.25,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.75,0.0,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.5,0.0,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 8:
    gsInfo << "Magenta, Polygon[{";
    mesh.evaluate1Ring(0.5,0.25,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.75,0.0,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.75,0.25,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 9:
    gsInfo << "Brown, Polygon[{";
    mesh.evaluate1Ring(1,0.0,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.75,0.0,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.75,0.25,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 10:
    gsInfo << "Pink, Polygon[{";
    mesh.evaluate1Ring(0.5,0.25,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.5,0.5,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.75,0.25,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 11:
    gsInfo << "Purple, Polygon[{";
    mesh.evaluate1Ring(0.25,0.25,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.5,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.0,0.5,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 12:
    gsInfo << "LightRed, Polygon[{";
    mesh.evaluate1Ring(0.0,0.75,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.5,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.0,0.5,0.5,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 13:
    gsInfo << "LightGreen, Polygon[{";
    mesh.evaluate1Ring(0.0,0.75,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.0,1.0,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.75,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 14:
    gsInfo << "LightBlue, Polygon[{";
    mesh.evaluate1Ring(0.0,0.75,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.5,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.75,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 15:
    gsInfo << "LightPink, Polygon[{";
    mesh.evaluate1Ring(0.5,0.5,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.5,0.25,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << ",";
    mesh.evaluate1Ring(0.25,0.75,0.0,mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7), mesh.getHalfFaceByBoundary(mesh.getVertexAtId(6), mesh.getVertexAtId(3), mesh.getVertexAtId(7)), refmeshfaces);
    gsInfo << "}]";
    gsInfo << "}]";
    gsInfo << "\n";

*/
    //Test 2: Input mesh is the reference configuration - evaluate the triangle (2,2),(3,2),(2,1)


/*
    std::vector<gsVertexHandle> vert = mesh2.getVertices();
    for (std::vector<gsVertexHandle>::iterator
              it = vert.begin(); it!= vert.end(); ++it)
    {
        gsInfo << "Index: " << (*it)->getId() << " X: " << (*it)->x() << " Y: " << (*it)->y() << "\n \n";
    }

    gsInfo << "Triangle: X: " << mesh2.getVertexAtId(1)->x() << " Y: " << mesh2.getVertexAtId(1)->y() << "\n";
    gsInfo << "Triangle: X: " << mesh2.getVertexAtId(2)->x() << " Y: " << mesh2.getVertexAtId(2)->y() << "\n";
    gsInfo << "Triangle: X: " << mesh2.getVertexAtId(3)->x() << " Y: " << mesh2.getVertexAtId(3)->y() << "\n \n";
    gsInfo << "General triangle: " << "\n";
    mesh2.evaluate1Ring(1,0.0,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    mesh2.evaluate1Ring(0,1.0,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    mesh2.evaluate1Ring(0,0.0,1.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);



    mesh2.evaluate1Ring(1,0.0,0.0,mesh2.getVertexAtId(8), mesh2.getVertexAtId(2), mesh2.getVertexAtId(1), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(8), mesh2.getVertexAtId(2), mesh2.getVertexAtId(1)),refmeshfaces);
    mesh2.evaluate1Ring(0,1.0,0.0,mesh2.getVertexAtId(8), mesh2.getVertexAtId(2), mesh2.getVertexAtId(1), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(8), mesh2.getVertexAtId(2), mesh2.getVertexAtId(1)),refmeshfaces);
    mesh2.evaluate1Ring(0,0.0,1.0,mesh2.getVertexAtId(8), mesh2.getVertexAtId(2), mesh2.getVertexAtId(1), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(8), mesh2.getVertexAtId(2), mesh2.getVertexAtId(1)),refmeshfaces);

    mesh2.evaluate1Ring(1,0.0,0.0,mesh2.getVertexAtId(8), mesh2.getVertexAtId(1), mesh2.getVertexAtId(7), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(8), mesh2.getVertexAtId(1), mesh2.getVertexAtId(7)),refmeshfaces);
    mesh2.evaluate1Ring(0,1.0,0.0,mesh2.getVertexAtId(8), mesh2.getVertexAtId(1), mesh2.getVertexAtId(7), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(8), mesh2.getVertexAtId(1), mesh2.getVertexAtId(7)),refmeshfaces);
    mesh2.evaluate1Ring(0,0.0,1.0,mesh2.getVertexAtId(8), mesh2.getVertexAtId(1), mesh2.getVertexAtId(7), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(8), mesh2.getVertexAtId(1), mesh2.getVertexAtId(7)),refmeshfaces);

    mesh2.evaluate1Ring(1,0.0,0.0,mesh2.getVertexAtId(7), mesh2.getVertexAtId(1), mesh2.getVertexAtId(6), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(7), mesh2.getVertexAtId(1), mesh2.getVertexAtId(6)),refmeshfaces);
    mesh2.evaluate1Ring(0,1.0,0.0,mesh2.getVertexAtId(7), mesh2.getVertexAtId(1), mesh2.getVertexAtId(6), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(7), mesh2.getVertexAtId(1), mesh2.getVertexAtId(6)),refmeshfaces);
    mesh2.evaluate1Ring(0,0.0,1.0,mesh2.getVertexAtId(7), mesh2.getVertexAtId(1), mesh2.getVertexAtId(6), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(7), mesh2.getVertexAtId(1), mesh2.getVertexAtId(6)),refmeshfaces);

    mesh2.evaluate1Ring(1,0.0,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(5), mesh2.getVertexAtId(6), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(5), mesh2.getVertexAtId(6)),refmeshfaces);
    mesh2.evaluate1Ring(0,1.0,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(5), mesh2.getVertexAtId(6), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(5), mesh2.getVertexAtId(6)),refmeshfaces);
    mesh2.evaluate1Ring(0,0.0,1.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(5), mesh2.getVertexAtId(6), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(5), mesh2.getVertexAtId(6)),refmeshfaces);

    mesh2.evaluate1Ring(1,0.0,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(3), mesh2.getVertexAtId(5), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(3), mesh2.getVertexAtId(5)),refmeshfaces);
    mesh2.evaluate1Ring(0,1.0,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(3), mesh2.getVertexAtId(5), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(3), mesh2.getVertexAtId(5)),refmeshfaces);
    mesh2.evaluate1Ring(0,0.0,1.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(3), mesh2.getVertexAtId(5), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(3), mesh2.getVertexAtId(5)),refmeshfaces);







    //Inner triangle:
    gsInfo << "Graphics[{Red, Polygon[{{2, 2}, {2, 1}, {3, 2}}],";
    gsInfo << "Blue, Polygon[{";
    mesh2.evaluate1Ring(0.5,0.25,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.5,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.25,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 1:
    gsInfo << "Gray, Polygon[{";
    mesh2.evaluate1Ring(0.5,0.25,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.5,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.5,0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 2:
    gsInfo << "Green, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.25,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.0,0.75,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.25,0.75,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 3:
    gsInfo << "Yellow, Polygon[{";
    mesh2.evaluate1Ring(0.0,0.0,1,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.0,0.75,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.25,0.75,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 4:
    gsInfo << "Black, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.25,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.0,0.75,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.0,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 5:
    gsInfo << "Orange, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.25,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.25,0.75,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.5,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 6:
    gsInfo << "White, Polygon[{";
    mesh2.evaluate1Ring(0.5,0.25,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.25,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.0,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 7:
    gsInfo << "Cyan, Polygon[{";
    mesh2.evaluate1Ring(0.5,0.25,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75,0.0,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.0,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 8:
    gsInfo << "Magenta, Polygon[{";
    mesh2.evaluate1Ring(0.5,0.25,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75,0.0,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75,0.25,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 9:
    gsInfo << "Brown, Polygon[{";
    mesh2.evaluate1Ring(1,0.0,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75,0.0,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75,0.25,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 10:
    gsInfo << "Pink, Polygon[{";
    mesh2.evaluate1Ring(0.5,0.25,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.5,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75,0.25,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 11:
    gsInfo << "Purple, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.25,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.5,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.5,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 12:
    gsInfo << "LightRed, Polygon[{";
    mesh2.evaluate1Ring(0.0,0.75,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.5,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.5,0.5,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 13:
    gsInfo << "LightGreen, Polygon[{";
    mesh2.evaluate1Ring(0.0,0.75,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,1.0,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.75,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 14:
    gsInfo << "LightBlue, Polygon[{";
    mesh2.evaluate1Ring(0.0,0.75,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.5,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.75,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 15:
    gsInfo << "LightPink, Polygon[{";
    mesh2.evaluate1Ring(0.5,0.5,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.5,0.25,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.75,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}]";
    gsInfo << "}]";
    gsInfo << "\n";

*/

    //Test 3: Input mesh is the reference configuration - evaluate the triangle (2,2),(3,2),(2,1); goal is a lightly curved triangle
    //(-0.25,+0.125,+0.125)

/*
    std::vector<gsVertexHandle> vert = mesh2.getVertices();
    for (std::vector<gsVertexHandle>::iterator
              it = vert.begin(); it!= vert.end(); ++it)
    {
        gsInfo << "Index: " << (*it)->getId() << " X: " << (*it)->x() << " Y: " << (*it)->y() << "\n \n";
    }

    gsInfo << "Triangle: X: " << mesh2.getVertexAtId(1)->x() << " Y: " << mesh2.getVertexAtId(1)->y() << "\n";
    gsInfo << "Triangle: X: " << mesh2.getVertexAtId(2)->x() << " Y: " << mesh2.getVertexAtId(2)->y() << "\n";
    gsInfo << "Triangle: X: " << mesh2.getVertexAtId(3)->x() << " Y: " << mesh2.getVertexAtId(3)->y() << "\n \n";
    gsInfo << "General triangle: " << "\n";
    mesh2.evaluate1Ring(1,0.0,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    mesh2.evaluate1Ring(0.0,1.0,0.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    mesh2.evaluate1Ring(0.0,0.0,1.0,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);


    //Inner triangle:
    gsInfo << "Graphics[{Red, Polygon[{{2, 2}, {2, 1}, {3, 2}}],";
    gsInfo << "Blue, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.375,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.625,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.375,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 1:
    gsInfo << "Gray, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.375,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.625,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.625,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 2:
    gsInfo << "Green, Polygon[{";
    mesh2.evaluate1Ring(0.0,0.375,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.125,0.875,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(-0.25,0.375,0.875,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 3:
    gsInfo << "Yellow, Polygon[{";
    mesh2.evaluate1Ring(-0.25,0.125,1.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.125,0.875,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(-0.25,0.375,0.875,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 4:
    gsInfo << "Black, Polygon[{";
    mesh2.evaluate1Ring(0.0,0.375,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.125,0.875,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.125,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 5:
    gsInfo << "Orange, Polygon[{";
    mesh2.evaluate1Ring(0.0,0.375,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(-0.25,0.375,0.875,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(-0.25,0.625,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 6:
    gsInfo << "White, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.375,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.375,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.125,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 7:
    gsInfo << "Cyan, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.375,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.125,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.125,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 8:
    gsInfo << "Magenta, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.375,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.125,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.375,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 9:
    gsInfo << "Brown, Polygon[{";
    mesh2.evaluate1Ring(0.75,0.125,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.125,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.375,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 10:
    gsInfo << "Pink, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.375,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25,0.625,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5,0.375,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 11:
    gsInfo << "Purple, Polygon[{";
    mesh2.evaluate1Ring(0.0,0.375,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.625,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(-0.25,0.625,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 12:
    gsInfo << "LightRed, Polygon[{";
    mesh2.evaluate1Ring(-0.25,0.875,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.625,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(-0.25,0.625,0.625,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 13:
    gsInfo << "LightGreen, Polygon[{";
    mesh2.evaluate1Ring(-0.25,0.875,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(-0.25,1.125,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.875,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 14:
    gsInfo << "LightBlue, Polygon[{";
    mesh2.evaluate1Ring(-0.25,0.875,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.625,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.875,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 15:
    gsInfo << "LightPink, Polygon[{";
    mesh2.evaluate1Ring(0.25,0.625,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.625,0.375,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0,0.875,0.125,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}]";
    gsInfo << "}]";
    gsInfo << "\n";
*/

    //Test 4: Input mesh is the reference configuration - evaluate the triangle (2,2),(3,2),(2,1); goal are several triangle shifts
    //O1 (-0.25,0.125,0.125)
    //O2 (-0.5,0.25,0.25)
    //O3 (-0.25,0.15,0.1)
    //O4 (0.125,-0.25,0.125)
    //O5 (0.5,-1.0,0.5)
    //O6 (0.375,-0.75,0.375)
    //O7 (0.375,0.375,-0.75)
    //O8 (0.125,0.125,-0.25)

/*
    inc1 = 0.125;
    inc2 = 0.125;
    inc3 = -0.25;

    gsInfo << "Triangle: X: " << mesh2.getVertexAtId(1)->x() << " Y: " << mesh2.getVertexAtId(1)->y() << "\n";
    gsInfo << "Triangle: X: " << mesh2.getVertexAtId(2)->x() << " Y: " << mesh2.getVertexAtId(2)->y() << "\n";
    gsInfo << "Triangle: X: " << mesh2.getVertexAtId(3)->x() << " Y: " << mesh2.getVertexAtId(3)->y() << "\n \n";
    gsInfo << "General triangle: " << "\n";
    mesh2.evaluate1Ring(1+inc1,0.0+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    mesh2.evaluate1Ring(0+inc1,1.0+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    mesh2.evaluate1Ring(0+inc1,0.0+inc2,1.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);



    //Inner triangle:
    gsInfo << "Graphics[{Red, Polygon[{{2, 2}, {2, 1}, {3, 2}}],";
    gsInfo << "Blue, Polygon[{";
    mesh2.evaluate1Ring(0.5+inc1,0.25+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.5+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.25+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 1:
    gsInfo << "Gray, Polygon[{";
    mesh2.evaluate1Ring(0.5+inc1,0.25+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.5+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5+inc1,0.5+inc2,0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 2:
    gsInfo << "Green, Polygon[{";
    mesh2.evaluate1Ring(0.25+inc1,0.25+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.0+inc2,0.75+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0+inc1,0.25+inc2,0.75+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 3:
    gsInfo << "Yellow, Polygon[{";
    mesh2.evaluate1Ring(0.0+inc1,0.0+inc2,1+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.0+inc2,0.75+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0+inc1,0.25+inc2,0.75+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 4:
    gsInfo << "Black, Polygon[{";
    mesh2.evaluate1Ring(0.25+inc1,0.25+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.0+inc2,0.75+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5+inc1,0.0+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 5:
    gsInfo << "Orange, Polygon[{";
    mesh2.evaluate1Ring(0.25+inc1,0.25+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0+inc1,0.25+inc2,0.75+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0+inc1,0.5+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 6:
    gsInfo << "White, Polygon[{";
    mesh2.evaluate1Ring(0.5+inc1,0.25+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.25+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5+inc1,0.0+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 7:
    gsInfo << "Cyan, Polygon[{";
    mesh2.evaluate1Ring(0.5+inc1,0.25+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75+inc1,0.0+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5+inc1,0.0+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 8:
    gsInfo << "Magenta, Polygon[{";
    mesh2.evaluate1Ring(0.5+inc1,0.25+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75+inc1,0.0+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75+inc1,0.25+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 9:
    gsInfo << "Brown, Polygon[{";
    mesh2.evaluate1Ring(1+inc1,0.0+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75+inc1,0.0+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75+inc1,0.25+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 10:
    gsInfo << "Pink, Polygon[{";
    mesh2.evaluate1Ring(0.5+inc1,0.25+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.5+inc1,0.5+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.75+inc1,0.25+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 11:
    gsInfo << "Purple, Polygon[{";
    mesh2.evaluate1Ring(0.25+inc1,0.25+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.5+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0+inc1,0.5+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 12:
    gsInfo << "LightRed, Polygon[{";
    mesh2.evaluate1Ring(0.0+inc1,0.75+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.5+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0+inc1,0.5+inc2,0.5+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 13:
    gsInfo << "LightGreen, Polygon[{";
    mesh2.evaluate1Ring(0.0+inc1,0.75+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.0+inc1,1.0+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.75+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 14:
    gsInfo << "LightBlue, Polygon[{";
    mesh2.evaluate1Ring(0.0+inc1,0.75+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.5+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.75+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}],";
    gsInfo << "\n";

    //Triangle 15:
    gsInfo << "LightPink, Polygon[{";
    mesh2.evaluate1Ring(0.5+inc1,0.5+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.5+inc2,0.25+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << ",";
    mesh2.evaluate1Ring(0.25+inc1,0.75+inc2,0.0+inc3,mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3), mesh2.getHalfFaceByBoundary(mesh2.getVertexAtId(1), mesh2.getVertexAtId(2), mesh2.getVertexAtId(3)),refmeshfaces);
    gsInfo << "}]";
    gsInfo << "}]";
    gsInfo << "\n";

*/

    if ( plot )
    {
        // Output a paraview file and open it
        gsWriteParaview<real_t>(hm , "HeEdgeGraph");

      // Call paraview on exit
      char cmdParaview[100];
      strcpy(cmdParaview,"paraview HeEdgeGraph.vtp\0");
      strcat(cmdParaview," &");  
      return system(cmdParaview);
    }
    else
    {
      return 0;
    }

}
