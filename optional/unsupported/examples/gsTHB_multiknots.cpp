#include <gismo.h>
#include <gsIO/gsIOUtils.h>
#include <iostream>




using namespace gismo;

int main(int argc, char *argv[])
{

//    //measuring the computational time
//    //int clo=clock();
//    std::string filename;
//    gsHTensorBasis<2> * hbs = NULL;
//    try
//    {
//      gsCmdLine cmd("Hi, give me a file (.xml) with some hierarchical splines.");

//      gsArgValPlain<std::string> a1("filename","File containing hierarchical splines.",
//                    false,"", "string",cmd );
//      cmd.parse(argc,argv);
//      filename = a1.getValue();

//    } catch ( gsArgException& e )
//    { gsInfo << "Error: " << e.error() << " " << e.argId() << "\n"; }

//    if ( ! filename.empty() )
//    {
//      // The file data
//      gsFileData<>  data( filename );
//      hbs = data.getAnyFirst< gsHTensorBasis<2> >();
//    }
//    else
//    {
//      gsInfo<< "Running default example.\n";
//      // The file data

//      
//      filename= "basis_thbs.xml"; //default example
//      gsFileData<>  data( filename );
//      hbs = data.getFirst< gsTHBSplineBasis<2> >();
//      gsInfo<<"file: "<< filename  <<"\n";
//    }

//    if ( hbs != 0 )
//    {
//       gsInfo<< "  Got "<< *hbs << "\n";
//    }
//    else
//    {
//        gsInfo<< " Nothing interesting found in the file, quitting..\n";
//        return 0;
//    }


    int deg_x = 2;
    int deg_y = 2;
    int kn = 1;
    gsKnotVector<> T_KV (0, 1, kn , deg_x+1, 1 ) ;
    gsKnotVector<> T_KV1 (0, 1,kn,deg_y+1,1) ;
    gsInfo<<"Knot Vector"<<T_KV<<"\n";
    gsInfo<<"Knot Vector"<<T_KV1<<"\n";

//    for (int i = 0; i < rand_knot_num; i ++){
//        int kn = rand()%T_KV1.size();
//        if(T_KV1.multiplicity(T_KV1[kn]) < deg_y){
//            T_KV1.insert(T_KV1[kn]);
//        }
//    }

    gsTensorBSplineBasis<2, real_t> T_tbasis( T_KV, T_KV1 );
    std::vector<index_t> boxes;
    boxes.push_back(1);
    boxes.push_back(0);
    boxes.push_back(0);
    boxes.push_back(2);
    boxes.push_back(2);
    //helps to create the refinement

    gsHBSplineBasis<2, real_t>  THB(T_tbasis, boxes);
    THB.printCharMatrix();
    THB.increaseMultiplicity(1, 0, 0.25, 1);
    THB.printCharMatrix();
    gsWriteParaview(THB,"THB_multipl_test",500);


//    gsMesh<double> mesh;
//    makeMesh(THB, mesh);

//    gsWriteParaview(mesh, "aaaa");

//  delete hbs;

  return 0;
}

