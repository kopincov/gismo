
#include <gismo.h>

#include <iostream>



using namespace gismo;

int main(int argc, char *argv[])
{

    //measuring the computational time
    //int clo=clock();
    std::string filename = "basis_thbs.xml";
    gsHTensorBasis<2>::uPtr hbs;
    gsCmdLine cmd("Hi, give me a file (.xml) with some hierarchical splines.");
    cmd.addPlainString("filename", "File containing hierarchical splines (.xml).", filename);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // The file data
    gsFileData<>  data( filename );
    hbs = data.getAnyFirst< gsHTensorBasis<2> >();

    if ( hbs )
    {
       gsInfo<< "  Got "<< *hbs << "\n";
    }
    else
    {
        gsInfo<< " Nothing interesting found in the file, quitting..\n";
        return 0;
    }

    //hbs->printCharMatrix();

  gsInfo<<"\n"<< "Size of the basis "<<hbs->size()<<"\n";
  gsInfo<<"The tree has "<< hbs->treeSize() << " nodes.\n" << "\n";


  //gsMatrix<> grev  = HB.anchors();
  //gsInfo<<"\n"<< "The Greville points: "<< "\n";
  //gsInfo<<  grev << "\n";

  // Testing evaluation
  gsMatrix<> para  = hbs->support();
  gsInfo<<"\n"<< "The parameter range is: "<< "\n" << para <<"\n";

  gsInfo<<"Num. of leaves: "<< hbs->tree().leafSize() <<".\n";

  hbs->tree().printLeaves();

  gsVector<> c0 = para.col(0);
  gsVector<> c1 = para.col(1);
  gsMatrix<> pts = uniformPointGrid(c0,c1, 11) ;
  gsMatrix<>   ev  = hbs->eval( pts ) ;
  gsMatrix<index_t> act = hbs->active( pts ) ;
  gsInfo<<"points\n"<< pts   <<"\n";
  gsInfo<<"eval  \n"<<  ev    <<"\n";
  gsInfo<<"act   \n"<<  act   <<"\n";
  gsInfo<<"Sums   \n"<< ev.colwise().sum()   <<"\n"; //----checkMatrixPartitionOfUnitClose

  gsMatrix<>   evs  = hbs->evalSingle( 0, pts ) ;
  // for all i = 0..act.cols()
  // if act.col(i) contains 0 at position k
  //  check that evs(0,k) == ev(act(k,i),i)
  // gsInfo<<"eval  0: \n"<< evs    <<"\n"; // checkMatrixClose(ev,evs)

  return 0;
}

