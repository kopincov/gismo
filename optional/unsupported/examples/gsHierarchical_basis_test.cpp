
#include <iostream>
#include <set>
#include <map>

#include <gismo.h>
#include <gismo_dev.h>




using namespace gismo;

int main(int argc, char *argv[])
{

    //measuring the computational time
    //int clo=clock();
    std::string filename = "basis_thbs_01.xml";
    gsHTensorBasis<2>::uPtr hbs;
    
    gsCmdLine cmd("Hi, give me a file (.xml) with some hierarchical splines.");
    cmd.addPlainString("filename", "File containing hierarchical splines (.xml).", filename);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // The file data
    gsFileData<>  data( filename );
    if ( data.has< gsHBSplineBasis<2> >() )
    {
        hbs = data.getFirst< gsHBSplineBasis<2> >();
    }
    if ( data.has< gsTHBSplineBasis<2> >() )
    {
        hbs = data.getFirst< gsTHBSplineBasis<2> >();
    }
    gsInfo<< "  Got "<< *hbs << "\n";

    hbs->printCharMatrix();

    gsInfo<<"Size of the basis "<<hbs->size()<<"\n";
    gsInfo<<"The tree has "<< hbs->tree().size() << " nodes." << "\n";
    gsMatrix<> ins;
    ins.resize(2,4);
    ins(0,0) = 0.0;
    ins(1,0) = 0;
    ins(0,1) = 1;
    ins(1,1) = 1;

    ins(0,2) = 0.0;
    ins(1,2) = 0;
    ins(0,3) = 0.5;
    ins(1,3) = 0.5;
    hbs->refine(ins);

    hbs->printCharMatrix();

    gsInfo<<"Size of the basis "<<hbs->size()<<"\n";
    gsInfo<<"The tree has "<< hbs->tree().size() << " nodes." << "\n";



    gsMatrix<> grev  = hbs->anchors();
    gsInfo<<"\n"<< "The Greville points: "<< "\n";
    gsInfo<<  grev << "\n"<< "\n";

    // Testing evaluation
    gsMatrix<> para  = hbs->support();
    gsInfo<< "The parameter range is: "<< "\n" << para <<"\n";

    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    gsMatrix<> pts = uniformPointGrid(c0,c1, 11) ;
    gsMatrix<>   ev  = hbs->eval( pts ) ;
    gsMatrix<index_t>  act = hbs->active( pts ) ;
    gsInfo<<"Partition of unity test \n"<< ev.colwise().sum()   <<"\n";

    hbs->evalSingle_into( 0, pts, ev ) ;

    ////////THB spline tests///////////////////////

    return 0;
}
