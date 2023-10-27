
#include <iostream>
#include <set>
#include <map>
#include <math.h>
#include <ctime>

#include <gismo.h>





using namespace gismo;
//lvl- max level of domain refinement, size in x and y direction of the inserted boxes, r-radius of a circle
std::vector<index_t> curvilinear_refinement(int lvl, int size_x, int size_y,int r, gsTHBSplineBasis<2> *bas){
    gsVector<index_t> i1;
    gsVector<index_t> i2;
    i1.resize(2);
    i2.resize(2);
    gsVector<> v1;
    gsVector<> v2;
    v1.resize(2);
    v2.resize(2);
    std::vector<index_t> q;
    for (int i = 1; i <= lvl; i++) {
        if(unsigned (lvl)>bas->tree().getIndexLevel()){
            gsInfo<<"The insertion ended at level "<<lvl<<"."<<"\n";
            gsInfo<<"For insertion of higher levels change the limit in the THB-basis"<<"\n";
            break;
        }
        //for (unsigned int j = 0; j < bas->m_cvs[0][i].size()-size_x; j++) {
        for (size_t j = 0; j < bas->getBases()[i]->component(0).knots().size()-size_x; j++) {
            for (size_t k = 0; k < bas->getBases()[i]->component(1).knots().size()-size_y; k++) {
                v1[0] = bas->getBases()[i]->component(0).knots()[j]; //v1[0] = bas->m_cvs[0][i][j];//values of knots
                v1[1] = bas->getBases()[i]->component(1).knots()[k]; //v1[1] = bas->m_cvs[1][i][k];
                v2[0] = bas->getBases()[i]->component(0).knots()[j + size_x]; //v2[0] = bas->m_cvs[0][i][j+size_x];
                v2[1] = bas->getBases()[i]->component(1).knots()[k + size_y]; //v2[1] = bas->m_cvs[1][i][k+size_y];
                real_t a = (v1[0]*v1[0] + v1[1]*v1[1]) - r*r;
                real_t b = (v2[0]*v2[0] + v2[1]*v2[1]) - r*r;

                real_t c = (v1[0]*v1[0] + v2[1]*v2[1]) - r*r;
                real_t d = (v2[0]*v2[0] + v1[1]*v1[1]) - r*r;


                if((a*b<0)||(c*d<0)){
                    //insert box to tree->necesary to change to global indices
                    q.push_back(i);
                    q.push_back(j);
                    q.push_back(k);
                    q.push_back(j+size_x);
                    q.push_back(k+size_y);
                }
            }
        }
    }
    return q;
}


//lvl- max level to be inserted, sixe of the inserted boxes in x and y direction, definition of a line ax+by+c = 0->define a,b,c
std::vector<index_t> line_refinement(int lvl, int size_x, int size_y, std::vector<real_t> line,  gsTHBSplineBasis<2> *bas){
    gsVector<index_t> i1;
    gsVector<index_t> i2;
    i1.resize(2);
    i2.resize(2);
    gsVector<> v1;
    gsVector<> v2;
    v1.resize(2);
    v2.resize(2);
    std::vector<index_t> q;
    real_t a, b, c, d;
    // For over all potential functions; if it contains a diagonal, add it to the tree.
    for (int i = 1; i <= lvl; i++) {
        if(unsigned (lvl)>bas->tree().getIndexLevel()){
            gsInfo<<"The insertion ended at level "<<lvl<<"."<<"\n";
            gsInfo<<"For insertion of higher levels change the limit in the THB-basis"<<"\n";
            break;
        }
        for (size_t j = 0; j < bas->getBases()[i]->component(0).knots().size()-size_x; j++) {//for (unsigned int j = 0; j < bas->m_cvs[0][i].size()-size_x; j++) {
            for (size_t k = 0; k < bas->getBases()[i]->component(1).knots().size()-size_y; k++) {//for (unsigned int k = 0; k < bas->m_cvs[1][i].size()-size_y; k++) {
                v1[0] = bas->getBases()[i]->component(0).knots()[j]; //bas->m_cvs[0][i][j];//values of knots
                v1[1] = bas->getBases()[i]->component(1).knots()[k]; //bas->m_cvs[1][i][k];
                v2[0] = bas->getBases()[i]->component(0).knots()[j+size_x]; //bas->m_cvs[0][i][j+size_x];
                v2[1] = bas->getBases()[i]->component(1).knots()[k+size_y]; //bas->m_cvs[1][i][k+size_y];
                //Test whether the box is not degenerated.
                if( (math::abs(v1[0]-v2[0])!=0)&& (math::abs(v1[1]-v2[1])!=0) )
                {
                    a = line[0]*v1[0] + line[1]*v1[1] + line[2];
                    b = line[0]*v2[0] + line[1]*v2[1] + line[2];

                    c = line[0]*v1[0] + line[1]*v2[1] + line[2];
                    d = line[0]*v2[0] + line[1]*v1[1] + line[2];
                    //test if the box cros the line
                    if((a*b<0)||(c*d<0)){
                        //insert box to tree->necesary to change to global indices
                        q.push_back(i);
                        q.push_back(j);
                        q.push_back(k);
                        q.push_back(j+size_x);
                        q.push_back(k+size_y);
                    }
                }
            }
        }
    }
    return q;
}

int main()
{
  gsVector<index_t> i1;
  gsVector<index_t> i2;
  i1.resize(2);
  i2.resize(2);



  ////////THB spline tests///////////////////////


    gsInfo << " ----- Test THBSplineBasis ----------  \n";
    gsKnotVector<> T_KV (0, 1, 13,3, 1 ) ;
    gsInfo<<"Knot Vector"<<T_KV<<"\n";
    gsBSplineBasis<> T_bsp( T_KV );
    gsTensorBSplineBasis<2,real_t > T_tbasis( new gsBSplineBasis<>(T_bsp) , new gsBSplineBasis<>(T_bsp) ) ;

    //empty basis used for creating the refinement
    gsTHBSplineBasis<2>  TT( T_tbasis );
    std::vector<real_t> line;
    line.resize(3);
    line[0]=1;
    line[1]=-1;
    line[2] =0;
    gsTHBSplineBasis<2>  THB( T_tbasis , line_refinement(2, 4,4,line, &TT) ) ;

    gsInfo<<"The tree has "<< THB.tree().size() << " nodes.\n" << "\n";
    //set the matrices
    //THB.initialize();
    gsMatrix<>  T_para  = THB.support();//wrong if sKnotVector<> T_KV (0, 1, 3,1, 1 ) ;
    gsInfo<<"\n"<< "The parameter range is: "<< "\n" << T_para <<"\n";

    gsVector<> T_c0 = T_para.col(0);
    gsVector<> T_c1 = T_para.col(1);
    gsMatrix<> T_pts = uniformPointGrid(T_c0,T_c1, 8) ;
    gsInfo<<"points\n"<< T_pts   <<"\n";

    gsMatrix<>  T_ev  = THB.eval( T_pts) ;
    gsInfo<<"eval  \n"<< T_ev    <<"\n";
    gsInfo<<"Sums   \n"<< T_ev.colwise().sum()   <<"\n";
  return 0;

}

