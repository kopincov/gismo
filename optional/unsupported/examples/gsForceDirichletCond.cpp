

#include <gismo.h>

#include <gsAssembler/gsForceDirichletConditions.h>

#include <iostream>
#include <vector>



using namespace gismo;

using std::flush;
;

#define TEST(a , b)\
    tmp2=(a);\
    passed = (passed && tmp2);\
    gsInfo << (tmp2? "TEST OK\n":"TEST FAIL\n");\
    if ( !tmp2 ){b;};

bool tmp2;
bool passed=true;

bool operator== (const gsSparseMatrix<> &a, const gsSparseMatrix<> &b )
{
    bool eq = (a.cols()==b.cols() && a.rows()==b.rows() );
    if (!eq) return false;
    for (int i=0; i<a.rows(); ++i)
        for (int j=0; j<a.cols(); ++j)
        {
            eq = eq && a(i,j)==b(i,j);
        }
    return eq;
}

template<int i>
void test();
template<int i>
void test2();


int main ()
{
    gsInfo<<"ColMajor\n";
    test<Eigen::ColMajor>();
    test2<Eigen::ColMajor>();
    gsInfo<<"RowMajor\n";
    test<Eigen::RowMajor>();
    test2<Eigen::RowMajor>();
    return passed? 0:1;
}

template <int S>
void test()
{

    gsSparseMatrix<real_t,S>  matrix, cmatrix;
    gsVector<>        vector, cvector;

    gsMatrix<unsigned>   dirichlet_dofs;
    gsMatrix<>           dirichlet_val;

    matrix.resize(5,5);
    cmatrix.resize(5,5);
    vector.setZero(5);
    cvector.setZero(5);

    for (int i=0; i<5; i++)
    {
        for (int j=0; j<5; j++)
            matrix(i,j)=i+5*j;
        vector(i)=i;
    }
    cmatrix=matrix;
    cvector=vector;
    for (int i=2; i<4; i++)
        for (int j=0; j<5; j++)
        {
            cmatrix(i,j)=0;
            cmatrix(j,i)=0;
        }
    cmatrix(2,2)=1;
    cmatrix(3,3)=1;
    cvector(2)=14;
    cvector(3)=15;
    cvector(0)-=365;
    cvector(1)-=394;
    cvector(2)=14;
    cvector(3)=15;
    cvector(4)-=481;

    dirichlet_dofs.resize(2,1);
    dirichlet_dofs<<2,3;
    dirichlet_val.resize(2,1);
    dirichlet_val<<14,15;
    gsForceDirichletConditions(matrix,vector,dirichlet_dofs,dirichlet_val);

    gsInfo<<"Test vector: ";
    TEST(vector==cvector, gsInfo<<vector<<"\n"<<cvector<<"\n");
    gsInfo<<"Test matrix: ";
    TEST(matrix==cmatrix, gsInfo<<matrix<<"\n"<<cmatrix<<"\n");

}


template <int S>
void test2()
{

    gsSparseMatrix<real_t,S>  matrix, cmatrix;
    gsVector<>          vector, cvector;

    gsMatrix<unsigned>   dirichlet_dofs;
    gsMatrix<>           dirichlet_val;

    matrix.resize(5,5);
    cmatrix.resize(5,5);
    vector.setZero(5);
    cvector.setZero(5);

    for (int i=0; i<5; i++)
    {
        for (int j=0; j<5; j++)
        {
            if(i!=j)
            {
                matrix(i,j)=i+5*j;
            }
        }
        vector(i)=i;
    }
    cmatrix=matrix;
    cvector=vector;
    for (int i=2; i<4; i++)
        for (int j=0; j<5; j++)
        {
            cmatrix(i,j)=0;
            cmatrix(j,i)=0;
        }
    cmatrix(2,2)=1;
    cmatrix(3,3)=1;
    cvector(2)=14;
    cvector(3)=15;
    cvector(0)-=365;
    cvector(1)-=394;
    cvector(2)=14;
    cvector(3)=15;
    cvector(4)-=481;

    dirichlet_dofs.resize(2,1);
    dirichlet_dofs<<2,3;
    dirichlet_val.resize(2,1);
    dirichlet_val<<14,15;
    gsForceDirichletConditions(matrix,vector,dirichlet_dofs,dirichlet_val);

    gsInfo<<"Test vector: ";
    TEST(vector==cvector, gsInfo<<vector<<"\n"<<cvector<<"\n");
    gsInfo<<"Test matrix: ";
    TEST(matrix==cmatrix, gsInfo<<matrix<<"\n"<<cmatrix<<"\n");

}
