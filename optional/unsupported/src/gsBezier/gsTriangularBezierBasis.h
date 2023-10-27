/** @file gsTriangularBezierBasis.h

    @brief Provides declaration of TriangularBezierBasis class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): G. Kiss
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBasisFun.h>

#include <gsCore/gsDebug.h>

#include <gsCore/gsBoundary.h>

#include <gsCore/gsBasis.h>

#include <gsUtils/gsCombinatorics.h>//
 
namespace gismo
{

template<short_t d,class T>
class gsTriangularBezierBasis : public gsBasis<T>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
    /// Shared pointer for gsTriangularBezierBasis
    typedef memory::shared_ptr< gsTriangularBezierBasis > Ptr;

    /// Unique pointer for gsTriangularBezierBasis
    typedef memory::unique_ptr< gsTriangularBezierBasis > uPtr;

    typedef T Scalar_t;

    static const bool IsRational = false;

    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;

    /// Dimension of the parameter domain
    static const int Dim = d;

private:

    //unsigned int d;//dimension
    short_t m_degree;
    gsMatrix<T,3,2> m_domain;
    std::vector<gsVector<unsigned int> > m_compositions;

public:
    gsTriangularBezierBasis() {  }

    gsTriangularBezierBasis(unsigned int p) {
        // d = 2;
        m_degree = p;
        //m_domain.resize(3,2);
        m_domain(0,0) = 0;
        m_domain(0,1) = 0;
        m_domain(1,0) = 1;
        m_domain(1,1) = 0;
        m_domain(2,0) = 0;
        m_domain(2,1) = 1;
        getCompositions(m_compositions);
    }

    gsTriangularBezierBasis(unsigned int p, const gsMatrix<T,3,2> & domain) 
    {
        // d = 2;
        m_degree = p;
        //m_domain.resize(3,2);
        m_domain(0,0) = domain(0,0);
        m_domain(0,1) = domain(0,1);
        m_domain(1,0) = domain(1,0);
        m_domain(1,1) = domain(1,1);
        m_domain(2,0) = domain(2,0);
        m_domain(2,1) = domain(2,1);
        getCompositions(m_compositions);
    }

    /// Destructor
    ~gsTriangularBezierBasis() { };

public:

    // Look at gsBasis class for a description
    short_t domainDim() const;

    // Look at gsBasis class for a description
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
    { 
        const int sz = size();
        result.resize( sz, u.cols() );
        for ( int i = 0; i< sz; ++i )
            result.row(i).setConstant(i);
    }

    // Look at gsBasis class for a description
    gsMatrix<T> support() const;

    // Look at gsBasis class for a description
    gsMatrix<T> support(const index_t & i) const;

    // Look at gsBasis class for a description
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    // Look at gsBasis class for a description
    void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;


    // Look at gsBasis class for a description
    void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    GISMO_CLONE_FUNCTION(gsTriangularBezierBasis)

    // Look at gsBasis class for a description
    memory::unique_ptr<gsGeometry<T> > makeGeometry(gsMatrix<T> coefs ) const
    { GISMO_NO_IMPLEMENTATION }

    // Look at gsBasis class for a description
    domainIter makeDomainIterator() const
    { GISMO_NO_IMPLEMENTATION }


    // Look at gsBasis class for a description
    std::ostream &print(std::ostream &os) const
    {
        os <<"Triangular basis of degree "<< m_degree <<"\n";
        return os;
    }

    /// Prints the object as a string with extended details.
    std::string detail() const
    {
        // By default just uses print(..)
        std::ostringstream os;
        print(os);
        return os.str();
    }

    /*
      Member functions that may be implemented or not in the derived class
    */

    /// The number of basis functions in this basis.
    index_t size() const;
    
    short_t totalDegree() const
    { return m_degree; }

    short_t degree(short_t i = 0) const
    { return m_degree; }

//private:
    gsVector<T> getBaricentricCoordinate(gsVector<T>param)const;

    //check if the baricentric corrdinates in coord are all positive
    bool isInDomain (gsVector<T> coord)const;

    //create compositions for basis functions for m_degree
    void getCompositions( std::vector<gsVector<unsigned int> > & compos)const;

    //create compositions for basis functions for arbitrary degree
    void getCompositions( std::vector<gsVector<unsigned int> > & compos, unsigned degree)const;

    //rt derivative
    void rThDerivSingle(unsigned int r,unsigned int i, const gsMatrix<T> & dir, const gsMatrix<T> & u, gsMatrix<T>& result)const;

    ///find p in \a compos and return the index
    int findIndex(gsVector<unsigned int> const p, std::vector<gsVector<unsigned int> >const compos )const;

    //evaluation of single function by giving the baricentric coordinates of a point
    //u contains the baricentric coordinates
    void evalSingle_into_baric(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    //retunrs degree
    int get_degree()const{
        return m_degree;
    }

    gsMatrix<T,3,2> get_domain()const{
        return m_domain;
    }
}; // class gsTriangularBezierBasis


template<short_t d,class T>
short_t gsTriangularBezierBasis<d,T>::domainDim() const{ return 2; }

template<short_t d,class T>
index_t gsTriangularBezierBasis<d,T>::size() const{
    return  m_compositions.size();
}

template<short_t d,class T>
gsMatrix<T> gsTriangularBezierBasis<d,T>::support(const index_t & i) const{
    gsMatrix<T> res(d,2);

    for (unsigned j = 0; j < d; ++j)
    {
        T min= 1000000;
        T max = -10000000;
        for(int k = 0; k < m_domain.rows();k++)
        {
            if(m_domain(k,j) > max){
                max = m_domain(k,j);
            }
            if(m_domain(k,j) < min){
                min = m_domain(k,j);
            }
        }
        res(j,0) = min;
        res(j,1) = max;
    }
// res.row(i) =  m_bases[i]->support();
    return res;
}

template<short_t d,class T>
gsMatrix<T> gsTriangularBezierBasis<d,T>:: support() const{
    return support(0);
}
/////////////evaluations/////////////////
template<short_t d,class T>
void gsTriangularBezierBasis<d,T>::evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const{
    gsVector<T>  baric;
    gsVector<T>  point;
    result.setZero(1,u.cols() );
    for(int j = 0; j < u.cols(); j++)
    {
        point =  u.col(j);
        baric = getBaricentricCoordinate(point);
        if(isInDomain(baric))
        {
            T denom = 1.0;
            T temp = 1.0;
            for(int k = 0; k < m_compositions[i].size(); k++)
            {
                denom = denom * factorial(m_compositions[i][k]);
                temp = temp * math::pow(baric[k], static_cast<int>(m_compositions[i][k]));
            }
            result(0,j) =  (factorial(m_degree)/denom) * temp;
        }
    }
}

template<short_t d,class T>
void gsTriangularBezierBasis<d,T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const{
    result.setZero(m_compositions.size(), u.cols());
    gsMatrix<T> single_result;
    for(unsigned int i = 0; i < m_compositions.size(); i++)
    {
        evalSingle_into(i, u, single_result);
        result.row(i) = single_result.row(0);
    }
}

// de Casteljau alg. for triangular patches
//template<short_t d,class T>
//void gsTriangularBezierBasis<d,T>::evalFunc_into(const gsMatrix<T> &u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const 
// {
//
//  }

///////////firts derivatives/////////////////////
template<short_t d,class T>
void gsTriangularBezierBasis<d,T>::derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const{
    gsMatrix<T> coordinates;
    gsMatrix<T> tempresult;
    gsMatrix<T> directions(d,2);
    int debugme = -1;
    coordinates.setZero(d+1,d);
    result.resize(d,u.cols());
    //create matrix where each column is a poin to compute the directional derivatives in the u and v direction

    for(int j = 1; j < coordinates.rows(); j++){
        coordinates(j,j-1) = 1;
    }
    if(debugme > 0){
        std::cout<<"coordinates for the directions are: "<<std::endl<<coordinates<<std::endl;
    }
    if(debugme > 5){
        std::cout<<"poits are "<<std::endl<<u<<std::endl;
    }
    directions.row(0) = coordinates.row(0);
    for(unsigned int j = 0; j < d; j++){
        directions.row(1) = coordinates.row(j+1);
        if(debugme > 0){
            std::cout<<"directions for the first derivative are: "<<std::endl<<directions<<std::endl;
        }
        rThDerivSingle(1,i ,directions,u,tempresult);
        if(debugme > 5){
            std::cout<<"first derivatives are: "<<std::endl<<tempresult<<std::endl;
        }
        result.row(j) = tempresult.row(0);
    }
    //rThDerivSingle(1,i,u,result);

}

template<short_t d,class T>
void gsTriangularBezierBasis<d,T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const{
    result.setZero(m_compositions.size()*this->dim(), u.cols());
    gsMatrix<T> single_result;
    for(unsigned int i = 0; i < m_compositions.size(); i++)
    {
        derivSingle_into(i, u, single_result);
        for(int j = 0; j < single_result.rows();j++){
            result.row(2*i+j) = single_result.row(j);
        }
    }
}


//////////////second derivatives//////////////
template<short_t d,class T>
void gsTriangularBezierBasis<d,T>::deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{

    gsMatrix<T> coordinates;
    gsMatrix<T> tempresult;
    gsMatrix<T> directions(d,2);
    int debugme = -1;
    coordinates.setZero(((d * (d + 1)) / 2)+1, d);
    result.resize(((d * (d + 1)) / 2), u.cols());
    //create matrix where each column is a poin to compute the directional derivatives in the u and v direction

    for(unsigned j = 1; j < d+1; j++){
        coordinates(j,j-1) = 1;
    }

    for(int j = d+1; j < coordinates.rows(); j++){
        coordinates(j,(j-1)%coordinates.cols()) = 1;
        coordinates(j,j%coordinates.cols()) = 1;
    }

    if(debugme > 0){
        std::cout<<"coordinates for the directions are: "<<std::endl<<coordinates<<std::endl;
    }
    if(debugme > 5){
        std::cout<<"poits are "<<std::endl<<u<<std::endl;
    }
    directions.row(0) = coordinates.row(0);
    for(unsigned int j = 0; j < ((d * (d + 1)) / 2); j++){
        directions.row(1) = coordinates.row(j+1);
        if(debugme > 0){
            std::cout<<"directions for the second derivative are: "<<std::endl<<directions<<std::endl;
        }
        rThDerivSingle(2,i ,directions,u,tempresult);
        if(debugme > 5){
            std::cout<<"second derivatives are: "<<std::endl<<tempresult<<std::endl;
        }
        result.row(j) = tempresult.row(0);
    }
}

template<short_t d,class T>
void gsTriangularBezierBasis<d,T>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    result.setZero(m_compositions.size()*((d * (d + 1)) / 2), u.cols());
    gsMatrix<T> single_result;
    for(unsigned int i = 0; i < m_compositions.size(); i++)
    {
        deriv2Single_into(i, u, single_result);
        for(int j = 0; j < single_result.rows();j++){
            result.row(single_result.rows()*i+j) = single_result.row(j);
        }
    }
}




//private:
template<short_t d,class T>
gsVector<T>  gsTriangularBezierBasis<d,T>::getBaricentricCoordinate(gsVector<T>param)const{
    gsVector<T> result(3);

    gsMatrix<T> dd(3,3);
    dd(0,0) = 1;
    dd(0,1) = 1;
    dd(0,2) = 1;

    dd.row(1) = m_domain.col(0).transpose();
    dd.row(2) = m_domain.col(1).transpose();

    T det = dd.determinant();

    dd(1,0) = param[0];
    dd(2,0) = param[1];
    result[0] = dd.determinant()/det;

    dd(1,0) = m_domain(0,0);
    dd(2,0) = m_domain(0,1);
    dd(1,1) = param[0];
    dd(2,1) = param[1];
    result[1] = dd.determinant()/det;

    dd(1,1) = m_domain(1,0);
    dd(2,1) = m_domain(1,1);
    dd(1,2) = param[0];
    dd(2,2) = param[1];
    result[2] = dd.determinant()/det;
    return result;
}

template<short_t d,class T>
bool  gsTriangularBezierBasis<d,T>::isInDomain (gsVector<T> coord)const{
    for(int i = 0; i < coord.size(); i++)
        if(coord[i]<-0.00001)
            return false;
    //coord.array() < 0 ;
    return true;
}

template<short_t d,class T>
void gsTriangularBezierBasis<d,T>::getCompositions(std::vector<gsVector<unsigned int> > & compos)const{
    this->getCompositions(compos,m_degree);
}


template<short_t d,class T>
void gsTriangularBezierBasis<d,T>::getCompositions(std::vector<gsVector<unsigned int> > & compos, unsigned degree)const{
    int k = this->dim();//dimension
    int t = degree;
    //std::vector< gsVector<int> > result;
    gsVector<unsigned int,3> RR;
    RR.setZero();
    RR[0] = degree;
    compos.push_back(RR);
    while(RR[k]!=degree)
    {
        for (int i = 0; i <= k; i++)
        {
            if(RR[i]!=0)
            {
                t = RR[i];
                RR[i] = 0;
                RR[0] = t-1;
                RR[i+1] +=1;
                compos.push_back(RR);
                break;
            }
        }
    }
}

// based on derivatives from http://www.ann.jussieu.fr/~frey/papers/meshing/Farin%20G.,%20Triangular%20Bernstein%20Bezier%20patches.pdf
//page 11 equation 2.7
template<short_t d,class T>
void gsTriangularBezierBasis<d,T>::rThDerivSingle(unsigned int r, unsigned int i,const gsMatrix<T> & dir, const gsMatrix<T> & u, gsMatrix<T>& result)const{
    int debugme = -1;
    result.setZero(1,u.cols());//u,v directional derivative
    gsMatrix<T> temp_eval, temp1_eval;
    T sum_u = 0;
    gsMatrix<T> directions (3,1);//directions u and v for derivative
    gsTriangularBezierBasis<d,T> temp(r);
    std::vector<gsVector<unsigned int> > temp_compos;
    gsTriangularBezierBasis<d,T> temp1(m_degree-r);
    std::vector<gsVector<unsigned int> > temp1_compos;
    directions.col(0) = getBaricentricCoordinate(dir.row(1)) - getBaricentricCoordinate(dir.row(0));
    temp.getCompositions(temp_compos);//coposition for temp

    temp1.getCompositions(temp1_compos);//compositions for temp1
    for(int k = 0; k < u.cols(); k++){

        for(unsigned int j = 0; j < temp_compos.size();j++)
        {
            temp.evalSingle_into_baric(j,directions,temp_eval);

            int kk = findIndex(m_compositions[i]-temp_compos[j], temp1_compos);
            if(debugme>5){
                std::cout<<"m_composition[i] = "<<m_compositions[i].transpose()<<"\n  temp_compos[j] = "<< temp_compos[j].transpose()<<"\n   temp1_compos: "<<std::endl;//<<temp1_compos<<std::endl;
            }
            if(kk>=0){
                temp1.evalSingle_into(kk, u.col(k) , temp1_eval);
                sum_u += temp_eval(0,0)*temp1_eval(0,0);
            }
        }
        T cons = T(factorial(m_degree))/ T(factorial(m_degree-r));

        if(debugme>10){
            std::cout<< "constatnt: "<<cons<<std::endl;
        }
        result(0,k) = cons * sum_u;
        sum_u = 0;
    }
}

template<short_t d,class T>
int gsTriangularBezierBasis<d,T>::findIndex(gsVector<unsigned int> const p, std::vector<gsVector<unsigned int> >const compos )const{
    for(unsigned int i = 0; i < compos.size();i++){
        if(p==compos[i]){
            return i;
        }
    }
    return -1;
}

template<short_t d,class T>
void gsTriangularBezierBasis<d,T>::evalSingle_into_baric(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const{
    result.setZero(1,u.cols() );
    for(int j = 0; j < u.cols(); j++)
    {
        T denom = 1.0;
        T temp = 1.0;
        for(int k = 0; k < m_compositions[i].size(); k++)
        {
            denom = denom * factorial(m_compositions[i][k]);
            temp = temp * math::pow(u(k,j), static_cast<int>(m_compositions[i][k]));
        }
        result(0,j) =  (factorial(m_degree)/denom) * temp;
    }
}


}; // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTriangularBezierBasisXML.h)
#endif
