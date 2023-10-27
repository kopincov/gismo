/** @file gsTriangularBezier.h

    @brief Provides declaration of TriangularBezier class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): G. Kiss
*/

#pragma once
#include <gsUtils/gsPointGrid.h>
#include <algorithm>
#include <gsCore/gsGeometry.h>
#include <gsUtils/gsMesh/gsMesh.h>
#include <gsBezier/gsTriangularBezierBasis.h>

namespace gismo
{

/** \brief
    A triangular Bezier function, in \em d dimensions.

    This is the geometry type associated with gsTriangularBezierBasis.

    \tparam T is the coefficient type

    \ingroup geometry
*/

template<short_t d, class T>
class gsTriangularBezier : public gsGeoTraits<d,T>::GeometryBase
{
public:
    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    typedef gsTriangularBezierBasis<d,T> Basis;

    /// Shared pointer for gsTriangularBezier
    typedef memory::shared_ptr< gsTriangularBezier<2,T> > Ptr;

    /// Unique pointer for gsTriangularBezier
    typedef memory::unique_ptr< gsTriangularBezier > uPtr;

public:

    /// Construct Triangular Bezie patch from basis functions and coefficient matrix
    gsTriangularBezier( const Basis * basis, const gsMatrix<T> * coefs ) :
    Base( basis, coefs ) { }

    /// Construct Triangular Bezie patch from basis functions and coefficient matrix
    gsTriangularBezier( const Basis & basis, const gsMatrix<T> & coefs ) :
    Base( basis, coefs ) { }

    /// Copy constructor
    gsTriangularBezier( const gsTriangularBezier & other )
    {
        this->m_basis = other.basis().clone().release();
        this->m_coefs = other.coefs();
    }

    GISMO_CLONE_FUNCTION(gsTriangularBezier)

    GISMO_BASIS_ACCESSORS

public:

    //////////////////////////////////////////////////
    // Virtual member functions required by the base class
    //////////////////////////////////////////////////

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "Triangular Bezier.\n"; return os; };

    /// Returns the degree wrt direction i
    short_t degree(short_t i = 0) const
    { return this->basis().component(i).degree(i); };

    //TODO rename to eval_into
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    //subdivision of the triangular patch into 3 patches at a parameter value defined in u

    void subdivide(const gsVector<T>& u, gsMultiPatch<T> & result);


    void toMesh(gsMesh<T> & msh, int npoints) const;
//    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
}; // class gsTriangularBezier


}; // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////




namespace gismo
{

template<short_t d, class T>
void gsTriangularBezier<d,T>::toMesh(gsMesh<T> & msh, int npoints) const
{
    // compute number of coordinates in each axis
    const gsMatrix<T> param = this->parameterRange();
    const gsVector<T> a = param.col(0);
    const gsVector<T> b = param.col(1);
    gsVector<unsigned> np    = uniformSampleCount(a, b, npoints );
    np[0] = std::max(np[0],np[1]);

    std::vector<gsVector<unsigned int> > compos;
    this->basis().getCompositions(compos, np[0]);
    std::vector<gsVector<T> > baric;
    gsVector<T,3> temp;
    gsMatrix<T> cp;
    for(unsigned int i = 0; i < compos.size();i++){
        temp[0] = T(compos[i][0])/np[0];
        temp[1] = T(compos[i][1])/np[0];
        temp[2] = T(compos[i][2])/np[0];
        baric.push_back(temp);
    }

    gsMatrix<T,3,2> domain = this->basis().get_domain();
    gsVector<T,2>eucCoord;
    for(unsigned int i = 0; i < baric.size();i++){
        eucCoord[0] = baric[i][0] * domain(0,0) + baric[i][1] * domain(1,0) + baric[i][2] * domain(2,0);
        eucCoord[1] = baric[i][1] * domain(0,1) + baric[i][1] * domain(1,1) + baric[i][2] * domain(2,1);
        this->eval_into( eucCoord,cp);
        msh.addVertex( cp );
    }

    //insert all needed triangles joining the vertices
    int max = np[0];
    int shift = max;
    int rcount = 0;
    int iter = 0;//number of trianges used in the jth iteration of the deCast.alg.
    for(unsigned int k = 0; k <= np[0]+1 ;k++)
        iter += k;
    for(int k = 0; k < iter-2; k++){
        if(k!=max){
            msh.addFace(k  , k+1, k+shift+1);
            if((rcount>0) && (rcount<shift)){
                msh.addFace(k  , k+shift, k+shift+1);
            }
            rcount++;
        }else{
            max +=shift;
            shift--;
            rcount = 0;
        }
    }
}

template<short_t d, class T>
void gsTriangularBezier<d,T>:: eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const{
    result.resize(3,u.cols());//every column is an evaluated point
    for(int i = 0; i < u.cols(); i++){
        gsMatrix<T>  coef = this->m_coefs;
        gsVector<> bar = this->basis().getBaricentricCoordinate(u.col(i));//baricentric coordinates of a parameter value
        for(int j = 0; j < this->basis().degree(); j++){ //start of the deCastel. alg.
            int max = this->basis().degree() - j;
            int shift = max;
            int iter = 0;//number of trianges used in the jth iteration of the deCast.alg.
            int pos = 0;//position in the coefficient matrix
            for(int k = 0; k <= this->basis().degree()-j+1 ;k++)
                iter += k;
            for(int k = 0; k < iter-2; k++){
                if(k!=max){
                    coef(pos,0) = bar[0]*coef(k,0) + bar[1]*coef(k+1,0) +
                            bar[2]*coef(k+1+shift,0);//x-coordinate
                    coef(pos,1) = bar[0]*coef(k,1) + bar[1]*coef(k+1,1) +
                            bar[2]*coef(k+1+shift,1);//y-coordinate
                    coef(pos,2) = bar[0]*coef(k,2) + bar[1]*coef(k+1,2) +
                            bar[2]*coef(k+1+shift,2);//z-coordinate
                    pos++;
                }else{
                    max +=shift;
                    shift--;
                }
            }
        }
        result(0,i) =coef(0,0);
        result(1,i) =coef(0,1);
        result(2,i) =coef(0,2);
    }
}

template<short_t d, class T>
void gsTriangularBezier<d,T>::subdivide(const gsVector<T>& u, gsMultiPatch<T> & result)
{
    gsMatrix<T> g1, g2, g3;//coefficients for the resulting geometries
    g1.resize(this->m_coefs.rows(),this->m_coefs.cols());
    g2.resize(this->m_coefs.rows(),this->m_coefs.cols());
    g3.resize(this->m_coefs.rows(),this->m_coefs.cols());
    int pos = 0;//the index of the next CP to be added to g1,g2,g3
    gsMatrix<T>  coef = this->m_coefs;
    gsVector<> bar = this->basis().getBaricentricCoordinate(u);
    for(int j = 0; j < this->basis().degree(); j++){
        //here copy the control points to result
        int shiftG2 = this->basis().degree()+1-j;
        int shiftG3 = this->basis().degree()-j;
        int posG2 = shiftG2;
        int posG3 = shiftG3;

        g1.row(pos) = coef.row(0);
        g2.row(pos) = coef.row(0);
        g3.row(pos) = coef.row(posG3);
        posG3 += shiftG3;
        //shiftG3--;
        pos++;
        for(int i = 1; i <  this->basis().degree() - j+1;i++){
            g1.row(pos) = coef.row(i);
            g2.row(pos) = coef.row(posG2);
            g3.row(pos) = coef.row(posG3);
            shiftG2--;
            shiftG3--;
            posG2 += shiftG2;
            posG3 += shiftG3;
            pos++;

        }
        int max = this->basis().degree() - j;
        int shift = max;
        int iter = 0;
        int pos2 = 0;
        for(int k = 0; k <= this->basis().degree()-j+1 ;k++)
            iter += k;
        for(int k = 0; k < iter-2; k++){
            if(k!=max){
                coef(pos2,0) = bar[0]*coef(k,0) + bar[1]*coef(k+1,0) +
                        bar[2]*coef(k+1+shift,0);
                coef(pos2,1) = bar[0]*coef(k,1) + bar[1]*coef(k+1,1) +
                        bar[2]*coef(k+1+shift,1);
                coef(pos2,2) = bar[0]*coef(k,2) + bar[1]*coef(k+1,2) +
                        bar[2]*coef(k+1+shift,2);
                pos2++;
            }else{
                max +=shift;
                shift--;
            }
        }
    }
    //adding the result from the last iteration of the deCastl.alg.
    g1.row(pos) = coef.row(0);
    g2.row(pos) = coef.row(0);
    g3.row(pos) = coef.row(0);
    gsMatrix<T,3,2> domain = this->basis().get_domain();
    domain.row(2) = u;
    gsTriangularBezierBasis<2,T> basis_g1(this->basis().degree(), domain);
    domain.row(1) = this->basis().get_domain().row(2);
    gsTriangularBezierBasis<2,T> basis_g2(this->basis().degree(), domain);
    domain.row(0) = this->basis().get_domain().row(1);
    gsTriangularBezierBasis<2,T> basis_g3(this->basis().degree(), domain);

    typename gsTriangularBezier<2,T>::uPtr geom_g1(new gsTriangularBezier<2,T>(basis_g1,g1));
    typename gsTriangularBezier<2,T>::uPtr geom_g2(new gsTriangularBezier<2,T>(basis_g2,g2));
    typename gsTriangularBezier<2,T>::uPtr geom_g3(new gsTriangularBezier<2,T>(basis_g3,g3));
    result.addPatch(give(geom_g1));
    result.addPatch(give(geom_g2));
    result.addPatch(give(geom_g3));
}

}; // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTriangularBezierXML.h)
#endif
