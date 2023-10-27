
#pragma once 

#include <gsCore/gsMemory.h>
#include <gsCore/gsBoundary.h>
#include <gsCore/gsGeometry.h>

namespace gismo
{

template <short_t d, class T>
gsMatrix<T> gsBoxSplineBasis<d,T>::support() const
{
    gsMatrix<T> supp(d,2);                                 // corners of the bounding box of the support -> 2 points in R^d columnwise!
    gsMatrix<T,d> m_fullX = getFullX(m_X, m_n);            // get full matrix to compute the support

    // compute the support by summing up all positive numbers separately and all negative numbers separately in each component
    for (unsigned i = 0; i < d; i++){
        T sumneg = 0;                                      // the sum of all negative numbers
        T sumpos = 0;
        for (int j = 0; j < m_fullX.cols(); j++){
           if ( m_fullX(i,j) < 0 )         sumneg += m_fullX(i,j);
           else if ( m_fullX(i,j) > 0 )    sumpos += m_fullX(i,j);
        }
        supp(i,0) = sumneg;
        supp(i,1) = sumpos;
    }
    return supp;
}

template <short_t d, class T>
gsMatrix<T> gsBoxSplineBasis<d,T>::support(const index_t & i) const
{
    GISMO_UNUSED(i);
    return support();
}

template <short_t d, class T>
gsMatrix<T,d> gsBoxSplineBasis<d,T>::getFullX(const gsMatrix<T,d> & X, const gsVector<int> & n) const
{
    gsMatrix<T,d> m_fullX ;                                         // constructing a matrix of ALL directions
    m_fullX.resize(d, n.sum() );

    int nrows = n.rows();
    int r = 0;
    for ( int i = 0; i < nrows; i++ ){                              // over all rows in m_n
        for ( int j = 0; j < n(i); j++ ){                           // add it as many times as it is its multiplicity
            for (unsigned k = 0; k < d; k++ ) m_fullX( k, r ) = X( k, i );
            r++;
        }
    }
    return m_fullX;
}

template <short_t d, class T>
gsMatrix<T,d> gsBoxSplineBasis<d,T>::getX(const gsMatrix<T,d> & X) const
{
    gsMatrix<T,d> m_distX ;                             // constructing a matrix of DISTINCT directions
    m_distX.resize(d, 1);
    m_distX.col(0) = X.col(0);                          // add the first vector to m_X from full X

    bool appear = 0;

    for (int i = 1; i < X.cols(); i++) {                // over all directions in X
        appear = 0;                                     // assume it did not appear in m_X yet
        for (int j = 0; j < m_distX.cols(); j++) {      // over the directions in the m_distX so far
            for (unsigned l = 0; l < d; l++){           // check all d components
                if ( X(l,i) != m_distX(l,j) ) break;    // if one component differs break and continue with the next vector in m_X
                if ( l == d-1 ) appear = 1;             // all components are the same => it is included
            }
        }
        if ( appear == 0 ) {                                                            // if it did not appear yet
            m_distX.conservativeResize(d, m_distX.cols() + 1);                          // add it to m_distX
            for (unsigned k = 0; k < d; k++) m_distX(k, m_distX.cols()-1) = X(k,i);
        }
    }
    return m_distX;
}

template <short_t d, class T>
gsVector<int> gsBoxSplineBasis<d,T>::getn(const gsMatrix<T,d> & X) const
{
    // we count the multiplicities of the different directions as they appear in X
    gsVector<int> m_mult ;
    m_mult.resize(m_X.cols());
    m_mult.setZero();
    m_mult(0) = 1;

    bool flag = 0;                                  // flag indicating the presence of the vector in m_X

    for (int i = 1; i < X.cols(); i++) {            // over all directions in X (full matrix)
        flag = 0;                                   // assuming it is not yet included in m_X
        for (int j = 0; j < m_X.cols(); j++) {      // over the directions in the m_X so far
            for (unsigned l = 0; l < d; l++) {      // over all d components
                if ( X(l,i) != m_X(l,j) ) break;    // if there is one component different, break and continue to next vector in m_X
                if ( l == d-1 ) {                   // all the components are the same => vector already exists in m_X
                    m_mult(j) += 1;                 // increase the multiplicity
                    flag = 1;                       // raise flag it is alreday included
                }
            }
            if ( flag == 1 ) break;                 // if it is already included break and continue to the next vector in full X
        }
    }
    return m_mult ;
}

template <short_t d, class T>
void gsBoxSplineBasis<d,T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    evalSingle_into(0, u, result);
}


template <short_t d, class T>
void gsBoxSplineBasis<d,T>::evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    result.resize(1, u.cols() );
    int ncols = m_X.cols();

    gsMatrix<T> Y;                                                  // the matrix for least squares
    Y = m_X.transpose() * (m_X * m_X.transpose()).inverse();

    gsVector<int> m(ncols);                                         // vector for the delayed translation - the position in the recursion tree
    m.setZero();                                                    // set to zero

    // compute the normals beforehand
    gsVector<T> M(ncols);
    M.setZero();
    int step = d-1;
    int k = ncols;
    gsMatrix<T,d> N = _computeBoxNormals(step,k,M);

    gsMatrix<T> J(1,u.cols());                                      // construct a matrix of ones for the matrix evaluation
    J.setOnes();
    gsMatrix<T> t = Y * u;                                          // t is a matrix, columnwise containing least square represenation for each u.col(i)
    result = _boxSplineRec( m_n, m, N, Y, t, J, u );                // recursion

}

template <short_t d, class T>
gsMatrix<T> gsBoxSplineBasis<d,T>::evalresult(const gsMatrix<T> & u) const
{
    gsMatrix<T> result(1, u.cols() );
    int ncols = m_X.cols();

    gsMatrix<T> Y;                                                  // the matrix for least squares
    Y = m_X.transpose() * (m_X * m_X.transpose()).inverse();

    gsVector<int> m(ncols);                                         // vector for the delayed translation - the position in the recursion tree
    m.setZero();                                                    // set to zero

    // compute the normals beforehand
    gsVector<T> M(ncols);
    M.setZero();
    int step = d-1;
    int k = ncols;
    gsMatrix<T,d> N = _computeBoxNormals(step,k,M);

    gsMatrix<T> J(1,u.cols());                                      // construct a matrix of ones for the matrix evaluation
    J.setOnes();
    gsMatrix<T> t = Y * u;                                          // t is a matrix, columnwise containing least square represenation for each u.col(i)
    result = _boxSplineRec( m_n, m, N, Y, t, J, u );                // recursion
    return result;
}

template <short_t d, class T>
gsMatrix<T> gsBoxSplineBasis<d,T>::derivativeresult(const gsMatrix<T> & u, const gsVector<int> & r) const
{
    //old implementation
    //gsBasis<T>::derivSingle_into(i, u, result);// default implementation

    //implementation by Dominik

    gsMatrix <T> result(1, u.cols() );

    gsMatrix<T> newu = u;
    gsVector<int> n = m_n;
    bool one_time_control = false;
    bool equal_check = false;
    for(int j = 0; j < u.cols(); j++)
    {
        for(int l = 0; l < r.rows(); l++)
        {
            newu(l,j) = newu(l,j) - r(l);
        }
    }
    for(int k = 0; k < m_X.cols(); k++)
    {
        equal_check = false;
        for(int l = 0; l < m_X.rows(); l++)
        {
            if(m_X(l,k) != r(l))
            {
                equal_check = true;
            }
        }
        if(!equal_check && one_time_control == false)
            n(k)=n(k)-1;
        else
            one_time_control = true;
    }
    //eval algorithm
    int ncols = m_X.cols();

    gsMatrix<T> Y;                                                  // the matrix for least squares
    Y = m_X.transpose() * (m_X * m_X.transpose()).inverse();

    gsVector<int> m(ncols);                                         // vector for the delayed translation - the position in the recursion tree
    m.setZero();                                                    // set to zero

    // compute the normals beforehand
    gsVector<T> M(ncols);
    M.setZero();
    int step = d-1;
    int k = ncols;
    gsMatrix<T,d> N = _computeBoxNormals(step,k,M);

    gsMatrix<T> J(1,u.cols());                                      // construct a matrix of ones for the matrix evaluation
    J.setOnes();
    gsMatrix<T> t = Y * u;                                          // t is a matrix, columnwise containing least square represenation for each u.col(i)
    gsMatrix<T> result1 = _boxSplineRec( n, m, N, Y, t, J, u );
    t = Y * newu;// recursion
    gsMatrix<T> result2 = _boxSplineRec( n, m, N, Y, t, J, newu );
    for(int j = 0; j < result1.cols(); j++)
    {
        for(int l = 0; l < result1.rows(); l++)
        {
            result(l,j)=result1(l,j)-result2(l,j);
        }
    }
    return result;
}

template <short_t d, class T> inline
void gsBoxSplineBasis<d,T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const 
{
    //implementation by Dominik
    derivSingle_into(0, u, result);

}

template <short_t d, class T>  inline
void gsBoxSplineBasis<d,T>::derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    //old implementation
    //gsBasis<T>::derivSingle_into(i, u, result);// default implementation

    //implementation by Dominik
    result.resize(1, u.cols() );

    //first for a fixed derivation r for testing
    gsVector<int> r(2);
    r << 1,0 ;

    gsMatrix<T> newu = u;
    gsVector<int> n = m_n;
    bool one_time_control = false;
    bool equal_check = false;
    for(int j = 0; j < u.cols(); j++)
    {
        for(int l = 0; l < r.rows(); l++)
        {
            newu(l,j) = newu(l,j) - r(l);
        }
    }
    for(int k = 0; k < m_X.cols(); k++)
    {
        equal_check = false;
        for(int l = 0; l < m_X.rows(); l++)
        {
            if(m_X(l,k) != r(l))
            {
                equal_check = true;
            }
        }
        if(!equal_check && one_time_control == false)
            n(k)=n(k)-1;
        else
            one_time_control = true;
    }
    //eval algorithm
    int ncols = m_X.cols();

    gsMatrix<T> Y;                                                  // the matrix for least squares
    Y = m_X.transpose() * (m_X * m_X.transpose()).inverse();

    gsVector<int> m(ncols);                                         // vector for the delayed translation - the position in the recursion tree
    m.setZero();                                                    // set to zero

    // compute the normals beforehand
    gsVector<T> M(ncols);
    M.setZero();
    int step = d-1;
    int k = ncols;
    gsMatrix<T,d> N = _computeBoxNormals(step,k,M);

    gsMatrix<T> J(1,u.cols());                                      // construct a matrix of ones for the matrix evaluation
    J.setOnes();
    gsMatrix<T> t = Y * u;                                          // t is a matrix, columnwise containing least square represenation for each u.col(i)
    gsMatrix<T> result1 = _boxSplineRec( n, m, N, Y, t, J, u );
    t = Y * newu;// recursion
    gsMatrix<T> result2 = _boxSplineRec( n, m, N, Y, t, J, newu );
    for(int j = 0; j < result1.cols(); j++)
    {
        for(int l = 0; l < result1.rows(); l++)
        {
            result(l,j)=result1(l,j)-result2(l,j);
        }
    }
}

template <short_t d, class T>
memory::unique_ptr<gsGeometry<T> > gsBoxSplineBasis<d,T>::makeGeometry(gsMatrix<T> coefs ) const
{
    return memory::unique_ptr<gsGeometry<T> >();
}

template <short_t d, class T>
void gsBoxSplineBasis<d,T>::uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots, int mul)
{
}


template <short_t d, class T>
void gsBoxSplineBasis<d,T>::uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, 
                                                       int numKnots, int mul)
{
    //gsSparseRows<T> trans;
}


// template <short_t d, class T>
// void gsBoxSplineBasis<d,T>::expandDirection(gsMatrix<T,d> & XX)
// {
//   XX.resize( d, m_mult.sum() );
//   XX.setOnes();
// }

// template <short_t d, class T>
// void gsBoxSplineBasis<d,T>::computeBoxNormals() { }


template <short_t d, class T>
std::ostream & 
gsBoxSplineBasis<d,T>::print(std::ostream &os) const
{
    os << "BoxSplineBasis function with \nfull direction matrix:\n" << this->getFullX(m_X,m_n)
       << "\nno-multi-direction matrix:\n" << m_X
       << "\nmultiplicity vector:\n" << m_n.transpose()
       << "\ndomain:\n"<< this->support();
    return os;
}

/* OLD RECURSION FUNCTION
 *
template <short_t d, class T>
T gsBoxSplineBasis<d,T>::_boxSplineRecBruteForce( const gsMatrix<T> & X, const gsMatrix<T> & u )
{
    assert( u.cols() == 1 ); // Just one point

    using std::abs;
    unsigned N = X.cols(); // number of all directions

    if ( N == d ) // Matrix is square, i.e. the base case
	{
	    const gsMatrix<T> Xip = X.inverse() * u;
	    for ( unsigned k= 0; k!=d; ++k)
        if ( (Xip(k,0) < 0) || (Xip(k,0) >= 1) )
		    return 0;
	    return T(1)/abs(X.determinant()) ;
	}
    else
	{
	    // Least norm vector s.t. Xz=u
        const gsMatrix<T> z = ( X.transpose() * ( X * X.transpose() ).inverse() ) * u;
	    //gsMatrix<T> z;
	    //z = X.jacobiSvd(Eigen::ComputeThinU | 
	    //                  Eigen::ComputeThinV).solve(u); 

	    T sum = 0;
        for ( unsigned k = 0; k != N; ++k )
		{
		    gsMatrix<T> Y(d,N-1);// X without X.col(k)
		    unsigned r = 0;
            for ( unsigned i = 0; i != N; ++i )
			if ( i == k )
			    continue;
			else
			    Y.col(r++) = X.col(i);

            if ( Y.colPivHouseholderQr().rank() < index_t(d) ) continue;
            else
                sum += z(k,0)    * box_spline_rec(Y,u          ) +
                    ( 1-z(k,0) ) * box_spline_rec(Y,u - X.col(k) ) ;
		}
        return sum / static_cast<T>(N-d) ;
	}
}
*/


// NEW RECURSION FUNCTION BY URSKA

template <short_t d, class T>
gsMatrix<T> gsBoxSplineBasis<d,T>::_boxSplineRec( const gsVector<int> & n,
                                                  const gsVector<int> & m,
                                                  const gsMatrix<T,d> & N,
                                                  const gsMatrix<T> & Y,
                                                  const gsMatrix<T> & t,
                                                  const gsMatrix<T,1> & J,
                                                  const gsMatrix<T> & u ) const
{
    using std::abs;

    // recursive case

    if (n.sum() > (index_t)(d))
    {
        gsMatrix<T> sum(1,u.cols());                        // sum in the recursion formula
        sum.setZero();
        int j = 0;                                          // for indicating position in t

        for ( int k = 0; k < m_X.cols(); ++k ) {            // over all distinct directions in m_X
            gsVector<int> nn = n ;
            nn(k) -= 1 ;                                    // reducing the multiplicity of the considered vector
            gsVector<int> mm = m ;
            mm(k) += 1 ;                                    // increasing the appearance of the considered vector in the delayed translation

            if (n(k) > 1) {                                 // if there are more vectors of this instance in remaining dir. matrix
                gsMatrix<T> res_1 = _boxSplineRec(nn, m , N, Y, t                       , J, u);    // recursion in the first summand
                gsMatrix<T> res_2 = _boxSplineRec(nn, mm, N, Y, t - (Y * m_X.col(k)) * J, J, u);    // recursion in the second summand

                for (int ll = 0; ll < u.cols(); ll++)
                    sum(0,ll) += t(j,ll) * res_1(0,ll) + ( n(k) - t(j,ll) ) * res_2(0,ll) ;
                j += 1 ;
            }
            else if (n(k) > 0) {                             // we reduced multiplicity of the vector to 0 => change matrix for least squares
                unsigned r = 0;
                gsMatrix<T,d> Z;
                for (int i = 0; i < m_X.cols(); i++) {       // rebuild a matrix for least squares, i.e., without the vector of which the multiplicity is reduced to 0
                    if (nn(i) > 0) {
                        Z.conservativeResize(d,r + 1);
                        Z.col(r++) = m_X.col(i);
                    }
                }
                if (Z.colPivHouseholderQr().rank() == (index_t)(d)) {                                             // if the matrix has full rank
                    gsMatrix<T> W = Z.transpose() * (Z * Z.transpose()).inverse();                              // new least squares matrix

                    gsMatrix<T> s = u;
                    for(int i = 0; i < m.rows(); i++) s -= static_cast<T>(m(i)) * (m_X.col(i) * J);             // for the second summand compute the delayed translation

                    gsMatrix<T> res_1 = _boxSplineRec(nn, m , N, W, W *  s                  , J, u);            // recursion in the first summand
                    gsMatrix<T> res_2 = _boxSplineRec(nn, mm, N, W, W * (s - m_X.col(k) * J), J, u);            // recursion in the second summand

                    for (int ll = 0; ll < u.cols(); ll++)
                       sum(0,ll) += t(j,ll) * res_1(0,ll) + ( n(k) - t(j,ll) ) * res_2(0,ll) ;
                }
                j += 1;
            }

        }
        return sum / static_cast<T>(n.sum() - d) ;
    }

    // the base case ( <=> matrix is square)
    else {
        gsMatrix<T> b(1,u.cols());
        b.setZero();
        unsigned r = 0;
        gsMatrix<T,d> Z;
        for (int i = 0; i < m_X.cols(); i++) {                                                  // rebuild a matrix
            if (n(i) > 0) {
                Z.conservativeResize(d,Z.cols()+1);
                Z.col(r++) = m_X.col(i);
            }
         }
        if (Z.colPivHouseholderQr().rank() == (index_t)(d)) {                                     // if the matrix has full rank
            b.setOnes();
            gsMatrix<T> s = u;
            for(int i = 0; i < m.rows(); i++) if (m(i) > 0) s -= static_cast<T>(m(i)) * (m_X.col(i) * J);

            // determine the indices of the remaining vectors from m_X
            gsVector<int> ind(d);
            int cind = 0;
            for (int l = 0; l < n.rows(); l++) if (n(l) > 0) ind(cind++) = l;

            // checking against all hyperplanes in the set of remaining vectors
            unsigned nb;
            for (unsigned ii = 0; ii < d; ii++) {

                // construct a bit number to find the normal vector
                nb = 0;
                for (int jj = 0; jj < n.rows(); jj++) 
                    if (jj != ind(ii) && n(jj) > 0) 
                        nb |= (1 << jj);

                gsVector<T> norm = N.col(nb);                              // normal vector to the hyperplane of included vectors wothout the p-th one
                T p = norm.dot(m_X.col(ind(ii)));                          // dot product with the p-th included vector

                T q;
                for (int kk = 0; kk < u.cols(); kk++){
                    q = norm.dot(s.col(kk));                               // dot product with the translation s = x - ksi_1 - ksi_2 - ...
                    b(0,kk) = math::min(b(0,kk), T(1.0) - 
                    ((p > 0 && q < 0) || (p < 0 && q >= 0)));
                }

                for (int kk = 0; kk < b.cols(); kk++){
                    q = norm.dot(s.col(kk) - m_X.col(ind(ii)));            // dot product with the translation s = x - ksi_1 - ksi_2 - ... - (p-th vector)
                    b(0,kk) = math::min(b(0,kk),T(1.0) - 
                    ((p > 0 && q >= 0) || (p < 0 && q < 0)));
                }
            }
            b = b / math::abs(Z.determinant()) ;
        }
        return b;
        }
}

template <short_t d, class T>
gsMatrix<T,d> gsBoxSplineBasis<d,T>::_computeBoxNormals(const int & r,
                                                        const int & k,
                                                        const gsVector<T> & M) const
{

    // construct a matrix that will store the normals; the size is (d x 2^k)
    int ncols = 1 << k;
    gsMatrix<T> N(d,ncols);
    N.setZero();

    // r tells how many vectors we still have to choose before reaching the base case
    // k tells which is the next direction to choose in m_X, going from the end to the beginning
    // if k < t then even if we reach the base case we will not have d-1 vectors for the hyperplane
    if (k >= r) {
        // not the base case yet, still have to choose
        if (r > 0) {
            // construct two matrices of half of the size of N and then combine them into N
            int ncolss = 1 << (k-1);

            // left one will contain normal vectors of the hyperplanes WITHOUT the k-th vector
            gsMatrix<T> N1(d,ncolss);
            N1 = _computeBoxNormals(r  ,k-1,M);

            // right one will contain normal vectors of the hyperplanes WITH the k-th vector
            gsMatrix<T> N2(d,ncolss);
            gsVector<T> M2 = M;
            M2(k-1) += 1;
            N2 = _computeBoxNormals(r-1,k-1,M2);

            N << N1,N2;

        }
        // base case, r = 0 => compute the normal vector to the hyperplane of the chosen vectors
        else {
            gsMatrix<T> S(d,d-1);
            S.setZero();
            int count = 0;
            for (int i = 0; i < m_X.cols(); i++) if (M(i) > 0) S.col(count++) = m_X.col(i); // fill S with vectors forming the hyperplane

            // works only for d = 2 or d = 3!!!
            switch (d) {
            case 2:
            {
                N.col(0) << -S(1,0),S(0,0);
                break;
            }
            case 3:
            {
                //N(0,0) = S(1,0)*S(2,1) - S(1,1)*S(2,0);
                //N(1,0) = S(2,0)*S(0,1) - S(0,0)*S(2,1);
                //N(2,0) = S(0,0)*S(1,1) - S(0,1)*S(1,0);
                N.col(0) = S.col(0).template segment<3>(0).cross( S.col(1).template segment<3>(0) );
                break;
            }
            default:
            {
                // Do a singular value decomposition on the matrix
                //Eigen::JacobiSVD<gismo::gsMatrix> svd(S, Eigen::ComputeFullV);

                // Get the V matrix
                //gsMatrix<T> V((int)svd.matrixV().rows(), (int)svd.matrixV().cols());
                //N.col(0) = svd.matrixV();

                break;
            }
            }

        }
    }
    return N;
}


} // namespace gismo
