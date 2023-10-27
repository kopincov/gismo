#pragma once

#include <gsCore/gsGeometry.h>
#include <gsPolynomial/gsMonomialBasis.h>

namespace gismo
{

/** \brief
    An univariate polynomial in monomial basis.

    This is the geometry type associated with gsMonomialBasis.

    \tparam T coefficient type

    \ingroup geometry
*/

// replaces appeareances of \a oldStr with \a newStr inside the string
// \a str
inline void stringReplace(std::string& str, 
                          const std::string& oldStr, 
                          const std::string& newStr)
{
    size_t pos = 0;
    while((pos = str.find(oldStr, pos)) != std::string::npos)
    {
        str.replace(pos, oldStr.length(), newStr);
        pos += newStr.length();
    }
}


template<class T>
class gsMonomialPoly : public gsGeoTraits<1,T>::GeometryBase
{
public:
    typedef gsMonomialBasis<T> Basis;

    typedef typename gsGeoTraits<1,T>::GeometryBase Base;

    /// Shared pointer for gsMonomialPoly
    typedef memory::shared_ptr< gsMonomialPoly > Ptr;

    /// Unique pointer for gsMonomialPoly
    typedef memory::unique_ptr< gsMonomialPoly > uPtr;

    using Base::m_coefs;
    
    // Default constructor
    //gsMonomialPoly() {}

    /// Constructs a monomial polynomial by coefficient matrix. degree
    /// is induced by the size of coefficients
    explicit gsMonomialPoly(const gsMatrix<T> & coefs, int p = -1) :
    Base(Basis(p==-1?coefs.rows()-1:p), coefs )
    { }
    
    explicit gsMonomialPoly(const std::string & str, std::string var = "x")
    {this->set_str(str,var);}

    /// Constructs a monomial polynomial by basis and coefficient matrix
    gsMonomialPoly(const Basis & basis, const gsMatrix<T> & coefs ) :
    Base( basis, coefs )
    { }
    
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        const int p = deg();
        if (isConstant())
        {
            os << "gsMonomialPoly (constant): "<< m_coefs.at(0);
            return os;
        }
        else
        {
            os << "gsMonomialPoly (deg="<<p<<"): ";
            for (index_t i = p; i!=0; --i)
            {
                if ( 0 != m_coefs.at(i) )
                {
                    if (i!=p) os << (m_coefs.at(i) > 0 ? "+" : "-");
                    const T a = math::abs(m_coefs.at(i)); 
                    if (1!=a) os << a <<"*";
                    if (1!=i) os<<"x^"<<i; else os<<"x";
                }
            }
            if ( 0 != m_coefs.at(0) )
                os <<(m_coefs.at(0)>0?"+":"")<<m_coefs.at(0);
            return os;
        }
    }

    GISMO_BASIS_ACCESSORS
    
    GISMO_CLONE_FUNCTION(gsMonomialPoly)

    int deg() const { return basis().deg();}
    
    /// returns true iff the polynomial is equal to zero
	bool isNull() const
    {
        return (m_coefs.array() == 0).all();
    }

    /// returns true iff the polynomial is constant
	bool isConstant() const
    {
        return (m_coefs.bottomRows(deg()).array() == 0).all();
    }

    /// returns true iff the polynomial is monic
	bool isMonic() const
    {
        return (1 == leadCoeff()).all();
    }

    /// returns the leading coefficient
	T leadCoeff()
    {
        GISMO_ASSERT(1==m_coefs.cols(), "to do");
        GISMO_ASSERT(!isNull(), "Poly is zero");
        int lead = deg();// right-most non-zero
        while ( 0==m_coefs.at(lead) ) --lead;
        return m_coefs.at(lead);
    }

    /// returns the trailing coefficient
	T trailCoeff()
    {
        GISMO_ASSERT(1==m_coefs.cols(), "to do");
        return m_coefs.at(0);
    }

    
    /**
       \brief Performs conversion from monomial basis to Bernstein
       basis (on [0,1]).
    */
    void asBezier(gsBezier<T> & bezier) const // convert functions ?
    {
        const int p = this->deg();
        const int n = this->geoDim();
        gsBernsteinBasis<T> bbasis(0, 1, p);
        
        gsMatrix<T> diff_table(p+1,p+1); // difference table
        gsMatrix<T> coefs_bezier(bbasis.size(), n);
        
        // Loop over the dimension of the coefficients
        const gsMatrix<T> & coefs = m_coefs;
        for(index_t j=0; j!=n; j++)
        {
            // Load monomial coefficients into the left column and scale them
            diff_table.col(0)=coefs.col(j);
            for(int i=1; i<=p-1; i++)
                diff_table(i,0) /= binomial(p, i);
            
            // Compute difference table backwards
            for(int column=p-1; column>=0; column--)
                for(int row=1; row<=p-column; row++)
                    diff_table(column, row)=diff_table(column+1,row-1)+diff_table(column,row-1);
            
            // Extract Bernstein coefficients from the the top row
            coefs_bezier.col(j)=diff_table.row(0).transpose();
        }
        
        // Generate Bernstein polynomial with the computed coefficients
        bezier = gsBezier<T>(bbasis, give(coefs_bezier));
    }

protected:

    void set_str(const std::string & str, std::string var = "x");
};

template<class T>
void gsMonomialPoly<T>::set_str(const std::string & str, std::string var)
{
    // step 1. Normalize
    std::string poly(" ");
    for ( std::string::const_iterator it=str.begin(); it!=str.end(); ++it)
    {
        if (*it==var[0]) // var=char..
        {
            if (*poly.rbegin() != '*' ) poly += "1*";
            poly += "x";
            if ( (it+1==str.end() || *(it+1)!='^') ) poly += "^1";
        }
        else if (*it=='+') poly += " ";
        else if (*it=='-') poly += " -";
        else if (*it!=' ') poly += *it;
    }
    
    // step 2. Read in
    std::string t;
    std::istringstream term, in(poly);
    T cf;
    std::vector<T> vcf(1,0);
    unsigned ex;
    std::vector<unsigned> vex(1,0);
    
    while (in >> t)
    {
        //gsInfo << "|"<< t <<"|\n";
        if (t.find("x")!=std::string::npos)
        {
            stringReplace(t, "*x^", " ");
            term.clear();term.str(t);
            if (!gsGetValue(term, cf)) gsWarn<<"Error parsing coefficient.\n";
            if (!gsGetInt  (term, ex)) gsWarn<<"Error parsing exponent.\n";
            vcf.push_back(cf);
            vex.push_back(ex);
        }
        else
        {
            term.clear();term.str(t);
            if (!gsGetValue(term, cf)) gsWarn<<"Error parsing constant coefficient.\n";
            vcf.push_back(cf);
            vex.push_back(0);
        }
    }
    
    // step 3. Write polynomial coefficients
    GISMO_ASSERT(NULL==this->m_basis, "to do");
    this->m_basis = new Basis(*std::max_element(vex.begin(), vex.end()));
    m_coefs.setZero(deg()+1, 1);
    for (size_t i = 0; i!= vex.size(); ++i)
        m_coefs.at( vex[i] ) += vcf[i];
}

} // namespace gismo
