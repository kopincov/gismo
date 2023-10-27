
#include <gsIO/gsXmlUtils.h>

namespace gismo
{

namespace internal
{


/// Get a TriangularBezier Basisfrom XML data
template<short_t d,class T>
class gsXml< gsTriangularBezierBasis<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTriangularBezierBasis<TMPLA2(d,T)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "gsTriangularBezierBasis"+ (d>1 ? to_string(d):""); }

    static gsTriangularBezierBasis<d,T> * get (gsXmlNode * node)
    {

        // Read degree
        unsigned deg = atoi( node->first_attribute("degree")->value() );
        
        gsTriangularBezierBasis<d,T> * TrBB = new gsTriangularBezierBasis<d,T>(deg);
        return TrBB;
    }

    static gsXmlNode * put (const gsTriangularBezierBasis<d,T> & obj,
                            gsXmlTree & data )
    {
        // Add a new node (without data)
        gsXmlNode* tp_node = internal::makeNode("Basis" , data);
        
        tp_node->append_attribute( makeAttribute("type",
                                                 internal::gsXml<gsTriangularBezierBasis<d,T> >::type().c_str(), data) );
        
        std::ostringstream str;
        str<< obj.get_degree();
        tp_node->append_attribute( makeAttribute( "degree", str.str(),data ));
        return tp_node;
    }
};


}// namespace internal

}// namespace gismo
