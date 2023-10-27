

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo
{

namespace internal
{

/// Get a Triangular Bezier from XML data
template<class T, short_t d>
class gsXml< gsTriangularBezier<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTriangularBezier<TMPLA2(d,T)>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "TriangularBezier"+to_string(d); }

    static gsTriangularBezier<d,T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsTriangularBezier<d,T> >(node);
    }

    static gsXmlNode * put (const gsTriangularBezier<d,T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml< gsTriangularBezier<d,T> >(obj,data);
    }
};


}// namespace internal

}// namespace gismo
