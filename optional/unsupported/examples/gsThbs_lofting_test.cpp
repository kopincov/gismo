/* Preparations for lofting with THB splines.
* Author: Dominik Mokris, dominik.mokris@jku.at
*/

#include <iostream>
#include <fstream>
#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;

int main()
{
    std::string filename = "thbs_00.xml";
    gsFileData<real_t> data(filename);
    gsTHBSpline<2> THBS;
    //Linking error: data.getFirst<gsTHBSpline<2> >(THBS);
    THBS = *data.getFirst< gsTHBSpline<2> >();

    gsInfo << THBS << "\n";

    gsMatrix<index_t> ll;
    gsMatrix<index_t> ur;
    gsVector<index_t> lvl;
    THBS.basis().tree().getBoxes(ll,ur,lvl);
    gsInfo << ll << ";\n" << ur <<";\n" << lvl << "\n";

    gsTensorBSpline<2,real_t > BSpline;
    THBS.convertToBSpline( BSpline );

    gsInfo << "AFTER" << "\n";
    //gsInfo << BSpline.coefs() ;
    gsInfo << BSpline.basis();

    gsFileData<> newdata;
    newdata << BSpline;
    newdata.dump( "thbs_00_converted_to_BSpline.xml" );

    return 0;
}
