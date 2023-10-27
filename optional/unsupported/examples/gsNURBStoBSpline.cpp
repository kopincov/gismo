/** @file gsNURBStoBSpline.cpp

    @brief Converts geometry patches from other formats to
    tensor-prdocut B-Spline format, approximately

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string fn("quarter_annulus.xml");
    gsCmdLine cmd("Translate NURBS 2 BSpline.");    

    cmd.addPlainString("filename", "File containing data to convert", fn);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr geo = gsReadFile<>(fn);

    if ( geo )
    {
        gsInfo<<"Got "<< *geo <<"\n";
    }
    else
    {
        gsInfo<<"No input. Quitting\n";
        return 0;
    }

    gsMatrix<> intPoints, values;
    gsFileData<> fd;

    gsGeometry<>::uPtr g = gsGeometry<>::uPtr();

    for (size_t n = 0; n != geo->nPatches(); ++n)
    {

        if ( gsTensorNurbs<2>* tn = dynamic_cast< gsTensorNurbs<2>*>(&geo->patch(n)) )
        {
            const gsTensorBSplineBasis<2> & tb = tn->basis().source();

            tb.anchors_into(intPoints);
            tn->eval_into(intPoints, values);

            g = tb.gsBasis<>::interpolateData(values, intPoints);

            fd << *g ;
            gsInfo<<"Converted gsTensorNurbs<2> to BSpline.\n";
        }

    
        if ( gsTensorNurbs<3>* tn = dynamic_cast< gsTensorNurbs<3>*>(&geo->patch(n)) )
        {
            const gsTensorBSplineBasis<3> & tb = tn->basis().source();

            tb.anchors_into(intPoints);
            tn->eval_into(intPoints, values);

            g = tb.gsBasis<>::interpolateData(values, intPoints);

            fd << *g ;
            gsInfo<<"Converted gsTensorNurbs<3> to BSpline.\n";
        }

        if ( gsHTensorBasis<2>* tn = dynamic_cast< gsHTensorBasis<2>*>(&geo->patch(n).basis()) )
        {
            const gsTensorBSplineBasis<2,real_t> & tb = tn->tensorLevel(0);

            tb.anchors_into(intPoints);
            geo->patch(n).eval_into(intPoints, values);

            g = tb.gsBasis<>::interpolateData(values, intPoints);

            fd << *g ;
            gsInfo<<"Converted gsHTensorBasis<2> to BSpline.\n";
        }

        if ( gsHTensorBasis<3>* tn = dynamic_cast< gsHTensorBasis<3>*>(&geo->patch(n).basis()) )
        {
            const gsTensorBSplineBasis<3,real_t > & tb = tn->tensorLevel(0);

            tb.anchors_into(intPoints);
            geo->patch(n).eval_into(intPoints, values);

            g = tb.gsBasis<>::interpolateData(values, intPoints);

            gsInfo<<"Converted to "<<*g<<"\n";

            fd << *g ;
            gsInfo<<"Converted gsHTensorBasis<3> to BSpline.\n";
        }
    
    }// end for all patches

    if ( fd.numData() > 0)
        fd.save("convertedtoBSpline");
    else
        gsInfo<<"Did not find a geometry to be converted in the file.\n";

    return 0;
}
