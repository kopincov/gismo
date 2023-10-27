/** @file tutorial_topology.cpp

    @brief Prints a topology of multi patch. 

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    std::string multipatch_file = MOTOR_DATA_DIR "jku/triangle_topology.xml";
    std::string output = "topology";
    
    gsCmdLine cmd("Outputs a topology from a given multi-patch");
    cmd.addPlainString("filename", "File containing a multi-patch", multipatch_file);
    cmd.addString("b", "boundary", "File containing boundary data", output);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Input Arguments: \n"
           << "Multipatch: " << multipatch_file << "\n"
           << "Output:     " << output << "\n"
           << "---------------------------------------------------------\n\n";

    gsFileData<> fd(multipatch_file);
    gsMultiPatch<>* multipatch = fd.getAnyFirst< gsMultiPatch<> >().release();
    
    for (gsMultiPatch<>::const_biterator it = multipatch->bBegin();
         it != multipatch->bEnd(); ++it)
    {
        const gsGeometry<>& geometry = multipatch->patch(it->patch);
        const gsGeometry<>::uPtr boundary = geometry.boundary(it->side());
	    
        std::stringstream ss;
        ss << output << "_boundary__p_" << util::to_string(it->patch)
           << "_" << it->side();
	
        std::string out = ss.str();
        gsInfo << "Writing: " << out << "\n";
        gsWriteParaview(*boundary, out);
    }

    for (gsMultiPatch<>::const_iiterator it = multipatch->iBegin();
         it != multipatch->iEnd(); ++it)
    {
        const gsGeometry<>& geometry = multipatch->patch(it->first().patch);
        const gsGeometry<>::uPtr interface = geometry.boundary(it->first().side());
	
        std::stringstream ss;
        ss << output << "_interface___p_" << util::to_string(it->first().patch)
           << "_" <<  it->first().side() << "___p_" << util::to_string(it->second().patch)
           << "_" << it->second().side();
	
        std::string out = ss.str();
        gsInfo << "Writing: " << out << "\n";
        gsWriteParaview(*interface, out);
    }

    return 0;
}


