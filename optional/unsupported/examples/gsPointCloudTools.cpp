/** @file gsPointCloudTools.h

    @brief It allows to convert from .par to Gismo's .xml format.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <gsFitting/gsPointCloudTools.h>
#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string input("");
    std::string output("");
    
    index_t flag   = 0; // User can decide which function to run
    
    real_t u_min    = 0;
    real_t u_max    = 1;
    real_t v_min    = 0;
    real_t v_max    = 1;
    real_t left     = 0;
    real_t right    = 1;
    bool dir        = true;
    real_t scale    = 1;
    

    // ----------------------------------------------------------------------

    gsCmdLine cmd("Converting from .par to .xml.");
    
    cmd.addString("i", "input", ".par file with parameters and points", input);
    cmd.addString("o", "output", "G+Smo .xml file having parameters and points", output);
    
    cmd.addInt("f", "flag", "Function to run: \n"
					     "0) points selection and scaling \n"
					     "1) file conversion from .xml to .par \n"
					     "2) file conversion from .par to .xml \n"
					     "3) Scaling of the Parameters u and v give in input in a .par file\n"
                         "4) Getting rid of values which are in a given interval\n"
                         "5) Transpose the parametric domain", flag);
    
    cmd.addReal("u", "u_min", "minimum u value (optional, default: 0)", u_min);
    cmd.addReal("U", "u_max", "maximum u value (optional, default: 1)", u_max);
    cmd.addReal("v", "v_min", "minimum v value (optional, default: 0)", v_min);
    cmd.addReal("V", "v_max", "maximum v value (optional, default: 1)", v_max);
    cmd.addReal("l", "left", "left extreme of the interval to get rid of: [left, right]", left);
    cmd.addReal("r", "right", "right extreme of the interval to get rid of: [left, right]", right);
    cmd.addSwitch("d", "direction", "direction of the interval to get rid of: u = 0, v = 1", dir);
    
    cmd.addReal("s", "scale", "scaling factor for length conversion", scale);
    
     
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    // ----------------------------------------------------------------------
    
    gsInfo << "----------------------------------------\n"
           << "INPUT DATA:\n\n" 
           << "input:                " << input << "\n\n"
	   << "output:               " << output << "\n\n"
           << "----------------------------------------" << std::endl;
    
    // ----------------------------------------------------------------------
        
    if (input == "")
    {
        gsInfo << "No input. Exiting...\n";
        return 0;
    }

    //restrictedPar(input, output, 0.01, 0.99, 0, 0.25, true, 1./1000); // fig.2(a) or fig.7 paper Adaptative fitting 
    //restrictedPar(input, output, 0.3, 0.7, 0, 0.25, true, 1); // fig.2(b) or fig.8 paper Adaptative fitting
    
    switch(flag){
        case 5:
            transpose_domain(input, output);
            break;
        case 4:
            cutting_interval(input, output, left, right, dir);
            gsInfo << "direction: " << dir << std::endl;
            break;
        case 3:
            //scaled_Parameters(input, output);
            opt_scale(input, output);
            break;
        case 2:
            parToXML(input, output);
            break;
        case 1:
            XMLToPar(input, output);
            break;
        case 0:
            restrictedPar(input, output, u_min, u_max, v_min, v_max, true, scale);
            break;
        default:
            std::cout << "Invalid function choice. " << "\n";
    }
  
  
  
  
    return 0;
}


