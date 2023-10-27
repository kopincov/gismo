

// stuff that does not fit here but I still do not know where to put it

#include <iostream>
#include <string>
#include <gsCore/gsLinearAlgebra.h>

#include "gsRemappedBasis.h"



namespace gismo {


GISMO_EXPORT void writeVTSdata     (std::ofstream &file, const std::string &name, const gsMatrix<real_t> &data);
GISMO_EXPORT void writeVTSpreamble      (std::ofstream &file, const gsMatrix<> &pts, const gsVector<unsigned> &np);
GISMO_EXPORT void writeVTSpostamble     (std::ofstream &file);
GISMO_EXPORT void printAllFunctions  (const gsRemappedBasis &space, const std::string &filename, size_t numPoints=1000);


} // gismo
