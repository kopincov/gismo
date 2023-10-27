/**
   @file gismo-dev.h is the main header file to be included in
   applications using gismo-dev.
*/

#pragma once

#include <gsCore/gsDevConfig.h>
#include <gsCore/gsGeometryEvaluator.h>

#include <gsBoxSplines/gsBoxSplineBasis.h>

#include <gsSegment/gsVolumeSegment.h>

#include <gsBezier/gsBernsteinBasis.h>
#include <gsBezier/gsBezier.h>

#include <gsBezier/gsTensorBezier.h>
#include <gsBezier/gsTensorBernsteinBasis.h>

#include <gsPolynomial/gsMonomialBasis.h>

#include <gsIO/gsGnuplot.h>

#include <gsUtils/gsNorms.h>

#include <gsAssembler/gsNormL2.h>
#include <gsAssembler/gsSeminormH1.h>
#include <gsAssembler/gsVisitorNeumann2.h>

namespace {
const bool gismo_dev_filedata_ok = gismo::gsFileManager::addSearchPaths(GISMO_DEV_DATA_DIR);
}
