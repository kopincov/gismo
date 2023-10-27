/** @file gsMotorIOUtils.h

    @brief ...

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Sajavicius, J. Speh
*/

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <gismo.h>
#include <gsIO/gsIOUtils.h>

namespace gismo
{

gsMultiPatch<>* readMultipatch(const std::string& inputFile);
void writeMultipatch(const gsMultiPatch<>& multipatch,
                     const std::string& outputFile,
                     unsigned npts = 1000);
gsGeometry<>* readGeometry(const std::string& inputFile);
gsGeometry<>* readGeometry(const std::string& inputFile,
                           const int id);
void writeGeometry(const gsGeometry<>& geom,
                   const std::string& outputFile);
void writeControlPoints(const gsMatrix<>& ctrlPts,
                        const std::string& outputFile);
//gsMatrix<> readMatrix(const std::string& inputFile);
void writeToTxtFile(const gsMatrix<>& X,
                    const gsMatrix<>& Y,
                    const std::string& outputFile);
void writeToTxtFile(const gsMatrix<>& X,
                    const std::string& outputFile);
void writeDisplacements(const gsMatrix<>& ctrlPts,
                        const gsMatrix<>& replacedCtrlPts,
                        const std::string& filename);
void writeMap(const gsGeometry<>& map,
              const std::string& outputPrefix);
void writeMap(const gsGeometry<>& map,
              const int id,
              const std::string& outputPrefix);
void writeKnots(const gsGeometry<>& map,
                const std::string& outputPrefix);
void writeKnots(const gsGeometry<>& map,
                const int id,
                const std::string& outputPrefix);

template <typename T>
void printMatrixSizes(const gsMatrix<T> &mat);

gsMultiPatch<>* readMultipatch(const std::string& inputFile)
{
    const gsFileData<> fd(inputFile);
    gsMultiPatch<>* mp = NULL;

    if ( fd.has< gsMultiPatch<> >() )
    {
        mp = fd.getAnyFirst< gsMultiPatch<> >().release();
    }
    else
    {
        gsInfo << "There is no multipatch in the file " << inputFile << std::endl;
    }

    return mp;
}

void writeMultipatch(const gsMultiPatch<>& multipatch,
                     const std::string& outputFile,
                     unsigned npts)
{
    gsInfo << "Writing " << outputFile << ".xml" << std::endl;
    gsFileData<> fd;
    fd << multipatch;
    fd.dump(outputFile);

    gsInfo << "Writing " << outputFile << ".pvd" << std::endl;
    gsWriteParaview(multipatch, outputFile, npts);
}

gsGeometry<>* readGeometry(const std::string& inputFile)
{
    gsFileData<> fd(inputFile);
    gsGeometry<>* geom = NULL;

    if ( fd.has< gsGeometry<> >() )
    {
        geom = fd.getAnyFirst< gsGeometry<> >().release();
    }
    else
    {
        gsInfo << "There is no geometry in the file " << inputFile << std::endl;
    }

    return geom;
}

gsGeometry<>* readGeometry(const std::string& inputFile,
                           const int id)
{
    gsFileData<> fd(inputFile);
    gsGeometry<>* geom = NULL;

    if ( fd.has< gsGeometry<> >() )
    {
        geom = fd.getId< gsGeometry<> >(id).release();
    }
    else
    {
        gsInfo << "There is no geometry in the file " << inputFile << std::endl;
    }

    return geom;
}


void writeGeometry(const gsGeometry<>& geom,
                   const std::string& outputFile)
{
    gsFileData<> fd;
    gsInfo << "Writing " << outputFile << ".xml" << std::endl;
    fd << geom;
    fd.dump(outputFile);
    gsInfo << "Writing " << outputFile << ".pvd" << std::endl;
    gsWriteParaview(geom, outputFile);
}

void writeControlPoints(const gsMatrix<>& ctrlPts,
                        const std::string& filename)
{
    const gsMatrix<>& points = ctrlPts.transpose();
    gsWriteParaviewPoints(points, filename);
}

//gsMatrix<> readMatrix(std::string& inputFile)
//{
//    std::ifstream input(inputFile);//("points.txt");
//    std::string line("");

//    index_t i = 0;
//    gsMatrix<> points(2, 2);
//    points.setZero();

//    while(getline(input, line)) {
//        if (!line.length() || line[0] == '#')
//        {
//            continue;
//        }
//        std::istringstream iss(line);
//        iss >> points(0, i) >> points(1, i);
//        i++;
//    }
//    return points;
//}

void writeToTxtFile(const gsMatrix<>& X,
                    const gsMatrix<>& Y,
                    const std::string& outputFile)
{
    std::ofstream file;

    gsInfo << "Writing " << outputFile << std::endl;
    file.open(outputFile.c_str());

    file << X;
    file << "\n";
    file << Y;
    file << "\n";
}

void writeToTxtFile(const gsMatrix<>& X,
                    const std::string& outputFile)
{
    std::ofstream file;

    gsInfo << "Writing " << outputFile << std::endl;
    file.open(outputFile.c_str());

    file << X;
    file << "\n";
}


void writeDisplacements(const gsMatrix<>& ctrlPts,
                        const gsMatrix<>& replacedCtrlPts,
                        const std::string& filename)
{
    gsKnotVector<> kv(0.0, 1.0, 0, 2);
    gsMultiPatch<> mp;

    for (index_t curve = 0; curve != ctrlPts.rows(); curve++)
    {
        gsMatrix<> coefs(2, 2);
        coefs.row(0) = ctrlPts.row(curve);
        coefs.row(1) = replacedCtrlPts.row(curve);

        mp.addPatch(gsBSpline<>(kv, coefs));
    }

    gsWriteParaview(mp, filename);

}

void writeMap(const gsGeometry<>& map,
              const std::string& outputPrefix)
{
    std::string out = outputPrefix + "_THBMap";
    //gsWriteParaview(map, out);
    gsWrite(map, out);
}


void writeMap(const gsGeometry<>& map,
              const int id,
              const std::string& outputPrefix)
{
    std::string out = outputPrefix + "_THBMap_" + util::to_string(id);
    gsWriteParaview(map, out);
    gsWrite(map, out);
}

void writeKnots(const gsGeometry<>& map,
                const std::string& outputPrefix)
{
    gsMesh<> mesh;
    makeMesh<>(map.basis(), mesh, 3);
    std::string out = outputPrefix + "_Knot_config";
    gsWriteParaview(mesh, out);

    map.evaluateMesh(mesh);
    out = outputPrefix + "_Knot_config_physical";
    gsWriteParaview(mesh, out);
}

void writeKnots(const gsGeometry<>& map,
                const int id,
                const std::string& outputPrefix)
{
    gsMesh<> mesh;
    makeMesh<>(map.basis(), mesh, 3);
    std::string out = outputPrefix + "_Knot_config_" + util::to_string(id);
    gsWriteParaview(mesh, out);

    map.evaluateMesh(mesh);
    out = outputPrefix + "_Knot_config_physical_" + util::to_string(id);
    gsWriteParaview(mesh, out);
}

template <typename T>
void printMatrixSizes(const gsMatrix<T>& mat)
{
    gsInfo << mat.rows() << " x " << mat.cols() << std::endl;
}

} // namespace gismo

