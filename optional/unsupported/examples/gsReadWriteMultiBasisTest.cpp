/** @file gsReadWriteMultiBasis.cpp

    @brief A test for reading/writing a multi basis from/to xml.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/
#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;

void readingCommandLineArguments(int argc,
                                 char* argv[],
                                 std::string& multiBasisFile)
{
    multiBasisFile = "basis3d/multiBasis.xml";

    gsCmdLine cmd("Testing the reading/writing of xml files of multi basis");

    cmd.addString("i", "input", "Input xml file containing multi basis", multiBasisFile);

    cmd.getValues(argc,argv);

    gsInfo << "---------------------------------------------------------\n\n"
              "Input Arguments: \n\n"
              "Input file: " << multiBasisFile << "\n\n"
              "---------------------------------------------------------\n\n";
}

bool areEqual(const gsBasis<>& basis1,
              const gsBasis<>& basis2)
{
    if (basis1.dim() != basis2.dim())
    {
        return false;
    }

    if (basis1.size() != basis2.size())
    {
        return false;
    }

    return true;
}

bool areEqual(const std::vector<patchSide> boundaries1,
              const std::vector<patchSide> boundaries2)
{
    for (size_t i = 0; i != boundaries1.size(); i++)
    {
        if (boundaries1[i] != boundaries2[i])
        {
            return false;
        }
    }
    return true;
}

bool areEqual(const std::vector<boundaryInterface> interfaces1,
              const std::vector<boundaryInterface> interfaces2)
{
    for (size_t i = 0; i != interfaces1.size(); i++)
    {
        if ( !(interfaces1[i] == interfaces2[i]) )
        {
            return false;
        }
    }

    return true;
}

bool areEqual(const gsBoxTopology& topology1,
              const gsBoxTopology& topology2)
{
    if (topology1.nInterfaces() != topology2.nInterfaces())
    {
        return false;
    }

    if (topology1.nBoundary() != topology2.nBoundary())
    {
        return false;
    }

    if (!areEqual(topology1.boundaries(), topology2.boundaries()))
    {
        return false;
    }

    if (!areEqual(topology1.interfaces(), topology2.interfaces()))
    {
        return false;
    }

    return true;
}


bool areEqual(const gsMultiBasis<>& mb1,
              const gsMultiBasis<>& mb2)
{
    if (mb1.nBases() != mb2.nBases())
    {
        return false;
    }

    for (size_t index = 0; index != mb1.nBases(); index++)
    {
        if (!areEqual(mb1.basis(index), mb2.basis(index)))
        {
            return false;
        }
    }

    if (!areEqual(mb1.topology(), mb2.topology()))
    {
        return false;
    }

    return true;
}


int main(int argc, char* argv[])
{
    std::string multiBasisFile;

    try { readingCommandLineArguments(argc, argv, multiBasisFile);  } catch (int rv) { return rv; }

    // read the multi basis
    gsFileData<> fd1(multiBasisFile);
    gsMultiBasis<>::uPtr mb1 = fd1.getAnyFirst< gsMultiBasis<> >();

    // write it back to file
    gsFileData<> fd2;
    fd2 << *mb1;
    fd2.dump("some_file");

    // read it again
    gsFileData<> fd3("some_file.xml");
    gsMultiBasis<>::uPtr mb3 = fd3.getAnyFirst< gsMultiBasis<> >();

    int failed;
    if (areEqual(*mb1, *mb3))
    {
        gsInfo << "Multi basis are equal" << "\n";
        failed = 0;
    }
    else
    {
        gsInfo << "Multi basis are not equal" << "\n";
        failed = 1;
    }

    return failed;
}
