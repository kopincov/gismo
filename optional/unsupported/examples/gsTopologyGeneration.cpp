/** @file gsTopologyGeneration.cpp

    @brief Testing the gsTopologyGraph functionality

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/


#include <gismo.h>
#include <gismo_dev.h>

#include <gsUtils/gsTopologyGraph.h>
#include <stdio.h>
#include <time.h>

using std::cout;
using std::endl;
using std::make_pair;
using std::string;
using namespace gismo;
typedef std::pair<unsigned,unsigned> Couple;
typedef std::pair<patchSide,patchSide> PsCouple;

void getGraphBuildingForSides(const std::vector<unsigned> &sides,GraphBuilding &g)
{
    unsigned total=0;
    for(size_t i = 0;i<sides.size();++i)
        total+=sides[i];
    g=GraphBuilding(total);
    unsigned id=0;
    std::vector<size_t> groups;
    for(size_t i = 0;i<sides.size();++i)
    {
        for(unsigned j = 0;j<sides[i];++j)
            groups.push_back(id++);
        if(sides[i]>1)
            g.setTypeIGroup(groups);
        groups.clear();
    }
}

void addCoupledSides(unsigned side0,unsigned side1,bool sameDir,
                     const std::vector<unsigned> &sides,std::vector<Couple> &couples)
{
    GISMO_ASSERT(side0<sides.size()&&side1<sides.size()&&sides[side0]==sides[side1]&&side0!=side1,
                 "cannot couple the two sides");
    if(side0>side1)
    {
        unsigned temp;
        temp=side1;
        side1=side0;
        side0=temp;
    }
    unsigned total=0,start0=0,start1=0;
    for(unsigned i = 0;i<sides.size();++i)
    {
        if(i<side0)
            start0+=sides[i];
        if(i<side1)
            start1+=sides[i];
        total+=sides[i];
    }
    unsigned end1=start1+sides[side1];
    for(unsigned i = 0;i<sides[side0];++i)
    {
        unsigned val0 = start0+i;
        unsigned val1 = sameDir ? start1+i : end1-1-i;
        if(val1==total)
            val1=0;
        if(val0!=val1)
        {
            Couple c = make_pair(val0,val1);
            couples.push_back(c);
        }
    }
}

void addCoupledSide(unsigned side,
                     const std::vector<unsigned> &sides,std::vector<Couple> &couples)
{
    GISMO_ASSERT(side<sides.size(),"cannot couple the side");
    unsigned start=0;
    for(unsigned i = 0;i<side;++i)
        start+=sides[i];
    unsigned end=start+sides[side];
    for(unsigned i = 0;i<sides[side]/2;++i)
    {
        Couple c = make_pair(start+i,end-1-i);
        couples.push_back(c);
    }
}

void get1Side(const std::vector<unsigned> &sides,
                   GraphBuilding &g,std::vector<std::vector<Couple> > &symmetricCouples)
{
    GISMO_ASSERT(sides.size()==1,"wrong amount of sides");
    symmetricCouples.clear();
    getGraphBuildingForSides(sides,g);
}

void get2Sides(const std::vector<unsigned> &sides,
                   GraphBuilding &g,std::vector<std::vector<Couple> > &symmetricCouples)
{
    GISMO_ASSERT(sides.size()==2,"wrong amount of sides");
    symmetricCouples.clear();
    getGraphBuildingForSides(sides,g);
    if(sides[0]==sides[1])
    {
        std::vector<Couple> symm;
        addCoupledSides(0,1,false,sides,symm);
        symmetricCouples.push_back(symm);
    }
}

void get3Sides(const std::vector<unsigned> &sides,
                 GraphBuilding &g,std::vector<std::vector<Couple> > &symmetricCouples)
{
    GISMO_ASSERT(sides.size()==3,"wrong amount of sides");
    symmetricCouples.clear();
    getGraphBuildingForSides(sides,g);
    if(sides[0]==sides[1])
    {
        std::vector<Couple> symm;
        addCoupledSides(0,1,false,sides,symm);
        addCoupledSide(2,sides,symm);
        symmetricCouples.push_back(symm);
    }
    if(sides[1]==sides[2])
    {
        std::vector<Couple> symm;
        addCoupledSides(1,2,false,sides,symm);
        addCoupledSide(0,sides,symm);
        symmetricCouples.push_back(symm);
    }
    if(sides[0]==sides[2])
    {
        std::vector<Couple> symm;
        addCoupledSides(0,2,false,sides,symm);
        addCoupledSide(1,sides,symm);
        symmetricCouples.push_back(symm);
    }
}

void get4Sides(const std::vector<unsigned> &sides,
                 GraphBuilding &g,std::vector<std::vector<Couple> > &symmetricCouples)
{
    GISMO_ASSERT(sides.size()==4,"wrong amount of sides");
    symmetricCouples.clear();
    getGraphBuildingForSides(sides,g);
    if(sides[0]==sides[2])//senkrechte
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(0,2,false,sides,symmetrie);
        addCoupledSide(1,sides,symmetrie);
        addCoupledSide(3,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[1]==sides[3])//waagrechte
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(3,1,false,sides,symmetrie);
        addCoupledSide(0,sides,symmetrie);
        addCoupledSide(2,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[0]==sides[1]&&sides[2]==sides[3])//schräge links unten rechts oben
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(0,1,false,sides,symmetrie);
        addCoupledSides(2,3,false,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[1]==sides[2]&&sides[0]==sides[3])//schräge links oben rechts unten
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(2,1,false,sides,symmetrie);
        addCoupledSides(0,3,false,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
}

void get5Sides(const std::vector<unsigned> &sides,
                   GraphBuilding &g,std::vector<std::vector<Couple> > &symmetricCouples)
{
    GISMO_ASSERT(sides.size()==5,"wrong amount of sides");
    symmetricCouples.clear();
    getGraphBuildingForSides(sides,g);
    if(sides[1]==sides[4]&&sides[2]==sides[3])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(1,4,false,sides,symmetrie);
        addCoupledSides(2,3,false,sides,symmetrie);
        addCoupledSide(0,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[0]==sides[2]&&sides[3]==sides[4])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(0,2,false,sides,symmetrie);
        addCoupledSides(3,4,false,sides,symmetrie);
        addCoupledSide(1,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[1]==sides[3]&&sides[4]==sides[0])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(1,3,false,sides,symmetrie);
        addCoupledSides(4,0,false,sides,symmetrie);
        addCoupledSide(2,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[2]==sides[4]&&sides[0]==sides[1])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(2,4,false,sides,symmetrie);
        addCoupledSides(0,1,false,sides,symmetrie);
        addCoupledSide(3,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[3]==sides[0]&&sides[1]==sides[2])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(3,0,false,sides,symmetrie);
        addCoupledSides(2,1,false,sides,symmetrie);
        addCoupledSide(4,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
}

void get6Sides(const std::vector<unsigned> &sides,
                   GraphBuilding &g,std::vector<std::vector<Couple> > &symmetricCouples)
{
    GISMO_ASSERT(sides.size()==6,"wrong amount of sides");
    symmetricCouples.clear();
    getGraphBuildingForSides(sides,g);
    if(sides[1]==sides[0]&&sides[2]==sides[5]&&sides[3]==sides[4])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(0,1,false,sides,symmetrie);
        addCoupledSides(2,5,false,sides,symmetrie);
        addCoupledSides(3,4,false,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[1]==sides[2]&&sides[3]==sides[0]&&sides[5]==sides[4])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(2,1,false,sides,symmetrie);
        addCoupledSides(3,0,false,sides,symmetrie);
        addCoupledSides(5,4,false,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[3]==sides[2]&&sides[4]==sides[1]&&sides[5]==sides[0])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(2,3,false,sides,symmetrie);
        addCoupledSides(4,1,false,sides,symmetrie);
        addCoupledSides(5,0,false,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[1]==sides[5]&&sides[2]==sides[4])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(5,1,false,sides,symmetrie);
        addCoupledSides(2,4,false,sides,symmetrie);
        addCoupledSide(0,sides,symmetrie);
        addCoupledSide(3,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[2]==sides[0]&&sides[3]==sides[5])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(5,3,false,sides,symmetrie);
        addCoupledSides(2,0,false,sides,symmetrie);
        addCoupledSide(1,sides,symmetrie);
        addCoupledSide(4,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
    if(sides[1]==sides[3]&&sides[0]==sides[4])
    {
        std::vector<Couple> symmetrie;
        addCoupledSides(3,1,false,sides,symmetrie);
        addCoupledSides(0,4,false,sides,symmetrie);
        addCoupledSide(5,sides,symmetrie);
        addCoupledSide(2,sides,symmetrie);
        symmetricCouples.push_back(symmetrie);
    }
}

string getNSides(const std::vector<unsigned> &sides,
                 GraphBuilding &g,std::vector<std::vector<Couple> > &symmetricCouples)
{
    std::stringstream ss;
    switch(sides.size())
    {
    case 0:
        GISMO_ERROR("cannot handle 0 sides.");
        break;
    case 1:
        ss << "OneSide";
        get1Side(sides,g,symmetricCouples);
        break;
    case 2:
        ss << "TwoSides";
        get2Sides(sides,g,symmetricCouples);
        break;
    case 3:
        ss << "ThreeSides";
        get3Sides(sides,g,symmetricCouples);
        break;
    case 4:
        ss << "FourSides";
        get4Sides(sides,g,symmetricCouples);
        break;
    case 5:
        ss << "FiveSides";
        get5Sides(sides,g,symmetricCouples);
        break;
    case 6:
        ss << "SixSides";
        get6Sides(sides,g,symmetricCouples);
        break;
    default:
        ss << sides.size() << "Sides";
        getGraphBuildingForSides(sides,g);
        std::cout << "symmetries not implemented for more than 6 sides." << std::endl;
    }
    ss << ": [" << sides[0];
    for(unsigned i=1;i<sides.size();++i)
        ss << "," << sides[i];
    ss << "]";
    return ss.str();
}

void printTimeStamp()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];
    time (&rawtime);
    timeinfo = localtime (&rawtime);

    strftime (buffer,80,"%F | %R",timeinfo);
    puts (buffer);
}

void createGraphMatrixFiles(int boundaries,int minFaces,int upToFaces,bool keepDuplicates=false)
{
    GraphBuilding g;
    std::vector<std::vector<Couple> > symmetries;
    std::vector<unsigned> sides;
    for(int i = 0;i<boundaries;++i)
        sides.push_back(1);
    string desc=getNSides(sides,g,symmetries);
    GraphBuilding::resetCounters();
    GraphBuilding::setKeepDuplicates(keepDuplicates);
    for(int faces = minFaces;faces<=upToFaces;++faces)
    {
        printTimeStamp();
        cout<<desc<<" with "<<faces<<" faces"<<endl;
        std::vector<Graph> graphs = g.findAllPossibleGraphs(faces);
        if(graphs.size()>0)
        {
            std::stringstream ss;
            if(!keepDuplicates)
                ss<<"Top"<<boundaries<<".xml";
            else
                ss<<"Top"<<boundaries<<"Dup.xml";
            gsInfo << "Writing output to " << ss.str() << "\n";
            Graph::printGraphListToFile(ss.str(),graphs);
        }
    }
}

unsigned countGraphs(int boundaries,int minFaces,int upToFaces, bool countDuplicates=false)
{
    GraphBuilding g;
    std::vector<std::vector<Couple> > symmetries;
    std::vector<unsigned> sides;
    unsigned totalNumber=0;
    for(int i = 0;i<boundaries;++i)
        sides.push_back(1);
    string desc=getNSides(sides,g,symmetries);
    std::vector<unsigned>numbers;
    GraphBuilding::resetCounters();
    GraphBuilding::setKeepDuplicates(countDuplicates);
    for(int faces = minFaces;faces<=upToFaces;++faces)
    {
        printTimeStamp();
        cout<<desc<<" with "<<faces<<" faces"<<endl;
        std::vector<Graph> graphs = g.findAllPossibleGraphs(faces);
        numbers.push_back(graphs.size());
    }
    std::cout<<"Summary:\n";
    if(countDuplicates)
        std::cout<<"DUPLICATES counted also!"<<std::endl;
    unsigned index=0;
    for(int faces=minFaces;faces<=upToFaces;++faces)
    {
        std::cout<<"Faces: " << faces << " Count: " << numbers[index]<<"\n";
        totalNumber+=numbers[index];
        index++;
    }
    std::cout<<std::endl;
    return totalNumber;
}

bool readGraphs(int boundarySize,int faces,std::vector<Graph>& graphList,bool getDuplicates=false)
{
    bool success=true;
    std::stringstream ss;
    if(!getDuplicates)
        ss<<"TopologyMatrizes/"<<"Top"<<boundarySize<<".xml";
    else
        ss<<"TopologyMatrizes/"<<"Top"<<boundarySize<<"Dup.xml";
    std::string file = gsFileManager::find(ss.str());

    std::vector<Graph> graphs;
    Graph::readGraphListFromFile(file,graphs);
    std::vector<std::vector<Graph> > graphPartitions;
    Graph::getPartitionsOfGraphVector(graphs,graphPartitions,faces+boundarySize,faces+boundarySize);

    graphList.clear();
    if(graphPartitions.size()==1)
        graphList=graphPartitions[0];
    else
        success=false;

    return success;
}

void sortOutNonUniqueGraphs(std::vector<Graph>& graphList)
{
    Graph g;
    bool found;
    std::vector<unsigned> uniqueIndexList,duplicateIndexList;
    for(unsigned i = 0;i<graphList.size();++i)
    {
        g=graphList[i];
        found=false;
        for(index_t j = 0;j<g.getNrBoundaryVertices();++j)
        {
            for(unsigned k = 0;k<uniqueIndexList.size();++k)
            {
                if(Graph::equal(graphList[uniqueIndexList[k]],g))
                    found=true;
                if(found)
                    break;
            }
            if(found)
                break;
            g=g.rotateBoundaryByOne();
            g=g.normalizeGraph();
        }
        if(!found)
            uniqueIndexList.push_back(i);
        else
            duplicateIndexList.push_back(i);
    }
    for(int i = duplicateIndexList.size()-1;i>=0;i--)
    {
        graphList.erase(graphList.begin()+duplicateIndexList[i]);
    }
}

bool isNotTowEdgeConnected(Graph g) {return !(g.innerGraphIsTwoEdgeConnected());}

void sortOutNonTwoEdgeConnectedGraphs(std::vector<Graph>& graphList)
{
    graphList.erase(std::remove_if(graphList.begin(),graphList.end(),isNotTowEdgeConnected),graphList.end());
}

bool hasValenceHigherThanFive(Graph g)
{
    gsBoxTopology topol = ConstructMultiPatchFromGraph::getBoxTopology(g);
    return topol.getMaxValence()>5;
}

void sortOutValenceSmallerEqualFive(std::vector<Graph>& graphList)
{
    graphList.erase(std::remove_if(graphList.begin(),graphList.end(),hasValenceHigherThanFive),graphList.end());
}

int main()
{
    int success=0;
    int action=0; // select a run action
    bool workWithDuplicates=false; // select true if duplicate graphs should be considered
    switch(action)
    {
    case 0: // standard test
    {
        const unsigned c = countGraphs(2,1,6,workWithDuplicates);
        if(c==5)
        {
            gsInfo<<"countGraphs(2,1,6) = "<< c <<", OK.\n.";
        }
        else
        {
            gsInfo<<"countGraphs(2,1,6) was "<< c <<" instead of "<<5 <<".\n.";
            success = 1;
        }
        break;
    }
    case 1: // create graph matrix files
    {
//        createGraphMatrixFiles(2,1,8,workWithDuplicates);
//        createGraphMatrixFiles(4,1,8,workWithDuplicates);
//        createGraphMatrixFiles(6,1,8,workWithDuplicates);
//        createGraphMatrixFiles(8,1,7,workWithDuplicates);
//        createGraphMatrixFiles(10,1,7,workWithDuplicates);
//        createGraphMatrixFiles(12,1,7,workWithDuplicates);
        createGraphMatrixFiles(14,1,9,workWithDuplicates);
        createGraphMatrixFiles(16,1,9,workWithDuplicates);
        createGraphMatrixFiles(18,1,9,workWithDuplicates);
        createGraphMatrixFiles(20,1,9,workWithDuplicates);
        createGraphMatrixFiles(22,1,9,workWithDuplicates);
        createGraphMatrixFiles(24,1,9,workWithDuplicates);
        break;
    }
    case 2: // count graphs function
    {
        countGraphs(2,8,9,workWithDuplicates);
        countGraphs(4,8,9,workWithDuplicates);
        countGraphs(6,8,9,workWithDuplicates);
        countGraphs(8,1,7,workWithDuplicates);
        countGraphs(10,1,7,workWithDuplicates);
        countGraphs(12,1,7,workWithDuplicates);
        break;
    }
    case 3: // count uniqueGraphs only
    {
        bool valenceFiveOrLower=false; // if true only valence five or lower is considered
        gsInfo<<"unique graphs";
        if(valenceFiveOrLower)
            gsInfo<<" with valence 5 or lower";
        gsInfo<<": "<<std::endl;
        bool successfulRead;
        std::vector<Graph> graphList;
        for(unsigned boundarySides=2;boundarySides<=12;boundarySides=boundarySides+2)
        {
            gsInfo<<"Boundary Sides: "<<boundarySides<<std::endl<<std::endl;
            for(unsigned faces = 1;faces<=9;++faces)
            {
                graphList.clear();
                successfulRead=readGraphs(boundarySides,faces,graphList,workWithDuplicates);
                if(!successfulRead)
                {
                    gsInfo<<"Faces: " << faces << " Count: 0"<<std::endl;
                    continue;
                }
                sortOutNonUniqueGraphs(graphList);
                if(valenceFiveOrLower)
                    sortOutValenceSmallerEqualFive(graphList);
                gsInfo<<"Faces: " << faces << " Count: "<<graphList.size()<<std::endl;
            }
            gsInfo<<std::endl;
        }
        break;
    }
    case 4: // count unique graphs which are twoEdgeConnected
    {
        bool valenceFiveOrLower=false; // if true only valence five or lower is considered
        gsInfo<<"unique, two-edge-connected graphs";
        if(valenceFiveOrLower)
            gsInfo<<" with valence 5 or lower";
        gsInfo<<": "<<std::endl;
        bool successfulRead;
        std::vector<Graph> graphList;
        for(unsigned boundarySides=2;boundarySides<=12;boundarySides=boundarySides+2)
        {
            gsInfo<<"Boundary Sides: "<<boundarySides<<std::endl<<std::endl;
            for(unsigned faces = 1;faces<=9;++faces)
            {
                graphList.clear();
                successfulRead=readGraphs(boundarySides,faces,graphList,workWithDuplicates);
                if(!successfulRead)
                {
                    gsInfo<<"Faces: " << faces << " Count: 0"<<std::endl;
                    continue;
                }
                sortOutNonUniqueGraphs(graphList);
                sortOutNonTwoEdgeConnectedGraphs(graphList);
                if(valenceFiveOrLower)
                    sortOutValenceSmallerEqualFive(graphList);
                gsInfo<<"Faces: " << faces << " Count: "<<graphList.size()<<std::endl;
            }
            gsInfo<<std::endl;
        }
        break;
    }
    default:
    {
        gsInfo<<"No Correct option was selected."<<std::endl;
        success = 1;
        break;
    }
    }
    return success;
}
