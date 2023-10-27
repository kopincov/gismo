/** @file gsTopologyGraph.cpp

    @brief Implementation File

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger, J. Vogl
*/

#include <gsCore/gsBoundary.h>
#include <gsIO/gsFileData.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsCore/gsMultiPatch.h>
#include<gsUtils/gsTopologyGraph.h>
#include <gsCore/gsBoxTopology.h>

#include <fstream>

namespace gismo {

size_t GraphBuilding::m_countRecursives=0;
size_t GraphBuilding::m_countTrees=0;
bool GraphBuilding::m_keepDuplicateGraphs=false;

void Graph::addEdge(size_t v1, size_t v2, short_t pos1, short_t pos2)
{
    patchSide ps1(v1, posToSide(pos1));
    patchSide ps2(v2, posToSide(pos2));
    PsCouple c = std::make_pair(ps1, ps2);
    GraphBuilding::insertNonDuplicatesOrdered(m_edges, c);
}
void Graph::addEdge(patchSide ps1,patchSide ps2)
{
    PsCouple c=std::make_pair(ps1,ps2);
    GraphBuilding::insertNonDuplicatesOrdered(m_edges,c);
}

bool Graph::removeEdge(PsCouple ps1)
{
    for (size_t i = 0; i < m_edges.size(); ++i)
        if (m_edges[i] == ps1)
        {
            m_edges.erase(m_edges.begin() + i);
            return true;
        }
    return false;
}

std::vector<size_t> Graph::vertexConnections(size_t vertex) const
{
    std::vector<size_t> result;
    for(size_t i = 0;i<m_edges.size();++i)
    {
        if(m_edges[i].first.patch == static_cast<index_t>(vertex))
            result.push_back(m_edges[i].second.patch);
        else if(m_edges[i].second.patch== static_cast<index_t>(vertex))
            result.push_back(m_edges[i].first.patch);
    }
    return result;
}

bool Graph::equal(const Graph & g1, Graph const & g2)
{
    if(g1.m_vertices!=g2.m_vertices)
        return false;
    bool equal = true;
    for (size_t i = 0; i < g1.m_edges.size(); ++i)
        if(g1.m_edges[i].first.patch!=g2.m_edges[i].first.patch ||
                g1.m_edges[i].second.patch!=g2.m_edges[i].second.patch)
            equal=false;
    return equal;
}

bool Graph::similar(Graph const & g1, Graph const & g2)
{
    if(g1.m_vertices!=g2.m_vertices)
        return false;
    Graph newG1 = g1.normalizeGraph();
    Graph newG2 = g2.normalizeGraph();
    return equal(newG1,newG2);
}

Graph Graph::normalizeGraph() const
{
    Graph oldG = *this;
    Graph newG(m_vertices, m_boundaryVertices);
    newG.putEdgesInGraph(oldG.m_edges, m_boundaryVertices - 1);
    for (index_t i = m_boundaryVertices; i < newG.m_vertices; ++i)
    {
        size_t smallest = oldG.findSmallestVertex(i);
        swapInVertices(oldG.m_edges, smallest, i);
        newG.putEdgesInGraph(oldG.m_edges, i);
    }
    if (oldG.m_edges.size() > 0)
        GISMO_ERROR("not all edges transfered to the new graph.");
    return newG;
}

bool Graph::innerGraphIsTwoEdgeConnected() const
{
    std::vector<PsCouple> edges = this->getEdges();
    Graph newG;
    for (size_t i = 0; i < edges.size(); ++i)
    {
        if (edges[i].first.patch < m_boundaryVertices ||
            edges[i].second.patch < m_boundaryVertices)
            continue;
        newG = *this;
        newG.removeEdge(edges[i]);
        if (!newG.innerGraphIsConnected())
            return false;
    }
    return true;
}

void Graph::travelGraph(size_t vertex,std::vector<bool>& visited) const
{
    visited[vertex]=true;
    std::vector<size_t> connects = this->vertexConnections(vertex);
    if(connects.size()==0)
        return;
    for(size_t i = 0;i<connects.size();++i)
        if(!visited[connects[i]])
            travelGraph(connects[i],visited);
}


bool Graph::innerGraphIsConnected() const
{
    std::vector<bool> visited;
    for(index_t i = 0;i<this->getNrVertices();++i)
        visited.push_back(i<m_boundaryVertices);
    size_t startVertex=m_boundaryVertices;
    travelGraph(startVertex,visited);
    for(size_t i = 0;i<visited.size();++i)
        if(!visited[i])
            return false;
    return true;
}

Graph Graph::rotateBoundaryByOne() const
{
    Graph newG=*this;
    for(index_t i = 1;i<m_boundaryVertices;++i)
    {
        swapInVertices(newG.m_edges,0,i);
    }
    return newG;
}

bool Graph::isSymmetric(std::vector<Couple> symmetricVertices) const
{
    Graph gNew=*this;
    for(size_t i = 0;i<symmetricVertices.size();++i)
        swapInVertices(gNew.m_edges,symmetricVertices[i].first,symmetricVertices[i].second);
    return similar(*this,gNew);
}

void Graph::levelMatrixToGraphs(const gsMatrix<index_t>& mat,std::vector<Graph>& graphs)
{
    graphs.clear();
    GISMO_ASSERT(mat.rows()>0&&mat.cols()>0,"Matrix should have positive size.");
    int vertices = mat(0,0);
    int boundVertices = mat(0,1);
    for(index_t i = 0;i<mat.rows();++i)
    {
        GISMO_ASSERT(mat(i,0)==vertices,"All graphs should have same number of vertices.");
        GISMO_ASSERT(mat(i,1)==boundVertices,"All graphs should have same number of boundary vertices.");
        Graph g(vertices,boundVertices);
        for(index_t j = 2;j<mat.cols();j=j+4)
        {
            if(mat(i,j)!=-1)
            {
                patchSide ps1(mat(i,j),mat(i,j+1));
                patchSide ps2(mat(i,j+2),mat(i,j+3));
                g.addEdge(ps1,ps2);
            }
        }
        graphs.push_back(g);
    }
}

void Graph::graphsToLevelMatrix(std::vector<Graph>& graphs,gsMatrix<index_t>& mat)
{
    GISMO_ASSERT(graphs.size()>0,"Graph vector should have positive size.");
    int vertices=-1,max=-1,edges_i,vertices_i,boundary_i;
    std::vector<PsCouple> edges;
    for(size_t i = 0;i<graphs.size();++i)
    {
        edges_i = graphs[i].getEdges().size();
        max=edges_i*4+2>max ? edges_i*4+2 : max;
        vertices_i = graphs[i].getNrVertices();
        vertices=vertices_i>vertices ? vertices_i : vertices;
    }
    mat.resize(graphs.size(),max);
    for(size_t i = 0;i<graphs.size();++i)
    {
        vertices_i = graphs[i].getNrVertices();
        boundary_i = graphs[i].getNrBoundaryVertices();
        GISMO_ASSERT(vertices==vertices_i,"All Graphs should have the same number of vertices.");
        edges = graphs[i].getEdges();
        mat(i,0)=vertices_i;
        mat(i,1)=boundary_i;
        size_t j=0;
        for(;j<edges.size();++j)
        {
            mat(i,2+j*4)   = edges[j].first.patch;
            mat(i,2+j*4+1) = edges[j].first.side();
            mat(i,2+j*4+2) = edges[j].second.patch;
            mat(i,2+j*4+3) = edges[j].second.side();
        }
        for(index_t k=j*4+2;k<mat.cols();++k)
            mat(i,k)=-1;
    }
}

void Graph::getPartitionsOfGraphVector(std::vector<Graph>& graphs,std::vector<std::vector<Graph> >& graphPartitions,
                                index_t minFaces, index_t maxFaces)
{
    graphPartitions.clear();
    std::vector<size_t> limits;
    index_t vertices=-1,vertices_i;
    size_t i=0;
    for(;i<graphs.size();++i)
    {
        vertices_i=graphs[i].getNrVertices();
        if(vertices<vertices_i)
        {
            if(vertices_i>maxFaces)
                break;
            if(minFaces<=vertices_i)
                limits.push_back(i);
            vertices=vertices_i;
        }
    }
    limits.push_back(i);
    std::vector<Graph> temp;
    for(i = 1;i<limits.size();++i)
    {
        temp.clear();
        for(size_t j = limits[i-1];j<limits[i];++j)
        {
            temp.push_back(graphs[j]);
        }
        graphPartitions.push_back(temp);
    }
}

void Graph::readGraphListFromFile(std::string filename,std::vector<Graph>& graphs)
{
    graphs.clear();
    gsFileData<real_t> graphMatrizes( filename );
    int numberOfLevels = graphMatrizes.count<gsMatrix<index_t> >();
    gsMatrix<index_t>::uPtr levelMatrix;
    std::vector<Graph> graphsOfLevel;
    for(int i = 0;i<numberOfLevels;++i)
    {
        levelMatrix = graphMatrizes.getId<gsMatrix<index_t> >(i);
        levelMatrixToGraphs(*(levelMatrix.release()),graphsOfLevel);
        graphs.insert( graphs.end(), graphsOfLevel.begin(), graphsOfLevel.end() );
    }
}

void Graph::printGraphListToFile(std::string filename,std::vector<Graph>& graphs)
{
    std::vector<std::vector<Graph> > oldGraphPartitions,newGraphPartitions,graphPartitions;
    if( std::ifstream(filename.c_str()).good() )
    {
        std::vector<Graph> oldGraphs;
        readGraphListFromFile(filename,oldGraphs);
        getPartitionsOfGraphVector(oldGraphs,oldGraphPartitions);
    }

    getPartitionsOfGraphVector(graphs,newGraphPartitions);

    size_t i = 0,j = 0;
    std::vector<Graph> temp;
    while(i<oldGraphPartitions.size()||j<newGraphPartitions.size())
    {
        if(j<newGraphPartitions.size() &&
                (i>=oldGraphPartitions.size()||oldGraphPartitions[i][0].getNrVertices()>=newGraphPartitions[j][0].getNrVertices()))
        {
            temp.clear();
            for(size_t k = 0;k<newGraphPartitions[j].size();++k)
                temp.push_back(newGraphPartitions[j][k]);
            graphPartitions.push_back(temp);
            if(i<oldGraphPartitions.size()&&oldGraphPartitions[i][0].getNrVertices()==newGraphPartitions[j][0].getNrVertices())
                i++;
            j++;
        }
        else if(j>=newGraphPartitions.size()||oldGraphPartitions[i][0].getNrVertices()<newGraphPartitions[j][0].getNrVertices())
        {
            temp.clear();
            for(size_t k = 0;k<oldGraphPartitions[i].size();++k)
                temp.push_back(oldGraphPartitions[i][k]);
            graphPartitions.push_back(temp);
            i++;
        }
    }

    gsFileData<real_t> output;
    gsMatrix<index_t> mat;
    for(size_t k = 0;k<graphPartitions.size();++k)
    {
        graphsToLevelMatrix(graphPartitions[k],mat);
        output << mat;
    }
    output.dump(filename);
}

short_t Graph::posToSide(short_t pos)
{
    switch(pos)
    {
        case 1: return 3;
        case 2: return 1;
        case 3: return 4;
        case 4: return 2;
        default: return 0;
    }
}

size_t Graph::findSmallestVertex(index_t start) const
{
    std::vector<size_t> lowestConnectedVertex(m_vertices-start,m_vertices+1);
    std::vector<size_t> highestNumberOfConnectedVertices(m_vertices-start,0);
    std::vector<size_t> lowestSumConnectedVertices(m_vertices-start,0);
    for(size_t i=0;i<m_edges.size();++i)
    {
        if(m_edges[i].first.patch>=start)
        {
            if(lowestConnectedVertex[m_edges[i].first.patch-start]>(size_t)m_edges[i].second.patch)
                lowestConnectedVertex[m_edges[i].first.patch-start]=m_edges[i].second.patch;
            if(m_edges[i].second.patch<start)
            {
                highestNumberOfConnectedVertices[m_edges[i].first.patch-start]+=1;
                lowestSumConnectedVertices[m_edges[i].first.patch-start]+=m_edges[i].second.patch;
            }
        }
        if(m_edges[i].second.patch>=start)
        {
            //cout << "reach here";
            if(lowestConnectedVertex[m_edges[i].second.patch-start]>(size_t)m_edges[i].first.patch)
                lowestConnectedVertex[m_edges[i].second.patch-start]=m_edges[i].first.patch;
            if(m_edges[i].first.patch<start)
            {
                highestNumberOfConnectedVertices[m_edges[i].second.patch-start]+=1;
                lowestSumConnectedVertices[m_edges[i].second.patch-start]+=m_edges[i].first.patch;
            }
        }
    }
    std::vector<size_t> candidates;
    size_t min = *std::min_element(lowestConnectedVertex.begin(),lowestConnectedVertex.end());
    for(size_t i = 0;i<lowestConnectedVertex.size();++i)
        if(lowestConnectedVertex[i]==min)
            candidates.push_back(i);
    if(candidates.size()==1)
        return candidates[0]+start;

    std::vector<size_t> highestNumberCandidates;
    for(size_t i = 0;i<candidates.size();++i)
        highestNumberCandidates.push_back(highestNumberOfConnectedVertices[candidates[i]]);
    size_t max = *std::max_element(highestNumberCandidates.begin(),highestNumberCandidates.end());
    std::vector<size_t> maxCandidates;
    for(size_t i = 0;i<highestNumberCandidates.size();++i)
        if(highestNumberCandidates[i]==max)
            maxCandidates.push_back(candidates[i]);
    if(maxCandidates.size()==1)
        return maxCandidates[0]+start;

    std::vector<size_t> lowestSumCandidates;
    for(size_t i = 0;i<maxCandidates.size();++i)
        lowestSumCandidates.push_back(lowestSumConnectedVertices[maxCandidates[i]]);
    size_t minSum = *std::min_element(lowestSumCandidates.begin(),lowestSumCandidates.end());
    std::vector<size_t> minSumCandidates;
    for(size_t i = 0;i<lowestSumCandidates.size();++i)
        if(lowestSumCandidates[i]==minSum)
            minSumCandidates.push_back(maxCandidates[i]);
    if(minSumCandidates.size()==1)
        return minSumCandidates[0]+start;
    else
        return start;
        //GISMO_ERROR("could not find the smallest vertex.");
}

void Graph::swapInVertices(std::vector<PsCouple> &edges,index_t swap1,index_t swap2)
{
    if(swap1==swap2)
        return;
    for(size_t i = 0;i<edges.size();++i)
    {
        if(edges[i].first.patch==swap1)
            edges[i].first.patch=swap2;
        else if(edges[i].first.patch==swap2)
            edges[i].first.patch=swap1;
        if(edges[i].second.patch==swap1)
            edges[i].second.patch=swap2;
        else if(edges[i].second.patch==swap2)
            edges[i].second.patch=swap1;
    }
}

void Graph::putEdgesInGraph(std::vector<PsCouple> &edges, index_t limit)
{
    std::deque<size_t> removeList;
    for(size_t i=0;i<edges.size();++i)
    {
        if(edges[i].first.patch<=limit&&
                edges[i].second.patch<=limit)
        {
            addEdge(edges[i].first,edges[i].second);
            removeList.push_front(i);
        }
    }
    for (size_t i = 0; i < removeList.size(); ++i)
        edges.erase(edges.begin()+removeList[i]);
}

GraphBuilding::GraphBuilding(size_t nrEdges,std::vector<Couple> &typeI) :
    m_typeI(typeI),m_origLimit(nrEdges),m_maxElement(nrEdges-1),m_graph(nrEdges,nrEdges)
{
    for(size_t i =0;i<nrEdges;++i)
    {
        m_edges.push_back(i);
        m_origEdges=m_edges;
        std::vector<size_t> rel;
        rel.push_back(i);
        m_edgeToVertexRelations.push_back(rel);
    }
}

GraphBuilding::GraphBuilding(size_t nrEdges) :
    m_origLimit(nrEdges),m_maxElement(nrEdges-1),m_graph(nrEdges,nrEdges)
{
    for(size_t i =0;i<nrEdges;++i)
    {
        m_edges.push_back(i);
        m_origEdges=m_edges;
        std::vector<size_t> rel;
        rel.push_back(i);
        m_edgeToVertexRelations.push_back(rel);
    }
}

void GraphBuilding::setTypeIGroup(std::vector<size_t> &edges)
{
    for(size_t i = 0;i<edges.size();++i)
        for(size_t j = i+1;j<edges.size();++j)
            m_typeI.push_back(std::make_pair(edges[i],edges[j]));
}

void GraphBuilding::setTypeIIGroup(std::vector<size_t> &edges)
{
    for(size_t i = 0;i<edges.size();++i)
        for(size_t j = i+1;j<edges.size();++j)
            m_typeII.push_back(std::make_pair(edges[i],edges[j]));
}

std::vector<GraphBuilding::Couple> GraphBuilding::possibleUnifiedEdges() const
{
    std::vector<Couple> possibleEdges;
//        for(size_t i = 0;i<m_edges.size();++i)
//            for(size_t j = i+1;j<m_edges.size();++j)
//                possibleEdges.push_back(std::make_pair(m_edges[i],m_edges[j]));
//        for(size_t i = 0;i<m_edges.size()-1;++i)
//            possibleEdges.push_back(std::make_pair(m_edges[i],m_edges[i+1]));
//        possibleEdges.push_back(std::make_pair(m_edges[0],m_edges.back()));
    for (size_t i = 1; i < m_edges.size(); i += 2)
        possibleEdges.push_back(std::make_pair(m_edges[0], m_edges[i]));
    for (size_t i = possibleEdges.size() - 1; i < -1UL; --i)
        if (!possibleConnect(possibleEdges[i]))
            possibleEdges.erase(possibleEdges.begin() + i);
    return possibleEdges;
}

bool GraphBuilding::possibleConnect(Couple c) const
{
    // original edges cannot be connected
    if( c.first<m_origLimit && c.second<m_origLimit )
        return false;
    for(size_t i = 0;i<m_alreadyTried.size();++i)
        if(m_alreadyTried[i]==c || (m_alreadyTried[i].first==c.second && m_alreadyTried[i].second==c.first))
            return false;
    // check typeI and typeII connectivity
    Couple c2 = std::make_pair(c.second,c.first);
    for(size_t i = 0;i<m_typeI.size();++i)
        if(m_typeI[i]==c || m_typeI[i]==c2)
            return false;
    for(size_t i = 0;i<m_typeII.size();++i)
        if(m_typeII[i]==c || m_typeII[i]==c2)
            return false;

    bool found=false;
    size_t lower=0,upper=0;
    for(size_t i = 0;i<m_edges.size();++i)
    {
        if(c.first==m_edges[i] || c.second==m_edges[i])
        {
            if(found)
                upper=i;
            else
                lower=i;
            found=!found;
        }
    }
    GISMO_ASSERT(upper!=0 || lower!=0,"upper and lower not found");
    // even number of edges in each subgraph
    if((upper-lower)%2 == 0)
        return false;
    // connectivity criteria
    for(size_t i = lower+1;i<upper;++i)
        if(m_origEdges[lower]>=m_origEdges[i] || m_origEdges[i]>=m_origEdges[upper])
            return true;
    // if they are neighbours then its ok too
    return upper-lower<=1;
}

bool GraphBuilding::isActiveEdge(size_t e) const
{
    bool found = false;
    for(size_t i = 0;i<m_edges.size();++i)
        if(e==m_edges[i])
            found = true;
    return found;
}

std::pair<GraphBuilding,GraphBuilding> GraphBuilding::unifyEdge(Couple c) const
{
    //std::cout << "unify: (" << c.first << "," <<  c.second << ")" << endl;
    size_t firstIndex = find(m_edges.begin(),m_edges.end(),c.first) - m_edges.begin();
    size_t secondIndex = find(m_edges.begin(),m_edges.end(),c.second) - m_edges.begin();
    if(!(firstIndex<m_edges.size()&&secondIndex<m_edges.size()))
    {
        std::cout << "edges: ";
        for(size_t i = 0;i<m_edges.size();++i)
            std::cout << m_edges[i] << " ";
        std::cout << "try to unify: (" << c.first << "," << c.second << ")" <<std::endl;
        GISMO_ERROR("elements not found in edge-list");
    }
    size_t low = ( firstIndex>secondIndex ? secondIndex : firstIndex );
    size_t high = ( firstIndex>secondIndex ? firstIndex : secondIndex );

    GraphBuilding g1,g2;
    for(size_t i = 0;i<m_edges.size();++i)
    {
        if( i<low || i>high )
        {
            g1.m_edges.push_back(m_edges[i]);
            g1.m_origEdges.push_back(m_origEdges[i]);
        }
        else if(i>low && i<high)
        {
            g2.m_edges.push_back(m_edges[i]);
            g2.m_origEdges.push_back(m_origEdges[i]);
        }
    }
    std::vector<Couple> newTypeIIs;
    for(size_t i = 0;i<m_typeI.size();++i)
    {
        Couple typeI = m_typeI[i];
        size_t tFirstIndex = find(m_edges.begin(),m_edges.end(),typeI.first) - m_edges.begin();
        size_t tSecondIndex = find(m_edges.begin(),m_edges.end(),typeI.second) - m_edges.begin();
        if(!(tFirstIndex<m_edges.size()&&tSecondIndex<m_edges.size()))
        {
            std::cout << "edges: ";
            for(size_t k = 0;k<m_edges.size();++k)
                std::cout << m_edges[k] << " ";
            std::cout << "looking for typeI: (" << typeI.first << "," << typeI.second << ")" <<std::endl;
            GISMO_ERROR("elements not found in edge-list");
        }
        size_t lowTypeI = ( tFirstIndex>tSecondIndex ? tSecondIndex : tFirstIndex );
        size_t highTypeI = ( tFirstIndex>tSecondIndex ? tFirstIndex : tSecondIndex );
        if(( lowTypeI<low || lowTypeI>high ) && (highTypeI<low || highTypeI>high))
            g1.m_typeI.push_back(typeI);
        else if(lowTypeI>low && lowTypeI<high && highTypeI>low && highTypeI<high)
            g2.m_typeI.push_back(typeI);


        if(typeI.first==c.first||typeI.second==c.first||typeI.first==c.second||typeI.second==c.second)
        {
            size_t first = (typeI.first==c.first||typeI.first==c.second) ? typeI.second : typeI.first;
            size_t cPart = (typeI.first==c.first||typeI.second==c.first) ? c.second : c.first;
            for(size_t k = 0;k<m_typeI.size();++k)
            {
                if( m_typeI[k].first==cPart || m_typeI[k].second==cPart )
                {
                    size_t second = m_typeI[k].first==cPart ? m_typeI[k].second : m_typeI[k].first;
                    Couple newTypeII=std::make_pair(first,second);
                    newTypeIIs.push_back(newTypeII);
                }
            }
        }
    }
    for(size_t i = 0;i<m_typeII.size()+newTypeIIs.size();++i)
    {
        Couple typeII = i<m_typeII.size() ? m_typeII[i] : newTypeIIs[i-m_typeII.size()];
        size_t tFirstIndex = find(m_edges.begin(),m_edges.end(),typeII.first) - m_edges.begin();
        size_t tSecondIndex = find(m_edges.begin(),m_edges.end(),typeII.second) - m_edges.begin();
        if(!(tFirstIndex<m_edges.size()&&tSecondIndex<m_edges.size()))
        {
            std::cout << "edges: ";
            for(size_t k = 0;k<m_edges.size();++k)
                std::cout << m_edges[k] << " ";
            std::cout << "looking for typeII: (" << typeII.first << "," << typeII.second << ")" <<std::endl;
            GISMO_ERROR("elements not found in edge-list");
        }
        size_t lowTypeII = ( tFirstIndex>tSecondIndex ? tSecondIndex : tFirstIndex );
        size_t highTypeII = ( tFirstIndex>tSecondIndex ? tFirstIndex : tSecondIndex );
        if(( lowTypeII<low || lowTypeII>high ) && (highTypeII<low || highTypeII>high))
            g1.m_typeII.push_back(typeII);
        else if(lowTypeII>low && lowTypeII<high && highTypeII>low && highTypeII<high)
            g2.m_typeII.push_back(typeII);
    }
    g1.m_origLimit=m_origLimit;
    g2.m_origLimit=m_origLimit;
    g1.m_graph=m_graph;
    g2.m_graph=m_graph;
    g1.m_edgeToVertexRelations=m_edgeToVertexRelations;
    g2.m_edgeToVertexRelations=m_edgeToVertexRelations;
    g1.m_maxElement=m_maxElement;
    g2.m_maxElement=m_maxElement;
    g1.m_alreadyTried=m_alreadyTried;
    g2.m_alreadyTried=m_alreadyTried;
    size_t connectVert1 = 0,connectVert2 = 0, pos1, pos2;
#ifndef NDEBUG
    GISMO_ASSERT(findVertexInRelations(c.first, connectVert1, pos1),"connect vertex not found in relations list");
    GISMO_ASSERT(findVertexInRelations(c.second,connectVert2, pos2),"connect vertex not found in relations list");
#else
    findVertexInRelations(c.first, connectVert1, pos1);
    findVertexInRelations(c.second,connectVert2, pos2);
#endif
    g1.m_graph.addEdge(connectVert1,connectVert2,pos1+1,pos2+1);
    g2.m_graph.addEdge(connectVert1,connectVert2,pos1+1,pos2+1);
    return std::make_pair(g1,g2);
}

Graph GraphBuilding::combineGraphs(Graph g1,Graph g2)
{
    GISMO_ASSERT(g1.getNrVertices()==g2.getNrVertices(),"two graphs have to agree on the vertices");
    GISMO_ASSERT(g1.getNrBoundaryVertices()==g2.getNrBoundaryVertices(),"two graphs have to agree on the boundary vertices");
    Graph combined(g1.getNrVertices(),g1.getNrBoundaryVertices());
    std::vector<PsCouple> g1Edges=g1.getEdges();
    std::vector<PsCouple> g2Edges=g2.getEdges();
    for(size_t i = 0;i<g1Edges.size();++i)
        combined.addEdge(g1Edges[i].first,g1Edges[i].second);
    for(size_t i = 0;i<g2Edges.size();++i)
        combined.addEdge(g2Edges[i].first,g2Edges[i].second);
    return combined;
}

std::vector<size_t> GraphBuilding::possibleVertexLocations(size_t lastInserted) const
{
    std::vector<size_t> pos;
    if(lastInserted == -1UL)
        pos.push_back(m_edges.back());
    else
    {
        for(size_t i = 0;i<m_edges.size();++i)
            if(m_edges[i]>static_cast<size_t>(lastInserted))
                pos.push_back(m_edges[i]);
    }
    return pos;
}

GraphBuilding GraphBuilding::getGraphBuildingWithAddedVertex(size_t connect) const
{
    //std::cout << "add Vertex: " << connect << std::endl;
    bool found = false;
    size_t i = 0;
    for(;i<m_edges.size();++i)
        if(m_edges[i]==connect)
        {
            found=true;
            break;
        }
    if( !found )
    {
        std::cout << "edges: ";
        for(i = 0;i<m_edges.size();++i)
            std::cout << m_edges[i] << " ";
        std::cout << "try to add: " << connect <<std::endl;
        GISMO_ERROR("not possible to add the vertex");
    }

    GraphBuilding g = *this;
    std::vector<size_t> newEdges;
    newEdges.push_back(m_maxElement+1);
    newEdges.push_back(m_maxElement+2);
    newEdges.push_back(m_maxElement+3);
    g.m_edges.insert(g.m_edges.begin()+i+1,newEdges.begin(),newEdges.end());
    g.m_edges.erase(g.m_edges.begin()+i);
    g.m_origEdges.insert(g.m_origEdges.begin()+i+1,2,m_origEdges[i]);
    g.removeTypeII(connect);
    g.convertTypeItoTypeII(connect,newEdges);
    std::vector<size_t> typeIEdges;
    typeIEdges.push_back(newEdges[0]);
    typeIEdges.push_back(newEdges[1]);
    g.setTypeIGroup(typeIEdges);
    typeIEdges.clear();
    typeIEdges.push_back(newEdges[1]);
    typeIEdges.push_back(newEdges[2]);
    g.setTypeIGroup(typeIEdges);
    std::vector<size_t> typeIIEdges;
    typeIIEdges.push_back(newEdges[0]);
    typeIIEdges.push_back(newEdges[1]);
    typeIIEdges.push_back(newEdges[2]);
    g.setTypeIIGroup(typeIIEdges);

    g.m_graph.addVertices(1);
    g.m_maxElement+=3;
    size_t connectVert=0;
    size_t side;
#ifndef NDEBUG
    GISMO_ASSERT(findVertexInRelations(connect,connectVert, side),"connect vertex not found in relations list");
#else
    findVertexInRelations(connect,connectVert, side);
#endif
    g.m_graph.addEdge(connectVert,g.m_graph.getNrVertices()-1,side+1,4);
    newEdges.push_back(connect);
    g.m_edgeToVertexRelations.push_back(newEdges);
    return g;
}

bool GraphBuilding::findVertexInRelations(size_t vert,size_t &rel, size_t &side) const
{
    for (size_t i = 0; i < m_edgeToVertexRelations.size(); ++i)
        for (size_t j = 0; j < m_edgeToVertexRelations[i].size(); ++j)
            if (m_edgeToVertexRelations[i][j] == vert)
            {
                rel = i;
                side = j;
                return true;
            }
    return false;
}

void GraphBuilding::insertNonDuplicatesOrdered(std::vector<Couple> &couples,Couple insert)
{
    size_t big = insert.first>insert.second ? insert.first : insert.second;
    size_t small = insert.first>insert.second ? insert.second : insert.first;
    Couple c = std::make_pair(small,big);
    for(size_t pos=0;pos<couples.size();++pos)
    {
        if(couples[pos].first==small&&couples[pos].second==big)
            return;
    }
    couples.push_back(c);
    sortEdges(couples);
}

void GraphBuilding::insertNonDuplicatesOrdered(std::vector<PsCouple> &couples,PsCouple insert)
{
    patchSide big = insert.first.patch>insert.second.patch ? insert.first : insert.second;
    patchSide small = insert.first.patch>insert.second.patch ? insert.second : insert.first;
    PsCouple c = std::make_pair(small,big);
    for(size_t pos=0;pos<couples.size();++pos)
    {
        if(couples[pos].first==small&&couples[pos].second==big)
            return;
    }
    couples.push_back(c);
    sortEdges(couples);
}

void GraphBuilding::insertNonDuplicates(std::vector<Graph> &graphs, Graph graph)
{
    bool found=false;
    for(size_t l=0;l<graphs.size();++l)
        if(Graph::similar(graphs[l],graph))
        {
            found=true;
            break;
        }
    if(!found || GraphBuilding::m_keepDuplicateGraphs)
        graphs.push_back(graph);
}

std::vector<Graph> GraphBuilding::recursivelyUnifyEdges(size_t recursiveCall)
{
    m_countRecursives++;
    recursiveCall++;
    std::vector<size_t> positions = freeEdges();
    if(positions.size()==0)
    {
        std::vector<Graph> retGraph;
        retGraph.push_back(getGraph());
        return retGraph;
    }
    std::vector<Couple> possibleEdges = possibleUnifiedEdges();
    std::vector<Graph> combines;
    for(size_t i = 0;i<possibleEdges.size();++i)
    {
        std::pair<GraphBuilding,GraphBuilding> graphs = unifyEdge(possibleEdges[i]);
        std::vector<Graph> combines1 = graphs.first.recursivelyUnifyEdges(recursiveCall);
        std::vector<Graph> combines2 = graphs.second.recursivelyUnifyEdges(recursiveCall);
        if(combines1.size()!=0 && combines2.size()!=0)
        {
            for(size_t k = 0;k<combines1.size();++k)
                for(size_t j = 0;j<combines2.size();++j)
                {
                    Graph combine = GraphBuilding::combineGraphs(combines1[k],combines2[j]);
                    insertNonDuplicates(combines,combine);
                    //combines.push_back(combine);
                }
        }
        addAlreadyTried(possibleEdges[i]);
    }
    return combines;
}

bool GraphBuilding::nextInserts(gsVector<size_t> &inserts,size_t startSize)
{
    size_t size=startSize+2*(inserts.rows()-1);
    bool carryUp=true,carryUpHappened=false;
    for(index_t i = inserts.rows()-1;i>=0;i--)
    {
        if(carryUp)
        {
            inserts(i)+=1;
            carryUp=false;
        }
        if(inserts(i)>size-1&&i>0)
        {
            inserts(i)=0;
            carryUp=true;
            carryUpHappened=true;
        }
        size-=2;
    }
    if(inserts(0)>=startSize)
        return false;
    if(carryUpHappened)
        for(index_t i = inserts.rows()-1;i>=0;i--)
        {
            inserts(i)=inserts(i)>inserts(0)?inserts(i):inserts(0);
        }
    return true;
}

std::vector<Graph> GraphBuilding::findAllPossibleGraphs(size_t maxToInsertVertices, size_t lastInserted)
{
    if(maxToInsertVertices==0)
    {
        m_countTrees++;
        return recursivelyUnifyEdges();
    }
    std::vector<size_t> positions = possibleVertexLocations(lastInserted);
    GISMO_ASSERT(positions.size()>0,"GraphBuilding object with no free edges");
    std::vector<Graph> combines;
    for(size_t j = 0;j<positions.size();++j)
    {
        GraphBuilding gNew = getGraphBuildingWithAddedVertex(positions[j]);
        std::vector<Graph> combine = gNew.findAllPossibleGraphs(maxToInsertVertices-1,positions[j]);
        for(size_t k = 0;k<combine.size();++k)
        {
            Graph normalized = combine[k].normalizeGraph();
            insertNonDuplicates(combines,normalized);
        }
    }
    return combines;
}

void GraphBuilding::removeTypeII(size_t oldEdge)
{
    std::vector<size_t> toErase;
    for (size_t i = 0; i < m_typeII.size(); ++i)
    {
        if (m_typeII[i].first == oldEdge || m_typeII[i].second == oldEdge)
            toErase.push_back(i);
    }
    for (size_t i = toErase.size() - 1; i < -1UL; --i)
    {
        m_typeII.erase(m_typeII.begin() + toErase[i]);
    }
}

void GraphBuilding::convertTypeItoTypeII(size_t oldEdge,std::vector<size_t> &newEdges)
{
    std::vector<size_t> toErase;
    for (size_t i = 0; i < m_typeI.size(); ++i)
    {
        if (m_typeI[i].first == oldEdge)
        {
            for (size_t j = 0; j < newEdges.size(); ++j)
                m_typeII.push_back(std::make_pair(newEdges[j], m_typeI[i].second));
            toErase.push_back(i);
        }
        else if (m_typeI[i].second == oldEdge)
        {
            for (size_t j = 0; j < newEdges.size(); ++j)
                m_typeII.push_back(std::make_pair(m_typeI[i].first, newEdges[j]));
            toErase.push_back(i);
        }
    }
    for (size_t i = toErase.size() - 1; i < -1UL; --i)
    {
        m_typeI.erase(m_typeI.begin() + toErase[i]);
    }
}

}
