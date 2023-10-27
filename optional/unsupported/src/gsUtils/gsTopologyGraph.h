/** @file gsTopologyGraph.h

    @brief Classes for producing graphs related to multipatch layouts,
        details can be found in the paper
        F. Buchegger, B. Jüttler:
        Planar Multi-Patch Domain Parameterization via Patch Adjacency Graphs.
        NFN Technical Report No. 44
        http://www.gs.jku.at/pubs/NFNreport44.pdf

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsCore/gsBoundary.h>
#include <gsIO/gsFileData.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsCore/gsMultiPatch.h>

namespace gismo
{

class GISMO_EXPORT Graph {
public:
    typedef std::pair<size_t,size_t> Couple;
    typedef std::pair<patchSide,patchSide> PsCouple;

    Graph() : m_vertices(0),m_boundaryVertices(0) { }

    Graph(size_t vertices,size_t boundaryVert) :
        m_vertices(vertices),m_boundaryVertices(boundaryVert)
    { }

    index_t getNrVertices() const { return m_vertices; }

    index_t getNrBoundaryVertices() const { return m_boundaryVertices; }

    void addVertices(size_t vertices) { m_vertices+= vertices; }

    void addEdge(size_t v1,size_t v2,short_t pos1,short_t pos2);

    void addEdge(patchSide ps1,patchSide ps2);

    bool removeEdge(PsCouple ps1);

    std::vector<PsCouple> getEdges() const { return m_edges; }

    std::vector<size_t> vertexConnections(size_t vertex) const;

    void print() const
    {
        std::cout << "Graph: \nVertices: " << m_vertices << "\nEdges: ";
        for (size_t i = 0; i < m_edges.size(); ++i)
            std::cout << "(" << m_edges[i].first << "," << m_edges[i].second << ")";
        std::cout << std::endl;
    }

    static bool equal(Graph const & g1,Graph const & g2);

    static bool similar(const Graph & g1, const Graph & g2);

    // reordering of Graph: first vertex
    // is the one with connection to the lowest of smaller vertices, if
    // the same then to the one with highest amount of connections to smaller vertices,
    // if the same to the smaller sum of them
    Graph normalizeGraph() const;

    Graph rotateBoundaryByOne() const;

    bool innerGraphIsTwoEdgeConnected() const;

    bool innerGraphIsConnected() const;

    bool isSymmetric(std::vector<Couple> symmetricVertices) const;

    static void levelMatrixToGraphs(const gsMatrix<index_t>& mat,std::vector<Graph>& graphs);

    static void graphsToLevelMatrix(std::vector<Graph>& graphs,gsMatrix<index_t>& mat);

    static void getPartitionsOfGraphVector(std::vector<Graph>& graphs,std::vector<std::vector<Graph> >& graphPartitions,
                                    index_t minFaces = -1, index_t maxFaces = 1000000000);

    static void readGraphListFromFile(std::string filename,std::vector<Graph>& graphs);

    static void printGraphListToFile(std::string filename,std::vector<Graph>& graphs);

private:

    static short_t posToSide(short_t pos);

    size_t findSmallestVertex(index_t start) const;

    static void swapInVertices(std::vector<PsCouple> &edges,index_t swap1,index_t swap2);

    void putEdgesInGraph(std::vector<PsCouple> &edges,index_t limit);

    void travelGraph(size_t vertex,std::vector<bool>& visited) const;

private:

    index_t m_vertices;
    index_t m_boundaryVertices;
    std::vector<PsCouple> m_edges;
};

class GISMO_EXPORT GraphBuilding {
public:
    typedef Graph::Couple Couple;
    typedef Graph::PsCouple PsCouple;

    GraphBuilding() :
        m_origLimit(0),m_maxElement(0)
    {

    }

    GraphBuilding(size_t nrEdges,std::vector<Couple> &typeI);

    GraphBuilding(size_t nrEdges);

    void print() const
    {
        std::cout << "edges: ";
        for(size_t i = 0;i<m_edges.size();++i)
            std::cout << m_edges[i] << " ";
        std::cout << "\ntypeI: ";
        for(size_t i = 0;i<m_typeI.size();++i)
            std::cout << "(" << m_typeI[i].first << "," << m_typeI[i].second << ") ";
        std::cout << "\ntypeII: ";
        for(size_t i = 0;i<m_typeII.size();++i)
            std::cout << "(" << m_typeII[i].first << "," << m_typeII[i].second << ") ";
        std::cout << "\n";
    }

public :
    void setTypeIGroup(std::vector<size_t> &edges);

    void setTypeIIGroup(std::vector<size_t> &edges);

    void addAlreadyTried(Couple c) { insertNonDuplicatesOrdered(m_alreadyTried,c);  }

    // connecting edges (unifying)

    std::vector<Couple> possibleUnifiedEdges() const;

    bool possibleConnect(Couple c) const;

    bool isActiveEdge(size_t e) const;

    std::pair<GraphBuilding,GraphBuilding> unifyEdge(Couple c) const;

public:

    Graph getGraph() { return m_graph; }

    size_t getOrigLimit() const { return m_origLimit; }

    static Graph combineGraphs(Graph g1,Graph g2);

    // adding vertices:

    std::vector<size_t> freeEdges() const { return m_edges;  }

    std::vector<size_t> possibleVertexLocations(size_t lastInserted=-1) const;

    GraphBuilding getGraphBuildingWithAddedVertex(size_t connect) const;

    bool findVertexInRelations(size_t vert,size_t &rel, size_t &side) const;

    static void insertNonDuplicatesOrdered(std::vector<Couple> &couples,Couple insert);

    static void insertNonDuplicatesOrdered(std::vector<PsCouple> &couples,PsCouple insert);

    static void insertNonDuplicates(std::vector<Graph> &graphs, Graph graph);

    std::vector<Graph> recursivelyUnifyEdges(size_t recursiveCall=0);

    static bool nextInserts(gsVector<size_t> &inserts,size_t startSize);

    std::vector<Graph> findAllPossibleGraphs(size_t maxToInsertVertices,size_t lastInserted=-1);

    static void resetCounters()
    {
        m_countRecursives=0;
        m_countTrees=0;
    }

    static size_t getTreeCount() { return m_countTrees; }

    static size_t getRecursiveCount() { return m_countRecursives; }

    static void setKeepDuplicates(bool keepDuplicatesGraphs) { m_keepDuplicateGraphs=keepDuplicatesGraphs; }

private:

    void removeTypeII(size_t oldEdge);

    void convertTypeItoTypeII(size_t oldEdge,std::vector<size_t> &newEdges);

    static bool sortCouples (Couple i,Couple j) // to do: free functions as (private) members of class Graph/GraphBuilding
    {
        return ( i.first<j.first ||
                 (i.first==j.first && i.second<j.second) );
    }

    static void sortEdges(std::vector<PsCouple> &couples)
    {
        std::sort (couples.begin(), couples.end(), sortPsCouples);
    }

    static bool sortPsCouples (PsCouple i,PsCouple j)
    {
        return ( i.first.patch<j.first.patch ||
                 (i.first.patch==j.first.patch && i.second.patch<j.second.patch) );
    }

    static void sortEdges(std::vector<Couple> &couples)
    {
        std::sort (couples.begin(), couples.end(), sortCouples);
    }

private:
    std::vector<size_t> m_edges;
    std::vector<size_t> m_origEdges;
    std::vector<Couple> m_typeI;
    std::vector<Couple> m_typeII;
    std::vector<Couple> m_alreadyTried;
    size_t m_origLimit;
    size_t m_maxElement;
    Graph m_graph;
    std::vector<std::vector<size_t> > m_edgeToVertexRelations;

    static size_t m_countRecursives;
    static size_t m_countTrees;
    static bool m_keepDuplicateGraphs;
};

class ConstructMultiPatchFromGraph {

    typedef Graph::Couple Couple;
    typedef Graph::PsCouple PsCouple;

public:

    static gsBoxTopology getBoxTopology(Graph g)
    {
        size_t boundarySize = g.getNrBoundaryVertices();
        size_t boxes = g.getNrVertices()-boundarySize;
        std::vector<boundaryInterface> interfaces;
        for(size_t i = 0;i<g.getEdges().size();++i)
        {
            size_t patch1=g.getEdges()[i].first.patch;
            size_t patch2=g.getEdges()[i].second.patch;
            boxSide side1=g.getEdges()[i].first.side();
            boxSide side2=g.getEdges()[i].second.side();
            if(patch1>=boundarySize&&patch2>=boundarySize)
            {
                bool orientation = true;
                if( (side1==1&&side2==1) ||
                        (side1==1&&side2==4) ||
                        (side1==2&&side2==2) ||
                        (side1==2&&side2==3) ||
                        (side1==3&&side2==2) ||
                        (side1==3&&side2==3) ||
                        (side1==4&&side2==1) ||
                        (side1==4&&side2==4) )
                {
                    orientation=false;
                }
                patchSide ps1=g.getEdges()[i].first;
                ps1.patch-=boundarySize;
                patchSide ps2=g.getEdges()[i].second;
                ps2.patch-=boundarySize;
                boundaryInterface interface(ps1,ps2,orientation);
                interfaces.push_back(interface);
            }

        }
        std::vector<patchSide> boundaries;
        gsBoxTopology topol(2,boxes,boundaries,interfaces);
        topol.addAutoBoundaries();
        return topol;
    }

    static gsGeometry<real_t>* getPatch(size_t nrOfControlPointsInOneRow)
    {
        size_t deg=2;
        gsKnotVector<real_t> kv(0,1,nrOfControlPointsInOneRow-deg-1,deg+1);
        gsTensorBSplineBasis<2,real_t> basis = gsTensorBSplineBasis<2,real_t>(kv,kv);
        gsMatrix<real_t> coefs;
        coefs.setConstant(nrOfControlPointsInOneRow*nrOfControlPointsInOneRow,2,100);
        gsGeometry<real_t>* patch = new gsTensorBSpline<2,real_t>(basis,coefs);
        return patch;
    }

    static void setBoundaryCoefs(gsMatrix<real_t> &allCoefs,gsMatrix<real_t> bound,boxSide side)
    {
        switch(side)
        {
        case 1:
            for(index_t i = 0;i<bound.rows();++i)
            {
                    allCoefs(i*bound.rows(),0)=bound(bound.rows()-i-1,0);
                    allCoefs(i*bound.rows(),1)=bound(bound.rows()-i-1,1);
            }
            break;
        case 2:
            for(index_t i = 0;i<bound.rows();++i)
            {
                    allCoefs((i+1)*bound.rows()-1,0)=bound(i,0);
                    allCoefs((i+1)*bound.rows()-1,1)=bound(i,1);
            }
            break;
        case 3: allCoefs.block(0,0,bound.rows(),2)=bound;
                break;
        case 4:
            for(index_t i = 0;i<bound.rows();++i)
            {
                    allCoefs(allCoefs.rows()-bound.rows()+i,0)=bound(bound.rows()-i-1,0);
                    allCoefs(allCoefs.rows()-bound.rows()+i,1)=bound(bound.rows()-i-1,1);
            }
            break;
        default:
            GISMO_ERROR("Should have valid side.");
            break;
        }
    }

    static void setCornerCoefs(gsMatrix<real_t>& coefs,boxCorner corner,gsMatrix<real_t> cornerCoefs)
    {
        size_t rowlength = math::isqrt((size_t)coefs.rows()); //assume that same number of coefs in u and v
        switch(corner)
        {
        case 1:
            coefs.block(0,0,1,2) = cornerCoefs;
            break;
        case 2:
            coefs.block(rowlength-1,0,1,2) = cornerCoefs;
            break;
        case 3:
            coefs.block(coefs.rows()-rowlength,0,1,2) = cornerCoefs;
            break;
        case 4:
            coefs.block(coefs.rows()-1,0,1,2) = cornerCoefs;
            break;
           default:
            GISMO_ERROR("Should have valid corner.");
            break;
        }
    }

    static gsMatrix<real_t> extractCornerCoefs(gsMatrix<real_t>& coefs,boxCorner corner)
    {
        gsMatrix<real_t> result(1,2);
        size_t rowlength = math::isqrt((size_t)coefs.rows()); //assume that same number of coefs in u and v
        switch(corner)
        {
        case 1:
            result = coefs.block(0,0,1,2);
            break;
        case 2:
            result = coefs.block(rowlength-1,0,1,2);
            break;
        case 3:
            result = coefs.block(coefs.rows()-rowlength,0,1,2);
            break;
        case 4:
            result = coefs.block(coefs.rows()-1,0,1,2);
            break;
           default:
            GISMO_ERROR("Should have valid corner.");
            break;
        }
        return result;
    }

    static void setAllCornerCoefs(gsMultiPatch<real_t>& mp,patchCorner start,std::vector<patchCorner> cornerList)
    {
        gsMatrix<real_t> coefs = extractCornerCoefs(mp.patch(start.patch).coefs(),start);
        for(size_t i  = 0;i<cornerList.size();++i)
        {
            setCornerCoefs(mp.patch(cornerList[i].patch).coefs(),cornerList[i],coefs);
        }
    }

    static void graphToMP(Graph g,std::vector<gsMatrix<real_t> >& boundary,gsMultiPatch<real_t>& result)
    {
        gsBoxTopology topol=getBoxTopology(g);
        std::vector<gsGeometry<real_t>* > patches;
        for(index_t i = 0;i<topol.nBoxes();++i)
            patches.push_back(getPatch(boundary[0].rows()));
        for(gsBoxTopology::const_biterator b = topol.bBegin(); b!=topol.bEnd();b++)
        {
            size_t patch=b->patch;
            boxSide side=b->side();
            setBoundaryCoefs(patches[patch]->coefs(),boundary[patch+boundary.size()],side);
        }
        result=gsMultiPatch<real_t>(patches,topol.boundaries(),topol.interfaces());

        std::vector<patchSide> bounds = result.boundaries();
        for(size_t k = 0;k<bounds.size();++k)
        {
            std::vector<boxCorner> boxCorners;
            bounds[k].side().getContainedCorners(2,boxCorners);
            for(size_t j=0;j<boxCorners.size();++j)
            {
                patchCorner c(bounds[k].patch,boxCorners[j]);
                std::vector<patchCorner> localCornerList;
                result.getCornerList(c,localCornerList);
                setAllCornerCoefs(result,c,localCornerList);
            }
        }
        // go through all corner of the boundary and sanitize the coefs
    }

    static void graphsToMPs(std::vector<Graph> gs,std::vector<gsMatrix<real_t> >& boundary,std::vector<gsMultiPatch<real_t> >& results)
    {
        results.clear();
        for(size_t i = 0;i<gs.size();++i)
        {
            gsMultiPatch<real_t> mp;
            graphToMP(gs[i],boundary,mp);
            results.push_back(mp);
        }
    }

};





} // namespace gismo
