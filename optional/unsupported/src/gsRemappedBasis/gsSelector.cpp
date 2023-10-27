/** @file gsSelector.cpp

    @brief Implementation of non templated parts of gsSelector

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
*/

#include <gsRemappedBasis/gsSelector.h>

namespace gismo
{

std::ostream& operator<<  (std::ostream &out, const NodeData &node)
{
    out<<"as fork: (dir: "<<node.dir<<", par: "<<node.par<< "); "
      <<"as leaf: (space: "<<node.space<< ")\n";
    return out;
}


std::ostream& operator<<  (std::ostream &out, const gsBoxList &boxes)
{
    GISMO_ASSERT(boxes.domainDim()==2,"can print only 2d partitions");
    for(size_t b=0; b<boxes.size();++b)
    {
        gsMatrix<real_t> center=(boxes.box(b).col(0)+boxes.box(b).col(1))/2;
        out<<"\\draw[fill=black!"<< 10*boxes.basisId(b)%100<<"!white ] ("<<boxes.box(b)(0,0)<<','<<boxes.box(b)(1,0)
          <<") rectangle ("<<boxes.box(b)(0,1)<<','<<boxes.box(b)(1,1)<<") "
         <<"("<<center(0,0)<<','<<center(1,0)<<") node {\\tiny"<<boxes.basisId(b)<<"};\n";
    }
    return out;
}

// Initializers

void gsSelector::initFromBoxes (
        const gsBoxList        &boxes,
        const gsMatrix<real_t> &boundingBox,
        size_t patch)
{
    if (m_domainMap.size()<=patch)
        m_domainMap.resize(patch+1);
    m_domainMap[patch].initFromBoxes(boxes,boundingBox);
}

void gsSelector::initFromBoxesMax (
        const gsBoxList        &boxes,
        const gsMatrix<real_t> &boundingBox,
        size_t patch)
{
    if (m_domainMap.size()<=patch)
        m_domainMap.resize(patch+1);
    m_domainMap[patch].initFromBoxesMax(boxes,boundingBox);
}

void gsSelector::initFromBoxesMin (
        const gsBoxList        &boxes,
        const gsMatrix<real_t> &boundingBox,
        size_t patch)
{
    if (m_domainMap.size()<=patch)
        m_domainMap.resize(patch+1);
    m_domainMap[patch].initFromBoxesMin(boxes,boundingBox);
}

void gsSelector::addBoxesMax(
        const gsBoxList        &boxes,
        size_t patch
        )
{
    if (m_domainMap.size()<=patch)
        m_domainMap.resize(patch+1);
    m_domainMap[patch].addBoxesMax(boxes);
}

void gsSelector::addBoxesMin(
        const gsBoxList        &boxes,
        size_t patch
        )
{
    if (m_domainMap.size()<=patch)
        m_domainMap.resize(patch+1);
    m_domainMap[patch].addBoxesMin(boxes);
}


// Getters

gsDomainMap&       gsSelector::patch(size_t p)
{
    if (p>=m_domainMap.size())
        m_domainMap.resize(p+1);
    return m_domainMap[p];
}

const gsDomainMap& gsSelector::patch(size_t p) const
{
    return m_domainMap[p];
}

size_t gsSelector::patchNum () const
{
    return m_domainMap.size();
}


/**
     * Returns the basisId of the basis capable of depicting all the basis functions
     * active in the given points.
     * @param points , each given as a column of the matrix;
     * @param out returned basisId's (one per point).
     */
void gsSelector::getBasisAt ( const gsMatrix<real_t> &points,
                              std::vector<basisIdT>  &out,
                              patchIdT               patch ) const
{
    return m_domainMap[patch].getBasisAt(points,out);
}

// Exports

gsBoxList gsSelector::asBoxList(size_t  patch ) const
{
    return m_domainMap[patch].asBoxList();
}


void gsSelector::exportToTex(std::string filename) const
{
    const size_t numPatches=m_domainMap.size();
    std::ofstream fout;
    fout.open((filename + ".tex").c_str());
    if( fout.is_open() )
    {

        fout<<"\\documentclass[tikz,multi]{standalone}\n"
              "\\usepackage{tikz}\n"
              "\\begin{document}\n";
        for (size_t patch=0; patch<numPatches; ++patch)
            fout<<"\\begin{tikzpicture}[scale=20]\n"
               <<m_domainMap[patch].asBoxList()
              <<"\\end{tikzpicture}"
             <<(patch==numPatches-1 ? "\n\\end{document}\n" : "\n\\newpage\n\n" )
            <<std::endl;
        fout.close();
    }
    else
        gsWarn << "Error opening file, output not saved!";
}

} // namespace gismo
