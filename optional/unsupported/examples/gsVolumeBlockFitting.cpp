// Author: Jaka Speh
// Generates the volumes and impose the C0 continuity.

#include <cmath>

#include <iostream>
#include <fstream>

#include <gismo.h>
#include <gsIO/gsIOUtils.h>


#ifdef GISMO_WITH_PSOLID
#include <gsParasolid/gsWriteParasolid.h>
using gismo::extensions::gsWriteParasolid;
#endif

using namespace gismo;

static const real_t PI = 3.14159265358979323846;


// a control point match is a list of pairs, where first = id of patch, and
// second = index of control pt
//typedef std::vector<std::pair<int, int> > CpMatchType;


struct CpMatchEntry {
    int patchId;
    int ctrlPtId;
    boxSide bs;
};

typedef std::vector< CpMatchEntry > CpMatchType;

// auxilliary function that converts an integer to a string
std::string int2str(int number)
{
    std::ostringstream convert;
    convert << number;
    return convert.str();
}


int nthIntRoot(int base, int n)
{
    if(n < 1)
        return -1;
    int i = 1;
    int val = math::ipow(i,n);
    while(val < base)
    {
        i++;
        val = math::ipow(i,n);
    }
    if(val == base)
        return i;
    return -1;
}


/// Get the information about incident faces.
/// \param pairs consists of vectors of size four
///              [ [idVolume1, idFace1, idVolume2, idFace2],
///                [ ...]
///                ...
///              ]
///              which means, that volumes[idVolume1] idFace1-th face is
///              incident with volumes[idVolume2] idFace2-th face
/// \param fileName name of the file where the information about pairs is
///                 stored
void getIncidentFaces(std::vector<std::vector<unsigned> >& pairs,
                      const std::vector<std::string> files,
                      const std::string& fileName)
{
    std::ifstream pairFile;
    pairFile.open(fileName.c_str());

    // file name of volume1 and file name of volume2
    std::string f1, f2;

    int counterPair = 0;
    while (pairFile >> f1)
    {
        pairs.push_back(std::vector<unsigned>());

        f1 += ".xml"; // fix I am comparing file names (with extension

        unsigned idFace1;
        pairFile >> idFace1;

        unsigned idf1 = 100;
        for (size_t index = 0; index < files.size(); index++)
        {
            if (files[index] == f1)
            {
                idf1 = index;
                break;
            }
        }

        pairs[counterPair].push_back(idf1);
        pairs[counterPair].push_back(idFace1);

        pairFile >> f2;

        f2 += ".xml";

        unsigned idFace2;
        pairFile >> idFace2;

        unsigned idf2 = 100;
        for (size_t index = 0; index < files.size(); index++)
        {
            if (files[index] == f2)
            {
                idf2 = index;
                break;
            }
        }

        pairs[counterPair].push_back(idf2);
        pairs[counterPair].push_back(idFace2);

        counterPair++;
    }
}

/// Returns four corner indices of a bspline volume from the boxSide side.
///
/// \param indices vector of output indices
/// \param bspline input volume
/// \param side boundary side
void getCornerBoundaryIndices(std::vector<gsVector<index_t, 3> >& indices,
                              const gsTensorBSpline<3, real_t>& bspline,
                              const boxSide side)
{
    int fixedDirection = side.direction();
    unsigned fixedIndex = side.parameter() * (bspline.basis().size(fixedDirection) - 1);


    int direction1 = (fixedDirection == 0) ? 1 : 0;
    int direction2 = (fixedDirection == 2) ? 1 : 2;

    gsVector<index_t, 3> tensorIndex(3);
    tensorIndex.setZero();
    tensorIndex(fixedDirection) = fixedIndex;

    indices.push_back(tensorIndex);

    tensorIndex(direction1) = bspline.basis().size(direction1) - 1;

    indices.push_back(tensorIndex);

    tensorIndex(direction1) = 0;
    tensorIndex(direction2) = bspline.basis().size(direction2) - 1;

    indices.push_back(tensorIndex);

    tensorIndex(direction1) = bspline.basis().size(direction1) - 1;

    indices.push_back(tensorIndex);

}

/// Compute the absolute distance between two (1 x 3) row matrices
/// \param mat1 (1 x 3) matrix
/// \param mat2 (1 x 3) matrix
real_t absoluteDistance(const gsMatrix<> mat1,
                const gsMatrix<> mat2)
{
    real_t diff = 0;
    for (unsigned col = 0; col < 3; col++)
    {
        diff += math::abs(mat1(0, col) - mat2(0, col));
    }

    return diff;
}

/// Gets two vectors, it computes the difference between these two vectors and
/// normalize it. For example
/// start = (10, 5, 0)
/// end = (5, 10, 0)
///
/// end - start = (-5, 5, 0)
///
/// result = (-1, 1, 0)
///
/// \param start vector
/// \param end vector
gsVector<int, 3> getStandartDirectionalVector(const gsVector<index_t, 3>& start,
                                              const gsVector<index_t, 3>& end)
{
    gsVector<int, 3> vec;

    for (int index = 0; index < 3; index++)
    {
        int diff = static_cast<int>(end(index)) - static_cast<int>(start(index));

        if (diff != 0)
        {
            diff = (diff < 0) ? -1: 1;
        }

        vec(index) = diff;
    }

    return vec;
}


/// Cast int vector to unsigned vector.
///
/// \param vec int vector
/// \return unsigned vector
gsVector<unsigned, 3> intToUnsigned(gsVector<int, 3> vec)
{
    gsVector<unsigned, 3> newVec;

    for (int row = 0; row < 3; row++)
    {
        newVec(row) = static_cast<unsigned>(vec(row));
    }

    return newVec;
}


/// Cast unsigned vector to int vector.
///
/// \param vec unsigned vector
/// \return int vector
gsVector<int, 3> unsignedToInt(gsVector<unsigned, 3> vec)
{
    gsVector<int, 3> newVec;

    for (int row = 0; row < 3; row++)
    {
        newVec(row) = static_cast<int>(vec(row));
    }

    return newVec;
}

//if controlpoint is on the boundary of the surface, return 0 (boxSide 0...none)
boxSide setSide(CpMatchEntry &m, gsTensorBSpline<3, real_t> const &vol, boxSide side)
{
    index_t ni = vol.basis().size(0);
    index_t nj = vol.basis().size(1);
    index_t nk = vol.basis().size(2);

    gsVector<index_t,3> ind = vol.basis().tensorIndex(m.ctrlPtId);

    //if controlpoint is on west or east and on the edge of the surface, set side to 0
    if(side == 1 || side == 2)
    {
        if(ind[1] == 0 || ind[1] == nj-1 || ind[2] == 0 || ind[2] == nk-1)
            return 0;
    }
    //if controlpoint is on north or south and on the edge of the surface, set side to 0
    if(side == 3 || side == 4)
    {
        if(ind[0] == 0 || ind[0] == ni-1 || ind[2] == 0 || ind[2] == nk-1)
            return 0;
    }
    //if controlpoint is on front or back and on the edge of the surface, set side to 0
    if(side == 5 || side == 6)
    {
        if(ind[0] == 0 || ind[0] == ni-1 || ind[1] == 0 || ind[1] == nj-1)
            return 0;
    }
    return side;
}


void addMatch(std::vector<CpMatchType> &matches,
              int patch1, int cpIdx1,
              int patch2, int cpIdx2,
              boxSide side1, boxSide side2,
              gsTensorBSpline<3, real_t> const &vol1,
              gsTensorBSpline<3, real_t> const &vol2)
{
    // search for the points
    size_t nMatches = matches.size();
    int matchIdx1 = -1;
    int matchIdx2 = -1;
    for(size_t i = 0; i < nMatches; i++)
    {
        size_t nPts = matches[i].size();
        for(size_t j = 0; j < nPts; j++)
        {
            CpMatchEntry curMatch = matches[i][j];
            if(curMatch.patchId == patch1 && curMatch.ctrlPtId == cpIdx1)
            {
                GISMO_ASSERT(matchIdx1 == -1, "Point in multiple matches");
                matchIdx1 = (int)i;
            }
            if(curMatch.patchId == patch2 && curMatch.ctrlPtId == cpIdx2)
            {
                GISMO_ASSERT(matchIdx2 == -1, "Point in multiple matches");
                matchIdx2 = (int)i;
            }
        }
    }
    
    if(matchIdx1 != -1 && matchIdx2 != -1 && matchIdx1 != matchIdx2)
    {
        // both points seen, merge matches
        matches[matchIdx1].insert(matches[matchIdx1].end(),
                                  matches[matchIdx2].begin(),
                                  matches[matchIdx2].end());
        matches.erase(matches.begin() + matchIdx2);
    }
    else if(matchIdx1 != -1 && matchIdx2 != -1)
    {
        // points already seen, already matched
    }
    else if(matchIdx1 != -1)
    {
        // one point seen, add the other one
        CpMatchEntry m;
        m.patchId = patch2;
        m.ctrlPtId = cpIdx2;
        m.bs = side2;
        m.bs = setSide(m, vol2, m.bs);
        matches[matchIdx1].push_back(m);
    }
    else if(matchIdx2 != -1)
    {
        // one point seen, add the other one
        CpMatchEntry m;
        m.patchId = patch1;
        m.ctrlPtId = cpIdx1;
        m.bs = side1;
        m.bs = setSide(m, vol1, m.bs);
        matches[matchIdx2].push_back(m);
    }
    else
    {
        // neither point seen, make a new match containing both
        CpMatchType newMatch;
        CpMatchEntry m1,m2;

        m1.patchId = patch1;
        m1.ctrlPtId = cpIdx1;
        m1.bs = side1;
        m1.bs = setSide(m1, vol1, m1.bs);
        m2.patchId = patch2;
        m2.ctrlPtId = cpIdx2;
        m2.bs = side2;
        m2.bs = setSide(m2, vol2, m2.bs);

        newMatch.push_back(m1);
        newMatch.push_back(m2);

        matches.push_back(newMatch);
    }
}

/// Volume1 side1 almost (C0) matches volume2 side2. We move the control
/// points of both splines such that they (C0) match along the boundary.
///
/// \param volume1 first tensor bspline volume
/// \param side1 side from the volume1
/// \param volume2 second tensor bspline volume
/// \param side2 side from the volume2
void matchBoundaryOfVolumes(
                            int volIndex1,
                            gsTensorBSpline<3, real_t>& volume1,
                            boxSide side1,
                            int volIndex2,
                            gsTensorBSpline<3, real_t>& volume2,
                            boxSide side2,
                            std::vector<CpMatchType> &matches)
{

    std::vector<gsVector<index_t, 3> > ind1, tmpInd2, ind2;

    // getting 4 corner indices of side1 of volumei, i = 1, 2
    getCornerBoundaryIndices(ind1, volume1, side1);
    getCornerBoundaryIndices(tmpInd2, volume2, side2);

    // matching 4 corner indices
    for (unsigned i = 0; i < 4; i++)
    {
        unsigned coefIndex1 = volume1.basis().index(ind1[i]);

        gsMatrix<> coef1 = volume1.coef(coefIndex1);

        // distance between fixed coefficient from first volume to all four
        // boundary vertices from second volume
        gsVector<real_t, 4> distance;

        for (unsigned j = 0; j < 4; j++)
        {
            unsigned coefIndex2 = volume2.basis().index(tmpInd2[j]);
            gsMatrix<> coef2 = volume2.coef(coefIndex2);

            real_t dst = absoluteDistance(coef1, coef2);
            distance(j) = dst;
        }

        int minIndex;
        distance.minCoeff(&minIndex);

        ind2.push_back(tmpInd2[minIndex]);
    }


    // directional vectors, we can access all boundary indices via these vectors
    //
    // boundary1 point (i, j) has index
    // ind1[0] + i * vec1u + j * vec1v
    //
    // similar for boundary2

    gsVector<index_t, 3> vec1u = getStandartDirectionalVector(ind1[0], ind1[1]);
    gsVector<index_t, 3> vec1v = getStandartDirectionalVector(ind1[0], ind1[2]);

    gsVector<index_t, 3> vec2u = getStandartDirectionalVector(ind2[0], ind2[1]);
    gsVector<index_t, 3> vec2v = getStandartDirectionalVector(ind2[0], ind2[2]);


//    gsInfo << "vec1u: " << vec1u << "\n"
//              << "vec1v: " << vec1v << "\n\n"
//              << "vec2u: " << vec2u << "\n"
//              << "vec2v: " << vec2v << "\n";

    // FIX THIS, NOW IT IS WORKING BECAUSE ALL DIRECTIONS ARE THE SAME

    // number of basis function in direction vec1u
    int dir;
    vec1u.maxCoeff(&dir); // should have just one non-zero coefficient
    int size1u = volume1.basis().size(dir);
    vec1v.maxCoeff(&dir);
    int size1v = volume1.basis().size(dir);

    vec2u.maxCoeff(&dir);
    int size2u = volume2.basis().size(dir);
    vec2v.maxCoeff(&dir);
    int size2v = volume2.basis().size(dir);


    if (!(size1u == size2u) && (size1v == size2v))
    {
        GISMO_ERROR("Sizes of boundaries are not the same!");
    }

    real_t maxDst = 0;

    // we will change the coefficients
    gsMatrix<>& coefs1 = volume1.coefs();
    gsMatrix<>& coefs2 = volume2.coefs();

    for (int i = 0; i < size1u; i++)
    {
        for (int j = 0; j < size1v; j++)
        {
            // boundary index 1
            gsVector<index_t, 3> tmp1 = (ind1[0]) + i * vec1u + j * vec1v;
            gsVector<index_t, 3> tmp2 = (ind2[0]) + i * vec2u + j * vec2v;
            gsVector<index_t, 3> bIndex1 = (tmp1);
            gsVector<index_t, 3> bIndex2 = (tmp2);

            index_t coefIndex1 = volume1.basis().index(bIndex1);
            index_t coefIndex2 = volume2.basis().index(bIndex2);


            real_t dst = absoluteDistance(coefs1.row(coefIndex1),
                                          coefs2.row(coefIndex2));

            if (maxDst < dst)
            {
                maxDst = dst;
            }
            
            addMatch(matches, volIndex1, coefIndex1, volIndex2, coefIndex2, side1, side2, volume1, volume2);

            //gsMatrix<T> mid = coefs1.row(coefIndex1);
            // not good
            // (coefs1.row(coefIndex1) + coefs2.row(coefIndex2)) / 2;

            //coefs1.row(coefIndex1) = mid;
            //coefs2.row(coefIndex2) = mid;
       
        }
    }

    gsInfo << "Maximal distance (error): \n"
                 "    - before matching: "
              << maxDst << "\n";
}


/// Saves the output in various formats.
///
/// \param volumes the volumes we want to save
/// \param output where we want to save the volumes
void saveVolumes(const std::vector<gsTensorBSpline<3, real_t> >& volumes,
                 const std::string& output,
                 const bool writeGoToolsFile = false)
{

    gsMultiPatch<> mp;

    for (size_t index = 0; index < volumes.size(); index++)
    {
        gsGeometry<>::uPtr geom(new gsTensorBSpline<3, real_t>(volumes[index]));
        mp.addPatch(give(geom));


        // writing each volume indidividually

        std::string prefix = output + util::to_string(index);
        // std::string strVol = prefix + "volume";
        // gsWriteParaview(volumes[index], strVol);


        // writing xml file for each volume
         std::string strXml = prefix + "volume.xml";
         gsFileData<> fileVolume;
         fileVolume << volumes[index];
         fileVolume.dump(strXml);


        // writing mesh for each volume
         gsMesh<> mesh;
         makeMesh(volumes[index].basis(), mesh, 10);
         volumes[index].evaluateMesh(mesh);
         std::string strMesh = prefix + "mesh";
         gsWriteParaview(mesh, strMesh);
    }

    gsWriteParaview(mp, output);

    gsFileData<> fileVolume;
    fileVolume << mp;
    fileVolume.dump(output);

    if (writeGoToolsFile)
    {
        std::string strGoTools = output + "GoTools.g2";

        std::ofstream file;
        file.open(strGoTools.c_str());

        if (!file.is_open())
        {
            gsWarn << "Can not open file: " << strGoTools << "\n"
                      "Aborting ...";
            return;
        }

        for (size_t index = 0; index < volumes.size(); index++)
        {
            gsWriteGoToolsSpline<3, real_t>(volumes[index], file);
        }

        file.close();
    }
}

/// Saves the output in various formats.
///
/// \param volumes the volumes we want to save
/// \param output where we want to save the volumes
void saveFaces(const std::vector<gsTensorBSpline<2, real_t> >& faces,
                 const std::string& output)
{

    gsMultiPatch<> mp;

    for (size_t index = 0; index < faces.size(); index++)
    {
        if(faces[index].coefsSize() != 0)
        {
            gsGeometry<>::uPtr geom(new gsTensorBSpline<2, real_t>(faces[index]));
            mp.addPatch(give(geom));


            // writing each topFace indidividually

            std::string prefix = output + util::to_string(index);

            // writing xml file for each volume
             std::string strXml = prefix + "TopFace.xml";
             gsFileData<> fileFace;
             fileFace << faces[index];
             fileFace.dump(strXml);


        // writing mesh for each volume
//         gsMesh<> mesh;
//         makeMesh(faces[index].basis(), mesh, 10);
//         faces[index].evaluateMesh(mesh);
//         std::string strMesh = prefix + "mesh";
//         gsWriteParaview(mesh, strMesh);
        }
    }

    //gsWriteParaview(mp, output);

    gsFileData<> fileFace;
    fileFace << mp;
    fileFace.dump(output+"TopFace");
}


/// From coefs that represents a single leg of the chair stand example, we
/// construct all other 4 legs with rotation of the coefficients.
///
/// \param coefs input coefficient for the chair stand example
/// \param i which leg we want to compute
/// \param newCoefs coefficients of the i-th leg
///
/// TO DO: rewrite this code with G+SMO geometry methods
void transformCoefs(const gsMatrix<>& coefs,
                    const int i,
                    gsMatrix<>& newCoefs)
{
    newCoefs.resize(coefs.rows(), coefs.cols());

    gsVector<> vec(3);
    vec << -0.2, 0.5, 0;

    gsMatrix<> mat(3, 3);
    mat << math::cos(2 * i * PI / 5), math::sin(2 * i * PI / 5), 0,
           -math::sin(2 * i * PI / 5), math::cos(2 * i * PI / 5), 0,
            0, 0, 1;

    for (int row = 0; row < coefs.rows(); row++)
    {
        gsVector<> point = coefs.row(row);
        newCoefs.row(row) = (mat * (point - vec)) + vec;
    }
}


/// get all the files from the pairFile
std::vector<std::string> getAllFileNames(std::string& pairFileName)
{
    std::vector<std::string> fileNames;

    std::string pairFilePath = gsFileManager::find(pairFileName);
    if (pairFilePath.empty())
    {
        gsWarn << "Cannot find file " << pairFileName << ".\n";
        return fileNames;
    }
    
    std::ifstream pairFile;
    pairFile.open(pairFilePath.c_str());

    // file name of volume1 and file name of volume2
    std::string f1, f2;


    std::vector<std::string> tmp;

    while (pairFile >> f1)
    {
        unsigned dontNeed;
        pairFile >> dontNeed;

        pairFile >> f2;
        pairFile >> dontNeed;

        f1 += ".xml";
        f2 += ".xml";

        bool alreadyInside = false;
        for (size_t index = 0; index < fileNames.size(); index++)
        {
            if (fileNames[index] == f1)
            {
                alreadyInside = true;
                break;
            }
        }

        if (!alreadyInside)
        {
            fileNames.push_back(f1);
        }

        tmp.push_back(f2);
    }

    for (size_t index1 = 0; index1 < tmp.size(); index1++)
    {
        bool alreadyInside = false;
        for (size_t index2 = 0; index2 < fileNames.size(); index2++)
        {
            if (tmp[index1] == fileNames[index2])
            {
                alreadyInside = true;
                break;
            }
        }

        if (!alreadyInside)
        {
            fileNames.push_back(tmp[index1]);
        }
    }

    return fileNames;
}


/// This function is usefull only in chair stand example (for Siemens people,
/// they want some crazy stuff).
gsTensorBSpline<3, real_t> constructVolume(gsTensorBSpline<3, real_t>& volume,
                                           const boxSide& side,
                                           int numInternalKnots,
                                           int degree,
                                           gsVector<>& axis,
                                           int fraction,
                                           bool turnAround = true)
{
    gsTensorBSpline<2, real_t>* surf =
            dynamic_cast<gsTensorBSpline<2, real_t>*> (volume.boundary(side).release());

    gsKnotVector<> kv(0, 1, numInternalKnots, degree + 1, 1);
    int numBasisFun = degree + 1 + numInternalKnots;

    const gsMatrix<>& coefs = surf->coefs();


    //const distance = (coefs.row(0).transpose() - axis).norm() * percent;

    gsMatrix<> newCoefs(coefs.rows() * numBasisFun, coefs.cols());

    int newRow = 0;
    for (int n = 0; n < numBasisFun; ++n)
    {
        const real_t factor = static_cast<real_t>(n) / (fraction * (numBasisFun - 1));


        for (int row = 0; row < coefs.rows(); ++row)
        {
            if (n == 0)
            {
                newCoefs.row(newRow) = coefs.row(row);
            }
            else
            {
                real_t coefX = coefs(row, 0);
                real_t coefY = coefs(row, 1);


                real_t x = coefX + factor * (axis(0) - coefX);
                real_t y = coefY + factor * (axis(1) - coefY);



//                if (row == 0)
//                {
//                gsInfo << "coefX: " << coefX << "  axis(0) " << axis(0)
//                          << "  x: " << x << "\n"
//                          << "coefX: " << coefY << "  axis(1) " << axis(1)
//                          << "  y: " << y << "\n" << "\n";
//                }

                newCoefs(newRow, 0) = x;
                newCoefs(newRow, 1) = y;
                newCoefs(newRow, 2) = coefs(row, 2);

            }

            ++newRow;
        }
    }

    gsMatrix<> mat(newCoefs.rows(), newCoefs.cols());
    const int maxX = surf->basis().component(0).size();
    const int maxY = surf->basis().component(1).size();

    if (turnAround)
    {
        int row = 0;
        const int surfPts = maxX * maxY; // points in surface

        for (int z = 0; z < numBasisFun; ++z)
        {
            for (int y = 0; y < maxY; ++y)
            {
                for (int x = 0; x < maxX; ++x)
                {
                    mat.row(row) = newCoefs.row((numBasisFun - 1 - z) * surfPts +
                                                y * maxX + x);
                    ++row;
                }

            }
        }
    }
    else
    {
        mat = newCoefs;
    }


    gsTensorBSpline<3, real_t> bspl(surf->basis().component(0).knots(),
                                    surf->basis().component(1).knots(),
                                    kv,
                                    give(mat));

    delete surf;

    return bspl;
}


/// Save the boundary of volumes into .g2 and .xmt_txt format. It takes all
/// faces of all volumes from vector volumes, but omits faces that are stored
/// in vector middleFaces. String output determines where the files will be stored.
///
/// \param volumes vector of volumes
/// \param middleFaces vectors of [volumes id, side] we must omit
/// \param output output string
void saveBoundary(const std::vector<gsTensorBSpline<3, real_t> >& volumes,
                  const std::vector<std::pair<unsigned, boxSide> >& middleFaces,
                  const std::string& output)
{

    const std::string prefix = output + "_boundary_";
    int counter = 0;

    gsMultiPatch<> multiPatch;
    for (unsigned index = 0; index < volumes.size(); index++)
    {
        for (boxSide side = boxSide::getFirst(3); boxSide::getEnd(3); ++side )
        {
            bool skip = false;
            for (size_t i = 0; i < middleFaces.size(); i++)
            {
                if (middleFaces[i].first == index)
                {
                    if (side == middleFaces[i].second)
                    {
                        skip = true;
                    }
                }
            }

            if (skip)
            {
                continue;
            }

            gsGeometry<>::uPtr g = volumes[index].boundary(side);

            //std::string name = prefix + util::to_string(counter);
            counter++;

            // gsWriteParaview(*g, name);
//            gsFileData<> fileData;
//            fileData << *g;
//            fileData.dump(name);

            multiPatch.addPatch(give(g)); // multi patch deletes the geometries
            //delete g;
        }

    }

#ifdef GISMO_WITH_PSOLID
    gsWriteParasolid(multiPatch, output + "_boundary");
#endif

    gsWriteGoTools(multiPatch, output + "_boundary");

    gsFileData<> fileData;
    fileData << multiPatch;
    fileData.dump(output + "_boundary");

    gsWriteParaview(multiPatch, output + "_boundary");

}



void checkTheJacobianDeterminant(
        const std::vector<gsTensorBSpline<3, real_t> >& volumes,
        const std::string& output)
{
    int numPositive = 0;

    for (size_t index = 0; index < volumes.size(); index++)
    {
        const gsTensorBSpline<3, real_t>& tbspl = volumes[index];

        gsMatrix<> para  = tbspl.support();

        gsVector<> c0 = para.col(0);
        gsVector<> c1 = para.col(1);
        gsMatrix<> pts = uniformPointGrid(c0, c1, 10000) ;

        bool minus = false;
        bool plus = false;
        bool zero = false;

        // here we will save points with positive Jacobian determinant
        // this are points we don't want
        gsMatrix<> mat(pts.rows(), pts.cols());
        mat.setZero();


        for (int col = 0; col < pts.cols(); col++)
        {
            real_t determinant = tbspl.jacobian(pts.col(col)).determinant();
            if (determinant < 0)
            {
                minus = true;
            }
            else if (determinant > 0)
            {
                gsMatrix<> tmp = tbspl.eval(pts.col(col));
                mat.col(col) = tmp.col(0);
                numPositive++;
                plus = true;
            }
            else
            {
                zero = true;
            }
        }

        if (plus)
        {
            std::string pointsString = output + "PositiveJacobianPoints" +
                                       util::to_string(index);

            gsWriteParaviewPoints(mat, pointsString);
        }

        gsInfo << "volume: " << index << "\n"
                  << "sign of Jac(tbspl): \n"
                  << "                    minus: " << minus << "\n"
                  << "                    plus: " << plus << "\n"
                  << "                    zero: " << zero << "\n"
                  << "                    number of sample points: "
                  << pts.cols() << "\n";
    }

    gsInfo << "Number of points where the Jacobian det is positive: "
              << numPositive << "\n";
}


void modifyVolumesForSiemens(const int param,
                             const int degree,
                             std::vector<gsTensorBSpline<3, real_t> >& volumes,
                             std::vector<std::pair<unsigned, boxSide> >& middleFaces
                             )
{

    gsVector<> axis1(2);
    axis1 << -0.2, 0.5;

    gsVector<> axis2(2);
    axis2 << 0.92, 0.5;

    if (param == 2)
    {
        gsTensorBSpline<3, real_t> bspl1 = constructVolume(volumes[1],
                boundary::north, 2, degree, axis1, 4);


        gsTensorBSpline<3, real_t> bspl2 = constructVolume(volumes[8],
                boundary::south, 2, degree, axis2, 2, false);

        gsTensorBSpline<3, real_t> bspl3 = constructVolume(volumes[6],
                boundary::front, 2, degree, axis2, 2);

        volumes.push_back(bspl1);
        volumes.push_back(bspl2);
        volumes.push_back(bspl3);

        middleFaces.push_back(std::make_pair(1, boundary::north));
        middleFaces.push_back(std::make_pair(6, boundary::front));
        middleFaces.push_back(std::make_pair(8, boundary::south));

        middleFaces.push_back(std::make_pair(9, boundary::back));

        middleFaces.push_back(std::make_pair(10, boundary::front));
        middleFaces.push_back(std::make_pair(10, boundary::south));
        middleFaces.push_back(std::make_pair(10, boundary::north));

        middleFaces.push_back(std::make_pair(11, boundary::back));
        middleFaces.push_back(std::make_pair(11, boundary::east));
        middleFaces.push_back(std::make_pair(11, boundary::west));
    }
    else if (param == 1)
    {

        gsTensorBSpline<3, real_t> bspl1 = constructVolume(volumes[3],
                boundary::back, 2, degree, axis1, 4, false);


        gsTensorBSpline<3, real_t> bspl2 = constructVolume(volumes[0],
                boundary::north, 2, degree, axis2, 2, true);

        gsTensorBSpline<3, real_t> bspl3 = constructVolume(volumes[5],
                boundary::south, 2, degree, axis2, 2, false);

        volumes.push_back(bspl1);
        volumes.push_back(bspl2);
        volumes.push_back(bspl3);

        middleFaces.push_back(std::make_pair(0, boundary::north));
        middleFaces.push_back(std::make_pair(3, boundary::back));
        middleFaces.push_back(std::make_pair(5, boundary::south));

        middleFaces.push_back(std::make_pair(7, boundary::front));

        middleFaces.push_back(std::make_pair(8, boundary::south));
        middleFaces.push_back(std::make_pair(8, boundary::north));
        middleFaces.push_back(std::make_pair(8, boundary::back));

        middleFaces.push_back(std::make_pair(9, boundary::south));
        middleFaces.push_back(std::make_pair(9, boundary::north));
        middleFaces.push_back(std::make_pair(9, boundary::front));
    }

}


void saveWholeChairStandExample(
        const std::vector<gsTensorBSpline<3, real_t> >& volumes,
        const std::string& output)
{
    for (int i = 1; i < 5; i++)
    {
        std::string outDirectory = output + "Part" +
                util::to_string(i) + "-";

        std::vector<gsTensorBSpline<3, real_t> > transformedVolumes;

        for (size_t index = 0; index < volumes.size(); index++)
        {
            const gsTensorBSplineBasis<3, real_t>& basis = volumes[index].basis();
            const gsMatrix<>& coefs = volumes[index].coefs();
            gsMatrix<> newCoefs;

            transformCoefs(coefs, static_cast<int>(i), newCoefs);

            gsTensorBSpline<3, real_t> newVolume(basis, newCoefs);

            transformedVolumes.push_back(newVolume);
        }

        saveVolumes(transformedVolumes, outDirectory);
    }
}

/// This function modifies volumes to get C0 continuity between pieces.
/// Input volumes are not watertight, function makes them watertight by
/// moving boundary control points.
///
/// \param[in] pairs vector of pairs [[Volume1 | Side1 | Volume2  | Side2], ...]
///                  Side1 of Volume1 must match Side2 of Volume2
/// \param[in] volBlock abstract presentation of volumes
/// \param[in, out] volumes this volumes will be modified
/// \param[out] middleFaces faces which match (Side1 of Volume1, Side2 of Volume2, ...)
void matchFaces(const std::vector<std::vector<unsigned> >& pairs,
                const std::vector<gsVolumeBlock<>* >& volBlocks,
                std::vector<gsTensorBSpline<3, real_t> >& volumes,
                std::vector<gsTensorBSpline<2, real_t> >& topFaces,
                std::vector<std::pair<unsigned, boxSide> >& middleFaces,
                bool randomizeGaps, real_t interfaceErrorSize,
                std::string interfaceErrors, std::string gapDim)
{
    gsInfo << "\n\nMatching:\n\n"
                 "Volume pairs that must match :\n"
              << "     | Vol1  | Side1 | Vol2  | Side2 " << "\n";

    for (size_t cPair = 0; cPair != pairs.size(); cPair++)
    {
        gsInfo << " " << cPair << ".  ";

        for (size_t index = 0; index != pairs[cPair].size(); index++)
        {
            gsInfo << "| " << pairs[cPair][index] << "     ";
        }
        gsInfo << "\n";
    }

    gsInfo << "\n";
    
    // A vector of all the "matches" between control points. Each call to
    // matchBoundaryOfVolumes will build up the matches in this vector.
    std::vector<CpMatchType> matches;

    for (size_t index = 0; index != pairs.size(); index++)
    {
        gsInfo << "Matching " << index << ". pair ...\n" << "\n";

        // volume block 1 and 2
        unsigned volBlockIndex1 = pairs[index][0];
        unsigned volBlockIndex2 = pairs[index][2];

        gsVolumeBlock<>* vb1 = volBlocks[volBlockIndex1];
        gsVolumeBlock<>* vb2 = volBlocks[volBlockIndex2];

        boxSide side1 = vb1->getSideOfHexahedron(pairs[index][1]);
        boxSide side2 = vb2->getSideOfHexahedron(pairs[index][3]);

        matchBoundaryOfVolumes(volBlockIndex1, volumes[volBlockIndex1], side1,
                               volBlockIndex2, volumes[volBlockIndex2], side2,
                               matches);

        middleFaces.push_back(std::make_pair(pairs[index][0], side1));
        middleFaces.push_back(std::make_pair(pairs[index][2], side2));
    }
    
    // perform the matches
    size_t nMatches = matches.size();
    for(size_t i = 0; i < nMatches; i++)
    {
        gsMatrix<> mid(1, 3);
        mid.setZero();
        
        // compute midpoint
        size_t nPts = matches[i].size();
        for(size_t j = 0; j < nPts; j++)
        {
            int thisCoefIndex = matches[i][j].ctrlPtId;
            mid += volumes[matches[i][j].patchId].coefs().row(thisCoefIndex);
        }
        mid /= nPts;
        
        // move all the relevant control points to the midpoint
        for(size_t j = 0; j < nPts; j++)
        {
            int thisCoefIndex = matches[i][j].ctrlPtId;
            volumes[matches[i][j].patchId].coefs().row(thisCoefIndex) = mid;
        }
    }

    // sanity check
    for(size_t i = 0; i < nMatches; i++)
    {
        gsMatrix<>& coefs1 = volumes[matches[i][0].patchId].coefs();
        int coefIndex1 = matches[i][0].ctrlPtId;
        gsMatrix<> testPoint = coefs1.row(coefIndex1);
        
        #ifndef NDEBUG
        size_t nPts = matches[i].size();
        for(size_t j = 0; j < nPts; j++)
        {
            GISMO_ASSERT(
            volumes[matches[i][j].patchId].coefs().row( matches[i][j].ctrlPtId ) == testPoint,
            "Failed sanity check, points did not match after they were forced to match");
        }
        #endif
    }
    
    gsInfo << "There are " << matches.size() << " matches.\n";


// now the edges and connection surfaces are watertight
// for testing, we need some examples, that are not watertight
// we need some plane surfaces and some volumes, with that property
// to get the surfaces, we extract the "top-level-faces" of the volumes
// this will only work, if the highest surfaces have a constant z-value
// afterwards, we move the corresponding controlpoints to get gaps or overlaps


// -------- extract top-level faces of chairStand model -----------------
    if(gapDim == "surface" || gapDim == "both")
    {
        // get top level-value (biggest z-component)
        real_t max = 0;
        real_t tol = 0.99;
        for(size_t i=0; i<volumes.size(); i++)
        {
            size_t NC = volumes[i].coefs().size()/3;
            for(size_t j=0; j<NC; j++)
            {
                if(volumes[i].coef(j)(2) > max)
                    max = volumes[i].coef(j)(2);
            }
        }

        // find volumes, that have faces on top-level
        // and extract them into surfaces
        std::vector<int> counts;
        std::vector<std::vector<gsVector<index_t, 3> > > indCoefs;
        for(size_t i=0; i<volumes.size(); i++)
        {
            int count = 0;
            int NC = volumes[i].coefs().size()/3;
            std::vector<gsVector<index_t, 3> > indices;
            for(int j=0; j<NC; j++)
            {
                if(volumes[i].coef(j)(2) > tol*max)
                {
                    count++;
                    gsVector<index_t, 3> index = volumes[i].basis().tensorIndex(j);
                    indices.push_back(index);
                }
            }
            counts.push_back(count);
            indCoefs.push_back(indices);
        }

        std::vector<std::pair<int,int> > constInd;

        //get CP-index of top face, depending on the orientation
        //of the volume it is either i, j, or k
        for(size_t i=0; i<volumes.size(); i++)
        {
            int NC = volumes[i].coefs().size()/3;
            int Nsquare = nthIntRoot(NC,3);
            Nsquare *= Nsquare;

            //if surface is not entirely on the top face
            if(counts[i] != Nsquare)
            {
                std::pair<int, int> p(-1,-1);
                constInd.push_back(p);
            }
            //compare the first and the last coef of the surface
            //two of the indices have to be different and one has to be the same
            else
            {
                if(indCoefs[i][0][0] == indCoefs[i][Nsquare-1][0])
                {
                    std::pair<int, int> p(0,indCoefs[i][0][0]);
                    constInd.push_back(p);
                }
                else if(indCoefs[i][0][1] == indCoefs[i][Nsquare-1][1])
                {
                    std::pair<int, int> p(1,indCoefs[i][0][1]);
                    constInd.push_back(p);
                }
                else if(indCoefs[i][0][2] == indCoefs[i][Nsquare-1][2])
                {
                    std::pair<int, int> p(2,indCoefs[i][0][2]);
                    constInd.push_back(p);
                }
            }
        }

        //construct surfaces and assign CPs of the surfaces
        std::vector<bool> onTop;
        for(size_t i=0; i<volumes.size(); i++)
        {
            //if the volume has a surface on the top-level
            if(constInd[i].first > -1)
            {
                int NC = volumes[i].coefs().size()/3;
                const int N = nthIntRoot(NC,3);
                const int Nsquare = N*N;

                gsKnotVector<> kv = volumes[i].basis().knots(0);
                gsMatrix<> coefs(Nsquare,3);
                gsTensorBSpline<2,real_t> face(kv, kv, coefs);

                for(int r=0; r<N; r++)
                {
                    for(int s=0; s<N; s++)
                    {
                        int indVol;
                        if(constInd[i].first == 0)
                            indVol = volumes[i].basis().index(constInd[i].second,r,s);
                        else if(constInd[i].first == 1)
                            indVol = volumes[i].basis().index(r,constInd[i].second,s);
                        else // if(constInd[i].first == 2)
                            indVol = volumes[i].basis().index(r,s,constInd[i].second);

                        const int indSur = face.basis().index(r,s);
                        face.coef(indSur) = volumes[i].coef(indVol);
                    }
                }
                topFaces.push_back(face);
                onTop.push_back(true);
            }
            //if the volume has no surface on the top-level, create a dummy-surface
            //that way, the indices of the surfaces and the volumes will be the same
            else
            {
                topFaces.push_back( gsTensorBSpline<2,real_t>() );
                onTop.push_back(false);
            }
        }


    // ------------- get matching CPs of TopLevel-surfaces ------------------
    //we need the CPs, that lie on more than one surface

        // A vector of all the "matches" between control points of the top-level surfaces.
        std::vector<CpMatchType> surMatches;

        for(size_t i=0; i<nMatches; i++)
        {
            int volID = matches[i][0].patchId;
            int ptID = matches[i][0].ctrlPtId;

            // if CP is on topFace
            if(volumes[volID].coef(ptID)(2) > tol*max)
            {
                CpMatchType t;
                int nCPs = 0;
                for(size_t j=0; j<matches[i].size(); j++)
                {
                    int patchID = matches[i][j].patchId;
                    int pointID = matches[i][j].ctrlPtId;

                    nCPs = topFaces[patchID].coefsSize();
                    //if surface is no dummy-surface, that is, if the volume has a top-level surface
                    if(nCPs>0)
                    {
                        int surID = -1;
                        int k = 0;
                        // get CP-index with respect to surface
                        while(k<nCPs && surID == -1)
                        {
                            if(topFaces[patchID].coef(k) == volumes[patchID].coef(pointID))
                                surID = k;
                            k++;
                        }
                        CpMatchEntry m;
                        m.patchId = patchID;
                        m.ctrlPtId = surID;

                        //get boxSide of CP, 0 if on a corner
                        int bs = 0;
                        const int n = nthIntRoot(nCPs,2);

                        if(surID > 0 && surID < n-1)
                            bs = 4;
                        else if(surID > n*(n-1) && surID < n*n-1)
                            bs = 3;
                        else
                        {
                            int l=1;
                            while(l<n-1 && bs == 0)
                            {
                                if(surID == l*n)
                                    bs = 1;
                                l++;
                            }
                            l=2;
                            while(l<n && bs == 0)
                            {
                                if(surID == l*n-1)
                                    bs = 2;
                                l++;
                            }
                        }
                        m.bs = bs;
                        t.push_back(m);
                    }
                }
                if(nCPs > 0)
                    surMatches.push_back(t);
            }
        }

    // ------------------- create gaps, overlaps or both for surfaces --------------------
        srand(static_cast<unsigned>( time(NULL)) );
        if(interfaceErrorSize > 0)
        {
            for(size_t i=0; i<surMatches.size(); i++)
            {
                size_t nPts = surMatches[i].size();
                for(size_t j=0; j<nPts; j++)
                {
                    boxSide side = surMatches[i][j].bs;
                    if(!(side == 0)) //if CP does not lie on a corner
                    {
                        int patchId = surMatches[i][j].patchId;
                        int ctrlPtId = surMatches[i][j].ctrlPtId;

                        //get inner neighbour ID
                        gsVector<index_t,2> ind = topFaces[patchId].basis().tensorIndex(ctrlPtId);
                        gsVector<index_t,2> indNB = ind;

                        //depending on which side we are, the inner neighbour will be computed
                        switch(side)
                        {
                        case 1: indNB[0]++;
                            break;
                        case 2: indNB[0]--;
                            break;
                        case 3: indNB[1]--;
                            break;
                        case 4: indNB[1]++;
                            break;
                        default:
                            GISMO_ERROR("unknown side "<<side);
                            break;
                        }
                        int innerNeighbourId = topFaces[patchId].basis().index(indNB[0], indNB[1]);

                        //get vector between point and inner neighbour
                        gsVector<real_t,3> point = topFaces[patchId].coef(ctrlPtId);
                        gsVector<real_t,3> neighbour = topFaces[patchId].coef(innerNeighbourId);

                        gsVector<real_t,3> dir = neighbour - point;
                        real_t norm = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
                        dir = dir / norm;

                        //move point towards inner neighbour
                        gsVector<real_t,3> movingVec = dir;
                        movingVec *= interfaceErrorSize;
                        if(randomizeGaps == true)
                        {
                            real_t randNum = rand() % 1000;
                            randNum /= 1000;
                            movingVec *= randNum;
                        }
                        int gapOverlap;
                        if(interfaceErrors == "gaps")
                            gapOverlap = 1;
                        else if(interfaceErrors == "overlaps")
                            gapOverlap = -1;
                        else
                        {
                            int randNum = rand() % 2;
                            gapOverlap = -1 + 2*randNum;
                        }

                        movingVec *= gapOverlap;

                        topFaces[patchId].coef(ctrlPtId) = point + movingVec;
                    }
                }
            }
        }
    }

// -------- create gaps, overlaps or both for volumes ------------------------------------------------

    if(gapDim == "volume" || gapDim == "both")
    {
        if(interfaceErrorSize > 0)
        {
            //srand(time(NULL));
            for(size_t i = 0; i < nMatches; i++)
            {
                size_t nPts = matches[i].size();
                if(nPts == 2) //if there are more than two controlpoints, they lie on an edge
                {
                    for(size_t j = 0; j < nPts; j++)
                    {
                        boxSide side = matches[i][j].bs;
                        if(!(side == 0)) //if controlpoint does not lie on an edge
                        {
                            //get inner neighbour ID
                            int patchId = matches[i][j].patchId;
                            int ctrlPtId = matches[i][j].ctrlPtId;
                            gsVector<index_t,3> ind = volumes[patchId].basis().tensorIndex(ctrlPtId);
                            gsVector<index_t,3> indNB = ind;

                            //depending on which side we are, the inner neighbour will be computed
                            switch(side)
                            {
                            case 1: indNB[0] ++;
                                break;
                            case 2: indNB[0] --;
                                break;
                            case 3: indNB[1] ++;
                                break;
                            case 4: indNB[1] --;
                                break;
                            case 5: indNB[2] ++;
                                break;
                            case 6: indNB[2] --;
                                break;
                            default:
                                GISMO_ERROR("unknown side "<<side);
                                break;
                            }

                            int innerNeighbourId = volumes[patchId].basis().index(indNB[0], indNB[1], indNB[2]);

                            //get vector between point and inner neighbour
                            gsVector<real_t,3> point = volumes[patchId].coef(ctrlPtId);
                            gsVector<real_t,3> neighbour = volumes[patchId].coef(innerNeighbourId);

                            gsVector<real_t,3> dir = neighbour - point;
                            real_t norm = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
                            dir = dir / norm;



                            //move point towards inner neighbour
                            gsVector<real_t,3> movingVec = dir;
                            movingVec *= interfaceErrorSize;
                            if(randomizeGaps == true)
                            {
                                real_t randNum = rand() % 1000;
                                randNum /= 1000;
                                movingVec *= randNum;
                            }
                            int gapOverlap;
                            if(interfaceErrors == "gaps")
                                gapOverlap = 1;
                            else if(interfaceErrors == "overlaps")
                                gapOverlap = -1;
                            else
                            {
                                int randNum = rand() % 2;
                                gapOverlap = -1 + 2*randNum;
                            }

                            movingVec *= gapOverlap;

                            volumes[patchId].coef(ctrlPtId) = point + movingVec;
                        }
                    }
                }
            }
        }
    }
}


void orderUp(gsTensorBSpline<3, real_t> &volume, bool expectNoChange)
{
    // change the coordinate system so that the third coordinate (in the
    // domain) points up.
    gsInfo << "Changing coordinates so that the third coordinate points approximately up.\n";
    
    gsMatrix<real_t> lookDirs(6, 3);
    lookDirs << 1, 0, 0,
            0, 1, 0,
            0, 0, 1,
            -1, 0, 0,
            0, -1, 0,
            0, 0, -1;
    
    // compute the Jacobian at a bunch of sample points and find the mean
    
    gsMatrix<real_t> para  = volume.support();
    gsVector<real_t> c0 = para.col(0);
    gsVector<real_t> c1 = para.col(1);
    
    gsMatrix<real_t> pts = uniformPointGrid(c0, c1, 125);
    gsMatrix<real_t> meanJ(3, 3);
    meanJ.setZero();
    
    for(int col = 0; col < pts.cols(); col++)
    {
        gsMatrix<real_t> j = volume.jacobian(pts.col(col));
        meanJ += j;
    }
    meanJ /= pts.cols();
    
    // multiply mean Jacobian by each of the lookDirs, and
    // choose the lookDir which results in a vector close to up
    gsMatrix<real_t> lookVals = meanJ * lookDirs.transpose();
    int mi = -1;
    real_t md = 0;
    gsMatrix<real_t> up(1, 3);
    up << 0, 0, 1;
    for(int col = 0; col < lookVals.cols(); col++)
    {
        gsMatrix<real_t> nd = up * lookVals.col(col);
        GISMO_ASSERT(nd.rows() == 1 && nd.cols() == 1, "Unexpected matrix size");
        if(nd(0, 0) > md)
        {
            mi = col;
            md = nd(0, 0);
        }
    }
    
    if(expectNoChange)
    {
        if(mi != 2)
        {
            gsInfo << "Sanity check failed: expected not to have to rearrange control points\n";
            GISMO_ASSERT(false, "Sanity check failed");
        }
        return;
    }
    
    // based on the best lookDir, reorder the control points of the volume
    
    gsTensorBSplineBasis<3, real_t> b = volume.basis();
    
    std::vector<int> sizes;
    sizes.push_back(b.size(0));
    sizes.push_back(b.size(1));
    sizes.push_back(b.size(2));
    
    GISMO_ASSERT(sizes[0] == sizes[1] && sizes[1] == sizes[2],
                 "Only implemented for volumes with the same # ctrl points on each dimension");
    
    std::vector<int> oldStrides;
    oldStrides.push_back(b.stride(0));
    oldStrides.push_back(b.stride(1));
    oldStrides.push_back(b.stride(2));
    
    gsMatrix<real_t> newCoefs(oldStrides[2] * sizes[2], 3);
    
    std::vector<int> newStrides;
    int newOffset = 0;
    int rot = mi % 3;
    bool refl = (mi >= 3);
    
    // rotate so that dimension rot maps to dimension 2
    newStrides.push_back(oldStrides[(2*rot + 2) % 3]);
    newStrides.push_back(oldStrides[(2*rot + 0) % 3]);
    newStrides.push_back(oldStrides[(2*rot + 1) % 3]);
    
    if(refl)
    {
        // flip two of the dimensions so that the third dimension maps up
        newOffset += newStrides[rot] * (sizes[rot] - 1);
        newStrides[rot] = -newStrides[rot];
        
        int rotNext = (rot + 1) % 3;
        newOffset += newStrides[rotNext] * (sizes[rotNext] - 1);
        newStrides[rotNext] = -newStrides[rotNext];
    }
    
    // now do the rearrangment of indices
    for(int i = 0; i < sizes[0]; i++)
        for(int j = 0; j < sizes[1]; j++)
            for(int k = 0; k < sizes[2]; k++)
            {
                int inpIdx = i * oldStrides[0] + j * oldStrides[1] + k * oldStrides[2];
                int outpIdx = i * newStrides[0] + j * newStrides[1] + k * newStrides[2] + newOffset;
                newCoefs.row(outpIdx) = volume.coefs().row(inpIdx);
            }
    
    gsMatrix<real_t> oldCoefs = volume.coefs();
    gsMatrix<real_t> coeffDiff = newCoefs - oldCoefs;
    volume.setCoefs(newCoefs);
    
    orderUp(volume, true);
}



int main(int argc, char* argv[])
{
    index_t n = 10;
    real_t lambda = 1e-9;
    //real_t eps = 0.0;
    index_t degree = 3;
    index_t internalKnots = 3;
    bool full = false;
    bool alignVertical = false;
    std::string pairFile( "solids/chairStand/"
                          "segmented_objStandFiveParts_pairs.txt");
    std::string directory("solids/chairStand/");
    std::string output("chairStand");
    bool randomizeGaps = false;
    std::string interfaceErrors("gaps");
    real_t interfaceErrorSize = 0.005;
    std::string gapDim("both");


    // "/michael/";
    // "/michael/facePairs.txt"
    
   // gsCmdLine cmd("Reading gsSolid");
    
    gsCmdLine cmd("Reading gsSolid");
    cmd.addSwitch("alignVertical", 
                  "Rearrange control points so that increasing z-parameter points "
                  "approximately upwards", alignVertical);
    cmd.addInt("n", "numPts", "Number of sample points", n);
    cmd.addReal("l", "lambda", "Smoothing parameter", lambda);
    cmd.addInt("d", "degree", "Degree of the B-Splines", degree);
    cmd.addInt("k", "numKnots", "Number of interior knots", internalKnots);
    cmd.addString("p", "pairFile", "Pair File", pairFile);
    cmd.addString("m", "directory",
                  "Directory where files are stored", directory);
    cmd.addString("o", "output", "Output prefix", output);
    cmd.addSwitch("full", "Outputs the hole chair stand.", full);
    cmd.addSwitch("randomize", "Randomize gaps and overlaps", randomizeGaps);
    cmd.addString("g", "gapDim", "extract top-level surfaces, "
                                 "and create 2D examples of gaps and overlaps "
                                 "(options are: volumes, surfaces, both)", gapDim);
    cmd.addString("i", "interfaceErrors",
                  "Introduce artificial interface errors "
                  "(options are: gaps, overlaps, mixed)",
                  interfaceErrors);
    cmd.addReal("s", "interfaceErrorSize",
                                       "Parameter that determines the size of "
                                       "the gap or overlap", interfaceErrorSize);

    //gsArgVal<> epsArg("e", "eps", "Precision for computing lengths of curves",
    //                           false, 1e-6, "double", cmd);
    //eps = epsArg.getValue();
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "------------------------------------------------------------"
                 "\nInput Arguments: \n\n"
                 "Ouput prefix: " << output << "\n\n"
                 "Number of sample points: " << n << "\n\n"
                 "Lambda: " << lambda << "\n\n"
                 "Degree: " << degree << "\n\n"
                 "Number of internal knots: " << internalKnots << "\n\n"
                 "Full chair stand example (only for chair stand example): " <<
                 full << "\n\n"
                 "Input directory: " << directory << "\n\n"
                 "Input file (how pieces match between each other): " <<
                 pairFile << "\n\n"
                 "Watertightnes-mode: "<< interfaceErrors << "\n\n"
                 "Surfaces or Volumes:" << gapDim << "\n\n"
                 "Randomized distances: " << randomizeGaps << "\n\n"
                 "Maximum distance: " << interfaceErrorSize << "\n\n"
                 "------------------------------------------------------------"
                 "\n\n";

    if(interfaceErrors != "gaps" && interfaceErrors != "overlaps" && interfaceErrors != "mixed")
    {
        std::cerr << "Error: wrong parameter for interfaceErrors" <<"\n";
        return -1;
    }
    if(interfaceErrorSize < 0)
    {
        std::cerr << "Error: negative parameter for interfaceErrorSize" <<"\n";
        return -1;
    }


    std::vector<std::string> allFiles = getAllFileNames(pairFile);


    for (size_t index = 0; index < allFiles.size(); index++)
    {
        //gsInfo << "File " << index << ": " << allFiles[index] << "\n";
        std::string f = allFiles[index];
        f = f.substr(0, f.size() - 4);
        gsInfo << index << " " << f << "\n";
    }

    // -------------------------------------------------------------------------
    // containers
    // -------------------------------------------------------------------------
    // file names of the individual pieces
    std::vector<std::string> files;

    // pointers to the boundary representation (trimmed surfaces) of the pieces
    std::vector<gsVolumeBlock<>* > volBlocks;

    // volumes of the pieces
    std::vector<gsTensorBSpline<3, real_t> > volumes;

    // faces incident to many pieces (they needs to be removed if we plot the
    // boundary / shell)
    std::vector<std::pair<unsigned, boxSide> > middleFaces;

    // topFaces of the pieces
    std::vector<gsTensorBSpline<2, real_t> > topFaces;

    // pointers to the boundary representation
    std::vector<gsSolid<>*> solids;
    // -------------------------------------------------------------------------


    // -------------------------------------------------------------------------
    // CREATING VOLUMES FROM BOUNDARY REPRESENTATION
    // -------------------------------------------------------------------------


    gsInfo << "\n\n\nFitting:\n\n";
    
    for (size_t i = 0; i < allFiles.size(); ++i)
    {
        std::string fileName(allFiles[i]);
        files.push_back(fileName);

        gsInfo << "(" << i << ") File: " << fileName << "\n"
                     "Volume fitting..." << "\n";

        gsSolid<>::uPtr pSolid;

        gsFileData<> data(directory + fileName);
        if (data.has< gsSolid<> >())
        {
            pSolid = data.getFirst< gsSolid<> >();
        }
        if(!pSolid)
        {
            std::cerr << "Either there was an error loading the file or no gsSolid was found.\n";
            return -1;
        }

        // input solids have just one volumeBlock
        volBlocks.push_back(pSolid->volume[0]);

        gsTensorBSpline<3, real_t> volume =
                volBlocks[i]->getTensorBSplineVolume(n, internalKnots, degree, lambda, 1e-3);
        

        volumes.push_back(volume);
        solids.push_back(pSolid.release());

        gsInfo << "\n";
    }


    // -------------------------------------------------------------------------
    // MAKING INTERFACES WATERTIGHT, MOVING CONTROL POINTS TO GET C0 CONTINUITY
    // BETWEEN PIECES
    // -------------------------------------------------------------------------


    std::vector<std::vector<unsigned> > pairs;
    getIncidentFaces(pairs, files, pairFile);
    
    matchFaces(pairs, volBlocks, volumes, topFaces, middleFaces, randomizeGaps,
               interfaceErrorSize, interfaceErrors, gapDim);
    
    
    // -------------------------------------------------------------------------
    // REORDER CONTROL POINTS SO THAT INCREASING Z-PARAMETER POINTS
    // APPROXIMATELY UPWARDS
    // -------------------------------------------------------------------------
    
    if(alignVertical)
    {
        for(size_t i = 0; i < volumes.size(); i++)
        {
            orderUp(volumes[i], false);
            orderUp(volumes[i], true); // test that it worked
        }
    }


    // -------------------------------------------------------------------------
    // SAVING VOLUMES, CHECKING JACOBIAN DETERMINANT, SIEMENS RELEVANT STUFF
    // -------------------------------------------------------------------------


    // siemens wants some additional volumes
    // modifyVolumesForSiemens(-1, degree, volumes, middleFaces);


    // uncomment if you want to check the sign of the Jacobian determinant
    // in 10000 sample points
    checkTheJacobianDeterminant(volumes, output);


    // saving

    gsInfo << "Saving volumes with prefix " << output << "\n";
    saveVolumes(volumes, output, true);

    gsInfo << "Saving topFaces with prefix " << output << "\n";
    saveFaces(topFaces, output);

    // uncomment if you want to save boundary (shell) of the volumes
    // saveBoundary(volumes, middleFaces, output);

    // relavant just for chair stand example
    // we can plot the whole example, not just one part
    if (full)
    {
        saveWholeChairStandExample(volumes, output);
    }


    for (size_t index = 0; index < solids.size(); index++)
    {
        delete solids[index];
    }

    return 0;
}
