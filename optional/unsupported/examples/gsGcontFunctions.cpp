/** @file gsGcontFunctions.cpp

    @brief Searches for G1 functions on a multipatch domain

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Birner, A. Mantzaflaris
*/

#include <gismo.h>
#include <gismo_dev.h>

// #include <boost/functional/hash.hpp>

using namespace gismo;

#define EigenSLV Eigen::FullPivLU<gsMatrix<>::Base>

#define useSparse 0
#if useSparse
#define MatType   gsSparseMatrix<real_t,ColMajor>
//#define MatType   gsSparseMatrix<real_t,RowMajor>
#else
#define MatType   gsMatrix<>
#endif


gsMultiBasis<> mb;

std::vector<index_t> cind, rind;
inline gsMatrix<> columns(const MatType & mat,
                          const gsVector<index_t> & ind)
{
#           if useSparse
            cind.assign(ind.begin(),ind.end());
            rind = mat.innerOf(cind);
            return mat.submatrix(rind,cind);
#           else
            gsMatrix<> sm;
            mat.submatrixCols(ind, sm);
            return sm;
#           endif
}

// Computes the matrix
void G1Matrix(const gsMultiPatch<> & mp,
              const index_t pdeg,
              MatType & mat, gsDofMapper & map, bool randomize = false);

void G1Matrix(const gsMultiPatch<> & mp,
              const index_t pdeg,
              const std::vector<gsGeometry<>::Ptr> & pert1,
              MatType & mat,
              gsDofMapper & map); //, bool randomize);


class gsGluingData;

void G1Matrix(const gsMultiPatch<> & mp,
              const index_t pdeg,
              const gsGluingData & gd,
              MatType & mat,
              gsDofMapper & map);

// Full geometric basis with non-homogeneous boundary
gsSparseMatrix<> FullBasisNonHom(const gsMultiPatch<> & mp,
                            // const index_t pdeg,
                            const gsGluingData & gd, gsDofMapper & map);


void PlotFunctions(const gsMatrix<> & ker, gsMultiPatch<> & mp,
                   const gsDofMapper & map, const int & ns, const bool & pm, const bool & pn);

gsVector<index_t> countSparsity(const gsMatrix<> & K);

gsMatrix<> eval_R(const gsKnotVector<> & m_knots,
                       const gsGeometry<> & patch,
                       const gsMatrix<> & u);

void localBasis_elwise(const gsMultiPatch<> & mp,
                       const gsMatrix<> & M,
                       const gsDofMapper & map,
                       gsMatrix<> & result);

template<int d>
void localBasis_incr(const gsMultiPatch<> & mp,
                     const index_t pdeg,
                     const std::vector<gsGeometry<>::Ptr> & pert,
                     gsDofMapper & map,
                     gsMatrix<> & result, bool bdr = false);

template<int d>
void localBasis_incr(const gsMultiPatch<> & mp,
                     const index_t pdeg,
                     gsDofMapper & map,
                     gsMatrix<> & result, bool bdr = false);

std::pair<index_t,index_t> computeDim(const gsMultiPatch<> & mp,  bool randomize = false);

size_t computeDimFormula(const gsMultiPatch<> & mp,  bool randomize = false);

void printSparsity(const gsMatrix<> & K, const gsDofMapper & map, const gsMultiPatch<> & mp);

std::string printSparsity(const gsMatrix<> & K);

template<int d> void
plotSparsity(const gsMatrix<> & K, const gsDofMapper & map, gsMultiPatch<> & mp);

void sortBySparsity(gsMatrix<> & K);

void selectRepresentatives(gsMatrix<> & K, gsVector<index_t> & count);

// true if the matrix is in RREF
bool isRref(const gsMatrix<> & M)
{
    const index_t nr = M.rows();
    const index_t nc = M.cols();
    index_t i, piv = 0;
    for (index_t r=0; r!=nr; ++r)
    {
        if (nc <= piv) return true;
        i = r;
        while (0 == M(i, piv))
        {
            ++i;
            if (nr == i)
            {
                i = r;
                ++piv;
                if (nc == piv) return true;
            }
        }

        if ( 1 != M(r, piv) ) return false;
    }
    return true;
}

// true if the matrix is in RREF
bool isRref(const gsSparseMatrix<real_t,RowMajor> & M)
{
    const index_t nr = M.rows();
    const index_t nc = M.cols();
    index_t i, piv = 0;
    for (index_t r=0; r!=nr; ++r)
    {
        if (nc <= piv) return true;
        i = r;
        while (0 == M(i, piv))
        {
            ++i;
            if (nr == i)
            {
                i = r;
                ++piv;
                if (nc == piv) return true;
            }
        }

        if ( 1 != M(r, piv) ) return false;
    }
    return true;
}

// Compute the nullspace dimension for a matrix in RREF
int rrefKerDim(const gsMatrix<> & M)
{
    GISMO_ASSERT( isRref(M), "Not in RREF form" );

    const index_t nr = M.rows();
    const index_t nc = M.cols();

    index_t k = 0, c = 0;
    for (index_t r=0; r<nr; ++r)
    {
        while ( c<nc && 0 == M(r,c) ) { ++c; ++k;}
        ++c;
    }
    if (c < nc ) k += nc - c;
    return k;
}

// Compute the nullspace dimension for a matrix in RREF
template<int rr>
int rrefKerDim(const gsSparseMatrix<real_t, rr> & M)
{
    //GISMO_ASSERT( isRref(M), "Not in RREF form" );

    const index_t nr = M.rows();
    const index_t nc = M.cols();

    index_t k = 0, c = 0;
    for (index_t r=0; r<nr; ++r)
    {
        while ( c<nc && 0 == M(r,c) ) { ++c; ++k;}
        ++c;
    }
    if (c < nc ) k += nc - c;
    return k;
}

// Compute a nullspace basis for a matrix in RREF
index_t rrefKerBasis(const gsMatrix<> & M, gsMatrix<> & ker,
                     std::vector<index_t> & nonPiv)
{
    GISMO_ASSERT( isRref(M), "Not in RREF form" );
    nonPiv.clear();

    std::vector<std::pair<index_t,index_t> > piv;
    piv.clear(); // (row,col)
    piv.reserve(std::max<index_t>(M.cols()-M.rows()+1,2));
    index_t i = 0, j = 0, k = 0;
    while(i<M.rows() && j<M.cols())
    {
        if ( 1 == M(i,j) )
        {
            piv.push_back( std::make_pair(i,j));
            ++i;
        }
        ++j;
    }

    if ( (size_t)M.cols() == piv.size() )
    {
        ker.setZero(M.cols(),1);
        return 0;
    }

    j = 0;// number of pivots at the left of col(i)
    ker.resize(M.cols(), M.cols()-piv.size());
    nonPiv.reserve(ker.cols());
    piv.push_back(std::make_pair(0,M.cols()));// clamp pivots
    for (i = 0; i!= M.cols(); ++i)
    {
        if (i==piv[j].second)
        {
            ++j;
            continue;
        }

        nonPiv.push_back(i);
        ker.col(k).setZero();
        ker(i,k) = -1;
        for (index_t l = 0; l!=j; ++l)
            ker(piv[l].second, k) = M(piv[l].first, i);
        ++k;
    }

    return ker.cols();
}

// Compute a nullspace basis for a matrix in RREF
index_t rrefKerBasis(const gsMatrix<> & M, gsMatrix<> & ker)
{
    GISMO_ASSERT( isRref(M), "Not in RREF form" );
    std::vector<index_t> npiv;
    return rrefKerBasis(M, ker, npiv);
}

bool isC1Boundary(const gsVector<index_t,3> & t, index_t lF)
{
    if( (t.head(2).array() > 1).all() &&
        (t.head(2).array() < lF).all() )
    {
        return false;
    }
    return true;
}

bool isBoundary(const gsVector<index_t,3> & t, index_t lF)
{
    if( (t.head(2).array() >= 1).all() &&
        (t.head(2).array() <= lF).all() )
    {
        return false;
    }
    return true;
}

bool isVertex(const gsVector<index_t,3> & t,  index_t lF) //Not all vertices yet!!!
{
    if( (t.head(2).array() < 2).all()  ||
        (t.head(2).array() >= lF).all() )
    {
        return true;
    }
    if( (t.head(2).array() < 2).any()  &&
        (t.head(2).array() >= lF).any() )
    {
        return true;
    }
    return false;

}

void perturbPatch(gsGeometry<> & mp, const real_t magnitude)
{
    gsMatrix<> randM;
//    std::srand((unsigned int) time(0));
    gsMatrix<> & cf = mp.coefs();
    randM = 0.5 * magnitude * gsMatrix<>::Random(cf.rows(), cf.cols());
    cf += randM;
}

gsGeometry<>::uPtr randomPoly3(int deg1, int deg2, const real_t magnitude)
{
    gsKnotVector<> KV0(0,1,0,deg1+1);
    gsKnotVector<> KV1(0,1,0,deg2+1);
    //gsKnotVector<> KV2(0,1,0,2); // degree = 1
    gsTensorBSplineBasis<2> b(KV0, KV1); //, KV2
//    std::srand((unsigned int) time(0));
    gsMatrix<> randM = 0.5 * magnitude * gsMatrix<>::Random(b.size(), 1);
    return b.makeGeometry(give(randM));
}

gsGeometry<>::uPtr randomPoly2(int deg1, const real_t magnitude)
{
    gsKnotVector<> KV0(0,1,0,deg1+1);
    //gsKnotVector<> KV1(0,1,0,2); // degree = 1
    gsBSplineBasis<> b(KV0); //, KV1
//    std::srand((unsigned int) time(0));
    gsMatrix<> randM = 0.5 * magnitude * gsMatrix<>::Random(b.size(), 1);
    return b.makeGeometry(give(randM));
}

void getPerturbation(std::vector<gsGeometry<>::Ptr> & pert1, const index_t d, const real_t tol = 10e-1)
{
    pert1.resize(d+1, gsGeometry<>::Ptr());
    if (3==d)
    {
        /*
        pert1[0] = randomPoly3(d  ,d-1, tol);
        pert1[1] = randomPoly3(d-1,d  , tol);
        pert1[2] = randomPoly3(d-1,d-1, tol);
        pert1[3] = randomPoly3(d-1,d-1, tol);
        */
        ///*
        pert1[0] = randomPoly3(d-2,d-3, tol);
        pert1[1] = randomPoly3(d-3,d-2, tol);
        pert1[2] = randomPoly3(d-3,d-3, tol);
        pert1[3] = randomPoly3(d-3,d-3, tol);
        //*/

        /* // Store the perturbation for future loading
          gsFileData<> fd;
          fd<< *pert1[0];
          fd<< *pert1[1];
          fd<< *pert1[2];
          fd<< *pert1[3];
          fd.save("ggdata");
        */
    }
    else
    {
        pert1[0] = randomPoly2(d  , tol);
        pert1[1] = randomPoly2(d-1, tol);
        pert1[2] = randomPoly2(d-1, tol);
    }
}

template<short_t d>
typename gsTensorBSplineBasis<d-1>::Self_t getGluingBasis(int deg1, int deg2 = -1)
{
    GISMO_ASSERT(deg2>=0, "Error");
    gsKnotVector<> KV0(0,1,0,deg1+1);
    gsKnotVector<> KV1(0,1,0,deg2+1);
    return gsTensorBSplineBasis<2>(KV0, KV1);
}

template<>
gsBSplineBasis<> getGluingBasis<2>(int deg, int deg2)
{
    GISMO_UNUSED(deg2);
    GISMO_ASSERT(-1==deg2, "Error");
    gsKnotVector<> KV0(0,1,0,deg+1);
    return gsBSplineBasis<>(KV0);
}

void perturbGeometry(gsMultiPatch<> & mp, const real_t magnitude, bool commonBdr)
{
    gsInfo<<"Perturb input\n";

    const index_t np = mp.size();
    const boundaryInterface & iFace = *mp.iBegin();

    gsMatrix<> randM;
//    std::srand((unsigned int) time(0));
    gsMatrix<index_t> bd;
    gsVector<index_t> p(2);
    p(0) = iFace.first() .patch;
    p(1) = iFace.second().patch;

    for (index_t i = 0; i!=np; ++i)
    {
        gsMatrix<> & cf = mp.patch(p(i)).coefs();
        randM = 0.5 * magnitude * gsMatrix<>::Random(cf.rows(), cf.cols());

        if (!commonBdr)
        {
            bd = mp.patch(i).basis().boundary(
                          p(i) == 0 ? iFace.first().side() : iFace.second().side());
            for ( index_t k = 0; k != bd.rows(); ++k)
                randM.row(bd.at(k)).setZero();
        }

        cf += randM;
    }

    if (commonBdr)
        mp.closeGaps(magnitude);
}

void genericGeometryDim(const gsMultiPatch<> & mp, const real_t magnitude, bool commonBdr, index_t n = 10)
{
    gsInfo<<"--------- Initial ---------\n";
    computeDimFormula(mp, false);

    gsMultiPatch<> mpp;
    typedef std::map<size_t,size_t> map_t;
    map_t counts;

    for (index_t i = 0; i!=n; ++i)
    {
        gsInfo<<"--------- Trial "<<i+1<<" ---------\n";
        mpp = mp;
        perturbGeometry(mpp, magnitude, commonBdr);
        const size_t cc = computeDimFormula(mp, false);
        ++counts[cc];
    }

    gsInfo<< "Statistics:\n";
    n = 0;
    for(map_t::const_iterator it = counts.begin(); it != counts.end(); ++it,++n)
        std::cout <<"Key "<< n << " found " << it->second << " times.\n";
}

void genericGluingDim(const gsMultiPatch<> & mp, const real_t magnitude, index_t n = 10)
{
    /*
    gsInfo<<"--------- Initial ---------\n";
    computeDimFormula(mp, false);
    //*/

    GISMO_UNUSED(magnitude); // (!) tol is hardcoded in G1matrix, for now

    typedef std::map<size_t,size_t> map_t;
    map_t counts;

    for (index_t i = 0; i!=n; ++i)
    {
        gsInfo<<"--------- Trial "<<i+1<<" ---------\n";
        const size_t cc = computeDimFormula(mp, true);
        ++counts[cc];
    }

    gsInfo<< "Statistics:\n";
    n = 0;
    for(map_t::const_iterator it = counts.begin(); it != counts.end(); ++it,++n)
        std::cout <<"Key "<< n << " found " << it->second << " times.\n";
}


/// Define the source function  -exp(x+z).*sin(y)
class gsGluingData
{
public:
    gsGluingData(const gsMultiPatch<> & mp,
                 // const index_t pdeg,
                 const std::vector<gsGeometry<>::Ptr> & pert1 = std::vector<gsGeometry<>::Ptr>())
    : d(mp.parDim())
    {
        if ( 3 == d )
            init<3>(mp, /*pdeg,*/ pert1);
        else
            init<2>(mp, /*pdeg,*/ pert1);
    }

    gsGluingData(gsMatrix<index_t> & degrees)
    {
        // degrees has dimension (d-1)x(d+1)
        GISMO_UNUSED(degrees);
        GISMO_NO_IMPLEMENTATION;
    }


    template<short_t _d> void
    init(const gsMultiPatch<> & mp,
         // const index_t pdeg,
         const std::vector<gsGeometry<>::Ptr> & pert1 )
    {
        const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
        const gsGeometry<> & P1 = mp.patch( iFace.first().patch );
        const gsGeometry<> & P2 = mp.patch( iFace.second().patch );
        const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));

        const index_t l = iFace.first()  .direction(); // "orthogonal" direction p1
        const index_t l2 = iFace.second().direction(); // "orthogonal" direction p2

        gsMatrix<> u2, ev1, ev2, glpt, gleval, tmp, D0(d,d);


        //(3,2), (2.3), (2,2), (2,2)
        std::vector<typename gsTensorBSplineBasis<_d-1>::Self_t > gb(d+1);

        if ( 3 == _d )
        {
            /*
            gb[0] = getGluingBasis<_d>(d  ,d-1);
            gb[1] = getGluingBasis<_d>(d-1,d  );
            gb[2] = getGluingBasis<_d>(d-1,d-1);
            gb[3] = getGluingBasis<_d>(d-1,d-1);
            */
            ///*
            gb[0] = getGluingBasis<_d>(d-2,d-3);
            gb[1] = getGluingBasis<_d>(d-3,d-2);
            gb[2] = getGluingBasis<_d>(d-3,d-3);
            gb[3] = getGluingBasis<_d>(d-3,d-3);
            //*/
        }
        else // 2 == _d
        {
            gb[0] = getGluingBasis<_d>(d  );
            gb[1] = getGluingBasis<_d>(d-1);
            gb[2] = getGluingBasis<_d>(d-1);
        }

        gdata.resize(d+1);

        int sgn = 1;
        gsMatrix<> u1;

        for (index_t k = 0; k<=d; ++k)
        {
            gb[k].anchors_into(tmp);
            glpt.resize(d, tmp.cols());
            glpt.topRows(l) = tmp.topRows(l);
            glpt.row(l).setConstant( iFace.first().parameter() ? 1 : 0); // l==2 and l2==2 usually
            glpt.bottomRows(d-l-1) = tmp.bottomRows(d-l-1);

            gleval.resize(1, tmp.cols()); // interpolation values

            for ( index_t i = 0; i != glpt.cols(); ++i )
            {
                u1 = glpt.col(i);
                ifaceMap.eval_into(u1, u2);
                P1.jacobian_into (u1, ev1);
                P2.jacobian_into (u2, ev2);

                if (d==k)
                {
                    gleval.at(i)  = sgn * ev1.determinant(); // J44
                }
                else
                {
                    // get submatrix
                    D0.col(d-1) = ev2.col(l2);
                    for (index_t t = 0; t!=d-1; ++t)
                        D0.col(t) = ev1.col(t + (index_t)(t>=k));

                    gleval.at(i)  = sgn * D0.determinant();  // J4k (k=0,1,2)
                }
            }
            gdata[k] = gb[k].interpolateAtAnchors(gleval);
            // Add perturbation...
            if ( ! pert1.empty() )
                gdata[k]->coefs() +=  sgn * pert1[k]->coefs(); // assuming the same pert. degree

            //sgn *= -1;
        }
    }


    void eval_into(const gsMatrix<> & u, gsMatrix<> & result) const
    {
        GISMO_ASSERT(u.cols()==1, "Expecting one point");
        result.resize(d+1, u.cols());

        gsMatrix<> tmp;
        for (index_t k = 0; k<=d; ++k)
        {
            gdata[k]->eval_into(u.topRows(d-1), tmp);
            result.at(k) = tmp.value();
        }
    }

    index_t parDim() const {return d;}

    const gsGeometry<> & get(const size_t k) const {return *gdata[k];}
    gsGeometry<> & get(const size_t k) {return *gdata[k];}

    const gsMultiBasis<> & basis() const {return mbasis;}
    gsMultiBasis<> & basis()             {return mbasis;}


    // load / save

public:

    index_t d;

    gsMultiBasis<> mbasis;

    // gluing data
    std::vector<gsGeometry<>::Ptr> gdata;
};


void storeBasis(const gsMatrix<> & ker, gsMultiPatch<> & mp,
                const gsDofMapper & map);


void basisDegree3(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result);

void basisDegree4(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result);

void basisPlnDegree3(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result);

void basisPlnDegree4(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result);

void basisLinDegree3(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result,
                  const std::vector<gsGeometry<>::Ptr>& pert);

void basisLinDegree4(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result,
                  const std::vector<gsGeometry<>::Ptr>& pert);

int main(int argc, char *argv[])
{
    const unsigned sd = (unsigned int)time(0);
    gsInfo <<"Random seed: "<< sd <<"\n";
    std::srand(sd); // random generator

    //gsInfo <<"GMP ver: " << gmp_version <<"\n";

    //std::string fn("planar/two_squares.xml");
    //std::string fn("planar/two_quadrilateralPatches.xml");
    //std::string fn("volumes/two_cuboids.xml");
    std::string fn("volumes/two_cuboids2.xml");
    //std::string fn("volumes/two_cubes.xml");
    //std::string fn("volumes/two_cuboidsPlanarBoundary.xml");
    index_t deg(3);
    index_t knot(0);
    index_t mult(-1);
    bool plot = false; // defaults to false
    index_t numSamples = 1000;
    bool plot_mesh = false;
    bool plot_net = false;
    bool c = false;
    bool getbasis = false, genericGeo = false, genGluing = true, pertInput = false, pertGluing = false;
    bool bdr     = false;

    gsCmdLine cmd("Hi, give me a file (eg: .xml) and I will try to draw it!");
    cmd.addPlainString("filename", "File containing the geometry", fn);
    cmd.addInt   ("d", "degree", "Degree to work with", deg);
    cmd.addInt   ("k", "knots", "Number of inserted knots", knot);
    cmd.addInt   ("m", "mult", "Multiplicity of the knots to insert", mult);
    cmd.addSwitch("plot"  , "Plot the result", plot);
    cmd.addSwitch("dim"  , "Compute dimension formulas", c);
    cmd.addSwitch("perturb"  , "Perturb the input", pertInput);
    cmd.addSwitch("pertGluing", "Perturb the gluing data", pertGluing);
    cmd.addSwitch("getBasis", "Computes a basis", getbasis);
    cmd.addSwitch("genericGeo"  , "Compute dimension of generic input", genericGeo);
    cmd.addSwitch("genericGluing"  , "Compute dimension for generic gluing data", genGluing);
    cmd.addSwitch("boundary"  , "Compute boundary functions", bdr);
    cmd.addInt   ("s", "samples", "Number of samples to use for viewing", numSamples);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (-1==mult) mult = deg-1; // if mult not given default to C^1

    gsFileData<>  filedata(fn);
    /*
    if ( filedata.has<gsSparseMatrix<> >() )
    {
        gsSparseMatrix<> gb;
        gsMultiBasis<> mb;
        filedata.getFirst(gb);
        gsInfo<< "Got matrix "<< gb.dim() <<".\n";
        filedata.getFirst(mb);
        gsInfo<< "Got "<< mb <<"\n";
        gsMultiBasis<> mb2 = mb;
        const gsBasis<> & b0 = mb.basis(0).component(0);
        const index_t nKnots = (b0.size() - b0.degree(0)-1) / (b0.degree(0)-1);
        mb2.uniformRefine(knot-nKnots, math::min(mult,deg));
        return 0;
    }
    */

    gsMultiPatch<> mp;
    filedata.getFirst(mp);
    mp.computeTopology();
    const boundaryInterface & iFace = *mp.iBegin();
    gsInfo<< "Got "<< mp <<"\n";
    gsInfo<< iFace <<"\n";

    // Get the degree of the input
    index_t pdeg = mp.patch(0).basis().maxDegree();

    mb = gsMultiBasis<>(mp);
    mb.degreeElevate(deg-pdeg);
    mb.uniformRefine(knot, math::min(mult,deg));
    //mb.uniformRefineComponent(0, knot, math::min(mult,deg));
    //mb.uniformRefineComponent(1, knot, math::min(mult,deg));
    // Here: degree-elevate the patches
    //mp.degreeElevate(deg-pdeg);
    // Here: knot insertion
    //mp.uniformRefine(knot, math::min(mult,deg));

    gsInfo<< "Initial configuration (patch 0): "<< mp[0] <<"\n";

    if(pertInput)
        perturbGeometry(mp, 0.001, true);

    if (c)
    {
        computeDimFormula(mp, false);
        return 0;
    }

    std::vector<gsGeometry<>::Ptr> pert;

    if (genGluing)
    {
        if (getbasis)
        {
            /* // Predefined
              gsFileData<> fd("ggdata.xml");
            pert.resize(4);
            pert[0] = fd.getId<gsGeometry<> >(0);
            pert[1] = fd.getId<gsGeometry<> >(1);
            pert[2] = fd.getId<gsGeometry<> >(2);
            pert[3] = fd.getId<gsGeometry<> >(3);
            */

            // Random
            getPerturbation( pert, mp.parDim() );
        }
        else
        {
            GISMO_ASSERT(0==knot, "Beware of initial knot insertion, got "<< knot);// table starts with zero knots...
            genericGluingDim(mp, 0, 20); // note: second arg. unused
            return 0;
        }
    }

    if (genericGeo)
    {
        GISMO_ASSERT(0==knot, "Beware of initial knot insertion, got "<< knot);// table starts with zero knots...
        // Note: deg==1  and pdeg==1 is the trilinear case
        genericGeometryDim(mp, 0.001, true, 10);
        return 0;
    }

    gsInfo << "Working basis:\n"<< mb[0] <<"\n";
    index_t kdim, odim;
    gsDofMapper map;
    gsMatrix<> Ker;
 //   gsVector<index_t> count;

    if (getbasis)
    {
  /*
        if ( 2 == mp.parDim() )
        localBasis_incr<2>(mp, pdeg, pert, map, Ker, bdr);
       // localBasis_incr<2>(mp, pdeg, map, Ker, bdr);
        else //3D
        localBasis_incr<3>(mp, pdeg, pert, map, Ker, bdr);
  */
 ///*
        gsSparseMatrix<> Bsparse;
        if ( 3 == mp.parDim() )
        {
            if (3==deg)
            {
                gsGluingData gd(mp/*, pdeg*/);
                Bsparse = FullBasisNonHom(mp, /*pdeg,*/ gd ,map);
            }
            //if (3==deg) basisLinDegree3(mp, pdeg, Bsparse, pert);
            //if (3==deg) basisDegree3(mp, pdeg, Bsparse);
            //if (4==deg) basisDegree4(mp, pdeg, Bsparse);
            if (4==deg) basisLinDegree4(mp, pdeg, Bsparse, pert);
        }
        else
        {
            if (3==deg) basisPlnDegree3(mp, pdeg, Bsparse);
            if (4==deg) basisPlnDegree4(mp, pdeg, Bsparse);
        }
        // Export XML
        gsFileData<> xml;
        xml << mb;
        xml << Bsparse;
        const index_t nKnots0 = (mb[0].component(0).size()-mb[0].component(0).degree(0)-1) / (mb[0].component(0).degree(0)-1);

        std::string name = "deg"+ util::to_string(deg) +"BasisLin_k" + util::to_string(nKnots0) + "_Map";
        //std::string name = "deg"+ util::to_string(deg) +"Basis_k" + util::to_string(nKnots0) + "_Map";

        xml.saveCompressed(name);
        gsInfo<< "dim="<<Bsparse.cols()<<", written to file "<<name<<".xml.gz \n";
        return 0;
// */

        kdim = Ker.cols();
        odim = map.freeSize();

        sortBySparsity(Ker);
        printSparsity(Ker, map, mp);
/*
        if ( kdim < 100 )
        {
            sortBySparsity(Ker);
            printSparsity(Ker, map, mp);
        }

        //selectRepresentatives(Ker, count);
        gsInfo << "\n ---- Summarized ---- \n\n";
        gsInfo << "Count: "<< count.transpose() << " (total "<<count.size()<<" kinds)\n";
        printSparsity(Ker, map, mp);
        gsInfo << "Coefs: "<< countSparsity(Ker).transpose() <<"\n";
        gsInfo << "Count: "<< count.transpose() << " (total "<<count.size()<<" kinds)\n";

        const index_t onInt = mp.patch(mp.iBegin()->first().patch).basis()
            .boundary(iFace.first().side()).size();
        gsInfo << "Rank            : "<< Ker.rows() + onInt - kdim <<"\n";

         Normalize and Store
         for( index_t i = 0; i!= Ker.cols(); ++i) //normalize rationally
         {
             Ker.col(i).array() /= Ker.col(i).array().abs().maxCoeff();
             if ( (Ker.col(i).array() <= 0 ).all() )
                 Ker.col(i).array() *= -1;
         }

        Ker.array() *= 2;
 */
       // storeBasis(Ker, mp, map);

    }
    else
    {
        std::pair<index_t,index_t> sdim = computeDim(mp, pertGluing);
        kdim = sdim.first;
        odim = sdim.second;
    }

    gsInfo << "Kernel dimension: "<< kdim <<"\n";
    gsInfo << "Final dimension : "<< kdim + odim <<"\n";

    if (plot)
    {
        GISMO_ENSURE(getbasis, "Basis not computed!");

        gsInfo << "Plotting in ParaView..\n";

        if (2==mp.parDim())
        {
            PlotFunctions(Ker, mp, map, numSamples, plot_mesh, plot_net);
            plotSparsity<2>(Ker, map, mp);
        }
        else
        {
            PlotFunctions(Ker, mp, map, numSamples, false, false);
            plotSparsity<3>(Ker, map, mp);
        }

        return system("paraview gcall.pvd &");
    }

    return 0;
}


void PlotFunctions(const gsMatrix<> & ker, gsMultiPatch<> & mp,
                   const gsDofMapper & map, const int & ns, const bool & pm, const bool & pn)
{
    gsFileData<> xml;
//    const boundaryInterface & iFace = *mp.iBegin();

//    const gsBasis<> & B1 = mb[0];
//    const gsBasis<> & B2 = mb[1];
//    const gsTensorBSplineBasis<3> & B1 =
//    static_cast<const gsTensorBSplineBasis<3> &>(mb[0]);
//    const gsTensorBSplineBasis<3> & B2 =
//    static_cast<const gsTensorBSplineBasis<3> &>(mb[1]);

    gsParaviewCollection col("gcall");
    const boundaryInterface & iFacePlot = *mp.iBegin();
    gsVector<index_t> p(2);
    p(0) = iFacePlot.first() .patch;
    p(1) = iFacePlot.second().patch;

    const index_t d = mp.parDim();
    mp.patch(0).embed(d+1);
    mp.patch(1).embed(d+1);
    for ( index_t k = 0; k != ker.cols(); ++k) // for all basis function
    {
        for ( size_t i = 0; i< mp.nPatches(); ++i) // for all patches
        {
            gsMatrix<> & coefs = mp.patch(p(i)).coefs();
            const index_t sz = coefs.rows();
            for ( index_t j = 0; j != sz; ++j)
            {
                const index_t jj = map.index(j,p(i));
                if ( map.is_boundary_index(jj) )// control point is in the matrix
                {
                    const int bjj = map.global_to_bindex(jj);
                    coefs(j,d) = ker(bjj,k);
                    //coefs(j,d) = ( 0!=ker(bjj,k) ? 1 : 0) ;
                    //coefs(j,d) = 1;
                    //For Debugging
                   /* if ( 0 != ker(bjj,k) )
                        gsInfo<< "patch="<<p(i)<<"j="<<j<<", ker="<<bjj<< ", Coeff=("<<
                                 (0==p(i) ? B1.tensorIndex(j).transpose() :
                                            B2.tensorIndex(j).transpose())<<
                                    ") is "<< ker(bjj,k) <<"\n";
                                    */
                }
            }
        }

        // Export XML
        xml << mp;

        // Export Paraview
        std::string gcBasisFct = "gcall" + util::to_string(k);
        gsWriteParaview(mp, gcBasisFct, ns, pm, pn);
        col.addTimestep(gcBasisFct, p(0), k, ".vts");
        col.addTimestep(gcBasisFct, p(1), k, ".vts");
        if(pm)
        {
            col.addTimestep(gcBasisFct, p(0), k, "_mesh.vtp");
            col.addTimestep(gcBasisFct, p(1), k, "_mesh.vtp");
        }
        if(pn)
        {
            col.addTimestep(gcBasisFct, p(0), k, "_cnet.vtp");
            col.addTimestep(gcBasisFct, p(1), k, "_cnet.vtp");
        }
    }
    col.save();

    xml.save("representative_basis");
}


void storeBasis(const gsMatrix<> & ker, gsMultiPatch<> & mp,
                const gsDofMapper & map)
{
    gsMultiBasis<> & newBasis = mb;
    const gsBasis<> & B0 = newBasis.basis(0);
    const gsBasis<> & B1 = newBasis.basis(1);
    const index_t nKnots0 = (B0.component(0).size() - B0.component(0).degree(0)-1) / (B0.component(0).degree(0)-1);

    gsSparseMatrix<> bMatrix(B0.size() + B1.size(), ker.cols());

    const boundaryInterface & iFacePlot = *mp.iBegin();
    gsVector<index_t> p(2);
    p(0) = iFacePlot.first() .patch;
    p(1) = iFacePlot.second().patch;

    for ( index_t k = 0; k != ker.cols(); ++k) // for all basis function
    {
        for ( size_t i = 0; i< mp.nPatches(); ++i) // for all patches
        {
            const index_t start = ( 0 == p(i) ? 0 : newBasis.basis(1).size() );
            const index_t sz    = newBasis.basis(p(i)).size();

            for ( index_t j = 0; j != sz; ++j)
            {
                const index_t jj = map.index(j,p(i));
                if ( map.is_boundary_index(jj) )// control point is in the matrix
                {
                    const int bjj = map.global_to_bindex(jj);
                    if ( ker(bjj,k) != 0 )
                        bMatrix.insert(start+j, k)   = ker(bjj,k);
                }
            }
        }

        // Export XML
        gsFileData<> xml;
        //gsMatrix<> test = bMatrix.toDense();
        //xml << test;
        xml << newBasis;
        xml << bMatrix;
        std::string name = "basis_d" + util::to_string(B0.component(0).degree(0))
            + "_k" + util::to_string(nKnots0);
        xml.save(name);
    }
}

void G1Matrix(const gsMultiPatch<> & mp,
              const index_t pdeg,
              const std::vector<gsGeometry<>::Ptr> & pert1,
              MatType & mat,
              gsDofMapper & map)
{
    const index_t d = mp.parDim();
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsGeometry<> & P1 = mp.patch( iFace.first().patch );
    const gsGeometry<> & P2 = mp.patch( iFace.second().patch );
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));

    gsMatrix<> pt1, gr, u2, ev1, ev2, D0(d,d);
    gsMatrix<index_t> act1, act2;

    gsVector<index_t> sz(2);
    sz[0] = B1.size();
    sz[1] = B2.size();
    map = gsDofMapper(sz);
    // glue interface
    B1.matchWith(iFace, B2, act1, act2);
    map.matchDofs(iFace.first().patch, act1, iFace.second().patch, act2);
    // mark dofs
    map.markBoundary(iFace.first().patch, act1);//interface
    act1 = B1.boundaryOffset(iFace.first() .side(), 1);
    map.markBoundary(iFace.first().patch, act1); //first adj. face
    act2 = B2.boundaryOffset(iFace.second().side(), 1);
    map.markBoundary(iFace.second().patch, act2);//second adj. face
    map.finalize();

    // the determinant lives in this basis
    const index_t l = iFace.first()  .direction(); // "orthogonal" direction p1
    const index_t l2 = iFace.second().direction(); // "orthogonal" direction p2
    gsBasis<>::uPtr detB = B1.boundaryBasis(iFace.first().side());
    detB->reduceContinuity(1);
    detB->degreeElevate(d*pdeg-1);
    detB->anchors_into(gr);
    pt1.resize(d, gr.cols());
    pt1.topRows(l) = gr.topRows(l);
    pt1.row(l).setConstant( iFace.first().parameter() ? 1 : 0); // l==2 usually
    pt1.bottomRows(d-l-1) = gr.bottomRows(d-l-1);
    //gsInfo<<"Points:\n"<< pt1<<"\n";

    const index_t N = gr.cols();
    const index_t M = map.boundarySize();

    #if useSparse
    mat.resize(N, M);
    index_t nn = 2 * detB->degree(0) + 1;
    nn *= 2 * detB->degree(1) + 1;
    mat.reservePerColumn(nn);
    #else
    mat.setZero(N, M);
    #endif

    gsMatrix<> u1, tmp;// temporary variable
    gr.resize(d+1,1); //gr[0]... gr[d]  --> values of alpha, beta, gamma, delta
    for ( index_t i = 0; i != N; ++i)
    {
        u1 = pt1.col(i);
        ifaceMap.eval_into(u1, u2);
        P1.jacobian_into (u1, ev1);
        P2.jacobian_into (u2, ev2);
        //gsInfo<<"Pair: ("<< u1.transpose() <<"), ("<<u2.transpose()<<")\n";

        // G1 = (g11, g12, w13)
        // G2 = (g21 ,g22, w23)
        //
        // | dxg11 dyg11 dxg21 |
        // | dxg12 dyg12 dxg22 | = 0
        // | dxw13 dyw13 dxw23 |
        // eg. x= l, y= !l

        // first case: coming from Geometry mapping
        int sgn = 1;
        D0.col(d-1) = ev2.col(l2);
        for (index_t k = 0; k!=d; ++k) //last row
        {
            for (index_t t = 0; t!=d-1; ++t)
                D0.col(t) = ev1.col(t + (index_t)(t>=k));
            gr.at(k)  = sgn * D0.determinant();  // J41, J42, J43
            sgn *= -1;
        }
        gr.at(d)  = sgn * ev1.determinant(); // J44

        // gsDebugVar(gr.transpose());
        // gd.eval_into(u1, gr2);
        // gsDebugVar(gr2.transpose());

        if ( !pert1.empty() )
        {
            gsDebugVar(*pert1.front());
            gsDebugVar(u1);
            //--- Perturbation with generic values
            for (index_t k = 0; k!=d; ++k) //last row
            {
                pert1[k]->eval_into(u1.topRows(d-1), tmp); // polynomial of same degree as a_k
                gr.at(k) += tmp(0,0); // using the first coordinate
            }

            pert1[d]->eval_into(u2.topRows(d-1), tmp);
            gr.at(d) += tmp(0,0);
        }

        B1.deriv_into (u1, ev1 );
        B1.active_into(u1, act1);
        map.localToGlobal(act1, iFace.first().patch, act1);
        B2.deriv_into (u2, ev2 );
        B2.active_into(u2, act2);
        map.localToGlobal(act2, iFace.second().patch, act2);

        for ( index_t j = 0; j != act1.rows(); ++j) // ==act2.rows()
        {
            const index_t jj1 = act1.at(j);
            if ( map.is_boundary_index(jj1) )
            {
                const index_t bjj = map.global_to_bindex(jj1);
#if useSparse
                mat.addTo(i, bjj, gr.col(0).topRows(d).dot( ev1.col(0).segment(d*j,d) ) );
#else
                mat(i, bjj) +=  gr.col(0).topRows(d).dot( ev1.col(0).segment(d*j,d) );
#endif
            }
            const index_t jj2 = act2.at(j);
            if ( map.is_boundary_index(jj2) )
            {
                const index_t bjj = map.global_to_bindex(jj2);
#if useSparse
                mat.addTo(i, bjj, ev2.at(d*j+l2) * gr.at(d) );
#else
                mat(i, bjj) += ev2.at(d*j+l2) * gr.at(d);
#endif
            }
        }
    }

    #if useSparse
    // Compress sparse matrix
    mat.makeCompressed();
    #endif
}

void G1Matrix(const gsMultiPatch<> & mp,
              const index_t pdeg,
              const gsGluingData & gd,
              MatType & mat,
              gsDofMapper & map)
{
    const index_t d = mp.parDim();
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));

    gsMatrix<> pt1, gr, u2, ev1, ev2;
    gsMatrix<index_t> act1, act2;

    gsVector<index_t> sz(2);
    sz[0] = B1.size();
    sz[1] = B2.size();
    map = gsDofMapper(sz);
    // glue interface
    B1.matchWith(iFace, B2, act1, act2);
    map.matchDofs(iFace.first().patch, act1, iFace.second().patch, act2);
    // mark dofs
    map.markBoundary(iFace.first().patch, act1);//interface
    act1 = B1.boundaryOffset(iFace.first() .side(), 1);
    map.markBoundary(iFace.first().patch, act1); //first adj. face
    act2 = B2.boundaryOffset(iFace.second().side(), 1);
    map.markBoundary(iFace.second().patch, act2);//second adj. face
    map.finalize();
    //map.print();

    // the determinant lives in this basis
    const index_t l = iFace.first()  .direction(); // "orthogonal" direction p1
    const index_t l2 = iFace.second().direction(); // "orthogonal" direction p2
    gsBasis<>::uPtr detB = B1.boundaryBasis(iFace.first().side());
    detB->reduceContinuity(1);
    detB->degreeElevate(d*pdeg-1);
    detB->anchors_into(gr);
    pt1.resize(d, gr.cols());
    pt1.topRows(l) = gr.topRows(l);
    pt1.row(l).setConstant( iFace.first().parameter() ? 1 : 0); // l==2 usually
    pt1.bottomRows(d-l-1) = gr.bottomRows(d-l-1);
    //gsInfo<<"Points:\n"<< pt1<<"\n";

    const index_t N = gr.cols();
    const index_t M = map.boundarySize();

    #if useSparse
    mat.resize(N, M);
    //mat.reservePerColumn(detB->totalDegree());// to improve using numActive()
    mat.reserve(gsVector<index_t>::Constant(mat.rows(), detB->totalDegree()));
    #else
    mat.setZero(N, M);
    #endif

    gsMatrix<> u1;
    gr.resize(d+1,1); //gr[0]... gr[d]  --> values of alpha, beta, gamma, delta
    for ( index_t i = 0; i != N; ++i)
    {
        u1 = pt1.col(i);
        gd.eval_into(u1, gr);
        B1.deriv_into (u1, ev1 );
        B1.active_into(u1, act1);
        map.localToGlobal(act1, iFace.first().patch, act1);
        ifaceMap.eval_into(u1, u2);
        B2.deriv_into (u2, ev2 );
        B2.active_into(u2, act2);
        map.localToGlobal(act2, iFace.second().patch, act2);

        for ( index_t j = 0; j != act1.rows(); ++j) // ==act2.rows()
        {
            const index_t jj1 = act1.at(j);
            if ( map.is_boundary_index(jj1) )
            {
                const index_t bjj = map.global_to_bindex(jj1);
#if useSparse
                mat.addTo(i, bjj, gr.col(0).topRows(d).dot( ev1.col(0).segment(d*j,d) ) );
#else
                mat(i, bjj) +=  gr.col(0).topRows(d).dot( ev1.col(0).segment(d*j,d) );
#endif
            }
            const index_t jj2 = act2.at(j);
            if ( map.is_boundary_index(jj2) )
            {
                const index_t bjj = map.global_to_bindex(jj2);
#if useSparse
                mat.addTo(i, bjj, ev2.at(d*j+l2) * gr.at(d) );
#else
                mat(i, bjj) += ev2.at(d*j+l2) * gr.at(d);
#endif
            }
        }
    }

    #if useSparse
    // Compress sparse matrix
    mat.makeCompressed();
    #endif
}

void G1Matrix(const gsMultiPatch<> & mp,
              const index_t pdeg,
              MatType & mat,
              gsDofMapper & map, bool randomize)
{
    const index_t d = mp.parDim();
    std::vector<gsGeometry<>::Ptr> pert1;
    real_t tol = 10e-1; // todo: add as option if needed
    if (randomize)
        getPerturbation(pert1, d, tol);

    //gsGluingData gd(mp, /*pdeg,*/ pert1);
    //G1Matrix(mp, pdeg, gd, mat, map);
    G1Matrix(mp, pdeg, pert1, mat, map);
}
gsMatrix<> eval_R(const gsKnotVector<> & m_knots,
                       const gsMatrix<> & u)
{
    const index_t m_p = m_knots.degree();
    STACK_ARRAY(real_t, left, m_p + 1);
    STACK_ARRAY(real_t, right, m_p + 1);
    gsMatrix<> result(m_p+1, u.cols() );

    for (index_t v = 0; v < u.cols(); ++v) // for all columns of u--gsMatrix
    {
        // Get span of absissae
        unsigned span = m_knots.iFind( u(0,v) ) - m_knots.begin() ;
        result(0,v)= 1;  // 0-th degree function value

        for(int j=1; j<= m_p; j++) // For all degrees
        {
            left[j]  = u(0,v) - m_knots[span+1-j];
            right[j] = m_knots[span+j] - u(0,v);
            real_t saved = 0;

            for(int r=0; r<j ; r++)
            {
                const real_t temp = result(r,v) / ( right[r+1] + left[j-r] );
                result(r,v)     = saved + (j==m_p?real_t(1):right[r+1]) * temp ;
                saved = (j==m_p?real_t(1):left[j-r]) * temp ;
            }
            result(j,v)     = saved;
        }

    }// end for all columns v
    return result;
}

gsSparseMatrix<> FullBasisNonHom(const gsMultiPatch<> & mp,
                            // const index_t pdeg,
                            const gsGluingData & gd,
                            gsDofMapper & map)
{
    gsSparseMatrix<> mat;
    gsMatrix<> b;

    const index_t d = mp.parDim();
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];

    gsMatrix<> pts1, pts2, gr, u2, ev1, ev2;
    gsMatrix<index_t> act1, act2;

    gsVector<index_t> sz(2);
    sz[0] = B1.size();
    sz[1] = B2.size();
    map = gsDofMapper(sz);
    // glue interface
    B1.matchWith(iFace, B2, act1, act2);
    map.matchDofs(iFace.first().patch, act1, iFace.second().patch, act2);
    // mark dofs
    map.markBoundary(iFace.first().patch, act1);//interface
    act1 = B1.boundaryOffset(iFace.first() .side(), 1);
    map.markBoundary(iFace.first().patch, act1); //first adj. face
    act2 = B2.boundaryOffset(iFace.second().side(), 1);
    map.markBoundary(iFace.second().patch, act2);//second adj. face
    map.finalize();

    gsMatrix<index_t> r1, r2, r3;
    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    // the determinant lives in this basis
    gsBasis<>::uPtr basis1 = B1.clone();
    gsTensorBSplineBasis<3> * bbasis1 = dynamic_cast<gsTensorBSplineBasis<3> *>(basis1.get());
    bbasis1->component(0).degreeReduce(2);
    bbasis1->component(1).degreeReduce(2);
    gsBasis<>::uPtr basis2 = B2.clone();
    gsTensorBSplineBasis<3> * bbasis2 = dynamic_cast<gsTensorBSplineBasis<3> *>(basis2.get());
    bbasis2->component(0).degreeReduce(2);
    bbasis2->component(1).degreeReduce(2);
    gsVector<index_t,3> ti;
    index_t t0sz = bbasis1->component(0).size();

    pts1.resize(d,r1.size()+r2.size());
    pts2.resize(d,r3.size());
    for(index_t k = 0; k!=r1.size(); ++k)
        pts1.col(k) = B1.anchor(r1.at(k));
    for(index_t k = 0; k!=r2.size(); ++k)
        pts1.col(k+r1.size()) = B1.anchor(r2.at(k));
    for(index_t k = 0; k!=r3.size(); ++k)
        pts2.col(k) = B2.anchor(r3.at(k));

    const index_t N = map.boundarySize(); // = r1.size()+r2.size()+r3.size()

    mat.resize(N, N);
    mat.reserve(gsVector<index_t>::Constant(mat.rows(), B1.totalDegree()));

    index_t M = bbasis1->boundary(iFace.first().side()).size();
    b.setZero(N, M);

    gsMatrix<> u1;
    gr.resize(d+1,1); //gr[0]... gr[d]  --> values of alpha, beta, gamma, delta

    gsMatrix<> corner;
    corner.resize(d,d);
    for(index_t k = 0; k!=d; ++k)
        corner.row(k) = (mp.patch(1).coefAtCorner(k+2)-mp.patch(1).coefAtCorner(1));

    // real_t vol = corner.determinant();

    const gsKnotVector<> & knots =//same knots in every direction
        dynamic_cast<const gsTensorBSplineBasis<3>&>(B1).knots(1);
    gsMatrix<> R = eval_R(knots, pts1);
    gsDebugVar(pts1);
    gsDebugVar(R);

    for ( index_t i = 0; i != pts1.cols(); ++i)
    {
        u1 = pts1.col(i);
        gd.eval_into(u1, gr);
        B1.eval_into (u1, ev1 );
        B1.active_into(u1, act1);
        map.localToGlobal(act1, iFace.first().patch, act1);

        // Matrix
        for ( index_t j = 0; j != act1.rows(); ++j) // ==act2.rows()
        {
            const index_t jj1 = act1.at(j);
            if ( map.is_boundary_index(jj1) )
            {
                const index_t bjj = map.global_to_bindex(jj1);
                mat.addTo(i, bjj, ev1.at(j) );
            }
        }

        bbasis1->eval_into  (u1, ev1);
        bbasis1->active_into(u1, act1);

        // Vector
        for ( index_t j = 0; j != act1.rows(); ++j) // ==act2.rows()
        {
            ti = bbasis1->tensorIndex(act1.at(j));
            if ( 1 == ti[2] )
                b(i, ti[0]*t0sz + ti[1] ) =  gr.at(2) * ev1.at(j);
        }
    }

    for ( index_t i = 0; i != pts2.cols(); ++i)
    {
        u2 = pts2.col(i);
        gd.eval_into(u2, gr);
        B2.eval_into (u2, ev2 );
        B2.active_into(u2, act2);
        map.localToGlobal(act2, iFace.second().patch, act2);

        // Matrix
        for ( index_t j = 0; j != act2.rows(); ++j) // ==act2.rows()
        {
            const index_t jj1 = act2.at(j);
            if ( map.is_boundary_index(jj1) )
            {
                const index_t bjj = map.global_to_bindex(jj1);
                mat.addTo(pts1.cols()+i, bjj, ev2.at(j) );
            }
        }

        bbasis2->eval_into  (u2, ev2);
        bbasis2->active_into(u2, act2);
        // Vector
        for ( index_t j = 0; j != act2.rows(); ++j) // ==act2.rows()
        {
            ti = bbasis2->tensorIndex(act2.at(j));
            if ( 1 == ti[2] )
                b(pts1.cols()+i, ti[0]*t0sz + ti[1] ) =  gr.at(3) * ev2.at(j) ;
        }
    }

    mat.makeCompressed();

    gsSparseSolver<>::QR slv;
    gsSparseMatrix<> result = slv.compute(mat).solve(b).sparseView();
    result.prune(1e-10,1);

    gsWrite(result, "phi_1");


    /*
    gsBasis<>::uPtr basis1 = B1.clone();
    gsTensorBSplineBasis<3> * bbasis1 = dynamic_cast<gsTensorBSplineBasis<3> *>(basis1.get());
    bbasis1->component(0).degreeReduce(2);
    bbasis1->component(1).degreeReduce(2);
    gsBasis<>::uPtr basis2 = B2.clone();
    gsTensorBSplineBasis<3> * bbasis2 = dynamic_cast<gsTensorBSplineBasis<3> *>(basis2.get());
    bbasis2->component(0).degreeReduce(2);
    bbasis2->component(1).degreeReduce(2);

    b.setZero(N, M);

    for ( index_t i = 0; i != pts1.cols(); ++i)
    {
        u1 = pts1.col(i);
        gd.eval_into(u1, gr);
        bbasis1->eval_into  (u1, ev1);
        bbasis1->active_into(u1, act1);

        // Vector
        for ( index_t j = 0; j != act1.rows(); ++j) // ==act2.rows()
        {
            ti = bbasis1->tensorIndex(act1.at(j));
            if ( 1 == ti[2] )
                b(i, ti[0]*t0sz + ti[1] ) =  gr.at(2) * ev1.at(j);
        }
    }

    for ( index_t i = 0; i != pts2.cols(); ++i)
    {
        u2 = pts2.col(i);
        gd.eval_into(u2, gr);
        bbasis2->eval_into  (u2, ev2);
        bbasis2->active_into(u2, act2);
        // Vector
        for ( index_t j = 0; j != act2.rows(); ++j) // ==act2.rows()
        {
            ti = bbasis2->tensorIndex(act2.at(j));
            if ( 1 == ti[2] )
                b(pts1.cols()+i, ti[0]*t0sz + ti[1] ) =  gr.at(3) * ev2.at(j) ;
        }
    }

    result = slv.solve(b).sparseView();
    result.prune(1e-10,1);

    gsWrite(result, "phi_2");
    */
    return result;
}


gsVector<index_t> countSparsity(const gsMatrix<> & K)
{
    gsVector<index_t> c(K.cols());
    for (index_t j = 0; j!=K.cols(); ++j)
        c[j] = (K.col(j).array()!=0).count();
    return c;
}

std::string printSparsity(const gsMatrix<> & M)
{
    std::ostringstream os;
    for (index_t i = 0; i!=M.rows(); ++i)
    {
        for (index_t j = 0; j!=M.cols(); ++j)
            os<< ( 0 == M(i,j) ? "\u00B7" : "x");
        os<<"  "<<(M.row(i).array()!=0).count()<<"\n";
    }
    return os.str();
}

void printSparsity(const gsMatrix<> & Ker,
                   const gsDofMapper & map,
                   const gsMultiPatch<> & mp)
{
    //const index_t d = mp.parDim();
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    gsVector<index_t> p(2), nz(3);
    p(0) = iFace.first() .patch;
    p(1) = iFace.second().patch;

    gsMatrix<index_t> r1, r2, r3;
    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    //const index_t v = ( 2 == B1.domainDim() ? 1 : B1.component(0).size() );
    const index_t d = B1.domainDim();
    const index_t v = B1.component(0).size();

    const index_t u = ( 2 == d ? 1 : v);

    gsMatrix<int> S(v, 3*u);
    for (index_t k = 0; k!=Ker.cols(); ++k)
    {
        index_t m = 0;
        for (index_t j = 0; j!=r1.size(); ++j)
            S.at(m++) = math::getSign( Ker(map.bindex(r1.at(j),p(0)), k) );
        for (index_t j = 0; j!=r2.size(); ++j)
            S.at(m++) = math::getSign( Ker(map.bindex(r2.at(j),p(0)), k) );
        for (index_t j = 0; j!=r3.size(); ++j)
            S.at(m++) = math::getSign( Ker(map.bindex(r3.at(j),p(1)), k) );

        nz[0] = (S.leftCols(u)    .array()!=0).count() ;
        nz[1] = (S.middleCols(u,u).array()!=0).count() ;
        nz[2] = (S.rightCols (u)  .array()!=0).count() ;
        if ( 3==d)
        {
            S.leftCols  (u)  .transposeInPlace();
            S.middleCols(u,u).transposeInPlace();
            S.rightCols (u)  .transposeInPlace();
        }

        gsInfo<<"Basis function "<< k <<"(nnz="<< nz.transpose() <<")\n";
        for (index_t i = 0; i!=S.rows(); ++i)
        {
            for (index_t j = 0; j!=S.cols(); ++j)
            {
                gsInfo << ( j!=0 && j%v == 0  ? " | " : "");
                //gsInfo << ( 0 == S(i,j) ? "\u00B7" : "x");
                gsInfo << ( 0 > S(i,j) ? "x" : ( 0 < S(i,j) ? "+" : "\u00B7") );
            }
            //gsInfo <<"  "<<(M.row(i).array()!=0).count()<<"\n";
            gsInfo << "\n";
        }
    }
}


template<int d>
void plotSparsity(const gsMatrix<> & Ker, const gsDofMapper & map, gsMultiPatch<> & mp)
{
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsTensorBasis<d> & B1 = static_cast<const gsTensorBasis<d>&>(mb[0]);
    const gsTensorBasis<d> & B2 = static_cast<const gsTensorBasis<d>&>(mb[1]);
    gsVector<index_t> p(2);
    p(0) = iFace.first() .patch;
    p(1) = iFace.second().patch;
    gsVector<index_t, d> tt;

    gsParaviewCollection col("sparsity");

    const index_t sz1 = B1.size();
    const index_t sz2 = B1.size();
    gsMatrix<> CP(4, sz1), CP2(4, sz1+sz2);

    for (index_t j = 0; j!= Ker.cols(); ++j)
    {
        for (index_t k = 0; k!= sz1; ++k)
        {
            tt = B1.tensorIndex(k);
            CP.col(k).head<d>() = tt.template cast<real_t>();
            const index_t jj = map.index(k,p(0));
            if ( map.is_boundary_index(jj)  && 0 != Ker(map.global_to_bindex(jj), j)  )
                CP(d,k) = 0 ;
            else
                CP(d,k) = 1;
        }

        CP2.leftCols(sz1) = CP;

        for (index_t k = 0; k!= sz2; ++k)
        {
            tt = B2.tensorIndex(k);
            CP.col(k).head<d>() = tt.template cast<real_t>();
            const index_t jj = map.index(k,p(1));
            if ( map.is_boundary_index(jj) && 0 != Ker(map.global_to_bindex(jj), j)  )
                CP(d,k) = 0 ;
            else
                CP(d,k) = -1;
        }

        CP.row(2).array() *= -1; // mirror wrz z-coord
        CP2.rightCols(sz2) = CP;
        const std::string cur = "sparsity" + util::to_string(j);
        gsWriteParaviewPoints(CP2, cur);

        col.addTimestep(cur, j, ".vtp");
    }
    col.save();
}


void localBasis_elwise(const gsMultiPatch<> & mp,
                       const gsMatrix<> & M,
                       const gsDofMapper & map,
                       gsMatrix<> & result)
{
    const index_t d = mp.parDim();//=2
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    const gsAffineFunction<> ifaceMap(mp.getMapForInterface(iFace));
    const index_t nc = map.boundarySize();

    const index_t l = iFace.first().direction(); // "orthogonal" direction
    gsMatrix<> pt1, gr, u2, tmp;
    static_cast<const gsTensorBSplineBasis<2>&>(B1).component(!l).knots().centers_into(gr);
    pt1.resize(d, gr.cols());// center points on the edge
    pt1.row(!l) = gr;
    pt1.row( l).setConstant( iFace.first().parameter() ? 1 : 0);

    gsMatrix<index_t> act1, act2;
    gsSortedVector<unsigned> C;
    index_t c = 0;

    result.setZero(nc, M.fullPivLu().dimensionOfKernel());
    EigenSLV slv;
    EigenSLV kslv;

    for (index_t i = 0; i != pt1.cols(); ++i) // for each center point
    {
        if (c == result.cols() ) return;

        const gsMatrix<> & u1 = pt1.col(i);

        ifaceMap.eval_into(u1, u2);

        B1.active_into(u1, act1);
        B2.active_into(u2, act2);
        map.localToGlobal(act1, iFace.first().patch , act1);
        map.localToGlobal(act2, iFace.second().patch, act2);

        C.clear();
        for ( index_t j = 0; j != act1.rows(); ++j) // ==act2.rows()
        {
            const index_t jj1 = act1.at(j);

            if ( map.is_boundary_index(jj1) )
            {
                const index_t bjj = map.global_to_bindex(jj1);
                C.push_sorted_unique(bjj);
            }
            const index_t jj2 = act2.at(j);
            if ( map.is_boundary_index(jj2) )
            {
                const index_t bjj = map.global_to_bindex(jj2);
                C.push_sorted_unique(bjj);
            }
        }

        // form submatrix
        M.submatrixCols( gsAsVector<unsigned>(C), tmp);
        kslv.compute(tmp);

        if ( 0 != kslv.dimensionOfKernel() ) // any new ones ?
        {
            tmp = kslv.kernel();

            if (0==c) // first one ?
            {
                for (index_t k=0; k!= tmp.cols(); ++k)
                    for (size_t t=0; t!= C.size(); ++t)
                        result(C[t],k) = tmp(t,k);
                c = tmp.cols();
            }
            else
            {
                for (index_t k=0; k!= tmp.cols(); ++k)
                {
                    result.col(c).setZero();
                    for (size_t t=0; t!= C.size(); ++t)
                        result(C[t],c) = tmp(t,k);

                    if (slv.compute(result.leftCols(c+1)).rank() == c+1 ) // l.i. ?
                        ++c; // keep this one
                }
            }
        }
    }

    if (c != result.cols() ) // in case we skipped any functions
    {
        gsInfo<<"Found "<< c <<", still needed "<< result.cols() - c <<"\n";

        tmp = kslv.compute(M).kernel();//full kernel
        slv.compute(result.leftCols(c));
        for (index_t k=0; k!= tmp.cols(); ++k)
        {
            result.col(c) = tmp.col(k);

            if (slv.compute(result.leftCols(c+1)).rank() == c+1 ) // l.i. ?
            {
                ++c;
                if (c == result.cols() ) return;
            }
        }
    }

    GISMO_ASSERT(c == result.cols(), "Still missing "<<result.cols()-c);
}

gsVector<index_t> asVector(const std::set<index_t> & s)
{
    const index_t sz = s.size();
    gsVector<index_t> v(sz);
    util::copy(s.begin(), s.end(), v.data());
    return v;
}


std::pair<index_t,index_t> computeDim(const gsMultiPatch<> & mp, bool randomize)
{
    index_t kDim;
    MatType G;
    gsDofMapper map;
    index_t pdeg = mb[0].maxDegree();

    // ----------
    std::vector<gsGeometry<>::Ptr> pert1;
    if (randomize)
        getPerturbation(pert1, mp.parDim() );

    // gsGluingData gd(mp, /*pdeg,*/ pert1);
    // G1Matrix(mp, pdeg, gd, G, map);

    //gsGluingData gd(mp, /*pdeg,*/ pert1);
    //gsSparseMatrix<> fb = FullBasisNonHom(mp, /*pdeg,*/ gd, map);
    //gsInfo<<"OK all done\n";exit(0);

    G1Matrix(mp, pdeg, G, map, randomize);

/*
    MatType Gq;
    G1Matrix(mp, pdeg, Gq, map, randomize);
    gsInfo<<"Diff\n"<< (G-Gq).cast<double>() <<"\n";
    gsDebugVar(G.cast<double>() );
    gsDebugVar(Gq.cast<double>());
    gsInfo<<"Diff norm: "<< (G-Gq).squaredNorm() <<"\n";
    G.swap(Gq);
*/

#ifdef GISMO_WITH_GMP
    G.rrefInPlace();
    kDim = rrefKerDim(G);
#else
#  if useSparse
    gsSparseSolver<>::QR slv;
    kDim = G.cols() - slv.compute(G).rank();
#  else
    kDim = G.fullPivLu().dimensionOfKernel();
#  endif
#endif
    return std::make_pair(kDim, map.freeSize());
}

size_t computeDimFormula(const gsMultiPatch<> & mp, bool randomize)
{
    const boundaryInterface & iFace = *mp.iBegin();

    index_t pdeg = mb[0].maxDegree();

    // Ensure the degree is at least 3
    gsMultiPatch<> _mp = mp;
    if (pdeg < 3) // p+2
    {
        _mp.degreeElevate(3-pdeg);
        pdeg = 3;
    }

    gsMatrix<index_t> kDim(3,3), rValues(3,3), fDim(3,3);
    gsMatrix<> val(2,9);
    gsDofMapper map;
    gsMultiPatch<> tmp;


    const index_t d = mp.parDim();
    std::vector<gsGeometry<>::Ptr> pert1;
    real_t tol = 10e-1; // todo: add as option if needed

    if (randomize)
        getPerturbation(pert1, d, tol);

    for (index_t p = 0; p!=3; ++p) //pdeg .. pdeg+p
    {
        for (index_t k = 0; k!=3; ++k) //--
        {
            gsInfo<<"p = "<<pdeg+p<<", k = "<<k<<"\n";
            tmp = _mp;
            tmp.degreeElevate(p);
            tmp.uniformRefine(k, pdeg+p-1);//mult=deg-1

            MatType G;

            //gsGluingData gd(tmp, /*pdeg,*/ pert1);
            //G1Matrix(tmp, pdeg, gd, G, map);
            G1Matrix(tmp, pdeg, pert1, G, map);

            #ifdef GISMO_WITH_GMP
            G.rrefInPlace();
            kDim(k,p)  = rrefKerDim(G);
#           else
            #if useSparse
            gsSparseSolver<>::QR slv;
            kDim(k,p)  = G.cols() - slv.compute(G).rank();
            #else
            kDim(k,p)  = G.fullPivLu().dimensionOfKernel();
            #endif
#           endif

            val(0,3*p+k) = pdeg+p;
            val(1,3*p+k) = k    ;
            // note: full dim is bi-cubic!
            //kDim(k,p) = G.fullPivLu().dimensionOfKernel() + map.freeSize();

            // rank values
            const index_t onInt = tmp.patch(tmp.iBegin()->first().patch).basis()
                                       .boundary(iFace.first().side()).size();
            rValues(k,p) = G.cols() + onInt - kDim(k,p);
            fDim   (k,p) = kDim(k,p) + map.freeSize();
        }
    }

    gsInfo<<"samples(k=0..2,p="<<pdeg<<".."<<pdeg+2<<"):\n"<< kDim <<"\n";
    gsInfo<<"ranks(k=0..2,p="<<pdeg<<".."<<pdeg+2<<"):\n"<< rValues <<"\n";
    gsInfo<<"final(k=0..2,p="<<pdeg<<".."<<pdeg+2<<"):\n"<< fDim <<"\n";
    gsMonomialBasis<real_t> dfu(2);
    gsGenericTensorBasis<2,real_t> df(dfu,dfu);
    gsSparseMatrix<> cm;
    df.collocationMatrix(val,cm);
    gsSparseSolver<>::LU slv;
    //gsMatrix<> cf = cm.toDense().fullPivLu().solve(kDim.asVector().cast<real_t>());
    gsMatrix<> cf = slv.compute(cm).solve(kDim.asVector().cast<real_t>());
    cf.resize(3,3);
#ifndef GISMO_WITH_GMP
    cf.removeNoise(10e-11);
#endif
    gsInfo<<"result(p,k):\n"<< cf <<"\n";

    return util::hash_range(kDim.data(), kDim.data()+9);
}

void sortByColWeight(gsMatrix<> & K, gsVector<index_t> & wgt)
{
    unsigned lastSwapDone = K.cols() - 1;
    unsigned lastCheckIdx = lastSwapDone;
    bool didSwap;

    gsVector<index_t> sz(K.cols()), prm = gsVector<index_t>::LinSpaced(K.cols(),0,lastSwapDone);
    for( index_t i=0; i != K.cols(); i++)
        sz[i] = (K.col(i).array()!=0).count();

    do
    {
        didSwap = false;
        lastCheckIdx = lastSwapDone;

        for( unsigned i=0; i < lastCheckIdx; i++)
            if( wgt[i] > wgt[i+1] )
            {
                std::swap(wgt[i] , wgt[i+1]);
                std::swap(prm[i], prm[i+1]);
                didSwap = true;
                lastSwapDone = i;
            }
    }
    while( didSwap );

    K =  K * prm.asPermutation();
}

void sortBySparsity(gsMatrix<> & K)
{
    gsVector<index_t> sz(K.cols()), prm = gsVector<index_t>::LinSpaced(K.cols(),0,K.cols() - 1);
    for( index_t i=0; i != K.cols(); i++)
        sz[i] = (K.col(i).array()!=0).count();
    sortByColWeight(K, sz);
}

void selectRepresentatives(gsMatrix<> & K, gsVector<index_t> & count)
{
    if ( 0==K.cols() ) return;

    sortBySparsity(K);

    gsMatrix<> Krep;
    gsVector<index_t> sz(K.cols());
    for( index_t i=0; i != K.cols(); i++)
        sz[i] = (K.col(i).array()!=0).count();

    std::vector<index_t> usz(sz.begin(), sz.end());
    std::vector<index_t>::iterator it = std::unique(usz.begin(), usz.end());
    usz.resize(std::distance(usz.begin(), it) );

    index_t m = 0;
    usz[m]    = 0;
    for (index_t j=1; j!=sz.size(); ++j)
        if (sz[j]!=sz[j-1]) usz[m++] = j;
    usz[m] = sz.size();

    count.resize(usz.size());
    count[0] = usz[0];
    for (size_t k=1; k!=usz.size(); ++k)
        count[k] = usz[k] - usz[k-1];

    gsAsVector<index_t>(usz).array() -= 1;
    K.submatrixCols(usz, Krep);
    Krep.swap(K);
    //sortByColWeight(K, count);
}


void order_boundary(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind,
    index_t s = 2,
    const bool trim = true)
{
    const index_t sz  = B1.component(0).size();
    s = std::min<index_t>(s, sz/2);
    const index_t numInt  =
        r1.size() - (B1.component(0).size()-2*s) * (B1.component(1).size()-2*s)
        - ( trim ? s * B1.component(0).size() : 0 );
    ind.resize( 3*numInt );
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz,sz);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz,sz);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz,sz);
    index_t m = 0;

    for(index_t i=s; i < sz; ++i)
        for(index_t k=0; k < s; ++k)
        {
            ind(m++) = map.bindex(rr2(i,k), iFace.first() .patch);
            ind(m++) = map.bindex(rr1(i,k), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,k), iFace.second().patch);
        }
    for(index_t j=s; j < sz; ++j)
        for(index_t k=0; k < s; ++k)
        {
            ind(m++) = map.bindex(rr2(sz-1-k,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr1(sz-1-k,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(sz-1-k,j), iFace.second().patch);
        }

    for(index_t i=sz-1-s; i >= (trim?s:0); --i)
        for(index_t k=0; k < s; ++k)
        {
            ind(m++) = map.bindex(rr2(i,sz-1-k), iFace.first() .patch);
            ind(m++) = map.bindex(rr1(i,sz-1-k), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,sz-1-k), iFace.second().patch);
        }
    if (!trim)
        for(index_t j=sz-1-s; j >= 0; --j)
            for(index_t k=0; k < s; ++k)
            {
                ind(m++) = map.bindex(rr2(k,j), iFace.first() .patch);
                ind(m++) = map.bindex(rr1(k,j), iFace.first() .patch);
                ind(m++) = map.bindex(rr3(k,j), iFace.second().patch);
            }
}

void order_patchwise(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind
    //,index_t s = 2
        )
{
    const index_t numInt = r2.size();
    ind.resize( 3*numInt );
    const index_t sz  = B1.component(0).size();
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz,sz);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz,sz);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz,sz);
    index_t m = 0;

    const index_t d = B1.component(0).maxDegree();
    const index_t k = ((index_t)std::sqrt(numInt)-(d+1)) / (d-1);
    //const index_t sqrtNumInt = (index_t)std::sqrt(numInt);
    for(index_t i=0; i<=d+k*(d-1); ++i)
    {
        m=((i+1)*(i+1)-1)*3;
        ind(m++) = map.bindex(rr1(i,i), iFace.first() .patch);
        ind(m++) = map.bindex(rr2(i,i), iFace.first() .patch);
        ind(m++) = map.bindex(rr3(i,i), iFace.second().patch);
    }

    for(index_t j=0; j<d+k*(d-1);++j)
    {
        for(index_t i=j+1; i<=d+k*(d-1);++i)
        {
            m=i*i*3+6*j;
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    }

    for(index_t i=0; i<d+k*(d-1);++i)
    {
        for(index_t j=i+1; j<=d+k*(d-1);++j)
        {
            m=(j*j+1)*3+6*i;
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }

    }

}
void order_pln(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind,
    const index_t st, const index_t e)
{
    const index_t sz = B1.component(0).size();
    index_t m = 0;
    ind.resize(3*(e-st));
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz,1);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz,1);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz,1);

    for(index_t i=st; i < e; ++i)
        {
            ind(m++) = map.bindex(rr1(i,0), iFace.first() .patch);
            ind(m++) = map.bindex(rr2(i,0), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,0), iFace.second().patch);
        }
}
void order_pln1_1(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind)
{
    const index_t sz = B1.component(0).size();
    index_t m = 0;
    ind.resize(4);

    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz,1);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz,1);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz,1);
            ind(m++) = map.bindex(rr1(0,0), iFace.first() .patch);
            ind(m++) = map.bindex(rr2(0,0), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(0,0), iFace.second().patch);
            ind(m++) = map.bindex(rr1(1,0), iFace.first() .patch);
}
void order_pln1_2(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind)
{
    const index_t sz = B1.component(0).size();
    index_t m = 0;
    ind.resize(4);
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz,1);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz,1);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz,1);
            ind(m++) = map.bindex(rr2(0,0), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(0,0), iFace.second().patch);
            ind(m++) = map.bindex(rr1(1,0), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(1,0), iFace.second().patch);

}
void order_pln1_3(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind)
{
    const index_t sz = B1.component(0).size();
    index_t m = 0;
    ind.resize(4);
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz,1);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz,1);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz,1);
            ind(m++) = map.bindex(rr1(sz-2,0), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(sz-2,0), iFace.second().patch);
            ind(m++) = map.bindex(rr1(sz-1,0), iFace.first() .patch);
            ind(m++) = map.bindex(rr2(sz-1,0), iFace.first() .patch);

}
void order_pln1_4(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind)
{
    const index_t sz = B1.component(0).size();
    index_t m = 0;
    ind.resize(4);
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz,1);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz,1);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz,1);
            ind(m++) = map.bindex(rr3(sz-2,0), iFace.second().patch);
            ind(m++) = map.bindex(rr1(sz-1,0), iFace.first() .patch);
            ind(m++) = map.bindex(rr2(sz-1,0), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(sz-1,0), iFace.second().patch);

}
void order_pln_if(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind,
    const index_t st, const index_t e)
{
    const index_t sz = B1.component(0).size();
    index_t m = 0;
    ind.resize(2*(e-st+1));
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz,1);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz,1);
    for(index_t i=st; i < e+1; ++i)
    {
        ind(m++) = map.bindex(rr1(i,0), iFace.first() .patch);
        ind(m++) = map.bindex(rr3(i,0), iFace.second().patch);
    }
}
void order_pln_oth(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind,
    const index_t st, const index_t e)
{
    const index_t sz = B1.component(0).size();
    index_t m = 0;
    ind.resize(3*(e-st-1)+2);
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz,1);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz,1);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz,1);

        ind(m++) = map.bindex(rr3(st,0), iFace.second() .patch);
    for(index_t i=st+1; i < e; ++i)
    {
        ind(m++) = map.bindex(rr1(i,0), iFace.first() .patch);
        ind(m++) = map.bindex(rr2(i,0), iFace.first() .patch);
        ind(m++) = map.bindex(rr3(i,0), iFace.second().patch);
    }
        ind(m++) = map.bindex(rr1(e,0), iFace.first().patch);

}

void order_tp(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t st1, const index_t e1,
    const index_t sz2, const index_t st2, const index_t e2)
{
    index_t m = 0;
    ind.resize(3 * (e1-st1) * (e2-st2));
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    for(index_t i=st1; i < e1; ++i)
        for(index_t j=st2; j < e2; ++j)
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);

    for(index_t i=st1; i < e1; ++i)
        for(index_t j=st2; j < e2; ++j)
        {
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
}
void order_ln1(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t sz2, const index_t st1, const index_t st2)
{
    index_t m = 0;
    ind.resize(2);
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    ind(m++) = map.bindex(rr1(st1,st2), iFace.first() .patch);
    ind(m++) = map.bindex(rr3(st1,st2), iFace.second().patch);
}

void order_ln2(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t st1, const index_t e1,
    const index_t sz2, const index_t st2, const index_t e2,
        const index_t stif1, const index_t eif1,
        const index_t stif2, const index_t eif2)
{
    index_t m = 0;
    ind.resize((e1-st1) * (e2-st2)-4 + (eif1-stif1) * (eif2-stif2));
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    for(index_t i=stif1; i < eif1; ++i)
        for(index_t j=stif2; j < eif2; ++j)
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);

    for(index_t i=st1+1; i < e1-1; ++i)
        for(index_t j=st2; j < st2+1; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=st1; i < e1; ++i)
        for(index_t j=st2+1; j < st2+4; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=st1; i < e1-1; ++i)
        for(index_t j=e2-2; j < e2-1; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second() .patch);
        }
    ind(m++) = map.bindex(rr1(e1-1,e2-2), iFace.first() .patch);

    for(index_t i=st1+1; i < e1-1; ++i)
        for(index_t j=e2-1; j < e2; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second() .patch);
        }
}
void order_ln3(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t st1, const index_t e1,
    const index_t sz2, const index_t st2, const index_t e2,
        const index_t stif1, const index_t eif1,
        const index_t stif2, const index_t eif2)
{
    index_t m = 0;
    ind.resize((e1-st1) * (e2-st2)-4 + (eif1-stif1) * (eif2-stif2));
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    for(index_t i=stif1; i < eif1; ++i)
        for(index_t j=stif2; j < eif2; ++j)
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);

    for(index_t i=st1+1; i < e1-1; ++i)
        for(index_t j=st2; j < st2+1; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=st1; i < e1; ++i)
        for(index_t j=st2+1; j < st2+2; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=st1; i < e1-1; ++i)
        for(index_t j=e2-2; j < e2-1; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second() .patch);
        }
    ind(m++) = map.bindex(rr1(e1-1,e2-2), iFace.first() .patch);

    for(index_t i=st1+1; i < e1-1; ++i)
        for(index_t j=e2-1; j < e2; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second() .patch);
        }
}
void order_ln4(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t st1, const index_t e1,
    const index_t sz2, const index_t st2, const index_t e2,
        const index_t stif1, const index_t eif1,
        const index_t stif2, const index_t eif2)
{
    index_t m = 0;
    ind.resize((e1-st1) * (e2-st2)-4 + (eif1-stif1) * (eif2-stif2));
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    for(index_t i=stif1; i < eif1; ++i)
        for(index_t j=stif2; j < eif2; ++j)
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);

    for(index_t i=st1+1; i < e1-1; ++i)
        for(index_t j=st2; j < st2+1; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=st1; i < e1; ++i)
        for(index_t j=st2+1; j < st2+3; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=st1; i < e1-1; ++i)
        for(index_t j=e2-2; j < e2-1; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second() .patch);
        }
    ind(m++) = map.bindex(rr1(e1-1,e2-2), iFace.first() .patch);

    for(index_t i=st1+1; i < e1-1; ++i)
        for(index_t j=e2-1; j < e2; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second() .patch);
        }
}
void order_ln5(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t st1, const index_t e1,
    const index_t sz2, const index_t st2, const index_t e2,
        const index_t stif1, const index_t eif1,
        const index_t stif2, const index_t eif2)
{
    index_t m = 0;
    ind.resize((e1-st1) * (e2-st2)-4 + (eif1-stif1) * (eif2-stif2));
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    for(index_t i=stif1; i < eif1; ++i)
        for(index_t j=stif2; j < eif2; ++j)
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);

    for(index_t i=st1+1; i < e1-1; ++i)
        for(index_t j=st2; j < st2+1; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=st1; i < e1; ++i)
        for(index_t j=st2+1; j < st2+2; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=st1; i < e1-1; ++i)
        for(index_t j=e2-2; j < e2-1; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second() .patch);
        }
    ind(m++) = map.bindex(rr1(e1-1,e2-2), iFace.first() .patch);

    for(index_t i=st1+1; i < e1-1; ++i)
        for(index_t j=e2-1; j < e2; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second() .patch);
        }
}
void order_ln6(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t st1, const index_t e1,
    const index_t sz2, const index_t st2, const index_t e2,
        const index_t stif1, const index_t eif1,
        const index_t stif2, const index_t eif2)
{
    index_t m = 0;
    ind.resize((e1-st1) * (e2-st2)-4 + (eif1-stif1) * (eif2-stif2));
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    for(index_t i=stif1; i < eif1; ++i)
        for(index_t j=stif2; j < eif2; ++j)
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);

    for(index_t i=st1+1; i < e1-1; ++i)
        for(index_t j=st2; j < st2+1; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=st1; i < e1; ++i)
        for(index_t j=st2+1; j < st2+3; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=st1; i < e1-1; ++i)
        for(index_t j=e2-2; j < e2-1; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second() .patch);
        }
    ind(m++) = map.bindex(rr1(e1-1,e2-2), iFace.first() .patch);

    for(index_t i=st1+1; i < e1-1; ++i)
        for(index_t j=e2-1; j < e2; ++j)
        {
            ind(m++) = map.bindex(rr3(i,j), iFace.second() .patch);
        }
}
void order_tp_if(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t st1, const index_t e1,
    const index_t sz2, const index_t st2, const index_t e2,
        const index_t stif1, const index_t eif1,
        const index_t stif2, const index_t eif2)
{
    index_t m = 0;
    ind.resize(2 * (e1-st1) * (e2-st2) + (eif1-stif1) * (eif2-stif2));
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    for(index_t i=stif1; i < eif1; ++i)
        for(index_t j=stif2; j < eif2; ++j)
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);

    for(index_t i=st1; i < e1; ++i)
        for(index_t j=st2; j < e2; ++j)
        {
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
}
void order_tp_if5(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t st1, const index_t e1,
    const index_t sz2, const index_t st2, const index_t e2,
        const index_t stif1, const index_t eif1,
        const index_t stif2, const index_t eif2)
{
    index_t m = 0;
    ind.resize(2 * (e1-st1) * (e2-st2) + (eif1-stif1) * (eif2-stif2)-2);
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    for(index_t i=stif1; i < eif1; ++i)
        for(index_t j=stif2; j < eif2; ++j)
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);

    for(index_t i=st1; i < e1-1; ++i)
        for(index_t j=st2; j < e2; ++j)
        {
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=e1-1; i < e1; ++i)
        for(index_t j=st2; j < e2-1; ++j)
        {
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
}
void order_tp_if6(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t st1, const index_t e1,
    const index_t sz2, const index_t st2, const index_t e2,
        const index_t stif1, const index_t eif1,
        const index_t stif2, const index_t eif2)
{
    index_t m = 0;
    ind.resize(2 * (e1-st1) * (e2-st2) + (eif1-stif1) * (eif2-stif2)-3);
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    for(index_t i=stif1; i < eif1; ++i)
        for(index_t j=stif2; j < eif2; ++j)
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);

    for(index_t i=st1; i < e1-1; ++i)
        for(index_t j=st2; j < e2; ++j)
        {
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=e1-1; i < e1; ++i)
        for(index_t j=st2; j < e2-2; ++j)
        {
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=e1-1; i < e1; ++i)
        for(index_t j=e2-2; j < e2-1; ++j)
        {
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
        }
}

void order_tp_if6p(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    gsVector<index_t> & ind,
    const index_t sz1, const index_t st1, const index_t e1,
    const index_t sz2, const index_t st2, const index_t e2,
        const index_t stif1, const index_t eif1,
        const index_t stif2, const index_t eif2)
{
    index_t m = 0;
    ind.resize(2 * (e1-st1) * (e2-st2) + (eif1-stif1) * (eif2-stif2)-4);
    gsAsConstMatrix<index_t> rr1 = r1.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr2 = r2.reshape(sz1,sz2);
    gsAsConstMatrix<index_t> rr3 = r3.reshape(sz1,sz2);

    for(index_t i=stif1; i < eif1; ++i)
        for(index_t j=stif2; j < eif2; ++j)
            ind(m++) = map.bindex(rr2(i,j), iFace.first() .patch);

    for(index_t i=st1; i < e1-1; ++i)
        for(index_t j=st2; j < e2; ++j)
        {
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=e1-1; i < e1; ++i)
        for(index_t j=st2; j < e2-3; ++j)
        {
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
            ind(m++) = map.bindex(rr3(i,j), iFace.second().patch);
        }
    for(index_t i=e1-1; i < e1; ++i)
        for(index_t j=e2-3; j < e2-1; ++j)
        {
            ind(m++) = map.bindex(rr1(i,j), iFace.first() .patch);
        }
}

void order_cropped(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind,
    const index_t s = 0)
{
    if (2 == B1.dim())
    {
        const index_t sz1 = B1.component(0).size();
        order_tp(iFace, r1, r2, r3, map, ind, sz1, s, sz1-s, 1, 0, 1);
    }
    else //3
    {
        const index_t sz1 = B1.component(0).size();
        const index_t sz2 = B1.component(1).size();
        order_tp(iFace, r1, r2, r3, map, ind, sz1, s, sz1-s, sz2, s, sz2-s);
    }
}

void order_trimmed(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind,
    const index_t s = 2)
{
    if (2 == B1.dim())
    {
        const index_t sz1 = B1.component(0).size();
        order_tp(iFace, r1, r2, r3, map, ind, sz1, 0, sz1-s, 1, 0, 1);
    }
    else //3
    {
        const index_t sz1 = B1.component(0).size();
        const index_t sz2 = B1.component(1).size();
        order_tp(iFace, r1, r2, r3, map, ind, sz1, 0, sz1, sz2, 0, sz2-s);
    }
}

void order_corner(
    const boundaryInterface & iFace,
    const gsMatrix<index_t> & r1,
    const gsMatrix<index_t> & r2,
    const gsMatrix<index_t> & r3,
    const gsDofMapper & map,
    const gsBasis<> & B1,
    gsVector<index_t> & ind,
    const index_t s)
{
    if (2 == B1.dim())
    {
        const index_t sz1 = B1.component(0).size();
        order_tp(iFace, r1, r2, r3, map, ind, sz1, 0, s, 1, 0, 1);
    }
    else //3
    {
        const index_t sz1 = B1.component(0).size();
        const index_t sz2 = B1.component(1).size();
        order_tp(iFace, r1, r2, r3, map, ind, sz1, 0, s, sz2, 0, s);
    }
}

void append_ind(gsVector<index_t> & ind, gsVector<index_t> & indNew)
{
    ind.conservativeResize(ind.size()+indNew.size());
    ind.bottomRows(indNew.size()) = indNew;
}

template<int d>
void localBasis_incr(const gsMultiPatch<> & mp,
                     const index_t pdeg,
                     gsDofMapper & map,
                     gsMatrix<> & result, bool bdr)
{
    MatType M;
    G1Matrix(mp, pdeg, M, map, false);// matrix and mapper

    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb.basis( iFace.first() .patch );
    const gsBasis<> & B2 = mb.basis( iFace.second().patch );

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;
    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    if ( !iFace.dirOrientation(iFace.first(), !iFace.first().direction()) )
    {
        gsInfo<< "Reversing (2D).\n";
        std::reverse(r3.data(), r3.data()+r3.size());
    }
    if ( 2 == mp.parDim() )
        order_pln(iFace, r1, r2, r3, map, B1, ind, 0, B1.component(0).size());
    else //3D
    {
//      order_trimmed (iFace, r1, r2, r3, map, B1, ind, 2);
//      order_boundary(iFace, r1, r2, r3, map, B1, ind, 3, true);
        if (bdr)
            order_cropped(iFace, r1, r2, r3, map, B1, ind, 0);
//          order_corner (iFace, r1, r2, r3, map, B1, ind, B1.component(0).size()-2);
        else
            order_cropped(iFace, r1, r2, r3, map, B1, ind, 2);
     }
//    order_patchwise(iFace, r1, r2, r3, map, B1, ind);
    gsMatrix<> Mperm, Ker;

    Mperm = columns(M,ind);

/*
    gsInfo << "Matrix ("<<M.dim()<<", nz="
#if useSparse
           <<double((M.array() != 0).count())/M.size()<<")\n";
#else
           <<double(M.nonZeros())/M.size()<<")\n";
#endif
*/

#ifdef GISMO_WITH_GMP
    Mperm.rrefInPlace();
    (void)rrefKerBasis(Mperm, Ker);
    Ker.transpose().gaussElim();
    // if (Ker.cols() < 100 )
    //     gsInfo<< "K=\n"<< printSparsity(Ker.transpose()) <<"\n";
#else
    #if !useSparse
    Ker = Mperm.fullPivLu().kernel();
#endif
    Ker.removeNoise(10e-11);
#endif
    result.setZero(3*r1.size(), Ker.cols() );
    for (index_t m=0; m!= Ker.rows(); ++m)
        result.row(ind[m]) = Ker.row(m);
}


template<int d>
void localBasis_incr(const gsMultiPatch<> & mp,
                     const index_t pdeg,
                     const std::vector<gsGeometry<>::Ptr> & pert,
                     gsDofMapper & map,
                     gsMatrix<> & result, bool bdr)
{
    MatType M;
    gsGluingData gd(mp, /*pdeg,*/ pert);
    G1Matrix(mp, pdeg, gd, M, map);
   // G1Matrix(mp, pdeg, pert, M, map);// matrix and mapper
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb.basis( iFace.first() .patch );
    const gsBasis<> & B2 = mb.basis( iFace.second().patch );

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;
    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    if ( !iFace.dirOrientation(iFace.first(), !iFace.first().direction()) )
    {
        gsInfo<< "Reversing (2D).\n";
        std::reverse(r3.data(), r3.data()+r3.size());
    }

//    order_trimmed (iFace, r1, r2, r3, map, B1, ind, 2);
    //order_boundary(iFace, r1, r2, r3, map, B1, ind, 3, true);
    if (bdr)
        order_corner (iFace, r1, r2, r3, map, B1, ind, B1.component(0).size()-2);
    else
        order_cropped(iFace, r1, r2, r3, map, B1, ind, 2);

//    order_patchwise(iFace, r1, r2, r3, map, B1, ind);
    gsMatrix<> Mperm, Ker;
    Mperm=columns(M,ind);

    gsInfo << "Matrix ("<<M.dim()<<", nz="
           <<double(M.nonZeros())/M.size()<<")\n";

#ifdef GISMO_WITH_GMP
    Mperm.rrefInPlace();
    (void)rrefKerBasis(Mperm, Ker);
    Ker.transpose().gaussElim();
    // if (Ker.cols() < 100 )
    //     gsInfo<< "K=\n"<< printSparsity(Ker.transpose()) <<"\n";
#else
    #if !useSparse
    Ker = Mperm.fullPivLu().kernel();
#endif
    Ker.removeNoise(10e-11);
#endif
    result.setZero(3*r1.size(), Ker.cols() );
    for (index_t m=0; m!= Ker.rows(); ++m)
        result.row(ind[m]) = Ker.row(m);
}
void basisDegree3(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result)
{
    GISMO_ASSERT(mb[0].size()>7, "Problem, need more knots per direction");

    MatType M;
    gsDofMapper map;
    G1Matrix(mp, pdeg, M, map, false);// matrix and mapper

    //const index_t d = mp.parDim();//=2
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    const index_t sz1 = B1.component(0).size();
    const index_t sz2 = B1.component(1).size();
    const index_t sP1 = B1.size();
    const index_t sP2 = B2.size();

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;

    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    const boundaryInterface & iFacePlot = *mp.iBegin();
    gsVector<index_t> p(2);
    p(0) = iFacePlot.first() .patch;
    p(1) = iFacePlot.second().patch;

    GISMO_ENSURE(3==B1.component(0).degree(0), "Degree is not 3");
    const index_t nKnots0 = (B1.component(0).size() - B1.component(0).degree(0)-1) / (B1.component(0).degree(0)-1);
    const index_t nKnots1 = (B1.component(1).size() - B1.component(1).degree(1)-1) / (B1.component(1).degree(1)-1);
    const index_t nKnots2 = (B1.component(2).size() - B1.component(2).degree(2)-1) / (B1.component(2).degree(2)-1);

    gsMatrix<> subm, subk;
    index_t c = 0;
    const index_t dgr=B1.component(0).degree(0);

    //result = gsSparseMatrix<>(sP1+sP2, (nKnots-2)*(nKnots-2) );
    result = gsSparseMatrix<>(sP1+sP2, (nKnots0-2)*(nKnots0-2) + 2*(dgr+1+nKnots0*(dgr-1)-4)*(dgr+1+nKnots1*(dgr-1)-4)*(dgr+1+nKnots2*(dgr-1)-4));
    std::vector<std::pair<index_t,index_t> > pi;

    for (index_t i = 0; i != nKnots0-2; ++i)
    {
        for (index_t j = 0; j != nKnots0-2; ++j)
        {
            order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 2+i*2, 8+i*2, sz2, 2+j*2, 8+j*2, 3+i*2, 7+i*2, 3+j*2, 7+j*2);

            subm = columns(M, ind);
            subm.rrefInPlace();
            (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

            for (index_t s = 0; s!=ind.size(); ++s)
            {
                map.preImage(ind[s]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(s, 0);
                }
            }

            c++;// next function
        }
    }
    for (index_t i = 2; i!= dgr+nKnots2*(dgr-1)-1; ++i)
    {
        gsMatrix<index_t> r4,r5;
        r4 = B1.boundaryOffset(iFace.first() .side(), i);
        r5 = B2.boundaryOffset(iFace.second() .side(), i);

        const index_t sz = B1.component(0).size();
        ind.resize(1);
        gsAsConstMatrix<index_t> rr4 = r4.reshape(sz,sz);
        gsAsConstMatrix<index_t> rr5 = r5.reshape(sz,sz);
        for(index_t j=2; j!=sz-2; ++j)
            for(index_t k=2; k!=sz-2; ++k)
            {
                ind(0) = map.bindex(rr4(j,k), iFace.first() .patch);
                map.preImage(ind[0]+ map.freeSize(), pi);
                if (pi[0].first == p(0))
                    result.insert(sP1+pi[0].second, c)
                          = 1;
                else
                    result.insert(pi[0].second, c)
                          = 1;
                c++;

                ind(0) = map.bindex(rr5(j,k), iFace.second() .patch);
                map.preImage(ind[0]+ map.freeSize(), pi);
                if (pi[0].first == p(0))
                    result.insert(sP1+pi[0].second, c)
                              = 1;
                else
                    result.insert(pi[0].second, c)
                              = 1;
                c++;
            }
     }
}
//Generic case degree 4 with arbitrary boundary surfaces
/*
void basisDegree4(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result)
{
    MatType M;
    gsDofMapper map;
    G1Matrix(mp, pdeg, M, map, false);// matrix and mapper

    //const index_t d = mp.parDim();//=2
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    const index_t sz1 = B1.component(0).size();
    const index_t sz2 = B1.component(1).size();
    const index_t sP1 = B1.size();
    const index_t sP2 = B2.size();

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;
    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    const boundaryInterface & iFacePlot = *mp.iBegin();
    gsVector<index_t> p(2);
    p(0) = iFacePlot.first() .patch;
    p(1) = iFacePlot.second().patch;

    GISMO_ENSURE(4==B1.component(0).degree(0), "Degree is not 4");
    const index_t nKnots0 = (B1.component(0).size() - B1.component(0).degree(0)-1) / (B1.component(0).degree(0)-1);
    const index_t nKnots1 = (B1.component(1).size() - B1.component(1).degree(1)-1) / (B1.component(1).degree(1)-1);
    const index_t nKnots2 = (B1.component(2).size() - B1.component(2).degree(2)-1) / (B1.component(2).degree(2)-1);

    gsMatrix<> subm, subk;
    index_t c = 0;
    const index_t dgr=B1.component(0).degree(0);

    //result = gsSparseMatrix<>(sP1+sP2, (5*(nKnots*nKnots))-(6*nKnots)+2 );
    //result = gsSparseMatrix<>(sP1+sP2, (5*(nKnots*nKnots))-(6*nKnots)+2 + 2*(dgr+1+nKnots*(dgr-1)-4)*(dgr+1+nKnots*(dgr-1)-4)*(dgr+1+nKnots*(dgr-1)-4));
    result = gsSparseMatrix<>(sP1+sP2, (5*(nKnots0*nKnots0))-(6*nKnots0)+2 + 2*(dgr+1+nKnots0*(dgr-1)-4)*(dgr+1+nKnots1*(dgr-1)-4)*(dgr+1+nKnots2*(dgr-1)-4));

    std::vector<std::pair<index_t,index_t> > pi;
    //Type 1

    for (index_t i = 0; i != 1; ++i)
    {
        for (index_t j = 0; j != nKnots1; ++j)
        {
            order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 6+i*3, sz2, 2+j*3, 6+j*3, 3+i*3, 5+i*3, 3+j*3, 5+j*3);
            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

            for (index_t i = 0; i!=ind.size(); ++i)
            {
                map.preImage(ind[i]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(i, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(i, 0);
                }
            }

            c++;// next function
        }
    }
    for (index_t j = 0; j != 1; ++j)
    {
        for (index_t i = 1; i != nKnots1; ++i)
        {
            order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 6+i*3, sz2, 2+j*3, 6+j*3, 3+i*3, 5+i*3, 3+j*3, 5+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

            for (index_t i = 0; i!=ind.size(); ++i)
            {
                map.preImage(ind[i]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(i, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(i, 0);
                }
            }

            c++;// next function
        }
    }

    if(nKnots0>1)
        //Type 2.1
    {
        for (index_t i = 0; i != 1; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 6+i*3, sz2, 3+j*3, 8+j*3, 3+i*3, 5+i*3, 4+j*3,7+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t i = 0; i!=ind.size(); ++i)
                {
                    map.preImage(ind[i]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(i, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(i, 0);
                    }
                }

                c++;// next function
            }
        }

        //Type 2.2
        for (index_t i = 0; i != nKnots0-1; ++i)
        {
            for (index_t j = 0; j != nKnots1; ++j)
            {
                order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 3+i*3, 8+i*3, sz2, 2+j*3, 6+j*3, 4+i*3, 7+i*3, 3+j*3, 5+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t i = 0; i!=ind.size(); ++i)
                {
                    map.preImage(ind[i]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(i, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(i, 0);
                    }
                }

                c++;// next function
            }
        }
        //Type 3
        for (index_t i = 0; i != nKnots0-1; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 3+i*3, 8+i*3, sz2, 3+j*3, 8+j*3, 4+i*3, 7+i*3, 4+j*3, 7+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t i = 0; i!=ind.size(); ++i)
                {
                    map.preImage(ind[i]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(i, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(i, 0);
                    }
                }

                c++;// next function
            }
        }
        //Type 4
        for (index_t i = 0; i != nKnots0-1; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 9+i*3, sz2, 2+j*3, 9+j*3, 0, 0, 0, 0);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t i = 0; i!=ind.size(); ++i)
                {
                    map.preImage(ind[i]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(i, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(i, 0);
                    }
                }

                c++;// next function
            }
        }
        //Type 5
        for (index_t i = 0; i != nKnots0-1; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if5(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 9+i*3, sz2, 2+j*3, 9+j*3, 6+i*3, 8+i*3, 6+j*3, 8+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t i = 0; i!=ind.size(); ++i)
                {
                    map.preImage(ind[i]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(i, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(i, 0);
                    }
                }

                c++;// next function
            }
        }
        //Type 6
        for (index_t i = 0; i != nKnots0-1; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if6(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 9+i*3, sz2, 2+j*3, 9+j*3, 6+i*3, 8+i*3, 4+j*3, 8+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t i = 0; i!=ind.size(); ++i)
                {
                    map.preImage(ind[i]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(i, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(i, 0);
                    }
                }

                c++;// next function
            }
        }
    }
    for (index_t i = 2; i!= dgr+nKnots2*(dgr-1)-1; ++i)
    {
        gsMatrix<index_t> r4,r5;
        r4 = B1.boundaryOffset(iFace.first() .side(), i);
        r5 = B2.boundaryOffset(iFace.second() .side(), i);

        const index_t sz = B1.component(0).size();
        ind.resize(1);
        gsAsConstMatrix<index_t> rr4 = r4.reshape(sz,sz);
        gsAsConstMatrix<index_t> rr5 = r5.reshape(sz,sz);
        for(index_t j=2; j!=sz-2; ++j)
            for(index_t k=2; k!=sz-2; ++k)
            {
                ind(0) = map.bindex(rr4(j,k), iFace.first() .patch);
                map.preImage(ind[0]+ map.freeSize(), pi);
                if (pi[0].first == p(0))
                    result.insert(sP1+pi[0].second, c)
                          = 1;
                else
                    result.insert(pi[0].second, c)
                          = 1;
                c++;

                ind(0) = map.bindex(rr5(j,k), iFace.second() .patch);
                map.preImage(ind[0]+ map.freeSize(), pi);
                if (pi[0].first == p(0))
                    result.insert(sP1+pi[0].second, c)
                              = 1;
                else
                    result.insert(pi[0].second, c)
                              = 1;
                c++;
            }
     }
}*/

//Generic case degree 4 with planar boundary surfaces but non-planar interface
void basisDegree4(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result)
{
    GISMO_ASSERT(mb[0].size()>10, "Problem, need more knots per direction");

    MatType M;
    gsDofMapper map;
    G1Matrix(mp, pdeg, M, map, false);// matrix and mapper

    //const index_t d = mp.parDim();//=2
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    const index_t sz1 = B1.component(0).size();
    const index_t sz2 = B1.component(1).size();
    const index_t sP1 = B1.size();
    const index_t sP2 = B2.size();

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;
    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    const boundaryInterface & iFacePlot = *mp.iBegin();
    gsVector<index_t> p(2);
    p(0) = iFacePlot.first() .patch;
    p(1) = iFacePlot.second().patch;

    GISMO_ENSURE(4==B1.component(0).degree(0), "Degree is not 4");
    const index_t nKnots0 = (B1.component(0).size() - B1.component(0).degree(0)-1) / (B1.component(0).degree(0)-1);
    const index_t nKnots1 = (B1.component(1).size() - B1.component(1).degree(1)-1) / (B1.component(1).degree(1)-1);
    const index_t nKnots2 = (B1.component(2).size() - B1.component(2).degree(2)-1) / (B1.component(2).degree(2)-1);

    gsMatrix<> subm, subk;
    index_t c = 0;
    const index_t dgr=B1.component(0).degree(0);

    //result = gsSparseMatrix<>(sP1+sP2, (5*(nKnots*nKnots))-(6*nKnots)+2 );
    //result = gsSparseMatrix<>(sP1+sP2, (5*(nKnots*nKnots))-(6*nKnots)+2 + 2*(dgr+1+nKnots*(dgr-1)-4)*(dgr+1+nKnots*(dgr-1)-4)*(dgr+1+nKnots*(dgr-1)-4));
    result = gsSparseMatrix<>(sP1+sP2, (5*(nKnots0*nKnots0))-(6*nKnots0)+2 + 2*(dgr+1+nKnots0*(dgr-1)-4)*(dgr+1+nKnots1*(dgr-1)-4)*(dgr+1+nKnots2*(dgr-1)-4));


    std::vector<std::pair<index_t,index_t> > pi;
    //Type 1

    for (index_t i = 0; i != 1; ++i)
    {
        for (index_t j = 0; j != nKnots1; ++j)
        {
            order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 6+i*3, sz2, 2+j*3, 6+j*3, 3+i*3, 5+i*3, 3+j*3, 5+j*3);
            // form submatrix
            subm = columns(M, ind);
            subm.rrefInPlace();
            (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

            for (index_t s = 0; s!=ind.size(); ++s)
            {
                map.preImage(ind[s]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(s, 0);
                }
            }

            c++;// next function
        }
    }
    for (index_t j = 0; j != 1; ++j)
    {
        for (index_t i = 1; i != nKnots1; ++i)
        {
            order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 6+i*3, sz2, 2+j*3, 6+j*3, 3+i*3, 5+i*3, 3+j*3, 5+j*3);

            // form submatrix
            subm = columns(M, ind);
            subm.rrefInPlace();
            (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

            for (index_t s = 0; s!=ind.size(); ++s)
            {
                map.preImage(ind[s]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(s, 0);
                }
            }

            c++;// next function
        }
    }

    if(nKnots1>1)
        //Type 2.1
    {
        for (index_t i = 0; i != 1; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 6+i*3, sz2, 3+j*3, 8+j*3, 3+i*3, 5+i*3, 4+j*3, 7+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t s = 0; s!=ind.size(); ++s)
                {
                    map.preImage(ind[s]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(s, 0);
                    }
                }

                c++;// next function
            }
        }

        //Type 2.2
        for (index_t i = 0; i != nKnots0-1; ++i)
        {
            for (index_t j = 0; j != nKnots1; ++j)
            {
                order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 3+i*3, 8+i*3, sz2, 2+j*3, 6+j*3, 4+i*3, 7+i*3, 3+j*3, 5+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t s = 0; s!=ind.size(); ++s)
                {
                    map.preImage(ind[s]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(s, 0);
                    }
                }

                c++;// next function
            }
        }
        //Type 3
        for (index_t i = 0; i != nKnots0-1; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 3+i*3, 8+i*3, sz2, 3+j*3, 8+j*3, 4+i*3, 7+i*3, 4+j*3, 7+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t s = 0; s!=ind.size(); ++s)
                {
                    map.preImage(ind[s]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(s, 0);
                    }
                }

                c++;// next function
            }
        }
        //Type 4
        for (index_t i = 0; i != nKnots0-1; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 9+i*3, sz2, 2+j*3, 9+j*3, 0, 0, 0, 0);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t s = 0; s!=ind.size(); ++s)
                {
                    map.preImage(ind[s]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(s, 0);
                    }
                }

                c++;// next function
            }
        }
        //Type 5
        for (index_t i = 0; i != nKnots0-1; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if5(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 9+i*3, sz2, 2+j*3, 9+j*3, 6+i*3, 8+i*3, 6+j*3, 8+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t s = 0; s!=ind.size(); ++s)
                {
                    map.preImage(ind[s]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(s, 0);
                    }
                }

                c++;// next function
            }
        }

        //Type 6.1
        for (index_t i = 0; i != 1; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if6p(iFace, r1, r2, r3, map, ind, sz1, 2+(nKnots1-2)*3, 9+(nKnots1-2)*3, sz2, 2+j*3, 9+j*3, 6+(nKnots1-2)*3, 8+(nKnots1-2)*3, 4+j*3, 8+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t s = 0; s!=ind.size(); ++s)
                {
                    map.preImage(ind[s]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(s, 0);
                    }
                }

                c++;// next function
            }
        }

        //Type 6.2
        for (index_t i = 0; i != nKnots0-2; ++i)
        {
            for (index_t j = 0; j != nKnots1-1; ++j)
            {
                order_tp_if6(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 9+i*3, sz2, 2+j*3, 9+j*3, 6+i*3, 8+i*3, 4+j*3, 8+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t  s = 0; s!=ind.size(); ++s)
                {
                    map.preImage(ind[s]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(s, 0);
                    }
                }

                c++;// next function
            }
        }
    }
    for (index_t i = 2; i!= dgr+nKnots2*(dgr-1)-1; ++i)
    {
        gsMatrix<index_t> r4,r5;
        r4 = B1.boundaryOffset(iFace.first() .side(), i);
        r5 = B2.boundaryOffset(iFace.second() .side(), i);

        const index_t sz = B1.component(0).size();
        ind.resize(1);
        gsAsConstMatrix<index_t> rr4 = r4.reshape(sz,sz);
        gsAsConstMatrix<index_t> rr5 = r5.reshape(sz,sz);
        for(index_t j=2; j!=sz-2; ++j)
            for(index_t k=2; k!=sz-2; ++k)
            {
                ind(0) = map.bindex(rr4(j,k), iFace.first() .patch);
                map.preImage(ind[0]+ map.freeSize(), pi);
                if (pi[0].first == p(0))
                    result.insert(sP1+pi[0].second, c)
                          = 1;
                else
                    result.insert(pi[0].second, c)
                          = 1;
                c++;

                ind(0) = map.bindex(rr5(j,k), iFace.second() .patch);
                map.preImage(ind[0]+ map.freeSize(), pi);
                if (pi[0].first == p(0))
                    result.insert(sP1+pi[0].second, c)
                              = 1;
                else
                    result.insert(pi[0].second, c)
                              = 1;
                c++;
            }
     }
}
void basisPlnDegree3(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result)
{
    MatType M;
    gsDofMapper map;
    G1Matrix(mp, pdeg, M, map, false);// matrix and mapper

    //const index_t d = mp.parDim();//=2
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    const index_t sP1 = B1.size();
    const index_t sP2 = B2.size();
    const index_t sz = B1.component(0).size();

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;
    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    const boundaryInterface & iFacePlot = *mp.iBegin();
    gsVector<index_t> p(2);
    p(0) = iFacePlot.first() .patch;
    p(1) = iFacePlot.second().patch;

    GISMO_ENSURE(3==B1.component(0).degree(0), "Degree is not 3");
    const index_t nKnots = (B1.component(1).size() - B1.component(0).degree(0)-1) /
        (B1.component(0).degree(0)-1);

    gsMatrix<> subm, subk;
    index_t c = 0;
    const index_t dgr=B1.component(0).degree(0);

    //result = gsSparseMatrix<>(sP1+sP2, 7+2*nKnots );
    result = gsSparseMatrix<>(sP1+sP2, 7+2*nKnots + 2*(dgr+1+nKnots*(dgr-1)-2)*(dgr+1+nKnots*(dgr-1)));

    std::vector<std::pair<index_t,index_t> > pi;

    //Type 1.1
    order_pln1_1(iFace, r1, r2, r3, map, B1, ind);
    // form submatrix
    subm = columns(M, ind);
    subm.rrefInPlace();
    (void)rrefKerBasis(subm, subk);

     GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

     for (index_t i = 0; i!=ind.size(); ++i)
     {
         map.preImage(ind[i]+ map.freeSize(), pi);
         for (size_t k = 0; k!= pi.size(); ++k)
         {
            if (pi[k].first == p(0)){

               result.insert(sP1+pi[k].second, c)
                         = subk(i, 0);}
            else
            {
               result.insert(pi[k].second, c)
                         = subk(i, 0);
            }
         }
     }
     c++;// next function


     //Type 1.2
     order_pln1_2(iFace, r1, r2, r3, map, B1, ind);
     // form submatrix
     subm = columns(M, ind);
     subm.rrefInPlace();
     (void)rrefKerBasis(subm, subk);

     GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

     for (index_t i = 0; i!=ind.size(); ++i)
     {
         map.preImage(ind[i]+ map.freeSize(), pi);
         for (size_t k = 0; k!= pi.size(); ++k)
         {
            if (pi[k].first == p(0))
             result.insert(sP1+pi[k].second, c)
                       = subk(i, 0);
            else
             result.insert(pi[k].second, c)
                       = subk(i, 0);
         }
     }
     c++;// next function

   //Type 1.3
     order_pln1_3(iFace, r1, r2, r3, map, B1, ind);
   // form submatrix
     subm = columns(M, ind);
     subm.rrefInPlace();
     (void)rrefKerBasis(subm, subk);

     GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

    for (index_t i = 0; i!=ind.size(); ++i)
     {
         map.preImage(ind[i]+ map.freeSize(), pi);
         for (size_t k = 0; k!= pi.size(); ++k)
         {
            if (pi[k].first == p(0))
             result.insert(sP1+pi[k].second, c)
                       = subk(i, 0);
            else
             result.insert(pi[k].second, c)
                       = subk(i, 0);
         }
     }
     c++;// next function

     //Type 1.4
       order_pln1_4(iFace, r1, r2, r3, map, B1, ind);
     // form submatrix
       subm = columns(M, ind);
       subm.rrefInPlace();
       (void)rrefKerBasis(subm, subk);

       GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

      for (index_t i = 0; i!=ind.size(); ++i)
       {
           map.preImage(ind[i]+ map.freeSize(), pi);
           for (size_t k = 0; k!= pi.size(); ++k)
           {
              if (pi[k].first == p(0))
               result.insert(sP1+pi[k].second, c)
                         = subk(i, 0);
              else
               result.insert(pi[k].second, c)
                         = subk(i, 0);
           }
       }
       c++;// next function

    //Type 2
    for (index_t i = 0; i != 2; ++i)
    {
        order_pln_if(iFace, r1,  r3, map, B1, ind, 1+i*(nKnots*2-1), 3+i*(nKnots*2-1));

        // form submatrix
        subm = columns(M, ind);
        subm.rrefInPlace();
        (void)rrefKerBasis(subm, subk);

        GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

        for (index_t s = 0; s!=ind.size(); ++s)
        {
            map.preImage(ind[s]+ map.freeSize(), pi);
            for (size_t k = 0; k!= pi.size(); ++k)
            {
                if (pi[k].first == p(0))
                  result.insert(sP1+pi[k].second, c)
                            = subk(s, 0);
                else
                  result.insert(pi[k].second, c)
                            = subk(s, 0);
            }
        }

        c++;// next function
    }

    //Type 3.1
      order_pln_oth(iFace, r1, r2, r3, map, B1, ind, 0, 3);
    // form submatrix
      subm = columns(M, ind);
      subm.rrefInPlace();
      (void)rrefKerBasis(subm, subk);

      GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

     for (index_t s = 0; s!=ind.size(); ++s)
      {
          map.preImage(ind[s]+ map.freeSize(), pi);
          for (size_t k = 0; k!= pi.size(); ++k)
          {
             if (pi[k].first == p(0))
              result.insert(sP1+pi[k].second, c)
                        = subk(s, 0);
             else
              result.insert(pi[k].second, c)
                        = subk(s, 0);
          }
      }
      c++;// next function

      //Type 3.2
      for (index_t i = 0; i != nKnots-1; ++i)
      {
          order_pln_if(iFace, r1,  r3, map, B1, ind, 2+i*2, 5+i*2);

          // form submatrix
          subm = columns(M, ind);
          subm.rrefInPlace();
          (void)rrefKerBasis(subm, subk);

          GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

          for (index_t s = 0; s!=ind.size(); ++s)
          {
              map.preImage(ind[s]+ map.freeSize(), pi);
              for (size_t k = 0; k!= pi.size(); ++k)
              {
                  if (pi[k].first == p(0))
                    result.insert(sP1+pi[k].second, c)
                              = subk(s, 0);
                  else
                    result.insert(pi[k].second, c)
                              = subk(s, 0);
              }
          }

          c++;// next function
      }

      //Type 3.1
        order_pln_oth(iFace, r1, r2, r3, map, B1, ind, sz-4, sz-1);
      // form submatrix
        subm = columns(M, ind);
        subm.rrefInPlace();
        (void)rrefKerBasis(subm, subk);

        GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

       for (index_t i = 0; i!=ind.size(); ++i)
        {
            map.preImage(ind[i]+ map.freeSize(), pi);
            for (size_t k = 0; k!= pi.size(); ++k)
            {
               if (pi[k].first == p(0))
                result.insert(sP1+pi[k].second, c)
                          = subk(i, 0);
               else
                result.insert(pi[k].second, c)
                          = subk(i, 0);
            }
        }
      c++;// next function

        //Type 4
      for (index_t i = 0; i != 2; ++i)
        {
            order_pln_oth(iFace, r1, r2, r3, map, B1, ind, 1+i*(2*nKnots-3), 5+i*(2*nKnots-3));

            // form submatrix
            subm = columns(M, ind);
            subm.rrefInPlace();
            (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

            for (index_t s = 0; s!=ind.size(); ++s)
            {
                map.preImage(ind[s]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(s, 0);
                }
            }

            c++;// next function
      }

      //Type 5
      for (index_t i = 0; i != nKnots-2; ++i)
      {
          order_pln_oth(iFace, r1, r2, r3, map, B1, ind, 2+i*2, 7+i*2);

          // form submatrix
          subm = columns(M, ind);
          subm.rrefInPlace();
          (void)rrefKerBasis(subm, subk);

          GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

          for (index_t s = 0; s!=ind.size(); ++s)
          {
              map.preImage(ind[s]+ map.freeSize(), pi);
              for (size_t k = 0; k!= pi.size(); ++k)
              {
                  if (pi[k].first == p(0))
                    result.insert(sP1+pi[k].second, c)
                              = subk(s, 0);
                  else
                    result.insert(pi[k].second, c)
                              = subk(s, 0);
              }
          }

          c++;// next function
      }

      for (index_t i = 2; i!= dgr+1+nKnots*(dgr-1); ++i)
      {
          gsMatrix<index_t> r4,r5;
          r4 = B1.boundaryOffset(iFace.first() .side(), i);
          r5 = B2.boundaryOffset(iFace.second() .side(), i);

          ind.resize(1);
          gsAsConstMatrix<index_t> rr4 = r4.reshape(sz,1);
          gsAsConstMatrix<index_t> rr5 = r5.reshape(sz,1);
          for(index_t j=0; j!=sz; ++j)
          {
              ind(0) = map.bindex(rr4(j,0), iFace.first() .patch);

              map.preImage(ind[0]+ map.freeSize(), pi);
              if (pi[0].first == p(0))
                  result.insert(sP1+pi[0].second, c)
                            = 1;
              else
                  result.insert(pi[0].second, c)
                            = 1;
              c++;
          }
          for(index_t j=0; j!=sz; ++j)
          {
              ind(0) = map.bindex(rr5(j,0), iFace.second() .patch);

              map.preImage(ind[0]+ map.freeSize(), pi);
              if (pi[0].first == p(0))
                  result.insert(sP1+pi[0].second, c)
                            = 1;
              else
                  result.insert(pi[0].second, c)
                            = 1;
              c++;
          }
       }
}

void basisPlnDegree4(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result)
{
    MatType M;
    gsDofMapper map;
    G1Matrix(mp, pdeg, M, map, false);// matrix and mapper

    //const index_t d = mp.parDim();//=2
    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    const index_t sP1 = B1.size();
    const index_t sP2 = B2.size();
    const index_t sz = B1.component(0).size();

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;
    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    const boundaryInterface & iFacePlot = *mp.iBegin();
    gsVector<index_t> p(2);
    p(0) = iFacePlot.first() .patch;
    p(1) = iFacePlot.second().patch;

    GISMO_ENSURE(4==B1.component(0).degree(0), "Degree is not 4");
    const index_t nKnots = (B1.component(1).size() - B1.component(0).degree(0)-1) / (B1.component(0).degree(0)-1);

    gsMatrix<> subm, subk;
    index_t c = 0;
    const index_t dgr=B1.component(0).degree(0);

    result = gsSparseMatrix<>(sP1+sP2, 9+4*nKnots + 2*(dgr+1+nKnots*(dgr-1)-2)*(dgr+1+nKnots*(dgr-1)));
    std::vector<std::pair<index_t,index_t> > pi;

    //Type 1.1
    order_pln1_1(iFace, r1, r2, r3, map, B1, ind);
    // form submatrix
    subm = columns(M, ind);
    subm.rrefInPlace();
    (void)rrefKerBasis(subm, subk);

     GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

     for (index_t i = 0; i!=ind.size(); ++i)
     {
         map.preImage(ind[i]+ map.freeSize(), pi);
         for (size_t k = 0; k!= pi.size(); ++k)
         {
            if (pi[k].first == p(0))
               result.insert(sP1+pi[k].second, c)
                         = subk(i, 0);
            else
               result.insert(pi[k].second, c)
                         = subk(i, 0);
         }
     }
     c++;// next function

     //Type 1.2
     order_pln1_2(iFace, r1, r2, r3, map, B1, ind);
     // form submatrix
     subm = columns(M, ind);
     subm.rrefInPlace();
     (void)rrefKerBasis(subm, subk);

     GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

     for (index_t i = 0; i!=ind.size(); ++i)
     {
         map.preImage(ind[i]+ map.freeSize(), pi);
         for (size_t k = 0; k!= pi.size(); ++k)
         {
            if (pi[k].first == p(0))
             result.insert(sP1+pi[k].second, c)
                       = subk(i, 0);
            else
             result.insert(pi[k].second, c)
                       = subk(i, 0);
         }
     }
     c++;// next function

     //Type 1.3
     for (index_t i = 0; i != 2; ++i)
     {
         order_pln_if(iFace, r1, r3, map, B1, ind, 1+i*(nKnots*3+1), 2+i*(nKnots*3+1));

         // form submatrix
         subm = columns(M, ind);
         subm.rrefInPlace();
         (void)rrefKerBasis(subm, subk);

         GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

         for (index_t s = 0; s!=ind.size(); ++s)
         {
             map.preImage(ind[s]+ map.freeSize(), pi);
             for (size_t k = 0; k!= pi.size(); ++k)
             {
                 if (pi[k].first == p(0))
                   result.insert(sP1+pi[k].second, c)
                             = subk(s, 0);
                 else
                   result.insert(pi[k].second, c)
                             = subk(s, 0);
             }
         }

         c++;// next function
     }

   //Type 1.4
     order_pln1_3(iFace, r1, r2, r3, map, B1, ind);
   // form submatrix
     subm = columns(M, ind);
     subm.rrefInPlace();
     (void)rrefKerBasis(subm, subk);

     GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

    for (index_t i = 0; i!=ind.size(); ++i)
     {
         map.preImage(ind[i]+ map.freeSize(), pi);
         for (size_t k = 0; k!= pi.size(); ++k)
         {
            if (pi[k].first == p(0))
             result.insert(sP1+pi[k].second, c)
                       = subk(i, 0);
            else
             result.insert(pi[k].second, c)
                       = subk(i, 0);
         }
     }
     c++;// next function

     //Type 1.5
       order_pln1_4(iFace, r1, r2, r3, map, B1, ind);
     // form submatrix
       subm = columns(M, ind);
       subm.rrefInPlace();
       (void)rrefKerBasis(subm, subk);

       GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

      for (index_t i = 0; i!=ind.size(); ++i)
       {
           map.preImage(ind[i]+ map.freeSize(), pi);
           for (size_t k = 0; k!= pi.size(); ++k)
           {
              if (pi[k].first == p(0))
               result.insert(sP1+pi[k].second, c)
                         = subk(i, 0);
              else
               result.insert(pi[k].second, c)
                         = subk(i, 0);
           }
       }
       c++;// next function

    //Type 2
    for (index_t i = 0; i != 2; ++i)
    {
        order_pln_oth(iFace, r1, r2, r3, map, B1, ind, 0+i*(nKnots*3+2), 2+i*(nKnots*3+2));

        // form submatrix
        subm = columns(M, ind);
        subm.rrefInPlace();
        (void)rrefKerBasis(subm, subk);

        GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

        for (index_t s = 0; s!=ind.size(); ++s)
        {
            map.preImage(ind[s]+ map.freeSize(), pi);
            for (size_t k = 0; k!= pi.size(); ++k)
            {
                if (pi[k].first == p(0))
                  result.insert(sP1+pi[k].second, c)
                            = subk(s, 0);
                else
                  result.insert(pi[k].second, c)
                            = subk(s, 0);
            }
        }

        c++;// next function
    }

      //Type 3
    for(index_t j = 0; j != nKnots;++j)
    {
      for (index_t i = 0; i != 2; ++i)
      {
          order_pln_if(iFace, r1,  r3, map, B1, ind, 2+i+j*3, 4+i+j*3);

          // form submatrix
          subm = columns(M, ind);
          subm.rrefInPlace();
          (void)rrefKerBasis(subm, subk);

          GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

          for (index_t s = 0; s!=ind.size(); ++s)
          {
              map.preImage(ind[s]+ map.freeSize(), pi);
              for (size_t k = 0; k!= pi.size(); ++k)
              {
                  if (pi[k].first == p(0))
                    result.insert(sP1+pi[k].second, c)
                              = subk(s, 0);
                  else
                    result.insert(pi[k].second, c)
                              = subk(s, 0);
              }
          }

          c++;// next function
      }
    }

        //Type 4
      for (index_t i = 0; i !=2; ++i)
        {
            order_pln_oth(iFace, r1, r2, r3, map, B1, ind, 1+i, 4+i);

            // form submatrix
            subm = columns(M, ind);
            subm.rrefInPlace();
            (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

            for (index_t s = 0; s!=ind.size(); ++s)
            {
                map.preImage(ind[s]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(s, 0);
                }
            }

            c++;// next function
      }
      for (index_t i = 0; i != nKnots-2; ++i)
        {
            order_pln_oth(iFace, r1, r2, r3, map, B1, ind, 5+3*i, 8+3*i);

            // form submatrix
            subm = columns(M, ind);
            subm.rrefInPlace();
            (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

            for (index_t s = 0; s!=ind.size(); ++s)
            {
                map.preImage(ind[s]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(s, 0);
                }
            }

            c++;// next function
      }
      for (index_t i = 0; i != 2; ++i)
        {
            order_pln_oth(iFace, r1, r2, r3, map, B1, ind, sz-6+i, sz-3+i);

            // form submatrix
            subm = columns(M, ind);
            subm.rrefInPlace();
            (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

            for (index_t s = 0; s!=ind.size(); ++s)
            {
                map.preImage(ind[s]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(s, 0);
                }
            }

            c++;// next function
      }

      //Type 5
      for (index_t i = 0; i != nKnots-1; ++i)
      {
          order_pln_oth(iFace, r1, r2, r3, map, B1, ind, 3+i*3, 7+i*3);

          // form submatrix
          subm = columns(M, ind);
          subm.rrefInPlace();
          (void)rrefKerBasis(subm, subk);

          GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

          for (index_t s = 0; s!=ind.size(); ++s)
          {
              map.preImage(ind[s]+ map.freeSize(), pi);
              for (size_t k = 0; k!= pi.size(); ++k)
              {
                  if (pi[k].first == p(0))
                    result.insert(sP1+pi[k].second, c)
                              = subk(s, 0);
                  else
                    result.insert(pi[k].second, c)
                              = subk(s, 0);
              }
          }

          c++;// next function
      }
      for (index_t i = 2; i!= dgr+1+nKnots*(dgr-1); ++i)
      {
          gsMatrix<index_t> r4,r5;
          r4 = B1.boundaryOffset(iFace.first() .side(), i);
          r5 = B2.boundaryOffset(iFace.second() .side(), i);

          ind.resize(1);
          gsAsConstMatrix<index_t> rr4 = r4.reshape(sz,1);
          gsAsConstMatrix<index_t> rr5 = r5.reshape(sz,1);
          for(index_t j=0; j!=sz; ++j)
          {
              ind(0) = map.bindex(rr4(j,0), iFace.first() .patch);

              map.preImage(ind[0]+ map.freeSize(), pi);
              if (pi[0].first == p(0))
                  result.insert(sP1+pi[0].second, c)
                            = 1;
              else
                  result.insert(pi[0].second, c)
                            = 1;
              c++;
          }
          for(index_t j=0; j!=sz; ++j)
          {
              ind(0) = map.bindex(rr5(j,0), iFace.second() .patch);

              map.preImage(ind[0]+ map.freeSize(), pi);
              if (pi[0].first == p(0))
                  result.insert(sP1+pi[0].second, c)
                            = 1;
              else
                  result.insert(pi[0].second, c)
                            = 1;
              c++;
          }
       }
}

void basisLinDegree3(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result,
                  const std::vector<gsGeometry<>::Ptr>& pert)
{
    GISMO_ASSERT(mb[0].size()>10, "Problem, need more knots per direction");

    MatType M;
    gsDofMapper map;
    gsGluingData gd(mp, /*pdeg,*/ pert);
    G1Matrix(mp, pdeg, gd, M, map);
    //G1Matrix(mp, pdeg, M, map, false);// matrix and mapper

    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    const index_t sz1 = B1.component(0).size();
    const index_t sz2 = B1.component(1).size();
    const index_t sP1 = B1.size();
    const index_t sP2 = B2.size();

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;
    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    const boundaryInterface & iFacePlot = *mp.iBegin();
    gsVector<index_t> p(2);
    p(0) = iFacePlot.first() .patch;
    p(1) = iFacePlot.second().patch;

    GISMO_ENSURE(3==B1.component(0).degree(0), "Degree is not 3");
    const index_t nKnots0 = (B1.component(0).size() - B1.component(0).degree(0)-1) / (B1.component(0).degree(0)-1);

    gsMatrix<> subm, subk;
    index_t c = 0;

    if(nKnots0<2)
        result = gsSparseMatrix<>(sP1+sP2, (2*nKnots0)*(2*nKnots0));
    if(nKnots0>1)
        result = gsSparseMatrix<>(sP1+sP2, (2*nKnots0)*(2*nKnots0)+(nKnots0-2)*(nKnots0-2));

    std::vector<std::pair<index_t,index_t> > pi;
    //Type 1

    M.rrefInPlace();
    (void)rrefKerBasis(subm, subk);

    for (index_t i = 0; i != 2*nKnots0; ++i)
    {
        for (index_t j = 0; j != 2*nKnots0; ++j)
        {
            order_ln1(iFace, r1, r3, map, ind, sz1, sz2, 2+i, 2+j);

            // form submatrix
            subm = columns(M, ind);
            subm.rrefInPlace();
            (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols() && subk.squaredNorm()!=0 ,
                         "Dimension should be 1");

            for (index_t s = 0; s!=ind.size(); ++s)
            {
                map.preImage(ind[s]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(s, 0);
                }
            }

            c++;// next function
        }
    }

    if(nKnots0>2)
        //Type 2
    {
        for (index_t i = 0; i != nKnots0-2; ++i)
        {
            for (index_t j = 0; j != nKnots0-2; ++j)
            {
                order_ln2(iFace, r1, r2, r3, map, ind, sz1, 2+i*2, 8+i*2, sz2, 2+j*2, 8+j*2, 3+i*2, 7+i*2, 3+j*2, 7+j*2);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t s = 0; s!=ind.size(); ++s)
                {
                    map.preImage(ind[s]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(s, 0);
                    }
                }

                c++;// next function
            }
        }
    }
}

void basisLinDegree4(const gsMultiPatch<> & mp,
                  const index_t pdeg,
                  gsSparseMatrix<> & result,
                  const std::vector<gsGeometry<>::Ptr>& pert)
{
    GISMO_ASSERT(mb[0].size()>10, "Problem, need more knots per direction");

    MatType M;
    gsDofMapper map;
    gsGluingData gd(mp, /*pdeg,*/ pert);
    G1Matrix(mp, pdeg, gd, M, map);

    const boundaryInterface & iFace = *mp.iBegin();// assume one single interface
    const gsBasis<> & B1 = mb[0];
    const gsBasis<> & B2 = mb[1];
    const index_t sz1 = B1.component(0).size();
    const index_t sz2 = B1.component(1).size();
    const index_t sP1 = B1.size();
    const index_t sP2 = B2.size();

    gsMatrix<index_t> r1, r2, r3;
    gsVector<index_t> ind;
    r1 = B1.boundaryOffset(iFace.first() .side(), 1);
    r2 = B1.boundaryOffset(iFace.first() .side(), 0);
    r3 = B2.boundaryOffset(iFace.second().side(), 1);

    const boundaryInterface & iFacePlot = *mp.iBegin();
    gsVector<index_t> p(2);
    p(0) = iFacePlot.first() .patch;
    p(1) = iFacePlot.second().patch;

    GISMO_ENSURE(4==B1.component(0).degree(0), "Degree is not 4");
    const index_t nKnots0 = (B1.component(0).size() - B1.component(0).degree(0)-1) / (B1.component(0).degree(0)-1);

    gsMatrix<> subm, subk;
    index_t c = 0;

    if(nKnots0<1)
        result = gsSparseMatrix<>(sP1+sP2, 1);
    if(nKnots0>0)
        result = gsSparseMatrix<>(sP1+sP2, 13*nKnots0*nKnots0+2*nKnots0+2);

    std::vector<std::pair<index_t,index_t> > pi;
    //Type 1

    M.rrefInPlace();
    (void)rrefKerBasis(subm, subk);

    for (index_t i = 0; i != 3*nKnots0+1; ++i)
    {
        for (index_t j = 0; j != 3*nKnots0+1; ++j)
        {
            order_ln1(iFace, r1, r3, map, ind, sz1, sz2, 2+i, 2+j);

            // form submatrix
            subm = columns(M, ind);
            subm.rrefInPlace();
            (void)rrefKerBasis(subm, subk);

            GISMO_ASSERT(1 == subk.cols() && subk.squaredNorm()!=0 ,
                         "Dimension should be 1");

            for (index_t s = 0; s!=ind.size(); ++s)
            {
                map.preImage(ind[s]+ map.freeSize(), pi);
                for (size_t k = 0; k!= pi.size(); ++k)
                {
                    if (pi[k].first == p(0))
                      result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                    else
                      result.insert(pi[k].second, c)
                                = subk(s, 0);
                }
            }

            c++;// next function
        }
    }

    if(nKnots0>0)
        //Type 2
    {
        for (index_t i = 0; i != nKnots0; ++i)
        {
            for (index_t j = 0; j != nKnots0; ++j)
            {
                order_ln3(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 6+i*3, sz2, 2+j*3, 6+j*3, 3+i*3, 5+i*3, 3+j*3, 5+j*3);

            // form submatrix
                subm = columns(M, ind);
                subm.rrefInPlace();
                (void)rrefKerBasis(subm, subk);

                GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                for (index_t s = 0; s!=ind.size(); ++s)
                {
                    map.preImage(ind[s]+ map.freeSize(), pi);
                    for (size_t k = 0; k!= pi.size(); ++k)
                    {
                        if (pi[k].first == p(0))
                        result.insert(sP1+pi[k].second, c)
                                = subk(s, 0);
                        else
                        result.insert(pi[k].second, c)
                                = subk(s, 0);
                    }
                }

                c++;// next function
            }
        }
    }
    if(nKnots0>1)
    {
        //Type 3
        {
            for (index_t i = 0; i != nKnots0; ++i)
            {
                for (index_t j = 0; j != nKnots0-1; ++j)
                {
                    order_ln4(iFace, r1, r2, r3, map, ind, sz1, 2+i*3, 6+i*3, sz2, 3+j*3, 8+j*3, 3+i*3, 5+i*3, 4+j*3, 7+j*3);

                // form submatrix
                    subm = columns(M, ind);
                    subm.rrefInPlace();
                    (void)rrefKerBasis(subm, subk);

                    GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                    for (index_t s = 0; s!=ind.size(); ++s)
                    {
                        map.preImage(ind[s]+ map.freeSize(), pi);
                        for (size_t k = 0; k!= pi.size(); ++k)
                        {
                            if (pi[k].first == p(0))
                            result.insert(sP1+pi[k].second, c)
                                    = subk(s, 0);
                            else
                            result.insert(pi[k].second, c)
                                    = subk(s, 0);
                        }
                    }

                    c++;// next function
                }
            }
        }
        //Type 4
        {
            for (index_t i = 0; i != nKnots0-1; ++i)
            {
                for (index_t j = 0; j != nKnots0; ++j)
                {
                    order_ln3(iFace, r1, r2, r3, map, ind, sz1, 3+i*3, 8+i*3, sz2, 2+j*3, 6+j*3, 4+i*3, 7+i*3, 3+j*3, 5+j*3);

                // form submatrix
                    subm = columns(M, ind);
                    subm.rrefInPlace();
                    (void)rrefKerBasis(subm, subk);

                    GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                    for (index_t s = 0; s!=ind.size(); ++s)
                    {
                        map.preImage(ind[s]+ map.freeSize(), pi);
                        for (size_t k = 0; k!= pi.size(); ++k)
                        {
                            if (pi[k].first == p(0))
                            result.insert(sP1+pi[k].second, c)
                                    = subk(s, 0);
                            else
                            result.insert(pi[k].second, c)
                                    = subk(s, 0);
                        }
                    }

                    c++;// next function
                }
            }
        }
        //Type 5
        {
            for (index_t i = 0; i != nKnots0-1; ++i)
            {
                for (index_t j = 0; j != nKnots0-1; ++j)
                {
                    order_ln4(iFace, r1, r2, r3, map, ind, sz1, 3+i*3, 8+i*3, sz2, 3+j*3, 8+j*3, 4+i*3, 7+i*3, 4+j*3, 7+j*3);

                // form submatrix
                    subm = columns(M, ind);
                    subm.rrefInPlace();
                    (void)rrefKerBasis(subm, subk);

                    GISMO_ASSERT(1 == subk.cols(), "Dimension should be 1");

                    for (index_t s = 0; s!=ind.size(); ++s)
                    {
                        map.preImage(ind[s]+ map.freeSize(), pi);
                        for (size_t k = 0; k!= pi.size(); ++k)
                        {
                            if (pi[k].first == p(0))
                            result.insert(sP1+pi[k].second, c)
                                    = subk(s, 0);
                            else
                            result.insert(pi[k].second, c)
                                    = subk(s, 0);
                        }
                    }

                    c++;// next function
                }
            }
        }
    }


}
