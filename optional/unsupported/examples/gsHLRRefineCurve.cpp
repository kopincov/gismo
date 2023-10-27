
#include <gismo.h>

#include <gsRemappedBasis/gsHLR.h>
#include <gsRemappedBasis/gsTHB.h>
#include <gsRemappedBasis/otherUtils.h>

using namespace gismo;



class gsMCurve
{
public:
    gsMatrix<> support() const
    {
        gsMatrix<> result(1,2);
        result<<0,1;
        return result;
    }

    gsMatrix<> eval (const gsMatrix<> & points) const
    {
        gsMatrix<> res;
        res.resize(2,points.cols());
        for (index_t p=0; p< points.cols();++p)
            res.col(p)=eval(points(0,p));
        return res;
    }

    gsMatrix<real_t,2,1> eval(real_t t)const
    {
        if(t<0) t=0;
        if(t>1) t=1;
        gsMatrix<real_t,2,1> res;

        real_t soms=.4*(1-t*t);
        res << soms * math::sin(13.0*t)+.5, soms * math::cos(13.0*t)+.5;
        return res;
    }

};

gsKnotVector<real_t> getKnotVector(int lvl,int deg)
{
    gsKnotVector<> knots(0,1,0,deg+1);
    for (int i=0;i<lvl;++i)
        knots.uniformRefine();
    return knots;
}

bool inside(const gsMatrix<> &p, const gsMatrix<> &box)
{
    return (p.array()>=box.col(0).array()).all() && (p.array()<=box.col(1).array()).all();
}

int roundUp(int numToRound, int multiple, int deg)
{
    numToRound-=deg;
    int remainder = numToRound % multiple;
    return (remainder ? numToRound + multiple - remainder  : numToRound )+deg;
}
int roundDown(int numToRound, int multiple, int deg)
{
    numToRound-=deg;
    int remainder = numToRound % multiple;
    return numToRound-remainder+deg;
}

gsBoxList getBoxes (int lvl,int deg)
{
    gsMCurve curve;

    gsKnotVector<real_t> knots=getKnotVector(lvl,deg);
    // trace curve and give the boxes at maximum lvl that intersect the curve
    // append boxes for grading
    gsMatrix<real_t> domain(2,2);
    domain<<0,1,0,1;
    gsMatrix<real_t> supp=curve.support();
    gsMatrix<real_t> curBox(2,2);
    curBox<<0,0,0,0;

    gsBoxList        result(2);
    gsMatrix<>       curParam=supp.col(0);
    gsMatrix<>       curPoint=curve.eval(curParam);
    gsMatrix<>       tmpPoint;
    gsMatrix<>       tmpParam;

    const index_t max=knots.size()-1;

    real_t tol=math::pow((real_t)(2),-lvl-1);
    gsMatrix<> parIncr(1,1);
    parIncr(0,0)=math::pow((real_t)(2),-lvl-2)*(supp.col(1)-supp.col(0)).norm();

    index_t minElements=deg+1;

    while (!inside(curPoint,domain) )
    {
        curParam+=parIncr;
    }
    while (inside(curPoint,domain))
    {
        index_t ind0l=indexOfLastLessOrEqual(knots,curPoint(0,0));
        index_t ind0u=indexOfFirstGreaterOrEqual(knots,curPoint(0,0));
        if (ind0l>=ind0u)
            ind0u=ind0l+1; // we are on the knot line
        index_t ind1l=indexOfLastLessOrEqual(knots,curPoint(1,0));
        index_t ind1u=indexOfFirstGreaterOrEqual(knots,curPoint(1,0));
        if (ind1l>=ind1u)
            ind1u=ind1l+1; // we are on the knot line

        index_t curInd0=ind0l;
        index_t curInd1=ind1l;

        // add bigger box so we have a function at this lvl
        ind0l=roundDown(math::max(static_cast<index_t>(0),   ind0l-(minElements+3)/4),2,deg);
        ind0u=roundUp  (math::min(max, ind0u+(minElements+3)/4),2,deg);
        ind1l=roundDown(math::max(static_cast<index_t>(0),   ind1l-(minElements+3)/4),2,deg);
        ind1u=roundUp  (math::min(max, ind1u+(minElements+3)/4),2,deg);

        curBox(0,0)=knots[ind0l];
        curBox(0,1)=knots[ind0u];
        curBox(1,0)=knots[ind1l];
        curBox(1,1)=knots[ind1u];

        result.append(curBox, static_cast<gsBoxList::basisIdT>(lvl));

        while (static_cast<size_t>(curInd0)==indexOfLastLessOrEqual(knots,curPoint(0,0))
               && static_cast<size_t>(curInd1)==indexOfLastLessOrEqual(knots,curPoint(1,0))) // find a point in a different cell
        {
            do // search for a point close enough to respect the tolerance
            {
                tmpParam=curParam+parIncr;
                if (tmpParam(0,0)>supp.col(1)(0))
                    goto end;
                tmpPoint=curve.eval(tmpParam);
                parIncr/=2;
            }
            while((tmpPoint-curPoint).norm()>tol && parIncr.norm()>std::numeric_limits<real_t>::epsilon());

            parIncr(0,0)=math::pow((real_t)(2),-lvl-1)*(supp.col(1)-supp.col(0)).norm();
            curPoint=tmpPoint;
            curParam=tmpParam;
        }
    }
end:
    return result;
}

gsMatrix<> computeInverseDist( const gsMatrix<> &pts/*, int lvl*/)
{
    static gsMatrix<> values(2,0);
    if (values.cols()< 8000 ) // ?
    {
        gsVector<> low(1); low<<0;
        gsVector<> upp(1); upp<<1;
        gsVector<unsigned> np(1); np<<8000;
        const gsMatrix<> tmp = gsPointGrid(low,upp,np);
        values=gsMCurve().eval(tmp);
    }

    gsMatrix<> res(1,pts.cols());
    for (index_t p=0; p<pts.cols();++p)
    {
        res(0,p)=math::min((real_t)(100),1/(values.colwise()-pts.col(p)).colwise().norm().minCoeff());
    }
    return res;
}


int main (int argn, char** args)
{
    index_t deg=2;
    index_t lvl=2;

    gsCmdLine cmd("do something");
    cmd.addInt    ("l", "lvl",   "",  lvl );
    cmd.addInt    ("d", "deg",  "",  deg );
    try { cmd.getValues (argn,args); } catch (int rv) { return rv; }

    gsBoxList boxes = getBoxes(lvl,deg);

    typedef gsTensorBSplineBasis<2,real_t> tensorBasisT;

    gsKnotVector<> knots(0,1,0,deg+1);
    tensorBasisT current(knots,knots);
    std::vector<gsBoxList::basisPtr> bases;

    bases.push_back( current.clone() );
    //bases.push_back(gsBoxList::basisPtr(current.clone().release()));
    for (int i=1; i<=lvl;++i)
    {
        current.uniformRefine();
        bases.push_back( current.clone() );
    }

    gsStopwatch clock;
    gsHLR<2> space(bases,boxes);
    gsInfo << "space built in : "<<clock.stop()<< "seconds\n"<<std::flush; clock.restart();

    /*
    space.exportDefinitionToTex("HLRcurveLevels");
    int pdfLatexErr=0;
    pdfLatexErr=system(
        "(pdflatex HLRcurveLevels.tex >/dev/null ");
    GISMO_UNUSED(pdfLatexErr);
    gsInfo << "domain exported in : "<<clock.stop()<< "seconds\n"<<std::flush; clock.restart();
    */

    gsBasis<>::domainIter it=space.makeDomainIterator();
    gsVector<index_t> np(2);
    np.setConstant(2,2);
    gsGaussRule<> rul(np);
    gsVector<> unused;
    gsMatrix<> pts;
    gsFuncData<> data(NEED_ACTIVE|NEED_VALUE|SAME_ELEMENT);
    gsMatrix<>   distInv = computeInverseDist(pts);

    std::vector<real_t> rhsData;
    gsSparseEntries<real_t>   triplets;
    triplets.reserve(5*space.size());
    index_t totalP=0;
    for( ; it->good(); it->next() )
    {
        rul.mapTo(it->lowerCorner(), it->upperCorner(), pts,unused);
        space.compute(pts,data);
        distInv = computeInverseDist(pts);

        rhsData.insert(rhsData.end(),distInv.data(),distInv.data()+distInv.size());
        for (index_t a=0; a<data.actives.rows();++a)
            for (index_t p=0; p<pts.cols();++p)
                triplets.add(totalP+p,data.actives(a,0),data.values[0](a,p));
        totalP+=pts.cols();
    }
    gsSparseMatrix<> pr(totalP,space.size());
    pr.setFromTriplets(triplets.begin(),triplets.end());
    pr.makeCompressed();

    gsInfo << "data computed in : "<<clock.stop()<< "seconds\n"<<std::flush; clock.restart();

    gsSparseMatrix<> sys=pr.transpose()*pr;
    gsMatrix<>       rhs=pr.transpose()*gsAsVector<>(rhsData);
    gsInfo << "system assembled in : "<<clock.stop()<< "seconds\n"<<std::flush; clock.restart();

    sys.makeCompressed();
    gsMatrix<> coefs =  gsSparseSolver<>::LU(sys).solve(rhs);
    gsInfo << "solution computed in : "<<clock.stop()<< "seconds\n"<<std::flush;
    gsRemappedBasis *sol=space.makeReduced(coefs, false);

    printAllFunctions(*sol,"HLRinterp",11000);

    delete sol;
    return 0;
}
