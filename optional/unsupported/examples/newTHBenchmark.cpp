
#include <string>
#include <exception>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsRemappedBasis/gsTHB.h>
#include <gsRemappedBasis/gsHLR.h>
#include <gsIO/gsIOUtils.h>

#include <iomanip>


using namespace gismo;


struct Timing
{
    double initTime;
    double evalTime;
    double smartEval;

    Timing()
        : initTime(NAN), evalTime(NAN), smartEval(NAN)
    {}
};
std::ostream& operator<< (std::ostream & out, const Timing t)
{
    std::ios::fmtflags old_settings = out.flags();
    std::streamsize old_precision = out.precision();
    out.precision(3);
    out.setf(std::ios_base::scientific);
    out<<std::setw(10)<< t.initTime<<"   "<<std::setw(10)<<t.evalTime<<"   "<<std::setw(10)<<t.smartEval;
    out.flags(old_settings);
    out.precision(old_precision);
    return out;
}

// utility functions

gsBoxList getBoxes ();
void evalSimple(Timing &time);
void evalSmart(Timing &time);

void compareValues(real_t &diff);

// methods of construction

typedef gsFunctionSet<>* (*constructor) (Timing &time, const gsBoxList &);

gsFunctionSet<> *getOldThb(Timing &time, const gsBoxList &boxes);
gsFunctionSet<> *getNewThb(Timing &time, const gsBoxList &boxes);


constructor methods[]={getNewThb,getOldThb};
std::string methodName[]={"new","old"};
int         methodNum=sizeof(methods) / sizeof(constructor);

// global data

index_t     grading=-1;
index_t     degree=2;
index_t     level=3;
index_t     numPoints=100000;
bool        truncOpt=true;


gsFunctionSet<> *basis;
gsFunctionSet<> *basisCMP;
std::string      geoName="curves2d/e-curve.xml";

int main (int argn, char** args)
{
    gsCmdLine cmd("test running time of a hierarchiacally refined basis");

    index_t       m=0;
    bool        cmp=false;

    cmd.addString ("c", "curve",   "name of a file containing the refining curve", geoName);
    cmd.addInt    ("m", "method",  "either 0 for new or 1 for old", m );
    cmd.addInt    ("g", "grading",    "",   grading );
    cmd.addInt    ("d", "degree",     "",   degree );
    cmd.addInt    ("l", "level",      "",   level );
    cmd.addInt    ("n", "number",     "",   numPoints );
    cmd.addSwitch ("D", "difference", "",   cmp);
    cmd.addSwitch ("T", "truncate",   "",   truncOpt);

    try { cmd.getValues(argn,args); } catch (int rv) { return rv; }

    if (grading==-1)
        grading=degree-1;

    gsBoxList boxes=getBoxes();

    if (cmp)
    {
        std::vector<Timing> times(methodNum);
        std::vector<gsFunctionSet<>*> bases(methodNum);

        for (index_t i=0;i<methodNum;++i)
        {
            bases[i]=methods[i](times[i],boxes);
        }

        for (index_t i=0;i<methodNum;++i)
        {
            basis=bases[i];
            evalSimple(times[i]);
            evalSmart(times[i]);
            gsInfo<<"\n"<<methodName[i]<<" "<<times[i];
        }

        real_t maxDiff;
        for (index_t i=1;i<methodNum;++i)
        {
            gsInfo<<"\nCOMPARISON OF "<<methodName[0]<<" WITH "<<methodName[i]<<"\n";
            basis=bases[0];
            basisCMP=bases[i];
            compareValues(maxDiff);
            gsInfo<<"maxDiff("<<methodName[0]<<", "<<methodName[i]<<") = "<<maxDiff<<"\n";
        }
        freeAll(bases);
    }
    else
    {
        Timing times;
        basis=methods[m](times,boxes);
        evalSimple(times);
        evalSmart (times);

        gsInfo<<times<<"\n";

        delete basis;
    }

    return 0;
}


bool inside(const gsMatrix<> &p, const gsMatrix<> &box)
{
    return (p.array()>=box.col(0).array()).all() && (p.array()<=box.col(1).array()).all();
}

index_t roundUp(index_t numToRound, index_t multiple)
{
    numToRound-=degree;
    index_t remainder = numToRound % multiple;
    return (remainder ? numToRound + multiple - remainder  : numToRound )+degree;
}
index_t roundDown(index_t numToRound, index_t multiple)
{
    numToRound-=degree;
    index_t remainder = numToRound % multiple;
    return numToRound-remainder+degree;
}

gsKnotVector<real_t> getKnotVector(index_t levelParam)
{
    gsKnotVector<> knots(0,1,0,degree+1);
    gsBSplineBasis<> tbasis(knots);

    gsTHBSplineBasis<1> oldB(tbasis);
    oldB.tensorLevel(levelParam); // force existence of the level
    // boxes
    return oldB.getBases()[levelParam]->knots(0);
}

gsBoxList getBoxes ()
{
    gsKnotVector<real_t> knots=getKnotVector(level);
    // trace curve and give the boxes at maximum level that intersect the curve
    // append boxes for grading
    gsBSpline<>::uPtr curve = gsReadFile<>(geoName);
    if (!curve || curve->targetDim()!=2 || curve->domainDim()!=1 )
        throw std::logic_error("The provided geometry is not a curve");

    gsMatrix<real_t> domain(2,2);
    domain<<0,1,0,1;
    gsMatrix<real_t> supp=curve->support();
    gsMatrix<real_t> curBox(2,2);
    curBox<<0,0,0,0;

    gsBoxList        result(2);
    gsMatrix<>       curParam=supp.col(0);
    gsMatrix<>       curPoint=curve->eval(curParam);
    gsMatrix<>       tmpPoint;
    gsMatrix<>       tmpParam;

    const index_t max=knots.size()-1;

    real_t tol=math::pow((real_t)(2),-level-1);
    gsMatrix<> parIncr(1,1);
    parIncr(0,0)=math::pow((real_t)(2),-level-2)*(supp.col(1)-supp.col(0)).norm();

    index_t minElements=degree+1;

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

        // add bigger box so we have a function at this level
        ind0l=roundDown(math::max(static_cast<index_t>(0), ind0l-(minElements+3)/4),2);
        ind0u=roundUp  (math::min(max                    , ind0u+(minElements+3)/4),2);
        ind1l=roundDown(math::max(static_cast<index_t>(0), ind1l-(minElements+3)/4),2);
        ind1u=roundUp  (math::min(max                    , ind1u+(minElements+3)/4),2);

        curBox(0,0)=knots[ind0l];
        curBox(0,1)=knots[ind0u];
        curBox(1,0)=knots[ind1l];
        curBox(1,1)=knots[ind1u];

        result.append(curBox, static_cast<gsBoxList::basisIdT>(level));

        for (gsBoxList::basisIdT l=static_cast<gsBoxList::basisIdT>(level)-1; l>0 ; --l)
        {
            const index_t factor=1<<(level-l);
            ind0l=ind0l-(grading+1)/2*factor;
            ind0u=ind0u+(grading+1)/2*factor;
            ind1l=ind1l-(grading+1)/2*factor;
            ind1u=ind1u+(grading+1)/2*factor;
            // round
            ind0l=math::max(static_cast<index_t>(0), roundDown(ind0l,2*factor));
            ind0u=math::min(max,                     roundUp  (ind0u,2*factor));
            ind1l=math::max(static_cast<index_t>(0), roundDown(ind1l,2*factor));
            ind1u=math::min(max,                     roundUp  (ind1u,2*factor));

            curBox(0,0)=knots[ind0l];
            curBox(0,1)=knots[ind0u];
            curBox(1,0)=knots[ind1l];
            curBox(1,1)=knots[ind1u];

            result.append(curBox, l);
        }

        while (static_cast<size_t>(curInd0)==indexOfLastLessOrEqual(knots,curPoint(0,0))
               && static_cast<size_t>(curInd1)==indexOfLastLessOrEqual(knots,curPoint(1,0))) // find a point in a different cell
        {
            do // search for a point close enough to respect the tolerance
            {
                tmpParam=curParam+parIncr;
                if (tmpParam(0,0)>supp.col(1)(0))
                    return result;
                tmpPoint=curve->eval(tmpParam);
                parIncr/=2;
            }
            while((tmpPoint-curPoint).norm()>tol && parIncr.norm()>std::numeric_limits<real_t>::epsilon());

            parIncr(0,0)=math::pow((real_t)(2),-level-1)*(supp.col(1)-supp.col(0)).norm();
            curPoint=tmpPoint;
            curParam=tmpParam;
        }
    }

    return result;
}

void evalSimple(Timing &time)
{
    gsMatrix<> bbox(2,2);
    bbox<<0,1,0,1;
    gsFuncData<> result(NEED_VALUE|NEED_DERIV);
    const gsMatrix<> points = gsPointGrid(bbox, numPoints);
    gsStopwatch clock;
    basis->compute(points,result);
    time.evalTime=clock.stop();
}

gsDomainIterator<> * getIterator()
{
    gsHTensorBasis<2> *thbold=dynamic_cast<gsHTensorBasis<2>*>(basis);
    if (thbold)
        return thbold->makeDomainIterator().release();
    gsTHB<2> * thbnew = dynamic_cast<gsTHB<2> *>(basis);
    if (thbnew)
        return thbnew->makeDomainIterator().release();
    throw std::logic_error("proper basis needed");
}


void evalSmart(Timing &time)
{
    // for each element
    gsDomainIterator<> *it=getIterator();
    gsStopwatch clock;
    gsVector<index_t> numQuad;
    numQuad.setConstant(2,degree+1);
    gsGaussRule<> quad(numQuad);

    gsMatrix<> nodes;
    gsVector<> weights;

    gsFuncData<> result(NEED_VALUE|NEED_DERIV|SAME_ELEMENT);
    for(;it->good();it->next())
    {
        quad.mapTo(it->lowerCorner(), it->upperCorner(),nodes, weights);
        basis->compute(nodes,result); // make sure this is not eliminated from the compiler
    }
    time.smartEval=clock.stop();
    delete it;
}

gsMatrix<> removeZeroRows(const gsMatrix<>& mat)
{
    gsMatrix<> m(mat.rows(), mat.cols());

    index_t counter = 0;
    for (index_t r = 0; r != mat.rows(); r++)
    {
        bool zeroRow = true;

        for (index_t c = 0; c != mat.cols(); c++)
        {
            if (mat(r, c) != 0.0)
            {
                zeroRow = false;
                break;
            }
        }
        if (!zeroRow)
        {
            m.row(counter) = mat.row(r);
            counter++;
        }
    }
    m.conservativeResize(counter, mat.cols());
    return m;
}




void compareValues(real_t &maxDiff)
{
    // for each element
    gsDomainIterator<> *it=getIterator();

    gsVector<index_t> numQuad;
    numQuad.setConstant(2,degree+1);
    gsGaussRule<> quad(numQuad);

    gsMatrix<> nodes;
    gsVector<> weights;

    gsFuncData<> resultN(NEED_ACTIVE|NEED_VALUE|NEED_DERIV);
    gsFuncData<> resultO(NEED_ACTIVE|NEED_VALUE|NEED_DERIV);

    maxDiff=0;
    for(;it->good();it->next())
    {
        quad.mapTo(it->lowerCorner(), it->upperCorner(),nodes, weights);
        basisCMP->compute(nodes, resultN);
        basis->compute(nodes, resultO);
        gsMatrix<> evalN = removeZeroRows(resultN.values[0]);
        gsMatrix<> evalO = removeZeroRows(resultO.values[0]);

        if (evalN.rows() == evalO.rows() && evalN.cols() == evalO.cols())
            maxDiff=math::max((evalN - evalO).norm(),maxDiff);
        else
        {
            gsInfo << "eval: " << evalN.rows() << " x " << evalN.cols() << " =?= "
                      << evalO.rows() << " x " << evalO.cols() << "   "
                      << resultN.actives.cols() << " =?= " << resultO.actives.cols() << "\n\n"
                      << evalN << "\n\n"
                      << evalO << "\n\n"
                      << resultN.actives << "\n\n"
                      << resultO.actives << "\n\n"
                      << evalN.colwise().sum()<< "\n\n"
                      << evalO.colwise().sum()<< "\n\n"
                      << ((it->lowerCorner()+it->upperCorner()).array()/2).transpose()<< "\n\n"
                      << "\n";
            maxDiff=std::numeric_limits<real_t>::infinity();
            break;
        }
    }
    delete it;
}

gsFunctionSet<> *getOldThb(Timing &time, const gsBoxList &boxes)
{
    gsKnotVector<> knots(0,1,0,degree+1);
    gsTensorBSplineBasis<2> tbasis(knots,knots);

    gsHTensorBasis<2> *oldB = truncOpt? static_cast<gsHTensorBasis<2> *>(new gsTHBSplineBasis<2>(tbasis)):static_cast<gsHTensorBasis<2> *>(new gsHBSplineBasis<2>(tbasis));
    oldB->tensorLevel(level); // force existence of the level

    std::vector<index_t> boxOld;
    boxes.toRefineElementFormat(oldB->getBases(),boxOld);

    gsStopwatch clock;
    // clock only the refinement as we somehow forsed unnecessary computations before
    oldB->refineElements(boxOld);
    time.initTime=clock.stop();

    gsMesh<real_t> msh;
    makeMesh<real_t>(*oldB, msh,8);
    gsWriteParaview(msh, "THBG", false);

    return oldB;
}


gsFunctionSet<> *getNewThb(Timing &time, const gsBoxList &boxes)
{
    typedef gsTensorBSplineBasis<2,real_t> tensorBasisT;
    gsTHB<2>* result;

    gsKnotVector<> knots(0,1,0,degree+1);
    tensorBasisT::uPtr current(new tensorBasisT(knots,knots));

    std::vector<gsBoxList::basisPtr> bases;
    bases.push_back(gsBoxList::basisPtr(current->clone().release()));
    for (index_t i=1; i<=level;++i)
    {
        current->uniformRefine();
        bases.push_back(gsBoxList::basisPtr(current->clone().release()));
    }

    gsStopwatch clock;
    result = new gsTHB<2>(bases,boxes,truncOpt);
    time.initTime=clock.stop();

    result->exportSelectorToTex("THB");
    int errNo=system("pdflatex  THB.tex >/dev/null 2>&1 &"); GISMO_UNUSED(errNo);

    return result;
}

