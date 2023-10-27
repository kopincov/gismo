
/*
* gsBSplineRootsTest.cpp created on 29.07.2014
*
* Author: Andrea Bressan
*
* This file is part of the G+SMO library
*/


#include <gismo.h>
#include <gismo_dev.h>
#include <gsNurbs/gsBSplineSolver.h>
#include <gsIO/gsWriteParaview.h>
#include <gsIO/gsFileData.h>


using namespace gismo;

int puzzlePiece (int argc, char** argv);


int main (int argc, char** argv)
{
    gsMatrix<> coefs(12,1);
    gsBSplineBasis<> my_basis(0,10,9,2);
    std::vector<Root<> > roots;
    gsVector<>  normal(1);
    normal(0)=1;
    unsigned count=0;

    coefs<< 0,1,2,3,0,0,0,0,1,2,3,4;
    gsBSpline<> my_curve(my_basis,coefs);
    count = findHyperPlaneIntersections<real_t>(my_curve, normal,(real_t) 0, 1000*math::limits::epsilon(),roots);
    if (
        count != 2
        || roots[0].type != odd
        || roots[0].parameter != 0

        || roots[1].type != even_interval
        || roots[1].begParameter != 4
        || roots[1].endParameter != 6
        )
        return 1;
    roots.resize(0);

    coefs<< 0,1,2,3,0,0,0,0,-1,-2,-3,-4;
    
    my_curve = gsBSpline<>(my_basis,coefs);
    count = findHyperPlaneIntersections<real_t>(my_curve, normal,(real_t) 0, 1000*math::limits::epsilon(),roots);
    if (
        count != 2
        || roots[0].type != odd
        || roots[0].parameter != 0

        || roots[1].type != odd_interval
        || roots[1].begParameter != 4
        || roots[1].endParameter != 6
        )
        return 1;
    roots.resize(0);

    coefs<< 0,0,0,3,0,0,0,0,-1,-2,-3,-4;
    my_curve = gsBSpline<>(my_basis,coefs);
    count = findHyperPlaneIntersections<real_t>(my_curve, normal,(real_t) 0, 1000*math::limits::epsilon(),roots);
    if (
        count != 2
        || roots[0].type != odd_interval
        || roots[0].begParameter != 0
        || roots[0].endParameter != 1

        || roots[1].type != odd_interval
        || roots[1].begParameter != 4
        || roots[1].endParameter != 6
        )
        return 1;
    roots.resize(0);




    coefs.resize(13,1);
    my_basis.insertKnot(5);

    coefs<< 1,1,1,1,1,1,0,1,1,1,1,1,2;
    my_curve = gsBSpline<>(my_basis,coefs);
    count = findHyperPlaneIntersections<real_t>(my_curve, normal,(real_t) 0, 1000*math::limits::epsilon(),roots);
    if (
        count != 1
        || roots[0].type != even
        || roots[0].parameter != 5
        )
        return 1;
    roots.resize(0);


    coefs<< 1,1,1,1,1,1,0,-1,-1,-1,-1,-1,-2;
    my_curve = gsBSpline<>(my_basis,coefs);
    count = findHyperPlaneIntersections<real_t>(my_curve, normal,(real_t) 0, 1000*math::limits::epsilon(),roots);
    if (
        count != 1
        || roots[0].type != odd
        || roots[0].parameter != 5
        )
        return 1;
    roots.resize(0);

    coefs<< 1,1,1,1,1,0,0,1,1,1,1,1,2;
    my_curve = gsBSpline<>(my_basis,coefs);
    count = findHyperPlaneIntersections<real_t>(my_curve, normal,(real_t) 0, 1000*math::limits::epsilon(),roots);
    if (
        count != 1
        || roots[0].type != even
        || roots[0].parameter != 5
        )
        return 1;
    roots.resize(0);

    coefs<< 1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-2;
    my_curve = gsBSpline<>(my_basis,coefs);
    count = findHyperPlaneIntersections<real_t>(my_curve, normal,(real_t) 0, 1000*math::limits::epsilon(),roots);
    if (
        count != 1
        || roots[0].type != odd
        || roots[0].parameter != 4
        )
        return 1;
    roots.resize(0);

    coefs<< 1,1,1,1,1,1,1,1,1,1,1,1,0;
    my_curve = gsBSpline<>(my_basis,coefs);
    count = findHyperPlaneIntersections<real_t>(my_curve, normal,(real_t) 0, 1000*math::limits::epsilon(),roots);
    if (
        count != 1
        || roots[0].type != odd
        || roots[0].parameter != 10
        )
        return 1;
    roots.resize(0);

    coefs<< 1,1,1,1,1,1,1,1,1,1,0,0,0;
    my_curve = gsBSpline<>(my_basis,coefs);
    count = findHyperPlaneIntersections<real_t>(my_curve, normal,(real_t) 0, 1000*math::limits::epsilon(),roots);
    if (
        count != 1
        || roots[0].type != odd_interval
        || roots[0].begParameter != 9
        || roots[0].endParameter != 10
        )
        return 1;
    roots.resize(0);

    return puzzlePiece (argc, argv);
}




int puzzlePiece (int argc, char** argv)
{
    typedef gsBSpline<real_t> gsCurveRR;


    std::string fn("planar/planarDomainPuzzle1.xml");

    gsPlanarDomain<>::uPtr pDomain;
    gsCurveRR::uPtr        curve;
    gsMatrix<>             boundingBox;
    index_t                    nSample = 1000;
    
    gsCmdLine cmd("Hi, give me a file (.txt, .axl) and I will try to draw it!");
    cmd.addString("I", "input", "Input planar domain", fn);
    cmd.addInt("s","samples", "Number of samples ordinates", nSample);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<>  file(fn);

    if (file.has< gsPlanarDomain<> >() && nSample>0)
    {
        pDomain = file.getFirst< gsPlanarDomain<> >();
    }
    else
    {
        gsWarn << "Input arguments are \"strange\". Aborting...\n";
        return 1;
    }


    gsCurveLoop<>& cLoop = pDomain->outer();
    curve  = memory::convert_ptr<gsCurveRR>(cLoop.singleCurve());
    if (!curve)
    {
        gsWarn << "The input does not have a BSpline curve boundary. Aborting...\n";
        return 1;
    }
    boundingBox = pDomain->boundingBox();


    real_t y;
    real_t ss=nSample-1;
    std::vector<Root<> >           xs;
    gsMatrix<>                     points;

    gsVector<>                     normal(2);
    normal(0)=.5;normal(1)=.4;

    std::vector<unsigned>   count(nSample);

    for (int i=0;i<nSample;++i)
    {
        real_t ii=i;
        y= ii/ss * boundingBox(1,0) + (ss-ii)/ss*boundingBox(1,1);
        count[i]=findHyperPlaneIntersections<real_t>(*curve,normal,y, 1000*math::limits::epsilon(), xs);
        gsInfo<<"\n"<<"line: "<<i<<"\n";
    }

    points.resize(2,xs.size());

    unsigned j=0, k=0;
    while (count[j]==0)
    {
        ++j;
    }

    for (unsigned i=0; i< xs.size(); ++i,++k)
    {
        switch (xs[i].type)
        {
            case odd:
            case even:
                points.col(i)=xs[i].point;
                break;
            default:
                points.col(i)=(xs[i].begPoint+xs[i].endPoint)/2;
        }
        if (k==count[j])
        {
            k=0;
            ++j;
        }
        while (count[j]==0)
        {
            ++j;
        }
    }

    gsWriteParaview(*curve, std::string("RootTestCurve"));
    gsWriteParaviewPoints(points,"RootTestCurveIntersections"); // columns in a matrix

    return 0;
}
