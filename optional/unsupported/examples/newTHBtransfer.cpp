

#include <gismo.h>

#include <gsRemappedBasis/gsTHB.h>

using namespace gismo;

#define dim 2

enum Action
{
    none,
    refine,
    addBox,
    prtBox,
    exitAct
};


void readBox (std::string &data, gsMatrix<real_t,dim,2> &box,  gsBoxList::basisIdT &lvl)
{
    std::istringstream input(data);

    if (!(input>>lvl) )
        throw std::logic_error("The first entry in '"+data+"' is not a valid level");
    for (int i=0;i<dim;++i)
    {
        if ( ! (input>>box(i,0)) )
            throw std::logic_error("Not a valid box begin in direction "+util::to_string(i)+"'"+input.str()+"'");
        if ( ! (input>>box(i,1)) )
            throw std::logic_error("Not a valid box begin in direction "+util::to_string(i)+"'"+input.str()+"'");
    }
}


Action readInput(std::string &cmd, gsMatrix<real_t,dim,2> &box,  gsBoxList::basisIdT &lvl)
{
    std::cout << "read the following line '"<<cmd<<"'"<<std::endl;
    static int line=0;
    ++line;
    Action result=none;

    size_t delim=cmd.find(' ');
    std::string action;
    std::string args;

    if (delim!=std::string::npos)
    {
        action = cmd.substr(0,delim);
        args   = cmd.substr(delim);
    }
    else
        action = cmd;

    if (action=="box")
    {
        readBox(args, box,lvl);
        result=addBox;
    }
    else if (action=="exit")
        result= exitAct;
    else if (action=="refine")
        result= refine;
    else if (action=="print")
        result= prtBox;
    return result;
}


void sanitizeRefinement( gsBoxList &boxes, gsTHB<dim> &sp1)
{
    for (size_t b=0; b<boxes.size(); ++b)
    {
        gsAsMatrix<real_t> box=boxes.box(b);
        const gsBoxList::basisIdT baId=boxes.basisId(b);
        const gsTHB<dim>::tensorBasisT &ba=*sp1.getTensorBases()[baId];
        index_t pos[dim];
        index_t end[dim];
        for (int dir =0; dir<dim; ++dir)
        {
            // align up to tolerance
            pos[dir]=std::lower_bound(ba.knots(dir).begin(),ba.knots(dir).end(), box(dir,0))-ba.knots(dir).begin();
            if (math::abs(ba.knots(dir)[pos[dir]]-box(dir,0))>1e-10 )
            {
                --pos[dir];
            }
            box(dir,0)=ba.knots(dir)[pos[dir]];
            end[dir]=std::lower_bound(ba.knots(dir).begin(),ba.knots(dir).end(), box(dir,1))-ba.knots(dir).begin();
            if (math::abs(ba.knots(dir)[end[dir]]-box(dir,1))>1e-10 )
            {
                --end[dir];
            }
            box(dir,1)=ba.knots(dir)[end[dir]];
        }
        if (baId<=1)
            continue;
        // insert coarser box
        const gsTHB<dim>::tensorBasisT &baC=*sp1.getTensorBases()[baId-1];
        gsMatrix<real_t> newBox(dim,2);
        for (int dir =0; dir<dim; ++dir)
        {
            // align up to tolerance
            pos[dir]=std::lower_bound(baC.knots(dir).begin(),baC.knots(dir).end(), box(dir,0))-baC.knots(dir).begin();
            if (math::abs(baC.knots(dir)[pos[dir]]-box(dir,0))>1e-10 )
            {
                --pos[dir];
            }
            newBox(dir,0)=baC.knots(dir)[pos[dir]];
            end[dir]=std::lower_bound(baC.knots(dir).begin(),baC.knots(dir).end(), box(dir,1))-baC.knots(dir).begin();
            if (math::abs(baC.knots(dir)[end[dir]-1]-box(dir,1))<1e-10 )
            {
                --end[dir];
            }
            newBox(dir,1)=baC.knots(dir)[end[dir]];
        }
        boxes.append(newBox,baId-1);
    }
}


void checkTransfer(gsTHB<dim> &sp1, gsTHBSplineBasis<dim> &sp2, gsBoxList &boxes)
{
    sanitizeRefinement(boxes, sp1);
    double newT,oldT;

    gsSparseMatrix<real_t> tr1S;
    gsSparseMatrix<real_t> tr1,tr2;

    std::vector<index_t> elements;
    boxes.toRefineElementFormat(sp1.getBases(),elements);

    gsStopwatch clock;
    sp1.refineElementsWithTransfer(elements,tr1S);
    newT=clock.stop();

    clock.restart();
    sp2.refineElements_withTransfer(elements,tr2);
    oldT=clock.stop();

    tr1 = tr1S;

    real_t diffNorm=(tr1-tr2).norm() ;
    if ( diffNorm > 1e-10 )
    {
        std::cout.unsetf ( std::ios::floatfield );
        std::cout.precision(2);
        std::cout.setf(std::ios::scientific);
        std::cout << "transfer new:\n"<<tr1<<"\ntransfer old:\n"<<tr2<<"\ndifference\n"<<(tr1-tr2)<<std::endl;
    }
    std::cout<<"diff norm is: "<<diffNorm<<"  times are: "<<newT<<"  "<<oldT<<std::endl;
}


std::vector<gsBoxList::basisPtr> getBases(gsTHBSplineBasis<dim>  &sp2)
{
    typedef gsTHBSplineBasis<dim>::tensorBasis tensorBasis;

    std::vector<gsBoxList::basisPtr>      result;
    const std::vector<tensorBasis*> &tBas=sp2.getBases();

    for (size_t i=0; i<tBas.size();++i)
        result.push_back(gsBoxList::basisPtr(tBas[i]->clone().release()));

    tensorBasis::uPtr tmp(tBas.back());
    for (size_t i=result.size();i<12;++i)
    {
        tmp->uniformRefine();
        result.push_back(gsBoxList::basisPtr(tmp->clone().release()));
    }
    return result;
}

int main(int argc, char *argv[])
{
    bool run = false;
    gsCmdLine cmd("test computation of transfer against old implementation");
    cmd.addSwitch ("R", "run", "",   run);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (!run)
    {
        gsInfo<<"I did nothing! :P"<<std::endl;
        return 0;
    }

    gsMatrix<real_t,dim,2> box;
    gsBoxList::basisIdT    lvl;

    gsBoxList              curRef(dim);

    gsKnotVector<real_t>   knot(0,1,0,3);
    gsTHBSplineBasis<dim>::tensorBasis lev0(std::vector<gsKnotVector<> >(2,knot));

    gsTHBSplineBasis<dim>  sp2(lev0);
    gsTHB<dim>             sp1(getBases(sp2), curRef);

    std::cout<<"The domain is \n"<<sp2.support()<<std::endl;

    std::string line;
    while (!std::cin.eof())
    {
        getline(std::cin,line,'\n');
        switch (readInput(line, box, lvl))
        {
        case none:
            if (line!="")
                std::cout<<"Invalid input '"<<line<<"\n"<<std::flush;
            continue;
            break;
        case refine:
            std::cout<<"refi"<<std::endl;
            checkTransfer(sp1,sp2,curRef);
            curRef=gsBoxList(2);
            break;
        case addBox:
            std::cout<<"add"<<std::endl;
            curRef.append(box,lvl);
            break;
        case prtBox:
            for (size_t b=0; b<curRef.size();++b)
                std::cout<<"set\n"<<curRef.box(b)<<"\n to lvl "<<curRef.basisId(b)<<"\n";
            std::cout<<std::endl;
            break;
        case exitAct:
            std::cout<<"exit"<<std::endl;
            goto out;
        default:
            return 1;
        };
    }
out:
    return 0;
}

