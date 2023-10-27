/** @file

    @brief gsStokesInfSupTest export matrices needed for computing the infsup in matlab

    Author(s): A. Bressan
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsRecipeAssembler/gsRecipeAssemblerStokes.h>
#include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>

#include <gsRecipeAssemblerAdaptive/gsAdaptiveSolver.h>
#include <gsRecipeAssemblerAdaptive/gsCompositeBasisSpaceRefiners.h>
#include <gsRecipeAssemblerAdaptive/gsTensorBasisSpaceRefiners.h>

#include <gsNurbs/gsNurbsCreator.h>
#include <gsUtils/gsExportMatrix.h>

#include <algorithm>

using namespace gismo;


bool mstrcmp(const char* a,const char* b)
{
    return !strcmp(a,b);
}


struct InfSupCommonOptions
{
    std::string     m_name;
    int             m_argn;
    char**          m_args;

    index_t         m_deg;
    index_t         m_smooth;

    index_t         m_level;
    index_t         m_levSep;
    index_t         m_numElLev0;
    index_t         m_boxSize;

    index_t         m_macroSize;

    index_t         m_method;
    bool            m_basisFromDomain;

    InfSupCommonOptions(std::string name, int argn, char** args)
        : m_name(name), m_argn(argn), m_args(args)
    {
        static const char* endArg[2]={"--", NULL};
        if (argn==0)
        {
            m_args=const_cast<char**>(endArg);
            m_argn=1;
        }
        else
        {
            m_argn = std::search(args+1,args+argn,endArg,endArg+1,mstrcmp)-args;
        }
        gsCmdLine   cmd ("InfSupTest common options");
        std::string methodName("TH");
        m_deg=2;
        m_smooth=-1;
        m_level=2;
        m_boxSize=-1;
        m_levSep=-1;
        m_numElLev0=1; // number of elements at level 0
        m_basisFromDomain=false;

        cmd.addString("m", "method",          "which element: TH, SG, RT?",    methodName);
        cmd.addInt   ("s", "smoothness",      "snoothness of the spaces",      m_smooth);
        cmd.addInt   ("d", "degree",          "Degree of the spaces",          m_deg);
        cmd.addInt   ("l", "level",           "level of refinemen",            m_level);
        cmd.addInt   ("S", "macroSize",       "size of the refined box in elements",m_boxSize);
        cmd.addInt   ("D", "levelDistance",   "distance between two levels",   m_levSep);
        cmd.addInt   ("e", "minElements",     "number of elements at level 0", m_numElLev0);
        cmd.addSwitch("B", "basisFromDomain", "make start basis from domain",  m_basisFromDomain);

        cmd.getValues(m_argn,m_args);

        if (m_smooth<0)
            m_smooth=m_deg-1;
        else if (m_smooth>=m_deg)
            GISMO_ERROR("smoothness must be smaller than degree");
        if ( !strcmp(methodName.c_str(),"TH") || !strcmp(methodName.c_str(),"TaylorHood" ) )
            m_method=gsStokesRefiner::TH;
        else if ( !strcmp(methodName.c_str(),"SG") || !strcmp(methodName.c_str(),"SubGrid" ) )
            m_method=gsStokesRefiner::SG;
        else if ( !strcmp(methodName.c_str(),"RT") || !strcmp(methodName.c_str(),"RaviartThomas" ) )
            m_method=gsStokesRefiner::RT;
        else GISMO_ERROR("Invalid method name");

        // default depending on degree and smoothness
        m_macroSize   = (m_smooth+m_deg+1)/(m_deg-m_smooth);
        m_boxSize     = m_boxSize>=0 ? m_boxSize : m_macroSize;
        m_levSep = m_levSep>=0 ? m_levSep : m_macroSize;
    }
};


std::vector<gsBasis<>*> makeInitialBasis(const gsMultiPatch<> &domain, const InfSupCommonOptions & opt)
{
    std::vector<gsBasis<>*> result(domain.size());
    if (opt.m_basisFromDomain)
    {
        for (int p=0; p<domain.size();++p)
            result[p]=&domain.basis(p);
    }
    else
    {
        gsKnotVector<> knots(0,1,opt.m_numElLev0-1, opt.m_deg+1, opt.m_deg-opt.m_smooth);
        gsTensorBSplineBasis<2> temp(knots,knots);
        std::fill(result.begin(),result.end(),new gsHBSplineBasis<2>(temp));
    }
    return result;
}


class InfSupTest : public InfSupCommonOptions, public gsStokesRefiner
{
protected:
    using gsStokesRefiner::m_method;
    using gsStokesRefiner::m_domain;
    using InfSupCommonOptions::m_smooth;
    int                    m_dim;
public:

    InfSupTest(const gsMultiPatch<> &domain, std::string name, int argn, char** args)
        :
          InfSupCommonOptions(name,argn,args),
          gsStokesRefiner(domain, makeInitialBasis(domain, InfSupCommonOptions(name,argn,args)), InfSupCommonOptions::m_method, InfSupCommonOptions::m_smooth)
    {
        m_dim=m_domain.dim();
    }

    virtual ~InfSupTest()
    {}

    const gsMultiPatch<> &domain() const {return m_domain;}

    void init()
    {
        for (size_t patch=0; patch < m_basis.nPatches(); ++patch )// for all patches
        {
            m_basis.refineElements( patch, getPatchBoxes(patch), false );
        }
        m_basis.repairPatches();
        m_basis.updateTopol();
        makeVelocityBasis();
    }

    virtual std::vector<index_t> getPatchBoxes(size_t p) = 0;

    virtual std::vector<index_t>  getEliminated()
    {
        std::vector<gsPhysicalSpace*> m_space=getSpaces();

        std::vector<index_t> result;
        gsMultiPatch<>::const_biterator iter=m_domain.bBegin();
        for ( ;iter!=m_domain.bEnd(); ++iter)
        {
            std::vector<index_t> patchDofs = m_space[velocity]->boundaryDofs(*iter,0);
            for (size_t i=0; i<patchDofs.size();++i )
                result.erase(std::remove(result.begin(),result.end(),patchDofs[i]),result.end());
            result.insert(result.end(),patchDofs.begin(),patchDofs.end()) ;
        }
        freeAll(m_space);
        return result;
    }


    void printMeshes(std::string fileName) const
    {
        gsParaviewCollection colP(fileName+"_Pmesh");
        gsParaviewCollection colV(fileName+"_Vmesh");
        for (int i=0; i<m_domain.size();++i)
        {
            std::stringstream nameP;
            nameP<<fileName<<"_P"<<i<<"_Pmesh";
            gsMesh<real_t> mshP;
            makeMesh<real_t>(m_basis.basis(i), mshP,8);
            m_domain[i].evaluateMesh(mshP);
            gsWriteParaview(mshP, nameP.str(), false);
            colP.addPart(nameP.str(),".vtp");

            std::stringstream nameV;
            nameV<<fileName<<"_P"<<i<<"_Vmesh";
            gsMesh<real_t> mshV;
            makeMesh<real_t>(m_basisV[0]->basis(i), mshV,8);
            m_domain[i].evaluateMesh(mshV);
            gsWriteParaview(mshV, nameV.str(), false);
            colV.addPart(nameV.str(),".vtp");
        }
        colP.save();
        colV.save();
    }

    virtual std::string describe()
    {
        static const char* methodName[]={"NONE","TH","SG","RT"};
        std::ostringstream out;
        out<<m_name<<"_"<<methodName[m_method]<<"d"<<m_deg<<"s"<<m_smooth<<"l"<<m_level;
        return out.str();
    }
};


class gsMatricesAssembler : public gsRecipeAssembler
{
protected:
    gsWeightMapper<real_t> *map[2];
public:
    gsMatricesAssembler(const gsMultiPatch<real_t> &domain)
        : gsRecipeAssembler(domain)
    {
        map[0]=NULL;
        map[1]=NULL;
        m_shiftsSource.resize(2);
        m_shiftsSource[0]=0;
        m_shiftsSource[1]=0;
    }

    ~gsMatricesAssembler()
    {
        freeAll(map, map+2);
    }

    struct TestMatrices
    {
        gsSparseMatrix<>   norm1Mat;
        gsSparseMatrix<>   norm2Mat;
        gsSparseMatrix<>   operatorMat;
        //        gsMatrix<unsigned> dirichletDofs;  // not needed as we can delete dofs here
    };

    TestMatrices m_matrices;

    const TestMatrices& matrices() const {return m_matrices;}

    void initMappers()
    {
        for (int i=0; i<2;++i)
        {
            map[i]=m_space[i]->getMapper();
            map[i]->optimize();
        }
        std::sort(m_eliminatedDofs[0].begin(),m_eliminatedDofs[0].end());
        map[0]->sourceToTarget(m_eliminatedDofs[0], m_eliminatedTarget);
        reorderMapperTarget(*map[0],m_eliminatedTarget,&m_permutation);
        map[0]->optimize();
        map[1]->optimize();
    }

    void initSystemMatrices()
    {
        const index_t Vsize=getFreeLimit();
        const index_t Psize=map[1]->getNrOfTargets();

        m_matrices.norm1Mat.resize(Vsize,Vsize);
        m_matrices.norm2Mat.resize(Psize,Psize);
        m_matrices.operatorMat.resize(Psize,Vsize);

        index_t estimatedOverlap=0;
        for (index_t patch=0; patch<m_domain.size();++patch)
            estimatedOverlap=std::max<index_t>(getPatchMaxOverlap(patch),estimatedOverlap);

        gsVector<index_t> perColReserve;
        perColReserve.setConstant(Vsize,math::min(estimatedOverlap,Vsize));
        m_matrices.norm1Mat.reserve(perColReserve);

        perColReserve.setConstant(Psize,math::min(estimatedOverlap,Psize));
        m_matrices.norm2Mat.reserve(perColReserve);
        perColReserve.setConstant(Vsize,math::min(estimatedOverlap,Psize));
        m_matrices.operatorMat.reserve(perColReserve);
    }

void postProcess()
{
    m_matrices.operatorMat.prune( (real_t)(0.1));
    m_matrices.norm1Mat.prune((real_t)(0.1));
    m_matrices.norm2Mat.prune((real_t)(0.1));

    m_matrices.operatorMat.makeCompressed();
    m_matrices.norm1Mat.makeCompressed();
    m_matrices.norm2Mat.makeCompressed();
}

    void init()
    {
        initMappers();
        initSystemMatrices();
    }

    virtual index_t getFreeLimit () const  {return map[0]->getNrOfTargets()-m_eliminatedTarget.size();}

    virtual gsRecipe<real_t>    getBoundaryRecipe  (patchSide s)
    {
        return gsRecipe<real_t>();
    }

    virtual gsRecipe<real_t> getPatchRecipe(index_t patch)
    {
        typedef gsMatAndRhsModWriter<SMatT,gsNullWriter<real_t> > MyW;

        gsRecipe<real_t> result;
        gsRecipeIngredient<real_t> ingr;

        ingr.setOperator(new gsGradGradOp<real_t>());
        ingr.setTestSpace(velocity);
        ingr.setUnknownSpace(velocity);
        ingr.setRule(new gsL2GMapper<MyW>(
                         *map[0],
                     *map[0],
                MyW(getFreeLimit(),getFreeLimit(),m_matrices.norm1Mat,gsNullWriter<real_t>::unique_instance)
                ));
        result.add(ingr);

        ingr.setOperator(new gsDivergenceOp<real_t>());
        ingr.setTestSpace( pressure);
        ingr.setUnknownSpace( velocity);
        ingr.setRule(new gsL2GMapper<MyW>(
                         *map[0],
                     *map[1],
                MyW(getFreeLimit(),getFreeLimit(),m_matrices.operatorMat,gsNullWriter<real_t>::unique_instance)
                ));
        result.add(ingr);


        ingr.setOperator(new gsL2ScalarOp<real_t>());
        ingr.setTestSpace(pressure);
        ingr.setUnknownSpace(pressure);
        ingr.setRule(new gsL2GMapper<MyW>(
                         *map[1],
                     *map[1],
                MyW(getFreeLimit(),getFreeLimit(),m_matrices.norm2Mat,gsNullWriter<real_t>::unique_instance)
                ));
        result.add(ingr);

        return result;
    }
};


// FACTORY METHOD


typedef InfSupTest* (*InfSupTestAllocator) (std::string name, int argn, char** args);

struct InfSupTestFactory
{
protected:
    std::map<std::string, InfSupTestAllocator> m_registry;
public:
    static InfSupTestFactory& getDefaultFactory() {static InfSupTestFactory f; return f;}

    static char addToDefaultFactory(InfSupTestAllocator allocator, std::string key)
    {
        return getDefaultFactory().add(allocator,key);
    }

    void printAll()
    {
        std::map<std::string, InfSupTestAllocator>::const_iterator it=m_registry.begin();
        std::map<std::string, InfSupTestAllocator>::const_iterator end=m_registry.end();
        for ( ; it!=end; ++it)
        {
            gsInfo<<it->first<<" "<<(it->second)<<"\n";
        }
    }

    bool isRegistered (const std::string &name)
    {
        std::map<std::string, InfSupTestAllocator>::const_iterator it, end;
        it = m_registry.find(name);
        end= m_registry.end();
        return it!=end;
    }

    InfSupTest* build (const std::string &name, int argn, char** args )
    {
        std::map<std::string, InfSupTestAllocator>::const_iterator it, end;
        it = m_registry.find(name);
        end= m_registry.end();
        if (it!=end)
            return it->second(it->first,argn,args);
        else return NULL;
    }

    char add(InfSupTestAllocator allocator, std::string key)
    {
        m_registry.insert(std::make_pair(key,allocator));
        return 1;
    }

};
template<typename T>
struct RegInfSupTest
{
    static char __dummy_infsup_test;
    static InfSupTest *build (std::string name, int argn, char** args)
    {
        return new T(name,argn, args);
    }
};
#define infsupTestRegister(A,B)        template <> char RegInfSupTest<A>::__dummy_infsup_test = InfSupTestFactory::addToDefaultFactory(RegInfSupTest<A>::build, B)
#define infsupTestRegisterOther(A,B,C) template <> char RegInfSupTest<A>::__dummy_infsup_test = InfSupTestFactory::addToDefaultFactory(C, B)


// preallocated domains
gsMultiPatch<> *domains[]=
{
    new gsMultiPatch<>(*gsNurbsCreator<real_t>::BSplineSquare()),
    ((gsMultiPatch<>::uPtr)gsReadFile<>("planar/thbs_multipatch_01.xml")).release(),
    new gsMultiPatch<real_t>(gsNurbsCreator<real_t>::BSplineSquareGrid(1, 1, 1)),
    new gsMultiPatch<real_t>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus()),
    ((gsMultiPatch<>::uPtr)gsReadFile<>("planar/thbPuzzlePiece.xml")).release()
};

// MAIN FUNCTION

int main(int argn, char *args[])
{
    // parse arguments
    static const char  *endArgs[]={"--", NULL};

    std::string testName="corner";
    std::string exportName="";
    char** newArgs = std::find_first_of(args,args+argn,endArgs,endArgs+1,mstrcmp);
    int    ourArgn = newArgs-args;
    int    newArgn = argn-ourArgn;
    bool   printHead=false;
    bool   meshFile =false;
    bool   avoidComp=false;

    gsCmdLine cmd  ("Assemble matrices required for infsup test for the divergence operator");
    cmd.addString  ("t", "test", "Input test name",             testName);
    cmd.addSwitch  ("H", "head", "print header in time table",  printHead);
    cmd.addSwitch  ("M", "mesh", "save mesh file",              meshFile);
    cmd.addSwitch  ("A", "dont", "do not compute",              avoidComp);
    cmd.addString  ("o", "out",  "outfile name",                exportName);


    try { cmd.getValues(ourArgn,args); } catch (int rv) { return rv; }

    InfSupTestFactory F=InfSupTestFactory::getDefaultFactory();
    if (! F.isRegistered(testName))
    {
        gsInfo<<"\nThe registered tests are:\n";
        F.printAll();
        gsInfo<<"\n";
        exit(EXIT_FAILURE);
    }

    double initTime, assembleTime, exportTime, totalTime;
    gsStopwatch clock;

    InfSupTest *mytest=F.build(testName,newArgn,newArgs);
    mytest->init();
    std::vector<index_t> eliminated=mytest->getEliminated();
    initTime     = clock.stop(); clock.restart();

    // compute
    gsMatricesAssembler assembler(mytest->domain());
    if (!avoidComp)
    {
        assembler.setSpace(mytest->getSpaces());
        assembler.eliminateDofs(eliminated,velocity);
        assembler.assemble();
    }
    assembleTime = clock.stop(); clock.restart();

    // export
    if(exportName.length()==0)
        exportName=mytest->describe();
    if (meshFile)
        mytest->printMeshes(exportName);
    if (!avoidComp)
    {
        exportMatrixToASCII(exportName+"_vN.dat",  assembler.matrices().norm1Mat);
        exportMatrixToASCII(exportName+"_pN.dat",  assembler.matrices().norm2Mat);
        exportMatrixToASCII(exportName+"_div.dat", assembler.matrices().operatorMat);
    }
    exportTime = clock.stop();
    totalTime  = exportTime + assembleTime + initTime;

    // output times
    gsInfo.setf(std::ios::left, std::ios::floatfield);
    gsInfo<<std::scientific;
    gsInfo.precision(3);
    if(printHead)
        gsInfo << std::setw(exportName.length())<<"test"
                  << std::setw(12)<<"init"
                  << std::setw(12)<<"assemble"
                  << std::setw(12)<<"export"
                  << std::setw(12)<<"total"
                  << "\n";
    gsInfo<<std::setw(exportName.length())<<exportName
            <<std::setw(12)<<initTime
           <<std::setw(12)<<assembleTime
          <<std::setw(12)<<exportTime
         <<std::setw(12)<<totalTime<<"\n";


    delete mytest;
    freeAll(domains,domains+sizeof(domains)/sizeof(gsMultiPatch<> *));
    return 0;
}



// test code

class StokesCornerRefinement : public  InfSupTest
{
public:
    StokesCornerRefinement(std::string name, int argn, char** args)
        : InfSupTest(*domains[0],name,argn,args)
    {}

protected:
    virtual std::vector<index_t> getPatchBoxes(size_t p)
    {
        std::vector<index_t> result;
        result.reserve( (2*m_dim+1)*m_level);
        for (int l=1; l<=m_level; ++l)
        {
            const int beg = 0;
            const int end = math::min(beg+m_boxSize, (m_numElLev0)*(1<<l));

            result.push_back(l);
            for ( int s = 0; s< m_dim; ++s)
                result.push_back(beg);
            for ( int s = 0; s< m_dim; ++s)
                result.push_back(end);
        }
        return result;
    }
};
infsupTestRegister(StokesCornerRefinement, "corner");

class StokesShiftedCornerRefinement : public InfSupTest
{
public:
    StokesShiftedCornerRefinement(std::string name, int argn, char** args)
        : InfSupTest(*domains[0],name,argn,args)
    {}

    virtual std::vector<index_t> getPatchBoxes(size_t p)
    {
        std::vector<index_t> result;
        result.reserve( (2*m_dim+1)*m_level);
        for (int l=1; l<=m_level; ++l)
        {
            const int beg = (2<<l);
            const int end = math::min(beg+m_boxSize, (m_numElLev0)*(1<<l));

            result.push_back(l);
            for ( int s = 0; s< m_dim; ++s)
                result.push_back(beg);
            for ( int s = 0; s< m_dim; ++s)
                result.push_back(end);
        }
        return result;
    }
};
infsupTestRegister(StokesShiftedCornerRefinement, "scorner");



class StokesTP : public InfSupTest
{
public:
    StokesTP(std::string name, int argn, char** args)
        : InfSupTest(*domains[0],name,argn,args)
    {}

    virtual std::vector<index_t> getPatchBoxes(size_t p)
    {
        std::vector<index_t> result;
        return result;
    }
};
infsupTestRegister(StokesTP, "tensor");


class StokesMacro : public StokesTP
{
public:
    StokesMacro(std::string name, int argn, char** args)
        : StokesTP(name,argn,args)
    {
        m_numElLev0=m_macroSize;
    }
    void catVector(gsMatrix<unsigned> &a, gsMatrix<unsigned> &b)
    {
        a.conservativeResize(a.rows()+b.rows(),1);
        a.bottomRows(b.rows())=b;
    }

    virtual std::vector<index_t>  getEliminated()
    {
        std::vector<gsPhysicalSpace*> m_space=getSpaces();

        std::vector<index_t> tmp, result;
        gsMultiPatch<>::const_biterator iter=m_domain.bBegin();
        for ( ;iter!=m_domain.bEnd(); ++iter)
        {
            std::vector<index_t> patchDofs = m_space[velocity]->boundaryDofs(*iter,m_smooth);
            for (size_t i=0; i<patchDofs.size();++i )
                tmp.erase(std::remove(tmp.begin(),tmp.end(),patchDofs[i]),tmp.end());
            tmp.insert(tmp.end(),patchDofs.begin(),patchDofs.end()) ;
        }
        freeAll(m_space);
        return result;
    }
};
infsupTestRegister(StokesMacro, "macro");



class StokesCenterRefinement : public InfSupTest
{
public:
    StokesCenterRefinement(std::string name, int argn, char** args)
        : InfSupTest(*domains[0],name,argn,args)
    {
        // change number of elements at level 0 to make space for m_level different levels
        std::vector<index_t> boxes=StokesCenterRefinement::getPatchBoxes(0);
        m_numElLev0 = (boxes[m_dim+2]+1)/2+m_levSep;
        initBasis(makeInitialBasis(domain(),*this));
    }

//    virtual std::vector<index_t>  getPatchBoxes(size_t p)
//    {
//        std::vector<index_t> result;
//        result.reserve( (2*m_dim+1)*m_level);
//        for (int l=1; l<=m_level; ++l)
//        {
//            const int mid  = (m_minElements)*static_cast<int>(1<<(l-1));
//            const int half = (m_boxSize+1)/2;
//            const int beg  = math::max(mid-half, 0);
//            const int end  = math::min(mid+half, static_cast<int>((1<<l)*(m_minElements)));

//            result.push_back(l);
//            for ( int s = 0; s< m_dim; ++s)
//                result.push_back(beg);
//            for ( int s = 0; s< m_dim; ++s)
//                result.push_back(end);
//        }
//        return result;
//    }

    virtual std::vector<index_t>  getPatchBoxes(size_t p)
    {
        const int bSize=2*m_dim+1;
        std::vector<index_t> result;
        result.resize(bSize*m_level);

        int beg=0;
        for (int l=0; l<m_level;++l)
        {
            beg=(beg+m_levSep)*2;
            result[l*bSize]=l+1;
            for ( int s = 0; s< m_dim; ++s)
                result[l*bSize+s+1]=beg;
        }

        int end=beg+m_boxSize;
        for ( int s = 0; s< m_dim; ++s)
            result[(m_level-1)*bSize+s+1+m_dim]=end;

        for (int l=m_level-2; l>=0;--l)
        {
            end=(end+1)/2+m_levSep;
            for ( int s = 0; s< m_dim; ++s)
                result[l*bSize+s+1+m_dim]=end;
        }

//        // insert max level box
//        if (m_level>0)
//        {
//            const int mid  = (m_numElLev0)*static_cast<int>(1<<(m_level-1));
//            const int half = (m_boxSize+1)/2;
//            beg = math::max(mid-half, 0);
//            end = math::min(mid+half, static_cast<int>((1<<m_level)*(m_numElLev0)));

//            result.push_back(m_level);
//            for ( int s = 0; s< m_dim; ++s)
//                result.push_back(beg);
//            for ( int s = 0; s< m_dim; ++s)
//                result.push_back(end);
//        }
//        // insert lower lvl boxes
//        for (int l=m_level-1; l>0; --l)
//        {
//            beg = math::max(beg/2-m_levSep, 0);
//            end = math::min((end+1)/2+m_levSep,static_cast<int>((1<<l)*(m_numElLev0)));

//            if (end-beg<m_boxSize)
//            {
//                const int move = (m_boxSize-end+beg+1)/2;
//                beg = math::max(beg-move, 0);
//                end = math::min(end+move, static_cast<int>((1<<l)*(m_numElLev0)));
//            }

//            result.push_back(l);
//            for ( int s = 0; s< m_dim; ++s)
//                result.push_back(beg);
//            for ( int s = 0; s< m_dim; ++s)
//                result.push_back(end);
//        }
        return result;
    }
};
infsupTestRegister(StokesCenterRefinement, "center");



class StokesEdgeRefinement : public InfSupTest
{
public:
    StokesEdgeRefinement(std::string name, int argn, char** args)
        : InfSupTest(*domains[0],name,argn,args)
    {}

    virtual unsigned minNumElements() { return m_numElLev0+=m_numElLev0%2; }
    virtual std::vector<index_t>  getPatchBoxes(size_t p)
    {
        const int m_basisSupp=(m_deg+1)/(m_deg-m_smooth);

        const int levBlock = (2*m_dim+1);
        std::vector<index_t> result;
        result.resize( levBlock*m_level);

        int beg, end, dir1beg,dir1end;

        const int mid  = (m_numElLev0)*static_cast<int>(1<<(m_level-1));
        const int half = (m_boxSize+1)/2;
        const int dist = m_levSep;

        beg = math::max<int>(mid-half, 0);
        beg-= beg%2;
        end = math::max<int>(mid+half, 0);
        end+= end%2;
        dir1beg = mid-m_boxSize;
        dir1beg -= dir1beg%2;
        dir1beg = math::max<int>(dir1beg, 0);
        dir1end = mid;

        for (int l=m_level-1; l>=0; --l)
        {
            const int levS=l*levBlock;
            const int levSB=levS+1;
            const int levSE=levSB+m_dim;
            // write box
            result[levS]=l+1;
            result[levSB]=dir1beg;
            for ( int s = 1; s< m_dim; ++s)
                result[levSB+s]=beg;
            result[levSE]=dir1end;
            for ( int s = 1; s< m_dim; ++s)
                result[levSE+s]=end;
            // update box
            beg  = (beg-2*dist)/2;
            beg -= beg%2;
            beg  = math::max<int>(beg, 0);
            end  = math::max<int>((end+2*dist+end%2)/2, beg+m_basisSupp);
            end += end%2;
            end  = math::min<int>(end, (1<<l)*(m_numElLev0));
            dir1end /= 2;
            dir1beg  = math::min<int>((dir1beg-2*dist)/2, dir1end-m_basisSupp);
            dir1beg -= dir1beg%2;
            dir1beg  = math::max<int>(dir1beg, 0);
        }
        return result;
    }
};
infsupTestRegister(StokesEdgeRefinement, "edge");
