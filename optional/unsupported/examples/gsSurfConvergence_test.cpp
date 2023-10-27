
#include <iostream>
#include <vector>

#include <gismo.h>
#include <gismo_dev.h>
#include <gsUtils/gsBVProblem.h>

#include <gsAssembler/gsGaussSurfaceAssembler.h>



using namespace gismo;


gsField<> solveOnSurface( gsBVProblem<> const & bvp,
                          std::vector< gsBasis<>* > & bases, 
                          bool conforming = true, bool elimination = true);

void eliminateDirichletDofs(gsBVProblem<> const & bvp,
                            gsDofMapper const & mapper,
                            gsMatrix<> & ddof,
                            std::vector< gsBasis<>* > const & bases);

//////// L2 - Norm /////////////////

real_t igaSurfL2Distance(const gsGeometry<>& patch, 
                         const gsGeometry<>& func, 
                         const gsFunction<>& v, 
                         bool v_isParam)
{
    gsGeometryEvaluator<>::uPtr geoEval(getEvaluator(NEED_VALUE | NEED_MEASURE , patch));
    // assuming real-valued function
    gsGeometryEvaluator<>::uPtr funcEval(getEvaluator(NEED_VALUE, func));

    // degree of the underlying Gauss rule to use
    gsVector<index_t> numNodes( func.parDim() );
    for (index_t i = 0; i < numNodes.size(); ++i)
        numNodes[i] = func.basis().degree(i) + 1;
 
    real_t sum(0);
    gsBasis<>::domainIter domIt = func.basis().makeDomainIterator();
    gsGaussRule<> quRule(numNodes);
    gsMatrix<> quNodes;
    gsVector<> quWeights;
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    for (; domIt->good(); domIt->next())
    {        
        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(quNodes);
        funcEval->evaluateAt(quNodes);
        const gsMatrix<> & func_vals = funcEval->values();
        
        // Evaluate function v
        gsMatrix<> v_val = v_isParam ? v.eval(quNodes)
            : v.eval( geoEval->values() ); 
        
        // perform the quadrature
        for (index_t k = 0; k < quRule.numNodes(); ++k) // loop over quadrature nodes
        {            
            const gsMatrix<> & Jk = geoEval->jacobians().block(0, k*2, 3, 2);
            gsMatrix<> FirstFund = Jk.transpose()*Jk;
            const real_t weight = quWeights[k] * math::sqrt(FirstFund.determinant());

            const real_t diff = func_vals(0,k) - v_val(0,k);
            sum += weight * diff * diff;
        }
    }
    return math::sqrt(sum);
}

real_t igaSurfFieldL2Distance(const gsField<>& u, const gsFunction<>& v, bool v_isParam = false)
{
    real_t dist(0);

    for (index_t i = 0; i < u.nPatches(); ++i)
    {
        const gsGeometry<> & func  = static_cast<const gsGeometry<> &>( u.function(i) );
        const real_t curDist = igaSurfL2Distance( u.patch(i), func, v, v_isParam);
        dist += curDist * curDist;
    }

    return math::sqrt(dist);
}

real_t igaSurfFieldL2Distance(const gsField<>& u, const gsField<>& v)
{
    real_t dist(0);

    for (index_t i = 0; i < u.nPatches(); ++i)
    {
        const gsGeometry<> & ufunc  = static_cast<const gsGeometry<> &>( u.function(i) );
        const gsGeometry<> & vfunc  = static_cast<const gsGeometry<> &>( v.function(i) );
        const real_t curDist = igaSurfL2Distance( u.patch(i), ufunc, vfunc, v.isParametrized() );
        dist += curDist * curDist;
    }

    return math::sqrt(dist);
}

// H1 norm ///////////////////////////////////////////////////////////////////////////

real_t igaSurfH1Distance(const gsGeometry<>& patch, 
                         const gsGeometry<>& func, 
                         const gsFunction<>& v, 
                         bool v_isParam)
{
    const int d = func.parDim();

    gsGeometryEvaluator<>::uPtr geoEval(getEvaluator(NEED_VALUE |
                                                     NEED_GRAD_TRANSFORM |
                                                     NEED_MEASURE, patch));
    // assuming real-valued function
    gsGeometryEvaluator<>::uPtr funcEval(getEvaluator(NEED_JACOBIAN, func));
    
    // degree of the underlying Gauss rule to use
    gsVector<index_t> numNodes( func.parDim() );
    for (index_t i = 0; i < numNodes.size(); ++i)
        numNodes[i] = func.basis().degree(i) + 1;
 
    real_t sum(0);
    gsBasis<>::domainIter domIt = func.basis().makeDomainIterator();
    gsGaussRule<> quRule(numNodes);
    gsMatrix<> quNodes;
    gsVector<> quWeights;
    quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    for (; domIt->good(); domIt->next())
    {        
        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(quNodes);
        funcEval->evaluateAt(quNodes);
        gsMatrix<> func_ders = funcEval->jacobians();// (!) coping 
        
        // get the gradients to columns
        func_ders.transposeInPlace();
        func_ders.resize(d, quRule.numNodes() );

        // Evaluate function v
        gsMatrix<> v_ders = v_isParam ? v.deriv(quNodes)
            : v.deriv( geoEval->values() ); 

        // get the gradients to columns
        v_ders.transposeInPlace();
        if ( v_isParam )
            v_ders.resize(d, quRule.numNodes() );
        else
            v_ders.resize(d+1, quRule.numNodes() );

        // perform the quadrature
        for (index_t k = 0; k < quRule.numNodes(); ++k) // loop over quadrature nodes
        {
            // Transform the gradients
            const gsMatrix<> & Jk = geoEval->jacobians().block(0, k*2, 3, 2);
            gsMatrix<> FirstFund = Jk.transpose()*Jk;
            const real_t weight = quWeights[k] * math::sqrt(FirstFund.determinant());
            FirstFund  = Jk* FirstFund.inverse();

            gsMatrix<> func_pders =  FirstFund* func_ders.col(k);

            gsMatrix<> v_der_k;
            if ( v_isParam )
                v_der_k =  FirstFund* v_ders.col(k);
            else
                v_der_k = v_ders.col(k);
            
            sum += weight * (func_pders - v_der_k).squaredNorm();
        }
    }
    return math::sqrt(sum);
}

real_t igaSurfFieldH1Distance(const gsField<>& u, const gsFunction<>& v, bool v_isParam = false)
{
    real_t dist(0);

    for (index_t i = 0; i < u.nPatches(); ++i)
    {
        const gsGeometry<> & func  = static_cast<const gsGeometry<> &>( u.function(i) );
        const real_t curDist = igaSurfH1Distance( u.patch(i), func, v, v_isParam);
        dist += curDist * curDist;
    }

    return math::sqrt(dist);
}

real_t igaSurfFieldH1Distance(const gsField<>& u, const gsField<>& v)
{
    real_t dist(0);

    for (index_t i = 0; i < u.nPatches(); ++i)
    {
        const gsGeometry<> & ufunc  = static_cast<const gsGeometry<> &>( u.function(i) );
        const gsGeometry<> & vfunc  = static_cast<const gsGeometry<> &>( v.function(i) );
        const real_t curDist = igaSurfH1Distance( u.patch(i), ufunc, vfunc, v.isParametrized());
        dist += curDist * curDist;
    }

    return math::sqrt(dist);
}

/////// DG norm ///////////////////////////////////////////////////////////

real_t igaSurfDGDistanceJump(const gsGeometry<>& patch1, const gsGeometry<>& patch2,
                             const gsGeometry<>& func1,  const gsGeometry<>& func2, // approximati solution
                             const gsFunction<>& v1, const gsFunction<>& v2,	// exact solution
                             const boundaryInterface & bi, // interface
                             const real_t mu,
                             bool v_isParam)
{
    const int d = func1.parDim();

    gsGeometryEvaluator<>::uPtr geoEval1(getEvaluator(NEED_VALUE | NEED_MEASURE, patch1));
        
    gsGeometryEvaluator<>::uPtr geoEval2(getEvaluator(NEED_VALUE | NEED_MEASURE, patch2));
    // assuming real-valued function
    gsGeometryEvaluator<>::uPtr funcEval1(getEvaluator(NEED_VALUE, func1));
    
    gsGeometryEvaluator<>::uPtr funcEval2(getEvaluator(NEED_VALUE, func2));
    
    const boxSide side1 = bi[0].side();
    const boxSide side2 = bi[1].side();

    // "DG method not implemented yet for non-matching interfaces");
    
    // assumes matching orientation
    // degree of the underlying Gauss rule to use
    gsVector<index_t> intNodes1 ( func1.basis().dim() );
    const int dir1 = bi.first().direction();
    for (int i = 0; i < dir1; ++i)
        intNodes1[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;
    intNodes1[dir1] = 1;
    for (int i = dir1+1; i < func1.basis().dim(); ++i)
        intNodes1[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;    
    
    gsVector<index_t> intNodes2 ( func2.basis().dim() );
    const int dir2 = bi.second().direction() ;
    for (int i = 0; i < dir2; ++i)
        intNodes2[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;
    intNodes2[dir2] = 1;
    for (int i = dir2+1; i < func1.basis().dim(); ++i)
        intNodes2[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;
    
    
    // Temporaries
    // gsMatrix<> grads_k_1, grads_k_2;!! (not needed)
    gsVector<> unormal(d);
    
    
    real_t sum(0);
    // iterator on grid cells on the "right"
    gsDomainIterator<>::uPtr domIter2= func2.basis().makeDomainIterator(side2);


    //int count = 0;
    // iterate over all boundary grid cells on the "left"
    for (gsDomainIterator<>::uPtr domIter1 = func1.basis().makeDomainIterator(side1); 
         domIter1->good(); domIter1->next())
    {
        // Compute the quadrature rule on both sides
        gsGaussRule<> quRule1( intNodes1 );
        gsGaussRule<> quRule2( intNodes2 );

        gsMatrix<> quNodes1;
        gsVector<> quWeights1;
        quRule1.mapTo(domIter1->lowerCorner(), domIter1->upperCorner(), quNodes1, quWeights1);

        gsMatrix<> quNodes2;
        gsVector<> quWeights2;
        quRule2.mapTo(domIter2->lowerCorner(), domIter2->upperCorner(), quNodes2, quWeights2);

        // compute image of Gauss nodes under geometry mapping as well
        // as Jacobians
        geoEval1->evaluateAt(quNodes1);
        geoEval2->evaluateAt(quNodes2);
	
        funcEval1->evaluateAt(quNodes1);
        funcEval2->evaluateAt(quNodes2);
	
        gsMatrix<> func1_vals = funcEval1->values();// (!) coping 
        gsMatrix<> func2_vals = funcEval2->values();// (!) coping 
	
        // exact solution
        gsMatrix<> v1_vals = v_isParam ? v1.eval(quNodes1)
            : v1.eval( geoEval1->values() ); 
			    
        gsMatrix<> v2_vals = v_isParam ? v2.eval(quNodes2)
            : v2.eval( geoEval2->values() ); 

        for (index_t k=0; k!= quRule1.numNodes(); ++k)
        {            
            const gsMatrix<> & jac1 = geoEval1->jacobians().block(0, k*d, d+1, d);

            // Compute the tangent vector from patch1
            const gsVector<real_t,3> tangent  =  jac1.block<3, 1>(0,!dir1);

            // Integral transformation and quadarature weight (patch1)
            // assumed the same on both sides
            const real_t fff = mu * quWeights1[k] * tangent.norm();
   
            const real_t diff = func1_vals(0,k) - v1_vals(0,k) 
                - func2_vals(0,k) + v2_vals(0,k);
            sum += fff * diff*diff;
        }

        domIter2->next();
    }
    return math::sqrt(sum);
}

/////////
real_t igaSurfFieldDGDistance(const gsField<>& u, const gsFunction<>& v, bool v_isParam= false)
{
    real_t dist = igaSurfFieldH1Distance(u, v, v_isParam);
    dist *= dist; // square
 
    gsMultiPatch<> mp = u.patches();
    // iterate over all interfaces
    for ( gsMultiPatch<>::const_iiterator it = mp.iBegin();
          it != mp.iEnd(); ++it ) // *it ---> interface
    {
	  
	  
        const gsGeometry<> & func1  = static_cast<const gsGeometry<> &>( u.function(it->first().patch) );
        const gsGeometry<> & func2  = static_cast<const gsGeometry<> &>( u.function(it->second().patch) );

        // get penalty parametr
        const real_t h = math::pow( (real_t) func1.basis().size(), -1.0 / func1.basis().dim() );
        const real_t bdeg = (real_t)func1.basis().degree(0);
        real_t mu = ( (bdeg+func1.basis().dim())* (bdeg+1) * 2.0 / h );
        const real_t curDist = igaSurfDGDistanceJump( mp.patch(it->first().patch),
                                                      mp.patch(it->second().patch),
                                                      func1,
                                                      func2,
                                                      v, v,
                                                      *it,
                                                      mu, // mu 
                                                      v_isParam);
	    dist += curDist * curDist;
	}
    return math::sqrt(dist);
}

real_t igaSurfFieldDGDistance(const gsField<>& u, const gsField<>& v)
{
    real_t dist = igaSurfFieldH1Distance(u, v);
    dist *= dist; // square
    bool v_isParam = v.isParametrized() ;
    
    gsMultiPatch<> mp = u.patches();
    // iterate over all interfaces
    for ( gsMultiPatch<>::const_iiterator it = mp.iBegin();
          it != mp.iEnd(); ++it ) // *it ---> interface
    {
	    
        const gsGeometry<> & ufunc1  = static_cast<const gsGeometry<> &>( u.function(it->first().patch) );
        const gsGeometry<> & ufunc2  = static_cast<const gsGeometry<> &>( u.function(it->second().patch) );
	  
        const gsGeometry<> & vfunc1  = static_cast<const gsGeometry<> &>( v.function(it->first().patch) );
        const gsGeometry<> & vfunc2  = static_cast<const gsGeometry<> &>( v.function(it->second().patch) );

        // get penalty parametr
        const real_t h = math::pow( (real_t) ufunc1.basis().size(), -1.0 / ufunc1.basis().dim() );
        const real_t bdeg = (real_t)ufunc1.basis().degree(0);
        real_t mu = ( (bdeg+ufunc1.basis().dim())* (bdeg+1) * 2.0 / h );
        const real_t curDist = igaSurfDGDistanceJump( mp.patch(it->first().patch),
                                                      mp.patch(it->second().patch),
                                                      ufunc1,
                                                      ufunc2,
                                                      vfunc1, vfunc2,
                                                      *it,
                                                      mu, // mu 
                                                      v_isParam);
	    dist += curDist * curDist;
	}
    return math::sqrt(dist);
}
//////
int main(int argc, char *argv[])
{   
    // Input options
    index_t numElevate  = 0;
    index_t numHref     = 0;
    index_t basisDegree = 0;
    bool plot       = false;
    bool compareFinest = false;

    // Multipatch object
    gsMultiPatch<> mp;

    memory::unique_ptr< gsPde<> > pde;

    int result = 0;
    std::string fn("surfaces/quartercylinder2.xml");
    std::string fn_pde("pde/surfacepoisson_quartercylinder.xml");
    
    gsCmdLine cmd("Testing a multipatch problem.");
    cmd.addInt("r","hRefine", 
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("p","degree", 
               "Degree of the basis functions to use for solving (will elevate or reduce the input)",
               basisDegree);
    cmd.addString("q","pde","File containing a poisson PDE (.xml)", fn_pde);
    cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)", fn);
    cmd.addInt("e","degreeElevation", 
               "Number of degree elevation steps to perform on the Geometry's basis before solving", 
               numElevate);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("finest", "Compare against finest level solution for conv. computation.", compareFinest);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    
    gsReadFile<>(fn, mp);

    pde = gsReadFile<>(fn_pde);
    gsFunctionExpr<>::Ptr solution = gsReadFile<>(fn_pde);

    GISMO_ASSERT( mp.parDim()==2 && mp.geoDim()==3, "surface only");

    // Create BVP
    gsBVProblem<> bvp( mp, pde.release() );

    gsFunctionExpr<> n_bc("0",2);
    // Create Dirichlet boundary conditions for all boundaries
    for (gsMultiPatch<>::const_biterator 
             bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        // No BCs if we add a mass term (hom. Neumann problem)
        bvp.addCondition( *bit, condition_type::dirichlet, solution.get() );
//        bvp.addCondition( *bit, condition_type::neumann, &n_bc );
    }
    
    // Set up and elevate discretization bases
    std::vector< gsBasis<>* > bases;
    for ( gsMultiPatch<>::const_iterator it = mp.begin();
          it != mp.end(); ++it ) 
        //bases.push_back( it->basis().makeNonRational() );
        bases.push_back( (*it)->basis().makeNonRational().release() );
    
    if (basisDegree)
        for (size_t i = 0; i < mp.nPatches(); ++i)
            bases[i]->setDegree(basisDegree);
    else if (numElevate)
        for (size_t i = 0; i < mp.nPatches(); ++i)
            bases[i]->degreeElevate(numElevate);
    
    gsInfo<<"Solving "<< bvp << "\n";  

    // Run tests
    gsMatrix<> testsL2(numHref+1,7);
    testsL2.setZero();

    gsMatrix<> testsH1(numHref+1,7);
    testsH1.setZero();


    std::vector<gsField<> *> sols, sols_weakBC, sols_DG;

    int i = 0;
    do
    {
        int dofs = 0;
        for (size_t j = 0; j < bases.size(); ++j)
            dofs += bases[j]->size();

        //gsInfo<<"Discr. Space for patch 0: \n"<< *bases[0] << "\n";
        gsInfo<< "Computing Conforming solution..\n";// conforming, eliminate
        gsField<> sol = solveOnSurface(bvp, bases, true, true);
        sols.push_back(new gsField<>(sol));

        gsField<> sol_weakBC;
        if ( ! mp.isClosed() )
        {
            gsInfo<< "Computing Conforming solution, weak BCs..\n";// conforming, !eliminate
            sol_weakBC = solveOnSurface(bvp, bases, true, false);
            sols_weakBC.push_back(new gsField<>(sol_weakBC));
        }

        gsInfo<< "Computing DG solution..\n";// !conforming, !eliminate
        gsField<> sol_dg = solveOnSurface(bvp, bases, false, false);
        sols_DG.push_back(new gsField<>(sol_dg));

        // Collect data
        testsL2(i,0) = dofs;
        testsH1(i,0) = dofs;
                
        for (size_t j = 0; j < bases.size(); ++j )
            bases[j]->uniformRefine();
    } 
    while ( i++ < numHref );

    i = 0;

    // Construct the exact solution.. to be improved (it may not be ther)
    gsField<> exact0( mp, *solution, false ) ;

    gsField<> * exact;

    if ( compareFinest || (!solution) )
    {
        exact = sols.back();
    }
    else
    {
        exact = &exact0;
    }

    do
    {
        testsL2(i,1)= igaSurfFieldL2Distance(*sols[i], *exact );
        
        if ( sols_weakBC.size() )
            testsL2(i,3)= igaSurfFieldL2Distance(*sols_weakBC[i], *exact );
        testsL2(i,5)= igaSurfFieldL2Distance(*sols_DG[i], *exact );    
        testsH1(i,1) = igaSurfFieldH1Distance(*sols[i], *exact );
        
        if ( sols_weakBC.size() )
            testsH1(i,3)= igaSurfFieldH1Distance(*sols_weakBC[i], *exact );
        testsH1(i,5)= igaSurfFieldDGDistance(*sols_DG[i], *exact );
//         testsH1(i,5)=-1;
            
        if (i > 0)
        {
            testsL2(i,2) = testsL2(i-1,1) / testsL2(i,1);

            if ( sols_weakBC.size() )
                testsL2(i,4) = testsL2(i-1,3) / testsL2(i,3);
            testsL2(i,6)= testsL2(i-1,5) / testsL2(i,5);
            testsH1(i,2) = testsH1(i-1,1) / testsH1(i,1);
            
            if ( sols_weakBC.size() )
                testsH1(i,4)= testsH1(i-1,3) / testsH1(i,3);
            
            testsH1(i,6)= testsH1(i-1,5) / testsH1(i,5);
        }
        
        if (plot && i== numHref)
        {
            // Write approximate and exact solution to paraview files
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>( *sols.back(), "poisson_problem", 1000);
            
            gsWriteParaview<>( *exact, "poisson_problem_exact", 1000);

            // Run paraview
            result = system("paraview poisson_problem.pvd &");
        }
    } 
    while ( i++ < numHref );
    

    freeAll(sols);
    freeAll(sols_DG);
    if ( sols_weakBC.size() )
        freeAll(sols_weakBC);


    for(i = 1; i<= numHref; ++i)
    {   // Compute convergence rates
        testsL2(i,2)= math::log(testsL2(i,2))/std::log(2.0);
        testsL2(i,4)= math::log(testsL2(i,4))/std::log(2.0);
        testsL2(i,6)= math::log(testsL2(i,6))/std::log(2.0);

        testsH1(i,2)= math::log(testsH1(i,2))/std::log(2.0);
        testsH1(i,4)= math::log(testsH1(i,4))/std::log(2.0);
        testsH1(i,6)= math::log(testsH1(i,6))/std::log(2.0);
    }
    
    if ( mp.isClosed() )
    {
        testsL2.col(3).setZero();
        testsL2.col(4).setZero();
        testsH1.col(3).setZero();
        testsH1.col(4).setZero();
    }

    gsInfo << "Summary:\n\n";
    gsInfo << " (deg= "<<bases[0]->minDegree() <<")  |  CONFORMING (elim. BCs) |  CONFORMING (weak BCs) |   Disc. Galerkin      \n";
    gsInfo << "    --------------------------------------------------------------------------------------\n";
    gsInfo << "    Dofs   |  L2 error | conv. rate|  L2 error |conv. rate |  L2 error | conv. rate    \n" 
           << testsL2  << "\n";

    gsInfo << "    --------------------------------------------------------------------------------------\n";
    gsInfo << "    Dofs   |  H1 error  | conv. rate|  H1 error | conv. rate|  DG error | conv. rate    \n" 
           << testsH1  << "\n";

    if ( mp.isClosed() )
        gsInfo << "Note: surface is closed, no boundary conditions applied.\n";

    freeAll(bases);

    return result;
}


gsField<> solveOnSurface( gsBVProblem<> const & bvp,
                          std::vector< gsBasis<>* > & bases,
                          bool conforming, bool elimination)
{
    // Make assembler object
    gsGaussSurfaceAssembler<real_t> assembler2;

    // Setup global mapping    gsDofMapper dofmapper( bases );

    gsBoundaryConditions<>  empty_bc;
    gsBoxTopology            empty_topology;
    const gsBoundaryConditions<> *bc      =&empty_bc;
    const gsBoxTopology           *topology=&empty_topology;
    // set interfaces to glue to the patch interfaces
    if ( conforming )
        topology = &(bvp.patches());
    // set boundary dofs to eliminate corresponding to Dirichlet sides
    if ( elimination )
        bc = &(bvp.boundaryConditions());

    std::vector< gsBasis<>* > temp;
    temp.resize(bases.size());
    cloneAll(bases.begin(),bases.end(),temp.begin());
    gsDofMapper dofmapper( gsMultiBasis<>(temp, *topology), *bc );
    dofmapper.finalize();

    assembler2.mapper = dofmapper;

    // If the Dirichlet strategy is elimination then precompute
    // Dirichlet dofs
    gsMatrix<> ddof;
    if ( elimination )
    {
        eliminateDirichletDofs(bvp, dofmapper, ddof,bases);
    }

    // Setup the global system
    const int dofs = dofmapper.freeSize();

    gsSparseMatrix<> * K = new gsSparseMatrix<>(dofs,dofs) ;
    gsVector<>       * b = new gsVector<>(dofs);
    b->setZero();

    for (index_t np=0; np < bvp.nPatches(); ++np )
    {
        assembler2.setGeometry( bvp.patch(np) );

        // assemble stiffness matrix and rhs for the local patch np
        gsSparseSystem0<real_t> patchSystem =
            assembler2.assemble( *bases[np], ddof, bvp.pde(), np );

        // add result to the global system (K,b)
        *K += *patchSystem.matrix();
        *b += *patchSystem.rhs();

        patchSystem.free();
    }

    K->makeCompressed();
    gsSparseSystem0<real_t> system(K, b);

    // Enforce Dirichlet boundary conditions by Nitsche's method
    if ( ! elimination )
    {
        for ( gsBVProblem<>::const_iterator 
                  it = bvp.dirichletBegin(); it != bvp.dirichletEnd(); ++it )
        {
            assembler2.setGeometry( bvp.patch(it->patch()) );
            assembler2.applyBoundary( *bases[it->patch()], *it, assembler2.mapper, system );
        }
    }

    // Enforce Neumann boundary conditions
    for (gsBVProblem<>::const_iterator it= bvp.neumannBegin();
         it != bvp.neumannEnd(); ++it)
    {
        assembler2.setGeometry( bvp.patch(it->patch()) );
        assembler2.applyBoundary( *bases[it->patch()], *it,
                                  assembler2.mapper, system );
    }

    // If we are in in dg (Discontinous Galerkin) mode add interface
    // contributions
    if ( ! conforming )
    {
        for ( gsMultiPatch<>::const_iiterator it = 
                  bvp.patches().iBegin();
              it != bvp.patches().iEnd(); ++it )
        {
            assembler2.setGeometry( bvp.patch( it->first().patch ) );
            assembler2.applyDG( *bases[it->first().patch], *bases[it->second().patch],
                                bvp.patch( it->second().patch ),
                                *it, assembler2.mapper, system );
        }
    }

    gsInfo<<"Total DoFs              : "<< dofmapper.freeSize() << "\n";  
    gsInfo<<"Penalty constant        : "<< assembler2.getMu(*bases[0]) << "\n";  
    gsInfo<<"Gauss nodes per direction: "<< 
        assembler2.getNumIntNodesFor((*bases[0])).transpose() << "\n";

    // Solve linear system
    gsVector<> res;


    gsSparseSolver<>::CGDiagonal solver;
    res = solver.compute( * system.matrix() ).solve ( * system.rhs() );
    gsInfo << "residual error: " << solver.error() << "\n";
    gsInfo << "    iterations: " << solver.iterations() << "\n";


    // gsSparseSolver<>::QR 
    //     solver( *system.matrix());
    // res = solver.solve ( * system.rhs() );

    // gsSparseSolver<>::BiCGSTABILUT solver( system.matrix()->selfadjointView<Eigen::Lower>() );
    // res = solver.solve ( * system.rhs() );



    gsInfo << "---------------------------------------\n";

    system.free();

    // Reconstruct patch-wise solutions
    gsPiecewiseFunction<> * sols = new gsPiecewiseFunction<>;
    
    const gsMultiPatch<> & mp = bvp.patches();
   
    for (size_t np=0; np < mp.nPatches(); ++np )
    {    
        const int sz  = bases[np]->size();
        const int fsz = dofmapper.freeSize();
        const int dim = bvp.pde().fieldDim();
        
        gsMatrix<> coeffs(sz, dim);
        for (index_t i = 0; i < sz; ++i)
        {
            if ( dofmapper.is_free(i, np) ) // internal or interface
            {
                for (int k = 0; k < dim; ++k)
                    coeffs(i,k) = res[ k * fsz + dofmapper.index(i, np) ];
            }
            else // eliminated Dof: fill with Dirichlet data
            {
                coeffs.row(i) = ddof.row( dofmapper.bindex(i, np) );
            }
        }
        
        sols->addPiecePointer( bases[np]->makeGeometry( coeffs ) );
    }
    
    return gsField<>(mp, gsPiecewiseFunction<>::Ptr(sols), true);
}



void eliminateDirichletDofs(gsBVProblem<> const & bvp,
                            gsDofMapper const & mapper,
                            gsMatrix<> & ddof,
                            std::vector< gsBasis<>* > const & bases)
{
    ddof.resize( mapper.boundarySize(), bvp.pde().fieldDim() );

    for ( gsBVProblem<>::const_iterator 
              it = bvp.dirichletBegin(); it != bvp.dirichletEnd(); ++it )
    {
        // Get dofs on this boundary
        gsMatrix<index_t> boundary = bases[it->patch()]->boundary(it->side()) ;

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t k=0; k!= boundary.size(); ++k)
            {
                const int ii= mapper.bindex( (boundary)(k) , it->patch() );
                ddof.row(ii).setZero();
            }
            continue;
        }

        // Get the side information
        int dir =  it->side().direction();
        index_t param = ( it->side().parameter()  ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<> > rr;
        rr.reserve( bvp.patches().parDim() );

        for ( int i=0; i < bvp.patches().parDim(); ++i)
        {
            if ( i==dir )
            {
                gsVector<> b(1); 
                b[0] = ( bases[it->patch()]->component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {   
                rr.push_back( bases[it->patch()]->component(i).anchors().transpose() );
            }
        }

        // Compute dirichlet values
        gsMatrix<> fpts = 
            it->function()->eval( bvp.patch(it->patch()).eval(  gsPointGrid<real_t>( rr ) ) );

        // Interpolate dirichlet boundary 
        gsBasis<>::uPtr h = bases[it->patch()]->boundaryBasis(it->side());
        gsGeometry<>::uPtr geo = h->interpolateAtAnchors( fpts );
        const gsMatrix<> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t k=0; k!= boundary.size(); ++k)
        {
            const int ii= mapper.bindex( (boundary)(k) , it->patch() );
            ddof.row(ii) = dVals.row(k);
        }
    }
}
