/**  gsIETIdGAssembler.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on:  2015-05-07
*/


#include <gsIETI/gsIETIdGAssembler.h>
#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>

#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsMultiGrid/gsMultiGridAdapter.h>
#include <gsMultiGrid/gsMassSmoother.h>
#include <gsSolver/gsTwoLevel.h>
#include <gsSolver/gsCompositePrecOp.h>

namespace gismo {

template<class T> class gsPoissonHeterogeneousAssembler;

template<class T>
gsIETIdGAssembler<T>::gsIETIdGAssembler(gsAssembler<T>& assembler)
    : Base(assembler)
{
    setOptions(defaultOptions());
}

template<class T>
gsOptionList gsIETIdGAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addSwitch("EnableMixedInterfaceCoupling","allows the use of mixed cG and cG coupling, requires setting of m_iCoupling.", false);

    return opt;

}

template<class T>
void gsIETIdGAssembler<T>::setOptions(const gsOptionList &opt)
{
    Base::setOptions(opt);
}


template<class T>
void gsIETIdGAssembler<T>::init()
{
    if(m_IETIoptions.opt.getSwitch("EnableMixedInterfaceCoupling"))
    {
        gsPoissonHeterogeneousAssembler<T>* pAss = dynamic_cast<gsPoissonHeterogeneousAssembler<T>*>(m_assembler);
        if(pAss==NULL)
            GISMO_ASSERT(m_iCoupling.size() == (size_t)m_patches.nInterfaces(), "coupling conditions on interfaces not set or have wrong size.");
        else if(m_iCoupling.size() != (size_t)m_patches.nInterfaces())
            m_iCoupling = pAss->getMixedInterfaceConditions();

    }
    else
    {
        m_iCoupling.clear();
        m_iCoupling.resize(m_patches.nInterfaces(),iFace::dg);
    }
    Base::init();
}

template<class T>
void gsIETIdGAssembler<T>::init_mappers()
{
    size_t cDim = m_basis.size();

    //Create different mappers;
    m_locDofsMapper.resize(m_patches.nPatches());
    for(size_t np = 0; np <m_patches.nPatches();np++)
        m_locDofsMapper[np].reserve(cDim);
    m_primalDofMapper.reserve(cDim);
    m_stdMapper.reserve(cDim);

    gsVector<index_t> sizes(m_patches.nPatches());
    for(size_t np = 0; np <m_patches.nPatches();np++)
        sizes[np]=m_basis.front().size(np)+m_numExtraBasisTotal[np];



    //create plain mappers
    for(size_t c =0; c<cDim;c++)
    {
        //-------------------------------------------------------------------------------
        for(size_t np = 0; np<m_patches.nPatches();np++)
        {
            //create mapper for one patch with extra dofs
            gsVector<index_t> s(1);
            s(0)=m_basis[c].size(np)+m_numExtraBasisTotal[np];
            gsDofMapper locdofMapper(s);

            //store it
            m_locDofsMapper[np].push_back(locdofMapper);
        }
        //-------------------------------------------------------------------------------
        //create the mapper for the multipatch
        gsDofMapper stdMapper(sizes);
        m_stdMapper.push_back(stdMapper);

        //-------------------------------------------------------------------------------
        //The primalDofMapper still needs to be finalized -> init_PrimalDofs
        gsDofMapper primalDofMapper(sizes);
        m_primalDofMapper.push_back(primalDofMapper);
    }

    //Eliminate Dirichlet Dofs
    if ( m_assembler->options().getInt("DirichletStrategy") == dirichlet::elimination)
    {
        for(size_t c =0; c<cDim;c++)
        {
            const gsDofMapper& tempMapper = m_assembler->system().colMapper(c); //only used temporary for idenifying eliminated dofs.

            //Create localDofMappers (where the boundary is eliminated)
            for(size_t np = 0; np<m_patches.nPatches();np++)
                //eliminate the dirichlet boundary on the patch
                for (typename gsBoundaryConditions<T>::const_iterator
                     it = m_bc_loc[np].dirichletBegin() ; it != m_bc_loc[np].dirichletEnd(); ++it )
                {
                    gsMatrix<index_t> b = m_basis[c].basis(np).boundary(it->ps.side());
                    m_locDofsMapper[np][c].markBoundary(0, b);

                }
            //eliminate the dirichlet boundary
            for (typename gsBoundaryConditions<T>::const_iterator
                 it = m_bConditions.dirichletBegin() ; it != m_bConditions.dirichletEnd(); ++it )
            {
                gsMatrix<index_t> b = m_basis[c].basis(it->ps.patch).boundary(it->ps.side());
                m_stdMapper[c].markBoundary(it->ps.patch, b);
                m_primalDofMapper[c].markBoundary(it->ps.patch, b);
                //tempMapper.markBoundary(it->ps.patch, *b);
            }

            //-------------------------------------------------------------------------------

            //tempMapper.finalize();//finalize it in order to use the check for elimination
            //-------------------------------------------------------------------------------
            //set elimination for extra dg Basis
            int ii=0;
            for(typename gsMultiPatch<T>::iiterator it  = m_patches.iBegin();it<m_patches.iEnd();it++)
            {
                if(m_iCoupling[ii++]!=iFace::dg)
                    continue;

                boundaryInterface iface = *it;
                gsMatrix<index_t> bound1 = m_boundaryBasis[iface.first().patch][iface.first().side()];
                gsMatrix<index_t> bound2 = m_boundaryBasis[iface.second().patch][iface.second().side()];

                unsigned n1 = m_basis[c].size(iface.first().patch);
                unsigned n2 = m_basis[c].size(iface.second().patch);
                for(int i =0; i<bound1.rows();i++)
                {
                    if(tempMapper.is_boundary(bound1(i,0),iface.first().patch))
                    {
                        index_t dof = n2+dgOffset(iface.second().patch,iface.second().side())+i;
                        m_locDofsMapper[iface.second().patch][c].eliminateDof(dof,0);
                        m_stdMapper[c].eliminateDof(dof,iface.second().patch);
                        m_primalDofMapper[c].eliminateDof(dof,iface.second().patch);

                        m_elimExtraBasisConnection.push_back(
                                    std::pair<patchDof,patchDof>(patchDof(iface.first().patch,bound1(i,0)),
                                                                 patchDof(iface.second().patch,dof))
                                    );
                    }
                }
                for(int i =0; i<bound2.rows();i++)
                {
                    if(tempMapper.is_boundary(bound2(i,0),iface.second().patch))
                    {
                        index_t dof = n1+dgOffset(iface.first().patch,iface.first().side())+i;
                        m_locDofsMapper[iface.first().patch][c].eliminateDof(dof,0);
                        m_stdMapper[c].eliminateDof(dof,iface.first().patch);
                        m_primalDofMapper[c].eliminateDof(dof,iface.first().patch);

                        m_elimExtraBasisConnection.push_back(
                                    std::pair<patchDof,patchDof>(patchDof(iface.second().patch,bound2(i,0)),
                                                                 patchDof(iface.first().patch,dof))
                                    );
                    }
                }
            }

        }


    }
    //"Match" the interface
    for(size_t c =0; c<cDim;c++)
    {
        int ii=0;
        for(typename gsMultiPatch<T>::iiterator it  = m_patches.iBegin();it<m_patches.iEnd();it++)
        {
            boundaryInterface iface = *it;
            gsMatrix<index_t>
                    b1= m_boundaryBasis[iface.first().patch][iface.first() .side()],
                    b2= m_boundaryBasis[iface.second().patch][iface.second() .side()];

            unsigned n1= m_basis[c].size(iface.first().patch);
            unsigned n2= m_basis[c].size(iface.second().patch);
            // match interface dofs
            if(m_iCoupling[ii]==iFace::dg)
            {
                for (int i = 0; i<b1.size(); ++i)
                    m_stdMapper[c].matchDof(iface.first().patch, (b1)(i,0), iface.second().patch,n2+dgOffset(iface.second().patch,iface.second().side())+i);
                for (int i = 0; i<b2.size(); ++i)
                    m_stdMapper[c].matchDof(iface.second().patch, (b2)(i,0), iface.first().patch,n1+dgOffset(iface.first().patch,iface.first().side())+i);
            }
            else if(m_iCoupling[ii]==iFace::conforming)
            {
                for (int i = 0; i<b1.size(); ++i)
                    m_stdMapper[c].matchDofs(iface.first().patch,b1, iface.second().patch,b2);
            }
            ii++;
        }
    }

    //Finilize them
    for(size_t c =0; c<cDim;c++)
    {
        for(size_t np = 0; np<m_patches.nPatches();np++)
            m_locDofsMapper[np][c].finalize();
        m_stdMapper[c].finalize();
    }

    //Set offsets for vector valued ones.
    unsigned offset = 0, offset2 = 0, offsetb = 0,offsetb2=0;
    for(size_t np = 0;np<m_patches.nPatches();np++)
    {
        offset = offsetb = 0;
        for(size_t c =0; c<cDim;c++)
        {
            m_locDofsMapper[np][c].setShift(offset);
            offset+=m_locDofsMapper[np][c].freeSize();

            m_locDofsMapper[np][c].setBoundaryShift(offsetb);
            offsetb+=m_locDofsMapper[np][c].boundarySize();
        }
    }
    for(size_t c =0; c<cDim;c++)
    {
        m_stdMapper[c].setShift(offset2);
        offset2+=m_stdMapper[c].freeSize();

        m_stdMapper[c].setBoundaryShift(offsetb2);
        offsetb2+=m_stdMapper[c].boundarySize();
    }

    //set mapper for edge averages
    if(m_patches.parDim() ==3 && m_IETIoptions.strat.contains(primalDofMethod::edges))
    {
        gsVector<index_t> entr(m_patches.nPatches());
        for(size_t np=0;np<m_patches.nPatches();np++)
            entr[np]=3*12;
        //Each patch in 3D consists of 12 edges
        m_edgeMapper = gsDofMapper(entr);

        // Set the static members of Edge
        Base::InstantiateEdge();
    }
}

template<class T>
void gsIETIdGAssembler<T>::init_patchInterfaceSides()
{
    m_patchISides.resize(m_patches.nPatches());

    m_numExtraBasis.resize(m_patches.nPatches());
    m_boundaryBasis.resize(m_patches.nPatches());
    m_numExtraBasisTotal.resize(m_patches.nPatches(),0);
    m_patchSideToExtraIndex.resize(m_patches.nPatches());

    for(size_t np=0; np<m_patches.nPatches();np++)
    {
        m_numExtraBasis[np].resize(2*m_patches.parDim()+1,0);
        m_boundaryBasis[np].resize(2*m_patches.parDim()+1);
    }

    int i=0;
    for(typename gsMultiPatch<T>::iiterator it  = m_patches.iBegin();it<m_patches.iEnd();it++)
    {
        patchSide side = (*it).first();
        patchSide side2 = (*it).second();

        m_patchISides[side.patch].push_back(side);
        m_patchISides[side2.patch].push_back(side2);

        m_boundaryBasis[side.patch][side.index()] = m_basis[0].basis(side.patch).boundary(side);
        m_boundaryBasis[side2.patch][side2.index()] = m_basis[0].basis(side2.patch).boundary(side2);


        if(m_iCoupling[i]==iFace::dg)
        {
            m_numExtraBasis[side.patch][side.index()]=m_boundaryBasis[side2.patch][side2.index()].rows();
            m_numExtraBasis[side2.patch][side2.index()]=m_boundaryBasis[side.patch][side.index()].rows();

            m_numExtraBasisTotal[side.patch]+=m_numExtraBasis[side.patch][side.index()];
            m_numExtraBasisTotal[side2.patch]+=m_numExtraBasis[side2.patch][side2.index()];

            m_patchSideToExtraIndex[side.patch][side2] =  m_patchISides[side.patch].size()-1;
            m_patchSideToExtraIndex[side2.patch][side] =  m_patchISides[side2.patch].size()-1;
        }
        else
        {
            m_numExtraBasis[side.patch][side.index()]=0;
            m_numExtraBasis[side2.patch][side2.index()]=0;
        }
        ++i;

    }
    i=0;
    //Do again a round over the interfaces to add also the additional interface of an patch AT THE END.
    for(typename gsMultiPatch<T>::iiterator it  = m_patches.iBegin();it<m_patches.iEnd();it++)
    {
        if(m_iCoupling[i]==iFace::dg)
        {
            m_patchISides[(*it).first().patch].push_back((*it).second());
            m_patchISides[(*it).second().patch].push_back((*it).first());


        }
        ++i;
    }

    m_patchAveragesSides= m_patchISides;

}


template<class T>
void gsIETIdGAssembler<T>::init_primalDofs()
{
    int dim = m_patches.parDim();
    size_t cDim = m_basis.size();
    m_primalVdofs.resize(m_patches.nPatches());
    m_primalEdges.resize(m_patches.nPatches());

    std::vector<boxCorner> corners, corners2;
    gsMatrix<index_t> b1,b2;
    gsVector<index_t>  bSize(dim-1);
    for(int i=0; i<dim-1;i++)
        bSize[i]=2;
    int ii=0;
    for(typename gsMultiPatch<T>::iiterator it  = m_patches.iBegin();it<m_patches.iEnd();it++)
    {
        boundaryInterface iface = (*it);
        patchSide side = iface.first();
        patchSide side2 = iface.second();

        //unsigned n2= m_basis.front().size(side2.patch);

        if(m_IETIoptions.strat.contains(primalDofMethod::vertices))
        {
            corners.clear();
            corners2.clear();
            side.getContainedCorners(dim,corners);
            side2.getContainedCorners(dim,corners2);

            b1.setZero(corners.size(),1);
            b2.setZero(corners2.size(),1);

            for(size_t c=0;c<cDim;c++)
            {
                for(size_t i=0;i<corners.size();i++)
                {
                    b1(i,0) =m_basis[c][side.patch].functionAtCorner(corners[i]);
                    b2(i,0) =m_basis[c][side2.patch].functionAtCorner(corners2[i]);
                    if(m_iCoupling[ii]==iFace::dg )
                    {
                        if(m_locDofsMapper[side.patch][c].is_free(b1(i,0) ))
                        {
                            m_primalVdofs[side.patch].push_sorted_unique(compCalc(side.patch,b1(i,0),c));
                            index_t extraV2 =dgFindCorrespondingExtraIndex(side2, side ,b1(i,0));

                            //Add the extra Vertex
                            m_primalVdofs[side2.patch].push_sorted_unique(compCalc(side2.patch,extraV2,c));

                            GISMO_ASSERT(m_stdMapper[c].index(b1(i,0),side.patch)==m_stdMapper[c].index(extraV2,side2.patch), "Corner points are not coupled correctly!");

                            //match them
                            m_primalDofMapper[c].matchDof(side.patch,b1(i,0),side2.patch,extraV2);
                        }
                        if(m_locDofsMapper[side2.patch][c].is_free(b2(i,0) ) )
                        {
                            m_primalVdofs[side2.patch].push_sorted_unique(compCalc(side2.patch,b2(i,0),c));
                            index_t extraV1 =dgFindCorrespondingExtraIndex(side, side2 ,b2(i,0));

                            //Add the extra Vertex
                            m_primalVdofs[side.patch].push_sorted_unique(compCalc(side.patch,extraV1,c));

                            GISMO_ASSERT(m_stdMapper[c].index(b2(i,0),side2.patch)==m_stdMapper[c].index(extraV1,side.patch), "Corner points are not coupled correctly!");

                            //match them
                            m_primalDofMapper[c].matchDof(side2.patch,b2(i,0),side.patch,extraV1);
                        }
                    }
                    else if(iFace::conforming)
                    {
                        if(m_locDofsMapper[side.patch][c].is_free(b1(i,0) )&& m_locDofsMapper[side2.patch][c].is_free(b2(i,0) ) )
                        {
                            m_primalVdofs[side.patch].push_sorted_unique(compCalc(side.patch,b1(i,0),c));
                            m_primalVdofs[side2.patch].push_sorted_unique(compCalc(side2.patch,b2(i,0),c));

                            m_primalDofMapper[c].matchDof(side.patch,b1(i,0),side2.patch,b2(i,0));
                        }
                    }

                }


            }

        }
        if(dim ==3 && m_IETIoptions.strat.contains(primalDofMethod::edges))
        {
            std::vector<Edge> e1, e2;
            e1.reserve(4);
            e2.reserve(4);

            Base::getEdgesFromBoundary(side,e1);
            Base::getEdgesFromBoundary(side2,e2);

            gsMatrix<index_t> num(4,1);
            num << 0,1,2,3;
            iface.reorderCorners(num);

            for(size_t i = 0; i<e1.size();i++)
            {
                GISMO_ASSERT(side.patch == e1[i].side3D.patch && side2.patch == e2[i].side3D.patch, "Edges not matching with Patch");
                if(!e1[num(i,0)].isEliminated && !e2[i].isEliminated && e1[num(i,0)].number>2 && e2[i].number>2)
                {
                    m_primalEdges[side.patch].push_sorted_unique(e1[num(i,0)]);
                    m_primalEdges[side2.patch].push_sorted_unique(e2[i]);

                    Edge e2_extra(e1[num(i,0)]);
                    if(m_iCoupling[ii]==iFace::dg)
                    {
                        e2_extra.edgeNumber = e2[i].edgeNumber +12;
                        if(side2.direction()== (int)math::max(e2[i].directions[0],e2[i].directions[1]))
                            e2_extra.edgeNumber +=12;

                        Edge e1_extra(e2[i]);
                        e1_extra.edgeNumber = e1[num(i,0)].edgeNumber +12;

                        if(side.direction()== (int)math::max(e1[num(i,0)].directions[0],e1[num(i,0)].directions[1]))
                            e1_extra.edgeNumber +=12;

                        m_primalEdges[side.patch].push_sorted_unique(e1_extra);
                        m_primalEdges[side2.patch].push_sorted_unique(e2_extra);


                        m_edgeMapper.matchDof(side.patch,e1[num(i,0)].edgeNumber,side2.patch,e2_extra.edgeNumber);
                        m_edgeMapper.matchDof(side.patch,e1_extra.edgeNumber,side2.patch,e2[i].edgeNumber);

                        // m_edgeMapper.matchDof(side.patch,e1[num(i,0)].edgeNumber,side2.patch,e2[i].edgeNumber); //changeEDGE
                        // m_edgeMapper.matchDof(side.patch,e1_extra.edgeNumber,side2.patch,e2_extra.edgeNumber);  //changeEDGE
                    }
                    else
                        m_edgeMapper.matchDof(side.patch,e1[num(i,0)].edgeNumber,side2.patch,e2[i].edgeNumber);


                }


            }

        }
        ii++;
    }

    unsigned off = 0;
    std::vector<index_t> pShifts;

    for(size_t c=0;c<cDim;c++)
    {
        pShifts.push_back(off);
        m_primalDofMapper[c].finalize();
        m_primalDofMapper[c].setShift(off);
        off+=m_primalDofMapper[c].freeSize();
    }

    m_edgeMapper.finalize();


    //populate the mappings for the Primal dofs (loc->glob)
    m_pDofsLoc2Glob.resize(m_patches.nPatches());

    if(m_IETIoptions.strat.contains(primalDofMethod::vertices))
    {
        for(size_t np=0; np<m_patches.nPatches();np++)
        {
            for(std::vector<index_t>::iterator it=m_primalVdofs[np].begin();  it!=m_primalVdofs[np].end();it++)
            {
                size_t c=getComp(np,*it);

                // the number of standard dofs (not coupled and not eliminated
                int offset =(m_primalDofMapper[c].freeSize()-m_primalDofMapper[c].coupledSize());
                int idx= m_primalDofMapper[c].index(compCalcBack(np,*it),np)-pShifts[c]-offset; //so the primal dofs have get a index starting from zero.

                //add the other components.
                for(size_t cc = 0;cc<c;cc++)
                    idx+= m_primalDofMapper[cc].coupledSize();

                m_pDofsLoc2Glob[np].push_back(idx);
            }
        }


    }

    // m_pDofsLoc2Glob[np][offset : end] also corresponds to patchISides[np][1:end]
    int offset = 0;
    for(size_t c=0;c<cDim;c++)
        offset+= m_primalDofMapper[c].coupledSize();

    if((dim == 2 && m_IETIoptions.strat.contains(primalDofMethod::edges)) || (dim == 3 && m_IETIoptions.strat.contains(primalDofMethod::faces) ) )
    {
        for(typename gsMultiPatch<T>::iiterator it  = m_patches.iBegin();it<m_patches.iEnd();it++)
        {
            patchSide side = (*it).first();
            patchSide side2 = (*it).second();

            for(size_t c = 0;c<cDim;c++)
            {
                m_pDofsLoc2Glob[side.patch].push_back(offset);
                m_pDofsLoc2Glob[side2.patch].push_back(offset);
                offset++;

            }
        }

        //Do again a round over the interfaces to add also the additional interface of an patch AT THE END. This corrects also the info structure
        //ChangeFACE
        /*
        for(typename gsMultiPatch<T>::iiterator it  = m_patches.iBegin();it<m_patches.iEnd();it++)
        {
            patchSide side = (*it).first();
            patchSide side2 = (*it).second();

            for(size_t c = 0;c<cDim;c++)
            {
                m_pDofsLoc2Glob[side.patch].push_back(offset);
                m_pDofsLoc2Glob[side2.patch].push_back(offset);
                offset++;
            }
        }
*/


    }


    if(dim ==3 && m_IETIoptions.strat.contains(primalDofMethod::edges))
    {
        for(size_t np=0; np<m_patches.nPatches();np++)
        {
            //unsigned count=0; //changeEDGE
            for(std::vector<Edge>::iterator it=m_primalEdges[np].begin();  it!=m_primalEdges[np].end();it++)
            {
                //the number of non matched edges (to start from 0)
                int offset2 =(m_edgeMapper.freeSize()-m_edgeMapper.coupledSize());
                // the index starts from the number of Vertex and Face dofs (==offset).
                int idx= m_edgeMapper.index((*it).edgeNumber,np)-offset2+offset;

                for(size_t c=0; c<cDim;c++)
                    m_pDofsLoc2Glob[np].push_back(idx+c*m_edgeMapper.coupledSize());

                // count++; //changeEDGE
                // if(count >= m_primalEdges[np].size()/2) //changeEDGE
                //    break;
            }
        }
    }

}


template<class T>
void gsIETIdGAssembler<T>::init_boundaryDofs()
{
    unsigned cDim = m_basis.size();

    m_boundDofs.resize(m_patches.nPatches());
    for ( typename gsMultiPatch<T>::iiterator it =
          m_patches.iBegin(); it != m_patches.iEnd(); ++it )
    {
        boundaryInterface iFace=(*it);

        patchSide p1 = iFace.first();
        patchSide p2 = iFace.second();


        for(size_t c=0; c<cDim;c++)
        {
            gsMatrix<index_t> b1 = m_boundaryBasis[p1.patch][p1.index()];
            gsMatrix<index_t> b2 = m_boundaryBasis[p2.patch][p2.index()];

            for(int i=0;i<b1.rows();i++)
            {
                if(m_locDofsMapper[p1.patch][c].is_free((b1)(i,0)))
                    m_boundDofs[p1.patch].push_sorted_unique(compCalc(p1.patch,(b1)(i,0),c));
            }
            for(int i=0;i<b2.rows();i++)
            {
                if(m_locDofsMapper[p2.patch][c].is_free((b2)(i,0)))
                    m_boundDofs[p2.patch].push_sorted_unique(compCalc(p2.patch,(b2)(i,0),c));
            }

            //The extra basis functions (the extra dofs might be eliminated)
            size_t n1 = m_basis[c].size(p1.patch);
            size_t n2 = m_basis[c].size(p2.patch);

            for(int i=0; i<m_numExtraBasis[p1.patch][p1.index()];i++)
                if(m_locDofsMapper[p1.patch][c].is_free(n1 +dgOffset(p1.patch,p1.side())+i))
                    m_boundDofs[p1.patch].push_sorted_unique(compCalc(p1.patch, n1 +dgOffset(p1.patch,p1.side())+i ,c));

            for(int i=0; i<m_numExtraBasis[p2.patch][p2.index()];i++)
                if(m_locDofsMapper[p2.patch][c].is_free(n2 +dgOffset(p2.patch,p2.side())+i))
                    m_boundDofs[p2.patch].push_sorted_unique(compCalc(p2.patch, n2 +dgOffset(p2.patch,p2.side())+i ,c));
        }
    }

    std::vector<index_t>::iterator it;
    unsigned iI,iR,sum,idx;;

    m_glob2BoundInteriorIndex.resize(m_patches.nPatches());
    m_globIsBoundIndex.resize(m_patches.nPatches());

    for(size_t np=0; np<m_patches.nPatches();np++)
    {
        iI = iR = 0;

        //and fill up m_glob2BoundIndex[np]
        sum=0;
        for(size_t c=0;c<cDim;c++)
            sum+=m_locDofsMapper[np][c].freeSize();

        m_globIsBoundIndex[np].resize(sum);
        m_glob2BoundInteriorIndex[np].resize(sum);

        for(size_t c=0; c<cDim;c++)
            for(int i = 0; i<m_locDofsMapper[np][c].size();i++)
            {
                idx = m_locDofsMapper[np][c].index(i);
                if(m_locDofsMapper[np][c].is_free_index(idx))
                {
                    it = std::find(m_boundDofs[np].begin(), m_boundDofs[np].end(), compCalc(np,i,c));

                    if(it!=m_boundDofs[np].end())
                    {
                        m_globIsBoundIndex[np][idx]=true;
                        m_glob2BoundInteriorIndex[np].at(idx)=std::distance(m_boundDofs[np].begin(), it);
                    }
                    else
                    {
                        m_globIsBoundIndex[np][idx]=false;
                        m_glob2BoundInteriorIndex[np].at(idx)=iI;
                        iI++;
                    }
                }
            }
    }
}

template<class T>
void gsIETIdGAssembler<T>::init_jumpOperatorData()
{
    size_t cDim = m_basis.size();
    std::set<index_t> freeElimDof;
    gsMatrix<index_t> glob1, glob2;
    for ( typename gsMultiPatch<T>::const_iiterator it =
          m_patches.iBegin(); it != m_patches.iEnd(); ++it )
    {
        boundaryInterface iFace=(*it);

        patchSide p1 = iFace.first();
        patchSide p2 = iFace.second();

        for(size_t c=0;c<cDim;c++)
        {
            gsMatrix<index_t> basis1 = m_boundaryBasis[p1.patch][p1.index()];
            gsMatrix<index_t> basis2 = m_boundaryBasis[p2.patch][p2.index()];

            gsMatrix<index_t> b1(basis1.rows()+m_numExtraBasis[p1.patch][p1.index()],1);
            gsMatrix<index_t> b2(basis2.rows()+m_numExtraBasis[p2.patch][p2.index()],1);

            b1.topRows(basis1.rows())=basis1;
            b2.bottomRows(basis2.rows())=basis2;

            for(int i=0;i<m_numExtraBasis[p1.patch][p1.index()];i++)//for(int i=0;i<basis2.rows();i++)
                b1(basis1.rows()+i,0)= m_basis[c].size(p1.patch)+dgOffset(p1.patch,p1.side())+i;
            for(int i=0;i<m_numExtraBasis[p2.patch][p2.index()];i++)//for(int i=0;i<basis1.rows();i++)
                b2(i,0)= m_basis[c].size(p2.patch)+dgOffset(p2.patch,p2.side())+i;

            m_stdMapper[c].localToGlobal(b1,p1.patch,glob1);
            m_stdMapper[c].localToGlobal(b2,p2.patch,glob2);

            GISMO_ASSERT(glob1 == glob2 , "The mapping is not matching, something went wrong");
            for(int i=0;i<glob1.rows();i++)
            {
                bool isfree1 = m_locDofsMapper[p1.patch][c].is_free((b1)(i,0));
                bool isfree2 = m_locDofsMapper[p2.patch][c].is_free((b2)(i,0));

                //if both are free just add them
                if(isfree1 && isfree2)
                {
                    m_globalConnectionTable[glob1(i,0)].insert(patchDof(p1.patch,compCalc(p1.patch,(b1)(i,0),c)) );
                    m_globalConnectionTable[glob1(i,0)].insert(patchDof(p2.patch,compCalc(p2.patch,(b2)(i,0),c)) );
                }
                else if(isfree1 || isfree2)
                {
                    //if one of the two are not free check if this global dof already has a dirichlet dof
                    if(freeElimDof.find(glob1(i,0))==freeElimDof.end())
                    {
                        //if not then add  it to the free-eliminated dofs and add both dof to the map
                        freeElimDof.insert(glob1(i,0));
                        m_globalConnectionTable[glob1(i,0)].insert(patchDof(p1.patch,compCalc(p1.patch,(b1)(i,0),c)) );
                        m_globalConnectionTable[glob1(i,0)].insert(patchDof(p2.patch,compCalc(p2.patch,(b2)(i,0),c)) );
                    }
                    else
                    {
                        //if a dirichlet dof was already there, then only add the free dof.
                        if(isfree1)
                            m_globalConnectionTable[glob1(i,0)].insert(patchDof(p1.patch,compCalc(p1.patch,(b1)(i,0),c)) );
                        if(isfree2)
                            m_globalConnectionTable[glob1(i,0)].insert(patchDof(p2.patch,compCalc(p2.patch,(b2)(i,0),c)) );
                    }
                }


            }

        }
    }

    m_lagrangeTable.reserve(m_globalConnectionTable.size());

    bool flag=false;

    bool isfreeA, isfreeB;
    for(typename std::map<index_t,std::set<patchDof> >::iterator it=m_globalConnectionTable.begin(); it!=m_globalConnectionTable.end(); ++it)
    {
        flag = false;
        if(freeElimDof.find((*it).first)!=freeElimDof.end())
            flag=true;
        for(typename std::set<patchDof>::iterator itA=(*it).second.begin(); itA!=(*it).second.end(); ++itA)
        {
            if(itA==(*it).second.end())
                break;

            //This can be improved
            if(flag)
                isfreeA = m_locDofsMapper[(*itA).first][0].is_free((*itA).second);
            else
                isfreeA = true;

            for(typename std::set<patchDof>::iterator itB=itA; itB!=(*it).second.end(); ++itB)
            {
                if(itB==itA)
                    continue;

                //This can be improved
                if(flag)
                    isfreeB = m_locDofsMapper[(*itB).first][0].is_free((*itB).second);
                else
                    isfreeB = true;

                if(isfreeA && isfreeB)
                    m_lagrangeTable.push_back(std::pair<patchDof,patchDof>(*itA,*itB));
                else if(isfreeA || isfreeB)
                {
                    m_lagrangeTable.push_back(std::pair<patchDof,patchDof>(*itA,*itB));
                    m_freeElimLagr.push_back(m_lagrangeTable.size()-1);

                }
            }


        }
    }

}

template<class T>
void gsIETIdGAssembler<T>::setInfo()
{
    info.numberPatches = m_patches.nPatches();
    info.dim = m_patches.parDim();
    info.cDim = m_basis.size();

    info.dofTotal=0;
    for(size_t c=0;c<info.cDim;c++)
        info.dofTotal += m_primalDofMapper[c].freeSize();
    info.dofTotalB=0;
    info.dofTotalI=0;

    info.dofsB.resize(m_patches.nPatches(),0);
    info.dofsI.resize(m_patches.nPatches(),0);
    info.dofsP.resize(m_patches.nPatches(),0);
    info.dofs.resize(m_patches.nPatches(),0);


    int diffMethods=3;
    info.dofTotalPtype.resize(diffMethods,0);
    info.dofsPtype.resize(m_patches.nPatches());

    //set size on patch np for dofs, boundaryDofs, interiorDofs
    for (size_t np=0; np < m_patches.nPatches(); ++np )
    {
        info.dofsPtype[np].resize(diffMethods);

        for(size_t c=0;c<info.cDim;c++)
            info.dofs[np]+=(m_locDofsMapper[np][c].freeSize());

        info.dofsB[np]=(m_boundDofs[np].size());
        info.dofsI[np]=(info.dofs[np]-info.dofsB[np]);
        info.dofTotalB+= info.dofsB[np];
        info.dofTotalI+= info.dofsI[np];
    }


    //set size for primal dofs
    info.dofTotalP=0;
    int method;
    if(m_IETIoptions.strat.contains( primalDofMethod::vertices))
    {
        method = primalDofMethod::getMethodIndex(primalDofMethod::vertices);

        for (size_t np=0; np < m_patches.nPatches(); ++np )
            info.dofsPtype[np][method] = m_primalVdofs[np].size();
        for(size_t c=0;c<info.cDim;c++)
            info.dofTotalPtype[method]+=m_primalDofMapper[c].coupledSize();

        info.dofTotalP +=  info.dofTotalPtype[method];
    }
    if((info.dim == 2 && m_IETIoptions.strat.contains( primalDofMethod::edges ) ) ||  (info.dim ==3 && m_IETIoptions.strat.contains(primalDofMethod::faces) ))
    {
        if(info.dim == 2)
            method =primalDofMethod::getMethodIndex(primalDofMethod::edges);
        else
            method =primalDofMethod::getMethodIndex(primalDofMethod::faces);

        for (size_t np=0; np < m_patches.nPatches(); ++np )
        {
            unsigned i=0;
            for( i=0;i<m_patchISides[np].size();++i)
                if((unsigned)m_patchISides[np][i].patch !=np)
                    break;
            info.dofsPtype[np][method] = info.cDim*i; // m_patchISides[np].size()/2 changeFACE <<<< --- here is a change!
        //  info.dofsPtype[np][method] = info.cDim*m_patchISides[np].size();
        }

        info.dofTotalPtype[method]= info.cDim*m_patches.nInterfaces(); //changeFACE
        //  info.dofTotalPtype[method] = info.cDim*m_patches.nInterfaces()*2;
        info.dofTotalP +=   info.dofTotalPtype[method];
    }
    if(info.dim ==3 && ( m_IETIoptions.strat.contains(primalDofMethod::edges)))
    {
        method = primalDofMethod::getMethodIndex(primalDofMethod::edges);
        for (size_t np=0; np < m_patches.nPatches(); ++np )
            // info.dofsPtype[np][method] = info.cDim*m_primalEdges[np].size()/2; //changeEDGE
            info.dofsPtype[np][method] = info.cDim*m_primalEdges[np].size();
        info.dofTotalPtype[method]= info.cDim*m_edgeMapper.coupledSize();
        info.dofTotalP +=   info.dofTotalPtype[method];
    }


    //set size for the remaining dofs and the sum of primal dofs on patch np
    info.dofsR.resize(m_patches.nPatches());
    for (size_t np=0; np < m_patches.nPatches(); ++np )
    {
        unsigned sum =0;
        for(int i = 0; i<diffMethods;i++)
            sum+= info.dofsPtype[np][i];
        info.dofsP[np]=sum;

        //remaining are all dofs minus primal vertex dofs (only for noMinimalenergy)
        info.dofsR[np]=info.dofs[np] - info.dofsPtype[np][primalDofMethod::getMethodIndex(primalDofMethod::vertices)];
    }

    //for minimal energy, we only consider boundary dofs
    if(!m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
        info.dofs = info.dofsB;

    //number of lagrange multiplier
    info.lagrangeMult = m_lagrangeTable.size();

    //number of global Interface dofs
    info.dofInterface= m_globalConnectionTable.size();
}

template<class T>
void gsIETIdGAssembler<T>::init_remainingDofs()
{
    size_t cDim = m_basis.size();
    std::vector<index_t>::iterator it;

    unsigned sum,idx,iR;
    m_glob2PrimalRemainingIndex.resize(m_patches.nPatches());
    m_globIsPrimalIndex.resize(m_patches.nPatches());

    for(size_t np=0; np<m_patches.nPatches();np++)
    {
        iR=sum=0;
        for(size_t c=0;c<cDim;c++)
            sum+=m_locDofsMapper[np][c].freeSize();

        m_globIsPrimalIndex[np].resize(sum);
        m_glob2PrimalRemainingIndex[np].resize(sum);

        for(size_t c=0; c<cDim;c++)
            for(int i = 0; i<m_locDofsMapper[np][c].size();i++)
            {
                idx = m_locDofsMapper[np][c].index(i);
                if(m_locDofsMapper[np][c].is_free_index(idx))
                {
                    it = std::find(m_primalVdofs[np].begin(), m_primalVdofs[np].end(), compCalc(np,i,c));
                    if(it!=m_primalVdofs[np].end())
                    {
                        m_globIsPrimalIndex[np][idx]=true;
                        m_glob2PrimalRemainingIndex[np][idx]=std::distance(m_primalVdofs[np].begin(), it);
                    }
                    else
                    {
                        m_globIsPrimalIndex[np][idx]=false;
                        m_glob2PrimalRemainingIndex[np][idx]=iR;
                        iR++;
                    }
                }
            }
    }

    m_remDofs.resize(info.numberPatches);
    m_remDofIsAlsoBoundDof.resize(info.numberPatches);

    for(size_t np=0; np<info.numberPatches;np++)
    {
        m_remDofs[np].reserve(info.dofsR[np]);
        m_remDofIsAlsoBoundDof[np].reserve(info.dofsR[np]);

        //loop over all dofs, included the dirichlet ones
        for(size_t c = 0; c< info.cDim;c++)
        {
            for(index_t d = 0; d<info.dofs[np]+m_locDofsMapper[np][c].boundarySize();d++)
            {
                idx = m_primalDofMapper[c].index(d,np);
                if(!m_primalDofMapper[c].is_coupled_index(idx) && !m_primalDofMapper[c].is_boundary_index(idx))
                {
                    m_remDofs[np].push_back(compCalc(np,d,c));

                    //also test if d is a interface dof, i.e. is in m_boundDofs[np]:
                    if(std::find(m_boundDofs[np].begin(), m_boundDofs[np].end(), compCalc(np,d,c))!=m_boundDofs[np].end())
                        m_remDofIsAlsoBoundDof[np].push_back(true);
                    else
                        m_remDofIsAlsoBoundDof[np].push_back(false);
                }
            }
        }
    }
}

template<class T>
void gsIETIdGAssembler<T>::assemble( const gsMultiPatch<T> &curSol)
{

    Base::assembleInit();
    bool assembleRHS = !m_IETIoptions.opt.getSwitch("NoAssemblyOfRHS");
    std::vector<gsSparseMatrix<T> > matrices(info.numberPatches);
    std::vector<gsSparseMatrix<T> > matrices2(info.numberPatches);
    std::vector<gsMatrix<T> > rhs_loc(info.numberPatches);
    std::vector<gsMatrix<T> > spp_loc(info.numberPatches);

    Eigen::setNbThreads(1);

    Base::reserveSpace(m_dirDofs,rhs_loc);

    gsStopwatch time, time2;
    if(!m_IETIoptions.opt.getSwitch("ExtraTiming")) //The Fast One
    {
#pragma omp parallel
        {

            gsAssembler<T>* A = m_assembler->create();

#pragma omp for schedule(static, 1) nowait
            for (size_t np=0; np < info.numberPatches; ++np )
            {

                Base::assembleLocal(A, matrices[np], rhs_loc[np],m_dirDofs[np],np,curSol);

                //gsInfo<<"K:\n\n"<<matrix.toDense()<<std::endl;

                //increase the size so that the extra dofs fit
                gsSparseMatrix<T>& matrix = matrices[np];
                matrix.uncompress();
                matrix.conservativeResize(m_locDofsMapper[np][0].freeSize(),m_locDofsMapper[np][0].freeSize());

                //reserve for more dofs
                int nonZerosPerCol = m_assembler->numColNz();
                matrix.reserve(gsVector<index_t>::Constant(matrix.cols(), nonZerosPerCol));
                //the old entries must not be deleted, so matrix.reservePerColumn(..) does not work here
            }

            delete A;

#pragma omp barrier
            Base::printTime(time,"Time for assemling all patch local matrices: ");
#pragma omp single
            {
                updateDDofs();

                time.restart();
                //assemble DG contributions
                assembleDgInterfaceContribution(matrices,rhs_loc);
                // gsInfo<<"K:\n\n"<<matrices[0].toDense()<<std::endl;
                // gsInfo<<"K:\n\n"<<matrices[1].toDense()<<std::endl;
                time.stop();
                gsInfo<<"Assemble dG Interface: "<<time<<std::endl;

                for (size_t np=0; np < info.numberPatches; ++np )
                    matrices[np].makeCompressed();
            }

#pragma omp for schedule(static, 1) nowait
            for (size_t np=0; np < info.numberPatches; ++np )
            {
                // -- Reordering
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                {
                    matrices2[np] = matrices[np];
                    Base::makeReorderingPrimalRem(matrices2[np],np);
                }
                Base::makeReorderingBoundInt(matrices[np],np);

                //For Preconditioning and processing the rhs and solution ICO minimum energy
                Base::assembleKiiKibKbb(matrices[np], np);

                // The main matrices for applying the system matrix
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                    Base::assembleKrrKrpKpp(matrices2[np],np);
                else
                {
                    assembleC(np);
                    Base::assembleLUofKC(matrices[np],np);
                }


                //assemble local Spp
                Base::assembleSppLoc(spp_loc[np], np);

                //assemble local rhs
                if(assembleRHS) Base::assembleRhs(rhs_loc, np);

            } //end-for
        }//end-parallel
    }//end if
    else // The slow one
    {

#pragma omp parallel
        {
            gsAssembler<T>* A = m_assembler->create();
#pragma omp for schedule(static, 1) nowait
            for (size_t np=0; np < info.numberPatches; ++np )
            {
                Base::assembleLocal(A, matrices[np], rhs_loc[np],m_dirDofs[np], np,curSol);

                gsSparseMatrix<T>& matrix = matrices[np];
                matrix.uncompress();
                //increase the size so that the extra dofs fit
                matrix.conservativeResize(m_locDofsMapper[np][0].freeSize(),m_locDofsMapper[np][0].freeSize());

                //reserve for more dofs
                int nonZerosPerCol = m_assembler->numColNz();
                matrix.reserve(gsVector<index_t>::Constant(matrix.cols(), nonZerosPerCol));
                //the old entries must not be deleted, so matrix.reservePerColumn(..) does not work here
            }
            delete A;
        }

        Base::printTime(time,"Time for assemling all patch local matrices: ");

        //fill up the remaining dirichlet dofs (only for the free-eliminated ones)
        updateDDofs();

        time.restart();
        //assemble DG contributions
        assembleDgInterfaceContribution(matrices,rhs_loc);
        Base::printTime(time,"Assemble dG Interface: ");

        /*
        int max=0;
        for(unsigned np=0; np<info.numberPatches;np++)
        {
          //  gsInfo<<"the matrix starts here: \n"<<matrices[np].toDense()<<"\n\n";
            max=0;
            for(int i=0; i<matrices[np].cols();++i)
            {
                if(max<matrices[np].col(i).nonZeros())
                    max = matrices[np].col(i).nonZeros();

            }
            gsInfo<<"Max number of entries: "<<max <<"\n"<<"usage: "<<(real_t)max/m_assembler->options().numColNz(m_basis[0][0])<<"\n";
        }*/

        time.restart();

#pragma omp parallel
        {
            // -- Reordering
#pragma omp for schedule(static, 1)
            for (size_t np=0; np < info.numberPatches; ++np )
            {
                matrices[np].makeCompressed();
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                {
                    matrices2[np] = matrices[np];
                    Base::makeReorderingPrimalRem(matrices2[np],np);
                }
                Base::makeReorderingBoundInt(matrices[np],np);
            }

            Base::printTime(time,"Time for reordering all patch local matrices: ");

            //------------------------------------------------------------------------------/
            //For Preconditioning and processing the rhs and solution ICO minimum energy


#pragma omp for schedule(static, 1)
            for (size_t np=0; np < info.numberPatches; ++np )
                Base::assembleKiiKibKbb(matrices[np], np);

            Base::printTime(time,"Time for calculating LU of Kii : ");


            // The main matrices for applying the system matrix
            if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
            {
#pragma omp for schedule(static, 1)
                for (size_t np=0; np < info.numberPatches; ++np )
                    Base::assembleKrrKrpKpp(matrices2[np],np);

                Base::printTime(time,"Time for computing LU factorization of Krr : ");
            }
            else
            {
#pragma omp for schedule(static, 1)
                for (size_t np=0; np < info.numberPatches; ++np )
                    assembleC(np);

                Base::printTime(time,"Time for computing C : ");

#pragma omp for schedule(static, 1)
                for (size_t np=0; np < info.numberPatches; ++np )
                    Base::assembleLUofKC(matrices[np],np);

                Base::printTime(time,"Time for computing LU factorization of KC : ");
            }


#pragma omp for schedule(static, 1)
            for (size_t np=0; np < info.numberPatches; ++np )
                Base::assembleSppLoc(spp_loc[np], np);

            Base::printTime(time,"Time for assemling local Spp: ");

#pragma omp for schedule(static, 1)
            for (size_t np=0; np < info.numberPatches; ++np )
                if(assembleRHS) Base::assembleRhs(rhs_loc, np);

            Base::printTime(time,"Time for assemling rhs: ");
        }
    }
    Base::printTime(time2,"Time for doing parallel assembling part: ");
    //------------Do Serial stuff--------------
    Eigen::setNbThreads(0);
    //1. right hand side
    if(assembleRHS) Base::assembleRhsFreeElim();

    //2. Spp
    Base::assembleSpp(spp_loc);

    Base::printTime(time2,"Time for doing serial assembling part: ");
}

template<class T>
void gsIETIdGAssembler<T>::updateDDofs()
{
    //fill up the remaining dirichlet dofs (only for the free-eliminated ones). Needed for assembleDgInterfaceContribution
    for(typename std::vector<std::pair<patchDof,patchDof> >::iterator  it=m_elimExtraBasisConnection.begin(); it!=m_elimExtraBasisConnection.end(); ++it)
    {
        patchDof p1= (*it).first; //the original one
        patchDof p2 = (*it).second; //the extra Basis one
        int c = getComp(p2.first,p2.second); // p1 and p2 must have the same component

        unsigned bIdxExtra = m_locDofsMapper[p2.first][c].bindex(compCalcBack(p2.first,p2.second));
        unsigned bIdx = m_locDofsMapper[p1.first][c].bindex(compCalcBack(p1.first,p1.second));
        m_dirDofs[p2.first].row(bIdxExtra)= m_dirDofs[p1.first].row(bIdx);
    }
}

template<class T>
void gsIETIdGAssembler<T>::assembleDgInterfaceContribution(std::vector<gsSparseMatrix<T> >& matrices, std::vector<gsMatrix<T> > & rhs) const
{
    gsMatrix<index_t> *actives1,*actives2;
    gsMatrix<index_t> activesExtra1, activesExtra2;

    gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
    gsVector<T> quWeights;         // Mapped weights

    unsigned evFlags(0);     // Evaluation flags for the Geometry map
    gsQuadRule<T> QuRule;

    //TODO: make this stuff more generic.
    gsPoissonHeterogeneousAssembler<T>* pAss = dynamic_cast< gsPoissonHeterogeneousAssembler<T>*>(m_assembler);
    GISMO_ASSERT(pAss !=NULL,
                 "In order to use dG-IETI, the assembler needs to be derived from gsPoissonHeterogeneousAssembler<T>");

    int ii=0;
    for ( typename gsMultiPatch<T>::const_iiterator it = m_patches.iBegin(); it != m_patches.iEnd(); ++it )
    {
        if(m_iCoupling[ii++]!=iFace::dg)
            continue;

        const boundaryInterface & bi =
                ( m_basis[0][it->first() .patch].numElements(it->first() .side() ) <
                m_basis[0][it->second().patch].numElements(it->second().side() ) ?
                    it->getInverse() : *it );

        gsVisitorDg2<T,1> dg =  pAss->visitorDg(bi); //makes a copy, but class is very small.

        const gsAffineFunction<T> interfaceMap(m_patches.getMapForInterface(bi));

        const int patch1      = bi.first().patch;
        const int patch2      = bi.second().patch;
        const gsBasis<T> & B1 = m_basis[0][patch1];// (!) unknown 0
        const gsBasis<T> & B2 = m_basis[0][patch2];

        const int bSize1      = B1.numElements( bi.first() .side() );
        const int bSize2      = B2.numElements( bi.second().side() );

        const int ratio = bSize1 / bSize2;
        GISMO_ASSERT(bSize1 >= bSize2 && bSize1%bSize2==0,
                     "DG assumes nested interfaces. Got bSize1="<<
                     bSize1<<", bSize2="<<bSize2<<"." );

        // Initialize
        dg.initialize(B1, B2, QuRule, evFlags);

        // Initialize geometry evaluators
        typename gsGeometryEvaluator<T>::uPtr geoEval1(getEvaluator(evFlags, m_patches[patch1]));
        typename gsGeometryEvaluator<T>::uPtr geoEval2(getEvaluator(evFlags, m_patches[patch2]));

        // Initialize domain element iterators
        typename gsBasis<T>::domainIter domIt1 = B1.makeDomainIterator( bi.first() .side() );
        typename gsBasis<T>::domainIter domIt2 = B2.makeDomainIterator( bi.second().side() );

        int count = 0;
        // iterate over all boundary grid cells on the "left"
        for (; domIt1->good(); domIt1->next() )
        {

            count++;
            // if(patch1 == 0 && patch2 == 2 && count ==2)
            //    gsInfo<<"Stop\n";
            // Compute the quadrature rule on both sides
            QuRule.mapTo( domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights);
            interfaceMap.eval_into(quNodes1,quNodes2);

            // Perform required evaluations on the quadrature nodes
            dg.evaluate(B1, *geoEval1, B2, *geoEval2, quNodes1, quNodes2);

            // Assemble on element
            dg.assemble(*domIt1,*domIt2, *geoEval1, *geoEval2, quWeights);

            //extract the actives and prepare them for the IETI_locToGlob map
            dg.getActives(actives1, actives2);
            prepareActives(bi,*actives1,*actives2,activesExtra1,activesExtra2 );

            //do the map
            dg.localToGlobalIETI(m_locDofsMapper[patch1][0],m_dirDofs[patch1],
                    activesExtra1, matrices[patch1], rhs[patch1]);

            dg.revert();

            dg.localToGlobalIETI(m_locDofsMapper[patch2][0],m_dirDofs[patch2],
                    activesExtra2, matrices[patch2], rhs[patch2]);

            if ( count % ratio == 0 ) // next master element ?
            {
                domIt2->next();
            }

        }
    }

}

template<class T>
void gsIETIdGAssembler<T>::prepareActives(const boundaryInterface & bi, gsMatrix<index_t>& actives1,gsMatrix<index_t>& actives2,
                                          gsMatrix<index_t>& activesExtra1,gsMatrix<index_t>& activesExtra2) const
{
    const int patch1      = bi.first().patch;
    const int patch2      = bi.second().patch;
    const int n1          = actives1.rows();
    const int n2          = actives2.rows();

    actives1.conservativeResize(n1+n2, actives1.cols());
    actives2.conservativeResize(n2+n1, actives2.cols());
    actives1.bottomRows(n2).setZero();
    actives2.bottomRows(n1).setZero();

    activesExtra1.resize(m_numExtraBasis[patch1][bi.first().index()],1);
    activesExtra2.resize(m_numExtraBasis[patch2][bi.second().index()],1);
    int iter1=0;
    for(int i=0; i<n2;i++)
    {
        int k=dgFindCorrespondingExtraIndex(bi.first(), bi.second() ,(actives2)(i,0));
        if(k>=0)
        {
            (actives1)(n1+i,0)= k;
            activesExtra1(iter1,0)=i;
            iter1++;
        }
    }

    int iter2=0;
    for(int i=0; i<n1;i++)
    {
        int k=dgFindCorrespondingExtraIndex(bi.second(), bi.first() ,(actives1)(i,0));
        if(k>=0)
        {
            (actives2)(n2+i,0)=k;
            activesExtra2(iter2,0)=i;
            iter2++;
        }
    }
    //shrink to the appropriate size, since not all extra dofs are active on an element!
    activesExtra1.conservativeResize(iter1, activesExtra1.cols());
    activesExtra2.conservativeResize(iter2, activesExtra2.cols());
}

template<class T>
void gsIETIdGAssembler<T>::assembleC(size_t np)
{
    gsMatrix<T> val;
    std::vector<Trip> tripletList;


    m_C[np].resize(info.dofsP[np],info.dofsB[np]);
    int estimated_num=0;

    if(m_IETIoptions.strat.contains(primalDofMethod::vertices))
        estimated_num+=info.dofsP[np];

    if((info.dim == 2 && m_IETIoptions.strat.contains(primalDofMethod::edges)) || (info.dim == 3 && m_IETIoptions.strat.contains(primalDofMethod::faces)))
        estimated_num+= info.cDim*m_patchISides[np].size() *
                cast<T,int>(floor(math::sqrt(T(info.dofTotal))));

    if(info.dim ==3 && m_IETIoptions.strat.contains(primalDofMethod::edges) && info.dofsPtype[np][primalDofMethod::getMethodIndex(primalDofMethod::edges)]!=0 )
        estimated_num+= info.dofsPtype[np][primalDofMethod::getMethodIndex(primalDofMethod::edges)]*m_primalEdges[np].front().number;

    tripletList.reserve(estimated_num);

    int offset=0;
    if(m_IETIoptions.strat.contains(primalDofMethod::vertices))
    {
        for(size_t i=0; i<m_primalVdofs[np].size(); ++i)
        {
            int c = getComp(np,m_primalVdofs[np][i]);
            int idx = m_glob2BoundInteriorIndex[np][m_locDofsMapper[np][c].index(compCalcBack(np,m_primalVdofs[np][i]))];
            tripletList.push_back(Trip(i,idx,1));
        }
        offset+=m_primalVdofs[np].size();
    }

    if((info.dim == 2 && m_IETIoptions.strat.contains(primalDofMethod::edges)) || (info.dim == 3 && m_IETIoptions.strat.contains(primalDofMethod::faces)))
    {
        int method;
        info.dim == 2 ? method = primalDofMethod::getMethodIndex(primalDofMethod::edges) : method = primalDofMethod::getMethodIndex(primalDofMethod::faces);
        for(size_t i = 0; i< m_patchISides[np].size();i++)
        {
            const patchSide& side = m_patchISides[np][i];

            for(size_t c=0;c<info.cDim;c++)
            {
                val.setZero(info.dofsB[np],1);
                Base::getInterfaceAverageValue(side, val,c);

                if(i<info.dofsPtype[np][method]/info.cDim) //m_patchISides[np].size()/2)
                {
                    for(int j = 0; j<val.rows();j++)
                        if(val(j,0)!=0 && m_stdMapper[c].is_free(m_boundDofs[np][j],np))
                            tripletList.push_back(Trip(offset+info.cDim*i+c,j,val(j,0)));
                }
                else
                {
                    int sideIdx =  getIndexOfExtraSide(np,side);//i%(m_patchISides[np].size()/2);
                    GISMO_ASSERT(sideIdx>=0,"the side you are looking for is either not a dG side or not on this patch.");
                    for(int j = 0; j<val.rows();j++)
                        if(val(j,0)!=0)
                        {
                            int idx = dgFindCorrespondingExtraIndex(m_patchISides[np][sideIdx],side,m_boundDofs[side.patch][j]);
                            if(m_stdMapper[c].is_free(idx,np))
                            {
                                int bidx = m_glob2BoundInteriorIndex[np][m_locDofsMapper[np][c].index(idx)];
                                //   tripletList.push_back(Trip(offset+info.cDim*i+c,bidx,val(j,0)));
                                tripletList.push_back(Trip(offset+info.cDim*sideIdx+c,bidx,val(j,0)));  //ChangeFACE
                            }
                        }
                }
            }
        }
        offset+=info.dofsPtype[np][method];
    }
    if(info.dim ==3 && m_IETIoptions.strat.contains(primalDofMethod::edges))
    {
        int method = m_IETIoptions.strat.contains(primalDofMethod::edges);

        int dir;
        unsigned entr;
        patchSide side;
        for(size_t e=0;e<m_primalEdges[np].size();e++)
        {
            const Edge& edge = m_primalEdges[np][e];
            Edge* edge_real=NULL;

            int patch = edge.side3D.patch;

            if(edge.edgeNumber >= 12)
            {
                for(size_t k=0; k<m_primalEdges[np].size();k++)
                {
                    GISMO_ASSERT(m_primalEdges[np][k].edgeNumber<12, "Corresponding edge not found, something went wrong.");
                    if(edge.edgeNumber%12 == m_primalEdges[np][k].edgeNumber)
                    {
                        edge_real = &m_primalEdges[np][k];
                        break;
                    }
                }

                if(edge.edgeNumber < 24)
                    dir = math::min(edge_real->directions[0],edge_real->directions[1]);
                else
                    dir = math::max(edge_real->directions[0],edge_real->directions[1]);

                if(dir == edge_real->side3D.direction())
                    side = edge_real->side3D;
                else
                    side = edge_real->side3D_2;
            }


            for(size_t c=0;c<info.cDim;c++)
            {
                int n = m_basis[c].basis(np).size();
                typename gsBasis<T>::uPtr bbasis = m_basis[c].basis(patch).boundaryBasis(edge.side3D);
                gsMatrix<index_t> bI = m_basis[c].basis(patch).boundary(edge.side3D);
                gsMatrix<index_t> bbI = bbasis->boundary(edge.side2D);

                val.setZero(bbI.rows(),1);

                Base::getFaceEdgeAverageValue(edge,val,c,false);
                for(int j = 0;j<bbI.rows();j++)
                {
                    if(edge.edgeNumber < 12)
                        entr= (bI)((bbI)(j,0),0);
                    else
                        entr = n+dgOffset(side.patch,side)+(bbI)(j,0);

                    if(m_locDofsMapper[np][c].is_free(entr))
                    {
                        int idx = m_glob2BoundInteriorIndex[np][m_locDofsMapper[np][c].index(entr)];
                        tripletList.push_back(Trip(offset+info.cDim*e+c,idx, val(j,0)));
                        // tripletList.push_back(Trip(offset+info.cDim*(e%(m_primalEdges[np].size()/2))+c,idx, val(j,0))); //changeEDGE
                    }
                }
            }
        }

        offset+=info.dofsPtype[np][method];
    }
    m_C[np].setFromTriplets(tripletList.begin(), tripletList.end());
    //   if(np==0)
    //       gsInfo<<m_C[np].toDense()<<"\n\n";
    GISMO_ASSERT(m_C[np].rows()<=m_C[np].cols(), "More primal varibles than interface dofs, cannot be linear independet. Refine your domain, or choose less primal dofs (different strategy)");
    // gsInfo<<"Patch: "<<np<<std::endl;
    // gsInfo<<m_C[np].toDense()<<std::endl<<std::endl;

}

template<class T>
void gsIETIdGAssembler<T>::setupMGforKC(size_t np)
{
    typename gsPreconditionerOp<T>::Ptr mg_ptr;

    std::vector< gsMultiBasis<> > MG_bases;
    gsBasis<T>& coarseBasis = (gsBasis<T>&)*m_patches[np].basis().clone().release();
    gsBoundaryConditions<T> bc = m_bc_loc[np];

    m_exC[np].resize(m_C[np].rows(),info.dofsB[np]+info.dofsI[np]);
    gsVector<index_t> nonZerosPerCol(info.dofsB[np]+info.dofsI[np]);
    for(index_t i=0; i<info.dofsB[np]+info.dofsI[np];i++)
    {
        int nC=0;
        if(m_globIsBoundIndex[np][i])
            nC=m_C[np].col(m_glob2BoundInteriorIndex[np][i]).nonZeros();
        nonZerosPerCol(i) = cast<double,int>(nC*1.333);
    }

    m_exC[np].reserve(nonZerosPerCol);
    for (index_t k=0; k<m_C[np].outerSize(); ++k)
        for (typename gsSparseMatrix<T>::InnerIterator it(m_C[np],k); it; ++it)
        {
            int c=getComp(np,it.col());
            index_t col = m_locDofsMapper[np][c].index(m_boundDofs[np][it.col()]);
            m_exC[np].insert(it.row(),col) = it.value();
        }



    if(false)
    {
       // int degreeDiff = m_basis[0][np].degree(0) - coarseBasis.degree(0);
       // coarseBasis.degreeIncrease(degreeDiff);
       // coarseBasis.uniformRefine(1);
       // std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices_KC;

        //  gsGridHierarchy<>::buildByRefinement(gsMultiBasis<>(coarseBasis), bc, m_assembler->options(), m_options.KC_optionsMG.numLevels, m_options.KC_optionsMG.numRefine, 1)
        //.moveMultiBasesTo(MG_bases)
        //.moveTransferMatricesTo(transferMatrices_KC)
        //.clear();

        if(m_bc_loc[np].dirichletSides().size()==0)
        {
            // gsInfo<<"Mass: \n"<<genAss.assembleMass().toDense()<<"\n \n";
            gsStopwatch time;
            gsSparseMatrix<T> M;

            assembleParameterMassForTensorProductSpace(m_basis[0][np],m_bc_loc[np],M);
            m_timings.assemblingMass[np]+=time.stop();
          mg_ptr= getMultiGrid<T>(gsMultiBasis<>(m_basis[0][np]),
                    m_K[np]->topLeftCorner(m_basis[0][np].size(), m_basis[0][np].size())+m_IETIoptions.opt.getReal("Regularization")*M,bc,m_IETIoptions.opt.getGroup("MG_KC"));
        }
        else
        {
            //TODO: fixme: s not defined, and copy of K is made,
            size_t s = 0 ; //transferMatrices_KC.back().rows();
            mg_ptr= getMultiGrid<T>(gsMultiBasis<>(m_basis[0][np]),
                    m_K[np]->topLeftCorner(s,s),bc,m_IETIoptions.opt.getGroup("MG_KC"));
        }


    }
    else
    {
        //TODO: use coarsening or set it outside IETI algorithm
        coarseBasis.degreeDecrease(coarseBasis.maxDegree()-1);
        coarseBasis.uniformRefine();
        std::vector< gsSparseMatrix<real_t, RowMajor> > transferMatrices_KC;
        gsGridHierarchy<T>::buildByRefinement(gsMultiBasis<>(coarseBasis), bc, m_IETIoptions.opt.getGroup("MG_KC"))
            .moveMultiBasesTo(MG_bases)
            .moveTransferMatricesTo(transferMatrices_KC)
            .clear();

        gsPde<T>* pde = m_assembler->pde().restrictToPatch(np);
        gsAssembler<T>* A = m_assembler->create();
        const gsBasis<T>& P1_basis = MG_bases.back()[0];
        pde->boundaryConditions() = bc;
        A->initialize(*pde,P1_basis,m_assembler->options());
        A->assemble();
        gsSparseMatrix<T> matrix = A->matrix();
        //TODO: Add dG-Terms
        //-------------------
        gsMatrix<index_t> *actives1,*actives2;
        gsMatrix<index_t> activesExtra1, activesExtra2;

        gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
        gsVector<T> quWeights;         // Mapped weights

        unsigned evFlags(0);     // Evaluation flags for the Geometry map
        gsQuadRule<T> QuRule;

        gsPoissonHeterogeneousAssembler<T>* pAss = static_cast< gsPoissonHeterogeneousAssembler<T>*>(m_assembler);

        for ( typename gsMultiPatch<T>::const_iiterator it = m_patches.iBegin(); it != m_patches.iEnd(); ++it )
        {
            const boundaryInterface & bi =
                    ( m_basis[0][it->first() .patch].numElements(it->first() .side() ) <
                    m_basis[0][it->second().patch].numElements(it->second().side() ) ?
                        it->getInverse() : *it );

            const size_t patch1      = bi.first().patch;
            const size_t patch2      = bi.second().patch;

            if(!(patch1 == np || patch2 == np))
                continue;
            const gsAffineFunction<T> interfaceMap(m_patches.getMapForInterface(bi));
            gsVisitorDg2<T,1> dg =  pAss->visitorDg(bi); //makes a copy, but class is very small.

            const gsBasis<T> & B1 = m_basis[0][patch1];// (!) unknown 0
            const gsBasis<T> & B2 = m_basis[0][patch2];

            const int bSize1      = B1.numElements( bi.first() .side() );
            const int bSize2      = B2.numElements( bi.second().side() );

            const int ratio = bSize1 / bSize2;
            GISMO_ASSERT(bSize1 >= bSize2 && bSize1%bSize2==0, "DG assumes nested interfaces. Got bSize1="<<bSize1<<", bSize2="<<bSize2<<"." );

            // Initialize
            if(patch1==np)
                dg.initialize(P1_basis, P1_basis, QuRule, evFlags);
            else
                dg.initialize(P1_basis, P1_basis, QuRule, evFlags);

            // Initialize geometry evaluators
            typename gsGeometryEvaluator<T>::uPtr geoEval1(getEvaluator(evFlags, m_patches[patch1]));
            typename gsGeometryEvaluator<T>::uPtr geoEval2(getEvaluator(evFlags, m_patches[patch2]));

            // Initialize domain element iterators
            typename gsBasis<T>::domainIter domIt1 = B1.makeDomainIterator( bi.first() .side() );
            typename gsBasis<T>::domainIter domIt2 = B2.makeDomainIterator( bi.second().side() );

            int count = 0;
            // iterate over all boundary grid cells on the "left"
            for (; domIt1->good(); domIt1->next() )
            {
                count++;
                // Compute the quadrature rule on both sides
                QuRule.mapTo( domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights);
                interfaceMap.eval_into(quNodes1,quNodes2);

                // Perform required evaluations on the quadrature nodes
                if(patch1==np)
                    dg.evaluate(P1_basis, *geoEval1, P1_basis, *geoEval2, quNodes1, quNodes2);
                else
                    dg.evaluate(P1_basis, *geoEval1, P1_basis, *geoEval2, quNodes1, quNodes2);

                dg.assemble(*domIt1,*domIt2, *geoEval1, *geoEval2, quWeights);
                dg.getActives(actives1, actives2);
                activesExtra1.setZero(0,0);
                activesExtra2.setZero(0,0);
                gsMatrix<T> dummy(matrix.rows(),1);
                if(patch1==np)
                    dg.localToGlobalIETI(A->system().colMapper(0),m_dirDofs[patch1],
                                         activesExtra1, matrix, dummy);
                else
                {
                    dg.revert();

                    dg.localToGlobalIETI(A->system().colMapper(0),m_dirDofs[patch2],
                                         activesExtra2, matrix, dummy);
                }

                if ( count % ratio == 0 ) // next master element ?
                    domIt2->next();
            }
        }

        //-------------------

        gsMultiGridOp<>::Ptr mg_P1;
        if(m_bc_loc[np].dirichletSides().size()==0)
        {
            // gsInfo<<"Mass: \n"<<genAss.assembleMass().toDense()<<"\n \n";
            gsStopwatch time;
            gsSparseMatrix<T> M;

            assembleParameterMassForTensorProductSpace(MG_bases.back()[0],m_bc_loc[np],M);
            m_timings.assemblingMass[np]+=time.stop();
            mg_P1= gsMultiGridOp<>::make(matrix+m_IETIoptions.opt.getReal("Regularization")*M, give(transferMatrices_KC));
        }
        else
        {
            //    size_t s = m_transferMatrices_KC[np].back().rows();
            mg_P1= gsMultiGridOp<>::make(matrix, give(transferMatrices_KC));
        }
        mg_P1->setOptions(m_IETIoptions.opt.getGroup("MG_KC"));
        gsKroneckerOp<>::Ptr massInv;
        assembleParameterMassInverseForTensorProductSpace(m_basis[0][np],bc,massInv);
        gsInfo << "done." << std::endl;

        gsInfo << "Assemble rectangular mass..." << std::flush;
        gsSparseMatrix<real_t> Mrect;
        assembleParameterMassForTensorProductSpace<>( m_basis[0][np], MG_bases.back()[0], bc, Mrect);
        gsInfo << "done, " << Mrect.rows() << "x" << Mrect.cols() << " dofs." << std::endl;

        unsigned m= massInv->rows();
        int elimT = 0;
        for(size_t i=0; i<m_patchISides[np].size()/2;++i)
        {
            patchSide ps = m_patchISides[np][i];
            int idx =0;
            for(int j=1; j< ps.index();j++)
                if(m_numExtraBasis[np][j]!=0) idx++;
            int elim = 0;
            for(int j=0; j<m_numExtraBasis[ps.patch][ps.index()];j++)
                if(!m_locDofsMapper[ps.patch][0].is_free(m_basis[0].size(np) +dgOffset(ps.patch,ps.side())+j))
                    elim++;
            elimT+=elim;
        }
        gsPreconditionerOp<>::Ptr ps;
        gsSparseMatrix<T> M;
        typename gsSparseMatrix<T>::Block K_block = m_K[np]->block(0,0,m,m);

        if(m_bc_loc[np].dirichletSides().size()==0)
        {
            assembleParameterMassForTensorProductSpace(m_basis[0][np],m_bc_loc[np],M);
            typename gsSparseMatrix<T>::Ptr Kreg= memory::make_shared(new gsSparseMatrix<T>(K_block+ m_IETIoptions.opt.getReal("Regularization")*M));
            ps= gsCompositePrecOp<>::make(
                        gsPreconditionerFromOp<>::make( makeMatrixOp(  Kreg ), makeSubspaceCorrectedMassSmootherOperator(m_basis[0][np], m_IETIoptions.opt.getReal("MG_KC.TL.Damping"), bc)),
                    makeGaussSeidelOp( Kreg )
                    );
            mg_ptr= gsTwoLevel::make(
                        makeMatrixOp(  Kreg ),
                        ps,
                        give(Mrect),
                        massInv,
                        mg_P1);
        }
        else
        {
            ps= gsCompositePrecOp<>::make(
                        gsPreconditionerFromOp<>::make( makeMatrixOp(  K_block ), makeSubspaceCorrectedMassSmootherOperator(m_basis[0][np], m_IETIoptions.opt.getReal("MG_KC.TL.Damping"), bc)),
                    makeGaussSeidelOp( K_block )
                    );
            mg_ptr= gsTwoLevel::make(
                        K_block,
                        ps,
                        give(Mrect),
                        massInv,
                        mg_P1);
        }

        static_cast<gsTwoLevel*>(mg_ptr.get())->setOptions(m_IETIoptions.opt.getGroup("MG_KC.TL"));

        delete A;
    }

    //setup Precond for dG-extra parts
    int elimT = 0;
    unsigned n= mg_ptr->rows();

    for(size_t i=0; i<m_patchISides[np].size()/2;++i)
    {
        patchSide ps = m_patchISides[np][i];
        int idx =0;
        for(int j=1; j< ps.index();j++)
            if(m_numExtraBasis[np][j]!=0) idx++;
        int elim = 0;
        for(int j=0; j<m_numExtraBasis[ps.patch][ps.index()];j++)
            if(!m_locDofsMapper[ps.patch][0].is_free(m_basis[0].size(np) +dgOffset(ps.patch,ps.side())+j))
                elim++;
        elimT+=elim;
    }
    /*
    {
        for(int k=0; k<ps.index();++k)
        {
            for(int j=0; j<m_numExtraBasis[ps.patch][k];j++)
                if(!m_locDofsMapper[ps.patch][0].is_free(m_basis[0].size(np) +dgOffset(np,k)+j))
                    elimc++;
        }
        int nL = m_numExtraBasis[ps.patch][ps.index()]-elim;
        gsSparseMatrix<T> mat = m_K[np].block(n + dgOffset(np,ps)-elimc,n + dgOffset(np,ps)-elimc,m_numExtraBasis[np][ps.side().index()]-elim,m_numExtraBasis[np][ps.side().index()]-elim);
        dGPrec->addOperator(idx,idx,makeSparseCholeskySolver(mat));
        if(np==0)
        {
            if(i==0)
                gsInfo<<"m_K:\n"<<m_K[np].toDense()<<"\n\n";
            gsInfo<<"1-1-loc block: \n"<<m_K[np].block(n + dgOffset(np,ps)-elimc,n + dgOffset(np,ps)-elimc,m_numExtraBasis[np][ps.side().index()]-elim,m_numExtraBasis[np][ps.side().index()]-elim).toDense()
                    <<"\n\n";
        }
    }
*/

    gsMatrix<T> schur =  m_K[np]->bottomRightCorner(m_numExtraBasisTotal[np]-elimT,m_numExtraBasisTotal[np]-elimT);
    if(false && np==5)
    {
        gsInfo<<"0-0 block:\n"<<m_K[np]->topLeftCorner(n,n).toDense()<<"\n";
        gsInfo<<"1-1 block:\n"<<m_K[np]->bottomRightCorner(m_numExtraBasisTotal[np]-elimT,m_numExtraBasisTotal[np]-elimT).toDense()<<"\n";
        gsInfo<<"1-0 block:\n"<<m_K[np]->block(n,0,m_numExtraBasisTotal[np]-elimT,n).toDense()<<"\n";
        gsInfo<<"0-1 block:\n"<<m_K[np]->block(0,n,n,m_numExtraBasisTotal[np]-elimT).toDense()<<"\n";
        gsInfo<<"C block:\n"<<m_exC[np]<<"\n";
    }
    gsMatrix<T> col;

    for(index_t j=0; j<m_numExtraBasisTotal[np]-elimT;++j)
    {
        col.setZero(n,1);
        for( int p=0; p<1;++p )
            mg_ptr->step(m_K[np]->block(0,n+j,n,1).toDense(), col);
        schur.col(j) -= m_K[np]->block(n,0,m_numExtraBasisTotal[np]-elimT,n)*col;
    }

    typename gsSparseMatrix<T>::Block Koff=m_K[np]->block(n,0,m_numExtraBasisTotal[np]-elimT,n);
    typename gsSparseMatrix<T>::Block KoffT= m_K[np]->block(0,n,n,m_numExtraBasisTotal[np]-elimT);

    typename gsBlockOp<T>::Ptr op1 = gsBlockOp<T>::make(2,2);
    op1->addOperator(0,0,gsIdentityOp<T>::make(n));
    op1->addOperator(1,1,gsIdentityOp<T>::make(m_numExtraBasisTotal[np]-elimT));
    op1->addOperator(0,1,gsProductOp<T>::make(gsScaledOp<T>::make(makeMatrixOp(KoffT),-1),mg_ptr));
    typename gsBlockOp<T>::Ptr op2 = gsBlockOp<T>::make(2,2);
    op2->addOperator(0,0,mg_ptr);
    op2->addOperator(1,1,makePartialPivLUSolver(schur));
    typename gsBlockOp<T>::Ptr op3 = gsBlockOp<T>::make(2,2);
    op3->addOperator(0,0,gsIdentityOp<T>::make(n));
    op3->addOperator(1,1,gsIdentityOp<T>::make(m_numExtraBasisTotal[np]-elimT));
    op3->addOperator(1,0,gsProductOp<T>::make(mg_ptr,gsScaledOp<T>::make(makeMatrixOp(Koff),-1)));


    typename gsProductOp<T>::Ptr Kprec = gsProductOp<T>::make(op3, op2,op1);

    m_inexact_KC[np] = gsBlockOp<T>::make(2,2);
    //  typename gsBlockOp<T>::Ptr Kprec = gsBlockOp<T>::make(2,2);
    //  Kprec->addOperator(0,0,mg_ptr);
    // Kprec->addOperator(1,1,dGPrec);//gsIdentityOp<T>::make(m_K[np].cols()-mg_ptr->cols()));
    //  Kprec->addOperator(1,1,makePartialPivLUSolver(schur));

    m_inexact_KC[np]->addOperator(0,0,Kprec);

    getSchur(m_S[np],np);
    m_inexact_KC[np]->addOperator(1,1,makePartialPivLUSolver(m_S[np]));

    delete &coarseBasis;

}

template<class T>
void gsIETIdGAssembler<T>::getSchur(gsMatrix<T> &schur, size_t np)
{
    /*
    gsMatrix<T> d(m_Kbb[np].rows(),1);
    d.array()= m_Kbb[np].diagonal().array();
    gsSparseMatrix<T> d_inv(m_Kbb[np].rows(),m_Kbb[np].cols());
    d_inv.setIdentity();
    for(index_t i =0;i<m_Kbb[np].rows();i++)
        d_inv(i,i)/=(d(i));

    schur = (m_C[np]*d_inv*m_C[np].transpose().toDense());
*/
    /*
    gsMatrix<T> d(m_K[np].rows(),1);
    d.array()= m_K[np].diagonal().array();
    gsSparseMatrix<T> d_inv(m_K[np].rows(),m_K[np].cols());
    d_inv.setIdentity();
    for(size_t i =0;i<m_K[np].rows();i++)
        d_inv(i,i)/=d(i);

    schur = (m_exC[np]*d_inv*m_exC[np].transpose().toDense());
*/
    /*
    sparseSPDfact Kbb_LU(m_Kbb[np]);
    schur = (m_C[np]*Kbb_LU.solve(m_C[np].transpose().toDense()));
*/

    /*
    gsMatrix<T> d(m_Kii[np].rows(),1);
    d.array()= m_Kii[np].diagonal().array();
    gsSparseMatrix<T> d_inv(m_Kii[np].rows(),m_Kii[np].cols());
    d_inv.setIdentity();
    for(size_t i =0;i<m_Kii[np].rows();i++)
        d_inv(i,i)/=d(i);
    sparseLUfact Kbb_LU(m_Kbb[np]-m_Kib[np].transpose()*d_inv*m_Kib[np]);
    schur = (m_C[np]*Kbb_LU.solve(m_C[np].transpose().toDense()));
*/

    /*
    sparseSPDfact kiiinv(m_Kii[np]);
           // gsGenericAssembler<T> genAss(m_patches.patch(np), m_basis.front()[np],gsAssemblerOptions(),&m_bc_loc[np]);
    Eigen::PartialPivLU<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > Sb(m_Kbb[np].toDense()-m_Kib[np].transpose().toDense()*kiiinv.solve(m_Kib[np].toDense()));
   schur = m_C[np]*Sb.solve(m_C[np].transpose().toDense());
*/
    /*
    gsGenericAssembler<T> genAss(m_patches.patch(np), m_basis.front()[np],gsAssemblerOptions(),&m_bc_loc[np]);
    sparseSPDfact kinv(m_K[np]+1.e-5*genAss.assembleMass());
    schur = m_exC[np]*kinv.solve(m_exC[np].transpose().toDense());
    */


    gsMatrix<T> col;
    schur.setZero(info.dofsP[np],info.dofsP[np]);
    gsLinearOperator<T>* BlockMg = m_inexact_KC[np]->getOperator(0,0).get();
    for(index_t i=0; i<info.dofsP[np];++i)
    {
        col.setZero(m_exC[np].cols(),1);
        BlockMg->apply(m_exC[np].row(i).transpose().toDense(), col);
        schur.col(i) = m_exC[np]*col;
    }

    gsInfo<<"Finished schur\n";

}


} // namespace gismo
