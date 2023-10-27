/** @file gsDecoupledBasis.cpp

    @brief Tests partition of unity property of the decoupled basis
    and tries several exporting functions.

    Furthermore, functionality for exporting figures for the revised
    version of @cite bm2016 are saved here.

    This file is a part of the G+Smo library, in particular of the
    ambitious gsRemappedBasis project.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Authors: Dominik Mokris, Andrea Bressan
*/


#include <gismo.h>

#include <gsRemappedBasis/gsDecoupledBasis.h>
#include <gsRemappedBasis/gsTHB.h>
#include <gsRemappedBasis/otherUtils.h>

using namespace gismo;

bool partitionOfUnityTest()
{
    const int deg=3;
    const int maxLvl=2;
    std::vector<gsKnotVector<> > knotsx(1,gsKnotVector<>(0,1,2,deg+1,1));
    std::vector<gsKnotVector<> > knotsy(1,gsKnotVector<>(0,1,2,deg+1,1));

    std::vector<gsBoxList::basisPtr> bases;
    for (int l=0;l<maxLvl+1;++l)
    {
        bases.push_back((gsBoxList::basisPtr(new gsTensorBSplineBasis<2>(knotsx.back(),knotsy.back()))));
        knotsx.push_back(knotsx.back());
        knotsx.back().uniformRefine();
        knotsy.push_back(knotsy.back());
        knotsy.back().uniformRefine();
    }


#define FromBack(a,b)( a[a.size()-(b)-1] )
    gsBoxList  boxes(2);
    gsMatrix<> box(2,2);


    // lvl 2 area
    box<<knotsx[0][deg+1],knotsx[0][deg+2], knotsy[0][deg],FromBack(knotsy[0],deg);
    boxes.append(box,2);

    box<<knotsx[0][deg],FromBack(knotsx[0],deg), knotsy[0][deg+1],knotsy[0][deg+2];
    boxes.append(box,2);



    gsDecoupledBasis<2> decBasis(bases,boxes);

    printAllFunctions(decBasis,"decoupled",11000);
    decBasis.debugWithTeX();
    decBasis.exportSelectorToTex("remapped");

    return decBasis.testReprMatrix();
}

// This function is for the export of figures for the revised version
// of the paper A. Bressan, D. Mokris: A versatile strategy for the
// implementation of adaptive splines. Please, do not meddle with it
// unless you are one of the authors.
bool forThePaper()
{
    const int deg=4;
    const int maxLvl=1;
    std::vector<gsKnotVector<> > knotsx(1,gsKnotVector<>(0,1,4,deg+1,1));
    std::vector<gsKnotVector<> > knotsy(1,gsKnotVector<>(0,1,4,deg+1,1));

    std::vector<gsBoxList::basisPtr> bases;
    for (int l=0;l<maxLvl+1;++l)
    {
        bases.push_back((gsBoxList::basisPtr(new gsTensorBSplineBasis<2>(knotsx.back(),knotsy.back()))));
        knotsx.push_back(knotsx.back());
        knotsx.back().uniformRefine(1,2);
        knotsy.push_back(knotsy.back());
        knotsy.back().uniformRefine(1,2);
    }
#define FromBack(a,b)( a[a.size()-(b)-1] )
    gsBoxList  boxes(2);
    gsMatrix<> box(2,2);

    box<<knotsx[0][deg+1],knotsx[0][deg+3], knotsy[0][deg],FromBack(knotsy[0],deg);
    boxes.append(box,1);

    box<<knotsx[0][deg],knotsx[0][deg+1],knotsy[0][deg+2],knotsy[0][deg+4];
    boxes.append(box,1);

    box<<knotsx[0][deg+3],knotsx[0][deg+5],knotsy[0][deg+1],knotsy[0][deg+3];
    boxes.append(box,1);

    gsDecoupledBasis<2> DHB(bases,boxes);
    gsTHB<2> THB(bases,boxes);
    gsTHB<2> HB(bases,boxes,false);
    
    printAllFunctions(DHB,"remappedDHB",11000);
    printAllFunctions(THB,"remappedTHB",11000);
    printAllFunctions(HB,"remappedHB",11000);
    DHB.debugWithTeX();
    DHB.exportSelectorToTex("remapped");

    /* Exporting in ParaView: open the file (remapped*.vts), Apply,
       Filters->Alphabetical->Warp by scalar and choose the function
       you want to plot (nr. 40 in HB and THB, nrs. 60-63 in
       DHB). Contour colouring can be achieved by selecting the same
       function for colouring instead of "Solid color".
       
       The camera is saved in the revision folder (different from the
       submission camera!), it can be loaded with the "Adjust camera"
       button right next to 2D/3D switch of the plot. And the window
       is made smaller from the left to avoid the margins. Finally,
       export scene as a pdf.*/
    
    return DHB.testReprMatrix();
}

int main(int argc, char *argv[])
{
    bool partionOfUnity = true;
    bool figures = false;

    gsCmdLine cmd("Tests partition of unity property of the decoupled basis "
        "and tries several exporting functions.");
    cmd.addSwitch("figures", "Export figures for the revised version of @cite bm2016.", figures);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    bool result = true;
    
    if (partionOfUnity)
    {
        gsInfo << "PartitionOfUnityTest..." << std::flush;
        result &= partitionOfUnityTest();
        gsInfo << "done.\n";
    }

    if (figures)
    {
        gsInfo << "ForThePaper..." << std::flush;
        result &= forThePaper();
        gsInfo << "done.\n";
    }
    else
    {
        gsInfo << "To export figures for the revised version of @cite bm2016, call --figures.\n";
    }
    
    return result ? 0 : 1;
}
