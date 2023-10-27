

#include <gismo.h>

#include <gsRemappedBasis/gsHLR.h>
#include <gsRemappedBasis/gsTHB.h>
#include <gsRemappedBasis/otherUtils.h>

using namespace gismo;


int main()
{
    const int deg=2;

    gsBoxList boxes(2);

    std::vector<gsBoxList::basisPtr> bases;

    gsKnotVector<> knots(0,1,3,deg+1,1);
    gsKnotVector<> knots2=knots;
    knots2.uniformRefine();

    bases.push_back(gsBoxList::basisPtr(new gsTensorBSplineBasis<2>(knots,knots)));
    bases.push_back(gsBoxList::basisPtr(new gsTensorBSplineBasis<2>(knots2,knots2)));

    gsMatrix<> box(2,2);
    box<<knots2[4],knots[5],
         knots2[0],knots2[5];
    boxes.append(box,1);

    gsHLR<2> space(bases,boxes);
    space.exportDefinitionToTex("defDom");
    space.exportRequestToTex("reqDom");
    space.exportSelectorToTex ("repDom");

    /* // Note: this requires a proper latex installation
    int pdfLatexErr=0;
    pdfLatexErr=system("pdflatex defDom.tex >/dev/null");
    pdfLatexErr=system("pdflatex reqDom.tex >/dev/null");
    pdfLatexErr=system("pdflatex repDom.tex >/dev/null");
    pdfLatexErr=system("pdfmerger defDom.pdf reqDom.pdf repDom.pdf domains.pdf");
    GISMO_UNUSED(pdfLatexErr);
    */

    printAllFunctions(space,"basisFuncHLR");

    gsInfo<<"OK\n";

    return 0;
}


