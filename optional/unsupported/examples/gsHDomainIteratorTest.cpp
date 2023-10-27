// Author Jaka Speh
// testing hierarchical domain iterator
// work in progress
// Usage
// 1.) to run all tests (see if the area of all elements is equal to 1)
// $ gsHDomainIteratorTest
//
// 2.) to run test from a file, prints some data from each element
// $ gsHDomainIteratorTest -b filnameOfBasis -d dimension



#include <iostream>
#include <string>

#include <gismo.h>
#include <gismo_dev.h>
#include <gsHSplines/gsHDomainIterator.h>



using namespace gismo;

// gsHDomainIterator::nextBox is now a private function
//void testBoxIteration(gsHDomainIterator<2, real_t>& domIter)
//{
//    int numBox = 0;
//    while (domIter.nextBox())
//    {
//        gsInfo << numBox << ". box\n"
//                  << "level: " << domIter.m_cnode->level << "\n"
//                  << "low: " << domIter.m_cnode->lower.transpose() << "\n"
//                  << "upp: " << domIter.m_cnode->upper.transpose() << "\n"
//                  << "\n";

//        numBox++;
//    }

//    numBox = 0;
//    while (domIter.nextBox())
//    {
//        gsInfo << numBox << ". box\n"
//                  << "level: " << domIter.m_cnode->level << "\n"
//                  << "low: " << domIter.m_cnode->lower.transpose() << "\n"
//                  << "upp: " << domIter.m_cnode->upper.transpose() << "\n"
//                  << "\n";

//        numBox++;
//    }
//}

template<typename domainIter>
void testNextIteration(domainIter& domIter)
{
    int numElement = 0;

    for (; domIter.good(); domIter.next(), numElement++)
    {
        gsInfo << numElement << ". element\n"
                  << "level: " << domIter.getLevel() << "\n"
                  << "low: " << domIter.lowerCorner().transpose() << "\n"
                  << "upp: " << domIter.upperCorner().transpose() << "\n"
                  << "center: " << domIter.center.transpose() << "\n"
                  << "\n";
    }
}


template <unsigned d>
int testingDomain(gsHDomainIterator<real_t, d>& domIter)
{
    real_t sum = 0;

    for (; domIter.good(); domIter.next())
    {
        gsVector<real_t, d> low; // lower corner
        gsVector<real_t, d> upp; // upper corner
        gsVector<real_t, d> ell; // upp - low


        low = domIter.lowerCorner();
        upp = domIter.upperCorner();
        ell = upp - low;

        sum += ell.prod();
    }

    if ( math::almostEqual<real_t>(1.0, sum) )
    {
        gsInfo << "Sum of all elements is " << sum << "\n" << "\n";
        return 0;
    }
    else
    {
        gsInfo << "ERROR: sum of all elements is " << sum
                  << " instead of 1.\n" << "\n";
        return 1;
    }
}


int main(int argc, char *argv[])
{

    std::string baseFile("");
    index_t d = 0;
    gsCmdLine cmd("Test for hiearchical domain iterator (gsHDomainIterator.h)");
    cmd.addString("b", "basis",
                  "Input file for (truncated) hierarchical basis.",
                  baseFile);
    cmd.addInt("d", "dimension", "Dimension of the given example", d);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (baseFile != "")
    {
        if (d == 2)
        {
            gsHBSplineBasis<2>::uPtr thbs = gsReadFile<>(baseFile);

            if (!thbs)
            {
                gsInfo << "Could not find HB-spline basis, trying to find "
                             "THB-spline basis." << "\n";

                gsTHBSplineBasis<2>::uPtr t = gsReadFile<>(baseFile);
                if (!t)
                {
                    gsInfo << "Invalid basis!";
                    return 0;
                }
                else
                {
                    gsHDomainIterator<real_t, 2> domIter(*t);
                    testNextIteration<gsHDomainIterator<real_t, 2> >(domIter);
                }
            }
            else
            {
                gsInfo<< *thbs ;
                gsHDomainIterator<real_t, 2> domIter(*thbs);
                testNextIteration<gsHDomainIterator<real_t, 2> >(domIter);
            }
        }
        else if (d == 3)
        {

            gsHBSplineBasis<3>::uPtr thbs = gsReadFile<>(baseFile);

            if (!thbs)
            {
                gsTHBSplineBasis<3>::uPtr t = gsReadFile<>(baseFile);
                if (!t)
                {
                    gsInfo << "Invalid basis!";
                    return 0;
                }
                else
                {
                    gsHDomainIterator<real_t, 3> domIter(*t);
                    testNextIteration<gsHDomainIterator<real_t, 3> >(domIter);
                }
            }
            else
            {
                gsHDomainIterator<real_t, 3> domIter(*thbs);
                testNextIteration<gsHDomainIterator<real_t, 3> >(domIter);
            }
        }
    }
    else
    {
        gsInfo << "*******************************************************\n"
                     "** 2d cases                                          **\n"
                     "*******************************************************\n\n";

        const char* files2D[] =
        {    "basis_thbs.xml",
             "basis_thbs_01.xml",
             "basis_thbs_02.xml",
             "basis_thbs_03.xml",             
             "basis_thbs_04.xml",
             "basis_thbs_05.xml",
             "basis_thbs_06.xml",
             "basis_thbs_07.xml",
        };

        int failed = 0;

        for (unsigned i = 0; i < sizeof(files2D) / sizeof(files2D[0]); ++i)
        {
            gsTHBSplineBasis<2>::uPtr thbs;

            std::string filename(files2D[i]);
            gsInfo << filename << "\n";
            gsFileData<> data(filename);
            thbs = data.getFirst< gsTHBSplineBasis<2> >();
            gsHDomainIterator<real_t, 2> domIter(*thbs);

            failed += testingDomain<2>(domIter);

        }

        gsInfo << "\n\n\n"
                     "*******************************************************\n"
                     "** 3d cases                                          **\n"
                     "*******************************************************\n\n";

        const char* files3D[] =
        {    //"basis3d/basis_thbs_01.xml",
             //"basis3d/basis_thbs_02.xml",
             "basis3d/basis_thbs_03.xml"
        };

        for (unsigned i = 0; i < sizeof(files3D) / sizeof(files3D[0]); ++i)
        {
            gsTHBSplineBasis<3>::uPtr thbs;

            std::string filename(files3D[i]);
            gsInfo << filename << "\n";
            gsFileData<> data(filename);
            thbs = data.getFirst< gsTHBSplineBasis<3> >();
            gsHDomainIterator<real_t, 3> domIter(*thbs);

            failed += testingDomain<3>(domIter);

        }


        return failed;
    }

    return 0;
}

