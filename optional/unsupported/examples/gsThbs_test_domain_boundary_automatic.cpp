/* Automatic testing of gsHTensorBasis::domain_boundaries(...) and its evil minions.
 * Seven test cases are provided and answer to each of them is compared to the known answer.
 */

#include <iostream>
#include <fstream>

#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;

bool examine_xml( std::string filename_xml, std::string corr_answ )
{
    /* Takes *.xml file, reads it, gets the polylines and compares
     * them with the corr_answ. */
    gsTHBSpline<2>::uPtr THB;
    std::string filename = "basis2d/" + filename_xml;
    gsInfo<<filename<<"\n";
    gsFileData< real_t > data(filename);
    THB = data.getFirst< gsTHBSpline<2> >();

    std::stringstream ans_str;
    std::vector< std::vector<std::vector< std::vector< real_t > > > > polylines;
    std::vector<std::vector< std::vector<index_t> > > aabb;
    aabb = THB -> basis().domainBoundariesParams(polylines);

    // Now put everything into stringstream.
    for(unsigned int i0 = 0; i0 < polylines.size();i0++)
    {
        ans_str << "Level:"<< i0<<":";
        for(unsigned int i1 = 0; i1 < polylines[i0].size();i1++)
	{
            ans_str<<"ConComp:"<<i1<<":";
            ans_str<<"BoundBox:("<< aabb[i0][i1][0]<<","<< aabb[i0][i1][1]<<","<< aabb[i0][i1][2]<<","<< aabb[i0][i1][3] <<");";
            for(unsigned int i2 = 0; i2 < polylines[i0][i1].size();i2++)
	    {
                ans_str<<"Line:"<<i2<<":";
                ans_str<< "(";
                for(unsigned int i3 = 0; i3 < polylines[i0][i1][i2].size();i3++)
		{
                    ans_str << polylines[i0][i1][i2][i3];
                    if( i3+1 != polylines[i0][i1][i2].size() )
                        ans_str << ",";
                }
                ans_str <<");";
            }
        }
    }

    // Check the stringstream with the corr_answ. If they are different, shout.
    std::string ans = ans_str.str();
    //gsInfo << answer << "\n";
    if( ans != corr_answ )
    {
        gsInfo << "Warning, incorrect answer for the file " << filename_xml << "." << "\n";
        //std::fstream fout;
        std::ostream & fout = gsInfo;
        //fout.open( "gsThbs_test_domain_boundary.log", std::fstream::in );
        fout << "Your answer:" << "\n";
        fout << ans << "\n" << "\n";
        fout << "Correct answer:" << "\n";
        fout << corr_answ << "\n";
        //fout.close();

        return false;
    }
    else
        return true;
}


int main(int argc, char *argv[])
{
    if(
            ( examine_xml( "test_domain_boundary_00.xml", "Level:0:ConComp:0:BoundBox:(0,0,24,24);Line:0:(1,0,1,0.3333);Line:1:(0,0,1,0);Line:2:(0,0,0,1);Line:3:(0,1,1,1);Line:4:(1,0.5,1,1);Line:5:(0.6667,0.5,1,0.5);Line:6:(0.6667,0.5,0.6667,0.8333);Line:7:(0.3333,0.8333,0.6667,0.8333);Line:8:(0.3333,0.6667,0.3333,0.8333);Line:9:(0.1667,0.6667,0.3333,0.6667);Line:10:(0.1667,0.3333,0.1667,0.6667);Line:11:(0.1667,0.3333,0.3333,0.3333);Line:12:(0.3333,0.1667,0.3333,0.3333);Line:13:(0.3333,0.1667,0.6667,0.1667);Line:14:(0.6667,0.1667,0.6667,0.3333);Line:15:(0.6667,0.3333,1,0.3333);Level:1:ConComp:0:BoundBox:(4,4,16,16);Line:0:(0.3333,0.1667,0.3333,0.3333);Line:1:(0.1667,0.3333,0.3333,0.3333);Line:2:(0.1667,0.3333,0.1667,0.6667);Line:3:(0.1667,0.6667,0.3333,0.6667);Line:4:(0.3333,0.5,0.3333,0.6667);Line:5:(0.3333,0.5,0.5,0.5);Line:6:(0.5,0.3333,0.5,0.5);Line:7:(0.5,0.3333,0.6667,0.3333);Line:8:(0.6667,0.1667,0.6667,0.3333);Line:9:(0.3333,0.1667,0.6667,0.1667);ConComp:1:BoundBox:(8,8,24,20);Line:0:(0.6667,0.3333,0.6667,0.5);Line:1:(0.5,0.5,0.6667,0.5);Line:2:(0.5,0.5,0.5,0.6667);Line:3:(0.3333,0.6667,0.5,0.6667);Line:4:(0.3333,0.6667,0.3333,0.8333);Line:5:(0.3333,0.8333,0.6667,0.8333);Line:6:(0.6667,0.5,0.6667,0.8333);Line:7:(0.6667,0.5,1,0.5);Line:8:(1,0.3333,1,0.5);Line:9:(0.6667,0.3333,1,0.3333);Level:2:ConComp:0:BoundBox:(8,8,16,16);Line:0:(0.5,0.3333,0.5,0.5);Line:1:(0.3333,0.5,0.5,0.5);Line:2:(0.3333,0.5,0.3333,0.6667);Line:3:(0.3333,0.6667,0.5,0.6667);Line:4:(0.5,0.5,0.5,0.6667);Line:5:(0.5,0.5,0.6667,0.5);Line:6:(0.6667,0.3333,0.6667,0.5);Line:7:(0.5,0.3333,0.6667,0.3333);" ) )
            &&
            ( examine_xml( "test_domain_boundary_01.xml", "Level:0:ConComp:0:BoundBox:(0,2,2,4);Line:0:(0,0.2,0,0.4);Line:1:(0,0.4,0.5,0.4);Line:2:(0.5,0.2,0.5,0.4);Line:3:(0,0.2,0.5,0.2);ConComp:1:BoundBox:(0,6,2,8);Line:0:(0,0.6,0,0.8);Line:1:(0,0.8,0.5,0.8);Line:2:(0.5,0.6,0.5,0.8);Line:3:(0,0.6,0.5,0.6);Level:1:ConComp:0:BoundBox:(0,0,4,10);Line:0:(0,0,0,0.2);Line:1:(0,0.2,0.5,0.2);Line:2:(0.5,0.2,0.5,0.4);Line:3:(0,0.4,0.5,0.4);Line:4:(0,0.4,0,0.6);Line:5:(0,0.6,0.5,0.6);Line:6:(0.5,0.6,0.5,0.8);Line:7:(0,0.8,0.5,0.8);Line:8:(0,0.8,0,1);Line:9:(0,1,1,1);Line:10:(1,0,1,1);Line:11:(0,0,1,0);" ) )
            &&
            ( examine_xml( "test_domain_boundary_02.xml", "Level:0:ConComp:0:BoundBox:(2,0,10,10);Line:0:(0.8,0,0.8,0.2);Line:1:(0.2,0.2,0.8,0.2);Line:2:(0.2,0.2,0.2,0.4);Line:3:(0.2,0.4,0.6,0.4);Line:4:(0.6,0.4,0.6,0.6);Line:5:(0.2,0.6,0.6,0.6);Line:6:(0.2,0.6,0.2,0.8);Line:7:(0.2,0.8,0.4,0.8);Line:8:(0.4,0.8,0.4,1);Line:9:(0.4,1,1,1);Line:10:(1,0,1,1);Line:11:(0.8,0,1,0);Level:1:ConComp:0:BoundBox:(0,0,8,10);Line:0:(0.2,0.2,0.2,0.4);Line:1:(0.2,0.4,0.6,0.4);Line:2:(0.6,0.4,0.6,0.6);Line:3:(0.2,0.6,0.6,0.6);Line:4:(0.2,0.6,0.2,0.8);Line:5:(0.2,0.8,0.4,0.8);Line:6:(0.4,0.8,0.4,1);Line:7:(0,1,0.4,1);Line:8:(0,0,0,1);Line:9:(0,0,0.8,0);Line:10:(0.8,0,0.8,0.2);Line:11:(0.2,0.2,0.8,0.2);" ) )
            && // Warning: some of the domains in tests 03 and 07 get slightly enlarged to fulfill the condition that Omega^{\ell+1} is a union of cells from Omega^\ell.
            //( examine_xml( "test_domain_boundary_03.xml", "Level:0:ConComp:0:BoundBox:(2,0,10,10);Line:0:(0.8,0,0.8,0.2);Line:1:(0.2,0.2,0.8,0.2);Line:2:(0.2,0.2,0.2,0.4);Line:3:(0.2,0.4,0.6,0.4);Line:4:(0.6,0.4,0.6,0.6);Line:5:(0.2,0.6,0.6,0.6);Line:6:(0.2,0.6,0.2,0.8);Line:7:(0.2,0.8,0.4,0.8);Line:8:(0.4,0.8,0.4,1);Line:9:(0.4,1,1,1);Line:10:(1,0,1,1);Line:11:(0.8,0,1,0);Level:1:ConComp:0:BoundBox:(0,0,16,20);Line:0:(0.2,0.2,0.2,0.4);Line:1:(0.2,0.4,0.6,0.4);Line:2:(0.6,0.4,0.6,0.6);Line:3:(0.2,0.6,0.6,0.6);Line:4:(0.2,0.6,0.2,0.8);Line:5:(0.2,0.8,0.4,0.8);Line:6:(0.4,0.8,0.4,1);Line:7:(0,1,0.4,1);Line:8:(0,0,0,1);Line:9:(0,0,0.8,0);Line:10:(0.8,0,0.8,0.2);Line:11:(0.2,0.2,0.8,0.2);" ) )
            //&& // Test 03 disabled for now (L-shaped domains have somehow crawled back).
            ( examine_xml( "test_domain_boundary_04.xml", "Level:0:ConComp:0:BoundBox:(8,8,16,16);Line:0:(0.6,0.4,0.6,0.6);Line:1:(0.4,0.6,0.6,0.6);Line:2:(0.4,0.6,0.4,0.8);Line:3:(0.4,0.8,0.6,0.8);Line:4:(0.6,0.6,0.6,0.8);Line:5:(0.6,0.6,0.8,0.6);Line:6:(0.8,0.4,0.8,0.6);Line:7:(0.6,0.4,0.8,0.4);ConComp:1:BoundBox:(0,0,20,20);Line:0:(0,0,0,1);Line:1:(0,1,0.2,1);Line:2:(0.2,0.2,0.2,1);Line:3:(0.2,0.2,1,0.2);Line:4:(1,0,1,0.2);Line:5:(0,0,1,0);ConComp:2:BoundBox:(16,16,20,20);Line:0:(0.8,0.8,0.8,1);Line:1:(0.8,1,1,1);Line:2:(1,0.8,1,1);Line:3:(0.8,0.8,1,0.8);Level:1:ConComp:0:BoundBox:(8,8,16,16);Line:0:(0.4,0.4,0.4,0.8);Line:1:(0.4,0.8,0.6,0.8);Line:2:(0.6,0.6,0.6,0.8);Line:3:(0.6,0.6,0.8,0.6);Line:4:(0.8,0.4,0.8,0.6);Line:5:(0.4,0.4,0.8,0.4);ConComp:1:BoundBox:(4,4,20,20);Line:0:(0.2,0.2,0.2,1);Line:1:(0.2,1,0.8,1);Line:2:(0.8,0.8,0.8,1);Line:3:(0.8,0.8,1,0.8);Line:4:(1,0.2,1,0.8);Line:5:(0.2,0.2,1,0.2);Level:2:ConComp:0:BoundBox:(8,8,12,12);Line:0:(0.4,0.4,0.4,0.6);Line:1:(0.4,0.6,0.6,0.6);Line:2:(0.6,0.4,0.6,0.6);Line:3:(0.4,0.4,0.6,0.4);" ) )
            &&
            ( examine_xml( "test_domain_boundary_05.xml", "Level:0:ConComp:0:BoundBox:(8,8,40,48);Line:0:(0.625,0.125,0.625,0.25);Line:1:(0.125,0.125,0.625,0.125);Line:2:(0.125,0.125,0.125,0.75);Line:3:(0.125,0.75,0.625,0.75);Line:4:(0.625,0.625,0.625,0.75);Line:5:(0.25,0.625,0.625,0.625);Line:6:(0.25,0.25,0.25,0.625);Line:7:(0.25,0.25,0.625,0.25);ConComp:1:BoundBox:(0,0,64,64);Line:0:(0.75,0.5,0.75,0.875);Line:1:(0.375,0.5,0.75,0.5);Line:2:(0.375,0.375,0.375,0.5);Line:3:(0.375,0.375,0.75,0.375);Line:4:(0.75,0,0.75,0.375);Line:5:(0,0,0.75,0);Line:6:(0,0,0,1);Line:7:(0,1,1,1);Line:8:(1,0.875,1,1);Line:9:(0.75,0.875,1,0.875);Level:1:ConComp:0:BoundBox:(8,8,40,48);Line:0:(0.625,0.125,0.625,0.25);Line:1:(0.125,0.125,0.625,0.125);Line:2:(0.125,0.125,0.125,0.75);Line:3:(0.125,0.75,0.625,0.75);Line:4:(0.625,0.625,0.625,0.75);Line:5:(0.25,0.625,0.625,0.625);Line:6:(0.25,0.25,0.25,0.625);Line:7:(0.25,0.25,0.625,0.25);Level:2:ConComp:0:BoundBox:(48,0,64,56);Line:0:(0.75,0,0.75,0.875);Line:1:(0.75,0.875,1,0.875);Line:2:(1,0,1,0.875);Line:3:(0.75,0,1,0);Level:3:ConComp:0:BoundBox:(24,24,48,32);Line:0:(0.375,0.375,0.375,0.5);Line:1:(0.375,0.5,0.75,0.5);Line:2:(0.75,0.375,0.75,0.5);Line:3:(0.375,0.375,0.75,0.375);" ) )
            &&
            ( examine_xml( "test_domain_boundary_06.xml", "Level:0:ConComp:0:BoundBox:(16,16,20,20);Line:0:(0.5,0.5,0.5,0.625);Line:1:(0.5,0.625,0.625,0.625);Line:2:(0.625,0.5,0.625,0.625);Line:3:(0.5,0.5,0.625,0.5);ConComp:1:BoundBox:(12,12,24,24);Line:0:(0.375,0.375,0.375,0.75);Line:1:(0.375,0.75,0.75,0.75);Line:2:(0.75,0.375,0.75,0.75);Line:3:(0.375,0.375,0.75,0.375);ConComp:2:BoundBox:(8,8,28,28);Line:0:(0.25,0.25,0.25,0.875);Line:1:(0.25,0.875,0.875,0.875);Line:2:(0.875,0.25,0.875,0.875);Line:3:(0.25,0.25,0.875,0.25);ConComp:3:BoundBox:(0,0,32,32);Line:0:(0,0,0,1);Line:1:(0,1,0.125,1);Line:2:(0.125,0.125,0.125,1);Line:3:(0.125,0.125,1,0.125);Line:4:(1,0,1,0.125);Line:5:(0,0,1,0);Level:1:ConComp:0:BoundBox:(16,16,20,20);Line:0:(0.5,0.5,0.5,0.625);Line:1:(0.5,0.625,0.625,0.625);Line:2:(0.625,0.5,0.625,0.625);Line:3:(0.5,0.5,0.625,0.5);ConComp:1:BoundBox:(12,12,24,24);Line:0:(0.375,0.375,0.375,0.75);Line:1:(0.375,0.75,0.75,0.75);Line:2:(0.75,0.375,0.75,0.75);Line:3:(0.375,0.375,0.75,0.375);Level:2:ConComp:0:BoundBox:(8,8,28,28);Line:0:(0.25,0.25,0.25,0.875);Line:1:(0.25,0.875,0.875,0.875);Line:2:(0.875,0.25,0.875,0.875);Line:3:(0.25,0.25,0.875,0.25);ConComp:1:BoundBox:(4,4,32,32);Line:0:(0.125,0.125,0.125,1);Line:1:(0.125,1,1,1);Line:2:(1,0.125,1,1);Line:3:(0.125,0.125,1,0.125);" ) )
            &&
            ( examine_xml( "test_domain_boundary_07.xml", "Level:0:ConComp:0:BoundBox:(0,0,8,20);Line:0:(0,0,0,1);Line:1:(0,1,0.2,1);Line:2:(0.2,0.2,0.2,1);Line:3:(0.2,0.2,0.4,0.2);Line:4:(0.4,0,0.4,0.2);Line:5:(0,0,0.4,0);ConComp:1:BoundBox:(8,8,16,16);Line:0:(0.6,0.4,0.6,0.6);Line:1:(0.4,0.6,0.6,0.6);Line:2:(0.4,0.6,0.4,0.8);Line:3:(0.4,0.8,0.6,0.8);Line:4:(0.6,0.6,0.6,0.8);Line:5:(0.6,0.6,0.8,0.6);Line:6:(0.8,0.4,0.8,0.6);Line:7:(0.6,0.4,0.8,0.4);ConComp:2:BoundBox:(12,0,20,4);Line:0:(0.6,0,0.6,0.2);Line:1:(0.6,0.2,1,0.2);Line:2:(1,0,1,0.2);Line:3:(0.6,0,1,0);ConComp:3:BoundBox:(16,16,20,20);Line:0:(0.8,0.8,0.8,1);Line:1:(0.8,1,1,1);Line:2:(1,0.8,1,1);Line:3:(0.8,0.8,1,0.8);Level:1:ConComp:0:BoundBox:(4,0,20,20);Line:0:(0.6,0,0.6,0.2);Line:1:(0.4,0,0.6,0);Line:2:(0.4,0,0.4,0.1);Line:3:(0.4,0.1,0.5,0.1);Line:4:(0.5,0.1,0.5,0.4);Line:5:(0.5,0.4,0.8,0.4);Line:6:(0.8,0.4,0.8,0.6);Line:7:(0.6,0.6,0.8,0.6);Line:8:(0.6,0.6,0.6,0.8);Line:9:(0.4,0.8,0.6,0.8);Line:10:(0.4,0.2,0.4,0.8);Line:11:(0.2,0.2,0.4,0.2);Line:12:(0.2,0.2,0.2,1);Line:13:(0.2,1,0.8,1);Line:14:(0.8,0.8,0.8,1);Line:15:(0.8,0.8,1,0.8);Line:16:(1,0.2,1,0.8);Line:17:(0.6,0.2,1,0.2);Level:2:ConComp:0:BoundBox:(8,2,12,12);Line:0:(0.5,0.1,0.5,0.4);Line:1:(0.4,0.1,0.5,0.1);Line:2:(0.4,0.1,0.4,0.6);Line:3:(0.4,0.6,0.6,0.6);Line:4:(0.6,0.4,0.6,0.6);Line:5:(0.5,0.4,0.6,0.4);")))
        return 0;
    else
        return -1;
}
