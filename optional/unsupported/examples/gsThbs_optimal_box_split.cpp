/*
 *Created by: Gabor Kiss
Testing the optimal and aproximated rectangular partition of domain boundaries
*/

#include <iostream>
#include <fstream>

#include <gismo.h>
#include <gismo_dev.h>


using namespace gismo;
typedef gsHDomain<2,index_t>::point Point;
//struct Point{
//    real_t x;
//    real_t y;
//};

struct Ordering_FunctorX
{
    bool operator()(const Point& a, const Point& b)
    {
        if (a[0] == b[0])
            return a[1] < b[1];
        else
            return a[0] < b[0];
    }
};

struct Ordering_FunctorY
{
    bool operator()(const Point& a, const Point& b)
    {
        if (a[1] == b[1])
            return a[0] < b[0];
        else
            return a[1] < b[1];
    }
};

//Input:oriented connected component
std::vector<Point> getPoints(std::vector<std::vector<Point> > cc){
//all points in a connected component- including the points in the holes
    std::vector<Point> points;
    int position = 0;
    //go through all polylines in one conected component
    for(unsigned int i = 0; i < cc.size();i++)
    {
        //add space for new polygon vertices in points
        points.resize(points.size() + cc[i].size());
        //for every polyline take every second edge to get every point on the boundary
        for(unsigned int j = 0; j < cc[i].size();j++)
        {
            //line segment
            points[position][0] = cc[i][j][0];
            points[position][1] = cc[i][j][1];
            position++;
        }
    }
    //std::sort(points.begin(),points.end(),Ordering_Functor());
    return points;
}


index_t findCorner(std::vector< std::vector<index_t> > poly){
    index_t minx = 1000000;
    index_t maxy = 0;
    index_t result = 0;
    for(size_t i = 0; i < poly.size();i++)
    {
        if(poly[i][0]<minx)
        {
            result = i;
            minx = poly[i][0];
            maxy = poly[i][1];
        }
        if(poly[i][0]==minx)
        {
            if(maxy<poly[i][1])
            {
                result = i;
                maxy = poly[i][1];
            }
        }
    }
    return result;
}

//in:ordered corners of the polygon- including vertices of holes
//return all suspected diagonals-
std::vector< std::vector<Point> > identifyPossibleDiagonals(std::vector<Point> poly){
    std::vector<Point> temp_line;
    std::vector< std::vector<Point> > result;
    temp_line.resize(2);
    unsigned int dir;
    if(poly[0][0]==poly[1][0])
        dir = 0;
    else
        dir = 1;
    for(unsigned int i = 0; i < poly.size()-1; i++)
    {
        if((poly[i][dir] == poly[i+1][dir]) && (i%2==1))
        {
            temp_line[0][0] = poly[i][0];
            temp_line[0][1] = poly[i][1];
            temp_line[1][0] = poly[i+1][0];
            temp_line[1][1] = poly[i+1][1];
            result.push_back(temp_line);
        }
    }
    return result;
}

std::vector<std::vector<Point> > getOrientedComponent(std::vector<std::vector< std::vector<index_t> > >  cc){
    std::vector<std::vector<Point> > result;
    int position = 0;
    result.resize(cc.size());
    //1. find upper left corner
    for(unsigned int i = 0; i < cc.size();i++)
    {
        result[i].resize(cc[i].size());
        unsigned int  corner = findCorner(cc[i]);//its always a horizontal line
        result[i][0][0] = cc[i][corner][0];
        result[i][0][1] = cc[i][corner][1];
        position = 1;
        //2. chose orientation for the boundary (go left from the upper left corner) and oposit orientation for the holes(go down from upper left corner)
        if(i==0)// if the polygone is a domain boundary- oriented clockwise
        {
            if(cc[i][(corner+1)%cc[i].size()][0]==cc[i][corner][2])//if the next line segment has the same x coordinate as the end point of the corner line
            {
                result[i][position][0] = cc[i][corner][2];
                result[i][position][1] = cc[i][corner][3];
                position++;
                for(unsigned int j = 1; j < cc[i].size()-1;j++)// go through all line segments
                {
                    if( (cc[i][(corner+j)%cc[i].size()][0]!=result[i][position-1][0]) || (cc[i][(corner+j)%cc[i].size()][1]!=result[i][position-1][1]) )
                    {
                        result[i][position][0] = cc[i][(corner+j)%cc[i].size()][0];
                        result[i][position][1] = cc[i][(corner+j)%cc[i].size()][1];
                        position++;
                    }else
                    {
                        result[i][position][0] = cc[i][(corner+j)%cc[i].size()][2];
                        result[i][position][1] = cc[i][(corner+j)%cc[i].size()][3];
                        position++;
                    }
                }
            }else
            {
                result[i][position][0] = cc[i][corner][2];
                result[i][position][1] = cc[i][corner][3];
                position++;
                corner +=cc[i].size();//TODO TEST
                for(unsigned int j = 1; j < cc[i].size()-1;j++)// go through all line segments
                {
                    if( (cc[i][(corner-j)%cc[i].size()][0]!=result[i][position-1][0]) || (cc[i][(corner-j)%cc[i].size()][1]!=result[i][position-1][1]) )
                    {

                        result[i][position][0] = cc[i][(corner-j)%cc[i].size()][0];
                        result[i][position][1] = cc[i][(corner-j)%cc[i].size()][1];
                        position++;
                    }else
                    {
                        result[i][position][0] = cc[i][(corner-j)%cc[i].size()][2];
                        result[i][position][1] = cc[i][(corner-j)%cc[i].size()][3];
                        position++;
                    }
                }
            }
        }else{//holes in the domain should be oriented counteclockwise
            if(cc[i][(corner+1)%cc[i].size()][0]==cc[i][corner][2])
            {
                corner +=cc[i].size();
                for(unsigned int j = 1; j < cc[i].size();j++)
                {
                    if( (cc[i][(corner-j)%cc[i].size()][0]!=result[i][position-1][0]) || (cc[i][(corner-j)%cc[i].size()][1]!=result[i][position-1][1]) )
                    {
                        result[i][position][0] = cc[i][(corner-j)%cc[i].size()][0];
                        result[i][position][1] = cc[i][(corner-j)%cc[i].size()][1];
                        position++;
                    }else
                    {
                        result[i][position][0] = cc[i][(corner-j)%cc[i].size()][2];
                        result[i][position][1] = cc[i][(corner-j)%cc[i].size()][3];
                        position++;
                    }
                }
            }else
            {
                for(unsigned int j = 1; j < cc[i].size();j++)// go through all line segments
                {
                    if( (cc[i][(corner+j)%cc[i].size()][0]!=result[i][position-1][0]) || (cc[i][(corner+j)%cc[i].size()][1]!=result[i][position-1][1]) )
                    {
                        result[i][position][0] = cc[i][(corner+j)%cc[i].size()][0];
                        result[i][position][1] = cc[i][(corner+j)%cc[i].size()][1];
                        position++;
                    }else
                    {
                        result[i][position][0] = cc[i][(corner+j)%cc[i].size()][2];
                        result[i][position][1] = cc[i][(corner+j)%cc[i].size()][3];
                        position++;
                    }
                }
            }
        }
    }
    //3. return polygon
    return result;
}

//test intersection of two perpendicular lines
bool isIntersecting(std::vector<Point> diag, Point edge0, Point edge1)
{
    if( (diag[0][0] > edge0[0])&&(diag[0][0]< edge1[0]) && (diag[0][1] < edge0[1]) && (diag[1][1]> edge0[1])  )
    {//diag is x constant and edge y constant
        return true;
    }else
    {//diag is y constant and edge x constant
        if( (diag[0][1] > edge0[1]) && (diag[0][1]< edge1[1]) && (diag[0][0]<edge0[0]) && (diag[1][0]>edge1[0]) )
        return true;
    }
    return false;
}

//input: diag- vector of two Points
//       cc - connected component
bool testDiagonal(std::vector<Point> diag, std::vector<std::vector<Point> >cc){
    // test direction of diagonal
    bool diagDir = 0;//vertical diagonal
    unsigned int direction = 0;
    if(diag[0][0]!=diag[1][0])
    {
        diagDir = 1;//horizontal diagonal
    }
    //for all polylines test direction of the first edge
    for(unsigned int i = 0; i < cc.size();i++)
    {
        if(i==0)//boundary of the domain
        {
            if(diagDir)
                direction = 1;//horizontal
            else
                direction = 0;
        }else//holes in the domain
        {
            if(diagDir)
                direction = 0;//vertical
            else
                direction = 1;
        }
        for(unsigned int j = direction; j < cc[i].size()-1; j = j+2)
        {
            if(isIntersecting(diag,cc[i][j],cc[i][j+1])){
                return 0;
            }
        }
    }
    //test intersections with perpendicular edges
    //return true if no intersection found
    return 1;
}

int getDirection(Point v1, Point v2, Point v3)
{
    return (int(v2[0]) - int(v1[0]) ) * (int(v3[1]) - int(v2[1]) )  -
            (int(v3[0]) - int(v2[0]) ) * (int(v2[1]) - int(v1[1]) );
}

//return 1 if it is a diagonal
bool testAngles(std::vector<Point> diag, std::vector<std::vector<Point> >cc){
    //find the first point of the diagonal on cc
    for(unsigned int i = 0; i < cc.size();i++)
    {
        for(unsigned int j = 0; j < cc[i].size();j++)
        {
            if(cc[i][j]==diag[0]){
                int temp;//TODO change to TT
                if(j==0)
                    temp = getDirection(cc[i][cc[i].size()-1],cc[i][j],cc[i][(j+1)%cc[i].size()]);
                else
                    temp = getDirection(cc[i][j-1],cc[i][j],cc[i][(j+1)%cc[i].size()]);
                if( temp < 0 )
                    return 0;//not a diagonal
            }
        }
    }
    //compute angle between the two adjacent adges
    //if boundary and agle ==90- diagonal
    //if hole and angle ==270- diaognal
    //else delete

    return 1;
}


int main(int argc, char *argv[])
{
    //creating empty THB geometry
    gsTHBSpline<2>::uPtr THB;
    //reading file with THB spline geometry
    std::string filename = "basis2d/";
    //filename += "test_domain_boundary_00.xml";
    //filename += "test_domain_boundary_01.xml";
    //filename += "test_domain_boundary_02.xml";
    //filename += "test_domain_boundary_03.xml";
    //filename += "test_domain_boundary_04.xml";
    //filename += "test_domain_boundary_05.xml";
    //filename += "test_domain_boundary_06.xml";
    //filename += "test_domain_boundary_07.xml";
    filename += "test_domain_boundary_08.xml";
    gsFileData< real_t > data(filename);
    //creating the THB geometry object read from the file
    THB = data.getFirst< gsTHBSpline<2> >();
    gsInfo<<"THB data read correctly"<<"\n";

    //creating output variables for the getBsplinePatches_trimming(...) function see: THBsplineBasis
    gsVector<index_t> level;
    std::vector<std::vector<std::vector< std::vector<index_t> > > > con_comp;

    THB->basis().getConnectedComponents(con_comp,level);
    gsInfo<<"Levels are:\n"<<level<<"\n";


    for(unsigned int j = 0; j < con_comp[0].size();j++){
        if(j==0)
            gsInfo<<"Boundary of connected component 0 is: "<<"\n";
        else
            gsInfo<<"Hole "<< j-1<<" in connected component 0 is: "<<"\n";
        for (unsigned int i = 0; i < con_comp[0][j].size(); i++){//level 0, first curve
            gsInfo<<"["<< con_comp[0][j][i][0]<<" , "<<con_comp[0][j][i][1]<<"],["<<con_comp[0][j][i][2]<<" , "<<con_comp[0][j][i][3]<<"]"<<"\n";
        }
    }


    //getDiagonals(con_comp[0], diag,1);
    //std::reverse(con_comp[0][1].begin(), con_comp[0][1].end());
    gsInfo<<"Oriented line is: "<<"\n";
    std::vector<std::vector<Point> >oriented_line = getOrientedComponent(con_comp[0]);
    for(unsigned int i =0; i < oriented_line[0].size(); i++){
        gsInfo<<"["<<oriented_line[0][i][0]<<" , "<<oriented_line[0][i][1]<<"]"<<"\n";
    }
    gsInfo<<"Oriented hole"<<"\n";
    for(unsigned int i =0; i < oriented_line[1].size(); i++){
        gsInfo<<"["<<oriented_line[1][i][0]<<" , "<<oriented_line[1][i][1]<<"]"<<"\n";
    }


    std::vector<Point> polyPoints;
    polyPoints = getPoints(oriented_line);
    std::sort(polyPoints.begin(),polyPoints.end(),Ordering_FunctorX());//sort by X
    gsInfo<<"\n"<<"Sorted points are: "<<"\n";
    for(unsigned int i =0; i < polyPoints.size(); i++){
        gsInfo<<"["<< polyPoints[i][0]<<" , "<<polyPoints[i][1]<<"]"<<"\n";
    }

    std::vector< std::vector<Point> > diag;
    diag = identifyPossibleDiagonals(polyPoints);
    gsInfo<<"\n"<<"Possible diagonals are: "<<"\n";
    for(unsigned int i =0; i < diag.size(); i++){
        gsInfo<<"["<< diag[i][0][0]<<" , "<<diag[i][0][1]<<"], "<<"["<< diag[i][1][0]<<" , "<<diag[i][1][1]<<"]"<<"\n";
    }
    for (unsigned int i = 0; i<diag.size();i++){
        if(testDiagonal(diag[i],oriented_line)){
            gsInfo<<"hura 1"<<"\n";
        }else{
            gsInfo<<"buuuu 1"<<"\n";
            diag.erase(diag.begin()+i);
            i--;
        }
    }
    for (unsigned int i = 0; i<diag.size();i++){
        if(testAngles(diag[i],oriented_line)){
            gsInfo<<"hura 2"<<"\n";
        }else{
            gsInfo<<"buuuu 2"<<"\n";
            diag.erase(diag.begin()+i);
            i--;
        }
    }

    return 0;
}




//backup for some old functions
//geting a conected component
//std::vector<Point> getSortedPoints(std::vector< std::vector< std::vector< unsigned int > > > cc){
//all points in a cnected component- including the points in the holes
//    std::vector<Point> points;
//    int position = 0;
//    //go through all polylines in one conected component
//    for(unsigned int i = 0; i < cc.size();i++){
//        //add space for new polygon vertices in points
//        points.resize(points.size() + cc[i].size());
//        //for every polyline take every second edge to get every point on the boundary
//        for(unsigned int j = 0; j < cc[i].size();j = j+2)
//        {
//            //line segment
//            points[position][0] = cc[i][j][0];
//            points[position][1] = cc[i][j][1];
//            points[position+1][0] = cc[i][j][2];
//            points[position+1][1] = cc[i][j][3];
//            position =  position+2;
//        }
//    }
//    std::sort(points.begin(),points.end(),Ordering_Functor());
//    return points;
//}
//in:ordered corners of the polygon- including vertices of holes
//return all suspected diagonals-
//std::vector< std::vector<Point> > identifyPossibleDiagonals(std::vector<Point> poly){
//    std::vector<Point> temp_line;
//    std::vector< std::vector<Point> > result;
//    temp_line.resize(2);
//    for(unsigned int i = 0; i < poly.size()-1; i++){
//        if((poly[i][0] == poly[i+1][0]) && (i%2==1)){
//            temp_line[0][0] = poly[i][0];
//            temp_line[0][1] = poly[i][1];
//            temp_line[1][0] = poly[i+1][0];
//            temp_line[1][1] = poly[i+1][1];
//            result.push_back(temp_line);
//        }
//    }
//    return result;
//}
//test intersection of two perpendicular lines
//bool isIntersecting(std::vector<Point> diag, std::vector<index_t> edge)
//{
//    if( (diag[0][0] >= edge[0])&&(diag[0][0]<= edge[2]) && (diag[0][1] <= edge[1]) && (diag[1][1]>= edge[1])  )
//    {//diag is x constant and edge y constant
//        return true;
//    }else
//    {//diag is y constant and edge x constant
//        if( (diag[0][1] >= edge[1]) && (diag[0][1]<= edge[3]) && (diag[0][0]<=edge[0]) && (diag[1][0]>=edge[2]) )
//        return true;
//    }
//    return false;
//}


//direction 0-> xconstant
//direction 1-> yconstant
//cc->coencted component
//temp_diag- possible diagonals
//std::vector< std::vector<Point> > getDiagonals(std::vector<std::vector< std::vector< unsigned int > > > cc,
//                                               std::vector< std::vector<Point> > temp_diag, int direction)
//{
//    std::vector< std::vector<Point> > real_diag;
//    unsigned int k;
//    for(unsigned int ii = 0; ii < temp_diag.size();ii++)
//    {
//        unsigned int counter = 0;
//        for(unsigned int i = 0; i< cc.size();i++)
//        {
//            if(cc[i][0][direction] == cc[i][0][direction+2]){//if the first edge is constant in "direction"
//                k = 0;
//            }else
//            {
//                k = 1;
//            }
//            for(unsigned int j = k; j < cc[i].size()-1; j = j+2)
//            {
//                if(isIntersecting(temp_diag[ii],cc[i][j]))
//                {
//                    counter++;//count how may intersection there are with perpendicular edges
//                }
//            }
//        }
//        if(counter==2)//just the two endpoints can intersect with a perpencicular edge in case of a real diagonal
//        {

//        }
//    }
//    return real_diag;
//}



