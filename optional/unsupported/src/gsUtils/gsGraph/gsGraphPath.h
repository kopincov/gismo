#pragma once

#include <iostream>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

/** 
    A path connecting two vertices in a graph
*/

  
template<class T>
class gsGraphPath
{
// data members
protected:
    unsigned source; 
    unsigned target;
    T length;
    gsVector<unsigned>::Ptr previous;
  
public:
    /// Accessors to data members  
    unsigned getSource() const {return source;}
    unsigned getTarget() const {return target;}
    T getLength() const {return length;}
    gsVector<unsigned>::Ptr getPrevMat() const {return previous;}
  
public:
    /// Default empty constructor
    gsGraphPath() {}
  
    gsGraphPath(int const & s, int const & t, gsVector<unsigned>::Ptr prev, T const & l=0){source = s; target=t;previous = prev;length = l;}

    ~gsGraphPath() { } //destructor  
  
public:
  
    /// Prints the object as a string.  
    std::ostream &print(std::ostream &ost) const
    {
        std::vector<unsigned> path = computePath(previous,target);
        ost << "The shortest path from "<< source <<" to " << target << " : \n";	
        if ( ! path.empty() )
        {
            //for (unsigned i=0;i!=path.size();i++) ost << " "<< path[path.size()-1-i] << " ";
            for (unsigned i=0;i!=path.size();i++) ost << " "<< path[i] << " ";
        }

        ost << std::endl << "With length: " << length<<std::endl;
        return ost;
    }
     
    void setPrevious(gsVector<unsigned>::Ptr prev){previous=prev;};     

public:

    /// compute the path from *source* to *destination* given by *prev*
    std::vector<unsigned> computePath(gsVector<unsigned>::Ptr prev, unsigned destination) const;
    std::vector<unsigned> computePath() const {return computePath(previous,target);}

    void reArrangePath(std::vector<unsigned> & loop)
    {
        unsigned source2 = *(loop.end()-1);
        loop.insert(loop.begin(),source2);
        loop.erase(loop.end()-1);
    }

// Data members
private:


}; // class gsGraphPath


/// Define the functionality of the operator << acting on a class
template<class T>
std::ostream &operator<<(std::ostream &os, const gsGraphPath<T>& b)
{return b.print(os); }


//////////////////////////////////////////////////
//////////////////////////////////////////////////

template <class T>
std::vector<unsigned> gsGraphPath<T>::computePath(gsVector<unsigned>::Ptr prev, unsigned destination) const
{
    /// recover the path connecting *source* and *destination* using *prev* 
    std::vector<unsigned> path1;  
    path1.push_back(destination);
    unsigned  prepoint = (*prev)(destination);   
    path1.push_back(prepoint);
    while (prepoint != source)
    {      
        prepoint = (*prev)(prepoint); 
        path1.push_back(prepoint);
    }
     
    return path1;
}


} // namespace gismo
