#pragma once

#include <iostream>

namespace gismo
{

  /** 
      Options for volume segmentation of a solid given by its trimmed surfaces (gsSolid). 
  */

  
template<class T>
class gsGraphOptions
{
  
// Data members
public:
  int InfRep; // approach for REPresent INFinity; values: -1 (Inf=-1);
  
  
public:

  /// Default empty constructor
  gsGraphOptions() 
  { 
      InfRep = -1;
  };  

  virtual ~gsGraphOptions() {}; //destructor
  
  void updateBase(std::string field,int value)
  {
    if (field=="InfRep")
      InfRep = value;
//     else
//     {
//       std::cout << std::endl << "The intered name does not match any of the members of the class gsGraphOptions" << std::endl;
//       exit(1);
//     };
      
  };
  
  /// Prints the object as a string
  virtual std::ostream &print(std::ostream &os) const 
  {
    os<<"gsGraphOptions:\n";
    os<<"InfRep: "<< InfRep << std::endl;
    return os;       
  };  
  


}; 


//////////////////////////////////////////////////
//////////////////////////////////////////////////
/// Define the functionality of the operator << acting on a class
template<class T>
std::ostream &operator<<(std::ostream &os, const gsGraphOptions<T>& b)
{return b.print(os); };




}; // namespace gismo
