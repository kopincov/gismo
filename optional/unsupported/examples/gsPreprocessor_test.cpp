/*
* gsPreprocessor_test.cpp created on 11.05.2017
*
* Author: J. Vogl
*
* This file is part of the G+SMO library
*/

//#include <gismo.h>
//#include <gsCore/gsFunctionSet.h>
#include <iostream>

// Helper macros for counting arguments, works till highest number in __RSEQ_N
// Call with __VA_NARG__(__VA_ARGS__)
#define EXPAND(x) x

#define __VA_NARG__(...) __VA_NARG_(_0, ## __VA_ARGS__, __RSEQ_N())
#define __VA_NARG_(...) EXPAND(__VA_ARG_N(__VA_ARGS__))
#define __VA_ARG_N(_1, _2, _3, _4, N,...) N
#define __RSEQ_N() 3, 2, 1, 0

#define PP_EXPAND(x) x
#define PP_ARG_N(_1,_2,_3,_4,N,...) N
#define PP_RSEQ_N() 4,3,2,1,0
#define PP_NARG_(...) PP_EXPAND(PP_ARG_N(__VA_ARGS__))
#define PP_COMMASEQ_N()  1,1,1,0,0
#define PP_COMMA(...) ,
#define PP_HASCOMMA(...) PP_NARG_(__VA_ARGS__,PP_COMMASEQ_N())
#define PP_NARG(...) PP_NARG_HELPER1(PP_HASCOMMA(__VA_ARGS__),PP_HASCOMMA(PP_COMMA __VA_ARGS__ ()),PP_NARG_(__VA_ARGS__, PP_RSEQ_N()))
#define PP_NARG_HELPER1(a,b,N) PP_NARG_HELPER2(a, b, N)
#define PP_NARG_HELPER2(a,b,N) PP_NARG_HELPER2 ## a(b,N)
#define PP_NARG_HELPER20(b,N) !b
#define PP_NARG_HELPER21(b,N) N

//using namespace gismo;

int main (int argc, char* argv[])
{
    int result1 = 0;
    std::cout << "__VA_NARG__()       " << __VA_NARG__() << std::endl
              << "__VA_NARG__(x)      " << __VA_NARG__(x) << std::endl
              << "__VA_NARG__(x,x)    " << __VA_NARG__(x,x) << std::endl
              << "__VA_NARG__(x,x,x)  " << __VA_NARG__(x,x,x) << std::endl;
    result1 |= (__VA_NARG__() != 0);
    result1 |= (__VA_NARG__(x) != 1);
    result1 |= (__VA_NARG__(x,x) != 2);
    result1 |= (__VA_NARG__(x,x,x) != 3);
    std::cout << "result __VA_NARG__: " << result1 << std::endl << std::endl;

    std::cout << "PP_NARG()           " << PP_NARG() << std::endl
              << "PP_NARG(x)          " << PP_NARG(x) << std::endl
              << "PP_NARG(x,x)        " << PP_NARG(x,x) << std::endl
              << "PP_NARG(x,x,x)      " << PP_NARG(x,x,x) << std::endl
              << "PP_NARG(x,x,x,x)    " << PP_NARG(x,x,x,x) << std::endl;
    int result2 = 0;
    result2 |= ((PP_NARG()) != 0);
    result2 |= ((PP_NARG(x)) != 1);
    result2 |= ((PP_NARG(x,x)) != 2);
    result2 |= ((PP_NARG(x,x,x)) != 3);
    std::cout << "result PP_NARG:     " << result2 << std::endl;

    if (!result1 || !result2) // OR is correct. __VA_NARG__() can fail with some compilers.
        return 0;
    return (result1 || result2);
}
