/** @file gsJITCompiler_test.cpp

    @brief Provides test of the JIT-compiler

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moeller
*/

#include <gismo.h>
#include <gismo_dev.h>
#include <gsCore/gsJITCompiler.h> // to be included later

using namespace gismo;

// Main program
int main()
{    
    // Create compiler configuration
    gsJITCompilerConfig CC;

    // *** 1. Test C-code kernel ***
    {
        // Load C compiler configuration
        CC.load(GISMO_CONFIG_DIR "jit.xml", gsJITLang::C);
    
        // Create JIT compiler
        gsJITCompiler compiler(CC);
        
        // Write demo compute kernel
        compiler << "#include \"stdlib.h\"\n"
                 << "int foo1( int i, int j, int k )\n"
                 << "{\n"
                 << "  return i+j+k;\n"
                 << "}\n"
                 << util::type<real_t>::name() << " foo2( "
                 << util::type<real_t>::name() << " i, "
                 << util::type<real_t>::name() << " j, "
                 << util::type<real_t>::name() << " k )\n"
                 << "{\n"
                 << "  return i+j+k;\n"
                 << "}\n";
        
        // Output compute kernel
        gsInfo << "Compute kernel code:\n" << compiler << "\n";
        
        typedef int    (dl_func1)(int,int,int);
        typedef real_t (dl_func2)(real_t,real_t,real_t);
        
        // Create dynamic library object
        try {
            gsDynamicLibrary dl = compiler.build(true); 
            dl_func1 * f1 = dl.getSymbol<dl_func1>("foo1");
            dl_func2 * f2 = dl.getSymbol<dl_func2>("foo2");
            
            // Call dynamically load
            std::cout << f1(1,2,3) << ", "
                      << f2(1.0,2.0,3.0) << std::endl;
        } catch (const std::exception& e) {
            gsInfo << e.what() <<"\n";
            return EXIT_FAILURE;
        }
    }
    
    // *** 2. Test C++-code kernel ***
    {
        // Load C++ compiler configuration
        CC.load(GISMO_CONFIG_DIR "jit.xml", gsJITLang::CXX);
    
        // Create JIT compiler
        gsJITCompiler compiler(CC);
        
        // Write demo compute kernel
        compiler << "#include \"iostream\"\n"
                 << "EXPORT int foo1( int i, int j, int k )\n"
                 << "{\n"
                 << "  return i+j+k;\n"
                 << "}\n"
                 << "EXPORT "
                 << util::type<real_t>::name() << " foo2( "
                 << util::type<real_t>::name() << " i, "
                 << util::type<real_t>::name() << " j, "
                 << util::type<real_t>::name() << " k )\n"
                 << "{\n"
                 << "  return i+j+k;\n"
                 << "}\n";
        
        // Output compute kernel
        gsInfo << "Compute kernel code:\n" << compiler << "\n";
        
        typedef int    (dl_func1)(int,int,int);
        typedef real_t (dl_func2)(real_t,real_t,real_t);
        
        // Create dynamic library object
        try {
            gsDynamicLibrary dl = compiler.build(true); 
            dl_func1 * f1 = dl.getSymbol<dl_func1>("foo1");
            dl_func2 * f2 = dl.getSymbol<dl_func2>("foo2");
            
            // Call dynamically load
            std::cout << f1(1,2,3) << ", "
                      << f2(1.0,2.0,3.0) << std::endl;
        } catch (const std::exception& e) {
            gsInfo << e.what() <<"\n";
            return EXIT_FAILURE;
        }
    }

    // *** 3. Test Fortran-code kernel ***
    {
        // Load Fortran compiler configuration
        CC.load(GISMO_CONFIG_DIR "jit.xml", gsJITLang::Fortran);
        
        // Create JIT compiler
        gsJITCompiler compiler(CC);
        
        // Write demo compute kernel
        compiler << "INTEGER FUNCTION foo1(i, j, k)\n"
                 << "IMPLICIT NONE\n"
                 << "INTEGER, INTENT(IN) :: i,j,k\n"
                 << "foo1 = i+j+k\n"
                 << "END FUNCTION foo1\n";
        
        // Output compute kernel
        gsInfo << "Compute kernel code:\n" << compiler << "\n";
        
        typedef int (dl_func1)(int*,int*,int*);
        
        // Create dynamic library object
        try {
            gsDynamicLibrary dl = compiler.build(true); 
            dl_func1 * f1 = dl.getSymbol<dl_func1>("foo1_");
            
            // Call dynamically load
            int i = 1;
            int j = 2;
            int k = 3;
            std::cout << f1(&i, &j, &k) << std::endl;
        } catch (const std::exception& e) {
            gsInfo << e.what() <<"\n";
            return EXIT_FAILURE;
        }
    }
    
    return EXIT_SUCCESS;
}
