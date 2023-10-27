/** @file gsMemory2Test.cpp
    @brief Provides utility function related to memory management.
    This file is part of the G+Smo library. 
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/


#include <gismo.h>
#include <gsCore/gsMemory2.h>

using namespace gismo;

struct A {
    A() { std::cout << "A()\n"; }
    virtual ~A() { std::cout << "~A()\n"; }
};

struct B : public A {
    B() { std::cout << "B()\n"; }
    virtual ~B() { std::cout << "~B()\n"; }
};



int main() {

    {
        std::cout << "First testcase.\n";
    
        std::cout << "memory::shared_ptr<A> ptr( new B );\n";
        memory::shared_ptr<A> ptr( new B );
        
        std::cout << "memory::shared_ptr<B> ptr2 = memory::static_pointer_cast<B>(ptr);\n";
        memory::shared_ptr<B> ptr2 = memory::static_pointer_cast<B>(ptr);

        std::cout << "ptr.reset();\n";
        ptr.reset();
        
        std::cout << "Go out of scope.\n";
    
    }
    std::cout << "Done.\n";
    
    {
        std::cout << "Second testcase.\n";
    
        std::cout << "memory::shared_ptr<A> ptr( new B );\n";
        memory::shared_ptr<A> ptr( new B );
        
        std::cout << "memory::shared_ptr<B> ptr2 = memory::static_pointer_cast<B>(ptr);\n";
        memory::shared_ptr<B> ptr2 = memory::static_pointer_cast<B>(ptr);

        std::cout << "ptr2.reset();\n";
        ptr2.reset();
        
        std::cout << "Go out of scope.\n";
    
    }
    std::cout << "Done.\n";
    
    {
        std::cout << "Third testcase.\n";
    
        std::cout << "memory::shared_ptr<A> ptr( new B );\n";
        memory::shared_ptr<A> ptr( new B );
        
        std::cout << "memory::shared_ptr<B> ptr2 = memory::static_pointer_cast<B>(ptr);\n";
        memory::shared_ptr<B> ptr2 = memory::static_pointer_cast<B>(ptr);

        std::cout << "ptr.reset();\n";
        ptr.reset();
        
        std::cout << "ptr2.reset();\n";
        ptr2.reset();
        
        std::cout << "Go out of scope.\n";
    
    }
    std::cout << "Done.\n";
    

    {
        std::cout << "Fourth testcase.\n";
    
        std::cout << "memory::shared_ptr<A> ptr( new B );\n";
        memory::shared_ptr<A> ptr( new B );
        
        std::cout << "memory::shared_ptr<B> ptr2 = memory::static_pointer_cast<B>(ptr);\n";
        memory::shared_ptr<B> ptr2 = memory::static_pointer_cast<B>(ptr);
        
        std::cout << "Go out of scope.\n";
    
    }
    std::cout << "Done.\n";

    return 0;
}

