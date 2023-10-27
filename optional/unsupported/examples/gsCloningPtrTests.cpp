/** @file gsCloningPtrTest.cpp

    @brief 

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/

#include <gismo.h>
#include <gsUtils/gsCloningPtr.h>

using namespace gismo;


class testSuperBase
{
public:
    virtual ~testSuperBase() {}
    virtual testSuperBase* clone() const = 0;
    virtual int get() = 0;
};

class testBase : public testSuperBase
{
public:
    virtual ~testBase() {}
    virtual testBase* clone() const = 0;
    virtual int get() = 0;
};

class test : public testBase
{
    int m_val;
    
public :
    
    test(int val=0)     : m_val(val)           { gsInfo << "test(" << val << ")      " << this << " NULL -> " << m_val << "\n";                                                   }
    test(const test& o) : m_val(o.m_val)       { gsInfo << "test(const test&o)       " << this << " NULL -> " << m_val << "\n";                                                   }
    test& operator=(const test& o)             { gsInfo << "operator=(const test& o) " << this << " " << m_val << " -> " << o.m_val << "\n";  m_val = o.m_val; return *this;      }
#if EIGEN_HAS_RVALUE_REFERENCES
    test(test&& o)      : m_val(give(o.m_val)) { gsInfo << "test(test&&o)            " << this << " NULL -> " << m_val << "\n";                                                   }
    test& operator=(test&& o)                  { gsInfo << "operator=(test&& o)      " << this << " " << m_val << " -> " << o.m_val << "\n"; return *this;                        }
#endif
    ~test()                                    { gsInfo << "~test()                  " << this << " " << m_val << " -> NULL\n";                                                   }
    virtual test* clone() const                { gsInfo << "clone()                  " << this << " " << m_val << "\n"; return new test(*this);                                   }
    virtual int get()                          { return m_val;                                                                                                                    }
 
};


#define TEST_ASSIGN( type1, type2 )                                                 \
{                                                                                   \
    gsInfo << "\n\n**************************************************\n";           \
    gsInfo << "  a -> " << #type1 << "   b -> " << #type2 << "\n";                  \
    gsInfo << "\nAllocating...\n";                                                  \
    type1  a( new test(1) );                                                        \
    type2  b;                                                                       \
    gsInfo << "\nTesting b = a...\n";                                               \
    b = a;                                                                          \
    gsInfo << "\n a.get()  = " << a.get()  << " b.get()  = " << b.get()  << "\n";   \
    gsInfo <<   " a->get() = " << a->get() << " b->get() = " << b->get() << "\n";   \
    gsInfo << "\nCleaning up...\n";                                                 \
}

#define TEST_ASSIGN_CLONE( type1, type2 )                                           \
{                                                                                   \
    gsInfo << "\n\n**************************************************\n";           \
    gsInfo << "  a -> " << #type1 << "   b -> " << #type2 << "\n";                  \
    gsInfo << "\nAllocating...\n";                                                  \
    type1  a( new test(1) );                                                        \
    type2  b;                                                                       \
    gsInfo << "\nTesting b = memory::clone(a)...\n";                                \
    b = memory::clone(a);                                                           \
    gsInfo << "\n a.get()  = " << a.get()  << " b.get()  = " << b.get()  << "\n";   \
    gsInfo <<   " a->get() = " << a->get() << " b->get() = " << b->get() << "\n";   \
    gsInfo << "\nCleaning up...\n";                                                 \
}

#define TEST_ASSIGN_GIVE( type1, type2 )                                            \
{                                                                                   \
    gsInfo << "\n\n**************************************************\n";           \
    gsInfo << "  a -> " << #type1 << "   b -> " << #type2 << "\n";                  \
    gsInfo << "\nAllocating...\n";                                                  \
    type1  a( new test(1) );                                                        \
    type2  b;                                                                       \
    gsInfo << "\nTesting b = give(a)...\n";                                         \
    b = give(a);                                                                    \
    gsInfo << "\n a.get()  = " << a.get()  << " b.get()  = " << b.get()  << "\n";   \
    gsInfo <<   " a->get() = NON-DEF            b->get() = " << b->get() << "\n";   \
    gsInfo << "\nCleaning up...\n";                                                 \
}

int main()
{
      TEST_ASSIGN( memory::cloning_ptr<testBase> , memory::cloning_ptr<testSuperBase>  )
    //TEST_ASSIGN( memory::unique_ptr<testBase>, memory::cloning_ptr<testSuperBase>  )
      gsInfo << "\n\n**************************************************\n"
        << "Is not possible to ASSIGN a UNIQUE pointer to a CLONE pointer variable because both pointers have UNIQUE ownership and we do not want an implicit CLONE.";
    //TEST_ASSIGN( memory::shared_ptr<testBase>, memory::cloning_ptr<testSuperBase>  )
      gsInfo << "\n\n**************************************************\n"
        << "Is not possible to ASSIGN a SHARED pointer to a CLONE pointer variable because CLONE pointers have UNIQUE ownership and we do not want an implicit CLONE.";
    
    //TEST_ASSIGN( memory::cloning_ptr<testBase>,  memory::unique_ptr<testSuperBase> )
      gsInfo << "\n\n**************************************************\n"
        << "Is not possible to ASSIGN a CLONE pointer to a UNIQUE pointer variable because both pointers have UNIQUE ownership and we do not want an implicit CLONE.";
    //TEST_ASSIGN( memory::unique_ptr<testBase>, memory::unique_ptr<testSuperBase> )
      gsInfo << "\n\n**************************************************\n"
        << "Is not possible to ASSIGN a UNIQUE pointer to a UNIQUE pointer variable (already due to std).";
    //TEST_ASSIGN( memory::shared_ptr<testBase>, memory::unique_ptr<testSuperBase> )
      gsInfo << "\n\n**************************************************\n"
        << "Is not possible to ASSIGN a SHARED pointer to a UNIQUE pointer variable (already due to std).";

    //TEST_ASSIGN( memory::cloning_ptr<testBase>,  memory::shared_ptr<testSuperBase> )
      gsInfo << "\n\n**************************************************\n"
        << "Is not possible to ASSIGN a CLONE pointer to a SHARED pointer variable because CLONE pointers have UNIQUE ownership and we do not want an implicit CLONE.";
    //TEST_ASSIGN( memory::unique_ptr<testBase>, memory::shared_ptr<testSuperBase> )
      gsInfo << "\n\n**************************************************\n"
        << "Is not possible to ASSIGN a UNIQUE pointer to a SHARED pointer variable (already due to std).";
      TEST_ASSIGN( memory::shared_ptr<testBase>, memory::shared_ptr<testSuperBase> )

      
      
      TEST_ASSIGN_CLONE( memory::cloning_ptr<testBase>, memory::cloning_ptr<testSuperBase>  )
      TEST_ASSIGN_CLONE( memory::unique_ptr<testBase>,  memory::cloning_ptr<testSuperBase>  )
      TEST_ASSIGN_CLONE( memory::shared_ptr<testBase>,  memory::cloning_ptr<testSuperBase>  )
    
      TEST_ASSIGN_CLONE( memory::cloning_ptr<testBase>, memory::unique_ptr<testSuperBase> )
      TEST_ASSIGN_CLONE( memory::unique_ptr<testBase>,  memory::unique_ptr<testSuperBase> )
      TEST_ASSIGN_CLONE( memory::shared_ptr<testBase>,  memory::unique_ptr<testSuperBase> )

      TEST_ASSIGN_CLONE( memory::cloning_ptr<testBase>, memory::shared_ptr<testSuperBase> )
      TEST_ASSIGN_CLONE( memory::unique_ptr<testBase>,  memory::shared_ptr<testSuperBase> )
      TEST_ASSIGN_CLONE( memory::shared_ptr<testBase>,  memory::shared_ptr<testSuperBase> )
      
      
      
      TEST_ASSIGN_GIVE( memory::cloning_ptr<testBase>, memory::cloning_ptr<testSuperBase>  )
      TEST_ASSIGN_GIVE( memory::unique_ptr<testBase>,  memory::cloning_ptr<testSuperBase>  )
    //TEST_ASSIGN_GIVE( memory::shared_ptr<testBase>,  memory::cloning_ptr<testSuperBase>  )
      gsInfo << "\n\n**************************************************\n"
        << "Note that it is not possible to GIVE a SHARED pointer to a CLONE pointer variable because the other SHARED owners might disagree with that!";
      
      
      TEST_ASSIGN_GIVE( memory::cloning_ptr<testBase>, memory::unique_ptr<testSuperBase> )
      TEST_ASSIGN_GIVE( memory::unique_ptr<testBase>,  memory::unique_ptr<testSuperBase> )
    //TEST_ASSIGN_GIVE( memory::shared_ptr<testBase>,  memory::unique_ptr<testSuperBase> )
      gsInfo << "\n\n**************************************************\n"
        << "Note that it is not possible to GIVE a SHARED pointer to a UNIQUE pointer variable because the other SHARED owners might disagree with that!";

      TEST_ASSIGN_GIVE( memory::cloning_ptr<testBase>, memory::shared_ptr<testSuperBase> )
      TEST_ASSIGN_GIVE( memory::unique_ptr<testBase>,  memory::shared_ptr<testSuperBase> )
      TEST_ASSIGN_GIVE( memory::shared_ptr<testBase>,  memory::shared_ptr<testSuperBase> )
      
      
      return 0;
      
}
