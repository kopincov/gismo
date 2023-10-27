/** @file gsCmdLineWithEnumSupport.h

    @brief Extends gsCmdLine with an interface to enums.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include<gsCore/gsForwardDeclarations.h>
#include<gsIO/gsCmdLine.h>

namespace gismo {

namespace util {
#if __cplusplus >= 201103
template <class T> struct maybe_enum { enum { value = std::is_enum<T>::value }; };
#else
template <class T> struct maybe_enum { enum { value = 1 }; };
#endif
}

template <class E> class gsCmdLineEnum;
class gsEnumHelper;


class GISMO_EXPORT gsCmdLineWithEnumSupport : public gsCmdLine
{

public:

    gsCmdLineWithEnumSupport(const std::string& message,
              const char delimiter = ' ',
              bool helpAndVersion = true) : gsCmdLine(message, delimiter, helpAndVersion) {}
    void getValues( int argc, char** argv );
    ~gsCmdLineWithEnumSupport();

public:

    /// @brief Register a enum option for the command line
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used
    /// @param name       Long form of the flag
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const variable, initialized
    ///                   with the default value. When \a getValues is
    ///                   invoked and the user provided a value, it is
    ///                   written to that variable.
    ///
    /// All possible options have to be registered using add subsequently:
    /// \code{.cpp}
    ///   namespace foo { enum value { bar = 0, baz = 1 }; }
    ///   
    ///   int main( int argc, char** argv )
    ///   {
    ///       foo::value f = foo::bar;
    ///
    ///       gsCmdLineWithEnumSupport cmd("My test file");
    ///
    ///       cmd.addEnum( flag, name, description, f )
    ///           .add(foo::bar, "bar", "the BAR-Option")
    ///           .add(foo::baz, "baz", "whatever baz means");
    ///
    ///       cmd.getValues(argc, argv);
    ///
    ///       return 0;
    ///   }
    /// \endcode
    ///
    template <class E>
    gsCmdLineEnum<E> addEnum( const std::string& flag, const std::string& name, const std::string& desc, E& e )
    {
        GISMO_STATIC_ASSERT( util::maybe_enum<E>::value&&sizeof(E)==sizeof(unsigned),
           "This function only works for enums that are based on int or unsigned int." );
        return gsCmdLineEnum<E>( addEnumImplementation( flag, name, desc, reinterpret_cast<unsigned&>(e) ) );
    }

private:
    gsEnumHelper* addEnumImplementation( const std::string& flag, const std::string& name, const std::string& desc, unsigned& e );
    std::vector< gsEnumHelper* > m_enumhelpers;

};


class GISMO_EXPORT gsEnumHelper
{
public:
    
    void add( const unsigned& val, const std::string & key, const std::string & desc );
    void writeKeyOfChosenOptTo( std::string & key );
    void writeDescOfChosenOptTo( std::string & desc );

private:
    gsEnumHelper( unsigned& val, const std::string & name ) : m_val(val), m_name(name), m_desc_ptr(NULL), m_key_ptr(NULL) {}
    void handle( bool abortOnError );

    friend class gsCmdLineWithEnumSupport;

    unsigned&                m_val;
    std::string              m_str;
    std::vector<std::string> m_keys;
    std::vector<std::string> m_descs;
    std::vector<unsigned>    m_vals;
    std::string              m_name;
    std::string             *m_desc_ptr;
    std::string             *m_key_ptr;
};


template <class E>
class gsCmdLineEnum
{
public:
    
    /// Register an element of the enum
    gsCmdLineEnum& add( const E& val, const std::string & key, const std::string & desc )
    { m_helper->add( reinterpret_cast<const unsigned&>(val), key, desc ); return *this; }
    
    /// The key of the chosen option will be written to key
    gsCmdLineEnum& writeKeyOfChosenOptTo( std::string & key )
    { m_helper->writeKeyOfChosenOptTo(key); return *this; }
    
    /// The description of the chosen option will be written to desc
    gsCmdLineEnum& writeDescOfChosenOptTo( std::string & key )
    { m_helper->writeDescOfChosenOptTo(key); return *this; }
    
private:
    friend class gsCmdLineWithEnumSupport;

    gsCmdLineEnum( gsEnumHelper* helper ) : m_helper(helper) {}
    gsCmdLineEnum( const gsCmdLineEnum& );
    gsCmdLineEnum& operator= ( const gsCmdLineEnum& );

    gsEnumHelper* m_helper;
};

} // namespace gismo

