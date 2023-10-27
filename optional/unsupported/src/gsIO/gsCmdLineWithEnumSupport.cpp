/** @file gsCmdLineWithEnumSupport.cpp

    @brief Extends gsCmdLine with an interface to enums.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gsIO/gsCmdLineWithEnumSupport.h>

namespace gismo {

void gsEnumHelper::writeKeyOfChosenOptTo( std::string & key )
{
    GISMO_ASSERT ( !m_key_ptr, "Cannot invoke gsEnumHelper::writeKeyOfChosenOptTo more than once." );
    m_key_ptr = &key;
}

void gsEnumHelper::writeDescOfChosenOptTo( std::string & desc )
{
    GISMO_ASSERT ( !m_desc_ptr, "Cannot invoke gsEnumHelper::writeDescOfChosenOptTo more than once." );
    m_desc_ptr = &desc;
}

void gsEnumHelper::add( const unsigned& val, const std::string & key, const std::string & desc )
{
    GISMO_ASSERT ( !key.empty(), "gsEnumHelper::add: The given key must not be empty." );
    m_vals.push_back(val);
    m_keys.push_back(key);
    m_descs.push_back(desc);
}

void gsEnumHelper::handle( bool abortOnError )
{
    const index_t sz = m_vals.size();
    GISMO_ASSERT ( sz>0, "gsCmdLine: No values have been registered for the enum \"" << m_name << "\"." );
    // The parameter was not set, so we use defualt
    if ( m_str.empty() )
    {
        for (index_t i=0; i<sz; ++i)
            if ( m_vals[i] == m_val )
            {
                m_str = m_keys[i];
                if ( m_key_ptr )
                    *m_key_ptr = m_keys[i];
                if ( m_desc_ptr )
                    *m_desc_ptr = m_descs[i];
            }
        return;
    }
    // If the parameter was set, we search for it...
    for (index_t i=0; i<sz; ++i)
        if ( m_keys[i] == m_str )
        {
            m_val = m_vals[i];
            if ( m_key_ptr )
                *m_key_ptr = m_keys[i];
            if ( m_desc_ptr )
                *m_desc_ptr = m_descs[i];
            return;
        }
    // If the description has been provided, we might also accept it...
    for (index_t i=0; i<sz; ++i)
        if ( !m_descs[i].empty() && m_descs[i] == m_str )
        {
            m_val = m_vals[i];
            m_str = m_keys[i];
            if ( m_key_ptr )
                *m_key_ptr = m_keys[i];
            if ( m_desc_ptr )
                *m_desc_ptr = m_descs[i];
            return;
        }
    // If the given option is unknown...
    if (abortOnError)
    {
        gsInfo << " ERROR: The option \"" << m_name << "\" cannot take the value \"" << m_str << "\".\n\n USAGE:\n  Allowed values are:\n\n";
        std::string* default_opt = NULL;
        for (index_t i=0; i<sz; ++i)
        {
            gsInfo << "    * "
                   << std::setw( 9)<<std::left<< ("\""+m_keys[i]+"\"")
                   << " " << m_descs[i] << "\n";
            if ( m_vals[i] == m_val ) default_opt = &(m_keys[i]);
        }
        gsInfo << "\n";
        if (default_opt)
            gsInfo << "  The value \"" << *default_opt << "\" is the default.\n\n";
        std::exit(-1);
    }
    else
        throw std::runtime_error( "An undefined value has been chosen for an enum option." );
}

void gsCmdLineWithEnumSupport::getValues( int argc, char** argv )
{
    gsCmdLine::getValues( argc, argv );

    index_t sz = m_enumhelpers.size();
    bool abortOnError = gsCmdLine::getExceptionHandling();
    for (index_t i = 0; i<sz; ++i)
        m_enumhelpers[i]->handle( abortOnError );
}

gsEnumHelper* gsCmdLineWithEnumSupport::addEnumImplementation( const std::string& flag, const std::string& name, const std::string& desc, unsigned& e )
{
    index_t sz = m_enumhelpers.size();
    m_enumhelpers.push_back( new gsEnumHelper(e,name) );
    // internally handle this as a string option
    gsCmdLine::addString( flag, name, desc, m_enumhelpers[sz]->m_str );
    return m_enumhelpers[sz];
}

gsCmdLineWithEnumSupport::~gsCmdLineWithEnumSupport()
{
    freeAll(m_enumhelpers);
}


} // namespace gismo
