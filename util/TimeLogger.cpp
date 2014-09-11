/*
 * The MIT License
 *
 * Copyright (c) 2014 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/**
 *  \file   TimeLogger.cpp
 *  \date   Sep 8, 2014
 *  \author "James C. Sutherland"
 */

//-----------------------------------------------------------------------------

#include <util/TimeLogger.h>

#include <boost/foreach.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

//-----------------------------------------------------------------------------

TimeLogger::TimeLogger( const std::string logFileName,
                        const Format format )
: logFileName_( logFileName ),
  haveFileIO_( true ),
  format_( format )
{
  std::ifstream infile( logFileName.c_str() );
  if( infile.good() ){
    // load the existing file from disk rather than simply appending it.
    // This will ensure that the tree structure is properly maintained.
    try{
      switch( format_ ){
        case JSON: boost::property_tree::json_parser::read_json( infile, pt_ ); break;
        case XML : boost::property_tree::xml_parser ::read_xml ( infile, pt_ ); break;
      }
    }
    catch( std::exception& err ){
      logFileName_ = logFileName + ".new";
      std::cout << "\n\nNOTE: there was an error reading the log file " << logFileName
          << "Details follow:\n"
          << err.what()
          << "\nTo avoid overwriting the log file, a different log file named \n\t"
          << logFileName_
          << "\nwill be used\n\n";
    }
  }
  // start the timer that measures the lifetime of the TimeLogger.
  totalTime_.start();
}

TimeLogger::TimeLogger( const Format format )
: haveFileIO_( false ),
  format_( format )
{}

//-----------------------------------------------------------------------------

TimeLogger::~TimeLogger()
{
  totalTime_.stop();
  write_entries();
}

//-----------------------------------------------------------------------------

void TimeLogger::start( const std::string& label ){
  Entries::iterator it = entries_.find( label );
  if( it == entries_.end() ){
    it = entries_.insert( entries_.begin(), make_pair(label, Timer()) );
  }
  it->second.start();
}

//-----------------------------------------------------------------------------

void TimeLogger::reset( const std::string& label ){
  Entries::iterator it = entries_.find( label );
  if( it == entries_.end() ){
    it = entries_.insert( entries_.begin(), make_pair(label, Timer()) );
  }
  it->second.reset();
}

//-----------------------------------------------------------------------------

double TimeLogger::stop( const std::string& label )
{
# ifndef NDEBUG
  assert( entries_.find( label ) != entries_.end() );
# endif
  return entries_[label].stop();
}

//-----------------------------------------------------------------------------

void TimeLogger::add_entry( const std::string& label,
                            const Timer& timer )
{
  entries_[label] = timer;
}

//-----------------------------------------------------------------------------


const Timer& TimeLogger::timer( const std::string& label ) const
{
  Entries::const_iterator it = entries_.find( label );
  assert( it != entries_.end() );
  return it->second;
}

//-----------------------------------------------------------------------------

double TimeLogger::total_time()
{
  totalTime_.stop();
  const double t = totalTime_.elapsed_time();
  totalTime_.start();
  return t;
}

//-----------------------------------------------------------------------------

void TimeLogger::write_entries()
{
  const boost::posix_time::ptime time = boost::posix_time::second_clock::local_time();
  const std::string timeStamp = boost::posix_time::to_simple_string(time);

  boost::property_tree::ptree pt;

  BOOST_FOREACH( const Entries::value_type& vt, entries_ ){
    const Timer& t = vt.second;
    const std::string& label = vt.first;
    pt.put( label, t.elapsed_time() );
  }
  pt.put( "total_time", totalTime_.elapsed_time() );

  pt_.add_child( timeStamp, pt );

  switch( format_ ){
    case JSON:
      if( haveFileIO_ ) boost::property_tree::json_parser::write_json( logFileName_, pt_ );
      else              boost::property_tree::json_parser::write_json( std::cout,    pt_ );
      break;
    case XML:
     if( haveFileIO_ ) boost::property_tree::xml_parser::write_xml( logFileName_, pt_ );
     else              boost::property_tree::xml_parser::write_xml( std::cout,    pt_ );
      break;
  }
}

//-----------------------------------------------------------------------------
