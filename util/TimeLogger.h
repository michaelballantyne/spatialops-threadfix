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


#ifndef TIMELOGGER_H_
#define TIMELOGGER_H_

#include <map>
#include <fstream>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/property_tree/ptree.hpp>


/**
 * \class Timer
 * \date Sep 8, 2014
 * \author James C. Sutherland
 * \brief Provides basic timing functionality
 *
 * This class provides stopwatch-like functionality.  You can repeatedly call
 * start/stop and the cumulative time will be reported, just like on a stopwatch.
 */
class Timer{
  boost::posix_time::ptime startTime_, stopTime_;
  double previous_;

public:

  Timer(){ reset(); }

  ~Timer(){}

  /** \brief (re)-start the timer */
  inline void start(){
    previous_ = elapsed_time();
    startTime_ = boost::posix_time::microsec_clock::universal_time();
  }

  /**
   *  \brief Stop the timer and return the time since the last call to start().
   *  This is different from elapsed_time().
   */
  inline double stop(){
    stopTime_ = boost::posix_time::microsec_clock::universal_time();
    return (stopTime_-startTime_).total_microseconds()*1e-6;
  }

  /** \brief reset the timer */
  inline void reset(){
    previous_ = 0;
    stopTime_ = startTime_ = boost::posix_time::microsec_clock::universal_time();
  }

  /**
   * \brief Return the total time for the timer.  This is the sum of all time
   *  intervals between calls to start/stop since the timer was constructed.
   *
   *  Note that if the timer is running (start() without a paired stop()),
   *  that time will not be incorporated into this call.
   */
  inline double elapsed_time() const{
    return previous_ + (stopTime_-startTime_).total_microseconds()*1e-6;
  }

};

/**
 *  \class  TimeLogger
 *  \date   Sep 8, 2014
 *  \author "James C. Sutherland"
 *
 *  The TimeLogger class provides a facility to log timings of calculations.
 *  It provides the ability to create timers labeled by strings for fine-grained
 *  timing. It also allows for log files to be written to disk and for existing
 *  log files to be appended, thus creating a time history of performance.
 *
 */
class TimeLogger
{
public:

  enum Format{
    XML,
    JSON
  };

  /**
   * @param logFileName the name of the logging file to use (existing files will be appended)
   * @param format either JSON or XML. JSON is default if no argument is provided
   */
  TimeLogger( const std::string logFileName,
              const Format format=JSON );

  /**
   * @brief Log to std::cout
   * @param format either JSON or XML. JSON is default if no argument is provided
   */
  TimeLogger( const Format format=JSON );

  ~TimeLogger();

  /** (re)start the given timer */
  void start( const std::string& label );

  /** \brief stop the given timer */
  double stop( const std::string& label );

  /** \brief reset the given timer */
  void reset( const std::string& label );

  /** \brief add a Timer to the logger */
  void add_entry( const std::string& label,
                  const Timer& timer );

  /** \brief Obtain the requested timer */
  const Timer& timer( const std::string& label ) const;

  /** \brief Obtain the elapsed time that this Logger has been in existence (time since construction) */
  double total_time();

private:
  std::string logFileName_;
  const bool haveFileIO_;
  const Format format_;
  typedef std::map<std::string,Timer>  Entries;
  Entries entries_;

  Timer totalTime_;

  boost::property_tree::ptree pt_;

  void write_entries();
};


#endif /* TIMELOGGER_H_ */
