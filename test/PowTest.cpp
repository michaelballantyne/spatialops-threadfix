#include <iostream>

//--- SpatialOps includes ---//
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/IntVec.h>
#include <spatialops/structured/MemoryWindow.h>
#include <spatialops/structured/stencil/FVStaggeredOperatorTypes.h>
#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/stencil/StencilBuilder.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/Nebo.h>

//-- boost includes ---//
#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace po = boost::program_options;

using namespace SpatialOps;

#define RUN_TEST(TEST,							\
		 TYPE)							\
  boost::posix_time::ptime start( boost::posix_time::microsec_clock::universal_time() ); \
  boost::posix_time::ptime end( boost::posix_time::microsec_clock::universal_time() ); \
  int ii = 0;								\
  start = boost::posix_time::microsec_clock::universal_time();		\
									\
  for(; ii < number_of_runs; ii++) {					\
    TEST;								\
  };									\
									\
  end = boost::posix_time::microsec_clock::universal_time();		\
  std::cout << TYPE;							\
  std::cout << " runs: ";						\
  std::cout << number_of_runs;						\
  std::cout << " result: ";						\
  std::cout << (end - start).total_microseconds()*1e-6;			\
  std::cout << std::endl;

template<typename FieldType>
inline void evaluate_square(FieldType & result,
                            FieldType const & phi,
                            int const number_of_runs) {
  RUN_TEST(result <<= square(phi),
           "square                  ");
};

template<typename FieldType>
inline void evaluate_pow_double(FieldType & result,
                                FieldType const & phi,
                                double const exp,
                                int const number_of_runs) {
  RUN_TEST(result <<= pow(phi, exp),
           "original pow with double");
};

template<typename FieldType>
inline void evaluate_pow_int(FieldType & result,
                             FieldType const & phi,
                             int const exp,
                             int const number_of_runs) {
  RUN_TEST(result <<= pow(phi, exp),
           "pow_int                 ");
};

int main(int iarg, char* carg[]) {
    typedef SVolField Field;

    bool test;
    int nx, ny, nz;
    int number_of_runs;
    int exp = 2;
#ifdef ENABLE_THREADS
  int thread_count;
#endif

    // parse the command line options input describing the problem
    {
        po::options_description desc("Supported Options");
	desc.add_options()
	  ( "help", "print help message" )
	  ( "nx", po::value<int>(&nx)->default_value(10), "Grid in x" )
	  ( "ny", po::value<int>(&ny)->default_value(10), "Grid in y" )
	  ( "nz", po::value<int>(&nz)->default_value(10), "Grid in z" )
#ifdef ENABLE_THREADS
          ( "tc", po::value<int>(&thread_count)->default_value(NTHREADS), "Number of threads for Nebo")
#endif
	  ( "runs", po::value<int>(&number_of_runs)->default_value(1), "Number of iterations of each test");

	po::variables_map args;
	po::store( po::parse_command_line(iarg,carg,desc), args );
	po::notify(args);

	if (args.count("help")) {
	    std::cout << desc << "\n";
	    return 1;
	}

#ifdef ENABLE_THREADS
    set_hard_thread_count(thread_count);
#endif
    }

    const GhostData ghost(1);
    const BoundaryCellInfo bc = BoundaryCellInfo::build<Field>(false,false,false);
    const MemoryWindow window( get_window_with_ghost(IntVec(nx,ny,nz),ghost,bc) );

    Field a( window, bc, ghost, NULL );
    Field r( window, bc, ghost, NULL );

    Field::iterator ia = a.begin();
    for(size_t kk = 0; kk < window.glob_dim(2); kk++)
      for(size_t jj = 0; jj < window.glob_dim(1); jj++)
        for(size_t ii = 0; ii < window.glob_dim(0); ii++, ++ia)
          *ia = ii + jj * 2 + kk * 4;

    evaluate_square(r, a, number_of_runs);
    evaluate_pow_int(r, a, exp, number_of_runs);
    evaluate_pow_double(r, a, exp, number_of_runs);

    return 0;
};
