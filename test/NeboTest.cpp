#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>

#include <spatialops/Nebo.h>
#include <test/TestHelper.h>
#include <test/FieldHelper.h>

#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace SpatialOps;
using namespace structured;

#define MANUAL(EXPR)                                            \
    {                                                           \
        const int xLength = window.extent(0);                   \
        const int yLength = window.extent(1);                   \
        const int zLength = window.extent(2);                   \
                                                                \
        for(int kk = 0; kk < zLength; kk++) {                   \
            for(int jj = 0; jj < yLength; jj++) {               \
                for(int ii = 0; ii < xLength; ii++) {           \
                    ref(ii, jj, kk) = EXPR;                     \
                };                                              \
            };                                                  \
        };                                                      \
    }

#define INPUT1 input1(ii, jj, kk)
#define INPUT2 input2(ii, jj, kk)
#define INPUT3 input3(ii, jj, kk)

int main( int iarg, char* carg[] )
{
    int nx, ny, nz;
    bool bcplus[] = { false, false, false };

    po::options_description desc("Supported Options");
    desc.add_options()
        ( "help", "print help message\n" )
        ( "nx",   po::value<int>(&nx)->default_value(11), "number of points in x-dir for base mesh" )
        ( "ny",   po::value<int>(&ny)->default_value(11), "number of points in y-dir for base mesh" )
        ( "nz",   po::value<int>(&nz)->default_value(11), "number of points in z-dir for base mesh" )
        ( "bcx",  "physical boundary on +x side?" )
        ( "bcy",  "physical boundary on +y side?" )
        ( "bcz",  "physical boundary on +z side?" );

    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    if( args.count("bcx") ) bcplus[0] = true;
    if( args.count("bcy") ) bcplus[1] = true;
    if( args.count("bcz") ) bcplus[2] = true;

    if( args.count("help") ){
      cout << desc << endl
           << "Examples:" << endl
           << " test_nebo --nx 5 --ny 10 --nz 3 --bcx" << endl
           << " test_nebo --bcx --bcy --bcz" << endl
           << " test_nebo --nx 50 --bcz" << endl
           << endl;
      return -1;
    }

    typedef SVolField Field;

    const MemoryWindow window( get_window_with_ghost<Field>(IntVec(nx,ny,nz),bcplus[0],bcplus[1],bcplus[2]) );

    Field input1( window, NULL );
    Field input2( window, NULL );
    Field input3( window, NULL );
    Field test( window, NULL );
    Field ref( window, NULL );

    const int total = nx * ny * nz;

    initialize_field(input1, 0.0);
    initialize_field(input1, total);
    initialize_field(input1, 2 * total);

    TestHelper status(true);

    test <<= 0.0;
    MANUAL(0.0);
    status( ref == test, "scalar assignment test" );

    test <<= max(input1, 0.0);
    MANUAL(std::max(INPUT1, 0.0));
    status( ref == test, "max test" );

    test <<= min(input1, 0.0);
    MANUAL(std::min(INPUT1, 0.0));
    status( ref == test, "min test" );

    if( status.ok() ) {
      cout << "ALL TESTS PASSED :)" << endl;
      return 0;
    }
    else {
        cout << "******************************" << endl
             << " At least one test FAILED! :(" << endl
             << "******************************" << endl;
        return -1;
    };
}
