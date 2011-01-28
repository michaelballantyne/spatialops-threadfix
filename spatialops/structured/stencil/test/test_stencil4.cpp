#include <spatialops/FieldExpressionsExtended.h>
#include <spatialops/FieldReductions.h>

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/Grid.h>
#include <spatialops/structured/stencil/Stencil4.h>

#include <test/TestHelper.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace SpatialOps;
using namespace structured;

#include <stdexcept>
using std::cout;
using std::endl;



int main( int iarg, char* carg[] )
{
  int nx, ny, nz;
  bool bcplus[] = { false, false, false };

  {
    po::options_description desc("Supported Options");
    desc.add_options()
      ( "help", "print help message\n" )
      ( "nx",   po::value<int>(&nx)->default_value(5), "number of points in x-dir for base mesh" )
      ( "ny",   po::value<int>(&ny)->default_value(5), "number of points in y-dir for base mesh" )
      ( "nz",   po::value<int>(&nz)->default_value(5), "number of points in z-dir for base mesh" )
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
           << " test_stencil --nx 5 --ny 10 --nz 3 --bcx" << endl
           << " test_stencil --bcx --bcy --bcz" << endl
           << " test_stencil --nx 50 --bcz" << endl
           << endl;
      return -1;
    }
  }


  TestHelper status( true );
  const IntVec npts(nx,ny,nz);

  {
    std::string bcx = bcplus[0] ? "ON" : "OFF";
    std::string bcy = bcplus[1] ? "ON" : "OFF";
    std::string bcz = bcplus[2] ? "ON" : "OFF";
    cout << "Run information: " << endl
         << "  bcx    : " << bcx << endl
         << "  bcy    : " << bcy << endl
         << "  bcz    : " << bcz << endl
         << "  domain : " << npts << endl
         << endl;
  }

}
