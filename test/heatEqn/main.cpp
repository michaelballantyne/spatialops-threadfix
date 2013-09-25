#include <iostream>
using std::cout;
using std::endl;

//--- SpatialOps includes ---//
#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/Nebo.h>
#include <spatialops/structured/Grid.h>

#include <spatialops/structured/stencil/FVStaggeredOperatorTypes.h>
#include <spatialops/structured/stencil/StencilBuilder.h>

typedef SpatialOps::structured::SVolField   CellField;
typedef SpatialOps::structured::SSurfXField XSideField;
typedef SpatialOps::structured::SSurfYField YSideField;
typedef SpatialOps::structured::SSurfZField ZSideField;

typedef SpatialOps::structured::BasicOpTypes<CellField>::GradX      GradX;
typedef SpatialOps::structured::BasicOpTypes<CellField>::InterpC2FX InterpX;
typedef SpatialOps::structured::BasicOpTypes<CellField>::DivX       DivX;

typedef SpatialOps::structured::BasicOpTypes<CellField>::GradY      GradY;
typedef SpatialOps::structured::BasicOpTypes<CellField>::InterpC2FY InterpY;
typedef SpatialOps::structured::BasicOpTypes<CellField>::DivY       DivY;

typedef SpatialOps::structured::BasicOpTypes<CellField>::GradZ      GradZ;
typedef SpatialOps::structured::BasicOpTypes<CellField>::InterpC2FZ InterpZ;
typedef SpatialOps::structured::BasicOpTypes<CellField>::DivZ       DivZ;


//--- local includes ---//
#include "tools.h"

//-- boost includes ---//
#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace po = boost::program_options;
namespace SS = SpatialOps::structured;

int main( int iarg, char* carg[] )
{
  size_t ntime;
  std::vector<int> npts(3,1);
  std::vector<double> length(3,1.0);
  std::vector<double> spacing(3,1.0);

  // parse the command line options input describing the problem
  {
    po::options_description desc("Supported Options");
    desc.add_options()
      ( "help", "print help message" )
      ( "ntime", po::value<size_t>(&ntime)->default_value(1000), "number of 'iterations'" )
      //( "ntime", po::value<size_t>(&ntime), "number of 'iterations'" )
      ( "nx", po::value<int>(&npts[0])->default_value(10), "Grid in x" )
      ( "ny", po::value<int>(&npts[1])->default_value(10), "Grid in y" )
      ( "nz", po::value<int>(&npts[2])->default_value(10), "Grid in z" )
      ( "Lx", po::value<double>(&length[0])->default_value(1.0),"Length in x")
      ( "Ly", po::value<double>(&length[1])->default_value(1.0),"Length in y")
      ( "Lz", po::value<double>(&length[2])->default_value(1.0),"Length in z");

    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    if (args.count("help")) {
      cout << desc << "\n";
      return 1;
    }
  }

  cout << " [nx,ny,nz] = [" << npts[0] << "," << npts[1] << "," << npts[2] << "]" << endl
       << " ntime = " << ntime << endl
#     ifdef ENABLE_THREADS
       << " NTHREADS = " << NTHREADS << endl
#     endif
       << endl;

  // set mesh spacing (uniform, structured mesh)
  for( size_t i=0; i<3; ++i )
    spacing[i] = length[i]/double(npts[i]);

  // set face areas
  std::vector<double> area(3,1.0);
  area[0] = spacing[1]*spacing[2];
  area[1] = spacing[0]*spacing[2];
  area[2] = spacing[0]*spacing[1];

  // build the spatial operators
  SpatialOps::OperatorDatabase sodb;
  SpatialOps::structured::build_stencils( npts[0],   npts[1],   npts[2],
                                          length[0], length[1], length[2],
                                          sodb );

  // grab pointers to the operators
  const GradX* const gradx = sodb.retrieve_operator<GradX>();
  const GradY* const grady = sodb.retrieve_operator<GradY>();
  const GradZ* const gradz = sodb.retrieve_operator<GradZ>();

  const DivX* const divx = sodb.retrieve_operator<DivX>();
  const DivY* const divy = sodb.retrieve_operator<DivY>();
  const DivZ* const divz = sodb.retrieve_operator<DivZ>();

  const InterpX* const interpx = sodb.retrieve_operator<InterpX>();
  const InterpY* const interpy = sodb.retrieve_operator<InterpY>();
  const InterpZ* const interpz = sodb.retrieve_operator<InterpZ>();

  // build fields
  const SS::GhostData ghost(1);

  const SS::BoundaryCellInfo cellBC = SS::BoundaryCellInfo::build< CellField>(true,true,true);

  const SS::MemoryWindow vwindow( SS::get_window_with_ghost(npts,ghost,cellBC) );

  CellField temperature( vwindow, cellBC, ghost, NULL );
  CellField thermCond  ( vwindow, cellBC, ghost, NULL );
  CellField rhoCp      ( vwindow, cellBC, ghost, NULL );
  CellField xcoord     ( vwindow, cellBC, ghost, NULL );
  CellField ycoord     ( vwindow, cellBC, ghost, NULL );
  CellField zcoord     ( vwindow, cellBC, ghost, NULL );
  CellField rhs        ( vwindow, cellBC, ghost, NULL );

  SS::Grid grid( npts, length );
  grid.set_coord<SpatialOps::XDIR>( xcoord );
  grid.set_coord<SpatialOps::YDIR>( ycoord );
  grid.set_coord<SpatialOps::ZDIR>( zcoord );

  const SS::BoundaryCellInfo xBC = SS::BoundaryCellInfo::build<XSideField>(true,true,true);
  const SS::BoundaryCellInfo yBC = SS::BoundaryCellInfo::build<YSideField>(true,true,true);
  const SS::BoundaryCellInfo zBC = SS::BoundaryCellInfo::build<ZSideField>(true,true,true);
  const SS::MemoryWindow xwindow( SS::get_window_with_ghost(npts,ghost,xBC) );
  const SS::MemoryWindow ywindow( SS::get_window_with_ghost(npts,ghost,yBC) );
  const SS::MemoryWindow zwindow( SS::get_window_with_ghost(npts,ghost,zBC) );

  XSideField xflux( xwindow, xBC, ghost, NULL );
  YSideField yflux( ywindow, yBC, ghost, NULL );
  ZSideField zflux( zwindow, zBC, ghost, NULL );

  {
    using namespace SpatialOps;
    temperature <<= sin( xcoord ) + cos( ycoord ) + sin( zcoord );
    thermCond <<= xcoord + ycoord + zcoord;
    rhoCp <<= 1.0;
  }

  try{
    cout << "beginning 'timestepping'" << endl;

    const boost::posix_time::ptime time_start( boost::posix_time::microsec_clock::universal_time() );

    // mimic the effects of solving this PDE in time.
    for( size_t itime=0; itime<ntime; ++itime ){

      using namespace SpatialOps;

      if( npts[0]>1 ) calculate_flux( *gradx, *interpx, temperature, thermCond, xflux );
      if( npts[1]>1 ) calculate_flux( *grady, *interpy, temperature, thermCond, yflux );
      if( npts[2]>1 ) calculate_flux( *gradz, *interpz, temperature, thermCond, zflux );

      rhs <<= 0.0;
      if( npts[0]>1 ) calculate_rhs( *divx, xflux, rhoCp, rhs );
      if( npts[1]>1 ) calculate_rhs( *divy, yflux, rhoCp, rhs );
      if( npts[2]>1 ) calculate_rhs( *divz, zflux, rhoCp, rhs );

      // ordinarily at this point we would use rhs to update
      // temperature.  however, there are some other complicating
      // factors like setting boundary conditions that we have neglected
      // here to simplify things.

      //    cout << itime+1 << " of " << ntime << endl;
    }

    const boost::posix_time::ptime time_end( boost::posix_time::microsec_clock::universal_time() );
    const boost::posix_time::time_duration time_dur = time_end - time_start;

    cout << "done" << endl << endl
         << "time taken: "
         << time_dur.total_microseconds()*1e-6
         << endl << endl;


    return 0;
  }
  catch( std::exception& err ){
    cout << "ERROR!" << endl;
  }
  return -1;
}
