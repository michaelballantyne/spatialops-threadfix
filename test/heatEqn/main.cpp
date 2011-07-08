#include <iostream>
using std::cout;
using std::endl;

//--- SpatialOps includes ---//
#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/FieldExpressions.h>
#include <spatialops/structured/Grid.h>

#ifdef LINALG_STENCIL
# include <spatialops/structured/stencil/FVStaggeredOperatorTypes.h>
# include <spatialops/structured/stencil/StencilBuilder.h>
#else
# include <spatialops/structured/matrix/FVStaggeredOperatorTypes.h>
# include <spatialops/structured/matrix/FVStaggeredInterpolant.h>
# include <spatialops/structured/matrix/FVStaggeredGradient.h>
# include <spatialops/structured/matrix/FVStaggeredDivergence.h>
#endif


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
#ifdef LINALG_STENCIL
  SpatialOps::structured::build_stencils( npts[0],   npts[1],   npts[2],
                                          length[0], length[1], length[2],
                                          sodb );
#else
  sodb.register_new_operator<GradX>( GradX::Assembler( spacing[0], npts, true, true, true ) );
  sodb.register_new_operator<GradY>( GradY::Assembler( spacing[1], npts, true, true, true ) );
  sodb.register_new_operator<GradZ>( GradZ::Assembler( spacing[2], npts, true, true, true ) );

  sodb.register_new_operator<DivX>( DivX::Assembler( npts, area[0], spacing[0], true, true, true ) );
  sodb.register_new_operator<DivY>( DivY::Assembler( npts, area[1], spacing[1], true, true, true ) );
  sodb.register_new_operator<DivZ>( DivZ::Assembler( npts, area[2], spacing[2], true, true, true ) );

  sodb.register_new_operator<InterpX>( InterpX::Assembler( npts, true, true, true ) );
  sodb.register_new_operator<InterpY>( InterpY::Assembler( npts, true, true, true ) );
  sodb.register_new_operator<InterpZ>( InterpZ::Assembler( npts, true, true, true ) );
#endif

  // grap pointers to the operators
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
  const SpatialOps::structured::MemoryWindow vwindow( SpatialOps::structured::get_window_with_ghost<CellField>(npts,true,true,true) );

  CellField temperature( vwindow, NULL );
  CellField thermCond  ( vwindow, NULL );
  CellField rhoCp      ( vwindow, NULL );
  CellField xcoord     ( vwindow, NULL );
  CellField ycoord     ( vwindow, NULL );
  CellField zcoord     ( vwindow, NULL );
  CellField rhs        ( vwindow, NULL );

  SpatialOps::structured::Grid grid( npts, length );
  grid.set_coord<SpatialOps::XDIR>( xcoord );
  grid.set_coord<SpatialOps::YDIR>( ycoord );
  grid.set_coord<SpatialOps::ZDIR>( zcoord );

  const SpatialOps::structured::MemoryWindow xwindow( SpatialOps::structured::get_window_with_ghost<XSideField>(npts,true,true,true) );
  const SpatialOps::structured::MemoryWindow ywindow( SpatialOps::structured::get_window_with_ghost<YSideField>(npts,true,true,true) );
  const SpatialOps::structured::MemoryWindow zwindow( SpatialOps::structured::get_window_with_ghost<ZSideField>(npts,true,true,true) );

  XSideField xflux( xwindow, NULL );
  YSideField yflux( ywindow, NULL );
  ZSideField zflux( zwindow, NULL );

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
