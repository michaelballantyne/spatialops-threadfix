#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;

//--- SpatialOps includes ---//
#include <spatialops/structured/FVStaggered.h>
#include <spatialops/structured/Grid.h>
#ifdef ENABLE_THREADS
#include <spatialops/ThreadPool.h>
#endif

#include <util/TimeLogger.h>

typedef SpatialOps::SVolField   CellField;
typedef SpatialOps::SSurfXField XSideField;
typedef SpatialOps::SSurfYField YSideField;
typedef SpatialOps::SSurfZField ZSideField;

typedef SpatialOps::BasicOpTypes<CellField>::GradX      GradX;
typedef SpatialOps::BasicOpTypes<CellField>::InterpC2FX InterpX;
typedef SpatialOps::BasicOpTypes<CellField>::DivX       DivX;

typedef SpatialOps::BasicOpTypes<CellField>::GradY      GradY;
typedef SpatialOps::BasicOpTypes<CellField>::InterpC2FY InterpY;
typedef SpatialOps::BasicOpTypes<CellField>::DivY       DivY;

typedef SpatialOps::BasicOpTypes<CellField>::GradZ      GradZ;
typedef SpatialOps::BasicOpTypes<CellField>::InterpC2FZ InterpZ;
typedef SpatialOps::BasicOpTypes<CellField>::DivZ       DivZ;


//-- boost includes ---//
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace SpatialOps;

//==============================================================================

void
initialize_mask_points( const Grid& grid,
                        std::vector<IntVec>& xminus,
                        std::vector<IntVec>& xplus,
                        std::vector<IntVec>& yminus,
                        std::vector<IntVec>& yplus,
                        std::vector<IntVec>& zminus,
                        std::vector<IntVec>& zplus )
{
  for( int k=-1; k<grid.extent(2); ++k ){
    for( int  j=-1; j<=grid.extent(1); ++j ){
      xminus.push_back( IntVec(-1, j, k) );
      xplus .push_back( IntVec(grid.extent(0), j, k) );
    }
  }
  for( int k=-1; k<grid.extent(2); ++k ){
    for( int i=-1; i<=grid.extent(0); ++i ){
      xminus.push_back( IntVec(i,-1,k) );
      xplus .push_back( IntVec(i,grid.extent(1), k) );
    }
  }
  for( int j=-1; j<grid.extent(1); ++j ){
    for( int i=-1; i<=grid.extent(0); ++i ){
      xminus.push_back( IntVec(i,j,-1) );
      xplus .push_back( IntVec(i,j,grid.extent(2)) );
    }
  }
}

//==============================================================================

int main( int iarg, char* carg[] )
{
  double dt, tend;
  int nthreads = NTHREADS;
  IntVec npts;
  DoubleVec length;
  std::string logFileName;

  // parse the command line options input describing the problem
  {
    po::options_description desc("Supported Options");
    desc.add_options()
      ( "help", "print help message" )
      ( "dt",   po::value<double>(&dt  )->default_value(2e-4), "timestep" )
      ( "tend", po::value<double>(&tend)->default_value(0.04 ), "end time" )
      ( "nx", po::value<int>(&npts[0])->default_value(32), "nx" )
      ( "ny", po::value<int>(&npts[1])->default_value(32), "ny" )
      ( "nz", po::value<int>(&npts[2])->default_value(32), "nz" )
      ( "Lx", po::value<double>(&length[0])->default_value(1.0),"Length in x")
      ( "Ly", po::value<double>(&length[1])->default_value(1.0),"Length in y")
      ( "Lz", po::value<double>(&length[2])->default_value(1.0),"Length in z")
      ( "nthreads", po::value<int>(&nthreads)->default_value(NTHREADS),"Number of threads (no effect unless configured with threads enabled)" )
      ( "timings-file-name", po::value<std::string>(&logFileName)->default_value("timings.log"), "Name for performance timings file" );

    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    if (args.count("help")) {
      cout << "\nSolves the 3D transient diffusion equation\n\n" << desc << "\n";
      return 1;
    }
  }

# ifdef ENABLE_THREADS
  set_soft_thread_count( nthreads );
# endif

  cout << " [nx,ny,nz] = [" << npts[0] << "," << npts[1] << "," << npts[2] << "]" << endl
       << " dt   = " << dt << endl
       << " tend = " << tend << endl
#      ifdef ENABLE_THREADS
       << " NTHREADS = " << get_soft_thread_count() << endl
#      endif
       << endl;

  // set mesh spacing (uniform, structured mesh)
  const DoubleVec spacing = length/npts;

  const Grid grid( npts, length );

  // build the spatial operators
  SpatialOps::OperatorDatabase sodb;
  SpatialOps::build_stencils( grid, sodb );

  const GradX& gradx = *sodb.retrieve_operator<GradX>();
  const GradY& grady = *sodb.retrieve_operator<GradY>();
  const GradZ& gradz = *sodb.retrieve_operator<GradZ>();

  const DivX& divx = *sodb.retrieve_operator<DivX>();
  const DivY& divy = *sodb.retrieve_operator<DivY>();
  const DivZ& divz = *sodb.retrieve_operator<DivZ>();

  const InterpX& interpx = *sodb.retrieve_operator<InterpX>();
  const InterpY& interpy = *sodb.retrieve_operator<InterpY>();
  const InterpZ& interpz = *sodb.retrieve_operator<InterpZ>();

  // build fields
  const GhostData ghost(1);

  const BoundaryCellInfo cellBC = BoundaryCellInfo::build< CellField>(true,true,true);

  const MemoryWindow vwindow( get_window_with_ghost(npts,ghost,cellBC) );

  CellField temperature( vwindow, cellBC, ghost, NULL );
  CellField thermCond  ( vwindow, cellBC, ghost, NULL );
  CellField rhoCp      ( vwindow, cellBC, ghost, NULL );
  CellField xcoord     ( vwindow, cellBC, ghost, NULL );
  CellField ycoord     ( vwindow, cellBC, ghost, NULL );
  CellField zcoord     ( vwindow, cellBC, ghost, NULL );
  CellField rhs        ( vwindow, cellBC, ghost, NULL );

  grid.set_coord<SpatialOps::XDIR>( xcoord );
  grid.set_coord<SpatialOps::YDIR>( ycoord );
  grid.set_coord<SpatialOps::ZDIR>( zcoord );

  const BoundaryCellInfo xBC = BoundaryCellInfo::build<XSideField>(true,true,true);
  const BoundaryCellInfo yBC = BoundaryCellInfo::build<YSideField>(true,true,true);
  const BoundaryCellInfo zBC = BoundaryCellInfo::build<ZSideField>(true,true,true);
  const MemoryWindow xwindow( get_window_with_ghost(npts,ghost,xBC) );
  const MemoryWindow ywindow( get_window_with_ghost(npts,ghost,yBC) );
  const MemoryWindow zwindow( get_window_with_ghost(npts,ghost,zBC) );

  XSideField xflux( xwindow, xBC, ghost, NULL );
  YSideField yflux( ywindow, yBC, ghost, NULL );
  ZSideField zflux( zwindow, zBC, ghost, NULL );

  //----------------------------------------------------------------------------
  // Build and initialize masks:
  std::vector<IntVec> xminusPts, xplusPts, yminusPts, yplusPts, zminusPts, zplusPts;

  initialize_mask_points( grid, xminusPts, xplusPts, yminusPts, yplusPts, zminusPts, zplusPts );

  SpatialMask<CellField> xminus( temperature, xminusPts );
  SpatialMask<CellField> xplus ( temperature, xplusPts  );
  SpatialMask<CellField> yminus( temperature, yminusPts );
  SpatialMask<CellField> yplus ( temperature, yplusPts  );
  SpatialMask<CellField> zminus( temperature, zminusPts );
  SpatialMask<CellField> zplus ( temperature, zplusPts  );

# ifdef ENABLE_CUDA
  // Masks are created on CPU so we need to explicitly transfer them to GPU
  xminus.add_consumer( GPU_INDEX );
  xplus .add_consumer( GPU_INDEX );
  yminus.add_consumer( GPU_INDEX );
  yplus .add_consumer( GPU_INDEX );
  zminus.add_consumer( GPU_INDEX );
  zplus .add_consumer( GPU_INDEX );
# endif
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // initial conditions
  thermCond   <<= xcoord + ycoord + zcoord;
  rhoCp       <<= 1.0;
  temperature <<= cond( xminus, 1.0 )( xplus, -1.0 )
                      ( yminus, 1.0 )( yplus, -1.0 )
                      ( zminus, 1.0 )( zplus, -1.0 )
                      ( sin( xcoord ) + cos( ycoord ) + sin( zcoord ) );
  //----------------------------------------------------------------------------

  try{
    cout << "beginning timestepping" << endl;

    TimeLogger logger(logFileName);

    // mimic the effects of solving this PDE in time.
    double t = 0;
    while( t < tend ){

      logger.start("flux");
      if( npts[0]>1 ) xflux <<= -interpx(thermCond) * gradx(temperature);
      if( npts[1]>1 ) yflux <<= -interpy(thermCond) * grady(temperature);
      if( npts[2]>1 ) zflux <<= -interpz(thermCond) * gradz(temperature);
      logger.stop("flux");

      rhs <<= 0.0;

      logger.start("rhs");
//      if( npts[0]>1 ) rhs <<= rhs + divx( interpx(thermCond) * gradx(temperature) ) / rhoCp;
//      if( npts[1]>1 ) rhs <<= rhs + divy( interpy(thermCond) * grady(temperature) ) / rhoCp;
//      if( npts[2]>1 ) rhs <<= rhs + divz( interpz(thermCond) * gradz(temperature) ) / rhoCp;
      if( npts[0]>1 ) rhs <<= rhs - divx( xflux ) / rhoCp;
      if( npts[1]>1 ) rhs <<= rhs - divy( yflux ) / rhoCp;
      if( npts[2]>1 ) rhs <<= rhs - divz( zflux ) / rhoCp;
      logger.stop("rhs");

      logger.start("update");
      temperature <<= temperature + dt * rhs;
      logger.stop("update");

      logger.start("bc");
      temperature <<= cond( xminus, 1.0 )( xplus, -1.0 )
                          ( yminus, 1.0 )( yplus, -1.0 )
                          ( zminus, 1.0 )( zplus, -1.0 )
                          ( temperature );
      logger.stop("bc");

      t += dt;
    }

    cout << "done" << endl << endl
         << "\n------------------------------------------------------------------\n"
         << "Timing Summary:\n"
         << "---------------\n"
         << std::setw(10) << "Flux"
         << std::setw(10)  << "RHS"
         << std::setw(10)  << "update"
         << std::setw(10)  << "bc"
         << std::setw(10) << "TOTAL"
         << std::endl
         << std::setw(10) << logger.timer("flux").elapsed_time()
         << std::setw(10) << logger.timer("rhs").elapsed_time()
         << std::setw(10) << logger.timer("update").elapsed_time()
         << std::setw(10) << logger.timer("bc").elapsed_time()
         << std::setw(10) << logger.total_time()
         << "\n------------------------------------------------------------------\n"
         << endl << endl;

    return 0;
  }
  catch( std::exception& err ){
    cout << "ERROR!" << endl;
  }
  return -1;
}
