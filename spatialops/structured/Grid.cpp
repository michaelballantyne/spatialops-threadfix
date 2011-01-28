#include "Grid.h"
#include <spatialops/structured/FVStaggeredFieldTypes.h>

namespace SpatialOps{
namespace structured{

  // the amount that the volume coordinate for the mesh is shifted
  template< typename CoordT, typename StagT >
  struct StagShift{
    static double value(){ return 0.0; }
  };

  template< typename DirT >
  struct StagShift<DirT,DirT>{
    static double value(){ return -0.5; }
  };

  template< typename CoordT, typename FaceT >
  struct FaceShift{
    static double value(){ return 0.5; }
  };

  template< typename DirT >
  struct FaceShift<DirT,DirT>{
    // e.g. x-coordinate on an x-face.
    static double value(){ return 0.0; }
  };

  template< typename CoordT, typename StagT, typename FaceT >
  struct multiplier
  {
    static double value(){ return FaceShift<CoordT,FaceT>::value() + StagShift<CoordT,StagT>::value(); }
  };


  template< typename CoordT, typename FieldT >
  double get_offset( const double spc )
  {
    typedef typename FieldT::Location::StagLoc StagLoc;
    typedef typename FieldT::Location::FaceDir FaceDir;
    return spc * multiplier<CoordT,StagLoc,FaceDir>::value();
  }


  Grid::Grid( const IntVec npts,
              const std::vector<double>& length )
    : npts_( npts ),
      length_( length )
  {
    assert( length.size() == 3 );
    for( size_t i=0; i<3; ++i ){
      spacing_.push_back( length[i] / double(npts[i]) );
    }
  }

  template<typename CoordT> unsigned int get_dir();
  template<> unsigned int get_dir<XDIR>(){ return 0; }
  template<> unsigned int get_dir<YDIR>(){ return 1; }
  template<> unsigned int get_dir<ZDIR>(){ return 2; }

  //------------------------------------------------------------------

  template< typename CoordT >
  double Grid::spacing() const
  {
    return spacing_[ get_dir<CoordT>() ];
  }

  //------------------------------------------------------------------

  template< typename CoordT, typename FieldT >
  void Grid::set_coord( FieldT& f ) const
  {
    const unsigned int dir = get_dir<CoordT>();
    const double offset = get_offset<CoordT,FieldT>( spacing_[dir] );

    typedef typename FieldT::iterator FieldIter;

    FieldIter iter=f.begin();

    const MemoryWindow& mwInterior = f.window_without_ghost();
    const MemoryWindow& mw         = f.window_with_ghost();

    const IntVec lo(0,0,0);
    const IntVec hi( mw.extent() );

    const int ixOffset = mwInterior.offset(dir);
    IntVec ix;
    for( ix[2]=lo[2]; ix[2]<hi[2]; ++ix[2] ){
      for( ix[1]=lo[1]; ix[1]<hi[1]; ++ix[1] ){
        for( ix[0]=lo[0]; ix[0]<hi[0]; ++ix[0] ){
          *iter = spacing_[dir] * (ix[dir]-ixOffset) + offset;
          ++iter;
        }
      }
    }
  }

  //==================================================================
  // Explicit template instantiation
  //
# define DECLARE_COORD( FIELD )                                 \
  template void Grid::set_coord< XDIR, FIELD >( FIELD& ) const;	\
  template void Grid::set_coord< YDIR, FIELD >( FIELD& ) const;	\
  template void Grid::set_coord< ZDIR, FIELD >( FIELD& ) const;

# define DECLARE( VOL )                                 \
  DECLARE_COORD( VOL );                                 \
  DECLARE_COORD( FaceTypes<VOL>::XFace );               \
  DECLARE_COORD( FaceTypes<VOL>::YFace );               \
  DECLARE_COORD( FaceTypes<VOL>::ZFace );

  DECLARE( SVolField );
  DECLARE( XVolField );
  DECLARE( YVolField );
  DECLARE( ZVolField );

  template double Grid::spacing<XDIR>() const;
  template double Grid::spacing<YDIR>() const;
  template double Grid::spacing<ZDIR>() const;
  //
  //==================================================================

} // namespace structured
} // namespace SpatialOps
