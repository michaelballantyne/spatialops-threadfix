#include "BoxFilter.h"
#include <spatialops/structured/FVStaggeredFieldTypes.h>

namespace SpatialOps{
namespace structured{


  /*
   *  This function builds a vector of fields that are windowed into
   *  the supplied "src" field such that they can be used in a filtering operator.
   *  This should work in 1D, 2D or 3D.
   */
  template<typename FieldT>
  void
  build_src_fields( const FieldT& src,
                    std::vector<FieldT>& fields )
  {
    const MemoryWindow& ws = src.window_with_ghost();
    const size_t ihi = ws.glob_dim(0)>1 ? 3 : 1;
    const size_t jhi = ws.glob_dim(1)>1 ? 3 : 1;
    const size_t khi = ws.glob_dim(2)>1 ? 3 : 1;
    const IntVec of( ihi>1 ? 2 : 0, jhi>1 ? 2 : 0, khi>1 ? 2 : 0 );
    for( size_t k=0; k<khi; ++k ){
      for( size_t j=0; j<jhi; ++j ){
        for( size_t i=0; i<ihi; ++i ){
          fields.push_back( FieldT( MemoryWindow( ws.glob_dim(),
                                                  ws.offset()+IntVec(i,j,k),
                                                  ws.extent()-of,
                                                  ws.has_bc(0), ws.has_bc(1), ws.has_bc(2) ),
                                    &src[0],
                                    ExternalStorage ) );
        }
      }
    }
    assert( fields.size() == ihi*jhi*khi );
  }

  //-----------------------------------------------------------------

  template<typename FieldT>
  void
  BoxFilter<FieldT>::apply_to_field( const FieldT& src, FieldT& dest ) const
  {
    const MemoryWindow& w_dest = dest.window_with_ghost();

    srcFields_.clear();  srcIters_.clear();
    build_src_fields<FieldT>( src, srcFields_ );
    for( typename std::vector<FieldT>::const_iterator is=srcFields_.begin(); is!=srcFields_.end(); ++is ){
      srcIters_.push_back( is->begin() );
    }

    IntVec of, ex;
    for( int i=0; i<3; ++i ){
      of[i] = w_dest.glob_dim(i)>1 ? 1 : 0;
      ex[i] = w_dest.glob_dim(i)>1 ? -2 : 0;
    }

    // create the destination field memory window
    FieldT d( MemoryWindow( w_dest.glob_dim(),
                            w_dest.offset()+of,
                            w_dest.extent()+ex,
                            w_dest.has_bc(0), w_dest.has_bc(1), w_dest.has_bc(2) ),
              &dest[0],
              ExternalStorage );

    const double fac = 1.0 / double(srcFields_.size());
    typename FieldT::iterator id=d.begin();
    const typename FieldT::iterator ide=d.end();
    for( ; id!=ide; ++id ){
      *id = 0.0;
      typename std::vector<ConstFieldIter>::iterator isi=srcIters_.begin();
      const typename std::vector<ConstFieldIter>::const_iterator isie=srcIters_.end();
      for( ; isi!=isie; ++isi ){
        *id += **isi;
        ++(*isi);  // increment this source iterator to the next point
      }
      *id *= fac;
    }
  }

  //==================================================================
  // Explicit template instantiation
  //
# define DECLARE_FILTER_VARIANTS( VOLT )                \
  template class BoxFilter<VOLT>;                       \
  template class BoxFilter<FaceTypes<VOLT>::XFace>;     \
  template class BoxFilter<FaceTypes<VOLT>::YFace>;     \
  template class BoxFilter<FaceTypes<VOLT>::ZFace>;

  DECLARE_FILTER_VARIANTS( SVolField )
  DECLARE_FILTER_VARIANTS( XVolField )
  DECLARE_FILTER_VARIANTS( YVolField )
  DECLARE_FILTER_VARIANTS( ZVolField )
  //
  //==================================================================

} // namespace structured
} // namespace SpatialOps