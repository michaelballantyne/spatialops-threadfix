#ifndef FVToolsTemplates_h
#define FVToolsTemplates_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>

#include <spatialops/structured/MemoryWindow.h>

#include <vector>
#include <set>

namespace SpatialOps{
namespace structured{

  /** @file FVToolsTemplates.h
   *  @brief Provides function templates useful for field and operator
   *  creation for structured meshes.  Specialize these templates to
   *  your needs.
   */

  //------------------------------------------------------------------

  template<typename FieldT>
  int get_nx_with_ghost( const int nxNoGhost, const bool hasPlusFaceX )
  {
    int nx = nxNoGhost;
    if( nxNoGhost>1 ) nx += 2*FieldT::Ghost::NGHOST;
    if( hasPlusFaceX && (0==FieldT::Location::FaceDir::value) ) nx+=1;
    return nx;
  }

  template<typename FieldT>
  int get_ny_with_ghost( const int nyNoGhost, const bool hasPlusFaceY )
  {
    int ny = nyNoGhost;
    if( nyNoGhost>1 ) ny += 2*FieldT::Ghost::NGHOST;
    if( hasPlusFaceY && (0==FieldT::Location::FaceDir::value) ) ny+=1;
    return ny;
  }

  template<typename FieldT>
  int get_nz_with_ghost( const int nzNoGhost, const bool hasPlusFaceZ )
  {
    int nz = nzNoGhost;
    if( nzNoGhost>1 ) nz += 2*FieldT::Ghost::NGHOST;
    if( hasPlusFaceZ && (0==FieldT::Location::FaceDir::value) ) nz+=1;
    return nz;
  }

  template<typename FieldT>
  IntVec get_dim_with_ghost( const IntVec& dimNoGhost,
                             const bool hasPlusFaceX, const bool hasPlusFaceY, const bool hasPlusFaceZ )
  {
    return IntVec( get_nx_with_ghost<FieldT>(dimNoGhost[0],hasPlusFaceX),
                   get_ny_with_ghost<FieldT>(dimNoGhost[1],hasPlusFaceY),
                   get_nz_with_ghost<FieldT>(dimNoGhost[2],hasPlusFaceZ) );
  }

  template<typename FieldT>
  int get_ntot_with_ghost( const IntVec& dimNoGhost,
                           const bool hasPlusFaceX, const bool hasPlusFaceY, const bool hasPlusFaceZ )
  {
    const IntVec dim( get_dim_with_ghost<FieldT>( dimNoGhost, hasPlusFaceX, hasPlusFaceY, hasPlusFaceZ ) );
    return dim[0]*dim[1]*dim[2];
  }


  /**
   *  \brief Obtain the memory window for a field on a patch that is a single, contiguous memory block
   *  \param localDim  the dimensionality of the field in consideration (without ghost cells)
   */
  template<typename FieldT>
  MemoryWindow
  get_window_with_ghost( const IntVec& localDim, const bool hasPlusFaceX, const bool hasPlusFaceY, const bool hasPlusFaceZ )
  {
    const IntVec dim( get_nx_with_ghost<FieldT>(localDim[0],hasPlusFaceX),
                      get_ny_with_ghost<FieldT>(localDim[1],hasPlusFaceY),
                      get_nz_with_ghost<FieldT>(localDim[2],hasPlusFaceZ) );
    return MemoryWindow( dim );
  }

  //------------------------------------------------------------------

  /**
   *  \brief Obtain the memory window for a field on a patch that is a subset of a larger memory block
   *
   *  \param globalDim the global dimensionality of the memory block (without ghost cells)
   *  \param localDim  the dimensionality of the field in consideration (without ghost cells)
   *  \param offset    the offset (start ijk index) of the local field in the global address space.
   */
  template<typename FieldT>
  MemoryWindow
  get_window_with_ghost( const IntVec& globalDim, const IntVec& localDim, const IntVec& offset,
                         const bool hasPlusFaceX, const bool hasPlusFaceY, const bool hasPlusFaceZ )
  {
    const IntVec dimLoc( get_nx_with_ghost<FieldT>(localDim[0],hasPlusFaceX),
                         get_ny_with_ghost<FieldT>(localDim[1],hasPlusFaceY),
                         get_nz_with_ghost<FieldT>(localDim[2],hasPlusFaceZ) );
    const IntVec dimGlob( get_nx_with_ghost<FieldT>(globalDim[0],hasPlusFaceX),
                          get_ny_with_ghost<FieldT>(globalDim[1],hasPlusFaceY),
                          get_nz_with_ghost<FieldT>(globalDim[2],hasPlusFaceZ) );
    return MemoryWindow( dimGlob, offset, dimLoc );
  }

  //------------------------------------------------------------------

  /**
   *  @brief obtain the number of points for this field type in the requested direction
   *
   *  @param dir The direction.  See XDIR, YDIR, ZDIR.
   *  @param dim The dimensions of the mesh in each direction.
   *  @param hasPlusFace bool indicating if the patch has a physical
   *         boundary in the dir direction.  In that case, the fields
   *         have an extra point.
   *
   *  \deprecated{ this should phase out, as it does not allow functionality with MemoryWindow}
   */
  template<typename FieldT>
  inline int npts( size_t dir, const std::vector<int>& dim, const bool hasPlusFace )
  {
    size_t n = dim[dir];
    if( n>1 ){
      n += 2*FieldT::Ghost::NGHOST;
      if( hasPlusFace && (int(dir) == FieldT::Location::FaceDir::value) ){
        ++n;
      }
    }
    return n;
  }

  template<> inline int npts<double>( size_t dir, const std::vector<int>& dim, const bool hasPlusFace )
  {
    return dim[dir];
  }


  /**
   * @brief get the total number of points (including ghost cells) for
   * a field in the x-direction.
   *
   * @param dim A vector containing the number of cells in each
   * coordinate direction.  This is a three-component vector.
   *
   * @param hasPlusXSideFaces A boolean flag to indicate if this patch
   * is on a +x side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.
   *
   *  \deprecated{ this should phase out, as it does not allow functionality with MemoryWindow}
   */
  template<typename FieldT>
  inline size_t get_nx( const std::vector<int>& dim,
                        const bool hasPlusXSideFaces )
  {
    return npts<FieldT>( XDIR::value, dim, hasPlusXSideFaces );
  }

  /**
   * @brief get the total number of points (including ghost cells) for
   * a field in the y-direction.
   *
   * @param dim A vector containing the number of cells in each
   * coordinate direction.  This is a three-component vector.
   *
   * @param hasPlusYSideFaces A boolean flag to indicate if this patch
   * is on a +y side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.
   *
   *  \deprecated{ this should phase out, as it does not allow functionality with MemoryWindow}
   */
  template<typename FieldT> inline size_t get_ny( const std::vector<int>& dim,
                                                  const bool hasPlusYSideFaces )
  {
    return npts<FieldT>( YDIR::value, dim, hasPlusYSideFaces );
  }


  /**
   * @brief get the total number of points (including ghost cells) for
   * a field in the z-direction.
   *
   * @param dim A vector containing the number of cells in each
   * coordinate direction.  This is a three-component vector.
   *
   * @param hasPlusZSideFaces A boolean flag to indicate if this patch
   * is on a +z side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.
   *
   *  \deprecated{ this should phase out, as it does not allow functionality with MemoryWindow}
   */
  template<typename FieldT> inline size_t get_nz( const std::vector<int>& dim,
                                                  const bool hasPlusZSideFaces )
  {
    return npts<FieldT>( ZDIR::value, dim, hasPlusZSideFaces );
  }

  //==================================================================

  /**
   * @brief get the total number of points in a field, including ghost
   * cells.
   *
   * @param dim A vector containing the number of cells in each
   * coordinate direction.  This is a three-component vector.
   *
   * @param hasPlusXSideFaces A boolean flag to indicate if this patch
   * is on a +x side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.
   *
   * @param hasPlusYSideFaces A boolean flag to indicate if this patch
   * is on a +y side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.

   * @param hasPlusZSideFaces A boolean flag to indicate if this patch
   * is on a +z side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.
   *
   * @todo Remove default values.  This is very dangerous for parallel
   * computations to have a default value for the + side information.
   *
   *  \deprecated{ this should phase out, as it does not allow functionality with MemoryWindow}
   */
  template<typename FieldT> int get_n_tot( const std::vector<int>& dim,
                                           const bool hasPlusXSideFaces=true,
                                           const bool hasPlusYSideFaces=true,
                                           const bool hasPlusZSideFaces=true )
  {
    return get_nx<FieldT>(dim,hasPlusXSideFaces)
         * get_ny<FieldT>(dim,hasPlusYSideFaces)
         * get_nz<FieldT>(dim,hasPlusZSideFaces);
  }

  //====================================================================

  // intended for local use only.
  void _ghost_set_( const int ng,
                    const int nxt, const int nyt, const int nzt,
                    const IntVec& dim,
                    const bool hasPlusXSideFaces,
                    const bool hasPlusYSideFaces,
                    const bool hasPlusZSideFaces,
                    size_t& ix,
                    std::set<size_t>& ghostSet );

  //==================================================================

  /**
   *  @brief Obtain the set of indices corresponding to ghost cells
   *  for this field.
   *
   * @param hasPlusXSideFaces A boolean flag to indicate if this patch
   * is on a +x side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.
   *
   * @param hasPlusYSideFaces A boolean flag to indicate if this patch
   * is on a +y side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.

   * @param hasPlusZSideFaces A boolean flag to indicate if this patch
   * is on a +z side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.
   *
   * @todo Remove default values.  This is very dangerous for parallel
   * computations to have a default value for the + side information.
   *
   *  \todo JCS: need to rework this.  It only works for the "old"
   *  memory view.  Need to integrate this with the MemoryWindow
   *  stuff.
   */
  template<typename FieldT> 
  const std::set<size_t> get_ghost_set( const IntVec& dim,
                                        const bool hasPlusXSideFaces=true,
                                        const bool hasPlusYSideFaces=true,
                                        const bool hasPlusZSideFaces=true )
  {
    std::set<size_t> ghostSet;
    ghostSet.clear();
    size_t ix=0;
    _ghost_set_( FieldT::Ghost::NGHOST,
                 get_nx_with_ghost<FieldT>(dim[0],hasPlusXSideFaces),
                 get_ny_with_ghost<FieldT>(dim[1],hasPlusYSideFaces),
                 get_nz_with_ghost<FieldT>(dim[2],hasPlusZSideFaces),
                 dim,
                 hasPlusXSideFaces, hasPlusYSideFaces, hasPlusZSideFaces,
                 ix,
                 ghostSet );
    return ghostSet;
  }

  //==================================================================

  /**
   *  @struct IndexTriplet
   *  @brief  Holds the ijk index.
   */
  struct IndexTriplet
  {
    IndexTriplet( const int ii, const int jj, const int kk ): i(ii), j(jj), k(kk){}
    IndexTriplet(){ i=j=k=-1; }
    IndexTriplet( const IndexTriplet& x ){ i=x.i; j=x.j; k=x.k; }
    IndexTriplet& operator=(const IndexTriplet& x){ i=x.i; j=x.j; k=x.k; return *this; }
    int& operator[](const int dim)      { switch(dim) { case 0: return i; case 1: return j; case 2: return k; } }
    int  operator[](const int dim) const{ switch(dim) { case 0: return i; case 1: return j; case 2: return k; } }
    int i,j,k;
  };

  //==================================================================

  /**
   *  @brief Use this to transform a flat index to i,j,k indices.
   */
  template<typename FieldT>
  struct flat2ijk
  {
    static IndexTriplet value( const IntVec& dim, const int ix,
                               const bool hasPlusXSideFaces=true,
                               const bool hasPlusYSideFaces=true,
                               const bool hasPlusZSideFaces=true );
  };

  //==================================================================

  /**
   *  @brief Use this to transform i,j,k indices to a flat index.
   */
  template<typename FieldT>
  struct ijk2flat
  {
    static int value( const IntVec& dim, const IndexTriplet& ixt,
                      const bool hasPlusXSideFaces=true,
                      const bool hasPlusYSideFaces=true,
                      const bool hasPlusZSideFaces=true );
  };

  //====================================================================

  template<typename FieldT>
  inline IndexTriplet
  flat2ijk<FieldT>::value( const IntVec& dim, const int ix,
                           const bool hasPlusXSideFaces, const bool hasPlusYSideFaces, const bool hasPlusZSideFaces )
  {
    IndexTriplet triplet;

    const int nxt = get_nx_with_ghost<FieldT>(dim[0],hasPlusXSideFaces);
    const int nyt = get_ny_with_ghost<FieldT>(dim[1],hasPlusXSideFaces);
    triplet.i = ix%nxt;
    triplet.j = ix/nxt % nyt;
    triplet.k = ix/(nxt*nyt);

    return triplet;
  }

  //==================================================================

  template<typename FieldT>
  inline int
  ijk2flat<FieldT>::value( const IntVec& dim, const IndexTriplet& triplet,
                           const bool hasPlusXSideFaces, const bool hasPlusYSideFaces, const bool hasPlusZSideFaces )
  {
    const int nxt = get_nx_with_ghost<FieldT>(dim[0],hasPlusXSideFaces);
    const int nyt = get_ny_with_ghost<FieldT>(dim[1],hasPlusYSideFaces);
      
    return
      triplet.i +
      triplet.j * nxt +
      triplet.k * nxt*nyt;
  }

  //==================================================================

}// namespace structured
}// namespace SpatialOps

#endif
