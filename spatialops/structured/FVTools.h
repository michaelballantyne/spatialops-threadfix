#ifndef FVToolsTemplates_h
#define FVToolsTemplates_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>
#include <spatialops/SpatialOpsTools.h>

#include <spatialops/structured/MemoryWindow.h>

#include <vector>
#include <set>

namespace SpatialOps{
namespace structured{

  /** \file FVTools.h
   *
   *  \addtogroup structured
   *  @{
   *  \addtogroup fields
   *  @{
   *
   *  \brief Provides function templates useful for field and operator
   *  creation for structured meshes.  Specialize these templates to
   *  your needs.
   */

  //------------------------------------------------------------------

  /**
   *  \fn int nextra()
   *
   *  \brief obtain the number of extra points in this field in the
   *         given direction, relative to a cell centered field
   *
   *  \param hasPlusFace true if the (+) side of this patch is a physical boundary.
   *
   *  \tparam FieldT the type of field in question
   *  \tparam DirT   the direction in question
   */
  template< typename FieldT, typename DirT >
  struct ExtraPoints{
    /**
     *  \brief obtain the number of extra points.
     *  \param hasPlusFace indicates if there is a physical boundary on the (+) side of this patch.
     */
    static unsigned int value( const bool hasPlusFace ){
      return IsSameType<typename FieldT::Location::RootLoc,DirT>::result && hasPlusFace ? 1 : 0;
    }
  };

  template< typename FieldT > struct ExtraPoints<FieldT,NODIR>; // invalid usage

  //------------------------------------------------------------------

  namespace detail{
    template< typename RootT, typename DirT > struct LocDirComp             { enum{ value=0 }; };
    template<                 typename DirT > struct LocDirComp<DirT ,DirT >{ enum{ value=-1}; };
    template<                               > struct LocDirComp<NODIR,NODIR>{ enum{ value=0 }; };
  }

  /**
   *  \struct IndexStagger
   *
   *  \brief Obtain the index value for how far the given field type
   *         is staggered relative to a scalar cell centered variable.
   *         Nominally 0 or -1.
   *
   *  \tparam FieldT the type of field in consideration.
   *  \tparam DirT   the direction we are interested in.
   *
   *  Usage:
   *  \code
   *    const int n = IndexStagger<FieldT,DirT>::value;
   *  \endcode
   */
  template< typename FieldT, typename DirT >
  struct IndexStagger{
    enum{ value = detail::LocDirComp<typename FieldT::Location::RootLoc,DirT>::value };
  };

  template< typename FieldT > struct IndexStagger<FieldT,NODIR>; // invalid usage

  //------------------------------------------------------------------

  /**
   *  \fn unsigned int stride( const MemoryWindow& mw )
   *
   *  \brief Obtain the stride (flat index) in the requested direction
   *         for the given field
   *
   *  \tparam DirT the direction of interest
   *  \param mw the MemoryWindow associated with the field
   */
  template< typename DirT > unsigned int stride( const MemoryWindow& mw );
  
  template<> inline unsigned int stride<XDIR>( const MemoryWindow& mw ){ return 1; }
  template<> inline unsigned int stride<YDIR>( const MemoryWindow& mw ){ return mw.extent(0); }
  template<> inline unsigned int stride<ZDIR>( const MemoryWindow& mw ){ return mw.extent(0)*mw.extent(1); }

  //------------------------------------------------------------------

  /**
   *  \fn int get_nx_with_ghost( const int, const bool )
   *
   *  \brief obtain the number of points in the x direction
   *
   *  \param nxNoGhost number of points in the x-direction excluding
   *    ghost cells
   *
   *  \param hasPlusFaceX - is the (+) side a physical domain
   *    boundary?  If so, x-face fields get an extra entry.
   *
   *  \result the number of points in the x-direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want NX for.
   */
  template<typename FieldT>
  int get_nx_with_ghost( const int nxNoGhost, const bool hasPlusFaceX )
  {
    int nx = nxNoGhost;
    if( nxNoGhost>1 ) nx += ExtraPoints<FieldT,XDIR>::value(hasPlusFaceX) + 2*FieldT::Ghost::NGHOST;
    return nx;
  }

  /**
   *  \fn int get_ny_with_ghost( const int, const bool )
   *
   *  \brief obtain the number of points in the y direction
   *
   *  \param nyNoGhost number of points in the y-direction excluding
   *    ghost cells
   *
   *  \param hasPlusFaceY - is the (+) side a physical domain
   *    boundary?  If so, y-face fields get an extra entry.
   *
   *  \result the number of points in the y-direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want NY for.
   */
  template<typename FieldT>
  int get_ny_with_ghost( const int nyNoGhost, const bool hasPlusFaceY )
  {
    int ny = nyNoGhost;
    if( nyNoGhost>1 ) ny += ExtraPoints<FieldT,YDIR>::value(hasPlusFaceY) + 2*FieldT::Ghost::NGHOST;
    return ny;
  }

  /**
   *  \fn int get_nz_with_ghost( const int, const bool )
   *
   *  \brief obtain the number of points in the z direction
   *
   *  \param nzNoGhost number of points in the z-direction excluding
   *    ghost cells
   *
   *  \param hasPlusFaceZ - is the (+) side a physical domain
   *    boundary?  If so, z-face fields get an extra entry.
   *
   *  \result the number of points in the z-direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want NZ for.
   */
  template<typename FieldT>
  int get_nz_with_ghost( const int nzNoGhost, const bool hasPlusFaceZ )
  {
    int nz = nzNoGhost;
    if( nzNoGhost>1 ) nz += ExtraPoints<FieldT,ZDIR>::value(hasPlusFaceZ) + 2*FieldT::Ghost::NGHOST;
    return nz;
  }

  /**
   *  \fn int get_dim_with_ghost( const IntVec&, const bool, const bool, const bool )
   *
   *  \brief obtain the number of points in each direction for the given field type
   *
   *  \param dimNoGhost number of points in each direction excluding
   *    ghost cells
   *
   *  \param hasPlusFaceX - is the (+) side a physical domain
   *    boundary?  If so, x-face fields get an extra entry.
   *
   *  \param hasPlusFaceY - is the (+) side a physical domain
   *    boundary?  If so, y-face fields get an extra entry.
   *
   *  \param hasPlusFaceZ - is the (+) side a physical domain
   *    boundary?  If so, z-face fields get an extra entry.
   *
   *  \result the number of points in each direction, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want (NX,NY,NZ) for.
   */
  template<typename FieldT>
  IntVec get_dim_with_ghost( const IntVec& dimNoGhost,
                             const bool hasPlusFaceX, const bool hasPlusFaceY, const bool hasPlusFaceZ )
  {
    return IntVec( get_nx_with_ghost<FieldT>(dimNoGhost[0],hasPlusFaceX),
                   get_ny_with_ghost<FieldT>(dimNoGhost[1],hasPlusFaceY),
                   get_nz_with_ghost<FieldT>(dimNoGhost[2],hasPlusFaceZ) );
  }

  /**
   *  \fn int get_ntot_with_ghost( const IntVec&, const bool, const bool, const bool )
   *  \brief get the total number of points including ghost cells
   *
   *  \param dimNoGhost number of points in each direction excluding
   *    ghost cells
   *
   *  \param hasPlusFaceX - is the (+) side a physical domain
   *    boundary?  If so, x-face fields get an extra entry.
   *
   *  \param hasPlusFaceY - is the (+) side a physical domain
   *    boundary?  If so, y-face fields get an extra entry.
   *
   *  \param hasPlusFaceZ - is the (+) side a physical domain
   *    boundary?  If so, z-face fields get an extra entry.
   *
   *  \result the total number of points in the field, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want (NX,NY,NZ) for.
   */
   template<typename FieldT>
  int get_ntot_with_ghost( const IntVec& dimNoGhost,
                           const bool hasPlusFaceX, const bool hasPlusFaceY, const bool hasPlusFaceZ )
  {
    const IntVec dim( get_dim_with_ghost<FieldT>( dimNoGhost, hasPlusFaceX, hasPlusFaceY, hasPlusFaceZ ) );
    return dim[0]*dim[1]*dim[2];
  }


  /**
   *  \fn MemoryWindow get_window_with_ghost( const IntVec&, const bool, const bool, const bool )
   *  \brief Obtain the memory window for a field on a patch that is a single, contiguous memory block
   *
   *  \param dimNoGhost number of points in each direction excluding
   *    ghost cells
   *
   *  \param hasPlusFaceX - is the (+) side a physical domain
   *    boundary?  If so, x-face fields get an extra entry.
   *
   *  \param hasPlusFaceY - is the (+) side a physical domain
   *    boundary?  If so, y-face fields get an extra entry.
   *
   *  \param hasPlusFaceZ - is the (+) side a physical domain
   *    boundary?  If so, z-face fields get an extra entry.
   *
   *  \result the total number of points in the field, including ghost cells.
   *
   *  \tparam FieldT the type of field that we want (NX,NY,NZ) for.
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

  template<> inline int get_nx_with_ghost<double>( const int nxNoGhost, const bool hasPlusFaceX ){ return nxNoGhost; }
  template<> inline int get_ny_with_ghost<double>( const int nyNoGhost, const bool hasPlusFaceY ){ return nyNoGhost; }
  template<> inline int get_nz_with_ghost<double>( const int nzNoGhost, const bool hasPlusFaceZ ){ return nzNoGhost; }

  //------------------------------------------------------------------

  /**
   *  \brief Obtain the memory window for a field on a patch that is a subset of a larger memory block
   *
   *  \param globalDim the global dimensionality of the memory block (without ghost cells)
   *  \param localDim  the dimensionality of the field in consideration (without ghost cells)
   *  \param offset    the offset (start ijk index) of the local field in the global address space.
   *
   *  \param hasPlusFaceX - is the (+) side a physical domain
   *    boundary?  If so, x-face fields get an extra entry.
   *
   *  \param hasPlusFaceY - is the (+) side a physical domain
   *    boundary?  If so, y-face fields get an extra entry.
   *
   *  \param hasPlusFaceZ - is the (+) side a physical domain
   *    boundary?  If so, z-face fields get an extra entry.
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
   *  \brief Obtain the set of indices corresponding to ghost cells
   *  for this field.
   *
   * \param hasPlusXSideFaces A boolean flag to indicate if this patch
   * is on a +x side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.
   *
   * \param hasPlusYSideFaces A boolean flag to indicate if this patch
   * is on a +y side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.

   * \param hasPlusZSideFaces A boolean flag to indicate if this patch
   * is on a +z side physical boundary.  If so, then it is assumed
   * that there is an extra face on that side of the domain, and face
   * variable dimensions will be modified accordingly.
   *
   * \todo Remove default values.  This is very dangerous for parallel
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
   *  \brief Use this to transform a flat index to i,j,k indices.
   */
  template<typename FieldT>
  struct flat2ijk
  {
    static IntVec value( const IntVec& dim, const int ix,
                         const bool hasPlusXSideFaces=true,
                         const bool hasPlusYSideFaces=true,
                         const bool hasPlusZSideFaces=true )
    {
      IntVec triplet(0,0,0);
      const int nxt = get_nx_with_ghost<FieldT>(dim[0],hasPlusXSideFaces);
      const int nyt = get_ny_with_ghost<FieldT>(dim[1],hasPlusYSideFaces);
      triplet[0] = ix%nxt;
      triplet[1] = ix/nxt % nyt;
      triplet[2] = ix/(nxt*nyt);
      return triplet;
    }
  };

  //==================================================================

  /**
   *  \brief Use this to transform i,j,k indices to a flat index.
   */
  template<typename FieldT>
  struct ijk2flat
  {
    static int value( const IntVec& dim,
                      const IntVec& triplet,
                      const bool hasPlusXSideFaces=true,
                      const bool hasPlusYSideFaces=true,
                      const bool hasPlusZSideFaces=true )
    {
      const int nxt = get_nx_with_ghost<FieldT>(dim[0],hasPlusXSideFaces);
      const int nyt = get_ny_with_ghost<FieldT>(dim[1],hasPlusYSideFaces);
      return
        triplet[0] +
        triplet[1] * nxt +
        triplet[2] * nxt*nyt;
    }
  };

  //==================================================================

  /**
   *  @}
   *  @}
   */

}// namespace structured
}// namespace SpatialOps

#endif
