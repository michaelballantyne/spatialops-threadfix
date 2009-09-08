#ifndef FVToolsTemplates_h
#define FVToolsTemplates_h

#include <spatialops/SpatialOpsConfigure.h>

#include <vector>
#include <set>

namespace SpatialOps{
namespace FVStaggered{

  /** @file FVToolsTemplates.h
   *  @brief Provides function templates useful for field and operator
   *  creation for structured meshes.  Specialize these templates to
   *  your needs.
   *
   *  @todo Move these tools out of the FVStaggered namespace.
   */

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
   */
  template<typename FieldT>
  inline int get_nx( const std::vector<int>& dim,
                     const bool hasPlusXSideFaces );

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
   */
  template<typename FieldT> inline int get_ny( const std::vector<int>& dim,
                                               const bool hasPlusYSideFaces );


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
   */
  template<typename FieldT> inline int get_nz( const std::vector<int>& dim,
                                               const bool hasPlusZSideFaces );

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
  void _ghost_set_( const int ngm, const int ngp,
                    const int nxt, const int nyt, const int nzt,
                    const std::vector<int>& dim,
                    const bool hasPlusXSideFaces,
                    const bool hasPlusYSideFaces,
                    const bool hasPlusZSideFaces,
                    int& ix,
                    std::set<int>& ghostSet );

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
   */
  template<typename FieldT> 
  const std::set<int> get_ghost_set( const std::vector<int>& dim,
                                     const bool hasPlusXSideFaces=true,
                                     const bool hasPlusYSideFaces=true,
                                     const bool hasPlusZSideFaces=true )
  {
    typedef typename FieldT::Ghost G;
    std::set<int> ghostSet;
    ghostSet.clear();
    int ix=0;
    _ghost_set_( G::NM, G::NP,
                 get_nx<FieldT>(dim,hasPlusXSideFaces),
                 get_ny<FieldT>(dim,hasPlusYSideFaces),
                 get_nz<FieldT>(dim,hasPlusZSideFaces),
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
    static IndexTriplet value( const std::vector<int>& dim, const int ix,
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
    static int value( const std::vector<int>& dim, const IndexTriplet& ixt,
                      const bool hasPlusXSideFaces=true,
                      const bool hasPlusYSideFaces=true,
                      const bool hasPlusZSideFaces=true );
  };

  //====================================================================

  template<typename FieldT>
  inline IndexTriplet
  flat2ijk<FieldT>::value( const std::vector<int>& dim, const int ix,
                           const bool hasPlusXSideFaces, const bool hasPlusYSideFaces, const bool hasPlusZSideFaces )
  {
    IndexTriplet triplet;

    const int nxt = get_nx<FieldT>(dim,hasPlusXSideFaces);
    const int nyt = get_ny<FieldT>(dim,hasPlusXSideFaces);
    triplet.i = ix%nxt;
    triplet.j = ix/nxt % nyt;
    triplet.k = ix/(nxt*nyt);

    return triplet;
  }

  //==================================================================

  template<typename FieldT>
  inline int
  ijk2flat<FieldT>::value( const std::vector<int>& dim, const IndexTriplet& triplet,
                           const bool hasPlusXSideFaces, const bool hasPlusYSideFaces, const bool hasPlusZSideFaces )
  {
    const int nxt = get_nx<FieldT>(dim,hasPlusXSideFaces);
    const int nyt = get_ny<FieldT>(dim,hasPlusYSideFaces);
      
    return
      triplet.i +
      triplet.j * nxt +
      triplet.k * nxt*nyt;
  }

  //==================================================================

}// namespace FVStaggered
}// namespace SpatialOps

#endif
