#ifndef FVStaggeredTools_h
#define FVStaggeredTools_h

#include <set>

#include <spatialops/FVStaggeredTypes.h>

namespace SpatialOps{
namespace FVStaggered{

  //==================================================================

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
  template<typename FieldT> inline int get_nx( const std::vector<int>& dim,
                                               const bool hasPlusXSideFaces );

  template<> inline int get_nx<SVolField>( const std::vector<int>& dim,
                                           const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 ) return 1;
    return dim[0] + SVolField::Ghost::NM + SVolField::Ghost::NP;
  }
  template<> inline int get_nx<SVolRHS>( const std::vector<int>& dim,
                                         const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 ) return 1;
    return dim[0] + SVolRHS::Ghost::NM + SVolRHS::Ghost::NP;
  }
  template<> inline int get_nx<SSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 ) return 1;
    int npts = dim[0] + SSurfXField::Ghost::NM + SSurfXField::Ghost::NP;
    if( hasPlusXSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_nx<SSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    return dim[0] + SSurfYField::Ghost::NM + SSurfYField::Ghost::NP;
  }
  template<> inline int get_nx<SSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    return dim[0] + SSurfZField::Ghost::NM + SSurfZField::Ghost::NP;
  }

  template<> inline int get_nx<XVolField>( const std::vector<int>& dim,
                                           const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 ) return 1;
    int npts = dim[0] + XVolField::Ghost::NM + XVolField::Ghost::NP;
    if( hasPlusXSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_nx<XVolRHS>( const std::vector<int>& dim,
                                         const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 ) return 1;
    int npts = dim[0] + XVolRHS::Ghost::NM + XVolRHS::Ghost::NP;
    if( hasPlusXSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_nx<XSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 ) return 1;
    return dim[0] + XSurfXField::Ghost::NM + XSurfXField::Ghost::NP;
  }
  template<> inline int get_nx<XSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    int npts = dim[0] + XSurfYField::Ghost::NM + XSurfYField::Ghost::NP;
    if( hasPlusXSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_nx<XSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    int npts = dim[0] + XSurfZField::Ghost::NM + XSurfZField::Ghost::NP;
    if( hasPlusXSideFaces ) ++npts;
    return npts;
  }

  template<> inline int get_nx<YVolField>( const std::vector<int>& dim,
                                           const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    return dim[0] + YVolField::Ghost::NM + YVolField::Ghost::NP;
  }
  template<> inline int get_nx<YVolRHS>( const std::vector<int>& dim,
                                           const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    return dim[0] + YVolRHS::Ghost::NM + YVolRHS::Ghost::NP;
  }
  template<> inline int get_nx<YSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    int npts = dim[0] + YSurfXField::Ghost::NM + YSurfXField::Ghost::NP;
    if( hasPlusXSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_nx<YSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    return dim[0] + YSurfYField::Ghost::NM + YSurfYField::Ghost::NP;
  }
  template<> inline int get_nx<YSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[0] + YSurfZField::Ghost::NM + YSurfZField::Ghost::NP;
  }

  template<> inline int get_nx<ZVolField>( const std::vector<int>& dim,
                                           const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    return dim[0] + ZVolField::Ghost::NM + ZVolField::Ghost::NP;
  }
  template<> inline int get_nx<ZVolRHS>( const std::vector<int>& dim,
                                           const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    return dim[0] + ZVolRHS::Ghost::NM + ZVolRHS::Ghost::NP;
  }
  template<> inline int get_nx<ZSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    int npts = dim[0] + ZSurfXField::Ghost::NM + ZSurfXField::Ghost::NP;
    if( hasPlusXSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_nx<ZSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[0] + ZSurfYField::Ghost::NM + ZSurfYField::Ghost::NP;
  }
  template<> inline int get_nx<ZSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusXSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    return dim[0] + ZSurfZField::Ghost::NM + ZSurfZField::Ghost::NP;
  }


  //==================================================================
 
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

  template<> inline int get_ny<SVolField>( const std::vector<int>& dim,
                                           const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 ) return 1;
    return dim[1] + SVolField::Ghost::NM + SVolField::Ghost::NP;
  }
  template<> inline int get_ny<SVolRHS>( const std::vector<int>& dim,
                                           const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 ) return 1;
    return dim[1] + SVolRHS::Ghost::NM + SVolRHS::Ghost::NP;
  }
  template<> inline int get_ny<SSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    return dim[1] + SSurfXField::Ghost::NM + SSurfXField::Ghost::NP;
  }
  template<> inline int get_ny<SSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 ) return 1;
    int npts = dim[1] + SSurfYField::Ghost::NM + SSurfYField::Ghost::NP;
    if( hasPlusYSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_ny<SSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[1] + SSurfZField::Ghost::NM + SSurfZField::Ghost::NP;
  }

  template<> inline int get_ny<XVolField>( const std::vector<int>& dim,
                                           const bool hasPlusYSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    return dim[1] + XVolField::Ghost::NM + XVolField::Ghost::NP;
  }
  template<> inline int get_ny<XVolRHS>( const std::vector<int>& dim,
                                           const bool hasPlusYSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    return dim[1] + XVolRHS::Ghost::NM + XVolRHS::Ghost::NP;
  }
  template<> inline int get_ny<XSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    return dim[1] + XSurfXField::Ghost::NM + XSurfXField::Ghost::NP;
  }
  template<> inline int get_ny<XSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    int npts = dim[1] + XSurfYField::Ghost::NM + XSurfYField::Ghost::NP;
    if( hasPlusYSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_ny<XSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[1] + XSurfZField::Ghost::NM + XSurfZField::Ghost::NP;
  }

  template<> inline int get_ny<YVolField>( const std::vector<int>& dim,
                                           const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 ) return 1;
    int npts = dim[1] + YVolField::Ghost::NM + YVolField::Ghost::NP;
    if( hasPlusYSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_ny<YVolRHS>( const std::vector<int>& dim,
                                           const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 ) return 1;
    int npts = dim[1] + YVolRHS::Ghost::NM + YVolRHS::Ghost::NP;
    if( hasPlusYSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_ny<YSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 ) return 1;
    int npts = dim[1] + YSurfXField::Ghost::NM + YSurfXField::Ghost::NP;
    if( hasPlusYSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_ny<YSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 ) return 1;
    return dim[1] + YSurfYField::Ghost::NM + YSurfYField::Ghost::NP;
  }
  template<> inline int get_ny<YSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    int npts = dim[1] + YSurfZField::Ghost::NM + YSurfZField::Ghost::NP;
    if( hasPlusYSideFaces ) ++npts;
    return npts;
  }


  template<> inline int get_ny<ZVolField>( const std::vector<int>& dim,
                                           const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[1] + ZVolField::Ghost::NM + ZVolField::Ghost::NP;
  }
  template<> inline int get_ny<ZVolRHS>( const std::vector<int>& dim,
                                           const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[1] + ZVolRHS::Ghost::NM + ZVolRHS::Ghost::NP;
  }
  template<> inline int get_ny<ZSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[1] + ZSurfXField::Ghost::NM + ZSurfXField::Ghost::NP;
  }
  template<> inline int get_ny<ZSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    int npts = dim[1] + ZSurfYField::Ghost::NM + ZSurfYField::Ghost::NP;
    if( hasPlusYSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_ny<ZSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusYSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[1] + ZSurfZField::Ghost::NM + ZSurfZField::Ghost::NP;
  }

  //==================================================================

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

  template<> inline int get_nz<SVolField>( const std::vector<int>& dim,
                                           const bool hasPlusZSideFaces )
  {
    if( dim[2]<=1 ) return 1;
    return dim[2] + SVolField::Ghost::NM + SVolField::Ghost::NP;
  }
  template<> inline int get_nz<SVolRHS>( const std::vector<int>& dim,
                                           const bool hasPlusZSideFaces )
  {
    if( dim[2]<=1 ) return 1;
    return dim[2] + SVolRHS::Ghost::NM + SVolRHS::Ghost::NP;
  }
  template<> inline int get_nz<SSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    return dim[2] + SSurfXField::Ghost::NM + SSurfXField::Ghost::NP;
  }
  template<> inline int get_nz<SSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[2] + SSurfYField::Ghost::NM + SSurfYField::Ghost::NP;
  }
  template<> inline int get_nz<SSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[2]<=1 ) return 1;
    int npts = dim[2] + SSurfZField::Ghost::NM + SSurfZField::Ghost::NP;
    if( hasPlusZSideFaces ) ++npts;
    return npts;
  }

  template<> inline int get_nz<XVolField>( const std::vector<int>& dim,
                                           const bool hasPlusZSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    return dim[2] + XVolField::Ghost::NM + XVolField::Ghost::NP;
  }
  template<> inline int get_nz<XVolRHS>( const std::vector<int>& dim,
                                           const bool hasPlusZSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    return dim[2] + XVolRHS::Ghost::NM + XVolRHS::Ghost::NP;
  }
  template<> inline int get_nz<XSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    return dim[2] + XSurfXField::Ghost::NM + XSurfXField::Ghost::NP;
  }
  template<> inline int get_nz<XSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[2] + XSurfYField::Ghost::NM + XSurfYField::Ghost::NP;
  }
  template<> inline int get_nz<XSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    int npts = dim[2] + XSurfZField::Ghost::NM + XSurfZField::Ghost::NP;
    if( hasPlusZSideFaces ) ++npts;
    return npts;
  }

  template<> inline int get_nz<YVolField>( const std::vector<int>& dim,
                                           const bool hasPlusZSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[2] + YVolField::Ghost::NM + YVolField::Ghost::NP;
  }
  template<> inline int get_nz<YVolRHS>( const std::vector<int>& dim,
                                           const bool hasPlusZSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[2] + YVolRHS::Ghost::NM + YVolRHS::Ghost::NP;
  }
  template<> inline int get_nz<YSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[0]<=1 || dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[2] + YSurfXField::Ghost::NM + YSurfXField::Ghost::NP;
  }
  template<> inline int get_nz<YSurfYField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    return dim[2] + YSurfYField::Ghost::NM + YSurfYField::Ghost::NP;
  }
  template<> inline int get_nz<YSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    int npts = dim[2] + YSurfZField::Ghost::NM + YSurfZField::Ghost::NP;
    if( hasPlusZSideFaces ) ++npts;
    return npts;
  }

  template<> inline int get_nz<ZVolField>( const std::vector<int>& dim,
                                           const bool hasPlusZSideFaces )
  {
    if( dim[2]<=1 ) return 1;
    int npts = dim[2] + ZVolField::Ghost::NM + ZVolField::Ghost::NP;
    if( hasPlusZSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_nz<ZVolRHS>( const std::vector<int>& dim,
                                           const bool hasPlusZSideFaces )
  {
    if( dim[2]<=1 ) return 1;
    int npts = dim[2] + ZVolRHS::Ghost::NM + ZVolRHS::Ghost::NP;
    if( hasPlusZSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_nz<ZSurfXField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[0]<=1 || dim[2]<=1 ) return 1;
    int npts = dim[2] + ZSurfXField::Ghost::NM + ZSurfXField::Ghost::NP;
    if( hasPlusZSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_nz<ZSurfYField>( const std::vector<int>& dim,
                                                    const bool hasPlusZSideFaces )
  {
    if( dim[1]<=1 || dim[2]<=1 ) return 1;
    int npts = dim[2] + ZSurfYField::Ghost::NM + ZSurfYField::Ghost::NP;
    if( hasPlusZSideFaces ) ++npts;
    return npts;
  }
  template<> inline int get_nz<ZSurfZField>( const std::vector<int>& dim,
                                             const bool hasPlusZSideFaces )
  {
    if( dim[2]<=1 ) return 1;
    return dim[2] + ZSurfZField::Ghost::NM + ZSurfZField::Ghost::NP;
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

  //==================================================================

  // intended for local use only.
  inline void _ghost_set_( const int ngm, const int ngp,
                           const int nxt, const int nyt, const int nzt,
                           const std::vector<int>& dim,
                           const bool hasPlusXSideFaces,
                           const bool hasPlusYSideFaces,
                           const bool hasPlusZSideFaces,
                           int& ix,
                           std::set<int>& ghostSet )
  {
    const int ngxm = dim[0]>1 ? ngm : 0;
    const int ngxp = dim[0]>1 ? ngp : 0;
    const int ngym = dim[1]>1 ? ngm : 0;
    const int ngyp = dim[1]>1 ? ngp : 0;
    const int ngzm = dim[2]>1 ? ngm : 0;
    const int ngzp = dim[2]>1 ? ngp : 0;

    // -z side ghost layer
    if( dim[2]>1 ){
      for( int kg=0; kg<ngzm; ++kg )
        for( int j=0; j<nyt; ++j )
          for( int i=0; i<nxt; ++i )
            ghostSet.insert(ix++);
    }

    // z interior
    for( int k=ngzm; k<nzt-ngzp; ++k ){

      // -y side ghost layer
      if( dim[1]>1 ){
        for( int i=0; i<nxt; ++i )
          for( int jg=0; jg<ngym; ++jg )
            ghostSet.insert(ix++);
      }

      // y interior
      for( int j=ngym; j<nyt-ngyp; ++j ){
        // -x side ghost layer
        if( dim[0]>1 ) for( int ig=0; ig<ngxm; ++ig ) ghostSet.insert(ix++);
        // x interior
        ix+=nxt-ngxm-ngxp;
        // +x side ghost layer
        if( dim[0]>1 ) for( int ig=0; ig<ngxp; ++ig ) ghostSet.insert(ix++);
      }

      // +y side ghost layer
      if( dim[1]>1 ){
        for( int i=0; i<nxt; ++i )
          for( int jg=0; jg<ngyp; ++jg )
            ghostSet.insert(ix++);
      }
    }

    // +z side ghost layer
    if( dim[2]>1 ){
      for( int kg=0; kg<ngzp; ++kg )
        for( int i=0; i<nxt; ++i )
          for( int j=0; j<nyt; ++j )
            ghostSet.insert(ix++);
    }
  }

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
