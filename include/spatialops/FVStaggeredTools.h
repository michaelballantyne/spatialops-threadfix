#ifndef FVStaggeredTools_h
#define FVStaggeredTools_h

#include <spatialops/SpatialOpsConfigure.h>

#include <set>

#include <spatialops/FVStaggeredTypes.h>
#include <spatialops/FVToolsTemplates.h>

namespace SpatialOps{
namespace FVStaggered{

  //==================================================================

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

}// namespace FVStaggered
}// namespace SpatialOps

#endif
