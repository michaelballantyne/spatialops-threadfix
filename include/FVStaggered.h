#ifndef FVStaggered_h
#define FVStaggered_h

#include <FVStaggeredTypes.h>
#include <FVStaggeredInterpolant.h>
#include <FVStaggeredGradient.h>
#include <FVStaggeredDivergence.h>
#include <FVStaggeredScratch.h>


/**
 *  @todo Need to build & test scratch operator assemblers.
 *  @todo Need to fix linear system interface for spatial operators.
 *  @todo get funky interpolant operators working.  This is going to be painful. 
 */

namespace SpatialOps{

  //------------------------------------------------------------------

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  // Implementation of *= operator for surface scalar to surface vector fields
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::SSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::SSurfX,FVStaggered::NoGhost>::
  operator*=( const FVStaggered::SSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::SSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld *= *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::XSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::XSurfX,FVStaggered::NoGhost>::
  operator*=( const FVStaggered::XSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::XSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld *= *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::YSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::YSurfX,FVStaggered::NoGhost>::
  operator*=( const FVStaggered::YSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::YSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld *= *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::ZSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::ZSurfX,FVStaggered::NoGhost>::
  operator*=( const FVStaggered::ZSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::ZSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld *= *irhs;
    return *this;
  }

  //------------------------------------------------------------------

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  // Implementation of /= operator for surface scalar to surface vector fields
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::SSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::SSurfX,FVStaggered::NoGhost>::
  operator/=( const FVStaggered::SSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::SSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld /= *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::XSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::XSurfX,FVStaggered::NoGhost>::
  operator/=( const FVStaggered::XSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::XSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld /= *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::YSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::YSurfX,FVStaggered::NoGhost>::
  operator/=( const FVStaggered::YSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::YSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld /= *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::ZSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::ZSurfX,FVStaggered::NoGhost>::
  operator/=( const FVStaggered::ZSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::ZSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld /= *irhs;
    return *this;
  }

  //------------------------------------------------------------------

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  // Implementation of += operator for surface scalar to surface vector fields
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::SSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::SSurfX,FVStaggered::NoGhost>::
  operator+=( const FVStaggered::SSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::SSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld += *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::XSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::XSurfX,FVStaggered::NoGhost>::
  operator+=( const FVStaggered::XSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::XSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld += *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::YSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::YSurfX,FVStaggered::NoGhost>::
  operator+=( const FVStaggered::YSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::YSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld += *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::ZSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::ZSurfX,FVStaggered::NoGhost>::
  operator+=( const FVStaggered::ZSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::ZSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld += *irhs;
    return *this;
  }

  //------------------------------------------------------------------

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  // Implementation of -= operator for surface scalar to surface vector fields
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::SSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::SSurfX,FVStaggered::NoGhost>::
  operator-=( const FVStaggered::SSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::SSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld -= *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::XSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::XSurfX,FVStaggered::NoGhost>::
  operator-=( const FVStaggered::XSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::XSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld -= *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::YSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::YSurfX,FVStaggered::NoGhost>::
  operator-=( const FVStaggered::YSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::YSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld -= *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::ZSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::ZSurfX,FVStaggered::NoGhost>::
  operator-=( const FVStaggered::ZSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::ZSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld -= *irhs;
    return *this;
  }

  //------------------------------------------------------------------

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  // Implementation of = operator for surface scalar to surface vector fields
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::SSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::SSurfX,FVStaggered::NoGhost>::
  operator=( const FVStaggered::SSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::SSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::SSurfY,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::SSurfY,FVStaggered::NoGhost>::
  operator=( const FVStaggered::SSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::SSurfField::const_iterator irhs = rhs.begin();
    if( this->get_entries_per_comp(0)>0 )
      irhs += rhs.get_entries_per_comp(0);
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::SSurfZ,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::SSurfZ,FVStaggered::NoGhost>::
  operator=( const FVStaggered::SSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::SSurfField::const_iterator irhs=rhs.begin() ;
    if( this->get_entries_per_comp(0)>0 )
      irhs += rhs.get_entries_per_comp(0) + rhs.get_entries_per_comp(1);
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  //------------------------------------------------------------------

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::XSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::XSurfX,FVStaggered::NoGhost>::
  operator=( const FVStaggered::XSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::XSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::XSurfY,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::XSurfY,FVStaggered::NoGhost>::
  operator=( const FVStaggered::XSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::XSurfField::const_iterator irhs = rhs.begin();
    if( this->get_entries_per_comp(0)>0 )
      irhs += rhs.get_entries_per_comp(0);
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::XSurfZ,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::XSurfZ,FVStaggered::NoGhost>::
  operator=( const FVStaggered::XSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::XSurfField::const_iterator irhs=rhs.begin() ;
    if( this->get_entries_per_comp(0)>0 )
      irhs += rhs.get_entries_per_comp(0) + rhs.get_entries_per_comp(1);
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  //------------------------------------------------------------------

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::YSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::YSurfX,FVStaggered::NoGhost>::
  operator=( const FVStaggered::YSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::YSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::YSurfY,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::YSurfY,FVStaggered::NoGhost>::
  operator=( const FVStaggered::YSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::YSurfField::const_iterator irhs = rhs.begin();
    if( this->get_entries_per_comp(0)>0 )
      irhs += rhs.get_entries_per_comp(0);
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::YSurfZ,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::YSurfZ,FVStaggered::NoGhost>::
  operator=( const FVStaggered::YSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::YSurfField::const_iterator irhs=rhs.begin() ;
    if( this->get_entries_per_comp(0)>0 )
      irhs += rhs.get_entries_per_comp(0) + rhs.get_entries_per_comp(1);
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  //------------------------------------------------------------------

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::ZSurfX,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::ZSurfX,FVStaggered::NoGhost>::
  operator=( const FVStaggered::ZSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::ZSurfField::const_iterator irhs=rhs.begin();
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::ZSurfY,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::ZSurfY,FVStaggered::NoGhost>::
  operator=( const FVStaggered::ZSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::ZSurfField::const_iterator irhs = rhs.begin();
    if( this->get_entries_per_comp(0)>0 )
      irhs += rhs.get_entries_per_comp(0);
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  template<>
  template<>
  inline SpatialField<LinAlgTrilinos,FVStaggered::ZSurfZ,FVStaggered::NoGhost>&
  SpatialField<LinAlgTrilinos,FVStaggered::ZSurfZ,FVStaggered::NoGhost>::
  operator=( const FVStaggered::ZSurfField& rhs )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    FVStaggered::ZSurfField::const_iterator irhs=rhs.begin() ;
    if( this->get_entries_per_comp(0)>0 )
      irhs += rhs.get_entries_per_comp(0) + rhs.get_entries_per_comp(1);
    for( ; ifld!=iflde; ++ifld, ++irhs )
      *ifld = *irhs;
    return *this;
  }

  //------------------------------------------------------------------

}// namespace SpatialOps
#endif
