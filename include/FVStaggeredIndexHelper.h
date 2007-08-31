#ifndef FVStaggeredIndexHelper_h
#define FVStaggeredIndexHelper_h


#include <FVStaggeredTypes.h>
#include <FVStaggeredTools.h>

namespace SpatialOps{
namespace FVStaggered{


  //==================================================================

  template<typename T1, typename T2, int IsSurf> struct DirSelector{};
  template<typename T1, typename T2> struct DirSelector<T1,T2,1>{ typedef T1 Type; };
  template<typename T1, typename T2> struct DirSelector<T1,T2,0>{ typedef T2 Type; };


  template< typename SrcFieldT, typename DestFieldT >
  class IndexHelper
  {
  public:
    IndexHelper( const std::vector<int>& dim );
    void get_cols( const int irow, std::vector<int>& cols ) const;
    int get_ncol() const;
    int get_nrow() const;

    int calculate_stride( const int irow, const int icol ) const;

  private:
    const std::vector<int>& dim_;
  };

  /**
   *  @todo Implement scalar volume -> staggered surface index helpers.
   */

  //==================================================================

  template<> class IndexHelper<SVolField,XSurfField>{};  // not yet implemented
  template<> class IndexHelper<SVolField,YSurfField>{};  // not yet implemented
  template<> class IndexHelper<SVolField,ZSurfField>{};  // not yet implemented

  //==================================================================

  inline bool is_valid_entry( const std::vector<int>& dim,
			      const IndexTriplet& ixt1,
			      const IndexTriplet& ixt2 );

  //==================================================================

  template<typename SrcFieldT, typename DestFieldT>
  bool is_in_bounds( const std::vector<int>& dim, const int ixs, const int ixd );

  //==================================================================

  template<typename SrcFieldT, typename DestFieldT>
  void shift_dest_index( const std::vector<int>& dim, const int irow, IndexTriplet& triplet );

  //==================================================================






  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //                          Implementation
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






  //==================================================================


  //------------------------------------------------------------------ 
  template<typename SrcFieldT, typename DestFieldT>
  IndexHelper<SrcFieldT,DestFieldT>::
  IndexHelper( const std::vector<int>& dim )
    : dim_( dim )
  {
  }
  //------------------------------------------------------------------ 
  template<typename SrcFieldT, typename DestFieldT>
  void
  IndexHelper<SrcFieldT,DestFieldT>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    // from the row (flat index of destination field), determine the
    // ijk indices for the destination field
    IndexTriplet ijk = flat2ijk<DestFieldT,DestFieldT::Location::IsSurface>::value( dim_, irow );

    shift_dest_index<SrcFieldT,DestFieldT>( dim_, irow, ijk );

    if( ijk.i<0 || ijk.j<0 || ijk.k<0 ) return;

    // from the ijk index, determine the flat index for the source field.
    const int icol = ijk2flat<SrcFieldT,SrcFieldT::Location::IsSurface>::value( dim_, ijk );

    const int stride = calculate_stride( irow, icol );

    if( is_in_bounds<SrcFieldT,DestFieldT>(dim_,icol,irow) &&
	is_in_bounds<SrcFieldT,DestFieldT>(dim_,icol+stride,irow) ){

      const IndexTriplet ijk2 = flat2ijk<SrcFieldT,SrcFieldT::Location::IsSurface>::value( dim_, icol+stride );

      if( is_valid_entry( dim_, ijk, ijk2 ) ){
	cols.push_back( icol        );
	cols.push_back( icol+stride );
      }

    }
  }
  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  int
  IndexHelper<SrcFieldT,DestFieldT>::
  get_ncol() const
  {
    return get_n_tot<SrcFieldT>(dim_);
  }
  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  int
  IndexHelper<SrcFieldT,DestFieldT>::
  get_nrow() const
  {
    return get_n_tot<DestFieldT>(dim_);
  }
  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  int
  IndexHelper<SrcFieldT,DestFieldT>::
  calculate_stride( const int irow, const int icol ) const
  {
    typedef typename DirSelector<SrcFieldT,DestFieldT,SrcFieldT::Location::IsSurface>::Type FieldDirSelect;

    int n=-1;
    switch( FieldDirSelect::Location::Dir::value ){
    case XDIR::value:
      n=1;
      break;
    case YDIR::value:
      n = get_nx<SrcFieldT>(dim_[0]);
      break;
    case ZDIR::value:
      n = get_nx<SrcFieldT>(dim_[0]) * get_ny<SrcFieldT>(dim_[1]);
      break;
    default:
      cout << "ERROR: field dir enum value=" << FieldDirSelect::Location::Dir::value << endl;
      assert(0);
    }
    return n;
  }
  template<> inline int IndexHelper<SVolField,XVolField>::calculate_stride( const int irow, const int icol ) const
  {
    return 1;
  }
  template<> inline int IndexHelper<SVolField,YVolField>::calculate_stride( const int irow, const int icol ) const
  {
    return get_nx<SVolField>(dim_[0]);
  }
  template<> inline int IndexHelper<SVolField,ZVolField>::calculate_stride( const int irow, const int icol ) const
  {
    return get_nx<SVolField>(dim_[0]) * get_ny<SVolField>(dim_[1]);
  }
  template<> inline int IndexHelper<SVolField,SSurfField>::calculate_stride( const int irow, const int icol ) const
  {
    const int nx1 = get_nx_x<SSurfField>(dim_);
    const int ny1 = get_ny_x<SSurfField>(dim_);
    const int nz1 = get_nz_x<SSurfField>(dim_);
    const int n1 = dim_[0]>1 ? nx1*ny1*nz1 : 0;

    if( irow<n1 ) return 1;

    const int nx2 = get_nx_y<SSurfField>(dim_);
    const int ny2 = get_ny_y<SSurfField>(dim_);
    const int nz2 = get_nz_y<SSurfField>(dim_);
    const int n2 = dim_[1]>1 ? n1 + nx2*ny2*nz2 : n1;

    if( irow<n2 ) return get_nx<SVolField>(dim_[0]);

    return get_nx<SVolField>(dim_[0])*get_nx<SVolField>(dim_[1]);
  }
  //------------------------------------------------------------------

  // not yet implemented:
  template<> inline int IndexHelper<XVolField,XSurfField>::calculate_stride( const int irow, const int icol ) const;
  template<> inline int IndexHelper<YVolField,YSurfField>::calculate_stride( const int irow, const int icol ) const;
  template<> inline int IndexHelper<ZVolField,ZSurfField>::calculate_stride( const int irow, const int icol ) const;


  //==================================================================

  bool is_valid_entry( const std::vector<int>& dim,
		       const IndexTriplet& ixt1,
		       const IndexTriplet& ixt2 )
  {
    if( std::abs(ixt2.i - ixt1.i) > 1 ) return false;
    if( std::abs(ixt2.j - ixt1.j) > 1 ) return false;
    if( std::abs(ixt2.k - ixt1.k) > 1 ) return false;
    return true;
  }
  
  //==================================================================

  template<typename SrcFieldT, typename DestFieldT>
  bool is_in_bounds( const std::vector<int>& dim, const int ixs, const int ixd )
  {
    const IndexTriplet it = flat2ijk<SrcFieldT,SrcFieldT::Location::IsSurface>::value(dim,ixs);

    if( it.i<0 || it.j<0 || it.k<0 ) return false;

    const int nx = get_nx<SrcFieldT>(dim[0]);
    const int ny = get_ny<SrcFieldT>(dim[1]);
    const int nz = get_nz<SrcFieldT>(dim[2]);

    if( it.i<nx && it.j<ny && it.k<nz )  return true;

    return false;
  }

  //==================================================================

  template<typename SrcFieldT, typename DestFieldT>
  void shift_dest_index( const std::vector<int>& dim, const int irow, IndexTriplet& triplet ){}

  template<> inline void shift_dest_index<SSurfXField,SVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (SSurfXField::Ghost::NM-SVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (SSurfXField::Ghost::NM-SVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (SSurfXField::Ghost::NM-SVolField::Ghost::NM);
  }
  template<> inline void shift_dest_index<SSurfYField,SVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (SSurfYField::Ghost::NM-SVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (SSurfYField::Ghost::NM-SVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (SSurfYField::Ghost::NM-SVolField::Ghost::NM);
  }
  template<> inline void shift_dest_index<SSurfZField,SVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (SSurfZField::Ghost::NM-SVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (SSurfZField::Ghost::NM-SVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (SSurfZField::Ghost::NM-SVolField::Ghost::NM);
  }

  template<> inline void shift_dest_index<XSurfXField,XVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (XSurfXField::Ghost::NM-XVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (XSurfXField::Ghost::NM-XVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (XSurfXField::Ghost::NM-XVolField::Ghost::NM);
  }
  template<> inline void shift_dest_index<XSurfYField,XVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (XSurfYField::Ghost::NM-XVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (XSurfYField::Ghost::NM-XVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (XSurfYField::Ghost::NM-XVolField::Ghost::NM);
  }
  template<> inline void shift_dest_index<XSurfZField,XVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (XSurfZField::Ghost::NM-XVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (XSurfZField::Ghost::NM-XVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (XSurfZField::Ghost::NM-XVolField::Ghost::NM);
  }

  template<> inline void shift_dest_index<YSurfXField,YVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (YSurfXField::Ghost::NM-YVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (YSurfXField::Ghost::NM-YVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (YSurfXField::Ghost::NM-YVolField::Ghost::NM);
  }
  template<> inline void shift_dest_index<YSurfYField,YVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (YSurfYField::Ghost::NM-YVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (YSurfYField::Ghost::NM-YVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (YSurfYField::Ghost::NM-YVolField::Ghost::NM);
  }
  template<> inline void shift_dest_index<YSurfZField,YVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (YSurfZField::Ghost::NM-YVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (YSurfZField::Ghost::NM-YVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (YSurfZField::Ghost::NM-YVolField::Ghost::NM);
  }

  template<> inline void shift_dest_index<ZSurfXField,ZVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (ZSurfXField::Ghost::NM-ZVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (ZSurfXField::Ghost::NM-ZVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (ZSurfXField::Ghost::NM-ZVolField::Ghost::NM);
  }
  template<> inline void shift_dest_index<ZSurfYField,ZVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (ZSurfYField::Ghost::NM-ZVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (ZSurfYField::Ghost::NM-ZVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (ZSurfYField::Ghost::NM-ZVolField::Ghost::NM);
  }
  template<> inline void shift_dest_index<ZSurfZField,ZVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (ZSurfZField::Ghost::NM-ZVolField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (ZSurfZField::Ghost::NM-ZVolField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (ZSurfZField::Ghost::NM-ZVolField::Ghost::NM);
  }

  template<> inline void shift_dest_index<SVolField,SSurfField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    const std::vector<int> entriesPerComp = get_entries_per_comp<SSurfField>(dim);
    const int n1 = entriesPerComp[0];
    const int n2 = entriesPerComp[1] + n1;

    if( irow < n1 ){
      if( dim[1]>1 ) triplet.j += (SVolField::Ghost::NM - SSurfField::Ghost::NM);
      if( dim[2]>1 ) triplet.k += (SVolField::Ghost::NM - SSurfField::Ghost::NM);
    }
    else if( irow < n2 ){
      if( dim[0]>1 ) triplet.i += (SVolField::Ghost::NM - SSurfField::Ghost::NM);
      if( dim[2]>1 ) triplet.k += (SVolField::Ghost::NM - SSurfField::Ghost::NM);
    }
    else{
      if( dim[0]>1 ) triplet.i += (SVolField::Ghost::NM - SSurfField::Ghost::NM);
      if( dim[1]>1 ) triplet.j += (SVolField::Ghost::NM - SSurfField::Ghost::NM);
    }
  }

  template<> inline void shift_dest_index<SVolField,SSurfXField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[1]>1 ) triplet.j += (SVolField::Ghost::NM - SSurfXField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (SVolField::Ghost::NM - SSurfXField::Ghost::NM);
  }
  template<> inline void shift_dest_index<SVolField,SSurfYField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (SVolField::Ghost::NM - SSurfYField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (SVolField::Ghost::NM - SSurfYField::Ghost::NM);
  }
  template<> inline void shift_dest_index<SVolField,SSurfZField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (SVolField::Ghost::NM - SSurfZField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (SVolField::Ghost::NM - SSurfZField::Ghost::NM);
  }

  template<> inline void shift_dest_index<SVolField,XVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  { --triplet.i; }

  template<> inline void shift_dest_index<SVolField,YVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  { --triplet.j; }

  template<> inline void shift_dest_index<SVolField,ZVolField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  { --triplet.k; }

  template<> inline void shift_dest_index<XVolField,XSurfXField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[1]>1 ) triplet.j += (XVolField::Ghost::NM - XSurfXField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (XVolField::Ghost::NM - XSurfXField::Ghost::NM);
  }
  template<> inline void shift_dest_index<XVolField,XSurfYField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (XVolField::Ghost::NM - XSurfYField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (XVolField::Ghost::NM - XSurfYField::Ghost::NM);
  }
  template<> inline void shift_dest_index<XVolField,XSurfZField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (XVolField::Ghost::NM - XSurfYField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (XVolField::Ghost::NM - XSurfYField::Ghost::NM);
  }

  template<> inline void shift_dest_index<YVolField,YSurfXField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[1]>1 ) triplet.j += (YVolField::Ghost::NM - YSurfXField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (YVolField::Ghost::NM - YSurfXField::Ghost::NM);
  }
  template<> inline void shift_dest_index<YVolField,YSurfYField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (YVolField::Ghost::NM - YSurfYField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (YVolField::Ghost::NM - YSurfYField::Ghost::NM);
  }
  template<> inline void shift_dest_index<YVolField,YSurfZField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (YVolField::Ghost::NM - YSurfYField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (YVolField::Ghost::NM - YSurfYField::Ghost::NM);
  }

  template<> inline void shift_dest_index<ZVolField,ZSurfXField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[1]>1 ) triplet.j += (ZVolField::Ghost::NM - ZSurfXField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (ZVolField::Ghost::NM - ZSurfXField::Ghost::NM);
  }
  template<> inline void shift_dest_index<ZVolField,ZSurfYField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (ZVolField::Ghost::NM - ZSurfYField::Ghost::NM);
    if( dim[2]>1 ) triplet.k += (ZVolField::Ghost::NM - ZSurfYField::Ghost::NM);
  }
  template<> inline void shift_dest_index<ZVolField,ZSurfZField>( const std::vector<int>& dim, const int irow, IndexTriplet& triplet )
  {
    if( dim[0]>1 ) triplet.i += (ZVolField::Ghost::NM - ZSurfYField::Ghost::NM);
    if( dim[1]>1 ) triplet.j += (ZVolField::Ghost::NM - ZSurfYField::Ghost::NM);
  }


}// namespace FVStaggered
}// namespace SpatialOps

#endif
