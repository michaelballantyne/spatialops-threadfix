#ifndef FVStaggeredIndexHelper_h
#define FVStaggeredIndexHelper_h

#include <spatialops/SpatialOpsConfigure.h>

#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/FVStaggeredTypes.h>


namespace SpatialOps{
namespace structured{


  //==================================================================

  /**
   *  @class IndexHelper
   *  @author James C. Sutherland
   *
   *  Note that this class will only work properly with the field
   *  types defined in FVStaggeredTypes.h.
   *
   *  @todo Check on all implementations of shift_dest_to_src for consistency w.r.t. ghosting.
   *  @todo Document all of this so it is comprehensible.
   */
  template< typename SrcFieldT, typename DestFieldT >
  class IndexHelper
  {
  public:

    /**
     *  @param dim A vector with three elements indicating the domain
     *  extent (number of cells) in each coordinate direction.
     *
     *  @param hasPlusXSideFaces Boolean flag to indicate if the
     *  operator is to be constructed including face cells on the +X
     *  side of the domain.
     *
     *  @param hasPlusYSideFaces Boolean flag to indicate if the
     *  operator is to be constructed including face cells on the +Y
     *  side of the domain.
     *
     *  @param hasPlusZSideFaces Boolean flag to indicate if the
     *  operator is to be constructed including face cells on the +Z
     *  side of the domain.
     */
    IndexHelper( const IntVec& dim,
                 const bool hasPlusXSideFaces,
                 const bool hasPlusYSideFaces,
                 const bool hasPlusZSideFaces );
    void get_cols( const int irow, std::vector<int>& cols ) const;
    int get_ncol() const;
    int get_nrow() const;
    /**
     *  \brief Calculates the stride in a given direction based on the source
     *  and destination fields of the IndexHelper instance. One can determine
     *  the face direction from the destination field type and based on that
     *  the function calculates how many cells to stride.
     */    
    int calculate_stride() const;


  private:
    const IntVec& dim_;
    const bool hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_;



    template<typename T1, typename T2, int IsSurf> struct DirSelector{};
    template<typename T1, typename T2> struct DirSelector<T1,T2,1>{ typedef T1 Type; };
    template<typename T1, typename T2> struct DirSelector<T1,T2,0>{ typedef T2 Type; };
  };

  //==================================================================

  inline bool is_valid_entry( const IntVec& dim,
                              const IntVec& ixt1,
                              const IntVec& ixt2 );

  //==================================================================

  /**
   *  @brief Determines if the row and column correspond to interior
   *  points for the given source/destination fields.
   *
   *  @param ixs Flat index for the source field
   *
   *  @param ixd Flat index for the destination field.
   *
   *  @param dim A vector with three elements indicating the domain
   *  extent (number of cells) in each coordinate direction.
   *
   *  @param hasPlusXSideFaces Boolean flag to indicate if the
   *  operator is to be constructed including face cells on the +X
   *  side of the domain.
   *
   *  @param hasPlusYSideFaces Boolean flag to indicate if the
   *  operator is to be constructed including face cells on the +Y
   *  side of the domain.
   *
   *  @param hasPlusZSideFaces Boolean flag to indicate if the
   *  operator is to be constructed including face cells on the +Z
   *  side of the domain.
   */
  template< typename SrcFieldT,
            typename DestFieldT >
  bool is_in_bounds( const IntVec& dim,
                     const int ixs,
                     const int ixd,
                     const bool hasPlusXSideFaces,
                     const bool hasPlusYSideFaces,
                     const bool hasPlusZSideFaces );

  //==================================================================






  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //                          Implementation
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  /**
   *  @brief Given the IntVec for the destination, convert it
   *  to an IntVec for the source field.
   *
   *  @param dim A vector with three elements indicating the domain
   *  extent (number of cells) in each coordinate direction.
   *
   *  @param triplet The ijk location. On entry, this contains the
   *  location for the destination field. On exit it contains the ijk
   *  location for the source field.
   */
  template< typename SrcT, typename DestT >
  void shift_dest_to_src( const IntVec& dim,
                          IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (DestT::Ghost::NGHOST-SrcT::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (DestT::Ghost::NGHOST-SrcT::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (DestT::Ghost::NGHOST-SrcT::Ghost::NGHOST);
  }

  // Scalar volume to surface operations (e.g. gradient, interpolant)
  template<> inline void shift_dest_to_src<SVolField,SSurfXField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SSurfXField::Ghost::NGHOST-SVolField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (SSurfXField::Ghost::NGHOST-SVolField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (SSurfXField::Ghost::NGHOST-SVolField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<SVolField,SSurfYField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SSurfYField::Ghost::NGHOST-SVolField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (SSurfYField::Ghost::NGHOST-SVolField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (SSurfYField::Ghost::NGHOST-SVolField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<SVolField,SSurfZField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SSurfZField::Ghost::NGHOST-SVolField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (SSurfZField::Ghost::NGHOST-SVolField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (SSurfZField::Ghost::NGHOST-SVolField::Ghost::NGHOST)-1;
  }


  // X volume to surface operations (e.g. gradient, interpolant)
  template<> inline void shift_dest_to_src<XVolField,XSurfXField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (XSurfXField::Ghost::NGHOST-XVolField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (XSurfXField::Ghost::NGHOST-XVolField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (XSurfXField::Ghost::NGHOST-XVolField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<XVolField,XSurfYField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (XSurfYField::Ghost::NGHOST-XVolField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (XSurfYField::Ghost::NGHOST-XVolField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (XSurfYField::Ghost::NGHOST-XVolField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<XVolField,XSurfZField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (XSurfZField::Ghost::NGHOST-XVolField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (XSurfZField::Ghost::NGHOST-XVolField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (XSurfZField::Ghost::NGHOST-XVolField::Ghost::NGHOST)-1;
  }


  // Y volume to surface operations (e.g. gradient, interpolant)
  template<> inline void shift_dest_to_src<YVolField,YSurfXField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (YSurfXField::Ghost::NGHOST-YVolField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (YSurfXField::Ghost::NGHOST-YVolField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (YSurfXField::Ghost::NGHOST-YVolField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<YVolField,YSurfYField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (YSurfYField::Ghost::NGHOST-YVolField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (YSurfYField::Ghost::NGHOST-YVolField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (YSurfYField::Ghost::NGHOST-YVolField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<YVolField,YSurfZField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (YSurfZField::Ghost::NGHOST-YVolField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (YSurfZField::Ghost::NGHOST-YVolField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (YSurfZField::Ghost::NGHOST-YVolField::Ghost::NGHOST)-1;
  }


  // Z volume to surface operations (e.g. gradient, interpolant)
  template<> inline void shift_dest_to_src<ZVolField,ZSurfXField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (ZSurfXField::Ghost::NGHOST-ZVolField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (ZSurfXField::Ghost::NGHOST-ZVolField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (ZSurfXField::Ghost::NGHOST-ZVolField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<ZVolField,ZSurfYField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (ZSurfYField::Ghost::NGHOST-ZVolField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (ZSurfYField::Ghost::NGHOST-ZVolField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (ZSurfYField::Ghost::NGHOST-ZVolField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<ZVolField,ZSurfZField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (ZSurfZField::Ghost::NGHOST-ZVolField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (ZSurfZField::Ghost::NGHOST-ZVolField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (ZSurfZField::Ghost::NGHOST-ZVolField::Ghost::NGHOST)-1;
  }


  // Scalar volume to staggered volume
  template<> inline void shift_dest_to_src<SVolField,XVolField>( const IntVec& dim, IntVec& triplet ){ shift_dest_to_src<SVolField,SSurfXField>( dim, triplet ); }
  template<> inline void shift_dest_to_src<SVolField,YVolField>( const IntVec& dim, IntVec& triplet ){ shift_dest_to_src<SVolField,SSurfYField>( dim, triplet ); }
  template<> inline void shift_dest_to_src<SVolField,ZVolField>( const IntVec& dim, IntVec& triplet ){ shift_dest_to_src<SVolField,SSurfZField>( dim, triplet ); }

  // Scalar volume to staggered surfaces
  template<> inline void shift_dest_to_src<SVolField,XSurfXField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SVolField::Ghost::NGHOST-XSurfXField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (SVolField::Ghost::NGHOST-XSurfXField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (SVolField::Ghost::NGHOST-XSurfXField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<SVolField,XSurfYField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SVolField::Ghost::NGHOST-XSurfYField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (SVolField::Ghost::NGHOST-XSurfYField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (SVolField::Ghost::NGHOST-XSurfYField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<SVolField,XSurfZField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SVolField::Ghost::NGHOST-XSurfZField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (SVolField::Ghost::NGHOST-XSurfZField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (SVolField::Ghost::NGHOST-XSurfZField::Ghost::NGHOST)-1;
  }
  template<> inline void shift_dest_to_src<SVolField,YSurfXField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SVolField::Ghost::NGHOST-YSurfXField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (SVolField::Ghost::NGHOST-YSurfXField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (SVolField::Ghost::NGHOST-YSurfXField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<SVolField,YSurfYField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SVolField::Ghost::NGHOST-YSurfYField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (SVolField::Ghost::NGHOST-YSurfYField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (SVolField::Ghost::NGHOST-YSurfYField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<SVolField,YSurfZField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SVolField::Ghost::NGHOST-YSurfZField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (SVolField::Ghost::NGHOST-YSurfZField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (SVolField::Ghost::NGHOST-YSurfZField::Ghost::NGHOST)-1;
  }
  template<> inline void shift_dest_to_src<SVolField,ZSurfXField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SVolField::Ghost::NGHOST-ZSurfXField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (SVolField::Ghost::NGHOST-ZSurfXField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (SVolField::Ghost::NGHOST-ZSurfXField::Ghost::NGHOST)-1;
  }
  template<> inline void shift_dest_to_src<SVolField,ZSurfYField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SVolField::Ghost::NGHOST-ZSurfYField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (SVolField::Ghost::NGHOST-ZSurfYField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (SVolField::Ghost::NGHOST-ZSurfYField::Ghost::NGHOST)-1;
  }
  template<> inline void shift_dest_to_src<SVolField,ZSurfZField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (SVolField::Ghost::NGHOST-ZSurfZField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (SVolField::Ghost::NGHOST-ZSurfZField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (SVolField::Ghost::NGHOST-ZSurfZField::Ghost::NGHOST)-1;
  }


  // jcs ....
  // x-staggered to y-staggered
  template<> inline void shift_dest_to_src<XVolField,YSurfXField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (XVolField::Ghost::NGHOST - YSurfXField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (XVolField::Ghost::NGHOST - YSurfXField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (XVolField::Ghost::NGHOST - YSurfXField::Ghost::NGHOST);
  }

  template<> inline void shift_dest_to_src<XVolField,ZSurfXField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (XVolField::Ghost::NGHOST - ZSurfXField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (XVolField::Ghost::NGHOST - ZSurfXField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (XVolField::Ghost::NGHOST - ZSurfXField::Ghost::NGHOST)-1;
  }

  template<> inline void shift_dest_to_src<YVolField,XSurfYField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (YVolField::Ghost::NGHOST - XSurfYField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (YVolField::Ghost::NGHOST - XSurfYField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (YVolField::Ghost::NGHOST - XSurfYField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<YVolField,ZSurfYField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (YVolField::Ghost::NGHOST - ZSurfYField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (YVolField::Ghost::NGHOST - ZSurfYField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (YVolField::Ghost::NGHOST - ZSurfYField::Ghost::NGHOST)-1;
  }

  template<> inline void shift_dest_to_src<ZVolField,XSurfZField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (ZVolField::Ghost::NGHOST - XSurfZField::Ghost::NGHOST)-1;
    if( dim[1]>1 ) triplet[1] += (ZVolField::Ghost::NGHOST - XSurfZField::Ghost::NGHOST);
    if( dim[2]>1 ) triplet[2] += (ZVolField::Ghost::NGHOST - XSurfZField::Ghost::NGHOST);
  }
  template<> inline void shift_dest_to_src<ZVolField,YSurfZField>( const IntVec& dim, IntVec& triplet )
  {
    if( dim[0]>1 ) triplet[0] += (ZVolField::Ghost::NGHOST - YSurfZField::Ghost::NGHOST);
    if( dim[1]>1 ) triplet[1] += (ZVolField::Ghost::NGHOST - YSurfZField::Ghost::NGHOST)-1;
    if( dim[2]>1 ) triplet[2] += (ZVolField::Ghost::NGHOST - YSurfZField::Ghost::NGHOST);
  }

  //------------------------------------------------------------------


  //==================================================================


  //------------------------------------------------------------------ 

  template<typename SrcFieldT, typename DestFieldT>
  IndexHelper<SrcFieldT,DestFieldT>::
  IndexHelper( const IntVec& dim,
               const bool hasPlusXSideFaces,
               const bool hasPlusYSideFaces,
               const bool hasPlusZSideFaces )
    : dim_( dim ),
      hasPlusXSideFaces_( hasPlusXSideFaces ),
      hasPlusYSideFaces_( hasPlusYSideFaces ),
      hasPlusZSideFaces_( hasPlusZSideFaces )
  {}
  //------------------------------------------------------------------ 
  template<typename SrcFieldT, typename DestFieldT>
  void
  IndexHelper<SrcFieldT,DestFieldT>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    // from the row (flat index of destination field), determine the
    // ijk indices for the destination field
    IntVec ijk = flat2ijk<DestFieldT>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    
    shift_dest_to_src<SrcFieldT,DestFieldT>( dim_, ijk );

    if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 ) return;

    // from the ijk index, determine the flat index for the source field.
    const int icol = ijk2flat<SrcFieldT>::value( dim_, ijk, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );

    const int stride = calculate_stride();

    if( is_in_bounds<SrcFieldT,DestFieldT>(dim_,icol,       irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SrcFieldT,DestFieldT>(dim_,icol+stride,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) ){

      const IntVec ijk2 = flat2ijk<SrcFieldT>::value( dim_, icol+stride, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );

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
    return get_ntot_with_ghost<SrcFieldT>(dim_,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_);
  }
  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  int
  IndexHelper<SrcFieldT,DestFieldT>::
  get_nrow() const
  {
    return get_ntot_with_ghost<DestFieldT>(dim_,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_);
  }
  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  int
  IndexHelper<SrcFieldT,DestFieldT>::
  calculate_stride() const
  {
    // jcs this is not robust at all...
    typedef typename DirSelector< SrcFieldT, DestFieldT, IsSameType<typename DestFieldT::Location::FaceDir,NODIR>::result >::Type FieldDirSelect;

    int n=-1;
    switch( FieldDirSelect::Location::FaceDir::value ){
    case XDIR::value:
      n=1;
      break;
    case YDIR::value:
      n = get_nx_with_ghost<SrcFieldT>(dim_[0],hasPlusXSideFaces_);
      break;
    case ZDIR::value:
      n = get_nx_with_ghost<SrcFieldT>(dim_[0],hasPlusXSideFaces_) * get_ny_with_ghost<SrcFieldT>(dim_[1],hasPlusYSideFaces_);
      break;
    default:
      std::cout << "ERROR: field location enum value=" << FieldDirSelect::Location::FaceDir::value << std::endl;
      assert(0);
    }
    return n;
  }
  template<> inline int IndexHelper<SVolField,XVolField>::calculate_stride() const
  {
    return 1;
  }
  template<> inline int IndexHelper<SVolField,YVolField>::calculate_stride() const
  {
    return get_nx_with_ghost<SVolField>(dim_[0],hasPlusXSideFaces_);
  }
  template<> inline int IndexHelper<SVolField,ZVolField>::calculate_stride() const
  {
    return get_nx_with_ghost<SVolField>(dim_[0],hasPlusXSideFaces_) * get_ny_with_ghost<SVolField>(dim_[1],hasPlusYSideFaces_);
  }
  template<> inline int IndexHelper<XVolField,YSurfXField>::calculate_stride() const
  {
    return get_nx_with_ghost<XVolField>(dim_[0],hasPlusXSideFaces_);
  }
  template<> inline int IndexHelper<XVolField,ZSurfXField>::calculate_stride() const
  {
    return get_nx_with_ghost<XVolField>(dim_[0],hasPlusXSideFaces_) * get_ny_with_ghost<XVolField>(dim_[0],hasPlusYSideFaces_);
  }
  template<> inline int IndexHelper<YVolField,XSurfYField>::calculate_stride() const
  {
    return 1;
  }
  template<> inline int IndexHelper<YVolField,ZSurfYField>::calculate_stride() const
  {
    return get_nx_with_ghost<YVolField>(dim_[0],hasPlusXSideFaces_) * get_ny_with_ghost<YVolField>(dim_[1],hasPlusYSideFaces_);
  }
  template<> inline int IndexHelper<ZVolField,XSurfZField>::calculate_stride() const
  {
    return 1;
  }
  template<> inline int IndexHelper<ZVolField,YSurfZField>::calculate_stride() const
  {
    return get_nx_with_ghost<ZVolField>(dim_[0],hasPlusXSideFaces_);
  }
  template<> inline int IndexHelper<XVolField,SVolField>::calculate_stride() const
  {
    return 1;
  }
  template<> inline int IndexHelper<YVolField,SVolField>::calculate_stride() const
  {
    return get_nx_with_ghost<YVolField>(dim_[0],hasPlusXSideFaces_);
  }
  template<> inline int IndexHelper<ZVolField,SVolField>::calculate_stride() const
  {
    return get_nx_with_ghost<ZVolField>(dim_[0],hasPlusXSideFaces_) * get_ny_with_ghost<ZVolField>(dim_[1],hasPlusYSideFaces_);
  }
  //------------------------------------------------------------------

  //==================================================================

  bool is_valid_entry( const IntVec& dim,
                       const IntVec& ixt1,
                       const IntVec& ixt2 )
  {
    if( std::abs(ixt2[0] - ixt1[0]) > 1 ) return false;
    if( std::abs(ixt2[1] - ixt1[1]) > 1 ) return false;
    if( std::abs(ixt2[2] - ixt1[2]) > 1 ) return false;
    return true;
  }
  
  //==================================================================

  template<typename SrcFieldT, typename DestFieldT>
  bool is_in_bounds( const IntVec& dim, const int ixs, const int ixd,
                     const bool hasPlusXSideFaces, const bool hasPlusYSideFaces, const bool hasPlusZSideFaces )
  {
    const IntVec it = flat2ijk<SrcFieldT>::value(dim,ixs,hasPlusXSideFaces,hasPlusYSideFaces,hasPlusZSideFaces);

    if( it[0]<0 || it[1]<0 || it[2]<0 ) return false;

    const int nx = get_nx_with_ghost<SrcFieldT>(dim[0],hasPlusXSideFaces);
    const int ny = get_ny_with_ghost<SrcFieldT>(dim[1],hasPlusYSideFaces);
    const int nz = get_nz_with_ghost<SrcFieldT>(dim[2],hasPlusZSideFaces);

    if( it[0]<nx && it[1]<ny && it[2]<nz )  return true;

    return false;
  }

  //==================================================================


  //------------------------------------------------------------------

  template<> inline void
  IndexHelper<SVolField,XSurfXField>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    IntVec ijk = flat2ijk<XSurfXField>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    shift_dest_to_src<SVolField,XSurfXField>( dim_, ijk );
    if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 ) return;
    const int icol = ijk2flat<SVolField>::value( dim_, ijk, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    if( is_in_bounds<SVolField,XSurfXField>(dim_,icol,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) ){
      cols.push_back( icol );
    }
  }

  template<> inline void
  IndexHelper<SVolField,XSurfYField>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    IntVec ijk = flat2ijk<XSurfYField>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    shift_dest_to_src<SVolField,XSurfYField>( dim_, ijk );
    if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 ) return;

    const int icol = ijk2flat<SVolField>::value( dim_, ijk, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    const int icol2 = icol + 1;
    const int icol3 = icol + get_nx_with_ghost<SVolField>(dim_[0],hasPlusXSideFaces_);
    const int icol4 = icol3 + 1;

    if( is_in_bounds<SVolField,XSurfYField>(dim_,icol ,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,XSurfYField>(dim_,icol2,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,XSurfYField>(dim_,icol3,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,XSurfYField>(dim_,icol4,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) ){
      cols.push_back( icol  );
      cols.push_back( icol2 );
      cols.push_back( icol3 );
      cols.push_back( icol4 );
    }
  }

  template<> inline void
  IndexHelper<SVolField,XSurfZField>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    IntVec ijk = flat2ijk<XSurfZField>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    shift_dest_to_src<SVolField,XSurfZField>( dim_, ijk );
    if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 ) return;

    const int icol = ijk2flat<SVolField>::value( dim_, ijk, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    const int icol2 = icol + 1;
    const int icol3 = icol + get_nx_with_ghost<SVolField>(dim_[0],hasPlusXSideFaces_)*get_ny_with_ghost<SVolField>(dim_[1],hasPlusYSideFaces_);
    const int icol4 = icol3 + 1;

    if( is_in_bounds<SVolField,XSurfZField>(dim_,icol ,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,XSurfZField>(dim_,icol2,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,XSurfZField>(dim_,icol3,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,XSurfZField>(dim_,icol4,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) ){
      cols.push_back( icol  );
      cols.push_back( icol2 );
      cols.push_back( icol3 );
      cols.push_back( icol4 );
    }
  }


  template<> inline void
  IndexHelper<SVolField,YSurfXField>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    IntVec ijk = flat2ijk<YSurfXField>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    shift_dest_to_src<SVolField,YSurfXField>( dim_, ijk );
    if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 ) return;

    const int icol = ijk2flat<SVolField>::value( dim_, ijk, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    const int icol2 = icol + 1;
    const int icol3 = icol + get_nx_with_ghost<SVolField>(dim_[0],hasPlusXSideFaces_);
    const int icol4 = icol3 + 1;

    if( is_in_bounds<SVolField,YSurfXField>(dim_,icol ,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,YSurfXField>(dim_,icol2,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,YSurfXField>(dim_,icol3,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,YSurfXField>(dim_,icol4,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) ){
      cols.push_back( icol  );
      cols.push_back( icol2 );
      cols.push_back( icol3 );
      cols.push_back( icol4 );
    }
  }

  template<> inline void
  IndexHelper<SVolField,YSurfYField>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    IntVec ijk = flat2ijk<YSurfYField>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    shift_dest_to_src<SVolField,YSurfYField>( dim_, ijk );
    if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 ) return;
    const int icol = ijk2flat<SVolField>::value( dim_, ijk, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    if( is_in_bounds<SVolField,YSurfYField>(dim_,icol,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) ){
      cols.push_back( icol );
    }
  }

  template<> inline void
  IndexHelper<SVolField,YSurfZField>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    IntVec ijk = flat2ijk<YSurfZField>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    shift_dest_to_src<SVolField,YSurfZField>( dim_, ijk );
    if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 ) return;

    const int nx = get_nx_with_ghost<SVolField>(dim_[0],hasPlusXSideFaces_);
    const int ny = get_ny_with_ghost<SVolField>(dim_[1],hasPlusYSideFaces_);

    const int icol = ijk2flat<SVolField>::value( dim_, ijk, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    const int icol2 = icol + nx;
    const int icol3 = icol + nx*ny;
    const int icol4 = icol3 + nx;

    if( is_in_bounds<SVolField,YSurfZField>(dim_,icol ,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,YSurfZField>(dim_,icol2,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,YSurfZField>(dim_,icol3,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,YSurfZField>(dim_,icol4,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) ){
      cols.push_back( icol  );
      cols.push_back( icol2 );
      cols.push_back( icol3 );
      cols.push_back( icol4 );
    }
  }


  template<> inline void
  IndexHelper<SVolField,ZSurfXField>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    IntVec ijk = flat2ijk<ZSurfXField>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    shift_dest_to_src<SVolField,ZSurfXField>( dim_, ijk );
    if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 ) return;

    const int icol = ijk2flat<SVolField>::value( dim_, ijk, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    const int icol2 = icol + 1;
    const int icol3 = icol + get_nx_with_ghost<SVolField>(dim_[0],hasPlusXSideFaces_)*get_ny_with_ghost<SVolField>(dim_[1],hasPlusYSideFaces_);
    const int icol4 = icol3 + 1;

    if( is_in_bounds<SVolField,ZSurfXField>(dim_,icol ,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,ZSurfXField>(dim_,icol2,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,ZSurfXField>(dim_,icol3,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,ZSurfXField>(dim_,icol4,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) ){
      cols.push_back( icol  );
      cols.push_back( icol2 );
      cols.push_back( icol3 );
      cols.push_back( icol4 );
    }
  }

  template<> inline void
  IndexHelper<SVolField,ZSurfYField>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    IntVec ijk = flat2ijk<ZSurfYField>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    shift_dest_to_src<SVolField,ZSurfYField>( dim_, ijk );
    if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 ) return;

    const int nx = get_nx_with_ghost<SVolField>(dim_[0],hasPlusXSideFaces_);
    const int ny = get_ny_with_ghost<SVolField>(dim_[1],hasPlusYSideFaces_);

    const int icol = ijk2flat<SVolField>::value( dim_, ijk, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    const int icol2 = icol + nx;
    const int icol3 = icol + nx*ny;
    const int icol4 = icol3 + nx;

    if( is_in_bounds<SVolField,ZSurfYField>(dim_,icol ,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,ZSurfYField>(dim_,icol2,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,ZSurfYField>(dim_,icol3,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) &&
        is_in_bounds<SVolField,ZSurfYField>(dim_,icol4,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) ){
      cols.push_back( icol  );
      cols.push_back( icol2 );
      cols.push_back( icol3 );
      cols.push_back( icol4 );
    }
  }

  template<> inline void
  IndexHelper<SVolField,ZSurfZField>::
  get_cols( const int irow, std::vector<int>& cols ) const
  {
    IntVec ijk = flat2ijk<ZSurfZField>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    shift_dest_to_src<SVolField,ZSurfZField>( dim_, ijk );
    if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 ) return;
    const int icol = ijk2flat<SVolField>::value( dim_, ijk, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    if( is_in_bounds<SVolField,ZSurfZField>(dim_,icol,irow,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) ){
      cols.push_back( icol );
    }
  }

  //------------------------------------------------------------------

}// namespace structured
}// namespace SpatialOps

#endif
