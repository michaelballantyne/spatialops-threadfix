#ifndef FVStaggeredInterpolant_h
#define FVStaggeredInterpolant_h

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/FVStaggeredIndexHelper.h>


namespace SpatialOps{


  // forward declaration.
  namespace structured{ template<typename T1,typename T2> class LinearInterpolantAssembler; }

  // this is required for the SpatialOperator class.  It specifies
  // that we should use the LinearInterpolantAssembler to construct
  // Interpolant objects.
  template< typename SrcFieldT, typename DestFieldT >
  struct OpAssemblerSelector< Interpolant, SrcFieldT, DestFieldT >
  {
    typedef structured::LinearInterpolantAssembler<SrcFieldT,DestFieldT>  Assembler;
  };


namespace structured{


  /**
   *  @class  LinearInterpolantAssembler
   *  @author James C. Sutherland
   *
   *  @brief Assembles linear interpolant operators for UNIFORM staggered meshes.
   *
   *  An assembler for a linear interpolant operator to be used to
   *  construct a SpatialOperator.  This conforms to the requirements
   *  for a SpatialOperator Assembler as defined in the documentation
   *  on the OpAssemblerSelector.
   *
   *  This implementation assumes a UNIFORM MESH spacing.  For
   *  nonuniform meshes, the sparsity pattern will remain unchanged,
   *  but the coefficient values will change.
   *
   *  @par Template Parameters
   *  <ul>
   *  <li> \b SrcFieldT The type of the source field for the interpolant operator
   *  <li> \b DestFieldT The type of the destination field for the interpolant operator
   *  <\ul>
   */
  template< typename SrcFieldT,
            typename DestFieldT >
  class LinearInterpolantAssembler
  {
  public:

    unsigned int num_nonzeros() const;

    /**
     *  @brief Construct a LinearInterpolantAssembler.
     *  @param dimExtent 3-component vector indicating the number of
     *  interior points in each coordinate direction.
     */
    LinearInterpolantAssembler( const IntVec& dimExtent,
                                const bool hasPlusXSideFaces = true,
                                const bool hasPlusYSideFaces = true,
                                const bool hasPlusZSideFaces = true );

    ~LinearInterpolantAssembler(){}

    int get_ncols() const;   ///< Return the number of columns in this operator
    int get_nrows() const;   ///< Return the number of rows in this operator

    /**
     *  @brief Obtain the nonzero values for this operator on the given row.
     *
     *  @param irow The index for the row to obtain entries for.
     *  @param vals A vector containing the values for the nonzero
     *         entries in this row.
     *  @param ixs A vector containing the column indices
     *         corresponding to each nonzero entry.
     */
    void get_row_entries( const int irow,
                          std::vector<double> & vals,
                          std::vector<int> & ixs ) const;

    /**
     *  @brief Obtain the set of column indices corresponding to ghost entries.
     */
    void get_ghost_cols( std::set<size_t>& ghostCols ) const
    {
      ghostCols = get_ghost_set<SrcFieldT>( dim_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    }

    /**
     *  @brief Obtain the set of row indices corresponding to ghost entries.
     */
    void get_ghost_rows( std::set<size_t>& ghostRows ) const
    {
      ghostRows = get_ghost_set<DestFieldT>( dim_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    }

  private:

    const IndexHelper<SrcFieldT,DestFieldT> indexHelper_;
    const IntVec dim_;
    const bool hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_;

  };

  //==================================================================






  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //                           Implementation
  //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






  //==================================================================


  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  LinearInterpolantAssembler<SrcFieldT,DestFieldT>::
  LinearInterpolantAssembler( const IntVec& dimExtent,
                              const bool hasPlusXSideFaces,
                              const bool hasPlusYSideFaces,
                              const bool hasPlusZSideFaces )
    : indexHelper_( dimExtent, hasPlusXSideFaces, hasPlusYSideFaces, hasPlusZSideFaces ),
      dim_( dimExtent ),
      hasPlusXSideFaces_( hasPlusXSideFaces ),
      hasPlusYSideFaces_( hasPlusYSideFaces ),
      hasPlusZSideFaces_( hasPlusZSideFaces )
  {
  }
  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  unsigned int
  LinearInterpolantAssembler<SrcFieldT,DestFieldT>::
  num_nonzeros() const
  {
    return 2;
  }
  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  int
  LinearInterpolantAssembler<SrcFieldT,DestFieldT>::
  get_ncols() const
  {
    return indexHelper_.get_ncol();
  }
  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  int
  LinearInterpolantAssembler<SrcFieldT,DestFieldT>::
  get_nrows() const
  {
    return indexHelper_.get_nrow();
  }
  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  void
  LinearInterpolantAssembler<SrcFieldT,DestFieldT>::
  get_row_entries( const int irow,
                   std::vector<double> & vals,
                   std::vector<int> & ixs ) const
  {
    switch( int(DestFieldT::Location::FaceDir::value) ){
    case XDIR::value:
      if( dim_[0]==1 ) return;
      break;
    case YDIR::value:
      if( dim_[1]==1 ) return;
      break;
    case ZDIR::value:
      if( dim_[2]==1 ) return;
      break;
    }
    switch( int(DestFieldT::Location::StagLoc::value) ){
    case XDIR::value:
      if( dim_[0]==1 ) return;
      break;
    case YDIR::value:
      if( dim_[1]==1 ) return;
      break;
    case ZDIR::value:
      if( dim_[2]==1 ) return;
      break;
    }
    indexHelper_.get_cols( irow, ixs );
    if( !ixs.empty() ){
      assert( ixs.size() == num_nonzeros() );
      const double val = 1.0/double(num_nonzeros());
      vals.assign( num_nonzeros(), val );
    }
  }

  //------------------------------------------------------------------
  // specializations for XVolField->SSurfXField
  template<> inline unsigned int LinearInterpolantAssembler<XVolField,SSurfXField>::num_nonzeros() const{ return 1; }

  template<> inline void LinearInterpolantAssembler<XVolField,SSurfXField>::get_row_entries( const int irow, std::vector<double> & vals, std::vector<int> & ixs ) const
  {
    if( dim_[0]==1 ) return;

    const int ghostDiff = XVolField::Ghost::NGHOST - SSurfXField::Ghost::NGHOST;
    int icol = irow + ghostDiff;
    if( dim_[1]>1 ){
      icol += ghostDiff * get_nx_with_ghost<XVolField>(dim_[0],hasPlusXSideFaces_);
    }
    if( dim_[2]>1 ){
      icol += (XVolField::Ghost::NGHOST-SSurfXField::Ghost::NGHOST)
        * get_nx_with_ghost<XVolField>(dim_[0],hasPlusXSideFaces_)
        * get_ny_with_ghost<XVolField>(dim_[1],hasPlusYSideFaces_);
    }
    ixs.push_back(icol);
    vals.push_back(1.0);
  }

  //------------------------------------------------------------------
  // specializations for YVolField->SSurfYField
  template<> inline unsigned int LinearInterpolantAssembler<YVolField,SSurfYField>::num_nonzeros() const{ return 1; }

  template<> inline void LinearInterpolantAssembler<YVolField,SSurfYField>::get_row_entries( const int irow, std::vector<double> & vals, std::vector<int> & ixs ) const
  {
    if( dim_[1]==1 ) return;

    const int ghostDiff = YVolField::Ghost::NGHOST - SSurfYField::Ghost::NGHOST;
    int icol = irow + ghostDiff;
    if( dim_[0]>1 ){
      icol += ghostDiff
        * get_nx_with_ghost<YVolField>(dim_[0],hasPlusXSideFaces_);
    }
    if( dim_[2]>1 ){
      icol += ghostDiff
        * get_nx_with_ghost<YVolField>(dim_[0],hasPlusXSideFaces_)
        * get_ny_with_ghost<YVolField>(dim_[1],hasPlusYSideFaces_);
    }
    ixs.push_back(icol);
    vals.push_back(1.0);
  }

  //------------------------------------------------------------------
  // specializations for ZVolField->SSurfZField
  template<> inline unsigned int LinearInterpolantAssembler<ZVolField,SSurfZField>::num_nonzeros() const{ return 1; }

  template<> inline void LinearInterpolantAssembler<ZVolField,SSurfZField>::get_row_entries( const int irow, std::vector<double> & vals, std::vector<int> & ixs ) const
  {
    if( dim_[2]==1 ) return;

    const int ghostDiff = ZVolField::Ghost::NGHOST - SSurfZField::Ghost::NGHOST;
    int icol = irow + ghostDiff;
    if( dim_[0]>1 ){
      icol += ghostDiff * get_nx_with_ghost<ZVolField>(dim_[0],hasPlusXSideFaces_);
    }
    if( dim_[1]>1 ){
      icol += ghostDiff
        * get_nx_with_ghost<ZVolField>(dim_[0],hasPlusXSideFaces_)
        * get_ny_with_ghost<ZVolField>(dim_[1],hasPlusYSideFaces_);
    }
    ixs.push_back(icol);
    vals.push_back(1.0);
  }

  //--------------------------------------------------------------------
  // specializations for SVolField -> XSurf fields
  template<> inline unsigned int LinearInterpolantAssembler<SVolField,XSurfXField>::num_nonzeros() const{ return 1; }
  template<> inline unsigned int LinearInterpolantAssembler<SVolField,XSurfYField>::num_nonzeros() const{ return 4; }
  template<> inline unsigned int LinearInterpolantAssembler<SVolField,XSurfZField>::num_nonzeros() const{ return 4; }

  //------------------------------------------------------------------
  // specializations for SVolField -> YSurf fields
  template<> inline unsigned int LinearInterpolantAssembler<SVolField,YSurfXField>::num_nonzeros() const{ return 4; }
  template<> inline unsigned int LinearInterpolantAssembler<SVolField,YSurfYField>::num_nonzeros() const{ return 1; }
  template<> inline unsigned int LinearInterpolantAssembler<SVolField,YSurfZField>::num_nonzeros() const{ return 4; }

  //------------------------------------------------------------------
  // specializations for SVolField -> ZSurf fields
  template<> inline unsigned int LinearInterpolantAssembler<SVolField,ZSurfXField>::num_nonzeros() const{ return 4; }
  template<> inline unsigned int LinearInterpolantAssembler<SVolField,ZSurfYField>::num_nonzeros() const{ return 4; }
  template<> inline unsigned int LinearInterpolantAssembler<SVolField,ZSurfZField>::num_nonzeros() const{ return 1; }

  //==================================================================


} // namespace structured
} // namespace SpatialOps

#endif
