#ifndef FVSS_GradientAssembler_h
#define FVSS_GradientAssembler_h

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/matrix/FVStaggeredIndexHelper.h>


namespace SpatialOps{


  // forward declaration.
  namespace structured{ template<typename T1,typename T2> class GradientAssembler; }

  template< typename SrcField, typename DestField >
  struct OpAssemblerSelector< Gradient, SrcField, DestField >
  {
    typedef structured::GradientAssembler<SrcField,DestField>  Assembler;
  };


namespace structured{


  //==================================================================


  /**
   *  @class  GradientAssembler
   *  @author James C. Sutherland
   *
   *  @brief Assembles gradient operators for UNIFORM staggered meshes.
   *
   *  An assembler for a gradient operator to be used to construct a
   *  SpatialOperator.  This conforms to the requirements for a
   *  SpatialOperator Assembler as defined in the documentation on the
   *  OpAssemblerSelector.
   *
   *  This implementation assumes a UNIFORM MESH spacing.  For
   *  nonuniform meshes, the sparsity pattern will remain unchanged,
   *  but the coefficient values will change.
   *
   *  @par Template Parameters
   *    \li \b SrcField The type for the source field.
   *    \li \b DestField The type for the destination field.
   *
   */
  template< typename SrcField,
            typename DestField >
  class GradientAssembler
  {
  public:

    /** @brief Return the number of nonzero entries for this operator. */
    unsigned int num_nonzeros() const{ return 2; }

    /**
     *  @brief Construct a GradientAssembler object.
     *
     *  @param meshSpacing The grid spacing (assumed uniform) in the
     *  direction that this gradient operator applies to.
     *
     *  @param dimExtent A vector with three elements indicating the
     *  domain extent (number of cells) in each coordinate direction.
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
    GradientAssembler( const double meshSpacing,
                       const IntVec dimExtent,
                       const bool hasPlusXSideFaces,
                       const bool hasPlusYSideFaces,
                       const bool hasPlusZSideFaces );

    ~GradientAssembler(){}

    int get_ncols() const;
    int get_nrows() const;

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
      ghostCols = get_ghost_set<SrcField>( dim_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    }

    /**
     *  @brief Obtain the set of row indices corresponding to ghost entries.
     */
    void get_ghost_rows( std::set<size_t>& ghostRows ) const
    {
      ghostRows = get_ghost_set<DestField>( dim_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    }

  private:

    const IntVec dim_;
    const IndexHelper<SrcField,DestField> indexHelper_;
    const double coef_;
    const bool hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_;
  };


  //==================================================================




  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //  Implementations
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




  //==================================================================


  //------------------------------------------------------------------
  template< typename SrcField, typename DestField >
  GradientAssembler<SrcField,DestField>::
  GradientAssembler( const double meshSpacing,
                     const IntVec dimExtent,
                     const bool hasPlusXSideFaces,
                     const bool hasPlusYSideFaces,
                     const bool hasPlusZSideFaces )
    : dim_( dimExtent ),
      indexHelper_( dimExtent, hasPlusXSideFaces, hasPlusYSideFaces, hasPlusZSideFaces ),
      coef_( 1.0/meshSpacing ),
      hasPlusXSideFaces_( hasPlusXSideFaces ),
      hasPlusYSideFaces_( hasPlusYSideFaces ),
      hasPlusZSideFaces_( hasPlusZSideFaces )
  {
  }
  //------------------------------------------------------------------
  template< typename SrcField, typename DestField >
  int
  GradientAssembler<SrcField,DestField>::
  get_ncols() const
  {
    int n=1;
    if( get_ntot_with_ghost<DestField>( dim_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ ) > 1 )
      n = indexHelper_.get_ncol();
    return n;
  }
  //------------------------------------------------------------------
  template< typename SrcField, typename DestField >
  int
  GradientAssembler<SrcField,DestField>::
  get_nrows() const
  {
    return indexHelper_.get_nrow();
  }
  //------------------------------------------------------------------
  template< typename SrcField, typename DestField >
  void
  GradientAssembler<SrcField,DestField>::
  get_row_entries( const int irow,
                   std::vector<double> & vals,
                   std::vector<int> & ixs ) const
  {
    switch( int(DestField::Location::FaceDir::value) ){
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
    switch( int(DestField::Location::StagLoc::value) ){
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
    if( ixs.size() == 2 ){
      vals.push_back( -coef_ );
      vals.push_back(  coef_ );
    }
  }
  //------------------------------------------------------------------

} // namespace structured
} // namespace SpatialOps

#endif
