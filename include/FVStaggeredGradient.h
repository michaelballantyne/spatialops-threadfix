#ifndef FVSS_GradientAssembler_h
#define FVSS_GradientAssembler_h

#include <FVStaggeredIndexHelper.h>

#include <SpatialField.h>
#include <SpatialOperator.h>


namespace SpatialOps{


  // forward declaration.
  namespace FVStaggered{ template<typename T1,typename T2> class GradientAssembler; }

  template< typename SrcField, typename DestField >
  struct OpAssemblerSelector< FVStaggered::Gradient, SrcField, DestField >
  {
    typedef FVStaggered::GradientAssembler<SrcField,DestField>  Assembler;
  };


namespace FVStaggered{


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
     */
    GradientAssembler( const double meshSpacing,
		       const std::vector<int>& dimExtent );

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
    void get_ghost_cols( std::set<int>& ghostCols ) const
    {
      ghostCols = get_ghost_set<SrcField>(dim_);
    }

    /**
     *  @brief Obtain the set of row indices corresponding to ghost entries.
     */
    void get_ghost_rows( std::set<int>& ghostRows ) const
    {
      ghostRows = get_ghost_set<DestField>(dim_);
    }

  private:

    const std::vector<int>& dim_;
    const IndexHelper<SrcField,DestField> indexHelper_;
    const double coef_;

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
		     const std::vector<int>& dimExtent )
    : dim_( dimExtent ),
      indexHelper_( dimExtent ),
      coef_( 1.0/meshSpacing )
  {
  }
  //------------------------------------------------------------------
  template< typename SrcField, typename DestField >
  int
  GradientAssembler<SrcField,DestField>::
  get_ncols() const
  {
    int n=1;
    if( get_n_tot<DestField>(dim_)>1 )
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
    switch( DestField::Location::Dir::value ){
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
    switch( DestField::Location::StagDir::value ){
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

}
}

#endif
