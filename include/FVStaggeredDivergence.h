#ifndef FVSS_DivergenceAssembler_h
#define FVSS_DivergenceAssembler_h

#include <FVStaggeredIndexHelper.h>

#include <SpatialField.h>
#include <SpatialOperator.h>


namespace SpatialOps{

  // forward declaration.
  namespace FVStaggered{ template<typename T1,typename T2> class DivergenceAssembler; }

  template< typename SrcField, typename DestField >
  struct OpAssemblerSelector< FVStaggered::Divergence, SrcField, DestField >
  {
    typedef FVStaggered::DivergenceAssembler<SrcField,DestField>  Assembler;
  };


namespace FVStaggered{

  /**
   *  @class  DivergenceAssembler
   *  @author James C. Sutherland
   *
   *  @brief Assembles divergence operators for UNIFORM staggered meshes.
   *
   *  An assembler for a divergence operator to be used to construct a
   *  SpatialOperator.  This conforms to the requirements for a
   *  SpatialOperator Assembler as defined in the documentation on the
   *  OpAssemblerSelector.
   *
   *  This implementation assumes a UNIFORM MESH spacing.  For
   *  nonuniform meshes, the sparsity pattern will remain unchanged,
   *  but the coefficient values will change.
   *
   *  @par Template Parameters
   *    \li \b SrcField The policy for the source field.
   *    \li \b DestField The policy for the destination field.
   *
   */
  template< typename SrcField,
	    typename DestField >
  class DivergenceAssembler
  {
  public:

    /** @brief Return the number of nonzero entries for this operator. */
    static int num_nonzeros(){return 2;}

    /**
     *  @brief Construct a DivergenceAssembler object.
     *
     *  @param dimExtent A vector with three elements indicating the
     *  domain extent (number of cells) in each coordinate direction.
     *
     *  @param cellFaceArea The cell area (assumed constant) for the
     *  face normal to the direction of this divergence operator.
     *
     *  @param cellVolume The volume of the cell (assumed constant).
     */
    DivergenceAssembler( const std::vector<int>& dimExtent,
			 const double cellFaceArea,
			 const double cellVolume );

    ~DivergenceAssembler(){}

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
    const std::vector<int> extent_;
    const double coefValue_;

  };


  //==================================================================




  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //                          Implementation
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




  //==================================================================


  //------------------------------------------------------------------
  template< typename SrcField, typename DestField >
  DivergenceAssembler<SrcField,DestField>::
  DivergenceAssembler( const std::vector<int>& dimExtent,
		       const double cellFaceArea,
		       const double cellVolume )
    : dim_        ( dimExtent ),
      indexHelper_( dimExtent ),
      extent_     ( dimExtent ),
      coefValue_  ( cellFaceArea/cellVolume )
  {
  }
  //--------------------------------------------------------------------
  template< typename SrcField, typename DestField >
  int
  DivergenceAssembler<SrcField,DestField>::
  get_ncols() const
  {
    return indexHelper_.get_ncol();
  }
  //--------------------------------------------------------------------
  template< typename SrcField, typename DestField >
  int
  DivergenceAssembler<SrcField,DestField>::
  get_nrows() const
  {
    int n=1;
    if( get_n_tot<SrcField>(dim_) > 1 ) n=indexHelper_.get_nrow();
    return n;
  }
  //--------------------------------------------------------------------
  template< typename SrcField, typename DestField >
  void
  DivergenceAssembler<SrcField,DestField>::
  get_row_entries( const int irow,
		   std::vector<double> & vals,
		   std::vector<int> & ixs ) const
  {
    switch( SrcField::Location::Dir::value ){
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
    switch( SrcField::Location::StagDir::value ){
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
      vals.push_back( -coefValue_ );
      vals.push_back(  coefValue_ );
    }
  }
  //------------------------------------------------------------------


}
}

#endif
