#ifndef FVSS_DivergenceAssembler_h
#define FVSS_DivergenceAssembler_h

#include <spatialops/SpatialOpsConfigure.h>

#include <spatialops/SpatialField.h>
#include <spatialops/SpatialOperator.h>

#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/FVStaggeredIndexHelper.h>


namespace SpatialOps{

  // forward declaration.
  namespace structured{ template<typename T1,typename T2> class DivergenceAssembler; }

  template< typename SrcField, typename DestField >
  struct OpAssemblerSelector< Divergence, SrcField, DestField >
  {
    typedef structured::DivergenceAssembler<SrcField,DestField>  Assembler;
  };


namespace structured{

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
    static unsigned int num_nonzeros(){return 2;}

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
    DivergenceAssembler( const std::vector<int>& dimExtent,
                         const double cellFaceArea,
                         const double cellVolume,
                         const bool hasPlusXSideFaces,
                         const bool hasPlusYSideFaces,
                         const bool hasPlusZSideFaces );

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

    const std::vector<int>& dim_;
    const IndexHelper<SrcField,DestField> indexHelper_;
    const std::vector<int> extent_;
    const double coefValue_;
    const bool hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_;
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
                       const double cellVolume,
                       const bool hasPlusXSideFaces,
                       const bool hasPlusYSideFaces,
                       const bool hasPlusZSideFaces )
    : dim_        ( dimExtent ),
      indexHelper_( dimExtent, hasPlusXSideFaces, hasPlusYSideFaces, hasPlusZSideFaces ),
      extent_     ( dimExtent ),
      coefValue_  ( cellFaceArea/cellVolume ),
      hasPlusXSideFaces_( hasPlusXSideFaces ),
      hasPlusYSideFaces_( hasPlusYSideFaces ),
      hasPlusZSideFaces_( hasPlusZSideFaces )
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
    if( get_n_tot<SrcField>(dim_,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_) > 1 )
      n=indexHelper_.get_nrow();
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
    switch( int(SrcField::Location::FaceDir::value) ){
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
    switch( int(SrcField::Location::StagLoc::value) ){
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


} // namespace staggered
} // namespace SpatialOps

#endif
