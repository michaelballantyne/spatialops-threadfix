#ifndef FVRestrictOperator_h
#define FVRestrictOperator_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/FVTools.h>
#include <spatialops/SpatialOperator.h>

namespace SpatialOps{

  // forward declaration.
  template<typename T1> class RestrictionAssembler;

  // this is required for the SpatialOperator class.  It specifies
  // that we should use the RestrictionAssembler class to construct
  // restriction objects.
  template< typename FieldT >
  struct OpAssemblerSelector< Restriction, FieldT, FieldT >
  {
    typedef RestrictionAssembler<FieldT>  Assembler;
  };


  /**
   *  @class  RestrictionAssembler
   *  @author James C. Sutherland
   *  @date   April, 2009
   *
   *  @brief Assembles a restriction operator to move a field from a
   *  coarse mesh to a fine one.  This uses direct-injection (no
   *  interpolation).  It also assumes that restriction occurs on a
   *  single dimension only.
   */
  template<typename FieldT>
  class RestrictionAssembler
  {
  public:

    /**
     *  Build an assembler for a restriction operator.
     *
     *  @param dimExtentSrc The vector indicating the number of
     *         points in each ordinal direction for the source
     *         field.
     *
     *  @param dimExtentDest The vector indicating the number of
     *         points in each ordinal direction for the destination
     *         field.  Note that this vector should only differ from
     *         the dimExtentSrc parameter in a single dimension.
     *
     *  @param hasPlusXSideFaces Indicates whether there is a
     *         physical boundary on the right side of the patch in
     *         the x-direction.
     *  @param hasPlusYSideFaces Indicates whether there is a
     *         physical boundary on the right side of the patch in
     *         the y-direction.
     *  @param hasPlusZSideFaces Indicates whether there is a
     *         physical boundary on the right side of the patch in
     *         the z-direction.
     *
     */
    RestrictionAssembler( const std::vector<int>& dimExtentSrc,
                          const std::vector<int>& dimExtentDest,
                          const bool hasPlusXSideFaces = true,
                          const bool hasPlusYSideFaces = true,
                          const bool hasPlusZSideFaces = true  );


    ~RestrictionAssembler(){}

    unsigned int num_nonzeros() const{ return 1; }

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
      ghostCols = FVStaggered::get_ghost_set<FieldT>( dimSrc_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    }

    /**
     *  @brief Obtain the set of row indices corresponding to ghost entries.
     */
    void get_ghost_rows( std::set<int>& ghostRows ) const
    {
      ghostRows = FVStaggered::get_ghost_set<FieldT>( dimDest_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    }

  private:

    static unsigned int get_active_dim( const std::vector<int>& dimSrc,
                                        const std::vector<int>& dimDest );

    const unsigned int activeDim_;
    const std::vector<int> dimSrc_, dimDest_;
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


  template<typename FieldT>
  RestrictionAssembler<FieldT>::
  RestrictionAssembler( const std::vector<int>& dimExtentSrc,
                        const std::vector<int>& dimExtentDest,
                        const bool hasPlusXSideFaces,
                        const bool hasPlusYSideFaces,
                        const bool hasPlusZSideFaces )
    : activeDim_( get_active_dim(dimExtentSrc,dimExtentDest) ),
        
      dimSrc_ ( dimExtentSrc  ),
      dimDest_( dimExtentDest ),
        
      hasPlusXSideFaces_( hasPlusXSideFaces ),
      hasPlusYSideFaces_( hasPlusYSideFaces ),
      hasPlusZSideFaces_( hasPlusZSideFaces )
  {
    int ndifferent = 0;
    for( int i=0; i<3; ++i )
      if( dimExtentDest[i] != dimExtentSrc[i] )  ++ndifferent;
    assert( ndifferent <= 1 );
    assert( dimExtentSrc[activeDim_] < dimExtentDest[activeDim_] );
  }

  //----------------------------------------------------------------

  template<typename FieldT>
  unsigned int
  RestrictionAssembler<FieldT>::
  get_active_dim( const std::vector<int>& dimSrc,
                  const std::vector<int>& dimDest )
  {
    unsigned int dim = 99;
    for( size_t i=0; i<dimSrc.size(); ++i ){
      if( dimSrc[i] != dimDest[i] ){
        dim = i;
        break;
      }
    }
    assert( dim < 3 );
    return dim;
  }

  //----------------------------------------------------------------

  template<typename FieldT>
  int
  RestrictionAssembler<FieldT>::
  get_ncols() const
  {
    return FVStaggered::get_n_tot<FieldT>( dimSrc_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
  }

  //----------------------------------------------------------------

  template<typename FieldT>
  int
  RestrictionAssembler<FieldT>::
  get_nrows() const
  {
    return FVStaggered::get_n_tot<FieldT>( dimDest_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
  }

  //----------------------------------------------------------------

  template<typename FieldT>
  void
  RestrictionAssembler<FieldT>::
  get_row_entries( const int irow,
                   std::vector<double> & vals,
                   std::vector<int> & ixs ) const
  {
    const FVStaggered::IndexTriplet ijkDest( FVStaggered::flat2ijk<FieldT>::value( dimDest_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ ) );

    const int ngM = FieldT::Ghost::NM; // number of ghost points
    const int ngP = FieldT::Ghost::NP; // number of ghost points

    const int nsrc  = dimSrc_ [activeDim_];
    const int ndest = dimDest_[activeDim_];
    const int idest = ijkDest [activeDim_] - ngM;

    FVStaggered::IndexTriplet ijkSrc = ijkDest;

    if( ijkDest[activeDim_] < ngM )
      ijkSrc[activeDim_] = 0;
    else if( ijkDest[activeDim_]-ngM >= dimDest_[activeDim_] )
      ijkSrc[activeDim_] = dimSrc_[activeDim_]+ngP;
    else
      ijkSrc[activeDim_] = int( float(idest*nsrc)/float(ndest) + ngM );

    const int index = FVStaggered::ijk2flat<FieldT>::value( dimSrc_, ijkSrc, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );

    ixs.push_back( index );
    vals.push_back( 1.0 );
  }

  //----------------------------------------------------------------

} // namespace SpatialOps

#endif  // FVRestrictOperator_h
