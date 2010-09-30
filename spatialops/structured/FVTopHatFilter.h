#ifndef FVStructuredFilter_h
#define FVStructuredFilter_h

#include <spatialops/SpatialOpsConfigure.h>

#include <spatialops/SpatialOperator.h>

#include <spatialops/structured/FVTools.h>

namespace SpatialOps{

  // forward declaration.
  namespace structured{ template<typename T1> class TopHatFilterAssembler; }

  // this is required for the SpatialOperator class.  It specifies
  // that we should use the TopHatFilterAssembler class to construct
  // filter objects.
  template< typename FieldT >
  struct OpAssemblerSelector< Filter, FieldT, FieldT >
  {
    typedef structured::TopHatFilterAssembler<FieldT>  Assembler;
  };


namespace structured{

  /**
   *  @class TopHatFilterAssembler
   *  @author James C. Sutherland
   *  @brief Assembles a tophat filter on a uniform mesh.  Coefficients are all constant.
   *
   *  NOTE: this filter will shrink symmetrically as it approaches a
   *  patch boundary.  This means that if only one layer of ghost
   *  cells are used, then the filter will be at most 3 points at
   *  the boundary.
   */
  template< typename FieldT >
  class TopHatFilterAssembler
  {
  public:

    /**
     *  @brief Build a filter with a specified number of points.
     *  @param npts The number of points to include in the filter stencil.
     */
    TopHatFilterAssembler( const int npts,
                           const IntVec& dimExtent,
                           const bool hasPlusXSideFaces = true,
                           const bool hasPlusYSideFaces = true,
                           const bool hasPlusZSideFaces = true  );

    /**
     *  @brief Build a filter with a specified width.
     *  @param width The width of the filter.
     */
    TopHatFilterAssembler( const double width,
                           const std::vector<double>& meshSpacing,
                           const IntVec& dimExtent,
                           const bool hasPlusXSideFaces = true,
                           const bool hasPlusYSideFaces = true,
                           const bool hasPlusZSideFaces = true  );

    ~TopHatFilterAssembler(){}

    unsigned int num_nonzeros() const;

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
      ghostCols = get_ghost_set<FieldT>( dim_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    }

    /**
     *  @brief Obtain the set of row indices corresponding to ghost entries.
     */
    void get_ghost_rows( std::set<size_t>& ghostRows ) const
    {
      ghostRows = get_ghost_set<FieldT>( dim_, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ );
    }

  private:

    void get_npts( const IntVec& loc,
                   int& ntot,
                   IntVec& nxyz ) const;

    const IntVec dim_;
    const bool hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_;
    std::vector<unsigned int> npts_;
    
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
  TopHatFilterAssembler<FieldT>::
  TopHatFilterAssembler( const double width,
                         const std::vector<double>& meshSpacing,
                         const IntVec& dimExtent,
                         const bool hasPlusXSideFaces,
                         const bool hasPlusYSideFaces,
                         const bool hasPlusZSideFaces )
    : dim_( dimExtent ),
      hasPlusXSideFaces_( hasPlusXSideFaces ),
      hasPlusYSideFaces_( hasPlusYSideFaces ),
      hasPlusZSideFaces_( hasPlusZSideFaces ),
      npts_(3,0)
  {
    for( size_t i=0; i<3; ++i ){
      npts_[i] = 1 + width / meshSpacing[i];
      if( npts_[i]%2 ) ++npts_[i];
    }
  }

  //------------------------------------------------------------------

  template<typename FieldT>
  TopHatFilterAssembler<FieldT>::
  TopHatFilterAssembler( const int npts,
                         const IntVec& dimExtent,
                         const bool hasPlusXSideFaces,
                         const bool hasPlusYSideFaces,
                         const bool hasPlusZSideFaces )
    : dim_( dimExtent ),
      hasPlusXSideFaces_( hasPlusXSideFaces ),
      hasPlusYSideFaces_( hasPlusYSideFaces ),
      hasPlusZSideFaces_( hasPlusZSideFaces ),
      npts_( 3, npts )
  {
  }

  //------------------------------------------------------------------

  template<typename FieldT>
  int
  TopHatFilterAssembler<FieldT>::
  get_ncols() const
  {
    return get_ntot_with_ghost<FieldT>(dim_,hasPlusXSideFaces_,hasPlusYSideFaces_,hasPlusZSideFaces_);
  }

  //------------------------------------------------------------------

  template<typename FieldT>
  int
  TopHatFilterAssembler<FieldT>::
  get_nrows() const
  {
    return get_ncols();
  }

  //------------------------------------------------------------------

  template<typename FieldT>
  unsigned int
  TopHatFilterAssembler<FieldT>::
  num_nonzeros() const
  {
    unsigned int ntot = 1;
    for( std::vector<unsigned int>::const_iterator i=npts_.begin(); i!=npts_.end(); ++i ){
      ntot *= *i;
    }
    return ntot;
  }

  //------------------------------------------------------------------

  template<typename FieldT>
  void
  TopHatFilterAssembler<FieldT>::
  get_row_entries( const int irow,
                   std::vector<double> & vals,
                   std::vector<int> & ixs ) const
  {
    const IntVec ijk( flat2ijk<FieldT>::value( dim_, irow, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ ) );

    int ntot=1;
    IntVec nxyz(1,1,1), nside(1,1,1);
    get_npts(ijk,ntot,nxyz);
    for( size_t i=0; i<3; ++i ) nside[i] = nxyz[i]/2;
    const double fac = 1.0/double( ntot );

    for( int ipts=-nside[0]; ipts<=nside[0]; ++ipts ){
      for( int jpts=-nside[1]; jpts<=nside[1]; ++jpts ){
        for( int kpts=-nside[2]; kpts<=nside[2]; ++kpts ){
          const IntVec ixt( ijk[0]+ipts, ijk[1]+jpts, ijk[2]+kpts );
          ixs.push_back( ijk2flat<FieldT>::value( dim_, ixt, hasPlusXSideFaces_, hasPlusYSideFaces_, hasPlusZSideFaces_ ) );
          vals.push_back( fac );
        }
      }
    }
  }

  //------------------------------------------------------------------

  template<typename FieldT>
  void
  TopHatFilterAssembler<FieldT>::
  get_npts( const IntVec& loc,
            int& ntot,
            IntVec& nxyz ) const
  {
    ntot = 1;
    if( dim_[0] > 1 ){
      const int nx = get_nx_with_ghost<FieldT>(dim_[0],hasPlusXSideFaces_);
      const int nside = npts_[0] / 2;
      nxyz[0] = 1 + 2 * std::min( std::min( nside, loc[0] ),
                                  std::min( nside, nx-loc[0]-1 ) );
      ntot *= nxyz[0];
    }

    if( dim_[1] > 1 ){
      const int ny = get_ny_with_ghost<FieldT>(dim_[1],hasPlusYSideFaces_);
      const int nside = npts_[1] / 2;
      nxyz[1] = 1 + 2 * std::min( std::min( nside, loc[1]),
                                  std::min( nside, ny-loc[1]-1) );
      ntot *= nxyz[1];
    }

    if( dim_[2] > 1 ){
      const int nz = get_nz_with_ghost<FieldT>(dim_[2],hasPlusZSideFaces_);
      const int nside = npts_[2] / 2;
      nxyz[2] = 1 + 2 * std::min( std::min( nside, loc[2]),
                                  std::min( nside, nz-loc[2]-1) );
      ntot *= nxyz[2];
    }
  }

  //------------------------------------------------------------------

} // namespace structured
} // namespace Spatialops

#endif