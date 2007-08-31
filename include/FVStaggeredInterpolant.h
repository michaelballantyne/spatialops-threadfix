#ifndef FVStaggeredInterpolant_h
#define FVStaggeredInterpolant_h

#include <FVStaggeredIndexHelper.h>
#include <FVStaggeredTypes.h>

#include <SpatialField.h>
#include <SpatialOperator.h>
#include <SpatialOpsDefs.h>

namespace SpatialOps{


  // forward declaration.
  namespace FVStaggered{ template<typename T1,typename T2> class LinearInterpolantAssembler; }

  // this is required for the SpatialOperator class.  It specifies
  // that we should use the LinearInterpolantAssembler to construct
  // Interpolant objects.
  template< typename SrcFieldT, typename DestFieldT >
  struct OpAssemblerSelector< FVStaggered::Interpolant, SrcFieldT, DestFieldT >
  {
    typedef FVStaggered::LinearInterpolantAssembler<SrcFieldT,DestFieldT>  Assembler;
  };


namespace FVStaggered{


  /**
   *  @class  LinearInterpolantAssembler
   *  @author James C. Sutherland
   *
   *  @brief Assembles linear interpolant operators for staggered meshes.
   *
   *  An assembler for a linear interpolant operator to be used to
   *  construct a SpatialOperator.  This conforms to the requirements
   *  for a SpatialOperator Assembler as defined in the documentation
   *  on the OpAssemblerSelector.
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

    int num_nonzeros() const;

    /**
     *  @brief Construct a LinearInterpolantAssembler.
     *  @param dimExtent 3-component vector indicating the number of
     *  interior points in each coordinate direction.
     */
    LinearInterpolantAssembler( const std::vector<int>& dimExtent );

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
    void get_ghost_cols( std::set<int>& ghostCols ) const
    {
      ghostCols = get_ghost_set<SrcFieldT>(dim_);
    }

    /**
     *  @brief Obtain the set of row indices corresponding to ghost entries.
     */
    void get_ghost_rows( std::set<int>& ghostRows ) const
    {
      ghostRows = get_ghost_set<DestFieldT>(dim_);
    }

  private:

    const IndexHelper<SrcFieldT,DestFieldT> indexHelper_;
    const std::vector<int>& dim_;

  };

  //==================================================================

  // still need to implement these:
  template<>  class LinearInterpolantAssembler<SVolField,XSurfField >{};
  template<>  class LinearInterpolantAssembler<SVolField,XSurfXField>{};
  template<>  class LinearInterpolantAssembler<SVolField,XSurfYField>{};
  template<>  class LinearInterpolantAssembler<SVolField,XSurfZField>{};

  template<>  class LinearInterpolantAssembler<SVolField,YSurfField >{};
  template<>  class LinearInterpolantAssembler<SVolField,YSurfXField>{};
  template<>  class LinearInterpolantAssembler<SVolField,YSurfYField>{};
  template<>  class LinearInterpolantAssembler<SVolField,YSurfZField>{};

  template<>  class LinearInterpolantAssembler<SVolField,ZSurfField >{};
  template<>  class LinearInterpolantAssembler<SVolField,ZSurfXField>{};
  template<>  class LinearInterpolantAssembler<SVolField,ZSurfYField>{};
  template<>  class LinearInterpolantAssembler<SVolField,ZSurfZField>{};

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
  LinearInterpolantAssembler( const std::vector<int>& dimExtent )
    : indexHelper_( dimExtent ),
      dim_( dimExtent )
  {
  }
  //------------------------------------------------------------------
  template<typename SrcFieldT, typename DestFieldT>
  int
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
    switch( DestFieldT::Location::Dir::value ){
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
    switch( DestFieldT::Location::StagDir::value ){
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
    for( std::vector<int>::size_type i=0; i<ixs.size(); ++i ){
      vals.push_back( 0.5 );
    }
  }
  //------------------------------------------------------------------


  // specializations for XVolField->SSurfXField

  template<>
  inline int LinearInterpolantAssembler<XVolField,SSurfXField>::
  get_ncols() const{ return get_n_tot<XVolField>(dim_); }

  template<>
  inline int LinearInterpolantAssembler<XVolField,SSurfXField>::
  get_nrows() const{ return get_n_tot<SSurfXField>(dim_); }

  template<>
  inline int LinearInterpolantAssembler<XVolField,SSurfXField>::
  num_nonzeros() const{ return 1; }

  template<>
  inline void
  LinearInterpolantAssembler<XVolField,SSurfXField>::
  get_row_entries( const int irow, std::vector<double> & vals, std::vector<int> & ixs ) const
  {
    ixs.push_back(irow);
    vals.push_back(1.0);
  }

  //------------------------------------------------------------------

  // specializations for XVolField->XSurfXField

  //------------------------------------------------------------------

  // specializations for YVolField->SSurfYField

  template<>
  inline int LinearInterpolantAssembler<YVolField,SSurfYField>::
  get_ncols() const{ return get_n_tot<YVolField>(dim_); }

  template<>
  inline int LinearInterpolantAssembler<YVolField,SSurfYField>::
  get_nrows() const{ return get_n_tot<SSurfYField>(dim_); }

  template<>
  inline int LinearInterpolantAssembler<YVolField,SSurfYField>::
  num_nonzeros() const{ return 1; }

  template<>
  inline void
  LinearInterpolantAssembler<YVolField,SSurfYField>::
  get_row_entries( const int irow, std::vector<double> & vals, std::vector<int> & ixs ) const
  {
    ixs.push_back(irow);
    vals.push_back(1.0);
  }

  //------------------------------------------------------------------

  // specializations for ZVolField->SSurfZField

  template<>
  inline int LinearInterpolantAssembler<ZVolField,SSurfZField>::
  get_ncols() const{ return get_n_tot<ZVolField>(dim_); }

  template<>
  inline int LinearInterpolantAssembler<ZVolField,SSurfZField>::
  get_nrows() const{ return get_n_tot<SSurfZField>(dim_); }

  template<>
  inline int LinearInterpolantAssembler<ZVolField,SSurfZField>::
  num_nonzeros() const{ return 1; }

  template<>
  inline void
  LinearInterpolantAssembler<ZVolField,SSurfZField>::
  get_row_entries( const int irow, std::vector<double> & vals, std::vector<int> & ixs ) const
  {
    ixs.push_back(irow);
    vals.push_back(1.0);
  }

  //==================================================================


} // namespace FVStaggered
} // namespace SpatialOps

#endif
