#ifndef FVStaggeredBCTools_h
#define FVStaggeredBCTools_h

#include <FVStaggeredTools.h>

#include <boost/static_assert.hpp>

namespace SpatialOps{
namespace FVStaggered{


  //------------------------------------------------------------------

  template< typename SrcFieldT, typename DestFieldT, typename Dir >
  inline void bc_adjust_ijk( IndexTriplet& st,
			     IndexTriplet& dt,
			     const std::vector<int>& dim )
  {
    switch( Dir::value ){

    case XDIR::value :
      if( dim[0]<=1 ) return;

      if( st.i==0 ){ --st.i; }
      else      { ++st.i; ++dt.i; }
      st.i +=  SrcFieldT::Ghost::NM;
      dt.i += DestFieldT::Ghost::NM;

      if( dim[1]>1 ){
	st.j +=  SrcFieldT::Ghost::NM;
	dt.j += DestFieldT::Ghost::NM;
      }

      if( dim[2]>1 ){
	st.k +=  SrcFieldT::Ghost::NM;
	dt.k += DestFieldT::Ghost::NM;
      }
      break;

    case YDIR::value :
      if( dim[1]<=1 ) return;

      if( st.j==0 ){ --st.j; }
      else      { ++st.j; ++dt.j; }
      st.j +=  SrcFieldT::Ghost::NM;
      dt.j += DestFieldT::Ghost::NM;

      if( dim[0]>1 ){
	st.i +=  SrcFieldT::Ghost::NM;
	dt.i += DestFieldT::Ghost::NM;
      }

      if( dim[2]>1 ){
	st.k +=  SrcFieldT::Ghost::NM;
	dt.k += DestFieldT::Ghost::NM;
      }
      break;

    case ZDIR::value :
      if( dim[2]<=1 ) return;

      if( st.k==0 ){ --st.k; }
      else         { ++st.k; ++dt.k; }
      st.k +=  SrcFieldT::Ghost::NM;
      dt.k += DestFieldT::Ghost::NM;

      if( dim[0]>1 ){
	st.i +=  SrcFieldT::Ghost::NM;
	dt.i += DestFieldT::Ghost::NM;
      }

      if( dim[1]>1 ){
	st.j +=  SrcFieldT::Ghost::NM;
	dt.j += DestFieldT::Ghost::NM;
      }
      break;

    default:
      assert(1);
    }
  }

  //------------------------------------------------------------------

  /**
   *  @brief Apply either a Dirichlet or Neumann BC (depending on the
   *  operator supplied) to the given field.
   *
   *  @param op The operator to use in applying the BC.  Supplying a
   *            gradient operator results in a Neumann BC; supplying
   *            an interpolant operator results in a Dirichlet BC.
   *
   *  @param i  The x-direction index for the cell we want to apply the BC to. Index is 0-based on patch interior.
   *  @param j  The y-direction index for the cell we want to apply the BC to. Index is 0-based on patch interior.
   *  @param k  The z-direction index for the cell we want to apply the BC to. Index is 0-based on patch interior.
   *
   *  @param bcVal The value for the boundary condition to set.
   *
   *  @param f  The field we want to set the BC on.
   *
   *  @par Template Parameters
   *   <ul>
   *
   *   <li> \b OpT Specifies the type of SpatialOperator that will be
   *   used to apply this BC.  This must define a few things:
   *     <ul>
   *     <li> \b DestFieldType The type for the destination field.
   *     <li> \b SrcFieldType  The type for the source field.
   *     </ul>
   *
   *   <li> \b Dir Specifies the direction (boundary normal) for this BC.
   *
   *   </ul>
   *
   *  @todo Implement and test approach for vol->surface interpolant
   *  operators.  This will be required for dirichlet conditions
   */
  template< typename OpT, typename Dir >
  void assign_bc_point( const OpT& op,
			const int i,
			const int j,
			const int k,
			const std::vector<int>& dim,
			const double bcVal,
			typename OpT::SrcFieldType& f )
  {
    typedef typename OpT::SrcFieldType  SrcFieldT;
    typedef typename OpT::DestFieldType DestFieldT;

    IndexTriplet st(i,j,k), dt(i,j,k);

    // adjust field indices to include ghost cells
    bc_adjust_ijk<SrcFieldT,DestFieldT,Dir>( st, dt, dim );

    const int ixf  = ijk2flat< SrcFieldT, SrcFieldT::Location::IsSurface    >::value( dim, st );
    int irow       = ijk2flat<DestFieldT,DestFieldT::Location::IsSurface,Dir>::value( dim, dt );

    int ncol; double*vals=0; int*ixs=0;
    op.get_linalg_mat().ExtractMyRowView( irow, ncol, vals, ixs );

    // NOTE: This will NOT work in the case where we have multiple ghost cells!
    assert( OpT::SrcGhost::NM == OpT::SrcGhost::NP );
    assert( OpT::SrcGhost::NM == 1 );

    double prodsum=0.0; double ghostcoeff=0.0;
    for( int icol=0; icol<ncol; icol++ ){
      if (*ixs == ixf)	ghostcoeff = *vals;
      else              prodsum += (*vals) *f[*ixs];
      ++ixs; ++vals;
    }

    assert( ghostcoeff != 0.0 );
    f[ixf] = (bcVal - prodsum) / ghostcoeff;
  }


}// namespace FVStaggered
}// namespace SpatialOps

#endif
