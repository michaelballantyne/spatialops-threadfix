#ifndef FVStaggeredBCTools_h
#define FVStaggeredBCTools_h

#include <boost/static_assert.hpp>

#include <FVStaggeredTools.h>
#include <SpatialOpsDefs.h>

#include <boost/static_assert.hpp>

namespace SpatialOps{
namespace FVStaggered{


  //------------------------------------------------------------------

  // given ijk indices that are zero based on the interior, this
  // produces the flat index that is 0-based in the ghost cell.
  template< typename FieldT, typename Dir >
  inline int get_ghost_flat_ix_src( const std::vector<int>& dim, int i, int j, int k )
  {
    // 1. shift index by number of ghost cells.
    // 2. shift index to first ghost cell (from interior cell)
    if( dim[0]>1 ){
      if( IsSameType<Dir,XDIR>::result )
	if( i==0 ) --i; else ++i;
      i += FieldT::Ghost::NM;
    }
    if( dim[1]>1 ){
      if( IsSameType<Dir,YDIR>::result )
	if( j==0 ) --j; else ++j;
      j += FieldT::Ghost::NM;
    }
    if( dim[2]>1 ){
      if( IsSameType<Dir,ZDIR>::result )
	if( k==0 ) --k; else ++k;
      k += FieldT::Ghost::NM;
    }

    return ijk2flat<FieldT,FieldT::Location::IsSurface>::value(dim,IndexTriplet(i,j,k));
  }

  //------------------------------------------------------------------

  template<typename FieldT, typename Dir>
  inline int get_ghost_flat_ix_dest( const std::vector<int>& dim, int i, int j, int k )
  {
    if( dim[0]>1 ){
      if( IsSameType<Dir,XDIR>::result ) if(i>0) ++i;
      i += FieldT::Ghost::NM;
    }

    if( dim[1]>1 ){
      if( IsSameType<Dir,YDIR>::result ) if(j>0) ++j;
      j += FieldT::Ghost::NM;
    }

    if( dim[2]>1 ){
      if( IsSameType<Dir,ZDIR>::result ) if(k>0) ++k;
      k += FieldT::Ghost::NM;
    }

    return ijk2flat<FieldT,FieldT::Location::IsSurface,Dir>::value(dim,IndexTriplet(i,j,k));
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

    const int ixf = get_ghost_flat_ix_src <SrcFieldT,Dir >( dim, i, j, k );
    int irow      = get_ghost_flat_ix_dest<DestFieldT,Dir>( dim, i, j, k );

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


  /**
   *  @author James C. Sutherland
   *  @date Jan, 2008
   *
   *  Imprints an operator with a given boundary condition.
   *
   *  @param bcOp The BC operator that we are using to apply the BC.
   *              For Dirichlet conditions, use an interpolant
   *              operator; for Neumann conditions use a gradient
   *              operator.
   *
   *  @param op   The operator to imprint with the boundary conditions.
   *
   *  @param i  The x-direction index for the cell we want to apply the BC to. Index is 0-based on patch interior.
   *  @param j  The y-direction index for the cell we want to apply the BC to. Index is 0-based on patch interior.
   *  @param k  The z-direction index for the cell we want to apply the BC to. Index is 0-based on patch interior.
   *
   *  @param bcVal The value for the boundary condition to set.
   *
   *  @par Template Parameters
   *   <ul>
   *
   *   <li> \b BCOpT Specifies the type of SpatialOperator that will be
   *   used to apply this BC.  This must define a few things:
   *     <ul>
   *     <li> \b DestFieldType The type for the destination field.
   *     <li> \b SrcFieldType  The type for the source field.
   *     </ul>
   *
   *   <li> \b Dir Specifies the direction (boundary normal) for this BC.
   *
   *   <li> \b OpT Specifies the type of SpatialOperator that will be
   *   imprinting with this BC.
   *
   *   </ul>
   */
  template< typename BCOpT, typename Dir, typename OpT >
  void imprint_bc_on_op( const BCOpT& bcOp,
			 const int i, const int j, const int k,
			 const std::vector<int>& dim,
			 const double bcVal,
			 const typename BCOpT::SrcFieldType& f,
			 OpT& op,
			 double& rhs )
  {
    static int ncolMax = 10; // maximum number of nonzero columns
    struct BCInfo{ int ix; double coef; };
    BCInfo bcinfo[ncolMax];
    int nbcinfo = 0;

    // get the index into the field value at this point.
    const int ixf  = get_ghost_flat_ix_src <typename BCOpT::SrcFieldType, Dir>( dim, i, j, k );
    int irow       = get_ghost_flat_ix_dest<typename BCOpT::DestFieldType,Dir>( dim, i, j, k );

    int ncol; double*vals=0; int*ixs=0;
    bcOp.get_linalg_mat().ExtractMyRowView( irow, ncol, vals, ixs );
    assert( ncol<ncolMax );
    assert( OpT::SrcGhost::NM == OpT::SrcGhost::NP );
    assert( OpT::SrcGhost::NM == 1 );

    double prodsum=0.0; double ghostcoeff=0.0;
    BCInfo* bci = &bcinfo[0];
    for( int icol=0; icol<ncol; ++icol ){
      if (*ixs == ixf){
	ghostcoeff = *vals;
      }
      else{
	prodsum += (*vals) * f[*ixs];
	bci->coef = *vals;
	bci->ix   = *ixs;
	++bci;
	++nbcinfo;
      }
      ++ixs; ++vals;
    }

    assert( ghostcoeff != 0.0 );
    //     const double ghostVal = (bcVal - prodsum) / ghostcoeff;



    // currently, we are basically assuming that the operator here is
    // going to a linear system, and that it is like a Laplacian...
    BOOST_STATIC_ASSERT( bool( IsSameType<typename OpT::SrcFieldType, typename OpT::DestFieldType>::result ) );

    //
    // now set the operator value.  We must potentially alter each coefficient in this row.
    //
    const int ig = get_ghost_flat_ix_src <typename OpT::SrcFieldType, Dir>( dim, i, j, k );
    irow         = get_ghost_flat_ix_dest<typename OpT::DestFieldType,NODIR>( dim, i, j, k );

    op.get_linalg_mat().ExtractMyRowView( irow, ncol, vals, ixs );
    assert( ncol<ncolMax );
    assert( OpT::SrcGhost::NM == OpT::SrcGhost::NP );
    assert( OpT::SrcGhost::NM == 1 );

    double Sg = 0.0;
    int* iix=ixs;
    double* ival=vals;
    for( int icol=0; icol<ncol; ++icol ){
      if( *iix==ig ){
	Sg = *ival;
	*ival = 0.0;
	break;
      }
      ++iix; ++ival;
    }

    //
    // augment the RHS value
    //
    rhs -= bcVal/ghostcoeff*Sg;

    // imprint the LHS.
    ival=vals;
    iix=ixs;
    for( int icol=0; icol<ncol; ++icol ){
      if( *iix!=ig ){
	for( BCInfo* bci=&bcinfo[0]; bci!=&bcinfo[nbcinfo]; ++bci ){
	  if( bci->ix == *iix ){
	    *ival -= bci->coef/ghostcoeff * Sg;
	    break;
	  }
	} // for
      } // if
      ++iix; ++ival;
    } // column loop
    
  }


}// namespace FVStaggered
}// namespace SpatialOps

#endif
