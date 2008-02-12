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
  inline int get_ghost_flat_ix_src( const std::vector<int>& dim,
				    const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ,
				    int i, int j, int k )
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
    return ijk2flat<FieldT>::value(dim,IndexTriplet(i,j,k),bcFlagX,bcFlagY,bcFlagZ);
  }

  //------------------------------------------------------------------

  template<typename FieldT, typename Dir>
  inline int get_ghost_flat_ix_dest( const std::vector<int>& dim,
				     const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ,
				     int i, int j, int k )
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
    return ijk2flat<FieldT>::value(dim,IndexTriplet(i,j,k),bcFlagX,bcFlagY,bcFlagZ);
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
   *  @param dim A vector containing the number of cells in each
   *  coordinate direction.  This is a three-component vector.
   *
   *  @param bcFlagX A boolean flag to indicate if this patch is on a
   *  +x side physical boundary.  If so, then it is assumed that there
   *  is an extra face on that side of the domain, and face variable
   *  dimensions will be modified accordingly.
   *
   *  @param bcFlagY A boolean flag to indicate if this patch is on a
   *  +y side physical boundary.  If so, then it is assumed that there
   *  is an extra face on that side of the domain, and face variable
   *  dimensions will be modified accordingly.
   *
   *  @param bcFlagZ A boolean flag to indicate if this patch is on a
   *  +z side physical boundary.  If so, then it is assumed that there
   *  is an extra face on that side of the domain, and face variable
   *  dimensions will be modified accordingly.
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
			const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ,
			const double bcVal,
			typename OpT::SrcFieldType& f )
  {
    typedef typename OpT::SrcFieldType  SrcFieldT;
    typedef typename OpT::DestFieldType DestFieldT;

    const int ixf = get_ghost_flat_ix_src <SrcFieldT,Dir >( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k );
    int irow      = get_ghost_flat_ix_dest<DestFieldT,Dir>( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k );

    // NOTE: This will NOT work in the case where we have multiple ghost cells!
    assert( OpT::SrcGhost::NM == OpT::SrcGhost::NP );
    assert( OpT::SrcGhost::NM == 1 );

    const typename OpT::MatrixRow row = op.get_row(irow);
    typename       OpT::const_column_iterator icol =row.begin();
    const typename OpT::const_column_iterator icole=row.end();

    double prodsum=0.0; double ghostcoeff=0.0;
    for( ; icol!=icole; ++icol ){
      if (icol.index() == uint(ixf) )
	ghostcoeff = *icol;
      else
	prodsum += *icol *f[icol.index()];
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
   *  @param dim A vector containing the number of cells in each
   *  coordinate direction.  This is a three-component vector.
   *
   *  @param bcFlagX A boolean flag to indicate if this patch is on a
   *  +x side physical boundary.  If so, then it is assumed that there
   *  is an extra face on that side of the domain, and face variable
   *  dimensions will be modified accordingly.
   *
   *  @param bcFlagY A boolean flag to indicate if this patch is on a
   *  +y side physical boundary.  If so, then it is assumed that there
   *  is an extra face on that side of the domain, and face variable
   *  dimensions will be modified accordingly.
   *
   *  @param bcFlagZ A boolean flag to indicate if this patch is on a
   *  +z side physical boundary.  If so, then it is assumed that there
   *  is an extra face on that side of the domain, and face variable
   *  dimensions will be modified accordingly.
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
			 const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ,
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
    const int ixf  = get_ghost_flat_ix_src <typename BCOpT::SrcFieldType, Dir>( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k );
    int irow       = get_ghost_flat_ix_dest<typename BCOpT::DestFieldType,Dir>( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k );


    const typename BCOpT::MatrixRow bcrow = bcOp.get_row(irow);
    typename       BCOpT::const_column_iterator icolbc = bcrow.begin();
    const typename BCOpT::const_column_iterator icolbce= bcrow.end();

    double prodsum=0.0; double ghostcoeff=0.0;
    BCInfo* bci = &bcinfo[0];
    for( ; icolbc!=icolbce; ++icolbc ){
      if( icolbc.index() == uint(ixf) ){
	ghostcoeff = *icolbc;
      }
      else{
	prodsum += *icolbc * f[icolbc.index()];
	bci->coef = *icolbc;
	bci->ix   = icolbc.index();
	++bci;
	++nbcinfo;
      }
    }

    assert( ghostcoeff != 0.0 );


    // currently, we are basically assuming that the operator here is
    // going to a linear system, and that it is like a Laplacian...
    BOOST_STATIC_ASSERT( bool( IsSameType<typename OpT::SrcFieldType, typename OpT::DestFieldType>::result ) );

    //
    // now set the operator value.  We must potentially alter each coefficient in this row.
    //
    const int ig = get_ghost_flat_ix_src <typename OpT::SrcFieldType, Dir  >( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k );
    irow         = get_ghost_flat_ix_dest<typename OpT::DestFieldType,NODIR>( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k );

    typename OpT::MatrixRow row = op.get_row(irow);
    typename OpT::column_iterator icol=row.begin();
    const typename OpT::column_iterator icole=row.end();

    double Sg = 0.0;
    for( ; icol!=icole; ++icol ){
      if( icol.index() == uint(ig) ){
	double& val = *icol;
	Sg = val;
	val = 0.0;
	break;
      }
    }

    //
    // augment the RHS value
    //
    rhs -= bcVal/ghostcoeff*Sg;

    // imprint the LHS.
    typename OpT::column_iterator ic=row.begin();
    for( ; ic!=icole; ++ic ){
      const int ix = ic.index();
      if( ix!=ig ){
	for( BCInfo* bci=&bcinfo[0]; bci!=&bcinfo[nbcinfo]; ++bci ){
	  if( bci->ix == ix ){
	    *ic -= bci->coef/ghostcoeff * Sg;
	    break;
	  }
	} // for
      } // if
    } // column loop
    
  }


}// namespace FVStaggered
}// namespace SpatialOps

#endif
