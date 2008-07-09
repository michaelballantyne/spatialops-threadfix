#ifndef FVStaggeredBCTools_h
#define FVStaggeredBCTools_h

#include <boost/static_assert.hpp>

#include <FVStaggeredTools.h>
#include <SpatialOpsDefs.h>

#include <boost/static_assert.hpp>

namespace SpatialOps{
namespace FVStaggered{

  //------------------------------------------------------------------

  template< typename FieldT, typename Dir >
  inline void shift_ghost_ix( const std::vector<int>& dim, IndexTriplet& ijk )
  {
    if( IsSameType<Dir,XDIR>::result ) if( dim[0]>1 ){ if(ijk.i==0) --ijk.i; else ++ijk.i; }
    if( IsSameType<Dir,YDIR>::result ) if( dim[1]>1 ){ if(ijk.j==0) --ijk.j; else ++ijk.j; }
    if( IsSameType<Dir,ZDIR>::result ) if( dim[2]>1 ){ if(ijk.k==0) --ijk.k; else ++ijk.k; }
  }

  template< typename FieldT, typename Dir >
  inline void shift_ghost_ix_dest( const std::vector<int>& dim, IndexTriplet& ijk )
  {
    if( IsSameType<Dir,XDIR>::result ) if( dim[0]>1 ){ if(ijk.i>0) ++ijk.i; }
    if( IsSameType<Dir,YDIR>::result ) if( dim[1]>1 ){ if(ijk.j>0) ++ijk.j; }
    if( IsSameType<Dir,ZDIR>::result ) if( dim[2]>1 ){ if(ijk.k>0) ++ijk.k; }
  }

  //------------------------------------------------------------------

  // given ijk indices that are zero based on the interior, this
  // produces the flat index that is 0-based in the ghost cell.
  template< typename FieldT >
  inline int get_ghost_flat_ix( const std::vector<int>& dim,
				const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ,
				int i, int j, int k )
  {
    if( dim[0]>1 ) i += FieldT::Ghost::NM;
    if( dim[1]>1 ) j += FieldT::Ghost::NM;
    if( dim[2]>1 ) k += FieldT::Ghost::NM;
    return ijk2flat<FieldT>::value(dim,IndexTriplet(i,j,k),bcFlagX,bcFlagY,bcFlagZ);
  }

  //==================================================================

  /**
   *  @class  BCPoint 
   *  @author James C. Sutherland
   *  @date   June, 2008
   *
   *  @brief A simple class to impose BCs on a structured mesh using
   *  operators.  Intended for use when the BC is not located at the
   *  storage location for the field values and we need to use ghost
   *  values to achieve the desired BC.
   *
   *  The purpose of this class is to provide a more efficient way to
   *  impose BCs.  BCPoint objects can be constructed and stored in a
   *  container.  Much of the cost associated with imposing the BC is
   *  associated with extracting coefficients from the operator.  This
   *  is done at construction of a BCPoint object and is not repeated
   *  when the BC is imposed, leading to a more efficient BC
   *  implementation.  However, this operation is still applied at a
   *  point.  If we were to apply it at a set of points, we could
   *  potentially improve efficiency further.
   *
   *  The approach here is to set ghost values in a field to acheive a
   *  desired BC.  For example, a Dirichlet condition on a field can
   *  be achieved by using the interpolant operator and setting the
   *  field value at the ghost location so that the desired BC is
   *  achieved at the interpolated point.
   *
   *  @par Template Parameters
   *   <ul>
   *   <li> \b OpT Specifies the type of SpatialOperator that will be
   *   used to apply this BC.  This must define a few things:
   *     <ul>
   *     <li> \b DestFieldType The type for the destination field.
   *     <li> \b SrcFieldType  The type for the source field.
   *     </ul>
   *   </ul>
   *
   *  @todo Need a way of implementing time-varying BCs.  Perhaps this
   *  could be achieved via a functor passed into the constructor that
   *  would calculate the bc value given the coordinates and time?
   *  Then we would also need to know the coordinates and time,
   *  however.  The coordinates could be supplied at construction
   *  time, but the time would need to be an additional argument to
   *  BCPoint::operator().
   */
  template< typename OpT >
  class BCPoint
  {
    typedef typename OpT::SrcFieldType  SrcFieldT;
    typedef typename OpT::DestFieldType DestFieldT;

    const double bcVal_;
    const int ixf_;

    double ghostCoef_;
    typedef std::pair<int,double> IxValPair;
    std::vector<IxValPair> ixVals_;

  public:

    /**
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
     *  NOTE: the field to apply the BC on is not supplied here, since
     *  that would force its address to remain constant.  Rather, it
     *  is supplied to the BCPoint::operator(...) method.  Not quite
     *  as efficient, but less restrictive.  If we were guaranteed
     *  that the address of the field were fixed, then we could bind
     *  an iterator to the BC value insertion point.
     */
    BCPoint( const OpT& op,
	     const int i,
	     const int j,
	     const int k,
	     const std::vector<int>& dim,
	     const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ,
	     const double bcVal );

    ~BCPoint(){}

    /**
     *  Impose the BC on the supplied field.
     */
    inline void operator()( typename OpT::SrcFieldType& f ) const;
  };

  //==================================================================

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
   *   </ul>
   *
   */
  template< typename OpT >
  void assign_bc_point( const OpT& op,
			const int i,
			const int j,
			const int k,
			const std::vector<int>& dim,
			const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ,
			const double bcVal,
			typename OpT::SrcFieldType& f )
  {
    BCPoint<OpT> bcp(op,i,j,k,dim,bcFlagX,bcFlagY,bcFlagZ,bcVal);
    bcp(f);

//     typedef typename OpT::SrcFieldType  SrcFieldT;
//     typedef typename OpT::DestFieldType DestFieldT;

//     const int ixf = get_ghost_flat_ix<SrcFieldT >( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k );
//     int irow      = get_ghost_flat_ix<DestFieldT>( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k );

//     // NOTE: This will NOT work in the case where we have multiple ghost cells!
//     assert( OpT::SrcGhost::NM == OpT::SrcGhost::NP );
//     assert( OpT::SrcGhost::NM == 1 );

//     const typename OpT::MatrixRow row = op.get_row(irow);
//     typename       OpT::const_column_iterator icol =row.begin();
//     const typename OpT::const_column_iterator icole=row.end();

//     double prodsum=0.0; double ghostcoeff=0.0;
//     for( ; icol!=icole; ++icol ){
//       if( icol.index() == size_t(ixf) )
// 	ghostcoeff = *icol;
//       else
// 	prodsum += *icol *f[icol.index()];
//     }

//     assert( ghostcoeff != 0.0 );
//     f[ixf] = (bcVal - prodsum) / ghostcoeff;
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
			 OpT& op,
			 double& rhs )
  {
    typedef typename BCOpT::SrcFieldType  SrcFieldT;
    typedef typename BCOpT::DestFieldType DestFieldT;

    static int ncolMax = 10; // maximum number of nonzero columns
    struct BCInfo{ int ix; double coef; };
    BCInfo bcinfo[ncolMax];
    int nbcinfo = 0;

    // get the index into the field value at this point.
    IndexTriplet ijks(i,j,k);
    shift_ghost_ix<SrcFieldT,Dir>( dim, ijks );
    const int ixf  = get_ghost_flat_ix<SrcFieldT >( dim, bcFlagX, bcFlagY, bcFlagZ, ijks.i, ijks.j, ijks.k );

    IndexTriplet ijkd(i,j,k);
    shift_ghost_ix_dest<DestFieldT,Dir>( dim, ijkd );
    int irow       = get_ghost_flat_ix<DestFieldT>( dim, bcFlagX, bcFlagY, bcFlagZ, ijkd.i, ijkd.j, ijkd.k );

    const typename BCOpT::MatrixRow bcrow = bcOp.get_row(irow);
    typename       BCOpT::const_column_iterator icolbc = bcrow.begin();
    const typename BCOpT::const_column_iterator icolbce= bcrow.end();

    double ghostcoeff=0.0;
    BCInfo* bci = &bcinfo[0];
    for( ; icolbc!=icolbce; ++icolbc ){
      if( icolbc.index() == size_t(ixf) ){
	ghostcoeff = *icolbc;
      }
      else{
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
    const int ig = get_ghost_flat_ix<typename OpT::SrcFieldType >( dim, bcFlagX, bcFlagY, bcFlagZ, ijks.i, ijks.j, ijks.k );
    //    irow         = get_ghost_flat_ix<typename OpT::DestFieldType>( dim, bcFlagX, bcFlagY, bcFlagZ, ijkd.i, ijkd.j, ijkd.k );
    irow         = get_ghost_flat_ix<typename OpT::DestFieldType>( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k );

    typename OpT::MatrixRow row = op.get_row(irow);
    typename OpT::column_iterator icol=row.begin();
    const typename OpT::column_iterator icole=row.end();

    double Sg = 0.0;
    for( ; icol!=icole; ++icol ){
      if( icol.index() == size_t(ig) ){
	double& val = *icol;
	Sg = val;
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



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//                           IMPLEMENTATIONS
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



template< typename OpT >
BCPoint<OpT>::
BCPoint( const OpT& op,
	 const int i,
	 const int j,
	 const int k,
	 const std::vector<int>& dim,
	 const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ,
	 const double bcVal )
  : bcVal_( bcVal ),
    ixf_( get_ghost_flat_ix<SrcFieldT>( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k ) )
{
  const int irow = get_ghost_flat_ix<DestFieldT>( dim, bcFlagX, bcFlagY, bcFlagZ, i, j, k );

  // NOTE: This will NOT work in the case where we have multiple ghost cells!
  assert( OpT::SrcGhost::NM == OpT::SrcGhost::NP );
  assert( OpT::SrcGhost::NM == 1 );

  const typename OpT::MatrixRow row = op.get_row(irow);
  typename       OpT::const_column_iterator icol =row.begin();
  const typename OpT::const_column_iterator icole=row.end();
  ghostCoef_ = 0.0;
  for( ; icol!=icole; ++icol ){
    if( icol.index() == size_t(ixf_) )
      ghostCoef_ = *icol;
    else{
      ixVals_.push_back( std::make_pair(icol.index(),*icol) );
    }
  }
#ifndef NDEBUG
  if( ghostCoef_ == 0.0 ){
    std::cout << "Error in BCPoint." << std::endl
	      << "trying to set bc value: " << bcVal_ << endl
	      << "(i,j,k)=("<<i<<","<<j<<","<<k<<")"<<endl
	      << "ixf_ = " << ixf_ << endl
	      << "op coefs: ";
    for( typename OpT::const_column_iterator i=row.begin(); i!=icole; ++i ){
      cout << "  (" << i.index() << "," << *i << ")";
    }
    cout << endl;
  }
  assert( ghostCoef_ != 0.0 );
#endif
}

//--------------------------------------------------------------------

template< typename OpT >
void
BCPoint<OpT>::
operator()( SrcFieldT& f ) const
{
  double prodsum=0.0;
  for( std::vector<IxValPair>::const_iterator ix=ixVals_.begin(); ix!=ixVals_.end(); ++ix ){
    prodsum += ix->second * f[ix->first];
  }
  const double val = (bcVal_ - prodsum) / ghostCoef_;
  f[ixf_] = val;
}


}// namespace FVStaggered
}// namespace SpatialOps

//--------------------------------------------------------------------

#endif
