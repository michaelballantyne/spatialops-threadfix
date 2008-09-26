#ifndef FVStaggeredBCTools_h
#define FVStaggeredBCTools_h

#include <SpatialOperator.h>
#include <OperatorDatabase.h>
#include <FVStaggeredTools.h>

namespace SpatialOps{
namespace FVStaggered{

  /**
   *  @enum BCSide
   *  @brief For use with FV schemes.  Specifies the boundary
   *         condition on a face to the indicated side of the volume.
   */
  enum BCSide{
    X_PLUS_SIDE,
    X_MINUS_SIDE,
    Y_PLUS_SIDE,
    Y_MINUS_SIDE,
    Z_PLUS_SIDE,
    Z_MINUS_SIDE,
    NO_SHIFT
  };

  /**
   *  @brief A convenient way to implement constant valued boundary
   *  conditions.
   */
  struct ConstValEval{
    ConstValEval( const double val ) : val_(val) {}
    inline double operator()() const{ return val_; }
  private:
    const double val_;
  };

  /**
   *  @brief Obtain the flat index including ghost cells given the
   *  IndexTriplet based on the interior.
   *
   *  @param dim <code>vector<int></code> contiaining the interior
   *         domain extent in each direction.
   *
   *  @param bcFlagX Flag indicating whether this patch resides on a
   *         physical domain boundary in the (+x) direction.
   *
   *  @param bcFlagY Flag indicating whether this patch resides on a
   *         physical domain boundary in the (+y) direction.
   *
   *  @param bcFlagZ Flag indicating whether this patch resides on a
   *         physical domain boundary in the (+z) direction.
   *
   *  @param ijk The IndexTriplet (i,j,k) for the point, 0-based on
   *         domain interior (excludes ghost cells).
   */
  template< typename FieldT >
  int get_index_with_ghost( const std::vector<int>& dim,
                            const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ, 
                            IndexTriplet ijk );

  /**
   *  @brief Modify the IndexTriplet as appropriate to obtain the
   *         proper IndexTriplet for this field in the context of
   *         setting bcs on it using an operator.
   *
   *  @param dim <code>vector<int></code> contiaining the interior
   *         domain extent in each direction.
   *
   *  @param side The BCSide indicating which side (face) of the given
   *         cell the bc is to be applied on.
   *
   *  @param ijk The IndexTriplet (i,j,k) for the point, 0-based on
   *         domain interior (excludes ghost cells).
   */
  template<typename OpT, typename FieldT>
  IndexTriplet shift_to_ghost_ix( const std::vector<int>& dim,
                                  const BCSide side,
                                  IndexTriplet ijk );

  /**
   *  @class  BoundaryCondition
   *  @author James C. Sutherland
   *  @date   July, 2008
   *
   *  @brief Set a boundary condition directly on a variable at a single point.
   *
   *  @par Template Parameters
   *  <ul>
   *   <li> \b FieldT The type of field to set the bounary condition on.

   *   <li> \b BCEval An object of this type is supplied to calculate
   *        the boundary condition.  It takes no arguments and returns
   *        a single value which is the boundary condition at the
   *        point.  By using functors and bound functions (e.g. via
   *        boost::function and boost::lambda), you can achieve very
   *        complicated functors here using function composition.  For
   *        example, you could bind /f$g() = rho(t)*u(t)/f$ by binding
   *        a time functor to a functor for \f$rho(t)\f$ and
   *        \f$u(t)\f$ and then combining these to obtain /f$g() =
   *        rho(t)*u(t)/f$.

   *  </ul>
   *
   *  @par Design Considerations
   *  \li BoundaryCondition objects should be destroyed and rebuilt
   *      when mesh changes occur such as mesh point addition/removal.
   *      This is because the index for which they were originally
   *      built is no longer valid.  This convention also eliminates
   *      concern over operators becoming invalidated.
   */
  template< typename FieldT,
            typename BCEval >
  class BoundaryCondition
  {
    const int index_;
    const BCEval bcEval_;

  public:

    /**
     *  @param point The IndexTriplet specifying the location to apply
     *         this BC.  0-based on patch interior.
     *
     *  @param dim <code>vector<int></code> of extents in each
     *         direction (excluding ghost cells)
     *
     *  @param bcPlusX Flag indicating whether this patch resides on a
     *         physical domain boundary in the (+x) direction.
     *
     *  @param bcPlusY Flag indicating whether this patch resides on a
     *         physical domain boundary in the (+y) direction.
     *
     *  @param bcPlusZ Flag indicating whether this patch resides on a
     *         physical domain boundary in the (+z) direction.
     *
     *  @param bcEval A functor providing a method to evaluate the bc.
     *         It should have the following signature:
     *         <code>double()</code>.
     */
    BoundaryCondition( const IndexTriplet point,
                       const std::vector<int> dim,
                       const bool bcPlusX,
                       const bool bcPlusY,
                       const bool bcPlusZ,
                       const BCEval bcEval );

    ~BoundaryCondition(){}

    /**
     *  Evaluates the boundary condition.
     *
     *  @param f The field that we want to set the BC on.
     */
    inline void operator()( FieldT& f ) const;
  };

  //==================================================================

  /**
   *  @class BoundaryConditionOp
   *  @author James C. Sutherland
   *  @date July, 2008
   *
   *  @brief Imposes BCs on a field from a structured mesh using
   *         operators.  Intended for use when the BC is not located
   *         at the storage location for the field values and we need
   *         to use ghost values to achieve the desired BC.
   *
   *  The basic approach here is that we set the boundary condition on
   *  a field using an operator.  Supplying an interpolant operator
   *  results in a Dirichlet condition, while a gradient operator
   *  results in a Neumann condition.  Specifically, if \f$a_i\f$ are
   *  the coefficients for the operator at location "1" (the first
   *  interior mesh point), then the ghost value \f$f_0\f$ that
   *  produces the desired BC value, \f$f_{bc}\f$, is given as
   *  \f[
   *      f_0 = \frac{1}{a_0} \left( f_{bc} + \sum_{i=1}^{n} a_i f_i \right)
   * \f]
   *
   *  At construction, as much computation as possible is performed.
   *  This maximizes efficiency of the evaluation phase.
   *
   *  @par Template Parameters
   *  <ul>
   *  <li> \b OpT The type of operator that we are using to set the
   *       BCs.  This must conform to the interface of a
   *       SpatialOperator.
   *
   *  <li> \b BCEval The type for the functor being used to evaluate
   *       the boundary condition.  Suggestion: consider using
   *       Boost:Function here.  Any conforming interface
   *       <code>double()</code> should work, however.
   *  </ul>
   *
   *  @par Design Considerations
   *  \li BoundaryConditionOp objects should be destroyed and rebuilt
   *      when mesh changes occur such as mesh point addition/removal.
   *      This is because the index for which they were originally
   *      built is no longer valid.  This convention also eliminates
   *      concern over operators becoming invalidated.
   */
  template< typename OpT,
            typename BCEval >
  class BoundaryConditionOp
  {
    const BCEval bcEval_;
    const int index_;
    double ghostCoef_;

    typedef std::pair<int,double> IxValPair;
    std::vector<IxValPair> ixVals_;

    typedef typename OpT::SrcFieldType   SrcFieldT;
    typedef typename OpT::DestFieldType  DestFieldT;

  public:
    /**
     *  @param point The i,j,k location at which we want to specify
     *         the boundary condition (based on scalar cell center
     *         index)
     *
     *  @param side What side of the given point should be BC be applied to (+/-)?
     *
     *  @param eval The evalautor to obtain the bc value at this point.
     *
     *  @param soDatabase The database for spatial operators. An
     *         operator of type OpT will be extracted from this
     *         database.
     */
    BoundaryConditionOp( const std::vector<int>& dim,
                         const bool bcPlusX,
                         const bool bcPlusY,
                         const bool bcPlusZ,
                         const IndexTriplet point,
                         const BCSide side,
                         const BCEval bceval,
                         const OperatorDatabase& soDatabase );

    ~BoundaryConditionOp(){}

    /**
     *  Impose the boundary condition on the supplied field.
     */
    void operator()( SrcFieldT& f ) const;
  };



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
   *  @param ijk The IndexTriplet for the cell we want to apply the BC
   *         to. Indices are 0-based on patch interior.
   *
   *  @param dim A vector containing the number of cells in each
   *         coordinate direction.  This is a three-component vector.
   *
   *  @param bcFlagX A boolean flag to indicate if this patch is on a
   *         +x side physical boundary.  If so, then it is assumed
   *         that there is an extra face on that side of the domain,
   *         and face variable dimensions will be modified
   *         accordingly.
   *
   *  @param bcFlagY A boolean flag to indicate if this patch is on a
   *         +y side physical boundary.  If so, then it is assumed
   *         that there is an extra face on that side of the domain,
   *         and face variable dimensions will be modified
   *         accordingly.
   *
   *  @param bcFlagZ A boolean flag to indicate if this patch is on a
   *         +z side physical boundary.  If so, then it is assumed
   *         that there is an extra face on that side of the domain,
   *         and face variable dimensions will be modified
   *         accordingly.
   *
   *  @param bcVal The value for the boundary condition to set.
   *
   *  @par Template Parameters
   *   <ul>
   *   <li> \b BCOpT Specifies the type of SpatialOperator that will be
   *        used to apply this BC.  
   *   <li> \b OpT Specifies the type of SpatialOperator that will be
   *        imprinting with this BC.
   *   </ul>
   *
   *  @todo Need to allow for time-varying BCs. This effects the RHS
   *  but not the LHS.  This needs to hook into the expression somewhere.
   */
  template< typename BCOpT, typename OpT >
  void imprint_bc_on_op( const BCOpT& bcOp,
                         const IndexTriplet ijk,
                         const std::vector<int>& dim,
                         const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ,
                         const double bcVal,
                         const BCSide side,
                         OpT& op,
                         double& rhs );



  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //                         Implementation
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  //------------------------------------------------------------------

  template< typename FieldT, typename BCEval >
  BoundaryCondition<FieldT,BCEval>::
  BoundaryCondition( const IndexTriplet point,
                     const std::vector<int> dim,
                     const bool bcPlusX,
                     const bool bcPlusY,
                     const bool bcPlusZ,
                     const BCEval bcEval )
    : index_( get_index_with_ghost<FieldT>( dim, bcPlusX, bcPlusY, bcPlusZ, point ) ),
      bcEval_( bcEval )
  {}      

  //------------------------------------------------------------------

  template< typename FieldT, typename BCEval >
  void
  BoundaryCondition<FieldT,BCEval>::
  operator()( FieldT& f ) const
  {
    f[index_] = bcEval_();
  }

  //------------------------------------------------------------------


  //==================================================================


  //------------------------------------------------------------------

  template< typename OpT, typename BCEval >
  BoundaryConditionOp<OpT,BCEval>::
  BoundaryConditionOp( const std::vector<int>& dim,
                       const bool bcPlusX,
                       const bool bcPlusY,
                       const bool bcPlusZ,
                       const IndexTriplet point,
                       const BCSide side,
                       const BCEval bcEval,
                       const OperatorDatabase& soDatabase )
    : bcEval_( bcEval ),
      index_( get_index_with_ghost<SrcFieldT>( dim, bcPlusX, bcPlusY, bcPlusZ,
                                               shift_to_ghost_ix<OpT,SrcFieldT>(dim,side,point) ) )
  {
    const OpT* op = soDatabase.retrieve_operator<OpT>();
    const int irow = get_index_with_ghost<DestFieldT>( dim, bcPlusX, bcPlusY, bcPlusZ,
                                                       shift_to_ghost_ix<OpT,DestFieldT>(dim,side,point) );

    const typename OpT::MatrixRow row = op->get_row(irow);
    typename       OpT::const_column_iterator icol =row.begin();
    const typename OpT::const_column_iterator icole=row.end();
    ghostCoef_ = 0.0;
    for( ; icol!=icole; ++icol ){
      if( icol.index() == size_t(index_) )
        ghostCoef_ = *icol;
      else{
        ixVals_.push_back( std::make_pair(icol.index(),*icol) );
      }
    }
#   ifndef NDEBUG
    if( ghostCoef_ == 0.0 ){
      std::cout << "Error in BoundaryConditionOp." << std::endl
                << "(i,j,k)=("<<point.i<<","<< point.j <<","<< point.k <<")"<<endl
                << "index_ = " << index_ << endl
                << "row = " << irow << endl
                << "op coefs: ";
      for( typename OpT::const_column_iterator icol=row.begin(); icol!=icole; ++icol ){
        cout << "  (" << icol.index() << "," << *icol << ")";
      }
      cout << endl;
    }
    assert( ghostCoef_ != 0.0 );
#   endif
  }

  //------------------------------------------------------------------

  template< typename OpT, typename BCEval >
  void
  BoundaryConditionOp<OpT,BCEval>::
  operator()( SrcFieldT& f ) const
  {
    double prodsum=0.0;
    for( std::vector<IxValPair>::const_iterator ix=ixVals_.begin(); ix!=ixVals_.end(); ++ix ){
      prodsum += ix->second * f[ix->first];
    }
    const double val = ( bcEval_() - prodsum) / ghostCoef_;
    f[index_] = val;
  }

  //--------------------------------------------------------------------

  template< typename FieldT >
  int get_index_with_ghost( const std::vector<int>& dim,
                            const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ, 
                            IndexTriplet index )
  {
    if( dim[0]>1 )  index.i += FieldT::Ghost::NM;
    if( dim[1]>1 )  index.j += FieldT::Ghost::NM;
    if( dim[2]>1 )  index.k += FieldT::Ghost::NM;
    return ijk2flat<FieldT>::value( dim, index, bcFlagX, bcFlagY, bcFlagZ );
  }

  template<typename OpT, typename FieldT>
  IndexTriplet shift_to_ghost_ix( const std::vector<int>& dim,
                                  const BCSide side,
                                  IndexTriplet ijk )
  {
    if( IsSameType<typename OpT::SrcFieldType,FieldT>::result ){
      switch(side){
      case X_MINUS_SIDE: if(dim[0]>1) --ijk.i; break;
      case Y_MINUS_SIDE: if(dim[1]>1) --ijk.j; break;
      case Z_MINUS_SIDE: if(dim[2]>1) --ijk.k; break;
      case X_PLUS_SIDE : if(dim[0]>1) ++ijk.i; break;
      case Y_PLUS_SIDE : if(dim[1]>1) ++ijk.j; break;
      case Z_PLUS_SIDE : if(dim[2]>1) ++ijk.k; break;
      case NO_SHIFT: assert(1); break;
      }
    }
    else{
      switch(side){
      case X_PLUS_SIDE : if(dim[0]>1) ++ijk.i; break;
      case Y_PLUS_SIDE : if(dim[1]>1) ++ijk.j; break;
      case Z_PLUS_SIDE : if(dim[2]>1) ++ijk.k; break;
      case X_MINUS_SIDE: break;
      case Y_MINUS_SIDE: break;
      case Z_MINUS_SIDE: break;
      case NO_SHIFT: assert(1); break;
      }
    }
    return ijk;
  }

  template<> IndexTriplet
  shift_to_ghost_ix<GradXVolXSurfX,XVolField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case X_MINUS_SIDE: if(dim[0]>1) --ijk.i; break;
    case X_PLUS_SIDE : if(dim[0]>1) ijk.i+=2; break;
    // error cases:
    case Y_MINUS_SIDE:  case Y_PLUS_SIDE:
    case Z_MINUS_SIDE:  case Z_PLUS_SIDE:
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }
  template<> IndexTriplet
  shift_to_ghost_ix<GradXVolXSurfX,XSurfXField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case X_MINUS_SIDE: if(dim[0]>1) --ijk.i; break;
    case X_PLUS_SIDE : if(dim[0]>1) ++ijk.i; break;
    // error cases:
    case Y_MINUS_SIDE:  case Y_PLUS_SIDE: 
    case Z_MINUS_SIDE:  case Z_PLUS_SIDE: 
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }
  template<> IndexTriplet
  shift_to_ghost_ix<InterpXVolXSurfX,XVolField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case X_MINUS_SIDE: if(dim[0]>1) --ijk.i; break;
    case X_PLUS_SIDE : if(dim[0]>1) ijk.i+=2; break;
    // error cases:
    case Y_MINUS_SIDE:  case Y_PLUS_SIDE:
    case Z_MINUS_SIDE:  case Z_PLUS_SIDE:
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }
  template<> IndexTriplet
  shift_to_ghost_ix<InterpXVolXSurfX,XSurfXField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case X_MINUS_SIDE: if(dim[0]>1) --ijk.i; break;
    case X_PLUS_SIDE : if(dim[0]>1) ++ijk.i; break;
    // error cases:
    case Y_MINUS_SIDE:  case Y_PLUS_SIDE: 
    case Z_MINUS_SIDE:  case Z_PLUS_SIDE: 
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }
  template<> IndexTriplet
  shift_to_ghost_ix<GradYVolYSurfY,YVolField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case Y_MINUS_SIDE: if(dim[1]>1) --ijk.j; break;
    case Y_PLUS_SIDE : if(dim[1]>1) ijk.j+=2; break;
    // error cases:
    case X_MINUS_SIDE:  case X_PLUS_SIDE:
    case Z_MINUS_SIDE:  case Z_PLUS_SIDE:
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }
  template<> IndexTriplet
  shift_to_ghost_ix<GradYVolYSurfY,YSurfYField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case Y_MINUS_SIDE: if(dim[1]>1) --ijk.j; break;
    case Y_PLUS_SIDE : if(dim[1]>1) ++ijk.j; break;
    // error cases:
    case X_MINUS_SIDE:  case X_PLUS_SIDE: 
    case Z_MINUS_SIDE:  case Z_PLUS_SIDE: 
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }
  template<> IndexTriplet
  shift_to_ghost_ix<InterpYVolYSurfY,YVolField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case Y_MINUS_SIDE: if(dim[1]>1) --ijk.j; break;
    case Y_PLUS_SIDE : if(dim[1]>1) ijk.j+=2; break;
    // error cases:
    case X_MINUS_SIDE:  case X_PLUS_SIDE:
    case Z_MINUS_SIDE:  case Z_PLUS_SIDE:
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }
  template<> IndexTriplet
  shift_to_ghost_ix<InterpYVolYSurfY,YSurfYField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case Y_MINUS_SIDE: if(dim[1]>1) --ijk.j; break;
    case Y_PLUS_SIDE : if(dim[1]>1) ++ijk.j; break;
    // error cases:
    case X_MINUS_SIDE:  case X_PLUS_SIDE: 
    case Z_MINUS_SIDE:  case Z_PLUS_SIDE: 
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }

  template<> IndexTriplet
  shift_to_ghost_ix<GradZVolZSurfZ,ZVolField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case Z_MINUS_SIDE: if(dim[2]>1) --ijk.k; break;
    case Z_PLUS_SIDE : if(dim[2]>1) ijk.k+=2; break;
    // error cases:
    case X_MINUS_SIDE:  case X_PLUS_SIDE:
    case Y_MINUS_SIDE:  case Y_PLUS_SIDE:
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }
  template<> IndexTriplet
  shift_to_ghost_ix<GradZVolZSurfZ,ZSurfZField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case Z_MINUS_SIDE: if(dim[2]>1) --ijk.k; break;
    case Z_PLUS_SIDE : if(dim[2]>1) ++ijk.k; break;
    // error cases:
    case X_MINUS_SIDE:  case X_PLUS_SIDE: 
    case Y_MINUS_SIDE:  case Y_PLUS_SIDE: 
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }
  template<> IndexTriplet
  shift_to_ghost_ix<InterpZVolZSurfZ,ZVolField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case Z_MINUS_SIDE: if(dim[2]>1) --ijk.k; break;
    case Z_PLUS_SIDE : if(dim[2]>1) ijk.k+=2; break;
    // error cases:
    case X_MINUS_SIDE:  case X_PLUS_SIDE:
    case Y_MINUS_SIDE:  case Y_PLUS_SIDE:
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }
  template<> IndexTriplet
  shift_to_ghost_ix<InterpZVolZSurfZ,ZSurfZField>( const std::vector<int>& dim, const BCSide side, IndexTriplet ijk )
  {
    switch(side){
    case Z_MINUS_SIDE: if(dim[2]>1) --ijk.k; break;
    case Z_PLUS_SIDE : if(dim[2]>1) ++ijk.k; break;
    // error cases:
    case X_MINUS_SIDE:  case X_PLUS_SIDE: 
    case Y_MINUS_SIDE:  case Y_PLUS_SIDE: 
    case NO_SHIFT: assert(1); break;
    }
    return ijk;
  }


  //------------------------------------------------------------------

  template< typename BCOpT, typename OpT >
  void imprint_bc_on_op( const BCOpT& bcOp,
                         const IndexTriplet ijk,
                         const std::vector<int>& dim,
                         const bool bcFlagX, const bool bcFlagY, const bool bcFlagZ,
                         const double bcVal,
                         const BCSide side,
                         OpT& op,
                         double& rhs )
  {
    static int ncolMax = 10; // maximum number of nonzero columns
    struct BCInfo{ int ix; double coef; };
    BCInfo bcinfo[ncolMax];
    int nbcinfo = 0;

    // get the index into the field value at this point.
    IndexTriplet ijks( shift_to_ghost_ix<BCOpT,typename BCOpT::SrcFieldType >(dim,side,ijk) );
    IndexTriplet ijkd( shift_to_ghost_ix<BCOpT,typename BCOpT::DestFieldType>(dim,side,ijk) );
    const int ixf = get_index_with_ghost<typename BCOpT::SrcFieldType >( dim, bcFlagX, bcFlagY, bcFlagZ, ijks );
    int irow      = get_index_with_ghost<typename BCOpT::DestFieldType>( dim, bcFlagX, bcFlagY, bcFlagZ, ijkd );

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
    BOOST_STATIC_ASSERT( bool( IsSameType<typename   OpT::SrcFieldType, typename OpT::DestFieldType>::result ) );
    BOOST_STATIC_ASSERT( bool( IsSameType<typename BCOpT::SrcFieldType, typename OpT::SrcFieldType >::result ) );

    //
    // now set the operator value.  We must potentially alter each coefficient in this row.
    //
    const int ig = ixf;// = get_index_with_ghost<typename OpT::SrcFieldType >( dim, bcFlagX, bcFlagY, bcFlagZ, ijks );
    irow         = get_index_with_ghost<typename OpT::DestFieldType>( dim, bcFlagX, bcFlagY, bcFlagZ, ijk );

    typename OpT::MatrixRow row = op.get_row(irow);
    typename OpT::column_iterator icol=row.begin();
    const typename OpT::column_iterator icole=row.end();

    double Sg = 0.0;
    for( ; icol!=icole; ++icol ){
      if( icol.index() == size_t(ig) ){
        Sg = *icol;
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

} // namespace FVStaggered
} // namespace SpatialOps


#endif  // Expr_BoundaryCondition_h
