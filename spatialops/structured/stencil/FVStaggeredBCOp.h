#ifndef SpatialOps_FVStaggeredStencilBCOp_h
#define SpatialOps_FVStaggeredStencilBCOp_h

#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/IndexTriplet.h>

namespace SpatialOps{
namespace structured{


  /**
   * \enum BCSide
   * \brief Allows identification of whether we are setting the BC
   *        on the right or left side when using an operator.
   */
  enum BCSide{
    MINUS_SIDE,  ///< Minus side
    PLUS_SIDE    ///< Plus side
  };

  /**
   * \class BoundaryConditionOp
   * \brief Provides a simple interface to set a boundary condition via an operator.
   * \tparam OpT the operator type for use in setting a BC
   * \tparam BCEval a functor for obtaining the BC value
   *
   *  NOTE: The BoundaryConditionOp class should only be used with
   *        operators that involve the scalar volume.
   */
  template< typename OpT,
            typename BCEval >
  class BoundaryConditionOp
  {
    typedef typename OpT::SrcFieldType    SrcFieldT;
    typedef typename OpT::DestFieldType   DestFieldT;

    typedef typename  SrcFieldT::Location::Offset  SO;
    typedef typename DestFieldT::Location::Offset  DO;

    typedef typename Subtract<SO,DO>::result SOMinusDO;

    typedef typename UnitTriplet< typename GetNonzeroDir<SOMinusDO>::DirT >::type UnitVec;

    typedef typename Multiply< DO, UnitVec >::result  S1Shift;
    typedef typename Add< S1Shift, UnitVec >::result  S2Shift;

    const BCEval bcEval_;  ///< functor to set the value of the BC
    const IntVec apoint_;  ///< the index for the value in the source field we will set
    const IntVec bpoint_;  ///< the index for the value in the source field we use to obtain the value we want to set.
    double ca_, cb_;       ///< high and low coefficients for the operator

    BoundaryConditionOp& operator=( const BoundaryConditionOp& ); // no assignment
    BoundaryConditionOp();                                        // no default constructor

  public:

    typedef BCEval BCEvalT;  ///< Expose the BCEval type.
    
    /**
     *  \param destIndex The i,j,k location at which we want to specify
     *         the boundary condition.  This is indexed 0-based on
     *         the interior (neglecting ghost cells), and refers to
     *         the index in the "destination" field of the operator.
     *
     *  \param side The side of the cell (MINUS_SIDE or PLUS_SIDE) that
     *         this BC is to be applied on.
     *
     *  \param eval The evaluator to obtain the bc value at this point.
     *
     *  \param opdb The database for spatial operators. An operator of
     *         type OpT will be extracted from this database.
     */
    BoundaryConditionOp( const IntVec& destIndex,
                         const BCSide side,
                         const BCEval bceval,
                         const OperatorDatabase& opdb );

    ~BoundaryConditionOp(){}

    /**
     *  \brief Impose the boundary condition on the supplied field.
     */
    void operator()( SrcFieldT& f ) const;

    /**
     *  \brief Impose the boundary condition on the supplied fields.
     */
    void operator()( std::vector<SrcFieldT*>& f ) const;

  }; // class BoundaryConditionOp


  // ================================================================
  //
  //                            Implementation
  //
  // ================================================================

  template< typename OpT, typename BCEval >
  BoundaryConditionOp<OpT,BCEval>::
  BoundaryConditionOp( const IntVec& destPoint,
                       const BCSide side,
                       const BCEval bceval,
                       const OperatorDatabase& soDatabase )
      : bcEval_( bceval ),
        apoint_( destPoint + ( (side==MINUS_SIDE) ? S1Shift::int_vec() : S2Shift::int_vec() ) ),
        bpoint_( destPoint + ( (side==MINUS_SIDE) ? S2Shift::int_vec() : S1Shift::int_vec() ) )
  {
    // let phi_a be the ghost value, phi_b be the internal value, and phi_bc be the boundary condition,
    //   phi_bc = a*phi_a + b*phi_b
    // then
    //   phi_a = (phi_bc - b*phi_b) / a
    //
    const OpT* const op = soDatabase.retrieve_operator<OpT>();
    ca_ = (side==MINUS_SIDE ? op->get_minus_coef() : op->get_plus_coef()  );
    cb_ = (side==MINUS_SIDE ? op->get_plus_coef()  : op->get_minus_coef() );
  }

  //------------------------------------------------------------------

  template< typename OpT, typename BCEval >
  void
  BoundaryConditionOp<OpT,BCEval>::
  operator()( SrcFieldT& f ) const
  {
    // jcs: this is not very efficient (indexing slowness) but I
    //      am not sure that we can do any better at this point.
    const unsigned int ia = f.window_without_ghost().flat_index(apoint_);
    const unsigned int ib = f.window_without_ghost().flat_index(bpoint_);
    f[ia] = ( bcEval_() - cb_*f[ib] ) / ca_;
  }

  //------------------------------------------------------------------

} // namespace 
} // namespace SpatialOps

#endif // SpatialOps_FVStaggeredStencilBCOp_h
