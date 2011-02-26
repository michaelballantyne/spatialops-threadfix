#ifndef SpatialOps_FVStaggeredStencilBCOp_h
#define SpatialOps_FVStaggeredStencilBCOp_h

#include <spatialops/structured/MemoryWindow.h>
#include <spatialops/structured/FVStaggeredBCTools.h>

namespace SpatialOps{
namespace structured{

  template< typename OpT,
            typename BCEval >
  class BoundaryConditionOp
  {
    const BCEval bcEval_;
    IntVec apoint_, bpoint_;
    double ca_, cb_;

    typedef typename OpT::SrcFieldType   SrcFieldT;

    BoundaryConditionOp& operator=( const BoundaryConditionOp& ); // no assignment
    BoundaryConditionOp(); // no default constructor

  public:
    
    /**
     *  Expose the bcevaluator type.
     */
    typedef BCEval BCEvalT;
    
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
    BoundaryConditionOp( const IntVec& dim,
                         const bool bcPlusX,
                         const bool bcPlusY,
                         const bool bcPlusZ,
                         const IntVec point,
                         const BCSide side,
                         const BCEval bceval,
                         const OperatorDatabase& soDatabase );

    ~BoundaryConditionOp(){}

    /**
     *  Impose the boundary condition on the supplied field.
     */
    void operator()( SrcFieldT& f ) const;

    void operator()( std::vector<SrcFieldT*>& f ) const;

  }; // class BoundaryConditionOp



  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //                         Implementation
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  template< typename OpT, typename BCEval >
  BoundaryConditionOp<OpT,BCEval>::
  BoundaryConditionOp( const IntVec& dim,
                       const bool bcPlusX,
                       const bool bcPlusY,
                       const bool bcPlusZ,
                       const IntVec point,
                       const BCSide side,
                       const BCEval bceval,
                       const OperatorDatabase& soDatabase )
    : bcEval_( bceval )
  {
    // let phi_a be the ghost value, phi_b be the internal value, and phi_bc be the boundary condition,
    //   phi_bc = a*phi_a + b*phi_b
    // then
    //   phi_a = (phi_bc - b*phi_b) / a
    //
    IntVec iashift(0,0,0), ibshift(0,0,0);
    const OpT* const op = soDatabase.retrieve_operator<OpT>();
    switch(side){
    case X_MINUS_SIDE:
      ca_ = op->get_minus_coef();
      cb_ = op->get_plus_coef();
      iashift[0] = -1;
      ibshift[0] = 0;
      break;
    case Y_MINUS_SIDE:
      ca_ = op->get_minus_coef();
      cb_ = op->get_plus_coef();
      iashift[1] = -1;
      ibshift[1] = 0;
      break;
    case Z_MINUS_SIDE:
      ca_ = op->get_minus_coef();
      cb_ = op->get_plus_coef();
      iashift[2] = -1;
      ibshift[2] = 0;
      break;
    case X_PLUS_SIDE:
      cb_ = op->get_minus_coef();
      ca_ = op->get_plus_coef();
      iashift[0] = 0;
      ibshift[0] = -1;
      break;
    case Y_PLUS_SIDE:
      cb_ = op->get_minus_coef();
      ca_ = op->get_plus_coef();
      iashift[1] = 0;
      ibshift[1] = -1;
      break;
    case Z_PLUS_SIDE:
      cb_ = op->get_minus_coef();
      ca_ = op->get_plus_coef();
      iashift[2] = 0;
      ibshift[2] = -1;
      break;
    default:
      throw std::runtime_error("Invalid BC face specification");
    }
    apoint_ = point + iashift;
    bpoint_ = point + ibshift;
  }

  //------------------------------------------------------------------

  template< typename OpT, typename BCEval >
  void
  BoundaryConditionOp<OpT,BCEval>::
  operator()( SrcFieldT& f ) const
  {
    const unsigned int ia = f.window_without_ghost().flat_index(apoint_);
    const unsigned int ib = f.window_without_ghost().flat_index(bpoint_);
    //    std::cout << apoint_ << ", " << bpoint_ << " : " << ia << "," << ib << " : " << f.window_without_ghost() << " : "<< f.window_with_ghost() << endl;
    f[ia] = ( bcEval_() - cb_*f[ib] ) / ca_;
  }

  //------------------------------------------------------------------

} // namespace 
} // namespace SpatialOps

#endif SpatialOps_FVStaggeredStencilBCOp_h
