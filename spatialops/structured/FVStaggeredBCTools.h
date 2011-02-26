#ifndef FVStaggeredBCTools_h
#define FVStaggeredBCTools_h

#include <spatialops/SpatialOpsConfigure.h>

#include <spatialops/structured/FVTools.h>

namespace SpatialOps{
namespace structured{

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

  struct ConstValVec{
    ConstValVec( const std::vector<double>& vec ) : val_( vec ){}
    inline const std::vector<double>& operator()() const{ return val_; }
  private:
    const std::vector<double> val_;
  };

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
    const IntVec point_;
    const BCEval bcEval_;

  public:

    /**
     *  @param point The IntVec specifying the location to apply
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
    BoundaryCondition( const IntVec point,
                       const BCEval bcEval );

    ~BoundaryCondition(){}

    /**
     *  Evaluates the boundary condition.
     *
     *  @param f The field that we want to set the BC on.
     */
    inline void operator()( FieldT& f ) const;
  };



  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //                         Implementation
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  //------------------------------------------------------------------

  template< typename FieldT, typename BCEval >
  BoundaryCondition<FieldT,BCEval>::
  BoundaryCondition( const IntVec point,
                     const BCEval bcEval )
    : point_( point ),
      bcEval_( bcEval )
  {}      

  //------------------------------------------------------------------

  template< typename FieldT, typename BCEval >
  void
  BoundaryCondition<FieldT,BCEval>::
  operator()( FieldT& f ) const
  {
    const unsigned int ix = f.window_without_ghost().flat_index( point_ );
    f[ix] = bcEval_();
  }

} // namespace structured
} // namespace SpatialOps


#endif  // Expr_BoundaryCondition_h
