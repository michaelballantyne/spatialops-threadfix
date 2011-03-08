#ifndef SpatialOps_Structured_Stencil4_h
#define SpatialOps_Structured_Stencil4_h

namespace SpatialOps{
namespace structured{

  class MemoryWindow;
  class IntVec;

  /**
   *  \struct Stencil4
   *  \author James C. Sutherland
   *
   *  \brief Intended for use on 4-point stencils on structured, uniform meshes
   *
   *  \tparam OpT the type of operator (Interpolant, Gradient)
   *  \tparam SrcFieldT the type of field we apply the operator to
   *  \tparam DestFieldT the type of field produced by the operator
   *
   *  Examples of 4-point stencils:
   *   X-volume to y-surface or z-surface (advecting velocity)
   *   Y-volume to x-surface or z-surface (advecting velocity)
   *   Z-volume to x-surface or y-surface (advecting velocity)
   */
  template< typename OpT, typename SrcFieldT, typename DestFieldT >
  struct Stencil4
  {
    Stencil4( const double coef1,
              const double coef2,
              const double coef3,
              const double coef4 );

    void apply_to_field( const SrcFieldT& src, DestFieldT& dest ) const;

  private:
    const double coef1_, coef2_, coef3_, coef4_;
  };

  /**
   *  \struct Stencil4Helper
   *  \author James C. Sutherland
   *  \brief Provides information to customize the behavior of a Stencil4.
   */
  template< typename SrcFieldT, typename DestFieldT >
  struct Stencil4Helper;
  /*
  {
    Stencil4Helper( const MemoryWindow& wsrc,
                    const MemoryWindow& wdest );

    unsigned int src_offset_1() const;  ///< offset for the first source iterator
    unsigned int src_offset_2() const;  ///< offset for the second source iterator
    unsigned int src_offset_3() const;  ///< offset for the third source iterator
    unsigned int src_offset_4() const;  ///< offset for the fourth source iterator

    unsigned int dest_offset() const;  ///< offset for the destination iterator

    IntVec src_increment()  const;  ///< how far to increment the source iterators after each directional loop
    IntVec dest_increment() const;  ///< how far to increment the destination iterator after each directional loop

    IntVec low()  const; ///< the low index bounds for the x, y, and z loops
    IntVec high() const; ///< the high index bounds for the x, y, and z loops
  };
  */

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Structured_Stencil4_h
