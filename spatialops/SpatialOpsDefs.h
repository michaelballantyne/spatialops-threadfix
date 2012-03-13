#ifndef UT_SpatialOpsDefs_h
#define UT_SpatialOpsDefs_h

#include <spatialops/SpatialOpsConfigure.h>

namespace SpatialOps{

  //==================================================================

  /**
   * @defgroup DirectionDefinitions
   * @{
   */

  /**
   *  @struct XDIR
   *  @brief Defines a type for the x-direction.
   */
  struct XDIR{ enum{value=0}; };

  /**
   *  @struct YDIR
   *  @brief Defines a type for the y-direction.
   */
  struct YDIR{ enum{value=1}; };

  /**
   *  @struct ZDIR
   *  @brief Defines a type for the z-direction.
   */
  struct ZDIR{ enum{value=2}; };

  /**
   *  @struct NODIR
   *  @brief Defines a type to represent no direction
   */
  struct NODIR{ enum{value=-10}; };

  /** @} */  // end of Direction group.




  /**
   *  @defgroup Operator Types
   *  @{
   */

  /**
   *  @struct Interpolant
   *  @brief  Defines a type for Interpolant operators.
   */
  struct Interpolant{};

  /**
   *  @struct Gradient
   *  @brief  Defines a type for Gradient operators.
   */
  struct Gradient{};

  /**
   *  @struct Divergence
   *  @brief  Defines a type for Divergence operators.
   */
  struct Divergence{};

  /**
   *  @struct Scratch
   *  @brief  Defines a type for Scratch operators.
   */
  struct Scratch{};

  /**
   *  @struct Filter
   *  @brief  Defines a type for Filter operators.
   */
  struct Filter{};

  /**
   *  @struct Restriction
   *  @brief  Defines a type for Restriction operators.
   */
  struct Restriction{};


  /**
   * @struct InterpolantX
   * @brief X-interpolant for use with FD operations whose src and dest fields are the same type
   */
  struct InterpolantX{ typedef XDIR DirT; };
  /**
   * @struct InterpolantY
   * @brief Y-interpolant for use with FD operations whose src and dest fields are the same type
   */
  struct InterpolantY{ typedef YDIR DirT; };

  /**
   * @struct InterpolantZ
   * @brief Z-interpolant for use with FD operations whose src and dest fields are the same type
   */
  struct InterpolantZ{ typedef ZDIR DirT; };

  /**
   * @struct GradientX
   * @brief X-interpolant for use with FD operations whose src and dest fields are the same type
   */
  struct GradientX{ typedef XDIR DirT; };
  /**
   * @struct GradientY
   * @brief Y-interpolant for use with FD operations whose src and dest fields are the same type
   */
  struct GradientY{ typedef YDIR DirT; };

  /**
   * @struct GradientZ
   * @brief Z-interpolant for use with FD operations whose src and dest fields are the same type
   */
  struct GradientZ{ typedef ZDIR DirT; };


/** @} */  // end of Operator Types group

  //==================================================================

}

#endif
