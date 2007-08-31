#ifndef UT_SpatialOpsDefs_h
#define UT_SpatialOpsDefs_h

namespace SpatialOps{

  //==================================================================

  /**
   * @defgroup DirectionDefinitions
   *
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


  //==================================================================

}

#endif
