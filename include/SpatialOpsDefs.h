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
  struct XDIR{};

  /**
   *  @struct YDIR
   *  @brief Defines a type for the y-direction.
   */
  struct YDIR{};

  /**
   *  @struct ZDIR
   *  @brief Defines a type for the z-direction.
   */
  struct ZDIR{};

  /**
   *  @struct NODIR
   *  @brief Defines a type to represent no direction
   */
  struct NODIR{};


  /** @} */  // end of Direction group.


  //==================================================================


  /**
   * @defgroup SideDefinitions
   *
   * @{
   */
  /** @struct SideMinus
   *  @brief  Defines a type for the (-) side of the domain.
   */
  struct SideMinus{};
  /** @struct SideMinus
   *  @brief  Defines a type for the (+) side of the domain.
   */
  struct SidePlus{};
  /** @struct SideMinus
   *  @brief  Defines a type for no specification of a side of the domain.
   */
  struct NoSide{};

  /** @} */  // end of Side Definition group


  //==================================================================

}

#endif
