#ifndef FVStaggeredTypes_h
#define FVStaggeredTypes_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/structured/SpatialField.h>
#include <spatialops/SpatialOpsDefs.h>

#if defined(LINALG_UBLAS)
# include <spatialops/LinAlgUBlas.h>
  typedef SpatialOps::LinAlgUBlas LinAlg;
#elif defined(LINALG_TRILINOS)
# include <spatialops/LinAlgTrilinos.h>
  typedef SpatialOps::LinAlgTrilinos LinAlg;
#elif defined(LINALG_STENCIL)
  struct LinAlg{
    typedef int VecType;
    VecType& setup_vector( const int, double* ){ static VecType vt=0; return vt; }
  };
#else
#  error No LinAlg typedef was made!
#endif

/**
 *  \file FVStaggeredFieldTypes.h
 *
 *  \addtogroup structured
 *  @{
 *  \addtogroup fields
 *  @{
 *
 */

namespace SpatialOps{
namespace structured{

  // FaceDir: The direction relative to its volume field that this field is staggered.
  // StagLoc: The direction relative to the scalar volume field that this field's volume field is staggered.

  /**
   *  \struct SVol
   *  \brief Type traits for a scalar volume field
   *
   *  \struct SSurfX
   *  \brief Type traits for an x-surface field on a scalar volume
   *
   *  \struct SSurfY
   *  \brief Type traits for an y-surface field on a scalar volume
   *
   *  \struct SSurfZ
   *  \brief Type traits for an z-surface field on a scalar volume
   */
  struct SVol  { typedef NODIR FaceDir;  typedef NODIR StagLoc;  typedef NODIR RootLoc; };
  struct SSurfX{ typedef XDIR  FaceDir;  typedef NODIR StagLoc;  typedef XDIR  RootLoc; };
  struct SSurfY{ typedef YDIR  FaceDir;  typedef NODIR StagLoc;  typedef YDIR  RootLoc; };
  struct SSurfZ{ typedef ZDIR  FaceDir;  typedef NODIR StagLoc;  typedef ZDIR  RootLoc; };

  /**
   *  \struct XVol
   *  \brief Type traits for a x-staggered volume field
   *
   *  \struct XSurfX
   *  \brief Type traits for an x-surface field on a x-staggered volume
   *
   *  \struct XSurfY
   *  \brief Type traits for an y-surface field on a x-staggered volume
   *
   *  \struct XSurfZ
   *  \brief Type traits for an z-surface field on a x-staggered volume
   */
  struct XVol  { typedef NODIR FaceDir;  typedef XDIR StagLoc;  typedef NODIR RootLoc; }; //   typedef XDIR  RootLoc; };
  struct XSurfX{ typedef XDIR  FaceDir;  typedef XDIR StagLoc;  typedef XDIR  RootLoc; }; //   typedef NODIR RootLoc; };
  struct XSurfY{ typedef YDIR  FaceDir;  typedef XDIR StagLoc;  typedef YDIR  RootLoc; }; //   typedef YDIR  RootLoc; };
  struct XSurfZ{ typedef ZDIR  FaceDir;  typedef XDIR StagLoc;  typedef ZDIR  RootLoc; }; //   typedef ZDIR  RootLoc; };

  /**
   *  \struct YVol
   *  \brief Type traits for a y-staggered volume field
   *
   *  \struct YSurfX
   *  \brief Type traits for an x-surface field on a y-staggered volume
   *
   *  \struct YSurfY
   *  \brief Type traits for an y-surface field on a y-staggered volume
   *
   *  \struct YSurfZ
   *  \brief Type traits for an z-surface field on a y-staggered volume
   */
  struct YVol  { typedef NODIR FaceDir;  typedef YDIR  StagLoc;  typedef NODIR RootLoc; }; //  typedef YDIR  RootLoc; };
  struct YSurfX{ typedef XDIR  FaceDir;  typedef YDIR  StagLoc;  typedef XDIR  RootLoc; }; //  typedef XDIR  RootLoc; };
  struct YSurfY{ typedef YDIR  FaceDir;  typedef YDIR  StagLoc;  typedef YDIR  RootLoc; }; //  typedef NODIR RootLoc; };
  struct YSurfZ{ typedef ZDIR  FaceDir;  typedef YDIR  StagLoc;  typedef ZDIR  RootLoc; }; //  typedef ZDIR  RootLoc; };

  /**
   *  \struct ZVol
   *  \brief Type traits for a z-staggered volume field
   *
   *  \struct ZSurfX
   *  \brief Type traits for an x-surface field on a z-staggered volume
   *
   *  \struct ZSurfY
   *  \brief Type traits for an y-surface field on a z-staggered volume
   *
   *  \struct ZSurfZ
   *  \brief Type traits for an z-surface field on a z-staggered volume
   */
  struct ZVol  { typedef NODIR FaceDir;  typedef ZDIR  StagLoc;  typedef NODIR RootLoc; };  //  typedef ZDIR  RootLoc; };
  struct ZSurfX{ typedef XDIR  FaceDir;  typedef ZDIR  StagLoc;  typedef XDIR  RootLoc; };  //  typedef XDIR  RootLoc; };
  struct ZSurfY{ typedef YDIR  FaceDir;  typedef ZDIR  StagLoc;  typedef YDIR  RootLoc; };  //  typedef YDIR  RootLoc; };
  struct ZSurfZ{ typedef ZDIR  FaceDir;  typedef ZDIR  StagLoc;  typedef ZDIR  RootLoc; };  //  typedef NODIR RootLoc; };

  /**
   *  \struct DefaultGhost
   *  \brief define the default number of ghosts (1)
   *
   *  \struct NoGhost
   *  \brief define a type that has no ghost information
   */
  struct DefaultGhost{ enum{ NGHOST=1 }; };
  struct NoGhost     { enum{ NGHOST=0 }; };


  //-- Field Types --//

  /**
   *  \typedef typedef SpatialField< LinAlg, SVol,   DefaultGhost > SVolField;
   *  \brief defines a volume field on the scalar volume.
   *
   *  \typedef typedef SpatialField< LinAlg, SSurfX, DefaultGhost > SSurfXField;
   *  \brief defines a x-surface field on the scalar volume
   *
   *  \typedef typedef SpatialField< LinAlg, SSurfY, DefaultGhost > SSurfYField;
   *  \brief defines a y-surface field on the scalar volume
   *
   *  \typedef typedef SpatialField< LinAlg, SSurfZ, DefaultGhost > SSurfZField;
   *  \brief defines a z-surface field on the scalar volume
   */
  typedef SpatialField< LinAlg, SVol,   DefaultGhost > SVolField;
  typedef SpatialField< LinAlg, SSurfX, DefaultGhost > SSurfXField;
  typedef SpatialField< LinAlg, SSurfY, DefaultGhost > SSurfYField;
  typedef SpatialField< LinAlg, SSurfZ, DefaultGhost > SSurfZField;


  /**
   *  \typedef typedef SpatialField< LinAlg, XVol,   DefaultGhost > XVolField;
   *  \brief defines a volume field on the x-staggered volume
   *
   *  \typedef typedef SpatialField< LinAlg, XSurfX, DefaultGhost > XSurfXField;
   *  \brief defines a x-surface field on the x-staggered volume
   *
   *  \typedef typedef SpatialField< LinAlg, XSurfY, DefaultGhost > XSurfYField;
   *  \brief defines a y-surface field on the x-staggered volume
   *
   *  \typedef typedef SpatialField< LinAlg, XSurfZ, DefaultGhost > XSurfZField;
   *  \brief defines a z-surface field on the x-staggered volume
   */
  typedef SpatialField< LinAlg, XVol,   DefaultGhost > XVolField;
  typedef SpatialField< LinAlg, XSurfX, DefaultGhost > XSurfXField;
  typedef SpatialField< LinAlg, XSurfY, DefaultGhost > XSurfYField;
  typedef SpatialField< LinAlg, XSurfZ, DefaultGhost > XSurfZField;


  /**
   *  \typedef typedef SpatialField< LinAlg, YVol,   DefaultGhost > YVolField;
   *  \brief defines a volume field on the y-staggered volume
   *
   *  \typedef typedef SpatialField< LinAlg, YSurfX, DefaultGhost > YSurfXField;
   *  \brief defines a x-surface field on the y-staggered volume
   *
   *  \typedef typedef SpatialField< LinAlg, YSurfY, DefaultGhost > YSurfYField;
   *  \brief defines a y-surface field on the y-staggered volume
   *
   *  \typedef typedef SpatialField< LinAlg, YSurfZ, DefaultGhost > YSurfZField;
   *  \brief defines a z-surface field on the y-staggered volume
   */
  typedef SpatialField< LinAlg, YVol,   DefaultGhost > YVolField;
  typedef SpatialField< LinAlg, YSurfX, DefaultGhost > YSurfXField;
  typedef SpatialField< LinAlg, YSurfY, DefaultGhost > YSurfYField;
  typedef SpatialField< LinAlg, YSurfZ, DefaultGhost > YSurfZField;


  /**
   *  \typedef typedef SpatialField< LinAlg, ZVol,   DefaultGhost > ZVolField;
   *  \brief defines a volume field on the z-staggered volume
   *
   *  \typedef typedef SpatialField< LinAlg, ZSurfX, DefaultGhost > ZSurfXField;
   *  \brief defines a x-surface field on the z-staggered volume
   *
   *  \typedef typedef SpatialField< LinAlg, ZSurfY, DefaultGhost > ZSurfYField;
   *  \brief defines a y-surface field on the z-staggered volume
   *
   *  \typedef typedef SpatialField< LinAlg, ZSurfZ, DefaultGhost > ZSurfZField;
   *  \brief defines a z-surface field on the z-staggered volume
   */
  typedef SpatialField< LinAlg, ZVol,   DefaultGhost > ZVolField;
  typedef SpatialField< LinAlg, ZSurfX, DefaultGhost > ZSurfXField;
  typedef SpatialField< LinAlg, ZSurfY, DefaultGhost > ZSurfYField;
  typedef SpatialField< LinAlg, ZSurfZ, DefaultGhost > ZSurfZField;


  /**
   *  \struct FaceTypes
   *  \brief Define Face field types in terms of a cell field type.
   *
   *  Class template specializations exist for the following field types:
   *   - SVolField
   *   - XVolField
   *   - YVolField
   *   - ZVolField
   *
   *  Specializations of this struct define the following typedefs:
   *   - \c XFace - the type of the field on the x-face
   *   - \c YFace - the type of the field on the yface
   *   - \c ZFace - the type of the field on the z-face
   *
   *  Example usage:
   *  \code
   *  typedef FaceTypes< CellT >::XFace XFaceT;
   *  typedef FaceTypes< CellT >::YFace YFaceT;
   *  typedef FaceTypes< CellT >::ZFace ZFaceT;
   *  \endcode
   */
  template< typename CellFieldT > struct FaceTypes;

  template<> struct FaceTypes<SVolField>
  {
    typedef SpatialOps::structured::SSurfXField XFace;
    typedef SpatialOps::structured::SSurfYField YFace;
    typedef SpatialOps::structured::SSurfZField ZFace;
  };

  template<> struct FaceTypes<XVolField>
  {
    typedef SpatialOps::structured::XSurfXField XFace;
    typedef SpatialOps::structured::XSurfYField YFace;
    typedef SpatialOps::structured::XSurfZField ZFace;
  };

  template<> struct FaceTypes<YVolField>
  {
    typedef SpatialOps::structured::YSurfXField XFace;
    typedef SpatialOps::structured::YSurfYField YFace;
    typedef SpatialOps::structured::YSurfZField ZFace;
  };

  template<> struct FaceTypes<ZVolField>
  {
    typedef SpatialOps::structured::ZSurfXField XFace;
    typedef SpatialOps::structured::ZSurfYField YFace;
    typedef SpatialOps::structured::ZSurfZField ZFace;
  };


  /**
   *  \struct FaceTypes
   *  \brief Define Volume field types in terms of a face field type.
   *
   *  Class template specializations exist for the following field types:
   *   - SSurfXField
   *   - SSurfYField
   *   - SSurfZField
   *   - XSurfXField
   *   - XSurfYField
   *   - XSurfZField
   *   - YSurfXField
   *   - YSurfYField
   *   - YSurfZField
   *   - ZSurfXField
   *   - ZSurfYField
   *   - ZSurfZField
   *
   *  Example usage:
   *  \code
   *  typedef FaceTypes< FaceT       >::VolField FieldT;
   *  typedef FaceTypes< SSurfZField >::VolField FieldT;
   *  \endcode
   */
  template<typename FaceT> struct VolType;

  template<> struct VolType<SSurfXField>{ typedef SpatialOps::structured::SVolField VolField; };
  template<> struct VolType<SSurfYField>{ typedef SpatialOps::structured::SVolField VolField; };
  template<> struct VolType<SSurfZField>{ typedef SpatialOps::structured::SVolField VolField; };

  template<> struct VolType<XSurfXField>{ typedef SpatialOps::structured::XVolField VolField; };
  template<> struct VolType<XSurfYField>{ typedef SpatialOps::structured::XVolField VolField; };
  template<> struct VolType<XSurfZField>{ typedef SpatialOps::structured::XVolField VolField; };

  template<> struct VolType<YSurfXField>{ typedef SpatialOps::structured::YVolField VolField; };
  template<> struct VolType<YSurfYField>{ typedef SpatialOps::structured::YVolField VolField; };
  template<> struct VolType<YSurfZField>{ typedef SpatialOps::structured::YVolField VolField; };

  template<> struct VolType<ZSurfXField>{ typedef SpatialOps::structured::ZVolField VolField; };
  template<> struct VolType<ZSurfYField>{ typedef SpatialOps::structured::ZVolField VolField; };
  template<> struct VolType<ZSurfZField>{ typedef SpatialOps::structured::ZVolField VolField; };

  /**
   *  @}
   *  @}
   */


}// namespace structured
}// namespace SpatialOps

#endif
