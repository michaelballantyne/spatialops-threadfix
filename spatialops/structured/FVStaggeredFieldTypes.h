#ifndef FVStaggeredTypes_h
#define FVStaggeredTypes_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>
#include <spatialops/structured/SpatialField.h>
#include <spatialops/structured/IndexTriplet.h>

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
  //
  // Offset : The offset of this field relative to a scalar volume center.
  //
  // BCExtra: If a physical BC is present in the given direction, this
  //          typedef provides a way to determine the number of extra
  //          points in this field relative to a scalar volume field.

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
  struct SVol  { typedef NODIR FaceDir;  typedef IndexTriplet< 0, 0, 0> Offset;  typedef IndexTriplet<0,0,0>  BCExtra; };
  struct SSurfX{ typedef XDIR  FaceDir;  typedef IndexTriplet<-1, 0, 0> Offset;  typedef IndexTriplet<1,0,0>  BCExtra; };
  struct SSurfY{ typedef YDIR  FaceDir;  typedef IndexTriplet< 0,-1, 0> Offset;  typedef IndexTriplet<0,1,0>  BCExtra; };
  struct SSurfZ{ typedef ZDIR  FaceDir;  typedef IndexTriplet< 0, 0,-1> Offset;  typedef IndexTriplet<0,0,1>  BCExtra; };

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
  struct XVol  { typedef NODIR FaceDir;  typedef IndexTriplet<-1, 0, 0> Offset;  typedef IndexTriplet<1,0,0> BCExtra; };
  struct XSurfX{ typedef XDIR  FaceDir;  typedef IndexTriplet< 0, 0, 0> Offset;  typedef IndexTriplet<0,0,0> BCExtra; };
  struct XSurfY{ typedef YDIR  FaceDir;  typedef IndexTriplet<-1,-1, 0> Offset;  typedef IndexTriplet<0,1,0> BCExtra; };
  struct XSurfZ{ typedef ZDIR  FaceDir;  typedef IndexTriplet<-1, 0,-1> Offset;  typedef IndexTriplet<0,0,1> BCExtra; };

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
  struct YVol  { typedef NODIR FaceDir;  typedef IndexTriplet< 0,-1, 0> Offset;  typedef IndexTriplet<0,1,0> BCExtra; };
  struct YSurfX{ typedef XDIR  FaceDir;  typedef IndexTriplet<-1,-1, 0> Offset;  typedef IndexTriplet<1,0,0> BCExtra; };
  struct YSurfY{ typedef YDIR  FaceDir;  typedef IndexTriplet< 0, 0, 0> Offset;  typedef IndexTriplet<0,0,0> BCExtra; };
  struct YSurfZ{ typedef ZDIR  FaceDir;  typedef IndexTriplet< 0,-1,-1> Offset;  typedef IndexTriplet<0,0,1> BCExtra; };

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
  struct ZVol  { typedef NODIR FaceDir;  typedef IndexTriplet< 0, 0,-1> Offset;  typedef IndexTriplet<0,0,1> BCExtra; };
  struct ZSurfX{ typedef XDIR  FaceDir;  typedef IndexTriplet<-1, 0,-1> Offset;  typedef IndexTriplet<1,0,0> BCExtra; };
  struct ZSurfY{ typedef YDIR  FaceDir;  typedef IndexTriplet< 0,-1,-1> Offset;  typedef IndexTriplet<0,1,0> BCExtra; };
  struct ZSurfZ{ typedef ZDIR  FaceDir;  typedef IndexTriplet< 0, 0, 0> Offset;  typedef IndexTriplet<0,0,0> BCExtra; };

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
   *  \struct VolType
   *  \brief Define face field types in terms of a volume field type.
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
   *  typedef VolType< FaceT       >::VolField FieldT;
   *  typedef VolType< SSurfZField >::VolField FieldT;
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


}// namespace structured
}// namespace SpatialOps

/**
 *  @}
 *  @}
 */

#endif
