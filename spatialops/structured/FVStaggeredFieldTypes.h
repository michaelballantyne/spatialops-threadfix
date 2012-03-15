/*
 * Copyright (c) 2011 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef FVStaggeredTypes_h
#define FVStaggeredTypes_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>
#include <spatialops/structured/SpatialField.h>
#include <spatialops/structured/IndexTriplet.h>


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
   *  \struct OneGhost
   *  \brief define the one ghost cell in each direction on each face
   *
   *  \struct NoGhost
   *  \brief define a type that has no ghost information
   */
  struct TwoGhost{ typedef IndexTriplet<2,2,2> NGhostMinus;  typedef IndexTriplet<2,2,2> NGhostPlus; };
  struct OneGhost{ typedef IndexTriplet<1,1,1> NGhostMinus;  typedef IndexTriplet<1,1,1> NGhostPlus; };
  struct NoGhost { typedef IndexTriplet<0,0,0> NGhostMinus;  typedef IndexTriplet<0,0,0> NGhostPlus; };


  //-- Field Types --//

  /**
   *  \typedef typedef SpatialField< SVol,   DefaultGhost > SVolField;
   *  \brief defines a volume field on the scalar volume.
   *
   *  \typedef typedef SpatialField< SSurfX, DefaultGhost > SSurfXField;
   *  \brief defines a x-surface field on the scalar volume
   *
   *  \typedef typedef SpatialField< SSurfY, DefaultGhost > SSurfYField;
   *  \brief defines a y-surface field on the scalar volume
   *
   *  \typedef typedef SpatialField< SSurfZ, DefaultGhost > SSurfZField;
   *  \brief defines a z-surface field on the scalar volume
   */
  typedef SpatialField< SVol,   OneGhost > SVolField;
  typedef SpatialField< SSurfX, OneGhost > SSurfXField;
  typedef SpatialField< SSurfY, OneGhost > SSurfYField;
  typedef SpatialField< SSurfZ, OneGhost > SSurfZField;


  /**
   *  \typedef typedef SpatialField< XVol,   DefaultGhost > XVolField;
   *  \brief defines a volume field on the x-staggered volume
   *
   *  \typedef typedef SpatialField< XSurfX, DefaultGhost > XSurfXField;
   *  \brief defines a x-surface field on the x-staggered volume
   *
   *  \typedef typedef SpatialField< XSurfY, DefaultGhost > XSurfYField;
   *  \brief defines a y-surface field on the x-staggered volume
   *
   *  \typedef typedef SpatialField< XSurfZ, DefaultGhost > XSurfZField;
   *  \brief defines a z-surface field on the x-staggered volume
   */
  typedef SpatialField< XVol,   OneGhost > XVolField;
  typedef SpatialField< XSurfX, OneGhost > XSurfXField;
  typedef SpatialField< XSurfY, OneGhost > XSurfYField;
  typedef SpatialField< XSurfZ, OneGhost > XSurfZField;


  /**
   *  \typedef typedef SpatialField< YVol,   DefaultGhost > YVolField;
   *  \brief defines a volume field on the y-staggered volume
   *
   *  \typedef typedef SpatialField< YSurfX, DefaultGhost > YSurfXField;
   *  \brief defines a x-surface field on the y-staggered volume
   *
   *  \typedef typedef SpatialField< YSurfY, DefaultGhost > YSurfYField;
   *  \brief defines a y-surface field on the y-staggered volume
   *
   *  \typedef typedef SpatialField< YSurfZ, DefaultGhost > YSurfZField;
   *  \brief defines a z-surface field on the y-staggered volume
   */
  typedef SpatialField< YVol,   OneGhost > YVolField;
  typedef SpatialField< YSurfX, OneGhost > YSurfXField;
  typedef SpatialField< YSurfY, OneGhost > YSurfYField;
  typedef SpatialField< YSurfZ, OneGhost > YSurfZField;


  /**
   *  \typedef typedef SpatialField< ZVol,   DefaultGhost > ZVolField;
   *  \brief defines a volume field on the z-staggered volume
   *
   *  \typedef typedef SpatialField< ZSurfX, DefaultGhost > ZSurfXField;
   *  \brief defines a x-surface field on the z-staggered volume
   *
   *  \typedef typedef SpatialField< ZSurfY, DefaultGhost > ZSurfYField;
   *  \brief defines a y-surface field on the z-staggered volume
   *
   *  \typedef typedef SpatialField< ZSurfZ, DefaultGhost > ZSurfZField;
   *  \brief defines a z-surface field on the z-staggered volume
   */
  typedef SpatialField< ZVol,   OneGhost > ZVolField;
  typedef SpatialField< ZSurfX, OneGhost > ZSurfXField;
  typedef SpatialField< ZSurfY, OneGhost > ZSurfYField;
  typedef SpatialField< ZSurfZ, OneGhost > ZSurfZField;


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
