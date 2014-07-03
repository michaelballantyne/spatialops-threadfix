/*
 * Copyright (c) 2014 The University of Utah
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
#include <spatialops/structured/FVStaggeredLocationTypes.h>
#include <spatialops/structured/SpatialField.h>


/**
 *  \file FVStaggeredFieldTypes.h
 */

namespace SpatialOps{
namespace structured{

  /**
 *
 *  \addtogroup fieldtypes
 *  @{
 *
   *  \typedef typedef SpatialField< SVol > SVolField;
   *  \brief defines a volume field on the scalar volume.
   *
   *  \typedef typedef SpatialField< SSurfX > SSurfXField;
   *  \brief defines a x-surface field on the scalar volume
   *
   *  \typedef typedef SpatialField< SSurfY > SSurfYField;
   *  \brief defines a y-surface field on the scalar volume
   *
   *  \typedef typedef SpatialField< SSurfZ > SSurfZField;
   *  \brief defines a z-surface field on the scalar volume
   */
  typedef SpatialField< SVol   > SVolField;
  typedef SpatialField< SSurfX > SSurfXField;
  typedef SpatialField< SSurfY > SSurfYField;
  typedef SpatialField< SSurfZ > SSurfZField;


  /**
   *  \typedef typedef SpatialField< XVol > XVolField;
   *  \brief defines a volume field on the x-staggered volume
   *
   *  \typedef typedef SpatialField< XSurfX > XSurfXField;
   *  \brief defines a x-surface field on the x-staggered volume
   *
   *  \typedef typedef SpatialField< XSurfY > XSurfYField;
   *  \brief defines a y-surface field on the x-staggered volume
   *
   *  \typedef typedef SpatialField< XSurfZ > XSurfZField;
   *  \brief defines a z-surface field on the x-staggered volume
   */
  typedef SpatialField< XVol   > XVolField;
  typedef SpatialField< XSurfX > XSurfXField;
  typedef SpatialField< XSurfY > XSurfYField;
  typedef SpatialField< XSurfZ > XSurfZField;


  /**
   *  \typedef typedef SpatialField< YVol > YVolField;
   *  \brief defines a volume field on the y-staggered volume
   *
   *  \typedef typedef SpatialField< YSurfX > YSurfXField;
   *  \brief defines a x-surface field on the y-staggered volume
   *
   *  \typedef typedef SpatialField< YSurfY > YSurfYField;
   *  \brief defines a y-surface field on the y-staggered volume
   *
   *  \typedef typedef SpatialField< YSurfZ > YSurfZField;
   *  \brief defines a z-surface field on the y-staggered volume
   */
  typedef SpatialField< YVol   > YVolField;
  typedef SpatialField< YSurfX > YSurfXField;
  typedef SpatialField< YSurfY > YSurfYField;
  typedef SpatialField< YSurfZ > YSurfZField;


  /**
   *  \typedef typedef SpatialField< ZVol > ZVolField;
   *  \brief defines a volume field on the z-staggered volume
   *
   *  \typedef typedef SpatialField< ZSurfX > ZSurfXField;
   *  \brief defines a x-surface field on the z-staggered volume
   *
   *  \typedef typedef SpatialField< ZSurfY > ZSurfYField;
   *  \brief defines a y-surface field on the z-staggered volume
   *
   *  \typedef typedef SpatialField< ZSurfZ > ZSurfZField;
   *  \brief defines a z-surface field on the z-staggered volume
   */
  typedef SpatialField< ZVol   > ZVolField;
  typedef SpatialField< ZSurfX > ZSurfXField;
  typedef SpatialField< ZSurfY > ZSurfYField;
  typedef SpatialField< ZSurfZ > ZSurfZField;


  /**
   *  \typedef typedef SpatialField< SingleValue > SingleValueField;
   *  \brief defines a single value field
   *
   *  Warning: SingleValueFields should ONLY be built with MemoryWindows of size 1x1x1!
   */
  typedef SpatialField< SingleValue > SingleValueField;


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

/**
 *  @} // fieldtypes group
 */

}// namespace structured
}// namespace SpatialOps


#endif
