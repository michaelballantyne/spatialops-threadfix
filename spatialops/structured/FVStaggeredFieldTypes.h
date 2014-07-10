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

  /**
   *  \ingroup fieldtypes
   *  @{
   *
   *  \typedef typedef SVolField;
   *  \brief defines a volume field on the scalar volume.
   *
   *  \typedef typedef SSurfXField;
   *  \brief defines a x-surface field on the scalar volume
   *
   *  \typedef typedef SSurfYField;
   *  \brief defines a y-surface field on the scalar volume
   *
   *  \typedef typedef SSurfZField;
   *  \brief defines a z-surface field on the scalar volume
   */
  typedef SpatialField< SVol   > SVolField;
  typedef SpatialField< SSurfX > SSurfXField;
  typedef SpatialField< SSurfY > SSurfYField;
  typedef SpatialField< SSurfZ > SSurfZField;


  /**
   *  \typedef typedef XVolField;
   *  \brief defines a volume field on the x-staggered volume
   *
   *  \typedef typedef XSurfXField;
   *  \brief defines a x-surface field on the x-staggered volume
   *
   *  \typedef typedef XSurfYField;
   *  \brief defines a y-surface field on the x-staggered volume
   *
   *  \typedef typedef XSurfZField;
   *  \brief defines a z-surface field on the x-staggered volume
   */
  typedef SpatialField< XVol   > XVolField;
  typedef SpatialField< XSurfX > XSurfXField;
  typedef SpatialField< XSurfY > XSurfYField;
  typedef SpatialField< XSurfZ > XSurfZField;


  /**
   *  \typedef typedef YVolField;
   *  \brief defines a volume field on the y-staggered volume
   *
   *  \typedef typedef YSurfXField;
   *  \brief defines a x-surface field on the y-staggered volume
   *
   *  \typedef typedef YSurfYField;
   *  \brief defines a y-surface field on the y-staggered volume
   *
   *  \typedef typedef YSurfZField;
   *  \brief defines a z-surface field on the y-staggered volume
   */
  typedef SpatialField< YVol   > YVolField;
  typedef SpatialField< YSurfX > YSurfXField;
  typedef SpatialField< YSurfY > YSurfYField;
  typedef SpatialField< YSurfZ > YSurfZField;


  /**
   *  \typedef typedef ZVolField;
   *  \brief defines a volume field on the z-staggered volume
   *
   *  \typedef typedef ZSurfXField;
   *  \brief defines a x-surface field on the z-staggered volume
   *
   *  \typedef typedef ZSurfYField;
   *  \brief defines a y-surface field on the z-staggered volume
   *
   *  \typedef typedef ZSurfZField;
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
   *   - \link SpatialOps::SVolField  SVolField \endlink
   *   - \link SpatialOps::XVolField  XVolField \endlink
   *   - \link SpatialOps::YVolField  YVolField \endlink
   *   - \link SpatialOps::ZVolField  ZVolField \endlink
   *
   *  Specializations of this struct define the following typedefs:
   *   - \c XFace - the type of the field on the x-face
   *   - \c YFace - the type of the field on the yface
   *   - \c ZFace - the type of the field on the z-face
   *
   *  Example usage:
   *  \code{.cpp}
   *  typedef FaceTypes< CellT >::XFace XFaceT;
   *  typedef FaceTypes< CellT >::YFace YFaceT;
   *  typedef FaceTypes< CellT >::ZFace ZFaceT;
   *  \endcode
   *
   *  See also \ref VolType.
   */
  template< typename CellFieldT > struct FaceTypes;

  template<> struct FaceTypes<SVolField>
  {
    typedef SSurfXField XFace;
    typedef SSurfYField YFace;
    typedef SSurfZField ZFace;
  };

  template<> struct FaceTypes<XVolField>
  {
    typedef XSurfXField XFace;
    typedef XSurfYField YFace;
    typedef XSurfZField ZFace;
  };

  template<> struct FaceTypes<YVolField>
  {
    typedef YSurfXField XFace;
    typedef YSurfYField YFace;
    typedef YSurfZField ZFace;
  };

  template<> struct FaceTypes<ZVolField>
  {
    typedef ZSurfXField XFace;
    typedef ZSurfYField YFace;
    typedef ZSurfZField ZFace;
  };


  /**
   * @{
   *  \struct VolType
   *  \brief Define face field types in terms of a volume field type.
   *
   *  Class template specializations exist for the following field types:
   *   - \link SpatialOps::SSurfXField SSurfXField \endlink
   *   - \link SpatialOps::SSurfYField SSurfYField \endlink
   *   - \link SpatialOps::SSurfZField SSurfZField \endlink
   *   - \link SpatialOps::XSurfXField XSurfXField \endlink
   *   - \link SpatialOps::XSurfYField XSurfYField \endlink
   *   - \link SpatialOps::XSurfZField XSurfZField \endlink
   *   - \link SpatialOps::YSurfXField YSurfXField \endlink
   *   - \link SpatialOps::YSurfYField YSurfYField \endlink
   *   - \link SpatialOps::YSurfZField YSurfZField \endlink
   *   - \link SpatialOps::ZSurfXField ZSurfXField \endlink
   *   - \link SpatialOps::ZSurfYField ZSurfYField \endlink
   *   - \link SpatialOps::ZSurfZField ZSurfZField \endlink
   *
   *  Example usage:
   *  \code
   *  typedef VolType< FaceT       >::VolField FieldT;
   *  typedef VolType< SSurfZField >::VolField FieldT;
   *  \endcode
   *
   *  See also \ref FaceTypes.
   */
  template<typename FaceT> struct VolType;

  template<> struct VolType<SSurfXField>{ typedef SVolField VolField; };
  template<> struct VolType<SSurfYField>{ typedef SVolField VolField; };
  template<> struct VolType<SSurfZField>{ typedef SVolField VolField; };

  template<> struct VolType<XSurfXField>{ typedef XVolField VolField; };
  template<> struct VolType<XSurfYField>{ typedef XVolField VolField; };
  template<> struct VolType<XSurfZField>{ typedef XVolField VolField; };

  template<> struct VolType<YSurfXField>{ typedef YVolField VolField; };
  template<> struct VolType<YSurfYField>{ typedef YVolField VolField; };
  template<> struct VolType<YSurfZField>{ typedef YVolField VolField; };

  template<> struct VolType<ZSurfXField>{ typedef ZVolField VolField; };
  template<> struct VolType<ZSurfYField>{ typedef ZVolField VolField; };
  template<> struct VolType<ZSurfZField>{ typedef ZVolField VolField; };
  /**@}*/

/**
 *  @} // fieldtypes group
 */

}// namespace SpatialOps


#endif
