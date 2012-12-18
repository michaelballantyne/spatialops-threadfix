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

#ifndef SpatialOps_Reference_Stencil_h
#define SpatialOps_Reference_Stencil_h

#include <spatialops/structured/IndexTriplet.h>

namespace SpatialOps {
namespace structured {

    /**
     * \brief apply reference (2-point) stencil on given fields with given coefficients
     * \tparam SrcType - the type of field the operator is applied to
     * \tparam DestType - the type of field the operator produces
     * \param coefLo the coefficient to multiply the (-) side field by
     * \param coefHi the coefficient to multiply the (+) side field by
     * \param src the field that the operator is applied to
     * \param dest the resulting field.
     */
    template<typename SrcType, typename DestType>
    inline void ref_stencil2_apply_to_field( const double coefLo,
                                             const double coefHi,
                                             const SrcType  & src,
                                                   DestType & dest );

    /**
     * \brief apply reference (4-point) stencil on given fields with given coefficients
     * \tparam SrcType - the type of field the operator is applied to
     * \tparam DestType - the type of field the operator produces
     * \param coef1 the coefficient to multiply the first side field by
     * \param coef2 the coefficient to multiply the second side field by
     * \param coef3 the coefficient to multiply the third side field by
     * \param coef4 the coefficient to multiply the fourth side field by
     * \param src the field that the operator is applied to
     * \param dest the resulting field.
     */
    template<typename SrcType, typename DestType>
    inline void ref_stencil4_apply_to_field( const double coef1,
                                             const double coef2,
                                             const double coef3,
                                             const double coef4,
                                             const SrcType  & src,
                                                   DestType & dest );

    /**
     * \brief apply reference (FD 2-point) stencil on given fields with given coefficients
     * \tparam OpType - type of operator
     * \tparam FieldType - the type of field the operator is applied to and the type of field the operator produces
     * \param coefLo the coefficient to multiply the (-) side field by
     * \param coefHi the coefficient to multiply the (+) side field by
     * \param src the field that the operator is applied to
     * \param dest the resulting field.
     */
    template<typename OpType, typename FieldType>
    inline void ref_fd_stencil2_apply_to_field( const double coefLo,
                                                const double coefHi,
                                                const FieldType & src,
                                                      FieldType & dest );

    /*******************************************************************
     *
     * NOTE: all information in the RefStencil2Detail namespace is meant only for
     *       use within the ref_stencil2_apply_to_field function and should not be used elsewhere!
     *
     ******************************************************************/
    namespace RefStencil2Detail {

        /**
         *  \struct ActiveDir
         *
         *  \brief determines the direction that the operator is acting along
         *
         *  \tparam SrcT the source field type for this 2-point stencil operator
         *  \tparam DestT the destination field type for this 2-point stencil operator
         *
         *  Example:
         *  \code typedef ActiveDir<SrcT,DestT>::type OpDir; \endcode
         */
        template<typename SrcT, typename DestT>
        struct ActiveDir
        {
        private:
            typedef typename Subtract< typename  SrcT::Location::Offset,
                                       typename DestT::Location::Offset >::result  Difference;
        public:
            typedef typename GetNonzeroDir< Difference >::DirT type; ///< the direction that the operator acts in.
        };

        template<typename T> struct ActiveDir<T,T>;  ///< invalid - all ReferenceStencil2 must do something.

        /**
         * \struct ExtentsAndOffsets
         * \author James C. Sutherland
         * \brief Provides typedefs for dealing with extents and offsets for ReferenceStencil2 operators.
         */
        template<typename SrcT, typename DestT>
        struct ExtentsAndOffsets
        {
        private:
            typedef typename  SrcT::Location::Offset                SFO;            ///< Offset information for Source field
            typedef typename  SrcT::Location::BCExtra               SFBCExtra;      ///< Extra cell information for Source field
            typedef typename DestT::Location::Offset                DFO;            ///< Offset information for Destination field
            typedef typename DestT::Location::BCExtra               DFBCExtra;      ///< Extra cell information for Destination field

        public:
            typedef typename ActiveDir<SrcT,DestT>::type            Dir;            ///< The direction that this operator acts in
            typedef typename UnitTriplet<Dir>::type                 DirUnitVec;     ///< unit vector for the direction that this operator acts in.
            typedef IndexTriplet<0,0,0>                             Src1Offset;     ///< offset for the first source field
            typedef typename DirUnitVec::Negate                     Src1Extent;     ///< extent modification for the first source field
            typedef typename Subtract<SFBCExtra,
                DirUnitVec>::result::PositiveOrZero::Negate         Src1ExtentBC;   ///< amount to augment source extents by if a BC is present
            typedef typename Add<DirUnitVec,Src1Offset>::result     Src2Offset;     ///< offset for the second source field
            typedef Src1Extent                                      Src2Extent;     ///< extent modification for the second source field
            typedef Src1ExtentBC                                    Src2ExtentBC;   ///< additional extent modification if a BC is present
            typedef typename Multiply< DirUnitVec,
                typename Subtract<
                    typename Multiply<DirUnitVec,SFO>::result,
                    typename Multiply<DirUnitVec,DFO>::result
                    >::result
                >::result::PositiveOrZero                           DestOffset;     ///< the offset for the destination field
            typedef Src1Extent                                      DestExtent;     ///< the extent for the destination field
            typedef typename Subtract<
                typename Multiply<DFO,DFBCExtra>::result,
                    typename Multiply< DirUnitVec,
                        typename Multiply<SFO,SFBCExtra>::result
                    >::result
                >::result                                           DestExtentBC;   ///< amount to augment destination extents by if a BC is present
        };
    } // namespace RefStencil2Detail

    /*******************************************************************
     *
     * NOTE: all information in the RefStencil4Detail namespace is meant only for
     *       use within the ref_stencil4_apply_to_field function and should not be used elsewhere!
     *
     ******************************************************************/
    namespace RefStencil4Detail{
        typedef IndexTriplet<0,0,0>  ZeroTriplet;

        template< typename IT1, typename IT2, typename IT3> struct SelectNonzeroTriplet;

        template< typename IT1, typename IT2 > struct SelectNonzeroTriplet< IT1, IT2, ZeroTriplet >{ typedef IT1 type1;  typedef IT2 type2; };
        template< typename IT1, typename IT3 > struct SelectNonzeroTriplet< IT1, ZeroTriplet, IT3 >{ typedef IT1 type1;  typedef IT3 type2; };
        template< typename IT2, typename IT3 > struct SelectNonzeroTriplet< ZeroTriplet, IT2, IT3 >{ typedef IT2 type1;  typedef IT3 type2; };

        template< typename SrcT, typename DestT >
        struct ActiveDirs
        {
            typedef typename Subtract<
                typename  SrcT::Location::Offset,
                typename DestT::Location::Offset >::result                    Difference;
        private:
            typedef typename Multiply< typename UnitTriplet<XDIR>::type,
                                       Difference >::result::Abs              ActiveX;     // {1,0,0} if active
            typedef typename Multiply< typename UnitTriplet<YDIR>::type,
                                       Difference >::result::Abs              ActiveY;     // {0,1,0} if active
            typedef typename Multiply< typename UnitTriplet<ZDIR>::type,
                                       Difference >::result::Abs              ActiveZ;     // {0,0,1} if active

        public:
            typedef typename SelectNonzeroTriplet<ActiveX,ActiveY,ActiveZ>::type1  Dir1Vec;  ///< the first active direction UnitTriplet
            typedef typename SelectNonzeroTriplet<ActiveX,ActiveY,ActiveZ>::type2  Dir2Vec;  ///< the second active direction UnitTriplet
            typedef typename GetNonzeroDir< Dir1Vec >::DirT                        Dir1;     ///< the first active direction
            typedef typename GetNonzeroDir< Dir2Vec >::DirT                        Dir2;     ///< the second active direction
        };

        template< typename T > struct ActiveDirs<T,T>;  // invalid - all Stencil4 must do something.

        /**
         * \struct ExtentsAndOffsets
         * \author James C. Sutherland
         * \brief  Information about the extents and offsets for Stencil4
         */
        template< typename SrcT, typename DestT >
        struct ExtentsAndOffsets
        {
        private:
            typedef typename  SrcT::Location::Offset              SFO;            ///< Offset information for Source field
            typedef typename  SrcT::Location::BCExtra             SFBCExtra;      ///< Extra cell information for Source field
            typedef typename DestT::Location::Offset              DFO;            ///< Offset information for Destination field
            typedef typename DestT::Location::BCExtra             DFBCExtra;      ///< Extra cell information for Destination field

        public:
            typedef typename ActiveDirs<SrcT,DestT>::Dir1         Dir1;            ///< The first active direction
            typedef typename ActiveDirs<SrcT,DestT>::Dir2         Dir2;            ///< The second active direction
            typedef typename ActiveDirs<SrcT,DestT>::Dir1Vec      Dir1Vec;         ///< The first active direction
            typedef typename ActiveDirs<SrcT,DestT>::Dir2Vec      Dir2Vec;         ///< The second active direction
            typedef typename Multiply<SFBCExtra,Dir1Vec>::result  UpperLoopBCAug1; ///< shift for dir1 upper bounds when BC is present
            typedef typename Multiply<SFBCExtra,Dir2Vec>::result  UpperLoopBCAug2; ///< shift for dir2 upper bounds when BC is present
            typedef IndexTriplet<0,0,0>                           Src1Offset;      ///< offset for the first source field
            typedef typename Add<Dir1Vec,Dir2Vec>::result::Negate Src1Extent;      ///< extent modification for the first source field
            typedef typename Subtract<SFBCExtra,
                typename Add<Dir1Vec,Dir2Vec>::result
                >::result::PositiveOrZero::Negate                 Src1ExtentBC;    ///< amount to augment source1 extent by if a BC is present
            typedef typename Add<Dir1Vec,Src1Offset>::result      Src2Offset;      ///< offset for the second source field
            typedef Src1Extent                                    Src2Extent;      ///< extent modification for the second source field
            typedef Src1ExtentBC                                  Src2ExtentBC;    ///< amount to augment source2 extent by if a BC is present
            typedef typename Add<Dir2Vec,Src1Offset>::result      Src3Offset;      ///< offset for the third source field
            typedef Src1Extent                                    Src3Extent;      ///< extent modification for the third source field
            typedef Src1ExtentBC                                  Src3ExtentBC;    ///< amount to augment source3 extent by if a BC is present
            typedef typename Add<Dir1Vec,Src3Offset>::result      Src4Offset;      ///< offset for the fourth source field
            typedef Src1Extent                                    Src4Extent;      ///< extent modification for the fourth source field
            typedef Src3ExtentBC                                  Src4ExtentBC;    ///< amount to augment source4 extent by if a BC is present
            typedef typename Add<typename Multiply< Dir1Vec,
                    typename Subtract<
                        typename Multiply<Dir1Vec,SFO>::result,
                        typename Multiply<Dir1Vec,DFO>::result
                        >::result
                    >::result::PositiveOrZero,
                    typename Multiply< Dir2Vec,
                        typename Subtract<typename Multiply<Dir2Vec,SFO>::result,
                        typename Multiply<Dir2Vec,DFO>::result
                        >::result
                    >::result::PositiveOrZero
                >::result                                         DestOffset;      ///< the offset for the destination field
            typedef Src1Extent                                    DestExtent;      ///< the extent for the destination field
            typedef typename Subtract<typename Multiply<DFO,DFBCExtra>::result,
                typename Multiply< typename Add<Dir1Vec,Dir2Vec>::result,
                    typename Multiply<SFO,SFBCExtra>::result
                    >::result
                >::result                                         DestExtentBC;
        };
    } // namespace RefStencil4Detail

    //------------------------------------------------------------------

    template<typename SrcType, typename DestType>
    inline void ref_stencil2_apply_to_field( const double coefLo,
                                             const double coefHi,
                                             const SrcType  & src,
                                                   DestType & dest )
    {
        RefStencil2Detail::ExtentsAndOffsets<SrcType, DestType> typedef Extents;

        const MemoryWindow & wdest = dest.window_with_ghost();
        const MemoryWindow & wsrc = src.window_with_ghost();

        const MemoryWindow wd(wdest.glob_dim(),
                              wdest.offset() + Extents::DestOffset::int_vec(),
                              wdest.extent() + Extents::DestExtent::int_vec() + wdest.has_bc() * Extents::DestExtentBC::int_vec(),
                              wdest.has_bc(0),
                              wdest.has_bc(1),
                              wdest.has_bc(2));
        const MemoryWindow ws1(wsrc.glob_dim(),
                               wsrc.offset() + Extents::Src1Offset::int_vec(),
                               wsrc.extent() + Extents::Src1Extent::int_vec() + wsrc.has_bc() * Extents::Src1ExtentBC::int_vec(),
                               wsrc.has_bc(0),
                               wsrc.has_bc(1),
                               wsrc.has_bc(2));
        const MemoryWindow ws2(wsrc.glob_dim(),
                               wsrc.offset() + Extents::Src2Offset::int_vec(),
                               wsrc.extent() + Extents::Src2Extent::int_vec() + wsrc.has_bc() * Extents::Src2ExtentBC::int_vec(),
                               wsrc.has_bc(0),
                               wsrc.has_bc(1),
                               wsrc.has_bc(2));

#       ifndef NDEBUG
            assert(ws1.extent() == ws2.extent() && ws1.extent() == wd.extent());
#       endif

              DestType d(wd, dest.field_values(), ExternalStorage);
        const SrcType s1(ws1, src.field_values(), ExternalStorage);
        const SrcType s2(ws2, src.field_values(), ExternalStorage);

        typename DestType::iterator id = d.begin();
        typename DestType::iterator ide = d.end();
        typename SrcType::const_iterator is1 = s1.begin();
        typename SrcType::const_iterator is2 = s2.begin();

        for(; id != ide; ++id, ++is1, ++is2) {
            *id = *is1 * coefLo + *is2 * coefHi;
        };
    }

    //------------------------------------------------------------------

    template<typename SrcType, typename DestType>
    inline void ref_stencil4_apply_to_field( const double coef1,
                                             const double coef2,
                                             const double coef3,
                                             const double coef4,
                                             const SrcType  & src,
                                                   DestType & dest )
    {
        RefStencil4Detail::ExtentsAndOffsets<SrcType, DestType> typedef Extents;

        const MemoryWindow & wsrc = src.window_with_ghost();
        const MemoryWindow & wdest = dest.window_with_ghost();

        const MemoryWindow wd(wdest.glob_dim(),
                              wdest.offset() + Extents::DestOffset::int_vec(),
                              wdest.extent() + Extents::DestExtent::int_vec() + wdest.has_bc() * Extents::DestExtentBC::int_vec(),
                              wdest.has_bc(0),
                              wdest.has_bc(1),
                              wdest.has_bc(2));
        const MemoryWindow ws1(wsrc.glob_dim(),
                               wsrc.offset() + Extents::Src1Offset::int_vec(),
                               wsrc.extent() + Extents::Src1Extent::int_vec() + wsrc.has_bc() * Extents::Src1ExtentBC::int_vec(),
                               wsrc.has_bc(0),
                               wsrc.has_bc(1),
                               wsrc.has_bc(2));
        const MemoryWindow ws2(wsrc.glob_dim(),
                               wsrc.offset() + Extents::Src2Offset::int_vec(),
                               wsrc.extent() + Extents::Src2Extent::int_vec() + wsrc.has_bc() * Extents::Src2ExtentBC::int_vec(),
                               wsrc.has_bc(0),
                               wsrc.has_bc(1),
                               wsrc.has_bc(2));
        const MemoryWindow ws3(wsrc.glob_dim(),
                               wsrc.offset() + Extents::Src3Offset::int_vec(),
                               wsrc.extent() + Extents::Src3Extent::int_vec() + wsrc.has_bc() * Extents::Src3ExtentBC::int_vec(),
                               wsrc.has_bc(0),
                               wsrc.has_bc(1),
                               wsrc.has_bc(2));
        const MemoryWindow ws4(wsrc.glob_dim(),
                               wsrc.offset() + Extents::Src4Offset::int_vec(),
                               wsrc.extent() + Extents::Src4Extent::int_vec() + wsrc.has_bc() * Extents::Src4ExtentBC::int_vec(),
                               wsrc.has_bc(0),
                               wsrc.has_bc(1),
                               wsrc.has_bc(2));

#       ifndef NDEBUG
            assert(ws1.extent() == ws2.extent() && ws1.extent() == ws3.extent() && ws1.extent() == ws4.extent() && ws1.extent() == wd.extent());
#       endif

        DestType d(wd, dest.field_values(), ExternalStorage);
        SrcType s1(ws1, src.field_values(), ExternalStorage);
        SrcType s2(ws2, src.field_values(), ExternalStorage);
        SrcType s3(ws3, src.field_values(), ExternalStorage);
        SrcType s4(ws4, src.field_values(), ExternalStorage);

        typename DestType::iterator id = d.begin();
        typename DestType::iterator ide = d.end();
        typename SrcType::const_iterator is1 = s1.begin();
        typename SrcType::const_iterator is2 = s2.begin();
        typename SrcType::const_iterator is3 = s3.begin();
        typename SrcType::const_iterator is4 = s4.begin();

        for(; id != ide; ++id, ++is1, ++is2, ++is3, ++is4) {
            *id = *is1 * coef1 + *is2 * coef2 + *is3 * coef3 + *is4 * coef4;
        };
    };

    //------------------------------------------------------------------

    template<typename OpType, typename FieldType>
    inline void ref_fd_stencil2_apply_to_field( const double coefLo,
                                                const double coefHi,
                                                const FieldType & src,
                                                      FieldType & dest )
    {
        const MemoryWindow & w = src.window_with_ghost();
        typedef typename UnitTriplet<typename OpType::DirT>::type DirVec;
        const IntVec shift = DirVec::int_vec() + DirVec::int_vec();

        const MemoryWindow wd(w.glob_dim(),
                              w.offset() + DirVec::int_vec(),
                              w.extent() - shift,
                              w.has_bc(0),
                              w.has_bc(1),
                              w.has_bc(2));
        const MemoryWindow ws1(w.glob_dim(),
                               w.offset(),
                               w.extent() - shift,
                               w.has_bc(0),
                               w.has_bc(1),
                               w.has_bc(2));
        const MemoryWindow ws2(w.glob_dim(),
                               w.offset() + shift,
                               w.extent() - shift,
                               w.has_bc(0),
                               w.has_bc(1),
                               w.has_bc(2));

        FieldType d(wd, dest.field_values(), ExternalStorage);
        FieldType s1(ws1, src.field_values(), ExternalStorage);
        FieldType s2(ws2, src.field_values(), ExternalStorage);

        typename FieldType::iterator id = d.begin();
        typename FieldType::iterator ide = d.end();
        typename FieldType::const_iterator is1 = s1.begin();
        typename FieldType::const_iterator is2 = s2.begin();

        for(; id != ide; ++id, ++is1, ++is2) {
            *id = *is1 * coefLo + *is2 * coefHi;
        };
    };

  }// namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Reference_Stencil_h
