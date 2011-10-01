#ifndef SpatialOps_Structured_Stencil_h
#define SpatialOps_Structured_Stencil_h

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>

namespace SpatialOps {
  namespace structured {

    /**
     *  \class Stencil2
     *  \author James C. Sutherland
     *
     *  \brief Support for implementing simple two-point stencils in
     *         one-dimension on structured meshes.
     *
     *  \tparam OpT - the type of operator
     *  \tparam SrcT - the type of field the operator is applied to
     *  \tparam DestT - the type of field the operator produces
     *
     *  See also Stencil2Helper
     */
    template<typename OperatorT, typename SrcFieldT, typename DestFieldT>
    class Stencil2
    {
      const double coefLo_, coefHi_;
    public:

      typedef OperatorT OpT;
      typedef SrcFieldT SrcFieldType;
      typedef DestFieldT DestFieldType;

      /**
       *  \brief construct a stencil with the specified coefficients
       *  \param coefLo the coefficient to multiply the (-) side field by
       *  \param coefHi the coefficient to multiply the (+) side field by
       */
      Stencil2( const double coefLo, const double coefHi );

      ~Stencil2();

      void apply_to_field( const SrcFieldType& src, DestFieldType& dest ) const;

      double get_minus_coef() const{ return coefLo_; }
      double  get_plus_coef() const{ return coefHi_; }
    };

    /*******************************************************************
     *
     * NOTE: all information in the s2detail namespace is meant only for
     *       use within the Stencil2Helper class and should not be used
     *       elsewhere!
     *
     ******************************************************************/
    namespace s2detail {
//
//      /**
//       *  \fn unsigned int stride( const MemoryWindow& mw )
//       *
//       *  \brief Obtain the stride (flat index) in the requested direction
//       *         for the given field.  This is a "local" stride, not the global one.
//       *
//       *  \tparam DirT the direction of interest
//       *  \param mw the MemoryWindow associated with the field
//       */
//      template< typename DirT > unsigned int stride( const MemoryWindow& mw );
//
//      template<> inline unsigned int stride<XDIR>( const MemoryWindow& mw ){ return 1; }
//      template<> inline unsigned int stride<YDIR>( const MemoryWindow& mw ){ return mw.extent(0); }
//      template<> inline unsigned int stride<ZDIR>( const MemoryWindow& mw ){ return mw.extent(0)*mw.extent(1); }

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

      template<typename T> struct ActiveDir<T,T>;  ///< invalid - all Stencil2 must do something.


      /**
       * \struct ExtentsAndOffsets
       * \brief Provides typedefs for dealing with extents and offsets for Stencil2 operators.
       */
      template<typename SrcT, typename DestT>
      struct ExtentsAndOffsets
      {
      private:
        typedef typename  SrcT::Location::Offset                SFO;            ///< Offset information for Source field
        typedef typename  SrcT::Location::BCExtra               SFBCExtra;      ///< Extra cell information for Source field
        typedef typename DestT::Location::Offset                DFO;            ///< Offset information for Destination field
        typedef typename DestT::Location::BCExtra               DFBCExtra;      ///< Extra cell information for Destination field
        typedef typename ActiveDir<SrcT,DestT>::type            Dir;            ///< the direction that this operator acts in

        typedef typename UnitTriplet<Dir>::type                 DirUnitVec;     ///< unit vector for the direction that this operator acts in.

      public:

        // jcs these two should be private, after we fix the rtests.
        typedef typename DirUnitVec::Negate                     UpperLoopShift; ///< upper bounds on loop

        typedef typename Multiply<SFBCExtra,DirUnitVec>::result UpperLoopBCAug; ///< shift for upper bounds when BC is present

        typedef IndexTriplet<0,0,0>                             Src1Offset;     ///< offset for the first source field
        typedef UpperLoopShift                                  Src1Extent;     ///< extent modification for the first source field
        typedef typename Subtract<SFBCExtra,
            DFBCExtra>::result::PositiveOrZero::Negate          Src1ExtentBC;   ///< amount to augment source extents by if a BC is present

        typedef typename Add<DirUnitVec,Src1Offset>::result     Src2Offset;     ///< offset for the second source field
        typedef Src1Extent                                      Src2Extent;     ///< extent modification for the second source field
        typedef Src1ExtentBC                                    Src2ExtentBC;   ///< additional extent modification if a BC is present

        typedef typename Multiply< DirUnitVec,
            typename Subtract<SFO,DFO>::result
            >::result::PositiveOrZero                           DestOffset;     ///< the offset for the destination field
        typedef UpperLoopShift                                  DestExtent;     ///< the extent for the destination field
        typedef typename Subtract<DFBCExtra,
            SFBCExtra>::result::PositiveOrZero::Negate          DestExtentBC;   ///< amount to augment destination extents by if a BC is present

        //  \todo rip these out and modify corresponding regression test types...
        typedef IndexTriplet< 1,
            Kronecker<XDIR::value,Dir::value>::value * Abs<SFBCExtra::X-DFBCExtra::X>::result,
            Kronecker<YDIR::value,Dir::value>::value * Abs<SFBCExtra::Y-DFBCExtra::Y>::result
            >                                                   SrcIterInc;     ///< The increment for the source field iterators

        typedef SrcIterInc                                      DestIterInc;    ///< The increment for the destination field iterators

        typedef IndexTriplet< 0,
                              SFO::X * DFO::X,
                              SFO::Y * DFO::Y >                 SrcIterBCAug;

        typedef IndexTriplet< 0,
                              SFO::X - DFO::X,
                              SFO::Y - DFO::Y >                 DestIterBCAug;
      };

    } // namespace s2detail

  }// namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Structured_Stencil_h
