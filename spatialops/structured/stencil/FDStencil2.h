#ifndef FDSTENCIL2_H_
#define FDSTENCIL2_H_

#include <spatialops/structured/IndexTriplet.h>

namespace SpatialOps{
  namespace structured{

    template< typename OpT, typename FieldT, typename DirT >
    class FDStencil2{
      const double coefLo_, coefHi_;

      typedef typename UnitTriplet<DirT>::type  DirVec;

    public:
      typedef OpT       Type;
      typedef FieldT    SrcFieldType;
      typedef FieldT    DestFieldType;

      FDStencil2( const double coefLo, const double coefHi );
      ~FDStencil2();

      void apply_to_field( const SrcFieldType& src, DestFieldType& dest ) const;

      inline double get_minus_coef() const{ return coefLo_; } ///< get the (-) coefficient
      inline double  get_plus_coef() const{ return coefHi_; } ///< get the (+) coefficient
    };

  }
}


#endif /* FDSTENCIL2_H_ */
