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
    };

  }
}


#endif /* FDSTENCIL2_H_ */
