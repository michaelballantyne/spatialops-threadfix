#ifndef SpatialOps_Structured_Filter_h
#define SpatialOps_Structured_Filter_h

#include <spatialops/structured/IndexTriplet.h>

namespace SpatialOps{
namespace structured{

  /**
   *  \class Filter
   *  \date October, 2011
   *  \author James C. Sutherland
   *  \brief Provides a box filter operator for uniform structured meshes.
   */
  template< typename FieldT >
  class BoxFilter{
    typedef typename FieldT::const_iterator ConstFieldIter;

    mutable std::vector<FieldT> srcFields_;
    mutable std::vector<ConstFieldIter> srcIters_;

  public:
    BoxFilter(){};
    ~BoxFilter(){};

    void apply_to_field( const FieldT& src, FieldT& dest ) const;
  };

} // namespace structured
} // namespace SpatialOps

#endif /* SpatialOps_Structured_Filter_h */
