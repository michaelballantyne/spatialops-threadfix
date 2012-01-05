#ifndef SpatialOps_structured_Grid_h
#define SpatialOps_structured_Grid_h

#include <vector>
#include <spatialops/structured/MemoryWindow.h>

namespace SpatialOps{
namespace structured{

  /**
   *  \class Grid
   *  \author James C. Sutherland
   *  \brief Provides a simple interface to set coordinate fields.
   */
  class Grid
  {
    const IntVec npts_;
    const std::vector<double> length_;
    std::vector<double> spacing_;
  public:
    /**
     *  \param npts the number of scalar volume cells in each direction
     *  \param length the domain length in each direction
     */
    Grid( const IntVec npts,
          const std::vector<double>& length );

    /**
     *  \brief obtain the grid spacing in the requested direction
     */
    template< typename CoordT >
    double spacing() const;

    inline const IntVec& extent() const{ return npts_; }

    inline int extent( const size_t i ) const{ return npts_[i]; }

    inline const std::vector<double>& length() const{ return length_; }
    inline double length( const size_t i ) const{ return length_[i]; }

    /**
     *  \brief set the coordinates on the given field.
     *
     *  Examples:
     *  \code
     *  SVolField xsvol, ysvol;
     *  grid.set_coord<XDIR>( xsvol );
     *  grid.set_coord<YDIR>( ysvol );
     *
     *  XSurfYField xxsy, zxsy;
     *  grid.set_coord<XDIR>( xxsy );
     *  grid.set_coord<ZDIR>( zxsy );
     *  \endcode
     */
    template< typename CoordT, typename FieldT >
    void set_coord( FieldT& f ) const;
  };

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_structured_Grid_h
