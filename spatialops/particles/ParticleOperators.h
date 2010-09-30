#ifndef ParticleOperators_h
#define ParticleOperators_h

#include <spatialops/particles/ParticleFieldTypes.h>

namespace SpatialOps{
namespace Particle{

//   /**
//    *  Interpolates a particle field onto an underlying mesh field.
//    *
//    *  JCS: it isn't clear how to do this appropriately since we may
//    *  have multiple particles in a single cell...
//    */
//   template< typename CellField >
//   class ParticleToCell
//   {
//   public:
//     typedef ParticleField SrcFieldType;
//     typedef CellField     DestFieldType;

//     /**
//      *  @param particleCoord Field of coordinates for all particles
//      *  @param meshCoord Field of coordinates for the underlying mesh
//      */
//     ParticleToCell( const ParticleField& particleCoord,
//                     const CellField& meshCoord );

//     void apply_to_field( const SrcFieldType&,
//                          DestFieldType& ) const;

//   private:
//     const ParticleField& particleCoord_;
//     const CellField& cellCoord_;
//     const double dx_, xo_;
//   };


  //==================================================================


  /**
   *  @class CellToParticle
   *  @brief Operator to interpolate a mesh field onto a particle.
   *  @author James C. Sutherland
   */
  template< typename CellField >
  class CellToParticle
  {
  public:
    typedef CellField     SrcFieldType;
    typedef ParticleField DestFieldType;

    /**
     *  @param particleCoord Field of coordinates for all particles
     *  @param meshCoord Field of coordinates for the underlying mesh
     */
    CellToParticle( const CellField& meshCoord );

    void apply_to_field( const ParticleField& particleCoord,
                         const SrcFieldType&,
                         DestFieldType& ) const;
  private:
    const CellField& cCoord_;
    const double dx_, xo_;
  };



  // =================================================================
  //
  //                           Implementation
  //
  // =================================================================



  // template< typename CellField >
//   ParticleToCell<CellField>::
//   ParticleToCell( const ParticleField& particleCoord,
//                   const CellField& meshCoord )
//     : particleCoord_( particleCoord ),
//       cellCoord_    ( meshCoord     ),
//       dx_( meshCoord[1]-meshCoord[0] ),
//       xo_( meshCoord[0] )
//   {
//     const double TOL = 1e-6;

//     // note that this assumes 1D
//     bool isUniform = true;
//     typename CellField::const_iterator ix2=meshCoord.begin();
//     typename CellField::const_iterator ix = ix2++;
//     for( ; ix2!=meshCoord.end(); ++ix, ++ix2 ){
//       if( std::abs( dx_ - (*ix2-*ix) )/dx_ > TOL ){
//         isUniform = false;
//       }
//     }
//     assert( isUniform );
//   }

//   //------------------------------------------------------------------

//   template< typename CellField >
//   void
//   ParticleToCell<CellField>::
//   apply_to_field( const SrcFieldType& src,
//                   DestFieldType& dest ) const
//   {
//     dest = 0.0;
//     ParticleField::const_iterator ipx=particleCoord_.begin();
//     ParticleField::const_iterator isrc = src.begin();

//     const double halfwidth = 0.5*dx_;

//     for( ; ipx!=particleCoord_.end(); ++ipx, ++isrc ){
//       // given the current particle coordinate,
//       // determine what cell it is located in.
//       const size_t cellIx = (*ipx-halfwidth-xo_) / dx_ -1;
//       dest[cellIx] = *isrc;
//     }
//   }

  //==================================================================

  template< typename CellField >
  CellToParticle<CellField>::
  CellToParticle( const CellField& meshCoord )
    : cCoord_( meshCoord ),
      dx_( meshCoord[1]-meshCoord[0] ),
      xo_( meshCoord[0] )
  {
    const double TOL = 1e-6;

    // note that this assumes 1D
    bool isUniform = true;
    typename CellField::const_iterator ix2=meshCoord.begin();
    typename CellField::const_iterator ix = ix2; ++ix2;
    for( ; ix2!=meshCoord.end(); ++ix, ++ix2 ){
      if( std::abs( dx_ - (*ix2-*ix) )/dx_ > TOL ){
        isUniform = false;
      }
    }
    assert( isUniform );
  }

  //------------------------------------------------------------------

  template<typename CellField>
  void
  CellToParticle<CellField>::
  apply_to_field( const ParticleField& particleCoord,
                  const SrcFieldType& src,
                  DestFieldType& dest ) const
  {
    dest = 0.0;
    ParticleField::const_iterator ipx= particleCoord.begin();
    ParticleField::iterator idest = dest.begin();

    const double halfwidth = 0.5*dx_;

    for( ; ipx!=particleCoord.end(); ++ipx, ++idest ){
      // given the current particle coordinate, determine what cell it
      // is located in.  Then interpolate the src values (from the
      // mesh) onto the particle using linear interpolation
      const double xp = *ipx;
      const size_t i1 = size_t( (xp-halfwidth-xo_) / dx_ -1 );
      const double x1 = cCoord_[i1];
      const size_t i2 = ( xp > x1 ) ? i1+1 : i1-1;
      const double x2 = cCoord_[i2];
      const double c1 = (xp-x2)/(x1-x2);
      const double c2 = (xp-x1)/(x2-x1);
      *idest = src[i1]*c1 + src[i2]*c2;
    }
  }

  //------------------------------------------------------------------

} // namespace Particle
} // namespace SpatialOps

#endif // ParticleOperators_h
