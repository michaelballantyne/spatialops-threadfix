#ifndef ParticleOperators_h
#define ParticleOperators_h

#include <spatialops/particles/ParticleFieldTypes.h>

#include <stdexcept>
#include <math.h>

namespace SpatialOps{
namespace Particle{

   /**
    *  Interpolates a particle field onto an underlying mesh field.      
    */
   template< typename CellField >
   class ParticleToCell
   {
   public:
     typedef ParticleField SrcFieldType;
     typedef CellField     DestFieldType;

     /**
      *  @param meshCoord Vector of coordinates for the underlying mesh
      */
     ParticleToCell( const CellField& meshCoord );

     /**
    * @param particleCoord Field of coordinates for all particles (ParticleField)
    * @param src source field from which values are interpolated to partciles (ParticleField)
    * @param dest destination field to which values are interpolated (CellField)
    */
     void apply_to_field( const ParticleField& particleCoord,
                          const ParticleField& particleSize,
                          const SrcFieldType& src,
                          DestFieldType& dest ) const;

   private:
     const CellField& coordVec_;
     const double dx_, xo_;
   };


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
     *  @param meshCoord Field of coordinates for the underlying mesh
     */   
    CellToParticle( const CellField& meshCoord ); 

    /**
    * @param particleCoord Field of coordinates for all particles (ParticleField)
    * @param src source field from which values are interpolated to partciles (VolField)
    * @param dest destination field to which values are interpolated (ParticleField)
    */
    void apply_to_field( const ParticleField& particleCoord,
                         const SrcFieldType& src,
                         DestFieldType& dest ) const;
   

  private:
    const CellField& coordVec_;
    const double dx_, xo_;
  };



  // =================================================================
  //
  //                           Implementation
  //
  // =================================================================



   template< typename CellField >
   ParticleToCell<CellField>::
   ParticleToCell( const CellField& meshCoord )
     : coordVec_( meshCoord ),
       dx_( meshCoord[1]-meshCoord[0] ),
       xo_( meshCoord[0] )
   {
     const double TOL = 1e-6;

     // note that this assumes 1D
     bool isUniform = true;
     typename CellField::const_iterator ix2=meshCoord.begin();
     typename CellField::const_iterator ix = ix2;
     ++ix2;
     for( ; ix2!=meshCoord.end(); ++ix, ++ix2 ){
       if( fabs( dx_ - (*ix2-*ix) )/dx_ > TOL ){
         isUniform = false;
       }
     }
     if( !isUniform )
       throw std::runtime_error( "Particle operators require uniform mesh spacing" );
   }

   //------------------------------------------------------------------

   template< typename CellField >
   void
   ParticleToCell<CellField>::
   apply_to_field( const ParticleField& particleCoord,
                   const ParticleField& particleSize,
                   const SrcFieldType& src,
                   DestFieldType& dest ) const
   {
     dest = 0.0;
     ParticleField::const_iterator plociter = particleCoord.begin();
     //ParticleField::const_iterator psizeiter = particleSize.begin();
     ParticleField::const_iterator isrc = src.begin();

     const double halfwidth = 0.5*dx_;

     for( ; plociter!=particleCoord.end(); ++plociter, ++isrc ){
       // given the current particle coordinate,
       // determine what cell it is located in.
       const size_t cellIx1 = size_t((*plociter-halfwidth-xo_) / dx_);
       const size_t cellIx2 = cellIx1 +1 ;
       const double leftloc = coordVec_[cellIx1];
       const double rightloc = coordVec_[cellIx2];
       //std::cout<<" cellIx1 : "<<cellIx1<<"  cellIx1 : "<<cellIx2<<"  leftloc  : "<< leftloc<<"  rightloc : "<<rightloc<<std::endl;
       if( fabs( *plociter - leftloc) <= fabs( *plociter - rightloc ))
         dest[cellIx1] += *isrc;
       else
         dest[cellIx2] += *isrc;
     }
   }

 //==================================================================

  template< typename CellField >
  CellToParticle<CellField>::
  CellToParticle( const CellField& meshCoord )
    : coordVec_( meshCoord ),
      dx_( meshCoord[1]-meshCoord[0] ),
      xo_( meshCoord[0] )
  {
    const double TOL = 1e-6;

    // note that this assumes 1D
    bool isUniform = true;
    typename CellField::const_iterator ix2=meshCoord.begin();
    typename CellField::const_iterator ix = ix2;
    ++ix2;
    for( ; ix2!=meshCoord.end(); ++ix, ++ix2 ){     
      if( fabs( dx_ - (*ix2-*ix) )/dx_ > TOL ){
        isUniform = false;
      }
    }
    if( !isUniform )
      throw std::runtime_error( "Particle operators require uniform mesh spacing" );
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
      const size_t i1 = size_t( (xp-halfwidth-xo_) / dx_ );     
      const double x1 = coordVec_[i1];     
      const size_t i2 = ( xp > x1 ) ? i1+1 : i1-1;     
      const double x2 = coordVec_[i2];      
      const double c1 = (xp-x2)/(x1-x2);
      const double c2 = (xp-x1)/(x2-x1);
      *idest = src[i1]*c1 + src[i2]*c2;

    }
  }

} // namespace Particle
} // namespace SpatialOps

#endif // ParticleOperators_h
