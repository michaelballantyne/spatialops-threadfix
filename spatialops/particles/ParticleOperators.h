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
#ifndef ParticleOperators_h
#define ParticleOperators_h

#include <spatialops/particles/ParticleFieldTypes.h>
#include <spatialops/Nebo.h>

#include <stdexcept>
#include <cmath>

namespace SpatialOps{
namespace Particle{

   /**
    *  @class  ParticleToCell
    *  @author James C. Sutherland
    *  @brief Interpolates a particle field onto an underlying mesh field.
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
    *  @param particleCoord Field of coordinates for all particles (ParticleField)
    *  @param particleSize Field of size for all particles (ParticleField)
    *  @param src source field from which values are interpolated to partciles (ParticleField)
    *  @param dest destination field to which values are interpolated (CellField)
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
    * @param particleSize Field of size for all particles (ParticleField)
    * @param src source field from which values are interpolated to partciles (VolField)
    * @param dest destination field to which values are interpolated (ParticleField)
    */
    void apply_to_field( const ParticleField& particleCoord,
                         const ParticleField& particleSize,
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
     dest <<= 0.0;

     //const double halfwidth = 0.5*dx_;

#   ifndef NDEBUG
    const int nmax = dest.window_with_ghost().local_npts();
#   endif

    ParticleField::const_iterator plociter = particleCoord.begin();
    ParticleField::const_iterator psizeiter = particleSize.begin();
    ParticleField::const_iterator isrc = src.begin();
    for( ; plociter!=particleCoord.end(); ++plociter, ++isrc, ++psizeiter ){
      //particle location is the position where particle center is located
      //identify the location of the particle left boundary
      double leftloc = *plociter-( *psizeiter / 2) ;
      //identify the location of the particle right boundary
      const double rightloc = *plociter+( 0.5 * *psizeiter );
      //identify the cell index where particle left boundary located
      int leftcellIx1 = int((*plociter-( 0.5 * *psizeiter )) / dx_) + 1;
#     ifndef NDEBUG
      if( (leftcellIx1 >= nmax) || (leftcellIx1<0) ){
        throw std::runtime_error( "Particle is outside of the domain!" );
      }
#     endif
      //std::cout<<"leftloc  :  "<<leftloc<<"  rightloc  :  "<<rightloc<<"   leftcellIx1 : "<<leftcellIx1<<std::endl;
      //loop through all cells effected by the particle and add the source term
      while(leftloc < rightloc){
        double rb = coordVec_[leftcellIx1] + dx_/2 ;
        if(rb > rightloc)
          rb = rightloc ;
        dest[leftcellIx1] += *isrc * (rb-leftloc) / *psizeiter ;
        leftcellIx1 += 1 ;
        leftloc = rb ;
      }
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
                  const ParticleField& particleSize,
                  const SrcFieldType& src,
                  DestFieldType& dest ) const
  {
    dest <<= 0.0;

    //const double halfwidth = 0.5*dx_;
#   ifndef NDEBUG
    const int nmax = src.window_with_ghost().local_npts();
#   endif
    ParticleField::const_iterator plociter = particleCoord.begin();
    ParticleField::const_iterator psizeiter = particleSize.begin();
    ParticleField::iterator destiter = dest.begin();
    for( ; plociter!=particleCoord.end(); ++plociter, ++destiter, ++psizeiter ){
      double leftloc = *plociter-( *psizeiter / 2.0) ;
      const double rightloc = *plociter+( *psizeiter / 2.0) ;
      int leftcellIx1 = int((*plociter-( *psizeiter / 2.0)) / dx_) + 1;
#     ifndef NDEBUG
      if( leftcellIx1 >= nmax || leftcellIx1<0 ){
        throw std::runtime_error( "Particle is outside of the domain!" );
      }
#     endif
      //std::cout<<"leftloc  :  "<<leftloc<<"  rightloc  :  "<<rightloc<<"   leftcellIx1 : "<<leftcellIx1<<std::endl;
      while(leftloc < rightloc){
        double rb = coordVec_[leftcellIx1] + dx_/2.0 ;
        if(rb > rightloc)
          rb = rightloc ;
        *destiter += src[leftcellIx1] * (rb-leftloc) / *psizeiter ;
        leftcellIx1 += 1 ;
        leftloc = rb ;
      }
    }
  }

} // namespace Particle
} // namespace SpatialOps

#endif // ParticleOperators_h
