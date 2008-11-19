#ifndef Test_Functions_h
#define Test_Functions_h

#include <spatialops/FieldFunctions.h>

//====================================================================


template<typename FieldT>
class LinearFunction : public SpatialOps::FieldFunction3D<FieldT>
{
public:
  LinearFunction( const FieldT& x, const FieldT& y, const FieldT& z )
    : SpatialOps::FieldFunction3D<FieldT>(x,y,z)
  {}
  ~LinearFunction(){}
  void evaluate( FieldT& phi ) const
  {
    typename FieldT::const_iterator ix=this->get_x().begin();
    typename FieldT::const_iterator iy=this->get_y().begin();
    typename FieldT::const_iterator iz=this->get_z().begin();
    for( typename FieldT::iterator iphi=phi.begin(); iphi!=phi.end(); ++iphi,++ix,++iy,++iz ){
      *iphi = 2*(*ix)+3*(*iy)+4*(*iz);
    }
  }

  void dx( FieldT& gradPhi ) const{ gradPhi = 2.0; }
  void dy( FieldT& gradPhi ) const{ gradPhi = 3.0; }
  void dz( FieldT& gradPhi ) const{ gradPhi = 4.0; }

  void d2x( FieldT& d2phi ) const{ d2phi=0.0; }
  void d2y( FieldT& d2phi ) const{ d2phi=0.0; }
  void d2z( FieldT& d2phi ) const{ d2phi=0.0; }

private:
};


//====================================================================


template<typename FieldT>
class SinFun : public SpatialOps::FieldFunction3D<FieldT>
{
public:

  SinFun( const FieldT& x, const FieldT& y, const FieldT& z )
    : SpatialOps::FieldFunction3D<FieldT>(x,y,z),
      pi( 3.141592653589793 )
  {}
  ~SinFun(){}
  void evaluate( FieldT& phi ) const
  {
    typename FieldT::const_iterator ix=this->get_x().begin();
    typename FieldT::const_iterator iy=this->get_y().begin();
    typename FieldT::const_iterator iz=this->get_z().begin();
    for( typename FieldT::iterator iphi=phi.begin(); iphi!=phi.end(); ++iphi,++ix,++iy,++iz ){
      *iphi = std::sin(*ix*pi) + std::sin(0.5**iy*pi) + std::cos(*iz*pi) + 1;
    }
  }

  void dx( FieldT& gradPhi ) const
  {
    typename FieldT::const_iterator ix=this->get_x().begin();
    for( typename FieldT::iterator igrad=gradPhi.begin(); igrad!=gradPhi.end(); ++igrad,++ix ){
      *igrad = pi*std::cos(*ix*pi);
    }
  }

  void dy( FieldT& grad ) const
  {
    typename FieldT::const_iterator iy=this->get_y().begin();
    for( typename FieldT::iterator igrad=grad.begin(); igrad!=grad.end(); ++igrad,++iy ){
      *igrad = 0.5*pi*std::cos(0.5**iy*pi);
    }
  }

  void dz( FieldT& grad ) const
  {
    typename FieldT::const_iterator iz=this->get_z().begin();
    for( typename FieldT::iterator igrad=grad.begin(); igrad!=grad.end(); ++igrad,++iz ){
      *igrad =  -pi*std::sin(*iz*pi);
    }
  }


  void d2x( FieldT& d2phi ) const
  {
    typename FieldT::const_iterator ix=this->get_x().begin();
    for( typename FieldT::iterator iphi=d2phi.begin(); iphi!=d2phi.end(); ++iphi,++ix ){
      *iphi = -pi*pi*std::sin(*ix*pi);
    }
  }


  void d2y( FieldT& d2phi ) const
  {
    typename FieldT::const_iterator iy=this->get_y().begin();
    for( typename FieldT::iterator iphi=d2phi.begin(); iphi!=d2phi.end(); ++iphi,++iy ){
      *iphi = -0.25*pi*pi*std::sin(0.5**iy*pi);
    }
  }

  void d2z( FieldT& d2phi ) const
  {
    typename FieldT::const_iterator iz=this->get_z().begin();
    for( typename FieldT::iterator iphi=d2phi.begin(); iphi!=d2phi.end(); ++iphi,++iz ){
      *iphi = -pi*pi*std::cos(*iz*pi);
    }
  }

private:
  const double pi;
};



#endif
