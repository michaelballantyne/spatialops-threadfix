#include <vector>
#include <iostream>

#include <spatialops/OperatorDatabase.h>
#include <spatialops/structured/FVStaggered.h>

#include "Functions.h"

using namespace std;
using namespace SpatialOps;
using namespace structured;

//====================================================================

bool is_result_valid( const SVolField& f1,
                      const SVolField& f2 )
{
  bool isvalid = true;
  for( SVolField::const_iterator if1=f1.begin(), if2=f2.begin(); if1!=f1.end(); ++if1, ++if2 ){
    isvalid = (*if1 == *if2);
    if( !isvalid ){
      cout << "invalid: " << *if1 << ", " << *if2 << endl;
      break;
    }
  }
  return isvalid;
}

//====================================================================

void build_mesh( const IntVec& dim, SVolField& x, SVolField& y, SVolField& z )
{
  vector<double> spacing(3,1.0);
  for( int i=0; i<3; ++i ){
    spacing[i] = 1.0 / dim[i];
  }

  const int ihi = get_nx_with_ghost<SVolField>(dim[0],true);
  const int jhi = get_ny_with_ghost<SVolField>(dim[1],true);
  const int khi = get_nz_with_ghost<SVolField>(dim[2],true);

  const int ngx = dim[0]>1 ? SVolField::Ghost::NGHOST : 0;
  const int ngy = dim[1]>1 ? SVolField::Ghost::NGHOST : 0;
  const int ngz = dim[2]>1 ? SVolField::Ghost::NGHOST : 0;

  for( int k=0; k<khi; ++k ){
    const double zval = spacing[2]*(double(k)+0.5-ngz);
    for( int j=0; j<jhi; ++j ){
      const double yval = spacing[1]*(double(j)+0.5-ngy);
      for( int i=0; i<ihi; ++i ){
        const double xval = spacing[0]*(double(i)+0.5-ngx);
        x(i,j,k) = xval;
        y(i,j,k) = yval;
        z(i,j,k) = zval;
      }
    }
  }
}

//====================================================================

bool restrict_test( const int activeDim,
                    const IntVec& dimSrc, 
                    const IntVec& dimDest )
{
  OperatorDatabase opDB;

  RestrictionAssembler<SVolField> ra( dimSrc, dimDest, true, true, true );
  opDB.register_new_operator< RestrictSVol >( new RestrictSVol(ra) );

  SVolField coarse( get_dim_with_ghost<SVolField>(dimSrc, true,true,true), NULL );
  SVolField fine  ( get_dim_with_ghost<SVolField>(dimDest,true,true,true), NULL );
  SVolField fine2 ( get_dim_with_ghost<SVolField>(dimDest,true,true,true), NULL );

  SVolField xcoarse( get_dim_with_ghost<SVolField>(dimSrc,true,true,true), NULL );
  SVolField ycoarse( get_dim_with_ghost<SVolField>(dimSrc,true,true,true), NULL );
  SVolField zcoarse( get_dim_with_ghost<SVolField>(dimSrc,true,true,true), NULL );

  SVolField xfine( get_dim_with_ghost<SVolField>(dimDest,true,true,true), NULL );
  SVolField yfine( get_dim_with_ghost<SVolField>(dimDest,true,true,true), NULL );
  SVolField zfine( get_dim_with_ghost<SVolField>(dimDest,true,true,true), NULL );

  build_mesh( dimSrc, xcoarse, ycoarse, zcoarse );

  SinFun<SVolField> func( xcoarse, ycoarse, zcoarse );
  //LinearFunction<SVolField> func( xcoarse, ycoarse, zcoarse );
  func.evaluate( coarse );

  const RestrictSVol* const restrict = opDB.retrieve_operator<RestrictSVol>();
  restrict->apply_to_field( coarse, fine );

  // now restrict the coordinate fields.  Then see if the function
  // evaluated at the restricted coordinates is equal to the
  // restricted field.  It should be.
  restrict->apply_to_field( xcoarse, xfine );
  restrict->apply_to_field( ycoarse, yfine );
  restrict->apply_to_field( zcoarse, zfine );

  SinFun<SVolField> f2( xfine, yfine, zfine );
  f2.evaluate( fine2 );

  //#define WRITE_OUTPUT
# ifdef WRITE_OUTPUT
  coarse.write_matlab("coarse",true);
  fine.write_matlab("fine",true);
  xcoarse.write_matlab("xc",true);
  ycoarse.write_matlab("yc",true);
  zcoarse.write_matlab("zc",true);
  xfine.write_matlab("xf",true);
  yfine.write_matlab("yf",true);
  zfine.write_matlab("zf",true);
  restrict->write_matlab("R");
  fine2.write_matlab("f2",true);
# endif

  return is_result_valid( fine, fine2 );
}

//====================================================================

int main()
{
  IntVec dimSrc(1,1,1), dimDest(1,1,1);

  dimSrc [0] = 3;    dimSrc [1] = 10;    dimSrc [2] = 5;
  dimDest[0] = 16;   dimDest[1] = 10;    dimDest[2] = 5;
  const bool xpass = restrict_test( 0, dimSrc, dimDest );
  cout << "x: ";
  if( xpass )  cout << "PASS" << endl;
  else         cout << "FAIL" << endl;

  dimSrc [0] = 14;   dimSrc [1] = 10;    dimSrc [2] = 6;
  dimDest[0] = 14;   dimDest[1] = 21;    dimDest[2] = 6;
  const bool ypass = restrict_test( 1, dimSrc, dimDest );
  cout << "y: ";
  if( ypass )  cout << "PASS" << endl;
  else         cout << "FAIL" << endl;

  dimSrc [0] = 14;   dimSrc [1] = 11;    dimSrc [2] = 23;
  dimDest[0] = 14;   dimDest[1] = 11;    dimDest[2] = 61;
  const bool zpass = restrict_test( 2, dimSrc, dimDest );
  cout << "z: ";
  if( zpass )  cout << "PASS" << endl;
  else         cout << "FAIL" << endl;

  if( xpass && ypass && zpass ) return 0;
  else return -1;
}
