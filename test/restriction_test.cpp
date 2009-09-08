#include <vector>
#include <iostream>

#include <spatialops/FVStaggered.h>
#include <spatialops/OperatorDatabase.h>

#include <Functions.h>

using namespace std;
using namespace SpatialOps;
using namespace FVStaggered;

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

void build_mesh( const vector<int>& dim, SVolField& x, SVolField& y, SVolField& z )
{
  vector<double> spacing(3,1.0);
  for( int i=0; i<3; ++i ){
    spacing[i] = 1.0 / dim[i];
  }

  SVolField::iterator isvx=x.begin();
  SVolField::iterator isvy=y.begin();
  SVolField::iterator isvz=z.begin();

  const int ihi = get_nx<SVolField>(dim,true);
  const int jhi = get_ny<SVolField>(dim,true);
  const int khi = get_nz<SVolField>(dim,true);

  const int ngxm = dim[0]>1 ? SVolField::Ghost::NM : 0;
  const int ngym = dim[1]>1 ? SVolField::Ghost::NM : 0;
  const int ngzm = dim[2]>1 ? SVolField::Ghost::NM : 0;

  for( int k=0; k<khi; ++k ){
    const double z = spacing[2]*(double(k)+0.5-ngzm);
    for( int j=0; j<jhi; ++j ){
      const double y = spacing[1]*(double(j)+0.5-ngym);
      for( int i=0; i<ihi; ++i ){
        const double x = spacing[0]*(double(i)+0.5-ngxm);
        *isvx = x;   *isvy = y;   *isvz = z;
        ++isvx;      ++isvy;      ++isvz;
      }
    }
  }
}

//====================================================================

bool restrict_test( const int activeDim,
                    const vector<int>& dimSrc, 
                    const vector<int>& dimDest )
{
  OperatorDatabase opDB;

  RestrictionAssembler<SVolField> ra( dimSrc, dimDest, true, true, true );
  opDB.register_new_operator< RestrictSVol >( new RestrictSVol(ra) );

  SVolField coarse( get_n_tot<SVolField>(dimSrc, true,true,true), get_ghost_set<SVolField>(dimSrc, true,true,true), NULL );
  SVolField fine  ( get_n_tot<SVolField>(dimDest,true,true,true), get_ghost_set<SVolField>(dimDest,true,true,true), NULL );
  SVolField fine2 ( get_n_tot<SVolField>(dimDest,true,true,true), get_ghost_set<SVolField>(dimDest,true,true,true), NULL );

  SVolField xcoarse( get_n_tot<SVolField>(dimSrc,true,true,true), get_ghost_set<SVolField>(dimSrc,true,true,true), NULL );
  SVolField ycoarse( get_n_tot<SVolField>(dimSrc,true,true,true), get_ghost_set<SVolField>(dimSrc,true,true,true), NULL );
  SVolField zcoarse( get_n_tot<SVolField>(dimSrc,true,true,true), get_ghost_set<SVolField>(dimSrc,true,true,true), NULL );

  SVolField xfine( get_n_tot<SVolField>(dimDest,true,true,true), get_ghost_set<SVolField>(dimDest,true,true,true), NULL );
  SVolField yfine( get_n_tot<SVolField>(dimDest,true,true,true), get_ghost_set<SVolField>(dimDest,true,true,true), NULL );
  SVolField zfine( get_n_tot<SVolField>(dimDest,true,true,true), get_ghost_set<SVolField>(dimDest,true,true,true), NULL );

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
  vector<int> dimSrc(3,1), dimDest(3,1);

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
