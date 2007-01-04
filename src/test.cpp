#include <SpatialOperator.h>
#include <SpatialField.h>
#include <FVStaggeredSpatialOps.h>
#include <LinearSystem.h>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_VectorOut.h>

#include <vector>
#include <iostream>
#include <cmath>

using namespace std;
using namespace SpatialOps;

#define PI 3.1415

void setup_geom( const std::vector<int> & dim,
		 const int nghost,
		 int & nn, int&nx, int&ny, int&nz,
		 std::vector<double> & spacing,
		 std::vector<double> & area,
		 double & volume )
{

  // set up array dimensions for "source" arrays
  nx = dim[0]+2*nghost;
  ny = dim[1]>1 ? dim[1]+2*nghost : 1;
  nz = dim[2]>1 ? dim[2]+2*nghost : 1;
  nn = nx*ny*nz;

  for( int i=0; i<3; ++i )
    spacing[i] = (dim[i]==1) ? 1.0 : PI/(dim[i]-1);
  area[0] = spacing[1]*spacing[2];
  area[1] = spacing[0]*spacing[2];
  area[2] = spacing[0]*spacing[1];

  volume = spacing[0]*spacing[1]*spacing[2];
}

//====================================================================

bool test_linsys()
{
  using namespace FVStaggeredUniform;

  std::vector<int> dim(3,1);
  dim[0] = 10;
  dim[1] = 2;
  dim[2] = 1;

  std::vector<double> spacing(3), area(3);
  double volume;
  int nn, nx,ny,nz;
  const int nghost = 1;
  setup_geom( dim, nghost, nn,nx,ny,nz, spacing, area, volume );

  // create a laplacian
  const Gradient2ndOrder xGrad( spacing, dim, X_DIR );
  const Divergence2ndOrder xDiv( area, volume, dim, X_DIR );
  ScratchOperator xLaplacian( dim, 3, X_DIR );

  xDiv.apply( xGrad, xLaplacian );

//   EpetraExt::RowMatrixToMatrixMarketFile( "xDiv.mm", xDiv.epetra_mat(), "", "" );
//   EpetraExt::RowMatrixToMatrixMarketFile( "xGrad.mm", xGrad.epetra_mat(), "", "" );
//   EpetraExt::RowMatrixToMatrixMarketFile( "xLaplacian.mm", xLaplacian.epetra_mat(), "", "" );

  LinearSystem linSys( dim );
  linSys.reset();

  linSys.get_lhs().add_contribution( xLaplacian );

  EpetraExt::RowMatrixToMatrixMarketFile( "xLaplacian.mm", linSys.get_lhs().epetra_mat(), "", "" );

  SpatialField rhsField ( dim, nghost, NULL, SpatialField::InternalStorage );

  rhsField = 1.0;

  RHS & rhs = linSys.get_rhs();

  rhs.reset();
  rhs.add_contribution( rhsField );

  linSys.solve();

  SpatialField tmp( dim, 0, linSys.get_rhs().get_ptr(), SpatialField::ExternalStorage );
  EpetraExt::VectorToMatrixMarketFile( "rhs.mm",   tmp.epetra_vec(), "", "" );
  EpetraExt::VectorToMatrixMarketFile( "soln.mm",  linSys.get_soln_field_epetra_vec(), "", "" );

  //  linSys.get_lhs().Print(std::cout);

  return true;
}

//====================================================================

bool test_field()
{
  std::vector<int> dim(3,1);
  dim[0] = 10;

  const int nghost = 0;

  int nn=1;
  for( int i=0; i<3; ++i )
    nn *= dim[i]+nghost;

  vector<double> fptr(nn,0.0);
  vector<double> gptr(nn,0.0);

  SpatialField f( dim, nghost, &fptr[0], SpatialField::ExternalStorage );
  SpatialField g( dim, nghost, &gptr[0], SpatialField::ExternalStorage );
  SpatialField h( dim, nghost, NULL,     SpatialField::InternalStorage );

  for( int i=0; i<nn; ++i ){
    fptr[i] = 2*i;
    gptr[i] = nn+i;
  }

  h += g;
  h *= f;
  h -= g;

  bool ok = true;

  for( int i=0; i<nn; ++i ){
    const double ans = gptr[i]*fptr[i]-gptr[i];
    if( std::abs( h.get_ptr()[i]-ans ) > 0.0001 )
      ok = false;
  }

  return ok;
}

//====================================================================

bool test_spatial_ops_x()
{
  using namespace FVStaggeredUniform;

  std::vector<int> dim(3,1);
  dim[0] = 10;
  dim[1] = 5;
  dim[2] = 3;

  // set up array dimensions for "source" arrays
  const int nghost=1;
  int nx, ny, nz, nn;
  std::vector<double> spacing(3,0.0);
  std::vector<double>    area(3,0.0);
  double volume;

  setup_geom( dim, 1, nn, nx, ny, nz, spacing, area, volume );

  // build some operators and stash them in the database.
  SpatialOpDatabase & SODatabase = SpatialOpDatabase::self();
  {
    SpatialOperator * xinterp = new LinearInterpolant( dim, X_DIR );
    SpatialOperator * xGrad   = new Gradient2ndOrder( spacing, dim, X_DIR );
    SpatialOperator * xDiv = new Divergence2ndOrder( area, volume, dim, X_DIR );
    SpatialOperator * scratchOp = new ScratchOperator( dim, 3, X_DIR );

    SODatabase.register_new_operator( SpatialOpDatabase::INTERPOLANT_X, xinterp, "X-Interpolant Second Order Staggered" );
    SODatabase.register_new_operator( SpatialOpDatabase::DIVERGENCE_X,  xDiv,    "X-Divergence Second Order Staggered" );
    SODatabase.register_new_operator( SpatialOpDatabase::GRADIENT_X,    xGrad,   "X-Gradient Second Order Staggered" );
    SODatabase.register_new_operator( SpatialOpDatabase::SCRATCH_X,     scratchOp, "Scratch X Second Order" );
  }

  SpatialOperator *& xinterp  = SODatabase.retrieve_operator( SpatialOpDatabase::INTERPOLANT_X );
  SpatialOperator *& xDiv     = SODatabase.retrieve_operator( SpatialOpDatabase::DIVERGENCE_X  );
  SpatialOperator *& xGrad    = SODatabase.retrieve_operator( SpatialOpDatabase::GRADIENT_X    );
  SpatialOperator *& scratchOp= SODatabase.retrieve_operator( "Scratch X Second Order"         );

  /*
  EpetraExt::RowMatrixToMatrixMarketFile( "Int_x.mm",
					  xinterp->epetra_mat(), 
					  "xinterp",
					  "1-D interpolant in x-direction" );
  */
  EpetraExt::RowMatrixToMatrixMarketFile( "Grad_x.mm",
					  xGrad->epetra_mat(), 
					  "xGrad",
					  "1-D gradient in x-direction" );
  EpetraExt::RowMatrixToMatrixMarketFile( "Div_x.mm",
					  xDiv->epetra_mat(), 
					  "xDiv",
					  "1-D divergence in x-direction" );

  vector<double> fptr  (nn,0.0);
  vector<double> dfdx  (nn,0.0);
  vector<double> d2fdx2(nn,0.0);
  vector<double> xptr  (nn,0.0);

  // set values in fptr and xptr
  int ix=0;
  for( int k=0; k<nz; ++k ){
    for( int j=0; j<ny; ++j){
      for( int i=0; i<nx; ++i ){
	xptr[ix] = double(i)*spacing[0];
	fptr[ix] = std::sin(xptr[ix]);
	dfdx[ix] = std::cos(xptr[ix]);
	++ix;
      }
    }
  }

  SpatialField  f( dim, nghost, &fptr[0], SpatialField::ExternalStorage );
  SpatialField  x( dim, nghost, &xptr[0], SpatialField::ExternalStorage );
  SpatialField  g( dim, nghost, NULL,     SpatialField::InternalStorage );
  SpatialField df( dim, nghost, NULL,     SpatialField::InternalStorage );
  SpatialField xg( dim, nghost, NULL,     SpatialField::InternalStorage );
  SpatialField d2f(dim, nghost, NULL,     SpatialField::InternalStorage );

  xinterp->apply( f, g );
  xinterp->apply( x, xg);
  xGrad  ->apply( f, df );

  // form the laplacian
  xDiv->apply( *xGrad, *scratchOp );
  scratchOp->apply( f, d2f );

  EpetraExt::RowMatrixToMatrixMarketFile( "Laplace_x.mm", scratchOp->epetra_mat(),  "laplace x", "1-D laplacian in x-direction" );
  EpetraExt::VectorToMatrixMarketFile( "fx.mm",  f.epetra_vec(), "f", "original data" );
  EpetraExt::VectorToMatrixMarketFile(  "x.mm",   x.epetra_vec(), "x", "original x" );
  EpetraExt::VectorToMatrixMarketFile( "gx.mm",  g.epetra_vec(), "g", "intpolated data" );
  EpetraExt::VectorToMatrixMarketFile( "dfdx.mm", df.epetra_vec(), "dfdx", "grad x" );
  EpetraExt::VectorToMatrixMarketFile( "xg.mm", xg.epetra_vec(), "xg", "intpolated x" );
  EpetraExt::VectorToMatrixMarketFile( "d2fdx2.mm", d2f.epetra_vec(), "d2f", "" );
 
  /*
    scratchOp.reset_entries(0.0);
    scratchOp += xinterp;
    
    EpetraExt::RowMatrixToMatrixMarketFile( "I2.mm", scratchOp->epetra_mat(), "", "" );
    EpetraExt::RowMatrixToMatrixMarketFile( "I1.mm",   xinterp->epetra_mat(), "", "" );
  */

  return true;
}

//====================================================================

bool test_spatial_ops_y()
{
  using namespace FVStaggeredUniform;

  std::vector<int> dim(3,1);
  dim[0] = 2;
  dim[1] = 10;
  dim[2] = 2;

  const int nghost=1;
  int nn,nx,ny,nz;
  std::vector<double> spacing(3,0.0);
  std::vector<double>    area(3,0.0);
  double volume;
  setup_geom( dim, 1, nn,nx,ny,nz, spacing, area, volume );

  LinearInterpolant yinterp( dim, Y_DIR );
  Gradient2ndOrder yGrad( spacing, dim, Y_DIR );
  /*
  EpetraExt::RowMatrixToMatrixMarketFile( "Int_y.mm",
					  yinterp.epetra_mat(), 
					  "yinterp",
					  "1-D interpolant in y-direction" );

  EpetraExt::RowMatrixToMatrixMarketFile( "Grad_y.mm",
					  yGrad.epetra_mat(), 
					  "yGrad",
					  "1-D gradient in y-direction" );
  */
  vector<double> fptr(nn,0.0);
  vector<double> dfdx(nn,0.0);
  vector<double> yptr(nn,0.0);

  // set values in fptr and yptr
  int ix=0;
  for( int k=0; k<nz; ++k ){
    for( int j=0; j<ny; ++j){
      for( int i=0; i<nx; ++i ){
	yptr[ix] = double(j)*spacing[1];
	fptr[ix] = std::sin(yptr[ix]);
	++ix;
      }
    }
  }

  SpatialField  f( dim, nghost, &fptr[0], SpatialField::ExternalStorage );
  SpatialField  y( dim, nghost, &yptr[0], SpatialField::ExternalStorage );
  SpatialField df( dim, nghost, NULL,     SpatialField::InternalStorage );
  SpatialField  g( dim, nghost, NULL,     SpatialField::InternalStorage );
  SpatialField yg( dim, nghost, NULL,     SpatialField::InternalStorage );
  
  yinterp.apply( f, g );
  yinterp.apply( y, yg);
  yGrad.apply( f, df );
  
  EpetraExt::VectorToMatrixMarketFile( "fy.mm",  f.epetra_vec(),  "f", "original data"   );
  EpetraExt::VectorToMatrixMarketFile( "y.mm",   y.epetra_vec(),  "y", "original y"      );
  EpetraExt::VectorToMatrixMarketFile( "dfdy.mm",  df.epetra_vec(),  "dfdy", "gradient y" );
  EpetraExt::VectorToMatrixMarketFile( "gy.mm",  g.epetra_vec(),  "g", "intpolated data" );
  EpetraExt::VectorToMatrixMarketFile( "yg.mm", yg.epetra_vec(), "yg", "intpolated y"    );

  bool ok = true;
  return ok;
}

//====================================================================

bool test_spatial_ops_z()
{
  using namespace FVStaggeredUniform;

  std::vector<int> dim(3,1);
  dim[0] = 2;
  dim[1] = 2;
  dim[2] = 10;

  // set up array dimensions for "source" arrays
  const int nghost=1;
  const int nx = dim[0]+2*nghost;
  const int ny = dim[1]>1 ? dim[1]+2*nghost : 1;
  const int nz = dim[2]>1 ? dim[2]+2*nghost : 1;
  const int nn = nx*ny*nz;

  std::vector<double> spacing(3,0.0);
  for( int i=0; i<3; ++i ) spacing[i] = PI/(dim[i]-1);

  LinearInterpolant zinterp( dim, Z_DIR );
  Gradient2ndOrder zGrad( spacing, dim, Z_DIR );
  /*
  EpetraExt::RowMatrixToMatrixMarketFile( "Int_z.mm",
					  zinterp.epetra_mat(), 
					  "zinterp",
					  "1-D interpolant in z-direction" );

  EpetraExt::RowMatrixToMatrixMarketFile( "Grad_z.mm",
					  zGrad.epetra_mat(), 
					  "zinterp",
					  "1-D interpolant in z-direction" );
  */
  vector<double> fptr(nn,0.0);
  vector<double> zptr(nn,0.0);

  // set values in fptr and zptr
  int ix=0;
  for( int k=0; k<nz; ++k ){
    for( int j=0; j<ny; ++j){
      for( int i=0; i<nx; ++i ){
	zptr[ix] = double(k)*spacing[2];
	fptr[ix] = std::sin(zptr[ix]);
	++ix;
      }
    }
  }

  SpatialField  f( dim, nghost, &fptr[0], SpatialField::ExternalStorage );
  SpatialField  z( dim, nghost, &zptr[0], SpatialField::ExternalStorage );
  SpatialField df( dim, nghost, NULL,     SpatialField::InternalStorage );
  SpatialField  g( dim, nghost, NULL,     SpatialField::InternalStorage );
  SpatialField zg( dim, nghost, NULL,     SpatialField::InternalStorage );
  
  zinterp.apply( f, g );
  zinterp.apply( z, zg);
  zGrad.apply( f, df );
  
  EpetraExt::VectorToMatrixMarketFile( "fz.mm",  f.epetra_vec(),  "f", "original data"   );
  EpetraExt::VectorToMatrixMarketFile( "z.mm",   z.epetra_vec(),  "z", "original z"      );
  EpetraExt::VectorToMatrixMarketFile( "dfdz.mm",  df.epetra_vec(),  "dfdz", "gradient z" );
  EpetraExt::VectorToMatrixMarketFile( "gz.mm",  g.epetra_vec(),  "g", "intpolated data" );
  EpetraExt::VectorToMatrixMarketFile( "zg.mm", zg.epetra_vec(), "zg", "intpolated z"    );

  bool ok = true;
  return ok;
}

//====================================================================

bool test_spatial_ops_algebra()
{
  using namespace FVStaggeredUniform;

  bool ok = true;

  std::vector<int> dim(3,1);
  dim[0] = 2;
  dim[1] = 2;
  dim[2] = 3;

  // set up array dimensions for "source" arrays
  const int nghost=1;
  const int nx = dim[0]+2*nghost;
  const int ny = dim[1]>1 ? dim[1]+2*nghost : 1;
  const int nz = dim[2]>1 ? dim[2]+2*nghost : 1;
  const int nn = nx*ny*nz;

  LinearInterpolant zinterp( dim, Z_DIR );
  LinearInterpolant z2( dim, Z_DIR );

  vector<double> f1(nn,0.0);
  vector<double> f2(nn,0.0);

  for( int ix=0; ix<nn; ++ix ){
    f1[ix] = 2.0;
    f2[ix] = ix;
  }

  SpatialField field1( dim, nghost, &f1[0], SpatialField::InternalStorage );
  SpatialField field2( dim, nghost, &f2[0], SpatialField::InternalStorage );

  zinterp.right_scale( field1 );
  for( int i=0; i<zinterp.nrows(); ++i ){
    double*vals;    int*ixs;    int nentries;
    zinterp.epetra_mat().ExtractMyRowView( i, nentries, vals, ixs );
    assert( nentries == 2 );

    for( int k=0; k<nentries; ++k ){
      if( *vals++ != 1.0 ) ok=false;
    }
  }

  zinterp.right_scale( field2 );
  for( int i=0; i<zinterp.nrows(); ++i ){
    double*vals;    int*ixs;    int nentries;
    zinterp.epetra_mat().ExtractMyRowView( i, nentries, vals, ixs );
    assert( nentries == 2 );

    for( int k=0; k<nentries; ++k ){
      if( *vals++ != ixs[k] ){
	ok=false;
 	cout << *(vals-1) << "  " << ixs[k] << std::endl;
      }
    }    
  }

  return ok;
}

//====================================================================

int main()
{
  bool ok;
  ok = test_linsys();
  if( ok ) cout << "   linsys test:   PASS" << endl;
  else     cout << "   linsys test:   FAIL" << endl;

  ok = test_field();
  if( ok ) cout << "   field ops test:   PASS" << endl;
  else     cout << "   field ops test:   FAIL" << endl;

  ok = test_spatial_ops_x();
  if( ok ) cout << "   spatial ops X test: PASS" << endl;
  else     cout << "   spatial ops X test: FAIL" << endl;

  ok = test_spatial_ops_y();
  if( ok ) cout << "   spatial ops Y test: PASS" << endl;
  else     cout << "   spatial ops Y test: FAIL" << endl;

  ok = test_spatial_ops_z();
  if( ok ) cout << "   spatial ops Z test: PASS" << endl;
  else     cout << "   spatial ops Z test: FAIL" << endl;

  ok = test_spatial_ops_algebra();
  cout << "   spatial ops algebra test: ";
  if( ok ) cout << "PASS" << std::endl;
  else     cout << "FAIL" << std::endl;

  return 0;
}

//====================================================================
