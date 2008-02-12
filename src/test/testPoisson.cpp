#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;


#include <FVStaggered.h>
#include <FVStaggeredBCTools.h>
#include <LinearSystem.h>
#include <FieldFunctions.h>

#include <Grid.h>
#include <Functions.h>
#include <buildOps.h>

using namespace SpatialOps;
using namespace FVStaggered;

enum BCType
{
  DIRICHLET,
  NEUMANN
};

//--------------------------------------------------------------------

template<typename FieldT>
class QuadFun : public FieldFunction3D<FieldT>
{
public:
  QuadFun( const FieldT& x,
	   const FieldT& y,
	   const FieldT& z,
	   const std::vector<int>& dim,
	   const std::vector<bool>& bcFlag )
    : FieldFunction3D<FieldT>( x, y, z ),
      a(2.0), b(2.0), c(2.0), d(1.0),
      dim_(dim),
      tmp( get_n_tot<FieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
	   get_ghost_set<FieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
	   NULL )
  {}

  ~QuadFun(){}

  void evaluate( FieldT& phi ) const
  {
    phi = d;
    if( dim_[0]>1 ){
      tmp  = this->get_x();
      tmp *= this->get_x();
      tmp *= a;
      phi += tmp;
    }
    if( dim_[1]>1 ){
      tmp  = this->get_y();
      tmp *= this->get_y();
      tmp *= b;
      phi += tmp;
    }
    if( dim_[2]>1 ){
      tmp  = this->get_z();
      tmp *= this->get_z();
      tmp *= c;
      phi += tmp;
    }
  }

  void dx( FieldT& df ) const{
    df = 0.0;
    if( dim_[0]>1 ){
      df = this->get_x();  df*=2*a;
    }
  }

  void dy( FieldT& df ) const{
    df = 0.0;
    if( dim_[1]>1 ){
      df = this->get_y();
      df*=2*b;
    }
  }

  void dz( FieldT& df ) const{
    df = 0.0;
    if( dim_[2]>1 ){
      df = this->get_z();
      df*=2*c;
    }
  }

  void d2x( FieldT& d2f ) const{
    d2f=0.0;
    if( dim_[0]>1 ) d2f=2*a;
  }

  void d2y( FieldT& d2f ) const{
    d2f=0.0;
    if( dim_[1]>1 ) d2f=2*b;
  }

  void d2z( FieldT& d2f ) const{
    d2f=0.0;
    if( dim_[2]>1 ) d2f=2*c;
  }

private:
  const double a, b, c, d;
  const std::vector<int> dim_;
  mutable FieldT tmp;
};

//--------------------------------------------------------------------

double test_poisson( const Grid& grid,
		     const vector<int>& dim,
		     const vector<bool>& bcFlag,
		     const BCType bcType )
{
  ScratchSVol& Lx = *SpatialOpDatabase<ScratchSVol>::self().retrieve_operator( 1 );
  ScratchSVol& Ly = *SpatialOpDatabase<ScratchSVol>::self().retrieve_operator( 2 );
  ScratchSVol& Lz = *SpatialOpDatabase<ScratchSVol>::self().retrieve_operator( 3 );

  LinSysInfo lsi( dim, bcFlag[0], bcFlag[1], bcFlag[2] );
  LinearSystem& linsys = LinSysFactory::self().get_linsys( lsi );
  RHS& rhs = linsys.get_rhs();
  LHS& lhs = linsys.get_lhs();
  lhs.reset();

  const GradSVolSSurfX& Gx = *SpatialOpDatabase<GradSVolSSurfX>::self().retrieve_operator();
  const GradSVolSSurfY& Gy = *SpatialOpDatabase<GradSVolSSurfY>::self().retrieve_operator();
  const GradSVolSSurfZ& Gz = *SpatialOpDatabase<GradSVolSSurfZ>::self().retrieve_operator();

  const DivSSurfXSVol& Dx = *SpatialOpDatabase<DivSSurfXSVol>::self().retrieve_operator();  
  const DivSSurfYSVol& Dy = *SpatialOpDatabase<DivSSurfYSVol>::self().retrieve_operator();  
  const DivSSurfZSVol& Dz = *SpatialOpDatabase<DivSSurfZSVol>::self().retrieve_operator();  

  const InterpSVolSSurfX& Rx = *SpatialOpDatabase<InterpSVolSSurfX>::self().retrieve_operator();
  const InterpSVolSSurfY& Ry = *SpatialOpDatabase<InterpSVolSSurfY>::self().retrieve_operator();
  const InterpSVolSSurfZ& Rz = *SpatialOpDatabase<InterpSVolSSurfZ>::self().retrieve_operator();

  //
  // set up the Laplacian operator in each direction and assemble the
  // linear system to be solved.
  //
  if( dim[0]>1 )  Dx.apply_to_op( Gx, Lx );
  if( dim[1]>1 )  Dy.apply_to_op( Gy, Ly );
  if( dim[2]>1 )  Dz.apply_to_op( Gz, Lz );

//   Gx.write_matlab("Gx"); Dx.write_matlab("Dx"); Lx.write_matlab("Lxpbc");
//   Gy.write_matlab("Gy"); Dy.write_matlab("Dy");
//   Gz.write_matlab("Gz"); Dz.write_matlab("Dz");

  //
  // set the RHS field
  //
  SVolField rhsField( get_n_tot<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
		      get_ghost_set<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
		      NULL );
  SVolField tmpField( get_n_tot<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
		      get_ghost_set<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
		      NULL );
  QuadFun<SVolField> quadFunVol( grid.xcoord_svol(), grid.ycoord_svol(), grid.zcoord_svol(), dim, bcFlag );
  rhsField = 0.0;
  if( dim[0]>1 ){ quadFunVol.d2x(rhsField); }
  if( dim[1]>1 ){ quadFunVol.d2y(tmpField); rhsField += tmpField; }
  if( dim[2]>1 ){ quadFunVol.d2z(tmpField); rhsField += tmpField; }

  rhs.reset();
  rhs.add_field_contribution( rhsField );

  // set bcs on scratch matrix that will form the linear system.
  {
    QuadFun<SSurfXField> bcFunX( grid.xcoord_sxsurf(), grid.ycoord_sxsurf(), grid.zcoord_sxsurf(), dim, bcFlag );
    QuadFun<SSurfYField> bcFunY( grid.xcoord_sysurf(), grid.ycoord_sysurf(), grid.zcoord_sysurf(), dim, bcFlag );
    QuadFun<SSurfZField> bcFunZ( grid.xcoord_szsurf(), grid.ycoord_szsurf(), grid.zcoord_szsurf(), dim, bcFlag );

    SSurfXField bcValX( get_n_tot<SSurfXField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), get_ghost_set<SSurfXField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), NULL );
    SSurfYField bcValY( get_n_tot<SSurfYField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), get_ghost_set<SSurfYField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), NULL );
    SSurfZField bcValZ( get_n_tot<SSurfZField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), get_ghost_set<SSurfZField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), NULL );

    switch ( bcType ){
    case DIRICHLET:
      bcFunX.evaluate( bcValX );
      bcFunY.evaluate( bcValY );
      bcFunZ.evaluate( bcValZ );
      break;
    case NEUMANN:
      bcFunX.dx( bcValX );
      bcFunY.dy( bcValY );
      bcFunZ.dz( bcValZ );
      break;

    default:
      assert(0);
    }

    //
    // set the boundary conditions
    //
    const int ighost = dim[0]>1 ? SSurfXField::Ghost::NM : 0;
    const int jghost = dim[1]>1 ? SSurfYField::Ghost::NM : 0;
    const int kghost = dim[2]>1 ? SSurfZField::Ghost::NM : 0;

    // set bcs: x faces
    if( dim[0]>1 ){
      for( int ix=0; ix<2; ++ix ){
	int i=0, ii=0;
	if( ix!=0 ){
	  // set index on bf field, which has different dimensionality than the rhs field
	  ii = get_nx<SSurfXField>(dim,bcFlag[0]) - SSurfXField::Ghost::NP - SSurfXField::Ghost::NM - 1;
	  // index for the rhs field.
	  i = dim[0]-1;
	}
	for( int j=0; j<dim[1]; ++j ){
	  for( int k=0; k<dim[2]; ++k ){
	    // determine index for x field.
 	    const IndexTriplet ijk( ii+ighost, j+jghost, k+kghost );
 	    const int iix = ijk2flat<SSurfXField>::value(dim,ijk,bcFlag[0],bcFlag[1],bcFlag[2]);
	    const double bcval = bcValX[iix];
	    const int irow = i + j*dim[0] + k*dim[0]*dim[1];
	    switch ( bcType ){
	    case DIRICHLET:
	      imprint_bc_on_op<InterpSVolSSurfX,XDIR,ScratchSVol>( Rx, i, j, k, dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, rhsField, Lx, rhs[irow] );
	      break;
	    case NEUMANN:
	      imprint_bc_on_op<GradSVolSSurfX,  XDIR,ScratchSVol>( Gx, i, j, k, dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, rhsField, Lx, rhs[irow] );
	      break;
	    }
	  }
	}
      }
    }

    // set bcs: y faces
    if( dim[1]>1 ){
      for( int iy=0; iy<2; ++iy ){
	int j=0, jj=0;
	if( iy!=0 ){
	  j = dim[1]-1;
	  jj = get_ny<SSurfYField>(dim,bcFlag[1]) - SSurfYField::Ghost::NP - SSurfYField::Ghost::NM - 1;
	}
	for( int i=0; i<dim[0]; ++i ){
	  for( int k=0; k<dim[2]; ++k ){

	    // obtain the BC value
	    const IndexTriplet ijk( i+ighost, jj+jghost, k+kghost );
	    const int iix = ijk2flat<SSurfYField>::value(dim,ijk,bcFlag[0],bcFlag[1],bcFlag[2]);
	    const int irow = i + j*dim[0] + k*dim[0]*dim[1];
	    const double bcval = bcValY[iix];

	    // set the BC value:
	    switch ( bcType ){
	    case DIRICHLET:
	      imprint_bc_on_op<InterpSVolSSurfY,YDIR,ScratchSVol>( Ry, i, j, k, dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, rhsField, Ly, rhs[irow] );
	      break;
	    case NEUMANN:
	      imprint_bc_on_op<GradSVolSSurfY,  YDIR,ScratchSVol>( Gy, i, j, k, dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, rhsField, Ly, rhs[irow] );
	      break;
	    }
	  }
	}
      }
    }

    // set bcs: z faces
    if( dim[2]>1 ){
      for( int iz=0; iz<2; ++iz ){
	int k=0, kk=0;
	if( iz!=0 ){
	  k = dim[2]-1;
	  kk = get_nz<SSurfZField>(dim,bcFlag[2]) - SSurfZField::Ghost::NM - SSurfZField::Ghost::NP - 1;
	}
	for( int i=0; i<dim[0]; ++i ){
	  for( int j=0; j<dim[1]; ++j ){
	    const IndexTriplet ijk( i+ighost, j+jghost, kk+kghost );
	    const int iix = ijk2flat<SSurfZField>::value(dim,ijk,bcFlag[0],bcFlag[1],bcFlag[2]);
	    const double bcval = bcValZ[iix];
	    const int irow = i + j*dim[0] + k*dim[0]*dim[1];
	    switch ( bcType ){
	    case DIRICHLET:
	      imprint_bc_on_op<InterpSVolSSurfZ,ZDIR,ScratchSVol>( Rz, i, j, k, dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, rhsField, Lz, rhs[irow] );
	      break;
	    case NEUMANN:
	      imprint_bc_on_op<GradSVolSSurfZ,  ZDIR,ScratchSVol>( Gz, i, j, k, dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, rhsField, Lz, rhs[irow] );
	      break;
	    }
	  }
	}
      }
    }
  }

  if( dim[0]>1 )  lhs.add_op_contribution( Lx );
  if( dim[1]>1 )  lhs.add_op_contribution( Ly );
  if( dim[2]>1 )  lhs.add_op_contribution( Lz );

//   Lx.write_matlab("Lx");
//   Ly.write_matlab("Ly");
//   Lz.write_matlab("Lz");

  //
  // Solve the linear system for the solution.
  //
  linsys.solve();


  {
    //
    // examine the solution to determine error
    //
    SVolField phi( get_n_tot<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), get_ghost_set<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), NULL );
    quadFunVol.evaluate( phi );


    const SOLN& soln = linsys.get_soln_field();

    SOLN::const_interior_iterator isoln = soln.interior_begin();
    SVolField::interior_iterator iphi = phi.interior_begin();
    const SVolField::interior_iterator iphie = phi.interior_end();

    double maxAbsErr=0.0, maxRelErr=0.0;
    double avgAbsErr=0.0, avgRelErr=0.0;
    int nrel=0, nabs=0;
    for( ; iphi!=iphie; ++isoln, ++iphi ){
      const double err = abs(*iphi-*isoln);
      avgAbsErr += err;
      maxAbsErr = max(err,maxAbsErr);
      ++nabs;
      if( abs(*iphi)>1e-10 ){
	const double relErr = abs(err / *iphi);
	avgRelErr += relErr;
	maxRelErr = max(relErr,maxRelErr);
	++nrel;
      }
    }
    avgRelErr /= double(nrel);
    avgAbsErr /= double(nabs);


    /*
      cout << endl
      << scientific << setprecision(2)
      << "  max abs error: " << setw(8) << maxAbsErr << "  max rel err: " << setw(8) << maxRelErr << endl
      << "  avg abs error: " << setw(8) << avgAbsErr << "  avg rel err: " << setw(8) << avgRelErr << endl << endl;
      cout.unsetf(ios::scientific | ios::floatfield );

      SVolRHS tmp( get_n_tot<SVolRHS>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), get_ghost_set<SVolRHS>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), NULL );

      tmp = rhs;  tmp.write_matlab("rhs");

      rhs.reset();
      rhs.add_field_contribution( quadFunVol.get_x() );
      tmp = rhs;
      tmp.write_matlab("x");
      
      rhs.reset();
      rhs.add_field_contribution( quadFunVol.get_y() );
      tmp = rhs;
      tmp.write_matlab("y");

      rhs.reset();
      rhs.add_field_contribution( phi );
      tmp = rhs;
      tmp.write_matlab("phiExact");    

      tmp = linsys.get_soln_field();
      tmp.write_matlab("phi");

      linsys.get_lhs().write_matlab("A");
    */

    return maxAbsErr;
  }
}


//--------------------------------------------------------------------

void driver( const std::vector<int>& dim,
	     const std::vector<double>& length,
	     const double tol1, const double tol2 )
{
  std::vector<double> spacing(3,1.0);
  for( int i=0; i<3; ++i ){
    if( dim[i]>1 ) spacing[i] = length[i]/dim[i];
  }

  std::vector<bool> bcFlag(3,true);
  build_ops( dim, spacing, bcFlag );
  const Grid grid( dim, spacing, bcFlag );

  const double err1 = test_poisson( grid, dim, bcFlag, DIRICHLET );
  const double err2 = test_poisson( grid, dim, bcFlag, NEUMANN   );

  cout << "Poisson eqn: Dirichlet BC ... ";
  if( err1 > tol1 ){
    cout << endl << scientific << "FAIL (" << err1 << ")" << endl;
    cout.unsetf(ios::scientific);
  }
  else{
    cout << "PASS" << endl;
  }

  cout << "Poisson eqn: Neumann   BC ... ";
  if( err2 > tol1 ){
    cout << endl << scientific << "FAIL (" << err2 << ")" << endl;
    cout.unsetf(ios::scientific);
  }
  else{
    cout << "PASS" << endl;
  }

}

//--------------------------------------------------------------------

int main()
{
  vector<int> dim(3,1);
  dim[0] = 21 ;
  dim[1] = 46 ;
  dim[2] = 12;

  std::vector<double> length(3,1);
  length[1]=2;
  length[2]=0.5;

  driver( dim, length, 1.12e-3, 9.4e-4 );

//   dim[0]=30;      dim[1]=30;   dim[2]=30;
//   length[0]=1; length[1]=1; length[2]=1;
//   driver( dim, length, 1e-5, 1e-5 );
}
