#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;


#include <spatialops/structured/FVStaggered.h>
#include <spatialops/structured/FVStaggeredBCTools.h>
#include <spatialops/LinearSystem.h>
#include <spatialops/FieldFunctions.h>

#include "Grid.h"
#include "Functions.h"
#include "buildOps.h"

using namespace SpatialOps;
using namespace structured;

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
           const IntVec& dim,
           const std::vector<bool>& bcFlag )
    : FieldFunction3D<FieldT>( x, y, z ),
      a(2.0), b(2.0), c(2.0), d(1.0),
      dim_(dim),
      tmp( get_window_with_ghost<FieldT>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
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
  const IntVec dim_;
  mutable FieldT tmp;
};

//--------------------------------------------------------------------

double test_poisson( const OperatorDatabase& opDB,
                     const Grid& grid,
                     const IntVec& dim,
                     const vector<bool>& bcFlag,
                     const BCType bcType )
{
  ScratchSVol& Lx = *opDB.retrieve_operator<ScratchSVol>( 1 );
  ScratchSVol& Ly = *opDB.retrieve_operator<ScratchSVol>( 2 );
  ScratchSVol& Lz = *opDB.retrieve_operator<ScratchSVol>( 3 );

  LinSysInfo lsi( dim, bcFlag[0], bcFlag[1], bcFlag[2] );
  LinearSystem& linsys = LinSysFactory::self().get_linsys( lsi );
  RHS& rhs = linsys.get_rhs();
  LHS& lhs = linsys.get_lhs();
  lhs.reset();

  const GradSVolSSurfX& Gx = *opDB.retrieve_operator<GradSVolSSurfX>();
  const GradSVolSSurfY& Gy = *opDB.retrieve_operator<GradSVolSSurfY>();
  const GradSVolSSurfZ& Gz = *opDB.retrieve_operator<GradSVolSSurfZ>();

  const DivSSurfXSVol& Dx = *opDB.retrieve_operator<DivSSurfXSVol>();  
  const DivSSurfYSVol& Dy = *opDB.retrieve_operator<DivSSurfYSVol>();  
  const DivSSurfZSVol& Dz = *opDB.retrieve_operator<DivSSurfZSVol>();  

  const InterpSVolSSurfX& Rx = *opDB.retrieve_operator<InterpSVolSSurfX>();
  const InterpSVolSSurfY& Ry = *opDB.retrieve_operator<InterpSVolSSurfY>();
  const InterpSVolSSurfZ& Rz = *opDB.retrieve_operator<InterpSVolSSurfZ>();

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
  SVolField rhsField( get_window_with_ghost<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
                      NULL );
  SVolField tmpField( get_window_with_ghost<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]),
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

    SSurfXField bcValX( get_window_with_ghost<SSurfXField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), NULL );
    SSurfYField bcValY( get_window_with_ghost<SSurfYField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), NULL );
    SSurfZField bcValZ( get_window_with_ghost<SSurfZField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), NULL );

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

    // set bcs: x faces
    if( dim[0]>1 ){
      for( int ix=0; ix<2; ++ix ){
        const int i = ix==0 ? 0 : dim[0]-1;
        const BCSide side = ix==0 ? X_MINUS_SIDE : X_PLUS_SIDE;
        for( int j=0; j<dim[1]; ++j ){
          for( int k=0; k<dim[2]; ++k ){
            // determine index for x field.
            const IntVec ijk = shift_to_ghost_ix<InterpSVolSSurfX,ScratchSVol>( dim, side, IntVec(i,j,k) );
            const int iix = get_index_with_ghost<SSurfXField>(dim,bcFlag[0],bcFlag[1],bcFlag[2],ijk);
            const int irhs = i + j*dim[0] + k*dim[0]*dim[1];
            const double bcval = bcValX[iix];
            switch ( bcType ){
            case DIRICHLET:
              imprint_bc_on_op<InterpSVolSSurfX,ScratchSVol>( Rx, IntVec(i,j,k), dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, side, Lx, rhs[irhs] );
              break;
            case NEUMANN:
              imprint_bc_on_op<GradSVolSSurfX,  ScratchSVol>( Gx, IntVec(i,j,k), dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, side, Lx, rhs[irhs] );
              break;
            }
          }
        }
      }
    }

    // set bcs: y faces
    if( dim[1]>1 ){
      for( int iy=0; iy<2; ++iy ){
        const int j = iy==0 ? 0 : dim[1]-1;
        const BCSide side = iy==0 ? Y_MINUS_SIDE : Y_PLUS_SIDE;
        for( int i=0; i<dim[0]; ++i ){
          for( int k=0; k<dim[2]; ++k ){
            // obtain the BC value
            const IntVec ijk = shift_to_ghost_ix<InterpSVolSSurfY,ScratchSVol>( dim, side, IntVec(i,j,k) );
            const int iix = get_index_with_ghost<SSurfYField>(dim,bcFlag[0],bcFlag[1],bcFlag[2],ijk);
            const int irow = i + j*dim[0] + k*dim[0]*dim[1];
            const double bcval = bcValY[iix];
            // set the BC value:
            switch ( bcType ){
            case DIRICHLET:
              imprint_bc_on_op<InterpSVolSSurfY,ScratchSVol>( Ry, IntVec(i,j,k), dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, side, Ly, rhs[irow] );
              break;
            case NEUMANN:
              imprint_bc_on_op<GradSVolSSurfY, ScratchSVol>( Gy, IntVec(i,j,k), dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, side, Ly, rhs[irow] );
              break;
            }
          }
        }
      }
    }

    // set bcs: z faces
    if( dim[2]>1 ){
      for( int iz=0; iz<2; ++iz ){
        const int k = iz==0 ? 0 : dim[2]-1;
        const BCSide side = iz==0 ? Z_MINUS_SIDE : Z_PLUS_SIDE;
        for( int i=0; i<dim[0]; ++i ){
          for( int j=0; j<dim[1]; ++j ){
            const IntVec ijk = shift_to_ghost_ix<InterpSVolSSurfZ,ScratchSVol>( dim, side, IntVec(i,j,k) );
            const int iix = get_index_with_ghost<SSurfZField>(dim,bcFlag[0],bcFlag[1],bcFlag[2],ijk);
            const int irow = i + j*dim[0] + k*dim[0]*dim[1];
            const double bcval = bcValZ[iix];
            switch ( bcType ){
            case DIRICHLET:
              imprint_bc_on_op<InterpSVolSSurfZ,ScratchSVol>( Rz, IntVec(i,j,k), dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, side, Lz, rhs[irow] );
              break;
            case NEUMANN:
              imprint_bc_on_op<GradSVolSSurfZ,  ScratchSVol>( Gz, IntVec(i,j,k), dim, bcFlag[0], bcFlag[1], bcFlag[2], bcval, side, Lz, rhs[irow] );
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
    SVolField phi( get_window_with_ghost<SVolField>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), NULL );
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

      SVolRHS tmp( get_window_with_ghost<SVolRHS>(dim,bcFlag[0],bcFlag[1],bcFlag[2]), NULL );

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

bool driver( const std::vector<int>& dim,
             const std::vector<double>& length,
             const double tol1, const double tol2 )
{
  std::vector<double> spacing(3,1.0);
  for( int i=0; i<3; ++i ){
    if( dim[i]>1 ) spacing[i] = length[i]/dim[i];
  }

  std::vector<bool> bcFlag(3,true);
  OperatorDatabase opDB;
  build_ops( dim, spacing, bcFlag, opDB );
  const Grid grid( dim, spacing, bcFlag, opDB );

  const double err1 = test_poisson( opDB, grid, dim, bcFlag, DIRICHLET );
  cout << "Poisson eqn: Dirichlet BC ... ";
  if( err1 > tol1 ){
    cout << endl << scientific << "*** FAIL *** (" << err1 << ")" << endl;
    cout.unsetf(ios::scientific);
  }
  else{
    cout << "PASS" << endl;
  }

  bool isFailed = true;
  const double err2 = test_poisson( opDB, grid, dim, bcFlag, NEUMANN   );
  cout << "Poisson eqn: Neumann   BC ... ";
  if( err2 > tol1 ){
    cout << endl << scientific << "*** FAIL *** (" << err2 << ")" << endl;
    cout.unsetf(ios::scientific);
  }
  else{
    cout << "PASS" << endl;
    isFailed = false;
  }
  return isFailed;
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

  return driver( dim, length, 1.12e-3, 9.4e-4 );

//   dim[0]=30;      dim[1]=30;   dim[2]=30;
//   length[0]=1; length[1]=1; length[2]=1;
//   driver( dim, length, 1e-5, 1e-5 );
}
