#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

#include <FV2ndOrderTypes.h> 

#define TOL 1.e-6
#define MAX_MAXRELERR 1.e-3
#define MAX_MEANRELERR 1.e-3

#define X_DIR 0
#define Y_DIR 1
#define Z_DIR 2
 
using namespace SpatialOps;
using namespace FVStaggeredUniform;

#include <test_functions.h>
#include <grid.h>
#include <analytical_class.h>
#include <utilities.h>


// *****************************************************************************************************************************	
// Test Linear Interpolant Operator C2F
// *****************************************************************************************************************************	
bool test_linear_interpolant_C2F(const grid_class& grid, analytical_class& funct, std::vector<double> &mean_rel_err, std::vector<double> &max_rel_err)
{
  using namespace FVStaggeredUniform;
  std::cout.setf(ios::scientific);

  bool test_x, test_y, test_z, ok;
	
  std::cout << "Test Function: " << funct.test.description() << endl;
		
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Build the operators
  // -----------------------------------------------------------------------------------------------------------------------------	
  InterpXC2F::Assembler Rx_C2F_Assembler(grid.dim);				// Assembler Linear Interpolant X - Cells To Faces
  InterpYC2F::Assembler Ry_C2F_Assembler(grid.dim);				// Assembler Linear Interpolant X - Cells To Faces
  InterpZC2F::Assembler Rz_C2F_Assembler(grid.dim);				// Assembler Linear Interpolant X - Cells To Faces
		
  InterpXC2F Rx_C2F( Rx_C2F_Assembler );
  InterpYC2F Ry_C2F( Ry_C2F_Assembler );
  InterpZC2F Rz_C2F( Rz_C2F_Assembler );
	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Spatial Fields: Numerical Interpolations
  // -----------------------------------------------------------------------------------------------------------------------------	
  XSideField          Rfx( grid.dim, NULL, InternalStorage  );
  YSideField          Rfy( grid.dim, NULL, InternalStorage  );
  ZSideField          Rfz( grid.dim, NULL, InternalStorage  );
   	
  Rx_C2F.apply_to_field( funct.f, Rfx  );
  Ry_C2F.apply_to_field( funct.f, Rfy  );
  Rz_C2F.apply_to_field( funct.f, Rfz  );
  	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Check Equality
  // -----------------------------------------------------------------------------------------------------------------------------	
  std:: cout << " Test x: ";
  check_equality_C2F(grid, funct.f_faceX, Rfx, test_x, mean_rel_err[X_DIR], max_rel_err[X_DIR]);
  std:: cout << endl;
  			
  std:: cout << " Test y: ";
  check_equality_C2F(grid, funct.f_faceY, Rfy, test_y, mean_rel_err[Y_DIR], max_rel_err[Y_DIR]);
  std:: cout << endl;

  std:: cout << " Test z: ";
  check_equality_C2F(grid, funct.f_faceZ, Rfz, test_z, mean_rel_err[Z_DIR], max_rel_err[Z_DIR]);
  std:: cout << endl;
	
  ok = true;
  if (test_x==false || test_y==false || test_z==false)
    ok = false;
 
  return ok;
}

// *****************************************************************************************************************************	
// Test Linear Interpolant Operator F2C
// *****************************************************************************************************************************	
bool test_linear_interpolant_F2C(const grid_class& grid, analytical_class& funct, std::vector<double> &mean_rel_err, std::vector<double> &max_rel_err)
{
  using namespace FVStaggeredUniform;
  std::cout.setf(ios::scientific);

  bool test_x, test_y, test_z, ok;

  std::cout << "Test Function: " << funct.test.description() << endl;
	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Build the operators
  // -----------------------------------------------------------------------------------------------------------------------------	
  InterpXF2C::Assembler Rx_F2C_Assembler(grid.dim);				// Assembler Linear Interpolant X - Cells To Faces
  InterpYF2C::Assembler Ry_F2C_Assembler(grid.dim);				// Assembler Linear Interpolant X - Cells To Faces
  InterpZF2C::Assembler Rz_F2C_Assembler(grid.dim);				// Assembler Linear Interpolant X - Cells To Faces
	
  InterpXF2C Rx_F2C( Rx_F2C_Assembler );
  InterpYF2C Ry_F2C( Ry_F2C_Assembler );
  InterpZF2C Rz_F2C( Rz_F2C_Assembler );
 
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Spatial Fields: Numerical Interpolations
  // -----------------------------------------------------------------------------------------------------------------------------	
  CellField          Rfx( grid.dim, NULL, InternalStorage  );
  CellField          Rfy( grid.dim, NULL, InternalStorage  );
  CellField          Rfz( grid.dim, NULL, InternalStorage  );
    
  Rx_F2C.apply_to_field( funct.f_faceX, Rfx );
  Ry_F2C.apply_to_field( funct.f_faceY, Rfy );
  Rz_F2C.apply_to_field( funct.f_faceZ, Rfz );

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Check Equality
  // -----------------------------------------------------------------------------------------------------------------------------	
  std:: cout << " Test x: ";
  check_equality_F2C( grid, funct.f, Rfx, test_x, mean_rel_err[X_DIR], max_rel_err[X_DIR]);
  std:: cout << endl;
  			
  std:: cout << " Test y: ";
  check_equality_F2C(grid, funct.f, Rfy, test_y, mean_rel_err[Y_DIR], max_rel_err[Y_DIR]);
  std:: cout << endl;

  std:: cout << " Test z: ";
  check_equality_F2C(grid, funct.f, Rfz, test_z, mean_rel_err[Z_DIR], max_rel_err[Z_DIR]);
  std:: cout << endl;
	
  ok = true;
  if (test_x==false || test_y==false || test_z==false)
    ok = false;
 
  return ok;
}


// *****************************************************************************************************************************	
// Test Gradient C2F
// *****************************************************************************************************************************	
bool test_gradient_C2F(const grid_class& grid, analytical_class& funct, std::vector<double> &mean_rel_err, std::vector<double> &max_rel_err)
{
  using namespace FVStaggeredUniform;
  std::cout.setf(ios::scientific);

  std::cout << "Test Function: " << funct.test.description() << endl;
		
  bool test_x, test_y, test_z, ok;
	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Build the operators
  // -----------------------------------------------------------------------------------------------------------------------------	
  GradXC2F::Assembler xGrad_C2F_Assembler( grid.spacing, grid.dim );
  GradYC2F::Assembler yGrad_C2F_Assembler( grid.spacing, grid.dim );
  GradZC2F::Assembler zGrad_C2F_Assembler( grid.spacing, grid.dim );

  GradXC2F	xGrad_C2F( xGrad_C2F_Assembler);	
  GradYC2F	yGrad_C2F( yGrad_C2F_Assembler);	
  GradZC2F	zGrad_C2F( zGrad_C2F_Assembler);	

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Spatial Fields: Numerical Gradients
  // -----------------------------------------------------------------------------------------------------------------------------	
  XSideField      Gx( grid.dim, NULL,  InternalStorage  );
  YSideField      Gy( grid.dim, NULL,  InternalStorage  );
  ZSideField      Gz( grid.dim, NULL,  InternalStorage  );
 
  xGrad_C2F.apply_to_field( funct.f, Gx );	// Numerical Gradient X
  yGrad_C2F.apply_to_field( funct.f, Gy );	// Numerical Gradient Y
  zGrad_C2F.apply_to_field( funct.f, Gz );	// Numerical Gradient Z
  
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Check Equality
  // -----------------------------------------------------------------------------------------------------------------------------
  std:: cout << " Test x: ";
  check_equality_C2F(grid, funct.dfdx_faceX, Gx, test_x, mean_rel_err[X_DIR], max_rel_err[X_DIR]);
  std:: cout << endl;
  		
  std:: cout << " Test y: ";
  check_equality_C2F(grid, funct.dfdy_faceY, Gy, test_y, mean_rel_err[Y_DIR], max_rel_err[Y_DIR]);	
  std:: cout << endl;

  std:: cout << " Test z: ";
  check_equality_C2F(grid, funct.dfdz_faceZ, Gz, test_z, mean_rel_err[Z_DIR], max_rel_err[Z_DIR]);  	
  std:: cout << endl;
  	
  ok = true;
  if (test_x==false || test_y==false || test_z==false)
    ok = false;
 	
  return ok;
}

// *****************************************************************************************************************************	
// Test Gradient F2C
// *****************************************************************************************************************************	
bool test_gradient_F2C(const grid_class& grid, analytical_class& funct, std::vector<double> &mean_rel_err, std::vector<double> &max_rel_err)
{
  using namespace FVStaggeredUniform;
  std::cout.setf(ios::scientific);

  bool test_x, test_y, test_z, ok;

  std::cout << "Test Function: " << funct.test.description() << endl;	
		
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Build the operators
  // -----------------------------------------------------------------------------------------------------------------------------	
  GradXF2C::Assembler xGrad_F2C_Assembler( grid.spacing, grid.dim );
  GradYF2C::Assembler yGrad_F2C_Assembler( grid.spacing, grid.dim );
  GradZF2C::Assembler zGrad_F2C_Assembler( grid.spacing, grid.dim );

  GradXF2C	xGrad_F2C( xGrad_F2C_Assembler);	
  GradYF2C	yGrad_F2C( yGrad_F2C_Assembler);	
  GradZF2C	zGrad_F2C( zGrad_F2C_Assembler);	
	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Spatial Fields: Numerical Gradients
  // -----------------------------------------------------------------------------------------------------------------------------	
  CellField      Gx( grid.dim, NULL,  InternalStorage  );
  CellField      Gy( grid.dim, NULL,  InternalStorage  );
  CellField      Gz( grid.dim, NULL,  InternalStorage  );
 
  xGrad_F2C.apply_to_field( funct.f_faceX, Gx );	// Numerical Gradient X
  yGrad_F2C.apply_to_field( funct.f_faceY, Gy );	// Numerical Gradient Y
  zGrad_F2C.apply_to_field( funct.f_faceZ, Gz );	// Numerical Gradient Z
  
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Check Equality
  // -----------------------------------------------------------------------------------------------------------------------------
  std:: cout << " Test x: ";
  check_equality_F2C(grid, funct.dfdx, Gx, test_x, mean_rel_err[X_DIR], max_rel_err[X_DIR]);
  std:: cout << endl;
  		
  std:: cout << " Test y: ";
  check_equality_F2C( grid, funct.dfdy, Gy, test_y, mean_rel_err[Y_DIR], max_rel_err[Y_DIR]);	
  std:: cout << endl;

  std:: cout << " Test z: ";
  check_equality_F2C( grid, funct.dfdz, Gz, test_z, mean_rel_err[Z_DIR], max_rel_err[Z_DIR]);  	
  std:: cout << endl;
	
  ok = true;
  if (test_x==false || test_y==false || test_z==false)
    ok = false;
 	
  return ok;
}

// *****************************************************************************************************************************	
// Test Divergence C2F
// *****************************************************************************************************************************	
bool test_divergence_F2C(const grid_class& grid, analytical_class& funct, std::vector<double> &mean_rel_err, std::vector<double> &max_rel_err)
{
  using namespace FVStaggeredUniform;
  std::cout.setf(ios::scientific);

  bool test_x, test_y, test_z, ok;

  std::cout << "Test Function: " << funct.test.description() << endl;		
	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Build the operators
  // -----------------------------------------------------------------------------------------------------------------------------	
  DivXF2C::Assembler xDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
  DivYF2C::Assembler yDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
  DivZF2C::Assembler zDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
  
  GradXC2F::Assembler xGrad_C2F_Assembler( grid.spacing, grid.dim );
  GradYC2F::Assembler yGrad_C2F_Assembler( grid.spacing, grid.dim );
  GradZC2F::Assembler zGrad_C2F_Assembler( grid.spacing, grid.dim );
	
  DivXF2C    xDiv_F2C( xDiv_F2C_assembler );
  DivYF2C    yDiv_F2C( yDiv_F2C_assembler );
  DivZF2C    zDiv_F2C( zDiv_F2C_assembler );

  GradXC2F	xGrad_C2F( xGrad_C2F_Assembler);	
  GradYC2F	yGrad_C2F( yGrad_C2F_Assembler);	
  GradZC2F	zGrad_C2F( zGrad_C2F_Assembler);	
  	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Spatial Fields: Numerical Gradients
  // -----------------------------------------------------------------------------------------------------------------------------	
  XSideField   Gx( grid.dim, NULL,  InternalStorage  );
  YSideField   Gy( grid.dim, NULL,  InternalStorage  );
  ZSideField   Gz( grid.dim, NULL,  InternalStorage  );
 
  CellField    Dx( grid.dim, NULL, InternalStorage );	// Numerical  Divergence X
  CellField    Dy( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Y
  CellField    Dz( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Z

  xGrad_C2F.apply_to_field( funct.f, Gx );	// Numerical Gradient X
  yGrad_C2F.apply_to_field( funct.f, Gy );	// Numerical Gradient Y
  zGrad_C2F.apply_to_field( funct.f, Gz );	// Numerical Gradient Z
  		
  xDiv_F2C.apply_to_field( Gx, Dx );  // Numerical Divergence X
  yDiv_F2C.apply_to_field( Gy, Dy );  // Numerical Divergence Y
  zDiv_F2C.apply_to_field( Gz, Dz );  // Numerical Divergence Z
 
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Check Equality
  // -----------------------------------------------------------------------------------------------------------------------------
  std:: cout << " Test x: ";
  check_equality_F2C(grid, funct.d2fdx2, Dx, test_x, mean_rel_err[X_DIR], max_rel_err[X_DIR]);
  std:: cout << endl;
		
  std:: cout << " Test y: ";
  check_equality_F2C(grid, funct.d2fdy2, Dy, test_y, mean_rel_err[Y_DIR], max_rel_err[Y_DIR]);
  std:: cout << endl;

  std:: cout << " Test z: ";
  check_equality_F2C(grid, funct.d2fdz2, Dz, test_z, mean_rel_err[Z_DIR], max_rel_err[Z_DIR]);
  std:: cout << endl;
 	
  ok = true;
  if (test_x==false || test_y==false || test_z==false)
    ok = false;
 	
  return ok;
}

// *****************************************************************************************************************************	
// Test Divergence C2F
// *****************************************************************************************************************************	
bool test_divergence_C2F(const grid_class& grid, analytical_class& funct, std::vector<double> &mean_rel_err, std::vector<double> &max_rel_err)
{
  using namespace FVStaggeredUniform;
  std::cout.setf(ios::scientific);

  bool test_x, test_y, test_z, ok;

  std::cout << "Test Function: " << funct.test.description() << endl;
	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Build the operators
  // -----------------------------------------------------------------------------------------------------------------------------	
  DivXC2F::Assembler xDiv_C2F_assembler( grid.dim, grid.area, grid.volume );
  DivYC2F::Assembler yDiv_C2F_assembler( grid.dim, grid.area, grid.volume );
  DivZC2F::Assembler zDiv_C2F_assembler( grid.dim, grid.area, grid.volume );
  
  GradXF2C::Assembler xGrad_F2C_Assembler( grid.spacing, grid.dim );
  GradYF2C::Assembler yGrad_F2C_Assembler( grid.spacing, grid.dim );
  GradZF2C::Assembler zGrad_F2C_Assembler( grid.spacing, grid.dim );
	
  DivXC2F    xDiv_C2F( xDiv_C2F_assembler );
  DivYC2F    yDiv_C2F( yDiv_C2F_assembler );
  DivZC2F    zDiv_C2F( zDiv_C2F_assembler );

  GradXF2C	xGrad_F2C( xGrad_F2C_Assembler);	
  GradYF2C	yGrad_F2C( yGrad_F2C_Assembler);	
  GradZF2C	zGrad_F2C( zGrad_F2C_Assembler);	
  		
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Spatial Fields: Numerical Gradients
  // -----------------------------------------------------------------------------------------------------------------------------	
  CellField   Gx( grid.dim, NULL,  InternalStorage  );
  CellField   Gy( grid.dim, NULL,  InternalStorage  );
  CellField   Gz( grid.dim, NULL,  InternalStorage  );
 
  XSideField    Dx( grid.dim, NULL, InternalStorage );	// Numerical  Divergence X
  YSideField    Dy( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Y
  ZSideField    Dz( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Z

  xGrad_F2C.apply_to_field( funct.f_faceX, Gx );	// Numerical Gradient X
  yGrad_F2C.apply_to_field( funct.f_faceY, Gy );	// Numerical Gradient Y
  zGrad_F2C.apply_to_field( funct.f_faceZ, Gz );	// Numerical Gradient Z
   		
  xDiv_C2F.apply_to_field( Gx, Dx );  // Numerical Divergence X
  yDiv_C2F.apply_to_field( Gy, Dy );  // Numerical Divergence Y
  zDiv_C2F.apply_to_field( Gz, Dz );  // Numerical Divergence Z
 
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Check Equality
  // -----------------------------------------------------------------------------------------------------------------------------
  std:: cout << " Test x: ";
  check_equality_C2F(grid, funct.d2fdx2_faceX, Dx, test_x, mean_rel_err[X_DIR], max_rel_err[X_DIR]);
  std:: cout << endl;
		
  std:: cout << " Test y: ";
  check_equality_C2F( grid, funct.d2fdy2_faceY, Dy, test_y, mean_rel_err[Y_DIR], max_rel_err[Y_DIR]);
  std:: cout << endl;

  std:: cout << " Test z: ";
  check_equality_C2F( grid, funct.d2fdz2_faceZ, Dz, test_z, mean_rel_err[Z_DIR], max_rel_err[Z_DIR]);
  std:: cout << endl;
 	
  ok = true;
  if (test_x==false || test_y==false || test_z==false)
    ok = false;
 	
  return ok;
}


// *****************************************************************************************************************************	
// Test Scratch Operator - Laplacian
// *****************************************************************************************************************************	
bool test_scratch(const grid_class& grid, analytical_class& funct, std::vector<double> &mean_rel_err, std::vector<double> &max_rel_err)
{
  using namespace FVStaggeredUniform;
  std::cout.setf(ios::scientific);

  bool test_x, test_y, test_z, ok;

  std::cout << "Test Function: " << funct.test.description() << endl;

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Build the operators
  // -----------------------------------------------------------------------------------------------------------------------------	
  DivXF2C::Assembler xDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
  DivYF2C::Assembler yDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
  DivZF2C::Assembler zDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
  
  GradXC2F::Assembler xGrad_C2F_Assembler( grid.spacing, grid.dim );
  GradYC2F::Assembler yGrad_C2F_Assembler( grid.spacing, grid.dim );
  GradZC2F::Assembler zGrad_C2F_Assembler( grid.spacing, grid.dim );
	
  SxCell::Assembler  xLaplacian_Assembler(grid.dim);
  SyCell::Assembler  yLaplacian_Assembler(grid.dim);
  SzCell::Assembler  zLaplacian_Assembler(grid.dim);
	
  DivXF2C    xDiv_F2C( xDiv_F2C_assembler );
  DivYF2C    yDiv_F2C( yDiv_F2C_assembler );
  DivZF2C    zDiv_F2C( zDiv_F2C_assembler );

  GradXC2F	xGrad_C2F( xGrad_C2F_Assembler);	
  GradYC2F	yGrad_C2F( yGrad_C2F_Assembler);	
  GradZC2F	zGrad_C2F( zGrad_C2F_Assembler);	
	
  SxCell     xLaplacian(xLaplacian_Assembler);
  SyCell     yLaplacian(yLaplacian_Assembler);	
  SzCell     zLaplacian(zLaplacian_Assembler);	
	
  xDiv_F2C.apply_to_op(xGrad_C2F, xLaplacian);
  yDiv_F2C.apply_to_op(yGrad_C2F, yLaplacian);
  zDiv_F2C.apply_to_op(zGrad_C2F, zLaplacian);

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Laplacian Application
  // -----------------------------------------------------------------------------------------------------------------------------		
  CellField    Dx( grid.dim, NULL, InternalStorage );	// Numerical  Divergence X
  CellField    Dy( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Y
  CellField    Dz( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Z
  		
  xLaplacian.apply_to_field( funct.f, Dx );  // Numerical Divergence X
  yLaplacian.apply_to_field( funct.f, Dy );  // Numerical Divergence Y
  zLaplacian.apply_to_field( funct.f, Dz );  // Numerical Divergence Z		  		

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Check Equality - Laplacian
  // -----------------------------------------------------------------------------------------------------------------------------
  std:: cout << " Test x: ";
  check_equality_F2C(grid, funct.d2fdx2, Dx, test_x, mean_rel_err[X_DIR], max_rel_err[X_DIR]);
  std:: cout << endl;
		
  std:: cout << " Test y: ";
  check_equality_F2C(grid, funct.d2fdy2, Dy, test_y, mean_rel_err[Y_DIR], max_rel_err[Y_DIR]);
  std:: cout << endl;

  std:: cout << " Test z: ";
  check_equality_F2C(grid, funct.d2fdz2, Dz, test_z, mean_rel_err[Z_DIR], max_rel_err[Z_DIR]);
  std:: cout << endl;
	  	  	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Comparison with Div(Grad(F))
  // -----------------------------------------------------------------------------------------------------------------------------	
  {	
    XSideField   _Gx( grid.dim, NULL, InternalStorage );
    YSideField   _Gy( grid.dim, NULL, InternalStorage );
    ZSideField   _Gz( grid.dim, NULL, InternalStorage );
 	
    CellField    _Dx( grid.dim, NULL, InternalStorage );	// Numerical  Divergence X
    CellField    _Dy( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Y
    CellField    _Dz( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Z
  			
    xGrad_C2F.apply_to_field( funct.f, _Gx );	// Numerical Gradient X
    yGrad_C2F.apply_to_field( funct.f, _Gy );	// Numerical Gradient Y
    zGrad_C2F.apply_to_field( funct.f, _Gz );	// Numerical Gradient Z
  		
    xDiv_F2C.apply_to_field( _Gx, _Dx );  // Numerical Divergence X
    yDiv_F2C.apply_to_field( _Gy, _Dy );  // Numerical Divergence Y
    zDiv_F2C.apply_to_field( _Gz, _Dz );  // Numerical Divergence Z
	
    _Dx-=Dx;
    _Dy-=Dy;
    _Dz-=Dz;
 
    bool comparison_x = true;
    bool comparison_y = true;
    bool comparison_z = true;
  		
    for( int i=0; i<_Dx.get_ntotal(); ++i )
      {
	if (_Dx.begin()[i]<=1.e-16 && Dx.begin()[i]<=1.e-16) comparison_x = true;
	else if ( _Dx.begin()[i]/(Dx.begin()[i]+TOL) >= 1.e-10) comparison_x = false;
      }
    for( int i=0; i<_Dy.get_ntotal(); ++i )
      {	
	if (_Dy.begin()[i]<=1.e-16 && Dy.begin()[i]<=1.e-16) comparison_y = true;
	else if ( _Dy.begin()[i]/(Dy.begin()[i]+TOL) >= 1.e-10) comparison_y = false;
      }
    for( int i=0; i<_Dz.get_ntotal(); ++i )
      {
	if (_Dz.begin()[i]<=1.e-16 && Dz.begin()[i]<=1.e-16) comparison_z = true;
	else if ( _Dz.begin()[i]/(Dz.begin()[i]+TOL) >= 1.e-10) comparison_z = false;
      }
		
    if (comparison_x == true && comparison_y == true && comparison_z == true)
      std::cout << "Comparison with Div(Grad(F)): PASS" << endl;
    else
      std::cout << "Comparison with Div(Grad(F)): FAIL " << comparison_x << " " << comparison_y << " " << comparison_z << " " << endl;
  }

  ok = true;
  if (test_x==false || test_y==false || test_z==false)
    ok = false;
 	
  return ok;
}

// *****************************************************************************************************************************	
// Mixed Derivatives
// *****************************************************************************************************************************	
bool test_mixed_derivatives(const grid_class& grid, analytical_class& funct, std::vector<double> &mean_rel_err, std::vector<double> &max_rel_err)
{
  using namespace FVStaggeredUniform;
  std::cout.setf(ios::scientific);
	
  bool ok;
  bool test_xy, test_xz;
  bool test_yx, test_yz;
  bool test_zx, test_zy;

  double mean_rel_err_xy, mean_rel_err_xz;
  double mean_rel_err_yx, mean_rel_err_yz;
  double mean_rel_err_zx, mean_rel_err_zy;

  double max_rel_err_xy, max_rel_err_xz;
  double max_rel_err_yx, max_rel_err_yz;
  double max_rel_err_zx, max_rel_err_zy;
	
  std::cout << "Test Function: " << funct.test.description() << endl;

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Build the operators
  // -----------------------------------------------------------------------------------------------------------------------------	
  InterpXC2F::Assembler Rx_C2F_Assembler(grid.dim);				// Assembler Linear Interpolant X - Cells To Faces
  InterpYC2F::Assembler Ry_C2F_Assembler(grid.dim);				// Assembler Linear Interpolant Y - Cells To Faces
  InterpZC2F::Assembler Rz_C2F_Assembler(grid.dim);				// Assembler Linear Interpolant Z - Cells To Faces

  InterpXF2C::Assembler Rx_F2C_Assembler(grid.dim);				// Assembler Linear Interpolant X - Faces To Cells
  InterpYF2C::Assembler Ry_F2C_Assembler(grid.dim);				// Assembler Linear Interpolant Y - Faces To Cells
  InterpZF2C::Assembler Rz_F2C_Assembler(grid.dim);				// Assembler Linear Interpolant Z - Faces To Cells
	
  DivXF2C::Assembler xDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
  DivYF2C::Assembler yDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
  DivZF2C::Assembler zDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
  
  GradXC2F::Assembler xGrad_C2F_Assembler( grid.spacing, grid.dim );
  GradYC2F::Assembler yGrad_C2F_Assembler( grid.spacing, grid.dim );
  GradZC2F::Assembler zGrad_C2F_Assembler( grid.spacing, grid.dim );
				
  InterpXC2F Rx_C2F = InterpXC2F( Rx_C2F_Assembler );
  InterpYC2F Ry_C2F = InterpYC2F( Ry_C2F_Assembler );
  InterpZC2F Rz_C2F = InterpZC2F( Rz_C2F_Assembler );

  InterpXF2C Rx_F2C = InterpXF2C( Rx_F2C_Assembler );
  InterpYF2C Ry_F2C = InterpYF2C( Ry_F2C_Assembler );
  InterpZF2C Rz_F2C = InterpZF2C( Rz_F2C_Assembler );
		
  DivXF2C    xDiv_F2C( xDiv_F2C_assembler );
  DivYF2C    yDiv_F2C( yDiv_F2C_assembler );
  DivZF2C    zDiv_F2C( zDiv_F2C_assembler );

  GradXC2F	xGrad_C2F( xGrad_C2F_Assembler);	
  GradYC2F	yGrad_C2F( yGrad_C2F_Assembler);	
  GradZC2F	zGrad_C2F( zGrad_C2F_Assembler);	
		 	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Scratch Operator Application
  // -----------------------------------------------------------------------------------------------------------------------------	
	
  XSideField      FaceX( grid.dim, NULL,  InternalStorage  );
  YSideField      FaceY( grid.dim, NULL,  InternalStorage  );
  ZSideField      FaceZ( grid.dim, NULL,  InternalStorage  );

  CellField    D2DxDy( grid.dim, NULL, InternalStorage );	// Numerical  Divergence X
  CellField    D2DxDz( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Y
  CellField    D2DyDx( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Z
  CellField    D2DyDz( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Z
  CellField    D2DzDx( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Z
  CellField    D2DzDy( grid.dim, NULL, InternalStorage );	// Numerical  Divergence Z
	 
  // xy Derivatives
  yGrad_C2F.apply_to_field(funct.f,FaceY);	
  Ry_F2C.apply_to_field(FaceY,D2DxDy);
  Rx_C2F.apply_to_field(D2DxDy, FaceX);
  xDiv_F2C.apply_to_field(FaceX,D2DxDy);	

  // yx Derivatives	
  xGrad_C2F.apply_to_field(funct.f,FaceX);	
  Rx_F2C.apply_to_field(FaceX,D2DyDx);
  Ry_C2F.apply_to_field(D2DyDx, FaceY);
  yDiv_F2C.apply_to_field(FaceY,D2DyDx);	
		
  // xz Derivatives
  zGrad_C2F.apply_to_field(funct.f,FaceZ);	
  Rz_F2C.apply_to_field(FaceZ,D2DxDz);
  Rx_C2F.apply_to_field(D2DxDz, FaceX);
  xDiv_F2C.apply_to_field(FaceX,D2DxDz);	

  // zx Derivatives
  xGrad_C2F.apply_to_field(funct.f,FaceX);	
  Rx_F2C.apply_to_field(FaceX,D2DzDx);
  Rz_C2F.apply_to_field(D2DzDx, FaceZ);
  zDiv_F2C.apply_to_field(FaceZ,D2DzDx);	

  // xz Derivatives
  zGrad_C2F.apply_to_field(funct.f,FaceZ);	
  Rz_F2C.apply_to_field(FaceZ,D2DyDz);
  Ry_C2F.apply_to_field(D2DyDz, FaceY);
  yDiv_F2C.apply_to_field(FaceY,D2DyDz);	

  // zx Derivatives
  yGrad_C2F.apply_to_field(funct.f,FaceY);	
  Ry_F2C.apply_to_field(FaceY,D2DzDy);
  Rz_C2F.apply_to_field(D2DzDy, FaceZ);
  zDiv_F2C.apply_to_field(FaceZ,D2DzDy);	
				
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Check Equality - Scratch Operator
  // -----------------------------------------------------------------------------------------------------------------------------
  std:: cout << " Test x-y: ";
  check_equality_F2C(grid, funct.d2fdxdy, D2DxDy, test_xy, mean_rel_err_xy, max_rel_err_xy);
  std:: cout << endl;

  std:: cout << " Test x-z: ";
  check_equality_F2C(grid, funct.d2fdxdz, D2DxDz, test_xz, mean_rel_err_xz, max_rel_err_xz);
  std:: cout << endl;

  std:: cout << " Test y-x: ";
  check_equality_F2C(grid, funct.d2fdydx, D2DyDx, test_yx, mean_rel_err_yx, max_rel_err_yx);
  std:: cout << endl;

  std:: cout << " Test y-z: ";
  check_equality_F2C(grid, funct.d2fdydz, D2DyDz, test_yz, mean_rel_err_yz, max_rel_err_yz);
  std:: cout << endl;
  	
  std:: cout << " Test z-x: ";
  check_equality_F2C(grid, funct.d2fdzdx, D2DzDx, test_zx, mean_rel_err_zx, max_rel_err_zx);
  std:: cout << endl;

  std:: cout << " Test z-y: ";
  check_equality_F2C(grid, funct.d2fdzdy, D2DzDy, test_zy, mean_rel_err_zy, max_rel_err_zy);
  std:: cout << endl;
	 
  mean_rel_err[X_DIR] = mean_rel_err_xy;
  mean_rel_err[Y_DIR] = mean_rel_err_yz;
  mean_rel_err[Z_DIR] = mean_rel_err_zx;
	
  max_rel_err[X_DIR] = max_rel_err_xy;
  max_rel_err[Y_DIR] = max_rel_err_yz;
  max_rel_err[Z_DIR] = max_rel_err_zx;

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Symmetry Test
  // -----------------------------------------------------------------------------------------------------------------------------	
  {	  		
    bool comparison_xy = true;
    bool comparison_xz = true;
    bool comparison_yz = true;
  		
    for(unsigned int i=0; i<grid.index_bothsides_cell.size(); i++)
      {    
	int point = grid.index_bothsides_cell[i];
	if ( fabs( D2DxDy.begin()[point] - D2DyDx.begin()[point]) / (fabs(D2DxDy.begin()[point])+TOL) >= 1.e-4 )
	  comparison_xy = false; 		
      }	
		
    for(unsigned int i=0; i<grid.index_bothsides_cell.size(); i++)
      {    
	int point = grid.index_bothsides_cell[i];
	if ( fabs( D2DxDz.begin()[point] - D2DzDx.begin()[point]) / (fabs(D2DxDz.begin()[point])+TOL) >= 1.e-4 )
	  comparison_xz = false; 		
      }	
		
    for(unsigned int i=0; i<grid.index_bothsides_cell.size(); i++)
      {    
	int point = grid.index_bothsides_cell[i];
	if ( fabs( D2DyDz.begin()[point] - D2DzDy.begin()[point]) / (fabs(D2DyDz.begin()[point])+TOL) >= 1.e-4 )
	  comparison_yz = false; 		
      }	
		
    if (comparison_xy == true && comparison_xz == true && comparison_yz == true)
      std::cout << "Symmetry test: PASS" << endl;
    else
      std::cout << "Symmetry test: FAIL " << comparison_xy << " " << comparison_xz << " " << comparison_yz << " " << endl;
  }
	
  	
  ok = true;
  if (test_yz==false || test_xy==false || test_xz==false)
    ok = false;
 	
  return ok;
}



// *****************************************************************************************************************************	
// Test Scaling C2F
// *****************************************************************************************************************************	
void test_scale_C2F(string kind, std::vector<int> dim, std::vector<bool> &ok)
{
  using namespace FVStaggeredUniform;
  std::cout.setf(ios::scientific);
	
  bool ok_linear_interpolant = false;
  bool ok_gradient = false;
  bool ok_divergence = false;
		
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Geometry Construction
  // -----------------------------------------------------------------------------------------------------------------------------	
  grid_class grid(dim);

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Spatial Fields
  // -----------------------------------------------------------------------------------------------------------------------------	
  vector<double> f1(grid.ntot_cell);
  vector<double> f2(grid.ntot_cell);

  for( int ix=0; ix<grid.ntot_cell; ++ix )
    {
      f1[ix] = double(ix);							 
      f2[ix] = double(ix)*double(2*ix) + sqrt(double(ix)); 
    }
		
  XSideField      field1X( grid.dim, &f1[0],  	ExternalStorage );
  XSideField      field2X( grid.dim, &f2[0],  	ExternalStorage );
  YSideField      field1Y( grid.dim, &f1[0],  	ExternalStorage );
  YSideField      field2Y( grid.dim, &f2[0],  	ExternalStorage );
  ZSideField      field1Z( grid.dim, &f1[0],  	ExternalStorage );
  ZSideField      field2Z( grid.dim, &f2[0],  	ExternalStorage );
  CellField    	 field1( grid.dim, &f1[0],  	ExternalStorage );
  CellField      	 field2( grid.dim, &f2[0],  	ExternalStorage );

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Linear Interpolant Operator
  // -----------------------------------------------------------------------------------------------------------------------------	
  {
    InterpXC2F::Assembler Rx_C2F_Assembler(grid.dim);		// Assembler Linear Interpolant X - Cells To Faces
    InterpYC2F::Assembler Ry_C2F_Assembler(grid.dim);		// Assembler Linear Interpolant Y - Cells To Faces
    InterpZC2F::Assembler Rz_C2F_Assembler(grid.dim);			// Assembler Linear Interpolant Z - Cells To Faces
				
    InterpXC2F RxOld_C2F = InterpXC2F( Rx_C2F_Assembler );
    InterpYC2F RyOld_C2F = InterpYC2F( Ry_C2F_Assembler );
    InterpZC2F RzOld_C2F = InterpZC2F( Rz_C2F_Assembler );

    InterpXC2F Rx1_C2F = InterpXC2F( Rx_C2F_Assembler );
    InterpYC2F Ry1_C2F = InterpYC2F( Ry_C2F_Assembler );
    InterpZC2F Rz1_C2F = InterpZC2F( Rz_C2F_Assembler );

    InterpXC2F Rx2_C2F = InterpXC2F( Rx_C2F_Assembler );
    InterpYC2F Ry2_C2F = InterpYC2F( Ry_C2F_Assembler );
    InterpZC2F Rz2_C2F = InterpZC2F( Rz_C2F_Assembler );
	
    if (kind == "LEFT")
      {	
	Rx1_C2F.left_scale( field1X );
	Ry1_C2F.left_scale( field1Y);
	Rz1_C2F.left_scale( field1Z);
  	
	Rx2_C2F.left_scale( field2X );
	Ry2_C2F.left_scale( field2Y);
	Rz2_C2F.left_scale( field2Z );
      }
    else if (kind == "RIGHT")
      {	
	Rx1_C2F.right_scale( field1 );
	Ry1_C2F.right_scale( field1 );
	Rz1_C2F.right_scale( field1 );
  	
	Rx2_C2F.right_scale( field2 );
	Ry2_C2F.right_scale( field2 );
	Rz2_C2F.right_scale( field2 );
      }
  	
    ok_linear_interpolant = compare_scaling(kind, Rx1_C2F.get_linalg_mat(), RxOld_C2F.get_linalg_mat(), field1.begin());
    ok_linear_interpolant = compare_scaling(kind, Ry1_C2F.get_linalg_mat(), RyOld_C2F.get_linalg_mat(), field1.begin());
    ok_linear_interpolant = compare_scaling(kind, Rz1_C2F.get_linalg_mat(), RzOld_C2F.get_linalg_mat(), field1.begin());

    ok_linear_interpolant = compare_scaling(kind, Rx2_C2F.get_linalg_mat(), RxOld_C2F.get_linalg_mat(), field2.begin());
    ok_linear_interpolant = compare_scaling(kind, Ry2_C2F.get_linalg_mat(), RyOld_C2F.get_linalg_mat(), field2.begin());
    ok_linear_interpolant = compare_scaling(kind, Rz2_C2F.get_linalg_mat(), RzOld_C2F.get_linalg_mat(), field2.begin());
  }
	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Gradient Operator
  // -----------------------------------------------------------------------------------------------------------------------------	
  {
    GradXC2F::Assembler xGrad_C2F_Assembler( grid.spacing, grid.dim );
    GradYC2F::Assembler yGrad_C2F_Assembler( grid.spacing, grid.dim );
    GradZC2F::Assembler  zGrad_C2F_Assembler( grid.spacing, grid.dim );
				
    GradXC2F xGradOld_C2F = GradXC2F( xGrad_C2F_Assembler);	
    GradYC2F yGradOld_C2F = GradYC2F( yGrad_C2F_Assembler);	
    GradZC2F zGradOld_C2F = GradZC2F( zGrad_C2F_Assembler);	
	
    GradXC2F xGrad1_C2F = GradXC2F( xGrad_C2F_Assembler);	
    GradYC2F yGrad1_C2F = GradYC2F( yGrad_C2F_Assembler);	
    GradZC2F zGrad1_C2F = GradZC2F( zGrad_C2F_Assembler);	
	
    GradXC2F xGrad2_C2F = GradXC2F( xGrad_C2F_Assembler);	
    GradYC2F yGrad2_C2F = GradYC2F( yGrad_C2F_Assembler);	
    GradZC2F zGrad2_C2F = GradZC2F( zGrad_C2F_Assembler);	
	
    if (kind == "LEFT")
      {
	xGrad1_C2F.left_scale( field1X );
	yGrad1_C2F.left_scale( field1Y);
	zGrad1_C2F.left_scale( field1Z );
	
	xGrad2_C2F.left_scale( field2X );
	yGrad2_C2F.left_scale( field2Y );
	zGrad2_C2F.left_scale( field2Z );
      }
    else if (kind == "RIGHT")
      {	
	xGrad1_C2F.right_scale( field1 );
	yGrad1_C2F.right_scale( field1 );
	zGrad1_C2F.right_scale( field1 );
  	
	xGrad2_C2F.right_scale( field2 );
	yGrad2_C2F.right_scale( field2 );
	zGrad2_C2F.right_scale( field2 );
      }

    ok_gradient = compare_scaling(kind, xGrad1_C2F.get_linalg_mat(), xGradOld_C2F.get_linalg_mat(), field1.begin());		
    ok_gradient = compare_scaling(kind, yGrad1_C2F.get_linalg_mat(), yGradOld_C2F.get_linalg_mat(), field1.begin());		
    ok_gradient = compare_scaling(kind, zGrad1_C2F.get_linalg_mat(), zGradOld_C2F.get_linalg_mat(), field1.begin());

    ok_gradient = compare_scaling(kind, xGrad2_C2F.get_linalg_mat(), xGradOld_C2F.get_linalg_mat(), field2.begin());	
    ok_gradient = compare_scaling(kind, yGrad2_C2F.get_linalg_mat(), yGradOld_C2F.get_linalg_mat(), field2.begin());
    ok_gradient = compare_scaling(kind, yGrad2_C2F.get_linalg_mat(), yGradOld_C2F.get_linalg_mat(), field2.begin());
  }
  	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Divergence Operator
  // -----------------------------------------------------------------------------------------------------------------------------	
  {	
    DivXC2F::Assembler xDiv_C2F_assembler( grid.dim, grid.area, grid.volume );
    DivYC2F::Assembler yDiv_C2F_assembler( grid.dim, grid.area, grid.volume );
    DivZC2F::Assembler zDiv_C2F_assembler( grid.dim, grid.area, grid.volume );
  				
    DivXC2F xDivOld_C2F = DivXC2F( xDiv_C2F_assembler );
    DivYC2F yDivOld_C2F = DivYC2F( yDiv_C2F_assembler );
    DivZC2F zDivOld_C2F = DivZC2F( zDiv_C2F_assembler );

    DivXC2F xDiv1_C2F = DivXC2F( xDiv_C2F_assembler );
    DivYC2F yDiv1_C2F = DivYC2F( yDiv_C2F_assembler );
    DivZC2F zDiv1_C2F = DivZC2F( zDiv_C2F_assembler );

    DivXC2F xDiv2_C2F = DivXC2F( xDiv_C2F_assembler );
    DivYC2F yDiv2_C2F = DivYC2F( yDiv_C2F_assembler );
    DivZC2F zDiv2_C2F = DivZC2F( zDiv_C2F_assembler );
				
    if (kind == "LEFT")
      {
	xDiv1_C2F.left_scale( field1X );
	yDiv1_C2F.left_scale( field1Y );
	zDiv1_C2F.left_scale( field1Z );

	xDiv2_C2F.left_scale( field2X );
	yDiv2_C2F.left_scale( field2Y );
	zDiv2_C2F.left_scale( field2Z );
      }
    else if (kind == "RIGHT")
      {	
	xDiv1_C2F.right_scale( field1 );
	yDiv1_C2F.right_scale( field1 );
	zDiv1_C2F.right_scale( field1 );
  	
	xDiv2_C2F.right_scale( field2 );
	yDiv2_C2F.right_scale( field2 );
	zDiv2_C2F.right_scale( field2 );
      }

    ok_divergence = compare_scaling(kind, xDiv1_C2F.get_linalg_mat(), xDivOld_C2F.get_linalg_mat(), field1.begin());		
    ok_divergence = compare_scaling(kind, yDiv1_C2F.get_linalg_mat(), yDivOld_C2F.get_linalg_mat(), field1.begin());		
    ok_divergence = compare_scaling(kind, zDiv1_C2F.get_linalg_mat(), zDivOld_C2F.get_linalg_mat(), field1.begin());		
	    
    ok_divergence = compare_scaling(kind, xDiv2_C2F.get_linalg_mat(), xDivOld_C2F.get_linalg_mat(), field2.begin());		
    ok_divergence = compare_scaling(kind, yDiv2_C2F.get_linalg_mat(), yDivOld_C2F.get_linalg_mat(), field2.begin());		
    ok_divergence = compare_scaling(kind, zDiv2_C2F.get_linalg_mat(), zDivOld_C2F.get_linalg_mat(), field2.begin());		
  }
	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Final Results
  // -----------------------------------------------------------------------------------------------------------------------------	
  ok[0] = ok_linear_interpolant;      
  ok[1] = ok_gradient;
  ok[2] = ok_divergence;
}

// *****************************************************************************************************************************	
// Test Scaling F2C
// *****************************************************************************************************************************	
void test_scale_F2C(string kind, std::vector<int> dim, std::vector<bool> &ok)
{
  using namespace FVStaggeredUniform;
  std::cout.setf(ios::scientific);
	
  bool ok_linear_interpolant = false;
  bool ok_gradient = false;
  bool ok_divergence = false;
		
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Geometry Construction
  // -----------------------------------------------------------------------------------------------------------------------------	
  grid_class grid(dim);

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Spatial Fields
  // -----------------------------------------------------------------------------------------------------------------------------	
  vector<double> f1(grid.ntot_cell);
  vector<double> f2(grid.ntot_cell);

  for( int ix=0; ix<grid.ntot_cell; ++ix )
    {
      f1[ix] = double(ix);							 
      f2[ix] = double(ix)*double(2*ix) + sqrt(double(ix)); 
    }
		
  XSideField      field1X( grid.dim, &f1[0],  	ExternalStorage );
  XSideField      field2X( grid.dim, &f2[0],  	ExternalStorage );
  YSideField      field1Y( grid.dim, &f1[0],  	ExternalStorage );
  YSideField      field2Y( grid.dim, &f2[0],  	ExternalStorage );
  ZSideField      field1Z( grid.dim, &f1[0],  	ExternalStorage );
  ZSideField      field2Z( grid.dim, &f2[0],  	ExternalStorage );
  CellField    	 field1( grid.dim, &f1[0],  	ExternalStorage );
  CellField      	 field2( grid.dim, &f2[0],  	ExternalStorage );

  // -----------------------------------------------------------------------------------------------------------------------------	
  // Linear Interpolant Operator
  // -----------------------------------------------------------------------------------------------------------------------------	
  {
    InterpXF2C::Assembler Rx_F2C_Assembler(grid.dim);		// Assembler Linear Interpolant X - Cells To Faces
    InterpYF2C::Assembler Ry_F2C_Assembler(grid.dim);		// Assembler Linear Interpolant Y - Cells To Faces
    InterpZF2C::Assembler Rz_F2C_Assembler(grid.dim);			// Assembler Linear Interpolant Z - Cells To Faces
				
    InterpXF2C RxOld_F2C = InterpXF2C( Rx_F2C_Assembler );
    InterpYF2C RyOld_F2C = InterpYF2C( Ry_F2C_Assembler );
    InterpZF2C RzOld_F2C = InterpZF2C( Rz_F2C_Assembler );

    InterpXF2C Rx1_F2C = InterpXF2C( Rx_F2C_Assembler );
    InterpYF2C Ry1_F2C = InterpYF2C( Ry_F2C_Assembler );
    InterpZF2C Rz1_F2C = InterpZF2C( Rz_F2C_Assembler );

    InterpXF2C Rx2_F2C = InterpXF2C( Rx_F2C_Assembler );
    InterpYF2C Ry2_F2C = InterpYF2C( Ry_F2C_Assembler );
    InterpZF2C Rz2_F2C = InterpZF2C( Rz_F2C_Assembler );
	
    if (kind == "LEFT")
      {	
	Rx1_F2C.left_scale( field1 );
	Ry1_F2C.left_scale( field1);
	Rz1_F2C.left_scale( field1);
  	
	Rx2_F2C.left_scale( field2 );
	Ry2_F2C.left_scale( field2);
	Rz2_F2C.left_scale( field2 );
      }
    else if (kind == "RIGHT")
      {	
	Rx1_F2C.right_scale( field1X );
	Ry1_F2C.right_scale( field1Y );
	Rz1_F2C.right_scale( field1Z );
  	
	Rx2_F2C.right_scale( field2X );
	Ry2_F2C.right_scale( field2Y );
	Rz2_F2C.right_scale( field2Z );
      }
  	
    ok_linear_interpolant = compare_scaling(kind, Rx1_F2C.get_linalg_mat(), RxOld_F2C.get_linalg_mat(), field1.begin());
    ok_linear_interpolant = compare_scaling(kind, Ry1_F2C.get_linalg_mat(), RyOld_F2C.get_linalg_mat(), field1.begin());
    ok_linear_interpolant = compare_scaling(kind, Rz1_F2C.get_linalg_mat(), RzOld_F2C.get_linalg_mat(), field1.begin());

    ok_linear_interpolant = compare_scaling(kind, Rx2_F2C.get_linalg_mat(), RxOld_F2C.get_linalg_mat(), field2.begin());
    ok_linear_interpolant = compare_scaling(kind, Ry2_F2C.get_linalg_mat(), RyOld_F2C.get_linalg_mat(), field2.begin());
    ok_linear_interpolant = compare_scaling(kind, Rz2_F2C.get_linalg_mat(), RzOld_F2C.get_linalg_mat(), field2.begin());
  }
	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Gradient Operator
  // -----------------------------------------------------------------------------------------------------------------------------	
  {
    GradXF2C::Assembler xGrad_F2C_Assembler( grid.spacing, grid.dim );
    GradYF2C::Assembler yGrad_F2C_Assembler( grid.spacing, grid.dim );
    GradZF2C::Assembler  zGrad_F2C_Assembler( grid.spacing, grid.dim );
				
    GradXF2C xGradOld_F2C = GradXF2C( xGrad_F2C_Assembler);	
    GradYF2C yGradOld_F2C = GradYF2C( yGrad_F2C_Assembler);	
    GradZF2C zGradOld_F2C = GradZF2C( zGrad_F2C_Assembler);	
	
    GradXF2C xGrad1_F2C = GradXF2C( xGrad_F2C_Assembler);	
    GradYF2C yGrad1_F2C = GradYF2C( yGrad_F2C_Assembler);	
    GradZF2C zGrad1_F2C = GradZF2C( zGrad_F2C_Assembler);	
	
    GradXF2C xGrad2_F2C = GradXF2C( xGrad_F2C_Assembler);	
    GradYF2C yGrad2_F2C = GradYF2C( yGrad_F2C_Assembler);	
    GradZF2C zGrad2_F2C = GradZF2C( zGrad_F2C_Assembler);	
	
    if (kind == "LEFT")
      {
	xGrad1_F2C.left_scale( field1 );
	yGrad1_F2C.left_scale( field1);
	zGrad1_F2C.left_scale( field1 );
	
	xGrad2_F2C.left_scale( field2 );
	yGrad2_F2C.left_scale( field2 );
	zGrad2_F2C.left_scale( field2 );
      }
    else if (kind == "RIGHT")
      {	
	xGrad1_F2C.right_scale( field1X );
	yGrad1_F2C.right_scale( field1Y );
	zGrad1_F2C.right_scale( field1Z );
  	
	xGrad2_F2C.right_scale( field2X );
	yGrad2_F2C.right_scale( field2Y );
	zGrad2_F2C.right_scale( field2Z );
      }

    ok_gradient = compare_scaling(kind, xGrad1_F2C.get_linalg_mat(), xGradOld_F2C.get_linalg_mat(), field1.begin());		
    ok_gradient = compare_scaling(kind, yGrad1_F2C.get_linalg_mat(), yGradOld_F2C.get_linalg_mat(), field1.begin());		
    ok_gradient = compare_scaling(kind, zGrad1_F2C.get_linalg_mat(), zGradOld_F2C.get_linalg_mat(), field1.begin());

    ok_gradient = compare_scaling(kind, xGrad2_F2C.get_linalg_mat(), xGradOld_F2C.get_linalg_mat(), field2.begin());	
    ok_gradient = compare_scaling(kind, yGrad2_F2C.get_linalg_mat(), yGradOld_F2C.get_linalg_mat(), field2.begin());
    ok_gradient = compare_scaling(kind, yGrad2_F2C.get_linalg_mat(), yGradOld_F2C.get_linalg_mat(), field2.begin());
  }
  	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Divergence Operator
  // -----------------------------------------------------------------------------------------------------------------------------	
  {	
    DivXF2C::Assembler xDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
    DivYF2C::Assembler yDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
    DivZF2C::Assembler zDiv_F2C_assembler( grid.dim, grid.area, grid.volume );
  				
    DivXF2C xDivOld_F2C = DivXF2C( xDiv_F2C_assembler );
    DivYF2C yDivOld_F2C = DivYF2C( yDiv_F2C_assembler );
    DivZF2C zDivOld_F2C = DivZF2C( zDiv_F2C_assembler );

    DivXF2C xDiv1_F2C = DivXF2C( xDiv_F2C_assembler );
    DivYF2C yDiv1_F2C = DivYF2C( yDiv_F2C_assembler );
    DivZF2C zDiv1_F2C = DivZF2C( zDiv_F2C_assembler );

    DivXF2C xDiv2_F2C = DivXF2C( xDiv_F2C_assembler );
    DivYF2C yDiv2_F2C = DivYF2C( yDiv_F2C_assembler );
    DivZF2C zDiv2_F2C = DivZF2C( zDiv_F2C_assembler );
				
    if (kind == "LEFT")
      {
	xDiv1_F2C.left_scale( field1 );
	yDiv1_F2C.left_scale( field1 );
	zDiv1_F2C.left_scale( field1 );

	xDiv2_F2C.left_scale( field2 );
	yDiv2_F2C.left_scale( field2 );
	zDiv2_F2C.left_scale( field2 );
      }
    else if (kind == "RIGHT")
      {	
	xDiv1_F2C.right_scale( field1X );
	yDiv1_F2C.right_scale( field1Y );
	zDiv1_F2C.right_scale( field1Z );
  	
	xDiv2_F2C.right_scale( field2X );
	yDiv2_F2C.right_scale( field2Y );
	zDiv2_F2C.right_scale( field2Z );
      }

    ok_divergence = compare_scaling(kind, xDiv1_F2C.get_linalg_mat(), xDivOld_F2C.get_linalg_mat(), field1.begin());		
    ok_divergence = compare_scaling(kind, yDiv1_F2C.get_linalg_mat(), yDivOld_F2C.get_linalg_mat(), field1.begin());		
    ok_divergence = compare_scaling(kind, zDiv1_F2C.get_linalg_mat(), zDivOld_F2C.get_linalg_mat(), field1.begin());		
	    
    ok_divergence = compare_scaling(kind, xDiv2_F2C.get_linalg_mat(), xDivOld_F2C.get_linalg_mat(), field2.begin());		
    ok_divergence = compare_scaling(kind, yDiv2_F2C.get_linalg_mat(), yDivOld_F2C.get_linalg_mat(), field2.begin());		
    ok_divergence = compare_scaling(kind, zDiv2_F2C.get_linalg_mat(), zDivOld_F2C.get_linalg_mat(), field2.begin());		
  }
	
  // -----------------------------------------------------------------------------------------------------------------------------	
  // Final Results
  // -----------------------------------------------------------------------------------------------------------------------------	
  ok[0] = ok_linear_interpolant;      
  ok[1] = ok_gradient;
  ok[2] = ok_divergence;
}


// OPERATOR TEST 
void operator_test(int kind);

// CONVERGENCE TEST
void convergence_test(int kind);
 
// LEFT SCALE
void left_scale_test();

// RIGHT SCALE
void right_scale_test();

int main()
{
  const int size = 300;
  char comment[size];
	
  int iLinearInterpolant_C2F_Test;
  int iLinearInterpolant_F2C_Test;
  int iGradient_C2F_Test;
  int iGradient_F2C_Test;
  int iDivergence_C2F_Test;
  int iDivergence_F2C_Test; 	
  int iScratch_Laplacian_Test;
  int iMixedDerivatives_Test;
  int iLinearInterpolant_Convergence_C2F_Test;
  int iLinearInterpolant_Convergence_F2C_Test; 	
  int iGradient_Convergence_C2F_Test;
  int iGradient_Convergence_F2C_Test; 	
  int iDivergence_Convergence_C2F_Test;
  int iDivergence_Convergence_F2C_Test; 	
  int iScratch_Laplacian_Convergence_Test;
  int iMixedDerivatives_Convergence_Test;
  int iLeft_Scale_Test;
  int iRight_Scale_Test;
	
  ifstream iFile("src/test/test.input", ios::in);
		
  // 1. Operator - TESTS
  // -----------------------------------------------------------------------------------
  iFile.getline(comment, size);
  iFile.getline(comment, size);
  iFile.getline(comment, size);	
	
  iFile.getline(comment, size);	
  iFile >> iLinearInterpolant_C2F_Test;	iFile.getline(comment, size);
  iFile >> iLinearInterpolant_F2C_Test;	iFile.getline(comment, size);
  iFile.getline(comment, size);	
	
  iFile.getline(comment, size);		
  iFile >> iGradient_C2F_Test;			iFile.getline(comment, size);
  iFile >> iGradient_F2C_Test;			iFile.getline(comment, size);
  iFile.getline(comment, size);		
		
  iFile.getline(comment, size);		
  iFile >> iDivergence_C2F_Test;		iFile.getline(comment, size);
  iFile >> iDivergence_F2C_Test;		iFile.getline(comment, size); 	
  iFile.getline(comment, size);	

  iFile.getline(comment, size);	
  iFile >> iScratch_Laplacian_Test;		iFile.getline(comment, size);
  iFile.getline(comment, size);	

  iFile.getline(comment, size);		  	
  iFile >> iMixedDerivatives_Test;		iFile.getline(comment, size);
  iFile.getline(comment, size);	
	
  // 2. Operators - CONVERGENCE TESTS 
  // -----------------------------------------------------------------------------------	
  iFile.getline(comment, size);
  iFile.getline(comment, size);
  iFile.getline(comment, size);  	
	
  iFile.getline(comment, size);	
  iFile >> iLinearInterpolant_Convergence_C2F_Test;		iFile.getline(comment, size);
  iFile >> iLinearInterpolant_Convergence_F2C_Test; 		iFile.getline(comment, size);
  iFile.getline(comment, size);	

  iFile.getline(comment, size);	
  iFile >> iGradient_Convergence_C2F_Test;				iFile.getline(comment, size);
  iFile >> iGradient_Convergence_F2C_Test; 			iFile.getline(comment, size);
  iFile.getline(comment, size);	

  iFile.getline(comment, size);	
  iFile >> iDivergence_Convergence_C2F_Test;			iFile.getline(comment, size);
  iFile >> iDivergence_Convergence_F2C_Test; 			iFile.getline(comment, size);
  iFile.getline(comment, size);	
	
  iFile.getline(comment, size);	
  iFile >> iScratch_Laplacian_Convergence_Test;			iFile.getline(comment, size);
  iFile.getline(comment, size);	
	
  iFile.getline(comment, size);	  	
  iFile >> iMixedDerivatives_Convergence_Test;			iFile.getline(comment, size);
  iFile.getline(comment, size);	

  // 3. Spatial Fields - LEFT AND RIGHT SCALING
  // ----------------------------------------------------------------------------------- 	
  iFile.getline(comment, size);
  iFile.getline(comment, size);
  iFile.getline(comment, size);
  iFile >> iLeft_Scale_Test;								iFile.getline(comment, size);
  iFile >> iRight_Scale_Test;							iFile.getline(comment, size);
	
  iFile >> comment;
  iFile.close();
	
  if (strcmp(comment, "END"))
    {
      cout << "Error in file test.input!" << endl;
      exit(1);
    }

  std::cout << endl << endl;	
 
  // ------------------------------------------------------------------------
  // 1. OPERATORS - TESTS
  // ------------------------------------------------------------------------
  // LINEAR INTERPOLANT OPERATOR TEST
  if (iLinearInterpolant_C2F_Test==1)	operator_test(1);
  if (iLinearInterpolant_F2C_Test==1)  	operator_test(6);
						
  // GRADIENT OPERATOR TEST
  if (iGradient_C2F_Test==1)			operator_test(2);
  if (iGradient_F2C_Test==1)			operator_test(7);
		  	
  // DIVERGENCE OPERATOR TEST
  if (iDivergence_C2F_Test==1)		operator_test(3);
  if (iDivergence_F2C_Test==1)		operator_test(8);

  // SCRATCH OPERATOR TEST
  if (iScratch_Laplacian_Test==1)		operator_test(4);

  // MIXED DERIVATIVES TEST
  if (iMixedDerivatives_Test==1)		operator_test(5);
		 
  // ------------------------------------------------------------------------
  // 2. OPERATORS - CONVERGENCE TESTS
  // ------------------------------------------------------------------------
  // LINEAR INTERPOLANT OPERATOR TEST - CONVERGENCE
  if (iLinearInterpolant_Convergence_C2F_Test==1)	convergence_test(1);
  if (iLinearInterpolant_Convergence_F2C_Test==1)	convergence_test(6);

  // GRADIENT OPERATOR TEST - CONVERGENCE
  if (iGradient_Convergence_C2F_Test==1)			convergence_test(2);  	
  if (iGradient_Convergence_F2C_Test==1)			convergence_test(7);
		
  // DIVERGENCE OPERATOR TEST - CONVERGENCE
  if (iDivergence_Convergence_C2F_Test==1)		convergence_test(3);
  if (iDivergence_Convergence_F2C_Test==1)		convergence_test(8);
		  	
  // SCRATCH (LAPLACIAN) OPERATOR TEST - CONVERGENCE
  if (iScratch_Laplacian_Convergence_Test==1)		convergence_test(4);
	
  // MIXED DERIVATIVES OPERATOR TEST - CONVERGENCE
  if (iMixedDerivatives_Convergence_Test==1)		convergence_test(5);	
 
  // ------------------------------------------------------------------------
  // 3. SPATIAL FIELDS - LEFT AND RIGHT SCALING
  // ------------------------------------------------------------------------		 		 		
  // LEFT SCALE
  if (iLeft_Scale_Test==1)	left_scale_test();

  // RIGHT SCALE
  if (iRight_Scale_Test==1)	right_scale_test();
	
  // Return	
  return 0;
}

// 1. OPERATORS - TESTS
void operator_test(int kind)
{
  std::cout << "---------------------------------------------------------------"			<< endl;
  if (kind == 1) std::cout << " LINEAR INTERPOLANT OPERATOR - C2F"	<< endl;	
  if (kind == 2) std::cout << " GRADIENT OPERATOR - C2F"			<< endl;	
  if (kind == 3) std::cout << " DIVERGENCE OPERATOR - C2F"			<< endl;	
  if (kind == 4) std::cout << " SCRATCH LAPLACIAN OPERATOR"		<< endl;	
  if (kind == 5) std::cout << " MIXED DERIVATIVES OPERATOR"		<< endl;	
  if (kind == 6) std::cout << " LINEAR INTERPOLANT OPERATOR - F2C"	<< endl;	
  if (kind == 7) std::cout << " GRADIENT OPERATOR - F2C"			<< endl;	
  if (kind == 8) std::cout << " DIVERGENCE OPERATOR - F2C"			<< endl;	
  std::cout << "---------------------------------------------------------------" 			<< endl 	<< endl;
  	   
  bool ok = false;
  std::vector<int> dim(3,1);
  std::vector<double> mean_rel_err(3, 0.);
  std::vector<double> max_rel_err(3, 0.);
  	
  for (int nGridPoints=19, i=0; i<2; nGridPoints*=2, i++)
    {  	
      std::cout << "Grid: " << nGridPoints << " x " << nGridPoints << " x "<< nGridPoints << endl;
      dim[X_DIR] = nGridPoints; dim[Y_DIR] = nGridPoints; dim[Z_DIR] = nGridPoints;
  	 
      for (int iTestFunction=1; iTestFunction <=8; iTestFunction++)
	{
	  grid_class grid(dim);
	  analytical_class funct(grid);
	  funct.assign_function_name(iTestFunction);
	  funct.assign_function(iTestFunction);
	    
	  std::cout << "Test Function " << iTestFunction << ": " << funct.test.description() << endl;
	  std::cout << "---------------------------------------------------------------------------------" << endl;	
			
	  if (kind == 1) ok = test_linear_interpolant_C2F( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 2) ok = test_gradient_C2F          ( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 3) ok = test_divergence_C2F        ( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 4) ok = test_scratch               ( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 5) ok = test_mixed_derivatives     ( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 6) ok = test_linear_interpolant_F2C( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 7) ok = test_gradient_F2C          ( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 8) ok = test_divergence_F2C        ( grid, funct, mean_rel_err, max_rel_err );
			
	  std::cout << endl;
  	
	  if (kind == 1)
	    {
	      if( ok ) std::cout << "Linear Interpolant Operator C2F Test:   PASS" << endl;
	      else     std::cout << "Linear Interpolant Operator C2F Test:   FAIL" << endl;
	    }
	  if (kind == 2)
	    {
	      if( ok )  std::cout << "Gradient C2F Test:   PASS" << endl;
	      else     std::cout << "Gradient C2F test:   FAIL" << endl;
	    }
	  if (kind == 3)
	    {
	      if( ok ) std::cout << "Divergence Operator C2F Test:   PASS" << endl;
	      else     std::cout << "Divergence Operator C2F test:   FAIL" << endl;
	    }
	  if (kind == 4)
	    {
	      if( ok ) std::cout << "Scratch Laplacian Operator Test:   PASS" << endl;
	      else     std::cout << "Scratch Laplacian Operator test:   FAIL" << endl;
	    }		
	  if (kind == 5)
	    {
	      if( ok ) std::cout << "Mixed Derivatives Test:   PASS" << endl;
	      else     std::cout << "Mixed Derivatives test:   FAIL" << endl;
	    }		
	  if (kind == 6)
	    {
	      if( ok ) std::cout << "Linear Interpolant F2C Operator Test:   PASS" << endl;
	      else     std::cout << "Linear Interpolant F2C Operator Test:   FAIL" << endl;
	    }		
	  if (kind == 7)
	    {
	      if( ok )  std::cout << "Gradient F2C Test:   PASS" << endl;
	      else     std::cout << "Gradient F2C test:   FAIL" << endl;
	    }		
	  if (kind == 8)
	    {
	      if( ok )  std::cout << "Divergence Operator F2C Test:   PASS" << endl;
	      else     std::cout << "Divergence Operator F2C test:   FAIL" << endl;
	    }		
   			
	  std::cout << endl;
	}
    }
}

// OPERATORS - CONVERGENCE TESTS
void convergence_test(int kind)
{
  std::cout << "---------------------------------------------------------------"								<< endl;
  if (kind == 1) std::cout << " LINEAR INTERPOLANT OPERATOR - C2F - Convergence test"	<< endl;	
  if (kind == 2) std::cout << " GRADIENT OPERATOR - C2F - Convergence test"				<< endl;	
  if (kind == 3) std::cout << " DIVERGENCE OPERATOR - C2F - Convergence test"			<< endl;	
  if (kind == 4) std::cout << " SCRATCH LAPLACIAN OPERATOR - Convergence test"			<< endl;	
  if (kind == 5) std::cout << " SCRATCH MIXED DERIVATIVES - Convergence test"			<< endl;	
  if (kind == 6) std::cout << " LINEAR INTERPOLANT OPERATOR - F2C - Convergence test"	<< endl;	
  if (kind == 7) std::cout << " GRADIENT OPERATOR - F2C - Convergence test"				<< endl;	
  if (kind == 8) std::cout << " DIVERGENCE OPERATOR - F2C - Convergence test"			<< endl;	
  std::cout << "---------------------------------------------------------------"								<< endl	<< endl;
	 
  bool ok;
  std::vector<int> dim(3,1);
  std::vector<double> mean_rel_err(3, 0.);
  std::vector<double> max_rel_err(3, 0.);
  std::string fileName;
  std::string path; 
	
  if (kind == 1) path = "testfiles/convergence/linear_interpolant_C2F/TestFunction_";
  if (kind == 2) path = "testfiles/convergence/gradient_C2F/TestFunction_";
  if (kind == 3) path = "testfiles/convergence/divergence_C2F/TestFunction_";
  if (kind == 4) path = "testfiles/convergence/scratch_laplacian/TestFunction_";
  if (kind == 5) path = "testfiles/convergence/mixed_derivatives/TestFunction_";
  if (kind == 6) path = "testfiles/convergence/linear_interpolant_F2C/TestFunction_";
  if (kind == 7) path = "testfiles/convergence/gradient_F2C/TestFunction_";
  if (kind == 8) path = "testfiles/convergence/divergence_F2C/TestFunction_";
			  	    	
  for (int iTestFunction=1; iTestFunction <=8; iTestFunction++)
    {
      std::cout << "Test Function: " << iTestFunction << endl;
      std::cout << "---------------------------------------------------------------------------------" << endl;
			
      ofstream oFile;
      fileName = IntToString (iTestFunction);
      fileName = path + fileName;
      oFile.open(fileName.c_str(), ios::out);
      oFile.setf(ios::scientific);
      oFile << "Points" << "\t";
      oFile << "Mean Rel X" << "\t\t";
      oFile << "Mean Rel Y" << "\t\t";
      oFile << "Mean Rel Z" << "\t\t";
      oFile << "Max Rel X" << "\t\t";
      oFile << "Max Rel Y" << "\t\t";
      oFile << "Max Rel Z" << "\t\t";
      oFile << endl << endl;
      
      for (int nGridPoints=6, i=1; i<5; nGridPoints*=2, i++)
	{  			 
	  std::cout << "Grid: " << nGridPoints << " x " << nGridPoints << " x "<< nGridPoints << endl;
  		  	
	  dim[X_DIR] = nGridPoints; dim[Y_DIR] = nGridPoints; dim[Z_DIR] = nGridPoints;

	  grid_class grid( dim);
	  analytical_class funct( grid );
	  funct.assign_function_name(iTestFunction);
	  funct.assign_function(iTestFunction);
	
	  if (kind == 1) ok = test_linear_interpolant_C2F( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 2) ok = test_gradient_C2F          ( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 3) ok = test_divergence_C2F        ( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 4) ok = test_scratch               ( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 5) ok = test_mixed_derivatives     ( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 6) ok = test_linear_interpolant_F2C( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 7) ok = test_gradient_F2C          ( grid, funct, mean_rel_err, max_rel_err );
	  if (kind == 8) ok = test_divergence_F2C        ( grid, funct, mean_rel_err, max_rel_err );
  	
	  std::cout << endl;
  				
	  for (int k=0; k <3; k++)
	    {
	      if (mean_rel_err[k]==0.) mean_rel_err[k]=1.e-17;
	      if (max_rel_err[k]==0.) max_rel_err[k]=1.e-17;
	    }	
	  oFile << nGridPoints << "\t";
	  oFile << "\t" <<  mean_rel_err[X_DIR] <<  "\t" << mean_rel_err[Y_DIR] <<  "\t" << mean_rel_err[Z_DIR]; 
	  oFile << "\t" <<  max_rel_err[X_DIR]  <<  "\t" << max_rel_err[Y_DIR]  <<  "\t" << max_rel_err[Z_DIR];
	  oFile << endl;
	}
			
      oFile.close();
    }
}

// LEFT SCALE
void left_scale_test()
{
  std::cout << "---------------------------------------------------------------------" << endl;
  std::cout << " LEFT SCALING" << endl;	
  std::cout << "---------------------------------------------------------------------" << endl << endl;
	
  std::vector<bool> ok(3,true);
  std::vector<string> ok_string(3,"FAIL");
  std::vector<int>  dim(3,1);
	  
  ofstream oFile;
  oFile.open("testfiles/scaling/Left_Scale", ios::out);
  oFile.setf(ios::scientific);

  oFile << "#x\t" << "#y\t" << "#z\t" << "Interp.\t" << "Grad.\t" << "Diver.\t" << endl << endl;

  for (int nX=3, i=1; i<=3; nX*=2, i++)
    for (int nY=2, j=1; j<=3; nY*=2, j++)
      for (int nZ=4, k=1; k<=3; nZ*=2, k++)
	{  			 
	  std::cout << "Grid: " << nX << " x " << nY << " x "<< nZ << endl;
	  dim[X_DIR] = nX; dim[Y_DIR] = nY; dim[Z_DIR] = nZ;
	  			  						
	  test_scale_C2F("LEFT", dim, ok);
  			
	  if (ok[0] == true)	 	{std::cout << "Left Scaling C2F - Linear Interpolant: PASS" << endl; ok_string[0] = "PASS";}
	  if (ok[0] == false)		 std::cout << "Left Scaling C2F - Linear Interpolant: FAIL" << endl;

	  if (ok[1] == true)	 	{std::cout << "Left Scaling C2F - Gradient: PASS" << endl; ok_string[1] = "PASS";}
	  if (ok[1] == false)		 std::cout << "Left Scaling C2F - Gradient: FAIL" << endl;

	  if (ok[2] == true)	 	{std::cout << "Left Scaling C2F - Divergence: PASS" << endl; ok_string[2] = "PASS";}
	  if (ok[2] == false)		 std::cout << "Left Scaling C2F - Divergence: FAIL" << endl;
  		
	  std::cout << endl;
  		
	  oFile << nX << "\t" << nY << "\t" << nZ << "\t";
	  oFile << ok_string[0] << "\t" << ok_string[1] << "\t" << ok_string[2] << "\t";
	  oFile << endl;
				  				
	  test_scale_F2C("LEFT", dim, ok);
  			
	  if (ok[0] == true)	 	{std::cout << "Left Scaling F2C - Linear Interpolant: PASS" << endl; ok_string[0] = "PASS";}
	  if (ok[0] == false)		 std::cout << "Left Scaling F2C - Linear Interpolant: FAIL" << endl;

	  if (ok[1] == true)	 	{std::cout << "Left Scaling F2C - Gradient: PASS" << endl; ok_string[1] = "PASS";}
	  if (ok[1] == false)		 std::cout << "Left Scaling F2C - Gradient: FAIL" << endl;

	  if (ok[2] == true)	 	{std::cout << "Left Scaling F2C - Divergence: PASS" << endl; ok_string[2] = "PASS";}
	  if (ok[2] == false)		 std::cout << "Left Scaling F2C - Divergence: FAIL" << endl;
  		
	  std::cout << endl;
  		
	  oFile << nX << "\t" << nY << "\t" << nZ << "\t";
	  oFile << ok_string[0] << "\t" << ok_string[1] << "\t" << ok_string[2] << "\t";
	  oFile << endl;
	}
			
  oFile << endl;
  oFile.close();
}

// RIGHT SCALE
void right_scale_test()
{
  std::cout << "---------------------------------------------------------------------" << endl;
  std::cout << " RIGHT SCALING" << endl;	
  std::cout << "---------------------------------------------------------------------" << endl << endl;
	
  std::vector<bool> ok(3,true);
  std::vector<string> ok_string(3,"FAIL");
  std::vector<int>  dim(3,1);
	  
  ofstream oFile;
  oFile.open("testfiles/scaling/Right_Scale", ios::out);
  oFile.setf(ios::scientific);
	
  oFile << "#x\t" << "#y\t" << "#z\t" << "Interp.\t" << "Grad.\t" << "Diver.\t" << endl << endl;

  for (int nX=3, i=1; i<=3; nX*=2, i++)
    for (int nY=2, j=1; j<=3; nY*=2, j++)
      for (int nZ=4, k=1; k<=3; nZ*=2, k++)
	{  			 
	  std::cout << "Grid: " << nX << " x " << nY << " x "<< nZ << endl;
	  dim[X_DIR] = nX; dim[Y_DIR] = nY; dim[Z_DIR] = nZ;
	  			  						
	  test_scale_C2F("RIGHT", dim, ok);
  			
	  if (ok[0] == true)	 	{std::cout << "Right Scaling C2F - Linear Interpolant: PASS" << endl; ok_string[0] = "PASS";}
	  if (ok[0] == false)		 std::cout << "Right Scaling C2F - Linear Interpolant: FAIL" << endl;

	  if (ok[1] == true)	 	{std::cout << "Right Scaling C2F - Gradient: PASS" << endl; ok_string[1] = "PASS";}
	  if (ok[1] == false)		 std::cout << "Right Scaling C2F - Gradient: FAIL" << endl;

	  if (ok[2] == true)	 	{std::cout << "Right Scaling C2F - Divergence: PASS" << endl; ok_string[2] = "PASS";}
	  if (ok[2] == false)		 std::cout << "Right Scaling C2F - Divergence: FAIL" << endl;
  		
	  std::cout << endl;
  		
	  oFile << nX << "\t" << nY << "\t" << nZ << "\t";
	  oFile << ok_string[0] << "\t" << ok_string[1] << "\t" << ok_string[2] << "\t";
	  oFile << endl;
				
	  test_scale_F2C("RIGHT", dim, ok);
  			
	  if (ok[0] == true)	 	{std::cout << "Right Scaling F2C - Linear Interpolant: PASS" << endl; ok_string[0] = "PASS";}
	  if (ok[0] == false)		 std::cout << "Right Scaling F2C - Linear Interpolant: FAIL" << endl;

	  if (ok[1] == true)	 	{std::cout << "Right Scaling F2C - Gradient: PASS" << endl; ok_string[1] = "PASS";}
	  if (ok[1] == false)		 std::cout << "Right Scaling F2C - Gradient: FAIL" << endl;

	  if (ok[2] == true)	 	{std::cout << "Right Scaling F2C - Divergence: PASS" << endl; ok_string[2] = "PASS";}
	  if (ok[2] == false)		 std::cout << "Right Scaling F2C- Divergence: FAIL" << endl;
  		
	  std::cout << endl;
  		
	  oFile << nX << "\t" << nY << "\t" << nZ << "\t";
	  oFile << ok_string[0] << "\t" << ok_string[1] << "\t" << ok_string[2] << "\t";
	  oFile << endl;				
	}
			
  oFile << endl;
  oFile.close();
}



