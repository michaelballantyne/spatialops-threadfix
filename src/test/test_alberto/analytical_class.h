#include <FV2ndOrderTypes.h>

// *****************************************************************************************************************************	
// Analytical Function
// *****************************************************************************************************************************	
class analytical_class
{
public:

  // -----------------------------------------------------------------------------------		
  // Cell Center 
  // -----------------------------------------------------------------------------------
  CellField f;   			// Function
  	
  CellField dfdx; 			// Analitycal Gradient X
  CellField dfdy; 			// Analitycal Gradient Y
  CellField dfdz;			// Analitycal Gradient Z

  CellField d2fdx2; 			// Analitycal Laplacian XX
  CellField d2fdy2; 			// Analitycal Laplacian YY
  CellField d2fdz2;			// Analitycal Laplacian ZZ

  CellField d2fdxdy; 		// Analitycal Mixed Derivative XY
  CellField d2fdxdz; 		// Analitycal Mixed Derivative XZ
  CellField d2fdydz;			// Analitycal Mixed Derivative YZ

  // -----------------------------------------------------------------------------------
  // Cell Face
  // -----------------------------------------------------------------------------------
  XSideField f_xint; 			// Interpolated function X
  YSideField f_yint; 			// Interpolated function Y
  ZSideField f_zint;			// Interpolated function Z

  XSideField dfdx_int; 		// Analitycal Gradient X
  YSideField dfdy_int; 		// Analitycal Gradient Y
  ZSideField dfdz_int;		// Analitycal Gradient Z

  XSideField d2fdx2_xint; 			// Analitycal Laplacian XX
  YSideField d2fdy2_yint; 			// Analitycal Laplacian YY
  ZSideField d2fdz2_zint;			// Analitycal Laplacian ZZ
	
  // -----------------------------------------------------------------------------------
  // Public Functions
  // -----------------------------------------------------------------------------------	
  analytical_class( const grid_class & grid );

  void assign_function(int kind);
  void assign_function_name(int kind);

  AnalyticalFunction test;

private:
  const grid_class &grid_;
	
  TestFunction01 testFunction01;
  TestFunction02 testFunction02;
  TestFunction03 testFunction03;
  TestFunction04 testFunction04;
  TestFunction05 testFunction05;
  TestFunction06 testFunction06;
  TestFunction07 testFunction07;
  TestFunction08 testFunction08;
};

analytical_class::analytical_class( const grid_class& grid )
  : f   ( grid.dim, NULL, InternalStorage ),

    dfdx( grid.dim, NULL, InternalStorage ),
    dfdy( grid.dim, NULL, InternalStorage ),
    dfdz( grid.dim, NULL, InternalStorage ),

    d2fdx2( grid.dim, NULL, InternalStorage ),
    d2fdy2( grid.dim, NULL, InternalStorage ),
    d2fdz2( grid.dim, NULL, InternalStorage ),

    d2fdxdy( grid.dim, NULL, InternalStorage ),
    d2fdxdz( grid.dim, NULL, InternalStorage ),
    d2fdydz( grid.dim, NULL, InternalStorage ),

    f_xint( grid.dim, NULL, InternalStorage ),
    f_yint( grid.dim, NULL, InternalStorage ),
    f_zint( grid.dim, NULL, InternalStorage ),

    dfdx_int( grid.dim, NULL, InternalStorage ),
    dfdy_int( grid.dim, NULL, InternalStorage ),
    dfdz_int( grid.dim, NULL, InternalStorage ),

    d2fdx2_xint( grid.dim, NULL, InternalStorage ),
    d2fdy2_yint( grid.dim, NULL, InternalStorage ),
    d2fdz2_zint( grid.dim, NULL, InternalStorage ),

    grid_( grid )
{
}


void analytical_class::assign_function_name(int kind)
{		
  if (kind == 1)	test.setup(&testFunction01);	
  else if (kind == 2)	test.setup(&testFunction02);	
  else if (kind == 3)	test.setup(&testFunction03);	
  else if (kind == 4)	test.setup(&testFunction04);	
  else if (kind == 5)	test.setup(&testFunction05);	
  else if (kind == 6)	test.setup(&testFunction06);	
  else if (kind == 7)	test.setup(&testFunction07);	
  else if (kind == 8)	test.setup(&testFunction08);
}

void analytical_class::assign_function(int kind)
{		
  assign_function_name(kind);	
  	
  int point = 0;  	
  for( int k=0; k<grid_.nz; ++k )
    {
      for( int j=0; j<grid_.ny; ++j)
    	{
	  for( int i=0; i<grid_.nx; ++i )
	    {	
	      f[point] = test.f(grid_.x[point], grid_.y[point], grid_.z[point]);

	      f_xint[point] = test.f(grid_.x_int[point], grid_.y[point], grid_.z[point]);
	      f_yint[point] = test.f(grid_.x[point], grid_.y_int[point], grid_.z[point]);
	      f_zint[point] = test.f(grid_.x[point], grid_.y[point], grid_.z_int[point]);				

	      dfdx[point] = test.dfdx(grid_.x[point], grid_.y[point], grid_.z[point]);
	      dfdy[point] = test.dfdy(grid_.x[point], grid_.y[point], grid_.z[point]);
	      dfdz[point] = test.dfdz(grid_.x[point], grid_.y[point], grid_.z[point]);

	      dfdx_int[point] = test.dfdx(grid_.x_int[point], grid_.y[point], grid_.z[point]);
	      dfdy_int[point] = test.dfdy(grid_.x[point], grid_.y_int[point], grid_.z[point]);
	      dfdz_int[point] = test.dfdz(grid_.x[point], grid_.y[point], grid_.z_int[point]);

	      d2fdx2[point] = test.d2fdx2(grid_.x[point], grid_.y[point], grid_.z[point]);
	      d2fdy2[point] = test.d2fdy2(grid_.x[point], grid_.y[point], grid_.z[point]);
	      d2fdz2[point] = test.d2fdz2(grid_.x[point], grid_.y[point], grid_.z[point]);

	      d2fdx2_xint[point] = test.d2fdx2(grid_.x_int[point], grid_.y[point], grid_.z[point]);
	      d2fdy2_yint[point] = test.d2fdy2(grid_.x[point], grid_.y_int[point], grid_.z[point]);
	      d2fdz2_zint[point] = test.d2fdz2(grid_.x[point], grid_.y[point], grid_.z_int[point]);
				
	      d2fdxdy[point] = test.d2fdxdy(grid_.x[point], grid_.y[point], grid_.z[point]);
	      d2fdxdz[point] = test.d2fdxdz(grid_.x[point], grid_.y[point], grid_.z[point]);
	      d2fdydz[point] = test.d2fdydz(grid_.x[point], grid_.y[point], grid_.z[point]);
												
	      ++point;
	    }
    	} 
    }
}
