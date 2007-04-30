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

  CellField d2fdx2; 		// Analitycal Laplacian XX
  CellField d2fdy2; 		// Analitycal Laplacian YY
  CellField d2fdz2;		// Analitycal Laplacian ZZ

  CellField d2fdxdy; 		// Analitycal Mixed Derivative XY
  CellField d2fdxdz; 		// Analitycal Mixed Derivative XZ
  CellField d2fdydz;		// Analitycal Mixed Derivative YZ
  CellField d2fdydx; 		// Analitycal Mixed Derivative YX
  CellField d2fdzdx; 		// Analitycal Mixed Derivative ZX
  CellField d2fdzdy;		// Analitycal Mixed Derivative ZY

  // -----------------------------------------------------------------------------------
  // Cell Face
  // -----------------------------------------------------------------------------------
  XSideField f_faceX; 			// Interpolated function X
  YSideField f_faceY; 			// Interpolated function Y
  ZSideField f_faceZ;			// Interpolated function Z

  XSideField dfdx_faceX; 			// Analitycal Gradient X
  YSideField dfdy_faceY; 			// Analitycal Gradient Y
  ZSideField dfdz_faceZ;			// Analitycal Gradient Z

  XSideField d2fdx2_faceX; 		// Analitycal Laplacian XX
  YSideField d2fdy2_faceY; 		// Analitycal Laplacian YY
  ZSideField d2fdz2_faceZ;		// Analitycal Laplacian ZZ
 
	
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

analytical_class::analytical_class( const grid_class& grid )  : 

	f   ( grid.dim, NULL, InternalStorage ),

    	dfdx( grid.dim, NULL, InternalStorage ),
    	dfdy( grid.dim, NULL, InternalStorage ),
    	dfdz( grid.dim, NULL, InternalStorage ),

    	d2fdx2( grid.dim, NULL, InternalStorage ),
    	d2fdy2( grid.dim, NULL, InternalStorage ),
    	d2fdz2( grid.dim, NULL, InternalStorage ),

    	d2fdxdy( grid.dim, NULL, InternalStorage ),
    	d2fdxdz( grid.dim, NULL, InternalStorage ),
    	d2fdydz( grid.dim, NULL, InternalStorage ),
    	d2fdydx( grid.dim, NULL, InternalStorage ),
    	d2fdzdx( grid.dim, NULL, InternalStorage ),
    	d2fdzdy( grid.dim, NULL, InternalStorage ),

    	f_faceX( grid.dim, NULL, InternalStorage ),
    	f_faceY( grid.dim, NULL, InternalStorage ),
    	f_faceZ( grid.dim, NULL, InternalStorage ),

    	dfdx_faceX( grid.dim, NULL, InternalStorage ),
   	dfdy_faceY( grid.dim, NULL, InternalStorage ),
    	dfdz_faceZ( grid.dim, NULL, InternalStorage ),

    	d2fdx2_faceX( grid.dim, NULL, InternalStorage ),
    	d2fdy2_faceY( grid.dim, NULL, InternalStorage ),
    	d2fdz2_faceZ( grid.dim, NULL, InternalStorage ),

    	grid_( grid )
{

}


void analytical_class::assign_function_name(int kind)
{		
  if (kind == 1)	       test.setup(&testFunction01);	
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
	int point;
	assign_function_name(kind);	
  	
	// Cell Fields
	point = 0;  	
  	for( int k=0; k<grid_.nz_cell; ++k )
		for( int j=0; j<grid_.ny_cell; ++j)
			for( int i=0; i<grid_.nx_cell; ++i )
	   		 {	
	      			f[point] = test.f(grid_.x_cell[point], grid_.y_cell[point], grid_.z_cell[point]);

	      			dfdx[point] = test.dfdx(grid_.x_cell[point], grid_.y_cell[point], grid_.z_cell[point]);
	      			dfdy[point] = test.dfdy(grid_.x_cell[point], grid_.y_cell[point], grid_.z_cell[point]);
	      			dfdz[point] = test.dfdz(grid_.x_cell[point], grid_.y_cell[point], grid_.z_cell[point]);
								
			      	d2fdx2[point] = test.d2fdx2(grid_.x_cell[point], grid_.y_cell[point], grid_.z_cell[point]);
			      	d2fdy2[point] = test.d2fdy2(grid_.x_cell[point], grid_.y_cell[point], grid_.z_cell[point]);	
			     	d2fdz2[point] = test.d2fdz2(grid_.x_cell[point], grid_.y_cell[point], grid_.z_cell[point]);					
				
			     	d2fdxdy[point] = test.d2fdxdy(grid_.x_cell[point], grid_.y_cell[point], grid_.z_cell[point]);
			      	d2fdxdz[point] = test.d2fdxdz(grid_.x_cell[point], grid_.y_cell[point], grid_.z_cell[point]);	
			      	d2fdydz[point] = test.d2fdydz(grid_.x_cell[point], grid_.y_cell[point], grid_.z_cell[point]);

			     	d2fdydx[point] = d2fdxdy[point];
			      	d2fdzdx[point] = d2fdxdz[point];	
			      	d2fdzdy[point] = d2fdydz[point];
																
				++point;
	    		}
				
	point = 0;  	
  	for( int k=0; k<grid_.nz_faceX; ++k )
		for( int j=0; j<grid_.ny_faceX; ++j)
			for( int i=0; i<grid_.nx_faceX; ++i )
	   		 {	
			 	f_faceX[point] 		= test.f(grid_.x_faceX[point], grid_.y_faceX[point], grid_.z_faceX[point]);
	      			dfdx_faceX[point] 	= test.dfdx(grid_.x_faceX[point], grid_.y_faceX[point], grid_.z_faceX[point]);
			      	d2fdx2_faceX[point]	= test.d2fdx2(grid_.x_faceX[point], grid_.y_faceX[point], grid_.z_faceX[point]);
														
				++point;
	    		}
			
	point = 0;  	
  	for( int k=0; k<grid_.nz_faceY; ++k )
		for( int j=0; j<grid_.ny_faceY; ++j)
			for( int i=0; i<grid_.nx_faceY; ++i )
	   		 {	
	      			f_faceY[point] = test.f(grid_.x_faceY[point], grid_.y_faceY[point], grid_.z_faceY[point]);			 
	      			dfdy_faceY[point] = test.dfdy(grid_.x_faceY[point], grid_.y_faceY[point], grid_.z_faceY[point]);
			        d2fdy2_faceY[point] = test.d2fdy2(grid_.x_faceY[point], grid_.y_faceY[point], grid_.z_faceY[point]);
												
				++point;
	    		}
			
	point = 0;  	
  	for( int k=0; k<grid_.nz_faceZ; ++k )
		for( int j=0; j<grid_.ny_faceZ; ++j)
			for( int i=0; i<grid_.nx_faceZ; ++i )
	   		 {	
	      			f_faceZ[point] = test.f(grid_.x_faceZ[point], grid_.y_faceZ[point], grid_.z_faceZ[point]);							 
	      			dfdz_faceZ[point] = test.dfdz(grid_.x_faceZ[point], grid_.y_faceZ[point], grid_.z_faceZ[point]);
			       d2fdz2_faceZ[point] = test.d2fdz2(grid_.x_faceZ[point], grid_.y_faceZ[point], grid_.z_faceZ[point]);
								
				++point;
	    		}				
}
