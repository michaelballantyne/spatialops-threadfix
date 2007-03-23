// *****************************************************************************************************************************	
// Analytical Function
// *****************************************************************************************************************************	
class analytical_class
{
public:

	// -----------------------------------------------------------------------------------		
	// Cell Center 
	// -----------------------------------------------------------------------------------
  	vector<double> f;   			// Function
  	
  	vector<double> dfdx; 			// Analitycal Gradient X
	vector<double> dfdy; 			// Analitycal Gradient Y
	vector<double> dfdz;			// Analitycal Gradient Z

  	vector<double> d2fdx2; 			// Analitycal Laplacian XX
	vector<double> d2fdy2; 			// Analitycal Laplacian YY
	vector<double> d2fdz2;			// Analitycal Laplacian ZZ

  	vector<double> d2fdxdy; 		// Analitycal Mixed Derivative XY
	vector<double> d2fdxdz; 		// Analitycal Mixed Derivative XZ
	vector<double> d2fdydz;			// Analitycal Mixed Derivative YZ

	// -----------------------------------------------------------------------------------
	// Cell Face
	// -----------------------------------------------------------------------------------
  	vector<double> f_xint; 			// Interpolated function X
	vector<double> f_yint; 			// Interpolated function Y
	vector<double> f_zint;			// Interpolated function Z

  	vector<double> dfdx_int; 		// Analitycal Gradient X
	vector<double> dfdy_int; 		// Analitycal Gradient Y
	vector<double> dfdz_int;		// Analitycal Gradient Z

  	vector<double> d2fdx2_xint; 			// Analitycal Laplacian XX
	vector<double> d2fdy2_yint; 			// Analitycal Laplacian YY
	vector<double> d2fdz2_zint;			// Analitycal Laplacian ZZ
	
	// -----------------------------------------------------------------------------------
	// Public Functions
	// -----------------------------------------------------------------------------------	
	void setup(grid_class *grid);
	void assign_function(int kind);
	void assign_function_name(int kind);

	AnalyticalFunction test;

private:
	grid_class *ptGrid;
	
	TestFunction01 testFunction01;
	TestFunction02 testFunction02;
	TestFunction03 testFunction03;
	TestFunction04 testFunction04;
	TestFunction05 testFunction05;
	TestFunction06 testFunction06;
	TestFunction07 testFunction07;
	TestFunction08 testFunction08;
};

void analytical_class::setup(grid_class *grid)
{
	ptGrid = grid;
  	
  	f.resize(grid->nn,0.0);				// Function - Cell Center
  	dfdx.resize(grid->nn,0.0);			// Analitycal Gradient X - Cell Center
	dfdy.resize(grid->nn,0.0);			// Analitycal Gradient Y - Cell Center
	dfdz.resize(grid->nn,0.0);			// Analitycal Gradient Z - Cell Center	

  	d2fdx2.resize(grid->nn,0.0);		// Analitycal Laplacian X - Cell Center
	d2fdy2.resize(grid->nn,0.0);		// Analitycal Laplacian Y - Cell Center
	d2fdz2.resize(grid->nn,0.0);		// Analitycal Laplacian Z - Cell Center 
	
	d2fdxdy.resize(grid->nn,0.0);		// Analitycal Mixed Derivative X - Cell Center
	d2fdxdz.resize(grid->nn,0.0);		// Analitycal Mixed Derivative Y - Cell Center
	d2fdydz.resize(grid->nn,0.0);		// Analitycal Mixed Derivative Z - Cell Center 

  	dfdx_int.resize(grid->nn,0.0);			// Analitycal Gradient X - Cell Face
	dfdy_int.resize(grid->nn,0.0);			// Analitycal Gradient Y - Cell Face
	dfdz_int.resize(grid->nn,0.0);			// Analitycal Gradient Z - Cell Face	

  	f_xint.resize(grid->nn,0.0);				// Interpolated function X - Cell Face
	f_yint.resize(grid->nn,0.0);				// Interpolated function Y - Cell Face
	f_zint.resize(grid->nn,0.0);				// Interpolated function Z - Cell Face	
	
  	d2fdx2_xint.resize(grid->nn,0.0);		// Analitycal Laplacian X - Cell Center
	d2fdy2_yint.resize(grid->nn,0.0);		// Analitycal Laplacian Y - Cell Center
	d2fdz2_zint.resize(grid->nn,0.0);		// Analitycal Laplacian Z - Cell Center 
	
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
  	for( int k=0; k<ptGrid->nz; ++k )
  	{
    	for( int j=0; j<ptGrid->ny; ++j)
    	{
    		for( int i=0; i<ptGrid->nx; ++i )
    		{	
				f[point] = test.f(ptGrid->x[point], ptGrid->y[point], ptGrid->z[point]);

				f_xint[point] = test.f(ptGrid->x_int[point], ptGrid->y[point], ptGrid->z[point]);
				f_yint[point] = test.f(ptGrid->x[point], ptGrid->y_int[point], ptGrid->z[point]);
				f_zint[point] = test.f(ptGrid->x[point], ptGrid->y[point], ptGrid->z_int[point]);				

				dfdx[point] = test.dfdx(ptGrid->x[point], ptGrid->y[point], ptGrid->z[point]);
				dfdy[point] = test.dfdy(ptGrid->x[point], ptGrid->y[point], ptGrid->z[point]);
				dfdz[point] = test.dfdz(ptGrid->x[point], ptGrid->y[point], ptGrid->z[point]);

				dfdx_int[point] = test.dfdx(ptGrid->x_int[point], ptGrid->y[point], ptGrid->z[point]);
				dfdy_int[point] = test.dfdy(ptGrid->x[point], ptGrid->y_int[point], ptGrid->z[point]);
				dfdz_int[point] = test.dfdz(ptGrid->x[point], ptGrid->y[point], ptGrid->z_int[point]);

				d2fdx2[point] = test.d2fdx2(ptGrid->x[point], ptGrid->y[point], ptGrid->z[point]);
				d2fdy2[point] = test.d2fdy2(ptGrid->x[point], ptGrid->y[point], ptGrid->z[point]);
				d2fdz2[point] = test.d2fdz2(ptGrid->x[point], ptGrid->y[point], ptGrid->z[point]);

				d2fdx2_xint[point] = test.d2fdx2(ptGrid->x_int[point], ptGrid->y[point], ptGrid->z[point]);
				d2fdy2_yint[point] = test.d2fdy2(ptGrid->x[point], ptGrid->y_int[point], ptGrid->z[point]);
				d2fdz2_zint[point] = test.d2fdz2(ptGrid->x[point], ptGrid->y[point], ptGrid->z_int[point]);
				
				d2fdxdy[point] = test.d2fdxdy(ptGrid->x[point], ptGrid->y[point], ptGrid->z[point]);
				d2fdxdz[point] = test.d2fdxdz(ptGrid->x[point], ptGrid->y[point], ptGrid->z[point]);
				d2fdydz[point] = test.d2fdydz(ptGrid->x[point], ptGrid->y[point], ptGrid->z[point]);
												
				++point;
      		}
    	} 
  	}
}
