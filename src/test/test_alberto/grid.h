// *****************************************************************************************************************************	
// Grid Class
// *****************************************************************************************************************************	
double LA = 1.;				
double LB = 2.;

class grid_class
{
public:
	int nx;
	int ny;
	int nz;
	int nn;
	int ni;
	std::vector<int> 	dim;
	std::vector<double> spacing;
	std::vector<double> area;
	std::vector<int> 	nghostCell;
	
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	
	std::vector<double> x_int;
	std::vector<double> y_int;
	std::vector<double> z_int;
	
	double volume;	
	
	std::vector<int> index_bothsides;
	std::vector<int> index_bothsides_x;
	std::vector<int> index_bothsides_y;
	std::vector<int> index_bothsides_z;
	std::vector<int> index_noside;
	std::vector<int> index_noside_x;
	std::vector<int> index_noside_y;
	std::vector<int> index_noside_z;
	
	void setup_geom(std::vector<int> _dim, std::vector<int> _nghostCell);
	
private:

};

void grid_class::setup_geom(std::vector<int> _dim, std::vector<int> _nghostCell)
{
	// Staggered grid from 0 to PI in each direction

	// #point           0           1           2                       n          n+1
	   
	// Grid:      |-----O-----|-----O-----|-----O-----|-----O-----|-----O-----|-----O-----|
	
	// #face      0           1           2                       n          n+1
	
	// coord.                 0          PI/n       2PI/n                     PI        
	

  	// -----------------------------------------------------------------------------------------------------------------------------	
  	// Resizing vectors
  	// -----------------------------------------------------------------------------------------------------------------------------			
  	dim.resize(3,1);				// number of points in each direction
  	spacing.resize(3,0.0);			// spacing in each direction
  	area.resize(3,0.0);			// area in each direction
  	nghostCell.resize(6,1);		// number of ghost cells in each direction (X- X+  Y- Y+  Z- Z+)
  
  	
  	// -----------------------------------------------------------------------------------------------------------------------------	
  	// Set up array dimensions for "source" arrays
  	// -----------------------------------------------------------------------------------------------------------------------------	
  	dim[X_DIR] = _dim[0];
 	dim[Y_DIR] = _dim[1];
 	dim[Z_DIR] = _dim[2];
 	for( int i=0; i<5; ++i )
 		nghostCell[i] = _nghostCell[i];
 		
 	nx =                		      dim[X_DIR] + nghostCell[0] + nghostCell[1];
  	ny = dim[Y_DIR]>1 ? dim[Y_DIR] + nghostCell[2] + nghostCell[3] : 1;
 	nz = dim[Z_DIR]>1 ? dim[Z_DIR] + nghostCell[4] + nghostCell[5] : 1;
  	nn = nx*ny*nz;
  	ni = dim[X_DIR]*dim[Y_DIR]*dim[Z_DIR];
  	

  	for( int i=0; i<3; ++i )
    	spacing[i] = (dim[i]==1) ? 1.0 : (LB-LA)/(dim[i]);

  	area[X_DIR] = spacing[Y_DIR]*spacing[Z_DIR];
  	area[Y_DIR] = spacing[X_DIR]*spacing[Z_DIR];
  	area[Z_DIR] = spacing[X_DIR]*spacing[Y_DIR];

  	volume = spacing[X_DIR]*spacing[Y_DIR]*spacing[Z_DIR];
  	
  	// -----------------------------------------------------------------------------------------------------------------------------	
  	// Grid Coordinates
  	// -----------------------------------------------------------------------------------------------------------------------------	
  	x.resize(nn,0.0);		// xCoordinates
  	y.resize(nn,0.0);		// yCoordinates
  	z.resize(nn,0.0);		// zCoordinates
  	x_int.resize(nn,0.0);	// xCoordinates - interpolated
  	y_int.resize(nn,0.0);	// yCoordinates - interpolated
  	z_int.resize(nn,0.0);	// zCoordinates - interpolated

  	int point=0;	
  	for( int k=0; k<nz; ++k )
    	for( int j=0; j<ny; ++j)
    		for( int i=0; i<nx; ++i )
    		{
				x[point] = double(i-0.5 + (1-nghostCell[0]))*spacing[0] + LA;	// xCoordinates 				
   				y[point] = double(j-0.5 + (1-nghostCell[2]))*spacing[1] + LA;	// yCoordinates 
   				z[point] = double(k-0.5 + (1-nghostCell[4]))*spacing[2] + LA;	// zCoordinates
   				
   				x_int[point] = double(i-1. + (1-nghostCell[0]))*spacing[0] + LA;	// xCoordinates 				
   				y_int[point] = double(j-1. + (1-nghostCell[2]))*spacing[1] + LA;	// yCoordinates 
   				z_int[point] = double(k-1. + (1-nghostCell[4]))*spacing[2] + LA;	// zCoordinates
   				
				++point;
      		}
      		
	// skip BOTH sides
	for( int k=1; k<=nz; ++k )
    	for( int j=1; j<=ny; ++j)
      		for( int i=1; i<=nx; ++i )
      		{
      			if (i>nghostCell[0] && i<=(nx-nghostCell[1]))
      				if (j>nghostCell[2] && j<=(ny-nghostCell[3]))
      					if (k>nghostCell[4] && k<=(nz-nghostCell[5]))
      					{	
      						index_bothsides.push_back((i-1)+(j-1)*nx+(k-1)*nx*ny); 
      						index_bothsides_x.push_back(i); index_bothsides_y.push_back(j); index_bothsides_z.push_back(k);
      					}									
      		}
		
    	// no skipping
	for( int k=1; k<=nz; ++k )
    	for( int j=1; j<=ny; ++j)
      		for( int i=1; i<=nx; ++i )  	
      		{
      			index_noside.push_back((i-1)+(j-1)*nx+(k-1)*nx*ny); 
      			index_noside_x.push_back(i); index_noside_y.push_back(j); index_noside_z.push_back(k);
      		}
}







