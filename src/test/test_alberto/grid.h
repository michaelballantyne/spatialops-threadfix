#include <vector>

// *****************************************************************************************************************************	
// Grid Class
// *****************************************************************************************************************************	
double LA = 1.;				
double LB = 2.;

class grid_class
{
public:
	std::vector<int> 	dim;

  	int nx_cell, ny_cell, nz_cell;
  	int nx_faceX, ny_faceX, nz_faceX;
  	int nx_faceY, ny_faceY, nz_faceY;
  	int nx_faceZ, ny_faceZ, nz_faceZ;
	
  	int ntot_cell, nint_cell;
  	int ntot_faceX, nint_faceX; 
  	int ntot_faceY, nint_faceY; 
  	int ntot_faceZ, nint_faceZ; 
	 	
	std::vector<double> spacing;
  	std::vector<double> area;
	
  	std::vector<double> x_cell;
  	std::vector<double> y_cell;
  	std::vector<double> z_cell;
	
  	std::vector<double> x_faceX;
  	std::vector<double> y_faceX;
  	std::vector<double> z_faceX;
  	std::vector<double> x_faceY;
  	std::vector<double> y_faceY;
  	std::vector<double> z_faceY;
  	std::vector<double> x_faceZ;
  	std::vector<double> y_faceZ;
  	std::vector<double> z_faceZ;
	
  	double volume;	
	
  	std::vector<int> index_bothsides_cell;
  	std::vector<int> index_bothsides_faceX;
  	std::vector<int> index_bothsides_faceY;
  	std::vector<int> index_bothsides_faceZ;
  	
	std::vector<int> index_bothsides_x_cell;
  	std::vector<int> index_bothsides_y_cell;
  	std::vector<int> index_bothsides_z_cell;
	std::vector<int> index_bothsides_x_faceX;
  	std::vector<int> index_bothsides_y_faceX;
  	std::vector<int> index_bothsides_z_faceX;
	std::vector<int> index_bothsides_x_faceY;
  	std::vector<int> index_bothsides_y_faceY;
  	std::vector<int> index_bothsides_z_faceY;
	std::vector<int> index_bothsides_x_faceZ;
  	std::vector<int> index_bothsides_y_faceZ;
  	std::vector<int> index_bothsides_z_faceZ;
	
  	std::vector<int> index_noside_cell;
  	std::vector<int> index_noside_faceX;
  	std::vector<int> index_noside_faceY;
  	std::vector<int> index_noside_faceZ;
	
	std::vector<int> index_noside_x_cell;
  	std::vector<int> index_noside_y_cell;
  	std::vector<int> index_noside_z_cell;
	std::vector<int> index_noside_x_faceX;
  	std::vector<int> index_noside_y_faceX;
  	std::vector<int> index_noside_z_faceX;
	std::vector<int> index_noside_x_faceY;
  	std::vector<int> index_noside_y_faceY;
  	std::vector<int> index_noside_z_faceY;
	std::vector<int> index_noside_x_faceZ;
  	std::vector<int> index_noside_y_faceZ;
  	std::vector<int> index_noside_z_faceZ;
	
  	grid_class( const std::vector<int>& _dim);
	
private:

};

grid_class::grid_class( const std::vector<int> & _dim)
 
{
	int point;

	dim = _dim;
	
	// Cell Field
	nx_cell =  dim[X_DIR] + 2;
	ny_cell =  dim[Y_DIR]>1 ? dim[Y_DIR] + 2 : 1;
	nz_cell =  dim[Z_DIR]>1 ? dim[Z_DIR] + 2 : 1;
	
	ntot_cell = nx_cell*ny_cell*nz_cell;
        nint_cell  = dim[X_DIR]*dim[Y_DIR]*dim[Z_DIR];
	
	// FaceX Field
	nx_faceX = nx_cell + 1;
	ny_faceX = ny_cell;
	nz_faceX = nz_cell;
	ntot_faceX = nx_faceX * ny_faceX * nz_faceX;
	nint_faceX = dim[X_DIR]*dim[Y_DIR]*dim[Z_DIR];
	
	// FaceY Field
	nx_faceY = nx_cell;
	ny_faceY = ny_cell + 1;
	nz_faceY = nz_cell;
	ntot_faceY = nx_faceY * ny_faceY * nz_faceY;
	nint_faceY = dim[X_DIR]*dim[Y_DIR]*dim[Z_DIR];

	// FaceZ Field
	nx_faceZ = nx_cell;
	ny_faceZ = ny_cell;
	nz_faceZ = nz_cell+1;
	ntot_faceZ = nx_faceZ * ny_faceZ * nz_faceZ;
	nint_faceZ = dim[X_DIR]*dim[Y_DIR]*dim[Z_DIR];
	
	// Staggered grid from 0 to PI in each direction
	// #point           0           1           2                       n          n+1
	// Grid:      |-----O-----|-----O-----|-----O-----|-----O-----|-----O-----|-----O-----|
	// #face      0           1           2                       n          n+1
	// coord.                 0          PI/n       2PI/n                     PI        
	
	// -----------------------------------------------------------------------------------------------------------------------------	
  	// Resizing vectors
  	// -----------------------------------------------------------------------------------------------------------------------------			
  	spacing.resize(3,0.0);			// spacing in each direction
  	area.resize(3,0.0);			// area in each direction
  
  	for( int i=0; i<3; ++i )
    		spacing[i] = (dim[i]==1) ? 1.0 : (LB-LA)/(dim[i]);
	
	area[X_DIR] = spacing[Y_DIR]*spacing[Z_DIR];
	area[Y_DIR] = spacing[X_DIR]*spacing[Z_DIR];
  	area[Z_DIR] = spacing[X_DIR]*spacing[Y_DIR];

  	volume = spacing[X_DIR]*spacing[Y_DIR]*spacing[Z_DIR];
	
	// -----------------------------------------------------------------------------------------------------------------------------	
	// Grid Coordinates
	// -----------------------------------------------------------------------------------------------------------------------------	
  	x_cell.resize(ntot_cell,0.0);		// xCoordinates - Cell Centered
  	y_cell.resize(ntot_cell,0.0);		// yCoordinates - Cell Centered
  	z_cell.resize(ntot_cell,0.0);		// zCoordinates - Cell Centered
  	
	x_faceX.resize(ntot_faceX,0.0);		// xCoordinates - faceX Centered
  	y_faceX.resize(ntot_faceX,0.0);		// yCoordinates - faceX Centered
 	z_faceX.resize(ntot_faceX,0.0);		// zCoordinates - faceX Centered

	x_faceY.resize(ntot_faceY,0.0);		// xCoordinates - faceY Centered
  	y_faceY.resize(ntot_faceY,0.0);		// yCoordinates - faceY Centered
 	z_faceY.resize(ntot_faceY,0.0);		// zCoordinates - faceY Centered

	x_faceZ.resize(ntot_faceZ,0.0);		// xCoordinates - faceY Centered
  	y_faceZ.resize(ntot_faceZ,0.0);		// yCoordinates - faceY Centered
 	z_faceZ.resize(ntot_faceZ,0.0);		// zCoordinates - faceY Centered
		
	// Cell Centered Coordinates
	point=0;	
  	for( int k=0; k<nz_cell; ++k )
    		for( int j=0; j<ny_cell; ++j)
      			for( int i=0; i<nx_cell; ++i )
			{
	  			x_cell[point] = double(i-0.5)*spacing[0] + LA;		// xCoordinates 				
	  			y_cell[point] = double(j-0.5)*spacing[1] + LA;		// yCoordinates 
	  			z_cell[point] = double(k-0.5)*spacing[2]+ LA;		// zCoordinates
	  			++point;
			}
			
	// FaceX Coordinates
	point=0;	
  	for( int k=0; k<nz_faceX; ++k )
    		for( int j=0; j<ny_faceX; ++j)
      			for( int i=0; i<nx_faceX; ++i )
			{
	  			x_faceX[point] = double(i-1.0)*spacing[0] + LA;		// xCoordinates 				
	  			y_faceX[point] = double(j-0.5)*spacing[1] + LA;		// yCoordinates 
	  			z_faceX[point] = double(k-0.5)*spacing[2]+ LA;		// zCoordinates
	  			++point;
			}		;   				

	// FaceY Coordinates
	point=0;	
  	for( int k=0; k<nz_faceY; ++k )
    		for( int j=0; j<ny_faceY; ++j)
      			for( int i=0; i<nx_faceY; ++i )
			{
	  			x_faceY[point] = double(i-0.5)*spacing[0] + LA;		// xCoordinates 				
	  			y_faceY[point] = double(j-1.0)*spacing[1] + LA;		// yCoordinates 
	  			z_faceY[point] = double(k-0.5)*spacing[2]+ LA;		// zCoordinates
	  			++point;
			}		;   				
			
	// FaceZ Coordinates
	point=0;	
  	for( int k=0; k<nz_faceZ; ++k )
    		for( int j=0; j<ny_faceZ; ++j)
      			for( int i=0; i<nx_faceZ; ++i )
			{
	  			x_faceZ[point] = double(i-0.50)*spacing[0] + LA;		// xCoordinates 				
	  			y_faceZ[point] = double(j-0.50)*spacing[1] + LA;		// yCoordinates 
	  			z_faceZ[point] = double(k-1.0)*spacing[2]+ LA;		// zCoordinates
	  			++point;
			}		;   				
						
	// skip BOTH sides - Cell Centered
  	for( int k=1; k<=nz_cell; ++k )
    		for( int j=1; j<=ny_cell; ++j)
      			for( int i=1; i<=nx_cell; ++i )
			{
	  			if (i>1 && i<=(nx_cell-1))
	    				if (j>1 && j<=(ny_cell-1))
	      					if (k>1 && k<=(nz_cell-1))
						{	
		  					index_bothsides_cell.push_back((i-1)+(j-1)*nx_cell+(k-1)*nx_cell*ny_cell); 
		 	 				index_bothsides_x_cell.push_back(i);
							 index_bothsides_y_cell.push_back(j); 
							 index_bothsides_z_cell.push_back(k);
						}									
			}

	// skip BOTH sides - FaceX Centered
  	for( int k=1; k<=nz_faceX; ++k )
    		for( int j=1; j<=ny_faceX; ++j)
      			for( int i=1; i<=nx_faceX; ++i )
			{
	  			if (i>1 && i<=(nx_faceX-1))
	    				if (j>1 && j<=(ny_faceX-1))
	      					if (k>1 && k<=(nz_faceX-1))
						{	
		  					index_bothsides_faceX.push_back((i-1)+(j-1)*nx_faceX+(k-1)*nx_faceX*ny_faceX); 
		 	 				index_bothsides_x_faceX.push_back(i);
							 index_bothsides_y_faceX.push_back(j); 
							 index_bothsides_z_faceX.push_back(k);
						}									
			}
			
	// skip BOTH sides - FaceY Centered
  	for( int k=1; k<=nz_faceY; ++k )
    		for( int j=1; j<=ny_faceY; ++j)
      			for( int i=1; i<=nx_faceY; ++i )
			{
	  			if (i>1 && i<=(nx_faceY-1))
	    				if (j>1 && j<=(ny_faceY-1))
	      					if (k>1 && k<=(nz_faceY-1))
						{	
		  					index_bothsides_faceY.push_back((i-1)+(j-1)*nx_faceY+(k-1)*nx_faceY*ny_faceY); 
		 	 				index_bothsides_x_faceY.push_back(i);
							 index_bothsides_y_faceY.push_back(j); 
							 index_bothsides_z_faceY.push_back(k);
						}									
			}			
	
	// skip BOTH sides - FaceZ Centered
  	for( int k=1; k<=nz_faceZ; ++k )
    		for( int j=1; j<=ny_faceZ; ++j)
      			for( int i=1; i<=nx_faceZ; ++i )
			{
	  			if (i>1 && i<=(nx_faceZ-1))
	    				if (j>1 && j<=(ny_faceZ-1))
	      					if (k>1 && k<=(nz_faceZ-1))
						{	
		  					index_bothsides_faceZ.push_back((i-1)+(j-1)*nx_faceZ+(k-1)*nx_faceZ*ny_faceZ); 
		 	 				index_bothsides_x_faceZ.push_back(i);
							 index_bothsides_y_faceZ.push_back(j); 
							 index_bothsides_z_faceZ.push_back(k);
						}									
			}	
										
	// no skipping - Cell Centered
  	for( int k=1; k<=nz_cell; ++k )
    		for( int j=1; j<=ny_cell; ++j)
      			for( int i=1; i<=nx_cell; ++i )  	
			{
	 	 		index_noside_cell.push_back((i-1)+(j-1)*nx_cell+(k-1)*nx_cell*ny_cell); 
	  			index_noside_x_cell.push_back(i); 
				index_noside_y_cell.push_back(j); 
				index_noside_z_cell.push_back(k);
			}
			
	// no skipping - FaceX Centered
  	for( int k=1; k<=nz_faceX; ++k )
    		for( int j=1; j<=ny_faceX; ++j)
      			for( int i=1; i<=nx_faceX; ++i )  	
			{
	 	 		index_noside_faceX.push_back((i-1)+(j-1)*nx_faceX+(k-1)*nx_faceX*ny_faceX); 
	  			index_noside_x_faceX.push_back(i); 
				index_noside_y_faceX.push_back(j); 
				index_noside_z_faceX.push_back(k);
			}
			
	// no skipping - FaceY Centered
  	for( int k=1; k<=nz_faceY; ++k )
    		for( int j=1; j<=ny_faceY; ++j)
      			for( int i=1; i<=nx_faceY; ++i )  	
			{
	 	 		index_noside_faceY.push_back((i-1)+(j-1)*nx_faceY+(k-1)*nx_faceY*ny_faceY); 
	  			index_noside_x_faceY.push_back(i); 
				index_noside_y_faceY.push_back(j); 
				index_noside_z_faceY.push_back(k);
			}	
			
	// no skipping - FaceZ Centered
  	for( int k=1; k<=nz_faceZ; ++k )
    		for( int j=1; j<=ny_faceZ; ++j)
      			for( int i=1; i<=nx_faceZ; ++i )  	
			{
	 	 		index_noside_faceZ.push_back((i-1)+(j-1)*nx_faceZ+(k-1)*nx_faceZ*ny_faceZ); 
	  			index_noside_x_faceZ.push_back(i); 
				index_noside_y_faceZ.push_back(j); 
				index_noside_z_faceZ.push_back(k);
			}			
			
	// Information on Video
	//cout << "Number of grid points - Cells: " << nx_cell		<< " x " << ny_cell 	<< " x " << nz_cell << " = " << ntot_cell << endl;
	//cout << "Number of grid points - FaceX: " << nx_faceX << " x " << ny_faceX 	<< " x " << nz_faceX << " = " << ntot_faceX << endl;
	//cout << "Number of grid points - FaceY: " << nx_faceY << " x " << ny_faceY 	<< " x " << nz_faceY << " = " << ntot_faceY << endl;
	//cout << "Number of grid points - FaceZ: " << nx_faceZ << " x " << ny_faceZ 	<< " x " << nz_faceZ << " = " << ntot_faceZ<< endl;
}
