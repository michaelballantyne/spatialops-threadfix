#include <FVStaggeredSpatialOps.h>

// *****************************************************************************************************************************	
// Check Equality
// *****************************************************************************************************************************	
void compare(const double fAnalytical, const double fNumerical, double &mean_rel_err, double &max_rel_err)
{
  double rel_err;
  const double difference = std::fabs(fAnalytical-fNumerical);
  const double minimum = min(fabs(fAnalytical),fabs(fNumerical));
	
  // In this case the relative error could be wrong
  if (difference<=1.e-16 && minimum<=1.e-16) rel_err = 1.e-17;
  // This is usual definition of relative error
  else rel_err = (minimum == 0) ? 1.e-17 : difference / ( minimum + TOL );
   	
  mean_rel_err += rel_err;
  if (rel_err > max_rel_err)
    max_rel_err = rel_err;
}

template< typename Dir >
void indices_for_SideFields(grid_class &grid, std::vector<int> &index)
{
  typedef SpatialOps::FVStaggeredUniform::DefaultSideGhosting<Dir>  GhostTraits;

  const int nXMinus = GhostTraits::template get<SpatialOps::XDIR,SpatialOps::SideMinus>();
  const int nYMinus = GhostTraits::template get<SpatialOps::YDIR,SpatialOps::SideMinus>();
  const int nZMinus = GhostTraits::template get<SpatialOps::ZDIR,SpatialOps::SideMinus>();
  const int nXPlus  = GhostTraits::template get<SpatialOps::XDIR,SpatialOps::SidePlus>();
  const int nYPlus  = GhostTraits::template get<SpatialOps::YDIR,SpatialOps::SidePlus>();
  const int nZPlus  = GhostTraits::template get<SpatialOps::ZDIR,SpatialOps::SidePlus>();
		
  const int nX = grid.nx + (nXMinus-1) + (nXPlus-1);
  const int nY = grid.ny + (nYMinus-1) + (nYPlus-1);
  const int nZ = grid.nz + (nZMinus-1) + (nZPlus-1);
		
  // skip BOTH sides
  const int klo=nZMinus;  const int khi = nZ-nZPlus;
  const int jlo=nYMinus;  const int jhi = nY-nYPlus;
  const int ilo=nXMinus;  const int ihi = nX-nXPlus;
  for( int k=klo; k<khi; ++k )
    for( int j=jlo; j<jhi; ++j )
      for( int i=ilo; i<ihi; ++i )
	index.push_back( i + j*nX+k*nX*nY );
}


void check_equality_C2F(grid_class &grid, double *analytical_int, XSideField &numerical, bool &test, double &mean_rel_err, double &max_rel_err)
{
  // ----------------------------------------------------------------------------------------------------------------------------------
  //	Error checking
  //-----------------------------------------------------------------------------------------------------------------------------------
  mean_rel_err = 0.; 
  max_rel_err = 0.; 
	
  // skip BOTH sides
  std::vector<int> index;

  indices_for_SideFields<SpatialOps::XDIR>(grid, index);
			
  for(unsigned int i=0; i<grid.index_bothsides.size(); i++)
    {    
      int point = grid.index_bothsides[i];
      int point_interpolation = index[i];
      compare( analytical_int[point], numerical.get_ptr()[point_interpolation], mean_rel_err, max_rel_err);
      // cout << endl << analytical_int[point] << " " << numerical.get_ptr()[point_interpolation] ;
    }	
	
  //for(unsigned int i=0; i<grid.nx*grid.ny*grid.nz; i++)
  //	cout << endl << analytical_int[i] << " " << numerical.get_ptr()[i] ;
		
  mean_rel_err /= double(grid.index_bothsides.size());

  std::cout << "  Mean Relative Error: " << mean_rel_err  << "  ";
  std::cout << "  Max Relative Error: "  << max_rel_err   << "  ";
	  	
  if (max_rel_err > MAX_MAXRELERR || mean_rel_err > MAX_MEANRELERR)
    test = false;
  else test = true;
}

void check_equality_C2F(grid_class &grid, double *analytical_int, YSideField &numerical, bool &test, double &mean_rel_err, double &max_rel_err)
{
  // ----------------------------------------------------------------------------------------------------------------------------------
  //	Error checking
  //-----------------------------------------------------------------------------------------------------------------------------------
  mean_rel_err = 0.; 
  max_rel_err = 0.; 
	
  // skip BOTH sides
  std::vector<int> index;

  indices_for_SideFields<SpatialOps::YDIR>(grid, index);
			
  for(unsigned int i=0; i<grid.index_bothsides.size(); i++)
    {    
      int point = grid.index_bothsides[i];
      int point_interpolation = index[i];
      compare( analytical_int[point], numerical.get_ptr()[point_interpolation], mean_rel_err, max_rel_err);
      // cout << endl << analytical_int[point] << " " << numerical.get_ptr()[point_interpolation] ;
    }	
	
  //for(unsigned int i=0; i<grid.nx*grid.ny*grid.nz; i++)
  //	cout << endl << analytical_int[i] << " " << numerical.get_ptr()[i] ;
		
  mean_rel_err /= double(grid.index_bothsides.size());

  std::cout << "  Mean Relative Error: " << mean_rel_err  << "  ";
  std::cout << "  Max Relative Error: "  << max_rel_err   << "  ";
	  	
  if (max_rel_err > MAX_MAXRELERR || mean_rel_err > MAX_MEANRELERR)
    test = false;
  else test = true;
}

void check_equality_C2F(grid_class &grid, double *analytical_int, ZSideField &numerical, bool &test, double &mean_rel_err, double &max_rel_err)
{
  // ----------------------------------------------------------------------------------------------------------------------------------
  //	Error checking
  //-----------------------------------------------------------------------------------------------------------------------------------
  mean_rel_err = 0.; 
  max_rel_err = 0.; 
	
  // skip BOTH sides
  std::vector<int> index;

  indices_for_SideFields<SpatialOps::ZDIR>(grid, index);
			
  for(unsigned int i=0; i<grid.index_bothsides.size(); i++)
    {    
      int point = grid.index_bothsides[i];
      int point_interpolation = index[i];
      compare( analytical_int[point], numerical.get_ptr()[point_interpolation], mean_rel_err, max_rel_err);
      //cout << endl << analytical_int[point] << " " << numerical.get_ptr()[point_interpolation] ;
    }	
	
  //for(unsigned int i=0; i<grid.nx*grid.ny*grid.nz; i++)
  //	cout << endl << analytical_int[i] << " " << numerical.get_ptr()[i] ;
		
  mean_rel_err /= double(grid.index_bothsides.size());

  std::cout << "  Mean Relative Error: " << mean_rel_err  << "  ";
  std::cout << "  Max Relative Error: "  << max_rel_err   << "  ";
	  	
  if (max_rel_err > MAX_MAXRELERR || mean_rel_err > MAX_MEANRELERR)
    test = false;
  else test = true;
}

void check_equality_F2C(grid_class &grid, double *analytical_int, CellField &numerical, bool &test, double &mean_rel_err, double &max_rel_err)
{
  // ----------------------------------------------------------------------------------------------------------------------------------
  //	Error checking
  //-----------------------------------------------------------------------------------------------------------------------------------
  mean_rel_err = 0.; 
  max_rel_err = 0.; 
	   	   	
  for(unsigned int i=0; i<grid.index_bothsides.size(); i++)
    {    
      int point = grid.index_bothsides[i];
      compare( analytical_int[point], numerical.get_ptr()[point], mean_rel_err, max_rel_err);
      //cout << endl << analytical_int[point] << " " << numerical[point] ;
    }	
	
  mean_rel_err /= double(grid.index_bothsides.size());

  std::cout << "  Mean Relative Error: " << mean_rel_err  << "  ";
  std::cout << "  Max Relative Error: "  << max_rel_err   << "  ";
	  	
  if (max_rel_err > MAX_MAXRELERR || mean_rel_err > MAX_MEANRELERR)
    test = false;
  else test = true;
}

bool compare_scaling(string kind, Epetra_CrsMatrix &SpOp,  Epetra_CrsMatrix &SpOpOld, double *field)
{
  bool ok = true;
	
  for( int i=0; i<SpOp.NumGlobalRows(); ++i )
    {
      double *vals;    	int *ixs;    	int nentries;
      double *valsOld;    int *ixsOld;    	int nentriesOld;
    	
      SpOp.ExtractMyRowView( i, nentries, vals, ixs );
      SpOpOld.ExtractMyRowView( i, nentriesOld, valsOld, ixsOld );
    		
      assert( nentries == nentriesOld );
		
      if (kind == "LEFT")
	for( int k=0; k<nentries; ++k )
	  {
	    if( *vals++ != (*valsOld++)*field[i] ) ok=false;
	  }
			
      if (kind == "RIGHT")
	for( int k=0; k<nentries; ++k )
	  {
	    if( *vals++ != (*valsOld++)*field[*ixs++] ) ok=false;
	  }
    }
  	
  return ok;
}

// Conversion from Int to string
std::string IntToString(int number)
{
  std::ostringstream oss;

  // Works just like cout
  oss<< number;

  // Return the underlying string
  return oss.str();
}
