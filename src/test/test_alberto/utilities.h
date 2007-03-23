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

void indeces_for_SideFields(char C, grid_class &grid, std::vector<int> &index)
{
  int nX, nY, nZ;
  int nXMinus = 0;
  int  nXPlus = 0;
  int  nYMinus =0;
  int  nYPlus =0;
  int  nZMinus =0;
  int nZPlus = 0;
	
  nX = grid.nx;
  nY = grid.ny;
  nZ = grid.nz;
	
  if (C=='X')
    {
      nXMinus = XSideFieldTraits::GhostTraits::get<XDIR,SideMinus>();
      nYMinus = XSideFieldTraits::GhostTraits::get<YDIR,SidePlus>();
      nZMinus = XSideFieldTraits::GhostTraits::get<ZDIR,SidePlus>();
      nXPlus   = XSideFieldTraits::GhostTraits::get<XDIR,SidePlus>();
      nYPlus   = XSideFieldTraits::GhostTraits::get<YDIR,SidePlus>();
      nZPlus   = XSideFieldTraits::GhostTraits::get<ZDIR,SidePlus>();
    }

  else if (C=='Y')
    {
      nXMinus = YSideFieldTraits::GhostTraits::get<XDIR,SideMinus>();
      nYMinus = YSideFieldTraits::GhostTraits::get<YDIR,SidePlus>();
      nZMinus = YSideFieldTraits::GhostTraits::get<ZDIR,SidePlus>();
      nXPlus   = YSideFieldTraits::GhostTraits::get<XDIR,SidePlus>();
      nYPlus   = YSideFieldTraits::GhostTraits::get<YDIR,SidePlus>();
      nZPlus   = YSideFieldTraits::GhostTraits::get<ZDIR,SidePlus>();
    }
	
  else if (C=='Z')
    {
      nXMinus = ZSideFieldTraits::GhostTraits::get<XDIR,SideMinus>();
      nYMinus = ZSideFieldTraits::GhostTraits::get<YDIR,SidePlus>();
      nZMinus = ZSideFieldTraits::GhostTraits::get<ZDIR,SidePlus>();
      nXPlus   = ZSideFieldTraits::GhostTraits::get<XDIR,SidePlus>();
      nYPlus   = ZSideFieldTraits::GhostTraits::get<YDIR,SidePlus>();
      nZPlus   = ZSideFieldTraits::GhostTraits::get<ZDIR,SidePlus>();
    }

		
  nX += (nXMinus-1);
  nY += (nYMinus-1);
  nZ += (nZMinus-1);
		
  // skip BOTH sides
  for( int k=1; k<=nZ; ++k )
    for( int j=1; j<=nY; ++j )
      for( int i=1; i<=nX; ++i )
	{
	  if (i>nXMinus && i<=(nX-nXPlus))
	    if (j>nYMinus && j<=(nY-nYPlus))
	      if (k>nZMinus && k<=(nZ-nZPlus))
		index.push_back( (i-1)+(j-1)*nX+(k-1)*nX*nY ); 
	}  
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

  indeces_for_SideFields('X', grid, index);
			
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

  indeces_for_SideFields('Y', grid, index);
			
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

  indeces_for_SideFields('Z', grid, index);
			
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
