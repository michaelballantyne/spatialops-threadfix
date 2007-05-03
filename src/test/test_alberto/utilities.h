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


void check_equality_C2F(const grid_class &grid, XSideField &analytical, XSideField &numerical, bool &test, double &mean_rel_err, double &max_rel_err)
{
  // ----------------------------------------------------------------------------------------------------------------------------------
  //	Error checking
  //-----------------------------------------------------------------------------------------------------------------------------------
  mean_rel_err = 0.; 
  max_rel_err = 0.; 

  for(unsigned int i=0; i<grid.index_bothsides_faceX.size(); i++)
    {    
      int point = grid.index_bothsides_faceX[i];
      compare( analytical.begin()[point], numerical.begin()[point], mean_rel_err, max_rel_err);
    }	
			
  mean_rel_err /= double(grid.index_bothsides_faceX.size());

  std::cout << "  Mean Relative Error: " << mean_rel_err  << "  ";
  std::cout << "  Max Relative Error: "  << max_rel_err   << "  ";
	  	
  if (max_rel_err > MAX_MAXRELERR || mean_rel_err > MAX_MEANRELERR)
    test = false;
  else test = true;
}

void check_equality_C2F(const grid_class &grid, YSideField &analytical, YSideField &numerical, bool &test, double &mean_rel_err, double &max_rel_err)
{
  // ----------------------------------------------------------------------------------------------------------------------------------
  //	Error checking
  //-----------------------------------------------------------------------------------------------------------------------------------
  mean_rel_err = 0.; 
  max_rel_err = 0.; 
				
  for(unsigned int i=0; i<grid.index_bothsides_faceY.size(); i++)
    {    
      int point = grid.index_bothsides_faceY[i];
      compare( analytical.begin()[point], numerical.begin()[point], mean_rel_err, max_rel_err);
    }	
			
  mean_rel_err /= double(grid.index_bothsides_faceY.size());

  std::cout << "  Mean Relative Error: " << mean_rel_err  << "  ";
  std::cout << "  Max Relative Error: "  << max_rel_err   << "  ";
	  	
  if (max_rel_err > MAX_MAXRELERR || mean_rel_err > MAX_MEANRELERR)
    test = false;
  else test = true;
}

void check_equality_C2F(const grid_class &grid, ZSideField &analytical, ZSideField &numerical, bool &test, double &mean_rel_err, double &max_rel_err)
{
  // ----------------------------------------------------------------------------------------------------------------------------------
  //	Error checking
  //-----------------------------------------------------------------------------------------------------------------------------------
  mean_rel_err = 0.; 
  max_rel_err = 0.; 
	
  for(unsigned int i=0; i<grid.index_bothsides_faceZ.size(); i++)
    {    
      int point = grid.index_bothsides_faceZ[i];
      compare( analytical.begin()[point], numerical.begin()[point], mean_rel_err, max_rel_err);
    }	
			
  mean_rel_err /= double(grid.index_bothsides_faceZ.size());

  std::cout << "  Mean Relative Error: " << mean_rel_err  << "  ";
  std::cout << "  Max Relative Error: "  << max_rel_err   << "  ";
	  	
  if (max_rel_err > MAX_MAXRELERR || mean_rel_err > MAX_MEANRELERR)
    test = false;
  else test = true;
}

void check_equality_F2C(const grid_class &grid, CellField &analytical, CellField &numerical, bool &test, double &mean_rel_err, double &max_rel_err)
{
  // ----------------------------------------------------------------------------------------------------------------------------------
  //	Error checking
  //-----------------------------------------------------------------------------------------------------------------------------------
  mean_rel_err = 0.; 
  max_rel_err = 0.; 
	   	   	
  for(unsigned int i=0; i<grid.index_bothsides_cell.size(); i++)
    {    
      int point = grid.index_bothsides_cell[i];
      compare( analytical[point], numerical[point], mean_rel_err, max_rel_err);
    }	
    	
  mean_rel_err /= double(grid.index_bothsides_cell.size());

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
