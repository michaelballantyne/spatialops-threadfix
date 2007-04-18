#ifndef UT_FVStaggeredSpatialOps_h
#define UT_FVStaggeredSpatialOps_h

#include <SpatialOperator.h>
#include <SpatialField.h>
#include <SpatialOpsDefs.h>

namespace SpatialOps{
namespace FVStaggeredUniform{

  struct Side{};  // provides type information for face fields
  struct Cell{};  // provides type information for cell fields


  /** Policy for a field with no ghost cells */
  struct NoGhosting
  {
    template<typename Dir, typename SideType >
    static int get(){return 0;}
  };


  // Policy defining ghosting for a cell field
  struct DefaultCellGhosting
  {
    template<typename Dir, typename SideType >
    static int get(){return 1;}
  };


  template< typename Dir >
  struct DefaultSideGhosting
  {
    // note: specialized template member functions exist!
    template<typename Direction, typename SideType>
    static int get(){return 1;}
  };


  //==================================================================

  template< class Location,
	    class GhostPolicy >
  struct FieldTraits
  {
    typedef GhostPolicy GhostTraits;      // required by SpatialField policy for FieldTraits
    typedef Location    StorageLocation;  // required by SpatialField policy for FieldTraits
  };

  //==================================================================
  

  class OpInfo
  {
  public:

    OpInfo( const std::vector<int> & extent,
	    const std::vector<int> & ngs,
	    const std::vector<int> & ngd );

    ~OpInfo();

    void get_bounds( const int irow,
		     const bool staggerLeft,	
		     const int nentries,
		     bool* const inBoundsLo,
		     bool* const inBoundsHi,
		     bool* const inBounds,
		     int & icol ) const;

    inline const std::vector<int>& extent () const{ return extent_; }
    inline const std::vector<int>& ng_src () const{ return ngs_; }
    inline const std::vector<int>& ng_dest() const{ return ngd_; }
    inline int get_nxd() const{return nxd_;}
    inline int get_nyd() const{return nyd_;}
    inline int get_nzd() const{return nzd_;}
    inline int get_nxs() const{return nxs_;}
    inline int get_nys() const{return nys_;}
    inline int get_nzs() const{return nzs_;}

    inline int nrows() const{ return nrows_; }
    inline int ncols() const{ return ncols_;  }

  private:
    const std::vector<int> extent_, ngs_, ngd_;
    const int nxd_, nyd_, nzd_;
    const int nxs_, nys_, nzs_;
    const int nrows_, ncols_;
  };


  //==================================================================


  template<typename GhostPolicy>
  inline int entrycount( const std::vector<int> & dimExtent );


  //==================================================================


  template< typename Direction,
	    typename Location,        // location of source field => location of operator.
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  class LinearInterpolantAssembler
  {
  public:

    typedef Location  OpLocation;
    typedef Direction DirType;

    static int num_nonzeros(){ return 2; }

    LinearInterpolantAssembler( const std::vector<int> & dimExtent );
    
    ~LinearInterpolantAssembler(){};

    int get_ncols() const;
    int get_nrows() const;

    inline const std::vector<int>& get_extent() const{ return opInfo_.extent(); }

    void get_row_entries( const int irow,
			  std::vector<double> & vals,
			  std::vector<int> & ixs );

  private:

    const OpInfo opInfo_;
    bool inBoundsLo_[3], inBoundsHi_[3], inBounds_[3];

  };


  //==================================================================


  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  class DivergenceAssembler
  {
  public:

    typedef Direction DirType;
    typedef Location OpLocation;

    static int num_nonzeros(){ return 2; }

    DivergenceAssembler( const std::vector<int> & dimExtent,
			 const std::vector<double> & cellFaceArea,
			 const double cellVolume );

    ~DivergenceAssembler(){}

    int get_ncols() const;
    int get_nrows() const;

    inline const std::vector<int>& get_extent() const{ return opInfo_.extent(); }

    void get_row_entries( const int irow,
			  std::vector<double> & vals,
			  std::vector<int> & ixs );

  private:

    const OpInfo opInfo_;
    const std::vector<double> faceArea_;
    const double cellVol_;
    bool inBoundsLo_[3], inBoundsHi_[3], inBounds_[3];

  };


  //==================================================================


  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  class GradientAssembler
  {
  public:

    typedef Direction DirType;
    typedef Location OpLocation;

    static int num_nonzeros(){ return 2; }

    GradientAssembler( const std::vector<double> & meshSpacing,
		       const std::vector<int> & dimExtent );

    ~GradientAssembler(){}

    int get_ncols() const;
    int get_nrows() const;

    inline const std::vector<int>& get_extent() const{ return opInfo_.extent(); }

    void get_row_entries( const int irow,
			  std::vector<double> & vals,
			  std::vector<int> & ixs );

  private:

    const OpInfo opInfo_;
    const std::vector<double> spacing_;
    bool inBoundsLo_[3], inBoundsHi_[3], inBounds_[3];

  };


  //==================================================================


  template< int NumNonZero,
	    typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  class ScratchAssembler
  {
  public:

    typedef Direction DirType;
    typedef Location  OpLocation;

    static int num_nonzeros(){ return NumNonZero; }

    ScratchAssembler( const std::vector<int> & dimExtent );

    ~ScratchAssembler(){}

    int get_ncols() const;
    int get_nrows() const;

    inline const std::vector<int>& get_extent() const{ return opInfo_.extent(); }

    void get_row_entries( const int irow,
			  std::vector<double> & vals,
			  std::vector<int> & ixs );
    
  private:

    const OpInfo opInfo_;
    bool inBoundsLo_[3], inBoundsHi_[3], inBounds_[3];

  };
  //------------------------------------------------------------------


  //==================================================================


  // helper functions...
  template<typename Location>
  bool stagger_left();


  template<typename Direction>
  void interp_bounds_helper( const bool bnds[3], const bool loBnds[3], const bool hiBnds[3],
			     const OpInfo& opInfo,
			     double & val, int & ifac );

  template<typename Direction>
  void grad_bounds_helper( const bool bnds[3], const bool loBnds[3], const bool hiBnds[3],
			   const std::vector<double> & spacing,
			   const OpInfo& opInfo,
			   double & val, int & ifac, int & shift );

  template<typename Direction>
  void div_bounds_helper( const bool bnds[3], const bool loBnds[3], const bool hiBnds[3],
			  const std::vector<double> & faceArea,
			  const double cellVol,
			  const OpInfo& opInfo,
			  double & val, int & ifac, int& shift );

  template<typename Direction>
  void scratch_bounds_helper( const bool loBnds[3],
			      const OpInfo& opInfo,
			      int & skip,
			      int & shift );

  //==================================================================







  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //  Implementations
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







  //==================================================================


  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  LinearInterpolantAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  LinearInterpolantAssembler( const std::vector<int> & dimExtent )
    : opInfo_( dimExtent,
	       ghost_vec<SrcGhostPolicy>(),
	       ghost_vec<DestGhostPolicy>() )
  {
  }
  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  void
  LinearInterpolantAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_row_entries( const int irow,
		   std::vector<double> & vals,
		   std::vector<int> & ixs )
  {
    const bool staggerLeft = stagger_left<Location>();

    int icol;
    opInfo_.get_bounds( irow, staggerLeft, 2, inBoundsLo_, inBoundsHi_, inBounds_, icol );
    double val = 0.5;
    int ifac = (staggerLeft ? -1 : 1);

    interp_bounds_helper<Direction>( inBounds_, inBoundsLo_, inBoundsHi_, opInfo_, val, ifac );
    const int shift = ifac;

    vals.push_back( val );  ixs.push_back( icol       );
    vals.push_back( val );  ixs.push_back( icol+shift );
  }
  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  int
  LinearInterpolantAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_ncols() const
  {
    return entrycount<SrcGhostPolicy>( opInfo_.extent() );
  }
  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  int
  LinearInterpolantAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_nrows() const
  {
    return entrycount<DestGhostPolicy>( opInfo_.extent() );
  }
  //------------------------------------------------------------------


  //==================================================================


  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  DivergenceAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  DivergenceAssembler( const std::vector<int> & dimExtent,
		       const std::vector<double> & cellFaceArea,
		       const double cellVolume )
    : opInfo_( dimExtent,
	       ghost_vec<SrcGhostPolicy>(),
	       ghost_vec<DestGhostPolicy>() ),
      faceArea_( cellFaceArea ),
      cellVol_ ( cellVolume   )
  {
  }
  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  void
  DivergenceAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_row_entries( const int irow,
		   std::vector<double> & vals,
		   std::vector<int> & ixs )
  {
    const bool staggerLeft = stagger_left<Location>();

    int icol;
    opInfo_.get_bounds( irow, staggerLeft, 2, inBoundsLo_, inBoundsHi_, inBounds_, icol );
    int ifac = (staggerLeft ? -1 : 1);

    int shift;
    double fac;
    div_bounds_helper<Direction>( inBounds_, inBoundsLo_, inBoundsHi_, faceArea_, cellVol_, opInfo_,
				  fac, ifac, shift );

    vals.push_back( -ifac*fac );  ixs.push_back( icol       );
    vals.push_back(  ifac*fac );  ixs.push_back( icol+shift );
  }
  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  int
  DivergenceAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_ncols() const
  {
    return entrycount<SrcGhostPolicy>( opInfo_.extent() );
  }
  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  int
  DivergenceAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_nrows() const
  {
    return entrycount<DestGhostPolicy>( opInfo_.extent() );
  }
  //------------------------------------------------------------------
  

  //==================================================================


  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  GradientAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  GradientAssembler( const std::vector<double> & meshSpacing,
		     const std::vector<int> & dimExtent )
    : opInfo_( dimExtent,
	       ghost_vec<SrcGhostPolicy>(),
	       ghost_vec<DestGhostPolicy>() ),
      spacing_( meshSpacing )
  {
  }
  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  void
  GradientAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_row_entries( const int irow,
		   std::vector<double> & vals,
		   std::vector<int> & ixs )
  {
    const bool staggerLeft = stagger_left<Location>();

    int icol;
    opInfo_.get_bounds( irow, staggerLeft, 2, inBoundsLo_, inBoundsHi_, inBounds_, icol );
    int ifac = (staggerLeft ? -1 : 1);

    double fac = 0;
    int shift;
    grad_bounds_helper<Direction>( inBounds_, inBoundsLo_, inBoundsHi_, spacing_, opInfo_, fac, ifac, shift );

    vals.push_back( -ifac*fac );  ixs.push_back( icol       );
    vals.push_back(  ifac*fac );  ixs.push_back( icol+shift );
  }
  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  int
  GradientAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_ncols() const
  {
    return entrycount<SrcGhostPolicy>( opInfo_.extent() );
  }
  //------------------------------------------------------------------
  template< typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  int
  GradientAssembler<Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_nrows() const
  {
    return entrycount<DestGhostPolicy>( opInfo_.extent() );
  }
  //------------------------------------------------------------------
  

  //==================================================================


  //------------------------------------------------------------------
  template< int NumNonZero,
	    typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  ScratchAssembler<NumNonZero,Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  ScratchAssembler::ScratchAssembler( const std::vector<int> & dimExtent )
    : opInfo_( dimExtent,
	       ghost_vec<SrcGhostPolicy >(),
	       ghost_vec<DestGhostPolicy>() )
  {
  }
  //------------------------------------------------------------------
  template< int NumNonZero,
	    typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  void
  ScratchAssembler<NumNonZero,Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_row_entries( const int irow,
		   std::vector<double> & vals,
		   std::vector<int> & ixs )
  {
    static const bool staggerLeft = stagger_left<Location>();

    const double val = 0.0;
    int icol;
    opInfo_.get_bounds( irow, staggerLeft, NumNonZero, inBoundsLo_, inBoundsHi_, inBounds_, icol );

    const int ncols = opInfo_.ncols();

    int nShift = -( int(std::ceil( float(NumNonZero)/2.0 )) - 1 );
    static const int nOffDiagP = NumNonZero/2;
    if( NumNonZero%2 == 0 ){
      if( staggerLeft )  --nShift;
      else               ++nShift;
    }

    int skip = 1;
    scratch_bounds_helper<Direction>( inBoundsLo_, opInfo_, skip, nShift );

    int ix = nShift*skip;
    while( ix+icol <0 ) ix+=skip;
    while( ix+icol + nOffDiagP*skip > ncols ) ix-=skip;
    for( int i=0; i<NumNonZero; ++i ){
      const int ipos = icol + ix;
      if( ipos >= 0 && ipos < ncols ){
	vals.push_back( val );  ixs.push_back( ipos );
      }
      ix += skip;
    }
  }
  //------------------------------------------------------------------
  template< int NumNonZero,
	    typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  int
  ScratchAssembler<NumNonZero,Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_ncols() const
  {
    return entrycount<SrcGhostPolicy>( opInfo_.extent() );
  }
  //------------------------------------------------------------------
  template< int NumNonZero,
	    typename Direction,
	    typename Location,
	    class SrcGhostPolicy,
	    class DestGhostPolicy >
  int
  ScratchAssembler<NumNonZero,Direction,Location,SrcGhostPolicy,DestGhostPolicy>::
  get_nrows() const
  {
    return entrycount<DestGhostPolicy>( opInfo_.extent() );
  }
  //------------------------------------------------------------------

  //==================================================================

  template<typename GhostPolicy>
  int entrycount( const std::vector<int> & dims )
  {
    int nn;
    int nrows=1;

    nn = dims[0] + GhostPolicy::template get<XDIR,SideMinus>() + GhostPolicy::template get<XDIR,SidePlus>();
    nrows *= nn;

    if( dims[1] > 1 ){
      nn = dims[1] + GhostPolicy::template get<YDIR,SideMinus>() + GhostPolicy::template get<YDIR,SidePlus>();
      nrows *= nn;
    }
    if( dims[2] > 1 ){
      nn = dims[2] + GhostPolicy::template get<ZDIR,SideMinus>() + GhostPolicy::template get<ZDIR,SidePlus>();
      nrows *= nn;
    }
    return nrows;
  }

  //==================================================================


  template<>
  template<>
  inline int DefaultSideGhosting<XDIR>::get<XDIR,SidePlus>(){return 2;}

  template<>
  template<>
  inline int DefaultSideGhosting<YDIR>::get<YDIR,SidePlus>(){return 2;}

  template<>
  template<>
  inline int DefaultSideGhosting<ZDIR>::get<ZDIR,SidePlus>(){return 2;}


  //====================================================================


} // namespace FVStaggeredUniform
} // namespace SpatialOps

#endif
