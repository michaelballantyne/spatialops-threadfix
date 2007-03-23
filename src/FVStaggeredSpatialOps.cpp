#include <FVStaggeredSpatialOps.h>

#include <sstream>
#include <stdexcept>
#include <numeric>
#include <algorithm>

using std::vector;

namespace SpatialOps{
namespace FVStaggeredUniform{

  //==================================================================


  OpInfo::OpInfo( const std::vector<int> & extent,
		  const std::vector<int> & ngs,
		  const std::vector<int> & ngd )
    : extent_( extent ),
      ngs_( ngs ),
      ngd_( ngd ),

      nxd_( extent_[0]+ngd_[0]+ngd_[1] ),
      nyd_( extent_[1] > 1 ? extent_[1]+ngd_[2]+ngd_[3] : 1 ),
      nzd_( extent_[2] > 1 ? extent_[2]+ngd_[4]+ngd_[5] : 1 ),

      nxs_( extent_[0]+ngs_[0]+ngs_[1] ),
      nys_( extent_[1] > 1 ? extent_[1]+ngs_[2]+ngs_[3] : 1 ),
      nzs_( extent_[2] > 1 ? extent_[2]+ngs_[4]+ngs_[5] : 1 ),

      nrows_( nxd_*nyd_*nzd_ ),
      ncols_( nxs_*nys_*nzs_ )
  {
  }
  //------------------------------------------------------------------
  OpInfo::~OpInfo()
  {
  }
  //------------------------------------------------------------------
  void
  OpInfo::get_bounds( const int irow,
		      const bool staggerLeft,
		      const int nentries,
		      bool* const inBoundsLo,
		      bool* const inBoundsHi,
		      bool* const inBounds,
		      int & icol ) const
  {
    // get the (ijk) index for the dest array
    const int idest = irow%nxd_;
    const int jdest = irow/nxd_ % nyd_;
    const int kdest = irow/(nxd_*nyd_);

    // get the (ijk) index for the src array
    int isrc = idest - ngd_[0] + ngs_[0];
    int jsrc = jdest - ngd_[2] + ngs_[2];
    int ksrc = kdest - ngd_[4] + ngs_[4];

    int nShift = -(int(std::ceil( float(nentries)/2.0 )) - 1);
    if( nentries%2 == 0 )
      if( !staggerLeft ) nShift--;
    for(int i=0; i<3; ++i){ inBoundsLo[i]=inBoundsHi[i]=true; }

    if( isrc <     1+nShift ){ inBoundsLo[0]=false; isrc=0;      }
    if( isrc >= nxs_+nShift ){ inBoundsHi[0]=false; isrc=nxs_-1; }
    if( extent_[1] > 1 ){
      if( jsrc <     1+nShift ){ inBoundsLo[1]=false; jsrc=0;      }
      if( jsrc >= nys_+nShift ){ inBoundsHi[1]=false; jsrc=nys_-1; }
    }
    if( extent_[2] > 1 ){
      if( ksrc <     1+nShift ){ inBoundsLo[2]=false; ksrc=0;      }
      if( ksrc >= nzs_+nShift ){ inBoundsHi[2]=false; ksrc=nzs_-1; }
    }

    icol = ksrc*(nxs_*nys_) + jsrc*nxs_ + isrc;

    for( int i=0; i<3; ++i ){
      inBounds[i] = ( inBoundsLo[i] && inBoundsHi[i] );
    }
  }
  //------------------------------------------------------------------


  //==================================================================


  template<>
  bool stagger_left<Cell>(){ return true; }

  template<>
  bool stagger_left<Side>(){ return false; }


  //==================================================================


  template<>
  void interp_bounds_helper<XDIR>( const bool bnds[], const bool loBnds[], const bool hiBnds[],
				   const OpInfo& opInfo,
				   double & val,
				   int & ifac )
  {
    if( !bnds[0]   ) val = 0.0;
    if( !loBnds[0] ) ifac= 1;
    if( !hiBnds[0] ) ifac=-1;
  }
  template<>
  void interp_bounds_helper<YDIR>( const bool bnds[], const bool loBnds[], const bool hiBnds[],
				   const OpInfo& opInfo,
				   double & val,
				   int & ifac )
  {
    if( !bnds[1]   ) val = 0.0;
    if( !loBnds[1] ) ifac= 1;
    if( !hiBnds[1] ) ifac=-1;
    ifac *= opInfo.get_nxs();
  }
  template<>
  void interp_bounds_helper<ZDIR>( const bool bnds[], const bool loBnds[], const bool hiBnds[],
				   const OpInfo& opInfo,
				   double & val,
				   int & ifac )
  {
    if( !bnds[2]   ) val = 0.0;
    if( !loBnds[2] ) ifac= 1;
    if( !hiBnds[2] ) ifac=-1;
    ifac *= opInfo.get_nxs()*opInfo.get_nys();
  }


  //==================================================================


  template<>
  void grad_bounds_helper<XDIR>( const bool bnds[3], const bool loBnds[3], const bool hiBnds[3],
				 const std::vector<double> & spacing,
				 const OpInfo& opInfo,
				 double & val, int & ifac, int & shift )
  {
    val = bnds[0] ? 1.0/spacing[0] : 0.0;
    if( !loBnds[0] ) ifac= 1;
    if( !hiBnds[0] ) ifac=-1;
    shift = ifac;
  }
  template<>
  void grad_bounds_helper<YDIR>( const bool bnds[3], const bool loBnds[3], const bool hiBnds[3],
				 const std::vector<double> & spacing,
				 const OpInfo& opInfo,
				 double & val, int & ifac, int & shift )
  {
    val = bnds[1] ? 1.0/spacing[1] : 0.0;
    if( !loBnds[1] ) ifac= 1;
    if( !hiBnds[1] ) ifac=-1;
    shift = ifac*opInfo.get_nxs();
  }
  template<>
  void grad_bounds_helper<ZDIR>( const bool bnds[3], const bool loBnds[3], const bool hiBnds[3],
				 const std::vector<double> & spacing,
				 const OpInfo& opInfo,
				 double & val, int & ifac, int & shift )
  {
    val = bnds[2] ? 1.0/spacing[2] : 0.0;
    if( !loBnds[2] ) ifac= 1;
    if( !hiBnds[2] ) ifac=-1;
    shift = ifac * opInfo.get_nxs()*opInfo.get_nys();
  }


  //==================================================================


  template<>
  void div_bounds_helper<XDIR>( const bool bnds[3], const bool loBnds[3], const bool hiBnds[3],
				const std::vector<double> & faceArea,
				const double cellVol,
				const OpInfo& opInfo,
				double & val, int & ifac, int& shift )
  {
    val = bnds[0] ? faceArea[0]/cellVol : 0.0;
    if( !loBnds[0] ) ifac= 1;
    if( !hiBnds[0] ) ifac=-1;
    shift = ifac;
  }
  template<>
  void div_bounds_helper<YDIR>( const bool bnds[3], const bool loBnds[3], const bool hiBnds[3],
				const std::vector<double> & faceArea,
				const double cellVol,
				const OpInfo& opInfo,
				double & val, int & ifac, int& shift )
  {
    val = bnds[1] ? faceArea[1]/cellVol : 0.0;
    if( !loBnds[1] ) ifac= 1;
    if( !hiBnds[1] ) ifac=-1;
    shift = ifac * opInfo.get_nxs();
  }
  template<>
  void div_bounds_helper<ZDIR>( const bool bnds[3], const bool loBnds[3], const bool hiBnds[3],
				const std::vector<double> & faceArea,
				const double cellVol,
				const OpInfo& opInfo,
				double & val, int & ifac, int& shift )
  {
    val = bnds[2] ? faceArea[2]/cellVol : 0.0;
    if( !loBnds[2] ) ifac= 1;
    if( !hiBnds[2] ) ifac=-1;
    shift = ifac * opInfo.get_nxs()*opInfo.get_nys();
  }


  //==================================================================

  template<>
  void scratch_bounds_helper<XDIR>( const bool loBnds[3],
				    const OpInfo& opInfo,
				    int & skip,
				    int & shift )
  {
    skip = 1;
    if( !loBnds[0] ) ++shift;
  }
  template<>
  void scratch_bounds_helper<YDIR>( const bool loBnds[3],
				    const OpInfo& opInfo,
				    int & skip,
				    int & shift )
  {
    skip = opInfo.get_nxs();
    if( !loBnds[1] ) ++shift;
  }
  template<>
  void scratch_bounds_helper<ZDIR>( const bool loBnds[3],
				    const OpInfo& opInfo,
				    int & skip,
				    int & shift )
  {
    skip = opInfo.get_nxs()*opInfo.get_nys();
    if( !loBnds[2] ) ++shift;
  }

  //==================================================================


  template<>
  template<>
  int DefaultSideGhosting<XDIR>::get<XDIR,SidePlus>(){return 2;}

  template<>
  template<>
  int DefaultSideGhosting<YDIR>::get<YDIR,SidePlus>(){return 2;}

  template<>
  template<>
  int DefaultSideGhosting<ZDIR>::get<ZDIR,SidePlus>(){return 2;}


  //==================================================================


} // namespace FVStaggeredUniform
} // namespace SpatialOps
