#ifndef FVStaggeredTools_h
#define FVStaggeredTools_h

#include <boost/static_assert.hpp>

#include <FVStaggeredTypes.h>
#include <SpatialOpsTools.h> // for IsSameType

namespace SpatialOps{
namespace FVStaggered{


  //==================================================================

  /**
   * @brief get the total number of points (including ghost cells) for
   * a field in the x-direction.
   */
  template<typename FieldT> inline int get_nx( const int nxInterior )
  {
    typedef typename FieldT::Ghost G;
    int npts=1;
    if( nxInterior>1 ){
      npts = nxInterior + G::NM + G::NP;
      if( IsSameType<typename FieldT::Location, SSurfX>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, XVol  >::result ) ++npts;
      if( IsSameType<typename FieldT::Location, XSurfX>::result ) npts+=2;
      if( IsSameType<typename FieldT::Location, XSurfY>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, XSurfZ>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, YSurfX>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, ZSurfX>::result ) ++npts;
   }
    return npts;
  }

  //==================================================================
 
  /**
   * @brief get the total number of points (including ghost cells) for
   * a field in the y-direction.
   */
  template<typename FieldT> inline int get_ny( const int nyInterior )
  {
    typedef typename FieldT::Ghost G;
    int npts=1;
    if( nyInterior>1 ){
      npts = nyInterior + G::NM + G::NP;
      if( IsSameType<typename FieldT::Location, SSurfY>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, XSurfY>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, YVol  >::result ) ++npts;
      if( IsSameType<typename FieldT::Location, YSurfX>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, YSurfY>::result ) npts+=2;
      if( IsSameType<typename FieldT::Location, YSurfZ>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, ZSurfY>::result ) ++npts;
   }
    return npts;
  }

  //==================================================================

  /**
   * @brief get the total number of points (including ghost cells) for
   * a field in the z-direction.
   */
  template<typename FieldT> inline int get_nz( const int nzInterior )
  {
    typedef typename FieldT::Ghost G;
    int npts=1;
    if( nzInterior>1 ){
      npts = nzInterior + G::NM + G::NP;
      if( IsSameType<typename FieldT::Location, SSurfZ>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, XSurfZ>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, YSurfZ>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, ZVol  >::result ) ++npts;
      if( IsSameType<typename FieldT::Location, ZSurfX>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, ZSurfY>::result ) ++npts;
      if( IsSameType<typename FieldT::Location, ZSurfZ>::result ) npts+=2;
    }
    return npts;
  }

  //==================================================================

  /**
   * @brief get the total number of points in the x-direction
   * (including ghost cells) for a x-surface field.
   */
  template<typename FieldT> inline int get_nx_x( const std::vector<int>& dim )
  {
    BOOST_STATIC_ASSERT( FieldT::Location::IsSurface );
    switch( FieldT::Location::StagDir::value ){
    case YDIR::value: if( dim[1]==1 ) return 1; break;
    case ZDIR::value: if( dim[2]==1 ) return 1;
    }
    return get_nx<FieldT>(dim[0]);
  }

  // explicitly declare functions without implementations for field
  // types that this doesn't make sense for.
  template<> int get_nx_x<SSurfYField>( const std::vector<int>& dim );
  template<> int get_nx_x<SSurfZField>( const std::vector<int>& dim );
  template<> int get_nx_x<XSurfYField>( const std::vector<int>& dim );
  template<> int get_nx_x<XSurfZField>( const std::vector<int>& dim );
  template<> int get_nx_x<YSurfYField>( const std::vector<int>& dim );
  template<> int get_nx_x<YSurfZField>( const std::vector<int>& dim );
  template<> int get_nx_x<ZSurfYField>( const std::vector<int>& dim );
  template<> int get_nx_x<ZSurfZField>( const std::vector<int>& dim );


  //==================================================================

  /**
   * @brief get the total number of points in the x-direction
   * (including ghost cells) for a y-surface field.
   */
  template<typename FieldT> inline int get_nx_y( const std::vector<int>& dim )
  {
    BOOST_STATIC_ASSERT( FieldT::Location::IsSurface );
    switch( FieldT::Location::StagDir::value ){
    case XDIR::value: if( dim[0]==1 ) return 1; break;
    case ZDIR::value: if( dim[2]==1 ) return 1;
    }
    int npts = 1;
    if( dim[0]>1 ) npts = get_nx<FieldT>(dim[0]);
    return npts;
  }

  // explicitly declare functions without implementations for field
  // types that this doesn't make sense for.
  template<> int get_nx_y<SSurfXField>( const std::vector<int>& dim );
  template<> int get_nx_y<SSurfZField>( const std::vector<int>& dim );
  template<> int get_nx_y<XSurfXField>( const std::vector<int>& dim );
  template<> int get_nx_y<XSurfZField>( const std::vector<int>& dim );
  template<> int get_nx_y<YSurfXField>( const std::vector<int>& dim );
  template<> int get_nx_y<YSurfZField>( const std::vector<int>& dim );
  template<> int get_nx_y<ZSurfXField>( const std::vector<int>& dim );
  template<> int get_nx_y<ZSurfZField>( const std::vector<int>& dim );

  //==================================================================

  /**
   * @brief get the total number of points in the x-direction
   * (including ghost cells) for a z-surface field.
   */
  template<typename FieldT> inline int get_nx_z( const std::vector<int>& dim )
  {
    BOOST_STATIC_ASSERT( FieldT::Location::IsSurface );
    switch( FieldT::Location::StagDir::value ){
    case XDIR::value: if( dim[0]==1 ) return 1; break;
    case YDIR::value: if( dim[1]==1 ) return 1;
    }
    int npts = 1;
    if( dim[0]>1 ) npts = get_nx<FieldT>(dim[0]);
    return npts;
  }

  // explicitly declare functions without implementations for field
  // types that this doesn't make sense for.
  template<> int get_nx_z<SSurfXField>( const std::vector<int>& dim );
  template<> int get_nx_z<SSurfYField>( const std::vector<int>& dim );
  template<> int get_nx_z<XSurfXField>( const std::vector<int>& dim );
  template<> int get_nx_z<XSurfYField>( const std::vector<int>& dim );
  template<> int get_nx_z<YSurfXField>( const std::vector<int>& dim );
  template<> int get_nx_z<YSurfYField>( const std::vector<int>& dim );
  template<> int get_nx_z<ZSurfXField>( const std::vector<int>& dim );
  template<> int get_nx_z<ZSurfYField>( const std::vector<int>& dim );

  //==================================================================

  /**
   * @brief get the total number of points in the y-direction
   * (including ghost cells) for a x-surface field.
   */
  template<typename FieldT> inline int get_ny_x( const std::vector<int>& dim )
  {
    BOOST_STATIC_ASSERT( FieldT::Location::IsSurface );
    switch( FieldT::Location::StagDir::value ){
    case YDIR::value: if( dim[1]==1 ) return 1; break;
    case ZDIR::value: if( dim[2]==1 ) return 1;
    }
    int npts = 1;
    if( dim[1]>1 ) npts = get_ny<FieldT>(dim[1]);
    return npts;
  }

  // explicitly declare functions without implementations for field
  // types that this doesn't make sense for.
  template<> int get_ny_x<SSurfYField>( const std::vector<int>& dim );
  template<> int get_ny_x<SSurfZField>( const std::vector<int>& dim );
  template<> int get_ny_x<XSurfYField>( const std::vector<int>& dim );
  template<> int get_ny_x<XSurfZField>( const std::vector<int>& dim );
  template<> int get_ny_x<YSurfYField>( const std::vector<int>& dim );
  template<> int get_ny_x<YSurfZField>( const std::vector<int>& dim );
  template<> int get_ny_x<ZSurfYField>( const std::vector<int>& dim );
  template<> int get_ny_x<ZSurfZField>( const std::vector<int>& dim );

  //==================================================================

  /**
   * @brief get the total number of points in the y-direction
   * (including ghost cells) for a y-surface field.
   */
  template<typename FieldT> inline int get_ny_y( const std::vector<int>& dim )
  {
    BOOST_STATIC_ASSERT( FieldT::Location::IsSurface );
    switch( FieldT::Location::StagDir::value ){
    case XDIR::value: if( dim[0]==1 ) return 1; break;
    case ZDIR::value: if( dim[2]==1 ) return 1;
    }
    return dim[1]>1 ? get_ny<FieldT>(dim[1]) : 1;
  }

  // explicitly declare functions without implementations for field
  // types that this doesn't make sense for.
  template<> int get_ny_y<SSurfXField>( const std::vector<int>& dim );
  template<> int get_ny_y<SSurfZField>( const std::vector<int>& dim );
  template<> int get_ny_y<XSurfXField>( const std::vector<int>& dim );
  template<> int get_ny_y<XSurfZField>( const std::vector<int>& dim );
  template<> int get_ny_y<YSurfXField>( const std::vector<int>& dim );
  template<> int get_ny_y<YSurfZField>( const std::vector<int>& dim );
  template<> int get_ny_y<ZSurfXField>( const std::vector<int>& dim );
  template<> int get_ny_y<ZSurfZField>( const std::vector<int>& dim );

  //==================================================================

  /**
   * @brief get the total number of points in the y-direction
   * (including ghost cells) for a z-surface field.
   */
  template<typename FieldT> inline int get_ny_z( const std::vector<int>& dim )
  {
    BOOST_STATIC_ASSERT( FieldT::Location::IsSurface );
    switch( FieldT::Location::StagDir::value ){
    case XDIR::value: if( dim[0]==1 ) return 1; break;
    case YDIR::value: if( dim[1]==1 ) return 1;
    }
    int npts = 1;
    if( dim[1]>1 ) npts = get_ny<FieldT>(dim[1]);
    return npts;
  }

  // explicitly declare functions without implementations for field
  // types that this doesn't make sense for.
  template<> int get_ny_z<SSurfXField>( const std::vector<int>& dim );
  template<> int get_ny_z<SSurfYField>( const std::vector<int>& dim );
  template<> int get_ny_z<XSurfXField>( const std::vector<int>& dim );
  template<> int get_ny_z<XSurfYField>( const std::vector<int>& dim );
  template<> int get_ny_z<YSurfXField>( const std::vector<int>& dim );
  template<> int get_ny_z<YSurfYField>( const std::vector<int>& dim );
  template<> int get_ny_z<ZSurfXField>( const std::vector<int>& dim );
  template<> int get_ny_z<ZSurfYField>( const std::vector<int>& dim );

  //==================================================================

  /**
   * @brief get the total number of points in the z-direction
   * (including ghost cells) for a x-surface field.
   */
  template<typename FieldT> inline int get_nz_x( const std::vector<int>& dim )
  {
    BOOST_STATIC_ASSERT( FieldT::Location::IsSurface );
    switch( FieldT::Location::StagDir::value ){
    case YDIR::value: if( dim[1]==1 ) return 1; break;
    case ZDIR::value: if( dim[2]==1 ) return 1;
    }
    int npts = 1;
    if( dim[2]>1 ) npts = get_nz<FieldT>(dim[2]);
    return npts;
  }

  // explicitly declare functions without implementations for field
  // types that this doesn't make sense for.
  template<> int get_nz_x<SSurfYField>( const std::vector<int>& dim );
  template<> int get_nz_x<SSurfZField>( const std::vector<int>& dim );
  template<> int get_nz_x<XSurfYField>( const std::vector<int>& dim );
  template<> int get_nz_x<XSurfZField>( const std::vector<int>& dim );
  template<> int get_nz_x<YSurfYField>( const std::vector<int>& dim );
  template<> int get_nz_x<YSurfZField>( const std::vector<int>& dim );
  template<> int get_nz_x<ZSurfYField>( const std::vector<int>& dim );
  template<> int get_nz_x<ZSurfZField>( const std::vector<int>& dim );

  //==================================================================

  /**
   * @brief get the total number of points in the z-direction
   * (including ghost cells) for a y-surface field.
   */
  template<typename FieldT> inline int get_nz_y( const std::vector<int>& dim )
  {
    BOOST_STATIC_ASSERT( FieldT::Location::IsSurface );
    switch( FieldT::Location::StagDir::value ){
    case XDIR::value: if( dim[0]==1 ) return 1;  break;
    case ZDIR::value: if( dim[2]==1 ) return 1;
    }
    int npts = 1;
    if( dim[2]>1 ) npts = get_nz<FieldT>(dim[2]);
    return npts;
  }

  // explicitly declare functions without implementations for field
  // types that this doesn't make sense for.
  template<> int get_nz_y<SSurfXField>( const std::vector<int>& dim );
  template<> int get_nz_y<SSurfZField>( const std::vector<int>& dim );
  template<> int get_nz_y<XSurfXField>( const std::vector<int>& dim );
  template<> int get_nz_y<XSurfZField>( const std::vector<int>& dim );
  template<> int get_nz_y<YSurfXField>( const std::vector<int>& dim );
  template<> int get_nz_y<YSurfZField>( const std::vector<int>& dim );
  template<> int get_nz_y<ZSurfXField>( const std::vector<int>& dim );
  template<> int get_nz_y<ZSurfZField>( const std::vector<int>& dim );

  //==================================================================

  /**
   * @brief get the total number of points in the z-direction
   * (including ghost cells) for a z-surface field.
   */
  template<typename FieldT> inline int get_nz_z( const std::vector<int>& dim )
  {
    BOOST_STATIC_ASSERT( FieldT::Location::IsSurface );
    switch( FieldT::Location::StagDir::value ){
    case XDIR::value: if( dim[0]==1 ) return 1;  break;
    case YDIR::value: if( dim[1]==1 ) return 1;
    }
    return dim[2]>1 ? get_nz<FieldT>(dim[2]) : 1;
  }

  // explicitly declare functions without implementations for field
  // types that this doesn't make sense for.
  template<> int get_nz_z<SSurfXField>( const std::vector<int>& dim );
  template<> int get_nz_z<SSurfYField>( const std::vector<int>& dim );
  template<> int get_nz_z<XSurfXField>( const std::vector<int>& dim );
  template<> int get_nz_z<XSurfYField>( const std::vector<int>& dim );
  template<> int get_nz_z<YSurfXField>( const std::vector<int>& dim );
  template<> int get_nz_z<YSurfYField>( const std::vector<int>& dim );
  template<> int get_nz_z<ZSurfXField>( const std::vector<int>& dim );
  template<> int get_nz_z<ZSurfYField>( const std::vector<int>& dim );

  //==================================================================

  /**
   * @brief get the total number of points in a field, including ghost
   * cells.
   */
  template<typename FieldT> int get_n_tot( const std::vector<int>& dim )
  {
    BOOST_STATIC_ASSERT( FieldT::Location::IsSurface );

    switch( FieldT::Location::Dir::value ){
    case NODIR::value:  // surface fields (scalars located at all faces)
      {
	int n=0;
	if( dim[0]>1 )  // number of points for "x-face" field
	  n += (            get_nx_x<FieldT>(dim)     )
	    *  ( dim[1]>1 ? get_ny_x<FieldT>(dim) : 1 )
	    *  ( dim[2]>1 ? get_nz_x<FieldT>(dim) : 1 );
	if( dim[1]>1 )  // number of points for 'y-face" field
	  n += ( dim[0]>1 ? get_nx_y<FieldT>(dim) : 1 )
	    *  (            get_ny_y<FieldT>(dim)     )
	    *  ( dim[2]>1 ? get_nz_y<FieldT>(dim) : 1 );
	if( dim[2]>1 )  // number of points for "z-face" field
	  n += ( dim[0]>1 ? get_nx_z<FieldT>(dim) : 1 )
	    *  ( dim[1]>1 ? get_ny_z<FieldT>(dim) : 1 )
	    *  (            get_nz_z<FieldT>(dim)     );
	return n;
      }
    case XDIR::value:  // x-face vector component
      {
	int n=1;
	if( dim[0]>1 )
	  n = get_nx_x<FieldT>(dim) * get_ny_x<FieldT>(dim) * get_nz_x<FieldT>(dim);
	return n;
      }
    case YDIR::value:  // y-face vector component
      {
	int n=1;
	if( dim[1]>1 )
	  n = get_nx_y<FieldT>(dim) * get_ny_y<FieldT>(dim) * get_nz_y<FieldT>(dim);
	return n;
      }
    case ZDIR::value:  // z-face vector component
      {
	int n=1;
	if( dim[2]>1 )
	  n = get_nx_z<FieldT>(dim) * get_ny_z<FieldT>(dim) * get_nz_z<FieldT>(dim);
	return n;
      }
    }
    // should never get here.
    return -1;
  }

  template<> inline int get_n_tot<SVolField>( const std::vector<int>& dim )
  {
    return get_nx<SVolField>(dim[0]) * get_ny<SVolField>(dim[1]) * get_nz<SVolField>(dim[2]);
  }

  template<> inline int get_n_tot<SVolRHS>( const std::vector<int>& dim )
  {
    return get_nx<SVolRHS>(dim[0]) * get_ny<SVolRHS>(dim[1]) * get_nz<SVolRHS>(dim[2]);
  }

  template<> inline int get_n_tot<XVolField>( const std::vector<int>& dim )
  {
    int n=1;
    if( dim[0]>1 )
      n = get_nx<XVolField>(dim[0]) * get_ny<XVolField>(dim[1]) * get_nz<XVolField>(dim[2]);
    return n;
  }

  template<> inline int get_n_tot<YVolField>( const std::vector<int>& dim )
  {
    int n=1;
    if( dim[1]>1 )
      n = get_nx<YVolField>(dim[0]) * get_ny<YVolField>(dim[1]) * get_nz<YVolField>(dim[2]);
    return n;
  }

  template<> inline int get_n_tot<ZVolField>( const std::vector<int>& dim )
  {
    int n=1;
    if( dim[2]>1 )
      n = get_nx<ZVolField>(dim[0]) * get_ny<ZVolField>(dim[1]) * get_nz<ZVolField>(dim[2]);
    return n;
  }

  //==================================================================

  inline void _ghost_set_( const int ngm, const int ngp,
			   const int nxt, const int nyt, const int nzt,
			   const std::vector<int>& dim,
			   int& ix,
			   std::set<int>& ghostSet )
  {
    const int ngxm = dim[0]>1 ? ngm : 0;
    const int ngxp = dim[0]>1 ? ngp : 0;
    const int ngym = dim[1]>1 ? ngm : 0;
    const int ngyp = dim[1]>1 ? ngp : 0;
    const int ngzm = dim[2]>1 ? ngm : 0;
    const int ngzp = dim[2]>1 ? ngp : 0;

    // -z side ghost layer
    if( dim[2]>1 ){
      for( int kg=0; kg<ngzm; ++kg )
	for( int j=0; j<nyt; ++j )
	  for( int i=0; i<nxt; ++i )
	    ghostSet.insert(ix++);
    }

    // z interior
    for( int k=ngzm; k<nzt-ngzp; ++k ){

      // -y side ghost layer
      if( dim[1]>1 ){
	for( int i=0; i<nxt; ++i )
	  for( int jg=0; jg<ngym; ++jg )
	    ghostSet.insert(ix++);
      }

      // y interior
      for( int j=ngym; j<nyt-ngyp; ++j ){
	// -x side ghost layer
	if( dim[0]>1 ) for( int ig=0; ig<ngxm; ++ig ) ghostSet.insert(ix++);
	// x interior
	ix+=nxt-ngxm-ngxp;
	// +x side ghost layer
	if( dim[0]>1 ) for( int ig=0; ig<ngxp; ++ig ) ghostSet.insert(ix++);
      }

      // +y side ghost layer
      if( dim[1]>1 ){
	for( int i=0; i<nxt; ++i )
	  for( int jg=0; jg<ngyp; ++jg )
	    ghostSet.insert(ix++);
      }
    }

    // +z side ghost layer
    if( dim[2]>1 ){
      for( int kg=0; kg<ngzp; ++kg )
	for( int i=0; i<nxt; ++i )
	  for( int j=0; j<nyt; ++j )
	    ghostSet.insert(ix++);
    }
  }

  //==================================================================

  /**
   *  @brief Obtain the set of indices corresponding to ghost cells
   *  for this field.
   */
  template<typename FieldT> 
  const std::set<int>& get_ghost_set( const std::vector<int>& dim )
  {
    typedef typename FieldT::Ghost G;
    static std::set<int> ghostSet;
    ghostSet.clear();
    int ix=0;
    _ghost_set_( G::NM, G::NP,
		 get_nx<FieldT>(dim[0]),
		 get_ny<FieldT>(dim[1]),
		 get_nz<FieldT>(dim[2]),
		 dim,
		 ix,
		 ghostSet );
    return ghostSet;
  }

  //==================================================================

  struct IndexTriplet
  {
    IndexTriplet( const int ii, const int jj, const int kk ): i(ii), j(jj), k(kk){}
    IndexTriplet(){ i=j=k=-1; }
    IndexTriplet( const IndexTriplet& x ){ i=x.i; j=x.j; k=x.k; }
    IndexTriplet& operator=(const IndexTriplet& x){ i=x.i; j=x.j; k=x.k; return *this; }
    int i,j,k;
  };

  //==================================================================

  /**
   *  @brief Use this to transform a flat index to i,j,k indices.
   */
  template<typename FieldT, int IsSurfField, typename CompType=NODIR>
  struct flat2ijk
  {
    static IndexTriplet value( const std::vector<int>& dim, const int ix );
  };
  template<typename FieldT>
  struct flat2ijk<FieldT,0>
  {
    static IndexTriplet value( const std::vector<int>& dim, const int ix );
  };
  template<typename FieldT>
  struct flat2ijk<FieldT,1>
  {
    static IndexTriplet value( const std::vector<int>& dim, const int ix );
  };

  //==================================================================

  /**
   *  @brief Use this to transform i,j,k indices to a flat index.
   */
  template<typename FieldT, int IsSurfField, typename CompType=NODIR>
  struct ijk2flat
  {
    static int value( const std::vector<int>& dim, const IndexTriplet& ixt );
  };
  template<typename FieldT>
  struct ijk2flat<FieldT,0>
  {
    static int value( const std::vector<int>& dim, const IndexTriplet& ixt );
  };
  template<typename FieldT,typename CompType>
  struct ijk2flat<FieldT,1,CompType>
  {
    static int value( const std::vector<int>& dim, const IndexTriplet& ixt );
  };


  // implementation for all surface fields:
  template<typename FieldT>
  inline IndexTriplet
  flat2ijk<FieldT,1>::value( const std::vector<int>& dim, const int ix )
  {
    IndexTriplet triplet;

    int nxt=-1, nyt=-1;
    switch( FieldT::Location::Dir::value ){
    case XDIR::value:
      nxt = get_nx_x<FieldT>(dim);
      nyt = get_ny_x<FieldT>(dim);
      break;
    case YDIR::value:
      nxt = get_nx_y<FieldT>(dim);
      nyt = get_ny_y<FieldT>(dim);
      break;
    case ZDIR::value:
      nxt = get_nx_z<FieldT>(dim);
      nyt = get_ny_z<FieldT>(dim);
      break;
    case NODIR::value:
      {
	const int nx1 = get_nx_x<FieldT>(dim);
	const int ny1 = get_ny_x<FieldT>(dim);
	const int nz1 = get_nz_x<FieldT>(dim);
	const int nx2 = get_nx_y<FieldT>(dim);
	const int ny2 = get_ny_y<FieldT>(dim);
	const int nz2 = get_nz_y<FieldT>(dim);
	const int nx3 = get_nx_z<FieldT>(dim);
	const int ny3 = get_ny_z<FieldT>(dim);
	const int nz3 = get_nz_z<FieldT>(dim);
	const int n1 = dim[0]>1 ? nx1*ny1*nz1 : 0;
	const int n2 = dim[1]>1 ? n1 + nx2*ny2*nz2 : n1;
	const int n3 = dim[2]>1 ? n2 + nx3*ny3*nz3 : n2;

	if( ix<n1 ){
	  triplet.i = ix%nx1;
	  triplet.j = ix/nx1 % ny1;
	  triplet.k = ix/(nx1*ny1);
	}
	else if( ix<n2 ){
	  triplet.i = (nx2>1) ? (ix-n1)%nx2 : 0;
	  triplet.j = (ix-n1)/nx2 % ny2;
	  triplet.k = (ix-n1)/(nx2*ny2);
	}
	else if( ix<n3 ){
	  triplet.i = (nx3>1) ? (ix-n2)%nx3 : 0;
	  triplet.j = (ny3>1) ? (ix-n2)/nx3 % ny3 : 0;
	  triplet.k = (ix-n2)/(nx3*ny3);
	}
	else{
	  assert(0);  // should never get here.
	}
	return triplet;
      }
    default:
      assert(0);
    }
    triplet.i = ix%nxt;
    triplet.j = ix/nxt % nyt;
    triplet.k = ix/(nxt*nyt);

    return triplet;
  }

  // implementation for all volume fields
  template<typename FieldT>
  inline IndexTriplet
  flat2ijk<FieldT,0>::value( const std::vector<int>& dim, const int ix )
  {
    IndexTriplet triplet;

    const int nxt = get_nx<FieldT>(dim[0]);
    const int nyt = get_ny<FieldT>(dim[1]);

    triplet.i = ix%nxt;
    triplet.j = ix/nxt % nyt;
    triplet.k = ix/(nxt*nyt);

    return triplet;
  }

  //==================================================================

  // implementation for all surface fields
  template<typename FieldT,typename CompType>
  inline int
  ijk2flat<FieldT,1,CompType>::value( const std::vector<int>& dim, const IndexTriplet& triplet )
  {
    int nxt=-1, nyt=-1;
    switch( FieldT::Location::Dir::value ){
    case XDIR::value:
      nxt = get_nx_x<FieldT>(dim);
      nyt = get_ny_x<FieldT>(dim);
      break;
    case YDIR::value:
      nxt = get_nx_y<FieldT>(dim);
      nyt = get_ny_y<FieldT>(dim);
      break;
    case ZDIR::value:
      nxt = get_nx_z<FieldT>(dim);
      nyt = get_ny_z<FieldT>(dim);
      break;
    case NODIR::value:
      {
	switch( CompType::value ){
	case XDIR::value:
	  nxt = get_nx_x<FieldT>(dim);
	  nyt = get_ny_x<FieldT>(dim);
	  break;
	case YDIR::value:
	  nxt = get_nx_y<FieldT>(dim);
	  nyt = get_ny_y<FieldT>(dim);
	  break;
	case ZDIR::value:
	  nxt = get_nx_z<FieldT>(dim);
	  nyt = get_ny_z<FieldT>(dim);
	  break;
	}
	break;
      }
    }

    int flat = triplet.i + nxt*triplet.j + nxt*nyt*triplet.k;

    return flat;
  }

  // implementation for all volume fields
  template<typename FieldT>
  inline int
  ijk2flat<FieldT,0>::value( const std::vector<int>& dim, const IndexTriplet& triplet )
  {
    const int nxt = get_nx<FieldT>(dim[0]);
    const int nyt = get_ny<FieldT>(dim[1]);
      
    return
      triplet.i +
      triplet.j * nxt +
      triplet.k * nxt*nyt;
  }

  //==================================================================


}// namespace FVStaggered
}// namespace SpatialOps

#endif
