#ifndef SO_BoundaryCondition_h
#define SO_BoundaryCondition_h

#include <boost/static_assert.hpp>

#include <SpatialField.h>
#include <SpatialOperator.h>


namespace SpatialOps{
namespace FVStaggeredUniform{

  //==================================================================

  /**
   *  Given a flat index on some field with ghost cells, this function
   *  computes the (i,j,k) indices for the interior of the domain.
   *
   *  If the flat index resides in a ghost cell, an assertion will
   *  result.
   *
   *  @param ix   The flat index on the given field with the GhostTraits
   *              specified by the template parameter.
   *  @param nx   The number of interior points in the x-direction
   *  @param ny   The number of interior points in the y-direction
   *
   *  @param i    The x-direction index for the interior of the domain.
   *  @param j    The y-direction index for the interior of the domain.
   *  @param k    The k-direction index for the interior of the domain.
   *
   *  @par Template Parameters
   *
   * \li \b GhostTraits specifies the ghost strategy for the field
   * whose indices we are manipulating.  If the field has no ghost
   * cells, then this reduces to a simple map between 3D and 1D
   * indexing.
   */
  template< typename GhostTraits >
  inline void flatix_to_ijk_interior( const int ix,
				      const std::vector<int>& nxyz,
				      int& i, int& j, int& k )
  {
    static const int ngxm = GhostTraits::template get<XDIR,SideMinus>();
    static const int ngym = GhostTraits::template get<YDIR,SideMinus>();
    static const int ngzm = GhostTraits::template get<ZDIR,SideMinus>();
    static const int ngx  = ngxm + GhostTraits::template get<XDIR,SidePlus >();
    static const int ngy  = ngym + GhostTraits::template get<YDIR,SidePlus >();

    const int nxt = nxyz[0]+ngx;
    const int nyt = nxyz[1]>1 ? nxyz[1]+ngy : 1;

    i =             ix % nxt       - ngxm;      assert( i>=0 );  assert( i<nxyz[0] );
    j = nxyz[1]>1 ? ix / nxt % nyt - ngym : 0;  assert( j>=0 );  assert( j<nxyz[1] );
    k = nxyz[2]>1 ? ix /(nxt * nyt)- ngzm : 0;  assert( k>=0 );  assert( k<nxyz[2] );
  }

  //==================================================================

  /**
   *  Calculates a flat index in a field given the interior ijk indices.
   *
   *  @param i    The x-direction index for the interior of the domain.
   *  @param j    The y-direction index for the interior of the domain.
   *  @param k    The k-direction index for the interior of the domain.
   *  @param nx   The number of interior points in the x-direction
   *  @param ny   The number of interior points in the y-direction
   *
   *  @param ix   The flat index on the given field with the GhostTraits
   *              specified by the template parameter.
   *
   *  @par Template Parameters
   *
   * \li \b GhostTraits specifies the ghost strategy for the field
   * whose indices we are manipulating.  If the field has no ghost
   * cells, then this reduces to a simple map between 3D and 1D
   * indexing.
   */
  template< typename GhostTraits >
  inline void ijk_interior_to_flatix( const int i,  const int j, const int k,
				      const std::vector<int>& nxyz,
				      int& ix )
  {
    static const int ngxm = GhostTraits::template get<XDIR,SideMinus>();
    static const int ngym = GhostTraits::template get<YDIR,SideMinus>();
    static const int ngzm = GhostTraits::template get<ZDIR,SideMinus>();
    static const int ngx  = ngxm + GhostTraits::template get<XDIR,SidePlus >();
    static const int ngy  = ngym + GhostTraits::template get<YDIR,SidePlus >();

    // add in x-dir with ghost offests
    ix = i + ngxm;

    // add in y-dir, with ghosts offsets
    const int nxt = nxyz[0]+ngx;
    if( nxyz[1]>1 ) ix += (j+ngym)*nxt;

    const int nyt = nxyz[1]+ngy;
    // add in z-dir (if ny>1)  with ghost offsets
    if( nxyz[2]>1 ) ix += (k+ngzm)*(nxt*nyt);
  }

  //==================================================================

  namespace BCToolsLocal{

    /**
     *  @struct Shifter
     *  @author James C. Sutherland
     *  @date   June, 2007
     *
     *  @brief Helper class for shifting field indices by one toward
     *  the side indicated and in the direction indicated.
     */
    template<typename GhostTraits,typename Dir,typename Side>
    struct Shifter
    {
      /** This method is only implemented in specialized versions of Shifter */
      static void shift_index( const std::vector<int>& nxyz, int& ix );
    };
    //----------------------------------------------------------------
    /**
     *  @brief Shift the index to the (-) x.
     */
    template<typename GhostTraits>
    struct Shifter<GhostTraits,XDIR,SideMinus>
    {
      static void shift_index( const std::vector<int>& nxyz, int& ix )
      {
	--ix;
      }
    };
    //----------------------------------------------------------------
    /**
     *  @brief Shift the index to the (+) x.
     */
    template<typename GhostTraits>
    struct Shifter<GhostTraits,XDIR,SidePlus>
    {
      static void shift_index( const std::vector<int>& nxyz, int& ix )
      {
	++ix;
      }
    };
    //----------------------------------------------------------------
    /**
     *  @brief Shift the index to the (-) y.
     */
    template<typename GhostTraits>
    struct Shifter<GhostTraits,YDIR,SideMinus>
    {
      static void shift_index( const std::vector<int>& nxyz, int& ix )
      {
	ix -= nxyz[0] + GhostTraits::template get<XDIR,SideMinus>() + GhostTraits::template get<XDIR,SidePlus>();
      }
    };
    //----------------------------------------------------------------
    /**
     *  @brief Shift the index to the (+) y.
     */
    template<typename GhostTraits>
    struct Shifter<GhostTraits,YDIR,SidePlus>
    {
      static void shift_index( const std::vector<int>& nxyz, int& ix )
      {
	ix += nxyz[0] + GhostTraits::template get<XDIR,SideMinus>() + GhostTraits::template get<XDIR,SidePlus>();
      }
    };
    //----------------------------------------------------------------
    /**
     *  @brief Shift the index to the (-) z.
     */
    template<typename GhostTraits>
    struct Shifter<GhostTraits,ZDIR,SideMinus>
    {
      static void shift_index(const std::vector<int>& nxyz, int& ix)
      {
	const int nxt = nxyz[0] + GhostTraits::template get<XDIR,SideMinus>() + GhostTraits::template get<XDIR,SidePlus>();
	const int nyt = nxyz[1] + GhostTraits::template get<YDIR,SideMinus>() + GhostTraits::template get<YDIR,SidePlus>();
	ix -= nxt*nyt;
      }
    };
    //---------------------------------------------------------------- 
    /**
     *  @brief Shift the index to the (+) z.
     */
   template<typename GhostTraits>
    struct Shifter<GhostTraits,ZDIR,SidePlus>
    {
      static void shift_index(const std::vector<int>& nxyz, int& ix)
      {
	const int nxt = nxyz[0] + GhostTraits::template get<XDIR,SideMinus>() + GhostTraits::template get<XDIR,SidePlus>();
	const int nyt = nxyz[1] + GhostTraits::template get<YDIR,SideMinus>() + GhostTraits::template get<YDIR,SidePlus>();
	ix += nxt*nyt;
      }
    };
    //----------------------------------------------------------------

    //================================================================

    /**
     *  @struct RowIndexShifter
     *  @author James C. Sutherland
     *  @date   June, 2007
     *
     *  Given an interior index, on the left side of the domain, that
     *  is the correct index for the row of the operator.  On the
     *  right side of the domain, it is not correct, and must be
     *  shifted.  RowIndexShifter has one static method, which is only
     *  implemented in specialized versions of the struct.
     */
    template<typename GhostTraits, typename Dir, typename Side>
    struct RowIndexShifter{
      static void apply( const std::vector<int>& nxyz, int& ix );
    };
    //----------------------------------------------------------------
    template<typename GhostTraits>
    struct RowIndexShifter<GhostTraits,XDIR,SideMinus>{
      static void apply( const std::vector<int>& nxyz, int& ix ){}
    };
    //----------------------------------------------------------------
    template<typename GhostTraits>
    struct RowIndexShifter<GhostTraits,XDIR,SidePlus>{
      static void apply( const std::vector<int>& nxyz, int& ix )
      {
	++ix;
      }
    };
    //----------------------------------------------------------------
    template<typename GhostTraits>
    struct RowIndexShifter<GhostTraits,YDIR,SideMinus>{
      static void apply( const std::vector<int>& nxyz, int& ix ){}
    };
    //----------------------------------------------------------------
    template<typename GhostTraits>
    struct RowIndexShifter<GhostTraits,YDIR,SidePlus>{
      static void apply( const std::vector<int>& nxyz, int& ix )
      {
	static const int ngx = GhostTraits::template get<XDIR,SideMinus>() + GhostTraits::template get<XDIR,SidePlus>();
	ix += nxyz[0] + ngx;
      }
    };
    //----------------------------------------------------------------
    template<typename GhostTraits>
    struct RowIndexShifter<GhostTraits,ZDIR,SideMinus>{
      static void apply( const std::vector<int>& nxyz, int& ix ){}
    };
    //----------------------------------------------------------------
    template<typename GhostTraits>
    struct RowIndexShifter<GhostTraits,ZDIR,SidePlus>{
      static void apply( const std::vector<int>& nxyz, int& ix )
      {
	static const int ngx = GhostTraits::template get<XDIR,SideMinus>() + GhostTraits::template get<XDIR,SidePlus>();
	static const int ngy = GhostTraits::template get<YDIR,SideMinus>() + GhostTraits::template get<YDIR,SidePlus>();
	ix += (nxyz[0] + ngx) * (nxyz[1] + ngy);
      }
    };
    //----------------------------------------------------------------

  } // namespace BCToolsLocal

  //==================================================================

  /**
   *  @brief Apply either a Dirichlet or Neumann BC (depending on the
   *  operator supplied) to the given field.
   *
   *  @param op The operator to use in applying the BC.  Supplying a
   *            gradient operator results in a Neumann BC; supplying
   *            an interpolant operator results in a Dirichlet BC.
   *
   *  @param i  The x-direction index for the cell we want to apply the BC to.
   *  @param j  The y-direction index for the cell we want to apply the BC to.
   *  @param k  The z-direction index for the cell we want to apply the BC to.
   *
   *  @param bcVal The value for the boundary condition to set.
   *
   *  @param f  The field we want to set the BC on.
   *
   *  @par Template Parameters
   *   <ul>
   *
   *   <li> \b OpT Specifies the type of SpatialOperator that will be
   *   used to apply this BC.  This must define a few things:
   *     <ul>
   *     <li> \b DestGhost The ghost traits for the destination field.
   *     <li> \b SrcGhost  The ghost traits for the source field.
   *     <li> \b DirType   The direction that this operator applies in. See XDIR, YDIR, ZDIR.
   *     </ul>
   *
   *   <li> \b FieldT specifies the type of SpatialField that we are
   *        setting the BC on.  Note that FieldT must be compatible with
   *        the source field for the operator.  This must define a few types:
   *     <ul>
   *     <li> \b Ghost the ghosting strategy for this field.
   *     <li> \b Location the location for this field.
   *     </ul>
   *
   *   <li> \b %Side specifies the side (+/-) to specify the BC on.  See SideMinus SidePlus
   *
   *   </ul>
   */
  template< typename OpT, typename FieldT, typename Side >
  void assign_bc_point( const OpT& op, const int i, const int j, const int k, const double bcVal, FieldT& f )
  {
    // ensure that at a minimum the source field for the operator is
    // ghosted the same as the field we are applying the BC on.
    BOOST_STATIC_ASSERT( bool(IsSameType< typename OpT::SrcGhost,    typename FieldT::Ghost    >::result) );
    BOOST_STATIC_ASSERT( bool(IsSameType< typename OpT::SrcLocation, typename FieldT::Location >::result) );


    const std::vector<int>& extent = f.get_extent();

    int irow=-1;     int ixf=-1;
    ijk_interior_to_flatix<typename OpT::SrcGhost >( i,j,k, extent, ixf  );
    ijk_interior_to_flatix<typename OpT::DestGhost>( i,j,k, extent, irow );

    BCToolsLocal::RowIndexShifter<typename OpT::DestGhost,typename OpT::DirType,Side>::apply( extent, irow );
    BCToolsLocal::Shifter<typename OpT::SrcGhost,typename OpT::DirType,Side>::shift_index( extent, ixf );

    int ncol; double*vals=0; int*ixs=0;
    op.get_linalg_mat().ExtractMyRowView( irow, ncol, vals, ixs );

    double prodsum=0.0; double ghostcoeff=0.0;
    for (int icol = 0; icol < ncol; icol++){
      if (*ixs == ixf)	ghostcoeff = *vals;
      else              prodsum += (*vals) *f[*ixs];
      ++ixs; ++vals;
    }

    assert( ghostcoeff != 0.0 );
    f[ixf] = (bcVal - prodsum) / ghostcoeff;
  }

  //====================================================================


} // namespace FVStaggeredUniform
} // namespace SpatialOps


#endif
