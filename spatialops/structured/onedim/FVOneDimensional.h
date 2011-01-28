#ifndef FVStaggeredOneDimStretch_h
#define FVStaggeredOneDimStretch_h

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOperator.h>

#include "LagrangePoly.h"
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>


#include <vector>
#include <set>

namespace SpatialOps{

//====================================================================

  // forward declaration.
  namespace structured{
    class OneDimInterpolantAssembler;
    class OneDimGradientAssembler;
    class OneDimDivergenceAssembler;
  }

//====================================================================

  // this is required for the SpatialOperator class.  It specifies
  // that we should use the InterpolantAssembler to construct
  // Interpolant objects.
  template< typename T1, typename T2 >
  struct OpAssemblerSelector< Interpolant, T1, T2 >
  {
    typedef structured::OneDimInterpolantAssembler  Assembler;
  };

//====================================================================

  template< typename T1, typename T2 >
  struct OpAssemblerSelector< Gradient, T1, T2 >
  {
    typedef structured::OneDimGradientAssembler  Assembler;
  };

//====================================================================

  template<>
  struct OpAssemblerSelector< Divergence, structured::SSurfXField, structured::SVolField >
  {
    typedef structured::OneDimDivergenceAssembler  Assembler;
  };

//====================================================================

namespace structured{

  /**
   * @todo Need to rework setting ghost cell locations so that
   *       periodic BCs work.  This affects get_x_src and get_x_dest
   */

  std::vector<double>
  get_x_src( const int nGhost,
             const std::vector<double>& xsrc,
             const std::vector<double>& xdest );


  std::vector<double>
  get_x_dest( const int nGhost,
              const std::vector<double>& xsrc,
              const std::vector<double>& xdest );

//====================================================================

  /**
   *  @class OneDimInterpolantAssembler
   *  @author James C. Sutherland
   *  @date   July, 2008
   *  @brief Provides interpolation for a one-dimensional mesh with
   *         arbitrary mesh spacing.
   */
  class OneDimInterpolantAssembler
  {
  public:

    /**
     *  \param polynomialOrder The order of polynomial to be used for
     *         the interpolant.
     *  \param x The coordinates for the source field locations (locations
     *         we want to interpolate from)
     *  \param xs The coordinates for the destination field locations
     *         (locations where we want to interpolate to)
     */
    OneDimInterpolantAssembler( const int polynomialOrder,
                                const int nGhost,
                                const std::vector<double>& xsrc,
                                const std::vector<double>& xdest );

    inline unsigned int num_nonzeros() const{ return numNonzero_; }

    int get_ncols() const;
    int get_nrows() const;

    void get_row_entries( const int irow,
                          std::vector<double> & vals,
                          std::vector<int> & ixs ) const;

    void get_ghost_cols( std::set<size_t>& ghostCols ) const;
    void get_ghost_rows( std::set<size_t>& ghostRows ) const;

  private:

    const int ncol_, nrow_;
    int polyOrder_;
    const std::vector<double> xdest_;
    const LagrangeCoefficients lagrange_;
    const int nx_;
    unsigned int numNonzero_;
  }; // class InterpolantAssembler


  //==================================================================


  /**
   *  @class OneDimInterpolantAssembler
   *  @author James C. Sutherland
   *  @date   July, 2008
   *  @brief Provides volume->surface gradients for a
   *         one-dimensional mesh with arbitrary mesh spacing.
   */
  class OneDimGradientAssembler
  {
  public:

    /**
     *  \param polynomialOrder The order of polynomial to be used for
     *         the interpolant.
     *  \param x The coordinates for the volume locations (locations
     *         we want to interpolate from)
     *  \param xs The coordinates for the surface locations (locations
     *         where we want to interpolate to)
     */
    OneDimGradientAssembler( const int polynomialOrder,
                             const int nGhost,
                             const std::vector<double>& x,
                             const std::vector<double>& xs );

    inline unsigned int num_nonzeros() const{ return numNonzero_; }

    int get_ncols() const;
    int get_nrows() const;

    void get_row_entries( const int irow,
                          std::vector<double> & vals,
                          std::vector<int> & ixs ) const;

    void get_ghost_cols( std::set<size_t>& ghostCols ) const;
    void get_ghost_rows( std::set<size_t>& ghostRows ) const;

  private:

    const int ncol_, nrow_;
    int polyOrder_;
    const std::vector<double> xdest_;
    const LagrangeCoefficients lagrange_;
    const int nx_;
    unsigned int numNonzero_;

  }; // class OneDimGradientAssembler


  //==================================================================


  class OneDimDivergenceAssembler
  {
  public:


    /**
     *  \param x The coordinates for the volume locations (locations
     *         we want to interpolate from)
     *  \param xs The coordinates for the surface locations (locations
     *         where we want to interpolate to).
     *
     *  Note that x must be at the midpoints of xs.
     */
    OneDimDivergenceAssembler( const int nGhost,
                               const std::vector<double>& x,
                               const std::vector<double>& xs );

    inline unsigned int num_nonzeros() const{ return numNonzero_; }

    int get_ncols() const{ return ncol_; }
    int get_nrows() const{ return nrow_; }

    void get_row_entries( const int irow,
                          std::vector<double> & vals,
                          std::vector<int> & ixs ) const;

    void get_ghost_cols( std::set<size_t>& ghostCols ) const;
    void get_ghost_rows( std::set<size_t>& ghostRows ) const;

  private:

    const int ncol_, nrow_, nx_;
    const unsigned int numNonzero_;
    std::vector<double> dxinv_;

  }; // class OneDimDivergenceAssembler


  //==================================================================

} // namespace structured
} // namespace SpatialOps

#endif
