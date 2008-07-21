#ifndef FVStaggeredOneDimStretch_h
#define FVStaggeredOneDimStretch_h

#include <FVStaggeredTypes.h>
#include <LagrangePoly.h>
#include <SpatialField.h>
#include <SpatialOperator.h>
#include <FVStaggeredTools.h>


#include <vector>
#include <set>

namespace SpatialOps{

//====================================================================

  // forward declaration.
  namespace FVStaggered{
    class OneDimInterpolantAssembler;
    class OneDimGradientAssembler;
    class OneDimDivergenceAssembler;
  }

//====================================================================

  // this is required for the SpatialOperator class.  It specifies
  // that we should use the InterpolantAssembler to construct
  // Interpolant objects.
  template< typename SVolField, typename SSurfXField >
  struct OpAssemblerSelector< FVStaggered::Interpolant, SVolField, SSurfXField >
  {
    typedef FVStaggered::OneDimInterpolantAssembler  Assembler;
  };

//====================================================================

  template< typename SVolField, typename SSurfXField >
  struct OpAssemblerSelector< FVStaggered::Gradient, SVolField, SSurfXField >
  {
    typedef FVStaggered::OneDimGradientAssembler  Assembler;
  };

//====================================================================

  template< typename SVolField, typename SSurfXField >
  struct OpAssemblerSelector< FVStaggered::Divergence, SVolField, SSurfXField >
  {
    typedef FVStaggered::OneDimDivergenceAssembler  Assembler;
  };

//====================================================================

namespace FVStaggered{

//====================================================================

  /**
   *  @class OneDimInterpolantAssembler
   *  @author James C. Sutherland
   *  @date   July, 2008
   *  @brief Provides volume->surface interpolation for a
   *         one-dimensional mesh with arbitrary mesh spacing.
   */
  class OneDimInterpolantAssembler
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
    OneDimInterpolantAssembler( const int polynomialOrder,
                                const std::vector<double>& x,
                                const std::vector<double>& xs );

    inline unsigned int num_nonzeros() const{ return numNonzero_; }

    int get_ncols() const;
    int get_nrows() const;

    void get_row_entries( const int irow,
                          std::vector<double> & vals,
                          std::vector<int> & ixs ) const;

    void get_ghost_cols( std::set<int>& ghostCols ) const;
    void get_ghost_rows( std::set<int>& ghostRows ) const;

    static std::vector<double>
    get_x( const std::vector<double>& x,
           const std::vector<double>& xs );

    static std::vector<double>
    get_xs( const std::vector<double>& x,
            const std::vector<double>& xs );

  private:

    int polyOrder_;
    const std::vector<double> xsurf_;
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
                             const std::vector<double>& x,
                             const std::vector<double>& xs );

    inline unsigned int num_nonzeros() const{ return numNonzero_; }

    int get_ncols() const;
    int get_nrows() const;

    void get_row_entries( const int irow,
                          std::vector<double> & vals,
                          std::vector<int> & ixs ) const;

    void get_ghost_cols( std::set<int>& ghostCols ) const;
    void get_ghost_rows( std::set<int>& ghostRows ) const;

    inline std::vector<double>
    get_x( const std::vector<double>& x,
           const std::vector<double>& xs ) const{ return OneDimInterpolantAssembler::get_x(x,xs); }

    inline std::vector<double>
    get_xs( const std::vector<double>& x,
            const std::vector<double>& xs ) const{ return OneDimInterpolantAssembler::get_xs(x,xs); }

  private:

    int polyOrder_;
    const std::vector<double> xsurf_;
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
    OneDimDivergenceAssembler( const std::vector<double>& x,
                               const std::vector<double>& xs );

    inline unsigned int num_nonzeros() const{ return numNonzero_; }

    int get_ncols() const{ return nx_+3; }
    int get_nrows() const{ return nx_+2; }

    void get_row_entries( const int irow,
                          std::vector<double> & vals,
                          std::vector<int> & ixs ) const;

    void get_ghost_cols( std::set<int>& ghostCols ) const;
    void get_ghost_rows( std::set<int>& ghostRows ) const;

  private:

    const int nx_;
    const unsigned int numNonzero_;
    std::vector<double> dxinv_;

  }; // class OneDimDivergenceAssembler


  //==================================================================


} // namespace FVStaggered
} // namespace SpatialOps

#endif
