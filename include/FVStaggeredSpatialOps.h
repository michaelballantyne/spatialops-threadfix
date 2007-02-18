#ifndef UT_FVStaggeredSpatialOps_h
#define UT_FVStaggeredSpatialOps_h

#include <SpatialOperator.h>
#include <SpatialField.h>

namespace SpatialOps{
namespace FVStaggeredUniform{

  enum OpType{
    CellToFace,
    FaceToCell
  };

  int rowcount( const std::vector<int> & dimExtent,
		const std::vector<int> & nghost );
  int colcount( const std::vector<int> & dimExtent,
		const std::vector<int> & nghost );

//====================================================================

  /**
   *  @class  ScratchOperator
   *  @author James C. Sutherland
   *  @date   December, 2006
   *
   *  provides a "scratch" operator for cell fields.  This may be
   *  useful when building operators from other operators & fields.
   */
  class ScratchOperator : public SpatialOperator
  {
  public:
    ScratchOperator( const std::vector<int> & dimExtent,
		     const std::vector<int> & nghostSrc,
		     const std::vector<int> & nghostDest,
		     const int entriesPerRow,
		     const Direction dir );

    ~ScratchOperator();

    void setup_matrix();

  private:

    void get_row_entries( const int irow,
			  std::vector<double> & vals,
			  std::vector<int> & ixs ) const;

    const int entriesPerRow_;
    const Direction dir_;
    int ndim_;
  };

//====================================================================

  /**
   *  @class  LinearInterpolant
   *  @author James C. Sutherland
   *  @date   December, 2006
   */
  class LinearInterpolant : public SpatialOperator
  {
  public:

    LinearInterpolant( const std::vector<int> & dimExtent,
		       const std::vector<int> & nghostSrc,
		       const std::vector<int> & nghostDest,
		       const OpType opType,
		       const Direction dir );
    
    ~LinearInterpolant();

    void setup_matrix();

  private:

    static int entries_per_row(){ return 2; }

    void get_row_entries( const int irow,
			  std::vector<double> & vals,
			  std::vector<int> & ixs ) const;

    const OpType opType_;
    const Direction dir_;
    int ndim_;
  };

//====================================================================

  /**
   *  @class  Gradient2ndOrder
   *  @author James C. Sutherland
   *  @date   December, 2006
   */
  class Gradient2ndOrder : public SpatialOperator
  {
  public:

    Gradient2ndOrder( const std::vector<double> & meshSpacing,
		      const std::vector<int> & dimExtent,
		      const std::vector<int> & nghostSrc,
		      const std::vector<int> & nghostDest,
		      const OpType opType,
		      const Direction dir );

    ~Gradient2ndOrder();

  private:

    static int entries_per_row(){ return 2; }

    void setup_matrix();

    void get_row_entries( const int irow,
			  std::vector<double> & vals,
			  std::vector<int> & ixs ) const;

    const std::vector<double> spacing_;
    const OpType opType_;
    const Direction dir_;
    int ndim_;
  };

//====================================================================

  /**
   *  @class  Divergence2ndOrder
   *  @author James C. Sutherland
   *  @date   December, 2006
   */
  class Divergence2ndOrder : public SpatialOperator
  {
  public:
    Divergence2ndOrder( const std::vector<double> & cellFaceArea,
			const double cellVolume,
			const std::vector<int> & dimExtent,
			const std::vector<int> & nghostSrc,
			const std::vector<int> & nghostDest,
			const OpType opType,
			const Direction dir );

    ~Divergence2ndOrder();

  private:

    static int entries_per_row(){ return 2; }

    void setup_matrix();
    void get_row_entries( const int irow,

			  std::vector<double> & vals,
			  std::vector<int> & ixs ) const;

    const std::vector<double> faceArea_;
    const double cellVol_;
    const OpType opType_;
    const Direction dir_;
    int ndim_;
  };

//====================================================================

} // namespace FVStaggeredUniform
} // namespace SpatialOps

#endif
