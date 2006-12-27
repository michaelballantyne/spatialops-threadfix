#ifndef UT_FVStaggeredSpatialOps_h
#define UT_FVStaggeredSpatialOps_h

#include <SpatialOperator.h>
#include <SpatialField.h>

namespace SpatialOps{
namespace FVStaggeredUniform{

//====================================================================

  /**
   *  @class  ScratchOperator
   *  @author James C. Sutherland
   *  @date   December, 2006
   */
  class ScratchOperator : public SpatialOperator
  {
  public:
    ScratchOperator( const std::vector<int> & dimExtent,
		     const int entriesPerRow,
		     const Direction dir );

    ~ScratchOperator();

    void setup_matrix();

  private:
    static int rowcount( const std::vector<int> & dimExtent, const int nghost );
    static int colcount( const std::vector<int> & dimExtent, const int nghost );
    static int my_nghost(){return 1;}

    void get_row_entries( const int irow,
			  std::vector<double> & vals,
			  std::vector<int> & ixs ) const;

    const int entriesPerRow_;
    const Direction dir_;
    int ndim_;
  };

//====================================================================

  /**
   *  @class  ScratchOperator
   *  @author James C. Sutherland
   *  @date   December, 2006
   */
  class LinearInterpolant : public SpatialOperator
  {
  public:

    LinearInterpolant( const std::vector<int> & dimExtent,
		       const Direction dir );

    ~LinearInterpolant();

    void setup_matrix();

  private:

    static int colcount( const std::vector<int>& dims );

    static int rowcount( const std::vector<int>& dims );

    static int entries_per_row(){ return 2; }

    static int my_nghost(){ return 1; }

    void get_row_entries( const int irow,
			  std::vector<double> & vals,
			  std::vector<int> & ixs ) const;

    const Direction dir_;
    int ndim_;
  };

//====================================================================

  /**
   *  @class  ScratchOperator
   *  @author James C. Sutherland
   *  @date   December, 2006
   */
  class Gradient2ndOrder : public SpatialOperator
  {
  public:

    Gradient2ndOrder( const std::vector<double> & meshSpacing,
		      const std::vector<int> & dimExtent,
		      const Direction dir );

    ~Gradient2ndOrder();

  private:

    static int rowcount( const std::vector<int> & extent );
    static int colcount( const std::vector<int> & extent );
    static int entries_per_row(){ return 2; }
    static int my_nghost(){ return 1; }

    void setup_matrix();

    void get_row_entries( const int irow,
			  std::vector<double> & vals,
			  std::vector<int> & ixs ) const;

    const std::vector<double> spacing_;
    const Direction dir_;
    int ndim_;
  };

//====================================================================

  /**
   *  @class  ScratchOperator
   *  @author James C. Sutherland
   *  @date   December, 2006
   */
  class Divergence2ndOrder : public SpatialOperator
  {
  public:
    Divergence2ndOrder( const std::vector<double> & cellFaceArea,
			const double cellVolume,
			const std::vector<int> & dimExtent,
			const Direction dir );

    ~Divergence2ndOrder();

  private:

    static int rowcount( const std::vector<int> & extent );
    static int colcount( const std::vector<int> & extent );
    static int entries_per_row(){ return 2; }
    static int my_nghost(){ return 1; }

    void setup_matrix();
    void get_row_entries( const int irow,

			  std::vector<double> & vals,
			  std::vector<int> & ixs ) const;

    const std::vector<double> faceArea_;
    const double cellVol_;
    const Direction dir_;
    int ndim_;
  };

//====================================================================

} // namespace FVStaggeredUniform
} // namespace SpatialOps

#endif
