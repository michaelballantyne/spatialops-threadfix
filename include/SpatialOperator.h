#ifndef UT_SpatialOperator_h
#define UT_SpatialOperator_h

#include <vector>
#include <map>
#include <string>

//-------------------------------------
// trilinos class forward declarations:
class Epetra_CrsMatrix;
class Epetra_LocalMap;
class Epetra_SerialComm;
//-------------------------------------


namespace SpatialOps{


// forward declaration
class SpatialField;
class SpatialOperator;

//====================================================================


/**
 *  @class  SpatialOpDatabase
 *  @author James C. Sutherland
 *  @date   December, 2006
 *
 *  Factory for SpatialOperator objects.  These objects are registered
 *  here (created externally, ownership is transferred) and can be
 *  recalled for later use.
 *
 *  Note that one can have multiple operators registered for a given
 *  type, and activate them as needed.  This allows the potential for
 *  dynamic operator switching.
 *
 */
class SpatialOpDatabase
{
public:

  enum OperatorType{

    CELL_DIVERGENCE_X,     FACE_DIVERGENCE_X,
    CELL_DIVERGENCE_Y,     FACE_DIVERGENCE_Y,
    CELL_DIVERGENCE_Z,     FACE_DIVERGENCE_Z,

    CELL_GRADIENT_X,       FACE_GRADIENT_X,
    CELL_GRADIENT_Y,       FACE_GRADIENT_Y,
    CELL_GRADIENT_Z,       FACE_GRADIENT_Z,

    CELL_INTERPOLANT_X,    FACE_INTERPOLANT_X,
    CELL_INTERPOLANT_Y,    FACE_INTERPOLANT_Y,
    CELL_INTERPOLANT_Z,    FACE_INTERPOLANT_Z,

    CELL_SCRATCH_X,        FACE_SCRATCH_X,
    CELL_SCRATCH_Y,        FACE_SCRATCH_Y,
    CELL_SCRATCH_Z,        FACE_SCRATCH_Z

  };

  static SpatialOpDatabase& self();

  /**
   *  Registers a new operator.
   *
   *  @param opType : The type of operator.

   *  @param op : The operator itself.  Ownership is transfered.  This
   *  must be a heap-allocated object (build via "new")
   *
   *  @param opName : The name for this operator.  Must be a unique
   *  name.  Duplicate names will result in an exception.
   *
   *  @param makeDefault : [true] By default, registration of a new
   *  operator makes it the default operator.  If this flag is "false"
   *  then it will not replace the current default operator, unless
   *  one does not yet exist.
   */
  void register_new_operator( const OperatorType opType,
			      SpatialOperator * const op,
			      const std::string& opName,
			      const bool makeDefault = true );

  /**
   *  Reset the default operator to the one with the given name.
   */
  void set_default_operator( const OperatorType opType,
			     const std::string & opName,
			     const std::vector<int> & nxyz,
			     const std::vector<int> & nghostSrc,
			     const std::vector<int> & nghostDest );

  /**
   *  Obtain the spatial operator with the given type and shape.
   *  Throws an exception if no match is found.
   *
   *  This pointer reference can change if the default operator for
   *  this type is changed via a call to
   *  <code>set_default_operator</code> or
   *  <code>register_new_operator</code>.
   */
  SpatialOperator*& retrieve_operator( const OperatorType opType,
				       const std::vector<int> & nxyz,
				       const std::vector<int> & nghostSrc,
				       const std::vector<int> & nghostDest );

  /**
   *  Obtain the spatial operator with the given name and shape.
   *  Throws an exception if no match is found.
   *
   *  This returns a pointer reference that will never change.
   */
  SpatialOperator*& retrieve_operator( const std::string & opName,
				       const std::vector<int> & nxyz,
				       const std::vector<int> & nghostSrc,
				       const std::vector<int> & nghostDest );

private:

  SpatialOpDatabase();
  ~SpatialOpDatabase();

  SpatialOpDatabase(const SpatialOpDatabase&);
  SpatialOpDatabase&  operator=(const SpatialOpDatabase&);

  struct Shape{
    Shape( const std::vector<int> & extent,
	   const std::vector<int> ghostSrc,
	   const std::vector<int> ghostDest );
    std::vector<int> nxyz, ngS, ngD;
    bool operator ==( const Shape& s ) const;
    bool operator < ( const Shape& s ) const;
  };

  typedef std::map< Shape,SpatialOperator* > ShapeOpMap;

  typedef std::map< OperatorType, ShapeOpMap > TypeShapeMap;
  typedef std::map< std::string,  ShapeOpMap > NameShapeMap;

  NameShapeMap nameMap_;
  TypeShapeMap typeMap_;

};


//====================================================================


/**
 *  @class  SpatialOperator
 *  @author James C. Sutherland
 *  @date   December, 2006
 *
 *  Provides support for discrete spatial operators.
 *
 *  This is intended to be used as a base class only, and provides
 *  most of the functionality required.  Derived classes mainly are
 *  responsible for populating each row of the matrix.
 *
 *  Several rules apply:
 *
 *   - Application of a SpatialOperator must not involve parallel
 *   communication.  This is to ensure efficiency.
 *
 *   - 
 *   
 *
 */
class SpatialOperator
{

public:

  enum Dimension{
    XDIM = 0,
    YDIM = 1,
    ZDIM = 2
  };

  enum Side{
    MINUS = 0,
    PLUS  = 1
  };


  /**
   *  Construct a SpatialOperator.
   *
   *  @param nrows: The number of rows in this matrix
   *
   *  @param ncols : The number of columns in this matrix
   *
   *  @param nghost : The number of ghost cells on the (-) and (+)
   *  side of the patch in each coordinate direction.  For example,
   *  [ngxl, ngxr, ngyl, ngyr, ngzl, ngzr]
   *
   *  @param entriesPerRow : The number of nonzero entries on each row of this matrix operator.
   *
   *  @param extent : The number of points in each direction
   *  (excluding ghost cells) for the domain extent.
   */
  SpatialOperator( const int nrows,
		   const int ncols,
		   const std::vector<int> & nghostSrc,
		   const std::vector<int> & nghostDest,
		   const int entriesPerRow,
		   const std::vector<int> & extent );
		   
  virtual ~SpatialOperator();

  /**
   *  Apply this SpatialOperator to a field.  The source and
   *  destination fields must have compatible dimension and ghosting
   *  for use with this operator.
   */
  void apply( const SpatialField & srcField,
	      SpatialField & destField ) const;

  /**
   *  performs matrix-matrix multiplication, resulting in another
   *  SpatialOperator
   */
  void apply( const SpatialOperator & srcOp,
	      SpatialOperator & destOp ) const;


  /**
   *  returns true if this operator is ready for use
   */
  inline bool ready() const{ return isFinalized_; }

  
  SpatialOperator& operator  = ( const SpatialOperator& );

  SpatialOperator& operator += ( const SpatialOperator& );
  SpatialOperator& operator -= ( const SpatialOperator& );

  /** sum the given field into the diagonal */
  SpatialOperator& operator += ( const SpatialField& );

  /** subtract the given field from the diagonal */
  SpatialOperator& operator -= ( const SpatialField& );


  /**
   *  Scale the matrix such that A(i,j) = x(i)*A(i,j) where i denotes
   *  the row number of A and j denotes the column number of A.
   */
  void  left_scale( const SpatialField& );

  /**
   *  Scale the matrix such that A(i,j) = x(j)*A(i,j) where i denotes
   *  the global row number of A and j denotes the global column
   *  number of A.
   */
  void right_scale( const SpatialField& );


  /** reset the coefficients in the matrix */
  void reset_entries( const double val = 0 );

  //@{ /** Obtain a reference to the underlying Epetra_CrsMatrix object */

        Epetra_CrsMatrix & epetra_mat();
  const Epetra_CrsMatrix & epetra_mat() const;

  //}@

  /** Obtain the number of rows in this operator */
  inline int nrows() const{ return nrows_; }

  /** Obtain the number of columns in this operator */
  inline int ncols() const{ return ncols_; }



  inline const std::vector<int>& nghost_src()  const{return nghostSrc_ ;}
  inline const std::vector<int>& nghost_dest() const{return nghostDest_;}

  inline const int nghost_src( const Dimension dim, const Side side ) const
  {
    const int ix = 2*int(dim) + int(side);
    return nghostSrc_[ix];
  }

  inline const int nghost_dest( const Dimension dim, const Side side ) const
  {
    const int ix = 2*int(dim) + int(side);
    return nghostDest_[ix];
  }


  struct IndexTriplet{
    int index[3];
    inline int& operator[](const int i){return index[i];}
  };

  /** if provided, the IndexTriplet is populated with the interior ijk index for non-ghost entries.*/
  bool is_col_ghost( const int colnum, IndexTriplet* const ix=NULL ) const;
  bool is_row_ghost( const int rownum, IndexTriplet* const ix=NULL ) const;

  const std::vector<int> & get_extent() const{ return extent_; }

  void Print( std::ostream & ) const;

protected:


  /**
   *  Derived classes must implement this and should call it from
   *  their constructors.  After setting the matrix coefficients, call
   *  SpatialOperator::finalize()
   */
  virtual void setup_matrix() = 0;


  /**
   *  Check compatibility of this operator applied on the one provided
   *  as an argument.  The base class implementation simply performs
   *  dimensional compatibility checks.  Derived classes may
   *  specialize this method further.
   */
  virtual bool compatibility_check( const SpatialOperator& op,
				    const bool isResultOp ) const;


  enum FieldType{ SOURCE_FIELD, DEST_FIELD };

  /**
   * Check compatability of a field with this operator.  The base
   * class method simply checks for dimensional compatibility.
   * Derived classes may/should specialize this method to ensure
   * proper ghost availability.
  */
  virtual bool compatibility_check( const SpatialField & field,
				    const FieldType ) const;

  /** Insert a row into the matrix */
  void insert_row_entry( const int rownum,
			 std::vector<double> & rowValues,
			 std::vector<int> & rowIndices );

  /** Sum a row into the matrix */
  void sum_into_row( const int rownum,
		     std::vector<double> & rowValues,
		     std::vector<int> & rowIndices );

  /** Finalize the assembly of the matrix.  Must be called prior to using it. */
  void finalize();

  const Epetra_LocalMap *rowMap_, *colMap_;
  const std::vector<int> extent_;

private:

  const int nrows_, ncols_, entriesPerRow_;
  const std::vector<int> nghostSrc_, nghostDest_;
  bool isFinalized_;
  Epetra_CrsMatrix * mat_;

};


//====================================================================

/**
 *  @class  MapFactory
 *  @author James C. Sutherland
 *  @date   December, 2006
 */
class MapFactory
{
public:
  static MapFactory& self();
  
  const Epetra_LocalMap & get_map( const int npts );

private:
  MapFactory();
  ~MapFactory();

  Epetra_SerialComm * const com_;
  typedef std::map<int,Epetra_LocalMap*> InfoEpetraMap;
  InfoEpetraMap map_;
};


//====================================================================


} // namespace SpatialOps


#endif
