#ifndef UT_SpatialOperator_h
#define UT_SpatialOperator_h

#include <vector>
#include <map>

//-------------------------------------
// trilinos class forward declarations:
class Epetra_CrsMatrix;
class Epetra_LocalMap;
class Epetra_SerialComm;
//-------------------------------------


namespace SpatialOps{


// forward declaration
class SpatialField;


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
 *   communication.
 *
 *   - 
 *   .
 *
 */
class SpatialOperator
{

public:

  /**
   *  Construct a SpatialOperator.
   *
   *  @param nrows: The number of rows in this matrix
   *
   *  @param ncols : The number of columns in this matrix
   *
   *  @param entriesPerRow : The number of nonzero entries on each row of this matrix operator.
   *
   *  @param extent : The number of points in each direction
   *  (excluding ghost cells) for the domain extent.
   */
  SpatialOperator( const int nrows,
		   const int ncols,
		   const int nghost,
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


  /** zero out the coefficients in the matrix */
  void zero_entries();

  //@{ /** Obtain the underlying Epetra_CrsMatrix object */

        Epetra_CrsMatrix & epetra_mat();
  const Epetra_CrsMatrix & epetra_mat() const;

  //}@

  /** Obtain the number of rows in this operator */
  inline int nrows() const{ return nrows_; }

  /** Obtain the number of columns in this operator */
  inline int ncols() const{ return ncols_; }

  inline int nghost() const{ return nghost_; }

  const std::vector<int> & get_extent() const{ return extent_; }

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
  virtual bool compatibility_check( const SpatialOperator& op  ) const;


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

  const int nrows_, ncols_, nghost_, entriesPerRow_;
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
