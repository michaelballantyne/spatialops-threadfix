#ifndef UT_SpatialField_h
#define UT_SpatialField_h

#include <vector>

//-------------------------------------
// trilinos class forward declarations:
class Epetra_Vector;
//-------------------------------------


class RHS; // forward declaration



namespace SpatialOps{

enum Direction{ X_DIR=0, Y_DIR=1, Z_DIR=2 };


/**
 *  @class SpatialField
 *  @author James C. Sutherland
 *  @date   December, 2006
 *
 *  Base class for SpatialFields defined on a logically rectangular
 *  domain.
 */
class SpatialField
{
public:

  enum StorageMode{
    InternalStorage,
    ExternalStorage
  };

  /**
   *  Construct a SpatialField.
   *
   *  @param fieldDims: The number of points (excluding any ghost
   *  cells) for the domain in each of the three ordinal directions.
   *
   *  @param fieldValues : Pointer to the field values.  Behavior is
   *  dictated by the choice of <code>StorageMode</code>
   *
   *  @param mode : Storage options.  If InternalStorage then the
   *  fieldValues will be copied into an internal buffer.  If
   *  ExternalStorage then the fieldValues will be stored externally.
   *  Efficiency suggests that ExternalStorage is best, since it will
   *  avoid excessive copies.
   */
  SpatialField( const std::vector<int> & fieldDims,
		const int nghost,
		double * const fieldValues,
		const StorageMode mode );

  virtual ~SpatialField();


  /**
   *  @param npts : number of points (including ghost cells)
   *  @param values : array of values to overwrite with.
   */
  void reset_values( const int npts,
		     const double* const values );


  //@{ /** Obtain the underlying Trilinos Epetra_Vector object. */

        Epetra_Vector& epetra_vec();
  const Epetra_Vector& epetra_vec() const;

  //}@


  //@{  /** Operators for SpatialField objects */

  //@{
  /**
   *  These operators preserve/update ghost cell information. 
   */
  SpatialField& operator  =(const SpatialField&);
  SpatialField& operator +=(const SpatialField&);
  SpatialField& operator -=(const SpatialField&);
  SpatialField& operator *=(const SpatialField&);
  SpatialField& operator /=(const SpatialField&);

  SpatialField& operator  =(const double);
  SpatialField& operator *=(const double);
  SpatialField& operator +=(const double);
  //}@

  //@{
  /**
   *  Note that these operators will invalidate ghost cell information
   *  on this SpatialField, since RHS fields don't have ghost
   *  information.
   */
  SpatialField& operator  =(const RHS&);
  SpatialField& operator +=(const RHS&);
  SpatialField& operator -=(const RHS&);
  //}@

  //}@

  
  /** get the total number of points (including ghost layers) */
  inline int get_ntotal() const{ return npts_; }

  inline int nghost() const{ return nghost_; }

  inline const std::vector<int>& get_extent() const{return extent_;}

  /** obtain a pointer to the underlying field - this should be used carefully! */
  inline       double* get_ptr()      { return fieldValues_; }
  inline const double* get_ptr() const{return fieldValues_; }

  void Print( std::ostream& ) const;

protected:

  virtual bool consistency_check( const SpatialField& s ) const;

private:

  static int get_npts( const std::vector<int> & extent,
		       const int nghost );

  const std::vector<int> extent_;
  const int npts_;
  const int nghost_;
  const StorageMode storageMode_;
  double * const fieldValues_;

  Epetra_Vector * vec_;

  SpatialField( const SpatialField& );
  SpatialField();

};

//====================================================================

} // namespace SpatialOps

#endif
