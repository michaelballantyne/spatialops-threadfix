#include <map>
#include <iostream>

#include <LinAlgTrilinos.h>

//--------------------------------------
// Trilinos Includes
#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_SerialComm.h>
#include <EpetraExt_MatrixMatrix.h>
// Trilinos Includes
//--------------------------------------

namespace SpatialOps{


  //====================================================================



  /**
   *  @class  MapFactory
   *  @author James C. Sutherland
   *  @date   December, 2006
   */
  class MapFactory
  {

  public:

    static MapFactory& self()
    {
      static MapFactory s;
      return s;
    }

    const Epetra_LocalMap & get_map( const int npts )
    {
      InfoEpetraMap::const_iterator ii = map_.find( npts );
      if( ii == map_.end() ){
        Epetra_LocalMap * myMap = new Epetra_LocalMap( npts, 0, *com_ );
        map_.insert( std::make_pair(npts,myMap) );
        return *myMap;
      }
      return *(ii->second);
    }

  private:
    MapFactory()
      : com_( new Epetra_SerialComm() )
    {}

    ~MapFactory()
    {
      for( InfoEpetraMap::iterator ii=map_.begin(); ii!=map_.end(); ++ii ){
        delete ii->second;
      }
    }

    Epetra_SerialComm * const com_;
    typedef std::map<int,Epetra_LocalMap*> InfoEpetraMap;
    InfoEpetraMap map_;
   };


  //==================================================================












  //--------------------------------------------------------------------
  LinAlgTrilinos::LinAlgTrilinos()
  {
    mat_ = NULL;
    vec_ = NULL;
  }
  //------------------------------------------------------------------
  LinAlgTrilinos::~LinAlgTrilinos()
  {
    destroy_matrix();
    destroy_vector();
  }
  //------------------------------------------------------------------
  void
  LinAlgTrilinos::destroy_matrix()
  {
    if( mat_ != NULL ){
      delete mat_;
      mat_ = NULL;
    }
  }
  //------------------------------------------------------------------
  void
  LinAlgTrilinos::destroy_vector()
  {
    if( vec_ != NULL ){
      delete vec_;
      vec_=NULL;
    }
  }
  //------------------------------------------------------------------
  Epetra_Vector&
  LinAlgTrilinos::setup_vector( const int npts,
                                double* const fieldValues )
  {
    const Epetra_LocalMap & epetraMap = MapFactory::self().get_map( npts );
    vec_ = new Epetra_Vector( View, epetraMap, fieldValues );
    return *vec_;
  }
  //------------------------------------------------------------------
  Epetra_CrsMatrix&
  LinAlgTrilinos::setup_matrix( const int nrows,
                                const int ncols,
                                const int entriesPerRow )
  {
    rowMap_ = &MapFactory::self().get_map( nrows );
    colMap_ = &MapFactory::self().get_map( ncols ),
    mat_ = new Epetra_CrsMatrix( Copy,
                                 *rowMap_,
                                 *colMap_,
                                 entriesPerRow,
                                 true );
    return *mat_;
  }
  //------------------------------------------------------------------
  void
  LinAlgTrilinos::insert_row_values( const int rownum,
                                     std::vector<double> & rowValues,
                                     std::vector<int> & rowIndices )
  {
    using namespace std;

    assert( rowValues.size() == rowIndices.size() );

    const int flag = mat_->InsertMyValues( rownum,
                                           rowValues.size(),
                                           &rowValues[0],
                                           &rowIndices[0] );
    if( flag!=0 ){
      cout << flag << endl
           << "Error inserting values into row: " << rownum << endl
           << " nonzero column indices: [ ";
      for( vector<int>::const_iterator ii=rowIndices.begin(); ii!=rowIndices.end(); ++ii ){
        cout << *ii << " ";
      }
      cout << "]" << endl << endl;
    }
    assert( flag==0 );
  }
  //------------------------------------------------------------------
  void
  LinAlgTrilinos::multiply( const Epetra_CrsMatrix & B,
                            Epetra_CrsMatrix & C ) const
  {
    const bool useTranspose = false;
    const int flag =
      EpetraExt::MatrixMatrix::Multiply( *mat_, useTranspose,
                                         B, useTranspose,
                                         C );
    if( flag!=0 )
      std::cout << std::endl
                << "ERROR!  Flag=" << flag
                << " returned from EpetraExt::MatrixMatrix::Multiply()." << std::endl
                << "        This likely indicates incompatible matrices for multiplication." << std::endl
                << "        Check matrix sparsity patterns and dimensions for compatibility."
                << std::endl << std::endl;
  }
  //------------------------------------------------------------------
  void
  LinAlgTrilinos::multiply( const Epetra_Vector& src, Epetra_Vector& dest ) const
  {
    mat_->Multiply( false, src, dest );
  }
  //------------------------------------------------------------------
  void
  LinAlgTrilinos::left_scale( const VecType& vec )
  {
#ifdef NDEBUG
    mat_->LeftScale( vec );
#else
    const int flag = mat_->LeftScale( vec );
    assert( flag == 0 );
#endif
  }
  //------------------------------------------------------------------
  void
  LinAlgTrilinos::right_scale( const VecType& vec )
  { 
#ifdef NDEBUG
    mat_->RightScale( vec );
#else
    const int flag = mat_->RightScale( vec );
    assert( flag == 0 );
#endif
  }
  //------------------------------------------------------------------
  void
  LinAlgTrilinos::finalize()
  {
    mat_->FillComplete( *colMap_, *rowMap_ );
    mat_->OptimizeStorage();
  }
  //--------------------------------------------------------------------
  void
  LinAlgTrilinos::reset_entries( const double val )
  {
    mat_->PutScalar(val);
    return;
  }
  //--------------------------------------------------------------------
  void
  LinAlgTrilinos::print_vec( std::ostream& s ) const
  {
    vec_->Print(s);
  }
  //--------------------------------------------------------------------
  void
  LinAlgTrilinos::print_mat( std::ostream& s ) const
  {
    mat_->Print(s);
  }
  //--------------------------------------------------------------------





  //====================================================================

}
