#ifndef SpatialOps_SpatialField_h
#define SpatialOps_SpatialField_h

#include <iostream>
#include <cassert>

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/structured/MemoryWindow.h>

#ifdef SOPS_BOOST_SERIALIZATION
# include <boost/serialization/serialization.hpp>
# include <boost/serialization/split_member.hpp>
# include <boost/serialization/binary_object.hpp>
#endif

class RHS;

namespace SpatialOps{
namespace structured{

  enum StorageMode
    {
      InternalStorage,
      ExternalStorage
    };

  /**
   *  \class SpatialField
   *  \ingroup structured
   *  \ingroup fields
   *  \author James C. Sutherland
   *
   *  \brief Abstracts a field.
   *
   *  \tparam VecOps - the scheme to use in performing linear algebra
   *    operations on this field.  This must be compatible with the
   *    scheme used for operators that operate on this field.
   *
   *  \tparam FieldLocation - type traits to describe the location of
   *    this field.  On staggered meshes, this will describe the mesh
   *    this field is associated with.  It also defines whether this
   *    field is on a volume or surface.
   *
   *  \tparam GhostTraits - The ghost type traits.  Must define an
   *    enum value \c NGHOST that gives the number of ghost cells for
   *    this field.
   *
   *  \tparam T - the underlying datatype (defaults to \c double)
   *
   *  \par Related classes:
   *   - \ref MemoryWindow
   *   - \ref SpatialFieldStore
   *   - \ref SpatialOperator
   *
   *  \par Public Typedefs
   *   - \c field_type - this field's type
   *   - \c Ghost - the ghost type traits
   *   - \c Location - the location type traits
   *   - \c value_type  - the type of underlying data being stored in this SpatialField
   *   - \c iterator, \c const_iterator - iterators to the elements in this field
   *   - \c interior_iterator, \c const_interior_iterator - iterators to the interior elements in this field (excludes ghost cells).
   */
  template< typename VecOps,
            typename FieldLocation,
            typename GhostTraits,
            typename T=double >
  class SpatialField
  {
    typedef typename VecOps::VecType VecType;

    typedef SpatialField<VecOps,FieldLocation,GhostTraits,T> MyType;
    
    const MemoryWindow fieldWindow_;
    MemoryWindow interiorFieldWindow_;

    T* fieldValues_;

    const bool builtField_;

    VecOps linAlg_;
    VecType& vec_;

    inline void reset_values( const T* values );

#ifdef SOPS_BOOST_SERIALIZATION
    friend class boost::serialization::access;

    template<typename Archive>
    void save( Archive& ar, const unsigned int version ) const
    {
      ar << fieldWindow_;
      ar << interiorFieldWindow_; 
      ar << builtField_;

      const size_t npts = fieldWindow_.glob_dim(0)
                        * fieldWindow_.glob_dim(1)
                        * fieldWindow_.glob_dim(2);
      ar << npts;
      for( size_t i=0; i<npts; ++i ){
        ar << fieldValues_[i];
      }
      
    }

    template<typename Archive>
    void load( Archive& ar, const unsigned int version )
    {
      ar >> const_cast<MemoryWindow&>(fieldWindow_);
      ar >> interiorFieldWindow_;
      ar >> const_cast<bool&>(builtField_);

      size_t npts;
      ar >> npts;
      if( builtField_ ){
        fieldValues_ = new T[ npts ];
      }

      for( size_t i=0; i<npts; ++i ){
        ar >> fieldValues_[i];
      }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  public:

    typedef GhostTraits   Ghost;
    typedef FieldLocation Location;
    typedef T             AtomicT;

    typedef SpatialField<VecOps,FieldLocation,GhostTraits,T> field_type;
    typedef T value_type;

    typedef FieldIterator<T>       iterator;
    typedef FieldIterator<T>       interior_iterator;

    typedef ConstFieldIterator<T> const_iterator;
    typedef ConstFieldIterator<T> const_interior_iterator;

    /**
     *  \brief Construct a SpatialField
     *  \param window - the MemoryWindow that specifies this field
     *         including ghost cells.
     *  \param mode Storage options.  If InternalStorage then the
     *         fieldValues will be copied into an internal buffer.  If
     *         ExternalStorage then the fieldValues will be stored
     *         externally.  Efficiency suggests that ExternalStorage
     *         is best, since it will avoid excessive copies.  Safety
     *         suggests that InternalStorage is best, since it
     *         protects against memory corruption and inadvertant
     *         deletion of the field's underlying memory.
     */
    SpatialField( const MemoryWindow window,
                  T* const fieldValues,
                  const StorageMode mode = InternalStorage );

    SpatialField( const IntVec npts,
                  T* const fieldValues,
                  const StorageMode mode = InternalStorage );

    SpatialField( const SpatialField& other );

    virtual ~SpatialField();

    // warning: slow
    T& operator()( const size_t i, const size_t j, const size_t k );
    T& operator()( const size_t i, const size_t j, const size_t k ) const;

    // warning: slow
    T& operator[]( const size_t i );
    T& operator[]( const size_t i ) const;

    inline const_iterator begin() const{ return const_iterator(fieldValues_,fieldWindow_.flat_index(IntVec(0,0,0)),fieldWindow_); }
    inline       iterator begin()      { return       iterator(fieldValues_,fieldWindow_.flat_index(IntVec(0,0,0)),fieldWindow_); }

    inline const_iterator end() const;
    inline       iterator end();

    inline const_interior_iterator interior_begin() const{ return const_interior_iterator(fieldValues_,interiorFieldWindow_.flat_index(IntVec(0,0,0)),interiorFieldWindow_); }
    inline       interior_iterator interior_begin()      { return       interior_iterator(fieldValues_,interiorFieldWindow_.flat_index(IntVec(0,0,0)),interiorFieldWindow_); }

    inline const_interior_iterator interior_end() const;
    inline       interior_iterator interior_end();

    inline MyType& operator =(const MyType&);
    inline MyType& operator+=(const MyType&);
    inline MyType& operator-=(const MyType&);
    inline MyType& operator*=(const MyType&);
    inline MyType& operator/=(const MyType&);

    inline MyType& operator =(const T);
    inline MyType& operator+=(const T);
    inline MyType& operator-=(const T);
    inline MyType& operator*=(const T);
    inline MyType& operator/=(const T);

    inline SpatialField& operator= (const RHS&);  ///< Assign a RHS to this field (doesn't affect ghosts)
    inline SpatialField& operator+=(const RHS&);  ///< Add a RHS to this field (doesn't affect ghosts)
    inline SpatialField& operator-=(const RHS&);  ///< Subtract a RHS from this field (doesn't affect ghosts)

    bool operator!=(const MyType&) const;
    bool operator==(const MyType&) const;

    /**
     * @name
     * Obtain the underlying VecType object that corresponds to
     * the LinAlg strategy.
     */
    //@{
    inline       VecType & get_linalg_vec()      { assert( fieldWindow_.glob_dim() == fieldWindow_.extent() ); return vec_; }
    inline const VecType & get_linalg_vec() const{ assert( fieldWindow_.glob_dim() == fieldWindow_.extent() ); return vec_; }
    //@}

    const MemoryWindow& window_without_ghost() const{ return interiorFieldWindow_; }
    const MemoryWindow& window_with_ghost() const{ return fieldWindow_; }

  };

  //==================================================================
  //
  //                          Implementation
  //
  //==================================================================

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>::
  SpatialField( const MemoryWindow window,
                T* const fieldValues,
                const StorageMode mode )
    : fieldWindow_( window ),
      interiorFieldWindow_( window ), // reset with correct info later
      fieldValues_( (mode==ExternalStorage)
                    ? fieldValues
                    : new T[ window.glob_dim(0) * window.glob_dim(1) * window.glob_dim(2) ] ),
      builtField_( mode==InternalStorage ),
      vec_( linAlg_.setup_vector( fieldWindow_.npts(), fieldValues_ ) )
  {
    IntVec ext = window.extent();
    IntVec ofs = window.offset();
    for( size_t i=0; i<3; ++i ){
      if( ext[i]>1 ){
        ext[i] -= 2*GhostTraits::NGHOST;
        ofs[i] +=   GhostTraits::NGHOST;
      }
    }
    interiorFieldWindow_ = MemoryWindow( window.glob_dim(), ofs, ext );
    if( mode==InternalStorage )  reset_values( fieldValues );
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>::
  SpatialField( const IntVec npts,
                T* const fieldValues,
                const StorageMode mode )
    : fieldWindow_( npts, IntVec(0,0,0), npts ),
      interiorFieldWindow_( IntVec(0,0,0) ), // reset with correct info later
      fieldValues_( (mode==ExternalStorage)
                    ? fieldValues
                    : new T[npts[0]*npts[1]*npts[2]] ),
      builtField_( mode==InternalStorage ),
      vec_( linAlg_.setup_vector( fieldWindow_.npts(), fieldValues_ ) )
  {
    IntVec ext = fieldWindow_.extent();
    IntVec ofs = fieldWindow_.offset();
    for( size_t i=0; i<3; ++i ){
      if( ext[i]>1 ){
        ext[i] -= 2*GhostTraits::NGHOST;
        ofs[i] +=   GhostTraits::NGHOST;
      }
    }
    interiorFieldWindow_ = MemoryWindow( fieldWindow_.glob_dim(), ofs, ext );
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>::
  SpatialField( const SpatialField& other )
    : fieldWindow_( other.fieldWindow_ ),
      interiorFieldWindow_( other.interiorFieldWindow_ ),
      fieldValues_( new T[ fieldWindow_.glob_dim(0) * fieldWindow_.glob_dim(1) * fieldWindow_.glob_dim(2) ] ),
      builtField_( true ),
      vec_( linAlg_.setup_vector( fieldWindow_.npts(), fieldValues_ ) )
  {
    (*this) = other;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>::
  ~SpatialField()
  {
    if( builtField_ ) delete [] fieldValues_;
  }

  //------------------------------------------------------------------

  template< class VecOps, typename FieldLocation, typename GhostTraits, typename T >
  void
  SpatialField<VecOps,FieldLocation,GhostTraits,T>::
  reset_values( const T* values )
  {
    iterator ifld=begin();
    const iterator iflde=end();
    if( NULL == values ){
      for( ; ifld!=iflde; ++ifld ) *ifld = 0.0;
    }
    else{
      for( ; ifld!=iflde; ++ifld, ++values )  *ifld = *values;
    }
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>::const_iterator
  SpatialField<VecOps,Location,GhostTraits,T>::end() const
  {
    IntVec ijk = fieldWindow_.extent();
    for( size_t i=0; i<3; ++i ) ijk[i] -= 1;
    const size_t n = fieldWindow_.flat_index( ijk );
    const_iterator i(fieldValues_, n, fieldWindow_);
    return ++i;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>::iterator
  SpatialField<VecOps,Location,GhostTraits,T>::end()
  {
    IntVec ijk = fieldWindow_.extent();
    for( size_t i=0; i<3; ++i ) ijk[i] -= 1;
    const size_t n = fieldWindow_.flat_index( ijk );
    iterator i(fieldValues_, n, fieldWindow_);
    return ++i;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>::const_interior_iterator
  SpatialField<VecOps,Location,GhostTraits,T>::interior_end() const
  {
    IntVec ijk = interiorFieldWindow_.extent();
    for( size_t i=0; i<3; ++i ) ijk[i] -= 1;
    const_interior_iterator i( fieldValues_, interiorFieldWindow_.flat_index( ijk ), interiorFieldWindow_ );
    return ++i;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>::interior_iterator
  SpatialField<VecOps,Location,GhostTraits,T>::interior_end()
  {
    IntVec ijk = interiorFieldWindow_.extent();
    for( size_t i=0; i<3; ++i ) ijk[i] -= 1;
    interior_iterator i( fieldValues_, interiorFieldWindow_.flat_index( ijk ), interiorFieldWindow_ );
    return ++i;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  T&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator()( const size_t i, const size_t j, const size_t k )
  {
#   ifndef NDEBUG
    assert( i < fieldWindow_.extent(0) );
    assert( j < fieldWindow_.extent(1) );
    assert( k < fieldWindow_.extent(2) );
#   endif
    return fieldValues_[ fieldWindow_.flat_index(IntVec(i,j,k)) ];
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  T&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator()( const size_t i, const size_t j, const size_t k ) const
  {
#   ifndef NDEBUG
    assert( i < fieldWindow_.extent(0) );
    assert( j < fieldWindow_.extent(1) );
    assert( k < fieldWindow_.extent(2) );
#   endif
    return fieldValues_[ fieldWindow_.flat_index(IntVec(i,j,k)) ];
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  T&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator[]( const size_t i )
  {
    return fieldValues_[ fieldWindow_.flat_index( fieldWindow_.ijk_index_from_local(i) ) ];
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  T&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator[]( const size_t i ) const
  {
    return fieldValues_[ fieldWindow_.flat_index( fieldWindow_.ijk_index_from_local(i) ) ];
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator=(const MyType& other )
  {
    const_iterator iother=other.begin();
    const iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld, ++iother ){
      *ifld = *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator+=(const MyType& other )
  {
    const_iterator iother=other.begin();
    const iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld, ++iother ){
      *ifld += *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator-=(const MyType& other )
  {
    const_iterator iother=other.begin();
    const iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld, ++iother ){
      *ifld -= *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator*=(const MyType& other )
  {
    const_iterator iother=other.begin();
    const iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld, ++iother ){
      *ifld *= *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator/=(const MyType& other )
  {
    const_iterator iother=other.begin();
    const iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld, ++iother ){
      *ifld /= *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  bool
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator!=(const MyType& other) const
  {
    const_iterator iother=other.begin();
    const_iterator iend=this->end();
    for( const_iterator ifld=this->begin(); ifld!=iend; ++ifld, ++iother ){
      if( *ifld == *iother ) return false;
    }
    return true;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  bool
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator==(const MyType& other) const
  {
    const_iterator iother=other.begin();
    const_iterator iend=this->end();
    for( const_iterator ifld=this->begin(); ifld!=iend; ++ifld, ++iother ){
      if( *ifld != *iother ) return false;
    }
    return true;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator=(const T a)
  {
    const iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld ) *ifld = a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator+=(const T a)
  {
    const iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld ) *ifld += a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator-=(const T a)
  {
    const iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld ) *ifld -= a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator*=(const T a)
  {
    const iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld ) *ifld *= a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator/=(const T a)
  {
    const iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld ) *ifld /= a;
    return *this;
  }

  //------------------------------------------------------------------

} // namespace structured

//#include <spatialops/FieldOperations.h>

} // namespace SpatialOps

#endif // SpatialOps_SpatialField_h
