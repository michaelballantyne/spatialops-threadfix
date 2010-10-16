#ifndef SpatialOps_SpatialField_h
#define SpatialOps_SpatialField_h

#include <spatialops/structured/MemoryWindow.h>

namespace SpatialOps{
namespace structured{

  enum StorageMode
    {
      InternalStorage,
      ExternalStorage
    };

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

    T* const fieldValues_;

    const bool builtField_;

    VecOps linAlg_;
    VecType& vec_;

    inline void reset_values( const T* values );

  public:

    typedef GhostTraits Ghost;
    typedef FieldLocation Location;

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

    virtual ~SpatialField();

    // warning: slow
    T& operator()( const int i, const int j, const int k );
    T  operator()( const int i, const int j, const int k ) const;

    // warning: slow
    T& operator[]( const size_t i );
    T  operator[]( const size_t i ) const;

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

    bool operator!=(const MyType&) const;
    bool operator==(const MyType&) const;

    /**
     * @name
     * Obtain the underlying VecType object that corresponds to
     * the LinAlg strategy.
     */
    //@{
    inline       VecType & get_linalg_vec()      { return vec_; }
    inline const VecType & get_linalg_vec() const{ return vec_; }
    //@}

    const MemoryWindow& window() const{ return fieldWindow_; }

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
      interiorFieldWindow_( window ),
      fieldValues_( (mode==ExternalStorage)
                    ? fieldValues
                    : new T[ window.glob_dim(0) * window.glob_dim(1) * window.glob_dim(2) ] ),
      builtField_( mode==InternalStorage ),
      vec_( linAlg_.setup_vector( fieldWindow_.npts(), fieldValues_ ) )
  {
    for( size_t i=0; i<3; ++i ){
      if( window.extent(i)>1 ){
        interiorFieldWindow_.extent(i) -= 2*GhostTraits::NGHOST;
        interiorFieldWindow_.offset(i) +=   GhostTraits::NGHOST;
      }
    }
    if( mode==InternalStorage )  reset_values( fieldValues );
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>::
  SpatialField( const IntVec npts,
                T* const fieldValues,
                const StorageMode mode )
    : fieldWindow_( npts, IntVec(0,0,0), npts ),
      interiorFieldWindow_( npts ),
      fieldValues_( (mode==ExternalStorage)
                    ? fieldValues
                    : new T[npts[0]*npts[1]*npts[2]] ),
      builtField_( mode==InternalStorage ),
      vec_( linAlg_.setup_vector( fieldWindow_.npts(), fieldValues_ ) )
  {
    for( size_t i=0; i<3; ++i ){
      if( npts[i]>1 ){
        interiorFieldWindow_.extent(i) -= 2*GhostTraits::NGHOST;
        interiorFieldWindow_.offset(i) +=   GhostTraits::NGHOST;
      }
    }
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
    const size_t n = 1+fieldWindow_.flat_index( ijk );
    return const_iterator(fieldValues_, n, fieldWindow_);
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>::iterator
  SpatialField<VecOps,Location,GhostTraits,T>::end()
  {
    IntVec ijk = fieldWindow_.extent();
    for( size_t i=0; i<3; ++i ) ijk[i] -= 1;
    const size_t n = 1+fieldWindow_.flat_index( ijk );
    return iterator(fieldValues_, n, fieldWindow_);
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>::const_interior_iterator
  SpatialField<VecOps,Location,GhostTraits,T>::interior_end() const
  {
    IntVec ijk = interiorFieldWindow_.extent();
    for( size_t i=0; i<3; ++i ) ijk[i] -= 1;
    const size_t n = 1+interiorFieldWindow_.flat_index( ijk );
    return const_interior_iterator( fieldValues_, n, interiorFieldWindow_ );
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>::interior_iterator
  SpatialField<VecOps,Location,GhostTraits,T>::interior_end()
  {
    IntVec ijk = interiorFieldWindow_.extent();
    for( size_t i=0; i<3; ++i ) ijk[i] -= 1;
    const size_t n = 1+interiorFieldWindow_.flat_index( ijk );
    return interior_iterator( fieldValues_, n, interiorFieldWindow_ );
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  T&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator()( const int i, const int j, const int k )
  {
    assert( i < fieldWindow_.extent(0) );
    assert( j < fieldWindow_.extent(1) );
    assert( k < fieldWindow_.extent(2) );
    return fieldValues_[ fieldWindow_.flat_index(IntVec(i,j,k)) ];
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  T
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator()( const int i, const int j, const int k ) const
  {
    assert( i < fieldWindow_.extent(0) );
    assert( j < fieldWindow_.extent(1) );
    assert( k < fieldWindow_.extent(2) );
    return fieldValues_[ fieldWindow_.flat_index(IntVec(i,j,k)) ];
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  T&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator[]( const size_t i )
  {
    return fieldValues_[ fieldWindow_.flat_index( fieldWindow_.ijk_index(i) ) ];
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  T
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator[]( const size_t i ) const
  {
    return fieldValues_[ fieldWindow_.flat_index( fieldWindow_.ijk_index(i) ) ];
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
} // namespace SpatialOps

#endif // SpatialOps_SpatialField_h
