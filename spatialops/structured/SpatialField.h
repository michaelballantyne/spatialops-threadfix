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

  public:

    typedef GhostTraits Ghost;
    typedef FieldLocation Location;

    typedef FieldIterator<      T*>       iterator;
    typedef FieldIterator<const T*> const_iterator;

    typedef FieldIterator<      T*>       interior_iterator;
    typedef FieldIterator<const T*> const_interior_iterator;

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

    inline const_iterator begin() const{ return const_iterator(fieldValues_+fieldWindow_.flat_index(IntVec(0,0,0)),fieldWindow_); }
    inline       iterator begin()      { return       iterator(fieldValues_+fieldWindow_.flat_index(IntVec(0,0,0)),fieldWindow_); }

    inline const_iterator end() const{ return const_iterator(fieldValues_ + fieldWindow_.npts(), fieldWindow_); }
    inline       iterator end()      { return       iterator(fieldValues_ + fieldWindow_.npts(), fieldWindow_); }

    inline const_interior_iterator interior_begin() const{ return const_interior_iterator(fieldValues_+interiorFieldWindow_.flat_index(IntVec(0,0,0)),interiorFieldWindow_); }
    inline       interior_iterator interior_begin()      { return       interior_iterator(fieldValues_+interiorFieldWindow_.flat_index(IntVec(0,0,0)),interiorFieldWindow_); }

    inline const_interior_iterator interior_end() const{ return const_interior_iterator(fieldValues_ + interiorFieldWindow_.npts(), interiorFieldWindow_); }
    inline       interior_iterator interior_end()      { return       interior_iterator(fieldValues_ + interiorFieldWindow_.npts(), interiorFieldWindow_); }

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
      builtField_( mode==InternalStorage )
  {
    for( size_t i=0; i<3; ++i ){
      if( window.extent(i)>1 ){
        interiorFieldWindow_.extent(i) -= 2*GhostTraits::NGHOST;
        interiorFieldWindow_.offset(i) +=   GhostTraits::NGHOST;
      }
    }
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
      builtField_( mode==InternalStorage )
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

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  T&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator()( const int i, const int j, const int k )
  {
    return fieldValues_[ fieldWindow_.flat_index(IntVec(i,j,k)) ];
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  T
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator()( const int i, const int j, const int k ) const
  {
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
    iterator iend=this->end();
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
    for( iterator ifld=this->begin(); ifld!=this->end(); ++ifld, ++iother ){
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
    for( iterator ifld=this->begin(); ifld!=this->end(); ++ifld, ++iother ){
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
    for( iterator ifld=this->begin(); ifld!=this->end(); ++ifld, ++iother ){
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
    for( iterator ifld=this->begin(); ifld!=this->end(); ++ifld, ++iother ){
      *ifld /= *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator=(const T a)
  {
    iterator iend=this->end();
    for( iterator ifld=this->begin(); ifld!=iend; ++ifld ) *ifld = a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator+=(const T a)
  {
    for( iterator ifld=this->begin(); ifld!=this->end(); ++ifld ) *ifld += a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator-=(const T a)
  {
    for( iterator ifld=this->begin(); ifld!=this->end(); ++ifld ) *ifld -= a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator*=(const T a)
  {
    for( iterator ifld=this->begin(); ifld!=this->end(); ++ifld ) *ifld *= a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator/=(const T a)
  {
    for( iterator ifld=this->begin(); ifld!=this->end(); ++ifld ) *ifld /= a;
    return *this;
  }

  //------------------------------------------------------------------

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_SpatialField_h
