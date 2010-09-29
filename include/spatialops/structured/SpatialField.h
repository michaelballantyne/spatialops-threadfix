#ifndef SpatialOps_SpatialField_h
#define SpatialOps_SpatialField_h

#include <MemoryWindow.h>

namespace SpatialOps{

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

    typedef SpatialField<VecOps,FieldLocation,GhostTraits> MyType;

    friend class SpatialFieldStore<MyType>;

    const MemoryWindow fieldWindow_;
    MemoryWindow interiorFieldWindow_;

    T* const fieldValues_;

  public:

    typedef GhostTraits Ghost;
    typedef FieldLocation Location;

    typedef FieldIterator<      T*>       iterator;
    typedef FieldIterator<const T*> const_iterator;

    typedef FieldIterator<      T*>       interior_iterator;
    typedef FieldIterator>const T*> const_interior_iterator;

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

    SpatialField( const int[3] npts,
                  T* const fieldValues,
                  const StorageMode mode = InternalStorage );

    virtual ~SpatialField();


    inline SpatialField& operator=(const SpatialField&);
    inline SpatialField& operator+=(const SpatialField&);
    inline SpatialField& operator-=(const SpatialField&);
    inline SpatialField& operator*=(const SpatialField&);
    inline SpatialField& operator/=(const SpatialField&);

    inline const_iterator begin() const{ return const_iterator(fieldWindow_); }
    inline       iterator begin()      { return       iterator(fieldWindow_); }

    inline const_iterator end() const{ return const_iterator( fieldValues_ + fieldWindow_.flat_index( fieldWindow_.extent() ) ); }
    inline       iterator end() const{ return       iterator( fieldValues_ + fieldWindow_.flat_index( fieldWindow_.extent() ) ); }

    inline const_interior_iterator interior_begin() const{ return const_interior_iterator(interiorFieldWindow_;) }
    inline       interior_iterator interior_begin()      { return       interior_iterator(interiorFieldWindow_;) }

    inline const_interior_iterator interior_end() const{ return const_interior_iterator( fieldValues_ + interiorFieldWindow_.flat_index( interiorFieldWindow_.extent() ) ); }
    inline       interior_iterator interior_end()      { return       interior_iterator( fieldValues_ + interiorFieldWindow_.flat_index( interiorFieldWindow_.extent() ) ); }

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
    : fieldWindow_( window )
      interiorFieldWindow_( window ),
      fieldValues_( (mode==ExternalStorage)
                    ? fieldValues
                    : new T[ window.glob_dim(0) * window.glob_dim(1) * window.glob_dim(2) ] )
  {
    for( size_t i=0; i<3; ++i ){
      if( window.extent[i]>1 ){
        interiorFieldWindow_.extent(i) -= GhostTraits::NGHOST;
        interiorFieldWindow_.offset(i) += GhostTraits::NGHOST;
      }
    }
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  SpatialField<VecOps,Location,GhostTraits,T>::
  SpatialField( const int[3] npts,
                T* const fieldValues,
                const StorageMode mode )
    : fieldWindow_( npts, {0,0,0}, npts )
      interiorFieldWindow_( window ),
      fieldValues_( (mode==ExternalStorage)
                    ? fieldValues
                    : new T[npts[0]*npts[1]*npts[2]] )
  {
    for( size_t i=0; i<3; ++i ){
      if( window.extent[i]>1 ){
        interiorFieldWindow_.extent(i) -= GhostTraits::NGHOST;
        interiorFieldWindow_.offset(i) += GhostTraits::NGHOST;
      }
    }
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator=(const MyType& other )
  {
    for( iterator ifld=this->begin(), const_iterator iother=other.begin();
         ifld!=this->end(); ++ifld, ++iother ){
      *ifld = *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator+=(const MyType& other )
  {
    for( iterator ifld=this->begin(), const_iterator iother=other.begin();
         ifld!=this->end(); ++ifld, ++iother ){
      *ifld += *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator-=(const MyType& other )
  {
    for( iterator ifld=this->begin(), const_iterator iother=other.begin();
         ifld!=this->end(); ++ifld, ++iother ){
      *ifld -= *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator*=(const MyType& other )
  {
    for( iterator ifld=this->begin(), const_iterator iother=other.begin();
         ifld!=this->end(); ++ifld, ++iother ){
      *ifld *= *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator/=(const MyType& other )
  {
    for( iterator ifld=this->begin(), const_iterator iother=other.begin();
         ifld!=this->end(); ++ifld, ++iother ){
      *ifld /= *iother;
    }
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator=(const T a)
  {
    for( iterator ifld=this->begin(), ifld!=this->end(); ++ifld ) *ifld = a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator+=(const T a)
  {
    for( iterator ifld=this->begin(), ifld!=this->end(); ++ifld ) *ifld += a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator-=(const T a)
  {
    for( iterator ifld=this->begin(), ifld!=this->end(); ++ifld ) *ifld -= a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator*=(const T a)
  {
    for( iterator ifld=this->begin(), ifld!=this->end(); ++ifld ) *ifld *= a;
    return *this;
  }

  //------------------------------------------------------------------

  template< typename VecOps, typename Location, typename GhostTraits, typename T >
  typename SpatialField<VecOps,Location,GhostTraits,T>&
  SpatialField<VecOps,Location,GhostTraits,T>::
  operator/=(const T a)
  {
    for( iterator ifld=this->begin(), ifld!=this->end(); ++ifld ) *ifld /= a;
    return *this;
  }

  //------------------------------------------------------------------

} // namespace SpatialOps

#endif // SpatialOps_SpatialField_h
