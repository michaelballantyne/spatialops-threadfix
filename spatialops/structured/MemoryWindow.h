#ifndef SpatialOps_MemoryWindow_h
#define SpatialOps_MemoryWindow_h

#include <cassert>
#include <vector>
#include <iterator>

#include <boost/serialization/serialization.hpp>

namespace SpatialOps{
namespace structured{

  /**
   *  \class IntVec
   *  \brief provides a lightweight class to deal with a 3D vector of integers.
   */
  class IntVec
  {
    friend std::ostream& operator<<(std::ostream&, const IntVec&);

    int ijk[3];

    friend class boost::serialization::access;
    template<typename Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
      ar & ijk;
    }

  public:
    IntVec(){ ijk[0]=0; ijk[1]=0; ijk[2]=0; }
    inline IntVec( const int i, const int j, const int k )
    {
      ijk[0]=i; ijk[1]=j; ijk[2]=k;
    }
    inline IntVec( const int vec[3] )
    {
      ijk[0]=vec[0];  ijk[1]=vec[1];  ijk[2]=vec[2];
    }
    IntVec( const std::vector<int>& vec )
    {
      ijk[0]=vec[0]; ijk[1]=vec[1]; ijk[2]=vec[2];
    }
    inline IntVec( const IntVec& x )
    {
      ijk[0]=x.ijk[0];  ijk[1]=x.ijk[1];  ijk[2]=x.ijk[2];
    }

    inline const int  operator[](const size_t i) const{ return ijk[i]; }
    inline       int& operator[](const size_t i)      { return ijk[i]; }

    IntVec& operator=(const IntVec& x)
    {
      for( size_t i=0; i<3; ++i ) ijk[i] = x.ijk[i];
      return *this;
    }

    inline bool operator==(const IntVec& v) const
    {
      return (ijk[0]==v.ijk[0]) & (ijk[1]==v.ijk[1]) & (ijk[2]==v.ijk[2]);
    }
    inline bool operator!=(const IntVec& v) const
    {
      return (ijk[0]!=v.ijk[0]) | (ijk[1]!=v.ijk[1]) | (ijk[2]!=v.ijk[2]);
    }
  };

  /**
   *  \class MemoryWindow
   *  \author James C. Sutherland
   *  \date September 2010
   *
   *  \brief Provides tools to index into a sub-block of memory.
   *
   *  Given a block of memory, [Nx,Ny,Nz], assume that we want to deal
   *  with a sub-block of size [nx,ny,nz] that starts at [i,j,k] in
   *  the larger block.  The MemoryWindow class provides basic tools
   *  to help with this.
   */
  class MemoryWindow{

    friend std::ostream& operator<<( std::ostream&, const MemoryWindow& );

    IntVec nptsGlob_, offset_, extent_;

    friend class boost::serialization::access;

    template<typename Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
      ar & nptsGlob_;
      ar & offset_;
      ar & extent_;
    }

    template<typename Archive>
    void save_construct_data( Archive& ar, const MemoryWindow* w, const unsigned int version )
    {
      ar << w->nptsGlob_ << w->offset_ << w->extent_;
    }

    template<typename Archive>
    void load_construct_data( Archive& ar, const MemoryWindow* w, const unsigned int version )
    {
      IntVec npg, ofs, ext;
      ar >> npg >> ofs >> ext;
      ::new(w)MemoryWindow( npg, ofs, ext );
    }

  public:

    /**
     *  \brief construct a MemoryWindow object
     *  \param npts the total (global) number of points in each direction
     *  \param offset the offset into the memory
     *  \param extent the size of the block that we are considering
     */
    MemoryWindow( const int npts[3],
                  const int offset[3],
                  const int extent[3] );
		
    MemoryWindow( const IntVec& npts,
                  const IntVec& offset,
                  const IntVec& extent );

    /**
     *  \brief construct a MemoryWindow object
     *  \param npts the total (global) number of points in each direction
     */
    MemoryWindow( const int npts[3] );
    MemoryWindow( const IntVec npts );

    MemoryWindow( const MemoryWindow& other );
		
    ~MemoryWindow();

    /**
     *  \brief given the local ijk location (0-based on the local
     *         window), obtain the flat index in the global memory
     *         space.
     */
    inline int flat_index( IntVec loc ) const
    {
#     ifndef NDEBUG
      if( extent_[0]>1 ) assert( loc[0] < extent_[0] );
      if( extent_[1]>1 ) assert( loc[1] < extent_[1] );
      if( extent_[2]>1 ) assert( loc[2] < extent_[2] );
#     endif
      for( size_t i=0; i<3; ++i ) loc[i] += offset_[i];
      return loc[0] + loc[1]*nptsGlob_[0] + loc[2]*nptsGlob_[0]*nptsGlob_[1];
    }

    /**
     *  \brief given the local flat location (0-based on the local
     *         window), obtain the ijk index in the global memory
     *         space.
     */
    inline IntVec ijk_index( const int loc ) const
    {
      IntVec ijk( 0,0,0 );
      ijk[0] = loc % extent_[0]              + offset_[0];
      ijk[1] = loc / extent_[0] % extent_[1] + offset_[1];
      ijk[2] = loc / (extent_[0]*extent_[1]) + offset_[2];
      return ijk;
    }

    /**
     *  \brief obtain the number of points in the field.  Note that
     *  this is not necessarily contiguous memory
     */
    inline size_t npts() const{ return extent_[0] * extent_[1] * extent_[2]; }

    inline size_t glob_dim( const size_t i ) const{ return size_t(nptsGlob_[i]); }
    inline size_t offset  ( const size_t i ) const{ return size_t(offset_[i]); }
    inline size_t extent  ( const size_t i ) const{ return size_t(extent_[i]); }

    inline IntVec extent  () const{ return extent_; }
    inline IntVec offset  () const{ return offset_; }
    inline IntVec glob_dim() const{ return nptsGlob_; }

    inline IntVec& extent  (){ return extent_; }
    inline IntVec& offset  (){ return offset_; }
    inline IntVec& glob_dim(){ return nptsGlob_; }

    /**
     *  \brief obtain the stride in the requested direction.
     */
    inline int stride( const size_t i ) const{
      const int n = 1 + nptsGlob_[i] - extent_[i];
      return n;
    }

    inline bool operator==( const MemoryWindow& w ) const
    {
      return ( (nptsGlob_==w.nptsGlob_) & (extent_==w.extent_) & (offset_==w.offset_) );
    }

    inline bool operator!=( const MemoryWindow& w ) const
    {
      return (nptsGlob_!=w.nptsGlob_) | (extent_!=w.extent_) | (offset_!=w.offset_);
    }

  };

  template<int Dir> size_t stride( const MemoryWindow& mw );

  template<> inline size_t stride<0>( const MemoryWindow& mw ){ return 1; }
  template<> inline size_t stride<1>( const MemoryWindow& mw ){ return stride<0>(mw) + mw.glob_dim(0)-mw.extent(0); }
  template<> inline size_t stride<2>( const MemoryWindow& mw )
  {
    return stride<0>(mw)
      + ( mw.glob_dim(0)-mw.extent(0) ) * mw.glob_dim(1)
      + ( mw.glob_dim(1)-mw.extent(1) );
  }


  template<typename T> class ConstFieldIterator; // forward

  /**
   *  \class FieldIterator
   *  \author James C. Sutherland
   *  \date September, 2010
   *
   *  \brief Provides a forward iterator for a field that is
   *         associated with a MemoryWindow, allowing one to iterate
   *         over the "local" portion of that field as defined by the
   *         MemoryWindow.
   */
  template<typename T>
  class FieldIterator
  {
    friend class ConstFieldIterator<T>;
    T* current_;
    T* first_;
    const MemoryWindow& window_;
    IntVec stride_;
    size_t i_,j_,k_;

  public:
    typedef FieldIterator<T> self;
    typedef typename std::iterator_traits<T*>::value_type      value_type;
    typedef typename std::iterator_traits<T*>::reference       reference;
    typedef typename std::iterator_traits<T*>::pointer         pointer;
    typedef typename std::iterator_traits<T*>::difference_type difference_type;
    typedef          std::forward_iterator_tag                iterator_category;

    FieldIterator( const self& other )
      : current_( other.current_ ),
        first_  ( other.first_   ),
        window_ ( other.window_  )
    {
      stride_ = other.stride_;
      i_=other.i_; j_=other.j_; k_=other.k_;
    }
    
    FieldIterator( T* t, const size_t offset, const MemoryWindow& window )
      : current_( t+offset ),
        first_  ( t        ),
        window_ ( window   )
    {
      stride_[0] = stride<0>(window);
      stride_[1] = stride<1>(window);
      stride_[2] = stride<2>(window);
      i_ = j_ = k_ = 0;
    }

    inline self& operator++()
    {
      ++i_;
      if( i_<window_.extent(0) ){
        current_ += stride_[0];
      }
      else{
        i_=0;
        ++j_;
        if( j_ < window_.extent(1) )  current_ += stride_[1];
        else{
          j_=0;
          ++k_;
          if( k_ < window_.extent(2) ) current_ += stride_[2];
          else{
            IntVec ijkend = window_.extent();
            --ijkend[0]; --ijkend[1]; --ijkend[2];
            const size_t last = window_.flat_index( ijkend ) + 1;
            current_ = first_ + last;
          }
        }
      }
      return *this;
    }

    inline self& operator+( const size_t n )
    {
      for( size_t i=0; i<n; ++i )  ++(*this);
      return *this;
    }

    inline self& operator+=( const size_t n )
    {
      for( size_t i=0; i<n; ++i )  ++(*this);
      return *this;
    }

    inline bool operator==( const self& other ) const{ return current_==other.current_; }

    inline bool operator!=( const self& other ) const{ return current_!=other.current_; }

    inline self& operator=( const self& other )
    {
      current_ = other.current_;
      first_   = other.first_;
      i_       = other.i_;
      j_       = other.j_;
      k_       = other.k_;
      return *this;
    }

    inline reference operator*()
    {
#     ifndef NDEBUG
      if( window_.extent(2) > 1 )  assert( k_ < window_.extent(2) );
      if( window_.extent(1) > 1 )  assert( j_ < window_.extent(1) );
      if( window_.extent(0) > 1 )  assert( i_ < window_.extent(0) );
#     endif
      return *current_;
    }

    inline const reference operator*() const
    {
#     ifndef NDEBUG
      if( window_.extent(2) > 1 )  assert( k_ < window_.extent(2) );
      if( window_.extent(1) > 1 )  assert( j_ < window_.extent(1) );
      if( window_.extent(0) > 1 )  assert( i_ < window_.extent(0) );
#     endif
      return *current_;
    }

    inline size_t i(){ return i_; }
    inline size_t j(){ return j_; }
    inline size_t k(){ return k_; }
  };


  //==================================================================


  /**
   *  \class ConstFieldIterator
   *  \author James C. Sutherland
   *  \date September, 2010
   *
   *  \brief Provides a forward iterator for a field that is
   *         associated with a MemoryWindow, allowing one to iterate
   *         over the "local" portion of that field as defined by the
   *         MemoryWindow.
   */
  template<typename T>
  class ConstFieldIterator
  {
    const T* current_;
    const T* first_;
    const MemoryWindow& window_;
    IntVec stride_;
    size_t i_,j_,k_;

  public:
    typedef ConstFieldIterator<T> self;
    typedef typename std::iterator_traits<const T*>::value_type      value_type;
    typedef typename std::iterator_traits<const T*>::reference       reference;
    typedef typename std::iterator_traits<const T*>::pointer         pointer;
    typedef typename std::iterator_traits<const T*>::difference_type difference_type;
    typedef          std::forward_iterator_tag                       iterator_category;

    ConstFieldIterator( const self& other )
      : current_( other.current_ ),
        first_  ( other.first_   ),
        window_ ( other.window_  )
    {
      stride_ = other.stride_;
      i_=other.i_; j_=other.j_; k_=other.k_;
    }
    
    ConstFieldIterator( const T* t, const size_t offset, const MemoryWindow& window )
      : current_( t+offset ),
        first_  ( t        ),
        window_ ( window   )
    {
      stride_[0] = stride<0>(window);
      stride_[1] = stride<1>(window);
      stride_[2] = stride<2>(window);
      i_ = j_ = k_ = 0;
    }

    ConstFieldIterator( const FieldIterator<T> t )
      : current_( t.current_ ),
        first_  ( t.first_   ),
        window_ ( t.window_  )
    {
      stride_ = t.stride_;
      i_ = t.i_;
      j_ = t.j_;
      k_ = t.k_;
    }

    inline self& operator++()
    {
      ++i_;
      if( i_<window_.extent(0) )  current_ += stride_[0];
      else{
        i_=0;
        ++j_;
        if( j_ < window_.extent(1) )  current_ += stride_[1];
        else{
          j_=0;
          ++k_;
          if( k_ < window_.extent(2) )  current_ += stride_[2];
          else{
            IntVec ijkend = window_.extent();
            --ijkend[0]; --ijkend[1]; --ijkend[2];
            const size_t last = window_.flat_index( ijkend ) + 1;
            current_ = first_ + last;
          }
        }
      }
      return *this;
    }

    inline self& operator+( const size_t n )
    {
      for( size_t i=0; i<n; ++i )  ++(*this);
      return *this;
    }

    inline self& operator+=( const size_t n )
    {
      for( size_t i=0; i<n; ++i )  ++(*this);
      return *this;
    }

    inline bool operator==( const self& other ) const{ return current_==other.current_; }

    inline bool operator!=( const self& other ) const{ return current_!=other.current_; }

    inline self& operator=( const self& other )
    {
      current_ = other.current_;
      first_   = other.first_;
      i_       = other.i_;
      j_       = other.j_;
      k_       = other.k_;
      return *this;
    }

    inline const reference operator*() const
    {
#     ifndef NDEBUG
      if( window_.extent(2) > 1 )  assert( k_ < window_.extent(2) );
      if( window_.extent(1) > 1 )  assert( j_ < window_.extent(1) );
      if( window_.extent(0) > 1 )  assert( i_ < window_.extent(0) );
#     endif
      return *current_;
    }

  };

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_MemoryWindow_h
