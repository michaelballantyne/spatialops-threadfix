#ifndef SpatialOps_MemoryWindow_h
#define SpatialOps_MemoryWindow_h

#include <cassert>
#include <iterator>

#include <boost/type_traits.hpp>

namespace SpatialOps{
namespace structured{

  // james: consider re-implementing SpatialField here.
  //  spatial field would take a MemoryWindow object to its constructor?  That way we could define appropriate sizes as well as iterators.
  //  special case for interior_begin would use a MemoryWindow that excluded the ghost cells.  This would enable interior iterators without any problem.
  //  general iterator would use 
  //
  // implement some of the current member functions as free functions:
  //    write_matlab - take begin and end iterators as arguments.
  //    Print

  /**
   *  \class IntVec
   *  \brief provides a lightweight class to deal with a 3D vector of integers.
   */
  class IntVec
  {
    int ijk[3];
  public:
    inline IntVec( const int i, const int j, const int k )
    {
      ijk[0]=i; ijk[1]=j; ijk[2]=k;
    }
    inline IntVec( const int vec[3] )
    {
      ijk[0]=vec[0];  ijk[1]=vec[1];  ijk[2]=vec[2];
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

    IntVec nptsGlob_, offset_, extent_;

  public:

    /**
     *  \brief construct a MemoryWindow object
     *  \param npts the total (global) number of points in each direction
     *  \param offset the offset into the memory
     *  \param extent the size of the block that we are considering
     */
    MemoryWindow( const int npts[3],
                  const int offset[3],
                  const int extent[3] )
      : nptsGlob_( npts ), offset_( offset ), extent_( extent )
    {}
    MemoryWindow( const IntVec& npts, const IntVec& offset, const IntVec& extent )
      : nptsGlob_( npts ), offset_( offset ), extent_( extent )
    {}

    /**
     *  \brief construct a MemoryWindow object
     *  \param npts the total (global) number of points in each direction
     */
    MemoryWindow( const int npts[3] )
      : nptsGlob_( npts ), offset_(0,0,0), extent_( npts )
    {}
    MemoryWindow( const IntVec npts )
      : nptsGlob_( npts ), offset_(0,0,0), extent_( npts )
    {}

    MemoryWindow( const MemoryWindow& other )
      : nptsGlob_( other.nptsGlob_ ), offset_( other.offset_ ), extent_( other.extent_ )
    {}

    ~MemoryWindow(){}

    /**
     *  \brief given the local ijk location (0-based on the local
     *         window), obtain the flat index in the global memory
     *         space.
     */
    inline int flat_index( IntVec loc ) const
    {
      assert( loc[0] < extent_[0] );
      assert( loc[1] < extent_[1] );
      assert( loc[2] < extent_[2] );
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

    inline int npts() const{ return extent_[0] * extent_[1] * extent_[2]; }

    inline int glob_dim( const size_t i ) const{ return nptsGlob_[i]; }
    inline int offset  ( const size_t i ) const{ return offset_[i]; }
    inline int extent  ( const size_t i ) const{ return extent_[i]; }

    inline int& glob_dim( const size_t i ){ return nptsGlob_[i]; }
    inline int& offset  ( const size_t i ){ return offset_[i]; }
    inline int& extent  ( const size_t i ){ return extent_[i]; }

    inline IntVec extent  () const{ return extent_; }
    inline IntVec offset  () const{ return offset_; }
    inline IntVec glob_dim() const{ return nptsGlob_; }

    inline IntVec& extent  (){ return extent_; }
    inline IntVec& offset  (){ return offset_; }
    inline IntVec& glob_dim(){ return nptsGlob_; }


    /**
     *  \brief obtain the stride in the requested direction.
     */
    inline int stride( const size_t i ) const{ return nptsGlob_[i] - extent_[i]; }

  };

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
    const MemoryWindow& window_;
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
        window_( other.window_ )
    {
      i_=other.i_; j_=other.j_; k_=other.k_;
    }
    
    FieldIterator( T* t, const MemoryWindow& window )
      : current_( t ),
        window_( window )
    {
      i_ = j_ = k_ = 0;
    }

    inline self& operator++()
    {
      if( i_<window_.extent(0) ){
        ++i_;
        ++current_;
      }
      else{
        i_=0;
        ++j_;
        current_ += 1+window_.stride(0);
      }
      if( j_>=window_.extent(1) ){
        j_=0;
        ++k_;
        current_ += window_.stride(0) * window_.stride(1);
      }
      assert( k_ < window_.extent(2) );
      return *this;
    }

    inline bool operator==( const self& other ) const{ return current_==other.current_; }

    inline bool operator!=( const self& other ) const{ return current_!=other.current_; }

    inline self& operator=( const self& other )
    {
      current_ = other.current_;
      i_ = other.i_;
      j_ = other.j_;
      k_ = other.k_;
      return *this;
    }

    inline       reference operator*()      { return *current_; }
    inline const reference operator*() const{ return *current_; }
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
    const MemoryWindow& window_;
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
        window_( other.window_ )
    {
      i_=other.i_; j_=other.j_; k_=other.k_;
    }
    
    ConstFieldIterator( const T* t, const MemoryWindow& window )
      : current_( t ),
        window_( window )
    {
      i_ = j_ = k_ = 0;
    }

    ConstFieldIterator( const FieldIterator<T> t )
      : current_( t.current_ ), window_( t.window_ )
    {
      i_ = t.i_;
      j_ = t.j_;
      k_ = t.k_;
    }

    inline self& operator++()
    {
      if( i_<window_.extent(0) ){
        ++i_;
        ++current_;
      }
      else{
        i_=0;
        ++j_;
        current_ += 1+window_.stride(0);
      }
      if( j_>=window_.extent(1) ){
        j_=0;
        ++k_;
        current_ += window_.stride(0) * window_.stride(1);
      }
      assert( k_ < window_.extent(2) );
      return *this;
    }

    inline bool operator==( const self& other ) const{ return current_==other.current_; }

    inline bool operator!=( const self& other ) const{ return current_!=other.current_; }

    inline self& operator=( const self& other )
    {
      current_ = other.current_;
      i_ = other.i_;
      j_ = other.j_;
      k_ = other.k_;
      return *this;
    }

    inline       reference operator*()      { return *current_; }
    inline const reference operator*() const{ return *current_; }

  };

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_MemoryWindow_h
