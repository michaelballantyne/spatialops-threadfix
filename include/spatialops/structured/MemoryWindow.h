#ifndef SpatialOps_MemoryWindow_h
#define SpatialOps_MemoryWindow_h

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

    const int[3] nptsGlob_;
    const int[3] offset_;
    const int[3] extent_;

  public:

    /**
     *  \brief construct a MemoryWindow object
     *  \param npts the total (global) number of points in each direction
     *  \param offset the offset into the memory
     *  \param extent the size of the block that we are considering
     */
    MemoryWindow( const int[3] npts,
                  const int[3] offset,
                  const int[3] extent )
      : nptsGlob_( npts ), offset_( offset ), extent_( extent )
    {}

    /**
     *  \brief construct a MemoryWindow object
     *  \param npts the total (global) number of points in each direction
     */
    MemoryWindow( const int[3] npts )
      : nptsGlob_( npts ), offset_( npts ), extent_( npts )
    {}

    MemoryWindow( const MemoryWindow& other )
      : nptsGlob_( other.nptsGlob_ ),
        offset_( other.offset_ ),
        extent_( other.extent_ )
    {}

    ~MemoryWindow(){}

    /**
     *  \brief given the local ijk location (0-based on the local
     *         window), obtain the flat index in the global memory
     *         space.
     */
    inline int flat_index( int[3] loc ) const
    {
      for( size_t i=0; i<3; ++i ) loc[i]+=offset[i];
      return loc[0] + loc[1]*npts[1] + loc[2]*npts[1]*npts[2];
    }

    /**
     *  \brief given the local flat location (0-based on the local
     *         window), obtain the ijk index in the global memory
     *         space.
     */
    inline int[3] ijk_index( const int loc ) const
    {
      int[3] ijk={0 0 0};
      ijk[0] = loc % extent_[0]              + offset_[0];
      ijk[1] = loc / extent_[0] % extent_[1] + offset_[1];
      ijk[2] = loc / (extent_[0]*extent_[1]) + offset_[2];
      return ijk;
    }

    inline int glob_dim( const size_t i ) const{ return nptsGlob_[i]; }
    inline int offset  ( const size_t i ) const{ return offset_[i]; }
    inline int extent  ( const size_t i ) const{ return extent_[i]; }

    inline int& glob_dim( const size_t i ){ return nptsGlob_[i]; }
    inline int& offset  ( const size_t i ){ return offset_[i]; }
    inline int& extent  ( const size_t i ){ return extent_[i]; }

    inline int[3] extent  () const{ return extent_; }
    inline int[3] offset  () const{ return offset_; }
    inline int[3] glob_dim() const{ return nptsGlob_; }

    inline int[3]& extent  (){ return extent_; }
    inline int[3]& offset  (){ return offset_; }
    inline int[3]& glob_dim(){ return nptsGlob_; }


    /**
     *  \brief obtain the stride in the requested direction.
     */
    inline int stride( const size_t i ) const{ return nptsGlob_[i] - extent_[i]; }
  };


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
    T current_;
    const MemoryWindow window_;
    size_t i_,j_,k_;

  public:
    typedef FieldIterator<T> self;
    typedef typename std::iterator_traits<T>::value_type   value_type;
    typedef typename std::iterator_traits<T>::reference       reference;
    typedef typename std::iterator_traits<T>::pointer         pointer;
    typedef typename std::iterator_traits<T>::difference_type difference_type;
    typedef          std::forward_iterator_tag                              iterator_category;

    FieldIterator( T first, const MemoryWindow window )
      : current_( first ),
        window_( window ),
    {
      i_ = j_ = k_ = 0;

      // shift current to the beginning
      current_ += window_.flat_index( {0,0,0} );
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
        current_ += window_.stride(0);
      }
      if( j_>=window_.extent(1) ){
        j_=0;
        ++k_;
        current_ += window_.stride(0) * window_.stride(1);
      }
      assert( k_<window_.extent(2) );
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
  };

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_MemoryWindow_h
