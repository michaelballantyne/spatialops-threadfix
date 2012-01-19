#ifndef SpatialOps_MemoryWindow_h
#define SpatialOps_MemoryWindow_h

#include <vector>
#include <iterator>

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>

#include <spatialops/structured/IntVec.h>

#ifdef SOPS_BOOST_SERIALIZATION
# include <boost/serialization/serialization.hpp>
#endif

#ifndef NDEBUG
# include <cassert>
# include <sstream>
# include <stdexcept>
#endif

/**
 * \file MemoryWindow.h
 * \addtogroup structured
 * @{
 *
 */

namespace SpatialOps{
namespace structured{

  /**
   *  \class MemoryWindow
   *  \author James C. Sutherland
   *  \date September 2010
   *
   *  \ingroup structured
   *  \brief Provides tools to index into a sub-block of memory.
   *
   *  Given a block of memory, [Nx,Ny,Nz], assume that we want to deal
   *  with a sub-block of size [nx,ny,nz] that starts at [i,j,k] in
   *  the larger block.  The MemoryWindow class provides basic tools
   *  to help with this.
   */
  class MemoryWindow{

    friend std::ostream& operator<<( std::ostream&, const MemoryWindow& );

    IntVec nptsGlob_;   ///< The global number of points
    IntVec offset_;     ///< The offset for this window
    IntVec extent_;     ///< The extent of this window
    IntVec bc_;         ///< Indicates presence of a physical boundary on the + side of the window

#   ifdef SOPS_BOOST_SERIALIZATION
    friend class boost::serialization::access;

    template<typename Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
      ar & nptsGlob_;
      ar & offset_;
      ar & extent_;
      ar & bc_;
    }

    template<typename Archive>
    void save_construct_data( Archive& ar, const MemoryWindow* w, const unsigned int version )
    {
      ar << w->nptsGlob_ << w->offset_ << w->extent_ << w->bc_;
    }

    template<typename Archive>
    void load_construct_data( Archive& ar, const MemoryWindow* w, const unsigned int version )
    {
      IntVec npg, ofs, ext, bc;
      ar >> npg >> ofs >> ext >> bc;
      ::new(w)MemoryWindow( npg, ofs, ext, bc[0], bc[1], bc[2] );
    }
#   endif

  public:

    /**
     *  \brief construct a MemoryWindow object
     *  \param npts the total (global) number of points in each direction
     *  \param offset the offset into the memory
     *  \param extent the size of the block that we are considering
     *  \param bcx true if a physical boundary is present in the (+x) direction
     *  \param bcy true if a physical boundary is present in the (+y) direction
     *  \param bcz true if a physical boundary is present in the (+z) direction
     */
    MemoryWindow( const int npts[3],
                  const int offset[3],
                  const int extent[3],
                  const bool bcx,
                  const bool bcy,
                  const bool bcz );

    /**
     *  \brief construct a MemoryWindow object
     *
     *  \param npts the total (global) number of points in each direction
     *  \param offset the offset into the memory
     *  \param extent the size of the block that we are considering
     *  \param bcx true if a physical boundary is present in the (+x) direction
     *  \param bcy true if a physical boundary is present in the (+y) direction
     *  \param bcz true if a physical boundary is present in the (+z) direction
     */
    MemoryWindow( const IntVec& npts,
                  const IntVec& offset,
                  const IntVec& extent,
                  const bool bcx,
                  const bool bcy,
                  const bool bcz );

    /**
     *  \brief construct a MemoryWindow object where there is no "window"
     *
     *  \param npts the total (global) number of points in each direction
     *  \param bcx (optional - default false) true if a physical boundary is present in the (+x) direction
     *  \param bcy (optional - default false) true if a physical boundary is present in the (+y) direction
     *  \param bcz (optional - default false) true if a physical boundary is present in the (+z) direction
     */
    MemoryWindow( const int npts[3],
                  const bool bcx=false,
                  const bool bcy=false,
                  const bool bcz=false );

    /**
     *  \brief construct a MemoryWindow object where there is no "window"
     *
     *  \param npts the total (global) number of points in each direction
     *  \param bcx (optional - default false) true if a physical boundary is present in the (+x) direction
     *  \param bcy (optional - default false) true if a physical boundary is present in the (+y) direction
     *  \param bcz (optional - default false) true if a physical boundary is present in the (+z) direction
     */
    MemoryWindow( const IntVec& npts,
                  const bool bcx=false,
                  const bool bcy=false,
                  const bool bcz=false );

    MemoryWindow( const MemoryWindow& other );

    ~MemoryWindow();

    /**
     *  \brief Splits the MemoryWindow into a series of child windows ordered
     *         as a vector varying in x then y then z.
     *
     *  \param splitPattern the number of partitions to make in each ordinate direction.
     *
     *  \param npad (default [0,0,0]) Specifies the number of cells in
     *    each direction that we want to pad the children window with.
     *    Note that you should ensure that the parent window is also
     *    padded with the same number of cells to avoid problems.
     *
     *  \return vector<MemoryWindow> containing the child windows.
     */
    std::vector<MemoryWindow> split( const IntVec splitPattern,
                                     const IntVec npad=IntVec(0,0,0),
                                     const IntVec bcExtents=IntVec(0,0,0) ) const;

    /**
     *  \brief given the local ijk location (0-based on the local
     *         window), obtain the flat index in the global memory
     *         space.
     */
    inline int flat_index( IntVec loc ) const
    {
#     ifndef NDEBUG
      if( extent_[0]>1 ) assert( loc[0] < nptsGlob_[0] );
      if( extent_[1]>1 ) assert( loc[1] < nptsGlob_[1] );
      if( extent_[2]>1 ) assert( loc[2] < nptsGlob_[2] );
#     endif
      loc[0] = nptsGlob_[0] > 1 ? loc[0]+offset_[0] : 0;
      loc[1] = nptsGlob_[1] > 1 ? loc[1]+offset_[1] : 0;
      loc[2] = nptsGlob_[2] > 1 ? loc[2]+offset_[2] : 0;
      return loc[0] + nptsGlob_[0] * (loc[1] + loc[2]*nptsGlob_[1]);
    }

    /**
     *  \brief given the local flat location (0-based on the global
     *         field), obtain the ijk index in the global memory
     *         space.
     */
    inline IntVec ijk_index_from_global( const int loc ) const
    {
      return IntVec( loc % nptsGlob_[0],
                     loc / nptsGlob_[0] % nptsGlob_[1],
                     loc /(nptsGlob_[0] * nptsGlob_[1]) );
    }

    /**
     *  \brief given the local flat location (0-based on the local
     *         window), obtain the ijk index in the global memory
     *         space.
     */
    inline IntVec ijk_index_from_local( const int loc ) const
    {
      return IntVec( loc % extent_[0]               + offset_[0],
                     loc / extent_[0] % extent_[1]  + offset_[1],
                     loc /(extent_[0] * extent_[1]) + offset_[2] );
    }

    /**
     *  \brief obtain the global number of points in the field.  Note
     *  that this is not necessarily contiguous memory
     */
    inline size_t glob_npts() const{ return nptsGlob_[0] * nptsGlob_[1] * nptsGlob_[2]; }

    /**
     *  \brief obtain the local number of points in the field.  Note
     *  that this is not necessarily contiguous memory
     */
    inline size_t local_npts() const{ return extent_[0] * extent_[1] * extent_[2]; }

    inline size_t glob_dim( const size_t i ) const{ return size_t(nptsGlob_[i]); }
    inline size_t offset  ( const size_t i ) const{ return size_t(offset_[i]  ); }
    inline size_t extent  ( const size_t i ) const{ return size_t(extent_[i]  ); }

    inline const IntVec& extent  () const{ return extent_;   }
    inline const IntVec& offset  () const{ return offset_;   }
    inline const IntVec& glob_dim() const{ return nptsGlob_; }

    /**
     * \brief Query if there is a physical BC on the + side of this window in the given direction.
     * @param i The direction: (x,y,z) = (0,1,2)
     * @return true if a BC is present on the + side of this window, false otherwise.
     */
    inline bool has_bc( const size_t i ) const{ return bc_[i]; }

    /**
     * \brief Obtain an IntVec indicating the presence of a physical BC
     *        (1=present, 0=not present) on each + side of this window.
     * @return IntVec indicating the presence of a BC on the + side of this window.
     */
    inline const IntVec& has_bc() const{ return bc_; }

    /**
     * \brief compare two MemoryWindows for equality
     */
    inline bool operator==( const MemoryWindow& w ) const{
      return (nptsGlob_ == w.nptsGlob_) &&
             (extent_   == w.extent_  ) &&
             (offset_   == w.offset_  ) &&
             (bc_       == w.bc_      );
    }

    inline bool operator!=( const MemoryWindow& w ) const{
      return (nptsGlob_ != w.nptsGlob_) ||
             (extent_   != w.extent_  ) ||
             (offset_   != w.offset_  ) ||
             (bc_       != w.bc_      );
    }

  };

  /**
   *  \class IteratorIncrementor
   *  \brief provides increment/decrement facilities for iterators
   *
   *  \tparam DirT - the direction to increment the iterator in.
   */
  template< typename DirT > class IteratorIncrementor;

  template<>
  struct IteratorIncrementor<XDIR>
  {
  private:

    static inline size_t stride_x( const MemoryWindow* const mw ){ return 1; }
    static inline size_t stride_y( const MemoryWindow* const mw ){ return stride_x(mw) + mw->glob_dim(0)-mw->extent(0); }
    static inline size_t stride_z( const MemoryWindow* const mw ){ return stride_y(mw) + ( mw->glob_dim(0) ) * ( mw->glob_dim(1)-mw->extent(1) ); }

  public:

    /**
     *  \brief returns the increment for the pointer
     */
    static size_t increment( const MemoryWindow* const w,
                             size_t& i, size_t& j, size_t& k )
    {
      size_t inc=0;
      ++i;
      if( i < w->extent(0) ){
        inc += stride_x(w);
      }
      else{
        i=0;
        ++j;
        if( j < w->extent(1) ){
          inc += stride_y(w);
        }
        else{
          j=0;
          ++k;
          inc += stride_z(w);
        }
      }
      return inc;
    }

    /**
     *  \brief returns the decrement for the pointer.  If zero, then reset the pointer.
     */
    static size_t decrement( const MemoryWindow* const w,
                             size_t& i, size_t& j, size_t& k )
    {
      size_t inc=0;
      if( i > 0 ){
        --i;
        ++inc;
      }
      else{
        i = w->extent(0)-1;
        if( j > 0 ){
          --j;
          inc += stride_y(w);
        }
        else{
          j = w->extent(1)-1;
          if( k > 0 ){
            --k;
            inc += stride_z(w);
          }
          else{
            i = j = k = 0;
          }
        }
      }
      return inc;
    }

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
  template< typename FieldT >
  class FieldIterator
  {
    typedef typename FieldT::AtomicT  T;
    friend class ConstFieldIterator<FieldT>;
    T* current_;   ///< The current pointer that this iterator refers to
    T* first_;     ///< The first position in memory for the field this iterator is associated with
    const MemoryWindow* window_;  ///< The MemoryWindow associated with this field and iterator
    size_t i_,j_,k_;

  public:
    typedef FieldIterator<FieldT> self;
    typedef typename std::iterator_traits<T*>::value_type      value_type;
    typedef typename std::iterator_traits<T*>::reference       reference;
    typedef typename std::iterator_traits<T*>::pointer         pointer;
    typedef typename std::iterator_traits<T*>::difference_type difference_type;
    typedef          std::forward_iterator_tag                iterator_category;

    FieldIterator()
      : current_( NULL ), first_( NULL ), window_( NULL )
    {
      i_ = j_ = k_ = 0;
    }

    /**
     *  \brief Copy constructor
     */
    FieldIterator( const self& other )
      : current_( other.current_ ),
        first_  ( other.first_   ),
        window_ ( other.window_  )
    {
      i_=other.i_; j_=other.j_; k_=other.k_;
    }

    /**
     *  \brief Construct a FieldIterator
     *  \param t the raw pointer to the begin location of the field
     *  \param offset the global offset into the field where we want the iterator to be.
     *  \param window the MemoryWindow describing the portion of this field that we have access to.
     */
    FieldIterator( T* t, const size_t offset, const MemoryWindow* const window )
      : current_( t+offset ),
        first_  ( t        ),
        window_ ( window   )
    {
      const IntVec& ijk = window->ijk_index_from_global( offset );
      i_ = ijk[0] - window->offset(0);
      j_ = ijk[1] - window->offset(1);
      k_ = ijk[2] - window->offset(2);
    }

    /**
     *  \brief increment the iterator
     */
    inline self& operator++(){
      current_ += IteratorIncrementor<XDIR>::increment( window_, i_, j_, k_ );
      return *this;
    }

    /**
     *  \brief decrement operator
     */
    inline self& operator--(){
      const size_t dec = IteratorIncrementor<XDIR>::decrement( window_, i_, j_, k_ );
      if( dec==0 ) current_ = first_;
      else current_ -= dec;
      return *this;
    }

    inline self operator+( const size_t n ) const{
      self iter(*this);
      iter+=n;
      return iter;
    }

    inline self operator-( const size_t n ) const{
      self iter(*this);
      iter -= n;
      return iter;
    }

    inline self& operator+=( const size_t n ){
      for( size_t i=0; i<n; ++i )  ++(*this);
      return *this;
    }

    inline self& operator-=( const size_t n ){
      for( size_t i=0; i<n; ++i )  --(*this);
      return *this;
    }

    inline bool operator==( const self& other ) const{ return current_==other.current_; }

    inline bool operator!=( const self& other ) const{ return current_!=other.current_; }

    inline self& operator=( const self& other ){
      current_ = other.current_;
      first_   = other.first_;
      i_       = other.i_;
      j_       = other.j_;
      k_       = other.k_;
      return *this;
    }

    inline reference operator*(){
#     ifndef NDEBUG
      if( i_ >= window_->extent(0) + window_->offset(0) ||
          j_ >= window_->extent(1) + window_->offset(1) ||
          k_ >= window_->extent(2) + window_->offset(2) ){
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
            << "iterator is in an invalid state for dereference";
        throw std::runtime_error( msg.str() );
      }
#     endif
      return *current_;
    }

    inline reference operator*() const{
#     ifndef NDEBUG
      if( i_ >= window_->extent(0) + window_->offset(0) ||
          j_ >= window_->extent(1) + window_->offset(1) ||
          k_ >= window_->extent(2) + window_->offset(2) ){
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
            << "iterator is in an invalid state for dereference";
        throw std::runtime_error( msg.str() );
      }
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
   *
   *  See the documentation for FieldIterator
   */
  template< typename FieldT >
  class ConstFieldIterator
  {
    typedef typename FieldT::AtomicT  T;
    const T* current_;
    const T* first_;
    const MemoryWindow* window_;
    size_t i_, j_, k_;

  public:
    typedef ConstFieldIterator<FieldT> self;
    typedef typename std::iterator_traits<const T*>::value_type      value_type;
    typedef typename std::iterator_traits<const T*>::reference       reference;
    typedef typename std::iterator_traits<const T*>::pointer         pointer;
    typedef typename std::iterator_traits<const T*>::difference_type difference_type;
    typedef          std::forward_iterator_tag                       iterator_category;

    ConstFieldIterator()
      : current_( NULL ), first_( NULL ), window_( NULL )
    {
      i_ = j_ = k_ = 0;
    }

    ConstFieldIterator( const self& other )
      : current_( other.current_ ),
        first_  ( other.first_   ),
        window_ ( other.window_  )
    {
      i_=other.i_; j_=other.j_; k_=other.k_;
    }

    ConstFieldIterator( const T* t, const size_t offset, const MemoryWindow* const window )
      : current_( t+offset ),
        first_  ( t        ),
        window_ ( window   )
    {
      const IntVec ijk = window->ijk_index_from_global( offset );
      i_ = ijk[0] - window->offset(0);
      j_ = ijk[1] - window->offset(1);
      k_ = ijk[2] - window->offset(2);
    }

    /**
     *  \brief Copy constructor to promote a FieldIterator to a ConstFieldIterator
     */
    ConstFieldIterator( const FieldIterator<FieldT> t )
      : current_( t.current_ ),
        first_  ( t.first_   ),
        window_ ( t.window_  )
    {
      i_ = t.i_;
      j_ = t.j_;
      k_ = t.k_;
    }

    inline self& operator++(){
      current_ += IteratorIncrementor<XDIR>::increment( window_, i_, j_, k_ );
      return *this;
    }

    inline self& operator--(){
      const size_t dec = IteratorIncrementor<XDIR>::decrement( window_, i_, j_, k_ );
      if( dec==0 ) current_ = first_;
      else current_ -= dec;
      return *this;
    }

    inline self operator+( const size_t n ) const{
      self iter(*this);
      iter += n;
      return iter;
    }

    inline self operator-( const size_t n ) const{
      self iter(*this);
      iter -= n;
      return iter;
    }

    inline self& operator+=( const size_t n ){
      for( size_t i=0; i<n; ++i )  ++(*this);
      return *this;
    }

    inline self& operator-=( const size_t n ){
      for( size_t i=0; i<n; ++i )  --(*this);
      return *this;
    }

    inline bool operator==( const self& other ) const{ return current_==other.current_; }

    inline bool operator!=( const self& other ) const{ return current_!=other.current_; }

    inline self& operator=( const self& other ){
      current_ = other.current_;
      first_   = other.first_;
      i_       = other.i_;
      j_       = other.j_;
      k_       = other.k_;
      return *this;
    }

    inline reference operator*() const{
#     ifndef NDEBUG
      if( i_ >= window_->extent(0) + window_->offset(0) ||
          j_ >= window_->extent(1) + window_->offset(1) ||
          k_ >= window_->extent(2) + window_->offset(2) ){
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
            << "iterator is in an invalid state for dereference";
        throw std::runtime_error( msg.str() );
      }
#     endif
      return *current_;
    }

  };

} // namespace structured
} // namespace SpatialOps

/**
 * @}
 */

#endif // SpatialOps_MemoryWindow_h
