/*
 * Copyright (c) 2014 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef SpatialOps_MemoryWindow_h
#define SpatialOps_MemoryWindow_h

#include <cstddef>
#include <vector>
#include <iterator>
# include <sstream>
# include <iostream>

#include <spatialops/SpatialOpsConfigure.h>

#include <spatialops/structured/IntVec.h>
#include <spatialops/structured/GhostData.h>
#include <spatialops/structured/IndexTriplet.h> // jcs remove...

#ifndef NDEBUG
# include <cassert>
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
                  const int extent[3] );

    /**
     *  \brief construct a MemoryWindow object
     *
     *  \param npts the total (global) number of points in each direction
     *  \param offset the offset into the memory
     *  \param extent the size of the block that we are considering
     */
    MemoryWindow( const IntVec& npts,
                  const IntVec& offset,
                  const IntVec& extent );

    /**
     *  \brief construct a MemoryWindow object where there is no "window"
     *
     *  \param npts the total (global) number of points in each direction
     */
    MemoryWindow( const int npts[3] );

    /**
     *  \brief construct a MemoryWindow object where there is no "window"
     *
     *  \param npts the total (global) number of points in each direction
     *  \param bcx (optional - default false) true if a physical boundary is present in the (+x) direction
     *  \param bcy (optional - default false) true if a physical boundary is present in the (+y) direction
     *  \param bcz (optional - default false) true if a physical boundary is present in the (+z) direction
     */
    MemoryWindow( const IntVec& npts );

    // inline to improve performance since this is frequently called from nebo
    MemoryWindow( const MemoryWindow& other )
    : nptsGlob_( other.nptsGlob_ ),
      offset_  ( other.offset_   ),
      extent_  ( other.extent_   )
    {}

    MemoryWindow& operator=( const MemoryWindow& other );

    ~MemoryWindow();

    /**
     *  \brief given the local ijk location (0-based on the local
     *         window), obtain the flat index in the global memory
     *         space.
     */
    inline int flat_index( IntVec loc ) const{
#     ifndef NDEBUG
      if( extent_[0]>1 ) assert( loc[0] < nptsGlob_[0] );
      if( extent_[1]>1 ) assert( loc[1] < nptsGlob_[1] );
      if( extent_[2]>1 ) assert( loc[2] < nptsGlob_[2] );
#     endif
      loc[0] = nptsGlob_[0] > 1 ? loc[0]+offset_[0] : 0;
      loc[1] = nptsGlob_[1] > 1 ? loc[1]+offset_[1] : 0;
      loc[2] = nptsGlob_[2] > 1 ? loc[2]+offset_[2] : 0;
      return ijk_to_flat(nptsGlob_,loc);
    }

    /**
     *  \brief given the local flat location (0-based on the global
     *         field), obtain the ijk index in the global memory
     *         space.
     */
    inline IntVec ijk_index_from_global( const int loc ) const{
      return flat_to_ijk(nptsGlob_,loc);
    }

    /**
     *  \brief given the local flat location (0-based on the local
     *         window), obtain the ijk index in the global memory
     *         space.
     */
    inline IntVec ijk_index_from_local( const int loc ) const{
      return flat_to_ijk(extent_,loc);
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
     * \brief compare two MemoryWindows for equality
     */
    inline bool operator==( const MemoryWindow& w ) const{
      return (nptsGlob_ == w.nptsGlob_) &&
             (extent_   == w.extent_  ) &&
             (offset_   == w.offset_  );
    }

    inline bool operator!=( const MemoryWindow& w ) const{
      return (nptsGlob_ != w.nptsGlob_) ||
             (extent_   != w.extent_  ) ||
             (offset_   != w.offset_  );
    }

    /**
     * \brief return the current memory window with oldGhosts removed and newGhosts added
     *
     * \param oldGhosts number of ghost cells to remove from window
     * \param newGhosts number of ghost cells to add to window
     */
    MemoryWindow reset_ghosts( GhostData const & oldGhosts,
                               GhostData const & newGhosts ) const {
      IntVec oldTotal = oldGhosts.get_minus() + oldGhosts.get_plus();
      IntVec newTotal = newGhosts.get_minus() + newGhosts.get_plus();

      return MemoryWindow( glob_dim(),
                           offset() + oldGhosts.get_minus() - newGhosts.get_minus(),
                           extent() - oldTotal + newTotal );
    }

    /**
     * \brief return the current memory window with given ghosts removed
     *
     * \param ghosts number of ghost cells to remove from window
     */
    MemoryWindow remove_ghosts( GhostData const & ghosts ) const {
      return (*this).reset_ghosts( ghosts, GhostData(0,0,0,0,0,0) );
    }

    /**
     * \brief return if the current window fits in given window
     *
     * \param outer window into which current window fits
     *
     * If current and outer are equal, returns true.
     *
     * If global dimensions are different, returns false.
     */
    bool fits_in( MemoryWindow const & outer ) const {
      return ( glob_dim() == outer.glob_dim() &&
               offset()   >= outer.offset()   &&
               extent()   <= outer.extent()   );
    }

    /**
     * \brief return if current window fits between the
     *  two given windows
     *
     * \param inner window to fit in checked
     * \param outer window into which checked fits
     *
     * If current, inner, and outer are equal, returns true.
     *
     * If global dimensions are different, returns false.
     *
     * Debug mode asserts that inner fits in outer.
     */
    bool fits_between( MemoryWindow const & inner,
                       MemoryWindow const & outer ) const {
#ifndef NDEBUG
      assert( fit_in(inner, outer) );
#endif
      return ( inner.fits_in(*this)   &&
               (*this).fits_in(outer) );
    }

    /**
     * \brief Writes the internals of MemoryWindow to a string
     * \return a string value representing the internals of MemoryWindow.
     */
    inline std::string print() const {
      std::stringstream s;
      s << "Offset: " << offset_ << std::endl
        << "Extent: " << extent_ << std::endl;
      return s.str();
    }

    /**
     * \brief performs basic sanity checks to see if there is anything obviously wrong with this window.
     */
    bool sanity_check() const;

  };

  template<typename FieldType>
  class ConstFieldIterator;  // forward

  template<typename FieldType>
  class FieldIterator : public std::iterator<std::random_access_iterator_tag, typename FieldType::value_type> {
    friend class ConstFieldIterator<FieldType>;
    typedef FieldIterator<FieldType> Self;
    typedef typename FieldType::value_type AtomicType;

  public:
    FieldIterator()
    : current_(NULL),
      count_(0),
      xIndex_(0), yIndex_(0), zIndex_(0),
      yStep_(0), zStep_(0),
      xExtent_(0), yExtent_(0), zExtent_(0),
      xyExtent_(0)
    {}

    FieldIterator( AtomicType * field_values,
                   const MemoryWindow & w )
    : current_(field_values +
               w.offset(0) +
               w.offset(1) * w.glob_dim(0) +
               w.offset(2) * w.glob_dim(0) * w.glob_dim(1)),
      count_(0),
      xIndex_(0), yIndex_(0), zIndex_(0),
      yStep_(w.glob_dim(0) - w.extent(0)),
      zStep_((w.glob_dim(1) - w.extent(1)) * w.glob_dim(0)),
      xExtent_(w.extent(0)),
      yExtent_(w.extent(1)),
      zExtent_(w.extent(2)),
      xyExtent_(w.extent(0) * w.extent(1))
    {}

    //mutable dereference
    inline AtomicType & operator*() {
#     ifndef NDEBUG
      if( count_ != (xIndex_ +
          yIndex_ * xExtent_ +
          zIndex_ * xyExtent_) ){
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
            << "iterator's internal count is off";
        throw std::runtime_error(msg.str());
      }
      if( xIndex_ >= xExtent_ ||
          yIndex_ >= yExtent_ ||
          zIndex_ >= zExtent_ ||
          xIndex_ < 0 ||
          yIndex_ < 0 ||
          zIndex_ < 0 ){
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
            << "iterator is in an invalid state for dereference";
        throw std::runtime_error(msg.str());
      }
#     endif
      return *current_;
    }

    //immutable dereference
    inline AtomicType const & operator*() const {
#     ifndef NDEBUG
      if(count_ != (xIndex_ +
          yIndex_ * xExtent_ +
          zIndex_ * xyExtent_) ){
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
            << "iterator's internal count is off";
        throw std::runtime_error(msg.str());
      }
      if( xIndex_ >= xExtent_ ||
          yIndex_ >= yExtent_ ||
          zIndex_ >= zExtent_ ||
          xIndex_ < 0 ||
          yIndex_ < 0 ||
          zIndex_ < 0 ){
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
            << "iterator is in an invalid state for dereference";
        throw std::runtime_error(msg.str());
      }
#     endif
      return *current_;
    }

    //increment
    inline Self & operator++() {
      current_++; //xStep
      count_++;
      xIndex_++;
      if(xIndex_ == xExtent_) {
        current_ += yStep_; //yStep
        xIndex_ = 0;
        yIndex_++;
        if(yIndex_ == yExtent_) {
          current_ += zStep_; //zStep
          yIndex_ = 0;
          zIndex_++;
        }
      }
      return *this;
    }
    inline Self operator++(int) { Self result = *this; ++(*this); return result; }

    //decrement
    inline Self & operator--() {
      current_--; //xStep
      count_--;
      xIndex_--;
      if( xIndex_ == -1 ){
        current_ -= yStep_; //yStep
        xIndex_ = xExtent_ - 1;
        yIndex_--;
        if( yIndex_ == -1 ){
          current_ -= zStep_; //zStep
          yIndex_ = yExtent_ - 1;
          zIndex_--;
        }
      }
      return *this;
    }
    inline Self operator--(int) { Self result = *this; --(*this); return result; }

    //compound assignment
    inline Self & operator+=(int change) {
      //small change (only changes xIndex_)
      if( (change > 0 && //positive change
           change < xExtent_ - xIndex_) ||
          (change < 0 && //negative change
           -change < xIndex_) ){
        current_ += change;
        xIndex_ += change;
        count_ += change;
      }
      //bigger change (changes yIndex_ and/or zIndex_)
      else {
        current_ += (change + //xStep
            yStep_ * (((count_ + change) / xExtent_ ) - (count_ / xExtent_)) +
            zStep_ * (((count_ + change) / xyExtent_) - (count_ /xyExtent_)));
        count_ += change;
        xIndex_ = count_ % xExtent_;
        yIndex_ = (count_ % xyExtent_) / xExtent_;
        zIndex_ = count_ / xyExtent_;
      }
      return *this;
    }
    inline Self & operator-=(int change) { return *this += -change; }

    //addition/subtraction
    inline Self operator+ (int change) const { Self result = *this; result += change; return result; }
    inline Self operator- (int change) const { return *this + (-change); }

    //pointer subtraction
    inline ptrdiff_t operator- (Self const & other) const { return count_ - other.count_; }

    //offset dereference
    inline AtomicType & operator[](int change) { Self result = *this; result += change; return *result; }

    //comparisons
    inline bool operator==(Self const & other) const { return current_ == other.current_; }
    inline bool operator!=(Self const & other) const { return current_ != other.current_; }
    inline bool operator< (Self const & other) const { return current_ <  other.current_; }
    inline bool operator> (Self const & other) const { return current_ >  other.current_; }
    inline bool operator<=(Self const & other) const { return current_ <= other.current_; }
    inline bool operator>=(Self const & other) const { return current_ >= other.current_; }

    IntVec location() const{ return IntVec(xIndex_,yIndex_,zIndex_); }

  private:
    AtomicType * current_;
    int count_;
    int xIndex_, yIndex_, zIndex_;
    int yStep_, zStep_;
    int xExtent_, yExtent_, zExtent_, xyExtent_;
  };

  template<typename FieldType>
  class ConstFieldIterator : public std::iterator<std::random_access_iterator_tag, typename FieldType::value_type> {
    typedef ConstFieldIterator<FieldType> Self;
    typedef typename FieldType::value_type AtomicType;

  public:
    ConstFieldIterator()
    : current_(NULL),
      count_(0),
      xIndex_(0), yIndex_(0), zIndex_(0),
      yStep_(0), zStep_(0),
      xExtent_(0), yExtent_(0), zExtent_(0),
      xyExtent_(0)
  {}

    ConstFieldIterator(const AtomicType * const field_values,
                       const MemoryWindow & w)
    : current_(field_values +
               w.offset(0) * 1 +
               w.offset(1) * w.glob_dim(0) +
               w.offset(2) * w.glob_dim(0) * w.glob_dim(1)),
      count_(0),
      xIndex_(0), yIndex_(0), zIndex_(0),
      yStep_(w.glob_dim(0) - w.extent(0)),
      zStep_((w.glob_dim(1) - w.extent(1)) * w.glob_dim(0)),
      xExtent_(w.extent(0)),
      yExtent_(w.extent(1)),
      zExtent_(w.extent(2)),
      xyExtent_(w.extent(0) * w.extent(1))
    {}

    ConstFieldIterator(const FieldIterator<FieldType> it)
    : current_(it.current_),
      count_(it.count_),
      xIndex_(it.xIndex_),
      yIndex_(it.yIndex_),
      zIndex_(it.zIndex_),
      yStep_(it.yStep_),
      zStep_(it.zStep_),
      xExtent_(it.xExtent_),
      yExtent_(it.yExtent_),
      zExtent_(it.zExtent_),
      xyExtent_(it.xyExtent_)
    {}

    //immutable dereference
    inline AtomicType const & operator*() const
    {
#     ifndef NDEBUG
      if( count_ != (xIndex_ +
          yIndex_ * xExtent_ +
          zIndex_ * xyExtent_) ){
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
            << "iterator's internal count is off";
        throw std::runtime_error(msg.str());
      }
      if( xIndex_ >= xExtent_ ||
          yIndex_ >= yExtent_ ||
          zIndex_ >= zExtent_ ||
          xIndex_ < 0 ||
          yIndex_ < 0 ||
          zIndex_ < 0 ){
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__ << std::endl
            << "iterator is in an invalid state for dereference";
        throw std::runtime_error(msg.str());
      }
#         endif
      return *current_;
    }

    //increment
    inline Self & operator++() {
      current_++; //xStep
      count_++;
      xIndex_++;
      if( xIndex_ == xExtent_ ){
        current_ += yStep_; //yStep
        xIndex_ = 0;
        yIndex_++;
        if( yIndex_ == yExtent_ ){
          current_ += zStep_; //zStep
          yIndex_ = 0;
          zIndex_++;
        }
      }
      return *this;
    }
    inline Self operator++(int) { Self result = *this; ++(*this); return result; }

    //decrement
    inline Self & operator--() {
      --current_; //xStep
      --count_;
      --xIndex_;
      if( xIndex_ == -1 ){
        current_ -= yStep_; //yStep
        xIndex_ = xExtent_ - 1;
        yIndex_--;
        if( yIndex_ == -1 ){
          current_ -= zStep_; //zStep
          yIndex_ = yExtent_ - 1;
          --zIndex_;
        }
      }
      return *this;
    }
    inline Self operator--(int) { Self result = *this; --(*this); return result; }

    //compound assignment
    inline Self & operator+=(int change) {
      //small change (only changes xIndex_)
      if( (change > 0 && //positive change
           change < xExtent_ - xIndex_) ||
          (change < 0 && //negative change
              - change < xIndex_) ){
        current_ += change;
        xIndex_  += change;
        count_   += change;
      }
      //bigger change (changes yIndex_ and/or zIndex_)
      else {
        int new_count = count_ + change;
        int old_count = count_;
        current_ += (change + //xStep
            yStep_ * ((new_count / xExtent_) - (old_count / xExtent_)) +
            zStep_ * ((new_count / xyExtent_) - (old_count /xyExtent_)));
        count_ += change;
        xIndex_ = count_ % xExtent_;
        yIndex_ = (count_ % xyExtent_) / xExtent_;
        zIndex_ = count_ / xyExtent_;
      }
      return *this;
    }
    inline Self & operator-=(int change) { return *this += -change; }

    //addition/subtraction
    inline Self operator+ (int change) const { Self result = *this; result += change; return result; }
    inline Self operator- (int change) const { return *this + (-change); }

    //iterator subtraction
    inline ptrdiff_t operator- (Self const & other) const { return count_ - other.count_; }

    //offset dereference
    inline AtomicType & operator[](int change) { Self result = *this; result += change; return *result; }

    IntVec location() const{ return IntVec(xIndex_,yIndex_,zIndex_); }

    //comparisons
    inline bool operator==(Self const & other) const { return current_ == other.current_; }
    inline bool operator!=(Self const & other) const { return current_ != other.current_; }
    inline bool operator< (Self const & other) const { return current_ <  other.current_; }
    inline bool operator> (Self const & other) const { return current_ >  other.current_; }
    inline bool operator<=(Self const & other) const { return current_ <= other.current_; }
    inline bool operator>=(Self const & other) const { return current_ >= other.current_; }

  private:
    const AtomicType * current_;
    int count_;
    int xIndex_, yIndex_, zIndex_;
    int yStep_, zStep_;
    int xExtent_, yExtent_, zExtent_, xyExtent_;
  };

} // namespace structured
} // namespace SpatialOps

/**
 * @}
 */

#endif // SpatialOps_MemoryWindow_h
