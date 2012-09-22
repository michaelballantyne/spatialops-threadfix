/*
 * Copyright (c) 2011 The University of Utah
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

#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/SpatialOpsDefs.h>

#include <spatialops/structured/IntVec.h>
#include <spatialops/structured/IndexTriplet.h>

#ifdef SOPS_BOOST_SERIALIZATION
# include <boost/serialization/serialization.hpp>
#endif

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
                                     const IntVec nGhostMinus=IntVec(0,0,0),
                                     const IntVec nGhostPlus=IntVec(0,0,0),
                                     const IntVec bcExtents=IntVec(0,0,0) ) const;

    /**
     *  \brief Refines/reduces the MemoryWindow into a MemoryWindow within the original.
     *
     *  \param splitPattern the number of partitions made in each ordinate direction.
     *  \param location the specific MemoryWindow to return.
     *
     *  \return MemoryWindow containing the new subwindow.
     */
    inline MemoryWindow
    refine(const IntVec splitPattern,
           const IntVec location) const {
#       ifndef NDEBUG
        for( size_t i=0; i<3; ++i ){
            assert( extent_[i] >= splitPattern[i] );
            assert( splitPattern[i] > 0 );
            assert( extent_[i] + offset_[i] <= nptsGlob_[i] );
            assert( location[i] < splitPattern[i] );
            assert( location[i] >= 0 );
        }
#       endif
        // try to make windows close to same size
        const IntVec nextra = IntVec(extent_[0] % splitPattern[0],
                                     extent_[1] % splitPattern[1],
                                     extent_[2] % splitPattern[2]);
        const IntVec stdExtent = IntVec(extent_[0] / splitPattern[0],
                                        extent_[1] / splitPattern[1],
                                        extent_[2] / splitPattern[2]);
        const IntVec baseOffset = IntVec(offset_[0] + stdExtent[0] * location[0] + (location[0] < nextra[0] ? location[0] : nextra[0]),
                                         offset_[1] + stdExtent[1] * location[1] + (location[1] < nextra[1] ? location[1] : nextra[1]),
                                         offset_[2] + stdExtent[2] * location[2] + (location[2] < nextra[2] ? location[2] : nextra[2]));
        const IntVec baseExtent = IntVec(stdExtent[0] + (location[0] < nextra[0] ? 1 : 0),
                                         stdExtent[1] + (location[1] < nextra[1] ? 1 : 0),
                                         stdExtent[2] + (location[2] < nextra[2] ? 1 : 0));

        return MemoryWindow(nptsGlob_, baseOffset, baseExtent, 0, 0, 0);
    };

    /**
     *  \brief Resizes/reduces the MemoryWindow to given number of ghost cells.
     *
     *  \param oldGhost is a GhostData that specifices how many ghost cells are currently on each face
     *  \param newGhost is a GhostData that specifices how many ghost cells are to be on each face
     *
     *  \return new MemoryWindow reduced from having oldNeg/oldPos ghost cells (and extra cells along boundary conditions) to having newNeg/newPos ghost cells on each face.
     */
    template<typename OldGhost, typename NewGhost>
    inline MemoryWindow
    resize_ghost() const {
        return MemoryWindow(nptsGlob_,
                            offset_ + OldGhost::offset()    - NewGhost::offset(),
                            extent_ - OldGhost::extent(bc_) + NewGhost::extent(bc_),
                            bc_[0],
                            bc_[1],
                            bc_[2]);
    }

    /**
     *  \brief shifts/moves the MemoryWindow by given amounts.
     *
     *  \param size is an IndexTriplet that specifies how much to shift/move the memory window.
     *
     *  \return new MemoryWindow moved by size.
     */
    template<typename IndexTriplet>
    inline MemoryWindow shift() const {
        return MemoryWindow(nptsGlob_,
                            offset_ + IndexTriplet::int_vec(),
                            extent_,
                            bc_[0],
                            bc_[1],
                            bc_[2]);
    }

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
    inline IntVec ijk_index_from_local( const int loc ) const
    {
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

    /**
     * \brief Writes the internals of MemoryWindow to a string
     * \return a string value representing the internals of MemoryWindow.
     */
    inline std::string print() const {
      std::stringstream s;
      s << "Offset: " << offset_ << std::endl
        << "Extent: " << extent_ << std::endl
        << "BC :    " << bc_     << std::endl;
      return s.str();
    }

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
      {};

      FieldIterator(AtomicType * field_values,
                    const MemoryWindow & w)
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
      {};

      //mutable dereference
      inline AtomicType & operator*() {
#         ifndef NDEBUG
          if(count_ != (xIndex_ +
                        yIndex_ * xExtent_ +
                        zIndex_ * xyExtent_)) {
              std::ostringstream msg;
              msg << __FILE__ << " : " << __LINE__ << std::endl
                  << "iterator's internal count is off";
              throw std::runtime_error(msg.str());
          };
          if(xIndex_ >= xExtent_ ||
             yIndex_ >= yExtent_ ||
             zIndex_ >= zExtent_ ||
             xIndex_ < 0 ||
             yIndex_ < 0 ||
             zIndex_ < 0) {
              std::ostringstream msg;
              msg << __FILE__ << " : " << __LINE__ << std::endl
                  << "iterator is in an invalid state for dereference";
              throw std::runtime_error(msg.str());
          };
#         endif
          return *current_;
      };

      //immutable dereference
      inline AtomicType const & operator*() const {
#         ifndef NDEBUG
          if(count_ != (xIndex_ +
                        yIndex_ * xExtent_ +
                        zIndex_ * xyExtent_)) {
              std::ostringstream msg;
              msg << __FILE__ << " : " << __LINE__ << std::endl
                  << "iterator's internal count is off";
              throw std::runtime_error(msg.str());
          };
          if(xIndex_ >= xExtent_ ||
             yIndex_ >= yExtent_ ||
             zIndex_ >= zExtent_ ||
             xIndex_ < 0 ||
             yIndex_ < 0 ||
             zIndex_ < 0) {
              std::ostringstream msg;
              msg << __FILE__ << " : " << __LINE__ << std::endl
                  << "iterator is in an invalid state for dereference";
              throw std::runtime_error(msg.str());
          };
#         endif
          return *current_;
      };

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
              };
          };
          return *this;
      };
      inline Self operator++(int) { Self result = *this; ++(*this); return result; };

      //decrement
      inline Self & operator--() {
          current_--; //xStep
          count_--;
          xIndex_--;
          if(xIndex_ == -1) {
              current_ -= yStep_; //yStep
              xIndex_ = xExtent_ - 1;
              yIndex_--;
              if(yIndex_ == -1) {
                  current_ -= zStep_; //zStep
                  yIndex_ = yExtent_ - 1;
                  zIndex_--;
              };
          };
          return *this;
      };
      inline Self operator--(int) { Self result = *this; --(*this); return result; };

      //compound assignment
      inline Self & operator+=(int change) {
          //small change (only changes xIndex_)
          if((change > 0 && //positive change
              change < xExtent_ - xIndex_) ||
             (change < 0 && //negative change
              - change < xIndex_)) {
              current_ += change;
              xIndex_ += change;
              count_ += change;
          }
          //bigger change (changes yIndex_ and/or zIndex_)
          else {
              current_ += (change + //xStep
                           yStep_ * (((count_ + change) / xExtent_) - (count_ / xExtent_)) +
                           zStep_ * (((count_ + change) / xyExtent_) - (count_ /xyExtent_)));
              count_ += change;
              xIndex_ = count_ % xExtent_;
              yIndex_ = (count_ % xyExtent_) / xExtent_;
              zIndex_ = count_ / xyExtent_;
          };
          return *this;
      };
      inline Self & operator-=(int change) { return *this += -change; };

      //addition/subtraction
      inline Self operator+ (int change) const { Self result = *this; result += change; return result; };
      inline Self operator- (int change) const { return *this + (-change); };

      //pointer subtraction
      inline ptrdiff_t operator- (Self const & other) const { return count_ - other.count_; };

      //offset dereference
      inline AtomicType & operator[](int change) { Self result = *this; result += change; return *result; };

      //comparisons
      inline bool operator==(Self const & other) const { return current_ == other.current_; };
      inline bool operator!=(Self const & other) const { return current_ != other.current_; };
      inline bool operator< (Self const & other) const { return current_ <  other.current_; };
      inline bool operator> (Self const & other) const { return current_ >  other.current_; };
      inline bool operator<=(Self const & other) const { return current_ <= other.current_; };
      inline bool operator>=(Self const & other) const { return current_ >= other.current_; };

  private:
      AtomicType * current_;
      int count_;
      int xIndex_;
      int yIndex_;
      int zIndex_;
      int yStep_;
      int zStep_;
      int xExtent_;
      int yExtent_;
      int zExtent_;
      int xyExtent_;
  };

  template<typename FieldType>
  inline FieldIterator<FieldType> operator+(int change,
                                            FieldIterator<FieldType> const & iterator) {
      return iterator + change;
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
      {};

      ConstFieldIterator(AtomicType * field_values,
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
      {};

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
      {};

      //immutable dereference
      inline AtomicType const & operator*() const {
#         ifndef NDEBUG
          if(count_ != (xIndex_ +
                        yIndex_ * xExtent_ +
                        zIndex_ * xyExtent_)) {
              std::ostringstream msg;
              msg << __FILE__ << " : " << __LINE__ << std::endl
                  << "iterator's internal count is off";
              throw std::runtime_error(msg.str());
          };
          if(xIndex_ >= xExtent_ ||
             yIndex_ >= yExtent_ ||
             zIndex_ >= zExtent_ ||
             xIndex_ < 0 ||
             yIndex_ < 0 ||
             zIndex_ < 0) {
              std::ostringstream msg;
              msg << __FILE__ << " : " << __LINE__ << std::endl
                  << "iterator is in an invalid state for dereference";
              throw std::runtime_error(msg.str());
          };
#         endif
          return *current_;
      };

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
              };
          };
          return *this;
      };
      inline Self operator++(int) { Self result = *this; ++(*this); return result; };

      //decrement
      inline Self & operator--() {
          current_--; //xStep
          count_--;
          xIndex_--;
          if(xIndex_ == -1) {
              current_ -= yStep_; //yStep
              xIndex_ = xExtent_ - 1;
              yIndex_--;
              if(yIndex_ == -1) {
                  current_ -= zStep_; //zStep
                  yIndex_ = yExtent_ - 1;
                  zIndex_--;
              };
          };
          return *this;
      };
      inline Self operator--(int) { Self result = *this; --(*this); return result; };

      //compound assignment
      inline Self & operator+=(int change) {
          //small change (only changes xIndex_)
          if((change > 0 && //positive change
              change < xExtent_ - xIndex_) ||
             (change < 0 && //negative change
              - change < xIndex_)) {
              current_ += change;
              xIndex_ += change;
              count_ += change;
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
          };
          return *this;
      };
      inline Self & operator-=(int change) { return *this += -change; };

      //addition/subtraction
      inline Self operator+ (int change) const { Self result = *this; result += change; return result; };
      inline Self operator- (int change) const { return *this + (-change); };

      //iterator subtraction
      inline ptrdiff_t operator- (Self const & other) const { return count_ - other.count_; };

      //offset dereference
      inline AtomicType & operator[](int change) { Self result = *this; result += change; return *result; };

      //comparisons
      inline bool operator==(Self const & other) const { return current_ == other.current_; };
      inline bool operator!=(Self const & other) const { return current_ != other.current_; };
      inline bool operator< (Self const & other) const { return current_ <  other.current_; };
      inline bool operator> (Self const & other) const { return current_ >  other.current_; };
      inline bool operator<=(Self const & other) const { return current_ <= other.current_; };
      inline bool operator>=(Self const & other) const { return current_ >= other.current_; };

  private:
      AtomicType * current_;
      int count_;
      int xIndex_;
      int yIndex_;
      int zIndex_;
      int yStep_;
      int zStep_;
      int xExtent_;
      int yExtent_;
      int zExtent_;
      int xyExtent_;
  };

  template<typename FieldType>
  inline ConstFieldIterator<FieldType> operator+(int change,
                                                 ConstFieldIterator<FieldType> const & iterator) {
      return iterator + change;
  };

} // namespace structured
} // namespace SpatialOps

/**
 * @}
 */

#endif // SpatialOps_MemoryWindow_h
