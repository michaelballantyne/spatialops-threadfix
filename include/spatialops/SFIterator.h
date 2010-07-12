#ifndef SFIterator_h
#define SFIterator_h

#include <spatialops/SpatialOpsConfigure.h>

#include <iterator>


namespace SpatialOps{

/**
 *  @class  InteriorIterator
 *  @author James C. Sutherland
 *
 *  @brief Provides a forward iterator for SpatialField that omits
 *  ghost cells.
 */
template<typename ForwardIterator>
class InteriorIterator
{
public:

  typedef InteriorIterator<ForwardIterator> self;

  // required typedefs for stl compliance
  typedef typename std::iterator_traits<ForwardIterator>::value_type      value_type;
  typedef typename std::iterator_traits<ForwardIterator>::reference       reference;
  typedef typename std::iterator_traits<ForwardIterator>::pointer         pointer;
  typedef typename std::iterator_traits<ForwardIterator>::difference_type difference_type;
  typedef          std::forward_iterator_tag                              iterator_category;


  InteriorIterator( ForwardIterator first, const size_t ix,
                    std::set<size_t>::const_iterator ighost,
                    const std::set<size_t>::const_iterator ighostend );
  InteriorIterator( const self& other );

  // required operators for a forward iterator: ++ (prefix and postfix), *, ==, !=
  //@{
  inline self& operator++();
  inline self  operator++(int);

  inline       reference operator*();
  inline const reference operator*() const;

  inline bool operator==(const self& other) const;
  inline bool operator!=(const self& other) const;

  inline self& operator=(const self& other);
  //@}

private:
  ForwardIterator current_;
  size_t ix_;
  std::set<size_t>::const_iterator ig_;
  std::set<size_t>::const_iterator ige_; // Can't declare const because of the assignment operator.
};


//====================================================================


/**
 *  @class  GhostIterator
 *  @author James C. Sutherland
 *
 *  @brief Provides a forward iterator for SpatialField that omits
 *  interior cells, cycling only through ghost cells.
 *
 *  This has not been thoroughly tested...
 */
template<typename ForwardIterator>
class GhostIterator
{
public:

  typedef GhostIterator<ForwardIterator> self;

  // required typedefs
  typedef typename std::iterator_traits<ForwardIterator>::value_type      value_type;
  typedef typename std::iterator_traits<ForwardIterator>::reference       reference;
  typedef typename std::iterator_traits<ForwardIterator>::pointer         pointer;
  typedef typename std::iterator_traits<ForwardIterator>::difference_type difference_type;
  typedef          std::forward_iterator_tag                              iterator_category;


  GhostIterator( ForwardIterator first, const size_t ix,
                 std::set<size_t>::const_iterator ig,
                 const std::set<size_t>::const_iterator ige );
  GhostIterator( const GhostIterator& other );

  // required operators for a forward iterator: ++ (prefix and postfix), *, ==, !=
  //@{
  inline self& operator++();
  inline self  operator++(int);

  inline reference operator*();

  inline bool operator==(const self& other) const;
  inline bool operator!=(const self& other) const;

  inline self& operator=(const self& other);
  //@}

private:
  ForwardIterator current_;
  size_t ix_;
  std::set<size_t>::const_iterator ig_;
  const std::set<size_t>::const_iterator ige_;
};


//--------------------------------------------------------------------
template<typename ForwardIterator>
InteriorIterator<ForwardIterator>::
InteriorIterator( ForwardIterator fi,
                  const size_t ix,
                  std::set<size_t>::const_iterator ig,
                  const std::set<size_t>::const_iterator ige )
  : current_ ( fi ),
    ix_      ( ix ),
    ig_      ( ig ),
    ige_     ( ige)
{
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
InteriorIterator<ForwardIterator>::
InteriorIterator( const self& other )
  : current_ ( other.current_  ),
    ix_      ( other.ix_       ),
    ig_      ( other.ig_       ),
    ige_     ( other.ige_      )
{
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
InteriorIterator<ForwardIterator>&
InteriorIterator<ForwardIterator>::
operator=(const self& other)
{
  current_ = other.current_;
  ix_      = other.ix_;
  ig_      = other.ig_;
  ige_     = other.ige_;
  return *this;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
InteriorIterator<ForwardIterator>&
InteriorIterator<ForwardIterator>::
operator++()
{
  ++ix_; ++current_;
  while( *ig_ == ix_ && ig_!=ige_ ){ ++ix_; ++ig_; ++current_; }
  return *this;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
InteriorIterator<ForwardIterator>
InteriorIterator<ForwardIterator>::
operator++(int)
{
  self tmp( *this );
  ++(*this);
  return tmp;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
typename InteriorIterator<ForwardIterator>::reference
InteriorIterator<ForwardIterator>::
operator*()
{
  return *current_;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
const typename InteriorIterator<ForwardIterator>::reference
InteriorIterator<ForwardIterator>::
operator*() const
{
  return *current_;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
bool
InteriorIterator<ForwardIterator>::
operator!=(const self& other) const
{
  return current_ != other.current_;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
bool
InteriorIterator<ForwardIterator>::
operator==(const self& other) const
{
  return current_ == other.current_;
}
//--------------------------------------------------------------------


//====================================================================


//--------------------------------------------------------------------
template<typename ForwardIterator>
GhostIterator<ForwardIterator>::
GhostIterator( ForwardIterator fi,
               const size_t ix,
               std::set<size_t>::const_iterator ig,
               const std::set<size_t>::const_iterator ige )
  : current_ ( fi ),
    ix_      ( ix ),
    ig_      ( ig ),
    ige_     ( ige)
{
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
GhostIterator<ForwardIterator>::
GhostIterator( const self& other )
  : current_ ( other.current_  ),
    ix_      ( other.ix_       ),
    ig_      ( other.ig_       ),
    ige_     ( other.ige_      )
{
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
GhostIterator<ForwardIterator>&
GhostIterator<ForwardIterator>::
operator=(const self& other)
{
  current_ = other.current_;
  ix_      = other.ix_;
  ig_      = other.ig_;
  ige_     = other.ige_;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
GhostIterator<ForwardIterator>&
GhostIterator<ForwardIterator>::
operator++()
{
  // jcs this may be wrong:
  ++ix_; ++current_;
  while( *ig_ != ix_ && ig_!=ige_ ){ ++ix_; ++ig_; ++current_; }
  return *this;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
GhostIterator<ForwardIterator>
GhostIterator<ForwardIterator>::
operator++(int)
{
  self tmp( *this );
  ++(*this);
  return tmp;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
typename GhostIterator<ForwardIterator>::reference
GhostIterator<ForwardIterator>::
operator*()
{
  return *current_;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
bool
GhostIterator<ForwardIterator>::
operator!=(const self& other) const
{
  return current_ != other.current_;
}
//--------------------------------------------------------------------
template<typename ForwardIterator>
bool
GhostIterator<ForwardIterator>::
operator==(const self& other) const
{
  return current_ == other.current_;
}
//--------------------------------------------------------------------



}


#endif
