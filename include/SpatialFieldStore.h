#ifndef UT_SpatialFieldStore_h
#define UT_SpatialFieldStore_h

#include <SpatialField.h>

//
//  JCS: NOTE - this is not thread safe!  We could have concurrent
//  threads accessing the same memory!  I don't know how to beat this
//  elegantly.
//

namespace SpatialOps{

  // forward declaration
  template< typename T1,
	    typename T2,
	    typename T3 >     class SpatialField;


  /**
   *  @class SpatialFieldStore
   *  @author James C. Sutherland
   *  @date   May, 2007
   *
   *  The <code>SpatialFieldStore</code> class provides a mechanism to
   *  generate temporary <code>SpatialField</code> objects for use in
   *  internal operatations in the <code>SpatialField</code> class.
   *  This prevents multiple allocation/deallocation of such objects
   *  that would be required otherwise.
   *
   *  NOTE: this could get us into big trouble if we have threads
   *  running concurrently, since we would not be able to guarantee
   *  that two threads didn't use the same memory.
   *
   */
  template< typename T1,
	    typename T2,
	    typename T3 >
  class SpatialFieldStore
  {

  public:
    static SpatialFieldStore& self();

    SpatialField<T1,T2,T3>& get( const SpatialField<T1,T2,T3>& field,
				 const size_t fieldNum );

  private:

    SpatialFieldStore(){};
    ~SpatialFieldStore();

    typedef std::vector<SpatialField<T1,T2,T3>*> FieldStore;
    typedef std::map<int,FieldStore> FieldStoreMap;

    FieldStoreMap fsmap_;
  };


  //==================================================================






  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  //  Implementation
  //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





  //=================================================================


  //------------------------------------------------------------------
  template<typename T1,typename T2,typename T3>
  SpatialFieldStore<T1,T2,T3>::~SpatialFieldStore()
  {
    for( typename FieldStoreMap::iterator ii=fsmap_.begin(); ii!=fsmap_.end(); ++ii ){
      FieldStore& fs = ii->second;
      for( typename FieldStore::iterator jj=fs.begin(); jj!=fs.end(); ++jj ){
	delete *jj;
      }
    }
  }
  //------------------------------------------------------------------
  template<typename T1,typename T2,typename T3>
  SpatialFieldStore<T1,T2,T3>&
  SpatialFieldStore<T1,T2,T3>::self()
  {
    static SpatialFieldStore<T1,T2,T3> s;
    return s;
  }
  //------------------------------------------------------------------
  template<typename T1,typename T2,typename T3>
  SpatialField<T1,T2,T3>&
  SpatialFieldStore<T1,T2,T3>::get( const SpatialField<T1,T2,T3>& field,
				    const size_t fieldNum )
  {
    const int npts = field.get_ntotal();

    // find the proper map
    FieldStore& fs = fsmap_[npts];

    // check for existance of this field
    if( fieldNum >= fs.size() ){
      fs.push_back( new SpatialField<T1,T2,T3>(field) );
    }
    return *(fs[fieldNum]);
  }
  //------------------------------------------------------------------


  //====================================================================


} // namespace SpatialOps

#endif
