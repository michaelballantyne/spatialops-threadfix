#ifndef SpatialOpDatabase_h
#define SpatialOpDatabase_h

#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>

#include <map>
#include <list>
#include <stdexcept>
#include <sstream>

namespace SpatialOps{

/**
 *  @class OperatorDatabase
 *  @author James C. Sutherland
 *  @date  July, 2008
 *  @brief Provides a database to hold operators of any type.
 */
class OperatorDatabase
{

  template<typename OpT> class OpStore
  {
    typedef std::map<int,OpT*> IDOpMap;
    IDOpMap opMap_;
    int idCounter_;
  public:
    OpStore();
    ~OpStore();
    int add( OpT* const op );
    OpT* get( const int id=-1 ) const;
  }; // class OpStore

  typedef std::list< boost::any > StoreList;
  StoreList storeList_;

  template<typename OpT>
  OpStore<OpT>* find() const;

public:

  OperatorDatabase(){}
  ~OperatorDatabase(){}

  /**
   *  \brief Register a new operator of the given type
   *  \param op The operator to register.  Ownership is transfered.
   *         This should be heap-allocated via "new".
   *
   *  Example: <code>opDB.register_new_operator( new MyOpType(myOpAssembler) );</code>
   *
   *  \todo Consider passing the operator assembler here rather than
   *        the operator itself.  Then we could build the operator
   *        internally.  Note that then we would have to provide the
   *        template argument when calling this method, whereas it can
   *        be deduced currently.
   */
  template<typename OpT>
  inline int register_new_operator( OpT* const op );

  /**
   *  \brief Retrieve an operator of the requested type.
   *  \param id In the case where more than one operator of this type
   *         has been registered, this specifies the identifier for
   *         the desired operator.  This argument may be ommitted.
   *         However, if omitted and multiple operators of this type
   *         are registered, an exception will be thrown.
   */
  template<typename OpT>
  inline OpT* retrieve_operator( const int id=-1 ) const;

  /**
   *  \brief Remove all operators of all types
   */
  inline void empty_database();

}; // class OperatorDatabase


//====================================================================





// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//                          Implementation
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





//====================================================================


template<typename OpT>
typename OperatorDatabase::OpStore<OpT>*
OperatorDatabase::find() const
{
  typedef typename boost::shared_ptr< OpStore<OpT> > PtrT;
  for( typename StoreList::const_iterator istor=storeList_.begin(); istor!=storeList_.end(); ++istor ){
    if( typeid(PtrT) == istor->type() ){
      return boost::any_cast<PtrT>(*istor).get();
    }
  }
  return NULL;
}

//--------------------------------------------------------------------

template<typename OpT>
int
OperatorDatabase::register_new_operator( OpT* const op )
{
  // see if we have a store built for this yet.
  OpStore<OpT>* store = this->find<OpT>();
  if( store==NULL ){
    typedef boost::shared_ptr< OpStore<OpT> > PtrT;
    store = new OpStore<OpT>();
    storeList_.push_back( PtrT(store) );
  }
  // add the operator 
  return store->add(op);
}

//--------------------------------------------------------------------

template<typename OpT>
OpT*
OperatorDatabase::retrieve_operator( const int id ) const
{
  const OpStore<OpT>* const store = this->find<OpT>();
  if( store==NULL ){
    std::ostringstream msg;
    msg << "ERROR!  Attempted to retrieve an operator that does not exist." << std::endl
        << "        Operator type name: " << typeid(OpT).name() << std::endl;
    throw std::runtime_error( msg.str() );
  }
  return store->get(id);
}

//--------------------------------------------------------------------

void
OperatorDatabase::empty_database()
{
  storeList_.clear();
}

//--------------------------------------------------------------------

template<typename OpT>
OperatorDatabase::OpStore<OpT>::OpStore()
  : idCounter_(0)
{}

template<typename OpT>
OperatorDatabase::OpStore<OpT>::~OpStore()
{
  for( typename IDOpMap::iterator i=opMap_.begin(); i!=opMap_.end(); ++i ){
    delete i->second;
  }
}

//--------------------------------------------------------------------

template<typename OpT>
int
OperatorDatabase::OpStore<OpT>::add( OpT* const op )
{
  opMap_.insert( make_pair( ++idCounter_, op ) );
  return idCounter_;
}

//--------------------------------------------------------------------

template<typename OpT>
OpT*
OperatorDatabase::OpStore<OpT>::get( const int id ) const
{
  typename IDOpMap::const_iterator iop = opMap_.begin();
  if( id==-1 ){
    if( opMap_.size() > 1 ){
      std::ostringstream msg;
      msg << "ERROR!  You must provide a unique identifier, since multiple operators have been registered." << std::endl
          << "        registered ids:" << std::endl;
      for( iop=opMap_.begin(); iop!=opMap_.end(); ++iop ){
        msg << "     " << iop->first << std::endl;
      }
      throw std::runtime_error(msg.str());
    }
  }
  else{
    iop = opMap_.find( id );
  }
  if( iop == opMap_.end() ){
    std::ostringstream msg;
    msg << "ERROR!  Attempted to retrieve an operator that does not exist." << std::endl
        << "        Operator type name: " << typeid(OpT).name() << std::endl
        << "   Registered ids:" << std::endl;
    for( iop=opMap_.begin(); iop!=opMap_.end(); ++iop ){
      msg << "     " << iop->first <<  std::endl;
    }
    throw std::runtime_error( msg.str() );
  }
  return iop->second;
}

//--------------------------------------------------------------------

} // Namespace SpatialOps

#endif
