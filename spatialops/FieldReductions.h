#ifndef SpatialOps_FieldReductions_h
#define SpatialOps_FieldReductionss_h

#include <spatialops/FieldOperations.h>
#include <spatialops/FieldOperationDefinitions.h>

namespace SpatialOps{

  /**
   *  @author Christopher Earl
   *  @date January, 2011
   *  
   *  @brief Definition of fold for FieldExpression. (FieldExpression case.)
   *  
   *  @relates FieldExpression
   *  
   *  \tparam ResultType Type returned from \c proc.
   *  \tparam ExprType A FieldExpression-style Standard Type.
   *  \tparam FieldType Field type.
   *  
   *  \param proc A function pointer with type \code ResultType (*)(ResultType const &, AtomicType const &) \endcode where \c AtomicType is \c typename \c FieldType::value_type.
   *  \param initialValue A ResultType to start the fold.
   *  \param fexpr A FieldExpression object to run fold over.
   *
   *  \return result of folding \c proc over \c fexpr using \c initialValue.
   *  
   *  Fold works as follows:
   *  
   *  1) \c initialValue is used as the first intermediate result.
   *  
   *  2) An iterator to the first element of \c fexpr is set up.
   *  
   *  3) \c proc is called with the itermediate result and the current element of \c fexpr.
   *  
   *  4) The result of this call is the new intermediate result.
   *  
   *  5a) If there is an element after the current one in \c fexpr, it becomes the new current element.
   *  (The iterator is incremented.)
   *  Repeat steps 3-5.
   *  
   *  5b) If the iterator is at the end of \c fexpr, the intermediate result is returned as the final one.
   */
  template<typename ResultType, typename ExprType, typename FieldType>
    inline ResultType fieldFold(ResultType const & (*proc)(ResultType const &,
							   typename FieldType::value_type const &),
				ResultType const & initialValue,
				FieldExpression<ExprType,FieldType> & fexpr) {
    //initialize:
    ResultType result = initialValue;
    fexpr.init();
    
    //run fold:
    while(!fexpr.at_end()) {
      result = proc(result,
		    fexpr.eval());
      fexpr.next();
    };
    
    return result;
  };
  
  /**
   *  @author Christopher Earl
   *  @date January, 2011
   *  
   *  @brief Definition of fold for FieldExpression. (Field case.)
   *  
   *  @relates FieldExpression
   *  
   *  \tparam ResultType Type returned from \c proc.
   *  \tparam FieldType Field type.
   *  
   *  \param proc A function pointer with type \code ResultType (*)(ResultType const &, AtomicType const &) \endcode where \c AtomicType is \c typename \c FieldType::value_type.
   *  \param initialValue A ResultType to start the fold.
   *  \param field A Field object.
   *
   *  \return result of folding \c proc over \c field using \c initialValue.
   *  
   *  This instance of \c fold wraps \c field in a FieldExpression and then calls itself.
   */
  template<typename ResultType, typename FieldType>
    inline ResultType fieldFold(ResultType const & (*proc)(ResultType const &,
							   typename FieldType::value_type const &),
				ResultType const & initialValue,
				FieldType const & field) {
    FieldForm<FieldType> typedef ExprType;
    
    return fieldFold(proc,
		     initialValue,
		     FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /**
   *  @author Christopher Earl
   *  @date January, 2011
   *  
   *  @brief Definition of reduce for FieldExpression. (FieldExpression case.)
   *  
   *  @relates FieldExpression
   *  
   *  \tparam ExprType A FieldExpression-style Standard Type.
   *  \tparam FieldType Field type.
   *  
   *  \param proc A function pointer with type \code AtomicType (*)(AtomicType const &, AtomicType const &) \endcode where \c AtomicType is \c typename \c FieldType::value_type.
   *  \param fexpr A FieldExpression object to run reduce over.
   *
   *  \return result of folding \c proc over \c fexpr.
   *  
   *  Reduce works similarly to fold.
   *  Reduce works as follows:
   *  
   *  1) The first value in \c fexpr is used as the first intermediate result.
   *  
   *  2) An iterator to the second element of \c fexpr is set up.
   *  
   *  3) \c proc is called with the itermediate result and the current element of \c fexpr.
   *  
   *  4) The result of this call is the new intermediate result.
   *  
   *  5a) If there is an element after the current one in \c fexpr, it becomes the new current element.
   *  (The iterator is incremented.)
   *  Repeat steps 3-5.
   *  
   *  5b) If the iterator is at the end of \c fexpr, the intermediate result is returned as the final one.
   *
   *  \warning Assumes fexpr has at least one element. (Undefined behavior if fexpr has zero elements.)
   */
  template<typename ExprType, typename FieldType>
    inline typename FieldType::value_type fieldReduce(typename FieldType::value_type const & (*proc)(typename FieldType::value_type const &,
												     typename FieldType::value_type const &),
						      FieldExpression<ExprType,FieldType> & fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //initialize:
    fexpr.init();
    
    //set up first value:
    AtomicType result = fexpr.eval();
    fexpr.next();
    
    //run reduce:
    while(!fexpr.at_end()) {
      result = proc(result,
		    fexpr.eval());
      fexpr.next();
    };
    
    return result;
  };
  
  template<typename ExprType, typename FieldType>
    inline typename FieldType::value_type fieldReduce(typename FieldType::value_type (*proc)(typename FieldType::value_type const &,
											     typename FieldType::value_type const &),
						      FieldExpression<ExprType,FieldType> & fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //initialize:
    fexpr.init();
    
    //set up first value:
    AtomicType result = fexpr.eval();
    fexpr.next();
    
    //run reduce:
    while(!fexpr.at_end()) {
      result = proc(result,
		    fexpr.eval());
      fexpr.next();
    };
    
    return result;
  };
  
  /* fieldReduce for FieldType */
  template<typename FieldType>
    inline typename FieldType::value_type fieldReduce(typename FieldType::value_type const & (*proc)(typename FieldType::value_type const &,
												     typename FieldType::value_type const &),
						      FieldType const & field) {
    
    FieldForm<FieldType> typedef ExprType;
    
    return fieldReduce(proc,
		       FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* fieldReduce for FieldType */
  template<typename FieldType>
    inline typename FieldType::value_type fieldReduce(typename FieldType::value_type (*proc)(typename FieldType::value_type const &,
											     typename FieldType::value_type const &),
						      FieldType const & field) {
    
    FieldForm<FieldType> typedef ExprType;
    
    return fieldReduce(proc,
		       FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* Field version of max */
  template<typename ExprType, typename FieldType>
    inline typename FieldType::value_type fieldMax(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run max in terms of reduce:
    return fieldReduce(std::max, fexpr);
  };
  
  /* Field version of max */
  template<typename FieldType>
    inline typename FieldType::value_type fieldMax(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return fieldMax(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* Field version of min */
  template<typename ExprType, typename FieldType>
    inline typename FieldType::value_type fieldMin(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run min in terms of reduce:
    return fieldReduce(std::min, fexpr);
  };
  
  /* Field version of min */
  template<typename FieldType>
    inline typename FieldType::value_type fieldMin(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return fieldMin(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  template<typename AtomicType>
    inline AtomicType sum(AtomicType const & a,
			  AtomicType const & b) {
    return a + b;
  };
  
  /* Field version of sum */
  template<typename ExprType, typename FieldType>
    inline typename FieldType::value_type fieldSum(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run sum in terms of reduce:
    return fieldReduce(sum<AtomicType>, fexpr);
  };
  
  /* Field version of sum */
  template<typename FieldType>
    inline typename FieldType::value_type fieldSum(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return fieldSum(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* Field version of norm */
  template<typename ExprType, typename FieldType>
    inline typename FieldType::value_type fieldNorm(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run norm in terms of reduce:
    return std::sqrt(fieldSum(pow(fexpr,2)));
  };
  
  /* Field version of norm */
  template<typename FieldType>
    inline typename FieldType::value_type fieldNorm(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return fieldNorm(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
} // namespace SpatialOps

#endif // SpatialOps_FieldReductions_h
