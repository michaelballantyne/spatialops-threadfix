#ifndef SpatialOps_FieldReductions_h
#define SpatialOps_FieldReductions_h

#include <spatialops/FieldExpressions.h>
#include <spatialops/FieldExpressionsExtended.h>

//cwearl basic marcros:
#define I inline
#define S static
#define SI S I

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
    I ResultType field_fold(ResultType const & (*proc)(ResultType const &,
						       typename FieldType::value_type const &),
			    ResultType const & initialValue,
			    FieldExpression<ExprType,FieldType> & fexpr) {
    //initialize:
    ResultType result = initialValue;
    typename ExprType::template FullState<UseWholeIterator> expr = fexpr.expression().template init<UseWholeIterator>();
    //fexpr.init();
    
    //run fold:
    while(!fexpr.at_end()) {
      result = proc(result,
		    expr.eval());
      expr.next();
    };
    
    return result;
  };
  
  template<typename ResultType, typename ExprType, typename FieldType>
    I ResultType field_fold_interior(ResultType const & (*proc)(ResultType const &,
								typename FieldType::value_type const &),
				     ResultType const & initialValue,
				     FieldExpression<ExprType,FieldType> & fexpr) {
    //initialize:
    ResultType result = initialValue;
    typename ExprType::template FullState<UseInteriorIterator> expr = fexpr.expression().template init<UseInteriorIterator>();
    //fexpr.init();
    
    //run fold:
    while(!fexpr.at_end()) {
      result = proc(result,
		    expr.eval());
      expr.next();
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
    I ResultType field_fold(ResultType const & (*proc)(ResultType const &,
						       typename FieldType::value_type const &),
			    ResultType const & initialValue,
			    FieldType const & field) {
    FieldForm<FieldType> typedef ExprType;
    
    return field_fold(proc,
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
    I typename FieldType::value_type field_reduce(typename FieldType::value_type const & (*proc)(typename FieldType::value_type const &,
												 typename FieldType::value_type const &),
						  FieldExpression<ExprType,FieldType> & fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //initialize:
    typename ExprType::template FullState<UseWholeIterator> expr = fexpr.expression().template init<UseWholeIterator>();
    //fexpr.init();
    
    //set up first value:
    AtomicType result = expr.eval();
    expr.next();
    
    //run reduce:
    while(!expr.at_end()) {
      result = proc(result,
		    expr.eval());
      expr.next();
    };
    
    return result;
  };
  
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_reduce_iterior(typename FieldType::value_type const & (*proc)(typename FieldType::value_type const &,
													 typename FieldType::value_type const &),
							  FieldExpression<ExprType,FieldType> & fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //initialize:
    typename ExprType::template FullState<UseInteriorIterator> expr = fexpr.expression().template init<UseInteriorIterator>();
    //fexpr.init();
    
    //set up first value:
    AtomicType result = expr.eval();
    expr.next();
    
    //run reduce:
    while(!expr.at_end()) {
      result = proc(result,
		    expr.eval());
      expr.next();
    };
    
    return result;
  };
  
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_reduce(typename FieldType::value_type (*proc)(typename FieldType::value_type const &,
											 typename FieldType::value_type const &),
						  FieldExpression<ExprType,FieldType> & fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //initialize:
    typename ExprType::template FullState<UseWholeIterator> expr = fexpr.expression().template init<UseWholeIterator>();
    //fexpr.init();
    
    //set up first value:
    AtomicType result = expr.eval();
    expr.next();
    
    //run reduce:
    while(!expr.at_end()) {
      result = proc(result,
		    expr.eval());
      expr.next();
    };
    
    return result;
  };
  
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_reduce_interior(typename FieldType::value_type (*proc)(typename FieldType::value_type const &,
												  typename FieldType::value_type const &),
							   FieldExpression<ExprType,FieldType> & fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //initialize:
    typename ExprType::template FullState<UseInteriorIterator> expr = fexpr.expression().template init<UseInteriorIterator>();
    //fexpr.init();
    
    //set up first value:
    AtomicType result = expr.eval();
    expr.next();
    
    //run reduce:
    while(!expr.at_end()) {
      result = proc(result,
		    expr.eval());
      expr.next();
    };
    
    return result;
  };
  
  /* field_reduce for FieldType */
  template<typename FieldType>
    I typename FieldType::value_type field_reduce(typename FieldType::value_type const & (*proc)(typename FieldType::value_type const &,
												 typename FieldType::value_type const &),
						  FieldType const & field) {
    
    FieldForm<FieldType> typedef ExprType;
    
    return field_reduce(proc,
			FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* field_reduce for FieldType */
  template<typename FieldType>
    I typename FieldType::value_type field_reduce_interior(typename FieldType::value_type const & (*proc)(typename FieldType::value_type const &,
													  typename FieldType::value_type const &),
							   FieldType const & field) {
    
    FieldForm<FieldType> typedef ExprType;
    
    return field_reduce_interior(proc,
				 FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* field_reduce for FieldType */
  template<typename FieldType>
    I typename FieldType::value_type field_reduce(typename FieldType::value_type (*proc)(typename FieldType::value_type const &,
											 typename FieldType::value_type const &),
						  FieldType const & field) {
    
    FieldForm<FieldType> typedef ExprType;
    
    return field_reduce(proc,
			FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* field_reduce for FieldType */
  template<typename FieldType>
    I typename FieldType::value_type field_reduce_interior(typename FieldType::value_type (*proc)(typename FieldType::value_type const &,
												  typename FieldType::value_type const &),
							   FieldType const & field) {
    
    FieldForm<FieldType> typedef ExprType;
    
    return field_reduce_interior(proc,
				 FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* Field version of max */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_max(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run max in terms of reduce:
    return field_reduce(std::max, fexpr);
  };
  
  /* Field version of max */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_max_interior(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run max in terms of reduce:
    return field_reduce_interior(std::max, fexpr);
  };
  
  /* Field version of max */
  template<typename FieldType>
    I typename FieldType::value_type field_max(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return field_max(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* Field version of max */
  template<typename FieldType>
    I typename FieldType::value_type field_max_interior(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return field_max_interior(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* Field version of min */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_min(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run min in terms of reduce:
    return field_reduce(std::min, fexpr);
  };
  
  /* Field version of min */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_min_interior(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run min in terms of reduce:
    return field_reduce_interior(std::min, fexpr);
  };
  
  /* Field version of min */
  template<typename FieldType>
    I typename FieldType::value_type field_min(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return field_min(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* Field version of min */
  template<typename FieldType>
    I typename FieldType::value_type field_min_interior(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return field_min_interior(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  template<typename AtomicType>
    I AtomicType sum(AtomicType const & a,
		     AtomicType const & b) {
    return a + b;
  };
  
  /* Field version of sum */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_sum(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run sum in terms of reduce:
    return field_reduce(sum<AtomicType>, fexpr);
  };
  
  /* Field version of sum */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_sum_interior(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run sum in terms of reduce:
    return field_reduce_interior(sum<AtomicType>, fexpr);
  };
  
  /* Field version of sum */
  template<typename FieldType>
    I typename FieldType::value_type field_sum(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return field_sum(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* Field version of sum */
  template<typename FieldType>
    I typename FieldType::value_type field_sum_interior(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return field_sum_interior(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* Field version of norm */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_norm(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run norm in terms of reduce:
    return std::sqrt(field_sum(pow(fexpr,2)));
  };
  
  /* Field version of norm */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_norm_interior(FieldExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;
    
    //run norm in terms of reduce:
    return std::sqrt(field_sum_interior(pow(fexpr,2)));
  };
  
  /* Field version of norm */
  template<typename FieldType>
    I typename FieldType::value_type field_norm(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return field_norm(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
  /* Field version of norm */
  template<typename FieldType>
    I typename FieldType::value_type field_norm_interior(FieldType field) {
    FieldForm<FieldType> typedef ExprType;
    
    return field_norm_interior(FieldExpression<ExprType,FieldType>(ExprType(field)));
  };
  
} // namespace SpatialOps

//cwearl basic marcros:
#undef I
#undef S
#undef SI

#endif // SpatialOps_FieldReductions_h
