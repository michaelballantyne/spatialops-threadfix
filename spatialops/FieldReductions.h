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
   *  @brief Definition of fold for NeboExpression. (NeboExpression case.)
   *
   *  @relates NeboExpression
   *
   *  \tparam ResultType Type returned from \c proc.
   *  \tparam ExprType A NeboExpression-style Standard Type.
   *  \tparam FieldType Field type.
   *
   *  \param proc A function pointer with type \code ResultType (*)(ResultType const &, AtomicType const &) \endcode where \c AtomicType is \c typename \c FieldType::value_type.
   *  \param initialValue A ResultType to start the fold.
   *  \param fexpr A NeboExpression object to run fold over.
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
			    NeboExpression<ExprType,FieldType> & fexpr) {
    //initialize:
    ResultType result = initialValue;
    typename ExprType::template Iterator<UseWholeIterator>::SeqWalkType expr = fexpr.expr().template init<UseWholeIterator>();
    //typename ExprType::template FullState<UseWholeIterator> expr = fexpr.expression().template init<UseWholeIterator>();
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
				     NeboExpression<ExprType,FieldType> & fexpr) {
    //initialize:
    ResultType result = initialValue;
    typename ExprType::template Iterator<UseInteriorIterator>::SeqWalkType expr = fexpr.expr().template init<UseInteriorIterator>();
    //typename ExprType::template FullState<UseInteriorIterator> expr = fexpr.expression().template init<UseInteriorIterator>();
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
   *  @brief Definition of fold for NeboExpression. (NeboField case.)
   *
   *  @relates NeboExpression
   *
   *  \tparam ResultType Type returned from \c proc.
   *  \tparam FieldType Field type.
   *
   *  \param proc A function pointer with type \code ResultType (*)(ResultType const &, AtomicType const &) \endcode where \c AtomicType is \c typename \c FieldType::value_type.
   *  \param initialValue A ResultType to start the fold.
   *  \param field A NeboField object.
   *
   *  \return result of folding \c proc over \c field using \c initialValue.
   *
   *  This instance of \c fold wraps \c field in a NeboExpression and then calls itself.
   */
  template<typename ResultType, typename FieldType>
    I ResultType field_fold(ResultType const & (*proc)(ResultType const &,
						       typename FieldType::value_type const &),
			    ResultType const & initialValue,
			    FieldType const & field) {
    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_fold(proc,
                      initialValue,
                      NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /**
   *  @author Christopher Earl
   *  @date January, 2011
   *
   *  @brief Definition of reduce for NeboExpression. (NeboExpression case.)
   *
   *  @relates NeboExpression
   *
   *  \tparam ExprType A NeboExpression-style Standard Type.
   *  \tparam FieldType Field type.
   *
   *  \param proc A function pointer with type \code AtomicType (*)(AtomicType const &, AtomicType const &) \endcode where \c AtomicType is \c typename \c FieldType::value_type.
   *  \param fexpr A NeboExpression object to run reduce over.
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
						  NeboExpression<ExprType,FieldType> & fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //initialize:
    typename ExprType::template Iterator<UseWholeIterator>::SeqWalkType expr = fexpr.expr().template init<UseWholeIterator>();
    //typename ExprType::template FullState<UseWholeIterator> expr = fexpr.expression().template init<UseWholeIterator>();
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
    I typename FieldType::value_type field_reduce_interior(typename FieldType::value_type const & (*proc)(typename FieldType::value_type const &,
                                                                                                          typename FieldType::value_type const &),
                                                           NeboExpression<ExprType,FieldType> & fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //initialize:
    typename ExprType::template Iterator<UseInteriorIterator>::SeqWalkType expr = fexpr.expr().template init<UseInteriorIterator>();
    //typename ExprType::template FullState<UseInteriorIterator> expr = fexpr.expression().template init<UseInteriorIterator>();
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
						  NeboExpression<ExprType,FieldType> & fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //initialize:
    typename ExprType::template Iterator<UseWholeIterator>::SeqWalkType expr = fexpr.expr().template init<UseWholeIterator>();
    //typename ExprType::template FullState<UseWholeIterator> expr = fexpr.expression().template init<UseWholeIterator>();
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
							   NeboExpression<ExprType,FieldType> & fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //initialize:
    typename ExprType::template Iterator<UseInteriorIterator>::SeqWalkType expr = fexpr.expr().template init<UseInteriorIterator>();
    //typename ExprType::template FullState<UseInteriorIterator> expr = fexpr.expression().template init<UseInteriorIterator>();
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

    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_reduce(proc,
			NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /* field_reduce for FieldType */
  template<typename FieldType>
    I typename FieldType::value_type field_reduce_interior(typename FieldType::value_type const & (*proc)(typename FieldType::value_type const &,
													  typename FieldType::value_type const &),
							   FieldType const & field) {

    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_reduce_interior(proc,
				 NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /* field_reduce for FieldType */
  template<typename FieldType>
    I typename FieldType::value_type field_reduce(typename FieldType::value_type (*proc)(typename FieldType::value_type const &,
											 typename FieldType::value_type const &),
						  FieldType const & field) {

    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_reduce(proc,
			NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /* field_reduce for FieldType */
  template<typename FieldType>
    I typename FieldType::value_type field_reduce_interior(typename FieldType::value_type (*proc)(typename FieldType::value_type const &,
												  typename FieldType::value_type const &),
							   FieldType const & field) {

    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_reduce_interior(proc,
				 NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /* Field version of max */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_max(NeboExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //run max in terms of reduce:
    return field_reduce(std::max, fexpr);
  };

  /* Field version of max */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_max_interior(NeboExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //run max in terms of reduce:
    return field_reduce_interior(std::max, fexpr);
  };

  /* Field version of max */
  template<typename FieldType>
    I typename FieldType::value_type field_max(FieldType field) {
    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_max(NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /* Field version of max */
  template<typename FieldType>
    I typename FieldType::value_type field_max_interior(FieldType field) {
    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_max_interior(NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /* Field version of min */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_min(NeboExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //run min in terms of reduce:
    return field_reduce(std::min, fexpr);
  };

  /* Field version of min */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_min_interior(NeboExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //run min in terms of reduce:
    return field_reduce_interior(std::min, fexpr);
  };

  /* Field version of min */
  template<typename FieldType>
    I typename FieldType::value_type field_min(FieldType field) {
    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_min(NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /* Field version of min */
  template<typename FieldType>
    I typename FieldType::value_type field_min_interior(FieldType field) {
    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_min_interior(NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  template<typename AtomicType>
    I AtomicType sum(AtomicType const & a,
		     AtomicType const & b) {
    return a + b;
  };

  /* Field version of sum */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_sum(NeboExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //run sum in terms of reduce:
    return field_reduce(sum<AtomicType>, fexpr);
  };

  /* Field version of sum */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_sum_interior(NeboExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //run sum in terms of reduce:
    return field_reduce_interior(sum<AtomicType>, fexpr);
  };

  /* Field version of sum */
  template<typename FieldType>
    I typename FieldType::value_type field_sum(FieldType field) {
    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_sum(NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /* Field version of sum */
  template<typename FieldType>
    I typename FieldType::value_type field_sum_interior(FieldType field) {
    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_sum_interior(NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /* Field version of norm */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_norm(NeboExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //run norm in terms of reduce:
    return std::sqrt(field_sum(pow(fexpr,2)));
  };

  /* Field version of norm */
  template<typename ExprType, typename FieldType>
    I typename FieldType::value_type field_norm_interior(NeboExpression<ExprType,FieldType> fexpr) {
    typename FieldType::value_type typedef AtomicType;

    //run norm in terms of reduce:
    return std::sqrt(field_sum_interior(pow(fexpr,2)));
  };

  /* Field version of norm */
  template<typename FieldType>
    I typename FieldType::value_type field_norm(FieldType field) {
    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_norm(NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

  /* Field version of norm */
  template<typename FieldType>
    I typename FieldType::value_type field_norm_interior(FieldType field) {
    NeboConstField<Initial, FieldType> typedef ExprType;

    return field_norm_interior(NeboExpression<ExprType,FieldType>(ExprType(field)));
  };

} // namespace SpatialOps

//cwearl basic marcros:
#undef I
#undef S
#undef SI

#endif // SpatialOps_FieldReductions_h
