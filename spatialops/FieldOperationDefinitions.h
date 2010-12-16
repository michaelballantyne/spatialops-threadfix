#ifndef SpatialOps_FieldOperationDefinitions_h
#define SpatialOps_FieldOperationDefinitions_h

#include "FieldOperations.h"
#include <cmath>

namespace SpatialOps{

  /**
   *  @file FieldOperationsDefinitions.h
   */

  /**
   *  @page expressiontemplates
   *  @ingroup ExpressionTemplates
   *  @par Supported operations
   *
   *  The following expressions are supported.  Here, \c a \c b and \c
   *  c can be other valid expressions.
   *  \code
   *    a <<= sin(b); 
   *    a <<= cos(b);
   *    a <<= tan(c);
   *    a <<= tanh(c);
   *    a <<= exp(c);
   *  \endcode
   */
  
  /**
   *  @struct SumOp
   *  @author Christopher Earl
   *  @date October, 2010
   *  @ingroup ExpressionTemplatesDetail
   *  
   *  @brief Specialized FieldExpression-style representation of the addition operator, \c operator \c +.
   *  
   *  \tparam Operand1 First operand's type.
   *  \tparam Operand2 Second operand's type.
   *  \tparam FieldType Field type.
   *  
   *  SumOp allows the use of \c operator \c + in FieldExpression objects for addition/summation.
   *  That is, for any legal FieldExpression objects, \c operand1 and \c operand2, the following is a legal FieldExpression:
   *  \code
   *  operand1 + operand2
   *  \endcode
   *  When \c operator \c <<= is used to compute/assign the above SumOp, the successive elements of \c operand1 and \c operand2 will be computed, added together, and returned/assigned without a temporary FieldType object being used as an intermediary in the computation.
   *  
   *  \warning SumOp, along with the other SpecificBinOp operators, uses the standard C++ operator precedence.
   *  Thus, the two following FieldExpression objects return different results (assuming \c a is a FieldType object and using ProdOp):
   *  \code
   *  a + 4 * 6
   *  \endcode
   *  and
   *  \code
   *  (a + 4) * 6
   *  \endcode
   *  
   *  
   *  
   *  @typedef FieldType::value_type typedef SumOp::AtomicType
   *  @brief Typedef for AtomicType.
   *  
   *  Typedef for AtomicType based on FieldType's value_type.
   *  
   *  
   *  
   *  @fn SumOp::SumOp(Operand1 op1, Operand2 op2)
   *  @brief Constructor for SumOp.
   *  
   *  \param op1 an Operand1 object.
   *  \param op2 an Operand2 object.
   *  \return a new SumOp containing op1 and op2.
   *  
   *  Construct a new SumOp containing op1 and op2.
   *  
   *  
   *  
   *  @fn SumOp::first() const
   *  @brief Operand1's accessor
   *  
   *  \return constant reference to op1 (first operand passed to constructor).
   *  
   *  Returns constant reference to first operand.
   *  
   *  
   *  
   *  @fn SumOp::second() const
   *  @brief Operand2's accessor
   *  
   *  \return constant reference to op2 (second operand passed to constructor).
   *  
   *  Returns constant reference to second operand.
   *  
   *  
   *  
   *  @fn void SumOp::init()
   *  @brief Initializer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Initializes operands.
   *  
   *  
   *  
   *  @fn void SumOp::next()
   *  @brief Incrementer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Increments operands.
   *  
   *  
   *  
   *  @fn AtomicType SumOp::eval() const
   *  @brief Computes and returns current element.
   *  
   *  \return current element: sum of operands' current elements.
   *  
   *  Computes and returns sum of operands' current elements.
   *  
   *  
   *  
   *  @fn bool SumOp::at_end() const
   *  @brief Predicate: Current position is the end?
   *  
   *  \return Boolean; true, if currently at end; false; if not.
   *  
   *  Returns whether or not either operands' current position is the end/last position.
   *  
   *  
   *  
   *  @fn bool SumOp::has_length() const
   *  @brief Predicate: Can reach end position?
   *  
   *  \return Boolean; true, if one or both operands can reach the end.
   *  
   *  Returns whether or not either operands can reach the end position.
   */
  BUILD_BINARY_OPERATOR(SumOp, +, operator +);
  
  /**
   *  @struct DiffOp
   *  @author Christopher Earl
   *  @date October, 2010
   *  @ingroup ExpressionTemplatesDetail
   *  
   *  @brief Specialized FieldExpression-style representation of the subtraction operator, \c operator \c -.
   *  
   *  \tparam Operand1 First operand's type.
   *  \tparam Operand2 Second operand's type.
   *  \tparam FieldType Field type.
   *  
   *  DiffOp allows the use of \c operator \c - in FieldExpression objects for subtraction.
   *  That is, for any legal FieldExpression objects, \c operand1 and \c operand2, the following is a legal FieldExpression:
   *  \code
   *  operand1 - operand2
   *  \endcode
   *  When \c operator \c <<= is used to compute/assign the above DiffOp, the successive elements of \c operand1 and \c operand2 will be computed, subtracted, and returned/assigned without a temporary FieldType object being used as an intermediary in the computation.
   *  
   *  \warning DiffOp, along with the other SpecificBinOp operators, uses the standard C++ operator precedence.
   *  Thus, the two following FieldExpression objects return different results (assuming \c a is a FieldType object and using ProdOp):
   *  \code
   *  a - 4 * 6
   *  \endcode
   *  and
   *  \code
   *  (a - 4) * 6
   *  \endcode
   *  
   *  
   *  
   *  @typedef FieldType::value_type typedef DiffOp::AtomicType
   *  @brief Typedef for AtomicType.
   *  
   *  Typedef for AtomicType based on FieldType's value_type.
   *  
   *  
   *  
   *  @fn DiffOp::DiffOp(Operand1 op1, Operand2 op2)
   *  @brief Constructor for DiffOp.
   *  
   *  \param op1 an Operand1 object.
   *  \param op2 an Operand2 object.
   *  \return a new DiffOp containing op1 and op2.
   *  
   *  Construct a new DiffOp containing op1 and op2.
   *  
   *  
   *  
   *  @fn DiffOp::first() const
   *  @brief Operand1's accessor
   *  
   *  \return constant reference to op1 (first operand passed to constructor).
   *  
   *  Returns constant reference to first operand.
   *  
   *  
   *  
   *  @fn DiffOp::second() const
   *  @brief Operand2's accessor
   *  
   *  \return constant reference to op2 (second operand passed to constructor).
   *  
   *  Returns constant reference to second operand.
   *  
   *  
   *  
   *  @fn void DiffOp::init()
   *  @brief Initializer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Initializes operands.
   *  
   *  
   *  
   *  @fn void DiffOp::next()
   *  @brief Incrementer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Increments operands.
   *  
   *  
   *  
   *  @fn AtomicType DiffOp::eval() const
   *  @brief Computes and returns current element.
   *  
   *  \return current element: difference of operands' current elements.
   *  
   *  Computes and returns difference of operands' current elements.
   *  
   *  
   *  
   *  @fn bool DiffOp::at_end() const
   *  @brief Predicate: Current position is the end?
   *  
   *  \return Boolean; true, if currently at end; false; if not.
   *  
   *  Returns whether or not either operands' current position is the end/last position.
   *  
   *  
   *  
   *  @fn bool DiffOp::has_length() const
   *  @brief Predicate: Can reach end position?
   *  
   *  \return Boolean; true, if one or both operands can reach the end.
   *  
   *  Returns whether or not either operands can reach the end position.
   */
  BUILD_BINARY_OPERATOR(DiffOp, -, operator -);

  /**
   *  @struct ProdOp
   *  @author Christopher Earl
   *  @date October, 2010
   *  @ingroup ExpressionTemplatesDetail
   *  
   *  @brief Specialized FieldExpression-style representation of the multiplication operator, \c operator \c *.
   *  
   *  \tparam Operand1 First operand's type.
   *  \tparam Operand2 Second operand's type.
   *  \tparam FieldType Field type.
   *  
   *  @par
   *  
   *  ProdOp allows the use of \c operator \c * in FieldExpression objects for multiplication.
   *  That is, for any legal FieldExpression objects, \c operand1 and \c operand2, the following is a legal FieldExpression:
   *  \code
   *  operand1 * operand2
   *  \endcode
   *  When \c operator \c <<= is used to compute/assign the above ProdOp, the successive elements of \c operand1 and \c operand2 will be computed, multiplied together, and returned/assigned without a temporary FieldType object being used as an intermediary in the computation.
   *  
   *  \warning ProdOp, along with the other SpecificBinOp operators, uses the standard C++ operator precedence.
   *  Thus, the two following FieldExpression objects return different results (assuming \c a is a FieldType object and using DiffOp):
   *  \code
   *  a - 4 * 6
   *  \endcode
   *  and
   *  \code
   *  (a - 4) * 6
   *  \endcode
   *  
   *  
   *  
   *  @typedef FieldType::value_type typedef ProdOp::AtomicType
   *  @brief Typedef for AtomicType.
   *  
   *  Typedef for AtomicType based on FieldType's value_type.
   *  
   *  
   *  
   *  @fn ProdOp::ProdOp(Operand1 op1, Operand2 op2)
   *  @brief Constructor for ProdOp.
   *  
   *  \param op1 an Operand1 object.
   *  \param op2 an Operand2 object.
   *  \return a new ProdOp containing op1 and op2.
   *  
   *  Construct a new ProdOp containing op1 and op2.
   *  
   *  
   *  
   *  @fn ProdOp::first() const
   *  @brief Operand1's accessor
   *  
   *  \return constant reference to op1 (first operand passed to constructor).
   *  
   *  Returns constant reference to first operand.
   *  
   *  
   *  
   *  @fn ProdOp::second() const
   *  @brief Operand2's accessor
   *  
   *  \return constant reference to op2 (second operand passed to constructor).
   *  
   *  Returns constant reference to second operand.
   *  
   *  
   *  
   *  @fn void ProdOp::init()
   *  @brief Initializer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Initializes operands.
   *  
   *  
   *  
   *  @fn void ProdOp::next()
   *  @brief Incrementer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Increments operands.
   *  
   *  
   *  
   *  @fn AtomicType ProdOp::eval() const
   *  @brief Computes and returns current element.
   *  
   *  \return current element: product of operands' current elements.
   *  
   *  Computes and returns product of operands' current elements.
   *  
   *  
   *  
   *  @fn bool ProdOp::at_end() const
   *  @brief Predicate: Current position is the end?
   *  
   *  \return Boolean; true, if currently at end; false; if not.
   *  
   *  Returns whether or not either operands' current position is the end/last position.
   *  
   *  
   *  
   *  @fn bool ProdOp::has_length() const
   *  @brief Predicate: Can reach end position?
   *  
   *  \return Boolean; true, if one or both operands can reach the end.
   *  
   *  Returns whether or not either operands can reach the end position.
   */
  BUILD_BINARY_OPERATOR(ProdOp, *, operator *);

  /**
   *  @struct DivOp
   *  @author Christopher Earl
   *  @date October, 2010
   *  @ingroup ExpressionTemplatesDetail
   *  
   *  @brief Specialized FieldExpression-style representation of the division operator, \c operator \c /.
   *  
   *  \tparam Operand1 First operand's type.
   *  \tparam Operand2 Second operand's type.
   *  \tparam FieldType Field type.
   *  
   *  @par
   *  
   *  DivOp allows the use of \c operator \c / in FieldExpression objects for division.
   *  That is, for any legal FieldExpression objects, \c operand1 and \c operand2, the following is a legal FieldExpression:
   *  \code
   *  operand1 / operand2
   *  \endcode
   *  When \c operator \c <<= is used to compute/assign the above DivOp, the successive elements of \c operand1 and \c operand2 will be computed, divided one from the other, and returned/assigned without a temporary FieldType object being used as an intermediary in the computation.
   *  
   *  \warning DivOp, along with the other SpecificBinOp operators, uses the standard C++ operator precedence.
   *  Thus, the two following FieldExpression objects return different results (assuming \c a is a FieldType object and using DiffOp):
   *  \code
   *  a - 4 / 6
   *  \endcode
   *  and
   *  \code
   *  (a - 4) / 6
   *  \endcode
   *  
   *  
   *  
   *  @typedef FieldType::value_type typedef DivOp::AtomicType
   *  @brief Typedef for AtomicType.
   *  
   *  Typedef for AtomicType based on FieldType's value_type.
   *  
   *  
   *  
   *  @fn DivOp::DivOp(Operand1 op1, Operand2 op2)
   *  @brief Constructor for DivOp.
   *  
   *  \param op1 an Operand1 object.
   *  \param op2 an Operand2 object.
   *  \return a new DivOp containing op1 and op2.
   *  
   *  Construct a new DivOp containing op1 and op2.
   *  
   *  
   *  
   *  @fn DivOp::first() const
   *  @brief Operand1's accessor
   *  
   *  \return constant reference to op1 (first operand passed to constructor).
   *  
   *  Returns constant reference to first operand.
   *  
   *  
   *  
   *  @fn DivOp::second() const
   *  @brief Operand2's accessor
   *  
   *  \return constant reference to op2 (second operand passed to constructor).
   *  
   *  Returns constant reference to second operand.
   *  
   *  
   *  
   *  @fn void DivOp::init()
   *  @brief Initializer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Initializes operands.
   *  
   *  
   *  
   *  @fn void DivOp::next()
   *  @brief Incrementer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Increments operands.
   *  
   *  
   *  
   *  @fn AtomicType DivOp::eval() const
   *  @brief Computes and returns current element.
   *  
   *  \return current element: division of operands' current elements.
   *  
   *  Computes and returns division of operands' current elements.
   *  
   *  
   *  
   *  @fn bool DivOp::at_end() const
   *  @brief Predicate: Current position is the end?
   *  
   *  \return Boolean; true, if currently at end; false; if not.
   *  
   *  Returns whether or not either operands' current position is the end/last position.
   *  
   *  
   *  
   *  @fn bool DivOp::has_length() const
   *  @brief Predicate: Can reach end position?
   *  
   *  \return Boolean; true, if one or both operands can reach the end.
   *  
   *  Returns whether or not either operands can reach the end position.
   */
  BUILD_BINARY_OPERATOR(DivOp, /, operator /);
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Addition via function notation (rather than using the addition operator).
   *  
   *  \param first a double.
   *  \param second a double.
   *  \return a double (addition of first and second).
   *  
   *  @par
   *  Addition via function notation (rather than using the addition operator).
   *  Provided for use with SumFcn.
   *  
   *  \warning Using this function hardcodes AtomicType (alternatively FieldType::value_type) to be double.
   */
  inline double add (double first, double second) {
    return first + second;
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Subtraction via function notation (rather than using the subtraction operator).
   *  
   *  \param first a double.
   *  \param second a double.
   *  \return a double (subtraction of first and second).
   *  
   *  @par
   *  Subtraction via function notation (rather than using the subtraction operator).
   *  Provided for use with DiffFcn.
   *  
   *  \warning Using this function hardcodes AtomicType (alternatively FieldType::value_type) to be double.
   */
  inline double subt (double first, double second) {
    return first - second;
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Multiplication via function notation (rather than using the multiplication operator).
   *  
   *  \param first a double.
   *  \param second a double.
   *  \return a double (multiplication of first and second).
   *  
   *  @par
   *  Multiplication via function notation (rather than using the multiplication operator).
   *  Provided for use with MultFcn.
   *  
   *  \warning Using this function hardcodes AtomicType (alternatively FieldType::value_type) to be double.
   */
  inline double mult (double first, double second) {
    return first * second;
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Division via function notation (rather than using the division operator).
   *  
   *  \param first a double.
   *  \param second a double.
   *  \return a double (division of first and second).
   *  
   *  @par
   *  Division via function notation (rather than using the division operator).
   *  Provided for use with DuvFcn.
   *  
   *  \warning Using this function hardcodes AtomicType (alternatively FieldType::value_type) to be double.
   */
  inline double div (double first, double second) {
    return first / second;
  };
  
  /**
   *  @struct SumFcn
   *  @author Christopher Earl
   *  @date October, 2010
   *  @ingroup ExpressionTemplatesDetail
   *  
   *  @brief Specialized FieldExpression-style representation of the add function.
   *  
   *  \tparam Operand1 First operand's type.
   *  \tparam Operand2 Second operand's type.
   *  \tparam FieldType Field type.
   *  
   *  SumFcn allows the use of the function \c add in FieldExpression objects.
   *  That is, for any legal FieldExpression objects, \c operand1 and \c operand2, the following is a legal FieldExpression:
   *  \code
   *  add(operand1, operand2)
   *  \endcode
   *  When operator <<= is used to compute/assign the above SumFcn, the successive elements of \c operand1 and \c operand2 will be computed, applied to add, and returned/assigned without the use of a temporary FieldType object being used as an intermediary in the computation.
   *  
   *  
   *   \note SumFcn duplicates the functionality of SumOp.  SumFcn was provided as an example of what BUILD_BINARY_FUNCTION can do.
   *  
   *  @typedef FieldType::value_type typedef SumFcn::AtomicType
   *  @brief Typedef for AtomicType.
   *  
   *  Typedef for AtomicType based on FieldType's value_type.
   *  
   *  
   *  
   *  @fn SumFcn::SumFcn(Operand1 op1, Operand2 op2)
   *  @brief Constructor for SumFcn.
   *  
   *  \param op1 an Operand1 object.
   *  \param op2 an Operand2 object.
   *  \return a new SumFcn containing op1 and op2.
   *  
   *  Construct a new SumFcn containing op1 and op2.
   *  
   *  
   *  
   *  @fn SumFcn::first() const
   *  @brief Operand1's accessor
   *  
   *  \return constant reference to op1 (first operand passed to constructor).
   *  
   *  Returns constant reference to first operand.
   *  
   *  
   *  
   *  @fn SumFcn::second() const
   *  @brief Operand2's accessor
   *  
   *  \return constant reference to op2 (second operand passed to constructor).
   *  
   *  Returns constant reference to second operand.
   *  
   *  
   *  
   *  @fn void SumFcn::init()
   *  @brief Initializer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Initializes operands.
   *  
   *  
   *  
   *  @fn void SumFcn::next()
   *  @brief Incrementer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Increments operands.
   *  
   *  
   *  
   *  @fn AtomicType SumFcn::eval() const
   *  @brief Computes and returns current element.
   *  
   *  \return current element: the result of applying operands' current elements to add.
   *  
   *  Computes and returns the result of applying operands' current elements to add.
   *  
   *  
   *  
   *  @fn bool sumFcn::at_end() const
   *  @brief Predicate: Current position is the end?
   *  
   *  \return Boolean; true, if currently at end; false; if not.
   *  
   *  Returns whether or not either operands' current position is the end/last position.
   *  
   *  
   *  
   *  @fn bool sumFcn::has_length() const
   *  @brief Predicate: Can reach end position?
   *  
   *  \return Boolean; true, if one or both operands can reach the end.
   *  
   *  Returns whether or not either operands can reach the end position.
   */
  BUILD_BINARY_FUNCTION(SumFcn, add, add);

  /**
   *  @struct DiffFcn
   *  @author Christopher Earl
   *  @date October, 2010
   *  @ingroup ExpressionTemplatesDetail
   *  
   *  @brief Specialized FieldExpression-style representation of the subt function.
   *  
   *  \tparam Operand1 First operand's type.
   *  \tparam Operand2 Second operand's type.
   *  \tparam FieldType Field type.
   *  
   *  DiffFcn allows the use of the function \c subt in FieldExpression objects.
   *  That is, for any legal FieldExpression objects, \c operand1 and \c operand2, the following is a legal FieldExpression:
   *  \code
   *  subt(operand1, operand2)
   *  \endcode
   *  When operator <<= is used to compute/assign the above DiffFcn, the successive elements of \c operand1 and \c operand2 will be computed, applied to subt, and returned/assigned without the use of a temporary FieldType object being used as an intermediary in the computation.
   *  
   *  
   *   \note DiffFcn duplicates the functionality of DiffOp.  DiffFcn was provided as an example of what BUILD_BINARY_FUNCTION can do.
   *  
   *  @typedef FieldType::value_type typedef DiffFcn::AtomicType
   *  @brief Typedef for AtomicType.
   *  
   *  Typedef for AtomicType based on FieldType's value_type.
   *  
   *  
   *  
   *  @fn DiffFcn::DiffFcn(Operand1 op1, Operand2 op2)
   *  @brief Constructor for DiffFcn.
   *  
   *  \param op1 an Operand1 object.
   *  \param op2 an Operand2 object.
   *  \return a new DiffFcn containing op1 and op2.
   *  
   *  Construct a new DiffFcn containing op1 and op2.
   *  
   *  
   *  
   *  @fn DiffFcn::first() const
   *  @brief Operand1's accessor
   *  
   *  \return constant reference to op1 (first operand passed to constructor).
   *  
   *  Returns constant reference to first operand.
   *  
   *  
   *  
   *  @fn DiffFcn::second() const
   *  @brief Operand2's accessor
   *  
   *  \return constant reference to op2 (second operand passed to constructor).
   *  
   *  Returns constant reference to second operand.
   *  
   *  
   *  
   *  @fn void DiffFcn::init()
   *  @brief Initializer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Initializes operands.
   *  
   *  
   *  
   *  @fn void DiffFcn::next()
   *  @brief Incrementer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Increments operands.
   *  
   *  
   *  
   *  @fn AtomicType DiffFcn::eval() const
   *  @brief Computes and returns current element.
   *  
   *  \return current element: the result of applying operands' current elements to subt.
   *  
   *  Computes and returns the result of applying operands' current elements to subt.
   *  
   *  
   *  
   *  @fn bool DiffFcn::at_end() const
   *  @brief Predicate: Current position is the end?
   *  
   *  \return Boolean; true, if currently at end; false; if not.
   *  
   *  Returns whether or not either operands' current position is the end/last position.
   *  
   *  
   *  
   *  @fn bool DiffFcn::has_length() const
   *  @brief Predicate: Can reach end position?
   *  
   *  \return Boolean; true, if one or both operands can reach the end.
   *  
   *  Returns whether or not either operands can reach the end position.
   */
  BUILD_BINARY_FUNCTION(DiffFcn, subt, subt);

  /**
   *  @struct MultFcn
   *  @author Christopher Earl
   *  @date October, 2010
   *  @ingroup ExpressionTemplatesDetail
   *  
   *  @brief Specialized FieldExpression-style representation of the mult function.
   *  
   *  \tparam Operand1 First operand's type.
   *  \tparam Operand2 Second operand's type.
   *  \tparam FieldType Field type.
   *  
   *  MultFcn allows the use of the function \c mult in FieldExpression objects.
   *  That is, for any legal FieldExpression objects, \c operand1 and \c operand2, the following is a legal FieldExpression:
   *  \code
   *  mult(operand1, operand2)
   *  \endcode
   *  When operator <<= is used to compute/assign the above MultFcn, the successive elements of \c operand1 and \c operand2 will be computed, applied to mult, and returned/assigned without the use of a temporary FieldType object being used as an intermediary in the computation.
   *  
   *  
   *   \note MultFcn duplicates the functionality of MultOp.  MultFcn was provided as an example of what BUILD_BINARY_FUNCTION can do.
   *  
   *  @typedef FieldType::value_type typedef MultFcn::AtomicType
   *  @brief Typedef for AtomicType.
   *  
   *  Typedef for AtomicType based on FieldType's value_type.
   *  
   *  
   *  
   *  @fn MultFcn::MultFcn(Operand1 op1, Operand2 op2)
   *  @brief Constructor for MultFcn.
   *  
   *  \param op1 an Operand1 object.
   *  \param op2 an Operand2 object.
   *  \return a new MultFcn containing op1 and op2.
   *  
   *  Construct a new MultFcn containing op1 and op2.
   *  
   *  
   *  
   *  @fn MultFcn::first() const
   *  @brief Operand1's accessor
   *  
   *  \return constant reference to op1 (first operand passed to constructor).
   *  
   *  Returns constant reference to first operand.
   *  
   *  
   *  
   *  @fn MultFcn::second() const
   *  @brief Operand2's accessor
   *  
   *  \return constant reference to op2 (second operand passed to constructor).
   *  
   *  Returns constant reference to second operand.
   *  
   *  
   *  
   *  @fn void MultFcn::init()
   *  @brief Initializer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Initializes operands.
   *  
   *  
   *  
   *  @fn void MultFcn::next()
   *  @brief Incrementer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Increments operands.
   *  
   *  
   *  
   *  @fn AtomicType MultFcn::eval() const
   *  @brief Computes and returns current element.
   *  
   *  \return current element: the result of applying operands' current elements to mult.
   *  
   *  Computes and returns the result of applying operands' current elements to mult.
   *  
   *  
   *  
   *  @fn bool MultFcn::at_end() const
   *  @brief Predicate: Current position is the end?
   *  
   *  \return Boolean; true, if currently at end; false; if not.
   *  
   *  Returns whether or not either operands' current position is the end/last position.
   *  
   *  
   *  
   *  @fn bool MultFcn::has_length() const
   *  @brief Predicate: Can reach end position?
   *  
   *  \return Boolean; true, if one or both operands can reach the end.
   *  
   *  Returns whether or not either operands can reach the end position.
   */
  BUILD_BINARY_FUNCTION(MultFcn, mult, mult);

  /**
   *  @struct DivFcn
   *  @author Christopher Earl
   *  @date October, 2010
   *  @ingroup ExpressionTemplatesDetail
   *  
   *  @brief Specialized FieldExpression-style representation of the div function.
   *  
   *  \tparam Operand1 First operand's type.
   *  \tparam Operand2 Second operand's type.
   *  \tparam FieldType Field type.
   *  
   *  DivFcn allows the use of the function \c div in FieldExpression objects.
   *  That is, for any legal FieldExpression objects, \c operand1 and \c operand2, the following is a legal FieldExpression:
   *  \code
   *  div(operand1, operand2)
   *  \endcode
   *  When operator <<= is used to compute/assign the above DivFcn, the successive elements of \c operand1 and \c operand2 will be computed, applied to div, and returned/assigned without the use of a temporary FieldType object being used as an intermediary in the computation.
   *  
   *  
   *   \note DivFcn duplicates the functionality of DivOp.  DivFcn was provided as an example of what BUILD_BINARY_FUNCTION can do.
   *  
   *  @typedef FieldType::value_type typedef DivFcn::AtomicType
   *  @brief Typedef for AtomicType.
   *  
   *  Typedef for AtomicType based on FieldType's value_type.
   *  
   *  
   *  
   *  @fn DivFcn::DivFcn(Operand1 op1, Operand2 op2)
   *  @brief Constructor for DivFcn.
   *  
   *  \param op1 an Operand1 object.
   *  \param op2 an Operand2 object.
   *  \return a new DivFcn containing op1 and op2.
   *  
   *  Construct a new DivFcn containing op1 and op2.
   *  
   *  
   *  
   *  @fn DivFcn::first() const
   *  @brief Operand1's accessor
   *  
   *  \return constant reference to op1 (first operand passed to constructor).
   *  
   *  Returns constant reference to first operand.
   *  
   *  
   *  
   *  @fn DivFcn::second() const
   *  @brief Operand2's accessor
   *  
   *  \return constant reference to op2 (second operand passed to constructor).
   *  
   *  Returns constant reference to second operand.
   *  
   *  
   *  
   *  @fn void DivFcn::init()
   *  @brief Initializer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Initializes operands.
   *  
   *  
   *  
   *  @fn void DivFcn::next()
   *  @brief Incrementer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Increments operands.
   *  
   *  
   *  
   *  @fn AtomicType DivFcn::eval() const
   *  @brief Computes and returns current element.
   *  
   *  \return current element: the result of applying operands' current elements to div.
   *  
   *  Computes and returns the result of applying operands' current elements to div.
   *  
   *  
   *  
   *  @fn bool DivFcn::at_end() const
   *  @brief Predicate: Current position is the end?
   *  
   *  \return Boolean; true, if currently at end; false; if not.
   *  
   *  Returns whether or not either operands' current position is the end/last position.
   *  
   *  
   *  
   *  @fn bool DivFcn::has_length() const
   *  @brief Predicate: Can reach end position?
   *  
   *  \return Boolean; true, if one or both operands can reach the end.
   *  
   *  Returns whether or not either operands can reach the end position.
   */
  BUILD_BINARY_FUNCTION(DivFcn, div, div);
  
  /**
   *  @struct SinFcn
   *  @author Christopher Earl
   *  @date October, 2010
   *  @ingroup ExpressionTemplatesDetail
   *  
   *  @brief Specialized FieldExpression-style representation of the sin function.
   *  
   *  \tparam Operand First operand's type.
   *  \tparam FieldType Field type.
   *  
   *  SinFcn allows the use of the function \c sin in FieldExpression objects.
   *  That is, for any legal FieldExpression object, \c operand, the following is a legal FieldExpression:
   *  \code
   *  sin(operand)
   *  \endcode
   *  When operator <<= is used to compute/assign the above SinFcn, the successive elements of \c operand will be computed, applied to sin, and returned/assigned without the use of a temporary FieldType object being used as an intermediary in the computation.
   *  
   *  
   *  
   *  @typedef FieldType::value_type typedef SinFcn::AtomicType
   *  @brief Typedef for AtomicType.
   *  
   *  Typedef for AtomicType based on FieldType's value_type.
   *  
   *  
   *  
   *  @fn SinFcn::SinFcn(Operand oper)
   *  @brief Constructor for SinFcn.
   *  
   *  \param oper an Operand object.
   *  \return a new SinFcn containing oper.
   *  
   *  Construct a new SinFcn containing oper.
   *  
   *  
   *  
   *  @fn SinFcn::first() const
   *  @brief Operand1's accessor
   *  
   *  \return constant reference to op1 (first operand passed to constructor).
   *  
   *  Returns constant reference to first operand.
   *  
   *  
   *  
   *  @fn SinFcn::second() const
   *  @brief Operand2's accessor
   *  
   *  \return constant reference to op2 (second operand passed to constructor).
   *  
   *  Returns constant reference to second operand.
   *  
   *  
   *  
   *  @fn void SinFcn::init()
   *  @brief Initializer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Initializes operand.
   *  
   *  
   *  
   *  @fn void SinFcn::next()
   *  @brief Incrementer.
   *  
   *  \return nothing. Called for side-effects.
   *  
   *  Increments operand.
   *  
   *  
   *  
   *  @fn AtomicType SinFcn::eval() const
   *  @brief Computes and returns current element.
   *  
   *  \return current element: the result of applying operand's current element to sin.
   *  
   *  Computes and returns the result of applying operand's current element to sin.
   *  
   *  
   *  
   *  @fn bool SinFcn::at_end() const
   *  @brief Predicate: Current position is the end?
   *  
   *  \return Boolean; true, if currently at end; false; if not.
   *  
   *  Returns whether or not operand's current position is the end/last position.
   *  
   *  
   *  
   *  @fn bool SinFcn::has_length() const
   *  @brief Predicate: Can reach end position?
   *  
   *  \return Boolean; true, if operand can reach the end.
   *  
   *  Returns whether or not operand can reach the end position.
   */
  BUILD_UNARY_FUNCTION(SinFcn, std::sin, sin);
  BUILD_UNARY_FUNCTION(CosFcn, std::cos, cos);
  BUILD_UNARY_FUNCTION(TanFcn, std::tan, tan);

  BUILD_UNARY_FUNCTION(ExpFcn, std::exp, exp);

  BUILD_UNARY_FUNCTION(TanhFcn, std::tanh, tanh);

  BUILD_UNARY_FUNCTION(AbsFcn, std::abs, abs);

} // namespace SpatialOps

#endif // SpatialOps_FieldOperationDefinitions_h
