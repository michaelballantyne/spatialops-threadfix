#ifndef SpatialOps_FieldOperations_h
#define SpatialOps_FieldOperations_h

namespace SpatialOps{
  /**
   *  @file FieldOperations.h
   */


  /**
   *  @page expressiontemplates
   *  @ingroup ExpressionTemplates
   *  @par Introduction
   *   The FieldExpression structure and \c operator \c <<= allow assignment-style statements that compute a given expression involving fields and assign the results to another given field.
   *   \c operator \c <<= uses the given FieldExpression object to construct a loop using forward iterators to evaluate the given expression for each element of the fields, one at a time.
   *   The field_type itself is an inferred parameter to FieldExpression.
   *   (The requirements for the field type are given below.)
   *   
   *  @par Example
   *   \code
   *    extended_vector<int> a, b, c;
   *    ...
   *    c <<= a + sin(b);
   *   \endcode
   *   The final statement in this example expands to a loop which evaluates the expression the equivalent of the following:
   *   \code
   *    extended_vector::iterator c_iter = c.begin();
   *    extended_vector::const_iterator a_iter = a.begin();
   *    extended_vector::const_iterator b_iter = b.begin();
   *    while(c_iter != c.end()) {
   *       *c_iter = a_iter + sin(b_iter);
   *       ++c_iter;
   *       ++a_iter;
   *       ++b_iter;
   *    };
   *   \endcode
   *  
   *  @par \c operator \c <<= Syntax
   *   The syntax for a \c operator \c <<= assignment statement is as follows:
   *   \code
   *    field <<= fieldexpr;
   *   \endcode
   *   where \c field is a mutable field (a field_type object) and \c fieldexpr is a FieldExpression object.
   *  
   *  @par FieldExpression Syntax
   *   The syntax for building a FieldExpression object was to be straightforward yet extensible.
   *   
   *   \c operator \c <<= is overloaded so that scalars and fields can be used as FieldExpressions.  For example, if \c a and \c b are fields, both
   *   \code
   *    b <<= 42;
   *   \endcode
   *   and
   *   \code
   *    b <<= a;
   *   \endcode
   *   will compile.
   *   (The first will assign the constant 42 to every element of \c b, and the second will copy the elements of \c a into the elements of \c b.)
   *   
   *   A true FieldExpression object is an expression and therefore represents a function/operator application, whose operands may themselves be expressions.
   *  
   *  @par Field Type Requirements
   */

  template<typename AtomicType>
    struct Scalar;

  template<typename FieldType>
    struct FieldForm;

  template <typename Operand1, typename Operand2, typename FieldType>
    struct BinOp;
  
#define BUILD_BINARY_TYPE_PROTOTYPE(OBJECT_NAME)			\
  template <typename Operand1, typename Operand2, typename FieldType>	\
    struct OBJECT_NAME

  template <typename Operand, typename FieldType>
    struct UnOp;

#define BUILD_UNARY_TYPE_PROTOTYPE(OBJECT_NAME)		\
  template <typename Operand, typename FieldType>	\
    struct OBJECT_NAME

  template<int Num, typename FieldType>
    struct ArgForm;
  
  template<typename ExprType, typename FieldType>
    struct FieldExpression;

  template<typename ExprType, int CrtNum, typename Max, typename FieldType>
    struct FcnForm;
  
  template<int Num>
    struct Int;

  template<typename Max1, typename Max2>
    struct MetaFcnMax;

  template<typename Input, typename FieldType>
    struct StandardizeTerm;

  template<typename NewType, typename OldTerm, typename FieldType>
    struct LiftTerm;

  template<typename NewType, typename OldTerm1, typename OldTerm2, typename FieldType>
    struct CombineTerms;

  template<typename Arg, typename FieldType>
    struct StandardizeArg;

  template<typename BeginType, int CurrentArg, typename ArgType>
    struct ArgReplace;

  template<typename StandardType, int CrtNum, typename Max, typename FieldType>
    struct ArgResultTerm;
  
  /**
   *  @struct Scalar
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief FieldExpression-style representation of an element (AtomicType).
   *  
   *  \tparam AtomicType Basic element type.
   *  
   *  Scalar is a simple container structure that generalizes the use of AtomicType objects for the FieldExpression framework.
   *
   *  \note
   *  Current implementation keeps a copy of element contained rather than a pointer to the original.
   *  This was done with the assumption that AtomicType is the same size as a pointer.
   *  (Thus keeping a pointer does not reduce the size of Scalar and does increase the look-up time.)
   */
  template<typename AtomicType>
    struct Scalar {
    private:
      AtomicType val;
      
    public:
      
      /**
       *  @brief Constructor for Scalar.
       *  
       *  \param v an AtomicType (template parameter to Scalar).
       *  \return a new Scalar with internal value set to v.
       *  
       *  Constructs a new Scalar with internal value set to the parameter v.
       */
    Scalar(AtomicType const & v)
    : val(v)
      {};
      
      /**
       *  @brief Initializer.
       *  
       *  \return nothing. Called for side-effects (none here).
       *  
       *  A no-op. Included for regularity of complex FieldExpression objects.
       *  Since this function is inlined, it disappears in the executed code.
       */
      inline void init() const {};
      
      /**
       *  @brief Incrementer.
       *  
       *  \return nothing. Called for side-effects (none here).
       *  
       *  A no-op. Included for regularity of complex FieldExpression objects.
       *  Since this function is inlined, it disappears in the executed code.
       */
      inline void next() const {};
      
      /**
       *  @brief Reference to current value.
       *  
       *  \return a constant reference to the current element (an AtomicType/FieldType::value_type object).
       *  
       *  Returns a constant reference to the current element, which is constant.
       */
      inline AtomicType const & eval() const {
	return val;
      };
    };

  /**
   *  @struct FieldForm
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief FieldExpression-style representation of a field (FieldType).
   *  
   *  \tparam FieldType Field type.
   *  
   *  FieldType is a simple container structure that generalizes the use of FieldType objects for the FieldExpression framework.
   *
   *  A FieldType object must provide the following typedefs:
   *  
   *   \li \c field_type A FieldType object's own type (i.e. FieldType). For example, if \c SpatialField is a FieldType, 
   *       then SpatialField must contain:
   *       \code
   *   typedef SpatialField field_type
   *       \endcode
   *  
   *   \li \c value_type A FieldType object's element/atomic type.  Similar to \c field_type, except
   *       \c value_type is the type of the elements contained within FieldType.  Also known as
   *       AtomicType elsewhere in the code and this documentation.
   *  
   *   \li \c iterator Commonly a typedef for: \c value_type*.
   *  
   *   \li \c const_iterator Commonly a typedef for: \c value_type \c const*.
   *  
   *  A FieldType object must provide the following methods:
   *  
   *   \li \code iterator begin() \endcode
   *       Returns a pointer to the first element in the current FieldType object.
   *  
   *   \li \code iterator end() \endcode
   *       Returns a pointer to the final element in the current FieldType object.
   *  
   *   \li \code iterator operator ++ (iterator &) \endcode
   *       and
   *       \code const_iterator operator ++ (const_iterator &) \endcode
   *       Increments argument iterator/constant iterator.
   *   \note 
   *       Return value is not used; only the side-effect is important.
   */
  template<typename FieldType>
    struct FieldForm {
    private:
      FieldType const *fptr;
      typename FieldType::const_iterator iter;
      
    public:
      
      /**
       *  @brief Constructor for FieldForm.
       *  
       *  \param field a FieldType (template parameter to FieldForm).
       *  \return a new FieldForm containing field.
       *  
       *  Constructs a new Fieldform containing field.
       *
       *  Current implementation initializes field (may be redundant with init() function).
       */
    FieldForm(FieldType const & field)
    : fptr(&field), iter(field.begin())
      {};
      
      /**
       *  @brief Initializer.
       *  
       *  \return nothing. Called for side-effects.
       *  
       *  Initializes internal iterator to beginning the contained field.
       */
      inline void init() {
	iter = fptr->begin();
      };
      
      /**
       *  @brief Incrementer.
       *  
       *  \return nothing. Called for side-effects.
       *  
       *  Increments internal iterator to next element in the contained field.
       */
      inline void next() {
	++iter;
      };
      
      /**
       *  @brief Reference to current element.
       *  
       *  \return a constant reference to the current element (an AtomicType/FieldType::value_type object) in the field.
       *  
       *  Returns a constant reference to the current element in the contained field.
       */
      inline typename FieldType::value_type const & eval() const {
	return *iter;
      };
    };

  /**
   *  @struct BinOp/SpecificBinOp
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief General/specialized FieldExpression-style representation of a binary operation/function.
   *  
   *  \tparam Operand1 First operand's type.
   *  \tparam Operand2 Second operand's type.
   *  \tparam FieldType Field type.
   *  
   *  @par BinOp
   *   BinOp is the general, non-optimized representation of a binary function and therefore requires the function be passed to it.
   *   Briefly, to use a function, \c fcn, that does not have a SpecificBinOp defined for it, with operands \c op1 and \c op2, the \c app function is called:
   *   \code
   *   app(fcn, op1, op2)
   *   \endcode
   *   The signature for \c app is:
   *   \code
   *   BinOp<Operand1, Operand2, FieldType> app (typename FieldType::value_type (*)(typename FieldType::value_type, typename FieldType::value_type),
   *                                             Operand1 const &,
   *                                             Operand2 const &)
   *   \endcode
   *  
   *  @par SpecificBinOp
   *   Commonly used binary operators and functions have been given individual representations, for both optimization and ease-of-use reasons.
   *   (These optimized BinOp-like structures are generated with macros by the preprocessor and so do not show up directly in this documentation.)
   *   To use these sturctures, the name given to the macro is used - usually identical to the operator or function itself.
   *   For example, for addition of two operands, \c op1 and \c op2, use the '+' symbol with infix notation:
   *   \code
   *   op1 + op2
   *   \endcode
   *   The macros, \c BUILD_BINARY_OPERATOR and \c BUILD_BINARY_FUNCTION, define binary operators and functions respectively.
   *   For each member function in BinOp, there is an equivalent member function in each SpecificBinOp (without the function parameter).
   *   Note that usual C/C++ order of operations applies to binary operators defined in this way.
   */
  template <typename Operand1, typename Operand2, typename FieldType>
    struct BinOp {
      
      /**
       *  @brief Typedef for AtomicType.
       *  
       *  Typedef for AtomicType based on FieldType's value_type.
       */
      typename FieldType::value_type typedef AtomicType;
      
    private:
      typedef AtomicType (*OperationType)(AtomicType, AtomicType);
      OperationType op;
      Operand1 operand1;
      Operand2 operand2;
      
    public:
      
      /**
       *  @brief Constructor for BinOp.
       *  
       *  \param operation a function that takes two AtomicType objects and returns an AtomicType object.
       *  \param op1 an Operand1 object.
       *  \param op2 an Operand2 object.
       *  \return a new BinOp containing operation, op1, and op2.
       *  
       *  Constructs a new BinOp containing operation, op1, and op2.
       */
    BinOp(OperationType operation,//AtomicType (*operation)(AtomicType, AtomicType),
	  Operand1 op1,
	  Operand2 op2)
    : op(operation), operand1(op1), operand2(op2)
      {};
      
      /**
       *  @brief Operand1's accessor
       *  
       *  \return constant reference to op1 (first operand passed to constructor).
       *  
       *  Returns constant reference to first operand.
       */
      inline Operand1 const & first () const {
	return operand1;
      };
      
      /**
       *  @brief Operand2's accessor
       *  
       *  \return constant reference to op2 (second operand passed to constructor).
       *  
       *  Returns constant reference to second operand.
       */
      inline Operand2 const & second () const {
	return operand2;
      };
      
      /**
       *  @brief Operator's accessor
       *  
       *  \return constant reference to internal function (function passed to constructor).
       *  
       *  Returns constant reference to function.
       */
      inline OperationType const & fcn () const {
	return op;
      };
      
      /**
       *  @brief Initializer.
       *  
       *  \return nothing. Called for side-effects.
       *  
       *  Initializes operands.
       */
      inline void init() {
	operand1.init();
	operand2.init();
      };
      
      /**
       *  @brief Incrementer.
       *  
       *  \return nothing. Called for side-effects.
       *  
       *  Increments operands.
       */
      inline void next() {
	operand1.next();
	operand2.next();
      };
      
      /**
       *  @brief Computes and returns current element.
       *  
       *  \return current element based on current operands and operation.
       *  
       *  Computes and returns current element base on operands' current elements and the operation function.
       */
      inline AtomicType eval() const {
	return op(operand1.eval(),
		  operand2.eval());
      };
    };
  
#define BUILD_BINARY_OPERATOR_STRUCT(OBJECT_NAME, INTERNAL_NAME)	\
  template <typename Operand1, typename Operand2, typename FieldType>	\
    struct OBJECT_NAME {						\
									\
    typename FieldType::value_type typedef AtomicType;			\
    									\
    private:								\
    Operand1 operand1;							\
    Operand2 operand2;							\
    									\
    public:								\
    									\
    OBJECT_NAME(Operand1 op1,						\
		Operand2 op2)						\
    : operand1(op1), operand2(op2)					\
      {};								\
									\
    inline Operand1 const & first () const {				\
      return operand1;							\
    };									\
    									\
    inline Operand2 const & second () const {				\
      return operand2;							\
    };									\
    									\
    inline void init() {						\
      operand1.init();							\
      operand2.init();							\
    };									\
    									\
    inline void next() {						\
      operand1.next();							\
      operand2.next();							\
    };									\
    									\
    inline AtomicType eval() const {					\
      return operand1.eval() INTERNAL_NAME				\
	operand2.eval();						\
    };									\
    }
  
#define BUILD_BINARY_FUNCTION_STRUCT(OBJECT_NAME, INTERNAL_NAME)	\
  template <typename Operand1, typename Operand2, typename FieldType>	\
    struct OBJECT_NAME {						\
    									\
    typename FieldType::value_type typedef AtomicType;			\
    									\
    private:								\
    Operand1 operand1;							\
    Operand2 operand2;							\
    									\
    public:								\
    									\
    OBJECT_NAME(Operand1 op1,						\
		Operand2 op2)						\
    : operand1(op1), operand2(op2)					\
      {};								\
    									\
    inline Operand1 const & first () const {				\
      return operand1;							\
    };									\
    									\
    inline Operand2 const & second () const {				\
      return operand2;							\
    };									\
    									\
    inline void init() {						\
      operand1.init();							\
      operand2.init();							\
    };									\
    									\
    inline void next() {						\
      operand1.next();							\
      operand2.next();							\
    };									\
    									\
    inline AtomicType eval() const {					\
      return INTERNAL_NAME(operand1.eval(),				\
			   operand2.eval());				\
    };									\
    }
  
  /**
   *  @struct UnOp/SpecificUnOp
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief General/specialized FieldExpression-style representation of a unary function.
   *  
   *  \tparam Operand Operand's type.
   *  \tparam FieldType Field type.
   *  
   *  @par UnOp
   *   UnOp is the non-optimized representation of a unary function and therefore requires the function be passed to it.
   *  Briefly, to use a function, \c fcn, that does not have a SpecificUnOp defined for it, with operand \c op, the \c app function is called:
   *  \code
   *  app(fcn, op)
   *  \endcode
   *  The signature for \c app is:
   *  \code
   *  UnOp<Operand, FieldType> app (typename FieldType::value_type (*)(typename FieldType::value_type),
   *                                Operand const &)
   *  \endcode
   *  
   *  @par SpecificUnOp
   *  Commonly used unary functions have been given individual representations, for both optimization and ease-of-use reasons.
   *  (These optimized UnOp-like structures are generated with macros by the preprocessor and so do not show up directly in this documentation.)
   *  To use these sturctures, the name given to the macro is used - usually identical to the function itself.
   *  For example, for sin of the operand, \c op, usage is identical to applying sin to a double:
   *  \code
   *  sin(op)
   *  \endcode
   *  The macro, \c BUILD_UNARY_FUNCTION, defines unary functions.
   *   For each member function in UnOp, there is an equivalent member function in each SpecificUnOp (without the function parameter).
   *  Note that usual C/C++ unary operators can be defined by macros very similar to \c BUILD_UNARY_FUNCTION; however, this macro has not been implemented, because there is no immediate use.
   */
  template <typename Operand, typename FieldType>
    struct UnOp {
      
      /**
       *  @brief Typedef for AtomicType.
       *  
       *  Typedef for AtomicType based on FieldType's value_type.
       */
      typename FieldType::value_type typedef AtomicType;
      
    private:
      typedef AtomicType (*OperationType)(AtomicType, AtomicType);
      OperationType op;
      Operand operand;
      
    public:
      
      /**
       *  @brief Constructor for UnOp.
       *  
       *  \param operation a function that takes an AtomicType object and returns an AtomicType object.
       *  \param oper an Operand object.
       *  \return a new UnOp containing operation and oper.
       *  
       *  Constructs a new UnOp containing operation and oper.
       */
    UnOp(AtomicType (*operation)(AtomicType, AtomicType),
	 Operand oper)
    : op(operation), operand(oper)
      {};
      
      /**
       *  @brief Operand's accessor
       *  
       *  \return constant reference to oper (operand passed to constructor).
       *  
       *  Returns constant reference to operand.
       */
      inline Operand const & oper () const {
	return operand;
      };
      
      /**
       *  @brief Operator's accessor
       *  
       *  \return constant reference to internal function (function passed to constructor).
       *  
       *  Returns constant reference to function.
       */
      inline OperationType const & fcn () const {
	return op;
      };
      
      /**
       *  @brief Initializer.
       *  
       *  \return nothing. Called for side-effects.
       *  
       *  Initializes operand.
       */
      inline void init() {
	operand.init();
      };
      
      /**
       *  @brief Incrementer.
       *  
       *  \return nothing. Called for side-effects.
       *  
       *  Increments operand.
       */
      inline void next() {
	operand.next();
      };
      
      /**
       *  @brief Computes and returns current element.
       *  
       *  \return current element based on current operand and operation.
       *  
       *  Computes and returns current element base on the operand's current element and the operation function.
       */
      inline AtomicType eval() const {
	return op(operand.eval());
      };
    };

#define BUILD_UNARY_STRUCT(OBJECT_NAME, INTERNAL_NAME)		\
  template <typename Operand, typename FieldType>		\
    struct OBJECT_NAME {					\
    								\
    typename FieldType::value_type typedef AtomicType;		\
    								\
    private:							\
    Operand operand;						\
								\
    public:							\
								\
    OBJECT_NAME(Operand oper)					\
    : operand(oper)						\
      {};							\
    								\
    inline Operand const & oper() const {			\
      return operand;						\
    };								\
								\
    inline void init() {					\
      operand.init();						\
    };								\
								\
    inline void next() {					\
      operand.next();						\
    };								\
    								\
    inline AtomicType eval() const {				\
      return INTERNAL_NAME(operand.eval());			\
    };								\
    }
  
  /**
   *  @struct ArgForm
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief FieldExpression-style representation of an anonymous argument.
   *  
   *  \tparam Num Argument number.
   *  \tparam FieldType Field type.
   *  
   *   Anonymous arguments are used to create anonymous functions (FcnForm).
   *   (Recursive definition: An anonymous function is a FieldExpression containing at least one anonymous argument.)
   *   When expressions are applied to an anonymous function, each in order is associated with an integer, beginning at 0, for the argument it replaces.
   *   So in the following FieldExpression,
   *   \code
   *   (4 + ($1 - $0))(a, b)
   *   \endcode
   *   The expression \c a is associated with 0, and the expression \c b with 1.
   *   When applying an argument, everywhere the ArgForm with the associated number appears, it is replaced by the argument.
   *   Reusing the above example,
   *   \code
   *   (4 + ($1 - $0))(a, b)
   *   \endcode
   *   becomes
   *   \code
   *   (4 + ($1 - a))(b)
   *   \endcode
   *   which itself becomes
   *   \code
   *   (4 + (b - a))
   *   \endcode
   *   For more information on ArgForm, see FcnForm.
   *   
   *  @todo Add application functionality to simple/bare ArgForm.
   *   With this it will be possible to apply parameters to ArgForm directly.
   *   For example,
   *   \code
   *   ArgForm<0,FieldType>(4)
   *   \endcode
   *   will cause a compilation error rather than returning a Scalar, with a value of 4.
   */
  template<int Num, typename FieldType>
    struct ArgForm {
      FieldType typedef field_type;
  
      ArgForm()
      {};
    };
  
  /**
   *  @struct FieldExpression
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief FieldExpression-style representation of an expression to evaluate.
   *  
   *  \tparam ExprType The type representing the actual expression represented.
   *  \tparam FieldType Field type.
   *  
   *  A FieldExpression represents an expression to evaluate.
   *  A FieldExpression is either:
   *   \li A Scalar (based on an AtomicType, \c typename \c FieldType::value_type),
   *   \li A FieldForm (based on a FieldType),
   *   \li A BinOp of two expressions,
   *   \li A SpecificBinOp (see BinOp) of two expressions,
   *   \li A UnOp of an expression,
   *   \li A SpecificUnOp (see UnOp) of an expression, or
   *   \li The result of a FcnForm after all its arguments have been applied (see FcnForm for more information).
   *  
   *  @par
   *  
   *  The FieldExpression structure itself is a container structure and does not contain itself, even though two FieldExpressions may be used as operands to build a new FieldExpression.
   *  Equivalently, the FieldExpression structure only appears at the top of FieldExpression and never appears in a subexpression.
   *  For an example, consider the following situation:
   *   \li \c FT is the current FieldType,
   *   \li \c AT is the current AtomicType (\c typename \c FieldType::value_type),
   *   \li \c a has type \c FT,
   *   \li \c 4 has type \c AT,
   *   \li \c operator \c + has the return type \c SumOp defined for it using the macro \c BUILD_BINARY_OPERATOR (see BinOp for more information), and
   *   \li the expression
   *  \code
   *  a + 4
   *  \endcode
   *   appears in the program, and so needs to be evaluated.
   *  
   *  To determine the type of this expression, the compiler goes through a recursive process.
   *  First, it finds the type of the first subexpression \c a.
   *  This examination finds the following return type:
   *  \code
   *  FieldExpression<FieldForm<FT>, FT>
   *  \endcode
   *  Next, the compiler finds the type of the second subexpression \c 4.
   *  This examination finds the following return type:
   *  \code
   *  FieldExpression<Scalar<AT>, FT>
   *  \endcode
   *  Finally, the two subexpressions are stripped of their FieldExpression containers and used in a \c SumOp type; so the expression
   *  \code
   *  a + 4
   *  \endcode
   *  returns the type:
   *  \code
   *  FieldExpression<SumOp<FieldForm<FT>,
   *                   Scalar<AT>,
   *                   FT>,
   *             FT>
   *  \endcode
   *
   *  \note There is currently no (direct) way to have a function that takes three or more arguments in a FieldExpression.
   *  This can be added very easily (the implementation will be extremely similar to that of BinOp).
   *  There is no current use, so it was not implemented.
   */
  template<typename Operand, typename FieldType>
    struct FieldExpression {
    public:
      FieldType typedef field_type;
      
    private:
      Operand expr;
  
    public:
    /**
     *  @brief Constructor for FieldExpression.
     *  
     *  \param given an ExprType that is the internal expression type.
     *  \return a new FieldExpression containing given.
     *  
     *  Constructs a new FieldExpression containing given.
     */
    FieldExpression(Operand const & given)
    : expr(given)
      {};
      
      /**
       *  @brief Accessor for internal expression.
       *
       *  \return constant reference to internal expression.
       *
       *  Returns constant reference to internal expression.
       */
      Operand const & expression() const {
	return expr;
      };
    };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief FieldExpression-style representation of an anonymous function.
   *  
   *  \tparam ExprType The type representing the function's body/actual expression.
   *  \tparam CrtNum Integer representing the next argument to be applied.
   *  \tparam MaxNum Integer representing the total number of arguments that need to be applied.
   *  \tparam FieldType Field type.
   *  
   *  An anonymous function is a FieldExpression that contains at least one anonymous argument (see ArgForm).
   *  A FcnForm takes the place of FieldExpression as the container structure used to standardize type interactions.
   *  A FcnForm also controls/regulates application of arguments.
   *  
   *  @par Application
   *   Application of arguments to an anonymous function simulates textual substitution.
   *   That is, when a FieldExpression is applied to a FcnForm, wherever the next argument appears in the ExprType it is replaced by the applied FieldExpression.
   *   The next argument is determined by the value in CrtNum, which is incremented after every application.
   *   (See the extended example below for a high-level demonstration.)
   *   Because of this (quasi-)textual substitution, the difference between running a fully-applied anonymous function and the expression it represents is very small.
   *   For a brief example, evaluating (at run-time) the anonymous function,
   *  \code
   *  ($1 + $0)(4, a)
   *  \endcode
   *   and the expression
   *  \code
   *  (a + 4)
   *  \endcode
   *   takes roughly the same amount of time, with only slightly higher overhead for the anonymous function (which is independent from the size or complexity of the FieldType).
   *  
   *  @par Reversion to FieldExpression
   *   Once all the arguments that appear in an anonymous function are replaced by expressions via application, the container structure FcnForm is replaced by the general container structure FieldExpression.
   *   This means that any fully-applied anonymous function can be treated like a FieldExpression and used wherever a FieldExpression is appropriate/acceptable.
   *   The extended example below demonstrates this replacement.
   *  
   *  @par Currying
   *   Currying is the idea that not all of a functions parameters must be applied at the same time.
   *   The current implementation actually curries each argument separately.
   *   So in the running example, the anonymous function
   *  \code
   *  ($1 + $0)(4, a)
   *  \endcode
   *   is evaluated as if each argument was applied individually:
   *  \code
   *  ($1 + $0)(4)(a)
   *  \endcode
   *  
   *  @par Extended example
   *   Consider the running example:
   *  \code
   *  ($1 + $0)(4, a)
   *  \endcode
   *   Before any application occurs, the anonymous function itself,
   *  \code
   *  ($1 + $0)
   *  \endcode
   *   has the following type:
   *  \code
   *  FcnForm<SumOp<ArgForm<1, FT>,
   *                ArgForm<0, FT>,
   *                FT>,
   *          0,
   *          2,
   *          FT>
   *  \endcode
   *  When the arguments are applied they are applied individually, so \c 4 is applied as the first argument without considering the second argument \c a yet.
   *  After \c 4 is applied, the code is equivalent to:
   *  \code
   *  ($1 + 4)(a)
   *  \endcode
   *   The anonymous function itself,
   *  \code
   *  ($1 + 4)
   *  \endcode
   *   now has the following type:
   *  \code
   *  FcnForm<SumOp<ArgForm<1, FT>,
   *                Scalar<AT>,
   *                FT>,
   *          1,
   *          2,
   *          FT>
   *  \endcode
   *   All references to \c ArgForm<0, \c FT> have been replaced by the FieldExpression, \c Scalar<AT>.
   *   Also, the CrtNum template parameter has been incremented, signifying that the second parameter is the next parameter to be applied.
   *  After \c a is applied, the code is equivalent to:
   *  \code
   *  (a + 4)
   *  \endcode
   *   The anonymous function itself has disappeared (no anonymous arguments left) and is replaced by a FieldExpression with the following type:
   *  \code
   *  FieldExpression<SumOp<FieldForm<FT>,
   *                        Scalar<AT>,
   *                        FT>,
   *                  FT>
   *  \endcode
   *   So the original anonymous function,
   *  \code
   *  ($1 + $0)(4, a)
   *  \endcode
   *   after both its arguments have been applied returns a FieldExpression, namely:
   *  \code
   *  FieldExpression<SumOp<FieldForm<FT>,
   *                        Scalar<AT>,
   *                        FT>,
   *                  FT>
   *  \endcode
   *
   *   \warning The current implementation supports applying at most three arguments at a time.
   *   But there can be any number of different applications.
   *   For example, the following anonymous function will not compile:
   *  \code
   *  ($0 + $1 + $2 + $3)(4, a, b, 3)
   *  \endcode
   *   However, the following equivalent anonymous functions will compile:
   *  \code
   *  ($0 + $1 + $2 + $3)(4, a, b)(3)
   *  \endcode
   *  \code
   *  ($0 + $1 + $2 + $3)(4, a)(b, 3)
   *  \endcode
   *  \code
   *  ($0 + $1 + $2 + $3)(4)(a)(b)(3)
   *  \endcode
   */
  template<typename ExprType, int CrtNum, int MaxNum, typename FieldType>
    struct FcnForm<ExprType,CrtNum,Int<MaxNum>,FieldType> {
  public:
    FieldType typedef field_type;
    ExprType typedef expr_type;
    
  private:
    ExprType expr;
    
  public:
    /**
     *  @brief Constructor for FcnForm.
     *  
     *  \param given an ExprType that is the internal expression type.
     *  \return a new FcnForm containing given.
     *  
     *  Constructs a new FcnForm containing given.
     */
  FcnForm(ExprType given)
    :expr(given)
    {};
    
    /**
     *  @brief Accessor for internal expression.
     *
     *  \return constant reference to internal expression.
     *
     *  Returns constant reference to internal expression.
     */
    ExprType const & expression() const {
      return expr;
    };
    
    /**
     *  @brief Single argument function application.
     *
     *  \param arg an ArgType object.
     *  \return either a FcnForm or a FieldExpression, depending on the number of remaining arguments to apply.
     *  
     *  The return type is based on type-juggling to make sure everything matches up.
     *  If CrtNum is two less than MaxNum (one less after application), then there are no more arguments to apply, and a FieldExpression object is returned.
     *  However, if CrtNum is more than two less than MaxNum, then there are more arguments to apply, and a FcnForm object is returned.
     */
    template<typename ArgType>
      typename ArgResultTerm<typename ArgReplace<ExprType,
      CrtNum,
      typename StandardizeArg<ArgType,
      FieldType>::StandardType>::ResultType,
      CrtNum,
      //CrtNum + 1,
      Int<MaxNum>,
      FieldType>::ResultType operator () (ArgType const & arg) {
    
      /* Actual type of expression after application of (standardized) ArgType. */
      ArgReplace<ExprType,
	CrtNum,
	typename StandardizeArg<ArgType,
	FieldType>::StandardType> typedef ActualNextType;
      
      /* Wrapper type - if all arguments have been bound, wrapper is FieldExpression; otherwise, wrapper is FcnForm. */
      typename ArgResultTerm<typename ActualNextType::ResultType,
	CrtNum,
	Int<MaxNum>,
	FieldType>::ResultType typedef WrapperNextType;
      
      /* Actual code that runs: Call correct apply function followed by a typecast. */
      return WrapperNextType(ActualNextType::apply(expr,
						   StandardizeArg<ArgType,
						   FieldType>::standardType(arg)));
    };
  
    /*
     * Multiple arguments/currying the uncurried:
     */
    
    /**
     *  @brief Double argument function application.
     *
     *  \param arg1 an Arg1 object.
     *  \param arg2 an Arg2 object.
     *  \return either a FcnForm or a FieldExpression, depending on the number of remaining arguments to apply.
     *  
     *  The return type is based on type-juggling to make sure everything matches up.
     *  If CrtNum is three less than MaxNum (one less after application), then there are no more arguments to apply, and a FieldExpression object is returned.
     *  However, if CrtNum is more than three less than MaxNum, then there are more arguments to apply, and a FcnForm object is returned.
     */
    template<typename Arg1, typename Arg2>
      typename ArgResultTerm<typename ArgReplace<typename ArgReplace<ExprType,
      CrtNum,
      typename StandardizeArg<Arg1,
      FieldType>::StandardType>::ResultType,
      CrtNum + 1,
      typename StandardizeArg<Arg2,
      FieldType>::StandardType>::ResultType,
      CrtNum + 1,
      Int<MaxNum>,
      FieldType>::ResultType operator () (Arg1 const & arg1,
					  Arg2 const & arg2) {
      return this -> operator ()
	(arg1)
	(arg2);
    };
    
    /**
     *  @brief Triple argument function application.
     *
     *  \param arg1 an Arg1 object.
     *  \param arg2 an Arg2 object.
     *  \param arg3 an Arg3 object.
     *  \return either a FcnForm or a FieldExpression, depending on the number of remaining arguments to apply.
     *  
     *  The return type is based on type-juggling to make sure everything matches up.
     *  If CrtNum is four less than MaxNum (one less after application), then there are no more arguments to apply, and a FieldExpression object is returned.
     *  However, if CrtNum is more than four less than MaxNum, then there are more arguments to apply, and a FcnForm object is returned.
     */
    template<typename Arg1, typename Arg2, typename Arg3>
      typename ArgResultTerm<typename ArgReplace<typename ArgReplace<typename ArgReplace<ExprType,
      CrtNum,
      typename StandardizeArg<Arg1,
      FieldType>::StandardType>::ResultType,
      CrtNum + 1,
      typename StandardizeArg<Arg2,
      FieldType>::StandardType>::ResultType,
      CrtNum + 2,
      typename StandardizeArg<Arg3,
      FieldType>::StandardType>::ResultType,
      CrtNum + 2,
      Int<MaxNum>,
      FieldType>::ResultType operator () (Arg1 const & arg1,
					  Arg2 const & arg2,
					  Arg3 const & arg3) {
      return this -> operator ()
	(arg1)
	(arg2)
	(arg3);
    };
  };
  
  /**
   *  @struct Int
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Wrapper/container structure to make an integer available at compile-time.
   *  
   *  \tparam Num An integer.
   *  
   *  Int contains an integer as a template parameter to make that integer available at compile-time for various uses.
   *  Int should never be instantiated; its sole purpose is making an integer into a type.
   */
  template<int Num>
    struct Int {};

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the meta-function, max.
   *  
   *  \tparam Num1 An integer.
   *  \tparam Num2 An integer.
   *  
   *  MetaFcnMax encodes the max function (returns maximum of arguments) as a structure so that the result is available at compile-time.
   *  
   *  The arguments for max are the template parameters: Num1 and Num2, wrapped in an Int.
   *  The result is in the typedef'ed Max.
   *  The computation of the maximum argument is preformed in conjunction with InternalMetaFcnMax.
   *  
   *  MetaFcnMax should never be instantiated; its sole purpose is computing the maximum of two integers at compile-time.
   */
  template<int Num1, int Num2>
    struct MetaFcnMax<Int<Num1>,Int<Num2> > {
  
  private:
    template<typename TrueAnswer, typename FalseAnswer, bool answer>
      struct InternalMetaFcnMax;
  
  public:
    /**
     *  @brief Typedef for Max, the result of the meta-function, max.
     *  
     *  Typedef for Max, which is either \c Int<Num1> or \c Int<Num2>, whichever is larger.
     */
    typename InternalMetaFcnMax<Int<Num1>,Int<Num2>,(Num1 > Num2)>::Max typedef Max;
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the internals of meta-function, max. (True branch)
   *  
   *  \tparam Num1 An integer (from MetaFcnMax).
   *  \tparam Num2 An integer (from MetaFcnMax).
   *  \tparam TrueAnswer A type (Int<Num1>).
   *  \tparam FalseAnswer A type (Int<Num2>).
   *  \tparam answer A boolean (here \c true).
   *  
   *  InternalMetaFcnMax encodes the internals of the meta-function max (returns maximum of arguments) as a structure so that the result is available at compile-time.
   *  
   *  The arguments for max are the first two template parameters: Num1 and Num2, wrapped in an Int.
   *  The boolean result of the conditional (is the first argument larger than the second) is the third template parameter, here \c true.
   *  The specializations of InternalMetaFcnMax each cover a different branch of the conditional's result.
   *  The result is in the typedef'ed Max, here \c TrueAnswer (\c Int<Num1>).
   *  The computation of the maximum argument is preformed in conjunction with MetaFcnMax.
   *  
   *  InternalMetaFcnMax should never be instantiated; its sole purpose is computing the maximum of two integers at compile-time.
   */
  template<int Num1, int Num2>
    template<typename TrueAnswer, typename FalseAnswer>
    struct MetaFcnMax<Int<Num1>,Int<Num2> >::InternalMetaFcnMax<TrueAnswer,FalseAnswer,true> {
    TrueAnswer typedef Max;
  };
    
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the internals of meta-function, max. (False branch)
   *  
   *  \tparam Num1 An integer (from MetaFcnMax).
   *  \tparam Num2 An integer (from MetaFcnMax).
   *  \tparam TrueAnswer A type (Int<Num1>).
   *  \tparam FalseAnswer A type (Int<Num2>).
   *  \tparam answer A boolean (here \c false).
   *  
   *  InternalMetaFcnMax encodes the internals of the meta-function max (returns maximum of arguments) as a structure so that the result is available at compile-time.
   *  
   *  The arguments for max are the first two template parameters: Num1 and Num2, wrapped in an Int.
   *  The boolean result of the conditional (is the first argument larger than the second) is the third template parameter, here \c false.
   *  The specializations of InternalMetaFcnMax each cover a different branch of the conditional's result.
   *  The result is in the typedef'ed Max, here \c FalseAnswer (\c Int<Num2>).
   *  The computation of the maximum argument is preformed in conjunction with MetaFcnMax.
   *  
   *  InternalMetaFcnMax should never be instantiated; its sole purpose is computing the maximum of two integers at compile-time.
   */
  template<int Num1, int Num2>
    template<typename TrueAnswer, typename FalseAnswer>
    struct MetaFcnMax<Int<Num1>,Int<Num2> >::InternalMetaFcnMax<TrueAnswer,FalseAnswer,false> {
    FalseAnswer typedef Max;
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the StandardizeTerm type conversions (type juggling) for intermediate FieldExpression structures. (FieldType case).
   *  
   *  \tparam Input A type (here FieldType).
   *  \tparam FieldType Field type.
   *  
   *  StandardizeTerm provides the standard FieldExpression terms and types for the base types that underlie FieldExpression objects.
   *  
   *  StandardizeTerm preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> Input </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> FieldType </td>
   *  <td> FieldExpression<FieldForm<...>, ...> </td>
   *  <td> FieldForm<...> </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<...> </td>
   *  <td> FcnForm<ArgForm<...>, ...> </td>
   *  <td> ArgForm<...> </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType, ...> </td>
   *  <td> FcnForm<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  </tr>
   *  </table>
   *  
   *  A StandardTerm is either a FieldExpression or a FcnForm, based on whether the expression represented does not or does contain anonymous arguements, respectively.
   *  That is: If an underlying expression contains any anonymous arguments (ArgForm), its StandardTerm is a FcnForm.
   *  Conversely, if an underlying expression contains no anonymous arguments (ArgForm), its StandardTerm is a FieldExpression.
   *  A StandardType is always the underlying expression in either a FieldExpression or FcnForm.
   *  A FieldType (as here) is converted into its FieldExpression-style representation, FieldForm.
   *  
   *  Note that Scalar objects are not handled by StandardizeTerm.  Scalars and their underlying representation (developer/user defined) do not contain FieldType, so they are handled differently.
   *  
   *  StandardizeTerm should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename FieldType>
    struct StandardizeTerm<FieldType,FieldType> {
    
    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here FieldForm<FieldType>.
     */
    FieldForm<FieldType> typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FieldExpression<StandardType,FieldType>.
     */
    FieldExpression<StandardType, FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given An Input object (here a FieldType object).
     *  \return a new StandardType object (here a FieldForm object).
     *  
     *  Converts and returns a StandardType object based on given, (here a FieldForm object).
     */
    static inline StandardType standardType (FieldType const & given) {
      return StandardType(given);
    };
    
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given An Input object (here a FieldType object).
     *  \return a new StandardTerm (here a FieldExpression object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FieldExpression<FieldForm<...>, ...> object).
     */
    static inline StandardTerm standardTerm (FieldType const & given) {
      return StandardTerm(StandardType(given));
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the StandardizeTerm type conversions (type juggling) for intermediate FieldExpression structures. (FieldExpression case).
   *  
   *  \tparam Input A type (here FieldExpression).
   *  \tparam FieldType Field type.
   *  
   *  StandardizeTerm provides the standard FieldExpression terms and types for the base types that underlie FieldExpression objects.
   *  
   *  StandardizeTerm preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> Input </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> FieldType </td>
   *  <td> FieldExpression<FieldForm<...>, ...> </td>
   *  <td> FieldForm<...> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<...> </td>
   *  <td> FcnForm<ArgForm<...>, ...> </td>
   *  <td> ArgForm<...> </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType, ...> </td>
   *  <td> FcnForm<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  </tr>
   *  </table>
   *  
   *  A StandardTerm is either a FieldExpression or a FcnForm, based on whether the expression represented does not or does contain anonymous arguements, respectively.
   *  That is: If an underlying expression contains any anonymous arguments (ArgForm), its StandardTerm is a FcnForm.
   *  Conversely, if an underlying expression contains no anonymous arguments (ArgForm), its StandardTerm is a FieldExpression.
   *  A StandardType is always the underlying expression in either a FieldExpression or FcnForm.
   *  A FieldType is converted into its FieldExpression-style representation, FieldForm.
   *  
   *  Note that Scalar objects are not handled by StandardizeTerm.  Scalars and their underlying representation (developer/user defined) do not contain FieldType, so they are handled differently.
   *  
   *  StandardizeTerm should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename ExprType, typename FieldType>
    struct StandardizeTerm<FieldExpression<ExprType,FieldType>,FieldType> {
    
    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here ExprType.
     */
    ExprType typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FieldExpression<StandardType,FieldType>.
     */
    FieldExpression<StandardType,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given An Input object (here a FieldExpression<ExprType, ...> object).
     *  \return a reference to a StandardType object (here a ExprType object).
     *  
     *  Converts and returns a reference to a StandardType object based on given, (here a ExprType object).
     */
    static inline StandardType const & standardType (FieldExpression<ExprType,FieldType> const & given) {
      return given.expression();
    };
    
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given An Input object (here a FieldExpression<ExprType, ...> object).
     *  \return a reference to a StandardTerm object (here a FieldExpression<ExprType, ...> object).
     *  
     *  Converts and returns a reference to a StandardTerm object, given, (here a FieldExpression<ExprType, ...> object).
     */
    static inline StandardTerm const & standardTerm (FieldExpression<ExprType,FieldType> const & given) {
      return given;
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the StandardizeTerm type conversions (type juggling) for intermediate FieldExpression structures. (ArgForm case).
   *  
   *  \tparam Input A type (here ArgForm).
   *  \tparam FieldType Field type.
   *  
   *  StandardizeTerm provides the standard FieldExpression terms and types for the base types that underlie FieldExpression objects.
   *  
   *  StandardizeTerm preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> Input </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> FieldType </td>
   *  <td> FieldExpression<FieldForm<...>, ...> </td>
   *  <td> FieldForm<...> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<...> </td>
   *  <td> FcnForm<ArgForm<...>, ...> </td>
   *  <td> ArgForm<...> </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType, ...> </td>
   *  <td> FcnForm<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  </tr>
   *  </table>
   *  
   *  A StandardTerm is either a FieldExpression or a FcnForm, based on whether the expression represented does not or does contain anonymous arguements, respectively.
   *  That is: If an underlying expression contains any anonymous arguments (ArgForm), its StandardTerm is a FcnForm.
   *  Conversely, if an underlying expression contains no anonymous arguments (ArgForm), its StandardTerm is a FieldExpression.
   *  A StandardType is always the underlying expression in either a FieldExpression or FcnForm.
   *  A FieldType is converted into its FieldExpression-style representation, FieldForm.
   *  
   *  Note that Scalar objects are not handled by StandardizeTerm.  Scalars and their underlying representation (developer/user defined) do not contain FieldType, so they are handled differently.
   *  
   *  StandardizeTerm should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<int Num, typename FieldType>
    struct StandardizeTerm<ArgForm<Num,FieldType>,FieldType> {
    
    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here ArgForm<...>.
     */
    ArgForm<Num,FieldType> typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FcnForm<StandardType,...>.
     */
    FcnForm<StandardType,0,Int<Num + 1>,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given An Input object (here a ArgForm<...> object).
     *  \return a reference to a StandardType object (here a ArgForm<...> object).
     *  
     *  Converts and returns a reference to a StandardType object, given, (here a ExprType object).
     */
    static inline StandardType const & standardType (ArgForm<Num,FieldType> const & given) {
      return given;
    };
    
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given An Input object (here a ArgForm<...> object).
     *  \return a new StandardTerm (here a FcnForm<ArgForm<...>, ...> object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FcnForm<ArgForm<...>, ...> object).
     */
    static inline StandardTerm standardTerm (ArgForm<Num,FieldType> const & given) {
      return StandardTerm(given);
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the StandardizeTerm type conversions (type juggling) for intermediate FieldExpression structures. (FcnForm case).
   *  
   *  \tparam Input A type (here FcnForm).
   *  \tparam FieldType Field type.
   *  
   *  StandardizeTerm provides the standard FieldExpression terms and types for the base types that underlie FieldExpression objects.
   *  
   *  StandardizeTerm preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> Input </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> FieldType </td>
   *  <td> FieldExpression<FieldForm<...>, ...> </td>
   *  <td> FieldForm<...> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<...> </td>
   *  <td> FcnForm<ArgForm<...>, ...> </td>
   *  <td> ArgForm<...> </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType, ...> </td>
   *  <td> FcnForm<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  </table>
   *  
   *  A StandardTerm is either a FieldExpression or a FcnForm, based on whether the expression represented does not or does contain anonymous arguements, respectively.
   *  That is: If an underlying expression contains any anonymous arguments (ArgForm), its StandardTerm is a FcnForm.
   *  Conversely, if an underlying expression contains no anonymous arguments (ArgForm), its StandardTerm is a FieldExpression.
   *  A StandardType is always the underlying expression in either a FieldExpression or FcnForm.
   *  A FieldType is converted into its FieldExpression-style representation, FieldForm.
   *  
   *  Note that Scalar objects are not handled by StandardizeTerm.  Scalars and their underlying representation (developer/user defined) do not contain FieldType, so they are handled differently.
   *  
   *  StandardizeTerm should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename ExprType, typename Max, typename FieldType>
    struct StandardizeTerm<FcnForm<ExprType,0,Max,FieldType>,FieldType> {
    
    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here ExprType.
     */
    ExprType typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FcnForm<StandardType, ...>.
     */
    FcnForm<StandardType,0,Max,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given An Input object (here a FcnForm<ExprType, ...> object).
     *  \return a reference to a StandardType object (here a ExprType object).
     *  
     *  Converts and returns a reference to a StandardType object based on given, (here a ExprType object).
     */
    static inline StandardType const & standardType (FcnForm<ExprType,0,Max,FieldType> const & given) {
      return given.expression();
    };
    
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given An Input object (here a FcnForm<ExprType, ...> object).
     *  \return a reference to a StandardTerm object (here a FcnForm<ExprType, ...> object).
     *  
     *  Converts and returns a reference to a StandardTerm object, given, (here a FcnForm<ExprType, ...> object).
     */
    static inline StandardTerm const & standardTerm (FcnForm<ExprType,0,Max,FieldType> const & given) {
      return given;
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the LiftTerm type conversions (type juggling) for intermediate FieldExpression structures. (FieldExpression case).
   *  
   *  \tparam NewType A type.
   *  \tparam OldTerm A Term (here FieldExpression).
   *  \tparam FieldType Field type.
   *  
   *  LiftTerm creates the StandardTerm that matches OldTerm (here a FieldExpression), but contains NewType as its ExprType.
   *  
   *  LiftTerm preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> OldTerm </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<OldType, ...> </td>
   *  <td> FieldExpression<NewType, ...> </td>
   *  <td> NewType </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<OldType, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  </table>
   *  
   *  A StandardTerm is either a FieldExpression or a FcnForm, based on whether the expression represented does not or does contain anonymous arguements, respectively.
   *  That is: If an underlying expression contains any anonymous arguments (ArgForm), its StandardTerm is a FcnForm.
   *  Conversely, if an underlying expression contains no anonymous arguments (ArgForm), its StandardTerm is a FieldExpression.
   *  A StandardType is always the underlying expression in either a FieldExpression or FcnForm.
   *  
   *  Note that LiftTerm does not check that these conventions are maintained.  It is assumed that NewType is somehow based upon OldType, and therefore these conventions hold true.
   *  
   *  LiftTerm should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename NewType, typename OldType, typename FieldType>
    struct LiftTerm<NewType,
    FieldExpression<OldType,FieldType>,
    FieldType> {
    
    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here NewType.
     */
    NewType typedef StandardType;

    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FieldExpression<StandardType, ...>.
     */
    FieldExpression<NewType,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given A NewType object.
     *  \return a reference to a StandardType object.
     *  
     *  Converts and returns a reference to a StandardType object, given.
     */
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given A NewType object.
     *  \return a StandardTerm object (here a FieldExpression<NewType, ...> object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FieldExpression<NewType, ...> object).
     */
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the LiftTerm type conversions (type juggling) for intermediate FieldExpression structures. (FcnForm case).
   *  
   *  \tparam NewType A type.
   *  \tparam OldTerm A Term (here FcnForm).
   *  \tparam FieldType Field type.
   *  
   *  LiftTerm creates the StandardTerm that matches OldTerm (here a FcnForm), but contains NewType as its ExprType.
   *  
   *  LiftTerm preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> OldTerm </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<OldType, ...> </td>
   *  <td> FieldExpression<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<OldType, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  </table>
   *  
   *  A StandardTerm is either a FieldExpression or a FcnForm, based on whether the expression represented does not or does contain anonymous arguements, respectively.
   *  That is: If an underlying expression contains any anonymous arguments (ArgForm), its StandardTerm is a FcnForm.
   *  Conversely, if an underlying expression contains no anonymous arguments (ArgForm), its StandardTerm is a FieldExpression.
   *  A StandardType is always the underlying expression in either a FieldExpression or FcnForm.
   *  
   *  Note that LiftTerm does not check that these conventions are maintained.  It is assumed that NewType is somehow based upon OldType, and therefore these conventions hold true.
   *  
   *  LiftTerm should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename NewType, typename OldType, typename Max, typename FieldType>
    struct LiftTerm<NewType,
    FcnForm<OldType,0,Max,FieldType>,
    FieldType> {
    
    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here NewType.
     */
    NewType typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FcnForm<StandardType, ...>.
     */
    FcnForm<NewType,0,Max,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given A NewType object.
     *  \return a reference to a StandardType object.
     *  
     *  Converts and returns a reference to a StandardType object, given.
     */
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given A NewType object.
     *  \return a StandardTerm object (here a FcnForm<NewType, ...> object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FcnForm<NewType, ...> object).
     */
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the CombineTerms type conversions (type juggling) for intermediate FieldExpression structures. (two FieldExpression types case.)
   *  
   *  \tparam NewType A type.
   *  \tparam OldTerm1 A Term (here FieldExpression).
   *  \tparam OldTerm2 A Term (here FieldExpression).
   *  \tparam FieldType Field type.
   *  
   *  CombineTerms creates the StandardTerm that matches the combination of OldTerm1 and OldTerm2 (here both are FieldExpression types), but contains NewType as its ExprType.
   *  CombineTerms is a binary version of LiftTerm: In LiftTerm, NewType is based upon a single ExprType; whereas in CombineTerms, NewType is based upon two ExprTypes.
   *  
   *  CombineTerms preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> OldTerm1 </th>
   *  <th> OldTerm2 </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType1, ...> </td>
   *  <td> FieldExpression<ExprType2, ...> </td>
   *  <td> FieldExpression<NewType, ...> </td>
   *  <td> NewType </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType1, ...> </td>
   *  <td> FcnForm<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType1, ...> </td>
   *  <td> FieldExpression<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType1, ...> </td>
   *  <td> FcnForm<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  </table>
   *  
   *  A StandardTerm is either a FieldExpression or a FcnForm, based on whether the expression represented does not or does contain anonymous arguements, respectively.
   *  That is: If an underlying expression contains any anonymous arguments (ArgForm), its StandardTerm is a FcnForm.
   *  Conversely, if an underlying expression contains no anonymous arguments (ArgForm), its StandardTerm is a FieldExpression.
   *  A StandardType is always the underlying expression in either a FieldExpression or FcnForm.
   *  
   *  Note that CombineTerms does not check that these conventions are maintained.  It is assumed that NewType is somehow based upon ExprType1 and ExprType2, and therefore these conventions hold true.
   *  
   *  CombineTerms should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename NewType, typename ExprType1, typename ExprType2, typename FieldType>
    struct CombineTerms<NewType,
    FieldExpression<ExprType1,FieldType>,
    FieldExpression<ExprType2,FieldType>,
    FieldType> {

    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here NewType.
     */
    NewType typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FieldExpression<StandardType, ...>.
     */
    FieldExpression<NewType,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given A NewType object.
     *  \return a reference to a StandardType object.
     *  
     *  Converts and returns a reference to a StandardType object, given.
     */
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
    
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given A NewType object.
     *  \return a StandardTerm object (here a FieldExpression<NewType, ...> object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FieldExpression<NewType, ...> object).
     */
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the CombineTerms type conversions (type juggling) for intermediate FieldExpression structures. (FieldExpression and FcnForm case.)
   *  
   *  \tparam NewType A type.
   *  \tparam OldTerm1 A Term (here FieldExpression).
   *  \tparam OldTerm2 A Term (here FcnForm).
   *  \tparam FieldType Field type.
   *  
   *  CombineTerms creates the StandardTerm that matches the combination of OldTerm1 and OldTerm2 (here a FieldExpression and a FcnForm, respectively), but contains NewType as its ExprType.
   *  CombineTerms is a binary version of LiftTerm: In LiftTerm, NewType is based upon a single ExprType; whereas in CombineTerms, NewType is based upon two ExprTypes.
   *  
   *  CombineTerms preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> OldTerm1 </th>
   *  <th> OldTerm2 </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType1, ...> </td>
   *  <td> FieldExpression<ExprType2, ...> </td>
   *  <td> FieldExpression<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType1, ...> </td>
   *  <td> FcnForm<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType1, ...> </td>
   *  <td> FieldExpression<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType1, ...> </td>
   *  <td> FcnForm<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  </table>
   *  
   *  A StandardTerm is either a FieldExpression or a FcnForm, based on whether the expression represented does not or does contain anonymous arguements, respectively.
   *  That is: If an underlying expression contains any anonymous arguments (ArgForm), its StandardTerm is a FcnForm.
   *  Conversely, if an underlying expression contains no anonymous arguments (ArgForm), its StandardTerm is a FieldExpression.
   *  A StandardType is always the underlying expression in either a FieldExpression or FcnForm.
   *  
   *  Note that CombineTerms does not check that these conventions are maintained.  It is assumed that NewType is somehow based upon ExprType1 and ExprType2, and therefore these conventions hold true.
   *  
   *  CombineTerms should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename NewType, typename ExprType1, typename ExprType2, typename Max2, typename FieldType>
    struct CombineTerms<NewType,
    FieldExpression<ExprType1,FieldType>,
    FcnForm<ExprType2,0,Max2,FieldType>,
    FieldType> {

    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here NewType.
     */
    NewType typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FcnForm<StandardType, ...>.
     */
    FcnForm<NewType,0,Max2,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given A NewType object.
     *  \return a reference to a StandardType object.
     *  
     *  Converts and returns a reference to a StandardType object, given.
     */
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
    
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given A NewType object.
     *  \return a StandardTerm object (here a FcnForm<NewType, ...> object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FcnForm<NewType, ...> object).
     */
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the CombineTerms type conversions (type juggling) for intermediate FieldExpression structures. (FcnForm and FieldExpression case.)
   *  
   *  \tparam NewType A type.
   *  \tparam OldTerm1 A Term (here FcnForm).
   *  \tparam OldTerm2 A Term (here FieldExpression).
   *  \tparam FieldType Field type.
   *  
   *  CombineTerms creates the StandardTerm that matches the combination of OldTerm1 and OldTerm2 (here a FcnForm and a FieldExpression, respectively), but contains NewType as its ExprType.
   *  CombineTerms is a binary version of LiftTerm: In LiftTerm, NewType is based upon a single ExprType; whereas in CombineTerms, NewType is based upon two ExprTypes.
   *  
   *  CombineTerms preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> OldTerm1 </th>
   *  <th> OldTerm2 </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType1, ...> </td>
   *  <td> FieldExpression<ExprType2, ...> </td>
   *  <td> FieldExpression<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType1, ...> </td>
   *  <td> FcnForm<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType1, ...> </td>
   *  <td> FieldExpression<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType1, ...> </td>
   *  <td> FcnForm<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  </table>
   *  
   *  A StandardTerm is either a FieldExpression or a FcnForm, based on whether the expression represented does not or does contain anonymous arguements, respectively.
   *  That is: If an underlying expression contains any anonymous arguments (ArgForm), its StandardTerm is a FcnForm.
   *  Conversely, if an underlying expression contains no anonymous arguments (ArgForm), its StandardTerm is a FieldExpression.
   *  A StandardType is always the underlying expression in either a FieldExpression or FcnForm.
   *  
   *  Note that CombineTerms does not check that these conventions are maintained.  It is assumed that NewType is somehow based upon ExprType1 and ExprType2, and therefore these conventions hold true.
   *  
   *  CombineTerms should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename NewType, typename ExprType1, typename Max1, typename ExprType2, typename FieldType>
    struct CombineTerms<NewType,
    FcnForm<ExprType1,0,Max1,FieldType>,
    FieldExpression<ExprType2,FieldType>,
    FieldType> {

    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here NewType.
     */
    NewType typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FcnForm<StandardType, ...>.
     */
    FcnForm<NewType,0,Max1,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given A NewType object.
     *  \return a reference to a StandardType object.
     *  
     *  Converts and returns a reference to a StandardType object, given.
     */
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given A NewType object.
     *  \return a StandardTerm object (here a FcnForm<NewType, ...> object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FcnForm<NewType, ...> object).
     */
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the CombineTerms type conversions (type juggling) for intermediate FieldExpression structures. (Two FcnForm types case.)
   *  
   *  \tparam NewType A type.
   *  \tparam OldTerm1 A Term (here FcnForm).
   *  \tparam OldTerm2 A Term (here FcnForm).
   *  \tparam FieldType Field type.
   *  
   *  CombineTerms creates the StandardTerm that matches the combination of OldTerm1 and OldTerm2 (here both are FcnForm types), but contains NewType as its ExprType.
   *  CombineTerms is a binary version of LiftTerm: In LiftTerm, NewType is based upon a single ExprType; whereas in CombineTerms, NewType is based upon two ExprTypes.
   *  
   *  CombineTerms preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> OldTerm1 </th>
   *  <th> OldTerm2 </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType1, ...> </td>
   *  <td> FieldExpression<ExprType2, ...> </td>
   *  <td> FieldExpression<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType1, ...> </td>
   *  <td> FcnForm<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType1, ...> </td>
   *  <td> FieldExpression<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  </tr>
   *  <tr>
   *  <td> FcnForm<ExprType1, ...> </td>
   *  <td> FcnForm<ExprType2, ...> </td>
   *  <td> FcnForm<NewType, ...> </td>
   *  <td> NewType </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  </table>
   *  
   *  A StandardTerm is either a FieldExpression or a FcnForm, based on whether the expression represented does not or does contain anonymous arguements, respectively.
   *  That is: If an underlying expression contains any anonymous arguments (ArgForm), its StandardTerm is a FcnForm.
   *  Conversely, if an underlying expression contains no anonymous arguments (ArgForm), its StandardTerm is a FieldExpression.
   *  A StandardType is always the underlying expression in either a FieldExpression or FcnForm.
   *  
   *  Note that CombineTerms does not check that these conventions are maintained.  It is assumed that NewType is somehow based upon ExprType1 and ExprType2, and therefore these conventions hold true.
   *  
   *  CombineTerms should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename NewType, typename ExprType1, typename Max1, typename ExprType2, typename Max2, typename FieldType>
    struct CombineTerms<NewType,
    FcnForm<ExprType1,0,Max1,FieldType>,
    FcnForm<ExprType2,0,Max2,FieldType>,
    FieldType> {

    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here NewType.
     */
    NewType typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FcnForm<StandardType, ...>.
     */
    FcnForm<NewType,
      0,
      typename MetaFcnMax<Max1,Max2>::Max,
      FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given A NewType object.
     *  \return a reference to a StandardType object.
     *  
     *  Converts and returns a reference to a StandardType object, given.
     */
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
    
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given A NewType object.
     *  \return a StandardTerm object (here a FcnForm<NewType, ...> object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FcnForm<NewType, ...> object).
     */
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the StandardizeArg type conversions (type juggling) for intermediate FieldExpression structures. (Scalar case.)
   *  
   *  \tparam Arg A basic type.
   *  \tparam FieldType Field type.
   *  
   *  StandardizeArg creates the StandardTerm (FieldExpression) that matches Arg type.
   *  StandardizeArg is used to convert the actual arguments/parameters of anonymous arguments (ArgForm) to their FieldExpression representations.
   *  
   *  StandardizeArg preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> Arg </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> AtomicType / value_type </td>
   *  <td> FieldExpression<Scalar<...>, ...> </td>
   *  <td> Scalar<...> </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> FieldType </td>
   *  <td> FieldExpression<FieldForm<...>, ...> </td>
   *  <td> FieldForm<...> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  </tr>
   *  </table>
   *  
   *  While a StandardTerm can be either a FieldExpression or a FcnForm, StandardizeArg standardizes the potential arguments to FcnForm objects.
   *  Since FcnForm objects can only take FieldExpression objects as actual arguments, StandardizeArg only produces FieldExpression objects.
   *  
   *  StandardizeArg should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename FieldType>
    struct StandardizeArg<typename FieldType::value_type,FieldType> {
    typename FieldType::value_type typedef AtomicType;

    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here Scalar<...>.
     */
    Scalar<AtomicType> typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FieldExpression<StandardType, ...>.
     */
    FieldExpression<StandardType,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given An AtomicType object.
     *  \return a new StandardType object.
     *  
     *  Builds and returns a StandardType object based on given.
     */
    static inline StandardType standardType (AtomicType const & given) {
      return StandardType(given);
    };
    
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given An AtomicType object.
     *  \return a StandardTerm object (here a FieldExpression<Scalar<...>, ...> object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FieldExpression<Scalar<...>, ...> object).
     */
    static inline StandardTerm standardTerm (AtomicType const & given) {
      return StandardTerm(StandardType(given));
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the StandardizeArg type conversions (type juggling) for intermediate FieldExpression structures. (FieldForm case.)
   *  
   *  \tparam Arg A basic type.
   *  \tparam FieldType Field type.
   *  
   *  StandardizeArg creates the StandardTerm (FieldExpression) that matches Arg type.
   *  StandardizeArg is used to convert the actual arguments/parameters of anonymous arguments (ArgForm) to their FieldExpression representations.
   *  
   *  StandardizeArg preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> Arg </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> AtomicType / value_type </td>
   *  <td> FieldExpression<Scalar<...>, ...> </td>
   *  <td> Scalar<...> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldType </td>
   *  <td> FieldExpression<FieldForm<...>, ...> </td>
   *  <td> FieldForm<...> </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  </tr>
   *  </table>
   *  
   *  While a StandardTerm can be either a FieldExpression or a FcnForm, StandardizeArg standardizes the potential arguments to FcnForm objects.
   *  Since FcnForm objects can only take FieldExpression objects as actual arguments, StandardizeArg only produces FieldExpression objects.
   *  
   *  StandardizeArg should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename FieldType>
    struct StandardizeArg<FieldType,FieldType> {

    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here FieldForm<...>.
     */
    FieldForm<FieldType> typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FieldExpression<StandardType, ...>.
     */
    FieldExpression<StandardType,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given A FieldType object.
     *  \return a new StandardType object.
     *  
     *  Builds and returns a StandardType object based on given.
     */
    static inline StandardType standardType (FieldType const & given) {
      return StandardType(given);
    };
    
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given A FieldType object.
     *  \return a StandardTerm object (here a FieldExpression<FieldForm<...>, ...> object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FieldExpression<FieldForm<...>, ...> object).
     */
    static inline StandardTerm standardTerm (FieldType const & given) {
      return StandardTerm(StandardType(given));
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the StandardizeArg type conversions (type juggling) for intermediate FieldExpression structures. (FieldExpression case.)
   *  
   *  \tparam Arg A basic type.
   *  \tparam FieldType Field type.
   *  
   *  StandardizeArg creates the StandardTerm (FieldExpression) that matches Arg type.
   *  StandardizeArg is used to convert the actual arguments/parameters of anonymous arguments (ArgForm) to their FieldExpression representations.
   *  
   *  StandardizeArg preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> Arg </th>
   *  <th> Standard Term </th>
   *  <th> Standard Type </th>
   *  </tr>
   *  <tr>
   *  <td> AtomicType / value_type </td>
   *  <td> FieldExpression<Scalar<...>, ...> </td>
   *  <td> Scalar<...> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldType </td>
   *  <td> FieldExpression<FieldForm<...>, ...> </td>
   *  <td> FieldForm<...> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> FieldExpression<ExprType, ...> </td>
   *  <td> ExprType </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  </table>
   *  
   *  While a StandardTerm can be either a FieldExpression or a FcnForm, StandardizeArg standardizes the potential arguments to FcnForm objects.
   *  Since FcnForm objects can only take FieldExpression objects as actual arguments, StandardizeArg only produces FieldExpression objects.
   *  
   *  StandardizeArg should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename ExprType, typename FieldType>
    struct StandardizeArg<FieldExpression<ExprType,FieldType>,FieldType> {

    /**
     *  @brief Typedef for StandardType.
     *  
     *  Typedef for StandardType, here ExprType.
     */
    ExprType typedef StandardType;
    
    /**
     *  @brief Typedef for StandardTerm.
     *  
     *  Typedef for StandardTerm, here FieldExpression<StandardType, ...>.
     */
    FieldExpression<StandardType,FieldType> typedef StandardTerm;
    
    /**
     *  @brief Converts to StandardType.
     *  
     *  \param given A FieldExpression object.
     *  \return a new StandardType object.
     *  
     *  Builds and returns a StandardType object based on given.
     */
    static inline StandardType standardType (FieldExpression<ExprType,FieldType> const & given) {
      return given.expression();
    };
    
    /**
     *  @brief Converts to StandardTerm.
     *  
     *  \param given A FieldExpression object.
     *  \return a StandardTerm object (here a FieldExpression<ExprType, ...> object).
     *  
     *  Converts and returns a StandardTerm object based on given, (here a FieldExpression<ExprType, ...> object).
     */
    static inline StandardTerm standardTerm (FieldExpression<ExprType,FieldType> const & given) {
      return given;
    };
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the ArgReplace type conversions (type juggling) for applying arguments to a FcnForm. (Scalar case.)
   *  
   *  \tparam BeginType A FieldExpression-style type, which is the representation of an expression/anonymous function; here, Scalar<AtomicType>.
   *  \tparam CurrentArg An integer, representing the current argument to replace.
   *  \tparam ArgType A FieldExpression-style type, which is the actual argument to replace any ArgForm<CurrentArg, ...>.
   *  
   *  ArgReplace replaces ArgForm<CurrentArg, ...> with ArgType in BeginType.
   *  ArgReplace recurses down through the syntax tree/FieldExpression-style representation of an expression to find all applicable anonymous arguments.
   *  
   *  ArgReplace preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> BeginType </th>
   *  <th> Result Type </th>
   *  </tr>
   *  <tr>
   *  <td> Scalar<AtomicType> </td>
   *  <td> Scalar<AtomicType> </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> FieldForm<FieldType> </td>
   *  <td> FieldForm<FieldType> </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<CurrentArg, ...> </td>
   *  <td> ArgType </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> BinOp<Operand1, Operand2, ...> </td>
   *  <td> BinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificBinOp<Operand1, Operand2, ...> </td>
   *  <td> SpecificBinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> UnOp<Operand, ...> </td>
   *  <td> UnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificUnOp<Operand, ...> </td>
   *  <td> SpecificUnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  </table>
   *  
   *  ArgReplace should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<int CurrentArg, typename ArgType, typename AtomicType>
    struct ArgReplace<Scalar<AtomicType>,CurrentArg,ArgType> {
    
    /**
     *  @brief Typedef for ResultType.
     *  
     *  Typedef for ResultType, here Scalar<AtomicType>.
     */
    Scalar<AtomicType> typedef ResultType;
    
    /**
     *  @brief Converts to ResultType.
     *  
     *  \param state A Scalar<AtomicType> object.
     *  \param arg An ArgType object.
     *  \return a StandardType object (state).
     *  
     *  Returns state.
     */
    static inline ResultType const & apply(Scalar<AtomicType> const & state,
					   ArgType const & arg) {
      return state;
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the ArgReplace type conversions (type juggling) for applying arguments to a FcnForm. (FieldForm case.)
   *  
   *  \tparam BeginType A FieldExpression-style type, which is the representation of an expression/anonymous function; here, FieldForm<FieldType>.
   *  \tparam CurrentArg An integer, representing the current argument to replace.
   *  \tparam ArgType A FieldExpression-style type, which is the actual argument to replace any ArgForm<CurrentArg, ...>.
   *  
   *  ArgReplace replaces ArgForm<CurrentArg, ...> with ArgType in BeginType.
   *  ArgReplace recurses down through the syntax tree/FieldExpression-style representation of an expression to find all applicable anonymous arguments.
   *  
   *  ArgReplace preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> BeginType </th>
   *  <th> Result Type </th>
   *  </tr>
   *  <tr>
   *  <td> Scalar<AtomicType> </td>
   *  <td> Scalar<AtomicType> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldForm<FieldType> </td>
   *  <td> FieldForm<FieldType> </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<CurrentArg, ...> </td>
   *  <td> ArgType </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> BinOp<Operand1, Operand2, ...> </td>
   *  <td> BinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificBinOp<Operand1, Operand2, ...> </td>
   *  <td> SpecificBinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> UnOp<Operand, ...> </td>
   *  <td> UnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificUnOp<Operand, ...> </td>
   *  <td> SpecificUnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  </table>
   *  
   *  ArgReplace should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<int CurrentArg, typename ArgType, typename FieldType>
    struct ArgReplace<FieldForm<FieldType>,CurrentArg,ArgType> {
    
    /**
     *  @brief Typedef for ResultType.
     *  
     *  Typedef for ResultType, here FieldForm<FieldType>.
     */
    FieldForm<FieldType> typedef ResultType;
    
    /**
     *  @brief Converts to ResultType.
     *  
     *  \param state A FieldForm<FieldType> object.
     *  \param arg An ArgType object.
     *  \return a StandardType object (state).
     *  
     *  Returns state.
     */
    static inline ResultType const & apply(FieldForm<FieldType> const & state,
					   ArgType const & arg) {
      return state;
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the ArgReplace type conversions (type juggling) for applying arguments to a FcnForm. (Matching ArgForm case.)
   *  
   *  \tparam BeginType A FieldExpression-style type, which is the representation of an expression/anonymous function; here, ArgForm<CurrentArg, ...>.
   *  \tparam CurrentArg An integer, representing the current argument to replace.
   *  \tparam ArgType A FieldExpression-style type, which is the actual argument to replace any ArgForm<CurrentArg, ..>.
   *  
   *  ArgReplace replaces ArgForm<CurrentArg, ...> with ArgType in BeginType.
   *  ArgReplace recurses down through the syntax tree/FieldExpression-style representation of an expression to find all applicable anonymous arguments.
   *  
   *  ArgReplace preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> BeginType </th>
   *  <th> Result Type </th>
   *  </tr>
   *  <tr>
   *  <td> Scalar<AtomicType> </td>
   *  <td> Scalar<AtomicType> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldForm<FieldType> </td>
   *  <td> FieldForm<FieldType> </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<CurrentArg, ...> </td>
   *  <td> ArgType </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> BinOp<Operand1, Operand2, ...> </td>
   *  <td> BinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificBinOp<Operand1, Operand2, ...> </td>
   *  <td> SpecificBinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> UnOp<Operand, ...> </td>
   *  <td> UnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificUnOp<Operand, ...> </td>
   *  <td> SpecificUnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  </table>
   *  
   *  ArgReplace should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<int CurrentArg, typename ArgType, typename FieldType>
    struct ArgReplace<ArgForm<CurrentArg,FieldType>,CurrentArg,ArgType> {
    
    /**
     *  @brief Typedef for ResultType.
     *  
     *  Typedef for ResultType, here ArgType.
     */
    ArgType typedef ResultType;
    
    /**
     *  @brief Converts to ResultType.
     *  
     *  \param state An ArgForm<CurrentArg, ...> object.
     *  \param arg An ArgType object.
     *  \return an ArgType object (arg).
     *  
     *  Returns arg.
     */
    static inline ResultType const & apply(ArgForm<CurrentArg,FieldType> const & state,
					   ArgType const & arg) {
      return arg;
    };
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the ArgReplace type conversions (type juggling) for applying arguments to a FcnForm. (Non-matching ArgForm case.)
   *  
   *  \tparam BeginType A FieldExpression-style type, which is the representation of an expression/anonymous function; here, ArgForm<NotCurrentArg, ...>.
   *  \tparam CurrentArg An integer, representing the current argument to replace.
   *  \tparam ArgType A FieldExpression-style type, which is the actual argument to replace any ArgForm<CurrentArg, ..>.
   *  
   *  ArgReplace replaces ArgForm<CurrentArg, ...> with ArgType in BeginType.
   *  ArgReplace recurses down through the syntax tree/FieldExpression-style representation of an expression to find all applicable anonymous arguments.
   *  
   *  ArgReplace preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> BeginType </th>
   *  <th> Result Type </th>
   *  </tr>
   *  <tr>
   *  <td> Scalar<AtomicType> </td>
   *  <td> Scalar<AtomicType> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldForm<FieldType> </td>
   *  <td> FieldForm<FieldType> </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<CurrentArg, ...> </td>
   *  <td> ArgType </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> BinOp<Operand1, Operand2, ...> </td>
   *  <td> BinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificBinOp<Operand1, Operand2, ...> </td>
   *  <td> SpecificBinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> UnOp<Operand, ...> </td>
   *  <td> UnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificUnOp<Operand, ...> </td>
   *  <td> SpecificUnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  </table>
   *  
   *  ArgReplace should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<int NotCurrentArg, int CurrentArg, typename ArgType, typename FieldType>
    struct ArgReplace<ArgForm<NotCurrentArg,FieldType>,CurrentArg,ArgType> {
    
    /**
     *  @brief Typedef for ResultType.
     *  
     *  Typedef for ResultType, here ArgForm<NotCurrentArg, ...>.
     */
    ArgForm<NotCurrentArg,FieldType> typedef ResultType;
    
    /**
     *  @brief Converts to ResultType.
     *  
     *  \param state An ArgForm<NotCurrentArg, ...> object.
     *  \param arg An ArgType object.
     *  \return an ArgForm<NotCurrentArg, ...> object (state).
     *  
     *  Returns state.
     */
    static inline ResultType const & apply(ArgForm<NotCurrentArg,FieldType> const & state,
					   ArgType const & arg) {
      return state;
    };
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the ArgReplace type conversions (type juggling) for applying arguments to a FcnForm. (BinOp case.)
   *  
   *  \tparam BeginType A FieldExpression-style type, which is the representation of an expression/anonymous function; here, BinOp<Operand1, Operand2, ...>.
   *  \tparam CurrentArg An integer, representing the current argument to replace.
   *  \tparam ArgType A FieldExpression-style type, which is the actual argument to replace any ArgForm<CurrentArg, ..>.
   *  
   *  ArgReplace replaces ArgForm<CurrentArg, ...> with ArgType in BeginType.
   *  ArgReplace recurses down through the syntax tree/FieldExpression-style representation of an expression to find all applicable anonymous arguments.
   *  
   *  ArgReplace preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> BeginType </th>
   *  <th> Result Type </th>
   *  </tr>
   *  <tr>
   *  <td> Scalar<AtomicType> </td>
   *  <td> Scalar<AtomicType> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldForm<FieldType> </td>
   *  <td> FieldForm<FieldType> </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<CurrentArg, ...> </td>
   *  <td> ArgType </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> BinOp<Operand1, Operand2, ...> </td>
   *  <td> BinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificBinOp<Operand1, Operand2, ...> </td>
   *  <td> SpecificBinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> UnOp<Operand, ...> </td>
   *  <td> UnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificUnOp<Operand, ...> </td>
   *  <td> SpecificUnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  </table>
   *  
   *  ArgReplace should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename Operand1, typename Operand2, int CurrentArg, typename ArgType, typename FieldType> 
    struct ArgReplace<BinOp<Operand1,Operand2,FieldType>,CurrentArg,ArgType> {
    
    /**
     *  @brief Typedef for ResultType.
     *  
     *  Typedef for ResultType, here BinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...>.
     */
    BinOp<typename ArgReplace<Operand1,CurrentArg,ArgType>::ResultType,
      typename ArgReplace<Operand2,CurrentArg,ArgType>::ResultType,
      FieldType> typedef ResultType;
    
    /**
     *  @brief Converts to ResultType.
     *  
     *  \param state A BinOp<Operand1, Operand2, ...> object.
     *  \param arg An ArgType object.
     *  \return a new BinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> object.
     *  
     *  Builds and returns a new ResultType object based upon the ResultTypes of the ArgReplace'ed operand types.
     */
    static inline ResultType apply(BinOp<Operand1,Operand2,FieldType> const & state,
				   ArgType const & arg) {
      return ResultType(state.fcn(),
			ArgReplace<Operand1,CurrentArg,ArgType>::apply(state.first(),
								     arg),
			ArgReplace<Operand2,CurrentArg,ArgType>::apply(state.second(),
								     arg));
    };
  };

  /* No doxygen comments for the code generated by this macro.  This is a design-limitation of doxygen. */
  /* (For marco-generated types, doxygen documentation must refer to the type by name.) */
  /* (For specialized templates, doxygen documentation must refer to the type by proximity to the specialization.) */
  /* (This macro defines a specialization, so these two requirements conflict.) */
#define BUILD_BINARY_ARGUMENT_REPLACE(OBJECT_NAME)			\
  template <typename Operand1, typename Operand2, int CurrentArg, typename ArgType, typename FieldType> \
    struct ArgReplace<OBJECT_NAME<Operand1,Operand2,FieldType>,CurrentArg,ArgType> { \
    /* ResultType is  SpecificBinOp<ArgReplace<Operand1, ...>::ResultType, */ \
    /*                              ArgReplace<Operand2, ...>::ResultType, ...> */ \
    OBJECT_NAME<typename ArgReplace<Operand1,CurrentArg,ArgType>::ResultType, \
      typename ArgReplace<Operand2,CurrentArg,ArgType>::ResultType,	\
      FieldType> typedef ResultType;					\
    									\
    /* Returns ResultType build from state's operands. */		\
    static inline ResultType apply(OBJECT_NAME<Operand1,Operand2,FieldType> const & state, \
				   ArgType const & arg) {		\
      return ResultType(ArgReplace<Operand1,CurrentArg,ArgType>::apply(state.first(), \
								     arg), \
			ArgReplace<Operand2,CurrentArg,ArgType>::apply(state.second(), \
								     arg)); \
    };									\
  }
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the ArgReplace type conversions (type juggling) for applying arguments to a FcnForm. (UnOp case.)
   *  
   *  \tparam BeginType A FieldExpression-style type, which is the representation of an expression/anonymous function; here, UnOp<Operand, ...>.
   *  \tparam CurrentArg An integer, representing the current argument to replace.
   *  \tparam ArgType A FieldExpression-style type, which is the actual argument to replace any ArgForm<CurrentArg, ..>.
   *  
   *  ArgReplace replaces ArgForm<CurrentArg, ...> with ArgType in BeginType.
   *  ArgReplace recurses down through the syntax tree/FieldExpression-style representation of an expression to find all applicable anonymous arguments.
   *  
   *  ArgReplace preforms the following conversions:
   *  
   *  <table>
   *  <tr>
   *  <th> BeginType </th>
   *  <th> Result Type </th>
   *  </tr>
   *  <tr>
   *  <td> Scalar<AtomicType> </td>
   *  <td> Scalar<AtomicType> </td>
   *  </tr>
   *  <tr>
   *  <td> FieldForm<FieldType> </td>
   *  <td> FieldForm<FieldType> </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<CurrentArg, ...> </td>
   *  <td> ArgType </td>
   *  </tr>
   *  <tr>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  <td> ArgForm<NotCurrentArg, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> BinOp<Operand1, Operand2, ...> </td>
   *  <td> BinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificBinOp<Operand1, Operand2, ...> </td>
   *  <td> SpecificBinOp<ArgReplace<Operand1, ...>::ResultType, ArgReplace<Operand2, ...>::ResultType, ...> </td>
   *  </tr>
   *  <tr>
   *  <td> UnOp<Operand, ...> </td>
   *  <td> UnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  <td> (Current case.) </td>
   *  </tr>
   *  <tr>
   *  <td> SpecificUnOp<Operand, ...> </td>
   *  <td> SpecificUnOp<ArgReplace<Operand, ...>::ResultType, ...> </td>
   *  </tr>
   *  </table>
   *  
   *  ArgReplace should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename Operand, int CurrentArg, typename ArgType, typename FieldType>
    struct ArgReplace<UnOp<Operand,FieldType>,CurrentArg,ArgType> {
    
    /**
     *  @brief Typedef for ResultType.
     *  
     *  Typedef for ResultType, here UnOp<ArgReplace<Operand, ...>::ResultType, ...>.
     */
    UnOp<typename ArgReplace<Operand,CurrentArg,ArgType>::ResultType,FieldType> typedef ResultType;
    
    /**
     *  @brief Converts to ResultType.
     *  
     *  \param state A UnOp<Operand, ...> object.
     *  \param arg An ArgType object.
     *  \return a new UnOp<ArgReplace<Operand, ...>::ResultType, ...> object.
     *  
     *  Builds and returns a new ResultType object based upon the ResultTypes of the ArgReplace'ed Operand type.
     */
    static inline ResultType apply(UnOp<Operand,FieldType> const & state,
				   ArgType const & arg) {
      return ResultType(state.op,
			ArgReplace<Operand,CurrentArg,ArgType>::apply(state.operand,
								    arg));
    };
  };

  /* No doxygen comments for the code generated by this macro.  This is a design-limitation of doxygen. */
  /* (For marco-generated types, doxygen documentation must refer to the type by name.) */
  /* (For specialized templates, doxygen documentation must refer to the type by proximity to the specialization.) */
  /* (This macro defines a specialization, so these two requirements conflict.) */
#define BUILD_UNARY_ARGUMENT_REPLACE(OBJECT_NAME)			\
  template <typename Operand, int CurrentArg, typename ArgType, typename FieldType> \
    struct ArgReplace<OBJECT_NAME<Operand,FieldType>,CurrentArg,ArgType> { \
    /* ResultType is  SpecificBinOp<ArgReplace<Operand1, ...>::ResultType, ...> */ \
    OBJECT_NAME<typename ArgReplace<Operand,CurrentArg,ArgType>::ResultType,FieldType> typedef ResultType; \
									\
    /* Returns ResultType build from state's operand. */		\
    static inline ResultType apply(OBJECT_NAME<Operand,FieldType> const & state, \
				   ArgType const & arg) {		\
      return ResultType(ArgReplace<Operand,CurrentArg,ArgType>::apply(state.operand, \
								    arg)); \
    };									\
  }
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the ArgResultTerm type conversions (type juggling) for the result of applying arguments to a FcnForm. (All argumets applied case.)
   *  
   *  \tparam StandardType A FieldExpression-style type, which is the representation of an expression/anonymous function.
   *  \tparam CrtNum An integer, representing the current argument to replace.
   *  \tparam Max A type-represenstation for the maximum/final argument number (Int<Num>, for some integer Num; here, Num = CrtNum + 1).
   *  \tparam FieldType Field type.
   *  
   *  ArgResultTerm examines the current argument number (CrtNum) with the total (maximum) number of arguments.
   *  Arguments are indexed from 0, and the total number of arguments is not indexed from 0.
   *  Thus, if CrtNum is the last argument, then the representation of the total number of arguments (Max) is Int<CrtNum + 1>, as here.
   *  Otherwise, if CrtNum is not the last argument, then there is no (constant) relationship between CrtNum and Max.
   *
   *  If there are no more arguments to apply (i.e. the last argument has just been applied), then StandardType contains no more anonyous arguments and is no longer an anonymous function.
   *  Thus, if the last argument has just been applied, as here, the StandardTerm for StandardType is a FieldExpression.
   *  In this case, ResultType is FieldExpression<StandardType, FieldType>, as here.
   *  
   *  Otherwise, there are arguments still to apply, then there are anonymous arguments left in StandardType, and StandardType is still an anonymous function.
   *  Thus, if there are still anonymous arguments, the StandardTerm for StandardType is a FcnForm.
   *  In this case, ResultType is FcnForm<StandardType, CrtNum + 1, Int<MaxNum>, FieldType>.
   *  
   *  ArgResultTerm should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename StandardType, int CrtNum, typename FieldType>
    struct ArgResultTerm<StandardType,CrtNum,Int<CrtNum + 1>,FieldType> {
    
    /**
     *  @brief Typedef for ResultType.
     *  
     *  Typedef for ResultType, here FieldExpression<ArgReplace<StandardType, ...>::ResultType, ...>.
     */
    FieldExpression<StandardType,FieldType> typedef ResultType;
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Structure encoding the ArgResultTerm type conversions (type juggling) for the result of applying arguments to a FcnForm. (Some argumets remain unapplied case.)
   *  
   *  \tparam StandardType A FieldExpression-style type, which is the representation of an expression/anonymous function.
   *  \tparam CrtNum An integer, representing the current argument to replace.
   *  \tparam Max A type-represenstation for the maximum/final argument number (Int<MaxNum>, for some integer MaxNum; here, MaxNum with no relation to CrtNum).
   *  \tparam FieldType Field type.
   *  
   *  ArgResultTerm examines the current argument number (CrtNum) with the total (maximum) number of arguments.
   *  Arguments are indexed from 0, and the total number of arguments is not indexed from 0.
   *  Thus, if CrtNum is the last argument, then the representation of the total number of arguments (Max) is Int<CrtNum + 1>.
   *  Otherwise, if CrtNum is not the last argument, then there is no (constant) relationship between CrtNum and Max, as here.
   *
   *  If there are no more arguments to apply (i.e. the last argument has just been applied), then StandardType contains no more anonyous arguments and is no longer an anonymous function.
   *  Thus, if the last argument has just been applied, the StandardTerm for StandardType is a FieldExpression.
   *  In this case, ResultType is FieldExpression<StandardType, FieldType>.
   *  
   *  Otherwise, there are arguments still to apply, then there are anonymous arguments left in StandardType, and StandardType is still an anonymous function.
   *  Thus, if there are still anonymous arguments, as here, the StandardTerm for StandardType is a FcnForm.
   *  In this case, ResultType is FcnForm<StandardType, CrtNum + 1, Int<MaxNum>, FieldType>, as here.
   *  
   *  ArgResultTerm should never be instantiated; its sole purpose is type conversion/juggling.
   */
  template<typename StandardType, int CrtNum, int MaxNum, typename FieldType>
    struct ArgResultTerm<StandardType,CrtNum,Int<MaxNum>,FieldType> {
    
    /**
     *  @brief Typedef for ResultType.
     *  
     *  Typedef for ResultType, here FcnForm<ArgReplace<StandardType, ...>::ResultType, ...>.
     */
    FcnForm<StandardType,
      CrtNum + 1,
      Int<MaxNum>,
      FieldType> typedef ResultType;
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Definition of app for arguments SubExpr and SubExpr.
   *  
   *  @relates BinOp
   *  
   *  \tparam SubExpr1 First operand's type.
   *  \tparam SubExpr2 Second operand's type.
   *
   *  \param first A SubExpr1 object.
   *  \param second A SubExpr2 object.
   *
   *  \return the Standard Term that is the result of creating a BinOp with first and second as its operands.
   *  
   *  A SubExpr can be a FieldType, a FieldExpression, an ArgForm, or an FcnForm.
   *  Basically, a SubExpr can be any FieldExpression except for Scalar<AtomicType>.
   *  Scalar<AtomicType> is a special argument because it does not define field_type.
   *  (Every SubExpr defines field_type.
   *   Since Scalar<AtomicType> does not, it cannot be used as a SubExpr.)
   *  Scalar<AtomicType> as an argument is treated by special cases (not here).
   *  
   *  All app does is build the correct StandardTerm type (and the correct StandardTerm object) for a BinOp object with the given arguments as operands.
   *  Since a SubExpr can be any one of many different types, StandardizeTerm is used to find and build the correct types and objects.
   */
  template<typename SubExpr1, typename SubExpr2>
    typename CombineTerms<BinOp<typename StandardizeTerm<SubExpr1,
    typename SubExpr1::field_type>::StandardType,
    typename StandardizeTerm<SubExpr2,
    typename SubExpr1::field_type>::StandardType,
    typename SubExpr1::field_type>,
    typename StandardizeTerm<SubExpr1,
    typename SubExpr1::field_type>::StandardTerm,
    typename StandardizeTerm<SubExpr2,
    typename SubExpr1::field_type>::StandardTerm,
    typename SubExpr1::field_type>::StandardTerm
    app (typename SubExpr1::field_type::value_type (*fcn)(typename SubExpr1::field_type::value_type,
							  typename SubExpr1::field_type::value_type),
	 SubExpr1 const & first,
	 SubExpr2 const & second) {
  
    typename SubExpr1::field_type typedef FieldType;
  
    typename StandardizeTerm<SubExpr1,FieldType>::StandardTerm typedef Term1;
    typename StandardizeTerm<SubExpr2,FieldType>::StandardTerm typedef Term2;
  
    typename StandardizeTerm<SubExpr1,FieldType>::StandardType typedef Type1;
    typename StandardizeTerm<SubExpr2,FieldType>::StandardType typedef Type2;
  
    BinOp<Type1,Type2,FieldType> typedef ReturnType;
    typename CombineTerms<ReturnType,
      Term1,
      Term2,
      FieldType>::StandardTerm typedef ReturnTerm;
  
    return ReturnTerm(ReturnType(fcn,
				 StandardizeTerm<SubExpr1,FieldType>::standardType(first),
				 StandardizeTerm<SubExpr2,FieldType>::standardType(second)));
  };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Definition of app for arguments SubExpr and Scalar.
   *  
   *  @relates BinOp
   *  
   *  \tparam SubExpr First operand's type.
   *
   *  \param first A SubExpr object.
   *  \param second A Scalar object.
   *
   *  \return the Standard Term that is the result of creating a BinOp with first and second as its operands.
   *  
   *  A SubExpr can be a FieldType, a FieldExpression, an ArgForm, or an FcnForm.
   *  Basically, a SubExpr can be any FieldExpression except for Scalar<AtomicType>.
   *  Scalar<AtomicType> is a special argument because it does not define field_type.
   *  (Every SubExpr defines field_type.
   *   Since Scalar<AtomicType> does not, it cannot be used as a SubExpr.)
   *  Scalar<AtomicType> as an argument is treated by special cases, as here.
   *  
   *  All app does is build the correct StandardTerm type (and the correct StandardTerm object) for a BinOp object with the given arguments as operands.
   *  Since a SubExpr can be any one of many different types, StandardizeTerm is used to find and build the correct types and objects.
   */
  template<typename SubExpr>
    typename CombineTerms<BinOp<typename StandardizeTerm<SubExpr,
    typename SubExpr::field_type>::StandardType,
    Scalar<typename SubExpr::field_type::value_type>,
    typename SubExpr::field_type>,
    typename StandardizeTerm<SubExpr,
    typename SubExpr::field_type>::StandardTerm,
    FieldExpression<Scalar<typename SubExpr::field_type::value_type>,
    typename SubExpr::field_type>,
    typename SubExpr::field_type>::StandardTerm
    app (typename SubExpr::field_type::value_type (*fcn)(typename SubExpr::field_type::value_type,
							 typename SubExpr::field_type::value_type),
	 SubExpr const & first,
	 typename SubExpr::field_type::value_type const & second) {
  
    typename SubExpr::field_type typedef FieldType;
    typename FieldType::value_type typedef AtomicType;
  
    typename StandardizeTerm<SubExpr,FieldType>::StandardTerm typedef Term1;
    FieldExpression<Scalar<AtomicType>,FieldType> typedef Term2;
  
    typename StandardizeTerm<SubExpr,FieldType>::StandardType typedef Type1;
    Scalar<AtomicType> typedef Type2;
  
    BinOp<Type1,Type2,FieldType> typedef ReturnType;
    typename CombineTerms<ReturnType,
      Term1,
      Term2,
      FieldType>::StandardTerm typedef ReturnTerm;
  
    return ReturnTerm(ReturnType(fcn,
				 StandardizeTerm<SubExpr,FieldType>::standardType(first),
				 Type2(second)));
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Definition of app for arguments Scalar and SubExpr.
   *  
   *  @relates BinOp
   *  
   *  \tparam SubExpr Second operand's type.
   *
   *  \param second A Scalar object.
   *  \param first A SubExpr1 object.
   *
   *  \return the Standard Term that is the result of creating a BinOp with first and second as its operands.
   *  
   *  A SubExpr can be a FieldType, a FieldExpression, an ArgForm, or an FcnForm.
   *  Basically, a SubExpr can be any FieldExpression except for Scalar<AtomicType>.
   *  Scalar<AtomicType> is a special argument because it does not define field_type.
   *  (Every SubExpr defines field_type.
   *   Since Scalar<AtomicType> does not, it cannot be used as a SubExpr.)
   *  Scalar<AtomicType> as an argument is treated by special cases, as here.
   *  
   *  All app does is build the correct StandardTerm type (and the correct StandardTerm object) for a BinOp object with the given arguments as operands.
   *  Since a SubExpr can be any one of many different types, StandardizeTerm is used to find and build the correct types and objects.
   */
  template<typename SubExpr>
    typename CombineTerms<BinOp<Scalar<typename SubExpr::field_type::value_type>,
    typename StandardizeTerm<SubExpr,
    typename SubExpr::field_type>::StandardType,
    typename SubExpr::field_type>,
    FieldExpression<Scalar<typename SubExpr::field_type::value_type>,
    typename SubExpr::field_type>,
    typename StandardizeTerm<SubExpr,
    typename SubExpr::field_type>::StandardTerm,
    typename SubExpr::field_type>::StandardTerm
    app (typename SubExpr::field_type::value_type (*fcn)(typename SubExpr::field_type::value_type,
							 typename SubExpr::field_type::value_type),
	 typename SubExpr::field_type::value_type const & first,
	 SubExpr const & second) {
  
    typename SubExpr::field_type typedef FieldType;
    typename FieldType::value_type typedef AtomicType;
  
    FieldExpression<Scalar<AtomicType>,FieldType> typedef Term1;
    typename StandardizeTerm<SubExpr,FieldType>::StandardTerm typedef Term2;
  
    Scalar<AtomicType> typedef Type1;
    typename StandardizeTerm<SubExpr,FieldType>::StandardType typedef Type2;
  
    BinOp<Type1,Type2,FieldType> typedef ReturnType;
    typename CombineTerms<ReturnType,
      Term1,
      Term2,
      FieldType>::StandardTerm typedef ReturnTerm;
  
    return ReturnTerm(ReturnType(fcn,
				 Type1(first),
				 StandardizeTerm<SubExpr,FieldType>::standardType(second)));
  };

  /* No doxygen comments for the code generated by this macro.  This is a design-limitation of doxygen. */
  /* (For marco-generated types/functions, doxygen documentation must refer to the type/function by name.) */
  /* (For specialized templates, doxygen documentation must refer to the type/function by proximity to the specialization.) */
  /* (This macro defines a specialization, so these two requirements conflict.) */
#define BUILD_BINARY_INTERFACE(OBJECT_NAME, EXTERNAL_NAME)		\
  /* SubExpr X SubExpr: */						\
    template<typename SubExpr1, typename SubExpr2>			\
      typename CombineTerms<OBJECT_NAME<typename StandardizeTerm<SubExpr1, \
      typename SubExpr1::field_type>::StandardType,			\
      typename StandardizeTerm<SubExpr2,				\
      typename SubExpr1::field_type>::StandardType,			\
      typename SubExpr1::field_type>,					\
      typename StandardizeTerm<SubExpr1,				\
      typename SubExpr1::field_type>::StandardTerm,			\
      typename StandardizeTerm<SubExpr2,				\
      typename SubExpr1::field_type>::StandardTerm,			\
      typename SubExpr1::field_type>::StandardTerm			\
      EXTERNAL_NAME (SubExpr1 const & first,				\
		     SubExpr2 const & second) {				\
									\
      typename SubExpr1::field_type typedef FieldType;			\
									\
      typename StandardizeTerm<SubExpr1,FieldType>::StandardTerm typedef Term1; \
      typename StandardizeTerm<SubExpr2,FieldType>::StandardTerm typedef Term2; \
									\
      typename StandardizeTerm<SubExpr1,FieldType>::StandardType typedef Type1; \
      typename StandardizeTerm<SubExpr2,FieldType>::StandardType typedef Type2; \
									\
      OBJECT_NAME<Type1,Type2,FieldType> typedef ReturnType;		\
      typename CombineTerms<ReturnType,					\
	Term1,								\
	Term2,								\
	FieldType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(StandardizeTerm<SubExpr1,FieldType>::standardType(first), \
				   StandardizeTerm<SubExpr2,FieldType>::standardType(second))); \
    };									\
  									\
    /* SubExpr X Scalar: */						\
    template<typename SubExpr>						\
      typename CombineTerms<OBJECT_NAME<typename StandardizeTerm<SubExpr, \
      typename SubExpr::field_type>::StandardType,			\
      Scalar<typename SubExpr::field_type::value_type>,			\
      typename SubExpr::field_type>,					\
      typename StandardizeTerm<SubExpr,					\
      typename SubExpr::field_type>::StandardTerm,			\
      FieldExpression<Scalar<typename SubExpr::field_type::value_type>,	\
      typename SubExpr::field_type>,					\
      typename SubExpr::field_type>::StandardTerm			\
      EXTERNAL_NAME (SubExpr const & first,				\
		     typename SubExpr::field_type::value_type const & second) { \
      									\
      typename SubExpr::field_type typedef FieldType;			\
      typename FieldType::value_type typedef AtomicType;		\
									\
      typename StandardizeTerm<SubExpr,FieldType>::StandardTerm typedef Term1; \
      FieldExpression<Scalar<AtomicType>,FieldType> typedef Term2;		\
									\
      typename StandardizeTerm<SubExpr,FieldType>::StandardType typedef Type1; \
      Scalar<AtomicType> typedef Type2;					\
    									\
      OBJECT_NAME<Type1,Type2,FieldType> typedef ReturnType;		\
      typename CombineTerms<ReturnType,					\
	Term1,								\
	Term2,								\
	FieldType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(StandardizeTerm<SubExpr,FieldType>::standardType(first), \
				   Type2(second)));			\
    };									\
									\
    /* Scalar X SubExpr: */						\
    template<typename SubExpr>						\
      typename CombineTerms<OBJECT_NAME<Scalar<typename SubExpr::field_type::value_type>, \
      typename StandardizeTerm<SubExpr,					\
      typename SubExpr::field_type>::StandardType,			\
      typename SubExpr::field_type>,					\
      FieldExpression<Scalar<typename SubExpr::field_type::value_type>,	\
      typename SubExpr::field_type>,					\
      typename StandardizeTerm<SubExpr,					\
      typename SubExpr::field_type>::StandardTerm,			\
      typename SubExpr::field_type>::StandardTerm			\
      EXTERNAL_NAME (typename SubExpr::field_type::value_type const & first, \
		     SubExpr const & second) {				\
    									\
      typename SubExpr::field_type typedef FieldType;			\
      typename FieldType::value_type typedef AtomicType;		\
									\
      FieldExpression<Scalar<AtomicType>,FieldType> typedef Term1;		\
      typename StandardizeTerm<SubExpr,FieldType>::StandardTerm typedef Term2; \
									\
      Scalar<AtomicType> typedef Type1;					\
      typename StandardizeTerm<SubExpr,FieldType>::StandardType typedef Type2; \
									\
      OBJECT_NAME<Type1,Type2,FieldType> typedef ReturnType;		\
      typename CombineTerms<ReturnType,					\
	Term1,								\
	Term2,								\
	FieldType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(Type1(first),			\
				   StandardizeTerm<SubExpr,FieldType>::standardType(second))); \
    }
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Definition of app for argument SubExpr.
   *  
   *  @relates UnOp
   *  
   *  \tparam SubExpr Operand's type.
   *  
   *  \param argument A SubExpr1 object.
   *
   *  \return the Standard Term that is the result of creating an UnOp with argument as its operand.
   *  
   *  A SubExpr can be a FieldType, a FieldExpression, an ArgForm, or an FcnForm.
   *  Basically, a SubExpr can be any FieldExpression except for Scalar<AtomicType>.
   *  Scalar<AtomicType> is a special argument because it does not define field_type.
   *  (Every SubExpr defines field_type.
   *   Since Scalar<AtomicType> does not, it cannot be used as a SubExpr.)
   *  Scalar<AtomicType> cannot be used as an argument type.
   *  (However, where \c app(fcn, \c value) is desired, \c fcn(value) will work.)
   *  
   *  All app does is build the correct StandardTerm type (and the correct StandardTerm object) for an UnOp object with argument as operand.
   *  Since a SubExpr can be any one of many different types, StandardizeTerm is used to find and build the correct types and objects.
   */
  template<typename SubExpr>
    typename LiftTerm<UnOp<typename StandardizeTerm<SubExpr,
    typename SubExpr::field_type>::StandardType,
    typename SubExpr::field_type>,
    typename StandardizeTerm<SubExpr,
    typename SubExpr::field_type>::StandardTerm,
    typename SubExpr::field_type>::StandardTerm
    app (typename SubExpr::field_type::value_type (*fcn)(typename SubExpr::field_type::value_type),
	 SubExpr const & argument) {
  
    typename SubExpr::field_type typedef FieldType;
  
    typename StandardizeTerm<SubExpr,FieldType>::StandardTerm typedef Term;
  
    typename StandardizeTerm<SubExpr,FieldType>::StandardType typedef Type;
  
    UnOp<Type,FieldType> typedef ReturnType;
    typename LiftTerm<ReturnType,
      Term,
      FieldType>::StandardTerm typedef ReturnTerm;
  
    return ReturnTerm(ReturnType(fcn,
				 StandardizeTerm<SubExpr,FieldType>::standardType(argument)));
  };

  /* No doxygen comments for the code generated by this macro.  This is a design-limitation of doxygen. */
  /* (For marco-generated types/functions, doxygen documentation must refer to the type/function by name.) */
  /* (For specialized templates, doxygen documentation must refer to the type/function by proximity to the specialization.) */
  /* (This macro defines a specialization, so these two requirements conflict.) */
#define BUILD_UNARY_INTERFACE(OBJECT_NAME, EXTERNAL_NAME)		\
  /* SubExpr: */							\
    template<typename SubExpr>						\
      typename LiftTerm<OBJECT_NAME<typename StandardizeTerm<SubExpr,	\
      typename SubExpr::field_type>::StandardType,			\
      typename SubExpr::field_type>,					\
      typename StandardizeTerm<SubExpr,					\
      typename SubExpr::field_type>::StandardTerm,			\
      typename SubExpr::field_type>::StandardTerm			\
      EXTERNAL_NAME (SubExpr const & argument) {			\
    									\
      typename SubExpr::field_type typedef FieldType;			\
									\
      typename StandardizeTerm<SubExpr,FieldType>::StandardTerm typedef Term; \
									\
      typename StandardizeTerm<SubExpr,FieldType>::StandardType typedef Type; \
									\
      OBJECT_NAME<Type,FieldType> typedef ReturnType;			\
      typename LiftTerm<ReturnType,					\
	Term,								\
	FieldType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(StandardizeTerm<SubExpr,FieldType>::standardType(argument))); \
    }

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Definition of operator <<= (assignment) for argument FieldExpression. (AtomicType case.)
   *  
   *  @relates FieldExpression
   *  
   *  \tparam FieldType Field type.
   *  
   *  \param lhs A FieldType object.
   *  \param rhs A AtomicType (FieldType::value_type) object.
   *
   *  \return a reference to the newly-assigned lhs.
   *  
   *  This instance of operator <<= wraps rhs in a FieldExpression and then calls itself.
   */
  template <typename FieldType>
    static inline FieldType const & operator <<= (FieldType & lhs,
						  typename FieldType::value_type const & rhs) {
    Scalar<typename FieldType::value_type> typedef ExprType;
    
    return lhs <<= FieldExpression<ExprType,FieldType>(ExprType(rhs));
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Definition of operator <<= (assignment) for argument FieldExpression. (FieldType case.)
   *  
   *  @relates FieldExpression
   *  
   *  \tparam FieldType Field type.
   *  
   *  \param lhs A FieldType object.
   *  \param rhs A FieldType object.
   *
   *  \return a reference to the newly-assigned lhs.
   *  
   *  This instance of operator <<= wraps rhs in a FieldExpression and then calls itself.
   */
  template <typename FieldType>
    static inline FieldType const & operator <<= (FieldType & lhs,
						  FieldType const & rhs) {
    FieldForm<FieldType> typedef ExprType;
    
    return lhs <<= FieldExpression<ExprType,FieldType>(ExprType(rhs));
  };
  
  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Definition of operator <<= (assignment) for argument FieldExpression. (FieldExpression case.)
   *  
   *  @relates FieldExpression
   *  
   *  \tparam ExprType A FieldExpression-style Standard Type.
   *  \tparam FieldType Field type.
   *  
   *  \param lhs A FieldType object.
   *  \param rhs A FieldExpression object.
   *
   *  \return a reference to the newly-assigned lhs.
   *  
   *  This instance of operator <<= does the actual assignment.
   *  The while loop that iterates over the elements in rhs is explicit here.
   *  The body of the loop evaluates the current element of the rhs and assigns it to the current element of the lhs.
   *  
   *  The evaluation of an element is via the eval from the ExprType object within rhs.
   *  
   *  \note There are no checks here to make sure that lhs and rhs are the same size.
   *  (The number of iterations of the loop is based on the size of lhs.)
   */
  template <typename ExprType, typename FieldType>
    static inline FieldType const & operator <<= (FieldType & lhs,
						  FieldExpression<ExprType,
						  FieldType> rhs) {
    //initialize:
    typename FieldType::iterator iter = lhs.begin();
    //initialize rhs:
    ExprType expr = rhs.expression();
    expr.init();
    
    while(iter != lhs.end()) {
      *iter = expr.eval();
      //increment:
      ++iter;
      expr.next();
    };
    
    return lhs;
  };

#define BUILD_BINARY_OPERATOR(OBJECT_NAME, INTERNAL_NAME, EXTERNAL_NAME) \
  BUILD_BINARY_TYPE_PROTOTYPE(OBJECT_NAME);				\
  BUILD_BINARY_OPERATOR_STRUCT(OBJECT_NAME, INTERNAL_NAME);		\
  BUILD_BINARY_ARGUMENT_REPLACE(OBJECT_NAME);				\
  BUILD_BINARY_INTERFACE(OBJECT_NAME, EXTERNAL_NAME)

#define BUILD_BINARY_FUNCTION(OBJECT_NAME, INTERNAL_NAME, EXTERNAL_NAME) \
  BUILD_BINARY_TYPE_PROTOTYPE(OBJECT_NAME);				\
  BUILD_BINARY_FUNCTION_STRUCT(OBJECT_NAME, INTERNAL_NAME);		\
  BUILD_BINARY_ARGUMENT_REPLACE(OBJECT_NAME);				\
  BUILD_BINARY_INTERFACE(OBJECT_NAME, EXTERNAL_NAME)

#define BUILD_UNARY_FUNCTION(OBJECT_NAME, INTERNAL_NAME, EXTERNAL_NAME) \
  BUILD_UNARY_TYPE_PROTOTYPE(OBJECT_NAME);				\
  BUILD_UNARY_STRUCT(OBJECT_NAME, INTERNAL_NAME);			\
  BUILD_UNARY_ARGUMENT_REPLACE(OBJECT_NAME);				\
  BUILD_UNARY_INTERFACE(OBJECT_NAME, EXTERNAL_NAME)
  
#define DEFINE_ANONYMOUS_ARGUMENT(NUMBER, TYPE, NAME)	\
  ArgForm<NUMBER,TYPE> NAME
  
#define DEFINE_SET_OF_ANONYMOUS_ARGUMENTS(PRE, POST, TYPE) \
  DEFINE_ANONYMOUS_ARGUMENT(0, TYPE, PRE##0##POST);	   \
  DEFINE_ANONYMOUS_ARGUMENT(1, TYPE, PRE##1##POST);	   \
  DEFINE_ANONYMOUS_ARGUMENT(2, TYPE, PRE##2##POST);	   \
  DEFINE_ANONYMOUS_ARGUMENT(3, TYPE, PRE##3##POST);	   \
  DEFINE_ANONYMOUS_ARGUMENT(4, TYPE, PRE##4##POST);	   \
  DEFINE_ANONYMOUS_ARGUMENT(5, TYPE, PRE##5##POST);	   \
  DEFINE_ANONYMOUS_ARGUMENT(6, TYPE, PRE##6##POST);	   \
  DEFINE_ANONYMOUS_ARGUMENT(7, TYPE, PRE##7##POST);	   \
  DEFINE_ANONYMOUS_ARGUMENT(8, TYPE, PRE##8##POST);	   \
  DEFINE_ANONYMOUS_ARGUMENT(9, TYPE, PRE##9##POST)
  
#define DEFINE_$_ANONYMOUS_ARGUMENTS(TYPE) \
  DEFINE_SET_OF_ANONYMOUS_ARGUMENTS($, , TYPE)
  
  } // namespace SpatialOps

#endif // SpatialOps_FieldOperations_h
