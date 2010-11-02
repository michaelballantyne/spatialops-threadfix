#ifndef SpatialOps_FieldOperations_h
#define SpatialOps_FieldOperations_h

namespace SpatialOps{
  
  /*
   * representation structure prototypes:
   */
  
  /**
   *  @struct 
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief 
   *  
   *  @par Template Parameters
   *   \li \b 
   *  
   *  
   */

  /**
   *  @struct Scalar
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Meta-computation representation of an element (AtomicType).
   *  
   *  @par Template Parameters
   *   \li \b AtomicType Basic element type.
   *  
   *  Scalar is a simple container structure that generalizes the use of AtomicType's in meta-computation.
   */
  template<typename AtomicType>
    struct Scalar;

  /**
   *  @struct FieldForm
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Meta-computation representation of a field (FieldType).
   *  
   *  @par Template Parameters
   *   \li \b FieldType Field type.
   *  
   *  @par
   *  
   *  FieldType is a simple container structure that generalizes the use of FieldType's in meta-computation.
   *
   *  FieldType must provide the following typedefs:
   *  
   *   \li \c field_type FieldType's own type. For example, if \c SpatialField is a FieldType, 
   *       then SpatialField must contain: \code typedef SpatialField field_type \endcode
   *  
   *   \li \c value_type FieldType's element/atomic type.  Similar to \c field_type, except
   *       \c value_type is the type of the elements contained within FieldType.  Usually
   *       referred to as AtomicType elsewhere in this documentation.
   *  
   *   \li \c iterator Type: \c value_type*.
   *  
   *   \li \c const_iterator Type: \c value_type \c const*.
   *  
   *  FieldType must provide the following methods:
   *  
   *   \li \code iterator begin() \endcode Returns a pointer to the first element in current FieldType.
   *  
   *   \li \code iterator end() \endcode Returns a pointer to the final element in current FieldType.
   *  
   *   \li \code iterator operator ++ (iterator &) \endcode and
   *       \code const_iterator operator ++ (const_iterator &) \endcode Increments given iterator/const_iterator. Return value is not used; only the side-effect is important.
   */
  template<typename FieldType>
    struct FieldForm;

  /**
   *  @struct BinOp/SpecificBinOp
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief General/specialized meta-computation representation of a binary operation/function.
   *  
   *  @par Template Parameters
   *   \li \b Operand1 First operand's type.
   *   \li \b Operand2 Second operand's type.
   *   \li \b FieldType Field type.
   *  
   *  @par BinOp
   *   BinOp is the non-optimized representation of a binary function and therefore requires the function be passed to it.
   *   Briefly, to use a function, \c fcn, that does not have a SpecificBinOp defined for it, with operands \c op1 and \c op2, the \c app function is called:
   *   \code app(fcn, op1, op2) \endcode.  The signature for \c app is:
   *   \code
   *   BinOp<Operand1, Operand2, FieldType> app (typename FieldType::value_type (*)(typename FieldType::value_type, typename FieldType::value_type),
   *                                             Operand1 const &,
   *                                             Operand2 const &)
   *   \endcode
   *  
   *  @par SpecificBinOp
   *   Commonly used binary operators and functions have been given individual representations, for both optimization and ease-of-use reasons.
   *   (These optimized BinOp-like structures are generated with macros by the preprocessor and so do not show up directly in the documentation.)
   *   To use these sturctures, the name given to the macro is used - usually identical to the operator or function itself.
   *   For example, for addition of two operands, \c op1 and \c op2, use the '+' symbol with infix notation:
   *   \code
   *   op1 + op2
   *   \endcode
   *   The macros, \c BUILD_BINARY_OPERATOR and \c BUILD_BINARY_FUNCTION, define binary operators and functions respectively.
   *   Note that usual C/C++ order of operations applies to binary operators defined in this way.
   */
  template <typename Operand1, typename Operand2, typename FieldType>
    struct BinOp;
  
#define BUILD_BINARY_TYPE_PROTOTYPE(NAME)				\
  template <typename Operand1, typename Operand2, typename FieldType>	\
    struct NAME;

  /**
   *  @struct UnOp/SpecificUnOp
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief General/specialized meta-computation representation of a unary function.
   *  
   *  @par Template Parameters
   *   \li \b Operand Operand's type.
   *   \li \b FieldType Field type.
   *  
   *  @par UnOp UnOp is the non-optimized representation of a unary function and therefore requires the function be passed to it.
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
   *  @par SpecificUnOp Commonly used unary functions have been given individual representations, for both optimization and ease-of-use reasons.
   *  (These optimized UnOp-like structures are generated with macros by the preprocessor and so do not show up directly in the documentation.)
   *  To use these sturctures, the name given to the macro is used - usually identical to the function itself.
   *  For example, for sin of the operand, \c op, usage is identical to applying sin to a double:
   *  \code
   *  sin(op)
   *  \endcode
   *  The macro, \c BUILD_UNARY_FUNCTION, defines unary functions.
   *  Note that usual C/C++ unary operators can be defined by macros very similar to \c BUILD_UNARY_FUNCTION; however, this macro has not been coded, because there is no immediate use.
   */
  template <typename Operand, typename FieldType>
    struct UnOp;

#define BUILD_UNARY_TYPE_PROTOTYPE(NAME)		\
  template <typename Operand, typename FieldType>	\
    struct NAME;

  /**
   *  @struct ArgForm
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Meta-computation representation of an anonymous argument.
   *  
   *  @par Template Parameters
   *   \li \b Num Argument number.
   *   \li \b FieldType Field type.
   *  
   *  @par
   *   
   *   Anonymous arguments are used to create anonymous functions (FcnForm).
   *   (Recursive definition: An anonymous function is an Expression containing at least one anonymous argument.)
   *   When expressions are applied to an anonymous function, each in order is associated with an integer, beginning at 0, for the argument it replaces.
   *   So in the following Expression,
   *   \code
   *   (4 + ($1 - $0))(a, b)
   *   \endcode
   *   the expression \c a is associated with 0, and the expression \c b with 1.
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
   *   will return a Scalar, with a value of 4, rather than causing a compilation error.
   */
  template<int Num, typename FieldType>
    struct ArgForm;
  
  /*
   * Container/wrapper structure prototypes:
   */
  
  /**
   *  @struct Expression
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Meta-computation representation of an expression to evaluate.
   *  
   *  @par Template Parameters
   *   \li \b ExprType The type representing the actual expression represented.
   *   \li \b FieldType Field type.
   *  
   *  @par
   *  
   *  An Expression represents an expression to evaluate.
   *  An Expression is either:
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
   *  The Expression structure itself is a container structure and does not contain itself, even though two Expressions may be used as operands to build a new Expression.
   *  Equivalently, the Expression structure only appears at the top of Expression and never appears in a subexpression.
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
   *  Expression<FieldForm<FT>, FT>
   *  \endcode
   *  Next, the compiler finds the type of the second subexpression \c 4.
   *  This examination finds the following return type:
   *  \code
   *  Expression<Scalar<AT>, FT>
   *  \endcode
   *  Finally, the two subexpressions are stripped of their Expression containers and used in a \c SumOp type; so the expression
   *  \code
   *  a + 4
   *  \endcode
   *  returns the type:
   *  \code
   *  Expression<SumOp<FieldForm<FT>,
   *                   Scalar<AT>,
   *                   FT>,
   *             FT>
   *  \endcode
   *
   *  \note There is currently no (direct) way to have a function that takes three or more arguments in an Expression.
   *  This can be added very easily (the implementation will be extremely similar to that of BinOp).
   *  There is no current use, so it was not implemented.
   */
  template<typename ExprType, typename FieldType>
    struct Expression;

  /* doxygen description below at partial specification */
  template<typename ExprType, int CrtNum, typename Max, typename FieldType>
    struct FcnForm;
  
  /*
   * Auxiliary/type-juggling structure prototypes:
   */
  
  template<int Num>
    struct ArgNum;

  template<typename Max1, typename Max2>
    struct CompareMaxArg;

  template<int Num1, int Num2, bool answer>
    struct InternalCompareMaxArg;

  template<typename Input, typename FieldType>
    struct StandardizeTerm;

  template<typename NewType, typename OldType, typename FieldType>
    struct LiftTerm;

  template<typename NewType, typename Type1, typename Type2, typename FieldType>
    struct CombineTerms;

  template<typename Arg, typename FieldType>
    struct StandardizeArg;

  template<typename BeginType, int CurrentArg, typename ArgType>
    struct ArgApply;

  template<typename BeginType, int CrtNum, typename Max, typename ArgType, typename FieldType>
    struct AppResultFinder;
  
  /*
   * Inliner structure prototype:
   */
  
  template <typename StructTypeTemplate>
    struct Inline;
  

  /*
   * representation structure definitions:
   */

  /* Internal representation of constant doubles */
  template<typename AtomicType>
    struct Scalar {
      AtomicType value;
  
    Scalar(AtomicType const & v)
    : value(v)
      {};
    };

  /* Internal representation of vector/SpatialField */
  template<typename FieldType>
    struct FieldForm {
      FieldType const *fptr;
      typename FieldType::const_iterator iter;
  
      //assumption: const_iterator is a pointer
      //TOCONSIDER: should iter be initialized?
    FieldForm(FieldType const & field)
    : fptr(&field), iter(field.begin())
      {};
    };

  /* Internal representation of generic binary operation */
  template <typename Operand1, typename Operand2, typename FieldType>
    struct BinOp {
      typename FieldType::value_type typedef AtomicType;
  
      AtomicType (*op)(AtomicType, AtomicType);
      Operand1 operand1;
      Operand2 operand2;
  
    BinOp(AtomicType (*operation)(AtomicType, AtomicType),
	  Operand1 op1,
	  Operand2 op2)
    : op(operation), operand1(op1), operand2(op2)
      {};
    };

  /* Internal representation of specialized binary operation */
#define BUILD_BINARY_STRUCT(NAME)					\
  template <typename Operand1, typename Operand2, typename FieldType>	\
    struct NAME {							\
    typename FieldType::value_type typedef AtomicType;			\
    									\
    Operand1 operand1;							\
    Operand2 operand2;							\
    									\
    NAME(Operand1 op1,							\
	 Operand2 op2)							\
    : operand1(op1), operand2(op2)					\
      {};								\
    };
  
  /* Internal representation of generic unary operation */
  template <typename Operand, typename FieldType>
    struct UnOp {
      typename FieldType::value_type typedef AtomicType;
  
      AtomicType (*op)(AtomicType);
      Operand operand;
  
    UnOp(AtomicType (*operation)(AtomicType),
	 Operand oper)
    : op(operation), operand(oper)
      {};
    };

  /* Internal representation of specialized unary operation */
#define BUILD_UNARY_STRUCT(NAME)			\
  template <typename Operand, typename FieldType>	\
    struct NAME {					\
    typename FieldType::value_type typedef AtomicType;	\
    							\
    Operand operand;					\
    							\
    NAME(Operand op)					\
    : operand(op)					\
      {};						\
    };

  /* Internal represenation of anonymous argument */
  template<int Num, typename FieldType>
    struct ArgForm {
      FieldType typedef field_type;
  
      ArgForm()
      {};
    };
  
  
  /*
   * Container/wrapper structure definitions:
   */
  
  /* Wrapper/Container for passing expressions around (no anonymous arguments). */
  template<typename Operand, typename FieldType>
    struct Expression {
      FieldType typedef field_type;
  
      Operand expr;
  
    Expression(Operand const & given)
    : expr(given)
      {};
    };

  /**
   *  @author Christopher Earl
   *  @date October, 2010
   *  
   *  @brief Meta-computation representation of an anonymous function.
   *  
   *  @par Template Parameters
   *   \li \b ExprType The type representing the function's body/actual expression represented.
   *   \li \b CrtNum Integer representing the next argument to be applied.
   *   \li \b MaxNum Integer representing the total number of arguments that need to be applied.
   *   \li \b FieldType Field type.
   *  
   *  An anonymous function is an Expression that contains at least one anonymous argument (see ArgForm).
   *  A FcnForm takes the place of Expression as the container structure used to standardize type interactions.
   *  A FcnForm also controls/regulates application of arguments.
   *  
   *  @par Application
   *   Application of arguments to an anonymous function simulates textual substitution.
   *   That is, when an Expression is applied to a FcnForm, wherever the next argument appears in the ExprType it is replaced by the applied Expression.
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
   *  @par Reversion back to Expression
   *   Once all the arguments that appear in an anonymous function are replaced by expressions via application, the container structure FcnForm is replaced by the general container structure Expression.
   *   This means that any fully-applied anonymous function can be treated like an Expression and used wherever an Expression is appropriate.
   *   The extended example below demonstrates this replacement.
   *  
   *  @par Currying
   *   Currying is the idea that
   *
   *  @par Extended example
   *   Consider the
   */
  template<typename ExprType, int CrtNum, int MaxNum, typename FieldType>
    struct FcnForm<ExprType,CrtNum,ArgNum<MaxNum>,FieldType> {
    FieldType typedef field_type;
    ExprType typedef expr_type;
  
    ExprType expr;
  
  FcnForm(ExprType given)
    :expr(given)
    {};
  
    /* operator () definition (one argument version). */
    /* Most of the code here is to make sure the type rules match up. */
    template<typename ArgType>
      typename AppResultFinder<ExprType,
      CrtNum,
      ArgNum<MaxNum>,
      typename StandardizeArg<ArgType,
      FieldType>::StandardType,
      FieldType>::ResultType operator () (ArgType const & arg) {
    
      /* Wrapper type - if all arguments have been bound, wrapper is Expression; otherwise, wrapper is FcnForm. */
      typename AppResultFinder<ExprType,
	CrtNum,
	ArgNum<MaxNum>,
	typename StandardizeArg<ArgType,
	FieldType>::StandardType,
	FieldType>::ResultType typedef WrapperNextType;
  
      /* Actual type of expression after application of (standardized) ArgType. */
      ArgApply<ExprType,
	CrtNum,
	typename StandardizeArg<ArgType,
	FieldType>::StandardType> typedef ActualNextType;
    
      /* Actual code that runs: Call correct apply function followed by a typecast. */
      return WrapperNextType(ActualNextType::apply(expr,
						   StandardizeArg<ArgType,
						   FieldType>::standardType(arg)));
    };
  
    /*
     * Multiple arguments/currying the uncurried:
     */
  
    /* Two arguments. */
    template<typename Arg1, typename Arg2>
      typename AppResultFinder<typename ArgApply<ExprType,
      CrtNum,
      typename StandardizeArg<Arg1,
      FieldType>::StandardType>::ReturnType,
      CrtNum + 1,
      ArgNum<MaxNum>,
      typename StandardizeArg<Arg2,
      FieldType>::StandardType,
      FieldType>::ResultType operator () (Arg1 const & arg1,
					 Arg2 const & arg2) {
      return this -> operator ()
	(arg1)
	(arg2);
    };
  
    /* Three arguments. */
    template<typename Arg1, typename Arg2, typename Arg3>
      typename AppResultFinder<typename ArgApply<typename ArgApply<ExprType,
      CrtNum,
      typename StandardizeArg<Arg1,
      FieldType>::StandardType>::ReturnType,
      CrtNum + 1,
      typename StandardizeArg<Arg2,
      FieldType>::StandardType>::ReturnType,
      CrtNum + 2,
      ArgNum<MaxNum>,
      typename StandardizeArg<Arg3,
      FieldType>::StandardType,
      FieldType>::ResultType operator () (Arg1 const & arg1,
					 Arg2 const & arg2,
					 Arg3 const & arg3) {
      return this -> operator ()
	(StandardizeArg<Arg1,
	 FieldType>::standardType(arg1))
	(StandardizeArg<Arg2,
	 FieldType>::standardType(arg2))
	(StandardizeArg<Arg3,
	 FieldType>::standardType(arg3));
    };
  };
  
  
  /*
   * Auxiliary/type-juggling structure prototypes:
   */
  
  /* An interger (number of argument) in type form. */
  template<int Num>
    struct ArgNum {};

  /* Calls InternalCompareMaxArg to find larger argument number. */
  /*  Also, forces templated types to be type-form of integers. */
  template<int Num1, int Num2>
    struct CompareMaxArg<ArgNum<Num1>,ArgNum<Num2> > {
    typename InternalCompareMaxArg<Num1,Num2,(Num1 > Num2)>::Max typedef Max;
  };
  
  /* If comparison is true, return first number. */
  template<int Num1, int Num2>
    struct InternalCompareMaxArg<Num1,Num2,true> {
    ArgNum<Num1> typedef Max;
  };

  /* If comparison is false, return second number. */
  template<int Num1, int Num2>
    struct InternalCompareMaxArg<Num1,Num2,false> {
    ArgNum<Num2> typedef Max;
  };

  /* Standardize FieldType into Expression. */
  template<typename FieldType>
    struct StandardizeTerm<FieldType,FieldType> {
    FieldForm<FieldType> typedef StandardType;
    Expression<StandardType, FieldType> typedef StandardTerm;
  
    static inline StandardType standardType (FieldType const & given) {
      return StandardType(given);
    };
  
    static inline StandardTerm standardTerm (FieldType const & given) {
      return StandardTerm(StandardType(given));
    };
  
  };

  /* Standardize Expression into itself. */
  template<typename ExprType, typename FieldType>
    struct StandardizeTerm<Expression<ExprType,FieldType>,FieldType> {
    ExprType typedef StandardType;
    Expression<StandardType,FieldType> typedef StandardTerm;
  
    static inline StandardType const & standardType (Expression<ExprType,FieldType> const & given) {
      return given.expr;
    };
  
    static inline StandardTerm const & standardTerm (Expression<ExprType,FieldType> const & given) {
      return given;
    };
  
  };

  /* Standardize ArgForm into FcnForm. */
  template<int Num, typename FieldType>
    struct StandardizeTerm<ArgForm<Num,FieldType>,FieldType> {
    ArgForm<Num,FieldType> typedef StandardType;
    FcnForm<StandardType,0,ArgNum<Num + 1>,FieldType> typedef StandardTerm;
  
    static inline StandardType const & standardType (ArgForm<Num,FieldType> const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (ArgForm<Num,FieldType> const & given) {
      return StandardTerm(given);
    };
  
  };

  /* Standardize FcnForm into itself. */
  template<typename ExprType, typename Max, typename FieldType>
    struct StandardizeTerm<FcnForm<ExprType,0,Max,FieldType>,FieldType> {
    ExprType typedef StandardType;
    FcnForm<StandardType,0,Max,FieldType> typedef StandardTerm;
  
    static inline StandardType const & standardType (FcnForm<ExprType,0,Max,FieldType> const & given) {
      return given.expr;
    };
  
    static inline StandardTerm const & standardTerm (FcnForm<ExprType,0,Max,FieldType> const & given) {
      return given;
    };
  
  };

  /* Lift an Expresssion around a new type. */
  template<typename NewType, typename OldType, typename FieldType>
    struct LiftTerm<NewType,
    Expression<OldType,FieldType>,
    FieldType> {
    NewType typedef StandardType;
    Expression<NewType,FieldType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  };

  /* Lift a FcnForm around a new type. */
  template<typename NewType, typename OldType, typename Max, typename FieldType>
    struct LiftTerm<NewType,
    FcnForm<OldType,0,Max,FieldType>,
    FieldType> {
    NewType typedef StandardType;
    FcnForm<NewType,0,Max,FieldType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  };

  /* Combine two Expressions into a single Expression. */
  template<typename NewType, typename ExprType1, typename ExprType2, typename FieldType>
    struct CombineTerms<NewType,
    Expression<ExprType1,FieldType>,
    Expression<ExprType2,FieldType>,
    FieldType> {
    NewType typedef StandardType;
    Expression<NewType,FieldType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  
  };

  /* Combine Expression and FcnForm into a FcnForm. */
  template<typename NewType, typename ExprType1, typename ExprType2, typename Max2, typename FieldType>
    struct CombineTerms<NewType,
    Expression<ExprType1,FieldType>,
    FcnForm<ExprType2,0,Max2,FieldType>,
    FieldType> {
    NewType typedef StandardType;
    FcnForm<NewType,0,Max2,FieldType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  
  };

  /* Combine FcnForm and Expression into a FcnForm. */
  template<typename NewType, typename ExprType1, typename Max1, typename ExprType2, typename FieldType>
    struct CombineTerms<NewType,
    FcnForm<ExprType1,0,Max1,FieldType>,
    Expression<ExprType2,FieldType>,
    FieldType> {
    NewType typedef StandardType;
    FcnForm<NewType,0,Max1,FieldType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  
  };

  /* Combine two FcnForms into a single FcnForm. */
  template<typename NewType, typename ExprType1, typename Max1, typename ExprType2, typename Max2, typename FieldType>
    struct CombineTerms<NewType,
    FcnForm<ExprType1,0,Max1,FieldType>,
    FcnForm<ExprType2,0,Max2,FieldType>,
    FieldType> {
    NewType typedef StandardType;
    FcnForm<NewType,
      0,
      typename CompareMaxArg<Max1,Max2>::Max,
      FieldType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  
  };


  /* Standardize value_type into Scalar. */
  template<typename FieldType>
    struct StandardizeArg<typename FieldType::value_type,FieldType> {
    typename FieldType::value_type typedef AtomicType;
    Scalar<AtomicType> typedef StandardType;
    Expression<StandardType,FieldType> typedef StandardTerm;
  
    static inline StandardType standardType (AtomicType const & given) {
      return StandardType(given);
    };
  
    static inline StandardTerm standardTerm (AtomicType const & given) {
      return StandardTerm(StandardType(given));
    };
  
  };

  /* Standardize FieldType into FieldForm. */
  template<typename FieldType>
    struct StandardizeArg<FieldType,FieldType> {
    FieldForm<FieldType> typedef StandardType;
    Expression<StandardType,FieldType> typedef StandardTerm;
  
    static inline StandardType standardType (FieldType const & given) {
      return StandardType(given);
    };
  
    static inline StandardTerm standardTerm (FieldType const & given) {
      return StandardTerm(StandardType(given));
    };
  
  };

  /* Standardize Expression into itself. */
  template<typename ExprType, typename FieldType>
    struct StandardizeArg<Expression<ExprType,FieldType>,FieldType> {
    ExprType typedef StandardType;
    Expression<StandardType,FieldType> typedef StandardTerm;
  
    static inline StandardType standardType (Expression<ExprType,FieldType> const & given) {
      return given.expr;
    };
  
    static inline StandardTerm standardTerm (Expression<ExprType,FieldType> const & given) {
      return given;
    };
  
  };
  
  
  /* ArgApply returns the type/state of applying an argument to an expression with anonymous arguments. */
  /* Applying an parameter to a Scalar changes nothing. */
  template<int CurrentArg, typename ArgType, typename AtomicType>
    struct ArgApply<Scalar<AtomicType>,CurrentArg,ArgType> {
    Scalar<AtomicType> typedef ReturnType;
  
    static inline ReturnType const & apply(Scalar<AtomicType> const & state,
					   ArgType const & arg) {
      return state;
    };
  };

  /* Applying an parameter to a FieldForm changes nothing. */
  template<int CurrentArg, typename ArgType, typename FieldType>
    struct ArgApply<FieldForm<FieldType>,CurrentArg,ArgType> {
    FieldForm<FieldType> typedef ReturnType;
  
    static inline ReturnType const & apply(FieldForm<FieldType> const & state,
					   ArgType const & arg) {
      return state;
    };
  };

  /* Applying an parameter to the correct argument returns the parameter. */
  template<int Num, typename ArgType, typename FieldType>
    struct ArgApply<ArgForm<Num,FieldType>,Num,ArgType> {
    ArgType typedef ReturnType;
  
    static inline ReturnType const & apply(ArgForm<Num,FieldType> const & state,
					   ArgType const & arg) {
      return arg;
    };
  };

  /* Applying an parameter to the wrong argument changes nothing. */
  template<int Num, int CurrentArg, typename ArgType, typename FieldType>
    struct ArgApply<ArgForm<Num,FieldType>,CurrentArg,ArgType> {
    ArgForm<Num,FieldType> typedef ReturnType;
  
    static inline ReturnType const & apply(ArgForm<Num,FieldType> const & state,
					   ArgType const & arg) {
      return state;
    };
  };

  /* Applying an parameter to a binary expression recurses on the subexpressions. */
  template<typename Operand1, typename Operand2, int CurrentArg, typename ArgType, typename FieldType> 
    struct ArgApply<BinOp<Operand1,Operand2,FieldType>,CurrentArg,ArgType> {
    BinOp<typename ArgApply<Operand1,CurrentArg,ArgType>::ReturnType,
      typename ArgApply<Operand2,CurrentArg,ArgType>::ReturnType,
      FieldType> typedef ReturnType;
  
    static inline ReturnType apply(BinOp<Operand1,Operand2,FieldType> const & state,
				   ArgType const & arg) {
      return ReturnType(state.op,
			ArgApply<Operand1,CurrentArg,ArgType>::apply(state.operand1,
								     arg),
			ArgApply<Operand2,CurrentArg,ArgType>::apply(state.operand2,
								     arg));
    };
  };

  /* Applying an parameter to a binary expression recurses on the subexpressions. */
#define BUILD_BINARY_ARGUMENT_APPLY(NAME)				\
  template <typename Operand1, typename Operand2, int CurrentArg, typename ArgType, typename FieldType> \
    struct ArgApply<NAME<Operand1,Operand2,FieldType>,CurrentArg,ArgType> { \
    NAME<typename ArgApply<Operand1,CurrentArg,ArgType>::ReturnType,	\
      typename ArgApply<Operand2,CurrentArg,ArgType>::ReturnType,	\
      FieldType> typedef ReturnType;					\
    									\
    static inline ReturnType apply(NAME<Operand1,Operand2,FieldType> const & state, \
				   ArgType const & arg) {		\
      return ReturnType(ArgApply<Operand1,CurrentArg,ArgType>::apply(state.operand1, \
								     arg), \
			ArgApply<Operand2,CurrentArg,ArgType>::apply(state.operand2, \
								     arg)); \
    };									\
  };

  /* Applying an parameter to an unary expression recurses on the subexpression. */
  template<typename Operand, int CurrentArg, typename ArgType, typename FieldType>
    struct ArgApply<UnOp<Operand,FieldType>,CurrentArg,ArgType> {
    UnOp<typename ArgApply<Operand,CurrentArg,ArgType>::ReturnType,FieldType> typedef ReturnType;
  
    static inline ReturnType apply(UnOp<Operand,FieldType> const & state,
				   ArgType const & arg) {
      return ReturnType(state.op,
			ArgApply<Operand,CurrentArg,ArgType>::apply(state.operand,
								    arg));
    };
  };

  /* Applying an parameter to an unary expression recurses on the subexpression. */
#define BUILD_UNARY_ARGUMENT_APPLY(NAME)				\
  template <typename Operand, int CurrentArg, typename ArgType, typename FieldType> \
    struct ArgApply<NAME<Operand,FieldType>,CurrentArg,ArgType> {	\
    NAME<typename ArgApply<Operand,CurrentArg,ArgType>::ReturnType,FieldType> typedef ReturnType; \
									\
    static inline ReturnType apply(NAME<Operand,FieldType> const & state, \
				   ArgType const & arg) {		\
      return ReturnType(ArgApply<Operand,CurrentArg,ArgType>::apply(state.operand, \
								    arg)); \
    };									\
  };


  /* AppResultFinder returns final result: Either wrapped in a FcnForm or not (remaining arguments or not). */
  /* Final argument is applied: No FcnForm wrapper. */
  template<typename BeginType, int CrtNum, typename ArgType, typename FieldType>
    struct AppResultFinder<BeginType,CrtNum,ArgNum<CrtNum + 1>,ArgType,FieldType> {
    Expression<typename ArgApply<BeginType,CrtNum,ArgType>::ReturnType,FieldType> typedef ResultType;
  };

  /* Final argument has not been applied: Need FcnForm wrapper. */
  template<typename BeginType, int CrtNum, int MaxNum, typename ArgType, typename FieldType>
    struct AppResultFinder<BeginType,CrtNum,ArgNum<MaxNum>,ArgType,FieldType> {
    FcnForm<typename ArgApply<BeginType,CrtNum,ArgType>::ReturnType,
      CrtNum + 1,
      ArgNum<MaxNum>,
      FieldType> typedef ResultType;
  };

  /*
   * Inliner structure definitions:
   */
  
  /* Inline inlines the correct functions and operators where appropriate. */
  /* Inline Scalars: */
  template<typename AtomicType>
    struct Inline<Scalar<AtomicType> > {
    static inline void init (Scalar<AtomicType> const & scalar) {
    };
  
    static inline void next (Scalar<AtomicType> const & scalar) {
    };
  
    static inline AtomicType const & eval (Scalar<AtomicType> const & scalar) {
      return scalar.value;
    };
  };

  /* Inline FieldForms: */
  template <typename FieldType>
    struct Inline<FieldForm<FieldType> > {
    typename FieldType::value_type typedef AtomicType;
  
    static inline void init (FieldForm<FieldType> & field) {
      field.iter = field.fptr->begin();
    };
  
    static inline void next (FieldForm<FieldType> & field) {
      ++field.iter;
    };
  
    static inline AtomicType const & eval (FieldForm<FieldType> const & field) {
      return *field.iter;
    };
  };

  /* Inline Expressions: */
  template <typename ExprType, typename FieldType>
    struct Inline<Expression<ExprType,FieldType> > {
    typename FieldType::value_type typedef AtomicType;
  
    static inline void init (ExprType & expr) {
      Inline<ExprType>::init(expr);
    };
  
    static inline void next (ExprType & expr) {
      Inline<ExprType>::next(expr);
    };
  
    static inline AtomicType const & eval (ExprType const & expr) {
      return Inline<ExprType>::eval(expr);
    };
  };

  /* Inline BinOps: */
  template <typename Operand1, typename Operand2, typename FieldType> 
    struct Inline<BinOp<Operand1,Operand2,FieldType> > {
    typename FieldType::value_type typedef AtomicType;
  
    static inline void init (BinOp<Operand1,Operand2,FieldType> & opStruct) {
      Inline<Operand1>::init(opStruct.operand1);
      Inline<Operand2>::init(opStruct.operand2);
    };
  
    static inline void next (BinOp<Operand1,Operand2,FieldType> & opStruct) {
      Inline<Operand1>::next(opStruct.operand1);
      Inline<Operand2>::next(opStruct.operand2);
    };
  
    static inline AtomicType eval (BinOp<Operand1,Operand2,FieldType> const & opStruct) {
      return opStruct.op(Inline<Operand1>::eval(opStruct.operand1),
			 Inline<Operand2>::eval(opStruct.operand2));
    };
  };

  /* Inline Binary Operators: */
#define BUILD_BINARY_OPERATOR_INLINER(TYPE, OPERATOR)			\
									\
  template <typename Operand1, typename Operand2, typename FieldType>	\
    struct Inline<TYPE<Operand1,Operand2,FieldType> > {			\
    typename FieldType::value_type typedef AtomicType;			\
									\
    static inline void init (TYPE<Operand1,Operand2,FieldType> & opStruct) { \
      Inline<Operand1>::init(opStruct.operand1);			\
      Inline<Operand2>::init(opStruct.operand2);			\
    };									\
									\
    static inline void next (TYPE<Operand1,Operand2,FieldType> & opStruct) { \
      Inline<Operand1>::next(opStruct.operand1);			\
      Inline<Operand2>::next(opStruct.operand2);			\
    };									\
    									\
    static inline AtomicType eval (TYPE<Operand1,Operand2,FieldType> const & opStruct) { \
      return Inline<Operand1>::eval(opStruct.operand1) OPERATOR		\
	Inline<Operand2>::eval(opStruct.operand2);			\
    };									\
  };

  /* Inline Binary Functions: */
#define BUILD_BINARY_FUNCTION_INLINER(TYPE, FUNCTION)			\
									\
  template <typename Operand1, typename Operand2, typename FieldType>	\
    struct Inline<TYPE<Operand1,Operand2,FieldType> > {			\
    typename FieldType::value_type typedef AtomicType;			\
									\
    static inline void init (TYPE<Operand1,Operand2,FieldType> & opStruct) { \
      Inline<Operand1>::init(opStruct.operand1);			\
      Inline<Operand2>::init(opStruct.operand2);			\
    };									\
									\
    static inline void next (TYPE<Operand1,Operand2,FieldType> & opStruct) { \
      Inline<Operand1>::next(opStruct.operand1);			\
      Inline<Operand2>::next(opStruct.operand2);			\
    };									\
    									\
    static inline AtomicType eval (TYPE<Operand1,Operand2,FieldType> const & opStruct) { \
      return FUNCTION(Inline<Operand1>::eval(opStruct.operand1),	\
		      Inline<Operand2>::eval(opStruct.operand2));	\
    };									\
  };

  /* Inline UnOps: */
  template <typename Operand, typename FieldType>
    struct Inline<UnOp<Operand,FieldType> > {
    typename FieldType::value_type typedef AtomicType;
  
    static inline void init (UnOp<Operand,FieldType> & opStruct) {
      Inline<Operand>::init(opStruct.operand);
    };
  
    static inline void next (UnOp<Operand,FieldType> & opStruct) {
      Inline<Operand>::next(opStruct.operand);
    };
  
    static inline AtomicType eval (UnOp<Operand,FieldType> const & opStruct) {
      return opStruct.op(Inline<Operand>::eval(opStruct.operand));
    };
  };

  /* Inline Unary Functions: */
#define BUILD_UNARY_FUNCTION_INLINER(TYPE, FUNCTION)			\
									\
  template <typename Operand, typename FieldType>			\
    struct Inline<TYPE<Operand,FieldType> > {				\
    typename FieldType::value_type typedef AtomicType;			\
									\
    static inline void init (TYPE<Operand,FieldType> & opStruct) {	\
      Inline<Operand>::init(opStruct.operand);				\
    }									\
    									\
    static inline void next (TYPE<Operand,FieldType> & opStruct) {	\
      Inline<Operand>::next(opStruct.operand);				\
    }									\
    									\
    static inline AtomicType eval (TYPE<Operand,FieldType> const & opStruct) { \
      return FUNCTION(Inline<Operand>::eval(opStruct.operand));		\
    };									\
  };


  /*
   * Binary Operator/Function interface:
   */

  /*
   * 
   * Cases for building binary operations:
   * 
   * Input X Input
   * Input X Scalar
   * Scalar X Input
   * Scalar X Scalar
   *
   *
   * Input \in {FieldType, Expression, ArgForm, FcnForm}
   *
   */

  /* Input X Input: */
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

  /* Input X Scalar: */
  template<typename SubExpr>
    typename CombineTerms<BinOp<typename StandardizeTerm<SubExpr,
    typename SubExpr::field_type>::StandardType,
    Scalar<typename SubExpr::field_type::value_type>,
    typename SubExpr::field_type>,
    typename StandardizeTerm<SubExpr,
    typename SubExpr::field_type>::StandardTerm,
    Expression<Scalar<typename SubExpr::field_type::value_type>,
    typename SubExpr::field_type>,
    typename SubExpr::field_type>::StandardTerm
    app (typename SubExpr::field_type::value_type (*fcn)(typename SubExpr::field_type::value_type,
							 typename SubExpr::field_type::value_type),
	 SubExpr const & first,
	 typename SubExpr::field_type::value_type const & second) {
  
    typename SubExpr::field_type typedef FieldType;
    typename FieldType::value_type typedef AtomicType;
  
    typename StandardizeTerm<SubExpr,FieldType>::StandardTerm typedef Term1;
    Expression<Scalar<AtomicType>,FieldType> typedef Term2;
  
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

  /* Scalar X Input: */
  template<typename SubExpr>
    typename CombineTerms<BinOp<Scalar<typename SubExpr::field_type::value_type>,
    typename StandardizeTerm<SubExpr,
    typename SubExpr::field_type>::StandardType,
    typename SubExpr::field_type>,
    Expression<Scalar<typename SubExpr::field_type::value_type>,
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
  
    Expression<Scalar<AtomicType>,FieldType> typedef Term1;
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


#define BUILD_BINARY_INTERFACE(RETURN_TYPE, FUNCTION_NAME)		\
  /* Input X Input: */							\
    template<typename SubExpr1, typename SubExpr2>			\
      typename CombineTerms<RETURN_TYPE<typename StandardizeTerm<SubExpr1, \
      typename SubExpr1::field_type>::StandardType,			\
      typename StandardizeTerm<SubExpr2,				\
      typename SubExpr1::field_type>::StandardType,			\
      typename SubExpr1::field_type>,					\
      typename StandardizeTerm<SubExpr1,				\
      typename SubExpr1::field_type>::StandardTerm,			\
      typename StandardizeTerm<SubExpr2,				\
      typename SubExpr1::field_type>::StandardTerm,			\
      typename SubExpr1::field_type>::StandardTerm			\
      FUNCTION_NAME (SubExpr1 const & first,				\
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
      RETURN_TYPE<Type1,Type2,FieldType> typedef ReturnType;		\
      typename CombineTerms<ReturnType,					\
	Term1,								\
	Term2,								\
	FieldType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(StandardizeTerm<SubExpr1,FieldType>::standardType(first), \
				   StandardizeTerm<SubExpr2,FieldType>::standardType(second))); \
    };									\
  									\
    /* Input X Scalar: */						\
    template<typename SubExpr>						\
      typename CombineTerms<RETURN_TYPE<typename StandardizeTerm<SubExpr, \
      typename SubExpr::field_type>::StandardType,			\
      Scalar<typename SubExpr::field_type::value_type>,			\
      typename SubExpr::field_type>,					\
      typename StandardizeTerm<SubExpr,					\
      typename SubExpr::field_type>::StandardTerm,			\
      Expression<Scalar<typename SubExpr::field_type::value_type>,	\
      typename SubExpr::field_type>,					\
      typename SubExpr::field_type>::StandardTerm			\
      FUNCTION_NAME (SubExpr const & first,				\
		     typename SubExpr::field_type::value_type const & second) { \
									\
      typename SubExpr::field_type typedef FieldType;			\
      typename FieldType::value_type typedef AtomicType;			\
									\
      typename StandardizeTerm<SubExpr,FieldType>::StandardTerm typedef Term1; \
      Expression<Scalar<AtomicType>,FieldType> typedef Term2;		\
									\
      typename StandardizeTerm<SubExpr,FieldType>::StandardType typedef Type1; \
      Scalar<AtomicType> typedef Type2;					\
    									\
      RETURN_TYPE<Type1,Type2,FieldType> typedef ReturnType;		\
      typename CombineTerms<ReturnType,					\
	Term1,								\
	Term2,								\
	FieldType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(StandardizeTerm<SubExpr,FieldType>::standardType(first), \
				   Type2(second)));			\
    };									\
									\
    /* Scalar X Input: */						\
    template<typename SubExpr>						\
      typename CombineTerms<RETURN_TYPE<Scalar<typename SubExpr::field_type::value_type>, \
      typename StandardizeTerm<SubExpr,					\
      typename SubExpr::field_type>::StandardType,			\
      typename SubExpr::field_type>,					\
      Expression<Scalar<typename SubExpr::field_type::value_type>,	\
      typename SubExpr::field_type>,					\
      typename StandardizeTerm<SubExpr,					\
      typename SubExpr::field_type>::StandardTerm,			\
      typename SubExpr::field_type>::StandardTerm			\
      FUNCTION_NAME (typename SubExpr::field_type::value_type const & first, \
		     SubExpr const & second) {				\
    									\
      typename SubExpr::field_type typedef FieldType;			\
      typename FieldType::value_type typedef AtomicType;			\
									\
      Expression<Scalar<AtomicType>,FieldType> typedef Term1;		\
      typename StandardizeTerm<SubExpr,FieldType>::StandardTerm typedef Term2; \
									\
      Scalar<AtomicType> typedef Type1;					\
      typename StandardizeTerm<SubExpr,FieldType>::StandardType typedef Type2; \
									\
      RETURN_TYPE<Type1,Type2,FieldType> typedef ReturnType;		\
      typename CombineTerms<ReturnType,					\
	Term1,								\
	Term2,								\
	FieldType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(Type1(first),			\
				   StandardizeTerm<SubExpr,FieldType>::standardType(second))); \
    };
  
  
  /*
   * UnaryFunction interface:
   */

  /*
   * 
   * Case for building unary operations:
   * 
   * Input
   *
   *
   * Input \in {FieldType, Expression, ArgForm, FcnForm}
   *
   */

  /* Input: */
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

#define BUILD_UNARY_INTERFACE(RETURN_TYPE, FUNCTION_NAME)		\
  /* Input: */								\
    template<typename SubExpr>						\
      typename LiftTerm<RETURN_TYPE<typename StandardizeTerm<SubExpr,	\
      typename SubExpr::field_type>::StandardType,			\
      typename SubExpr::field_type>,					\
      typename StandardizeTerm<SubExpr,					\
      typename SubExpr::field_type>::StandardTerm,			\
      typename SubExpr::field_type>::StandardTerm			\
      FUNCTION_NAME (SubExpr const & argument) {			\
    									\
      typename SubExpr::field_type typedef FieldType;			\
									\
      typename StandardizeTerm<SubExpr,FieldType>::StandardTerm typedef Term; \
									\
      typename StandardizeTerm<SubExpr,FieldType>::StandardType typedef Type; \
									\
      RETURN_TYPE<Type,FieldType> typedef ReturnType;			\
      typename LiftTerm<ReturnType,					\
	Term,								\
	FieldType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(StandardizeTerm<SubExpr,FieldType>::standardType(argument))); \
    };


  /*
   * Assignment defintions/interface:
   */

  /* Assign a value_type to a FieldType. */
  template <typename FieldType>
    static inline void assign (FieldType & lhs,
			       typename FieldType::value_type const & given) {
    typename FieldType::value_type typedef AtomicType;
    Scalar<AtomicType> typedef ExprType;
  
    ExprType state = ExprType(given);
  
    typename FieldType::iterator iter = lhs.begin();
  
    Inline<ExprType>::init(state);
  
    for(Inline<ExprType>::init(state);
	iter != lhs.end();
	++iter,
	  Inline<ExprType>::next(state))
      *iter = Inline<ExprType>::eval(state);
  
  };

  /* Assign a FieldType to a FieldType. */
  template <typename FieldType>
    static inline void assign (FieldType & lhs,
			       FieldType const & given) {
    FieldForm<FieldType> typedef ExprType;
  
    ExprType state = ExprType(given);
  
    typename FieldType::iterator iter = lhs.begin();
  
    Inline<ExprType>::init(state);
  
    for(Inline<ExprType>::init(state);
	iter != lhs.end();
	++iter,
	  Inline<ExprType>::next(state))
      *iter = Inline<ExprType>::eval(state);
  
  };

  /* Assign an Expression to a FieldType. */
  template <typename ExprType, typename FieldType>
    static inline void assign (FieldType & lhs,
			       Expression<ExprType,
			       FieldType> state) {
    typename FieldType::iterator iter = lhs.begin();
  
    Inline<ExprType>::init(state.expr);
  
    for(Inline<ExprType>::init(state.expr);
	iter != lhs.end();
	++iter,
	  Inline<ExprType>::next(state.expr))
      *iter = Inline<ExprType>::eval(state.expr);
  
  };


  /*
   * Macro's for binary operators and binary/unary functions:
   */

#define BUILD_BINARY_OPERATOR(NAME, OPERATOR)		\
  BUILD_BINARY_TYPE_PROTOTYPE(NAME)			\
    BUILD_BINARY_STRUCT(NAME)				\
    BUILD_BINARY_ARGUMENT_APPLY(NAME)			\
    BUILD_BINARY_OPERATOR_INLINER(NAME, OPERATOR)	\
    BUILD_BINARY_INTERFACE(NAME, operator OPERATOR)
  //Binary Operators:
  BUILD_BINARY_OPERATOR(SumOp, +)
    BUILD_BINARY_OPERATOR(DiffOp, -)
    BUILD_BINARY_OPERATOR(ProdOp, *)
    BUILD_BINARY_OPERATOR(DivOp, /)


  inline double add (double first, double second) {
    return first + second;
  };

  inline double subt (double first, double second) {
    return first - second;
  };

  inline double mult (double first, double second) {
    return first * second;
  };

  inline double div (double first, double second) {
    return first / second;
  };

#define BUILD_BINARY_FUNCTION(NAME, FUNCTION)		\
  BUILD_BINARY_TYPE_PROTOTYPE(NAME)			\
    BUILD_BINARY_STRUCT(NAME)				\
    BUILD_BINARY_ARGUMENT_APPLY(NAME)			\
    BUILD_BINARY_FUNCTION_INLINER(NAME, FUNCTION)	\
    BUILD_BINARY_INTERFACE(NAME, FUNCTION)
  //Binary Functions:
  BUILD_BINARY_FUNCTION(SumFcn,add)
    BUILD_BINARY_FUNCTION(SubtFcn,subt)
    BUILD_BINARY_FUNCTION(MultFcn,mult)
    BUILD_BINARY_FUNCTION(DivFcn,div)

#define BUILD_UNARY_FUNCTION(NAME, FUNCTION)		\
    BUILD_UNARY_TYPE_PROTOTYPE(NAME)			\
      BUILD_UNARY_STRUCT(NAME)				\
      BUILD_UNARY_ARGUMENT_APPLY(NAME)			\
      BUILD_UNARY_FUNCTION_INLINER(NAME, FUNCTION)	\
      BUILD_UNARY_INTERFACE(NAME, FUNCTION)
    //Unary Functions:
    BUILD_UNARY_FUNCTION(SinOp, sin)
    
} // namespace SpatialOps

#endif // SpatialOps_FieldOperations_h
