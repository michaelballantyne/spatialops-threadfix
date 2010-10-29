#ifndef SpatialOps_FieldOperations_h
#define SpatialOps_FieldOperations_h

namespace SpatialOps{
  
  /*
   * structure prototypes:
   */
  
  template<int ArgNum, typename UserType>
    struct ArgForm;

  template<typename AtomicType>
    struct Scalar;

  template<typename UserType>
    struct FieldForm;

  template <typename Operand1, typename Operand2, typename UserType>
    struct BinOp;

#define BUILD_BINARY_TYPE_PROTOTYPE(NAME)				\
  template <typename Operand1, typename Operand2, typename UserType>	\
    struct NAME;

  template <typename Operand, typename UserType>
    struct UnOp;

#define BUILD_UNARY_TYPE_PROTOTYPE(NAME)		\
  template <typename Operand, typename UserType>	\
    struct NAME;

  template<typename ExprType, typename UserType>
    struct Expression;

  template<int Num>
    struct ArgNum;

  template<typename Max1, typename Max2>
    struct CompareMaxArg;

  template<int Num1, int Num2, bool answer>
    struct InternalCompareMaxArg;

  template<typename Input, typename UserType>
    struct StandardizeTerm;

  template<typename NewType, typename OldType, typename UserType>
    struct LiftTerm;

  template<typename NewType, typename Type1, typename Type2, typename UserType>
    struct CombineTerms;

  template<typename Arg, typename UserType>
    struct StandardizeArg;

  template<typename BeginType, int CurrentArg, typename ArgType>
    struct ArgApply;

  template<typename BeginType, int CrtArgNum, typename Max, typename ArgType, typename UserType>
    struct AppResultFinder;

  template<typename ExprType, int CrtArgNum, typename Max, typename UserType>
    struct FcnForm;
  
  /* General definition. */
  template <typename StructTypeTemplate>
    struct Inline;
  

  /*
   * structure definitions:
   */

  /* Internal represenation of anonymous argument */
  template<int ArgNum, typename UserType>
    struct ArgForm {
      UserType typedef field_type;
  
      ArgForm()
      {};
    };

  /* Internal representation of constant doubles */
  template<typename AtomicType>
    struct Scalar {
      AtomicType value;
  
    Scalar(AtomicType const & v)
    : value(v)
      {};
    };

  /* Internal representation of vector/SpatialField */
  template<typename UserType>
    struct FieldForm {
      UserType const *fptr;
      typename UserType::const_iterator iter;
  
      //assumption: const_iterator is a pointer
      //TOCONSIDER: should iter be initialized?
    FieldForm(UserType const & field)
    : fptr(&field), iter(field.begin())
      {};
    };

  /* Internal representation of generic binary operation */
  template <typename Operand1, typename Operand2, typename UserType>
    struct BinOp {
      typename UserType::value_type typedef AtomicType;
  
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
  template <typename Operand1, typename Operand2, typename UserType>	\
    struct NAME {							\
    typename UserType::value_type typedef AtomicType;			\
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
  template <typename Operand, typename UserType>
    struct UnOp {
      typename UserType::value_type typedef AtomicType;
  
      AtomicType (*op)(AtomicType);
      Operand operand;
  
    UnOp(AtomicType (*operation)(AtomicType),
	 Operand oper)
    : op(operation), operand(oper)
      {};
    };

  /* Internal representation of specialized unary operation */
#define BUILD_UNARY_STRUCT(NAME)			\
  template <typename Operand, typename UserType>	\
    struct NAME {					\
    typename UserType::value_type typedef AtomicType;	\
    							\
    Operand operand;					\
    							\
    NAME(Operand op)					\
    : operand(op)					\
      {};						\
    };


  /* Wrapper/Container for passing AtomicType and UserType to expressions */
  template<typename Operand, typename UserType>
    struct Expression {
      UserType typedef field_type;
  
      Operand expr;
  
    Expression(Operand const & given)
    : expr(given)
      {};
    };

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

  /* Standardize UserType into Expression. */
  template<typename UserType>
    struct StandardizeTerm<UserType,UserType> {
    FieldForm<UserType> typedef StandardType;
    Expression<StandardType, UserType> typedef StandardTerm;
  
    static inline StandardType standardType (UserType const & given) {
      return StandardType(given);
    };
  
    static inline StandardTerm standardTerm (UserType const & given) {
      return StandardTerm(StandardType(given));
    };
  
  };

  /* Standardize Expression into itself. */
  template<typename ExprType, typename UserType>
    struct StandardizeTerm<Expression<ExprType,UserType>,UserType> {
    ExprType typedef StandardType;
    Expression<StandardType,UserType> typedef StandardTerm;
  
    static inline StandardType const & standardType (Expression<ExprType,UserType> const & given) {
      return given.expr;
    };
  
    static inline StandardTerm const & standardTerm (Expression<ExprType,UserType> const & given) {
      return given;
    };
  
  };

  /* Standardize ArgForm into FcnForm. */
  template<int Num, typename UserType>
    struct StandardizeTerm<ArgForm<Num,UserType>,UserType> {
    ArgForm<Num,UserType> typedef StandardType;
    FcnForm<StandardType,0,ArgNum<Num + 1>,UserType> typedef StandardTerm;
  
    static inline StandardType const & standardType (ArgForm<Num,UserType> const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (ArgForm<Num,UserType> const & given) {
      return StandardTerm(given);
    };
  
  };

  /* Standardize FcnForm into itself. */
  template<typename ExprType, typename Max, typename UserType>
    struct StandardizeTerm<FcnForm<ExprType,0,Max,UserType>,UserType> {
    ExprType typedef StandardType;
    FcnForm<StandardType,0,Max,UserType> typedef StandardTerm;
  
    static inline StandardType const & standardType (FcnForm<ExprType,0,Max,UserType> const & given) {
      return given.expr;
    };
  
    static inline StandardTerm const & standardTerm (FcnForm<ExprType,0,Max,UserType> const & given) {
      return given;
    };
  
  };

  /* Lift an Expresssion around a new type. */
  template<typename NewType, typename OldType, typename UserType>
    struct LiftTerm<NewType,
    Expression<OldType,UserType>,
    UserType> {
    NewType typedef StandardType;
    Expression<NewType,UserType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  };

  /* Lift a FcnForm around a new type. */
  template<typename NewType, typename OldType, typename Max, typename UserType>
    struct LiftTerm<NewType,
    FcnForm<OldType,0,Max,UserType>,
    UserType> {
    NewType typedef StandardType;
    FcnForm<NewType,0,Max,UserType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  };

  /* Combine two Expressions into a single Expression. */
  template<typename NewType, typename ExprType1, typename ExprType2, typename UserType>
    struct CombineTerms<NewType,
    Expression<ExprType1,UserType>,
    Expression<ExprType2,UserType>,
    UserType> {
    NewType typedef StandardType;
    Expression<NewType,UserType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  
  };

  /* Combine Expression and FcnForm into a FcnForm. */
  template<typename NewType, typename ExprType1, typename ExprType2, typename Max2, typename UserType>
    struct CombineTerms<NewType,
    Expression<ExprType1,UserType>,
    FcnForm<ExprType2,0,Max2,UserType>,
    UserType> {
    NewType typedef StandardType;
    FcnForm<NewType,0,Max2,UserType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  
  };

  /* Combine FcnForm and Expression into a FcnForm. */
  template<typename NewType, typename ExprType1, typename Max1, typename ExprType2, typename UserType>
    struct CombineTerms<NewType,
    FcnForm<ExprType1,0,Max1,UserType>,
    Expression<ExprType2,UserType>,
    UserType> {
    NewType typedef StandardType;
    FcnForm<NewType,0,Max1,UserType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  
  };

  /* Combine two FcnForms into a single FcnForm. */
  template<typename NewType, typename ExprType1, typename Max1, typename ExprType2, typename Max2, typename UserType>
    struct CombineTerms<NewType,
    FcnForm<ExprType1,0,Max1,UserType>,
    FcnForm<ExprType2,0,Max2,UserType>,
    UserType> {
    NewType typedef StandardType;
    FcnForm<NewType,
      0,
      typename CompareMaxArg<Max1,Max2>::Max,
      UserType> typedef StandardTerm;
  
    static inline StandardType const & standardType (NewType const & given) {
      return given;
    };
  
    static inline StandardTerm standardTerm (NewType const & given) {
      return StandardTerm(given);
    };
  
  };


  /* Standardize value_type into Scalar. */
  template<typename UserType>
    struct StandardizeArg<typename UserType::value_type,UserType> {
    typename UserType::value_type typedef AtomicType;
    Scalar<AtomicType> typedef StandardType;
    Expression<StandardType,UserType> typedef StandardTerm;
  
    static inline StandardType standardType (AtomicType const & given) {
      return StandardType(given);
    };
  
    static inline StandardTerm standardTerm (AtomicType const & given) {
      return StandardTerm(StandardType(given));
    };
  
  };

  /* Standardize UserType into FieldForm. */
  template<typename UserType>
    struct StandardizeArg<UserType,UserType> {
    FieldForm<UserType> typedef StandardType;
    Expression<StandardType,UserType> typedef StandardTerm;
  
    static inline StandardType standardType (UserType const & given) {
      return StandardType(given);
    };
  
    static inline StandardTerm standardTerm (UserType const & given) {
      return StandardTerm(StandardType(given));
    };
  
  };

  /* Standardize Expression into itself. */
  template<typename ExprType, typename UserType>
    struct StandardizeArg<Expression<ExprType,UserType>,UserType> {
    ExprType typedef StandardType;
    Expression<StandardType,UserType> typedef StandardTerm;
  
    static inline StandardType standardType (Expression<ExprType,UserType> const & given) {
      return given.expr;
    };
  
    static inline StandardTerm standardTerm (Expression<ExprType,UserType> const & given) {
      return given;
    };
  
  };
  
  
  /* 
   * ArgApply returns the type/state of applying an argument to an expression with anonymous arguments.
   */

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
  template<int CurrentArg, typename ArgType, typename UserType>
    struct ArgApply<FieldForm<UserType>,CurrentArg,ArgType> {
    FieldForm<UserType> typedef ReturnType;
  
    static inline ReturnType const & apply(FieldForm<UserType> const & state,
					   ArgType const & arg) {
      return state;
    };
  };

  /* Applying an parameter to the correct argument returns the parameter. */
  template<int ArgNum, typename ArgType, typename UserType>
    struct ArgApply<ArgForm<ArgNum,UserType>,ArgNum,ArgType> {
    ArgType typedef ReturnType;
  
    static inline ReturnType const & apply(ArgForm<ArgNum,UserType> const & state,
					   ArgType const & arg) {
      return arg;
    };
  };

  /* Applying an parameter to the wrong argument changes nothing. */
  template<int ArgNum, int CurrentArg, typename ArgType, typename UserType>
    struct ArgApply<ArgForm<ArgNum,UserType>,CurrentArg,ArgType> {
    ArgForm<ArgNum,UserType> typedef ReturnType;
  
    static inline ReturnType const & apply(ArgForm<ArgNum,UserType> const & state,
					   ArgType const & arg) {
      return state;
    };
  };

  /* Applying an parameter to a binary expression recurses on the subexpressions. */
  template<typename Operand1, typename Operand2, int CurrentArg, typename ArgType, typename UserType> 
    struct ArgApply<BinOp<Operand1,Operand2,UserType>,CurrentArg,ArgType> {
    BinOp<typename ArgApply<Operand1,CurrentArg,ArgType>::ReturnType,
      typename ArgApply<Operand2,CurrentArg,ArgType>::ReturnType,
      UserType> typedef ReturnType;
  
    static inline ReturnType apply(BinOp<Operand1,Operand2,UserType> const & state,
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
  template <typename Operand1, typename Operand2, int CurrentArg, typename ArgType, typename UserType> \
    struct ArgApply<NAME<Operand1,Operand2,UserType>,CurrentArg,ArgType> { \
    NAME<typename ArgApply<Operand1,CurrentArg,ArgType>::ReturnType,	\
      typename ArgApply<Operand2,CurrentArg,ArgType>::ReturnType,	\
      UserType> typedef ReturnType;					\
    									\
    static inline ReturnType apply(NAME<Operand1,Operand2,UserType> const & state, \
				   ArgType const & arg) {		\
      return ReturnType(ArgApply<Operand1,CurrentArg,ArgType>::apply(state.operand1, \
								     arg), \
			ArgApply<Operand2,CurrentArg,ArgType>::apply(state.operand2, \
								     arg)); \
    };									\
  };

  /* Applying an parameter to an unary expression recurses on the subexpression. */
  template<typename Operand, int CurrentArg, typename ArgType, typename UserType>
    struct ArgApply<UnOp<Operand,UserType>,CurrentArg,ArgType> {
    UnOp<typename ArgApply<Operand,CurrentArg,ArgType>::ReturnType,UserType> typedef ReturnType;
  
    static inline ReturnType apply(UnOp<Operand,UserType> const & state,
				   ArgType const & arg) {
      return ReturnType(state.op,
			ArgApply<Operand,CurrentArg,ArgType>::apply(state.operand,
								    arg));
    };
  };

  /* Applying an parameter to an unary expression recurses on the subexpression. */
#define BUILD_UNARY_ARGUMENT_APPLY(NAME)				\
  template <typename Operand, int CurrentArg, typename ArgType, typename UserType> \
    struct ArgApply<NAME<Operand,UserType>,CurrentArg,ArgType> {	\
    NAME<typename ArgApply<Operand,CurrentArg,ArgType>::ReturnType,UserType> typedef ReturnType; \
									\
    static inline ReturnType apply(NAME<Operand,UserType> const & state, \
				   ArgType const & arg) {		\
      return ReturnType(ArgApply<Operand,CurrentArg,ArgType>::apply(state.operand, \
								    arg)); \
    };									\
  };


  /*
   * AppResultFinder returns final result: Either wrapped in a FcnForm or not (remaining arguments or not).
   */

  /* Final argument is applied: No FcnForm wrapper. */
  template<typename BeginType, int CrtArgNum, typename ArgType, typename UserType>
    struct AppResultFinder<BeginType,CrtArgNum,ArgNum<CrtArgNum + 1>,ArgType,UserType> {
    Expression<typename ArgApply<BeginType,CrtArgNum,ArgType>::ReturnType,UserType> typedef ResultType;
  };

  /* Final argument has not been applied: Need FcnForm wrapper. */
  template<typename BeginType, int CrtArgNum, int MaxArgNum, typename ArgType, typename UserType>
    struct AppResultFinder<BeginType,CrtArgNum,ArgNum<MaxArgNum>,ArgType,UserType> {
    FcnForm<typename ArgApply<BeginType,CrtArgNum,ArgType>::ReturnType,
      CrtArgNum + 1,
      ArgNum<MaxArgNum>,
      UserType> typedef ResultType;
  };

  /*
   * FcnForm wraps expressions that have arguments still in them.
   */

  /* General definition. */
  template<typename ExprType, int CrtArgNum, int MaxArgNum, typename UserType>
    struct FcnForm<ExprType,CrtArgNum,ArgNum<MaxArgNum>,UserType> {
    UserType typedef field_type;
    ExprType typedef expr_type;
  
    ExprType expr;
  
  FcnForm(ExprType given)
    :expr(given)
    {};
  
    /* operator () definition (one argument version). */
    /* Most of the code here is to make sure the type rules match up. */
    template<typename ArgType>
      typename AppResultFinder<ExprType,
      CrtArgNum,
      ArgNum<MaxArgNum>,
      typename StandardizeArg<ArgType,
      UserType>::StandardType,
      UserType>::ResultType operator () (ArgType const & arg) {
    
      /* Wrapper type - if all arguments have been bound, wrapper is Expression; otherwise, wrapper is FcnForm. */
      typename AppResultFinder<ExprType,
	CrtArgNum,
	ArgNum<MaxArgNum>,
	typename StandardizeArg<ArgType,
	UserType>::StandardType,
	UserType>::ResultType typedef WrapperNextType;
  
      /* Actual type of expression after application of (standardized) ArgType. */
      ArgApply<ExprType,
	CrtArgNum,
	typename StandardizeArg<ArgType,
	UserType>::StandardType> typedef ActualNextType;
    
      /* Actual code that runs: Call correct apply function followed by a typecast. */
      return WrapperNextType(ActualNextType::apply(expr,
						   StandardizeArg<ArgType,
						   UserType>::standardType(arg)));
    };
  
    /*
     * Multiple arguments/currying the uncurried:
     */
  
    /* Two arguments. */
    template<typename Arg1, typename Arg2>
      typename AppResultFinder<typename ArgApply<ExprType,
      CrtArgNum,
      typename StandardizeArg<Arg1,
      UserType>::StandardType>::ReturnType,
      CrtArgNum + 1,
      ArgNum<MaxArgNum>,
      typename StandardizeArg<Arg2,
      UserType>::StandardType,
      UserType>::ResultType operator () (Arg1 const & arg1,
					 Arg2 const & arg2) {
      return this -> operator ()
	(arg1)
	(arg2);
    };
  
    /* Three arguments. */
    template<typename Arg1, typename Arg2, typename Arg3>
      typename AppResultFinder<typename ArgApply<typename ArgApply<ExprType,
      CrtArgNum,
      typename StandardizeArg<Arg1,
      UserType>::StandardType>::ReturnType,
      CrtArgNum + 1,
      typename StandardizeArg<Arg2,
      UserType>::StandardType>::ReturnType,
      CrtArgNum + 2,
      ArgNum<MaxArgNum>,
      typename StandardizeArg<Arg3,
      UserType>::StandardType,
      UserType>::ResultType operator () (Arg1 const & arg1,
					 Arg2 const & arg2,
					 Arg3 const & arg3) {
      return this -> operator ()
	(StandardizeArg<Arg1,
	 UserType>::standardType(arg1))
	(StandardizeArg<Arg2,
	 UserType>::standardType(arg2))
	(StandardizeArg<Arg3,
	 UserType>::standardType(arg3));
    };
  };


  /*
   * Inline inlines the correct functions and operators where appropriate.
   */

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
  template <typename UserType>
    struct Inline<FieldForm<UserType> > {
    typename UserType::value_type typedef AtomicType;
  
    static inline void init (FieldForm<UserType> & field) {
      field.iter = field.fptr->begin();
    };
  
    static inline void next (FieldForm<UserType> & field) {
      ++field.iter;
    };
  
    static inline AtomicType const & eval (FieldForm<UserType> const & field) {
      return *field.iter;
    };
  };

  /* Inline Expressions: */
  template <typename ExprType, typename UserType>
    struct Inline<Expression<ExprType,UserType> > {
    typename UserType::value_type typedef AtomicType;
  
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
  template <typename Operand1, typename Operand2, typename UserType> 
    struct Inline<BinOp<Operand1,Operand2,UserType> > {
    typename UserType::value_type typedef AtomicType;
  
    static inline void init (BinOp<Operand1,Operand2,UserType> & opStruct) {
      Inline<Operand1>::init(opStruct.operand1);
      Inline<Operand2>::init(opStruct.operand2);
    };
  
    static inline void next (BinOp<Operand1,Operand2,UserType> & opStruct) {
      Inline<Operand1>::next(opStruct.operand1);
      Inline<Operand2>::next(opStruct.operand2);
    };
  
    static inline AtomicType eval (BinOp<Operand1,Operand2,UserType> const & opStruct) {
      return opStruct.op(Inline<Operand1>::eval(opStruct.operand1),
			 Inline<Operand2>::eval(opStruct.operand2));
    };
  };

  /* Inline Binary Operators: */
#define BUILD_BINARY_OPERATOR_INLINER(TYPE, OPERATOR)			\
									\
  template <typename Operand1, typename Operand2, typename UserType>	\
    struct Inline<TYPE<Operand1,Operand2,UserType> > {			\
    typename UserType::value_type typedef AtomicType;			\
									\
    static inline void init (TYPE<Operand1,Operand2,UserType> & opStruct) { \
      Inline<Operand1>::init(opStruct.operand1);			\
      Inline<Operand2>::init(opStruct.operand2);			\
    };									\
									\
    static inline void next (TYPE<Operand1,Operand2,UserType> & opStruct) { \
      Inline<Operand1>::next(opStruct.operand1);			\
      Inline<Operand2>::next(opStruct.operand2);			\
    };									\
    									\
    static inline AtomicType eval (TYPE<Operand1,Operand2,UserType> const & opStruct) { \
      return Inline<Operand1>::eval(opStruct.operand1) OPERATOR		\
	Inline<Operand2>::eval(opStruct.operand2);			\
    };									\
  };

  /* Inline Binary Functions: */
#define BUILD_BINARY_FUNCTION_INLINER(TYPE, FUNCTION)			\
									\
  template <typename Operand1, typename Operand2, typename UserType>	\
    struct Inline<TYPE<Operand1,Operand2,UserType> > {			\
    typename UserType::value_type typedef AtomicType;			\
									\
    static inline void init (TYPE<Operand1,Operand2,UserType> & opStruct) { \
      Inline<Operand1>::init(opStruct.operand1);			\
      Inline<Operand2>::init(opStruct.operand2);			\
    };									\
									\
    static inline void next (TYPE<Operand1,Operand2,UserType> & opStruct) { \
      Inline<Operand1>::next(opStruct.operand1);			\
      Inline<Operand2>::next(opStruct.operand2);			\
    };									\
    									\
    static inline AtomicType eval (TYPE<Operand1,Operand2,UserType> const & opStruct) { \
      return FUNCTION(Inline<Operand1>::eval(opStruct.operand1),	\
		      Inline<Operand2>::eval(opStruct.operand2));	\
    };									\
  };

  /* Inline UnOps: */
  template <typename Operand, typename UserType>
    struct Inline<UnOp<Operand,UserType> > {
    typename UserType::value_type typedef AtomicType;
  
    static inline void init (UnOp<Operand,UserType> & opStruct) {
      Inline<Operand>::init(opStruct.operand);
    };
  
    static inline void next (UnOp<Operand,UserType> & opStruct) {
      Inline<Operand>::next(opStruct.operand);
    };
  
    static inline AtomicType eval (UnOp<Operand,UserType> const & opStruct) {
      return opStruct.op(Inline<Operand>::eval(opStruct.operand));
    };
  };

  /* Inline Unary Functions: */
#define BUILD_UNARY_FUNCTION_INLINER(TYPE, FUNCTION)			\
									\
  template <typename Operand, typename UserType>			\
    struct Inline<TYPE<Operand,UserType> > {				\
    typename UserType::value_type typedef AtomicType;			\
									\
    static inline void init (TYPE<Operand,UserType> & opStruct) {	\
      Inline<Operand>::init(opStruct.operand);				\
    }									\
    									\
    static inline void next (TYPE<Operand,UserType> & opStruct) {	\
      Inline<Operand>::next(opStruct.operand);				\
    }									\
    									\
    static inline AtomicType eval (TYPE<Operand,UserType> const & opStruct) { \
      return FUNCTION(Inline<Operand>::eval(opStruct.operand));		\
    };									\
  };


  /*
   * State/Expression Binary Operation interface:
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
  
    typename SubExpr1::field_type typedef UserType;
  
    typename StandardizeTerm<SubExpr1,UserType>::StandardTerm typedef Term1;
    typename StandardizeTerm<SubExpr2,UserType>::StandardTerm typedef Term2;
  
    typename StandardizeTerm<SubExpr1,UserType>::StandardType typedef Type1;
    typename StandardizeTerm<SubExpr2,UserType>::StandardType typedef Type2;
  
    BinOp<Type1,Type2,UserType> typedef ReturnType;
    typename CombineTerms<ReturnType,
      Term1,
      Term2,
      UserType>::StandardTerm typedef ReturnTerm;
  
    return ReturnTerm(ReturnType(fcn,
				 StandardizeTerm<SubExpr1,UserType>::standardType(first),
				 StandardizeTerm<SubExpr2,UserType>::standardType(second)));
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
  
    typename SubExpr::field_type typedef UserType;
    typename UserType::value_type typedef AtomicType;
  
    typename StandardizeTerm<SubExpr,UserType>::StandardTerm typedef Term1;
    Expression<Scalar<AtomicType>,UserType> typedef Term2;
  
    typename StandardizeTerm<SubExpr,UserType>::StandardType typedef Type1;
    Scalar<AtomicType> typedef Type2;
  
    BinOp<Type1,Type2,UserType> typedef ReturnType;
    typename CombineTerms<ReturnType,
      Term1,
      Term2,
      UserType>::StandardTerm typedef ReturnTerm;
  
    return ReturnTerm(ReturnType(fcn,
				 StandardizeTerm<SubExpr,UserType>::standardType(first),
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
  
    typename SubExpr::field_type typedef UserType;
    typename UserType::value_type typedef AtomicType;
  
    Expression<Scalar<AtomicType>,UserType> typedef Term1;
    typename StandardizeTerm<SubExpr,UserType>::StandardTerm typedef Term2;
  
    Scalar<AtomicType> typedef Type1;
    typename StandardizeTerm<SubExpr,UserType>::StandardType typedef Type2;
  
    BinOp<Type1,Type2,UserType> typedef ReturnType;
    typename CombineTerms<ReturnType,
      Term1,
      Term2,
      UserType>::StandardTerm typedef ReturnTerm;
  
    return ReturnTerm(ReturnType(fcn,
				 Type1(first),
				 StandardizeTerm<SubExpr,UserType>::standardType(second)));
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
      typename SubExpr1::field_type typedef UserType;			\
									\
      typename StandardizeTerm<SubExpr1,UserType>::StandardTerm typedef Term1; \
      typename StandardizeTerm<SubExpr2,UserType>::StandardTerm typedef Term2; \
									\
      typename StandardizeTerm<SubExpr1,UserType>::StandardType typedef Type1; \
      typename StandardizeTerm<SubExpr2,UserType>::StandardType typedef Type2; \
									\
      RETURN_TYPE<Type1,Type2,UserType> typedef ReturnType;		\
      typename CombineTerms<ReturnType,					\
	Term1,								\
	Term2,								\
	UserType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(StandardizeTerm<SubExpr1,UserType>::standardType(first), \
				   StandardizeTerm<SubExpr2,UserType>::standardType(second))); \
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
      typename SubExpr::field_type typedef UserType;			\
      typename UserType::value_type typedef AtomicType;			\
									\
      typename StandardizeTerm<SubExpr,UserType>::StandardTerm typedef Term1; \
      Expression<Scalar<AtomicType>,UserType> typedef Term2;		\
									\
      typename StandardizeTerm<SubExpr,UserType>::StandardType typedef Type1; \
      Scalar<AtomicType> typedef Type2;					\
    									\
      RETURN_TYPE<Type1,Type2,UserType> typedef ReturnType;		\
      typename CombineTerms<ReturnType,					\
	Term1,								\
	Term2,								\
	UserType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(StandardizeTerm<SubExpr,UserType>::standardType(first), \
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
      typename SubExpr::field_type typedef UserType;			\
      typename UserType::value_type typedef AtomicType;			\
									\
      Expression<Scalar<AtomicType>,UserType> typedef Term1;		\
      typename StandardizeTerm<SubExpr,UserType>::StandardTerm typedef Term2; \
									\
      Scalar<AtomicType> typedef Type1;					\
      typename StandardizeTerm<SubExpr,UserType>::StandardType typedef Type2; \
									\
      RETURN_TYPE<Type1,Type2,UserType> typedef ReturnType;		\
      typename CombineTerms<ReturnType,					\
	Term1,								\
	Term2,								\
	UserType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(Type1(first),			\
				   StandardizeTerm<SubExpr,UserType>::standardType(second))); \
    };
  
  
  /*
   * Unary Function Application
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
  
    typename SubExpr::field_type typedef UserType;
  
    typename StandardizeTerm<SubExpr,UserType>::StandardTerm typedef Term;
  
    typename StandardizeTerm<SubExpr,UserType>::StandardType typedef Type;
  
    UnOp<Type,UserType> typedef ReturnType;
    typename LiftTerm<ReturnType,
      Term,
      UserType>::StandardTerm typedef ReturnTerm;
  
    return ReturnTerm(ReturnType(fcn,
				 StandardizeTerm<SubExpr,UserType>::standardType(argument)));
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
      typename SubExpr::field_type typedef UserType;			\
									\
      typename StandardizeTerm<SubExpr,UserType>::StandardTerm typedef Term; \
									\
      typename StandardizeTerm<SubExpr,UserType>::StandardType typedef Type; \
									\
      RETURN_TYPE<Type,UserType> typedef ReturnType;			\
      typename LiftTerm<ReturnType,					\
	Term,								\
	UserType>::StandardTerm typedef ReturnTerm;			\
									\
      return ReturnTerm(ReturnType(StandardizeTerm<SubExpr,UserType>::standardType(argument))); \
    };


  /*
   * Assignment interface:
   */

  /* Assign a value_type to a UserType. */
  template <typename UserType>
    static inline void assign (UserType & lhs,
			       typename UserType::value_type const & given) {
    typename UserType::value_type typedef AtomicType;
    Scalar<AtomicType> typedef ExprType;
  
    ExprType state = ExprType(given);
  
    typename UserType::iterator iter = lhs.begin();
  
    Inline<ExprType>::init(state);
  
    for(Inline<ExprType>::init(state);
	iter != lhs.end();
	++iter,
	  Inline<ExprType>::next(state))
      *iter = Inline<ExprType>::eval(state);
  
  };

  /* Assign a UserType to a UserType. */
  template <typename UserType>
    static inline void assign (UserType & lhs,
			       UserType const & given) {
    FieldForm<UserType> typedef ExprType;
  
    ExprType state = ExprType(given);
  
    typename UserType::iterator iter = lhs.begin();
  
    Inline<ExprType>::init(state);
  
    for(Inline<ExprType>::init(state);
	iter != lhs.end();
	++iter,
	  Inline<ExprType>::next(state))
      *iter = Inline<ExprType>::eval(state);
  
  };

  /* Assign an Expression to a UserType. */
  template <typename ExprType, typename UserType>
    static inline void assign (UserType & lhs,
			       Expression<ExprType,
			       UserType> state) {
    typename UserType::iterator iter = lhs.begin();
  
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


    double add (double first, double second) {
    return first + second;
  };

  double subt (double first, double second) {
    return first - second;
  };

  double mult (double first, double second) {
    return first * second;
  };

  double div (double first, double second) {
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
