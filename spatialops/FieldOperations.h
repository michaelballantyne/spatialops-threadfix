#ifndef SpatialOps_FieldOperations_h
#define SpatialOps_FieldOperations_h

namespace SpatialOps{

  /*
   *
   * cwearl's additions:
   * SpatialField<VecOps,FieldLocation,GhostTraits>
   */
  
  typedef double D;
  
  /*
   * State/Expression prototypes:
   */

  struct Scalar;

  template<typename VecOps, typename FieldLocation, typename GhostTraits>
    struct FieldForm;

  template<int ArgNum>
    struct ArgForm;

  template <typename Operand1, typename Operand2> 
    struct BinOp;

#define BUILD_BINARY_TYPE_PROTOTYPE(NAME)		\
  template <typename Operand1, typename Operand2>	\
    struct NAME;
  
  template <typename Operand> 
    struct UnOp;
  
#define BUILD_UNARY_TYPE_PROTOTYPE(NAME)	\
  template <typename Operand>			\
    struct NAME;

  
  /*
   * Anonymous Argument prototypes:
   */

  template<typename Given>
    struct Convert;

  template<typename CurrentType, int CurrentArg, typename ArgType>
    struct ArgApply;

  template<int Num>
    struct ArgNum;

  template<typename Max1, typename Max2>
    struct CompareMaxArg;

  template<int Num1, int Num2, bool answer>
    struct InternalCompareMaxArg;

  template<typename BeginType, int CrtArgNum, typename max, typename ArgType>
    struct AppResultFinder;

  template<typename ActualType, int CrtArgNum, typename max>
    struct FcnForm;


  /*
   * State/Expression definitions:
   */

  /* Internal representation of constant doubles */
  struct Scalar {
    D value;
  
  Scalar(D v)
  : value(v)
    {};
  };
  
  /* Internal representation of vector/SpatialField */
  template<typename VecOps, typename FieldLocation, typename GhostTraits>
    struct FieldForm {
      SpatialField<VecOps,FieldLocation,GhostTraits>*fptr;
      typename SpatialField<VecOps,FieldLocation,GhostTraits>::const_iterator iter;
      
      //assumption: const_iterator is a pointer
      //TOCONSIDER: should iter be initialized?
    FieldForm(SpatialField<VecOps,FieldLocation,GhostTraits>* field)
    : fptr(field), iter(field->begin())
      {};
    };
  
  /* Internal represenation of anonymous argument */
  template<int ArgNum>
    struct ArgForm {
      ArgForm() {};
    };

  /* Internal representation of generic binary operation */
  template <typename Operand1, typename Operand2> 
    struct BinOp {
      D (*op)(D, D);
      Operand1 operand1;
      Operand2 operand2;
  
    BinOp(D (*operation)(D, D),
	  Operand1 op1,
	  Operand2 op2)
    : op(operation), operand1(op1), operand2(op2)
      {};
    };

  /* Internal representation of specialized binary operation */
#define BUILD_BINARY_STRUCT(NAME)			\
  template <typename Operand1, typename Operand2>	\
    struct NAME {					\
    Operand1 operand1;					\
    Operand2 operand2;					\
							\
    NAME(Operand1 op1,					\
	 Operand2 op2)					\
    : operand1(op1), operand2(op2)			\
      {};						\
    };

  /* Internal representation of generic unary operation */
  template <typename Operand>
    struct UnOp {
      D (*op)(D);
      Operand operand;
  
    UnOp(D (*operation)(D),
	 Operand oper)
    : op(operation), operand(oper)
      {}
    };

  /* Internal representation of specialized unary operation */
#define BUILD_UNARY_STRUCT(NAME)		\
  template <typename Operand>			\
    struct NAME {				\
    Operand operand;				\
						\
    NAME(Operand op)				\
    : operand(op)				\
      {}					\
    };


  /*
   * Anonymous Argument Definitions:
   */

  /*
   * Convert converts given type to internal representation.
   */

  /* Convert general type to itself. */
  template<typename Given>
    struct Convert {
      Given typedef ResultType;
  
      static Given convert(Given given) { return given; };
    };

  /* Convert double into Scalar. */
  template<>
    struct Convert<D> {
    Scalar typedef ResultType;
  
    static Scalar convert(D given) { return Scalar(given); };
  };

  /* Convert Field into FieldForm. */
  template<typename VecOps, typename FieldLocation, typename GhostTraits>
    struct Convert<SpatialField<VecOps,FieldLocation,GhostTraits> > {
    FieldForm<VecOps,FieldLocation,GhostTraits> typedef ResultType;
  
    static FieldForm<VecOps,FieldLocation,GhostTraits> convert(SpatialField<VecOps,FieldLocation,GhostTraits>& given) { 
      return FieldForm<VecOps,FieldLocation,GhostTraits>(& given);
    };
  };

  /* 
   * ArgApply returns the type/state of applying an argument to an expression with anonymous arguments.
   */

  /* General definition. */
  template<typename CurrentType, int CurrentArg, typename ArgType>
    struct ArgApply {};

  /* Applying an parameter to a Scalar changes nothing. */
  template<int CurrentArg, typename ArgType>
    struct ArgApply<Scalar,CurrentArg,ArgType> {
    typedef Scalar ReturnType;
  
    static inline ReturnType apply(Scalar state,
				   ArgType arg) {
      return state;
    };
  };

  /* Applying an parameter to a FieldForm changes nothing. */
  template<typename VecOps, typename FieldLocation, typename GhostTraits, int CurrentArg, typename ArgType>
    struct ArgApply<FieldForm<VecOps,FieldLocation,GhostTraits>,CurrentArg,ArgType> {
    typedef FieldForm<VecOps,FieldLocation,GhostTraits> ReturnType;
  
    static inline ReturnType apply(FieldForm<VecOps,FieldLocation,GhostTraits> state,
				   ArgType arg) {
      return state;
    };
  };

  /* Applying an parameter to the correct argument returns the parameter. */
  template<int ArgNum, typename ArgType>
    struct ArgApply<ArgForm<ArgNum>,ArgNum,ArgType> {
    typedef ArgType ReturnType;
  
    static inline ReturnType apply(ArgForm<ArgNum> state,
				   ArgType arg) {
      return arg;
    };
  };

  /* Applying an parameter to the wrong argument changes nothing. */
  template<int ArgNum, int CurrentArg, typename ArgType>
    struct ArgApply<ArgForm<ArgNum>,CurrentArg,ArgType> {
    typedef ArgForm<ArgNum> ReturnType;
  
    static inline ReturnType apply(ArgForm<ArgNum> state,
				   ArgType arg) {
      return state;
    };
  };

  /* Applying an parameter to a binary expression recurses on the subexpressions. */
  template<typename Operand1, typename Operand2, int CurrentArg, typename ArgType> 
    struct ArgApply<BinOp<Operand1,Operand2>,CurrentArg,ArgType> {
    BinOp<typename ArgApply<Operand1,CurrentArg,ArgType>::ReturnType,
      typename ArgApply<Operand2,CurrentArg,ArgType>::ReturnType> typedef ReturnType;
  
    static inline ReturnType apply(BinOp<Operand1,Operand2> state, ArgType arg) {
      return ReturnType(state.op,
			ArgApply<Operand1,CurrentArg,ArgType>::apply(state.operand1,
								     arg),
			ArgApply<Operand2,CurrentArg,ArgType>::apply(state.operand2,
								     arg));
    };
  };

  /* Applying an parameter to a binary expression recurses on the subexpressions. */
#define BUILD_BINARY_ARGUMENT_APPLY(NAME)				\
  template <typename Operand1, typename Operand2, int CurrentArg, typename ArgType> \
    struct ArgApply<NAME<Operand1,Operand2>,CurrentArg,ArgType> {	\
    NAME<typename ArgApply<Operand1,CurrentArg,ArgType>::ReturnType,	\
      typename ArgApply<Operand2,CurrentArg,ArgType>::ReturnType> typedef ReturnType; \
    									\
    static inline ReturnType apply(NAME<Operand1,Operand2> state, ArgType arg) { \
      return ReturnType(ArgApply<Operand1,CurrentArg,ArgType>::apply(state.operand1, \
								     arg), \
			ArgApply<Operand2,CurrentArg,ArgType>::apply(state.operand2, \
								     arg)); \
    };									\
  };

  /* Applying an parameter to an unary expression recurses on the subexpression. */
  template<typename Operand, int CurrentArg, typename ArgType> 
    struct ArgApply<UnOp<Operand>,CurrentArg,ArgType> {
    UnOp<typename ArgApply<Operand,CurrentArg,ArgType>::ReturnType> typedef ReturnType;
  
    static inline ReturnType apply(UnOp<Operand> state, ArgType arg) {
      return ReturnType(state.op,
			ArgApply<Operand,CurrentArg,ArgType>::apply(state.operand,
								    arg));
    };
  };

  /* Applying an parameter to an unary expression recurses on the subexpression. */
#define BUILD_UNARY_ARGUMENT_APPLY(NAME)				\
  template <typename Operand, int CurrentArg, typename ArgType>		\
    struct ArgApply<NAME<Operand>,CurrentArg,ArgType> {			\
    NAME<typename ArgApply<Operand,CurrentArg,ArgType>::ReturnType> typedef ReturnType; \
									\
    static inline ReturnType apply(NAME<Operand> state, ArgType arg) {	\
      return ReturnType(ArgApply<Operand,CurrentArg,ArgType>::apply(state.operand, \
								    arg)); \
    };									\
  };

  /* An interger (number of argument) in type form. */
  template<int Num>
    struct ArgNum {};

  /*
   * CompareMaxArg determines maximum argument number (in type form).
   */

  /* General definition. */
  template<typename Max1, typename Max2>
    struct CompareMaxArg {};

  /* Calls InternalCompareMaxArg to find larger argument number. */
  /*  Also, forces templated types to be type-form of integers. */
  template<int Num1, int Num2>
    struct CompareMaxArg<ArgNum<Num1>,ArgNum<Num2> > {
    typedef typename InternalCompareMaxArg<Num1,Num2,(Num1 > Num2)>::max max;
  };

  /* General definition. */
  template<int Num1, int Num2, bool answer>
    struct InternalCompareMaxArg {};

  /* If comparison is true, return first number. */
  template<int Num1, int Num2>
    struct InternalCompareMaxArg<Num1,Num2,true> {
    typedef ArgNum<Num1> max;
  };

  /* If comparison is false, return second number. */
  template<int Num1, int Num2>
    struct InternalCompareMaxArg<Num1,Num2,false> {
    typedef ArgNum<Num2> max;
  };

  /*
   * AppResultFinder returns final result: Either wrapped in a FcnForm or not (remaining arguments or not).
   */

  /* General defintion. */
  template<typename BeginType, int CrtArgNum, typename max, typename ArgType>
    struct AppResultFinder {};

  /* Final argument is applied: No FcnForm wrapper. */
  template<typename BeginType, int CrtArgNum, typename ArgType>
    struct AppResultFinder<BeginType,CrtArgNum,ArgNum<CrtArgNum + 1>,ArgType> {
    typename ArgApply<BeginType,CrtArgNum,ArgType>::ReturnType typedef ResultType;
  };

  /* Final argument has not been applied: Need FcnForm wrapper. */
  template<typename BeginType, int CrtArgNum, int MaxArgNum, typename ArgType>
    struct AppResultFinder<BeginType,CrtArgNum,ArgNum<MaxArgNum>,ArgType> {
    FcnForm<typename ArgApply<BeginType,CrtArgNum,ArgType>::ReturnType,
      CrtArgNum + 1,
      ArgNum<MaxArgNum> > typedef ResultType;
  };

  /*
   * FcnForm wraps expressions that have arguments still in them.
   */

  /* General defintion. */
  template<typename ActualType, int CrtArgNum, int MaxArgNum>
    struct FcnForm<ActualType,CrtArgNum,ArgNum<MaxArgNum> > {
    ActualType state;
  
  FcnForm(ActualType given)
    :state(given)
    {};
    
    /* operator () definition (one argument version). */
    typename AppResultFinder<ActualType,
      CrtArgNum,
      ArgNum<MaxArgNum>,
      Scalar>::ResultType operator () (D arg) {
    
      /* Wrapper type - if no wrapper (FcnForm) is needed, this is identical to ActualNextType. */
      typename AppResultFinder<ActualType,
	CrtArgNum,
	ArgNum<MaxArgNum>,
	Scalar>::ResultType typedef WrapperNextType;
    
      /* Actual type of expression after application of (converted) ArgType. */
      ArgApply<ActualType,
	CrtArgNum,
	Scalar> typedef ActualNextType;
    
      /* Actual code that runs: Call correct apply function followed by a typecast. */
      return WrapperNextType(ActualNextType::apply(state,
						   Scalar(arg)));
      };

    /* Most of the code here is to make sure the type rules match up. */
    template<typename ArgType>
      typename AppResultFinder<ActualType,
      CrtArgNum,
      ArgNum<MaxArgNum>,
      typename Convert<ArgType>::ResultType>::ResultType operator () (ArgType & arg) {
    
	/* Wrapper type - if no wrapper (FcnForm) is needed, this is identical to ActualNextType. */
	typename AppResultFinder<ActualType,
	  CrtArgNum,
	  ArgNum<MaxArgNum>,
	  typename Convert<ArgType>::ResultType>::ResultType typedef WrapperNextType;
    
	/* Actual type of expression after application of (converted) ArgType. */
	ArgApply<ActualType,
	  CrtArgNum,
	  typename Convert<ArgType>::ResultType> typedef ActualNextType;
    
	/* Actual code that runs: Call correct apply function followed by a typecast. */
	return WrapperNextType(ActualNextType::apply(state,
						     Convert<ArgType>::convert(arg)));
      };
  
    /*
     * Multiple arguments/currying the uncurried:
     */
  
    /* Two arguments. */
    template<typename Arg1, typename Arg2>
      typename AppResultFinder<typename ArgApply<ActualType,
      CrtArgNum,
      typename Convert<Arg1>::ResultType>::ReturnType,
      CrtArgNum + 1,
      ArgNum<MaxArgNum>,
      typename Convert<Arg2>::ResultType>::ResultType operator () (Arg1 arg1, Arg2 arg2) {
	return this -> operator ()
	  (Convert<Arg1>::convert(arg1))
	  (Convert<Arg2>::convert(arg2));
      };
  
    /* Three arguments. */
    template<typename Arg1, typename Arg2, typename Arg3>
      typename AppResultFinder<typename ArgApply<typename ArgApply<ActualType,
      CrtArgNum,
      typename Convert<Arg1>::ResultType>::ReturnType,
      CrtArgNum + 1,
      typename Convert<Arg2>::ResultType>::ReturnType,
      CrtArgNum + 2,
      ArgNum<MaxArgNum>,
      typename Convert<Arg2>::ResultType>::ResultType operator () (Arg1 arg1, Arg2 arg2, Arg3 arg3) {
	return this -> operator ()
	  (Convert<Arg1>::convert(arg1))
	  (Convert<Arg2>::convert(arg2))
	  (Convert<Arg3>::convert(arg3));
      };
  };


  /*
   * State/Expression inlining:
   */

  /*
   * Inline inlines the correct functions and operators where appropriate.
   */

  /* General defintion. */
  template <typename StructTypeTemplate>
    struct Inline {};

  /* Inline Scalars: */
  template<>
    struct Inline<Scalar> {
    static inline void init (Scalar& scalar) { }
  
    static inline void next (Scalar& scalar) { }
  
    static inline D eval (Scalar& scalar) { return scalar.value; }
  };

  /* Inline FieldForms: */
  template<typename VecOps, typename FieldLocation, typename GhostTraits>
    struct Inline<FieldForm<VecOps,FieldLocation,GhostTraits> > {
    static inline void init (FieldForm<VecOps,FieldLocation,GhostTraits>& field) {
      field.iter = field.fptr->begin();
    }
  
    static inline void next (FieldForm<VecOps,FieldLocation,GhostTraits>& field) {
      ++field.iter;
    }
  
    static inline D eval (FieldForm<VecOps,FieldLocation,GhostTraits>& field) {
      return *field.iter;
    }
  };


  /* Inline BinOps: */
  template <typename Operand1, typename Operand2> 
    struct Inline<BinOp<Operand1,Operand2> > {
    static inline void init (BinOp<Operand1,Operand2>& opStruct) {
      Inline<Operand1>::init(opStruct.operand1);
      Inline<Operand2>::init(opStruct.operand2);
    }
  
    static inline void next (BinOp<Operand1,Operand2>& opStruct) {
      Inline<Operand1>::next(opStruct.operand1);
      Inline<Operand2>::next(opStruct.operand2);
    }
  
    static inline D eval (BinOp<Operand1,Operand2>& opStruct) {
      return opStruct.op(Inline<Operand1>::eval(opStruct.operand1),
			 Inline<Operand2>::eval(opStruct.operand2));
    }
  };

  /* Inline Binary Operators: */
#define BUILD_BINARY_OPERATOR_INLINER(TYPE, OPERATOR)			\
									\
  template <typename Operand1, typename Operand2>			\
    struct Inline<TYPE<Operand1,Operand2> > {				\
    static inline void init (TYPE<Operand1,Operand2>& opStruct) {	\
      Inline<Operand1>::init(opStruct.operand1);			\
      Inline<Operand2>::init(opStruct.operand2);			\
    }									\
									\
    static inline void next (TYPE<Operand1,Operand2>& opStruct) {	\
      Inline<Operand1>::next(opStruct.operand1);			\
      Inline<Operand2>::next(opStruct.operand2);			\
    }									\
									\
    static inline D eval (TYPE<Operand1,Operand2>& opStruct) {		\
      return Inline<Operand1>::eval(opStruct.operand1) OPERATOR		\
	Inline<Operand2>::eval(opStruct.operand2);			\
    }									\
  };

  /* Inline Binary Functions: */
#define BUILD_BINARY_FUNCTION_INLINER(TYPE, FUNCTION)			\
									\
  template <typename Operand1, typename Operand2>			\
    struct Inline<TYPE<Operand1,Operand2> > {				\
    static inline void init (TYPE<Operand1,Operand2>& opStruct) {	\
      Inline<Operand1>::init(opStruct.operand1);			\
      Inline<Operand2>::init(opStruct.operand2);			\
    }									\
									\
    static inline void next (TYPE<Operand1,Operand2>& opStruct) {	\
      Inline<Operand1>::next(opStruct.operand1);			\
      Inline<Operand2>::next(opStruct.operand2);			\
    }									\
									\
    static inline D eval (TYPE<Operand1,Operand2>& opStruct) {		\
      return FUNCTION(Inline<Operand1>::eval(opStruct.operand1),	\
		      Inline<Operand2>::eval(opStruct.operand2));	\
    }									\
  };

  /* Inline UnOps: */
  template <typename Operand>
    struct Inline<UnOp<Operand> > {
    static inline void init (UnOp<Operand>& opStruct) {
      Inline<Operand>::init(opStruct.operand);
    }
  
    static inline void next (UnOp<Operand>& opStruct) {
      Inline<Operand>::next(opStruct.operand);
    }
  
    static inline D eval (UnOp<Operand>& opStruct) {
      return opStruct.op(Inline<Operand>::eval(opStruct.operand));
    }
  };

  /* Inline Unary Functions: */
#define BUILD_UNARY_FUNCTION_INLINER(TYPE, FUNCTION)		\
								\
  template <typename Operand>					\
    struct Inline<TYPE<Operand> > {				\
    static inline void init (TYPE<Operand>& opStruct) {		\
      Inline<Operand>::init(opStruct.operand);			\
    }								\
								\
    static inline void next (TYPE<Operand>& opStruct) {		\
      Inline<Operand>::next(opStruct.operand);			\
    }								\
								\
    static inline D eval (TYPE<Operand>& opStruct) {		\
      return FUNCTION(Inline<Operand>::eval(opStruct.operand)); \
    }								\
  };


  /*
   * State/Expression Binary Operation interface:
   */

  /*
   * 
   * Cases for building binary operations:
   * 
   * Scalar X Scalar - doesn't work for doubles (maybe)
   * Scalar X Field
   * Scalar X StructType
   * Scalar X AnonArg
   * Scalar X AnonFcn
   * Field X Scalar
   * Field X Field
   * Field X StructType
   * Field X AnonArg
   * Field X AnonFcn
   * StructType X Scalar
   * StructType X Field
   * StructType X StructType
   * StructType X AnonArg
   * StructType X AnonFcn
   * AnonArg X Scalar
   * AnonArg X Field
   * AnonArg X StructType
   * AnonArg X AnonArg
   * AnonArg X AnonFcn
   * AnonFcn X Scalar
   * AnonFcn X Field
   * AnonFcn X StructType
   * AnonFcn X AnonArg
   * AnonFcn X AnonFcn
   * 
   */

#define SCALAR_TYPE Scalar
#define SCALAR_PRIMATIVE_TYPE D
#define SCALAR_REFERENCE 

#define FIELD_TYPE FieldForm<VecOps,FieldLocation,GhostTraits>
#define FIELD_PRIMATIVE_TYPE SpatialField<VecOps,FieldLocation,GhostTraits>&
#define FIELD_REFERENCE &

  /*
   * Scalar X TYPE 
   */

  /* Scalar X Field */
  template<typename VecOps, typename FieldLocation, typename GhostTraits>
  BinOp<SCALAR_TYPE,FIELD_TYPE> app (D (*fcn)(D, D),
				     SCALAR_PRIMATIVE_TYPE first,
				     FIELD_PRIMATIVE_TYPE second) {
    return BinOp<SCALAR_TYPE,FIELD_TYPE> (fcn,
					  SCALAR_TYPE(SCALAR_REFERENCE first),
					  FIELD_TYPE(FIELD_REFERENCE second));
  };
 
  /* Scalar X StructType */
  template<typename Operand>
    BinOp<SCALAR_TYPE,Operand> app (D (*fcn)(D, D),
				    SCALAR_PRIMATIVE_TYPE first,
				    Operand second) {
    return BinOp<SCALAR_TYPE,Operand> (fcn,
				       SCALAR_TYPE(SCALAR_REFERENCE first),
				       second);
  };

  /* Scalar X AnonArg */
  template<int NumOfArg>
    FcnForm<BinOp<SCALAR_TYPE,ArgForm<NumOfArg> >,
    0,
    ArgNum<NumOfArg + 1> > app (D (* fcn)(D, D),
				SCALAR_PRIMATIVE_TYPE first,
				ArgForm<NumOfArg> second) {
    return FcnForm<BinOp<SCALAR_TYPE,ArgForm<NumOfArg> >,
      0,
      ArgNum<NumOfArg + 1> >(BinOp<SCALAR_TYPE,ArgForm<NumOfArg> >(fcn,
								   SCALAR_TYPE(SCALAR_REFERENCE first),
								   second));
  };

  /* Scalar X AnonFcn */
  template<typename ArgExp, typename max>
    FcnForm<BinOp<SCALAR_TYPE,ArgExp>,0,max> app (D (* fcn)(D, D),
						  SCALAR_PRIMATIVE_TYPE first,
						  FcnForm<ArgExp,0,max> second) {
    return FcnForm<BinOp<SCALAR_TYPE,ArgExp>,0,max>(BinOp<SCALAR_TYPE,ArgExp>(fcn,
									      SCALAR_TYPE(SCALAR_REFERENCE first),
									      second.state));
  };


  /*
   * Field X TYPE 
   */

  /* Field X Scalar */
  template<typename VecOps, typename FieldLocation, typename GhostTraits>
  BinOp<FIELD_TYPE,SCALAR_TYPE> app (D (*fcn)(D, D),
				     FIELD_PRIMATIVE_TYPE first,
				     SCALAR_PRIMATIVE_TYPE second) {
    return BinOp<FIELD_TYPE,SCALAR_TYPE> (fcn,
					  FIELD_TYPE(FIELD_REFERENCE first),
					  SCALAR_TYPE(SCALAR_REFERENCE second));
  };

  /* Field X Field */
  template<typename VecOps, typename FieldLocation, typename GhostTraits>
  BinOp<FIELD_TYPE,FIELD_TYPE> app (D (*fcn)(D, D),
				    FIELD_PRIMATIVE_TYPE first,
				    FIELD_PRIMATIVE_TYPE second) {
    return BinOp<FIELD_TYPE,FIELD_TYPE> (fcn,
					 FIELD_TYPE(FIELD_REFERENCE first),
					 FIELD_TYPE(FIELD_REFERENCE second));
  };

  /* Field X StructType */
  template<typename VecOps, typename FieldLocation, typename GhostTraits, typename Operand>
    BinOp<FIELD_TYPE,Operand> app (D (*fcn)(D, D),
				   FIELD_PRIMATIVE_TYPE first,
				   Operand second) {
    return BinOp<FIELD_TYPE,Operand> (fcn,
				      FIELD_TYPE(FIELD_REFERENCE first),
				      second);
  };

  /* Field X AnonArg */
  template<typename VecOps, typename FieldLocation, typename GhostTraits, int NumOfArg>
    FcnForm<BinOp<FIELD_TYPE,ArgForm<NumOfArg> >,
    0,
    ArgNum<NumOfArg + 1> > app (D (* fcn)(D, D),
				FIELD_PRIMATIVE_TYPE first,
				ArgForm<NumOfArg> second) {
    return FcnForm<BinOp<FIELD_TYPE,ArgForm<NumOfArg> >,
      0,
      ArgNum<NumOfArg + 1> >(BinOp<FIELD_TYPE,ArgForm<NumOfArg> >(fcn,
								  FIELD_TYPE(FIELD_REFERENCE first),
								  second));
  };

  /* Field X AnonFcn */
  template<typename VecOps, typename FieldLocation, typename GhostTraits, typename ArgExp, typename max>
    FcnForm<BinOp<FIELD_TYPE,ArgExp>,0,max> app (D (* fcn)(D, D),
						 FIELD_PRIMATIVE_TYPE first,
						 FcnForm<ArgExp,0,max> second) {
    return FcnForm<BinOp<FIELD_TYPE,ArgExp>,0,max>(BinOp<FIELD_TYPE,ArgExp>(fcn,
									    FIELD_TYPE(FIELD_REFERENCE first),
									    second.state));
  };


  /*
   * StructType X TYPE 
   */

  /* StructType X Scalar */
  template<typename Operand>
    BinOp<Operand,SCALAR_TYPE> app (D (*fcn)(D, D),
				    Operand first,
				    SCALAR_PRIMATIVE_TYPE second) {
    return BinOp<Operand,SCALAR_TYPE> (fcn,
				       first,
				       SCALAR_TYPE(SCALAR_REFERENCE second));
  };

  /* StructType X Field */
  template<typename VecOps, typename FieldLocation, typename GhostTraits, typename Operand>
    BinOp<Operand,FIELD_TYPE> app (D (*fcn)(D, D),
				   Operand first,
				   FIELD_PRIMATIVE_TYPE second) {
    return BinOp<Operand,FIELD_TYPE> (fcn,
				      first,
				      FIELD_TYPE(FIELD_REFERENCE second));
  };

  /* StructType X StructType */
  template<typename Operand1, typename Operand2>
    BinOp<Operand1,Operand2> app (D (*fcn)(D, D),
				  Operand1 first,
				  Operand2 second) {
    return BinOp<Operand1,Operand2> (fcn,
				     first,
				     second);
  };

  /* StructType X AnonArg */
  template<typename Operand, int NumOfArg>
    FcnForm<BinOp<Operand,ArgForm<NumOfArg> >,
    0,
    ArgNum<NumOfArg + 1> > app (D (*fcn)(D, D),
				Operand first,
				ArgForm<NumOfArg> second) {
    return FcnForm<BinOp<Operand,ArgForm<NumOfArg> >,
      0,
      ArgNum<NumOfArg + 1> >(BinOp<Operand,ArgForm<NumOfArg> >(fcn,
							       first,
							       second));
  };

  /* StructType X AnonFcn */
  template<typename Operand, typename ArgExp, typename max>
    FcnForm<BinOp<Operand,ArgExp>,0,max> app (D (*fcn)(D, D),
					      Operand first,
					      FcnForm<ArgExp,0,max> second) {
    return FcnForm<BinOp<Operand,ArgExp>,0,max>(BinOp<Operand,ArgExp>(fcn,
								      first,
								      second.state));
  };


  /*
   * AnonArg X TYPE 
   */

  /* AnonArg X Scalar */
  template<int NumOfArg>
    FcnForm<BinOp<ArgForm<NumOfArg>,SCALAR_TYPE>,
    0,
    ArgNum<NumOfArg + 1> > app (D (*fcn)(D, D),
				ArgForm<NumOfArg> first,
				SCALAR_PRIMATIVE_TYPE second) {
    return FcnForm<BinOp<ArgForm<NumOfArg>,SCALAR_TYPE>,
      0,
      ArgNum<NumOfArg + 1> > (BinOp<ArgForm<NumOfArg>,SCALAR_TYPE>(fcn,
								   first,
								   SCALAR_TYPE(SCALAR_REFERENCE second)));
  };

  /* AnonArg X Field */
  template<typename VecOps, typename FieldLocation, typename GhostTraits, int NumOfArg>
    FcnForm<BinOp<ArgForm<NumOfArg>,FIELD_TYPE>,
    0,
    ArgNum<NumOfArg + 1> > app (D (*fcn)(D, D),
				ArgForm<NumOfArg> first,
				FIELD_PRIMATIVE_TYPE second) {
    return FcnForm<BinOp<ArgForm<NumOfArg>,FIELD_TYPE>,
      0,
      ArgNum<NumOfArg + 1> >(BinOp<ArgForm<NumOfArg>,FIELD_TYPE> (fcn,
								  first,
								  FIELD_TYPE(FIELD_REFERENCE second)));
  };

  /* AnonArg X StructType */
  template<int NumOfArg, typename Operand>
    FcnForm<BinOp<ArgForm<NumOfArg>,Operand>,
    0,
    ArgNum<NumOfArg + 1> > app (D (*fcn)(D, D),
				ArgForm<NumOfArg> first,
				Operand second) {
    return FcnForm<BinOp<ArgForm<NumOfArg>,Operand>,
      0,
      ArgNum<NumOfArg + 1> >(BinOp<ArgForm<NumOfArg>,Operand> (fcn,
							       first,
							       second));
  };

  /* AnonArg X AnonArg */
  template<int NumOfArg1, int NumOfArg2>
    FcnForm<BinOp<ArgForm<NumOfArg1>,
    ArgForm<NumOfArg2> >,
    0,
    typename CompareMaxArg<ArgNum<NumOfArg1 + 1>,
    ArgNum<NumOfArg2 + 1> >::max > app (D (*fcn)(D, D),
					ArgForm<NumOfArg1> first,
					ArgForm<NumOfArg2> second) {
    return FcnForm<BinOp<ArgForm<NumOfArg1>,
      ArgForm<NumOfArg2> >,
      0,
      typename CompareMaxArg<ArgNum<NumOfArg1 + 1>,
      ArgNum<NumOfArg2 + 1> >::max >(BinOp<ArgForm<NumOfArg1>,
				     ArgForm<NumOfArg2> > (fcn,
							   first,
							   second));
  };

  /* AnonArg X AnonFcn */
  template<int NumOfArg1, typename ArgExp, int NumOfArg2>
    FcnForm<BinOp<ArgForm<NumOfArg1>,
    ArgExp>,
    0,
    typename CompareMaxArg<ArgNum<NumOfArg1 + 1>,
    ArgNum<NumOfArg2> >::max > app (D (*fcn)(D, D),
				    ArgForm<NumOfArg1> first,
				    FcnForm<ArgExp,0,ArgNum<NumOfArg2> > second) {
    return FcnForm<BinOp<ArgForm<NumOfArg1>,
      ArgExp>,
      0,
      typename CompareMaxArg<ArgNum<NumOfArg1 + 1>,
      ArgNum<NumOfArg2> >::max >(BinOp<ArgForm<NumOfArg1>,
				 ArgExp> (fcn,
					  first,
					  second.state));
  };

  /*
   * AnonFcn X TYPE
   */

  /* AnonFcn X Scalar */
  template<typename ArgExp, typename max>
    FcnForm<BinOp<ArgExp,SCALAR_TYPE>,0,max> app (D (*fcn)(D, D),
						  FcnForm<ArgExp,0,max> first,
						  SCALAR_PRIMATIVE_TYPE second) {
    return FcnForm<BinOp<ArgExp,SCALAR_TYPE>,0,max> (BinOp<ArgExp,SCALAR_TYPE>(fcn,
									       first.state,
									       SCALAR_TYPE(SCALAR_REFERENCE second)));
  };

  /* AnonFcn X Field */
  template<typename VecOps, typename FieldLocation, typename GhostTraits, typename ArgExp, typename max>
    FcnForm<BinOp<ArgExp,FIELD_TYPE>,0,max> app (D (*fcn)(D, D),
						 FcnForm<ArgExp,0,max> first,
						 FIELD_PRIMATIVE_TYPE second) {
    return FcnForm<BinOp<ArgExp,FIELD_TYPE>,0,max>(BinOp<ArgExp,FIELD_TYPE> (fcn,
									     first.state,
									     FIELD_TYPE(FIELD_REFERENCE second)));
  };

  /* AnonFcn X StructType */
  template<typename ArgExp, typename max, typename Operand>
    FcnForm<BinOp<ArgExp,Operand>,0,max> app (D (*fcn)(D, D),
					      FcnForm<ArgExp,0,max> first,
					      Operand second) {
    return FcnForm<BinOp<ArgExp,Operand>,0,max>(BinOp<ArgExp,Operand> (fcn,
								       first.state,
								       second));
  };

  /* AnonFcn X AnonArg */
  template<typename ArgExp, int NumOfArg1, int NumOfArg2>
    FcnForm<BinOp<ArgExp,
    ArgForm<NumOfArg2> >,
    0,
    typename CompareMaxArg<ArgNum<NumOfArg1>,
    ArgNum<NumOfArg2 + 1> >::max > app (D (*fcn)(D, D),
					FcnForm<ArgExp,0,ArgNum<NumOfArg1> > first,
					ArgForm<NumOfArg2> second) {
    return FcnForm<BinOp<ArgExp,
      ArgForm<NumOfArg2> >,
      0,
      typename CompareMaxArg<ArgNum<NumOfArg1>,
      ArgNum<NumOfArg2 + 1> >::max >(BinOp<ArgExp,
				     ArgForm<NumOfArg2> > (fcn,
							   first.state,
							   second));
  };

  /* AnonFcn X AnonFcn */
  template<typename ArgExp1, typename max1, typename ArgExp2, typename max2>
    FcnForm<BinOp<ArgExp1,ArgExp2>,0, typename CompareMaxArg<max1,max2>::max > app (D (*fcn)(D, D),
										    FcnForm<ArgExp1,0,max1> first,
										    FcnForm<ArgExp2,0,max2> second) {
    return FcnForm<BinOp<ArgExp1,ArgExp2>,0,typename CompareMaxArg<max1,max2>::max >(BinOp<ArgExp1,ArgExp2> (fcn,
													     first.state,
													     second.state));
  };


#define BUILD_BINARY_INTERFACE(RETURN_TYPE, FUNCTION_NAME)		\
  /*									\
   * Scalar X TYPE							\
   */									\
									\
    /* Scalar X Field */						\
    template<typename VecOps, typename FieldLocation, typename GhostTraits> \
      RETURN_TYPE<SCALAR_TYPE,FIELD_TYPE> FUNCTION_NAME (SCALAR_PRIMATIVE_TYPE first, \
							 FIELD_PRIMATIVE_TYPE second) { \
      return RETURN_TYPE<SCALAR_TYPE,FIELD_TYPE> (SCALAR_TYPE(SCALAR_REFERENCE first), \
						  FIELD_TYPE(FIELD_REFERENCE second)); \
    };									\
    									\
    /* Scalar X StructType */						\
    template<typename Operand>						\
      RETURN_TYPE<SCALAR_TYPE,Operand> FUNCTION_NAME (SCALAR_PRIMATIVE_TYPE first, \
						      Operand second) {	\
      return RETURN_TYPE<SCALAR_TYPE,Operand> (SCALAR_TYPE(SCALAR_REFERENCE first), \
					       second);			\
    };									\
  									\
    /* Scalar X AnonArg */						\
    template<int NumOfArg>						\
      FcnForm<RETURN_TYPE<SCALAR_TYPE,ArgForm<NumOfArg> >,		\
      0,								\
      ArgNum<NumOfArg + 1> > FUNCTION_NAME (SCALAR_PRIMATIVE_TYPE first, \
					    ArgForm<NumOfArg> second) { \
      return FcnForm<RETURN_TYPE<SCALAR_TYPE,ArgForm<NumOfArg> >,	\
	0,								\
	ArgNum<NumOfArg + 1> >(RETURN_TYPE<SCALAR_TYPE,ArgForm<NumOfArg> >(SCALAR_TYPE(SCALAR_REFERENCE first), \
									   second)); \
    };									\
  									\
    /* Scalar X AnonFcn */						\
    template<typename ArgExp, typename max>				\
      FcnForm<RETURN_TYPE<SCALAR_TYPE,ArgExp>,				\
      0,								\
      max> FUNCTION_NAME (SCALAR_PRIMATIVE_TYPE first,			\
			  FcnForm<ArgExp,0,max> second) {		\
      return FcnForm<RETURN_TYPE<SCALAR_TYPE,ArgExp>,			\
	0,								\
	max>(RETURN_TYPE<SCALAR_TYPE,ArgExp>(SCALAR_TYPE(SCALAR_REFERENCE first), \
					     second.state));		\
    };									\
									\
    /*									\
     * Field X TYPE							\
     */									\
									\
    /* Field X Scalar */						\
    template<typename VecOps, typename FieldLocation, typename GhostTraits> \
      RETURN_TYPE<FIELD_TYPE,SCALAR_TYPE> FUNCTION_NAME (FIELD_PRIMATIVE_TYPE first, \
							 SCALAR_PRIMATIVE_TYPE second) { \
      return RETURN_TYPE<FIELD_TYPE,SCALAR_TYPE> (FIELD_TYPE(FIELD_REFERENCE first), \
						  SCALAR_TYPE(SCALAR_REFERENCE second)); \
    };									\
									\
    /* Field X Field */							\
    template<typename VecOps, typename FieldLocation, typename GhostTraits> \
    RETURN_TYPE<FIELD_TYPE,FIELD_TYPE> FUNCTION_NAME (FIELD_PRIMATIVE_TYPE first, \
						      FIELD_PRIMATIVE_TYPE second) { \
      return RETURN_TYPE<FIELD_TYPE,FIELD_TYPE> (FIELD_TYPE(FIELD_REFERENCE first), \
						 FIELD_TYPE(FIELD_REFERENCE second)); \
    };									\
									\
    /* Field X StructType */						\
    template<typename VecOps, typename FieldLocation, typename GhostTraits, typename Operand> \
      RETURN_TYPE<FIELD_TYPE,Operand> FUNCTION_NAME (FIELD_PRIMATIVE_TYPE first, \
						     Operand second) {	\
      return RETURN_TYPE<FIELD_TYPE,Operand> (FIELD_TYPE(FIELD_REFERENCE first), \
					      second);			\
    };									\
									\
									\
    /* Field X AnonArg */						\
    template<typename VecOps, typename FieldLocation, typename GhostTraits, int NumOfArg> \
      FcnForm<RETURN_TYPE<FIELD_TYPE,ArgForm<NumOfArg> >,		\
      0,								\
      ArgNum<NumOfArg + 1> > FUNCTION_NAME (FIELD_PRIMATIVE_TYPE first, \
					    ArgForm<NumOfArg> second) { \
      return FcnForm<RETURN_TYPE<FIELD_TYPE,ArgForm<NumOfArg> >,	\
	0,								\
	ArgNum<NumOfArg + 1> >(RETURN_TYPE<FIELD_TYPE,ArgForm<NumOfArg> >(FIELD_TYPE(FIELD_REFERENCE first), \
									  second)); \
    };									\
									\
    /* Field X AnonFcn */						\
    template<typename VecOps, typename FieldLocation, typename GhostTraits, typename ArgExp, typename max> \
      FcnForm<RETURN_TYPE<FIELD_TYPE,ArgExp>,				\
      0,								\
      max> FUNCTION_NAME (FIELD_PRIMATIVE_TYPE first,			\
			  FcnForm<ArgExp,0,max> second) {		\
      return FcnForm<RETURN_TYPE<FIELD_TYPE,ArgExp>,			\
	0,								\
	max>(RETURN_TYPE<FIELD_TYPE,ArgExp>(FIELD_TYPE(FIELD_REFERENCE first), \
					    second.state));		\
    };									\
									\
    /*									\
     * StructType X TYPE						\
     */									\
									\
    /* StructType X Scalar */						\
    template<typename Operand>						\
      RETURN_TYPE<Operand,SCALAR_TYPE> FUNCTION_NAME (Operand first,	\
						      SCALAR_PRIMATIVE_TYPE second) { \
      return RETURN_TYPE<Operand,SCALAR_TYPE> (first,			\
					       SCALAR_TYPE(SCALAR_REFERENCE second)); \
    };									\
									\
    /* StructType X Field */						\
    template<typename VecOps, typename FieldLocation, typename GhostTraits, typename Operand> \
      RETURN_TYPE<Operand,FIELD_TYPE> FUNCTION_NAME (Operand first,	\
						     FIELD_PRIMATIVE_TYPE second) { \
      return RETURN_TYPE<Operand,FIELD_TYPE> (first,			\
					      FIELD_TYPE(FIELD_REFERENCE second)); \
    };									\
									\
    /* StructType X StructType */					\
    template<typename Operand1, typename Operand2>			\
      RETURN_TYPE<Operand1,Operand2> FUNCTION_NAME (Operand1 first,	\
						    Operand2 second) {	\
      return RETURN_TYPE<Operand1,Operand2> (first,			\
					     second);			\
    };									\
									\
    /* StructType X AnonArg */						\
    template<typename Operand, int NumOfArg>				\
      FcnForm<RETURN_TYPE<Operand,ArgForm<NumOfArg> >,			\
      0,								\
      ArgNum<NumOfArg + 1> > FUNCTION_NAME (Operand first,		\
					    ArgForm<NumOfArg> second) { \
      return FcnForm<RETURN_TYPE<Operand,ArgForm<NumOfArg> >,		\
	0,								\
	ArgNum<NumOfArg + 1> >(RETURN_TYPE<Operand,ArgForm<NumOfArg> >(first, \
								       second)); \
    };									\
									\
    /* StructType X AnonFcn */						\
    template<typename Operand, typename ArgExp, typename max>		\
      FcnForm<RETURN_TYPE<Operand,ArgExp>,				\
      0,								\
      max> FUNCTION_NAME (Operand first,				\
			  FcnForm<ArgExp,0,max> second) {		\
      return FcnForm<RETURN_TYPE<Operand,ArgExp>,			\
	0,								\
	max>(RETURN_TYPE<Operand,ArgExp>(first,				\
					 second.state));		\
    };									\
  									\
    /*									\
     * AnonArg X TYPE							\
     */									\
  									\
    /* AnonArg X Scalar */						\
    template<int NumOfArg>						\
      FcnForm<RETURN_TYPE<ArgForm<NumOfArg>,SCALAR_TYPE>,		\
      0,								\
      ArgNum<NumOfArg + 1> > FUNCTION_NAME (ArgForm<NumOfArg> first,	\
					    SCALAR_PRIMATIVE_TYPE second) { \
      return FcnForm<RETURN_TYPE<ArgForm<NumOfArg>,SCALAR_TYPE>,	\
	0,								\
	ArgNum<NumOfArg + 1> > (RETURN_TYPE<ArgForm<NumOfArg>,SCALAR_TYPE>(first, \
									   SCALAR_TYPE(SCALAR_REFERENCE second))); \
    };									\
									\
    /* AnonArg X Field */						\
    template<typename VecOps, typename FieldLocation, typename GhostTraits, int NumOfArg> \
      FcnForm<RETURN_TYPE<ArgForm<NumOfArg>,FIELD_TYPE>,		\
      0,								\
      ArgNum<NumOfArg + 1> > FUNCTION_NAME (ArgForm<NumOfArg> first,	\
					    FIELD_PRIMATIVE_TYPE second) { \
      return FcnForm<RETURN_TYPE<ArgForm<NumOfArg>,FIELD_TYPE>,		\
	0,								\
	ArgNum<NumOfArg + 1> >(RETURN_TYPE<ArgForm<NumOfArg>,FIELD_TYPE> (first, \
									  FIELD_TYPE(FIELD_REFERENCE second))); \
    };									\
  									\
    /* AnonArg X StructType */						\
    template<int NumOfArg, typename Operand>				\
      FcnForm<RETURN_TYPE<ArgForm<NumOfArg>,Operand>,			\
      0,								\
      ArgNum<NumOfArg + 1> > FUNCTION_NAME (ArgForm<NumOfArg> first,	\
					    Operand second) {		\
      return FcnForm<RETURN_TYPE<ArgForm<NumOfArg>,Operand>,		\
	0,								\
	ArgNum<NumOfArg + 1> >(RETURN_TYPE<ArgForm<NumOfArg>,Operand> (first, \
								       second)); \
    };									\
  									\
    /* AnonArg X AnonArg */						\
    template<int NumOfArg1, int NumOfArg2>				\
      FcnForm<RETURN_TYPE<ArgForm<NumOfArg1>,				\
      ArgForm<NumOfArg2> >,						\
      0,								\
      typename CompareMaxArg<ArgNum<NumOfArg1 + 1>,			\
      ArgNum<NumOfArg2 + 1> >::max > FUNCTION_NAME (ArgForm<NumOfArg1> first, \
						    ArgForm<NumOfArg2> second) { \
      return FcnForm<RETURN_TYPE<ArgForm<NumOfArg1>,			\
	ArgForm<NumOfArg2> >,						\
	0,								\
	typename CompareMaxArg<ArgNum<NumOfArg1 + 1>,			\
	ArgNum<NumOfArg2 + 1> >::max >(RETURN_TYPE<ArgForm<NumOfArg1>,ArgForm<NumOfArg2> > (first, \
											    second)); \
    };									\
  									\
    /* AnonArg X AnonFcn */						\
    template<int NumOfArg1, typename ArgExp, int NumOfArg2>		\
      FcnForm<RETURN_TYPE<ArgForm<NumOfArg1>,				\
      ArgExp>,								\
      0,								\
      typename CompareMaxArg<ArgNum<NumOfArg1 + 1>,			\
      ArgNum<NumOfArg2> >::max > FUNCTION_NAME (ArgForm<NumOfArg1> first, \
						FcnForm<ArgExp,0,ArgNum<NumOfArg2> > second) { \
      return FcnForm<RETURN_TYPE<ArgForm<NumOfArg1>,			\
	ArgExp>,							\
	0,								\
	typename CompareMaxArg<ArgNum<NumOfArg1 + 1>,			\
	ArgNum<NumOfArg2> >::max >(RETURN_TYPE<ArgForm<NumOfArg1>,	\
				   ArgExp> (first,			\
					    second.state));		\
    };									\
  									\
    /*									\
     * AnonFcn X TYPE							\
     */									\
  									\
    /* AnonFcn X Scalar */						\
    template<typename ArgExp, typename max>				\
      FcnForm<RETURN_TYPE<ArgExp,SCALAR_TYPE>,				\
      0,								\
      max> FUNCTION_NAME (FcnForm<ArgExp,0,max> first,			\
			  SCALAR_PRIMATIVE_TYPE second) {		\
      return FcnForm<RETURN_TYPE<ArgExp,SCALAR_TYPE>,			\
	0,								\
	max> (RETURN_TYPE<ArgExp,SCALAR_TYPE>(first.state,		\
					      SCALAR_TYPE(SCALAR_REFERENCE second))); \
    };									\
									\
    /* AnonFcn X Field */						\
    template<typename VecOps, typename FieldLocation, typename GhostTraits, typename ArgExp, typename max> \
      FcnForm<RETURN_TYPE<ArgExp,FIELD_TYPE>,				\
      0,								\
      max> FUNCTION_NAME (FcnForm<ArgExp,0,max> first,			\
			  FIELD_PRIMATIVE_TYPE second) {		\
      return FcnForm<RETURN_TYPE<ArgExp,FIELD_TYPE>,			\
	0,								\
	max>(RETURN_TYPE<ArgExp,FIELD_TYPE> (first.state,		\
					     FIELD_TYPE(FIELD_REFERENCE second))); \
    };									\
  									\
    /* AnonFcn X StructType */						\
    template<typename ArgExp, typename max, typename Operand>		\
      FcnForm<RETURN_TYPE<ArgExp,Operand>,				\
      0,								\
      max> FUNCTION_NAME (FcnForm<ArgExp,0,max> first,			\
			  Operand second) {				\
      return FcnForm<RETURN_TYPE<ArgExp,Operand>,			\
	0,								\
	max>(RETURN_TYPE<ArgExp,Operand> (first.state,			\
					  second));			\
    };									\
  									\
    /* AnonFcn X AnonArg */						\
    template<typename ArgExp, int NumOfArg1, int NumOfArg2>		\
      FcnForm<RETURN_TYPE<ArgExp,					\
      ArgForm<NumOfArg2> >,						\
      0,								\
      typename CompareMaxArg<ArgNum<NumOfArg1>,				\
      ArgNum<NumOfArg2 + 1> >::max > FUNCTION_NAME (FcnForm<ArgExp,0,ArgNum<NumOfArg1> > first, \
						    ArgForm<NumOfArg2> second) { \
      return FcnForm<RETURN_TYPE<ArgExp,				\
	ArgForm<NumOfArg2> >,						\
	0,								\
	typename CompareMaxArg<ArgNum<NumOfArg1>,			\
	ArgNum<NumOfArg2 + 1> >::max >(RETURN_TYPE<ArgExp,		\
				       ArgForm<NumOfArg2> > (first.state, \
							     second));	\
    };									\
									\
    /* AnonFcn X AnonFcn */						\
    template<typename ArgExp1, typename max1, typename ArgExp2, typename max2> \
      FcnForm<RETURN_TYPE<ArgExp1,ArgExp2>,				\
      0,								\
      typename CompareMaxArg<max1,max2>::max > FUNCTION_NAME (FcnForm<ArgExp1,0,max1> first, \
							      FcnForm<ArgExp2,0,max2> second) { \
      return FcnForm<RETURN_TYPE<ArgExp1,ArgExp2>,			\
	0,								\
	typename CompareMaxArg<max1,max2>::max >(RETURN_TYPE<ArgExp1,ArgExp2> (first.state, \
									       second.state)); \
    };

  /*
   * Unary Function Application
   */

  /* Scalar */
/*   UnOp<SCALAR_TYPE> app (SCALAR_PRIMATIVE_TYPE (*fcn)(SCALAR_PRIMATIVE_TYPE), */
/* 			 SCALAR_PRIMATIVE_TYPE operand) { */
/*     return UnOp<SCALAR_TYPE> (fcn, */
/* 			      SCALAR_TYPE(SCALAR_REFERENCE operand)); */
/*   }; */
  template<typename D>
  UnOp<SCALAR_TYPE> app (D (*fcn)(D),
			 D operand) {
    return UnOp<SCALAR_TYPE> (fcn,
			      SCALAR_TYPE(SCALAR_REFERENCE operand));
  };

  /* Field */
  template<typename VecOps, typename FieldLocation, typename GhostTraits>
  UnOp<FIELD_TYPE> app (SCALAR_PRIMATIVE_TYPE (*fcn)(SCALAR_PRIMATIVE_TYPE),
			FIELD_PRIMATIVE_TYPE operand) {
    return UnOp<FIELD_TYPE> (fcn,
			     FIELD_TYPE(FIELD_REFERENCE operand));
  };

  /* StructType */
  template<typename Operand>
    UnOp<Operand> app (SCALAR_PRIMATIVE_TYPE (*fcn)(SCALAR_PRIMATIVE_TYPE),
		       Operand operand) {
    return UnOp<Operand> (fcn,
			  operand);
  };

  /* AnonArg */
  template<int NumOfArg>
    FcnForm<UnOp<ArgForm<NumOfArg> >,
    0,
    ArgNum<NumOfArg + 1> > app (SCALAR_PRIMATIVE_TYPE (*fcn)(SCALAR_PRIMATIVE_TYPE),
				ArgForm<NumOfArg> operand) {
    return FcnForm<UnOp<ArgForm<NumOfArg> >,
      0,
      ArgNum<NumOfArg + 1> >(UnOp<ArgForm<NumOfArg> >(fcn,
						      operand));
  };

  /* AnonFcn */
  template<typename ArgExp, typename max>
    FcnForm<UnOp<ArgExp>,
    0,
    max> app (SCALAR_PRIMATIVE_TYPE (*fcn)(SCALAR_PRIMATIVE_TYPE),
	      FcnForm<ArgExp,0,max> operand) {
    return FcnForm<UnOp<ArgExp>,
      0,
      max>(UnOp<ArgExp>(fcn,
			operand.state));
  };

#define BUILD_UNARY_INTERFACE(RETURN_TYPE, FUNCTION_NAME)		\
									\
  /* Field */								\
    template<typename VecOps, typename FieldLocation, typename GhostTraits> \
      RETURN_TYPE<FIELD_TYPE> FUNCTION_NAME (FIELD_PRIMATIVE_TYPE operand) { \
      return RETURN_TYPE<FIELD_TYPE> (FIELD_TYPE(FIELD_REFERENCE operand)); \
    };									\
									\
    /* StructType */							\
    template<typename Operand>						\
      RETURN_TYPE<Operand> FUNCTION_NAME (Operand operand) {		\
      return RETURN_TYPE<Operand> (operand);				\
    };									\
									\
    /* AnonArg */							\
    template<int NumOfArg>						\
      FcnForm<RETURN_TYPE<ArgForm<NumOfArg> >,				\
      0,								\
      ArgNum<NumOfArg + 1> > FUNCTION_NAME (ArgForm<NumOfArg> operand) { \
      return FcnForm<RETURN_TYPE<ArgForm<NumOfArg> >,			\
	0,								\
	ArgNum<NumOfArg + 1> >(RETURN_TYPE<ArgForm<NumOfArg> >(operand)); \
    };									\
									\
    /* AnonFcn */							\
    template<typename ArgExp, typename max>				\
      FcnForm<RETURN_TYPE<ArgExp>,					\
      0,								\
      max> FUNCTION_NAME (FcnForm<ArgExp,0,max> operand) {		\
      return FcnForm<RETURN_TYPE<ArgExp>,				\
	0,								\
	max>(RETURN_TYPE<ArgExp>(operand.state));			\
    };


  /*
   * Assignment interface:
   */

  template<typename VecOps, typename FieldLocation, typename GhostTraits, typename StructTypeTemplate>
    static inline void assign (SpatialField<VecOps,FieldLocation,GhostTraits>& lhs, StructTypeTemplate state) {
    typename SpatialField<VecOps,FieldLocation,GhostTraits>::iterator iter = lhs.begin();
  
/*     typename Convert<StructTypeTemplate>::ResultType typedef ActualType; */
/*     ActualType realState = Convert<StructTypeTemplate>::convert(state); */
  
/*     Inline<ActualType>::init(realState); */
  
/*     for(Inline<ActualType>::init(realState); */
/* 	iter != lhs.end(); */
/* 	++iter, */
/* 	  Inline<ActualType>::next(realState)) */
/*       *iter = Inline<ActualType>::eval(realState); */
  
      Inline<StructTypeTemplate>::init(state);
  
      for(Inline<StructTypeTemplate>::init(state);
          iter != lhs.end();
          ++iter,
    	Inline<StructTypeTemplate>::next(state))
        *iter = Inline<StructTypeTemplate>::eval(state);
  
  }

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
  //BUILD_BINARY_OPERATOR(SumOp, +)
  //BUILD_BINARY_OPERATOR(DiffOp, -)
  //BUILD_BINARY_OPERATOR(ProdOp, *)
  //BUILD_BINARY_OPERATOR(DivOp, /)

  template<typename Return>
    D add (D first, D second) {
    return first + second;
  };

  template<typename Return>
    Return subt (D first, D second) {
    return first - second;
  };
  
  template<typename Return>
    Return mult (D first, D second) {
    return first * second;
  };
  
  template<typename Return>
    Return div (D first, D second) {
    return first / second;
  };
  
#define BUILD_BINARY_FUNCTION(NAME, FUNCTION, FUNCTION_NAME)	\
  BUILD_BINARY_TYPE_PROTOTYPE(NAME)				\
    BUILD_BINARY_STRUCT(NAME)					\
    BUILD_BINARY_ARGUMENT_APPLY(NAME)				\
    BUILD_BINARY_FUNCTION_INLINER(NAME, FUNCTION)		\
    BUILD_BINARY_INTERFACE(NAME, FUNCTION_NAME)
  //Binary Functions:
  BUILD_BINARY_FUNCTION(SumFcn,add<D>,add)
    BUILD_BINARY_FUNCTION(SubtFcn,subt<D>,subt)
    BUILD_BINARY_FUNCTION(MultFcn,mult<D>,mult)
    BUILD_BINARY_FUNCTION(DivFcn,div<D>,div)

#define BUILD_UNARY_FUNCTION(NAME, FUNCTION)		\
    BUILD_UNARY_TYPE_PROTOTYPE(NAME)			\
      BUILD_UNARY_STRUCT(NAME)				\
      BUILD_UNARY_ARGUMENT_APPLY(NAME)			\
      BUILD_UNARY_FUNCTION_INLINER(NAME, FUNCTION)	\
      BUILD_UNARY_INTERFACE(NAME, FUNCTION)
    //Unary Functions:
    BUILD_UNARY_FUNCTION(SinOp, sin)
    
    
  /*
   * end of cwearl's additions
   */

} // namespace SpatialOps

#endif // SpatialOps_FieldOperations_h

