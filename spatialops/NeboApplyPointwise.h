#ifndef NEBO_APPLY_POINTWISE_H
   #define NEBO_APPLY_POINTWISE_H
   namespace SpatialOps {

    template<typename SubExpr, typename FieldType>
    inline SubExpr normalize(NeboExpression<SubExpr, FieldType> arg) {
        return arg.expr();
    }


    inline NeboScalar<Initial, double> normalize(double arg) {
        NeboScalar<Initial, double> typedef ReturnType;

        return ReturnType(arg);
    }

    template<typename FieldType>
    inline NeboConstField<Initial, FieldType> normalize(FieldType const & arg) {
        NeboConstField<Initial, FieldType> typedef ReturnType;

        return ReturnType(arg);
    }


    template<typename Arg>
    struct FinalType {
        NeboConstField<Initial, Arg> typedef Result;
    };

    template<>
    struct FinalType<double> {
        NeboScalar<Initial, double> typedef Result;
    };

    template<typename SubExpr, typename FieldType>
    struct FinalType<NeboExpression<SubExpr, FieldType> > {
        SubExpr typedef Result;
    };


    template<typename ValueType>
    class GenericFieldType {
        typedef ValueType value_type;
    };


    template<typename Arg>
    struct FindFieldType {
        Arg typedef Result;
    };

    template<>
    struct FindFieldType<double> {
        GenericFieldType<double> typedef Result;
    };

    template<typename SubExpr, typename FieldType>
    struct FindFieldType<NeboExpression<SubExpr, FieldType> > {
        FieldType typedef Result;
    };


    // Need something that only compiles if value types match, but returns the (non-generic) field type
    template<typename ValueType1, typename ValueType2, typename FieldType>
    struct ValueTypeCheck;

    template<typename ValueType, typename FieldType>
    struct ValueTypeCheck<ValueType, ValueType, FieldType> {
        FieldType typedef Result;
    };

    template<typename FieldType1, typename FieldType2>
    struct RefineFieldType;

    // Covers Generic + Generic (if matching value_types) and Type + Type
    template<typename FieldType>
    struct RefineFieldType<FieldType, FieldType>
    {
        FieldType typedef Result;
    };

    template<typename ValueType>
    struct RefineFieldType<GenericFieldType<ValueType>, GenericFieldType<ValueType> > {
        GenericFieldType<ValueType> typedef Result;
    };

    template<typename FieldType, typename ValueType>
    struct RefineFieldType<GenericFieldType<ValueType>, FieldType> {
        typename ValueTypeCheck<ValueType, typename FieldType::value_type, FieldType>::Result typedef Result;
    };

    template<typename FieldType, typename ValueType>
    struct RefineFieldType<FieldType, GenericFieldType<ValueType> > {
        typename ValueTypeCheck<ValueType, typename FieldType::value_type, FieldType>::Result typedef Result;
    };

/*    template<typename Functor, typename SubExpr1, typename FieldType1, typename SubExpr2, typename FieldType2, typename SubExpr3, typename FieldType3>
    double test_scalar(NeboExpression<SubExpr1, FieldType1> expr1, NeboExpression<SubExpr2, FieldType2> expr2, NeboExpression<SubExpr3, FieldType3> expr3) {
        typename RefineFieldType3<FieldType1, FieldType2, FieldType2>::Result typedef result_field_type;
        Functor f;
        return f(1, 2, 3);
    }
    */

   } /* SpatialOps */

#endif
/* FULMAR_APPLY_POINTWISE_H */
