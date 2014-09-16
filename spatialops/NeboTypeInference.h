/*
 * Copyright (c) 2014 The University of Utah
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

#ifndef NEBO_TYPE_INFERENCE_H
#define NEBO_TYPE_INFERENCE_H
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
  struct RefineFieldType<FieldType, FieldType> {
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

} /* SpatialOps */

#endif
/* NEBO_TYPE_INFERENCE_H */
