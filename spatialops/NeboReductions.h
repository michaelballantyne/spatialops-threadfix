/* This file was generated by fulmar version 0.9.0. */

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

#ifndef NEBO_REDUCTIONS_H
#  define NEBO_REDUCTIONS_H

#  define field_fold nebo_fold

#  define field_fold_interior nebo_fold_interior

#  define field_reduce nebo_reduce

#  define field_reduce_interior nebo_reduce_interior

#  define field_max nebo_max

#  define field_max_interior nebo_max_interior

#  define field_min nebo_min

#  define field_min_interior nebo_min_interior

#  define field_sum nebo_sum

#  define field_sum_interior nebo_sum_interior

#  define field_norm nebo_norm

#  define field_norm_interior nebo_norm_interior

   namespace SpatialOps {
      template<typename ResultType, typename ExprType, typename FieldType>
       inline ResultType nebo_fold(ResultType const & (*proc)(ResultType const &,
                                                              typename FieldType::
                                                              value_type const &),
                                   ResultType const & initialValue,
                                   NeboExpression<ExprType, FieldType> const &
                                   fexpr) {
          structured::GhostData ghosts = fexpr.expr().possible_ghosts();

          ResultType result = initialValue;

          typename ExprType::ReductionType expr = fexpr.expr().reduce_init();

          while(!(expr.at_end())) {
             result = proc(result, expr.eval());

             expr.next();
          };

          return result;
       };

      template<typename ResultType, typename FieldType>
       inline ResultType nebo_fold(ResultType const & (*proc)(ResultType const &,
                                                              typename FieldType::
                                                              value_type const &),
                                   ResultType const & initialValue,
                                   FieldType const & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_fold(proc,
                           initialValue,
                           NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ResultType, typename ExprType, typename FieldType>
       inline ResultType nebo_fold_interior(ResultType const & (*proc)(ResultType
                                                                       const &,
                                                                       typename
                                                                       FieldType::
                                                                       value_type
                                                                       const &),
                                            ResultType const & initialValue,
                                            NeboExpression<ExprType, FieldType>
                                            const & fexpr) {
          structured::GhostData ghosts(0);

          ResultType result = initialValue;

          typename ExprType::ReductionType expr = fexpr.expr().reduce_init();

          while(!(expr.at_end())) {
             result = proc(result, expr.eval());

             expr.next();
          };

          return result;
       };

      template<typename ResultType, typename FieldType>
       inline ResultType nebo_fold_interior(ResultType const & (*proc)(ResultType
                                                                       const &,
                                                                       typename
                                                                       FieldType::
                                                                       value_type
                                                                       const &),
                                            ResultType const & initialValue,
                                            FieldType const & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_fold_interior(proc,
                                    initialValue,
                                    NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_reduce(typename FieldType::
                                                         value_type const & (*
                                                                             proc)(typename
                                                                                   FieldType::
                                                                                   value_type
                                                                                   const
                                                                                   &,
                                                                                   typename
                                                                                   FieldType::
                                                                                   value_type
                                                                                   const
                                                                                   &),
                                                         NeboExpression<ExprType,
                                                                        FieldType>
                                                         const & fexpr) {
          structured::GhostData ghosts = fexpr.expr().possible_ghosts();

          typename ExprType::ReductionType expr = fexpr.expr().reduce_init();

          typename FieldType::value_type result = expr.eval();

          expr.next();

          while(!(expr.at_end())) {
             result = proc(result, expr.eval());

             expr.next();
          };

          return result;
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_reduce(typename FieldType::
                                                         value_type const & (*
                                                                             proc)(typename
                                                                                   FieldType::
                                                                                   value_type
                                                                                   const
                                                                                   &,
                                                                                   typename
                                                                                   FieldType::
                                                                                   value_type
                                                                                   const
                                                                                   &),
                                                         FieldType const & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_reduce(proc, NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_reduce_interior(typename
                                                                  FieldType::
                                                                  value_type
                                                                  const & (*proc)(typename
                                                                                  FieldType::
                                                                                  value_type
                                                                                  const
                                                                                  &,
                                                                                  typename
                                                                                  FieldType::
                                                                                  value_type
                                                                                  const
                                                                                  &),
                                                                  NeboExpression<ExprType,
                                                                                 FieldType>
                                                                  const & fexpr) {
          structured::GhostData ghosts(0);

          typename ExprType::ReductionType expr = fexpr.expr().reduce_init();

          typename FieldType::value_type result = expr.eval();

          expr.next();

          while(!(expr.at_end())) {
             result = proc(result, expr.eval());

             expr.next();
          };

          return result;
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_reduce_interior(typename
                                                                  FieldType::
                                                                  value_type
                                                                  const & (*proc)(typename
                                                                                  FieldType::
                                                                                  value_type
                                                                                  const
                                                                                  &,
                                                                                  typename
                                                                                  FieldType::
                                                                                  value_type
                                                                                  const
                                                                                  &),
                                                                  FieldType
                                                                  const & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_reduce_interior(proc, NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_reduce(typename FieldType::
                                                         value_type (*proc)(typename
                                                                            FieldType::
                                                                            value_type
                                                                            const
                                                                            &,
                                                                            typename
                                                                            FieldType::
                                                                            value_type
                                                                            const
                                                                            &),
                                                         NeboExpression<ExprType,
                                                                        FieldType>
                                                         const & fexpr) {
          structured::GhostData ghosts = fexpr.expr().possible_ghosts();

          typename ExprType::ReductionType expr = fexpr.expr().reduce_init();

          typename FieldType::value_type result = expr.eval();

          expr.next();

          while(!(expr.at_end())) {
             result = proc(result, expr.eval());

             expr.next();
          };

          return result;
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_reduce(typename FieldType::
                                                         value_type (*proc)(typename
                                                                            FieldType::
                                                                            value_type
                                                                            const
                                                                            &,
                                                                            typename
                                                                            FieldType::
                                                                            value_type
                                                                            const
                                                                            &),
                                                         FieldType const & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_reduce(proc, NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_reduce_interior(typename
                                                                  FieldType::
                                                                  value_type (*
                                                                              proc)(typename
                                                                                    FieldType::
                                                                                    value_type
                                                                                    const
                                                                                    &,
                                                                                    typename
                                                                                    FieldType::
                                                                                    value_type
                                                                                    const
                                                                                    &),
                                                                  NeboExpression<ExprType,
                                                                                 FieldType>
                                                                  const & fexpr) {
          structured::GhostData ghosts(0);

          typename ExprType::ReductionType expr = fexpr.expr().reduce_init();

          typename FieldType::value_type result = expr.eval();

          expr.next();

          while(!(expr.at_end())) {
             result = proc(result, expr.eval());

             expr.next();
          };

          return result;
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_reduce_interior(typename
                                                                  FieldType::
                                                                  value_type (*
                                                                              proc)(typename
                                                                                    FieldType::
                                                                                    value_type
                                                                                    const
                                                                                    &,
                                                                                    typename
                                                                                    FieldType::
                                                                                    value_type
                                                                                    const
                                                                                    &),
                                                                  FieldType
                                                                  const & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_reduce_interior(proc, NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_max(NeboExpression<ExprType,
                                                                     FieldType>
                                                      const & fexpr) {
          return nebo_reduce(std::max, fexpr);
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_max(FieldType const & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_max(NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_max_interior(NeboExpression<ExprType,
                                                                              FieldType>
                                                               const & fexpr) {
          return nebo_reduce_interior(std::max, fexpr);
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_max_interior(FieldType const &
                                                               field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_max_interior(NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_min(NeboExpression<ExprType,
                                                                     FieldType>
                                                      const & fexpr) {
          return nebo_reduce(std::min, fexpr);
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_min(FieldType const & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_min(NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_min_interior(NeboExpression<ExprType,
                                                                              FieldType>
                                                               const & fexpr) {
          return nebo_reduce_interior(std::min, fexpr);
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_min_interior(FieldType const &
                                                               field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_min_interior(NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename AtomicType>
       inline AtomicType sum(AtomicType const & a, AtomicType const & b) {
          return a + b;
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_sum(NeboExpression<ExprType,
                                                                     FieldType>
                                                      const & fexpr) {
          return nebo_reduce(sum<typename FieldType::value_type>, fexpr);
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_sum(FieldType const & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_sum(NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_sum_interior(NeboExpression<ExprType,
                                                                              FieldType>
                                                               const & fexpr) {
          return nebo_reduce_interior(sum<typename FieldType::value_type>, fexpr);
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_sum_interior(FieldType const &
                                                               field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_sum_interior(NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_norm(NeboExpression<ExprType,
                                                                      FieldType>
                                                       const & fexpr) {
          return std::sqrt(nebo_sum(pow(fexpr, 2)));
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_norm(FieldType const & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_norm(NeboExpression<ExprType, FieldType>(ExprType(field)));
       };

      template<typename ExprType, typename FieldType>
       inline typename FieldType::value_type nebo_norm_interior(NeboExpression<ExprType,
                                                                               FieldType>
                                                                const & fexpr) {
          return std::sqrt(nebo_sum_interior(pow(fexpr, 2)));
       };

      template<typename FieldType>
       inline typename FieldType::value_type nebo_norm_interior(FieldType const
                                                                & field) {
          NeboConstField<Initial, FieldType> typedef ExprType;

          return nebo_norm_interior(NeboExpression<ExprType, FieldType>(ExprType(field)));
       };
   } /* SpatialOps */

#endif
/* NEBO_REDUCTIONS_H */
