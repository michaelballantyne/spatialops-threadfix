#ifndef SpatialOps_FieldExpressions_h
#  define SpatialOps_FieldExpressions_h

#  include <spatialops/SpatialOpsConfigure.h>
#  include <spatialops/structured/SpatialField.h>
#  include <cmath>

   /*#include <iostream> */

#  ifdef FIELD_EXPRESSION_THREADS
#     include <vector>
#     include <boost/bind.hpp>
#     include <boost/ref.hpp>
#     include <spatialops/ThreadPool.h>
#     include <spatialops/structured/IntVec.h>
#     include <boost/interprocess/sync/interprocess_semaphore.hpp>
      namespace BI = boost::interprocess;
#  endif //FIELD_EXPRESSION_THREADS

   namespace SpatialOps {

      /* Meta-programming compiler flags */
      struct UseWholeIterator;
      struct UseInteriorIterator;

      template<typename Use, typename FieldType>
       struct IteratorStyle;

      /* UseWholeIterator */
      template<typename FieldType>
       struct IteratorStyle<UseWholeIterator, FieldType> {

         typename FieldType::memory_window typedef MemoryWindow;

         static inline MemoryWindow const & memory_window (FieldType const & field) {
            return field.window_with_ghost();
         };
      };

      /* UseInteriorIterator */
      template<typename FieldType>
       struct IteratorStyle<UseInteriorIterator, FieldType> {

         typename FieldType::memory_window typedef MemoryWindow;

         static inline MemoryWindow const & memory_window (FieldType const & field) {
            return field.window_without_ghost();
         };
      };

      /* Modes: */
      struct Initial;
      template<typename IteratorType, typename SourceType>
       struct ResizePrep;
      template<typename IteratorType, typename SourceType>
       struct Resize;
      template<typename IteratorType, typename SourceType>
       struct SeqWalk;

      template<typename CurrentMode, typename FieldType>
       struct Scalar;

      template<typename FieldType>
       struct Scalar<Initial, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          Scalar<Initial, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             Scalar<ResizePrep<IteratorType, MyType>, FieldType> typedef ResizePrepType;

             Scalar<SeqWalk<IteratorType, MyType>, FieldType> typedef SeqWalkType;
          };
          Scalar(AtomicType const & v)
          : value_(v)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline AtomicType const & value (void) const { return value_; };

         private:
          AtomicType const & value_;
      };

      template<typename IteratorType, typename SourceType, typename FieldType>
       struct Scalar<ResizePrep<IteratorType, SourceType>, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          Scalar<ResizePrep<IteratorType, SourceType>, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          Scalar<Resize<IteratorType, MyType>, FieldType> typedef ResizeType;
          Scalar(SourceType const & source)
          : value_(source.value())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline AtomicType const & value (void) const { return value_; };

         private:
          AtomicType const & value_;
      };

      template<typename IteratorType, typename SourceType, typename FieldType>
       struct Scalar<Resize<IteratorType, SourceType>, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          Scalar<Resize<IteratorType, SourceType>, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          Scalar<SeqWalk<IteratorType, MyType>, FieldType> typedef SeqWalkType;
          Scalar(MemoryWindow const & size, SourceType const & source)
          : value_(source.value())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline AtomicType const & value (void) const { return value_; };

         private:
          AtomicType const & value_;
      };

      template<typename IteratorType, typename SourceType, typename FieldType>
       struct Scalar<SeqWalk<IteratorType, SourceType>, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          Scalar<SeqWalk<IteratorType, SourceType>, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          Scalar(SourceType const & source)
          : value_(source.value())
          {};
          inline void next (void) {};
          inline bool at_end (void) const { return false; };
          inline bool has_length (void) const { return false; };
          inline AtomicType const & eval (void) const { return value_; };

         private:
          AtomicType const & value_;
      };

      template<typename CurrentMode, typename FieldType>
       struct ConstField;

      template<typename FieldType>
       struct ConstField<Initial, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ConstField<Initial, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             ConstField<ResizePrep<IteratorType, MyType>, FieldType> typedef ResizePrepType;

             ConstField<SeqWalk<IteratorType, MyType>, FieldType> typedef SeqWalkType;
          };
          ConstField(FieldType const & f)
          : field_(f)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline FieldType const & field (void) const { return field_; };

         private:
          FieldType const field_;
      };

      template<typename IteratorType, typename SourceType, typename FieldType>
       struct ConstField<ResizePrep<IteratorType, SourceType>, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ConstField<ResizePrep<IteratorType, SourceType>, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          ConstField<Resize<IteratorType, MyType>, FieldType> typedef ResizeType;
          ConstField(SourceType const & source)
          : field_(source.field())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline FieldType const & field (void) const { return field_; };

         private:
          FieldType const field_;
      };

      template<typename IteratorType, typename SourceType, typename FieldType>
       struct ConstField<Resize<IteratorType, SourceType>, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ConstField<Resize<IteratorType, SourceType>, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          ConstField<SeqWalk<IteratorType, MyType>, FieldType> typedef SeqWalkType;
          ConstField(MemoryWindow const & size, SourceType const & source)
          : field_(size, source.field().field_values(), structured::ExternalStorage)
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline FieldType const & field (void) const { return field_; };

         private:
          FieldType const field_;
      };

      template<typename IteratorType, typename SourceType, typename FieldType>
       struct ConstField<SeqWalk<IteratorType, SourceType>, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ConstField<SeqWalk<IteratorType, SourceType>, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          ConstField(SourceType const & source)
          : iter_(source.field().begin()), end_(source.field().end())
          {};
          inline void next (void) { ++iter_; };
          inline bool at_end (void) const { return (iter_ == end_); };
          inline bool has_length (void) const { return true; };
          inline AtomicType const & eval (void) const { return *iter_; };

         private:
          typename FieldType::const_iterator iter_;
          typename FieldType::const_iterator const end_;
      };

      template<typename CurrentMode, typename FieldType>
       struct Field;

      template<typename FieldType>
       struct Field<Initial, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          Field<Initial, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             Field<ResizePrep<IteratorType, MyType>, FieldType> typedef ResizePrepType;

             Field<SeqWalk<IteratorType, MyType>, FieldType> typedef SeqWalkType;
          };
          Field(FieldType & f)
          : field_(f)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline FieldType & field (void) { return field_; };

         private:
          FieldType field_;
      };

      template<typename IteratorType, typename SourceType, typename FieldType>
       struct Field<ResizePrep<IteratorType, SourceType>, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          Field<ResizePrep<IteratorType, SourceType>, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          Field<Resize<IteratorType, MyType>, FieldType> typedef ResizeType;
          Field(SourceType & source)
          : field_(source.field())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) {
             return ResizeType(newSize, *this);
          };
          inline FieldType & field (void) { return field_; };

         private:
          FieldType field_;
      };

      template<typename IteratorType, typename SourceType, typename FieldType>
       struct Field<Resize<IteratorType, SourceType>, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          Field<Resize<IteratorType, SourceType>, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          Field<SeqWalk<IteratorType, MyType>, FieldType> typedef SeqWalkType;
          Field(MemoryWindow const & size, SourceType & source)
          : field_(size, source.field().field_values(), structured::ExternalStorage)
          {};
          inline SeqWalkType init (void) { return SeqWalkType(*this); };
          inline FieldType & field (void) { return field_; };

         private:
          FieldType field_;
      };

      template<typename IteratorType, typename SourceType, typename FieldType>
       struct Field<SeqWalk<IteratorType, SourceType>, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          Field<SeqWalk<IteratorType, SourceType>, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          Field(SourceType & source)
          : iter_(source.field().begin()), end_(source.field().end())
          {};
          inline void next (void) { ++iter_; };
          inline bool at_end (void) const { return (iter_ == end_); };
          inline bool has_length (void) const { return true; };
          inline AtomicType & ref (void) { return *iter_; };
          inline AtomicType * ptr (void) { return &(*iter_); };

         private:
          typename FieldType::iterator iter_;
          typename FieldType::iterator const end_;
      };

      template<typename Operand, typename FieldType>
       struct FieldExpression {

         public:
          FieldType typedef field_type;
          FieldExpression(Operand const & given)
          : expr_(given)
          {};
          inline Operand const & expr (void) const { return expr_; };

         private:
          Operand expr_;
      };

      template<typename Input, typename FieldType>
       struct Standardize;

      template<typename FieldType>
       struct Standardize<FieldType, FieldType> {

         ConstField<Initial, FieldType> typedef StandardType;

         static inline StandardType standardType (FieldType const & given) {
            return StandardType(given);
         };
      };

      template<typename ExprType, typename FieldType>
       struct Standardize<FieldExpression<ExprType, FieldType>, FieldType> {

         ExprType typedef StandardType;

         static inline StandardType standardType (FieldExpression<ExprType, FieldType> const & given) {
            return given.expr();
         };
      };

      template<typename CurrentMode, typename Operand1, typename Operand2, typename FieldType>
       struct SumOp;

      template<typename Operand1, typename Operand2, typename FieldType>
       struct SumOp<Initial, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          SumOp<Initial, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             SumOp<ResizePrep<IteratorType, MyType>,
                   typename Operand1::template Iterator<IteratorType>::ResizePrepType,
                   typename Operand2::template Iterator<IteratorType>::ResizePrepType,
                   FieldType> typedef ResizePrepType;

             SumOp<SeqWalk<IteratorType, MyType>,
                   typename Operand1::template Iterator<IteratorType>::SeqWalkType,
                   typename Operand2::template Iterator<IteratorType>::SeqWalkType,
                   FieldType> typedef SeqWalkType;
          };
          SumOp(Operand1 const & op1, Operand2 const & op2)
          : operand1_(op1), operand2_(op2)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct SumOp<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          SumOp<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          SumOp<Resize<IteratorType, MyType>,
                typename Operand1::ResizeType,
                typename Operand2::ResizeType,
                FieldType> typedef ResizeType;
          SumOp(SourceType const & source)
          : operand1_(source.operand1()), operand2_(source.operand2())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct SumOp<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          SumOp<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          SumOp<SeqWalk<IteratorType, MyType>,
                typename Operand1::SeqWalkType,
                typename Operand2::SeqWalkType,
                FieldType> typedef SeqWalkType;
          SumOp(MemoryWindow const & size, SourceType const & source)
          : operand1_(size, source.operand1()), operand2_(size, source.operand2())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct SumOp<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          SumOp<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          SumOp(SourceType const & source)
          : operand1_(source.operand1()), operand2_(source.operand2())
          {};
          inline void next (void) { operand1_.next(); operand2_.next(); };
          inline bool at_end (void) const { return (operand1_.at_end() || operand2_.at_end()); };
          inline bool has_length (void) const {
             return (operand1_.has_length() || operand2_.has_length());
          };
          inline AtomicType eval (void) const { return (operand1_.eval() + operand2_.eval()); };

         private:
          Operand1 operand1_;
          Operand2 operand2_;
      };

      /* SubExpr X SubExpr */
      template<typename SubExpr1, typename SubExpr2>
       inline FieldExpression<SumOp<Initial,
                                    typename Standardize<SubExpr1, typename SubExpr1::field_type>::
                                    StandardType,
                                    typename Standardize<SubExpr2, typename SubExpr1::field_type>::
                                    StandardType,
                                    typename SubExpr1::field_type>,
                              typename SubExpr1::field_type> operator + (SubExpr1 const & arg1,
                                                                         SubExpr2 const & arg2) {

          typename SubExpr1::field_type typedef FieldType;

          typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef Type1;

          typename Standardize<SubExpr2, typename SubExpr1::field_type>::StandardType typedef Type2;

          SumOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)),
                                       Type2(Standardize<SubExpr2, FieldType>::standardType(arg2))));
       };

      /* SubExpr X Scalar */
      template<typename SubExpr1>
       inline FieldExpression<SumOp<Initial,
                                    typename Standardize<SubExpr1, typename SubExpr1::field_type>::
                                    StandardType,
                                    Scalar<Initial, typename SubExpr1::field_type>,
                                    typename SubExpr1::field_type>,
                              typename SubExpr1::field_type> operator + (SubExpr1 const & arg1,
                                                                         typename SubExpr1::
                                                                         field_type::value_type
                                                                         const & arg2) {

          typename SubExpr1::field_type typedef FieldType;

          typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef Type1;

          Scalar<Initial, typename SubExpr1::field_type> typedef Type2;

          SumOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)),
                                       Type2(Type2(arg2))));
       };

      /* Scalar X SubExpr */
      template<typename SubExpr2>
       inline FieldExpression<SumOp<Initial,
                                    Scalar<Initial, typename SubExpr2::field_type>,
                                    typename Standardize<SubExpr2, typename SubExpr2::field_type>::
                                    StandardType,
                                    typename SubExpr2::field_type>,
                              typename SubExpr2::field_type> operator + (typename SubExpr2::
                                                                         field_type::value_type
                                                                         const & arg1,
                                                                         SubExpr2 const & arg2) {

          typename SubExpr2::field_type typedef FieldType;

          Scalar<Initial, typename SubExpr2::field_type> typedef Type1;

          typename Standardize<SubExpr2, typename SubExpr2::field_type>::StandardType typedef Type2;

          SumOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Type1(arg1)),
                                       Type2(Standardize<SubExpr2, FieldType>::standardType(arg2))));
       };

      template<typename CurrentMode, typename Operand1, typename Operand2, typename FieldType>
       struct DiffOp;

      template<typename Operand1, typename Operand2, typename FieldType>
       struct DiffOp<Initial, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          DiffOp<Initial, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             DiffOp<ResizePrep<IteratorType, MyType>,
                    typename Operand1::template Iterator<IteratorType>::ResizePrepType,
                    typename Operand2::template Iterator<IteratorType>::ResizePrepType,
                    FieldType> typedef ResizePrepType;

             DiffOp<SeqWalk<IteratorType, MyType>,
                    typename Operand1::template Iterator<IteratorType>::SeqWalkType,
                    typename Operand2::template Iterator<IteratorType>::SeqWalkType,
                    FieldType> typedef SeqWalkType;
          };
          DiffOp(Operand1 const & op1, Operand2 const & op2)
          : operand1_(op1), operand2_(op2)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct DiffOp<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          DiffOp<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          DiffOp<Resize<IteratorType, MyType>,
                 typename Operand1::ResizeType,
                 typename Operand2::ResizeType,
                 FieldType> typedef ResizeType;
          DiffOp(SourceType const & source)
          : operand1_(source.operand1()), operand2_(source.operand2())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct DiffOp<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          DiffOp<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          DiffOp<SeqWalk<IteratorType, MyType>,
                 typename Operand1::SeqWalkType,
                 typename Operand2::SeqWalkType,
                 FieldType> typedef SeqWalkType;
          DiffOp(MemoryWindow const & size, SourceType const & source)
          : operand1_(size, source.operand1()), operand2_(size, source.operand2())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct DiffOp<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          DiffOp<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          DiffOp(SourceType const & source)
          : operand1_(source.operand1()), operand2_(source.operand2())
          {};
          inline void next (void) { operand1_.next(); operand2_.next(); };
          inline bool at_end (void) const { return (operand1_.at_end() || operand2_.at_end()); };
          inline bool has_length (void) const {
             return (operand1_.has_length() || operand2_.has_length());
          };
          inline AtomicType eval (void) const { return (operand1_.eval() - operand2_.eval()); };

         private:
          Operand1 operand1_;
          Operand2 operand2_;
      };

      /* SubExpr X SubExpr */
      template<typename SubExpr1, typename SubExpr2>
       inline FieldExpression<DiffOp<Initial,
                                     typename Standardize<SubExpr1, typename SubExpr1::field_type>::
                                     StandardType,
                                     typename Standardize<SubExpr2, typename SubExpr1::field_type>::
                                     StandardType,
                                     typename SubExpr1::field_type>,
                              typename SubExpr1::field_type> operator - (SubExpr1 const & arg1,
                                                                         SubExpr2 const & arg2) {

          typename SubExpr1::field_type typedef FieldType;

          typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef Type1;

          typename Standardize<SubExpr2, typename SubExpr1::field_type>::StandardType typedef Type2;

          DiffOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)),
                                       Type2(Standardize<SubExpr2, FieldType>::standardType(arg2))));
       };

      /* SubExpr X Scalar */
      template<typename SubExpr1>
       inline FieldExpression<DiffOp<Initial,
                                     typename Standardize<SubExpr1, typename SubExpr1::field_type>::
                                     StandardType,
                                     Scalar<Initial, typename SubExpr1::field_type>,
                                     typename SubExpr1::field_type>,
                              typename SubExpr1::field_type> operator - (SubExpr1 const & arg1,
                                                                         typename SubExpr1::
                                                                         field_type::value_type
                                                                         const & arg2) {

          typename SubExpr1::field_type typedef FieldType;

          typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef Type1;

          Scalar<Initial, typename SubExpr1::field_type> typedef Type2;

          DiffOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)),
                                       Type2(Type2(arg2))));
       };

      /* Scalar X SubExpr */
      template<typename SubExpr2>
       inline FieldExpression<DiffOp<Initial,
                                     Scalar<Initial, typename SubExpr2::field_type>,
                                     typename Standardize<SubExpr2, typename SubExpr2::field_type>::
                                     StandardType,
                                     typename SubExpr2::field_type>,
                              typename SubExpr2::field_type> operator - (typename SubExpr2::
                                                                         field_type::value_type
                                                                         const & arg1,
                                                                         SubExpr2 const & arg2) {

          typename SubExpr2::field_type typedef FieldType;

          Scalar<Initial, typename SubExpr2::field_type> typedef Type1;

          typename Standardize<SubExpr2, typename SubExpr2::field_type>::StandardType typedef Type2;

          DiffOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Type1(arg1)),
                                       Type2(Standardize<SubExpr2, FieldType>::standardType(arg2))));
       };

      template<typename CurrentMode, typename Operand1, typename Operand2, typename FieldType>
       struct ProdOp;

      template<typename Operand1, typename Operand2, typename FieldType>
       struct ProdOp<Initial, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ProdOp<Initial, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             ProdOp<ResizePrep<IteratorType, MyType>,
                    typename Operand1::template Iterator<IteratorType>::ResizePrepType,
                    typename Operand2::template Iterator<IteratorType>::ResizePrepType,
                    FieldType> typedef ResizePrepType;

             ProdOp<SeqWalk<IteratorType, MyType>,
                    typename Operand1::template Iterator<IteratorType>::SeqWalkType,
                    typename Operand2::template Iterator<IteratorType>::SeqWalkType,
                    FieldType> typedef SeqWalkType;
          };
          ProdOp(Operand1 const & op1, Operand2 const & op2)
          : operand1_(op1), operand2_(op2)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct ProdOp<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ProdOp<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          ProdOp<Resize<IteratorType, MyType>,
                 typename Operand1::ResizeType,
                 typename Operand2::ResizeType,
                 FieldType> typedef ResizeType;
          ProdOp(SourceType const & source)
          : operand1_(source.operand1()), operand2_(source.operand2())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct ProdOp<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ProdOp<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          ProdOp<SeqWalk<IteratorType, MyType>,
                 typename Operand1::SeqWalkType,
                 typename Operand2::SeqWalkType,
                 FieldType> typedef SeqWalkType;
          ProdOp(MemoryWindow const & size, SourceType const & source)
          : operand1_(size, source.operand1()), operand2_(size, source.operand2())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct ProdOp<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ProdOp<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          ProdOp(SourceType const & source)
          : operand1_(source.operand1()), operand2_(source.operand2())
          {};
          inline void next (void) { operand1_.next(); operand2_.next(); };
          inline bool at_end (void) const { return (operand1_.at_end() || operand2_.at_end()); };
          inline bool has_length (void) const {
             return (operand1_.has_length() || operand2_.has_length());
          };
          inline AtomicType eval (void) const { return (operand1_.eval() * operand2_.eval()); };

         private:
          Operand1 operand1_;
          Operand2 operand2_;
      };

      /* SubExpr X SubExpr */
      template<typename SubExpr1, typename SubExpr2>
       inline FieldExpression<ProdOp<Initial,
                                     typename Standardize<SubExpr1, typename SubExpr1::field_type>::
                                     StandardType,
                                     typename Standardize<SubExpr2, typename SubExpr1::field_type>::
                                     StandardType,
                                     typename SubExpr1::field_type>,
                              typename SubExpr1::field_type> operator * (SubExpr1 const & arg1,
                                                                         SubExpr2 const & arg2) {

          typename SubExpr1::field_type typedef FieldType;

          typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef Type1;

          typename Standardize<SubExpr2, typename SubExpr1::field_type>::StandardType typedef Type2;

          ProdOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)),
                                       Type2(Standardize<SubExpr2, FieldType>::standardType(arg2))));
       };

      /* SubExpr X Scalar */
      template<typename SubExpr1>
       inline FieldExpression<ProdOp<Initial,
                                     typename Standardize<SubExpr1, typename SubExpr1::field_type>::
                                     StandardType,
                                     Scalar<Initial, typename SubExpr1::field_type>,
                                     typename SubExpr1::field_type>,
                              typename SubExpr1::field_type> operator * (SubExpr1 const & arg1,
                                                                         typename SubExpr1::
                                                                         field_type::value_type
                                                                         const & arg2) {

          typename SubExpr1::field_type typedef FieldType;

          typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef Type1;

          Scalar<Initial, typename SubExpr1::field_type> typedef Type2;

          ProdOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)),
                                       Type2(Type2(arg2))));
       };

      /* Scalar X SubExpr */
      template<typename SubExpr2>
       inline FieldExpression<ProdOp<Initial,
                                     Scalar<Initial, typename SubExpr2::field_type>,
                                     typename Standardize<SubExpr2, typename SubExpr2::field_type>::
                                     StandardType,
                                     typename SubExpr2::field_type>,
                              typename SubExpr2::field_type> operator * (typename SubExpr2::
                                                                         field_type::value_type
                                                                         const & arg1,
                                                                         SubExpr2 const & arg2) {

          typename SubExpr2::field_type typedef FieldType;

          Scalar<Initial, typename SubExpr2::field_type> typedef Type1;

          typename Standardize<SubExpr2, typename SubExpr2::field_type>::StandardType typedef Type2;

          ProdOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Type1(arg1)),
                                       Type2(Standardize<SubExpr2, FieldType>::standardType(arg2))));
       };

      template<typename CurrentMode, typename Operand1, typename Operand2, typename FieldType>
       struct DivOp;

      template<typename Operand1, typename Operand2, typename FieldType>
       struct DivOp<Initial, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          DivOp<Initial, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             DivOp<ResizePrep<IteratorType, MyType>,
                   typename Operand1::template Iterator<IteratorType>::ResizePrepType,
                   typename Operand2::template Iterator<IteratorType>::ResizePrepType,
                   FieldType> typedef ResizePrepType;

             DivOp<SeqWalk<IteratorType, MyType>,
                   typename Operand1::template Iterator<IteratorType>::SeqWalkType,
                   typename Operand2::template Iterator<IteratorType>::SeqWalkType,
                   FieldType> typedef SeqWalkType;
          };
          DivOp(Operand1 const & op1, Operand2 const & op2)
          : operand1_(op1), operand2_(op2)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct DivOp<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          DivOp<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          DivOp<Resize<IteratorType, MyType>,
                typename Operand1::ResizeType,
                typename Operand2::ResizeType,
                FieldType> typedef ResizeType;
          DivOp(SourceType const & source)
          : operand1_(source.operand1()), operand2_(source.operand2())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct DivOp<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          DivOp<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          DivOp<SeqWalk<IteratorType, MyType>,
                typename Operand1::SeqWalkType,
                typename Operand2::SeqWalkType,
                FieldType> typedef SeqWalkType;
          DivOp(MemoryWindow const & size, SourceType const & source)
          : operand1_(size, source.operand1()), operand2_(size, source.operand2())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct DivOp<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          DivOp<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          DivOp(SourceType const & source)
          : operand1_(source.operand1()), operand2_(source.operand2())
          {};
          inline void next (void) { operand1_.next(); operand2_.next(); };
          inline bool at_end (void) const { return (operand1_.at_end() || operand2_.at_end()); };
          inline bool has_length (void) const {
             return (operand1_.has_length() || operand2_.has_length());
          };
          inline AtomicType eval (void) const { return (operand1_.eval() / operand2_.eval()); };

         private:
          Operand1 operand1_;
          Operand2 operand2_;
      };

      /* SubExpr X SubExpr */
      template<typename SubExpr1, typename SubExpr2>
       inline FieldExpression<DivOp<Initial,
                                    typename Standardize<SubExpr1, typename SubExpr1::field_type>::
                                    StandardType,
                                    typename Standardize<SubExpr2, typename SubExpr1::field_type>::
                                    StandardType,
                                    typename SubExpr1::field_type>,
                              typename SubExpr1::field_type> operator / (SubExpr1 const & arg1,
                                                                         SubExpr2 const & arg2) {

          typename SubExpr1::field_type typedef FieldType;

          typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef Type1;

          typename Standardize<SubExpr2, typename SubExpr1::field_type>::StandardType typedef Type2;

          DivOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)),
                                       Type2(Standardize<SubExpr2, FieldType>::standardType(arg2))));
       };

      /* SubExpr X Scalar */
      template<typename SubExpr1>
       inline FieldExpression<DivOp<Initial,
                                    typename Standardize<SubExpr1, typename SubExpr1::field_type>::
                                    StandardType,
                                    Scalar<Initial, typename SubExpr1::field_type>,
                                    typename SubExpr1::field_type>,
                              typename SubExpr1::field_type> operator / (SubExpr1 const & arg1,
                                                                         typename SubExpr1::
                                                                         field_type::value_type
                                                                         const & arg2) {

          typename SubExpr1::field_type typedef FieldType;

          typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef Type1;

          Scalar<Initial, typename SubExpr1::field_type> typedef Type2;

          DivOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)),
                                       Type2(Type2(arg2))));
       };

      /* Scalar X SubExpr */
      template<typename SubExpr2>
       inline FieldExpression<DivOp<Initial,
                                    Scalar<Initial, typename SubExpr2::field_type>,
                                    typename Standardize<SubExpr2, typename SubExpr2::field_type>::
                                    StandardType,
                                    typename SubExpr2::field_type>,
                              typename SubExpr2::field_type> operator / (typename SubExpr2::
                                                                         field_type::value_type
                                                                         const & arg1,
                                                                         SubExpr2 const & arg2) {

          typename SubExpr2::field_type typedef FieldType;

          Scalar<Initial, typename SubExpr2::field_type> typedef Type1;

          typename Standardize<SubExpr2, typename SubExpr2::field_type>::StandardType typedef Type2;

          DivOp<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Type1(arg1)),
                                       Type2(Standardize<SubExpr2, FieldType>::standardType(arg2))));
       };

      template<typename CurrentMode, typename Operand, typename FieldType>
       struct SinFcn;

      template<typename Operand, typename FieldType>
       struct SinFcn<Initial, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          SinFcn<Initial, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             SinFcn<ResizePrep<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::ResizePrepType,
                    FieldType> typedef ResizePrepType;

             SinFcn<SeqWalk<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::SeqWalkType,
                    FieldType> typedef SeqWalkType;
          };
          SinFcn(Operand const & op)
          : operand_(op)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct SinFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          SinFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          SinFcn<Resize<IteratorType, MyType>, typename Operand::ResizeType, FieldType> typedef
          ResizeType;
          SinFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct SinFcn<Resize<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          SinFcn<Resize<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          SinFcn<SeqWalk<IteratorType, MyType>, typename Operand::SeqWalkType, FieldType> typedef
          SeqWalkType;
          SinFcn(MemoryWindow const & size, SourceType const & source)
          : operand_(size, source.operand())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct SinFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          SinFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          SinFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline void next (void) { operand_.next(); };
          inline bool at_end (void) const { return (operand_.at_end()); };
          inline bool has_length (void) const { return (operand_.has_length()); };
          inline AtomicType eval (void) const { return std::sin(operand_.eval()); };

         private:
          Operand operand_;
      };

      /* SubExpr */
      template<typename SubExpr>
       inline FieldExpression<SinFcn<Initial,
                                     typename Standardize<SubExpr, typename SubExpr::field_type>::
                                     StandardType,
                                     typename SubExpr::field_type>,
                              typename SubExpr::field_type> sin (SubExpr const & arg) {

          typename SubExpr::field_type typedef FieldType;

          typename Standardize<SubExpr, typename SubExpr::field_type>::StandardType typedef Type;

          SinFcn<Initial, Type, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type(Standardize<SubExpr, FieldType>::standardType(arg))));
       };

      template<typename CurrentMode, typename Operand, typename FieldType>
       struct CosFcn;

      template<typename Operand, typename FieldType>
       struct CosFcn<Initial, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          CosFcn<Initial, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             CosFcn<ResizePrep<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::ResizePrepType,
                    FieldType> typedef ResizePrepType;

             CosFcn<SeqWalk<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::SeqWalkType,
                    FieldType> typedef SeqWalkType;
          };
          CosFcn(Operand const & op)
          : operand_(op)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct CosFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          CosFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          CosFcn<Resize<IteratorType, MyType>, typename Operand::ResizeType, FieldType> typedef
          ResizeType;
          CosFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct CosFcn<Resize<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          CosFcn<Resize<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          CosFcn<SeqWalk<IteratorType, MyType>, typename Operand::SeqWalkType, FieldType> typedef
          SeqWalkType;
          CosFcn(MemoryWindow const & size, SourceType const & source)
          : operand_(size, source.operand())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct CosFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          CosFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          CosFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline void next (void) { operand_.next(); };
          inline bool at_end (void) const { return (operand_.at_end()); };
          inline bool has_length (void) const { return (operand_.has_length()); };
          inline AtomicType eval (void) const { return std::cos(operand_.eval()); };

         private:
          Operand operand_;
      };

      /* SubExpr */
      template<typename SubExpr>
       inline FieldExpression<CosFcn<Initial,
                                     typename Standardize<SubExpr, typename SubExpr::field_type>::
                                     StandardType,
                                     typename SubExpr::field_type>,
                              typename SubExpr::field_type> cos (SubExpr const & arg) {

          typename SubExpr::field_type typedef FieldType;

          typename Standardize<SubExpr, typename SubExpr::field_type>::StandardType typedef Type;

          CosFcn<Initial, Type, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type(Standardize<SubExpr, FieldType>::standardType(arg))));
       };

      template<typename CurrentMode, typename Operand, typename FieldType>
       struct TanFcn;

      template<typename Operand, typename FieldType>
       struct TanFcn<Initial, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          TanFcn<Initial, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             TanFcn<ResizePrep<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::ResizePrepType,
                    FieldType> typedef ResizePrepType;

             TanFcn<SeqWalk<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::SeqWalkType,
                    FieldType> typedef SeqWalkType;
          };
          TanFcn(Operand const & op)
          : operand_(op)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct TanFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          TanFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          TanFcn<Resize<IteratorType, MyType>, typename Operand::ResizeType, FieldType> typedef
          ResizeType;
          TanFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct TanFcn<Resize<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          TanFcn<Resize<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          TanFcn<SeqWalk<IteratorType, MyType>, typename Operand::SeqWalkType, FieldType> typedef
          SeqWalkType;
          TanFcn(MemoryWindow const & size, SourceType const & source)
          : operand_(size, source.operand())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct TanFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          TanFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          TanFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline void next (void) { operand_.next(); };
          inline bool at_end (void) const { return (operand_.at_end()); };
          inline bool has_length (void) const { return (operand_.has_length()); };
          inline AtomicType eval (void) const { return std::tan(operand_.eval()); };

         private:
          Operand operand_;
      };

      /* SubExpr */
      template<typename SubExpr>
       inline FieldExpression<TanFcn<Initial,
                                     typename Standardize<SubExpr, typename SubExpr::field_type>::
                                     StandardType,
                                     typename SubExpr::field_type>,
                              typename SubExpr::field_type> tan (SubExpr const & arg) {

          typename SubExpr::field_type typedef FieldType;

          typename Standardize<SubExpr, typename SubExpr::field_type>::StandardType typedef Type;

          TanFcn<Initial, Type, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type(Standardize<SubExpr, FieldType>::standardType(arg))));
       };

      template<typename CurrentMode, typename Operand, typename FieldType>
       struct ExpFcn;

      template<typename Operand, typename FieldType>
       struct ExpFcn<Initial, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ExpFcn<Initial, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             ExpFcn<ResizePrep<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::ResizePrepType,
                    FieldType> typedef ResizePrepType;

             ExpFcn<SeqWalk<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::SeqWalkType,
                    FieldType> typedef SeqWalkType;
          };
          ExpFcn(Operand const & op)
          : operand_(op)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct ExpFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ExpFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          ExpFcn<Resize<IteratorType, MyType>, typename Operand::ResizeType, FieldType> typedef
          ResizeType;
          ExpFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct ExpFcn<Resize<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ExpFcn<Resize<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          ExpFcn<SeqWalk<IteratorType, MyType>, typename Operand::SeqWalkType, FieldType> typedef
          SeqWalkType;
          ExpFcn(MemoryWindow const & size, SourceType const & source)
          : operand_(size, source.operand())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct ExpFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          ExpFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          ExpFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline void next (void) { operand_.next(); };
          inline bool at_end (void) const { return (operand_.at_end()); };
          inline bool has_length (void) const { return (operand_.has_length()); };
          inline AtomicType eval (void) const { return std::exp(operand_.eval()); };

         private:
          Operand operand_;
      };

      /* SubExpr */
      template<typename SubExpr>
       inline FieldExpression<ExpFcn<Initial,
                                     typename Standardize<SubExpr, typename SubExpr::field_type>::
                                     StandardType,
                                     typename SubExpr::field_type>,
                              typename SubExpr::field_type> exp (SubExpr const & arg) {

          typename SubExpr::field_type typedef FieldType;

          typename Standardize<SubExpr, typename SubExpr::field_type>::StandardType typedef Type;

          ExpFcn<Initial, Type, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type(Standardize<SubExpr, FieldType>::standardType(arg))));
       };

      template<typename CurrentMode, typename Operand, typename FieldType>
       struct TanhFcn;

      template<typename Operand, typename FieldType>
       struct TanhFcn<Initial, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          TanhFcn<Initial, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             TanhFcn<ResizePrep<IteratorType, MyType>,
                     typename Operand::template Iterator<IteratorType>::ResizePrepType,
                     FieldType> typedef ResizePrepType;

             TanhFcn<SeqWalk<IteratorType, MyType>,
                     typename Operand::template Iterator<IteratorType>::SeqWalkType,
                     FieldType> typedef SeqWalkType;
          };
          TanhFcn(Operand const & op)
          : operand_(op)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct TanhFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          TanhFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          TanhFcn<Resize<IteratorType, MyType>, typename Operand::ResizeType, FieldType> typedef
          ResizeType;
          TanhFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct TanhFcn<Resize<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          TanhFcn<Resize<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          TanhFcn<SeqWalk<IteratorType, MyType>, typename Operand::SeqWalkType, FieldType> typedef
          SeqWalkType;
          TanhFcn(MemoryWindow const & size, SourceType const & source)
          : operand_(size, source.operand())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct TanhFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          TanhFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          TanhFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline void next (void) { operand_.next(); };
          inline bool at_end (void) const { return (operand_.at_end()); };
          inline bool has_length (void) const { return (operand_.has_length()); };
          inline AtomicType eval (void) const { return std::tanh(operand_.eval()); };

         private:
          Operand operand_;
      };

      /* SubExpr */
      template<typename SubExpr>
       inline FieldExpression<TanhFcn<Initial,
                                      typename Standardize<SubExpr, typename SubExpr::field_type>::
                                      StandardType,
                                      typename SubExpr::field_type>,
                              typename SubExpr::field_type> tanh (SubExpr const & arg) {

          typename SubExpr::field_type typedef FieldType;

          typename Standardize<SubExpr, typename SubExpr::field_type>::StandardType typedef Type;

          TanhFcn<Initial, Type, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type(Standardize<SubExpr, FieldType>::standardType(arg))));
       };

      template<typename CurrentMode, typename Operand, typename FieldType>
       struct AbsFcn;

      template<typename Operand, typename FieldType>
       struct AbsFcn<Initial, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          AbsFcn<Initial, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             AbsFcn<ResizePrep<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::ResizePrepType,
                    FieldType> typedef ResizePrepType;

             AbsFcn<SeqWalk<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::SeqWalkType,
                    FieldType> typedef SeqWalkType;
          };
          AbsFcn(Operand const & op)
          : operand_(op)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct AbsFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          AbsFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          AbsFcn<Resize<IteratorType, MyType>, typename Operand::ResizeType, FieldType> typedef
          ResizeType;
          AbsFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct AbsFcn<Resize<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          AbsFcn<Resize<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          AbsFcn<SeqWalk<IteratorType, MyType>, typename Operand::SeqWalkType, FieldType> typedef
          SeqWalkType;
          AbsFcn(MemoryWindow const & size, SourceType const & source)
          : operand_(size, source.operand())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct AbsFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          AbsFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          AbsFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline void next (void) { operand_.next(); };
          inline bool at_end (void) const { return (operand_.at_end()); };
          inline bool has_length (void) const { return (operand_.has_length()); };
          inline AtomicType eval (void) const { return std::abs(operand_.eval()); };

         private:
          Operand operand_;
      };

      /* SubExpr */
      template<typename SubExpr>
       inline FieldExpression<AbsFcn<Initial,
                                     typename Standardize<SubExpr, typename SubExpr::field_type>::
                                     StandardType,
                                     typename SubExpr::field_type>,
                              typename SubExpr::field_type> abs (SubExpr const & arg) {

          typename SubExpr::field_type typedef FieldType;

          typename Standardize<SubExpr, typename SubExpr::field_type>::StandardType typedef Type;

          AbsFcn<Initial, Type, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type(Standardize<SubExpr, FieldType>::standardType(arg))));
       };

      template<typename CurrentMode, typename Operand, typename FieldType>
       struct NegFcn;

      template<typename Operand, typename FieldType>
       struct NegFcn<Initial, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          NegFcn<Initial, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             NegFcn<ResizePrep<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::ResizePrepType,
                    FieldType> typedef ResizePrepType;

             NegFcn<SeqWalk<IteratorType, MyType>,
                    typename Operand::template Iterator<IteratorType>::SeqWalkType,
                    FieldType> typedef SeqWalkType;
          };
          NegFcn(Operand const & op)
          : operand_(op)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct NegFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          NegFcn<ResizePrep<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          NegFcn<Resize<IteratorType, MyType>, typename Operand::ResizeType, FieldType> typedef
          ResizeType;
          NegFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct NegFcn<Resize<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          NegFcn<Resize<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          NegFcn<SeqWalk<IteratorType, MyType>, typename Operand::SeqWalkType, FieldType> typedef
          SeqWalkType;
          NegFcn(MemoryWindow const & size, SourceType const & source)
          : operand_(size, source.operand())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand const & operand (void) const { return operand_; };

         private:
          Operand const operand_;
      };

      template<typename IteratorType, typename SourceType, typename Operand, typename FieldType>
       struct NegFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          NegFcn<SeqWalk<IteratorType, SourceType>, Operand, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          NegFcn(SourceType const & source)
          : operand_(source.operand())
          {};
          inline void next (void) { operand_.next(); };
          inline bool at_end (void) const { return (operand_.at_end()); };
          inline bool has_length (void) const { return (operand_.has_length()); };
          inline AtomicType eval (void) const { return -(operand_.eval()); };

         private:
          Operand operand_;
      };

      /* SubExpr */
      template<typename SubExpr>
       inline FieldExpression<NegFcn<Initial,
                                     typename Standardize<SubExpr, typename SubExpr::field_type>::
                                     StandardType,
                                     typename SubExpr::field_type>,
                              typename SubExpr::field_type> operator - (SubExpr const & arg) {

          typename SubExpr::field_type typedef FieldType;

          typename Standardize<SubExpr, typename SubExpr::field_type>::StandardType typedef Type;

          NegFcn<Initial, Type, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type(Standardize<SubExpr, FieldType>::standardType(arg))));
       };

      template<typename CurrentMode, typename Operand1, typename Operand2, typename FieldType>
       struct PowFcn;

      template<typename Operand1, typename Operand2, typename FieldType>
       struct PowFcn<Initial, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          PowFcn<Initial, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          template<typename IteratorType>
           struct Iterator {

             PowFcn<ResizePrep<IteratorType, MyType>,
                    typename Operand1::template Iterator<IteratorType>::ResizePrepType,
                    typename Operand2::template Iterator<IteratorType>::ResizePrepType,
                    FieldType> typedef ResizePrepType;

             PowFcn<SeqWalk<IteratorType, MyType>,
                    typename Operand1::template Iterator<IteratorType>::SeqWalkType,
                    typename Operand2::template Iterator<IteratorType>::SeqWalkType,
                    FieldType> typedef SeqWalkType;
          };
          PowFcn(Operand1 const & op1, Operand2 const & op2)
          : operand1_(op1), operand2_(op2)
          {};
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::SeqWalkType init (void) const {
              return typename Iterator<IteratorType>::SeqWalkType(*this);
           };
          template<typename IteratorType>
           inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {
              return typename Iterator<IteratorType>::ResizePrepType (*this);
           };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct PowFcn<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          PowFcn<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          PowFcn<Resize<IteratorType, MyType>,
                 typename Operand1::ResizeType,
                 typename Operand2::ResizeType,
                 FieldType> typedef ResizeType;
          PowFcn(SourceType const & source)
          : operand1_(source.operand1()), operand2_(source.operand2())
          {};
          inline ResizeType resize ( MemoryWindow const & newSize) const {
             return ResizeType(newSize, *this);
          };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct PowFcn<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          PowFcn<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          PowFcn<SeqWalk<IteratorType, MyType>,
                 typename Operand1::SeqWalkType,
                 typename Operand2::SeqWalkType,
                 FieldType> typedef SeqWalkType;
          PowFcn(MemoryWindow const & size, SourceType const & source)
          : operand1_(size, source.operand1()), operand2_(size, source.operand2())
          {};
          inline SeqWalkType init (void) const { return SeqWalkType(*this); };
          inline Operand1 const & operand1 (void) const { return operand1_; };
          inline Operand2 const & operand2 (void) const { return operand2_; };

         private:
          Operand1 const operand1_;
          Operand2 const operand2_;
      };

      template<typename IteratorType,
               typename SourceType,
               typename Operand1,
               typename Operand2,
               typename FieldType>
       struct PowFcn<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> {

         public:
          typename FieldType::value_type typedef AtomicType;
          PowFcn<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef MyType;
          typename FieldType::memory_window typedef MemoryWindow;
          PowFcn(SourceType const & source)
          : operand1_(source.operand1()), operand2_(source.operand2())
          {};
          inline void next (void) { operand1_.next(); operand2_.next(); };
          inline bool at_end (void) const { return (operand1_.at_end() || operand2_.at_end()); };
          inline bool has_length (void) const {
             return (operand1_.has_length() || operand2_.has_length());
          };
          inline AtomicType eval (void) const {
             return std::pow(operand1_.eval(), operand2_.eval());
          };

         private:
          Operand1 operand1_;
          Operand2 operand2_;
      };

      /* SubExpr X SubExpr */
      template<typename SubExpr1, typename SubExpr2>
       inline FieldExpression<PowFcn<Initial,
                                     typename Standardize<SubExpr1, typename SubExpr1::field_type>::
                                     StandardType,
                                     typename Standardize<SubExpr2, typename SubExpr1::field_type>::
                                     StandardType,
                                     typename SubExpr1::field_type>,
                              typename SubExpr1::field_type> pow (SubExpr1 const & arg1,
                                                                  SubExpr2 const & arg2) {

          typename SubExpr1::field_type typedef FieldType;

          typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef Type1;

          typename Standardize<SubExpr2, typename SubExpr1::field_type>::StandardType typedef Type2;

          PowFcn<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)),
                                       Type2(Standardize<SubExpr2, FieldType>::standardType(arg2))));
       };

      /* SubExpr X Scalar */
      template<typename SubExpr1>
       inline FieldExpression<PowFcn<Initial,
                                     typename Standardize<SubExpr1, typename SubExpr1::field_type>::
                                     StandardType,
                                     Scalar<Initial, typename SubExpr1::field_type>,
                                     typename SubExpr1::field_type>,
                              typename SubExpr1::field_type> pow (SubExpr1 const & arg1,
                                                                  typename SubExpr1::field_type::
                                                                  value_type const & arg2) {

          typename SubExpr1::field_type typedef FieldType;

          typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef Type1;

          Scalar<Initial, typename SubExpr1::field_type> typedef Type2;

          PowFcn<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)),
                                       Type2(Type2(arg2))));
       };

      /* Scalar X SubExpr */
      template<typename SubExpr2>
       inline FieldExpression<PowFcn<Initial,
                                     Scalar<Initial, typename SubExpr2::field_type>,
                                     typename Standardize<SubExpr2, typename SubExpr2::field_type>::
                                     StandardType,
                                     typename SubExpr2::field_type>,
                              typename SubExpr2::field_type> pow (typename SubExpr2::field_type::
                                                                  value_type const & arg1,
                                                                  SubExpr2 const & arg2) {

          typename SubExpr2::field_type typedef FieldType;

          Scalar<Initial, typename SubExpr2::field_type> typedef Type1;

          typename Standardize<SubExpr2, typename SubExpr2::field_type>::StandardType typedef Type2;

          PowFcn<Initial, Type1, Type2, FieldType> typedef ReturnType;

          FieldExpression<ReturnType, FieldType> typedef ReturnTerm;

          return ReturnTerm(ReturnType(Type1(Type1(arg1)),
                                       Type2(Standardize<SubExpr2, FieldType>::standardType(arg2))));
       };

#     define BUILD_BINARY_FUNCTION(OBJECT_NAME, INTERNAL_NAME, EXTERNAL_NAME)                      \
         template<typename CurrentMode, typename Operand1, typename Operand2, typename FieldType>  \
          struct OBJECT_NAME;                                                                      \
                                                                                                   \
         template<typename Operand1, typename Operand2, typename FieldType>                        \
          struct OBJECT_NAME<Initial, Operand1, Operand2, FieldType> {                             \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<Initial, Operand1, Operand2, FieldType> typedef MyType;                   \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             template<typename IteratorType>                                                       \
              struct Iterator {                                                                    \
                                                                                                   \
                OBJECT_NAME<ResizePrep<IteratorType, MyType>,                                      \
                            typename Operand1::template Iterator<IteratorType>::ResizePrepType,    \
                            typename Operand2::template Iterator<IteratorType>::ResizePrepType,    \
                            FieldType> typedef ResizePrepType;                                     \
                                                                                                   \
                OBJECT_NAME<SeqWalk<IteratorType, MyType>,                                         \
                            typename Operand1::template Iterator<IteratorType>::SeqWalkType,       \
                            typename Operand2::template Iterator<IteratorType>::SeqWalkType,       \
                            FieldType> typedef SeqWalkType;                                        \
             };                                                                                    \
             OBJECT_NAME(Operand1 const & op1, Operand2 const & op2)                               \
             : operand1_(op1), operand2_(op2)                                                      \
             {};                                                                                   \
             template<typename IteratorType>                                                       \
              inline typename Iterator<IteratorType>::SeqWalkType init (void) const {              \
                 return typename Iterator<IteratorType>::SeqWalkType(*this);                       \
              };                                                                                   \
             template<typename IteratorType>                                                       \
              inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {    \
                 return typename Iterator<IteratorType>::ResizePrepType (*this);                   \
              };                                                                                   \
             inline Operand1 const & operand1 (void) const { return operand1_; };                  \
             inline Operand2 const & operand2 (void) const { return operand2_; };                  \
                                                                                                   \
            private:                                                                               \
             Operand1 const operand1_;                                                             \
             Operand2 const operand2_;                                                             \
         };                                                                                        \
                                                                                                   \
         template<typename IteratorType,                                                           \
                  typename SourceType,                                                             \
                  typename Operand1,                                                               \
                  typename Operand2,                                                               \
                  typename FieldType>                                                              \
          struct OBJECT_NAME<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> { \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType>      \
             typedef MyType;                                                                       \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             OBJECT_NAME<Resize<IteratorType, MyType>,                                             \
                         typename Operand1::ResizeType,                                            \
                         typename Operand2::ResizeType,                                            \
                         FieldType> typedef ResizeType;                                            \
             OBJECT_NAME(SourceType const & source)                                                \
             : operand1_(source.operand1()), operand2_(source.operand2())                          \
             {};                                                                                   \
             inline ResizeType resize ( MemoryWindow const & newSize) const {                      \
                return ResizeType(newSize, *this);                                                 \
             };                                                                                    \
             inline Operand1 const & operand1 (void) const { return operand1_; };                  \
             inline Operand2 const & operand2 (void) const { return operand2_; };                  \
                                                                                                   \
            private:                                                                               \
             Operand1 const operand1_;                                                             \
             Operand2 const operand2_;                                                             \
         };                                                                                        \
                                                                                                   \
         template<typename IteratorType,                                                           \
                  typename SourceType,                                                             \
                  typename Operand1,                                                               \
                  typename Operand2,                                                               \
                  typename FieldType>                                                              \
          struct OBJECT_NAME<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> {    \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef  \
             MyType;                                                                               \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             OBJECT_NAME<SeqWalk<IteratorType, MyType>,                                            \
                         typename Operand1::SeqWalkType,                                           \
                         typename Operand2::SeqWalkType,                                           \
                         FieldType> typedef SeqWalkType;                                           \
             OBJECT_NAME(MemoryWindow const & size, SourceType const & source)                     \
             : operand1_(size, source.operand1()), operand2_(size, source.operand2())              \
             {};                                                                                   \
             inline SeqWalkType init (void) const { return SeqWalkType(*this); };                  \
             inline Operand1 const & operand1 (void) const { return operand1_; };                  \
             inline Operand2 const & operand2 (void) const { return operand2_; };                  \
                                                                                                   \
            private:                                                                               \
             Operand1 const operand1_;                                                             \
             Operand2 const operand2_;                                                             \
         };                                                                                        \
                                                                                                   \
         template<typename IteratorType,                                                           \
                  typename SourceType,                                                             \
                  typename Operand1,                                                               \
                  typename Operand2,                                                               \
                  typename FieldType>                                                              \
          struct OBJECT_NAME<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> {   \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef \
             MyType;                                                                               \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             OBJECT_NAME(SourceType const & source)                                                \
             : operand1_(source.operand1()), operand2_(source.operand2())                          \
             {};                                                                                   \
             inline void next (void) { operand1_.next(); operand2_.next(); };                      \
             inline bool at_end (void) const { return (operand1_.at_end() || operand2_.at_end()); }; \
             inline bool has_length (void) const {                                                 \
                return (operand1_.has_length() || operand2_.has_length());                         \
             };                                                                                    \
             inline AtomicType eval (void) const {                                                 \
                return INTERNAL_NAME(operand1_.eval(), operand2_.eval());                          \
             };                                                                                    \
                                                                                                   \
            private:                                                                               \
             Operand1 operand1_;                                                                   \
             Operand2 operand2_;                                                                   \
         };                                                                                        \
                                                                                                   \
         /* SubExpr X SubExpr */                                                                   \
         template<typename SubExpr1, typename SubExpr2>                                            \
          inline FieldExpression<OBJECT_NAME<Initial,                                              \
                                             typename Standardize<SubExpr1,                        \
                                                                  typename SubExpr1::field_type>:: \
                                             StandardType,                                         \
                                             typename Standardize<SubExpr2,                        \
                                                                  typename SubExpr1::field_type>:: \
                                             StandardType,                                         \
                                             typename SubExpr1::field_type>,                       \
                                 typename SubExpr1::field_type> EXTERNAL_NAME (SubExpr1 const & arg1, \
                                                                               SubExpr2 const & arg2) { \
                                                                                                   \
             typename SubExpr1::field_type typedef FieldType;                                      \
                                                                                                   \
             typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef   \
             Type1;                                                                                \
                                                                                                   \
             typename Standardize<SubExpr2, typename SubExpr1::field_type>::StandardType typedef   \
             Type2;                                                                                \
                                                                                                   \
             OBJECT_NAME<Initial, Type1, Type2, FieldType> typedef ReturnType;                     \
                                                                                                   \
             FieldExpression<ReturnType, FieldType> typedef ReturnTerm;                            \
                                                                                                   \
             return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)), \
                                          Type2(Standardize<SubExpr2, FieldType>::standardType(arg2)))); \
          };                                                                                       \
                                                                                                   \
         /* SubExpr X Scalar */                                                                    \
         template<typename SubExpr1>                                                               \
          inline FieldExpression<OBJECT_NAME<Initial,                                              \
                                             typename Standardize<SubExpr1,                        \
                                                                  typename SubExpr1::field_type>:: \
                                             StandardType,                                         \
                                             Scalar<Initial, typename SubExpr1::field_type>,       \
                                             typename SubExpr1::field_type>,                       \
                                 typename SubExpr1::field_type> EXTERNAL_NAME (SubExpr1 const & arg1, \
                                                                               typename SubExpr1:: \
                                                                               field_type::        \
                                                                               value_type const &  \
                                                                               arg2) {             \
                                                                                                   \
             typename SubExpr1::field_type typedef FieldType;                                      \
                                                                                                   \
             typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef   \
             Type1;                                                                                \
                                                                                                   \
             Scalar<Initial, typename SubExpr1::field_type> typedef Type2;                         \
                                                                                                   \
             OBJECT_NAME<Initial, Type1, Type2, FieldType> typedef ReturnType;                     \
                                                                                                   \
             FieldExpression<ReturnType, FieldType> typedef ReturnTerm;                            \
                                                                                                   \
             return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)), \
                                          Type2(Type2(arg2))));                                    \
          };                                                                                       \
                                                                                                   \
         /* Scalar X SubExpr */                                                                    \
         template<typename SubExpr2>                                                               \
          inline FieldExpression<OBJECT_NAME<Initial,                                              \
                                             Scalar<Initial, typename SubExpr2::field_type>,       \
                                             typename Standardize<SubExpr2,                        \
                                                                  typename SubExpr2::field_type>:: \
                                             StandardType,                                         \
                                             typename SubExpr2::field_type>,                       \
                                 typename SubExpr2::field_type> EXTERNAL_NAME (typename SubExpr2:: \
                                                                               field_type::        \
                                                                               value_type const &  \
                                                                               arg1,               \
                                                                               SubExpr2 const & arg2) { \
                                                                                                   \
             typename SubExpr2::field_type typedef FieldType;                                      \
                                                                                                   \
             Scalar<Initial, typename SubExpr2::field_type> typedef Type1;                         \
                                                                                                   \
             typename Standardize<SubExpr2, typename SubExpr2::field_type>::StandardType typedef   \
             Type2;                                                                                \
                                                                                                   \
             OBJECT_NAME<Initial, Type1, Type2, FieldType> typedef ReturnType;                     \
                                                                                                   \
             FieldExpression<ReturnType, FieldType> typedef ReturnTerm;                            \
                                                                                                   \
             return ReturnTerm(ReturnType(Type1(Type1(arg1)),                                      \
                                          Type2(Standardize<SubExpr2, FieldType>::standardType(arg2)))); \
          };

#     define BUILD_BINARY_OPERATOR(OBJECT_NAME, INTERNAL_NAME, EXTERNAL_NAME)                      \
         template<typename CurrentMode, typename Operand1, typename Operand2, typename FieldType>  \
          struct OBJECT_NAME;                                                                      \
                                                                                                   \
         template<typename Operand1, typename Operand2, typename FieldType>                        \
          struct OBJECT_NAME<Initial, Operand1, Operand2, FieldType> {                             \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<Initial, Operand1, Operand2, FieldType> typedef MyType;                   \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             template<typename IteratorType>                                                       \
              struct Iterator {                                                                    \
                                                                                                   \
                OBJECT_NAME<ResizePrep<IteratorType, MyType>,                                      \
                            typename Operand1::template Iterator<IteratorType>::ResizePrepType,    \
                            typename Operand2::template Iterator<IteratorType>::ResizePrepType,    \
                            FieldType> typedef ResizePrepType;                                     \
                                                                                                   \
                OBJECT_NAME<SeqWalk<IteratorType, MyType>,                                         \
                            typename Operand1::template Iterator<IteratorType>::SeqWalkType,       \
                            typename Operand2::template Iterator<IteratorType>::SeqWalkType,       \
                            FieldType> typedef SeqWalkType;                                        \
             };                                                                                    \
             OBJECT_NAME(Operand1 const & op1, Operand2 const & op2)                               \
             : operand1_(op1), operand2_(op2)                                                      \
             {};                                                                                   \
             template<typename IteratorType>                                                       \
              inline typename Iterator<IteratorType>::SeqWalkType init (void) const {              \
                 return typename Iterator<IteratorType>::SeqWalkType(*this);                       \
              };                                                                                   \
             template<typename IteratorType>                                                       \
              inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {    \
                 return typename Iterator<IteratorType>::ResizePrepType (*this);                   \
              };                                                                                   \
             inline Operand1 const & operand1 (void) const { return operand1_; };                  \
             inline Operand2 const & operand2 (void) const { return operand2_; };                  \
                                                                                                   \
            private:                                                                               \
             Operand1 const operand1_;                                                             \
             Operand2 const operand2_;                                                             \
         };                                                                                        \
                                                                                                   \
         template<typename IteratorType,                                                           \
                  typename SourceType,                                                             \
                  typename Operand1,                                                               \
                  typename Operand2,                                                               \
                  typename FieldType>                                                              \
          struct OBJECT_NAME<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType> { \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<ResizePrep<IteratorType, SourceType>, Operand1, Operand2, FieldType>      \
             typedef MyType;                                                                       \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             OBJECT_NAME<Resize<IteratorType, MyType>,                                             \
                         typename Operand1::ResizeType,                                            \
                         typename Operand2::ResizeType,                                            \
                         FieldType> typedef ResizeType;                                            \
             OBJECT_NAME(SourceType const & source)                                                \
             : operand1_(source.operand1()), operand2_(source.operand2())                          \
             {};                                                                                   \
             inline ResizeType resize ( MemoryWindow const & newSize) const {                      \
                return ResizeType(newSize, *this);                                                 \
             };                                                                                    \
             inline Operand1 const & operand1 (void) const { return operand1_; };                  \
             inline Operand2 const & operand2 (void) const { return operand2_; };                  \
                                                                                                   \
            private:                                                                               \
             Operand1 const operand1_;                                                             \
             Operand2 const operand2_;                                                             \
         };                                                                                        \
                                                                                                   \
         template<typename IteratorType,                                                           \
                  typename SourceType,                                                             \
                  typename Operand1,                                                               \
                  typename Operand2,                                                               \
                  typename FieldType>                                                              \
          struct OBJECT_NAME<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> {    \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<Resize<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef  \
             MyType;                                                                               \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             OBJECT_NAME<SeqWalk<IteratorType, MyType>,                                            \
                         typename Operand1::SeqWalkType,                                           \
                         typename Operand2::SeqWalkType,                                           \
                         FieldType> typedef SeqWalkType;                                           \
             OBJECT_NAME(MemoryWindow const & size, SourceType const & source)                     \
             : operand1_(size, source.operand1()), operand2_(size, source.operand2())              \
             {};                                                                                   \
             inline SeqWalkType init (void) const { return SeqWalkType(*this); };                  \
             inline Operand1 const & operand1 (void) const { return operand1_; };                  \
             inline Operand2 const & operand2 (void) const { return operand2_; };                  \
                                                                                                   \
            private:                                                                               \
             Operand1 const operand1_;                                                             \
             Operand2 const operand2_;                                                             \
         };                                                                                        \
                                                                                                   \
         template<typename IteratorType,                                                           \
                  typename SourceType,                                                             \
                  typename Operand1,                                                               \
                  typename Operand2,                                                               \
                  typename FieldType>                                                              \
          struct OBJECT_NAME<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> {   \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<SeqWalk<IteratorType, SourceType>, Operand1, Operand2, FieldType> typedef \
             MyType;                                                                               \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             OBJECT_NAME(SourceType const & source)                                                \
             : operand1_(source.operand1()), operand2_(source.operand2())                          \
             {};                                                                                   \
             inline void next (void) { operand1_.next(); operand2_.next(); };                      \
             inline bool at_end (void) const { return (operand1_.at_end() || operand2_.at_end()); }; \
             inline bool has_length (void) const {                                                 \
                return (operand1_.has_length() || operand2_.has_length());                         \
             };                                                                                    \
             inline AtomicType eval (void) const {                                                 \
                return (operand1_.eval() INTERNAL_NAME operand2_.eval());                          \
             };                                                                                    \
                                                                                                   \
            private:                                                                               \
             Operand1 operand1_;                                                                   \
             Operand2 operand2_;                                                                   \
         };                                                                                        \
                                                                                                   \
         /* SubExpr X SubExpr */                                                                   \
         template<typename SubExpr1, typename SubExpr2>                                            \
          inline FieldExpression<OBJECT_NAME<Initial,                                              \
                                             typename Standardize<SubExpr1,                        \
                                                                  typename SubExpr1::field_type>:: \
                                             StandardType,                                         \
                                             typename Standardize<SubExpr2,                        \
                                                                  typename SubExpr1::field_type>:: \
                                             StandardType,                                         \
                                             typename SubExpr1::field_type>,                       \
                                 typename SubExpr1::field_type> EXTERNAL_NAME (SubExpr1 const & arg1, \
                                                                               SubExpr2 const & arg2) { \
                                                                                                   \
             typename SubExpr1::field_type typedef FieldType;                                      \
                                                                                                   \
             typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef   \
             Type1;                                                                                \
                                                                                                   \
             typename Standardize<SubExpr2, typename SubExpr1::field_type>::StandardType typedef   \
             Type2;                                                                                \
                                                                                                   \
             OBJECT_NAME<Initial, Type1, Type2, FieldType> typedef ReturnType;                     \
                                                                                                   \
             FieldExpression<ReturnType, FieldType> typedef ReturnTerm;                            \
                                                                                                   \
             return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)), \
                                          Type2(Standardize<SubExpr2, FieldType>::standardType(arg2)))); \
          };                                                                                       \
                                                                                                   \
         /* SubExpr X Scalar */                                                                    \
         template<typename SubExpr1>                                                               \
          inline FieldExpression<OBJECT_NAME<Initial,                                              \
                                             typename Standardize<SubExpr1,                        \
                                                                  typename SubExpr1::field_type>:: \
                                             StandardType,                                         \
                                             Scalar<Initial, typename SubExpr1::field_type>,       \
                                             typename SubExpr1::field_type>,                       \
                                 typename SubExpr1::field_type> EXTERNAL_NAME (SubExpr1 const & arg1, \
                                                                               typename SubExpr1:: \
                                                                               field_type::        \
                                                                               value_type const &  \
                                                                               arg2) {             \
                                                                                                   \
             typename SubExpr1::field_type typedef FieldType;                                      \
                                                                                                   \
             typename Standardize<SubExpr1, typename SubExpr1::field_type>::StandardType typedef   \
             Type1;                                                                                \
                                                                                                   \
             Scalar<Initial, typename SubExpr1::field_type> typedef Type2;                         \
                                                                                                   \
             OBJECT_NAME<Initial, Type1, Type2, FieldType> typedef ReturnType;                     \
                                                                                                   \
             FieldExpression<ReturnType, FieldType> typedef ReturnTerm;                            \
                                                                                                   \
             return ReturnTerm(ReturnType(Type1(Standardize<SubExpr1, FieldType>::standardType(arg1)), \
                                          Type2(Type2(arg2))));                                    \
          };                                                                                       \
                                                                                                   \
         /* Scalar X SubExpr */                                                                    \
         template<typename SubExpr2>                                                               \
          inline FieldExpression<OBJECT_NAME<Initial,                                              \
                                             Scalar<Initial, typename SubExpr2::field_type>,       \
                                             typename Standardize<SubExpr2,                        \
                                                                  typename SubExpr2::field_type>:: \
                                             StandardType,                                         \
                                             typename SubExpr2::field_type>,                       \
                                 typename SubExpr2::field_type> EXTERNAL_NAME (typename SubExpr2:: \
                                                                               field_type::        \
                                                                               value_type const &  \
                                                                               arg1,               \
                                                                               SubExpr2 const & arg2) { \
                                                                                                   \
             typename SubExpr2::field_type typedef FieldType;                                      \
                                                                                                   \
             Scalar<Initial, typename SubExpr2::field_type> typedef Type1;                         \
                                                                                                   \
             typename Standardize<SubExpr2, typename SubExpr2::field_type>::StandardType typedef   \
             Type2;                                                                                \
                                                                                                   \
             OBJECT_NAME<Initial, Type1, Type2, FieldType> typedef ReturnType;                     \
                                                                                                   \
             FieldExpression<ReturnType, FieldType> typedef ReturnTerm;                            \
                                                                                                   \
             return ReturnTerm(ReturnType(Type1(Type1(arg1)),                                      \
                                          Type2(Standardize<SubExpr2, FieldType>::standardType(arg2)))); \
          };

#     define BUILD_UNARY_FUNCTION(OBJECT_NAME, INTERNAL_NAME, EXTERNAL_NAME)                       \
         template<typename CurrentMode, typename Operand, typename FieldType>                      \
          struct OBJECT_NAME;                                                                      \
                                                                                                   \
         template<typename Operand, typename FieldType>                                            \
          struct OBJECT_NAME<Initial, Operand, FieldType> {                                        \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<Initial, Operand, FieldType> typedef MyType;                              \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             template<typename IteratorType>                                                       \
              struct Iterator {                                                                    \
                                                                                                   \
                OBJECT_NAME<ResizePrep<IteratorType, MyType>,                                      \
                            typename Operand::template Iterator<IteratorType>::ResizePrepType,     \
                            FieldType> typedef ResizePrepType;                                     \
                                                                                                   \
                OBJECT_NAME<SeqWalk<IteratorType, MyType>,                                         \
                            typename Operand::template Iterator<IteratorType>::SeqWalkType,        \
                            FieldType> typedef SeqWalkType;                                        \
             };                                                                                    \
             OBJECT_NAME(Operand const & op)                                                       \
             : operand_(op)                                                                        \
             {};                                                                                   \
             template<typename IteratorType>                                                       \
              inline typename Iterator<IteratorType>::SeqWalkType init (void) const {              \
                 return typename Iterator<IteratorType>::SeqWalkType(*this);                       \
              };                                                                                   \
             template<typename IteratorType>                                                       \
              inline typename Iterator<IteratorType>::ResizePrepType resize_prep (void) const {    \
                 return typename Iterator<IteratorType>::ResizePrepType (*this);                   \
              };                                                                                   \
             inline Operand const & operand (void) const { return operand_; };                     \
                                                                                                   \
            private:                                                                               \
             Operand const operand_;                                                               \
         };                                                                                        \
                                                                                                   \
         template<typename IteratorType, typename SourceType, typename Operand, typename FieldType> \
          struct OBJECT_NAME<ResizePrep<IteratorType, SourceType>, Operand, FieldType> {           \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<ResizePrep<IteratorType, SourceType>, Operand, FieldType> typedef MyType; \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             OBJECT_NAME<Resize<IteratorType, MyType>, typename Operand::ResizeType, FieldType>    \
             typedef ResizeType;                                                                   \
             OBJECT_NAME(SourceType const & source)                                                \
             : operand_(source.operand())                                                          \
             {};                                                                                   \
             inline ResizeType resize ( MemoryWindow const & newSize) const {                      \
                return ResizeType(newSize, *this);                                                 \
             };                                                                                    \
             inline Operand const & operand (void) const { return operand_; };                     \
                                                                                                   \
            private:                                                                               \
             Operand const operand_;                                                               \
         };                                                                                        \
                                                                                                   \
         template<typename IteratorType, typename SourceType, typename Operand, typename FieldType> \
          struct OBJECT_NAME<Resize<IteratorType, SourceType>, Operand, FieldType> {               \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<Resize<IteratorType, SourceType>, Operand, FieldType> typedef MyType;     \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             OBJECT_NAME<SeqWalk<IteratorType, MyType>, typename Operand::SeqWalkType, FieldType>  \
             typedef SeqWalkType;                                                                  \
             OBJECT_NAME(MemoryWindow const & size, SourceType const & source)                     \
             : operand_(size, source.operand())                                                    \
             {};                                                                                   \
             inline SeqWalkType init (void) const { return SeqWalkType(*this); };                  \
             inline Operand const & operand (void) const { return operand_; };                     \
                                                                                                   \
            private:                                                                               \
             Operand const operand_;                                                               \
         };                                                                                        \
                                                                                                   \
         template<typename IteratorType, typename SourceType, typename Operand, typename FieldType> \
          struct OBJECT_NAME<SeqWalk<IteratorType, SourceType>, Operand, FieldType> {              \
                                                                                                   \
            public:                                                                                \
             typename FieldType::value_type typedef AtomicType;                                    \
             OBJECT_NAME<SeqWalk<IteratorType, SourceType>, Operand, FieldType> typedef MyType;    \
             typename FieldType::memory_window typedef MemoryWindow;                               \
             OBJECT_NAME(SourceType const & source)                                                \
             : operand_(source.operand())                                                          \
             {};                                                                                   \
             inline void next (void) { operand_.next(); };                                         \
             inline bool at_end (void) const { return (operand_.at_end()); };                      \
             inline bool has_length (void) const { return (operand_.has_length()); };              \
             inline AtomicType eval (void) const { return INTERNAL_NAME(operand_.eval()); };       \
                                                                                                   \
            private:                                                                               \
             Operand operand_;                                                                     \
         };                                                                                        \
                                                                                                   \
         /* SubExpr */                                                                             \
         template<typename SubExpr>                                                                \
          inline FieldExpression<OBJECT_NAME<Initial,                                              \
                                             typename Standardize<SubExpr,                         \
                                                                  typename SubExpr::field_type>::  \
                                             StandardType,                                         \
                                             typename SubExpr::field_type>,                        \
                                 typename SubExpr::field_type> EXTERNAL_NAME (SubExpr const & arg) { \
                                                                                                   \
             typename SubExpr::field_type typedef FieldType;                                       \
                                                                                                   \
             typename Standardize<SubExpr, typename SubExpr::field_type>::StandardType typedef Type; \
                                                                                                   \
             OBJECT_NAME<Initial, Type, FieldType> typedef ReturnType;                             \
                                                                                                   \
             FieldExpression<ReturnType, FieldType> typedef ReturnTerm;                            \
                                                                                                   \
             return ReturnTerm(ReturnType(Type(Standardize<SubExpr, FieldType>::standardType(arg)))); \
          };

      template<typename LhsType, typename RhsType>
       inline void field_expression_sequential_execute_internal (LhsType lhs, RhsType rhs) {
          while(!lhs.at_end()){ lhs.ref() = rhs.eval(); lhs.next(); rhs.next(); };
       };

      template<typename CallStyle, typename ExprType, typename FieldType>
       inline FieldType const & field_expression_sequential_execute (FieldType & initial_lhs,
                                                                     FieldExpression<ExprType,
                                                                                     FieldType>
                                                                     const & initial_rhs) {

          field_expression_sequential_execute_internal<typename Field<Initial, FieldType>::template
                                                       Iterator<CallStyle>::SeqWalkType,
                                                       typename ExprType::template Iterator<CallStyle>
                                                       ::SeqWalkType>(Field<Initial, FieldType>(initial_lhs).template
                                                                                                             init<CallStyle>(),
                                                                      initial_rhs.expr().template
                                                                                         init<CallStyle>());

          return initial_lhs;
       };

#     ifdef FIELD_EXPRESSION_THREADS
         template<typename CallStyle,
                  typename ResizeLhsType,
                  typename ResizeRhsType,
                  typename FieldType>
          inline void field_expression_thread_parallel_execute_internal (ResizeLhsType & lhs,
                                                                         ResizeRhsType const & rhs,
                                                                         typename FieldType::
                                                                         memory_window const &
                                                                         window,
                                                                         BI::interprocess_semaphore
                                                                         * semaphore) {

             field_expression_sequential_execute_internal<typename ResizeLhsType::ResizeType::
                                                          SeqWalkType,
                                                          typename ResizeRhsType::ResizeType::
                                                          SeqWalkType>(lhs.resize(window).init(),
                                                                       rhs.resize(window).init());

             semaphore->post();
          }
#     endif //FIELD_EXPRESSION_THREADS;

#     ifdef FIELD_EXPRESSION_THREADS
         template<typename CallStyle, typename ExprType, typename FieldType>
          inline FieldType const & field_expression_thread_parallel_execute (FieldType & initial_lhs,
                                                                             FieldExpression<ExprType,
                                                                                             FieldType>
                                                                             const & initial_rhs,
                                                                             int const
                                                                             number_of_partitions) {

             typename Field<Initial, FieldType>::template Iterator<CallStyle>::ResizePrepType
             typedef LhsType;

             typename ExprType::template Iterator<CallStyle>::ResizePrepType typedef RhsType;

             typename FieldType::memory_window typedef MemoryWindow;

             MemoryWindow window = IteratorStyle<CallStyle, FieldType>::memory_window(initial_lhs);

             int x = 1;
             int y = 1;
             int z = 1;

             if(number_of_partitions <= window.extent(2)){ z = number_of_partitions; }
             else if(number_of_partitions <= window.extent(1)){ y = number_of_partitions; }
             else if(number_of_partitions <= window.extent(0)){ x = number_of_partitions; };

             std::vector<typename FieldType::memory_window> vec_window = window.split(structured::
                                                                                      IntVec(x, y, z));

             std::vector<BI::interprocess_semaphore *> vec_semaphore;

             typename std::vector<typename FieldType::memory_window>::const_iterator window_iterator
             = vec_window.begin();

             typename std::vector<typename FieldType::memory_window>::const_iterator window_end =
             vec_window.end();

             for(; window_iterator != window_end; ++window_iterator){

                vec_semaphore.push_back(new BI::interprocess_semaphore(0));

                ThreadPoolFIFO::self().schedule(boost::bind(&
                                                            field_expression_thread_parallel_execute_internal<CallStyle,
                                                                                                              LhsType,
                                                                                                              RhsType,
                                                                                                              FieldType>,
                                                            Field<Initial, FieldType>(initial_lhs).template
                                                                                                   resize_prep<CallStyle>(),
                                                            initial_rhs.expr().template resize_prep<CallStyle>(),
                                                            *window_iterator,
                                                            vec_semaphore.back()));
             };

             std::vector<BI::interprocess_semaphore *>::iterator isem = vec_semaphore.begin();

             std::vector<BI::interprocess_semaphore *>::iterator const esem = vec_semaphore.end();

             for(; isem != esem; ++isem){ (*isem)->wait(); };

             for(isem = vec_semaphore.begin(); isem != esem; ++isem){ delete *isem; };

             return initial_lhs;
          }
#     endif //FIELD_EXPRESSION_THREADS;

      template<typename CallStyle, typename ExprType, typename FieldType>
       inline FieldType const & field_expression_general_execute (FieldType & initial_lhs,
                                                                  FieldExpression<ExprType,
                                                                                  FieldType> const &
                                                                  initial_rhs) {

          return
#                ifdef FIELD_EXPRESSION_THREADS
                    field_expression_thread_parallel_execute<CallStyle, ExprType, FieldType>(initial_lhs,
                                                                                             initial_rhs,
                                                                                             NTHREADS)
#                else
                    field_expression_sequential_execute<CallStyle, ExprType, FieldType>(initial_lhs,
                                                                                        initial_rhs)
#                endif //FIELD_EXPRESSION_THREADS
                 ;;
       };

      template<typename FieldType>
       inline FieldType const & operator <<= (FieldType & lhs,
                                              typename FieldType::value_type const & rhs) {

          Scalar<Initial, FieldType> typedef ExprType;

          return (lhs <<= FieldExpression<ExprType, FieldType>(ExprType(rhs)));
       };

      template<typename FieldType>
       inline FieldType const & operator <<= (FieldType & lhs, FieldType const & rhs) {

          ConstField<Initial, FieldType> typedef ExprType;

          return (lhs <<= FieldExpression<ExprType, FieldType>(ExprType(rhs)));
       };

      template<typename ExprType, typename FieldType>
       inline FieldType const & operator <<= (FieldType & lhs,
                                              FieldExpression<ExprType, FieldType> const & rhs) {
          return field_expression_general_execute<UseWholeIterator, ExprType, FieldType>(lhs, rhs);
       };

      template<typename FieldType>
       inline FieldType const & interior_assign (FieldType & lhs,
                                                 typename FieldType::value_type const & rhs) {

          Scalar<Initial, FieldType> typedef ExprType;

          return interior_assign(lhs, FieldExpression<ExprType, FieldType>(ExprType(rhs)));
       };

      template<typename FieldType>
       inline FieldType const & interior_assign (FieldType & lhs, FieldType const & rhs) {

          ConstField<Initial, FieldType> typedef ExprType;

          return interior_assign(lhs, FieldExpression<ExprType, FieldType>(ExprType(rhs)));
       };

      template<typename ExprType, typename FieldType>
       inline FieldType const & interior_assign (FieldType & lhs,
                                                 FieldExpression<ExprType, FieldType> const & rhs) {
          return field_expression_general_execute<UseInteriorIterator, ExprType, FieldType>(lhs, rhs);
       };
   } //SpatialOps;

#endif //SpatialOps_FieldExpressions_h
