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

#ifndef SpatialOps_FieldTill_h
#define SpatialOps_FieldTill_h

#include <vector>
//#include <iostream>

#include <spatialops/SpatialOpsConfigure.h>

#include <spatialops/FieldExpressionsMacros.h>
#include <spatialops/CoreMacrosDefine.h>
#include <spatialops/FieldExpressions.h>

namespace SpatialOps{
  //namespace structured{
  /**
   *  @file FieldExpressions.h
   */

  struct CallByVal;
  struct CallByRef;
  struct CallByPtr;

  //currently, TillMutableField is FieldForm (provides both ptr() and ref() interfaces)

  template<typename Mode, typename FieldType>
    struct TillVec_ConstField;

  BUILD_IMMUT_FIELD_EXPR_NO_SUBEXPR(/*Structure's Name*/
				    TillVec_ConstField
				    ,/*Template Parameters WITH COMMA*/
				    NOTHING
				    ,/*Template Arguments WITH COMMA*/
				    NOTHING
				    ,/*Initial's Typedefs*/
				    std::vector<ConstFieldType const *>              typedef vec_InputConstFieldType;
				    std::vector<   NCFieldType const *>              typedef vec_InputNCFieldType;
				    std::vector<ConstFieldType        >              typedef vec_ConstFieldType;
				    template<typename IteratorType>
				    struct Iterator2 {
				      ConstFieldInfo<IteratorType QC ConstFieldType> typedef     ConstFieldInfoType;
				      std::vector<ConstFieldInfoType>                typedef vec_ConstFieldInfoType;
				    }
				    ,/*Initial's Data Members*/
				    vec_ConstFieldType fields_;
				    ,/*Initial's Private Functions*/
				    NOTHING
				    ,/*Initial's Constructor Parameters*/
				    vec_InputConstFieldType const & fields
				    ,/*Initial's Constructor Assignments*/
				    fields_()
				    ,/*Initial's Constructor Computation*/
				    typename vec_InputConstFieldType::const_iterator i = fields.begin();
				    typename vec_InputConstFieldType::const_iterator e = fields.end  ();
				    for(; i != e; ++i)
				      fields_.push_back(*(*i));
				    ,/*Initial's External Functions*/
				    TillVec_ConstField(vec_InputNCFieldType const & fs) {
				      for(typename vec_InputNCFieldType::const_iterator i = fs.begin();
					  i != fs.end();
					  ++i)
					fields_.push_back(ConstFieldType(*(*i)));
				    };
				    template<typename IteratorType>
				    I typename Iterator2<IteratorType>::vec_ConstFieldInfoType const info() const {
				      typename Iterator2<IteratorType>::vec_ConstFieldInfoType vec;
				      for(typename vec_ConstFieldType::const_iterator i = fields_.begin();
					  i != fields_.end();
					  ++i)
					vec.push_back(typename Iterator2<IteratorType>::ConstFieldInfoType(*i));
				      return vec;
				    };
				    I vec_ConstFieldType const & fields() const { return fields_; }
				    ,/*ResizePrep's Typedefs*/
				    ConstFieldInfo<IteratorType QC ConstFieldType> typedef     ConstFieldInfoType;
				    std::vector<ConstFieldInfoType>                typedef vec_ConstFieldInfoType
				    ,/*ResizePrep's Data Members*/
				    vec_ConstFieldInfoType current_info;
				    ,/*ResizePrep's Private Functions*/
				    NOTHING
				    ,/*ResizePrep's Constructor Assignments*/
				    current_info(source.template info<IteratorType>())
				    ,/*ResizePrep's Constructor Computation*/
				    NOTHING
				    ,/*ResizePrep's public functions*/
				    template<typename Type>
				    I vec_ConstFieldInfoType const & info() const { return current_info; }
				    ,/*Resize's Typedefs*/
				    ConstFieldInfo<IteratorType QC ConstFieldType> typedef     ConstFieldInfoType;
				    std::vector<ConstFieldInfoType>                typedef vec_ConstFieldInfoType;
				    std::vector<ConstFieldType>                    typedef vec_ConstFieldType
				    ,/*Resize's Data Members*/
				    vec_ConstFieldType fields_;
				    ,/*Resize's Private Functions*/
				    NOTHING
				    ,/*Resize's Constructor Assignments*/
				    fields_()
				    ,/*Resize's Constructor Computation*/
				    typename vec_ConstFieldInfoType::const_iterator       i = source.template info<IteratorType>().begin();
				    typename vec_ConstFieldInfoType::const_iterator const e = source.template info<IteratorType>().end  ();
				    for(;
					i != e;
					++i) {
				      fields_.push_back(i->resize(size).new_field());
				      /* current_info.push_back(i->resize(size)); */
				    }
				    ,/*Resize's public functions*/
				    I vec_ConstFieldType const & fields() const { return fields_; }
				    ,/*SeqWalk's Typedefs*/
				    std::vector<ConstFieldType>               typedef vec_ConstFieldType;
				    IterFcns<ConstFieldType>                  typedef     Iter;
				    typename Iter::const_iterator_type        typedef     IterType;
				    std::vector<IterType>                     typedef vec_IterType;
				    std::vector<AtomicType>                   typedef vec_AtomicType;
				    ,/*SeqWalk's Data Members*/
				    vec_IterType   iters;
				    vec_IterType   ends;
				    vec_AtomicType values;
				    ,/*SeqWalk's Private Functions*/
				    I void fill() {
				      values.clear();
				      for(typename vec_IterType::const_iterator i = iters.begin();
					  i != iters.end();
					  ++i)
					values.push_back(*(*i));
				    };
				    I void increment() {
				      for(typename vec_IterType::iterator i = iters.begin();
					  i != iters.end();
					  ++i)
					++(*i);
				    };
				    ,/*SeqWalk's Constructor Assignments*/
				    iters () QC
				    ends  () QC
				    values()
				    ,/*SeqWalk's Constructor Computation*/
				    typename vec_ConstFieldType::const_iterator ifs = source.fields().begin();
				    typename vec_ConstFieldType::const_iterator efs = source.fields().end  ();
				    for(;
					ifs != efs;
					++ifs) {
				      iters.push_back(ifs->begin());
				      ends .push_back(ifs->end  ());
				    };
				    fill()
				    ,/*SeqWalk's public functions*/
				    NOTHING
				    ,/*SeqWalk's next Function Computation*/
				    increment();
				    fill()
				    ,/*SeqWalk's Accessors*/
				    NOTHING
				    ,/*SeqWalk's eval Function Return Type*/
				    vec_AtomicType
				    ,/*SeqWalk's eval Function Return Type Qualifiers*/
				    const &
				    ,/*SeqWalk's eval Function Expression*/
				    values
				    ,/*SeqWalk's at_end Function Expression*/
				    has_length() ? false : (iters.front() == ends.front())
				    ,/*SeqWalk's has_length Function Expression*/
				    !(ends.empty())
				    );

  template<typename Mode, typename FieldType>
    struct TillVec_Field;

  BUILD_MUT_FIELD_EXPR(/*Structure's Name*/
		       TillVec_Field
		       ,/*Template Parameters WITH COMMA*/
		       NOTHING
		       ,/*Template Arguments WITH COMMA*/
		       NOTHING
		       ,/*Initial's Typedefs*/
		       std::vector<FieldType *>              typedef vec_Field;
		       template<typename IteratorType>
		       struct Iterator2 {
			 FieldInfo<IteratorType QC FieldType> typedef FieldInfoType;
			 std::vector<FieldInfoType>           typedef vec_FieldInfoType;
		       }
		       ,/*Initial's Data Members*/
		       vec_Field * const fptrs;
		       ,/*Initial's Private Functions*/
		       NOTHING
		       ,/*Initial's Constructor Parameters*/
		       vec_Field & fields
		       ,/*Initial's Constructor Assignments*/
		       fptrs(&fields)
		       ,/*Initial's Constructor Computation*/
		       NOTHING
		       ,/*Initial's External Functions*/
		       template<typename IteratorType>
		       I typename Iterator2<IteratorType>::vec_FieldInfoType const info() const {
			 typename Iterator2<IteratorType>::vec_FieldInfoType vec;
			 for(typename vec_Field::iterator fptr = fptrs->begin();
			     fptr != fptrs->end();
			     ++fptr)
			   vec.push_back(typename Iterator2<IteratorType>::FieldInfoType(*(*fptr)));
			 return vec;
		       }
		       ,/*ResizePrep's Typedefs*/
		       FieldInfo<IteratorType QC FieldType> typedef FieldInfoType;
		       std::vector<FieldInfoType>           typedef vec_FieldInfoType
		       ,/*ResizePrep's Data Members*/
		       vec_FieldInfoType current_info;
		       ,/*ResizePrep's Private Functions*/
		       NOTHING
		       ,/*ResizePrep's Constructor Assignments*/
		       current_info(source.template info<IteratorType>())
		       ,/*ResizePrep's Constructor Computation*/
		       NOTHING
		       ,/*ResizePrep's public functions*/
		       template<typename Type>
		       I vec_FieldInfoType const & info() const { return current_info; }
		       ,/*Resize's Typedefs*/
		       FieldInfo<IteratorType QC FieldType> typedef FieldInfoType;
		       std::vector<FieldInfoType>           typedef vec_FieldInfoType
		       ,/*Resize's Data Members*/
		       vec_FieldInfoType current_info
		       ,/*Resize's Private Functions*/
		       NOTHING
		       ,/*Resize's Constructor Assignments*/
		       current_info(source.template info<IteratorType>().resize(size))
		       ,/*Resize's Constructor Computation*/
		       NOTHING
		       ,/*Resize's public functions*/
		       template<typename Type>
		       I vec_FieldInfoType & info() const { return current_info; }
		       ,/*SeqWalk's Typedefs*/
		       FieldInfo<IteratorType QC FieldType> typedef FieldInfoType;
		       std::vector<FieldInfoType>           typedef vec_FieldInfoType;
		       IterFcns<FieldType>                  typedef Iter;
		       typename Iter::iterator_type         typedef IterType;
		       std::vector<IterType>                typedef vec_IterType;
		       std::vector<AtomicType *>            typedef vec_AtomicType;
		       ,/*SeqWalk's Data Members*/
		       vec_IterType iters;
		       vec_IterType ends;
		       vec_AtomicType values;
		       ,/*SeqWalk's Private Functions*/
		       I void fill() {
			 values.clear();
			 for(typename vec_IterType::const_iterator i = iters.begin();
			     i != iters.end();
			     ++i)
			   values.push_back(&(**i));
		       };
		       I void increment() {
			 for(typename vec_IterType::iterator i = iters.begin();
			     i != iters.end();
			     ++i)
			   ++(*i);
		       };
		       ,/*SeqWalk's Constructor Assignments*/
		       iters () QC
		       ends  () QC
		       values()
		       ,/*SeqWalk's Constructor Computation*/
		       vec_FieldInfoType is = source.template info<IteratorType>();
		       vec_FieldInfoType es = source.template info<IteratorType>();
		       for(typename vec_FieldInfoType::const_iterator i = is.begin();
			   i != is.end();
			   ++i) {
			 iters.push_back(i->begin());
			 ends.push_back (i->end());
		       };
		       fill()
		       ,/*SeqWalk's public functions*/
		       NOTHING
		       ,/*SeqWalk's next Function Computation*/
		       increment();
		       fill()
		       ,/*SeqWalk's Accessors*/
		       NOTHING
		       ,/*SeqWalk's Standard Accessors' Return Type*/
		       vec_AtomicType
		       ,/*SeqWalk's ref Function Reference*/
		       values
		       ,/*SeqWalk's ptr Function Pointer*/
		       values
		       ,/*SeqWalk's at_end Function Expression*/
		       has_length() ? false : (iters.front() == ends.front())
		       ,/*SeqWalk's has_length Function Expression*/
		       !(ends.empty())
		       );

  /*sequential loop*/
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      till_2_sequential_execute_internal
		      ,/*Template Parameters*/
		      typename ProcType  QC
		      typename FirstType QC
		      typename SecondType
		      ,/*Function Parameters*/
		      ProcType   const proc   QC
		      FirstType        first  QC
		      SecondType       second
		      ,/*Function's Body*/
		      while(!first.at_end()) {            /*is there a better method?*/
			proc(first.ref(), second.eval()); /*call proc*/
			first.next();                     /*increment*/
			second.next();		          /*increment*/
		      }
		      );

  /*sequential external call*/
  BUILD_VOID_FIELD_FUNCTION(/*Function's Name*/
			    till_2_sequential_execute
			    ,/*Template Parameters WITH COMMA*/
			    typename CallStyle  QC
			    ,/*Function Parameters*/
			    void (*proc)(typename FieldType::value_type                    & QC
					 std::vector<typename FieldType::value_type> const &   )  QC
			    FieldType                            & initial_first                  QC
			    std::vector<FieldType const *> const & initial_second
			    ,/*Function's Body*/
			    till_2_sequential_execute_internal<
			    void (*)(typename FieldType::value_type                    & QC
				     std::vector<typename FieldType::value_type> const &   )                             QC
			    typename FieldForm         <Initial QC FieldType>::template Iterator<CallStyle>::SeqWalkType QC
			    typename TillVec_ConstField<Initial QC FieldType>::template Iterator<CallStyle>::SeqWalkType
			    >(proc                                                                                 QC
			      FieldForm         <Initial QC FieldType>(initial_first) .template init<CallStyle>()  QC
			      TillVec_ConstField<Initial QC FieldType>(initial_second).template init<CallStyle>()   )
			    );

#ifdef TILL_THREADS
  /*set-up for sequential loop in a single thread*/
  BUILD_VOID_FIELD_FUNCTION(/*Function's Name*/
			    till_2_thread_parallel_execute_internal
			    ,/*Template Parameters WITH COMMA*/
			    typename CallStyle        QC
			    typename ResizeFirstType  QC
			    typename ResizeSecondType QC
			    ,/*Function Parameters*/
			    void (*proc)(typename FieldType::value_type                    & QC
					 std::vector<typename FieldType::value_type> const &  )  QC
			    ResizeFirstType        & first                                QC
			    ResizeSecondType const & second                               QC
			    NewPartition     const & partition
			    ,/*Function's Body*/
			    till_2_sequential_execute_internal<
			    void (*)(typename FieldType::value_type                    & QC
				     std::vector<typename FieldType::value_type> const &   ) QC
			    typename ResizeFirstType ::ResizeType::SeqWalkType               QC
			    typename ResizeSecondType::ResizeType::SeqWalkType
			    >(proc                            QC
			      first .resize(partition).init() QC
			      second.resize(partition).init()  )
			    );
#endif

#ifdef TILL_THREADS
  /*set-up for all threads*/
  BUILD_VOID_FIELD_FUNCTION(/*Function's Name*/
			    till_2_thread_parallel_execute
			    ,/*Template Parameters WITH COMMA*/
			    typename CallStyle QC
			    ,/*Function Parameters*/
			    void (*proc)(typename FieldType::value_type                    & QC
					 std::vector<typename FieldType::value_type> const &   )  QC
			    FieldType                       & initial_first                       QC
			    std::vector<FieldType const *>  const & initial_second                QC
			    int                             const   number_of_threads
			    ,/*Function's Body*/
			    typename FieldForm         <Initial QC FieldType>::template Iterator<CallStyle>::ResizePrepType typedef FirstType;
			    typename TillVec_ConstField<Initial QC FieldType>::template Iterator<CallStyle>::ResizePrepType typedef SecondType;
			    GeneratePartitions<CallStyle QC FieldType> gen_par(initial_first QC number_of_threads);
			    for(; !gen_par.at_last(); ++gen_par) {
			      ThreadPoolFIFO::self().schedule(boost::bind(&till_2_thread_parallel_execute_internal<
								      CallStyle  QC
								      FirstType  QC
								      SecondType QC
								      FieldType>                                                                                 QC
								      proc                                                                                       QC
								      FieldForm         <Initial QC FieldType>(initial_first) .template resize_prep<CallStyle>() QC
								      TillVec_ConstField<Initial QC FieldType>(initial_second).template resize_prep<CallStyle>() QC
								      gen_par.partition()                                                                         )); };
			    ThreadPoolFIFO::self().wait()
			    );
#endif

  BUILD_VOID_FIELD_FUNCTION(/*Function's Name*/
			    till_2_general_execute
			    ,/*Template Parameters WITH COMMA*/
			    typename CallStyle QC
			    ,/*Function Parameters*/
			    void (*proc)(typename FieldType::value_type                    & QC
					 std::vector<typename FieldType::value_type> const &  ) QC
			    FieldType                            & initial_first                QC
			    std::vector<FieldType const *> const & initial_second
			    ,/*Function's Body*/
#ifdef TILL_THREADS
			    till_2_thread_parallel_execute<CallStyle QC FieldType>(proc           QC
										   initial_first  QC
										   initial_second QC
										   NTHREADS        )
#else
			    till_2_sequential_execute     <CallStyle QC FieldType>(proc           QC
										   initial_first  QC
										   initial_second  )
#endif
			    );



  template<typename FieldType>
    I void till(void (*proc)(typename FieldType::value_type &,
			     std::vector<typename FieldType::value_type> const &),
		FieldType                            & arg1                      ,
		std::vector<FieldType const *> const & arg2                      ) {

    till_2_general_execute<UseWholeIterator>(proc, arg1, arg2);
  };

  /* template<typename FieldType> */
  /* I void till(void (*proc)(typename FieldType::value_type &), */
  /* 	    FieldType & arg1) { */
  /*   typename TillMutableField<FieldType>::template FullState<UseWholeIterator,CallByRef> field = (TillMutableField<FieldType>(arg1)).template init<UseWholeIterator,CallByRef>(); */

  /*   while(!field.at_end()) { */
  /*     proc(field.ref()); */
  /*     field.next(); */
  /*   }; */
  /* }; */

  /* template<typename FieldType> */
  /* I void till(void (*proc)(typename FieldType::value_type * const), */
  /* 	    FieldType & arg1) { */
  /*   typename TillMutableField<FieldType>::template FullState<UseWholeIterator,CallByPtr> field = (TillMutableField<FieldType>(arg1)).template init<UseWholeIterator,CallByPtr>(); */

  /*   while(!field.at_end()) { */
  /*     proc(field.ptr()); */
  /*     field.next(); */
  /*   }; */
  /* }; */

  /* template<typename FieldType> */
  /* I void till(void (*proc)(typename FieldType::value_type &), */
  /* 	    FieldType * const arg1) { */
  /*   typename TillMutableField<FieldType>::template FullState<UseWholeIterator,CallByRef> field = (TillMutableField<FieldType>(arg1)).template init<UseWholeIterator,CallByRef>(); */

  /*   while(!field.at_end()) { */
  /*     proc(field.ref()); */
  /*     field.next(); */
  /*   }; */
  /* }; */

  /* template<typename FieldType> */
  /* I void till(void (*proc)(typename FieldType::value_type * const), */
  /* 	    FieldType * const arg1) { */
  /*   typename TillMutableField<FieldType>::template FullState<UseWholeIterator,CallByPtr> field = (TillMutableField<FieldType>(arg1)).template init<UseWholeIterator,CallByPtr>(); */

  /*   while(!field.at_end()) { */
  /*     proc(field.ptr()); */
  /*     field.next(); */
  /*   }; */
  /* }; */

  //} // namespace structured
} // namespace SpatialOps

#include <spatialops/CoreMacrosUndefine.h>

#endif // SpatialOps_FieldTill_h
