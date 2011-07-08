#ifndef SpatialOps_FieldExpressionsStencils_h
#define SpatialOps_FieldExpressionsStencils_h

//#include <iostream>

#include <boost/bind.hpp>

#include <spatialops/SpatialOpsConfigure.h>

#ifdef STENCIL_THREADS
#  include <spatialops/ThreadPool.h>
#endif

#include <spatialops/FieldExpressionsMacros.h>
#include <spatialops/CoreMacrosDefine.h>

#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/stencil/Stencil2.h>
#include <spatialops/structured/stencil/Stencil4.h>

namespace SpatialOps{
  //namespace structured{
  /**
   *  @file FieldExpressionsStencils.h
   */
  
  template<typename IteratorType, typename SrcType, typename DestType>
    struct StencilGeneratePartitions {
      GeneratePartitions<IteratorType,SrcType>  typedef SrcGenerator;
      GeneratePartitions<IteratorType,DestType> typedef DestGenerator;
    private:
      SrcGenerator  srcPart;   //partition for Source      Field
      DestGenerator destPart;  //partition for Destination Field
      
    public:
      StencilGeneratePartitions(SrcType  const & src ,
				 DestType const & dest,
				 int      const   n    )
      {
	/*find axis*/
	int axis = -1;
	typename DestType::memory_window::dimension_type const extent = IterFcnsStyle<IteratorType,DestType>::memory_window(dest).extent();
	for(int j = 2; j >= 0 && axis == -1; j--)
	  if(extent[j] >= n)
	    axis = j;
	
	/*find any physcial boundaries*/
	int src_extent  = IterFcnsStyle<IteratorType,SrcType> ::memory_window(src) .extent(axis);
	int dest_extent = IterFcnsStyle<IteratorType,DestType>::memory_window(dest).extent(axis);
	if(  src_extent != dest_extent)
	  if(src_extent <  dest_extent)
	    {
	      /* dest has physical boundary */
	      srcPart  = SrcGenerator (axis, n, false);
	      destPart = DestGenerator(axis, n, true );
	    }
	  else
	    {
	      /* src has physical boundary */
	      srcPart  = SrcGenerator (axis, n, true );
	      destPart = DestGenerator(axis, n, false);
	    };
      };
      
      I SrcGenerator  const & srcGen () const { return srcPart;  };
      I DestGenerator const & destGen() const { return destPart; };
    };
  
  /*sequential apply_to_field call*/
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      stencil_2_apply_to_field_sequential_execute_internal
		      ,/*Template Parameters*/
		      typename OperatorT QC
		      typename SrcT      QC
		      typename DestT
		      ,/*Function Parameters*/
		      SrcT   const & src  QC
		      DestT        & dest QC
		      double const   low  QC
		      double const   high 
		      ,/*Function's Body*/
		      typedef typename structured::GetConstT<SrcT>::type ConstSrcT;
		      /* typedef typename SrcT::const_field_type ConstSrcT; */
		      
		      const structured::Stencil2Helper<ConstSrcT QC DestT> helper( src .window_with_ghost() QC
										   dest.window_with_ghost()  );
		      
		      const structured::IntVec sinc = helper.src_increment ();
		      const structured::IntVec dinc = helper.dest_increment();
		      
		      typename DestT::iterator       idest = dest.begin() + helper.dest_offset  ();
		      typename SrcT ::const_iterator isrcm = src .begin() + helper.src_offset_lo();
		      typename SrcT ::const_iterator isrcp = src .begin() + helper.src_offset_hi();
		      
		      const structured::IntVec lo = helper.low ();
		      const structured::IntVec hi = helper.high();
		      
		      for( int k=lo[2]; k<hi[2]; ++k ){
			for( int j=lo[1]; j<hi[1]; ++j ){
			  for( int i=lo[0]; i<hi[0]; ++i ){
			    *idest = low * *isrcm + high * *isrcp;
			    idest += dinc[0];
			    isrcm += sinc[0];
			    isrcp += sinc[0];
			  }
			  idest += dinc[1];
			  isrcm += sinc[1];
			  isrcp += sinc[1];
			}
			idest += dinc[2];
			isrcm += sinc[2];
			isrcp += sinc[2];
		      }
		      );
  
  /*sequential external call*/
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      stencil_2_apply_to_field_sequential_execute
		      ,/*Template Parameters*/
		      typename CallStyle QC
		      typename OperatorT QC
		      typename SrcT      QC
		      typename DestT
		      ,/*Function Parameters*/
		      SrcT   const & src  QC
		      DestT        & dest QC
		      double const   low  QC
		      double const   high
		      ,/*Function's Body*/
		      stencil_2_apply_to_field_sequential_execute_internal<
		      OperatorT QC
		      SrcT      QC
		      DestT      >
		      (src  QC
		       dest QC
		       low  QC
		       high  )
		      );
  
#ifdef STENCIL_THREADS
  /*set-up for sequential loop in a single thread*/
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      stencil_2_apply_to_field_thread_parallel_execute_internal
		      ,/*Template Parameters*/
		      typename CallStyle QC
		      typename OperatorT QC
		      typename SrcT      QC
		      typename DestT
		      ,/*Function Parameters*/
		      SrcT         const & src           QC
		      DestT              & dest          QC
		      double       const   low           QC
		      double       const   high          QC
		      NewPartition const & src_partition QC
		      NewPartition const & dest_partition
		      ,/*Function's Body*/
		      stencil_2_apply_to_field_sequential_execute_internal<
		      OperatorT                      QC
		      typename SrcT ::ConstFieldType QC
		      typename DestT::   NCFieldType  >
		      (src. resize(src_partition) .field() QC
		       dest.resize(dest_partition).field() QC
		       low                                 QC
		       high                                 )
		      );
#endif
  
#ifdef STENCIL_THREADS
  /*set-up for all threads*/
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      stencil_2_apply_to_field_thread_parallel_execute
		      ,/*Template Parameters*/
		      typename CallStyle QC
		      typename OperatorT QC
		      typename SrcT      QC
		      typename DestT
		      ,/*Function Parameters*/
		      SrcT   const & src               QC
		      DestT        & dest              QC
		      double const   low               QC
		      double const   high              QC
		      int    const   number_of_threads
		      ,/*Function's Body*/
		      typename ConstFieldForm<Initial QC SrcT> ::template Iterator<CallStyle>::ResizePrepType typedef FirstType;
		      typename      FieldForm<Initial QC DestT>::template Iterator<CallStyle>::ResizePrepType typedef SecondType;
		      StencilGeneratePartitions<CallStyle QC
		                                SrcT      QC
                                                DestT      >                                                  typedef SGPType;
		      typename SGPType::SrcGenerator                                                          typedef SrcGPType;
		      typename SGPType::DestGenerator                                                         typedef DestGPType;

		      SGPType    gen_gen_par    (src QC dest QC number_of_threads);
		      SrcGPType  src_gen_par  = gen_gen_par.srcGen ();
		      DestGPType dest_gen_par = gen_gen_par.destGen();
		      for(; !dest_gen_par.at_last(); ++src_gen_par, ++dest_gen_par) {
			ThreadPoolFIFO::self().schedule(boost::bind(&stencil_2_apply_to_field_thread_parallel_execute_internal<
								CallStyle  QC
								OperatorT  QC
								FirstType  QC
								SecondType  >                                                            QC
								ConstFieldForm<Initial QC SrcT> (src) .template resize_prep<CallStyle>() QC
								     FieldForm<Initial QC DestT>(dest).template resize_prep<CallStyle>() QC
								low                                                                      QC
								high                                                                     QC
								src_gen_par .partition                                                () QC
								dest_gen_par.partition                                                ()  )); };
		      ThreadPoolFIFO::self().wait()
		      );
#endif
  
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      stencil_2_apply_to_field_general_execute
		      ,/*Template Parameters*/
		      typename CallStyle QC
		      typename OperatorT QC
		      typename SrcT      QC
		      typename DestT
		      ,/*Function Parameters*/
		      SrcT   const & src  QC
		      DestT        & dest QC
		      double const   low  QC
		      double const   high
		      ,/*Function's Body*/
#ifdef STENCIL_THREADS
		      stencil_2_apply_to_field_thread_parallel_execute<
		      CallStyle QC
		      OperatorT QC
		      SrcT      QC
		      DestT      >
		      (src     ,
		       dest    ,
		       low     ,
		       high    ,
		       NTHREADS )
#else
		      stencil_2_apply_to_field_sequential_execute<
		      CallStyle QC
		      OperatorT QC
		      SrcT      QC
		      DestT      >
		      (src ,
		       dest,
		       low ,
		       high )
#endif
		      );
  
  /*sequential apply_to_field call*/
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      stencil_4_apply_to_field_sequential_execute_internal
		      ,/*Template Parameters*/
		      typename OperatorT QC
		      typename SrcT      QC
		      typename DestT
		      ,/*Function Parameters*/
		      SrcT   const & src   QC
		      DestT        & dest  QC
		      double const   coef1 QC
		      double const   coef2 QC
		      double const   coef3 QC
		      double const   coef4
		      ,/*Function's Body*/
		      typedef typename structured::GetConstT<SrcT>::type ConstSrcT;
		      const structured::Stencil4Helper<ConstSrcT QC DestT> helper( src.window_with_ghost(),
										   dest.window_with_ghost() );
		      
		      const structured::IntVec sinc = helper.src_increment();
		      const structured::IntVec dinc = helper.dest_increment();
		      
		      typename  DestT::iterator idest = dest.begin() + helper.dest_offset();
		      
		      typedef typename SrcT::const_iterator SrcIter;
		      SrcIter isrc1 = src.begin() + helper.src_offset_1();
		      SrcIter isrc2 = src.begin() + helper.src_offset_2();
		      SrcIter isrc3 = src.begin() + helper.src_offset_3();
		      SrcIter isrc4 = src.begin() + helper.src_offset_4();
		      
		      const structured::IntVec lo = helper.low ();
		      const structured::IntVec hi = helper.high();
		      
		      for( int k=lo[2]; k<hi[2]; ++k ){
			for( int j=lo[1]; j<hi[1]; ++j ){
			  for( int i=lo[0]; i<hi[0]; ++i ){
			    *idest = coef1 * *isrc1
			           + coef2 * *isrc2
			           + coef3 * *isrc3
			           + coef4 * *isrc4;
			    idest += dinc[0];
			    isrc1 += sinc[0];
			    isrc2 += sinc[0];
			    isrc3 += sinc[0];
			    isrc4 += sinc[0];
			  }
			  idest += dinc[1];
			  isrc1 += sinc[1];
			  isrc2 += sinc[1];
			  isrc3 += sinc[1];
			  isrc4 += sinc[1];
			}
			idest += dinc[2];
			isrc1 += sinc[2];
			isrc2 += sinc[2];
			isrc3 += sinc[2];
			isrc4 += sinc[2];
		      }
		      );
  
  /*sequential external call*/
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      stencil_4_apply_to_field_sequential_execute
		      ,/*Template Parameters*/
		      typename CallStyle QC
		      typename OperatorT QC
		      typename SrcT      QC
		      typename DestT
		      ,/*Function Parameters*/
		      SrcT   const & src   QC
		      DestT        & dest  QC
		      double const   coef1 QC
		      double const   coef2 QC
		      double const   coef3 QC
		      double const   coef4
		      ,/*Function's Body*/
		      stencil_4_apply_to_field_sequential_execute_internal<
		      OperatorT QC
		      SrcT      QC
		      DestT      >
		      (src   QC
		       dest  QC
		       coef1 QC
		       coef2 QC
		       coef3 QC
		       coef4  )
		      );
  
#ifdef STENCIL_THREADS
  /*set-up for sequential loop in a single thread*/
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      stencil_4_apply_to_field_thread_parallel_execute_internal
		      ,/*Template Parameters*/
		      typename CallStyle QC
		      typename OperatorT QC
		      typename SrcT      QC
		      typename DestT
		      ,/*Function Parameters*/
		      SrcT         const & src           QC
		      DestT              & dest          QC
		      double       const   coef1         QC
		      double       const   coef2         QC
		      double       const   coef3         QC
		      double       const   coef4         QC
		      NewPartition const & src_partition QC
		      NewPartition const & dest_partition
		      ,/*Function's Body*/
		      stencil_4_apply_to_field_sequential_execute_internal<
		      OperatorT                      QC
		      typename SrcT ::ConstFieldType QC
		      typename DestT::   NCFieldType  >
		      (src. resize(src_partition) .field() QC
		       dest.resize(dest_partition).field() QC
		       coef1                               QC
		       coef2                               QC
		       coef3                               QC
		       coef4                                )
		      );
#endif
  
#ifdef STENCIL_THREADS
  /*set-up for all threads*/
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      stencil_4_apply_to_field_thread_parallel_execute
		      ,/*Template Parameters*/
		      typename CallStyle QC
		      typename OperatorT QC
		      typename SrcT      QC
		      typename DestT
		      ,/*Function Parameters*/
		      SrcT   const & src               QC
		      DestT        & dest              QC
		      double const   coef1             QC
		      double const   coef2             QC
		      double const   coef3             QC
		      double const   coef4             QC
		      int    const   number_of_threads
		      ,/*Function's Body*/
		      typename ConstFieldForm<Initial QC SrcT> ::template Iterator<CallStyle>::ResizePrepType typedef FirstType;
		      typename      FieldForm<Initial QC DestT>::template Iterator<CallStyle>::ResizePrepType typedef SecondType;
		      StencilGeneratePartitions<CallStyle QC
		                                SrcT      QC
                                                DestT      >                                                  typedef SGPType;
		      typename SGPType::SrcGenerator                                                          typedef SrcGPType;
		      typename SGPType::DestGenerator                                                         typedef DestGPType;

		      SGPType    gen_gen_par    (src QC dest QC number_of_threads);
		      SrcGPType  src_gen_par  = gen_gen_par.srcGen ();
		      DestGPType dest_gen_par = gen_gen_par.destGen();
		      for(; !dest_gen_par.at_last(); ++src_gen_par, ++dest_gen_par) {
			ThreadPoolFIFO::self().schedule(boost::bind(&stencil_4_apply_to_field_thread_parallel_execute_internal<
								CallStyle  QC
								OperatorT  QC
								FirstType  QC
								SecondType  >                                                            QC
								ConstFieldForm<Initial QC SrcT> (src) .template resize_prep<CallStyle>() QC
								     FieldForm<Initial QC DestT>(dest).template resize_prep<CallStyle>() QC
								coef1                                                                    QC
								coef2                                                                    QC
								coef3                                                                    QC
								coef4                                                                    QC
								src_gen_par .partition                                                () QC
								dest_gen_par.partition                                                ()  )); };
		      ThreadPoolFIFO::self().wait()
		      );
#endif
  
  BUILD_VOID_FUNCTION(/*Function's Name*/
		      stencil_4_apply_to_field_general_execute
		      ,/*Template Parameters*/
		      typename CallStyle QC
		      typename OperatorT QC
		      typename SrcT      QC
		      typename DestT
		      ,/*Function Parameters*/
		      SrcT   const & src   QC
		      DestT        & dest  QC
		      double const   coef1 QC
		      double const   coef2 QC
		      double const   coef3 QC
		      double const   coef4
		      ,/*Function's Body*/
#ifdef STENCIL_THREADS
		      stencil_4_apply_to_field_thread_parallel_execute<
		      CallStyle QC
		      OperatorT QC
		      SrcT      QC
		      DestT      >
		      (src     ,
		       dest    ,
		       coef1   ,
		       coef2   ,
		       coef3   ,
		       coef4   ,
		       NTHREADS )
#else
		      stencil_4_apply_to_field_sequential_execute<
		      CallStyle QC
		      OperatorT QC
		      SrcT      QC
		      DestT      >
		      (src ,
		       dest,
		       coef1,
		       coef2,
		       coef3,
		       coef4 )
#endif
		      );
  
  //} // namespace structured
} // namespace SpatialOps

#include <spatialops/CoreMacrosUndefine.h>

#endif // SpatialOps_FieldExpressionsStencils_h
