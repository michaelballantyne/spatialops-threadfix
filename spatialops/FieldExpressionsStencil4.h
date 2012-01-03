#ifndef SpatialOps_FieldExpressionsStencil_4_h
#  define SpatialOps_FieldExpressionsStencil_4_h

#  include <spatialops/SpatialOpsConfigure.h>
#  include <spatialops/FieldExpressions.h>
#  include <spatialops/structured/stencil/Stencil4.h>

   /*#include <iostream> */

#  ifdef STENCIL_THREADS
#     include <vector>
#     include <boost/bind.hpp>
#     include <boost/ref.hpp>
#     include <spatialops/ThreadPool.h>
#     include <spatialops/structured/IntVec.h>
#     include <boost/interprocess/sync/interprocess_semaphore.hpp>
      namespace BI = boost::interprocess;
#  endif //STENCIL_THREADS

   namespace SpatialOps {

      namespace structured {

         template<typename OperatorType, typename SrcType, typename DestType>
          inline void stencil_4_apply_to_field_sequential_execute_internal (SrcType const & src,
                                                                            DestType & dest,
                                                                            double const coef1,
                                                                            double const coef2,
                                                                            double const coef3,
                                                                            double const coef4) {

             s4detail::ExtentsAndOffsets<SrcType, DestType> typedef Extents;

             const MemoryWindow & ws = src.window_with_ghost();

             const MemoryWindow ws1(ws.glob_dim(),
                                    ws.offset() + Extents::Src1Offset::int_vec(),
                                    ws.extent() + Extents::Src1Extent::int_vec() + ws.has_bc() *
                                    Extents::Src1ExtentBC::int_vec(),
                                    ws.has_bc(0),
                                    ws.has_bc(1),
                                    ws.has_bc(2));

             const MemoryWindow ws2(ws.glob_dim(),
                                    ws.offset() + Extents::Src2Offset::int_vec(),
                                    ws.extent() + Extents::Src2Extent::int_vec() + ws.has_bc() *
                                    Extents::Src2ExtentBC::int_vec(),
                                    ws.has_bc(0),
                                    ws.has_bc(1),
                                    ws.has_bc(2));

             const MemoryWindow ws3(ws.glob_dim(),
                                    ws.offset() + Extents::Src3Offset::int_vec(),
                                    ws.extent() + Extents::Src3Extent::int_vec() + ws.has_bc() *
                                    Extents::Src3ExtentBC::int_vec(),
                                    ws.has_bc(0),
                                    ws.has_bc(1),
                                    ws.has_bc(2));

             const MemoryWindow ws4(ws.glob_dim(),
                                    ws.offset() + Extents::Src4Offset::int_vec(),
                                    ws.extent() + Extents::Src4Extent::int_vec() + ws.has_bc() *
                                    Extents::Src4ExtentBC::int_vec(),
                                    ws.has_bc(0),
                                    ws.has_bc(1),
                                    ws.has_bc(2));

             const MemoryWindow & wdest = dest.window_with_ghost();

             const MemoryWindow wd(wdest.glob_dim(),
                                   wdest.offset() + Extents::DestOffset::int_vec(),
                                   wdest.extent() + Extents::DestExtent::int_vec() + wdest.has_bc()
                                   * Extents::DestExtentBC::int_vec(),
                                   wdest.has_bc(0),
                                   wdest.has_bc(1),
                                   wdest.has_bc(2));

#            ifndef NDEBUG
                assert(ws1.extent() == ws2.extent() && ws1.extent() == ws3.extent() && ws1.extent()
                       == ws4.extent() && ws1.extent() == wd.extent())
#            endif //NDEBUG;

             DestType d(wd, &dest[0], ExternalStorage);

             SrcType s1(ws1, &src[0], ExternalStorage);

             SrcType s2(ws2, &src[0], ExternalStorage);

             SrcType s3(ws3, &src[0], ExternalStorage);

             SrcType s4(ws4, &src[0], ExternalStorage);

             typename DestType::iterator id = d.begin();

             typename DestType::iterator ide = d.end();

             typename SrcType::const_iterator is1 = s1.begin();

             typename SrcType::const_iterator is2 = s2.begin();

             typename SrcType::const_iterator is3 = s3.begin();

             typename SrcType::const_iterator is4 = s4.begin();

             for(; id != ide; ++id, ++is1, ++is2, ++is3, ++is4){
                *id = *is1 * coef1 + *is2 * coef2 + *is3 * coef3 + *is4 * coef4;
             };
          };

         template<typename OperatorType, typename SrcType, typename DestType>
          inline void stencil_4_apply_to_field_sequential_execute (SrcType const & src,
                                                                   DestType & dest,
                                                                   double const coef1,
                                                                   double const coef2,
                                                                   double const coef3,
                                                                   double const coef4) {

             stencil_4_apply_to_field_sequential_execute_internal<OperatorType, SrcType, DestType>(src,
                                                                                                   dest,
                                                                                                   coef1,
                                                                                                   coef2,
                                                                                                   coef3,
                                                                                                   coef4);
          };

#        ifdef STENCIL_THREADS
            template<typename OperatorType, typename SrcType, typename DestType>
             inline void stencil_4_apply_to_field_thread_parallel_execute_internal (SrcType const &
                                                                                    src,
                                                                                    DestType & dest,
                                                                                    double const
                                                                                    coef1,
                                                                                    double const
                                                                                    coef2,
                                                                                    double const
                                                                                    coef3,
                                                                                    double const
                                                                                    coef4,
                                                                                    MemoryWindow
                                                                                    const & sw,
                                                                                    MemoryWindow
                                                                                    const & dw,
                                                                                    BI::
                                                                                    interprocess_semaphore
                                                                                    * sem) {

                stencil_4_apply_to_field_sequential_execute_internal<OperatorType,
                                                                     typename SrcType::field_type,
                                                                     typename DestType::field_type>(src.resize(sw).field(),
                                                                                                    dest.resize(dw).field(),
                                                                                                    coef1,
                                                                                                    coef2,
                                                                                                    coef3,
                                                                                                    coef4);

                sem->post();
             }
#        endif //STENCIL_THREADS;

#        ifdef STENCIL_THREADS
            template<typename OperatorType, typename SrcType, typename DestType>
             inline void stencil_4_apply_to_field_thread_parallel_execute (SrcType const & src,
                                                                           DestType & dest,
                                                                           double const coef1,
                                                                           double const coef2,
                                                                           double const coef3,
                                                                           double const coef4,
                                                                           int const
                                                                           number_of_partitions) {

                typename ConstField<Initial, SrcType>::template Iterator<UseWholeIterator>::
                ResizePrepType typedef SrcPmtrType;

                typename Field<Initial, DestType>::template Iterator<UseWholeIterator>::
                ResizePrepType typedef DestPmtrType;

                MemoryWindow sw = src.window_with_ghost();

                MemoryWindow dw = dest.window_with_ghost();

                int x = 1;
                int y = 1;
                int z = 1;
                int g = SrcType::Ghost::NGHOST;
                structured::IntVec gcs = structured::IntVec(g, g, g);

                if(number_of_partitions <= sw.extent(2)){ z = number_of_partitions; }
                else if(number_of_partitions <= sw.extent(1)){ y = number_of_partitions; }
                else if(number_of_partitions <= sw.extent(0)){ x = number_of_partitions; };

                std::vector<MemoryWindow> vec_sw = sw.split(structured::IntVec(x, y, z), gcs);

                std::vector<MemoryWindow> vec_dw = dw.split(structured::IntVec(x, y, z), gcs);

                std::vector<BI::interprocess_semaphore *> vec_semaphore;

                typename std::vector<MemoryWindow>::const_iterator ivec_sw = vec_sw.begin();

                typename std::vector<MemoryWindow>::const_iterator ivec_dw = vec_dw.begin();

                typename std::vector<MemoryWindow>::const_iterator evec_sw = vec_sw.end();

                for(; ivec_sw != evec_sw; ++ivec_sw, ++ivec_dw){

                   vec_semaphore.push_back(new BI::interprocess_semaphore(0));

                   ThreadPoolFIFO::self().schedule(boost::bind(&
                                                               stencil_4_apply_to_field_thread_parallel_execute_internal<OperatorType,
                                                                                                                         SrcPmtrType,
                                                                                                                         DestPmtrType>,
                                                               ConstField<Initial, SrcType>(src).template
                                                                                                 resize_prep<UseWholeIterator>(),
                                                               Field<Initial, DestType>(dest).template
                                                                                              resize_prep<UseWholeIterator>(),
                                                               coef1,
                                                               coef2,
                                                               coef3,
                                                               coef4,
                                                               *ivec_sw,
                                                               *ivec_dw,
                                                               vec_semaphore.back()));
                };

                std::vector<BI::interprocess_semaphore *>::iterator isem = vec_semaphore.begin();

                std::vector<BI::interprocess_semaphore *>::iterator const esem = vec_semaphore.end();

                for(; isem != esem; ++isem){ (*isem)->wait(); };

                for(isem = vec_semaphore.begin(); isem != esem; ++isem){ delete *isem; };
             }
#        endif //STENCIL_THREADS;

         template<typename OperatorType, typename SrcType, typename DestType>
          inline void stencil_4_apply_to_field_general_execute (SrcType const & src,
                                                                DestType & dest,
                                                                double const coef1,
                                                                double const coef2,
                                                                double const coef3,
                                                                double const coef4) {

#            ifdef STENCIL_THREADS
                stencil_4_apply_to_field_thread_parallel_execute<OperatorType, SrcType, DestType>(src,
                                                                                                  dest,
                                                                                                  coef1,
                                                                                                  coef2,
                                                                                                  coef3,
                                                                                                  coef4,
                                                                                                  NTHREADS)
#            else
                stencil_4_apply_to_field_sequential_execute<OperatorType, SrcType, DestType>(src,
                                                                                             dest,
                                                                                             coef1,
                                                                                             coef2,
                                                                                             coef3,
                                                                                             coef4)
#            endif //STENCIL_THREADS
             ;
          };
      } //structured;
   } //SpatialOps;

#endif //SpatialOps_FieldExpressionsStencil_4_h
