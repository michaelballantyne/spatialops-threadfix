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

#ifndef SpatialOps_FieldExpressionsStencil_4_h
#  define SpatialOps_FieldExpressionsStencil_4_h

#  include <spatialops/SpatialOpsConfigure.h>
#  include <spatialops/FieldExpressions.h>
#  include <spatialops/structured/stencil/Stencil4.h>

   /*#include <iostream> */

#  ifdef STENCIL_THREADS
#     include <vector>
#     include <boost/bind.hpp>
#     include <spatialops/ThreadPool.h>
#     include <spatialops/structured/IntVec.h>
#     include <boost/interprocess/sync/interprocess_semaphore.hpp>
      namespace BI = boost::interprocess;
#  endif
   /* STENCIL_THREADS */

   namespace SpatialOps {

      namespace structured {

         template<typename OperatorType, typename SrcType, typename DestType>
          inline void stencil_4_apply_to_field_sequential_execute_internal(SrcType const & src,
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
                       == ws4.extent() && ws1.extent() == wd.extent());
#            endif
             /* NDEBUG */;

             DestType d(wd, dest.field_values(), ExternalStorage);

             SrcType s1(ws1, src.field_values(), ExternalStorage);

             SrcType s2(ws2, src.field_values(), ExternalStorage);

             SrcType s3(ws3, src.field_values(), ExternalStorage);

             SrcType s4(ws4, src.field_values(), ExternalStorage);

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
          inline void stencil_4_apply_to_field_sequential_execute(SrcType const & src,
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
             inline void stencil_4_apply_to_field_thread_parallel_execute_internal(SrcType const &
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
#        endif
         /* STENCIL_THREADS */;

#        ifdef STENCIL_THREADS
            template<typename OperatorType, typename SrcType, typename DestType>
             inline void stencil_4_apply_to_field_thread_parallel_execute(SrcType const & src,
                                                                          DestType & dest,
                                                                          double const coef1,
                                                                          double const coef2,
                                                                          double const coef3,
                                                                          double const coef4,
                                                                          int const
                                                                          number_of_partitions) {

                typename structured::FromGhost<typename SrcType::Ghost>::result typedef ValidGhost;

                structured::IndexTriplet<0, 0, 0> typedef InitialShift;

                typename NeboConstField<Initial, SrcType>::template Iterator<ValidGhost,
                                                                             InitialShift>::
                ResizePrepType typedef SrcPmtrType;

                typename NeboField<Initial, DestType>::template Iterator<ValidGhost, InitialShift>::
                ResizePrepType typedef DestPmtrType;

                MemoryWindow sw = src.window_with_ghost();

                MemoryWindow dw = dest.window_with_ghost();

                int x = 1;
                int y = 1;
                int z = 1;

                if(number_of_partitions <= sw.extent(2)){ z = number_of_partitions; }
                else if(number_of_partitions <= sw.extent(1)){ y = number_of_partitions; }
                else if(number_of_partitions <= sw.extent(0)){ x = number_of_partitions; };

                typename SrcType::field_type::Location::BCExtra typedef SrcBCExtra;
                typename DestType::field_type::Location::BCExtra typedef DestBCExtra;
                structured::IntVec sBC = sw.has_bc() * SrcBCExtra::int_vec();
                structured::IntVec dBC = dw.has_bc() * DestBCExtra::int_vec();

                std::vector<MemoryWindow> vec_sw = sw.split(structured::IntVec(x, y, z),
                                                            SrcType::Ghost::NGhostMinus::int_vec(),
                                                            SrcType::Ghost::NGhostPlus::int_vec(),
                                                            sBC);

                std::vector<MemoryWindow> vec_dw = dw.split(structured::IntVec(x, y, z),
                                                            SrcType::Ghost::NGhostMinus::int_vec(),
                                                            SrcType::Ghost::NGhostPlus::int_vec(),
                                                            dBC);

                BI::interprocess_semaphore semaphore(0);

                typename std::vector<MemoryWindow>::const_iterator ivec_sw = vec_sw.begin();

                typename std::vector<MemoryWindow>::const_iterator ivec_dw = vec_dw.begin();

                typename std::vector<MemoryWindow>::const_iterator evec_sw = vec_sw.end();

                for(; ivec_sw != evec_sw; ++ivec_sw, ++ivec_dw){

                   ThreadPoolFIFO::self().schedule(boost::bind(&
                                                               stencil_4_apply_to_field_thread_parallel_execute_internal<OperatorType,
                                                                                                                         SrcPmtrType,
                                                                                                                         DestPmtrType>,
                                                               NeboConstField<Initial, SrcType>(src).template
                                                                                                     resize_prep<ValidGhost,
                                                                                                                 InitialShift>(),
                                                               NeboField<Initial, DestType>(dest).template
                                                                                                  resize_prep<ValidGhost,
                                                                                                              InitialShift>(),
                                                               coef1,
                                                               coef2,
                                                               coef3,
                                                               coef4,
                                                               *ivec_sw,
                                                               *ivec_dw,
                                                               &semaphore));
                };

                for(int ii = 0; ii < vec_sw.size(); ii++){ semaphore.wait(); };
             }
#        endif
         /* STENCIL_THREADS */;

         template<typename OperatorType, typename SrcType, typename DestType>
          inline void stencil_4_apply_to_field_general_execute(SrcType const & src,
                                                               DestType & dest,
                                                               double const coef1,
                                                               double const coef2,
                                                               double const coef3,
                                                               double const coef4) {

#            ifdef STENCIL_THREADS
                (is_thread_parallel() ? stencil_4_apply_to_field_thread_parallel_execute<OperatorType,
                                                                                              SrcType,
                                                                                              DestType>(src,
                                                                                                        dest,
                                                                                                        coef1,
                                                                                                        coef2,
                                                                                                        coef3,
                                                                                                        coef4,
                                                                                                        get_soft_thread_count())
                 : stencil_4_apply_to_field_sequential_execute<OperatorType, SrcType, DestType>(src,
                                                                                                dest,
                                                                                                coef1,
                                                                                                coef2,
                                                                                                coef3,
                                                                                                coef4))
#            else
                stencil_4_apply_to_field_sequential_execute<OperatorType, SrcType, DestType>(src,
                                                                                             dest,
                                                                                             coef1,
                                                                                             coef2,
                                                                                             coef3,
                                                                                             coef4)
#            endif
             /* STENCIL_THREADS */
             ;
          };
      } /* structured */;
   } /* SpatialOps */;

#endif
/* SpatialOps_FieldExpressionsStencil_4_h */
