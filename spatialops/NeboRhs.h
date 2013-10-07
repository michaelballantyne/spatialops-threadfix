/* This file was generated by fulmar version 0.8.0. */

/*
 * Copyright (c) 2013 The University of Utah
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

#ifndef NEBO_RHS_H
#  define NEBO_RHS_H

   namespace SpatialOps {
      template<typename CurrentMode, typename AtomicType>
       struct NeboScalar;
      template<typename AtomicType>
       struct NeboScalar<Initial, AtomicType> {
         public:
          AtomicType typedef value_type;

          NeboScalar<SeqWalk, AtomicType> typedef SeqWalkType;

#         ifdef FIELD_EXPRESSION_THREADS
             NeboScalar<Resize, AtomicType> typedef ResizeType;
#         endif
          /* FIELD_EXPRESSION_THREADS */

#         ifdef __CUDACC__
             NeboScalar<GPUWalk, AtomicType> typedef GPUWalkType;
#         endif
          /* __CUDACC__ */

          NeboScalar<Reduction, AtomicType> typedef ReductionType;

          NeboScalar(value_type const v)
          : value_(v)
          {}

          inline structured::GhostData possible_ghosts(void) const {
             return structured::GhostData(GHOST_MAX);
          }

          inline SeqWalkType init(structured::IntVec const & minus,
                                  structured::IntVec const & plus,
                                  structured::IntVec const & shift) const {
             return SeqWalkType(value_);
          }

#         ifdef FIELD_EXPRESSION_THREADS
             inline ResizeType resize(structured::IntVec const & minus,
                                      structured::IntVec const & plus) const {
                return ResizeType(value_);
             }
#         endif
          /* FIELD_EXPRESSION_THREADS */

#         ifdef __CUDACC__
             inline bool cpu_ready(void) const { return true; }

             inline bool gpu_ready(int const deviceIndex) const { return true; }

             inline GPUWalkType gpu_init(structured::IntVec const & minus,
                                         structured::IntVec const & plus,
                                         structured::IntVec const & shift,
                                         int const deviceIndex) const {
                return GPUWalkType(value_);
             }

#            ifdef NEBO_GPU_TEST
                inline void gpu_prep(int const deviceIndex) const {}
#            endif
             /* NEBO_GPU_TEST */
#         endif
          /* __CUDACC__ */

          inline ReductionType reduce_init(structured::IntVec const & minus,
                                           structured::IntVec const & plus,
                                           structured::IntVec const & shift) const {
             return ReductionType(value_);
          }

         private:
          value_type const value_;
      };
#     ifdef FIELD_EXPRESSION_THREADS
         template<typename AtomicType>
          struct NeboScalar<Resize, AtomicType> {
            public:
             AtomicType typedef value_type;

             NeboScalar<SeqWalk, AtomicType> typedef SeqWalkType;

             NeboScalar(value_type const value)
             : value_(value)
             {}

             inline SeqWalkType init(structured::IntVec const & shift,
                                     structured::IntVec const & split,
                                     structured::IntVec const & location) const {
                return SeqWalkType(value_);
             }

            private:
             value_type const value_;
         }
#     endif
      /* FIELD_EXPRESSION_THREADS */;
      template<typename AtomicType>
       struct NeboScalar<SeqWalk, AtomicType> {
         public:
          AtomicType typedef value_type;

          NeboScalar(value_type const value)
          : value_(value)
          {}

          inline void next(void) {}

          inline value_type eval(void) const { return value_; }

         private:
          value_type const value_;
      };
#     ifdef __CUDACC__
         template<typename AtomicType>
          struct NeboScalar<GPUWalk, AtomicType> {
            public:
             AtomicType typedef value_type;

             NeboScalar(value_type const value)
             : value_(value)
             {}

             __device__ inline void start(int x, int y) {}

             __device__ inline void next(void) {}

             __device__ inline value_type eval(void) const { return value_; }

            private:
             value_type const value_;
         }
#     endif
      /* __CUDACC__ */;
      template<typename AtomicType>
       struct NeboScalar<Reduction, AtomicType> {
         public:
          AtomicType typedef value_type;

          NeboScalar(value_type const value)
          : value_(value)
          {}

          inline void next(void) {}

          inline bool at_end(void) const { return false; }

          inline bool has_length(void) const { return false; }

          inline value_type eval(void) const { return value_; }

         private:
          value_type const value_;
      };

      template<typename CurrentMode, typename FieldType>
       struct NeboConstField;
      template<typename FieldType>
       struct NeboConstField<Initial, FieldType> {
         public:
          FieldType typedef field_type;

          NeboConstField<SeqWalk, FieldType> typedef SeqWalkType;

#         ifdef FIELD_EXPRESSION_THREADS
             NeboConstField<Resize, FieldType> typedef ResizeType;
#         endif
          /* FIELD_EXPRESSION_THREADS */

#         ifdef __CUDACC__
             NeboConstField<GPUWalk, FieldType> typedef GPUWalkType;
#         endif
          /* __CUDACC__ */

          NeboConstField<Reduction, FieldType> typedef ReductionType;

          NeboConstField(FieldType const & f)
          : field_(f)
          {}

          inline structured::GhostData possible_ghosts(void) const {
             return field_.get_valid_ghost_data() + point_to_ghost(field_.boundary_info().has_extra());
          }

          inline SeqWalkType init(structured::IntVec const & minus,
                                  structured::IntVec const & plus,
                                  structured::IntVec const & shift) const {
             return SeqWalkType(resize_ghost_and_shift_window(field_,
                                                              minus,
                                                              plus - field_.boundary_info().has_extra(),
                                                              shift));
          }

#         ifdef FIELD_EXPRESSION_THREADS
             inline ResizeType resize(structured::IntVec const & minus,
                                      structured::IntVec const & plus) const {
                return ResizeType(resize_ghost(field_, minus, plus - field_.boundary_info().has_extra()));
             }
#         endif
          /* FIELD_EXPRESSION_THREADS */

#         ifdef __CUDACC__
             inline bool cpu_ready(void) const {
                return field_.find_consumer(LOCAL_RAM, 0);
             }

             inline bool gpu_ready(int const deviceIndex) const {
                return field_.find_consumer(EXTERNAL_CUDA_GPU, deviceIndex);
             }

             inline GPUWalkType gpu_init(structured::IntVec const & minus,
                                         structured::IntVec const & plus,
                                         structured::IntVec const & shift,
                                         int const deviceIndex) const {
                return GPUWalkType(deviceIndex,
                                   resize_ghost_and_shift_window(field_,
                                                                 minus,
                                                                 plus - field_.boundary_info().has_extra(),
                                                                 shift));
             }

#            ifdef NEBO_GPU_TEST
                inline void gpu_prep(int const deviceIndex) const {
                   const_cast<FieldType *>(&field_)->add_consumer(EXTERNAL_CUDA_GPU,
                                                                  deviceIndex);
                }
#            endif
             /* NEBO_GPU_TEST */
#         endif
          /* __CUDACC__ */

          inline ReductionType reduce_init(structured::IntVec const & minus,
                                           structured::IntVec const & plus,
                                           structured::IntVec const & shift) const {
             return ReductionType(resize_ghost_and_shift_window(field_,
                                                                minus,
                                                                plus - field_.boundary_info().has_extra(),
                                                                shift));
          }

         private:
          FieldType const field_;
      };
#     ifdef FIELD_EXPRESSION_THREADS
         template<typename FieldType>
          struct NeboConstField<Resize, FieldType> {
            public:
             FieldType typedef field_type;

             NeboConstField<SeqWalk, FieldType> typedef SeqWalkType;

             NeboConstField(FieldType const & f)
             : field_(f)
             {}

             inline SeqWalkType init(structured::IntVec const & shift,
                                     structured::IntVec const & split,
                                     structured::IntVec const & location) const {
                return SeqWalkType(shift_window(FieldType(field_.window_with_ghost().refine(split,
                                                                                            location),
                                                          field_),
                                                shift));
             }

            private:
             FieldType const field_;
         }
#     endif
      /* FIELD_EXPRESSION_THREADS */;
      template<typename FieldType>
       struct NeboConstField<SeqWalk, FieldType> {
         public:
          FieldType typedef field_type;

          typename field_type::value_type typedef value_type;

          NeboConstField(FieldType const & f)
          : iter_(f.begin())
          {}

          inline void next(void) { iter_++; }

          inline value_type eval(void) const { return *iter_; }

         private:
          typename FieldType::const_iterator iter_;
      };
#     ifdef __CUDACC__
         template<typename FieldType>
          struct NeboConstField<GPUWalk, FieldType> {
            public:
             FieldType typedef field_type;

             typename field_type::value_type typedef value_type;

             NeboConstField(int const deviceIndex, FieldType const & f)
             : current_(f.field_values(EXTERNAL_CUDA_GPU, deviceIndex) + f.window_with_ghost().offset(0)
               + f.window_with_ghost().glob_dim(0) * (f.window_with_ghost().offset(1)
                                                      + (f.window_with_ghost().glob_dim(1)
                                                         * f.window_with_ghost().offset(2)))),
               xLength_(f.window_with_ghost().glob_dim(0)),
               step_(xLength_ * f.window_with_ghost().glob_dim(1))
             {}

             __device__ inline void start(int x, int y) {
                current_ += x + y * xLength_;
             }

             __device__ inline void next(void) { current_ += step_; }

             __device__ inline value_type eval(void) const { return *current_; }

            private:
             value_type const * current_;

             int const xLength_;

             int const step_;
         }
#     endif
      /* __CUDACC__ */;
      template<typename FieldType>
       struct NeboConstField<Reduction, FieldType> {
         public:
          FieldType typedef field_type;

          typename field_type::value_type typedef value_type;

          NeboConstField(FieldType const & f)
          : iter_(f.begin()), end_(f.end())
          {}

          inline void next(void) { iter_++; }

          inline bool at_end(void) const { return iter_ == end_; }

          inline bool has_length(void) const { return true; }

          inline value_type eval(void) const { return *iter_; }

         private:
          typename FieldType::const_iterator iter_;

          typename FieldType::const_iterator const end_;
      };

      template<typename CurrentMode, typename T>
       struct NeboConstSingleValueField;
      template<typename T>
       struct NeboConstSingleValueField<Initial, T> {
         public:
          SpatialOps::structured::SpatialField<SpatialOps::structured::
                                               SingleValue,
                                               T> typedef field_type;

          SpatialOps::structured::SpatialField<SpatialOps::structured::
                                               SingleValue,
                                               T> typedef SingleValueFieldType;

          NeboConstSingleValueField<SeqWalk, T> typedef SeqWalkType;

#         ifdef FIELD_EXPRESSION_THREADS
             NeboConstSingleValueField<Resize, T> typedef ResizeType;
#         endif
          /* FIELD_EXPRESSION_THREADS */

#         ifdef __CUDACC__
             NeboConstSingleValueField<GPUWalk, T> typedef GPUWalkType;
#         endif
          /* __CUDACC__ */

          NeboConstSingleValueField<Reduction, T> typedef ReductionType;

          NeboConstSingleValueField(SingleValueFieldType const & f)
          : field_(f)
          {}

          inline structured::GhostData possible_ghosts(void) const {
             return structured::GhostData(GHOST_MAX);
          }

          inline SeqWalkType init(structured::IntVec const & minus,
                                  structured::IntVec const & plus,
                                  structured::IntVec const & shift) const {
             return SeqWalkType(* field_.field_values(LOCAL_RAM));
          }

#         ifdef FIELD_EXPRESSION_THREADS
             inline ResizeType resize(structured::IntVec const & minus,
                                      structured::IntVec const & plus) const {
                return ResizeType(* field_.field_values(LOCAL_RAM));
             }
#         endif
          /* FIELD_EXPRESSION_THREADS */

#         ifdef __CUDACC__
             inline bool cpu_ready(void) const {
                return field_.find_consumer(LOCAL_RAM, 0);
             }

             inline bool gpu_ready(int const deviceIndex) const {
                return field_.find_consumer(EXTERNAL_CUDA_GPU, deviceIndex);
             }

             inline GPUWalkType gpu_init(structured::IntVec const & minus,
                                         structured::IntVec const & plus,
                                         structured::IntVec const & shift,
                                         int const deviceIndex) const {
                return GPUWalkType(deviceIndex, field_);
             }

#            ifdef NEBO_GPU_TEST
                inline void gpu_prep(int const deviceIndex) const {
                   const_cast<SingleValueFieldType *>(&field_)->add_consumer(EXTERNAL_CUDA_GPU,
                                                                             deviceIndex);
                }
#            endif
             /* NEBO_GPU_TEST */
#         endif
          /* __CUDACC__ */

          inline ReductionType reduce_init(structured::IntVec const & minus,
                                           structured::IntVec const & plus,
                                           structured::IntVec const & shift) const {
             return ReductionType(* field_.field_values(LOCAL_RAM));
          }

         private:
          SingleValueFieldType const field_;
      };
#     ifdef FIELD_EXPRESSION_THREADS
         template<typename T>
          struct NeboConstSingleValueField<Resize, T> {
            public:
             SpatialOps::structured::SpatialField<SpatialOps::structured::
                                                  SingleValue,
                                                  T> typedef field_type;

             NeboConstSingleValueField<SeqWalk, T> typedef SeqWalkType;

             NeboConstSingleValueField(double const & v)
             : value_(v)
             {}

             inline SeqWalkType init(structured::IntVec const & shift,
                                     structured::IntVec const & split,
                                     structured::IntVec const & location) const {
                return SeqWalkType(value_);
             }

            private:
             double const value_;
         }
#     endif
      /* FIELD_EXPRESSION_THREADS */;
      template<typename T>
       struct NeboConstSingleValueField<SeqWalk, T> {
         public:
          SpatialOps::structured::SpatialField<SpatialOps::structured::
                                               SingleValue,
                                               T> typedef field_type;

          typename field_type::value_type typedef value_type;

          NeboConstSingleValueField(double const & v)
          : value_(v)
          {}

          inline void next(void) {}

          inline value_type eval(void) const { return value_; }

         private:
          double value_;
      };
#     ifdef __CUDACC__
         template<typename T>
          struct NeboConstSingleValueField<GPUWalk, T> {
            public:
             SpatialOps::structured::SpatialField<SpatialOps::structured::
                                                  SingleValue,
                                                  T> typedef field_type;

             typename field_type::value_type typedef value_type;

             SpatialOps::structured::SpatialField<SpatialOps::structured::
                                                  SingleValue,
                                                  T> typedef
             SingleValueFieldType;

             NeboConstSingleValueField(int const deviceIndex,
                                       SingleValueFieldType const & f)
             : ptr_(f.field_values(EXTERNAL_CUDA_GPU, deviceIndex)), value_(0.0)
             {}

             __device__ inline void start(int x, int y) { value_ = *ptr_; }

             __device__ inline void next(void) {}

             __device__ inline value_type eval(void) const { return value_; }

            private:
             value_type const * const ptr_;

             value_type value_;
         }
#     endif
      /* __CUDACC__ */;
      template<typename T>
       struct NeboConstSingleValueField<Reduction, T> {
         public:
          SpatialOps::structured::SpatialField<SpatialOps::structured::
                                               SingleValue,
                                               T> typedef field_type;

          typename field_type::value_type typedef value_type;

          NeboConstSingleValueField(double const & v)
          : value_(v)
          {}

          inline void next(void) {}

          inline bool at_end(void) const { return false; }

          inline bool has_length(void) const { return false; }

          inline value_type eval(void) const { return value_; }

         private:
          double const value_;
      };
   } /* SpatialOps */

#endif
/* NEBO_RHS_H */
