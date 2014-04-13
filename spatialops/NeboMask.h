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

#ifndef NEBO_MASK_H
#  define NEBO_MASK_H

   namespace SpatialOps {
      template<typename CurrentMode, typename FieldType>
       struct NeboMask;
      template<typename FieldType>
       struct NeboMask<Initial, FieldType> {
         public:
          FieldType typedef field_type;

          NeboMask<SeqWalk, FieldType> typedef SeqWalkType;

#         ifdef FIELD_EXPRESSION_THREADS
             NeboMask<Resize, FieldType> typedef ResizeType;
#         endif
          /* FIELD_EXPRESSION_THREADS */

#         ifdef __CUDACC__
             NeboMask<GPUWalk, FieldType> typedef GPUWalkType;
#         endif
          /* __CUDACC__ */

          NeboMask<Reduction, FieldType> typedef ReductionType;

          NeboMask(structured::SpatialMask<FieldType> const & m)
          : mask_(m)
          {}

          inline structured::GhostData possible_ghosts(void) const {
             return mask_.get_valid_ghost_data() + point_to_ghost(mask_.boundary_info().has_extra());
          }

          inline SeqWalkType init(void) const { return SeqWalkType(mask_); }

#         ifdef FIELD_EXPRESSION_THREADS
             inline ResizeType resize(void) const { return ResizeType(mask_); }
#         endif
          /* FIELD_EXPRESSION_THREADS */

#         ifdef __CUDACC__
             inline bool cpu_ready(void) const {
                return mask_.find_consumer(LOCAL_RAM, 0);
             }

             inline bool gpu_ready(int const deviceIndex) const {
                return mask_.find_consumer(EXTERNAL_CUDA_GPU, deviceIndex);
             }

             inline GPUWalkType gpu_init(int const deviceIndex) const {
                return GPUWalkType(deviceIndex, mask_);
             }

#            ifdef NEBO_GPU_TEST
                inline void gpu_prep(int const deviceIndex) const {
                   const_cast<structured::SpatialMask<FieldType> *>(&mask_)->
                   add_consumer(EXTERNAL_CUDA_GPU, deviceIndex);
                }
#            endif
             /* NEBO_GPU_TEST */
#         endif
          /* __CUDACC__ */

          inline ReductionType reduce_init(void) const {
             return ReductionType(mask_);
          }

         private:
          structured::SpatialMask<FieldType> const mask_;
      };
#     ifdef FIELD_EXPRESSION_THREADS
         template<typename FieldType>
          struct NeboMask<Resize, FieldType> {
            public:
             FieldType typedef field_type;

             NeboMask<SeqWalk, FieldType> typedef SeqWalkType;

             NeboMask(structured::SpatialMask<FieldType> const & m)
             : mask_(m)
             {}

             inline SeqWalkType init(void) const { return SeqWalkType(mask_); }

            private:
             structured::SpatialMask<FieldType> const mask_;
         }
#     endif
      /* FIELD_EXPRESSION_THREADS */;
      template<typename FieldType>
       struct NeboMask<SeqWalk, FieldType> {
         public:
          FieldType typedef field_type;

          typename field_type::value_type typedef value_type;

          NeboMask(structured::SpatialMask<FieldType> const & m)
          : mask_(m)
          {}

          inline value_type eval(int const x, int const y, int const z) const {
             return mask_(x, y, z);
          }

         private:
          structured::SpatialMask<FieldType> const mask_;
      };
#     ifdef __CUDACC__
         template<typename FieldType>
          struct NeboMask<GPUWalk, FieldType> {
            public:
             FieldType typedef field_type;

             typename field_type::value_type typedef value_type;

             NeboMask(int const deviceIndex,
                      structured::SpatialMask<FieldType> const & m)
             : bitField_(m.mask_values(EXTERNAL_CUDA_GPU, deviceIndex)),
               xOffset_(m.window_with_ghost().offset(0) + m.get_ghost_data().get_minus(0)),
               yOffset_(m.window_with_ghost().offset(1) + m.get_ghost_data().get_minus(1)),
               zOffset_(m.window_with_ghost().offset(2) + m.get_ghost_data().get_minus(2)),
               xGlob_(m.window_with_ghost().glob_dim(0)),
               yGlob_(m.window_with_ghost().glob_dim(1))
             {}

             __device__ inline bool eval(int const x, int const y, int const z) const {
                return deref(x, y, z);
             }

            private:
             __device__ inline int find_position(int const x,
                                                 int const y,
                                                 int const z) const {
                const int newX = xOffset_ + x;

                const int newY = yOffset_ + y;

                const int newZ = zOffset_ + z;

                return newX + xGlob_ * (newY + yGlob_ * newZ);
             }

             __device__ inline int find_block(int const position) const {
                return position / NEBO_INT_BIT;
             }

             __device__ inline int find_bit_position(int const position) const {
                return position % NEBO_INT_BIT;
             }

             __device__ inline int deref(int const x, int const y, int const z) const {
                const int position = find_position(x, y, z);

                return !(!(*(bitField_ + find_block(position)) & (1 <<
                                                                  find_bit_position(position))));
             }

             unsigned int const * bitField_;

             int const xOffset_;

             int const yOffset_;

             int const zOffset_;

             int const xGlob_;

             int const yGlob_;
         }
#     endif
      /* __CUDACC__ */;
      template<typename FieldType>
       struct NeboMask<Reduction, FieldType> {
         public:
          FieldType typedef field_type;

          typename field_type::value_type typedef value_type;

          NeboMask(structured::SpatialMask<FieldType> const & m)
          : iter_(m.begin()), end_(m.end())
          {}

          inline void next(void) { iter_++; }

          inline bool at_end(void) const { return iter_ == end_; }

          inline bool has_length(void) const { return true; }

          inline value_type eval(void) const { return *iter_; }

         private:
          typename structured::SpatialMask<FieldType>::const_iterator iter_;

          typename structured::SpatialMask<FieldType>::const_iterator const end_
          ;
      };
   } /* SpatialOps */

#endif
/* NEBO_MASK_H */
