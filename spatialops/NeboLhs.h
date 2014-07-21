/* This file was generated by fulmar version 0.9.2. */

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

#ifndef NEBO_LHS_H
   #define NEBO_LHS_H

   namespace SpatialOps {
      #ifdef __CUDACC__
         template<typename LhsType, typename RhsType>
          __global__ void gpu_assign_kernel(LhsType lhs,
                                            RhsType rhs,
                                            int const xLow,
                                            int const xHigh,
                                            int const yLow,
                                            int const yHigh,
                                            int const zLow,
                                            int const zHigh) {
             lhs.gpuwalk_assign(rhs, xLow, xHigh, yLow, yHigh, zLow, zHigh);
          }
      #endif
      /* __CUDACC__ */;

      template<typename CurrentMode, typename FieldType>
       struct NeboField;
      template<typename FieldType>
       struct NeboField<Initial, FieldType> {
         public:
          FieldType typedef field_type;

          NeboField<SeqWalk, FieldType> typedef SeqWalkType;

          #ifdef ENABLE_THREADS
             NeboField<Resize, FieldType> typedef ResizeType;
          #endif
          /* ENABLE_THREADS */

          #ifdef __CUDACC__
             NeboField<GPUWalk, FieldType> typedef GPUWalkType;
          #endif
          /* __CUDACC__ */

          NeboField(FieldType f)
          : field_(f)
          {}

          template<typename RhsType>
           inline void assign(bool const useGhost, RhsType rhs) {
              GhostData const ghosts = calculate_actual_ghost(useGhost,
                                                              field_.get_ghost_data(),
                                                              field_.boundary_info(),
                                                              rhs.ghosts_with_bc());

              /* field_.reset_valid_ghosts(ghosts) */;

              IntVec const extents = field_.window_with_ghost().extent() -
              field_.get_valid_ghost_data().get_minus() - field_.get_valid_ghost_data().get_plus();

              IntVec const hasBC = field_.boundary_info().has_bc();

              const GhostData limits = GhostData(- ghosts.get_minus(0),
                                                 extents[0] + ghosts.get_plus(0),
                                                 - ghosts.get_minus(1),
                                                 extents[1] + ghosts.get_plus(1),
                                                 - ghosts.get_minus(2),
                                                 extents[2] + ghosts.get_plus(2));

              #ifdef __CUDACC__
                 #ifdef NEBO_GPU_TEST
                    gpu_test_assign<RhsType>(rhs, extents, ghosts, hasBC, limits)
                 #else
                    if(gpu_ready()) {
                       if(rhs.gpu_ready(gpu_device_index())) {
                          gpu_assign<RhsType>(rhs,
                                              extents,
                                              ghosts,
                                              hasBC,
                                              limits);
                       }
                       else {
                          std::ostringstream msg;
                          msg << "Nebo error in " << "Nebo Assignment" << ":\n";
                          msg << "Left-hand side of assignment allocated on ";
                          msg << "GPU but right-hand side is not ";
                          msg << "(completely) accessible on the same GPU";
                          msg << "\n";
                          msg << "\t - " << __FILE__ << " : " << __LINE__;
                          throw(std::runtime_error(msg.str()));
                       };
                    }
                    else {
                       if(cpu_ready()) {
                          if(rhs.cpu_ready()) {
                             cpu_assign<RhsType>(rhs,
                                                 extents,
                                                 ghosts,
                                                 hasBC,
                                                 limits);
                          }
                          else {
                             std::ostringstream msg;
                             msg << "Nebo error in " << "Nebo Assignment" <<
                             ":\n";
                             msg << "Left-hand side of assignment allocated on ";
                             msg << "CPU but right-hand side is not ";
                             msg << "(completely) accessible on the same CPU";
                             msg << "\n";
                             msg << "\t - " << __FILE__ << " : " << __LINE__;
                             throw(std::runtime_error(msg.str()));
                          };
                       }
                       else {
                          std::ostringstream msg;
                          msg << "Nebo error in " << "Nebo Assignment" << ":\n";
                          msg << "Left-hand side of assignment allocated on ";
                          msg << "unknown device - not on CPU or GPU";
                          msg << "\n";
                          msg << "\t - " << __FILE__ << " : " << __LINE__;
                          throw(std::runtime_error(msg.str()));
                       };
                    }
                 #endif
                 /* NEBO_GPU_TEST */
              #else
                 cpu_assign<RhsType>(rhs, extents, ghosts, hasBC, limits)
              #endif
              /* __CUDACC__ */;
           }

          template<typename RhsType>
           inline void masked_assign(SpatialMask<FieldType> const & mask,
                                     RhsType rhs) {
              #ifdef NEBO_REPORT_BACKEND
                 std::cout << "Starting Nebo masked assignment" << std::endl
              #endif
              /* NEBO_REPORT_BACKEND */;

              SeqWalkType lhs = init();

              typename RhsType::SeqWalkType expr = rhs.init(IntVec(0, 0, 0),
                                                            GhostData(0),
                                                            IntVec(0, 0, 0));

              std::vector<IntVec>::const_iterator ip = mask.points().begin();

              std::vector<IntVec>::const_iterator const ep = mask.points().end();

              for(; ip != ep; ip++) {
                 int const x = (*ip)[0];

                 int const y = (*ip)[1];

                 int const z = (*ip)[2];

                 lhs.ref(x, y, z) = expr.eval(x, y, z);
              };

              #ifdef __CUDACC__
                 if(gpu_ready()) {
                    std::ostringstream msg;
                    msg << "Nebo error in " << "Nebo Masked Assignment" << ":\n"
                    ;
                    msg << "Left-hand side of masked assignment allocated on ";
                    msg << "GPU and this backend does not support GPU execution"
                    ;
                    msg << "\n";
                    msg << "\t - " << __FILE__ << " : " << __LINE__;
                    throw(std::runtime_error(msg.str()));
                 }
                 else {
                    if(cpu_ready()) {
                       if(rhs.cpu_ready()) {
                          SeqWalkType lhs = init();

                          typename RhsType::SeqWalkType expr = rhs.init(IntVec(0,
                                                                               0,
                                                                               0),
                                                                        GhostData(0),
                                                                        IntVec(0,
                                                                               0,
                                                                               0));

                          std::vector<IntVec>::const_iterator ip = mask.points().begin();

                          std::vector<IntVec>::const_iterator const ep = mask.points().end();

                          for(; ip != ep; ip++) {
                             int const x = (*ip)[0];

                             int const y = (*ip)[1];

                             int const z = (*ip)[2];

                             lhs.ref(x, y, z) = expr.eval(x, y, z);
                          };
                       }
                       else {
                          std::ostringstream msg;
                          msg << "Nebo error in " << "Nebo Assignment" << ":\n";
                          msg << "Left-hand side of assignment allocated on ";
                          msg << "CPU but right-hand side is not ";
                          msg << "(completely) accessible on the same CPU";
                          msg << "\n";
                          msg << "\t - " << __FILE__ << " : " << __LINE__;
                          throw(std::runtime_error(msg.str()));
                       };
                    }
                    else {
                       std::ostringstream msg;
                       msg << "Nebo error in " << "Nebo Assignment" << ":\n";
                       msg << "Left-hand side of assignment allocated on ";
                       msg << "unknown device - not on CPU or GPU";
                       msg << "\n";
                       msg << "\t - " << __FILE__ << " : " << __LINE__;
                       throw(std::runtime_error(msg.str()));
                    };
                 }
              #else
                 {
                    SeqWalkType lhs = init();

                    typename RhsType::SeqWalkType expr = rhs.init(IntVec(0, 0, 0),
                                                                  GhostData(0),
                                                                  IntVec(0, 0, 0));

                    std::vector<IntVec>::const_iterator ip = mask.points().begin();

                    std::vector<IntVec>::const_iterator const ep = mask.points().end();

                    for(; ip != ep; ip++) {
                       int const x = (*ip)[0];

                       int const y = (*ip)[1];

                       int const z = (*ip)[2];

                       lhs.ref(x, y, z) = expr.eval(x, y, z);
                    };
                 }
              #endif
              /* __CUDACC__ */;

              #ifdef NEBO_REPORT_BACKEND
                 std::cout << "Finished Nebo masked assignment" << std::endl
              #endif
              /* NEBO_REPORT_BACKEND */;
           }

          inline SeqWalkType init(void) { return SeqWalkType(field_); }

         private:
          template<typename RhsType>
           inline void cpu_assign(RhsType rhs,
                                  IntVec const & extents,
                                  GhostData const & ghosts,
                                  IntVec const & hasBC,
                                  GhostData const limits) {
              #ifdef ENABLE_THREADS
                 if(is_thread_parallel()) {
                    thread_parallel_assign<RhsType>(rhs,
                                                    extents,
                                                    ghosts,
                                                    hasBC,
                                                    limits);
                 }
                 else {
                    sequential_assign<RhsType>(rhs,
                                               extents,
                                               ghosts,
                                               hasBC,
                                               limits);
                 }
              #else
                 sequential_assign<RhsType>(rhs, extents, ghosts, hasBC, limits)
              #endif
              /* ENABLE_THREADS */;
           }

          template<typename RhsType>
           inline void sequential_assign(RhsType rhs,
                                         IntVec const & extents,
                                         GhostData const & ghosts,
                                         IntVec const & hasBC,
                                         GhostData const limits) {
              #ifdef NEBO_REPORT_BACKEND
                 std::cout << "Starting Nebo sequential" << std::endl
              #endif
              /* NEBO_REPORT_BACKEND */;

              init().seqwalk_assign(rhs.init(extents, ghosts, hasBC), limits);

              #ifdef NEBO_REPORT_BACKEND
                 std::cout << "Finished Nebo sequential" << std::endl
              #endif
              /* NEBO_REPORT_BACKEND */;
           }

          #ifdef ENABLE_THREADS
             template<typename RhsType>
              inline void thread_parallel_assign(RhsType rhs,
                                                 IntVec const & extents,
                                                 GhostData const & ghosts,
                                                 IntVec const & hasBC,
                                                 GhostData const limits) {
                 #ifdef NEBO_REPORT_BACKEND
                    std::cout << "Starting Nebo thread parallel" << std::endl
                 #endif
                 /* NEBO_REPORT_BACKEND */;

                 Semaphore semaphore(0);

                 const int thread_count = field_.get_partition_count();

                 typename RhsType::ResizeType typedef RhsResizeType;

                 ResizeType new_lhs = resize();

                 RhsResizeType new_rhs = rhs.resize();

                 GhostData localLimits;

                 const IntVec split = nebo_find_partition(IntVec(limits.get_plus(0)
                                                                 - limits.get_minus(0),
                                                                 limits.get_plus(1)
                                                                 - limits.get_minus(1),
                                                                 limits.get_plus(2)
                                                                 - limits.get_minus(2)),
                                                          thread_count);

                 const int max = nebo_partition_count(split);

                 IntVec location = IntVec(0, 0, 0);

                 for(int count = 0; count < max; count++) {
                    nebo_set_up_extents(location, split, localLimits, limits);

                    ThreadPoolFIFO::self().schedule(boost::bind(&ResizeType::
                                                                template
                                                                resize_assign<RhsResizeType>,
                                                                new_lhs,
                                                                new_rhs,
                                                                extents,
                                                                ghosts,
                                                                hasBC,
                                                                localLimits,
                                                                &semaphore));

                    location = nebo_next_partition(location, split);
                 };

                 for(int ii = 0; ii < max; ii++) { semaphore.wait(); };

                 #ifdef NEBO_REPORT_BACKEND
                    std::cout << "Finished Nebo thread parallel" << std::endl
                 #endif
                 /* NEBO_REPORT_BACKEND */;
              }

             inline ResizeType resize(void) { return ResizeType(field_); }
          #endif
          /* ENABLE_THREADS */

          #ifdef __CUDACC__
             template<typename RhsType>
              inline void gpu_assign(RhsType rhs,
                                     IntVec const & extents,
                                     GhostData const & ghosts,
                                     IntVec const & hasBC,
                                     GhostData const limits) {
                 #ifdef NEBO_REPORT_BACKEND
                    std::cout << "Starting Nebo CUDA" << std::endl
                 #endif
                 /* NEBO_REPORT_BACKEND */;

                 typename RhsType::GPUWalkType typedef RhsGPUWalkType;

                 int xExtent = limits.get_plus(0) - limits.get_minus(0);

                 int yExtent = limits.get_plus(1) - limits.get_minus(1);

                 int blockDim = 16;

                 int xGDim = xExtent / blockDim + ((xExtent % blockDim) > 0 ? 1
                                                   : 0);

                 int yGDim = yExtent / blockDim + ((yExtent % blockDim) > 0 ? 1
                                                   : 0);

                 dim3 dimBlock(blockDim, blockDim);

                 dim3 dimGrid(xGDim, yGDim);

                 #ifndef NDEBUG
                    cudaError err;

                    if(cudaSuccess != (err = cudaStreamSynchronize(field_.get_stream())))
                    {
                       std::ostringstream msg;
                       msg << "Nebo error in " << "CUDA Kernel - before call" <<
                       ":\n";
                       msg << "	 - " << cudaGetErrorString(err);
                       msg << "\n";
                       msg << "\t - " << __FILE__ << " : " << __LINE__;
                       throw(std::runtime_error(msg.str()));;
                    }
                 #endif
                 /* NDEBUG */;

                 gpu_assign_kernel<GPUWalkType, RhsGPUWalkType><<<dimGrid,
                                                                  dimBlock,
                                                                  0,
                                                                  field_.get_stream()>>>(gpu_init(),
                                                                                         rhs.gpu_init(extents,
                                                                                                      ghosts,
                                                                                                      hasBC,
                                                                                                      gpu_device_index()),
                                                                                         limits.get_minus(0),
                                                                                         limits.get_plus(0),
                                                                                         limits.get_minus(1),
                                                                                         limits.get_plus(1),
                                                                                         limits.get_minus(2),
                                                                                         limits.get_plus(2));

                 #ifndef NDEBUG
                    if(cudaSuccess != (err = cudaStreamSynchronize(field_.get_stream())))
                    {
                       std::ostringstream msg;
                       msg << "Nebo error in " << "CUDA Kernel - after call" <<
                       ":\n";
                       msg << "	 - " << cudaGetErrorString(err);
                       msg << "\n";
                       msg << "\t - " << __FILE__ << " : " << __LINE__;
                       throw(std::runtime_error(msg.str()));;
                    }
                 #endif
                 /* NDEBUG */;

                 #ifdef NEBO_REPORT_BACKEND
                    std::cout << "Finished Nebo CUDA" << std::endl
                 #endif
                 /* NEBO_REPORT_BACKEND */;
              }

             inline bool cpu_ready(void) const {
                return IS_CPU_INDEX(field_.active_device_index());
             }

             inline bool gpu_ready(void) const {
                return IS_GPU_INDEX(field_.active_device_index());
             }

             inline int gpu_device_index(void) const {
                return field_.active_device_index();
             }

             inline GPUWalkType gpu_init(void) { return GPUWalkType(field_); }

             #ifdef NEBO_GPU_TEST
                template<typename RhsType>
                 inline void gpu_test_assign(RhsType rhs,
                                             IntVec const & extents,
                                             GhostData const & ghosts,
                                             IntVec const & hasBC,
                                             GhostData const limits) {
                    #ifdef NEBO_REPORT_BACKEND
                       std::cout << "Starting Nebo CUDA with Nebo copying" <<
                       std::endl
                    #endif
                    /* NEBO_REPORT_BACKEND */;

                    rhs.gpu_prep(0);

                    if(CPU_INDEX == field_.active_device_index()) {
                       FieldType gpu_field(field_.window_with_ghost(),
                                           field_.boundary_info(),
                                           field_.get_valid_ghost_data(),
                                           NULL,
                                           InternalStorage,
                                           GPU_INDEX);

                       NeboField<Initial, FieldType> gpu_lhs(gpu_field);

                       ema::cuda::CUDADeviceInterface & CDI = ema::cuda::
                       CUDADeviceInterface::self();

                       FieldType const & ftmp_ = field_;

                       CDI.memcpy_to(gpu_field.field_values(GPU_INDEX),
                                     ftmp_.field_values(),
                                     ftmp_.allocated_bytes(),
                                     0,
                                     ftmp_.get_stream());

                       gpu_lhs.template gpu_assign<RhsType>(rhs,
                                                            extents,
                                                            ghosts,
                                                            hasBC,
                                                            limits);

                       CDI.memcpy_from(field_.field_values(),
                                       gpu_field.field_values(GPU_INDEX),
                                       field_.allocated_bytes(),
                                       0,
                                       field_.get_stream());
                    }
                    else {
                       gpu_assign<RhsType>(rhs, extents, ghosts, hasBC, limits);
                    };

                    #ifdef NEBO_REPORT_BACKEND
                       std::cout << "Finished Nebo CUDA with Nebo copying" <<
                       std::endl
                    #endif
                    /* NEBO_REPORT_BACKEND */;
                 }
             #endif
             /* NEBO_GPU_TEST */
          #endif
          /* __CUDACC__ */

          FieldType field_;
      };
      #ifdef ENABLE_THREADS
         template<typename FieldType>
          struct NeboField<Resize, FieldType> {
            public:
             FieldType typedef field_type;

             NeboField<SeqWalk, FieldType> typedef SeqWalkType;

             NeboField(FieldType f)
             : field_(f)
             {}

             #ifdef ENABLE_THREADS
                template<typename RhsType>
                 inline void resize_assign(RhsType const & rhs,
                                           IntVec const & extents,
                                           GhostData const & ghosts,
                                           IntVec const & hasBC,
                                           GhostData const limits,
                                           Semaphore * semaphore) {
                    init().seqwalk_assign(rhs.init(extents, ghosts, hasBC),
                                          limits);

                    semaphore->post();
                 }
             #endif
             /* ENABLE_THREADS */

            private:
             inline SeqWalkType init(void) { return SeqWalkType(field_); }

             FieldType field_;
         }
      #endif
      /* ENABLE_THREADS */;
      template<typename FieldType>
       struct NeboField<SeqWalk, FieldType> {
         public:
          FieldType typedef field_type;

          typename field_type::value_type typedef value_type;

          NeboField(FieldType f)
          : xGlob_(f.window_with_ghost().glob_dim(0)),
            yGlob_(f.window_with_ghost().glob_dim(1)),
            base_(f.field_values(CPU_INDEX) + (f.window_with_ghost().offset(0) +
                                               f.get_valid_ghost_data().get_minus(0))
                  + (f.window_with_ghost().glob_dim(0) * ((f.window_with_ghost().offset(1)
                                                           + f.get_valid_ghost_data().get_minus(1))
                                                          + (f.window_with_ghost().glob_dim(1)
                                                             * (f.window_with_ghost().offset(2)
                                                                + f.get_valid_ghost_data().get_minus(2))))))
          {}

          template<typename RhsType>
           inline void seqwalk_assign(RhsType rhs, GhostData const limits) {
              for(int z = limits.get_minus(2); z < limits.get_plus(2); z++) {
                 for(int y = limits.get_minus(1); y < limits.get_plus(1); y++) {
                    for(int x = limits.get_minus(0); x < limits.get_plus(0); x++)
                    { ref(x, y, z) = rhs.eval(x, y, z); };
                 };
              };
           }

          inline value_type & ref(int const x, int const y, int const z) {
             return base_[x + xGlob_ * (y + (yGlob_ * z))];
          }

         private:
          int const xGlob_;

          int const yGlob_;

          value_type * base_;
      };
      #ifdef __CUDACC__
         template<typename FieldType>
          struct NeboField<GPUWalk, FieldType> {
            public:
             FieldType typedef field_type;

             typename field_type::value_type typedef value_type;

             NeboField(FieldType f)
             : base_(f.field_values(f.active_device_index()) + (f.window_with_ghost().offset(0)
                                                                + f.get_valid_ghost_data().get_minus(0))
                     + (f.window_with_ghost().glob_dim(0) * ((f.window_with_ghost().offset(1)
                                                              + f.get_valid_ghost_data().get_minus(1))
                                                             + (f.window_with_ghost().glob_dim(1)
                                                                * (f.window_with_ghost().offset(2)
                                                                   + f.get_valid_ghost_data().get_minus(2)))))),
               valid_(false),
               xGlob_(f.window_with_ghost().glob_dim(0)),
               yGlob_(f.window_with_ghost().glob_dim(1))
             {}

             template<typename RhsType>
              __device__ inline void gpuwalk_assign(RhsType rhs,
                                                    int const xLow,
                                                    int const xHigh,
                                                    int const yLow,
                                                    int const yHigh,
                                                    int const zLow,
                                                    int const zHigh) {
                 const int ii = blockIdx.x * blockDim.x + threadIdx.x;

                 const int jj = blockIdx.y * blockDim.y + threadIdx.y;

                 const int x = ii + xLow;

                 const int y = jj + yLow;

                 start(x, y, xHigh, yHigh);

                 for(int z = zLow; z < zHigh; z++) {
                    if(valid()) { ref(x, y, z) = rhs.eval(x, y, z); };
                 };
              }

            private:
             __device__ inline bool valid(void) { return valid_; }

             __device__ inline void start(int x,
                                          int y,
                                          int const xHigh,
                                          int const yHigh) {
                valid_ = (x < xHigh && y < yHigh);
             }

             __device__ inline value_type & ref(int const x,
                                                int const y,
                                                int const z) {
                return base_[x + xGlob_ * (y + (yGlob_ * z))];
             }

             value_type * base_;

             int valid_;

             int const xGlob_;

             int const yGlob_;
         }
      #endif
      /* __CUDACC__ */;
   } /* SpatialOps */

#endif
/* NEBO_LHS_H */
