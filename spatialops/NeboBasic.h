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

#ifndef NEBO_BASIC_H
#  define NEBO_BASIC_H

#  include <spatialops/SpatialOpsConfigure.h>
#  include <spatialops/structured/IndexTriplet.h>
#  include <spatialops/structured/GhostData.h>
#  include <spatialops/structured/SpatialField.h>
#  include <spatialops/structured/SpatialMask.h>
#  include <spatialops/structured/FVStaggeredFieldTypes.h>
#  include <cmath>
#  include <math.h>

#  ifdef NEBO_REPORT_BACKEND
#     include <iostream>
#  endif
   /* NEBO_REPORT_BACKEND */

#  ifdef FIELD_EXPRESSION_THREADS
#     include <spatialops/SpatialOpsTools.h>
#     include <vector>
#     include <boost/bind.hpp>
#     include <spatialops/ThreadPool.h>
#     include <spatialops/structured/IntVec.h>
#     include <spatialops/Semaphore.h>
#  endif
   /* FIELD_EXPRESSION_THREADS */

#  ifdef __CUDACC__
#     include <sstream>
#     include <spatialops/structured/MemoryTypes.h>
#  endif
   /* __CUDACC__ */

   namespace SpatialOps {
      /* Meta-programming compiler flags */
      struct All;
      struct InteriorOnly;

      inline structured::GhostData calculate_actual_ghost(bool const useGhost,
                                                          structured::GhostData const & lhs,
                                                          structured::BoundaryCellInfo const & bc,
                                                          structured::GhostData const & rhs) {
        if(bc.has_bc(0) && rhs.get_plus(0) < bc.has_extra(0)) {
          std::ostringstream msg;
          msg << "Nebo error in " << "Nebo Ghost Checking" << ":\n";
          msg << "Not enough valid extra cells to validate all interior ";
          msg << "cells in the X direction";
          msg << "\n";
          msg << "\t - " << __FILE__ << " : " << __LINE__;
          throw(std::runtime_error(msg.str()));;
        };

        if(bc.has_bc(1) && rhs.get_plus(1) < bc.has_extra(1)) {
          std::ostringstream msg;
          msg << "Nebo error in " << "Nebo Ghost Checking" << ":\n";
          msg << "Not enough valid extra cells to validate all interior ";
          msg << "cells in the Y direction";
          msg << "\n";
          msg << "\t - " << __FILE__ << " : " << __LINE__;
          throw(std::runtime_error(msg.str()));;
        };

        if(bc.has_bc(2) && rhs.get_plus(2) < bc.has_extra(2)) {
          std::ostringstream msg;
          msg << "Nebo error in " << "Nebo Ghost Checking" << ":\n";
          msg << "Not enough valid extra cells to validate all interior ";
          msg << "cells in the Z direction";
          msg << "\n";
          msg << "\t - " << __FILE__ << " : " << __LINE__;
          throw(std::runtime_error(msg.str()));;
        };

        return ((useGhost
                 ? min((lhs + point_to_ghost(bc.has_extra())), rhs)
                 : structured::GhostData(structured::IntVec(0, 0, 0),
                                         bc.has_extra()))
                - point_to_ghost(bc.has_extra()));
      };

      inline structured::GhostData calculate_limits(bool const useGhost,
                                                    structured::MemoryWindow const & lhsMemoryWindow,
                                                    structured::GhostData const & lhsCurrentGhosts,
                                                    structured::GhostData const & lhsPossibleGhosts,
                                                    structured::BoundaryCellInfo const & bc,
                                                    structured::GhostData const rhsPossibleGhosts) {
        structured::GhostData lhsActualGhosts = calculate_actual_ghost(useGhost,
                                                                       lhsPossibleGhosts,
                                                                       bc,
                                                                       rhsPossibleGhosts);
        return structured::GhostData(lhsActualGhosts.get_minus(),
                                     (lhsMemoryWindow.extent() - lhsCurrentGhosts.get_minus() - lhsCurrentGhosts.get_plus() +
                                      lhsActualGhosts.get_plus()));
                                      
      };

      template<typename FieldType>
       inline FieldType resize_ghost(FieldType const & field,
                                     structured::IntVec const & minus,
                                     structured::IntVec const & plus) {
          const structured::IntVec oldMinus = field.get_valid_ghost_data().get_minus();
          const structured::IntVec oldPlus = field.get_valid_ghost_data().get_plus();
          const structured::IntVec offsetChange = oldMinus - minus;
          const structured::IntVec extentChange = minus + plus - oldMinus - oldPlus;

          return field.reshape(extentChange, offsetChange);
       };

      template<typename Type1, typename Type2>
       struct NeboFieldCheck;

      template<typename Type>
       struct NeboFieldCheck<Type, Type> { Type typedef Result; };

      inline structured::IntVec nebo_find_partition(structured::IntVec const & extent,
                                                    int const thread_count) {
         int x = 1;
         int y = 1;
         int z = 1;

         if(thread_count <= extent[2]) { z = thread_count; }
         else if(thread_count <= extent[1]) { y = thread_count; }
         else if(thread_count <= extent[0]) { x = thread_count; };

         return structured::IntVec(x, y, z);
      };

      inline int nebo_partition_count(structured::IntVec const & split) {
         return split[0] * split[1] * split[2];
      };

      inline void nebo_set_up_extents(structured::IntVec const & current,
                                      structured::IntVec const & split,
                                      int & localXLow,
                                      int & localXHigh,
                                      int & localYLow,
                                      int & localYHigh,
                                      int & localZLow,
                                      int & localZHigh,
                                      int const xLow,
                                      int const xHigh,
                                      int const yLow,
                                      int const yHigh,
                                      int const zLow,
                                      int const zHigh) {
        using namespace structured;

        //full extent indexed from 0 rather than DLow (which is nonpositive - zero or below)
        IntVec const fullExtent(xHigh - xLow,
                                yHigh - yLow,
                                zHigh - zLow);

        //sanity checks
#       ifndef NDEBUG
          for( size_t i=0; i<3; ++i ){
            assert( fullExtent[i] >= split[i] );
            assert( split[i] > 0 );
            assert( current[i] < split[i] );
            assert( current[i] >= 0 );
          }
#       endif

        //extent of a partition
        IntVec const stdExtent = fullExtent / split;

        //number of partitions with an extra cell (to cover what is not covered by stdExtent)
        IntVec const nExtra(fullExtent[0] % split[0],
                            fullExtent[1] % split[1],
                            fullExtent[2] % split[2]);

        //number of previous paritions with an extra cell
        IntVec const pastExtra(current[0] < nExtra[0] ? current[0] : nExtra[0],
                               current[1] < nExtra[1] ? current[1] : nExtra[1],
                               current[2] < nExtra[2] ? current[2] : nExtra[2]);

        //does current partition have an extra cell
        IntVec const currentExtra(current[0] < nExtra[0] ? 1 : 0,
                                  current[1] < nExtra[1] ? 1 : 0,
                                  current[2] < nExtra[2] ? 1 : 0);

        //calculate current partition's low and high
        IntVec const low = stdExtent * current + pastExtra;
        IntVec const high = low + stdExtent + currentExtra;

        //shift back to indexing from DLow rather than zero
        localXLow = low[0] + xLow;
        localYLow = low[1] + yLow;
        localZLow = low[2] + zLow;
        localXHigh = high[0] + xLow;
        localYHigh = high[1] + yLow;
        localZHigh = high[2] + zLow;
      };

      inline structured::IntVec nebo_next_partition(structured::IntVec const & current,
                                                    structured::IntVec const & split) {
        structured::IntVec result;

        if(current[2] < split[2] - 1)
          result = structured::IntVec(current[0], current[1], 1 + current[2]);
        else if(current[1] < split[1] - 1)
          result = structured::IntVec(current[0], 1 + current[1], 0);
        else result = structured::IntVec(1 + current[0], 0, 0);

        return result;
      };

      template<typename Operand, typename FieldType>
       struct NeboExpression {
         public:
          FieldType typedef field_type;

          Operand typedef Expression;

          NeboExpression(Operand const & given)
          : expr_(given)
          {}

          inline Operand const & expr(void) const { return expr_; }

         private:
          Operand expr_;
      };

      template<typename Operand, typename T>
       struct NeboSingleValueExpression {
         public:
          SpatialOps::structured::SpatialField<SpatialOps::structured::SingleValue, T> typedef field_type;

          Operand typedef Expression;

          NeboSingleValueExpression(Operand const & given)
          : expr_(given)
          {}

          inline Operand const & expr(void) const { return expr_; }

         private:
          Operand expr_;
      };

      template<typename Operand, typename FieldType>
       struct NeboBooleanExpression {
         public:
          FieldType typedef field_type;

          Operand typedef Expression;

          NeboBooleanExpression(Operand const & given)
          : expr_(given)
          {}

          inline Operand const & expr(void) const { return expr_; }

         private:
          Operand expr_;
      };

      template<typename Operand, typename T>
       struct NeboBooleanSingleValueExpression {
         public:
          SpatialOps::structured::SpatialField<SpatialOps::structured::SingleValue, T> typedef field_type;

          Operand typedef Expression;

          NeboBooleanSingleValueExpression(Operand const & given)
          : expr_(given)
          {}

          inline Operand const & expr(void) const { return expr_; }

         private:
          Operand expr_;
      };

      /* Modes: */
      struct Initial;
#     ifdef FIELD_EXPRESSION_THREADS
         struct Resize;
#     endif
      /* FIELD_EXPRESSION_THREADS */
      struct SeqWalk;
#     ifdef __CUDACC__
         struct GPUWalk
#     endif
      /* __CUDACC__ */
      struct Reduction;
   } /* SpatialOps */

#endif
/* NEBO_BASIC_H */
