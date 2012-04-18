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

#include <spatialops/FieldExpressions.h>

namespace SpatialOps {

#ifdef FIELD_EXPRESSION_THREADS
    bool is_nebo_thread_parallel() { return ThreadPoolFIFO::self().is_nebo_thread_parallel(); };

    int get_nebo_soft_thread_count() { return ThreadPoolResourceManager::self().get_max_active_worker_count(ThreadPoolFIFO::self()); };
    int set_nebo_soft_thread_count(int thread_count) {
        ThreadPoolFIFO::self().set_nebo_thread_parallel(thread_count > 0);  //if thread_count is zero or less, use sequential nebo
        return ThreadPoolResourceManager::self().resize_active(ThreadPoolFIFO::self(), thread_count);
    };

    int get_nebo_hard_thread_count() { return ThreadPoolResourceManager::self().get_worker_count(ThreadPoolFIFO::self()); };
    int set_nebo_hard_thread_count(int thread_count) {
        ThreadPoolFIFO::self().set_nebo_thread_parallel(thread_count > 0);  //if thread_count is zero or less, use sequential nebo
        return ThreadPoolResourceManager::self().resize(ThreadPoolFIFO::self(), thread_count);
    };
#endif
} /* SpatialOps */

