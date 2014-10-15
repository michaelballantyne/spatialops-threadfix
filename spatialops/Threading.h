#ifndef NEBO_THREADING_H
#define NEBO_THREADING_H

#include <boost/atomic.hpp>
#include <boost/function.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <vector>
#include <stdexcept>
#include <string>
#include <sched.h>

namespace SpatialOps {
    class Threading {

        public:
            static void setup(int nthreads) {
                instance = new Threading(nthreads);
            }

            static Threading & getInstance() {
                return *instance;
            }

            void scheduleTasks(std::vector<boost::function0<void> > & tasks) {
                //if (tasks.size() != nthreads) {
                    //throw std::invalid_argument("ntasks must match nthreads");
                //}

                for (int i = 0; i < nthreads; i++) {
                    tasksArray[i]->store(&tasks[i], boost::memory_order_relaxed);
                }

                remaining.store(nthreads, boost::memory_order_release);

                while(remaining.load(boost::memory_order_acquire) > 0) {
                    // Spin until done
                }
            }

            const int nthreads;

        private:
            static Threading *instance;

            // Make sure the only way to construct one is via getInstance
            Threading(const Threading & other);

            boost::atomic<boost::function0<void> *> **tasksArray;
            boost::thread **threadsArray;
            boost::atomic<int> remaining;

            // (Private) constructor taking thread count.
            Threading(int nthreadsArg) : nthreads(nthreadsArg), remaining(0) {
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(0, &cpuset);
            //    sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);

                // Create atomics for communication with threads
                tasksArray = new boost::atomic<boost::function0<void> *> *[nthreads];
                for (int i = 0; i < nthreads; i++) {
                    tasksArray[i] = new boost::atomic<boost::function0<void> *>(NULL);
                }

                // Start threads
                threadsArray = new boost::thread *[nthreads];
                for (int i = 0; i < nthreads; i++) {
                    threadsArray[i] = new boost::thread(boost::bind(&Threading::threadBody, this, i));
                }
            }

            // I probably don't need this whole thing at all.
            ~Threading() {
                // Interrupt and threads and wait for them to die.
                for (int i = 0; i < nthreads; i++) {
                    threadsArray[i]->interrupt();
                }

                for (int i = 0; i < nthreads; i++) {
                    threadsArray[i]->join();
                    delete threadsArray[i];
                }

                delete threadsArray;

                // Then we delete the boost::atomics we allocated.
                for (int i = 0; i < nthreads; i++) {
                    delete tasksArray[i];
                }

                delete tasksArray;
            }

            void threadBody(int threadId) {
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET((threadId + 1) * 2, &cpuset);
                //sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);

//                pthread_setname_np(std::to_string(threadId).c_str());
                while (true) {
                    boost::function0<void> *task = NULL;
                    while((task = tasksArray[threadId]->load(boost::memory_order_relaxed)) == NULL) {
                        // Wait for work
                    }

                    boost::atomic_thread_fence(boost::memory_order_acquire);

                    // Execute task
                    (*task)();

                    // Report back
                    tasksArray[threadId]->store(NULL, boost::memory_order_relaxed);

                    remaining.fetch_sub(1, boost::memory_order_release);

                    //boost::this_thread::interruption_point();
                }
            }
    };
}

#endif

