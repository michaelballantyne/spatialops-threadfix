#ifndef Nebo_Semaphore_h
#define Nebo_Semaphore_h

#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>

/**
 * Implemented to replace boost::interprocess::interprocess_semaphore due to posix semaphores bug in glibc:
 * http://sourceware.org/bugzilla/show_bug.cgi?id=12674
 *
 * However, there are many more bugs in glibc! On some (non x86 / x86_64) platforms, this one could bite us:
 * http://sourceware.org/bugzilla/show_bug.cgi?id=13690
 */
namespace SpatialOps{
    struct Semaphore {
        public:
            Semaphore(int initial) {
                val = initial;
            }

            inline void post() {
                boost::lock_guard<boost::mutex> lock(mut);
                val++;
                cond.notify_one();
            }

            inline void wait() {
                boost::unique_lock<boost::mutex> lock(mut);
                while (val <= 0) {
                    cond.wait(lock);
                }
                val--;
            }

        private:
            boost::condition_variable cond;
            boost::mutex mut;
            int val;
    };
}
#endif // Nebo_Semaphore_h
