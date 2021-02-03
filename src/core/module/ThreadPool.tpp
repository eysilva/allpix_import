/**
 * @file
 * @brief Template implementation of thread pool for module multithreading
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <cassert>
#include <climits>

namespace allpix {
    template <typename T> ThreadPool::SafeQueue<T>::SafeQueue(unsigned int max_size) : max_size_(max_size) {}

    template <typename T> ThreadPool::SafeQueue<T>::~SafeQueue() { invalidate(); }

    /*
     * Block until a value is available if the wait parameter is set to true. The wait exits when the queue is invalidated.
     */
    template <typename T> bool ThreadPool::SafeQueue<T>::pop(T& out, const std::function<void()>& func) {
        // Lock the mutex
        std::unique_lock<std::mutex> lock{mutex_};
        if(!valid_) {
            return false;
        }

        unsigned int buffer = 0; // TODO: add as param
        assert(buffer <= max_size_);

        // Wait for one of the queues to be available
        bool pop_priority = !priority_queue_.empty() && priority_queue_.top().first == current_id_;
        bool pop_standard = !queue_.empty() && priority_queue_.size() + buffer < max_size_;
        while(!pop_priority && !pop_standard) {
            // Wait for new item in the queue (unlocks the mutex while waiting)
            pop_condition_.wait(lock);
            if(!valid_) {
                return false;
            }
            pop_priority = !priority_queue_.empty() && priority_queue_.top().first == current_id_;
            pop_standard = !queue_.empty() && priority_queue_.size() + buffer < max_size_;
        }

        // Pop the appropriate queue
        if(pop_priority) {
            // Priority queue is missing a pop returning a non-const reference, so need to apply a const_cast
            out = std::move(const_cast<PQValue&>(priority_queue_.top())).second; // NOLINT
            priority_queue_.pop();
            pop_condition_.notify_one();
        } else { // pop_standard
            out = std::move(queue_.front());
            queue_.pop();
        }

        // Optionally execute the mutex protected function
        if(func != nullptr) {
            func();
        }

        // Notify possible pusher waiting to fill the queue
        push_condition_.notify_one();
        return true;
    }

    template <typename T> bool ThreadPool::SafeQueue<T>::push(T value) {
        bool wait = false;
        // Lock the mutex
        std::unique_lock<std::mutex> lock{mutex_};

        // Check if the queue reached its full size
        if(queue_.size() >= max_size_) {
            // Wait until the queue is below the max size or it was invalidated(shutdown)
            wait = true;
            push_condition_.wait(lock, [this]() { return queue_.size() < max_size_ || !valid_; });
        }

        // Abort the push operation if conditions not met
        if(queue_.size() >= max_size_ || !valid_) {
            return wait;
        }

        // Push a new element to the queue and notify possible consumer
        queue_.push(std::move(value));
        pop_condition_.notify_one();
        return wait;
    }

    template <typename T> bool ThreadPool::SafeQueue<T>::push(uint64_t n, T value) {
        bool wait = false;
        // Lock the mutex
        std::unique_lock<std::mutex> lock{mutex_};
        assert(n >= current_id_);

        // Check if the queue reached its full size
        if(priority_queue_.size() >= max_size_) {
            // Wait until the queue is below the max size or it was invalidated(shutdown)
            wait = true;
            push_condition_.wait(lock, [this]() { return priority_queue_.size() < max_size_ || !valid_; });
        }

        // Abort the push operation if conditions not met
        if(priority_queue_.size() >= max_size_ || !valid_) {
            return wait;
        }

        // Push a new element to the queue and notify possible consumer
        priority_queue_.push(std::make_pair(n, std::move(value)));
        pop_condition_.notify_one();
        return wait;
    }

    template <typename T> void ThreadPool::SafeQueue<T>::complete(uint64_t n) {
        std::lock_guard<std::mutex> lock{mutex_};
        completed_ids_.insert(n);
        auto iter = completed_ids_.begin();
        // TODO reason about performance impact of this loop
        while(iter != completed_ids_.end()) {
            if(*iter != current_id_) {
                return;
            }
            completed_ids_.erase(iter);
            ++current_id_;
            pop_condition_.notify_one();
            ++iter;
        }
    }

    template <typename T> uint64_t ThreadPool::SafeQueue<T>::currentId() const {
        std::lock_guard<std::mutex> lock{mutex_};
        return current_id_;
    }

    template <typename T> bool ThreadPool::SafeQueue<T>::isValid() const {
        std::lock_guard<std::mutex> lock{mutex_};
        return valid_;
    }

    template <typename T> bool ThreadPool::SafeQueue<T>::empty() const {
        std::lock_guard<std::mutex> lock{mutex_};
        return !valid_ || (queue_.empty() && priority_queue_.empty());
    }

    template <typename T> size_t ThreadPool::SafeQueue<T>::size() const {
        std::lock_guard<std::mutex> lock{mutex_};
        return queue_.size() + priority_queue_.size();
    }

    template <typename T> size_t ThreadPool::SafeQueue<T>::prioritySize() const {
        std::lock_guard<std::mutex> lock{mutex_};
        return priority_queue_.size();
    }

    /*
     * Used to ensure no conditions are being waited for in pop when a thread or the application is trying to exit. The queue
     * is invalid after calling this method and it is an error to continue using a queue after this method has been called.
     */
    template <typename T> void ThreadPool::SafeQueue<T>::invalidate() {
        std::lock_guard<std::mutex> lock{mutex_};
        std::priority_queue<PQValue, std::vector<PQValue>, std::greater<PQValue>>().swap(priority_queue_);
        std::queue<T>().swap(queue_);
        valid_ = false;
        push_condition_.notify_all();
        pop_condition_.notify_all();
    }

    template <typename Func, typename... Args> auto ThreadPool::submit(Func&& func, Args&&... args) {
        return submit(UINT64_MAX, std::forward<Func>(func), std::forward<Args>(args)...);
    }

    template <typename Func, typename... Args> auto ThreadPool::submit(uint64_t n, Func&& func, Args&&... args) {
        // Bind the arguments to the tasks
        auto bound_task = std::bind(std::forward<Func>(func), std::forward<Args>(args)...);

        // Construct packaged task with correct return type
        using PackagedTask = std::packaged_task<decltype(bound_task())()>;
        PackagedTask task(bound_task);

        // Get future and wrapper to add to vector
        auto future = task.get_future().share();
        auto task_function = [task = std::move(task), future = future]() mutable {
            task();
            future.get();
        };
        if(threads_.empty()) {
            task_function();
        } else {
            if(n == UINT64_MAX) {
                queue_.push(std::make_unique<std::packaged_task<void()>>(std::move(task_function)));
            } else {
                queue_.push(n, std::make_unique<std::packaged_task<void()>>(std::move(task_function)));
            }
        }
        return future;
    }

} // namespace allpix
