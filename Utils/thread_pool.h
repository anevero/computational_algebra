// MIT License
//
// Copyright (c) 2017 Robert Vaser
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Source: https://github.com/rvaser/thread_pool

#ifndef UTILS_THREAD_POOL_H_
#define UTILS_THREAD_POOL_H_

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <queue>
#include <utility>
#include <vector>

class Semaphore;
std::unique_ptr<Semaphore> createSemaphore(uint32_t value);

class ThreadPool;
std::unique_ptr<ThreadPool> createThreadPool(
    uint32_t num_threads = std::thread::hardware_concurrency() / 2);

class Semaphore {
 public:
  ~Semaphore() = default;
  Semaphore(const Semaphore&) = delete;
  const Semaphore& operator=(const Semaphore&) = delete;

  uint32_t value() const {
    return value_;
  }

  void wait();
  void post();

  friend std::unique_ptr<Semaphore> createSemaphore(uint32_t value);

 private:
  explicit Semaphore(uint32_t value);

  std::mutex mutex_;
  std::condition_variable condition_;
  uint32_t value_;
};

class ThreadPool {
 public:
  ~ThreadPool();
  ThreadPool(const ThreadPool&) = delete;
  const ThreadPool& operator=(const ThreadPool&) = delete;

  uint32_t num_threads() const {
    return threads_.size();
  }

  const std::vector<std::thread::id>& thread_identifiers() const {
    return thread_identifiers_;
  }

  template<typename T, typename... Ts>
  auto submit(T&& routine, Ts&& ... params)
  -> std::future<typename std::result_of<T(Ts...)>::type> {
    auto task =
        std::make_shared<
            std::packaged_task<typename std::result_of<T(Ts...)>::type()>>(
            std::bind(std::forward<T>(routine), std::forward<Ts>(params)...));
    auto task_result = task->get_future();
    auto task_wrapper = [task]() {
      (*task)();
    };

    queue_sem_->wait();

    task_queue_.emplace(task_wrapper);

    queue_sem_->post();
    active_sem_->post();

    return task_result;
  }

  friend std::unique_ptr<ThreadPool> createThreadPool(
      uint32_t num_threads);

 private:
  explicit ThreadPool(uint32_t num_threads);

  static void worker_thread(ThreadPool* thread_pool);

  std::vector<std::thread> threads_;
  std::vector<std::thread::id> thread_identifiers_;

  std::queue<std::function<void()>> task_queue_;

  std::unique_ptr<Semaphore> queue_sem_;
  std::unique_ptr<Semaphore> active_sem_;

  std::atomic<bool> terminate_;
};

#endif  // UTILS_THREAD_POOL_H_
