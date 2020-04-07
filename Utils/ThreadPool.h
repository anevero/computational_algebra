// Copyright 2010-2018 Google LLC
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Source:
// https://github.com/google/or-tools/blob/v7.4/ortools/base/threadpool.h

#ifndef UTILS_THREADPOOL_H_
#define UTILS_THREADPOOL_H_

#include <condition_variable>
#include <functional>
#include <list>
#include <mutex>
#include <thread>
#include <vector>

namespace matrix::matrix_utils {

class ThreadPool {
 public:
  explicit ThreadPool(int num_threads);
  ~ThreadPool();

  void StartWorkers();
  void Schedule(const std::function<void()>& closure);
  std::function<void()> GetNextTask();
  void SetQueueCapacity(int capacity);

 private:
  const int num_workers_;

  std::list<std::function<void()>> tasks_;
  int queue_capacity_ = 2'000'000'000;

  std::mutex mutex_;
  std::condition_variable condition_;
  std::condition_variable capacity_condition_;

  bool waiting_to_finish_ = false;
  bool waiting_for_capacity_ = false;
  bool started_ = false;

  std::vector<std::thread> all_workers_;
};

}  // namespace matrix::matrix_utils

#endif  // UTILS_THREADPOOL_H_
