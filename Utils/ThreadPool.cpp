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
// https://github.com/google/or-tools/blob/v7.4/ortools/base/threadpool.cc

#include "ThreadPool.h"

namespace matrix::matrix_utils {

namespace {

void RunWorker(void* data) {
  auto* const thread_pool = reinterpret_cast<ThreadPool*>(data);
  std::function<void()> work = thread_pool->GetNextTask();
  while (work != nullptr) {
    work();
    work = thread_pool->GetNextTask();
  }
}

}  // namespace

ThreadPool::ThreadPool(int num_workers) : num_workers_(num_workers) {}

ThreadPool::~ThreadPool() {
  if (!started_) return;
  {
    std::lock_guard<std::mutex> mutex_lock(mutex_);
    waiting_to_finish_ = true;
    condition_.notify_all();
  }
  for (int i = 0; i < num_workers_; ++i) {
    all_workers_[i].join();
  }
}

void ThreadPool::SetQueueCapacity(int capacity) {
  queue_capacity_ = capacity;
}

void ThreadPool::StartWorkers() {
  started_ = true;
  for (int i = 0; i < num_workers_; ++i) {
    all_workers_.emplace_back(&RunWorker, this);
  }
}

std::function<void()> ThreadPool::GetNextTask() {
  std::unique_lock<std::mutex> lock(mutex_);
  while (true) {
    if (!tasks_.empty()) {
      std::function<void()> task = tasks_.front();
      tasks_.pop_front();
      if (tasks_.size() < queue_capacity_ && waiting_for_capacity_) {
        waiting_for_capacity_ = false;
        capacity_condition_.notify_all();
      }
      return task;
    }
    if (waiting_to_finish_) {
      return nullptr;
    }
    condition_.wait(lock);
  }
  return nullptr;
}

void ThreadPool::Schedule(const std::function<void()>& closure) {
  std::unique_lock<std::mutex> lock(mutex_);
  while (tasks_.size() >= queue_capacity_) {
    waiting_for_capacity_ = true;
    capacity_condition_.wait(lock);
  }
  tasks_.push_back(closure);
  if (started_) {
    condition_.notify_all();
  }
}

}  // namespace matrix::matrix_utils
