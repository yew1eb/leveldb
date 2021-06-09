// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#ifndef STORAGE_LEVELDB_UTIL_ARENA_H_
#define STORAGE_LEVELDB_UTIL_ARENA_H_

#include <cstddef>
#include <vector>
#include <assert.h>
#include <stdint.h>

namespace leveldb {
/// Arena每次按kBlockSize(4096)单位向系统申请内存，提供地址对齐的内存，记录内存使用。
/// 当memtable申请内存时，如果size不大于kBlockSize的四分之一，就在当前空闲的内存block中分配，
/// 否则，直接向系统申请（malloc）。这个策略是为了能更好的服务小内存的申请，避免个别大内存使用影响。
class Arena {
 public:
  Arena();
  ~Arena();

  // Return a pointer to a newly allocated memory block of "bytes" bytes.
  char* Allocate(size_t bytes);

  // Allocate memory with the normal alignment guarantees provided by malloc
  char* AllocateAligned(size_t bytes);

  // Returns an estimate of the total memory usage of data allocated
  // by the arena (including space allocated but not yet used for user
  // allocations).
  size_t MemoryUsage() const {
    return blocks_memory_ + blocks_.capacity() * sizeof(char*);
  }

 private:
  char* AllocateFallback(size_t bytes);
  char* AllocateNewBlock(size_t block_bytes);

  // Allocation state
  char* alloc_ptr_; /// 当前空闲内存block内的可用地址
  size_t alloc_bytes_remaining_; /// 当前空闲内存block内的可用大小

  // Array of new[] allocated memory blocks
  std::vector<char*> blocks_; /// 已经申请的内存block

  // Bytes of memory in blocks allocated so far
  /// 累计分配的内存大小 // 一个memtable对应一个Arena, // memtable内的数据量就用这个值表示
  size_t blocks_memory_;

  // No copying allowed
  Arena(const Arena&);
  void operator=(const Arena&);
};

inline char* Arena::Allocate(size_t bytes) {
  // The semantics of what to return are a bit messy if we allow
  // 0-byte allocations, so we disallow them here (we don't need
  // them for our internal use).
  assert(bytes > 0);
  if (bytes <= alloc_bytes_remaining_) {
    char* result = alloc_ptr_;
    alloc_ptr_ += bytes;
    alloc_bytes_remaining_ -= bytes;
    return result;
  }
  return AllocateFallback(bytes);
}

}  // namespace leveldb

#endif  // STORAGE_LEVELDB_UTIL_ARENA_H_
