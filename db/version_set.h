// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.
//
// The representation of a DBImpl consists of a set of Versions.  The
// newest version is called "current".  Older versions may be kept
// around to provide a consistent view to live iterators.
//
// Each Version keeps track of a set of Table files per level.  The
// entire set of versions is maintained in a VersionSet.
//
// Version,VersionSet are thread-compatible, but require external
// synchronization on all accesses.
//// 为了均衡读写的效率，sstable文件分层次（level）管理，db预定义了最大的level 值。
/// compact 进程负责level之间的均衡。
#ifndef STORAGE_LEVELDB_DB_VERSION_SET_H_
#define STORAGE_LEVELDB_DB_VERSION_SET_H_

#include <map>
#include <set>
#include <vector>
#include "db/dbformat.h"
#include "db/version_edit.h"
#include "port/port.h"

namespace leveldb {

namespace log { class Writer; }

class Compaction;
class Iterator;
class MemTable;
class TableBuilder;
class TableCache;
class Version;
class VersionSet;
class WritableFile;

// Return the smallest index i such that files[i]->largest >= key.
// Return files.size() if there is no such file.
// REQUIRES: "files" contains a sorted list of non-overlapping files.
extern int FindFile(const InternalKeyComparator& icmp,
                    const std::vector<FileMetaData*>& files,
                    const Slice& key);

// Returns true iff some file in "files" overlaps the user key range
// [*smallest,*largest].
// smallest==NULL represents a key smaller than all keys in the DB.
// largest==NULL represents a key largest than all keys in the DB.
// REQUIRES: If disjoint_sorted_files, files[] contains disjoint ranges
//           in sorted order.
extern bool SomeFileOverlapsRange(
    const InternalKeyComparator& icmp,
    bool disjoint_sorted_files,
    const std::vector<FileMetaData*>& files,
    const Slice* smallest_user_key,
    const Slice* largest_user_key);

/// 将每次compact后的最新数据状态定义为Version，也就是当前db元信息以及每个level上具有最新数据状态的sstable集合。
/// compact会在某个level上新加入或者删除一些sstable，但可能这个时候，那些要删除的sstable正在被读，为了处理这样的读写竞争情况，
/// 基于sstable文件一旦生成就不会改动的特点，每个Version加入引用计数，读以及解除读操作会将引用计数相应加减一。
/// 这样， db中可能有多个Version同时存在（提供服务），它们通过链表链接起来。
/// 当Version的引用计数为0并且不是当前最新的Version时，它会从链表中移除，对应的，
/// 该Version内的sstable就可以删除了（这些废弃的sstable会在下一次compact完成时被清理掉）。
class Version {
 public:
  // Append to *iters a sequence of iterators that will
  // yield the contents of this Version when merged together.
  // REQUIRES: This version has been saved (see VersionSet::SaveTo)
  void AddIterators(const ReadOptions&, std::vector<Iterator*>* iters);

  // Lookup the value for key.  If found, store it in *val and
  // return OK.  Else return a non-OK status.  Fills *stats.
  // REQUIRES: lock is not held
  struct GetStats {
    FileMetaData* seek_file;
    int seek_file_level;
  };
  Status Get(const ReadOptions&, const LookupKey& key, std::string* val,
             GetStats* stats);

  // Adds "stats" into the current state.  Returns true if a new
  // compaction may need to be triggered, false otherwise.
  // REQUIRES: lock is held
  bool UpdateStats(const GetStats& stats);

  // Reference count management (so Versions do not disappear out from
  // under live iterators)
  void Ref();
  void Unref();

  void GetOverlappingInputs(
      int level,
      const InternalKey* begin,         // NULL means before all keys
      const InternalKey* end,           // NULL means after all keys
      std::vector<FileMetaData*>* inputs);

  // Returns true iff some file in the specified level overlaps
  // some part of [*smallest_user_key,*largest_user_key].
  // smallest_user_key==NULL represents a key smaller than all keys in the DB.
  // largest_user_key==NULL represents a key largest than all keys in the DB.
  bool OverlapInLevel(int level,
                      const Slice* smallest_user_key,
                      const Slice* largest_user_key);

  // Return the level at which we should place a new memtable compaction
  // result that covers the range [smallest_user_key,largest_user_key].
  int PickLevelForMemTableOutput(const Slice& smallest_user_key,
                                 const Slice& largest_user_key);

  int NumFiles(int level) const { return files_[level].size(); }

  // Return a human readable string that describes this version's contents.
  std::string DebugString() const;

 private:
  friend class Compaction;
  friend class VersionSet;

  class LevelFileNumIterator;
  Iterator* NewConcatenatingIterator(const ReadOptions&, int level) const;
/// 属于的VersionSet
  VersionSet* vset_;            // VersionSet to which this Version belongs
/// 链表指针
  Version* next_;               // Next version in linked list
  Version* prev_;               // Previous version in linked list
/// 引用计数
  int refs_;                    // Number of live refs to this version

/// 每个level的所有sstable元信息。
/// files_[i]中的FileMetaData按照FileMetaData::smallest排序，
/// 这是在每次更新都保证的。（参见VersionSet::Builder::Save()）
  // List of files per level
  std::vector<FileMetaData*> files_[config::kNumLevels];

  //// 需要compact的文件（allowed_seeks用光）
  // Next file to compact based on seek stats.
  /// leveldb对单个sstable文件的IO也做了细化的优化，设计了一个巧妙的策略。
  /// 首先，一次查找如果对多于一个sstable进行了查找（对sstable进行了查找可以认为对其中的datablock进行了一次寻道seek），
  /// 说明处于低level上的sstable并没有提供高的hit比率，可以认为它处在不最优的情况，而我们认为compact后会倾向于均衡的状态，
  /// 所以在一个sstable的seek次数达到一定阈值后，主动对其进行compact是合理的。 这个具体seek次数阈值(allowed_seeks)的确定，
  /// 依赖于sas盘的IO性能：
  /// a. 一次磁盘寻道seek耗费10ms。
  /// b. 读或者写1M数据耗费10ms （按100M/s IO吞吐能力）。
  /// c. compact 1M的数据需要25M的IO：从level-n中读1M数据，从level-n+1中读10～12M数据，写入level-n+1中10～12M数据。
  /// 所以，compact 1M的数据的时间相当于做25次磁盘seek，反过来说就是，1次seek相当于compact 40k数据。
  /// 那么，可以得到seek阈值allowed_seeks=sstable_size / 40k。保守设置，当前实际的allowed_seeks = sstable_size / 16k。
  /// 每次compact完成，构造新的Version时（Builder::Apply()）,每个sstable的allowed_seeks会计算出来保存在FileMetaData。
  /// 在每次get操作的时候，如果有超过一个sstable文件进行了查找，会将第一个进行查找的sstable的allowed_seeks减一，
  /// 并检查其是否已经用光了allowed_seeks,若是，则将该sstable记录成当前Version的file_to_compact_,
  /// 并记录其所在的level(file_to_compact_level_)。
  FileMetaData* file_to_compact_;
  int file_to_compact_level_;

  // Level that should be compacted next and its compaction score.
  // Score < 1 means compaction is not strictly needed.  These fields
  // are initialized by Finalize().
  //// 需要compact的文件（allowed_seeks用光）
  /// leveldb中分level管理sstable，对于写，可以认为与sstable无关。
  /// 而基于get的流程（参见get流程），各level中的sstable的count，size以及range分布，会直接影响读的效率。
  /// 可以预想的最佳情形可能是level-0中最多有一个sstable，level-1以及之上的各level中key-range分布均匀，
  /// 期望更多的查找可以遍历最少的level即可定位到。 将这种预想的最佳状态定义成: level处于均衡的状态。
  /// 当采用具体的参数量化，也就量化了各个level的不均衡比重，即compact权重： score。score越大，表示该level越不均衡，
  /// 需要更优先进行compact。 每个level的具体均衡参数及比重计算策略如下：
  /// a. 因为level-0的sstable range可能overlap，所以如果level-0上有过多的sstable，在做查找时，会严重影响效率。
  /// 同时，因为level-0中的sstable由memtable直接dump得到，并不受kTargetFileSize（生成sstable的size）的控制，
  /// 所以sstable的count更有意义。基于此，对于level-0，均衡的状态需要满足：sstable 的count < kL0_CompactionTrigger。
  /// score = sstable的 count/ kL0_CompactionTrigger。 为了控制这个数量，
  /// 另外还有kL0_SlowdownWritesTrigger/kL0_StopWritesTrigger两个阈值来主动控制写的速率（参见put流程）。
  /// b. 对于level-1及以上的level，sstable均由compact过程产生，生成的sstable大小被kTargetFileSize控制，所以可以限定sstable总的size。
  /// 当前的策略是设置初始值kBaseLevelSize，然后以10的指数级按level增长。
  /// 每个level可以容纳的quota_size = kBaseLevelSize * 10^(level_number-1)。
  /// 所以level-1可以容纳总共kBaseLevelSize的sstable，level-2 允许kBaseLevelSize*10……
  /// 基于此，对于level-1及以上的level 均衡的状态需要满足：sstable的size < quota_size。 score = sstable的size / quota_size。
  /// 每次compact完成，生效新的Version时（VersionSet::Finalize()），都会根据上述的策略，
  /// 计算出每个level的score,取最大值作为当前Version的compaction_score_,同时记录对应的level(compaction_level_)。
  double compaction_score_;
  int compaction_level_;

  explicit Version(VersionSet* vset)
      : vset_(vset), next_(this), prev_(this), refs_(0),
        file_to_compact_(NULL),
        file_to_compact_level_(-1),
        compaction_score_(-1),
        compaction_level_(-1) {
  }

  ~Version();

  // No copying allowed
  Version(const Version&);
  void operator=(const Version&);
};

/// 整个db的当前状态被VersionSet管理着，其中有当前最新的Version以及其他正在服务的Version链表；
/// 全局的SequnceNumber，FileNumber；当前的manifest_file_number; 封装sstable的TableCache。
/// 每个level中下一次compact要选取的start_key等等。
class VersionSet {
 public:
  VersionSet(const std::string& dbname,
             const Options* options,
             TableCache* table_cache,
             const InternalKeyComparator*);
  ~VersionSet();

  // Apply *edit to the current version to form a new descriptor that
  // is both saved to persistent state and installed as the new
  // current version.  Will release *mu while actually writing to the file.
  // REQUIRES: *mu is held on entry.
  // REQUIRES: no other thread concurrently calls LogAndApply()
  Status LogAndApply(VersionEdit* edit, port::Mutex* mu);

  // Recover the last saved descriptor from persistent storage.
  Status Recover();

  // Return the current version.
  Version* current() const { return current_; }

  // Return the current manifest file number
  uint64_t ManifestFileNumber() const { return manifest_file_number_; }

  // Allocate and return a new file number
  uint64_t NewFileNumber() { return next_file_number_++; }

  // Return the number of Table files at the specified level.
  int NumLevelFiles(int level) const;

  // Return the combined file size of all files at the specified level.
  int64_t NumLevelBytes(int level) const;

  // Return the last sequence number.
  uint64_t LastSequence() const { return last_sequence_; }

  // Set the last sequence number to s.
  void SetLastSequence(uint64_t s) {
    assert(s >= last_sequence_);
    last_sequence_ = s;
  }

  // Mark the specified file number as used.
  void MarkFileNumberUsed(uint64_t number);

  // Return the current log file number.
  uint64_t LogNumber() const { return log_number_; }

  // Return the log file number for the log file that is currently
  // being compacted, or zero if there is no such log file.
  uint64_t PrevLogNumber() const { return prev_log_number_; }

  // Pick level and inputs for a new compaction.
  // Returns NULL if there is no compaction to be done.
  // Otherwise returns a pointer to a heap-allocated object that
  // describes the compaction.  Caller should delete the result.
  Compaction* PickCompaction();

  // Return a compaction object for compacting the range [begin,end] in
  // the specified level.  Returns NULL if there is nothing in that
  // level that overlaps the specified range.  Caller should delete
  // the result.
  Compaction* CompactRange(
      int level,
      const InternalKey* begin,
      const InternalKey* end);

  // Return the maximum overlapping data (in bytes) at next level for any
  // file at a level >= 1.
  int64_t MaxNextLevelOverlappingBytes();

  // Create an iterator that reads over the compaction inputs for "*c".
  // The caller should delete the iterator when no longer needed.
  Iterator* MakeInputIterator(Compaction* c);

  // Returns true iff some level needs a compaction.
  bool NeedsCompaction() const {
    Version* v = current_;
    return (v->compaction_score_ >= 1) || (v->file_to_compact_ != NULL);
  }

  // Add all files listed in any live version to *live.
  // May also mutate some internal state.
  void AddLiveFiles(std::set<uint64_t>* live);

  // Return the approximate offset in the database of the data for
  // "key" as of version "v".
  uint64_t ApproximateOffsetOf(Version* v, const InternalKey& key);

  // Return a human-readable short (single-line) summary of the number
  // of files per level.  Uses *scratch as backing store.
  struct LevelSummaryStorage {
    char buffer[100];
  };
  const char* LevelSummary(LevelSummaryStorage* scratch) const;

 private:
  class Builder;

  friend class Compaction;
  friend class Version;

  void Finalize(Version* v);

  void GetRange(const std::vector<FileMetaData*>& inputs,
                InternalKey* smallest,
                InternalKey* largest);

  void GetRange2(const std::vector<FileMetaData*>& inputs1,
                 const std::vector<FileMetaData*>& inputs2,
                 InternalKey* smallest,
                 InternalKey* largest);

  void SetupOtherInputs(Compaction* c);

  // Save current contents to *log
  Status WriteSnapshot(log::Writer* log);

  void AppendVersion(Version* v);
/// 实际的Env
  Env* const env_;
  const std::string dbname_; /// db的数据路径
  const Options* const options_; //// 传入的option
  TableCache* const table_cache_; /// 操作sstable的TableCache
  const InternalKeyComparator icmp_; /// comparator
  uint64_t next_file_number_; /// 下一个可用的FileNumber
  uint64_t manifest_file_number_; /// manifest文件的FileNumber
  uint64_t last_sequence_; /// 最后用过的SequnceNumber
  uint64_t log_number_; /// log文件的FileNumber
  /// 辅助log文件的FileNumber，在compact memtable时，置为0.
  uint64_t prev_log_number_;  // 0 or backing store for memtable being compacted

  // Opened lazily
  WritableFile* descriptor_file_; /// manifest文件的封装
  log::Writer* descriptor_log_; /// manifest文件的writer
  /// 正在服务的Version链表
  Version dummy_versions_;  // Head of circular doubly-linked list of versions.
  /// 当前最新的的Version
  Version* current_;        // == dummy_versions_.prev_

  /// 为了尽量均匀compact每个level，所以会将这一次compact的end-key作为
  /// 下一次compact的start-key。compactor_pointer_就保存着每个level
  /// 下一次compact的start-key.
  /// 除了current_外的Version，并不会做compact，所以这个值并不保存在Version中。
  // Per-level key at which the next compaction at that level should start.
  // Either an empty string, or a valid InternalKey.
  std::string compact_pointer_[config::kNumLevels];

  // No copying allowed
  VersionSet(const VersionSet&);
  void operator=(const VersionSet&);
};

// A Compaction encapsulates information about a compaction.
class Compaction {
 public:
  ~Compaction();

  // Return the level that is being compacted.  Inputs from "level"
  // and "level+1" will be merged to produce a set of "level+1" files.
  int level() const { return level_; }

  // Return the object that holds the edits to the descriptor done
  // by this compaction.
  VersionEdit* edit() { return &edit_; }

  // "which" must be either 0 or 1
  int num_input_files(int which) const { return inputs_[which].size(); }

  // Return the ith input file at "level()+which" ("which" must be 0 or 1).
  FileMetaData* input(int which, int i) const { return inputs_[which][i]; }

  // Maximum size of files to build during this compaction.
  uint64_t MaxOutputFileSize() const { return max_output_file_size_; }

  // Is this a trivial compaction that can be implemented by just
  // moving a single input file to the next level (no merging or splitting)
  bool IsTrivialMove() const;

  // Add all inputs to this compaction as delete operations to *edit.
  void AddInputDeletions(VersionEdit* edit);

  // Returns true if the information we have available guarantees that
  // the compaction is producing data in "level+1" for which no data exists
  // in levels greater than "level+1".
  bool IsBaseLevelForKey(const Slice& user_key);

  // Returns true iff we should stop building the current output
  // before processing "internal_key".
  bool ShouldStopBefore(const Slice& internal_key);

  // Release the input version for the compaction, once the compaction
  // is successful.
  void ReleaseInputs();

 private:
  friend class Version;
  friend class VersionSet;

  explicit Compaction(int level);

  int level_; //// 要compact的level
  uint64_t max_output_file_size_; //// 生成sstable的最大size (kTargetFileSize)
  Version* input_version_;  //// compact时当前的Version
  VersionEdit edit_;   //// 记录compact过程中的操作

  //// inputs_[0]为level-n的sstable文件信息，
  //// inputs_[1]为level-n+1的sstable文件信息
  // Each compaction reads inputs from "level_" and "level_+1"
  std::vector<FileMetaData*> inputs_[2];      // The two sets of inputs

  //// 位于level-n+2，并且与compact的key-range有overlap的sstable。
  //// 保存grandparents_是因为compact最终会生成一系列level-n+1的sstable，
  //// 而如果生成的sstable与level-n+2中有过多的overlap的话，当compact
  //// level-n+1时，会产生过多的merge，为了尽量避免这种情况，compact过程中
  //// 需要检查与level-n+2中产生overlap的size并与
  //// 阈值kMaxGrandParentOverlapBytes做比较，
  //// 以便提前中止compact。 不中止有什么问题？ 中止后怎么继续？
  // State used to check for number of of overlapping grandparent files
  // (parent == level_ + 1, grandparent == level_ + 2)
  std::vector<FileMetaData*> grandparents_;
  //// // 记录compact时grandparents_中已经overlap的index
  size_t grandparent_index_;  // Index in grandparent_starts_
  /// // 记录是否已经有key检查overlap // 如果是第一次检查，发现有overlap，也不会增加overlapped_bytes_. ？
  bool seen_key_;             // Some output key has been seen
  /// 记录已经overlap的累计size
  int64_t overlapped_bytes_;  // Bytes of overlap between current output
                              // and grandparent files

  /// State for implementing IsBaseLevelForKey
/// compact时，当key的ValueType是kTypeDeletion时，
/// 要检查其在level-n+1以上是否存在（IsBaseLevelForKey()）
/// 来决定是否丢弃掉该key。因为compact时，key的遍历是顺序的，
/// 所以每次检查从上一次检查结束的地方开始即可，
/// level_ptrs_[i]中就记录了input_version_->levels_[i]中，上一次比较结束的
  // level_ptrs_ holds indices into input_version_->levels_: our state
  // is that we are positioned at one of the file ranges for each
  // higher level than the ones involved in this compaction (i.e. for
  // all L >= level_ + 2).
  size_t level_ptrs_[config::kNumLevels];
};

}  // namespace leveldb

#endif  // STORAGE_LEVELDB_DB_VERSION_SET_H_
