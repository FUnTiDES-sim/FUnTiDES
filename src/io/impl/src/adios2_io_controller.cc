#include "adios2_io_controller.h"

#include <adios2.h>

#include <algorithm>
#include <array>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <exception>
#include <filesystem>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

namespace funtides::io {
namespace {

using Real = HostVectorReal::non_const_value_type;
using Buffer = std::vector<Real>;

constexpr char kVariableName[] = "field";
constexpr char kEngine[] = "BP5";

/// Same rule as the POSIX backend: `shot_id` becomes a path component.
/// TODO: move to a shared internal header with the POSIX copy.
void validateShotId(const std::string& shot_id) {
  for (const char c : shot_id) {
    const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '-';
    if (!ok) {
      throw std::invalid_argument("funtides::io: shot_id must contain only alphanumerics, '_' or '-', got \"" +
                                  shot_id + "\"");
    }
  }
}

std::string outputDir(const IOConfig& cfg) {
  return cfg.shot_id.empty() ? cfg.output_dir : cfg.output_dir + "/" + cfg.shot_id;
}

std::string bpPath(const IOConfig& cfg) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "_r%04d.bp", cfg.ctx.rank);
  return outputDir(cfg) + "/" + cfg.prefix + buf;
}

/// The staging copies are plain memcpy of `size()` values, which is only
/// correct for a contiguous view.
void requireContiguous(const HostVectorReal& field) {
  if (!field.span_is_contiguous()) {
    throw std::runtime_error("funtides::io: the ADIOS2 backend needs a contiguous view");
  }
}

/// ADIOS2 rejects zero-length array attributes, so empty vectors are skipped.
void defineDimsAttribute(adios2::IO& io, const std::string& name, const std::vector<std::size_t>& v) {
  if (v.empty()) return;
  const std::vector<std::uint64_t> data(v.begin(), v.end());
  io.DefineAttribute<std::uint64_t>(name, data.data(), data.size());
}

}  // namespace

// -----------------------------------------------------------------------------
// Impl: the state shared between the caller's thread and the I/O thread.
//
// One mutex `mu_` protects every member below it, one condition variable `cv_`
// signals every change (notify_all; at most two threads ever wait on it).
// Rule: ADIOS2 objects live only on the I/O thread's stack; the caller's
// thread never touches them.
// -----------------------------------------------------------------------------
class Adios2IOController::Impl {
 public:
  Impl(OpenMode mode, const IOConfig& config) : mode_(mode), cfg_(config), path_(bpPath(config)) {}

  /// Joins the thread. Must never throw: called from destructors.
  ~Impl() { stop(); }

  Impl(const Impl&) = delete;
  Impl& operator=(const Impl&) = delete;

  /// Starts the I/O thread and waits until it has opened the file.
  void start() {
    if (mode_ == OpenMode::kWrite) {
      pool_.resize(kWriteBuffers);
      for (Buffer& b : pool_) free_.push_back(&b);
      thread_ = std::thread([this] { runWriter(); });
    } else {
      thread_ = std::thread([this] { runReader(); });
    }
    std::unique_lock lk(mu_);
    cv_.wait(lk, [&] { return opened_ || error_; });
    rethrowIfFailed();
  }

  /// Asks the thread to finish and joins it. A writer first drains its queue;
  /// a reader drops its pending prefetches.
  void stop() noexcept {
    {
      std::lock_guard lk(mu_);
      stop_ = true;
    }
    cv_.notify_all();
    if (thread_.joinable()) thread_.join();
  }

  void rethrowIfFailed() const {
    if (error_) std::rethrow_exception(error_);
  }

  void checkFailed() {
    std::lock_guard lk(mu_);
    rethrowIfFailed();
  }

  // ---------------------------------------------------------------- write --

  void write(const Real* src, std::size_t n, bool wait_landed) {
    std::unique_lock lk(mu_);
    rethrowIfFailed();
    if (!nelem_) {
      nelem_ = n;
    } else if (*nelem_ != n) {
      throw std::runtime_error("funtides::io: snapshot of " + std::to_string(n) + " values, previous ones had " +
                               std::to_string(*nelem_));
    }

    // Back-pressure: wait for a free staging buffer.
    cv_.wait(lk, [&] { return !free_.empty() || error_; });
    rethrowIfFailed();
    Buffer* buf = free_.front();
    free_.pop_front();

    // The buffer is owned by this thread until it is queued: copy unlocked so
    // the I/O thread can keep writing meanwhile.
    lk.unlock();
    buf->assign(src, src + n);
    lk.lock();

    jobs_.push_back(buf);
    cv_.notify_all();
    if (wait_landed) waitIdle(lk);
  }

  void flushWrites() {
    std::unique_lock lk(mu_);
    waitIdle(lk);
  }

  // ----------------------------------------------------------------- read --

  void read(Real* dst, std::size_t n, std::size_t index) {
    std::unique_lock lk(mu_);
    rethrowIfFailed();
    if (index >= nsteps_) {
      throw std::runtime_error("funtides::io: no snapshot with index " + std::to_string(index) + " in " + path_ +
                               " (" + std::to_string(nsteps_) + " stored)");
    }
    if (n != *nelem_) {
      throw std::runtime_error("funtides::io: " + path_ + " holds " + std::to_string(*nelem_) +
                               " values per snapshot, the view has " + std::to_string(n));
    }

    Slot* s = findSlot(index);
    if (s == nullptr) {
      s = evictSlot();
      s->index = index;
      s->state = SlotState::kQueued;
      requests_.push_front(s);  // demand read goes ahead of prefetches
      cv_.notify_all();
    } else if (s->state == SlotState::kQueued) {
      // Prefetch not started yet: move it to the front of the queue.
      requests_.erase(std::find(requests_.begin(), requests_.end(), s));
      requests_.push_front(s);
    }
    cv_.wait(lk, [&] { return s->state == SlotState::kReady || error_; });
    rethrowIfFailed();

    // A kReady slot is touched by this thread only, so copy unlocked.
    lk.unlock();
    std::copy(s->data.begin(), s->data.end(), dst);
    lk.lock();
    s->state = SlotState::kEmpty;

    prefetchAfter(index);
  }

 private:
  enum class SlotState { kEmpty, kQueued, kLoading, kReady };

  /// One prefetch slot. `data` is sized once, after the file is opened.
  struct Slot {
    SlotState state{SlotState::kEmpty};
    std::size_t index{0};
    Buffer data;
  };
  // One slot for the demand read plus kPrefetchDepth ahead. At most one slot
  // is kLoading at any time (single I/O thread), so evictSlot() always finds
  // a victim when there are at least 2 slots.
  static constexpr std::size_t kSlots = kPrefetchDepth + 1;
  static_assert(kSlots >= 2);

  /// Waits until the write queue is empty and the thread is not writing.
  void waitIdle(std::unique_lock<std::mutex>& lk) {
    cv_.wait(lk, [&] { return (jobs_.empty() && !busy_) || error_; });
    rethrowIfFailed();
  }

  Slot* findSlot(std::size_t index) {
    for (Slot& s : slots_) {
      if (s.state != SlotState::kEmpty && s.index == index) return &s;
    }
    return nullptr;
  }

  /// Returns a slot for a demand read: an empty one if any, else a queued or
  /// ready prefetch that is dropped. Never the slot being loaded.
  Slot* evictSlot() {
    for (Slot& s : slots_) {
      if (s.state == SlotState::kEmpty) return &s;
    }
    for (Slot& s : slots_) {
      if (s.state == SlotState::kQueued) {
        requests_.erase(std::find(requests_.begin(), requests_.end(), &s));
        return &s;
      }
    }
    for (Slot& s : slots_) {
      if (s.state == SlotState::kReady) return &s;
    }
    throw std::logic_error("funtides::io: no evictable prefetch slot");  // unreachable, see kSlots
  }

  /// Queues the next kPrefetchDepth snapshots in the access direction, using
  /// empty slots only: a prefetch never evicts another one.
  void prefetchAfter(std::size_t index) {
    const bool backward = last_read_ && index < *last_read_;
    last_read_ = index;
    for (std::size_t k = 1; k <= kPrefetchDepth; ++k) {
      if (backward && k > index) break;
      const std::size_t next = backward ? index - k : index + k;
      if (next >= nsteps_) break;
      if (findSlot(next) != nullptr) continue;
      Slot* s = nullptr;
      for (Slot& c : slots_) {
        if (c.state == SlotState::kEmpty) {
          s = &c;
          break;
        }
      }
      if (s == nullptr) break;
      s->index = next;
      s->state = SlotState::kQueued;
      requests_.push_back(s);
    }
    cv_.notify_all();
  }

  /// Records the first error of the I/O thread and wakes every waiter.
  void fail(std::exception_ptr e) {
    {
      std::lock_guard lk(mu_);
      if (!error_) error_ = e;
    }
    cv_.notify_all();
  }

  void markOpened() {
    {
      std::lock_guard lk(mu_);
      opened_ = true;
    }
    cv_.notify_all();
  }

  // ------------------------------------------------------ I/O thread side --

  void runWriter() {
    std::optional<adios2::Engine> engine;
    try {
      adios2::ADIOS adios;
      adios2::IO io = adios.DeclareIO("funtides_snapshots");
      io.SetEngine(kEngine);
      io.DefineAttribute<std::int32_t>("rank", static_cast<std::int32_t>(cfg_.ctx.rank));
      defineDimsAttribute(io, "local_dims", cfg_.local_dims);
      defineDimsAttribute(io, "global_dims", cfg_.global_dims);
      defineDimsAttribute(io, "start_offsets", cfg_.start_offsets);
      engine = io.Open(path_, adios2::Mode::Write);
      markOpened();

      adios2::Variable<Real> var;
      for (;;) {
        Buffer* buf = nullptr;
        {
          std::unique_lock lk(mu_);
          cv_.wait(lk, [&] { return stop_ || !jobs_.empty(); });
          if (jobs_.empty()) break;  // stop_ requested and queue drained
          buf = jobs_.front();
          jobs_.pop_front();
          busy_ = true;
        }

        if (!var) {
          const std::size_t n = buf->size();
          var = io.DefineVariable<Real>(kVariableName, {n}, {0}, {n}, adios2::ConstantDims);
        }
        engine->BeginStep();
        // Deferred: ADIOS2 reads `buf` during EndStep, and `buf` stays ours
        // until we hand it back below, so no extra internal copy is needed.
        engine->Put(var, buf->data(), adios2::Mode::Deferred);
        engine->EndStep();

        {
          std::lock_guard lk(mu_);
          busy_ = false;
          free_.push_back(buf);
        }
        cv_.notify_all();
      }
      engine->Close();
    } catch (...) {
      fail(std::current_exception());
      // Best-effort close so the data already written stays readable.
      if (engine && *engine) {
        try {
          engine->Close();
        } catch (...) {
        }
      }
    }
  }

  void runReader() {
    try {
      adios2::ADIOS adios;
      adios2::IO io = adios.DeclareIO("funtides_snapshots");
      io.SetEngine(kEngine);
      adios2::Engine engine = io.Open(path_, adios2::Mode::ReadRandomAccess);
      adios2::Variable<Real> var = io.InquireVariable<Real>(kVariableName);
      if (!var) {
        throw std::runtime_error("funtides::io: " + path_ + " has no \"" + kVariableName +
                                 "\" variable of the expected scalar type");
      }
      const adios2::Dims shape = var.Shape();
      if (shape.size() != 1) {
        throw std::runtime_error("funtides::io: \"" + std::string(kVariableName) + "\" in " + path_ +
                                 " is not 1D");
      }
      {
        std::lock_guard lk(mu_);
        nsteps_ = var.Steps();
        nelem_ = shape[0];
        for (Slot& s : slots_) s.data.resize(shape[0]);
      }
      markOpened();

      for (;;) {
        Slot* s = nullptr;
        {
          std::unique_lock lk(mu_);
          cv_.wait(lk, [&] { return stop_ || !requests_.empty(); });
          if (stop_) break;
          s = requests_.front();
          requests_.pop_front();
          s->state = SlotState::kLoading;
        }

        var.SetStepSelection({s->index, 1});
        engine.Get(var, s->data.data(), adios2::Mode::Sync);

        {
          std::lock_guard lk(mu_);
          s->state = SlotState::kReady;
        }
        cv_.notify_all();
      }
      engine.Close();
    } catch (...) {
      fail(std::current_exception());
    }
  }

  const OpenMode mode_;
  const IOConfig cfg_;  // copy: the thread must not depend on the caller's object
  const std::string path_;
  std::thread thread_;

  std::mutex mu_;
  std::condition_variable cv_;
  bool opened_{false};
  bool stop_{false};
  std::exception_ptr error_;
  std::optional<std::size_t> nelem_;  ///< Values per snapshot, set by the first write or by open (read).

  // Write side.
  std::vector<Buffer> pool_;   ///< Staging buffers; never resized after start().
  std::deque<Buffer*> free_;   ///< Buffers the caller may fill.
  std::deque<Buffer*> jobs_;   ///< Buffers waiting to be written, in snapshot order.
  bool busy_{false};           ///< True while the thread writes a buffer.

  // Read side.
  std::size_t nsteps_{0};
  std::array<Slot, kSlots> slots_;
  std::deque<Slot*> requests_;  ///< Slots in state kQueued, in load order.
  std::optional<std::size_t> last_read_;
};

// -----------------------------------------------------------------------------

Adios2IOController::Adios2IOController(OpenMode mode, const IOConfig& config)
    : IOControllerBase(config), mode_(mode) {
  validateShotId(config.shot_id);
  const std::string path = bpPath(config);

  if (mode_ == OpenMode::kWrite) {
    // Same caveat as the POSIX backend on a parallel filesystem.
    const std::string dir = outputDir(config);
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) throw std::runtime_error("funtides::io: cannot create " + dir + ": " + ec.message());
  } else if (!std::filesystem::exists(path)) {
    throw std::runtime_error("funtides::io: no snapshot file " + path);
  }

  impl_ = std::make_unique<Impl>(mode, config);
  impl_->start();  // if this throws, ~Impl joins the thread
}

Adios2IOController::~Adios2IOController() {
  try {
    close();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "[funtides::io] error during close: %s\n", e.what());
  } catch (...) {
    std::fprintf(stderr, "[funtides::io] unknown error during close\n");
  }
}

void Adios2IOController::requireMode(OpenMode expected, const char* what) const {
  if (mode_ != expected) {
    throw std::runtime_error(std::string("funtides::io: ") + what +
                             " called on a controller opened in the opposite mode");
  }
}

void Adios2IOController::writeSnapshot(const HostVectorReal& field) {
  requireMode(OpenMode::kWrite, "writeSnapshot");
  if (closed_) throw std::runtime_error("funtides::io: writeSnapshot after close");
  requireContiguous(field);
  impl_->write(field.data(), field.size(), !config().async_snapshots);
}

void Adios2IOController::readSnapshot(const HostVectorReal& field, std::size_t index) {
  requireMode(OpenMode::kRead, "readSnapshot");
  if (closed_) throw std::runtime_error("funtides::io: readSnapshot after close");
  requireContiguous(field);
  impl_->read(field.data(), field.size(), index);
}

void Adios2IOController::flush() {
  if (closed_) return;
  if (mode_ == OpenMode::kWrite) {
    impl_->flushWrites();
  } else {
    impl_->checkFailed();
  }
}

void Adios2IOController::close() {
  if (closed_) return;
  closed_ = true;
  impl_->stop();  // a writer drains its queue and closes the file first
  impl_->checkFailed();
}

}  // namespace funtides::io
