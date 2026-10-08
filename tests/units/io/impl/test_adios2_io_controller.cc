#include <stdexcept>

#include "test_io_controller_common.h"

namespace funtides::io::test {

// Built in every configuration: without FUNTIDES_ENABLE_ADIOS2 the fixture
// skips its tests, and AdiosAvailability checks the factory error instead.
class Adios2IOControllerTest : public IOControllerFixture {
 protected:
  void SetUp() override {
    IOControllerFixture::SetUp();
    if (!isBackendAvailable(BackendKind::kAdios2)) {
      GTEST_SKIP() << "ADIOS2 backend not built in (FUNTIDES_ENABLE_ADIOS2=OFF)";
    }
  }

  std::unique_ptr<IOControllerBase> open(OpenMode mode, const IOConfig& cfg) const {
    return makeIOController(BackendKind::kAdios2, mode, cfg);
  }

  IOConfig makeAsyncConfig(std::vector<std::size_t> dims = {4, 5, 6}) const {
    IOConfig cfg = makeConfig(std::move(dims));
    cfg.async_snapshots = true;
    return cfg;
  }

  /// Writes `count` snapshots seeded 1000*k and returns them.
  std::vector<HostVectorReal> writeSeries(const IOConfig& cfg, int count, std::size_t n) const {
    std::vector<HostVectorReal> src;
    for (int k = 0; k < count; ++k) src.push_back(makeField(n, 1000.0f * k));
    auto io = open(OpenMode::kWrite, cfg);
    for (const HostVectorReal& s : src) io->writeSnapshot(s);
    io->close();
    return src;
  }
};

// ============================================================================
// Build configuration
// ============================================================================

// When the backend is not compiled in, asking for it must fail loudly rather
// than fall back to another backend.
TEST(AdiosAvailability, FactoryMatchesIsBackendAvailable) {
  EXPECT_TRUE(isBackendAvailable(BackendKind::kPosix));
  if (isBackendAvailable(BackendKind::kAdios2)) GTEST_SKIP() << "ADIOS2 backend is built in";
  EXPECT_THROW(makeIOController(BackendKind::kAdios2, OpenMode::kWrite, IOConfig{}), std::runtime_error);
}

// ============================================================================
// Snapshot round trip
// ============================================================================

TEST_F(Adios2IOControllerTest, SnapshotRoundTrip) {
  const IOConfig cfg = makeConfig();
  const std::size_t n = 4 * 5 * 6;
  const HostVectorReal src = makeField(n);

  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(src);
    io->close();
  }

  HostVectorReal dst("dst", n);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 0);
  io->close();

  expectEqual(src, dst);
}

TEST_F(Adios2IOControllerTest, MultipleSnapshotsKeepTheirOrder) {
  const IOConfig cfg = makeConfig();
  const std::size_t n = 4 * 5 * 6;
  const auto src = writeSeries(cfg, 3, n);

  auto io = open(OpenMode::kRead, cfg);
  for (int k : {2, 0, 1}) {
    HostVectorReal dst("dst", n);
    io->readSnapshot(dst, static_cast<std::size_t>(k));
    expectEqual(src[k], dst);
  }
  io->close();
}

TEST_F(Adios2IOControllerTest, SingleElementField) {
  const IOConfig cfg = makeConfig({1, 1, 1});
  const HostVectorReal src = makeField(1, 3.5f);

  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(src);
    io->close();
  }

  HostVectorReal dst("dst", 1);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 0);
  io->close();

  expectEqual(src, dst);
}

TEST_F(Adios2IOControllerTest, LargeField) {
  const IOConfig cfg = makeAsyncConfig({64, 64, 64});
  const std::size_t n = 64 * 64 * 64;
  const auto src = writeSeries(cfg, 4, n);

  auto io = open(OpenMode::kRead, cfg);
  HostVectorReal dst("dst", n);
  for (int k = 0; k < 4; ++k) {
    io->readSnapshot(dst, static_cast<std::size_t>(k));
    expectEqual(src[k], dst);
  }
  io->close();
}

// ============================================================================
// Asynchronous write path
// ============================================================================

// The controller copies the field before writeSnapshot() returns, so the
// caller may overwrite its host mirror right away. This is the property that
// lets a solver reuse one mirror every time step; a backend that kept a
// pointer to the caller's buffer would store the last values in every step.
TEST_F(Adios2IOControllerTest, CallerBufferIsReusableWhenWriteReturns) {
  const IOConfig cfg = makeAsyncConfig();
  const std::size_t n = 4 * 5 * 6;
  constexpr int kCount = 8;

  {
    auto io = open(OpenMode::kWrite, cfg);
    HostVectorReal mirror("mirror", n);
    for (int k = 0; k < kCount; ++k) {
      const HostVectorReal values = makeField(n, 1000.0f * k);
      for (std::size_t i = 0; i < n; ++i) mirror(i) = values(i);
      io->writeSnapshot(mirror);
    }
    io->close();
  }

  auto io = open(OpenMode::kRead, cfg);
  HostVectorReal dst("dst", n);
  for (int k = 0; k < kCount; ++k) {
    io->readSnapshot(dst, static_cast<std::size_t>(k));
    expectEqual(makeField(n, 1000.0f * k), dst);
  }
  io->close();
}

// More snapshots than staging buffers: writeSnapshot() must block on a free
// buffer, not drop or reorder snapshots.
TEST_F(Adios2IOControllerTest, ManyAsyncSnapshotsKeepTheirOrder) {
  const IOConfig cfg = makeAsyncConfig();
  const std::size_t n = 4 * 5 * 6;
  constexpr int kCount = 32;
  const auto src = writeSeries(cfg, kCount, n);

  auto io = open(OpenMode::kRead, cfg);
  HostVectorReal dst("dst", n);
  for (int k = 0; k < kCount; ++k) {
    io->readSnapshot(dst, static_cast<std::size_t>(k));
    expectEqual(src[k], dst);
  }
  io->close();
}

// After flush() every queued snapshot is on disk: a reader opened while the
// writer is still alive must see all of them.
TEST_F(Adios2IOControllerTest, FlushBetweenWrites) {
  const IOConfig cfg = makeAsyncConfig();
  const std::size_t n = 4 * 5 * 6;
  const HostVectorReal a = makeField(n, 1.0f);
  const HostVectorReal b = makeField(n, 2.0f);

  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(a);
    io->flush();
    io->writeSnapshot(b);
    io->close();
  }

  HostVectorReal dst("dst", n);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 0);
  expectEqual(a, dst);
  io->readSnapshot(dst, 1);
  expectEqual(b, dst);
  io->close();
}

// ============================================================================
// Prefetching read path
// ============================================================================

// Backward sweep, as in the adjoint pass of an RTM: prefetching must follow
// the access direction and never hand back a neighbouring snapshot.
TEST_F(Adios2IOControllerTest, BackwardSweepReturnsEverySnapshot) {
  const IOConfig cfg = makeAsyncConfig();
  const std::size_t n = 4 * 5 * 6;
  constexpr int kCount = 10;
  const auto src = writeSeries(cfg, kCount, n);

  auto io = open(OpenMode::kRead, cfg);
  HostVectorReal dst("dst", n);
  for (int k = kCount - 1; k >= 0; --k) {
    io->readSnapshot(dst, static_cast<std::size_t>(k));
    expectEqual(src[k], dst);
  }
  io->close();
}

// Jumps, repeats and direction changes invalidate the prefetched slots in
// every way the slot logic has to handle.
TEST_F(Adios2IOControllerTest, IrregularAccessPattern) {
  const IOConfig cfg = makeAsyncConfig();
  const std::size_t n = 4 * 5 * 6;
  constexpr int kCount = 10;
  const auto src = writeSeries(cfg, kCount, n);

  auto io = open(OpenMode::kRead, cfg);
  HostVectorReal dst("dst", n);
  for (int k : {0, 1, 2, 9, 8, 8, 3, 0, 9, 5, 6, 4, 4}) {
    io->readSnapshot(dst, static_cast<std::size_t>(k));
    expectEqual(src[k], dst);
  }
  io->close();
}

// ============================================================================
// Shot separation and file layout
// ============================================================================

TEST_F(Adios2IOControllerTest, ShotIdSeparatesOutput) {
  const std::size_t n = 4 * 5 * 6;
  const HostVectorReal src_a = makeField(n, 1.0f);
  const HostVectorReal src_b = makeField(n, 2.0f);

  IOConfig cfg_a = makeConfig();
  cfg_a.shot_id = "shot_042";
  IOConfig cfg_b = makeConfig();
  cfg_b.shot_id = "shot_043";

  {
    auto io = open(OpenMode::kWrite, cfg_a);
    io->writeSnapshot(src_a);
    io->close();
  }
  {
    auto io = open(OpenMode::kWrite, cfg_b);
    io->writeSnapshot(src_b);
    io->close();
  }

  // One BP file (a directory for BP5) per rank inside each shot directory.
  EXPECT_TRUE(std::filesystem::exists(dir_ / "shot_042" / "test_r0000.bp"));
  EXPECT_TRUE(std::filesystem::exists(dir_ / "shot_043" / "test_r0000.bp"));

  HostVectorReal dst("dst", n);
  auto io_a = open(OpenMode::kRead, cfg_a);
  io_a->readSnapshot(dst, 0);
  io_a->close();
  expectEqual(src_a, dst);

  auto io_b = open(OpenMode::kRead, cfg_b);
  io_b->readSnapshot(dst, 0);
  io_b->close();
  expectEqual(src_b, dst);
}

// The rank is part of the file name: two ranks sharing a directory must not
// overwrite each other.
TEST_F(Adios2IOControllerTest, RankSeparatesOutput) {
  const std::size_t n = 4 * 5 * 6;
  IOConfig cfg_0 = makeConfig();
  IOConfig cfg_1 = makeConfig();
  cfg_1.ctx.rank = 1;

  const HostVectorReal src_0 = makeField(n, 10.0f);
  const HostVectorReal src_1 = makeField(n, 20.0f);
  {
    auto io_0 = open(OpenMode::kWrite, cfg_0);
    auto io_1 = open(OpenMode::kWrite, cfg_1);
    io_0->writeSnapshot(src_0);
    io_1->writeSnapshot(src_1);
    io_0->close();
    io_1->close();
  }

  HostVectorReal dst("dst", n);
  auto io_1 = open(OpenMode::kRead, cfg_1);
  io_1->readSnapshot(dst, 0);
  io_1->close();
  expectEqual(src_1, dst);

  auto io_0 = open(OpenMode::kRead, cfg_0);
  io_0->readSnapshot(dst, 0);
  io_0->close();
  expectEqual(src_0, dst);
}

TEST_F(Adios2IOControllerTest, InvalidShotIdThrows) {
  for (const char* bad : {"../escape", "sub/dir", "with space", "dollar$"}) {
    IOConfig cfg = makeConfig();
    cfg.shot_id = bad;
    EXPECT_THROW(open(OpenMode::kWrite, cfg), std::invalid_argument) << "shot_id = " << bad;
  }
}

// ============================================================================
// Lifecycle
// ============================================================================

TEST_F(Adios2IOControllerTest, CloseIsIdempotent) {
  const IOConfig cfg = makeConfig();
  auto io = open(OpenMode::kWrite, cfg);
  io->writeSnapshot(makeField(4 * 5 * 6));
  io->close();
  EXPECT_NO_THROW(io->close());
}

// With async writes still queued, the destructor must drain them before
// closing: otherwise the last snapshots of a run are silently lost.
TEST_F(Adios2IOControllerTest, DestructorDrainsQueuedWrites) {
  const IOConfig cfg = makeAsyncConfig();
  const std::size_t n = 4 * 5 * 6;
  constexpr int kCount = 6;

  std::vector<HostVectorReal> src;
  for (int k = 0; k < kCount; ++k) src.push_back(makeField(n, 1000.0f * k));
  {
    auto io = open(OpenMode::kWrite, cfg);
    for (const HostVectorReal& s : src) io->writeSnapshot(s);
  }

  auto io = open(OpenMode::kRead, cfg);
  HostVectorReal dst("dst", n);
  for (int k = 0; k < kCount; ++k) {
    io->readSnapshot(dst, static_cast<std::size_t>(k));
    expectEqual(src[k], dst);
  }
  io->close();
}

// ============================================================================
// Error paths
// ============================================================================

TEST_F(Adios2IOControllerTest, ReadingAMissingSnapshotThrows) {
  const IOConfig cfg = makeConfig();
  writeSeries(cfg, 1, 4 * 5 * 6);

  HostVectorReal dst("dst", 4 * 5 * 6);
  auto io = open(OpenMode::kRead, cfg);
  EXPECT_THROW(io->readSnapshot(dst, 7), std::runtime_error);
  EXPECT_THROW(io->readSnapshot(dst, 1), std::runtime_error);  // one past the end
}

TEST_F(Adios2IOControllerTest, ReadingIntoAMissizedViewThrows) {
  const IOConfig cfg = makeConfig();
  writeSeries(cfg, 1, 4 * 5 * 6);

  HostVectorReal too_small("dst", 10);
  auto io = open(OpenMode::kRead, cfg);
  EXPECT_THROW(io->readSnapshot(too_small, 0), std::runtime_error);
}

// A BP variable has one shape for all its steps.
TEST_F(Adios2IOControllerTest, ChangingTheSnapshotSizeThrows) {
  const IOConfig cfg = makeAsyncConfig();
  auto io = open(OpenMode::kWrite, cfg);
  io->writeSnapshot(makeField(4 * 5 * 6));
  EXPECT_THROW(io->writeSnapshot(makeField(10)), std::runtime_error);
  EXPECT_NO_THROW(io->close());  // the valid snapshot is still written
}

TEST_F(Adios2IOControllerTest, WritingOnAReadControllerThrows) {
  const IOConfig cfg = makeConfig();
  writeSeries(cfg, 1, 4 * 5 * 6);

  auto io = open(OpenMode::kRead, cfg);
  EXPECT_THROW(io->writeSnapshot(makeField(4 * 5 * 6)), std::runtime_error);
}

TEST_F(Adios2IOControllerTest, ReadingOnAWriteControllerThrows) {
  const IOConfig cfg = makeConfig();
  auto io = open(OpenMode::kWrite, cfg);
  HostVectorReal dst("dst", 4 * 5 * 6);
  EXPECT_THROW(io->readSnapshot(dst, 0), std::runtime_error);
}

TEST_F(Adios2IOControllerTest, WritingAfterCloseThrows) {
  const IOConfig cfg = makeConfig();
  auto io = open(OpenMode::kWrite, cfg);
  io->close();
  EXPECT_THROW(io->writeSnapshot(makeField(4 * 5 * 6)), std::runtime_error);
}

TEST_F(Adios2IOControllerTest, ReadingAfterCloseThrows) {
  const IOConfig cfg = makeConfig();
  writeSeries(cfg, 1, 4 * 5 * 6);

  auto io = open(OpenMode::kRead, cfg);
  io->close();
  HostVectorReal dst("dst", 4 * 5 * 6);
  EXPECT_THROW(io->readSnapshot(dst, 0), std::runtime_error);
}

// Unlike the POSIX backend, which opens one file per snapshot and so can only
// fail in readSnapshot(), this backend opens its single file at construction.
TEST_F(Adios2IOControllerTest, OpeningAMissingFileForReadingThrows) {
  const IOConfig cfg = makeConfig();
  EXPECT_THROW(open(OpenMode::kRead, cfg), std::runtime_error);
}

}  // namespace funtides::io::test
