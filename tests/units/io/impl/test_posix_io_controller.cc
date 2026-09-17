#include <stdexcept>

#include "test_io_controller_common.h"

namespace funtides::io::test {

class PosixIOControllerTest : public IOControllerFixture {
 protected:
  std::unique_ptr<IOControllerBase> open(OpenMode mode, const IOConfig& cfg) const {
    return makeIOController(BackendKind::kPosix, mode, cfg);
  }
};

// ============================================================================
// Snapshot round trip
// ============================================================================

TEST_F(PosixIOControllerTest, SnapshotRoundTrip) {
  const IOConfig cfg = makeConfig();
  const std::size_t n = 4 * 5 * 6;
  const HostVectorReal src = makeField(n);

  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(src, 42, 0.42f);
    io->close();
  }

  HostVectorReal dst("dst", n);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 0);
  io->close();

  expectEqual(src, dst);
}

// Snapshots are addressed by ordinal, not by timestep. Reading three of them
// back out of order is what proves the indexing is right.
TEST_F(PosixIOControllerTest, MultipleSnapshotsKeepTheirOrder) {
  const IOConfig cfg = makeConfig();
  const std::size_t n = 4 * 5 * 6;

  std::vector<HostVectorReal> src;
  for (int k = 0; k < 3; ++k) src.push_back(makeField(n, 1000.0f * k));

  {
    auto io = open(OpenMode::kWrite, cfg);
    for (int k = 0; k < 3; ++k) io->writeSnapshot(src[k], k * 10, k * 0.1f);
    io->close();
  }

  auto io = open(OpenMode::kRead, cfg);
  for (int k : {2, 0, 1}) {
    HostVectorReal dst("dst", n);
    io->readSnapshot(dst, static_cast<std::size_t>(k));
    expectEqual(src[k], dst);
  }
  io->close();
}

TEST_F(PosixIOControllerTest, SingleElementField) {
  const IOConfig cfg = makeConfig({1, 1, 1});
  const HostVectorReal src = makeField(1, 3.5f);

  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(src, 0, 0.0f);
    io->close();
  }

  HostVectorReal dst("dst", 1);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 0);
  io->close();

  expectEqual(src, dst);
}

// Large enough to cross the usual stdio buffer boundaries.
TEST_F(PosixIOControllerTest, LargeField) {
  const IOConfig cfg = makeConfig({64, 64, 64});
  const std::size_t n = 64 * 64 * 64;
  const HostVectorReal src = makeField(n);

  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(src, 0, 0.0f);
    io->close();
  }

  HostVectorReal dst("dst", n);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 0);
  io->close();

  expectEqual(src, dst);
}

// ============================================================================
// Lifecycle
// ============================================================================

// close() must be idempotent: the destructor calls it too.
TEST_F(PosixIOControllerTest, CloseIsIdempotent) {
  const IOConfig cfg = makeConfig();
  auto io = open(OpenMode::kWrite, cfg);
  io->writeSnapshot(makeField(4 * 5 * 6), 0, 0.0f);
  io->close();
  EXPECT_NO_THROW(io->close());
}

// Destroying without an explicit close() must still leave a readable file.
// This is the case that silently loses the last snapshot of a run if the
// destructor forgets to flush.
TEST_F(PosixIOControllerTest, DestructorFlushes) {
  const IOConfig cfg = makeConfig();
  const std::size_t n = 4 * 5 * 6;
  const HostVectorReal src = makeField(n);

  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(src, 0, 0.0f);
  }

  HostVectorReal dst("dst", n);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 0);
  io->close();

  expectEqual(src, dst);
}

TEST_F(PosixIOControllerTest, FlushBetweenWrites) {
  const IOConfig cfg = makeConfig();
  const std::size_t n = 4 * 5 * 6;
  const HostVectorReal a = makeField(n, 1.0f);
  const HostVectorReal b = makeField(n, 2.0f);

  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(a, 0, 0.0f);
    io->flush();
    io->writeSnapshot(b, 1, 0.1f);
    io->close();
  }

  HostVectorReal dst("dst", n);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 1);
  io->close();

  expectEqual(b, dst);
}

// The buffer contract the interface promises: the caller may overwrite its view
// as soon as writeSnapshot() returns. Trivial here, but this is the test that
// will catch a missing staging copy in the asynchronous ADIOS2 backend.
TEST_F(PosixIOControllerTest, CallerMayReuseBufferImmediately) {
  IOConfig cfg = makeConfig();
  cfg.async_snapshots = true;
  const std::size_t n = 4 * 5 * 6;

  const HostVectorReal expected = makeField(n, 1.0f);
  HostVectorReal scratch = makeField(n, 1.0f);

  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(scratch, 0, 0.0f);
    for (std::size_t i = 0; i < n; ++i) scratch(i) = -999.0f;
    io->close();
  }

  HostVectorReal dst("dst", n);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 0);
  io->close();

  expectEqual(expected, dst);
}

// ============================================================================
// Error paths
// ============================================================================

TEST_F(PosixIOControllerTest, ReadingAMissingSnapshotThrows) {
  const IOConfig cfg = makeConfig();
  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(makeField(4 * 5 * 6), 0, 0.0f);
    io->close();
  }

  HostVectorReal dst("dst", 4 * 5 * 6);
  auto io = open(OpenMode::kRead, cfg);
  EXPECT_THROW(io->readSnapshot(dst, 7), std::runtime_error);
}

TEST_F(PosixIOControllerTest, ReadingIntoAMissizedViewThrows) {
  const IOConfig cfg = makeConfig();
  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(makeField(4 * 5 * 6), 0, 0.0f);
    io->close();
  }

  HostVectorReal too_small("dst", 10);
  auto io = open(OpenMode::kRead, cfg);
  EXPECT_THROW(io->readSnapshot(too_small, 0), std::runtime_error);
}

TEST_F(PosixIOControllerTest, WritingOnAReadControllerThrows) {
  const IOConfig cfg = makeConfig();
  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(makeField(4 * 5 * 6), 0, 0.0f);
    io->close();
  }

  auto io = open(OpenMode::kRead, cfg);
  EXPECT_THROW(io->writeSnapshot(makeField(4 * 5 * 6), 1, 0.1f), std::runtime_error);
}

TEST_F(PosixIOControllerTest, ReadingOnAWriteControllerThrows) {
  const IOConfig cfg = makeConfig();
  auto io = open(OpenMode::kWrite, cfg);
  HostVectorReal dst("dst", 4 * 5 * 6);
  EXPECT_THROW(io->readSnapshot(dst, 0), std::runtime_error);
}

TEST_F(PosixIOControllerTest, WritingAfterCloseThrows) {
  const IOConfig cfg = makeConfig();
  auto io = open(OpenMode::kWrite, cfg);
  io->close();
  EXPECT_THROW(io->writeSnapshot(makeField(4 * 5 * 6), 0, 0.0f), std::runtime_error);
}

TEST_F(PosixIOControllerTest, ReadingFromAnEmptyDirectoryThrows) {
  const IOConfig cfg = makeConfig();
  auto io = open(OpenMode::kRead, cfg);
  HostVectorReal dst("dst", 4 * 5 * 6);
  EXPECT_THROW(io->readSnapshot(dst, 0), std::runtime_error);
}

}  // namespace funtides::io::test
