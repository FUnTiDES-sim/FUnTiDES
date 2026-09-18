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
    io->writeSnapshot(src);
    io->close();
  }

  HostVectorReal dst("dst", n);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 0);
  io->close();

  expectEqual(src, dst);
}

// Snapshots are addressed by the order they were written. Reading three of them
// back out of order is what proves the indexing is right.
TEST_F(PosixIOControllerTest, MultipleSnapshotsKeepTheirOrder) {
  const IOConfig cfg = makeConfig();
  const std::size_t n = 4 * 5 * 6;

  std::vector<HostVectorReal> src;
  for (int k = 0; k < 3; ++k) src.push_back(makeField(n, 1000.0f * k));

  {
    auto io = open(OpenMode::kWrite, cfg);
    for (int k = 0; k < 3; ++k) io->writeSnapshot(src[k]);
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
    io->writeSnapshot(src);
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
    io->writeSnapshot(src);
    io->close();
  }

  HostVectorReal dst("dst", n);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 0);
  io->close();

  expectEqual(src, dst);
}

// ============================================================================
// Shot separation
// ============================================================================

// A shot id sends output to its own subdirectory, so two shots sharing an
// output_dir never overwrite each other.
TEST_F(PosixIOControllerTest, ShotIdSeparatesOutput) {
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

  EXPECT_TRUE(std::filesystem::is_directory(dir_ / "shot_042"));
  EXPECT_TRUE(std::filesystem::is_directory(dir_ / "shot_043"));

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

// shot_id becomes a path component, so anything that could escape the output
// directory is rejected rather than quietly sanitized.
TEST_F(PosixIOControllerTest, InvalidShotIdThrows) {
  for (const char* bad : {"../escape", "sub/dir", "with space", "dollar$"}) {
    IOConfig cfg = makeConfig();
    cfg.shot_id = bad;
    EXPECT_THROW(open(OpenMode::kWrite, cfg), std::invalid_argument) << "shot_id = " << bad;
  }
}

// ============================================================================
// Lifecycle
// ============================================================================

// close() must be idempotent: the destructor calls it too.
TEST_F(PosixIOControllerTest, CloseIsIdempotent) {
  const IOConfig cfg = makeConfig();
  auto io = open(OpenMode::kWrite, cfg);
  io->writeSnapshot(makeField(4 * 5 * 6));
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
    io->writeSnapshot(src);
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
    io->writeSnapshot(a);
    io->flush();
    io->writeSnapshot(b);
    io->close();
  }

  HostVectorReal dst("dst", n);
  auto io = open(OpenMode::kRead, cfg);
  io->readSnapshot(dst, 1);
  io->close();

  expectEqual(b, dst);
}

// ============================================================================
// Error paths
// ============================================================================

TEST_F(PosixIOControllerTest, ReadingAMissingSnapshotThrows) {
  const IOConfig cfg = makeConfig();
  {
    auto io = open(OpenMode::kWrite, cfg);
    io->writeSnapshot(makeField(4 * 5 * 6));
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
    io->writeSnapshot(makeField(4 * 5 * 6));
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
    io->writeSnapshot(makeField(4 * 5 * 6));
    io->close();
  }

  auto io = open(OpenMode::kRead, cfg);
  EXPECT_THROW(io->writeSnapshot(makeField(4 * 5 * 6)), std::runtime_error);
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
  EXPECT_THROW(io->writeSnapshot(makeField(4 * 5 * 6)), std::runtime_error);
}

TEST_F(PosixIOControllerTest, ReadingFromAnEmptyDirectoryThrows) {
  const IOConfig cfg = makeConfig();
  auto io = open(OpenMode::kRead, cfg);
  HostVectorReal dst("dst", 4 * 5 * 6);
  EXPECT_THROW(io->readSnapshot(dst, 0), std::runtime_error);
}

}  // namespace funtides::io::test
