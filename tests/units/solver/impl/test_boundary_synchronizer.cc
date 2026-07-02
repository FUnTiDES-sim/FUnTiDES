#include <gtest/gtest.h>

#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

#include "boundary_synchronizer.h"
#include "parallel_topology.h"

namespace solver {
namespace fe {
namespace test {

struct FloatView {
  std::vector<float> data;
  explicit FloatView(int n, float val = 0.0f) : data(n, val) {}
  float& operator()(int i) { return data[i]; }
  float operator()(int i) const { return data[i]; }
};

class WrongSizeBackend : public BoundarySynchronizer::Backend {
 public:
  void exchange(const std::map<int, std::vector<float>>& sendBuffers,
                std::map<int, std::vector<float>>& recvBuffers) override {
    recvBuffers.clear();
    for (const auto& [rank, data] : sendBuffers) recvBuffers[rank].resize(data.size() + 1, 0.0f);
  }
};

// ======================================================================
// Construction
// ======================================================================

TEST(BoundarySynchronizerTest, NullBackendThrows) {
  EXPECT_THROW(BoundarySynchronizer(nullptr), std::invalid_argument);
}

// ======================================================================
// SerialBackend
// ======================================================================

TEST(SerialBackendTest, ExchangeClearsRecvBuffers) {
  SerialBackend backend;
  std::map<int, std::vector<float>> send = {{1, {1.0f, 2.0f}}};
  std::map<int, std::vector<float>> recv = {{99, {9.0f}}};
  backend.exchange(send, recv);
  EXPECT_TRUE(recv.empty());
}

TEST(BoundarySynchronizerTest, SynchronizeNoop_SingleRank) {
  utils::ParallelTopology topo;
  topo.numRanks = 1;

  FloatView field(5, 3.0f);
  BoundarySynchronizer sync(std::make_unique<SerialBackend>());
  sync.synchronize(field, topo);

  for (int i = 0; i < 5; ++i) EXPECT_FLOAT_EQ(field(i), 3.0f);
}

// ======================================================================
// DebugBackend
// ======================================================================

TEST(DebugBackendTest, ExchangeReturnsZeroBuffers) {
  DebugBackend backend(0);
  std::map<int, std::vector<float>> send = {{1, {1.0f, 2.0f, 3.0f}}};
  std::map<int, std::vector<float>> recv;
  backend.exchange(send, recv);

  ASSERT_NE(recv.find(1), recv.end());
  EXPECT_EQ(recv[1].size(), 3u);
  for (float v : recv[1]) EXPECT_FLOAT_EQ(v, 0.0f);
}

TEST(BoundarySynchronizerTest, SynchronizeWithDebugBackend_AccumulatesZeros) {
  utils::ParallelTopology topo;
  topo.myRank = 0;
  topo.numRanks = 2;
  topo.sharedNodes[1] = {2, 4};

  FloatView field(6, 1.0f);
  BoundarySynchronizer sync(std::make_unique<DebugBackend>(0));
  sync.synchronize(field, topo);

  // DebugBackend returns zeros → accumulate adds 0 → values unchanged
  for (int i = 0; i < 6; ++i) EXPECT_FLOAT_EQ(field(i), 1.0f);
}

// ======================================================================
// Error paths (both wrapped as runtime_error by synchronize catch block)
// ======================================================================

TEST(BoundarySynchronizerTest, MissingNeighborDataThrows) {
  // SerialBackend clears recvBuffers; distributed topo expects neighbor data → throws
  utils::ParallelTopology topo;
  topo.myRank = 0;
  topo.numRanks = 2;
  topo.sharedNodes[1] = {0, 1};

  FloatView field(4, 0.0f);
  BoundarySynchronizer sync(std::make_unique<SerialBackend>());
  EXPECT_THROW(sync.synchronize(field, topo), std::runtime_error);
}

TEST(BoundarySynchronizerTest, WrongBufferSizeThrows) {
  // length_error in accumulate is caught and re-wrapped as runtime_error
  utils::ParallelTopology topo;
  topo.myRank = 0;
  topo.numRanks = 2;
  topo.sharedNodes[1] = {0, 1};

  FloatView field(4, 0.0f);
  BoundarySynchronizer sync(std::make_unique<WrongSizeBackend>());
  EXPECT_THROW(sync.synchronize(field, topo), std::runtime_error);
}

}  // namespace test
}  // namespace fe
}  // namespace solver
