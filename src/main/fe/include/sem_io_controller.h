#include <adios2.h>
#include <data_type.h>

#include <cstddef>

#include "adios2/common/ADIOSTypes.h"
#include "adios2/cxx11/Operator.h"

#define RECEIVERS_FILE "receivers.bp"
#define SNAPS_FILE "snapshots.bp"

class SemIOController
{
 private:
  adios2::ADIOS adios_;
  adios2::IO io_;
  adios2::IO async_io_;
  adios2::Engine receiver_writer_;
  adios2::Engine snaps_writer_;
  adios2::Variable<float> receivers_;
  adios2::Variable<float> receivers_coords_;
  adios2::Variable<float> iter_times_;
  adios2::Variable<float> pn_;
  adios2::Variable<int> timestep_;
  adios2::Operator compressor_op_;
  adios2::Operator receiver_op_;

  void initAdios() { adios_ = adios2::ADIOS(); }

  void configureIO()
  {
    io_ = adios_.DeclareIO("AccousticSEMOutput");
    io_.SetEngine("BP5");
    io_.SetParameter("Threads", "4");

    async_io_ = adios_.DeclareIO("AsyncAccousticSEMOutput");
    async_io_.SetEngine("BP5");
    async_io_.SetParameter("AsyncWrite", "On");
    async_io_.SetParameter("Threads", "4");
    async_io_.SetParameter("Profile", "On");
    async_io_.SetParameter("ProfileUnits", "Microseconds");

    // I/O Compression. See function AttachOperator
    // receiver_op_ = adios_.DefineOperator("bloscComp", "blosc");
    // compressor_op_ = adios_.DefineOperator("bloscSnap", "blosc");
  }

  void launchWriters()
  {
    receiver_writer_ = io_.Open(RECEIVERS_FILE, adios2::Mode::Write);
    snaps_writer_ = async_io_.Open(SNAPS_FILE, adios2::Mode::Write);
  }

  void defineVariable(const size_t nb_nodes, const size_t nb_iter, const size_t nb_receiver)
  {
    receivers_ = io_.DefineVariable<float>("AccousticReceiver", {nb_receiver, nb_iter}, {0, 0}, {nb_receiver, nb_iter});

    receivers_coords_ =
        io_.DefineVariable<float>("AccousticReceiverCoords", {nb_receiver, 3}, {0, 0}, {nb_receiver, 3});

    iter_times_ = io_.DefineVariable<float>("IterationTimes", {nb_iter}, {0}, {nb_iter});

    pn_ = async_io_.DefineVariable<float>("PressureField", {nb_nodes}, {0}, {nb_nodes});

    timestep_ = async_io_.DefineVariable<int>("TimeStep");
  }

  void attachOperator()
  {
    // NOTE: This is about I/O compression. This is currently diseable.
    // Alexis (01/10/2025): This could compress the I/O by a factor 2. However
    // this slow down the execution time by two when saving every 10 time step
    // iteration on my laptop. From 4.6 to 2.4 GB, and 23 to 33 secondes.
    // This is multithreaded.
    // Also it add blosc, libzstd and libblz4 as dependancies.
    // Don't forget to enable blosc2 and openmp in ADIOS2.cmake.
    //
    // pn_.AddOperation(compressor_op_, {{"compressor", "zstd"}, {"clevel",
    // "5"}, {"nthreads", "4"}}); receivers_.AddOperation(receiver_op_,
    // {{"compressor", "lz4"}, {"clevel", "3"}, {"nthreads", "4"}});
    // receivers_coords_.AddOperation(receiver_op_, {{"compressor", "lz4"},
    // {"clevel", "1"}, {"nthreads", "2"}});
  }

 public:
  SemIOController(const size_t nb_nodes, const size_t nb_iter, const size_t nb_receiver)
  {
    initAdios();
    configureIO();
    launchWriters();
    defineVariable(nb_nodes, nb_iter, nb_receiver);
    attachOperator();
  }

  SemIOController() = delete;

  ~SemIOController()
  {
    snaps_writer_.Close();
    receiver_writer_.Close();
  }

  void saveReceiver(vectorReal& receiver, const std::array<float, 3>& coords)
  {
    receiver_writer_.BeginStep();
    receiver_writer_.Put(receivers_, receiver.data());
    receiver_writer_.Put(receivers_coords_, coords.data());
    receiver_writer_.EndStep();
  }

  void saveSnapshot(const vectorReal& pnGlobal, const int timestep)
  {
    snaps_writer_.BeginStep();
    snaps_writer_.Put(timestep_, timestep);
    snaps_writer_.Put(pn_, pnGlobal.data());
    snaps_writer_.EndStep();
  }
};
