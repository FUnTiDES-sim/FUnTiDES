/**
 * @file sem_io_controller.h
 * @brief ADIOS2-based output of receiver traces and 3D field snapshots.
 */

#ifndef FUNTIDES_MAIN_FE_INCLUDE_SEM_IO_CONTROLLER_H_
#define FUNTIDES_MAIN_FE_INCLUDE_SEM_IO_CONTROLLER_H_

#include <adios2.h>
#include <data_type.h>

#include <cstddef>

#include "adios2/common/ADIOSTypes.h"
#include "adios2/cxx/Operator.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

/// Base name of the receiver trace output.
#define RECEIVERS_FILE "receivers"
/// Base name of the 3D snapshot output.
#define SNAPS_FILE "snapshots"

/**
 * @brief Writes receiver traces and 3D field snapshots to disk with ADIOS2.
 *
 * Receiver data goes through a synchronous IO, snapshots through an
 * asynchronous one. All files are opened by the constructor and closed by
 * the destructor. Snapshot data is float and is written as one chunk per MPI
 * rank into a global array.
 */
class SemIOController {
 private:
  adios2::ADIOS adios_;             ///< ADIOS2 instance (on MPI_COMM_WORLD when USE_MPI is defined).
  adios2::IO io_;                   ///< Synchronous IO, used for receivers.
  adios2::IO async_io_;             ///< Asynchronous IO, used for snapshots.
  adios2::Engine receiver_writer_;  ///< Writer engine for receivers.
  adios2::Engine snaps_writer_;     ///< Writer engine for snapshots.

  adios2::Variable<float> receivers_;         ///< Receiver traces, shape {nb_receiver, nb_iter}.
  adios2::Variable<float> receivers_coords_;  ///< Receiver coordinates, shape {nb_receiver, 3}.
  adios2::Variable<float> iter_times_;        ///< Iteration times, shape {nb_iter}.
  adios2::Variable<float> pn_;                ///< 3D field snapshot, global shape given at construction.
  adios2::Variable<int> timestep_;            ///< Time step index of the snapshot.

  adios2::Operator compressor_op_;  ///< Compression operator; not used yet.
  adios2::Operator receiver_op_;    ///< Receiver operator; not used yet.

  std::string rcv_file_{"rcv_not_set.bp"};    ///< Receiver output file name.
  std::string snap_file_{"snap_not_set.bp"};  ///< Snapshot output file name.

  /// Creates the ADIOS2 instance, on MPI_COMM_WORLD when USE_MPI is defined.
  void initAdios() {
#ifdef USE_MPI
    adios_ = adios2::ADIOS(MPI_COMM_WORLD);
#else
    adios_ = adios2::ADIOS();
#endif
  }

  /// Sets the output file names to RECEIVERS_FILE and SNAPS_FILE.
  void configureFilesName() {
    rcv_file_ = RECEIVERS_FILE;
    snap_file_ = SNAPS_FILE;
  }

  /// Declares the two IO objects (BP5 engine, 4 threads); the snapshot one writes asynchronously and profiles.
  void configureIO() {
    io_ = adios_.DeclareIO("AccousticSEMOutput");
    io_.SetEngine("BP5");
    io_.SetParameter("Threads", "4");

    async_io_ = adios_.DeclareIO("AsyncAccousticSEMOutput");
    async_io_.SetEngine("BP5");
    async_io_.SetParameter("AsyncWrite", "On");
    async_io_.SetParameter("Threads", "4");
    async_io_.SetParameter("Profile", "On");
    async_io_.SetParameter("ProfileUnits", "Microseconds");
  }

  /// Opens the receiver and snapshot files in write mode.
  void launchWriters() {
    receiver_writer_ = io_.Open(rcv_file_, adios2::Mode::Write);
    snaps_writer_ = async_io_.Open(snap_file_, adios2::Mode::Write);
  }

  /**
   * @brief Defines all ADIOS2 variables.
   * @param[in] global_dims Global shape of the snapshot array over all ranks.
   * @param[in] start_offsets Offset of this rank's chunk in the global array.
   * @param[in] local_dims Shape of this rank's chunk.
   * @param[in] nb_iter Number of time samples per receiver trace.
   * @param[in] nb_receiver Number of receivers.
   *
   * Receiver variables are declared with a zero offset and a local shape equal
   * to the global one.
   */
  void defineVariable(const std::vector<size_t>& global_dims, const std::vector<size_t>& start_offsets,
                      const std::vector<size_t>& local_dims, const size_t nb_iter, const size_t nb_receiver) {
    receivers_ = io_.DefineVariable<float>("AccousticReceiver", {nb_receiver, nb_iter}, {0, 0}, {nb_receiver, nb_iter});

    receivers_coords_ =
        io_.DefineVariable<float>("AccousticReceiverCoords", {nb_receiver, 3}, {0, 0}, {nb_receiver, 3});

    iter_times_ = io_.DefineVariable<float>("IterationTimes", {nb_iter}, {0}, {nb_iter});

    pn_ = async_io_.DefineVariable<float>("PressureField", global_dims, start_offsets, local_dims);

    timestep_ = async_io_.DefineVariable<int>("TimeStep", {1}, {0}, {1});
  }

  /// Placeholder for attaching compression operators; does nothing.
  void attachOperator() {}

 public:
  /**
   * @brief Creates the ADIOS2 pipeline and opens the receiver and snapshot files.
   * @param[in] global_dims Global shape of the snapshot array over all ranks.
   * @param[in] start_offsets Offset of this rank's chunk in the global array.
   * @param[in] local_dims Shape of this rank's chunk.
   * @param[in] nb_iter Number of time samples per receiver trace.
   * @param[in] nb_receiver Number of receivers.
   */
  SemIOController(const std::vector<size_t>& global_dims, const std::vector<size_t>& start_offsets,
                  const std::vector<size_t>& local_dims, const size_t nb_iter, const size_t nb_receiver) {
    initAdios();
    configureIO();
    configureFilesName();
    defineVariable(global_dims, start_offsets, local_dims, nb_iter, nb_receiver);
    launchWriters();
    attachOperator();
  }

  /// A controller cannot exist without its dimensions.
  SemIOController() = delete;

  /// Closes the snapshot and receiver files.
  ~SemIOController() {
    snaps_writer_.Close();
    receiver_writer_.Close();
  }

  /**
   * @brief Writes the receiver traces and coordinates as one ADIOS2 step.
   * @param[in] receiver Host array of the traces, expected to hold nb_receiver * nb_iter floats.
   * @param[in] coords Receiver coordinates in the order x, y, z.
   *
   * @todo VERIFY: what is the coordinate unit, and is the trace layout receiver-major (index r * nb_iter + it)?
   */
  void saveReceiver(vectorReal::host_mirror_type& receiver, const std::array<float, 3>& coords) {
    receiver_writer_.BeginStep();
    receiver_writer_.Put(receivers_, receiver.data());
    receiver_writer_.Put(receivers_coords_, coords.data());
    receiver_writer_.EndStep();
  }

  /**
   * @brief Writes one 3D field snapshot as one ADIOS2 step, asynchronously.
   * @param[in] pnGlobal Host array of this rank's chunk, of size local_dims[0] * local_dims[1] * local_dims[2].
   * @param[in] timestep Index of the time step, stored with the snapshot.
   */
  void saveSnapshot(const vectorReal::host_mirror_type& pnGlobal, const int timestep) {
    snaps_writer_.BeginStep();

    int ts_value[1] = {timestep};
    snaps_writer_.Put(timestep_, ts_value);
    snaps_writer_.Put(pn_, pnGlobal.data());
    snaps_writer_.EndStep();
  }
};

#endif  // FUNTIDES_MAIN_FE_INCLUDE_SEM_IO_CONTROLLER_H_
