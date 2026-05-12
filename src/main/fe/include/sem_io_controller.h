/**
 * @file sem_io_controller.h
 * @brief Controls input/output operations for the SEM proxy using the ADIOS2 library.
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

/// Default output filename base for receiver traces.
#define RECEIVERS_FILE "receivers"
/// Default output filename base for 3D field snapshots.
#define SNAPS_FILE "snapshots"

/**
 * @class SemIOController
 * @brief Manages writing simulation data (receivers and snapshots) to disk using ADIOS2.
 *
 * This class encapsulates the ADIOS2 API to handle both synchronous (for traces) and
 * asynchronous (for large 3D snapshots) I/O operations. It manages file engines,
 * variables, and MPI I/O communicators.
 */
class SemIOController {
 private:
  adios2::ADIOS adios_;             ///< Main ADIOS2 instance.
  adios2::IO io_;                   ///< IO handler for synchronous operations (receivers).
  adios2::IO async_io_;             ///< IO handler for asynchronous operations (snapshots).
  adios2::Engine receiver_writer_;  ///< File writing engine for receivers.
  adios2::Engine snaps_writer_;     ///< File writing engine for snapshots.

  adios2::Variable<float> receivers_;         ///< ADIOS2 variable for receiver signal data.
  adios2::Variable<float> receivers_coords_;  ///< ADIOS2 variable for receiver 3D coordinates.
  adios2::Variable<float> iter_times_;        ///< ADIOS2 variable for iteration timestamps.
  adios2::Variable<float> pn_;                ///< ADIOS2 variable for 3D domain snapshots (e.g., pressure).
  adios2::Variable<int> timestep_;            ///< ADIOS2 variable for the current simulation timestep.

  adios2::Operator compressor_op_;  ///< ADIOS2 operator for data compression (e.g., zstd, blosc).
  adios2::Operator receiver_op_;    ///< ADIOS2 operator for receiver data processing.

  std::string rcv_file_{"rcv_not_set.bp"};    ///< Name of the receiver output file.
  std::string snap_file_{"snap_not_set.bp"};  ///< Name of the snapshot output file.

  /**
   * @brief Initializes the ADIOS2 instance, hooking into MPI if enabled.
   */
  void initAdios() {
#ifdef USE_MPI
    adios_ = adios2::ADIOS(MPI_COMM_WORLD);
#else
    adios_ = adios2::ADIOS();
#endif
  }

  /**
   * @brief Configures the base names for the output files.
   */
  void configureFilesName() {
    rcv_file_ = RECEIVERS_FILE;
    snap_file_ = SNAPS_FILE;
  }

  /**
   * @brief Declares the ADIOS2 IO interfaces and configures their engine parameters.
   */
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

  /**
   * @brief Opens the actual output files via the ADIOS2 engines.
   */
  void launchWriters() {
    receiver_writer_ = io_.Open(rcv_file_, adios2::Mode::Write);
    snaps_writer_ = async_io_.Open(snap_file_, adios2::Mode::Write);
  }

  /**
   * @brief Defines the data structure, shape, and offsets for all ADIOS2 variables.
   * * @param global_dims Global dimensions of the 3D grid across all MPI ranks.
   * @param start_offsets Global offset of the local 3D grid chunk.
   * @param local_dims Dimensions of the local 3D grid chunk for this MPI rank.
   * @param nb_iter Total number of simulation iterations (time samples).
   * @param nb_receiver Total number of receivers.
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

  /**
   * @brief Attaches data compression or transformation operators (currently a placeholder).
   */
  void attachOperator() {}

 public:
  /**
   * @brief Constructs the IO Controller and fully configures the ADIOS2 pipeline.
   * * @param global_dims Global dimensions of the 3D mesh.
   * @param start_offsets Local MPI chunk offset within the global mesh.
   * @param local_dims Local MPI chunk dimensions.
   * @param nb_iter Total number of time iterations to be saved.
   * @param nb_receiver Total number of receivers recording data.
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

  /// Deleted default constructor to enforce proper initialization parameters.
  SemIOController() = delete;

  /**
   * @brief Destroys the IO Controller and ensures all file streams are properly closed.
   */
  ~SemIOController() {
    snaps_writer_.Close();
    receiver_writer_.Close();
  }

  /**
   * @brief Synchronously writes receiver signal traces and coordinates to disk.
   * * @param receiver The host mirror view containing the receiver time series.
   * @param coords The 3D coordinates [X, Y, Z] of the receiver.
   */
  void saveReceiver(vectorReal::host_mirror_type& receiver, const std::array<float, 3>& coords) {
    receiver_writer_.BeginStep();
    receiver_writer_.Put(receivers_, receiver.data());
    receiver_writer_.Put(receivers_coords_, coords.data());
    receiver_writer_.EndStep();
  }

  /**
   * @brief Asynchronously dispatches a 3D field snapshot write operation to disk.
   * * @param pnGlobal The host mirror view of the 3D field array (e.g., pressure or displacement).
   * @param timestep The integer index of the current time step.
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
