#include <adios2.h>
#include <data_type.h>

#include <cstddef>

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
  adios2::Variable<float> iter_times_;
  adios2::Variable<float> pn_;
  adios2::Variable<int> timestep_;

  void initAdios() { adios_ = adios2::ADIOS(); }

  void configureIO()
  {
    io_ = adios_.DeclareIO("AccousticSEMOutput");
    io_.SetEngine("BP5");
    io_.SetParameter("Threads", "4");
    io_.SetParameter("CompressionMethod", "zlib");

    async_io_ = adios_.DeclareIO("AsyncAccousticSEMOutput");
    async_io_.SetEngine("BP5");
    async_io_.SetParameter("AsyncWrite", "On");
    async_io_.SetParameter("DirectIO", "On");
    async_io_.SetParameter("Threads", "4");
    async_io_.SetParameter("Profile", "On");
    async_io_.SetParameter("ProfileUnits", "Microseconds");
    async_io_.SetParameter("Verbose", "5");
  }

  void launchWriters()
  {
    receiver_writer_ = io_.Open(RECEIVERS_FILE, adios2::Mode::Write);
    snaps_writer_ = async_io_.Open(SNAPS_FILE, adios2::Mode::Write);
  }

 public:
  SemIOController(const size_t nb_nodes, const size_t nb_iter,
                  const size_t nb_receiver)
  {
    initAdios();
    configureIO();
    launchWriters();

    receivers_ =
        io_.DefineVariable<float>("AccousticReceiver", {nb_receiver, nb_iter},
                                  {0, 0}, {nb_receiver, nb_iter});
    iter_times_ =
        io_.DefineVariable<float>("IterationTimes", {nb_iter}, {0}, {nb_iter});
    pn_ = async_io_.DefineVariable<float>("PressureField", {nb_nodes}, {0},
                                          {nb_nodes});
    timestep_ = async_io_.DefineVariable<int>("TimeStep");
  }

  SemIOController() = default;
  ~SemIOController()
  {
    snaps_writer_.Close();
    receiver_writer_.Close();
  }

  void saveReceivers(vectorReal& receivers)
  {
    receiver_writer_.BeginStep();
    receiver_writer_.Put(receivers_, receivers.data());
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
