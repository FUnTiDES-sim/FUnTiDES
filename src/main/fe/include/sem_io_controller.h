#include <adios2.h>
#include <data_type.h>

#include <cstddef>

#define RECEIVERS_FILE "receivers.bp"

class SemIOController
{
 private:
  adios2::ADIOS adios_;
  adios2::IO io_;

  adios2::Variable<float> receivers_;
  adios2::Variable<float> iter_times_;

  void initAdios() { adios_ = adios2::ADIOS(); }

  void configureIO()
  {
    io_ = adios_.DeclareIO("AccousticSEMOutput");
    io_.SetEngine("BP4");
    io_.SetParameter("Threads", "4");
    io_.SetParameter("CompressionMethod", "zlib");
  }

 public:
  SemIOController(const size_t nb_iter, const size_t nb_receiver)
  {
    initAdios();
    configureIO();

    receivers_ =
        io_.DefineVariable<float>("AccousticReceiver", {nb_receiver, nb_iter},
                                  {0, 0}, {nb_receiver, nb_iter});

    iter_times_ =
        io_.DefineVariable<float>("IterationTimes", {nb_iter}, {0}, {nb_iter});
  }

  SemIOController() = default;
  ~SemIOController() = default;

  void saveReceivers(vectorReal& receivers)
  {
    using ValueType = typename vectorReal::value_type;

    auto writer = io_.Open(RECEIVERS_FILE, adios2::Mode::Write);

    writer.BeginStep();
    writer.Put(receivers_, receivers.data());
    writer.EndStep();
    writer.Close();
  }
};
