#include "sem_proxy.h"

#include <boundary_synchronizer.h>
#include <cartesian_partitioner.h>
#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <source_and_receiver_utils.h>

#include <cmath>
#include <cxxopts.hpp>
#include <fstream>
#include <iomanip>

#ifdef COMPILE_DG
#include "dg_solver_data.h"
#endif
#include "rhs_acoustoelastic.h"
#ifdef USE_MPI
#include "mpi_backend.h"
#endif
#include "sem_solver.h"
#include "sem_solver_acoustoelastic.h"
#include "topology_factory.h"
#include "wavefield_acoustoelastic.h"

using namespace SourceAndReceiverUtils;
using namespace solver::fe;
using namespace solver::fe::solver_factory;
using namespace utils::enums;

SEMproxy::SEMproxy(const SemProxyOptions& opt) {
  int mpi_initialized = 0;
  InitMpi(&mpi_initialized);

  auto start_init = std::chrono::high_resolution_clock::now();

  InitSimParams(opt);
  InitMeshParams(opt);
  InitTopology();
  InitSync();
  InitTimeParams(opt);

  SetupSolver(opt);
  SetupAttenuation(opt);
  InitFiniteElem();
  SetupIo(opt);
  SetupDas(opt);

  time_init_ += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_init).count();
  DisplayInitMsg(opt);
}

void SEMproxy::SetupSolver(const SemProxyOptions& opt) {
  is_acousto_elastic_ = opt.isAcoustoElastic;
  const methodType method_type = GetMethod(opt.method);
  is_dg_ = (method_type == utils::enums::methodType::kDg);
  const implemType implem_type = GetImplem(opt.implem);
  const meshType mesh_type = GetMesh(opt.mesh);
  const modelLocationType model_location =
      opt.isModelOnNodes ? modelLocationType::kOnNodes : modelLocationType::kOnElements;
  const physicType physic_type = is_acousto_elastic_ ? physicType::kAcoustoElastic
                                                     : (opt.isElastic ? physicType::kElastic : physicType::kAcoustic);

  solver_ = createSolver(method_type, implem_type, mesh_type, model_location, physic_type, opt.order);

  const model::AnisotropyType anisotropy_type = GetAnisotropy(opt.anisotropy);
  use_sponge_ = opt.use_sponge;
  sponge_size_ = {opt.boundaries_size, opt.boundaries_size, opt.boundaries_size};
  surface_sponge_ = opt.surface_sponge;
  taper_delta_ = opt.taper_delta;

  if (opt.isElastic) {
    solver_->setAnisotropyType(anisotropy_type);
    if (anisotropy_type == model::AnisotropyType::kTTI && !opt.isModelOnNodes) {
      mesh_->initElasticityTensors(anisotropy_type);
    }
  }
}

void SEMproxy::SetupAttenuation(const SemProxyOptions& opt) {
  if (opt.qp > 0 || opt.qs > 0) {
    float qp = (opt.qp > 0) ? opt.qp : 1e9f;
    float qs = (opt.qs > 0) ? opt.qs : 1e9f;
    mesh_->setQualityFactors(qp, qs);
  }

  // Create temporary mirrors, fill on Host, deep_copy to Device
  auto ToView = [](const std::vector<float>& v, const std::string& name) {
    vectorReal d_view(name, v.size());
    auto h_view = Kokkos::create_mirror_view(d_view);
    for (size_t i = 0; i < v.size(); ++i) h_view(i) = v[i];
    Kokkos::deep_copy(d_view, h_view);
    return d_view;
  };

  if (!opt.sls_reference_angular_frequencies.empty()) {
    auto freq_view = ToView(opt.sls_reference_angular_frequencies, "slsFreqs");
    auto coeff_view = ToView(opt.sls_anelasticity_coefficients, "slsCoeffs");
    solver_->setSLSAttenuation(freq_view, coeff_view);
  } else if (opt.qp > 0 || opt.qs > 0) {
    float omega0 = 2.0f * M_PI * opt.f0;
    vectorReal freq_view("slsFreqAuto", 1);
    auto h_freq = Kokkos::create_mirror_view(freq_view);
    h_freq(0) = omega0;
    Kokkos::deep_copy(freq_view, h_freq);
    solver_->setSLSAttenuation(freq_view);
  }
}

void SEMproxy::SetupIo(const SemProxyOptions& opt) {
  size_t g_nx = static_cast<size_t>(opt.ex) * opt.order + 1;
  size_t g_ny = static_cast<size_t>(opt.ey) * opt.order + 1;
  size_t g_nz = static_cast<size_t>(opt.ez) * opt.order + 1;
  std::vector<size_t> global_dims = {g_nx, g_ny, g_nz};

  float dx = local_params_.global_lx / opt.ex;
  size_t elem_offset_x = static_cast<size_t>(std::round(local_params_.origin_x / dx));
  size_t node_offset_x = elem_offset_x * opt.order;
  std::vector<size_t> start_offsets = {node_offset_x, 0, 0};

  size_t l_nx = static_cast<size_t>(num_nodes_[0]);
  size_t l_ny = static_cast<size_t>(num_nodes_[1]);
  size_t l_nz = static_cast<size_t>(num_nodes_[2]);
  std::vector<size_t> local_dims = {l_nx, l_ny, l_nz};

  io_ctrl_ = std::make_shared<SemIOController>(global_dims, start_offsets, local_dims,
                                               static_cast<size_t>(num_samples_), static_cast<size_t>(1));

  is_snapshots_ = opt.snapshots;
  if (is_snapshots_) {
    snap_time_interval_ = opt.snap_time_interval;
  }
}

void SEMproxy::SetupDas(const SemProxyOptions& opt) {
  if (opt.das_type == "dipole") {
    das_type_ = SourceAndReceiverUtils::DASType::kDipole;
    das_num_samples_ = 2;
  } else if (opt.das_type == "strain") {
    das_type_ = SourceAndReceiverUtils::DASType::kStrainIntegration;
    das_num_samples_ = opt.das_samples;
  } else {
    das_type_ = SourceAndReceiverUtils::DASType::kNone;
  }

  if (das_type_ != SourceAndReceiverUtils::DASType::kNone) {
    if (!is_elastic_) throw std::runtime_error("DAS receivers require elastic simulation");
    das_gauge_length_ = opt.das_gauge_length;
    das_direction_ = SourceAndReceiverUtils::ComputeDASVector(opt.das_dip, opt.das_azimuth);
    das_vector_ = das_direction_;
    if (das_type_ == SourceAndReceiverUtils::DASType::kDipole) {
      for (int i = 0; i < 3; ++i) das_vector_[i] /= das_gauge_length_;
    }
  }
}

void SEMproxy::Run() {
  auto start_run_init = std::chrono::high_resolution_clock::now();
  solver_->setUseSponge(use_sponge_);
  solver_->computeFEInit(*mesh_, sponge_size_, surface_sponge_, taper_delta_);

  if (par_topology_.isDistributed()) {
    if (is_acousto_elastic_) {
      syncer_->synchronize(solver_->getMassMatrixAcoustic(), par_topology_);
      syncer_->synchronize(solver_->getMassMatrixElastic(), par_topology_);
    } else if (is_elastic_) {
      syncer_->synchronize(solver_->getMassMatrixElastic(), par_topology_);
    } else {
      syncer_->synchronize(solver_->getMassMatrixAcoustic(), par_topology_);
    }
    for (int c = 0; c < solver_->getNumComponents(); ++c) {
      syncer_->synchronize(solver_->getDampingMatrix(c), par_topology_);
    }
  }

  time_init_ += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_run_init).count();

  std::chrono::time_point<std::chrono::high_resolution_clock> start_compute_time, start_output_time;
  std::chrono::duration<double> total_compute_time(0), total_output_time(0);

#ifdef COMPILE_DG
  if (is_dg_) {
    DGWavefieldAcoustic wavefield(pn_dg_prev_, pn_dg_curr_);
    RhsAcoustic rhs(rhs_term_, rhs_element_, rhs_weights_);

    DGsolverDataAcoustic dgData(wavefield, rhs);

    for (int time_index = 0; time_index < num_samples_; time_index++) {
      start_compute_time = system_clock::now();
      solver_->computeOneStep(dt_, time_index, dgData);
      total_compute_time += system_clock::now() - start_compute_time;

      start_output_time = system_clock::now();

      auto d_field = dgData.getPreviousField(0);

      if (time_index % 50 == 0) {
        auto h_field = Kokkos::create_mirror_view(d_field);
        Kokkos::deep_copy(h_field, d_field);
        std::cout << "DG TimeStep=" << time_index << "  p(src_elem, dof=0)=" << h_field(source_element_, 0)
                  << "  p(rcv_elem, dof=0)=" << h_field(h_rhs_element_rcv_(0), 0) << std::endl;
      }

      if (is_snapshots_ && time_index % snap_time_interval_ == 0) {
        // Wait for previous to finish
        WaitSnapshots();
        Kokkos::fence();
        // Note: For optimal performance, h_field should be pre-allocated in InitArrays like the others
        auto h_field = Kokkos::create_mirror_view(d_field);
        Kokkos::deep_copy(h_field, d_field);

        auto io_task = [this, time_index, h_field]() {
          const int order = mesh_->getOrder();
          const int ex = num_elements_[0];
          const int ey = num_elements_[1];
          const int ez = num_elements_[2];
          const int zElem = ez / 2;
          const int n1d = order + 1;
          const int icZ = order / 2;
          std::ostringstream fname;
          fname << "slice_dg_" << std::setfill('0') << std::setw(5) << time_index << ".dat";
          std::ofstream fslice(fname.str());
          for (int ej_idx = 0; ej_idx < ey; ++ej_idx) {
            for (int ib = 0; ib < n1d; ++ib) {
              bool first = true;
              for (int ei_idx = 0; ei_idx < ex; ++ei_idx) {
                int const elem = ei_idx + ej_idx * ex + zElem * ex * ey;
                for (int ia = 0; ia < n1d; ++ia) {
                  int const dof = ia + ib * n1d + icZ * n1d * n1d;
                  if (!first) fslice << " ";
                  fslice << h_field(elem, dof);
                  first = false;
                }
              }
              fslice << "\n";
            }
          }
          fslice.close();
        };

        if (d_field.data() != h_field.data()) {
          snapshot_futures_.push_back(std::async(std::launch::async, io_task));
        } else {
          io_task();  // CPU-only fallback
        }
      }

      auto h_field = Kokkos::create_mirror_view(d_field);
      Kokkos::deep_copy(h_field, d_field);
      const int order = mesh_->getOrder();
      float varnp1 = 0.0f;
      for (int i = 0; i <= order; i++) {
        for (int j = 0; j <= order; j++) {
          for (int k = 0; k <= order; k++) {
            int localDof = i + j * (order + 1) + k * (order + 1) * (order + 1);
            varnp1 += h_field(h_rhs_element_rcv_(0), localDof) * h_rhs_weights_rcv_(0, localDof);
          }
        }
      }
      h_pn_at_receiver_(0, time_index) = varnp1;
      dgData.swapWavefields();
      total_output_time += system_clock::now() - start_output_time;
    }

    std::ofstream fout("receiver_trace.txt");
    fout << "# time pressure_at_receiver\n";
    for (int t = 0; t < num_samples_; ++t) fout << t * dt_ << " " << h_pn_at_receiver_(0, t) << "\n";
    fout.close();
  }
#endif  // COMPILE_DG
  if (is_acousto_elastic_) {
    WavefieldAcoustoElastic wavefield(pn_global_prev_, pn_global_curr_, uxn_global_prev_, uxn_global_curr_,
                                      uyn_global_prev_, uyn_global_curr_, uzn_global_prev_, uzn_global_curr_);
    RhsAcoustoElastic rhs(rhs_term_, rhs_element_, rhs_weights_, rhs_term_x_, rhs_term_y_, rhs_term_z_);
    SEMsolverDataAcoustoElastic solver_data(wavefield, rhs);

    for (int time_index = 0; time_index < num_samples_; time_index++) {
      start_compute_time = std::chrono::high_resolution_clock::now();
      solver_->computeOneStep(dt_, time_index, solver_data);
      total_compute_time += std::chrono::high_resolution_clock::now() - start_compute_time;

      start_output_time = std::chrono::high_resolution_clock::now();
      if (time_index % 50 == 0) {
        solver_->outputSolutionValues(time_index, h_rhs_element_(0), pn_global_prev_, "pnGlobal");
        solver_->outputSolutionValues(time_index, h_rhs_element_rcv_(0), uxn_global_prev_, "uxnGlobal");
        solver_->outputSolutionValues(time_index, h_rhs_element_rcv_(0), uyn_global_prev_, "uynGlobal");
        solver_->outputSolutionValues(time_index, h_rhs_element_rcv_(0), uzn_global_prev_, "uznGlobal");
      }

      if (is_snapshots_ && time_index % snap_time_interval_ == 0) {
        WaitSnapshots();
        Kokkos::deep_copy(h_pn_global_prev_, pn_global_prev_);
        Kokkos::deep_copy(h_uzn_global_prev_, uzn_global_prev_);

        auto io_task = [this, time_index, h_pn = h_pn_global_prev_, h_uz = h_uzn_global_prev_]() {
          io_ctrl_->saveSnapshot(h_pn, time_index);
          io_ctrl_->saveSnapshot(h_uz, time_index);
        };

        // Fallback checks if GPU data == CPU data (CPU-only build)
        if (pn_global_prev_.data() != h_pn_global_prev_.data()) {
          snapshot_futures_.push_back(std::async(std::launch::async, io_task));
        } else {
          io_task();
        }
      }

      Kokkos::deep_copy(h_pn_global_curr_, pn_global_curr_);
      const int order = mesh_->getOrder();
      float var_np1 = 0.0;
      for (int i = 0; i < order + 1; i++) {
        for (int j = 0; j < order + 1; j++) {
          for (int k = 0; k < order + 1; k++) {
            int node_idx = mesh_->globalNodeIndex(h_rhs_element_rcv_(0), i, j, k);
            int global_node_on_elem = i + j * (order + 1) + k * (order + 1) * (order + 1);
            var_np1 += h_pn_global_curr_(node_idx) * h_rhs_weights_rcv_(0, global_node_on_elem);
          }
        }
      }
      h_pn_at_receiver_(0, time_index) = var_np1;
      solver_data.swapWavefields();
      total_output_time += std::chrono::high_resolution_clock::now() - start_output_time;
    }

    start_output_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < h_pn_at_receiver_.extent(0); i++) {
      auto subview = Kokkos::subview(h_pn_at_receiver_, i, Kokkos::ALL());
      vectorReal::host_mirror_type subset("receiver_save", num_samples_);
      for (int j = 0; j < num_samples_; ++j) subset(j) = subview(j);
      io_ctrl_->saveReceiver(subset, src_coord_);
    }
    total_output_time += std::chrono::high_resolution_clock::now() - start_output_time;

  } else if (!is_elastic_) {
    WavefieldAcoustic wavefield(pn_global_prev_, pn_global_curr_);
    RhsAcoustic rhs(rhs_term_, rhs_element_, rhs_weights_);
    SEMsolverDataAcoustic solver_data(wavefield, rhs);

    for (int time_index = 0; time_index < num_samples_; time_index++) {
      start_compute_time = std::chrono::high_resolution_clock::now();
      solver_->computeForces(dt_, time_index, solver_data);
      if (par_topology_.isDistributed()) {
        for (int c = 0; c < solver_->getNumComponents(); ++c) {
          syncer_->synchronize(solver_->getForceVector(c), par_topology_);
        }
      }
      solver_->updateSolutionForward(dt_, solver_data);
      total_compute_time += std::chrono::high_resolution_clock::now() - start_compute_time;

      start_output_time = std::chrono::high_resolution_clock::now();

      if (time_index % 50 == 0) {
        solver_->outputSolutionValues(time_index, h_rhs_element_(0), pn_global_prev_, "pnGlobal");
      }

      if (is_snapshots_ && time_index % snap_time_interval_ == 0) {
#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        WaitSnapshots();
        Kokkos::deep_copy(h_pn_global_prev_, pn_global_prev_);

        // Offload both ADIOS2 and Slice formatting to background thread
        auto io_task = [this, time_index, h_field = h_pn_global_prev_]() {
          io_ctrl_->saveSnapshot(h_field, time_index);

          int nx = num_nodes_[0];
          int ny = num_nodes_[1];
          int nz = num_nodes_[2];
          int z_slice = nz / 2;
          std::ostringstream fname;
          fname << "slice_" << std::setfill('0') << std::setw(5) << time_index << ".dat";
          std::ofstream fslice(fname.str());
          for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
              int node_idx = ix + iy * nx + z_slice * nx * ny;
              fslice << h_field(node_idx);
              if (ix < nx - 1) fslice << " ";
            }
            fslice << "\n";
          }
          fslice.close();
        };

        if (pn_global_prev_.data() != h_pn_global_prev_.data()) {
          snapshot_futures_.push_back(std::async(std::launch::async, io_task));
        } else {
          io_task();
        }
      }

      // Bring wavefield array back to Host to compute Receiver Value
      Kokkos::deep_copy(h_pn_global_curr_, pn_global_curr_);
      const int order = mesh_->getOrder();
      float var_np1 = 0.0;
      for (int i = 0; i < order + 1; i++) {
        for (int j = 0; j < order + 1; j++) {
          for (int k = 0; k < order + 1; k++) {
            int node_idx = mesh_->globalNodeIndex(h_rhs_element_rcv_(0), i, j, k);
            int global_node_on_elem = i + j * (order + 1) + k * (order + 1) * (order + 1);
            var_np1 += h_pn_global_curr_(node_idx) * h_rhs_weights_rcv_(0, global_node_on_elem);
          }
        }
      }
      h_pn_at_receiver_(0, time_index) = var_np1;

      solver_data.swapWavefields();
      total_output_time += std::chrono::high_resolution_clock::now() - start_output_time;
    }

    start_output_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < h_pn_at_receiver_.extent(0); i++) {
      auto subview = Kokkos::subview(h_pn_at_receiver_, i, Kokkos::ALL());
      vectorReal::host_mirror_type subset("receiver_save", num_samples_);
      for (int j = 0; j < num_samples_; ++j) subset(j) = subview(j);
      io_ctrl_->saveReceiver(subset, src_coord_);
    }

    std::ofstream fout("receiver_trace.txt");
    fout << "# time pressure_at_receiver\n";
    for (int t = 0; t < num_samples_; ++t) fout << t * dt_ << " " << h_pn_at_receiver_(0, t) << "\n";
    fout.close();
    total_output_time += std::chrono::high_resolution_clock::now() - start_output_time;

  } else {  // acoustic
    WavefieldElastic wavefield(uxn_global_prev_, uxn_global_curr_, uyn_global_prev_, uyn_global_curr_, uzn_global_prev_,
                               uzn_global_curr_);
    RhsElastic rhs(rhs_term_x_, rhs_term_y_, rhs_term_z_, rhs_element_, rhs_weights_);
    SEMsolverDataElastic solver_data(wavefield, rhs);

    for (int time_index = 0; time_index < num_samples_; time_index++) {
      start_compute_time = std::chrono::high_resolution_clock::now();
      solver_->computeForces(dt_, time_index, solver_data);
      if (par_topology_.isDistributed()) {
        for (int c = 0; c < solver_->getNumComponents(); ++c)
          syncer_->synchronize(solver_->getForceVector(c), par_topology_);
      }
      solver_->updateSolutionForward(dt_, solver_data);
      total_compute_time += std::chrono::high_resolution_clock::now() - start_compute_time;

      start_output_time = std::chrono::high_resolution_clock::now();

      if (time_index % 50 == 0) {
        solver_->outputSolutionValues(time_index, h_rhs_element_(0), uxn_global_prev_, "uxnGlobal");
        solver_->outputSolutionValues(time_index, h_rhs_element_(0), uyn_global_prev_, "uynGlobal");
        solver_->outputSolutionValues(time_index, h_rhs_element_(0), uzn_global_prev_, "uznGlobal");
      }

      if (is_snapshots_ && time_index % snap_time_interval_ == 0) {
        WaitSnapshots();
        Kokkos::deep_copy(h_uxn_global_prev_, uxn_global_prev_);

        auto io_task = [this, time_index, h_ux = h_uxn_global_prev_]() { io_ctrl_->saveSnapshot(h_ux, time_index); };

        if (uxn_global_prev_.data() != h_uxn_global_prev_.data()) {
          snapshot_futures_.push_back(std::async(std::launch::async, io_task));
        } else {
          io_task();
        }
      }

      Kokkos::deep_copy(h_uxn_global_curr_, uxn_global_curr_);
      Kokkos::deep_copy(h_uyn_global_curr_, uyn_global_curr_);
      Kokkos::deep_copy(h_uzn_global_curr_, uzn_global_curr_);

      const int order = mesh_->getOrder();
      float var_ux_np1 = 0.0;
      float var_uy_np1 = 0.0;
      float var_uz_np1 = 0.0;
      for (int i = 0; i < order + 1; i++) {
        for (int j = 0; j < order + 1; j++) {
          for (int k = 0; k < order + 1; k++) {
            int node_idx = mesh_->globalNodeIndex(h_rhs_element_rcv_(0), i, j, k);
            int global_node_on_elem = i + j * (order + 1) + k * (order + 1) * (order + 1);
            var_ux_np1 += h_uxn_global_curr_(node_idx) * h_rhs_weights_rcv_(0, global_node_on_elem);
            var_uy_np1 += h_uyn_global_curr_(node_idx) * h_rhs_weights_rcv_(0, global_node_on_elem);
            var_uz_np1 += h_uzn_global_curr_(node_idx) * h_rhs_weights_rcv_(0, global_node_on_elem);
          }
        }
      }

      h_uxn_at_receiver_(0, time_index) = var_ux_np1;
      h_uyn_at_receiver_(0, time_index) = var_uy_np1;
      h_uzn_at_receiver_(0, time_index) = var_uz_np1;

      if (das_type_ != SourceAndReceiverUtils::DASType::kNone) {
        float das_val = 0.0f;
        const int total_das_nodes = static_cast<int>(das_node_ids_.size());
        for (int i_node = 0; i_node < total_das_nodes; ++i_node) {
          int n_id = das_node_ids_[i_node];
          if (n_id >= 0) {
            float w = das_weights_[i_node];
            das_val += (h_uxn_global_curr_(n_id) * das_vector_[0] + h_uyn_global_curr_(n_id) * das_vector_[1] +
                        h_uzn_global_curr_(n_id) * das_vector_[2]) *
                       w;
          }
        }
        h_das_signal_(time_index) = das_val;
      }

      solver_data.swapWavefields();
      total_output_time += std::chrono::high_resolution_clock::now() - start_output_time;
    }

    start_output_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < h_uxn_at_receiver_.extent(0); i++) {
      auto subview = Kokkos::subview(h_uxn_at_receiver_, i, Kokkos::ALL());
      vectorReal::host_mirror_type subset("receiver_save", num_samples_);
      for (int j = 0; j < num_samples_; ++j) subset(j) = subview(j);
      io_ctrl_->saveReceiver(subset, src_coord_);
    }

    if (das_type_ != SourceAndReceiverUtils::DASType::kNone) {
      std::ofstream f_das("das_trace.txt");
      f_das << "# time das_signal\n";
      for (int t = 0; t < num_samples_; ++t) f_das << t * dt_ << " " << h_das_signal_(t) << "\n";
      f_das.close();
    }
    total_output_time += std::chrono::high_resolution_clock::now() - start_output_time;
  }

  WaitSnapshots();
  time_compute_ = total_compute_time.count();
  time_io_ = total_output_time.count();
  DisplayPerfMsg();
}

void SEMproxy::InitArrays() {
  const auto n_nodes = mesh_->getNumberOfNodes();
  const auto n_elements = mesh_->getNumberOfElements();
  const auto n_pts_per_elem = mesh_->getNumberOfPointsPerElement();

  // Create Device Views
  rhs_element_ = vectorInt("rhsElement", num_rhs_);
  rhs_weights_ = arrayReal("RHSWeight", num_rhs_, n_pts_per_elem);
  rhs_element_rcv_ = vectorInt("rhsElementRcv", 1);
  rhs_weights_rcv_ = arrayReal("RHSWeightRcv", 1, n_pts_per_elem);

  if (is_acousto_elastic_) {
    rhs_term_ = allocateArray2D<arrayReal>(num_rhs_, num_samples_, "RHSTerm");
    rhs_term_x_ = allocateArray2D<arrayReal>(num_rhs_, num_samples_, "RHSTermx");
    rhs_term_y_ = allocateArray2D<arrayReal>(num_rhs_, num_samples_, "RHSTermy");
    rhs_term_z_ = allocateArray2D<arrayReal>(num_rhs_, num_samples_, "RHSTermz");
    pn_global_curr_ = allocateVector<vectorReal>(n_nodes, "pnGlobalCurr");
    pn_global_prev_ = allocateVector<vectorReal>(n_nodes, "pnGlobalPrev");
    pn_at_receiver_ = allocateArray2D<arrayReal>(1, num_samples_, "pn_at_receiver_");
    uxn_global_curr_ = allocateVector<vectorReal>(n_nodes, "uxnGlobalCurr");
    uyn_global_curr_ = allocateVector<vectorReal>(n_nodes, "uynGlobalCurr");
    uzn_global_curr_ = allocateVector<vectorReal>(n_nodes, "uznGlobalCurr");
    uxn_global_prev_ = allocateVector<vectorReal>(n_nodes, "uxnGlobalPrev");
    uyn_global_prev_ = allocateVector<vectorReal>(n_nodes, "uynGlobalPrev");
    uzn_global_prev_ = allocateVector<vectorReal>(n_nodes, "uznGlobalPrev");
  } else if (is_dg_) {
    rhs_term_ = allocateArray2D<arrayReal>(num_rhs_, num_samples_, "RHSTerm");
    pn_dg_prev_ = allocateArray2D<arrayReal>(n_elements, n_pts_per_elem, "pnDGPrev");
    pn_dg_curr_ = allocateArray2D<arrayReal>(n_elements, n_pts_per_elem, "pnDGCurr");
    pn_at_receiver_ = allocateArray2D<arrayReal>(1, num_samples_, "pn_at_receiver_");
  } else if (!is_elastic_) {
    rhs_term_ = allocateArray2D<arrayReal>(num_rhs_, num_samples_, "RHSTerm");
    pn_global_curr_ = allocateVector<vectorReal>(n_nodes, "pnGlobalCurr");
    pn_global_prev_ = allocateVector<vectorReal>(n_nodes, "pnGlobalPrev");
    pn_at_receiver_ = allocateArray2D<arrayReal>(1, num_samples_, "pn_at_receiver_");
  } else {
    rhs_term_x_ = arrayReal("RHSTermx", num_rhs_, num_samples_);
    rhs_term_y_ = arrayReal("RHSTermy", num_rhs_, num_samples_);
    rhs_term_z_ = arrayReal("RHSTermz", num_rhs_, num_samples_);
    uxn_global_curr_ = vectorReal("uxnGlobalCurr", n_nodes);
    uyn_global_curr_ = vectorReal("uynGlobalCurr", n_nodes);
    uzn_global_curr_ = vectorReal("uznGlobalCurr", n_nodes);
    uxn_global_prev_ = vectorReal("uxnGlobalPrev", n_nodes);
    uyn_global_prev_ = vectorReal("uynGlobalPrev", n_nodes);
    uzn_global_prev_ = vectorReal("uznGlobalPrev", n_nodes);
    uxn_at_receiver_ = arrayReal("uxnAtReceiver", 1, num_samples_);
    uyn_at_receiver_ = arrayReal("uynAtReceiver ", 1, num_samples_);
    uzn_at_receiver_ = arrayReal("uznAtReceiver", 1, num_samples_);
  }

  if (das_type_ != SourceAndReceiverUtils::DASType::kNone) {
    das_signal_ = vectorReal("dasSignal", num_samples_);
  }

  // Create Host Mirrors
  h_rhs_element_ = Kokkos::create_mirror_view(rhs_element_);
  h_rhs_element_rcv_ = Kokkos::create_mirror_view(rhs_element_rcv_);
  h_rhs_weights_ = Kokkos::create_mirror_view(rhs_weights_);
  h_rhs_weights_rcv_ = Kokkos::create_mirror_view(rhs_weights_rcv_);
  h_pn_global_curr_ = Kokkos::create_mirror_view(pn_global_curr_);
  h_pn_global_prev_ = Kokkos::create_mirror_view(pn_global_prev_);
  h_uxn_global_curr_ = Kokkos::create_mirror_view(uxn_global_curr_);
  h_uyn_global_curr_ = Kokkos::create_mirror_view(uyn_global_curr_);
  h_uzn_global_curr_ = Kokkos::create_mirror_view(uzn_global_curr_);
  h_uxn_global_prev_ = Kokkos::create_mirror_view(uxn_global_prev_);
  h_uyn_global_prev_ = Kokkos::create_mirror_view(uyn_global_prev_);
  h_uzn_global_prev_ = Kokkos::create_mirror_view(uzn_global_prev_);

  if (is_acousto_elastic_) {
    h_rhs_term_ = Kokkos::create_mirror_view(rhs_term_);
    h_rhs_term_x_ = Kokkos::create_mirror_view(rhs_term_x_);
    h_rhs_term_y_ = Kokkos::create_mirror_view(rhs_term_y_);
    h_rhs_term_z_ = Kokkos::create_mirror_view(rhs_term_z_);
    h_pn_at_receiver_ = Kokkos::create_mirror_view(pn_at_receiver_);
  } else if (!is_elastic_) {
    h_rhs_term_ = Kokkos::create_mirror_view(rhs_term_);
    h_pn_at_receiver_ = Kokkos::create_mirror_view(pn_at_receiver_);
  } else if (is_dg_) {
    h_rhs_term_ = Kokkos::create_mirror_view(rhs_term_);
    h_pn_at_receiver_ = Kokkos::create_mirror_view(pn_at_receiver_);
  } else {
    h_rhs_term_x_ = Kokkos::create_mirror_view(rhs_term_x_);
    h_rhs_term_y_ = Kokkos::create_mirror_view(rhs_term_y_);
    h_rhs_term_z_ = Kokkos::create_mirror_view(rhs_term_z_);
    h_uxn_at_receiver_ = Kokkos::create_mirror_view(uxn_at_receiver_);
    h_uyn_at_receiver_ = Kokkos::create_mirror_view(uyn_at_receiver_);
    h_uzn_at_receiver_ = Kokkos::create_mirror_view(uzn_at_receiver_);
  }

  if (das_type_ != SourceAndReceiverUtils::DASType::kNone) {
    h_das_signal_ = Kokkos::create_mirror_view(das_signal_);
  }
}

void SEMproxy::InitSource() {
  int ex = num_elements_[0];
  int ey = num_elements_[1];
  int ez = num_elements_[2];
  int lx = domain_size_[0];
  int ly = domain_size_[1];
  int lz = domain_size_[2];

  float rel_src_x = src_coord_[0] - local_params_.origin_x;
  float rel_src_y = src_coord_[1] - local_params_.origin_y;
  float rel_src_z = src_coord_[2] - local_params_.origin_z;

  bool source_on_rank = (rel_src_x >= 0 && rel_src_x < lx);
  int src_index = source_on_rank ? floor((rel_src_x * ex) / lx) + floor((rel_src_y * ey) / ly) * ex +
                                       floor((rel_src_z * ez) / lz) * ey * ex
                                 : 0;

  for (int i = 0; i < 1; i++) h_rhs_element_(i) = src_index;

  float corner_coords[8][3];
  int corner_iter = 0;
  int nodes_corner[2] = {0, mesh_->getOrder()};

  for (int k : nodes_corner) {
    for (int j : nodes_corner) {
      for (int i : nodes_corner) {
        int node_idx = mesh_->globalNodeIndex(h_rhs_element_(0), i, j, k);
        corner_coords[corner_iter][0] = mesh_->nodeCoord(node_idx, 0);
        corner_coords[corner_iter][2] = mesh_->nodeCoord(node_idx, 2);
        corner_coords[corner_iter][1] = mesh_->nodeCoord(node_idx, 1);
        corner_iter++;
      }
    }
  }

  std::vector<float> source_term = utils_.computeSourceTerm(num_samples_, dt_, f0_, ricker_order_, t_peak_);

  if (is_acousto_elastic_) {
    bool const source_in_fluid = (src_coord_[2] >= local_params_.acoustoElasticBoundaryZ);
    if (source_in_fluid) {
      for (int j = 0; j < num_samples_; j++) h_rhs_term_(0, j) = source_term[j];
    } else {
      for (int j = 0; j < num_samples_; j++) {
        h_rhs_term_x_(0, j) = source_term[j];
        h_rhs_term_y_(0, j) = source_term[j];
        h_rhs_term_z_(0, j) = source_term[j];
      }
    }
  } else if (!is_elastic_) {
    for (int j = 0; j < num_samples_; j++) h_rhs_term_(0, j) = source_term[j];
  } else {
    for (int j = 0; j < num_samples_; j++) {
      h_rhs_term_x_(0, j) = source_term[j];
      h_rhs_term_y_(0, j) = source_term[j];
      h_rhs_term_z_(0, j) = source_term[j];
    }
  }

  source_element_ = h_rhs_element_(0);
  int order = mesh_->getOrder();

  if (source_on_rank) {
    switch (order) {
      case 1:
        SourceAndReceiverUtils::ComputeRHSWeights<1>(corner_coords, src_coord_, h_rhs_weights_);
        break;
      case 2:
        SourceAndReceiverUtils::ComputeRHSWeights<2>(corner_coords, src_coord_, h_rhs_weights_);
        break;
      case 3:
        SourceAndReceiverUtils::ComputeRHSWeights<3>(corner_coords, src_coord_, h_rhs_weights_);
        break;
      case 4:
        SourceAndReceiverUtils::ComputeRHSWeights<4>(corner_coords, src_coord_, h_rhs_weights_);
        break;
      case 5:
        SourceAndReceiverUtils::ComputeRHSWeights<5>(corner_coords, src_coord_, h_rhs_weights_);
        break;
      default:
        throw std::runtime_error("Unsupported order: " + std::to_string(order));
    }
  } else {
    for (int k = 0; k < mesh_->getNumberOfPointsPerElement(); ++k) h_rhs_weights_(0, k) = 0.0f;
  }

  float rel_rcv_x = rcv_coord_[0] - local_params_.origin_x;
  int rcv_index =
      floor((rel_rcv_x * ex) / lx) + floor((rcv_coord_[1] * ey) / ly) * ex + floor((rcv_coord_[2] * ez) / lz) * ey * ex;

  if (rcv_index < 0 || rcv_index >= mesh_->getNumberOfElements()) rcv_index = 0;
  for (int i = 0; i < 1; i++) h_rhs_element_rcv_(i) = rcv_index;

  float corner_coords_rcv[8][3];
  corner_iter = 0;
  for (int k : nodes_corner) {
    for (int j : nodes_corner) {
      for (int i : nodes_corner) {
        int node_idx = mesh_->globalNodeIndex(h_rhs_element_rcv_(0), i, j, k);
        corner_coords_rcv[corner_iter][0] = mesh_->nodeCoord(node_idx, 0);
        corner_coords_rcv[corner_iter][2] = mesh_->nodeCoord(node_idx, 2);
        corner_coords_rcv[corner_iter][1] = mesh_->nodeCoord(node_idx, 1);
        corner_iter++;
      }
    }
  }

  switch (order) {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(corner_coords_rcv, rcv_coord_, h_rhs_weights_rcv_);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(corner_coords_rcv, rcv_coord_, h_rhs_weights_rcv_);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(corner_coords_rcv, rcv_coord_, h_rhs_weights_rcv_);
      break;
    case 4:
      SourceAndReceiverUtils::ComputeRHSWeights<4>(corner_coords_rcv, rcv_coord_, h_rhs_weights_rcv_);
      break;
    case 5:
      SourceAndReceiverUtils::ComputeRHSWeights<5>(corner_coords_rcv, rcv_coord_, h_rhs_weights_rcv_);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }

  if (das_type_ != SourceAndReceiverUtils::DASType::kNone) {
    const int npe = mesh_->getNumberOfPointsPerElement();
    const int total_das_nodes = das_num_samples_ * npe;
    das_node_ids_.resize(total_das_nodes, -1);
    das_weights_.resize(total_das_nodes, 0.0f);
    std::vector<float> sample_locs(das_num_samples_);
    if (das_num_samples_ == 1) {
      sample_locs[0] = 0.0f;
    } else {
      for (int i = 0; i < das_num_samples_; ++i)
        sample_locs[i] = -0.5f + static_cast<float>(i) / (das_num_samples_ - 1);
    }
    std::vector<float> integration_consts(das_num_samples_);
    if (das_type_ == SourceAndReceiverUtils::DASType::kDipole) {
      integration_consts[0] = -1.0f;
      integration_consts[1] = 1.0f;
    } else {
      for (int i = 0; i < das_num_samples_; ++i) integration_consts[i] = 1.0f / das_num_samples_;
    }

    for (int i_sample = 0; i_sample < das_num_samples_; ++i_sample) {
      std::array<float, 3> sample_coord;
      for (int d = 0; d < 3; ++d)
        sample_coord[d] = rcv_coord_[d] + das_direction_[d] * das_gauge_length_ * sample_locs[i_sample];
      float rel_x = sample_coord[0] - local_params_.origin_x;
      int sample_elem_idx = static_cast<int>(std::floor((rel_x * ex) / lx)) +
                            static_cast<int>(std::floor((sample_coord[1] * ey) / ly)) * ex +
                            static_cast<int>(std::floor((sample_coord[2] * ez) / lz)) * ey * ex;
      if (sample_elem_idx < 0 || sample_elem_idx >= mesh_->getNumberOfElements()) sample_elem_idx = 0;

      float sample_corners[8][3];
      int ci = 0;
      for (int kk : nodes_corner) {
        for (int jj : nodes_corner) {
          for (int ii : nodes_corner) {
            int n_idx = mesh_->globalNodeIndex(sample_elem_idx, ii, jj, kk);
            sample_corners[ci][0] = mesh_->nodeCoord(n_idx, 0);
            sample_corners[ci][1] = mesh_->nodeCoord(n_idx, 1);
            sample_corners[ci][2] = mesh_->nodeCoord(n_idx, 2);
            ci++;
          }
        }
      }

      int base_idx = i_sample * npe;
      for (int kk = 0; kk < order + 1; ++kk) {
        for (int jj = 0; jj < order + 1; ++jj) {
          for (int ii = 0; ii < order + 1; ++ii) {
            int local_node = ii + jj * (order + 1) + kk * (order + 1) * (order + 1);
            das_node_ids_[base_idx + local_node] = mesh_->globalNodeIndex(sample_elem_idx, ii, jj, kk);
          }
        }
      }

      switch (order) {
        case 1:
          SourceAndReceiverUtils::ComputeDASWeightsForSample<1>(sample_corners, sample_coord, das_direction_,
                                                                integration_consts[i_sample], das_type_,
                                                                &das_weights_[base_idx]);
          break;
        case 2:
          SourceAndReceiverUtils::ComputeDASWeightsForSample<2>(sample_corners, sample_coord, das_direction_,
                                                                integration_consts[i_sample], das_type_,
                                                                &das_weights_[base_idx]);
          break;
        case 3:
          SourceAndReceiverUtils::ComputeDASWeightsForSample<3>(sample_corners, sample_coord, das_direction_,
                                                                integration_consts[i_sample], das_type_,
                                                                &das_weights_[base_idx]);
          break;
        case 4:
          SourceAndReceiverUtils::ComputeDASWeightsForSample<4>(sample_corners, sample_coord, das_direction_,
                                                                integration_consts[i_sample], das_type_,
                                                                &das_weights_[base_idx]);
          break;
        case 5:
          SourceAndReceiverUtils::ComputeDASWeightsForSample<5>(sample_corners, sample_coord, das_direction_,
                                                                integration_consts[i_sample], das_type_,
                                                                &das_weights_[base_idx]);
          break;
        default:
          throw std::runtime_error("Unsupported order for DAS: " + std::to_string(order));
      }
    }
  }

  // --- DMA Copy everything to the GPU ---
  Kokkos::deep_copy(rhs_element_, h_rhs_element_);
  Kokkos::deep_copy(rhs_element_rcv_, h_rhs_element_rcv_);
  Kokkos::deep_copy(rhs_weights_, h_rhs_weights_);
  Kokkos::deep_copy(rhs_weights_rcv_, h_rhs_weights_rcv_);

  if (is_acousto_elastic_) {
    if (src_coord_[2] >= local_params_.acoustoElasticBoundaryZ)
      Kokkos::deep_copy(rhs_term_, h_rhs_term_);
    else {
      Kokkos::deep_copy(rhs_term_x_, h_rhs_term_x_);
      Kokkos::deep_copy(rhs_term_y_, h_rhs_term_y_);
      Kokkos::deep_copy(rhs_term_z_, h_rhs_term_z_);
    }
  } else if (!is_elastic_) {
    Kokkos::deep_copy(rhs_term_, h_rhs_term_);
  } else {
    Kokkos::deep_copy(rhs_term_x_, h_rhs_term_x_);
    Kokkos::deep_copy(rhs_term_y_, h_rhs_term_y_);
    Kokkos::deep_copy(rhs_term_z_, h_rhs_term_z_);
  }
}

void SEMproxy::SaveSnapshot(int timestep, const vectorReal& d_data, vectorReal::host_mirror_type& h_data) const {
  // Ultra-fast deep copy into the PRE-ALLOCATED host buffer. No malloc!
  Kokkos::deep_copy(h_data, d_data);

  // Hand directly to ADIOS2
  io_ctrl_->saveSnapshot(h_data, timestep);
}

implemType SEMproxy::GetImplem(std::string implem_arg) {
  if (implem_arg == "makutu") return implemType::kMakutu;
  throw std::invalid_argument("Implementation type invalid. Only 'makutu' is supported.");
};

meshType SEMproxy::GetMesh(std::string mesh_arg) {
  if (mesh_arg == "cartesian") return meshType::kStruct;
  if (mesh_arg == "ucartesian") return meshType::kUnstruct;
  throw std::invalid_argument("Mesh type invalid.");
};

methodType SEMproxy::GetMethod(std::string method_arg) {
  if (method_arg == "sem") return methodType::kSem;
  if (method_arg == "dg") return methodType::kDg;
  throw std::invalid_argument("Method type invalid.");
};

model::AnisotropyType SEMproxy::GetAnisotropy(std::string anisotropy_arg) {
  if (anisotropy_arg == "iso") return model::AnisotropyType::kIso;
  if (anisotropy_arg == "vti") return model::AnisotropyType::kVTI;
  if (anisotropy_arg == "tti") return model::AnisotropyType::kTTI;
  throw std::invalid_argument("Anisotropy type invalid.");
};

float SEMproxy::FindCflDt(float cfl_factor) {
  float sqrt_dim3 = 1.73f;
  float min_spacing = mesh_->getMinSpacing();
  float v_max = mesh_->getMaxSpeed();
  return cfl_factor * min_spacing / (sqrt_dim3 * v_max);
}
void SEMproxy::InitMpi(int* mpi_init) {
#ifdef USE_MPI
  MPI_Initialized(mpi_init);
  if (mpi_init) {
    MPI_Comm_rank(MPI_COMM_WORLD, &dist_ctx_.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &dist_ctx_.size);
  } else {
    dist_ctx_.rank = 0;
    dist_ctx_.size = 1;
  }
#else
  dist_ctx_.rank = 0;
  dist_ctx_.size = 1;
#endif
}
void SEMproxy::InitSimParams(const SemProxyOptions& opt) {
  model::CartesianParams<float, int> global_params(opt.order, opt.ex, opt.ey, opt.ez, opt.lx, opt.ly, opt.lz,
                                                   opt.isModelOnNodes, opt.isElastic);
  global_params.isAcoustoElastic = opt.isAcoustoElastic;
  global_params.acoustoElasticBoundaryZ = opt.acoustoElasticBoundaryZ;
  global_params.origin_x = 0;
  model::CartesianXPartitioner<float, int> partitioner;
  local_params_ = partitioner.partition(global_params, dist_ctx_.rank, dist_ctx_.size);
  local_params_.isAcoustoElastic = opt.isAcoustoElastic;
  local_params_.acoustoElasticBoundaryZ = opt.acoustoElasticBoundaryZ;
  local_params_.model_file = opt.model_file;

  num_elements_[0] = local_params_.ex;
  num_elements_[1] = local_params_.ey;
  num_elements_[2] = local_params_.ez;
  num_nodes_[0] = local_params_.ex * opt.order + 1;
  num_nodes_[1] = local_params_.ey * opt.order + 1;
  num_nodes_[2] = local_params_.ez * opt.order + 1;
  domain_size_[0] = local_params_.lx;
  domain_size_[1] = local_params_.ly;
  domain_size_[2] = local_params_.lz;
  src_coord_[0] = opt.srcx;
  src_coord_[1] = opt.srcy;
  src_coord_[2] = opt.srcz;
  t_peak_ = opt.tpeak;
  f0_ = opt.f0;
  ricker_order_ = opt.ricker_order;
  rcv_coord_[0] = opt.rcvx;
  rcv_coord_[1] = opt.rcvy;
  rcv_coord_[2] = opt.rcvz;
  is_elastic_ = opt.isElastic;
}
void SEMproxy::InitMeshParams(const SemProxyOptions& opt) {
  const meshType mesh_type = GetMesh(opt.mesh);
  if (mesh_type == meshType::kStruct) {
    switch (opt.order) {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> builder(
            local_params_.ex, local_params_.lx, local_params_.ey, local_params_.ly, local_params_.ez, local_params_.lz,
            opt.isModelOnNodes, opt.isElastic, local_params_.origin_x, local_params_.origin_y, local_params_.origin_z,
            -1.0f, -1.0f, -1.0f, 0.0f, 0.0f, 0.0f, opt.isAcoustoElastic, opt.acoustoElasticBoundaryZ, opt.model_file);
        mesh_ = builder.getModel(opt.free_surface);
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> builder(
            local_params_.ex, local_params_.lx, local_params_.ey, local_params_.ly, local_params_.ez, local_params_.lz,
            opt.isModelOnNodes, opt.isElastic, local_params_.origin_x, local_params_.origin_y, local_params_.origin_z,
            -1.0f, -1.0f, -1.0f, 0.0f, 0.0f, 0.0f, opt.isAcoustoElastic, opt.acoustoElasticBoundaryZ, opt.model_file);
        mesh_ = builder.getModel(opt.free_surface);
        break;
      }
      case 3: {
        model::CartesianStructBuilder<float, int, 3> builder(
            local_params_.ex, local_params_.lx, local_params_.ey, local_params_.ly, local_params_.ez, local_params_.lz,
            opt.isModelOnNodes, opt.isElastic, local_params_.origin_x, local_params_.origin_y, local_params_.origin_z,
            -1.0f, -1.0f, -1.0f, 0.0f, 0.0f, 0.0f, opt.isAcoustoElastic, opt.acoustoElasticBoundaryZ, opt.model_file);
        mesh_ = builder.getModel(opt.free_surface);
        break;
      }
      case 4: {
        model::CartesianStructBuilder<float, int, 4> builder(
            local_params_.ex, local_params_.lx, local_params_.ey, local_params_.ly, local_params_.ez, local_params_.lz,
            opt.isModelOnNodes, opt.isElastic, local_params_.origin_x, local_params_.origin_y, local_params_.origin_z,
            -1.0f, -1.0f, -1.0f, 0.0f, 0.0f, 0.0f, opt.isAcoustoElastic, opt.acoustoElasticBoundaryZ, opt.model_file);
        mesh_ = builder.getModel(opt.free_surface);
        break;
      }
      case 5: {
        model::CartesianStructBuilder<float, int, 5> builder(
            local_params_.ex, local_params_.lx, local_params_.ey, local_params_.ly, local_params_.ez, local_params_.lz,
            opt.isModelOnNodes, opt.isElastic, local_params_.origin_x, local_params_.origin_y, local_params_.origin_z,
            -1.0f, -1.0f, -1.0f, 0.0f, 0.0f, 0.0f, opt.isAcoustoElastic, opt.acoustoElasticBoundaryZ, opt.model_file);
        mesh_ = builder.getModel(opt.free_surface);
        break;
      }
      default:
        throw std::runtime_error("Order other than 1-5 is not supported");
    }
  } else if (mesh_type == meshType::kUnstruct) {
    model::CartesianUnstructBuilder<float, int> builder(local_params_);
    mesh_ = builder.getModel(opt.free_surface);
  }
}

void SEMproxy::InitTopology() {
  par_topology_ =
      TopologyFactory::createFromMesh(*mesh_, dist_ctx_.rank, dist_ctx_.size, local_params_.origin_x, local_params_.lx);
}

void SEMproxy::InitSync() {
#if USE_MPI
  if (dist_ctx_.size > 1) {
    syncer_ = std::make_unique<BoundarySynchronizer>(std::make_unique<solver::fe::MPIBackend>());
  } else {
    syncer_ = std::make_unique<BoundarySynchronizer>(std::make_unique<SerialBackend>());
  }
#else
  syncer_ = std::make_unique<BoundarySynchronizer>(std::make_unique<SerialBackend>());
#endif
}

void SEMproxy::InitTimeParams(const SemProxyOptions& opt) {
  if (opt.autodt) {
    float cfl_factor = (opt.order == 2) ? 0.5f : 0.7f;
    dt_ = FindCflDt(cfl_factor);
  } else {
    dt_ = opt.dt;
  }
  time_max_ = opt.timemax;
  num_samples_ = time_max_ / dt_;
}

void SEMproxy::DisplayInitMsg(const SemProxyOptions& opt) {
  std::cout << "\n==========================================================\n";
  std::cout << "               SEM PROXY LAUNCH PARAMETERS                  \n";
  std::cout << "==========================================================\n";
  std::cout << "[Mesh Configuration]\n";
  std::cout << "  - Number of Nodes:      " << mesh_->getNumberOfNodes() << "\n";
  std::cout << "  - Number of Elements:   " << mesh_->getNumberOfElements() << "\n";
  std::cout << "  - Polynomial Order:     " << opt.order << "\n";
  std::cout << "  - Mesh Type:            " << opt.mesh << "\n\n";

  std::cout << "[Physics & Solver]\n";
  std::cout << "  - Method:               " << opt.method << "\n";
  std::cout << "  - Implementation:       " << opt.implem << "\n";
  std::cout << "  - Primary Physics:      "
            << (is_acousto_elastic_ ? "Acousto-Elastic" : (is_elastic_ ? "Elastic" : "Acoustic")) << "\n";
  std::cout << "  - Model Locality:       " << (opt.isModelOnNodes ? "Nodes" : "Elements") << "\n\n";

  std::cout << "[Time Domain]\n";
  std::cout << "  - Time Step (dt):       " << std::fixed << std::setprecision(6) << dt_ << " s\n";
  std::cout << "  - Total Sim Time:       " << time_max_ << " s\n";
  std::cout << "  - Total Samples:        " << num_samples_ << "\n\n";

  std::cout << "[Output Directives]\n";
  if (is_snapshots_) {
    std::cout << "  - Snapshots:            Enabled (Every " << snap_time_interval_ << " iterations)\n";
  } else {
    std::cout << "  - Snapshots:            Disabled\n";
  }
  std::cout << "==========================================================\n\n";
}

void SEMproxy::DisplayPerfMsg() const {
  double total_time = time_init_ + time_compute_ + time_io_;
  double pct_init = (total_time > 0) ? (time_init_ / total_time) * 100.0 : 0.0;
  double pct_comp = (total_time > 0) ? (time_compute_ / total_time) * 100.0 : 0.0;
  double pct_io = (total_time > 0) ? (time_io_ / total_time) * 100.0 : 0.0;

  std::cout << "\n==========================================================\n";
  std::cout << "               SEM PROXY PERFORMANCE SUMMARY                \n";
  std::cout << "==========================================================\n";
  std::cout << "  - Initialization Time:  " << std::fixed << std::setprecision(4) << std::setw(8) << time_init_
            << " s (" << std::setprecision(1) << std::setw(4) << pct_init << "%)\n";
  std::cout << "  - Compute Kernel Time:  " << std::fixed << std::setprecision(4) << std::setw(8) << time_compute_
            << " s (" << std::setprecision(1) << std::setw(4) << pct_comp << "%)\n";
  std::cout << "  - I/O & Output Time:    " << std::fixed << std::setprecision(4) << std::setw(8) << time_io_ << " s ("
            << std::setprecision(1) << std::setw(4) << pct_io << "%)\n";
  std::cout << "----------------------------------------------------------\n";
  std::cout << "  - Total Execution Time: " << std::fixed << std::setprecision(4) << std::setw(8) << total_time
            << " s\n";
  std::cout << "==========================================================\n\n";
}

void SEMproxy::WaitSnapshots() {
  for (auto& f : snapshot_futures_) {
    if (f.valid()) f.wait();
  }
  snapshot_futures_.clear();
}
