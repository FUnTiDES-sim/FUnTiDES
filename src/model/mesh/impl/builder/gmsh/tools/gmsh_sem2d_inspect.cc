#include <fstream>
#include <iostream>
#include <sstream>

#include "gmsh_2d_reader.h"

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  int status = 0;
  try {
    if (argc < 3) throw std::invalid_argument("Usage: gmsh_sem2d_inspect mesh.msh materials.csv [SEM order]");
    std::ifstream in(argv[2]);
    std::string line;
    std::getline(in, line);
    std::map<int, model::GmshAcousticMaterial> materials;
    while (std::getline(in, line)) {
      std::replace(line.begin(), line.end(), ',', ' ');
      std::istringstream row(line);
      int tag;
      float vp, rho;
      if (!(row >> tag >> vp >> rho) || !materials.emplace(tag, model::GmshAcousticMaterial{vp, rho}).second)
        throw std::runtime_error("Invalid material table");
    }
    auto data = model::readGmsh2D(argv[1], argc > 3 ? std::stoi(argv[3]) : 2, materials,
                                  {{1, model::BoundaryFlag::Surface}, {2, model::BoundaryFlag::Damping}});
    model::ModelUnstruct<float, int> mesh(data);
    int surface = 0, damping = 0;
    for (int i = 0; i < data.n_node_; ++i) {
      surface += mesh.isFreeSurface(i);
      damping += mesh.boundaryType(i) == model::BoundaryFlag::Damping;
    }
    std::cout << "{\"elements\":" << data.n_element_ << ",\"GLL_nodes\":" << data.n_node_
              << ",\"order\":" << data.order_ << ",\"free_surface_nodes\":" << surface
              << ",\"absorbing_nodes\":" << damping << "}\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    status = 1;
  }
  Kokkos::finalize();
  return status;
}
