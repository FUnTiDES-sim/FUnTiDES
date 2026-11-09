#include "gmsh_3d_reader.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "gllpoints.h"

namespace model {
namespace {
using Tag = long long;
using Face = std::array<Tag, 4>;
struct Point {
  double x, y, z;
};
struct Hex {
  std::array<Tag, 8> vertices;
  int material;
};
constexpr int corners[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
constexpr int faces[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {3, 2, 6, 7}, {0, 3, 7, 4}};
Face face(Face f) {
  std::sort(f.begin(), f.end());
  return f;
}
Face hexFace(const Hex& h, int f) {
  Face v{};
  for (int i = 0; i < 4; ++i) v[i] = h.vertices[faces[f][i]];
  return face(v);
}
Point interpolate(const Hex& h, const std::map<Tag, Point>& points, double x, double y, double z,
                  double* determinant = nullptr) {
  double t[3] = {x, y, z}, J[3][3]{};
  Point p{};
  for (int v = 0; v < 8; ++v) {
    double a[3];
    for (int d = 0; d < 3; ++d) a[d] = (1 + (2 * corners[v][d] - 1) * t[d]) / 2;
    auto q = points.at(h.vertices[v]);
    double coord[3] = {q.x, q.y, q.z};
    p.x += a[0] * a[1] * a[2] * q.x;
    p.y += a[0] * a[1] * a[2] * q.y;
    p.z += a[0] * a[1] * a[2] * q.z;
    for (int d = 0; d < 3; ++d)
      for (int c = 0; c < 3; ++c) J[c][d] += coord[c] * (2 * corners[v][d] - 1) * 0.5 * a[(d + 1) % 3] * a[(d + 2) % 3];
  }
  if (determinant)
    *determinant = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                   J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                   J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
  return p;
}
void require(bool ok, const std::string& message) {
  if (!ok) throw std::runtime_error("Gmsh 3D: " + message);
}
}  // namespace

ModelUnstructData<float, int> readGmsh3D(const std::filesystem::path& file, int order,
                                         const std::map<int, GmshAcousticMaterial>& materials,
                                         const std::map<int, BoundaryFlag>& boundaries) {
  require(order >= 1 && order <= MAX_GLL_ORDER, "unsupported GLL order");
  for (const auto& item : materials)
    require(
        item.second.vp > 0 && item.second.rho > 0 && std::isfinite(item.second.vp) && std::isfinite(item.second.rho),
        "materials must have finite positive velocity and density");
  std::ifstream in(file);
  require(bool(in), "cannot open " + file.string());
  std::map<Tag, Point> points;
  std::vector<Hex> hexes;
  std::map<Face, int> boundary_faces;
  std::set<Tag> element_ids;
  bool format = false, have_nodes = false, have_elements = false;
  std::string section;
  while (in >> section) {
    if (section == "$MeshFormat") {
      double version;
      int binary, width;
      require(bool(in >> version >> binary >> width), "invalid mesh format");
      require(version == 2.2 && binary == 0 && width == 8,
              "expected MSH 2.2 ASCII; export with Mesh.MshFileVersion=2.2 and Mesh.Binary=0");
      std::string end;
      in >> end;
      require(end == "$EndMeshFormat", "missing format terminator");
      format = true;
    } else if (section == "$Nodes") {
      require(format && !have_nodes, "missing format or duplicate node section");
      int count;
      require(bool(in >> count) && count > 0, "invalid node count");
      for (int i = 0; i < count; ++i) {
        Tag tag;
        double x, y, z;
        require(bool(in >> tag >> x >> y >> z) && tag > 0 && std::isfinite(x + y + z), "invalid node record");
        require(points.emplace(tag, Point{x, y, z}).second, "duplicate node tag");
      }
      std::string end;
      in >> end;
      require(end == "$EndNodes", "missing node terminator");
      have_nodes = true;
    } else if (section == "$Elements") {
      require(have_nodes && !have_elements, "missing nodes or duplicate element section");
      int count;
      require(bool(in >> count) && count > 0, "invalid element count");
      std::string line;
      std::getline(in, line);
      for (int i = 0; i < count; ++i) {
        require(bool(std::getline(in, line)), "truncated element section");
        std::istringstream row(line);
        Tag id;
        int type, ntags;
        require(bool(row >> id >> type >> ntags) && id > 0 && ntags >= 0, "invalid element record");
        require(element_ids.insert(id).second, "duplicate element tag");
        int physical = 0;
        for (int j = 0; j < ntags; ++j) {
          int t;
          require(bool(row >> t), "truncated physical tags");
          if (j == 0) physical = t;
        }
        require(type == 1 || type == 3 || type == 5 || type == 15,
                "only Hex8 volumes, Q4 surfaces, lines and points are supported; found element type " +
                    std::to_string(type));
        int n = type == 5 ? 8 : type == 3 ? 4 : type == 1 ? 2 : 1;
        std::array<Tag, 8> v{};
        for (int j = 0; j < n; ++j)
          require(bool(row >> v[j]) && points.count(v[j]), "unknown or missing node in element");
        std::string extra;
        require(!(row >> extra), "unexpected element data");
        if (type == 5) {
          require(materials.count(physical), "unknown physical material tag " + std::to_string(physical));
          require(std::set<Tag>(v.begin(), v.end()).size() == 8, "repeated hex corner");
          Hex h{v, physical};
          for (int k = 0; k <= 2 * order; ++k)
            for (int j = 0; j <= 2 * order; ++j)
              for (int i = 0; i <= 2 * order; ++i) {
                double det;
                interpolate(h, points, -1.0 + double(i) / order, -1.0 + double(j) / order, -1.0 + double(k) / order,
                            &det);
                require(std::isfinite(det) && det > 0, "nonpositive Jacobian in hex " + std::to_string(id));
              }
          hexes.push_back(h);
        } else if (type == 3) {
          Face f{v[0], v[1], v[2], v[3]};
          require(std::set<Tag>(f.begin(), f.end()).size() == 4 && boundaries.count(physical),
                  "unknown or degenerate physical boundary surface");
          require(boundary_faces.emplace(face(f), physical).second, "duplicate boundary face");
        }
      }
      std::string end;
      in >> end;
      require(end == "$EndElements", "missing element terminator");
      have_elements = true;
    } else {
      require(section.size() > 1 && section[0] == '$' && section.rfind("$End", 0) != 0, "unexpected section");
      std::string end = "$End" + section.substr(1), line;
      bool found = false;
      while (std::getline(in, line)) {
        if (line == end || line == end + "\r") {
          found = true;
          break;
        }
      }
      require(found, "unterminated " + section);
    }
  }
  require(format && have_nodes && have_elements && !hexes.empty(), "incomplete or empty hexahedral mesh");
  std::map<Face, int> incidence;
  for (const auto& h : hexes)
    for (int f = 0; f < 6; ++f) ++incidence[hexFace(h, f)];
  for (const auto& item : incidence) {
    require(item.second <= 2, "nonmanifold face");
    require(item.second == 1 ? boundary_faces.count(item.first) == 1 : boundary_faces.count(item.first) == 0,
            "unmarked external face or boundary surface on internal face; check conformity");
  }
  for (const auto& item : boundary_faces)
    require(incidence.count(item.first), "boundary surface is not a hexahedron face");

  // Integer tensor weights identify a node independently of the local face
  // rotation/reflection, without merging geometrically close unrelated nodes.
  using NodeKey = std::vector<std::pair<Tag, int>>;
  std::map<NodeKey, int> shared;
  std::vector<Point> coordinates;
  std::vector<BoundaryFlag> flags;
  std::vector<std::vector<int>> connectivity;
  const int n = order + 1, npe = n * n * n;
  for (const auto& h : hexes) {
    std::vector<int> nodes(npe);
    std::array<int, 6> boundary_tags{};
    std::array<bool, 6> marked{};
    for (int f = 0; f < 6; ++f) {
      auto it = boundary_faces.find(hexFace(h, f));
      if (it != boundary_faces.end()) {
        boundary_tags[f] = it->second;
        marked[f] = true;
      }
    }
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
          int t[3] = {i, j, k};
          NodeKey key;
          for (int v = 0; v < 8; ++v) {
            int weight = 1;
            for (int d = 0; d < 3; ++d) weight *= corners[v][d] ? t[d] : order - t[d];
            if (weight) key.emplace_back(h.vertices[v], weight);
          }
          std::sort(key.begin(), key.end());
          auto p = interpolate(h, points, GLLPoints::get(order, i), GLLPoints::get(order, j), GLLPoints::get(order, k));
          auto entry = shared.emplace(key, coordinates.size());
          int id = entry.first->second;
          if (entry.second) {
            coordinates.push_back(p);
            flags.push_back(BoundaryFlag::InteriorNode);
          }
          nodes[i + n * j + n * n * k] = id;
          bool on[6] = {k == 0, k == order, j == 0, i == order, j == order, i == 0};
          for (int f = 0; f < 6; ++f) {
            // A physical tag can be zero when explicitly mapped by the caller.
            if (on[f] && marked[f]) {
              auto flag = boundaries.at(boundary_tags[f]);
              if (flags[id] != BoundaryFlag::Surface || flag == BoundaryFlag::Surface) flags[id] = flag;
            }
          }
        }
    connectivity.push_back(std::move(nodes));
  }
  ModelUnstructData<float, int> d;
  d.dim_ = 3;
  d.order_ = order;
  d.n_element_ = hexes.size();
  d.n_node_ = coordinates.size();
  double minx = coordinates[0].x, maxx = minx, minz = coordinates[0].z, maxz = minz;
  double miny = coordinates[0].y, maxy = miny;
  for (auto p : coordinates) {
    minx = std::min(minx, p.x);
    maxx = std::max(maxx, p.x);
    miny = std::min(miny, p.y);
    maxy = std::max(maxy, p.y);
    minz = std::min(minz, p.z);
    maxz = std::max(maxz, p.z);
  }
  d.ox_ = minx;
  d.oy_ = miny;
  d.oz_ = minz;
  d.lx_ = maxx - minx;
  d.ly_ = maxy - miny;
  d.lz_ = maxz - minz;
  d.isModelOnNodes_ = false;
  d.isElastic_ = false;
  d.global_node_index_ = allocateArray2D<dualArrayInt>(d.n_element_, npe, "Gmsh GLL connectivity");
  d.nodes_coords_x_ = allocateVector<dualVectorReal>(d.n_node_, "Gmsh x");
  d.nodes_coords_y_ = allocateVector<dualVectorReal>(d.n_node_, "Gmsh y");
  d.nodes_coords_z_ = allocateVector<dualVectorReal>(d.n_node_, "Gmsh z");
  d.boundaries_t_ = allocateVector<dualVectorInt>(d.n_node_, "Gmsh boundaries");
  d.model_vp_element_ = allocateVector<dualVectorReal>(d.n_element_, "Gmsh velocity");
  d.model_rho_element_ = allocateVector<dualVectorReal>(d.n_element_, "Gmsh density");
  for (int n = 0; n < d.n_node_; ++n) {
    d.nodes_coords_x_.view_host()(n) = coordinates[n].x;
    d.nodes_coords_y_.view_host()(n) = coordinates[n].y;
    d.nodes_coords_z_.view_host()(n) = coordinates[n].z;
    d.boundaries_t_.view_host()(n) = static_cast<int>(flags[n]);
  }
  for (int e = 0; e < d.n_element_; ++e) {
    for (int k = 0; k < npe; ++k) d.global_node_index_.view_host()(e, k) = connectivity[e][k];
    const auto& mat = materials.at(hexes[e].material);
    d.model_vp_element_.view_host()(e) = mat.vp;
    d.model_rho_element_.view_host()(e) = mat.rho;
  }
  for (auto* v :
       {&d.nodes_coords_x_, &d.nodes_coords_y_, &d.nodes_coords_z_, &d.model_vp_element_, &d.model_rho_element_})
    v->modify_host();
  d.global_node_index_.modify_host();
  d.boundaries_t_.modify_host();
  return d;
}
}  // namespace model
