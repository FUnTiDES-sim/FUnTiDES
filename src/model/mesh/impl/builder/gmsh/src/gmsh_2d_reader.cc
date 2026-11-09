#include "gmsh_2d_reader.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "gllpoints.h"

namespace model {
namespace {
using Tag = long long;
using Edge = std::pair<Tag, Tag>;
struct Point {
  double x, z;
};
struct Quad {
  std::array<Tag, 4> vertices;
  int material;
};
Edge edge(Tag a, Tag b) { return std::minmax(a, b); }
void require(bool ok, const std::string& message) {
  if (!ok) throw std::runtime_error("Gmsh 2D: " + message);
}
}  // namespace

ModelUnstructData<float, int> readGmsh2D(const std::filesystem::path& file, int order,
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
  std::vector<Quad> quads;
  std::map<Edge, int> boundary_edges;
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
        require(std::abs(y) < 1e-8, "mesh must lie in the XZ plane (y=0)");
        require(points.emplace(tag, Point{x, z}).second, "duplicate node tag");
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
        require(type == 1 || type == 3 || type == 15,
                "only linear quads, lines and points are supported; found element type " + std::to_string(type));
        int n = type == 3 ? 4 : type == 1 ? 2 : 1;
        std::array<Tag, 4> v{};
        for (int j = 0; j < n; ++j)
          require(bool(row >> v[j]) && points.count(v[j]), "unknown or missing node in element");
        std::string extra;
        require(!(row >> extra), "unexpected element data");
        if (type == 3) {
          require(materials.count(physical), "unknown physical material tag " + std::to_string(physical));
          require(std::set<Tag>(v.begin(), v.end()).size() == 4, "repeated quad corner");
          double area = 0;
          for (int j = 0; j < 4; ++j) {
            auto a = points.at(v[j]), b = points.at(v[(j + 1) % 4]);
            area += a.x * b.z - b.x * a.z;
          }
          require(std::abs(area) > 1e-12, "degenerate quad");
          if (area < 0) std::swap(v[1], v[3]);
          for (int j = 0; j < 4; ++j) {
            auto a = points.at(v[j]), b = points.at(v[(j + 1) % 4]), c = points.at(v[(j + 2) % 4]);
            require((b.x - a.x) * (c.z - b.z) - (b.z - a.z) * (c.x - b.x) > 1e-12, "nonconvex or folded quad");
          }
          quads.push_back({v, physical});
        } else if (type == 1) {
          require(v[0] != v[1] && boundaries.count(physical), "unknown or degenerate physical boundary line");
          require(boundary_edges.emplace(edge(v[0], v[1]), physical).second, "duplicate boundary edge");
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
  require(format && have_nodes && have_elements && !quads.empty(), "incomplete or empty quadrilateral mesh");
  std::map<Edge, int> incidence;
  for (const auto& q : quads)
    for (int j = 0; j < 4; ++j) ++incidence[edge(q.vertices[j], q.vertices[(j + 1) % 4])];
  for (const auto& item : incidence) {
    require(item.second <= 2, "nonmanifold edge");
    require(item.second == 1 ? boundary_edges.count(item.first) == 1 : boundary_edges.count(item.first) == 0,
            "unmarked external edge or boundary line on internal edge; check mesh conformity");
  }
  for (const auto& item : boundary_edges)
    require(incidence.count(item.first), "boundary line is not a quadrilateral edge");

  std::vector<Point> coordinates;
  std::map<Tag, int> vertices;
  std::map<std::tuple<Tag, Tag, int>, int> edge_nodes;
  std::vector<std::vector<int>> connectivity;
  auto vertex = [&](Tag tag) {
    auto it = vertices.find(tag);
    if (it != vertices.end()) return it->second;
    int n = coordinates.size();
    coordinates.push_back(points.at(tag));
    vertices[tag] = n;
    return n;
  };
  for (const auto& q : quads) {
    std::vector<int> nodes((order + 1) * (order + 1));
    for (int j = 0; j <= order; ++j)
      for (int i = 0; i <= order; ++i) {
        int corner = -1;
        if (i == 0 && j == 0)
          corner = 0;
        else if (i == order && j == 0)
          corner = 1;
        else if (i == order && j == order)
          corner = 2;
        else if (i == 0 && j == order)
          corner = 3;
        if (corner >= 0) {
          nodes[i + j * (order + 1)] = vertex(q.vertices[corner]);
          continue;
        }
        Tag a = 0, b = 0;
        int position = 0;
        if (j == 0) {
          a = q.vertices[0];
          b = q.vertices[1];
          position = i;
        } else if (i == order) {
          a = q.vertices[1];
          b = q.vertices[2];
          position = j;
        } else if (j == order) {
          a = q.vertices[3];
          b = q.vertices[2];
          position = i;
        } else if (i == 0) {
          a = q.vertices[0];
          b = q.vertices[3];
          position = j;
        }
        auto key = std::make_tuple(std::min(a, b), std::max(a, b), a < b ? position : order - position);
        if (a && edge_nodes.count(key)) {
          nodes[i + j * (order + 1)] = edge_nodes.at(key);
          continue;
        }
        double xi = GLLPoints::get(order, i), eta = GLLPoints::get(order, j);
        std::array<double, 4> weights{(1 - xi) * (1 - eta) / 4, (1 + xi) * (1 - eta) / 4, (1 + xi) * (1 + eta) / 4,
                                      (1 - xi) * (1 + eta) / 4};
        Point p{0, 0};
        for (int k = 0; k < 4; ++k) {
          p.x += weights[k] * points.at(q.vertices[k]).x;
          p.z += weights[k] * points.at(q.vertices[k]).z;
        }
        int n = coordinates.size();
        coordinates.push_back(p);
        nodes[i + j * (order + 1)] = n;
        if (a) edge_nodes[key] = n;
      }
    connectivity.push_back(std::move(nodes));
  }
  std::vector<BoundaryFlag> flags(coordinates.size(), BoundaryFlag::InteriorNode);
  for (const auto& item : boundary_edges) {
    auto flag = boundaries.at(item.second);
    auto apply = [&](int node) {
      if (flags[node] != BoundaryFlag::Surface || flag == BoundaryFlag::Surface) flags[node] = flag;
    };
    apply(vertices.at(item.first.first));
    apply(vertices.at(item.first.second));
    for (int k = 1; k < order; ++k) apply(edge_nodes.at(std::make_tuple(item.first.first, item.first.second, k)));
  }
  ModelUnstructData<float, int> d;
  d.dim_ = 2;
  d.order_ = order;
  d.n_element_ = quads.size();
  d.n_node_ = coordinates.size();
  double minx = coordinates[0].x, maxx = minx, minz = coordinates[0].z, maxz = minz;
  for (auto p : coordinates) {
    minx = std::min(minx, p.x);
    maxx = std::max(maxx, p.x);
    minz = std::min(minz, p.z);
    maxz = std::max(maxz, p.z);
  }
  d.ox_ = minx;
  d.oy_ = 0;
  d.oz_ = minz;
  d.lx_ = maxx - minx;
  d.ly_ = 0;
  d.lz_ = maxz - minz;
  d.isModelOnNodes_ = false;
  d.isElastic_ = false;
  d.global_node_index_ =
      allocateArray2D<dualArrayInt>(d.n_element_, (order + 1) * (order + 1), "Gmsh GLL connectivity");
  d.nodes_coords_x_ = allocateVector<dualVectorReal>(d.n_node_, "Gmsh x");
  d.nodes_coords_y_ = allocateVector<dualVectorReal>(d.n_node_, "Gmsh y");
  d.nodes_coords_z_ = allocateVector<dualVectorReal>(d.n_node_, "Gmsh z");
  d.boundaries_t_ = allocateVector<dualVectorInt>(d.n_node_, "Gmsh boundaries");
  d.model_vp_element_ = allocateVector<dualVectorReal>(d.n_element_, "Gmsh velocity");
  d.model_rho_element_ = allocateVector<dualVectorReal>(d.n_element_, "Gmsh density");
  for (int n = 0; n < d.n_node_; ++n) {
    d.nodes_coords_x_.view_host()(n) = coordinates[n].x;
    d.nodes_coords_z_.view_host()(n) = coordinates[n].z;
    d.boundaries_t_.view_host()(n) = static_cast<int>(flags[n]);
  }
  for (int e = 0; e < d.n_element_; ++e) {
    for (int k = 0; k < (order + 1) * (order + 1); ++k) d.global_node_index_.view_host()(e, k) = connectivity[e][k];
    const auto& mat = materials.at(quads[e].material);
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
