#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARAMS_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARAMS_H_
namespace model {
/**
 * @brief Parameters describing a Cartesian mesh (or one MPI subdomain of it) and its model source.
 * @tparam Coord Type of the coordinates and lengths.
 * @tparam Index Type of the element counts.
 */
template <typename Coord, typename Index>
struct CartesianParams {
  int order;                      ///< Polynomial order of the elements.
  Index ex, ey, ez;               ///< Number of elements along x, y, z (local subdomain).
  Coord lx, ly, lz;               ///< Extent along x, y, z (local subdomain).
  bool isModelOnNodes;            ///< True if model properties are given per node, false if per element.
  bool isElastic;                 ///< True for an elastic model, false for an acoustic one.
  std::string model_file;         ///< Path of the model file; empty if no file is used.
  bool isAcoustoElastic{false};   ///< True if the domain couples an acoustic and an elastic part.
  Coord acoustoElasticBoundaryZ{static_cast<Coord>(0)};  ///< @todo VERIFY: z coordinate of the acousto-elastic interface, in which unit and frame?
  Coord DgSemBoundaryZ{static_cast<Coord>(0)};           ///< @todo VERIFY: z coordinate of the DG/SEM interface, in which unit and frame?

  /// @name Global domain (for MPI decomposition)
  /// @{
  Coord global_lx{0}, global_ly{0}, global_lz{0};
  Coord global_origin_x{0}, global_origin_y{0}, global_origin_z{0};
  /// @}

  /// @name Origin of the local subdomain
  /// @{
  Coord origin_x{0}, origin_y{0}, origin_z{0};
  /// @}

  /// @brief Leaves the members without default initializer uninitialized.
  CartesianParams() = default;

  /**
   * @brief Sets the local mesh size and model kind; the model file is left empty.
   * @param[in] order_ Polynomial order.
   * @param[in] ex_,ey_,ez_ Number of elements along x, y, z.
   * @param[in] lx_,ly_,lz_ Extent along x, y, z.
   * @param[in] isModelOnNodes_ True if model properties are given per node.
   * @param[in] isElastic_ True for an elastic model.
   */
  CartesianParams(int order_, Index ex_, Index ey_, Index ez_, Coord lx_, Coord ly_, Coord lz_, bool isModelOnNodes_,
                  bool isElastic_)
      : order(order_),
        ex(ex_),
        ey(ey_),
        ez(ez_),
        lx(lx_),
        ly(ly_),
        lz(lz_),
        isModelOnNodes(isModelOnNodes_),
        isElastic(isElastic_),
        model_file{""} {}
};
}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARAMS_H_
