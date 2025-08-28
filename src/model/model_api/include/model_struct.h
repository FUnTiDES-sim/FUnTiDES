#ifndef SRC_MODEL_MODELAPI_INCLUDE_MODEL_STRUCT_H_
#define SRC_MODEL_MODELAPI_INCLUDE_MODEL_STRUCT_H_

#include "model_api.h"
#include "gllpoints.h"

namespace model {

template<typename FloatType, typename ScalarType>
struct ModelStructData {
    public:
    // GPU-compatible special member functions
    PROXY_HOST_DEVICE ModelStructData() = default;
    PROXY_HOST_DEVICE ~ModelStructData() = default;
    PROXY_HOST_DEVICE ModelStructData(const ModelStructData&) = default;
    PROXY_HOST_DEVICE ModelStructData& operator=(const ModelStructData&) = default;

    ScalarType ex_, ey_, ez_;
    FloatType dx_, dy_, dz_;
    int order_;
};

/**
 * @brief Abstract base class representing a structured 3D mesh.
 */
template <typename FloatType, typename ScalarType>
class ModelStruct: public ModelApi<FloatType, ScalarType>{
 public:
    /**
     * @brief Default constructor.
     */
    PROXY_HOST_DEVICE ModelStruct() = default;

    /**
     * @brief Constructor from ModelStructData.
     * @param data ModelStructData structure containing all the mesh data
     */
    PROXY_HOST_DEVICE ModelStruct(const ModelStructData<FloatType, ScalarType>& data) :
        ex_(data.ex_), ey_(data.ey_), ez_(data.ez_),
        hx_(data.dx_), hy_(data.dy_), hz_(data.dz_),
        order_(data.order_) {

        nx_ = order_ * ex_ + 1;
        ny_ = order_ * ey_ + 1;
        nz_ = order_ * ez_ + 1;

        lx_ = ex_ * hx_;
        ly_ = ey_ * hy_;
        lz_ = ez_ * hz_;
    }

    /**
     * @brief Copy constructor.
     */
    PROXY_HOST_DEVICE ModelStruct(const ModelStruct&) = default;

    /**
     * @brief Assignment operator.
     */
    PROXY_HOST_DEVICE ModelStruct& operator=(const ModelStruct&) = default;

    /**
     * @brief Destructor.
     */
    PROXY_HOST_DEVICE ~ModelStruct() = default;

    /**
     * @brief Get the coordinate of a global node in the given dimension.
     * @param dofGlobal Global node index
     * @param dim Dimension index (0 = x, 1 = y, 2 = z)
     * @return Coordinate value in the specified dimension
     */
    PROXY_HOST_DEVICE
    FloatType nodeCoord(ScalarType dofGlobal, int dim) const {
        // Calculate total number of nodes per dimension
        int nodesPerDim[3];
        nodesPerDim[0] = (ex_ * order_) + 1;
        nodesPerDim[1] = (ey_ * order_) + 1;
        nodesPerDim[2] = (ez_ * order_) + 1;

        // Convert global node index to 3D node indices (i, j, k)
        int k = dofGlobal / (nodesPerDim[0] * nodesPerDim[1]);
        int remainder = dofGlobal % (nodesPerDim[0] * nodesPerDim[1]);
        int j = remainder / nodesPerDim[0];
        int i = remainder % nodesPerDim[0];

        int nodeIdx[3] = {i, j, k};

        // Determine which element this node belongs to and local position within element
        int elemIdx = nodeIdx[dim] / order_;  // Element index in the requested dimension
        int localIdx = nodeIdx[dim] % order_; // Local node index within element (0 to order_)

        // Handle boundary case: if we're at the last node of an element (except the last element),
        // it's actually the first node of the next element
        if (localIdx == order_ && elemIdx < (dim == 0 ? ex_ : (dim == 1 ? ey_ : ez_)) - 1) {
            elemIdx++;
            localIdx = 0;
        }

        // Get the GLL point coordinate in reference element [-1, 1]
        FloatType gllPoint = GLLPoints::get(order_, localIdx);

        // Map from reference element to physical element
        FloatType elementSize = (dim == 0) ? hx_ : ((dim == 1) ? hy_ : hz_);
        FloatType elementStart = elemIdx * elementSize;

        // Transform from [-1, 1] to physical coordinates
        FloatType physicalCoord = elementStart + (gllPoint + 1.0) * elementSize * 0.5;

        return physicalCoord;
    }


    /**
     * @brief Get the global node index for a local element-node triplet.
     * @param e Element index
     * @param i Local i-index in the element
     * @param j Local j-index in the element
     * @param k Local k-index in the element
     * @return Global node index
     */
    PROXY_HOST_DEVICE
    ScalarType globalNodeIndex(ScalarType e,
                               int i,
                               int j,
                               int k) const {
        ScalarType elemZ = e / (ex_* ey_);
        ScalarType tmp   = e % (ex_* ey_);
        ScalarType elemY = tmp / ex_;
        ScalarType elemX = tmp % ex_;

        int ix = elemX * order_ + i;
        int iy = elemY * order_ + j;
        int iz = elemZ * order_ + k;

        return ix + iy * nx_ + iz * nx_ * ny_;
    }

    /**
     * @brief Get the P-wave velocity value at a global node.
     * @param n Global node index
     * @return Model P-wave velocity value at the node
     */
    PROXY_HOST_DEVICE
    FloatType getModelVpOnNodes(ScalarType n) const {
        // TODO: Not returning magic number
        return 1500;
    }

    /**
     * @brief Get the average P-wave velocity value on a given element.
     * @param e Element index
     * @return Model P-wave velocity value for the element
     */
    PROXY_HOST_DEVICE
    FloatType getModelVpOnElement(ScalarType e) const {
        // TODO: Not returning a magic number
        return 1500;
    }

    /**
     * @brief Get the density value at a global node.
     * @param n Global node index
     * @return Model density value at the node
     */
    PROXY_HOST_DEVICE
    FloatType getModelRhoOnNodes(ScalarType n) const {
        // TODO: Not returning a magic number
        return 1;
    }

    /**
     * @brief Get the average density value on a given element.
     * @param e Element index
     * @return Model density value for the element
     */
    PROXY_HOST_DEVICE
    FloatType getModelRhoOnElement(ScalarType e) const {
        // TODO: Not returning a magic number
        return 1;
    }

    /**
     * @brief Get the total number of elements in the mesh.
     * @return Total element count
     */
    PROXY_HOST_DEVICE
    ScalarType getNumberOfElements() const {
        return ex_ * ey_ * ez_;
    }

    /**
     * @brief Get the total number of global nodes in the mesh.
     * @return Total node count
     */
    PROXY_HOST_DEVICE
    ScalarType getNumberOfNodes() const {
        return (order_ * ex_ + 1) * (order_ * ey_ + 1) * (order_ * ez_ + 1);
    }

    /**
     * @brief Get the number of interpolation points per element.
     * @return Number of interpolation points in one element
     */
    PROXY_HOST_DEVICE
    int getNumberOfPointsPerElement() const {
        int n_nodes_per_elem = (order_ + 1)
            * (order_ + 1) * (order_ + 1);
        return n_nodes_per_elem;
    }

    /**
     * @brief Get the polynomial order of the elements.
     * @return ORDER
     */
    PROXY_HOST_DEVICE
    int getOrder() const {
        return order_;
    }

    /**
     * @brief Get the boundary type of a given node.
     * @param n Global node index
     * @return A combination of BoundaryFlag values
     */
    PROXY_HOST_DEVICE
    BoundaryFlag boundaryType(ScalarType n) const {
        // throw std::runtime_error("boundaryType not implemented (model struct)");
        return BoundaryFlag::InteriorNode;
    }

    /**
     * @brief Compute the outward unit normal vector of an element face.
     * @param e Element index
     * @param dir Axis direction (0 = x, 1 = y, 2 = z)
     * @param face 1 = negative side, 2 = positive side in that direction
     * @param[out] v Output array (size 3) holding the normal vector
     */
    PROXY_HOST_DEVICE
    void faceNormal(ScalarType e,
                    int dir,
                    int face,
                    FloatType v[3]) const {
        // throw std::runtime_error("faceNormal not implemented (model struct)");
        return;
    }

    /**
     * @brief Get the size of the domain in the specified dimension.
     * @param dim The dimension index (0 for X, 1 for Y, 2 for Z).
     * @return The size of the domain along the specified dimension.
     */
    PROXY_HOST_DEVICE
    FloatType domainSize(int dim) const {
        switch(dim) {
            case 0: return lx_;
            case 1: return ly_;
            case 2: return lz_;
            default: return -1; // throw std::runtime_error("Incorrect dim in domain size");
        }
    }

 private:
  ScalarType ex_, ey_, ez_; // Nb elements in each direction
  ScalarType nx_, ny_, nz_; // Nb nodes in each direction
  FloatType lx_, ly_, lz_;      // domain size
  FloatType hx_, hy_, hz_;      // element size
  int order_;

};

}  // namespace model
#endif  // SRC_MODEL_MODELAPI_INCLUDE_MODEL_STRUCT_H_
