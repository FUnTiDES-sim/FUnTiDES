#ifndef FUNTIDES_PARALLEL_INCLUDE_DISTRIBUTED_CTX_H_
#define FUNTIDES_PARALLEL_INCLUDE_DISTRIBUTED_CTX_H_
namespace utils {

/**
 * @brief Identifies this process within the set of parallel ranks.
 *
 * Holds the rank count and the current rank. Defaults describe a
 * single-process run.
 */
struct DistributedContext {
  int rank{0};  ///< Index of the current rank, in [0, size).
  int size{1};  ///< Total number of ranks.
};

}  // namespace utils
#endif  // FUNTIDES_PARALLEL_INCLUDE_DISTRIBUTED_CTX_H_
