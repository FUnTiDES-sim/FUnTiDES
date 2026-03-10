#ifndef FUNTIDES_UTILS_INCLUDE_DISTRIBUTED_CTX_H_
#define FUNTIDES_UTILS_INCLUDE_DISTRIBUTED_CTX_H_
namespace utils
{

/**
 * @brief Hold the distributed context.
 *
 * Discribs the distributed context, such as number of parallel ranks, current
 * rank and later on the MPI context or other kind of communicator.
 * */
struct DistributedContext
{
  int rank{0};  //< Current rank
  int size{1};  //< Total number of rank
};

}  // namespace utils
#endif  // FUNTIDES_UTILS_INCLUDE_DISTRIBUTED_CTX_H_
