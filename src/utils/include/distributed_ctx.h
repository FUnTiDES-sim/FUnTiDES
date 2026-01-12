#ifndef SRC_UTILS_INCLUDE_DISTRIBUTED_CXT_H
#define SRC_UTILS_INCLUDE_DISTRIBUTED_CXT_H

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

#endif  // SRC_UTILS_INCLUDE_DISTRIBUTED_CXT_H
