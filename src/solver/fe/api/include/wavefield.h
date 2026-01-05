#ifndef WAVEFIELD_H_
#define WAVEFIELD_H_

namespace solver
{
namespace fe
{
/**
 * @brief Base Wavefield data structure.
 */
struct Wavefield
{
  virtual ~Wavefield() = default;

  /**
   * @brief Swap data to advance the wavefield to the next time step.
   */
  virtual void advance() = 0;

  virtual void print() const = 0;
};
}  // namespace fe
}  // namespace solver

#endif  // WAVEFIELD_H_