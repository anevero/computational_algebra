#include "utils.h"

namespace matrix::matrix_utils {

// Constant seed is used to make the results of the tests more predictable.
std::mt19937_64 random_generator = std::mt19937_64{42};

void SetRandomSeed(unsigned long long seed) {
  random_generator.seed(seed);
}

}  // namespace matrix::matrix_utils
