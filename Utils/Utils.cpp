#include "Utils.h"

namespace {

// Constant seed is used to make the results of the tests more predictable.
std::mt19937_64 random_generator(42);

}  // namespace

long double Random(int modulo, bool force_positive) {
  long double
      value = 1.0L + random_generator() % modulo + 1.0L / random_generator();
  if (!force_positive && random_generator() % 2 == 1) {
    value = -value;
  }
  return value;
}
