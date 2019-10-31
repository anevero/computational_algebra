#include "Utils.h"

namespace {

std::mt19937_64 random_generator
    (std::chrono::system_clock::now().time_since_epoch().count());

}  // namespace

long double Random(int modulo, bool force_positive) {
  long double
      value = 1.0L + random_generator() % modulo + 1.0L / random_generator();
  if (!force_positive && random_generator() % 2 == 1) {
    value = -value;
  }
  return value;
}
