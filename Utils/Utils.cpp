#include "Utils.h"

namespace {

std::mt19937_64 random_generator
    (std::chrono::system_clock::now().time_since_epoch().count());

}  // namespace

long double Random() {
  return (1.0L + random_generator() % 5000
      + 1.0L / (random_generator() % 100000));
}
