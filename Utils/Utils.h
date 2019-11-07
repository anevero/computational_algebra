#ifndef UTILS_UTILS_H_
#define UTILS_UTILS_H_

#include <fstream>
#include <iomanip>
#include <random>

namespace matrix::matrix_utils {

void SetRandomSeed(unsigned long long seed);
long double Random(int modulo = 5000, bool force_positive = false);

}  // namespace matrix::matrix_utils

#endif  // UTILS_UTILS_H_
