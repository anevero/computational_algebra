#ifndef UTILS_DATA_OUTPUT_H_
#define UTILS_DATA_OUTPUT_H_

#include <fstream>
#include <iomanip>
#include <string>
#include <tuple>

#include "../tasks/task1.h"
#include "../tasks/task2.h"
#include "../tasks/task3.h"
#include "../tasks/task5.h"

namespace matrix::matrix_utils {

class DataOutput {
 public:
  static void Task1_RandomInput();
  static void Task1_SingleVsMultiThread();
  static void Task2_RandomBForSecondMatrix();
  static void Task3_RandomInput();
  static void Task5_DifferentOmegas();
};

}  // namespace matrix::matrix_utils

#endif  // UTILS_DATA_OUTPUT_H_
