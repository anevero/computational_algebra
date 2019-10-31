#ifndef UTILS_DATAOUTPUT_H_
#define UTILS_DATAOUTPUT_H_

#include <fstream>
#include <iomanip>
#include <tuple>

#include "../Tasks/Task1.h"
#include "../Tasks/Task2.h"
#include "../Tasks/Task3.h"

class DataOutput {
 public:
  static void Task1_RandomInput();
  static void Task1_SingleVsMultiThread();
  static void Task2_RandomBForSecondMatrix();
  static void Task3_RandomInput();
};

#endif  // UTILS_DATAOUTPUT_H_
