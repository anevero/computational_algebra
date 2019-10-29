#ifndef DATAOUTPUT_H_
#define DATAOUTPUT_H_

#include "../Tasks/Task1.h"
#include "../Tasks/Task2.h"
#include "../Tasks/Task3.h"

#include <fstream>
#include <iomanip>
#include <tuple>

class DataOutput {
 public:
  static void Task1_RandomInput();
  static void Task1_SingleVsMultiThread();
  static void Task2_RandomBForSecondMatrix();
  static void Task3_RandomInput();
};

#endif  // DATAOUTPUT_H_
