#include "DataOutput.h"

void DataOutput::Task1_RandomInput() {
  Task1<double> runner;
  auto results = runner.RandomInput();
  std::fstream out;
  out << std::fixed << std::setprecision(5);

  out.open("1_time_common_f.txt");
  for (const auto& item : results) {
    out << item.time_common_f << ", ";
  }
  out.close();

  out.open("1_time_optimized_f.txt");
  for (const auto& item : results) {
    out << item.time_optimized_f << ", ";
  }
  out.close();

  out.open("1_time_common_d.txt");
  for (const auto& item : results) {
    out << item.time_common_d << ", ";
  }
  out.close();
  out.open("1_time_optimized_d.txt");
  for (const auto& item : results) {
    out << item.time_optimized_d << ", ";
  }
  out.close();

  out.open("1_time_common_ld.txt");
  for (const auto& item : results) {
    out << item.time_common_ld << ", ";
  }
  out.close();
  out.open("1_time_optimized_ld.txt");
  for (const auto& item : results) {
    out << item.time_optimized_ld << ", ";
  }
  out.close();
}

void DataOutput::Task1_SingleVsMultiThread() {
  Task1<double> runner;
  auto results = runner.SingleVsMultiThread();
  std::fstream out;
  out << std::fixed << std::setprecision(5);

  out.open("1_time_thread_single.txt");
  for (auto item : std::get<0>(results)) {
    out << item << ", ";
  }
  out.close();

  out.open("1_time_thread_multi.txt");
  for (auto item : std::get<1>(results)) {
    out << item << ", ";
  }
  out.close();
}

void DataOutput::Task2_RandomBForSecondMatrix() {
  Task2<double> runner;
  auto results = runner.RandomBForSecondMatrix();
  std::fstream out;
  out << std::fixed << std::setprecision(10);

  out.open("2_fault_f.txt");
  for (auto item : std::get<0>(results)) {
    out << item << ", ";
  }
  out.close();

  out.open("2_fault_d.txt");
  for (auto item : std::get<1>(results)) {
    out << item << ", ";
  }
  out.close();

  out.open("2_fault_ld.txt");
  for (auto item : std::get<2>(results)) {
    out << item << ", ";
  }
  out.close();
}

void DataOutput::Task3_RandomInput() {
  Task3<double> runner;
  auto results = runner.RandomInput();
  std::fstream out;
  out << std::fixed << std::setprecision(5);

  out.open("3_time_tlu.txt");
  for (auto item : std::get<0>(results)) {
    out << item << ", ";
  }
  out.close();

  out.open("3_time_ldl.txt");
  for (auto item : std::get<1>(results)) {
    out << item << ", ";
  }
  out.close();
}
