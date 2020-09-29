#include "data_output.h"

namespace matrix::matrix_utils {

void DataOutput::Task1_RandomInput() {
  auto results = matrix_tasks::Task1<double>::RandomInput();
  std::ofstream out;
  out << std::fixed << std::setprecision(7);

  out.open("1_time_common_f.txt");
  for (const auto& item : results) {
    out << item.time_common_f << ", ";
  }
  out.close();

  out.open("1_time_tlu_f.txt");
  for (const auto& item : results) {
    out << item.time_tlu_f << ", ";
  }
  out.close();

  out.open("1_time_triangular_f.txt");
  for (const auto& item : results) {
    out << item.time_triangular_f << ", ";
  }
  out.close();

  out.open("1_time_common_d.txt");
  for (const auto& item : results) {
    out << item.time_common_d << ", ";
  }
  out.close();

  out.open("1_time_tlu_d.txt");
  for (const auto& item : results) {
    out << item.time_tlu_d << ", ";
  }
  out.close();

  out.open("1_time_triangular_d.txt");
  for (const auto& item : results) {
    out << item.time_triangular_d << ", ";
  }
  out.close();

  out.open("1_time_common_ld.txt");
  for (const auto& item : results) {
    out << item.time_common_ld << ", ";
  }
  out.close();

  out.open("1_time_tlu_ld.txt");
  for (const auto& item : results) {
    out << item.time_tlu_ld << ", ";
  }
  out.close();

  out.open("1_time_triangular_ld.txt");
  for (const auto& item : results) {
    out << item.time_triangular_ld << ", ";
  }
  out.close();
}

void DataOutput::Task1_SingleVsMultiThread() {
  auto results = matrix_tasks::Task1<double>::SingleVsMultiThread();
  std::ofstream out;
  out << std::fixed << std::setprecision(7);

  out.open("1_time_thread_single_tlu.txt");
  for (auto item : std::get<0>(results)) {
    out << item << ", ";
  }
  out.close();

  out.open("1_time_thread_multi_tlu.txt");
  for (auto item : std::get<1>(results)) {
    out << item << ", ";
  }
  out.close();

  out.open("1_time_thread_single_triangular.txt");
  for (auto item : std::get<2>(results)) {
    out << item << ", ";
  }
  out.close();

  out.open("1_time_thread_multi_triangular.txt");
  for (auto item : std::get<3>(results)) {
    out << item << ", ";
  }
  out.close();
}

void DataOutput::Task2_RandomBForSecondMatrix() {
  auto results = matrix_tasks::Task2<double>::RandomBForSecondMatrix();
  std::ofstream out;
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
  auto results = matrix_tasks::Task3<double>::RandomInput();
  std::ofstream out;
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

void DataOutput::Task5_DifferentOmegas() {
  auto result = matrix_tasks::Task5<double>::GetTime(0.8, 1.2);
  std::ofstream out;

  for (int i = 0; i < 3; ++i) {
    out.open("5_omega" + std::to_string(i) + ".txt");
    for (auto item : result[i]) {
      out << item << ", ";
    }
    out.close();
  }
}

}  // namespace matrix::matrix_utils
