#ifndef EXAMPLES_HPP
#define EXAMPLES_HPP

#include <charconv>
#include <stdexcept>
#include <string>
#include <system_error>

inline std::string sds_q_code_to_decimal(std::string q_code) {
  q_code.insert(1, ".");
  return q_code;
}

inline long long int sds_parse_integer_argument(
    const std::string &text, const char *name) {
  const char *begin = text.data();
  const char *const end = begin + text.size();
  if(begin != end && *begin == '+') ++begin;

  long long int value = 0;
  const auto result = std::from_chars(begin, end, value);
  if(begin == end || result.ec != std::errc{} || result.ptr != end) {
    throw std::invalid_argument(std::string(name) + " must be an integer");
  }
  return value;
}

void run_sourced_eqn(void);
void run_coupled_eqn(void);
void run_teukolsky_precise_eqn(void);
void run_sds_precise_eqn(void);
void run_sds_areal_scan(const std::string &q_code, long long int s,
                        long long int l, const std::string &beta_text);

#endif
