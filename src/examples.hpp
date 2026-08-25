#ifndef EXAMPLES_HPP
#define EXAMPLES_HPP

#include <string>

inline std::string sds_q_code_to_decimal(std::string q_code) {
  q_code.insert(1, ".");
  return q_code;
}

void run_sourced_eqn(void);
void run_coupled_eqn(void);
void run_teukolsky_precise_eqn(void);
void run_sds_precise_eqn(void);
void run_sds_areal_scan(const std::string &q_code, long long int s,
                        long long int l, const std::string &beta_text);

#endif
