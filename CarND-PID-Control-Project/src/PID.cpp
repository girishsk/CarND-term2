#include "PID.h"
#include <iostream>
#include <vector>

using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {}

PID::~PID() {}

void PID::Init(double Kp, double Ki, double Kd) {
  _Kp = Kp;
  _Ki = Ki;
  _Kd = Kd;
  cte_initialized = false;
  i_error = 0.0;
  vector<double> ctes ;


}

void PID::UpdateError(double cte) {

  if( !cte_initialized ){
    p_error = cte;
    cte_initialized = true;
  }
  d_error = cte - p_error;
  p_error = cte;

  i_error += cte;

}

double PID::TotalError() {
    return _Kp * p_error +_Kd * d_error + _Ki * i_error;
}

