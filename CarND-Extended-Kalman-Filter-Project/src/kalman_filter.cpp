#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
#include <math.h>


using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
#define PI 3.14159265

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  // std::cout << "Inside Predict" << std::endl;
  std::cout << "x_= " << x_ << std::endl;
  std::cout << "F_= " << F_ << std::endl;
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */

  VectorXd z_pred = VectorXd(3);
  float ms = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
  z_pred[0] = ms;
  z_pred[1] = atan2 (x_[1],x_[0]);
  cout<< "phi ::" << z_pred[1] << endl;
  if(ms == 0){
    z_pred[2] = 0;
  }
  else
    z_pred[2] = (x_[0]*x_[2] + x_[1]*x_[3])/ms;
  //VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
  while(y[1] > PI)
  {
    y[1] = y[1] - 2*PI;
  }
  while(y[1] < -1*PI)
  {
    y[1] = y[1] + 2*PI;
  }
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
