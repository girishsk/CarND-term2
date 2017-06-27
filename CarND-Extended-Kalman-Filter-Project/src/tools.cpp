#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size() == 0 || estimations.size() != ground_truth.size()) {
      cout<<"Vectors not equal or vector size zero"<<endl;
      return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        VectorXd diff = estimations[i] - ground_truth[i];
        rmse = rmse.array() + diff.array()*diff.array();
  }

  //calculate the mean
  rmse /= estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE

	//check division by zero

	if(px == 0 && py==0){
	  cout << "Divide by zero error" << endl;
	  return Hj;
	}

	float denominator = px*px + py*py;



	Hj<< px/pow(denominator,0.5),py/pow(denominator,0.5),0,0,
	     -1*py/denominator,px/denominator,0,0,
	     py*(vx*py-vy*px)/pow(denominator,3/2),px*(vy*px-vx*py)/pow(denominator,3/2),px/pow(denominator,0.5),py/pow(denominator,0.5);

	//compute the Jacobian matrix

	return Hj;
}
