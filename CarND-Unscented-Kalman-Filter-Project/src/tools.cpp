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
  VectorXd rmse(5);
  rmse << 0,0,0,0,0;

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
