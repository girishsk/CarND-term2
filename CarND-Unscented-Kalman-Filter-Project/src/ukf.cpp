#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  bool is_initialized_ = false
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a = 0.2;

  //Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd = 0.2;

  //define spreading parameter
  double lambda = 3 - n_aug;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    cout <<" Inside initialize" << endl;
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    ekf_.x_ = VectorXd(5);
    ekf_.x_ << 1, 1, 1, 1,1;
    ekf_.P_ = MatrixXd(5, 5);

    ekf_.F_ = MatrixXd(5, 5);

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      cout << "RADAR " << endl;
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      ekf_.x_[0] = rho*cos(phi);
      ekf_.x_[1] = rho*sin(phi);
      //TODO: initialize others

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "LASER " << endl;
      ekf_.x_[0] = measurement_pack.raw_measurements_[0];
      ekf_.x_[1] = measurement_pack.raw_measurements_[1];

    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  //TODO
  // predict and update
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  cout << "Time"<< dt << endl;
  int noise_ax = 9;
  int noise_ay = 9;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;


  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
        0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
        dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
        0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  if(dt !=0) {
    Prediction(dt);
  }
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

    UpdateRadar(measurement_pack.raw_measurements_);
  } else {
    // Laser updates

    UpdateLidar(measurement_pack.raw_measurements_);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 *
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {


  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  //create augmented mean state
  VectorXd x_aug_mean = VectorXd(7);
  x_aug_mean.head(5) = x;
  x_aug_mean[5]= 0;
  x_aug_mean[6]= 0;


  //create augmented covariance matrix
  P_aug.topLeftCorner( 5,5 ) = P;
  P_aug(5,5) = std_a * std_a;
  P_aug(6,6) = std_yawdd * std_yawdd;



  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug_mean;
  MatrixXd mul_a =  sqrt(n_aug + lambda)*A;
  for(int i = 0; i < n_aug; i ++){
    Xsig_aug.col(i+1) = x_aug_mean + mul_a.col(i) ;
    Xsig_aug.col(n_aug+i+1) = x_aug_mean - mul_a.col(i) ;
  }

  //write result
  *Xsig_out = Xsig_aug;

}



void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

  //calculate square root of P
  MatrixXd A = P.llt().matrixL();

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //your code goes here
  Xsig.col(0) = x;
  MatrixXd x_1_n = sqrt(lambda + n_x)*A;
//   std::cout << "x_1_n = " << std::endl << x_1_n << std::endl;
  for (int i =1 ; i <  n_x + 1; i++){
    Xsig.col(i) = x + x_1_n.col(i-1);
    Xsig.col(n_x+i) = x - x_1_n.col(i-1);
  }

  //calculate sigma points ...
  //set sigma points as columns of matrix Xsig

  //write result
  *Xsig_out = Xsig;

}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out,double delta_t,MatrixXd* Xsig_aug) {

  // input augmented sigma points

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  //predict sigma points
  //avoid division by zero
  for (int i = 0; i < 2 * n_aug + 1; i++) {
    if (Xsig_aug.col(i)[4] != 0) {
      VectorXd _x = Xsig_aug.col(i);
      VectorXd ad = VectorXd(5);
      ad[0] = _x[2] / _x[4] * (sin(_x[3] + _x[4] * delta_t) - sin(_x[3])) +
              0.5 * cos(_x[3]) * delta_t * delta_t * _x[5];
      ad[1] = _x[2] / _x[4] * (-1 * cos(_x[3] + _x[4] * delta_t) + cos(_x[3])) +
              0.5 * sin(_x[3]) * delta_t * delta_t * _x[5];
      ad[2] = delta_t * _x[5];
      ad[3] = _x[4] * delta_t + 0.5 * delta_t * delta_t * _x[6];
      ad[4] = delta_t * _x[6];

      Xsig_pred.col(i) = Xsig_aug.col(i).head(n_x) + ad;
    } else {
      VectorXd _x = Xsig_aug.col(i);
      VectorXd ad = VectorXd(5);
      ad[0] = _x[2] * sin(_x[3]) * delta_t + 0.5 * cos(_x[3]) * delta_t * delta_t * _x[5];
      ad[1] = _x[2] * -1 * cos(_x[3]) * delta_t + 0.5 * sin(_x[3]) * delta_t * delta_t * _x[5];
      ad[2] = delta_t * _x[5];
      ad[3] = _x[4] * delta_t + 0.5 * delta_t * delta_t * _x[6];
      ad[4] = delta_t * _x[6];
      Xsig_pred.col(i) = Xsig_aug.col(i).head(n_x) + ad;
    }
  }
}

  void UKF::PredictMeanAndCovariance(VectorXd *x_out, MatrixXd *P_out) {

    //create vector for weights
    VectorXd weights = VectorXd(2 * n_aug + 1);
    weights.fill(0);

    for (int i = 0; i < 2 * n_aug + 1; i++) {
      if (i == 0) {
        weights[i] = lambda / (n_aug + lambda);
      } else {
        weights[i] = 1 / (2 * (n_aug + lambda));
      }

      x += weights[i] * Xsig_pred.col(i);

    }
    for (int i = 0; i < 2 * n_aug + 1; i++) {
      P += weights[i] * ((Xsig_pred.col(i) - x) * ((Xsig_pred.col(i) - x).transpose()));
    }
    //set weights
    //predict state mean
    //predict state covariance matrix
  }

  void UKF::Prediction(double delta_t) {
    /**
    TODO:

    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */
    // TODO :


    // 1. Generate sigma points
    // GenerateSigmaPoints
    // AugmentedSigmaPoints

    // 2. Predict Sigma points
    // SigmaPointPrediction


    // 3. Predict mean and Covariance
    // PredictMeanAndCovariance


  }


  void UKF::PredictRadarMeasurement(VectorXd *z_out, MatrixXd *S_out) {

    MatrixXd R = MatrixXd(3, 3);
    R.fill(0);
    R(0, 0) = std_radr * std_radr;
    R(1, 1) = std_radphi * std_radphi;
    R(2, 2) = std_radrd * std_radrd;



    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug + 1; i++) {
      VectorXd _z = VectorXd(n_z);
      VectorXd _x = Xsig_pred.col(i);
      _z[0] = sqrt(_x[0] * _x[0] + _x[1] * _x[1]);
      _z[1] = atan2(_x[1], _x[0]);
      _z[2] = (_x[0] * cos(_x[3]) * _x[2] + _x[1] * sin(_x[3]) * _x[2]) / (sqrt(_x[0] * _x[0] + _x[1] * _x[1]));
      Zsig.col(i) = _z;
      z_pred += weights[i] * _z;

    }
    for (int i = 0; i < 2 * n_aug + 1; i++) {
      S += weights[i] * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();
    }
    S += R;
    //calculate mean predicted measurement

    //calculate measurement covariance matrix S

  }


  void UKF::UpdateState(VectorXd *x_out, MatrixXd *P_out) {

    //set vector for weights
    VectorXd weights = VectorXd(2 * n_aug + 1);
    double weight_0 = lambda / (lambda + n_aug);
    weights(0) = weight_0;
    for (int i = 1; i < 2 * n_aug + 1; i++) {
      double weight = 0.5 / (n_aug + lambda);
      weights(i) = weight;
    }


    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x, n_z);

    //calculate cross correlation matrix
    for (int i = 0; i < 2 * n_aug + 1; i++) {
      Tc += weights[i] * (Xsig_pred.col(i) - x) * (Zsig.col(i) - z_pred).transpose();
    }
    //calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    //update state mean and covariance matrix
    x = x + K * (z - z_pred);
    P = P - K * S * K.transpose();

  }

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
  void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */
  }

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
  void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */
  }

