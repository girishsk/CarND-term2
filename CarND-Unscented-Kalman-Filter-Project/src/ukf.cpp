#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
#define PI 3.14159265

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    bool is_initialized_ = false;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 1.5;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.57;

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
    n_x_ = 5;

    //set augmented dimension
    n_aug_ = 7;
    n_z_ = 3;

////  //Process noise standard deviation longitudinal acceleration in m/s^2
//  std_a_ = 0.2;
//    //Process noise standard deviation yaw acceleration in rad/s^2
//  std_yawdd_ = 0.2;

    // Process noise standard deviation longitudinal acceleration in m/s^2
//    std_a_ = 6;

    // Process noise standard deviation yaw acceleration in rad/s^2
//    std_yawdd_ = 0.5;
//
//    // Laser measurement noise standard deviation position1 in m
//    std_laspx_ = 0.15;
//
//    // Laser measurement noise standard deviation position2 in m
//    std_laspy_ = 0.15;


    S_ = MatrixXd(n_z_,n_z_);

    //define spreading parameter
    lambda_ = 3 - n_x_;
    weights_ = VectorXd(2 * n_aug_+ 1);
    z_pred_= VectorXd(n_z_);
    R_laser_ = MatrixXd(2, 2);
    H_laser_ = MatrixXd(2, 5);

//    //radar measurement noise standard deviation radius in m
//    std_radr_ = 0.3;
//
//    //radar measurement noise standard deviation angle in rad
//    std_radphi_ = 0.0175;
//
//    //radar measurement noise standard deviation radius change in m/s
//    std_radrd_ = 0.1;



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
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
    /**
    TODO:

    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
    if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        x_ << 0, 0, 0, 0,0;
        P_ = MatrixXd(5, 5);
        P_ = MatrixXd::Identity(n_x_, n_x_);
        H_laser_ << 1.0, 0.0, 0.0, 0.0,0.0,
                0.0, 1.0, 0.0, 0.0,0.0;


        previous_timestamp_ = measurement_pack.timestamp_;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            float rho = measurement_pack.raw_measurements_[0];
            float phi = measurement_pack.raw_measurements_[1];
            float rho_dot = measurement_pack.raw_measurements_[2];
            float vx = rho_dot * cos(phi);
            float vy = rho_dot * sin(phi);
            float v  = sqrt(vx * vx + vy * vy);
            x_[0] = rho*cos(phi);
            x_[1] = rho*sin(phi);
            x_[2] = v;
            //TODO: initialize others
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            x_[0] = measurement_pack.raw_measurements_[0];
            x_[1] = measurement_pack.raw_measurements_[1];

        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
    // predict and update
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

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
        VectorXd x = VectorXd(2);
        double rho = measurement_pack.raw_measurements_[0];
        double phi = measurement_pack.raw_measurements_[1];
        double  rho_dot = measurement_pack.raw_measurements_[2];
        x[0] = rho*cos(phi);
        x[1] = rho*sin(phi);


        UpdateRadar(measurement_pack);
    } else {
        // Laser updates

        UpdateLidar(measurement_pack);
    }
    cout<<"*******************************************************************************"<<endl;
}



void UKF::Prediction(double delta_t) {
    /**
    TODO:

    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */

    // GenerateSigmaPoints
    MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_+ 1);
    Xsig.fill(0);
    //calculate square root of P
    MatrixXd A = P_.llt().matrixL();

    Xsig.col(0) = x_;
    MatrixXd x_1_n = sqrt(lambda_ + n_aug_)*A;

    for (int i =1 ; i <  n_x_ + 1; i++){
        VectorXd _x = VectorXd(5);
        _x = x_ + x_1_n.col(i-1);
        Xsig.col(i) = _x;
        _x = x_ - x_1_n.col(i-1);
        Xsig.col(n_x_+i) = _x;
    }


    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_+ 1);
    Xsig_aug.fill(0);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);
    P_aug.fill(0);


    //create augmented mean state
    VectorXd x_aug_mean = VectorXd(7);
    x_aug_mean.fill(0);
    x_aug_mean.head(5) = x_;
    x_aug_mean[5]= 0;
    x_aug_mean[6]= 0;


    //create augmented covariance matrix
    P_aug.topLeftCorner( 5,5 ) = P_;
    P_aug(5,5) = std_a_ * std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;


    lambda_ = 3 - n_aug_;
    //create square root matrix
    MatrixXd A_aug = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0) = x_aug_mean;
    MatrixXd mul_a =  sqrt(n_aug_+ lambda_)*A_aug;
    for(int i = 0; i < n_aug_; i ++){
        Xsig_aug.col(i+1) = x_aug_mean + mul_a.col(i) ;
        Xsig_aug.col(n_aug_+i+1) = x_aug_mean - mul_a.col(i) ;
    }


    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_+ 1);
    Xsig_pred.fill(0);


    //predict sigma points
    //avoid division by zero
    for (int i = 0; i < 2 * n_aug_+ 1; i++) {
        if (fabs(Xsig_aug.col(i)[4]) > 0.001) {
            VectorXd _x = Xsig_aug.col(i);
            VectorXd ad = VectorXd(5);
            ad[0] = _x[2] / _x[4] * (sin(_x[3] + _x[4] * delta_t) - sin(_x[3])) +
                    0.5 * cos(_x[3]) * delta_t * delta_t * _x[5];
            ad[1] = _x[2] / _x[4] * (-1 * cos(_x[3] + _x[4] * delta_t) + cos(_x[3])) +
                    0.5 * sin(_x[3]) * delta_t * delta_t * _x[5];
            ad[2] = delta_t * _x[5];
            ad[3] = _x[4] * delta_t + 0.5 * delta_t * delta_t * _x[6];
            ad[4] = delta_t * _x[6];

            Xsig_pred.col(i) = Xsig_aug.col(i).head(n_x_) + ad;
        } else {
            VectorXd _x = Xsig_aug.col(i);
            VectorXd ad = VectorXd(5);
            ad[0] = _x[2] * sin(_x[3]) * delta_t + 0.5 * cos(_x[3]) * delta_t * delta_t * _x[5];
            ad[1] = _x[2] * -1 * cos(_x[3]) * delta_t + 0.5 * sin(_x[3]) * delta_t * delta_t * _x[5];
            ad[2] = delta_t * _x[5];
            ad[3] = _x[4] * delta_t + 0.5 * delta_t * delta_t * _x[6];
            ad[4] = delta_t * _x[6];
            Xsig_pred.col(i) = Xsig_aug.col(i).head(n_x_) + ad;
        }

    }

    Xsig_pred_ = Xsig_pred;


    // 3. Predict mean and Covariance
    // PredictMeanAndCovariance(Xsig_pred_);
    //create vector for weights
    MatrixXd P = MatrixXd(n_x_, n_x_);
    P_.fill(0);
    //create vector for weights
    VectorXd weights = VectorXd(2*n_aug_+1);
    weights.fill(0);

    //create vector for predicted state
    VectorXd x = VectorXd(n_x_);
    x_.fill(0);

    for (int i = 0; i < 2*n_aug_+1; i++){
        if(i == 0){
            weights[i] = lambda_/(n_aug_+lambda_);
        }
        else{
            weights[i] = 1/(2*(n_aug_+lambda_));
        }

        x_ += weights[i]*Xsig_pred.col(i);

    }
    for(int i = 0; i < 2*n_aug_+1; i ++){
        VectorXd _x = VectorXd(5);
        _x = Xsig_pred.col(i)-x_;
        while (_x[3] > PI) _x[3] -= 2. * PI;
        while (_x[3] < -PI) _x[3] += 2. * PI;
        P_ += weights[i]*_x *_x .transpose();
    }

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
    VectorXd z = meas_package.raw_measurements_;

    MatrixXd R_lidar_ = MatrixXd(2, 2);
    R_lidar_ << std_laspx_*std_laspx_,0,
            0,std_laspy_*std_laspy_;

    VectorXd z_pred = H_laser_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_laser_.transpose();
    MatrixXd S = H_laser_ * P_ * Ht + R_lidar_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    while (x_[3] > PI) x_[3] -= 2. * PI;
    while (x_[3] < -PI) x_[3] += 2. * PI;
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_laser_) * P_;


}


void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */
    int n_z = 3;
//    Xsig_pred_ <<
//              5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
//            1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
//            2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
//            0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
//            0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;


    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0);
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0);


    MatrixXd R = MatrixXd(3,3);
    R.fill(0);
    R(0,0) = std_radr_*std_radr_;
    R(1,1) = std_radphi_*std_radphi_;
    R(2,2) = std_radrd_*std_radrd_;

//    double lambda = 3 - n_aug_;
    VectorXd weights = VectorXd(2*n_aug_+1);
    weights.fill(0);
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {
        double weight = 0.5/(n_aug_+lambda_);
        weights(i) = weight;
    }

    //transform sigma points into measurement space
    for(int i =0; i < 2 * n_aug_ + 1 ; i++){
        VectorXd _z = VectorXd(n_z);
        VectorXd _x = Xsig_pred_.col(i);
        double _px = _x[0];
        double _py = _x[1];
        double _v  = _x[2];
        double _yaw = _x[3];
        double _vx = cos(_yaw)*_v;
        double _vy = sin(_yaw)*_v;
        _z[0] = sqrt(_px*_px + _py*_py);
        _z[1] = atan2(_py,_px);
        _z[2] = (_px*cos(_yaw)*_v + _py*sin(_yaw)*_v)/_z[0];
        Zsig.col(i) = _z;
        z_pred += weights[i]*_z;

    }

    for(int i = 0 ; i < 2 * n_aug_ + 1 ; i++ ) {
        VectorXd y = VectorXd(n_z);
        y = Zsig.col(i).array() - z_pred.array();
        while (y[1] > PI) y[1] -= 2. * PI;
        while (y[1] < -PI) y[1] += 2. * PI;

        S += weights[i]*y*y.transpose() ;
    }
    S+= R;



//    PredictRadarMeasurement(&z_pred,&S);
    VectorXd z = meas_package.raw_measurements_;
    VectorXd diff = VectorXd(n_z);
    diff = z_pred.array() - z.array();
    while(diff[1] > PI)
    {
        diff[1] = diff[1] - 2*PI;
    }
    while(diff[1] < -1*PI)
    {
        diff[1] = diff[1] + 2*PI;
    }

    cout<<"-------DIFF in radar---------"<<endl;
    cout<<diff<<endl;
    cout<<"----------------"<<endl;
    double nis = diff.transpose()*S.inverse()*diff;

    cout<<"NIS"<<nis<<endl;
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0);

    //calculate cross correlation matrix
    for(int i =0; i < 2*n_aug_ +1 ; i++){
        VectorXd y = VectorXd(n_z);
        y = Zsig.col(i) - z_pred;
        while(y[1] > PI)
        {
            y[1] = y[1] - 2*PI;
        }
        while(y[1] < -1*PI)
        {
            y[1] = y[1] + 2*PI;
        }
        VectorXd x = VectorXd(5);
        x = Xsig_pred_.col(i) - x_;

        while(x[3] > PI)
        {
            x[3] = x[3] - 2*PI;
        }
        while(x[3] < -1*PI)
        {
            x[3] = x[3] + 2*PI;
        }

        Tc += weights[i]*(x)*(y).transpose();
    }
    //calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    //update state mean and covariance matrix
    VectorXd y = VectorXd(n_z);
    y = z - z_pred;
    while(y[1] > PI)
    {
        y[1] = y[1] - 2*PI;
    }
    while(y[1] < -1*PI)
    {
        y[1] = y[1] + 2*PI;
    }
    x_ = x_ + K*(y);
    P_ = P_ - K*S*K.transpose();


}
