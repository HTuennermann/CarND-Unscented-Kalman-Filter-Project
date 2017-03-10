#include <iostream>
#include "ukf.h"
using namespace std;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  VectorXd diag(5);
  diag << 1,1,1000,M_PI_2,0.1; // we need some confidence in yawdd otherwise we have a divergence porblem
  P_ = diag.asDiagonal();
  cout << P_ << endl;

  n_x_ = 5;
  n_aug_ = 7;
  n_sig_ = 2*n_aug_+1;
  lambda_ = 3 - n_aug_;


  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  // Process noise standard deviation longitudinal acceleration in m/s^2

  std_a_ = 3; // needs change

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5; // was way to big




  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.085;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.085;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.5;

  // Weight initialization
  weights_ = VectorXd(n_sig_);
  weights_(0) =lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    weights_(i) = 0.5/(n_aug_+lambda_);
  }


  is_initialized_ = false;
  use_laser_ = true;
  use_radar_ = true;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  // Drop packages if not necessary
  if(meas_package.sensor_type_==MeasurementPackage::LASER && use_laser_==false)
    return;
  if(meas_package.sensor_type_==MeasurementPackage::RADAR && use_radar_==false)
    return;


  if(! is_initialized_) {
    if(meas_package.sensor_type_==MeasurementPackage::LASER) {
      cout << "laser" << endl;

      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      time_us_ = meas_package.timestamp_;
      cout <<   meas_package.raw_measurements_[0] << " " << meas_package.raw_measurements_[1] << endl;
    }
    else {
      cout << "radar" << endl;

      double  meas_rho =  meas_package.raw_measurements_[0];

      double meas_phi =  meas_package.raw_measurements_[1];
      cout <<   meas_rho * cos(meas_phi) << " " << meas_rho * sin(meas_phi) << endl;
      x_ << meas_rho * cos(meas_phi), meas_rho * sin(meas_phi), 0, 0, 0;
      time_us_ = meas_package.timestamp_;

    }
    is_initialized_=true;
  }
  else {


    double delta_t = (meas_package.timestamp_ - time_us_) /1000000.0;
    time_us_= meas_package.timestamp_;



    if(meas_package.sensor_type_==MeasurementPackage::LASER){
      Prediction(delta_t);
      UpdateLidar(meas_package);
    }
    else {
      Prediction(delta_t);
      UpdateRadar(meas_package);
    }


  }


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // Sigma points

  VectorXd x_aug(n_aug_);
  x_aug.fill(0);
  x_aug.head(n_x_) = x_;


  //MatrixXd(n_aug_, 2 * n_aug_ + 1);


  MatrixXd P_aug(n_aug_, n_aug_);
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;


  MatrixXd A = P_aug.llt().matrixL();


  double s = sqrt(lambda_ + n_aug_);


  MatrixXd Xsig_aug(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i<n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + s * A.col(i);
    Xsig_aug.col(n_aug_+i+1) = x_aug - s * A.col(i);
  }


  // sigma point prediction


  for(int i = 0; i<n_sig_; i++) {

    double p_x = Xsig_aug.col(i)(0);
    double p_y = Xsig_aug.col(i)(1);
    double v = Xsig_aug.col(i)(2);
    double yaw = Xsig_aug.col(i)(3);
    double yawd = Xsig_aug.col(i)(4);
    double nu_a = Xsig_aug.col(i)(5);
    double nu_yawdd = Xsig_aug.col(i)(6);

    VectorXd add(n_x_);
    if(abs(yawd) > 0.0001) {
    add << v / yawd * (sin(yaw+yawd*delta_t)-sin(yaw)) + 0.5*delta_t*delta_t*cos(yaw)*nu_a,
    v / yawd * (-cos(yaw+yawd*delta_t)+cos(yaw)) + 0.5*delta_t*delta_t*sin(yaw)*nu_a,
    0+ nu_a*delta_t,
    yawd*delta_t+0.5*delta_t*delta_t*nu_yawdd,
    0+delta_t*nu_yawdd;
    }
    else {
      add << v*delta_t*cos(yaw) + 0.5*delta_t*delta_t*cos(yaw)*nu_a,
      v*delta_t*sin(yaw) + 0.5*delta_t*delta_t*sin(yaw)*nu_a,
      0+ nu_a*delta_t,
      yawd*delta_t+0.5*delta_t*delta_t*nu_yawdd,
      0+delta_t*nu_yawdd;
    }

    Xsig_pred_.col(i) = Xsig_aug.col(i).head(5) + add;

  }


  //predicted state mean

  x_ = Xsig_pred_ * weights_;




  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization



    x_diff(3) = x_diff(3)-ceil((x_diff(3)-M_PI)/(2.*M_PI))*2.*M_PI;




    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();


  }



}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);

  // Just copy over the rows
  Zsig.row(0) = Xsig_pred_.row(0);
  Zsig.row(1) = Xsig_pred_.row(1);

  z_pred = Zsig * weights_;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0);
  for(int i = 0; i<n_sig_; i++) {
    VectorXd temp = Zsig.col(i) - z_pred;
    S = S + weights_(i)*temp*temp.transpose();
  }

  MatrixXd R(2,2);
  R << std_laspx_*std_laspx_,0,
          0,std_laspy_*std_laspy_;
  S = S + R;


  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0);
  for(int i=0; i<n_sig_; i++) {
    VectorXd tx = Xsig_pred_.col(i) - x_;
    VectorXd tz = Zsig.col(i) - z_pred;
    Tc = Tc + weights_(i)*tx*tz.transpose();
  }



  MatrixXd K = Tc*(S.inverse());


  VectorXd xo = x_;

  x_ = x_ + K*(meas_package.raw_measurements_-z_pred);
  NIS_laser_ = ((meas_package.raw_measurements_-z_pred).transpose())*S.inverse()*(meas_package.raw_measurements_-z_pred);

  P_ = P_ - K*S*K.transpose();


}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);


  for(int i = 0; i<n_sig_; i++) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yawd = Xsig_pred_(4,i);
    double rho = sqrt(px*px + py*py);
    double phi = atan2(py,px);
    double rhod;

    if(rho>0.01)
      rhod = (px*cos(yaw)*v+py*sin(yaw)*v)/rho;
    else
      rhod = 0;
    Zsig.col(i)<< rho, phi, rhod;

  }
  z_pred = Zsig * weights_;
  S.fill(0);
  for(int i = 0; i<n_sig_; i++) {
    VectorXd temp = Zsig.col(i) - z_pred;
    S = S + weights_(i)*temp*temp.transpose();
  }

  MatrixXd R(3,3);
  R << std_radr_*std_radr_,0,0,
  0,std_radphi_*std_radphi_,0,
  0,0,std_radrd_*std_radrd_;
  S = S + R;
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);


  Tc.fill(0);
  for(int i=0; i<n_sig_; i++) {
    VectorXd tx = Xsig_pred_.col(i) - x_;
    VectorXd tz = Zsig.col(i) - z_pred;
    Tc = Tc + weights_(i)*tx*tz.transpose();
  }



  MatrixXd K = Tc*(S.inverse());

  x_ = x_ + K*(meas_package.raw_measurements_-z_pred);
  NIS_radar_ = ((meas_package.raw_measurements_-z_pred).transpose())*S.inverse()*(meas_package.raw_measurements_-z_pred);


  P_ = P_ - K*S*K.transpose();


}
