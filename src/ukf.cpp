#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF(double std_a, double std_yaw) {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;

  n_aug_ = n_x_ + 2;

  n_sig_ = (2 * n_aug_) + 1;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = std_a;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = std_yaw;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
	P_ << 1,0,0,0,0,
        0,1,0,0,0,
        0,0,1,0,0,
        0,0,0,1,0,
        0,0,0,0,1;

  n_x_ = 5;
  lambda_ = 5.5 - n_aug_;
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
  laser_count_ = 0;
  radar_count_ = 0;
  weights_ = VectorXd(n_sig_);
  InitializeWeights();

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  std::cout<<"Process Measurements"<<std::endl;
  if(false == is_initialized_)
  {
    std::cout<<"First measurement"<<std::endl;
    x_ = VectorXd(5);
    if( MeasurementPackage::LASER == meas_package.sensor_type_)
    {
      std::cout<<"First measurement is LASER"<<std::endl;
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
      x_(2) = 0;
      x_(3) = 0;//atan2(x_(1),x_(0));
      x_(4) = 0;
    }
    else
    {
      std::cout<<"First measurement is RADAR"<<std::endl;
      x_(0) = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
      x_(1) = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));
      x_(2) = meas_package.raw_measurements_(2);
      x_(3) = meas_package.raw_measurements_(1);
      x_(4) = 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    std::cout<<"First x : "<<std::endl;
    std::cout<<x_;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<"First P : "<<std::endl;
    std::cout<<P_;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<"End of Measurement processing.\n\n\n\n\n\n\n\n";
    return;
  }
  double dt = ((double)meas_package.timestamp_ - (double)time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  std::cout<<"Going to perform prediction, time elapsed : "<<dt<<std::endl;
  Prediction(dt);

  if(MeasurementPackage::LASER == meas_package.sensor_type_)
  {
    std::cout<<"LASER update"<<std::endl;
    UpdateLidar(meas_package);
  }
  else
  {
    std::cout<<"RADAR update"<<std::endl;
    UpdateRadar(meas_package);
  }
  std::cout<<"End of Measurement processing.\n\n\n\n\n\n\n\n";

}


void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
  /* 1. Generate Sigma points */
  /* 2. Augment Sigma points with noise */
  std::cout<<"Calculating augmented sigma points"<<std::endl;
  MatrixXd x_aug = AugmentSigmaPoints();
  std::cout<<"Augmented X :"<<std::endl;
  std::cout<<x_aug;
  std::cout<<std::endl;
  std::cout<<std::endl;
  /* 3. Augment Covariance Matrix with noise */
  /* 4. Predict state updated sigma points */
  std::cout<<"Predicting updated sigma points"<<std::endl;

  PredictSigmaPoints(x_aug, delta_t);
  std::cout<<Xsig_pred_;
  std::cout<<std::endl;
  std::cout<<std::endl;
  /* 5. Get state vector and covariance from predicted sigma points */
  std::cout<<"Performing state update from predicted sigma points"<<std::endl;
  GetStateFromSigmaPoints();
  std::cout<<"State : "<<std::endl;
  std::cout<<x_;
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<"P Matrix : "<<std::endl;
  std::cout<<P_;
  std::cout<<std::endl;
  std::cout<<std::endl;
}
MatrixXd UKF::GetSigmaPoints() {
  MatrixXd sigmaPoints(n_x_, (2 * n_x_) + 1);
  lambda_ = 3 - n_x_;

  MatrixXd sqrtP = P_.llt().matrixL();

  sqrtP = sqrt(lambda_ + n_x_) * sqrtP;

  sigmaPoints.col(0) = x_;
  for(int i = 1; i < n_x_ + 1; i++)
  {
    sigmaPoints.col(i) = x_ + sqrtP.col(i-1);
    sigmaPoints.col(i + n_x_) = x_ - sqrtP.col(i-1);
  }
  return sigmaPoints;
}

MatrixXd UKF::AugmentSigmaPoints() {
  /* Local variable declarations */
  MatrixXd augSigmaPoints(n_aug_, n_sig_);
  VectorXd x_aug(n_aug_);

  MatrixXd P_aug(n_aug_, n_aug_);
  MatrixXd Q(2,2);

  /* Build the noise matrix */
  Q.fill(0.0);
  Q(0,0) = std_a_ * std_a_;
  Q(1,1) = std_yawdd_ * std_yawdd_;

  /* Get augmented Covariance matrix */
  P_aug.fill(0.0);
  P_aug.block(0,0,n_x_,n_x_) = P_;
  P_aug.block(n_x_,n_x_,2,2) = Q;

  /* Get Square root of covariance matrix */
  MatrixXd sqrtP = P_aug.llt().matrixL();
  sqrtP = sqrt(lambda_ + n_aug_) * sqrtP;

  /* Create Augmented state matrix */
  x_aug.segment(0,n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;

  augSigmaPoints.col(0) = x_aug;

  for(int i = 0; i < n_aug_; i++)
  {
    augSigmaPoints.col(i + 1) = x_aug + sqrtP.col(i);
    augSigmaPoints.col(i + n_aug_+ 1) = x_aug - sqrtP.col(i);
  }

  return augSigmaPoints;

}


void UKF::PredictSigmaPoints(MatrixXd x_aug, double dt) {

  for(int i = 0; i < n_sig_; i++)
  {
    double state_update_x;
    double state_update_y;
    double state_update_v;
    double state_update_rho;
    double state_update_rho_d;
    // State update Equations
    if(fabs(x_aug(4, i)) >= 0.001)
    {
      state_update_x     = x_aug(2 ,i)/x_aug(4,i) * ( sin(x_aug(3, i) + x_aug(4, i) * dt ) - sin(x_aug(3, i)) );
      state_update_y     = x_aug(2, i)/x_aug(4,i) * ( cos(x_aug(3, i)) - cos(x_aug(3, i) + x_aug(4, i) * dt ) );
      state_update_v     = 0;
      state_update_rho   = x_aug(4, i) * dt;
      state_update_rho_d = 0;
    }
    else
    {
      state_update_x     = x_aug(2, i) * cos(x_aug(3, i)) * dt;
      state_update_y     = x_aug(2, i) * sin(x_aug(3, i)) * dt;
      state_update_v     = 0;
      state_update_rho   = x_aug(4, i) * dt;
      state_update_rho_d = 0;
    }

    // Noise Update
    state_update_x     += 0.5 * (dt * dt) * cos(x_aug(3, i)) * x_aug(5, i);
    state_update_y     += 0.5 * (dt * dt) * sin(x_aug(3, i)) * x_aug(5, i);
    state_update_v     += dt * x_aug(5, i);
    state_update_rho   += 0.5 * (dt * dt) * x_aug(6, i);
    state_update_rho_d += dt * x_aug(6, i);
    // Final update to class variable 
    Xsig_pred_.col(i)(0) = x_aug(0,i) + state_update_x;
    Xsig_pred_.col(i)(1) = x_aug(1,i) + state_update_y;
    Xsig_pred_.col(i)(2) = x_aug(2,i) + state_update_v;
    Xsig_pred_.col(i)(3) = (x_aug(3,i) + state_update_rho);
    Xsig_pred_.col(i)(4) = x_aug(4,i) + state_update_rho_d;
#ifdef DEBUG
    if(Xsig_pred_.col(i)(3) - NormalizeAngle(x_aug(3,i) + state_update_rho) > 1e-8)
    {
     std::cout<<std::endl<<"Predicted Sigma psi : "<<Xsig_pred_.col(i)(3)<<std::endl;
     std::cout<<"Un-normalized Sigma psi : "<<NormalizeAngle(x_aug(3,i) + state_update_rho)<<std::endl;
    }
#endif
  }
}

void UKF::GetStateFromSigmaPoints()
{
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  x_(3) = NormalizeAngle(x_(3));

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 1; i < n_sig_; i++) {  //iterate over sigma points

    // state difference
    /* Trick to ensure Covariance matrix remains positive definite, since lambda is negative
     * this is not guaranteed without this trick. This trick causes some loss of accuracy.
     */
    VectorXd x_diff = Xsig_pred_.col(i) - Xsig_pred_.col(0);

    x_diff(3) = NormalizeAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  } 
}

MatrixXd UKF::LidarGetZSigPoints(int n_z)
{
  MatrixXd Zsig(n_z, n_sig_);
  for(int i = 0; i < n_sig_; i++)
  {
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
 
    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;

  }
  return Zsig;

}


void  UKF::InitializeWeights()
{
  VectorXd weights = VectorXd(n_sig_);
  double weight_0 = lambda_/(lambda_ + n_aug_);

  weights_(0) = weight_0;

  for (int i=1; i<n_sig_; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

}

VectorXd UKF::GetZPred(int n_z, MatrixXd Zsig, bool isRadar)
{
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_sig_; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  if(true == isRadar)
  {
    z_pred(1) = NormalizeAngle(z_pred(1));
  }
  return z_pred;
}

void UKF::GetSMatrix(int n_z, MatrixXd R, MatrixXd Zsig, VectorXd z_pred, bool isRadar, MatrixXd &S)
{
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) 
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    if(true == isRadar)
    {
      //angle normalization
      z_diff(1) = NormalizeAngle(z_diff(1));
    }

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R;
}


MatrixXd UKF::GetTMatrix(int n_z, MatrixXd Zsig, VectorXd z_pred, bool isRadar)
{
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    if(true == isRadar)
    {
      //angle normalization
      z_diff(1) = NormalizeAngle(z_diff(1));
    }
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    if(true == isRadar)
    {
      //angle normalization
      x_diff(3) = NormalizeAngle(x_diff(3));
    }

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  return Tc;

}

double UKF::NormalizeAngle(double input)
{
#if 0
  while(input >  M_PI) input -= (2.0 * M_PI);
  while(input < -M_PI) input += (2.0 * M_PI);
#else
  if(input < M_PI && input > -M_PI)
    return input;

  input = fmod(input, 2 * M_PI);

  if(input < -M_PI)
  {
    input += 2 * M_PI;
  }
  else if(input > M_PI)
  {
    input -= 2 * M_PI;
  }
#endif
  return input;
}
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z = 2;

#ifdef LIDAR_NOT_USING_UKF
  std::cout<<"LIDAR : Measurement points : "<<std::endl;
  std::cout<<meas_package.raw_measurements_;
  std::cout<<std::endl;
  std::cout<<std::endl;

  MatrixXd H(2,5);

  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;

  VectorXd z = meas_package.raw_measurements_;

  VectorXd z_diff = z - (H * x_);
  std::cout<<"LIDAR : z_diff = "<<std::endl;
  std::cout<<z_diff;
  std::cout<<std::endl;
  std::cout<<std::endl;

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_ * std_laspx_, 0                      , 
          0                      , std_laspy_ * std_laspy_;


  MatrixXd S = (H * P_ * H.transpose()) + R;
  std::cout<<"LIDAR : S = "<<std::endl;
  std::cout<<S;
  std::cout<<std::endl;
  std::cout<<std::endl;

  MatrixXd K = P_ * H.transpose() * S.inverse();
  std::cout<<"LIDAR : K = "<<std::endl;
  std::cout<<K;
  std::cout<<std::endl;
  std::cout<<std::endl;

  x_ = x_ + K * z_diff;
  x_(3) = NormalizeAngle(x_(3));

  P_ = (MatrixXd::Identity(5,5) - (K * H)) * P_;

  std::cout<<"LIDAR : x after update : "<<std::endl;
  std::cout<<x_;
  std::cout<<std::endl;
  std::cout<<std::endl;

  //P_ = P_ - K*S*K.transpose();
  std::cout<<"LIDAR : P after update : "<<std::endl;
  std::cout<<P_;
  std::cout<<std::endl;
  std::cout<<std::endl;
#else
  MatrixXd Zsig = LidarGetZSigPoints(n_z);
  std::cout<<"LIDAR : Zsig = "<<std::endl;
  std::cout<<Zsig;
  std::cout<<std::endl;
  std::cout<<std::endl;

  VectorXd z_pred = GetZPred(n_z, Zsig, false);
  std::cout<<"LIDAR : z_pred = "<<std::endl;
  std::cout<<z_pred;
  std::cout<<std::endl;
  std::cout<<std::endl;

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_ * std_laspx_, 0, 
          0                      , std_laspy_ * std_laspy_;

  MatrixXd S(n_z, n_z);
  GetSMatrix(n_z, R, Zsig, z_pred, false, S);
  std::cout<<"LIDAR : S = "<<std::endl;
  std::cout<<S;
  std::cout<<std::endl;
  std::cout<<std::endl;

  //create matrix for cross correlation Tc
  MatrixXd Tc = GetTMatrix(n_z, Zsig, z_pred, false);

  std::cout<<"LIDAR : T = "<<std::endl;
  std::cout<<Tc;
  std::cout<<std::endl;
  std::cout<<std::endl;
  
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  std::cout<<"LIDAR : K = "<<std::endl;
  std::cout<<K;
  std::cout<<std::endl;
  std::cout<<std::endl;

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  std::cout<<"LIDAR : z_diff = "<<std::endl;
  std::cout<<z_diff;
  std::cout<<std::endl;
  std::cout<<std::endl;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  //angle normalization
  x_(3) = NormalizeAngle(x_(3));
  std::cout<<"LIDAR : x after update : "<<std::endl;
  std::cout<<x_;
  std::cout<<std::endl;
  std::cout<<std::endl;

  P_ = P_ - K*S*K.transpose();
  std::cout<<"LIDAR : P after update : "<<std::endl;
  std::cout<<P_;
  std::cout<<std::endl;
  std::cout<<std::endl;
#endif

  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  laser_count_++;

  std::cout<<"LIDAR : NIS = "<<NIS_laser_<<std::endl<<std::endl;
}

MatrixXd UKF::RadarGetZSigPoints(int n_z)
{
  MatrixXd Zsig(n_z, n_sig_);

  for(int i = 0; i < n_sig_; i++)
  {
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    if(Zsig(0,i) < 0.00001)
    {
      p_x = 0.00001;
      p_y = 0.00001;
      Zsig(0,i) = 0.00001;
      Zsig(2,i) = 0;                   //r_dot
    }
    else
    {
      Zsig(2,i) = (p_x*v1 + p_y*v2 )/Zsig(0,i);                   //r_dot
    }
    Zsig(1,i) = atan2(p_y,p_x);                               //phi

  }
  return Zsig;

}


void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z = 3;
  //set vector for weights

  std::cout<<"RADAR : Measurement points : "<<std::endl;
  std::cout<<meas_package.raw_measurements_;
  std::cout<<std::endl;
  std::cout<<std::endl;

  MatrixXd Zsig = RadarGetZSigPoints(n_z);
  std::cout<<"RADAR : Zsig = "<<std::endl;
  std::cout<<Zsig;
  std::cout<<std::endl;
  std::cout<<std::endl;

  VectorXd z_pred = GetZPred(n_z, Zsig, true);
  std::cout<<"RADAR : z_pred = "<<std::endl;
  std::cout<<z_pred;
  std::cout<<std::endl;
  std::cout<<std::endl;

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_ * std_radr_, 0                        , 0,
          0                    , std_radphi_ * std_radphi_, 0,
          0                    , 0                        , std_radrd_ * std_radrd_;

  MatrixXd S (n_z, n_z);
  GetSMatrix(n_z, R, Zsig, z_pred, true, S);
  std::cout<<"RADAR : S = "<<std::endl;
  std::cout<<S;
  std::cout<<std::endl;
  std::cout<<std::endl;


  MatrixXd Tc = GetTMatrix(n_z, Zsig, z_pred, true);
  std::cout<<"RADAR : T = "<<std::endl;
  std::cout<<Tc;
  std::cout<<std::endl;
  std::cout<<std::endl;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  std::cout<<"RADAR : K = "<<std::endl;
  std::cout<<K;
  std::cout<<std::endl;
  std::cout<<std::endl;

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  //angle normalization
  z_diff(1) = NormalizeAngle(z_diff(1));
  std::cout<<"RADAR : z_diff = "<<std::endl;
  std::cout<<z_diff;
  std::cout<<std::endl;
  std::cout<<std::endl;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;

  //angle normalization
  x_(3) = NormalizeAngle(x_(3));
  std::cout<<"RADAR : x after update : "<<std::endl;
  std::cout<<x_;
  std::cout<<std::endl;
  std::cout<<std::endl;

  P_ = P_ - K*S*K.transpose();
  std::cout<<"RADAR : P after update : "<<std::endl;
  std::cout<<P_;
  std::cout<<std::endl;
  std::cout<<std::endl;

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  radar_count_++;
  std::cout<<"RADAR : NIS = "<<NIS_radar_<<std::endl<<std::endl;
}
