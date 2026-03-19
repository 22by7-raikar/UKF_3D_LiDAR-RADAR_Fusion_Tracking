#include "ukf.h"

#include <cmath>

#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

constexpr double kEpsilon = 1e-6;

// Keep all heading-related values in [-pi, pi] so the filter does not treat
// equivalent orientations as large residuals.
double NormalizeAngle(double angle) {
  while (angle > M_PI) {
    angle -= 2.0 * M_PI;
  }
  while (angle < -M_PI) {
    angle += 2.0 * M_PI;
  }
  return angle;
}

}  // namespace

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
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.setIdentity();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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
  // Start in the uninitialized state and wait for the first sensor reading to
  // anchor the state estimate.
  is_initialized_ = false;
  time_us_ = 0;

  // Standard CTRV UKF dimensions: 5 state variables plus 2 process-noise
  // variables for longitudinal acceleration and yaw acceleration.
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // Store the predicted sigma points so the radar update can reuse them.
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  // Precompute sigma-point weights once because they do not change across
  // prediction/update cycles.
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size(); ++i) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  // Use the first measurement only to initialize position. Velocity, yaw, and
  // yaw rate are not directly observable from a single measurement.
  if (!is_initialized_) {
    x_.fill(0.0);
    P_.setIdentity();

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Lidar already measures Cartesian position directly.
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);

      // The first lidar reading gives a better estimate of position than of
      // speed or turn rate, so start with tighter covariance on px/py only.
      P_(0, 0) = std_laspx_ * std_laspx_;
      P_(1, 1) = std_laspy_ * std_laspy_;
    } else {
      // Radar gives polar coordinates, so convert to Cartesian position.
      const double rho = meas_package.raw_measurements_(0);
      const double phi = meas_package.raw_measurements_(1);
      const double rho_dot = meas_package.raw_measurements_(2);

      x_(0) = rho * std::cos(phi);
      x_(1) = rho * std::sin(phi);
      x_(2) = std::fabs(rho_dot);

      P_(0, 0) = std_radr_ * std_radr_;
      P_(1, 1) = std_radr_ * std_radr_;
    }

    // Avoid exact zeros near the origin so later divisions stay numerically
    // well behaved during radar projection.
    if (std::fabs(x_(0)) < kEpsilon) {
      x_(0) = kEpsilon;
    }
    if (std::fabs(x_(1)) < kEpsilon) {
      x_(1) = kEpsilon;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  const double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Predict the CTRV state forward to the current measurement time before the
  // sensor-specific correction step. A zero-delta prediction is still useful
  // because it generates a valid sigma-point set for an immediate second
  // sensor update at the same timestamp.
  if (delta_t >= 0.0) {
    Prediction(delta_t);
  }

  // Apply the matching measurement model.
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR &&
             use_radar_) {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  // Augment the current state with process noise so the sigma points can model
  // uncertainty from acceleration and turn-rate changes.
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  // Square-root decomposition produces the spread directions for sigma-point
  // generation.
  MatrixXd L = P_aug.llt().matrixL();

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;

  // Create the symmetric sigma-point set around the augmented mean.
  const double sigma_scale = std::sqrt(lambda_ + n_aug_);
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug.col(i + 1) = x_aug + sigma_scale * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sigma_scale * L.col(i);
  }

  // Propagate each sigma point through the CTRV motion model and inject the
  // process noise terms.
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    const double p_x = Xsig_aug(0, i);
    const double p_y = Xsig_aug(1, i);
    const double v = Xsig_aug(2, i);
    const double yaw = Xsig_aug(3, i);
    const double yawd = Xsig_aug(4, i);
    const double nu_a = Xsig_aug(5, i);
    const double nu_yawdd = Xsig_aug(6, i);

    double px_p;
    double py_p;

    if (std::fabs(yawd) > kEpsilon) {
      px_p = p_x + v / yawd * (std::sin(yaw + yawd * delta_t) - std::sin(yaw));
      py_p = p_y + v / yawd * (-std::cos(yaw + yawd * delta_t) + std::cos(yaw));
    } else {
      px_p = p_x + v * delta_t * std::cos(yaw);
      py_p = p_y + v * delta_t * std::sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    px_p += 0.5 * nu_a * delta_t * delta_t * std::cos(yaw);
    py_p += 0.5 * nu_a * delta_t * delta_t * std::sin(yaw);
    v_p += nu_a * delta_t;
    yaw_p += 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p += nu_yawdd * delta_t;

    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // Recover the predicted mean state from the weighted sigma points.
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // Recover the predicted covariance, normalizing yaw residuals so wrapping at
  // +/-pi does not inflate the uncertainty.
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  const int n_z = 2;

  // Lidar is linear in the state space, so a standard Kalman update is enough.
  MatrixXd H = MatrixXd(n_z, n_x_);
  H.fill(0.0);
  H(0, 0) = 1.0;
  H(1, 1) = 1.0;

  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0.0,
       0.0, std_laspy_ * std_laspy_;

  const VectorXd z = meas_package.raw_measurements_;
  const VectorXd z_pred = H * x_;
  const VectorXd y = z - z_pred;
  const MatrixXd Ht = H.transpose();
  const MatrixXd S = H * P_ * Ht + R;
  const MatrixXd K = P_ * Ht * S.inverse();

  // Correct the predicted state with the position residual.
  x_ += K * y;

  const MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  const int n_z = 3;

  // Transform predicted sigma points into radar measurement space
  // [range, bearing, range_rate].
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    const double p_x = Xsig_pred_(0, i);
    const double p_y = Xsig_pred_(1, i);
    const double v = Xsig_pred_(2, i);
    const double yaw = Xsig_pred_(3, i);

    const double v1 = std::cos(yaw) * v;
    const double v2 = std::sin(yaw) * v;
    const double rho = std::sqrt(p_x * p_x + p_y * p_y);
    const double safe_rho = rho < kEpsilon ? kEpsilon : rho;

    Zsig(0, i) = safe_rho;
    Zsig(1, i) = std::atan2(p_y, p_x);
    Zsig(2, i) = (p_x * v1 + p_y * v2) / safe_rho;
  }

  // Predicted radar measurement mean.
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Innovation covariance in radar space, again normalizing the angle term.
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0.0, 0.0,
       0.0, std_radphi_ * std_radphi_, 0.0,
       0.0, 0.0, std_radrd_ * std_radrd_;
  S += R;

  // Cross-correlation between state space and radar measurement space gives us
  // the Kalman gain for the nonlinear update.
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  const MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  z_diff(1) = NormalizeAngle(z_diff(1));

  // Apply the nonlinear measurement correction.
  x_ += K * z_diff;
  x_(3) = NormalizeAngle(x_(3));
  P_ -= K * S * K.transpose();
}