#include "ukf.h"

#include <cmath>

#include <Eigen/Dense>

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

UKF::UKF() {
  x_ = VectorXd::Zero(5);
  P_ = MatrixXd::Identity(5, 5);

  // Process noise, tuned for the highway traffic.
  accel_noise_ = 1.5;
  yaw_accel_noise_ = 0.5;

  // Sensor noise, as specified for the lidar and radar.
  lidar_std_x_ = 0.15;
  lidar_std_y_ = 0.15;
  radar_range_std_ = 0.3;
  radar_bearing_std_ = 0.03;
  radar_rate_std_ = 0.3;

  // CTRV: 5 state variables, plus longitudinal and yaw acceleration noise.
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // The radar update reuses the sigma points of the last prediction.
  sigma_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size(); ++i) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }
}

void UKF::processMeasurement(const Measurement& m) {
  // Use the first measurement only to initialize position. Velocity, yaw, and
  // yaw rate are not directly observable from a single measurement.
  if (!is_initialized_) {
    x_.fill(0.0);
    P_.setIdentity();

    if (m.sensor == Measurement::Sensor::Lidar) {
      // Lidar already measures Cartesian position directly.
      x_(0) = m.z(0);
      x_(1) = m.z(1);

      // The first lidar reading gives a better estimate of position than of
      // speed or turn rate, so start with tighter covariance on px/py only.
      P_(0, 0) = lidar_std_x_ * lidar_std_x_;
      P_(1, 1) = lidar_std_y_ * lidar_std_y_;
    } else {
      // Radar gives polar coordinates, so convert to Cartesian position.
      const double rho = m.z(0);
      const double phi = m.z(1);
      const double rho_dot = m.z(2);

      x_(0) = rho * std::cos(phi);
      x_(1) = rho * std::sin(phi);
      x_(2) = std::fabs(rho_dot);

      P_(0, 0) = radar_range_std_ * radar_range_std_;
      P_(1, 1) = radar_range_std_ * radar_range_std_;
    }

    // Avoid exact zeros near the origin so later divisions stay numerically
    // well behaved during radar projection.
    if (std::fabs(x_(0)) < kEpsilon) {
      x_(0) = kEpsilon;
    }
    if (std::fabs(x_(1)) < kEpsilon) {
      x_(1) = kEpsilon;
    }

    time_us_ = m.timestamp_us;
    is_initialized_ = true;
    return;
  }

  const double delta_t = (m.timestamp_us - time_us_) / 1000000.0;
  time_us_ = m.timestamp_us;

  // Predict the CTRV state forward to the current measurement time before the
  // sensor-specific correction step. A zero-delta prediction is still useful
  // because it generates a valid sigma-point set for an immediate second
  // sensor update at the same timestamp.
  if (delta_t >= 0.0) {
    predict(delta_t);
  }

  // Apply the matching measurement model.
  if (m.sensor == Measurement::Sensor::Lidar && use_lidar_) {
    updateLidar(m);
  } else if (m.sensor == Measurement::Sensor::Radar &&
             use_radar_) {
    updateRadar(m);
  }
}

void UKF::predict(double delta_t) {
  // Augment the current state with process noise so the sigma points can model
  // uncertainty from acceleration and turn-rate changes.
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = accel_noise_ * accel_noise_;
  P_aug(n_x_ + 1, n_x_ + 1) = yaw_accel_noise_ * yaw_accel_noise_;

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

    sigma_pred_(0, i) = px_p;
    sigma_pred_(1, i) = py_p;
    sigma_pred_(2, i) = v_p;
    sigma_pred_(3, i) = yaw_p;
    sigma_pred_(4, i) = yawd_p;
  }

  // Recover the predicted mean state from the weighted sigma points.
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    x_ += weights_(i) * sigma_pred_.col(i);
  }

  // Recover the predicted covariance, normalizing yaw residuals so wrapping at
  // +/-pi does not inflate the uncertainty.
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_diff = sigma_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::updateLidar(const Measurement& m) {
  const int n_z = 2;

  // Lidar is linear in the state space, so a standard Kalman update is enough.
  MatrixXd H = MatrixXd(n_z, n_x_);
  H.fill(0.0);
  H(0, 0) = 1.0;
  H(1, 1) = 1.0;

  MatrixXd R = MatrixXd(n_z, n_z);
  R << lidar_std_x_ * lidar_std_x_, 0.0,
       0.0, lidar_std_y_ * lidar_std_y_;

  const VectorXd z = m.z;
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

void UKF::updateRadar(const Measurement& m) {
  const int n_z = 3;

  // Transform predicted sigma points into radar measurement space
  // [range, bearing, range_rate].
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    const double p_x = sigma_pred_(0, i);
    const double p_y = sigma_pred_(1, i);
    const double v = sigma_pred_(2, i);
    const double yaw = sigma_pred_(3, i);

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
  R << radar_range_std_ * radar_range_std_, 0.0, 0.0,
       0.0, radar_bearing_std_ * radar_bearing_std_, 0.0,
       0.0, 0.0, radar_rate_std_ * radar_rate_std_;
  S += R;

  // Cross-correlation between state space and radar measurement space gives us
  // the Kalman gain for the nonlinear update.
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));

    VectorXd x_diff = sigma_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  const MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = m.z - z_pred;
  z_diff(1) = NormalizeAngle(z_diff(1));

  // Apply the nonlinear measurement correction.
  x_ += K * z_diff;
  x_(3) = NormalizeAngle(x_(3));
  P_ -= K * S * K.transpose();
}