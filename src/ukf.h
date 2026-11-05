// Unscented Kalman filter with a CTRV (constant turn rate and velocity)
// motion model, fusing lidar and radar measurements of one vehicle.
#ifndef UKF_H
#define UKF_H

#include <Eigen/Dense>

#include "measurement.h"

class UKF {
 public:
  UKF();

  // Initialises on the first measurement; afterwards predicts to the
  // measurement time and applies the matching sensor update.
  void processMeasurement(const Measurement& m);

  // Propagates sigma points, state and covariance by dt seconds.
  void predict(double dt);

  // Linear Kalman update with a lidar position measurement.
  void updateLidar(const Measurement& m);

  // Unscented update with a radar (range, bearing, range rate) measurement.
  void updateRadar(const Measurement& m);

  bool is_initialized_ = false;
  bool use_lidar_ = true;  // if false, lidar is used only to initialise
  bool use_radar_ = true;  // if false, radar is used only to initialise

  Eigen::VectorXd x_;           // state [px, py, v, yaw, yaw_rate] (SI units, rad)
  Eigen::MatrixXd P_;           // state covariance
  Eigen::MatrixXd sigma_pred_;  // predicted sigma points, one per column
  long long time_us_ = 0;       // time of the current state [us]

  // Process noise.
  double accel_noise_;      // longitudinal acceleration std [m/s^2]
  double yaw_accel_noise_;  // yaw acceleration std [rad/s^2]

  // Measurement noise of the sensors.
  double lidar_std_x_;        // [m]
  double lidar_std_y_;        // [m]
  double radar_range_std_;    // [m]
  double radar_bearing_std_;  // [rad]
  double radar_rate_std_;     // [m/s]

  Eigen::VectorXd weights_;  // sigma point weights
  int n_x_;                  // state dimension
  int n_aug_;                // augmented state dimension (state + process noise)
  double lambda_;            // sigma point spreading parameter
};

#endif  // UKF_H
