// One sensor reading passed to the filter.
#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include <Eigen/Dense>

struct Measurement {
  enum class Sensor { Lidar, Radar };

  Sensor sensor = Sensor::Lidar;
  long long timestamp_us = 0;  // time of the reading [us]
  Eigen::VectorXd z;           // lidar: [px, py]; radar: [rho, phi, rho_dot]
};

#endif  // MEASUREMENT_H
