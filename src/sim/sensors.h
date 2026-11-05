// Simulated lidar and radar measurements of the traffic, track drawing,
// RMSE, and PCD file input/output.
#ifndef SENSORS_H
#define SENSORS_H

#include <string>
#include <vector>

#include <Eigen/Dense>
#include <pcl/io/pcd_io.h>

#include "../render/render.h"

struct LidarMarker {
  double x, y;
};

struct RadarMarker {
  double rho, phi, rho_dot;
};

class Sensors {
 public:
  std::vector<Eigen::VectorXd> estimations;   // [px, py, vx, vy] per tracked car and step
  std::vector<Eigen::VectorXd> ground_truth;  // matching true values

  // Gaussian noise sample; the seed makes runs repeatable.
  double noise(double stddev, long long seed);

  // Measures the car's position with noise, draws the marker, and feeds the car's UKF.
  LidarMarker lidarSense(Car& car, pcl::visualization::PCLVisualizer::Ptr& viewer,
                         long long timestamp, bool visualize);

  // Measures range, bearing and range rate from the ego car, draws them, and
  // feeds the car's UKF.
  RadarMarker radarSense(Car& car, const Car& ego, pcl::visualization::PCLVisualizer::Ptr& viewer,
                         long long timestamp, bool visualize);

  // Draws the UKF estimate, and optionally its predicted path over `time`
  // seconds in `steps` steps.
  void renderTrack(const Car& car, pcl::visualization::PCLVisualizer::Ptr& viewer, double time,
                   int steps);

  // Root-mean-square error of the estimations against the ground truth.
  static Eigen::VectorXd rmse(const std::vector<Eigen::VectorXd>& estimations,
                              const std::vector<Eigen::VectorXd>& ground_truth);

  static void savePcd(const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud, const std::string& file);
  static pcl::PointCloud<pcl::PointXYZ>::Ptr loadPcd(const std::string& file);
};

#endif  // SENSORS_H
