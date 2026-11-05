#include "sensors.h"

#include <cmath>
#include <iostream>
#include <random>

using Eigen::VectorXd;

double Sensors::noise(double stddev, long long seed) {
  std::mt19937 engine(static_cast<std::mt19937::result_type>(seed));
  return std::normal_distribution<double>(0.0, stddev)(engine);
}

LidarMarker Sensors::lidarSense(Car& car, pcl::visualization::PCLVisualizer::Ptr& viewer,
                                long long timestamp, bool visualize) {
  const LidarMarker marker{car.position.x + noise(0.15, timestamp),
                           car.position.y + noise(0.15, timestamp + 1)};
  if (visualize)
    viewer->addSphere(pcl::PointXYZ(marker.x, marker.y, 3.0), 0.5, 1, 0, 0,
                      car.name + "_lidar");

  Measurement m;
  m.sensor = Measurement::Sensor::Lidar;
  m.timestamp_us = timestamp;
  m.z = VectorXd(2);
  m.z << marker.x, marker.y;
  car.ukf.processMeasurement(m);
  return marker;
}

RadarMarker Sensors::radarSense(Car& car, const Car& ego,
                                pcl::visualization::PCLVisualizer::Ptr& viewer,
                                long long timestamp, bool visualize) {
  const double dx = car.position.x - ego.position.x;
  const double dy = car.position.y - ego.position.y;
  const double rho = std::sqrt(dx * dx + dy * dy);
  const double phi = std::atan2(dy, dx);
  const double rho_dot = (car.velocity * std::cos(car.angle) * rho * std::cos(phi) +
                          car.velocity * std::sin(car.angle) * rho * std::sin(phi)) / rho;

  const RadarMarker marker{rho + noise(0.3, timestamp + 2), phi + noise(0.03, timestamp + 3),
                           rho_dot + noise(0.3, timestamp + 4)};
  if (visualize) {
    const double hx = ego.position.x + marker.rho * std::cos(marker.phi);
    const double hy = ego.position.y + marker.rho * std::sin(marker.phi);
    viewer->addLine(pcl::PointXYZ(ego.position.x, ego.position.y, 3.0),
                    pcl::PointXYZ(hx, hy, 3.0), 1, 0, 1, car.name + "_rho");
    viewer->addArrow(pcl::PointXYZ(hx, hy, 3.0),
                     pcl::PointXYZ(hx + marker.rho_dot * std::cos(marker.phi),
                                   hy + marker.rho_dot * std::sin(marker.phi), 3.0),
                     1, 0, 1, car.name + "_rho_dot");
  }

  Measurement m;
  m.sensor = Measurement::Sensor::Radar;
  m.timestamp_us = timestamp;
  m.z = VectorXd(3);
  m.z << marker.rho, marker.phi, marker.rho_dot;
  car.ukf.processMeasurement(m);
  return marker;
}

void Sensors::renderTrack(const Car& car, pcl::visualization::PCLVisualizer::Ptr& viewer,
                          double time, int steps) {
  UKF ukf = car.ukf;
  viewer->addSphere(pcl::PointXYZ(ukf.x_[0], ukf.x_[1], 3.5), 0.5, 0, 1, 0, car.name + "_ukf");
  viewer->addArrow(pcl::PointXYZ(ukf.x_[0], ukf.x_[1], 3.5),
                   pcl::PointXYZ(ukf.x_[0] + ukf.x_[2] * std::cos(ukf.x_[3]),
                                 ukf.x_[1] + ukf.x_[2] * std::sin(ukf.x_[3]), 3.5),
                   0, 1, 0, car.name + "_ukf_vel");
  if (time <= 0) return;
  const double dt = time / steps;
  for (double t = dt; t <= time; t += dt) {
    ukf.predict(dt);
    const std::string id = car.name + "_ukf" + std::to_string(t);
    viewer->addSphere(pcl::PointXYZ(ukf.x_[0], ukf.x_[1], 3.5), 0.5, 0, 1, 0, id);
    viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY,
                                        1.0 - 0.8 * (t / time), id);
  }
}

VectorXd Sensors::rmse(const std::vector<VectorXd>& estimations,
                       const std::vector<VectorXd>& ground_truth) {
  VectorXd rmse = VectorXd::Zero(4);
  if (estimations.empty() || estimations.size() != ground_truth.size()) {
    std::cerr << "rmse: estimations and ground truth must be non-empty and of equal size\n";
    return rmse;
  }
  for (size_t i = 0; i < estimations.size(); ++i) {
    const VectorXd residual = estimations[i] - ground_truth[i];
    rmse += residual.cwiseProduct(residual);
  }
  return (rmse / static_cast<double>(estimations.size())).cwiseSqrt();
}

void Sensors::savePcd(const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud, const std::string& file) {
  pcl::io::savePCDFileASCII(file, *cloud);
  std::cerr << "Saved " << cloud->points.size() << " points to " << file << std::endl;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr Sensors::loadPcd(const std::string& file) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
  if (pcl::io::loadPCDFile<pcl::PointXYZ>(file, *cloud) == -1)
    std::cerr << "Couldn't read " << file << std::endl;
  return cloud;
}
