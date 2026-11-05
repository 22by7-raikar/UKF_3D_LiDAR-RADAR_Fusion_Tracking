// Three-lane highway scene: an ego car and three traffic cars that
// accelerate and change lanes. Each traffic car is tracked by its own UKF
// from simulated lidar and radar, and the running RMSE is checked against
// the accuracy targets.
#ifndef HIGHWAY_H
#define HIGHWAY_H

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "../render/render.h"
#include "lidar.h"
#include "sensors.h"

class Highway {
 public:
  std::vector<Car> traffic;
  Car egoCar;
  Sensors sensors;
  std::unique_ptr<Lidar> lidar;

  // RMSE targets for [px, py, vx, vy], checked after the first second.
  std::vector<double> rmseThreshold = {0.30, 0.16, 0.95, 0.70};
  std::vector<double> rmseFailLog = {0.0, 0.0, 0.0, 0.0};
  bool pass = true;

  // Parameters
  // --------------------------------
  // Which cars to track with a UKF
  std::vector<bool> trackCars = {true, true, true};
  // Draw the sensor measurements
  bool visualize_lidar = true;
  bool visualize_radar = true;
  // Draw the recorded lidar point clouds instead of the car models
  bool visualize_pcd = false;
  // Draw the UKF's predicted path over this many seconds, in this many steps
  double projectedTime = 0;
  int projectedSteps = 0;
  // --------------------------------

  explicit Highway(pcl::visualization::PCLVisualizer::Ptr& viewer) {
    egoCar = Car(Vect3(0, 0, 0), Vect3(4, 2, 2), Color(0, 1, 0), 0, 0, 2, "egoCar");

    // Each actuation is (time [us], acceleration [m/s^2], steering [rad]).
    Car car1(Vect3(-10, 4, 0), Vect3(4, 2, 2), Color(0, 0, 1), 5, 0, 2, "car1");
    car1.setInstructions({Actuation(0.5e6, 0.5, 0.0), Actuation(2.2e6, 0.0, -0.2),
                          Actuation(3.3e6, 0.0, 0.2), Actuation(4.4e6, -2.0, 0.0)});

    Car car2(Vect3(25, -4, 0), Vect3(4, 2, 2), Color(0, 0, 1), -6, 0, 2, "car2");
    car2.setInstructions({Actuation(4.0e6, 3.0, 0.0), Actuation(8.0e6, 0.0, 0.0)});

    Car car3(Vect3(-12, 0, 0), Vect3(4, 2, 2), Color(0, 0, 1), 1, 0, 2, "car3");
    car3.setInstructions({Actuation(0.5e6, 2.0, 1.0), Actuation(1.0e6, 2.5, 0.0),
                          Actuation(3.2e6, 0.0, -1.0), Actuation(3.3e6, 2.0, 0.0),
                          Actuation(4.5e6, 0.0, 0.0), Actuation(5.5e6, -2.0, 0.0),
                          Actuation(7.5e6, 0.0, 0.0)});

    traffic = {car1, car2, car3};
    for (size_t i = 0; i < traffic.size(); ++i)
      if (trackCars[i]) traffic[i].setUKF(UKF());

    lidar = std::make_unique<Lidar>(traffic, 0);

    renderHighway(0, viewer);
    egoCar.render(viewer);
    for (Car& car : traffic) car.render(viewer);
  }

  void stepHighway(double egoVelocity, long long timestamp, int frame_per_sec,
                   pcl::visualization::PCLVisualizer::Ptr& viewer) {
    if (visualize_pcd) {
      auto cloud = Sensors::loadPcd(std::string(PCD_DIR) + "/highway_" +
                                    std::to_string(timestamp) + ".pcd");
      renderPointCloud(viewer, cloud, "trafficCloud",
                       Color(184.0f / 256, 223.0f / 256, 252.0f / 256));
    }

    renderHighway(egoVelocity * timestamp / 1e6, viewer);
    egoCar.render(viewer);

    for (size_t i = 0; i < traffic.size(); ++i) {
      Car& car = traffic[i];
      car.move(1.0 / frame_per_sec, timestamp);
      if (!visualize_pcd) car.render(viewer);
      if (!trackCars[i]) continue;

      Eigen::VectorXd truth(4);
      truth << car.position.x, car.position.y, car.velocity * std::cos(car.angle),
          car.velocity * std::sin(car.angle);
      sensors.ground_truth.push_back(truth);
      sensors.lidarSense(car, viewer, timestamp, visualize_lidar);
      sensors.radarSense(car, egoCar, viewer, timestamp, visualize_radar);
      sensors.renderTrack(car, viewer, projectedTime, projectedSteps);

      const double v = car.ukf.x_(2);
      const double yaw = car.ukf.x_(3);
      Eigen::VectorXd estimate(4);
      estimate << car.ukf.x_(0), car.ukf.x_(1), std::cos(yaw) * v, std::sin(yaw) * v;
      sensors.estimations.push_back(estimate);
    }

    const Eigen::VectorXd rmse = Sensors::rmse(sensors.estimations, sensors.ground_truth);
    const char* labels[] = {" X: ", " Y: ", "Vx: ", "Vy: "};
    const char* ids[] = {"x", "y", "vx", "vy"};
    viewer->addText("Accuracy - RMSE:", 30, 300, 20, 1, 1, 1, "rmse");
    for (int k = 0; k < 4; ++k)
      viewer->addText(labels[k] + std::to_string(rmse[k]), 30, 275 - 25 * k, 20, 1, 1, 1,
                      std::string("rmse_") + ids[k]);

    if (timestamp > 1.0e6) {
      for (int k = 0; k < 4; ++k) {
        if (rmse[k] > rmseThreshold[k]) {
          rmseFailLog[k] = rmse[k];
          pass = false;
        }
      }
    }
    if (!pass) {
      viewer->addText("RMSE Failed Threshold", 30, 150, 20, 1, 0, 0, "rmse_fail");
      for (int k = 0; k < 4; ++k)
        if (rmseFailLog[k] > 0)
          viewer->addText(labels[k] + std::to_string(rmseFailLog[k]), 30, 125 - 25 * k, 20, 1,
                          0, 0, std::string("rmse_fail_") + ids[k]);
    }
  }
};

#endif  // HIGHWAY_H
