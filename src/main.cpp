// Runs the highway scene for 10 s at 30 frames per second, tracking each
// traffic car with a UKF and showing the result in a PCL viewer.

#include <iostream>

#include "sim/highway.h"

int main() {
  pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
  viewer->setBackgroundColor(0, 0, 0);
  viewer->initCameraParameters();
  viewer->setCameraPosition(-26, 0, 15.0, 25, 0, 0, 0, 0, 1);

  Highway highway(viewer);

  const int frame_per_sec = 30;
  const int sec_interval = 10;
  const double egoVelocity = 25;

  for (int frame = 0; frame < frame_per_sec * sec_interval; ++frame) {
    const long long time_us = 1000000LL * frame / frame_per_sec;
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();
    highway.stepHighway(egoVelocity, time_us, frame_per_sec, viewer);
    viewer->spinOnce(1000 / frame_per_sec);
  }
  const Eigen::VectorXd rmse = Sensors::rmse(highway.sensors.estimations,
                                             highway.sensors.ground_truth);
  std::cout << "RMSE [px py vx vy]: " << rmse.transpose() << "\n"
            << (highway.pass ? "within" : "outside") << " the RMSE targets after 1 s\n";
  return highway.pass ? 0 : 1;
}
