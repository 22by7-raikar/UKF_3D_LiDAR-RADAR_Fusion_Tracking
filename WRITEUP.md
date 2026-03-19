# SFND Unscented Kalman Filter Project Writeup

## Project Overview

This project implements an Unscented Kalman Filter (UKF) to track multiple vehicles on a highway using fused lidar and radar measurements. The tracker uses a Constant Turn Rate and Velocity (CTRV) motion model and estimates the state vector

$$x = [p_x, p_y, v, \psi, \dot{\psi}]^T$$

where $p_x$ and $p_y$ are position, $v$ is the velocity magnitude, $\psi$ is yaw, and $\dot{\psi}$ is yaw rate.

The project environment in [src/highway.h](src/highway.h) creates three traffic vehicles and feeds each one noisy lidar and radar measurements. RMSE is computed against ground truth in the order $[p_x, p_y, v_x, v_y]$ and must remain below $[0.30, 0.16, 0.95, 0.70]$ after one second.

## Build and Run

From the project root:

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`
5. `./ukf_highway`

The build configuration was updated in [CMakeLists.txt](CMakeLists.txt) to use a modern C++ standard required by the installed PCL version.

## UKF Implementation

The core implementation is in [src/ukf.cpp](src/ukf.cpp) and [src/ukf.h](src/ukf.h).

### 1. Initialization

The constructor initializes the filter dimensions, sigma-point weights, noise parameters, and covariance storage.

Important settings:

- State dimension: $n_x = 5$
- Augmented dimension: $n_{aug} = 7$
- Sigma point parameter: $\lambda = 3 - n_{aug}$
- Process noise tuning:
  - $\sigma_a = 1.5$
  - $\sigma_{\ddot{\psi}} = 0.5$

In [ProcessMeasurement](src/ukf.cpp#L96-L154), the first measurement initializes the filter state:

- lidar initializes $p_x$ and $p_y$ directly
- radar initializes $p_x$ and $p_y$ by converting from polar to Cartesian coordinates, and seeds $v$ from $|\dot{\rho}|$

The initial covariance matrix is also tightened for position using the corresponding sensor noise values.

### 2. Prediction Step

The prediction logic is implemented in [Prediction](src/ukf.cpp#L157-L250).

The algorithm performs these steps:

1. Create an augmented state vector and covariance matrix.
2. Generate augmented sigma points using the Cholesky decomposition.
3. Propagate each sigma point through the CTRV motion model.
4. Add process noise for longitudinal acceleration and yaw acceleration.
5. Recover the predicted mean state and covariance using the sigma-point weights.

Angle normalization is handled by `NormalizeAngle()` so wrapped headings do not produce large residuals.

A subtle but important fix was made in [ProcessMeasurement](src/ukf.cpp#L140-L146): prediction is also run for $\Delta t = 0$. This project can deliver a lidar update and a radar update at the same timestamp, so generating sigma points for the second update is still necessary even when no time has elapsed.

### 3. Lidar Update

The lidar update is implemented in [UpdateLidar](src/ukf.cpp#L252-L276).

Because lidar directly measures Cartesian position, the measurement model is linear:

$$z = Hx + v$$

with

$$
H = \begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0
\end{bmatrix}
$$

A standard Kalman update is used:

- compute residual $y = z - Hx$
- compute innovation covariance $S$
- compute Kalman gain $K$
- update state and covariance

### 4. Radar Update

The radar update is implemented in [UpdateRadar](src/ukf.cpp#L278-L345).

Radar is nonlinear, so predicted sigma points are transformed into measurement space:

$$z = [\rho, \phi, \dot{\rho}]^T$$

For each sigma point:

- $\rho = \sqrt{p_x^2 + p_y^2}$
- $\phi = \text{atan2}(p_y, p_x)$
- $\dot{\rho} = \frac{p_x v_x + p_y v_y}{\rho}$

The update then computes:

- predicted radar mean $\hat{z}$
- innovation covariance $S$
- cross-correlation matrix $T_c$
- Kalman gain $K = T_c S^{-1}$

Finally, the state and covariance are corrected using the radar residual. Angle normalization is applied both to sigma-point residuals and the final measurement residual.

## Tuning and Debugging

The main RMSE issue during development was a slight violation of the $v_y$ threshold. The final passing version came from two changes:

1. Running prediction even when $\Delta t = 0$ so same-timestamp sensor updates use valid predicted sigma points.
2. Tuning initialization and process noise:
   - reduced process noise to more realistic values
   - improved first-measurement covariance setup
   - seeded radar-based initial speed from measured radial velocity magnitude

These updates were enough to bring the RMSE within rubric limits.

## Highway Parameters

Useful visualization and testing parameters are defined in [src/highway.h](src/highway.h#L20-L29):

- `trackCars` toggles which traffic vehicles are tracked
- `visualize_lidar` and `visualize_radar` toggle sensor marker rendering
- `visualize_pcd` shows traffic from the lidar point-cloud perspective
- `projectedTime` and `projectedSteps` show future UKF predictions

These settings are useful for debugging individual targets and visualizing how the CTRV model extrapolates the path of each vehicle.

## Rubric Discussion

### Compiling and Testing

The project compiles successfully with `cmake` and `make`.

### Code Efficiency

The implementation avoids repeated recalculation where possible:

- sigma-point weights are precomputed once in the constructor
- predicted sigma points are stored and reused during radar updates
- matrix operations are kept local to each update step

The code prioritizes clarity and numerical stability while still following efficient UKF practices.

### Accuracy

The final implementation achieves RMSE below the rubric thresholds after the first second of simulation:

$$[p_x, p_y, v_x, v_y] \le [0.30, 0.16, 0.95, 0.70]$$

### Correct Algorithm

The filter follows the expected UKF workflow taught in the course:

1. initialize from the first measurement
2. predict forward using the CTRV model
3. apply a lidar or radar measurement update depending on sensor type
4. repeat for each incoming measurement

## Conclusion

The final project implements a working UKF for multi-vehicle highway tracking with lidar/radar fusion. It compiles successfully, follows the expected algorithmic flow, and meets the accuracy requirements in the project rubric.
