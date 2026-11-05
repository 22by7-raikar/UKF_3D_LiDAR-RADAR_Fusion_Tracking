# Design

The filter tracks each car with the CTRV (constant turn rate and velocity)
motion model and the state

$$x = [p_x, p_y, v, \psi, \dot{\psi}]^T$$

where $p_x, p_y$ are position, $v$ the speed, $\psi$ the yaw and
$\dot{\psi}$ the yaw rate.

The scene in [src/sim/highway.h](../src/sim/highway.h) has three traffic cars,
each measured by a simulated lidar and radar every frame. RMSE is computed
against ground truth in the order $[p_x, p_y, v_x, v_y]$ and must stay below
$[0.30, 0.16, 0.95, 0.70]$ after the first second.

## UKF Implementation

The filter is in [src/ukf.h](../src/ukf.h) and [src/ukf.cpp](../src/ukf.cpp).

### 1. Initialization

The constructor initializes the filter dimensions, sigma-point weights, noise parameters, and covariance storage.

Important settings:

- State dimension: $n_x = 5$
- Augmented dimension: $n_{aug} = 7$
- Sigma point parameter: $\lambda = 3 - n_{aug}$
- Process noise tuning:
  - $\sigma_a = 1.5$
  - $\sigma_{\ddot{\psi}} = 0.5$

In `processMeasurement`, the first measurement initializes the filter state:

- lidar initializes $p_x$ and $p_y$ directly
- radar initializes $p_x$ and $p_y$ by converting from polar to Cartesian coordinates, and seeds $v$ from $|\dot{\rho}|$

The initial covariance matrix is also tightened for position using the corresponding sensor noise values.

### 2. Prediction Step

The prediction is in `predict`.

The algorithm performs these steps:

1. Create an augmented state vector and covariance matrix.
2. Generate augmented sigma points using the Cholesky decomposition.
3. Propagate each sigma point through the CTRV motion model.
4. Add process noise for longitudinal acceleration and yaw acceleration.
5. Recover the predicted mean state and covariance using the sigma-point weights.

Angle normalization (`NormalizeAngle`) keeps wrapped headings from producing large residuals.

Prediction also runs for $\Delta t = 0$: the simulation can deliver a lidar update and a radar update at the same timestamp, so generating sigma points for the second update is still necessary even when no time has elapsed.

### 3. Lidar Update

The lidar update is in `updateLidar`.

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

The radar update is in `updateRadar`.

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

With these changes every RMSE stays within its target.

## Highway Parameters

Visualization and testing parameters are in [src/sim/highway.h](../src/sim/highway.h):

- `trackCars` toggles which traffic vehicles are tracked
- `visualize_lidar` and `visualize_radar` toggle sensor marker rendering
- `visualize_pcd` shows traffic from the lidar point-cloud perspective
- `projectedTime` and `projectedSteps` show future UKF predictions

These settings are useful for debugging individual targets and visualizing how the CTRV model extrapolates the path of each vehicle.

## Efficiency

- Sigma point weights are computed once, in the constructor.
- The predicted sigma points are stored and reused by the radar update.

## Results

RMSE over the 10 s run:

| | px [m] | py [m] | vx [m/s] | vy [m/s] |
|---|---|---|---|---|
| Target | 0.30 | 0.16 | 0.95 | 0.70 |
| UKF | 0.068 | 0.080 | 0.433 | 0.477 |
