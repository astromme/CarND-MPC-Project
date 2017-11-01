#include "MPC.h"
#include "Eigen-3.3/Eigen/Core"

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;


class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  MPC *mpc; // weak pointer
  json config;
  FG_eval(Eigen::VectorXd coeffs, MPC *mpc, json config) {
    this->coeffs = coeffs;
    this->mpc = mpc;
    this->config = config;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // The part of the cost based on the reference state.
    for (int t = 0; t < mpc->N; t++) {
      fg[0] += config["cost"]["cte"]*CppAD::pow(vars[mpc->cte_start + t], 2);
      fg[0] += config["cost"]["epsi"]*CppAD::pow(vars[mpc->epsi_start + t], 2);
      fg[0] += config["cost"]["v"]*CppAD::pow(vars[mpc->v_start + t] - mpc->ref_v, 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < mpc->N; t++) {
      fg[0] += config["cost"]["delta"]*CppAD::pow(vars[mpc->delta_start + t], 2);
      fg[0] += config["cost"]["a"]*CppAD::pow(vars[mpc->a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < mpc->N - 1; t++) {
      fg[0] += config["cost"]["delta_gap"]*CppAD::pow(vars[mpc->delta_start + t + 1] - vars[mpc->delta_start + t], 2);
      fg[0] += config["cost"]["a_gap"]*CppAD::pow(vars[mpc->a_start + t + 1] - vars[mpc->a_start + t], 2);
    }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + mpc->x_start] = vars[mpc->x_start];
    fg[1 + mpc->y_start] = vars[mpc->y_start];
    fg[1 + mpc->psi_start] = vars[mpc->psi_start];
    fg[1 + mpc->v_start] = vars[mpc->v_start];
    fg[1 + mpc->cte_start] = vars[mpc->cte_start];
    fg[1 + mpc->epsi_start] = vars[mpc->epsi_start];

    // The rest of the constraints
    for (int t = 1; t < mpc->N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[mpc->x_start + t];
      AD<double> y1 = vars[mpc->y_start + t];
      AD<double> psi1 = vars[mpc->psi_start + t];
      AD<double> v1 = vars[mpc->v_start + t];
      AD<double> cte1 = vars[mpc->cte_start + t];
      AD<double> epsi1 = vars[mpc->epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[mpc->x_start + t - 1];
      AD<double> y0 = vars[mpc->y_start + t - 1];
      AD<double> psi0 = vars[mpc->psi_start + t - 1];
      AD<double> v0 = vars[mpc->v_start + t - 1];
      AD<double> cte0 = vars[mpc->cte_start + t - 1];
      AD<double> epsi0 = vars[mpc->epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[mpc->delta_start + t - 1];
      AD<double> a0 = vars[mpc->a_start + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
      AD<double> psides0 = CppAD::atan(coeffs[1] + (2 * coeffs[2] * x0) + (3 * coeffs[3]* CppAD::pow(x0, 2)));

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      fg[1 + mpc->x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * mpc->dt);
      fg[1 + mpc->y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * mpc->dt);
      fg[1 + mpc->psi_start + t] = psi1 - (psi0 - v0 * delta0 / Lf * mpc->dt);
      fg[1 + mpc->v_start + t] = v1 - (v0 + a0 * mpc->dt);
      fg[1 + mpc->cte_start + t] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * mpc->dt));
      fg[1 + mpc->epsi_start + t] =
          epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * mpc->dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC(json config) {
  // We set the number of timesteps to 25
  // and the timestep evaluation frequency or evaluation
  // period to 0.05.
  this->N = config["N"];
  this->dt = config["dt"];
  this->ref_v = config["ref_v"];

  this->config = config;

  // The solver takes all the state variables and actuator
  // variables in a singular vector. Thus, we should to establish
  // when one variable starts and another ends to make our lifes easier.
  x_start = 0;
  y_start = x_start + N;
  psi_start = y_start + N;
  v_start = psi_start + N;
  cte_start = v_start + N;
  epsi_start = cte_start + N;
  delta_start = epsi_start + N;
  a_start = delta_start + N - 1;
}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  size_t num_states = 6;
  size_t num_inputs = 2;
  size_t n_vars = num_states * N + num_inputs * N;

  // TODO: Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }


  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }


  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }


  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);

  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }


  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;


  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, this, config);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);


  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  mpc_x_vals.clear();
  mpc_y_vals.clear();

  for (int i=0; i<N; i++) {
    double carx = solution.x[x_start+i];
    double cary = solution.x[y_start+i];

    mpc_x_vals.push_back(carx);
    mpc_y_vals.push_back(cary);
  }

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  return {solution.x[x_start + 1],   solution.x[y_start + 1],
          solution.x[psi_start + 1], solution.x[v_start + 1],
          solution.x[cte_start + 1], solution.x[epsi_start + 1],
          solution.x[delta_start],   solution.x[a_start]};
}
