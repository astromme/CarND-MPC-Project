#ifndef MPC_H
#define MPC_H

#include <vector>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "json.hpp"

// for convenience
using json = nlohmann::json;


using namespace std;
using CppAD::AD;

typedef CPPAD_TESTVECTOR(double) Dvector;

class MPC {
 public:
  MPC(json config);
  //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
  // the points in the simulator are connected by a Green line
  vector<double> mpc_x_vals;
  vector<double> mpc_y_vals;

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  json config;

  size_t N;
  double dt;
  double ref_v;

  // The solver takes all the state variables and actuator
  // variables in a singular vector. Thus, we should to establish
  // when one variable starts and another ends to make our lifes easier.
  size_t x_start;
  size_t y_start;
  size_t psi_start;
  size_t v_start;
  size_t cte_start;
  size_t epsi_start;
  size_t delta_start;
  size_t a_start;

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

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
