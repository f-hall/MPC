#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Tested different stepsizes and number of steps (N = 5..30) (d = 0.01..0.5)
// Get the best results with this values
size_t N = 20;
double dt = 0.05;

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

// ref values for cross track, epsilon and maximum speed.
double ref_cte = 0;
double ref_epsi = 0;
// You can chose even lower or higher speed (tested until 60 where it begins to slinger a bit on the street)
double ref_v = 45;

// Setting switch values where the vectors starts with theses values (example: v-values starts at 3*N = 60)
size_t switch_x = 0;
size_t switch_y = N;
size_t switch_psi = 2*N;
size_t switch_v = 3*N;
size_t switch_cte = 4*N;
size_t switch_epsi = 5*N;
size_t switch_delta = 6*N;
size_t switch_a = 7*N-1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {

    // Setting the vector for constraints
    fg[0] = 0;
    // The coeffs in the constraints (like 150, 8 and 250)
    // are used to ensure that we don't have to hectic movements
    // You can look it up on one of the last project-pages at udacity for more infos on this!
    for(int i = 0; i<N; i++)
        {
			fg[0] += CppAD::pow(vars[switch_cte + i] - ref_cte,2)
                  + 150*CppAD::pow(vars[switch_epsi + i] - ref_epsi,2)
			      + CppAD::pow(vars[switch_v + i] - ref_v,2);
            if(i<N-1)
            {
                fg[0] += 8*CppAD::pow(vars[switch_delta + i],2)
                      + 8*CppAD::pow(vars[switch_a + i],2);
            }
            if(i<N-2)
            {
                fg[0] += 250*CppAD::pow(vars[switch_delta + i] - vars[switch_delta + i + 1],2)
                      + 15*CppAD::pow(vars[switch_a + i] - vars[switch_a + i +1],2);
            }
		}
    fg[1 + switch_x] = vars[switch_x];
    fg[1 + switch_y] = vars[switch_y];
    fg[1 + switch_psi] = vars[switch_psi];
    fg[1 + switch_v] = vars[switch_v];
    fg[1 + switch_cte] = vars[switch_cte];
    fg[1 + switch_epsi] = vars[switch_epsi];

    AD<double> delta_help = vars[switch_delta];
    AD<double> a_help = vars[switch_a];

    // Going through all steps for the other contraints
    for (int i = 0; i < N - 1; i++)
        {

        AD<double> x_0 = vars[switch_x + i];
        AD<double> x_1 = vars[switch_x + i + 1];

        AD<double> y_0 = vars[switch_y + i];
        AD<double> y_1 = vars[switch_y + i + 1];

        AD<double> psi_0 = vars[switch_psi + i];
        AD<double> psi_1 = vars[switch_psi + i + 1];

        AD<double> v_0 = vars[switch_v + i];
        AD<double> v_1 = vars[switch_v + i + 1];

        AD<double> cte_0 = vars[switch_cte + i];
        AD<double> cte_1 = vars[switch_cte + i + 1];

        AD<double> epsi_0 = vars[switch_epsi + i];
        AD<double> epsi_1 = vars[switch_epsi + i + 1];

        AD<double> delta;
        AD<double> a;

        // I use this to avoid problems with the 100ms latency given in this project
        // So I don't use the i values for delta and a but the i+2 to get the timing right (seems to work really good)
        // The idea for this comes from a forum entry at the udacity-forums for this project - credits to them!
        int time_help = 0.01/dt;

        if (i >= time_help)
        {
            delta = vars[switch_delta + i - time_help];
            a = vars[switch_a + i - time_help];
        }
        else
        {
            delta = delta_help;
            a = a_help;
        }

        fg[2 + switch_x + i] = x_1 - (x_0 + v_0 * CppAD::cos(psi_0) * dt);
        fg[2 + switch_y + i] = y_1 - (y_0 + v_0 * CppAD::sin(psi_0) * dt);
        fg[2 + switch_psi + i] = psi_1 - (psi_0 + v_0/Lf * delta * dt);
        fg[2 + switch_v + i] = v_1 - (v_0 + a * dt);

        AD<double> f_0 = coeffs[0] + coeffs[1] * x_0 + coeffs[2] * CppAD::pow(x_0,2)+ coeffs[3] * CppAD::pow(x_0,3);
        AD<double> psi_0_help = CppAD::atan(coeffs[1]);

        fg[2 + switch_cte + i] = cte_1 - ((f_0 - y_0) + (v_0 * CppAD::sin(epsi_0) * dt));
        fg[2 + switch_epsi + i] = epsi_1 - ((psi_0 - psi_0_help) + v_0/Lf * delta * dt);

        }
    }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x_state = state[0];
  double y_state = state[1];
  double psi_state = state[2];
  double v_state = state[3];
  double cte_state = state[4];
  double epsi_state = state[5];

  double s_state = state[6];
  double t_state = state[7];

  // Setting number of model variables - the state vector has 6 elements and the
  // actuator vector 2 (Combined in state[0]-state[7]
  size_t n_vars = N*6 + (N-1)*2;

  // Setting number of constraints
  size_t n_constraints = N*6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.

  // Setting just the switch values to the state
  // values
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  vars[switch_x] = x_state;
  vars[switch_y] = y_state;
  vars[switch_psi] = psi_state;
  vars[switch_v] = v_state;

  vars[switch_cte] = cte_state;
  vars[switch_epsi] = epsi_state;

  vars[switch_delta] = s_state;
  vars[switch_a] = t_state;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Setting upper and lower bounds

  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  for (int i = 0; i < n_vars; i++)
  {
       if (0 <= i < switch_delta)
       {
            vars_lowerbound[i] = -42.0e15;
            vars_upperbound[i] = 42.0e15;
       }
       else if (i < switch_a)
       {
           vars_lowerbound[i] = -25.*M_PI/180.;
           vars_upperbound[i] = 25.*M_PI/180.;
       }
       else
       {
           vars_lowerbound[i] = -1.0;
           vars_upperbound[i] = 1.0;
       }
  }

  constraints_lowerbound[switch_x] = x_state;
  constraints_upperbound[switch_x] = x_state;

  constraints_lowerbound[switch_y] = y_state;
  constraints_upperbound[switch_y] = y_state;

  constraints_lowerbound[switch_psi] = psi_state;
  constraints_upperbound[switch_psi] = psi_state;

  constraints_lowerbound[switch_v] = v_state;
  constraints_upperbound[switch_v] = v_state;

  constraints_lowerbound[switch_cte] = cte_state;
  constraints_upperbound[switch_cte] = cte_state;

  constraints_lowerbound[switch_epsi] = epsi_state;
  constraints_upperbound[switch_epsi] = epsi_state;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

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

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  return {-1*solution.x[switch_delta], solution.x[switch_a],
		solution.x[switch_x], solution.x[switch_x+1],solution.x[switch_x+2],
		solution.x[switch_x+3],solution.x[switch_x+4], solution.x[switch_x+5],
		solution.x[switch_y],solution.x[switch_y+1],solution.x[switch_y+2],
		solution.x[switch_y+3],solution.x[switch_y+4], solution.x[switch_y+5]};
}
