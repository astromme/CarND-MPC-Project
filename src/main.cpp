#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  // Read config file
  std::ifstream in("config.json");
  json config;
  in >> config;
  assert(config.count("N") > 0);
  assert(config.count("dt") > 0);
  assert(config.count("ref_v") > 0);
  assert(config.count("cost") > 0);
  string costs[] = {"cte", "epsi", "v", "delta", "a", "delta_gap", "a_gap"};
  for (auto cost : costs) {
    assert(config["cost"].count(cost) > 0);
  }

  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc(config);

  h.onMessage([&mpc, config](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    if (config.count("debug_incoming_message") > 0 && config["debug_incoming_message"]) {
      cout << sdata << endl;
    }
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double a = j[1]["throttle"];
          double delta = j[1]["steering_angle"];

          Eigen::VectorXd ptsxEigen(ptsx.size());
          Eigen::VectorXd ptsyEigen(ptsy.size());

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          for (auto i=0; i<ptsx.size() ; ++i){
            double x = (ptsx[i] - px) * cos(psi) + (ptsy[i] - py) * sin(psi);
            double y = -(ptsx[i] - px) * sin(psi) + (ptsy[i] - py) * cos(psi);

            next_x_vals.push_back(x);
            next_y_vals.push_back(y);

            ptsxEigen(i) = x;
            ptsyEigen(i) = y;
          }

          auto coeffs = polyfit(ptsxEigen, ptsyEigen, 3);
          if (config.count("debug_coeffs") > 0 && config["debug_coeffs"]) {
            for (auto i=0; i<coeffs.size() ; ++i) { cout << coeffs[i] << ", "; }
            cout << endl;
          }

          // Apply latency
          double psides0 = atan(coeffs[1]);
          double latency_dt = 1/1000.0 * ((double) config["latency_ms"]);
          double x_latency = 0 + v * cos(0) * latency_dt;
          double y_latency = 0 + v * sin(0) * latency_dt;
          double psi_latency = 0 - v / mpc.Lf * delta * latency_dt;
          double v_latency = v + a * latency_dt;
          double epsi_latency = 0 - psides0 + v * delta / mpc.Lf * latency_dt;
          double cte_latency = polyeval(coeffs, x_latency) - y_latency;

          Eigen::VectorXd state(6);
          state << x_latency, y_latency, psi_latency, v_latency, cte_latency, epsi_latency; //px, py, psi, v, cte, epsi;

          auto vars = mpc.Solve(state, coeffs);

          if (config.count("debug_solution") > 0 && config["debug_solution"]) {
            int i = 0;
            for (auto type : {"x", "y", "psi", "v", "cte", "epsi", "delta", "a"}) {
              std::cout << std::setfill(' ') << setw(7) << type << ": ";
              for (auto j=0; j < mpc.N; j++) {
                if (mpc.solution.x[i] > 100) {
                  std::cout << std::fixed << std::setfill(' ') << setw(5) << mpc.solution.x[i];
                } else {
                  std::cout << std::fixed << std::setfill(' ') << setw(5) << std::setprecision(2) << mpc.solution.x[i];
                }
                i += 1;
                if (j < (mpc.N-1)) {
                  std::cout << " >";
                }
              }
              std::cout << std::endl;
            }
            std::cout << std::endl;
          }

          double steer_value = vars[6];
          double throttle_value = vars[7];
          if (config["disable_throttle"]) {
            throttle_value = 0;
          }

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value / deg2rad(25);
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory
          msgJson["mpc_x"] = mpc.mpc_x_vals;
          msgJson["mpc_y"] = mpc.mpc_y_vals;

          //Display the waypoints/reference line
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          if (config.count("debug_response") > 0 && config["debug_response"]) {
            std::cout << msg << std::endl;
          }
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(config["latency_ms"]));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
