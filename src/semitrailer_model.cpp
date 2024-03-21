#include "vehicle_model/semitrailer_model.hpp"

SemiTrailerModel::State SemiTrailerModel::stateFunction(const State& x, const Input& u)
{
  State x_dot;
  x_dot(X) = u(V) * std::cos(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) *
             std::cos(x(THETA));
  x_dot(Y) = u(V) * std::cos(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) *
             std::sin(x(THETA));
  x_dot(THETA) = u(V) * (std::sin(x(BETA)) / param_(L_2) -
                         param_(M) / (param_(L_1) * param_(L_2)) * std::cos(x(BETA)) * std::tan(u(ALPHA)));
  x_dot(BETA) = u(V) * (std::tan(u(ALPHA)) / param_(L_1) - std::sin(x(BETA)) / param_(L_2) +
                        param_(M) / (param_(L_1) * param_(L_2)) * std::cos(x(BETA)) * std::tan(u(ALPHA)));
  return x_dot;
}

SemiTrailerModel::Jacobian SemiTrailerModel::stateJacobian(const State& x, const Input& u)
{
  Jacobian A;
  A(X, X) = 0;
  A(X, Y) = 0;
  A(X, THETA) = -u(V) * std::cos(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) *
                std::sin(x(THETA));
  A(X, BETA) = u(V) *
               (-std::sin(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) +
                param_(M) / param_(L_1) / std::cos(x(BETA)) * std::tan(u(ALPHA))) *
               std::cos(x(THETA));
  A(Y, X) = 0;
  A(Y, Y) = 0;
  A(Y, THETA) = u(V) * std::cos(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) *
                std::cos(x(THETA));
  A(Y, BETA) = u(V) *
               (-std::sin(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) +
                param_(M) / param_(L_1) / std::cos(x(BETA)) * std::tan(u(ALPHA))) *
               std::sin(x(THETA));
  A(THETA, X) = 0;
  A(THETA, Y) = 0;
  A(THETA, THETA) = 0;
  A(THETA, BETA) = u(V) * (std::cos(x(BETA)) / param_(L_2) +
                           param_(M) / (param_(L_1) * param_(L_2)) * std::sin(x(BETA)) * std::tan(u(ALPHA)));
  A(BETA, X) = 0;
  A(BETA, Y) = 0;
  A(BETA, THETA) = 0;
  A(BETA, BETA) = u(V) * (-std::cos(x(BETA)) / param_(L_2) -
                          param_(M) / (param_(L_1) * param_(L_2)) * std::sin(x(BETA)) * std::tan(u(ALPHA)));
  return A;
}
