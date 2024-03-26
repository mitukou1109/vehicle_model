#include "vehicle_model/semitrailer_model.hpp"

namespace vehicle_model
{
SemitrailerModel::State SemitrailerModel::stateFunction(const State& x, const Input& u)
{
  State f;
  f(X) = u(V) * std::cos(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) *
         std::cos(x(THETA));
  f(Y) = u(V) * std::cos(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) *
         std::sin(x(THETA));
  f(THETA) = u(V) * (std::sin(x(BETA)) / param_(L_2) -
                     param_(M) / (param_(L_1) * param_(L_2)) * std::cos(x(BETA)) * std::tan(u(ALPHA)));
  f(BETA) = u(V) * (std::tan(u(ALPHA)) / param_(L_1) - std::sin(x(BETA)) / param_(L_2) +
                    param_(M) / (param_(L_1) * param_(L_2)) * std::cos(x(BETA)) * std::tan(u(ALPHA)));
  return f;
}

SemitrailerModel::StateJacobian SemitrailerModel::stateJacobian(const State& x, const Input& u)
{
  StateJacobian dfdx;
  dfdx(X, X) = 0;
  dfdx(X, Y) = 0;
  dfdx(X, THETA) = -u(V) * std::cos(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) *
                   std::sin(x(THETA));
  dfdx(X, BETA) = u(V) *
                  (-std::sin(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) +
                   param_(M) / param_(L_1) / std::cos(x(BETA)) * std::tan(u(ALPHA))) *
                  std::cos(x(THETA));
  dfdx(Y, X) = 0;
  dfdx(Y, Y) = 0;
  dfdx(Y, THETA) = u(V) * std::cos(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) *
                   std::cos(x(THETA));
  dfdx(Y, BETA) = u(V) *
                  (-std::sin(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) +
                   param_(M) / param_(L_1) / std::cos(x(BETA)) * std::tan(u(ALPHA))) *
                  std::sin(x(THETA));
  dfdx(THETA, X) = 0;
  dfdx(THETA, Y) = 0;
  dfdx(THETA, THETA) = 0;
  dfdx(THETA, BETA) = u(V) * (std::cos(x(BETA)) / param_(L_2) +
                              param_(M) / (param_(L_1) * param_(L_2)) * std::sin(x(BETA)) * std::tan(u(ALPHA)));
  dfdx(BETA, X) = 0;
  dfdx(BETA, Y) = 0;
  dfdx(BETA, THETA) = 0;
  dfdx(BETA, BETA) = u(V) * (-std::cos(x(BETA)) / param_(L_2) -
                             param_(M) / (param_(L_1) * param_(L_2)) * std::sin(x(BETA)) * std::tan(u(ALPHA)));
  return dfdx;
}

SemitrailerModel::InputJacobian SemitrailerModel::inputJacobian(const State& x, const Input& u)
{
  InputJacobian dfdu;
  dfdu(X, V) =
      std::cos(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) * std::cos(x(THETA));
  dfdu(X, ALPHA) = u(V) * std::cos(x(BETA)) * param_(M) / param_(L_1) * std::tan(x(BETA)) /
                   std::pow(std::cos(u(ALPHA)), 2) * std::cos(x(THETA));
  dfdu(Y, V) =
      std::cos(x(BETA)) * (1 + param_(M) / param_(L_1) * std::tan(x(BETA)) * std::tan(u(ALPHA))) * std::sin(x(THETA));
  dfdu(Y, ALPHA) = u(V) * std::cos(x(BETA)) * param_(M) / param_(L_1) * std::tan(x(BETA)) /
                   std::pow(std::cos(u(ALPHA)), 2) * std::sin(x(THETA));
  dfdu(THETA, V) = std::sin(x(BETA)) / param_(L_2) -
                   param_(M) / (param_(L_1) * param_(L_2)) * std::cos(x(BETA)) * std::tan(u(ALPHA));
  dfdu(THETA, ALPHA) =
      -u(V) * param_(M) / (param_(L_1) * param_(L_2)) * std::cos(x(BETA)) / std::pow(std::cos(u(ALPHA)), 2);
  dfdu(BETA, V) = std::tan(u(ALPHA)) / param_(L_1) - std::sin(x(BETA)) / param_(L_2) +
                  param_(M) / (param_(L_1) * param_(L_2)) * std::cos(x(BETA)) * std::tan(u(ALPHA));
  dfdu(BETA, ALPHA) = u(V) * (1 / param_(L_1) + param_(M) / (param_(L_1) * param_(L_2)) * std::cos(x(BETA))) /
                      std::pow(std::cos(u(ALPHA)), 2);
  return dfdu;
}
}  // namespace vehicle_model
