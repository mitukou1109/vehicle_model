#pragma once

#include <Eigen/Dense>

template <std::size_t num_of_states, std::size_t num_of_inputs, std::size_t num_of_params>
class VehicleModel
{
public:
  using State = Eigen::Matrix<double, num_of_states, 1>;
  using Input = Eigen::Matrix<double, num_of_inputs, 1>;
  using Jacobian = Eigen::Matrix<double, num_of_states, num_of_states>;
  using Param = Eigen::Matrix<double, num_of_params, 1>;

  VehicleModel(const Param& param) : param_(param){};

  virtual ~VehicleModel();

  virtual State stateFunction(const State& x, const Input& u) = 0;

  virtual Jacobian stateJacobian(const State& x, const Input& u) = 0;

protected:
  Param param_;
};
