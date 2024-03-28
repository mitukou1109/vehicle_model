#pragma once

#include "vehicle_model/vehicle_model.hpp"

namespace vehicle_model
{
class SemitrailerModel : public VehicleModel<4, 2, 3>
{
public:
  enum : Eigen::Index
  {
    X = 0,
    Y = 1,
    THETA = 2,
    BETA = 3,
    V = 0,
    ALPHA = 1,
    L_1 = 0,
    M = 1,
    L_2 = 2
  };

  SemitrailerModel(const Param& param) : VehicleModel(param){};

  SemitrailerModel(const double tractor_length, const double hitch_length, const double trailer_length)
    : VehicleModel((Param() << tractor_length, hitch_length, trailer_length).finished()){};

  ~SemitrailerModel(){};

  State stateFunction(const State& x, const Input& u) override;

  StateJacobian stateJacobian(const State& x, const Input& u) override;

  InputJacobian inputJacobian(const State& x, const Input& u) override;
};
}  // namespace vehicle_model