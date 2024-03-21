#pragma once

#include "vehicle_model/vehicle_model.hpp"

class SemiTrailerModel : public VehicleModel<4, 2, 3>
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

  SemiTrailerModel(const Param& param) : VehicleModel(param){};

  ~SemiTrailerModel() override{};

  State stateFunction(const State& x, const Input& u) override;

  Jacobian stateJacobian(const State& x, const Input& u) override;
};