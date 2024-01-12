#pragma once

#include "surface_tension.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
SurfaceTensionForce<DataDelegationType>::SurfaceTensionForce(BaseRelationType &base_relation)
    : ForcePrior(base_relation.getSPHBody(), "SurfaceTensionForce"), DataDelegationType(base_relation),
      rho_(this->particles_->rho_), Vol_(this->particles_->Vol_), force_prior_(this->particles_->force_prior_),
      color_gradient_(*this->particles_->getVariableByName<Vecd>("ColorGradient")),
      surface_tension_stress_(*this->particles_->getVariableByName<Matd>("SurfaceTensionStress")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
