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
      rho_(base_particles_.rho_), Vol_(base_particles_.Vol_),
      color_gradient_(*base_particles_.getVariableByName<Vecd>("ColorGradient")),
      surface_tension_stress_(*base_particles_.getVariableByName<Matd>("SurfaceTensionStress")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
