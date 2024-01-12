#pragma once

#include "viscous_dynamics.h"

namespace SPH
{
  namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
ViscousForce<DataDelegationType>::ViscousForce(BaseRelationType &base_relation)
    : ForcePrior(base_relation.getSPHBody(), "ViscousForce"), DataDelegationType(base_relation),
      rho_(this->particles_->rho_), Vol_(this->particles_->Vol_), vel_(this->particles_->vel_),
      mu_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceViscosity()),
      smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
