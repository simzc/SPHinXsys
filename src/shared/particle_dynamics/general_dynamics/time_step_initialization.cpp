#include "time_step_initialization.h"

namespace SPH
{
//=================================================================================================//
ForcePrior::ForcePrior(SPHBody &sph_body, const std::string &force_name)
    : LocalDynamics(sph_body), force_prior_(base_particles_.force_prior_),
      force_(*base_particles_.registerSharedVariable<Vecd>(force_name)),
      previous_force_(*base_particles_.registerSharedVariable<Vecd>("Previous" + force_name))
{
    base_particles_.addVariableToRestart<Vecd>("Previous" + force_name);
}
//=================================================================================================//
void ForcePrior::update(size_t index_i, Real dt)
{
    force_prior_[index_i] += force_[index_i] - previous_force_[index_i];
    previous_force_[index_i] = force_[index_i];
}
//=================================================================================================//
GravityForce::GravityForce(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
    : ForcePrior(sph_body, "GravityForce"),
      pos_(base_particles_.pos_), mass_(base_particles_.mass_),
      gravity_(gravity_ptr_keeper_.assignPtr(gravity_ptr)) {}
//=================================================================================================//
void GravityForce::update(size_t index_i, Real dt)
{
    force_[index_i] = mass_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
    ForcePrior::update(index_i, dt);
}
//=================================================================================================//
TimeStepInitialization::TimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
    : BaseTimeStepInitialization(sph_body, gravity_ptr), GeneralDataDelegateSimple(sph_body),
      pos_(particles_->pos_), force_prior_(particles_->force_prior_), mass_(particles_->mass_) {}
//=================================================================================================//
void TimeStepInitialization::update(size_t index_i, Real dt)
{
    force_prior_[index_i] = mass_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
}
//=================================================================================================//
RandomizeParticlePosition::RandomizeParticlePosition(SPHBody &sph_body)
    : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
      pos_(particles_->pos_), randomize_scale_(sph_body.sph_adaptation_->MinimumSpacing()) {}
//=================================================================================================//
void RandomizeParticlePosition::update(size_t index_i, Real dt)
{
    Vecd &pos_n_i = pos_[index_i];
    for (int k = 0; k < pos_n_i.size(); ++k)
    {
        pos_n_i[k] += dt * rand_uniform(-1.0, 1.0) * randomize_scale_;
    }
}
//=================================================================================================//
} // namespace SPH
