#include "force_prior.h"

namespace SPH
{
//=================================================================================================//
ForcePrior::ForcePrior(BaseParticles *base_particles, const std::string &force_name)
    : force_prior_(base_particles->force_prior_),
      force_(*base_particles->registerSharedVariable<Vecd>(force_name)),
      previous_force_(*base_particles->registerSharedVariable<Vecd>("Previous" + force_name))
{
    base_particles->addVariableToRestart<Vecd>("Previous" + force_name);
}
//=================================================================================================//
void ForcePrior::update(size_t index_i, Real dt)
{
    force_prior_[index_i] += force_[index_i] - previous_force_[index_i];
    previous_force_[index_i] = force_[index_i];
}
//=================================================================================================//
GravityForce::GravityForce(SPHBody &sph_body, Gravity &gravity)
    : LocalDynamics(sph_body), ForcePrior(&base_particles_, "GravityForce"),
      GeneralDataDelegateSimple(sph_body), gravity_(gravity),
      pos_(base_particles_.pos_), mass_(base_particles_.mass_) {}
//=================================================================================================//
void GravityForce::update(size_t index_i, Real dt)
{
    force_[index_i] = mass_[index_i] * gravity_.InducedAcceleration(pos_[index_i]);
    ForcePrior::update(index_i, dt);
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
