#include "bidirectional_buffer.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
BidirectionalBuffer::TagBufferParticles::
    TagBufferParticles(BodyAlignedBoxByCell &aligned_box_part, int axis)
    : BaseLocalDynamics<BodyPartByCell>(aligned_box_part),
      FluidDataSimple(aligned_box_part.getSPHBody()),
      aligned_box_(aligned_box_part.aligned_box_), axis_(axis), pos_(particles_->pos_),
      buffer_particle_indicator_(*particles_->registerSharedVariable<int>("BufferParticleIndicator"))
{
    particles_->registerSortableVariable<int>("BufferParticleIndicator");
}
//=================================================================================================//
void BidirectionalBuffer::TagBufferParticles::update(size_t index_i, Real dt = 0.0)
{
    buffer_particle_indicator_[index_i] = 0;
    if (aligned_box_.checkInBounds(axis_, pos_[index_i]))
    {
        buffer_particle_indicator_[index_i] = 1;
    }
}
//=================================================================================================//
void BidirectionalBuffer::Injection::update(size_t index_i, Real dt)
{
    if (buffer_particle_indicator_[index_i] == 1)
    {
        if (aligned_box_.checkUpperBound(axis_, pos_n_[index_i]))
        {
            mutex_switch_to_real_.lock();
            particle_buffer_.checkEnoughBuffer(*particles_);
            buffer_particle_indicator_[index_i] = 0;
            particles_->copyFromAnotherParticle(particles_->total_real_particles_, index_i);
            particles_->total_real_particles_ += 1;
            mutex_switch_to_real_.unlock();

            pos_n_[index_i] = aligned_box_.getUpperPeriodic(axis_, pos_n_[index_i]);
            p_[index_i] = buffer_.getTargetPressure(dt);
            rho_n_[index_i] = fluid_.DensityFromPressure(p_[index_i]);
            previous_surface_indicator_[index_i] = 1;
        }
    }
}
//=====================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
