/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file bidirectional_buffer.h
 * @brief Particle buffer for implementing velocity and pressure boundary condition in channel flows.
 * @details The buffer generates particles for inlet and delete particles for outlet.
 * Note that the buffer is able to generate and delete particles at the same buffer surface.
 * @author	Shuoguo Zhang and Xiangyu Hu
 */

#ifndef BIDIRECTIONAL_BUFFER_H
#define BIDIRECTIONAL_BUFFER_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
class BidirectionalBuffer
{
  protected:
    ConcurrentIndexVector buffer_particle_list_;
    virtual Real getTargetPressure(Real dt) = 0;

    class TagBufferParticles : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
    {
      public:
        TagBufferParticles(BodyAlignedBoxByCell &aligned_box_part, int axis);
        virtual ~TagBufferParticles(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
        AlignedBoxShape &aligned_box_;
        int axis_;
        StdLargeVec<Vecd> &pos_;
        StdLargeVec<int> &buffer_particle_indicator_;
    };

    class Injection : public BaseLocalDynamics<BodyPartByCell>, public FluidDataSimple
    {
      public:
        Injection(BodyAlignedBoxByCell &aligned_box_part, ParticleBuffer<Base> &particle_buffer, int axis)
            : BaseLocalDynamics<BodyPartByCell>(aligned_box_part),
              FluidDataSimple(aligned_box_part.getSPHBody()),
              aligned_box_(aligned_box_part.aligned_box_),
              particle_buffer_(particle_buffer), axis_(axis),
              fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
              pos_n_(particles_->pos_), rho_n_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
              previous_surface_indicator_(*particles_->getVariableByName<int>("PreviousSurfaceIndicator")),
              buffer_particle_indicator_(*particles_->getVariableByName<int>("BufferParticleIndicator")){};
        virtual ~Injection(){};

        /** This class is only implemented in sequential due to memory conflicts. */
        virtual void update(size_t index_i, Real dt = 0.0);

      protected:
        std::mutex mutex_switch_to_real_;
        AlignedBoxShape &aligned_box_;
        Fluid &fluid_;
        StdLargeVec<Vecd> &pos_n_;
        StdLargeVec<Real> &rho_n_, &p_;
        ParticleBuffer<Base> &particle_buffer_;
        const int axis_;
        StdLargeVec<int> &previous_surface_indicator_, &buffer_particle_indicator_;
    };

  public:
    BidirectionalBuffer(RealBody &real_body, SharedPtr<AlignedBoxShape> shape_ptr,
                        ParticleBuffer<Base> &particle_buffer, int axis_direction)
        : tag_buffer_particles(real_body, shape_ptr),
          injection(real_body, shape_ptr, particle_buffer, axis_direction)
    {
        particle_buffer.checkParticlesReserved();
    };
    virtual ~BidirectionalBuffer(){};

    SimpleDynamics<TagBufferParticles> tag_buffer_particles;
    SimpleDynamics<Injection> injection;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // BIDIRECTIONAL_BUFFER_H