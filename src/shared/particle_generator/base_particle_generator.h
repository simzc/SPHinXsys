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
 * @file base_particle_generator.h
 * @brief This is the base class of particle generator, which generates particles
 * with given positions and volumes. The direct generator simply generate
 * particle with given position and volume.
 * The particle generators naming are using template partial specialization,
 * with three keywords. The first indicates volume (default), surface and line particles,
 * the second generating methods, and the third the control parameters,
 * such as adaptive for adaptive resolution.
 * There are geometric dynamically allocating variables for particle generation.
 * After particle generation, their values are copied to corresponding
 * fix allcoated particle variables for numerical simulations.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_PARTICLE_GENERATOR_H
#define BASE_PARTICLE_GENERATOR_H

#include "base_data_package.h"
#include "large_data_containers.h"
#include "sph_data_containers.h"

namespace SPH
{

class SPHBody;
class BaseParticles;
//---------------------------------------------------------------------------
// Geometric types of particles. The default type is volume metric particles.
//---------------------------------------------------------------------------
class Surface;
class ThickSurface; // Surface thickness equal or larger than the particle spacing
class Observer;
class Reload;

template <typename... Parameters>
class ParticleGenerator;

template <typename... Parameters>
class GeneratingMethod;

template <> // default volume metric particle generator
class ParticleGenerator<Base>
{
  private:
    ValueAddressAssemble<StdLargeVec> geometric_data_;
    VariableAddressAssemble<StdLargeVec> geometric_variables_;
    VariableUniquePtrAssemble<StdLargeVec> geometric_variables_ptrs_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body);
    virtual ~ParticleGenerator(){};
    virtual void initializeGeometricVariables() = 0;
    void generateParticlesWithBasicVariables();

  protected:
    StdLargeVec<Vecd> &position_;
    StdLargeVec<Real> &volumetric_measure_;
    BaseParticles &base_particles_;
    Real particle_spacing_ref_;

    template <typename DataType>
    StdLargeVec<DataType> *registerGeometricVariable(const std::string &name);
    virtual void initializePosition(const Vecd &position);
    virtual void initializePositionAndVolumetricMeasure(const Vecd &position, Real volumetric_measure);
};

template <> // generate surface particles
class ParticleGenerator<Surface> : public ParticleGenerator<Base>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body);
    virtual ~ParticleGenerator(){};

  protected:
    //----------------------------------------------------------------------
    // Geometric dynamically allocating variables for particle generation
    //----------------------------------------------------------------------
    StdLargeVec<Vecd> &surface_normal_;
    StdLargeVec<Real> &surface_thickness_;
    virtual void initializeSurfaceProperties(const Vecd &surface_normal, Real thickness);
};

template <> // generate observer particles
class ParticleGenerator<Observer> : public ParticleGenerator<Base>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body)
        : ParticleGenerator<Base>(sph_body){};
    ParticleGenerator(SPHBody &sph_body, const StdVec<Vecd> &positions);
    virtual ~ParticleGenerator(){};
    virtual void initializeGeometricVariables() override;
};

template <> // generate particles by reloading dynamically relaxed particles
class ParticleGenerator<Reload> : public ParticleGenerator<Base>
{
    BaseMaterial &base_material_;
    std::string file_path_;

  public:
    ParticleGenerator(SPHBody &sph_body, const std::string &reload_body_name);
    virtual ~ParticleGenerator(){};
    virtual void initializeGeometricVariables() override;
};

} // namespace SPH
#endif // BASE_PARTICLE_GENERATOR_H
