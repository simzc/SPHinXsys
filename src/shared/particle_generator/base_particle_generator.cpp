#include "base_particle_generator.hpp"

#include "base_body.h"
#include "base_particles.h"
#include "io_all.h"

namespace SPH
{
//=================================================================================================//
ParticleGenerator<Base>::ParticleGenerator(SPHBody &sph_body)
    : position_(*registerGeometricVariable<Vecd>("Position")),
      volumetric_measure_(*registerGeometricVariable<Real>("VolumetricMeasure")),
      base_particles_(sph_body.getBaseParticles()),
      particle_spacing_ref_(sph_body.sph_adaptation_->ReferenceSpacing()) {}
//=================================================================================================//
void ParticleGenerator<Base>::initializePosition(const Vecd &position)
{
    position_.push_back(position);
    unsorted_id_.push_back(base_particles_.total_real_particles_);
    base_particles_.total_real_particles_++;
}
//=================================================================================================//
void ParticleGenerator<Base>::generateParticlesWithBasicVariables()
{
    initializeGeometricVariables();
    base_particles_.initializeAllParticlesBounds();
}
//=================================================================================================//
void ParticleGenerator<Base>::initializePositionAndVolumetricMeasure(
    const Vecd &position, Real volumetric_measure)
{
    initializePosition(position);
    volumetric_measure_.push_back(volumetric_measure);
}
//=================================================================================================//
ParticleGenerator<Surface>::ParticleGenerator(SPHBody &sph_body)
    : ParticleGenerator<Base>(sph_body),
      surface_normal_(*registerGeometricVariable<Vecd>("Normal")),
      surface_thickness_(*registerGeometricVariable<Real>("Thickness")),
      n_(*base_particles_.getVariableByName<Vecd>("NormalDirection")),
      thickness_(*base_particles_.getVariableByName<Real>("Thickness")) {}
//=================================================================================================//
void ParticleGenerator<Surface>::initializeSurfaceProperties(const Vecd &surface_normal, Real thickness)
{
    n_.push_back(surface_normal);
    thickness_.push_back(thickness);
}
//=================================================================================================//
ParticleGenerator<Observer>::ParticleGenerator(SPHBody &sph_body, const StdVec<Vecd> &positions)
    : ParticleGenerator<Base>(sph_body)
{
    for (size_t i = 0; i < positions.size(); ++i)
    {
        initializePositionAndVolumetricMeasure(positions[i], 0.0);
    }
}
//=================================================================================================//
void ParticleGenerator<Observer>::initializeGeometricVariables()
{
}
//=================================================================================================//
ParticleGenerator<Reload>::ParticleGenerator(SPHBody &sph_body, const std::string &reload_body_name)
    : ParticleGenerator<Base>(sph_body), base_material_(sph_body.getBaseMaterial())
{
    std::string reload_folder = sph_body.getSPHSystem().getIOEnvironment().reload_folder_;
    if (!fs::exists(reload_folder))
    {
        std::cout << "\n Error: the particle reload folder:" << reload_folder << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    file_path_ = reload_folder + "/" + reload_body_name + "_rld.xml";
}
//=================================================================================================//
void ParticleGenerator<Reload>::initializeGeometricVariables()
{
    base_material_.registerReloadLocalParameters(&base_particles_);
    base_particles_.readFromXmlForReloadParticle(file_path_);
}
//=================================================================================================//
} // namespace SPH
