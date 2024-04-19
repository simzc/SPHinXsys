#include "body_part_by_particle.h"

#include "base_particles.hpp"
namespace SPH
{
//=================================================================================================//
void BodyPartByParticle::tagParticles(TaggingParticleMethod &tagging_particle_method)
{
    for (size_t i = 0; i < base_particles_.total_real_particles_; ++i)
    {
        tagging_particle_method(i);
    }
};
//=================================================================================================//
BodyRegionByParticle::
    BodyRegionByParticle(SPHBody &sph_body, Shape &body_part_shape)
    : BodyPartByParticle(sph_body, body_part_shape.getName()),
      body_part_shape_(body_part_shape)
{
    TaggingParticleMethod tagging_particle_method = std::bind(&BodyRegionByParticle::tagByContain, this, _1);
    tagParticles(tagging_particle_method);
}
//=================================================================================================//
BodyRegionByParticle::BodyRegionByParticle(SPHBody &sph_body, SharedPtr<Shape> shape_ptr)
    : BodyRegionByParticle(sph_body, *shape_ptr.get())
{
    shape_ptr_keeper_.assignRef(shape_ptr);
}
//==
//=================================================================================================//
void BodyRegionByParticle::tagByContain(size_t particle_index)
{
    if (body_part_shape_.checkContain(base_particles_.pos_[particle_index]))
    {
        body_part_particles_.push_back(particle_index);
    }
}
//=================================================================================================//
BodySurface::BodySurface(SPHBody &sph_body)
    : BodyPartByParticle(sph_body, "BodySurface"),
      particle_spacing_min_(sph_body.sph_adaptation_->MinimumSpacing())
{
    TaggingParticleMethod tagging_particle_method = std::bind(&BodySurface::tagNearSurface, this, _1);
    tagParticles(tagging_particle_method);
    std::cout << "Number of surface particles : " << body_part_particles_.size() << std::endl;
}
//=================================================================================================//
void BodySurface::tagNearSurface(size_t particle_index)
{
    Real phi = sph_body_.getInitialShape().findSignedDistance(base_particles_.pos_[particle_index]);
    if (fabs(phi) < particle_spacing_min_)
        body_part_particles_.push_back(particle_index);
}
//=================================================================================================//
BodySurfaceLayer::BodySurfaceLayer(SPHBody &sph_body, Real layer_thickness)
    : BodyPartByParticle(sph_body, "InnerLayers"),
      thickness_threshold_(sph_body.sph_adaptation_->ReferenceSpacing() * layer_thickness)
{
    TaggingParticleMethod tagging_particle_method = std::bind(&BodySurfaceLayer::tagSurfaceLayer, this, _1);
    tagParticles(tagging_particle_method);
    std::cout << "Number of inner layers particles : " << body_part_particles_.size() << std::endl;
}
//=================================================================================================//
void BodySurfaceLayer::tagSurfaceLayer(size_t particle_index)
{
    Real distance = fabs(sph_body_.getInitialShape().findSignedDistance(base_particles_.pos_[particle_index]));
    if (distance < thickness_threshold_)
    {
        body_part_particles_.push_back(particle_index);
    }
}
//=================================================================================================//
} // namespace SPH
