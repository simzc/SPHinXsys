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
 * @file 	body_part_by_particle.h
 * @brief 	This is the base classes of body parts.
 * @details	There two main type of body parts. One is part by particle.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BODY_PART_BY_PARTICLE_H
#define BODY_PART_BY_PARTICLE_H

#include "body_part.h"

namespace SPH
{
/**
 * @class BodyPartByParticle
 * @brief A body part with a collection of particles.
 */
class BodyPartByParticle : public BodyPart
{
  public:
    IndexVector body_part_particles_; /**< Collection particle in this body part. */
    BaseParticles &getBaseParticles() { return base_particles_; };
    IndexVector &LoopRange() { return body_part_particles_; };
    size_t SizeOfLoopRange() { return body_part_particles_.size(); };

    BodyPartByParticle(SPHBody &sph_body, const std::string &body_part_name)
        : BodyPart(sph_body, body_part_name), base_particles_(sph_body.getBaseParticles()),
          body_part_bounds_(Vecd::Zero(), Vecd::Zero()), body_part_bounds_set_(false){};
    virtual ~BodyPartByParticle(){};

    void setBodyPartBounds(BoundingBox bbox)
    {
        body_part_bounds_ = bbox;
        body_part_bounds_set_ = true;
    };

    BoundingBox getBodyPartBounds()
    {
        if (!body_part_bounds_set_)
            std::cout << "WARNING: the body part bounds are not set for BodyPartByParticle." << std::endl;
        return body_part_bounds_;
    }

  protected:
    BaseParticles &base_particles_;
    BoundingBox body_part_bounds_;
    bool body_part_bounds_set_;

    typedef std::function<void(size_t)> TaggingParticleMethod;
    void tagParticles(TaggingParticleMethod &tagging_particle_method);
};

/**
 * @class BodyRegionByParticle
 * @brief A  body part with the collection of particles within by a prescribed shape.
 */
class BodyRegionByParticle : public BodyPartByParticle
{
  private:
    SharedPtrKeeper<Shape> shape_ptr_keeper_;

  public:
    BodyRegionByParticle(SPHBody &sph_body, Shape &body_part_shape);
    BodyRegionByParticle(SPHBody &sph_body, SharedPtr<Shape> shape_ptr);
    virtual ~BodyRegionByParticle(){};
    Shape &getBodyPartShape() { return body_part_shape_; };

  private:
    Shape &body_part_shape_;
    void tagByContain(size_t particle_index);
};

/**
 * @class BodySurface
 * @brief A  body part with the collection of particles at surface of a body
 */
class BodySurface : public BodyPartByParticle
{
  public:
    explicit BodySurface(SPHBody &sph_body);
    virtual ~BodySurface(){};

  private:
    Real particle_spacing_min_;
    void tagNearSurface(size_t particle_index);
};

/**
 * @class BodySurfaceLayer
 * @brief A  body part with the collection of particles within the surface layers of a body.
 */
class BodySurfaceLayer : public BodyPartByParticle
{
  public:
    explicit BodySurfaceLayer(SPHBody &sph_body, Real layer_thickness = 3.0);
    virtual ~BodySurfaceLayer(){};

  private:
    Real thickness_threshold_;
    void tagSurfaceLayer(size_t particle_index);
};

} // namespace SPH
#endif // BODY_PART_BY_PARTICLE_H
