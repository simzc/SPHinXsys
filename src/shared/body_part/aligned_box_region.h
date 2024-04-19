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
 * @file 	aligned_box_region.h
 * @brief 	This is the base classes of body parts.
 * @details	There two main type of body parts. One is part by particle.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef ALIGNED_BOX_REGION_H
#define ALIGNED_BOX_REGION_H

#include "body_part.h"
#include "body_part_by_cell.h"
#include "body_part_by_particle.h"

namespace SPH
{
/**
 * @class AlignedBoxRegion
 * @brief A template body part with the collection of particles within by an AlignedBoxShape.
 */
template <class BodyRegionType>
class AlignedBoxRegion : public BodyRegionType
{
  public:
    AlignedBoxShape &aligned_box_;

    AlignedBoxRegion(RealBody &real_body, SharedPtr<AlignedBoxShape> aligned_box_ptr)
        : BodyRegionType(real_body, aligned_box_ptr), aligned_box_(*aligned_box_ptr.get()){};
    virtual ~AlignedBoxRegion(){};
};

using BodyAlignedBoxByParticle = AlignedBoxRegion<BodyRegionByParticle>;
using BodyAlignedBoxByCell = AlignedBoxRegion<BodyRegionByCell>;
} // namespace SPH
#endif // ALIGNED_BOX_REGION_H
