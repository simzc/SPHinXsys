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
 * @file 	body_part.h
 * @brief 	This is the base classes of body parts.
 * @details	There two main type of body parts. One is part by particle.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BODY_PART_H
#define BODY_PART_H

#include "base_body.h"

namespace SPH
{
using namespace std::placeholders;

/**
 * @class BodyPart
 * @brief An auxillary class for SPHBody to indicate a part of the body.
 */
class BodyPart
{
  public:
    BodyPart(SPHBody &sph_body, const std::string &body_part_name)
        : sph_body_(sph_body), body_part_name_(body_part_name){};
    virtual ~BodyPart(){};

    SPHBody &getSPHBody() { return sph_body_; };
    std::string getName() { return body_part_name_; };

  protected:
    SPHBody &sph_body_;
    std::string body_part_name_;
};
} // namespace SPH
#endif // BODY_PART_H
