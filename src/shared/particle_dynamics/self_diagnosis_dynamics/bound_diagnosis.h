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
 * @file 	bound_diagnosis.h
 * @brief 	diagnosis related bounding check
 * @author	Xiangyu Hu
 */

#ifndef BOUND_DIAGNOSIS_H
#define BOUND_DIAGNOSIS_H

#include "general_dynamics.h"

namespace SPH
{
/**
 * @class VelocityBoundCheck
 * @brief  check whether particle velocity within a given bound
 */
class VelocityBoundCheck : public LocalDynamicsReduce<bool, ReduceOR>,
                           public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &vel_;
    Real velocity_bound_;

  public:
    VelocityBoundCheck(SPHBody &sph_body, Real velocity_bound);
    virtual ~VelocityBoundCheck(){};

    bool reduce(size_t index_i, Real dt = 0.0);
};
} // namespace SPH
#endif // BOUND_DIAGNOSIS_H