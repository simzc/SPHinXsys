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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
#ifndef DIFFUSION_DYNAMICS_IMPLICIT_H
#define DIFFUSION_DYNAMICS_IMPLICIT_H

#include "particle_dynamics_dissipation.h"

#include "general_operators.h"

namespace SPH
{
/**
 * @class BaseDiffusionSplittingInner
 * @brief A quantity damping by splitting scheme
 * this method modifies the quantity directly.
 * Note that, if periodic boundary condition is applied,
 * the parallelized version of the method requires the one using ghost particles
 * because the splitting partition only works in this case.
 */
template <typename DataType, class CoefficientType>
class BaseDiffusionSplittingInner : public OperatorInner<DataType, DataType, CoefficientType>
{
  public:
    template <typename CoefficientArg>
    BaseDiffusionSplittingInner(BaseInnerRelation &inner_relation,
                                const std::string &variable_name, const CoefficientArg &eta);
    virtual ~BaseDiffusionSplittingInner(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &Vol_, &mass_;
    StdLargeVec<DataType> &variable_;

    virtual ErrorAndParameters<DataType> computeErrorAndParameters(size_t index_i, Real dt);
    virtual void updateStates(size_t index_i, Real dt, const ErrorAndParameters<DataType> &error_and_parameters);
};

template <typename DataType>
using DiffusionSplittingInner = BaseDiffusionSplittingInner<DataType, ConstantCoefficient<Real>>;
template <typename DataType>
using DiffusionSplittingInnerCoefficientByParticle = BaseDiffusionSplittingInner<DataType, CoefficientByParticle<Real>>;

/**
 * @class BaseDiffusionSplittingWithWall
 * @brief Note that the values on the wall are constrained.
 */
template <typename DataType, class CoefficientType>
class BaseDiffusionSplittingWithWall : public BaseDiffusionSplittingInner<DataType, CoefficientType>, public DissipationDataWithWall
{
  public:
    template <typename CoefficientArg>
    BaseDiffusionSplittingWithWall(ComplexRelation &complex_wall_relation,
                                   const std::string &variable_name, const CoefficientArg &eta);
    virtual ~BaseDiffusionSplittingWithWall(){};

  protected:
    virtual ErrorAndParameters<DataType> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;

  private:
    StdVec<StdLargeVec<DataType> *> wall_variable_;
};

template <typename DataType>
using DiffusionSplittingWithWall = BaseDiffusionSplittingWithWall<DataType, ConstantCoefficient<Real>>;
template <typename DataType>
using DiffusionSplittingWithWallCoefficientByParticle = BaseDiffusionSplittingWithWall<DataType, CoefficientByParticle<Real>>;

} // namespace SPH
#endif // DIFFUSION_DYNAMICS_IMPLICIT_H