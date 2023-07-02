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
#ifndef GENERAL_DYNAMICS_EXTRA_H
#define GENERAL_DYNAMICS_EXTRA_H

#include "general_dynamics.h"

namespace SPH
{
/**
 * @class ValueAssignment
 * @brief set value for a discrete variable
 * @details Note that this class only prepares the data,
 * concrete realization wll be defined in application.
 */
template <typename DataType>
class ValueAssignment : public LocalDynamics, public GeneralDataDelegateSimple
{
  public:
    ValueAssignment(SPHBody &sph_body, const std::string &variable_name)
        : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
          variable_(*particles_->getVariableByName<DataType>(variable_name)){};
    virtual ~ValueAssignment(){};

  protected:
    StdLargeVec<DataType> &variable_;
};

/**
 * @class ImposingSourceTerm
 * @brief set source effect to a discrete variable
 */
template <typename DataType>
class ImposingSourceTerm : public LocalDynamics, public GeneralDataDelegateSimple
{
  public:
    ImposingSourceTerm(SPHBody &sph_body, const std::string &variable_name, const DataType &source_strength)
        : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
          variable_(*particles_->getVariableByName<DataType>(variable_name)),
          source_strength_(source_strength){};
    virtual ~ImposingSourceTerm(){};
    void setSourceStrength(Real source_strength) { source_strength_ = source_strength; };
    void update(size_t index_i, Real dt)
    {
        variable_[index_i] += source_strength_ * dt;
    };

  protected:
    StdLargeVec<DataType> &variable_;
    DataType source_strength_;
};

class ConstraintTotalScalarAmount : public LocalDynamics, public GeneralDataDelegateSimple
{
  public:
    ConstraintTotalScalarAmount(SPHBody &sph_body, const std::string &variable_name);
    virtual ~ConstraintTotalScalarAmount(){};
    void setupInitialScalarAmount();
    void setupDynamics(Real dt = 0.0) override;
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &variable_;
    ReduceDynamics<QuantityMoment<Real>> total_scalar_;
    bool is_initialized_;
    Real inital_total_;
    Real increment_;
};

template <typename Datatype>
inline Real getSquaredNorm(const Datatype &variable)
{
    return variable.squaredNorm();
};

template <>
inline Real getSquaredNorm<Real>(const Real &variable)
{
    return variable * variable;
};

/**
 * @class MaximumNorm
 * @brief  obtained the maximum norm of a variable
 */
template <typename DataType>
class MaximumNorm : public LocalDynamicsReduce<Real, ReduceMax>,
                    public GeneralDataDelegateSimple
{
  public:
    MaximumNorm(SPHBody &sph_body, const std::string &variable_name)
        : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
          GeneralDataDelegateSimple(sph_body),
          variable_(*particles_->getVariableByName<DataType>(variable_name)){};
    virtual ~MaximumNorm(){};

    virtual Real outputResult(Real reduced_value) override { return std::sqrt(reduced_value); }
    Real reduce(size_t index_i, Real dt = 0.0) { return getSquaredNorm(variable_[index_i]); };

  protected:
    StdLargeVec<DataType> &variable_;
};
} // namespace SPH
#endif // GENERAL_DYNAMICS_EXTRA_H
