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
  * @file relax_dynamics.h
  * @brief Particle relax is the process to produce body fitted
  * particle distribution with zero-order consistency.
  * @author	Bo Zhang, Xiangyu Hu
  */

#ifndef RELAX_DYNAMICS_H
#define RELAX_DYNAMICS_H

#include "base_relax_dynamics.h"
#include "general_constraint.h"

namespace SPH
{
class GemoteryShape;
class LevelSetShape;

namespace relax_dynamics
{
template <class DataDelegationType>
class BaseIntegration : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseIntegration(BaseRelationType &base_relation);
    virtual ~BaseIntegration(){};

  protected:
    SPHAdaptation *sph_adaptation_;
    StdLargeVec<Vecd> &residue_;
};

template <typename... RelaxtionTypes>
class RelaxationIntegration;

template <class KernelCorrectionType>
class RelaxationIntegration<Inner<>, KernelCorrectionType>
    : public BaseIntegration<RelaxDataDelegateInner>
{
  public:
    explicit RelaxationIntegration(BaseInnerRelation &inner_relation);
    RelaxationIntegration(BaseInnerRelation &inner_relation, std::string sub_shape_name);
    virtual ~RelaxationIntegration(){};
    Shape &getRelaxShape() { return relax_shape_; };
    void interaction(size_t index_i, Real dt = 0.0);
    void update£¨size_t index_i, Real dt = 0.0);

  protected:
    Shape &relax_shape_;
    KernelCorrectionType correction_;
};

template <class KernelCorrectionType>
class RelaxationIntegration<Inner<LevelSetCorrection>, KernelCorrectionType> : 
    public RelaxationIntegration<Inner<>, KernelCorrectionType>
{
  public:
    template <typename... Args>
    RelaxationIntegration(Args &&...args);
    template <typename BodyRelationType, typename FirstArg>
    explicit RelaxationIntegration(ConstructorArgs<BodyRelationType, FirstArg> parameters)
        : RelaxationIntegration(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~RelaxationIntegration(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_;
    LevelSetShape &level_set_shape_;
};

template <class KernelCorrectionType>
class RelaxationIntegration<Contact<>, KernelCorrectionType>
    : public BaseIntegration<RelaxDataDelegateContact>
{
  public:
    explicit RelaxationResidue(BaseContactRelation &contact_relation)
        : RelaxationResidue<Base, RelaxDataDelegateContact>(contact_relation){};
    virtual ~RelaxationResidue(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType correction_;
    StdLargeVecd<KernelCorrectionType> contact_corrections_;
};

using RelaxationInner = RelaxationIntegration<Inner<>, NoKernelCorrection>;
using RelaxationWithLevelSetInner = RelaxationIntegration<Inner<LevelSetCorrection>, NoKernelCorrection>;
using RelaxationComplex = ComplexInteraction<RelaxationIntegration<Inner<>, Contact<>>, NoKernelCorrection>;
using RelaxationWithLevelSetComplex = ComplexInteraction<RelaxationIntegration<Inner<LevelSetCorrection>, Contact<>>, NoKernelCorrection>;
using RelaxationCorrectionInner = RelaxationIntegration<Inner<>, KernelCorrection>;
using RelaxationCorrectionWithLevelSetInner= RelaxationIntegration<Inner<LevelSetCorrection>, KernelCorrection>;
using RelaxationCorrectionComplex = ComplexInteraction<RelaxationIntegration<Inner<>, Contact<>>, KernelCorrection>;
using RelaxationCorrectionWithLevelSetComplex = ComplexInteraction<RelaxationIntegration<Inner<LevelSetCorrection>, Contact<>>, KernelCorrection>;
}
#endif // RELAX_DYNAMICS_H