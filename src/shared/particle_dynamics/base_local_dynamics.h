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
 * @file    base_local_dynamics.h
 * @brief 	This is for the base classes of local particle dynamics, which describe the
 * 			dynamics of a particle and it neighbors.
 * @author	Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef BASE_LOCAL_DYNAMICS_H
#define BASE_LOCAL_DYNAMICS_H

#include "base_data_package.h"
#include "base_particle_dynamics.h"
#include "particle_functors.h"
#include "sph_data_containers.h"

namespace SPH
{
//----------------------------------------------------------------------
// Interaction type identifies
//----------------------------------------------------------------------
template <typename... InnerParameters>
class Inner {}; /**< Inner interaction: interaction within a body*/

template <typename... ContactParameters>
class Contact {}; /**< Contact interaction: interaction between a body with one or several another bodies */

class Boundary {};        /**< Interaction with boundary */
class Wall {};            /**< Interaction with wall boundary */
class Extended {};        /**< An extened method of an interaction type */
class SpatialTemporal {}; /**< A interaction considering spatial temporal correlations */
class Dynamic {};         /**< A dynamic interaction */

/**
 * @class BaseLocalDynamics
 * @brief The base class for all local particle dynamics.
 */
template <class DynamicsIdentifier>
class BaseLocalDynamics
{
  public:
    explicit BaseLocalDynamics(DynamicsIdentifier &identifier)
        : identifier_(identifier){};
    virtual ~BaseLocalDynamics(){};
    SPHBody &getSPHBody() { return sph_body_; };
    DynamicsIdentifier &getDynamicsIdentifier() { return identifier_; };
    virtual void setupDynamics(Real dt = 0.0){}; // setup global parameters
  protected:
    DynamicsIdentifier &identifier_;
    SPHBody &sph_body_ = identifier_.getSPHBody();
    BaseParticles &base_particles_ = sph_body_.getBaseParticles();
};
using LocalDynamics = BaseLocalDynamics<SPHBody>;

/**
 * @class BaseLocalDynamicsReduce
 * @brief The base class for all local particle dynamics for reducing.
 */
template <typename ReturnType, typename Operation, class DynamicsIdentifier>
class BaseLocalDynamicsReduce : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    BaseLocalDynamicsReduce(DynamicsIdentifier &identifier, ReturnType reference)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier), reference_(reference),
          quantity_name_("ReducedQuantity"){};
    virtual ~BaseLocalDynamicsReduce(){};

    using ReduceReturnType = ReturnType;
    ReturnType Reference() { return reference_; };
    std::string QuantityName() { return quantity_name_; };
    Operation &getOperation() { return operation_; };
    virtual ReturnType outputResult(ReturnType reduced_value) { return reduced_value; }

  protected:
    ReturnType reference_;
    Operation operation_;
    std::string quantity_name_;
};
template <typename ReturnType, typename Operation>
using LocalDynamicsReduce = BaseLocalDynamicsReduce<ReturnType, Operation, SPHBody>;

/**
 * @class Average
 * @brief Derives class for computing particle-wise averages.
 */
template <class ReduceSumType>
class Average : public ReduceSumType
{
  public:
    template <class DynamicsIdentifier, typename... Args>
    Average(DynamicsIdentifier &identifier, Args &&...args)
        : ReduceSumType(identifier, std::forward<Args>(args)...){};
    virtual ~Average(){};
    using ReturnType = typename ReduceSumType::ReduceReturnType;

    virtual ReturnType outputResult(ReturnType reduced_value)
    {
        ReturnType sum = ReduceSumType::outputResult(reduced_value);
        return sum / Real(this->getDynamicsIdentifier().SizeOfLoopRange());
    }
};

/**
 * @class ConstructorArgs
 * @brief Class template argument deduction (CTAD) for constructor arguments.
 */
template <typename BodyRelationType, typename... OtherArgs>
struct ConstructorArgs
{
    BodyRelationType &body_relation_;
    std::tuple<OtherArgs...> others_;
    SPHBody &getSPHBody() { return body_relation_.getSPHBody(); };
    ConstructorArgs(BodyRelationType &body_relation, OtherArgs... other_args)
        : body_relation_(body_relation),
          others_(other_args...){};
};

/**
 * @class ComplexInteraction
 * @brief A class that integrates multiple local dynamics.
 * Typically, it includes an inner interaction and one or
 * several contact interaction and boundary conditions.
 */
template <typename... T>
class ComplexInteraction;

template <typename... T>
class ComplexInteractionKernel;

template <typename... CommonParameters, template <typename... InteractionTypes> class LocalDynamicsName>
class ComplexInteractionKernel<LocalDynamicsName<CommonParameters...>>
{
  public:
    ComplexInteractionKernel(){};

    void interaction(size_t index_i, Real dt = 0.0){};
};

template <typename... CommonParameters, template <typename... InteractionTypes> class LocalDynamicsName,
          class FirstInteraction, class... OtherInteractions>
class ComplexInteractionKernel<LocalDynamicsName<CommonParameters...>, FirstInteraction, OtherInteractions...>
    : public decltype(LocalDynamicsName<FirstInteraction, CommonParameters...>::device_kernel)::KernelType
{
    using LocalDynamicsNameKernel = typename decltype(LocalDynamicsName<FirstInteraction, CommonParameters...>::device_kernel)::KernelType;

  protected:
//    LocalDynamicsNameKernel *first_interaction_;
    ComplexInteractionKernel<LocalDynamicsName<CommonParameters...>, OtherInteractions...> *other_interactions_;

  public:
    explicit ComplexInteractionKernel(LocalDynamicsNameKernel *first_interaction_kernel,
                                      ComplexInteractionKernel<LocalDynamicsName<CommonParameters...>, OtherInteractions...> *other_interactions_kernel)
        : LocalDynamicsNameKernel(*first_interaction_kernel),
          other_interactions_(other_interactions_kernel) {}

    void interaction(size_t index_i, Real dt = 0.0)
    {
//        first_interaction_->interaction(index_i, dt);
        LocalDynamicsNameKernel::interaction(index_i, dt);
        other_interactions_->interaction(index_i, dt);
    };
};

template <typename... CommonParameters, template <typename... InteractionTypes> class LocalDynamicsName>
class ComplexInteraction<LocalDynamicsName<>, CommonParameters...>
{
  public:
    ComplexInteraction(){};

    void interaction(size_t index_i, Real dt = 0.0){};

    execution::DeviceImplementation<ComplexInteractionKernel<LocalDynamicsName<CommonParameters...>>> device_kernel;
};

template <typename... CommonParameters, template <typename... InteractionTypes> class LocalDynamicsName,
          class FirstInteraction, class... OtherInteractions>
class ComplexInteraction<LocalDynamicsName<FirstInteraction, OtherInteractions...>, CommonParameters...>
    : public LocalDynamicsName<FirstInteraction, CommonParameters...>
{
  protected:
    ComplexInteraction<LocalDynamicsName<OtherInteractions...>, CommonParameters...> other_interactions_;

  public:
    template <class FirstParameterSet, typename... OtherParameterSets>
    explicit ComplexInteraction(FirstParameterSet &&first_parameter_set,
                                OtherParameterSets &&...other_parameter_sets)
        : LocalDynamicsName<FirstInteraction, CommonParameters...>(first_parameter_set),
          other_interactions_(std::forward<OtherParameterSets>(other_parameter_sets)...),
          device_kernel(LocalDynamicsName<FirstInteraction, CommonParameters...>::device_kernel.get_ptr(),
                        other_interactions_.device_kernel.get_ptr()) {};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        LocalDynamicsName<FirstInteraction, CommonParameters...>::interaction(index_i, dt);
        other_interactions_.interaction(index_i, dt);
    };

    execution::DeviceImplementation<
        ComplexInteractionKernel<LocalDynamicsName<CommonParameters...>,FirstInteraction, OtherInteractions...>> device_kernel;
};
} // namespace SPH
#endif // BASE_LOCAL_DYNAMICS_H
