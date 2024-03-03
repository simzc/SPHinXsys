#ifndef RELAX_DYNAMICS_HPP
#define RELAX_DYNAMICS_HPP

#include "relax_dynamics.h"

namespace SPH
{
namespace relax_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
BaseIntegration<DataDelegationType>::BaseIntegration(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      sph_adaptation_(sph_body_.sph_adaptation_),
      residue_(*this->particles_->template registerSharedVariable<Vecd>("ZeroOrderResidue")) {}
//=================================================================================================//
template <class KernelCorrectionType>
RelaxationIntegration<Inner<>, KernelCorrectionType>::
RelaxationIntegration(BaseInnerRelation &inner_relation, std::string sub_shape_name)
    :BaseIntegration<RelaxDataDelegateInner>(inner_relation),
    relax_shape_(*DynamicCast<ComplexShape>(this, *sph_body_.initial_shape_)
        .getSubShapeByName(sub_shape_name), correction_(particles_) {}
//=================================================================================================//
template <class KernelCorrectionType>
void RelaxationIntegration<Inner<>, KernelCorrectionType>::interaction(size_t index_i, Real dt)
{
    Vecd residue = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        residue -= 2.0 * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
    }
    residue_[index_i] = residue;
}
//=================================================================================================//
template <class KernelCorrectionType>
void RelaxationIntegration<Inner<>, KernelCorrectionType>::update(size_t index_i, Real dt)
{

}
//=================================================================================================//
template <class KernelCorrectionType>
template <typename... Args>
RelaxationIntegration<Inner<LevelSetCorrtion>>::RelaxationIntegration(Args &&...args)
:RelaxationIntegration<Inner<>>(std::forward<Args>(args)...), pos_(particles_->pos)
level_set_shape_(DynamicCast<LevelSetShape>(this, this->getRelaxShape())){};
//=================================================================================================//
template <class KernelCorrectionType>
RelaxationIntegration<Inner<LevelSetCorrection>>::
RelaxationIntegration(ConstructorArgs<BodyRelationType, FirstArg> parameters)
: RelaxationResidue(parameter.body_relation_, std::get<0>(parameters.others_)){};
//=================================================================================================//
template <class KernelCorrectionType>
void RelaxationIntegration<Inner<LevelSetCorrection>, KernelCorrectionType>::
interaction(size_t index_i, Real dt)
{
    RelaxationResidue<Inner<>>::interaction(index_i, dt);
    residue_[index_i] -= 2.0 * level_set_shape_.computeKernelGradientIntegral(
                                   pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
};
//=================================================================================================//
template <class KernelCorrectionType>
RelaxationIntegration<Contact<>, KernelCorrectionType>::
RelaxationIntegration(BaseContactRelation &contact_relation)
:BaseIntegration<RelaxDataDelegateContact>(contact_relation), correction_(this->particles)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_corrections_.push_back(KernelCorrectionType(this->contact_particles_[k]));
    }
};
//=================================================================================================//
template <class KernelCorrectionType>
void RelaxationIntegration<Contact<>, KernelCorrectionType>::interaction(size_t index_i, Real dt)
{
    Vecd residue = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            residue -= 2.0 * contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
        }
    }
    residue_[index_i] += residue;
};
//=================================================================================================//
} //namespace relax_dynamics
} //namespace SPH
#endif //RELAX_DYNAMICS_HPP