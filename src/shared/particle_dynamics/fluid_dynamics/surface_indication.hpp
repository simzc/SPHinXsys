#pragma once

#include "surface_indication.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
FreeSurfaceIndication<DataDelegationType>::FreeSurfaceIndication(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      indicator_(*this->particles_->template getVariableByName<int>("Indicator")),
      threshold_by_dimensions_(0.75 * Dimensions)
{
    this->particles_->registerVariable(pos_div_, "PositionDivergence");
}
//=================================================================================================//
template <class FreeSurfaceIdentificationType>
template <typename... ConstructorArgs>
SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIdentificationType>::
    SpatialTemporalFreeSurfaceIdentification(ConstructorArgs &&...args)
    : FreeSurfaceIdentificationType(std::forward<ConstructorArgs>(args)...)
{
    this->particles_->registerVariable(previous_surface_indicator_, "PreviousSurfaceIndicator", 1);
    this->particles_->template registerSortableVariable<int>("PreviousSurfaceIndicator");
}
//=================================================================================================//
template <class FreeSurfaceIdentificationType>
void SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIdentificationType>::
    interaction(size_t index_i, Real dt)
{
    FreeSurfaceIdentificationType::interaction(index_i, dt);

    if (this->pos_div_[index_i] < this->threshold_by_dimensions_ &&
        previous_surface_indicator_[index_i] != 1 &&
        !isNearPreviousFreeSurface(index_i))
        this->pos_div_[index_i] = 2.0 * this->threshold_by_dimensions_;
}
//=================================================================================================//
template <class FreeSurfaceIdentificationType>
bool SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIdentificationType>::
    isNearPreviousFreeSurface(size_t index_i)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        if (previous_surface_indicator_[inner_neighborhood.j_[n]] == 1)
        {
            is_near_surface = true;
            break;
        }
    }
    return is_near_surface;
}
//=================================================================================================//
template <class FreeSurfaceIdentificationType>
void SpatialTemporalFreeSurfaceIdentification<FreeSurfaceIdentificationType>::
    update(size_t index_i, Real dt)
{
    FreeSurfaceIdentificationType::update(index_i, dt);

    previous_surface_indicator_[index_i] = this->indicator_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
  //=================================================================================================//
} // namespace SPH