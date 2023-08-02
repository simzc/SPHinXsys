#include "bound_diagnosis.h"

namespace SPH
{
//=================================================================================================//
VelocityBoundCheck::
    VelocityBoundCheck(SPHBody &sph_body, Real velocity_bound)
    : LocalDynamicsReduce<bool, ReduceOR>(sph_body, false),
      GeneralDataDelegateSimple(sph_body),
      vel_(particles_->vel_), velocity_bound_(velocity_bound) {}
//=================================================================================================//
bool VelocityBoundCheck::reduce(size_t index_i, Real dt)
{
    return vel_[index_i].norm() > velocity_bound_;
}
//=================================================================================================//
} // namespace SPH
