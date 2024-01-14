#ifndef FLUID_STRUCTURE_INTERACTION_HPP
#define FLUID_STRUCTURE_INTERACTION_HPP

#include "fluid_structure_interaction.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType>
BasePressureForceFromFluid<RiemannSolverType>::
    BasePressureForceFromFluid(BaseContactRelation &contact_relation)
    : BasePressureForceFromFluid(true, contact_relation)
{
    particles_->registerVariable(force_from_fluid_, "PressureForceFromFluid");
}
//=================================================================================================//
template <class RiemannSolverType>
BasePressureForceFromFluid<RiemannSolverType>::
    BasePressureForceFromFluid(bool mostDerived, BaseContactRelation &contact_relation)
    : ForcePrior(contact_relation.getSPHBody(), "ForceFromFluid"),
      BaseForceFromFluid(contact_relation),
      mass_(particles_->mass_), vel_ave_(*particles_->AverageVelocity()),
      force_ave_(*particles_->AverageForce()), n_(particles_->n_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_rho_n_.push_back(&(contact_particles_[k]->rho_));
        contact_mass_.push_back(&(contact_particles_[k]->mass_));
        contact_vel_.push_back(&(contact_particles_[k]->vel_));
        contact_p_.push_back(contact_particles_[k]->template getVariableByName<Real>("Pressure"));
        contact_force_prior_.push_back(&(contact_particles_[k]->force_prior_));
        riemann_solvers_.push_back(RiemannSolverType(*contact_fluids_[k], *contact_fluids_[k]));
    }
}
//=================================================================================================//
template <class RiemannSolverType>
void BasePressureForceFromFluid<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Real> &rho_n_k = *(contact_rho_n_[k]);
        StdLargeVec<Real> &mass_k = *(contact_mass_[k]);
        StdLargeVec<Real> &p_k = *(contact_p_[k]);
        StdLargeVec<Vecd> &vel_k = *(contact_vel_[k]);
        StdLargeVec<Vecd> &force_prior_k = *(contact_force_prior_[k]);
        RiemannSolverType &riemann_solvers_k = riemann_solvers_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Vecd wall_acceleration = force_prior_k[index_j] / mass_k[index_j] -
                                     force_ave_[index_i] / particles_->mass_[index_i];
            Real p_in_wall = p_k[index_j] +
                             rho_n_k[index_j] * r_ij * SMAX(Real(0), wall_acceleration.dot(e_ij));
            Real u_jump = 2.0 * (vel_k[index_j] - vel_ave_[index_i]).dot(n_[index_i]);
            force += (riemann_solvers_k.DissipativePJump(u_jump) * n_[index_i] - (p_in_wall + p_k[index_j]) * e_ij) *
                     Vol_[index_i] * contact_neighborhood.dW_ijV_j_[n];
        }
    }
    force_from_fluid_[index_i] = force;
    force_[index_i] = force;
}
//=================================================================================================//
template <class PressureForceType>
template <class ViscousForceFromFluidType>
BaseAllForceFromFluid<PressureForceType>::
    BaseAllForceFromFluid(BaseContactRelation &contact_relation,
                          ViscousForceFromFluidType &viscous_force_from_fluid)
    : PressureForceType(false, contact_relation),
      viscous_force_from_fluid_(viscous_force_from_fluid.getForceFromFluid())
{
    this->particles_->registerVariable(this->force_from_fluid_, "AllForceFromFluid");
}
//=================================================================================================//
template <class PressureForceType>
void BaseAllForceFromFluid<PressureForceType>::interaction(size_t index_i, Real dt)
{
    PressureForceType::interaction(index_i, dt);
    this->force_from_fluid_[index_i] += viscous_force_from_fluid_[index_i];
    this->force_[index_i] += viscous_force_from_fluid_[index_i];
}
//=================================================================================================//
template <class ForceFromFluidDynamicsType>
TotalForceFromFluid::TotalForceFromFluid(ForceFromFluidDynamicsType &force_from_fluid_dynamics,
                                         const std::string &force_name)
    : LocalDynamicsReduce<Vecd, ReduceSum<Vecd>>(force_from_fluid_dynamics.getSPHBody(), Vecd::Zero()),
      force_from_fluid_dynamics_(force_from_fluid_dynamics),
      force_from_fluid_(force_from_fluid_dynamics.getForceFromFluid())
{
    quantity_name_ = force_name;
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
#endif // FLUID_STRUCTURE_INTERACTION_HPP
