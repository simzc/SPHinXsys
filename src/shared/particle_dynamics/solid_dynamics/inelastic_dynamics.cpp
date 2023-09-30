#include "inelastic_dynamics.h"

namespace SPH
{
//=====================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
PlasticIntegration1stHalf::
    PlasticIntegration1stHalf(BaseInnerRelation &inner_relation)
    : Integration1stHalf(inner_relation),
      plastic_solid_(DynamicCast<PlasticSolid>(this, elastic_solid_))
{
    numerical_dissipation_factor_ = 0.5;
}
//=================================================================================================//
void PlasticIntegration1stHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    rho_[index_i] = rho0_ / F_[index_i].determinant();

    // TODO: this naming is temporary, to be revised.
    stress_PK1_B_[index_i] = plastic_solid_.PlasticConstitutiveRelation(F_[index_i], index_i, dt);
}
//=================================================================================================//
DecomposedPlasticIntegration1stHalf::
    DecomposedPlasticIntegration1stHalf(BaseInnerRelation &inner_relation)
    : DecomposedIntegration1stHalf(inner_relation),
      plastic_solid_(DynamicCast<PlasticSolid>(this, elastic_solid_))
{
    particles_->registerVariable(scaling_matrix_, "ScalingMatrix");
}
//=================================================================================================//
void DecomposedPlasticIntegration1stHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    Real J = F_[index_i].determinant();
    Real one_over_J = 1.0 / J;
    rho_[index_i] = rho0_ * one_over_J;

    Matd be = plastic_solid_.ElasticLeftCauchy(F_[index_i], index_i, dt);
    Matd bt = F_[index_i] * F_[index_i].transpose();
    scaling_matrix_[index_i] = be * bt.inverse();

    Matd inverse_F_T = F_[index_i].inverse().transpose();
    Matd deviatoric_PK = plastic_solid_.DeviatoricKirchhoff(be - be.trace() * OneOverDimensions * Matd::Identity());
    stress_on_particle_[index_i] =
        inverse_F_T * elastic_solid_.VolumetricKirchhoff(J) +
        (deviatoric_PK + elastic_solid_.NumericalDampingLeftCauchy(F_[index_i], dF_dt_[index_i], smoothing_length_, index_i)) * inverse_F_T;
}
//=================================================================================================//
void DecomposedPlasticIntegration1stHalf::interaction(size_t index_i, Real dt)
{
    // including gravity and force from fluid
    Vecd acceleration = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Vecd e_ij = inner_neighborhood.e_ij_[n];
        Vecd pair_distance = pos_[index_i] - pos_[index_j];
        Vecd shear_force_ij = 100.0 * elastic_solid_.ShearModulus() *
                              (0.5 * (scaling_matrix_[index_i] + scaling_matrix_[index_j]) * e_ij).dot(e_ij) *
                              (pair_distance - 0.5 * (F_[index_i] + F_[index_j]) * r_ij * e_ij);
        acceleration += ((stress_on_particle_[index_i] + stress_on_particle_[index_j]) * e_ij + shear_force_ij) *
                        inner_neighborhood.dW_ijV_j_[n] * inv_rho0_;
    }
    acc_[index_i] = acceleration;
}
//=================================================================================================//
} // namespace solid_dynamics
  //=====================================================================================================//
} // namespace SPH
