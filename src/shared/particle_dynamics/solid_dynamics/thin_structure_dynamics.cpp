#include "thin_structure_dynamics.h"

namespace SPH
{
//=====================================================================================================//
namespace thin_structure_dynamics
{
//=================================================================================================//
ShellDynamicsInitialCondition::ShellDynamicsInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body), ShellDataSimple(sph_body),
      n0_(particles_->n0_), n_(particles_->n_), pseudo_n_(particles_->pseudo_n_),
      pos0_(particles_->pos0_), transformation_matrix_(particles_->transformation_matrix_) {}
//=================================================================================================//
UpdateShellNormalDirection::UpdateShellNormalDirection(SPHBody &sph_body)
    : LocalDynamics(sph_body), ShellDataSimple(sph_body),
      n_(particles_->n_), F_(particles_->F_),
      transformation_matrix_(particles_->transformation_matrix_) {}
//=========================================================================================//
void UpdateShellNormalDirection::update(size_t index_i, Real dt)
{
    /** Calculate the current normal direction of mid-surface. */
    n_[index_i] = transformation_matrix_[index_i].transpose() * getNormalFromDeformationGradientTensor(F_[index_i]);
}
//=================================================================================================//
ShellAcousticTimeStepSize::ShellAcousticTimeStepSize(SPHBody &sph_body, Real CFL)
    : LocalDynamicsReduce<Real, ReduceMin>(sph_body, MaxReal),
      ShellDataSimple(sph_body), CFL_(CFL), vel_(particles_->vel_), force_(particles_->force_),
      angular_vel_(particles_->angular_vel_), dangular_vel_dt_(particles_->dangular_vel_dt_),
      force_prior_(particles_->force_prior_),
      thickness_(particles_->thickness_), mass_(particles_->mass_),
      rho0_(particles_->elastic_solid_.ReferenceDensity()),
      E0_(particles_->elastic_solid_.YoungsModulus()),
      nu_(particles_->elastic_solid_.PoissonRatio()),
      c0_(particles_->elastic_solid_.ReferenceSoundSpeed()),
      smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
Real ShellAcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    // Since the particle does not change its configuration in pressure relaxation step,
    // I chose a time-step size according to Eulerian method.
    Real time_setp_0 = SMIN((Real)sqrt(smoothing_length_ / ((force_[index_i] + force_prior_[index_i]).norm() / mass_[index_i] + TinyReal)),
                            smoothing_length_ / (c0_ + vel_[index_i].norm()));
    Real time_setp_1 = SMIN((Real)sqrt(1.0 / (dangular_vel_dt_[index_i].norm() + TinyReal)),
                            Real(1.0) / (angular_vel_[index_i].norm() + TinyReal));
    Real time_setp_2 = smoothing_length_ * (Real)sqrt(rho0_ * (1.0 - nu_ * nu_) / E0_ /
                                                      (2.0 + (Pi * Pi / 12.0) * (1.0 - nu_) *
                                                                 (1.0 + 1.5 * pow(smoothing_length_ / thickness_[index_i], 2))));
    return CFL_ * SMIN(time_setp_0, time_setp_1, time_setp_2);
}
//=================================================================================================//
ShellCorrectConfiguration::
    ShellCorrectConfiguration(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ShellDataInner(inner_relation),
      B_(particles_->B_),
      n0_(particles_->n0_), transformation_matrix_(particles_->transformation_matrix_) {}
//=================================================================================================//
ShellDeformationGradientTensor::
    ShellDeformationGradientTensor(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ShellDataInner(inner_relation),
      pos_(particles_->pos_), pseudo_n_(particles_->pseudo_n_), n0_(particles_->n0_),
      B_(particles_->B_), F_(particles_->F_), F_bending_(particles_->F_bending_),
      transformation_matrix_(particles_->transformation_matrix_) {}
//=================================================================================================//
BaseShellRelaxation::BaseShellRelaxation(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ShellDataInner(inner_relation),
      rho_(particles_->rho_),
      thickness_(particles_->thickness_), mass_(particles_->mass_),
      pos_(particles_->pos_), vel_(particles_->vel_),
      force_(particles_->force_),
      force_prior_(particles_->force_prior_),
      n0_(particles_->n0_), pseudo_n_(particles_->pseudo_n_),
      dpseudo_n_dt_(particles_->dpseudo_n_dt_), dpseudo_n_d2t_(particles_->dpseudo_n_d2t_),
      rotation_(particles_->rotation_), angular_vel_(particles_->angular_vel_),
      dangular_vel_dt_(particles_->dangular_vel_dt_),
      B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_),
      F_bending_(particles_->F_bending_), dF_bending_dt_(particles_->dF_bending_dt_),
      transformation_matrix_(particles_->transformation_matrix_) {}
//=================================================================================================//
ShellStressRelaxationFirstHalf::
    ShellStressRelaxationFirstHalf(BaseInnerRelation &inner_relation,
                                   int number_of_gaussian_points, bool hourglass_control, Real hourglass_control_factor)
    : BaseShellRelaxation(inner_relation),
      elastic_solid_(particles_->elastic_solid_),
      global_stress_(particles_->global_stress_),
      global_moment_(particles_->global_moment_),
      mid_surface_cauchy_stress_(particles_->mid_surface_cauchy_stress_),
      numerical_damping_scaling_(particles_->numerical_damping_scaling_),
      global_shear_stress_(particles_->global_shear_stress_),
      rho0_(elastic_solid_.ReferenceDensity()),
      inv_rho0_(1.0 / rho0_),
      smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      E0_(elastic_solid_.YoungsModulus()),
      G0_(elastic_solid_.ShearModulus()),
      nu_(elastic_solid_.PoissonRatio()),
      hourglass_control_factor_(hourglass_control_factor),
      hourglass_control_(hourglass_control),
      number_of_gaussian_points_(number_of_gaussian_points)
{
    /** Note that, only three-point and five-point Gaussian quadrature rules are defined. */
    switch (number_of_gaussian_points)
    {
    case 1:
        gaussian_point_ = one_gaussian_point_;
        gaussian_weight_ = one_gaussian_weight_;
        break;
    case 5:
        gaussian_point_ = five_gaussian_points_;
        gaussian_weight_ = five_gaussian_weights_;
        break;
    default:
        gaussian_point_ = three_gaussian_points_;
        gaussian_weight_ = three_gaussian_weights_;
    }
}
//=================================================================================================//
void ShellStressRelaxationFirstHalf::initialization(size_t index_i, Real dt)
{
    // Note that F_[index_i], F_bending_[index_i], dF_dt_[index_i], dF_bending_dt_[index_i]
    // and rotation_[index_i], angular_vel_[index_i], dangular_vel_dt_[index_i], B_[index_i]
    // are defined in local coordinates, while others in global coordinates.
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
    pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;

    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;

    Real J = F_[index_i].determinant();
    Matd inverse_F = F_[index_i].inverse();

    rho_[index_i] = rho0_ / J;

    /** Get transformation matrix from global coordinates to current local coordinates. */
    Matd current_transformation_matrix = getTransformationMatrix(pseudo_n_[index_i]);

    Matd resultant_stress = Matd::Zero();
    Matd resultant_moment = Matd::Zero();
    Vecd resultant_shear_stress = Vecd::Zero();

    for (int i = 0; i != number_of_gaussian_points_; ++i)
    {
        Matd F_gaussian_point = F_[index_i] + gaussian_point_[i] * F_bending_[index_i] * thickness_[index_i] * 0.5;
        Matd dF_gaussian_point_dt = dF_dt_[index_i] + gaussian_point_[i] * dF_bending_dt_[index_i] * thickness_[index_i] * 0.5;
        Matd inverse_F_gaussian_point = F_gaussian_point.inverse();
        Matd current_local_almansi_strain = current_transformation_matrix * transformation_matrix_[index_i].transpose() * 0.5 *
                                            (Matd::Identity() - inverse_F_gaussian_point.transpose() * inverse_F_gaussian_point) *
                                            transformation_matrix_[index_i] * current_transformation_matrix.transpose();

        /** correct Almansi strain tensor according to plane stress problem. */
        current_local_almansi_strain = getCorrectedAlmansiStrain(current_local_almansi_strain, nu_);

        /** correct out-plane numerical damping. */
        Matd cauchy_stress = elastic_solid_.StressCauchy(current_local_almansi_strain, index_i) + current_transformation_matrix * transformation_matrix_[index_i].transpose() * F_gaussian_point *
                                                                                                      elastic_solid_.NumericalDampingRightCauchy(F_gaussian_point, dF_gaussian_point_dt, numerical_damping_scaling_[index_i], index_i) * F_gaussian_point.transpose() * transformation_matrix_[index_i] * current_transformation_matrix.transpose() / F_gaussian_point.determinant();

        /** Impose modeling assumptions. */
        cauchy_stress.col(Dimensions - 1) *= shear_correction_factor_;
        cauchy_stress.row(Dimensions - 1) *= shear_correction_factor_;
        cauchy_stress(Dimensions - 1, Dimensions - 1) = 0.0;

        if (i == 0)
        {
            mid_surface_cauchy_stress_[index_i] = cauchy_stress;
        }

        /** Integrate Cauchy stress along thickness. */
        resultant_stress +=
            0.5 * thickness_[index_i] * gaussian_weight_[i] * cauchy_stress;
        resultant_moment +=
            0.5 * thickness_[index_i] * gaussian_weight_[i] * (cauchy_stress * gaussian_point_[i] * thickness_[index_i] * 0.5);
        resultant_shear_stress -=
            0.5 * thickness_[index_i] * gaussian_weight_[i] * cauchy_stress.col(Dimensions - 1);

        resultant_stress.col(Dimensions - 1) = Vecd::Zero();
        resultant_moment.col(Dimensions - 1) = Vecd::Zero();
    }

    /** stress and moment in global coordinates for pair interaction */
    global_stress_[index_i] = J * current_transformation_matrix.transpose() * resultant_stress * current_transformation_matrix * transformation_matrix_[index_i].transpose() * inverse_F.transpose() * transformation_matrix_[index_i];
    global_moment_[index_i] = J * current_transformation_matrix.transpose() * resultant_moment * current_transformation_matrix * transformation_matrix_[index_i].transpose() * inverse_F.transpose() * transformation_matrix_[index_i];
    global_shear_stress_[index_i] = J * current_transformation_matrix.transpose() * resultant_shear_stress;
}
//=================================================================================================//
void ShellStressRelaxationFirstHalf::update(size_t index_i, Real dt)
{
    vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
    angular_vel_[index_i] += dangular_vel_dt_[index_i] * dt;
}
//=================================================================================================//
void ShellStressRelaxationSecondHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
    dpseudo_n_dt_[index_i] = transformation_matrix_[index_i].transpose() *
                             getVectorChangeRateAfterThinStructureRotation(local_pseudo_n_0, rotation_[index_i], angular_vel_[index_i]);
    pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void ShellStressRelaxationSecondHalf::update(size_t index_i, Real dt)
{
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
ConstrainShellBodyRegion::
    ConstrainShellBodyRegion(BodyPartByParticle &body_part)
    : BaseLocalDynamics<BodyPartByParticle>(body_part), ShellDataSimple(sph_body_),
      vel_(particles_->vel_), angular_vel_(particles_->angular_vel_) {}
//=================================================================================================//
void ConstrainShellBodyRegion::update(size_t index_i, Real dt)
{
    vel_[index_i] = Vecd::Zero();
    angular_vel_[index_i] = Vecd::Zero();
}
//=================================================================================================//
ConstrainShellBodyRegionAlongAxis::ConstrainShellBodyRegionAlongAxis(BodyPartByParticle &body_part, int axis)
    : BaseLocalDynamics<BodyPartByParticle>(body_part), ShellDataSimple(sph_body_),
      axis_(axis), pos_(particles_->pos_), pos0_(particles_->pos0_), vel_(particles_->vel_),
      force_(particles_->force_), rotation_(particles_->rotation_), angular_vel_(particles_->angular_vel_),
      dangular_vel_dt_(particles_->dangular_vel_dt_), mass_(particles_->mass_) {}
//=================================================================================================//
void ConstrainShellBodyRegionAlongAxis::update(size_t index_i, Real dt)
{
    vel_[index_i][axis_] = 0.0;
    vel_[index_i][2] = 0.0;
    force_[index_i][axis_] = 0.0;
    force_[index_i][2] = 0.0;

    angular_vel_[index_i][1 - axis_] = 0.0;
    dangular_vel_dt_[index_i][1 - axis_] = 0.0;
}
//=================================================================================================//
DistributingPointForcesToShell::
    DistributingPointForcesToShell(SPHBody &sph_body, std::vector<Vecd> point_forces,
                                   std::vector<Vecd> reference_positions, Real time_to_full_external_force,
                                   Real particle_spacing_ref, Real h_spacing_ratio)
    : LocalDynamics(sph_body), ShellDataSimple(sph_body),
      point_forces_(point_forces), reference_positions_(reference_positions),
      time_to_full_external_force_(time_to_full_external_force),
      particle_spacing_ref_(particle_spacing_ref), h_spacing_ratio_(h_spacing_ratio),
      pos0_(particles_->pos0_), force_prior_(particles_->force_prior_),
      thickness_(particles_->thickness_)
{
    for (size_t i = 0; i < point_forces_.size(); i++)
    {
        weight_.push_back(StdLargeVec<Real>(0.0));
        time_dependent_point_forces_.push_back(Vecd::Zero());
        sum_of_weight_.push_back(0.0);
        particles_->registerVariable(weight_[i], "Weight_" + std::to_string(i));
    }

    getWeight(); // TODO: should be revised and parallelized, using SimpleDynamics
}
//=================================================================================================//
void DistributingPointForcesToShell::getWeight()
{
    Kernel *kernel_ = sph_body_.sph_adaptation_->getKernel();
    Real reference_smoothing_length = sph_body_.sph_adaptation_->ReferenceSmoothingLength();
    Real smoothing_length = h_spacing_ratio_ * particle_spacing_ref_;
    Real h_ratio = reference_smoothing_length / smoothing_length;
    Real cutoff_radius_sqr = pow(2.0 * smoothing_length, 2);
    for (size_t i = 0; i < point_forces_.size(); ++i)
    {
        sum_of_weight_[i] = 0.0;
        for (size_t index = 0; index < particles_->total_real_particles_; ++index)
        {
            weight_[i][index] = 0.0;
            Vecd displacement = reference_positions_[i] - pos0_[index];
            if (displacement.squaredNorm() <= cutoff_radius_sqr)
            {
                weight_[i][index] = kernel_->W(h_ratio, displacement.norm(), displacement);
                sum_of_weight_[i] += weight_[i][index];
            }
        }
    }
}
//=================================================================================================//
void DistributingPointForcesToShell::setupDynamics(Real dt)
{
    Real current_time = GlobalStaticVariables::physical_time_;
    for (size_t i = 0; i < point_forces_.size(); ++i)
    {
        time_dependent_point_forces_[i] = current_time < time_to_full_external_force_
                                              ? current_time * point_forces_[i] / time_to_full_external_force_
                                              : point_forces_[i];
    }
}
//=================================================================================================//
void DistributingPointForcesToShell::update(size_t index_i, Real dt)
{
    force_prior_[index_i] = Vecd::Zero();
    for (size_t i = 0; i < point_forces_.size(); ++i)
    {
        Vecd force = weight_[i][index_i] / (sum_of_weight_[i] + TinyReal) * time_dependent_point_forces_[i];
        force_prior_[index_i] += force;
    }
}
//=================================================================================================//
ShellCurvature::ShellCurvature(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ShellDataInner(inner_relation),
      n0_(particles_->n0_), B_(particles_->B_), transformation_matrix_(particles_->transformation_matrix_),
      n_(particles_->n_), F_(particles_->F_), F_bending_(particles_->F_bending_),
      k1_(*particles_->registerSharedVariable<Real>("1stPrincipleCurvature")),
      k2_(*particles_->registerSharedVariable<Real>("2ndPrincipleCurvature"))
{
    particles_->registerVariable(dn_0_, "InitialNormalGradient");
};
//=================================================================================================//
void ShellCurvature::compute_initial_curvature()
{
    particle_for(
        par,
        particles_->total_real_particles_,
        [this](size_t index_i)
        {
            Matd dn_0_i = Matd::Zero();
            // transform initial local B_ to global B_
            const Matd B_global_i = transformation_matrix_[index_i].transpose() * B_[index_i] * transformation_matrix_[index_i];
            const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                const size_t index_j = inner_neighborhood.j_[n];
                const Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                dn_0_i -= (n0_[index_i] - n0_[index_j]) * gradW_ijV_j.transpose();
            }
            dn_0_[index_i] = dn_0_i * B_global_i;
            auto [k1, k2] = get_principle_curvatures(dn_0_[index_i]);
            k1_[index_i] = k1;
            k2_[index_i] = k2;
        });
}
//=================================================================================================//
void ShellCurvature::update(size_t index_i, Real)
{
    Matd dn_0_i = dn_0_[index_i] + transformation_matrix_[index_i].transpose() *
                                       F_bending_[index_i] * transformation_matrix_[index_i];
    Matd inverse_F = F_[index_i].inverse();
    Matd dn_i = dn_0_i * transformation_matrix_[index_i].transpose() * inverse_F * transformation_matrix_[index_i];
    auto [k1, k2] = get_principle_curvatures(dn_i);
    k1_[index_i] = k1;
    k2_[index_i] = k2;
}
//=================================================================================================//
AverageShellCurvature::AverageShellCurvature(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ShellDataInner(inner_relation),
      n_(particles_->n_),
      k1_ave_(*particles_->registerSharedVariable<Real>("Average1stPrincipleCurvature")),
      k2_ave_(*particles_->registerSharedVariable<Real>("Average2ndPrincipleCurvature")){};
//=================================================================================================//
void AverageShellCurvature::update(size_t index_i, Real)
{
    Matd dn_i = Matd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        const size_t index_j = inner_neighborhood.j_[n];
        const Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        dn_i -= (n_[index_i] - n_[index_j]) * gradW_ijV_j.transpose();
    }
    auto [k1, k2] = get_principle_curvatures(dn_i);
    k1_ave_[index_i] = k1;
    k2_ave_[index_i] = k2;
}
//=================================================================================================//
} // namespace thin_structure_dynamics
} // namespace SPH