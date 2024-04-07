#include "solid_particles.h"
#include "solid_particles_variable.h"

#include "base_body.h"
#include "elastic_solid.h"
#include "inelastic_solid.h"
#include "xml_engine.h"

namespace SPH
{
//=============================================================================================//
SolidParticles::SolidParticles(SPHBody &sph_body, Solid *solid)
    : BaseParticles(sph_body, solid), solid_(*solid) {}
//=================================================================================================//
void SolidParticles::initializeOtherVariables()
{
    BaseParticles::initializeOtherVariables();
}
//=============================================================================================//
ElasticSolidParticles::
    ElasticSolidParticles(SPHBody &sph_body, ElasticSolid *elastic_solid)
    : SolidParticles(sph_body, elastic_solid),
      elastic_solid_(*elastic_solid) {}
//=================================================================================================//
void ElasticSolidParticles::initializeOtherVariables()
{
    SolidParticles::initializeOtherVariables();

    registerSharedVariable<Matd>("DeformationGradient", Matd::Identity());

    addDerivedVariableToWrite<Displacement>();
    addDerivedVariableToWrite<VonMisesStrain>();
    addDerivedVariableToWrite<VonMisesStress>();

    addVariableToRestart<Matd>("DeformationGradient");
}
//=============================================================================================//
ShellParticles::ShellParticles(SPHBody &sph_body, ElasticSolid *elastic_solid)
    : ElasticSolidParticles(sph_body, elastic_solid), thickness_ref_(1.0)
{
    //----------------------------------------------------------------------
    //		modify kernel function for surface particles
    //----------------------------------------------------------------------
    sph_body.sph_adaptation_->getKernel()->reduceOnce();
    //----------------------------------------------------------------------
    //		register geometric data only
    //----------------------------------------------------------------------
    registerSharedVariable<Vecd>("NormalDirection");
    registerSharedVariable<Real>("Thickness");

    addVariableToWrite<Vecd>("NormalDirection");

    addVariableToList<Vecd>(variables_to_reload_, "NormalDirection");
    addVariableToList<Real>(variables_to_reload_, "Thickness");
}
//=================================================================================================//
void ShellParticles::initializeOtherVariables()
{
    BaseParticles::initializeOtherVariables();
    /**
     * register particle data
     */
    registerVariable(pos0_, "InitialPosition",
                     [&](size_t i) -> Vecd
                     { return pos_[i]; });
    registerVariable(n0_, "InitialNormalDirection",
                     [&](size_t i) -> Vecd
                     { return n_[i]; });
    registerVariable(transformation_matrix_, "TransformationMatrix");
    registerVariable(B_, "LinearGradientCorrectionMatrix", [&](size_t i) -> Matd
                     { return Matd::Identity(); });
    registerVariable(F_, "DeformationGradient", [&](size_t i) -> Matd
                     { return Matd::Identity(); });
    registerVariable(dF_dt_, "DeformationRate");
    registerVariable(pseudo_n_, "PseudoNormal",
                     [&](size_t i) -> Vecd
                     { return n_[i]; });
    registerVariable(dpseudo_n_dt_, "PseudoNormalChangeRate");
    registerVariable(dpseudo_n_d2t_, "PseudoNormal2ndOrderTimeDerivative");
    registerVariable(rotation_, "Rotation");
    registerVariable(angular_vel_, "AngularVelocity");
    registerVariable(dangular_vel_dt_, "AngularAcceleration");
    registerVariable(F_bending_, "BendingDeformationGradient");
    registerVariable(dF_bending_dt_, "BendingDeformationGradientChangeRate");
    registerVariable(global_shear_stress_, "GlobalShearStress");
    registerVariable(global_stress_, "GlobalStress");
    registerVariable(global_moment_, "GlobalMoment");
    registerVariable(mid_surface_cauchy_stress_, "MidSurfaceCauchyStress");
    registerVariable(numerical_damping_scaling_, "NumericalDampingScaling",
                     [&](size_t i) -> Matd
                     { return Matd::Identity() * sph_body_.sph_adaptation_->ReferenceSmoothingLength(); });
    /**
     * for FSI
     */
    registerVariable(vel_ave_, "AverageVelocity");
    registerVariable(force_ave_, "AverageForce");
    /**
     * for rotation.
     */
    addVariableToRestart<Matd>("DeformationGradient");
    addVariableToRestart<Vecd>("PseudoNormal");
    addVariableToRestart<Vecd>("Rotation");
    addVariableToRestart<Vecd>("AngularVelocity");
    /**
     * add basic output particle data
     */
    addVariableToWrite<Vecd>("NormalDirection");
    addDerivedVariableToWrite<Displacement>();
    addDerivedVariableToWrite<VonMisesStress>();
    addDerivedVariableToWrite<VonMisesStrain>();
    addVariableToRestart<Matd>("DeformationGradient");
    addVariableToWrite<Vecd>("Rotation");
    addDerivedVariableToWrite<MidSurfaceVonMisesStress>();
    /**
     * initialize transformation matrix
     */
    for (size_t i = 0; i != real_particles_bound_; ++i)
    {
        transformation_matrix_[i] = getTransformationMatrix(n_[i]);
        numerical_damping_scaling_[i](Dimensions - 1, Dimensions - 1) =
            thickness_[i] < sph_body_.sph_adaptation_->ReferenceSmoothingLength() ? thickness_[i] : sph_body_.sph_adaptation_->ReferenceSmoothingLength();
    }
}

} // namespace SPH
