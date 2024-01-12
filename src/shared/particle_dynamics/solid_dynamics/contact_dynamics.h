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
 * @file 	contact_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid contact dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef CONTACT_DYNAMICS_H
#define CONTACT_DYNAMICS_H

#include "force_prior.h"
#include "general_solid_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{
typedef DataDelegateContact<SolidParticles, SolidParticles> ContactDynamicsData;
typedef DataDelegateContact<SolidParticles, SolidParticles> ContactWithWallData;

class ContactDensityAccessor
{
  protected:
    ContactDensityAccessor(BaseParticles &particles, const std::string &variable_name)
        : contact_density_(*particles.registerSharedVariable<Real>(variable_name)){};
    ~ContactDensityAccessor() = default;
    StdLargeVec<Real> &contact_density_;
};

/**
 * @class SelfContactDensitySummation
 * @brief Computing the summation density due to solid self-contact model.
 */
class SelfContactDensitySummation : public ContactDensityAccessor, public LocalDynamics, public SolidDataInner
{
  public:
    explicit SelfContactDensitySummation(SelfSurfaceContactRelation &self_contact_relation);
    virtual ~SelfContactDensitySummation(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Real sigma = 0.0;
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            Real corrected_W_ij = std::max(inner_neighborhood.W_ij_[n] - offset_W_ij_, Real(0));
            sigma += corrected_W_ij * mass_[inner_neighborhood.j_[n]];
        }
        contact_density_[index_i] = sigma;
    };

  protected:
    StdLargeVec<Real> &mass_;
    Real offset_W_ij_;
};

/**
 * @class ContactDensitySummation
 * @brief Computing the summation density due to solid-solid contact model.
 */
class ContactDensitySummation : public ContactDensityAccessor, public LocalDynamics, public ContactDynamicsData
{
  public:
    explicit ContactDensitySummation(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~ContactDensitySummation(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        /** Contact interaction. */
        Real sigma = 0.0;
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            StdLargeVec<Real> &contact_mass_k = *(contact_mass_[k]);
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];

            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                Real corrected_W_ij = std::max(contact_neighborhood.W_ij_[n] - offset_W_ij_[k], Real(0));
                sigma += corrected_W_ij * contact_mass_k[contact_neighborhood.j_[n]];
            }
        }
        contact_density_[index_i] = sigma;
    };

  protected:
    StdLargeVec<Real> &mass_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
    StdVec<Real> offset_W_ij_;
};

/**
 * @class ShellContactDensity
 * @brief Computing the contact density due to shell contact using a
 * 		 surface integral being solved by Gauss-Legendre quadrature integration.
 */
class ShellContactDensity : public ContactDensityAccessor, public LocalDynamics, public ContactDynamicsData
{
  public:
    explicit ShellContactDensity(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~ShellContactDensity(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        /** shell contact interaction. */
        Real sigma = 0.0;
        Real contact_density_i = 0.0;

        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            StdLargeVec<Real> &contact_Vol_k = *(contact_Vol_[k]);
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                Real corrected_W_ij = std::max(contact_neighborhood.W_ij_[n] - offset_W_ij_[k], Real(0));
                sigma += corrected_W_ij * contact_Vol_k[contact_neighborhood.j_[n]];
            }
            constexpr Real heuristic_limiter = 0.1;
            // With heuristic_limiter, the maximum contact pressure is heuristic_limiter * K (Bulk modulus).
            // The contact pressure applied to fewer particles than on solids, yielding high acceleration locally,
            // which is one source of instability. Thus, we add a heuristic_limiter
            // to maintain enough contact pressure to prevent penetration while also maintaining stability.
            contact_density_i += heuristic_limiter * sigma * calibration_factor_[k];
        }
        contact_density_[index_i] = contact_density_i;
    };

  protected:
    Solid &solid_;
    Kernel *kernel_;

    Real particle_spacing_;
    StdVec<Real> calibration_factor_;
    StdVec<Real> offset_W_ij_;
    StdVec<StdLargeVec<Real> *> contact_Vol_;

    /** Abscissas and weights for Gauss-Legendre quadrature integration with n=3 nodes */
    const StdVec<Real> three_gaussian_points_ = {-0.7745966692414834, 0.0, 0.7745966692414834};
    const StdVec<Real> three_gaussian_weights_ = {0.5555555555555556, 0.8888888888888889, 0.5555555555555556};
};

/**
 * @class SelfContactForce
 * @brief Computing the self-contact force.
 */
class SelfContactForce : public ForcePrior, public SolidDataInner
{
  public:
    explicit SelfContactForce(SelfSurfaceContactRelation &self_contact_relation);
    virtual ~SelfContactForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Solid &solid_;
    StdLargeVec<Real> &mass_, &self_contact_density_, &Vol_;
    StdLargeVec<Vecd> &vel_;
    Real contact_impedance_;
};

/**
 * @class ContactForce
 * @brief Computing the contact force.
 */
class ContactForce : public ForcePrior, public ContactDynamicsData
{
  public:
    explicit ContactForce(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~ContactForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Solid &solid_;
    StdLargeVec<Real> &contact_density_, &Vol_, &mass_;
    StdVec<Solid *> contact_solids_;
    StdVec<StdLargeVec<Real> *> contact_contact_density_;
};

/**
 * @class ContactForceFromWall
 * @brief Computing the contact force from a rigid wall.
 *  Note that the body surface of the wall should be
 *  updated before computing the contact force.
 */
class ContactForceFromWall : public ForcePrior, public ContactWithWallData
{
  public:
    explicit ContactForceFromWall(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~ContactForceFromWall(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Solid &solid_;
    StdLargeVec<Real> &contact_density_, &Vol_;
};

/**
 * @class ContactForceToWall
 * @brief Computing contact force acting on a rigid wall.
 */
class ContactForceToWall : public ForcePrior, public ContactDynamicsData
{
  public:
    explicit ContactForceToWall(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~ContactForceToWall(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &Vol_;
    StdVec<Solid *> contact_solids_;
    StdVec<StdLargeVec<Real> *> contact_contact_density_;
};

/**
 * @class PairwiseFrictionFromWall
 * @brief Damping to wall by which the wall velocity is not updated
 * and the mass of wall particle is not considered.
 * Note that, currently, this class works only when the contact
 * bodies have the same resolution.
 */
class PairwiseFrictionFromWall : public LocalDynamics, public ContactWithWallData
{
  public:
    PairwiseFrictionFromWall(BaseContactRelation &contact_relation, Real eta);
    virtual ~PairwiseFrictionFromWall(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real eta_; /**< friction coefficient */
    StdLargeVec<Real> &Vol_, &mass_;
    StdLargeVec<Vecd> &vel_;
    StdVec<StdLargeVec<Vecd> *> wall_vel_n_, wall_n_;
};

/**
 * @class DynamicContactForceWithWall
 * @brief Computing the contact force with a rigid wall.
 *  Note that the body surface of the wall should be
 *  updated before computing the contact force.
 */
class DynamicContactForceWithWall : public ForcePrior, public ContactDynamicsData
{
  public:
    explicit DynamicContactForceWithWall(SurfaceContactRelation &solid_body_contact_relation, Real penalty_strength = 1.0);
    virtual ~DynamicContactForceWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Solid &solid_;
    StdLargeVec<Real> &Vol_;
    StdLargeVec<Vecd> &vel_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_, contact_n_;
    Real penalty_strength_;
    Real impedance_, reference_pressure_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // CONTACT_DYNAMICS_H
