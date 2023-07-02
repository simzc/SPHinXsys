#include "general_dynamics_extra.h"
#include <limits>

namespace SPH
{
//=================================================================================================//
ConstraintTotalScalarAmount::
    ConstraintTotalScalarAmount(SPHBody &sph_body, const std::string &variable_name)
    : LocalDynamics(sph_body), GeneralDataDelegateSimple(sph_body),
      variable_(*particles_->getVariableByName<Real>(variable_name)),
      total_scalar_(sph_body, variable_name), is_initialized_(false),
      inital_total_(0), increment_(0) {}
//=================================================================================================//
void ConstraintTotalScalarAmount::setupInitialScalarAmount()
{
    inital_total_ = total_scalar_.exec();
    is_initialized_ = true;
}
//=================================================================================================//
void ConstraintTotalScalarAmount::setupDynamics(Real dt)
{
    if (is_initialized_)
    {
        Real current_total = total_scalar_.exec();
        increment_ = (inital_total_ - current_total) / current_total;
    }
    else
    {
        std::cout << "\n Error: inital scalar amount is not defined!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
void ConstraintTotalScalarAmount::update(size_t index_i, Real dt)
{
    variable_[index_i] *= (1.0 + increment_);
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//