#ifndef DIFFUSION_DYNAMICS_IMPLICIT_HPP
#define DIFFUSION_DYNAMICS_IMPLICIT_HPP

#include "diffusion_dynamics_implicit.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType, class CoefficientType>
template <typename CoefficientArg>
BaseDiffusionSplittingInner<DataType, CoefficientType>::
    BaseDiffusionSplittingInner(BaseInnerRelation &inner_relation,
                                const std::string &variable_name, const CoefficientArg &eta)
    : OperatorInner<DataType, DataType, CoefficientType>(
          inner_relation, variable_name, variable_name, eta),
      Vol_(this->particles_->Vol_), mass_(this->particles_->mass_),
      variable_(this->in_variable_) {}
//=================================================================================================//
template <typename DataType, class CoefficientType>
ErrorAndParameters<DataType>
BaseDiffusionSplittingInner<DataType, CoefficientType>::computeErrorAndParameters(size_t index_i, Real dt)
{
    Real Vol_i = Vol_[index_i];
    Real mass_i = mass_[index_i];
    DataType &variable_i = variable_[index_i];
    ErrorAndParameters<DataType> error_and_parameters;
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        // linear projection
        DataType variable_derivative = (variable_i - variable_[index_j]);
        Real parameter_b = 2.0 * this->coefficient_(index_i, index_j) *
                           inner_neighborhood.dW_ijV_j_[n] * Vol_i * dt / inner_neighborhood.r_ij_[n];

        error_and_parameters.error_ -= variable_derivative * parameter_b;
        error_and_parameters.a_ += parameter_b;
        error_and_parameters.c_ += parameter_b * parameter_b;
    }
    error_and_parameters.a_ -= mass_i;
    return error_and_parameters;
}
//=================================================================================================//
template <typename DataType, class CoefficientType>
void BaseDiffusionSplittingInner<DataType, CoefficientType>::
    updateStates(size_t index_i, Real dt, const ErrorAndParameters<DataType> &error_and_parameters)
{
    Real parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
    DataType parameter_k = error_and_parameters.error_ / (parameter_l + TinyReal);
    variable_[index_i] += parameter_k * error_and_parameters.a_;

    Real Vol_i = Vol_[index_i];
    Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        Real parameter_b = 2.0 * this->coefficient_(index_i, index_j) *
                           inner_neighborhood.dW_ijV_j_[n] * Vol_i * dt / inner_neighborhood.r_ij_[n];

        // predicted quantity at particle j
        DataType variable_j = variable_[index_j] - parameter_k * parameter_b;
        DataType variable_derivative = (variable_[index_i] - variable_j);

        // exchange in conservation form
        variable_[index_j] -= variable_derivative * parameter_b / mass_[index_j];
    }
}
//=================================================================================================//
template <typename DataType, class CoefficientType>
void BaseDiffusionSplittingInner<DataType, CoefficientType>::interaction(size_t index_i, Real dt)
{
    ErrorAndParameters<DataType> error_and_parameters = computeErrorAndParameters(index_i, dt);
    updateStates(index_i, dt, error_and_parameters);
}
//=================================================================================================//
template <typename DataType, class CoefficientType>
template <typename CoefficientArg>
BaseDiffusionSplittingWithWall<DataType, CoefficientType>::
    BaseDiffusionSplittingWithWall(ComplexRelation &complex_wall_relation,
                                   const std::string &variable_name, const CoefficientArg &eta)
    : BaseDiffusionSplittingInner<DataType, CoefficientType>(
          complex_wall_relation.getInnerRelation(), variable_name, eta),
      DissipationDataWithWall(complex_wall_relation.getContactRelation())
{
    for (size_t k = 0; k != DissipationDataWithWall::contact_particles_.size(); ++k)
    {
        wall_variable_.push_back(contact_particles_[k]->template getVariableByName<DataType>(variable_name));
    }
}
//=================================================================================================//
template <typename DataType, class CoefficientType>
ErrorAndParameters<DataType>
BaseDiffusionSplittingWithWall<DataType, CoefficientType>::
    computeErrorAndParameters(size_t index_i, Real dt)
{
    ErrorAndParameters<DataType> error_and_parameters =
        BaseDiffusionSplittingInner<DataType, CoefficientType>::computeErrorAndParameters(index_i, dt);

    const DataType &variable_i = this->variable_[index_i];
    Real Vol_i = this->Vol_[index_i];
    /** Contact interaction. */
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<DataType> &variable_k = *(this->wall_variable_[k]);
        Neighborhood &contact_neighborhood = (*DissipationDataWithWall::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];

            // linear projection
            DataType variable_derivative = (variable_i - variable_k[index_j]);
            Real parameter_b = 2.0 * this->coefficient_(index_i) *
                               contact_neighborhood.dW_ijV_j_[n] * Vol_i * dt / contact_neighborhood.r_ij_[n];

            error_and_parameters.error_ -= variable_derivative * parameter_b;
            error_and_parameters.a_ += parameter_b;
        }
    }
    return error_and_parameters;
}
//=================================================================================================//
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_IMPLICIT_HPP