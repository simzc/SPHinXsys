#ifndef BASE_PARTICLE_GENERATOR_HPP
#define BASE_PARTICLE_GENERATOR_HPP

#include "base_particle_generator.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
StdLargeVec<DataType> *ParticleGenerator<Base>::registerGeometricVariable(const std::string &name)
{
    auto *variable = findVariableByName<DataType, StdLargeVec>(geometric_variables_, name);
    if (variable == nullptr)
    {
        return addVariableToAssemble<DataType, StdLargeVec>(geometric_data_, geometric_variables_,
                                                            geometric_variables_ptrs_, name);
    }
    else
    {
        return variable->ValueAddress();
    }
}
//=================================================================================================//
} // namespace SPH
#endif // BASE_PARTICLE_GENERATOR_HPP
