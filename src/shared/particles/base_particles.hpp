/**
 * @file 	base_particles.hpp
 * @brief 	This is the implementation of the template functions in base_particles.h
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_PARTICLES_HPP
#define BASE_PARTICLES_HPP

#include "base_particles.h"
#include "particle_dynamics_algorithms.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::registerSingleValueVariable(const std::string &name, DataType initial_value)
{
    SingleValueVariable<DataType> *variable = findVariableByName<DataType>(single_value_variables_, name);

    return variable != nullptr
               ? variable->ValueAddress()
               : addVariableToAssemble<DataType>(single_value_variables_,
                                                 single_value_variable_ptrs_, name, initial_value)
                     ->ValueAddress();
}
//=================================================================================================//
template <typename DataType>
DataType *BaseParticles::getSingleValueVariableByName(const std::string &name)
{
    SingleValueVariable<DataType> *variable = findVariableByName<DataType>(single_value_variables_, name);

    if (variable != nullptr)
    {
        return variable->ValueAddress();
    }

    std::cout << "\nError: the variable '" << name << "' is not registered!\n";
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    return nullptr;
}
//=================================================================================================//
template <typename DataType>
StdLargeVec<DataType> *BaseParticles::
    registerSharedVariable(const std::string &name, DataType initial_value)
{
    return registerSharedVariable<DataType>([initial_value](size_t i) -> DataType
                                            { return initial_value; },
                                            name);
}
//=================================================================================================//
template <typename DataType>
StdLargeVec<DataType> *BaseParticles::
    registerSharedVariable(const std::string &name, const std::string &target_name)
{
    DiscreteVariable<DataType> *target_variable =
        findVariableByName<DataType, DiscreteVariable>(all_discrete_variables_, target_name);

    if (target_variable != nullptr)
    {
        StdLargeVec<DataType> &variable_data = *target_variable->ValueAddress();
        return registerSharedVariable<DataType>([&variable_data](size_t i) -> DataType
                                                { return variable_data[i]; },
                                                name);
    }

    std::cout << "\nError: the target variable '" << target_name << "' is not registered!\n";
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    return nullptr;
}
//=================================================================================================//
template <typename DataType, class InitializationFunction>
StdLargeVec<DataType> *BaseParticles::
    registerSharedVariable(const InitializationFunction &initialization, const std::string &name)
{

    DiscreteVariable<DataType> *variable =
        findVariableByName<DataType, DiscreteVariable>(all_discrete_variables_, name);

    constexpr int type_index = DataTypeIndex<DataType>::value;
    if (variable == nullptr)
    {
        DiscreteVariable<DataType> *new_variable =
            addVariableToAssemble<DataType, DiscreteVariable>(all_discrete_variables_, all_discrete_variable_ptrs_, name);
        StdLargeVec<DataType> *variable_addrs = new_variable->ValueAddress();
        for (size_t i = 0; i != particles_bound_; ++i)
        {
            (*variable_addrs)[i] = initialization(i);
        }
        std::get<type_index>(all_particle_data_).push_back(variable_addrs);
        return variable_addrs;
    }
    else
    {
        return variable->ValueAddress();
    }
}
//=================================================================================================//
template <typename DataType>
StdLargeVec<DataType> *BaseParticles::getVariableByName(const std::string &name)
{
    DiscreteVariable<DataType> *variable =
        findVariableByName<DataType, DiscreteVariable>(all_discrete_variables_, name);

    if (variable != nullptr)
    {
        variable->ValueAddress();
    }

    std::cout << "\nError: the variable '" << name << "' is not registered!\n";
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    return nullptr;
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToList(ParticleVariables &variable_set, const std::string &name)
{
    DiscreteVariable<DataType> *variable =
        findVariableByName<DataType, DiscreteVariable>(all_discrete_variables_, name);

    if (variable != nullptr)
    {
        DiscreteVariable<DataType> *listed_variable =
            findVariableByName<DataType, DiscreteVariable>(variable_set, name);

        if (listed_variable == nullptr)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            std::get<type_index>(variable_set).push_back(variable);
        }
    }
    else
    {
        std::cout << "\n Error: the variable '" << name << "' to write is not particle data!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToWrite(const std::string &name)
{
    addVariableToList<DataType>(variables_to_write_, name);
}
//=================================================================================================//
template <class DerivedVariableMethod, class... Ts>
void BaseParticles::addDerivedVariableToWrite(Ts &&...args)
{
    SimpleDynamics<DerivedVariableMethod> *derived_data =
        derived_particle_data_.createPtr<SimpleDynamics<DerivedVariableMethod>>(sph_body_, std::forward<Ts>(args)...);
    derived_variables_.push_back(derived_data);
    using DerivedDataType = typename DerivedVariableMethod::DerivedDataType;
    addVariableToList<DerivedDataType>(variables_to_write_, derived_data->name_);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToRestart(const std::string &name)
{
    addVariableToList<DataType>(variables_to_restart_, name);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::addVariableToReload(const std::string &name)
{
    addVariableToList<DataType>(variables_to_reload_, name);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::registerSortableVariable(const std::string &name)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, name);

    if (variable != nullptr)
    {
        DiscreteVariable<DataType> *listed_variable = findVariableByName<DataType>(sortable_variables_, name);

        if (listed_variable == nullptr)
        {
            constexpr int type_index = DataTypeIndex<DataType>::value;
            std::get<type_index>(sortable_variables_).push_back(variable);
            StdLargeVec<DataType> *variable_data = std::get<type_index>(all_particle_data_)[variable->IndexInContainer()];
            std::get<type_index>(sortable_data_).push_back(variable_data);
        }
    }
    else
    {
        std::cout << "\n Error: the variable '" << name << "' to write is not particle data!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
template <typename SequenceMethod>
void BaseParticles::sortParticles(SequenceMethod &sequence_method)
{
    StdLargeVec<size_t> &sequence = sequence_method.computingSequence(*this);
    particle_sorting_.sortingParticleData(sequence.data(), total_real_particles_);
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::ResizeParticles::
operator()(DataContainerAddressKeeper<StdLargeVec<DataType>> &data_keeper, size_t new_size)
{
    for (size_t i = 0; i != data_keeper.size(); ++i)
    {
        data_keeper[i]->resize(new_size, ZeroData<DataType>::value);
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::CopyParticleData::
operator()(DataContainerAddressKeeper<StdLargeVec<DataType>> &data_keeper, size_t index, size_t another_index)
{
    for (size_t i = 0; i != data_keeper.size(); ++i)
    {
        (*data_keeper[i])[index] = (*data_keeper[i])[another_index];
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::WriteAParticleVariableToXml::
operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables, ParticleData &all_particle_data)
{
    for (size_t i = 0; i != variables.size(); ++i)
    {
        size_t index = 0;
        StdLargeVec<DataType> &variable_data = *variables[i]->ValueAddress();
        for (auto child = xml_parser_.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            xml_parser_.setAttributeToElement(child, variables[i]->Name(), variable_data[index]);
            index++;
        }
    }
}
//=================================================================================================//
template <typename DataType>
void BaseParticles::ReadAParticleVariableFromXml::
operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables, ParticleData &all_particle_data)
{
    for (size_t i = 0; i != variables.size(); ++i)
    {
        size_t index = 0;
        StdLargeVec<DataType> &variable_data = *variables[i]->ValueAddress();
        for (auto child = xml_parser_.first_element_->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            xml_parser_.queryAttributeValue(child, variables[i]->Name(), variable_data[index]);
            index++;
        }
    }
}
//=================================================================================================//
template <typename StreamType>
void BaseParticles::writeParticlesToVtk(StreamType &output_stream)
{
    size_t total_real_particles = total_real_particles_;

    // write sorted particles ID
    output_stream << "    <DataArray Name=\"SortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
    output_stream << "    ";
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        output_stream << i << " ";
    }
    output_stream << std::endl;
    output_stream << "    </DataArray>\n";

    // write unsorted particles ID
    output_stream << "    <DataArray Name=\"UnsortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
    output_stream << "    ";
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        output_stream << unsorted_id_[i] << " ";
    }
    output_stream << std::endl;
    output_stream << "    </DataArray>\n";

    // write integers
    constexpr int type_index_int = DataTypeIndex<int>::value;
    for (DiscreteVariable<int> *variable : std::get<type_index_int>(variables_to_write_))
    {
        StdLargeVec<int> &variable_data = *variable->ValueAddress();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Int32\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            output_stream << std::fixed << std::setprecision(9) << variable_data[i] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write scalars
    constexpr int type_index_Real = DataTypeIndex<Real>::value;
    for (DiscreteVariable<Real> *variable : std::get<type_index_Real>(variables_to_write_))
    {
        StdLargeVec<Real> &variable_data = *variable->ValueAddress();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            output_stream << std::fixed << std::setprecision(9) << variable_data[i] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write vectors
    constexpr int type_index_Vecd = DataTypeIndex<Vecd>::value;
    for (DiscreteVariable<Vecd> *variable : std::get<type_index_Vecd>(variables_to_write_))
    {
        StdLargeVec<Vecd> &variable_data = *variable->ValueAddress();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Vec3d vector_value = upgradeToVec3d(variable_data[i]);
            output_stream << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }

    // write matrices
    constexpr int type_index_Matd = DataTypeIndex<Matd>::value;
    for (DiscreteVariable<Matd> *variable : std::get<type_index_Matd>(variables_to_write_))
    {
        StdLargeVec<Matd> &variable_data = *variable->ValueAddress();
        output_stream << "    <DataArray Name=\"" << variable->Name() << "\" type= \"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
        output_stream << "    ";
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            Mat3d matrix_value = upgradeToMat3d(variable_data[i]);
            for (int k = 0; k != 3; ++k)
            {
                Vec3d col_vector = matrix_value.col(k);
                output_stream << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
            }
        }
        output_stream << std::endl;
        output_stream << "    </DataArray>\n";
    }
}
//=================================================================================================//
template <typename DataType>
BaseDerivedVariable<DataType>::
    BaseDerivedVariable(SPHBody &sph_body, const std::string &name)
    : name_(name),
      derived_variable_(sph_body.getBaseParticles().registerSharedVariable<DataType>(name)){};
//=================================================================================================//
} // namespace SPH
#endif // BASE_PARTICLES_HPP