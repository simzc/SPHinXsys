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
 * @file 	base_variable.h
 * @brief 	Here gives classes for the base variables used in simulation.
 * @details These variables are those discretized in spaces and time.
 * @author	Xiangyu Hu
 */

#ifndef BASE_VARIABLES_H
#define BASE_VARIABLES_H

#include "base_data_package.h"

namespace SPH
{
template <typename... Parameters>
class Variable;

class BaseVariable
{
  public:
    explicit BaseVariable(const std::string &name) : name_(name){};
    virtual ~BaseVariable(){};
    std::string Name() const { return name_; };

  private:
    const std::string name_;
};

template <typename DataType>
class Variable<StdLargeVec<DataType>> : public BaseVariable
{
  public:
    explicit Variable(const std::string &name) : BaseVariable(name){};
    virtual ~Variable(){};

    StdLargeVec<DataType> *ValueAddress() { return &value_; };

  private:
    StdLargeVec<DataType> value_;
};

template <typename ContainedValueType>
using ValueAddressKeeper = StdVec<ContainedValueType *>;

template <typename ContainedValueType>
using VariableAddressKeeper = StdVec<Variable<ContainedValueType> *>;

template <typename ContainedValueType>
using VariableUniquePtrsKeeper = UniquePtrsKeeper<Variable<ContainedValueType>>;

template <template <typename> typename KeeperType, template <typename> typename ContainerType>
using DataAssemble = std::tuple<KeeperType<ContainerType<Real>>,
                                KeeperType<ContainerType<Vec2d>>,
                                KeeperType<ContainerType<Vec3d>>,
                                KeeperType<ContainerType<Mat2d>>,
                                KeeperType<ContainerType<Mat3d>>,
                                KeeperType<ContainerType<int>>>;

template <template <typename> typename ContainerType>
using ValueAddressAssemble = DataAssemble<ValueAddressKeeper, ContainerType>;

template <template <typename> typename ContainerType>
using VariableAddressAssemble = DataAssemble<VariableAddressKeeper, ContainerType>;

template <template <typename> typename ContainerType>
using VariableUniquePtrAssemble = DataAssemble<VariableUniquePtrsKeeper, ContainerType>;

template <typename DataType, template <typename> class ContainerType>
Variable<ContainerType<DataType>> *findVariableByName(VariableAddressAssemble<ContainerType> &variable_addrs_assemble,
                                                      const std::string &name)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    auto &variables = std::get<type_index>(variable_addrs_assemble);
    auto result = std::find_if(variables.begin(), variables.end(),
                               [&](auto &variable) -> bool
                               { return variable->Name() == name; });

    return result != variables.end() ? *result : nullptr;
};

template <typename DataType, template <typename> class ContainerType, typename... Args>
ContainerType<DataType> *addVariableToAssemble(ValueAddressAssemble<ContainerType> &value_addrs_assemble,
                                               VariableAddressAssemble<ContainerType> &variable_addrs_assemble,
                                               VariableUniquePtrAssemble<ContainerType> &variable_ptr_assemble,
                                               const std::string &name, Args &&...args)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    UniquePtrsKeeper<Variable<ContainerType<DataType>>> &variable_ptrs = std::get<type_index>(variable_ptr_assemble);
    Variable<ContainerType<DataType>> *new_variable =
        variable_ptrs.template createPtr<Variable<ContainerType<DataType>>>(name, std::forward<Args>(args)...);
    std::get<type_index>(variable_addrs_assemble).push_back(new_variable);
    ContainerType<DataType> *new_value_addrs = new_variable->ValueAddress();
    std::get<type_index>(value_addrs_assemble).push_back(new_value_addrs);
    return new_value_addrs;
};

template <typename DataType>
class SingleVariable : public BaseVariable
{
  public:
    SingleVariable(const std::string &name, DataType &value)
        : BaseVariable(name), value_(value){};
    virtual ~SingleVariable(){};

    DataType *ValueAddress() { return &value_; };

  private:
    DataType value_;
};

template <typename DataType>
class DiscreteVariable : public BaseVariable
{
  public:
    DiscreteVariable(const std::string &name, size_t index)
        : BaseVariable(name), index_in_container_(index){};
    virtual ~DiscreteVariable(){};

    size_t IndexInContainer() const { return index_in_container_; };

  private:
    size_t index_in_container_;
};

template <typename DataType, template <typename VariableDataType> class VariableType>
VariableType<DataType> *findVariableByName(DataContainerAddressAssemble<VariableType> &assemble,
                                           const std::string &name)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    auto &variables = std::get<type_index>(assemble);
    auto result = std::find_if(variables.begin(), variables.end(),
                               [&](auto &variable) -> bool
                               { return variable->Name() == name; });

    return result != variables.end() ? *result : nullptr;
};

template <typename DataType, template <typename VariableDataType> class VariableType, typename... Args>
VariableType<DataType> *addVariableToAssemble(DataContainerAddressAssemble<VariableType> &assemble,
                                              DataContainerUniquePtrAssemble<VariableType> &ptr_assemble, Args &&...args)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    UniquePtrsKeeper<VariableType<DataType>> &variable_ptrs = std::get<type_index>(ptr_assemble);
    VariableType<DataType> *new_variable =
        variable_ptrs.template createPtr<VariableType<DataType>>(std::forward<Args>(args)...);
    std::get<type_index>(assemble).push_back(new_variable);
    return new_variable;
};
} // namespace SPH
#endif // BASE_VARIABLES_H
