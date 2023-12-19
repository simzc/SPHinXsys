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
 * @file 	host_device_type_interface.h
 * @brief 	This is the date type definition for SPHinXsys.
 * @author	Xiangyu Hu
 */

#ifndef HOST_DEVICE_TYPE_INTERFACE_H
#define HOST_DEVICE_TYPE_INTERFACE_H

#include "base_data_type.h"
#include "base_device_data_type.h"

namespace SPH
{
template <typename T>
struct HostType
{
    using type = T;
};

template <typename T>
using Host = typename HostType<T>::type;

template <typename T>
struct DeviceType
{
    using type = T;
};

template <>
struct DeviceType<Real>
{
    using type = DeviceReal;
};

template <>
struct DeviceType<Vec2d>
{
    using type = DeviceVec2d;
};

template <>
struct DeviceType<Vec3d>
{
    using type = DeviceVec3d;
};

template <>
struct DeviceType<Mat2d>
{
    using type = DeviceMat2d;
};

template <>
struct DeviceType<Mat3d>
{
    using type = DeviceMat3d;
};

template <typename T>
using Device = typename DeviceType<T>::type;

template <typename Type, class Enable = void>
struct DataTypeEquivalence
{
    static constexpr bool type_defined = false;
    static_assert("Type non recognized as host or device type.");
};

template <class CheckType, class Type1, class Type2>
using enable_if_is_either_t =
    std::enable_if_t<std::disjunction_v<std::is_same<CheckType, Type1>, std::is_same<CheckType, Type2>>>;

template <class TypeReal>
struct DataTypeEquivalence<TypeReal, enable_if_is_either_t<TypeReal, Real, DeviceReal>>
{
    static constexpr bool type_defined = true;
    using host_type = Real;
    using device_type = DeviceReal;
};

template <class TypeVec2d>
struct DataTypeEquivalence<TypeVec2d, enable_if_is_either_t<TypeVec2d, Vec2d, DeviceVec2d>>
{
    static constexpr bool type_defined = true;
    using host_type = Vec2d;
    using device_type = DeviceVec2d;
};

template <class TypeVec3d>
struct DataTypeEquivalence<TypeVec3d, enable_if_is_either_t<TypeVec3d, Vec3d, DeviceVec3d>>
{
    static constexpr bool type_defined = true;
    using host_type = Vec3d;
    using device_type = DeviceVec3d;
};

template <class TypeMat2d>
struct DataTypeEquivalence<TypeMat2d, enable_if_is_either_t<TypeMat2d, Vec2d, DeviceMat2d>>
{
    static constexpr bool type_defined = true;
    using host_type = Mat2d;
    using device_type = DeviceMat2d;
};

template <class TypeMat3d>
struct DataTypeEquivalence<TypeMat3d, enable_if_is_either_t<TypeMat3d, Mat3d, DeviceMat3d>>
{
    static constexpr bool type_defined = true;
    using host_type = Mat3d;
    using device_type = DeviceMat3d;
};

template <class CheckType, class HostOrDeviceType>
using enable_both_host_device_t =
    enable_if_is_either_t<CheckType, typename DataTypeEquivalence<HostOrDeviceType>::host_type,
                          typename DataTypeEquivalence<HostOrDeviceType>::device_type>;

} // namespace SPH

#endif // HOST_DEVICE_TYPE_INTERFACE_H
