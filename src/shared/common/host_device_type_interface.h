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
 * @brief 	This is the data type definition for SPHinXsys.
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

} // namespace SPH

#endif // HOST_DEVICE_TYPE_INTERFACE_H
