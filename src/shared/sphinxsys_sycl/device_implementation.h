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
 * @file    device_implementation.h
 * @brief   TBD.
 * @author  Alberto Guarnieri and Xiangyu Hu
 */

#ifndef SPHINXSYS_DEVICE_IMPLEMENTATION_H
#define SPHINXSYS_DEVICE_IMPLEMENTATION_H

#include "execution_policy.h"
#include "ownership.h"

#include <memory>
#include <sycl/sycl.hpp>

namespace SPH::execution
{
template <class Device>
class DeviceImplementation
{
  public:
    DeviceImplementation() = default;

    template <class... Args>
    explicit DeviceImplementation(Args &&...device_args)
        : device_(std::make_shared<Device>(std::forward<Args>(device_args)...)) {}

    using KernelType = Device;

    inline Device *get_ptr() const { return device_.get(); }

    inline sycl::buffer<Device> &get_buffer()
    {
        if (!device_buffer_)
            device_buffer_ = std::make_unique<sycl::buffer<Device>>(device_, 1);
        return *device_buffer_;
    }

  private:
    SharedPtr<Device> device_;
    SharedPtr<sycl::buffer<Device>> device_buffer_;
};
} // namespace SPH::execution

#endif // SPHINXSYS_DEVICE_IMPLEMENTATION_H
