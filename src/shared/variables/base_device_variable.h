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
 * @file 	base_device_variable.h
 * @brief 	Here gives classes for the base device variables used in simulation.
 * @details These variables are those discretized in spaces and time.
 * @author	Alberto Guarnieri and Xiangyu Hu
 */

#ifndef BASE_DEVICE_VARIABLE_H
#define BASE_DEVICE_VARIABLE_H

#include "base_variable.h"
#include "host_device_type_interface.h"
#include "memory_transfer.h"

namespace SPH
{
template <typename DataType>
class DeviceParticleDataContainer
{
  public: // type aliasing
    explicit DeviceParticleDataContainer(size_t size)
        : device_addr_(allocateDeviceData<Device<DataType>>(size)){};
    virtual ~DeviceParticleDataContainer()
    {
        freeDeviceData(device_addr_);
    };

    Device<DataType> *VariableAddress() { return device_addr_; }

  protected:
    Device<DataType> *device_addr_;
};
} // namespace SPH
#endif // BASE_DEVICE_VARIABLE_H
