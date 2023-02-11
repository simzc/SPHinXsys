/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file 	io_observation.h
 * @brief 	Classes for input and output with vtk (Paraview) files.
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#pragma once

#include "io_base.h"
#include "general_interpolation.h"

namespace SPH
{

	/**
	 * @class ObservingQuantity
	 * @brief Observing a variable from contact bodies.
	 */
	template <typename VariableType>
	class ObservingQuantity : public InteractionDynamics<BaseInterpolation<VariableType>>
	{
	public:
		template <class ContactRelationType>
		explicit ObservingQuantity(ObservingContact<ContactRelationType> &observing_relation, const std::string &variable_name)
			: InteractionDynamics<BaseInterpolation<VariableType>>(observing_relation, variable_name),
			  observing_relation_(observing_relation)
		{
			this->interpolated_quantities_ = registerObservedQuantity(variable_name);
		};
		virtual ~ObservingQuantity(){};

		/** The sequential function for executing the average operations on particles and their neighbors. */
		virtual void exec(Real dt = 0.0) override
		{
			observing_relation_.updateConfiguration();
			InteractionDynamics<BaseInterpolation<VariableType>>::exec(dt);
		};
		/** The parallel function for executing the average operations on particles and their neighbors. */
		virtual void parallel_exec(Real dt = 0.0) override
		{
			observing_relation_.updateConfiguration();
			InteractionDynamics<BaseInterpolation<VariableType>>::parallel_exec(dt);
		};

	protected:
		StdLargeVec<VariableType> observed_quantities_;
		BaseContactRelation &observing_relation_;

		/** Register the  observed variable if the variable name is new.
		 * If the variable is registered already, the registered variable will be returned. */
		StdLargeVec<VariableType> *registerObservedQuantity(const std::string &variable_name)
		{
			BaseParticles *particles = this->particles_;
			constexpr int type_index = DataTypeIndex<VariableType>::value;
			if (particles->all_variable_maps_[type_index].find(variable_name) == particles->all_variable_maps_[type_index].end())
			{
				particles->registerVariable(observed_quantities_, variable_name, ZeroData<VariableType>::value);
				return &observed_quantities_;
			}
			return particles->getVariableByName<VariableType>(variable_name);
		};
	};

	/**
	 * @class  BaseObservation
	 * @brief base class for sampling from simulation
	 */
	class BaseObservation : public BaseIO
	{
	protected:
		PltEngine plt_engine_;
		std::string dynamics_range_name_;
		const std::string quantity_name_;
		std::string filefullpath_output_;

	public:
		BaseObservation(IOEnvironment &io_environment, const std::string &dynamics_range_name,
						const std::string &quantity_name)
			: BaseIO(io_environment), plt_engine_(), dynamics_range_name_(dynamics_range_name),
			  quantity_name_(quantity_name)
		{
			filefullpath_output_ = io_environment_.output_folder_ + "/" + dynamics_range_name_ + "_" + quantity_name_ + ".dat";
		};
	};

	/**
	 * @class ObservedQuantityRecording
	 * @brief write files for observed quantity
	 */
	template <typename VariableType>
	class ObservedQuantityRecording : public ObservingQuantity<VariableType>,
									  public BaseObservation
	{
	protected:
		SPHBody &observer_;
		BaseParticles &base_particles_;

	public:
		VariableType type_indicator_; /*< this is an indicator to identify the variable type. */

		template <class ContactRelationType>
		ObservedQuantityRecording(const std::string &quantity_name, IOEnvironment &io_environment,
<<<<<<< HEAD
								  ObservingContact<ContactRelationType> &observing_relation)
			: ObservingQuantity<VariableType>(observing_relation, quantity_name),
			  BaseObservation(io_environment, observing_relation.sph_body_.getName(), quantity_name),
			  observer_(observing_relation.sph_body_), base_particles_(observer_.getBaseParticles())
=======
								  BaseContactRelation &contact_relation)
			: BodyStatesRecording(io_environment, contact_relation.getSPHBody()),
			  ObservingAQuantity<VariableType>(contact_relation, quantity_name),
			  observer_(contact_relation.getSPHBody()), plt_engine_(),
			  base_particles_(observer_.getBaseParticles()), 
			  dynamics_range_name_(contact_relation.getSPHBody().getName()),
			  quantity_name_(quantity_name)
>>>>>>> xiangyu/revise_base_dynamics
		{
			std::ofstream out_file(this->filefullpath_output_.c_str(), std::ios::app);
			out_file << "run_time"
					 << "   ";
			for (size_t i = 0; i != base_particles_.total_real_particles_; ++i)
			{
				std::string quantity_name_i = this->quantity_name_ + "[" + std::to_string(i) + "]";
				this->plt_engine_.writeAQuantityHeader(out_file, (*this->interpolated_quantities_)[i], quantity_name_i);
			}
			out_file << "\n";
			out_file.close();
		};
		virtual ~ObservedQuantityRecording(){};

		virtual void writeToFileByStep() override
		{
			if (this->checkToRecord())
			{
				this->parallel_exec();
				std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
				out_file << GlobalStaticVariables::physical_time_ << "   ";
				for (size_t i = 0; i != base_particles_.total_real_particles_; ++i)
				{
					plt_engine_.writeAQuantity(out_file, (*this->interpolated_quantities_)[i]);
				}
				out_file << "\n";
				out_file.close();
			}
		}

		StdLargeVec<VariableType> *getObservedQuantity()
		{
			return this->interpolated_quantities_;
		}
	};

	/**
	 * @class ReducedQuantityRecording
	 * @brief write reduced quantity of a body
	 */
	template <class ReduceMethodType>
	class ReducedQuantityRecording : public ReduceMethodType,
									 public BaseObservation
	{
	public:
		/*< deduce variable type from reduce method. */
		using VariableType = typename ReduceMethodType::ReduceReturnType;
		VariableType type_indicator_; /*< this is an indicator to identify the variable type. */

	public:
		template <typename... ConstructorArgs>
		ReducedQuantityRecording(IOEnvironment &io_environment, ConstructorArgs &&...args)
			: ReduceMethodType(std::forward<ConstructorArgs>(args)...),
			  BaseObservation(io_environment, this->DynamicsRangeName(), this->ReducedQuantityName())
		{
			std::ofstream out_file(this->filefullpath_output_.c_str(), std::ios::app);
			out_file << "\"run_time\""
					 << "   ";
			plt_engine_.writeAQuantityHeader(out_file, this->Reference(), this->quantity_name_);
			out_file << "\n";
			out_file.close();
		};
		virtual ~ReducedQuantityRecording(){};

		virtual void writeToFileByStep() override
		{
			if (this->checkToRecord())
			{
				std::ofstream out_file(this->filefullpath_output_.c_str(), std::ios::app);
				out_file << GlobalStaticVariables::physical_time_ << "   ";
				plt_engine_.writeAQuantity(out_file, this->parallel_exec());
				out_file << "\n";
				out_file.close();
			}
		};
	};
}
