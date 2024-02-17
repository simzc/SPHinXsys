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
 * @file base_particles.h
 * @brief This is the base class of SPH particles. The basic data of the particles
 * is saved in separated large vectors. Each derived class will introduce several extra
 * vectors for the new data. Note that there is no class of single particle.
 * @author	Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef BASE_PARTICLES_H
#define BASE_PARTICLES_H

#include "base_data_package.h"
#include "base_material.h"
#include "base_variable.h"
#include "particle_sorting.h"
#include "sph_data_containers.h"
#include "xml_parser.h"

#include <fstream>

namespace SPH
{

class SPHBody;
class BaseMaterial;
class BodySurface;
template <class ReturnType>
class BaseDynamics;

/**
 * @class BaseParticles
 * @brief Particles with essential (geometric and kinematic) data.
 * @details There are three types of particles， all particles of a same type are saved with continuous memory segments.
 * The first is real particle whose states are updated by particle dynamics.
 * The second is buffer particle whose states are not updated by particle dynamics.
 * Buffer particles are saved behind real particles.
 * The global value of total_real_particles_ separate the real and buffer particles.
 * Buffer particles may be switched from real particles or to real particles.
 * The total number of real particles and buffer particles gives the real_particles_bound_,
 * i.e. the maximum possible number of real particles.
 * As the memory for both particles are continuous, such switch is achieved at the memory boundary sequentially.
 * The basic idea is swap the data of the last real particle with the one will be switched,
 * and then switch this swapped last particle as buffer particle by decrease the total_real_particles_ by one.
 * Switch from buffer particle to real particle is easy. One just need to assign expect state to
 * the first buffer particle and increase total_real_particles_ by one.
 * The third is ghost particle which is created according to a corresponding real particle.
 * The states of ghost particle are updated according to boundary condition and those of the real particle.
 * The ghost particles are saved behind the buffer particles within one or more ghost bounds.
 * All particles are bounded by particle_bound_, which is the total number of particles in all types.
 * In SPHinXsys, the discrete variables (state of each particle) and single variables (state for all particles)
 * registered in general particle data (ParticleData).
 * Generally, all discrete variable should be owned by a BaseParticles object
 * so that other objects can use it by the function getVariableByName.
 * A shared discrete variable can also be defined by several objects
 * and a unique discrete variable is only used by the object who registers the variable.
 * According the memory span of the variable, the discrete variable can be classified into
 * sortable variable and derivable variables.
 * Sortable variables keep their memories from the beginning to the end of the simulation,
 * that is why they keeps their values after particle sorting.
 * Derivable variables are updated from sortable variables and other derivable variables.
 * Therefore, the memory span of derivable variables is only from its updating to the end of a single simulation step.
 * There are also (diagnostic) variables which are derived and not used for evolution of the system but only for output or visualization.
 */
class BaseParticles
{
  private:
    DataContainerUniquePtrAssemble<DiscreteVariable> all_discrete_variable_ptrs_;
    DataContainerUniquePtrAssemble<StdLargeVec> shared_particle_data_ptrs_;
    DataContainerUniquePtrAssemble<SingleVariable> all_global_variable_ptrs_;
    UniquePtrsKeeper<BaseDynamics<void>> diagnostic_particle_data_;

  public:
    explicit BaseParticles(SPHBody &sph_body, BaseMaterial *base_material);
    virtual ~BaseParticles(){};

    StdLargeVec<Vecd> pos_;         /**< Position */
    StdLargeVec<Vecd> vel_;         /**< Velocity */
    StdLargeVec<Vecd> force_;       /**< Force induced by pressure- or stress */
    StdLargeVec<Vecd> force_prior_; /**< Other, such as gravity and viscous, forces computed before force_ */

    StdLargeVec<Real> Vol_;      /**< Volumetric measure, also area and length of surface and linear particle */
    StdLargeVec<Real> rho_;      /**< Density */
    StdLargeVec<Real> mass_;     /**< Mass*/
    StdLargeVec<int> indicator_; /**< particle indicator: 0 for bulk, 1 for free surface indicator, other to be defined */
    //----------------------------------------------------------------------
    // Global information for defining particle groups
    //----------------------------------------------------------------------
    size_t total_real_particles_; /**< Global value of the present number of real particles. */
    size_t real_particles_bound_; /**< Maximum possible number of real particles. */
    size_t total_bound_;          /**< Total number of particles for all types. */

    SPHBody &getSPHBody() { return sph_body_; };
    BaseMaterial &getBaseMaterial() { return base_material_; };
    ParticleData &getAllParticleData() { return all_particle_data_; };
    /** initialize other variables after the particles are generated */
    virtual void initializeOtherVariables();
    //----------------------------------------------------------------------
    //		Generalized particle manipulation
    //----------------------------------------------------------------------
    void initializeAllParticlesBounds();
    void increaseAllParticlesBounds(size_t buffer_size);
    void copyFromAnotherParticle(size_t index, size_t another_index);
    void updateGhostParticle(size_t ghost_index, size_t index);
    void switchToBufferParticle(size_t index);
    //----------------------------------------------------------------------
    //		Parameterized management on generalized particle data
    //----------------------------------------------------------------------
    template <typename DataType>
    void registerVariable(StdLargeVec<DataType> &variable_addrs, const std::string &variable_name,
                          DataType initial_value = ZeroData<DataType>::value);
    template <typename DataType, class InitializationFunction>
    void registerVariable(StdLargeVec<DataType> &variable_addrs, const std::string &variable_name,
                          const InitializationFunction &initialization);
    template <typename DataType>
    StdLargeVec<DataType> *registerSharedVariable(
        const std::string &variable_name, const DataType &default_value = ZeroData<DataType>::value);
    template <typename DataType>
    StdLargeVec<DataType> *getVariableByName(const std::string &variable_name);
    ParticleVariables &AllDiscreteVariables() { return all_discrete_variables_; };

    template <typename DataType>
    DataType *registerSingleVariable(const std::string &variable_name,
                                     DataType initial_value = ZeroData<DataType>::value);
    template <typename DataType>
    DataType *getSingleVariableByName(const std::string &variable_name);
    //----------------------------------------------------------------------
    //		Manage subsets of particle variables
    //----------------------------------------------------------------------
    template <typename DataType>
    void addVariableToList(ParticleVariables &variable_set, const std::string &variable_name);
    template <typename DataType>
    void addVariableToWrite(const std::string &variable_name);
    template <typename DataType>
    void addVariableToRestart(const std::string &variable_name);
    inline const ParticleVariables &getVariablesToRestart() const { return variables_to_restart_; }
    template <typename DataType>
    void addVariableToReload(const std::string &variable_name);
    inline const ParticleVariables &getVariablesToReload() const { return variables_to_reload_; }

    template <class DiagnosticVariableMethod, class... Ts>
    void addDiagnosticVariableToWrite(Ts &&...);
    void computeDiagnosticVariables();
    //----------------------------------------------------------------------
    //		Particle data for sorting
    //----------------------------------------------------------------------
    StdLargeVec<size_t> unsorted_id_; /**< the ids assigned just after particle generated. */
    StdLargeVec<size_t> sorted_id_;   /**< the sorted particle ids of particles from unsorted ids. */
    StdLargeVec<size_t> sequence_;    /**< the sequence referred for sorting. */
    ParticleData sortable_data_;
    ParticleVariables sortable_variables_;
    ParticleSorting particle_sorting_;

    template <typename DataType>
    void registerSortableVariable(const std::string &variable_name);
    template <typename SequenceMethod>
    void sortParticles(SequenceMethod &sequence_method);
    //----------------------------------------------------------------------
    //		Particle data ouput functions
    //----------------------------------------------------------------------
    template <typename StreamType>
    void writeParticlesToVtk(StreamType &output_stream);
    void writeParticlesToPltFile(std::ofstream &output_file);
    virtual void writeSurfaceParticlesToVtuFile(std::ostream &output_file, BodySurface &surface_particles);
    void resizeXmlDocForParticles(XmlParser &xml_parser);
    void writeParticlesToXmlForRestart(std::string &filefullpath);
    void readParticleFromXmlForRestart(std::string &filefullpath);
    void writeToXmlForReloadParticle(std::string &filefullpath);
    void readFromXmlForReloadParticle(std::string &filefullpath);
    XmlParser *getReloadXmlParser() { return &reload_xml_parser_; };
    virtual BaseParticles *ThisObjectPtr() { return this; };
    //----------------------------------------------------------------------
    //		Relation relate volume, surface and linear particles
    //----------------------------------------------------------------------
    virtual Real ParticleVolume(size_t index) { return Vol_[index]; }
    virtual Real ParticleSpacing(size_t index) { return std::pow(Vol_[index], 1.0 / Real(Dimensions)); }

  protected:
    SPHBody &sph_body_;
    std::string body_name_;
    BaseMaterial &base_material_;
    XmlParser restart_xml_parser_;
    XmlParser reload_xml_parser_;
    ParticleData all_particle_data_;
    ParticleVariables all_discrete_variables_;
    SingleVariables all_single_variables_;
    ParticleVariables variables_to_write_;
    ParticleVariables variables_to_restart_;
    ParticleVariables variables_to_reload_;
    StdVec<BaseDynamics<void> *> diagnostic_variables_;

    virtual void writePltFileHeader(std::ofstream &output_file);
    virtual void writePltFileParticleData(std::ofstream &output_file, size_t index);
    void writeParticlesToXml(XmlParser &xml_parser, ParticleVariables &particle_variables);
    void readParticleFromXml(XmlParser &xml_parser, ParticleVariables &particle_variables);
    //----------------------------------------------------------------------
    //		Small structs for generalize particle operations
    //----------------------------------------------------------------------
    template <typename DataType>
    class resizeParticleData
    {
        ParticleData &all_particle_data_;

      public:
        resizeParticleData(ParticleData &all_particle_data);
        void operator()(size_t new_size) const;
    };

    /** Add a particle data with default value. */
    template <typename DataType>
    struct addParticleDataWithDefaultValue
    {
        void operator()(ParticleData &particle_data) const;
    };

    template <typename DataType>
    struct copyParticleData
    {
        void operator()(ParticleData &particle_data, size_t index, size_t another_index) const;
    };

  public:
    //----------------------------------------------------------------------
    //		Assemble based generalize particle operations
    //----------------------------------------------------------------------
    DataAssembleOperation<resizeParticleData> resize_particle_data_;
    DataAssembleOperation<addParticleDataWithDefaultValue> add_particle_data_with_default_value_;
    DataAssembleOperation<copyParticleData> copy_particle_data_;
};

/**
 * @struct WriteAParticleVariableToXml
 * @brief Define a operator for writing particle variable to XML format.
 */
struct WriteAParticleVariableToXml
{
    XmlParser &xml_parser_;
    size_t &total_real_particles_;
    WriteAParticleVariableToXml(XmlParser &xml_parser, size_t &total_real_particles)
        : xml_parser_(xml_parser), total_real_particles_(total_real_particles){};

    template <typename DataType>
    void operator()(const std::string &variable_name, StdLargeVec<DataType> &variable) const;
};

/**
 * @struct ReadAParticleVariableFromXml
 * @brief Define a operator for reading particle variable to XML format.
 */
struct ReadAParticleVariableFromXml
{
    XmlParser &xml_parser_;
    size_t &total_real_particles_;
    ReadAParticleVariableFromXml(XmlParser &xml_parser, size_t &total_real_particles)
        : xml_parser_(xml_parser), total_real_particles_(total_real_particles){};

    template <typename DataType>
    void operator()(const std::string &variable_name, StdLargeVec<DataType> &variable) const;
};

/**
 * @class DiagnosticVariable
 * @brief Variable whose value is derived from other variables
 * and not used for evolution of the system.
 */
template <typename DataType>
class DiagnosticVariable
{
  public:
    using DiagnosticDataType = DataType;
    std::string variable_name_;

    DiagnosticVariable(SPHBody &sph_body, const std::string &variable_name);
    virtual ~DiagnosticVariable(){};

  protected:
    StdLargeVec<DataType> diagnostic_variable_;
};
} // namespace SPH
#endif // BASE_PARTICLES_H
