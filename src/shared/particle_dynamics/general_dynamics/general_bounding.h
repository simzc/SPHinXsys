/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	general_bounding.h
 * @brief 	This is the particle dynamics for domain bounding
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GENERAL_BOUNDING_H
#define GENERAL_BOUNDING_H

#include "general_dynamics.h"

namespace SPH
{

/**
 * @class BoundingAlongAxis
 * @brief Bounding particle position in along axis.
 * The axis must be 0, 1 for 2d and 0, 1, 2 for 3d
 */
class BoundingAlongAxis : public LocalDynamics,
                          public GeneralDataDelegateSimple
{
  protected:
    const int axis_;              /**< the axis directions for bounding*/
    BoundingBox bounding_bounds_; /**< lower and upper bound for checking. */
    StdLargeVec<Vecd> &pos_;
    BaseCellLinkedList &cell_linked_list_;
    Real cut_off_radius_max_; /**< maximum cut off radius to avoid boundary particle depletion */

  public:
    BoundingAlongAxis(RealBody &real_body, BoundingBox bounding_bounds, int axis);
    virtual ~BoundingAlongAxis(){};

    template <class OperationType>
    void checkNearLowerBound(size_t index_i, const OperationType &operation)
    {
        if (pos_[index_i][axis_] > bounding_bounds_.first_[axis_] &&
            pos_[index_i][axis_] < (bounding_bounds_.first_[axis_] + cut_off_radius_max_))
        {
            operation(index_i);
        }
    };

    template <class OperationType>
    void checkNearUpperBound(size_t index_i, const OperationType &operation)
    {
        if (pos_[index_i][axis_] < bounding_bounds_.second_[axis_] &&
            pos_[index_i][axis_] > (bounding_bounds_.second_[axis_] - cut_off_radius_max_))
        {
            operation(index_i);
        }
    };
};

class PeriodicBoundary : public BoundingAlongAxis
{
  public:
    PeriodicBoundary(RealBody &real_body, BoundingBox bounding_bounds, int axis)
        : BoundingAlongAxis(real_body, bounding_bounds, axis),
          translation_(bounding_bounds.second_[axis] - bounding_bounds.first_[axis]){};

    void LowerBoundTranslation(size_t index_i)
    {
        pos_[index_i][axis_] += translation_;
    };

    void UpperBoundTranslation(size_t index_i)
    {
        pos_[index_i][axis_] -= translation_;
    };

    Vecd LowerBoundTranslation(Vecd pos)
    {
        pos[axis_] += translation_;
        return pos;
    };

    Vecd UpperBoundTranslation(Vecd &pos)
    {
        pos[axis_] -= translation_;
        return pos;
    };

  private:
    Real translation_;
};

class LeeEdwardsBoundary : public PeriodicBoundary
{
  public:
    LeeEdwardsBoundary(RealBody &real_body, const BoundingBox &bounding_bounds, int axis,
                       int shear_direction, Real shear_rate)
        : PeriodicBoundary(real_body, bounding_bounds, axis),
          vel_(particles_->vel_), shear_direction_(shear_direction),
          shear_center_(0.5 * (bounding_bounds.second_[shear_direction] + bounding_bounds.first_[shear_direction])),
          pos_increment_(bounding_bounds.second_[shear_direction] - bounding_bounds.first_[shear_direction]),
          vel_increment_(pos_increment_ * shear_rate){};

    void LowerBoundTranslation(size_t index_i)
    {
        PeriodicBoundary::LowerBoundTranslation(index_i);
        flipState(index_i);
    };

    void UpperBoundTranslation(size_t index_i)
    {
        PeriodicBoundary::UpperBoundTranslation(index_i);
        flipState(index_i);
    };

  protected:
    StdLargeVec<Vecd> &vel_;

  private:
    int shear_direction_; // upper bound moving toward positive side
    Real shear_center_;
    Real pos_increment_;
    Real vel_increment_;

    void flipState(size_t index_i)
    {
        bool is_positive_side = pos_[index_i][shear_direction_] > shear_center_;
        pos_[index_i][shear_direction_] += is_positive_side ? -pos_increment_ : pos_increment_;
        vel_[index_i][shear_direction_] += is_positive_side ? -vel_increment_ : vel_increment_;
    };
};

/**
 * @class BasePeriodicCondition
 * @brief Base class for two different type periodic boundary conditions.
 */
template <class PeriodicBoundaryType, class ExecutionPolicy>
class BasePeriodicCondition
{
  protected:
    StdVec<CellLists> bound_cells_data_;
    /**
     * @class PeriodicBounding
     * @brief Periodic bounding particle position in an axis direction
     */
    class PeriodicBounding : public BaseDynamics<void>, public PeriodicBoundaryType
    {
      protected:
        StdVec<CellLists> &bound_cells_data_;

        void checkLowerBound(size_t index_i)
        {
            if (this->pos[index_i][axis_] < this->bounding_bounds_.first_[axis_])
            {
                this->LowerBoundTranslation(index_i);
            }
        };

        void checkUpperBound(size_t index_i)
        {
            if (this->pos_[index_i][axis_] > this->bounding_bounds_.second_[axis_])
            {
                this->UpperBoundTranslation(index_i);
            }
        };

      public:
        template <typename... Args>
        PeriodicBounding(StdVec<CellLists> &bound_cells_data,
                         RealBody &real_body, BoundingBox bounding_bounds, int axis, Args &&...args)
            : BaseDynamics<void>(real_body),
              PeriodicBoundaryType(real_body, bounding_bounds, axis, std::forward<Args>(args)...),
              bound_cells_data_(bound_cells_data){};
        virtual ~PeriodicBounding(){};

        virtual void exec(Real dt = 0.0) override
        {
            particle_for(ExecutionPolicy(), bound_cells_data_[0].first,
                         [&](size_t i)
                         { checkLowerBound(i) });

            particle_for(ExecutionPolicy(), bound_cells_data_[1].first,
                         [&](size_t i)
                         { checkLowerBound(i) });
        };
    };

    /**
     * @class PeriodicImage
     * @brief Creating periodic image for boundary condition in an axis direction
     */
    class PeriodicImage : public BaseDynamics<void>, public PeriodicBoundaryType
    {
      protected:
        std::mutex mutex_cell_list_entry_; /**< mutex exclusion for memory conflict */
        StdVec<CellLists> &bound_cells_data_;
        virtual void checkLowerBound(ListDataVector &cell_list_data, Real dt = 0.0) = 0;
        virtual void checkUpperBound(ListDataVector &cell_list_data, Real dt = 0.0) = 0;

      public:
        template <typename... Args>
        PeriodicImage(StdVec<CellLists> &bound_cells_data,
                      RealBody &real_body, BoundingBox bounding_bounds, int axis, Args &&...args)
            : BaseDynamics<void>(real_body),
              PeriodicBoundaryType(real_body, bounding_bounds, axis, std::forward<Args>(args)...),
              bound_cells_data_(bound_cells_data){};
        virtual ~PeriodicImage(){};

        virtual void exec(Real dt = 0.0) override
        {
            this->setupDynamics(dt);

            particle_for(ExecutionPolicy(), bound_cells_data_[0].second,
                         [&](ListDataVector *cell_list)
                         { checkLowerBound(*cell_list, dt); });

            particle_for(ExecutionPolicy(), bound_cells_data_[1].second,
                         [&](ListDataVector *cell_list)
                         { checkUpperBound(*cell_list, dt); });
        };
    };

  public:
    BasePeriodicCondition(RealBody &real_body, BoundingBox bounding_bounds, int axis)
    {
        bound_cells_data_.resize(2);
        BaseCellLinkedList &cell_linked_list = real_body.getCellLinkedList();
        cell_linked_list.tagBoundingCells(bound_cells_data_, bounding_bounds, axis);
    };
    virtual ~BasePeriodicCondition(){};
};

/**
 * @class PeriodicConditionUsingCellLinkedList
 * @brief The method imposing periodic boundary condition in an axis direction.
 *	It includes two different steps, i.e. imposing periodic bounding and condition.
 *	The first step is carried out before update cell linked list and
 *	the second after the updating.
 *	If the exec or parallel_exec is called directly, error message will be given.
 */
class PeriodicConditionUsingCellLinkedList
    : public BasePeriodicCondition<PeriodicBoundary, execution::ParallelPolicy>
{
  protected:
    /**
     * @class PeriodicCellLinkedList
     * @brief Periodic boundary condition in an axis direction using cell linked list image
     */
    class PeriodicCellLinkedList : public PeriodicImage
    {
      protected:
        virtual void checkLowerBound(ListDataVector &cell_list_data, Real dt = 0.0) override;
        virtual void checkUpperBound(ListDataVector &cell_list_data, Real dt = 0.0) override;

      public:
        PeriodicCellLinkedList(StdVec<CellLists> &bound_cells_data,
                               RealBody &real_body, BoundingBox bounding_bounds, int axis)
            : PeriodicImage(bound_cells_data, real_body, bounding_bounds, axis){};
        virtual ~PeriodicCellLinkedList(){};
    };

  public:
    PeriodicConditionUsingCellLinkedList(RealBody &real_body, BoundingBox bounding_bounds, int axis)
        : BasePeriodicCondition<PeriodicBoundary, execution::ParallelPolicy>(real_body, bounding_bounds, axis),
          bounding_(bound_cells_data_, real_body, bounding_bounds, axis),
          update_cell_linked_list_(bound_cells_data_, real_body, bounding_bounds, axis)
    {
        real_body.addBeforeUpdateCellLinkedList(&bounding_);
        real_body.addAfterUpdateCellLinkedList(&update_cell_linked_list_);
    };
    virtual ~PeriodicConditionUsingCellLinkedList(){};

    PeriodicBounding bounding_;
    PeriodicCellLinkedList update_cell_linked_list_;
};

/**
 * @class PeriodicConditionUsingGhostParticles
 * @brief The method imposing periodic boundary condition in an axis direction by using ghost particles.
 *	It includes three different steps, i.e. imposing periodic bounding, creating ghosts and update ghost state.
 *	The first step is carried out before update cell linked list and
 *	the second and third after the updating.
 *	If the exec or parallel_exec is called directly, error message will be given.
 *  Note that, currently, this class is not for periodic condition in combined directions,
 *  such as periodic condition in both x and y directions.
 */
template <class PeriodicBoundaryType>
class PeriodicConditionUsingGhostParticles
    : public BasePeriodicCondition<PeriodicBoundaryType, execution::ParallelPolicy>
{
  protected:
    StdVec<IndexVector> ghost_particles_;

    /**
     * @class PeriodicGhost
     * @brief create ghost particles in an axis direction
     */
    class PeriodicGhost : public PeriodicImage
    {
      protected:
        StdVec<IndexVector> &ghost_particles_;
        virtual void setupDynamics(Real dt = 0.0) override
        {
            for (size_t i = 0; i != ghost_particles_.size(); ++i)
                ghost_particles_[i].clear();
        };
        virtual void checkLowerBound(ListDataVector &cell_list_data, Real dt = 0.0) override
        {
            for (size_t num = 0; num < cell_list_data.size(); ++num)
            {
                size_t index_i = std::get<0>(cell_list_data[num]);
                Vecd particle_position = std::get<1>(cell_list_data[num]);
                if (particle_position[axis_] > bounding_bounds_.first_[axis_] &&
                    particle_position[axis_] < (bounding_bounds_.first_[axis_] + cut_off_radius_max_))
                {
                    mutex_cell_list_entry_.lock();
                    size_t ghost_particle_index = this->particles_->insertAGhostParticle(index_i);
                    this->LowerBoundTranslation(size_t index_i);
                    /** insert ghost particle to cell linked list */
                    cell_linked_list_.InsertListDataEntry(ghost_particle_index,
                                                          this->pos_[ghost_particle_index],
                                                          std::get<2>(cell_list_data[num]));
                    mutex_cell_list_entry_.unlock();
                }
            }
        };
        virtual void checkUpperBound(ListDataVector &cell_list_data, Real dt = 0.0) override
        {
            for (size_t num = 0; num < cell_list_data.size(); ++num)
            {
                size_t index_i = std::get<0>(cell_list_data[num]);
                Vecd particle_position = std::get<1>(cell_list_data[num]);
                if (particle_position[axis_] < bounding_bounds_.second_[axis_] &&
                    particle_position[axis_] > (bounding_bounds_.second_[axis_] - cut_off_radius_max_))
                {
                    mutex_cell_list_entry_.lock();
                    size_t ghost_particle_index = particles_->insertAGhostParticle(index_i);
                    this->LowerBoundTranslation(size_t index_i);
                    /** insert ghost particle to cell linked list */
                    cell_linked_list_.InsertListDataEntry(ghost_particle_index,
                                                          this->pos_[ghost_particle_index],
                                                          std::get<2>(cell_list_data[num]));
                    mutex_cell_list_entry_.unlock();
                }
            }
        };

      public:
        template <typename... Args>
        PeriodicGhost(StdVec<IndexVector> &ghost_particles, StdVec<CellLists> &bound_cells_data,
                      RealBody &real_body, BoundingBox bounding_bounds, int axis, Args &&...args)
            : PeriodicImage(bound_cells_data, real_body, bounding_bounds, axis, std::forward<Args>(args)...),
              ghost_particles_(ghost_particles){};
        virtual ~PeriodicGhost(){};
    };

    /**
     * @class UpdatePeriodicGhost
     * @brief update ghost particles in an axis direction
     */
    class UpdatePeriodicGhost : public BaseDynamics<void>, public PeriodicBoundaryType
    {
      protected:
        StdVec<IndexVector> &ghost_particles_;
        void checkLowerBound(size_t index_i)
        {
            this->particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
            this->LowerBoundTranslation(index_i);
        };
        void checkUpperBound(size_t index_i)
        {
            this->particles_->updateFromAnotherParticle(index_i, sorted_id_[index_i]);
            this->UpperBoundTranslation(index_i);
        };

      public:
        template <typename... Args>
        UpdatePeriodicGhost(StdVec<IndexVector> &ghost_particles, RealBody &real_body,
                            BoundingBox bounding_bounds, int axis, Args &&...args)
            : BaseDynamics<void>(real_body),
              PeriodicBoundaryType(real_body, bounding_bounds, axis, std::forward<Args>(args)...),
              ghost_particles_(ghost_particles){};
        virtual ~UpdatePeriodicGhost(){};

        virtual void exec(Real dt = 0.0) override
        {
            particle_for(execution::ParallelPolicy(), ghost_particles_[0],
                         [&](size_t i)
                         { checkLowerBound(i); });

            particle_for(execution::ParallelPolicy(), ghost_particles_[1],
                         [&](size_t i)
                         { checkUpperBound(i); });
        };

      public:
        template <typename... Args>
        PeriodicConditionUsingGhostParticles(Args &&...args)
            : BasePeriodicCondition<PeriodicBoundaryType, execution::ParallelPolicy>(std::forward<Args>(args)...),
              bounding_(this->bound_cells_data_, std::forward<Args>(args)...),
              ghost_creation_(this->ghost_particles_, this->bound_cells_data_, std::forward<Args>(args)...),
              ghost_update_(this->ghost_particles_, std::forward<Args>(args)...)
        {
            ghost_particles_.resize(2);
            real_body.addBeforeUpdateCellLinkedList(&bounding_);
            real_body.addAfterUpdateCellLinkedList(&ghost_creation_);
        };

        virtual ~PeriodicConditionUsingGhostParticles(){};

        PeriodicBounding bounding_;
        PeriodicGhost ghost_creation_;
        UpdatePeriodicGhost ghost_update_;
    };
};
} // namespace SPH
#endif // GENERAL_BOUNDING_H
