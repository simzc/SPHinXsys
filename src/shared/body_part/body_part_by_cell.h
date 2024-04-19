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
 * @file 	body_part_by_cell.h
 * @brief 	This is the base classes of body parts.
 * @details	There two main type of body parts. One is part by particle.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BODY_PART_BY_CELL_H
#define BODY_PART_BY_CELL_H

#include "body_part.h"

namespace SPH
{
/**
 * @class BodyPartByCell
 * @brief A body part with a collection of cell lists.
 */
class BodyPartByCell : public BodyPart
{
  public:
    ConcurrentCellLists body_part_cells_; /**< Collection of cells to indicate the body part. */
    ConcurrentCellLists &LoopRange() { return body_part_cells_; };
    size_t SizeOfLoopRange();

    BodyPartByCell(RealBody &real_body, const std::string &body_part_name)
        : BodyPart(real_body, body_part_name), cell_linked_list_(real_body.getCellLinkedList()){};
    virtual ~BodyPartByCell(){};

  protected:
    BaseCellLinkedList &cell_linked_list_;
    typedef std::function<bool(Vecd, Real)> TaggingCellMethod;
    void tagCells(TaggingCellMethod &tagging_cell_method);
};

/**
 * @class BodyRegionByCell
 * @brief A body part with the cell lists within a prescribed shape.
 */
class BodyRegionByCell : public BodyPartByCell
{
  private:
    SharedPtrKeeper<Shape> shape_ptr_keeper_;

  public:
    BodyRegionByCell(RealBody &real_body, Shape &body_part_shape);
    BodyRegionByCell(RealBody &real_body, SharedPtr<Shape> shape_ptr);
    virtual ~BodyRegionByCell(){};
    Shape &getBodyPartShape() { return body_part_shape_; };

  private:
    Shape &body_part_shape_;
    bool checkNotFar(Vecd cell_position, Real threshold);
};

/**
 * @class NearShapeSurface
 * @brief A body part with the cell lists near the surface of a prescribed shape.
 * @ details The body part shape can be that of the body,
 * or a sub shape of the body shape, or a shape independent of the body shape.
 * Note that only cells near the surface of the body part shape are included.
 */
class NearShapeSurface : public BodyPartByCell
{
  private:
    UniquePtrKeeper<LevelSetShape> level_set_shape_keeper_;

  public:
    NearShapeSurface(RealBody &real_body, SharedPtr<Shape> shape_ptr);
    NearShapeSurface(RealBody &real_body, LevelSetShape &level_set_shape);
    explicit NearShapeSurface(RealBody &real_body);
    NearShapeSurface(RealBody &real_body, const std::string &sub_shape_name);
    virtual ~NearShapeSurface(){};
    LevelSetShape &getLevelSetShape() { return level_set_shape_; };

  private:
    LevelSetShape &level_set_shape_;
    bool checkNearSurface(Vecd cell_position, Real threshold);
};

} // namespace SPH
#endif // BODY_PART_BY_CELL_H
