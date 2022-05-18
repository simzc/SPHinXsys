/**
 * @file 	geometric_shape.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "geometric_shape.h"

namespace SPH
{
    //=================================================================================================//
    GeometricShapeRectangle::GeometricShapeRectangle(const Vec2d &halfsize, const std::string &shape_name)
        : Shape(shape_name), halfsize_(halfsize)
    {
        if (halfsize[0] < 0.0 || halfsize[1] < 0.0)
        {
            std::cout << "\n Error: the GeometricShapeRectangle half size must be positive! " << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
    }
    //=================================================================================================//
    bool GeometricShapeRectangle::checkContain(const Vec2d &pnt, bool BOUNDARY_INCLUDED)
    {
        return ABS(pnt[0]) < halfsize_[0] && ABS(pnt[1]) < halfsize_[1];
    }
    //=================================================================================================//
    Vec2d GeometricShapeRectangle::findClosestPoint(const Vec2d &pnt)
    {
        Real x = pnt[0] + halfsize_[0] > halfsize_[0] - pnt[0] ? halfsize_[0] : -halfsize_[0];
        Real y = pnt[1] + halfsize_[1] > halfsize_[1] - pnt[1] ? halfsize_[1] : -halfsize_[1];
        return Vec2d(x, y);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeRectangle::findBounds()
    {
        return BoundingBox(-halfsize_, halfsize_);
    }
    //=================================================================================================//
    GeometricShapeCircle::GeometricShapeCircle(const Vec2d &center, Real radius,
                                               const std::string &shape_name)
        : Shape(shape_name), center_(center), radius_(radius) {}
    //=================================================================================================//
    bool GeometricShapeCircle::checkContain(const Vec2d &pnt, bool BOUNDARY_INCLUDED)
    {
        return (pnt - center_).norm() < radius_;
    }
    //=================================================================================================//
    Vec2d GeometricShapeCircle::findClosestPoint(const Vec2d &pnt)
    {
        Vec2d displacement = pnt - center_;
        Real distance = displacement.norm();
        Real cosine = (SGN(displacement[0]) * (ABS(displacement[0])) + TinyReal) / (distance + TinyReal);
        Real sine = displacement[1] / (distance + TinyReal);
        return pnt + (radius_ - distance) * Vec2d(cosine, sine);
    }
    //=================================================================================================//
    BoundingBox GeometricShapeCircle::findBounds()
    {
        Vec2d shift = Vec2d(radius_, radius_);
        return BoundingBox(center_ - shift, center_ + shift);
    }
    //=================================================================================================//
    Real GeometricShapeCircle::findSignedDistance(const Vecd &pnt)
    {
        return (pnt - center_).norm() - radius_;
    }
    //=================================================================================================//
}