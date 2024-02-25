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
 * @file    particle_functors.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef PARTICLE_FUNCTORS_H
#define PARTICLE_FUNCTORS_H

#include "base_particles.h"

namespace SPH
{
//----------------------------------------------------------------------
// Particle group scope functors
//----------------------------------------------------------------------
class AllParticles
{
  public:
    explicit AllParticles(BaseParticles *base_particles){};
    bool operator()(size_t index_i)
    {
        return true;
    };
};

template <int INDICATOR>
class IndicatedParticles
{
    StdLargeVec<int> &indicator_;

  public:
    explicit IndicatedParticles(BaseParticles *base_particles)
        : indicator_(*base_particles->getVariableByName<int>("Indicator")){};
    bool operator()(size_t index_i)
    {
        return indicator_[index_i] == INDICATOR;
    };
};

using BulkParticles = IndicatedParticles<0>;

template <int INDICATOR>
class NotIndicatedParticles
{
    StdLargeVec<int> &indicator_;

  public:
    explicit NotIndicatedParticles(BaseParticles *base_particles)
        : indicator_(*base_particles->getVariableByName<int>("Indicator")){};
    bool operator()(size_t index_i)
    {
        return indicator_[index_i] != INDICATOR;
    };
};
//----------------------------------------------------------------------
// Sorted and unsorted particle functors
//----------------------------------------------------------------------
class UnSortedIndexVector : public IndexVector
{
    StdLargeVec<size_t> &sorted_id_;

  public:
    UnSortedIndexVector(BaseParticles *base_particles)
        : IndexVector(), sorted_id_(base_particles->sorted_id_){};

    size_t SortedID(size_t i) const
    {
        return sorted_id_[(*this)[i]];
    };
};
//----------------------------------------------------------------------
// Particle average functors
//----------------------------------------------------------------------
template <typename T>
class PairAverageFixed
{
    const T average_;

  public:
    PairAverageFixed(const T &c1, const T &c2)
        : average_(0.5 * (c1 + c2)){};
    explicit PairAverageFixed(const T &c)
        : PairAverageFixed(c, c){};
    T operator()(size_t index_i, size_t index_j)
    {
        return average_;
    };
};

class GeomAverage
{
  public:
    GeomAverage(){};

  protected:
    Real inverse(const Real &x) { return 1.0 / x; };

    template <typename MatrixType>
    MatrixType inverse(const MatrixType &x) { return x.inverse(); };
};

template <typename T>
class PairGeomAverageFixed : public GeomAverage
{
    const T geom_average_;

  public:
    PairGeomAverageFixed(const T &c1, const T &c2)
        : GeomAverage(), geom_average_(2.0 * c1 * c2 * inverse(c1 + c2)){};
    explicit PairGeomAverageFixed(const T &c)
        : PairGeomAverageFixed(c, c){};
    T operator()(size_t index_i, size_t index_j)
    {
        return geom_average_;
    };
};

template <typename T>
class PairAverageVariable
{
    StdLargeVec<T> &v1_, &v2_;

  public:
    PairAverageVariable(StdLargeVec<T> &v1, StdLargeVec<T> &v2)
        : v1_(v1), v2_(v2){};
    explicit PairAverageVariable(StdLargeVec<T> &v)
        : PairAverageVariable(v, v){};
    T operator()(size_t index_i, size_t index_j)
    {
        return 0.5 * (v1_[index_i] + v2_[index_j]);
    };
};

template <typename T>
class PairGeomAverageVariable : public GeomAverage
{
    StdLargeVec<T> &v1_, &v2_;

  public:
    PairGeomAverageVariable(StdLargeVec<T> &v1, StdLargeVec<T> &v2)
        : GeomAverage(), v1_(v1), v2_(v2){};
    explicit PairGeomAverageVariable(StdLargeVec<T> &v)
        : PairGeomAverageVariable(v, v){};

    T operator()(size_t index_i, size_t index_j)
    {
        return 2.0 * v1_[index_i] * v2_[index_j] * inverse(v1_[index_i] + v2_[index_j]);
    };
};
//----------------------------------------------------------------------
// Particle kernel functors
//----------------------------------------------------------------------
class NoKernelCorrection
{
  public:
    NoKernelCorrection(BaseParticles *particles){};
    Real operator()(size_t index_i)
    {
        return 1.0;
    };
};

class KernelCorrection
{
  public:
    KernelCorrection(BaseParticles *particles)
        : B_(*particles->getVariableByName<Matd>("KernelCorrectionMatrix")){};

    Matd operator()(size_t index_i)
    {
        return B_[index_i];
    };

  protected:
    StdLargeVec<Matd> &B_;
};

class SingleResolution
{
  public:
    SingleResolution(BaseParticles *particles){};
    Real operator()(size_t index_i)
    {
        return 1.0;
    };
};
//----------------------------------------------------------------------
// Particle adaptation functors
//----------------------------------------------------------------------
class AdaptiveResolution
{
  public:
    AdaptiveResolution(BaseParticles *particles)
        : h_ratio_(*particles->getVariableByName<Real>("SmoothingLengthRatio")){};

    Real operator()(size_t index_i)
    {
        return h_ratio_[index_i];
    };

  protected:
    StdLargeVec<Real> &h_ratio_;
};
//----------------------------------------------------------------------
// Particle reduce functors
//----------------------------------------------------------------------
template <class ReturnType>
struct ReduceSum
{
    ReturnType operator()(const ReturnType &x, const ReturnType &y) const { return x + y; };
};

struct ReduceMax
{
    Real operator()(Real x, Real y) const { return SMAX(x, y); };
};

struct ReduceMin
{
    Real operator()(Real x, Real y) const { return SMIN(x, y); };
};

struct ReduceOR
{
    bool operator()(bool x, bool y) const { return x || y; };
};

struct ReduceAND
{
    bool operator()(bool x, bool y) const { return x && y; };
};

struct ReduceLowerBound
{
    Vecd operator()(const Vecd &x, const Vecd &y) const
    {
        Vecd lower_bound;
        for (int i = 0; i < lower_bound.size(); ++i)
            lower_bound[i] = SMIN(x[i], y[i]);
        return lower_bound;
    };
};
struct ReduceUpperBound
{
    Vecd operator()(const Vecd &x, const Vecd &y) const
    {
        Vecd upper_bound;
        for (int i = 0; i < upper_bound.size(); ++i)
            upper_bound[i] = SMAX(x[i], y[i]);
        return upper_bound;
    };
};
} // namespace SPH
#endif // PARTICLE_FUNCTORS_H
