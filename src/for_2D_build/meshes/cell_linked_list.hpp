/**
 * @file 	cell_linked_list.hpp
 * @brief 	Here gives the classes for managing cell linked lists. This is the basic class
 * 			for building the particle configurations.
 * @details The cell linked list saves for each body a list of particles
 * 			located within the cell.
 * @author	Chi ZHang, Yongchuan and Xiangyu Hu
 */

#pragma once

#include "base_particles.h"
#include "cell_linked_list.h"
#include "mesh_iterators.hpp"
#include "particle_iterators.h"

namespace SPH
{
//=================================================================================================//
template <class DynamicsRange, typename GetSearchDepth, typename GetNeighborRelation>
void CellLinkedList::searchNeighborsByParticles(
    DynamicsRange &dynamics_range, ParticleConfiguration &particle_configuration,
    GetSearchDepth &get_search_depth, GetNeighborRelation &get_neighbor_relation)
{
    StdLargeVec<Vecd> &pos = dynamics_range.getBaseParticles().pos_;
    particle_for(execution::ParallelPolicy(), dynamics_range.LoopRange(),
                 [&](size_t index_i)
                 {
                     int search_depth = get_search_depth(index_i);
                     Array2i target_cell_index = CellIndexFromPosition(pos[index_i]);

                     Neighborhood &neighborhood = particle_configuration[index_i];
                     mesh_for_each(
                         Array2i::Zero().max(target_cell_index - search_depth * Array2i::Ones()),
                         all_cells_.min(target_cell_index + (search_depth + 1) * Array2i::Ones()),
                         [&](int l, int m)
                         {
                             ListDataVector &target_particles = cell_data_lists_[l][m];
                             for (const ListData &list_data : target_particles)
                             {
                                 get_neighbor_relation(neighborhood, pos[index_i], index_i, list_data);
                             }
                         });
                 });
}
//=================================================================================================//
template <class DynamicsRange, typename GetSearchDepth, typename GetNeighborRelation>
void CellLinkedListKernel::searchNeighborsByParticles(
    DynamicsRange &dynamics_range, NeighborhoodDevice *particle_configuration,
    GetSearchDepth &get_search_depth, GetNeighborRelation& get_neighbor_relation)
{
    auto *pos = dynamics_range.getBaseParticles().template getDeviceVariableByName<DeviceVec2d>("Position");
    auto get_neighbor_relation_kernel = *get_neighbor_relation.getDeviceProxy().getKernel();
    const int search_depth = get_search_depth(0);
    executionQueue.getQueue()
        .submit(
            [&, index_list=index_list_, index_head_list=index_head_list_, list_data_pos=list_data_pos_,
             list_data_Vol=list_data_Vol_, kernel=kernel_, mesh_lower_bound=mesh_lower_bound_,
             grid_spacing=grid_spacing_, all_grid_points=all_grid_points_, all_cells=all_cells_](sycl::handler &cgh) {
                cgh.parallel_for(dynamics_range.LoopRange(), [=](sycl::item<1> it) {
                                     const size_t index_i = it.get_id(0);
                                     const auto& pos_i = pos[index_i];
                                     auto &neighborhood = particle_configuration[index_i];
                                     const auto target_cell_index = CellIndexFromPosition(pos_i, mesh_lower_bound,
                                                                             grid_spacing, all_grid_points);
                                     mesh_for_each(
                                         sycl::max(DeviceArray2i{0}, target_cell_index - search_depth),
                                         sycl::min(all_cells, target_cell_index + search_depth + 1),
                                         [&](int l, int m) {
                                             const auto linear_cell_index = transferCellIndexTo1D({l,m}, all_cells);
                                             size_t index_j = index_head_list[linear_cell_index];
                                             // Cell list ends when index_j == 0, if index_j is already zero then cell is empty.
                                             while(index_j--) {  // abbreviates while(index_j != 0) { index_j -= 1; ... }
                                                 get_neighbor_relation_kernel(neighborhood, pos_i, index_i, index_j, list_data_pos[index_j],
                                                                              list_data_Vol[index_j], kernel);
                                                 index_j = index_list[index_j];
                                             }
                                         });
                                 });
            }).wait_and_throw();
}
//=================================================================================================//
} // namespace SPH
