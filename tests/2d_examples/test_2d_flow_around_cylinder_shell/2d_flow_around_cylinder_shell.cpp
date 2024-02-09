/**
 * @file 	2d_flow_around_cylinder.cpp
 * @brief 	This is the benchmark test for the wall modeling of viscous flow.
 * @details We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu
 */
#include "2d_flow_around_cylinder_shell.h"
#include "sphinxsys.h"

using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
// handle command line arguments
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(sph_system);
    ParameterizationIO &parameterization_io = io_environment.defineParameterizationIO();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineParticlesAndMaterial<BaseParticles, ParameterizedWaterMaterial>(parameterization_io, rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody cylinder(sph_system, makeShared<DefaultShape>("Cylinder"));
    cylinder.defineAdaptationRatios(1.15, resolution_ref / resolution_shell);
    cylinder.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(1, 1, 0);
    cylinder.generateParticles<CylinderParticleGenerator>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ParticleGeneratorObserver>(observation_locations);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelationToShell water_block_contact(water_block, {&cylinder}, {true});
    ContactRelationFromShell cylinder_contact(cylinder, {&water_block}, {true});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    ShellInnerRelationWithContactKernel cylinder_curvature_inner(cylinder, water_block);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which are only use for now for updating configurations.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> cylinder_average_curvature(cylinder_curvature_inner);
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_block, water_block.getBodyShapeBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_y(water_block, water_block.getBodyShapeBounds(), yAxis);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_velocity_correction(water_block_inner, water_block_contact);
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    BodyRegionByCell free_stream_buffer(water_block, makeShared<MultiPolygonShape>(createBufferShape()));
    SimpleDynamics<FreeStreamCondition> freestream_condition(free_stream_buffer);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_on_cylinder(cylinder_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid> pressure_force_on_cylinder(cylinder_contact);
    /** Computing viscous force acting on wall with wall model. */
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    cylinder.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    cylinder.addBodyStateForRecording<Vecd>("ViscousForceFromFluid");
    cylinder.addBodyStateForRecording<Vecd>("PressureForceFromFluid");
    BodyStatesRecordingToVtp write_real_body_states(sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<QuantitySummation<Vecd>>> write_total_viscous_force_from_fluid(cylinder, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_pressure_force_from_fluid(cylinder, "PressureForceFromFluid");
    ObservedQuantityRecording<Vecd>
        write_fluid_velocity("Velocity", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    periodic_condition_x.update_cell_linked_list_.exec();
    periodic_condition_y.update_cell_linked_list_.exec();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** initialize surface normal direction for the insert body. */
    cylinder_average_curvature.exec();
    water_block_complex.updateConfiguration();
    cylinder_contact.updateConfiguration();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 200.0;
    Real output_interval = end_time / 200.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;

        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();

            /** FSI for viscous force. */
            viscous_force_on_cylinder.exec();
            size_t inner_ite_dt = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation, first half. */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                pressure_force_on_cylinder.exec();
                /** Fluid pressure relaxation, second half. */
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                freestream_condition.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            periodic_condition_x.bounding_.exec();
            periodic_condition_y.bounding_.exec();
            water_block.updateCellLinkedListWithParticleSort(100);
            periodic_condition_x.update_cell_linked_list_.exec();
            periodic_condition_y.update_cell_linked_list_.exec();
            /** one need update configuration after periodic condition. */
            water_block_complex.updateConfiguration();
            cylinder_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
        write_total_pressure_force_from_fluid.writeToFile(number_of_iterations);
        fluid_observer_contact.updateConfiguration();
        write_fluid_velocity.writeToFile(number_of_iterations);

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    write_total_viscous_force_from_fluid.testResult();

    return 0;
}