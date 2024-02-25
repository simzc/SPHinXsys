/**
 * @file fsi2.cpp
 * @brief This is the benchmark test of fluid-structure interaction.
 * @details We consider a flow-induced vibration of an elastic beam behind a cylinder in 2D.
 * The case can be found in Chi Zhang, Massoud Rezavand, Xiangyu Hu,
 * Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics.
 * Journal of Computation Physics 404 (2020) 109135.
 * @author Chi Zhang and Xiangyu Hu
 */
#include "fsi2.h" // case file to setup the test case
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);  // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(false);        // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments

    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<Lattice>();

    SolidBody insert_body(sph_system, makeShared<Insert>("InsertedBody"));
    insert_body.defineAdaptationRatios(1.15, 2.0);
    insert_body.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    insert_body.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? insert_body.generateParticles<Reload>(insert_body.getName())
        : insert_body.generateParticles<Lattice>();

    ObserverBody beam_observer(sph_system, "BeamObserver");
    StdVec<Vecd> beam_observation_location = {0.5 * (BRT + BRB)};
    beam_observer.generateParticles<Observer>(beam_observation_location);
    ObserverBody fluid_observer(sph_system, "FluidObserver");
    auto fluid_observer_particle_generator = fluid_observer.makeSelfDefined<ParticleGeneratorFluidObserver>();
    fluid_observer.generateParticles(fluid_observer_particle_generator);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation insert_body_inner(insert_body);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(insert_body);
        RelaxationStepInner relaxation_step_inner(insert_body_inner);
        BodyStatesRecordingToVtp write_insert_body_to_vtp({&insert_body});
        ReloadParticleIO write_particle_reload_files({&insert_body});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_insert_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_insert_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_insert_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation insert_body_inner(insert_body);
    ContactRelation water_block_contact(water_block, RealBodyVector{&wall_boundary, &insert_body});
    ContactRelation insert_body_contact(insert_body, {&water_block});
    ContactRelation beam_observer_contact(beam_observer, {&insert_body});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // and only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation using verlet time stepping. */
    /** Here, we do not use Riemann solver for pressure as the flow is viscous. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    /** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_correction(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    /** Computing vorticity in the flow for visualization. */
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    /** Inflow boundary condition. */
    BodyAlignedBoxByCell inflow_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);
    /** Periodic BCs in x direction. */
    PeriodicAlongAxis periodic_along_x(water_block.getSPHBodyBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition(water_block, periodic_along_x);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> insert_body_normal_direction(insert_body);
    /** Corrected configuration for the elastic insert body. */
    InteractionWithUpdate<KernelCorrectionMatrixInner> insert_body_corrected_configuration(insert_body_inner);
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(insert_body_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(insert_body_contact);
    /** Compute the average velocity of the insert body. */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(insert_body);
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    /** Compute time step size of elastic solid. */
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> insert_body_computing_time_step_size(insert_body);
    /** Stress relaxation for the inserted body. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> insert_body_stress_relaxation_first_half(insert_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> insert_body_stress_relaxation_second_half(insert_body_inner);
    /** Constrain region of the inserted body. */
    BodyRegionByParticle beam_base(insert_body, makeShared<MultiPolygonShape>(createBeamBaseShape()));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_beam_base(beam_base);
    /** Update surface norm direction.*/
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> insert_body_update_normal(insert_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system.real_bodies_);
    RegressionTestTimeAverage<ReducedQuantityRecording<QuantitySummation<Vecd>>> write_total_viscous_force_from_fluid(insert_body, "ViscousForceFromFluid");
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_beam_tip_displacement("Position", beam_observer_contact);
    ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    periodic_condition.update_cell_linked_list_.exec();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing surface normal direction for the wall. */
    wall_boundary_normal_direction.exec();
    /** computing surface normal direction for the insert body. */
    insert_body_normal_direction.exec();
    /** computing linear reproducing configuration for the insert body. */
    insert_body_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
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
    write_beam_tip_displacement.writeToFile(number_of_iterations);
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
            transport_correction.exec();

            /** FSI for viscous force. */
            viscous_force_from_fluid.exec();
            /** Update normal direction on elastic body.*/
            insert_body_update_normal.exec();
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                pressure_force_from_fluid.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(insert_body_computing_time_step_size.exec(), dt - dt_s_sum);
                    insert_body_stress_relaxation_first_half.exec(dt_s);
                    constraint_beam_base.exec();
                    insert_body_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                parabolic_inflow.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            periodic_condition.bounding_.exec();

            water_block.updateCellLinkedListWithParticleSort(100);
            periodic_condition.update_cell_linked_list_.exec();
            water_block_complex.updateConfiguration();
            /** one need update configuration after periodic condition. */
            insert_body.updateCellLinkedList();
            insert_body_contact.updateConfiguration();
            /** write run-time observation into file */
            write_beam_tip_displacement.writeToFile(number_of_iterations);
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
        fluid_observer_contact.updateConfiguration();
        write_fluid_velocity.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        // The lift force at the cylinder is very small and not important in this case.
        write_total_viscous_force_from_fluid.generateDataBase({1.0e-2, 1.0e-2}, {1.0e-2, 1.0e-2});
        write_beam_tip_displacement.generateDataBase(1.0e-2);
    }
    else
    {
        write_total_viscous_force_from_fluid.testResult();
        write_beam_tip_displacement.testResult();
    }

    return 0;
}
