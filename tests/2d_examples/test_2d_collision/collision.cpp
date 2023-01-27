/**
 * @file 	collision.cpp
 * @brief 	two soft balls with and without internal damping bouncing within a confined boundary
 * @details This is the first case for test collision dynamics for
 * 			understanding SPH method for complex simulation.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 8.0;				  /**< box length. */
Real DH = 4.0;				  /**< box height. */
Real resolution_ref = 0.025;  /**< reference resolution. */
Real BW = resolution_ref * 4; /**< wall width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
Vec2d ball_center_1(2.0, 2.0);
Vec2d ball_center_2(6.0, 2.0);
Real ball_radius = 0.5;
// observer location
StdVec<Vecd> observation_location_1 = {ball_center_1};
StdVec<Vecd> observation_location_2 = {ball_center_2};
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real gravity_g = 1.0;
Real rho0_s = 1.0e3;
Real Youngs_modulus = 5.0e4;
Real poisson = 0.45;
Real physical_viscosity = 10000.0;
//----------------------------------------------------------------------
//	Geometric shapes
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		std::vector<Vecd> outer_wall_shape;
		outer_wall_shape.push_back(Vecd(-BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, -BW));

		std::vector<Vecd> inner_wall_shape;
		inner_wall_shape.push_back(Vecd(0.0, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, DH));
		inner_wall_shape.push_back(Vecd(DL, DH));
		inner_wall_shape.push_back(Vecd(DL, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, 0.0));

		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
class FreeBall : public MultiPolygonShape
{
public:
	explicit FreeBall(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addACircle(ball_center_1, ball_radius, 100, ShapeBooleanOps::add);
	}
};
class DampingBall : public MultiPolygonShape
{
public:
	explicit DampingBall(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addACircle(ball_center_2, ball_radius, 100, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** Tag for running particle relaxation for the initially body-fitted distribution */
	sph_system.setRunParticleRelaxation(true);
	/** Tag for starting with relaxed body-fitted particles distribution */
	sph_system.setReloadParticles(true);
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody free_ball(sph_system, makeShared<FreeBall>("FreeBall"));
	free_ball.defineBodyLevelSetShape();
	free_ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? free_ball.generateParticles<ParticleGeneratorReload>(io_environment, free_ball.getName())
		: free_ball.generateParticles<ParticleGeneratorLattice>();

	SolidBody damping_ball(sph_system, makeShared<DampingBall>("DampingBall"));
	damping_ball.defineBodyLevelSetShape();
	damping_ball.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? damping_ball.generateParticles<ParticleGeneratorReload>(io_environment, damping_ball.getName())
		: damping_ball.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	wall_boundary.generateParticles<ParticleGeneratorLattice>();

	ObserverBody free_ball_observer(sph_system, "FreeBallObserver");
	free_ball_observer.generateParticles<ObserverParticleGenerator>(observation_location_1);
	ObserverBody damping_ball_observer(sph_system, "DampingBallObserver");
	damping_ball_observer.generateParticles<ObserverParticleGenerator>(observation_location_2);
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Define body relation map used for particle relaxation.
		//----------------------------------------------------------------------
		InnerRelation free_ball_inner(free_ball);
		InnerRelation damping_ball_inner(damping_ball);
		//----------------------------------------------------------------------
		//	Define the methods for particle relaxation.
		//----------------------------------------------------------------------
		SimpleDynamics<RandomizeParticlePosition> free_ball_random_particles(free_ball);
		SimpleDynamics<RandomizeParticlePosition> damping_ball_random_particles(damping_ball);
		relax_dynamics::RelaxationStepInner free_ball_relaxation_step_inner(free_ball_inner);
		relax_dynamics::RelaxationStepInner damping_ball_relaxation_step_inner(damping_ball_inner);
		//----------------------------------------------------------------------
		//	Output for particle relaxation.
		//----------------------------------------------------------------------
		ReducedQuantityRecording<ReduceDynamics<QuantitySummation<Vecd>>> 
			free_ball_residue_force_recording(io_environment, free_ball, "Acceleration");
		ReducedQuantityRecording<ReduceDynamics<QuantitySummation<Vecd>>>
			damping_ball_residue_force_recording(io_environment, damping_ball, "Acceleration");
		BodyStatesRecordingToVtp write_ball_state(io_environment, sph_system.real_bodies_);
		ReloadParticleIO write_particle_reload_files(io_environment, {&free_ball, &damping_ball});
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		free_ball_random_particles.parallel_exec(0.25);
		damping_ball_random_particles.parallel_exec(0.25);
		sph_system.updateSystemCellLinkedLists();
		sph_system.updateSystemRelations();
		//----------------------------------------------------------------------
		//	First output before the simulation.
		//----------------------------------------------------------------------
		write_ball_state.writeToFileByStep();
		//----------------------------------------------------------------------
		//	From here iteration for particle relaxation begins.
		//----------------------------------------------------------------------
		int relax_step = 1000;
		while (sph_system.TotalSteps() < relax_step)
		{
			free_ball_relaxation_step_inner.exec();
			damping_ball_relaxation_step_inner.exec();
			sph_system.accumulateTotalSteps();

			free_ball_residue_force_recording.writeToFileByStep();
			damping_ball_residue_force_recording.writeToFileByStep();
			write_ball_state.writeToFileByStep();

			sph_system.monitorSteps("FreeBallResidueForce", free_ball_residue_force_recording.ResultValue().norm(),
									"DampingBallResidueForce", damping_ball_residue_force_recording.ResultValue().norm());

			sph_system.updateSystemCellLinkedLists();
			sph_system.updateSystemRelations();
		}
		std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
		write_particle_reload_files.writeToFileByStep();
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	TotalLagrangian<InnerRelation> free_ball_inner(free_ball);
	SurfaceContactRelation free_ball_contact(free_ball, {&wall_boundary});
	TotalLagrangian<InnerRelation> damping_ball_inner(damping_ball);
	SurfaceContactRelation damping_ball_contact(damping_ball, {&wall_boundary});

	TotalLagrangian<ObserverRelation> free_ball_observer_contact(free_ball_observer, RealBodyVector{&free_ball});
	TotalLagrangian<ObserverRelation> damping_all_observer_contact(damping_ball_observer, RealBodyVector{&damping_ball});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
	SimpleDynamics<TimeStepInitialization> free_ball_step_initialize(free_ball, gravity_ptr);
	SimpleDynamics<TimeStepInitialization> damping_ball_step_initialize(damping_ball, gravity_ptr);
	InteractionDynamics<solid_dynamics::CorrectConfiguration> free_ball_corrected_configuration(free_ball_inner);
	InteractionDynamics<solid_dynamics::CorrectConfiguration> damping_ball_corrected_configuration(damping_ball_inner);
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> free_ball_get_time_step_size(free_ball);
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> damping_ball_get_time_step_size(damping_ball);
	/** stress relaxation for the balls. */
	Dynamics1Level<solid_dynamics::Integration1stHalf> free_ball_stress_relaxation_first_half(free_ball_inner);
	Dynamics1Level<solid_dynamics::Integration2ndHalf> free_ball_stress_relaxation_second_half(free_ball_inner);
	Dynamics1Level<solid_dynamics::Integration1stHalf> damping_ball_stress_relaxation_first_half(damping_ball_inner);
	Dynamics1Level<solid_dynamics::Integration2ndHalf> damping_ball_stress_relaxation_second_half(damping_ball_inner);
	/** Algorithms for solid-solid contact. */
	InteractionDynamics<solid_dynamics::ContactDensitySummation, BodyPartByParticle> free_ball_update_contact_density(free_ball_contact);
	InteractionDynamics<solid_dynamics::ContactForceFromWall, BodyPartByParticle> free_ball_compute_solid_contact_forces(free_ball_contact);
	InteractionDynamics<solid_dynamics::ContactDensitySummation, BodyPartByParticle> damping_ball_update_contact_density(damping_ball_contact);
	InteractionDynamics<solid_dynamics::ContactForceFromWall, BodyPartByParticle> damping_ball_compute_solid_contact_forces(damping_ball_contact);
	/** Damping for one ball */
	DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
		damping(0.5, damping_ball_inner, "Velocity", physical_viscosity);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		free_ball_displacement_recording("Position", io_environment, free_ball_observer_contact);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		damping_ball_displacement_recording("Position", io_environment, damping_all_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.updateSystemCellLinkedLists();
	sph_system.updateSystemRelations();
	free_ball_corrected_configuration.parallel_exec();
	damping_ball_corrected_configuration.parallel_exec();
	//----------------------------------------------------------------------
	//	Initial states output.
	//----------------------------------------------------------------------
	body_states_recording.writeToFileByTime();
	free_ball_displacement_recording.writeToFileByStep();
	damping_ball_displacement_recording.writeToFileByStep();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	Real T0 = 10.0;
	Real end_time = T0;
	Real output_interval = 0.01 * T0;
	Real Dt = 0.1 * output_interval;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_interval)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_free = free_ball_get_time_step_size.parallel_exec();
				Real dt_damping = damping_ball_get_time_step_size.parallel_exec();
				Real dt = SMIN(dt_free, dt_damping);
				free_ball_step_initialize.parallel_exec();
				damping_ball_step_initialize.parallel_exec();

				free_ball_update_contact_density.parallel_exec();
				free_ball_compute_solid_contact_forces.parallel_exec();
				free_ball_stress_relaxation_first_half.parallel_exec(dt);
				free_ball_stress_relaxation_second_half.parallel_exec(dt);

				damping_ball_update_contact_density.parallel_exec();
				damping_ball_compute_solid_contact_forces.parallel_exec();
				damping_ball_stress_relaxation_first_half.parallel_exec(dt);
				damping.parallel_exec(dt);
				damping_ball_stress_relaxation_second_half.parallel_exec(dt);

				sph_system.accumulateTotalSteps();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				sph_system.monitorSteps("Time", GlobalStaticVariables::physical_time_, "elastic_dynamics_dt", dt);

				free_ball_displacement_recording.writeToFileByStep();
				damping_ball_displacement_recording.writeToFileByStep();

				/** Update cell linked list and configuration. */
				sph_system.updateSystemCellLinkedLists();
				sph_system.updateSystemRelations();
			}
		}
		tick_count t2 = tick_count::now();
		body_states_recording.writeToFileByTime();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	free_ball_displacement_recording.newResultTest();
	damping_ball_displacement_recording.newResultTest();

	return 0;
}
