/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: Lid Driven Cavity 2.5D                     *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// setup properties
Real particle_spacing = 0.02;
Real gravity_g = 0.0;
Real end_time = 30.0;

// material properties
Real rho = 1000.0;       // reference density
Real u_lid = 1.0;        // lid velocity
Real SOS = 10.0 * u_lid; // numerical speed of sound

// non-Newtonian properties
Real K = 1;     // consistency index
Real n = 0.5;   // power index
Real tau_y = 0; // yield stress

Real min_shear_rate = 1e-2; // cutoff low shear rate
Real max_shear_rate = 1e+3; // cutoff high shear rate

// mesh geometry data
std::string path_to_lid_boundary = "./input/lid_boundary.stl";
std::string path_to_no_slip_boundary = "./input/no_slip_boundary.stl";
std::string path_to_slip_boundary = "./input/slip_boundary.stl";
std::string path_to_fluid = "./input/fluid.stl";

//----------------------------------------------------------------------
//	Complex shapes for wall boundary
//----------------------------------------------------------------------
class Lid_Boundary : public ComplexShape
{
  public:
    explicit Lid_Boundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_lid_boundary, Vecd(0.0, 0.0, 0.0), 1.0);
    }
};
class No_Slip_Boundary : public ComplexShape
{
  public:
    explicit No_Slip_Boundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_no_slip_boundary, Vecd(0.0, 0.0, 0.0), 1.0);
    }
};
class Slip_Boundary : public ComplexShape
{
  public:
    explicit Slip_Boundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_slip_boundary, Vecd(0.0, 0.0, 0.0), 1.0);
    }
};
class FluidFilling : public ComplexShape
{
  public:
    explicit FluidFilling(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_fluid, Vecd(0.0, 0.0, 0.0), 1.0);
    }
};

class InitialVelocity
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    InitialVelocity(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          fluid_particles_(dynamic_cast<BaseParticles *>(&sph_body.getBaseParticles())){};

    void update(size_t index_i, Real dt)
    {
        /** initial velocity profile */
        vel_[index_i][1] = u_lid;
    }

  protected:
    BaseParticles *fluid_particles_;
};

int main(int ac, char *av[])
{
    //	Build up an SPHSystem
    BoundingBox system_domain_bounds(Vecd(-0.3, -0.7, -0.7), Vecd(0.3, 0.7, 0.7));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();

    //	Creating bodies with corresponding materials and particles
    FluidBody fluid(sph_system, makeShared<FluidFilling>("FluidBody"));
    fluid.defineParticlesAndMaterial<BaseParticles, HerschelBulkleyFluid>(rho, SOS, min_shear_rate, max_shear_rate, K, n, tau_y);
    fluid.generateParticles<ParticleGeneratorLattice>();

    SolidBody no_slip_boundary(sph_system, makeShared<No_Slip_Boundary>("NoSlipWall"));
    no_slip_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    no_slip_boundary.generateParticles<ParticleGeneratorLattice>();
    no_slip_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");

    SolidBody slip_boundary(sph_system, makeShared<Slip_Boundary>("SlipWall"));
    slip_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    slip_boundary.generateParticles<ParticleGeneratorLattice>();
    slip_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");

    SolidBody lid_boundary(sph_system, makeShared<Lid_Boundary>("LidWall"));
    lid_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    lid_boundary.generateParticles<ParticleGeneratorLattice>();
    lid_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");

    //	Define body relation map
    InnerRelation fluid_inner(fluid);
    ContactRelation fluid_all_walls(fluid, {&lid_boundary, &slip_boundary, &no_slip_boundary});
    ContactRelation fluid_no_slip(fluid, {&no_slip_boundary, &lid_boundary});

    ComplexRelation fluid_walls_complex(fluid_inner, fluid_all_walls);

    //	Define the numerical methods used in the simulation
    Gravity gravity(Vec3d(0.0, 0.0, gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(fluid, gravity);
    SimpleDynamics<InitialVelocity> initial_condition(lid_boundary);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(fluid_inner, fluid_all_walls);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(fluid_inner, fluid_all_walls);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(fluid_inner, fluid_all_walls);

    InteractionDynamics<fluid_dynamics::VelocityGradientWithWall> vel_grad_calculation(fluid_inner, fluid_no_slip);
    InteractionDynamics<fluid_dynamics::ShearRateDependentViscosity> shear_rate_calculation(fluid_inner);
    InteractionWithUpdate<fluid_dynamics::GeneralizedNewtonianViscousForceWithWall> viscous_acceleration(fluid_inner, fluid_no_slip);

    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_velocity_correction(fluid_inner, fluid_all_walls);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(fluid, u_lid);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_acoustic_time_step_size(fluid);
    ReduceDynamics<fluid_dynamics::SRDViscousTimeStepSize> get_viscous_time_step_size(fluid);

    //	Define the methods for I/O operations, observations
    fluid.addBodyStateForRecording<Real>("Pressure");
    BodyStatesRecordingToVtp write_fluid_states(sph_system.real_bodies_);

    //	Prepare the simulation
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    constant_gravity.exec();
    initial_condition.exec();

    //	Setup for time-stepping control
    // size_t number_of_iterations = sph_system.RestartStep();
    int nmbr_of_outputs = 100;
    Real output_interval = end_time / nmbr_of_outputs;
    Real dt = 0;
    Real Dt = 0;
    Real Dt_visc = 0;
    Real Dt_adv = 0;
    Real Dt_aco = 0;
    int iteration = 0;
    int output_counter = 1;

    //	First output before the main loop.
    write_fluid_states.writeToFile(0);
    TickCount t1 = TickCount::now();
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        iteration++;
        TimeInterval tt;
        TickCount t2 = TickCount::now();
        tt = t2 - t1;
        Dt_adv = get_fluid_advection_time_step_size.exec();
        Dt_visc = get_viscous_time_step_size.exec();
        Dt = SMIN(Dt_visc, Dt_adv);

        update_density_by_summation.exec(Dt);

        vel_grad_calculation.exec(Dt);
        shear_rate_calculation.exec(Dt);
        viscous_acceleration.exec(Dt);

        transport_velocity_correction.exec();

        Real relaxation_time = 0.0;
        while (relaxation_time < Dt)
        {
            Dt_aco = get_acoustic_time_step_size.exec();
            dt = SMIN(Dt_aco, Dt);
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);
            relaxation_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
        }

        if (iteration < 100 || output_counter * output_interval < GlobalStaticVariables::physical_time_)
        {
            std::cout << "Iteration: " << iteration << " | sim time in %: " << GlobalStaticVariables::physical_time_ / end_time * 100 << " | physical time in s: " << GlobalStaticVariables::physical_time_ << " | computation time in s: " << tt.seconds() << " | dt_adv: " << Dt_adv << " | dt_visc: " << Dt_visc << " | dt_aco: " << Dt_aco << "\r" << std::flush;
        }

        if (output_counter * output_interval < GlobalStaticVariables::physical_time_)
        {
            write_fluid_states.writeToFile();
            output_counter++;
        }
        fluid.updateCellLinkedListWithParticleSort(100);
        fluid_walls_complex.updateConfiguration();
        fluid_no_slip.updateConfiguration();
    }
    TickCount t3 = TickCount::now();
    TimeInterval te;
    te = t3 - t1;
    std::cout << "Done with iterations: " << iteration << " | Total computation time in s: " << (t3 - t1).seconds() << std::endl;
    return 0;
}