 /**
 * @file 	droplet_impact.cpp
 * @brief 	Impact of 2d water droplet in air onto a solid surface
 * @details This is the one of the basic test cases for understanding SPH method for multi-phase simulation.
 * @author 	Nilotpal Chakraborty
 */
#include "droplet_impact.h"
#include "sphinxsys.h"
using namespace SPH;

int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	FluidBody air_block(sph_system, makeShared<AirBlock>("AirBody"));
	air_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_a, c_f, mu_a);
	air_block.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexRelation water_air_complex(water_block, {&air_block});
	ContactRelation water_wall_contact(water_block, {&wall_boundary});
	ComplexRelation air_water_complex(air_block, {&water_block});
	ContactRelation air_wall_contact(air_block, {&wall_boundary});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	/** Initial condition with droplet impact velocity */
	// SimpleDynamics<DropletImpactSpeed> initial_condition(water_block);
	/** Initialize particle acceleration. */
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
	SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block, gravity_ptr);
	SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block, gravity_ptr);
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex>
		update_water_density_by_summation(water_wall_contact, water_air_complex);
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex>
		update_air_density_by_summation(air_wall_contact, air_water_complex);
	/** transport formulation for regularizing particle distribution. */
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>>
		air_transport_correction(air_wall_contact, air_water_complex, 2.0e-3);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>>
		water_transport_correction(water_air_complex, 2.0e-3);
	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step_size(water_block, U_max);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_air_advection_time_step_size(air_block, U_max);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_time_step_size(water_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_air_time_step_size(air_block);
	/** Time step size with considering sound speed. */
	ReduceDynamics<fluid_dynamics::SoundSpeedTimeStepSize> get_water_sound_time_step_size(water_block);
	ReduceDynamics<fluid_dynamics::SoundSpeedTimeStepSize> get_air_sound_time_step_size(air_block);
	/** Time step size with considering individual particle acceleration */
	ReduceDynamics<fluid_dynamics::ParticleAccelerationTimeStepSize> get_water_acc_time_step_size(water_block);
	ReduceDynamics<fluid_dynamics::ParticleAccelerationTimeStepSize> get_air_acc_time_step_size(air_block);
	/** Pressure relaxation for water by using position verlet time stepping. */
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemannWithWall>
        water_pressure_relaxation(water_wall_contact, water_air_complex);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
        water_density_relaxation(water_wall_contact, water_air_complex);
	/** Extend Pressure relaxation is used for air. */
	Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall>
		air_pressure_relaxation(air_wall_contact, air_water_complex, 2.0);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall>
		air_density_relaxation(air_wall_contact, air_water_complex);
	/** Viscous acceleration. */
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> air_viscous_acceleration(air_wall_contact, air_water_complex);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> water_viscous_acceleration(water_wall_contact, water_air_complex);
	/** Surface tension and wetting effects. */
	// InteractionWithUpdate<fluid_dynamics::FreeSurfaceIndicationComplex> surface_detection(water_air_complex.getInnerRelation(), water_wall_contact);
	// InteractionDynamics<fluid_dynamics::ColorFunctionGradientComplex> color_gradient(water_air_complex.getInnerRelation(), water_wall_contact);
	// InteractionDynamics<fluid_dynamics::ColorFunctionGradientInterpolationInner> color_gradient_interpolation(water_air_complex.getInnerRelation());
	// InteractionDynamics<fluid_dynamics::SurfaceTensionAccelerationInner> surface_tension_acceleration(water_air_complex.getInnerRelation(), tension_force);
	InteractionWithUpdate<fluid_dynamics::FreeSurfaceIndicationInner> surface_detection(water_air_complex.getInnerRelation());
	InteractionDynamics<fluid_dynamics::ColorFunctionGradientInner> color_gradient(water_air_complex.getInnerRelation());
	InteractionDynamics<fluid_dynamics::ColorFunctionGradientInterpolationInner> color_gradient_interpolation(water_air_complex.getInnerRelation());
	InteractionDynamics<fluid_dynamics::SurfaceTensionAccelerationInner> surface_tension_acceleration(water_air_complex.getInnerRelation(), tension_force);
	/** Wetting effects. */
	InteractionDynamics<fluid_dynamics::SurfaceNormWithWall> wetting_norm(water_wall_contact, contact_angle);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_boundary_normal_direction.exec();
	// initial_condition.exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 10;
	Real end_time = 0.05;				  /**< End time. */
	Real output_interval = end_time / 100; /**< Time stamps for output of body states. */
	Real delta_t1 = 0.0; //time step given by surface tension criterion
	Real delta_t2 = 0.0; //time step given by speed of sound criterion
	Real delta_t3 = 0.0; //time step given by particle acceleration criterion
	Real delta_t4 = 0.0; //time step given by viscous diffusion criterion
	Real dt = 0.0;						  /**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	TickCount t1 = TickCount::now();
	TickCount::interval_t interval;
	TickCount::interval_t interval_computing_time_step;
	TickCount::interval_t interval_computing_pressure_relaxation;
	TickCount::interval_t interval_updating_configuration;
	TickCount time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile();
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = TickCount::now();
			initialize_a_water_step.exec();
			initialize_a_air_step.exec();

			Real Dt_f = get_water_advection_time_step_size.exec();
			Real Dt_a = get_air_advection_time_step_size.exec();
			// Real Dt = SMIN(Dt_f, Dt_a);

			Real Dt_acc_f = get_water_acc_time_step_size.exec();
			Real Dt_acc_a = get_air_acc_time_step_size.exec();
			delta_t1 = SMIN(Dt_acc_f, Dt_acc_a);
			Real Dt_sound_f = get_water_sound_time_step_size.exec();
			Real Dt_sound_a = get_air_sound_time_step_size.exec();
			delta_t2 = SMIN(Dt_sound_f, Dt_sound_a);
			delta_t3 = 0.25*pow(1.3*particle_spacing_ref/gravity_g, 0.5);
			delta_t4 = 0.25*pow(0.1*rho0_f*pow(1.3*particle_spacing_ref,3)/(2*M_PI*tension_force),0.5);
			Real Dt = SMIN(Dt_f,Dt_a,delta_t1,delta_t2,delta_t3,delta_t4);

			update_water_density_by_summation.exec();
			update_air_density_by_summation.exec();
			air_transport_correction.exec();
			water_transport_correction.exec();

			air_viscous_acceleration.exec();
			water_viscous_acceleration.exec();

			surface_detection.exec();
			color_gradient.exec();
			color_gradient_interpolation.exec();
			wetting_norm.exec();
			surface_tension_acceleration.exec();

			interval_computing_time_step += TickCount::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = TickCount::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = get_water_time_step_size.exec();
				Real dt_a = get_air_time_step_size.exec();
				dt = SMIN(SMIN(dt_f, dt_a), Dt);

				water_pressure_relaxation.exec(dt);
				air_pressure_relaxation.exec(dt);

				water_density_relaxation.exec(dt);
				air_density_relaxation.exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_pressure_relaxation += TickCount::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = TickCount::now();

			water_block.updateCellLinkedListWithParticleSort(100);
			water_air_complex.updateConfiguration();
			water_wall_contact.updateConfiguration();

			air_block.updateCellLinkedListWithParticleSort(100);
			air_water_complex.updateConfiguration();
			air_wall_contact.updateConfiguration();

			interval_updating_configuration += TickCount::now() - time_instance;
		}

		TickCount t2 = TickCount::now();
		body_states_recording.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}

	TickCount t4 = TickCount::now();

	TickCount::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
			  << " seconds." << std::endl;
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
			  << interval_computing_time_step.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
			  << interval_computing_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
			  << interval_updating_configuration.seconds() << "\n";

	return 0;
}
