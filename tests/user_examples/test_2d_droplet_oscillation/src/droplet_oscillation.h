/**
 * @file 	wetting.h
 * @brief 	Numerical parameters and body definition for 2D two-phase wetting flow.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real insert_circle_radius = 0.2;			  /**< Radius of the droplet. */
Real DL = 0.5;						   /**< Domain half length. */
Real DH = 0.5;							   /**< Domain half height. */
Real particle_spacing_ref = insert_circle_radius / 20.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;	   /**< Extending width for BCs. */
// Real particle_spacing_ref = 2 / 40.0; /**< Initial reference particle spacing. */
// Real BW = particle_spacing_ref * 2;	   /**< Extending width for BCs. */
Vec2d insert_circle_center(0.0, 0.0);	  /**< Location of the droplet center. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL-BW, -DH-BW), Vec2d(DL+BW, DH+BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;								  /**< Reference density of water. */
Real rho0_a = 1.0e-3;							  /**< Reference density of air. */
Real gravity_g = 9.8;                             /**< Gravity. */
Real U_max = 0.3645;								  /**< Characteristic velocity (MATLAB calculation) */	
// Real U_max = 1.0;
Real c_f = 34.3;						  /**< Reference sound speed, (Mach no. = 0.1) */
// Real c_f = 10.0;
Real mu_f = 5.0e-2;								  /**< Water viscosity. */
Real mu_a = 5.0e-4;								  /**< Air viscosity. */
Real tension_force = 1;
//----------------------------------------------------------------------
//	Geometric elements used in the shape modeling.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
	// geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-0.2, -0.2));
	water_block_shape.push_back(Vecd(-0.2, 0.2));
	water_block_shape.push_back(Vecd(0.2, 0.2));
	water_block_shape.push_back(Vecd(0.2, -0.2));
	water_block_shape.push_back(Vecd(-0.2, -0.2));
	return water_block_shape;
}

std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-BW-DL, -BW-DH));
	outer_wall_shape.push_back(Vecd(-BW-DL, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW-DH));
	outer_wall_shape.push_back(Vecd(-BW-DL, -BW-DH));

	return outer_wall_shape;
}

std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL, -DH));
	inner_wall_shape.push_back(Vecd(-DL, DH));
	inner_wall_shape.push_back(Vecd(DL, DH));
	inner_wall_shape.push_back(Vecd(DL, -DH));
	inner_wall_shape.push_back(Vecd(-DL, -DH));

	return inner_wall_shape;
}
//----------------------------------------------------------------------
//	case-dependent geometric shape of water block.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
	}
};

// class WaterBlock : public MultiPolygonShape
// {
// public:
// 	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
// 	{
// 		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
// 	}
// };

//----------------------------------------------------------------------
//	case-dependent geometric shape of air block.
//----------------------------------------------------------------------
class AirBlock : public MultiPolygonShape
{
public:
	explicit AirBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);
		multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
	}
};

// class AirBlock : public MultiPolygonShape
// {
// public:
// 	explicit AirBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
// 	{
// 		multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);
// 		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::sub);
// 	}
// };
//----------------------------------------------------------------------
//	Wall boundary shape definition.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		std::vector<Vecd> outer_shape = createOuterWallShape();
		std::vector<Vecd> inner_shape = createInnerWallShape();
		multi_polygon_.addAPolygon(outer_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_shape, ShapeBooleanOps::sub);
	}
};
// observer location at centre of mass of upper 1/4th part of droplet
// StdVec<Vecd> observation_location = {Vecd(-4*insert_circle_radius/(3*M_PI), 4*insert_circle_radius/(3*M_PI))};
//	No-slip boundary condition
//----------------------------------------------------------------------
// struct NoSlipVelocity
// {
// 	Real u_ref_;

// 	template <class BoundaryConditionType>
// 	NoSlipVelocity(BoundaryConditionType &boundary_condition)
// 		: u_ref_(0.0) {}

// 	Vecd operator()(Vecd &position, Vecd &velocity)
// 	{
// 		Vecd target_velocity = Vecd::Zero();
// 		target_velocity[0] = u_ref_;
// 		target_velocity[0] = u_ref_;
// 		return target_velocity;
// 	}
// };
//	Application dependent initial condition
//----------------------------------------------------------------------
class DropletInitialSpeed
	: public fluid_dynamics::FluidInitialCondition
{
public:
	explicit DropletInitialSpeed(SPHBody &sph_body)
		: fluid_dynamics::FluidInitialCondition(sph_body){};

	void update(size_t index_i, Real dt)
	{
		/** initialize particle velocity */
		Real V_0 = 1.0;
		Real r_0 = 0.05;
		Real r = sqrt(pow(pos_[index_i][0],2) + pow(pos_[index_i][1],2));
		vel_[index_i][0] = V_0*(pos_[index_i][0]/r_0)*(1-pow(pos_[index_i][1],2)/(r_0*r))*exp(-r/r_0);
		vel_[index_i][1] = -V_0*(pos_[index_i][1]/r_0)*(1-pow(pos_[index_i][0],2)/(r_0*r))*exp(-r/r_0);
	}
};