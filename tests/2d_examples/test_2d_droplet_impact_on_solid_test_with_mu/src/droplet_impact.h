/**
 * @file 	wetting.h
 * @brief 	Numerical parameters and body definition for 2D two-phase wetting flow.
 * @author 	Nilotpal Chakraborty
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real insert_circle_radius = 1e-3;			  /**< Radius of the droplet. */
Real DL = 50 * insert_circle_radius;						   /**< Tank length. */
Real DH = 50 * insert_circle_radius;							   /**< Tank height. */
Real particle_spacing_ref = insert_circle_radius / 20.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;	   /**< Extending width for BCs. */
Vec2d insert_circle_center(0.5 * DL, insert_circle_radius);	  /**< Location of the cylinder center. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real gravity_g = 9.8;
Real rho0_f = 1000.0;								  /**< Reference density of water. */
Real rho0_a = 1.0;							  /**< Reference density of air. */
Real U_max = 2.65;								  /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;						  /**< Reference sound speed. */
Real mu_f = 0.2;								  /**< Water viscosity. */
Real mu_a = 2.0e-3;								  /**< Air viscosity. */
Real contact_angle = (90.0 / 180.0) * 3.1415926; /**< Contact angle with Wall. */
Real tension_force = 0.073;
//----------------------------------------------------------------------
//	Geometric elements used in the shape modeling.
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, -BW));

	return outer_wall_shape;
}

std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(0.0, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, DH));
	inner_wall_shape.push_back(Vecd(DL, DH));
	inner_wall_shape.push_back(Vecd(DL, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, 0.0));

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
//	Application dependent initial condition
//----------------------------------------------------------------------
class DropletImpactSpeed
	: public fluid_dynamics::FluidInitialCondition
{
public:
	explicit DropletImpactSpeed(SPHBody &sph_body)
		: fluid_dynamics::FluidInitialCondition(sph_body){};

	void update(size_t index_i, Real dt)
	{
		/** initialize particle velocity */
		vel_[index_i][0] = 0;
		vel_[index_i][1] = -1.0;
	}
};
