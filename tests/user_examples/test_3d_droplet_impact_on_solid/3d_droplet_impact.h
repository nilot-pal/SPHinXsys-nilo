/**
 * @file 	wetting.h
 * @brief 	Numerical parameters and body definition for 3D two-phase droplet impact.
 * @author 	Nilotpal Chakraborty
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real insert_circle_radius = 1.36e-3;			  /**< Radius of the droplet. */
// Real DL = 10.0 * insert_circle_radius;						   /**< Tank length. */
// Real DH = 10.0 * insert_circle_radius;							   /**< Tank height. */
// Real DW = 10.0 * insert_circle_radius;               // tank width
Real DL = 4.7 * insert_circle_radius;						   /**< Tank length. */
Real DH = 4.7 * insert_circle_radius;							   /**< Tank height. */
Real DW = 4.7 * insert_circle_radius;                   // tank width
Real particle_spacing_ref = 5e-5; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;	   /**< Extending width for BCs. */
Vecd circle_center(0.5 * DL, insert_circle_radius, 0.5 * DW);	  /**< Location of the cylinder center. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real gravity_g = 9.8;
Real rho0_f = 1000.0;								  /**< Reference density of water. */
Real rho0_a = 1.0;							  /**< Reference density of air. */
Real U_max = 4.0;								  /**< Characteristic velocity. */
// Real c_f = 10.0 * U_max;						  /**< Reference sound speed. */
Real c_f = 50.0;
Real mu_f = 5.0e-2;								  /**< Water viscosity. */
Real mu_a = 5.0e-5;								  /**< Air viscosity. */
Real contact_angle = (90.0 / 180.0) * 3.1415926; /**< Contact angle with Wall. */
Real tension_force = 0.073;
//----------------------------------------------------------------------
//	Geometric elements used in the shape modeling.
//----------------------------------------------------------------------
Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
Vecd halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
// std::vector<Vecd> createOuterWallShape()
// {
// 	std::vector<Vecd> outer_wall_shape;
// 	outer_wall_shape.push_back(Vecd(-BW, -BW, -BW));
// 	outer_wall_shape.push_back(Vecd(-BW, DH + BW, DW));
// 	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
// 	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
// 	outer_wall_shape.push_back(Vecd(-BW, -BW));

// 	return outer_wall_shape;
// }

// std::vector<Vecd> createInnerWallShape()
// {
// 	std::vector<Vecd> inner_wall_shape;
// 	inner_wall_shape.push_back(Vecd(0.0, 0.0));
// 	inner_wall_shape.push_back(Vecd(0.0, DH));
// 	inner_wall_shape.push_back(Vecd(DL, DH));
// 	inner_wall_shape.push_back(Vecd(DL, 0.0));
// 	inner_wall_shape.push_back(Vecd(0.0, 0.0));

// 	return inner_wall_shape;
// }
//----------------------------------------------------------------------
//	case-dependent geometric shape of water block.
//----------------------------------------------------------------------
// class WaterBlock : public ComplexShape
// {
// public:
// 	explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
// 	{
// 		Transform ball_centre_1(circle_center);
// 		add<TransformShape<GeometricShapeBall>>(Transform(ball_centre_1), circle_center, insert_circle_radius);
// 	}
// };
//----------------------------------------------------------------------
//	case-dependent geometric shape of air block.
//----------------------------------------------------------------------
class AirBlock : public ComplexShape
{
public:
	explicit AirBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Transform translation_wall_1(halfsize_inner);
		Transform ball_centre_2(circle_center);
		add<TransformShape<GeometricShapeBox>>(Transform(translation_wall_1), halfsize_inner);
		subtract<TransformShape<GeometricShapeBall>>(Transform(ball_centre_2), circle_center, insert_circle_radius);
	}
};
//----------------------------------------------------------------------
//	Wall boundary shape definition.
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
    		Transform translation_wall_2(halfsize_inner); // change of origin from Vecd(0,0,0) to halfsize_inner
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall_2), halfsize_outer);
        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_wall_2), halfsize_inner);
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
		vel_[index_i][1] = -2.64;
		vel_[index_i][2] = 0;
	}
};