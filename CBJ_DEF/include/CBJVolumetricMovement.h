#pragma once
#include "../include/SU2_DEF.hpp"
#include "../include/rbf.h"
#include <random>
#include <set>

class CbjVolumetricMovement : public CGridMovement
{
public:
	unsigned short nDim;		/*!< \brief Number of dimensions. */
	unsigned short nVar;		/*!< \brief Number of variables. */

	unsigned long nPoint;		/*!< \brief Number of points. */
	unsigned long nPointDomain;		/*!< \brief Number of points in the domain. */
	bool rbfControlPointSet;
	double allowErr;
	int maxPointNum;
	
	Wendland_C2 *rbf;

	CbjVolumetricMovement(CGeometry *geometry);
	void RbfInitialize(const double & supportRadius, const bool & rbfControlPointSet, const double & allowErr, const int & maxPointNum);
	~CbjVolumetricMovement();
	double Check_Grid(CGeometry* geometry);
	void SetVolume_Deformation(CGeometry *geometry, CConfig *config, bool UpdateGeo);
	/*!
	* \brief Compute the shape functions for hexahedron
	* \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	*/
	double GetHexa_Volume(double CoordCorners[8][3]);

	/*!
	* \brief Compute the shape functions for hexahedron
	* \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	*/
	double GetTetra_Volume(double CoordCorners[8][3]);

	/*!
	* \brief Compute the shape functions for hexahedron
	* \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	*/
	double GetWedge_Volume(double CoordCorners[8][3]);

	/*!
	* \brief Compute the shape functions for hexahedron
	* \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	*/
	double GetPyram_Volume(double CoordCorners[8][3]);

	/*!
	* \brief Compute the shape functions for hexahedron
	* \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	*/
	double GetTriangle_Area(double CoordCorners[8][3]);

	/*!
	* \brief Compute the shape functions for hexahedron
	* \param[in] CoordCorners - coordinates of the cornes of the hexahedron.
	*/
	double GetRectangle_Area(double CoordCorners[8][3]);

	static long random_select(long min, long max);

	double RbfSelectPoint();
	void RbfCalculateAlpha(vector<Point>& alpha, vector<Point>& control_point, vector<Point>& control_point_var);
	void RbfPointSeletion(vector<Point>& control_point, vector<Point>& control_point_var, vector<Point>& surface_point, vector<Point>& surface_point_var, vector<unsigned long>& globalIndex);
	void RbfApplyDeform(CGeometry *geometry, vector<Point>& control_point, vector<Point>& control_point_var);
};


void OutTime(double startTime);
