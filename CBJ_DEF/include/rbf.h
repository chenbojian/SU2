#ifndef RBF_TOOL_H
#define RBF_TOOL_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>
using namespace std;
struct Point
{
	double p[3];
	Point()
	{
		p[0] = 0.0;
		p[1] = 0.0;
		p[2] = 0.0;
	}
	Point(double a, double b) 
	{		
		p[0] = a;
		p[1] = b;
		p[2] = 0.0;

	}
	Point(double a, double b, double c) 
	{
		p[0] = a;
		p[1] = b;
		p[2] = c;
	}
	Point(double *p,int n)
	{
		if (n == 3)
		{
			copy_n(p, 3, this->p);
			
		}
		else if (n == 2)
		{
			copy_n(p, 2, this->p);
			this->p[2] = 0.0;
		}
		else
		{
			Point();
		}
	}
	double & operator[](int i)
	{
		return p[i];
	}
	double * begin()
	{
		return p;
	}
	double * end()
	{
		return p + 3;
	}
};
class RadialBasisFunction
{
private:
	virtual double basic_fun(double x) const = 0;

	void fill_matrix_rbf_part(Eigen::MatrixXd& A, std::vector<Point>& pos_i, std::vector<Point>& pos_j) const;
	void fill_matrix_linear_polynomial_part(Eigen::MatrixXd& A,  std::vector<Point>& pos_i, const unsigned nc) const;
	void fill_matrix_constrain_part(Eigen::MatrixXd& A, const unsigned ni, std::vector<Point>& pos_j) const;

protected:
	int nDim;
	double support_radius;

public:
	RadialBasisFunction(int n) : nDim(n), support_radius(1.0) { }
	RadialBasisFunction(int n, double r) : nDim(n), support_radius(r)
	{
//		const double eps = 1.e-6;
//		assertion(support_radius > eps, "The Support Radius Must be Positive");
	}

	void set_suprad(double r)
	{
//		const double eps = 1.e-6;
//		assertion(r > eps, "The Support Radius Must be Positive");
		support_radius = r;
	}

	double operator( ) (double x) const
	{
		x = fabs(x);
		return basic_fun(x);
	}
	double operator( ) ( Point & v) const
	{
		double x = 0.0;
		for (int iDim = 0; iDim<nDim; ++iDim) x += v[iDim] * v[iDim];
		x = sqrt(x);
		return basic_fun(x);
	}
	double operator( ) ( Point & xp, Point & xp0) const
	{
		Point v;
		for (int iDim = 0; iDim<nDim; ++iDim) v[iDim] = xp[iDim] - xp0[iDim];
		return this->operator()(v);
	}

	void set_solve_matrix(Eigen::MatrixXd& M,  std::vector<Point>& xc, bool add_poly) const;

};

class Wendland_C2 : public RadialBasisFunction
{
private:
	double basic_fun(double x) const
	{
		double r = x / support_radius;
		double phi = 0.0;
		const double eps = 1.e-10;
		if (r < (1.0 - eps)) phi = pow((1.0 - r), 4) * (4.0 * r + 1.0);
		return phi;
	}
public:
	Wendland_C2(int n, double r) : RadialBasisFunction(n, r) { }
};

#endif // RBF_TOOL_H

