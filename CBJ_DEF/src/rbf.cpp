#include "../include/rbf.h"
#include <random>


void RadialBasisFunction::set_solve_matrix(Eigen::MatrixXd& M, std::vector<Point>& xc, bool add_poly) const
{
	unsigned int nc = xc.size();
	unsigned int matrix_dim = nc;
	if (add_poly) matrix_dim += nDim + 1;
	M = Eigen::MatrixXd::Zero(matrix_dim, matrix_dim);
	fill_matrix_rbf_part(M, xc, xc);
	if (add_poly)
	{
		fill_matrix_linear_polynomial_part(M, xc, nc);
		fill_matrix_constrain_part(M, nc, xc);
	}
}

void RadialBasisFunction::fill_matrix_rbf_part(Eigen::MatrixXd& A,
											   std::vector<Point>& pos_i, std::vector<Point>& pos_j) const
{
	Point x_i;
	Point x_j;

	unsigned int imax = pos_i.size();
	unsigned int jmax = pos_j.size();

	for (unsigned j = 0; j < jmax; ++j)
	{
		x_j = pos_j[j];
		for (unsigned i = 0; i < imax; ++i)
		{
			x_i = pos_i[i];
			A(i, j) = (*this)(x_i, x_j);
		}
	}
}

void RadialBasisFunction::fill_matrix_linear_polynomial_part(Eigen::MatrixXd& A,
															 std::vector<Point>& pos_i, const unsigned nc) const
{
	unsigned int imax = pos_i.size();
	for (unsigned i = 0; i < imax; ++i)
		A(i, nc) = 1.0;
	for (int iDim = 0; iDim < nDim; ++iDim)
		for (unsigned i = 0; i < imax; ++i)
			A(i, nc + iDim + 1) = pos_i[iDim][i];
}


void RadialBasisFunction::fill_matrix_constrain_part(Eigen::MatrixXd& A,
													 const unsigned ni, std::vector<Point>& pos_j) const
{
	unsigned int jmax = pos_j.size();
	for (unsigned j = 0; j < jmax; ++j)
	{
		A(ni, j) = 1.0;
		for (int iDim = 0; iDim < nDim; ++iDim)
			A(ni + iDim + 1, j) = pos_j[iDim][j];
	}
}
