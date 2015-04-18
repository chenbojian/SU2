#include "../include/CBJVolumetricMovement.h"
using namespace std;


CbjVolumetricMovement::CbjVolumetricMovement(CGeometry* geometry)
{
	nDim = geometry->GetnDim();
}

void CbjVolumetricMovement::RbfInitialize(const double& supportRadius, const bool& rbfControlPointSet, const double& allowErr, const int& maxPointNum)
{
	rbf = new Wendland_C2(nDim, supportRadius);
	this->rbfControlPointSet = rbfControlPointSet;
	this->allowErr = allowErr;
	this->maxPointNum = maxPointNum;
}

CbjVolumetricMovement::~CbjVolumetricMovement()
{
	delete rbf;
}

double CbjVolumetricMovement::Check_Grid(CGeometry *geometry) {

	unsigned long iElem, ElemCounter = 0, PointCorners[8];
	double Area = 0.0, Volume = 0.0, MaxArea = -1E22, MaxVolume = -1E22, MinArea = 1E22, MinVolume = 1E22, CoordCorners[8][3];
	unsigned short nNodes = 0, iNodes, iDim;
	bool RightVol = true;

	int rank = MASTER_NODE;

#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	/*--- Load up each triangle and tetrahedron to check for negative volumes. ---*/

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     nNodes = 3;
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    nNodes = 4;
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  nNodes = 4;
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      nNodes = 5;
		if (geometry->elem[iElem]->GetVTK_Type() == WEDGE)        nNodes = 6;
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   nNodes = 8;

		for (iNodes = 0; iNodes < nNodes; iNodes++) {
			PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
			for (iDim = 0; iDim < nDim; iDim++) {
				CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
			}
		}

		/*--- Triangles ---*/

		if (nDim == 2) {

			if (nNodes == 3) Area = GetTriangle_Area(CoordCorners);
			if (nNodes == 4) Area = GetRectangle_Area(CoordCorners);

			if (Area >= -EPS) RightVol = true;
			else RightVol = false;;

			MaxArea = max(MaxArea, Area);
			MinArea = min(MinArea, Area);

		}

		/*--- Tetrahedra ---*/
		if (nDim == 3) {

			if (nNodes == 4) Volume = GetTetra_Volume(CoordCorners);
			if (nNodes == 5) Volume = GetPyram_Volume(CoordCorners);
			if (nNodes == 6) Volume = GetWedge_Volume(CoordCorners);
			if (nNodes == 8) Volume = GetHexa_Volume(CoordCorners);

			if (Volume >= -EPS) RightVol = true;
			else RightVol = false;;

			MaxVolume = max(MaxVolume, Volume);
			MinVolume = min(MinVolume, Volume);

		}

		if (!RightVol) ElemCounter++;

	}

#ifdef HAVE_MPI
	unsigned long ElemCounter_Local = ElemCounter; ElemCounter = 0;
	double MaxVolume_Local = MaxVolume; MaxVolume = 0.0;
	double MinVolume_Local = MinVolume; MinVolume = 0.0;
	MPI_Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&MaxVolume_Local, &MaxVolume, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&MinVolume_Local, &MinVolume, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

	if ((ElemCounter != 0) && (rank == MASTER_NODE))
		cout << "There are " << ElemCounter << " elements with negative volume.\n" << endl;

	if (nDim == 2) return MinArea;
	else return MinVolume;

}

void CbjVolumetricMovement::SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo)
{
	unsigned long IterLinSol = 0, Smoothing_Iter, iNonlinear_Iter;
	unsigned long iMarker, iVertex, iPoint, iDim, total_index;
	double MinVolume, NumError, Tol_Factor;
	bool Screen_Output;
	int rank = MASTER_NODE;
	int size = SINGLE_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
	/*--- Initialize the number of spatial dimensions, length of the state
	vector (same as spatial dimensions for grid deformation), and grid nodes. ---*/

	nDim = geometry->GetnDim();
	nVar = geometry->GetnDim();
	nPoint = geometry->GetnPoint();
	nPointDomain = geometry->GetnPointDomain();
	vector<double *> surface_point, surface_point_var;
	vector<Point> control_point, control_point_var;
	vector<unsigned long> globalIndex;
	vector<Point> all_surface_point, all_surface_point_var;
	vector<unsigned long> allGlobalIndex;
	double *coordVar;
	double *coord;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
	{
		if (config->GetMarker_All_KindBC(iMarker) != SYMMETRY_PLANE
			&& config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE
			&& config->GetMarker_All_KindBC(iMarker) != FAR_FIELD)
		{
			cout << geometry->GetnVertex(iMarker) << ':' << iMarker << endl;
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++)
			{
				coordVar = geometry->vertex[iMarker][iVertex]->GetVarCoord();
				coord = geometry->vertex[iMarker][iVertex]->GetCoord();
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				globalIndex.push_back(geometry->node[iPoint]->GetGlobalIndex());
				surface_point.push_back(coord);
				surface_point_var.push_back(coordVar);
			}
		}
	}
#ifdef HAVE_MPI
	//all send to master
	MPI_Barrier(MPI_COMM_WORLD);
	if (size > 1)
	{
		if (rank == MASTER_NODE)
		{
			Point sp, sp_var;
			for (iPoint = 0; iPoint < surface_point.size(); iPoint++)
			{
				for (iDim = 0; iDim < nDim; iDim++)
				{
					sp[iDim] = surface_point[iPoint][iDim];
					sp_var[iDim] = surface_point_var[iPoint][iDim];
				}
				allGlobalIndex.push_back(globalIndex[iPoint]);
				all_surface_point.push_back(sp);
				all_surface_point_var.push_back(sp_var);
			}
			for (int iRank = 1; iRank < size; iRank++)
			{
				unsigned long nLen;
				int source = iRank, tag = iRank;
				unsigned long * bufGlobalIndex;
				double * bufSurfacePoint, *bufSurfacePointVar;
				MPI_Status status;
				MPI_Recv(&nLen, 1, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, &status);
				bufGlobalIndex = new unsigned long[nLen];
				bufSurfacePoint = new double[nLen*nDim];
				bufSurfacePointVar = new double[nLen*nDim];

				MPI_Recv(bufGlobalIndex, nLen, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(bufSurfacePoint, nLen*nDim, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(bufSurfacePointVar, nLen*nDim, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
				for (unsigned long i = 0; i < nLen; i++)
				{
					for (iDim = 0; iDim < nDim; iDim++)
					{
						sp[iDim] = bufSurfacePoint[i*nDim + iDim];
						sp_var[iDim] = bufSurfacePointVar[i*nDim + iDim];
					}
					allGlobalIndex.push_back(bufGlobalIndex[i]);
					all_surface_point.push_back(sp);
					all_surface_point_var.push_back(sp_var);
				}
				delete bufSurfacePoint;
				delete bufGlobalIndex;
				delete bufSurfacePointVar;
			}
		}
		else
		{
			unsigned long nLen = globalIndex.size();
			int dest = MASTER_NODE, tag = rank;
			unsigned long * bufGlobalIndex;
			double * bufSurfacePoint, *bufSurfacePointVar;
			MPI_Ssend(&nLen, 1, MPI_UNSIGNED_LONG, dest, tag, MPI_COMM_WORLD);
			//send global index
			bufGlobalIndex = new unsigned long[nLen];
			bufSurfacePoint = new double[nLen*nDim];
			bufSurfacePointVar = new double[nLen*nDim];
			for (unsigned long i = 0; i < nLen; i++)
			{
				bufGlobalIndex[i] = globalIndex[i];
				for (iDim = 0; iDim < nDim; iDim++)
				{
					bufSurfacePoint[i*nDim + iDim] = surface_point[i][iDim];
					bufSurfacePointVar[i*nDim + iDim] = surface_point_var[i][iDim];
				}
			}
			MPI_Ssend(bufGlobalIndex, nLen, MPI_UNSIGNED_LONG, dest, tag, MPI_COMM_WORLD);
			MPI_Ssend(bufSurfacePoint, nLen*nDim, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			MPI_Ssend(bufSurfacePointVar, nLen*nDim, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);

			delete bufGlobalIndex;
			delete bufSurfacePoint;
			delete bufSurfacePointVar;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		//master to all
		if (rank == MASTER_NODE)
		{
			unsigned long nLen = allGlobalIndex.size();
			unsigned long * bufGlobalIndex;
			double * bufSurfacePoint, *bufSurfacePointVar;
			bufGlobalIndex = new unsigned long[nLen];
			bufSurfacePoint = new double[nLen*nDim];
			bufSurfacePointVar = new double[nLen*nDim];
			for (unsigned long i = 0; i < nLen; i++)
			{
				bufGlobalIndex[i] = allGlobalIndex[i];
				for (iDim = 0; iDim < nDim; iDim++)
				{
					bufSurfacePoint[i*nDim + iDim] = all_surface_point[i][iDim];
					bufSurfacePointVar[i*nDim + iDim] = all_surface_point_var[i][iDim];
				}
			}
			for (int iRank = 1; iRank < size; iRank++)
			{
				int dest = iRank, tag = iRank;
				MPI_Ssend(&nLen, 1, MPI_UNSIGNED_LONG, dest, tag, MPI_COMM_WORLD);
				MPI_Ssend(bufGlobalIndex, nLen, MPI_UNSIGNED_LONG, dest, tag, MPI_COMM_WORLD);
				MPI_Ssend(bufSurfacePoint, nLen*nDim, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
				MPI_Ssend(bufSurfacePointVar, nLen*nDim, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			}
			delete bufGlobalIndex;
			delete bufSurfacePoint;
			delete bufSurfacePointVar;
		}
		else
		{
			unsigned long nLen;
			int source = MASTER_NODE, tag = rank;
			unsigned long * bufGlobalIndex;
			double * bufSurfacePoint, *bufSurfacePointVar;
			MPI_Status status;
			MPI_Recv(&nLen, 1, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, &status);
			bufGlobalIndex = new unsigned long[nLen];
			bufSurfacePoint = new double[nLen*nDim];
			bufSurfacePointVar = new double[nLen*nDim];

			MPI_Recv(bufGlobalIndex, nLen, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, &status);
			MPI_Recv(bufSurfacePoint, nLen*nDim, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
			MPI_Recv(bufSurfacePointVar, nLen*nDim, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
			Point sp, sp_var;
			for (unsigned long i = 0; i < nLen; i++)
			{
				for (iDim = 0; iDim < nDim; iDim++)
				{
					sp[iDim] = bufSurfacePoint[i*nDim + iDim];
					sp_var[iDim] = bufSurfacePointVar[i*nDim + iDim];
				}
				allGlobalIndex.push_back(bufGlobalIndex[i]);
				all_surface_point.push_back(sp);
				all_surface_point_var.push_back(sp_var);
			}
			delete bufGlobalIndex;
			delete bufSurfacePoint;
			delete bufSurfacePointVar;
		}
	}
	else
	{
		Point sp, sp_var;
		for (iPoint = 0; iPoint < surface_point.size(); iPoint++)
		{
			for (iDim = 0; iDim < nDim; iDim++)
			{
				sp[iDim] = surface_point[iPoint][iDim];
				sp_var[iDim] = surface_point_var[iPoint][iDim];
			}
			allGlobalIndex.push_back(globalIndex[iPoint]);
			all_surface_point.push_back(sp);
			all_surface_point_var.push_back(sp_var);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

#endif
	if (rbfControlPointSet)
	{
		if (rank == MASTER_NODE)
		{
			cout << "begin Point selection" << endl;
			RbfPointSeletion(control_point, control_point_var, all_surface_point, all_surface_point_var, allGlobalIndex);

		}
	}
	else
	{
		if (rank == MASTER_NODE)
		{
			ifstream globalIndexFile;
			int idx;
			globalIndexFile.open("globalindex.dat", ios::in);
			if (globalIndexFile.fail())
			{
				exit(0);
			}
			while (globalIndexFile >> idx)
			{
				//can optimize by using unordered_map
				for (int i = 0; i < allGlobalIndex.size(); i++)
				{
					if (idx == allGlobalIndex[i])
					{
						Point &cp_var = all_surface_point_var[i];
						Point &cp = all_surface_point[i];
						control_point_var.push_back(cp_var);
						control_point.push_back(cp);
						break;
					}
				}
			}
			globalIndexFile.close();
			cout << "surface control point: " << control_point.size() << endl;
			unsigned long nLen;
			double* bufControlPoint;
			double* bufControlPointVar;
			nLen = control_point.size();
			bufControlPoint = new double[nLen*nDim];
			bufControlPointVar = new double[nLen*nDim];
			for (int i = 0; i < nLen; i++)
			{
				for (iDim = 0; iDim < nDim; iDim++)
				{
					bufControlPoint[i*nDim + iDim] = control_point[i][iDim];
					bufControlPointVar[i*nDim + iDim] = control_point_var[i][iDim];
				}
			}
			for (int iRank = 1; iRank < size; iRank++)
			{
				MPI_Ssend(&nLen, 1, MPI_UNSIGNED_LONG, iRank, iRank, MPI_COMM_WORLD);
				MPI_Ssend(bufControlPoint, nLen*nDim, MPI_DOUBLE, iRank, iRank, MPI_COMM_WORLD);
				MPI_Ssend(bufControlPointVar, nLen*nDim, MPI_DOUBLE, iRank, iRank, MPI_COMM_WORLD);
			}
			delete bufControlPoint;
			delete bufControlPointVar;
		}
		else
		{
			unsigned long nLen;
			double* bufControlPoint;
			double* bufControlPointVar;
			Point cp, cp_var;
			MPI_Status status;
			MPI_Recv(&nLen, 1, MPI_UNSIGNED_LONG, MASTER_NODE, rank, MPI_COMM_WORLD, &status);
			bufControlPoint = new double[nLen*nDim];
			bufControlPointVar = new double[nLen*nDim];
			MPI_Recv(bufControlPoint, nLen*nDim, MPI_DOUBLE, MASTER_NODE, rank, MPI_COMM_WORLD, &status);
			MPI_Recv(bufControlPointVar, nLen*nDim, MPI_DOUBLE, MASTER_NODE, rank, MPI_COMM_WORLD, &status);
			for (int i = 0; i < nLen; i++)
			{
				for (iDim = 0; iDim < nDim; iDim++)
				{
					cp[iDim] = bufControlPoint[i*nDim + iDim];
					cp_var[iDim] = bufControlPointVar[i*nDim + iDim];
				}
				control_point.push_back(cp);
				control_point_var.push_back(cp_var);
			}
			delete bufControlPoint;
			delete bufControlPointVar;
		}
		RbfApplyDeform(geometry, control_point, control_point_var);
	}
#ifdef HAVE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	MinVolume = Check_Grid(geometry);
	cout << "cbj--->MinVolume is :" << MinVolume << " on  " << rank << endl;

}

/*--- Deallocate vectors for the linear system. ---*/

double CbjVolumetricMovement::GetTriangle_Area(double CoordCorners[8][3]) {

	unsigned short iDim;
	double a[3] = { 0.0, 0.0, 0.0 }, b[3] = { 0.0, 0.0, 0.0 };
	double *Coord_0, *Coord_1, *Coord_2, Area;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[1];
	Coord_2 = CoordCorners[2];

	for (iDim = 0; iDim < nDim; iDim++) {
		a[iDim] = Coord_0[iDim] - Coord_2[iDim];
		b[iDim] = Coord_1[iDim] - Coord_2[iDim];
	}

	Area = 0.5*fabs(a[0] * b[1] - a[1] * b[0]);

	return Area;

}

double CbjVolumetricMovement::GetRectangle_Area(double CoordCorners[8][3]) {

	unsigned short iDim;
	double a[3] = { 0.0, 0.0, 0.0 }, b[3] = { 0.0, 0.0, 0.0 };
	double *Coord_0, *Coord_1, *Coord_2, Area;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[1];
	Coord_2 = CoordCorners[2];

	for (iDim = 0; iDim < nDim; iDim++) {
		a[iDim] = Coord_0[iDim] - Coord_2[iDim];
		b[iDim] = Coord_1[iDim] - Coord_2[iDim];
	}

	Area = 0.5*fabs(a[0] * b[1] - a[1] * b[0]);

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[2];
	Coord_2 = CoordCorners[3];

	for (iDim = 0; iDim < nDim; iDim++) {
		a[iDim] = Coord_0[iDim] - Coord_2[iDim];
		b[iDim] = Coord_1[iDim] - Coord_2[iDim];
	}

	Area += 0.5*fabs(a[0] * b[1] - a[1] * b[0]);

	return Area;

}

long CbjVolumetricMovement::random_select(long min, long max)
{
	default_random_engine dre(random_device{}());
	uniform_int_distribution<long> r(min, max);
	return r(dre);
}

void CbjVolumetricMovement::RbfCalculateAlpha(vector<Point>& alpha, vector<Point>& control_point, vector<Point>& control_point_var)
{
	Eigen::MatrixXd M, M_inv;
	rbf->set_solve_matrix(M, control_point, false);
	//	for (int i = 0; i < control_point.size(); i++)
	//	{
	//		for (int j = 0; j < control_point.size(); j++)
	//			cout << M(i, j) << ",";
	//		cout << endl;
	//	}
	//		
	M_inv = M.inverse();
	//	for (int i = 0; i < control_point.size(); i++)
	//	{
	//		for (int j = 0; j < control_point.size(); j++)
	//			cout << M_inv(i, j) << ",";
	//		cout << endl;
	//	}
	alpha.resize(control_point.size());
	for (int i = 0; i < control_point.size(); i++)
	{
		for (int j = 0; j < control_point.size(); j++)
		{
			for (int iDim = 0; iDim < nDim; iDim++)
			{
				alpha[i][iDim] += M_inv(i, j)*control_point_var[j][iDim];
			}
		}
	}

}

void CbjVolumetricMovement::RbfPointSeletion(vector<Point>& control_point, vector<Point>& control_point_var,
											 vector<Point>& surface_point, vector<Point>& surface_point_var,
											 vector<unsigned long>& globalIndex)
{
	int iDim;
	double allowErr = this->allowErr, maxSelectErr = 0.0;
	int maxSelectIdx = 0;
	const int initialSelectPoint = 4;
	set<unsigned long> selectedIndex;
	ofstream controlPointGlobalIndex;
	controlPointGlobalIndex.open("globalindex.dat", ios::out);
	for (int i = 0; i < initialSelectPoint; i++)
	{
		unsigned long idx = random_select(0, surface_point.size() - 1);
		if (selectedIndex.find(idx) != selectedIndex.end())
		{
			i--;
			continue;
		}
		control_point.push_back(surface_point[idx]);
		control_point_var.push_back(surface_point_var[idx]);
		controlPointGlobalIndex << globalIndex[idx] << endl;
		selectedIndex.insert(idx);
	}
	//start selection
	while (control_point.size() < this->maxPointNum)
	{
		maxSelectErr = 0.0;
		vector<Point> alpha;
		vector<Point> selection_point_var;
		//		cout << "1----" << endl;
		selection_point_var.resize(control_point.size(), Point(1.0, 1.0, 1.0));
		//		cout << "2----" << endl;
		RbfCalculateAlpha(alpha, control_point, selection_point_var);
		int imax = surface_point.size();
		int jmax = control_point.size();
		Point sp, cp;
		Point sp_var;
		for (unsigned long i = 0; i < imax; ++i)
		{
			for (iDim = 0; iDim < nDim; ++iDim)
				sp_var[iDim] = 0.0;
			for (iDim = 0; iDim < nDim; ++iDim)
				sp[iDim] = surface_point[i][iDim];
			for (int j = 0; j < jmax; ++j)
			{
				for (iDim = 0; iDim < nDim; ++iDim)
					cp[iDim] = control_point[j][iDim];
				for (iDim = 0; iDim < nDim; ++iDim)
				{
					sp_var[iDim] += (*rbf)(sp, cp)*alpha[j][iDim];
				}
			}
			//			cout << "3----"<<i << endl;

			double err = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				err += (sp_var[iDim] - 1.0)*(sp_var[iDim] - 1.0);
			err = sqrt(err);
			if (err > maxSelectErr)
			{
				if (selectedIndex.find(i) != selectedIndex.end())
					continue;
				maxSelectErr = err;
				maxSelectIdx = i;
			}
		}
		cout << "maxErr " << maxSelectErr << endl;
		cout << "maxindex " << maxSelectIdx << endl;
		if (maxSelectErr > allowErr)
		{
			control_point.push_back(surface_point[maxSelectIdx]);
			control_point_var.push_back(surface_point_var[maxSelectIdx]);
			controlPointGlobalIndex << globalIndex[maxSelectIdx] << endl;
			selectedIndex.insert(maxSelectIdx);
		}
		else
		{
			break;
		}
		//		cout << "Percent" << control_point.size() / 6 <<"% "<< '\r';
	}
	controlPointGlobalIndex.close();
	cout << "maxErr:" << maxSelectErr << " " << control_point.size() << endl;
	//真实变形误差
	maxSelectErr = 0.0;
	vector<Point> alpha;
	RbfCalculateAlpha(alpha, control_point, control_point_var);
	int imax = surface_point.size();
	int jmax = control_point.size();
	Point sp, cp;
	Point sp_var;
	for (int i = 0; i < imax; ++i)
	{
		fill_n(sp_var.begin(), nDim, 0.0);
		for (iDim = 0; iDim < nDim; ++iDim)
			sp[iDim] = surface_point[i][iDim];
		for (int j = 0; j < jmax; ++j)
		{
			for (iDim = 0; iDim < nDim; ++iDim)
				cp[iDim] = control_point[j][iDim];
			for (iDim = 0; iDim < nDim; ++iDim)
			{
				sp_var[iDim] += (*rbf)(sp, cp)*alpha[j][iDim];
			}
		}
		double err = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			err += (sp_var[iDim] - surface_point_var[i][iDim])*(sp_var[iDim] - surface_point_var[i][iDim]);
		err = sqrt(err);
		if (err > maxSelectErr)
		{
			maxSelectErr = err;
			maxSelectIdx = i;
		}

	}
	cout << "real deform max err is :" << maxSelectErr << endl
		<< "global index is :" << globalIndex[maxSelectIdx] << endl;
}

void CbjVolumetricMovement::RbfApplyDeform(CGeometry* geometry, vector<Point>& control_point, vector<Point>& control_point_var)
{
	int iDim;
	long nControlPoint = control_point.size();
	long nVolumePoint = geometry->GetnPoint();

	vector<Point> alpha;
	RbfCalculateAlpha(alpha, control_point, control_point_var);
	int imax = nVolumePoint;
	int jmax = nControlPoint;
	Point vp, cp;
	Point vp_var;
	for (int i = 0; i < imax; ++i)
	{
		for (iDim = 0; iDim < nDim; iDim++)
			vp_var[iDim] = 0.0;
		for (iDim = 0; iDim < nDim; ++iDim)
			vp[iDim] = geometry->node[i]->GetCoord(iDim);
		for (int j = 0; j < jmax; ++j)
		{
			for (iDim = 0; iDim < nDim; ++iDim)
				cp[iDim] = control_point[j][iDim];
			for (iDim = 0; iDim < nDim; ++iDim)
				vp_var[iDim] += (*rbf)(vp, cp)*alpha[j][iDim];
		}
		//deform
		for (iDim = 0; iDim < nDim; iDim++) {
			double new_coord = geometry->node[i]->GetCoord(iDim) + vp_var[iDim];
			if (fabs(new_coord) < EPS*EPS) new_coord = 0.0;
			geometry->node[i]->SetCoord(iDim, new_coord);
		}

	}
}

void OutTime(double startTime)
{
#ifdef HAVE_MPI
	double endTime = MPI_Wtime();
	cout << "inverse time is :" << endTime - startTime << endl;
#endif
}


double CbjVolumetricMovement::GetTetra_Volume(double CoordCorners[8][3]) {

	unsigned short iDim;
	double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
	double r1[3] = { 0.0, 0.0, 0.0 }, r2[3] = { 0.0, 0.0, 0.0 }, r3[3] = { 0.0, 0.0, 0.0 }, CrossProduct[3] = { 0.0, 0.0, 0.0 }, Volume;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[1];
	Coord_2 = CoordCorners[2];
	Coord_3 = CoordCorners[3];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	return Volume;

}

double CbjVolumetricMovement::GetPyram_Volume(double CoordCorners[8][3]) {

	unsigned short iDim;
	double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
	double r1[3] = { 0.0, 0.0, 0.0 }, r2[3] = { 0.0, 0.0, 0.0 }, r3[3] = { 0.0, 0.0, 0.0 }, CrossProduct[3] = { 0.0, 0.0, 0.0 }, Volume;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[1];
	Coord_2 = CoordCorners[2];
	Coord_3 = CoordCorners[4];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[2];
	Coord_2 = CoordCorners[3];
	Coord_3 = CoordCorners[4];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	return Volume;

}

double CbjVolumetricMovement::GetWedge_Volume(double CoordCorners[8][3]) {

	unsigned short iDim;
	double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
	double r1[3] = { 0.0, 0.0, 0.0 }, r2[3] = { 0.0, 0.0, 0.0 }, r3[3] = { 0.0, 0.0, 0.0 }, CrossProduct[3] = { 0.0, 0.0, 0.0 }, Volume;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[2];
	Coord_2 = CoordCorners[1];
	Coord_3 = CoordCorners[5];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[5];
	Coord_2 = CoordCorners[1];
	Coord_3 = CoordCorners[4];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[5];
	Coord_2 = CoordCorners[4];
	Coord_3 = CoordCorners[3];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	return Volume;

}

double CbjVolumetricMovement::GetHexa_Volume(double CoordCorners[8][3]) {

	unsigned short iDim;
	double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
	double r1[3] = { 0.0, 0.0, 0.0 }, r2[3] = { 0.0, 0.0, 0.0 }, r3[3] = { 0.0, 0.0, 0.0 }, CrossProduct[3] = { 0.0, 0.0, 0.0 }, Volume;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[1];
	Coord_2 = CoordCorners[2];
	Coord_3 = CoordCorners[5];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[2];
	Coord_2 = CoordCorners[7];
	Coord_3 = CoordCorners[5];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[2];
	Coord_2 = CoordCorners[3];
	Coord_3 = CoordCorners[7];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	Coord_0 = CoordCorners[0];
	Coord_1 = CoordCorners[5];
	Coord_2 = CoordCorners[7];
	Coord_3 = CoordCorners[4];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	Coord_0 = CoordCorners[2];
	Coord_1 = CoordCorners[7];
	Coord_2 = CoordCorners[5];
	Coord_3 = CoordCorners[6];

	for (iDim = 0; iDim < nDim; iDim++) {
		r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
		r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
		r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
	}

	CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1])*r3[0];
	CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2])*r3[1];
	CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0])*r3[2];

	Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

	return Volume;

}

