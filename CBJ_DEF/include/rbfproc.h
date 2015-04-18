#ifndef DEFPROC_H
#define DEFPROC_H

#include "../../zhongpeinan/include/rbf_tool.hpp"
#include "../../zhongpeinan/include/spatial_range.hpp"
#include "../../zhongpeinan/include/mesh_info.hpp"
#include "../../zhongpeinan/include/mesh_deform.hpp"

#include "../../zhongpeinan/include/preprocess.hpp"
#include "SU2_DEF.hpp"

#include <string>
#include <array>

class DeformProcedure
{
private:
	typedef CGridMove_Rbf GridMove;

	std::string              config_filename;
	CConfig*                 config;
	CGeometry*               geometry;
	int                      nDim;
	int                      curstep;
	int                      maxstep;
	GridMove*                grdmv;

	const RadialBasisFunction*     rbf;
	const SpatialRange*            sparg;

	enum { CPOINT, CPOINT_NEW, STRE, FORC, FORTAG, MVPT, MESH_OUT, N_MFILE };

	std::array<std::string, N_MFILE> message;

	static inline
		void make_tag(const std::string& tag_name, const std::string& content = "NONE")
	{
		std::ofstream tag(tag_name, std::ios::out);
		tag << content << std::endl;
		tag.close();
	}

	static inline
		void wait_for(const std::string& target)
	{
		std::cout << "Wait the file: " << target << std::endl;
		std::ifstream f(target, std::ios::in);
		while (!f.good())
		{
			f.clear(); f.close();
			f.open(target, std::ios::in);
		}
		std::cout << target << " exists." << std::endl;
	}

public:
	DeformProcedure();
	DeformProcedure(const std::string& aConfig_file);
	~DeformProcedure() { /*... delete grdmv; ...*/ }
	int dim() const { return nDim; }
	void set_maxstep(int n) { maxstep = n; }
	int get_maxstep() const { return maxstep; }
	int get_current_step() const { return curstep; }

	CGeometry* get_geometry() { return geometry; }
	CConfig*   get_config() { return config; }

	void preprocess(const RadialBasisFunction& aRbf, const SpatialRange& range);
	void set_rbf(const RadialBasisFunction& aRbf, const SpatialRange& range)
	{
		rbf = &aRbf;
		sparg = &range;
	}
	void set_default_message_file();

	bool deform();

	const string& get_new_mesh_file() const { return message[MESH_OUT]; }
	void set_new_mesh_file() { geometry->SetMeshFile(config, message[MESH_OUT]); }

};

#endif // DEFPROC_H
