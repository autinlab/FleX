
class Actine : public Scene
{
public:

	Actine(const char* name) : Scene(name) {}
	struct Branch
	{
		int npoints;
		std::vector<int> p_indices;
	};
	std::vector<Branch> all_branch;
	std::vector<Vec3> all_points;
	std::vector<int> ropes_mIndices;
	float main_scale = 1.0f / 100.0f;
	int iGroupCounter = 0;
	ofstream output;
	ofstream output_bin;
	
	virtual void CreatePersistenceRope(Rope& rope, int current, int persistence, float stiffness,
		int nfloat, float give, float D, float hardness)
	{
		for (int j = 1; j < persistence + 1; j++){
			//float r = Randf((D / 100.0f)*-1.0f, 0.0f);
			float r = Randf(-hardness, hardness);
			if (rope.mIndices.size()>j)
				CreateSpringInter(rope.mIndices[current - j], rope.mIndices[current], stiffness, give, (D*(float)j) + r);// *(float)(j - 1));
		}
	}

	virtual void MakeAdamGreatAgain(){
		for (int i = 0; i < all_branch.size(); i++){
			//data_curve is an array of float
			std::vector<Vec3> points;
			for (int j = 0; j < all_branch[i].npoints; j++){
				points.push_back(all_points[all_branch[i].p_indices[j]]);
			}
			float* data_curve = reinterpret_cast<float*>(points.data());
			int rope_phase = NvFlexMakePhase(iGroupCounter++, eNvFlexPhaseSelfCollide);
			Rope curve;

			CreateRopeFromData(curve, //rope
				main_scale, // scale
				1.0f, //stifness
				data_curve, //data
				0.0f, //length
				all_branch[i].npoints * 3,  //nfloat
				rope_phase,//phase
				0.0f,//spiral angle
				1.0f,//invmass
				0.0f,//give
				0,//extend_nb
				false,//extend
				false,//close
				g_params.radius*2.0f,
				0);
			curve.persistence = 0;
			printf("create Rope from data OK with %i points\n", curve.mIndices.size());
			g_ropes.push_back(curve);//instance
		}
	}

	virtual void BuildNetwork(){
		int begin = 0;
		int persistence = 5;
		float D = (7 / 1000.0f) / 2.0f;
		for (int i = 0; i < all_branch.size(); i++){
			//data_curve is an array of float
			Rope curve;
			curve.mIndices.push_back(all_branch[i].p_indices[0]);
			for (int j = 0; j < all_branch[i].npoints; j++){
				curve.mIndices.push_back(all_branch[i].p_indices[j]);
				//CreatePersistence(curve, begin, 0, 1.0f, 0, 0.0f, g_params.mRadius*2.0f);
				CreateSpringInter(all_branch[i].p_indices[j], all_branch[i].p_indices[j - 1], 1.0f);//or distance
				//CreatePersistence(curve, int current, int persistence, float stiffness, int nfloat, float give, float D)
				/*for (int p = 1; p < persistence + 1; p++){
				float r = Randf((D)*-1.0f, 0.0f);
				if (curve.mIndices.size()>p)
				CreateSpringInter(all_branch[i].p_indices[j - p], all_branch[i].p_indices[j], 1.0f);// , (D*(float)p) + r*(float)(p - 1));
				}*/
			}
			curve.persistence = persistence;
			printf("create Rope from data OK with %i points\n", curve.mIndices.size());
			g_ropes.push_back(curve);//instance
		}
	}

	virtual int checkIfExistInOtherBranch(int query_indice, int current){
		bool found = false;
		int branch_i = -1;
		for (int i = 0; i < all_branch.size(); i++)
		{
			if (i == current) continue;
			std::vector<int>::iterator iter = std::find(all_branch[i].p_indices.begin(), all_branch[i].p_indices.end(), query_indice); //used_indices.find(pindice);		  //already used ?
			if (iter != all_branch[i].p_indices.end())
			{
				found = true;
				branch_i = i;
				break;
			}
		}
		return branch_i;
	}

	// A BFS based function to check whether d is reachable from s.
	
	virtual void BuildParticleAndNetworkFromSimplifiedGraph(){
		//std::set<int> used_indices;
		ropes_mIndices.clear();
		std::vector<int> used_indices;
		std::vector<int> used_indices_gpart;
		std::set<int>::iterator it;

		int begin = 0;
		int persistence = 1;
		float D = g_params.radius;//distance between two points for collision
		int start = int(g_buffers->positions.size());
		float give = 0.0f;
		float stiffness = 1.0f;
		float r = 1.0f;//biased on the 1-3 spring
		int current = 0;
		//if closed do the last point ?
		int count = 0;
		int prev = 0;
		float hardness = 0.0f; //lead to D/1000.0f
		//this part doesnt necessary find the branching..
		for (int i = 0; i < all_branch.size(); i++)
		{
			int phase = NvFlexMakePhase(iGroupCounter++, eNvFlexPhaseSelfCollide);
			if ((i % 500) == 0) cout << " i " << i << endl;
			//if (count > 5) break;
			//if (i > 2050) break;
			Rope rope;
			//int phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);
			//printf("add a point %i %f %f %f \n",i/3,data[i]/1000.0f, data[i+1]/1000.0f, data[i+2]/1000.0f);
			int begin = int(g_buffers->positions.size());
			int previous = 0;
			//cout << i << endl;
			//cout << all_branch[i].npoints << endl;
			//cout << all_branch[i].p_indices.size() << endl;
			//cout << all_branch[i].p_indices[0] << endl;
			if (all_branch[i].npoints < 4) continue;
			Vec3 start = all_points[all_branch[i].p_indices[0]];
			for (int j = 0; j < all_branch[i].npoints - 1; j++){// all_branch[i].npoints; j++){
				int pindice = all_branch[i].p_indices[j]; //indice in all_pos
				int gindice = int(g_buffers->positions.size());    //indice in particle
				float dist = Length(all_points[pindice] * main_scale - start*main_scale);

				if ((j > 0) && (dist < g_params.radius)){
					continue;
				}
				else
				{
					start = all_points[pindice];
				}

				//this step is really long,here should use it only for first and last point
				if (j == 0){
					std::vector<int>::iterator iter = std::find(used_indices.begin(), used_indices.end(), pindice); //used_indices.find(pindice);		  //already used ?
					if (iter != used_indices.end())
					{
						//indice found can reuse it
						int k = std::distance(used_indices.begin(), iter);
						//cout << " k is " << k << " size is " << used_indices_gpart.size() << "query was " << pindice << " check " << used_indices[k] << endl;
						gindice = used_indices_gpart.at(k);
						count++;
					}
					else
					{
						used_indices.push_back(pindice);
						used_indices_gpart.push_back(gindice);
						g_buffers->positions.push_back(Vec4(all_points[pindice].x*main_scale, all_points[pindice].y*main_scale, all_points[pindice].z*main_scale, 1.0f));
						g_buffers->velocities.push_back(0.0f);
						g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));					
					}
				}
				else
				{
					used_indices.push_back(pindice);
					used_indices_gpart.push_back(gindice);
					g_buffers->positions.push_back(Vec4(all_points[pindice].x*main_scale, all_points[pindice].y*main_scale, all_points[pindice].z*main_scale, 1.0f));
					g_buffers->velocities.push_back(0.0f);
					g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
				}

				current = int(rope.mIndices.size());
				rope.mIndices.push_back(gindice);
				ropes_mIndices.push_back(gindice);
				CreatePersistenceRope(rope, current, persistence, stiffness, 0, give, D / 2.0f, hardness);
				//CreateSpringInter(rope.mIndices[current - j], rope.mIndices[current], stiffness, give);
				prev = begin;//int(g_buffers->positions.size())-1;
				current = begin;
				//
			}
			//add the last element
			int pindice = all_branch[i].p_indices[all_branch[i].npoints - 1]; //indice in all_pos
			int gindice = int(g_buffers->positions.size());    //indice in particle
			std::vector<int>::iterator iter = std::find(used_indices.begin(), used_indices.end(), pindice); //used_indices.find(pindice);		  //already used ?
			if (iter != used_indices.end())
			{
				//indice found can reuse it
				int k = std::distance(used_indices.begin(), iter);
				//cout << " k is " << k << " size is " << used_indices_gpart.size() << "query was " << pindice << " check " << used_indices[k] << endl;
				gindice = used_indices_gpart.at(k);
				count++;
			}
			else
			{
				used_indices.push_back(pindice);
				used_indices_gpart.push_back(gindice);
				g_buffers->positions.push_back(Vec4(all_points[pindice].x*main_scale, all_points[pindice].y*main_scale, all_points[pindice].z*main_scale, 1.0f));
				g_buffers->velocities.push_back(0.0f);
				g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
			}
			current = int(rope.mIndices.size());
			rope.mIndices.push_back(gindice);
			ropes_mIndices.push_back(gindice);
			CreatePersistenceRope(rope, current, persistence, stiffness, 0, give, D / 2.0f, hardness);

			float dist = Length(Vec3(g_buffers->positions[rope.mIndices[0]]) - Vec3(g_buffers->positions[rope.mIndices[rope.mIndices.size() - 1]]));
			float a = (D / 2.0f)*(float)rope.mIndices.size();// all_branch[i].npoints);
			float b = a;
			if (dist > a) {
				b = dist;
			}
			CreateSpringInter(rope.mIndices[0], rope.mIndices[rope.mIndices.size() - 1], stiffness, give, b);
			rope.persistence = persistence;
			g_ropes.push_back(rope);//instance
		}
	}


	virtual void ParseVTK(){
		/*
		write a cahce for all points/all lines
		parse branche.vtk
		the points are defined first after
		DATASET POLYDATA
		POINTS n dataType
		p 0x p 0y p 0z
		p 1x p 1y p 1z
		then the line indices
		LINES n size
		numPoints 0 , i 0 ,j 0 ,k 0 , ...
		numPoints 1 , i 1 ,j 1 ,k 1 , ...
		...
		numPoints n-1 , i n-1 ,j n-1 ,k n-1 , ...
		then some attributes
		SCALARS dataName dataType numComp
		LOOKUP_ TABLE tableName
		s 0
		s 1
		...
		s n-1
		*/
		//int phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);
		string datapath = "..\\..\\data\\cellpack\\";
		ifstream source(datapath + "branches.vtk");
		std::cout << "read in " << datapath + "branches.vtk" << endl;
		for (std::string line; std::getline(source, line);)   //read stream line by line
		{
			std::istringstream in(line);      //make a stream for the line itself

			std::string type;
			in >> type;                  //and read the first whitespace-separated token

			if (type == "POINTS")       //and check its value
			{
				//POINTS 828183 float
				int N;
				in >> N >> type;       //now read the whitespace-separated floats
				std::cout << " found " << N << " points " << endl;
				//should be as many lines as N points
				for (int i = 0; i < N; i++){
					std::getline(source, line);
					std::istringstream in(line);
					float x, y, z;
					in >> x >> y >> z;
					all_points.push_back(Vec3(x*main_scale, y*main_scale, z*main_scale));
					//std::cout << "i " << i << endl;
					//g_buffers->positions.push_back(Vec4(x*main_scale, y*main_scale, z*main_scale, 1.0f));
					//g_buffers->velocities.push_back(0.0f);
					//g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
				}
				std::cout << " parsed " << all_points.size() << " points " << endl;
			}
			else if (type == "LINES")
			{
				int N;
				int size;
				in >> N >> size;       //now read the whitespace-separated floats
				std::cout << " found " << N << " lines " << endl;
				for (int i = 0; i < N; i++){
					Branch b = Branch();
					std::getline(source, line);
					std::istringstream in(line);
					int np;
					in >> np;
					b.npoints = np;
					for (int j = 0; j < np; j++){
						int p_index;
						in >> p_index;
						b.p_indices.push_back(p_index);
					}
					all_branch.push_back(b);
				}
				std::cout << " parsed " << all_branch.size() << " lines " << endl;
				//LINES n size  LINES 26500 856121
				//numPoints 0, i 0, j 0, k 0, ...
			}
			else {
				continue;
			}
		}
	}

	virtual void readCache(std::string filename){
		std::ifstream source(filename, std::ios::binary);
		int N;
		source.read(reinterpret_cast<char*>(&N), sizeof(int));
		std::cout << "read in npoints " << N << endl;
		g_buffers->positions.resize(N);
		g_buffers->velocities.resize(N);
		g_buffers->phases.resize(N);
		source.read(reinterpret_cast<char*>(&g_buffers->positions[0]), sizeof(float) * 4 * N);
		source.read(reinterpret_cast<char*>(&g_buffers->velocities[0]), sizeof(float) * 3 * N);
		source.read(reinterpret_cast<char*>(&g_buffers->phases[0]), sizeof(int) * N);
		int Nspring;
		std::cout << "read in nspring " << Nspring << endl;
		source.read(reinterpret_cast<char*>(&Nspring), sizeof(int));
		g_buffers->springIndices.resize(Nspring * 2);
		g_buffers->springLengths.resize(Nspring);
		g_buffers->springStiffness.resize(Nspring);
		source.read(reinterpret_cast<char*>(&g_buffers->springIndices[0]), sizeof(int) * 2 * Nspring);
		source.read(reinterpret_cast<char*>(&g_buffers->springLengths[0]), sizeof(float) * Nspring);
		source.read(reinterpret_cast<char*>(&g_buffers->springStiffness[0]), sizeof(float) * Nspring);
		//the rope data
		int Nrope;
		source.read(reinterpret_cast<char*>(&Nrope), sizeof(int));
		std::cout << "read in n rope " << Nrope << endl;
		g_ropes.clear();
		g_ropes.resize(Nrope);
		ropes_mIndices.clear();
		for (int i = 0; i < Nrope; i++){
			//write Npoint,follow by points
			int Npoints;
			source.read(reinterpret_cast<char*>(&Npoints), sizeof(int));
			Rope r;
			r.mIndices.resize(Npoints);
			source.read(reinterpret_cast<char*>(&r.mIndices[0]), sizeof(int) * Npoints);
			g_ropes[i] = r;
			for (int j = 0; j < Npoints; j++) {
				ropes_mIndices.push_back(r.mIndices[j]);
			}
		}
		source.close();
		output_bin.open((filename + "_graph.txt").c_str(), ios::out);

		//int Nspring = int(g_buffers->springLengths.size());
		int j = 0;
		for (int i = 0; i < Nspring; i++){
			output_bin << g_buffers->springIndices[j] << " " << g_buffers->springIndices[j + 1] << endl;
			j += 2;
		}

		//save txt of g_position
		output_bin.close();
		output_bin.open((filename + "_points.txt").c_str(), ios::out);
		for (int i = 0; i < N; i++){
			output_bin << g_buffers->positions[i].x << " " << g_buffers->positions[i].y << " " << g_buffers->positions[i].z << " " << g_params.radius / 2.0f << endl;
		}
		output_bin.close();
	}

	virtual void writeCache(std::string filename){
		ofstream output_bin;
		output_bin.open(filename.c_str(), ios::out | ios::app | ios::binary);
		//first Npoints
		int N = int(g_buffers->positions.size());
		output_bin.write((char *)&N, sizeof(int));
		//write all the points
		output_bin.write((char *)&g_buffers->positions[0], sizeof(float) * 4 * N);
		output_bin.write((char *)&g_buffers->velocities[0], sizeof(float) * 3 * N);
		output_bin.write((char *)&g_buffers->phases[0], sizeof(int) * N);
		//write the springs
		int Nspring = int(g_buffers->springLengths.size());
		output_bin.write((char *)&Nspring, sizeof(int));
		output_bin.write((char *)&g_buffers->springIndices[0], sizeof(int) * 2 * Nspring);
		output_bin.write((char *)&g_buffers->springLengths[0], sizeof(float) *  Nspring);
		output_bin.write((char *)&g_buffers->springStiffness[0], sizeof(float) *  Nspring);
		//the rope data
		int Nrope = int(g_ropes.size());
		output_bin.write((char *)&Nrope, sizeof(int));
		for (int i = 0; i < Nrope; i++){
			//write Npoint,follow by points
			int Npoints = int(g_ropes[i].mIndices.size());
			output_bin.write((char *)&Npoints, sizeof(int));
			output_bin.write((char *)&g_ropes[i].mIndices[0], sizeof(int) * Npoints);
		}
		output_bin.close();
		//ofstream output_bin;
		output_bin.open((filename + "_graph.txt").c_str(), ios::out);

		//int Nspring = int(g_buffers->springLengths.size());
		int j = 0;
		for (int i = 0; i < Nspring; i++){
			output_bin << g_buffers->springIndices[j] << " " << g_buffers->springIndices[j + 1] << endl;
			j += 2;
		}

		//save txt of g_position
		output_bin.close();
		output_bin.open((filename + "_points.txt").c_str(), ios::out);
		for (int i = 0; i < N; i++){
			output_bin << g_buffers->positions[i].x << " " << g_buffers->positions[i].y << " " << g_buffers->positions[i].z << " " << g_params.radius / 2.0f << endl;
		}
		output_bin.close();

	}

	virtual void readCacheRaw(std::string filename){
		std::ifstream source(filename, std::ios::binary);
		int N;
		source.read(reinterpret_cast<char*>(&N), sizeof(int));
		std::cout << "read Points in " << N << endl;
		all_points.clear();
		all_points.resize(N);
		source.read(reinterpret_cast<char*>(&all_points[0]), sizeof(float) * 3 * N);
		//the branch data
		int Nbranch;
		source.read(reinterpret_cast<char*>(&Nbranch), sizeof(int));
		std::cout << "read Branhc in " << Nbranch << endl;
		all_branch.clear();
		all_branch.resize(Nbranch);
		for (int i = 0; i < Nbranch; i++){
			//write Npoint,follow by points
			int Npoints;
			source.read(reinterpret_cast<char*>(&Npoints), sizeof(int));
			if (i == 0)std::cout << "Branch " << i << " " << Npoints << endl;
			Branch r = Branch();
			r.p_indices.clear();
			r.p_indices.resize(Npoints);
			r.npoints = Npoints;
			source.read(reinterpret_cast<char*>(&r.p_indices[0]), sizeof(int) * Npoints);
			all_branch[i] = r;
			if (i == 0)std::cout << "Branch pt 0 " << r.p_indices[0] << endl;
		}
		source.close();
	}

	virtual void writeCacheRaw(std::string filename){
		ofstream output_bin;
		output_bin.open(filename.c_str(), ios::out | ios::app | ios::binary);
		//first Npoints
		int N = int(all_points.size());
		output_bin.write((char *)&N, sizeof(int));
		//write all the points
		output_bin.write((char *)&all_points[0], sizeof(float) * 3 * N);
		int Nbranch = int(all_branch.size());
		output_bin.write((char *)&Nbranch, sizeof(int));
		for (int i = 0; i < Nbranch; i++){
			//write Npoint,follow by points
			int Npoints = int(all_branch[i].npoints);
			output_bin.write((char *)&Npoints, sizeof(int));
			output_bin.write((char *)&all_branch[i].p_indices[0], sizeof(int) * Npoints);
		}
		output_bin.close();
	}

	virtual void Initialize()
	{
		//unit in file 1 unit=1angstrom?? or 1 unit = 10nm
		main_scale = 1.0f / 100.0f;
		//actin filaments (13.5 μm) persistance length 13000/7
		float actine_radius = 27.6f / 2.0f;// 7.0f; thast in nm rise is 27.6A
		float beads_radius = actine_radius*main_scale;//5.0f*main_scale;//0.08f;//5.0f*main_scale;//smallest sphere size
		g_params.radius = beads_radius * 2.0f;

		//check if file exist
		std::string datapath = "..\\..\\data\\";
		std::ifstream source(datapath + "cache.bin", std::ios::binary);
		bool force = false;
		if ((source.is_open()) && (!force)) {
			std::cout << "read in " << datapath + "cache.bin" << endl;
			source.close();
			readCache(datapath + "cache.bin");
		}
		else {
			std::ifstream source_raw(datapath + "cache_raw_new.bin", std::ios::binary);
			if (source_raw.is_open()) {
				source_raw.close();
				readCacheRaw(datapath + "cache_raw_new.bin");
			}
			else {
				ParseVTK();
				writeCacheRaw(datapath + "cache_raw.bin");
			}
			BuildParticleAndNetworkFromSimplifiedGraph();
			writeCache(datapath + "cache.bin");
		}

		g_numExtraParticles = 1024 * 1024;
		g_params.gravity[1] = 0.0f;//-9.f;
		g_params.damping = 1.0f;// 3.0f;
		g_params.radius = beads_radius * 2;

		g_windStrength = 0.0f;
		g_windFrequency = 0.0f;

		g_params.numIterations = 5;

		//g_params.fluid = false;
		g_params.vorticityConfinement = 0.0f;
		//g_params.mFluidRestDistance = g_params.mRadius*0.5;
		g_params.anisotropyScale = 2.5f / g_params.radius;
		g_params.smoothing = 0.5f;
		g_params.relaxationFactor = 1.f;
		g_params.restitution = 0.0f;
		g_params.collisionDistance = 0.01f;

		g_params.dynamicFriction = 0.25f;
		g_params.viscosity = 1.5f;
		g_params.cohesion = 1.1f;
		g_params.adhesion = 1.0f;
		g_params.surfaceTension = 1.0f;


		g_lightDistance *= 2.0f;

		// draw options		
		g_drawEllipsoids = true;
		g_params.maxSpeed = 100.0f;


		g_wireframe = true;
		g_pointScale = 1.0f;
		//        g_blur = 2.0f;
		g_pause = true;
		g_warmup = false;
		g_drawMesh = false;
		g_drawRopes = false;

		g_params.dynamicFriction = 0.4f;
		g_params.dissipation = 0.0f;

		g_params.particleCollisionMargin = g_params.radius*0.05f;
		g_params.drag = 0.0f;
		
		// better convergence with global relaxation factor
		g_params.relaxationMode = eNvFlexRelaxationGlobal;
		g_params.relaxationFactor = 0.25f;

		g_windStrength = 0.0f;

		g_numSubsteps = 5;	

		// draw options		
		g_drawPoints = true;
		g_drawSprings = false;
		g_drawCloth = false;
	}
};