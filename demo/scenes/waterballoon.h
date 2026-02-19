
class WaterBalloon : public Scene
{
public:
	float main_scale = 0.001f;
	bool write_output = false;
	bool force_write_frame = false;
	bool unique_name = false;
	float vlength_scale = 1.73205f;
	float vlength = 0.5f;
	float stiffness;


	// we need different length for penta and hexa. p6 and p5 neighboors length.
	WaterBalloon(const char* name) : Scene(name) {}

	virtual ~WaterBalloon()
	{
		for (size_t i = 0; i < mCloths.size(); ++i)
			NvFlexExtDestroyTearingCloth(mCloths[i].asset);
	}


	void AddSimpleInflatable(const Mesh* mesh, float overPressure, int phase)
	{
		const int numParticles = int(mesh->m_positions.size());
		const int maxParticles = numParticles * 2;
		std::cout << " create a cloth mesh using the global positions / indices " << numParticles << endl;

		Balloon balloon;
		balloon.particleOffset = g_buffers->positions.size();
		balloon.triangleOffset = g_buffers->triangles.size();
		// after you know numParticles and numTriangles
		balloon.vertToTris.clear();
		balloon.vertToTris.resize(numParticles);
		
		Vec3 mean(0, 0, 0);
		for (size_t i = 0; i < mesh->GetNumVertices(); ++i)
			mean += Vec3(mesh->m_positions[i]);
		mean /= float(mesh->GetNumVertices());

		// --- covariance matrix ---
		float C[3][3] = { {0,0,0},{0,0,0},{0,0,0} };

		const int startVertex = g_buffers->positions.size();
		const int startTriangle = g_buffers->triangles.size();

		// add mesh to system
		for (size_t i = 0; i < mesh->GetNumVertices(); ++i)
		{
			const Vec3 p3 = Vec3(mesh->m_positions[i]);

			g_buffers->positions.push_back(Vec4(p3.x, p3.y, p3.z, 1.0f));
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(phase);

			Vec3 p = Vec3(p3) - mean;
			C[0][0] += p.x*p.x; C[0][1] += p.x*p.y; C[0][2] += p.x*p.z;
			C[1][0] += p.y*p.x; C[1][1] += p.y*p.y; C[1][2] += p.y*p.z;
			C[2][0] += p.z*p.x; C[2][1] += p.z*p.y; C[2][2] += p.z*p.z;
		}
		// --- power iteration to get dominant eigenvector ---
		Vec3 v(1, 0, 0); // initial guess
		for (int it = 0; it < 20; ++it)
		{
			Vec3 v2(
				C[0][0] * v.x + C[0][1] * v.y + C[0][2] * v.z,
				C[1][0] * v.x + C[1][1] * v.y + C[1][2] * v.z,
				C[2][0] * v.x + C[2][1] * v.y + C[2][2] * v.z
			);
			float l = Length(v2);
			if (l < 1e-10f) break;
			v = v2 / l;
		}
		balloon.coneAxis = Normalize(v);
		int tLocal = 0;
		std::cout << " m_indices " << mesh->m_indices.size() << " " << mesh->GetNumFaces() << endl;

		int triOffset = g_buffers->triangles.size();
		int triCount = mesh->GetNumFaces();

		g_buffers->inflatableTriOffsets.push_back(triOffset / 3);
		g_buffers->inflatableTriCounts.push_back(mesh->GetNumFaces());
		g_buffers->inflatablePressures.push_back(overPressure);
		
		int tlocal = 0;

		for (size_t i = 0; i < mesh->m_indices.size(); i += 3, tlocal++)
		{
			int a = mesh->m_indices[i + 0];
			int b = mesh->m_indices[i + 1];
			int c = mesh->m_indices[i + 2];

			Vec3 n = -Normalize(Cross(mesh->m_positions[b] - mesh->m_positions[a], mesh->m_positions[c] - mesh->m_positions[a]));
			g_buffers->triangleNormals.push_back(n);

			g_buffers->triangles.push_back(a + startVertex);
			g_buffers->triangles.push_back(b + startVertex);
			g_buffers->triangles.push_back(c + startVertex);

			balloon.vertToTris[a].push_back(tLocal);
			balloon.vertToTris[b].push_back(tLocal);
			balloon.vertToTris[c].push_back(tLocal);
		}
		// why not use NvFlexExtCreateClothFromMesh here ? 

		// create a cloth mesh using the global positions / indices
		// 
		//ClothMesh* cloth = new ClothMesh(&g_buffers->positions[0], g_buffers->positions.size(), &g_buffers->triangles[triOffset], triCount * 3, 0.8f, 1.0f);
		//cloth->particleOffset = startVertex;
		//cloth->triangleOffset = startTriangle;

		//for (size_t i = 0; i < cloth->mConstraintIndices.size(); ++i)
		//	g_buffers->springIndices.push_back(cloth->mConstraintIndices[i]);

		//for (size_t i = 0; i < cloth->mConstraintCoefficients.size(); ++i)
		//	g_buffers->springStiffness.push_back(cloth->mConstraintCoefficients[i]);

		//for (size_t i = 0; i < cloth->mConstraintRestLengths.size(); ++i)
		//	g_buffers->springLengths.push_back(cloth->mConstraintRestLengths[i]);
		NvFlexExtAsset* cloth = NvFlexExtCreateClothFromMesh((float*)&g_buffers->positions[startVertex], startVertex, (int*)&mesh->m_indices[0], mesh->GetNumFaces(), 1.0f, 1.0f, 1.0f, 0.0f, 0.1f);
		balloon.asset = cloth;
		balloon.tearable = false;
		std::cout << " cloth " << cloth << endl;
		mCloths.push_back(balloon);

		// add inflatable params
		g_buffers->inflatableVolumes.push_back(cloth->inflatableVolume);
		g_buffers->inflatableCoefficients.push_back(cloth->inflatableStiffness);
	}

	void AddInflatable(const Mesh* mesh, float overPressure, float invMass, int phase)
	{
		// create a cloth mesh using the global positions / indices
		
		const int numParticles = int(mesh->m_positions.size());
		const int maxParticles = numParticles * 2;
		std::cout << " create a cloth mesh using the global positions / indices " << numParticles << endl;

		Balloon balloon;
		balloon.particleOffset = g_buffers->positions.size();
		balloon.triangleOffset = g_buffers->triangles.size();
		balloon.splitThreshold = 4.0f;
		// after you know numParticles and numTriangles
		balloon.vertToTris.clear();
		balloon.vertToTris.resize(numParticles);
		Vec3 mean(0, 0, 0);
		for (size_t i = 0; i < mesh->GetNumVertices(); ++i)
			mean += Vec3(mesh->m_positions[i]);
		mean /= float(mesh->GetNumVertices());

		// --- covariance matrix ---
		float C[3][3] = { {0,0,0},{0,0,0},{0,0,0} };

		// add particles to system
		for (size_t i = 0; i < mesh->GetNumVertices(); ++i)
		{
			const Vec3 p3 = Vec3(mesh->m_positions[i]);

			g_buffers->positions.push_back(Vec4(p3.x, p3.y, p3.z, invMass));
			g_buffers->restPositions.push_back(Vec4(p3.x, p3.y, p3.z, invMass));

			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(phase);

			Vec3 p = Vec3(p3) - mean;
			C[0][0] += p.x*p.x; C[0][1] += p.x*p.y; C[0][2] += p.x*p.z;
			C[1][0] += p.y*p.x; C[1][1] += p.y*p.y; C[1][2] += p.y*p.z;
			C[2][0] += p.z*p.x; C[2][1] += p.z*p.y; C[2][2] += p.z*p.z;
		}
		// --- power iteration to get dominant eigenvector ---
		Vec3 v(1, 0, 0); // initial guess
		for (int it = 0; it < 20; ++it)
		{
			Vec3 v2(
				C[0][0] * v.x + C[0][1] * v.y + C[0][2] * v.z,
				C[1][0] * v.x + C[1][1] * v.y + C[1][2] * v.z,
				C[2][0] * v.x + C[2][1] * v.y + C[2][2] * v.z
			);
			float l = Length(v2);
			if (l < 1e-10f) break;
			v = v2 / l;
		}
		balloon.coneAxis = Normalize(v);
		int tLocal = 0;
		std::cout << " m_indices " << mesh->m_indices.size() << " " << mesh->GetNumFaces() << endl;
		for (size_t i = 0; i < mesh->m_indices.size(); i += 3, ++tLocal)
		{
			int a = mesh->m_indices[i + 0];
			int b = mesh->m_indices[i + 1];
			int c = mesh->m_indices[i + 2];

			Vec3 n = -Normalize(Cross(mesh->m_positions[b] - mesh->m_positions[a], mesh->m_positions[c] - mesh->m_positions[a]));
			g_buffers->triangleNormals.push_back(n);

			g_buffers->triangles.push_back(a + balloon.particleOffset);
			g_buffers->triangles.push_back(b + balloon.particleOffset);
			g_buffers->triangles.push_back(c + balloon.particleOffset);


			balloon.vertToTris[a].push_back(tLocal);
			balloon.vertToTris[b].push_back(tLocal);
			balloon.vertToTris[c].push_back(tLocal);
		}

		// create tearing asset or NvFlexExtCreateClothFromMesh
		NvFlexExtAsset* cloth = NvFlexExtCreateTearingClothFromMesh(
			(float*)&g_buffers->positions[balloon.particleOffset],
			numParticles, maxParticles,
			(int*)&mesh->m_indices[0], mesh->GetNumFaces(),
			1.0f, 1.0f, 0.0f);
		//(const float* particles, int numVertices, const int* indices, int numTriangles, float stretchStiffness, float bendStiffness, float tetherStiffness, float tetherGive, float pressure)
		// NvFlexExtAsset* cloth = NvFlexExtCreateClothFromMesh((float*)&g_buffers->positions[balloon.particleOffset], numParticles, (int*)&mesh->m_indices[0], mesh->GetNumFaces(), 1.0f, 1.0f, 1.0f, 0.0f,0.1f);
		balloon.asset = cloth;
		balloon.tearable = true;
		std::cout << " cloth " << cloth << endl;
		mCloths.push_back(balloon);
		std::cout << " mCloths " << mCloths.size() <<" "<< balloon.asset->numParticles << endl;
	}
	
	inline void AddOneParticule(Vec3 p3, int phase, float invMass) {
		g_buffers->positions.push_back(Vec4(p3.x, p3.y, p3.z, invMass));
		g_buffers->restPositions.push_back(Vec4(p3.x, p3.y, p3.z, invMass));
		g_buffers->velocities.push_back(0.0f);
		g_buffers->phases.push_back(phase);
	}
	
	inline void AddSpring(int i, int j, float stiffness)
	{
		if (i == j) return;

		const Vec3 pi = Vec3(g_buffers->positions[i]);
		const Vec3 pj = Vec3(g_buffers->positions[j]);
		const float rest = Length(pj - pi);

		g_buffers->springIndices.push_back(i);
		g_buffers->springIndices.push_back(j);
		g_buffers->springLengths.push_back(rest);
		g_buffers->springStiffness.push_back(stiffness);
	}

	void AddInflatableCustom(const Mesh* mesh, float overPressure, float invMass, int phase)
	{
		const int startParticle = int(g_buffers->positions.size());

		const int numFaces = int(mesh->GetNumFaces());
		const int numParticles = numFaces * 3;

		std::cout << "AddInflatableCustom: " << numParticles << " particles (3 per tri)\n";

		// Map: original vertex -> all created particle indices that represent it
		std::vector<std::vector<int>> origVertToParticles(mesh->GetNumVertices());
		origVertToParticles.reserve(mesh->GetNumVertices());

		// Optional: store per-triangle created particle indices (for easy springs + faces)
		std::vector<int> triCornerParticle; // length = numFaces*3
		triCornerParticle.reserve(numParticles);

		int tLocal = 0;

		// -----------------------
		// 1) Create particles
		// -----------------------
		for (size_t idx = 0; idx < mesh->m_indices.size(); idx += 3, ++tLocal)
		{
			int a = mesh->m_indices[idx + 0];
			int b = mesh->m_indices[idx + 1];
			int c = mesh->m_indices[idx + 2];

			Vec3 va = mesh->m_positions[a];
			Vec3 vb = mesh->m_positions[b];
			Vec3 vc = mesh->m_positions[c];

			Vec3 vcenter = (va + vb + vc) / 3.0f;

			// pull corners slightly toward tri center (optional; keep if it helps)
			Vec3 ca = va + Normalize(vcenter - va) * g_params.radius;
			Vec3 cb = vb + Normalize(vcenter - vb) * g_params.radius;
			Vec3 cc = vc + Normalize(vcenter - vc) * g_params.radius;

			const int pa = startParticle + int(triCornerParticle.size()) + 0;
			const int pb = startParticle + int(triCornerParticle.size()) + 1;
			const int pc = startParticle + int(triCornerParticle.size()) + 2;

			AddOneParticule(ca, phase, invMass);
			AddOneParticule(cb, phase, invMass);
			AddOneParticule(cc, phase, invMass);

			triCornerParticle.push_back(pa);
			triCornerParticle.push_back(pb);
			triCornerParticle.push_back(pc);

			// Record duplicates for glue springs
			origVertToParticles[a].push_back(pa);
			origVertToParticles[b].push_back(pb);
			origVertToParticles[c].push_back(pc);
		}

		// -----------------------
		// 2) Build triangle index list for rendering / pressure / etc.
		//    (these triangles refer to the NEW particles, not original vertices)
		// -----------------------
		for (int t = 0; t < numFaces; ++t)
		{
			int pa = triCornerParticle[t * 3 + 0];
			int pb = triCornerParticle[t * 3 + 1];
			int pc = triCornerParticle[t * 3 + 2];

			//g_buffers->triangles.push_back(pa);
			//g_buffers->triangles.push_back(pb);
			//g_buffers->triangles.push_back(pc);

			// Triangle normal from particle positions (not mesh positions)
			Vec3 p0 = Vec3(g_buffers->positions[pa]);
			Vec3 p1 = Vec3(g_buffers->positions[pb]);
			Vec3 p2 = Vec3(g_buffers->positions[pc]);
			Vec3 n = -Normalize(Cross(p1 - p0, p2 - p0));
			//g_buffers->triangleNormals.push_back(n);
		}

		// -----------------------
		// 3) Springs
		// -----------------------
		const float kTri = 1.0f;  // TODO tune
		const float kGlue = 2.0f;  // TODO tune (stronger)

		// 3a) Intra-triangle springs (3 edges)
		for (int t = 0; t < numFaces; ++t)
		{
			int pa = triCornerParticle[t * 3 + 0];
			int pb = triCornerParticle[t * 3 + 1];
			int pc = triCornerParticle[t * 3 + 2];

			AddSpring(pa, pb, kTri);
			AddSpring(pb, pc, kTri);
			AddSpring(pc, pa, kTri);
		}

		// 3b) Glue springs: connect all duplicates belonging to same original vertex
		// Strong around vertices
		for (int v = 0; v < int(origVertToParticles.size()); ++v)
		{
			auto& lst = origVertToParticles[v];
			const int m = int(lst.size());
			if (m <= 1) continue;

			// Star connect to first (O(m)) instead of full clique (O(m^2))
			for (int j = 0; j < m; ++j)
			{
				int a = lst[j];
				int b = lst[(j + 1) % m];
				AddSpring(a, b, kGlue);
			}
			//int root = lst[0];
			//for (int j = 1; j < m; ++j)
			//	AddSpring(root, lst[j], kGlue);
		}
		// Done: you now have a “triangle soup” cloth where shared vertices are strongly glued.
	}

	static int CeilDiv(int a, int b) { return (a + b - 1) / b; }

	void Initialize()
	{
		mCloths.resize(0);
		
		vlength_scale = 1.73205f;
		vlength = 1.0f / vlength_scale;

		float minSize = 0.25f;
		float maxSize = 0.5f;
		float spacing = 4.0f;
		//float main_scale = 0.001f;
		// const int dim = 128;
		// convex rocks
		//for (int i = 0; i < 4; i++)
		//	for (int j = 0; j < 1; j++)
		//		AddRandomConvex(10, Vec3(i*maxSize*spacing, 0.0f, j*maxSize*spacing), minSize, maxSize, Vec3(0.0f, 1.0f, 0.0f), Randf(0.0f, k2Pi));
		// NvFlexDistanceFieldId sdf = CreateSDF(GetFilePathByPlatform("../../data/cellpack/boxZ.ply").c_str(), dim, 0.1f, 0.0f, 1.0f);
		// AddSDF(sdf, Vec3(0.0f, 0.0f, 0.0f), QuatFromAxisAngle(Vec3(0.0f, 0.0f, 0.0f), DegToRad(0.0f)), 1.0f);

		Mesh* box_mesh = ImportMesh(GetFilePathByPlatform("../../data/cellpack/boxZ.ply").c_str());
		// box_mesh->Normalize();
		box_mesh->Transform(ScaleMatrix(main_scale));
		NvFlexTriangleMeshId amesh = CreateTriangleMesh(box_mesh);
		AddTriangleMesh(amesh, Vec3(), Quat(), 1.0f);
		
		// hexagon-heaxgon distance is 51Ang * 0.001f
		// g_params.radius = (vlength * vlength_scale) / 2.0;
		float radius = 0.09f;// (vlength * vlength_scale) / 2.0; // 0.10f;
		int group = 0;

		g_numExtraParticles = 500000;
		g_numSubsteps = 3;

		g_params.radius = radius;
		g_params.dynamicFriction = 1.0f; //  0.125f;
		g_params.staticFriction = 1.0f;
		g_params.particleFriction = 1.0f;
		g_params.damping = 30.0f;
		g_params.dissipation = 0.0f;
		g_params.numIterations = 5;
		g_params.particleCollisionMargin = g_params.radius*0.05f;
		g_params.relaxationFactor = 1.0f;
		g_params.drag = 0.0f;
		g_params.smoothing = 1.f;
		g_params.maxSpeed = 0.5f*g_numSubsteps*radius / g_dt;
		g_params.gravity[1] *= 1.0f;
		g_params.collisionDistance = 0.01f;
		g_params.solidPressure = 0.0f;

		g_params.fluidRestDistance = radius*0.65f;
		g_params.viscosity = 0.0;
		g_params.adhesion = 0.0f;
		g_params.cohesion = 0.02f;


		// add inflatables
		std::vector<const char*> files = {
			"../../data/cellpack/rspi_0.ply",
			"../../data/cellpack/rspi_10.ply",
			"../../data/cellpack/rspi_20.ply",
			"../../data/cellpack/rspi_30.ply",
			"../../data/cellpack/rspi_40.ply",
			// "../../data/cellpack/rspi_50.ply",   // Problem
			"../../data/cellpack/rspi_60.ply",
			// "../../data/cellpack/rspi_70.ply",   // Problem
			"../../data/cellpack/rspi_80.ply",
			"../../data/cellpack/rspi_90.ply",
			"../../data/cellpack/rspi_97.ply",

			"../../data/cellpack/rspi_0_230.ply",
			"../../data/cellpack/rspi_10_230.ply",
			// "../../data/cellpack/rspi_20_230.ply",
			// "../../data/cellpack/rspi_50_230.ply",
			"../../data/cellpack/rspi_90_230.ply",
		};

		std::vector<Mesh*> inmesh;
		inmesh.reserve(files.size());

		for (const char* rel : files) {
			const std::string path = GetFilePathByPlatform(rel);
			inmesh.push_back(ImportMesh(path.c_str()));
		}

		// std::vector<Mesh*> listmesh;
		const int copiesPerMesh = 6;   // change to 2, 3, 4, 6, etc.
		const int meshCount = (int)inmesh.size();
		const int total = meshCount * copiesPerMesh;
		const int maxCols = 4;
		
		const int cols = std::min(maxCols, (int)std::ceil(std::sqrt((double)total)));              // never more than 4, and not > total
		const int rows = (total + cols - 1) / cols;             // ceil(total/cols)

		const float halfCols = (cols - 1) * 0.5f;
		const float halfRows = (rows - 1) * 0.5f;

		//const int cols = (int)std::ceil(std::sqrt((double)total));
		//const int rows = CeilDiv(total, cols);

		//const float halfCols = 0.5f * (cols - 1);
		//const float halfRows = 0.5f * (rows - 1);

		const float gspace = 1.5f;
		Point3 lower(0.0, 1.0, 0.0);
		for (int k = 0; k < total; ++k)
		{
			const int i = k / copiesPerMesh; // which mesh
			const int j = k % copiesPerMesh; // which copy of that mesh (if you need it)

			const int col = k % cols;
			const int row = k / cols;

			const float x = (col - halfCols) * gspace;
			const float y = (row - halfRows) * gspace + halfRows * gspace;
			const float z = 0.0f;
			std::cout << " k " << k << " mesh " << i << " " << inmesh.size() << endl;
			// IMPORTANT: you need a real copy/clone here (whatever your Mesh API is)
			Mesh* mesh = new Mesh(*inmesh[i]);  // or new Mesh(*inmesh[i]) etc.
			mesh->Transform(ScaleMatrix(main_scale));
			// mesh->Normalize();

			// base placement + grid offset in XZ
			Point3 p(lower.x + x, lower.y + y, lower.z + z);
			mesh->Transform(TranslationMatrix(p));
			// overPressure, invMass, group
			
			AddInflatable(mesh,
				1.0f, 0.25f,
				NvFlexMakePhase(group++,
					eNvFlexPhaseSelfCollide | eNvFlexPhaseSelfCollideFilter));
			delete mesh;  //?
			// else AddSimpleInflatable(mesh, 1.0f, NvFlexMakePhase(group++,
			// 	eNvFlexPhaseSelfCollide | eNvFlexPhaseSelfCollideFilter));
			// listmesh.push_back(mesh);
		}
		// add some extra manual particles
		
		/*Mesh* mesh = new Mesh(*inmesh[0]);
		mesh->Transform(ScaleMatrix(main_scale));
		// mesh->Normalize();

		// base placement + grid offset in XZ
		Point3 p(0.5f,0.5f,0.5f);
		mesh->Transform(TranslationMatrix(p));
		AddInflatableCustom(mesh, 1.0f, 1.0, NvFlexMakePhase(group++,
			eNvFlexPhaseSelfCollide | eNvFlexPhaseSelfCollideFilter));
		*/

		g_numSolidParticles = g_buffers->positions.size();
		g_numExtraParticles = g_buffers->positions.size();

		// fill inflatables with water
		//can we do it with SDF instead ? 
		//std::vector<Vec3> positions(10000);
		//int n = PoissonSample3D(0.45f, g_params.radius*0.42f, &positions[0], positions.size(), 10000);
		//int n = TightPack3D(0.45f, g_params.radius*0.42f, &positions[0], positions.size());
		mNumFluidParticles = 0;
		group++;
		for (size_t k = 0; k < total; ++k)
		{
			const int i = k / copiesPerMesh; // which mesh
			const int j = k % copiesPerMesh; // which copy of that mesh (if you need it)

			const int col = k % cols;
			const int row = k / cols;

			const float x = (col - halfCols) * gspace;
			const float y = (row - halfRows) * gspace + halfRows * gspace;
			const float z = 0.0f;

			Mesh* mesh = new Mesh(*inmesh[i]);
			// mesh->Normalize();
			mesh->Transform(ScaleMatrix(main_scale));
			// base placement + grid offset in XZ
			Vec3 p(lower.x + x, lower.y + y, lower.z + z);
			//mesh->Transform(TranslationMatrix(p));
			int vstart = int(g_buffers->positions.size());
			Vec3 mlower, mupper;
			mesh->GetBounds(mlower, mupper);
			// Vec3 center = (mlower + mupper) * 0.5f;
			// Vec3(center.x + x, center.y + y, center.z + z)
			// CreateParticleShape(const Mesh* srcMesh, Vec3 lower, Vec3 scale, float rotation, float spacing, Vec3 velocity, float invMass, bool rigid, float rigidStiffness, int phase, bool skin, float jitter = 0.005f, Vec3 skinOffset = 0.0f, float skinExpand = 0.0f, Vec4 color = Vec4(0.0f), float springStiffness = 0.0f)
			CreateParticleShape(mesh, mlower+p, 1.2f, 0.0f, g_params.radius*0.45f, Vec3(0.0f, 0.0f, 0.0f), 1.0f, false, 0.f, NvFlexMakePhase(group, eNvFlexPhaseSelfCollide | eNvFlexPhaseFluid), false, 0.0f);
			int vend = int(g_buffers->positions.size());
			mNumFluidParticles += (vend - vstart);
			//put a pasta inside ? 
			//Rope r;
			//float rlength = Length(mupper - mlower);
			//CreateRope(r, mlower + p, Vec3(0.0f, 1.0f, 0.0f), 0.25f, int(rlength / radius), rlength, NvFlexMakePhase(group++, eNvFlexPhaseSelfCollide));
			//g_ropes.push_back(r);

			/*
			const int vertStart = i*mesh->GetNumVertices();
			const int vertEnd = vertStart + mesh->GetNumVertices();

			const int phase = NvFlexMakePhase(group++, eNvFlexPhaseSelfCollide | eNvFlexPhaseFluid);

			Vec3 center;
			for (int v = vertStart; v < vertEnd; ++v)
				center += Vec3(g_buffers->positions[v]);

			center /= float(vertEnd - vertStart);
			// const Mesh* srcMesh, Vec3 lower, Vec3 scale, float rotation, float spacing, Vec3 velocity, float invMass, bool rigid, float rigidStiffness, int phase, bool skin, float jitter=0.005f, Vec3 skinOffset=0.0f, float skinExpand=0.0f, Vec4 color=Vec4(0.0f), float springStiffness=0.0f
			// why cant we use CreateParticleShape(mesh, Vec3(-2.0f + j*size, 3.0f + j*size, i*size), size, 0.0f, spacing, Vec3(0.0f, 0.0f, 0.0f), 1.0f, true, 1.f, phase, false, 0.0f); ?
			printf("%d, %d - %f %f %f\n", vertStart, vertEnd, center.x, center.y, center.z);

			for (int i = 0; i < n; ++i)
			{
				g_buffers->positions.push_back(Vec4(center + positions[i], 1.0f));
				g_buffers->restPositions.push_back(Vec4());
				g_buffers->velocities.push_back(0.0f);
				g_buffers->phases.push_back(phase);
			}
			mNumFluidParticles += n;
			*/
		}

		for (int i=0; i < inmesh.size(); i++)
			delete inmesh[i];

		g_drawPoints = false;
		g_drawEllipsoids = true;
		g_drawSprings = 0;
		g_drawCloth = false;
		g_warmup = false;
		g_pause = true;
	}

	void RebuildConstraints()
	{
		// update constraint data
		g_buffers->triangles.resize(0);
		g_buffers->springIndices.resize(0);
		g_buffers->springStiffness.resize(0);
		g_buffers->springLengths.resize(0);

		for (int c = 0; c < int(mCloths.size()); ++c)
		{
			Balloon& balloon = mCloths[c];
			balloon.vertToTris.clear();
			balloon.vertToTris.resize(balloon.asset->numParticles);
			int tLocal = 0;
			for (int i = 0; i < balloon.asset->numTriangles; ++i, tLocal++)
			{
				g_buffers->triangles.push_back(balloon.asset->triangleIndices[i * 3 + 0] + balloon.particleOffset);
				g_buffers->triangles.push_back(balloon.asset->triangleIndices[i * 3 + 1] + balloon.particleOffset);
				g_buffers->triangles.push_back(balloon.asset->triangleIndices[i * 3 + 2] + balloon.particleOffset);

				int a = balloon.asset->triangleIndices[i * 3 + 0];
				int b = balloon.asset->triangleIndices[i * 3 + 1];
				int c = balloon.asset->triangleIndices[i * 3 + 2];

				balloon.vertToTris[a].push_back(tLocal);
				balloon.vertToTris[b].push_back(tLocal);
				balloon.vertToTris[c].push_back(tLocal);
			}

			for (int i = 0; i < balloon.asset->numSprings * 2; ++i)
				g_buffers->springIndices.push_back(balloon.asset->springIndices[i] + balloon.particleOffset);


			for (int i = 0; i < balloon.asset->numSprings; ++i)
			{
				g_buffers->springStiffness.push_back(balloon.asset->springCoefficients[i]);
				g_buffers->springLengths.push_back(balloon.asset->springRestLengths[i]);
			}
		}
	}

	//use inflatable as well with inside volume instead of tearing water ballon ? 

	virtual void Sync()
	{
		// send new particle data to the GPU
		NvFlexSetRestParticles(g_solver, g_buffers->restPositions.buffer, NULL);

		// update solver
		NvFlexSetSprings(g_solver, g_buffers->springIndices.buffer, g_buffers->springLengths.buffer, g_buffers->springStiffness.buffer, g_buffers->springLengths.size());
		NvFlexSetDynamicTriangles(g_solver, g_buffers->triangles.buffer, g_buffers->triangleNormals.buffer, g_buffers->triangles.size() / 3);
		// NvFlexSetInflatables(g_solver, g_buffers->inflatableTriOffsets.buffer, g_buffers->inflatableTriCounts.buffer, g_buffers->inflatableVolumes.buffer, g_buffers->inflatablePressures.buffer, g_buffers->inflatableCoefficients.buffer, 1);
	}

	virtual void Update()
	{
		// temporarily restore the mouse particle's mass so that we can tear it
		if (g_mouseParticle != -1)
			g_buffers->positions[g_mouseParticle].w = g_mouseMass;

		// force larger radius for solid interactions to prevent interpenetration
		g_params.solidRestDistance = g_params.radius;

		// build new particle arrays
		std::vector<Vec4> newParticles;
		std::vector<Vec4> newParticlesRest;
		std::vector<Vec3> newVelocities;
		std::vector<int> newPhases;
		std::vector<Vec4> newNormals;

		for (int c = 0; c < int(mCloths.size()); ++c)
		{
			Balloon& balloon = mCloths[c];
			//redo the triangle 
			const int destOffset = newParticles.size();
			// std::cout << " destOffset " << destOffset  << endl;
			// std::cout << " balloon.asset->numParticles " << balloon.asset->numParticles << endl;
			// append existing particles
			for (int i = 0; i < balloon.asset->numParticles; ++i)
			{
				newParticles.push_back(g_buffers->positions[balloon.particleOffset + i]);
				newParticlesRest.push_back(g_buffers->restPositions[balloon.particleOffset + i]);
				newVelocities.push_back(g_buffers->velocities[balloon.particleOffset + i]);
				newPhases.push_back(g_buffers->phases[balloon.particleOffset + i]);
				newNormals.push_back(g_buffers->normals[balloon.particleOffset + i]);
			}

			// perform splitting
			const int maxCopies = 2048;
			const int maxEdits = 2048;

			NvFlexExtTearingParticleClone particleCopies[maxCopies];
			int numParticleCopies;

			NvFlexExtTearingMeshEdit triangleEdits[maxEdits];
			int numTriangleEdits;

			// update asset's copy of the particles
			memcpy(balloon.asset->particles, &g_buffers->positions[balloon.particleOffset], sizeof(Vec4)*balloon.asset->numParticles);

			// tear only if tearable
			if (balloon.tearable) NvFlexExtTearClothMesh(balloon.asset, balloon.splitThreshold, 1, particleCopies, &numParticleCopies, maxCopies, triangleEdits, &numTriangleEdits, maxEdits);

			// resize particle data arrays
			newParticles.resize(newParticles.size() + numParticleCopies);
			newParticlesRest.resize(newParticlesRest.size() + numParticleCopies);
			newVelocities.resize(newVelocities.size() + numParticleCopies);
			newPhases.resize(newPhases.size() + numParticleCopies);
			newNormals.resize(newNormals.size() + numParticleCopies);

			// copy particles
			for (int i = 0; i < numParticleCopies; ++i)
			{
				const int srcIndex = balloon.particleOffset + particleCopies[i].srcIndex;
				const int destIndex = destOffset + particleCopies[i].destIndex;

				newParticles[destIndex] = g_buffers->positions[srcIndex];
				newParticlesRest[destIndex] = g_buffers->restPositions[srcIndex];
				newVelocities[destIndex] = g_buffers->velocities[srcIndex];
				newPhases[destIndex] = g_buffers->phases[srcIndex];
				newNormals[destIndex]  = g_buffers->normals[srcIndex];
			}

			if (numParticleCopies)
			{
				// reduce split threshold for this balloon
				balloon.splitThreshold = 1.75f;
			}

			balloon.particleOffset = destOffset;
			balloon.asset->numParticles += numParticleCopies;
		}

		// append fluid particles
		const int fluidStart = g_numSolidParticles;
		const int fluidEnd = fluidStart + mNumFluidParticles;

		g_numSolidParticles = newParticles.size();

		for (int i = fluidStart; i < fluidEnd; ++i)
		{
			newParticles.push_back(g_buffers->positions[i]);
			newParticlesRest.push_back(Vec4());
			newVelocities.push_back(g_buffers->velocities[i]);
			newPhases.push_back(g_buffers->phases[i]);
			newNormals.push_back(g_buffers->normals[i]);
		}

		g_buffers->positions.assign(&newParticles[0], newParticles.size());
		g_buffers->restPositions.assign(&newParticlesRest[0], newParticlesRest.size());
		g_buffers->velocities.assign(&newVelocities[0], newVelocities.size());
		g_buffers->phases.assign(&newPhases[0], newPhases.size());
		g_buffers->normals.assign(&newNormals[0], newNormals.size());

		// build active indices list
		g_buffers->activeIndices.resize(g_buffers->positions.size());
		for (int i = 0; i < g_buffers->positions.size(); ++i)
			g_buffers->activeIndices[i] = i;

		// update constraint buffers
		RebuildConstraints();

		if (write_output || (force_write_frame && (g_frame % 5) == 0)) {
			if (!force_write_frame)
			{
				WriteOutput();
				WriteOutputOBJ();
			}
			WriteOutputPH();
			write_output = false;
			if (force_write_frame && (g_frame % 5) == 0) write_output = true;
		}

		// restore mouse mass		
		if (g_mouseParticle != -1)
			g_buffers->positions[g_mouseParticle].w = 0.0f;
	}

	virtual void Draw(int pass)
	{
		if (!g_drawMesh)
			return;

		for (size_t i = 0; i < mCloths.size(); ++i)
		{
			DrawCloth(&g_buffers->positions[0], &g_buffers->normals[0], NULL, &g_buffers->triangles[mCloths[i].triangleOffset], mCloths[i].asset->numTriangles, g_buffers->positions.size(), (i + 2) % 6);//, g_params.radius*0.25f);			
		}
	}

	virtual void DoGui()
	{
		if (imguiSlider("stiffness", &stiffness, 0.001f, 2.0f, 0.001f))
		{
			//for (int i = 0; i < int(g_buffers->springStiffness.size()); ++i)
		    //		g_buffers->springStiffness[i] = stiffness;
			for (int c = 0; c < int(mCloths.size()); ++c)
			{
				Balloon& balloon = mCloths[c];
				for (int i = 0; i < balloon.asset->numSprings; ++i)
				{
					balloon.asset->springCoefficients[i] = stiffness;
					//g_buffers->springStiffness.push_back(balloon.asset->springCoefficients[i]);
					//g_buffers->springLengths.push_back(balloon.asset->springRestLengths[i]);
				}
			}
		}
		if (imguiSlider("vlength_scale", &vlength_scale, 0.01f, 2.0f, 0.01f))
		{
			// g_params.radius = (vlength * vlength_scale) / 2.0;
			//this affect in the same way p5 and p6
			for (int c = 0; c < int(mCloths.size()); ++c)
			{
				Balloon& balloon = mCloths[c];
				for (int i = 0; i < balloon.asset->numSprings; ++i)
				{
					if ((i % 2) == 0) balloon.asset->springRestLengths[i] = vlength; // 0.09f radius even 
					else balloon.asset->springRestLengths[i] = vlength * vlength_scale;  // 0.15f  distance to next hexagon should be 11 odd
				}
			}
		}
		if (imguiSlider("vlength", &vlength, 0.001f, 0.2f, 0.001f))
		{
			// g_params.radius = (vlength * vlength_scale) / 2.0;
			for (int c = 0; c < int(mCloths.size()); ++c)
			{
				Balloon& balloon = mCloths[c];
				for (int i = 0; i < balloon.asset->numSprings; ++i)
				{
					if ((i % 2) == 0) balloon.asset->springRestLengths[i] = vlength; // 0.09f radius even 
					else balloon.asset->springRestLengths[i] = vlength * vlength_scale;  // 0.15f  distance to next hexagon should be 11 odd
				}
			}
		}
		if (imguiCheck("write output", write_output))
		{
			write_output = !write_output;
		}
		if (imguiCheck("write output every frame", force_write_frame))
		{
			force_write_frame = !force_write_frame;
		}
		if (imguiCheck("write output with unique name", unique_name))
		{
			unique_name = !unique_name;
		}
	}

	inline Vec3 RotateAroundAxis(const Vec3& v, const Vec3& axis, float angleRad)
	{
		// Rodrigues' rotation formula
		return v * cosf(angleRad)
			+ Cross(axis, v) * sinf(angleRad)
			+ axis * Dot(axis, v) * (1.0f - cosf(angleRad));
	}

	Vec3 SnapHexTangent60(
		const Vec3& p,
		const Vec3& n,
		const Vec3& t0,                    // initial tangent (unit, ⟂ n)
		const std::vector<int>& neighbors // global indices of neighbors (should be 6)
	)
	{
		const float sixty = 2.0f * M_PI / 6.0f;

		float bestScore = -1e30f;
		Vec3 bestT = t0;

		for (int k = 0; k < 6; ++k)
		{
			float angle = k * sixty;
			Vec3 t = RotateAroundAxis(t0, n, angle);

			float score = 0.0f;

			for (int gj : neighbors)
			{
				Vec3 d = Vec3(g_buffers->positions[gj]) - p;

				// project neighbor direction into tangent plane
				d -= Dot(d, n) * n;
				float dl = Length(d);
				if (dl < 1e-8f) continue;
				d /= dl;

				// hex lattice prefers ±t alignment
				score += fabsf(Dot(t, d));
			}

			if (score > bestScore)
			{
				bestScore = score;
				bestT = t;
			}
		}

		return Normalize(bestT);
	}

	virtual void WriteOutput() {
		// write current data in files
		// we want one file with all the vertices
		// another with the dual faces
		// std:string outPath = "../../data/cellpack/capsides_sim_" + std::to_string(g_frame) + ".xyz";
		std::string fn = (unique_name) ? "0" : std::to_string(g_frame);
		std::string outPath =
			"G:\\Dev\\cellPACKgpu.git\\Data\\lattices\\hiv\\capsides_sim_" + fn + ".xyz";
		std::ofstream out(outPath);
		if (!out.is_open())
		{
			std::cerr << "Failed to open " << outPath << " for writing.\n";
			return;
		}
		// Count total particles across all balloons/cloths
		// out << std::fixed << std::setprecision(6);

		for (int c = 0; c < int(mCloths.size()); ++c)
		{
			Balloon& balloon = mCloths[c];
			const int n = balloon.asset->numParticles;
			const int base = balloon.particleOffset;

			// --- XYZ "frame" header ---
			out << n << "\n";
			out << "capsid=" << c << " offset=" << base << " n=" << n << "\n";

			// --- points ---
			for (int i = 0; i < n; ++i)
			{
				char elem = (i < 12) ? char('P') : 'H';
				const Vec4 p = g_buffers->positions[base + i];
				const Vec4 n = g_buffers->normals[base + i];
				out << elem << " " << p.x / main_scale << " " << p.y / main_scale << " " << p.z / main_scale << " " << n.x << " " << n.y << " " << n.z << "\n";
			}
		}
		out.close();
		// std::cout << "Wrote " << mCloths.size() << " frames to " << outPath << "\n";
	}

	virtual void WriteOutputPH()
	{
		std::string fn = (unique_name) ? "0" : std::to_string(g_frame);
		//std::string outPath = "../../data/cellpack/ph_capsides_sim_" + std::to_string(g_frame) + ".xyz";
		std::string outPath =
			"G:\\Dev\\cellPACKgpu.git\\Data\\lattices\\hiv\\ph_capsides_sim_" + fn + ".xyz";
		
		std::ofstream out(outPath);
		if (!out.is_open())
		{
			std::cerr << "Failed to open " << outPath << " for writing.\n";
			return;
		}

		// Count total particles across all balloons/cloths
		size_t totalPts = 0;
		for (int bi = 0; bi < int(mCloths.size()); ++bi)
			totalPts += size_t(mCloths[bi].asset->numParticles);

		// XYZ header
		out << totalPts << "\n";
		out << "P/H capsid points with normal+tangent (frame=" << g_frame << ")\n";

		// Helpers
		auto safeNormalize = [](const Vec3& v, float eps = 1e-8f) -> Vec3 {
			float l = Length(v);
			if (l < eps) return Vec3(0.0f, 1.0f, 0.0f);
			return v / l;
		};

		auto projectToPlane = [](const Vec3& v, const Vec3& n) -> Vec3 {
			return v - Dot(v, n) * n;
		};

		for (int bi = 0; bi < int(mCloths.size()); ++bi)
		{
			Balloon& balloon = mCloths[bi];
			const int numVerts = balloon.asset->numParticles;
			const int base = balloon.particleOffset;

			for (int vi = 0; vi < numVerts; ++vi)
			{
				char elem = (vi < 12) ? 'P' : 'H';

				const int gi = base + vi;                 // global particle index
				const Vec4 pv4 = g_buffers->positions[gi];
				const Vec3 p = Vec3(pv4);

				// --- 1) compute vertex normal (area-weighted) ---
				Vec3 nsum(0, 0, 0);

				// Collect neighbors (global indices) from incident triangles
				// degree is small (5 or 6), so simple vector is fine
				int neighbors[16];
				int nbCount = 0;

				auto addNeighbor = [&](int gj) {
					if (gj == gi) return;
					for (int k = 0; k < nbCount; ++k)
						if (neighbors[k] == gj) return;
					if (nbCount < 16) neighbors[nbCount++] = gj;
				};

				for (int triLocal : balloon.vertToTris[vi])
				{
					int triIntStart = balloon.triangleOffset + triLocal * 3;

					int ga = g_buffers->triangles[triIntStart + 0];
					int gb = g_buffers->triangles[triIntStart + 1];
					int gc = g_buffers->triangles[triIntStart + 2];

					const Vec3 p0 = Vec3(g_buffers->positions[ga]);
					const Vec3 p1 = Vec3(g_buffers->positions[gb]);
					const Vec3 p2 = Vec3(g_buffers->positions[gc]);

					nsum += Cross(p1 - p0, p2 - p0); // area-weighted
					// std::cout << "edge distance " << Length(p1 - p0) << "\n" << std::endl; //0.06 0.059
					// add the other two vertices of this triangle as neighbors
					if (ga == gi) { addNeighbor(gb); addNeighbor(gc); }
					else if (gb == gi) { addNeighbor(ga); addNeighbor(gc); }
					else if (gc == gi) { addNeighbor(ga); addNeighbor(gb); }
				}

				Vec3 n = safeNormalize(nsum);

				// --- 2) choose a STABLE tangent from neighbors ---
				// Build a reference direction in tangent plane to break +/- ambiguity
				Vec3 ref = Cross(Vec3(0, 1, 0), n);
				if (Length(ref) < 1e-8f) ref = Cross(Vec3(1, 0, 0), n);
				ref = safeNormalize(ref);

				Vec3 bestT(0, 0, 0);
				float bestScore = -1e30f;

				for (int k = 0; k < nbCount; ++k)
				{
					int gj = neighbors[k];
					Vec3 d = Vec3(g_buffers->positions[gj]) - p;   // edge direction
					Vec3 tproj = projectToPlane(d, n);
					float tl = Length(tproj);
					if (tl < 1e-8f) continue;

					tproj /= tl;

					// choose neighbor direction most aligned with ref (stable selection)
					float score = Dot(tproj, ref);
					if (score > bestScore)
					{
						bestScore = score;
						bestT = tproj;
					}
				}

				// Fallback if something degenerate happens
				Vec3 t = bestT;
				if (Length(t) < 1e-8f)
					t = ref;

				// --- 3) write out ---
				out << elem << " "
					<< p.x / main_scale << " " << p.y / main_scale << " " << p.z / main_scale << " "
					<< n.x << " " << n.y << " " << n.z << " "
					<< t.x << " " << t.y << " " << t.z << "\n";
			}
		}

		out.close();
		std::cout << "Wrote " << totalPts << " points to " << outPath << "\n";
	}


	static void SortNeighborsAround(
		const Vec3& center,
		const Vec3& n,                  // desired plane normal (approx ok)
		int* neighbors, int nbCount
	) {
		// Build a local 2D basis (u,v) in the plane
		Vec3 u = Cross(Vec3(0, 1, 0), n);
		if (Length(u) < 1e-8f) u = Cross(Vec3(1, 0, 0), n);
		u = u / Length(u);
		Vec3 v = Cross(n, u); // already normalized

		// compute angle for each neighbor
		struct Item { int id; float ang; };
		Item items[16];

		for (int k = 0; k < nbCount; ++k)
		{
			Vec3 pk = Vec3(g_buffers->positions[neighbors[k]]) - center;
			// project to plane
			pk -= Dot(pk, n) * n;

			float x = Dot(pk, u);
			float y = Dot(pk, v);
			items[k] = { neighbors[k], atan2f(y, x) };
		}

		// sort by angle
		std::sort(items, items + nbCount, [](const Item& a, const Item& b) { return a.ang < b.ang; });

		for (int k = 0; k < nbCount; ++k)
			neighbors[k] = items[k].id;
	}

	virtual void WriteOutputPH1()
	{
		std::string outPath =
			"G:\\Dev\\cellPACKgpu.git\\Data\\lattices\\hiv\\ph_capsides_sim_" + std::to_string(g_frame) + ".xyz";

		std::ofstream out(outPath);
		if (!out.is_open())
		{
			std::cerr << "Failed to open " << outPath << " for writing.\n";
			return;
		}

		// Count total particles across all balloons/cloths
		size_t totalPts = 0;
		for (int bi = 0; bi < int(mCloths.size()); ++bi)
			totalPts += size_t(mCloths[bi].asset->numParticles);

		// XYZ header
		out << totalPts << "\n";
		out << "P/H capsid points with neighbor-plane normal + tangent (frame=" << g_frame << ")\n";

		auto safeNormalize = [](const Vec3& v, float eps = 1e-8f) -> Vec3 {
			float l = Length(v);
			if (l < eps) return Vec3(0.0f, 0.0f, 0.0f);
			return v / l;
		};

		auto projectToPlane = [](const Vec3& v, const Vec3& n) -> Vec3 {
			return v - Dot(v, n) * n;
		};

		for (int bi = 0; bi < int(mCloths.size()); ++bi)
		{
			Balloon& balloon = mCloths[bi];

			const int numVerts = balloon.asset->numParticles;
			const int base = balloon.particleOffset;

			// Cone axis is assumed precomputed and normalized
			const Vec3 coneAxis = safeNormalize(balloon.coneAxis);
			// A point near the cone centerline (mean of vertices)
			Vec3 axisPoint(0, 0, 0);
			for (int vi = 0; vi < numVerts; ++vi)
				axisPoint += Vec3(g_buffers->positions[base + vi]);
			axisPoint /= float(std::max(1, numVerts));

			for (int vi = 0; vi < numVerts; ++vi)
			{
				const char elem = (vi < 12) ? 'P' : 'H';

				const int gi = base + vi;     // global particle index
				const Vec3 pc = Vec3(g_buffers->positions[gi]);

				// ---- Collect neighbors from incident triangles ----
				int neighbors[16];
				int nbCount = 0;

				auto addNeighbor = [&](int gj) {
					if (gj == gi) return;
					for (int k = 0; k < nbCount; ++k)
						if (neighbors[k] == gj) return;
					if (nbCount < 16) neighbors[nbCount++] = gj;
				};

				for (int triLocal : balloon.vertToTris[vi])
				{
					const int triIntStart = balloon.triangleOffset + triLocal * 3;
					const int ga = g_buffers->triangles[triIntStart + 0];
					const int gb = g_buffers->triangles[triIntStart + 1];
					const int gc = g_buffers->triangles[triIntStart + 2];

					if (ga == gi) { addNeighbor(gb); addNeighbor(gc); }
					else if (gb == gi) { addNeighbor(ga); addNeighbor(gc); }
					else if (gc == gi) { addNeighbor(ga); addNeighbor(gb); }
				}

				if (nbCount < 3)
				{
					// Degenerate; write something sane and continue
					std::cout << "Degenerate; write something sane and continue\n" << endl;
					out << elem << " "
						<< pc.x / main_scale << " " << pc.y / main_scale << " " << pc.z / main_scale << " "
						<< 0 << " " << 1 << " " << 0 << " "
						<< 1 << " " << 0 << " " << 0 << "\n";
					continue;
				}
				std::cout << "nbCount " << nbCount << endl;
				
				// ---- Compute neighbor ring center ----
				Vec3 center(0, 0, 0);
				for (int k = 0; k < nbCount; ++k)
					center += Vec3(g_buffers->positions[neighbors[k]]);
				center /= float(nbCount);

				// ---- Compute plane normal from ordered neighbor ring ----

				// First, get a rough normal to define the sorting plane.
				// Using center->pc (or cone radial) is a good rough normal.
				Vec3 roughN = pc - axisPoint;
				roughN = roughN - Dot(roughN, coneAxis) * coneAxis;
				if (Length(roughN) < 1e-8f) roughN = coneAxis;
				roughN = safeNormalize(roughN);

				// Sort neighbors in a consistent loop around roughN
				SortNeighborsAround(center, roughN, neighbors, nbCount);

				// Now polygon normal using consecutive edges
				Vec3 nsum(0, 0, 0);
				for (int k = 0; k < nbCount; ++k)
				{
					int a = neighbors[k];
					int b = neighbors[(k + 1) % nbCount];

					Vec3 va = Vec3(g_buffers->positions[a]) - center;
					Vec3 vb = Vec3(g_buffers->positions[b]) - center;

					// project to plane of roughN to reduce noise
					va -= Dot(va, roughN) * roughN;
					vb -= Dot(vb, roughN) * roughN;

					nsum += Cross(va, vb);
				}

				Vec3 n = safeNormalize(nsum);
				if (Length(n) < 1e-8f)
					n = roughN;

				// ---- Tangent: center -> "first" neighbor (simple, not sorted) ----
				Vec3 t(0, 0, 0);
				{
					const Vec3 p0 = Vec3(g_buffers->positions[neighbors[0]]);
					t = p0 - center;

					// Always keep tangent in the plane (cheap + stable)
					// t = projectToPlane(t, n);
					t = safeNormalize(t);

					// if (Length(t) < 1e-8f)
					// {
						// fallback: use cone axis projected into plane
					// 	t = projectToPlane(coneAxis, n);
					// 	t = safeNormalize(t);
					// }
				}

				// ---- Output: position = ring center (your current choice) ----
				out << elem << " "
					<< center.x / main_scale << " " << center.y / main_scale << " " << center.z / main_scale << " "
					<< n.x << " " << n.y << " " << n.z << " "
					<< t.x << " " << t.y << " " << t.z << "\n";
			}
		}

		out.close();
		std::cout << "Wrote " << totalPts << " points to " << outPath << "\n";
	}


	virtual void WriteOutputOBJ() {
		// write current data in files
		// we want one file with all the vertices
		// another with the dual faces
		std::string fn = (unique_name) ? "0" : std::to_string(g_frame);
		std::string outPath =
			"G:\\Dev\\cellPACKgpu.git\\Assets\\Geoms\\HIV\\capsides_sim_" + fn + ".obj";
		std::ofstream out(outPath);
		if (!out.is_open())
		{
			std::cerr << "Failed to open " << outPath << " for writing.\n";
			return;
		}
		out << "# Flex balloons/cloths export (one object per balloon)\n";

		for (int c = 0; c < int(mCloths.size()); ++c)
		{
			Balloon& balloon = mCloths[c];
			const int n = balloon.asset->numParticles;
			const int base = balloon.particleOffset;
			out << "\n";
			out << "o balloon_" << c << "\n";
			
			// --- vertices ---
			for (int i = 0; i < n; ++i)
			{
				const Vec4 p = g_buffers->positions[base + i];
				out << "v " << p.x / main_scale << " " << p.y / main_scale << " " << p.z / main_scale << "\n";// << n.x << " " << n.y << " " << n.z << "\n";
			}
			// --- face indices ---
			const int f = balloon.asset->numTriangles;
			const int fbase = balloon.triangleOffset;

			for (int i = 0; i < f; i++)
			{
				int vc = g_buffers->triangles[fbase + i * 3 + 0];
				int vb = g_buffers->triangles[fbase + i * 3 + 1];
				int va = g_buffers->triangles[fbase + i * 3 + 2];
				//Vec4 fnormal = g_buffers->triangleNormals[fbase + f];
				out << "f " << (va + 1) << " " << (vb + 1) << " " << (vc + 1) << "\n";
			}
			out << "# balloon_" << c
				<< " verts=" << n
				<< " tris=" << f << "\n";
		}
		out.close();
		// std::cout << "Wrote balloons OBJ to " << outPath << "\n";
	}
	
	struct Balloon
	{
		NvFlexExtAsset* asset;

		int particleOffset;
		int triangleOffset;

		bool tearable;

		float splitThreshold;

		Vec3 coneAxis;

		// vertex -> list of local triangle indices
		std::vector<std::vector<int>> vertToTris;
	};

	int mNumFluidParticles;

	std::vector<Balloon> mCloths;
};
