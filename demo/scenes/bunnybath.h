extern std::map<NvFlexTriangleMeshId, GpuMesh*> g_meshes;

class BunnyBath : public Scene
{
public:

	BunnyBath(const char* name, bool dam) : Scene(name), mDam(dam) {}

	void Initialize()
	{
		float radius = 0.14f*2.0f;

		// deforming bunny
		float s = radius*0.5f;
		float m = 1.0f;
		int group = 1;
		//(const char* filename, Vec3 lower, Vec3 scale, float rotation, float spacing, Vec3 velocity, float invMass, bool rigid, float rigidStiffness, int phase, bool skin, float jitter=0.005f, Vec3 skinOffset=0.0f, float skinExpand=0.0f, Vec4 color=Vec4(0.0f), float springStiffness=0.0f)
		// CreateParticleShape(GetFilePathByPlatform("../../data/bunny.ply").c_str(), Vec3(4.0f, 0.0f, 0.0f), 0.5f, 0.0f, s, Vec3(0.0f, 0.0f, 0.0f), m, false, 1.0f, NvFlexMakePhase(group++, 0), true, 0.0f);
		// CreateParticleShape(GetFilePathByPlatform("../../data/box.ply").c_str(), Vec3(4.0f, 0.0f, 1.0f), 0.45f, 0.0f, s, Vec3(0.0f, 0.0f, 0.0f), m, true, 1.0f, NvFlexMakePhase(group++, 0), true, 0.0f);
		// CreateParticleShape(GetFilePathByPlatform("../../data/bunny.ply").c_str(), Vec3(3.0f, 0.0f, 0.0f), 0.5f, 0.0f, s, Vec3(0.0f, 0.0f, 0.0f), m, true, 1.0f, NvFlexMakePhase(group++, 0), true, 0.0f);
		// CreateParticleShape(GetFilePathByPlatform("../../data/sphere.ply").c_str(), Vec3(3.0f, 0.0f, 1.0f), 0.45f, 0.0f, s, Vec3(0.0f, 0.0f, 0.0f), m, true, 1.0f, NvFlexMakePhase(group++, 0), true, 0.0f);
		// CreateParticleShape(GetFilePathByPlatform("../../data/bunny.ply").c_str(), Vec3(2.0f, 0.0f, 1.0f), 0.5f, 0.0f, s, Vec3(0.0f, 0.0f, 0.0f), m, true, 1.0f, NvFlexMakePhase(group++, 0), true, 0.0f);
		// CreateParticleShape(GetFilePathByPlatform("../../data/box.ply").c_str(), Vec3(2.0f, 0.0f, 0.0f), 0.45f, 0.0f, s, Vec3(0.0f, 0.0f, 0.0f), m, true, 1.0f, NvFlexMakePhase(group++, 0), true, 0.0f);
		// CreateParticleShape(GetFilePathByPlatform("../../data/4STD.ply").c_str(), Vec3(1.0f, 0.0f, 0.0f), 0.5f, 0.0f, s, Vec3(0.0f, 0.0f, 0.0f), m, true, 1.0f, NvFlexMakePhase(group++, 0), true, 0.0f);
		const float meshScale = 0.1f;
		// Mesh* mesh1 = ImportMesh(GetFilePathByPlatform("../../data/SITES.ply").c_str());
		// mesh1->Transform(ScaleMatrix(meshScale));
		
		// AddTriangleMesh(CreateTriangleMesh(mesh1), Vec3(), Quat(), 1.0f);

		Mesh* mesh = ImportMesh(GetFilePathByPlatform("../../data/4STD.ply").c_str());
		// mesh->Normalize();
		
		mesh->Transform(ScaleMatrix(meshScale));
		mRenderMesh = mesh;
		mMesh = CreateTriangleMesh(mesh);
		g_meshes.erase(mMesh);
		LoadStickSitesFromPdb(GetFilePathByPlatform("../../data/sites.pdb").c_str(), meshScale);
		InitializeSiteProxyParticles();
		
		// AddTriangleMesh(mMesh, Vec3(), Quat(), 1.0f);

		g_numSolidParticles = 0; // g_buffers->positions.size();		

		float restDistance = radius*0.55f;

		if (mDam)
		{
			CreateParticleGrid(Vec3(0.0f, 0.0f, 0.0f), 32, 32, 32, restDistance, Vec3(0.0f), 1.0f, false, 0.0f, NvFlexMakePhase(0, eNvFlexPhaseSelfCollide | eNvFlexPhaseFluid), 0.005f);
			g_lightDistance *= 2.5f;
		}

		g_sceneLower = Vec3(0.0f);

		g_numSubsteps = 2;

		g_params.radius = radius;
		g_params.dynamicFriction = 0.01f; 
		g_params.viscosity = 2.0f;
		g_params.numIterations = 4;
		g_params.vorticityConfinement = 40.0f;
		g_params.fluidRestDistance = restDistance;
		g_params.solidPressure = 0.f;
		g_params.relaxationFactor = 0.0f;
		g_params.cohesion = 0.02f;
		g_params.collisionDistance = 0.01f;		

		g_maxDiffuseParticles = 64*1024;
		g_diffuseScale = 0.5f;

		g_fluidColor = Vec4(0.113f, 0.425f, 0.55f, 1.f);

		Emitter e1;
		e1.mDir = Vec3(1.0f, 0.0f, 0.0f);
		e1.mRight = Vec3(0.0f, 0.0f, -1.0f);
		e1.mPos = Vec3(radius, 1.f, 0.65f);
		e1.mSpeed = (restDistance/g_dt)*2.0f; // 2 particle layers per-frame
		e1.mEnabled = true;

		g_emitters.push_back(e1);

		g_numExtraParticles = 48*1024;

		g_lightDistance = 1.8f;

		g_clearColor = Vec3(1.0f);
		g_params.numPlanes = 5;

		g_waveFloorTilt = 0.0f;
		g_waveFrequency = 1.5f;
		g_waveAmplitude = 2.0f;
		
		g_warmup = true;

		// draw options
		g_drawPlane = false;
		g_drawPoints = false;
		g_drawEllipsoids = true;
		g_drawDiffuse = true;
	}

	void Update()
	{
		ClearShapes();

		if (!pauseMotion)
			mTime += g_dt;

		Vec3 pos, prevPos;
		Quat rot, prevRot;
		GetAnimatedPose(pos, prevPos, rot, prevRot);

		AddTriangleMesh(mMesh, pos, rot, 1.0f);

		UpdateSiteProxyParticles(pos, prevPos, rot);
		
		g_buffers->shapePrevPositions[0] = Vec4(prevPos, 0.0f);
		g_buffers->shapePrevRotations[0] = prevRot;

		UpdateShapes();
	}

	virtual void Draw(int pass)
	{
		if (!g_drawMesh || !mRenderMesh || pass != 0)
			return;

		Vec3 pos, prevPos;
		Quat rot, prevRot;
		GetAnimatedPose(pos, prevPos, rot, prevRot);

		Mesh m = *mRenderMesh;
		Matrix44 rotM = RotationMatrix(rot);
		Matrix44 xform = TranslationMatrix(Point3(pos))*rotM;

		for (int i = 0; i < int(m.m_positions.size()); ++i)
			m.m_positions[i] = Point3(xform*Vec4(Vec3(m.m_positions[i]), 1.0f));

		for (int i = 0; i < int(m.m_normals.size()); ++i)
			m.m_normals[i] = Normalize(Vec3(rotM*Vec4(m.m_normals[i], 0.0f)));

		DrawMesh(&m, Vec3(1.0f));
	}

	virtual void DoGui()
	{
		imguiSlider("x", &x, 0.0f, 5.0, 0.01f);
		imguiSlider("y", &y, 0.0f, 5.0, 0.01f);
		imguiSlider("z", &z, 0.0f, 5.0f, 0.01f);
		if (imguiCheck("Pause Motion", pauseMotion))
			pauseMotion = !pauseMotion;
		imguiSlider("Move Length", &moveLength, 0.0f, 5.0f, 0.01f);
		if (imguiCheck("Enable Sticking", enableSticking))
			enableSticking = !enableSticking;

		// imguiSlider("x", &ox, -1.0f, 1.0, 0.01f);
		// imguiSlider("y", &oy, -1.0f, 1.0, 0.01f);
		// imguiSlider("z", &oz, -1.0f, 1.0f, 0.01f);
		if (imguiCheck("Enable Site Proxy", enableSiteProxy))
			enableSiteProxy = !enableSiteProxy;
	}

	void GetAnimatedPose(Vec3& pos, Vec3& prevPos, Quat& rot, Quat& prevRot) const
	{
		float startTime = 1.0f;
		float time = Max(0.0f, mTime-startTime);
		float lastTime = Max(0.0f, time-g_dt);
		const float translationSpeed = 1.0f;

		pos = Vec3(x, y + moveLength*translationSpeed*(1.0f - cosf(time)), z);
		prevPos = Vec3(x, y + moveLength*translationSpeed*(1.0f - cosf(lastTime)), z);
		rot = QuatFromAxisAngle(Vec3(0.0f, 1.0f, 0.0f), 0.0f);
		prevRot = QuatFromAxisAngle(Vec3(0.0f, 1.0f, 0.0f), 0.0f);
	}

	void LoadStickSitesFromPdb(const char* path, float scale)
	{
		mStickSitesLocal.clear();
		mStickParticleForSite.clear();

		ifstream f(path);
		if (!f)
			return;

		string line;
		vector<Vec3> rawSites;

		while (getline(f, line))
		{
			if (line.size() < 54)
				continue;

			if (!(line.compare(0, 6, "HETATM") == 0 || line.compare(0, 4, "ATOM") == 0))
				continue;

			float px = 0.0f, py = 0.0f, pz = 0.0f;

			// Standard PDB columns: x[30:38], y[38:46], z[46:54]
			try
			{
				px = std::stof(line.substr(30, 8));
				py = std::stof(line.substr(38, 8));
				pz = std::stof(line.substr(46, 8));
			}
			catch (...)
			{
				continue;
			}

			rawSites.push_back(Vec3(px, py, pz));
		}

		if (!rawSites.empty())
		{
			Vec3 lower = rawSites[0];
			Vec3 upper = rawSites[0];

			for (int i = 1; i < int(rawSites.size()); ++i)
			{
				lower = Min(lower, rawSites[i]);
				upper = Max(upper, rawSites[i]);
			}

			Vec3 center = 0.5f * (lower + upper);
			mStickSitesLocal.reserve(rawSites.size());

			for (int i = 0; i < int(rawSites.size()); ++i)
				mStickSitesLocal.push_back((rawSites[i] - center) * scale);
		}

		mStickParticleForSite.assign(mStickSitesLocal.size(), -1);
	}

	void UpdateStuckParticles(const Vec3& pos, const Quat& rot)
	{
		if (!enableSticking || mStickSitesLocal.empty())
			return;

		if (pos.y < stickStartY)
			return;

		if (g_numSolidParticles >= int(g_buffers->positions.size()))
			return;

		const int firstFluid = int(g_numSolidParticles);
		const float captureRadiusSq = stickCaptureRadius * stickCaptureRadius;
		vector<Vec3> worldSites(mStickSitesLocal.size());
		vector<char> siteUsed(mStickSitesLocal.size(), 0);

		for (int s = 0; s < int(mStickSitesLocal.size()); ++s)
		{
			worldSites[s] = pos + Rotate(rot, mStickSitesLocal[s]);

			int particleIndex = mStickParticleForSite[s];
			if (particleIndex >= firstFluid && particleIndex < int(g_buffers->positions.size()))
			{
				siteUsed[s] = 1;
				g_buffers->positions[particleIndex] = Vec4(worldSites[s], 0.0f);
				if (particleIndex < int(g_buffers->velocities.size()))
					g_buffers->velocities[particleIndex] = Vec3(0.0f);
			}
		}

		for (int i = firstFluid; i < int(g_buffers->positions.size()); ++i)
		{
			bool alreadyAssigned = false;
			for (int s = 0; s < int(mStickParticleForSite.size()); ++s)
			{
				if (mStickParticleForSite[s] == i)
				{
					alreadyAssigned = true;
					break;
				}
			}
			if (alreadyAssigned)
				continue;

			Vec3 p = Vec3(g_buffers->positions[i]);
			float bestDistSq = captureRadiusSq;
			int bestSite = -1;

			for (int s = 0; s < int(worldSites.size()); ++s)
			{
				if (siteUsed[s])
					continue;

				float d2 = LengthSq(p - worldSites[s]);
				if (d2 < bestDistSq)
				{
					bestDistSq = d2;
					bestSite = s;
				}
			}

			if (bestSite != -1)
			{
				mStickParticleForSite[bestSite] = i;
				siteUsed[bestSite] = 1;
				g_buffers->positions[i] = Vec4(worldSites[bestSite], 0.0f);
				if (i < int(g_buffers->velocities.size()))
					g_buffers->velocities[i] = Vec3(0.0f);
			}
		}
	}

	void InitializeSiteProxyParticles()
	{
		mSiteProxyParticleIndices.clear();
		if (!enableSiteProxy || mStickSitesLocal.empty())
			return;

		Vec3 pos, prevPos;
		Quat rot, prevRot;
		GetAnimatedPose(pos, prevPos, rot, prevRot);

		int phase = NvFlexMakePhase(1, 0);

		for (int s = 0; s < int(mStickSitesLocal.size()); ++s)
		{
			Vec3 p = pos + Rotate(rot, mStickSitesLocal[s]);
			mSiteProxyParticleIndices.push_back(int(g_buffers->positions.size()));
			g_buffers->positions.push_back(Vec4(p, 0.0f)); // kinematic particle
			g_buffers->velocities.push_back(Vec3(0.0f));
			g_buffers->phases.push_back(phase);
		}
		// Build once (outside the loop / constructor init if possible)
		// static std::vector<Vec3> sphereOffsets;
		// if (sphereOffsets.empty())
		// {
		// 	const float proxyRadius  = g_params.radius * 2.0f;   // target sphere radius = 0.28
		// 	const float spacing      = g_params.radius;          // distance between proxy particles (tune this)
		// 	const float r2           = proxyRadius * proxyRadius;

		// 	const int n = int(std::ceil(proxyRadius / spacing));

		// 	for (int ix = -n; ix <= n; ++ix)
		// 	{
		// 		for (int iy = -n; iy <= n; ++iy)
		// 		{
		// 			for (int iz = -n; iz <= n; ++iz)
		// 			{
		// 				Vec3 d(float(ix) * spacing, float(iy) * spacing, float(iz) * spacing);
		// 				sphereOffsets.push_back(d);
		// 				// const float d2 = d.x * d.x + d.y * d.y + d.z * d.z;
		// 				// if (d2 <= r2)
		// 				// {
		// 				// 	sphereOffsets.push_back(d);
		// 				// }
		// 			}
		// 		}
		// 	}
		// }	
		// // For each stick site, stamp the sphere of particles
		// for (int s = 0; s < int(mStickSitesLocal.size()); ++s)
		// {
		// 	Vec3 center = pos + Rotate(rot, mStickSitesLocal[s]);

		// 	// (Optional) if you need to know where this site's proxies start:
		// 	// int startIndex = int(g_buffers->positions.size());

		// 	for (const Vec3& off : sphereOffsets)
		// 	{
		// 		Vec3 p = center + off;  // sphere is isotropic, no need to rotate 'off'

		// 		mSiteProxyParticleIndices.push_back(int(g_buffers->positions.size()));
		// 		g_buffers->positions.push_back(Vec4(p, 0.0f)); // kinematic particle
		// 		g_buffers->velocities.push_back(Vec3(0.0f));
		// 		g_buffers->phases.push_back(phase);
		// 	}

		// 	// (Optional) if downstream code assumes 1 proxy per site, store ranges instead:
		// 	// mSiteProxyStarts.push_back(startIndex);
		// 	// mSiteProxyCounts.push_back(int(g_buffers->positions.size()) - startIndex);
		// }		
	}

	void UpdateSiteProxyParticles(const Vec3& posin, const Vec3& prevPos, const Quat& rot)
	{
		if (!enableSiteProxy)
			return;
		Vec3 pos = posin;
		// Only follow the mesh while it moves upward.
		if (pos.y <= prevPos.y)
			pos = Vec3(pos.x, -2.0f, pos.z);

		Vec3 delta = pos - prevPos;
		Vec3 v = delta / g_dt;

		int n = mSiteProxyParticleIndices.size();//Min(int(mSiteProxyParticleIndices.size()), int(mStickSitesLocal.size()));
		for (int s = 0; s < n; ++s)
		{
			if (s < 0 || s >= int(g_buffers->positions.size()))
				continue;

			Vec3 p = pos + Rotate(rot, mStickSitesLocal[s] + Vec3(ox, oy, oz));
			g_buffers->positions[s] = Vec4(p, 0.0f);
			if (s < int(g_buffers->velocities.size()))
				g_buffers->velocities[s] = 0.0f;
		}
	}

	float mTime = 0.0f;
	int mType;
	float x = 2.5f;
	float y = 1.13f;
	float z = 2.5f;
	float ox = -0.75f;
	float oy = 0.31f;
	float oz = -0.17f;	
	bool pauseMotion = false;
	float moveLength = 1.66f;
	bool enableSticking = true;
	float stickCaptureRadius = 0.05f;
	float stickStartY = 1.6f;
	bool enableSiteProxy = true;
	bool mDam;
	Mesh* mRenderMesh = NULL;
	vector<Vec3> mStickSitesLocal;
	vector<int> mStickParticleForSite;
	vector<int> mSiteProxyParticleIndices;
	NvFlexTriangleMeshId mMesh;
};
