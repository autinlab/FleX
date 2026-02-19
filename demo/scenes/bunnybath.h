
class BunnyBath : public Scene
{
public:

	BunnyBath(const char* name, bool dam) : Scene(name), mDam(dam) {}

	void Initialize()
	{
		float radius = 0.14f;

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
		
		Mesh* mesh = ImportMesh(GetFilePathByPlatform("../../data/4STD.ply").c_str());
		// mesh->Normalize();
		mesh->Transform(ScaleMatrix(0.1f));
		mMesh = CreateTriangleMesh(mesh);
		
		// AddTriangleMesh(mMesh, Vec3(), Quat(), 1.0f);

		g_numSolidParticles = g_buffers->positions.size();		

		float restDistance = radius*0.55f;

		if (mDam)
		{
			CreateParticleGrid(Vec3(0.0f, 0.0f, 0.0f), 64, 64, 64, restDistance, Vec3(0.0f), 1.0f, false, 0.0f, NvFlexMakePhase(0, eNvFlexPhaseSelfCollide | eNvFlexPhaseFluid), 0.005f);
			g_lightDistance *= 0.5f;
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

		g_params.numPlanes = 5;

		g_waveFloorTilt = 0.0f;
		g_waveFrequency = 1.5f;
		g_waveAmplitude = 2.0f;
		
		g_warmup = true;

		// draw options		
		g_drawPoints = false;
		g_drawEllipsoids = true;
		g_drawDiffuse = true;
	}

	void Update()
	{
		ClearShapes();

		mTime += g_dt;

		// let cloth settle on object
		float startTime = 1.0f;

		float time = Max(0.0f, mTime-startTime);
		float lastTime = Max(0.0f, time-g_dt);

		const float rotationSpeed = 1.0f;
		const float translationSpeed = 1.0f;

		//Vec3 pos = Vec3(translationSpeed*(1.0f-cosf(time)), 0.5f, 0.0f);
		//Vec3 prevPos = Vec3(translationSpeed*(1.0f-cosf(lastTime)), 0.5f, 0.0f);

		Vec3 pos = Vec3(x, y + translationSpeed*(1.0f - cosf(time)), z);
		Vec3 prevPos = Vec3(x, y + translationSpeed*(1.0f - cosf(lastTime)), z);

		Quat rot = QuatFromAxisAngle(Vec3(0.0f, 1.0f, 0.0f), 0.0f);// 1.0f - cosf(rotationSpeed*time));
		Quat prevRot = QuatFromAxisAngle(Vec3(0.0f, 1.0f, 0.0f), 0.0f);//1.0f-cosf(rotationSpeed*lastTime));

		AddTriangleMesh(mMesh, pos, rot, 1.0f);
		
		g_buffers->shapePrevPositions[0] = Vec4(prevPos, 0.0f);
		g_buffers->shapePrevRotations[0] = prevRot;

		UpdateShapes();
	}

	virtual void DoGui()
	{
		imguiSlider("x", &x, 0.0f, 5.0, 0.01f);
		imguiSlider("y", &y, 0.0f, 5.0, 0.01f);
		imguiSlider("z", &z, 0.0f, 5.0f, 0.01f);
	}

	float mTime;
	int mType;
	float x = 2.5f;
	float y = 1.4f;
	float z = 2.5f;
	bool mDam;
	NvFlexTriangleMeshId mMesh;
};
