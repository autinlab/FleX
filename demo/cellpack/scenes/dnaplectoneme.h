
class DNAplectoneme : public Scene
{
public:

	DNAplectoneme(const char* name) : Scene(name) {}
	float main_scale = 1.0f / 100.0f;
	float scale_angle = 0.0f;
	float angle = 0.0f;
	float rotation = 0.0f;
	float rotationSpeed = 0.0f;
	float strength = 1.0f;
	float jitter_strength = 1.0f;

	bool savetxt = false;
	bool use_position = false;
	bool use_jitter = false;
	bool write_output = false;

	NvFlexExtForceField forcefield;
	NvFlexExtForceFieldCallback* callback;

	void Initialize()
	{
		g_buffers->positions.resize(0);
		g_buffers->velocities.resize(0);
		g_buffers->phases.resize(0);

		g_buffers->rigidTranslations.resize(0);
		g_buffers->rigidRotations.resize(0);
		g_buffers->rigidCoefficients.resize(0);
		g_buffers->rigidIndices.resize(0);
		g_buffers->rigidLocalPositions.resize(0);
		g_buffers->rigidOffsets.resize(0);

		rotation = 0.0f;
		rotationSpeed = 0.0f;
		angle = 0;
		float radius = 0.1f;
		//1beads is 10bp, 1 bp is 3.4A->34
		float beads_radius = 2.0f*11.85f*main_scale;//5.0f*main_scale;//0.08f;//5.0f*main_scale;//smallest sphere size
		g_params.radius = beads_radius * 2.0f;

		//g_meshcellPACK * cp;
		g_cp = new cellPACK();
		g_cp->use_rb = true;
		g_cp->main_radius = beads_radius;//((10.0f/100.0f)*(10.0f/100.0f));
		g_cp->main_scale = main_scale;
		g_cp->mGroupCounter = 0;
		g_cp->mInstances.resize(0);

		//experiment is
		g_cp->use_instances_mesh = false;//or true
		g_cp->loadRecipe("../../data/cellpack/recipes/DNAplectoneme.1.0.json");
		g_cp->placeFibers(true);//check partners
		g_cp->printBatchId();
		/*g_params.dynamicFriction = 0.6f;
		g_params.staticFriction = 0.35f;
		g_params.particleFriction = 0.25f;
		g_params.dissipation = 0.0f;
		g_params.numIterations = 2;
		g_params.viscosity = 0.0f;
		g_params.drag = 0.0f;
		g_params.lift = 0.0f;
		*/
		g_params.numPlanes = 0;
		g_params.gravity[1] = 0;// -9.f;
		g_params.gravity[2] = 0;// -10;// -9.f;
		g_params.damping = 0.0f;// 3.0f;
		// update collision planes to match flexs


		g_params.collisionDistance = radius*0.5f;
		g_params.particleCollisionMargin = radius*0.25f;

		g_params.numIterations = 6;
		g_numSubsteps = 6;

		g_numExtraParticles = 1024 * 1024;

		g_drawPoints = true;
		g_wireframe = false;
		g_drawMesh = false;
		g_drawRopes = false;
		g_pause = true;

		//g_numSubsteps = 2;

		g_lightDistance *= 3.0f;
		callback = NULL;
		//add soolvent ?
	}

	void rotateAroundI(int i, float arotation, int off=1) {
		//float strength = 1.0f;
		float mult = 1.0f / g_dt;// Length(g_buffers->velocities[i]);
		Vec3 p1 = Vec3(g_buffers->positions[i + off]);
		Vec3 p2 = Vec3(g_buffers->positions[i - off]);
		
		Vec3 v1 = Vec3(g_buffers->positions[i + off]) - Vec3(g_buffers->positions[i]);
		Vec3 v2 = Vec3(g_buffers->positions[i - off]) - Vec3(g_buffers->positions[i]);
		Vec3 axis = Normalize(v1 + v2);
		Quat q = QuatFromAxisAngle(axis, arotation);//should start at 0
		Vec3 center = Vec3(g_buffers->positions[i]);
		Vec3 newp0 = Rotate(q, v1);
		Vec3 newp1 = Rotate(q, v2);
		//use position
		if (use_position){
			g_buffers->positions[i + off] = Vec4(center + newp0, g_buffers->positions[i + off].w);
			g_buffers->positions[i - off] = Vec4(center + newp1, g_buffers->positions[i - off].w);
		}
		else {
			//use velocity
			g_buffers->velocities[i + off] = Normalize(Vec3(center + newp0) - p1) * strength * mult;
			g_buffers->velocities[i - off] = Normalize(Vec3(center + newp1) - p2) * strength * mult;
		}
	}

	void Update()
	{
		rotation = rotationSpeed;
		int nParticles = g_buffers->activeIndices.size();
		//rotate 2 particles
		if (angle!=0)
		{
			rotateAroundI(0+(int)angle, rotation, (int)angle);
			//g_buffers->positions[0].w = 0.0f;
			//g_buffers->velocities[0] = Vec3(0,0,0);
			//g_buffers->positions[1].w = 0.0f;
			//g_buffers->velocities[1] = Vec3(0, 0, 0);
			rotateAroundI((int)((float)nParticles*(1.0 / 3.0)), rotation, (int)angle);
			rotateAroundI((int)((float)nParticles*(2.0 / 3.0)), rotation);
			/*
			Vec3 v1 = Vec3(g_buffers->positions[2].x, g_buffers->positions[2].y, g_buffers->positions[2].z) - Vec3(g_buffers->positions[1].x, g_buffers->positions[1].y, g_buffers->positions[1].z);
			Vec3 v2 = Vec3(g_buffers->positions[0].x, g_buffers->positions[0].y, g_buffers->positions[0].z) - Vec3(g_buffers->positions[1].x, g_buffers->positions[1].y, g_buffers->positions[1].z);
			Vec3 axis = Normalize(v1 + v2);
			Quat q = QuatFromAxisAngle(axis, rotation);//should start at 0
			Vec3 center = Vec3(g_buffers->positions[1].x, g_buffers->positions[1].y, g_buffers->positions[1].z);
			Vec3 newp0 = Rotate(q, v1);
			Vec3 newp1 = Rotate(q, v2);
			g_buffers->positions[2] = Vec4(center + newp0, 1.0f);
			g_buffers->positions[0] = Vec4(center + newp1, 1.0f);
			*/
			//cout << " rot " << rotation << " " << endl;
			//cout << " q " << q.x << " " << q.y << " " << q.z << " " << q.w << endl;
			//cout << " v1 " << v1.x << " " << v1.y << " " << v1.z <<  endl;
			//cout << " newp1 " << newp0.x << " " << newp0.y << " " << newp0.z << endl;
		}
		if (use_jitter) {
			float time = g_frame*g_dt;

			(Vec3&)forcefield.mPosition = Vec3((sinf(time)), 0.5f, 0.0f);
			forcefield.mRadius = (sinf(time*1.5f)*0.5f + 0.5f);
			forcefield.mStrength = jitter_strength;
			forcefield.mMode = eNvFlexExtModeVelocityChange;// eNvFlexExtModeForce;
			forcefield.mLinearFalloff = true;

			NvFlexExtSetForceFields(callback, &forcefield, 1);
		}
		else {
			float time = g_frame*g_dt;

			(Vec3&)forcefield.mPosition = Vec3((sinf(time)), 0.5f, 0.0f);
			forcefield.mRadius = (sinf(time*1.5f)*0.5f + 0.5f);
			forcefield.mStrength = 0.0f;
			forcefield.mMode = eNvFlexExtModeVelocityChange;// eNvFlexExtModeForce;//eNvFlexExtModeVelocityChange
			forcefield.mLinearFalloff = true;

			//NvFlexExtSetForceFields(callback, &forcefield, 1);
		}
		if (savetxt){
			g_cp->writeFiberTxt("D:\\Data\\cellPACK_data\\DNAPlectoneme_control_points.txt");
			savetxt = false;
		}
		if (write_output) {
			std::string pfix = "";// "_" + std::to_string(current_Exp) + "_" + std::to_string(current_run);
			g_cp->writeToBinary(pfix);
			write_output = false;
		}
	}

	virtual void PostInitialize()
	{
		//Vec3 up = Normalize(Vec3(-g_waveFloorTilt, 1.0f, 0.0f));
		//(Vec4&)g_params.planes[0] = Vec4(up.x, up.y, up.z, -8.0f);
		//g_sceneLower = Vec3(-5.0f, 0.0f, 0.0f);
		//g_sceneUpper = g_sceneLower + Vec3(10.0f, 10.0f, 5.0f);
		// free previous callback, todo: destruction phase for tests
		if (callback)
			NvFlexExtDestroyForceFieldCallback(callback);

		// create new callback
		callback = NvFlexExtCreateForceFieldCallback(g_solver);
	}

	virtual void Sync()
	{
		//NvFlexSetRigids(g_flex, g_buffers->rigidOffsets.buffer, g_buffers->rigidIndices.buffer, g_buffers->rigidLocalPositions.buffer, g_buffers->rigidLocalNormals.buffer, g_buffers->rigidCoefficients.buffer, g_buffers->rigidRotations.buffer, g_buffers->rigidTranslations.buffer, g_buffers->rigidOffsets.size() - 1, g_buffers->rigidIndices.size());
		NvFlexSetRigids(g_solver, g_buffers->rigidOffsets.buffer,
			g_buffers->rigidIndices.buffer,
			g_buffers->rigidLocalPositions.buffer,
			g_buffers->rigidLocalNormals.buffer,
			g_buffers->rigidCoefficients.buffer,
			g_buffers->rigidPlasticThresholds.buffer,
			g_buffers->rigidPlasticCreeps.buffer,
			g_buffers->rigidRotations.buffer,
			g_buffers->rigidTranslations.buffer,
			g_buffers->rigidOffsets.size() - 1,
			g_buffers->rigidIndices.size());
	}

	virtual void KeyDown(int key)
	{
		if (key == 'B')
		{
			float bombStrength = 10.0f;

			Vec3 bombPos = g_emitters[0].mPos + g_emitters[0].mDir*5.0f;
			bombPos.y -= 5.0f;

			for (int i = 0; i < int(g_buffers->velocities.size()); ++i)
			{
				Vec3 dir = Vec3(g_buffers->positions[i]) - bombPos;

				g_buffers->velocities[i] += Normalize(dir)*bombStrength*Randf(0.0, 1.0f);
			}
		}

		if (key == '/')
		{
			//save the curve data as txt files
			//x,y,z
			savetxt = !savetxt;
		}
	}

	virtual void DoGui()
	{
		if (imguiCheck("write output", write_output))
		{
			write_output = !write_output;
		}
		imguiSlider("RotationAngle", &rotationSpeed, -3.14f, 3.14f, 0.1f);
		imguiSlider("RotationRange", &angle, 0, 10, 1);
		imguiSlider("RotationStrength", &strength, 0.0f, 5.0f, 0.1f);
		if (imguiCheck("Use position", use_position))
			use_position = !use_position;
		if (imguiCheck("Use jitter", use_jitter))
			use_jitter = !use_jitter;
		imguiSlider("JitterStrength", &jitter_strength, 0.0f, 20.0f, 0.1f);
	}


	void Draw(int pass)
	{
		//cp->DrawMeshInstance(pass);
	}
};