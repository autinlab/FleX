
class Mycoplasma : public Scene
{
public:

	Mycoplasma(const char* name) : Scene(name) {}
	float main_scale = 1.0f / 200.0f;
	float scale_angle = 0.0f;
	float angle = 1.0f;
	float rotation = 0.0f;
	float rotationSpeed = 0.0f;
	bool savetxt = false;
	bool write_output = false;
	bool force_write = false;
	bool dojitter_steared = false;
	bool dojitter = false;
	bool expand_loop = false;
	bool compress_comp = false;
	bool report_overlap = false;
	bool force_field = false;

	float current_radius;
	float current_sign = 1.0f;
	float extended_radius = 0.32f;

	float x = 0;
	float y = 0;
	float z = 0;
	float s = 5.0f;
	float is = 1.0f;
	float step = 0.01f;

	void Initialize()
	{
		main_scale = g_scale;
		//run with cd D:\Dev\FLEX_EXP\nvidia_flex.git\bin\win64\;NvFlexDemoReleaseCUDA_x64.exe -force_not_center=1 -dna_persistence=1 -lod_to_use=1 -shared_memory -cellpackrecipe D:\Dev\FLEX_EXP\nvidia_flex.git\data\cellpack\\recipe_serialized.json -cellpackmodel D:\Dev\FLEX_EXP\nvidia_flex.git\data\cellpack\\model.bin
		//cd D:\Dev\FLEX_EXP\nvidia_flex.git\bin\win64\;NvFlexDemoReleaseCUDA_x64.exe -dna_persistence=10 -lod_to_use=1 -shared_memory -cellpackrecipe D:\Dev\FLEX_EXP\nvidia_flex.git\data\cellpack\\recipe_serialized.json -cellpackmodel D:\Dev\FLEX_EXP\nvidia_flex.git\data\cellpack\\model.bin
		g_buffers->rigidTranslations.resize(0);
		g_buffers->rigidRotations.resize(0);
		g_buffers->rigidCoefficients.resize(0);
		g_buffers->rigidIndices.resize(0);
		g_buffers->rigidLocalPositions.resize(0);
		g_buffers->rigidOffsets.resize(0);

		float radius = 0.1f;
		float beads_radius = 11.85f*main_scale;//5.0f*main_scale;//0.08f;//5.0f*main_scale;//smallest sphere size
		beads_radius = 34.0f/2.0f * main_scale;
		//beads_radius = ( ((34.0f / 2.0f)-5.0f) * main_scale);
		//g_params.radius = beads_radius * 2.0f;
		rotation = 0.0f;
		rotationSpeed = 0.0f;

		//g_meshcellPACK * cp;
		g_cp = new cellPACK();
		g_cp->use_rb = true;
		//g_cp->main_radius = beads_radius;//((10.0f/100.0f)*(10.0f/100.0f));
		g_cp->main_scale = main_scale;
		g_cp->mGroupCounter = 0;
		g_cp->mInstances.resize(0);
		g_cp->lodproxy_to_use = g_lod_to_use;
		g_cp->dna_persistence = g_dna_persistence;
		g_cp->force_not_center = g_force_not_center;
		g_cp->use_partners_properties = (g_use_partners_properties != 0);
		//mycoplasma experiment is
		//g_cp->loadRecipe((g_cp->mainpath + "recipes\\Mycoplasma1.7.json").c_str());
		//g_cp->loadRecipe("D:\\Data\\cellPACK_data\\Mycoplasma\\Mpn_1.0.json");
		g_cp->use_instances_mesh = false;//or true

		//g_cp->loadRecipe("../../data/cellpack/recipes/rootmyco.json");
		//g_cp->loadResults("../../data/cellpack/recipes/results.json");
		g_cp->main_radius = beads_radius;
		if (g_radius != 0.0f) {
			cout << "use radius" << g_radius * main_scale << endl;
			g_cp->overwrite_radius = false;
			g_cp->main_radius = g_radius * main_scale;
			beads_radius = g_cp->main_radius;
		}
		if (g_cp_recipe != "") {
			cout << "use recipe file " << g_cp_recipe << endl;
			g_cp->loadRecipeSerialized(g_cp_recipe, true);
		}
		else
		{
			g_cp->loadRecipeSerialized("../../data/cellpack/recipes/recipe_serialized.json", false);
		}

		beads_radius = g_cp->main_radius;
		g_params.radius = beads_radius * 2.0f;
		extended_radius = g_params.radius * 1.5f;
		//update the radius from the recipe
		//g_cp->loadResultsBinary("../../data/cellpack/recipes/results_serialized.bin");
		if (g_cp_model != "") {
			cout << "use model file " << g_cp_model << endl;
			g_cp->loadResultsBinary(g_cp_model,false, g_ignore_comp);//mpnresults_serialized.bin
		}
		else
		{
			g_cp->loadResultsBinary("../../data/cellpack/recipes/model.bin");//mpnresults_serialized.bin
		}
		//g_cp->compactInstances();
		//g_cp->placeFibers(false);//force Circle ? or check the plectnonem

		//g_cp->printBatchId();

		//g_cp->finalize();
		//mpn_results_rapid_tr.json
		//loadResults((mainpath + "recipes\\Mycoplasma1.5_mixed_pdb_fixed.json").c_str());
		
		//g_cp->haltondistribute(30);
		//cp->compactInstances();

		//for (int i = 0; i < g_cp->pnames.size(); i++)
			//cout << i << " " << g_cp->pnames[i] << endl;

		//g_params.dynamicFriction = 0.6f;
		//g_params.staticFriction = 0.35f;
		//g_params.particleFriction = 0.25f;
		g_params.dissipation = 0.0f;
		g_params.numIterations = 6;
		g_numSubsteps = 12;
		g_params.viscosity = 0.0f;
		g_params.drag = 0.0f;
		g_params.lift = 0.0f;
		g_params.numPlanes = 0;
		g_params.gravity[1] = 0.0f;//-9.f;
		g_params.damping = 30.0f;// 3.0f;

		g_params.collisionDistance = radius;// *0.5f;
		g_params.particleCollisionMargin = radius*0.25f;

		g_numExtraParticles = 1024 * 1024;

		g_drawPoints = true;
		g_wireframe = true;
		g_drawMesh = false;
		g_drawRopes = false;
		g_drawSprings = true;
		g_pause = true;
		g_warmup = false;
		//g_numSubsteps = 2;

		g_lightDistance *= 3.0f;
		
	}

	void Update()
	{
		//can we change the compartment mesh here ? and do compress/compact experiment ?
		float time = g_frame * g_dt;

		(Vec3&)forcefield.mPosition = Vec3((sinf(time)), 0.5f, 0.0f);
		forcefield.mRadius = (sinf(time*1.5f)*0.5f + 0.5f);
		forcefield.mStrength = 100.0f;
		forcefield.mMode = eNvFlexExtModeForce;
		forcefield.mLinearFalloff = true;
		
		//if (force_field && (callback != nullptr))NvFlexExtSetForceFields(callback, &forcefield, 1);

		if (compress_comp) {
			//start with scaled comp
			//then scale down
			//need to scale and stay in place
			for (int i = 0; i < g_buffers->shapeFlags.size(); ++i) {
				const int flags = g_buffers->shapeFlags[i];
				NvFlexCollisionGeometry& geo = g_buffers->shapeGeometry[i];
				int type = int(flags&eNvFlexShapeFlagTypeMask);
				Vec3 centerBeforeScaling = Vec3(0.0f);
				Vec3 centerAfterScaling = Vec3(0.0f);
				Vec3 scaleChange = Vec3(0.0f);
				Vec3 positionAdjustment = Vec3(0.0f);
				Vec3 prevPosition = Vec3(g_buffers->shapePositions[i]);
				Vec3 prev = Vec3(g_buffers->shapePositions[i]);
				if (type == eNvFlexShapeSDF) {
					geo.sdf.scale = geo.sdf.scale - step;
					std::cout << "compress? " << i << " " << geo.sdf.scale << endl;
					if (geo.sdf.scale <= 1.0f) {
						geo.sdf.scale = 1.0f;
					}
				}
				else if (type == eNvFlexShapeTriangleMesh) {
					Mesh* mesh = g_cp->comp_tri[i];
					Vector3 minExtents, maxExtents;
					mesh->GetBounds(minExtents, maxExtents);
					centerBeforeScaling = ScaleMatrix(geo.triMesh.scale[0]) * (minExtents + maxExtents) * 0.5f;
					geo.triMesh.scale[0] = geo.triMesh.scale[0] - step;
					geo.triMesh.scale[1] = geo.triMesh.scale[1] - step;
					geo.triMesh.scale[2] = geo.triMesh.scale[2] - step;
					if (geo.triMesh.scale[0] < is) {
						geo.triMesh.scale[0] = is;
						geo.triMesh.scale[1] = is;
						geo.triMesh.scale[2] = is;
					}
					centerAfterScaling = ScaleMatrix(geo.triMesh.scale[0]) * (minExtents + maxExtents) * 0.5f;
					positionAdjustment = centerBeforeScaling - centerAfterScaling;
				}
				// Apply the position adjustment to ensure the mesh remains centered.
				g_buffers->shapePrevPositions[i] = g_buffers->shapePositions[i];
				g_buffers->shapePositions[i] = Vec4(prevPosition + positionAdjustment, 0.0f);
			}
			UpdateShapes();
		}

		if (expand_loop) {
			//max
			std::cout << " distance is " << g_params.solidRestDistance << " sign " << current_sign << " radius " << g_params.radius << endl;
			if (current_radius >= extended_radius) {
				current_sign = -1.0f;
			}
			//min
			if (current_radius <= 0.10f) {
				current_sign = 1.0f;
			}
			current_radius = current_radius + 0.005f * current_sign;
			g_params.solidRestDistance = current_radius;
			g_params.radius = current_radius;
		}
		else {
		
		}
		//slowly raise radius and then decrease it ?
		// copy transforms out
		/*rotation = rotationSpeed;
		for (int i = 0; i < int(g_cp->mInstances.size()); ++i)
		{
			int ingrIndex = g_cp->mInstances[i].mMeshIndex;
			Vec3 ingrpcpal = Vec3(g_cp->iBatches[ingrIndex].pcpalVectorx, g_cp->iBatches[ingrIndex].pcpalVectory, g_cp->iBatches[ingrIndex].pcpalVectorz);
			Quat q = QuatFromAxisAngle(ingrpcpal, rotation);//should start at 0
			//g_cp->mInstances[i].mTranslation = g_buffers->rigidTranslations[i];
			//g_cp->mInstances[i].mRotation = g_buffers->rigidRotations[i]*q;
			//g_buffers->rigidRotations[i] = g_cp->mInstances[0].mRotation;
			//rotate the particles instead ?
			NvFlexExtAsset* asset = g_cp->iBatches[ingrIndex].mAsset;
			for (int j = 0; j < asset->numParticles; ++j)
			{
				Vec3 localPos = Vec3(&asset->particles[j * 4]) - Vec3(&asset->shapeCenters[0]);//local position of the proxy
				Vec3 rpos = Rotate(q, localPos);
				rpos = Rotate(g_buffers->rigidRotations[i], rpos);
				g_buffers->positions[g_cp->mInstances[i].mParticleOffset + j] = Vec4(g_buffers->rigidTranslations[i] + rpos, 1.0f);
				//g_buffers->velocities[g_cp->mInstances[i].mParticleOffset + j] = Vec3(0, 0, 0);
			}
		}*/
		/*if ((g_frame % 5) == 0)
		{
			//dojitter_strength
			//dojitter_steared
			//dojitter_biased_strength
			if (dojitter) g_cp->jitter(dojitter_steared);
		}*/
		if (report_overlap) {
			int count = g_cp->countOverlap();
			std::cout << " total count is " << count  << endl;
		}
		if (dojitter) {
			g_cp->jitter(dojitter_steared, force_field);
			//g_cp->RandomDiffusion();
		}
		if (savetxt){
			std::string pfix = "_" + std::to_string(g_frame);
			g_cp->writeModelPDB(pfix,1.0/100.0);
			savetxt = false;
		}
		if (write_output) {
			std::string pfix = "_" + std::to_string(g_frame);// "_" + std::to_string(current_Exp) + "_" + std::to_string(current_run);
			g_cp->writeToBinary(pfix);
			if (!force_write) write_output = false;
		}
	}

	virtual void PostInitialize()
	{
		cout <<  "PostInitialize begin " << endl;
		// free previous callback, todo: destruction phase for tests
		if (callback != nullptr) {
			cout << "NvFlexExtDestroyForceFieldCallback" << endl;
			try {
				//NvFlexExtDestroyForceFieldCallback(callback);
			}
			catch (std::system_error& e) {
				std::cerr << e.code().message() << std::endl;
			}
		}

		// create new callback
		cout << "create new callback" << endl;
		//callback = NvFlexExtCreateForceFieldCallback(g_solver);
		//g_sceneLower = Vec3(-5.0f, 0.0f, 0.0f);
		//g_sceneUpper = g_sceneLower + Vec3(10.0f, 10.0f, 5.0f);
		cout << "PostInitialize end" << endl;
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
		if (key == 'x')
		{
			std::cout << "key x" << endl;
			std:string pfix = "";// "_" + std::to_string(current_Exp) + "_" + std::to_string(current_run);
			g_cp->writeToBinary(pfix);
		}
		if (key == 'j')
		{
			std::cout << "key j" << endl;
			g_cp->RandomDiffusion();
		}
	}

	virtual void DoGui()
	{
		imguiSlider("expand_radius", &extended_radius, 0.0f, g_cp->main_radius * 4.0f, 0.01f);
		
		if (imguiCheck("force_field", force_field))
		{
			force_field = !force_field;
		}
		if (imguiCheck("compress_comp", compress_comp)) {
			if (!compress_comp) {
				//start compression
				Vec3 centerBeforeScaling = Vec3(0.0f);
				Vec3 centerAfterScaling = Vec3(0.0f);
				Vec3 scaleChange = Vec3(0.0f);
				Vec3 positionAdjustment = Vec3(0.0f);
				// should correspond to comp_mesh
				for (int i = 0; i < g_buffers->shapeFlags.size(); ++i) {
					
					const int flags = g_buffers->shapeFlags[i];
					NvFlexCollisionGeometry& geo = g_buffers->shapeGeometry[i];
					Vec3 tocenter = Vec3(0.0f);
					Vec3 center1 = Vec3(0.0f);
					Vec3 center2 = Vec3(0.0f);
					Vec3 localLower;
					Vec3 localUpper;
					Vec3 prev = Vec3(g_buffers->shapePositions[i]);
					int type = int(flags&eNvFlexShapeFlagTypeMask);
					Vec3 prevPosition = Vec3(g_buffers->shapePositions[i]);
					std::cout << i << " prevPosition " << prevPosition[0] << " " << prevPosition[1] << " " << prevPosition[2] << endl;
					if (type == eNvFlexShapeSDF) {
						geo.sdf.scale = 5.0f;
					}
					if (type == eNvFlexShapeTriangleMesh) {
						Mesh* mesh = g_cp->comp_tri[i];
						Vector3 minExtents, maxExtents;
						mesh->GetBounds(minExtents, maxExtents);
						centerBeforeScaling = ScaleMatrix(geo.triMesh.scale[0]) * (minExtents + maxExtents) * 0.5f;
						std::cout << i << " " << s << " centerBeforeScaling " << centerBeforeScaling[0] << " " << centerBeforeScaling[1] << " " << centerBeforeScaling[2] << endl;

						geo.triMesh.scale[0] = s;
						geo.triMesh.scale[1] = s;
						geo.triMesh.scale[2] = s;
						
						//NvFlexGetTriangleMeshBounds(g_flexLib, geo.triMesh.mesh, centerAfterScaling, scaleChange); // Reuse to find the new bounds and center.
						//centerAfterScaling = (scaleChange - centerAfterScaling) * 0.5f + centerAfterScaling;
						centerAfterScaling = ScaleMatrix(s) * centerBeforeScaling;
						std::cout << i << " centerAfterScaling " << centerAfterScaling[0] << " " << centerAfterScaling[1] << " " << centerAfterScaling[2] << endl;

						//NvFlexGetTriangleMeshBounds(g_flexLib, geo.triMesh.mesh, localLower, localUpper);
						//center1 = (localUpper - localLower)*0.5f;//local
						positionAdjustment = (centerBeforeScaling - centerAfterScaling);
						std::cout << i << " positionAdjustment " << positionAdjustment[0] << " " << positionAdjustment[1] << " " << positionAdjustment[2] << endl;
						//localLower *= Vec3(geo.triMesh.scale);
						//localUpper *= Vec3(geo.triMesh.scale);
						//center2 = (localUpper - localLower) * 0.5f - center1;
						//std::cout << i <<" center2 " << center2[0] << " " << center2[1] << " " << center2[2] << endl;
						//x = 0.0f;//center2[0];// *3.2f;//for segmented shape
						//y = 0.0f;//-center2[1];// *1.5f;
						//z = 0.0f;// -center2[2];
					}
					g_buffers->shapePositions[i] = Vec4(prevPosition + positionAdjustment, 0.0f);
					g_buffers->shapePrevPositions[i] = g_buffers->shapePositions[i];
					// g_buffers->shapePrevPositions[i] = Vec4(prevPosition, 0.0f);
					//std::cout << i << " a pos " << g_buffers->shapePositions[i][0] << " " << g_buffers->shapePositions[i][1] << " " << g_buffers->shapePositions[i][2] << endl;
					//g_buffers->shapePositions[i] = Vec4(Vec3(x,y,z), 0.0f);
					//g_buffers->shapePrevPositions[i] = Vec4(prev, 0.0f); /// Vec3(geo.triMesh.scale)
					std::cout << i << " b pos " << g_buffers->shapePositions[i][0] << " " << g_buffers->shapePositions[i][1] << " " << g_buffers->shapePositions[i][2] << endl;
				}
				UpdateShapes();
			}
			else {
				//stop compress
				for (int i = 0; i < g_buffers->shapeFlags.size(); ++i) {
					const int flags = g_buffers->shapeFlags[i];
					NvFlexCollisionGeometry& geo = g_buffers->shapeGeometry[i];
					int type = int(flags&eNvFlexShapeFlagTypeMask);
					if (type == eNvFlexShapeSDF) {
						geo.sdf.scale = main_scale;
					}
					if (type == eNvFlexShapeTriangleMesh) {
						geo.triMesh.scale[0] = 1.0f;
						geo.triMesh.scale[1] = 1.0f;
						geo.triMesh.scale[2] = 1.0f;
					}
					g_buffers->shapePositions[i] = Vec4(0.0f);
					g_buffers->shapePrevPositions[i] = Vec4(0.0f);
				}
				UpdateShapes();
			}
			compress_comp = !compress_comp;
		}
		if (imguiCheck("expand_loop", expand_loop))
		{
			float mradius = 0.50f;
			expand_loop = !expand_loop;
			if (expand_loop) {
				current_sign = 1.0f;
				current_radius = g_cp->main_radius * 2.0f;
				//g_params.radius = mradius;
				//g_params.collisionDistance = mradius;
				//g_params.solidRestDistance = mradius;
			}
			else {
				current_radius = g_cp->main_radius * 2.0f;
				g_params.radius = g_cp->main_radius * 2.0f;
				// g_params.collisionDistance = g_cp->main_radius * 2.0f;
				g_params.solidRestDistance = g_cp->main_radius * 2.0f;
			}
		}
		if (imguiCheck("motion", dojitter))
		{
			dojitter = !dojitter;
		}
		if (imguiCheck("allparticle", dojitter_steared))
		{
			dojitter_steared = !dojitter_steared;
		}
		if (imguiCheck("toward center", g_cp->dojitter_biased))
		{
			g_cp->dojitter_biased = !g_cp->dojitter_biased;
		}
		
		if (imguiCheck("report overlap", report_overlap)) {
			report_overlap = !report_overlap;
		}
		imguiSlider("miniscale", &is, 0.1f, 1.0, 0.01f);
		imguiSlider("maxiscale", &s, 1.0f, 10.0, 0.01f);
		imguiSlider("step", &step, 0.0001f, 1.0f, 0.001f);
		/*
		imguiSlider("x off", &x, -150.0f, 150.0, 0.1f);
		imguiSlider("y off", &y, -150.0f, 150.0, 0.1f);
		imguiSlider("z off", &z, -150.0f, 150.0, 0.1f);
		for (int i = 0; i < g_buffers->shapeFlags.size(); ++i) {
			g_buffers->shapePrevPositions[i] = g_buffers->shapePositions[i];
			g_buffers->shapePositions[i] = Vec4(x,y,z, 0.0f);
		}
		*/
		imguiSlider("motion strength", &g_cp->dojitter_strength, 0.0f, 1.0, 0.001f);
		imguiSlider("motion bias strength", &g_cp->dojitter_biased_strength, 0.0f, 2.0, 0.0f);
		imguiSlider("comp_radius", &g_cp->comp_radius, 0.1f, 10.0f, 0.1f);

		if (imguiCheck("write output", write_output))
		{
			write_output = !write_output;
		}
		if (imguiCheck("force write output", force_write))
		{
			force_write = !force_write;
		}
		if (imguiCheck("write output txt", savetxt))
		{
			savetxt = !savetxt;
		}
		imguiSlider("RotationSpeed", &rotationSpeed, -20.0f, 20.0f, 0.1f);
		imguiSlider("RotationAngle", &angle, -20.0f, 20.0f, 0.1f);
	}

	void Draw(int pass)
	{
		g_cp->DrawMeshInstance(pass);
	}

	NvFlexExtForceField forcefield;

	NvFlexExtForceFieldCallback* callback;
};
