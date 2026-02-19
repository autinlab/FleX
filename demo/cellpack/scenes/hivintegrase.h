class HIVIntegrase : public Scene
{
public:

	HIVIntegrase(const char* name) : Scene(name) {}
	
	cellPACK * cp; 
	float main_scale;
	struct Integrase
	{
		int nb_site_occupied;
		int ccd_occupied[2];
		int ccd_ctdind[2];
		int ctd_occupied[2];
		int ctd_ccdind[2];
		int spring_offset[2];
		float mini_distance[2];
		int mini_indices[2];
	};

	std::vector<Integrase> mIntegrase;
	std::vector<int> free_mIntegrase_CCD;
	std::vector<int> free_mIntegrase_CTD;

	std::vector<std::vector<int> > ccd_ctd;
	std::vector<std::vector<int> > ctd_ccd;
	std::vector<int> ccd_ctd_1;
	std::vector<int> ccd_ctd_2;
	std::vector<int> ctd_ccd_1;
	std::vector<int> ctd_ccd_2;
	std::vector<int> strings_indices;
	std::vector<float> prev_strings_distances;//at 1 ms
	std::vector<std::vector<int> > binding_in;
	std::vector<std::vector<int> > binding_out;

	bool write_output = false;
	bool dojitter_steared = true;
	bool dojitter = true;
	bool dojitter_biased = true;
	bool use_threshold_binding = true;
	float dojitter_strength = 1.0f;
	float dojitter_biased_strength = 1.0f;
	float threshold_binding;
	bool dosimulation = false;

	std::vector<int> record_free;
	std::vector<float> record_size;

	int stop_criterion_time = 100; // how many last frame to store and check
	int count_stop = 0;
	int nexp = 0;
	int current_Exp = 0; // how many experiment?
	int nrun = 5;
	int current_run = 0;
	std::map <int, std::vector<float>> exp_params;//carry other the experiments parameters
	struct Node
	{
		int nb_site_occupied;
		std::vector<int> binding_in_occupied;
		std::vector<int> binding_in_ind;
		std::vector<int> binding_out_occupied;
		std::vector<int> binding_out_ind;
		std::vector<int> spring_offset;
		std::vector<float> mini_distance;
		std::vector<int> mini_indices;
	}; 
	std::vector<Node> mNodes;
	//float dojitter_strength = 1.0f;
	//float dojitter_biased_strength = 1.0f;
	//float threshold_binding;

	virtual void DoGui()
	{
		if (imguiCheck("write output", write_output))
		{
			write_output = !write_output;
		}

		if (imguiCheck("dosimulation", dosimulation))
		{
			dosimulation = !dosimulation;
		}
		if (imguiCheck("motion", dojitter))
		{
			dojitter = !dojitter;
		}
		if (imguiCheck("steared", dojitter_steared))
		{
			dojitter_steared = !dojitter_steared;
			//if (dojitter_steared) dojitter_strength = 1.0f;
			//else dojitter_strength = 0.01f;
		}
		if (imguiCheck("center biased", dojitter_biased))
		{
			dojitter_biased = !dojitter_biased;
		}

		if (imguiCheck("use binding thr", use_threshold_binding))
		{
			use_threshold_binding = !use_threshold_binding;
			if (!use_threshold_binding)
				threshold_binding = 999999.9;
		}
		imguiSlider("motion strength", &dojitter_strength, 0.0f, 2.0, 0.1f);
		imguiSlider("center biased strength", &dojitter_biased_strength, 0.0f, 2.0f, 0.1f);
		imguiSlider("binding distance", &threshold_binding, 0.0f, 20.0f, 0.01f);

		if (imguiButton("Export node as PDB")){
			exportPDBRigidTransform(current_Exp, current_run);
			//writeSoftTransform(current_Exp, current_run, 0, false);	//#the beads
			//writeRigidTransform(current_Exp, current_run, false);		//the shape rigid bodu
		}
		//if (imguiButton("Export PDB all beads"))
		//	exportPDBRigidTransform(current_Exp, current_run);
	}

	virtual void KeyDown(int key)
	{
		if (key == 'n')
		{
			std::cout << "key n" << endl;
			addOneRandomBinding();
		}

		if (key == 'b')
		{
			std::cout << "key b" << endl;
			reportBinding();
			checkDistance(false);
			checkRadiusAggregate();
		}
		if (key == 'S')
		{
			std::cout << "key S" << endl;
			//writeRigidTransform();
			//writeSoftTransform(0,0,0);
			//exportPDBRigidTransform(current_Exp, current_run);
		}
	}

	virtual void addOneRandomBinding(){
		bool doit = true;
		if (free_mIntegrase_CCD.size() >= 1){
			if (free_mIntegrase_CTD.size() >= 1) {
				int mind1 = Rand() % free_mIntegrase_CCD.size();
				int mind2 = Rand() % free_mIntegrase_CTD.size();
				int id1 = checkOccupationCCD(free_mIntegrase_CCD[mind1]);
				int id2 = checkOccupationCTD(free_mIntegrase_CTD[mind2]);
				//cout << id1 << " " << id2 << endl;
				if (free_mIntegrase_CCD[mind1] == free_mIntegrase_CTD[mind2]) doit = false;
				if (id1 == -1){
					doit = false;
					free_mIntegrase_CCD.erase(free_mIntegrase_CCD.begin() + mind1);
				}
				if (id2 == -1){
					doit = false;
					free_mIntegrase_CTD.erase(free_mIntegrase_CTD.begin() + mind2);
				}
				if (doit) {
					//cout << "ADD " << mind1 << " " << free_mIntegrase_CCD[mind1] << " " << id1 << " " << mind2 << " " << free_mIntegrase_CTD[mind2] << " " << id2 << endl;
					boundTwoIds(free_mIntegrase_CCD[mind1], free_mIntegrase_CTD[mind2], id1, id2);
				}
			}
		}
	}

	virtual void addOneBinding(int mind1, int mind2, int ccd_id1, int ctd_id2){
		bool doit = true;
		//int id1 = checkOccupationCCD(free_mIntegrase_CCD[mind1]);// Rand() % 2;
		//int id2 = checkOccupationCTD(free_mIntegrase_CTD[mind2]);
		if ((mIntegrase[mind1].ccd_occupied[ccd_id1] != -1) || (mIntegrase[mind2].ctd_occupied[ctd_id2] != -1)) doit = false;
		if (doit) {
			//cout << "ADD " << mind1 << " " << free_mIntegrase_CCD[mind1] << " " << id1 << " " << mind2 << " " << free_mIntegrase_CTD[mind2] << " " << id2 << endl;
			boundTwoIds(free_mIntegrase_CCD[mind1], free_mIntegrase_CTD[mind2], ccd_id1, ctd_id2);
		}
	}

	virtual void boundTwo(int i, int j){
		//randomly pick one acceptor//one donor
		int id1 = Rand() % 2;
		int id2 = Rand() % 2;
		boundTwoIds(i, j, id1, id2);
	}

	int checkOccupationCCD(int i){
		//cout << "check " << i << " " << mIntegrase[i].ccd_occupied[0] << " " << mIntegrase[i].ccd_occupied[1] << endl;
		if (mIntegrase[i].ccd_occupied[0] != -1 && mIntegrase[i].ccd_occupied[1] != -1)
		{
			return -1;
		}
		else if (mIntegrase[i].ccd_occupied[0] != -1){
			return 1;
		}
		else if (mIntegrase[i].ccd_occupied[1] != -1){
			return 0;
		}
		else //if ((!mIntegrase[i].ccd_occupied[0]) && (!mIntegrase[i].ccd_occupied[1]))
		{
			return Rand() % 2;
		}
	}

	int getOccupationCCD(int i){
		//cout << "check " << i << " " << mIntegrase[i].ccd_occupied[0] << " " << mIntegrase[i].ccd_occupied[1] << endl;
		int count = 0;
		if (mIntegrase[i].ccd_occupied[0] != -1)
			count++;
		if (mIntegrase[i].ccd_occupied[1] != -1)
			count++;
		if (mIntegrase[i].ctd_occupied[0] != -1)
			count++;
		if (mIntegrase[i].ctd_occupied[1] != -1)
			count++;
		return count;
	}

	int checkOccupationCTD(int i){
		if ((mIntegrase[i].ctd_occupied[0] != -1) && (mIntegrase[i].ctd_occupied[1] != -1))
		{
			return -1;
		}
		else if (mIntegrase[i].ctd_occupied[0] != -1){
			return 1;
		}
		else if (mIntegrase[i].ctd_occupied[1] != -1){
			return 0;
		}
		else //if ((!mIntegrase[i].ctd_occupied[0]) && (!mIntegrase[i].ctd_occupied[1]))
		{
			return Rand() % 2;
		}
	}

	virtual int boundTwoIds(int i, int j, int id1, int id2){
		int poffseti = cp->mInstances[i].mParticleOffset;
		int poffsetj = cp->mInstances[j].mParticleOffset;
		int sp_offset = strings_indices.size();
		if (mIntegrase[i].ccd_occupied[id1] != -1) return -1;
		if (mIntegrase[j].ctd_occupied[id2] != -1) return -1;
		//if (mIntegrase[i].ccd_occupied[1 - id1] == j || mIntegrase[j].ctd_occupied[1 - id2] == i) return -1;
		for (int n = 0; n < ccd_ctd[id1].size(); n++) {
			for (int m = 0; m < ctd_ccd[id2].size(); m++) {
				//cout << "bound " << id1 << " " << n << " " << poffseti + ccd_ctd[id1][n] << " " << id2 << " " << m << " " << poffsetj + ctd_ccd[id2][m] << " " << g_params.mRadius << " " << g_springLengths.size() << " " g_springIndices << endl;
				strings_indices.push_back(g_buffers->springLengths.size());
				CreateSpringInter(poffseti + ccd_ctd[id1][n], poffsetj + ctd_ccd[id2][m], 0.01f, 0.0f, g_params.radius); //float give = 0.0f, float length = 0.0f
			}
		}
		mIntegrase[i].ccd_occupied[id1] = j;
		mIntegrase[j].ctd_occupied[id2] = i;
		mIntegrase[i].ccd_ctdind[id1] = id2;
		mIntegrase[j].ctd_ccdind[id2] = id1;
		mIntegrase[i].spring_offset[id1] = sp_offset;
		//flexSetSprings(g_flex, &g_springIndices[0], &g_springLengths[0], &g_springStiffness[0], g_springLengths.size(), eFlexMemoryHost);
		//cout << "bound " << i << " " << j << " " << id1 << " " << id2 << " " << mIntegrase[i].ccd_occupied[id1] << " " << mIntegrase[j].ctd_occupied[id2] << endl;
		return 1;
	}

	virtual void setupBoundsAndLinker(){
		//go over all instances, bind to every other instances
		//indicesBounds
		for (int i = 0; i < int(cp->mInstances.size()); ++i)
		{
			//IngredientInstance& inst = cp->mInstances[i];
			Integrase integrase = Integrase();
			for (int i = 0; i < 2; i++){
				integrase.ccd_ctdind[i] = -1;
				integrase.ccd_occupied[i] = -1;
				integrase.ctd_occupied[i] = -1;
				integrase.ccd_ctdind[i] = -1;
				integrase.ctd_ccdind[i] = -1;
				integrase.spring_offset[i] = -1;
				integrase.mini_distance[i] = 99999.0f;
				integrase.mini_indices[i] = -1;
			}
			integrase.nb_site_occupied = 0;
			mIntegrase.push_back(integrase);
			free_mIntegrase_CCD.push_back(i);
			free_mIntegrase_CTD.push_back(i);
		}
		cout << " nb integrase " << mIntegrase.size() << endl;
	}

	virtual void reportBinding(bool report = true){
		std::vector<int> reports_count = { 0, 0, 0 };
		bool bounded = false;
		bool founded = false;
		int bound_order = 0;
		int count_free = 0;
		int unsatisfied = 0;
		for (int i = 0; i < int(cp->mInstances.size()) - 1; ++i)
		{
			NvFlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			int poffseti = cp->mInstances[i].mParticleOffset;
			if (report) {
				//cout << i << " CCD " << mIntegrase[i].ccd_occupied[0] << " " << mIntegrase[i].ccd_occupied[1];//miniD >?
				//cout << " CTD " << mIntegrase[i].ctd_occupied[0] << " " << mIntegrase[i].ctd_occupied[1] << endl;
			}
			if (mIntegrase[i].ccd_occupied[0] == -1) count_free++;
			if (mIntegrase[i].ccd_occupied[1] == -1) count_free++;
			continue;
			bound_order = 0;
			for (int j = i + 1; j < int(cp->mInstances.size()); j++) {
				//build spring ccd-ctd
				if (bound_order == 2) continue;
				int poffsetj = cp->mInstances[j].mParticleOffset;
				founded = false;
				for (int o = 0; o < 2; o++){
					founded = false;
					for (int p = 0; p < 2; p++){
						if (founded) continue;
						for (int n = 0; n < ccd_ctd[o].size(); n++) {
							if (founded) continue;
							for (int m = 0; m < ctd_ccd[p].size(); m++) {
								if (founded) continue;
								float D = Length(Vec3(g_buffers->positions[poffseti + ccd_ctd[o][n]]) - Vec3(g_buffers->positions[poffsetj + ctd_ccd[p][m]]));
								//cout << D << endl;
								if (D <= g_params.radius) {
									//cout << " bounded " << D << " " << i << " " << j << " " << o << " " << n << " " << p << " " << m << endl;
									founded = true;
								}
							}
						}
					}
					if (founded) {
						bound_order += 1;
						founded = false;
					}
				}
			}
			if (bound_order>0){
				for (int k = 0; k < bound_order; k++)
					reports_count[k] += 1;
			}
			else {
				reports_count[2] += 1;
			}
		}
		record_free.push_back(count_free);
		if (report)
			cout << count_free << " " << reports_count[0] << " " << reports_count[1] << " " << reports_count[2] << endl;
	}

	virtual float getMinDistances(int i, int j, int ccd, int ctd){
		float miniD = 999999.9f;
		int poffseti = cp->mInstances[i].mParticleOffset;
		int poffsetj = cp->mInstances[j].mParticleOffset;
		for (int n = 0; n < ccd_ctd[ccd].size(); n++) {
			for (int m = 0; m < ctd_ccd[ctd].size(); m++) {
				float D = Length(Vec3(g_buffers->positions[poffseti + ccd_ctd[ccd][n]]) - Vec3(g_buffers->positions[poffsetj + ctd_ccd[ctd][m]]));
				if (D < miniD) {
					miniD = D;
				}
			}
		}
		return miniD;
	}

	virtual int checkDistance(bool break_spring = false){
		//for every occupied pair, do we satsify the distance.  If not  remove the spring ?
		//should have store the spring indice for each site as well, or just the offset indice as we know the number of spring.
		int unsatisfied_count = 0;
		int nspring = ccd_ctd[0].size()*ctd_ccd[0].size() + 1;
		for (int i = 0; i < mIntegrase.size(); i++){
			//check if the ccd occupied are actually close
			for (int o = 0; o< 2; o++){
				if (mIntegrase[i].ccd_occupied[o] != -1){
					//mesure distance with the other ctd
					int j = mIntegrase[i].ccd_occupied[o];
					int p = mIntegrase[i].ccd_ctdind[o];
					float D = getMinDistances(i, j, o, p);
					if (D < mIntegrase[i].mini_distance[o]) {
						mIntegrase[i].mini_distance[o] = D;
						mIntegrase[i].mini_indices[o] = j;
					}
					if (D > g_params.radius*2.0f){
						unsatisfied_count++;
						if (break_spring){
							//remove those strings_indices from g_springIndices
							int indice_sp = mIntegrase[i].spring_offset[o];
							if (strings_indices[indice_sp] == -1) continue;
							cout << "remove " << indice_sp << " " << strings_indices[indice_sp] << " " << strings_indices[indice_sp] * 2 << endl;
							//g_buffers->springIndices.erase(g_springIndices.begin() + strings_indices[indice_sp] * 2, g_springIndices.begin() + strings_indices[indice_sp] * 2 + nspring * 2);
							//g_buffers->springLengths.erase(g_springLengths.begin() + strings_indices[indice_sp], g_springLengths.begin() + strings_indices[indice_sp] + nspring);
							//g_buffers->springStiffness.erase(g_springStiffness.begin() + strings_indices[indice_sp], g_springStiffness.begin() + strings_indices[indice_sp] + nspring);
							strings_indices[indice_sp] = -1;
							//change the indices ?
							for (int i = indice_sp + 1; i < strings_indices.size(); i++){
								strings_indices[i] -= nspring;
							}
							//strings_indices.erase(strings_indices.begin()+indice_sp);
							//unbind
							mIntegrase[i].ccd_occupied[o] = -1;
							mIntegrase[i].ccd_ctdind[o] = -1;
							mIntegrase[j].ctd_occupied[p] = -1;
						}
					}
				}
			}

		}
		NvFlexSetSprings(g_solver, g_buffers->springIndices.buffer, g_buffers->springLengths.buffer, g_buffers->springStiffness.buffer, g_buffers->springLengths.size());
		cout << "unsatisfied " << unsatisfied_count << endl;
		//break spring ?
		return unsatisfied_count;
	}

	virtual int oneDomain(const float* all_particles, int numParticles, std::vector<int> indices, std::vector<int>& outSpringIndices,
		std::vector<float>& outSpringLengths, std::vector<float>& outSpringStiffness, float radius, float stiffness = 1.0f, int offset = 0)
	{

		int count = 0;
		std::vector<Vec3> particles;
		for (int i = 0; i < numParticles; ++i)
		{
			Vec3 localPos = Vec3(&all_particles[(indices[i] + offset) * 4]);
			particles.push_back(localPos);
		}
		for (int l = 0; l < indices.size() - 1; l++){
			for (int m = l + 1; m < indices.size(); m++){
				outSpringIndices.push_back(indices[l] + offset);
				outSpringIndices.push_back(indices[m] + offset);
				outSpringLengths.push_back(Length(Vec3(particles[indices[l] + offset]) - Vec3(particles[indices[m] + offset])));
				outSpringStiffness.push_back(stiffness);
				count++;
			}
		}
		return count;
	}

	virtual void overwriteSpringsBackup(){
		/* using domain and linker information overwrite spring network of the assets.
		#domain CCD1  - 6,  7, 11, 12, 14, 15, 19, 22, 24, 26, 33, 35, 37, 40, 47
		#Linkera1     - 17,21,28,32,40
		#domain CTD1  - 1,2,3,4,5,9,10,16,18,31,39,43,44
		#Linkerb1     - 45,46
		#domain NTD1  - 8,13,20,23,25,27,29,30,34,36,38,40,41,42,48,49
		#domain CCD2  - 56,57,61,62,64,65,69,72,74,76,83,85,87,90,97
		#linkera2+50  - 67,71,78,82,90
		#domain CTD2  - 51, 52, 53, 54, 55, 59, 60, 66, 68, 81, 89, 93, 94
		#linkerb2+50  - 95,96
		#domain NTD2  - 58, 63, 70, 73, 75, 77, 79, 80, 84, 86, 88, 90, 91, 92, 98, 99
		*/
		float radius = g_params.radius*0.5f;
		float stiffness = 1.0f;
		std::vector<int> CTD1 = { 1, 2, 3, 4, 5, 9, 10, 16, 18, 28, 31, 39, 43, 44 };
		std::vector<int> CCD1 = { 0, 6, 7, 11, 12, 14, 15, 19, 22, 24, 26, 33, 35, 37, 40, 46, 47 };
		std::vector<int> NTD1 = { 8, 13, 20, 23, 25, 27, 29, 30, 34, 36, 38, 41, 42, 45, 48, 49 };
		std::vector<int> LinkerA = { 32, 28 };/// { 40, 21, 17, 32, 28 };
		std::vector<int> LinkerB = { 32 + 50, 28 + 50 };// { 45, 46 };
		std::vector<int> dimer = { 8, 75, 25, 58 };// { 11, 12, 14, 15, 22, 26, 33, 47, 8, 23, 25, 34, 49 };
		std::vector<int> CCD = { 0, 6, 7, 11, 12, 14, 15, 19, 22, 24, 26, 33, 35, 37, 40, 46, 47,
			50, 56, 57, 61, 62, 64, 65, 69, 72, 74, 76, 83, 85, 87, 90, 96, 97 };
		std::vector<int> NTD = { 8, 13, 20, 23, 25, 27, 29, 30, 34, 36, 38, 41, 42, 45, 48, 49,
			58, 63, 70, 73, 75, 77, 79, 80, 84, 86, 88, 91, 92, 95, 98, 99 };
		std::vector<int> CCDNTD = { 0, 6, 7, 11, 12, 14, 15, 17, 19, 21, 22, 24, 26, 32, 33, 35, 37, 40, 46, 47,
			50, 56, 57, 61, 62, 64, 65, 69, 72, 74, 76, 83, 85, 87, 90, 96, 97, 8, 13, 20, 23, 25, 27, 29, 30, 34, 36, 38, 41, 42, 45, 48, 49,
			58, 63, 70, 73, 75, 77, 79, 80, 84, 86, 88, 91, 92, 95, 98, 99, 71, 67, 82
		};

		NvFlexExtAsset* asset = cp->iBatches[0].mAsset;
		std::vector<int> springIndices;
		std::vector<float> springLengths;
		std::vector<float> springStiffness;
		int numLinks = 0;
		std::vector<int> springIndicestmp;
		std::vector<float> springLengthstmp;
		std::vector<float> springStiffnesstmp;
		// create links between particles

		int numLinksTmp = CreateLinksLocal(asset->particles, CTD1.size(),
			CTD1,
			springIndices, springLengths, springStiffness,
			radius, stiffness);
		numLinks += numLinksTmp;
		/*
		numLinksTmp = CreateLinksLocal(asset->mParticles, CCD.size(), CCD,
		springIndices, springLengths, springStiffness,
		radius, stiffness);
		numLinks += numLinksTmp;
		numLinksTmp = CreateLinksLocal(asset->mParticles, NTD1.size(), NTD1,
		springIndices, springLengths, springStiffness,
		radius, stiffness);
		numLinks += numLinksTmp;
		*/
		numLinksTmp = CreateLinksLocal(asset->particles, 
			CCDNTD.size(), CCDNTD,
			springIndices, springLengths, springStiffness,
			radius, stiffness);
		numLinks += numLinksTmp;


		//vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
		//create spring for Linker 
		for (int l = 0; l < LinkerA.size() - 1; l++){
			springIndices.push_back(LinkerA[l]);
			springIndices.push_back(LinkerA[l + 1]);
			springLengths.push_back(Length(Vec3(&asset->particles[LinkerA[l] * 4]) - Vec3(&asset->particles[LinkerA[l + 1] * 4])));
			springStiffness.push_back(stiffness);
			numLinks++;
		}
		/*
		for (int l = 0; l < LinkerB.size() - 1; l++){
		springIndices.push_back(LinkerB[l]);
		springIndices.push_back(LinkerB[l + 1]);
		springLengths.push_back(Length(Vec3(&asset->mParticles[LinkerB[l] * 4]) - Vec3(&asset->mParticles[LinkerB[l + 1] * 4])));
		springStiffness.push_back(stiffness);
		numLinks++;
		}
		*/
		numLinksTmp = CreateLinksLocal(asset->particles, CTD1.size(), CTD1,
			springIndices, springLengths, springStiffness,
			radius, stiffness, 50);
		numLinks += numLinksTmp;

		for (int l = 0; l < LinkerA.size() - 1; l++){
			springIndices.push_back(LinkerA[l] + 50);
			springIndices.push_back(LinkerA[l + 1] + 50);
			springLengths.push_back(Length(Vec3(&asset->particles[(LinkerA[l] + 50) * 4]) - Vec3(&asset->particles[(LinkerA[l + 1] + 50) * 4])));
			springStiffness.push_back(stiffness);
			numLinks++;
		}
		/*numLinksTmp = CreateLinksLocal(asset->mParticles, NTD1.size(), NTD1,
		springIndices, springLengths, springStiffness,
		radius, stiffness, 50);
		numLinks += numLinksTmp;
		numLinksTmp = CreateLinksLocal(asset->mParticles, CCD1.size(), CCD1,
		springIndices, springLengths, springStiffness,
		radius, stiffness, 50);
		numLinks += numLinksTmp;
		for (int l = 0; l < LinkerB.size() - 1; l++){
		springIndices.push_back(LinkerB[l] + 50);
		springIndices.push_back(LinkerB[l + 1] + 50);
		springLengths.push_back(Length(Vec3(&asset->mParticles[(LinkerB[l] + 50) * 4]) - Vec3(&asset->mParticles[(LinkerB[l + 1] + 50) * 4])));
		springStiffness.push_back(stiffness);
		numLinks++;

		for (int l = 0; l < dimer.size() - 1; l++){
		springIndices.push_back(dimer[l]);
		springIndices.push_back(dimer[l+1]);
		springLengths.push_back(Length(Vec3(&asset->mParticles[(dimer[l]) * 4]) - Vec3(&asset->mParticles[(dimer[l+1]) * 4])));
		springStiffness.push_back(0.5f);
		numLinks++;
		l++;
		}

		}
		*/



		//DIMER LINK CCD/CCD and NTD/NTD


		// assign links
		if (numLinks)
		{
			asset->springIndices = new int[numLinks * 2];
			memcpy(asset->springIndices, &springIndices[0], sizeof(int)*springIndices.size());

			asset->springCoefficients = new float[numLinks];
			memcpy(asset->springCoefficients, &springStiffness[0], sizeof(float)*numLinks);

			asset->springRestLengths = new float[numLinks];
			memcpy(asset->springRestLengths, &springLengths[0], sizeof(float)*numLinks);

			asset->numSprings = numLinks;
		}

		//overwrite teh shapes
		int numClusters = 5;

		/*asset->mNumShapeIndices = 3;
		asset->mShapeIndices
		asset->mShapeOffsets
		asset->mShapeCenters
		asset->mShapeCoefficients
		*/
		std::vector<int> clusterIndices;
		std::vector<int> clusterOffsets;
		std::vector<Vec3> clusterPositions;
		std::vector<float> clusterCoefficients;
		float clusterStiffness = 1.0f;
		//CCDNTD
		//CTD1
		//CTD2
		Vec3 center = Vec3(0, 0, 0);
		for (int i = 0; i < CCDNTD.size(); ++i)
		{
			clusterIndices.push_back(CCDNTD[i]);
			center += Vec3(asset->particles[CCDNTD[i] * 4 + 0], asset->particles[CCDNTD[i] * 4 + 1], asset->particles[CCDNTD[i] * 4 + 2]);
		}
		clusterOffsets.push_back(int(clusterIndices.size()));
		clusterPositions.push_back(center / CCDNTD.size());

		center = Vec3(0, 0, 0);
		for (int i = 0; i < CTD1.size(); ++i)
		{
			clusterIndices.push_back(CTD1[i]);
			center += Vec3(asset->particles[CTD1[i] * 4 + 0], asset->particles[CTD1[i] * 4 + 1], asset->particles[CTD1[i] * 4 + 2]);
		}
		clusterPositions.push_back(center / CCDNTD.size());
		clusterOffsets.push_back(int(clusterIndices.size()));
		center = Vec3(0, 0, 0);
		for (int i = 0; i < CTD1.size(); ++i){
			clusterIndices.push_back(CTD1[i] + 50);
			center += Vec3(asset->particles[(CTD1[i] + 50) * 4 + 0], asset->particles[(CTD1[i] + 50) * 4 + 1], asset->particles[(CTD1[i] + 50) * 4 + 2]);
		}
		clusterPositions.push_back(center / CCDNTD.size());
		clusterOffsets.push_back(int(clusterIndices.size()));


		center = Vec3(0, 0, 0);
		for (int i = 0; i < LinkerA.size(); ++i){
			clusterIndices.push_back(LinkerA[i]);
			center += Vec3(asset->particles[(LinkerA[i]) * 4 + 0], 
				asset->particles[(LinkerA[i]) * 4 + 1], 
				asset->particles[(LinkerA[i]) * 4 + 2]);
		}
		clusterPositions.push_back(center / LinkerA.size());
		clusterOffsets.push_back(int(clusterIndices.size()));

		center = Vec3(0, 0, 0);
		for (int i = 0; i < LinkerB.size(); ++i){
			clusterIndices.push_back(LinkerB[i]);
			center += Vec3(asset->particles[(LinkerB[i]) * 4 + 0], 
				asset->particles[(LinkerB[i]) * 4 + 1], 
				asset->particles[(LinkerB[i]) * 4 + 2]);
		}
		clusterPositions.push_back(center / LinkerA.size());
		clusterOffsets.push_back(int(clusterIndices.size()));


		// assign shapes
		clusterCoefficients.resize(numClusters, clusterStiffness);

		asset->shapeIndices = new int[clusterIndices.size()];
		memcpy(asset->shapeIndices, &clusterIndices[0], sizeof(int)*clusterIndices.size());

		asset->shapeOffsets = new int[numClusters];
		memcpy(asset->shapeOffsets, &clusterOffsets[0], sizeof(int)*numClusters);

		asset->shapeCenters = new float[numClusters * 3];
		memcpy(asset->shapeCenters, &clusterPositions[0], sizeof(float)*numClusters * 3);

		asset->shapeCoefficients = new float[numClusters];
		memcpy(asset->shapeCoefficients, &clusterCoefficients[0], sizeof(float)*numClusters);

		asset->numShapeIndices = int(clusterIndices.size());
		asset->numShapes = numClusters;
		cout << "numCluster " << numClusters << " " << asset->numShapeIndices << endl;
	}

	virtual void overwriteSprings(){
		float radius = g_params.radius*0.75f;
		float stiffness = 1.0f;
		//body come from the mapping in the sphere 

		cellPACK::IngredientSphereTree ing_spheres = cp->mIngrSphereTree[0];
		NvFlexExtAsset* asset = cp->iBatches[0].mAsset;
		std::vector<int> springIndices;
		std::vector<float> springLengths;
		std::vector<float> springStiffness;
		int numLinks = 0;
		std::vector<int> springIndicestmp;
		std::vector<float> springLengthstmp;
		std::vector<float> springStiffnesstmp;
		// create links between particles

		int numClusters = ing_spheres.LevelCounts[1];

		std::vector<int> clusterIndices;
		std::vector<int> clusterOffsets;
		std::vector<Vec3> clusterPositions;
		std::vector<float> clusterCoefficients;
		float clusterStiffness = 1.0f;

		//assume 4 LOD atom,lvl0,lvl1,lvl2 lvl0 used in flex. lvl2 used for the softbody
		int start = ing_spheres.LevelMappingStarts[2];
		for (int i = 0; i < ing_spheres.LevelCounts[1]; i++){
			//use the mapping to get the id.
			int map_start = ing_spheres.LevelMapping[start + i].x;
			int map_count = ing_spheres.LevelMapping[start + i].y;
			std::vector<int> mapping;
			Vec3 center = Vec3(0, 0, 0);
			for (int m = map_start; m < (map_start + map_count); m++){
				mapping.push_back(m);
				clusterIndices.push_back(m);
				//cout << "cluster " << i << " ptid " << m << endl;
				center += Vec3(asset->particles[m * 4 + 0],
					asset->particles[m * 4 + 1],
					asset->particles[m * 4 + 2]);
			}
			//cout << "start " << map_start << " " << map_count << " " << ing_spheres.LevelStarts[1] << clusterIndices.size() << endl;
			//create link for all domain except 0 and 3 
			if ((i != 0)&(i != 3)){
				int numLinksTmp = CreateLinksLocal(asset->particles, mapping.size(),
					mapping,
					springIndices, springLengths, springStiffness,
					radius, stiffness);
				numLinks += numLinksTmp;
			}
			clusterOffsets.push_back(int(clusterIndices.size()));
			clusterPositions.push_back(center / map_count);// ing_spheres.LevelPoints[ing_spheres.LevelStarts[1]]);
		}

		//dimer connect the two central domain 0 and 3
		int map_start1 = ing_spheres.LevelMapping[start + 0].x;
		int map_count1 = ing_spheres.LevelMapping[start + 0].y;
		int map_start2 = ing_spheres.LevelMapping[start + 3].x;
		int map_count2 = ing_spheres.LevelMapping[start + 3].y;
		std::vector<int> mapping;
		Vec3 center = Vec3(0, 0, 0);
		for (int m = map_start1; m < (map_start1 + map_count1); m++){
			mapping.push_back(m);
			//clusterIndices.push_back(m);
			center += Vec3(asset->particles[m * 4 + 0],
				asset->particles[m * 4 + 1],
				asset->particles[m * 4 + 2]);
		}
		for (int m = map_start2; m < (map_start2 + map_count2); m++){
			mapping.push_back(m);
			//clusterIndices.push_back(m);
			center += Vec3(asset->particles[m * 4 + 0],
				asset->particles[m * 4 + 1],
				asset->particles[m * 4 + 2]);
		}
		int numLinksTmp = CreateLinksLocal(asset->particles, mapping.size(), mapping,
			springIndices, springLengths, springStiffness,
			radius + 0.01f, stiffness);
		numLinks += numLinksTmp;

		//Linker between body, 
		for (int i = 0; i < ing_spheres.LevelCounts[1] - 1; i++){
			//link to next one
			//skip the middle one because homodimer
			if (i == 2) continue;
			int last = ing_spheres.LevelMapping[start + i].x + ing_spheres.LevelMapping[start + i].y - 1;
			int next = ing_spheres.LevelMapping[start + i + 1].x;
			springIndices.push_back(last);
			springIndices.push_back(next);
			springLengths.push_back(Length(Vec3(&asset->particles[last * 4]) - Vec3(&asset->particles[next * 4])));
			springStiffness.push_back(stiffness);
			numLinks++;
		}

		// assign links
		if (numLinks)
		{
			asset->springIndices = new int[numLinks * 2];
			memcpy(asset->springIndices, &springIndices[0], sizeof(int)*springIndices.size());

			asset->springCoefficients = new float[numLinks];
			memcpy(asset->springCoefficients, &springStiffness[0], sizeof(float)*numLinks);

			asset->springRestLengths = new float[numLinks];
			memcpy(asset->springRestLengths, &springLengths[0], sizeof(float)*numLinks);

			asset->numSprings = numLinks;
		}

		// assign shapes
		clusterCoefficients.resize(numClusters, clusterStiffness);

		asset->shapeIndices = new int[clusterIndices.size()];
		memcpy(asset->shapeIndices, &clusterIndices[0], sizeof(int)*clusterIndices.size());

		asset->shapeOffsets = new int[numClusters];
		memcpy(asset->shapeOffsets, &clusterOffsets[0], sizeof(int)*numClusters);

		asset->shapeCenters = new float[numClusters * 3];
		memcpy(asset->shapeCenters, &clusterPositions[0], sizeof(float)*numClusters * 3);

		asset->shapeCoefficients = new float[numClusters];
		memcpy(asset->shapeCoefficients, &clusterCoefficients[0], sizeof(float)*numClusters);

		asset->numShapeIndices = int(clusterIndices.size());
		asset->numShapes = numClusters;
		cout << "numCluster " << numClusters << endl;
	}

	virtual void Initialize()
	{
		bool ignore_comp = true;
		/* integrase initialization */
		std::string recipe = "recipes/HIV_IN_XP.json";
		std::string results = "recipes/INT_50_random.json";
		std::string wrkDir = "../../data/cellpack/";
		ccd_ctd_1 = { 24, 25, 79 };
		ccd_ctd_2 = { 74, 75, 29 };
		ccd_ctd.push_back(ccd_ctd_1);
		ccd_ctd.push_back(ccd_ctd_2);

		ctd_ccd_1 = { 40, 41, 42, 43, 49 };
		ctd_ccd_2 = { 90, 91, 92, 93, 99 };
		ctd_ccd.push_back(ctd_ccd_1);
		ctd_ccd.push_back(ctd_ccd_2);


		main_scale = 1.0f / 100.0f;

		float beads_radius = 6.5f*main_scale;
		g_params.radius = beads_radius * 2.0f;
		threshold_binding = g_params.radius*10.0f;
		
		threshold_binding = g_params.radius*2.0f;
		//setup experiments
		exp_params[0] = { 1.0f, 1.0f, g_params.radius*2.0f };
		exp_params[1] = { 100.0f, 0.0f, g_params.radius*100.0f };//no bias toward center but higher distance and random force
		exp_params[2] = { 100.0f, 0.0f, g_params.radius*10000.0f };//higher distance or infinite distance ?
		nexp = 3;
		nrun = 3;
		//g_mesh

		cp = new cellPACK();
		cp->use_rb = false;

		cp->main_radius = beads_radius;//((10.0f/100.0f)*(10.0f/100.0f));
		cp->main_scale = main_scale;
		cp->mainpath = wrkDir;// HIV_IN_XP.json"
		cp->datapath = wrkDir;// HIV_IN_XP.json"
		cp->loadRecipe(wrkDir + recipe, ignore_comp);
		//cp->loadRecipe((cp->mainpath + "HIV_IN_XP.1.0.json").c_str(), ignore_comp);
		overwriteSprings();

		cp->loadResults(wrkDir + results);
		//if (!ignore_comp)cp->buildMembrane(1.0f, 1.0f, flexMakePhase(99999, 0));

		setupBoundsAndLinker();

		int group = cp->iGroupCounter;
		g_numExtraParticles = 1024 * 1024;
		g_params.numPlanes = 0;
		g_params.gravity[1] = 0.0f;//-9.f;
		g_params.damping = 3.0f;// 3.0f;
		g_params.radius = beads_radius * 2;
		//g_params.mSolidRestDistance = g_params.mRadius*2.0f;
		g_windStrength = 0.0f;
		g_windFrequency = 0.0f;
		//        g_params.mDynamicFriction = 0.00f;
		//        g_params.mFluid = true;
		//        g_params.mViscosity = 0.0f;
		// g_params.fluid = false;
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

		g_params.gravity[1] = 0.0f;

		g_lightDistance *= 2.0f;

		// draw options		
		g_drawEllipsoids = true;
		g_params.maxSpeed = 100.0f;

		g_wireframe = true;
		g_pointScale = 1.0f;
		//        g_blur = 2.0f;
		g_pause = true;
		g_warmup = false;
		g_drawMesh = true;
		g_drawRopes = false;
		g_drawPoints = true;
		g_params.dynamicFriction = 0.4f;

		g_params.dissipation = 0.0f;
		g_params.numIterations = 4;
		g_params.viscosity = 0.0f;
		g_params.drag = 0.0f;
		g_params.lift = 0.0f;

		g_params.particleCollisionMargin = g_params.radius*0.05f;
		g_params.drag = 0.0f;
		g_params.collisionDistance = 0.01f;

		// better convergence with global relaxation factor
		g_params.relaxationMode = eNvFlexRelaxationGlobal;
		g_params.relaxationFactor = 0.25f;

		g_windStrength = 0.0f;

		g_numSubsteps = 5;
		g_params.numIterations = 5;

		//mSplitThreshold.resize(mCloths.size(), 45.0f);

		// draw options		
		g_drawSprings = 1;
		g_drawCloth = false;
	}

	virtual void setup_Exp(){
		//exp id define parameters 
		//current_Exp

		dojitter_strength = exp_params[current_Exp][0];
		dojitter_biased_strength = exp_params[current_Exp][1];
		threshold_binding = exp_params[current_Exp][2];
		//reset everything to initial position?
		record_free.clear();
		record_size.clear();
		//remove all spring
		/*
		g_buffers->springIndices.erase(g_buffers->springIndices.begin() + strings_indices[0] * 2, g_buffers->springIndices.begin() + strings_indices[strings_indices.size() - 1] * 2 + 2);
		g_buffers->springLengths.erase(g_buffers->springLengths.begin() + strings_indices[0], g_buffers->springLengths.begin() + strings_indices[strings_indices.size() - 1] + 1);
		g_buffers->springStiffness.erase(g_buffers->springStiffness.begin() + strings_indices[0], g_buffers->springStiffness.begin() + strings_indices[strings_indices.size() - 1] + 1);
		*/
		g_buffers->springIndices.resize(0);
		g_buffers->springLengths.resize(0);
		g_buffers->springStiffness.resize(0);
		NvFlexSetSprings(g_solver, g_buffers->springIndices.buffer, g_buffers->springLengths.buffer, g_buffers->springStiffness.buffer, g_buffers->springLengths.size());

		strings_indices.clear();
		mIntegrase.clear();
		free_mIntegrase_CCD.clear();
		free_mIntegrase_CTD.clear();
		setupBoundsAndLinker();
		//replace object to original position
		//cp->loadResults("INT_50_random.json", true);
		//nrun ?
		current_run++;
		if (current_run >= nrun){
			current_run = 0;
			current_Exp++;
			if (current_Exp >= nexp){
				//stop
				current_Exp = -1;
				g_pause = false;
			}
		}
		g_pause = false;
		if (current_Exp == -1) {
			g_pause = true;
			return;
		}
	}

	void showNetowrk(){
		//create particle for each instance
		//connect according binding.
		//one instance pos = one particle (random ?)
		//one connection ccd->ctd
	}

	void checkDistanceBinding(){
		//find the closest pair ctd-ccd free
		//use a list ?
		int mini_ccd_ind1 = -1;
		int mini_ctd_ind2 = -1;
		int mini_instance_id1 = -1;
		int mini_instance_id2 = -1;
		float mini_distance = 999999.0f;

		//should I use free_mIntegrase_CCD?
		for (int i = 0; i < int(cp->mInstances.size()); ++i)
		{
			NvFlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			int poffseti = cp->mInstances[i].mParticleOffset;

			mIntegrase[i].mini_distance[0] = 9999.0f;
			mIntegrase[i].mini_distance[1] = 9999.0f;

			for (int j = 0; j < int(cp->mInstances.size()); j++) {
				if (i == j) continue;
				int poffsetj = cp->mInstances[j].mParticleOffset;
				for (int o = 0; o < 2; o++){ //ccd
					for (int p = 0; p < 2; p++){ //ctd
						float D = getMinDistances(i, j, o, p);
						if ((mIntegrase[i].ccd_occupied[o] != -1) || (mIntegrase[j].ctd_occupied[p] != -1)){ continue; }
						//we can constraint only one interaction per pair of dimer or not
						//if ((mIntegrase[i].ccd_occupied[1 - o] == j) || (mIntegrase[j].ctd_occupied[1 - p] == i)){ continue; } //crab
						//if ((mIntegrase[i].ctd_occupied[1 - o] == j) || (mIntegrase[j].ccd_occupied[1 - p] == i)){ continue; } //head-tail-cross
						//if ((mIntegrase[i].ctd_occupied[o] == j) || (mIntegrase[j].ccd_occupied[p] == i)){ continue; }		 //head-tail
						//for (int n = 0; n < ccd_ctd[o].size(); n++) {//beads ccd
						//	for (int m = 0; m < ctd_ccd[p].size(); m++) {//beads ctd
						//		float D = Length(Vec3(g_positions[poffseti + ccd_ctd[o][n]]) - Vec3(g_positions[poffsetj + ctd_ccd[p][m]]));
						if (D < mini_distance){
							if (D < threshold_binding){
								//check if free
								//if ((mIntegrase[i].ccd_occupied[o]==-1) && (mIntegrase[j].ctd_occupied[p]==-1)){
								//	if ((mIntegrase[i].ccd_occupied[1 - o] != j) && (mIntegrase[j].ctd_occupied[1 - p] != i)){
								mini_distance = D;
								mini_ccd_ind1 = o;
								mini_ctd_ind2 = p;
								mini_instance_id1 = i;
								mini_instance_id2 = j;
								//	}
								//}
							}
						}
						//	}
						//}
					}
				}
			}
			/*if ((mIntegrase[i].ccd_occupied[0] == -1) || (mIntegrase[i].ccd_occupied[1] == -1)){
			cout << endl;
			cout << mIntegrase[i].mini_distance[0] << " " << mIntegrase[i].mini_distance[1] << " " << mIntegrase[i].mini_indices[0] << " " << mIntegrase[i].mini_indices[1] << " " << mini_distance << endl;
			cout << mIntegrase[i].ccd_occupied[0] << " " << mIntegrase[i].ccd_occupied[1] << " " << mIntegrase[i].ctd_occupied[0] << " " << mIntegrase[i].ctd_occupied[1] << endl;
			cout << mIntegrase[i].ccd_ctdind[0] << " " << mIntegrase[i].ccd_ctdind[1] << " " << mIntegrase[i].ctd_ccdind[0] << " " << mIntegrase[i].ctd_ccdind[1] << endl;
			}*/
		}
		if (mini_instance_id1 != -1){
			boundTwoIds(mini_instance_id1, mini_instance_id2, mini_ccd_ind1, mini_ctd_ind2);
			//free_mIntegrase_CCD.erase(free_mIntegrase_CCD.begin() + mini_instance_id1);
			//free_mIntegrase_CTD.erase(free_mIntegrase_CTD.begin() + mini_instance_id2);
		}
		else {
			//cout << " mini distance found is " << mini_distance << endl;
			//check the mini_distances..same as no distance threshold
			return;
			int res1 = 0;
			int res2 = 0;
			for (int i = 0; i < int(cp->mInstances.size()) - 1; ++i)
			{

				if ((mIntegrase[i].ccd_occupied[0] == -1)){
					//cout << mIntegrase[i].mini_distance[0] << " " << mIntegrase[i].mini_distance[1] << " " << mIntegrase[i].mini_indices[0] << " " << mIntegrase[i].mini_indices[1] << " " << mini_distance << endl;
					res1 = boundTwoIds(i, mIntegrase[i].mini_indices[0], 0, 0);
					res2 = boundTwoIds(i, mIntegrase[i].mini_indices[1], 0, 1);
					if (res1 != -1 || res2 != -1)
						break;
				}
				if (mIntegrase[i].ccd_occupied[1] == -1) {
					//cout << mIntegrase[i].mini_distance[0] << " " << mIntegrase[i].mini_distance[1] << " " << mIntegrase[i].mini_indices[0] << " " << mIntegrase[i].mini_indices[1] << " " << mini_distance << endl;
					res1 = boundTwoIds(i, mIntegrase[i].mini_indices[0], 1, 0);
					res2 = boundTwoIds(i, mIntegrase[i].mini_indices[1], 1, 1);
					if (res1 != -1 || res2 != -1)
						break;
				}
			}
			//if (res1 == -1 && res2 == -1)
			//	cout << "cant find another binding" << endl;
		}
	}

	void checkRadiusAggregate(){
		std::vector<Vec3> all_pos;
		Vec3 agg_center = Vec3(0, 0, 0);
		for (int i = 0; i < int(cp->mInstances.size()); ++i){
			int poffseti = cp->mInstances[i].mParticleOffset;
			NvFlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			Vec3 center = Vec3(0, 0, 0);
			for (int j = 0; j < asset->numParticles; j++){
				center += Vec3(g_buffers->positions[poffseti + j]);
			}
			center /= asset->numParticles;
			all_pos.push_back(center);
			agg_center += center;
		}
		agg_center /= all_pos.size();
		float maxiD = 0.0f;
		for (int i = 0; i < all_pos.size(); i++){
			float D = Length(all_pos[i] - agg_center);
			if (D > maxiD)
				maxiD = D;
		}
		//cout << " max distance is " << maxiD;
		record_size.push_back(maxiD);
	}

	void writeReport(int exp, int run, int free, bool append = false){
		//write information about this model
		//radius, nb free binding site ccd and ctd
		//networ of interaction 1 node can be attach to 4 other object
		//report integrase struct
		/*	struct Integrase
		{
		int nb_site_occupied;
		int ccd_occupied[2];
		int ccd_ctdind[2];
		int ctd_occupied[2];
		int ctd_ccdind[2];
		int spring_offset[2];
		float mini_distance[2];
		int mini_indices[2];
		};*/
		int  unsatisfied_count = checkDistance(false);
		ofstream of;
		std::string fname = "../../data/pack_data_soft";
		int nintegrase = mIntegrase.size();
		if (append) {
			fname = fname + ".txt";
			of.open(fname.c_str(), ios::out | ios::app);
			of << "# " << nintegrase << " " << exp << " " << run << " " << free << " " << unsatisfied_count << endl;
		}
		else {
			fname = fname + "_" + std::to_string(exp) + "_" + std::to_string(run) + "_" + std::to_string(free) + ".txt";
			of.open(fname.c_str(), ios::out);
		}

		//node / instance Id neiboorccd 1 neiborccd 2 neiboorctd 1 neiborctd 2
		for (int i = 0; i < nintegrase; i++){
			of << i << " " << mIntegrase[i].ccd_occupied[0] << " " << mIntegrase[i].ccd_occupied[1];// << " " << mIntegrase[i].ctd_occupied[0] << " " << mIntegrase[i].ctd_occupied[1];
			of << " " << mIntegrase[i].ccd_ctdind[0] << " " << mIntegrase[i].ccd_ctdind[1];// << " " << mIntegrase[i].ctd_ccdind[0] << " " << mIntegrase[i].ctd_ccdind[1];
			of << " " << mIntegrase[i].mini_distance[0] << " " << mIntegrase[i].mini_distance[1] << endl;// << " " << mIntegrase[i].ctd_ccdind[0] << " " << mIntegrase[i].ctd_ccdind[1] << endl;
		}
		of.close();
	}

	void Update()
	{
		//NvFlexTimers atimers;
		//float ** p;
		//float ** v;
		//int ** ph;
		//float ** no;
		//int a;
		//distance 71-3
		//distance center->farest dimer
		//reportBinding();

		//if (cp->comp_shape.size())
		//	flexExtSetRigidTransformation(cp->fcontainer, cp->mInstances.size());
		if ((g_frame % 5) == 0)
		{
			if (dojitter) cp->jitter(!dojitter_steared,false);
			//cout << "free " << free_mIntegrase_CCD.size() << " " << free_mIntegrase_CTD.size() << endl;
			//if (free_mIntegrase_CCD.size() != 0) addOneRandomBinding();
		}
		checkDistanceBinding();
		
		//WeightSpringStifness(strings_indices);
		//flexSetSprings(g_flex, &g_springIndices[0], &g_springLengths[0], &g_springStiffness[0], g_springLengths.size(), eFlexMemoryHost);

		//save packing leftHand
		if ((g_frame % 1) == 0)
		{
			//if (server)
			//	sendToClient();
			if (write_output){
			    std::string pfix = "";// "_" + std::to_string(current_Exp) + "_" + std::to_string(current_run);
				cp->writeToBinary(pfix);
			}
		}
		checkRadiusAggregate();
		reportBinding(false);
		//check stoping criteria
		if (dosimulation){
			if (record_free.size() == 100 && current_Exp != -1){
				//test..compare average and last value
				bool stop = false;
				int free = 0;
				float radius_agg = 0;
				for (int i = 0; i < 100; i++){
					free += record_free[i];
					radius_agg += record_size[i];
				}
				free /= 100;
				radius_agg /= 100;
				int dtfree = free - record_free[99];
				float dtrad = radius_agg - record_size[99];
				//remove first entry
				if (dtfree == 0 && dtrad < 0.05f)  count_stop++;
				record_free.erase(record_free.begin());
				record_size.erase(record_size.begin());
				if (count_stop > 100)
				{
					count_stop = 0;
					g_pause = true;
					//writeSoftTransform(current_Exp, current_run, free, true);	//#the beads
					//writeRigidTransform(current_Exp, current_run, true);		//the shape rigid bodu
					exportPDBRigidTransform(current_Exp, current_run);			//shape 1 as PDB dummy atom
					//cout << "1 free " << dtfree << " " << dtrad << " " << current_Exp << " " << current_run << endl;
					writeReport(current_Exp, current_run, free, true);			//connectivity
					cout << "2 free " << dtfree << " " << dtrad << " " << current_Exp << " " << current_run << endl;
					//init next experiment
					setup_Exp();
					cout << "free " << dtfree << " " << dtrad << " " << current_Exp << " " << current_run << endl;
				}

			}
		}
		if (current_Exp == -1) current_Exp = 0;
		//if ((g_frame % 5) && (grow_fiber))
		//	cp->growFiber((int)Nsub);
	}

	void exportPDBRigidTransform(int exp, int run)
	{
		NvFlexVector<Quat> quat(g_flexLib);
		NvFlexVector<Vec3> pos(g_flexLib);
		quat.map();
		pos.map();
		int totalNbBody = 0;
		int nInst = cp->mInstances.size();
		if (cp->use_rb){
			pos.resize(cp->mInstances.size());
			quat.resize(cp->mInstances.size());
			//pos = new Vec3[cp->mInstances.size()];
			//quat = new Quat[cp->mInstances.size()];
			totalNbBody = nInst;
		}
		else {
			for (int i = 0; i < cp->mInstances.size(); i++){
				NvFlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
				totalNbBody += asset->numShapes;
			}
			//pos = new Vec3[totalNbBody];
			//quat = new Quat[totalNbBody];
			pos.resize(totalNbBody);
			quat.resize(totalNbBody);
		}
		NvFlexGetRigids(g_solver, NULL, NULL, NULL, NULL, NULL, NULL, NULL, quat.buffer, pos.buffer);
		std::string name = "../../data/pack_result" + std::to_string(exp) + "_" + std::to_string(run);
		name = name + ".pdb";
		FILE *of;

		fopen_s(&of, name.c_str(), "w");
		//write position
		for (int i = 0; i<nInst; i++){
			//cout << "  ? " << i << endl;
			//output << -pos[i].x*(1.0f/main_scale) << " " <<pos[i].y*(1.0f/main_scale)<< " " <<pos[i].z*(1.0f/main_scale)<< " " << cp->mInstances[i].mMeshIndex << "\n";
			//int occup = getOccupationCCD(i);
			float bf = (mNodes[i].mini_distance[0] + mNodes[i].mini_distance[1]) / 2.0f;
			if (bf > 500.0f) bf = 0.0f;
			float ind = (float)cp->mInstances[i].mMeshIndex;//instance 
			NvFlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			int nbody = asset->numShapes;
			int bodyindice = i*nbody;
			float p[4] = { pos[bodyindice].x*(1.0f / main_scale)*0.05f,
				pos[bodyindice].y*(1.0f / main_scale)*0.05f,
				pos[bodyindice].z*(1.0f / main_scale)*0.05f, ind };
			fprintf(of, "%6s%5d%3s %1s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
				"ATOM  ", i, "C", " ", "ALA", "A", i, " ", p[0], p[1], p[2], (float)0, bf, "C", "");
		}
		//connectivity
		//8 -1 -1 -1 -1 9999 9999
		//7 49 9 1 0 0.110658 0.112724
		int nNodes = mNodes.size();
		for (int i = 0; i < nNodes; i++){
			if (mNodes[i].nb_site_occupied>0){
				fprintf(of, "CONECT%5d", i);//node i is connected to (can be connected to 3 other node)
				for (int b = 0; b < binding_in.size(); b++){
					if (mNodes[i].binding_in_occupied[b] != -1)
						fprintf(of, "%5d", mNodes[i].binding_in_occupied[b]);
					fprintf(of, "\n");
				}
			}
		}
		fclose(of);
	}


	virtual void Sync()
	{
		NvFlexSetRigids(g_solver,	g_buffers->rigidOffsets.buffer, 
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
};
