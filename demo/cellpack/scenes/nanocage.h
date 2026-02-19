class NanoCage : public Scene
{
public:
	cellPACK * cp;
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

	float main_scale = 1.0f / 100.0f;
	std::vector<std::vector<int> > binding_in;
	std::vector<std::vector<int> > binding_out;
	std::vector<int> strings_indices;
	int group;
	int stop_criterion_time = 100; // how many last frame to store and check
	int count_stop = 0;
	int nexp = 0;
	int current_Exp = 0; // how many experiment?
	int nrun = 5;
	int current_run = 0;
	std::map <int, std::vector<float>> exp_params;//carry other the experiments parameters
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

	std::vector<Integrase> mIntegrase;
	std::vector<int> free_mIntegrase_CCD;
	std::vector<int> free_mIntegrase_CTD;
	std::string wrkDir;
	bool bindTwo = false;
	bool next_step = false;
	int current_step;
	int current_node_id=0;
	int current_edge_id = 0;
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
		int particle_index;
		int node_index; 
	};
	std::vector<Node> mNodes;

	Graph G;
	GraphAttributes GA;
	List<edge> edges;//one edge is two nodes, source and target
	List<node> nodes;
	//NodeArray<bool> visited(G, false);
	//EdgeArray<bool> edge_visited(H, false);
	int poffset = 0;
	NanoCage(const char* name) : Scene(name) {}

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
			//NvFlexGetParticles(g_flex, g_buffers->positions.buffer, g_buffers->positions.size());
			bindTwo = false;
			next_step = true;
		}

		if (key == 'b')
		{
			std::cout << "key b" << endl;
			//reportBinding();
			//checkDistance(false);
			//checkRadiusAggregate();
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
		//randomly pick one acceptor//one donor
		int i = (int) (Rand() % cp->mInstances.size());
		int j = (int) (Rand() % cp->mInstances.size());
		int id1 = Rand() % 2;
		int id2 = Rand() % 2;
		boundTwoIds(i, j, id1, id2);
	}

	virtual float getMinDistances(int i, int j, int ccd, int ctd){
		float miniD = 999999.9f;
		int poffseti = cp->mInstances[i].mParticleOffset;
		int poffsetj = cp->mInstances[j].mParticleOffset;
		for (int n = 0; n < binding_in[ccd].size(); n++) {
			for (int m = 0; m < binding_out[ctd].size(); m++) {
				float D = Length(Vec3(g_buffers->positions[poffseti + binding_in[ccd][n]]) - Vec3(g_buffers->positions[poffsetj + binding_out[ctd][m]]));
				if (D < miniD) {
					miniD = D;
				}
			}
		}
		return miniD;
	}

	virtual int checkOccupation(int i, int j, int b){
		for (int n = 0; n < binding_in.size(); n++) {
			if (n == b) continue;
			if (mNodes[i].binding_in_occupied[n] == j) return n;
			if (mNodes[j].binding_in_occupied[n] == i) return n;
			if (mNodes[i].binding_out_occupied[n] == j) return n;
			if (mNodes[j].binding_out_occupied[n] == i) return n;
		}
		return -1;
	}

	virtual int boundTwoIds(int i, int j, int id1, int id2){
		int poffseti = cp->mInstances[i].mParticleOffset;
		int poffsetj = cp->mInstances[j].mParticleOffset;
		int sp_offset = strings_indices.size();
		if (mNodes[i].binding_in_occupied[id1] != -1) return -1;
		if (mNodes[j].binding_out_occupied[id2] != -1) return -1;
		cout << "bind " << i << " " << id1 << " to " << j << " " << id2 << endl;
		cout << "size " << binding_in[id1].size() << endl;//5 
		cout << "size " << binding_out[id2].size() << endl;//0 4 ?
		for (int n = 0; n < binding_in[id1].size(); n++) {
			for (int m = 0; m < binding_out[id2].size(); m++) {
				float stiff = 1.0f;
				float length = g_params.radius; // cp->mIngrSphereTree[0].DistancesMatrix[n*binding_in[id1].size() + m] * cp->main_scale;
				//cout << "bound id1 " << id1 << " " << n << " " << poffseti << " " <<binding_in[id1][n] << endl;
				//cout << "to  id2 " << id2 << " " << m << " " << poffsetj << " " << binding_out[id2][m] << endl;
				//cout << "length " << length << endl;
				strings_indices.push_back(g_buffers->springLengths.size());
				int inda = poffseti + binding_in[id1][n];
				int indb = poffsetj + binding_out[id2][m];
				//cout << g_buffers->positions.size() << endl;
				//cout << inda << " x " << g_buffers->positions[inda].x << endl;
				//cout << indb << " x " << g_buffers->positions[indb].x << endl;
				CreateSpringInter(poffseti + binding_in[id1][n], poffsetj + binding_out[id2][m], stiff, 0.0f, length); //float give = 0.0f, float length = 0.0f
			}
		}
		mNodes[i].binding_in_occupied[id1] = j;
		mNodes[j].binding_out_occupied[id2] = i;
		mNodes[j].binding_in_occupied[id2] = i;
		mNodes[i].binding_in_ind[id1] = id2;
		mNodes[j].binding_in_ind[id2] = id1;
		mNodes[i].spring_offset[id1] = sp_offset;
		mNodes[i].nb_site_occupied++;
		//bound the node 
		//flexSetSprings(g_flex, &g_springIndices[0], &g_buffers->springLengths[0], &g_springStiffness[0], g_buffers->springLengths.size(), eFlexMemoryHost);
		cout << "bound done " << i << " " << j << " " << id1 << " " << id2 << " " << mNodes[i].binding_in_occupied[id1] << " " << mNodes[j].binding_out_occupied[id2] << endl;
		return 1;
	}

	virtual void setupGraphx(){
		//nanocage graph is 12 node
		//each node has 3 possible binding 0,1,2. indexing clockwise
		//0 [0 - 1]
		//0 [1 - 2]
		//0 [2 - 3]
		//1 [0 - 2]
		//1 [1 - 0]
		//1 [2 - 4]
		//2 [0 - 0]
		//2 [1 - 1]
		//2 [2 - 5]
		//3 [0 - 6]
		//3 [1 - 7]
		//3 [2 - 0]
		//4 [0 - 8]
		//4 [1 - 9]
		//4 [2 - 1]
		//5 [0 - 10]
		//5 [1 - 11]
		//5 [2 - 2]
		//6 [0 - 7]
		//6 [1 - 3]
		//6 [2 - 11]	
		//7 [0 - 7]
		//7 [1 - 3]
		//7 [2 - 8]
		//8 [0 - 9]
		//8 [1 - 4]
		//8 [2 - 7]
		//9 [0 - 8]
		//9 [1 - 4]
		//9 [2 - 10]
		//10 [0 - 11]
		//10 [1 - 5]
		//10 [2 - 9]
		//11 [0 - 10]
		//11 [1 - 5]
		//11 [2 - 7]
	}

	virtual void OneNodeInstance(int index){
		//IngredientInstance& inst = cp->mInstances[i];
		Node node = Node();
		for (int i = 0; i < binding_in.size(); i++){
			node.binding_in_occupied.push_back(-1);
			node.binding_in_ind.push_back(-1);
			node.binding_out_occupied.push_back(-1);
			node.binding_out_ind.push_back(-1);
			node.spring_offset.push_back(-1);
			node.mini_distance.push_back(99999.9f);
			node.mini_indices.push_back(-1);
		}
		node.nb_site_occupied = 0;
		node.node_index = index;
		/*int phase = NvFlexMakePhase(9999999, 0);
		int particleOffset = g_buffers->positions.size();
		g_buffers->positions[particleOffset] = Vec4(RandomUnitVector(), 1.0f);
		g_buffers->velocities[particleOffset] = Vec3(0.0f);
		g_buffers->phases[particleOffset] = phase;
		node.particle_index = particleOffset;
		*/
		mNodes[index] = node;
	}

	virtual void setupBoundsAndLinker(){
		//go over all instances, bind to every other instances
		//indicesBounds
		for (int i = 0; i < int(cp->mInstances.size()); ++i)
		{
			//IngredientInstance& inst = cp->mInstances[i];
			Node node = Node();
			for (int i = 0; i < binding_in.size(); i++){
				node.binding_in_occupied.push_back(-1);
				node.binding_in_ind.push_back(-1);
				node.binding_out_occupied.push_back(-1);
				node.binding_out_ind.push_back(-1);
				node.spring_offset.push_back(-1);
				node.mini_distance.push_back(99999.9f);
				node.mini_indices.push_back(-1);
			}
			node.nb_site_occupied = 0;
			mNodes.push_back(node);
			//free_mIntegrase_CCD.push_back(i);
			//free_mIntegrase_CTD.push_back(i);
		}
		//cout << " nb integrase " << mIntegrase.size() << endl;
	}

	virtual void setupBinding(){
		cellPACK::IngredientSphereTree sph = cp->mIngrSphereTree[0];
		for (int i = 0; i < sph.BindingStarts.size(); i += 2){
			int start = sph.BindingStarts[i];
			int count = sph.BindingStarts[i + 1];
			std:vector<int> lb;
			cout << i << " start " << start << " count " << count << endl;
			for (int j = 0; j < count; j++){
				int id = sph.BindingSites[start + j];
				lb.push_back(id);
				cout << "add binding id " << id << endl;
			}
			binding_in.push_back(lb);
		}
		binding_out.push_back(binding_in[0]);
		binding_out.push_back(binding_in[1]);
		binding_out.push_back(binding_in[2]);
		binding_out.push_back(binding_in[2]);
	}

	virtual void overwriteSprings(){
		float radius = g_params.radius*0.75f;
		float stiffness = 1.0f;
		//body come from the mapping in the sphere tree file

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

	virtual void setupGraph(){
		string filename = wrkDir + "nanocage.gml";
		GA = GraphAttributes(G, GA.nodeLabel | GA.nodeId | GA.edgeLabel);
		GraphIO::readGML(GA, G, filename);
		G.allEdges(edges);
		G.allNodes(nodes);
		cp->mInstances.clear();
		cp->mInstances.resize(G.numberOfNodes());
		mNodes.resize(G.numberOfNodes());
	}

	virtual void addNodesParticles(int particleOffset){
		for (int i = 0; i < cp->iBatches[0].nInstances; i++){
			//also create the beads for it ?
			int phase = NvFlexMakePhase(cp->iGroupCounter++, eNvFlexPhaseSelfCollide | eNvFlexPhaseSelfCollideFilter);
			g_buffers->positions[particleOffset] = Vec4(RandomUnitVector(), 1.0f);
			g_buffers->velocities[particleOffset] = Vec3(0.0f);
			g_buffers->phases[particleOffset] = phase;
			mNodes[i].particle_index = particleOffset;
			particleOffset++;
		}
		g_buffers->activeIndices.resize(particleOffset);
		for (int i = 0; i < particleOffset; ++i)
			g_buffers->activeIndices[i] = i;
		
		//CreateSpringInter(mNodes[n1->index()].particle_index, mNodes[n2->index()].particle_index, 1.0, 0.0f, 2.0*g_params.radius); //float give = 0.0f, float length = 0.0f
		for (auto e : edges){
			node n1 = e->source();
			node n2 = e->target();
			int btype = atoi(GA.label(e).c_str());
			float length = g_params.radius;
			if (btype == 0) {
				length = g_params.radius;
			}
			CreateSpringInter(mNodes[n1->index()].particle_index, mNodes[n2->index()].particle_index, 1.0, 0.0f,length); //float give = 0.0f, float length = 0.0f
		}
	}

	virtual void do_nextStep(){
		//get next node
		if (current_edge_id >= edges.size())
		{
			addNodesParticles(poffset);
			return;
		}
		edge e = *edges.get(current_edge_id);
		node n1 = e->source();
		node n2 = e->target();
		
		int btype = atoi(GA.label(e).c_str());
		//type is 0 or 2 for nanocage
		cout << " connect " << n1->index() << " " << GA.label(n1) << endl;
		cout << " to " << n2->index() << " " << GA.label(n2) << endl;
		cout << " type " << btype << endl;
		cout << " ninstance " << cp->iBatches[0].nInstances << endl;
		if (cp->iBatches[0].nInstances <= n1->index()){
			
			//poffset = cp->createInstanceIngredientAt(0, RandomUnitVector(), Quat(), n1->index(), poffset);
			poffset = cp->oneInstanceAt(0, n1->index(), RandomUnitVector(), Quat(), poffset);
			OneNodeInstance(n1->index());
		}
		if (cp->iBatches[0].nInstances <= n2->index()){
			//poffset = cp->createInstanceIngredientAt(0, RandomUnitVector(), Quat(), n2->index(), poffset);
			poffset = cp->oneInstanceAt(0, n2->index(), RandomUnitVector(), Quat(), poffset);
			OneNodeInstance(n2->index());
		}
		boundTwoIds(n1->index(), n2->index(), btype, btype+1);
		current_edge_id++;
		cp->compactInstances();
		//addNodesParticles(poffset);
		//CreateSpringInter(mNodes[n1->index()].particle_index, mNodes[n2->index()].particle_index, 1.0, 0.0f, 2.0*g_params.radius); //float give = 0.0f, float length = 0.0f
	}

	virtual void Initialize()
	{
		current_edge_id = 0;
		poffset = 0;
		g_buffers->rigidTranslations.resize(0);
		g_buffers->rigidRotations.resize(0);
		g_buffers->rigidCoefficients.resize(0);
		g_buffers->rigidIndices.resize(0);
		g_buffers->rigidLocalPositions.resize(0);
		g_buffers->rigidOffsets.resize(0);
		mNodes.clear();
		bool ignore_comp = true;
		/* validation experiement setup */
		std::string recipe = "NanoCage.json";
		std::string results = "nanoCage_results.json";
		wrkDir = "../../../ProteinNanocage/";

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
		cp->use_rb = true;

		cp->main_radius = beads_radius;//((10.0f/100.0f)*(10.0f/100.0f));
		cp->main_scale = main_scale;
		cp->mainpath = wrkDir;// HIV_IN_XP.json"
		cp->datapath = wrkDir;// HIV_IN_XP.json"
		cp->loadRecipe(wrkDir + recipe, ignore_comp);

		cp->output_bin.open("../../data/pack_result.bin", ios::out | ios::binary);
		cp->output_bin.close();
		
		setupGraph();

		setupBinding();

		//overwriteSprings();

		//cp->loadResults(wrkDir + results);
		
		//setupBoundsAndLinker();
		//first step 
		//do_nextStep();
		//cp->compactInstances();
		next_step = true;

		group = cp->iGroupCounter;
		g_numExtraParticles = 1024 * 1024;
		g_params.numPlanes = 0;
		g_params.gravity[1] = 0.0f;//-9.f;
		g_params.damping = 10.0f;// 3.0f;
		g_params.radius = beads_radius * 2;
		//g_params.mSolidRestDistance = g_params.mRadius*2.0f;
		g_windStrength = 0.0f;
		g_windFrequency = 0.0f;
		//        g_params.mDynamicFriction = 0.00f;
		//        g_params.mFluid = true;
		//        g_params.mViscosity = 0.0f;
		g_params.fluid = false;
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
		g_params.dynamicFriction = 1.4f;

		g_params.dissipation = 0.0f;
		g_params.viscosity = 0.0f;
		g_params.drag = 1.0f;
		g_params.lift =1.0f;

		g_params.particleCollisionMargin = g_params.radius*0.05f;

		// better convergence with global relaxation factor
		g_params.relaxationMode = eNvFlexRelaxationGlobal;
		g_params.relaxationFactor = 0.25f;

		g_windStrength = 0.0f;

		g_numSubsteps = 5;
		g_params.numIterations = 5;

		// draw options		
		g_drawSprings = 1;
		g_drawCloth = false;
	}

	void checkDistanceBinding(){
		//find the closest pair ctd-ccd free
		//use a list ?
		int mini_binding_in = -1;
		int mini_binding_out = -1;
		int mini_binding_in_id = -1;
		int mini_binding_out_id = -1;
		float mini_distance = 999999.0f;

		//go other all instances
		for (int i = 0; i < int(cp->mInstances.size()); ++i)
		{
			NvFlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			int poffseti = cp->mInstances[i].mParticleOffset;
			for (int i = 0; i < binding_in.size(); i++){
				mNodes[i].mini_distance[i] = 9999.0f;
			}
			for (int j = 0; j < int(cp->mInstances.size()); j++) {
				if (i == j) continue;
				int poffsetj = cp->mInstances[j].mParticleOffset;
				for (int o = 0; o < binding_in.size(); o++){ //ccd
					//for (int p = 0; p < binding_out.size(); p++){ //ctd
					int p = o;
					float D = getMinDistances(i, j, o, p);
					if ((mNodes[i].binding_in_occupied[o] != -1) || (mNodes[j].binding_out_occupied[p] != -1)){ continue; }
					//need to avoid 1-1 and 2-2 so check if i already bound to j in the other binding site
					if (checkOccupation(i, j, o) != -1) continue;
					if (D < mini_distance){
						if (D < threshold_binding){
							//check if free
							//if ((mIntegrase[i].ccd_occupied[o]==-1) && (mIntegrase[j].ctd_occupied[p]==-1)){
							//	if ((mIntegrase[i].ccd_occupied[1 - o] != j) && (mIntegrase[j].ctd_occupied[1 - p] != i)){
							mini_distance = D;
							mini_binding_in = o;
							mini_binding_out = p;
							mini_binding_in_id = i;
							mini_binding_out_id = j;
							//	}
							//}
						}
					}
					//	}
					//}
					//}
				}
			}
			/*if ((mIntegrase[i].ccd_occupied[0] == -1) || (mIntegrase[i].ccd_occupied[1] == -1)){
			cout << endl;
			cout << mIntegrase[i].mini_distance[0] << " " << mIntegrase[i].mini_distance[1] << " " << mIntegrase[i].mini_indices[0] << " " << mIntegrase[i].mini_indices[1] << " " << mini_distance << endl;
			cout << mIntegrase[i].ccd_occupied[0] << " " << mIntegrase[i].ccd_occupied[1] << " " << mIntegrase[i].ctd_occupied[0] << " " << mIntegrase[i].ctd_occupied[1] << endl;
			cout << mIntegrase[i].ccd_ctdind[0] << " " << mIntegrase[i].ccd_ctdind[1] << " " << mIntegrase[i].ctd_ccdind[0] << " " << mIntegrase[i].ctd_ccdind[1] << endl;
			}*/
		}
		if (mini_binding_in_id != -1){
			boundTwoIds(mini_binding_in_id, mini_binding_out_id, mini_binding_in, mini_binding_out);
			//free_mIntegrase_CCD.erase(free_mIntegrase_CCD.begin() + mini_instance_id1);
			//free_mIntegrase_CTD.erase(free_mIntegrase_CTD.begin() + mini_instance_id2);
		}
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
		NvFlexGetRigidTransforms(g_flex, quat.buffer, pos.buffer);
		std:string name = "../../data/pack_result" + std::to_string(exp) + "_" + std::to_string(run);
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
		NvFlexSetSprings(g_flex, g_buffers->springIndices.buffer, g_buffers->springLengths.buffer, g_buffers->springStiffness.buffer, g_buffers->springLengths.size());
		
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
	
	virtual void Sync()
	{
		NvFlexSetSprings(g_flex, g_buffers->springIndices.buffer, g_buffers->springLengths.buffer, g_buffers->springStiffness.buffer, g_buffers->springLengths.size());
		NvFlexSetRigids(g_flex, g_buffers->rigidOffsets.buffer, g_buffers->rigidIndices.buffer, g_buffers->rigidLocalPositions.buffer, g_buffers->rigidLocalNormals.buffer, g_buffers->rigidCoefficients.buffer, g_buffers->rigidRotations.buffer, g_buffers->rigidTranslations.buffer, g_buffers->rigidOffsets.size() - 1, g_buffers->rigidIndices.size());
	}

	virtual void PostInitialize()
	{
		g_sceneLower = Vec3(-5.0f, 0.0f, 0.0f);
		g_sceneUpper = g_sceneLower + Vec3(10.0f, 10.0f, 5.0f);
	}

	void Update()
	{
		if (next_step){
			do_nextStep();
			
			next_step = false;
		}
		if (bindTwo) {
			cout << " 290 " << g_buffers->positions[290].x << endl;
			addOneRandomBinding();
			bindTwo = false;
			//NvFlexSetSprings(g_flex, g_buffers->springIndices.buffer, g_buffers->springLengths.buffer, g_buffers->springStiffness.buffer, g_buffers->springLengths.size());
		}
		return;
		NvFlexTimers atimers;
		float ** p;
		float ** v;
		int ** ph;
		float ** no;
		int a;
		//distance 71-3
		//distance center->farest dimer
		//reportBinding();

		//if (cp->comp_shape.size())
		//		flexExtSetRigidTransformation(cp->fcontainer, cp->mInstances.size());
		if ((g_frame % 5) == 0)
		{
			if (dojitter) cp->jitter(!dojitter_steared);
			//cout << "free " << free_mIntegrase_CCD.size() << " " << free_mIntegrase_CTD.size() << endl;
			//if (free_mIntegrase_CCD.size() != 0) addOneRandomBinding();
		}
		checkDistanceBinding();
		NvFlexSetSprings(g_flex, g_buffers->springIndices.buffer, g_buffers->springLengths.buffer, g_buffers->springStiffness.buffer, g_buffers->springLengths.size());
		//save packing leftHand
		if ((g_frame % 1) == 0)
		{
			if (write_output){
				std:string pfix = "";// "_" + std::to_string(current_Exp) + "_" + std::to_string(current_run);
				cp->writeToBinary(pfix);
			}
		}
		//checkRadiusAggregate();
		//reportBinding(false);
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
					//writeReport(current_Exp, current_run, free, true);			//connectivity
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
};
