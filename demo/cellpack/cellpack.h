#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sys/types.h> 
#include <string>
#include "halton.hpp"

#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#define ChangeDir _chdir
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#define ChangeDir chdir
#endif

#include "include/json/json.h"

/*
//GRAPH stuff
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/layered/DfsAcyclicSubgraph.h>
#include <ogdf/basic/List.h>
#include <ogdf/graphalg/ShortestPathAlgorithms.h>
#include <ogdf\graphalg\Dijkstra.h>
#include <ogdf\graphalg\ShortestPathWithBFM.h>
#include <ogdf\graphalg\Voronoi.h>
#include <ogdf\decomposition\BCTree.h>
#include <ogdf/fileformats/GmlParser.h>

using namespace ogdf;
*/

char cCurrentPath[FILENAME_MAX];
const int MAXALPH = 26;
string getRandomString(int n)
{
	char alphabet[MAXALPH] = { 'A', 'B', 'C', 'D', 'E', 'F', 'G',
						  'H', 'I', 'J', 'K', 'L', 'M', 'N',
						  'O', 'P', 'Q', 'R', 'S', 'T', 'U',
						  'V', 'W', 'X', 'Y', 'Z' };

	string res = "";
	for (int i = 0; i < n; i++)
		res = res + alphabet[rand() % MAXALPH];
	return res;
}

std::vector<std::string> getThreeLetterString(int amount)
{
	char alphabet[MAXALPH] = { 'A', 'B', 'C', 'D', 'E', 'F', 'G',
						  'H', 'I', 'J', 'K', 'L', 'M', 'N',
						  'O', 'P', 'Q', 'R', 'S', 'T', 'U',
						  'V', 'W', 'X', 'Y', 'Z' };
	std::vector<std::string> allcombination;
	for (char a = 'A'; a <= 'Z'; ++a) {
		for (char b = 'A'; b <= 'Z'; ++b) {
			//if (b == a) continue;
			for (char c = 'A'; c <= 'Z'; ++c) {
				//if (c == a) continue;
				//if (c == b) continue;
				//std::cout << a << b << c << '\n';
				std::stringstream ss;
				ss << a << b << c;
				allcombination.push_back(ss.str());
				if (allcombination.size() == amount) {
					return allcombination;
				}
			}
		}
	}
	return allcombination;
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = str.find(from);
	if (start_pos == std::string::npos)
		return false;
	str.replace(start_pos, from.length(), to);
	return true;
}


float SampleSDFX(float* sdf, int dim, int x, int y, int z)
{
	if (x >= dim || x < 0) return 9999;
	if (y >= dim || y < 0) return 9999;
	if (z >= dim || z < 0) return 9999;
	return sdf[z*dim*dim + y*dim + x];
}



NvFlexDistanceFieldId CreateSDFFromMeshF(const char* meshFile, Mesh* mesh, int dim, float margin = 0.1f, float expand = 0.0f, float scale = 1.0f)
{
	//Mesh* mesh = ImportMesh(meshFile);
	// include small margin to ensure valid gradients near the boundary
	//mesh->Normalize(1.0f - margin);
	//mesh->Transform(TranslationMatrix(Point3(margin, margin, margin)*0.5f));

	Vec3 lower(0.0f);
	Vec3 upper(1.0f);
	mesh->GetBounds(lower, upper);
	
	// try and load the sdf from disc if it exists
	// Begin Add Android Support
#ifdef ANDROID
	string sdfFile = string(meshFile, strlen(meshFile) - strlen(strrchr(meshFile, '.'))) + ".pfm";
#else
	string sdfFile = string(meshFile, strrchr(meshFile, '.')) + ".pfm";
#endif
	// End Add Android Support

	PfmImage pfm;
	if (!PfmLoad(sdfFile.c_str(), pfm))
	{
		pfm.m_width = dim;
		pfm.m_height = dim;
		pfm.m_depth = dim;
		pfm.m_data = new float[dim*dim*dim];

		printf("Cooking SDF: %s - dim: %d^3\n", sdfFile.c_str(), dim);

		CreateSDF(mesh, dim, lower, upper, pfm.m_data);

		PfmSave(sdfFile.c_str(), pfm);
	}

	//printf("Loaded SDF, %d\n", pfm.m_width);

	assert(pfm.m_width == pfm.m_height && pfm.m_width == pfm.m_depth);

	// cheap collision offset
	int numVoxels = int(pfm.m_width*pfm.m_height*pfm.m_depth);
	for (int i = 0; i < numVoxels; ++i)
	{
		pfm.m_data[i] += expand;
		pfm.m_data[i] *= scale;
	}

	NvFlexVector<float> field(g_flexLib);
	field.assign(pfm.m_data, pfm.m_width*pfm.m_height*pfm.m_depth);
	field.unmap();

	// set up flex collision shape
	NvFlexDistanceFieldId sdf = NvFlexCreateDistanceField(g_flexLib);
	NvFlexUpdateDistanceField(g_flexLib, sdf, dim, dim, dim, field.buffer);

	// entry in the collision->render map
	g_fields[sdf] = CreateGpuMesh(mesh);

	delete mesh;
	delete[] pfm.m_data;

	return sdf;
}


PfmImage CreateSDFFromMesh(const char* meshFile, Mesh* mesh, int dim, float margin = 0.1f, float expand = 0.0f)
{
	Vec3 lower(0.0f);
	Vec3 upper(1.0f);
	mesh->GetBounds(lower, upper);
	// try and load the sdf from disc if it exists
	// Begin Add Android Support
#ifdef ANDROID
	string sdfFile = string(meshFile, strlen(meshFile) - strlen(strrchr(meshFile, '.'))) + ".pfm";
#else
	string sdfFile = string(meshFile, strrchr(meshFile, '.')) + ".pfm";
#endif
	// End Add Android Support

	PfmImage pfm;
	if (!PfmLoad(sdfFile.c_str(), pfm))
	{
		pfm.m_width = dim;
		pfm.m_height = dim;
		pfm.m_depth = dim;
		pfm.m_data = new float[dim*dim*dim];

		printf("Cooking SDF: %s - dim: %d^3\n", sdfFile.c_str(), dim);

		CreateSDF(mesh, dim, lower, upper, pfm.m_data);

		PfmSave(sdfFile.c_str(), pfm);
	}

	//printf("Loaded SDF, %d\n", pfm.m_width);

	assert(pfm.m_width == pfm.m_height && pfm.m_width == pfm.m_depth);

	// cheap collision offset
	int numVoxels = int(pfm.m_width*pfm.m_height*pfm.m_depth);
	for (int i = 0; i < numVoxels; ++i)
		pfm.m_data[i] += expand;

	//NvFlexVector<float> field(g_solverLib);
	//field.assign(pfm.m_data, pfm.m_width*pfm.m_height*pfm.m_depth);
	//field.unmap();

	// set up flex collision shape
	//NvFlexDistanceFieldId sdf = NvFlexCreateDistanceField(g_solverLib);
	//NvFlexUpdateDistanceField(g_solverLib, sdf, dim, dim, dim, field.buffer);

	// entry in the collision->render map
	//g_fields[sdf] = CreateGpuMesh(mesh);

	//delete mesh;

	return pfm;
}

NvFlexExtAsset* flexExtCreateRigidFromPoints(std::vector<Vec3> points, bool to_center = true)
{
	Vec3 center = Vec3(0.0f,0.0f,0.0f);
	std::vector<Vec4> particles;

	NvFlexExtAsset* asset = new NvFlexExtAsset();
	memset(asset, 0, sizeof(*asset));

	if (points.size())
	{
		for (int i = 0; i < points.size(); i++)
		{
			center += points[i];
			particles.push_back(Vec4(points[i].x, points[i].y, points[i].z, 1.0f));
		}
		const int numParticles = int(particles.size());
		asset->numParticles = numParticles;
		asset->maxParticles = numParticles;

		asset->particles = new float[numParticles * 4];
		memcpy(asset->particles, &particles[0], sizeof(Vec4)*numParticles);

		// store center of mass should actually be 0,0,0
		center /= float(numParticles);
		//center = Vec3(0.00001f,0.00001f,0.00001f);
		//if (!to_center) center = Vec3(0.00001f, 0.00001f, 0.00001f);
		cout << " center was " << center.x << " " << center.y << " " << center.z << endl;
		asset->numShapes = 1;
		asset->numShapeIndices = numParticles;

		// for rigids we just reference all particles in the shape
		asset->shapeIndices = new int[numParticles];

		for (int i = 0; i < numParticles; ++i)
			asset->shapeIndices[i] = i;

		asset->shapeCenters = new float[4];
		asset->shapeCenters[0] = center.x;
		asset->shapeCenters[1] = center.y;
		asset->shapeCenters[2] = center.z;
		//asset->shapeCenters[3] = 0;

		asset->shapePlasticThresholds = NULL;
		asset->shapePlasticCreeps = NULL;

		asset->shapeCoefficients = new float[1];
		asset->shapeCoefficients[0] = 1.0f;

		asset->shapeOffsets = new int[1];
		asset->shapeOffsets[0] = numParticles;
	}
	return asset;
}


void RemoveSpring(int i, set<int> toRemove)
{
	//copy
	std::vector<float> spL;
	std::vector<float> spS;
	std::vector<int>   spI;
	int sp_counter = 0;
	//remove all spring of length Length , should be the Persistence Length 
	for (int m = 0; m < g_buffers->springLengths.size(); m++){
		if ((toRemove.find(g_buffers->springIndices[sp_counter]) != toRemove.end()) || (toRemove.find(g_buffers->springIndices[sp_counter + 1]) != toRemove.end()))
		{
			sp_counter += 2;
		}
		else {
			spL.push_back(g_buffers->springLengths[m]);
			spS.push_back(g_buffers->springStiffness[m]);
			spI.push_back(g_buffers->springIndices[sp_counter]);
			spI.push_back(g_buffers->springIndices[sp_counter + 1]);
			sp_counter += 2;
		}
	}
	g_buffers->springIndices.assign(&spI[0], spI.size());
	g_buffers->springStiffness.assign(&spS[0], spS.size());
	g_buffers->springLengths.assign(&spL[0], spL.size());
}

bool foundPair(Vec2 pair, std::set<Vec2> toRemove){
	//cout << "test " << pair.x << " " << pair.y << endl;
	for (auto n : toRemove) {
		if ((n.x == pair.x) && (n.y == pair.y))
			return true;
	}
	return false;
}

void RemoveSpringPair(std::set<Vec2> toRemove)
{
	//copy
	std::vector<float> spL;
	std::vector<float> spS;
	std::vector<int>   spI;
	int sp_counter = 0;
	int found = 0;
	//remove all spring of length Length , should be the Persistence Length 
	for (int m = 0; m < g_buffers->springLengths.size(); m++){
		Vec2 pair = Vec2( g_buffers->springIndices[sp_counter], g_buffers->springIndices[sp_counter + 1] );
		if (foundPair(pair,toRemove))
		{
			found++;
			cout << "found one" << found << endl;
			sp_counter += 2;
		}
		else {
			spL.push_back(g_buffers->springLengths[m]);
			spS.push_back(g_buffers->springStiffness[m]);
			spI.push_back(g_buffers->springIndices[sp_counter]);
			spI.push_back(g_buffers->springIndices[sp_counter + 1]);
			sp_counter += 2;
		}
	}
	cout << "found " << found << endl;
	g_buffers->springIndices.assign(&spI[0], spI.size());
	g_buffers->springStiffness.assign(&spS[0], spS.size());
	g_buffers->springLengths.assign(&spL[0], spL.size());
}

void ReplaceSpring(int i, int j, int k, int l){
	//go through all spring and find i,j ?
	// or generate a new array storing pair information id
	int sp_counter = 0;
	//printf ("search for %i and %i to be replace by %i and %i\n",i,j,k,l);
	for (int m = 0; m<g_buffers->springLengths.size(); m++){
		if (g_buffers->springIndices[sp_counter] == i){
			//printf ("found i  j is %i\n",g_springIndices[sp_counter+1]);
			if (g_buffers->springStiffness[sp_counter + 1] == j){
				//printf ("found the pair i %i j %i\n",i,j);
				g_buffers->springIndices[sp_counter] = k;
				g_buffers->springIndices[sp_counter + 1] = l;
				break;
			}
		}
		/*else {
		if (g_springIndices[sp_counter] == j){
		printf ("found j but i is %i\n",g_springIndices[sp_counter+1]);
		if (g_springIndices[sp_counter+1] == i){
		printf ("found the pair j %i i %i\n",j,i);
		g_springIndices[sp_counter] = l;
		g_springIndices[sp_counter+1] = k;
		break;
		}
		}
		}*/
		sp_counter += 2;
	}
}

void CreateSpringInter(int i, int j, float stiffness, float give = 0.0f, float length = 0.0f)
{
	//printf(" %i and %i \n", i, j);
	g_buffers->springIndices.push_back(i);
	g_buffers->springIndices.push_back(j);
	if (length != 0.0f) g_buffers->springLengths.push_back(length);
	else g_buffers->springLengths.push_back((1.0f + give)*Length(Vec3(g_buffers->positions[i]) - Vec3(g_buffers->positions[j])));
	g_buffers->springStiffness.push_back(stiffness);
}

void CreatePersistence(Rope& rope, int current, int persistence, float stiffness, int nfloat, float give, float D)
{
	//D = 150.0f * 1.0f / 200.0f;
	float st = stiffness;
	for (int j = 1; j < persistence + 1; j++){
		float r = Randf((D / 100.0f)*-1.0f, 0.0f);
		if (j == 1) st = -1.0f;
		else st = stiffness;
		// cout << "create link " << current - j << " " << current << " " << rope.mIndices.size() << endl;
		if (rope.mIndices.size() > j) {
			CreateSpringInter(current - j, current, st, give, (D*(float)j));// +r*(float)(j - 1));
		}
	}
}

void CreatePersistenceDisc(Rope& rope, int current, int persistence, float stiffness, int nfloat, float give, float D)
{
	//do it within i, i+persistence
	if (persistence == 0) return;
	if (rope.mIndices.size() <= persistence) return;
	if (current < persistence) return;
	float r = Randf(-D/10.0f,D/10.0f);// Randf((D / 100.0f)*-1.0f, 0.0f);
	CreateSpringInter(current - persistence, current, stiffness, give, (D*(float)persistence) + r);
}


void CreateClosedPersistence(Rope& rope, int start, int end, int persistence, float stiffness, int nfloat, float give, float D)
{
	//Lv1 = i+1
	//Lv2 = i+2
	//Lv3 = i+3
	//Lv4 = i+4
	for (int l = 0; l < persistence; l++){
		float r = Randf((D / 100.0f)*-1.0f, 0.0f);;// float r = Randf(-1.0f, 1.0f) / 2000.0f;
		for (int k = l, i = 0; k >= 0, i < l + 1; k--, i++)
		{
			CreateSpringInter(start + i, end - k, stiffness, give, (D*(float)(l + 1)) + r*(float)(l));
			cout << l << "start " << i << " end " << k << " " << l << " " << D << " " << D*(float)(l + 1) << endl;
		}
	}
}

void CreateClosedPersistenceDisc(Rope& rope, int start, int end, int persistence, float stiffness, int nfloat, float give, float D)
{
	if (persistence == 0) return;
	int k = persistence;
	for (int  i = 0; i < k;  i++)
	{
		float r = Randf(-D / 10.0f, D / 10.0f);// Randf((D / 100.0f)*-1.0f, 0.0f);
		CreateSpringInter(end - (k-1) + i , start + i, stiffness, give, (D*(float)persistence) + r);
		cout << "start " << i << " end " << k << " " << end - (k - 1) + i << " " << start + i << " " << D << " " << (D*(float)persistence) + r << endl;
	}
}


void CreateInsertedPersistence(Rope rope, int ip, float stiffness, float give, float D)
{
	int a, b, k, l;
	//cout << "fix spring around " << ip << " " << rope.persistence << endl;
	for (int lvl = 0; lvl < rope.persistence; lvl++){
		float r = Randf((D / 100.0f)*-1.0f, 0.0f); //float r = Randf(-1.0f, 1.0f) / 2000.0f;
		for (int i = 0; i < lvl + 1; i++)
		{
			a = rope.mIndices[ip - (lvl - i)];
			b = rope.mIndices[ip + (i + 2)];
			k = a;
			l = rope.mIndices[ip + (i + 1)];
			ReplaceSpring(a, b, k, l);
			//cout << a << " " << b << " " << k << " " << l << endl;
			//cout << lvl << " replace " << ip - (lvl - i) << " / " << ip + (i + 2) << " by " << ip - (lvl - i) << " / " << ip + (i + 1) << endl;
		}
		//cout << lvl << "create  " << ip + 1 << " / " << ip + (lvl + 2) << " D " << (D*(float)(lvl+1)) + r*(float)(lvl) << endl;
		CreateSpringInter(rope.mIndices[ip + 1], rope.mIndices[ip + (lvl + 2)], stiffness, give, (D*(float)(lvl + 1)) + r*(float)(lvl));
	}
}

void CreateRopeFromData(Rope& rope, float scale, float stiffness, float *data,
	float length, int nfloat, int phase, float spiralAngle = 0.0f,
	float invmass = 1.0f, float give = 0.075f, int extend_nb = 57,
	bool extend = false, bool closed = true, float D = 0.0f, int persistence = 2)
{
	//stiffness = 1.0f;
	int start = int(g_buffers->positions.size());
	if (give != 0.0f) give = 0.0f;
	float r = 1.0f;//biased on the 1-3 spring
	int current = 0;
	//if closed do the last point ?
	int longP = 0;// 40;
	float longstiffness = 0;
	cout << "nfloat " << nfloat << " extend " << extend<< " " << extend_nb <<" invmass " << invmass <<endl;
	if (extend_nb == 0)extend_nb = 1;
	for (int i = 0; i < nfloat; i += 3)
	{
		//if (i/3>20) break;
	    //printf("add a point %i %f %f %f \n",i/3,data[i] * scale, data[i+1] * scale, data[i+2] * scale);
		int begin = int(g_buffers->positions.size());
		int prev = begin;//int(g_positions.size())-1;

		if (!extend){
			current = begin;
			rope.mIndices.push_back(begin);
			g_buffers->positions.push_back(Vec4(data[i] * scale, data[i + 1] * scale, data[i + 2] * scale, invmass));
			g_buffers->normals.push_back(Vec4(0.0f, 0.0f, 1.0f, 1.0f));
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(phase);//int(g_positions.size()));
			// cout << " should bind to next point ? " << begin << " D " << D << endl;
			CreatePersistence(rope, begin, persistence, stiffness, nfloat, give, D);
			CreatePersistenceDisc(rope, begin, longP, longstiffness, nfloat, give, D);//dna specfici test
		}
		else {
			Vec3 point = Vec3(data[i] * scale, data[i + 1] * scale, data[i + 2] * scale);
			Vec4 next_point = Vec4(0, 0, 0, 0);
			Vec3 dirtopt = Vec3(0, 0, 0);
			if (i + 5 < nfloat) {
				next_point = Vec4(data[i + 3] * scale, data[i + 4] * scale, data[i + 5] * scale, 1.0f);
			}
			else {
				next_point = Vec4(data[0] * scale, data[1] * scale, data[2] * scale, 1.0f);
			}
			dirtopt = Vec3(next_point.x - point.x, next_point.y - point.y, next_point.z - point.z);

			for (int j = 0; j<extend_nb; j++){
				current = int(g_buffers->positions.size());
				float perc = ((float)j / (float)extend_nb);
				rope.mIndices.push_back(int(g_buffers->positions.size()));

				g_buffers->positions.push_back(Vec4(point.x + dirtopt.x*perc,
					point.y + dirtopt.y*perc,
					point.z + dirtopt.z*perc,
					invmass));
				g_buffers->normals.push_back(Vec4(0.0f, 0.0f, 1.0f, 1.0f));
				g_buffers->velocities.push_back(0.0f);
				g_buffers->phases.push_back(phase);//int(g_positions.size()));
				//printf ("subdivide how many %i %i %i %i\n",i,j,extend_nb,current);
				CreatePersistence(rope, current, persistence, stiffness, nfloat, give, D);
				CreatePersistenceDisc(rope, current, longP, longstiffness, nfloat, give, D);//dna specfici test
			}
		}
	}
	if (closed) {
		//cout << "nb poitns " <<rope.mIndices.size() << endl;
		r = Randf(-1.0f, 1.0f) / 2000.0f;
		int startindex = rope.mIndices[0];
		int endindex = int(g_buffers->positions.size() - 1);// rope.mIndices[rope.mIndices.size() - 1];
		CreateClosedPersistence(rope, startindex, endindex, persistence, stiffness, nfloat, give, D);
		CreateClosedPersistenceDisc(rope, startindex, endindex, longP, longstiffness, nfloat, give, D);
	}
	rope.coarseIndices = rope.mIndices;
}

void CreateRopeFromData(Rope& rope, float scale, float stiffness, std::vector<float> data,
	float length, int nfloat, int phase, float spiralAngle = 0.0f,
	float invmass = 1.0f, float give = 0.075f, int extend_nb = 57,
	bool extend = false, bool closed = true, float D = 0.0f, int persistence = 2)
{
	//stiffness = 1.0f;
	int start = int(g_buffers->positions.size());
	if (give != 0.0f) give = 0.0f;
	float r = 1.0f;//biased on the 1-3 spring
	int current = 0;
	//if closed do the last point ?
	int longP = 0;// 40;
	float longstiffness = 0;
	//cout << "nfloat " << nfloat << " extend " << extend << " " << extend_nb << endl;
	if (extend_nb == 0)extend_nb = 1;
	for (int i = 0; i < nfloat; i += 3)
	{
		//if (i/3>20) break;
		//printf("add a point %i %f %f %f \n", i / 3, data[i] * scale, data[i + 1] * scale, data[i + 2] * scale);
		int begin = int(g_buffers->positions.size());
		int prev = begin;//int(g_positions.size())-1;

		if (!extend) {
			current = begin;
			rope.mIndices.push_back(begin);
			g_buffers->positions.push_back(Vec4(data[i] * scale, data[i + 1] * scale, data[i + 2] * scale, 1.0f));
			g_buffers->normals.push_back(Vec4(0.0f, 0.0f, 1.0f, 1.0f));
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(phase);//int(g_positions.size()));
			CreatePersistence(rope, begin, persistence, stiffness, nfloat, give, D);
			CreatePersistenceDisc(rope, begin, longP, longstiffness, nfloat, give, D);//dna specfici test
		}
		else {
			Vec3 point = Vec3(data[i] * scale, data[i + 1] * scale, data[i + 2] * scale);
			Vec4 next_point = Vec4(0, 0, 0, 0);
			Vec3 dirtopt = Vec3(0, 0, 0);
			if (i + 5 < nfloat) {
				next_point = Vec4(data[i + 3] * scale, data[i + 4] * scale, data[i + 5] * scale, 1.0f);
			}
			else {
				next_point = Vec4(data[0] * scale, data[1] * scale, data[2] * scale, 1.0f);
			}
			dirtopt = Vec3(next_point.x - point.x, next_point.y - point.y, next_point.z - point.z);

			for (int j = 0; j<extend_nb; j++) {
				current = int(g_buffers->positions.size());
				float perc = ((float)j / (float)extend_nb);
				rope.mIndices.push_back(int(g_buffers->positions.size()));

				g_buffers->positions.push_back(Vec4(point.x + dirtopt.x*perc,
					point.y + dirtopt.y*perc,
					point.z + dirtopt.z*perc,
					1.0f));
				g_buffers->normals.push_back(Vec4(0.0f, 0.0f, 1.0f, 1.0f));
				g_buffers->velocities.push_back(0.0f);
				g_buffers->phases.push_back(phase);//int(g_positions.size()));
												   //printf ("subdivide how many %i %i %i %i\n",i,j,extend_nb,current);
				CreatePersistence(rope, current, persistence, stiffness, nfloat, give, D);
				CreatePersistenceDisc(rope, current, longP, longstiffness, nfloat, give, D);//dna specfici test
			}
		}
	}
	if (closed) {
		//cout << "nb poitns " << rope.mIndices.size() << endl;
		r = Randf(-1.0f, 1.0f) / 2000.0f;
		int startindex = rope.mIndices[0];
		int endindex = int(g_buffers->positions.size() - 1);// rope.mIndices[rope.mIndices.size() - 1];
		CreateClosedPersistence(rope, startindex, endindex, persistence, stiffness, nfloat, give, D);
		CreateClosedPersistenceDisc(rope, startindex, endindex, longP, longstiffness, nfloat, give, D);
	}
	rope.coarseIndices = rope.mIndices;
}

int HaltonSample3D(float* radius, float separation, Vec3* points, int maxPoints){
	int N = maxPoints;
	int DIM_MAX = 3;
	int* base = new int[DIM_MAX];
	int i;
	int j;
	int* leap = new int[DIM_MAX];
	double* r = new double[DIM_MAX*N];
	int* seed = new int[DIM_MAX];
	int step;
	halton_dim_num_set(DIM_MAX);
	step = 0;
	halton_step_set(step);
	for (i = 0; i < DIM_MAX; i++)
	{
		seed[i] = 0;
	}
	halton_seed_set(seed);
	for (i = 0; i < DIM_MAX; i++)
	{
		base[i] = prime(i + 1);
	}
	halton_base_set(base);

	i4vec_transpose_print(DIM_MAX, seed, "  SEED = ");
	i4vec_transpose_print(DIM_MAX, base, "  BASE = ");
	halton_sequence(N, r);
	for (j = 0; j < N; j++)
	{
		for (i = 0; i < DIM_MAX; i++)
		{
			points[j][i] = r[i + j*DIM_MAX] * radius[i];//scale ?
		}
	}
	return 0;
}

class cellPACK
{
public:
	struct Instance
	{
		Vec3 mTranslation;
		Quat mRotation;
		float mLifetime;

		int mGroup;
		int mParticleOffset;

		int mMeshIndex;

		int bounded;
	};

	struct MeshBatch
	{
		GpuMesh* mMesh;
		NvFlexExtAsset* mAsset;

		std::vector<Matrix44> mInstanceTransforms;
	};

	struct MeshAsset
	{
		const char* file;
		float scale;
	};

	struct AssetBatch
	{
		GpuMesh* mMesh;
		NvFlexTriangleMeshId nvMesh;
		NvFlexExtAsset* mAsset;
		float offsetx;
		float offsety;
		float offsetz;
		float pcpalVectorx;
		float pcpalVectory;
		float pcpalVectorz;
		int compId;
		int nInstances;
		const char* ingr_name;
		std::vector<Matrix44> mInstanceTransforms;
	};

	struct IngredientSphereTree{
		int nbLevel;
		//radius..
		std::vector<int> LevelCounts;
		std::vector<Vec3> LevelPoints;
		std::vector<int> LevelStarts;
		std::vector<Vec2> LevelMapping;//start,count
		std::vector<int> LevelMappingStarts;
		int nBinding;
		std::vector<int> BindingStarts;//start,count
		std::vector<int> BindingSites;
		std::vector<float> DistancesMatrix;
	};

	struct IngredientPartner
	{
		int nPartner;
		std::vector<Json::Value> partner_nodes;
		std::vector<int> batchs_id;
		std::vector<int> instances_id;
	};

	std::vector<IngredientPartner> iPartners;
	std::vector<IngredientPartner> iPartnersFibers;
	
	bool use_partners_properties = false;
	bool use_instances_mesh;
	bool overwrite_radius = true;

	Json::Value book_json;
	Json::Value results_json;
	int lodproxy_to_use=0;
	int dna_persistence=5;
	int force_not_center = 0;
	float main_scale = 1.0f;
	float main_radius = 10.0f;
	string mainpath = "..\\..\\data\\cellpack\\";// "D:\\Data\\cellPAC_data\\cellPACK_database_1.1.0\\";//D:\Data\cellPAC_data\cellPACK_database_1.1.0
	string datapath = "..\\..\\data\\cellpack\\beads\\";//"D:\\Data\\cellPAC_data\\cellPACK_database_1.1.0\\other\\";
	string geompath = "..\\..\\data\\cellpack\\geoms\\";//"D:\\Data\\cellPAC_data\\cellPACK_database_1.1.0\\geometries\\";
	int iGroupCounter=0;
	Vec3 minExtents, maxExtents; //boudningBox
	int maxParticles;
	int numParticles = 0;
	std::vector<int> mask;
	std::vector<int> maks_protein;
	std::vector<int> maks_fiber;
	std::vector<int> mask_membrane;

	std::vector<int> fiber_togrow;
	std::vector<int> fiber_togrow_stpos;
	std::vector<int> fiber_togrow_length;

	std::vector<string> pnames;
	std::vector<string> pnames_fiber;
	std::vector<Json::Value> pnames_fiber_nodes;

	std::vector<AssetBatch> iBatches;
	std::vector<AssetBatch> iBatchesFiber;
	std::vector<IngredientSphereTree> mIngrSphereTree;
	std::vector<Instance> mInstances;
	std::vector<Instance> mInstancesFiber;

	std::vector<Json::Value> proteins_nodes;

	std::vector<NvFlexTriangleMeshId> comp_mesh;
	std::vector<Mesh*> comp_tri;

	float comp_radius = 1600.0f;
	bool write_output = false;
	bool dojitter_steared = true;
	bool dojitter = true;
	bool dojitter_biased = false;
	bool use_threshold_binding = true;
	float dojitter_strength = 0.01f;
	float dojitter_biased_strength = 1.0f;
	float threshold_binding;
	bool dosimulation = false;

	float Nsub = 10.0f;

	ofstream output;
	ofstream output_bin;

	NvFlexDistanceFieldId sdf;
	PfmImage sdfdata;
	int dim = 128;
	int totalNBMol = 0;
	bool use_rb = false;

	cellPACK(){
		GetCurrentDir(cCurrentPath, sizeof(cCurrentPath));
		cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */
		string path = string(cCurrentPath) + "\\..\\..\\data\\cellpack\\";
		mainpath = path;
		datapath = mainpath + "\\beads\\";
		geompath = mainpath + "\\geoms\\";
	};
	

	float getValues(float i, float j, float k) {
		float u = (k * dim * dim) + (j * dim) + i;
		if (u > dim * dim * dim) u = dim * dim * dim;//or 0?
		if (u < 0) u = 0;
		return sdfdata.m_data[(int)u];
	}


	float trilinearInterpolation(Vec3 apoint, Vec3 bbmin, float grid_unit) {
		// Find the x, y and z values of the
		// 8 vertices of the cube that surrounds the point
		Vec3 p = (apoint - bbmin)*(1.0 / grid_unit);//+/-0.5 ?
		return SampleSDFX(sdfdata.m_data, dim, (int)p.x, (int)p.y, (int)p.z);
		/*float x0 = floor(p.x);
		float x1 = floor(p.x) + 1.0;
		float y0 = floor(p.y);
		float y1 = floor(p.y) + 1.0;
		float z0 = floor(p.z);
		float z1 = floor(p.z) + 1.0;
		// Look up the values of the 8 points surroundng the cube
		// Find the weights for each dimension
		float x = (p.x - x0);
		float y = (p.y - y0);
		float z = (p.z - z0);
		float V000 = getValues(x0, y0, z0);
		float V100 = getValues(x1, y0, z0);
		float V010 = getValues(x0, y1, z0);
		float V001 = getValues(x0, y0, z1);
		float V101 = getValues(x1, y0, z1);
		float V011 = getValues(x0, y1, z1);
		float V110 = getValues(x1, y1, z0);
		float V111 = getValues(x1, y1, z1);
		float Vxyz = V000 * (1.0 - x)*(1.0 - y)*(1.0 - z) +
			V100 * x*(1.0 - y)*(1.0 - z) +
			V010 * (1.0 - x)*y*(1.0 - z) +
			V001 * (1.0 - x)*(1.0 - y)*z +
			V101 * x*(1.0 - y)*z +
			V011 * (1.0 - x)*y*z +
			V110 * x*y*(1.0 - z) +
			V111 * x*y*z;
		return Vxyz;*/
	}

	Vec3 CalculateSurfaceNormal(Vec3 p, Vec3 bbmin, float grid_unit)
	{
		float H = grid_unit; //1.0f/grid_unit;// 0.001f;
		float dx = trilinearInterpolation(p + Vec3(H, 0.0f, 0.0f), bbmin, grid_unit) - trilinearInterpolation(p - Vec3(H, 0.0f, 0.0f), bbmin, grid_unit);
		float dy = trilinearInterpolation(p + Vec3(0.0f, H, 0.0f), bbmin, grid_unit) - trilinearInterpolation(p - Vec3(0.0f, H, 0.0f), bbmin, grid_unit);
		float dz = trilinearInterpolation(p + Vec3(0.0f, 0.0f, H), bbmin, grid_unit) - trilinearInterpolation(p - Vec3(0.0f, 0.0f, H), bbmin, grid_unit);

		return Normalize(Vec3(dx, dy, dz));
	}

	Vec3 GetAwayFromSurface(Vec3 p, int compId) {
		Vec3 lower = comp_tri[compId]->m_minExtents;// (0.0f);
		Vec3 upper = comp_tri[compId]->m_maxExtents;//(1.0f);
		// comp_tri[compId]->GetBounds(lower, upper);
		Vec3 botsdf = Vec3(lower[0], lower[1], lower[2]);
		float spacing = (upper[0] - lower[0]) / (float)dim;
		return CalculateSurfaceNormal(p, botsdf, spacing);
	}

	float GetSurfaceDistance(Vec3 p, int compId) {
		Vec3 lower = comp_tri[compId]->m_minExtents;// (0.0f);
		Vec3 upper = comp_tri[compId]->m_maxExtents;//(1.0f);
		// comp_tri[compId]->GetBounds(lower, upper);
		// cout << " bounds ? " << lower[0] << " " << lower[1] << " " << lower[2] << endl;
		Vec3 botsdf = Vec3(lower[0], lower[1], lower[2]);
		float spacing = (upper[0] - lower[0]) / (float)dim;
		return trilinearInterpolation(p, botsdf, spacing);
	}

	bool isInside(Vec3 p, int compId){
		//find inside point using sdf
		Vec3 lower = comp_tri[compId]->m_minExtents;// (0.0f);
		Vec3 upper = comp_tri[compId]->m_maxExtents;//(1.0f);
		// comp_tri[compId]->GetBounds(lower, upper);
		float spacing = (upper[0] - lower[0]) / (float)dim;
		Vec3 botsdf = Vec3(lower[0], lower[1], lower[2]);
		Vec3 topsdf = Vec3(upper[0], upper[1], upper[2]);
		Vec3 pquery = Vec3((p.x - botsdf.x) / spacing, (p.y - botsdf.y) / spacing, (p.z - botsdf.z) / spacing);
		Vec3 pquery_index = Vec3((int)pquery.x, (int)pquery.y, (int)pquery.z);
		float D = SampleSDFX(sdfdata.m_data, dim, (int)pquery.x, (int)pquery.y, (int)pquery.z);
		return (D < 0);
	}

	IngredientSphereTree parseProxy(string filename)
	{
		IngredientSphereTree ingr_spheres;
		//check the extension type
		cout << "parsing " << filename << endl;
		if (filename.substr(filename.find_last_of(".") + 1) == "bin") {
			cout << "parsing bin " << filename << endl;
			std::ifstream ifs(filename, std::ios::binary);
			if (ifs.is_open()) {
				int l;
				ifs.read(reinterpret_cast<char*>(&l), sizeof(int));
				ingr_spheres.nbLevel = l;
				cout << "found nbLevel " << ingr_spheres.nbLevel << endl;
				for (int i = 0; i < ingr_spheres.nbLevel; i++){
					int tmp;
					ifs.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
					ingr_spheres.LevelCounts.push_back(tmp);
					cout << "parsing level "<< i << " Cluster Counts " << tmp << endl;
				}
				//do we have binding info
				int b;
				ifs.read(reinterpret_cast<char*>(&b), sizeof(int));
				ingr_spheres.nBinding = b;
				cout << "found nBinding " << ingr_spheres.nBinding << endl;
				//gather binding site indices
				int start = 0;
				for (int i = 0; i < ingr_spheres.nBinding; i++){
					int count;
					ifs.read(reinterpret_cast<char*>(&count), sizeof(count));
					ingr_spheres.BindingStarts.push_back(start);
					ingr_spheres.BindingStarts.push_back(count);
					start += count;
				}
				//gather coordinates per level
				for (int i = 0; i < ingr_spheres.nbLevel; i++){
					ingr_spheres.LevelStarts.push_back(ingr_spheres.LevelPoints.size());
					for (int j = 0; j < ingr_spheres.LevelCounts[i]; j++){
						Vec4 tmp;
						ifs.read(reinterpret_cast<char*>(&tmp[0]), sizeof(float) * 4);
						ingr_spheres.LevelPoints.push_back(Vec3(tmp.x, tmp.y, tmp.z));
						//cout << tmp.x << " " << tmp.y << " " << tmp.z << endl;
						//cout << "parsing " << ingr_spheres.LevelCounts[i] << " points nb " << ingr_spheres.LevelPoints.size() << endl;
					}
					cout << "parsing level " << i << " with " << ingr_spheres.LevelCounts[i] << " cluster and points nb " << ingr_spheres.LevelPoints.size() << endl;
				}
				//gather mapping
				/*
				mapping.extend(mappingL1_atom_order)#L1_size *2
				mapping.extend(mappingL2_atom_order)#L2_size *2
				mapping.extend(mappingL2_L1_order)#L2_size *2
				*/
				for (int i = 0; i < ingr_spheres.nbLevel - 1; i++){
					//get atom mapping
					cout << "gather mapping for level " << i << endl;
					ingr_spheres.LevelMappingStarts.push_back(ingr_spheres.LevelMapping.size());
					for (int j = 0; j < ingr_spheres.LevelCounts[i]; j++){
						int start;
						ifs.read(reinterpret_cast<char*>(&start), sizeof(int));
						int count;
						ifs.read(reinterpret_cast<char*>(&count), sizeof(int));
						ingr_spheres.LevelMapping.push_back(Vec2(start, count));
						cout << "atom_map " << i << " " << j << " start " << start << " count " << count << endl;
					}
					//check with previous
					for (int j = i - 1; j < ingr_spheres.nbLevel - 2; j++){
						if (j < 0) continue;
						if (i == j) continue;
						//get mapping
						ingr_spheres.LevelMappingStarts.push_back(ingr_spheres.LevelMapping.size());
						for (int k = 0; k < ingr_spheres.LevelCounts[i]; k++){
							int start;
							ifs.read(reinterpret_cast<char*>(&start), sizeof(int));
							int count;
							ifs.read(reinterpret_cast<char*>(&count), sizeof(int));
							ingr_spheres.LevelMapping.push_back(Vec2(start, count));
							cout <<"level_map " << i << " " << j << " " << k << " start " << start << " count " << count << endl;
						}
					}
				}
				//gather bindingsite info if any
				if (ingr_spheres.nBinding != 0){
					int total = 1;
					for (int i = 0; i < ingr_spheres.BindingStarts.size(); i+=2){
						int start = ingr_spheres.BindingStarts[i];
						int count = ingr_spheres.BindingStarts[i + 1];
						cout << "start " << start << " count " << count << endl;
						//int start = ingr_spheres.BindingStarts[i*ingr_spheres.nBinding + 0];
						for (int j = 0; j < count; j++){
							int tmp;
							ifs.read(reinterpret_cast<char*>(&tmp), sizeof(int));
							ingr_spheres.BindingSites.push_back(tmp);
						}
						total *= count;
					}
					cout << " total distance number " << total << endl;
					//gather the distance matrix, mapping ?
					for (int i = 0; i < total; i++){
						float tmp;
						ifs.read(reinterpret_cast<char*>(&tmp), sizeof(float));
						ingr_spheres.DistancesMatrix.push_back(tmp);
						cout << " i " << i << " d " << tmp << endl;
					}
				}
			}
			else {
				std::cout << "Error opening file";
			}
		}
		else {
			std::ifstream ifs(filename);
			if (ifs.is_open()) {
				std::cout << "file opened and closed";
				for (std::string line; std::getline(ifs, line);)   //read stream line by line
				{
					std::istringstream in(line);      //make a stream for the line itself
					float x, y, z;
					in >> x >> y >> z;       //now read the whitespace-separated floats
					ingr_spheres.LevelPoints.push_back(Vec3(x, y, z));
				}
				ingr_spheres.LevelCounts.push_back(ingr_spheres.LevelPoints.size());
			}
			else {
				std::cout << "Error opening file";
			}
		}
		return ingr_spheres;
	}
	
	Json::Value getJsonIngredientsSerialized(Json::Value  comp, string name) {
		//std::cout << "getIngredientFromCompartment " << comp["name"].asString() << endl;
		if (comp["IngredientGroups"].size() != 0) {
			//std::cout << "comp[IngredientGroups].size() " << comp["IngredientGroups"].size() << endl;
			Json::Value igroup = comp["IngredientGroups"][0];
			if (igroup["Ingredients"].size() != 0) {
				Json::Value ingredients = igroup["Ingredients"];
				//std::cout << "compartment should have n ingredients " << ingredients.size() << endl;
				Json::Value ingr_node_name;
				for (Json::ValueIterator itr = ingredients.begin(); itr != ingredients.end(); itr++) {
					ingr_node_name = *itr;
					string ingr_name = ingr_node_name["name"].asString();
					if (ingr_name == name)
						return ingr_node_name;
				}
			}
		}
		return 0;
	}

	Json::Value getIngredientsSerialized(string name) {
		Json::Value  root = book_json;
		Json::Value ingr_node_name;
		ingr_node_name = getJsonIngredientsSerialized(root, name);
		if (ingr_node_name != 0)
			return ingr_node_name;
		Json::Value comp = book_json["Compartments"];
		int i = 0;
		if (comp.size() != 0)
		{
			//std::cout << "find compartments " << comp.size() << endl;
			for (int i = 0; i < comp.size(); i++)
			{
				Json::Value comp_name = comp[i];
				//std::cout << "compartment should have several childs " << comp_name.size() << " " << comp_name["name"].asString() << endl;
				if (comp_name.size() == 0) continue;
				Json::Value  comp_childs = comp_name["Compartments"];
				for (int j = 0; j < comp_childs.size(); j++) {
					if (comp_childs[j]["name"].asString() == "surface") {
						ingr_node_name = getJsonIngredientsSerialized(comp_childs[j], name);
						if (ingr_node_name != 0)
							return ingr_node_name;
					}
					else if (comp_childs[i]["name"].asString() == "interior") {
						ingr_node_name = getJsonIngredientsSerialized(comp_childs[j], name);
						if (ingr_node_name != 0)
							return ingr_node_name;
					}
					else {
						ingr_node_name = getJsonIngredientsSerialized(comp_childs[j], name);
						if (ingr_node_name != 0)
							return ingr_node_name;
					}
				}
			}
		}
		return 0;
	}

	Json::Value getProteinNode(int ingrId) {
		//cout << "getProteinNode ingrId " << ingrId << endl;
		if (ingrId < 0 || ingrId >= pnames.size()) return NULL;
		//cout << "getFiberNode pnames_fiber " << pnames_fiber[fiberId] << endl;
		return proteins_nodes[ingrId];
	}

	Json::Value getJsonIngredients(Json::Value  ingr_nodes, string name){
		Json::Value ingr_node_name;
		for (Json::ValueIterator itr = ingr_nodes.begin(); itr != ingr_nodes.end(); itr++) {
			ingr_node_name = *itr;
			string ingr_name = ingr_node_name["name"].asString();
			if (ingr_name == name)
				return ingr_node_name;
		}
		return 0;
	}

	Json::Value getIngredients(string name){
		Json::Value  cyto = book_json["cytoplasme"];
		Json::Value ingr_node_name;
		if (cyto != 0)
		{
			Json::Value cyto_ingredients = cyto["ingredients"];
			if (cyto_ingredients != 0) {
				ingr_node_name = getJsonIngredients(cyto_ingredients, name);
				if (ingr_node_name != 0)
					return ingr_node_name;
			}

		}

		Json::Value  comp = book_json["compartments"];
		if (comp != 0)
		{
			for (Json::ValueIterator itr = comp.begin(); itr != comp.end(); itr++)
			{
				Json::Value comp_name = *itr;
				if (comp_name.size() == 0) continue;
				Json::Value comp_surface = comp_name["surface"];//at 0 ?
				if (comp_surface != 0)
				{
					Json::Value surf_ingredients = comp_surface["ingredients"];
					if (surf_ingredients != 0)
					{
						ingr_node_name = getJsonIngredients(surf_ingredients, name);
						if (ingr_node_name != 0)
							return ingr_node_name;
					}
				}
				if (comp_name.size() == 1) continue;
				Json::Value comp_interior = comp_name["interior"];//at 0 ?
				if (comp_interior != 0)
				{
					Json::Value int_ingredients = comp_interior["ingredients"];
					if (int_ingredients != 0) {
						ingr_node_name = getJsonIngredients(int_ingredients, name);
						if (ingr_node_name != 0)
							return ingr_node_name;
					}
				}
			}
		}
		return 0;
	}

	bool testFiberName(string aname){
		// return true;
		if (aname.find("DNA") != std::string::npos) return true;
		else if (aname.find("RNA") != std::string::npos) return true;
		else if (aname.find("PRO") != std::string::npos) return true;
		else if (aname.find("peptide") != std::string::npos) return true;
		else if (aname.find("lypoglycane") != std::string::npos) return true;
		else return false;
	}

	int getFiberId(string name) {
		for (int i = 0; i < pnames_fiber.size(); i++)
		{
			if (pnames_fiber[i] == name) return i;
		}
		return -1;
	}

	Json::Value getFiberNode(int fiberId) {
		//cout << "getFiberNode fiberId " << fiberId << endl;
		if (fiberId < 0 || fiberId >= pnames_fiber.size()) return NULL;
		//cout << "getFiberNode pnames_fiber " << pnames_fiber[fiberId] << endl;
		return pnames_fiber_nodes[fiberId];
	}

	int getFiberIdFromPointId(int part_id) {
		for (int i = 0; i < g_ropes.size(); i++) {
			int ropetype = g_ropes[i].ropeType;
			int startindex = g_ropes[i].mIndices[0];
			int count = g_ropes[i].mIndices.size();
			int endindex = g_ropes[i].mIndices[count - 1];
			if (part_id >= startindex && part_id <= endindex) {
				return ropetype;
			}
		}
		return -1;
	}

	int getPartnerId(string pname, Json::Value ingr_node) {
		if (ingr_node["partners_properties"].size() != 0) {
			for (int i = 0; i < ingr_node["partners_properties"].size(); i++) {
				cout << i << " " << pname << " " << ingr_node["partners_properties"][i]["partner_name"].asString() << endl;
				if (pname == ingr_node["partners_properties"][i]["partner_name"].asString())
				{
					return i;
				}
			}
		}
		return -1;
	}

	void checkPartnerAlongCurvePoints(int rope_id, Json::Value ingr_node)
	{
		Rope rope = g_ropes[rope_id];
		set<int> toRemove; //i,j or just i
		set<Vec2> toRemovePair;
		if (ingr_node["partners_name"].empty())
			return;
		int nPartner = ingr_node["partners_name"].size();
		int nPoints = rope.mIndices.size();//0 ? 
		cout << "ok " << nPoints << " " << nPartner << " " << ingr_node["name"].asString() << endl;
		std::vector<Json::Value> partner_nodes;
		std::vector<int> batchs_id;

		for (int j = 0; j < nPartner; j++)
		{
			//if (ingr_node["partners_name"][j] != "mpn529") continue;
			cout << "partner " << ingr_node["partners_name"][j].asString() << endl;
			Json::Value ingr_partner_node = getIngredients(ingr_node["partners_name"][j].asString());
			partner_nodes.push_back(ingr_partner_node);
			int batchid = getIngredientBatchId(ingr_node["partners_name"][j].asString());
			batchs_id.push_back(batchid);
			
			int nbMol = ingr_partner_node["nbMol"].asInt();
			float proba = (float)nbMol / (float)nPoints;
			cout << " nbmol " << ingr_partner_node["nbMol"].asInt() << " proba is " << proba << endl;
		}
		int koffset = 1;//actual persistence
		int previously_use_point = 0;
		for (int curveI = 3; curveI < nPoints - 3; curveI++)
		{
			if (nPartner <= 0)
				break;
			//pick a partner
			int idPartner = Rand() % nPartner;
			Json::Value ingr_partner_node = partner_nodes[idPartner];// getIngredients(ingr_node["partners_name"][idPartner].asString());
			int batchid = getIngredientBatchId(ingr_node["partners_name"][idPartner].asString());//batchs_id[idPartner];// 
			//proba to bind
			int nbMol = ingr_partner_node["nbMol"].asInt();
			float proba = (float)nbMol / (float)nPoints;
			float r = Random(0.0f, 1.0f);// Randf(0.0f, 1.0f);

			if (r > proba)
				continue;

			//partner number is total or per fiber??
			if (iBatches[batchid].nInstances >= nbMol){
				partner_nodes.erase(partner_nodes.begin() + idPartner);
				nPartner--;
				if (nPartner <= 0)
					break;
				continue;
			}
			//get position
			cout << "ok " << curveI << " " << proba << " " << r << " " << idPartner << " " << ingr_partner_node["nbMol"].asFloat() << " " << nPartner << " " << iBatches[batchid].nInstances << endl;
			int curve_ln = ingr_partner_node["properties"]["range"][0].asInt();//75?
			if (curve_ln > nPoints)
			{
				curve_ln = nPoints / 2;
			}
			int curve_var = int(Randf(-1.0f, 1.0f)* (float)ingr_partner_node["properties"]["range"][1].asInt()); //Rand(-1, 1) % ingrproperties["range"][1].asInt();
			if (curveI + curve_ln + curve_var >= nPoints)
			{
				curve_ln = (curveI + curve_ln + curve_var) - nPoints;
			}
			int middleIndex = curveI;
			if (!ingr_partner_node["properties"]["pairs"].empty())
			{
				middleIndex = curveI + (ingr_partner_node["properties"]["pairs"].size() / 2) - 1;
			}

			if (middleIndex - previously_use_point < 3)
				//make sure not two ingredient place to close
				continue;

			Vec3 to1 = g_buffers->positions[rope.mIndices[middleIndex + 1]] - g_buffers->positions[rope.mIndices[middleIndex]];
			Vec3 to2 = g_buffers->positions[rope.mIndices[middleIndex + 2]] - g_buffers->positions[rope.mIndices[middleIndex + 1]];
			Vec3 up = Cross(Normalize(to1), Vec3(0, 0, 1));

			Vec3 ingrpcpal = Vec3(ingr_partner_node["principalVector"][0].asFloat()*1.0f,
				ingr_partner_node["principalVector"][1].asFloat()*1.0f,
				ingr_partner_node["principalVector"][2].asFloat())*1.0f;

			Vec3 offsetPos = Vec3(ingr_partner_node["offset"][0].asFloat()*main_scale,
				ingr_partner_node["offset"][1].asFloat()*main_scale,
				ingr_partner_node["offset"][2].asFloat())*main_scale;

			//align pcpal to up
			Quat rot = AlignVec3s(up, ingrpcpal);
			//cout << "palce ingredient " << rot.x << " " << rot.y << " " << rot.z << endl;
			Vec3 pos = Vec3(g_buffers->positions[rope.mIndices[middleIndex]]);// +Rotate(rot, offsetPos*2.0f);// +to1*1.0f / 2.0f) + offset;
			pos += to1 / 2.0f;

			//cout << "before adding elem " << g_positions.size() << endl;

			//create instance
			int offset = createInstanceIngredient(batchid, pos, Quat(0, 0, 0, 1));

			previously_use_point = middleIndex;
			if (!ingr_partner_node["properties"]["pairs"].empty())
			{
				cout << " pairs binding " << endl;
				for (int i = 0; i < ingr_partner_node["properties"]["pairs"].size(); i++){
					for (int j = 0; j < ingr_partner_node["properties"]["pairs"][i].size(); j++){
						int beads = ingr_partner_node["properties"]["pairs"][i][j].asInt();
						cout << "add1 " << curveI + i << " " << rope.mIndices[curveI + i] << " " << offset + beads << endl;
						if (curveI + i < nPoints)
							CreateSpringInter(rope.mIndices[curveI + i], offset + beads, 1.0f, 0.0f, main_radius);
						//we should remove the 1-40 for theses
						//else 
						//	CreateSpring(rope.mIndices[(curveI + i) - rope.mIndices.size()], offset + beads, 1.0f, 0.0f, main_radius*2.0f);
					}
					//remove?
					//if (curveI + i - 40 >= 0) toRemovePair.insert(Vec2(curveI + i - 40, curveI + i));
					//if (curveI + i + 40 < rope.mIndices.size()) toRemovePair.insert(Vec2(curveI + i, curveI + i + 40));
				}
				cout << curveI << " from " << middleIndex - koffset << " to " << middleIndex + koffset << endl;
				//for (int ak = middleIndex - koffset; ak < middleIndex; ak++){
				//	if (ak + 40 < rope.mIndices.size()) toRemovePair.insert(Vec2(ak, ak+40));
				//}
				for (int ak = middleIndex + 1; ak <= middleIndex + 20; ak++){
					if (ak - 20 >= 0) toRemovePair.insert(Vec2(ak - 20, ak));
				}
			}
			else if (!ingr_partner_node["properties"]["beadsin"].empty())
			{
				cout << " beadsin binding " << endl;
				for (int i = 0; i < ingr_partner_node["properties"]["beadsin"].size(); i++)
				{
					int beadsin = ingr_partner_node["properties"]["beadsin"][i].asInt();
					CreateSpringInter(rope.mIndices[curveI], offset + beadsin, 1, 0.0f, main_radius*2.0f);
					if (curveI + 1 < nPoints)
					{
						CreateSpringInter(rope.mIndices[curveI + 1], offset + beadsin, 0.5f, 0.0f, main_radius*2.0f);
					}
					if (curveI - 1 > 0)
					{
						CreateSpringInter(rope.mIndices[curveI - 1], offset + beadsin, 0.5f, 0.0f, main_radius*2.0f);
					}
					cout << "addx " << curveI + 1 << " " << rope.mIndices[curveI] << " " << offset + beadsin << endl;
				}
				//if (curveI + 1 - 40 >= 0) toRemovePair.insert(Vec2(curveI + 1, curveI + 1 - 40));
				//if (curveI + 1 + 40 < rope.mIndices.size()) toRemovePair.insert(Vec2(curveI + 1 + 40, curveI + 1));
				//if (curveI - 1 - 40 >= 0) toRemovePair.insert(Vec2(curveI -1, curveI - 1 - 40));
				//if (curveI - 1 + 40 < rope.mIndices.size()) toRemovePair.insert(Vec2(curveI - 1 + 40, curveI - 1));

				for (int i = 0; i < ingr_partner_node["properties"]["beadsout"].size(); i++)
				{
					int beadsout = ingr_partner_node["properties"]["beadsout"][i].asInt();
					int add_id = curveI + curve_ln + curve_var;
					if (add_id >= nPoints)
					{
						add_id = add_id - nPoints;
					}
					if (add_id < nPoints)
					{
						CreateSpringInter(rope.mIndices[add_id], offset + beadsout, 1, 0.0f, main_radius*2.0f);
						cout << i << " " << beadsout << "addy " << add_id << " " << rope.mIndices[add_id] << " " << offset + beadsout << endl;
						//if (rope.mIndices[add_id] - 40 >= 0) toRemovePair.insert(Vec2(rope.mIndices[add_id], rope.mIndices[add_id] - 40));
						//if (rope.mIndices[add_id] + 40 < rope.mIndices.size()) toRemovePair.insert(Vec2(rope.mIndices[add_id] + 40, rope.mIndices[add_id]));
					}
					if (add_id + 1 < nPoints)
					{
						CreateSpringInter(rope.mIndices[add_id + 1], offset + beadsout, 0.5f, 0.0f, main_radius*2.0f);
						cout << i << " " << beadsout << "addy " << add_id + 1 << " " << rope.mIndices[add_id + 1] << " " << offset + beadsout << endl;
						//if (rope.mIndices[add_id + 1] - 40 >= 0) toRemovePair.insert(Vec2(rope.mIndices[add_id + 1], rope.mIndices[add_id + 1] - 40));
						//if (rope.mIndices[add_id + 1] + 40 < rope.mIndices.size()) toRemovePair.insert(Vec2(rope.mIndices[add_id + 1] + 40, rope.mIndices[add_id + 1]));
					}
					if (add_id - 1 > 0)
					{
						CreateSpringInter(rope.mIndices[add_id - 1], offset + beadsout, 0.5f, 0.0f, main_radius*2.0f);
						cout << i << " " << beadsout << "addy " << add_id - 1 << " " << rope.mIndices[add_id - 1] << " " << offset + beadsout << endl;
						//if (rope.mIndices[add_id - 1] - 40 >= 0) toRemovePair.insert(Vec2(rope.mIndices[add_id - 1], rope.mIndices[add_id - 1] - 40));
						//if (rope.mIndices[add_id - 1] + 40 < rope.mIndices.size()) toRemovePair.insert(Vec2(rope.mIndices[add_id - 1] + 40, rope.mIndices[add_id - 1]));
					}

				}
			}
			else {
				
			}
			//from curveI to curveI + pairs.size() +- 5 beads ?


			if (!ingr_partner_node["properties"]["st_ingr"].empty())
			{
				Json::Value ingr_fiber_spawn = getIngredients(ingr_partner_node["properties"]["st_ingr"].asString());
				int fiber_id = getFiberId(ingr_partner_node["properties"]["st_ingr"].asString());

				//local coordinate or beads ? nead to atacch to
				Vec3 stpt1 = Vec3(ingr_partner_node["properties"]["st_pt1"][0].asFloat()*1.0f,
					ingr_partner_node["properties"]["st_pt1"][1].asFloat()*1.0f,
					ingr_partner_node["properties"]["st_pt1"][2].asFloat())*1.0f;

				Vec3 stpt2 = Vec3(ingr_partner_node["properties"]["st_pt2"][0].asFloat()*1.0f,
					ingr_partner_node["properties"]["st_pt2"][1].asFloat()*1.0f,
					ingr_partner_node["properties"]["st_pt2"][2].asFloat())*1.0f;
				//need to start a new  rope of kind 

				if (!ingr_partner_node["properties"]["b_pt1"].empty()){
					int beadsid = ingr_partner_node["properties"]["b_pt1"][0].asInt();
					//attach the new rope starting beads to the above beads.
					//spawn event to grow a new fiber
					fiber_togrow.push_back(fiber_id);
					fiber_togrow_stpos.push_back(offset + beadsid);
					//get the length based on position on rope input ?
					if (!ingr_node["properties"]["lnSeg"].empty()){
						int N = rope.mIndices.size();//total nbbeads
						int lengthSoFar = -1;
						//closest starting indices
						float lnSeg = ingr_node["properties"]["lnSeg"].asFloat();
						if (lnSeg < 0) {
							//count from end
							float ratio = ingr_node["properties"]["ratioSeg"].asFloat();
							lengthSoFar = (int)((float)(N - middleIndex)*ratio); //size from end
						}
						else {
							int nSegments = (int)((lnSeg*main_scale) / g_params.radius);
							lengthSoFar = middleIndex%nSegments;
							//int startSegments = middleIndex - lengthSoFar;
						}
						//cout << "found " << nSegments << " " << middleIndex << " " << lengthSoFar << endl;
						fiber_togrow_length.push_back(lengthSoFar);
					}
					else
					{
						fiber_togrow_length.push_back(-1);
					}
					//rope.mIndices[middleIndex]
					//fiber_togrow_length.push_back();//nbBeads or actual length ?
					//Rope curve;
					//curve.phase = NvFlexMakePhase(iGroupCounter++, eNvFlexPhaseSelfCollide);
					//g_ropes.push_back(curve);//instance
					//mask.push_back(-1);
					//maks_protein.push_back(-(fiber_id + 1));//it ptype is 0 ?
					//maks_fiber.push_back(fiber_id);//proteinType
				}
			}
		}
		cout << "remove " << toRemovePair.size() << endl;
		RemoveSpringPair(toRemovePair);
	}

	void growOneCurveFromSurface(int pType, int compId, int vIndex, float uLength, float length){
		cout << "one curve " << compId - 1 << " " << comp_tri.size() << endl;
		Vec3 start = Vec3(comp_tri[compId - 1]->m_positions[vIndex].x, comp_tri[compId - 1]->m_positions[vIndex].y, comp_tri[compId - 1]->m_positions[vIndex].z);
		Vec3 normal = Vec3(comp_tri[compId - 1]->m_normals[vIndex].x, comp_tri[compId - 1]->m_normals[vIndex].y, comp_tri[compId - 1]->m_normals[vIndex].z);
		std::vector<Vec3> startend;
		startend.push_back(start);
		startend.push_back(start + (normal / Length(normal))*length);
		int curve_phase = NvFlexMakePhase(iGroupCounter++, eNvFlexPhaseSelfCollide);

		float* datac = reinterpret_cast<float*>(startend.data());
		int nbfloat = 6;
		Rope curve;
		int subdivid = (int)length / (main_radius*2.0);

		CreateRopeFromData(curve, //rope
			1.0f, // scale
			1.0f, //stifness
			datac, //data
			uLength, //length
			nbfloat,  //nfloat
			curve_phase,//phase
			0.0f,//spiral angle
			1.0f,//invmass
			0.0f,//give
			subdivid,//extend_nb
			true,//extend
			false,//close
			main_radius*2.0f);
		curve.phase = curve_phase;
		//we still need to fix start to surface, we could attach the particle using sdf, or mass ?
		// put huge infiinte ? mass
		g_buffers->positions[curve.mIndices[0]][3] = 0.0f; //first points attached or both first points
		g_ropes.push_back(curve);//instance

		mask.push_back(-1);
		maks_protein.push_back(-(pType + 1));//it ptype is 0 ?
		maks_fiber.push_back(pType);//proteinType
	}

	void placeOneFiberSurface(int fiber_id){
		Json::Value ingr_node_name = getIngredients(pnames_fiber[fiber_id]);
		int compId = iBatchesFiber[fiber_id].compId;
		cout << "build rope " << endl;
		int ncomp = comp_tri.size();
		int nVertices = comp_tri[compId - 1]->GetNumVertices();
		//nb of curve
		int nbCurve = ingr_node_name["nbMol"].asInt();
		//length of one curve
		float length = ingr_node_name["length"].asFloat()*main_scale;
		float uLength = ingr_node_name["uLength"].asFloat();
		//place xParticles in direction of the normal, let flex relax
		cout << nbCurve << " nbCurve of length " << length << endl;
		nbCurve = 200;
		for (int j = 0; j < nbCurve; j++){
			//pick a random vertices, start from here
			int vIndex = Rand() % nVertices;
			growOneCurveFromSurface(fiber_id, compId, vIndex, uLength, length);
		}
	}

	void placeOneFiber(int fiber_id)
	{
		int nfloat = 0;
		std::vector<char>  result_curve;
		Json::Value ingr_node_name = getIngredients(pnames_fiber[fiber_id]);
		int compId = iBatchesFiber[fiber_id].compId;
		cout << "build rope " << fiber_id << endl;
		int nbCurve = ingr_node_name["nbMol"].asInt();
		float L = ingr_node_name["length"].asFloat()*main_scale;//in angstrom 30000.0*main_scale;// 
		//L = 100000.0*main_scale;
		int npoints = (int)(L / g_params.radius);//or main_radius *2 
		float *data_curve;
		//starting configuration is either a circle, a line or a given file
		string mode = ingr_node_name["startingMode"].asString();
		//mode = "";
		if (pnames_fiber[fiber_id] == "DNA") {
			if (mode == "") 
				mode = "pdb_file";//"file"
		}
		
		bool close = ingr_node_name["closed"].asBool();//in angstrom
		float D = g_params.radius;//(g_params.mRadius/4.0f)*1000.0f;//could be the persitence length, which is 450bp 34A is one turn 10bp. Distance between spheres
		int subdivid = 1;// (int)((L / D) / (float)N);
		
		int persistence = 2;// (int)round(((500.0f*main_scale) / (main_radius*2.0f)) / 4.0f);
		bool extend = false;
		//mode = "line";
		//close = false;
		L = 7000.0f*3.4f*main_scale;
		if (mode == "circle")
		{

			float radius_circle = (L*50) / k2Pi; 
			npoints = (int)(L / (g_params.radius*2.0f));
			nfloat = npoints * 3;
			data_curve = new float[nfloat];
			int count = 0;
			//a go from 0 to k2Pi
			// circle 1 238 8 37.8789 1 0.1185
			cout << "circle " << nbCurve << " " << L << " " << npoints << " " << radius_circle << " " << close << " " << main_radius << " " << g_params.radius << " " << (float)npoints  * (g_params.radius*2.0f) * 1.0f/main_scale<< endl;
			float a = 0.0f;
			for (int i = 0; i < npoints; i++){
				data_curve[count] = radius_circle * cos(a);
				data_curve[count + 1] = radius_circle * sin(a) + radius_circle + D*20;//offset 
				data_curve[count + 2] = Randf(-D*10.0f, D*10.0f);// 0.0f;//randomize z
				count += 3;
				a += k2Pi / (float)npoints;
			}
			extend = false;
		}
		else if (mode == "file")
		{
			//use thebinary file...dna
			cout << "build rope " << datapath<< " " << "tps_path.bin" << endl;
			ifstream ifs2(datapath + "tps_path.bin", ios::binary | ios::ate);
			ifstream::pos_type pos2 = ifs2.tellg();

			result_curve.resize(pos2);

			ifs2.seekg(0, ios::beg);
			ifs2.read(&result_curve[0], pos2);

			printf("read bytes %zu\n", result_curve.size());
			//now convert to float and place a particle there
			int step = 3;

			nfloat = result_curve.size() / sizeof(float);
			printf("which are %i float %i %zu\n", nfloat, nfloat / step, result_curve.size() / sizeof(float));

			int N = nfloat / step;
			float PL = 500.0f;//thats 25 spheres
			subdivid = (int)((L / D) / (float)N);
			data_curve = reinterpret_cast<float*>(result_curve.data());
			//g_numExtraParticles = int(L/D);
			printf(" use N=%i L=%f D=%f total excpected=%i\n", N, L, D, npoints);
			printf("create Rope from data with sibdivd %i and D %f and %i nExtraParticle\n", subdivid, D, g_numExtraParticles);
			extend = true;
		}
		else if (mode == "pdb_file")
		{
			//use the pdb file from David Goodsell
			cout << "build rope " << datapath << " " << "myco_relax_jpc.pdb" << endl;
			std::ifstream ifs(datapath + "myco_relax_jpc.pdb");
			std::vector<float> tmp_data;
			if (ifs.is_open()) {
				std::cout << "file opened and closed";
				for (std::string line; std::getline(ifs, line);)   //read stream line by line
				{
					//MODEL        1
					//ATOM      1  CA  GLY A   2      45.609 - 4.405   7.659
					if (line.substr(0, 4) == "ATOM") 
					{
						stringstream in(line);
						string header;
						int atnum;
						string atname;
						string resname;
						string chname;
						int resnum;
						float x, y, z;
						in >> header;
						in >> atnum;
						in >> atname;
						in >> resname;
						in >> chname;
						in >> resnum;
						in >> x;
						in >> y;
						in >> z;
						tmp_data.push_back(x * 10.0f);
						tmp_data.push_back(y * 10.0f);
						tmp_data.push_back(z * 10.0f);
						//cout << header << " " <<atnum << " " << atname << " " << resname << " " << chname << " " << resnum << " " << x << " " << y << " " << z << endl;
					}
				}
				nfloat = tmp_data.size();
			}
			else {
				std::cout << "Error opening file";
			}
			
			int step = 3;
			int N = nfloat / step;
			npoints = N;
			float PL = 50.0f;//thats 25 spheres
			subdivid = 1;// (int)((L / D) / (float)N);
			//data_curve = reinterpret_cast<float*>(tmp_data.data());
			//g_numExtraParticles = int(L/D);
			printf(" use N=%i L=%f D=%f total excpected=%i\n", N, L, D, npoints);
			printf("create Rope from data with sibdivd %i and D %f and %i nExtraParticle\n", subdivid, D, g_numExtraParticles);
			extend = false;
			int rope_phase = NvFlexMakePhase(iGroupCounter++, eNvFlexPhaseSelfCollide);

			Rope curve;
			CreateRopeFromData(curve, //rope
				main_scale, // scale
				1.0f, //stifness
				tmp_data, //data
				L, //length
				nfloat,  //nfloat
				rope_phase,//phase
				0.0f,//spiral angle
				1.0f,//invmass
				0.0f,//give
				subdivid,//extend_nb
				extend,//extend
				close,//close
				main_radius*2.0f,
				persistence);
			curve.persistence = persistence;
			printf("create Rope from data OK with %zu points\n", curve.mIndices.size());
			g_ropes.push_back(curve);//instance
			mask.push_back(-1);
			maks_protein.push_back(-(fiber_id + 1));//it ptype is 0 ?
			maks_fiber.push_back(fiber_id);//proteinType
			printf("create Rope from data OK with %i  %i\n", -(fiber_id + 1), fiber_id);
			return;
		}
		else if (mode == "line"){
			//build a straight line
			npoints = (int)(L / (g_params.radius*2.0));
			nfloat = npoints * 3;
			data_curve = new float[nfloat];
			int count = 0;
			for (int i = 0; i < npoints; i++){
				data_curve[count] = 100.0f*g_params.radius*(float)i;//x
				data_curve[count + 1] = 0.0f;//offset 
				data_curve[count + 2] = 0.0f;// 0.0f;//randomize z
				count += 3;
			}
			extend = false;
		}
		else
		{
			//line
			//from a random starting point generate a line with Npoint
			Vec3 start = Vec3(Randf(0.0f, 1.0f), Randf(0.0f, 1.0f), Randf(0.0f, 1.0f));
			std::vector<Vec3> startend;
			startend.push_back(start);
			startend.push_back(start + Normalize(start)*L);
			data_curve = reinterpret_cast<float*>(startend.data());
			subdivid = npoints;
			nfloat = 6;
			extend = true;
		}
		int rope_phase = NvFlexMakePhase(iGroupCounter++, eNvFlexPhaseSelfCollide);

		Rope curve;
		CreateRopeFromData(curve, //rope
			main_scale, // scale
			1.0f, //stifness
			data_curve, //data
			L, //length
			nfloat,  //nfloat
			rope_phase,//phase
			0.0f,//spiral angle
			1.0f,//invmass
			0.0f,//give
			subdivid,//extend_nb
			extend,//extend
			close,//close
			main_radius*2.0f,
			persistence);
		curve.persistence = persistence;
		printf("create Rope from data OK with %zu points\n", curve.mIndices.size());
		g_ropes.push_back(curve);//instance
		mask.push_back(-1);
		maks_protein.push_back(-(fiber_id + 1));//it ptype is 0 ?
		maks_fiber.push_back(fiber_id);//proteinType
		printf("create Rope from data OK with %i  %i\n", -(fiber_id + 1), fiber_id);
	}

	void placeFibers(bool check_partners){
		//loop over fiber ingredients
		for (int i = 0; i<pnames_fiber.size(); i++)
		{
			if (pnames_fiber[i] != "DNA") continue;
			cout << "place fibers " << pnames_fiber[i] << " " << i << endl;
			Json::Value ingr_node_name = getIngredients(pnames_fiber[i]);
			//get partner if any
			
			IngredientPartner partners;
			partners.nPartner = ingr_node_name["partners_name"].size();
			if (partners.nPartner > 0)
			{
				cout << "place fibers " << pnames_fiber[i] << " with " << partners.nPartner << "partners" << endl;
				for (int j = 0; j < partners.nPartner; j++)
				{
					Json::Value ingr_partner_node = getIngredients(ingr_node_name["partners_name"][j].asString());
					partners.partner_nodes.push_back(ingr_partner_node);
					int batchid = getIngredientBatchId(ingr_node_name["partners_name"][j].asString());
					partners.batchs_id.push_back(batchid);
				}
			}
			iPartnersFibers.push_back(partners);
			 
			int compId = iBatchesFiber[i].compId;
			//compid should gave us the mesh
			if (compId > 0) {
				 placeOneFiberSurface(i);
			}
			else {
				placeOneFiber(i);
				//break;
			}
		}
		if (!check_partners) return;
		for (int i = 0; i < g_ropes.size(); i++){
			Json::Value ingr_node_name = getIngredients(pnames_fiber[maks_fiber[i]]);
			checkPartnerAlongCurvePoints(i, ingr_node_name);
		}
		int rope_id = g_ropes.size();
		int rope_nb = fiber_togrow.size();
		
		fiberToGrow();

		cout << rope_id << " rope id " << rope_nb << " rope_nb " << g_ropes.size() << endl;

		//checkpartner on theses new ropes
		for (int i = rope_id; i < rope_id + rope_nb; i++){
			Json::Value ingr_node_name = getIngredients(pnames_fiber[maks_fiber[i]]);
			int nPartner = ingr_node_name["partners_name"].size();
			int nPoints = g_ropes[i].mIndices.size();//0 ? 
			cout << "check " << nPoints << " " << nPartner << " " << ingr_node_name["name"].asString() << endl;
			checkPartnerAlongCurvePoints(i, ingr_node_name);
		}
		//fiber To grow ?
		fiberToGrow();
	}

	void duplicateMeshMembrane(string path)
	{
		// membrane thickness is 47A
		float mb = g_mb / 2.0f * main_scale;
		cout << "mb " << mb << endl;
		Mesh* upshape = ImportMesh(GetFilePathByPlatform((path).c_str()).c_str());
		// Calculate the original center of the mesh
		Vector3 minExtents, maxExtents;
		upshape->GetBounds(minExtents, maxExtents);
		Vector3 originalCenter = (minExtents + maxExtents) * 0.5f;

		upshape->Transform(ScaleMatrix(main_scale));
		upshape->Transform(ScaleMatrix(1 + mb));

		// Re-calculate the center after scaling (assuming the mesh might be modified by Transform)
		Vector3 newMinExtents, newMaxExtents;
		upshape->GetBounds(newMinExtents, newMaxExtents);
		Vector3 newCenter = (newMinExtents + newMaxExtents) * 0.5f;

		// Calculate the translation required to re-center the mesh
		Vector3 centerTranslation = originalCenter - newCenter;
		upshape->Transform(TranslationMatrix(Point3(centerTranslation)));

		if (g_winding) {
			for (int i = 0; i < int(upshape->GetNumFaces()); ++i)
				swap(upshape->m_indices[i * 3 + 0], upshape->m_indices[i * 3 + 1]);
		}
		upshape->CalculateNormals();
		NvFlexTriangleMeshId upmesh = CreateTriangleMesh(upshape);
		AddTriangleMesh(upmesh, Vec3(), Quat(), 1.0f);
		comp_mesh.push_back(upmesh);
		comp_tri.push_back(upshape);

		Mesh* downshape = ImportMesh(GetFilePathByPlatform((path).c_str()).c_str());
		downshape->Transform(ScaleMatrix(main_scale));
		downshape->Transform(ScaleMatrix(1 - mb));
		// Re-calculate the center after scaling (assuming the mesh might be modified by Transform)
		downshape->GetBounds(newMinExtents, newMaxExtents);
		newCenter = (newMinExtents + newMaxExtents) * 0.5f;

		// Calculate the translation required to re-center the mesh
		centerTranslation = originalCenter - newCenter;
		downshape->Transform(TranslationMatrix(Point3(centerTranslation)));

		if (!g_winding) {
			for (int i = 0; i < int(downshape->GetNumFaces()); ++i)
				swap(downshape->m_indices[i * 3 + 0], downshape->m_indices[i * 3 + 1]);
		}
		downshape->CalculateNormals();
		NvFlexTriangleMeshId downmesh = CreateTriangleMesh(downshape);
		AddTriangleMesh(downmesh, Vec3(), Quat(), 1.0f);
		comp_mesh.push_back(downmesh);
		comp_tri.push_back(downshape);
	}

	void duplicateMeshMembrane(int nv, int nf, Vec3* vertices, Vec3* normals, int* faces)
	{
		// membrane thickness is 47A
		float mb = g_mb / 2.0f * main_scale;
		cout << "mb " << mb << endl;
		Mesh* upshape = createMesh(nv, nf, vertices, normals, faces);
		upshape->Transform(ScaleMatrix(main_scale));
		Vector3 minExtents, maxExtents;
		upshape->GetBounds(minExtents, maxExtents);
		Vector3 originalCenter = (minExtents + maxExtents) * 0.5f;

		upshape->Transform(ScaleMatrix(1 + mb));
		// Re-calculate the center after scaling (assuming the mesh might be modified by Transform)
		Vector3 newMinExtents, newMaxExtents;
		upshape->GetBounds(newMinExtents, newMaxExtents);
		Vector3 newCenter = (newMinExtents + newMaxExtents) * 0.5f;

		// Calculate the translation required to re-center the mesh
		Vector3 centerTranslation = originalCenter - newCenter;
		upshape->Transform(TranslationMatrix(Point3(centerTranslation)));

		if (g_winding) {
			for (int i = 0; i < int(upshape->GetNumFaces()); ++i)
				swap(upshape->m_indices[i * 3 + 0], upshape->m_indices[i * 3 + 1]);
		};
		upshape->CalculateNormals();
		NvFlexTriangleMeshId upmesh = CreateTriangleMesh(upshape);
		AddTriangleMesh(upmesh, Vec3(), Quat(), 1.0f);
		comp_mesh.push_back(upmesh);
		comp_tri.push_back(upshape);
		Mesh* downshape = createMesh(nv, nf, vertices, normals, faces);
		downshape->Transform(ScaleMatrix(main_scale));
		downshape->Transform(ScaleMatrix(1 - mb));
		
		// Re-calculate the center after scaling (assuming the mesh might be modified by Transform)
		downshape->GetBounds(newMinExtents, newMaxExtents);
		newCenter = (newMinExtents + newMaxExtents) * 0.5f;

		// Calculate the translation required to re-center the mesh
		centerTranslation = originalCenter - newCenter;
		downshape->Transform(TranslationMatrix(Point3(centerTranslation)));

		if (!g_winding) {
			for (int i = 0; i < int(downshape->GetNumFaces()); ++i)
				swap(downshape->m_indices[i * 3 + 0], downshape->m_indices[i * 3 + 1]);
		}
		downshape->CalculateNormals();
		NvFlexTriangleMeshId downmesh = CreateTriangleMesh(downshape);
		AddTriangleMesh(downmesh, Vec3(), Quat(), 1.0f);
		comp_mesh.push_back(downmesh);
		comp_tri.push_back(downshape);
	}


	void duplicateMeshMembraneNode(Json::Value node)
	{
		Mesh* upshape = getMeshFromJsonNode(node);
		upshape->Transform(ScaleMatrix(main_scale));
		upshape->Transform(ScaleMatrix(1 + g_params.radius / 100));
		upshape->CalculateNormals();
		NvFlexTriangleMeshId upmesh = CreateTriangleMesh(upshape);
		AddTriangleMesh(upmesh, Vec3(), Quat(), 1.0f);
		comp_mesh.push_back(upmesh);
		Mesh* downshape = getMeshFromJsonNode(node);
		downshape->Transform(ScaleMatrix(main_scale));
		downshape->Transform(ScaleMatrix(1 - g_params.radius / 100));
		for (int i = 0; i < int(downshape->GetNumFaces()); ++i)
			swap(downshape->m_indices[i * 3 + 0], downshape->m_indices[i * 3 + 1]);
		downshape->CalculateNormals();
		NvFlexTriangleMeshId downmesh = CreateTriangleMesh(downshape);
		AddTriangleMesh(downmesh, Vec3(), Quat(), 1.0f);
		comp_mesh.push_back(downmesh);
	}

	Mesh* createMesh(int nv, int nf, Vec3* vertices, Vec3* normals, int* faces){
		int numVertices = nv;
		int numFaces = nf;
		Mesh* mesh = new Mesh;
		mesh->m_positions.resize(numVertices);
		mesh->m_normals.resize(numVertices);
		mesh->m_colours.resize(numVertices, Colour(1.0f, 1.0f, 1.0f, 1.0f));
		mesh->m_indices.reserve(nf);
		for (uint32_t v = 0; v < numVertices; ++v)
		{
			mesh->m_positions[v] = Point3(-vertices[v].x/0.065f, vertices[v].y / 0.065f, vertices[v].z / 0.065f);
			mesh->m_normals[v] = Vector3(-normals[v].x, normals[v].y, normals[v].z);
		}
		for (uint32_t f = 0; f < numFaces; ++f)
		{
			mesh->m_indices.push_back(faces[f]);
		}
		//for (uint32_t i = 0; i < numVertices; ++i)
		//{
		//	mesh->m_normals[i] = SafeNormalize(mesh->m_normals[i], Vector3(0.0f, 1.0f, 0.0f));
		//}
		return mesh;
	}

	Mesh* getMeshFromJsonNode(Json::Value node) {
		int numVertices = node["verts"].size() / 3;
		int numFaces = node["faces"].size() / 3;
		Mesh* mesh = new Mesh;
		mesh->m_positions.resize(numVertices);
		mesh->m_normals.resize(numVertices);
		mesh->m_colours.resize(numVertices, Colour(1.0f, 1.0f, 1.0f, 1.0f));
		mesh->m_indices.reserve(numFaces * 3);
		for (uint32_t v = 0; v < numVertices; ++v)
		{
			mesh->m_positions[v] = Point3(node["verts"][v * 3].asFloat(),
				node["verts"][v * 3 + 1].asFloat(),
				node["verts"][v * 3 + 2].asFloat());
			mesh->m_normals[v] = Vector3(node["normals"][v * 3].asFloat(),
				node["normals"][v * 3 + 1].asFloat(),
				node["normals"][v * 3 + 2].asFloat());
			//mesh->m_normals[v] = Vector3(0.0f, 0.0f, 0.0f);
		}
		for (uint32_t f = 0; f < numFaces; ++f)
		{
			mesh->m_indices.push_back(node["faces"][f * 3].asUInt());
			mesh->m_indices.push_back(node["faces"][f * 3 + 1].asUInt());
			mesh->m_indices.push_back(node["faces"][f * 3 + 2].asUInt());
			/*
			// calculate vertex normals as we go
			Point3& v0 = mesh->m_positions[indices[0]];
			Point3& v1 = mesh->m_positions[indices[1]];
			Point3& v2 = mesh->m_positions[indices[2]];

			Vector3 n = SafeNormalize(Cross(v1-v0, v2-v0), Vector3(0.0f, 1.0f, 0.0f));

			for (uint32_t i=0; i < numIndices; ++i)
			{
				mesh->m_normals[indices[i]] += n;
			}
			*/
		}
		for (uint32_t i = 0; i < numVertices; ++i)
		{
			mesh->m_normals[i] = SafeNormalize(mesh->m_normals[i], Vector3(0.0f, 1.0f, 0.0f));
		}
		return mesh;
	}

	void compartmentMesh(Json::Value comp){
		Json::Value comp_geom_type = comp["geom_type"];//at 0 ?
		Json::Value comp_geom = comp["geom"];//at 0 ?
		Json::Value comp_geom_filename = comp["filename"];//at 0 ?
		Json::Value comp_name = comp["name"];//at 0 ?
		// if (comp_geom == 0) return;
		if (comp_geom_type.asString() == "raw") {
			Mesh* mesh = getMeshFromJsonNode(comp_geom);
			mesh->Transform(ScaleMatrix(main_scale));
			//for (int i = 0; i < int(shape->GetNumFaces()); ++i)
			//	swap(shape->m_indices[i * 3 + 0], shape->m_indices[i * 3 + 1]);
			//duplicate and inverse?
			duplicateMeshMembraneNode(comp_geom);
			Vec3 lower = Vec3(0.0f);
			Vec3 upper = Vec3(1.0f);
			mesh->GetBounds(lower, upper);
			comp_tri.push_back(mesh);
			//sdfdata = CreateSDFFromMesh(GetFilePathByPlatform((geompath + comp_name.asString()).c_str()).c_str(), mesh, dim);
		}
		if (comp_geom_type.asString() == "sphere") {
			comp_radius = 1600.0f*main_scale;
			Mesh* shape = CreateSphere(20, 20, 1600.0f*main_scale);//1512.38230769f
			Vec3 lower = Vec3(0.0f);
			Vec3 upper = Vec3(1.0f);
			shape->GetBounds(lower, upper);
			comp_tri.push_back(shape);
			//sdfdata = CreateSDFFromMesh(GetFilePathByPlatform((geompath + comp_name.asString()).c_str()).c_str(), shape, dim);

			Mesh* upshape = CreateSphere(20, 20, 1448.4f*main_scale);
			upshape->Transform(ScaleMatrix(1 + g_params.radius / 100));
			upshape->CalculateNormals();
			NvFlexTriangleMeshId upmesh = CreateTriangleMesh(upshape);
			AddTriangleMesh(upmesh, Vec3(), Quat(), 1.0f);
			comp_mesh.push_back(upmesh);
			Mesh* downshape = CreateSphere(20, 20, 1448.4f*main_scale);
			downshape->Transform(ScaleMatrix(1 - g_params.radius / 100));
			for (int i = 0; i < int(downshape->GetNumFaces()); ++i)
				swap(downshape->m_indices[i * 3 + 0], downshape->m_indices[i * 3 + 1]);
			downshape->CalculateNormals();
			NvFlexTriangleMeshId downmesh = CreateTriangleMesh(downshape);
			AddTriangleMesh(downmesh, Vec3(), Quat(), 1.0f);
			comp_mesh.push_back(downmesh);
		}
		if (comp_geom_type.asString() == "mb" && (comp_geom_filename.asString()== "" || comp_geom_filename.asString() == "None")) {
			Json::Value mbs = comp["mb"];
			comp_radius = 1600.0f*main_scale;
			float radius = mbs["radii"][0].asFloat();
			std::cout << "comp_geom_type mb  : " << radius << endl;
			//take the first sphere for now
			float x = mbs["positions"][0].asFloat()*main_scale;
			float y = mbs["positions"][1].asFloat()*main_scale;
			float z = mbs["positions"][2].asFloat()*main_scale;
			std::cout << "com pos " << x << " " << y << " " << z << endl;
			Matrix44 xform = TranslationMatrix(Point3(x,y,z));
			Mesh* shape = CreateSphere(20, 20, radius*main_scale);//1512.38230769f
			comp_tri.push_back(shape);
			//sdfdata = CreateSDFFromMesh(GetFilePathByPlatform((geompath + comp_name.asString()).c_str()).c_str(), shape, dim);
			Mesh* upshape = CreateSphere(20, 20, radius*main_scale);
			//upshape->Transform(xform);
			upshape->Transform(xform);//*ScaleMatrix(1 + g_params.radius / 100)
			upshape->CalculateNormals();
			NvFlexTriangleMeshId upmesh = CreateTriangleMesh(upshape);
			AddTriangleMesh(upmesh, Vec3(), Quat(), 1.0f);
			comp_mesh.push_back(upmesh);
			Mesh* downshape = CreateSphere(20, 20, radius*main_scale);
			//downshape->Transform(xform);
			downshape->Transform(xform);//*ScaleMatrix(1 - g_params.radius / 100)
			for (int i = 0; i < int(downshape->GetNumFaces()); ++i)
				swap(downshape->m_indices[i * 3 + 0], downshape->m_indices[i * 3 + 1]);
			downshape->CalculateNormals();
			NvFlexTriangleMeshId downmesh = CreateTriangleMesh(downshape);
			AddTriangleMesh(downmesh, Vec3(), Quat(), 1.0f);
			comp_mesh.push_back(downmesh);
		}
		if (comp_geom_type.asString() == "file" && comp_geom_filename.asString() != "" && comp_geom_filename.asString() != "None") {
			string path = "";
			if (comp_geom != 0) path = comp_geom.asString();//"mpn_surface_center.obj";// comp["geom"].asString();// "MMycoideHD.dae"; //comp["geom"].asString();
			if (comp_geom_filename != 0) path = comp_geom_filename.asString();
			std::cout << "path  : " << path << " compId " << comp_geom_filename << " " << comp_geom << endl;
			if (path == "") return;
			string objpath = string(path.c_str(), strrchr(path.c_str(), '.')) + ".obj";
			Mesh* shape = ImportMesh(GetFilePathByPlatform((geompath + objpath).c_str()).c_str());
			shape->Transform(ScaleMatrix(main_scale));
			for (int i = 0; i < int(shape->GetNumFaces()); ++i)
				swap(shape->m_indices[i * 3 + 0], shape->m_indices[i * 3 + 1]);
			duplicateMeshMembrane(geompath + objpath);
			/*
			// invert box faces
			for (int i = 0; i < int(shape->GetNumFaces()); ++i)
				swap(shape->m_indices[i * 3 + 0], shape->m_indices[i * 3 + 1]);
			shape->CalculateNormals();
			NvFlexTriangleMeshId mesh = CreateTriangleMesh(shape);
			AddTriangleMesh(mesh, Vec3(), Quat(), 1.0f);
			comp_mesh.push_back(mesh);
			*/
			Vec3 lower = Vec3(0.0f);
			Vec3 upper = Vec3(1.0f);
			shape->GetBounds(lower, upper);
			comp_tri.push_back(shape);
			//get the sdf

			//sdfdata = CreateSDFFromMesh(GetFilePathByPlatform((geompath + objpath).c_str()).c_str(), shape, dim);
			//add the sdf ? 
			//NvFlexDistanceFieldId sdf = CreateSDF(GetFilePathByPlatform((geompath + objpath).c_str()).c_str(), dim, 0.0f, 0.0f,-1.0f);
			//inverse field
			//AddSDF(sdf, Vec3(), Quat(), main_scale);
		}

	}

	void OneCompartmentMesh(bool outside = true) {
		//comp_radius = 1600.0f*main_scale;
		//take the first sphere for now
		std::cout << "OneCompartmentMesh  radius " << comp_radius << " scale " << main_scale << endl;
		Mesh* shape = CreateSphere(20, 20, comp_radius);//1512.38230769f
		comp_tri.push_back(shape);
		//sdfdata = CreateSDFFromMesh(GetFilePathByPlatform((geompath + "comp").c_str()).c_str(), shape, dim);
		Mesh* upshape = CreateSphere(20, 20, comp_radius);
		upshape->Transform(ScaleMatrix(1 + g_params.radius / 100));
		upshape->CalculateNormals();
		NvFlexTriangleMeshId upmesh = CreateTriangleMesh(upshape);
		AddTriangleMesh(upmesh, Vec3(), Quat(), 1.0f);
		comp_mesh.push_back(upmesh);
		Mesh* downshape = CreateSphere(20, 20, comp_radius);
		downshape->Transform(ScaleMatrix(1 - g_params.radius / 100));
		for (int i = 0; i < int(downshape->GetNumFaces()); ++i)
			swap(downshape->m_indices[i * 3 + 0], downshape->m_indices[i * 3 + 1]);
		downshape->CalculateNormals();
		NvFlexTriangleMeshId downmesh = CreateTriangleMesh(downshape);
		AddTriangleMesh(downmesh, Vec3(), Quat(), 1.0f);
		comp_mesh.push_back(downmesh);
	}

	int  getPcpalAxis(Vec3 pcpalVec){
		float m = maxf(maxf(fabsf(pcpalVec.x), fabsf(pcpalVec.y)), fabsf(pcpalVec.z));
		if (m == pcpalVec.x) return 0;
		else if (m == pcpalVec.y) return 1;
		else if (m == pcpalVec.z) return 2;
		return -1;
	}

	std::vector<Vec3> filterSurfaceBeads(Vec3 offset, Vec3 pcpalVec, std::vector<Vec3> points)
	{
		//skip all beads distance > beads_radius from offset
		int axis = getPcpalAxis(pcpalVec);
		cout << "pcpal Axis " << pcpalVec.x << " " << pcpalVec.y << " " << pcpalVec.z << " " << axis << endl;
		std::vector<Vec3> fpoints;
		for (int i = 0; i < points.size(); i++) {
			Vec3 tp = points[i] + offset;
			cout << " D " << tp.x << " " << tp.y << " " << tp.z << " " << fabsf(tp[axis]) << " " << g_params.radius  << endl;
			if (fabsf(tp[axis]) <= g_params.radius) {
				//offset to avoid the membrane
				float sign = 1;
				if (tp[axis] < 0)
					sign = -1;
				Vec3 p = Vec3(points[i].x, points[i].y, points[i].z);
				p[axis] = (p[axis] - tp[axis] )+ ((g_params.radius/1.5f))*sign;
				fpoints.push_back(p);
			}
			else 
				fpoints.push_back(points[i]);
		}
		return fpoints;
	}

	void IngredientToFixSDF(Json::Value ingr_node, Vec3 ipos){
		NvFlexDistanceFieldId sdf = CreateSDF(GetFilePathByPlatform("../../data/sphere.ply").c_str(), 32, 0.0f, 0.0f);
		Mesh* sphere_shape = CreateSphere(20, 20, 1.0);
		// sphere_shape->Transform(ScaleMatrix(1 + g_params.radius / 100));
		sphere_shape->CalculateNormals();
		NvFlexTriangleMeshId upmesh = CreateTriangleMesh(sphere_shape);
		int p = lodproxy_to_use;
		Json::Value jsonpos = ingr_node["positions"][p]["coords"];
		Json::Value jsonrad = ingr_node["radii_lod"][p]["radii"];
		for (int i = 0; i < jsonpos.size() / 3; i++) {
			float r = jsonrad[i].asFloat() * main_scale;
			Vec3 beadp = Vec3(-jsonpos[i * 3].asFloat(),
				jsonpos[i * 3 + 1].asFloat(),
				jsonpos[i * 3 + 2].asFloat()) * main_scale;
			if (isnan(beadp.x) || isnan(beadp.y) || isnan(beadp.z)) {
				beadp = Vec3(std::stof(jsonpos[i * 3].asString()),
					std::stof(jsonpos[i * 3 + 1].asString()),
					std::stof(jsonpos[i * 3 + 2].asString()));
			}
	
			//sdf, position, rotation, width/scale
			//float radius = 2000.0f * main_scale;
			// AddSDF(sdf, beadp + ipos + Vec3(-r/ 2.0f, -r/ 2.0f, -r/ 2.0f), Quat(), r);
			AddTriangleMesh(upmesh, beadp + ipos, Quat(), r);
			comp_mesh.push_back(upmesh);
		}
	}

	//create a flex asset from the ingredient data
	void ingredientToFlexAsset(Json::Value ingr_node, int compId, bool rb = false)
	{
		//need to transform with pcp and offset. using Quat.Rotate - Vec3 rpos = Rotate(rotation, localPos);
		string proxyname;
		Json::Value ingr_source = ingr_node["source"];
		string ingr_source_pdb = ingr_source["pdb"].asString();

		proxyname = datapath + ingr_source_pdb + "_cl.txt";
		std::cout << "ingr proxy  : " << proxyname << endl;

		Vec3 pcpalVector = Vec3(ingr_node["principalVector"][0].asFloat(),
			ingr_node["principalVector"][1].asFloat(),
			ingr_node["principalVector"][2].asFloat());
		if (pcpalVector == Vec3(0, 0, 0)) pcpalVector = Vec3(0, 0, 1);
		Json::Value  offsetnode = ingr_source["transform"]["offset"];
		Vec3 offset = Vec3(0.0f, 0.0f, 0.0f);//should be ingr_node["source"]["transform"]["offset"] if exist
		if (offsetnode != 0) {
			offset = Vec3(offsetnode[0].asFloat(),
				offsetnode[1].asFloat(),
				offsetnode[2].asFloat());
			//if surface use -offset
			if (compId <= 0) {
				offset = -offset;
			}
		}
		std::cout << "ingr pcpalVector  : " << pcpalVector.x << " " << pcpalVector.y << " " << pcpalVector.z << endl;
		std::cout << "ingr offset  : " << offset.x << " " << offset.y << " " << offset.z << endl;
		Json::Value m = ingr_node["meshFile"];
		Json::Value mesh_type = ingr_node["meshType"];
		string meshFile;
		Mesh* mesh;
		bool use_mesh = false;
		if (m != NULL)
		{
			//meshFile can also vert+face
			if (ingr_node["meshFile"].isString())
			{
				meshFile = ingr_node["meshFile"].asString();
				//std::cout << "meshFile  : " << meshFile << endl;

				if ((!meshFile.empty()) && (use_instances_mesh)) {
					meshFile.replace(meshFile.size() - 3, 3, "obj"); // (6) REMOVE DAE TO OBJ
					std::cout << "meshFile after : " << meshFile << endl;
					mesh = ImportMesh(GetFilePathByPlatform((geompath + meshFile).c_str()).c_str());
					//mesh = ImportMesh(GetFilePathByPlatform(meshFile.c_str()).c_str());
					mesh->Transform(ScaleMatrix(main_scale));
					for (int i = 0; i < int(mesh->GetNumFaces()); ++i)
						swap(mesh->m_indices[i * 3 + 0], mesh->m_indices[i * 3 + 1]);
					mesh->CalculateNormals(); //??
					use_mesh = true;
				}
			}
			else {
				//mesh define directly with "verts","faces" and "normals" from ,"meshType"=="raw"
				mesh = getMeshFromJsonNode(ingr_node["meshFile"]);
				mesh->Transform(ScaleMatrix(main_scale));
				//for (int i = 0; i < int(mesh->GetNumFaces()); ++i)
				//	swap(mesh->m_indices[i * 3 + 0], mesh->m_indices[i * 3 + 1]);
				//mesh->CalculateNormals(); //??
				use_mesh = true;
			}
		}
		std::cout << "parse the proxy and create a rigid body asset " << ingr_node["ingtype"].asString() << " "<< lodproxy_to_use << endl;
		//parse the proxy and create a rigid body asset
		//check if "Type":"Grow"
		IngredientSphereTree ingr_spheres;
		if (ingr_node["Type"].asString() == "Grow" || ingr_node["ingtype"].asString()=="fiber")
		{
			std::cout << "fiber" << endl;
			//model as rope in flex
			std::vector<Vec3> points_to_use;
			pnames_fiber.push_back(ingr_node["name"].asString());
			pnames_fiber_nodes.push_back(ingr_node);
			std::cout << "ok 2 " << endl;
			AssetBatch b;//for ribid body protein
			int p = lodproxy_to_use;
			std::cout << "lodproxy_to_use "<< lodproxy_to_use << endl;
			Json::Value jsonpos = ingr_node["positions"][p]["coords"];
			Quat raxe = AlignVec3s(Vec3(0, 0, 1), pcpalVector);
			if (pcpalVector == Vec3(0, 0, 1) || pcpalVector == Vec3(0, 0, -1))
				raxe = get_rotation_between(pcpalVector, Vec3(0, 0, 1));
			if (isnan(raxe.x) || isnan(raxe.y) || isnan(raxe.z)) {
				std::cout << "raxe " << raxe.x << " " << raxe.y << " " << raxe.z << endl;
				std::cout << "ingr pcpalVector  : " << pcpalVector.x << " " << pcpalVector.y << " " << pcpalVector.z << endl;
				raxe = Quat(0, 0, 0, 1);
			}
			std::cout << "ingr jsonpos.size()  : " << jsonpos.size() << endl;
			for (int i = 0; i < jsonpos.size() / 3; i++) {//1
				if (jsonpos[i * 3].asString() == "NaN" ||
					jsonpos[i * 3 + 1].asString() == "NaN" ||
					jsonpos[i * 3 + 2].asString() == "NaN") {
					// std::cout << " i is Nan" << i << endl;
					continue;
				}
				Vec3 beadp = Vec3(jsonpos[i * 3].asFloat(),
					jsonpos[i * 3 + 1].asFloat(),
					jsonpos[i * 3 + 2].asFloat());
				if (isnan(beadp.x) || isnan(beadp.y) || isnan(beadp.z)) {
					beadp = Vec3(std::stof(jsonpos[i * 3].asString()),
						std::stof(jsonpos[i * 3 + 1].asString()),
						std::stof(jsonpos[i * 3 + 2].asString()));
					std::cout << i << " " << typeid(jsonpos[i * 3]).name() << endl;
					std::cout << i << " " << beadp.x << " " << beadp.y << " " << beadp.z << " " << jsonpos[i * 3] << " " << jsonpos[i * 3 + 1] << " " << jsonpos[i * 3 + 2] << endl;
				}
				if (isnan(offset.x) || isnan(offset.y) || isnan(offset.z)) {
					std::cout << i << " offset " << offset.x << " " << offset.y << " " << offset.z << endl;
				}
				if (ingr_node["source"]["transform"]["center"].asBool() && (force_not_center == 0)) {
					beadp = Rotate(raxe, (beadp + offset));
				}
				ingr_spheres.LevelPoints.push_back(beadp*main_scale);
				if (isnan(beadp.x) || isnan(beadp.y) || isnan(beadp.z)) {
					std::cout << typeid(jsonpos[i * 3]).name() << " " << std::stof(jsonpos[i * 3 + 2].asString()) << endl;
					std::cout << beadp.x << " " << beadp.y << " " << beadp.z << " " << jsonpos[i * 3] << " " << jsonpos[i * 3 + 1] << " " << jsonpos[i * 3 + 2] << endl;
				}
			}
			std::cout << "fiber ingr positions  nb " << ingr_spheres.LevelPoints.size() << endl;
			points_to_use = ingr_spheres.LevelPoints;
			b.mAsset = flexExtCreateRigidFromPoints(points_to_use);
			b.compId = compId;
			b.offsetx = offset.x* main_scale;
			b.offsety = offset.y* main_scale;
			b.offsetz = offset.z* main_scale;
			b.pcpalVectorx = pcpalVector.x;
			b.pcpalVectory = pcpalVector.y;
			b.pcpalVectorz = pcpalVector.z;
			b.nInstances = 0;
			b.ingr_name = ingr_node["name"].asString().c_str();
			iBatchesFiber.push_back(b);
			std::cout << "fiber ingr name : " << ingr_node["name"].asString() << " " << pnames_fiber.size()-1 << endl;
		}
		else
		{
			//use proxy file or positions
			std::vector<Vec3> points_to_use;
			Json::Value  sphereFile = ingr_node["sphereFile"];

			//std::cout << "try loading " << ingr_node["sphereFile"] << " " << sphereFile.asString() << " " << datapath + sphereFile.asString() << endl;
			string filename = sphereFile.asString();

			if (filename.substr(filename.find_last_of(".") + 1) == "sph")
				filename = ingr_source_pdb + ".pdb_kmeans15.txt";

			if (!sphereFile.empty())
			{
				//filename = ingr_source_pdb + ".pdb_kmeans15.txt";
				//if (compId > 0){
				//	replace(ingr_source_pdb, "_mb", "");
				//}
				//filename = ingr_source_pdb + "_cl.txt";
				//filename = sphereFile;
				//std::cout << "load " << datapath + filename << endl;
				ingr_spheres = parseProxy(datapath + filename);
				//use lvl 0
				for (int i = 0; i < ingr_spheres.LevelCounts[0]; i++){
					points_to_use.push_back(ingr_spheres.LevelPoints[i] * main_scale);
				}
				std::cout << "sphereFile ingr proxy nb " << points_to_use.size() << endl;
			}
			else {
				std::cout << "load from positions " << ingr_node["positions"].size() << " " << ingr_node["positions"][0].size() << endl;
				int p = lodproxy_to_use;
				//if (ingr_node["positions"].size() > 1)
				//		p = 1;
				Json::Value jsonpos = ingr_node["positions"][p]["coords"];
				if (overwrite_radius)
				{
					std::cout << "overwrite_radius" << endl;
					Json::Value jsonrad = ingr_node["radii_lod"][p]["radii"];
					std::cout << jsonrad << endl;
					std::cout << jsonrad[0] << endl;
					std::cout << std::stof(jsonrad[0].asString()) << endl;
					std::cout << std::stof(jsonrad[0].asString())*main_scale << endl;
					main_radius = std::stof(jsonrad[0].asString())*main_scale;
					overwrite_radius = false;
				}
				Quat raxe = AlignVec3s(Vec3(0, 0, 1), pcpalVector);
				if (pcpalVector == Vec3(0, 0, 1) || pcpalVector == Vec3(0, 0, -1))
					raxe = get_rotation_between(pcpalVector, Vec3(0, 0, 1));
				// if (pcpalVector == Vec3(0, 0, 1)) raxe = Quat(0, 0, 0, 1);
				if (isnan(raxe.x) || isnan(raxe.y) || isnan(raxe.z)) {
					std::cout << "raxe " << raxe.x << " " << raxe.y << " " << raxe.z << endl;
					std::cout << "ingr pcpalVector  : " << pcpalVector.x << " " << pcpalVector.y << " " << pcpalVector.z << endl;
					raxe = Quat(0, 0, 0, 1);
				}
				std::cout << "raxe " << raxe.x << " " << raxe.y << " " << raxe.z << endl;
				std::cout << "jsonpos.size()/3 " << (jsonpos.size() / 3) << endl;
				for (int i = 0; i < jsonpos.size()/3; i++){//1
					if (jsonpos[i * 3].asString() == "NaN" ||
						jsonpos[i * 3 + 1].asString() == "NaN" ||
						jsonpos[i * 3 + 2].asString() == "NaN") {
						//std::cout << " i is Nan" << i << endl;
						continue;
					}
					Vec3 beadp = Vec3(jsonpos[i * 3].asFloat(),
						jsonpos[i * 3 + 1].asFloat(),
						jsonpos[i * 3 + 2].asFloat());
					if (isnan(beadp.x) || isnan(beadp.y) || isnan(beadp.z)) {
						beadp = Vec3(std::stof(jsonpos[i * 3].asString()),
							std::stof(jsonpos[i * 3 + 1].asString()),
								std::stof(jsonpos[i * 3 + 2].asString()));
						std::cout << i << " " <<  typeid(jsonpos[i * 3]).name() << endl;
						std::cout << i << " " << beadp.x << " " << beadp.y << " " << beadp.z << " " << jsonpos[i * 3] << " " << jsonpos[i * 3 + 1] << " " << jsonpos[i * 3 + 2] << endl;
					}
					if (isnan(offset.x) || isnan(offset.y) || isnan(offset.z)) {
						std::cout << i << " offset " << offset.x << " " << offset.y << " " << offset.z << endl;
					}

					//if (ingr_node["source"]["transform"]["center"].asBool())
					//{
					//beadp = raxe * (beadp + offset);
					// std::cout << "force_not_center " << force_not_center << endl;
					if (ingr_node["source"]["transform"]["center"].asBool() && (force_not_center == 0) )
					{
						// std::cout << "force_not_center " << force_not_center << endl;
						beadp = Rotate(raxe, (beadp + offset));
					}
					
					//std::cout  << jsonpos[i * 3] << " " << jsonpos[i * 3 + 1] << " " << jsonpos[i * 3 + 2] << endl;
					//std::cout << beadp.x << " " << beadp.y << " " << beadp.z << endl;
					//}
					ingr_spheres.LevelPoints.push_back(beadp*main_scale);
					if (isnan(beadp.x) || isnan(beadp.y) || isnan(beadp.z)) {
						std::cout << typeid(jsonpos[i * 3]).name() << " " << std::stof(jsonpos[i * 3 + 2].asString()) << endl;
						std::cout << beadp.x<<" "<< beadp.y<< " " << beadp.z << " " << jsonpos[i * 3] << " " << jsonpos[i * 3 + 1] << " " << jsonpos[i * 3 + 2] << endl;
					}
					//std::cout << jsonpos[i*3].asFloat() << endl;
				}
				std::cout << "ingr positions  nb " << ingr_spheres.LevelPoints.size() << endl;
				points_to_use = ingr_spheres.LevelPoints;
			}

			if (points_to_use.size() == 0) {
				std::cout << "ingr proxy  0 return " << points_to_use.size() << endl;
				//try to do it from the mesh ?
				//return;
			}
			std::cout << "points center ?  : " << ingr_node["source"]["transform"]["center"].asBool() << endl;
			std::cout << "compId " << compId << endl;
			//before creating the asset, if surface, remove beads at intersection of the enveloppe.
			use_mesh = false;
			AssetBatch b;
			if (use_rb) {
				//if (compId > 0 ) 
				//	points_to_use = filterSurfaceBeads(offset, pcpalVector, points_to_use);
				//const float spacing = radius*0.5f;
				//b.mAsset = NvFlexExtCreateRigidFromMesh((float*)&mesh->m_positions[0], 
				//			int(mesh->m_positions.size()), (int*)&mesh->m_indices[0], 
				//			mesh->m_indices.size(), spacing, -spacing*0.5f);
				if (use_mesh) {
					const float spacing = g_params.radius *0.5f;
					b.mAsset = NvFlexExtCreateRigidFromMesh((float*)&mesh->m_positions[0],
						int(mesh->m_positions.size()), (int*)&mesh->m_indices[0], 
						mesh->m_indices.size(), spacing, -spacing*0.5f);
				}
				else {
					if (points_to_use.size() != 0) {
						b.mAsset = flexExtCreateRigidFromPoints(points_to_use, ingr_node["source"]["transform"]["center"].asBool());
					}
				}
				//b.mMesh = CreateGpuMesh(mesh);
				if (!meshFile.empty()) 
					if (use_instances_mesh)
						b.nvMesh = CreateTriangleMesh(mesh);
			}
			//const float* points, int numParticles, float clusterSpacing, float clusterRadius, float clusterStiffness, float linkRadius, float linkStiffness
			else {
				//b.mAsset = flexExtCreateSoftFromPoints((float*)&points_to_use[0], points_to_use.size(),
				//0.10f,
				//0.10f,
				//0.05f,
				//g_params.radius,// / 2.0f,
				//1.00f);
				//overwrite the string
			}
			//b.mMesh = nullptr;//CreateGpuMesh(mesh);
			b.compId = compId;
			b.offsetx = offset.x* main_scale;
			b.offsety = offset.y* main_scale;
			b.offsetz = offset.z* main_scale;
			b.pcpalVectorx = pcpalVector.x;
			b.pcpalVectory = pcpalVector.y;
			b.pcpalVectorz = pcpalVector.z;
			b.nInstances = 0;
			b.ingr_name = ingr_node["name"].asString().c_str();
			iBatches.push_back(b);
			pnames.push_back(ingr_node["name"].asString());
			proteins_nodes.push_back(ingr_node);
			std::cout << "ingr name : " << ingr_node["name"].asString() << " " << proteins_nodes.size() - 1 << " ingrIndex " << iBatches.size() - 1 << endl;
			if (ingr_node["name"].asString() == "MG_404_MONOMER") {
				std::cout << "particle count " << b.mAsset->numParticles << " " << points_to_use.size() << endl;
			}
			//std::cout << "ingr index  " << iBatches.size() - 1 << " " << pnames.size() - 1 << " " << ingr_node["name"].asString() << endl;
			//std::cout << "particle count " << b.mAsset->numParticles << endl;
			//std::cout << "nShape " << b.mAsset->numShapes << endl;
		}
		mIngrSphereTree.push_back(ingr_spheres);
	}

	//parse ingredient in a compartment
	void parseJsonIngredients(Json::Value  ingr_nodes, int comp){
		//iteration is using alphabetic order?
		//JsonCpp keeps its values in a std::map<CZString, Value>, which is always sorted by the CZString comparison,
		for (Json::ValueIterator itr = ingr_nodes.begin(); itr != ingr_nodes.end(); itr++) {
			Json::Value ingr_node_name = *itr;
			ingredientToFlexAsset(ingr_node_name, comp);
			string ingr_name = ingr_node_name["name"].asString();
			Json::Value ingr_node_results = ingr_node_name["results"];
			totalNBMol += ingr_node_name["nbMol"].asInt();
			if (ingr_node_results.size() != 0)
			{
				std::cout << "array size ? #nb of instances" << ingr_node_results.size() << endl;
			}
		}
	}

	//parse ingredient in a compartment
	void parseJsonIngredientsSerialized(Json::Value  ingr_nodes, int comp) {
		//iteration is using alphabetic order?
		//JsonCpp keeps its values in a std::map<CZString, Value>, which is always sorted by the CZString comparison,
		int nIngredients = ingr_nodes.size();
		for (int i = 0; i < nIngredients; i++) {
			Json::Value ingr_node_name = ingr_nodes[i];
			ingredientToFlexAsset(ingr_node_name, comp);
			string ingr_name = ingr_node_name["name"].asString();
			Json::Value ingr_node_results = ingr_node_name["results"];
			totalNBMol += ingr_node_name["nbMol"].asInt();
			if (ingr_node_results.size() != 0)
			{
				std::cout << "array size ? #nb of instances" << ingr_node_results.size() << endl;
			}
			//"buildtype":"supercell"
			if (ingr_node_name["buildtype"].asString() == "supercell" ) {
				std::cout << ingr_node_name["name"] << "is a supercell ingredient skip and replace with big sphere" << endl;
				//sphere radius is 1.0
				//use a margin that offset the dimension
				//default margin =0.1f
				//NvFlexDistanceFieldId sdf = CreateSDF(GetFilePathByPlatform("../../data/sphere.ply").c_str(), dim, 0.0f,0.0f);
				//sdf, position, rotation, width/scale
				//float radius = 2000.0f * main_scale;
				//AddSDF(sdf, Vec3(-radius/2.0f, -radius/2.0f, -radius/2.0f), Quat(0.0f, 0.0f, 0.0f, 1.0f), radius);
			}
		}
	}
	

	void printBatchId(){
		for (int i = 0; i < pnames.size(); i++)
		{
			std::cout << "ingredient " << pnames[i] << " id " << i << endl;
		}
		for (int i = 0; i < pnames_fiber.size(); i++){
			std::cout << "ingredient fiber " << pnames_fiber[i] << " id " << i << endl;
		}
	}
	
	int getIngredientBatchId(string name){
		int index = 0;
		bool found = false;
		for (int i = 0; i <pnames.size(); i++)
		{
			if (pnames[i] == name) {
				index = i;
				found = true;
				break;
			}
		}
		if (!found)
		{
			for (int i = 0; i <pnames_fiber.size(); i++)
			{
				if (pnames_fiber[i] == name) {
					index = i;
					found = true;
					break;
				}
			}
		}
		return index;
	}


	void parseJsonIngredientsResults(Json::Value  ingr_nodes, int comp, bool redo = false){
		for (Json::ValueIterator itr = ingr_nodes.begin(); itr != ingr_nodes.end(); itr++) {
			Json::Value ingr_node_name = *itr;
			string ingr_name = ingr_node_name["name"].asString();
			Json::Value ingr_node_results = ingr_node_name["results"];
			//JSONNode ingr_node_results_data =  ingr_node_results.as_array();//array or array
			//if (ingr_name != "DNA-binding protein HU")
			//	continue;
			if (ingr_node_results.size() != 0)
			{//do something
				//continue;
				std::cout << ingr_name << " array size ? #nb of instances " << ingr_node_results.size() << endl;
				//poopulate
				int ingrIndex = getIngredientBatchId(ingr_name);
				int nResult = ingr_node_results.size();
				for (int i = 0; i < nResult; i++)
				{
					Vec3 pos = Vec3(ingr_node_results[i][0][0].asFloat(),
						ingr_node_results[i][0][1].asFloat(),
						ingr_node_results[i][0][2].asFloat())*main_scale;
					Quat quat = Quat(ingr_node_results[i][1][0].asFloat(),
						ingr_node_results[i][1][1].asFloat(),
						ingr_node_results[i][1][2].asFloat(),
						ingr_node_results[i][1][3].asFloat());

					if (!redo) createInstanceIngredient(ingrIndex, pos, quat);
					else updateInstanceIngredient(ingrIndex, pos, quat, i);
					//if (i == 2) return;
					//break;
				}
			}//if curve it will have the controls_points as "curves"[[],[],[]]
			//nbCurve and curve$id
			//list of Vector3

			Json::Value ingr_curve_results = ingr_node_name["nbCurve"];
			if (ingr_curve_results != NULL) {
				int persistence = 2;
				int fiber_id = getFiberId(ingr_name);
				if (fiber_id==-1) continue;
				int compId = iBatchesFiber[fiber_id].compId;
				cout << "build rope " << ingr_name<< " id " << fiber_id << endl;
				float L = ingr_node_name["length"].asFloat()*main_scale;
				bool close = ingr_node_name["closed"].asBool();
				float D = g_params.radius;				
				int nCurve = ingr_curve_results.asInt();
				for (int i = 0; i < nCurve; i++) {
					Json::Value curve = ingr_node_name["curve"+ to_string(i)];
					int npoints = curve.size();
					float *data_curve = new float[npoints * 3];
					for (int i = 0; i < npoints; i++){
						data_curve[i*3] = curve[i][0].asFloat();
						data_curve[i*3+1] = curve[i][1].asFloat();
						data_curve[i*3+2] = curve[i][2].asFloat();
					}
					//extend = false;
					int rope_phase = NvFlexMakePhase(iGroupCounter++, eNvFlexPhaseSelfCollide);
					Rope rcurve;
					CreateRopeFromData(rcurve, //rope
						main_scale, // scale
						1.0f, //stifness
						data_curve, //data
						L, //length
						npoints * 3,  //nfloat
						rope_phase,//phase
						0.0f,//spiral angle
						1.0f,//invmass
						0.0f,//give
						0,//extend_nb
						false,//extend
						close,//close
						main_radius*2.0f,
						persistence);
					rcurve.persistence = persistence;
					printf("create Rope from data OK with %i points\n", rcurve.mIndices.size());
					g_ropes.push_back(rcurve);//instance
					mask.push_back(-1);
					maks_protein.push_back(-(fiber_id + 1));//it ptype is 0 ?
					maks_fiber.push_back(fiber_id);//proteinType
					printf("create Rope from data OK with %i  %i\n", -(fiber_id + 1), fiber_id);
				}
			}
		}
	}

	void updateInstanceIngredient(int ingrIndex, Vec3 position, Quat rotation, int i){
		Instance& inst = mInstances[i];
		NvFlexExtAsset* asset = iBatches[inst.mMeshIndex].mAsset;
		for (int j = 0; j < asset->numParticles; ++j)
		{
			Vec3 localPos = Vec3(&asset->particles[j * 4]) - Vec3(&asset->shapeCenters[0]);//local position of the proxy
			Vec3 rpos = Rotate(rotation, localPos);
			g_buffers->positions[inst.mParticleOffset + j] = Vec4(position + rpos, 1.0f);
			g_buffers->velocities[inst.mParticleOffset + j] = Vec3(0, 0, 0);
		}
	}

	void UpdateInstanceTransform(int instId, Vec3 pos, Quat quat)
	{
		Instance& inst = mInstances[instId];
		NvFlexExtAsset* asset = iBatches[inst.mMeshIndex].mAsset;
		for (int j = 0; j < asset->numParticles; ++j)
		{
			Vec3 localPos = Vec3(&asset->particles[j * 4]);//- Vec3(&asset->shapeCenters[0]);local position of the proxy
			Vec3 rpos = Rotate(quat, localPos);
			g_buffers->positions[inst.mParticleOffset + j] = Vec4(pos + rpos, 1.0f);
		}
		g_buffers->rigidTranslations[instId].Set(pos.x, pos.y, pos.z);
		g_buffers->rigidRotations[instId].Set(quat.x, quat.y, quat.z, quat.w);
	}

	int createInstanceIngredient(int ingrIndex, Vec3 position, Quat rotation, bool rb = false, float mass = 1.0f)
	{

		const int particleOffset = g_buffers->positions.size();
		const int indexOffset = g_buffers->rigidOffsets.back();

		NvFlexExtAsset* asset = iBatches[ingrIndex].mAsset;
		//iBatches[inst.mMeshIndex].mAsset
		Json::Value ingr_node = getProteinNode(ingrIndex);// (iBatches[ingrIndex].ingr_name);
		//"buildtype":"supercell"
		if (ingr_node["buildtype"].asString() == "supercell" ) {
			std::cout << iBatches[ingrIndex].ingr_name << "is a supercell ingredient skip and replace with big sphere" << endl;
			//sphere radius is 1.0
			//NvFlexDistanceFieldId sdf = CreateSDF(GetFilePathByPlatform("../../data/sphere.ply").c_str(), dim);
			IngredientToFixSDF(ingr_node, position);
			//sdf, position, rotation, width/scale
			//float radius = 1000.0f * main_scale;
			//AddSDF(sdf, Vec3(0.0f, 0.0f, 0.0f), Quat(0.0f, 0.0f, 0.0f, 1.0f), radius);
			//create sphere sdf or mesh 
			//return -1;
		}
		
		if (use_instances_mesh) {
			std::cout << " use_instances_mesh " << use_instances_mesh << endl;
			const int dim = 128; 
			if (iBatches[ingrIndex].nvMesh!=NULL) AddTriangleMesh(iBatches[ingrIndex].nvMesh, position, rotation, 1.0f);//add colliding mesh ? but not the particles ?
			return 0;
		}
		else {
			Quat q = rotation;// Quat();

			Instance inst;

			inst.mLifetime = 1000.0f;// Randf(5.0f, 60.0f);
			inst.mParticleOffset = particleOffset;
			inst.mRotation = q;// q;// Quat(0, 0, 0, 1);// 
			inst.bounded = 0;
			inst.mTranslation = Vec3(position.x, position.y, position.z);

			inst.mMeshIndex = ingrIndex;

			Vec3 linearVelocity = Vec3(0.0f);//g_emitters[0].mDir*15.0f;//*Randf(5.0f, 10.0f);//Vec3(Randf(0.0f, 10.0f), 0.0f, 0.0f);
			Vec3 angularVelocity = Vec3(0.0f);//Vec3(UniformSampleSphere()*Randf()*k2Pi);

			mask.push_back(iBatches[ingrIndex].compId);
			maks_protein.push_back(ingrIndex);
			inst.mGroup = iGroupCounter++;

			//int mm = 2 ^ 24;
			if (g_buffers->rigidIndices.empty())
				g_buffers->rigidOffsets.push_back(0);

			int phase = NvFlexMakePhase(inst.mGroup, 0);
			//NvFlexMakePhaseWithChannels(i, eNvFlexPhaseSelfCollide, eNvFlexPhaseShapeChannel0<<i)//compId 
			//if (!use_rb) phase = NvFlexMakePhase(inst.mGroup, eNvFlexPhaseSelfCollide | eNvFlexPhaseSelfCollideFilter);
			// std::cout << " instance for " << ingrIndex << " compid " << iBatches[ingrIndex].compId << endl;inst1.mMeshIndex == 12
			mass = 1;
			// if (ingr_node["name"].asString() == "PDH") mass = 0;
			//special case
			// ingr_node["partners_properties"].size() != 0
			
			// if (ingr_node["name"].asString() == "mRNA") mass = 0;
			// if (ingr_node["name"].asString() == "POL_CA") mass = 0;
			// if (ingr_node["name"].asString() == "NP_CA") mass = 0;
			// hard coded mass == 0 ? 
			// or use the info binary data ?
			if (iBatches[ingrIndex].compId > 0 || ingr_node["buildtype"].asString() == "supercell" ) {
				// we need a more advanced filtering here for the masking. In unity we have a field for masking that should be used here...
				// std::cout << " mass 0 for ingrIndex " << ingrIndex << " " << ingr_node["name"].asString() << " compid " << iBatches[ingrIndex].compId << endl;
				mass = 0;
			}
			for (int j = 0; j < asset->numParticles; j++)
			{
				Vec3 localPos = Vec3(&asset->particles[j * 4]) - Vec3(&asset->shapeCenters[0]);//cluster position
				g_buffers->rigidLocalPositions.push_back(localPos);
				//What about rigidLocalNormals nx,ny,nz,distance
				g_buffers->rigidIndices.push_back(int(g_buffers->positions.size()));
				Vec3 rpos = Rotate(inst.mRotation, Vec3(&asset->particles[j * 4]));//localPos);// 
				//test if part inside
				rpos = inst.mTranslation + rpos;
				// std::cout << j <<" "<< localPos.x << " " << localPos.y << " " << localPos.z << endl;
				//if (iBatches[ingrIndex].compId < 0) 
				//{
					//if (!isInside(rpos, 0))
				Vec3 vel = Vec3(0.0f);
				float doff = Length(rpos);
				if (iBatches[ingrIndex].compId < 0 && doff  >= (comp_radius - g_params.radius*2)) {//add the membrane thickness to it
					//std::cout << j << " " << iBatches[ingrIndex].ingr_name << " " << ingrIndex << " outside " << Length(rpos) << " radius " << comp_radius << endl;
					rpos = (Normalize(rpos) * (comp_radius - g_params.radius * 2) );
					//std::cout << " after " << Length(rpos) << endl;
					vel = -Normalize(rpos) * g_params.radius / g_dt;
					//inst.mTranslation = inst.mTranslation - Normalize(rpos) * ((doff - comp_radius) * 2.0f);
				}
				//}
				g_buffers->positions.push_back(Vec4(rpos, mass));//Vec4(localPos,mass));// Vec4(inst.mTranslation + inst.mRotation * localPos, mass));// Vec4(inst.mTranslation + rpos, 1.0f));//inst.mRotation*
				g_buffers->velocities.push_back(vel);//linearVelocity + Cross(angularVelocity, localPos);
				g_buffers->phases.push_back(phase);
			}
			//what if fix one only ?
			g_buffers->rigidCoefficients.push_back(1.0f);//0.15?
			g_buffers->rigidTranslations.push_back(inst.mTranslation);//0?
			g_buffers->rigidRotations.push_back(inst.mRotation);//);

			g_buffers->rigidOffsets.push_back(int(g_buffers->rigidIndices.size()));

			const int numRigids = g_buffers->rigidOffsets.size() - 1;

			// calculate local rest space positions
			//g_buffers->rigidLocalPositions.resize(g_buffers->rigidOffsets.back());
			//CalculateRigidLocalPositions(&g_buffers->positions[0], g_buffers->positions.size(),
			//	&g_buffers->rigidOffsets[0], &g_buffers->rigidIndices[0],
			//	numRigids, &g_buffers->rigidLocalPositions[0]);


			//particleOffset += asset->numParticles;
			mInstances.push_back(inst);
			iBatches[ingrIndex].nInstances++;
			// Draw transform
			Matrix44 xform = TranslationMatrix(Point3(inst.mTranslation - inst.mRotation*Vec3(asset->shapeCenters)))*RotationMatrix(inst.mRotation);
			iBatches[ingrIndex].mInstanceTransforms.push_back(xform);
			return inst.mParticleOffset;
		}
		

		//if (use_rb) UpdateInstanceTransform(mInstances.size() - 1, position, rotation);
		
	}

	int createInstanceIngredientAt(int ingrIndex, Vec3 position, Quat rotation, int instance_id, int offset, bool rb = false)
	{

		int particleOffset = offset;
		const int indexOffset = g_buffers->rigidOffsets.back();
		NvFlexExtAsset* asset = iBatches[ingrIndex].mAsset;

		Quat q = rotation;// Quat();

		Instance inst;

		inst.mLifetime = 1000.0f;// Randf(5.0f, 60.0f);
		inst.mParticleOffset = particleOffset;
		inst.mRotation = q;// q;// Quat(0, 0, 0, 1);// 

		inst.mTranslation = Vec3(position.x, position.y, position.z);

		inst.mMeshIndex = ingrIndex;

		Vec3 linearVelocity = Vec3(0.0f);//g_emitters[0].mDir*15.0f;//*Randf(5.0f, 10.0f);//Vec3(Randf(0.0f, 10.0f), 0.0f, 0.0f);
		Vec3 angularVelocity = Vec3(0.0f);//Vec3(UniformSampleSphere()*Randf()*k2Pi);

		mask.push_back(iBatches[ingrIndex].compId);
		maks_protein.push_back(ingrIndex);
		inst.mGroup = iGroupCounter++;

		//int mm = 2 ^ 24;
		
		int phase = NvFlexMakePhase(inst.mGroup, 0);
		if (!use_rb) phase = NvFlexMakePhase(inst.mGroup, eNvFlexPhaseSelfCollide | eNvFlexPhaseSelfCollideFilter);
		for (int j = 0; j < asset->numParticles; ++j)
		{
			Vec3 localPos = Vec3(&asset->particles[j * 4]) - Vec3(&asset->shapeCenters[0]);
			//Vec3 rpos = Rotate(inst.mRotation, localPos);
			g_buffers->positions[inst.mParticleOffset + j] = Vec4(inst.mTranslation + inst.mRotation*localPos, 1.0f);
			g_buffers->velocities[inst.mParticleOffset + j] = Vec3(0.0f);
			g_buffers->phases[inst.mParticleOffset + j] = phase;
			// rebuild active indices
			
			//for (int i = 0; i < numActive; ++i)
			//	g_buffers->activeIndices[i] = i;
			//g_buffers->positions.push_back(Vec4(inst.mTranslation + inst.mRotation*localPos, 1.0f));// Vec4(inst.mTranslation + rpos, 1.0f));//inst.mRotation*
			//g_buffers->velocities.push_back(Vec3(0.0f));//linearVelocity + Cross(angularVelocity, localPos);
			//g_buffers->phases.push_back(phase);
		}


		g_buffers->rigidCoefficients.push_back(1.0f);
		g_buffers->rigidTranslations.push_back(inst.mTranslation);
		g_buffers->rigidRotations.push_back(inst.mRotation);

		for (int j = 0; j < asset->numShapeIndices; ++j)
		{
			g_buffers->rigidLocalPositions.push_back(Vec3(&asset->particles[j * 4]) - Vec3(&asset->shapeCenters[0]));
			g_buffers->rigidIndices.push_back(asset->shapeIndices[j] + numActive);
		}

		g_buffers->rigidOffsets.push_back(g_buffers->rigidIndices.size());
		//}
		const int numRigids = g_buffers->rigidOffsets.size() - 1;

		// calculate local rest space positions
		//g_buffers->rigidLocalPositions.resize(g_buffers->rigidOffsets.back());
		//CalculateRigidLocalPositions(&g_buffers->positions[0], g_buffers->positions.size(), 
		//							&g_buffers->rigidOffsets[0], &g_buffers->rigidIndices[0], 
		//							numRigids, &g_buffers->rigidLocalPositions[0]);
		
		//Matrix44 xform = TranslationMatrix(Point3(inst.mTranslation - inst.mRotation*Vec3(asset->shapeCenters)))*RotationMatrix(inst.mRotation);
		//iBatches[inst.mMeshIndex].mInstanceTransforms.push_back(xform);

		g_buffers->activeIndices.resize(particleOffset + asset->numParticles);
		for (int i = 0; i < g_buffers->activeIndices.size(); ++i)
			g_buffers->activeIndices[i] = i;

		//particleOffset += asset->mNumParticles;
		particleOffset += asset->numParticles;
		mInstances[instance_id]=inst;
		iBatches[ingrIndex].nInstances++;
		cout << " created instance " << instance_id << " for ingredient " << ingrIndex << endl;
		//if (use_rb) UpdateInstanceTransform(mInstances.size() - 1, position, rotation);
		return particleOffset;
	}

	int createInstanceIngredientFiber(int ingrIndex, Vec3 position, Quat rotation, float mass = 1.0f, bool rb = false)
	{
		std::cout << " instance for " << ingrIndex << endl;
		const int particleOffset = g_buffers->positions.size();
		const int indexOffset = g_buffers->rigidOffsets.back();

		NvFlexExtAsset* asset = iBatchesFiber[ingrIndex].mAsset;

		Json::Value ingr_node = getFiberNode(ingrIndex);// (iBatches[ingrIndex].ingr_name);
		// std::cout << " instance for node " << ingr_node << endl;

		if (use_instances_mesh) {
			std::cout << " use_instances_mesh " << use_instances_mesh << endl;
			const int dim = 128;
			if (iBatchesFiber[ingrIndex].nvMesh != NULL) AddTriangleMesh(iBatchesFiber[ingrIndex].nvMesh, position, rotation, 1.0f);//add colliding mesh ? but not the particles ?
			return 0;
		}
		else {
			Quat q = rotation;// Quat();

			Instance inst;

			inst.mLifetime = 1000.0f;// Randf(5.0f, 60.0f);
			inst.mParticleOffset = particleOffset;
			inst.mRotation = q;// q;// Quat(0, 0, 0, 1);// 

			inst.mTranslation = Vec3(position.x, position.y, position.z);

			inst.mMeshIndex = ingrIndex;

			Vec3 linearVelocity = Vec3(0.0f);//g_emitters[0].mDir*15.0f;//*Randf(5.0f, 10.0f);//Vec3(Randf(0.0f, 10.0f), 0.0f, 0.0f);
			Vec3 angularVelocity = Vec3(0.0f);//Vec3(UniformSampleSphere()*Randf()*k2Pi);

			mask.push_back(iBatchesFiber[ingrIndex].compId);
			// maks_fiber.push_back(ingrIndex);
			inst.mGroup = iGroupCounter++;

			//int mm = 2 ^ 24;
			if (g_buffers->rigidIndices.empty())
				g_buffers->rigidOffsets.push_back(0);

			int phase = NvFlexMakePhase(inst.mGroup, 0);
			//NvFlexMakePhaseWithChannels(i, eNvFlexPhaseSelfCollide, eNvFlexPhaseShapeChannel0<<i)//compId 
			//if (!use_rb) phase = NvFlexMakePhase(inst.mGroup, eNvFlexPhaseSelfCollide | eNvFlexPhaseSelfCollideFilter);
			std::cout << " instance for " << ingrIndex << " compid " << iBatchesFiber[ingrIndex].compId << endl;
			//std::cout << iBatches[ingrIndex].ingr_name  << " pos " << inst.mTranslation.x << " d " << Length(inst.mTranslation) << " radius " << comp_radius << " " << asset->numParticles  << endl;
			for (int j = 0; j < asset->numParticles; j++)
			{
				Vec3 localPos = Vec3(&asset->particles[j * 4]) - Vec3(&asset->shapeCenters[0]);//cluster position
				g_buffers->rigidLocalPositions.push_back(localPos);
				//What about rigidLocalNormals nx,ny,nz,distance
				g_buffers->rigidIndices.push_back(int(g_buffers->positions.size()));
				Vec3 rpos = Rotate(inst.mRotation, Vec3(&asset->particles[j * 4]));//localPos);// 
				//mass = 1;
				if (iBatchesFiber[ingrIndex].compId > 0 || ingr_node["buildtype"].asString() == "supercell") mass = 0;
				//test if part inside
				rpos = inst.mTranslation + rpos;
				// std::cout << j <<" "<< localPos.x << " " << localPos.y << " " << localPos.z << endl;
				//if (iBatches[ingrIndex].compId < 0) 
				//{
					//if (!isInside(rpos, 0))
				Vec3 vel = Vec3(0.0f);
				float doff = Length(rpos);
				if (iBatchesFiber[ingrIndex].compId < 0 && doff >= (comp_radius - g_params.radius * 2)) {//add the membrane thickness to it
					//std::cout << j << " " << iBatches[ingrIndex].ingr_name << " " << ingrIndex << " outside " << Length(rpos) << " radius " << comp_radius << endl;
					rpos = (Normalize(rpos) * (comp_radius - g_params.radius * 2));
					//std::cout << " after " << Length(rpos) << endl;
					vel = -Normalize(rpos) * g_params.radius / g_dt;
					//inst.mTranslation = inst.mTranslation - Normalize(rpos) * ((doff - comp_radius) * 2.0f);
				}
				//}
				g_buffers->positions.push_back(Vec4(rpos, mass));//Vec4(localPos,mass));// Vec4(inst.mTranslation + inst.mRotation * localPos, mass));// Vec4(inst.mTranslation + rpos, 1.0f));//inst.mRotation*
				g_buffers->velocities.push_back(vel);//linearVelocity + Cross(angularVelocity, localPos);
				g_buffers->phases.push_back(phase);
			}
			//what if fix one only ?
			g_buffers->rigidCoefficients.push_back(1.0f);//0.15?
			g_buffers->rigidTranslations.push_back(inst.mTranslation);//0?
			g_buffers->rigidRotations.push_back(inst.mRotation);//);

			g_buffers->rigidOffsets.push_back(int(g_buffers->rigidIndices.size()));

			const int numRigids = g_buffers->rigidOffsets.size() - 1;

			// calculate local rest space positions
			//g_buffers->rigidLocalPositions.resize(g_buffers->rigidOffsets.back());
			//CalculateRigidLocalPositions(&g_buffers->positions[0], g_buffers->positions.size(),
			//	&g_buffers->rigidOffsets[0], &g_buffers->rigidIndices[0],
			//	numRigids, &g_buffers->rigidLocalPositions[0]);


			//particleOffset += asset->numParticles;
			mInstancesFiber.push_back(inst);
			iBatchesFiber[ingrIndex].nInstances++;
			// Draw transform
			Matrix44 xform = TranslationMatrix(Point3(inst.mTranslation - inst.mRotation*Vec3(asset->shapeCenters)))*RotationMatrix(inst.mRotation);
			iBatchesFiber[ingrIndex].mInstanceTransforms.push_back(xform);
			return inst.mParticleOffset;
		}
		//if (use_rb) UpdateInstanceTransform(mInstances.size() - 1, position, rotation);

	}

	void CreateInstanceChainFromData(Rope& rope, float* data,
		int nfloat, int phase, float mass = 1.0f,
		bool closed = true, int persistence = 2)
	{
		//stiffness = 1.0f;
		int start = int(g_buffers->positions.size());
		float r = 1.0f;//biased on the 1-3 spring
		int current = 0;
		//if closed do the last point ?
		int longP = 0;// 40;
		float longstiffness = 0;
		int count = 0;
		int pid0 = 0;
		int pid1 = 0;
		//cout << "nfloat " << nfloat << " extend " << extend << " " << extend_nb << endl;
		int ninstance = nfloat / 3;
		int prev_moff=0;
		rope.ropeInstanceIdStart = mInstancesFiber.size();

		for (int i = 0; i < nfloat - 3; i += 3)
		{
			//if (i/3>20) break;
			//printf("add a point %i %f %f %f \n", i / 3, data[i] * scale, data[i + 1] * scale, data[i + 2] * scale);
			int begin = int(g_buffers->positions.size());
			int prev = begin;//int(g_positions.size())-1;
			Vec3 p = Vec3(data[i],data[i + 1],data[i + 2]) * main_scale;
			Quat q = Quat(0,0,0,1);//rotation to align on the path
			if (count < ninstance -1){
				Vec3 np = Vec3(data[i + 3], data[i + 4], data[i + 5]) * main_scale;
				Vec3 along = np - p;
				q = AlignVec3s(along, Vec3(0, 0, 1));
			}
			else {
				Vec3 pp = Vec3(data[i - 3], data[i - 2], data[i - 1]) * main_scale;
				Vec3 along = p - pp;
				q = AlignVec3s(along, Vec3(0, 0, 1));
			}
			std::cout << " createInstanceIngredientFiber "  << endl;
			int off = g_buffers->positions.size();
			int moff = createInstanceIngredientFiber(rope.ropeType, p, q, mass, false);
			rope.mIndices.push_back(off);
			int pid1 = off;
			std::vector<int> beads1;
			std::vector<float> stiffness;
			std::vector<int> beads2;
			//attach to previous one
			if (count > 0){
				float cut = 34.0f;
				int found = 0;
				int cfound = 5;
				NvFlexExtAsset* asset = iBatchesFiber[rope.ropeType].mAsset;
				for (int k = 0; k < asset->numParticles - 1; k++) {
					found = 0;
					for (int l = k; l < asset->numParticles; l++) {
						if (found > cfound) continue;
						Vec3 p1 = Vec3(g_buffers->positions[prev_moff + k]);
						Vec3 p2 = Vec3(g_buffers->positions[moff + l]);
						float distance = Length(p1 - p2);	
						if (distance < cut * main_scale) {
							beads1.push_back(prev_moff + k);
							stiffness.push_back(1.0f);
							beads2.push_back(moff + l);
							found = found+1;
						}			
					}
				}
				for (int j = 0; j < beads1.size(); j++) {
					//std::cout << "create spring1 " << beads1[j] << " " << beads2[j] << endl;
					CreateSpringInter(beads1[j], beads2[j], stiffness[j], 0.0f, main_radius*2.0f);
				}
				//cout << "link pid " << pid0 << " to " << pid1 << endl;
				//CreateSpringInter(pid0, pid1, 1.0f, 0.0f, 0.0f);
			}
			//CreatePersistence(rope, begin, persistence, stiffness, nfloat, give, D);
			//CreatePersistenceDisc(rope, begin, longP, longstiffness, nfloat, give, D);//dna specfici test
			pid0 = pid1;
			prev_moff = moff;
			count++;
		}
		if (closed) {
			//cout << "nb poitns " << rope.mIndices.size() << endl;
			r = Randf(-1.0f, 1.0f) / 2000.0f;
			int startindex = rope.mIndices[0];
			int endindex = int(g_buffers->positions.size() - 1);// rope.mIndices[rope.mIndices.size() - 1];
			//CreateClosedPersistence(rope, startindex, endindex, persistence, stiffness, nfloat, give, D);
			//CreateClosedPersistenceDisc(rope, startindex, endindex, longP, longstiffness, nfloat, give, D);
		}
		rope.coarseIndices = rope.mIndices;
	}



	void loadRecipe(string filename, bool ignore_comp = false)
	{
		totalNBMol = 0;
		std::cout << filename << endl;
		std::ifstream ifs(filename);
		if (ifs.is_open())
		{
			ifs >> book_json;
			std::cout << "file opened and closed";
		}
		else
		{
			std::cout << "Error opening file";
		}
		Json::Value bbox = book_json["options"]["boundingBox"];
		minExtents = Vec3(bbox[0][0].asFloat(), bbox[0][1].asFloat(), bbox[0][2].asFloat())*main_scale;
		maxExtents = Vec3(bbox[1][0].asFloat(), bbox[1][1].asFloat(), bbox[1][2].asFloat())*main_scale;

		Json::Value  cyto = book_json["cytoplasme"];
		if (cyto != 0)
		{
			std::cout << "find compartments cyto " << cyto.size() << endl;
			Json::Value cyto_ingredients = cyto["ingredients"];
			std::cout << "compartment cyto should have n ingredients " << cyto_ingredients.size() << endl;
			if (cyto_ingredients != 0)
			{
				parseJsonIngredients(cyto_ingredients, 0);
			}
		}
		Json::Value  comp = book_json["compartments"];
		int i = 0;
		if (comp != 0)
		{
			std::cout << "find compartments " << comp.size() << endl;

			for (Json::ValueIterator itr = comp.begin(); itr != comp.end(); itr++)
			{
				//if (i==0){
				//skip first compartments for now
				//    i++;
				//    continue;
				// }
				Json::Value comp_name = *itr;
				//CompMask* cm = new CompMask();
				std::cout << "compartment should have several childs " << comp_name.size() << " " << comp_name["name"].asString() << endl;
				if (comp_name.size() == 0) continue;
				Json::Value comp_geom = comp_name["geom"];//at 0 ?
				std::cout << "compartmentMesh" << comp_geom.asString() << endl;
				if (comp_geom != 0 || comp_name["filename"] != 0 )
				{
					std::cout << "compartmentMesh" << endl;
					compartmentMesh(comp_name);
					//if (!ignore_comp)compartmentsSDF(comp_name);
				}
				Json::Value comp_surface = comp_name["surface"];//at 0 ?
				if (comp_surface != 0)
				{
					Json::Value surf_ingredients = comp_surface["ingredients"];
					std::cout << "compartment should have n ingredients " << surf_ingredients.size() << endl;
					if (surf_ingredients != 0)
					{
						parseJsonIngredients(surf_ingredients, i + 1);
						//cm->mask.push_back(i);
					}
				}
				if (comp_name.size() == 1) continue;
				Json::Value comp_interior = comp_name["interior"];//at 0 ?
				if (comp_interior != 0)
				{
					Json::Value int_ingredients = comp_interior["ingredients"];
					std::cout << "compartment should have n ingredients " << int_ingredients.size() << endl;
					if (int_ingredients != 0) {
						parseJsonIngredients(int_ingredients, -(i + 1));
						//cm->mask.push_back(-i);
					}
				}
				//comp_mask.push_back(cm);
				i++;
			}
		}
		iGroupCounter = 0;
		maxParticles = 1024 * 1024;
	}

	void loadIngredientFromCompartment(Json::Value comp, int compid) {
		std::cout << "load IngredientFromCompartment " << comp["name"].asString() << " " << compid << endl;
		if (comp["IngredientGroups"].size() != 0) {
			std::cout << "comp[IngredientGroups].size() " << comp["IngredientGroups"].size() << endl;
			Json::Value igroup = comp["IngredientGroups"][0];
			if (igroup["Ingredients"].size() != 0) {
				Json::Value ingredients = igroup["Ingredients"];
				std::cout << "compartment should have n ingredients " << ingredients.size() << endl;
				parseJsonIngredientsSerialized(ingredients, compid);
			}
		}
	}

	int loadOneCompartmentSerialized(Json::Value comp, int compid) {
		Json::Value comp_name = comp;
		//CompMask* cm = new CompMask();
		std::cout << "1-compartment should have several childs " << comp_name.size() << " " << comp_name["name"].asString() << " " << comp_name["geom_type"].asString() << " " << comp_name["filename"].asString() << " " << compid << endl;
		if (comp_name.size() == 0) return 0;
		Json::Value comp_geom = comp_name["geom_type"];//at 0 ?
		if (comp_geom.asString() != "None" || comp_name["filename"] != 0)
		{
			compartmentMesh(comp_name);
			//if (!ignore_comp)compartmentsSDF(comp_name);
		}
		//check if compartments child
		Json::Value  comp_childs = comp_name["Compartments"];
		int count_comp = 0;
		//need to load in order surface, interior, compartment
		Json::Value surface;
		Json::Value interior;
		for (int j = 0; j < comp_childs.size(); j++) {
			if (comp_childs[j]["name"].asString() == "surface") {
				surface = comp_childs[j];
				//loadIngredientFromCompartment(comp_childs[j], compid);
			}
			else if (comp_childs[j]["name"].asString() == "interior") {
				interior = comp_childs[j];
				//loadIngredientFromCompartment(comp_childs[j],-compid);
			}
			//else {
			//	count_comp++;
			//	count_comp+=loadOneCompartmentSerialized(comp_childs[j], compid + count_comp);
			//}
		}
		if (surface != NULL) loadIngredientFromCompartment(surface, compid);
		if (interior != NULL) loadIngredientFromCompartment(interior, -compid);
		for (int j = 0; j < comp_childs.size(); j++) {
			if (comp_childs[j]["name"].asString() == "surface") {
				continue;
			}
			else if (comp_childs[j]["name"].asString() == "interior") {
				continue;
			}
			else {
				count_comp++;
				count_comp += loadOneCompartmentSerialized(comp_childs[j], compid + count_comp);
			}
		}
		return count_comp;
	}

	void loadRecipeSerialized(string filename, bool ignore_comp = false)
	{
		totalNBMol = 0;
		std::cout << filename << endl;
		std::ifstream ifs;
		// Set exceptions to be thrown on failure
		ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		std::cout << "try to open" << endl;
		try {
			ifs.open(filename);
		}
		catch (std::system_error& e) {
			std::cerr << e.code().message() << std::endl;
		}
		std::cout << "is_open?" << ifs.is_open() << endl;
		if (ifs.is_open())
		{
			ifs >> book_json;
			std::cout << "file opened and closed";
		}
		else
		{
			std::cout << "Error opening file" << endl;
			return;
		}
		std::cout << "is_open-true" << endl;
		//Json::Value bbox = book_json["options"]["boundingBox"];
		//minExtents = Vec3(bbox[0][0].asFloat(), bbox[0][1].asFloat(), bbox[0][2].asFloat())*main_scale;
		//maxExtents = Vec3(bbox[1][0].asFloat(), bbox[1][1].asFloat(), bbox[1][2].asFloat())*main_scale;
		//book_json is root
		Json::Value  root = book_json;
		//check if any ingredient in the outisde
		loadIngredientFromCompartment(root,0);
		std::cout << "start loading" << endl;
		Json::Value  comp = book_json["Compartments"];
		int i = 0;
		int compid = 1;
		if (comp.size() != 0)
		{
			std::cout << "find compartments " << comp.size() << endl;
			for (int i=0;i< comp.size();i++) 
			{
				compid = compid + i;
				Json::Value comp_name = comp[i];
				//CompMask* cm = new CompMask();
				std::cout << compid << " compartment should have several childs " << comp_name.size() << " " << comp_name["name"].asString() << endl;
				if (comp_name.size() == 0) continue;
				Json::Value comp_geom = comp_name["geom_type"];//at 0 ?
				std::cout << "comp_geom " << comp_geom.asString() << endl;
				if (comp_geom.asString() != "None" || comp_name["filename"] != 0)
				{
					if (!ignore_comp) compartmentMesh(comp_name);
					//if (!ignore_comp)compartmentsSDF(comp_name);
				}
				//check if compartments child
				Json::Value  comp_childs = comp_name["Compartments"];
				Json::Value surface;
				Json::Value interior;
				for (int j = 0; j < comp_childs.size(); j++) {
					if (comp_childs[j]["name"].asString() == "surface") {
						surface = comp_childs[j];
					}
					else if (comp_childs[j]["name"].asString() == "interior") {
						interior = comp_childs[j];
					}
				}
				if (surface != NULL) loadIngredientFromCompartment(surface, compid);
				if (interior != NULL) loadIngredientFromCompartment(interior, -compid);
				for (int j = 0; j < comp_childs.size(); j++) {
					if (comp_childs[j]["name"].asString() == "surface") {
						continue;
					}
					else if (comp_childs[j]["name"].asString() == "interior") {
						continue;
					}
					else {
						int nchild = loadOneCompartmentSerialized(comp_childs[j], compid + 1);
						compid += nchild;
					}
				}
			}
		}
		iGroupCounter = 0;
		maxParticles = 1024 * 1024;
	}
	

	void finalize(){
		g_buffers->activeIndices.resize(g_buffers->positions.size());
		for (int i = 0; i < g_buffers->activeIndices.size(); i++)
			g_buffers->activeIndices[i] = i;
	}

	void loadResults(string filename, bool redo = false)
	{
		std::ifstream ifs(filename);
		if (ifs.is_open())
		{
			ifs >> results_json;
			std::cout << "file opened and closed";
		}
		else
		{
			std::cout << "Error opening file " + filename;
		}

		Json::Value  cyto = results_json["cytoplasme"];
		if (cyto != 0)
		{
			std::cout << "find compartments cyto " << cyto.size() << endl;
			Json::Value cyto_ingredients = cyto["ingredients"];
			std::cout << "compartment cyto should have n ingredients " << cyto_ingredients.size() << endl;
			if (cyto_ingredients != 0)
			{
				parseJsonIngredientsResults(cyto_ingredients, 0,redo);
			}
		}
		Json::Value  comp = results_json["compartments"];
		int i = 0;
		if (comp != 0)
		{
			std::cout << "find compartments " << comp.size() << endl;

			for (Json::ValueIterator itr = comp.begin(); itr != comp.end(); itr++)
			{
				Json::Value comp_name = *itr;
				//CompMask* cm = new CompMask();
				std::cout << "2-compartment should have several childs " << comp_name.size() << endl;
				if (comp_name.size() == 0) continue;
				Json::Value comp_surface = comp_name["surface"];//at 0 ?
				if (comp_surface != 0)
				{
					Json::Value surf_ingredients = comp_surface["ingredients"];
					std::cout << "3-compartment should have n ingredients " << surf_ingredients.size() << endl;
					if (surf_ingredients != 0)
					{
						parseJsonIngredientsResults(surf_ingredients, i + 1, redo);
						//cm->mask.push_back(i);
					}
				}
				if (comp_name.size() == 1) continue;
				Json::Value comp_interior = comp_name["interior"];//at 0 ?
				if (comp_interior != 0)
				{
					Json::Value int_ingredients = comp_interior["ingredients"];
					std::cout << "compartment should have n ingredients " << int_ingredients.size() << endl;
					if (int_ingredients != 0) {
						parseJsonIngredientsResults(int_ingredients, -(i + 1),redo);
						//cm->mask.push_back(-i);
					}
				}
				//comp_mask.push_back(cm);
				i++;
			}
		}
		//repair for soft body
		//CompactObjects();
	}


	void updateFromBinaryBuffer() {
		//need some clearing first ?
	}




	void loadResultsBinary(string filename, bool redo = false, bool ignore_comp = false)
	{
		std::ifstream ifs(filename, std::ios::binary);
		if (ifs.is_open())
		{
			std::cout << "file opened" << endl;
		}
		else
		{
			std::cout << "Error opening file " + filename << endl;
			return;
		}
		int ninstance;
		int ncurvepts;
		int ncurves;
		ifs.read(reinterpret_cast<char*>(&ninstance), sizeof(int));
		ifs.read(reinterpret_cast<char*>(&ncurvepts), sizeof(int));
		ifs.read(reinterpret_cast<char*>(&ncurves), sizeof(int));
		std::cout << "found ninstance "  << ninstance  << " ncurvepts " << ncurvepts << " ncurve " << ncurves << " for this many names "<< pnames.size() << endl; // ?00
		if (ninstance != 0)
		{
			Vec4* pos = new Vec4[ninstance];
			Vec4* rot = new Vec4[ninstance];
			ifs.read(reinterpret_cast<char*>(&pos[0]), sizeof(float) * 4 * ninstance);
			ifs.read(reinterpret_cast<char*>(&rot[0]), sizeof(float) * 4 * ninstance);
			for (int i=0;i<ninstance;i++){
				if (std::isnan(pos[i][0]) || std::isnan(pos[i][1]) || std::isnan(pos[i][2]) || std::isnan(pos[i][3])) {
					std::cout << "IS NAN instance " << i << " with pid " << pos[i][3] << " " << pnames[pos[i][3]] << endl;
				}
				Vec3 p = Vec3(-pos[i][0],pos[i][1],pos[i][2])*main_scale;
				Quat q = Quat(-rot[i][0],rot[i][1],rot[i][2],-rot[i][3]);
				int ingrIndex = (int) pos[i][3];
				if (ingrIndex >= pnames.size() || ingrIndex >= iBatches.size())
				{
					std::cout << "pnames " << pnames.size() << " " << iBatches.size() << endl;
					std::cout << "instance " << i << " with pid " << ingrIndex << " " << pnames[ingrIndex] << " " << iBatches[ingrIndex].ingr_name << endl;
				}
				int offset_id = createInstanceIngredient(ingrIndex, p, q);	
				// Instance inst1 = mInstances[i];
				// std::cout << ingrIndex << " instance " << i << " with pid " << ingrIndex << " " << pnames[ingrIndex] << " inst1.mMeshIndex " << inst1.mMeshIndex << endl;
			}
		}
		int curve_offset = g_buffers->positions.count;
		int curve_rb = 0;
		if (ncurves!=0 && ncurvepts != 0)
		{
			int persistence = 5;
			float D = g_params.radius;			
			Vec4* positions = new Vec4[ncurvepts];
			Vec4* normals = new Vec4[ncurvepts];
			Vec4* infos = new Vec4[ncurvepts];//curveId, count, startid, mask
			Vec4* curves = new Vec4[ncurves];//curveType,rs,rc,tu
			ifs.read(reinterpret_cast<char*>(&positions[0]), sizeof(float) * 4 * ncurvepts);
			ifs.read(reinterpret_cast<char*>(&normals[0]), sizeof(float) * 4 * ncurvepts);
			ifs.read(reinterpret_cast<char*>(&infos[0]), sizeof(float) * 4 * ncurvepts);
			ifs.read(reinterpret_cast<char*>(&curves[0]), sizeof(float) * 4 * ncurves);
			int npoints = 0;
			int current_curveId = (int) infos[0][0];
			int current_curveType = (int)curves[current_curveId][0]; //(int) infos[0][1];
			std::vector<float> points;
			for (int i=0;i<ncurvepts;i++){
				int curveId = (int) infos[i][0];
				int curveType = (int)curves[curveId][0]; //(int) infos[i][1];
				float L = infos[i][3];
				float mass = (infos[i][3] == -1.0f) ? 0.0f : 1.0f;
				// std::cout << "curveType " << curveType << " " << curveId << " infos ? " << infos[i][1] << " " << mass << " " << infos[i][3] << endl;
				if (curveId != current_curveId) {
					//get length and if its closed and fiber id
					//extend = false;
					std::cout << i << " curveType " << curveType << " curveId " << curveId << " current_curveType " << current_curveType << " current_curveId " << current_curveId << " pnames_fiber[current_curveType] " << pnames_fiber[current_curveType] << " size "<< points.size() << " " << main_radius << " " << mass << " " << infos[i][3] << endl;
					Json::Value fiber = getFiberNode(current_curveType);
					bool close = false;// fiber["closed"].asBool();
					L = fiber["length"].asFloat()*main_scale;
					L = 150.0f*main_scale;
					int rope_phase = NvFlexMakePhase(iGroupCounter++, eNvFlexPhaseSelfCollide);// eNvFlexPhaseFluid);// eNvFlexPhaseSelfCollide);
					Rope rcurve;
					//persistence should be a parameter per fiber
					//(aname.find("DNA") != std::string::npos)
					if (current_curveType == 0 && pnames_fiber[current_curveType].find("DNA") != std::string::npos) {
						//dna
						close = true;
						persistence = dna_persistence;
					}
					else {
						close = false;
						persistence = 1;
					}
					persistence = dna_persistence;
					rcurve.persistence = persistence;
					rcurve.ropeType = current_curveType;
					if (!testFiberName(pnames_fiber[current_curveType])){
						// create instance along control point
						rcurve.ropeOffsetIdStart = g_buffers->rigidTranslations.size();
						CreateInstanceChainFromData(rcurve,
											points.data(),points.size(),
											rope_phase,mass,close,persistence);
						curve_rb++;
					}
					else {
						rcurve.ropeOffsetIdStart = -1;
						rcurve.ropeInstanceIdStart = -1;
						CreateRopeFromData(rcurve, //rope
							main_scale, // scale
							1.0f, //stifness
							points.data(), //data
							L, //length
							points.size(),  //nfloat
							rope_phase,//phase
							0.0f,//spiral angle
							mass,//invmass
							0.0f,//give
							0,//extend_nb
							false,//extend
							close,//close
							main_radius*2.0f,
							persistence);
					}
					//printf("create Rope from data OK with %i points type %i\n", rcurve.mIndices.size(), current_curveType);
					//std::cout << " start at " << rcurve.mIndices[0] << " and finish at " << rcurve.mIndices[rcurve.mIndices.size()-1] << "  " << rcurve.mIndices.size() << endl;
					g_ropes.push_back(rcurve);//instance
					mask.push_back(-1);
					maks_protein.push_back(-(current_curveType + 1));//it ptype is 0 ?
					maks_fiber.push_back(current_curveType);//proteinType
					//printf("create Rope from data OK with %i  %i\n", -(current_curveType + 1), current_curveType);
					points.clear();
				}
				points.push_back(-positions[i][0]);
				points.push_back(positions[i][1]);
				points.push_back(positions[i][2]);
				current_curveId = curveId;
				current_curveType = curveType;
			};
			//do the last one
			Json::Value fiber = getFiberNode(current_curveType);
			float L = fiber["length"].asFloat()*main_scale; 
			std::cout << "last current_curveType " << current_curveType << " " << current_curveType << " " << L << endl;
			bool close = false;// fiber["closed"].asBool();
			
			int rope_phase = NvFlexMakePhase(iGroupCounter++, eNvFlexPhaseSelfCollide);// eNvFlexPhaseSelfCollideeNvFlexPhaseFluid);
			Rope rcurve;
			/*if (current_curveType == 0 && pnames_fiber[current_curveType].find("DNA") != std::string::npos) {
				close = true;
				persistence = dna_persistence;
			}
			else {
				close = false;
				persistence = 1;
			}*/
			persistence = dna_persistence;
			rcurve.persistence = persistence;
			rcurve.ropeType = current_curveType;
			float m = (infos[ncurvepts-2][3] == -1.0f) ? 0.0f : 1.0f;
			if (!testFiberName(pnames_fiber[current_curveType])){
				std::cout << "testFiberName " << endl;
				rcurve.ropeOffsetIdStart = g_buffers->rigidTranslations.size();
				//create instance along control point
				CreateInstanceChainFromData(rcurve,
									points.data(),points.size(),
									rope_phase,m,close,persistence);
			}
			else {
				rcurve.ropeOffsetIdStart = -1;
				rcurve.ropeInstanceIdStart = -1;
				CreateRopeFromData(rcurve, //rope
					main_scale, // scale
					1.0f, //stifness
					points.data(), //data
					L, //length
					points.size(),  //nfloat
					rope_phase,//phase
					0.0f,//spiral angle
					m,//invmass
					0.0f,//give
					0,//extend_nb
					false,//extend
					close,//close
					main_radius*2.0f,
					persistence);
			}
			
			//printf("create Rope from data OK with %i points\n", rcurve.mIndices.size());
			g_ropes.push_back(rcurve);//instance
			mask.push_back(-1);
			maks_protein.push_back(-(current_curveType + 1));//it ptype is 0 ?
			maks_fiber.push_back(current_curveType);//proteinType
			//printf("create Rope from data OK with %i  %i\n", -(current_curveType + 1), current_curveType);
		}
		//list of bound object
		//return;//debug for jitter
		//ifs.close(); return;
		if (!ifs.eof()) {
			//Id1 Id2ignore_comp
			//if id is negative, curve instance Id
			int nlinks;
			ifs.read(reinterpret_cast<char*>(&nlinks), sizeof(int));
			std::cout << "found nlinks " << nlinks/4  << endl;
			if (nlinks > 0)
			{
				int* links = new int[nlinks];
				ifs.read(reinterpret_cast<char*>(&links[0]), sizeof(int) * nlinks);
				//loop over the linkgs and create springs
				//mInstances and 
				//nlinks = 0;

				for (int i = 0; i < nlinks / 4; i++) {
					//Type ID Type ID
					int type1 = links[i * 4];
					int id1 = links[i * 4 + 1];
					int type2 = links[i * 4 + 2];
					int id2 = links[i * 4 + 3];
					int pid1 = 0;
					int pid2 = 0;
					std::vector<int> beads1;
					std::vector<float> stiffness;
					std::vector<int> beads2;
					std::cout << curve_offset << " found nlinks " << i << " " << type1 << " " << id1 << " " << type2 <<" "<< id2 << " " << mInstances.size() << endl;
					if (id1 < 0 && id2 >= 0) //control point
					{
						pid1 = curve_offset + (-id1);
						Instance inst2 = mInstances[id2];
						Json::Value ingr_node = getProteinNode(inst2.mMeshIndex);// getIngredientsSerialized(pnames[inst2.mMeshIndex]);
						//get the beads
					}
					else if (id1 >=0 && id2 < 0) {
						pid2 = curve_offset + (-id2) -1;
						Instance& inst1 = mInstances[id1];
						pid1 = inst1.mParticleOffset;
						string name = pnames[inst1.mMeshIndex];
						std::cout << pid1 <<" x "<< id1 <<" found mMeshIndex " << inst1.mMeshIndex << " " << name <<" for type1 " << type1 << " " << pid2 << " " << type2 << endl;
						
						int fiber_type = type2;// getFiberIdFromPointId(pid2);
						
						string fibname = pnames_fiber[fiber_type];
						//std::cout << "found partner id? " << fiber_type << " " << fibname << endl;
						if (inst1.mMeshIndex != type1) {
							std::cout << i << " pb with id1 " << id1 << " mMeshIndex " << inst1.mMeshIndex << " type1 " << type1 << " " << name << " id2 " << id2 << " type2 " << type2 << " " << fibname << endl;
						}
						//std::cout << "found fibname " << fibname << endl;
						//if (name == "NC_capsid")continue;

						Json::Value ingr_node = getProteinNode(type1);// inst1.mMeshIndex); //getIngredientsSerialized(name);
						
						if (ingr_node == NULL) {
							std::cout << "not found " << name << endl;
							continue;
						}
						//std::cout << "found node? " << ingr_node["name"].asString() << endl;
						//partner id ?
						int partner_id = getPartnerId(fibname, ingr_node) + inst1.bounded;
						std::cout << "partner_id " << partner_id << " found partner name? " << fibname << " " << ingr_node["partners_properties"].size() <<" " << use_partners_properties << " " << inst1.bounded << endl;
						if (ingr_node["partners_properties"].size() != 0 && use_partners_properties){
							if (ingr_node["partners_properties"][0]["binding_site_lod"].size() == 2){
								Json::Value b1 = ingr_node["partners_properties"][partner_id]["binding_site_lod"][0]["binding_site"];
								Json::Value b2 = ingr_node["partners_properties"][partner_id]["binding_site_lod"][1]["binding_site"];
								std::cout << "b1.size " << b1.size() << " b2.size " << b2.size() << endl;
								for (int j = 0; j < b1.size(); j++)
								{
									beads1.push_back(pid1+b1[j].asInt());
									stiffness.push_back(1.0f);
								}
								for (int j = 0; j < b2.size(); j++)
								{
									beads2.push_back(pid2+b2[j].asInt());	
								}		
								inst1.bounded = 1;
							}
							std::cout << "size properties " << ingr_node["partners_properties"][0]["binding_site_lod"].size() << " " << inst1.bounded << endl;

						}
						if (beads1.size() == 0 && beads2.size() == 0) {
							//get the closest from inst1.mParticleOffset-> particle count
							//std::cout << "1-automatic binding using closest bead " << main_radius * 2.0f << endl;
							NvFlexExtAsset* asset = iBatches[inst1.mMeshIndex].mAsset;
							float cut = main_radius * 2.0f; // 34.0f;
							int cfound = 1;
							int found = 0;
							float miniD = 9999999.9f;
							int mini_pid = inst1.mParticleOffset;
							for (int j = 0; j < asset->numParticles; ++j)
							{
								if (found > cfound) break;
								Vec3 p = Vec3(g_buffers->positions[inst1.mParticleOffset + j]);
								float distance = Length(Vec3(g_buffers->positions[inst1.mParticleOffset + j]) - Vec3(g_buffers->positions[pid2]));
								//use mini distance
								if (distance < miniD) {
									miniD = distance;
									mini_pid = inst1.mParticleOffset + j;
								}
								/*if (distance < 500.0f * main_scale) {
									// std::cout << "distance is " << distance << " cutoff " << cut * main_scale << endl;
									mini_pid = inst1.mParticleOffset + j;
									beads1.push_back(mini_pid);
									stiffness.push_back(1.0f);
									beads2.push_back(pid2);
									found = found + 1;
									//change the phase
								}*/
							}
							beads1.push_back(mini_pid);
							stiffness.push_back(1.0f);
							beads2.push_back(pid2);
							// g_buffers->phases[pid2] = g_buffers->phases[inst1.mParticleOffset];
						}
						//ing.partners_properties[0].binding_site_lod[0].binding_site.Count
					}
					else if (id1 >=0 && id2 >= 0) {
						Instance inst1 = mInstances[id1];
						NvFlexExtAsset* asset1 = iBatches[inst1.mMeshIndex].mAsset;
						Json::Value ingr_node1 = getProteinNode(inst1.mMeshIndex); //getIngredientsSerialized(pnames[inst1.mMeshIndex]);
						Instance inst2 = mInstances[id2];
						NvFlexExtAsset* asset2 = iBatches[inst2.mMeshIndex].mAsset;
						Json::Value ingr_node2 = getProteinNode(inst2.mMeshIndex); //getIngredientsSerialized(pnames[inst2.mMeshIndex]);
						std::cout << "2-found inst1 " << pnames[inst1.mMeshIndex] << " inst2 " << pnames[inst2.mMeshIndex] << endl;
						//closest pair with cutoff
						float cut = main_radius * 2.0f; //34.0f;
						int cfound = 2;
						int found = 0;
						for (int k = 0; k < asset1->numParticles - 1; k++) {
							found = 0;
							for (int l = k; l < asset2->numParticles; l++) {
								if (found > cfound) continue;
								Vec3 p1 = Vec3(g_buffers->positions[inst1.mParticleOffset + k]);
								Vec3 p2 = Vec3(g_buffers->positions[inst2.mParticleOffset + l]);
								float distance = Length(p1 - p2);
								// std::cout << " distance "<< k <<" " << l << " " << distance << " " << (cut * main_scale) << endl;
								if (distance < cut * main_scale) {
									beads1.push_back(inst1.mParticleOffset + k);
									stiffness.push_back(1.0f);
									beads2.push_back(inst2.mParticleOffset + l);
									found = found+1;
								}
								// g_buffers->phases[inst2.mParticleOffset + l] = g_buffers->phases[inst1.mParticleOffset + k];
							}
						}
					}
					else if (id1 < 0 && id2 < 0) {
						pid1 = curve_offset + (-id1);
						pid2 = curve_offset + (-id2);
						//tether ?
						std::cout << "inter curve " << pid1 << " " << pid2 << endl;
						CreateSpringInter(pid1, pid2, 1.0f, 0.0f, main_radius*2.0f);
						//beads1.push_back(pid1);
						//beads2.push_back(pid2);
					}
					//std::cout << "found pair pid1 " << pid1 << " pid2 " << pid2 << endl;
					if (beads1.size() == beads2.size()){
						//1-1 links
						//std::cout << "found beads equal size " << beads1.size() << " " << beads2.size() << endl;
						for (int j = 0; j < beads1.size(); j++) {
							//std::cout << "create spring1 " << beads1[j] << " " << beads2[j] << endl;
							CreateSpringInter(beads1[j], beads2[j], stiffness[j], 0.0f, main_radius*2.0f);
						}
					} else {
						//1-all
						//std::cout << "found beads of size " << beads1.size() << " " << beads2.size() << endl;
						for (int j = 0; j < beads1.size(); j++) 
						{
							for (int k = 0; k < beads2.size(); k++) {
								//std::cout << "create spring2 " << beads1[j] << " " << beads2[k] << endl;
								CreateSpringInter(beads1[j], beads2[k], stiffness[j], 0.0f, main_radius*2.0f);
							}
						}
					}
					//if (i > 77) break;
				}
			}
			if (!ifs.eof()) {
				std::cout << "search ncomps " << endl;
				int ncomps;
				ifs.read(reinterpret_cast<char*>(&ncomps), sizeof(int));
				std::cout << "found ncomps " << ncomps << endl;
				if (ncomps > 0 && ncomps < 1000) {
					std::cout << "found ncomps " << ncomps << endl;
					Vec4* comps = new Vec4[ncomps];
					ifs.read(reinterpret_cast<char*>(&comps[0]), sizeof(float) * 4 * ncomps);
					//use first one for now
					comp_radius = comps[0][3] * main_scale;
					std::cout << "found radius compartment " << comps[0][3] << " " << comp_radius << endl;
					//update the compartment mesh
					//OneCompartmentMesh();
				}
				if (!ifs.eof()) {
					int nmesh;
					ifs.read(reinterpret_cast<char*>(&nmesh), sizeof(int));
					std::cout << "found nmesh " << nmesh << endl;
					if (nmesh > 0 && !ignore_comp) {
						// std::cout << "found nmesh " << nmesh << endl;
						int nv;
						int nf;
						ifs.read(reinterpret_cast<char*>(&nv), sizeof(int));
						ifs.read(reinterpret_cast<char*>(&nf), sizeof(int));
						Vec3* vertices = new Vec3[nv];
						Vec3* normals = new Vec3[nv];
						int* faces = new int[nf];
						ifs.read(reinterpret_cast<char*>(&vertices[0]), sizeof(float) * 3 * nv);
						ifs.read(reinterpret_cast<char*>(&normals[0]), sizeof(float) * 3 * nv);
						ifs.read(reinterpret_cast<char*>(&faces[0]), sizeof(int) * nf);
						//rebuild the mesh !
						int numVertices = nv;
						int numFaces = nf;
						Mesh* mesh = createMesh(nv, nf, vertices, normals, faces);
						//center 
						Vec3 lower = Vec3(0.0f);
						Vec3 upper = Vec3(1.0f);
						
						Vector3 minExtents, maxExtents;
						mesh->GetBounds(minExtents, maxExtents);
						Vector3 originalCenter = (minExtents + maxExtents) * 0.5f;
						std::cout << " originalCenter " << originalCenter[0] << " " << originalCenter[1] << " " << originalCenter[2] << endl;
						mesh->Transform(ScaleMatrix(main_scale));
						
						// Re-calculate the center after scaling (assuming the mesh might be modified by Transform)
						Vector3 newMinExtents, newMaxExtents;
						mesh->GetBounds(newMinExtents, newMaxExtents);
						Vector3 newCenter = (newMinExtents + newMaxExtents) * 0.5f;
						std::cout << " newCenter " << newCenter[0] << " " << newCenter[1] << " " << newCenter[2] << endl;
						// Calculate the translation required to re-center the mesh
						//Vector3 centerTranslation = originalCenter - newCenter;
						// mesh->Transform(TranslationMatrix(Point3(centerTranslation)));
						//NvFlexDistanceFieldId sdf = CreateSDFFromMeshF(GetFilePathByPlatform("../../data/cellpack.ply").c_str(), mesh, 128);
						//AddSDF(sdf, (upper - lower)*0.5f, Quat(), 1.0f);
						duplicateMeshMembrane(nv, nf, vertices, normals, faces);
						//mesh->Transform(TranslationMatrix(Point3((upper - lower))*0.5f));
						//NvFlexTriangleMeshId upmesh = CreateTriangleMesh(mesh);
						//AddTriangleMesh(upmesh, Vec3(), Quat(), 1.0f);
						//comp_tri.push_back(mesh);
					}
				}
			}
		}
		ifs.close();
		//force recenter vRNP
		/*for (int p = 0; p < mInstances.size(); p++) {
			Instance inst1 = mInstances[p];
			string astr; astr.assign(iBatches[inst1.mMeshIndex].ingr_name);
			std::cout << " aname is " << iBatches[inst1.mMeshIndex].ingr_name <<" astr " << astr  << " "<< inst1.mMeshIndex << endl;
			if (inst1.mMeshIndex == 11 || inst1.mMeshIndex == 12) {
				UpdateInstanceTransform(p, Vec3(0.0f), inst1.mRotation);
			}
		} //and the rna ? 
		for (int f = 0; f < g_ropes.size(); f++) {
			Rope r = g_ropes[f];
			if (r.ropeType == 0) {
				for (int p = 0; p < r.mIndices.size(); p++) {
					g_buffers->positions[r.mIndices[p]] = Vec4(0.0f,0.0f,0.0f, 1.0f);
				}
			}
		}*/
		std::cout << "loading done!" << endl;
	}

	void setupCellPACK(){
		use_rb = true;
		main_scale = 1.0f / 100.0f;
		main_radius = 11.85f*main_scale;//((10.0f/100.0f)*(10.0f/100.0f));
		//mycoplasma experiment is
		loadRecipe((mainpath + "recipes\\Mycoplasma1.6_full.json").c_str());
		//loadResults((mainpath + "recipes\\Mycoplasma1.5_mixed_pdb_fixed.json").c_str());
		haltondistribute(100);
		compactInstances();
		placeFibers(false);
	}

	void DrawMeshInstance(int pass)
	{
		if (!g_drawMesh)
			return;

		for (int b = 0; b < int(iBatches.size()); ++b)
		{
			if (iBatches[b].nvMesh == NULL) continue;
			if (iBatches[b].mInstanceTransforms.size())
			{
				extern Colour g_colors[];
				DrawGpuMeshInstances(g_meshes[iBatches[b].nvMesh], &iBatches[b].mInstanceTransforms[0], iBatches[b].mInstanceTransforms.size(), Vec3(g_colors[b % 8]));
			}
		}
		//what about triangle collision ?
	}
	
	bool oneInstance(int index){
		int particleOffset = NvFlexGetActiveCount(g_solver);
		NvFlexExtAsset* asset = iBatches[index].mAsset;
		// check we can fit in the container
		if (int(g_buffers->positions.size()) - particleOffset < asset->numParticles)
			return false;
		Instance inst;
		inst.mLifetime = 1000.0f;// Randf(5.0f, 60.0f);
		inst.mParticleOffset = particleOffset;
		inst.mRotation = QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);

		float spread = 0.2f;
		inst.mTranslation = g_emitters[0].mPos + Vec3(Randf(-spread, spread), Randf(-spread, spread), 0.0f);
		inst.mMeshIndex = index;
		Vec3 linearVelocity = g_emitters[0].mDir*15.0f;//*Randf(5.0f, 10.0f);//Vec3(Randf(0.0f, 10.0f), 0.0f, 0.0f);
		Vec3 angularVelocity = Vec3(UniformSampleSphere()*Randf()*k2Pi);

		inst.mGroup = mGroupCounter++;

		const int phase = NvFlexMakePhase(inst.mGroup, 0);

		// generate initial particle positions
		for (int j = 0; j < asset->numParticles; ++j)
		{
			Vec3 localPos = Vec3(&asset->particles[j * 4]);// -Vec3(&asset->shapeCenters[0]);

			g_buffers->positions[inst.mParticleOffset + j] = Vec4(inst.mTranslation + inst.mRotation*localPos, 1.0f);
			g_buffers->velocities[inst.mParticleOffset + j] = linearVelocity + Cross(angularVelocity, localPos);
			g_buffers->phases[inst.mParticleOffset + j] = phase;
		}

		particleOffset += asset->numParticles;

		mInstances.push_back(inst);
		return true;
	}

	int oneInstanceAt(int ingr_index, int instance_index, Vec3 pos, Quat quat, int particleOffset){
		//can only be called after initialization
		cout << "ingr " << ingr_index << endl;
		//int particleOffset = NvFlexGetActiveCount(g_solver);
		cout << "offset " << particleOffset << endl;
		NvFlexExtAsset* asset = iBatches[ingr_index].mAsset;
		// check we can fit in the container
		
		Instance inst;
		inst.mLifetime = 1000.0f;// Randf(5.0f, 60.0f);
		inst.mParticleOffset = particleOffset;
		inst.mRotation = quat;// QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);

		float spread = 0.2f;
		inst.mTranslation = pos;// g_emitters[0].mPos + Vec3(Randf(-spread, spread), Randf(-spread, spread), 0.0f);
		inst.mMeshIndex = ingr_index;
		Vec3 linearVelocity = Vec3();// g_emitters[0].mDir*15.0f;//*Randf(5.0f, 10.0f);//Vec3(Randf(0.0f, 10.0f), 0.0f, 0.0f);
		Vec3 angularVelocity = Vec3();// UniformSampleSphere()*Randf()*k2Pi);

		inst.mGroup = iGroupCounter++;

		const int phase = NvFlexMakePhase(iGroupCounter++, 0);

		// generate initial particle positions
		for (int j = 0; j < asset->numParticles; ++j)
		{
			Vec3 localPos = Vec3(&asset->particles[j * 4])  - Vec3(&asset->shapeCenters[0]);

			g_buffers->positions[inst.mParticleOffset + j] = Vec4(inst.mTranslation + inst.mRotation*localPos, 1.0f);
			g_buffers->velocities[inst.mParticleOffset + j] = Vec3(0.0f);
			g_buffers->phases[inst.mParticleOffset + j] = phase;
		}
		particleOffset += asset->numParticles;
		mInstances[instance_index] = inst;
		iBatches[ingr_index].nInstances++;
		return particleOffset;
	}

	void compactInstances(){
		// compact instances 
		static std::vector<Vec4> particles(g_buffers->positions.size());
		static std::vector<Vec3> velocities(g_buffers->velocities.size());
		static std::vector<int> phases(g_buffers->phases.size());
		
		std::vector<Vec4>::iterator itpos;
		std::vector<Vec3>::iterator itvel;
		std::vector<int>::iterator itindices;

		g_buffers->rigidTranslations.resize(0);
		g_buffers->rigidRotations.resize(0);
		g_buffers->rigidCoefficients.resize(0);
		g_buffers->rigidIndices.resize(0);
		g_buffers->rigidLocalPositions.resize(0);
		g_buffers->rigidOffsets.resize(0);

		// start index
		g_buffers->rigidOffsets.push_back(0);

		// clear mesh batches
		for (int i = 0; i < int(iBatches.size()); ++i)
			iBatches[i].mInstanceTransforms.resize(0);
		int rsize = 0;
		for (int i = 0; i < g_ropes.size(); i++)
		{
			rsize += g_ropes[i].mIndices.size();
		}
		g_buffers->positions.copyto(&particles[0], rsize);
		g_buffers->velocities.copyto(&velocities[0], rsize);
		g_buffers->phases.copyto(&phases[0], rsize);

		numActive = rsize;
		//treat the rope data
		/*for (int i = 0; i < g_ropes.size(); i++)
		{
			itindices = g_ropes[i].mIndices.begin();


			cout << " numActive " << numActive << " " << g_ropes[i].mIndices.size() << " " << g_buffers->positions.size() << " " << particles.size() << " " << velocities.size() << endl;
			//particles.assign()
			//for (int j = 0; j < g_ropes[i].mIndices.size(); j++)
			//{
				//if (g_ropes[i].mIndices[j] > g_buffers->positions.size()) cout << g_ropes[i].mIndices[j] << endl;
			//	particles[numActive + j] = g_buffers->positions[g_ropes[i].mIndices[j]];
			//	velocities[numActive + j] = Vec3(0.0f);// g_buffers->velocities[inst.mParticleOffset + j];
			//	phases[numActive + j] = g_buffers->phases[g_ropes[i].mIndices[j]];
			//}
			numActive += g_ropes[i].mIndices.size();
			//cout << " numActive " << numActive << endl;
		}
		*/
		for (int i = 0; i < int(mInstances.size()); ++i)
		{
			Instance& inst = mInstances[i];
			if (inst.mMeshIndex < 0) continue;
			//cout << " compact instance of " << inst.mMeshIndex << " size batch " << iBatches.size() << endl;
			NvFlexExtAsset* asset = iBatches[inst.mMeshIndex].mAsset;

			for (int j = 0; j < asset->numParticles; ++j)
			{
				particles[numActive + j] = g_buffers->positions[inst.mParticleOffset + j];
				velocities[numActive + j] = Vec3(0.0f);// g_buffers->velocities[inst.mParticleOffset + j];
				phases[numActive + j] = g_buffers->phases[inst.mParticleOffset + j];
			}

			g_buffers->rigidCoefficients.push_back(1.0f);
			g_buffers->rigidTranslations.push_back(inst.mTranslation);
			g_buffers->rigidRotations.push_back(inst.mRotation);

			for (int j = 0; j < asset->numShapeIndices; ++j)
			{
				g_buffers->rigidLocalPositions.push_back(Vec3(&asset->particles[j * 4]) - Vec3(&asset->shapeCenters[0]));
				g_buffers->rigidIndices.push_back(asset->shapeIndices[j] + numActive);
			}

			g_buffers->rigidOffsets.push_back(g_buffers->rigidIndices.size());

			//mInstances[i].mParticleOffset = numActive;

			// Draw transform
			Matrix44 xform = TranslationMatrix(Point3(inst.mTranslation - inst.mRotation*Vec3(asset->shapeCenters)))*RotationMatrix(inst.mRotation);
			iBatches[inst.mMeshIndex].mInstanceTransforms.push_back(xform);

			numActive += asset->numParticles;
		}

		// update particle buffers
		g_buffers->positions.assign(&particles[0], particles.size());
		g_buffers->velocities.assign(&velocities[0], velocities.size());
		g_buffers->phases.assign(&phases[0], phases.size());

		// rebuild active indices
		g_buffers->activeIndices.resize(numActive);
		for (int i = 0; i < numActive; ++i)
			g_buffers->activeIndices[i] = i;
		if (g_buffers->rigidOffsets.size())
		{
			assert(g_buffers->rigidOffsets.size() > 1);

			const int numRigids = g_buffers->rigidOffsets.size() - 1;
			// If the centers of mass for the rigids are not yet computed, this is done here
			// (If the CreateParticleShape method is used instead of the NvFlexExt methods, the centers of mass will be calculated here)
			if (g_buffers->rigidTranslations.size() == 0)
			{
				g_buffers->rigidTranslations.resize(g_buffers->rigidOffsets.size() - 1, Vec3());
				CalculateRigidCentersOfMass(&g_buffers->positions[0], 
					g_buffers->positions.size(), 
					&g_buffers->rigidOffsets[0], 
					&g_buffers->rigidTranslations[0], 
					&g_buffers->rigidIndices[0], numRigids);
			}

			// calculate local rest space positions
			g_buffers->rigidLocalPositions.resize(g_buffers->rigidOffsets.back());
			CalculateRigidLocalPositions(&g_buffers->positions[0], 
				&g_buffers->rigidOffsets[0], 
				&g_buffers->rigidTranslations[0], 
				&g_buffers->rigidIndices[0], 
				numRigids, 
				&g_buffers->rigidLocalPositions[0]);
			
			// set rigidRotations to correct length, probably NULL up until here
			g_buffers->rigidRotations.resize(g_buffers->rigidOffsets.size() - 1, Quat());
			//g_buffers->rigidRotations.resize(g_buffers->rigidOffsets.size() - 1, Quat());
			//g_buffers->rigidTranslations.resize(g_buffers->rigidOffsets.size() - 1, Vec3());

		}
	}


	void nbMolhaltondistribute(){
		//distribute until everything placed
		if (iBatches.size() == 0) return;
		Vec3 top = maxExtents;
		Vec3 bot = minExtents;
		std::vector<Vec3> halton_positions(totalNBMol);
		Vec3 adim = (top - bot)*0.75f;//Vec3(5, 5, 5);
		Vec3 center = (-bot)*0.75f;//dim/2.0f;
		float scale_dim[3] = { adim.x, adim.y, adim.z };//scale on x y z should be the bounding box size
		//int n = PoissonSample3D(2.45f, radius*0.42f, &positions[0], positions.size(), 2);
		int n = HaltonSample3D(scale_dim, main_radius, &halton_positions[0], halton_positions.size());
		int nToPlace = totalNBMol;
		//for each of this position place a solid object ?
		int particleOffset = g_buffers->positions.size();
		int count = 0;
		int safety_count = 0;
		int big_safety_count = 0;
		int hpos = 0;

		while (nToPlace > 0){
			if (big_safety_count > 200000) break;
			int ingrIndex = Rand() % iBatches.size();
			
			Json::Value ingr_node = getIngredients(iBatches[ingrIndex].ingr_name);
			//test n toward nbmol
			if (iBatches[ingrIndex].nInstances >= ingr_node["nbMol"].asInt()){
				big_safety_count++;
				continue;
			}

			int compId = iBatches[ingrIndex].compId;
			Vec3 p = Vec3(halton_positions[hpos].x,
				halton_positions[hpos].y,
				halton_positions[hpos].z) - center;
			Quat r = QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);
			if (compId < 0) //inside 
			{
				if (!isInside(p, 0))
				{
					safety_count++;
					hpos++;
					if (safety_count > 10000)
						break;
					if (hpos > totalNBMol)
						break;
					continue;
				}
				else {
					safety_count = 0;
				}
			}
			else if (compId > 0)//surface
			{
				//pick a surface points
				int ncomp = comp_tri.size();
				int nVertices = comp_tri[compId - 1]->GetNumVertices();
				int vIndex = Rand() % nVertices;
				p = Vec3(comp_tri[compId - 1]->m_positions[vIndex].x, comp_tri[compId - 1]->m_positions[vIndex].y, comp_tri[compId - 1]->m_positions[vIndex].z);
				Vec3 normal = Vec3(comp_tri[compId - 1]->m_normals[vIndex].x, comp_tri[compId - 1]->m_normals[vIndex].y, comp_tri[compId - 1]->m_normals[vIndex].z);
				//align pcpalVector to normal
				Vec3 ingrpcpal = Vec3(iBatches[ingrIndex].pcpalVectorx, iBatches[ingrIndex].pcpalVectory, iBatches[ingrIndex].pcpalVectorz);
				Vec3 offsetPos = Vec3(iBatches[ingrIndex].offsetx, iBatches[ingrIndex].offsety, iBatches[ingrIndex].offsetz);
				//align pcpal to up
				r = AlignVec3s(normal, ingrpcpal);
				p = p + r*offsetPos;		
			}
			particleOffset = createInstanceIngredient(ingrIndex, p, r, true, 1.0f);
			count++;
			hpos++;
			nToPlace--;
		}

	}

	void haltondistribute(int N)
	{
		//how to insure inside compartments
		if (iBatches.size() == 0) return;
		Vec3 top = maxExtents;
		Vec3 bot = minExtents;
		std::vector<Vec3> halton_positions(N);
		Vec3 adim = (top - bot)*0.75f;//Vec3(5, 5, 5);
		Vec3 center = (-bot)*0.75f;//dim/2.0f;
		float scale_dim[3] = { adim.x, adim.y, adim.z };//scale on x y z should be the bounding box size
		//int n = PoissonSample3D(2.45f, radius*0.42f, &positions[0], positions.size(), 2);
		int n = HaltonSample3D(scale_dim, main_radius, &halton_positions[0], halton_positions.size());

		//for each of this position place a solid object ?
		int particleOffset = g_buffers->positions.size();
		int count = 0;
		int safety_count = 0;
		int hpos = 0;
		
		while (count < N){
			float mass = 1.0f;
			int ingrIndex = Rand() % iBatches.size();
			int compId = iBatches[ingrIndex].compId;
			Vec3 p = Vec3(halton_positions[hpos].x, 
						  halton_positions[hpos].y, 
						  halton_positions[hpos].z) - center;
			Quat r = QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);
			if (compId < 0) //inside 
			{
				/*//find inside point using sdf
				Vec3 lower = bot;
				Vec3 upper = top;
				//comp_tri[0]->GetBounds(lower, upper);
				
				float spacing = (upper[0] - lower[0]) / (float)dim;
				Vec3 botsdf = Vec3(lower[0], lower[1], lower[2]);
				Vec3 topsdf = Vec3(upper[0], upper[1], upper[2]);
				Vec3 pquery = Vec3((p.x - botsdf.x) / spacing, (p.y - botsdf.y) / spacing, (p.z - botsdf.z) / spacing);
				Vec3 pquery_index = Vec3((int)pquery.x, (int)pquery.y, (int)pquery.z);
				float D = SampleSDFX(sdfdata.m_data, dim, (int)pquery.x, (int)pquery.y, (int)pquery.z);
				//continue;
				cout << " D is " << D << endl;
				cout << isInside(p, 0) << endl;
				*/
				if (!isInside(p, 0))
				{
					safety_count++;
					hpos++;
					if (safety_count > 10000)
						break;
					if (hpos > N)
						break;
					continue;
				}
				else {
					safety_count = 0;
				}
			}
			else if (compId > 0)//surface
			{
				//pick a surface points
				int ncomp = comp_tri.size();
				int nVertices = comp_tri[compId - 1]->GetNumVertices();
				int vIndex = Rand() % nVertices;
				p = Vec3(comp_tri[compId - 1]->m_positions[vIndex].x, comp_tri[compId - 1]->m_positions[vIndex].y, comp_tri[compId - 1]->m_positions[vIndex].z);
				Vec3 normal = Vec3(comp_tri[compId - 1]->m_normals[vIndex].x, comp_tri[compId - 1]->m_normals[vIndex].y, comp_tri[compId - 1]->m_normals[vIndex].z);
				//align pcpalVector to normal
				Vec3 ingrpcpal = Vec3(iBatches[ingrIndex].pcpalVectorx, iBatches[ingrIndex].pcpalVectory, iBatches[ingrIndex].pcpalVectorz);
				Vec3 offsetPos = Vec3(iBatches[ingrIndex].offsetx, iBatches[ingrIndex].offsety, iBatches[ingrIndex].offsetz);
				//align pcpal to up
				r = AlignVec3s(normal, ingrpcpal);
				p = p + r*offsetPos;
				mass = 0.0f;//fix but shouldnt do it like this, better to used the membrane point
			}
			particleOffset = createInstanceIngredient(ingrIndex, p, r, true, mass);
			count++;
			hpos++;
		}
	}

	void randomDistribute(int N)
	{
		if (iBatches.size() == 0) return;
		Vec3 top = maxExtents;
		Vec3 bot = minExtents;
		std::vector<Vec3> halton_positions(N);
		Vec3 dim = top - bot;//Vec3(5, 5, 5);
		Vec3 center = -bot;//dim/2.0f;
		float scale_dim[3] = { dim.x, dim.y, dim.z };//scale on x y z should be the bounding box size
		//int n = PoissonSample3D(2.45f, radius*0.42f, &positions[0], positions.size(), 2);
		int n = HaltonSample3D(scale_dim, main_radius, &halton_positions[0], halton_positions.size());

		//for each of this position place a solid object ?
		for (int i = 0; i<N; i++)
		{
			const int ingrIndex = Rand() % iBatches.size();
			Quat q = Quat(0, 0, 0, 1);// QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);
			createInstanceIngredient(ingrIndex, 
				Vec3(halton_positions[i].x, halton_positions[i].y, halton_positions[i].z) - center, q);
		}
	}

	/*
	void writeToSharedMemoryRB(int n, int start) {
		for (size_t i = 0; i < shared_buffer_size_proteins; i++)
		{
			if ((i + start) >= mInstances.size())
				return;

		}
	}
	void writeToSharedMemoryCurve() {
	}

	void writeToSharedMemory() {
		//update info
		int nInst = mInstances.size();
		sharedMemInfo[0] = (float)nInst; //g_buffers->rigidTranslations.size();
		sharedMemInfo[1] = (float)g_ropes.size();
		sharedMemInfo[2] = 0;
		sharedMemInfo[3] = 0;
		NvFlexVector<Quat> quat(g_flexLib);
		NvFlexVector<Vec3> pos(g_flexLib);
		quat.map();
		pos.map();
		pos.resize(nInst);
		quat.resize(nInst);
		int totalNbBody = nInst;
		quat.unmap();
		pos.unmap();
		NvFlexGetRigids(g_solver, NULL, NULL, NULL, NULL, NULL, NULL, NULL, quat.buffer, pos.buffer);
		quat.map();
		pos.map();
		for (size_t i = 0; i < num_shared_buffers_proteins; i++)
		{
			int start = i * shared_buffer_size_proteins;
			if (start >= nInst)
				break;
			for (size_t j = 0; j < shared_buffer_size_proteins; j++)
			{
				if ((j + start) >= nInst)
					return;
				int indice = j + start;
				float ind = (float)mInstances[indice].mMeshIndex;
				NvFlexExtAsset* asset = iBatches[mInstances[indice].mMeshIndex].mAsset;
				Vec3 center = Vec3(asset->shapeCenters[0], asset->shapeCenters[1], asset->shapeCenters[2]);
				center = Rotate(quat[indice], center);
				float p[4] = { -(pos[indice].x - center.x)*(1.0f / main_scale), (pos[indice].y - center.y)*(1.0f / main_scale), (pos[indice].z - center.z)*(1.0f / main_scale), ind };
				float q[4] = { quat[indice].x, -quat[indice].y, -quat[indice].z, quat[indice].w };
				for (size_t k = 0; k < 4; k++)
				{
					sharedMemBuffer_proteins[i][j * 4 + k] = p[j];
					sharedMemBuffer_rotation[i][j * 4 + k] = q[j];
					sharedMemBuffer_info[i][j * 4 + k] = 0;
					if (j == 3) sharedMemBuffer_proteins[i][j * 4 + k] = ind;
				}
				sharedMemBuffer_info[i][j * 4 + 0] = ind;
			}

		}
		UpdateSharedMemSpring();
	}
	*/

	void writeToBinary(std::string postfix = "", bool regular = true)
	{
		NvFlexVector<Quat> quat(g_flexLib);
		NvFlexVector<Vec3> pos(g_flexLib);
		quat.map();
		pos.map();

		int nInst = mInstances.size();
		int totalNbBody = 0;
		if (use_rb){
			pos.resize(mInstances.size());
			quat.resize(mInstances.size());
			totalNbBody = nInst;
		}
		else {

			for (int i = 0; i < mInstances.size(); i++){
				NvFlexExtAsset* asset =  iBatches[ mInstances[i].mMeshIndex].mAsset;
				totalNbBody += asset->numShapes;
			}
			pos.resize(totalNbBody);
			quat.resize(totalNbBody);
		}
		//USE g_buffers->rigidRotations.buffer AND g_buffers->rigidTranslations.buffer ?
		//NvFlexGetRigids(NvFlexSolver* solver, NvFlexBuffer* offsets, NvFlexBuffer* indices, NvFlexBuffer* restPositions, NvFlexBuffer* restNormals, NvFlexBuffer* stiffness, NvFlexBuffer* thresholds, NvFlexBuffer* creeps, NvFlexBuffer* rotations, NvFlexBuffer* translations);
		//NvFlexGetRigidTransforms(g_solver, quat.buffer, pos.buffer);
		//readback rigid transforms
		quat.unmap();
		pos.unmap();
		NvFlexGetRigids(g_solver, NULL, NULL, NULL, NULL, NULL, NULL, NULL, quat.buffer, pos.buffer);
		quat.map();
		pos.map();

		std:string fname = "../../data/cellpack/traj/pack_result" + postfix + ".bin";

		output_bin.open(fname.c_str(), ios::out  | ios::binary);//| ios::app

		//number of instances
		output_bin.write((char *)&nInst, sizeof(nInst));
		//number of controle points total
		int nptsTotal = 0;
		for (int r = 0; r < g_ropes.size(); r++){
			nptsTotal += g_ropes[r].mIndices.size();
		}
		output_bin.write((char *)&nptsTotal, sizeof(nptsTotal));
		//number of curve
		int nCurve = g_ropes.size();
		output_bin.write((char *)&nCurve, sizeof(nCurve));
		//number of rigidbody
		if (!regular) output_bin.write((char *)&totalNbBody, sizeof(totalNbBody));
		//number of membrane points total
		int npointmembrane = 0;
		if (mask_membrane.size() != 0)
		{
			for (int m = 0; m <  mask_membrane.size() / 2; m += 2)
			{
				npointmembrane +=  mask_membrane[m + 1];
			}
		}
		if (!regular) output_bin.write((char *)&npointmembrane, sizeof(npointmembrane));//nb of membrane point

		//write position
		if ( use_rb ){
			for (int i = 0; i < mInstances.size(); i++){
				//output << -pos[i].x*(1.0f/main_scale) << " " <<pos[i].y*(1.0f/main_scale)<< " " <<pos[i].z*(1.0f/main_scale)<< " " << cp->mInstances[i].mMeshIndex << "\n";
				float ind = (float) mInstances[i].mMeshIndex;
				NvFlexExtAsset* asset = iBatches[mInstances[i].mMeshIndex].mAsset;
				Vec3 center = Vec3 (asset->shapeCenters[0], asset->shapeCenters[1], asset->shapeCenters[2]);
				center = Rotate(quat[i], center);
				//subtract or add the center
				float p[4] = { -(pos[i].x- center.x)*(1.0f / main_scale), (pos[i].y - center.y)*(1.0f / main_scale), (pos[i].z - center.z)*(1.0f / main_scale), ind };
				output_bin.write((char *)&p, sizeof(float) * 4);
			}
			for (int i = 0; i < mInstances.size(); i++) {
				float p[4] = { quat[i].x, -quat[i].y, -quat[i].z, quat[i].w };
				output_bin.write((char *)&p, sizeof(float) * 4);
			}
			//output_bin.write((char *)&quat[0], sizeof(float) * 4 * nInst);
		}
		else {
			int count = 0;
			for (int i = 0; i <  mInstances.size(); i++){
				float ind = (float) mInstances[i].mMeshIndex;
				NvFlexExtAsset* asset =  iBatches[ mInstances[i].mMeshIndex].mAsset;
				for (int j = 0; j < asset->numShapes; ++j)
				{
					float p[4] = { -pos[count].x*(1.0f / main_scale), pos[count].y*(1.0f / main_scale), pos[count].z*(1.0f / main_scale), ind };
					output_bin.write((char *)&p, sizeof(float) * 4);
					count++;
				}
			}
			output_bin.write((char *)&quat[0], sizeof(float) * 4 * totalNbBody);
		}
		//g_buffers->positions.map();
		//g_buffers->normals.map();
		//write control point
		int nRope = g_ropes.size();
		std::vector<Vec4> positions;
		std::vector<Vec4> normals;
		std::vector<Vec4> infos;
		//should the guffer be map ?
		std::vector<Vec4> curve_infos;
		for (int r = 0; r<g_ropes.size(); r++){
			int npts = g_ropes[r].mIndices.size();
			float cType = (float) maks_fiber[r];
			float cId = (float)r;
			curve_infos.push_back(Vec4(cType,0,0,0)); ////ing_id,rs,rc,tu
			//cout << "write rope " << r << " cType " << cType << " cId " << cId << " npts " << npts << endl;
			for (int i = 0; i<npts; i++){
				Vec4 pos = g_buffers->positions[g_ropes[r].mIndices[i]];
				//unsigned int  typeind = (cType << 24) | (unsigned int) r;///combine cType and r
				//float p[4] = { -pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale),0.0f };
				positions.push_back(Vec4(-pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale), main_radius / main_scale));
				//output_bin.write((char *)&p, sizeof(float) * 4);
				Vec4 norm = g_buffers->normals[g_ropes[r].mIndices[i]];
				//float n[4] = { -norm.x, norm.y, norm.z,1.0f };
				normals.push_back(Vec4(-norm.x, norm.y, norm.z, 1.0f));
				//output_bin.write((char *)&n, sizeof(float) * 4);
				//info
				infos.push_back(Vec4(cId, cType, 0.0f, 0.0f));//curveId, curveType, maxAngle, ulength
				//float info[4] = { cId, cType, 0.0f,0.0f };//curveId, curveType, maxAngle, ulength
				//output_bin.write((char *)&info, sizeof(float) * 4);
				//output_bin.write((char *)&cType, sizeof(float));
				//output_bin.write((char *)&cId, sizeof(float));
			}
		}
		output_bin.write((char *)&positions[0], sizeof(float) * 4 * nptsTotal);
		output_bin.write((char *)&normals[0], sizeof(float) * 4 * nptsTotal);
		output_bin.write((char *)&infos[0], sizeof(float) * 4 * nptsTotal);
		output_bin.write((char *)&curve_infos[0], sizeof(float) * 4 * nRope);
		//curve infos
		//write membrane point
		if (!regular) {
			int typemb = 0;
			for (int m = 0; m < mask_membrane.size() / 2; m += 2)
			{
				int start = mask_membrane[m];
				int npoints = mask_membrane[m + 1];
				//cout << "mb npoints " << start << " " << npoints << endl;
				for (int i = start; i < start + npoints; i++)
				{
					Vec4 pos = g_buffers->positions[i];
					float p[3] = { -pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale) };
					output_bin.write((char *)&p, sizeof(float) * 3);
					unsigned int  typem = typemb;
					output_bin.write((char *)&typem, sizeof(unsigned int));
				}
				typemb += 1;
			}
		}
		output_bin.close();
	}

	void jitter(bool allparticles, bool diffuse){
		//pick one particle per instance, apply some velocity change
		//stay in vicinity of sphere
		float weight = 1;
		for (int i = 0; i < int( mInstances.size()); ++i)
		{
			if (!allparticles){
				int poffseti =  mInstances[i].mParticleOffset;
				Vec3 toward_center = -Vec3(g_buffers->positions[poffseti]);
				float d = Length(toward_center);
				if (d < 10) weight = 0;
				else weight = dojitter_biased_strength;//if (Length(toward_center) > 8.0f) 
				//Vec3 p = Lerp(Vec3(g_positions[poffseti]), Vec3(g_positions[poffseti]) + Normalize(toward_center)*weight + RandomUnitVector() * 20, 0.8f);
				Vec3 p = Lerp(Vec3(g_buffers->positions[poffseti]), Vec3(g_buffers->positions[poffseti]) + RandomUnitVector()*dojitter_strength, 0.8f);
				if (dojitter_biased){
					p = Lerp(Vec3(g_buffers->positions[poffseti]), Vec3(g_buffers->positions[poffseti]) + Normalize(toward_center)*weight + RandomUnitVector()*dojitter_strength, 0.8f);
				}
				Vec3 delta = Vec3(0,0,0);// p - Vec3(g_buffers->positions[poffseti]);
				float doff = Length(Vec3(g_buffers->positions[poffseti]));
				if (iBatches[mInstances[i].mMeshIndex].compId < 0) {
					//inside 
					NvFlexExtAsset* asset = iBatches[mInstances[i].mMeshIndex].mAsset;
					Vec3 ppos = Vec3(g_buffers->positions[mInstances[i].mParticleOffset]); // asset->numParticles
					Vec3 ppos2 = Vec3(g_buffers->positions[mInstances[i].mParticleOffset + asset->numParticles-1]);
					int cid = iBatches[mInstances[i].mMeshIndex].compId;
					// cout << "query distance with  " << cid << " " << poffseti << " instance " << i << endl;
					float D1 = GetSurfaceDistance(ppos, 0);
					float D2 = GetSurfaceDistance(ppos2, 0);
					if ((D1 >= 0 || D2 >= 0) && cid < 0) {
						//cout << "distance is " << D1 << " " << D2 << " " << i << endl;
						Vec3 normal = GetAwayFromSurface(ppos, 0);
						// cout << "normal is " << normal.x << " " << normal.y << " " << normal.z << endl;
						delta = -normal * dojitter_strength;
						//g_buffers->velocities[mInstances[i].mParticleOffset].x = delta.x;
						//g_buffers->velocities[mInstances[i].mParticleOffset].y = delta.y;
						//g_buffers->velocities[mInstances[i].mParticleOffset].z = delta.z;
						g_buffers->positions[mInstances[i].mParticleOffset].x += delta.x;
						g_buffers->positions[mInstances[i].mParticleOffset].y += delta.y;
						g_buffers->positions[mInstances[i].mParticleOffset].z += delta.z;
					}
					/*
					//delta = Normalize(toward_center)*dojitter_strength;
					if (doff >= (comp_radius - g_params.radius * 2)) {
						p = Normalize(Vec3(g_buffers->positions[poffseti])) * (comp_radius - g_params.radius * 2);
						delta = Normalize(toward_center)*dojitter_strength;
						g_buffers->positions[poffseti].x = p.x;
						g_buffers->positions[poffseti].y = p.y;
						g_buffers->positions[poffseti].z = p.z;
					}
					else {
					
					}
					*/
				}


				/*g_buffers->velocities[poffseti].x = delta.x / g_dt;
				g_buffers->velocities[poffseti].y = delta.y / g_dt;
				g_buffers->velocities[poffseti].z = delta.z / g_dt;
				*/
			}
			else {
				int poffseti =  mInstances[i].mParticleOffset;
				NvFlexExtAsset* asset =  iBatches[ mInstances[i].mMeshIndex].mAsset;
				Vec3 toward_center = -Vec3(g_buffers->positions[poffseti]);
				float d = Length(toward_center);
				if (d < 1.5f) weight = 0;
				else weight = dojitter_biased_strength / 100.0f;//if (Length(toward_center) > 8.0f) 
				Vec3 p = Lerp(Vec3(g_buffers->positions[poffseti]), Vec3(g_buffers->positions[poffseti]) + RandomUnitVector()*(dojitter_strength / 100.0f), 0.8f);
				if (dojitter_biased){
					p = Lerp(Vec3(g_buffers->positions[poffseti]), Vec3(g_buffers->positions[poffseti]) + Normalize(toward_center)*weight + RandomUnitVector()*(dojitter_strength / 100.0f), 0.8f);
				}
				Vec3 delta = p - Vec3(g_buffers->positions[poffseti]);
				Quat q = QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);
				//Rotate(g_rigidRotations[rigidIndex], localPos)
				//rotate the point and apply the forice
				//can we access g_buffers->rigidTranslations[i]
				Vec3 roff = Vec3(0.0);
				/*if (iBatches[mInstances[i].mMeshIndex].compId > 0) {
					float dm = Length(g_buffers->rigidTranslations[i]);
					if (abs(dm - comp_radius) > 0.1) 
					{
						//x--d---o---d
						//std::cout <<" distance "<< dm << " " << g_buffers->rigidTranslations[i].x << " " << g_buffers->rigidTranslations[i].y << " " << g_buffers->rigidTranslations[i].z << endl;
						//roff = Normalize(g_buffers->rigidTranslations[i]) * (comp_radius - dm);
						//in that case apply an offset on all particles ?
					}
				}
				*/
				if (diffuse && iBatches[mInstances[i].mMeshIndex].compId < 0) {
					Json::Value ingr_node1 = getProteinNode(mInstances[i].mMeshIndex);
					Vec3 direction_random = UniformSampleSphere();
					Vec3 rotation_random = UniformSampleSphere();
					int ri = (int) Random(0.0f, 0.9f) * asset->numParticles;
					float temperature = 298.0f;
					float Radius = maxf(ingr_node1["encapsulatingRadius"].asFloat(), 10.0f);//nm
					float r = Random(0.0f, 1.0f);
					float dtheta = 0.0f;
					if (r < 0.5f)
						dtheta = 1.0f;
					else
						dtheta = -1.0f;
					//? ? 25 0.000193045 0.000193045 0.000240098
					//? ? 0 0 0
					direction_random = direction_random * sqrt((6.0f * 0.245f * g_realdt) / Radius)*main_scale;//in nm
					//dtheta = dtheta * sqrt((2.0f * 0.184f * g_realdt) / (Radius * Radius * Radius));//in rad
					//float angle = (180.0f / M_PI) * dtheta * sqrt(temperature / 298.0f) * 10.0f;
					Vec3 new_pos = direction_random * sqrt(temperature / 298.0f)*200.0f;// *dojitter_strength;// 1.0f / uScale;g_buffers->rigidTranslations[i] +
					//Quat new_rot = QuatFromAxisAngle(rotation_random, angle); //g_buffers->rigidRotations[i] *
					//if (i < 1) {
					//	std::cout << " direction_random " << Radius << " " << direction_random.x << " " << direction_random.y << " " << direction_random.z << endl;
					//	std::cout << " new_pos " << Radius << " " << new_pos.x << " " << new_pos.y << " " << new_pos.z << endl;
					//}
					//Vec3 localPos = g_buffers->positions[poffseti];//- Vec3(&asset->shapeCenters[0]);local position of the proxy
					//Vec3 rp = new_pos + localPos;
					//if (i < 1) std::cout << " rp " << rp.x << " " << rp.y << " " << rp.z << endl;
					g_buffers->velocities[poffseti+ri].x = new_pos.x;
					g_buffers->velocities[poffseti + ri].y = new_pos.y;
					g_buffers->velocities[poffseti + ri].z = new_pos.z;
					g_buffers->positions[poffseti + ri].x += new_pos.x;
					g_buffers->positions[poffseti + ri].y += new_pos.y;
					g_buffers->positions[poffseti + ri].z += new_pos.z;
				}
				for (int j = 0; j < asset->numParticles; j++){
					Vec3 p_rot = Rotate(q, Vec3(g_buffers->positions[poffseti + j])) - Vec3(g_buffers->positions[poffseti + j]);
					Vec3 delta_rot = Vec3(0, 0, 0);// Lerp(Vec3(g_positions[poffseti + j]), p_rot, 0.8f)-Vec3(g_positions[poffseti + j]);
					float doff = Length(Vec3(g_buffers->positions[poffseti + j]));
					Vec3 new_pos;
					Quat new_rot;
					if (iBatches[mInstances[i].mMeshIndex].compId < 0) {
						//delta = Normalize(toward_center)*dojitter_strength;
						p = Normalize(Vec3(g_buffers->positions[poffseti + j])) * (comp_radius - g_params.radius * 2);
						if (doff >= (comp_radius - g_params.radius * 2)) {
							//Vec3 delta = p - Vec3(g_buffers->positions[poffseti]);
							g_buffers->velocities[poffseti + j].x = dojitter_strength * -p.x / g_dt;
							g_buffers->velocities[poffseti + j].y = dojitter_strength * -p.y / g_dt;
							g_buffers->velocities[poffseti + j].z = dojitter_strength * -p.z / g_dt;
						}
						if (doff >= comp_radius) {
							g_buffers->positions[poffseti + j].x = p.x;
							g_buffers->positions[poffseti + j].y = p.y;
							g_buffers->positions[poffseti + j].z = p.z;
						}
						//if (diffuse) {
						//	Vec3 localPos = g_buffers->positions[poffseti + j];//- Vec3(&asset->shapeCenters[0]);local position of the proxy
						//	Vec3 rp = new_pos + localPos;
						//	if (i < 1) std::cout << " rp " << rp.x << " " << rp.y << " " << rp.z << endl;
						//	g_buffers->velocities[poffseti + j].x = new_pos.x / g_dt;
						//	g_buffers->velocities[poffseti + j].y = new_pos.y / g_dt;
						//	g_buffers->velocities[poffseti + j].z = new_pos.z / g_dt;
						//}
					}
					else if (iBatches[mInstances[i].mMeshIndex].compId > 0) {
						//surface ingredient
						//go in localspace by removing te instance position ?
						//g_buffers->positions[poffseti + j].x += roff.x;
						//g_buffers->positions[poffseti + j].y += roff.y;
						//g_buffers->positions[poffseti + j].z += roff.z;
					}
					//g_buffers->velocities[poffseti + j].x = (delta.x + delta_rot.x) / g_dt;
					//g_buffers->velocities[poffseti + j].y = (delta.y + delta_rot.y) / g_dt;
					//g_buffers->velocities[poffseti + j].z = (delta.z + delta_rot.z) / g_dt;
				}

			}
		} // do the curve
		int serial_id = 0;
		for (int r = 0; r < g_ropes.size(); r++) {
			int npts = g_ropes[r].mIndices.size();
			float cType = (float)maks_fiber[r];
			float cId = (float)r;
			//cout << "write rope " << r << " cType " << cType << " cId " << cId << " npts " << npts << endl;
			for (int i = 0; i < npts; i++) {
				serial_id = g_ropes[r].mIndices[i];
				Vec4 pos = g_buffers->positions[serial_id];
				float doff = Length(Vec3(g_buffers->positions[serial_id]));
				Vec3 p = Normalize(Vec3(g_buffers->positions[serial_id])) * (comp_radius - g_params.radius * 2);
				/*if (doff >= (comp_radius - g_params.radius * 2)) {
					//Vec3 delta = p - Vec3(g_buffers->positions[serial_id]);
					g_buffers->velocities[serial_id].x = dojitter_strength * -p.x / g_dt;
					g_buffers->velocities[serial_id].y = dojitter_strength * -p.y / g_dt;
					g_buffers->velocities[serial_id].z = dojitter_strength * -p.z / g_dt;
				}
				if (doff >= comp_radius) {
					g_buffers->positions[serial_id].x = p.x;
					g_buffers->positions[serial_id].y = p.y;
					g_buffers->positions[serial_id].z = p.z;
				}*/
				// cout << "query distance with  " << cid << " " << poffseti << " instance " << i << endl;
				float D1 = GetSurfaceDistance(Vec3(g_buffers->positions[serial_id]), 0);
				if (D1 >= 0  && cId < 0) {
					//cout << "rope distance is " << D1 << " " << serial_id << endl;
					Vec3 normal = GetAwayFromSurface(Vec3(g_buffers->positions[serial_id]), 0);
					// cout << "normal is " << normal.x << " " << normal.y << " " << normal.z << endl;
					Vec3 delta = -normal * dojitter_strength;
					//g_buffers->velocities[serial_id].x = delta.x;
					//g_buffers->velocities[serial_id].y = delta.y;
					//g_buffers->velocities[serial_id].z = delta.z;
					g_buffers->positions[serial_id].x += delta.x;
					g_buffers->positions[serial_id].y += delta.y;
					g_buffers->positions[serial_id].z += delta.z;
				}
			}
		}
	}

	//increase damping to max and do diffusion based on size
	void RandomDiffusion() {
		// MapBuffers(g_buffers);
		for (int i = 0; i < int(mInstances.size()); ++i)
		{
			Instance& inst = mInstances[i];
			if (iBatches[mInstances[i].mMeshIndex].compId >= 0) 
			{
				continue;
			}
			if (i > 100) continue;
			//Quat q = QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);
			Json::Value ingr_node1 = getProteinNode(inst.mMeshIndex);
			Vec3 direction_random = UniformSampleSphere();
			Vec3 rotation_random = UniformSampleSphere();
			float temperature = 298.0f;
			float Radius = ingr_node1["encapsulatingRadius"].asFloat();//nm
			float r = Random(0.0f, 1.0f);
			float dtheta = 0.0f;
			if (r < 0.5f)
				dtheta = 1.0f;
			else
				dtheta = -1.0f;
			direction_random = direction_random * sqrt((6.0f * 0.245f * g_realdt) / Radius)*main_scale;//in nm
			dtheta = dtheta * sqrt((2.0f * 0.184f * g_realdt) / (Radius * Radius * Radius));//in rad
			float angle = (180.0f / M_PI) * dtheta * sqrt(temperature / 298.0f) * 10.0f;
			Vec3 new_pos =  direction_random * sqrt(temperature / 298.0f);// 1.0f / uScale;g_buffers->rigidTranslations[i] +
			Quat new_rot =  QuatFromAxisAngle(rotation_random, angle); //g_buffers->rigidRotations[i] *
			// cout << direction_random.x << endl;
			NvFlexExtAsset* asset = iBatches[inst.mMeshIndex].mAsset;
			for (int j = 0; j < asset->numParticles; ++j)
			{
				Vec3 localPos = Vec3(&asset->particles[j * 4]);//- Vec3(&asset->shapeCenters[0]);local position of the proxy
				Vec3 rpos = Rotate(new_rot, localPos);
				Vec3 rp = new_pos + rpos;
				// cout << rp.x << endl;
				// continue;
				// g_buffers->positions[inst.mParticleOffset + j] = Vec4(new_pos + rpos, 1.0f);
				//g_buffers->positions[inst.mParticleOffset + j].x = rp.x;
				//g_buffers->positions[inst.mParticleOffset + j].y = rp.y;
				//g_buffers->positions[inst.mParticleOffset + j].z = rp.z;
				g_buffers->velocities[inst.mParticleOffset + j].x += rp.x ;
				g_buffers->velocities[inst.mParticleOffset + j].y += rp.y ;
				g_buffers->velocities[inst.mParticleOffset + j].z += rp.z ;
			}
			//g_buffers->rigidTranslations[i].Set(new_pos.x, new_pos.y, new_pos.z);
			//g_buffers->rigidRotations[i].Set(new_rot.x, new_rot.y, new_rot.z, new_rot.w);
		}
		// UnmapBuffers(g_buffers);
	}

	void stayInside() {
		for (int i = 0; i < int(mInstances.size()); ++i)
		{
			int poffseti = mInstances[i].mParticleOffset;
			//distance to surface
			//Length(toward_center)
			Vec3 toward_center = -Vec3(g_buffers->positions[poffseti]);
			float d = Length(toward_center);
		}
	}

	virtual void insertPointRope(int rope_id, int ipair, float D, float give,
		float stiffness, int aphase, int current)
	{
		Rope& arope = g_ropes[rope_id];
		//std::cout << current << " insertPoint at : " << ipair << " particle id in rope " << arope.mIndices[ipair] << " " << D << endl;
		if (arope.mIndices[ipair] == 0) return;
		Vec4 point = Vec4(g_buffers->positions[arope.mIndices[ipair]]);
		Vec4 next_point = Vec4(g_buffers->positions[arope.mIndices[ipair + 1]]);
		Vec3 dirtopt = Vec3(next_point.x - point.x, next_point.y - point.y, next_point.z - point.z);

		g_buffers->positions[current] = Vec4(point.x + dirtopt.x*0.5f,
			point.y + dirtopt.y*0.5f,
			point.z + dirtopt.z*0.5f,
			1.0f);
		g_buffers->velocities[current] = 0.0f;
		g_buffers->phases[current] = aphase;

		int p_id = arope.mIndices[ipair];
		//insert new point after ipair
		std::vector<int>::iterator it = arope.coarseIndices.begin();
		arope.coarseIndices.insert(it + ipair + 1, current);
		it = arope.mIndices.begin();
		arope.mIndices.insert(it + ipair + 1, current);

		//std::cout << current << " insertedPoint at : " << ipair + 1 << " " << arope.mIndices[ipair + 1] << " particle id in rope " << arope.mIndices[ipair] << " " << ipair << endl;

		//restort the spring for the given persistence.
		//arope.persistence = 2;
		CreateInsertedPersistence(arope, ipair, stiffness, give, D);
		/*ReplaceSpring(p_id,arope.mIndices[ipair+1],p_id,current);
		CreateSpring(current, arope.mIndices[ipair+1] ,stiffness, give,D);//current +1 doesnt exist yet
		//1-3
		ReplaceSpring(arope.mIndices[ipair-1],arope.mIndices[ipair+1],arope.mIndices[ipair-1],current);//replace
		ReplaceSpring(p_id,arope.mIndices[ipair+2],p_id,arope.mIndices[ipair+1]);
		float r=Randf(-1.0f,1.0f)/2000.0f;
		CreateSpring(current, arope.mIndices[ipair+2], stiffness, give,(D*2.0f)+r);//*0.5f*/
	}
	
	int countOverlap() {
		//maxParticles*maxParticleNeighbors ints
		//maxParticles ints in length
		//buffer of indices ints
		NvFlexVector<int> neighborsBuffer(g_flexLib);
		NvFlexVector<int> countsBuffer(g_flexLib);
		NvFlexVector<int> apiToInternalBuffer(g_flexLib);
		NvFlexVector<int> internalToApiBuffer(g_flexLib);
		neighborsBuffer.map();
		countsBuffer.map();
		apiToInternalBuffer.map();
		internalToApiBuffer.map();
		neighborsBuffer.resize(maxParticles*g_maxNeighborsPerParticle);
		countsBuffer.resize(maxParticles);
		apiToInternalBuffer.resize(maxParticles);
		internalToApiBuffer.resize(maxParticles);
		neighborsBuffer.unmap();
		countsBuffer.unmap();
		apiToInternalBuffer.unmap();
		internalToApiBuffer.unmap();
		NvFlexGetNeighbors(g_solver, neighborsBuffer.buffer, countsBuffer.buffer, apiToInternalBuffer.buffer, internalToApiBuffer.buffer);
		//neighborsBuffer.map();
		//countsBuffer.map();
		//apiToInternalBuffer.map();
		//internalToApiBuffer.map();
		int* neighbors = (int*)NvFlexMap(neighborsBuffer.buffer, 0);
		int* counts = (int*)NvFlexMap(countsBuffer.buffer, 0);
		int* apiToInternal = (int*)NvFlexMap(apiToInternalBuffer.buffer, 0);
		int* internalToApi = (int*)NvFlexMap(internalToApiBuffer.buffer, 0);
		// neighbors are stored in a strided format so that the first neighbor
		// of each particle is stored sequentially, then the second, and so on
		int stride = maxParticles;
		int total_overlap = 0;
		float cutoff = main_radius * 2.0f;// g_params.radius;
		for (int i = 0; i < maxParticles; ++i)
		{
			// find offset in the neighbors buffer
			int offset = apiToInternal[i];
			int count = counts[offset];
			int phase = g_buffers->phases[offset];
			for (int c = 0; c < count; ++c)
			{
				//int neighbor = remap[neighbors[c*stride + offset]];
				int phase_neighbor = g_buffers->phases[offset];
				int neighbor = internalToApi[neighbors[c*stride + offset]];
				float distance = Length(Vec3(g_buffers->positions[offset]) - Vec3(g_buffers->positions[neighbor]));
				//check the phase ?
				//printf("All Particle %d's phase %d neighbor %d is particle %d phase %d at distance %f\n", i, phase, c, neighbor, phase_neighbor, distance);
				if (distance < cutoff) {
					printf("Pass cutoff Particle %d's phase %d neighbor %d is particle %d phase %d at distance %f\n", i, phase, c, neighbor, phase_neighbor, distance);
					total_overlap++;
				}
			}
		}
		NvFlexUnmap(neighborsBuffer.buffer);
		NvFlexUnmap(countsBuffer.buffer);
		NvFlexUnmap(apiToInternalBuffer.buffer);
		NvFlexUnmap(internalToApiBuffer.buffer);
		return total_overlap;
	}

	virtual void addPointRope(int rope_id, int ipair, float D, float give,
		float stiffness, int aphase, int current)
	{
		Rope& arope = g_ropes[rope_id];
		//std::cout << current << " insertPoint at : " << ipair << " particle id in rope " << arope.mIndices[ipair] << " " << D << endl;
		if (arope.mIndices[ipair] == 0) return;
		Vec4 point = Vec4(g_buffers->positions[arope.mIndices[ipair]]);
		Vec4 next_point = Vec4(g_buffers->positions[arope.mIndices[ipair + 1]]);
		Vec3 dirtopt = Vec3(next_point.x - point.x, next_point.y - point.y, next_point.z - point.z);

		g_buffers->positions[current] = Vec4(point.x + dirtopt.x*0.5f,
			point.y + dirtopt.y*0.5f,
			point.z + dirtopt.z*0.5f,
			1.0f);
		g_buffers->velocities[current] = 0.0f;
		g_buffers->phases[current] = aphase;

		int p_id = arope.mIndices[ipair];
		//insert new point after ipair
		std::vector<int>::iterator it = arope.coarseIndices.begin();
		arope.coarseIndices.insert(it + ipair + 1, current);
		it = arope.mIndices.begin();
		arope.mIndices.insert(it + ipair + 1, current);

		//std::cout << current << " insertedPoint at : " << ipair + 1 << " " << arope.mIndices[ipair + 1] << " particle id in rope " << arope.mIndices[ipair] << " " << ipair << endl;

		//restort the spring for the given persistence.
		//arope.persistence = 2;
		CreateInsertedPersistence(arope, ipair, stiffness, give, D);
		/*ReplaceSpring(p_id,arope.mIndices[ipair+1],p_id,current);
		CreateSpring(current, arope.mIndices[ipair+1] ,stiffness, give,D);//current +1 doesnt exist yet
		//1-3
		ReplaceSpring(arope.mIndices[ipair-1],arope.mIndices[ipair+1],arope.mIndices[ipair-1],current);//replace
		ReplaceSpring(p_id,arope.mIndices[ipair+2],p_id,arope.mIndices[ipair+1]);
		float r=Randf(-1.0f,1.0f)/2000.0f;
		CreateSpring(current, arope.mIndices[ipair+2], stiffness, give,(D*2.0f)+r);//*0.5f*/
	}

	int proceduralRope(int rope_id, int Nsub, bool insert = false){
		int rope_phase = g_ropes[rope_id].phase;
		float D = g_params.radius;// main_radius*2.0f;//g_params.mRadius;//*1000.0f;//could be the persitence length, which is 450bp 34A is one turn 10bp. Distance between spheres
		//D = D + D/10.0f;
		float give = 0.0f;
		float stiffness = 1.0f;
		//initialized if not
		//update frequency
		//update nb
		//every n frame add m particle to rope r
		int current = NvFlexGetActiveCount(g_solver);
		int ipair = int(Randf(0.75f, 0.90f)*g_ropes[0].mIndices.size());
		//std::cout << current << " insert at : " << ipair << " particle id in rope " << g_ropes[0].mIndices[ipair] << " " << D << " " << rope_phase << " " << g_ropes[0].mIndices.size() << endl;
		for (int i = 0; i < Nsub; i++){
			//insert or add
			if (insert) insertPointRope(rope_id, ipair, D, give, stiffness, rope_phase, current);
			else addPointRope(rope_id, ipair, D, give, stiffness, rope_phase, current);
			current++;
		}
		//insert Hue ?
		float h = Randf(0.0f, 1.0f);
		//if (h < 0.045390661f)
		//    makeHue(ipair);//need to record the position of hue
		//0.045390661
		NvFlexCopyDesc copyDesc;
		copyDesc.dstOffset = 0;
		copyDesc.srcOffset = 0;
		copyDesc.elementCount = g_buffers->positions.size();

		NvFlexSetSprings(g_solver, g_buffers->springIndices.buffer, g_buffers->springLengths.buffer, 
								   g_buffers->springStiffness.buffer, g_buffers->springLengths.size());
		NvFlexSetParticles(g_solver, g_buffers->positions.buffer, &copyDesc);
		NvFlexSetVelocities(g_solver, g_buffers->velocities.buffer, &copyDesc);
		NvFlexSetPhases(g_solver, g_buffers->phases.buffer, &copyDesc);
		g_buffers->activeIndices.resize(current);// g_positions.size());
		for (size_t i = 0; i < current; i++)// g_activeIndices.size(); ++i)
			g_buffers->activeIndices[i] = i;
		NvFlexSetActive(g_solver, g_buffers->activeIndices.buffer, &copyDesc);
		return ipair;
	}

	void proceduralGrowRope(int rope_id, int Nsub, int start, bool insert = false){
		int rope_phase = g_ropes[rope_id].phase;
		float D = g_params.radius;// main_radius*2.0f;//g_params.mRadius;//*1000.0f;//could be the persitence length, which is 450bp 34A is one turn 10bp. Distance between spheres
		//D = D + D/10.0f;
		float give = 0.0f;
		float stiffness = 1.0f;
		int persistence = 1;//part of the ingreident
		Rope arope = g_ropes[rope_id];
		int current = NvFlexGetActiveCount(g_solver);
		arope.mIndices.push_back(current);
		arope.coarseIndices.push_back(current);
		Vec4 point_start = Vec4(g_buffers->positions[start]);
		g_buffers->positions[current] = point_start;
		g_buffers->velocities[current] = 0.0f;
		g_buffers->phases[current] = rope_phase;
		CreateSpringInter(start, current, 1.0f, 0.0f, main_radius);
		current++;
		for (int i = 0; i < Nsub; i++){		
			Vec4 next_point = point_start + Vec4( UniformSampleSphere()*D,0);
			g_buffers->positions[current + i] = next_point;
			g_buffers->velocities[current+i] = 0.0f;
			g_buffers->phases[current+i] = rope_phase;
			arope.mIndices.push_back(current+i);
			arope.coarseIndices.push_back(current+i);
			CreatePersistence(arope, current + i, persistence, stiffness, 6, give, D);
		}
	}

	void proceduralGrowRopePush(int rope_id, int Nsub, int start, bool insert = false, int compid = 0){
		int rope_phase = g_ropes[rope_id].phase;
		float D = g_params.radius;// main_radius*2.0f;//g_params.mRadius;//*1000.0f;//could be the persitence length, which is 450bp 34A is one turn 10bp. Distance between spheres
		//D = D + D/10.0f;
		float give = 0.0f;
		float stiffness = 1.0f;
		int persistence = 2;//part of the ingredent
		Rope& arope = g_ropes[rope_id];
		int current = g_buffers->positions.size();
		arope.mIndices.push_back(current);
		arope.coarseIndices.push_back(current);
		Vec4 point_start = Vec4(g_buffers->positions[start]);
		g_buffers->positions.push_back(point_start);
		g_buffers->velocities.push_back(0.0f);
		g_buffers->phases.push_back(rope_phase);
		CreateSpringInter(start, current, 1.0f, 0.0f, D);
		Vec4 direction = Vec4(UniformSampleSphere(), 0);
		current++;
		Vec4 previous = point_start;
		for (int i = 0; i < Nsub; i++){
			Vec4 next_point = previous + direction*D;
			//check if insideif compId < 0 
			//cout << " comp " <<  compid << endl;
			if (compid<0) 
			{
				Vec4 temp_ = previous + direction*D*4.0f;
				//cout << "inside ? " << !isInside(Vec3(next_point.x, next_point.y, next_point.z), 0) << endl;
				if (!isInside(Vec3(temp_.x, temp_.y, temp_.z), 0))  {
					//cout << " not inside " << endl;
					for (int j = 0; j < 10; j++) {
						//find one inside
						direction = Vec4(UniformSampleSphere(), 0);
						next_point = previous + direction*D;
						temp_ = previous + direction*D*4.0f;
						//cout << temp_.x << " " << temp_.y << " " << temp_.z << endl;
						if (isInside(Vec3(temp_.x, temp_.y, temp_.z), 0)){
							break;
						}
					}
				}
			}
			g_buffers->positions.push_back(next_point);
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(rope_phase);
			arope.mIndices.push_back(current + i);
			arope.coarseIndices.push_back(current + i);
			CreatePersistence(arope, current + i, persistence, stiffness, 6, give, D);
			previous = next_point;
		}
	}

	float getFiberLength(int rope_id)
	{
		float L1 = 0.0f;
		//float L2=0.0f;
		for (int i = 0; i<g_ropes[rope_id].coarseIndices.size() - 2; i++){
			float length = Length(g_buffers->positions[g_ropes[rope_id].coarseIndices[i + 1]] - g_buffers->positions[g_ropes[rope_id].coarseIndices[i]]);
			L1 += length;
			//L2+=length/4.0f;
		}
		printf("length is currently %i %f %f\n", g_ropes[rope_id].coarseIndices.size(), L1*1.0f / main_scale, (L1 / 4.0f)*1.0f / main_scale);
		return L1 / main_scale;
	}

	void growFiber(int Nsub)
	{
		for (int i = 0; i<maks_fiber.size(); i++)
		{
			int fiber_id = maks_fiber[i];
			int compId = iBatchesFiber[fiber_id].compId;
			if (compId > 0) {
				//placeOneFiberSurface(i);
			}
			//how many of this is there already
			//how long they are
			//binder ?
			//compid should gave us the mesh
			else {
				Json::Value ingr_node_name = getIngredients(pnames_fiber[fiber_id]);
				cout << "place fibers " << pnames_fiber[fiber_id] << " rope id " << fiber_id << " " << i << " compId " << compId << endl;
				float target_length = ingr_node_name["length"].asFloat();
				float current_length = getFiberLength(i);
				if (current_length < target_length)
				{
					int curve_index = proceduralRope(i, Nsub);
					//decide if attach a partner here
					cout << "partners ? " << iPartnersFibers[fiber_id].nPartner << endl;
					//if (iPartnersFibers[fiber_id].nPartner > 0)
					//	{
					//	pickPlacePartner(fiber_id, i, curve_index, ingr_node_name);
					//}
				}
				break;
			}
		}
	}
	
	void fiberToGrow(){
		if (fiber_togrow.size() == 0) return;
		for (int i = 0; i < fiber_togrow.size(); i++)
		{
			
			int fiber_id = fiber_togrow[i]; //maks_fiber[rope_id];
			int start_bead = fiber_togrow_stpos[i];
			int compId = iBatchesFiber[fiber_id].compId;
			int rope_id = g_ropes.size();
			Rope curve;
			curve.phase = NvFlexMakePhase(iGroupCounter++, eNvFlexPhaseSelfCollide);
			g_ropes.push_back(curve);//instance
			mask.push_back(-1);
			maks_protein.push_back(-(fiber_id + 1));//it ptype is 0 ?
			maks_fiber.push_back(fiber_id);//proteinType

			if (compId > 0) {
				//placeOneFiberSurface(i);
			}
			else 
			{
				Json::Value ingr_node_name = getIngredients(pnames_fiber[fiber_id]);
				cout << "place fibers " << pnames_fiber[fiber_id] << " rope id " << rope_id << " compId " << compId << " starting at " << start_bead << endl;
				float target_length = ingr_node_name["length"].asFloat();
				//float current_length = getFiberLength(rope_id);
				int rope_phase = g_ropes[rope_id].phase;
				float D = g_params.radius;// main_radius*2.0f;//g_params.mRadius;//*1000.0f;//could be the persitence length, which is 450bp 34A is one turn 10bp. Distance between spheres
				//D = D + D/10.0f;
				int Nsub = (int) ((target_length*main_scale) / D);//280-300 for rna, but also depends on position on dna gene
				if (fiber_togrow_length[i] != -1) {
					Nsub = fiber_togrow_length[i];
				}
				//cout << " will start with length " << Nsub << " " << fiber_togrow_length[i] <<  endl;
				float give = 0.0f;
				float stiffness = 1.0f;
				proceduralGrowRopePush(rope_id, Nsub, start_bead, false, compId);
				//cout << "npoint is " << g_ropes[rope_id].mIndices.size() << endl;
				//checkPartnerAlongCurvePoints(g_ropes[rope_id], ingr_node_name);
			}
		}
		fiber_togrow.resize(0);
		fiber_togrow_stpos.resize(0);
		fiber_togrow_length.resize(0);
		//place ribosome
	}

	void exportPDBSoftTransform(int exp, int run)
	{
		ofstream of;
		string name = "../../data/pack_result_soft";
		bool append = false;
		if (append){
			name = name + ".txt";
			of.open(name.c_str(), ios::out | ios::app);
			of << "# " << mInstances.size() << " " << exp << " " << run << " " << endl;
		}
		else {
			name = name + "_" + std::to_string(exp) + "_" + std::to_string(run) + ".pdb";
			of.open(name.c_str(), ios::out);
		}

		//of.open("C:\\Users\\ludov\\OneDrive\\Documents\\cellVIEW\ -\ i\\Assets\\Data\\pack_result_soft.txt", ios::out);

		int nInst = mInstances.size();
		for (int i = 0; i <  mInstances.size(); i++){
			//grab the particles positions for each isntances
			int poffseti =  mInstances[i].mParticleOffset;
			NvFlexExtAsset* asset = iBatches[mInstances[i].mMeshIndex].mAsset;
			int nPart = asset->numParticles;
			for (int j = 0; j < nPart; j++){
				of << i << " " << j << " " << g_buffers->positions[poffseti + j].x*(1.0f / main_scale) << " " << g_buffers->positions[poffseti + j].y*(1.0f / main_scale) << " " << g_buffers->positions[poffseti + j].z*(1.0f / main_scale) << endl;
			}
		}
		of.close();
	}

	void writeFiberBinary(string fname){
		//std:string fname = "C:/Users/ludov/OneDrive/Documents/OnlinePacking_Tobias/cellVIEW-OP/Data/pack_result_fiber.bin";// "../../data/pack_result_fiber.bin";
		output_bin.open(fname.c_str(), ios::out | ios::app | ios::binary);
		//write control point
		int nRope = g_ropes.size();
		//number of controle points total
		int nptsTotal = 0;
		for (int r = 0; r < g_ropes.size(); r++){
			nptsTotal += g_ropes[r].mIndices.size();
		}
		output_bin.write((char *)&nptsTotal, sizeof(nptsTotal));
		for (int r = 0; r<g_ropes.size(); r++){
			int npts = g_ropes[r].mIndices.size();
			float cType = (float) maks_fiber[r];
			float cId = (float)r;
			for (int i = 0; i<npts; i++){
				Vec4 pos = g_buffers->positions[g_ropes[r].mIndices[i]];
				//unsigned int  typeind = (cType << 24) | (unsigned int) r;///combine cType and r
				float p[3] = { pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale) };
				output_bin.write((char *)&p, sizeof(float) * 3);
				output_bin.write((char *)&cType, sizeof(float));//fiber_type
				output_bin.write((char *)&cId, sizeof(float));//fiber_id
			}
		}
		output_bin.close();
	}

	void writeFiberTxt(string fname){
		//std:string fname = "C:/Users/ludov/OneDrive/Documents/OnlinePacking_Tobias/cellVIEW-OP/Data/pack_result_fiber.bin";// "../../data/pack_result_fiber.bin";
		ofstream output;
		output.open(fname.c_str(), ios::out);// | ios::app);
		//write control point
		int nRope = g_ropes.size();
		//number of controle points total
		int nptsTotal = 0;
		for (int r = 0; r < g_ropes.size(); r++){
			nptsTotal += g_ropes[r].mIndices.size();
		}
		output << nRope << " " << nptsTotal << " " << main_radius*(1.0f / main_scale) << endl;
		for (int r = 0; r<g_ropes.size(); r++){
			int npts = g_ropes[r].mIndices.size();
			float cType = (float)maks_fiber[r];
			float cId = (float)r;
			for (int i = 0; i<npts; i++){
				Vec4 pos = g_buffers->positions[g_ropes[r].mIndices[i]];
				//unsigned int  typeind = (cType << 24) | (unsigned int) r;///combine cType and r
				float p[3] = { -pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale) };
				output << p[0] << " " << p[1] << " " << p[2] << " " << cType << " " << cId << endl;
			}
		}
		output_bin.close();
	}

	void writeModelPDB(string postfix, float oscale) {
		//'ATOM  %5d %-4s %3s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s';
		//sprintf(formatString, serial, atomname, a.resname, defaults(a.chainname, ' '), a.resno, a.x, a.y, a.z, defaults(a.occupancy, 1.0), defaults(a.bfactor, 0.0), '', // segid
		//defaults(a.element, ''))
		int nTypes = iBatches.size();//number of ingredient
		int fTypes = pnames_fiber.size();
		std::vector<std::string> combs = getThreeLetterString(nTypes + fTypes);
		std::string fname = "../../data/cellpack/ascii/pack_result" + postfix + ".pdb";
		//output.open(fname.c_str(), ios::out);// | ios::app);
		FILE *fp;
		fp = fopen(fname.c_str(), "w");
		int serial_id = 0;
		int resno = 0;
		cout << nTypes << " instance " << fTypes << " fibers " << combs.size() << "letters " << endl;
		for (int i = 0; i < int(mInstances.size()); ++i)
		{
			resno = i;
			int ingType = mInstances[i].mMeshIndex;
			int poffseti = mInstances[i].mParticleOffset;
			NvFlexExtAsset* asset = iBatches[ingType].mAsset;
			//cout << i << " id " << ingType << " offset " << poffseti << " num beads "<< asset->numParticles << endl;
			for (int j = 0; j < asset->numParticles; j++) {
				serial_id = poffseti + j;
				Vec4 pos = g_buffers->positions[poffseti + j]*(1.0f / main_scale);
				if (serial_id >= 99999) serial_id = 99999;
				if (resno >= 99999) resno = 99999;
				fprintf(fp, "HETATM%5d %-4s %3s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n", 
					serial_id, "CA", combs[ingType], ' ', resno, pos.x*oscale, pos.y*oscale, pos.z*oscale,0.0,0.0,"","");
			}
		}
		cout <<"FIBERS"<< endl;
		for (int r = 0; r < g_ropes.size(); r++) {
			resno = 0;
			int npts = g_ropes[r].mIndices.size();
			float cType = (float)maks_fiber[r];
			float cId = (float)r;
			cout << "write rope " << r << " cType " << cType << " cId " << cId << " npts " << npts << endl;
			for (int i = 0; i < npts; i++) {
				serial_id = g_ropes[r].mIndices[i];
				if (serial_id >= 99999) serial_id = 99999;
				if (resno >= 99999) resno = 99999;
				Vec4 pos = g_buffers->positions[g_ropes[r].mIndices[i]] * (1.0f / main_scale);
				//positions.push_back(Vec4(-pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale), main_radius / main_scale));
				//infos.push_back(Vec4(cId, cType, 0.0f, 0.0f));//curveId, curveType, maxAngle, ulength
				//fprintf(fp, "%8d %4d %4d %10.3f %10.3f %10.3f\n", g_ropes[r].mIndices[i], (int)cType, (int)cId, pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale));
				fprintf(fp, "ATOM  %5d %-4s %3s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n",
					serial_id, "CA", combs[nTypes+ cType], ' ', resno, pos.x*oscale, pos.y*oscale, pos.z*oscale, 0.0, 0.0, to_string((int)cId), "");
				resno++;
			}
		}
		fclose(fp);
	}

	void writeModelTxt(string postfix) {
		//use the PDB format ?
		//'ATOM  %5d %-4s %3s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s';
		//currently "%8d %4d %4d %10.3f %10.3f %10.3f\n" Particle ID, Ingredient Type, Instance ID, X , Y , Z
		//write david's format.
		//string s = str( format("%2% %2% %1%\n") % "world" % "hello" );
		//cout << format("%1% %2% %|40t|%3%\n") % first[i] % last[i] % tel[i];
		//ofstream output;
		std::string fname = "../../data/cellpack/ascii/pack_result" + postfix + ".txt";
		//output.open(fname.c_str(), ios::out);// | ios::app);
		FILE *fp;
		fp = fopen(fname.c_str(), "w");
		for (int i = 0; i < int(mInstances.size()); ++i)
		{
			int poffseti = mInstances[i].mParticleOffset;
			NvFlexExtAsset* asset = iBatches[mInstances[i].mMeshIndex].mAsset;
			for (int j = 0; j < asset->numParticles; j++) {
				Vec3 pos = Vec3(g_buffers->positions[poffseti + j]);
				fprintf(fp, "%8d %4d %4d %10.3f %10.3f %10.3f\n", poffseti + j, mInstances[i].mMeshIndex, i, pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale));
			}
		}
		for (int r = 0; r < g_ropes.size(); r++) {
			int npts = g_ropes[r].mIndices.size();
			float cType = (float)maks_fiber[r];
			float cId = (float)r;
			cout << "write rope " << r << " cType " << cType << " cId " << cId << " npts " << npts << endl;
			for (int i = 0; i < npts; i++) {
				Vec4 pos = g_buffers->positions[g_ropes[r].mIndices[i]];
				//positions.push_back(Vec4(-pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale), main_radius / main_scale));
				//infos.push_back(Vec4(cId, cType, 0.0f, 0.0f));//curveId, curveType, maxAngle, ulength
				fprintf(fp, "%8d %4d %4d %10.3f %10.3f %10.3f\n", g_ropes[r].mIndices[i], (int) cType, (int) cId, pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale));
			}
		}
		fclose(fp);
	}

	float mAttractForce;	

	std::vector<MeshBatch> mBatches;
	int numActive;
	int mGroupCounter;

};
