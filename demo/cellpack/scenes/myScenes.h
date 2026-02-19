//
// Created by ludo on 2/10/16.
//

#ifndef FLEXPACK_MYSCENES_H
#define FLEXPACK_MYSCENES_H

#endif //FLEXPACK_MYSCENES_H

#include "halton.hpp"
#include "json/json.h"
#include "mmtf_parser.h"

//#include "include/flex.h"
//#include "include/flexExt.h"
//#include "../include/flexExt.h"

#include <iostream>
#include <fstream>

#include <stdlib.h> 


//GRAPH stuff
/*
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
#include <stdio.h>
#include <sys/types.h> 
/*
#if _WIN32
#undef UNICODE
//#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment(lib, "Ws2_32.lib")
#else
#include <sys/socket.h>
#include <netinet/in.h>
#endif
//#include <Python.h>
*/

#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#define ChangeDir _chdir
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#define ChangeDir chdir
#endif

char cCurrentPath[FILENAME_MAX];

#include <algorithm>

#define DEFAULT_BUFLEN 512
#define DEFAULT_PORT "12000"

#define BACKLOG 10

#define TAKE_N_BITS_FROM(b, p, n) ((b) >> (p)) & ((1 << (n)) - 1)


// basic SAP based acceleration structure for point cloud queries
//taken from flexExtSoft.cpp
struct SweepAndPrune
{
	struct Entry
	{
		Entry(Vec3 p, int i) : point(p), index(i) {}

		Vec3 point;
		int index;
	};

	SweepAndPrune(const Vec3* points, int n)
	{
		entries.reserve(n);
		for (int i = 0; i < n; ++i)
			entries.push_back(Entry(points[i], i));

		struct SortOnAxis
		{
			int axis;

			SortOnAxis(int axis) : axis(axis) {}

			bool operator()(const Entry& lhs, const Entry& rhs) const
			{
				return lhs.point[axis] < rhs.point[axis];
			}
		};

		// calculate particle bounds and longest axis
		Vec3 lower(FLT_MAX), upper(-FLT_MAX);
		for (int i = 0; i < n; ++i)
		{
			lower = Min(points[i], lower);
			upper = Max(points[i], upper);
		}

		Vec3 edges = upper - lower;

		if (edges.x > edges.y && edges.x > edges.z)
			longestAxis = 0;
		else if (edges.y > edges.z)
			longestAxis = 1;
		else
			longestAxis = 2;

		std::sort(entries.begin(), entries.end(), SortOnAxis(longestAxis));
	}

	void QuerySphere(Vec3 center, float radius, std::vector<int>& indices)
	{
		// find start point in the array
		int low = 0;
		int high = int(entries.size());

		// the point we are trying to find
		float queryLower = center[longestAxis] - radius;
		float queryUpper = center[longestAxis] + radius;

		// binary search to find the start point in the sorted entries array
		while (low < high)
		{
			const int mid = (high + low) / 2;

			if (queryLower > entries[mid].point[longestAxis])
				low = mid + 1;
			else
				high = mid;
		}

		// scan forward over potential overlaps
		float radiusSq = radius*radius;

		for (int i = low; i < int(entries.size()); ++i)
		{
			Vec3 p = entries[i].point;

			if (LengthSq(p - center) < radiusSq)
			{
				indices.push_back(entries[i].index);
			}
			else if (entries[i].point[longestAxis] > queryUpper)
			{
				// early out if ther are no more possible candidates
				break;
			}
		}
	}

	int longestAxis;	// [0,2] -> x,y,z

	std::vector<Entry> entries;
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


// creates distance constraints between particles within some radius
int CreateLinksLocal(const float* all_particles, int numParticles, std::vector<int> indices, std::vector<int>& outSpringIndices,
	std::vector<float>& outSpringLengths, std::vector<float>& outSpringStiffness, float radius, float stiffness = 1.0f, int offset = 0)
{
	int count = 0;

	std::vector<Vec3> particles;
	for (int i = 0; i < numParticles; ++i)
	{
		Vec3 localPos = Vec3(&all_particles[(indices[i]+offset) * 4]);
		particles.push_back(localPos);
	}

	std::vector<int> neighbors;
	SweepAndPrune sap(&particles[0], numParticles);

	for (int i = 0; i < numParticles; ++i)
	{
		neighbors.resize(0);

		sap.QuerySphere(Vec3(particles[i]), radius, neighbors);

		for (int j = 0; j < int(neighbors.size()); ++j)
		{
			const int nj = neighbors[j];

			if (nj != i)
			{
				outSpringIndices.push_back(indices[i] + offset);
				outSpringIndices.push_back(indices[nj] + offset);
				outSpringLengths.push_back(Length(Vec3(particles[i]) - Vec3(particles[nj])));
				outSpringStiffness.push_back(stiffness);

				++count;
			}
		}
	}
	return count;
}


struct SDF
{
	float mLower[3];			//!< Shape AABB lower bounds in world space
	float mUpper[3];			//!< Shape AABB upper bounds in world space
	float mInvEdgeLength[3];	//!< 1/(mUpper-mLower)

	unsigned int mWidth;		//!< Field x dimension in voxels
	unsigned int mHeight;		//!< Field y dimension in voxels
	unsigned int mDepth;		//!< Field z dimension in voxels

	float* mField;		//!< SDF voxel data, must be mWidth*mHeight*mDepth in length
	NvFlexDistanceFieldId fsdf;		//pointer to flexSDF
};



Mesh* GetMesh(const char* meshFile, float scale)
{
    Mesh* mesh = ImportMesh(meshFile);
    if (!mesh)
    {
        printf("Could not open mesh for reading: %s\n", meshFile);
        return NULL;
    }
    else
    {
        mesh->Transform(ScaleMatrix(scale));
        mesh->m_colours.resize(0);
    }
    return mesh;
}

/*
SDF* CreateSDFfromMesh(const char* sdfFilein, Mesh* mesh, float scale, Vec3 lower,
                           Vec3 minExtents,Vec3 maxExtents, float expand=0.0f)
{
	const int dim = 256;
    // voxelize mesh
    if (!mesh)
    {
        printf("Could not process mesh \n", sdfFilein);
		return NULL;
    }
    Vec3 edges;
    //mesh->GetBounds(minExtents, maxExtents);
    mesh->m_colours.resize(0);

    // square extents
    edges = maxExtents - minExtents;
    float longestAxis = max(max(edges.x, edges.y), edges.z);
    edges = longestAxis;

    minExtents = minExtents - edges * 0.1f;
    maxExtents = minExtents + edges * 1.1f;
    edges = maxExtents - minExtents;

    string sdfFile = string(sdfFilein, strrchr(sdfFilein, '.')) + ".pfm";

	PfmImage pfm;
	if (!PfmLoad(sdfFile.c_str(), pfm))
    {
        

		pfm.m_width = dim;
		pfm.m_height = dim;
		pfm.m_depth = dim;
		pfm.m_data = new float[dim*dim*dim];

        printf("Cooking SDF: %s - dim: %d^3\n", sdfFile.c_str(), dim);

		CreateSDF(mesh, dim, minExtents, maxExtents, pfm.m_data);

		PfmSave(sdfFile.c_str(), pfm);
    }

	printf("Loaded SDF, %d\n", pfm.m_width);

	assert(pfm.m_width == pfm.m_height && pfm.m_width == pfm.m_depth);

    // cheap collision offset
	int numVoxels = int(pfm.m_width*pfm.m_height*pfm.m_depth);
    for (int i=0; i < numVoxels; ++i)
    {
		pfm.m_data[i] += expand;
		//pfm.m_data[i] *= -1.0f;
    }

	// set up flex collision shape
	FlexSDF* sdf = flexCreateSDF();
	flexUpdateSDF(sdf, dim, dim, dim, pfm.m_data, eFlexMemoryHost);

	// entry in the collision->render map
	//g_fields[sdf] = CreateGpuMesh(mesh);

	//delete mesh;
	SDF* asdf = new SDF();

    // set up flex collision shape
	asdf->mWidth = pfm.m_width;
	asdf->mHeight = pfm.m_height;
	asdf->mDepth = pfm.m_depth;
	(Vec3&)asdf->mLower = minExtents;
	(Vec3&)asdf->mUpper = maxExtents;
	(Vec3&)asdf->mInvEdgeLength = Vec3(1.0f / edges.x, 1.0f / edges.y, 1.0f / edges.z);
	asdf->mField = pfm.m_data;
	asdf->fsdf = sdf;
	return asdf;
}

FlexExtAsset* flexExtCreateRigidFromPoints(std::vector<Vec3> points)
{
    Vec3 center;
    std::vector<Vec4> particles;

    FlexExtAsset* asset = new FlexExtAsset();
    memset(asset, 0, sizeof(*asset));

    if (points.size())
    {
        for (int i=0; i < points.size(); i++)
        {
            center += points[i];
            particles.push_back(Vec4(points[i].x, points[i].y, points[i].z, 1.0f));
        }
        const int numParticles = int(particles.size());
        asset->mNumParticles = numParticles;

        asset->mParticles = new float[numParticles*4];
        memcpy(asset->mParticles, &particles[0], sizeof(Vec4)*numParticles);

        // store center of mass
        center /= float(numParticles);
        center*=0.0;
        asset->mNumShapes = 1;
        asset->mNumShapeIndices = numParticles;

        // for rigids we just reference all particles in the shape
        asset->mShapeIndices = new int[numParticles];

        for (int i = 0; i < numParticles; ++i)
            asset->mShapeIndices[i] = i;

        asset->mShapeCenters = new float[4];
        asset->mShapeCenters[0] = center.x;
        asset->mShapeCenters[1] = center.y;
        asset->mShapeCenters[2] = center.z;

        asset->mShapeCoefficients = new float[1];
        asset->mShapeCoefficients[0] = 1.0f;

        asset->mShapeOffsets = new int[1];
        asset->mShapeOffsets[0] = numParticles;
    }
    return asset;
}
*/

void CreateSpringInter(int i, int j, float stiffness, float give=0.0f, float length=0.0f)
{
	//printf(" %i and %i \n", i, j);
    g_buffers->springIndices.push_back(i);
    g_buffers->springIndices.push_back(j);
    if (length!= 0.0f) g_buffers->springLengths.push_back(length);
    else g_buffers->springLengths.push_back((1.0f+give)*Length(Vec3(g_buffers->positions[i])-Vec3(g_buffers->positions[j])));
    g_buffers->springStiffness.push_back(stiffness);
}


void WeightSpringStifness(std::vector<int> indices_spring){
	//as distance decrease increase stiffness
	//use lennard jones
	//VLJ(r) = 4 e [ (t / r)12 - (t / r)6 ] 
	float sigma = g_params.radius;// finite distance at which the inter-particle potential is zero
	float e = 1.0f;    // depth of the potential well
	int i = 0;
	int j = 0;
	int s = 0;
	float r = 0.0f;
	for (int k = 0; k < indices_spring.size(); k++){
		s = indices_spring[k];
		if (s == -1) continue;
		i = g_buffers->springIndices[s * 2];
		j = g_buffers->springIndices[(s * 2) + 1];
		r = Length(Vec3(g_buffers->positions[i]) - Vec3(g_buffers->positions[j]));
		//g_buffers->springLengths[s] = r; 
		//cout << k << " " << s << " " << i << " " << j << " " << r ;//0 192000 1283 5114 6.5972
		//sigma = g_buffers->springLengths[s];
		//g_buffers->springStiffness[s] = 4.0f*e*(pow((sigma/r),12.0f)-pow((sigma/r),6.0f));
		//if (r > 2.5f*sigma)
		//	r = 2.5f*sigma;
		if (r == 0.0)
			r = sigma;
		if (r < sigma)
			r = sigma;
		g_buffers->springStiffness[s] = fabsf(e*(Pow((sigma / r), 12.0f) - 2.0f*Pow((sigma / r), 6.0f)));
		if (r > 2.5f*sigma)
			g_buffers->springStiffness[s] = 0.0050f; //0.00025f;

		//cout << " stiffness " << g_buffers->springStiffness[s] << " " << r << " " << sigma << " " << g_buffers->springLengths[s] << endl;
		//if (r < sigma) 
		//	cout << "stiffness " << g_buffers->springStiffness[s] << " " << r << " " << sigma << endl;
	}
}

//replace ?
void ReplaceSpring(int i, int j, int k, int l){
	//go through all spring and find i,j ?
	// or generate a new array storing pair information id
	int sp_counter = 0;
	//printf ("search for %i and %i to be replace by %i and %i\n",i,j,k,l);
	for (int m = 0; m<g_buffers->springLengths.size(); m++){
		if (g_buffers->springIndices[sp_counter] == i){
			//printf ("found i  j is %i\n",g_buffers->springIndices[sp_counter+1]);
			if (g_buffers->springIndices[sp_counter + 1] == j){
				//printf ("found the pair i %i j %i\n",i,j);
				g_buffers->springIndices[sp_counter] = k;
				g_buffers->springIndices[sp_counter + 1] = l;
				break;
			}
		}
		/*else {
		if (g_buffers->springIndices[sp_counter] == j){
		printf ("found j but i is %i\n",g_buffers->springIndices[sp_counter+1]);
		if (g_buffers->springIndices[sp_counter+1] == i){
		printf ("found the pair j %i i %i\n",j,i);
		g_buffers->springIndices[sp_counter] = l;
		g_buffers->springIndices[sp_counter+1] = k;
		break;
		}
		}
		}*/
		sp_counter += 2;
	}
}

void addPointToRope(Rope& rope,int posid){
    //how to insert a point a pos posid in the rope ?
    //rope.mIndices is the list g_position Id for the rope
    //we need to break spring posid->posid+1
    //and create posid->newPoint, newPoint->posid+1
}

//this assumed orderer  indice and not reused one
void CreatePersistence(Rope& rope, int current, int persistence, float stiffness, int nfloat, float give, float D)
{
	for (int j = 1; j < persistence+1; j++){
		float r = Randf((D/100.0f)*-1.0f, 0.0f);
		if (rope.mIndices.size()>j)
			CreateSpringInter(current - j, current, stiffness, give, (D*(float)j) + r*(float)(j - 1));
	}
}

//this is using indices in rope indices not in particles indices
void CreatePersistenceRope(Rope& rope, int current, int persistence, float stiffness, 
	int nfloat, float give, float D, float hardness)
{
	for (int j = 1; j < persistence + 1; j++){
		//float r = Randf((D / 100.0f)*-1.0f, 0.0f);
		float r = Randf(-hardness, hardness);
		if (rope.mIndices.size()>j)
			CreateSpringInter(rope.mIndices[current - j], rope.mIndices[current], stiffness, give, (D*(float)j) + r);// *(float)(j - 1));
	}
}

void CreateClosedPersistence(Rope& rope, int start,int end, int persistence, float stiffness, int nfloat, float give, float D)
{
	//Lv1 = i+1
	//Lv2 = i+2
	//Lv3 = i+3
	//Lv4 = i+4
	for (int l = 0; l < persistence; l++){
		float r = Randf((D / 100.0f)*-1.0f, 0.0f);;// float r = Randf(-1.0f, 1.0f) / 2000.0f;
		for (int k = l ,  i = 0; k >= 0, i < l + 1 ; k--, i++)
		{
			CreateSpringInter(start + i, end - k, stiffness, give, (D*(float)(l + 1)) + r*(float)(l));
			cout << l << "start " << i << " end " << k << " " << l << " " << D << " " << D*(float)(l+1 ) << endl;
		}
	}
}

void CreateInsertedPersistence(Rope rope, int ip,  float stiffness, float give, float D)
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
			ReplaceSpring(a,b,k,l);
			//cout << a << " " << b << " " << k << " " << l << endl;
			//cout << lvl << " replace " << ip - (lvl - i) << " / " << ip + (i + 2) << " by " << ip - (lvl - i) << " / " << ip + (i + 1) << endl;
		}
		//cout << lvl << "create  " << ip + 1 << " / " << ip + (lvl + 2) << " D " << (D*(float)(lvl+1)) + r*(float)(lvl) << endl;
		CreateSpringInter(rope.mIndices[ip + 1], rope.mIndices[ip + (lvl + 2)], stiffness, give, (D*(float)(lvl + 1)) + r*(float)(lvl));
	}
}

void CreateRopeFromData(Rope& rope, float scale, float stiffness, float *data,
                        float length, int nfloat, int phase, float spiralAngle=0.0f,
                        float invmass=1.0f, float give=0.075f, int extend_nb =57,
                        bool extend = false, bool closed = true,float D=0.0f,int persistence = 2)
{
    int start = int(g_buffers->positions.size());
    if (give != 0.0f) give=0.0f;
    float r=1.0f;//biased on the 1-3 spring
    int current=0;
    //if closed do the last point ?
    for (int i=0; i < nfloat; i+=3)
    {
        //if (i/3>20) break;
        //printf("add a point %i %f %f %f \n",i/3,data[i]/1000.0f, data[i+1]/1000.0f, data[i+2]/1000.0f);
        int begin = int(g_buffers->positions.size());
        int prev = begin;//int(g_buffers->positions.size())-1;

        if (!extend){
            current = begin;
            rope.mIndices.push_back(begin);
            g_buffers->positions.push_back(Vec4(data[i]*scale, data[i+1]*scale, data[i+2]*scale, 1.0f));
            g_buffers->velocities.push_back(0.0f);
            g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));

			CreatePersistence(rope, begin, persistence, stiffness, nfloat, give, D);

            //if (rope.mIndices.size()>1)
            //    CreateSpring(begin-1, begin, stiffness, give,D);
            //else
            //    CreateSpring(begin, start, stiffness, give,D);
            //if (rope.mIndices.size() > 2){
            //    r=Randf(-1.0f,1.0f)/2000.0f;
            //    if ((i+3 < nfloat)) {
            //        CreateSpring(begin-2, begin, stiffness, give,(D*2.0f)+r);//*0.5f
            //    }
            // }
        }
        else {
            Vec3 point = Vec3(data[i]*scale, data[i+1]*scale, data[i+2]*scale);
            Vec4 next_point = Vec4(0,0,0,0);
            Vec3 dirtopt = Vec3(0,0,0);
            if (i+5 < nfloat) {
                next_point = Vec4(data[i+3]*scale, data[i+4]*scale, data[i+5]*scale,1.0f);
            }
            else {
                next_point = Vec4(data[0]*scale, data[1]*scale, data[2]*scale,1.0f);
            }
            dirtopt = Vec3 (next_point.x - point.x,next_point.y - point.y,next_point.z - point.z);

            for (int j=0;j<extend_nb;j++){
                current = int(g_buffers->positions.size());
                float perc = ((float)j/(float)extend_nb);
                rope.mIndices.push_back(int(g_buffers->positions.size()));

                g_buffers->positions.push_back(Vec4( point.x+dirtopt.x*perc,
                                            point.y+dirtopt.y*perc,
                                            point.z+dirtopt.z*perc,
                                            1.0f));
                g_buffers->velocities.push_back(0.0f);
                g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
                //printf ("subdivide how many %i %i %i %i\n",i,j,extend_nb,current);
				CreatePersistence(rope, current, persistence, stiffness, nfloat, give, D);
				/*
                if (rope.mIndices.size() > 1) {
                    CreateSpring(current-1, current, stiffness, give,D);//current +1 doesnt exist yet
                }
                //printf ("i>3? %i %i\n",current,current-2);
                if (rope.mIndices.size() > 2){
                    r=Randf(-1.0f,1.0f)/2000.0f;
                    //printf ("bias is %f\n",r);
                    //if ((j > 2)) {
                    CreateSpring(current-2, current, stiffness, give,(D*2.0f)+r);
                    //}s
                }*/
            }
            //should fixed he point with a str
        }
        // close
        // tether
        //if (i > 0 && i%4 == 0)
        //CreateSpring(prev-3, prev+1, -0.25f);

        // bending spring
        //if (i > 0)
        //	CreateSpring(prev-1, prev+1, stiffness*0.5f, give,D*2.0);

    }
    if (closed) {
		r = Randf(-1.0f, 1.0f) / 2000.0f;
        //printf ("bias is %f\n",r);
		int startindex = rope.mIndices[0];
		int endindex = int(g_buffers->positions.size()-1);// rope.mIndices[rope.mIndices.size() - 1];
		//CreateSpring(endindex, startindex, stiffness, give, D);
		CreateClosedPersistence(rope, startindex, endindex, persistence, stiffness, nfloat, give, D);		
		//CreateSpring(endindex, startindex + 1, stiffness, give, (D*2.0f) + r);
		//CreateSpring(endindex - 1, startindex, stiffness, give, (D*2.0f) + r);
    }
    rope.coarseIndices = rope.mIndices;
    //CreateSpring(int(g_buffers->positions.size())-1,0, stiffness*0.5f, give);
    //CreateSpring(int(g_buffers->positions.size()),0, stiffness*0.5f, give);
}


Quat AlignVec3s(Vec3 fixed, Vec3 moving)
{
	//return rotation aligning moving to fixed
	Vec3 axis = Cross(moving / Length(moving), fixed / Length(fixed));
	float dot = Dot(moving / Length(moving), fixed / Length(fixed));
	float angle = ACos(dot);
	return QuatFromAxisAngle(axis, angle);
}

int HaltonSample3D(float* radius, float separation, Vec3* points, int maxPoints){
    int N=maxPoints;
    int DIM_MAX = 3;
    int* base = new int [DIM_MAX];
    int i;
    int j;
    int* leap = new int [DIM_MAX];
    int n;
    int dim_num;
    double* r = new double[DIM_MAX*N];
    int* seed = new int[DIM_MAX];
    int step;
    halton_dim_num_set ( DIM_MAX );
    step = 0;
    halton_step_set ( step );
    for ( i = 0; i < DIM_MAX; i++ )
    {
        seed[i] = 0;
    }
    halton_seed_set ( seed );
    for ( i = 0; i < DIM_MAX; i++ )
    {
        base[i] = prime ( i + 1 );
    }
    halton_base_set ( base );

    i4vec_transpose_print ( DIM_MAX, seed, "  SEED = " );
    i4vec_transpose_print ( DIM_MAX, base, "  BASE = " );
    halton_sequence ( N, r );
    for ( j = 0; j < N; j++ )
    {
        for ( i = 0; i < DIM_MAX; i++ )
        {
            points[j][i] = r[i+j*DIM_MAX]*radius[i];//scale ?
        }
    }
	return 0;
}

class cellPACK
{
public:
    Json::Value book_json;
	Json::Value results_json;
	//const string datapath="/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/other/";
    //const string geompath="/opt/data/dev/cellPACK/cellPACK_data_git/cellPACK_database_1.1.0/geometries/";
	string mainpath = "..\\..\\data\\cellpack\\";// "D:\\Data\\cellPAC_data\\cellPACK_database_1.1.0\\";//D:\Data\cellPAC_data\cellPACK_database_1.1.0
	string datapath = "..\\..\\data\\cellpack\\beads\\";//"D:\\Data\\cellPAC_data\\cellPACK_database_1.1.0\\other\\";
	string geompath = "..\\..\\data\\cellpack\\geoms\\";//"D:\\Data\\cellPAC_data\\cellPACK_database_1.1.0\\geometries\\";

	float main_scale = 1.0f;
    float main_radius=10.0f;
    int iGroupCounter;
    Vec3 minExtents, maxExtents; //boudningBox
	NvFlexExtContainer* fcontainer;
    int maxParticles;	 
    int numParticles = 0;
    std::vector<int> mask;
    std::vector<int> maks_protein;
    std::vector<int> maks_fiber;
	std::vector<int> mask_membrane;
    float rope_unit_length ;//= g_params.mRadius;
    int rope_phase;
    std::vector<char>  result_curve;
    int nfloat;
    float* data_curve ;
    Rope r_dna;
    
	float mClusterSpacing = 5.75f;
	float mClusterRadius = 10.0f;
	float mClusterStiffness = 0.15f;
	float mLinkRadius = 10.0f;
	float mLinkStiffness = 1.0f;

    std::vector<string> pnames;
    std::vector<string> pnames_fiber;
    
	bool use_rb=false;

    struct CompMask
    {
        std::vector<int> mask;
    };

    struct IngredientInstance
    {
        Vec3 mTranslation;
        Quat mRotation;
        float mLifetime;

        int mGroup;
        int mParticleOffset;
        int mMeshIndex;
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

	std::vector<NvFlexExtAsset> iBatches;
	std::vector<NvFlexExtAsset> iBatchesFiber;
	std::vector<IngredientSphereTree> mIngrSphereTree;
    std::vector<IngredientInstance> mInstances;
	std::vector<NvFlexDistanceFieldId*> comp_shape;//SDF

    std::vector<Mesh*> comp_mesh;
	std::vector<ClothMesh*> comp_cloth;
    std::vector<CompMask*> comp_mask;
    CompMask* compmask;
	NvFlexExtForceField* forceFields;


	cellPACK(){
		GetCurrentDir(cCurrentPath, sizeof(cCurrentPath));
		cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */
		string path = string(cCurrentPath) + "\\..\\..\\data\\cellpack\\";
		//std::cout << "The current working directory is " << cCurrentPath << " " << path << endl;
		//ChangeDir( path.c_str());
		//GetCurrentDir(cCurrentPath, sizeof(cCurrentPath));
		//std::cout << "The current working directory is " << cCurrentPath << endl;
		mainpath = path;
		datapath = mainpath + "\\beads\\";
		geompath = mainpath + "\\geoms\\";
	};

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
				cout << "found nbLevel " << ingr_spheres.nbLevel << " " << l << endl;
				for (int i = 0; i < ingr_spheres.nbLevel; i++){
					int tmp;
					ifs.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
					ingr_spheres.LevelCounts.push_back(tmp);
					cout << "parsing " << tmp << endl;
				}
				//do we have binding info
				int b;
				ifs.read(reinterpret_cast<char*>(&b), sizeof(int));
				ingr_spheres.nBinding = b;
				cout << "found nBinding " << ingr_spheres.nBinding << " " << b << endl;
				int start = 0;
				for (int i = 0; i < ingr_spheres.nBinding; i++){
					int count;
					ifs.read(reinterpret_cast<char*>(&count), sizeof(count));
					ingr_spheres.BindingStarts.push_back(start);
					ingr_spheres.BindingStarts.push_back(count);
					start += count;
				}
				//gather coordinates
				for (int i = 0; i < ingr_spheres.nbLevel; i++){
					ingr_spheres.LevelStarts.push_back(ingr_spheres.LevelPoints.size());
					for (int j = 0; j < ingr_spheres.LevelCounts[i]; j++){
						Vec4 tmp;
						ifs.read(reinterpret_cast<char*>(&tmp[0]), sizeof(float) * 4);
						ingr_spheres.LevelPoints.push_back(Vec3(tmp.x,tmp.y,tmp.z));
						cout << tmp.x << " " << tmp.y << " " << tmp.z << endl;
						//cout << "parsing " << ingr_spheres.LevelCounts[i] << " points nb " << ingr_spheres.LevelPoints.size() << endl;
					}
					cout << "parsing " << ingr_spheres.LevelCounts[i] << " points nb " << ingr_spheres.LevelPoints.size() << endl;
				}
				//gather mapping
				/*
				mapping.extend(mappingL1_atom_order)#L1_size *2
				mapping.extend(mappingL2_atom_order)#L2_size *2
				mapping.extend(mappingL2_L1_order)#L2_size *2
				*/
				for (int i = 0; i < ingr_spheres.nbLevel - 1; i++){
					//get atom mapping
					ingr_spheres.LevelMappingStarts.push_back(ingr_spheres.LevelMapping.size());
					for (int j = 0; j < ingr_spheres.LevelCounts[i]; j++){
						int start;
						ifs.read(reinterpret_cast<char*>(&start), sizeof(int));
						int count;
						ifs.read(reinterpret_cast<char*>(&count), sizeof(int));
						ingr_spheres.LevelMapping.push_back(Vec2(start,count));
						cout << i <<" " << j << " start " << start << " count " << count << endl;
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
							cout << i << " j " << j << " k " << k << " start " << start << " count " << count << endl;
						}
					}
				}
				//gather bindingsite info if any
				if (ingr_spheres.nBinding != 0){
					int total = 1;
					for (int i = 0; i < ingr_spheres.nBinding; i++){
						int count = ingr_spheres.BindingStarts[i*ingr_spheres.nBinding + 1];
						cout << " count " << count << endl;
						//int start = ingr_spheres.BindingStarts[i*ingr_spheres.nBinding + 0];
						for (int j = 0; j < count; j++){
							int tmp;
							ifs.read(reinterpret_cast<char*>(&tmp), sizeof(int));
							ingr_spheres.BindingSites.push_back(tmp);
						}
						total *= count;
					}
					cout << " total distance number " << total << endl;
					//gather the distance matrix
					for (int i = 0; i < total; i++){
						float tmp;
						ifs.read(reinterpret_cast<char*>(&tmp), sizeof(float));
						ingr_spheres.DistancesMatrix.push_back(tmp);
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
	
	void buildMembrane(float stiffness, float mass, int phase, float L =0.0f)
	{
		//AddInflatable(cp->comp_mesh[0], 1.0f, cp->iGroupCounter++);
		// add particles to system
		int ph = NvFlexMakePhase(phase, eNvFlexPhaseSelfCollide);
		int start = g_buffers->positions.size();
		mask_membrane.push_back(start);
		for (size_t i = 0; i < comp_mesh[0]->GetNumVertices(); ++i)
		{
			const Vec3 p = Vec3(comp_mesh[0]->m_positions[i]) + Vec3(comp_mesh[0]->m_normals[i])*11.0f;
			
			g_buffers->positions.push_back(Vec4(p.x, p.y, p.z, 1.00f / mass));
			g_buffers->restPositions.push_back(Vec4(p.x, p.y, p.z, 0.0));

			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(ph); 
		}

		for (size_t i = 0; i <  comp_mesh[0]->m_indices.size(); i += 3)
		{
			int a =  comp_mesh[0]->m_indices[i + 0];
			int b =  comp_mesh[0]->m_indices[i + 1];
			int c =  comp_mesh[0]->m_indices[i + 2];

			float La = Length(Vec3(g_buffers->positions[start + a]) - Vec3(g_buffers->positions[start + b]));
			float Lb = Length(Vec3(g_buffers->positions[start + b]) - Vec3(g_buffers->positions[start + c]));
			float Lc = Length(Vec3(g_buffers->positions[start + c]) - Vec3(g_buffers->positions[start + a]));

			if (L == 0.0f)
				L = La;

			CreateSpringInter(start + a, start + b, stiffness, 0.0f, (L == 0.0f) ? La : L);// main_radius*2.0f);
			CreateSpringInter(start + b, start + c, stiffness, 0.0f, (L == 0.0f) ? Lb : L);
			CreateSpringInter(start + c, start + a, stiffness, 0.0f, (L == 0.0f) ? Lc : L);

		}

		start = g_buffers->positions.size();
		ph = NvFlexMakePhase(phase+1, eNvFlexPhaseSelfCollide);
		for (size_t i = 0; i < comp_mesh[0]->GetNumVertices(); ++i)
		{
			const Vec3 p = Vec3(comp_mesh[0]->m_positions[i]) - Vec3(comp_mesh[0]->m_normals[i])*11.0f;

			g_buffers->positions.push_back(Vec4(p.x, p.y, p.z, 1.00f / mass));
			g_buffers->restPositions.push_back(Vec4(p.x, p.y, p.z, 0.0));

			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(ph);
		}

		for (size_t i = 0; i < comp_mesh[0]->m_indices.size(); i += 3)
		{
			int a = comp_mesh[0]->m_indices[i + 0];
			int b = comp_mesh[0]->m_indices[i + 1];
			int c = comp_mesh[0]->m_indices[i + 2];
			CreateSpringInter(start + a, start + b, stiffness, 0.0f, Length(Vec3(g_buffers->positions[start + a]) - Vec3(g_buffers->positions[start + b])));// main_radius*2.0f);
			CreateSpringInter(start + b, start + c, stiffness, 0.0f, Length(Vec3(g_buffers->positions[start + b]) - Vec3(g_buffers->positions[start + c])));
			CreateSpringInter(start + c, start + a, stiffness, 0.0f, Length(Vec3(g_buffers->positions[start + c]) - Vec3(g_buffers->positions[start + a])));
		}

		mask_membrane.push_back(comp_mesh[0]->GetNumVertices()*2);
	}

    void compartmentsSDF(Json::Value comp){
        Json::Value comp_geom = comp["geom"].asString();//at 0 ?
        Json::Value comp_name = comp["name"].asString();//at 0 ?
		//comp_geom = "MMycoideHD.dae";
		if (comp_geom == 0) return; 
		string path = "Mmycoides_PackingSurf_1_Med.obj";// comp["geom"].asString();// "MMycoideHD.dae"; //comp["geom"].asString();
		std::cout << "path  : " << path << " compId " << comp_name << endl;
        string objpath = string(path.c_str(), strrchr(path.c_str(), '.')) + ".obj";
		Mesh* mesh = GetMesh(GetFilePathByPlatform((geompath+objpath).c_str()).c_str(),main_scale);
		//SDF* sdf = CreateSDFfromMesh(GetFilePathByPlatform((geompath + objpath).c_str()).c_str(), mesh,
        //          main_scale, Vec3(0.0f, 0.0f, 0.0f),minExtents , maxExtents,0.0f);//main_radius/5.0f);
		const int dim = 128;
		NvFlexDistanceFieldId sdf = CreateSDF(GetFilePathByPlatform((geompath + objpath).c_str()).c_str(), dim);
		AddSDF(sdf, Vec3(-1.f, 0.0f, 0.0f), QuatFromAxisAngle(Vec3(0.0f, 1.0f, 0.0f), DegToRad(-45.0f)), 0.5f);

        //comp_shape.push_back(sdf); 
        comp_mesh.push_back(mesh);
		
		//FlexTriangleMesh* fmesh = CreateTriangleMesh(mesh);
		//AddTriangleMesh(fmesh, Vec3(), Quat(), 1.0f);

		//AddSDF(sdf->fsdf, Vec3(0.0f, 0.0f, 0.0f), Quat(), main_scale);
        //make a forcefield from it ?
		//create cloths rigid ? MMycoideHD
    }

	void compartmentsSDF(string objpath){
		Mesh* mesh = GetMesh(GetFilePathByPlatform((geompath + objpath).c_str()).c_str(), main_scale);
		mesh->GetBounds(minExtents, maxExtents);
		cout << main_scale << endl;
		cout << minExtents.x << " " << minExtents.y << " " << minExtents.z << endl;
		cout << maxExtents.x << " " << maxExtents.y << " " << maxExtents.z << endl;
		const int dim = 128;
		NvFlexDistanceFieldId sdf = CreateSDF(GetFilePathByPlatform((geompath + objpath).c_str()).c_str(), dim);
		AddSDF(sdf, Vec3(-1.f, 0.0f, 0.0f), QuatFromAxisAngle(Vec3(0.0f, 1.0f, 0.0f), DegToRad(-45.0f)), 0.5f);

		//SDF* sdf = CreateSDFfromMesh(GetFilePathByPlatform((geompath + objpath).c_str()).c_str(), mesh,
		//	main_scale, Vec3(0.0f, 0.0f, 0.0f), minExtents, maxExtents, 0.0f);//main_radius/5.0f);
		//comp_shape.push_back(sdf);
		comp_mesh.push_back(mesh);
	}

    void sdfToForceField()
    {
		forceFields = new NvFlexExtForceField[comp_shape.size()];
		/*
        //for every surface ingredient use sdf to snap
        //flexSetFields(g_flex, &cp->comp_shape[0], cp->comp_shape.size());
        for (int i=0;i<comp_shape.size();i++){
			NvFlexDistanceFieldId * fsdf = comp_shape[i];
            //forceFields[i].sdf = fsdf;
            forceFields[i].dim = fsdf->mDepth;
            //forceFields[i].nInstance = mask.size();
            //forceFields[i].mGroupIndices = &mask[0];
            forceFields[i].spacing = (fsdf->mUpper[0] - fsdf->mLower[0]) / (float)fsdf->mWidth;
            forceFields[i].mField = fsdf->mField;
            forceFields[i].mRadius =(float)mask.size();
            forceFields[i].mPosition[0] = fsdf->mLower[0];
            forceFields[i].mPosition[1] = fsdf->mLower[1];
            forceFields[i].mPosition[2] = fsdf->mLower[2];
            forceFields[i].mPosition2[0] = fsdf->mUpper[0];
            forceFields[i].mPosition2[1] = fsdf->mUpper[1];
            forceFields[i].mPosition2[2] = fsdf->mUpper[2];
			forceFields[i].mStrength = 10.0f;
            //cout << fsdf->mUpper << endl;
            //cout << fsdf->mLower << endl;
            //cout << fsdf->mWidth << endl;
            cout << " ? spacing ? " << forceFields[i].spacing << endl;
			cout << " ? mField[0] ? " << forceFields[i].mField[0] << endl;
        }
        //FlexSDF sdf[comp_shape.size()];
        //sdf[0] = *comp_shape[0];
        //sdf[1] = *comp_shape[1];
        //fcontainer->nInstance = (int)mask.size();
        //fcontainer->mGroupIndices = &mask[0];
        //also pass the iBatches
        //for (int i=0;i<maks_protein.size();i++){
        //    cout << "insta " << maks_protein[i] << "comp " << iBatches[maks_protein[i]].compId << endl;
        //    }
        */
    }

	void setForceField()
	{
		//need another function to setup with only membrane
		//flexExtSetForceFields(fcontainer, forceFields, &iBatches[0], iBatches.size(), &iBatchesFiber[0], iBatchesFiber.size(),
		//	&maks_protein[0], comp_shape.size(), eFlexMemoryHost);
	}

    void ingredientToFlexAsset(Json::Value ingr_node, int compId, bool rb=false)
    {
		/*
        string proxyname;
        Json::Value ingr_source = ingr_node["source"];

        std::cout << "ingr name : " << ingr_node["name"].asString() << endl;

        string ingr_source_pdb = ingr_source["pdb"].asString();
        std::cout << "ingr pdb  : " << ingr_source_pdb << " compId " << compId << endl;
        if (ingr_source_pdb.size() == 4 )
        {
            proxyname = datapath+ingr_source_pdb+".pdb_kmeans3.txt";
        }
        else
        {
            proxyname = datapath+ingr_source_pdb+".pdb_kmeans3.txt";
        }
        std::cout << "ingr proxy  : " << proxyname << endl;

        Vec3 pcpalVector = Vec3(ingr_node["principalVector"][0].asFloat(),
                                ingr_node["principalVector"][1].asFloat(),
                                ingr_node["principalVector"][2].asFloat());
		Json::Value  offsetnode = ingr_source["transform"]["offset"];
		Vec3 offset = Vec3(0.0f,0.0f,0.0f);//should be ingr_node["source"]["transform"]["offset"] if exist
		if (offsetnode != 0) {
			offset = Vec3(offsetnode[0].asFloat()*main_scale,
				offsetnode[1].asFloat()*main_scale,
				offsetnode[2].asFloat()*main_scale);
		}

        //parse the proxy and create a rigid body asset
        //check if "Type":"Grow"
		IngredientSphereTree ingr_spheres;
        if (ingr_node["Type"].asString() == "Grow")
        {
            //model as rope in flex
            pnames_fiber.push_back(ingr_node["name"].asString());
			NvFlexExtAsset b;//for ribid body protein
            //b.mAsset = flexExtCreateRigidFromPoints(points);
            //b.mMesh = nullptr;//CreateGpuMesh(mesh);
            //b.compId = compId;
            //b.offsetx = offset.x;
			//b.offsety = offset.y;
			//b.offsetz = offset.z;
            //b.pcpalVectorx = pcpalVector.x;
            //b.pcpalVectory = pcpalVector.y;
            //b.pcpalVectorz = pcpalVector.z;
            //b.nInstances = 0;
            //b.ingr_name = ingr_node["name"].asString().c_str();
            iBatchesFiber.push_back(b);  
        }
        else {
			//use proxy file or positions
			std::vector<Vec3> points_to_use;
			Json::Value  sphereFile = ingr_node["sphereFile"];
			
			std::cout << "try loading " << ingr_node["sphereFile"] <<" "<< sphereFile.asString() << " " << datapath + sphereFile.asString() << endl;
			string filename = sphereFile.asString();

			if (filename.substr(filename.find_last_of(".") + 1) == "sph")
				filename = ingr_source_pdb + ".pdb_kmeans3.txt";
			if (!sphereFile.empty())
			{
				std::cout << "load " << filename << endl;
				ingr_spheres = parseProxy(datapath + filename);
				//use lvl0
				for (int i = 0; i < ingr_spheres.LevelCounts[0];i++){
					points_to_use.push_back(ingr_spheres.LevelPoints[i] * main_scale);
				}
				std::cout << "ingr proxy nb " << points_to_use.size() << endl;
			}
			else {
				std::cout << "load from positions" << ingr_node["positions"].size() << ingr_node["positions"][0].size() << endl;
				Json::Value jsonpos = ingr_node["positions"][1];
				for (int i = 0; i < jsonpos.size(); i++){
					ingr_spheres.LevelPoints.push_back(Vec3(jsonpos[i][0].asFloat(), jsonpos[i][1].asFloat(), jsonpos[i][2].asFloat())*main_scale);
					//std::cout << points[i].x << " " << points[i].y << " " << points[i].z << endl;
				}
				std::cout << "ingr positions  nb " << ingr_spheres.LevelPoints.size() << endl;
				points_to_use = ingr_spheres.LevelPoints;
			}

			if (points_to_use.size() == 0) {
				std::cout << "ingr proxy  0 " << points_to_use.size() << endl;
                return;
            }

			NvFlexExtAsset b;
			if (use_rb) b.mAsset = flexExtCreateRigidFromPoints(points_to_use);
			//const float* points, int numParticles, float clusterSpacing, float clusterRadius, float clusterStiffness, float linkRadius, float linkStiffness
			else {
				b.mAsset = flexExtCreateSoftFromPoints((float*)&points_to_use[0], points_to_use.size(),
					0.10f,
					0.10f,
					0.05f,
					g_params.mRadius,// / 2.0f,
					1.00f);
				//overwrite the string
			}

            //b.mMesh = nullptr;//CreateGpuMesh(mesh);	

            b.compId = compId;
            b.offsetx = offset.x;
            b.offsety = offset.y;
            b.offsetz = offset.z;
            b.pcpalVectorx = pcpalVector.x;
            b.pcpalVectory = pcpalVector.y;
            b.pcpalVectorz = pcpalVector.z;
            b.nInstances = 0;
            b.ingr_name = ingr_node["name"].asString().c_str();
            iBatches.push_back(b);
            pnames.push_back(ingr_node["name"].asString());
            
            std::cout << "ingr index  " << iBatches.size()-1 << " " << b.ingr_name << endl;
			std::cout << "particle count " << b.mAsset->mNumParticles << endl;
			std::cout << "nShape " << b.mAsset->mNumShapes << endl;
            }
		mIngrSphereTree.push_back(ingr_spheres);
		*/
    }

    void parseJsonIngredients(Json::Value  ingr_nodes, int comp){
        for( Json::ValueIterator itr = ingr_nodes.begin() ; itr != ingr_nodes.end() ; itr++ ) {
            Json::Value ingr_node_name = *itr;
			//if (ingr_node_name["name"].asString() != "Ves203MeshsParent")
			//		continue;
            ingredientToFlexAsset(ingr_node_name,comp);
            string ingr_name = ingr_node_name["name"].asString();
            Json::Value ingr_node_results = ingr_node_name["results"];
            //JSONNode ingr_node_results_data =  ingr_node_results.as_array();//array or array
            if (ingr_node_results.size()!=0)
            {//do something
                std::cout << "array size ? #nb of instances" << ingr_node_results.size()<< endl;
            }
        }
    }

	void parseJsonIngredientsResults(Json::Value  ingr_nodes, int comp,bool redo=false){
		for (Json::ValueIterator itr = ingr_nodes.begin(); itr != ingr_nodes.end(); itr++) {
			Json::Value ingr_node_name = *itr;
			string ingr_name = ingr_node_name["name"].asString();
			Json::Value ingr_node_results = ingr_node_name["results"];
			//JSONNode ingr_node_results_data =  ingr_node_results.as_array();//array or array
			//if (ingr_name != "DNA-binding protein HU")
			//	continue;
			if (ingr_node_results.size() != 0)
			{//do something
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
					else updateInstanceIngredient(ingrIndex, pos, quat,i);
					//if (i == 2) return;
					//break;
				}
			}
		}
	}

    Json::Value getJsonIngredients(Json::Value  ingr_nodes, string name){
        Json::Value ingr_node_name;
        for( Json::ValueIterator itr = ingr_nodes.begin() ; itr != ingr_nodes.end() ; itr++ ) {
            ingr_node_name = *itr;
            string ingr_name = ingr_node_name["name"].asString();
            if (ingr_name == name) 
                return ingr_node_name;
        }
        return 0;
    }
    
    int getIngredientBatchId(string name){
        int index = 0;
		bool found = false;
        for (int i=0;i <pnames.size();i++)
        {
            if ( pnames[i] == name ) {
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
    
	NvFlexExtAsset getIngredientBatch(string name){
		int index = 0;
		bool found = false;
		NvFlexExtAsset b;
		
		for (int i = 0; i <pnames.size(); i++)
		{
			if (pnames[i] == name) {
				index = i;
				found = true;
				b = iBatches[i];
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
					b = iBatchesFiber[i];
					break;
				}
			}
		}
		return b;
	}

    Json::Value getIngredients(string name){
        Json::Value  cyto = book_json["cytoplasme"];
        Json::Value ingr_node_name;
        if (cyto != 0)
        {
            Json::Value cyto_ingredients = cyto["ingredients"];
            if (cyto_ingredients != 0) {
                ingr_node_name = getJsonIngredients(cyto_ingredients,name);
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
                if (comp_name.size() == 0 ) continue;
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
                if (comp_name.size() == 1 ) continue;
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

    void loadRecipe(string filename, bool ignore_comp=false )
    {

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
        minExtents = Vec3(bbox[0][0].asFloat(),bbox[0][1].asFloat(),bbox[0][2].asFloat())*main_scale;
        maxExtents = Vec3(bbox[1][0].asFloat(),bbox[1][1].asFloat(),bbox[1][2].asFloat())*main_scale;

        Json::Value  cyto = book_json["cytoplasme"];
        if (cyto != 0)
        {
            std::cout << "find compartments cyto " << cyto.size() << endl;
            Json::Value cyto_ingredients = cyto["ingredients"];
            std::cout << "compartment cyto should have n ingredients " << cyto_ingredients.size() << endl;
			if (cyto_ingredients != 0)
			{
				 parseJsonIngredients(cyto_ingredients,0);
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
                if (comp_name.size() == 0 ) continue;
                Json::Value comp_geom = comp_name["geom"];//at 0 ?
                if (comp_geom != 0)
                {
					if (!ignore_comp)compartmentsSDF(comp_name);
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
                if (comp_name.size() == 1 ) continue;
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
        iGroupCounter=0;
        maxParticles = 1024*1024;
    }

    void createFlexInstanceIngredient(int ingrIndex, Vec3 position, Quat rotation)
    {
        //int particleOffset = g_buffers->positions.size();
        /*FlexExtAsset* asset = iBatches[ingrIndex].mAsset;

        Matrix33 rot = Matrix33(rotation);
        Matrix44 transform = Matrix44(Vec4(rot.cols[0]),Vec4(rot.cols[1]),Vec4(rot.cols[2]),Vec4(position));

        //transform.SetTranslation(Point3(position.x,position.y,position.z));
        int phase = NvFlexMakePhase( iGroupCounter, 0);
        iGroupCounter++;
        float invMassScale = 1.0f;

        FlexExtInstance* ingr_inst = flexExtCreateInstance(fcontainer, asset , transform, 0.0f, 0.0f, 0.0f, phase,invMassScale );
        //ingr_inst->mUserData = &phase;

        numParticles+=asset->mNumParticles;
		*/
    }

	void updateInstanceIngredient(int ingrIndex, Vec3 position, Quat rotation, int i){
		/*
		IngredientInstance& inst = mInstances[i];
		FlexExtAsset* asset = iBatches[inst.mMeshIndex].mAsset;
		for (int j = 0; j < asset->mNumParticles; ++j)
		{
			Vec3 localPos = Vec3(&asset->mParticles[j * 4]);//local position of the proxy
			Vec3 rpos = Rotate(rotation, localPos);
			g_buffers->positions[inst.mParticleOffset + j] = Vec4(position + rpos, 1.0f);
			g_buffers->velocities[inst.mParticleOffset + j] = Vec3(0, 0, 0);
		}
		*/
	}

	int createInstanceIngredient(int ingrIndex, Vec3 position, Quat rotation, bool rb = false)
    {
		/*
		const int particleOffset = g_buffers->positions.size();
		const int indexOffset = g_rigidOffsets.back();

		FlexExtAsset* asset = iBatches[ingrIndex].mAsset;
		// check we can fit in the container
		//if (int(g_buffers->positions.size()) - particleOffset < asset->mNumParticles)
		//    break;

		Quat q = rotation;// Quat();

		IngredientInstance inst;

		inst.mLifetime = 1000.0f;// Randf(5.0f, 60.0f);
		inst.mParticleOffset = particleOffset;
		inst.mRotation = Quat(0, 0, 0, 1);// q;// Quat(0, 0, 0, 1);// 

		inst.mTranslation = Vec3(position.x, position.y, position.z);
		
		inst.mMeshIndex = ingrIndex;

		Vec3 linearVelocity = Vec3(0.0f);//g_emitters[0].mDir*15.0f;//*Randf(5.0f, 10.0f);//Vec3(Randf(0.0f, 10.0f), 0.0f, 0.0f);
		Vec3 angularVelocity = Vec3(0.0f);//Vec3(UniformSampleSphere()*Randf()*k2Pi);

		mask.push_back(iBatches[ingrIndex].compId);
		maks_protein.push_back(ingrIndex);
		inst.mGroup = iGroupCounter++;

		//int mm = 2 ^ 24;
		
		int phase = NvFlexMakePhase(inst.mGroup, 0);
		if (!use_rb) phase = NvFlexMakePhase(inst.mGroup, eFlexPhaseSelfCollide | eFlexPhaseSelfCollideFilter);
		for (int j = 0; j < asset->mNumParticles; ++j)
		{
			Vec3 localPos = Vec3(&asset->mParticles[j * 4]);// - Vec3(&asset->mShapeCenters[0]);
			Vec3 rpos = Rotate(inst.mRotation, localPos);
			
			g_buffers->positions.push_back(Vec4(inst.mTranslation + rpos, 1.0f));//inst.mRotation*
			g_buffers->velocities.push_back(Vec3(0.0f));//linearVelocity + Cross(angularVelocity, localPos);
			g_buffers->phases.push_back(phase);
		}
		
		if (!use_rb){
			// add shape data to solver
			//cout << "create soft body with " << asset->mNumShapeIndices << " mNumShapeIndices" << endl; //38 ? should 100
			for (int i = 0; i < asset->mNumShapeIndices; ++i){
				//cout << "create soft body with " << asset->mShapeIndices[i] + particleOffset << " particleOffset" << endl;
				g_rigidIndices.push_back(asset->mShapeIndices[i] + particleOffset);// particle indices
			}
			 
			//cout << "create soft body with " << asset->mNumShapes << " mNumShapes" << endl;
			for (int i = 0; i < asset->mNumShapes; ++i)
			{
				//cout << "create soft body with " << asset->mShapeOffsets[i] + indexOffset << " indexOffset" << endl;
				g_rigidOffsets.push_back(asset->mShapeOffsets[i] + indexOffset);//rigid body indices
				g_rigidTranslations.push_back(Vec3(&asset->mShapeCenters[i * 3]) + inst.mTranslation);
				g_rigidRotations.push_back(inst.mRotation);
				g_rigidCoefficients.push_back(asset->mShapeCoefficients[i]);
			}

			// add link data to the solver 
			for (int i = 0; i < asset->mNumSprings; ++i)
			{
				g_buffers->springIndices.push_back(asset->mSpringIndices[i * 2 + 0] + particleOffset);
				g_buffers->springIndices.push_back(asset->mSpringIndices[i * 2 + 1] + particleOffset);

				g_buffers->springStiffness.push_back(asset->mSpringCoefficients[i]);
				g_buffers->springLengths.push_back(asset->mSpringRestLengths[i]);
			}
		}
		else {
		
			g_rigidCoefficients.push_back(0.15f);
			g_rigidTranslations.push_back(inst.mTranslation);
			g_rigidRotations.push_back(inst.mRotation);

			for (int j = 0; j < asset->mNumShapeIndices; ++j)
			{
				g_rigidLocalPositions.push_back(Vec3(&asset->mParticles[j * 4]));// -Vec3(&asset->mShapeCenters[0]));
				g_rigidIndices.push_back(asset->mShapeIndices[j] + particleOffset);
			}

			g_rigidOffsets.push_back(g_rigidIndices.size());
		}

		//particleOffset += asset->mNumParticles;

		mInstances.push_back(inst);
		iBatches[ingrIndex].nInstances++;
		if (use_rb) UpdateInstanceTransform(mInstances.size() - 1, position, rotation);
		return particleOffset;
		*/
    }

	/*void CompactObjects(bool rb=false)
	{
		//PyGILState_Release(gstate);
		// compact instances
		static std::vector<Vec4> particles(g_buffers->positions.size());
		static std::vector<Vec3> velocities(g_buffers->velocities.size());
		static std::vector<int> phases(g_buffers->phases.size());

		g_rigidTranslations.resize(0);
		g_rigidRotations.resize(0);
		g_rigidCoefficients.resize(0);
		g_rigidIndices.resize(0);
		g_rigidLocalPositions.resize(0);
		g_rigidOffsets.resize(0);

		// start index
		g_rigidOffsets.push_back(0);

		int numActive = 0;
		int particleOffset = g_buffers->positions.size();
		int indexOffset = g_rigidOffsets.back();


		for (int i = 0; i < int(mInstances.size()); ++i)
		{
			IngredientInstance& inst = mInstances[i];

			FlexExtAsset* asset = iBatches[inst.mMeshIndex].mAsset;
			particleOffset = g_buffers->positions.size();

			for (int j = 0; j < asset->mNumParticles; ++j)
			{
				particles[numActive + j] = g_buffers->positions[inst.mParticleOffset + j];
				velocities[numActive+j] = g_buffers->velocities[inst.mParticleOffset+j];
				phases[numActive + j] = g_buffers->phases[inst.mParticleOffset + j];
			}

			if (use_rb){
				g_rigidCoefficients.push_back(1.0f);
				g_rigidTranslations.push_back(Vec3(0,0,0));// inst.mTranslation);
				g_rigidRotations.push_back(Quat(0,0,0,1));// inst.mRotation);

				for (int j = 0; j < asset->mNumShapeIndices; ++j)
				{
					g_rigidLocalPositions.push_back(Vec3(&asset->mParticles[j * 4]) - Vec3(&asset->mShapeCenters[0]));//should we rotate ?
					g_rigidIndices.push_back(asset->mShapeIndices[j] + numActive);
				}

				g_rigidOffsets.push_back(g_rigidIndices.size());
			}
			else {
				indexOffset = g_rigidOffsets.back();
				for (int i = 0; i < asset->mNumShapeIndices; ++i)
					g_rigidIndices.push_back(asset->mShapeIndices[i] + particleOffset);

				for (int i = 0; i < asset->mNumShapes; ++i)
				{
					g_rigidOffsets.push_back(asset->mShapeOffsets[i] + indexOffset);
					g_rigidTranslations.push_back(Vec3(&asset->mShapeCenters[i * 3]) + inst.mTranslation);
					g_rigidRotations.push_back(inst.mRotation);
					g_rigidCoefficients.push_back(asset->mShapeCoefficients[i]);
				}
			}
			mInstances[i].mParticleOffset = numActive;

			// Draw transform
			//Matrix44 xform = TranslationMatrix(Point3(inst.mTranslation));// - inst.mRotation*Vec3(asset->mShapeCenters)))*RotationMatrix(inst.mRotation);
			//iBatches[inst.mMeshIndex].mInstanceTransforms.push_back(xform);

			numActive += asset->mNumParticles;
			//if instance has partner bind them
		}

		// update solver
		swap(g_buffers->positions, particles);
		swap(g_buffers->velocities, velocities);
		swap(g_buffers->phases, phases);
	}*/

    void randomDistribute(int N)
    {
        if (iBatches.size() == 0) return;
        Vec3 top = maxExtents;
        Vec3 bot = minExtents;
        std::vector<Vec3> halton_positions(N);
        Vec3 dim = top-bot;//Vec3(5, 5, 5);
        Vec3 center =  -bot;//dim/2.0f;
        float scale_dim[3] = { dim.x, dim.y, dim.z};//scale on x y z should be the bounding box size
        //int n = PoissonSample3D(2.45f, radius*0.42f, &positions[0], positions.size(), 2);
        int n = HaltonSample3D(scale_dim, main_radius, &halton_positions[0], halton_positions.size());

        //for each of this position place a solid object ?
        for ( int i=0;i<N;i++)
        {
            const int ingrIndex = Rand()%iBatches.size();
            Quat q = Quat(0,0,0,1);// QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);
            createFlexInstanceIngredient(ingrIndex, Vec3(halton_positions[i].x,halton_positions[i].y,halton_positions[i].z)-center, q);
            //createFlexInstanceIngredient(ingrIndex, Vec3(0,0,0), q);
        }
    }

	void UpdateInstanceTransform(int instId, Vec3 pos, Quat quat)
	{
		IngredientInstance& inst = mInstances[instId];
		FlexExtAsset* asset = iBatches[inst.mMeshIndex].mAsset;
		for (int j = 0; j < asset->mNumParticles; ++j)
		{
			Vec3 localPos = Vec3(&asset->mParticles[j * 4]);//local position of the proxy
			Vec3 rpos = Rotate(quat, localPos);
			g_buffers->positions[inst.mParticleOffset + j]=Vec4(pos + rpos, 1.0f);
		}
		g_rigidTranslations[instId].Set(pos.x, pos.y, pos.z);
		g_rigidRotations[instId].Set(quat.x, quat.y, quat.z, quat.w);
	}

	void checkPartnerAlongCurvePoints(Rope& rope, Json::Value ingr_node)
	{
		if (ingr_node["partners_name"].empty())
			return;
		int nPartner = ingr_node["partners_name"].size();
		int nPoints = rope.mIndices.size();
		cout << "ok " << nPoints << " " << nPartner << endl;
		std::vector<Json::Value> partner_nodes;
		std::vector<int> batchs_id;
		for (int j = 0; j < nPartner; j++)
		{
			Json::Value ingr_partner_node = getIngredients(ingr_node["partners_name"][j].asString());
			partner_nodes.push_back(ingr_partner_node);
			int batchid = getIngredientBatchId(ingr_node["partners_name"][j].asString());
			batchs_id.push_back(batchid);
		}
		int previously_use_point = 0;
		for (int curveI = 3; curveI < nPoints-3; curveI++)
		{
			//pick a partner
			int idPartner = Rand() % nPartner;
			Json::Value ingr_partner_node = partner_nodes[idPartner];// getIngredients(ingr_node["partners_name"][idPartner].asString());
			int batchid = batchs_id[idPartner];// getIngredientBatchId(ingr_node["partners_name"][idPartner].asString());
			//proba to bind
			int nbMol = ingr_partner_node["nbMol"].asInt();
			float proba = (float)nbMol / (float)nPoints;
			float r = Random(0.0f,1.0f);// Randf(0.0f, 1.0f);
			
			if (r > proba)
				continue;
			
			if (iBatches[batchid].nInstances >= nbMol)
				continue;

			//get position
			cout << "ok " << curveI<<" "<<proba << " " << r << " " << idPartner << " " << ingr_partner_node["nbMol"].asFloat() << endl;
			int curve_ln = ingr_partner_node["properties"]["range"][0].asInt();
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
			
			//cout << "before adding elem " << g_buffers->positions.size() << endl;

			//create instance
			int offset = createInstanceIngredient(batchid, pos, Quat(0,0,0,1));

			previously_use_point = middleIndex;
			if (!ingr_partner_node["properties"]["pairs"].empty())
			{
				cout << " pairs " << endl;
				for (int i = 0; i < ingr_partner_node["properties"]["pairs"].size(); i++){
					for (int j = 0; j < ingr_partner_node["properties"]["pairs"][i].size(); j++){
						int beads = ingr_partner_node["properties"]["pairs"][i][j].asInt();
						cout << "add1 " << curveI + i << " " << rope.mIndices[curveI + i] << " " << offset + beads << endl;
						if (curveI + i < nPoints)
							CreateSpringInter(rope.mIndices[curveI + i], offset + beads, 1.0f, 0.0f, main_radius*2.0f);
						//else 
						//	CreateSpring(rope.mIndices[(curveI + i) - rope.mIndices.size()], offset + beads, 1.0f, 0.0f, main_radius*2.0f);
					}
				}
			}
			else
			{
				for (int i = 0; i < ingr_partner_node["properties"]["beadsin"].size(); i++)
				{
					int beadsin = ingr_partner_node["properties"]["beadsin"][i].asInt();
					CreateSpringInter(rope.mIndices[curveI], offset + beadsin, 1, 0.0f, main_radius*2.0f);
					if (curveI + 1 < nPoints)
						CreateSpringInter(rope.mIndices[curveI + 1], offset + beadsin, 0.5f, 0.0f, main_radius*2.0f);
					if (curveI - 1 > 0 )
						CreateSpringInter(rope.mIndices[curveI - 1], offset + beadsin, 0.5f, 0.0f, main_radius*2.0f);
					cout << "addx " << curveI + 1 << " " << rope.mIndices[curveI] << " " << offset + beadsin << endl;
				}
				for (int i = 0; i < ingr_partner_node["properties"]["beadsout"].size(); i++)
				{
					int beadsout = ingr_partner_node["properties"]["beadsout"][i].asInt();
					int add_id = curveI + curve_ln + curve_var;
					if (add_id >=  nPoints)
					{
						add_id = add_id - nPoints;
					}
					if (add_id < nPoints)
					{
						CreateSpringInter(rope.mIndices[add_id], offset + beadsout, 1, 0.0f, main_radius*2.0f);
						cout << i << " " << beadsout << "addy " << add_id << " " << rope.mIndices[add_id] << " " << offset + beadsout << endl;
					}
					if (add_id + 1 < nPoints)
					{
						CreateSpringInter(rope.mIndices[add_id + 1], offset + beadsout, 0.5f, 0.0f, main_radius*2.0f);
						cout << i << " " << beadsout << "addy " << add_id+1 << " " << rope.mIndices[add_id+1] << " " << offset + beadsout << endl;
					}
					if (add_id - 1 > 0)
					{
						CreateSpringInter(rope.mIndices[add_id - 1], offset + beadsout, 0.5f, 0.0f, main_radius*2.0f);
						cout << i << " " << beadsout << "addy " << add_id -1<< " " << rope.mIndices[add_id-1] << " " << offset + beadsout << endl;
					}
					
				}
			}

		}
	
	}

	void placePartner()
	{
		printf("place partner %i\n", mInstances.size());
		for (int i = 0; i < mInstances.size(); i++) 
		{
			IngredientInstance& inst = mInstances[i];
			AssetBatch batch = iBatches[(int)inst.mMeshIndex];
			Json::Value ingr_node = getIngredients(pnames[(int)inst.mMeshIndex]);
			cout << "OK" << ingr_node["name"] << endl;
			Json::Value ingr_node_partners = ingr_node["partners_name"];
			//cout << "OK" << endl;
			cout << ingr_node_partners.empty() << " " << ingr_node_partners.isArray() << endl;
			if (!ingr_node_partners.empty()){
				for (int j = 0; j < ingr_node_partners.size(); j++)
				{
					cout << "place " << ingr_node_partners[j].asString() << endl;
					Json::Value ingr_partner = getIngredients(ingr_node_partners[j].asString());
					//AssetBatch aBatch = getIngredientBatch(ingr_node_partners[i].asString());
					//cout << "type is " << ingr_partner["Type"].asString() << endl;
					if (ingr_partner["Type"].asString() == "Grow"){
						//find a rope of this type
						Rope curve = g_ropes[0];
						placePartnerOnCurve(curve, ingr_node, i);
					}
				}
			}
		}
	}

	void placePartnerOnCurve(Rope& rope, Json::Value ingr, int instId)
	{
		//how many partner
		cout << "palce ingredient " << ingr["name"].asString() << "instance id "<< instId << endl;
		//get the ingredients property
		Json::Value ingrproperties = ingr["properties"];
		//"beadsin": 0,
		//	"beadsout" : 7,
		//	"range" : [70, 5]
		//use the folowing beads on the partner
		//int beadsin = ingrproperties["beadsin"].asInt();
		//int beadsout = ingrproperties["beadsout"].asInt();
		//distance between curve beads
		int curve_ln = ingrproperties["range"][0].asInt();
		int curve_var = int(Randf(-1.0f, 1.0f)* (float) ingrproperties["range"][1].asInt()); //Rand(-1, 1) % ingrproperties["range"][1].asInt();

		int batchid = getIngredientBatchId(ingr["name"].asString());
		//instance id mInstances
		
		IngredientInstance& inst = mInstances[instId];

		//pick a random position to attach
		int curve_index = Rand() % rope.mIndices.size();
		//what happen if I change the translation/rotation
		//g_rigidRotations[instId].Set
		//g_rigidTranslations[instId].Set

		//also insure curve_index not already use ?
		//FlexExtAsset* asset = iBatches[batchid].mAsset;
		//int npart = asset->mNumParticles;
		int offset = inst.mParticleOffset;
		int middleIndex = curve_index;
		//cout << curve_index << " " << curve_var << " " << beadsin << " " << beadsout << " offset " << offset << endl;
		//create spring between offset+beadsin  and curve_index
		//create spring between offset+beadsout  and curve_index+(curve_ln+/-curve_var)
		if (!ingr["properties"]["pairs"].empty())
		{
			cout << " paris " << endl;
			for (int i = 0; i < ingr["properties"]["pairs"].size(); i++){
				for (int j = 0; j < ingr["properties"]["pairs"][i].size(); j++){
					int beads = ingr["properties"]["pairs"][i][j].asInt();
					CreateSpringInter(rope.mIndices[curve_index + i], offset + beads, 1.0f, 0.0f, main_radius*2.0f);
				}
			}

			middleIndex = curve_index + ( ingr["properties"]["pairs"].size() / 2 ) - 1;
			
			if (middleIndex >= rope.mIndices.size())
			{
				middleIndex = middleIndex - rope.mIndices.size();
			}
			cout << middleIndex << " " << ingr["properties"]["pairs"].size() / 2 << " " << curve_index << endl;
		}
		else
		{
			for (int i = 0; i < ingrproperties["beadsin"].size(); i++)
			{
				int beadsin = ingrproperties["beadsin"][i].asInt();
				CreateSpringInter(rope.mIndices[curve_index + 1], offset + beadsin, 1, 0.0f, main_radius*2.0f);
				CreateSpringInter(rope.mIndices[curve_index], offset + beadsin, 1, 0.0f, main_radius*2.0f);
				CreateSpringInter(rope.mIndices[curve_index - 1], offset + beadsin, 1, 0.0f, main_radius*2.0f);
			}
			for (int i = 0; i < ingrproperties["beadsout"].size(); i++)
			{
				int beadsout = ingrproperties["beadsout"][i].asInt();
				CreateSpringInter(rope.mIndices[curve_index + curve_ln + curve_var + 1], offset + beadsout, 1, 0.0f, main_radius*2.0f);
				CreateSpringInter(rope.mIndices[curve_index + curve_ln + curve_var], offset + beadsout, 1, 0.0f, main_radius*2.0f);
				CreateSpringInter(rope.mIndices[curve_index + curve_ln + curve_var - 1], offset + beadsout, 1, 0.0f, main_radius*2.0f);
			}
		}

		Vec3 to1 = g_buffers->positions[rope.mIndices[middleIndex + 1]] - g_buffers->positions[rope.mIndices[middleIndex]];
		Vec3 to2 = g_buffers->positions[rope.mIndices[middleIndex + 2]] - g_buffers->positions[rope.mIndices[middleIndex + 1]];
		Vec3 up = Cross(Normalize(to1), Vec3(0, 0, 1));
		cout << "palce ingredient " << up.x << " " <<
			up.y << " " <<
			up.z << endl;

		Vec3 ingrpcpal = Vec3(ingr["principalVector"][0].asFloat()*1.0f,
			ingr["principalVector"][1].asFloat()*1.0f,
			ingr["principalVector"][2].asFloat())*1.0f;
		cout << "palce ingredient " << ingrpcpal.x << " " << ingrpcpal.y << " " << ingrpcpal.z << endl;
		Vec3 offsetPos = Vec3(ingr["offset"][0].asFloat()*main_scale,
			ingr["offset"][1].asFloat()*main_scale,
			ingr["offset"][2].asFloat())*main_scale;
		//align pcpal to up
		Quat rot = AlignVec3s(up, ingrpcpal);
		cout << "palce ingredient " << rot.x << " " << rot.y << " " << rot.z << endl;
		Vec3 pos = g_buffers->positions[rope.mIndices[middleIndex]];// +to1*1.0f / 2.0f) + offset;
		pos += to1 / 2.0f;
		//pos += Rotate(rot, offsetPos);
		//update rigidposition of instance
		UpdateInstanceTransform(instId, pos, rot);

		/*
		CreateSpring(rope.mIndices[curve_index], offset + beadsin, 1, 0.0f, main_radius*2.0f);
		CreateSpring(rope.mIndices[curve_index + curve_ln + curve_var], offset + beadsout,1, 0.0f, main_radius*2.0f);

		CreateSpring(rope.mIndices[curve_index-1], offset + beadsin, 1, 0.0f, main_radius*2.0f);
		CreateSpring(rope.mIndices[curve_index + curve_ln + curve_var-1], offset + beadsout, 1, 0.0f, main_radius*2.0f);


		CreateSpring(rope.mIndices[curve_index], offset + beadsin +1 , 1, 0.0f, main_radius*4.0f);
		CreateSpring(rope.mIndices[curve_index + curve_ln + curve_var], offset + beadsout -1 , 1, 0.0f, main_radius*4.0f);

		CreateSpring(rope.mIndices[curve_index - 1], offset + beadsin+1, 1, 0.0f, main_radius*4.0f);
		CreateSpring(rope.mIndices[curve_index + curve_ln + curve_var - 1], offset + beadsout - 1, 1, 0.0f, main_radius*4.0f);
		*/
	}

	void placeOneFiberSurface(int fiber_id){
		Json::Value ingr_node_name = getIngredients(pnames_fiber[fiber_id]);
		int compId = iBatchesFiber[fiber_id].compId;
		cout << "build rope " << endl;
		int ncomp = comp_mesh.size();
		int nVertices = comp_mesh[compId - 1]->GetNumVertices();
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

	void placeOneFiber(int fiber_id, bool partners=false)
	{
		Json::Value ingr_node_name = getIngredients(pnames_fiber[fiber_id]);
		int compId = iBatchesFiber[fiber_id].compId;
		cout << "build rope " << fiber_id << endl;
		int nbCurve = ingr_node_name["nbMol"].asInt();
		float L = ingr_node_name["length"].asFloat()*main_scale;//in angstrom
		int npoints = (int)(L) / (main_radius*2.0f);
		float *data_curve;
		//starting configuration is either a circle, a line or a given file
		string mode = ingr_node_name["startingMode"].asString();
		bool close = ingr_node_name["closed"].asBool();//in angstrom
		float D = g_params.mRadius;//(g_params.mRadius/4.0f)*1000.0f;//could be the persitence length, which is 450bp 34A is one turn 10bp. Distance between spheres
		int subdivid = 1;// (int)((L / D) / (float)N);
		int persistence = 3;// (int)round(((500.0f*main_scale) / (main_radius*2.0f)) / 4.0f);
		bool extend = false;
		if (mode == "circle")
		{
			float radius_circle = (L * 50.0f) / k2Pi;
			nfloat = npoints * 3;
			data_curve = new float[nfloat];
			int count = 0;
			//a go from 0 to k2Pi
			// circle 1 238 8 37.8789 1 0.1185
			cout << "circle " << nbCurve << " " << L << " " << npoints << " " << radius_circle << " " << close << " " << main_radius << endl;
			float a = 0.0f;
			for (int i = 0; i < npoints; i++){
				data_curve[count] = radius_circle * cos(a);
				data_curve[count + 1] = radius_circle * sin(a);
				data_curve[count + 2] = 0.0f;
				count += 3;
				a += k2Pi / (float)npoints;
			}
			extend = false;
		}
		else if (mode == "file")
		{
			//use thebinary file...dna
			ifstream ifs2(datapath + "tps_path.bin", ios::binary | ios::ate);
			ifstream::pos_type pos2 = ifs2.tellg();

			result_curve.resize(pos2);

			ifs2.seekg(0, ios::beg);
			ifs2.read(&result_curve[0], pos2);

			printf("read bytes %i\n", result_curve.size());
			//now convert to float and place a particle there
			int step = 3;

			nfloat = result_curve.size() / sizeof(float);
			printf("which are %i float %i %i\n", nfloat, nfloat / step, result_curve.size() / sizeof(float));

			int N = nfloat / step;
			float PL = 500.0f;//thats 25 spheres
			subdivid = (int)((L / D) / (float)N);
			data_curve = reinterpret_cast<float*>(result_curve.data());
			//g_numExtraParticles = int(L/D);
			printf(" use N=%i L=%f D=%f total excpected=%i\n", N, L, D, npoints);
			printf("create Rope from data with sibdivd %i and D %f and %i nExtraParticle\n", subdivid, D, g_numExtraParticles);
			extend = true;
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
		rope_phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);

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
		printf("create Rope from data OK with %i points\n", curve.mIndices.size());

		// liste of couple Fixed
		/*std::vector<int> couple = { 1, 50, 11, 40, //21, 30,
		71, 120, 81, 110, //91, 100,
		141, 190, 151, 180, //161, 170,
		211, 260, 221, 250, //231, 240,
		281, 330, 291, 320, //301, 310,
		60, 130, 200, 340
		};*/
		// liste of random couple for loop making
		//npoints
		/*
		int nloop = 1;
		count = 0;
		std::vector<int> couple = { 0, 0 };
		for (int i = 0; i < nloop; i++)
		{
		couple[count] = int(Randf(0.0f, 0.5f)*npoints);
		couple[count + 1] = int(Randf(0.5f, 1.00f) * npoints);
		//couple[count + 2] = couple[count] + 10;
		//couple[count + 3] = couple[count + 1] + 10;
		count += 2;
		}
		count = 0;
		for (int i = 0; i < couple.size() / 2; i++){
		//CreateSpring(r_dna.mIndices[couple[count]], r_dna.mIndices[couple[count+1]], 0.01f, 0.0f, main_radius*2.0f);
		count += 2;
		}
		*/
		g_ropes.push_back(curve);//instance
		mask.push_back(-1);
		maks_protein.push_back(-(fiber_id + 1));//it ptype is 0 ?
		maks_fiber.push_back(fiber_id);//proteinType
		printf("create Rope from data OK with %i  %i\n", -(fiber_id + 1), fiber_id);
		if (partners) 
		{
			checkPartnerAlongCurvePoints(curve, ingr_node_name);
			//CompactObjects();
		}
	}

	void placeFibers(bool check_partners){
        //loop over fiber ingredients
        for (int i=0;i<pnames_fiber.size();i++)
        {

			cout << "place fibers " << pnames_fiber[i] << " " << i << endl;
			if (pnames_fiber[i] != "DNA") continue;
            Json::Value ingr_node_name = getIngredients(pnames_fiber[i]);
			//get partner if any
			IngredientPartner partners;
			partners.nPartner = ingr_node_name["partners_name"].size();
			if (partners.nPartner > 0)
			{
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
				placeOneFiber(i, check_partners);
				break;
			}

        } 
    }
    
    void growOneCurveFromSurface(int pType, int compId, int vIndex, float uLength, float length){
        cout << "one curve "  << compId-1 << " " << comp_mesh.size() << endl; 
        Vec3 start = Vec3(comp_mesh[compId-1]->m_positions[vIndex].x,comp_mesh[compId-1]->m_positions[vIndex].y,comp_mesh[compId-1]->m_positions[vIndex].z);
        Vec3 normal = Vec3(comp_mesh[compId-1]->m_normals[vIndex].x,comp_mesh[compId-1]->m_normals[vIndex].y,comp_mesh[compId-1]->m_normals[vIndex].z);
        std::vector<Vec3> startend;
        startend.push_back(start);
        startend.push_back(start+(normal/Length(normal))*length);
        int curve_phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);
        
        float* datac = reinterpret_cast<float*>(startend.data());
        int nbfloat = 6;
        Rope curve;
        int subdivid = (int) length/(main_radius*2.0);
        
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
        //we still need to fix start to surface, we could attach the particle using sdf, or mass ?
        // put huge infiinte ? mass
        g_buffers->positions[curve.mIndices[0]][3] = 0.0f; //first points attached or both first points
        g_ropes.push_back(curve);//instance

		mask.push_back(-1);
        maks_protein.push_back(-(pType+1));//it ptype is 0 ?
        maks_fiber.push_back(pType);//proteinType
    }
    
	void loadResults( string filename, bool redo=false )
	{
		std::ifstream ifs(filename);
		if (ifs.is_open())
		{
			ifs >> results_json;
			std::cout << "file opened and closed";
		}
		else
		{
			std::cout << "Error opening file "+filename;
		}
		
		Json::Value  cyto = results_json["cytoplasme"];
		if (cyto != 0)
		{
			std::cout << "find compartments cyto " << cyto.size() << endl;
			Json::Value cyto_ingredients = cyto["ingredients"];
			std::cout << "compartment cyto should have n ingredients " << cyto_ingredients.size() << endl;
			if (cyto_ingredients != 0)
			{
				//parseJsonIngredientsResults(cyto_ingredients, 0,redo);
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
				std::cout << "compartment should have several childs " << comp_name.size() << endl;
				if (comp_name.size() == 0) continue;
				Json::Value comp_surface = comp_name["surface"];//at 0 ?
				if (comp_surface != 0)
				{
					Json::Value surf_ingredients = comp_surface["ingredients"];
					std::cout << "compartment should have n ingredients " << surf_ingredients.size() << endl;
					if (surf_ingredients != 0)
					{
						parseJsonIngredientsResults(surf_ingredients, i + 1,redo);
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
						//parseJsonIngredientsResults(int_ingredients, -(i + 1),redo);
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

    void haltondistribute(int N)
    {
        //PyGILState_STATE gstate;
        //gstate = PyGILState_Ensure();
    
        //PyObject *pName, *pModule, *pDict, *pFunc;
        //PyObject *pArgs, *pValue;
        
        //PyObject *m = PyImport_AddModule("__main__");
        
        //PyObject *main_module = PyImport_ImportModule("__main__");
        //PyObject *main_dict   = PyModule_GetDict(main_module);
        
        //pFunc = PyMapping_GetItemString(main_dict, "pickIngredient");
        
        //pFunc = PyObject_GetAttrString(m,"pickIngredient");
        
        if (iBatches.size() == 0) return;
		        
		Vec3 top = maxExtents;
        Vec3 bot = minExtents;
        std::vector<Vec3> halton_positions(N);
        Vec3 dim = (top-bot)*0.75f;//Vec3(5, 5, 5);
        Vec3 center =  (-bot)*0.75f;//dim/2.0f;
        float scale_dim[3] = { dim.x, dim.y, dim.z};//scale on x y z should be the bounding box size
        //int n = PoissonSample3D(2.45f, radius*0.42f, &positions[0], positions.size(), 2);
        int n = HaltonSample3D(scale_dim, main_radius, &halton_positions[0], halton_positions.size());

        //for each of this position place a solid object ?
        int particleOffset = g_buffers->positions.size();
        //comp mask
        //compmask = new CompMask();
		//while not complete
		int count = 0;
		int safety_count = 0;
		int hpos = 0;
		while (count < N){
        //for ( int i=0;i<N;i++){
            // choose a random mesh to emit
            //check where is the points
            //
            //pValue = PyObject_CallFunction(pFunc, NULL);
            //PyObject* objectsRepresentation = PyObject_Repr(pValue);
            //const char* name = PyString_AsString(objectsRepresentation);
            //Py_DECREF(objectsRepresentation);
            //cout << "ingredient name is : " << name << endl;
            //getIngredients(string name)
            //int ingrIndex = getIngredientBatchId(name);

			//check if i inside the compartments?

			int ingrIndex =  Rand() % iBatches.size();
			
			int compId = iBatches[ingrIndex].compId;
			Vec3 p = Vec3(halton_positions[hpos].x, halton_positions[hpos].y, halton_positions[hpos].z) - center;
			Quat r = QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);
			SDF * fsdf = comp_shape[abs(compId) - 1];
			float spacing =  (fsdf->mUpper[0] - fsdf->mLower[0]) / (float)fsdf->mWidth;
			Vec3 botsdf = Vec3(fsdf->mLower[0], fsdf->mLower[1], fsdf->mLower[2]);
			Vec3 topsdf = Vec3(fsdf->mUpper[0], fsdf->mUpper[1], fsdf->mUpper[2]);
			float3 pquery = { (p.x - botsdf.x) / spacing, (p.y - botsdf.y) / spacing, (p.z - botsdf.z) / spacing };
			int3 pquery_index = { (int)pquery.x, (int)pquery.y, (int) pquery.z };
			int dimFF = fsdf->mDepth;
			float D = SampleSDFX(fsdf->mField, dimFF, (int)pquery.x, (int)pquery.y, (int)pquery.z);
			//cout << abs(compId) - 1 << " dimf " << dimFF << " index " << (int)pquery.x << " " << (int)pquery.y << " " << (int)pquery.z << " " << count << " " << D << endl;
			
			if (compId < 0) //inside
			{
				//continue;
				if (D > 0.0f)
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
				int ncomp = comp_mesh.size();
				int nVertices = comp_mesh[compId - 1]->GetNumVertices();
				int vIndex = Rand() % nVertices;
				p = Vec3(comp_mesh[compId - 1]->m_positions[vIndex].x, comp_mesh[compId - 1]->m_positions[vIndex].y, comp_mesh[compId - 1]->m_positions[vIndex].z);
				Vec3 normal = Vec3(comp_mesh[compId - 1]->m_normals[vIndex].x, comp_mesh[compId - 1]->m_normals[vIndex].y, comp_mesh[compId - 1]->m_normals[vIndex].z);
				//align pcpalVector to normal
				Vec3 ingrpcpal = Vec3(iBatches[ingrIndex].pcpalVectorx, iBatches[ingrIndex].pcpalVectory, iBatches[ingrIndex].pcpalVectorz);
				Vec3 offsetPos = Vec3(iBatches[ingrIndex].offsetx, iBatches[ingrIndex].offsety, iBatches[ingrIndex].offsetz);
				//align pcpal to up
				r = AlignVec3s(normal, ingrpcpal);
				p = p + r*offsetPos;
			}

			FlexExtAsset* asset = iBatches[ingrIndex].mAsset;
            // check we can fit in the container
            //if (int(g_buffers->positions.size()) - particleOffset < asset->mNumParticles)
            //    break;

            IngredientInstance inst;
            inst.mLifetime = 1000.0f;// Randf(5.0f, 60.0f);
            inst.mParticleOffset = particleOffset;
			inst.mRotation = r;// QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);
			
			inst.mTranslation = p;//g_emitters[0].mPos + Vec3(Randf(-spread, spread), Randf(-spread, spread), 0.0f);
            //inst.mTranslation = Vec3(0.0f);//g_emitters[0].mPos + Vec3(Randf(-spread, spread), Randf(-spread, spread), 0.0f);

            inst.mMeshIndex = ingrIndex;

            Vec3 linearVelocity = Vec3(0.0f);//g_emitters[0].mDir*15.0f;//*Randf(5.0f, 10.0f);//Vec3(Randf(0.0f, 10.0f), 0.0f, 0.0f);
            Vec3 angularVelocity = Vec3(0.0f);//Vec3(UniformSampleSphere()*Randf()*k2Pi);

            mask.push_back(iBatches[ingrIndex].compId);
            maks_protein.push_back(ingrIndex);
            inst.mGroup = iGroupCounter++;

            int mm = 2^24;
            const int phase = NvFlexMakePhase(inst.mGroup, 0);
            int g = TAKE_N_BITS_FROM(phase,0,24);
            //cout << "phase is " << phase << " group " << inst.mGroup << " g " << g << " index "<<ingrIndex << " comp " << iBatches[ingrIndex].compId << " size " << iBatches.size() << endl;
            // generate initial particle positions
            for (int j=0; j < asset->mNumParticles; ++j)
            {
                Vec3 localPos = Vec3(&asset->mParticles[j*4]);// - Vec3(&asset->mShapeCenters[0]);

				g_buffers->positions.push_back(Vec4(inst.mTranslation + localPos, 1.0f));//inst.mRotation*
                g_buffers->velocities.push_back(Vec3(0.0f));//linearVelocity + Cross(angularVelocity, localPos);
                g_buffers->phases.push_back(phase);
            }

            particleOffset += asset->mNumParticles;

            mInstances.push_back(inst);
            iBatches[ingrIndex].nInstances++;
			count++;
			hpos++;
			//iBatches[ingrIndex].nInstances
        }
        
        //PyGILState_Release(gstate);
        // compact instances
        static std::vector<Vec4> particles(g_buffers->positions.size());
        static std::vector<Vec3> velocities(g_buffers->velocities.size());
        static std::vector<int> phases(g_buffers->phases.size());

        g_rigidTranslations.resize(0);
        g_rigidRotations.resize(0);
        g_rigidCoefficients.resize(0);
        g_rigidIndices.resize(0);
        g_rigidLocalPositions.resize(0);
        g_rigidOffsets.resize(0);

        // start index
        g_rigidOffsets.push_back(0);

        // clear mesh batches
        //for (int i=0; i < int(iBatches.size()); ++i)
        //    iBatches[i].mInstanceTransforms.resize(0);

        int numActive = 0;

        for (int i=0; i < int(mInstances.size()); ++i)
        {
            IngredientInstance& inst = mInstances[i];

            FlexExtAsset* asset = iBatches[inst.mMeshIndex].mAsset;

            for (int j=0; j < asset->mNumParticles; ++j)
            {
                particles[numActive+j] = g_buffers->positions[inst.mParticleOffset+j];
                //velocities[numActive+j] = g_buffers->velocities[inst.mParticleOffset+j];
                phases[numActive+j] = g_buffers->phases[inst.mParticleOffset+j];
            }

            g_rigidCoefficients.push_back(1.0f);
            g_rigidTranslations.push_back(inst.mTranslation);
			g_rigidRotations.push_back(inst.mRotation);

            for (int j=0; j < asset->mNumShapeIndices; ++j)
            {
                g_rigidLocalPositions.push_back(Vec3(&asset->mParticles[j*4]) - Vec3(&asset->mShapeCenters[0]));
                g_rigidIndices.push_back(asset->mShapeIndices[j] + numActive);
            }

            g_rigidOffsets.push_back(g_rigidIndices.size());

            mInstances[i].mParticleOffset = numActive;

            // Draw transform
            //Matrix44 xform = TranslationMatrix(Point3(inst.mTranslation));// - inst.mRotation*Vec3(asset->mShapeCenters)))*RotationMatrix(inst.mRotation);
            //iBatches[inst.mMeshIndex].mInstanceTransforms.push_back(xform);

            numActive += asset->mNumParticles;
			//if instance has partner bind them
        }

        // update solver
        swap(g_buffers->positions, particles);
        swap(g_buffers->velocities, velocities);
        swap(g_buffers->phases, phases);

    }

    void initProceduralRope(){
        //initialized if not
        //update frequency
        //update nb
        //every n frame add m particle to rope r
        int step = 7;
        rope_unit_length = main_radius;
        rope_phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);
        ifstream ifs2(datapath+"tps_path.bin", ios::binary|ios::ate);
        ifstream::pos_type pos2 = ifs2.tellg();

        result_curve.resize(pos2);

        ifs2.seekg(0, ios::beg);
        ifs2.read(&result_curve[0], pos2);

        printf("read bytes %i\n",result_curve.size());
        //now convert to float and place a particle there
        step = 3;

        nfloat = result_curve.size()/sizeof(float);
        printf("which are %i float %i %i\n",nfloat,nfloat/step,result_curve.size()/sizeof(float));

        int N = nfloat/step;
        float L = 4119790.0f*main_scale;
        //20? 0.08f
        //distance between point. should be the spring distance
        //different from the radius of interaction and the rest radius , radius should be less ?
        //if radius is 34
        //D is the length of the spring between particle for the rope
        float D = rope_unit_length*2.0f;//(g_params.mRadius/4.0f)*1000.0f;//could be the persitence length, which is 450bp 34A is one turn 10bp. Distance between spheres
        //D = D + D/10.0f;
        float PL = 500.0f;//thats 25 spheres
        //D=4*D;
		int subdivid =  (int)((L / D) / (float)N);
        data_curve = reinterpret_cast<float*>(result_curve.data());
        //g_numExtraParticles = int(L/D);
        printf (" use N=%i L=%f D=%f g_numExtraParticles=%i\n",N,L,D,g_numExtraParticles);
        printf("create Rope from data with sibdivd %i and D %f and %i nExtraParticle\n",subdivid,D,g_numExtraParticles);

        //CreateRopeFromData(r, 1.0f, data_curve, 110.0f, nfloat, rope_phase,0.0f,1.0f,0.0f,51,true,true,(D/2.0f)/1000.0f);
        //CreateDualFromData(r, 1.0f, data_curve, 110.0f, nfloat, rope_phase,0.0f,1.0f,0.0f,subdivid,true,true,D/1000.0f);
        //CreateTriangleFromData(r, 1.0f, data_curve, 110.0f, nfloat, rope_phase,0.0f,1.0f,0.0f,subdivid,true,true,(D/2.0f)*main_scale);
		//nfloat = 30 * 4;
        CreateRopeFromData(r_dna, //rope
                           main_scale, // scale
                           1.0f, //stifness
                           data_curve, //data
                           110.0f*2.0f, //length
                           nfloat,  //nfloat
                           rope_phase,//phase
                           0.0f,//spiral angle
                           1.0f,//invmass
                           0.0f,//give
                           subdivid,//extend_nb
                           false,//extend
                           true,//closew
						   main_radius*2.0f);

        printf("create Rope from data OK with %i points\n",r_dna.mIndices.size());
		
		int index1 = Rand() % r_dna.mIndices.size();
		int index2 = Rand() % r_dna.mIndices.size();

		//CreateSpring(r_dna.mIndices[index1], r_dna.mIndices[index2], 1.0f, 0.0f, main_radius*2.0f);
        g_ropes.push_back(r_dna);//instance

        mask.push_back(-1);
        maks_protein.push_back(-2);
        maks_fiber.push_back(0);//proteinType

    }

	void initDNACircleExperiment(){
		//initialized if not
		//update frequency
		//update nb
		//every n frame add m particle to rope rs
		int nb = 100;
		//or 100 ?
		//int npoints = 350;//should be 55 beads...
		//generate coordinate circle
		float L = (float) nb * 3.4f;
		float r = L;//Length in Angstrom circunmference ?
		int npoints = (int)((L*main_scale) / (g_params.mRadius/2.0f));
		nfloat = npoints * 3;
		float* data_circle = new float[nfloat];
		int count = 0;
		//a go from 0 to k2Pi
		
		float a = 0.0f;
		for (int i = 0; i < npoints; i++){
			data_circle[count] = r * cos(a);
			data_circle[count + 1] = r * sin(a);
			data_circle[count + 2] = 0.0f;
			count += 3;
			a += k2Pi / (float)npoints;
		}

		int step = 7;
		rope_unit_length = main_radius;
		rope_phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);
		
		float D = rope_unit_length*2.0f;//(g_params.mRadius/4.0f)*1000.0f;//could be the persitence length, which is 450bp 34A is one turn 10bp. Distance between spheres
		int subdivid = 1;// (int)((L / D) / (float)N);
		int persistence = 2;// (int)round(((500.0f*main_scale) / (main_radius*2.0f)) / 4.0f);
		CreateRopeFromData(r_dna, //rope
			main_scale, // scale
			1.0f, //stifness
			data_circle, //data
			110.0f*2.0f, //length
			nfloat,  //nfloat
			rope_phase,//phase
			0.0f,//spiral angle
			1.0f,//invmass
			0.0f,//give
			subdivid,//extend_nb
			false,//extend
			true,//close
			main_radius*2.0f,
			persistence);
		r_dna.persistence = persistence;
		printf("create Rope from data OK with %i points\n", r_dna.mIndices.size());
		// liste of couple Fixed
		/*std::vector<int> couple = { 1, 50, 11, 40, //21, 30, 
			71, 120, 81, 110, //91, 100, 
			141, 190, 151, 180, //161, 170,
			211, 260, 221, 250, //231, 240, 
			281, 330, 291, 320, //301, 310, 
			60, 130, 200, 340
		};*/
		// liste of random couple for loop making
		//npoints
		int nloop = 1;
		count = 0;
		std::vector<int> couple = { 0, 0 };
		for (int i = 0; i < nloop; i++)
		{
			couple[count] = int(Randf(0.0f, 0.5f)*npoints);
			couple[count + 1] = int(Randf(0.5f, 1.00f) * npoints);
			//couple[count + 2] = couple[count] + 10;
			//couple[count + 3] = couple[count + 1] + 10;
			count += 2;
		}
		count = 0;
		for (int i = 0; i < couple.size()/2; i++){
			//CreateSpring(r_dna.mIndices[couple[count]], r_dna.mIndices[couple[count+1]], 0.01f, 0.0f, main_radius*2.0f);
			count += 2;
		}
		
		g_ropes.push_back(r_dna);//instance

		mask.push_back(-1);
		maks_protein.push_back(-2);
		maks_fiber.push_back(0);//proteinType

	}

    virtual void insertPointRope(int rope_id,int ipair,float D, float give, 
								float stiffness, int aphase, int current)
	{
        Rope& arope = g_ropes[rope_id];
        //std::cout << current << " insertPoint at : " << ipair << " particle id in rope " << arope.mIndices[ipair] << " " << D << endl;
        if (arope.mIndices[ipair] == 0) return;
        Vec4 point = Vec4(g_buffers->positions[arope.mIndices[ipair]]);
        Vec4 next_point = Vec4(g_buffers->positions[arope.mIndices[ipair+1]]);
        Vec3 dirtopt = Vec3 (next_point.x - point.x,next_point.y - point.y,next_point.z - point.z);

        g_buffers->positions[current]=Vec4(  point.x+dirtopt.x*0.5f,
                                    point.y+dirtopt.y*0.5f,
                                    point.z+dirtopt.z*0.5f,
                                    1.0f);
        g_buffers->velocities[current]=0.0f;
        g_buffers->phases[current]=aphase;

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

    int proceduralRope(int rope_id,int Nsub){
		float D = g_params.mRadius;// main_radius*2.0f;//g_params.mRadius;//*1000.0f;//could be the persitence length, which is 450bp 34A is one turn 10bp. Distance between spheres
        //D = D + D/10.0f;
        float give =0.0f;
        float stiffness = 1.0f;
        //initialized if not
        //update frequency
        //update nb
        //every n frame add m particle to rope r
        int current = flexGetActiveCount(g_flex);
		int ipair = int(Randf(0.75f, 0.90f)*g_ropes[0].mIndices.size());
		//std::cout << current << " insert at : " << ipair << " particle id in rope " << g_ropes[0].mIndices[ipair] << " " << D << " " << rope_phase << " " << g_ropes[0].mIndices.size() << endl;
        for (int i=0;i<Nsub;i++){
			insertPointRope(rope_id, ipair, D, give, stiffness, rope_phase, current);
            current++;
        }
        //insert Hue ?

        float h = Randf(0.0f,1.0f);
        //if (h < 0.045390661f)
        //    makeHue(ipair);//need to record the position of hue
        //0.045390661

        flexSetSprings(g_flex, &g_buffers->springIndices[0], &g_buffers->springLengths[0], &g_buffers->springStiffness[0], g_buffers->springLengths.size(), eFlexMemoryHost);
        flexSetParticles(g_flex, &g_buffers->positions[0].x, g_buffers->positions.size(), eFlexMemoryHost);
        flexSetVelocities(g_flex, &g_buffers->velocities[0].x, g_buffers->velocities.size(), eFlexMemoryHost);
        flexSetPhases(g_flex, &g_buffers->phases[0], g_buffers->phases.size(), eFlexMemoryHost);
		g_activeIndices.resize(current);// g_buffers->positions.size());
		for (size_t i = 0; i < current;i++)// g_activeIndices.size(); ++i)
            g_activeIndices[i] = i;
        flexSetActive(g_flex, &g_activeIndices[0], g_activeIndices.size(), eFlexMemoryHost);
		return ipair;
		/*
        Vec4 point = Vec4(g_buffers->positions[r_dna.mIndices[ipair]]);
        Vec4 next_point = Vec4(g_buffers->positions[r_dna.mIndices[ipair+1]]);
        Vec3 dirtopt = Vec3 (next_point.x - point.x,next_point.y - point.y,next_point.z - point.z);
        //int current = int(g_buffers->positions.size());


        g_buffers->positions[current]=Vec4(  point.x+dirtopt.x*0.5f,
                                    point.y+dirtopt.y*0.5f,
                                    point.z+dirtopt.z*0.5f,
                                    1.0f);
        g_buffers->velocities[current]=0.0f;
        g_buffers->phases[current]=rope_phase;//int(g_buffers->positions.size()));
        //replace 3 spring and create one
        //printf ("subdivide how many %i %i %i %i\n",i,j,extend_nb,current);
        int p_id = r_dna.mIndices[ipair];//before insertion
        ReplaceSpring(p_id,p_id+1,p_id,current);
        CreateSpring(current, p_id+1 ,stiffness, give,D);//current +1 doesnt exist yet
        //1-3
        ReplaceSpring(p_id-1,p_id+1,p_id-1,current);//replace
        ReplaceSpring(p_id,p_id+2,p_id,p_id+1);
        float r=Randf(-1.0f,1.0f)/1000.0f;
        CreateSpring(current, p_id+2, stiffness*0.5f, give,(D*2.0f)+r);
        //CreateSpring(ipair, current ,stiffness, give,D);//current +1 doesnt exist yet

        //printf ("bias is %f\n",r);
        //CreateSpring(ipair-1, current, stiffness*0.5f, give,(D*2.0f)+r);
        //CreateSpring(current, ipair+2, stiffness*0.5f, give,(D*2.0f)+r);
        //we need to remove the spring for
        std::vector<int>::iterator it = r_dna.coarseIndices.begin();
        r_dna.coarseIndices.insert (it+ipair+1,current);
        it = r_dna.mIndices.begin();
        //std::cout << "indices  " << current <<" " <<r_dna.mIndices[ipair-1] << " " <<r_dna.mIndices[ipair] <<" " <<r_dna.mIndices[ipair+1] <<" "<< r_dna.mIndices[ipair+2] << endl;
        r_dna.mIndices.insert (it+ipair+1,current);
        //std::cout << "indices  " <<r_dna.mIndices[ipair-1] << " " <<r_dna.mIndices[ipair] <<" " <<r_dna.mIndices[ipair+1] <<" "<< r_dna.mIndices[ipair+2] << endl;

        flexSetSprings(g_flex, &g_buffers->springIndices[0], &g_buffers->springLengths[0], &g_buffers->springStiffness[0], g_buffers->springLengths.size(), eFlexMemoryHost);
        flexSetParticles(g_flex, &g_buffers->positions[0].x, g_buffers->positions.size(), eFlexMemoryHost);
        flexSetVelocities(g_flex, &g_buffers->velocities[0].x, g_buffers->velocities.size(), eFlexMemoryHost);
        flexSetPhases(g_flex, &g_buffers->phases[0], g_buffers->phases.size(), eFlexMemoryHost);
        flexSetActive(g_flex, &g_activeIndices[0], current+1, eFlexMemoryHost);
		*/
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
		return L1/main_scale;
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
				cout << "place fibers " << pnames_fiber[fiber_id] << " rope id " << i << " compId " << compId << endl;
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

	int getInstanceOfType(int batchId, std::vector<int> used_instances_id){
		int instance_id = -1;
		for (int i = 0; i < mInstances.size(); i++)
		{
			if (std::find(used_instances_id.begin(), used_instances_id.end(), i) != used_instances_id.end())
				continue;
			IngredientInstance& inst = mInstances[i];

			if ((int)inst.mMeshIndex == batchId)
			{
				instance_id = i;
				break;
			}
		}
		return instance_id;
	}

	void pickPlacePartner(int fiber_id, int rope_id, int curveI, Json::Value ingr_node_name)
	{
		Rope rope = g_ropes[rope_id];
		int current_nPoints = rope.mIndices.size();

		float L = ingr_node_name["length"].asFloat()*main_scale;//in angstrom
		int target_nPoints = (int)(L*main_scale) / (main_radius*2.0f);

		int idPartner = Rand() % iPartnersFibers[fiber_id].nPartner;

		cout << "OK pick " << idPartner << " " << iPartnersFibers[fiber_id].partner_nodes.size() << endl;
		cout << "name" << iPartnersFibers[fiber_id].partner_nodes[idPartner]["name"] << endl;
		Json::Value ingr_partner_node = iPartnersFibers[fiber_id].partner_nodes[idPartner];// getIngredients(ingr_node["partners_name"][idPartner].asString());
		int batchid = iPartnersFibers[fiber_id].batchs_id[idPartner];// getIngredientBatchId(ingr_node["partners_name"][idPartner].asString());
		//proba to bind
		int nbMol = ingr_partner_node["nbMol"].asInt();
		float proba = (float)nbMol / (float)target_nPoints;
		float r = Random(0.0f, 1.0f);// Randf(0.0f, 1.0f);
		cout << "proba " << nbMol << " " << target_nPoints << " " << proba << " " << r << " " << (r > proba) << endl;
		if (r > proba)
			return;

		//if (iBatches[batchid].nInstances >= nbMol)
		//	return;

		//get position
		cout << "ok " << curveI << " " << proba << " " << r << " " << idPartner << " " << ingr_partner_node["nbMol"].asFloat() << endl;
		int curve_ln = ingr_partner_node["properties"]["range"][0].asInt();
		if (curve_ln > current_nPoints)
		{
			curve_ln = current_nPoints / 2;
		}
		int curve_var = int(Randf(-1.0f, 1.0f)* (float)ingr_partner_node["properties"]["range"][1].asInt()); //Rand(-1, 1) % ingrproperties["range"][1].asInt();
		if (curveI + curve_ln + curve_var >= current_nPoints)
		{
			curve_ln = (curveI + curve_ln + curve_var) - current_nPoints;
		}
		int middleIndex = curveI;
		if (!ingr_partner_node["properties"]["pairs"].empty())
		{
			middleIndex = curveI + (ingr_partner_node["properties"]["pairs"].size() / 2) - 1;
		}

		//if (middleIndex - previously_use_point < 3)
		//	//make sure not two ingredient place to close
		//	continue;

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

		//cout << "before adding elem " << g_buffers->positions.size() << endl;

		//get instance
		//maks_protein
		//instances
		int instance_id = getInstanceOfType(batchid, iPartnersFibers[fiber_id].instances_id);
		if (instance_id == -1)
			return;
		iPartnersFibers[fiber_id].instances_id.push_back(instance_id);
		//reposition the instances
		UpdateInstanceTransform(instance_id, pos, rot);
		int offset = mInstances[instance_id].mParticleOffset;
		//int offset = createInstanceIngredient(batchid, pos, Quat(0, 0, 0, 1));

		//previously_use_point = middleIndex;
		if (!ingr_partner_node["properties"]["pairs"].empty())
		{
			cout << " pairs " << endl;
			for (int i = 0; i < ingr_partner_node["properties"]["pairs"].size(); i++){
				for (int j = 0; j < ingr_partner_node["properties"]["pairs"][i].size(); j++){
					int beads = ingr_partner_node["properties"]["pairs"][i][j].asInt();
					cout << "add1 " << curveI + i << " " << rope.mIndices[curveI + i] << " " << offset + beads << endl;
					if (curveI + i < current_nPoints)
						CreateSpringInter(rope.mIndices[curveI + i], offset + beads, 1.0f, 0.0f, g_params.mRadius);
					//else 
					//	CreateSpring(rope.mIndices[(curveI + i) - rope.mIndices.size()], offset + beads, 1.0f, 0.0f, main_radius*2.0f);
				}//(int i, int j, float stiffness, float give=0.0f, float length=0.0f)
			}
		}
		else
		{
			for (int i = 0; i < ingr_partner_node["properties"]["beadsin"].size(); i++)
			{
				int beadsin = ingr_partner_node["properties"]["beadsin"][i].asInt();
				CreateSpringInter(rope.mIndices[curveI], offset + beadsin, 1, 0.0f, g_params.mRadius);
				if (curveI + 1 < current_nPoints)
					CreateSpringInter(rope.mIndices[curveI + 1], offset + beadsin, 0.5f, 0.0f, g_params.mRadius);
				if (curveI - 1 > 0)
					CreateSpringInter(rope.mIndices[curveI - 1], offset + beadsin, 0.5f, 0.0f, g_params.mRadius);
				cout << "addx " << curveI + 1 << " " << rope.mIndices[curveI] << " " << offset + beadsin << endl;
			}
			for (int i = 0; i < ingr_partner_node["properties"]["beadsout"].size(); i++)
			{
				int beadsout = ingr_partner_node["properties"]["beadsout"][i].asInt();
				int add_id = curveI + curve_ln + curve_var;
				if (add_id >= current_nPoints)
				{
					add_id = add_id - current_nPoints;
				}
				if (add_id < current_nPoints)
				{
					CreateSpringInter(rope.mIndices[add_id], offset + beadsout, 1, 0.0f, g_params.mRadius);
					cout << i << " " << beadsout << "addy " << add_id << " " << rope.mIndices[add_id] << " " << offset + beadsout << endl;
				}
				if (add_id + 1 < current_nPoints)
				{
					CreateSpringInter(rope.mIndices[add_id + 1], offset + beadsout, 0.5f, 0.0f, g_params.mRadius);
					cout << i << " " << beadsout << "addy " << add_id + 1 << " " << rope.mIndices[add_id + 1] << " " << offset + beadsout << endl;
				}
				if (add_id - 1 > 0)
				{
					CreateSpringInter(rope.mIndices[add_id - 1], offset + beadsout, 0.5f, 0.0f, g_params.mRadius);
					cout << i << " " << beadsout << "addy " << add_id - 1 << " " << rope.mIndices[add_id - 1] << " " << offset + beadsout << endl;
				}

			}
		}
	}
};

class Mycoplasma : public Scene
{
public:

	Mycoplasma(const char* name) : Scene(name) {}
    cellPACK * cp;
	flexServer *fserver;

    int base;
    int width;
    int height;
    int depth;
    float main_scale = 1.0f/100.0f;//0.2f;//1.0f/100.0f;//0.5 ? radius is 0.05
    bool shrink = false;
	
	bool server = false;
	bool write_output = false;
	bool grow_fiber = false;
	bool dojitter_steared = true;
	bool dojitter = true;
	bool dojitter_biased = true;
	bool use_threshold_binding = true;
	float dojitter_strength = 1.0f;
	float dojitter_biased_strength = 1.0f;
	float threshold_binding;
    int N = 10;
	float Nsub = 10.0f;
    int group;
	bool dosimulation = false;
    ofstream output;
    ofstream output_bin;

    struct Instance
    {
        Vec3 mTranslation;
        Quat mRotation;
        float mLifetime;

        int mGroup;
        int mParticleOffset;

        int mMeshIndex;
    };

    struct MeshBatch
    {
        GpuMesh* mMesh;
        FlexExtAsset* mAsset;

        std::vector<Matrix44> mInstanceTransforms;
    };

    struct MeshAsset
    {
        string file;
        float scale;
    };

    float mAttractForce;
    float mesh_scale = 0.1f;
    std::vector<MeshBatch> mBatches;

    int mGroupCounter;
    std::vector<Instance> mInstances;

	float mPressure;
	float strength=10.0f;

	std::vector<ClothMesh*> mCloths;
	std::vector<float> mRestVolume;
	std::vector<int> mTriOffset;
	std::vector<int> mTriCount;
	std::vector<float> mOverPressure;
	std::vector<float> mConstraintScale;
	std::vector<float> mSplitThreshold;


    virtual void testFlex(){
        //nt group = 0;
        int phase = NvFlexMakePhase(group++,0);//eFlexPhaseSelfCollide
        //in order to test flex mRadius, solidRadius, collision Distance and spring Length

        //Vec4 center1 = Vec4(0.0f,0.0f,0.0f,1.0f);

        g_buffers->positions.push_back(Vec4(0.0f,0.0f,0.0f,1.0f));
        g_buffers->velocities.push_back(0.0f);
        g_buffers->phases.push_back(NvFlexMakePhase(group,0));

        g_buffers->positions.push_back(Vec4(1.0f,0.0f,0.0f,1.0f));
        g_buffers->velocities.push_back(0.0f);
        g_buffers->phases.push_back(NvFlexMakePhase(group,0));

        g_buffers->positions.push_back(Vec4(0.0f,1.0f,0.0f,1.0f));
        g_buffers->velocities.push_back(0.0f);
        g_buffers->phases.push_back(NvFlexMakePhase(group,0));

        g_buffers->positions.push_back(Vec4(1.0f,1.0f,0.0f,1.0f));
        g_buffers->velocities.push_back(0.0f);
        g_buffers->phases.push_back(NvFlexMakePhase(group,0));

        float d = Length(g_buffers->positions[1]-g_buffers->positions[0]);
        cout << d << endl;
//
//        g_buffers->positions.push_back(Vec4(g_params.mCollisionDistance,g_params.mRadius*2.0f,0.0f,1.0f));
//        g_buffers->velocities.push_back(0.0f);
//        g_buffers->phases.push_back(NvFlexMakePhase(group++,0));
//
//        //same system but with spring
//        int p1 = int(g_buffers->positions.size());
//        g_buffers->positions.push_back(Vec4(0.0f,-g_params.mRadius,0.0f,1.0f));
//        g_buffers->velocities.push_back(0.0f);
//        g_buffers->phases.push_back(NvFlexMakePhase(group++,0));
//
//        int p2 = int(g_buffers->positions.size());
//        g_buffers->positions.push_back(Vec4(g_params.mRadius,-g_params.mRadius,0.0f,1.0f));
//        g_buffers->velocities.push_back(0.0f);
//        g_buffers->phases.push_back(NvFlexMakePhase(group++,0));
//        CreateSpring(p1,p2, 0.1f,0.0f,g_params.mRadius);
//        //CreateSpring(previous[k], current[k], stiffness, give,D)
//        //CreateSpring(previous[k], current[k], stiffness, give,D)
//
//         flexSetParticles(g_flex, &g_buffers->positions[0].x, g_buffers->positions.size(), eFlexMemoryHost);
//         flexSetVelocities(g_flex, &g_buffers->velocities[0].x, g_buffers->velocities.size(), eFlexMemoryHost);
//         flexSetPhases(g_flex, &g_buffers->phases[0], g_buffers->phases.size(), eFlexMemoryHost);
    }

	void AddInflatable(const Mesh* mesh, float overPressure, int phase)
	{
		const int startVertex = g_buffers->positions.size();

		// add mesh to system
		for (size_t i = 0; i < mesh->GetNumVertices(); ++i)
		{
			const Vec3 p = Vec3(mesh->m_positions[i]);

			g_buffers->positions.push_back(Vec4(p.x, p.y, p.z, 1.0f));
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(phase);
		}

		int triOffset = g_triangles.size();
		int triCount = mesh->GetNumFaces();

		mTriOffset.push_back(triOffset / 3);
		mTriCount.push_back(mesh->GetNumFaces());
		mOverPressure.push_back(overPressure);

		for (size_t i = 0; i < mesh->m_indices.size(); i += 3)
		{
			int a = mesh->m_indices[i + 0];
			int b = mesh->m_indices[i + 1];
			int c = mesh->m_indices[i + 2];

			Vec3 n = -Normalize(Cross(mesh->m_positions[b] - mesh->m_positions[a], mesh->m_positions[c] - mesh->m_positions[a]));
			g_triangleNormals.push_back(n);

			g_triangles.push_back(a + startVertex);
			g_triangles.push_back(b + startVertex);
			g_triangles.push_back(c + startVertex);
		}

		// create a cloth mesh using the global positions / indices
		ClothMesh* cloth = new ClothMesh(&g_buffers->positions[0], g_buffers->positions.size(), &g_triangles[triOffset], triCount * 3, 0.8f, 1.0f);

		for (size_t i = 0; i < cloth->mConstraintIndices.size(); ++i)
			g_buffers->springIndices.push_back(cloth->mConstraintIndices[i]);

		g_buffers->springStiffness.insert(g_buffers->springStiffness.end(), cloth->mConstraintCoefficients.begin(), cloth->mConstraintCoefficients.end());
		g_buffers->springLengths.insert(g_buffers->springLengths.end(), cloth->mConstraintRestLengths.begin(), cloth->mConstraintRestLengths.end());

		mCloths.push_back(cloth);

		// add inflatable params
		mRestVolume.push_back(cloth->mRestVolume);
		mConstraintScale.push_back(cloth->mConstraintScale);
	}

    virtual void Initialize()
    {

		g_rigidTranslations.resize(0);
		g_rigidRotations.resize(0);
		g_rigidCoefficients.resize(0);
		g_rigidIndices.resize(0);
		g_rigidLocalPositions.resize(0);
		g_rigidOffsets.resize(0);

		// start index
		g_rigidOffsets.push_back(0);

        main_scale =1.0f/100.0f;
		//or 68 2turn of DNA
		float beads_radius = 11.85f*main_scale;//5.0f*main_scale;//0.08f;//5.0f*main_scale;//smallest sphere size
		g_params.mRadius = beads_radius * 2.0f;

		output_bin.open("../../data/pack_result.bin", ios::out | ios::binary);
		output_bin.close();
			
        //g_meshcellPACK * cp;
        cp = new cellPACK();
		cp->use_rb = true;
        cp->main_radius = beads_radius;//((10.0f/100.0f)*(10.0f/100.0f));
        cp->main_scale = main_scale;
		
		//mycoplasma experiment is
		cp->loadRecipe((cp->mainpath + "recipes\\Mycoplasma1.6_full.json").c_str());
		//cp->loadResults((cp->mainpath + "recipes\\Mycoplasma1.6_full_result_1.json").c_str());
		
		cp->loadResults((cp->mainpath + "recipes\\Mycoplasma1.5_mixed_pdb_fixed.json").c_str());
		cp->placeFibers(false); 
		
		cp->buildMembrane(0.01f, 0.5f, NvFlexMakePhase(99999, 0) );
		//cp->buildMembrane(1.0f, 0.0f, NvFlexMakePhase(99999, 0));
		

		//HIV ? need recipe and result ?
		//cp->loadRecipe((cp->mainpath + "\\results\\BloodHIV1.0_mixed_fixed_nc1.cpr").c_str());
		//cp->loadResults((cp->mainpath + "\\results\\BloodHIV1.0_mixed_fixed_nc1.cpr").c_str());
		//cp->placeFibers(false);
		//cp->buildMembrane(0.05f, 0.5f, NvFlexMakePhase(99999, 0));

		//create spring on every triangle ?
		//david experiment is 
		//cp->loadRecipe((cp->mainpath + "\\recipes\\DNAplectoneme.1.0.json").c_str());
		//cp->placeFiber();

		//maxParticles = cp->maxParticles;
		//cp->haltondistribute(10);// 37000);//carefull not to reach max num particle
         //cp->randomDistribute(5);

        group = cp->iGroupCounter;

         //cp->initProceduralRope();
		 //cp->initDNACircleExperiment();
         

		 //cp->placePartner();
		unsigned int test = (10 << 24) | 15;
		unsigned int testshift = (10 << 24);
		unsigned int r1 = TAKE_N_BITS_FROM(test, 0, 24);
		unsigned int r2 = TAKE_N_BITS_FROM(test, 24, 32);
        printf ("after first24 %u last8 %u shift %u\n",r1,r2,testshift); 
        // output.open ("../../data/pack_result.txt");
        // output << "#nbingredients "<< cp->iBatches.size() << "\n" ;//total nb of ingredients
        // for (int i=0;i<cp->iBatches.size();i++){
        //    output <<"#object " << i << " " << cp->pnames[i]<<" " << cp->iBatches[i].nInstances << "\n" ;
        //    std::cout << i << " " << cp->pnames[i] <<" " << cp->iBatches[i].nInstances << "\n" ;
        // }
        
        //curve ?
        // output.close();
        
        //testFlex();
		
		g_numExtraParticles = 1024 * 1024;
        g_params.mNumPlanes = 0;
        //InitMeshes(radius);
        //populateInstance(radius);

        //Mesh* mesh = ImportMesh(GetFilePathByPlatform("../../data/MMycoideHD.obj").c_str());
        //mesh->Transform(ScaleMatrix(1.0f/500.0f));
        //mesh->CalculateNormals();
        //float main_scale = 1.0f/100.0f;

//        Mesh* mesh = ImportMesh(GetFilePathByPlatform("../../data/MMycoideHD.obj").c_str());
//        mesh->Transform(ScaleMatrix(1.0f/1000.0f));
//        mesh->CalculateNormals();
//        g_staticMesh = mesh;

//        g_staticMesh = ImportMesh(GetFilePathByPlatform("../../data/sphere.ply").c_str());
//        g_staticMesh->Normalize(2.0f); 
//        g_staticMesh->CalculateNormals();

        g_params.mGravity[1] = 0.0f;//-9.f;
		g_params.mDamping = 1.0f;// 3.0f;
        g_params.mRadius = beads_radius * 2;
        //g_params.mSolidRestDistance = g_params.mRadius*2.0f;
        g_windStrength = 0.0f;
        g_windFrequency = 0.0f;

//        g_params.mDynamicFriction = 0.00f;
//        g_params.mFluid = true;
//        g_params.mViscosity = 0.0f;
        g_params.mNumIterations = 5;
//        g_params.mVorticityConfinement = 0.0f;
//        g_params.mAnisotropyScale = 50.0f;
//        g_params.mSmoothing = 1.f;
//        //g_params.mFluidRestDistance = restDistance;
//        g_params.mNumPlanes = 0;
//        g_params.mCohesion = 0.0025f;
//        g_params.mSurfaceTension = 0.0f;
//        g_params.mCollisionDistance = 0.001f;
//        g_params.mRestitution = 0.0f;
//
//        g_params.mRelasxationFactor = 1.0f;
//

		g_params.mFluid = false;
		g_params.mVorticityConfinement = 0.0f;
		//g_params.mFluidRestDistance = g_params.mRadius*0.5;
		g_params.mAnisotropyScale = 2.5f / g_params.mRadius;
		g_params.mSmoothing = 0.5f;
		g_params.mRelaxationFactor = 1.f;
		g_params.mRestitution = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		g_params.mDynamicFriction = 0.25f;
		g_params.mViscosity = 1.5f;
		g_params.mCohesion = 1.1f;
		g_params.mAdhesion = 1.0f;
		g_params.mSurfaceTension = 1.0f;

		g_params.mGravity[1] = 0.0f;

		
		g_lightDistance *= 2.0f;

		// draw options		
		g_drawEllipsoids = true;
        g_params.mMaxSpeed = 100.0f;
          //g_params.mMaxSpeed = 0.5f*g_params.mRadius*g_numSubsteps/g_dt;
//
//        g_numSubsteps = 2;
//
//        g_fluidColor = Vec4(0.2f, 0.6f, 0.9f, 1.0f);
//        g_clearColor = Vec3(0.5f,0.0f,0.0f);
//        g_lightDistance *= 0.85f;
//
//        // draw options
//        g_drawDensity = true;
//        g_drawDiffuse = true;
//        g_drawEllipsoids = false;
//        g_drawPoints = true;
//
          g_wireframe = true;
          g_pointScale = 1.0f;
//        g_blur = 2.0f;
          g_pause = true;
          g_warmup = false;
          g_drawMesh = true;
          g_drawRopes = false;

		  g_params.mDynamicFriction = 0.4f;
		  g_params.mDissipation = 0.0f;
		  
		  g_params.mParticleCollisionMargin = g_params.mRadius*0.05f;
		  g_params.mDrag = 0.0f;
		  g_params.mCollisionDistance = 0.01f;

		  // better convergence with global relaxation factor
		  g_params.mRelaxationMode = eFlexRelaxationGlobal;
		  g_params.mRelaxationFactor = 0.25f;

		  g_windStrength = 0.0f;

		  g_numSubsteps = 5;
		  g_params.mNumIterations = 5;

		  mSplitThreshold.resize(mCloths.size(), 45.0f);

		  // draw options		
		  g_drawPoints = true;
		  g_drawSprings = 0;
		  g_drawCloth = false;
    }

	virtual void startServer()
	{
		fserver = new flexServer();
		fserver->InitializeWinsock();	
	}

	virtual void stopServer()
	{
		fserver->stopServer();
	}

    virtual void DoGui()
    {
        if (imguiCheck("shrink", shrink))
        {
            shrink = !shrink;
        }
        imguiSlider("MaxSpeed", &g_params.mMaxSpeed, 0.0f, 50000.0f, 0.1f);
		imguiSlider("Nsubidiv", &Nsub, 0.0f, 50.0f, 1.0f);
		if (imguiCheck("start server", server))
		{
			server = !server;
			if (server)
				startServer();
			else
				stopServer();
		}
		if (imguiCheck("write output", write_output))
		{
			write_output = !write_output;
		}
		if (imguiCheck("grow fiber update", grow_fiber))
		{
			grow_fiber = !grow_fiber;
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
		imguiSlider("motion strength", &dojitter_strength, 0.0f, 100.0f, 0.1f);
		imguiSlider("center biased strength", &dojitter_biased_strength, 0.0f, 100.0f, 0.1f);
		imguiSlider("binding distance", &threshold_binding, 0.0f, 20.0f, 0.01f);
		imguiSlider("membrane strength", &strength, 0.0f, 200.0f, 1.0f);
		//flexExtSetFFStrength(cp->fcontainer, 0, strength);
    }

    virtual void KeyDown(int key)
    {
        if (key == 'B')
        {
            cp->randomDistribute(5);
        }
		  
        if (key == 'l')
        {
            std::cout << "add point" << endl;
			cp->growFiber((int)Nsub);
        }
		if (key == 'm')
		{
			std::cout << "key m" << endl;
			std::cout << cp->mInstances.size() << endl;
			if (server) sendToClient();
		}
		//if (key == 'n')
		//{
		//	std::cout << "key n" << endl;
			//std::cout << cp->mInstances.size() << endl;
		//	if (server) sendToClientParticles();
		//}
		if (key == 'n')
		{
			writeFiberBinary();
		}
    }

	void sendToClientParticles(){
		//int current = int(g_buffers->positions.size());
		int current = flexGetActiveCount(g_flex);
		//flexGetActive(g_flex, &g_activeIndices[0], eFlexMemoryHost);
		//flexGetParticles(g_flex, &g_buffers->positions[0].x, cp->maxParticles, eFlexMemoryHost);

		printf("current is %d..\n", current);
		//int(g_buffers->positions.size()
		float* P = new float[current];
		int fcount = 0;
		for (int i = 0; i < current; i++){
			P[fcount] = -g_buffers->positions[i].x*(1.0f / main_scale);
			P[fcount + 1] = g_buffers->positions[i].y*(1.0f / main_scale);
			P[fcount + 2] = g_buffers->positions[i].z*(1.0f / main_scale);
			P[fcount + 3] = g_buffers->positions[i].w;
			fcount += 4;
		}
		int size = sizeof(float) * (current * 4);
		printf("size is %d..\n", size);
		char * test = (char*)&P[0];
		std::cout << "send nInst" << endl;
		fserver->sendToClient(test, size);
		// flexGetParticles(g_flex, &g_buffers->positions[0].x, cp->maxParticles, eFlexMemoryHost);
	}

	void sendToClient(){
		
		int nInst = (int)cp->mInstances.size();
		
		//printf("nInst is %d..\n", nInst);
		
		int nptsTotal = 0;
		for (int r = 0; r<g_ropes.size(); r++){
			nptsTotal += g_ropes[r].coarseIndices.size();
		}
		//printf("nInst is %d..\n", nInst);
		//printf("nptsTotal is %d..\n", nptsTotal);
		//uint32_t* header = new uint32_t[2];
		//header[0] = (uint32_t) nInst;
		//header[1] = (uint32_t) nptsTotal;
		
		Vec3* pos = new Vec3[nInst];//offset from surface
		Quat* quat = new Quat[nInst];//
		flexGetRigidTransforms(g_flex, (float*)&quat[0], (float*)&pos[0], eFlexMemoryHost);
		
		int fsize = 4 * nInst * 2 + 2 + nptsTotal * 4;
		float* P = new float[fsize];
		P[0] = (float) nInst;
		P[1] = (float) nptsTotal;
		int fcount = 2;
		for (int i = 0; i<nInst; i++){
			//output << -pos[i].x*(1.0f/main_scale) << " " <<pos[i].y*(1.0f/main_scale)<< " " <<pos[i].z*(1.0f/main_scale)<< " " << cp->mInstances[i].mMeshIndex << "\n";
			float ind = (float)cp->mInstances[i].mMeshIndex;
			P[fcount] = -pos[i].x*(1.0f / main_scale); 
			P[fcount+1] = pos[i].y*(1.0f / main_scale);
			P[fcount+2] = pos[i].z*(1.0f / main_scale);
			P[fcount+3] = ind;
			fcount += 4;
		}
		for (int i = 0; i<nInst; i++){
			//output << quat[i].x << " " << quat[i].y << " " << quat[i].z << " " << quat[i].w << "\n";
			//float ind = (float)cp->mInstances[i].mMeshIndex;
			P[fcount] =     quat[i].x;
			P[fcount + 1] = quat[i].y;
			P[fcount + 2] = quat[i].z;
			P[fcount + 3] = quat[i].w;
			fcount += 4;
		}

		int nRope = g_ropes.size();
		cout << "#nbrope " << g_ropes.size() << "\n"; //total nb instances for this frame
		//output_bin.write((char *)&nRope, sizeof(nRope));
		Vec4* p = new Vec4[nptsTotal];
		int count = 0;
		for (int r = 0; r < nRope; r++){
			int npts = g_ropes[r].coarseIndices.size();
			int cType = cp->maks_fiber[r];
			cout << "#ropesize" << g_ropes[r].coarseIndices.size() << "\n"; //total nb instances for this frame
			//output_bin.write((char *)&npts, sizeof(npts));
			for (int i = 0; i < g_ropes[r].coarseIndices.size(); i++){
				p[count] = g_buffers->positions[g_ropes[r].coarseIndices[i]];
				//float typeind = concatenateInt(r,cType);//(float) (cType << 8 | r);///combine cType and r
				unsigned int  typeind = (cType << 24) | (unsigned int)r;///combine cType and r
				//float p[3] = { -pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale) };
				p[count].w = (float)typeind;
				P[fcount] = -p[count].x*(1.0f / main_scale);
				P[fcount + 1] = p[count].y*(1.0f / main_scale);
				P[fcount + 2] = p[count].z*(1.0f / main_scale);
				P[fcount + 3] = p[count].w;
				fcount += 4;
				count++;
				//output_bin.write((char *)&p, sizeof(float) * 3);
				//output_bin.write((char *)&typeind, sizeof(unsigned int));
			}
		}

		int size = sizeof(float) * fsize;
		printf("size is %d..\n", size);
		//char* buffer = new char[size];
		//memcpy(buffer, &P, size);
		char * test = (char*)&P[0];

		//int header[2] = { nInst, nptsTotal };
		//send(fserver->ClientSocket, (char *)nInst, sizeof(int), 0);
		std::cout << "send nInst" <<  endl;
		//std::string buf;
		//buf.reserve(sizeof(nInst) + sizeof(nptsTotal) + sizeof(float) * 4 * nInst + sizeof(float) * 4 * nInst);
		/*buf = (char *)&nInst;
		buf.append((char *)&nptsTotal);
		buf.append((char *)&pos[0]);
		buf.append((char *)&quat[0]);
		buf.append((char *)&p[0]);*/

		//char *achar = (char*)malloc(strlen((char *)header)  + 1);
		//strcpy(achar, (char *)nInst);
		//strcat(achar, (char *)nptsTotal);
		//strcat(achar, (char *)&pos[0]);
		//strcat(achar, (char *)&quat[0]);
		
		//buf.append((char *)&nptsTotal);
		//buf.append((char *)&header);
		//buf.append((char *)&pos[0]);
		//buf.append((char *)&quat[0]);

		//fserver->sendToClient((char *)&nInst, sizeof(int));
		//fserver->sendToClient((char *)&header, sizeof(int) * 2);

		//buf.append("");
		//std::cout << buf << endl;
		//4+4?
		//fserver->sendToClient(buffer,size);// +sizeof(float) * 4 * nInst + sizeof(float) * 4 * nInst);
		fserver->sendToClient(test,size );
		//fserver->sendToClient(buf.c_str(), size);
		//std::cout << "send nBytes" << sizeof(nInst) + sizeof(nptsTotal) ;// +sizeof(float) * 4 * nInst + sizeof(float) * 4 * nInst << endl;
		//free(achar);
		/*
		fserver->sendToClient("start", 5);
		fserver->sendToClient((char *)&header, sizeof(int) * 2);
		//fserver->sendToClient((char *)&nptsTotal, sizeof(int));
		fserver->sendToClient((char *)&pos[0], sizeof(float) * 4 * nInst);
		fserver->sendToClient((char *)&quat[0], sizeof(float) * 4 * nInst);
		fserver->sendToClient("end", 3);
		*/
		//output_bin.write((char *)&nptsTotal, sizeof(nptsTotal));
	}

    void Warmup()
    {
        printf("Warming up sim..\n");
        //uint32_t numParticles = g_buffers->positions.size();
        //uint32_t maxParticles = numParticles + g_numExtraParticles*g_numExtraMultiplier;
        float ** p;
        float ** v;
        int ** ph;
        float ** no;
        int a;
        // warm it up (relax positions to reach rest density without affecting velocity)send
        FlexParams copy = g_params;
        copy.mNumIterations = 4;

        flexSetParams(g_flex, &copy);

        const int kWarmupIterations = 100;

        for (int i = 0; i < kWarmupIterations; ++i)
        {
            FlexTimers atimers;
            flexExtTickContainer(cp->fcontainer,  0.0001f, 1, &atimers);
            //flexUpdateSolver(g_flex, 0.0001f, 1, &atimers);
            //flexSetVelocities(g_flex, (float*)&g_buffers->velocities[0], maxParticles, eFlexMemoryHost);
        }

        // udpate host copy
        int active_count = cp->numParticles;//flexExtGetParticleCount(cp->fcontainer);

        g_buffers->positions.resize(cp->maxParticles);
        g_buffers->velocities.resize(cp->maxParticles);
        g_buffers->phases.resize(cp->maxParticles);
        g_activeIndices.resize(cp->maxParticles);//flexGetActiveCount(g_flex));
        flexGetActive(g_flex, &g_activeIndices[0], eFlexMemoryHost);
        flexGetParticles(g_flex, &g_buffers->positions[0].x, cp->maxParticles, eFlexMemoryHost);
        //g_buffers->positions;
        //g_activeIndices;
        flexGetVelocities(g_flex, &g_buffers->velocities[0].x, cp->maxParticles, eFlexMemoryHost);
        flexGetPhases(g_flex, &g_buffers->phases[0], cp->maxParticles, eFlexMemoryHost);

        a=10;
        printf("Finished warm up, %d \n",active_count);
    }

    float concatenateInt(int a, int b){
        float c;
        c = (float)b;
        while( c >= 1.0f ) c *= 0.1f; //moving the decimal point (.) to left most
        c = (float)a + c;
        return c; 
    }

	void writeRigidTransform(int exp, int run, bool append = false)
	{
		Vec3* pos;
		Quat* quat;
		int totalNbBody = 0;
		int nInst = cp->mInstances.size();
		if (cp->use_rb){
			pos = new Vec3[cp->mInstances.size()];
			quat = new Quat[cp->mInstances.size()];
			totalNbBody = nInst;
		}
		else {
			for (int i = 0; i < cp->mInstances.size(); i++){
				FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
				totalNbBody += asset->mNumShapes;
			}
			pos = new Vec3[totalNbBody];
			quat = new Quat[totalNbBody];
		}

		flexGetRigidTransforms(g_flex, (float*)&quat[0], (float*)&pos[0], eFlexMemoryHost);

		ofstream of;
		std:string name = "../../data/pack_result";
		if (append){
			name = name + ".txt";
			of.open(name.c_str(), ios::out | ios::app);
			of << "# " << cp->mInstances.size() << " " << exp << " " << run << " " << endl;
		}
		else {
			name = name + "_" + std::to_string(exp) + "_" + std::to_string(run)  + ".txt";
			of.open(name.c_str(), ios::out);
		}
	
		//write position
		for (int i = 0; i<totalNbBody; i++){
			//output << -pos[i].x*(1.0f/main_scale) << " " <<pos[i].y*(1.0f/main_scale)<< " " <<pos[i].z*(1.0f/main_scale)<< " " << cp->mInstances[i].mMeshIndex << "\n";
			int ind = cp->mInstances[0].mMeshIndex;//instance 
			float p[4] = { pos[i].x*(1.0f / main_scale), pos[i].y*(1.0f / main_scale), pos[i].z*(1.0f / main_scale), ind };
			of << ind << " " << p[0] << " " << p[1] << " " << p[2] << " " << quat[i].x << " " << quat[i].y << " " << quat[i].z << " " <<  quat[i].w << endl;
		}
		//write rotation
		of.close();
	}
	
	void writeSoftTransform(int exp, int run, int free, bool append = false)
	{
		ofstream of;
		std:string name = "../../data/pack_result_soft";
		if (append){
			name = name + ".txt";
			of.open(name.c_str(), ios::out | ios::app); 
			of << "# " << cp->mInstances.size() << " " << exp << " " << run << " " <<endl;
		}
		else { 
			name = name + "_" + std::to_string(exp) + "_" + std::to_string(run) + "_" + std::to_string(free) + ".txt";
			of.open(name.c_str(), ios::out); 
		}
		
		//of.open("C:\\Users\\ludov\\OneDrive\\Documents\\cellVIEW\ -\ i\\Assets\\Data\\pack_result_soft.txt", ios::out);
		
		int nInst = cp->mInstances.size();
		for (int i = 0; i < cp->mInstances.size(); i++){
			//grab the particles positions for each isntances
			int poffseti = cp->mInstances[i].mParticleOffset;
			FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			int nPart = asset->mNumParticles;
			for (int j = 0; j < nPart; j++){
				of << i << " " << j << " " << g_buffers->positions[poffseti + j].x*(1.0f / main_scale) << " " << g_buffers->positions[poffseti + j].y*(1.0f / main_scale) << " " << g_buffers->positions[poffseti + j].z*(1.0f / main_scale) << endl;
			}
		}
		of.close();
	}

	void exportPDBSoftTransform(int exp, int run)
	{
		ofstream of;
		std:string name = "../../data/pack_result_soft";
		bool append = false;
		if (append){
			name = name + ".txt";
			of.open(name.c_str(), ios::out | ios::app);
			of << "# " << cp->mInstances.size() << " " << exp << " " << run << " " << endl;
		}
		else {
			name = name + "_" + std::to_string(exp) + "_" + std::to_string(run) + ".pdb";
			of.open(name.c_str(), ios::out);
		}

		//of.open("C:\\Users\\ludov\\OneDrive\\Documents\\cellVIEW\ -\ i\\Assets\\Data\\pack_result_soft.txt", ios::out);

		int nInst = cp->mInstances.size();
		for (int i = 0; i < cp->mInstances.size(); i++){
			//grab the particles positions for each isntances
			int poffseti = cp->mInstances[i].mParticleOffset;
			FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			int nPart = asset->mNumParticles;
			for (int j = 0; j < nPart; j++){
				of << i << " " << j << " " << g_buffers->positions[poffseti + j].x*(1.0f / main_scale) << " " << g_buffers->positions[poffseti + j].y*(1.0f / main_scale) << " " << g_buffers->positions[poffseti + j].z*(1.0f / main_scale) << endl;
			}
		}
		of.close();
	}

	void writeToBinary(std::string postfix="")
	{
		Vec3* pos;
		Quat* quat;
		int nInst = cp->mInstances.size();
		int totalNbBody = 0;
		if (cp->use_rb){
			pos = new Vec3[cp->mInstances.size()];
			quat = new Quat[cp->mInstances.size()];
			totalNbBody = nInst;
		}
		else {
			
			for (int i = 0; i < cp->mInstances.size(); i++){
				FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
				totalNbBody += asset->mNumShapes;
			}
			pos = new Vec3[totalNbBody];
			quat = new Quat[totalNbBody];
		}

		flexGetRigidTransforms(g_flex, (float*)&quat[0], (float*)&pos[0], eFlexMemoryHost);
		std:string fname = "../../data/pack_result" + postfix + ".bin";	
		
		output_bin.open(fname.c_str(), ios::out | ios::app | ios::binary);

		//number of instances
		output_bin.write((char *)&nInst, sizeof(nInst));
		//number of rigidbody
		output_bin.write((char *)&totalNbBody, sizeof(totalNbBody));
		//number of controle points total
		int nptsTotal = 0;
		for (int r = 0; r < g_ropes.size(); r++){
			nptsTotal += g_ropes[r].mIndices.size();
		}
		output_bin.write((char *)&nptsTotal, sizeof(nptsTotal));
		//number of membrane points total
		int npointmembrane = 0;
		if (cp->mask_membrane.size() != 0)
		{
			for (int m = 0; m < cp->mask_membrane.size() / 2; m += 2)
			{
				npointmembrane += cp->mask_membrane[m + 1];
			}
		}
		output_bin.write((char *)&npointmembrane, sizeof(npointmembrane));//nb of membrane point

		//write position
		if (cp->use_rb){
			for (int i = 0; i < cp->mInstances.size(); i++){
				//output << -pos[i].x*(1.0f/main_scale) << " " <<pos[i].y*(1.0f/main_scale)<< " " <<pos[i].z*(1.0f/main_scale)<< " " << cp->mInstances[i].mMeshIndex << "\n";
				float ind = (float)cp->mInstances[i].mMeshIndex;
				float p[4] = { -pos[i].x*(1.0f / main_scale), pos[i].y*(1.0f / main_scale), pos[i].z*(1.0f / main_scale), ind };
				output_bin.write((char *)&p, sizeof(float) * 4);
			}
			output_bin.write((char *)&quat[0], sizeof(float) * 4 * nInst);
		}
		else {
			int count = 0;
			for (int i = 0; i < cp->mInstances.size(); i++){
				float ind = (float)cp->mInstances[i].mMeshIndex;
				FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
				for (int j = 0; j < asset->mNumShapes; ++j)
				{
					float p[4] = { -pos[count].x*(1.0f / main_scale), pos[count].y*(1.0f / main_scale), pos[count].z*(1.0f / main_scale), ind };
					output_bin.write((char *)&p, sizeof(float) * 4);
					count++;
				}
			}
			output_bin.write((char *)&quat[0], sizeof(float) * 4 * totalNbBody);
		}

		//write control point
		int nRope = g_ropes.size();

		for (int r = 0; r<g_ropes.size(); r++){
			int npts = g_ropes[r].mIndices.size();
			float cType = (float)cp->maks_fiber[r];
			float cId = (float)r;
			for (int i = 0; i<npts; i++){
				Vec4 pos = g_buffers->positions[g_ropes[r].mIndices[i]];
				//unsigned int  typeind = (cType << 24) | (unsigned int) r;///combine cType and r
				float p[3] = { -pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale) };
				output_bin.write((char *)&p, sizeof(float) * 3);
				output_bin.write((char *)&cType, sizeof(float));
				output_bin.write((char *)&cId, sizeof(float));
			}
		}

		//write membrane point
		int typemb = 0;
		for (int m = 0; m < cp->mask_membrane.size() / 2; m += 2)
		{
			int start = cp->mask_membrane[m];
			int npoints = cp->mask_membrane[m + 1];
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
		output_bin.close();
	}

	void writeFiberBinary(){
	std:string fname = "C:/Users/ludov/OneDrive/Documents/OnlinePacking_Tobias/cellVIEW-OP/Data/pack_result_fiber.bin";// "../../data/pack_result_fiber.bin";
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
			float cType = (float)cp->maks_fiber[r];
			float cId = (float)r;
			for (int i = 0; i<npts; i++){
				Vec4 pos = g_buffers->positions[g_ropes[r].mIndices[i]];
				//unsigned int  typeind = (cType << 24) | (unsigned int) r;///combine cType and r
				float p[3] = { -pos.x*(1.0f / main_scale), pos.y*(1.0f / main_scale), pos.z*(1.0f / main_scale) };
				output_bin.write((char *)&p, sizeof(float) * 3);
				output_bin.write((char *)&cType, sizeof(float));//fiber_type
				output_bin.write((char *)&cId, sizeof(float));//fiber_id
			}
		}
		output_bin.close();
	}

	void constraintToSurface()
	{
		flexGetRigidTransforms(g_flex, (float*)&g_rigidRotations[0], (float*)&g_rigidTranslations[0], eFlexMemoryHost);
		//loop over instance surface
		for (int i = 0; i < cp->mInstances.size(); i++)
		{
			
			if (cp->mask[i] > 0) //surface
			{
				//cout << "constraint object " << i << " " << cp->mask[i] << endl;
				int ingrIndex = cp->mInstances[i].mMeshIndex;
				
				Vec3 p = g_rigidTranslations[i];
				Quat r = g_rigidRotations[i];
				SDF * fsdf = cp->comp_shape[abs(cp->mask[i]) - 1];
				float spacing = (fsdf->mUpper[0] - fsdf->mLower[0]) / (float)fsdf->mWidth;
				Vec3 botsdf = Vec3(fsdf->mLower[0], fsdf->mLower[1], fsdf->mLower[2]);
				Vec3 topsdf = Vec3(fsdf->mUpper[0], fsdf->mUpper[1], fsdf->mUpper[2]);
				float3 pquery = { (p.x - botsdf.x) / spacing, (p.y - botsdf.y) / spacing, (p.z - botsdf.z) / spacing };
				int3 pquery_index = { (int)pquery.x, (int)pquery.y, (int)pquery.z };
				int dimFF = fsdf->mDepth;
				
				float drigid = SampleSDFX(fsdf->mField, dimFF, (int)pquery.x, (int)pquery.y, (int)pquery.z)*100.0f;
				Vec3 vsdf = SampleSDFGradX(fsdf->mField, dimFF, (int)pquery.x, (int)pquery.y, (int)pquery.z);
				//cout << "constraint object " << cp->iBatches[ingrIndex].ingr_name << " " << drigid << endl;
				Vec3 pal = Vec3(cp->iBatches[ingrIndex].pcpalVectorx, cp->iBatches[ingrIndex].pcpalVectory, cp->iBatches[ingrIndex].pcpalVectorz);
				Vec3 offset = Vec3(cp->iBatches[ingrIndex].offsetx, cp->iBatches[ingrIndex].offsety, cp->iBatches[ingrIndex].offsetz);

				Vec3 rpal = Rotate(r, pal);
				Vec3 roffset = Rotate(r, offset);

				Vec3 axis = Cross(rpal / Length(rpal), vsdf / Length(vsdf));
				float dot = Dot(rpal / Length(rpal), vsdf / Length(vsdf));
				float angle = ACos(dot);
				Quat q = QuatFromAxisAngle(axis, angle);

				Vec3 delta_pos = (vsdf / Length(vsdf))*drigid;
				cp->UpdateInstanceTransform(i, p - delta_pos, Quat(0,0,0,1));
			}
		}
		flexSetParticles(g_flex, &g_buffers->positions[0].x, g_buffers->positions.size(), eFlexMemoryHost);
		//flexSetRigids(g_flex, &g_rigidOffsets[0], &g_rigidIndices[0], (float*)&g_rigidLocalPositions[0], 
		//	g_rigidLocalNormals.size() ? (float*)&g_rigidLocalNormals[0] : NULL, &g_rigidCoefficients[0], 
		//	(float*)&g_rigidRotations[0], (float*)&g_rigidTranslations[0], cp->mInstances.size(), eFlexMemoryHost);
	}

    void Update()
    {
        FlexTimers atimers;
        float ** p;
        float ** v;
        int ** ph;
        float ** no;
        int a;

        //if (!g_pause) {
        //g_drawRopes = false;
        //Vec3* pos = new Vec3[cp->mInstances.size()];//offset from surface
        //Vec3* pcpalVector;//direction rotated to be align
        //Quat* quat= new Quat[cp->mInstances.size()];//
		//constraintToSurface();
		
		//flexSetSprings(g_flex, &g_buffers->springIndices[0], &g_buffers->springLengths[0], &g_buffers->springStiffness[0], g_buffers->springLengths.size(), eFlexMemoryHost);

		if (cp->comp_shape.size()) 
			flexExtSetRigidTransformation(cp->fcontainer, cp->mInstances.size());
         //save packing leftHand
		if ((g_frame % 1) == 0)
		{
			if (server)
				sendToClient();
			if (write_output)
				writeToBinary();
			//		 sendToClientParticles();
			//sendToClient();
		}
		 if ((g_frame % 5) && (grow_fiber)) 
			 cp->growFiber((int)Nsub);

		 //if (((g_frame % 1) == 0) && (write_output)){
		 //	writeToBinary();
         }

    void PostInitialize()
    {
		if (cp->comp_shape.size()) 
		{
			cp->fcontainer = flexExtCreateContainer(g_flex, cp->maxParticles);
			cp->sdfToForceField();
			cp->setForceField();
			flexExtSetRigidTransformation(cp->fcontainer, cp->mInstances.size());
		}
		//distribute
        //distribute
        //cp->randomDistribute(200);
       // FlexSDF* sdf = new FlexSDF[cp->comp_shape.size()];
       // for (int i=0;i<cp->comp_shape.size();i++)
       //     sdf[i] = *cp->comp_shape[0];
        //sdf[1] = *cp->comp_shape[1];
        //flexSetFields(g_flex, sdf, cp->comp_shape.size());
		//flexSetRigids(g_flex, &g_rigidOffsets[0], &g_rigidIndices[0], (float*)&g_rigidLocalPositions[0], 
		//	(float*)&g_rigidLocalNormals[0], &g_rigidCoefficients[0], (float*)&g_rigidRotations[0], 
		//	(float*)&g_rigidTranslations[0], g_rigidOffsets.size() - 1, eFlexMemoryHost);
		
        
        //overwrite the plane
        g_params.mNumPlanes = 0;
        (Vec4&)g_params.mPlanes[0] = Vec4(0.0f, 1.0f, 0.0f, -cp->minExtents.y);
		(Vec4&)g_params.mPlanes[1] = Vec4(0.0f, 0.0f, 1.0f, -cp->minExtents.z);// -cp->minExtents.z);
        (Vec4&)g_params.mPlanes[2] = Vec4(1.0f, 0.0f, 0.0f, -cp->minExtents.x);
        (Vec4&)g_params.mPlanes[3] = Vec4(-1.0f, 0.0f, 0.0f, cp->maxExtents.x);
		(Vec4&)g_params.mPlanes[4] = Vec4(0.0f, 0.0f, -1.0f, cp->maxExtents.z);// cp->maxExtents.z);
        (Vec4&)g_params.mPlanes[5] = Vec4(0.0f, -1.0f, 0.0f, cp->maxExtents.y);
        //gColors;
        //cp->proceduralRope(10);
		//sendToClient();
		//flexSetInflatables(g_flex, &mTriOffset[0], &mTriCount[0], &mRestVolume[0], 
		//	&mOverPressure[0], &mConstraintScale[0], mCloths.size(), eFlexMemoryHost);

    }

    void Draw(int pass)
    {
        if (!g_drawMesh)
            return;
        for (int i=0;i < cp->comp_mesh.size();i++)
        {
            DrawMesh(cp->comp_mesh[i], g_meshColor);
        }
    }

};

class DNAplectoneme : public Mycoplasma
{
public:

	DNAplectoneme(const char* name) : Mycoplasma(name) {}

	virtual void Initialize()
	{
		g_rigidTranslations.resize(0);
		g_rigidRotations.resize(0);
		g_rigidCoefficients.resize(0);
		g_rigidIndices.resize(0);
		g_rigidLocalPositions.resize(0);
		g_rigidOffsets.resize(0);

		// start index
		g_rigidOffsets.push_back(0);

		main_scale = 1.0f / 100.0f;
		//or 68 2turn of DNA
		float beads_radius = 11.85f*main_scale;//5.0f*main_scale;//0.08f;//5.0f*main_scale;//smallest sphere size
		g_params.mRadius = beads_radius * 2.0f;

		output_bin.open("../../data/pack_result.bin", ios::out | ios::binary);
		output_bin.close();

		//g_mesh
		cp = new cellPACK();

		cp->main_radius = beads_radius;//((10.0f/100.0f)*(10.0f/100.0f));
		cp->main_scale = main_scale;

		cp->loadRecipe((cp->mainpath + "\\recipes\\DNAplectoneme.1.0.json").c_str());
		cp->placeFibers(true);
		//cp->initProceduralRope();
		//cp->initDNACircleExperiment();
		//cp->placePartner();

		group = cp->iGroupCounter;
		g_numExtraParticles = 1024 * 1024;
		g_params.mNumPlanes = 0;
		g_params.mGravity[1] = 0.0f;//-9.f;
		g_params.mDamping = 1.0f;// 3.0f;
		g_params.mRadius = beads_radius * 2;
		//g_params.mSolidRestDistance = g_params.mRadius*2.0f;
		g_windStrength = 0.0f;
		g_windFrequency = 0.0f;

		//        g_params.mDynamicFriction = 0.00f;
		//        g_params.mFluid = true;
		//        g_params.mViscosity = 0.0f;
		g_params.mNumIterations = 5;
		//        g_params.mVorticityConfinement = 0.0f;
		//        g_params.mAnisotropyScale = 50.0f;
		//        g_params.mSmoothing = 1.f;
		//        //g_params.mFluidRestDistance = restDistance;
		//        g_params.mNumPlanes = 0;
		//        g_params.mCohesion = 0.0025f;
		//        g_params.mSurfaceTension = 0.0f;
		//        g_params.mCollisionDistance = 0.001f;
		//        g_params.mRestitution = 0.0f;
		//
		//        g_params.mRelasxationFactor = 1.0f;
		//

		g_params.mFluid = false;
		g_params.mVorticityConfinement = 0.0f;
		//g_params.mFluidRestDistance = g_params.mRadius*0.5;
		g_params.mAnisotropyScale = 2.5f / g_params.mRadius;
		g_params.mSmoothing = 0.5f;
		g_params.mRelaxationFactor = 1.f;
		g_params.mRestitution = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		g_params.mDynamicFriction = 0.25f;
		g_params.mViscosity = 1.5f;
		g_params.mCohesion = 1.1f;
		g_params.mAdhesion = 1.0f;
		g_params.mSurfaceTension = 1.0f;

		g_params.mGravity[1] = 0.0f;


		g_lightDistance *= 2.0f;

		// draw options		
		g_drawEllipsoids = true;
		g_params.mMaxSpeed = 100.0f;
		//g_params.mMaxSpeed = 0.5f*g_params.mRadius*g_numSubsteps/g_dt;
		//
		//        g_numSubsteps = 2;
		//
		//        g_fluidColor = Vec4(0.2f, 0.6f, 0.9f, 1.0f);
		//        g_clearColor = Vec3(0.5f,0.0f,0.0f);
		//        g_lightDistance *= 0.85f;
		//
		//        // draw options
		//        g_drawDensity = true;
		//        g_drawDiffuse = true;
		//        g_drawEllipsoids = false;
		//        g_drawPoints = true;
		//
		g_wireframe = true;
		g_pointScale = 1.0f;
		//        g_blur = 2.0f;
		g_pause = true;
		g_warmup = false;
		g_drawMesh = true;
		g_drawRopes = false;

		g_params.mDynamicFriction = 0.4f;
		g_params.mDissipation = 0.0f;

		g_params.mParticleCollisionMargin = g_params.mRadius*0.05f;
		g_params.mDrag = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		// better convergence with global relaxation factor
		g_params.mRelaxationMode = eFlexRelaxationGlobal;
		g_params.mRelaxationFactor = 0.25f;

		g_windStrength = 0.0f;

		g_numSubsteps = 5;
		g_params.mNumIterations = 5;

		mSplitThreshold.resize(mCloths.size(), 45.0f);

		// draw options		
		g_drawPoints = false;
		g_drawSprings = 0;
		g_drawCloth = false;
	}
};

class HIVIntegrase : public Mycoplasma
{
public:

	HIVIntegrase(const char* name) : Mycoplasma(name) {}
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

	std::vector<int> record_free;
	std::vector<float> record_size;

	int stop_criterion_time=100; // how many last frame to store and check
	int count_stop=0;
	int nexp = 0;
	int current_Exp = 0; // how many experiment?
	int nrun = 5;
	int current_run=0;
	std::map <int, std::vector<float>> exp_params;//carry other the experiments parameters
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
			writeSoftTransform(current_Exp, current_run, 0,false);	//#the beads
			writeRigidTransform(current_Exp, current_run, false);		//the shape rigid bodu
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
		if ((mIntegrase[mind1].ccd_occupied[ccd_id1]!=-1) || (mIntegrase[mind2].ctd_occupied[ctd_id2]!=-1)) doit = false;
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
				//cout << "bound " << id1 << " " << n << " " << poffseti + ccd_ctd[id1][n] << " " << id2 << " " << m << " " << poffsetj + ctd_ccd[id2][m] << " " << g_params.mRadius << " " << g_buffers->springLengths.size() << " " g_buffers->springIndices << endl;
				strings_indices.push_back(g_buffers->springLengths.size());
				CreateSpringInter(poffseti + ccd_ctd[id1][n], poffsetj + ctd_ccd[id2][m], 0.01f, 0.0f, g_params.mRadius); //float give = 0.0f, float length = 0.0f
			}
		}
		mIntegrase[i].ccd_occupied[id1] = j;
		mIntegrase[j].ctd_occupied[id2] = i;
		mIntegrase[i].ccd_ctdind[id1] = id2;
		mIntegrase[j].ctd_ccdind[id2] = id1;
		mIntegrase[i].spring_offset[id1] = sp_offset;
		//flexSetSprings(g_flex, &g_buffers->springIndices[0], &g_buffers->springLengths[0], &g_buffers->springStiffness[0], g_buffers->springLengths.size(), eFlexMemoryHost);
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
				integrase.ccd_ctdind[i]=-1;
				integrase.ctd_ccdind[i]=-1;
				integrase.spring_offset[i]=-1;
				integrase.mini_distance[i]=99999.0f;
				integrase.mini_indices[i]=-1;
			}
			integrase.nb_site_occupied = 0;
			mIntegrase.push_back(integrase);
			free_mIntegrase_CCD.push_back(i);
			free_mIntegrase_CTD.push_back(i);
		}
		cout << " nb integrase " << mIntegrase.size() << endl;
	}
	
	virtual void reportBinding(bool report=true){
		std::vector<int> reports_count = {0,0,0};
		bool bounded = false;
		bool founded = false;
		int bound_order = 0;
		int count_free = 0;
		int unsatisfied = 0;
		for (int i = 0; i < int(cp->mInstances.size()) - 1; ++i)
		{
			FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			int poffseti = cp->mInstances[i].mParticleOffset;
			if (report) {
				//cout << i << " CCD " << mIntegrase[i].ccd_occupied[0] << " " << mIntegrase[i].ccd_occupied[1];//miniD >?
				//cout << " CTD " << mIntegrase[i].ctd_occupied[0] << " " << mIntegrase[i].ctd_occupied[1] << endl;
			}
			if (mIntegrase[i].ccd_occupied[0]==-1) count_free++;
			if (mIntegrase[i].ccd_occupied[1]==-1) count_free++;
			continue;
			bound_order = 0;
			for (int j = i + 1; j < int(cp->mInstances.size()); j++) {
				//build spring ccd-ctd
				if (bound_order==2) continue;
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
								if (D <= g_params.mRadius) {
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
				for (int k = 0 ; k < bound_order;k++)
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
		int poffsetj= cp->mInstances[j].mParticleOffset;
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

	virtual int checkDistance(bool break_spring=false){
		//for every occupied pair, do we satsify the distance.  If not  remove the spring ?
		//should have store the spring indice for each site as well, or just the offset indice as we know the number of spring.
		int unsatisfied_count = 0;
		int nspring = ccd_ctd[0].size()*ctd_ccd[0].size()+1;
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
					if (D > g_params.mRadius*2.0f){
						unsatisfied_count++;
						if (break_spring){
							//remove those strings_indices from g_buffers->springIndices
							int indice_sp = mIntegrase[i].spring_offset[o];
							if (strings_indices[indice_sp] == -1) continue;
							cout << "remove " << indice_sp << " " << strings_indices[indice_sp] << " " << strings_indices[indice_sp] * 2 << endl;
							g_buffers->springIndices.erase(g_buffers->springIndices.begin() + strings_indices[indice_sp]*2, g_buffers->springIndices.begin() + strings_indices[indice_sp]*2 + nspring * 2);
							g_buffers->springLengths.erase(g_buffers->springLengths.begin() + strings_indices[indice_sp], g_buffers->springLengths.begin() + strings_indices[indice_sp] + nspring);
							g_buffers->springStiffness.erase(g_buffers->springStiffness.begin() + strings_indices[indice_sp], g_buffers->springStiffness.begin() + strings_indices[indice_sp] + nspring);
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
		flexSetSprings(g_flex, &g_buffers->springIndices[0], &g_buffers->springLengths[0], &g_buffers->springStiffness[0], g_buffers->springLengths.size(), eFlexMemoryHost);
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
		float radius = g_params.mRadius*0.5f;
		float stiffness = 1.0f;
		std::vector<int> CTD1 = {1, 2, 3, 4, 5, 9, 10, 16, 18, 28, 31, 39, 43, 44 };
		std::vector<int> CCD1 = {0, 6, 7, 11, 12, 14, 15, 19, 22, 24, 26, 33, 35, 37, 40, 46, 47 };
		std::vector<int> NTD1 = { 8, 13, 20, 23, 25, 27, 29, 30, 34, 36, 38, 41, 42, 45, 48, 49 };
		std::vector<int> LinkerA = { 32, 28 };/// { 40, 21, 17, 32, 28 };
		std::vector<int> LinkerB = { 32 + 50, 28 + 50 };// { 45, 46 };
		std::vector<int> dimer = { 8, 75, 25, 58 };// { 11, 12, 14, 15, 22, 26, 33, 47, 8, 23, 25, 34, 49 };
		std::vector<int> CCD = { 0, 6, 7, 11, 12, 14, 15, 19, 22, 24, 26, 33, 35, 37, 40, 46, 47,
			50, 56, 57, 61, 62, 64, 65, 69, 72, 74, 76, 83, 85, 87, 90, 96, 97 };
		std::vector<int> NTD = { 8, 13, 20, 23, 25, 27, 29, 30, 34, 36, 38, 41, 42, 45, 48, 49 ,
			58, 63, 70, 73, 75, 77, 79, 80, 84, 86, 88, 91, 92, 95, 98, 99 };
		std::vector<int> CCDNTD = { 0, 6, 7, 11, 12, 14, 15, 17, 19, 21,22, 24, 26, 32, 33, 35, 37, 40, 46, 47,
			50, 56, 57, 61, 62, 64, 65, 69, 72, 74, 76, 83, 85, 87, 90, 96, 97, 8, 13, 20, 23, 25, 27, 29, 30, 34, 36, 38, 41, 42, 45, 48, 49,
			58, 63, 70, 73, 75, 77, 79, 80, 84, 86, 88, 91, 92, 95, 98, 99, 71, 67, 82
		};

		FlexExtAsset* asset = cp->iBatches[0].mAsset;
		std::vector<int> springIndices;
		std::vector<float> springLengths;
		std::vector<float> springStiffness;
		int numLinks = 0;
		std::vector<int> springIndicestmp;
		std::vector<float> springLengthstmp;
		std::vector<float> springStiffnesstmp;
		// create links between particles

		int numLinksTmp = CreateLinksLocal(asset->mParticles, CTD1.size(), CTD1,
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
		numLinksTmp = CreateLinksLocal(asset->mParticles, CCDNTD.size(), CCDNTD,
			springIndices, springLengths, springStiffness,
			radius, stiffness);
		numLinks += numLinksTmp;
		

		//vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
		//create spring for Linker 
		for (int l = 0; l < LinkerA.size()-1; l++){
			springIndices.push_back(LinkerA[l]);
			springIndices.push_back(LinkerA[l + 1]);
			springLengths.push_back(Length(Vec3(&asset->mParticles[LinkerA[l] * 4]) - Vec3(&asset->mParticles[LinkerA[l + 1] * 4])));
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
		numLinksTmp = CreateLinksLocal(asset->mParticles, CTD1.size(), CTD1,
			springIndices, springLengths, springStiffness,
			radius, stiffness, 50);
		numLinks += numLinksTmp;

		for (int l = 0; l < LinkerA.size() - 1; l++){
			springIndices.push_back(LinkerA[l] + 50);
			springIndices.push_back(LinkerA[l + 1] + 50);
			springLengths.push_back(Length(Vec3(&asset->mParticles[(LinkerA[l] + 50) * 4]) - Vec3(&asset->mParticles[(LinkerA[l + 1] + 50) * 4])));
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
			asset->mSpringIndices = new int[numLinks * 2];
			memcpy(asset->mSpringIndices, &springIndices[0], sizeof(int)*springIndices.size());

			asset->mSpringCoefficients = new float[numLinks];
			memcpy(asset->mSpringCoefficients, &springStiffness[0], sizeof(float)*numLinks);

			asset->mSpringRestLengths = new float[numLinks];
			memcpy(asset->mSpringRestLengths, &springLengths[0], sizeof(float)*numLinks);

			asset->mNumSprings = numLinks;
		}

		//overwrite teh shapes
		int numClusters=5;

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
			center += Vec3(asset->mParticles[CCDNTD[i] * 4 + 0], asset->mParticles[CCDNTD[i] * 4 + 1], asset->mParticles[CCDNTD[i] * 4 + 2]);
		}
		clusterOffsets.push_back(int(clusterIndices.size()));
		clusterPositions.push_back(center / CCDNTD.size());

		center = Vec3(0, 0, 0);
		for (int i = 0; i < CTD1.size(); ++i)
		{
			clusterIndices.push_back(CTD1[i]);
			center += Vec3(asset->mParticles[CTD1[i] * 4 + 0], asset->mParticles[CTD1[i] * 4 + 1], asset->mParticles[CTD1[i] * 4 + 2]);
		}
		clusterPositions.push_back(center / CCDNTD.size());
		clusterOffsets.push_back(int(clusterIndices.size()));
		center = Vec3(0, 0, 0);
		for (int i = 0; i < CTD1.size(); ++i){
			clusterIndices.push_back(CTD1[i] + 50);
			center += Vec3(asset->mParticles[(CTD1[i] + 50) * 4 + 0], asset->mParticles[(CTD1[i] + 50) * 4 + 1], asset->mParticles[(CTD1[i] + 50) * 4 + 2]);
		}
		clusterPositions.push_back(center / CCDNTD.size());
		clusterOffsets.push_back(int(clusterIndices.size()));
		
		
		center = Vec3(0, 0, 0);
		for (int i = 0; i < LinkerA.size(); ++i){
			clusterIndices.push_back(LinkerA[i]);
			center += Vec3(asset->mParticles[(LinkerA[i]) * 4 + 0], asset->mParticles[(LinkerA[i]) * 4 + 1], asset->mParticles[(LinkerA[i]) * 4 + 2]);
		}
		clusterPositions.push_back(center / LinkerA.size());
		clusterOffsets.push_back(int(clusterIndices.size()));

		center = Vec3(0, 0, 0);
		for (int i = 0; i < LinkerB.size(); ++i){
			clusterIndices.push_back(LinkerB[i]);
			center += Vec3(asset->mParticles[(LinkerB[i]) * 4 + 0], asset->mParticles[(LinkerB[i]) * 4 + 1], asset->mParticles[(LinkerB[i]) * 4 + 2]);
		}
		clusterPositions.push_back(center / LinkerA.size());
		clusterOffsets.push_back(int(clusterIndices.size()));


		// assign shapes
		clusterCoefficients.resize(numClusters, clusterStiffness);

		asset->mShapeIndices = new int[clusterIndices.size()];
		memcpy(asset->mShapeIndices, &clusterIndices[0], sizeof(int)*clusterIndices.size());

		asset->mShapeOffsets = new int[numClusters];
		memcpy(asset->mShapeOffsets, &clusterOffsets[0], sizeof(int)*numClusters);

		asset->mShapeCenters = new float[numClusters * 3];
		memcpy(asset->mShapeCenters, &clusterPositions[0], sizeof(float)*numClusters * 3);

		asset->mShapeCoefficients = new float[numClusters];
		memcpy(asset->mShapeCoefficients, &clusterCoefficients[0], sizeof(float)*numClusters);

		asset->mNumShapeIndices = int(clusterIndices.size());
		asset->mNumShapes = numClusters;
		cout << "numCluster " << numClusters << " " << asset->mNumShapeIndices << endl;
	}

	virtual void overwriteSprings(){
		float radius = g_params.mRadius*0.75f;
		float stiffness = 1.0f;
		//body come from the mapping in the sphere tree file

		IngredientSphereTree ing_spheres = cp->mIngrSphereTree[0];
		FlexExtAsset* asset = cp->iBatches[0].mAsset;
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
			for (int m = map_start; m < (map_start+map_count); m++){
				mapping.push_back(m);
				clusterIndices.push_back(m);
				//cout << "cluster " << i << " ptid " << m << endl;
				center += Vec3(asset->mParticles[m * 4 + 0], asset->mParticles[m * 4 + 1], asset->mParticles[m * 4 + 2]);
			}
			//cout << "start " << map_start << " " << map_count << " " << ing_spheres.LevelStarts[1] << clusterIndices.size() << endl;
			//create link for all domain except 0 and 3 
			if ((i != 0)&(i != 3)){
				int numLinksTmp = CreateLinksLocal(asset->mParticles, mapping.size(), mapping,
					springIndices, springLengths, springStiffness,
					radius, stiffness);
				numLinks += numLinksTmp;
			}
			clusterOffsets.push_back(int(clusterIndices.size()));
			clusterPositions.push_back(center / map_count);// ing_spheres.LevelPoints[ing_spheres.LevelStarts[1]]);
		}

		/*integrase dimer attachement */
		
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
			center += Vec3(asset->mParticles[m * 4 + 0], asset->mParticles[m * 4 + 1], asset->mParticles[m * 4 + 2]);
		}
		for (int m = map_start2; m < (map_start2 + map_count2); m++){
			mapping.push_back(m);
			//clusterIndices.push_back(m);
			center += Vec3(asset->mParticles[m * 4 + 0], asset->mParticles[m * 4 + 1], asset->mParticles[m * 4 + 2]);
		}
		int numLinksTmp = CreateLinksLocal(asset->mParticles, mapping.size(), mapping,
			springIndices, springLengths, springStiffness,
			radius+0.01f, stiffness);
		numLinks += numLinksTmp;
		

		//Linker between body, 
		for (int i = 0; i < ing_spheres.LevelCounts[1]-1; i++){
			//link to next one
			//skip the middle one because homodimer
			if (i == 2) continue;
			int last = ing_spheres.LevelMapping[start + i].x + ing_spheres.LevelMapping[start + i].y - 1;
			int next = ing_spheres.LevelMapping[start + i + 1].x;
			springIndices.push_back(last);
			springIndices.push_back(next);
			springLengths.push_back(Length(Vec3(&asset->mParticles[last * 4]) - Vec3(&asset->mParticles[next * 4])));
			springStiffness.push_back(stiffness);
			numLinks++;
		}

		// assign links
		if (numLinks)
		{
			asset->mSpringIndices = new int[numLinks * 2];
			memcpy(asset->mSpringIndices, &springIndices[0], sizeof(int)*springIndices.size());

			asset->mSpringCoefficients = new float[numLinks];
			memcpy(asset->mSpringCoefficients, &springStiffness[0], sizeof(float)*numLinks);

			asset->mSpringRestLengths = new float[numLinks];
			memcpy(asset->mSpringRestLengths, &springLengths[0], sizeof(float)*numLinks);

			asset->mNumSprings = numLinks;
		}

		// assign shapes
		clusterCoefficients.resize(numClusters, clusterStiffness);

		asset->mShapeIndices = new int[clusterIndices.size()];
		memcpy(asset->mShapeIndices, &clusterIndices[0], sizeof(int)*clusterIndices.size());

		asset->mShapeOffsets = new int[numClusters];
		memcpy(asset->mShapeOffsets, &clusterOffsets[0], sizeof(int)*numClusters);

		asset->mShapeCenters = new float[numClusters * 3];
		memcpy(asset->mShapeCenters, &clusterPositions[0], sizeof(float)*numClusters * 3);

		asset->mShapeCoefficients = new float[numClusters];
		memcpy(asset->mShapeCoefficients, &clusterCoefficients[0], sizeof(float)*numClusters);

		asset->mNumShapeIndices = int(clusterIndices.size());
		asset->mNumShapes = numClusters;
		cout << "numCluster " << numClusters << endl;
	}
	
	virtual void Initialize()
	{
		/*
		#new tobind are indices based on symmetry:
		#CCD->CTD1 - 11,12,83
		#CCD->CTD2+50 - 61,62,33
		#CTD->CCD1 - 3,4,10,39,43
		#CTD->CCD2+50 - 53,54,60,89,93
		ccd_ctd_1 = { 83, 71, 56, 45, 11 };//41
		ccd_ctd_2 = { 77, 36, 35, 18, 8 };
		ctd_ccd_1 = { 3, 58, 64 };
		ctd_ccd_2 = { 30, 88, 99 };
		*/

		bool ignore_comp = true;
		/* integrase initialization */
		std::string recipe = "HIV_IN_XP.json";
		std::string results = "INT_50_random.json";
		std::string wrkDir = "../../data/";
		ccd_ctd_1 = { 24, 25, 79 };
		ccd_ctd_2 = { 74, 75, 29 };
		ccd_ctd.push_back(ccd_ctd_1);
		ccd_ctd.push_back(ccd_ctd_2);
		
		ctd_ccd_1 = { 40, 41, 42, 43, 49 };
		ctd_ccd_2 = { 90, 91, 92, 93, 99 };
		ctd_ccd.push_back(ctd_ccd_1);
		ctd_ccd.push_back(ctd_ccd_2);
		
		/* validation experiement setup */
		/*
		std::string recipe = "NanoCage.json";
		std::string results = "nanoCage_results.json";
		std::string wrkDir = "../../../ProteinNanocage/";
		ccd_ctd.push_back({ 33, 31, 30, 29, 28 });
		ccd_ctd.push_back({ 6, 7, 16, 18 });
		ccd_ctd.push_back({ 2, 5, 15 });

		ctd_ccd.push_back({ 33, 31, 30, 29, 28 });
		ctd_ccd.push_back({ 2, 5, 15 });
		ctd_ccd.push_back({ 6, 7, 16, 18 });
		*/

		g_rigidTranslations.resize(0);
		g_rigidRotations.resize(0);
		g_rigidCoefficients.resize(0);
		g_rigidIndices.resize(0);
		g_rigidLocalPositions.resize(0);
		g_rigidOffsets.resize(0);

		// start index
		g_rigidOffsets.push_back(0);

		main_scale = 1.0f / 100.0f;
		
		float beads_radius = 6.5f*main_scale;
		g_params.mRadius = beads_radius * 2.0f;
		threshold_binding = g_params.mRadius*10.0f;
		output_bin.open("../../data/pack_result.bin", ios::out | ios::binary);
		output_bin.close();
		threshold_binding = g_params.mRadius*2.0f;
		//setup experiments
		exp_params[0] = { 1.0f, 1.0f, g_params.mRadius*2.0f };
		exp_params[1] = { 100.0f, 0.0f, g_params.mRadius*100.0f };//no bias toward center but higher distance and random force
		exp_params[2] = { 100.0f, 0.0f, g_params.mRadius*10000.0f };//higher distance or infinite distance ?
		nexp = 3;
		nrun = 3;
		//g_mesh
		
		cp = new cellPACK();
		cp->use_rb = false;

		cp->main_radius = beads_radius;//((10.0f/100.0f)*(10.0f/100.0f));
		cp->main_scale = main_scale;
		cp->mainpath = wrkDir;// HIV_IN_XP.json"
		cp->datapath = wrkDir;// HIV_IN_XP.json"
		cp->loadRecipe(wrkDir+recipe, ignore_comp);
		//cp->loadRecipe((cp->mainpath + "HIV_IN_XP.1.0.json").c_str(), ignore_comp);
		overwriteSprings();

		cp->loadResults(wrkDir+results);
		if (!ignore_comp)cp->buildMembrane(1.0f, 1.0f, NvFlexMakePhase(99999, 0));

		setupBoundsAndLinker();

		group = cp->iGroupCounter;
		g_numExtraParticles = 1024 * 1024;
		g_params.mNumPlanes = 0;
		g_params.mGravity[1] = 0.0f;//-9.f;
		g_params.mDamping = 3.0f;// 3.0f;
		g_params.mRadius = beads_radius * 2;
		//g_params.mSolidRestDistance = g_params.mRadius*2.0f;
		g_windStrength = 0.0f;
		g_windFrequency = 0.0f;
		//        g_params.mDynamicFriction = 0.00f;
		//        g_params.mFluid = true;
		//        g_params.mViscosity = 0.0f;
		g_params.mFluid = false;
		g_params.mVorticityConfinement = 0.0f;
		//g_params.mFluidRestDistance = g_params.mRadius*0.5;
		g_params.mAnisotropyScale = 2.5f / g_params.mRadius;
		g_params.mSmoothing = 0.5f;
		g_params.mRelaxationFactor = 1.f;
		g_params.mRestitution = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		g_params.mDynamicFriction = 0.25f;
		g_params.mViscosity = 1.5f;
		g_params.mCohesion = 1.1f;
		g_params.mAdhesion = 1.0f;
		g_params.mSurfaceTension = 1.0f;

		g_params.mGravity[1] = 0.0f;

		g_lightDistance *= 2.0f;

		// draw options		
		g_drawEllipsoids = true;
		g_params.mMaxSpeed = 100.0f;

		g_wireframe = true;
		g_pointScale = 1.0f;
		//        g_blur = 2.0f;
		g_pause = true;
		g_warmup = false;
		g_drawMesh = true;
		g_drawRopes = false;
		g_drawPoints = true;
		g_params.mDynamicFriction = 0.4f;
		
		g_params.mDissipation = 0.0f;
		g_params.mNumIterations = 4;
		g_params.mViscosity = 0.0f;
		g_params.mDrag = 0.0f;
		g_params.mLift = 0.0f;

		g_params.mParticleCollisionMargin = g_params.mRadius*0.05f;
		g_params.mDrag = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		// better convergence with global relaxation factor
		g_params.mRelaxationMode = eFlexRelaxationGlobal;
		g_params.mRelaxationFactor = 0.25f;

		g_windStrength = 0.0f;

		g_numSubsteps = 5;
		g_params.mNumIterations = 5;

		mSplitThreshold.resize(mCloths.size(), 45.0f);

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
		g_buffers->springIndices.erase(g_buffers->springIndices.begin() + strings_indices[0]*2, g_buffers->springIndices.begin() + strings_indices[strings_indices.size()-1] * 2 + 2 );
		g_buffers->springLengths.erase(g_buffers->springLengths.begin() + strings_indices[0], g_buffers->springLengths.begin() + strings_indices[strings_indices.size()-1]+1 );
		g_buffers->springStiffness.erase(g_buffers->springStiffness.begin() + strings_indices[0], g_buffers->springStiffness.begin() + strings_indices[strings_indices.size() - 1]+1 );
		flexSetSprings(g_flex, &g_buffers->springIndices[0], &g_buffers->springLengths[0], &g_buffers->springStiffness[0], g_buffers->springLengths.size(), eFlexMemoryHost);

		strings_indices.clear();
		mIntegrase.clear();
		free_mIntegrase_CCD.clear();
		free_mIntegrase_CTD.clear();
		setupBoundsAndLinker();
		//replace object to original position
		cp->loadResults("INT_50_random.json",true);
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

	void jitter(bool allparticles){
		//pick one particle per instance, apply some velocity change
		//stay in vicinity of sphere
		float weight = 1;
		for (int i = 0; i < int(cp->mInstances.size()); ++i)
		{
			if (!allparticles){
				//FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
				int poffseti = cp->mInstances[i].mParticleOffset + 19;

				//g_buffers->velocities[poffseti] += RandomUnitVector()*100.0f;
				//g_buffers->velocities[poffseti] -= Vec3(g_buffers->positions[poffseti]);
				Vec3 toward_center = -Vec3(g_buffers->positions[poffseti]);
				if (Length(toward_center) < 1.5f) weight = 0;
				else weight = dojitter_biased_strength;//if (Length(toward_center) > 8.0f) 
				//Vec3 p = Lerp(Vec3(g_buffers->positions[poffseti]), Vec3(g_buffers->positions[poffseti]) + Normalize(toward_center)*weight + RandomUnitVector() * 20, 0.8f);
				Vec3 p = Lerp(Vec3(g_buffers->positions[poffseti]), Vec3(g_buffers->positions[poffseti]) + RandomUnitVector()*dojitter_strength, 0.8f);
				if (dojitter_biased){
					p = Lerp(Vec3(g_buffers->positions[poffseti]), Vec3(g_buffers->positions[poffseti]) + Normalize(toward_center)*weight + RandomUnitVector()*dojitter_strength, 0.8f);
				}
				Vec3 delta = p - Vec3(g_buffers->positions[poffseti]);
				/*
				g_buffers->positions[poffseti].x = p.x;
				g_buffers->positions[poffseti].y = p.y;
				g_buffers->positions[poffseti].z = p.z;
				*/
				g_buffers->velocities[poffseti].x = delta.x / g_dt;
				g_buffers->velocities[poffseti].y = delta.y / g_dt;
				g_buffers->velocities[poffseti].z = delta.z / g_dt;
			}
			else {
				int poffseti = cp->mInstances[i].mParticleOffset;
				FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
				Vec3 toward_center = -Vec3(g_buffers->positions[poffseti]);
				if (Length(toward_center) < 1.5f) weight = 0;
				else weight = dojitter_biased_strength/100.0f;//if (Length(toward_center) > 8.0f) 
				Vec3 p = Lerp(Vec3(g_buffers->positions[poffseti]), Vec3(g_buffers->positions[poffseti]) + RandomUnitVector()*(dojitter_strength/100.0f), 0.8f);
				if (dojitter_biased){
					p = Lerp(Vec3(g_buffers->positions[poffseti]), Vec3(g_buffers->positions[poffseti]) + Normalize(toward_center)*weight + RandomUnitVector()*(dojitter_strength / 100.0f), 0.8f);
				}
				Vec3 delta = p - Vec3(g_buffers->positions[poffseti]);
				Quat q = QuatFromAxisAngle(UniformSampleSphere(), Randf()*k2Pi);
				//Rotate(g_rigidRotations[rigidIndex], localPos)
				//rotate the point and apply the forice
				for (int j = 0; j < asset->mNumParticles; j++){
					Vec3 p_rot = Rotate(q, Vec3(g_buffers->positions[poffseti + j])) - Vec3(g_buffers->positions[poffseti+j]);
					Vec3 delta_rot = Vec3(0, 0, 0);// Lerp(Vec3(g_buffers->positions[poffseti + j]), p_rot, 0.8f)-Vec3(g_buffers->positions[poffseti + j]);
					g_buffers->velocities[poffseti + j].x = (delta.x + delta_rot.x) / g_dt;
					g_buffers->velocities[poffseti + j].y = (delta.y + delta_rot.y) / g_dt;
					g_buffers->velocities[poffseti + j].z = (delta.z + delta_rot.z) / g_dt;
				}
			
			}
		}
		//flexSetVelocities(g_flex, &g_buffers->velocities[0].x, g_buffers->velocities.size(), eFlexMemoryHost);
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
			FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
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
						//		float D = Length(Vec3(g_buffers->positions[poffseti + ccd_ctd[o][n]]) - Vec3(g_buffers->positions[poffsetj + ctd_ccd[p][m]]));
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

				if ((mIntegrase[i].ccd_occupied[0] == -1) ){
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
			FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			Vec3 center = Vec3(0, 0, 0);
			for (int j = 0; j < asset->mNumParticles; j++){
				center += Vec3(g_buffers->positions[poffseti + j]);
			}
			center /= asset->mNumParticles;
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
		std:string fname = "../../data/pack_data_soft";
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
		FlexTimers atimers;
		float ** p;
		float ** v;
		int ** ph;
		float ** no;
		int a;
		//distance 71-3
		//distance center->farest dimer
		//reportBinding();
		
		if (cp->comp_shape.size())
			flexExtSetRigidTransformation(cp->fcontainer, cp->mInstances.size());
		if ((g_frame % 5) == 0)
		{
			if (dojitter) jitter(!dojitter_steared);
			//cout << "free " << free_mIntegrase_CCD.size() << " " << free_mIntegrase_CTD.size() << endl;
			//if (free_mIntegrase_CCD.size() != 0) addOneRandomBinding();
		}
		checkDistanceBinding();
		WeightSpringStifness(strings_indices);
		flexSetSprings(g_flex, &g_buffers->springIndices[0], &g_buffers->springLengths[0], &g_buffers->springStiffness[0], g_buffers->springLengths.size(), eFlexMemoryHost);

		//save packing leftHand
		if ((g_frame % 1) == 0)
		{
			if (server)
				sendToClient();
			if (write_output){
				std:string pfix = "";// "_" + std::to_string(current_Exp) + "_" + std::to_string(current_run);
				writeToBinary(pfix);
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
					writeSoftTransform(current_Exp, current_run, free, true);	//#the beads
					writeRigidTransform(current_Exp, current_run, true);		//the shape rigid bodu
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
		if ((g_frame % 5) && (grow_fiber))
			cp->growFiber((int)Nsub);
	}

	void exportPDBRigidTransform(int exp, int run)
	{
		Vec3* pos;
		Quat* quat;
		int totalNbBody = 0;
		int nInst = cp->mInstances.size();
		if (cp->use_rb){
			pos = new Vec3[cp->mInstances.size()];
			quat = new Quat[cp->mInstances.size()];
			totalNbBody = nInst;
		}
		else {
			for (int i = 0; i < cp->mInstances.size(); i++){
				FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
				totalNbBody += asset->mNumShapes;
			}
			pos = new Vec3[totalNbBody];
			quat = new Quat[totalNbBody];
		}

		flexGetRigidTransforms(g_flex, (float*)&quat[0], (float*)&pos[0], eFlexMemoryHost);

		std:string name = "../../data/pack_result"+std::to_string(exp) + "_" + std::to_string(run);
		name = name + ".pdb";
		FILE *of;

		/*
		ofstream of;
		if (append){
		name = name + ".pdb";
		of.open(name.c_str(), ios::out | ios::app);
		of << "# " << cp->mInstances.size() << " " << exp << " " << run << " " << endl;
		}
		else {
		name = name + "_" + std::to_string(exp) + "_" + std::to_string(run) + ".txt";
		of.open(name.c_str(), ios::out);
		}*/

		fopen_s(&of, name.c_str(), "w");
		//write position
		for (int i = 0; i<nInst; i++){
			//cout << "  ? " << i << endl;
			//output << -pos[i].x*(1.0f/main_scale) << " " <<pos[i].y*(1.0f/main_scale)<< " " <<pos[i].z*(1.0f/main_scale)<< " " << cp->mInstances[i].mMeshIndex << "\n";
			int occup = getOccupationCCD(i);
			float bf = (mIntegrase[i].mini_distance[0] + mIntegrase[i].mini_distance[1])/2.0f;
			if (bf > 500.0f) bf = 0.0f;
			float ind = (float)cp->mInstances[i].mMeshIndex;//instance 
			FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			int nbody = asset->mNumShapes;
			int bodyindice = i*nbody;
			float p[4] = {  pos[bodyindice].x*(1.0f / main_scale)*0.05f, 
							pos[bodyindice].y*(1.0f / main_scale)*0.05f, 
							pos[bodyindice].z*(1.0f / main_scale)*0.05f, ind };
			/*
			cout << "  ? " << i << p[0] << " " << p[1] << " " << p[2] << endl;
			printf("%6s", "HETATM");
			printf("%5d", i);
			printf("%4s", "FE");
			printf("%1s", " ");
			printf("%3s", "INT");
			printf("%1s", "A");
			printf("%4d", 0);
			printf("%1s", " ");
			printf("%8.3f", p[0]);
			printf("%8.3f", p[1]);
			printf("%8.3f", p[2]);
			printf("%6.2f", 1.0f);
			printf("%6.2f", 1.0f);
			printf("%2s", "FE");
			printf("%2s\n", "");
			//{:6s}{:5d}{:^4s}{:1s} {:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}
			printf("%6s%5d%3s %1s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
				"HETATM", i, "FE", " ", "INT", "A", i, " ", p[0], p[1], p[2], 1.0, 1.0, "FE", "");
				*/
			fprintf(of, "%6s%5d%3s %1s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
				"ATOM  ", i, "C", " ", "ALA", "A", i, " ", p[0], p[1], p[2], (float)occup, bf, "C", "");
			//of << ind << " " << p[0] << " " << p[1] << " " << p[2] << " " << quat[i].x << " " << quat[i].y << " " << quat[i].z << " " << quat[i].w << endl;
			//skip the other one
			//"HETATM", 0, "FE"," ","INT","A", 0," ", coords[0][0],coords[0][1],coords[0][2],1.0,1.0, "FE","1.0"
		}
		//connectivity
		//8 -1 -1 -1 -1 9999 9999
		//7 49 9 1 0 0.110658 0.112724
		int nintegrase = mIntegrase.size();
		for (int i = 0; i < nintegrase; i++){
			if (mIntegrase[i].ccd_occupied[0] != -1 || mIntegrase[i].ccd_occupied[1] != -1){
				fprintf(of, "CONECT%5d", i);
				if (mIntegrase[i].ccd_occupied[0] != -1)
					fprintf(of, "%5d", mIntegrase[i].ccd_occupied[0]);
				if (mIntegrase[i].ccd_occupied[1] != -1)
					fprintf(of, "%5d", mIntegrase[i].ccd_occupied[1]);
				fprintf(of, "\n");
			}
		}
		fclose(of);
	}

};

class NanoCage : public HIVIntegrase
{
public:
	std::vector<std::vector<int> > binding_in;
	std::vector<std::vector<int> > binding_out;
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
	NanoCage(const char* name) : HIVIntegrase(name) {}

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

	virtual int checkOccupation(int i,int j,int b){
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
		//if (mIntegrase[i].ccd_occupied[1 - id1] == j || mIntegrase[j].ctd_occupied[1 - id2] == i) return -1;
		for (int n = 0; n < binding_in[id1].size(); n++) {
			for (int m = 0; m < binding_out[id2].size(); m++) {
				float stiff = 1.0f;
				float length = cp->mIngrSphereTree[0].DistancesMatrix[n*binding_in[id1].size()+m];
				//cout << "bound " << id1 << " " << n << " " << poffseti + ccd_ctd[id1][n] << " " << id2 << " " << m << " " << poffsetj + ctd_ccd[id2][m] << " " << g_params.mRadius << " " << g_buffers->springLengths.size() << " " g_buffers->springIndices << endl;
				strings_indices.push_back(g_buffers->springLengths.size());
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
		//flexSetSprings(g_flex, &g_buffers->springIndices[0], &g_buffers->springLengths[0], &g_buffers->springStiffness[0], g_buffers->springLengths.size(), eFlexMemoryHost);
		cout << "bound " << i << " " << j << " " << id1 << " " << id2 << " " << mNodes[i].binding_in_occupied[id1] << " " << mNodes[j].binding_out_occupied[id2] << endl;
		return 1;
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
		cout << " nb integrase " << mIntegrase.size() << endl;
	}

	virtual void setupBinding(){
		/*
		binding_in.push_back({ 58, 59, 63 });
		binding_in.push_back({ 9, 10, 33 });
		binding_in.push_back({ 1, 27, 34 });

		binding_out.push_back({ 58, 59, 63 });
		binding_out.push_back({ 1, 16, 18 });
		binding_out.push_back({ 1, 27, 34 });
		*/
		IngredientSphereTree sph = cp->mIngrSphereTree[0];
		for (int i = 0; i < sph.nBinding; i++){
			int start = sph.BindingStarts[i*sph.nBinding + 0];
			int count = sph.BindingStarts[i*sph.nBinding + 1];
			std:vector<int> lb;
			for (int j = 0; j < count; j++){
				int id = sph.BindingSites[start + j];
				lb.push_back(id);
				cout << "aadd " << id << endl;
			}
			binding_in.push_back(lb);
		}
		//bindingout is the invert
		for (int i = sph.nBinding; i >= 0; i--){
			binding_out.push_back(binding_in[i]);
		}
	}

	virtual void Initialize()
	{
		bool ignore_comp = true;
		/* validation experiement setup */
		std::string recipe = "NanoCage.json";
		std::string results = "nanoCage_results.json";
		std::string wrkDir = "../../../ProteinNanocage/";

		/*binding_in.push_back({ 33, 31, 30, 29, 28 });
		binding_in.push_back({ 0,3, 6, 7, 16, 18 });
		binding_in.push_back({ 2, 5, 15 });

		binding_out.push_back({ 33, 31, 30, 29, 28 });
		binding_out.push_back({ 2, 5, 15 });
		binding_out.push_back({ 0,3,6, 7, 16, 18 });
		[58, 59, 63]
		[58, 59, 63]
		[9, 10, 33]
		[1, 27, 34]
		[9, 10, 33]
		[1, 27, 34]

		
		come from the sphTree file. or ingredients dictionary.

		*/

		g_rigidTranslations.resize(0);
		g_rigidRotations.resize(0);
		g_rigidCoefficients.resize(0);
		g_rigidIndices.resize(0);
		g_rigidLocalPositions.resize(0);
		g_rigidOffsets.resize(0);

		// start index
		g_rigidOffsets.push_back(0);

		main_scale = 1.0f / 100.0f;

		float beads_radius = 6.5f*main_scale;
		g_params.mRadius = beads_radius * 2.0f;
		threshold_binding = g_params.mRadius*10.0f;
		output_bin.open("../../data/pack_result.bin", ios::out | ios::binary);
		output_bin.close();
		threshold_binding = g_params.mRadius*2.0f;
		//setup experiments
		exp_params[0] = { 1.0f, 1.0f, g_params.mRadius*2.0f };
		exp_params[1] = { 100.0f, 0.0f, g_params.mRadius*100.0f };//no bias toward center but higher distance and random force
		exp_params[2] = { 100.0f, 0.0f, g_params.mRadius*10000.0f };//higher distance or infinite distance ?
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
		//cp->loadRecipe((cp->mainpath + "HIV_IN_XP.1.0.json").c_str(), ignore_comp);
		//set up binding_in and the distance matrice.
		setupBinding();
		overwriteSprings();

		cp->loadResults(wrkDir + results);
		if (!ignore_comp)cp->buildMembrane(1.0f, 1.0f, NvFlexMakePhase(99999, 0));

		setupBoundsAndLinker();

		group = cp->iGroupCounter;
		g_numExtraParticles = 1024 * 1024;
		g_params.mNumPlanes = 0;
		g_params.mGravity[1] = 0.0f;//-9.f;
		g_params.mDamping = 3.0f;// 3.0f;
		g_params.mRadius = beads_radius * 2;
		//g_params.mSolidRestDistance = g_params.mRadius*2.0f;
		g_windStrength = 0.0f;
		g_windFrequency = 0.0f;
		//        g_params.mDynamicFriction = 0.00f;
		//        g_params.mFluid = true;
		//        g_params.mViscosity = 0.0f;
		g_params.mFluid = false;
		g_params.mVorticityConfinement = 0.0f;
		//g_params.mFluidRestDistance = g_params.mRadius*0.5;
		g_params.mAnisotropyScale = 2.5f / g_params.mRadius;
		g_params.mSmoothing = 0.5f;
		g_params.mRelaxationFactor = 1.f;
		g_params.mRestitution = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		g_params.mDynamicFriction = 0.25f;
		g_params.mViscosity = 1.5f;
		g_params.mCohesion = 1.1f;
		g_params.mAdhesion = 1.0f;
		g_params.mSurfaceTension = 1.0f;

		g_params.mGravity[1] = 0.0f;

		g_lightDistance *= 2.0f;

		// draw options		
		g_drawEllipsoids = true;
		g_params.mMaxSpeed = 100.0f;

		g_wireframe = true;
		g_pointScale = 1.0f;
		//        g_blur = 2.0f;
		g_pause = true;
		g_warmup = false;
		g_drawMesh = true;
		g_drawRopes = false;
		g_drawPoints = true;
		g_params.mDynamicFriction = 0.4f;

		g_params.mDissipation = 0.0f;
		g_params.mNumIterations = 4;
		g_params.mViscosity = 0.0f;
		g_params.mDrag = 0.0f;
		g_params.mLift = 0.0f;

		g_params.mParticleCollisionMargin = g_params.mRadius*0.05f;
		g_params.mDrag = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		// better convergence with global relaxation factor
		g_params.mRelaxationMode = eFlexRelaxationGlobal;
		g_params.mRelaxationFactor = 0.25f;

		g_windStrength = 0.0f;

		g_numSubsteps = 5;
		g_params.mNumIterations = 5;

		mSplitThreshold.resize(mCloths.size(), 45.0f);

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
			FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
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
		Vec3* pos;
		Quat* quat;
		int totalNbBody = 0;
		int nInst = cp->mInstances.size();
		if (cp->use_rb){
			pos = new Vec3[cp->mInstances.size()];
			quat = new Quat[cp->mInstances.size()];
			totalNbBody = nInst;
		}
		else {
			for (int i = 0; i < cp->mInstances.size(); i++){
				FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
				totalNbBody += asset->mNumShapes;
			}
			pos = new Vec3[totalNbBody];
			quat = new Quat[totalNbBody];
		}

		flexGetRigidTransforms(g_flex, (float*)&quat[0], (float*)&pos[0], eFlexMemoryHost);

		std:string name = "../../data/pack_result" + std::to_string(exp) + "_" + std::to_string(run);
		name = name + ".pdb";
		FILE *of;

		/*
		ofstream of;
		if (append){
		name = name + ".pdb";
		of.open(name.c_str(), ios::out | ios::app);
		of << "# " << cp->mInstances.size() << " " << exp << " " << run << " " << endl;
		}
		else {
		name = name + "_" + std::to_string(exp) + "_" + std::to_string(run) + ".txt";
		of.open(name.c_str(), ios::out);
		}*/

		fopen_s(&of, name.c_str(), "w");
		//write position
		for (int i = 0; i<nInst; i++){
			//cout << "  ? " << i << endl;
			//output << -pos[i].x*(1.0f/main_scale) << " " <<pos[i].y*(1.0f/main_scale)<< " " <<pos[i].z*(1.0f/main_scale)<< " " << cp->mInstances[i].mMeshIndex << "\n";
			//int occup = getOccupationCCD(i);
			float bf = (mNodes[i].mini_distance[0] + mNodes[i].mini_distance[1]) / 2.0f;
			if (bf > 500.0f) bf = 0.0f;
			float ind = (float)cp->mInstances[i].mMeshIndex;//instance 
			FlexExtAsset* asset = cp->iBatches[cp->mInstances[i].mMeshIndex].mAsset;
			int nbody = asset->mNumShapes;
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

	void Update()
	{
		FlexTimers atimers;
		float ** p;
		float ** v;
		int ** ph;
		float ** no;
		int a;
		//distance 71-3
		//distance center->farest dimer
		//reportBinding();

		if (cp->comp_shape.size())
			flexExtSetRigidTransformation(cp->fcontainer, cp->mInstances.size());
		if ((g_frame % 5) == 0)
		{
			if (dojitter) jitter(!dojitter_steared);
			//cout << "free " << free_mIntegrase_CCD.size() << " " << free_mIntegrase_CTD.size() << endl;
			//if (free_mIntegrase_CCD.size() != 0) addOneRandomBinding();
		}
		checkDistanceBinding();
		//WeightSpringStifness(strings_indices);
		flexSetSprings(g_flex, &g_buffers->springIndices[0], &g_buffers->springLengths[0], &g_buffers->springStiffness[0], g_buffers->springLengths.size(), eFlexMemoryHost);

		//save packing leftHand
		if ((g_frame % 1) == 0)
		{
			if (server)
				sendToClient();
			if (write_output){
			std:string pfix = "";// "_" + std::to_string(current_Exp) + "_" + std::to_string(current_run);
				writeToBinary(pfix);
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
					writeSoftTransform(current_Exp, current_run, free, true);	//#the beads
					writeRigidTransform(current_Exp, current_run, true);		//the shape rigid bodu
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
		if ((g_frame % 5) && (grow_fiber))
			cp->growFiber((int)Nsub);
	}
};


class ActineBranch : public Mycoplasma
{
public:

	ActineBranch(const char* name) : Mycoplasma(name) {}
	struct Branch
	{
		int npoints;
		std::vector<int> p_indices;
	};
	std::vector<Branch> all_branch;
	std::vector<Vec3> all_points;
	//float main_scale = 1.0f / 1000.0f;
	int iGroupCounter=0;
	Graph H;
	Graph G;
	GraphAttributes HA;
	GraphAttributes GA;

	virtual void MakeAdamGreatAgain(){
		for (int i = 0; i < all_branch.size(); i++){
			//data_curve is an array of float
			std::vector<Vec3> points;
			for (int j = 0; j < all_branch[i].npoints; j++){
				points.push_back(all_points[all_branch[i].p_indices[j]]);
			}
			float* data_curve = reinterpret_cast<float*>(points.data());
			int rope_phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);
			Rope curve;

			CreateRopeFromData(curve, //rope
				main_scale, // scale
				1.0f, //stifness
				data_curve, //data
				0.0f, //length
				all_branch[i].npoints*3,  //nfloat
				rope_phase,//phase
				0.0f,//spiral angle
				1.0f,//invmass
				0.0f,//give
				0,//extend_nb
				false,//extend
				false,//close
				g_params.mRadius*2.0f,
				0);
			curve.persistence = 0;
			printf("create Rope from data OK with %i points\n", curve.mIndices.size());
			g_ropes.push_back(curve);//instance
		}
	}

	virtual void BuildNetwork(){
		int begin = 0;
		int persistence = 5;
		float D = (7 / 1000.0f)/2.0f;
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
		int branch_i=-1;
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
	bool isReachable(node s, node d, Graph G)
	{
		// Base case
		if (s == d)
			return true;

		// Mark all the vertices as not visited
		NodeArray<bool> visited(G, false);
		// Create a queue for BFS
		SListPure<node> bfs_queue;
		bfs_queue.pushBack(s);
		visited[s] = true;

		// it will be used to get all adjacent vertices of a vertex
		//list<int>::iterator i;
		node w;
		node adj;
		edge e;

		while (!bfs_queue.empty())
		{

			// Dequeue a vertex from queue and print it
			s = bfs_queue.popFrontRet();

			// Get all adjacent vertices of the dequeued vertex s
			// If a adjacent has not been visited, then mark it visited
			// and enqueue it
			forall_adj_edges(e, w)
			//for (i = adj[s].begin(); i != adj[s].end(); ++i)
			{
				adj = e->opposite(w);
				// If this adjacent node is the destination node, then return true
				if (adj == d)
					return true;
				if (adj->degree()>2){
					if (!visited[adj])
					{
						visited[adj] = true;
					}
					continue;
				}
				// Else, continue to do BFS
				if (!visited[adj])
				{
					visited[adj] = true;
					bfs_queue.pushBack(adj);
				}
			}
		}
		return false;
	}

	List<node> pathFromTo(node s, node d, Graph G, NodeArray<bool> visited)
	{
		//std::cout << " pathFromTo " << endl;
		List<node> bfs_queue;
		bfs_queue.pushBack(s);
		// Base case
		if (s == d){
			//std::cout << "already found ? " << endl;
			bfs_queue.pushBack(d);
			return bfs_queue;
		}
		// Mark all the vertices as not visited
		// NodeArray<bool> visited(G, false);
		// Create a queue for BFS
		
		visited[s] = true;

		// it will be used to get all adjacent vertices of a vertex
		//list<int>::iterator i;
		node w;
		node adj;
		edge e;

		//std::cout << " ? " << endl;

		std:vector<List<node>> res;
		forall_adj_edges(e, s){
			adj = e->opposite(s);
			if (adj == d) {
				//::cout << "already found  " << endl;
				bfs_queue.pushBack(d);
				return bfs_queue;
			}
			List<node> queue;
			//queue.pushBack(s);
			//if (adj->degree()>25)
			//{
				//two points segments ?
				//remove the branch ?
				//std::cout << " two points segments" << adj << " " << adj->index() << " " << adj->degree() << endl;
				//if (adj == d) {
				//	std::cout << "already found ? " << endl;
				//	bfs_queue.pushBack(d);
				//	return bfs_queue;
				//}
				//else 
				//	continue;//should continue the for loop

			//}
			//else {
				queue.pushBack(adj);
				res.push_back(queue);
			//}
			//std::cout << "add to explore branch starting at " << e << " " << adj->index() << " " << GA.label(adj) << endl;
		}
		bool found=false;
		int found_branch=-1;
		for (int i = 0; i < res.size(); i++){
			if (found) break;
			SListPure<node> aqueue;
			aqueue.pushBack(res[i].front());
			//cout << " test queue " << i << " " <<  endl;
			while (!aqueue.empty())
			{

				// Dequeue a vertex from queue and print it
				w = aqueue.popFrontRet();
				if (w == d){
					res[i].pushFront(s);
					res[i].pushBack(w);
					found = true;
					found_branch = i;
					//std::cout << " found " << i << " " << adj->index() << " " << res[i].size() << endl;
					aqueue.clear();
					break;
				}
				//cout << "test " << w << endl;
				//if (visited[w])
				//{
				//	cout << "visited " << w << endl;
				//	continue;
				//}
				// Get all adjacent vertices of the dequeued vertex s
				// If a adjacent has not been visited, then mark it visited
				// and enqueue it
				forall_adj_edges(e, w)
				//for (i = adj[s].begin(); i != adj[s].end(); ++i)
				{
					//cout << "edge " << e << endl;
					adj = e->opposite(w);//opposite? could it be start ?
					//cout << "Test " << adj << endl;
					if (visited[adj]) continue;
					if (adj == s) {
						//special case
						visited[w] = true;
						res[i].pushBack(w);
						//aqueue.clear();
						//cout << "back to start " << w << " " << s << " " << e->source() << " " << e->target() << " " << aqueue.empty() << endl;
						continue;//what does it break
					}
					// If this adjacent node is the destination node, then return true
					else if (adj == d){
						res[i].pushFront(s);
						res[i].pushBack(adj);
						found = true;
						found_branch = i;
						//std::cout << " found " << i << " " << adj->index() << " " << res[i].size() << endl;
						//aqueue.clear();
						break;
					}
					else if (adj->degree()>2){
						
						if (!visited[adj])
						{
							visited[adj] = true;
						}
						//empty the list
						//aqueue.clear();
						//remove the branch ?
						//std::cout << " end with other node " << i << " " << adj << " " << adj->degree() << " " << aqueue.empty() <<  endl;
						continue;//should continue the for loop
					}
					// Else, continue to do BFS
					else if (!visited[adj])
					{
						//std::cout << "visit " << adj << endl;
						visited[adj] = true;
						aqueue.pushBack(adj);//or pushFront?
						res[i].pushBack(adj);
					}
					else {
					//already visit
						//std::cout << " else ??  " << adj << endl;

					}
				}
			}
		}
		if (found) {
			return res[found_branch];
		}
		else{
			std::cout << "not found" << s << " " << d << endl;
			return bfs_queue;
		}
	}

	std::vector<std::map<node, node>> _bidirectional_pred_succ(Graph G, node source, node target, node& f)
	{
		//"""Bidirectional shortest path helper.
		//
		//Returns(pred, succ, w) where
		//pred is a dictionary of predecessors from w to the source, and
		//succ is a dictionary of successors from w to the target.
		//"""
		//# does BFS from both source and target and meets in the middle
		std::map<node, node> mymapS;
		mymapS[source] = source;
		std::map<node, node> mymapT;
		mymapS[target] = target;
		std:vector<std::map<node, node>> results;
		
		if (target == source){
			results.push_back(mymapS);
			results.push_back(mymapT);
			return results;
		}
		//# handle either directed or undirected
		//Gpred = G.neighbors_iter
		//Gsucc = G.neighbors_iter

		//# predecesssor and successors in search
		std::map<node, node> pred; 
		pred[source] = NULL;
		std::map<node, node> succ; 
		succ[target] = NULL;

		//# initialize fringes, start with forward
		SListPure<node> forward_fringe; forward_fringe.pushBack(source);
		SListPure<node> reverse_fringe; reverse_fringe.pushBack(target);
		
		edge e;
		node v;
		node w;

		while ((!forward_fringe.empty()) && (!reverse_fringe.empty()))
		{
			if (forward_fringe.size() < reverse_fringe.size()){
				SListPure<node> this_level = forward_fringe;
				forward_fringe.clear();
				for (auto v : this_level){
					forall_adj_edges(e, v){
						w = e->opposite(v);
						if (pred.find(w) == pred.end()) {
							forward_fringe.pushBack(w);
							pred[w] = v;
						}
						if (succ.find(w) != succ.end()){
							f = w;
							results.push_back(pred);
							results.push_back(succ);
							return results;
						}
					}
				}
			}
			else {
				SListPure<node> this_level = reverse_fringe;
				reverse_fringe.clear();
				for (auto v : this_level){
					forall_adj_edges(e, v){
						w = e->opposite(v);
						if (succ.find(w) == succ.end()) {
							succ[w] = v;
							reverse_fringe.pushBack(w);
						}
						if (pred.find(w) != pred.end()){
							f = w;
							results.push_back(pred);
							results.push_back(succ);
							return results;
						}
					}
				}
			}
		}
		std::cout << "No path between %s and %s." << endl;
		results.push_back(mymapS);
		results.push_back(mymapT);
		return results;
	}

	List<node> bidirectional_shortest_path(Graph G, node source, node target)
	{
		node w;
		std::vector<std::map<node, node>> results = _bidirectional_pred_succ(G, source, target, w);
		node start = w;
		//std::cout << " found " << w << " " << start <<  endl;

		//# build path from pred + w + succ
		List<node> path;
		
		//# from source to w
		while (results[0].find(w) != results[0].end())
		{
			path.pushBack(w);
			w = results[0][w];//pred
		}
		path.reverse();
		//# from w to target
		w = results[1][start];
		while (results[1].find(w) != results[1].end())
		{
			path.pushBack(w);
			w = results[1][w];//succ
		}
		return path;
	}

	virtual void BuildParticleAndNetworkFromSimplifiedGraph(){
		//std::set<int> used_indices;
		ropes_mIndices.clear();
		std::vector<int> used_indices;
		std::vector<int> used_indices_gpart;
		std::set<int>::iterator it;

		int begin = 0;
		int persistence = 1;
		float D = g_params.mRadius;//distance between two points for collision
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
			int phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);
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
				
				if ((j > 0) && (dist < g_params.mRadius)){
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
			CreatePersistenceRope(rope, current, persistence, stiffness, 0, give, D/2.0f, hardness);

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

	virtual void BuildParticleAndNetworkFromGraph(){
		//std::set<int> used_indices;
		std::vector<int> used_indices;
		std::vector<int> used_indices_gpart;
		std::set<int>::iterator it;
		int begin = 0;
		int persistence = 1;
		float D = g_params.mRadius;//distance between two points for collision
		int start = int(g_buffers->positions.size());
		float give = 0.0f;
		float stiffness = 1.0f;
		float r = 1.0f;//biased on the 1-3 spring
		int current = 0;
		//if closed do the last point ?
		int count = 0;
		int prev = 0;
		float hardness = 0.0f; //lead to D/1000.0f
		//#segments = []
		//#search = set()
		//#for ed in H.edge :
		//#    edgs = H.edge[ed].keys()
		//#    for i in range(len(edgs)) :
		//#        if str(int(edgs[i]))+"_"+str(int(ed)) in search :
		//#            continue
		//#        li = nx.shortest_path(G, int(ed), int(edgs[i]))#give back all the point for this segment
		//#        new_points = interpolate(li, xyz, 27.6)
		//#        segments.append(new_points)
		//#        search.add(str(int(ed)) + "_" + str(int(edgs[i])))

		//go through all edge of the graph H
		List<edge> edgs;//one edge is two nodes, source and target
		H.allEdges(edgs);
		List<node> nodes;
		G.allNodes(nodes);
		NodeArray<bool> visited(G, false);
		EdgeArray<bool> edge_visited(H, false);
		int countedge = 0;
		//for (int i = 0; i < edgs.size(); i++){//edgs.size()
		for (auto e : edgs){
			countedge++;
			//edge e = (*edgs.get(i));
			if (countedge > 50) break;
			if (edge_visited[e]) continue;
			node a = e->source();
			node b = e->target();
			//need the node from G
			int ai = atoi(HA.label(a).c_str());
			int bi = atoi(HA.label(b).c_str());
			//get the node in G
			node ga = (*nodes.get(ai));
			node gb = (*nodes.get(bi));
			cout << "G" << " " << ga << " " << gb << endl;
			//cout << a->index() << " " << ai << " " << b->index() << " " << bi << endl;
			List<node> test = pathFromTo(ga, gb, G, visited);
			if (test.size() == 1) break;
			//List<node> test = bidirectional_shortest_path(G, ga, gb);
			//for (int i = 0; i < test.size(); i++){
			//	node n = (*test.get(i));
			//		int labelid = atoi(GA.label(n).c_str());
			//	cout << n->index() << " ";
			//}
			//that is a branch
			int phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);
			Rope rope;
			node n = *test.begin();
			Vec3 start = all_points[atoi(GA.label(n).c_str())];
			if (test.size() == 2){
				continue;
				int pindice = atoi(GA.label(ga).c_str()); //indice in all_pos
				int gindice = int(g_buffers->positions.size());    //indice in particle
				std::vector<int>::iterator iter = std::find(used_indices.begin(), used_indices.end(), pindice); //used_indices.find(pindice);		  //already used ?
				if (iter != used_indices.end())
				{
					//indice found can reuse it
					int k = std::distance(used_indices.begin(), iter);
					gindice = used_indices_gpart.at(k);
					count++;
				}
				else
				{
					used_indices.push_back(pindice);
					used_indices_gpart.push_back(gindice);
					g_buffers->positions.push_back(Vec4(all_points[pindice].x, all_points[pindice].y, all_points[pindice].z, 1.0f));
					g_buffers->velocities.push_back(0.0f);
					g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
				}
				pindice = atoi(GA.label(gb).c_str()); //indice in all_pos
				int gindice2 = int(g_buffers->positions.size());    //indice in particle
				iter = std::find(used_indices.begin(), used_indices.end(), pindice); //used_indices.find(pindice);		  //already used ?
				if (iter != used_indices.end())
				{
					//indice found can reuse it
					int k = std::distance(used_indices.begin(), iter);
					gindice2 = used_indices_gpart.at(k);
					count++;
				}
				else
				{
					used_indices.push_back(pindice);
					used_indices_gpart.push_back(gindice2);
					g_buffers->positions.push_back(Vec4(all_points[pindice].x, all_points[pindice].y, all_points[pindice].z, 1.0f));
					g_buffers->velocities.push_back(0.0f);
					g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
				}
				CreateSpringInter(gindice, gindice2, stiffness, give, D / 2.0f);
			}
			else {
				//generate one segments
				for (int j = 0; j < test.size()-1; j++){
					n = (*test.get(j));
					int pindice = atoi(GA.label(n).c_str()); //indice in all_pos
					int gindice = int(g_buffers->positions.size());    //indice in particle
					float dist = Length(all_points[pindice] - start);
					if ((j > 0) && (dist < g_params.mRadius)){
						continue;
					}
					else {
						start = all_points[pindice];
					}
					std::vector<int>::iterator iter = std::find(used_indices.begin(), used_indices.end(), pindice); //used_indices.find(pindice);		  //already used ?
					if (iter != used_indices.end())
					{
						//indice found can reuse it
						int k = std::distance(used_indices.begin(), iter);
						gindice = used_indices_gpart.at(k);
						count++;
					}
					else
					{
						used_indices.push_back(pindice);
						used_indices_gpart.push_back(gindice);
						g_buffers->positions.push_back(Vec4(all_points[pindice].x, all_points[pindice].y, all_points[pindice].z, 1.0f));
						g_buffers->velocities.push_back(0.0f);
						g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
					}
					current = int(rope.mIndices.size());
					rope.mIndices.push_back(gindice);
					CreatePersistenceRope(rope, current, persistence, stiffness, 0, give, D / 2.0f, hardness);
				}//end for loop

				//add the last element
				n = (*test.get(test.size() - 1));
				int pindice = atoi(GA.label(n).c_str()); //indice in all_pos
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
					g_buffers->positions.push_back(Vec4(all_points[pindice].x, all_points[pindice].y, all_points[pindice].z, 1.0f));
					g_buffers->velocities.push_back(0.0f);
					g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
				}
				current = int(rope.mIndices.size());
				rope.mIndices.push_back(gindice);
				//CreatePersistenceRope(rope, current, persistence, stiffness, 0, give, D/2.0f, hardness);

				float dist = Length(Vec3(g_buffers->positions[rope.mIndices[0]]) - Vec3(g_buffers->positions[rope.mIndices[rope.mIndices.size() - 1]]));
				float a = (D / 2.0f)*(float)rope.mIndices.size();// all_branch[i].npoints);
				float b = a;
				if (dist > a) {
					b = dist;
				}
				//CreateSpringInter(rope.mIndices[0], rope.mIndices[rope.mIndices.size() - 1], stiffness, give, b);
				rope.persistence = persistence;
				g_ropes.push_back(rope);//instance
			}//end else
			cout << test.size() << " " << countedge << endl;
			edge_visited[e] = true;
			
		}
	}

	virtual void BuildParticleAndNetwork(){
		//std::set<int> used_indices;
		std::vector<int> used_indices;
		std::vector<int> used_indices_gpart;
		std::set<int>::iterator it;
			
		int begin = 0;
		int persistence = 1;
		float D = g_params.mRadius;//distance between two points for collision
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
			int phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);
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

			Vec3 start = all_points[all_branch[i].p_indices[0]];
			for (int j = 0; j < all_branch[i].npoints-1; j ++){// all_branch[i].npoints; j++){
				int pindice = all_branch[i].p_indices[j]; //indice in all_pos
				int gindice = int(g_buffers->positions.size());    //indice in particle
				float dist = Length(all_points[pindice] - start); 
				/*
				if ((j > 0) && (dist < g_params.mRadius)){
					int branchid = checkIfExistInOtherBranch(pindice, i);
					if (branchid == -1) { continue; }//store indice ?
					else {
						start = all_points[pindice];
					}
				}
				else {
					start = all_points[pindice];
				}
				*/
				
				std::vector<int>::iterator iter = std::find(used_indices.begin(), used_indices.end(),  pindice); //used_indices.find(pindice);		  //already used ?
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
					g_buffers->positions.push_back(Vec4(all_points[pindice].x, all_points[pindice].y, all_points[pindice].z, 1.0f));
					g_buffers->velocities.push_back(0.0f);
					g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
				}

				current = int(rope.mIndices.size());
				rope.mIndices.push_back(gindice);
				CreatePersistenceRope(rope, current, persistence, stiffness, 0, give, D/2.0f,hardness);
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
				g_buffers->positions.push_back(Vec4(all_points[pindice].x, all_points[pindice].y, all_points[pindice].z, 1.0f));
				g_buffers->velocities.push_back(0.0f);
				g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
			}
			current = int(rope.mIndices.size());
			//rope.mIndices.push_back(gindice);
			//CreatePersistenceRope(rope, current, persistence, stiffness, 0, give, D/2.0f, hardness);
			
			float dist = Length(Vec3(g_buffers->positions[rope.mIndices[0]]) - Vec3(g_buffers->positions[rope.mIndices[rope.mIndices.size()-1]]));
			float a = (D / 2.0f)*(float)rope.mIndices.size();// all_branch[i].npoints);
			float b = a;
			if (dist > a) {
				b = dist;
			}
			//CreateSpringInter(rope.mIndices[0], rope.mIndices[rope.mIndices.size() - 1], stiffness, give, b);
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
		string datapath = "..\\..\\data\\";
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
			output_bin << g_buffers->positions[i].x << " " << g_buffers->positions[i].y << " " << g_buffers->positions[i].z << " " << g_params.mRadius / 2.0f << endl;
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
		output_bin.open((filename+"_graph.txt").c_str(), ios::out );

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
			output_bin << g_buffers->positions[i].x << " " << g_buffers->positions[i].y << " " << g_buffers->positions[i].z << " " << g_params.mRadius/2.0f << endl;
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
			if (i==0)std::cout << "Branch pt 0 " << r.p_indices[0] << endl;
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
		g_rigidTranslations.resize(0);
		g_rigidRotations.resize(0);
		g_rigidCoefficients.resize(0);
		g_rigidIndices.resize(0);
		g_rigidLocalPositions.resize(0);
		g_rigidOffsets.resize(0);

		// start index
		g_rigidOffsets.push_back(0);
		//unit in file 1 unit=1angstrom?? or 1 unit = 10nm
		main_scale = 1.0f / 100.0f;
		//actin filaments (13.5 μm) persistance length 13000/7
		float actine_radius = 27.6f/2.0f;// 7.0f; thast in nm rise is 27.6A
		float beads_radius = actine_radius*main_scale;//5.0f*main_scale;//0.08f;//5.0f*main_scale;//smallest sphere size
		g_params.mRadius = beads_radius * 2.0f;

		cp = new cellPACK();
		cp->use_rb = true;
		cp->main_radius = beads_radius;//((10.0f/100.0f)*(10.0f/100.0f));
		cp->main_scale = main_scale;
		cp->mask.push_back(-1);
		cp->maks_protein.push_back(-(0 + 1));//it ptype is 0 ?
		cp->maks_fiber.push_back(0);//proteinType

		//filement mesh
		/*string geompath = "C:\\Dev\\flexpack_dev_1.0\\data\\";
		cout << " read mesh 1 " << geompath + "filament.obj" << endl;
		Mesh* mesh1 = GetMesh(GetFilePathByPlatform((geompath + "filament.obj").c_str()).c_str(), main_scale);
		Quat q = QuatFromAxisAngle(Vec3(1, 0, 0), 3.14f / 2.0f);
		mesh1->Transform(RotationMatrix(q));
		mesh1->CalculateNormals();
		cout << " ok mesh 1 " << geompath + "filament.obj" << endl;

		FlexTriangleMesh* mesh = CreateTriangleMesh(mesh1);
		AddTriangleMesh(mesh, Vec3(), Quat(), 1.0f);
		*/
		//read in the point and connection for all fiber
		//string geompath = "C:\\Dev\\flexpack_dev_1.0\\data\\";
		/*cout << " read mesh 1 " << geompath + "membrane_red.obj" << endl;
		Mesh* mesh0 = GetMesh(GetFilePathByPlatform((geompath + "membrane_red.obj").c_str()).c_str(), main_scale);
		//Quat q = QuatFromAxisAngle(Vec3(1, 0, 0), 3.14f / 2.0f);
		mesh0->Transform(RotationMatrix(q));
		mesh0->CalculateNormals();
		cout << " ok mesh 1 " << geompath + "membrane_red.obj" << endl;

		FlexTriangleMesh* mesh00 = CreateTriangleMesh(mesh0);
		AddTriangleMesh(mesh00, Vec3(), Quat(), 1.0f);
		*/
		/*

		Mesh* mesh2 = GetMesh(GetFilePathByPlatform((geompath + "erf_red.obj").c_str()).c_str(), main_scale);
		mesh2->Transform(RotationMatrix(q));
		mesh2->CalculateNormals();
		FlexTriangleMesh* meshc1 = CreateTriangleMesh(mesh2);
		AddTriangleMesh(meshc1, Vec3(), Quat(), 1.0f);
		
		Mesh* mesh3 = GetMesh(GetFilePathByPlatform((geompath + "psd_red.obj").c_str()).c_str(), main_scale);
		mesh3->Transform(RotationMatrix(q));
		mesh3->CalculateNormals();
		FlexTriangleMesh* meshc2 = CreateTriangleMesh(mesh3);
		AddTriangleMesh(meshc2, Vec3(), Quat(), 1.0f);
		*/

		//check if file exist
		std::string datapath = "..\\..\\data\\";
		std::ifstream source(datapath + "cache.bin", std::ios::binary);
		bool force = false;
		if ((source.is_open())&&(!force)) {
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
			//read the grap

			/*
			string filename = datapath + "graph_node.gml";
			HA = GraphAttributes(H, HA.nodeLabel|HA.nodeId);
			GraphIO::readGML(HA, H, filename);

			string filenameG = datapath + "graph_node_all.gml";
			GA = GraphAttributes(G, GA.nodeLabel | GA.nodeId);
			GraphIO::readGML(GA, G, filenameG);
			*/
			//GraphIO::readGML(&G, datapath + "graph_node_all.gml")
			BuildParticleAndNetworkFromSimplifiedGraph();
			//BuildParticleAndNetwork();
			writeCache(datapath + "cache.bin");
		}
		//load fomr cache if exist
		
		//Mesh* level = ImportMesh(GetFilePathByPlatform("../../data/testzone.obj").c_str());
		//level->Normalize(100.0f);
		//level->Transform(TranslationMatrix(Point3(0.0f, -5.0f, 0.0f)));
		//level->CalculateNormals();

		//Vec3 lower, upper;
		//level->GetBounds(lower, upper);
		//Vec3 center = (lower + upper)*0.5f;


		
		g_numExtraParticles = 1024 * 1024;
		g_params.mNumPlanes = 0;
		
		g_params.mGravity[1] = 0.0f;//-9.f;
		g_params.mDamping = 1.0f;// 3.0f;
		g_params.mRadius = beads_radius * 2;
		
		g_windStrength = 0.0f;
		g_windFrequency = 0.0f;
		
		g_params.mNumIterations = 5;
		
		g_params.mFluid = false;
		g_params.mVorticityConfinement = 0.0f;
		//g_params.mFluidRestDistance = g_params.mRadius*0.5;
		g_params.mAnisotropyScale = 2.5f / g_params.mRadius;
		g_params.mSmoothing = 0.5f;
		g_params.mRelaxationFactor = 1.f;
		g_params.mRestitution = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		g_params.mDynamicFriction = 0.25f;
		g_params.mViscosity = 1.5f;
		g_params.mCohesion = 1.1f;
		g_params.mAdhesion = 1.0f;
		g_params.mSurfaceTension = 1.0f;

		g_params.mGravity[1] = 0.0f;


		g_lightDistance *= 2.0f;

		// draw options		
		g_drawEllipsoids = true;
		g_params.mMaxSpeed = 100.0f;

		
		g_wireframe = true;
		g_pointScale = 1.0f;
		//        g_blur = 2.0f;
		g_pause = true;
		g_warmup = false;
		g_drawMesh = false;
		g_drawRopes = false;

		g_params.mDynamicFriction = 0.4f;
		g_params.mDissipation = 0.0f;

		g_params.mParticleCollisionMargin = g_params.mRadius*0.05f;
		g_params.mDrag = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		// better convergence with global relaxation factor
		g_params.mRelaxationMode = eFlexRelaxationGlobal;
		g_params.mRelaxationFactor = 0.25f;

		g_windStrength = 0.0f;

		g_numSubsteps = 5;
		g_params.mNumIterations = 5;

		// draw options		
		g_drawPoints = true;
		g_drawSprings = false;
		g_drawCloth = false;
	}
};

class MeshTest : public Scene
{
public:

	MeshTest(const char* name) : Scene(name) {}
	float main_scale;
	float mPressure;
	float radius;
	std::vector<ClothMesh*> mCloths;
	std::vector<float> mRestVolume;
	std::vector<int> mTriOffset;
	std::vector<int> mTriCount;
	std::vector<float> mOverPressure;
	std::vector<float> mConstraintScale;
	std::vector<float> mSplitThreshold;
	cellPACK * cp;

	void AddInflatable(const Mesh* mesh, float overPressure, int phase, float L=0.0f)
	{
		const int startVertex = g_buffers->positions.size();

		// add mesh to system
		for (size_t i = 0; i < mesh->GetNumVertices(); ++i)
		{
			const Vec3 p = Vec3(mesh->m_positions[i]);

			g_buffers->positions.push_back(Vec4(p.x, p.y, p.z, 1.0f));
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(phase);
		}

		int triOffset = g_triangles.size();
		int triCount = mesh->GetNumFaces();

		mTriOffset.push_back(triOffset / 3);
		mTriCount.push_back(mesh->GetNumFaces());
		mOverPressure.push_back(overPressure);

		for (size_t i = 0; i < mesh->m_indices.size(); i += 3)
		{
			int a = mesh->m_indices[i + 0];
			int b = mesh->m_indices[i + 1];
			int c = mesh->m_indices[i + 2];

			Vec3 n = -Normalize(Cross(mesh->m_positions[b] - mesh->m_positions[a], mesh->m_positions[c] - mesh->m_positions[a]));
			g_triangleNormals.push_back(n);

			g_triangles.push_back(a + startVertex);
			g_triangles.push_back(b + startVertex);
			g_triangles.push_back(c + startVertex);
		}

		// create a cloth mesh using the global positions / indices
		ClothMesh* cloth = new ClothMesh(&g_buffers->positions[0], g_buffers->positions.size(), &g_triangles[triOffset], triCount * 3, 0.8f, 1.0f);

		for (size_t i = 0; i < cloth->mConstraintIndices.size(); ++i)
			g_buffers->springIndices.push_back(cloth->mConstraintIndices[i]);

		g_buffers->springStiffness.insert(g_buffers->springStiffness.end(), cloth->mConstraintCoefficients.begin(), cloth->mConstraintCoefficients.end());
		//change Rest Length
		if (L != 0.0f){
			for (int i = 0; i < cloth->mConstraintRestLengths.size(); i++){
				cloth->mConstraintRestLengths[i] = L;
			}
		}
		g_buffers->springLengths.insert(g_buffers->springLengths.end(), cloth->mConstraintRestLengths.begin(), cloth->mConstraintRestLengths.end());

		mCloths.push_back(cloth);

		// add inflatable params
		mRestVolume.push_back(cloth->mRestVolume);
		mConstraintScale.push_back(cloth->mConstraintScale);
	}

	void buildMassSpringQuad(Mesh* mesh, float stiffness, float mass, int phase, float L=0.0f)
	{
		//AddInflatable(cp->comp_mesh[0], 1.0f, cp->iGroupCounter++);
		// add particles to system
		int ph = NvFlexMakePhase(phase, eFlexPhaseSelfCollide);
		int start = g_buffers->positions.size();
		
		for (size_t i = 0; i < mesh->GetNumVertices(); ++i)
		{
			const Vec3 p = Vec3(mesh->m_positions[i]);// +Vec3(mesh->m_normals[i])*11.0f;

			g_buffers->positions.push_back(Vec4(p.x, p.y, p.z, 1.00f / mass));
			g_buffers->restPositions.push_back(Vec4(p.x, p.y, p.z, 0.0));
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(ph);
		}

		for (size_t i = 0; i < mesh->m_indices.size(); i += 4)
		{
			int a = mesh->m_indices[i + 0];
			int b = mesh->m_indices[i + 1];
			int c = mesh->m_indices[i + 2];
			int d = mesh->m_indices[i + 3];
			float La = Length(Vec3(g_buffers->positions[start + a]) - Vec3(g_buffers->positions[start + b]));
			float Lb = Length(Vec3(g_buffers->positions[start + b]) - Vec3(g_buffers->positions[start + c]));
			float Lc = Length(Vec3(g_buffers->positions[start + c]) - Vec3(g_buffers->positions[start + d]));
			float Ld = Length(Vec3(g_buffers->positions[start + d]) - Vec3(g_buffers->positions[start + a]));

			CreateSpringInter(start + a, start + b, stiffness, 0.0f, (L == 0.0f) ? La : L);// main_radius*2.0f);
			CreateSpringInter(start + b, start + c, stiffness, 0.0f, (L == 0.0f) ? Lb : L);
			CreateSpringInter(start + c, start + d, stiffness, 0.0f, (L == 0.0f) ? Lc : L);
			CreateSpringInter(start + d, start + a, stiffness, 0.0f, (L == 0.0f) ? Ld : L);
		}
	}

	void buildMassSpringTri(Mesh* mesh, float stiffness, float mass, int phase, float L = 0.0f)
	{
		//AddInflatable(cp->comp_mesh[0], 1.0f, cp->iGroupCounter++);
		// add particles to system

		int ph = NvFlexMakePhase(phase, eFlexPhaseSelfCollide);
		int start = g_buffers->positions.size();
		cp->mask_membrane.push_back(start);

		for (size_t i = 0; i < mesh->GetNumVertices(); ++i)
		{
			const Vec3 p = Vec3(mesh->m_positions[i]) + Vec3(mesh->m_normals[i])*g_params.mRadius;

			g_buffers->positions.push_back(Vec4(p.x, p.y, p.z, 1.00f / mass));
			g_buffers->restPositions.push_back(Vec4(p.x, p.y, p.z, 0.0));
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(ph);
		}

		for (size_t i = 0; i < mesh->m_indices.size(); i += 3)
		{
			int a = mesh->m_indices[i + 0];
			int b = mesh->m_indices[i + 1];
			int c = mesh->m_indices[i + 2];
			
			float La = Length(Vec3(g_buffers->positions[start + a]) - Vec3(g_buffers->positions[start + b]));
			float Lb = Length(Vec3(g_buffers->positions[start + b]) - Vec3(g_buffers->positions[start + c]));
			float Lc = Length(Vec3(g_buffers->positions[start + c]) - Vec3(g_buffers->positions[start + a]));
			
			if (L == 0.0f){
				L = La;
				cout << "length is " << L << endl;
			}
			CreateSpringInter(start + a, start + b, stiffness, 0.0f, (L == 0.0f) ? La : L);// main_radius*2.0f);
			CreateSpringInter(start + b, start + c, stiffness, 0.0f, (L == 0.0f) ? Lb : L);
			CreateSpringInter(start + c, start + a, stiffness, 0.0f, (L == 0.0f) ? Lc : L);
		}
	}

	virtual void Initialize()
	{
		g_rigidTranslations.resize(0);
		g_rigidRotations.resize(0);
		g_rigidCoefficients.resize(0);
		g_rigidIndices.resize(0);
		g_rigidLocalPositions.resize(0);
		g_rigidOffsets.resize(0);

		// start index
		g_rigidOffsets.push_back(0);
		main_scale = 1.0f / 100.0f;
		cp = new cellPACK();
		cp->use_rb = true;
		cp->main_radius = 1.0;//((10.0f/100.0f)*(10.0f/100.0f));
		cp->main_scale = main_scale;
		cp->maxParticles = 1024 * 1024;
		//actin filaments (13.5 μm) persistance length
		float atom_radius = 10.0f;
		float beads_radius = atom_radius*main_scale;//5.0f*main_scale;//0.08f;//5.0f*main_scale;//smallest sphere size
		g_params.mRadius = beads_radius * 2.0f;

		//read in the mesh - quad
		//build a network with vertices->beads
		//edge->spring with uniform length
		//relax

		string geompath = "C:\\Dev\\flexpack_dev_1.0\\data\\";
		string gname = "Mmycoides_PackingSurf_3_quad.obj";// "Mmycoides_PackingSurf_3.obj";
		cout << " read mesh 1 " << geompath + gname << endl;
		//cp->compartmentsSDF(gname);

		Mesh* mesh = GetMesh(GetFilePathByPlatform((geompath + gname).c_str()).c_str(), main_scale);
		//Quat q = QuatFromAxisAngle(Vec3(1, 0, 0), 3.14f / 2.0f);
		//mesh->Transform(RotationMatrix(q));
		//mesh->CalculateNormals();
		//cout << " ok mesh 1 " << geompath + gname << endl;
		//
		//FlexTriangleMesh* fmesh = CreateTriangleMesh(mesh);
		//AddTriangleMesh(fmesh, Vec3(), Quat(), 1.0f);
		Vec3 minExtents = Vec3(0, 0, 0);
		Vec3 maxExtents = Vec3(0, 0, 0);
		mesh->GetBounds(minExtents, maxExtents);
		cout << minExtents.x << " " << minExtents.y << " " << minExtents.z << endl;
		/*
		const int dim = 128;
		FlexSDF* sdf = CreateSDF(GetFilePathByPlatform("../../data/bunny.ply").c_str(), dim); 
		AddSDF(sdf, Vec3(-1.f, 0.0f, 0.0f), QuatFromAxisAngle(Vec3(0.0f, 1.0f, 0.0f), DegToRad(-45.0f)), 0.5f);
		//*/
		//SDF* sdf = CreateSDFfromMesh(GetFilePathByPlatform((geompath + gname).c_str()).c_str(), cp->comp_mesh[0],
		//	main_scale, Vec3(0.0f, 0.0f, 0.0f), minExtents, maxExtents, 0.0f);//main_radius/5.0f);
		float margin = 0.1f;
		FlexSDF* sdf = CreateSDF(GetFilePathByPlatform("../../data/Mmycoides_PackingSurf_3_quad.obj").c_str(), 256);//normalized s/maxedge
		AddSDF(sdf, minExtents - Vec3(margin, margin, margin)*0.5f* (maxExtents.x - minExtents.x) / 0.9f, QuatFromAxisAngle(Vec3(0.0f, 1.0f, 0.0f), DegToRad(0.0f)), (maxExtents.x - minExtents.x) / 0.9f);//maxEdge/(1.0-margin)

		//sdf process does scaling and trans
		//mesh->Normalize(1.0f - margin);
		//mesh->Transform(TranslationMatrix(Point3(margin, margin, margin)*0.5f));

		//g_fields[cp->comp_shape[0]->fsdf] = CreateGpuMesh(cp->comp_mesh[0]);
		//FlexSDF* sdf = CreateSDF(GetFilePathByPlatform("../../data/bunny.ply").c_str(), dim);
		//AddSDF(cp->comp_shape[0]->fsdf, Vec3(), Quat(), 1);
		//Vec3 lower, upper;
		//TransformBounds(Vec3(0.0f), Vec3(1.0f), Vec4(), Quat(), main_scale, lower, upper);
		/*
		FlexCollisionGeometry geo;
		geo.mSDF.mField = cp->comp_shape[0]->fsdf;
		geo.mSDF.mScale = 3;

		g_shapeStarts.push_back(g_shapeGeometry.size());
		g_shapeAabbMin.push_back(Vec4(minExtents, 0.0f));
		g_shapeAabbMax.push_back(Vec4(maxExtents, 0.0f));
		g_shapePositions.push_back(Vec4());
		g_shapeRotations.push_back(Quat());
		g_shapePrevPositions.push_back(Vec4());
		g_shapePrevRotations.push_back(Quat());
		g_shapeGeometry.push_back(geo);
		g_shapeFlags.push_back(flexMakeShapeFlags(eFlexShapeSDF, false));
		*/
		int phase = NvFlexMakePhase(99999, eFlexPhaseSelfCollide | eFlexPhaseSelfCollideFilter);//eFlexPhaseSelfCollide | eFlexPhaseSelfCollideFilter
		float L = 0.35f;// g_params.mRadius*2.0f;
		//buildMassSpringTri(cp->comp_mesh[0], 1.0f, 1.0f, phase, L);
		buildMassSpringTri(mesh, 1.0f, 1.0f, phase, L);

		//AddInflatable(mesh, 1, phase, L);
		//mSplitThreshold.resize(mCloths.size(), 45.0f);

		g_numExtraParticles = 1024 * 1024;
		g_params.mNumPlanes =0;
		//overwrite the plane
		(Vec4&)g_params.mPlanes[0] = Vec4(0.0f, 1.0f, 0.0f, -500.0);
		(Vec4&)g_params.mPlanes[1] = Vec4(0.0f, 0.0f, 1.0f, -cp->minExtents.z);// -cp->minExtents.z);
		(Vec4&)g_params.mPlanes[2] = Vec4(1.0f, 0.0f, 0.0f, -cp->minExtents.x);
		(Vec4&)g_params.mPlanes[3] = Vec4(-1.0f, 0.0f, 0.0f, cp->maxExtents.x);
		(Vec4&)g_params.mPlanes[4] = Vec4(0.0f, 0.0f, -1.0f, cp->maxExtents.z);// cp->maxExtents.z);
		(Vec4&)g_params.mPlanes[5] = Vec4(0.0f, -1.0f, 0.0f, cp->maxExtents.y);

		g_params.mGravity[1] = 0.0f;//-9.f;
		g_params.mDamping = 1.0f;// 3.0f;
		g_params.mRadius = beads_radius * 2;

		g_windStrength = 0.0f;
		g_windFrequency = 0.0f;

		g_params.mNumIterations = 5;

		g_params.mFluid = false;
		g_params.mVorticityConfinement = 0.0f;
		//g_params.mFluidRestDistance = g_params.mRadius*0.5;
		g_params.mAnisotropyScale = 2.5f / g_params.mRadius;
		g_params.mSmoothing = 0.5f;
		g_params.mRelaxationFactor = 1.f;
		g_params.mRestitution = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		g_params.mDynamicFriction = 0.25f;
		g_params.mViscosity = 1.5f;
		g_params.mCohesion = 1.1f;
		g_params.mAdhesion = 1.0f;
		g_params.mSurfaceTension = 1.0f;

		g_params.mGravity[1] = 0.0f;


		g_lightDistance *= 2.0f;

		// draw options		
		g_drawEllipsoids = true;
		g_params.mMaxSpeed = 100.0f;


		g_wireframe = true;
		g_pointScale = 1.0f;
		//        g_blur = 2.0f;
		g_pause = true;
		g_warmup = false;
		g_drawMesh = false;
		g_drawRopes = false;

		g_params.mDynamicFriction = 0.4f;
		g_params.mDissipation = 0.0f;

		g_params.mParticleCollisionMargin = g_params.mRadius*0.05f;
		g_params.mDrag = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		// better convergence with global relaxation factor
		g_params.mRelaxationMode = eFlexRelaxationGlobal;
		g_params.mRelaxationFactor = 0.25f;

		g_windStrength = 0.0f;

		g_numSubsteps = 5;
		g_params.mNumIterations = 5;

		// draw options		
		g_drawPoints = true;
		g_drawSprings = true;
		g_drawCloth = false;
	}

	void PostInitialize()
	{
		return;
		if (cp->comp_shape.size())
		{
			cp->fcontainer = flexExtCreateContainer(g_flex, cp->maxParticles);
			cp->sdfToForceField();
			cp->setForceField();
			//flexExtSetRigidTransformation(cp->fcontainer, cp->mInstances.size());
		}
	}

};




class MMTFTest : public Scene
{
public:

	MMTFTest(const char* name) : Scene(name) {}
	int iGroupCounter = 0;
	
	struct Branch
	{
		int npoints;
		std::vector<int> p_indices;
	};
	
	std::vector<Branch> all_branch;
	std::vector<Vec3> all_points;
	float main_scale = 1.0f / 100.0f;

	virtual void ParseMMTF(){
		//see the header and the api of mmtf in doc/index.html
		int phase = NvFlexMakePhase(iGroupCounter++, eFlexPhaseSelfCollide);
		string datapath = "..\\..\\data\\";
		string filename = datapath + "1HTQ.mmtf";
		//string filename = "C:\\Dev\\flexpack_dev_1.0\\data\\3j3q";
		MMTF_container* example = MMTF_container_new();
		MMTF_unpack_from_file(filename.c_str(), example);
		cout << example->mmtfVersion << endl;
		cout << "numAtromts" << example->numAtoms << endl;
		int N = example->numAtoms < 1024 * 1024 ? example->numAtoms : 1024 * 1024;
		for (int i = 0; i < N; i++){
			g_buffers->positions.push_back(Vec4(example->xCoordList[i] * main_scale, example->yCoordList[i] * main_scale, example->zCoordList[i] * main_scale, 1.0f));
			g_buffers->velocities.push_back(0.0f);
			g_buffers->phases.push_back(phase);//int(g_buffers->positions.size()));
		}
		MMTF_container_free(example);
	}

	virtual void Initialize()
	{
		g_rigidTranslations.resize(0);
		g_rigidRotations.resize(0);
		g_rigidCoefficients.resize(0);
		g_rigidIndices.resize(0);
		g_rigidLocalPositions.resize(0);
		g_rigidOffsets.resize(0);

		// start index
		g_rigidOffsets.push_back(0);

		main_scale = 1.0f / 100.0f;
		//actin filaments (13.5 μm) persistance length
		float atom_radius = 1.0f;
		float beads_radius = atom_radius*main_scale;//5.0f*main_scale;//0.08f;//5.0f*main_scale;//smallest sphere size
		g_params.mRadius = beads_radius * 2.0f;

		//read in the point and connection for all fiber

		ParseMMTF();
		
		g_numExtraParticles = 1024 * 1024;
		g_params.mNumPlanes = 0;

		g_params.mGravity[1] = 0.0f;//-9.f;
		g_params.mDamping = 1.0f;// 3.0f;
		g_params.mRadius = beads_radius * 2;

		g_windStrength = 0.0f;
		g_windFrequency = 0.0f;

		g_params.mNumIterations = 5;

		g_params.mFluid = false;
		g_params.mVorticityConfinement = 0.0f;
		//g_params.mFluidRestDistance = g_params.mRadius*0.5;
		g_params.mAnisotropyScale = 2.5f / g_params.mRadius;
		g_params.mSmoothing = 0.5f;
		g_params.mRelaxationFactor = 1.f;
		g_params.mRestitution = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		g_params.mDynamicFriction = 0.25f;
		g_params.mViscosity = 1.5f;
		g_params.mCohesion = 1.1f;
		g_params.mAdhesion = 1.0f;
		g_params.mSurfaceTension = 1.0f;

		g_params.mGravity[1] = 0.0f;


		g_lightDistance *= 2.0f;

		// draw options		
		g_drawEllipsoids = true;
		g_params.mMaxSpeed = 100.0f;


		g_wireframe = true;
		g_pointScale = 1.0f;
		//        g_blur = 2.0f;
		g_pause = true;
		g_warmup = false;
		g_drawMesh = false;
		g_drawRopes = false;

		g_params.mDynamicFriction = 0.4f;
		g_params.mDissipation = 0.0f;

		g_params.mParticleCollisionMargin = g_params.mRadius*0.05f;
		g_params.mDrag = 0.0f;
		g_params.mCollisionDistance = 0.01f;

		// better convergence with global relaxation factor
		g_params.mRelaxationMode = eFlexRelaxationGlobal;
		g_params.mRelaxationFactor = 0.25f;

		g_windStrength = 0.0f;

		g_numSubsteps = 5;
		g_params.mNumIterations = 5;

		// draw options		
		g_drawPoints = true;
		g_drawSprings = true;
		g_drawCloth = false;
	}




};
