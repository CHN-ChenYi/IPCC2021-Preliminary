// SLIC.cpp: implementation of the SLIC class.
//===========================================================================
// This code implements the zero parameter superpixel segmentation technique
// described in:
//
//
//
// "SLIC Superpixels Compared to State-of-the-art Superpixel Methods"
//
// Radhakrishna Achanta, Appu Shaji, Kevin Smith, Aurelien Lucchi, Pascal Fua,
// and Sabine Susstrunk,
//
// IEEE TPAMI, Volume 34, Issue 11, Pages 2274-2282, November 2012.
//
// https://www.epfl.ch/labs/ivrl/research/slic-superpixels/
//===========================================================================
// Copyright (c) 2013 Radhakrishna Achanta.
//
// For commercial use please contact the author:
//
// Email: firstname.lastname@epfl.ch
//===========================================================================

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include "SLIC.h"
#include <chrono>

#include <omp.h>
#include <mpi.h>

typedef chrono::high_resolution_clock Clock;

// For superpixels
const int dx4[4] = {-1,  0,  1,  0};
const int dy4[4] = { 0, -1,  0,  1};
//const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
//const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

// For supervoxels
const int dx10[10] = {-1,  0,  1,  0, -1,  1,  1, -1,  0, 0};
const int dy10[10] = { 0, -1,  0,  1, -1, -1,  1,  1,  0, 0};
const int dz10[10] = { 0,  0,  0,  0,  0,  0,  0,  0, -1, 1};

#ifdef MYMPI
// For mpi
int world_rank;
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SLIC::SLIC()
{
	// m_lvec = NULL;
	// m_avec = NULL;
	// m_bvec = NULL;

	m_vec=NULL;

	m_lvecvec = NULL;
	m_avecvec = NULL;
	m_bvecvec = NULL;
}

SLIC::~SLIC()
{
	if(m_vec) delete [] m_vec;

	if(m_lvecvec)
	{
		for( int d = 0; d < m_depth; d++ ) delete [] m_lvecvec[d];
		delete [] m_lvecvec;
	}
	if(m_avecvec)
	{
		for( int d = 0; d < m_depth; d++ ) delete [] m_avecvec[d];
		delete [] m_avecvec;
	}
	if(m_bvecvec)
	{
		for( int d = 0; d < m_depth; d++ ) delete [] m_bvecvec[d];
		delete [] m_bvecvec;
	}
}

//==============================================================================
///	RGB2XYZ
///
/// sRGB (D65 illuninant assumption) to XYZ conversion
//==============================================================================
void SLIC::RGB2XYZ(
	const int&		sR,
	const int&		sG,
	const int&		sB,
	double&			X,
	double&			Y,
	double&			Z)
{
	double R = sR/255.0;
	double G = sG/255.0;
	double B = sB/255.0;

	double r, g, b;

	if(R <= 0.04045)	r = R/12.92;
	else				r = pow((R+0.055)/1.055,2.4);
	if(G <= 0.04045)	g = G/12.92;
	else				g = pow((G+0.055)/1.055,2.4);
	if(B <= 0.04045)	b = B/12.92;
	else				b = pow((B+0.055)/1.055,2.4);

	X = r*0.4124564 + g*0.3575761 + b*0.1804375;
	Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
	Z = r*0.0193339 + g*0.1191920 + b*0.9503041;
}

//===========================================================================
///	RGB2LAB
//===========================================================================
void SLIC::RGB2LAB(const int& sR, const int& sG, const int& sB, double& lval, double& aval, double& bval)
{
	//------------------------
	// sRGB to XYZ conversion
	//------------------------
	double X, Y, Z;
	RGB2XYZ(sR, sG, sB, X, Y, Z);

	//------------------------
	// XYZ to LAB conversion
	//------------------------
	double epsilon = 0.008856;	//actual CIE standard
	double kappa   = 903.3;		//actual CIE standard

	double Xr = 0.950456;	//reference white
	double Yr = 1.0;		//reference white
	double Zr = 1.088754;	//reference white

	double xr = X/Xr;
	double yr = Y/Yr;
	double zr = Z/Zr;

	double fx, fy, fz;
	if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
	else				fx = (kappa*xr + 16.0)/116.0;
	if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
	else				fy = (kappa*yr + 16.0)/116.0;
	if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
	else				fz = (kappa*zr + 16.0)/116.0;

	lval = 116.0*fy-16.0;
	aval = 500.0*(fx-fy);
	bval = 200.0*(fy-fz);
}

//===========================================================================
///	DoRGBtoLABConversion
///
///	For whole image: overlaoded floating point version
//===========================================================================
void SLIC::DoRGBtoLABConversion(
	const unsigned int*&		ubuff,
	lab_wrapper*&				vec)
{
	int sz = m_width*m_height;
	vec=new lab_wrapper[sz];

	#pragma omp parallel for
	for( int j = 0; j < sz; j++ )
	{
		int r = (ubuff[j] >> 16) & 0xFF;
		int g = (ubuff[j] >>  8) & 0xFF;
		int b = (ubuff[j]      ) & 0xFF;

		RGB2LAB( r, g, b, vec[j].l, vec[j].a, vec[j].b);
	}
	/*
	int mid = sz / 2;
	MPI_Status status;
	// MPI_Request request;
	if(world_rank==0)
	{
		#pragma omp parallel for
		for( int j = 0; j < mid; j++ )
		{
			int r = (ubuff[j] >> 16) & 0xFF;
			int g = (ubuff[j] >>  8) & 0xFF;
			int b = (ubuff[j]      ) & 0xFF;

			RGB2LAB( r, g, b, vec[j].l, vec[j].a, vec[j].b );
		}
		// MPI_Isend(lvec,mid,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&request[0]);
		// MPI_Isend(avec,mid,MPI_DOUBLE,1,2,MPI_COMM_WORLD,&request[1]);
		// MPI_Isend(bvec,mid,MPI_DOUBLE,1,4,MPI_COMM_WORLD,&request[2]);
		// MPI_Irecv(lvec+mid,sz-mid,MPI_DOUBLE,1,1,MPI_COMM_WORLD,&request[3]);
		// MPI_Irecv(avec+mid,sz-mid,MPI_DOUBLE,1,3,MPI_COMM_WORLD,&request[4]);
		// MPI_Irecv(bvec+mid,sz-mid,MPI_DOUBLE,1,5,MPI_COMM_WORLD,&request[5]);
		// MPI_Sendrecv(lvec,mid,MPI_DOUBLE,1,0,lvec+mid,sz-mid,MPI_DOUBLE,1,1,MPI_COMM_WORLD,&status[0]);
		// MPI_Sendrecv(avec,mid,MPI_DOUBLE,1,2,avec+mid,sz-mid,MPI_DOUBLE,1,3,MPI_COMM_WORLD,&status[1]);
		// MPI_Sendrecv(bvec,mid,MPI_DOUBLE,1,4,bvec+mid,sz-mid,MPI_DOUBLE,1,5,MPI_COMM_WORLD,&status[2]);
		MPI_Sendrecv(vec,mid*3,MPI_DOUBLE,1,0,vec+mid,(sz-mid)*3,MPI_DOUBLE,1,1,MPI_COMM_WORLD,&status);
	}
	else
	{
		#pragma omp parallel for
		for( int j = mid; j < sz; j++ )
		{
			int r = (ubuff[j] >> 16) & 0xFF;
			int g = (ubuff[j] >>  8) & 0xFF;
			int b = (ubuff[j]      ) & 0xFF;

			RGB2LAB( r, g, b, vec[j].l, vec[j].a, vec[j].b );
		}
		// MPI_Isend(lvec+mid,sz-mid,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&request[0]);
		// MPI_Isend(avec+mid,sz-mid,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&request[1]);
		// MPI_Isend(bvec+mid,sz-mid,MPI_DOUBLE,0,5,MPI_COMM_WORLD,&request[2]);
		// MPI_Irecv(lvec,mid,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&request[3]);
		// MPI_Irecv(avec,mid,MPI_DOUBLE,0,2,MPI_COMM_WORLD,&request[4]);
		// MPI_Irecv(bvec,mid,MPI_DOUBLE,0,4,MPI_COMM_WORLD,&request[5]);
		// MPI_Sendrecv(lvec+mid,sz-mid,MPI_DOUBLE,0,1,lvec,mid,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status[0]);
		// MPI_Sendrecv(avec+mid,sz-mid,MPI_DOUBLE,0,3,avec,mid,MPI_DOUBLE,0,2,MPI_COMM_WORLD,&status[1]);
		// MPI_Sendrecv(bvec+mid,sz-mid,MPI_DOUBLE,0,5,bvec,mid,MPI_DOUBLE,0,4,MPI_COMM_WORLD,&status[2]);
		MPI_Sendrecv(vec+mid,(sz-mid)*3,MPI_DOUBLE,0,1,vec,mid*3,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
	}
	// MPI_Waitall(6,request,status);
	*/
}

//==============================================================================
///	DetectLabEdges
//==============================================================================
void SLIC::DetectLabEdges(
	const lab_wrapper*			vec,
	const int&					width,
	const int&					height,
	vector<double>&				edges)
{
	int sz = width*height;

	edges.resize(sz,0);
	#pragma omp parallel for
	for( int j = 1; j < height-1; j++ )
	{
		for( int k = 1; k < width-1; k++ )
		{
			int i = j*width+k;

			double dx = (vec[i-1].l-vec[i+1].l)*(vec[i-1].l-vec[i+1].l) +
						(vec[i-1].a-vec[i+1].a)*(vec[i-1].a-vec[i+1].a) +
						(vec[i-1].b-vec[i+1].b)*(vec[i-1].b-vec[i+1].b);

			double dy = (vec[i-width].l-vec[i+width].l)*(vec[i-width].l-vec[i+width].l) +
						(vec[i-width].a-vec[i+width].a)*(vec[i-width].a-vec[i+width].a) +
						(vec[i-width].b-vec[i+width].b)*(vec[i-width].b-vec[i+width].b);

			//edges[i] = (sqrt(dx) + sqrt(dy));
			edges[i] = (dx + dy);
		}
	}
}

//===========================================================================
///	PerturbSeeds
//===========================================================================
void SLIC::PerturbSeeds(
	vector<labxy_wrapper>&		kseeds,
	const vector<double>&		edges)
{
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	
	int numseeds = kseeds.size();

	for( int n = 0; n < numseeds; n++ )
	{
		int ox = kseeds[n].x;//original x
		int oy = kseeds[n].y;//original y
		int oind = oy*m_width + ox;

		int storeind = oind;
		for( int i = 0; i < 8; i++ )
		{
			int nx = ox+dx8[i];//new x
			int ny = oy+dy8[i];//new y

			if( nx >= 0 && nx < m_width && ny >= 0 && ny < m_height)
			{
				int nind = ny*m_width + nx;
				if( edges[nind] < edges[storeind])
				{
					storeind = nind;
				}
			}
		}
		if(storeind != oind)
		{
			kseeds[n].x = storeind%m_width;
			kseeds[n].y = storeind/m_width;
			kseeds[n].l = m_vec[storeind].l;
			kseeds[n].a = m_vec[storeind].a;
			kseeds[n].b = m_vec[storeind].b;
		}
	}
}

//===========================================================================
///	GetLABXYSeeds_ForGivenK
///
/// The k seed values are taken as uniform spatial pixel samples.
//===========================================================================
void SLIC::GetLABXYSeeds_ForGivenK(
	vector<labxy_wrapper>&		kseeds,
	const int&					K,
	const bool&					perturbseeds,
	const vector<double>&		edgemag)
{
	int sz = m_width*m_height;
	double step = sqrt(double(sz)/double(K));
	int T = step;
	int xoff = step/2;
	int yoff = step/2;
	
	int n(0);int r(0);
	for( int y = 0; y < m_height; y++ )
	{
		int Y = y*step + yoff;
		if( Y > m_height-1 ) break;

		for( int x = 0; x < m_width; x++ )
		{
			//int X = x*step + xoff;//square grid
			int X = x*step + (xoff<<(r&0x1));//hex grid
			if(X > m_width-1) break;

			int i = Y*m_width + X;

			//_ASSERT(n < K);
			
			//kseeds[n] = m_vec[i].l;
			//kseeds[n] = m_vec[i].a;
			//kseeds[n] = m_vec[i].b;
			//kseeds[n] = X;
			//kseeds[n] = Y;
			kseeds.push_back(labxy_wrapper(m_vec[i].l,m_vec[i].a,m_vec[i].b,X,Y));
			n++;
		}
		r++;
	}

	if(perturbseeds)
	{
		PerturbSeeds(kseeds, edgemag);
	}
}

//===========================================================================
///	PerformSuperpixelSegmentation_VariableSandM
///
///	Magic SLIC - no parameters
///
///	Performs k mean segmentation. It is fast because it looks locally, not
/// over the entire image.
/// This function picks the maximum value of color distance as compact factor
/// M and maximum pixel distance as grid step size S from each cluster (13 April 2011).
/// So no need to input a constant value of M and S. There are two clear
/// advantages:
///
/// [1] The algorithm now better handles both textured and non-textured regions
/// [2] There is not need to set any parameters!!!
///
/// SLICO (or SLIC Zero) dynamically varies only the compactness factor S,
/// not the step size S.
//===========================================================================
void SLIC::PerformSuperpixelSegmentation_VariableSandM(
	vector<labxy_wrapper>&		kseeds,
	int*						klabels,
	const int&					STEP,
	const int&					NUMITR)
{
	int sz = m_width*m_height;
	const int numk = kseeds.size();
	//double cumerr(99999.9);
	int numitr(0);

	//----------------
	int offset = STEP;
	if(STEP < 10) offset = STEP*1.5;
	//----------------


	vector<labxy_wrapper> sigma(numk);
	vector<int> clustersize(numk, 0);
	vector<double> inv(numk, 0);//to store 1/clustersize[k] values
	vector<double> distxy(sz, DBL_MAX);
	vector<double> distlab(sz, DBL_MAX);
	vector<double> distvec(sz, DBL_MAX);
	vector<double> maxlab(numk, 10*10);//THIS IS THE VARIABLE VALUE OF M, just start with 10
	vector<double> maxxy(numk, STEP*STEP);//THIS IS THE VARIABLE VALUE OF M, just start with 10

	double invxywt = 1.0/(STEP*STEP);//NOTE: this is different from how usual SLIC/LKM works

	// vector<int> distidx(sz, -1);
	int *distidx = new int[sz];
	double *_maxlab[OMP_NUM_THREADS];
	double *_maxxy[OMP_NUM_THREADS];
	// double *_sigmal[OMP_NUM_THREADS];
	// double *_sigmaa[OMP_NUM_THREADS];
	// double *_sigmab[OMP_NUM_THREADS];
	// double *_sigmax[OMP_NUM_THREADS];
	// double *_sigmay[OMP_NUM_THREADS];
	labxy_wrapper *_sigma[OMP_NUM_THREADS];
	int *_clustersize[OMP_NUM_THREADS];

	// memset(distidx, 0, sizeof(int) * sz);
	#pragma omp parallel for
	for (int i = 0; i < OMP_NUM_THREADS; ++i) {
		_maxlab[i] = new double[numk];
		_maxxy[i] = new double[numk];
		_sigma[i] = new labxy_wrapper[numk];
		_clustersize[i] = new int[numk];
	}
#ifdef PROF
	auto cost1 = 0LL;
	auto cost2 = 0LL;
	auto cost3 = 0LL;
	long long ccnt = 0;
	chrono::time_point<std::chrono::high_resolution_clock> startTime, endTime;
	chrono::microseconds compTime;
#endif
	while( numitr < NUMITR )
	{
		//------
		//cumerr = 0;
		numitr++;
		//------

#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
		startTime = Clock::now();
	}
#endif
		#pragma omp parallel for
		for (int i = 0; i < sz; ++i) {
			distvec[i] = DBL_MAX;
		}

		for( int n = 0; n < numk; n++ )
		{
			const double _kseedsl = kseeds[n].l;
			const double _kseedsa = kseeds[n].a;
			const double _kseedsb = kseeds[n].b;
			const double _kseedsx = kseeds[n].x;
			const double _kseedsy = kseeds[n].y;
			const int y1 = max(0,			(int)(_kseedsy-offset));
			const int y2 = min(m_height,	(int)(_kseedsy+offset));
			const int x1 = max(0,			(int)(_kseedsx-offset));
			const int x2 = min(m_width,	(int)(_kseedsx+offset));

			#pragma omp parallel for
			for( int y = y1; y < y2; y++ )
			{
				for( int x = x1; x < x2; x++ )
				{
					const int i = y*m_width + x;
					//_ASSERT( y < m_height && x < m_width && y >= 0 && x >= 0 );

					const double l = m_vec[i].l;
					const double a = m_vec[i].a;
					const double b = m_vec[i].b;

					const double _distlab =	(l - _kseedsl)*(l - _kseedsl) +
											(a - _kseedsa)*(a - _kseedsa) +
											(b - _kseedsb)*(b - _kseedsb);

					const double _distxy =	(x - _kseedsx)*(x - _kseedsx) +
											(y - _kseedsy)*(y - _kseedsy);

					//------------------------------------------------------------------------
					const double dist = _distlab/maxlab[n] + _distxy*invxywt;//only varying m, prettier superpixels
					//double dist = distlab[i]/maxlab[n] + distxy[i]/maxxy[n];//varying both m and S
					//------------------------------------------------------------------------

					distlab[i] = _distlab;
					distxy[i] = _distxy;
					
					if( dist < distvec[i] )
					{
						distvec[i] = dist;
						klabels[i]  = n;
					}
				}
			}
		}

#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
    	endTime = Clock::now();
    	compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
		
   		cost1 += compTime.count();
	}
#endif
		//-----------------------------------------------------------------
		// Assign the max color distance for a cluster
		//-----------------------------------------------------------------
		// if(0 == numitr)	// I think it won't be executed
		// {
		// 	maxlab.assign(numk,1);
		// 	maxxy.assign(numk,1);
		// }
#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
		startTime = Clock::now();
	}
#endif
		#pragma omp parallel for
		for (int i = 0; i < OMP_NUM_THREADS; ++i) {
			memset(_maxlab[i], 0, sizeof(double) * numk);
			memset(_maxxy[i], 0, sizeof(double) * numk);
		}
		#pragma omp parallel for
		for( int i = 0; i < sz; i++ ) {
			int idx = omp_get_thread_num();
			if(_maxlab[idx][klabels[i]] < distlab[i]) _maxlab[idx][klabels[i]] = distlab[i];
			if(_maxxy[idx][klabels[i]] < distxy[i]) _maxxy[idx][klabels[i]] = distxy[i];
		}
		#pragma omp parallel for
		for (int i = 0; i < numk; ++i)
			for (int j = 0; j < OMP_NUM_THREADS; ++j) {
				if(maxlab[i] < _maxlab[j][i]) maxlab[i] = _maxlab[j][i];
				if(maxxy[i] < _maxxy[j][i]) maxxy[i] = _maxxy[j][i];
			}
#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
    	endTime = Clock::now();
    	compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
		
   		cost2 += compTime.count();
	}
#endif
		//-----------------------------------------------------------------
		// Recalculate the centroid and store in the seed values
		//-----------------------------------------------------------------
#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
		startTime = Clock::now();
	}
#endif
		sigma.assign(numk,labxy_wrapper());
		clustersize.assign(numk, 0);

		
		#pragma omp parallel for
		for (int i = 0; i < OMP_NUM_THREADS; ++i) {
			memset(_sigma[i], 0, 5 * sizeof(double) * numk);
			memset(_clustersize[i], 0, sizeof(int) * numk);
		}
		#pragma omp parallel
		{
		int idx = omp_get_thread_num();
		#pragma omp for
		for (int j = 0; j < sz; ++j) {
			_sigma[idx][klabels[j]].l += m_vec[j].l;
			_sigma[idx][klabels[j]].a += m_vec[j].a;
			_sigma[idx][klabels[j]].b += m_vec[j].b;
			_sigma[idx][klabels[j]].x += (j%m_width);
			_sigma[idx][klabels[j]].y += (j/m_width);

			_clustersize[idx][klabels[j]]++;
		}
		}
		#pragma omp parallel for
		for (int i = 0; i < numk; ++i)
			for (int j = 0; j < OMP_NUM_THREADS; ++j) {
				sigma[i].l += _sigma[j][i].l;
				sigma[i].a += _sigma[j][i].a;
				sigma[i].b += _sigma[j][i].b;
				sigma[i].x += _sigma[j][i].x;
				sigma[i].y += _sigma[j][i].y;
				clustersize[i] += _clustersize[j][i];
			}

#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
    	endTime = Clock::now();
    	compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
		
   		cost3 += compTime.count();
	}
#endif

		#pragma omp parallel
		{
		#pragma omp for
		for( int k = 0; k < numk; k++ )
		{
			//_ASSERT(clustersize[k] > 0);
			if( clustersize[k] <= 0 ) clustersize[k] = 1;
			inv[k] = 1.0/double(clustersize[k]);//computing inverse now to multiply, than divide later
		}
		}
		
		#pragma omp parallel
		{
		#pragma omp for
		for( int k = 0; k < numk; k++ )
		{
			kseeds[k].l = sigma[k].l*inv[k];
			kseeds[k].a = sigma[k].a*inv[k];
			kseeds[k].b = sigma[k].b*inv[k];
			kseeds[k].x = sigma[k].x*inv[k];
			kseeds[k].y = sigma[k].y*inv[k];
		}
		}
	}

#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
    cout << "Part1 time=" << cost1/1000 << " ms" << endl;
    cout << "Part2 time=" << cost2/1000 << " ms" << endl;
    cout << "Part3 time=" << cost3/1000 << " ms" << endl;
	printf("sz = %d\n", sz);
	printf("numk = %d\n", numk);
	printf("ccnt = %lld\n", ccnt);
	}
#endif

	delete [] distidx;
	#pragma omp parallel for
	for (int i = 0; i < OMP_NUM_THREADS; ++i) {
		delete [] _maxlab[i];
		delete [] _maxxy[i];
		delete [] _sigma[i];
		delete [] _clustersize[i];
	}
}

//===========================================================================
///	SaveSuperpixelLabels2PGM
///
///	Save labels to PGM in raster scan order.
//===========================================================================
void SLIC::SaveSuperpixelLabels2PPM(
	char*                           filename, 
	int *                           labels, 
	const int                       width, 
	const int                       height)
{
    FILE* fp;
    char header[20];
 
    fp = fopen(filename, "wb");
 
    // write the PPM header info, such as type, width, height and maximum
    fprintf(fp,"P6\n%d %d\n255\n", width, height);
 
    // write the RGB data
    unsigned char *rgb = new unsigned char [ (width)*(height)*3 ];
    int k = 0;
	unsigned char c = 0;
    for ( int i = 0; i < (height); i++ ) {
        for ( int j = 0; j < (width); j++ ) {
			c = (unsigned char)(labels[k]);
            rgb[i*(width)*3 + j*3 + 2] = labels[k] >> 16 & 0xff;  // r
            rgb[i*(width)*3 + j*3 + 1] = labels[k] >> 8  & 0xff;  // g
            rgb[i*(width)*3 + j*3 + 0] = labels[k]       & 0xff;  // b

			// rgb[i*(width) + j + 0] = c;
            k++;
        }
    }
    fwrite(rgb, width*height*3, 1, fp);

    delete [] rgb;
 
    fclose(fp);

}

//===========================================================================
///	EnforceLabelConnectivity
///
///		1. finding an adjacent label for each new component at the start
///		2. if a certain component is too small, assigning the previously found
///		    adjacent label to this component, and not incrementing the label.
//===========================================================================
void SLIC::EnforceLabelConnectivity(
	const int*					labels,//input labels that need to be corrected to remove stray labels
	const int&					width,
	const int&					height,
	int*						nlabels,//new labels
	int&						numlabels,//the number of labels changes in the end if segments are removed
	const int&					K) //the number of superpixels desired by the user
{
//	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
//	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

	const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};

	const int sz = width*height;
	const int SUPSZ = sz/K;
	//nlabels.resize(sz, -1);
	memset(nlabels, -1, sizeof(int) * sz);
	// for( int i = 0; i < sz; i++ ) nlabels[i] = -1;
	int label(0);
	int* xvec = new int[sz];
	int* yvec = new int[sz];
	int oindex(0);
	int adjlabel(0);//adjacent label
	for( int j = 0; j < height; j++ )
	{
		for( int k = 0; k < width; k++ )
		{
			if( 0 > nlabels[oindex] )
			{
				nlabels[oindex] = label;
				//--------------------
				// Start a new segment
				//--------------------
				xvec[0] = k;
				yvec[0] = j;
				//-------------------------------------------------------
				// Quickly find an adjacent label for use later if needed
				//-------------------------------------------------------
				{for( int n = 0; n < 4; n++ )
				{
					int x = xvec[0] + dx4[n];
					int y = yvec[0] + dy4[n];
					if( (x >= 0 && x < width) && (y >= 0 && y < height) )
					{
						int nindex = y*width + x;
						if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
					}
				}}

				int count(1);
				for( int c = 0; c < count; c++ )
				{
					for( int n = 0; n < 4; n++ )
					{
						int x = xvec[c] + dx4[n];
						int y = yvec[c] + dy4[n];

						if( (x >= 0 && x < width) && (y >= 0 && y < height) )
						{
							int nindex = y*width + x;

							if( 0 > nlabels[nindex] && labels[oindex] == labels[nindex] )
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;
							}
						}

					}
				}
				//-------------------------------------------------------
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				//-------------------------------------------------------
				if(count <= SUPSZ >> 2)
				{
					for( int c = 0; c < count; c++ )
					{
						int ind = yvec[c]*width+xvec[c];
						nlabels[ind] = adjlabel;
					}
					label--;
				}
				label++;
			}
			oindex++;
		}
	}
	numlabels = label;

	if(xvec) delete [] xvec;
	if(yvec) delete [] yvec;
}

//===========================================================================
///	PerformSLICO_ForGivenK
///
/// Zero parameter SLIC algorithm for a given number K of superpixels.
//===========================================================================
void SLIC::PerformSLICO_ForGivenK(
	const unsigned int*			ubuff,
	const int					width,
	const int					height,
	int*						klabels,
	int&						numlabels,
	const int&					K,//required number of superpixels
	const double&				m)//weight given to spatial distance
{
	// vector<double> kseedsl(0);
	// vector<double> kseedsa(0);
	// vector<double> kseedsb(0);
	// vector<double> kseedsx(0);
	// vector<double> kseedsy(0);

	vector<labxy_wrapper> kseeds(0);

	//--------------------------------------------------
	m_width  = width;
	m_height = height;
	int sz = m_width*m_height;
	//--------------------------------------------------
	//if(0 == klabels) klabels = new int[sz];
	memset(klabels, -1, sizeof(int) * sz);
	// for( int s = 0; s < sz; s++ ) klabels[s] = -1;
	//--------------------------------------------------

#ifdef PROF
	chrono::time_point<std::chrono::high_resolution_clock> startTime, endTime;
	chrono::microseconds compTime;
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
    startTime = Clock::now();
	}
#endif
	if(1)//LAB
	{
		DoRGBtoLABConversion(ubuff, m_vec);
	}
	else//RGB
	{
		// m_lvec = new double[sz]; m_avec = new double[sz]; m_bvec = new double[sz];
		m_vec = new lab_wrapper[sz];
		for( int i = 0; i < sz; i++ )
		{
			m_vec[i].l = ubuff[i] >> 16 & 0xff;
			m_vec[i].a = ubuff[i] >>  8 & 0xff;
			m_vec[i].b = ubuff[i]       & 0xff;
		}
	}
#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
    endTime = Clock::now();
    compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
    cout <<  "DoRGBtoLABConversion time=" << compTime.count()/1000 << " ms" << endl;
	}
#endif
	//--------------------------------------------------

	bool perturbseeds(true);
	vector<double> edgemag(0);
#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
    startTime = Clock::now();
	}
#endif
	if(perturbseeds) DetectLabEdges(m_vec, m_width, m_height, edgemag);
	GetLABXYSeeds_ForGivenK(kseeds, K, perturbseeds, edgemag);
#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
    endTime = Clock::now();
    compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
    cout <<  "GetLABXYSeeds_ForGivenK time=" << compTime.count()/1000 << " ms" << endl;
	}
#endif

#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
    startTime = Clock::now();
	}
#endif
	int STEP = sqrt(double(sz)/double(K)) + 2.0;//adding a small value in the even the STEP size is too small.
	PerformSuperpixelSegmentation_VariableSandM(kseeds,klabels,STEP,10);
	numlabels = kseeds.size();
#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
    endTime = Clock::now();
    compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
    cout <<  "PerformSuperpixelSegmentation_VariableSandM time=" << compTime.count()/1000 << " ms" << endl;
	}
#endif

	int* nlabels = new int[sz];
#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
	startTime = Clock::now();
	}
#endif
	EnforceLabelConnectivity(klabels, m_width, m_height, nlabels, numlabels, K);
	{
		memcpy(klabels, nlabels, sizeof(int) * sz);
		// for(int i = 0; i < sz; i++) klabels[i] = nlabels[i];
	}
#ifdef PROF
#ifdef MYMPI
	if (world_rank == 0)
#endif
	{
	endTime = Clock::now();
	compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
	cout <<  "EnforceLabelConnectivity time=" << compTime.count()/1000 << " ms" << endl;
	}
#endif
	if(nlabels) delete [] nlabels;
}

//===========================================================================
/// Load PPM file
///
/// 
//===========================================================================
void LoadPPM(char* filename, unsigned int** data, int* width, int* height)
{
    char header[1024];
    FILE* fp = NULL;
    int line = 0;
 
    fp = fopen(filename, "rb");
 
    // read the image type, such as: P6
    // skip the comment lines
    while (line < 2) {    
        fgets(header, 1024, fp);
        if (header[0] != '#') {
            ++line;
        }
    }
    // read width and height
    sscanf(header,"%d %d\n", width, height);
 
    // read the maximum of pixels
    fgets(header, 20, fp);
 
    // get rgb data
    unsigned char *rgb = new unsigned char [ (*width)*(*height)*3 ];
    fread(rgb, (*width)*(*height)*3, 1, fp);

    *data = new unsigned int [ (*width)*(*height)*4 ];
    int k = 0;
    for ( int i = 0; i < (*height); i++ ) {
        for ( int j = 0; j < (*width); j++ ) {
            unsigned char *p = rgb + i*(*width)*3 + j*3;
                                      // a ( skipped )
            (*data)[k]  = p[2] << 16; // r
            (*data)[k] |= p[1] << 8;  // g
            (*data)[k] |= p[0];       // b
            k++;
        }
    }

    // ofc, later, you'll have to cleanup
    delete [] rgb;
 
    fclose(fp);
}

//===========================================================================
/// Load PPM file
///
/// 
//===========================================================================
int CheckLabelswithPPM(char* filename, int* labels, int width, int height)
{
    char header[1024];
    FILE* fp = NULL;
    int line = 0, ground = 0;
 
    fp = fopen(filename, "rb");
 
    // read the image type, such as: P6
    // skip the comment lines
    while (line < 2) {    
        fgets(header, 1024, fp);
        if (header[0] != '#') {
            ++line;
        }
    }
    // read width and height
	int w(0);
	int h(0);
    sscanf(header,"%d %d\n", &w, &h);
	if (w != width || h != height) return -1;
 
    // read the maximum of pixels
    fgets(header, 20, fp);
 
    // get rgb data
    unsigned char *rgb = new unsigned char [ (w)*(h)*3 ];
    fread(rgb, (w)*(h)*3, 1, fp);

    int num = 0, k = 0;
    for ( int i = 0; i < (h); i++ ) {
        for ( int j = 0; j < (w); j++ ) {
            unsigned char *p = rgb + i*(w)*3 + j*3;
                                  // a ( skipped )
            ground  = p[2] << 16; // r
            ground |= p[1] << 8;  // g
            ground |= p[0];       // b
            
			if (ground != labels[k])
				num++;

			k++;
        }
    }

    // ofc, later, you'll have to cleanup
    delete [] rgb;
 
    fclose(fp);

	return num;
}

//===========================================================================
///	The main function
//===========================================================================
int main(int argc, char **argv)
{
#ifdef MYMPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif
	unsigned int *img = NULL;
	int width(0);
	int height(0);
#ifdef SLCT
	//===========================case_switch=================================
	const int name_len = 100;
	char input_flag = '1';
	if (argc == 2)
		input_flag = argv[1][0];
	int m_spcount;
	switch (input_flag)
	{
	case '1':
		m_spcount = 200;
		break;
	case '2':
		m_spcount = 400;
		break;
	case '3':
		m_spcount = 150;
		break;
	default:
		fprintf(stderr, "argument error\n");
		exit(0);
		break;
	}
	char input_name[name_len], check_name[name_len], output_name[name_len];
	sprintf(input_name, "case%c/input_image.ppm", input_flag);
	sprintf(check_name, "case%c/check.ppm", input_flag);
	sprintf(output_name, "case%c/output_labels.ppm", input_flag);
	//==========================case_switch_end==============================
	LoadPPM(input_name, &img, &width, &height);
#else
	int m_spcount = 200;
	LoadPPM((char *)"input_image.ppm", &img, &width, &height);
#endif

	if (width == 0 || height == 0)
		return -1;

	int sz = width * height;
	int *labels = new int[sz];
	int numlabels(0);
	SLIC slic;
	// int m_spcount;
	double m_compactness;
	// m_spcount = 200;
	m_compactness = 10.0;
	auto startTime = Clock::now();

	slic.PerformSLICO_ForGivenK(img, width, height, labels, numlabels, m_spcount, m_compactness); //for a given number K of superpixels

#ifdef MYMPI
	if(world_rank==0)
	{
#endif
    auto endTime = Clock::now();
    auto compTime = chrono::duration_cast<chrono::microseconds>(endTime - startTime);
    cout <<  "Computing time=" << compTime.count()/1000 << " ms" << endl;

#ifdef SLCT
	int num = CheckLabelswithPPM(check_name, labels, width, height);
#else
	int num = CheckLabelswithPPM((char *)"check.ppm", labels, width, height);
#endif

	if (num < 0)
	{
		cout << "The result for labels is different from output_labels.ppm." << endl;
	}
	else
	{
		cout << "There are " << num << " points' labels are different from original file." << endl;
	}

#ifdef SLCT
	slic.SaveSuperpixelLabels2PPM(output_name, labels, width, height);
#else
	slic.SaveSuperpixelLabels2PPM((char *)"output_labels.ppm", labels, width, height);
#endif

#ifdef MYMPI
	}
#endif

	if (labels)
		delete[] labels;
	if (img)
		delete[] img;

#ifdef MYMPI
	MPI_Finalize();
#endif
	return 0;
}