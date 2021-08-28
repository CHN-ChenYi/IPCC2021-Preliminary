// SLIC.h: interface for the SLIC class.
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
//
//===========================================================================
// Copyright (c) 2013 Radhakrishna Achanta.
//
// For commercial use please contact the author:
//
// Email: firstname.lastname@epfl.ch
//===========================================================================

#if !defined(_SLIC_H_INCLUDED_)
#define _SLIC_H_INCLUDED_


#include <omp.h>

#include <vector>
#include <string>
#include <algorithm>
using namespace std;

#define REDUCE_NUM_THREADS 64
#define OMP_NUM_THREADS 64

struct lab_wrapper
{
	double l,a,b;
};

struct labxy_wrapper
{
	labxy_wrapper(double l=0, double a=0, double b=0, double x=0, double y=0):l(l),a(a),b(b),x(x),y(y) {}
	double l,a,b,x,y;
};

class SLIC  
{
public:
	SLIC();
	virtual ~SLIC();

	//============================================================================
	// Superpixel segmentation for a given number of superpixels
	//============================================================================
	void PerformSLICO_ForGivenK(
		const unsigned int*			ubuff,//Each 32 bit unsigned int contains ARGB pixel values.
		const int					width,
		const int					height,
		int*						klabels,
		int&						numlabels,
		const int&					K,
		const double&				m);

	//============================================================================
	// Save superpixel labels to pgm in raster scan order
	//============================================================================
	void SaveSuperpixelLabels2PPM(
		char*                       filename, 
		int *                       labels, 
		const int                   width, 
		const int                   height);

private:

	//============================================================================
	// Magic SLIC. No need to set M (compactness factor) and S (step size).
	// SLICO (SLIC Zero) varies only M dynamicaly, not S.
	//============================================================================
	void PerformSuperpixelSegmentation_VariableSandM(
		vector<labxy_wrapper>&		kseeds,
		int*						klabels,
		const int&					STEP,
		const int&					NUMITR);

	//============================================================================
	// Pick seeds for superpixels when number of superpixels is input.
	//============================================================================
	void GetLABXYSeeds_ForGivenK(
		vector<labxy_wrapper>&		kseeds,
		const int&					STEP,
		const bool&					perturbseeds,
		const vector<double>&		edges);

	//============================================================================
	// Move the seeds to low gradient positions to avoid putting seeds at region boundaries.
	//============================================================================
	void PerturbSeeds(
		vector<labxy_wrapper>&		kseeds,
		const vector<double>&		edges);
	
	//============================================================================
	// Detect color edges, to help PerturbSeeds()
	//============================================================================
	void DetectLabEdges(
		const lab_wrapper*			vec,
		const int&					width,
		const int&					height,
		vector<double>&				edges);
	
	//============================================================================
	// xRGB to XYZ conversion; helper for RGB2LAB()
	//============================================================================
	void RGB2XYZ(
		const int&					sR,
		const int&					sG,
		const int&					sB,
		double&						X,
		double&						Y,
		double&						Z);
	
	//============================================================================
	// sRGB to CIELAB conversion
	//============================================================================
	void RGB2LAB(
		const int&					sR,
		const int&					sG,
		const int&					sB,
		double&						lval,
		double&						aval,
		double&						bval);
	
	//============================================================================
	// sRGB to CIELAB conversion for 2-D images
	//============================================================================
	void DoRGBtoLABConversion(
		const unsigned int*&		ubuff,
		lab_wrapper*&				vec);

	//============================================================================
	// Post-processing of SLIC segmentation, to avoid stray labels.
	//============================================================================
	void EnforceLabelConnectivity(
		const int*					labels,
		const int&					width,
		const int&					height,
		int*						nlabels,//input labels that need to be corrected to remove stray labels
		int&						numlabels,//the number of labels changes in the end if segments are removed
		const int&					K); //the number of superpixels desired by the user

private:
	int										m_width;
	int										m_height;
	int										m_depth;

	// double*									m_lvec;
	// double*									m_avec;
	// double*									m_bvec;

	lab_wrapper * m_vec;

	double**								m_lvecvec;
	double**								m_avecvec;
	double**								m_bvecvec;

};

#endif // !defined(_SLIC_H_INCLUDED_)
