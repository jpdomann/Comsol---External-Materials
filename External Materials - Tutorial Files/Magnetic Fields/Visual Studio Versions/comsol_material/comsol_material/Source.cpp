/******************************************************************************************/
/******************************************************************************************/
/*					Interface to a magnetic saturation model.							  */
/******************************************************************************************/
/******************************************************************************************/


/******************************************************************************************/
/******************************************************************************************/
/*								Header Section											  */
/******************************************************************************************/
/******************************************************************************************/
/* Load libraries */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <limits>

//THIS IS DIFFERENT THAN WHAT IS IN THE ORIGINAL COMSOL FILES!!!
//THIS NEEDS TO BE UPDATED TO AVOID NAME MANGLING
#ifdef _MSC_VER									//Check if code is being compiled with visual studio 
#define EXPORT extern "C" __declspec(dllexport)	//Define command to export symbols to a .dll
#else
#define EXPORT
#endif

/* Function Declarations */
void interp1d(const double* Hg, const double* Bg, const double*  H, double* B, double* Jac, int const Ng);


/******************************************************************************************/
/******************************************************************************************/
/*								Main "eval" function									  */
/******************************************************************************************/
/******************************************************************************************/

/** Interface to a general B(H) relation with a variable number of material
*   model parameters and a variable number of states. The example code implements
*   an interpolation function representative of soft iron. */  
EXPORT int eval(double *oldH,      // Magnetic field, previous converged step, 3-vector, input
	double *H,         // Magnetic field, 3-vector, input
	double *B,         // Magnetic flux density, 3-vector, state
	double *Jac,       // Jacobian of B with respect to H, 3-by-3 matrix in row-major order, output
	int *nPar,         // Number of material model parameters, scalar, input
	double *par,       // Material model parameters, nPar-vector, input
	int *nStates,      // Number of states, scalar, input
	double *states)	   // States, nStates-vector
{  
	
	// Interpolation data
	int const nInt = 16;
	double const intH[16] = { 0.0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3, 20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310 };
	double const intB[16] = { 0.0, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4 };

	//pointers to start of interpolation data
	const double* pH = &intH[0];	
	const double* pB = &intB[0];

	// Check inputs (these are the only error codes Comsol natively understands)
	if (nPar[0] != 0)
		return 1;	              // error code 1 = "Wrong number of parameters"
	if (nStates[0] != 0)
		return 2;                 // error code 2 = "Wrong number of states"

	// Interpolate
	interp1d( pH, pB, H, B, Jac, nInt);
	// note that B is passed as a pointer



	return 0;  // Return value 0 is success, any other value triggers an exception
}//eval (main function called by Comsol)


/******************************************************************************************/
/******************************************************************************************/
/*									Function Definitions								  */
/******************************************************************************************/
/******************************************************************************************/

void interp1d(const double* Hg,			//pointer to known H points - input
	const double* Bg,					//pointer to known B points - input
	const double* H,					//pointer to H points where you want to know B - input
	double *B,							//pointer to B - output
	double *Jac,						//pointer to 3x3 jacobian matrix for B Jac[i][j] = dB[i]/dH[j]  - output
	int const nInt){					//number of points in Hg and Bg - input
	
	/* 1 Dimensional linear interoplation function*/

	int i, j, pos;
	const double eps = 3e-16;	
	double normH2, normH, normEpsH2, normEpsH, normB, dNormB, s;

	// Compute input norm
	normH2 = H[0] * H[0] + H[1] * H[1] + H[2] * H[2];

	if (normH2 > 1000){
		int debug_time = 1;
	}
	normEpsH2 = normH2 + eps;
	normH = sqrt(normH2);
	normEpsH = sqrt(normEpsH2);
	// Find the interpolation interval 
	pos = 0;
	while (pos < nInt && Hg[pos] < normH) {
		++pos;
	}
	// Handle extrapolation
	if (pos == 0) {
		pos = 1;
	}
	else if (pos == nInt) {
		pos = nInt - 1;
	}
	// Interpolate
	s = (normH - Hg[pos - 1]) / (Hg[pos] - Hg[pos - 1]);		//relative dist. between H points 
	normB = s*Bg[pos] + (1.0 - s)*Bg[pos - 1];
	dNormB = (Bg[pos] - Bg[pos - 1]) / (Hg[pos] - Hg[pos - 1]);
	
	// Write B vector
	if (normH > 0) {
		for (i = 0; i < 3; ++i) {
			B[i] = normB*H[i] / normEpsH;
		}
	}
	else {
		for (i = 0; i<3; ++i) {
			B[i] = dNormB*H[i];
		}
	}
	// Write Jacobian of B with respect to H
	if (normH>0) {
		for (i = 0; i < 3; ++i) {
			for (j = 0; j < 3; ++j) {
				Jac[i * 3 + j] = dNormB*H[i] * H[j] / (normH*normEpsH) + normB*(-H[i] * H[j] / (normEpsH2*normEpsH) +
					(i == j ? 1.0 / normEpsH : 0.0));
			}
		}
	}
	else {
		for (i = 0; i < 3; ++i) {
			for (j = 0; j < 3; ++j) {
				Jac[i * 3 + j] = (i == j ? dNormB : 0.0);
			}
		}
	}
}//interp1d


void keep_window_open(){
	std::cout << "Press ENTER to continue..." << std::flush;
	std::cin.ignore(std::numeric_limits<std::streamsize> ::max(), '\n');
}