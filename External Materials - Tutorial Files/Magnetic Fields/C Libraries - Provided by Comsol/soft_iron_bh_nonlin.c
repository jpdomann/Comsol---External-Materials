/** Interface to a magnetic saturation model. */

/** You are allowed to use, modify, and publish this External Material File and your modifications of it subject
*   to the terms and conditions of the COMSOL Software License Agreement (www.comsol.com/sla). */

/** Copyright © 2015 by COMSOL. */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#ifdef _MSC_VER								//Check if code is being compiled with visual studio 
#define EXPORT __declspec(dllexport)		//Define command to export symbol to a .dll
#else
#define EXPORT
#endif

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
				double *states) {  // States, nStates-vector
  int i, j, pos;
  const double eps = 3e-16;
  double normH2, normH, normEpsH2, normEpsH, normB, dNormB, s;
  // Interpolation data
  const int nInt = 17;
  double intH[17] = {0.0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3, 20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310, 5.00E+05};
  double intB[17] = {0.0, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.6283};


  // Check inputs
  if (nPar[0]!=0)
    return 1;	              // error code 1 = "Wrong number of parameters"
  if (nStates[0]!=0)
	return 2;                 // error code 2 = "Wrong number of states"

  // Compute input norm
  normH2 = H[0]*H[0]+H[1]*H[1]+H[2]*H[2];
  normEpsH2 = normH2+eps;
  normH = sqrt(normH2);
  normEpsH = sqrt(normEpsH2);
  // Find the interpolation interval 
  pos = 0;
  while (pos<nInt && intH[pos]<normH) {
    ++pos;
  }
  // Handle extrapolation
  if (pos==0) {
    pos = 1;
  } else if (pos==nInt) {
	pos = nInt-1;
  }
  // Interpolate
  s = (normH-intH[pos-1])/(intH[pos]-intH[pos-1]);
  normB = s*intB[pos]+(1.0-s)*intB[pos-1];
  dNormB = (intB[pos]-intB[pos-1])/(intH[pos]-intH[pos-1]);
  // Write B vector
  if (normH>0) {
    for (i=0; i<3; ++i) {
      B[i] = normB*H[i]/normEpsH;
    }
  } else {
    for (i=0; i<3; ++i) {
	  B[i] = dNormB*H[i];
	}
  }
  // Write Jacobian of B with respect to H
  if (normH>0) {
    for (i=0; i<3; ++i) {
      for (j=0; j<3; ++j) {
	    Jac[i*3+j] = dNormB*H[i]*H[j]/(normH*normEpsH)+normB*(-H[i]*H[j]/(normEpsH2*normEpsH)+(i==j ? 1.0/normEpsH : 0.0));
      }
    }
  } else {
    for (i=0; i<3; ++i) {
      for (j=0; j<3; ++j) {
	    Jac[i*3+j] = (i==j ? dNormB : 0.0);
      }
    }  
  }
  return 0;  // Return value 0 is success, any other value triggers an exception
}
