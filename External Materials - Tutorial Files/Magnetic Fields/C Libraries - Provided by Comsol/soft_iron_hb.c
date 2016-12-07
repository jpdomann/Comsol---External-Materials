/** Interface to a magnetic saturation model. */

/** You are allowed to use, modify, and publish this External Material File and your modifications of it subject
*   to the terms and conditions of the COMSOL Software License Agreement (www.comsol.com/sla). */

/** Copyright © 2015 by COMSOL. */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#ifdef _MSC_VER
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

/** Interface to a general H(B) relation with a variable number of material
*   model parameters and a variable number of states. The example code implements
*   an interpolation function representative of soft iron.*/
EXPORT int eval(double *oldB,      // Magnetic flux density, previous converged step, 3-vector, input
                double *B,         // Magnetic flux density, 3-vector, input
                double *H,         // Magnetic field, 3-vector, state
				double *Jac,       // Jacobian of H with respect to B, 3-by-3 matrix in row-major order, output
				int *nPar,         // Number of material model parameters, scalar, input
                double *par,       // Material model parameters, nPar-vector, input
				int *nStates,      // Number of states, scalar, input
				double *states) {  // States, nStates-vector
  int i, j, pos;
  const double eps = 3e-16;
  double normB2, normB, normEpsB2, normEpsB, normH, dNormH, s;
  // Interpolation data
  const int nInt = 16;
  double intB[16] = {0.0, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};
  double intH[16] = {0.0, 663.146, 1067.5, 1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3, 20462.8, 32169.6, 61213.4, 111408, 175070, 261469, 318310};
  
  // Check inputs
  if (nPar[0]!=0)
    return 1;	              // error code 1 = "Wrong number of parameters"
  if (nStates[0]!=0)
	return 2;                 // error code 2 = "Wrong number of states"

  // Compute input norm
  normB2 = B[0]*B[0]+B[1]*B[1]+B[2]*B[2];
  normEpsB2 = normB2+eps;
  normB = sqrt(normB2);
  normEpsB = sqrt(normEpsB2);
  // Find the interpolation interval 
  pos = 0;
  while (pos<nInt && intB[pos]<normB) {
    ++pos;
  }
  // Handle extrapolation
  if (pos==0) {
    pos = 1;
  } else if (pos==nInt) {
	pos = nInt-1;
  }
  // Interpolate
  s = (normB-intB[pos-1])/(intB[pos]-intB[pos-1]);
  normH = s*intH[pos]+(1.0-s)*intH[pos-1];
  dNormH = (intH[pos]-intH[pos-1])/(intB[pos]-intB[pos-1]);
  // Write H vector
  if (normB>0) {
    for (i=0; i<3; ++i) {
      H[i] = normH*B[i]/normEpsB;
    }
  } else {
    for (i=0; i<3; ++i) {
	  H[i] = dNormH*B[i];
	}
  }
  // Write Jacobian of H with respect to B
  if (normB>0) {
    for (i=0; i<3; ++i) {
      for (j=0; j<3; ++j) {
	    Jac[i*3+j] = dNormH*B[i]*B[j]/(normB*normEpsB)+normH*(-B[i]*B[j]/(normEpsB2*normEpsB)+(i==j ? 1.0/normEpsB : 0.0));
      }
    }
  } else {
    for (i=0; i<3; ++i) {
      for (j=0; j<3; ++j) {
	    Jac[i*3+j] = (i==j ? dNormH : 0.0);
      }
    }  
  }
  return 0;  // Return value 0 is success, any other value triggers an exception
}
