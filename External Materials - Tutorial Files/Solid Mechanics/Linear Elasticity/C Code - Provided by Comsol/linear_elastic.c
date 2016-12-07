/** Interface to an isotropic linear elastic solid with two parameters E and nu. 
*   Example code implements linear elastic behaviour. */

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

EXPORT int eval(double e[6],         // Input: Green-Lagrange strain tensor components in Voigt order (xx,yy,zz,yz,zx,xy)
                double s[6],         // Output: Second Piola-Kirchhoff stress components in Voigt order (xx,yy,zz,yz,zx,xy)
                double D[6][6],         // Output: Jacobian of stress with respect to strain, 6-by-6 matrix in row-major order
                int *nPar,         // Input: Number of material model parameters, scalar
                double *par,       // Input: Parameters: par[0] = E, par[1] = nu
                int *nStates,      // Input: Number of states, scalar        
                double *states) {  // States, nStates-vector
				
  int i,j;
  double E, nu;
 
  // Check inputs
  if (nPar[0] != 2)        // only two parameters needed, E and nu
    return 1;              // error code 1 = "Wrong number of parameters"
  if (nStates[0] !=0)      // simple linear elastic material, no states needed
  return 2;                // error code 2 = "Wrong number of states" 
 
  // Read input parameters from parameter vector      
  E = par[0];
  nu = par[1];

  // Check values of input parameters
  if (E <= 0.0) return -1;
  if (nu >= 0.5 || nu <= -1.0) return -1;
  
  // Set up Jacobian matrix
  // 6x6 Elasticity matrix D based on E and nu, initially filled by zeros  
  for (i = 0; i < 6; i++){
    for (j = 0; j < 6; j++) {
        D[i][j] = 0.0;
    }
  }
  D[0][0] = D[1][1] = D[2][2] = E * (1.0 - nu) / (1.0 + nu)/(1.0 - 2.0 * nu);                        // upper diagonal
  D[3][3] = D[4][4] = D[5][5] = E / (1.0 + nu);                                                      // lower diagonal, note that this equals 2G
  D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu);// off diagonal, equals lambda Lame
  
  // Compute Hooke's law: s = D*e
  for (i = 0; i < 6; i++){
    s[i] = 0.0;
    for (j = 0; j < 6; j++) {
      s[i] += D[i][j]*e[j];
    }
  }
  
  return 0;  // Return value 0 if success, any other value trigger an exception
}