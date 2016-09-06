// ==========================================================================
// $Id: Rectify.C $
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2010-, UBC
// 
// Creator: bradleyd ()
// Email:   bradleyd@cs.ubc.ca
// ==========================================================================

#ifndef CAMERACALIB_RECTIFY_H
#define CAMERACALIB_RECTIFY_H

#include "Rectify.hh"
#include "Camera.hh"
#include "MDA/LinearAlgebra/LinAlg.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


  /******************************************************************************
   * Some functions are not implemented anywhere that I can find them, 
   * but I need them here, so I'll put them here for now.
   ******************************************************************************/
  void matInverse3x3(Matrix M, Matrix& Mi)
  {
    double a,b,c,d,e,f,g,h,i;
    a = M[0][0]; b=M[0][1]; c=M[0][2];
    d = M[1][0]; e=M[1][1]; f=M[1][2];
    g = M[2][0]; h=M[2][1]; i=M[2][2];
    
    double det = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
    
    Mi[0][0] = e*i-f*h;
    Mi[0][1] = c*h-b*i;
    Mi[0][2] = b*f-c*e;
    Mi[1][0] = f*g-d*i;
    Mi[1][1] = a*i-c*g;
    Mi[1][2] = c*d-a*f;
    Mi[2][0] = d*h-e*g;
    Mi[2][1] = b*g-a*h;
    Mi[2][2] = a*e-b*d;
    
    for (int row=0;row<3;row++)
      for(int col=0; col<3; col++)
	Mi[row][col]/=det;
  }

  void crossProduct(const double& v1x, const double& v1y, const double& v1z,
		    const double& v2x, const double& v2y, const double& v2z,
		    double *vx, double *vy, double *vz)
  {
    *vx = (v1y*v2z - v1z*v2y);
    *vy = (v1z*v2x - v1x*v2z);
    *vz = (v1x*v2y - v1y*v2x);
  }

  /******************************************************************************
   * end of functions that I needed
   ******************************************************************************/


  /** Internal helper function to create rectifying transformations
   *  (converted from RectifKitv2.2 rectify.m matlab code) */
  void createRectTrans(Matrix K1, Matrix R1, Vector t1, Matrix K2, Matrix R2, Vector t2, Matrix& Rect1, Matrix& Rect2, Vector& d1, Vector& d2)
  {
    if (fabs(d1[1] - d2[1]) > 0.000001)
      {
	printf("Rectifying Error: left and right vertical displacements must be the same\n");
	return;
      }

    // build projection matrices
    Matrix P1(3,4);
    Matrix P2(3,4);
    Matrix Rt1(3,4);
    Matrix Rt2(3,4);
    for (int i=0; i<3; i++)
      {
	for (int j=0; j<3; j++)
	  {
	    Rt1[i][j] = R1[i][j];
	    Rt2[i][j] = R2[i][j];
	  }
	Rt1[i][3] = t1[i];
	Rt2[i][3] = t2[i];
      }
    multMatrixMatrix(K1, Rt1, P1);
    multMatrixMatrix(K2, Rt2, P2);
  
    // optical centers (unchanged)
    Vector c1(3);
    Vector c2(3);
    Matrix R1t = R1.getTranspose();
    Matrix R2t = R2.getTranspose();
    Matrix K1inv(3,3); // K1 inverse
    Matrix K2inv(3,3); // K2 inverse
    matInverse3x3(K1, K1inv);
    matInverse3x3(K2, K2inv);
    for (int i=0; i<3; i++)
      {
	c1[i] = P1[i][3];
	c2[i] = P2[i][3];
      }

    Vector c1tmp(3);
    Vector c2tmp(3);
    K1inv.rightMultiply(c1, c1tmp);
    c1=c1tmp;
    R1t.rightMultiply(c1, c1tmp);
    c1=c1tmp;
    K2inv.rightMultiply(c2, c2tmp);
    c2=c2tmp;
    R2t.rightMultiply(c2, c2tmp);
    c2=c2tmp;
    for (int i=0; i<3; i++)
      {
	c1[i] *= -1;
	c2[i] *= -1;
      }

    // new x axis (baseline from c1 to c2)
    Vector v1(3);
    for (int i=0; i<3; i++)
      {
	v1[i] = c2[i]-c1[i];
      }

    // new y axis (orthogonal to old z and new x)
    Vector v2(3);
    crossProduct(R1[2][0], R1[2][1], R1[2][2], v1[0], v1[1], v1[2], &(v2[0]), &(v2[1]), &(v2[2]));

    // new z axes (no choice, orthogonal to baseline and y)
    Vector v3(3);
    crossProduct(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], &(v3[0]), &(v3[1]), &(v3[2]));

    // new extrinsic (translation unchanged)
    v1.normalize();
    v2.normalize();
    v3.normalize();

    Matrix R(3,3);
    for(int i=0; i<3; i++)
      {
	R[0][i] = v1[i];
	R[1][i] = v2[i];
	R[2][i] = v3[i];
      }

    // new intrinsic (arbitrarily choose K2)
    Matrix Kn1 = K2;
    Matrix Kn2 = K2;
    Kn1[0][1] = 0;
    Kn2[0][1] = 0;

    // translate image centers
    Kn1[0][2] = Kn1[0][2] + d1[0];
    Kn1[1][2] = Kn1[1][2] + d1[1];
    Kn2[0][2] = Kn2[0][2] + d2[0];
    Kn2[1][2] = Kn2[1][2] + d2[1];

    // new projection matrices
    Matrix Rtn1(3,4);
    Matrix Rtn2(3,4);
    R.rightMultiply(c1,c1tmp);
    c1=c1tmp;
    R.rightMultiply(c2,c2tmp);
    c2=c2tmp;
    for (int i=0; i<3; i++)
      {
	c1[i] *= -1;
	c2[i] *= -1;
      }
    for (int i=0; i<3; i++)
      {
	for (int j=0; j<3; j++)
	  {
	    Rtn1[i][j] = R[i][j];
	    Rtn2[i][j] = R[i][j];
	  }
	Rtn1[i][3] = c1[i];
	Rtn2[i][3] = c2[i];
      }
    Matrix pl(3,4);
    Matrix pr(3,4);
    multMatrixMatrix(Kn1, Rtn1, pl);
    multMatrixMatrix(Kn2, Rtn2, pr);

    // rectifying image transformation
    Matrix Pn1(3,3); Matrix Po1(3,3); Matrix Po1inv(3,3);
    Matrix Pn2(3,3); Matrix Po2(3,3); Matrix Po2inv(3,3);
    for (int i=0; i<3; i++)
      {
	for (int j=0; j<3; j++)
	  {
	    Pn1[i][j] = pl[i][j];
	    Pn2[i][j] = pr[i][j];
	    Po1[i][j] = P1[i][j];
	    Po2[i][j] = P2[i][j];
	  }
      }
    matInverse3x3(Po1, Po1inv);
    matInverse3x3(Po2, Po2inv);
    multMatrixMatrix(Pn1, Po1inv, Rect1);
    multMatrixMatrix(Pn2, Po2inv, Rect2);
}
  


  /** compute the 3x3 rectifying transformations for a pair of cameras */
  void getRectifyTransforms(const Camera C1, const Camera C2, Matrix& Rect1, Matrix& Rect2)
  {
    // create rectifying transformations without centering    
    Vector centerL(2,0.0);
    Vector centerR(2,0.0);
    createRectTrans(C1.K, C1.R, C1.T, C2.K, C2.R, C2.T, Rect1, Rect2, centerL, centerR); 
    
    // center left
    Vector tmp(3), tmp2(3);
    tmp[0] = C1.K[0][2];
    tmp[1] = C1.K[1][2];
    tmp[2] = 1;
    Rect1.rightMultiply(tmp,tmp2);
    tmp=tmp2;
    tmp[0]/=tmp[2];
    tmp[1]/=tmp[2];
    centerL[0] = C1.K[0][2] - tmp[0];
    centerL[1] = C1.K[1][2] - tmp[1];

    // center right 
    tmp[0] = C2.K[0][2];
    tmp[1] = C2.K[1][2];
    tmp[2] = 1;
    Rect2.rightMultiply(tmp,tmp2);
    tmp=tmp2;
    tmp[0]/=tmp[2];
    tmp[1]/=tmp[2];
    centerR[0] = C2.K[0][2] - tmp[0];
    centerR[1] = C2.K[1][2] - tmp[1];

    // vertical displacements must be the same
    centerL[1] = centerR[1];

    // create final rectifying transformations, with centering
    createRectTrans(C1.K, C1.R, C1.T, C2.K, C2.R, C2.T, Rect1, Rect2, centerL, centerR);

  }
  
} /* namespace */

#endif /* CAMERACALIB_RECTIFY_H */
