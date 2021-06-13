#include "../sm/Elements/Beams/nlbeam_internalcontact.h"
#include "material.h"
#include "crosssection.h"
#include "node.h"

#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "engngm.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(NlBeamInternalContact);


  NlBeamInternalContact :: NlBeamInternalContact (int n, Domain *aDomain) : NlBeam_SM(n, aDomain)
{
}


// ====================================================
// AUXILIARY FUNCTIONS, TO BE REPLACED BY STANDARD ONES
// ====================================================

 IRResultType
NlBeamInternalContact :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    // first call parent
    NlBeam_SM :: initializeFrom(ir);    
    return IRRT_OK;
}


  
double
NlBeamInternalContact :: signum(double x)
{
  if (x>0.)
    return 1.;
  else if (x<0.)
    return -1.;
  return 0.;
}

double
NlBeamInternalContact :: L2norm(double x, double y, double z)
{
  return sqrt(x*x+y*y+z*z);
}

/*
Solution of a set of up to four linear equations by Gauss elimination.
In Matlab just x=A\b.
*/
void
NlBeamInternalContact :: solve(double A[4][4], double b[4], double x[4], int n)
{
 int i,j,k;
 for (i=0; i<n; i++) x[i] = b[i];
 // forward
 for (i=0; i<n-1; i++) for(j=i+1; j<n; j++){
  double aux = A[j][i]/A[i][i];
  for (k=i; k<n; k++) A[j][k] -= aux*A[i][k];
  x[j] -= aux*x[i];
  }
 // backward
 for (j=n-1; j>0; j--){
    x[j] /= A[j][j];
    for (i=0; i<j; i++) x[i] -= A[i][j]*x[j];
 }
 x[0] /= A[0][0];
}

/*
Solution of a set of three linear equations by Gauss elimination.
In Matlab just x=A\b.
*/
void
NlBeamInternalContact :: solve3(double A[3][3], double b[3], double x[3])
{
 int i,j,k;
 for (i=0; i<3; i++) x[i] = b[i];
 // forward
 for (i=0; i<2; i++) for(j=i+1; j<3; j++){
  double aux = A[j][i]/A[i][i];
  for (k=i; k<3; k++) A[j][k] -= aux*A[i][k];
  x[j] -= aux*x[i];
  }
 // backward
 for (j=2; j>0; j--){
    x[j] /= A[j][j];
    for (i=0; i<j; i++) x[i] -= A[i][j]*x[j];
 }
 x[0] /= A[0][0];
}

/*
Inversion of a 3x3 matrix.
*/
void
NlBeamInternalContact :: invert(double A[3][3], double Ainv[3][3])
{
  int i,j,k,l;
 for (i=0; i<3; i++){
   for (j=0; j<3; j++) Ainv[i][j] = 0.;
   Ainv[i][i] = 1.;
 }
 // forward
 for (i=0; i<2; i++) for(j=i+1;j<3;j++){
  double aux = A[j][i]/A[i][i];
  for (k=i; k<3; k++) A[j][k] -= aux*A[i][k];
  for (l=0; l<3; l++) Ainv[j][l] -= aux*Ainv[i][l];
  }
 // backward
 for (j=2; j>0; j--){
    for (l=0; l<3; l++) Ainv[j][l] /= A[j][j];
    for (i=0; i<j; i++) for (l=0; l<3; l++) Ainv[i][l] -= A[i][j]*Ainv[j][l];
 }
 for (l=0; l<3; l++) Ainv[0][l] /= A[0][0];
}

/*
Inversion of a 4x4 matrix.
*/
void
NlBeamInternalContact :: invert4(double A[4][4], double Ainv[4][4])
{
  int i,j,k,l;
 for (i=0; i<=3; i++){
   for (j=0; j<=3; j++) Ainv[i][j] = 0.;
   Ainv[i][i] = 1.;
 }
 // forward
 for (i=0; i<3; i++) for(j=i+1;j<=3;j++){
  double aux = A[j][i]/A[i][i];
  for (k=i; k<=3; k++) A[j][k] -= aux*A[i][k];
  for (l=0; l<=3; l++) Ainv[j][l] -= aux*Ainv[i][l];
  }
 // backward
 for (j=3; j>0; j--){
    for (l=0; l<=3; l++) Ainv[j][l] /= A[j][j];
    for (i=0; i<j; i++) for (l=0; l<=3; l++) Ainv[i][l] -= A[i][j]*Ainv[j][l];
 }
 for (l=0; l<=3; l++) Ainv[0][l] /= A[0][0];
}

// ============================================================================================
// FUNCTIONS RELATED TO THE MOMENT-CURVATURE RELATION (THE SIMPLEST LINEAR ONE IS ASSUMED HERE)
// ============================================================================================

double
NlBeamInternalContact :: computeCurvatureFromMoment(double M)
{
  return M/EI;
}

double
NlBeamInternalContact :: computeDerMomentFromCurvature(double kappa)
{
  return EI;
}

// ============================================================================================
// FUNCTIONS RELATED TO THE CONTACT MODEL WITH DRY COULOMB FRICTION
// ============================================================================================

/*
Evaluation of the contact loading function fc from the given normal force and shear force 
at the contact point. Negative values of fc indicate sticking and zero value indicates sliding.
The shear force is actually transmitted in the direction normal to the contact surface and
the normal force in the direction tangential to the contact surface. Whether positive shear
force means tension (inadmissible) or compression in the contact surface depends on which
segment is above and which one is below. Positive shear force means that the left segment
acts on the right segment upwards.

Input:
Nc ... normal force at contact point
Qc ... shear force at contact point
cmode ... contact mode

Return value:
fc ... contact loading function
*/   
double
NlBeamInternalContact :: evalContactLoadingFunction(double Nc, double Qc, contactModeType cmode)
{
  double fc;
  if (cmode==AA_cmode || cmode==AB_cmode || cmode==AC_cmode || cmode==CA_cmode) // right segment above the left one
      fc = fabs(Nc) - friction*Qc; 
  else // right segment below the left one
      fc = fabs(Nc) + friction*Qc; 
  return fc;
}

/*
Evaluation of the differential (linearized increment) of the contact loading function fc.

Input:
Nc ... normal force at contact point
dN ... increment of normal force at contact point
dQ ... increment of shear force at contact point
cmode ... contact mode

Return value:
dfc ... linearized increment of the contact loading function
*/   
double
NlBeamInternalContact :: evalDerContactLoadingFunction(double Nc, double dN, double dQ, contactModeType cmode)
{
  // contribution of normal force increment
  double dfc = signum(Nc)*dN;
  // contribution of shear force increment
  if (cmode==AA_cmode || cmode==AB_cmode || cmode==AC_cmode || cmode==CA_cmode) // right segment above the left one
      dfc -= friction*dQ; 
  else // right segment below the left one
      dfc += friction*dQ; 
  return dfc;
}

// ========================================================================
// KEY ALGORITHM - INTEGRATION ALONG A BEAM SEGMENT
// ========================================================================

/*
  Numerical integration of differential equations describing a beam segment.
  Left-end displacements set to zero, left-end forces specified as input.
  Computes displacements at the end of the segment of given length Lb.
  If an inflection point c is detected, computes its distance from the left end, Lc, and the displacements at c.
  Also computes two Jacobi matrices that contain derivatives of displacements at b or c with respect to the left-end forces
  and with respect to the spatial coordinate. 
  If the user specifies it, integration stops at the first inflection point and the results that refer to the segment end are then meaningless.
  Everything is done in local coordinates aligned with the left beam end.

  Input variables:
  * fab ... horizontal force, vertical force and moment at the left end
  * deltaPhi ... prescribed displacement jump at the inflexion point (useful for inverted modes, otherwise zero) 
  * Lb ... total length of the segment
  * inflection_only ... if this is set to "true", the integration stops as soon as an inflection point is found
  * outfile ... file into which the results should be printed

  Output variables:
  * ub ... horizontal displacement, vertical displacement and rotation at the end of the segment
  * jac_b ... derivatives of ub with respect to fab and x (evaluated at x=Lb)
  * Lc ... coordinate of inflection point
  * uc ... horizontal displacement, vertical displacement and rotation at the inflection point
  * jac_c ... derivatives of uc with respect to fab and x (evaluated at x=Lc)

  Return value:
  * true or false, indicating whether an inflection point has been detected (if not, the output Lc, uc and jac_c is meaningless)
*/
bool
NlBeamInternalContact :: integrateAlongSegment(double fab[3], double deltaPhi, double Lb, double ub[3], double jac_b[3][4], double* Lc, double uc[3], double jac_c[3][4], bool inflection_only, FILE* outfile)
{
  bool printflag = (outfile!=NULL);
  bool inflection_detected = false;
  int i, j, istep = 0;
  double aux, jacobi[3][3], jacobi_prev[3][3], jac_s[3][3];
  
  // initialization at the left end
  double Xab = fab[0], Zab = fab[1], Mab = fab[2];
  double M = -Mab;
  double dM[3];
  dM[0] = 0.; dM[1] = 0.; dM[2] = -1.;
  double kappa = computeCurvatureFromMoment(M);
  double dMdkappa = computeDerMomentFromCurvature(kappa);
  double dkappa[3];
  for (j=0; j<3; j++)
    dkappa[j] = dM[j]/dMdkappa;
  double x = 0., u[3], u_prev[3];
  for (i=0; i<3; i++){
    u[i] = 0.;
    for (j=0; j<3; j++)
      jacobi[i][j] = 0.;
  }
  if (printflag) fprintf(outfile,"%g %g %g %g\n",x,u[0],u[1],u[2]);

  // basic loop over spatial steps
  
  int nstep = ceil(Lb/DX); // the spatial step size is fixed and the segment length does not need to be its integer multiple
  for (istep=1; istep<=nstep; istep++){
    x += DX;
    for (i=0; i<3; i++)
      u_prev[i] = u[i];
    // rotation at midstep and its derivatives with respect to the left-end forces
    double phi_mid = u_prev[2]+kappa*DX/2.;
    double dphi_mid[3];
    for (j=0; j<3; j++)
      dphi_mid[j] = jacobi[2][j]+dkappa[j]*DX/2.;
    // normal force at midstep and its derivatives with respect to the left-end forces
    double N_mid = -Xab*cos(phi_mid)+Zab*sin(phi_mid);
    double dN_mid[3];
    for (j=0; j<3; j++)
      dN_mid[j] = (Xab*sin(phi_mid)+Zab*cos(phi_mid))*dphi_mid[j];
    dN_mid[0] -= cos(phi_mid); dN_mid[1] += sin(phi_mid);
    // horizontal displacement at the end of the step
    u[0] = u_prev[0]+DX*((1.+N_mid/EA)*cos(phi_mid)-1.);
   // vertical displacement at the end of the step
    u[1] = u_prev[1]-DX*(1.+N_mid/EA)*sin(phi_mid);
    // bending moment and curvature at the end of the step and their derivatives with respect to the left-end forces
    double M_prev = M; // (store the moment at the beginning of the step, needed for the inflection check)
    M = -Mab+Xab*u[1]-Zab*(x+u[0]);
    for (j=0; j<3; j++)
      dM[j] = Xab*jacobi[1][j]-Zab*jacobi[0][j];
    dM[0] += u[1];
    dM[1] += -(x+u[0]);
    dM[2] += -1.;
    kappa = computeCurvatureFromMoment(M);
    dMdkappa = computeDerMomentFromCurvature(kappa);
    for (j=0; j<3; j++)
      dkappa[j] = dM[j]/dMdkappa;
    // rotation at the end of the step
    u[2] = phi_mid+kappa*DX/2.;

    // test whether Jacobi at the beginning of the step needs to be stored
    if (istep==nstep || (!inflection_detected && M*M_prev<=0. && M_prev!=M)){
      for (i=0; i<3; i++)
	for (j=0; j<3; j++)
	  jacobi_prev[i][j] = jacobi[i][j];
    }

    // update Jacobi matrix
    for (j=0; j<3; j++){
      jacobi[0][j] += DX*(dN_mid[j]/EA)*cos(phi_mid) - DX*(1.+N_mid/EA)*sin(phi_mid)*dphi_mid[j];
      jacobi[1][j] -= DX*(dN_mid[j]/EA)*sin(phi_mid) + DX*(1.+N_mid/EA)*cos(phi_mid)*dphi_mid[j];
      jacobi[2][j] = dphi_mid[j]+dkappa[j]*DX/2.;
    }
     
    // test whether inflection occurs (the first occurence is considered)
    if (!inflection_detected && M*M_prev<=0. && M_prev!=M){
      // inflection point is detected - output parameters Lc, uc and jac_c will be computed
      inflection_detected = true;
      aux = M_prev / (M_prev-M);
      *Lc = x - DX + aux*DX;
      // displacements, rotation and their derivatives at Lc by linear interpolation within the current step
      for (i=0; i<3; i++){
	uc[i] = u_prev[i] + aux*(u[i]-u_prev[i]);
	for (j=0; j<3; j++)
	  jac_c[i][j] = aux*jacobi[i][j] + (1.-aux)*jacobi_prev[i][j];
      }
      // derivatives with respect to the spatial coordinate
      for (i=0; i<3; i++)
	  jac_c[i][3] = (u[i]-u_prev[i])/DX;
      // return if the user cares only about the inflection point and not about the segment end
      if (inflection_only)
	return true;
      // change the rotation by deltaPhi (useful for analysis of inverted contact modes)
      u[2] += deltaPhi;
    }

    // test whether the end of the segment has been reached
    if (istep==nstep){
      // displacements, rotation and their derivatives at Lb by linear interpolation within the last step
      aux = Lb/DX - (nstep-1);
      for (i=0; i<3; i++){
	ub[i] = u_prev[i] + aux*(u[i]-u_prev[i]);
	for (j=0; j<3; j++)
	  jac_b[i][j] = aux*jacobi[i][j] + (1.-aux)*jacobi_prev[i][j];
      }
      // derivatives with respect to the spatial coordinate
      for (i=0; i<3; i++)
	  jac_b[i][3] = (u[i]-u_prev[i])/DX;
      if (printflag) fprintf(outfile,"%g %g %g %g\n",Lb,ub[0],ub[1],ub[2]);
    } else {
      if (printflag) fprintf(outfile,"%g %g %g %g\n",x,u[0],u[1],u[2]);
    }
  } // end of loop over spatial steps
  
  return inflection_detected;
}

// ========================================================================
// AUXILIARY ALGORITHMS FOR THE BEAM ELEMENT WITH INTERNAL CONTACT
// ========================================================================

/*
  This method transforms the relative displacements and rotation of the right end with respect to the left end
  (with displacement components in coordinate system aligned with the left end) 
  into the relative displacements and rotation of the left end with respect to the right end
  (with displacement components in coordinate system aligned with the right end).

  Input: ub ... right-end relative displacements and rotation
  Output: ua ... left-end relative displacements and rotation
*/ 
void
NlBeamInternalContact :: transform_ub2ua(double ub[3], double ua[3])
{
  ua[0] = (beamLength+ub[0])*cos(ub[2]) - ub[1]*sin(ub[2]) - beamLength;
  ua[1] = (beamLength+ub[0])*sin(ub[2]) + ub[1]*cos(ub[2]);
  ua[2] = -ub[2];
}

/*
  This method computes the right-end forces and moment from the left-end forces and moment
  and from the relative displacements and rotation of the right end, based on equilibrium.
  
  Input: 
  ub ... right-end relative displacements and rotation
  fab ... left-end forces and moment
  
  Output: 
  fba ... right-end forces and moment
*/ 
void
NlBeamInternalContact :: transform_fab2fba(double ub[3], double fab[3], double fba[3])
{
  fba[0] = fab[0]*cos(ub[2]) - fab[1]*sin(ub[2]);
  fba[1] = fab[0]*sin(ub[2]) + fab[1]*cos(ub[2]);
  fba[2] = fab[0]*ub[1] - fab[1]*(beamLength+ub[0]) - fab[2];
}

/*
  This method computes the right-end forces and moment from the left-end forces and moment
  and from the relative displacements and rotation of the right end, based on equilibrium.

  Input: 
  ub ... right-end relative displacements and rotation
  fba ... right-end forces and moment
  
  Output: 
  fab ... left-end forces and moment
*/ 
void
NlBeamInternalContact :: transform_fba2fab(double ub[3], double fba[3], double fab[3])
{
  double ua[3];
  transform_ub2ua(ub, ua);
  transform_fab2fba(ua, fba, fab);
}

bool
NlBeamInternalContact :: checkLeftSegmentLengthAdmissibility(double Lc)
{
  return (Lc>=0 && Lc<=leftSegmentLength);
}

bool
NlBeamInternalContact :: checkRightSegmentLengthAdmissibility(double Lc)
{
  return (Lc>=0 && Lc<=rightSegmentLength);
}

/*
  For the given contact mode (which must be one of the tip modes),
  this method checks whether the given rotation jump at the contact point would be admissible.

  Input:
  deltaPhi ... rotation jump
  cmode ... contact mode (one of the tip modes)
  
  Return value:
  'true' if the rotation jump is admissible
*/
bool
NlBeamInternalContact :: checkRotationAdmissibility(double deltaPhi, contactModeType cmode)
{
  if (cmode==AC_cmode || cmode==CB_cmode) // right segment must rotate counterclockwise wrt to the left one
    return (sin(deltaPhi)>=0.); 
  else // right segment must rotate clockwise
    return (sin(deltaPhi)<=0.); 
}

/////////////////////////////////////////////////////////////////////
// TIP CONTACT
/////////////////////////////////////////////////////////////////////

/*
 Iterative search for the left-end forces and variables that characterize the contact behavior.
 Everything is done in local coordinates aligned with the left beam end.
 The tip of the right segment is assumed to be in contact with some point on the left segment.
 For modes AC and BC, this assumption is satisfied naturally.
 For modes CA and CB, it is necessary to swap left and right, which is done before the present method is invoked.
 It is specified whether "sticking" is assumed, in which case the position of the contact point is prescribed 
 and does not change, or "sliding" is assumed, in which case this position needs to be determined.
 Note that this function does NOT check admissibility of the solution. It only makes sure that
 three conditions are satisfied. Two of them are compatibility conditions (in the deformed
 configuration, the position of the contact point on the left segment coincides with the
 position of the tip of the right segment). The third condition is "zero moment at the (fixed)
 contact point" for the sticking mode and "zero contact loading function at the contact point"
 for the sliding mode. The loading function is computed from the normal and shear forces and
 its zero value indicates the critical combination that leads to frictional sliding.

 Input variables: 
 * ub_target ... displacements and rotation of the right end wrt to the fixed left end that should be achieved
 * fab ... initial guess of the left-end forces 
 * Lac ... initial active length of the left segment (can change if sliding occurs)
 * Lb ... total length of the right segment (is fixed, tip contact assumed)
 * cmode ... contact mode (one of the tip modes)
 * process ... type of process (0=sticking, 1=sliding)
 * printflag ... flag indicating whether detailed results should be printed

 Output variables: 
 * fab ... left-end forces at the end of the step
 * Lac ... final active length of the left segment (different from initial value if sliding occured)
 * Nca ... normal force in the left segment at the contact point (= inflection point)
 * Qca ... shear force in the left segment at the contact point
 * deltaPhi ... relative rotation of the right segment wrt to the left segment at the contact point
 
 Return value:
 * true or false, indicates whether a converged solution has been found (admissibility is not checked here)
*/
bool
NlBeamInternalContact :: findLeftEndForcesLocal_Tip_SoS(double ub_target[3], double fab[3], double* Lac, double Lb, contactModeType cmode, int process, double* Nca, double* Qca, double* deltaPhi, bool printflag)
{
  bool converged = false;
  double res[3], fba[3], uac[3], jacobi_ac[3][4], ubc[3], jacobi_bc[3][4], jacobi[3][3], dforces[3];
  double Mca, error;
  double phib = ub_target[2];
  double c = cos(phib);
  double s = sin(phib);
  // transformation matrix 
  double F[3][3];
  F[0][0] = c; F[0][1] = -s; F[0][2] = 0.;
  F[1][0] = s; F[1][1] =  c; F[1][2] = 0.;
  F[2][0] = ub_target[1]; F[2][1] = -beamLength-ub_target[0]; F[2][2] = -1.;
  int iter = 0, i, j;
  // weights for conversion to dimensionless quantities
  double weight_disp = 1./beamLength;
  double weight_mom = beamLength/EI; 
  double weight_force = beamLength*weight_mom; 

  while (iter++<MAXIT_BEAM){// main iterative loop of the Newton-Raphson method
    
    // integrate along the left segment
    double L_dummy, u_dummy[3], jac_dummy[3][4];
    if (process==0){ // sticking - we do not care about the inflection point here (dummy variables used)
      trialProcess = Stick_proc;
      integrateAlongSegment(fab, 0., *Lac, uac, jacobi_ac, &L_dummy, u_dummy, jac_dummy, false, NULL);
    } else { // sliding
      trialProcess = Slid_proc;
      double La = 10.*leftSegmentLength; // a very safe estimate - we only want to find the inflection point
      bool inflex = integrateAlongSegment(fab, 0., La, u_dummy, jac_dummy, Lac, uac, jacobi_ac, true, NULL);
      if (!inflex)
	return false; // no inflexion point found, contact condition cannot be enforced
      if (printflag) printf("left-end displacements: %g %g %g\n",uac[0],uac[1],uac[2]);
    }
    // internal forces at the contact point
    *Nca = -fab[0]*cos(uac[2])+fab[1]*sin(uac[2]);
    *Qca = -fab[0]*sin(uac[2])-fab[1]*cos(uac[2]);
    Mca = -fab[2]+fab[0]*uac[1]-fab[1]*(*Lac+uac[0]);
    // from equilibrium, we calculate the corresponding initial guess of the right-end forces
    transform_fab2fba(ub_target, fab, fba);
    if (printflag) printf("right end forces: %g %g %g\n",fba[0],fba[1],fba[2]);   
    // integrate along the right segment - we do not care about the inflection point here
    integrateAlongSegment(fba, 0., Lb, ubc, jacobi_bc, &L_dummy, u_dummy, jac_dummy, false, NULL);
    if (printflag) printf("right-end displacements: %g %g %g\n",ubc[0],ubc[1],ubc[2]);   
    // evaluate the residual
    // the first two components are always displacement differences
    res[0] = *Lac + uac[0] - beamLength - ub_target[0] + (Lb+ubc[0])*c + ubc[1]*s;
    res[1] = uac[1] - ub_target[1] - (Lb+ubc[0])*s + ubc[1]*c;
    // the third component depends on the type of process
    if (process==0){ // sticking
      res[2] = Mca;
      error = L2norm(weight_disp*res[0], weight_disp*res[1], weight_mom*res[2]);
    } else { // sliding
      res[2] = evalContactLoadingFunction(*Nca, *Qca, cmode);
      error = L2norm(weight_disp*res[0], weight_disp*res[1], weight_force*res[2]);
    }
    if (printflag) printf("residual: %g %g %g\nerror: %g\n",res[0],res[1],res[2],error);
    // check the convergence criterion
    if (error<TOL_BEAM){// converged solution found
      // evaluate the rotation jump at the contact point
      *deltaPhi = ub_target[2] + ubc[2] - uac[2];
      if (!stiffEvalMode) return true;
      converged = true; // return postponed, since Jacobi should be updated
    }
    
    // set up the Jacobi matrix (derivatives of the residual wrt to the left-end forces)

    // prepare auxiliary variables
    double du[3], dw[3], dN[3], dQ[3]; 
    for (j=0; j<3; j++){
      du[j] = 0.;
      dw[j] = 0.;
      for (i=0; i<3; i++){
	du[j] += jacobi_bc[0][i]*F[i][j];
	dw[j] += jacobi_bc[1][i]*F[i][j];
      }
    }
    if (process==1){ // sliding  
      // prepare derivatives of internal forces at the contact point
      for (j=0; j<3; j++){
	dN[j] = (fab[0]*sin(uac[2])+fab[1]*cos(uac[2])) * jacobi_ac[2][j];
 	dQ[j] = (-fab[0]*cos(uac[2])+fab[1]*sin(uac[2])) * jacobi_ac[2][j];
      }
      dN[0] += -cos(uac[2]);
      dN[1] += sin(uac[2]);
      dQ[0] += -sin(uac[2]);
      dQ[1] += -cos(uac[2]);
    }
    // evaluate the part of Jacobi matrix that has a regular structure (additional terms come later)
    for (j=0; j<3; j++){
      jacobi[0][j] = jacobi_ac[0][j] + du[j]*c + dw[j]*s;
      jacobi[1][j] = jacobi_ac[1][j] - du[j]*s + dw[j]*c;
      if (process==0) // sticking
	jacobi[2][j] = fab[0]*jacobi_ac[1][j]-fab[1]*jacobi_ac[0][j]; 
      else // sliding
	jacobi[2][j] = evalDerContactLoadingFunction(*Nca, dN[j], dQ[j], cmode);
    }
    if (process==0){ // additional terms for sticking (correction of the third row)
      jacobi[2][0] += uac[1];
      jacobi[2][1] += -(*Lac+uac[0]);
      jacobi[2][2] += -1.;
    } else { // additional terms for sliding (correction of the first two rows)
      // this piece of code takes into account the fact that changes of end forces lead to changes in the position of inflection point
      double dMdx = fab[0]*jacobi_ac[1][3]-fab[1]*(1.+jacobi_ac[0][3]); 
      double dMc[3], dLac[3];
      dMc[0] = uac[1];
      dMc[1] = -(*Lac+uac[0]);
      dMc[2] = -1.;
      for (j=0; j<3; j++){
	dMc[j] += fab[0]*jacobi_ac[1][j]-fab[1]*jacobi_ac[0][j];
	dLac[j] = -dMc[j]/dMdx; // derivatives of the position of inflection point wrt the left-end forces
      }
      for (j=0; j<3; j++){
	jacobi[0][j] += (1.+jacobi_ac[0][3])*dLac[j];
	jacobi[1][j] += jacobi_ac[1][3]*dLac[j];
      }  
    }

    if (converged && stiffEvalMode){ // prepare a block of the stiffness matrix
      invert(jacobi, Kblock);
      double c1 = (Lb+ubc[0])*s - ubc[1]*c; 
      double c2 = (Lb+ubc[0])*c + ubc[1]*s;
      // the last column is a linear combination of the first two columns
      for (i=0; i<3; i++)
	Kblock[i][2] = c1*Kblock[i][0] + c2*Kblock[i][1];
      return true;
    }
    
    // solve the linearized problem
    solve3(jacobi,res,dforces);
    if (printflag){
      printf("corrections of left-end forces: ");
      for (j=0; j<3; j++) {printf("%g ",-dforces[j]);}// note the negative sign
      printf("\n");
    }
    // correct the left-end forces
    for (i=0; i<3; i++)
      fab[i] -= dforces[i]; // note the negative sign
    if (printflag) printf("left end forces: %g %g %g\n",fab[0],fab[1],fab[2]);
  } // end of while loop
  
  // no convergence in MAXIT_BEAM iterations
  if (printflag){
    printf("No convergence in findLeftEndForcesLocal_Tip_SoS for process = \n");
    if (process==0)
      printf("sticking\n");
    else
      printf("sliding\n");
  }
  return false;
}

/* 
   Auxiliary function: the argument 'x' is shifted by an integer multiple of 2.*PI to get a number between -PI and PI.
*/
double
NlBeamInternalContact :: shiftToIntervalFromMinusPiToPi(double x)
{
  double answer = x - 2.*PI*floor((x+PI)/(2.*PI));
  return answer;
}

/*
  The most likely contact mode at the end of the step is determined based on the information
  on the results obtained with the assumption that the mode is equal to a given tip mode.
  The information provided is the rotation jump at the contact point and the result of
  segment length admissibility evaluation.

  Input:
  deltaPhi ... rotation jump at the contact point
  seglength_ok ... 'true' if the active length of both segments is between zero and the total segment length
  current_cmode ... the mode for which the solution has been computed (must be one of the tip modes)

  Return value:
  The contact mode that is expected to occur at the end of the step.
*/
contactModeType
NlBeamInternalContact :: suggestedMode(double deltaPhi, bool seglength_ok, contactModeType current_cmode)
{
  // if the tip in contact has slipped beyond the other tip, it is possible that the other tip gets activated
  if (!seglength_ok){
    switch (current_cmode){
    case AC_cmode: return CB_cmode;
    case BC_cmode: return CA_cmode;
    case CA_cmode: return BC_cmode;
    case CB_cmode: return AC_cmode;
    }
  }
  // now we know that the tip in contact has NOT slipped beyond the other tip
  // the rotation jump is first shifted into the interval from -PI to PI 
  deltaPhi = shiftToIntervalFromMinusPiToPi(deltaPhi);
  // if the rotation jump is admissible, the current_cmode should be the correct one,
  // otherwise we can expect a change into a smooth mode, depending on which boundary
  // of the admissible interval has been crossed
  switch (current_cmode){
  case AC_cmode:
    if (deltaPhi>=0.)
      return AC_cmode;
    else if (deltaPhi>-PI/2.)
      return AA_cmode;
    else
      return AB_cmode;
  case BC_cmode:
    if (deltaPhi<=0.)
      return BC_cmode;
    else if (deltaPhi<PI/2.)
      return BB_cmode;
    else
      return BA_cmode;
  case CA_cmode: 
    if (deltaPhi<=0.)
      return CA_cmode;
    else if (deltaPhi<PI/2.)
      return AA_cmode;
    else
      return BA_cmode;
  case CB_cmode: 
    if (deltaPhi>=0.)
      return CB_cmode;
    if (deltaPhi>-PI/2.)
      return BB_cmode;
    else
      return AB_cmode;
  }
  return current_cmode;
}

/*
 Evaluation of the left-end forces and moment for given relative displacements and rotation
 of the right end, assuming tip contact.
 Everything is done in local coordinates aligned with the left beam end.
 The tip of one segment is assumed to be in contact with some point on the other segment.
 For mode AC, the tip of the right segment touches the left segment from "above", for mode BC from "below".
 For mode CB, the tip of the left segment touches the right segment from "above", for mode CA from "below".
 First, "sticking" is assumed, and if this solution does not satisfy the loading condition
 at the contact point, then "sliding" is assumed.
 Admissibility of the solution is checked: 
 * The relative rotation angle at the contact point must be in the appropriate range, otherwise a mode change occurs. 
 * The active length of the segment must not exceed its total length, otherwise a mode change occurs. 
 The return value indicates whether the solution has been found and whether it is admissible. 

 Input variables: 
 * ub_target ... displacements and rotation of the right end wrt to the fixed left end that should be achieved
 * fab ... initial guess of the left-end forces 
 * cmode ... assumed contact mode (must be one of the tip modes)
 * printflag ... flag indicating whether detailed results should be printed

 Output variables: 
 * fab ... left-end forces at the end of the step

 Return value:
 * true or false, indicates whether a converged solution has been found
*/
contactModeType
NlBeamInternalContact :: findLeftEndForcesLocal_Tip(double ub_target[3], double fab[3], contactModeType cmode, bool printflag)
{
  double Nca, Qca, deltaPhi, fc, Lac;
  bool  success, rotation_ok;
  // store the initial guess for possible restart with sliding assumption
  double f_init[3], fba[3], ua_target[3];
  int i;
  
  bool tipIsOnRightSegment;
  switch(cmode){
  case AC_cmode: case BC_cmode:
    tipIsOnRightSegment = true;
    break;
  case CA_cmode: case CB_cmode:
    tipIsOnRightSegment = false;
    break;
  default:
    printf("findLeftEndForcesLocal_Tip cannot be used for cmode = %d\n",cmode);
    exit(1);
  }

  if (Process!=Slid_proc){// if it is known that the process is sliding, proceed directly to that case
  // sticking is assumed first
    
    trialProcess = Stick_proc;
    if (tipIsOnRightSegment){
      Lac = trialLeftActiveSegmentLength;
      for (i=0; i<3; i++)
	f_init[i] = fab[i]; // fab stored for possible restart with sliding
      success = findLeftEndForcesLocal_Tip_SoS(ub_target, fab, &Lac, rightSegmentLength, cmode, 0, &Nca, &Qca, &deltaPhi, printflag);  
    } else {
      Lac = trialRightActiveSegmentLength;
      // from equilibrium, we calculate the corresponding initial guess of the right-end forces
      transform_fab2fba(ub_target, fab, fba);
      for (i=0; i<3; i++)
	f_init[i] = fba[i]; // fba stored for possible restart with sliding
      transform_ub2ua(ub_target, ua_target);
      success = findLeftEndForcesLocal_Tip_SoS(ua_target, fba, &Lac, leftSegmentLength, cmode, 0, &Nca, &Qca, &deltaPhi, printflag);
      deltaPhi = -deltaPhi; // compensating for swapped segments
      transform_fba2fab(ub_target, fba, fab);
    }
    if (!success){
      if (printflag)
	printf("findLeftEndForcesLocal_Tip failed for sticking\n");
      // give sliding a chance - make sure that restart takes place with the original guess of fab or fba
      if (tipIsOnRightSegment)
	for (i=0; i<3; i++) fab[i] = f_init[i];
      else
	for (i=0; i<3; i++) fba[i] = f_init[i];
    } else { // a solution for sticking has been found but its admissibility still needs to be checked
      fc = evalContactLoadingFunction(Nca, Qca, cmode);
      rotation_ok = checkRotationAdmissibility(deltaPhi, cmode);
      if (printflag){
	printf("findLeftEndForcesLocal_Tip converged in sticking\n");
	printf(" Normal force: %g\n", Nca);
	printf(" Shear force: %g\n", Qca);
	printf(" Rotation jump: %g\n",deltaPhi);
	printf(" Loading function: %g\n",fc);
	if (fc>0.)
	  printf(" The solution violates the contact conditions\n");
	if (!rotation_ok)
	  printf(" The rotation is not in the admissible range\n");
	if (fc<=0. && rotation_ok)
	  printf("The solution is admissible\n");
	else
	  printf("The solution is not admissible\n");     
      }
      if (fc<=0. && rotation_ok){ // the solution with sticking assumed is admissible
	if (tipIsOnRightSegment){
	  trialRightActiveSegmentLength = rightSegmentLength;
	} else {
	  trialLeftActiveSegmentLength = leftSegmentLength;
	}
	return cmode;
      }
    }
  }

  // sliding is assumed
  trialProcess = Slid_proc;
  if (tipIsOnRightSegment)
    success = findLeftEndForcesLocal_Tip_SoS(ub_target, fab, &Lac, rightSegmentLength, cmode, 1, &Nca, &Qca, &deltaPhi, printflag);
  else {
    success = findLeftEndForcesLocal_Tip_SoS(ua_target, fba, &Lac, leftSegmentLength, cmode, 1, &Nca, &Qca, &deltaPhi, printflag);
    deltaPhi = -deltaPhi; // compensating for swapped segments
  }
  if (!success){
    if (printflag)
      printf("findLeftEndForcesLocal_Tip failed for sliding\n");
    trialLeftActiveSegmentLength = leftSegmentLength;
    trialRightActiveSegmentLength = rightSegmentLength;
    return N_cmode;
  }
  // a solution for sliding has been found but its admissibility still needs to be checked
  rotation_ok = checkRotationAdmissibility(deltaPhi, cmode);
  bool seglength_ok;
  if (tipIsOnRightSegment)
    seglength_ok = checkLeftSegmentLengthAdmissibility(Lac);
  else
    seglength_ok = checkRightSegmentLengthAdmissibility(Lac);
  if (printflag){
    fc = evalContactLoadingFunction(Nca, Qca, cmode);
    printf("findLeftEndForcesLocal_Tip converged in sliding\n");
    printf(" Normal force: %g\n", Nca);
    printf(" Shear force: %g\n", Qca);
    printf(" Rotation jump: %g\n",deltaPhi);
    printf(" Loading function: %g\n",fc);
    if (!rotation_ok)
      printf(" The rotation is not in the admissible range\n");
    if (!seglength_ok)
      printf(" The segment length after sliding is not in the admissible range\n");
    if (rotation_ok && seglength_ok)
      printf("The solution is admissible\n");
    else
     printf("The solution is not admissible\n");     
  }
  if (seglength_ok){ 
    if (tipIsOnRightSegment){
      trialLeftActiveSegmentLength = Lac;
      trialRightActiveSegmentLength = rightSegmentLength;
    } else {
      transform_fba2fab(ub_target, fba, fab); 
      trialLeftActiveSegmentLength = leftSegmentLength;
      trialRightActiveSegmentLength = Lac;
    }
    if (rotation_ok) return cmode;// the solution with sliding assumed is admissible
  } else {
    trialLeftActiveSegmentLength = leftSegmentLength;
    trialRightActiveSegmentLength = rightSegmentLength;
  }
  return suggestedMode(deltaPhi, seglength_ok, cmode);
}

/////////////////////////////////////////////////////////////////////
// SMOOTH CONTACT
/////////////////////////////////////////////////////////////////////

/*
Find forces and moment at the left end that lead to given displacements and rotations
at the right end, assuming that the displacement and rotation at the left end are zero
and that the process corresponds to ROLLING, i.e., the sum of the active segment lengths
remains constant during the step.
The solution is computed iteratively, using the Newton-Raphson method. 
Everything is done here in the local coordinate system aligned with the left end of the beam.
The contact point is located at the inflexion point and its distances from the left and right ends  
(i.e., the active segment lengths) are considered as input as well as output parameters. 
On input, what matters is only their sum. The deformed shape is computed by integrating
over the whole beam because the total length is known. The left-end forces are adjusted
such that the displacements and rotation at the right end are correct. However, if no inflexion
point is detected for the converged solution, the contact conditions cannot be satisfied
and the solution is not admissible. 

Input:
ub_target ... displacements and rotations of the right end with respect to the left end
fab ... initial guess of the left-end forces and moment
deltaPhi ... rotation jump at the contact point (zero for regular modes, +-PI for the inverted modes)
Lac ... initial value of the active length of the left segment
Lbc ... initial value of the active length of the right segment
printflag ... flag indicating whether detailed results should be printed

Output:
fab ... left-end forces and moment at the end of the step
Lac ... active length of the left segment at the end of the step
Lbc ... active length of the right segment at the end of the step
Nc ... normal force at the contact point
Qc ... shear force at the contact point

Return value:
indicates success or failure (failure means that the iterative process does not converge or the deformed
beam does not have an inflexion point) 
*/
bool
NlBeamInternalContact :: findLeftEndForcesLocal_Smooth_Rolling(double ub_target[3], double fab[3], double deltaPhi, double* Lac, double* Lbc, double* Nc, double* Qc, bool printflag)
{
  double res[3], dforces[3], ub_loc[3], jac[3][3], jac_b[3][4], uc[3], jac_c_dummy[3][4]; 
  double ub, wb, phib;
  int iter = 0, i, j;
  double L = *Lac + *Lbc;
  // weight for conversion to dimensionless quantities
  double weight_disp = 1./beamLength;
 
  while (iter++<MAXIT_BEAM){
    // integrate along beam of length L
     bool inflection_found = integrateAlongSegment(fab, deltaPhi, L, ub_loc, jac_b, Lac, uc, jac_c_dummy, false, NULL);
     for (i=0; i<3; i++) // the left 3x3 block of jac_b is stored in Jacobi, potentially used later for stiffness evaluation 
      for (j=0; j<3; j++)
	Jacobi[i][j] = jac_b[i][j];
    *Nc = -fab[0]*cos(uc[2])+fab[1]*sin(uc[2]);
    *Qc = -fab[0]*sin(uc[2])-fab[1]*cos(uc[2]);

    if (printflag){
      printf("\niteration %d\ndisplacements: %g %g %g\nJacobi matrix:\n",iter,ub_loc[0],ub_loc[1],ub_loc[2]);
      for (i=0; i<3; i++){
	for (j=0; j<3; j++) {printf("%g ",Jacobi[i][j]);}
	printf("\n");
      }
    }
    // evaluate the residual and check the convergence criterion
    for (i=0; i<3; i++)
      res[i] = ub_target[i]-ub_loc[i]; 
    // take into account the fact that, due to previous sliding,
    // the current beam length may differ from the initial distance between joints
    // (displacement ub_target[0] is taken with respect to the initial length of the beam)
    res[0] += beamLength - L;
    double error = L2norm(weight_disp*res[0], weight_disp*res[1], res[2]);
    if (printflag) printf("error: %g\n",error);
    if (error<TOL_BEAM){// converged
      if (!inflection_found)
	return false;// converged but no inflexion detected
      else{
	*Lbc = L - *Lac;
	return true;// converged and inflexion detected
      }
    }
    // compute the iterative correction of left-end forces
    for (i=0; i<3; i++)
      for (j=0; j<3; j++) // we make a copy of Jacobi because solve3 will destroy it
	jac[i][j] = Jacobi[i][j];
    solve3(jac,res,dforces);
    if (printflag){
      printf("corrections of left-end forces: ");
      for (j=0; j<3; j++) {printf("%g ",dforces[j]);}
      printf("\n");
    }
    for (i=0; i<3; i++)
      fab[i] += dforces[i]; 
    if (printflag) printf("left end forces: %g %g %g\n",fab[0],fab[1],fab[2]);
  }
  // no convergence in MAXIT_BEAM iterations
  if (printflag) printf("No convergence in findLeftEndForcesLocal (sticking)\n");
  return false;
}

/*
Find forces and moment at the left end that lead to given displacements and rotations
at the right end, assuming that the displacement and rotation at the left end are zero
and that the process corresponds to SLIDING, i.e., the sum of the active segment lengths
varies during the step.
Everything is done here in the local coordinate system aligned with the left end of the beam.
The solution is computed iteratively, using the Newton-Raphson method. 
In addition to the left-end forces and moment, the fourth basic unknown is the total beam length
(sum of the active segment lengths). In addition to the conditions dealing with prescribed
displacements and rotation at the right end, the fourth equation is the contact condition
at the contact point. Since sliding is assumed, the combination of the normal force and shear
force must be such that the contact loading function fc vanishes. 

The contact point is located at the inflexion point and its distances from the left and right ends  
(i.e., the active segment lengths) are considered as input as well as output parameters. 
On input, what matters is only their sum. The deformed shape is computed by integrating
over the whole beam because the total length is known. The left-end forces are adjusted
such that the displacements and rotation at the right end are correct. However, if no inflexion
point is detected for the converged solution, the contact conditions cannot be satisfied
and the solution is not admissible. 

Input:
ub_target ... displacements and rotations of the right end with respect to the left end
fab ... initial guess of the left-end forces and moment
deltaPhi ... rotation jump at the contact point (zero for regular modes, +-PI for the inverted modes)
Lac ... initial value of the active length of the left segment
Lbc ... initial value of the active length of the right segment
printflag ... flag indicating whether detailed results should be printed

Output:
fab ... left-end forces and moment at the end of the step
Lac ... active length of the left segment at the end of the step
Lbc ... active length of the right segment at the end of the step
Nc ... normal force at the contact point
Qc ... shear force at the contact point

Return value:
indicates success or failure (failure means that the iterative process does not converge or the deformed
beam does not have an inflexion point) 
*/

bool
NlBeamInternalContact :: findLeftEndForcesLocal_Smooth_Sliding(double ub_target[3], double fab[3], double deltaPhi, double* Lac, double* Lbc, contactModeType cmode, bool printflag)
{
  bool converged = false;
  double ub, wb, phib, res[4], dforces[4], ub_loc[3], jac_b[3][4], uc[3], jac_c[3][4], jac[4][4];
  double Nc, Qc, dN[3], dQ[3];
  int iter = 0, i, j;
  double L = *Lac + *Lbc;
  // weights for conversion to dimensionless quantities
  double weight_disp = 1./beamLength;
  double weight_force = beamLength*beamLength/EI; 

  while (iter++<MAXIT_BEAM){
    bool inflection_found = integrateAlongSegment(fab, deltaPhi, L, ub_loc, jac_b, Lac, uc, jac_c, false, NULL); 
    if (!inflection_found) return false; // no inflexion point - contact condition cannot be enforced
    // internal forces at the contact point
    Nc = -fab[0]*cos(uc[2])+fab[1]*sin(uc[2]);
    Qc = -fab[0]*sin(uc[2])-fab[1]*cos(uc[2]);
    // evaluate the residual and check the convergence criterion
    for (i=0; i<3; i++)
      res[i] = ub_loc[i]-ub_target[i];
    res[3] = evalContactLoadingFunction(Nc, Qc, cmode);
    // take into account the fact that, due to sliding,
    // the current beam length may differ from the initial distance between joints
    // (displacement ub_target[0] is taken with respect to the initial length of the beam)
    res[0] += L - beamLength;
    double error = L2norm(weight_disp*sqrt(res[0]*res[0]+res[1]*res[1]), res[2], weight_force*res[3]);
    if (printflag) printf("error: %g\n",error);
    if (error<TOL_BEAM)
      converged = true; // the actual return is a bit postponed, to make update of Jacobi possible if needed

    // prepare the Jacobi matrix for the set of 4 equations
    if (!converged || stiffEvalMode){
      // evaluate the first three rows
      for (i=0; i<3; i++)
	for (j=0; j<4; j++)
	  Jacobi44[i][j] = jac_b[i][j];
      // correct one entry by a term that reflects the line res[0] += L - beamLength;
      Jacobi44[0][3] += 1.;
      // derivatives of internal forces at the contact point
      for (j=0; j<3; j++){
	dN[j] = (fab[0]*sin(uc[2])+fab[1]*cos(uc[2])) * jac_c[2][j];
	dQ[j] = (-fab[0]*cos(uc[2])+fab[1]*sin(uc[2])) * jac_c[2][j];
      }
      dN[0] += -cos(uc[2]);
      dN[1] += sin(uc[2]);
      dQ[0] += -sin(uc[2]);
      dQ[1] += -cos(uc[2]);
      // evaluate the fourth row of Jacobi matrix 
      for (j=0; j<3; j++)
	Jacobi44[3][j] = evalDerContactLoadingFunction(Nc, dN[j], dQ[j], cmode);
      Jacobi44[3][3] = 0.; // the total length does not affect the inflection point
      if (printflag){
	printf("\niteration %d\ndisplacements: %g %g %g\nJacobi matrix:\n",iter,ub_loc[0],ub_loc[1],ub_loc[2]);
	for (i=0; i<4; i++){
	  for (j=0; j<4; j++) printf("%g ",Jacobi44[i][j]);
	  printf("\n");
	}
      } 
    }
    // this is the postponed return after convergence
    if (converged){
      *Lbc = L - *Lac;
      return true;
    }
   // compute the iterative correction of left-end forces and of active beam length
    for (i=0; i<4; i++)
      for (j=0; j<4; j++)
	jac[i][j] = Jacobi44[i][j];

    solve(jac,res,dforces,4);
    if (printflag){
      printf("corrections of left-end forces and of the length: ");
      for (j=0; j<4; j++) printf("%g ",-dforces[j]);
    }
    for (i=0; i<3; i++)
      fab[i] -= dforces[i];
    L -= dforces[3];
    if (printflag) printf("left-end forces and length: %g %g %g %g\n",fab[0],fab[1],fab[2],L);
  }
  // no convergence in MAXIT_BEAM iterations
  if (printflag) printf("No convergence in findLeftEndForcesLocal (sliding)\n");
  return false;
}

/*
Find forces and moment at the left end and the length of the beam that lead to given displacements and rotations
at the right end, assuming that the displacement and rotation at the left end are zero, and to the satisfaction
of the conditions describing the sliding unit inside the beam.
The solution is computed iteratively, using the Newton-Raphson method. 
Everything is done here in the local coordinate system aligned with the left end of the beam.
The contact point is located at the inflexion point and its distances from the left and right ends  
(i.e., the active segment lengths) are considered as input as well as output parameters. 

Input:
ub_target ... displacements and rotations of the right end with respect to the left end
fab ... initial guess of the left-end forces and moment
cmode ... assumed contact mode (must be one of the smooth ones)
printflag ... flag indicating whether detailed results should be printed

Output:
fab ... left-end forces and moment at the end of the step

Return value:
the same as the assumed mode (input variable cmode) if the solution is admissible, otherwise a guess of the likely contact mode
*/
contactModeType
NlBeamInternalContact :: findLeftEndForcesLocal_Smooth(double ub_target[3], double fab[3], contactModeType cmode, bool printflag)
{
  double Nc, Qc;
  // store the initial guess for possible restart with sliding assumption
  double fab_init[3], deltaPhi, Lac, Lbc;
  bool success;
  int i;
  for (i=0; i<3; i++) fab_init[i] = fab[i];

  if (Process!=Slid_proc){// if it is known that the process is sliding, proceed directly to that case

  // first assume that no additional sliding occurs during this step (i.e., fix Lac+Lbc and iterate on fab)
    trialProcess = Roll_proc;  
  deltaPhi = 0.;
  Lac = trialLeftActiveSegmentLength;
  Lbc = trialRightActiveSegmentLength;
  if (cmode==AB_cmode)
    deltaPhi = PI;
  else if (cmode==BA_cmode)
    deltaPhi = -PI;
  bool success = findLeftEndForcesLocal_Smooth_Rolling(ub_target, fab, deltaPhi, &Lac, &Lbc, &Nc, &Qc, printflag);
  // the result is 'false' if the left-end forces do not converge, or if they converge but the beam has no inflection point 
  if (!success){
    if (printflag) printf("findLeftEndForcesLocal_Smooth failed for sticking\n");
    // give sliding a chance - make sure that the solution that assumes sliding will later start from the original guess of left-end forces 
    for (i=0; i<3; i++) fab[i] = fab_init[i];
  } else { // sticking solution has been found, now let us see whether it is admissible
    double fc = evalContactLoadingFunction(Nc, Qc, cmode);
    if (printflag){
      printf("findLeftEndForcesLocal_Smooth converged in sticking\n");
      printf(" Normal force: %g\n", Nc);
      printf(" Shear force: %g\n", Qc);
      printf(" Left active segment length: %g\n",Lac);
      printf(" Right active segment length: %g\n",Lbc);
      printf(" Loading function: %g\n",fc);
      if (fc>0.)
	printf(" The solution violates the contact conditions\n");
      else
	printf("The solution is admissible\n");     
    }
    if (fc<=0.){ // admissible solution
      trialLeftActiveSegmentLength = Lac;
      trialRightActiveSegmentLength = Lbc;
      return cmode;
    }
  }
  }

  // now assume that sliding occurs during this step (i.e., iterate on fab as well as on Lac+Lbc)
  trialProcess = Slid_proc;  
  Lac = trialLeftActiveSegmentLength;
  Lbc = trialRightActiveSegmentLength;
  success = findLeftEndForcesLocal_Smooth_Sliding(ub_target, fab, deltaPhi, &Lac, &Lbc, cmode, printflag);
  if (success){ // the active segment lengths change due to sliding
    if ((Lac+Lbc)>leftSegmentLength+rightSegmentLength)
      return N_cmode; // most likely jumps to no contact
    else if (Lac>leftSegmentLength){
      trialLeftActiveSegmentLength = leftSegmentLength;
      trialRightActiveSegmentLength = Lbc;
      if (cmode==AA_cmode)
	return CA_cmode; // most likely changes into tip contact
      else
	return CB_cmode;
    }
    else if (Lbc>rightSegmentLength){
      trialLeftActiveSegmentLength = Lac;
      trialRightActiveSegmentLength = rightSegmentLength;
      if (cmode==AA_cmode)
	return AC_cmode; // most likely changes into tip contact 
      else
	return BC_cmode;
    }
    // admissible solution
    trialLeftActiveSegmentLength = Lac;
    trialRightActiveSegmentLength = Lbc;
    return cmode;
  }
  else
    return N_cmode; // most likely smoothly changes into no contact
}

/*
  This function checks whether the no-contact mode is possible for the given relative displacements.
  If it is, the return value is N_cmode.
  Otherwise the function tries to guess the most likely contact mode but the result is not always
  correct because the behavior is path-dependent.
  The function is used in situation when the mode at the beginning of the step is NOT the no-contact mode.
  Based on the behavior (e.g., lack of convergence or inadmissible solutions for assumed contact modes),
  it is concluded that contact is probably lost, but this assumption needs to be verified, which is
  the purpose of the present function. 
*/
contactModeType
NlBeamInternalContact :: predictContactMode(double ub[3])
{
  // assuming that both segments are straight, find their intersection
  
  double phib = ub[2];
  if (sin(phib)==0.) // if both segments are parallel, the no-contact mode is always possible
    return N_cmode;
  double Lbc = -ub[1]/sin(phib);
  if (Lbc<0. || Lbc>rightSegmentLength) return N_cmode;
  double Lac = beamLength+ub[0]-Lbc*cos(phib);
  if (Lac<0. || Lac>leftSegmentLength) return N_cmode;
  if (rightSegmentLength-Lbc < leftSegmentLength-Lac){ // right tip probably touches the left segment first 
    if (ub[1]>0.) return BC_cmode;
    else return AC_cmode;
  }
  // left tip probably touches the right segment first
  if (ub[1]>0.) return CA_cmode;
  return CB_cmode;
}

/*
  This method finds the time at which the tip of a segment hits the x-axis.
  It is assumed that the opposite end of the segment is initially located at (u_prev[0],u_prev[1])
  and rotated by u_prev[2], and that its displacement and rotation linearly increases during the step.
  The 'time' is a dimensionless parameter running from 0 at the beginning to 1 at the end of the step.
  Linear increase of rotation of the opposite end means that the tip travels along a curved trajectory.
  A nonlinear equation describing the condition that the current value of the second coordinate vanishes
  is solved iteratively. 

  Input variables:
  u_prev[3] ... 2 coordinates and rotation at the end of the previous step,
  du[3] ... increments of 2 coordinates and rotation during the step

  Return value:
  dimensionless time at which the tip passes through the x-axis
  (if no converged solution is found, the return value is -1., which is outside the admissible interval) 
 */
double
NlBeamInternalContact :: findTipContactTime(double u_prev[3], double du[3], double segLength)
{
  double t;
  // simple linearized estimate, used as the initial guess for the iterative scheme
  double db = du[1]+segLength*(sin(u_prev[2]+du[2])-sin(u_prev[2])); 
  if (db!=0.){
    t = (-u_prev[1]-segLength*sin(u_prev[2]))/db;
  } else {
    t = 0.5;
  }
  // iterative solution of nonlinear equation
  for (int iter=1; iter<=MAXIT_CONTACT_TIME; iter++){ 
      double f = u_prev[1] + t*du[1] + segLength*sin(u_prev[2]+t*du[2]);
      if (fabs(f)<TOL_CONTACT_TIME*segLength) 
	return t;// converged solution
      double fp = du[1] + segLength*cos(u_prev[2]+t*du[2])*du[2];
      if (fp!=0.)
	t -= f/fp;
  }
    // converged solution not found
  t = -1.;
}

/*
  This function checks whether, during the step from ub_prev to ub, contact would occur. 
  If it does, the return value indicates which contact mode would arise, 
  otherwise the result is the N_cmode.
  It is assumed that the state at ub_prev is a no-contact state.
  For contact to occur during the step, the tip of one segment must pass through the other segment.
  The tip that hits the other segment first determines the contact mode.

  Input variables:
  ub[3] ... relative displacements and rotation at the end of the current step
  ub_prev[3] ... relative displacements and rotation at the end of the previous step,

  Return value:
  most likely type of contact mode at the end of the step
*/
contactModeType
NlBeamInternalContact :: predictContactMode(double ub[3], double ub_prev[3])
{
  // assuming that both segments are straight,
  //        check at which time 'tb' and location 'xb' the right tip hits the left segment
  // (the 'time' is a dimensionless parameter running from 0 at the beginning to 1 at the end of the step)
  int i;
  double xa, xb;
  double du[3];
  for (i=0; i<3; i++)
    du[i] = ub[i] - ub_prev[i];
  double tb = findTipContactTime(ub_prev, du, rightSegmentLength);
  if (tb>=0. && tb<=1.) // find location of the hit event
    xb = beamLength + ub_prev[0] + tb*du[0] - rightSegmentLength*cos(ub_prev[2]+tb*du[2]);
   
  // now check at which time 'ta' and location 'xa' the left tip hits the right segment
  double ua[3], ua_prev[3];
  transform_ub2ua(ub, ua);
  transform_ub2ua(ub_prev, ua_prev);
  for (i=0; i<3; i++)
    du[i] = ua[i] - ua_prev[i];
  double ta = findTipContactTime(ua_prev, du, leftSegmentLength);
  if (ta>=0. && ta<=1.) // find location of the hit event
    xa = beamLength + ua_prev[0] + ta*du[0] - leftSegmentLength*cos(ua_prev[2]+ta*du[2]);

  // determine which hit event occurs at an admissible time and location
  bool left_hits = (ta>=0. && ta<=1. && xa>=0. && xa<=rightSegmentLength);
  bool right_hits = (tb>=0. && tb<=1. && xb>=0. && xb<=leftSegmentLength);

  if (!left_hits && !right_hits) return N_cmode;
  if (left_hits && right_hits)
    left_hits = ta<tb;
  if (left_hits){
    trialLeftActiveSegmentLength = leftSegmentLength;
    trialRightActiveSegmentLength = xa;
    if (ua[1]-ua_prev[1]+leftSegmentLength*(sin(ua[2])-sin(ua_prev[2])) > 0.)
      return CA_cmode;
    else 
      return CB_cmode;
  }
  trialLeftActiveSegmentLength = xb;
  trialRightActiveSegmentLength = rightSegmentLength;
  if (ub[1]-ub_prev[1]+rightSegmentLength*(sin(ub[2])-sin(ub_prev[2])) > 0.)
    return AC_cmode;
  else
    return BC_cmode;
}

/*
This method attempts to find (iteratively) the left-end forces and moment that lead to
the given right-end displacements and rotation, taken as relative with respect to the left end.
Everything is done in the local coordinate system attached to the left end.

Input variables:
ub[3] ... relative displacements and rotation at the end of the current step
ub_prev[3] ... relative displacements and rotation at the end of the previous step,
fab[3] ...  initial guess of the left-end forces and moment
printflag ... flag indicating whether detailed results should be printed

Output variables:
fab[3] ...  left-end forces and moment (values of initial guess are rewritten)

Return value:
boolean, indicating whether an admissible solution has been found

The method first assumes that the contact mode remains the
same as in the previous step and attempts to compute the
corresponding solution. If this attempt fails, the assumption
is modified and the whole process is repeated until an admissible
solution is found or until the number of attempts exceeds
a given limit (currently set to 5). If an admissible solution
is found, the key internal variables (the contact mode and the 
active lengths of both segments) are updated, the left-end forces 
at the end of the step are passed as an output variable, and the method
returns Boolean value 'true' to indicate success.

If the assumed mode is no-contact, it is not necessary to solve
a special set of equations, just to check whether the cantilevers
would intersect if they remain straight. If an intersection is 
found, the most likely mode is guessed (based on certain rules),
otherwise the no-contact solution is accepted, which of course
implies that the end forces and element stiffness are set to zero.
The check for intersection of two straight cantilevers is performed 
using the predictContactMode method, which has two different versions,
depending on whether the contact mode at the beginning of the step
is or is not the no-contact mode.
*/
bool
NlBeamInternalContact :: findLeftEndForcesLocal(double ub[3], double ub_prev[3], double fab[3], bool printflag)
{
  contactModeType assumed_cmode = contactMode;
  trialLeftActiveSegmentLength = leftActiveSegmentLength;
  trialRightActiveSegmentLength = rightActiveSegmentLength;

  for (int attempt=1; attempt<=MAXIT_CONTACT_MODE; attempt++){
    switch (assumed_cmode){
      // smooth modes
    case AA_cmode:
    case BB_cmode:
    case AB_cmode:
    case BA_cmode:
      trialContactMode = findLeftEndForcesLocal_Smooth(ub, fab, assumed_cmode, printflag);
      break;
      // tip modes
    case AC_cmode:
    case BC_cmode:
    case CA_cmode:
    case CB_cmode:
      trialContactMode = findLeftEndForcesLocal_Tip(ub, fab, assumed_cmode, printflag);
      break;
      // no-contact mode
    case N_cmode:
      fab[0] = fab[1] = fab[2] = 0.;
      if (contactMode!=N_cmode)
	trialContactMode = predictContactMode(ub);
      else 
	trialContactMode = predictContactMode(ub, ub_prev);
      break;
      // this should never happen
    default:
      printf("Mode %d has not been implemented yet in findLeftEndForcesLocal\n",assumed_cmode);
      exit(1);
    }
    if (trialContactMode==assumed_cmode){// solution accepted, update internal variables
      Process = trialProcess;
      return true;
    }
    // solution not accepted, let us try another one with modified assumption which mode occurs
    assumed_cmode = trialContactMode;
  }
  // no admissible solution found after MAXIT_CONTACT_MODE attempts
  return false;
}


/*
This method sets up the 3x3 stiffness matrix that contains derivatives of the left-end
forces and moment with respect to the right-end displacements and rotation taken relative 
to the left-end displacements and rotation. Everything is done in the local coordinate system 
attached to the left end.
This matrix is later transformed into the full 6x6 element stiffness matrix in global coordinates. 

Input variables:
ub[3] ... relative displacements and rotation at the end of the current step
ub_prev[3] ... relative displacements and rotation at the end of the previous step,
fab[3] ...  a guess (or converged value) of the left-end forces and moment
printflag ... flag indicating whether detailed results should be printed

Output variables:
fab[3] ...  left-end forces and moment (values of initial guess are rewritten)
Kloc[3][3] ... a block of the stiffness matrix in local coordinates

Return value:
boolean, indicating whether everything has run smoothly
*/

//bool evaluateLocalStiffness(double ub[3], double ub_prev[3], double fab[3], double Kloc[3][3], bool printflag)
//{
  /*
In most cases, the method will be called when the values of fab have already been computed.
However, to make sure that these are really the correct values, they are first recomputed.
If the input contains the correct converged values, the process will now converge after
one iteration. Otherwise, the full iterative process is used and it may not converge.
  */

void
NlBeamInternalContact :: construct_T(double T[3][3], double phia)
{
  T[0][0] = T[1][1] = cos(alpha-phia);
  T[0][1] = sin(alpha-phia);
  T[1][0] = -T[0][1];
  T[2][2] = 1.;
  T[0][2] = T[1][2] = T[2][0] = T[2][1] = 0.;
}

/*
Auxiliary matrix for transformation of stiffness.
It is obtained by differentiating T with respect to phia.
*/
void
NlBeamInternalContact :: construct_Tprime(double T[3][3], double phia)
{
  T[0][0] = T[1][1] = sin(alpha-phia);
  T[0][1] = -cos(alpha-phia);
  T[1][0] = -T[0][1];
  T[2][2] = 0.;
  T[0][2] = T [1][2] = T[2][0] = T[2][1] = 0.;
}

void
NlBeamInternalContact :: construct_l(double l[3], double phia)
{
  l[0] = beamLength * (cos(phia)-1.);
  l[1] = beamLength * sin(phia);
  l[2] = 0.;
}

void
NlBeamInternalContact :: construct_l_IC(double l[3], double phia, double L)
{
  l[0] = beamLength * cos(phia) - L;
  l[1] = beamLength * sin(phia);
  l[2] = 0.;
}

void
NlBeamInternalContact :: construct_lprime(double l[3], double phia)
{
  l[0] = -beamLength * sin(phia);
  l[1] = beamLength * cos(phia);
  l[2] = 0.;
}

/*
Find forces and moment at the left end that lead to given displacements and rotations
at both ends.
The input displacements and output forces are taken with respect to the global coordinate system.
Note that the transformation matrix T is affected by angle alpha that specifies the initial beam geometry. 
*/
bool
NlBeamInternalContact :: findLeftEndForces(double u[6], double u_prev[6], double fab[3])
{
  double ub_loc[3], ub_prev_loc[3], T[3][3], T_prev[3][3], fab_loc[3];
  int i, j;
  double phia = u[2];
  double phia_prev = u_prev[2];
 
  // compute displacements of the right end with respect to the auxiliary coordinate system
  construct_l(ub_loc, phia);
  construct_l(ub_prev_loc, phia_prev);
  construct_T(T, phia); 
  construct_T(T_prev, phia_prev); 
  for (i=0; i<=2; i++)
    for (j=0; j<=2; j++){
      ub_loc[i] += T[i][j]*(u[3+j]-u[j]);
      ub_prev_loc[i] += T_prev[i][j]*(u_prev[3+j]-u_prev[j]);
    }

  // transform initial guess to local coordinates
  for (i=0; i<=2; i++){
    fab_loc[i] = 0.;
    for (j=0; j<=2; j++)
      fab_loc[i] += T[i][j]*fab[j];
  }

  // compute the corresponding left-end forces and moment in the auxiliary coordinate system
  if(!findLeftEndForcesLocal(ub_loc, ub_prev_loc, fab_loc, false))
    return false;
 
  // transform forces to global coordinates
  for (i=0; i<=2; i++){
    fab[i] = 0.;
    for (j=0; j<=2; j++)
      fab[i] += T[j][i]*fab_loc[j];
  }
  return true;
}


void
NlBeamInternalContact :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{

  // solution vector
  answer.resize(6);
  FloatArray vu, vui;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, vu);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Incremental, tStep, vui);
  
  
  double u[6], u_prev[6], f[6];
  for(int i = 1; i <= 6; i++) {
    u[i-1] = vu.at(i);
    u_prev[i-1]  = vu.at(i) - vui.at(i);
  }
  this->findEndForces(u, u_prev, f);
  for(int i = 1; i <= 6; i++) {
    answer.at(i) = f[i-1];  
  }

}



/*
Find end forces and moments that lead to given displacements and rotations
at both ends.
The input displacements and output forces are taken with respect to the global coordinate system.
*/
bool
NlBeamInternalContact :: findEndForces(double u[6], double u_prev[6], double f[6])
{
/* When evaluating the end forces "from scratch", we set the process status to unknown.
The algorithm always tries sticking first and if the solution is not admissible or not found,
it proceeds to the sliding assumption. This is always the case when we call findEndForces.
However, it is also possible to call directly findLeftEndForces and then, if the process
status has been previously detected and stored, the algorithm skips the sticking/rolling part
in cases when it is know that the actual process is sliding. This is used when the stiffness
is evaluated after previous evaluation of internal forces for the same iterated increment.
*/
  Process = trialProcess = Unknown_proc;
  findLeftEndForces(u, u_prev, f); // only the first three entries of f are computed
  f[3] = -f[0];
  f[4] = -f[1];
  double c1 = beamLength*sin(alpha) + u[4] - u[1];
  double c2 = -beamLength*cos(alpha) - u[3] + u[0];
  f[5] = c1*f[0] + c2*f[1] - f[2];
}

/*
Evaluate the tangent stiffness matrix based on given displacements at both ends
and given end forces at the left end (must be provided).
It is assumed that findLeftEndForces has been run first and that the 3x3 Jacobi matrix
has been stored in Jacobi[3][3] or the 4x4 Jacobi matrix has been stored in Jacobi44[4][4].

Input variables:
u[6] ... end displacements and rotations at the end of the current step, in global coordinates
u_prev[6] ... end displacements and rotations at the end of the previous step, in global coordinates
fab[3] ...  a guess (or converged value) of the left-end forces and moment, in local coordinates

Output variables:
fab[3] ...  left-end forces and moment (values of initial guess are rewritten)
K[6][6] ... element stiffness matrix in global coordinates

Return value:
boolean, indicating whether everything has run smoothly
*/

void
NlBeamInternalContact :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
  FloatArray vu, vui;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, vu);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Incremental, tStep, vui);
  
  
  double u[6], u_prev[6], f[6];
  for(int i = 1; i <= 6; i++) {
    u[i-1] = vu.at(i);
    u_prev[i-1]  = vu.at(i) - vui.at(i);
  }
  double K[6][6], fab[3];
  this->evalStiffnessMatrix(u, u_prev, fab,K);
  for(int i = 1; i <=6; i++) {
    for(int j = 1; j <=6; j++) {
      answer.at(i,j) = K[i-1][j-1];
    }
  }
    
}

bool
NlBeamInternalContact :: evalStiffnessMatrix(double u[6], double u_prev[6], double fab[3], double K[6][6])
{
  int i, j, k;
  double L, L_prev;
  double T[3][3], Tprime[3][3], Ginv[3][3], Ginv44[4][4], TtGinv[3][3], lprime[3], fab_loc[3], ub_loc[3], Tu[3];

  stiffEvalMode = true;
  bool success = findLeftEndForces(u, u_prev, fab);
  stiffEvalMode = false;
  if (!success)
    return false;

  // compute auxiliary matrices
  double phia = u[2];
  construct_T(T, phia); 
  construct_Tprime(Tprime, phia); 
  construct_lprime(lprime, phia);

  // transform left-end forces to local coordinates
  for (i=0; i<=2; i++){
    fab_loc[i] = 0.;
    for(j=0; j<=2; j++)
      fab_loc[i] += T[i][j]*fab[j];
  }

  // get Jacobi matrix in local coordinates and invert it
  switch (trialContactMode){
    // smooth modes
  case AA_cmode:
  case BB_cmode:
  case AB_cmode:
  case BA_cmode:
    L_prev = leftActiveSegmentLength + rightActiveSegmentLength;
    L = trialLeftActiveSegmentLength + trialRightActiveSegmentLength;
    if (fabs(L-L_prev)>TOL_LENGTH*beamLength){ // sliding
      invert4(Jacobi44, Ginv44);
      for (i=0; i<3; i++)
	for (j=0; j<3; j++)
	  Ginv[i][j] = Ginv44[i][j];
    } else { // rolling
      invert(Jacobi, Ginv);
    }
    break;
    // tip modes
  case AC_cmode:
  case BC_cmode:
      for (i=0; i<3; i++)
	for (j=0; j<3; j++)
	  Ginv[i][j] = Kblock[i][j];
      break;
    // no-contact mode
  case N_cmode:
    for (i=0; i<6; i++)
      for (j=0; j<6; j++)
	K[i][j] = 0.;
    return true;
    // this should never happen
  default:
    printf("Mode %d has not been implemented yet in evalStiffnessMatrix\n",trialContactMode);
    exit(1);
  }

  // compute product Ttransposed*Ginverse
  for (i=0; i<=2; i++)
    for(j=0; j<=2; j++){
      double aux = 0.;
      for(k=0; k<=2; k++) aux += T[k][i]*Ginv[k][j];
      TtGinv[i][j] = aux;
    }
  
  // compute product Ttransposed*Ginverse*T and store it in the upper stiffness block 
  for (i=0; i<=2; i++)
    for(j=0; j<=2; j++){
      double aux = 0.;
      for(k=0; k<=2; k++) aux += TtGinv[i][k]*T[k][j];
      K[i][j] = -aux;
      K[i][3+j] = aux;
    }

  // compute Tprime*(ub-ua)
  for (i=0; i<=2; i++){
    Tu[i] = 0.;
    for(j=0; j<=2; j++)
      Tu[i] += Tprime[i][j]*(u[3+j]-u[j]);
  }

  // compute additional stiffness terms in the third column (associated with phia)
  for (i=0; i<=2; i++)
    for(j=0; j<=2; j++)
      K[i][2] += TtGinv[i][j]*(Tu[j]+lprime[j]) + Tprime[j][i]*fab_loc[j];

  // copy minus first row into fourth row and minus second row into fifth row
  for(j=0; j<=5; j++){
    K[3][j] = -K[0][j];
    K[4][j] = -K[1][j];
  }

  // construct the sixth row (complete formula)
  double c1 = beamLength*sin(alpha) + u[4] - u[1];
  double c2 = -beamLength*cos(alpha) - u[3] + u[0];
  K[5][0] = fab[1];
  K[5][1] = -fab[0];
  K[5][2] = 0.;
  K[5][3] = -fab[1];
  K[5][4] = fab[0];
  K[5][5] = 0.;
  for(j=0; j<=5; j++)
    K[5][j] += c1*K[0][j] + c2*K[1][j] - K[2][j];

  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// All functions below serve only for the purpose of testing and plotting the results.
// They will not be needed by OOFEM.
/////////////////////////////////////////////////////////////////////////////////////////////

/*
  An auxiliary function for plotting of the deformed segment shape.
  The Jacobi matrix is not evaluated.
*/
void NlBeamInternalContact ::  integrateAlongSegmentAndPlot(double fab[3], double Lb, double segmentLength, double u0[2], double T[2][2], FILE* outfile)
{
  if (outfile==NULL) return;
  double ub[3], Lc, uc[3], uplot[2];
  bool inflection_detected = false;
  int i, j, istep = 0;
  double aux;
  
  // initialization at the left end
  double Xab = fab[0], Zab = fab[1], Mab = fab[2];
  double M = -Mab;
  double kappa = computeCurvatureFromMoment(M);
  double x = 0., u[3], u_prev[3];
  for (i=0; i<3; i++)
    u[i] = 0.;
  
  for (i=0; i<2; i++){
    uplot[i] = u0[i] + T[i][0]*(x+u[0]) + T[i][1]*u[1];
  }
  fprintf(outfile,"%g %g \n",uplot[0],uplot[1]);

  // basic loop over spatial steps
  
  int nstep = ceil(Lb/DX); // the spatial step size is fixed and the segment length does not need to be its integer multiple
  for (istep=1; istep<=nstep; istep++){
    x += DX;
    for (i=0; i<3; i++)
      u_prev[i] = u[i];
    // rotation at midstep 
    double phi_mid = u_prev[2]+kappa*DX/2.;
    // normal force at midstep 
    double N_mid = -Xab*cos(phi_mid)+Zab*sin(phi_mid);
    // horizontal displacement at the end of the step
    u[0] = u_prev[0]+DX*((1.+N_mid/EA)*cos(phi_mid)-1.);
   // vertical displacement at the end of the step
    u[1] = u_prev[1]-DX*(1.+N_mid/EA)*sin(phi_mid);
    // bending moment and curvature at the end of the step 
    double M_prev = M; // (store the moment at the beginning of the step, needed for the inflection check)
    M = -Mab+Xab*u[1]-Zab*(x+u[0]);
    kappa = computeCurvatureFromMoment(M);
    // rotation at the end of the step
    u[2] = phi_mid+kappa*DX/2.;

    // test whether inflection occurs (the first occurence is considered)
    if (!inflection_detected && M*M_prev<=0. && M_prev!=M){
      // inflection point is detected - output parameter Lc, uc will be computed
      inflection_detected = true;
      //aux = M_prev / (M_prev-M);
      //*Lc = x - DX + aux*DX;
    }

    // test whether the end of the segment has been reached
    if (istep==nstep){
      // displacements and rotation at Lb by linear interpolation within the last step
      aux = Lb/DX - (nstep-1);
      for (i=0; i<3; i++)
	ub[i] = u_prev[i] + aux*(u[i]-u_prev[i]);
  
      for (i=0; i<2; i++)
	uplot[i] = u0[i] + T[i][0]*(Lb+ub[0]) + T[i][1]*ub[1];   
      fprintf(outfile,"%g %g \n",uplot[0],uplot[1]);
      // plot also the straight segment behind the contact point, if it exists
      double Lstraight = segmentLength - Lb;
      if (Lstraight>0.){
	for (i=0; i<2; i++)
	  uplot[i] = u0[i] + T[i][0]*(Lb+ub[0]+Lstraight*cos(ub[2])) + T[i][1]*(ub[1]-Lstraight*sin(ub[2]));   
	fprintf(outfile,"%g %g \n",uplot[0],uplot[1]);
      }
    } else {
      for (i=0; i<2; i++)
	uplot[i] = u0[i] + T[i][0]*(x+u[0]) + T[i][1]*u[1];   
      fprintf(outfile,"%g %g \n",uplot[0],uplot[1]);
    }
  } // end of loop over spatial steps
  fprintf(outfile,"\n");
}

void NlBeamInternalContact :: plotSegment(double fab[3], double ub[3], bool isLeftSegment, FILE* outfile)
{
  // transformation vector and matrix for conversion of displacements
  double u[2], T[2][2], f[3];
  if (isLeftSegment){
    u[0] = u[1] = 0.;
    T[0][0] = T[1][1] = 1.;
    T[0][1] = T[1][0] = 0.;
  } else {
    u[0] = beamLength+ub[0]; u[1] = ub[1];
    double alphab = PI + ub[2];
    double c = cos(alphab);
    double s = sin(alphab);
    T[0][0] = T[1][1] = c;
    T[0][1] = s; T[1][0] = -s;
    // transformation of forces and moment
    transform_fab2fba(ub, fab, f);
  }
  if (isLeftSegment)
    integrateAlongSegmentAndPlot(fab, leftActiveSegmentLength, leftSegmentLength, u, T, outfile);
  else
    integrateAlongSegmentAndPlot(f, rightActiveSegmentLength, rightSegmentLength, u, T, outfile);
}

#define MAXSTAGE 10

void
NlBeamInternalContact :: plotResponse(int nstage, int nstep[MAXSTAGE], double ustep[3][MAXSTAGE], int iplot[MAXSTAGE])
{
  double ub[3], ub_prev[3], fab[3];
  ub[0] = ub[1] = ub[2] = 0.;
  fab[0] = fab[1] = fab[2] = 0.;
  
  double L = beamLength;
  int i, istage, istep, jplot=0;

  for (istage=0; istage<nstage; istage++){
    for (istep=0; istep<nstep[istage]; istep++){
      for (i=0; i<3; i++){
	ub_prev[i] = ub[i];
	ub[i] += ustep[i][istage];
      }
      
      bool success = findLeftEndForcesLocal(ub, ub_prev, fab, false);
      if (success){
	contactMode = trialContactMode;
	leftActiveSegmentLength = trialLeftActiveSegmentLength;
	rightActiveSegmentLength = trialRightActiveSegmentLength;
      } else 
	fab[0] = fab[1] = fab[2] = 0.;
      printf("%g %g %g %g %g %g %d %g %g\n",ub[0],ub[1],ub[2],fab[0],fab[1],fab[2],contactMode,leftActiveSegmentLength,rightActiveSegmentLength);
      if ((istep+1)%iplot[istage]==0){
	jplot++;
	//FILE* outfile_left = open_file(jplot, true);
	//FILE* outfile_right = open_file(jplot, false);
	//plotSegment(fab, ub, true, outfile_left); // plot the left segment
	//plotSegment(fab, ub, false, outfile_right); // plot the right segment
	//fclose(outfile_left);
 	//fclose(outfile_right);
      }
    }
  }
}

}
