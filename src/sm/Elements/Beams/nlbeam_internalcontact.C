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


  NlBeamInternalContact :: NlBeamInternalContact (int n, Domain *aDomain) : NlBeam_SM(n, aDomain), Jacobi(3,3), Jacobi44(4,4), Kblock(3,3)
{
  this->internalForces.resize(3);
}


 IRResultType
NlBeamInternalContact :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    // first call parent
    NlBeam_SM :: initializeFrom(ir);


    IR_GIVE_FIELD(ir, leftSegmentLength, _IFT_NlBeamInternalContact_LeftSegmentLength );
    IR_GIVE_FIELD(ir, rightSegmentLength, _IFT_NlBeamInternalContact_RightSegmentLength );
    IR_GIVE_FIELD(ir, leftActiveSegmentLength, _IFT_NlBeamInternalContact_LeftActiveSegmentLength );
    IR_GIVE_FIELD(ir, rightActiveSegmentLength, _IFT_NlBeamInternalContact_RightActiveSegmentLength );
    int cmode;
    IR_GIVE_FIELD(ir, cmode, _IFT_NlBeamInternalContact_ContactMode);
    contactMode = (ContactModeType) cmode;

    IR_GIVE_FIELD(ir, DX, _IFT_NlBeamInternalContact_dx );
    IR_GIVE_FIELD(ir, friction, _IFT_NlBeamInternalContact_Friction );


    auto nodeA  = this->giveNode(1);
    auto nodeB  = this->giveNode(2);
    FloatArray pointA({nodeA->giveCoordinate(1), nodeA->giveCoordinate(3)});
    FloatArray pointB({nodeB->giveCoordinate(1), nodeB->giveCoordinate(3)});
    
    cosAlpha = (pointB.at(1) - pointA.at(1))/beamLength;
    sinAlpha = (pointB.at(2) - pointA.at(2))/beamLength;
    
    return IRRT_OK;


}



double
NlBeamInternalContact :: L2norm(double x, double y, double z)
{
  return sqrt(x*x+y*y+z*z);
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

double
NlBeamInternalContact :: evalContactLoadingFunction(double Nc, double Qc, ContactModeType cmode)
{
  double fc;  
  if (cmode == CMT_AA || cmode == CMT_AB || cmode == CMT_AC || cmode==CMT_CA) {
    // right segment above the left one
      fc = fabs(Nc) - friction*Qc; 
  } else {
    // right segment below the left one
      fc = fabs(Nc) + friction*Qc;
  }
  return fc;
}

double
NlBeamInternalContact :: evalDerContactLoadingFunction(double Nc, double dN, double dQ, ContactModeType cmode)
{
  // contribution of normal force increment
  double dfc = signum(Nc)*dN;
  // contribution of shear force increment
  if (cmode==CMT_AA || cmode==CMT_AB || cmode==CMT_AC || cmode==CMT_CA) {
    // right segment above the left one
    dfc -= friction*dQ; 
  } else { // right segment below the left one
      dfc += friction*dQ;
  }
  return dfc;
}


bool
NlBeamInternalContact :: integrateAlongSegment(FloatArray &fab, double deltaPhi, double Lb, FloatArray &ub, FloatMatrix &jac_b, double &Lc, FloatArray &uc, FloatMatrix &jac_c, bool inflection_only)
{
  
  bool inflection_detected = false;
  FloatMatrix jacobi(3,3), jacobi_prev(3,3), jac_s(3,3);
  // initialization at the left end
  double Xab = fab.at(1);
  double Zab = fab.at(2);
  double Mab = fab.at(3);
  double M = -Mab;
  FloatArray dM(3);
  dM.zero();
  dM.at(3) = -1;
  double kappa = computeCurvatureFromMoment(M);
  double dMdkappa = computeDerMomentFromCurvature(kappa);
  FloatArray dkappa;
  dkappa = dM;
  dkappa.times(1./dMdkappa);
  FloatArray u(3), u_prev;
  u_prev = ub;
 
  // basic loop over spatial steps  
  int nstep = ceil(Lb/DX); // the spatial step size is fixed and the segment length does not need to be its integer multiple
  double x = 0;
  for (int istep = 1; istep <= nstep; istep++){
    x += DX;
    ///@todo: check the following line
    u_prev = u;
    // rotation at midstep and its derivatives with respect to the left-end forces
    double phi_mid = u_prev.at(3) + kappa * DX / 2.;
    FloatArray dphi_mid(3);
    for (int j = 1; j <= 3; j++)
      dphi_mid.at(j) = jacobi.at(3,j)+dkappa.at(j)*DX/2.;
    // normal force at midstep and its derivatives with respect to the left-end forces
    double N_mid = -Xab*cos(phi_mid)+Zab*sin(phi_mid);
    FloatArray dN_mid(3);
    dN_mid = dphi_mid;
    dN_mid.times(Xab*sin(phi_mid)+Zab*cos(phi_mid));
    dN_mid.at(1) -= cos(phi_mid);
    dN_mid.at(2) += sin(phi_mid);
    // horizontal displacement at the end of the step
    u.at(1) = u_prev.at(1) + DX * ((1.+N_mid/EA)*cos(phi_mid)-1.);
   // vertical displacement at the end of the step
    u.at(2) = u_prev.at(2) - DX * ( 1. + N_mid / EA ) * sin(phi_mid);
    // bending moment and curvature at the end of the step and their derivatives with respect to the left-end forces
    double M_prev = M; // (store the moment at the beginning of the step, needed for the inflection check)
    M = -Mab+Xab*u.at(2)-Zab*(x+u.at(1));
    for (int j = 1; j <= 3; j++) {
      dM.at(j) = Xab*jacobi.at(2,j)-Zab*jacobi.at(1,j);
    }
    dM.at(1) += u.at(2);
    dM.at(2) += -(x+u.at(1));
    dM.at(3) += -1.;
    kappa = computeCurvatureFromMoment(M);
    dMdkappa = computeDerMomentFromCurvature(kappa);
    dkappa = dM;
    dkappa.times(1./dMdkappa);
    // rotation at the end of the step
    u.at(3) = phi_mid+kappa*DX/2.;
    // test whether Jacobi at the beginning of the step needs to be stored
    if (istep==nstep || (!inflection_detected && M*M_prev<=0. && M_prev!=M)){
      jacobi_prev = jacobi;
    }

    // update Jacobi matrix
    for (int j = 1; j <= 3; j++){
      jacobi.at(1,j) += DX*(dN_mid.at(j)/EA)*cos(phi_mid) - DX*(1.+N_mid/EA)*sin(phi_mid)*dphi_mid.at(j);
      jacobi.at(2,j) -= DX*(dN_mid.at(j)/EA)*sin(phi_mid) + DX*(1.+N_mid/EA)*cos(phi_mid)*dphi_mid.at(j);
      jacobi.at(3,j) = dphi_mid.at(j)+dkappa.at(j)*DX/2.;
    }
     
    // test whether inflection occurs (the first occurence is considered)
    if (!inflection_detected && M*M_prev<=0. && M_prev!=M){
      // inflection point is detected - output parameters Lc, uc and jac_c will be computed
      inflection_detected = true;
      double aux = M_prev / (M_prev-M);
      //@todo:check this
      Lc = x - DX + aux * DX;
      // displacements, rotation and their derivatives at Lc by linear interpolation within the current step
      for (int i = 1; i <= 3; i++){
	uc.at(i) = u_prev.at(i) + aux*(u.at(i)-u_prev.at(i));
	for (int j = 1; j <= 3; j++) {
	  jac_c.at(i,j) = aux*jacobi.at(i,j) + (1.-aux)*jacobi_prev.at(i,j);
	}
      }
      // derivatives with respect to the spatial coordinate
      for (int i = 1; i <= 3; i++) {
	jac_c.at(i,4) = (u.at(i)-u_prev.at(i))/DX;
      }
      // return if the user cares only about the inflection point and not about the segment end
      if (inflection_only)
	return true;
      // change the rotation by deltaPhi (useful for analysis of inverted contact modes)
      u.at(3) += deltaPhi;
    }

    // test whether the end of the segment has been reached
    if (istep==nstep){
      // displacements, rotation and their derivatives at Lb by linear interpolation within the last step
      double aux = Lb/DX - (nstep-1);
      for (int i = 1; i <= 3; i++){
	ub.at(i) = u_prev.at(i) + aux*(u.at(i)-u_prev.at(i));
	for (int j = 1; j <= 3; j++) {
	  jac_b.at(i,j) = aux*jacobi.at(i,j) + (1.-aux)*jacobi_prev.at(i,j);
	}
      }
      // derivatives with respect to the spatial coordinate
      for (int i = 1; i <= 3; i++) {
	jac_b.at(i,4) = (u.at(i)-u_prev.at(i))/DX;
      }
    } 
  } // end of loop over spatial steps
  
  return inflection_detected;
}

void
NlBeamInternalContact :: transform_ub2ua(const FloatArray &ub, FloatArray &answer)
{
  answer.resize(3);
  answer.at(1) = (beamLength+ub.at(1))*cos(ub.at(3)) - ub.at(2)*sin(ub.at(3)) - beamLength;
  answer.at(2) = (beamLength+ub.at(1))*sin(ub.at(3)) + ub.at(2)*cos(ub.at(3));
  answer.at(3) = -ub.at(3);
}

void
NlBeamInternalContact :: transform_fab2fba(const FloatArray &ub, const FloatArray &fab, FloatArray &answer)
{
  answer.resize(3);
  answer.at(1) = fab.at(1)*cos(ub.at(3)) - fab.at(2)*sin(ub.at(3));
  answer.at(2) = fab.at(1)*sin(ub.at(3)) + fab.at(2)*cos(ub.at(3));
  answer.at(3) = fab.at(1)*ub.at(2) - fab.at(2)*(beamLength+ub.at(1)) - fab.at(3);
}


void
NlBeamInternalContact :: transform_fba2fab(const FloatArray &ub, const FloatArray &fba, FloatArray &answer)
{
  FloatArray ua;
  transform_ub2ua(ub, ua);
  transform_fab2fba(ua, fba, answer);
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
NlBeamInternalContact :: checkRotationAdmissibility(double deltaPhi, ContactModeType cmode)
{
  if (cmode == CMT_AC || cmode == CMT_CB) // right segment must rotate counterclockwise wrt to the left one
    return (sin(deltaPhi)>=0.); 
  else // right segment must rotate clockwise
    return (sin(deltaPhi)<=0.); 
}

/////////////////////////////////////////////////////////////////////
// TIP CONTACT
/////////////////////////////////////////////////////////////////////


bool
NlBeamInternalContact :: findLeftEndForcesLocal_Tip_SoS(const FloatArray &ub_target, FloatArray &fab, double &Lac, double Lb, ContactModeType cmode, int process, double &Nca, double &Qca, double &deltaPhi)
{
  bool converged = false;
  FloatArray res(3), fba(3), uac(3),  dforces(3), ubc(3);
  FloatMatrix jacobi_ac(3,4), jacobi_bc(3,4), jacobi(3,3);
  double Mca, error;
  double phib = ub_target.at(3);
  double c = cos(phib);
  double s = sin(phib);
  // transformation matrix 
  FloatMatrix F(3,3);
  F = {{c, s, ub_target.at(2)},{-s,c, -beamLength-ub_target.at(1) },{0,0,-1}};
  // weights for conversion to dimensionless quantities
  double weight_disp = 1./beamLength;
  double weight_mom = beamLength/EI; 
  double weight_force = beamLength*weight_mom; 

  int iter = 0;
  while (iter++ < MAXIT_BEAM){// main iterative loop of the Newton-Raphson method
    
    // integrate along the left segment
    double L_dummy;
    FloatArray u_dummy(3);
    FloatMatrix jac_dummy(3,4);
    if (process==0){ // sticking - we do not care about the inflection point here (dummy variables used)
      trialProcess = PT_Stick;
      this->integrateAlongSegment(fab, 0., Lac, uac, jacobi_ac, L_dummy, u_dummy, jac_dummy, false);
    } else { // sliding
      trialProcess = PT_Slide;
      double La = 10.*leftSegmentLength; // a very safe estimate - we only want to find the inflection point
      bool inflex = integrateAlongSegment(fab, 0., La, u_dummy, jac_dummy, Lac, uac, jacobi_ac, true);
      if (!inflex) {
	return false; // no inflexion point found, contact condition cannot be enforced
      }
    }
    // internal forces at the contact point
    Nca = -fab.at(1) * cos(uac.at(3))+fab.at(2)*sin(uac.at(3));
    Qca = -fab.at(1) * sin(uac.at(3))-fab.at(2)*cos(uac.at(3));
    Mca = -fab.at(3) + fab.at(1) * uac.at(2) - fab.at(2) * (Lac + uac.at(1));
    // from equilibrium, we calculate the corresponding initial guess of the right-end forces
    transform_fab2fba(ub_target, fab, fba);
    // integrate along the right segment - we do not care about the inflection point here
    integrateAlongSegment(fba, 0., Lb, ubc, jacobi_bc, L_dummy, u_dummy, jac_dummy, false);
    // evaluate the residual
    // the first two components are always displacement differences
    res.at(1) = Lac + uac.at(1) - beamLength - ub_target.at(1) + (Lb+ubc.at(1))*c + ubc.at(2)*s;
    res.at(2) = uac.at(2) - ub_target.at(2) - (Lb+ubc.at(1))*s + ubc.at(2)*c;
    // the third component depends on the type of process
    if (process==0){ // sticking
      res.at(3) = Mca;
      //@todo
      FloatArray re(res);
      re.at(1) *= weight_disp;
      re.at(2) *= weight_disp;
      re.at(3) *= weight_mom;
      error = re.computeNorm();
    } else { // sliding
      res.at(3) = evalContactLoadingFunction(Nca, Qca, cmode);
      FloatArray re(res);
      re.at(1) *= weight_disp;
      re.at(2) *= weight_disp;
      re.at(3) *= weight_force;
      error = re.computeNorm();
    }
    // check the convergence criterion
    if (error<TOL_BEAM){// converged solution found
      // evaluate the rotation jump at the contact point
      deltaPhi = ub_target.at(3) + ubc.at(3) - uac.at(3);
      if (!stiffEvalMode){
	return true;
      }
      converged = true; // return postponed, since Jacobi should be updated
    }
    
    // set up the Jacobi matrix (derivatives of the residual wrt to the left-end forces)

    // prepare auxiliary variables
    FloatArray du(3), dw(3), dN(3), dQ(3);
    for (int j = 1; j <= 3; j++){
      for (int i = 1; i <= 3; i++){
	du.at(j) += jacobi_bc.at(1,i)*F.at(i,j);
	dw.at(j) += jacobi_bc.at(2,i)*F.at(i,j);
      }
    }
    if (process==1){ // sliding  
      // prepare derivatives of internal forces at the contact point
      for ( int j = 1; j <= 3; j++){
	//@todo:check
	dN.at(j) = (fab.at(1)*sin(uac.at(3))+fab.at(2)*cos(uac.at(3))) * jacobi_ac.at(3,j);
 	dQ.at(j) = (-fab.at(1)*cos(uac.at(3))+fab.at(2)*sin(uac.at(3))) * jacobi_ac.at(3,j);
      }
      dN.at(1) += -cos(uac.at(3));
      dN.at(2) +=  sin(uac.at(3));
      dQ.at(1) += -sin(uac.at(3));
      dQ.at(2) += -cos(uac.at(3));
    }
    // evaluate the part of Jacobi matrix that has a regular structure (additional terms come later)
    for (int j = 1; j <= 3; j++){
      jacobi.at(1,j) = jacobi_ac.at(1,j) + du.at(j)*c + dw.at(j)*s;
      jacobi.at(2,j) = jacobi_ac.at(2,j) - du.at(j)*s + dw.at(j)*c;
      if (process==0) // sticking
	jacobi.at(3,j) = fab.at(1)*jacobi_ac.at(2,j)-fab.at(2)*jacobi_ac.at(1,j); 
      else // sliding
	jacobi.at(3,j) = evalDerContactLoadingFunction(Nca, dN.at(j), dQ.at(j), cmode);
    }
    if (process==0){ // additional terms for sticking (correction of the third row)
      jacobi.at(3,1) += uac.at(2);
      jacobi.at(3,2) += -(Lac+uac.at(1));
      jacobi.at(3,3) += -1.;
    } else { // additional terms for sliding (correction of the first two rows)
      // this piece of code takes into account the fact that changes of end forces lead to changes in the position of inflection point
      double dMdx = fab.at(1)*jacobi_ac.at(2,4)-fab.at(2)*(1.+jacobi_ac.at(1,4)); 
      FloatArray dMc(3), dLac(3);
      dMc.at(1) = uac.at(2);
      dMc.at(2) = -(Lac+uac.at(1));
      dMc.at(3) = -1.;
      for (int j = 1; j <= 3; j++){
	dMc.at(j) += fab.at(1)*jacobi_ac.at(2,j)-fab.at(2)*jacobi_ac.at(1,j);
	dLac.at(j) = -dMc.at(j)/dMdx; // derivatives of the position of inflection point wrt the left-end forces
      }
      for (int j = 1; j <= 3; j++){
	jacobi.at(1,j) += (1.+jacobi_ac.at(1,4))*dLac.at(j);
	jacobi.at(2,j) += jacobi_ac.at(2,4)*dLac.at(j);
      }  
    }

    if (converged && stiffEvalMode){ // prepare a block of the stiffness matrix
      //@todo: check this
      Kblock.beInverseOf(jacobi);
      double c1 = (Lb+ubc.at(1))*s - ubc.at(2)*c; 
      double c2 = (Lb+ubc.at(1))*c + ubc.at(2)*s;
      // the last column is a linear combination of the first two columns
      for (int i = 1; i <= 3; i++)
	Kblock.at(i,3) = c1*Kblock.at(i,1) + c2*Kblock.at(i,2);
      return true;
    }
    
    // solve the linearized problem
    jacobi.solveForRhs(res, dforces);
    // correct the left-end forces
    fab.subtract(dforces); 
  } // end of while loop
  return false;
}


double
NlBeamInternalContact :: shiftToIntervalFromMinusPiToPi(double x)
{
  double answer = x - 2.*M_PI*floor((x+M_PI)/(2.* M_PI));
  return answer;
}


ContactModeType
NlBeamInternalContact :: suggestedMode(double deltaPhi, bool seglength_ok, ContactModeType current_cmode)
{
  // if the tip in contact has slipped beyond the other tip, it is possible that the other tip gets activated
  if (!seglength_ok){
    switch (current_cmode){
    case CMT_AC:
      return CMT_CB;
    case CMT_BC:
      return CMT_CA;
    case CMT_CA:
      return CMT_BC;
    case CMT_CB:
      return CMT_AC;
    }
  }
  // now we know that the tip in contact has NOT slipped beyond the other tip
  // the rotation jump is first shifted into the interval from -PI to PI 
  deltaPhi = shiftToIntervalFromMinusPiToPi(deltaPhi);
  // if the rotation jump is admissible, the current_cmode should be the correct one,
  // otherwise we can expect a change into a smooth mode, depending on which boundary
  // of the admissible interval has been crossed
  switch (current_cmode){
  case CMT_AC:
    if (deltaPhi>=0.)
      return CMT_AC;
    else if (deltaPhi>-M_PI/2.)
      return CMT_AA;
    else
      return CMT_AB;
  case CMT_BC:
    if (deltaPhi<=0.)
      return CMT_BC;
    else if (deltaPhi<M_PI/2.)
      return CMT_BB;
    else
      return CMT_BA;
  case CMT_CA: 
    if (deltaPhi<=0.)
      return CMT_CA;
    else if (deltaPhi<M_PI/2.)
      return CMT_AA;
    else
      return CMT_BA;
  case CMT_CB: 
    if (deltaPhi>=0.)
      return CMT_CB;
    if (deltaPhi>-M_PI/2.)
      return CMT_BB;
    else
      return CMT_AB;
  }
  return current_cmode;
}


ContactModeType
NlBeamInternalContact :: findLeftEndForcesLocal_Tip(FloatArray &ub_target, FloatArray &fab, ContactModeType cmode)
{
  double Nca, Qca, deltaPhi, fc, Lac;
  bool  success, rotation_ok;
  // store the initial guess for possible restart with sliding assumption
  FloatArray f_init(3), fba(3), ua_target(3);
  
  bool tipIsOnRightSegment;
  switch(cmode){
  case CMT_AC: case CMT_BC:
    tipIsOnRightSegment = true;
    break;
  case CMT_CA: case CMT_CB:
    tipIsOnRightSegment = false;
    transform_fab2fba(ub_target, fab, fba);
    break;
  default:
    printf("findLeftEndForcesLocal_Tip cannot be used for cmode = %d\n",cmode);
    exit(1);
  }

  if (Process!=PT_Slide){// if it is known that the process is sliding, proceed directly to that case
  // sticking is assumed first
    
    trialProcess = PT_Stick;
    if (tipIsOnRightSegment){
      Lac = trialLeftActiveSegmentLength;
      f_init = fab; // fab stored for possible restart with sliding
      success = findLeftEndForcesLocal_Tip_SoS(ub_target, fab, Lac, rightSegmentLength, cmode, 0, Nca, Qca, deltaPhi);  
    } else {
      Lac = trialRightActiveSegmentLength;
      // from equilibrium, we calculate the corresponding initial guess of the right-end forces
      f_init = fba; // fba stored for possible restart with sliding
      transform_ub2ua(ub_target, ua_target);
      success = findLeftEndForcesLocal_Tip_SoS(ua_target, fba, Lac, leftSegmentLength, cmode, 0, Nca, Qca, deltaPhi);
      deltaPhi = -deltaPhi; // compensating for swapped segments
      transform_fba2fab(ub_target, fba, fab);
    }
    if (!success){
      // give sliding a chance - make sure that restart takes place with the original guess of fab or fba
      if (tipIsOnRightSegment) {
	fab = f_init;
      } else {
	fba = f_init;
      }
    } else { // a solution for sticking has been found but its admissibility still needs to be checked
      fc = evalContactLoadingFunction(Nca, Qca, cmode);
      rotation_ok = checkRotationAdmissibility(deltaPhi, cmode);
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
  trialProcess = PT_Slide;
  if (tipIsOnRightSegment) {
    Lac = trialLeftActiveSegmentLength;
    success = findLeftEndForcesLocal_Tip_SoS(ub_target, fab, Lac, rightSegmentLength, cmode, 1, Nca, Qca, deltaPhi);
  }  else {
    Lac = trialRightActiveSegmentLength;
    success = findLeftEndForcesLocal_Tip_SoS(ua_target, fba, Lac, leftSegmentLength, cmode, 1, Nca, Qca, deltaPhi);
    deltaPhi = -deltaPhi; // compensating for swapped segments
  }
  if (!success){
    trialLeftActiveSegmentLength = leftSegmentLength;
    trialRightActiveSegmentLength = rightSegmentLength;
    return CMT_N;
  }
  // a solution for sliding has been found but its admissibility still needs to be checked
  rotation_ok = checkRotationAdmissibility(deltaPhi, cmode);
  bool seglength_ok;
  if (tipIsOnRightSegment) {
    seglength_ok = checkLeftSegmentLengthAdmissibility(Lac);
  } else {
    seglength_ok = checkRightSegmentLengthAdmissibility(Lac);
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


bool
NlBeamInternalContact :: findLeftEndForcesLocal_Smooth_Rolling(FloatArray &ub_target, FloatArray &fab, double deltaPhi, double &Lac, double &Lbc, double &Nc, double &Qc)
{

  FloatArray res(3), dforces(3), ub_loc(3), uc(3);
  FloatMatrix jac_c_dummy(3,4), jac(3,3), jac_b(3,4);
 
  double L = Lac + Lbc;
  // weight for conversion to dimensionless quantities
  double weight_disp = 1./beamLength;

  int iter = 0;
  while (iter++<MAXIT_BEAM){
    // integrate along beam of length L
     bool inflection_found = integrateAlongSegment(fab, deltaPhi, L, ub_loc, jac_b, Lac, uc, jac_c_dummy, false);
     Jacobi.beSubMatrixOf(jac_b,1,3,1,3);
     
     Nc = -fab.at(1)*cos(uc.at(3))+fab.at(2)*sin(uc.at(3));
     Qc = -fab.at(1)*sin(uc.at(3))-fab.at(2)*cos(uc.at(3));

     // evaluate the residual and check the convergence criterion
    for (int i = 1; i <= 3; i++)
      res.at(i) = ub_target.at(i)-ub_loc.at(i); 
    // take into account the fact that, due to previous sliding,
    // the current beam length may differ from the initial distance between joints
    // (displacement ub_target[0] is taken with respect to the initial length of the beam)
    res.at(1) += beamLength - L;
    double error = L2norm(weight_disp*res.at(1), weight_disp*res.at(2), res.at(3));
    if (error<TOL_BEAM){// converged
      if (!inflection_found)
	return false;// converged but no inflexion detected
      else{
	Lbc = L - Lac;
	return true;// converged and inflexion detected
      }
    }
    // compute the iterative correction of left-end forces
    jac = Jacobi;
    jac.solveForRhs(res,dforces);
    fab.add(dforces);
  }
  // no convergence in MAXIT_BEAM iterations
    return false;
}



bool
NlBeamInternalContact :: findLeftEndForcesLocal_Smooth_Sliding(FloatArray &ub_target, FloatArray &fab, double deltaPhi, double &Lac, double &Lbc, ContactModeType cmode)
{
  bool converged = false;
  FloatArray res(4), dforces(4), ub_loc(3), uc(3);
  FloatMatrix jac_b(3,4), jac_c(3,4), jac(4,4);
  double Nc, Qc;
  FloatArray dN(3), dQ(3);
  double L = Lac + Lbc;
  // weights for conversion to dimensionless quantities
  double weight_disp = 1./beamLength;
  double weight_force = beamLength*beamLength/EI; 

  int iter = 0;
  while (iter++<MAXIT_BEAM){
    bool inflection_found = integrateAlongSegment(fab, deltaPhi, L, ub_loc, jac_b, Lac, uc, jac_c,false); 
    if (!inflection_found) {
      return false; // no inflexion point - contact condition cannot be enforced
    }
    // internal forces at the contact point
    Nc = -fab.at(1)*cos(uc.at(3))+fab.at(2)*sin(uc.at(3));
    Qc = -fab.at(1)*sin(uc.at(3))-fab.at(2)*cos(uc.at(3));
    // evaluate the residual and check the convergence criterion
    res = ub_loc;
    res.subtract(ub_target);
    res.resizeWithValues(4);
    res.at(4) = evalContactLoadingFunction(Nc, Qc, cmode);
    // take into account the fact that, due to sliding,
    // the current beam length may differ from the initial distance between joints
    // (displacement ub_target[0] is taken with respect to the initial length of the beam)
    res.at(1) += L - beamLength;
    double error = L2norm(weight_disp*sqrt(res.at(1)*res.at(1)+res.at(2)*res.at(2)), res.at(3), weight_force*res.at(4));
   
    if (error < TOL_BEAM) {
      // the actual return is a bit postponed, to make update of Jacobi possible if needed
      converged = true; 
    }

    // prepare the Jacobi matrix for the set of 4 equations
    if (!converged || stiffEvalMode){
      // evaluate the first three rows
      for(int i = 1; i <= 3; i++) {
	for (int j = 1; j <= 4; j++) {
	  Jacobi44.at(i,j) = jac_b.at(i,j);
	}
      }
      // correct one entry by a term that reflects the line res[0] += L - beamLength;
      Jacobi44.at(1,4) += 1.;
      // derivatives of internal forces at the contact point
      for (int j = 1; j <= 3; j++){
	dN.at(j) = (fab.at(1)*sin(uc.at(3))+fab.at(2)*cos(uc.at(3))) * jac_c.at(3,j);
	dQ.at(j) = (-fab.at(1)*cos(uc.at(3))+fab.at(2)*sin(uc.at(3))) * jac_c.at(3,j);
      }
      dN.at(1) += -cos(uc.at(3));
      dN.at(2) += sin(uc.at(3));
      dQ.at(1) += -sin(uc.at(3));
      dQ.at(2) += -cos(uc.at(3));
      // evaluate the fourth row of Jacobi matrix 
      for (int j = 1; j <= 3; j++) {
	Jacobi44.at(4,j) = evalDerContactLoadingFunction(Nc, dN.at(j), dQ.at(j), cmode);
      }
      //Jacobi44.at(4,4) = 0.; // the total length does not affect the inflection point
    }
    // this is the postponed return after convergence
    if (converged){
      Lbc = L - Lac;
      return true;
    }
   // compute the iterative correction of left-end forces and of active beam length
    jac = Jacobi44;
    jac.solveForRhs(res,dforces);
    for(int i = 1; i <= 3; i++) {
      fab.at(i) -= dforces.at(i);
    }
    
    L -= dforces.at(4);
  }
  // no convergence in MAXIT_BEAM iterations
  return false;
}


ContactModeType
NlBeamInternalContact :: findLeftEndForcesLocal_Smooth(FloatArray &ub_target, FloatArray &fab, ContactModeType cmode)
{
  double Nc, Qc;
  // store the initial guess for possible restart with sliding assumption
  FloatArray fab_init(3);
  double deltaPhi, Lac, Lbc;
  bool success;
  fab_init = fab;

   deltaPhi = 0.;
   if (cmode==CMT_AB) {
     deltaPhi = M_PI;
   } else if (cmode==CMT_BA) {
     deltaPhi = -M_PI;
   }
   
  if (Process!=PT_Slide){// if it is known that the process is sliding, proceed directly to that case
    // first assume that no additional sliding occurs during this step (i.e., fix Lac+Lbc and iterate on fab)
    trialProcess = PT_Roll;  
    Lac = trialLeftActiveSegmentLength;
    Lbc = trialRightActiveSegmentLength;
    bool success = findLeftEndForcesLocal_Smooth_Rolling(ub_target, fab, deltaPhi, Lac, Lbc, Nc, Qc);
    // the result is 'false' if the left-end forces do not converge, or if they converge but the beam has no inflection point 
    if (!success) {
      // give sliding a chance - make sure that the solution that assumes sliding will later start from the original guess of left-end forces 
      fab = fab_init;
    } else { // sticking solution has been found, now let us see whether it is admissible
      double fc = evalContactLoadingFunction(Nc, Qc, cmode);
      if (fc<=0.){ // admissible solution
	trialLeftActiveSegmentLength = Lac;
	trialRightActiveSegmentLength = Lbc;
	return cmode;
      }
    }
  }
  // now assume that sliding occurs during this step (i.e., iterate on fab as well as on Lac+Lbc)
  trialProcess = PT_Slide;  
  Lac = trialLeftActiveSegmentLength;
  Lbc = trialRightActiveSegmentLength;
  success = findLeftEndForcesLocal_Smooth_Sliding(ub_target, fab, deltaPhi, Lac, Lbc, cmode);
  if (success){ // the active segment lengths change due to sliding
    if ((Lac+Lbc)>leftSegmentLength+rightSegmentLength)
      return CMT_N; // most likely jumps to no contact
    else if (Lac > leftSegmentLength){
      trialLeftActiveSegmentLength = leftSegmentLength;
      trialRightActiveSegmentLength = Lbc;
      if (cmode == CMT_AA)
	return CMT_CA; // most likely changes into tip contact
      else
	return CMT_CB;
    }
    else if (Lbc > rightSegmentLength) {
      trialLeftActiveSegmentLength = Lac;
      trialRightActiveSegmentLength = rightSegmentLength;
      if (cmode == CMT_AA)
	return CMT_AC; // most likely changes into tip contact 
      else
	return CMT_BC;
    }
    // admissible solution
    trialLeftActiveSegmentLength = Lac;
    trialRightActiveSegmentLength = Lbc;
    return cmode;
  }
  else
    return CMT_N; // most likely smoothly changes into no contact
}


ContactModeType
NlBeamInternalContact :: predictContactMode(FloatArray &ub)
{
  // assuming that both segments are straight, find their intersection
  
  double phib = ub.at(3);
  if (sin(phib)==0.) // if both segments are parallel, the no-contact mode is always possible
    return CMT_N;
  double Lbc = -ub.at(2)/sin(phib);
  if (Lbc<0. || Lbc>rightSegmentLength) {
    return CMT_N;
  }
  double Lac = beamLength+ub.at(1)-Lbc*cos(phib);
  if (Lac < 0. || Lac > leftSegmentLength) {
    return CMT_N;
  }
  if (rightSegmentLength - Lbc < leftSegmentLength - Lac) { // right tip probably touches the left segment first 
    if (ub.at(2)>0.) {
      return CMT_BC;
    } else {
      return CMT_AC;
    }
  }
  // left tip probably touches the right segment first
  if (ub.at(2) > 0.) {
    return CMT_CA;
  }
  return CMT_CB;
}


double
NlBeamInternalContact :: findTipContactTime(FloatArray &u_prev, FloatArray &du, double segLength)
{
  double t;
  // simple linearized estimate, used as the initial guess for the iterative scheme
  double db = du.at(2) + segLength*(sin(u_prev.at(3)+du.at(3))-sin(u_prev.at(3))); 
  if (db!=0.){
    t = (-u_prev.at(2)-segLength*sin(u_prev.at(3)))/db;
  } else {
    t = 0.5;
  }
  // iterative solution of nonlinear equation
  for (int iter = 1; iter <= MAXIT_CONTACT_TIME; iter++){ 
    double f = u_prev.at(2) + t*du.at(2) + segLength*sin(u_prev.at(3)+t*du.at(3));
    if (fabs(f) < TOL_CONTACT_TIME * segLength) {
	return t;// converged solution
    }
    double fp = du.at(2) + segLength*cos(u_prev.at(3)+t*du.at(3))*du.at(3);
    if (fp!=0.) {
	t -= f/fp;
    }
  }
  // converged solution not found
  return -1;
}


ContactModeType
NlBeamInternalContact :: predictContactMode(FloatArray ub, FloatArray ub_prev)
{
  // assuming that both segments are straight,
  //        check at which time 'tb' and location 'xb' the right tip hits the left segment
  // (the 'time' is a dimensionless parameter running from 0 at the beginning to 1 at the end of the step)
  double xa, xb;
  FloatArray du(3);
  du = ub;
  du.subtract(ub_prev);
  double tb = findTipContactTime(ub_prev, du, rightSegmentLength);
  if (tb>=0. && tb<=1.) {
    // find location of the hit event
    xb = beamLength + ub_prev.at(1) + tb*du.at(1) - rightSegmentLength*cos(ub_prev.at(3)+tb*du.at(3));
  }
  // now check at which time 'ta' and location 'xa' the left tip hits the right segment
  FloatArray ua(3), ua_prev(3);
  transform_ub2ua(ub, ua);
  transform_ub2ua(ub_prev, ua_prev);
  du = ua;
  du.subtract(ua_prev);
  double ta = findTipContactTime(ua_prev, du, leftSegmentLength);
  if (ta >= 0. && ta <= 1.) {
    // find location of the hit event
    xa = beamLength + ua_prev.at(1) + ta*du.at(1) - leftSegmentLength*cos(ua_prev.at(3)+ta*du.at(3));
  }
  // determine which hit event occurs at an admissible time and location
  bool left_hits = ( ta >= 0. && ta <= 1. && xa >= 0. && xa <= rightSegmentLength);
  bool right_hits = (tb >= 0. && tb <= 1. && xb >= 0. && xb <= leftSegmentLength);

  if (!left_hits && !right_hits) {
    return CMT_N;
  }
  if (left_hits && right_hits) {
    left_hits = ta < tb;
  }
  if (left_hits) {
    trialLeftActiveSegmentLength = leftSegmentLength;
    trialRightActiveSegmentLength = xa;
    if (ua.at(2)-ua_prev.at(2)+leftSegmentLength*(sin(ua.at(3))-sin(ua_prev.at(3))) > 0.) {
      return CMT_CA;
    } else {
      return CMT_CB;
    }
  }
  trialLeftActiveSegmentLength = xb;
  trialRightActiveSegmentLength = rightSegmentLength;
  if (ub.at(2)-ub_prev.at(2)+rightSegmentLength*(sin(ub.at(3))-sin(ub_prev.at(3))) > 0.) {
    return CMT_AC;
  } else {
    return CMT_BC;
  }
}


bool
NlBeamInternalContact :: findLeftEndForcesLocal(FloatArray &ub, FloatArray &ub_prev, FloatArray &fab)
{
  ContactModeType assumed_cmode = contactMode;
  trialLeftActiveSegmentLength = leftActiveSegmentLength;
  trialRightActiveSegmentLength = rightActiveSegmentLength;

  for (int attempt = 1; attempt <= MAXIT_CONTACT_MODE; attempt++) {
    switch (assumed_cmode){
      // smooth modes
    case CMT_AA:
    case CMT_BB:
    case CMT_AB:
    case CMT_BA:
      trialContactMode = findLeftEndForcesLocal_Smooth(ub, fab, assumed_cmode);
      break;
      // tip modes
    case CMT_AC:
    case CMT_BC:
    case CMT_CA:
    case CMT_CB:
      trialContactMode = findLeftEndForcesLocal_Tip(ub, fab, assumed_cmode);
      break;
      // no-contact mode
    case CMT_N:
      fab.zero();
      if (contactMode != CMT_N) {
	trialContactMode = predictContactMode(ub);
      } else {
	trialContactMode = predictContactMode(ub, ub_prev);
      }
      break;
      // this should never happen
    default:
      OOFEM_ERROR("Mode %d has not been implemented yet in findLeftEndForcesLocal\n",assumed_cmode);
    }
    if ( trialContactMode == assumed_cmode) {// solution accepted, update internal variables
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

void
NlBeamInternalContact :: construct_T(FloatMatrix &T, double phia)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = cosAlpha*cos(phia) + sinAlpha*sin(phia);
  T.at(1,2) =  sinAlpha *cos(phia) - cosAlpha * sin(phia);
  T.at(2,1) = - T.at(1,2);
  T.at(3,3) = 1.;
  /*
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = cos(alpha-phia);
  T.at(1,2) = sin(alpha-phia);
  T.at(2,1) = -T.at(1,2);
  T.at(3,3) = 1.;
  */
}

void
NlBeamInternalContact :: construct_T(FloatMatrix &T, double phia, double rot)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = (cosAlpha*cos(rot) + sinAlpha*sin(rot))*cos(phia) + (sinAlpha*cos(rot)-cosAlpha*sin(rot))*sin(phia);
  T.at(1,2) =  (sinAlpha*cos(rot)-cosAlpha*sin(rot)) *cos(phia) - (cosAlpha*cos(rot)+sinAlpha*sin(rot)) * sin(phia);
  T.at(2,1) = - T.at(1,2);
  T.at(3,3) = 1.;
}


void
NlBeamInternalContact :: construct_Tprime(FloatMatrix &T, double phia)
{
    T.resize(3,3);
  T.at(1,1) = T.at(2,2) = sinAlpha *cos(phia) - cosAlpha*sin(phia);
  T.at(1,2) = -cosAlpha*cos(phia) - sinAlpha*sin(phia);
  T.at(2,1) = -T.at(1,2);
  /*
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = sin(alpha-phia);
  T.at(1,2) = -cos(alpha-phia);
  T.at(2,1) = -T.at(1,2);
  */

}

void
NlBeamInternalContact :: construct_l(FloatArray &l, double phia)
{
  l.resize(3);
  l.at(1) = beamLength * (cos(phia)-1.);
  l.at(2) = beamLength * sin(phia);
}

void
NlBeamInternalContact :: construct_l(FloatArray &l, double phia, double L)
{
  l.resize(3);
  l.at(1) = L * (cos(phia) - 1.);
  l.at(2) = L *  sin(phia);
}



void
NlBeamInternalContact :: construct_l_IC(FloatArray &l, double phia, double L)
{
  l.at(1) = beamLength * cos(phia) - L;
  l.at(2) = beamLength * sin(phia);
}

void
NlBeamInternalContact :: construct_lprime(FloatArray &l, double phia)
{
  l.resize(3);
  l.at(1) = -beamLength * sin(phia);
  l.at(2) = beamLength * cos(phia);
}


bool
NlBeamInternalContact :: findLeftEndForces(const FloatArray &u, const FloatArray &u_prev, FloatArray &fab)
{
  FloatArray ub_loc(3), ub_prev_loc(3), fab_loc(3);
  FloatMatrix T(3,3), T_prev(3,3);
  double phia = u.at(3);
  double phia_prev = u_prev.at(3);
  Process = trialProcess = PT_Unknown;
  // compute displacements of the right end with respect to the auxiliary coordinate system
  FloatArray l_ub_loc, l_ub_prev_loc;
  construct_l(l_ub_loc, phia);
  construct_l(l_ub_prev_loc, phia_prev);
  construct_T(T, phia); 
  construct_T(T_prev, phia_prev);


  FloatArray u_a,u_b, u_ba, uprev_ba, uprev_a, uprev_b;
  u_a.beSubArrayOf(u, {1,2,3});
  u_b.beSubArrayOf(u, {4,5,6});
  u_ba.beDifferenceOf(u_b, u_a);
  
  uprev_a.beSubArrayOf(u_prev, {1,2,3});
  uprev_b.beSubArrayOf(u_prev, {4,5,6});
  uprev_ba.beDifferenceOf(uprev_b, uprev_a);

  ub_loc.beProductOf(T, u_ba);
  ub_loc.add(l_ub_loc);
  ub_prev_loc.beProductOf(T_prev,uprev_ba);
  ub_prev_loc.add(l_ub_prev_loc);
  fab_loc.beProductOf(T,fab);
  // compute the corresponding left-end forces and moment in the auxiliary coordinate system
  if(!findLeftEndForcesLocal(ub_loc, ub_prev_loc, fab_loc)) {
    return false;
  }
  //@todo: comment the following 3 lines???
  /* contactMode = trialContactMode;
  leftActiveSegmentLength = trialLeftActiveSegmentLength;
  rightActiveSegmentLength = trialRightActiveSegmentLength;
  */
  
  fab.beTProductOf(T,fab_loc);
  return true;
}


void
NlBeamInternalContact :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  /* When evaluating the end forces "from scratch", we set the process status to unknown.
     The algorithm always tries sticking first and if the solution is not admissible or not found,
     it proceeds to the sliding assumption. This is always the case when we call findEndForces.
     However, it is also possible to call directly findLeftEndForces and then, if the process
     status has been previously detected and stored, the algorithm skips the sticking/rolling part
     in cases when it is know that the actual process is sliding. This is used when the stiffness
     is evaluated after previous evaluation of internal forces for the same iterated increment.
  */
  answer.resize(6);
  FloatArray u;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
  if(!useUpdatedGpRecord) {
    // solution vector
    FloatArray ui, u_prev, f(3);
    this->computeVectorOf({D_u, D_w, R_v}, VM_Incremental, tStep, ui);
    u_prev = u - ui;
    Process = trialProcess = PT_Unknown;
    this->findLeftEndForces(u, u_prev, this->internalForces); // only the first three entries of f are computed
    //
  }
    answer.at(1) = this->internalForces.at(1);
    answer.at(2) = this->internalForces.at(2);
    answer.at(3) = this->internalForces.at(3);  
    answer.at(4) = -answer.at(1);
    answer.at(5) = -answer.at(2);
    double c1 = beamLength*sin(pitch) + u.at(5) - u.at(2);
    double c2 = -beamLength*cos(pitch) - u.at(4) + u.at(1);
    answer.at(6) = c1 * answer.at(1) + c2 * answer.at(2) - answer.at(3);
  
}


void
NlBeamInternalContact :: giveInternalForcesVector_fromU(FloatArray &answer, TimeStep *tStep, FloatArray &u, FloatArray &u_prev)
{
  /* When evaluating the end forces "from scratch", we set the process status to unknown.
     The algorithm always tries sticking first and if the solution is not admissible or not found,
     it proceeds to the sliding assumption. This is always the case when we call findEndForces.
     However, it is also possible to call directly findLeftEndForces and then, if the process
     status has been previously detected and stored, the algorithm skips the sticking/rolling part
     in cases when it is know that the actual process is sliding. This is used when the stiffness
     is evaluated after previous evaluation of internal forces for the same iterated increment.
  */
  answer.resize(6);
  FloatArray f(3);
  this->findLeftEndForces(u, u_prev, f); // only the first three entries of f are computed
    answer.at(1) = f.at(1);
    answer.at(2) = f.at(2);
    answer.at(3) = f.at(3);  
    answer.at(4) = -answer.at(1);
    answer.at(5) = -answer.at(2);
    double c1 = beamLength*sin(pitch) + u.at(5) - u.at(2);
    double c2 = -beamLength*cos(pitch) - u.at(4) + u.at(1);
    answer.at(6) = c1 * answer.at(1) + c2 * answer.at(2) - answer.at(3);
  
}


void
NlBeamInternalContact :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
  answer.resize(6,6);
  FloatArray u(6), ui(6), u_prev(6);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Incremental, tStep, ui);
  u_prev  = u - ui;
  double L, L_prev;
  FloatArray lprime(3), fab_loc(3), ub_loc(3), Tu(3), fab(3);
  FloatMatrix T(3,3), Tprime(3,3), Ginv(3,3), Ginv44(4,4), TtGinv(3,3);
  //  FloatArray f(6);
  // @todo: check this
  fab = this->internalForces;
  //  bool s = findLeftEndForces(u, u_prev, fab);
  // compute auxiliary matrices
  double phia = u.at(3);
  construct_T(T, phia); 
  construct_Tprime(Tprime, phia); 
  construct_lprime(lprime, phia);
  //  
  fab_loc.beProductOf(T,fab);
  // get Jacobi matrix in local coordinates and invert it
  switch (trialContactMode) {
    // smooth modes
  case CMT_AA:
  case CMT_BB:
  case CMT_AB:
  case CMT_BA:
    L_prev = leftActiveSegmentLength + rightActiveSegmentLength;
    L = trialLeftActiveSegmentLength + trialRightActiveSegmentLength;
    if (fabs(L - L_prev) > TOL_LENGTH * beamLength) { // sliding
      Ginv44.beInverseOf(Jacobi44);
      Ginv.beSubMatrixOf(Ginv44,1,3,1,3);
    } else { // rolling
      Ginv.beInverseOf(Jacobi);
    }
    break;
    // tip modes
  case CMT_AC:
  case CMT_BC:
    Ginv.beSubMatrixOf(Kblock,1,3,1,3);
    break;
    // no-contact mode
  case CMT_N:
    answer.zero();
    return;
    // this should never happen
  default:
    Ginv.beSubMatrixOf(Kblock,1,3,1,3);
    //comment OOFEM_ERROR("Mode %d has not been implemented yet in evalStiffnessMatrix\n",trialContactMode);
  }

  TtGinv.beTProductOf(T, Ginv);

  
  // compute product Ttransposed*Ginverse*T and store it in the upper stiffness block
  FloatMatrix aux;
  aux.beProductOf(TtGinv, T);
  answer.setSubMatrix(aux, 1,1);
  answer.times(-1);
  answer.setSubMatrix(aux, 1,4);

  FloatArray u_a,u_b, u_ba;
  u_a.beSubArrayOf(u, {1,2,3});
  u_b.beSubArrayOf(u, {4,5,6});
  u_ba.beDifferenceOf(u_b, u_a);

  // compute additional stiffness terms in the third column (associated with phia)
  Tu.beProductOf(Tprime, u_ba);
  FloatArray col3_1,col3_2;
  Tu.add(lprime);
  col3_1.beProductOf(TtGinv,Tu);
  col3_2.beTProductOf(Tprime, fab_loc);
  answer.addSubVectorCol(col3_1, 1,3);
  answer.addSubVectorCol(col3_2, 1,3);
			      
  // construct the sixth row (complete formula)
  double c1 = beamLength*sinAlpha + u.at(5) - u.at(2);
  double c2 = -beamLength*cosAlpha - u.at(4) + u.at(1);
  answer.at(6,1) =  fab.at(2);
  answer.at(6,2) = -fab.at(1);
  answer.at(6,4) = -fab.at(2);
  answer.at(6,5) =  fab.at(1);
  for(int j = 1; j <= 6; j++) {
    answer.at(4,j) = - answer.at(1,j);
    answer.at(5,j) = - answer.at(2,j);
    answer.at(6,j) += c1 * answer.at(1,j) + c2 * answer.at(2,j) - answer.at(3,j);
  }

  //compute the stiffness matrix numerically
  /*
  double du = 1.e-6; // small perturbation for numerical evaluation of stiffness
  FloatArray f;
  f = this->internalForces;
  FloatMatrix f_pert(6,6), Knum(6,6);

  for (int i=1; i<=6; i++){
    u.at(i) += du;
    this->giveInternalForcesVector_fromU(f, tStep, u, u_prev);
    u.at(i) -= du;
    for (int j=1; j<=6; j++){
      f_pert.at(j,i) = f.at(j);
    }
  }
  this->giveInternalForcesVector_fromU(f, tStep, u, u_prev);

  double stifferr = 0;
  double stiffdiagsum = 0;
  for (int i=1; i<=6; i++){
    for (int j=1; j<=6; j++){
      Knum.at(i,j) = (f_pert.at(i,j)-f.at(i))/du;
      stifferr += fabs(Knum.at(i,j)-answer.at(i,j));
      if (i==j) {
	stiffdiagsum += answer.at(i,j);
      }
    }
  }
  double ei = stifferr/stiffdiagsum;
  */

}

void
NlBeamInternalContact :: updateYourself(TimeStep *tStep)
{
    NlBeam_SM :: updateYourself(tStep);
    contactMode = trialContactMode;
    leftActiveSegmentLength = trialLeftActiveSegmentLength;
    rightActiveSegmentLength = trialRightActiveSegmentLength;
}


Interface *NlBeamInternalContact :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case VTKXMLExportModuleElementInterfaceType:
      return static_cast< VTKXMLExportModuleElementInterface * >( this );
    default:
      return StructuralElement :: giveInterface(it);
    }
}

void
NlBeamInternalContact :: printOutputAt(FILE *file, TimeStep *tStep)
{
  NlBeam_SM :: printOutputAt(file, tStep);
  fprintf(file, " Contact Mode %d  ", this->contactMode);
}


void
NlBeamInternalContact :: giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
  //two segments
  vtkPieces.resize(2);
  
  FloatArray ul, xLeft, xRight;
  FloatMatrix uMatrixLeft, uMatrixRight, T;
  FloatArray uab, ua(3),ub(3);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, uab);
  ua.beSubArrayOf(uab,{1,2,3});
  ub.beSubArrayOf(uab,{4,5,6});
  //
  /*To be removed*/
  double phia = uab.at(3);
  FloatArray ub_loc, ub_ua(ub), loc, fab_loc, u0(2);
  ub_ua.subtract(ua);  
  // compute displacements of the right end with respect to the auxiliary coordinate system
  construct_l(loc, phia);
  construct_T(T, phia); 
  ub_loc.beProductOf(T, ub_ua);
  ub_loc.add(loc);
  fab_loc.beProductOf(T, this->internalForces);
  FloatMatrix Tn = {{1,0},{0,1}};  
  this->computeSegmentDisplacements(uMatrixLeft,fab_loc, leftActiveSegmentLength, leftSegmentLength, u0, Tn, uab );
  // 
  ////
  FloatMatrix Tb,Ttb;
  FloatArray ur(2);
  Node *nodeB = this->giveNode(2);
  ur.at(1) = beamLength+ub_loc.at(1);
  ur.at(2) = ub_loc.at(2);
  double alphab = M_PI + ub_loc.at(3);
  construct_T(Tb, alphab);
  Ttb.beTranspositionOf(Tb);
  // transformation of forces and moment
  FloatArray f;
  transform_fab2fba(ub_loc, fab_loc, f);
  this->computeSegmentDisplacements(uMatrixRight,f, rightActiveSegmentLength, rightSegmentLength, ur, Ttb, uab );



  
  /*
  ///
  u.at(1) = ua.at(1);
  u.at(2) = ua.at(2);
  //
  Node *nodeA = this->giveNode(1);
  ul.resize(2);
  FloatMatrix Tt = {{1,0},{0,1}};
  this->computeSegmentDisplacements(uMatrixLeft,this->internalForces, leftActiveSegmentLength, leftSegmentLength, u, Tt );  
  ///
  FloatMatrix Tb,Ttb;
  Node *nodeB = this->giveNode(2);
  u.at(1) = nodeB->giveCoordinate(1)+ub.at(1);
  u.at(2) = nodeB->giveCoordinate(3)+ub.at(2);
  double alphab = M_PI + ub.at(3);
  //construct_T(T, alphab, M_PI);
  construct_T(Tb, alphab);
  Ttb.beTranspositionOf(Tb);
  // transformation of forces and moment
  FloatArray f;
  transform_fab2fba(ub, this->internalForces, f);
  this->computeSegmentDisplacements(uMatrixRight,f, rightActiveSegmentLength, rightSegmentLength, u, Ttb );
  */
  //@tobedeleted
  Node *nodeA = this->giveNode(1);
  //
  const int numCellNodes  = 2; // linear line
  //
  int numCellsLeft = ceil(leftActiveSegmentLength/DX);
  int nNodesLeft = numCellsLeft * numCellNodes;
  vtkPieces.at(0).setNumberOfCells(numCellsLeft);
  vtkPieces.at(0).setNumberOfNodes(nNodesLeft);
  //
  int numCellsRight = ceil(rightActiveSegmentLength/DX);
  int nNodesRight = numCellsRight * numCellNodes;
  vtkPieces.at(1).setNumberOfCells(numCellsRight);
  vtkPieces.at(1).setNumberOfNodes(nNodesRight);
  //
  /// Piece 0
  int val    = 1;
  int offset = 0;
  IntArray nodes(numCellNodes);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, uab);
  int nodeNum = 1;
  FloatArray nodeCoords(3);
  IntArray connectivity(2);
  for ( int iElement = 1; iElement <= numCellsLeft; iElement++ ) {
    for (int iNode = 1; iNode <= numCellNodes; iNode++) {
      //double L = (iElement-1) * DX + (iNode -1) * DX;
      nodeCoords.at(1) = nodeA->giveCoordinate(1);// + L * cosAlpha;
      nodeCoords.at(2) = 0;
      nodeCoords.at(3) = nodeA->giveCoordinate(3);// + L * sinAlpha;
      //
      vtkPieces.at(0).setNodeCoords(nodeNum, nodeCoords);
      nodeNum++;
      connectivity.at(iNode) = val++;
    }
    vtkPieces.at(0).setConnectivity(iElement, connectivity);
    offset += 2;
    vtkPieces.at(0).setOffset(iElement, offset);
    vtkPieces.at(0).setCellType(iElement, 3);
  }
  /// Piece 1
  nodeNum = 1;
  val    = 1;
  offset = 0;
  for ( int iElement = 1; iElement <= numCellsRight; iElement++ ) {
    for (int iNode = 1; iNode <= numCellNodes; iNode++) {
      //double L = (iElement-1) * DX + (iNode -1) * DX;
      nodeCoords.at(1) = 0;//nodeB->giveCoordinate(1) + L * cos(pitch + M_PI);
      nodeCoords.at(2) = 0;
      nodeCoords.at(3) = 0;//nodeB->giveCoordinate(3) + L * sin(pitch + M_PI);
      //
      vtkPieces.at(1).setNodeCoords(nodeNum, nodeCoords);
      nodeNum++;
      connectivity.at(iNode) = val++;
    }
    vtkPieces.at(1).setConnectivity(iElement, connectivity);
    offset += 2;
    vtkPieces.at(1).setOffset(iElement, offset);
    vtkPieces.at(1).setCellType(iElement, 3);
  }


  

    int n = primaryVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfPrimaryVarsToExport(n, nNodesLeft);
    vtkPieces [ 1 ].setNumberOfPrimaryVarsToExport(n, nNodesRight);
    for ( int i = 1; i <= n; i++ ) {
        UnknownType utype = ( UnknownType ) primaryVarsToExport.at(i);
        if ( utype == DisplacementVector ) {
	  for ( int nN = 1; nN <= nNodesLeft; nN++ ) {
	    int lN = nN % 2;
	    int iNode;
	    if(lN == 0) {
	      iNode = nN/2 + 1;
	    } else {
	      iNode = (nN + 1) / 2;
	    }

	    FloatArray u(3);
	    u.at(1) = uMatrixLeft.at(iNode, 1);
	    u.at(3) = uMatrixLeft.at(iNode, 2);
	    vtkPieces.at(0).setPrimaryVarInNode(i, nN, u);
	  }
	  for ( int nN = 1; nN <= nNodesRight; nN++ ) {
	    int lN = nN % 2;
	    int iNode;
	    if(lN == 0) {
	      iNode = nN/2 + 1;
	    } else {
	      iNode = (nN + 1) / 2;
	    }

	    FloatArray u(3);
	    u.at(1) = uMatrixRight.at(iNode, 1);
	    u.at(3) = uMatrixRight.at(iNode, 2);

	    vtkPieces.at(1).setPrimaryVarInNode(i, nN, u);
	  }
        }
    }
    /*
    FILE *FIDL, *FIDR;

    char fext[100];
    sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );
    std :: string fileName, functionname, temp;
    fileName = this->giveDomain()->giveEngngModel()->giveOutputBaseFileName();
    size_t foundDot;
    foundDot = fileName.rfind(".");
    
    while (foundDot != std :: string :: npos) {
      fileName.replace(foundDot, 1, "_");
      foundDot = fileName.rfind(".");   
    }
    
    fileName += fext;
    temp = fileName;
    size_t backslash = temp.rfind("/");
    if (backslash != std :: string :: npos ) {
      functionname = temp.substr(backslash+1, std :: string :: npos);
    } else {
      functionname = temp;
    }

    std :: string fileNameL(fileName), fileNameR(fileName);    
    fileNameL += "Left.out";
    fileNameR += "Right.out";
    if ( ( FIDL = fopen(fileNameL.c_str(), "w") ) == NULL ) {
      OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }
    if ( ( FIDR = fopen(fileNameR.c_str(), "w") ) == NULL ) {
      OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }
    uMatrixLeft.printYourself(FIDL);
    uMatrixRight.printYourself(FIDR);

    fclose(FIDL);
    fclose(FIDR);  
    */

}

void
NlBeamInternalContact ::computeSegmentDisplacements(FloatMatrix &uMatrix, const FloatArray &fab, double Lb, double segmentLength, const FloatArray &u0, const FloatMatrix &T, const FloatArray &uab)
{
  FloatArray ub, uc, uplot;
  bool inflection_detected = false;
  double aux;
  
  // initialization at the left end
  double Xab = fab.at(1);
  double Zab = fab.at(2);
  double Mab = fab.at(3);
  double M = -Mab;
  double kappa = computeCurvatureFromMoment(M);

  double x = 0.;
  FloatArray u(3), u_prev(3);
  // basic loop over spatial steps
  int nstep = ceil(Lb/DX); // the spatial step size is fixed and the segment length does not need to be its integer multiple
  uMatrix.resize(nstep+1, 2);


  FloatArray ul(3), l;
  FloatMatrix Tg;

  ul.at(1) = u0.at(1) + T.at(1,1)*(x+u.at(1)) + T.at(1,2)*u.at(2);
  ul.at(2) = u0.at(2) + T.at(2,1)*(x+u.at(1)) + T.at(2,2)*u.at(2);

  this->construct_l(l, uab.at(3), x);
  this->construct_T(Tg, uab.at(3));
  FloatArray ug;
  ul.subtract(l);
  ug.beTProductOf(Tg, ul);
  uMatrix.at(1,1) = ug.at(1) + uab.at(1);
  uMatrix.at(1,2) = ug.at(2) + uab.at(2);


  
  for (int istep=1; istep<=nstep; istep++){
    x += DX;
    u_prev = u;
    // rotation at midstep and its derivatives with respect to the left-end forces
    double phi_mid = u_prev.at(3) + kappa * DX / 2.;
    // normal force at midstep 
    double N_mid = -Xab*cos(phi_mid)+Zab*sin(phi_mid);
    // horizontal displacement at the end of the step
    u.at(1) = u_prev.at(1) + DX * ((1.+N_mid/EA)*cos(phi_mid)-1.);
   // vertical displacement at the end of the step
    u.at(2) = u_prev.at(2) - DX * ( 1. + N_mid / EA ) * sin(phi_mid);
    // bending moment and curvature at the end of the step 
    double M_prev = M; // (store the moment at the beginning of the step, needed for the inflection check)
    M = -Mab+Xab*u.at(2)-Zab*(x+u.at(1));
    kappa = computeCurvatureFromMoment(M);
    // rotation at the end of the step
    u.at(3) = phi_mid+kappa*DX/2.;
    // test whether inflection occurs (the first occurence is considered)
    if (!inflection_detected && M*M_prev<=0. && M_prev!=M){
      // inflection point is detected - output parameter Lc, uc will be computed      inflection_detected = true;
    }

    // test whether the end of the segment has been reached
    if (istep==nstep){
      // displacements and rotation at Lb by linear interpolation within the last step
      aux = Lb/DX - (nstep-1);
      ub.resize(3);
      for (int i=1; i<=3; i++) {
	ub.at(i) = u_prev.at(i) + aux*(u.at(i)-u_prev.at(i));
      }
      /*      uMatrix.at(istep+1, 1) = u0.at(1) + T.at(1,1)*(Lb+ub.at(1)) + T.at(1,2)*ub.at(2);
      uMatrix.at(istep+1, 2) = u0.at(2) + T.at(2,1)*(Lb+ub.at(1)) + T.at(2,2)*ub.at(2);
      */
      ul.at(1) = u0.at(1) + T.at(1,1)*(Lb+ub.at(1)) + T.at(1,2)*ub.at(2);
      ul.at(2) = u0.at(2) + T.at(2,1)*(Lb+ub.at(1)) + T.at(2,2)*ub.at(2);
      // plot also the straight segment behind the contact point, if it exists
      double Lstraight = segmentLength - Lb;

      if (Lstraight>0.){
	/*uMatrix.at(istep+1, 1) = u0.at(1) + T.at(1,1)*(Lb+ub.at(1)+Lstraight*cos(ub.at(3))) + T.at(1,2)*(ub.at(2)-Lstraight*sin(ub.at(3))) ;
	uMatrix.at(istep+1, 2) = u0.at(2) + T.at(2,1)*(Lb+ub.at(1)+Lstraight*cos(ub.at(3))) + T.at(2,2)*(ub.at(2)-Lstraight*sin(ub.at(3)));
	*/
	ul.at(1) = u0.at(1) + T.at(1,1)*(Lb+ub.at(1)+Lstraight*cos(ub.at(3))) + T.at(1,2)*(ub.at(2)-Lstraight*sin(ub.at(3))) ;
	ul.at(2) = u0.at(2) + T.at(2,1)*(Lb+ub.at(1)+Lstraight*cos(ub.at(3))) + T.at(2,2)*(ub.at(2)-Lstraight*sin(ub.at(3)));
      }
    } else {
      ul.at(1) = u0.at(1) + T.at(1,1)*(x+u.at(1)) + T.at(1,2)*u.at(2);
      ul.at(2) = u0.at(2) + T.at(2,1)*(x+u.at(1)) + T.at(2,2)*u.at(2);
      /*uMatrix.at(istep+1, 1) = u0.at(1) + T.at(1,1)*(x+u.at(1)) + T.at(1,2)*u.at(2);
      uMatrix.at(istep+1, 2) = u0.at(2) + T.at(2,1)*(x+u.at(1)) + T.at(2,2)*u.at(2);
      */     
    }
    this->construct_l(l, uab.at(3), x);
    this->construct_T(Tg, uab.at(3));
    ul.subtract(l);
    ug.beTProductOf(Tg, ul);
    uMatrix.at(istep+1,1) = ug.at(1) + uab.at(1);
    uMatrix.at(istep+1,2) = ug.at(2) + uab.at(2);

  } // end of loop over spatial steps
  
}

  
}

