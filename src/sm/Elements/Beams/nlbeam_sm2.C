#include "../sm/Elements/Beams/nlbeam_sm2.h"
#include "material.h"
#include "crosssection.h"
#include "node.h"

#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "engngm.h"
#include "mathfem.h"
#include "classfactory.h"
#include "function.h"

namespace oofem {
REGISTER_Element(NlBeam_SM2);


NlBeam_SM2 :: NlBeam_SM2(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    this->internalForces.resize(3);

}


void
NlBeam_SM2 :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_w, R_v
    };
}



double
NlBeam_SM2 :: computeLength()
// Returns the length of the receiver.
{
    double dx, dy;
    Node *nodeA, *nodeB;

    if ( beamLength == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        beamLength  = sqrt(dx * dx + dy * dy);
    }

    return beamLength;
}
  

 double
NlBeam_SM2 :: givePitch()
// Returns the pitch of the receiver.
{
    double xA, xB, yA, yB;
    Node *nodeA, *nodeB;

    if ( pitch == 10. ) {             // 10. : dummy initialization value
        nodeA  = this->giveNode(1);
        nodeB  = this->giveNode(2);
        xA     = nodeA->giveCoordinate(1);
        xB     = nodeB->giveCoordinate(1);
        yA     = nodeA->giveCoordinate(3);
        yB     = nodeB->giveCoordinate(3);
	pitch  = atan2(yB + this->eval_w0(curvedbeamLength) - yA - this->eval_w0(0) , xB + this->eval_u0(curvedbeamLength)  - xA - this->eval_u0(0) );
    }

    return pitch;
}

 IRResultType
NlBeam_SM2 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    // first call parent
    StructuralElement :: initializeFrom(ir);

    // Numerical parameters
    // 1. number of segments for numerical integration along the beam, default value 100
    IR_GIVE_OPTIONAL_FIELD(ir, NIP, _IFT_NlBeam_SM2_NIP);
    IR_GIVE_FIELD(ir, EA, _IFT_NlBeam_SM2_EA);
    IR_GIVE_FIELD(ir, EI, _IFT_NlBeam_SM2_EI);
    IR_GIVE_OPTIONAL_FIELD(ir, RADIUS, _IFT_NlBeam_SM2_RADIUS);
    IR_GIVE_OPTIONAL_FIELD(ir, DEPTH, _IFT_NlBeam_SM2_DEPTH);
    
    this->computeLength();

    cosBeta = 1;
    sinBeta = 0;

    auto nodeA  = this->giveNode(1);
    auto nodeB  = this->giveNode(2);
    FloatArray pointA({nodeA->giveCoordinate(1), nodeA->giveCoordinate(3)});
    FloatArray pointB({nodeB->giveCoordinate(1), nodeB->giveCoordinate(3)});
       
    tangentPoint.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, tangentPoint, _IFT_NlBeam_SM2_tangentPoint);    
    
    if(tangentPoint.giveSize()) {
      IR_GIVE_FIELD(ir, u_0, _IFT_NlBeam_SM2_u0);
      IR_GIVE_FIELD(ir, w_0, _IFT_NlBeam_SM2_w0);
      IR_GIVE_FIELD(ir, phi_0, _IFT_NlBeam_SM2_phi0);
      IR_GIVE_FIELD(ir, kappa_0, _IFT_NlBeam_SM2_kappa0);
      FloatArray tangentLine(tangentPoint);
      tangentLine.subtract(pointA);
      curvedbeamLength = tangentLine.computeNorm();
      cosBeta = (curvedbeamLength+ this->eval_u0(curvedbeamLength))/beamLength;
      sinBeta = this->eval_w0(curvedbeamLength)/beamLength;
    }
   
    cosAlpha = (pointB.at(1) - pointA.at(1))/beamLength;
    sinAlpha = (pointB.at(2) - pointA.at(2))/beamLength;

    //check consistency
    /*    double cosGamma = (tangentPoint.at(1)-pointA.at(1))/curvedbeamLength;
    double sinGamma = (tangentPoint.at(2)-pointA.at(2))/curvedbeamLength;
    double x_x0 = tangentPoint.at(1) + this->eval_u0(curvedbeamLength) * cosGamma + this->eval_w0(curvedbeamLength) * sinGamma;
    double y_y0 = tangentPoint.at(2) - this->eval_u0(curvedbeamLength) * sinGamma + this->eval_w0(curvedbeamLength) * cosGamma;
    */
    
    /*    // relative tolerance for iterations at the section level
    IR_GIVE_OPTIONAL_FIELD(ir, section_tol, _IFT_FbarElementExtensionInterface_fbarflag);
    // maximum number of iterations at the section level
    IR_GIVE_OPTIONAL_FIELD(ir, section_maxit, _IFT_FbarElementExtensionInterface_fbarflag);
    */
    
    // relative tolerance for iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_tol, _IFT_NlBeam_SM2_Beam_Tolerance);
    // maximum number of iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_maxit, _IFT_NlBeam_SM2_Beam_MaxIteration);

    
    this->x.resize(NIP+1);
    this->u.resize(NIP+1);
    this->w.resize(NIP+1);
    this->phi.resize(NIP+1);

    this->vN.resize(NIP+1);
    this->vV.resize(NIP+1);
    this->vM.resize(NIP+1);

    
    return IRRT_OK;
}


 

 

double
NlBeam_SM2 :: computeMomentFromCurvature(double kappa)
{
   return EI*kappa;
}

double
NlBeam_SM2 :: computeDerMomentFromCurvature(double kappa)
{
  return EI;
}

double
NlBeam_SM2 :: computeCurvatureFromMoment(double M)
{
  return M/EI;
}



/*
Functions defining the initial stress-free shape
*/
double
NlBeam_SM2 :: eval_kappa0(double x)
{  
  if(kappa_0.isDefined()) {
    return (& kappa_0)->eval( { { "x", x } }, this->giveDomain() );
  } else {
    return 0;
  }
  
}


  
double
NlBeam_SM2 :: eval_phi0(double x)
{  
  if(phi_0.isDefined()) {
    return (& phi_0)->eval( { { "x", x } }, this->giveDomain() );
  } else {
    return 0;
  }
  
}
double
NlBeam_SM2 :: eval_u0(double x)
{
  if(u_0.isDefined()) {
    return (&u_0)->eval( { { "x", x } }, this->giveDomain() );
  } else {
    return 0;
  }
}
double
NlBeam_SM2 :: eval_w0(double x)
{
  if(w_0.isDefined()) {
    return (&w_0)->eval( { { "x", x } }, this->giveDomain() );
  } else {
    return 0;
  }
}
/*
Functions defining the relations between internal forces and deformation variables
*/
double 
NlBeam_SM2 :: computeDeltaCurvatureFromInternalForces(double M, double N)
{
  double eps = (N+M/RADIUS)/EA;
  double delta_kappa = eps/RADIUS + M/(EI*(1.+0.15*DEPTH*DEPTH/(RADIUS*RADIUS)));
  return delta_kappa;
}
double 
NlBeam_SM2 :: computeDerCurvatureMoment(double M, double N)
{
  return 1./(RADIUS*RADIUS*EA) + 1./(EI*(1.+0.15*DEPTH*DEPTH/(RADIUS*RADIUS)));
}
double 
NlBeam_SM2 :: computeDerCurvatureNormalForce(double M, double N)
{
  return 1./(RADIUS*EA);
}
double 
NlBeam_SM2 :: computeCenterlineStrainFromInternalForces(double M, double N)
{
  return (N+M/RADIUS)/EA;
}
double 
NlBeam_SM2 :: computeDerStrainMoment(double M, double N)
{
  return 1./(RADIUS*EA);
}
double 
NlBeam_SM2 :: computeDerStrainNormalForce(double M, double N)
{
  return 1./EA;
}


void
NlBeam_SM2 :: integrateAlongBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi)
{
  if(tangentPoint.giveSize()) {
    this->integrateAlongCurvedBeamAndGetJacobi(fab, ub, jacobi);
   } else {
    this->integrateAlongStraightBeamAndGetJacobi(fab, ub, jacobi);
   }
}

void
NlBeam_SM2 :: integrateAlongCurvedBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi)
{
  // initialization at the left end
  ub.resize(3);
  double Xab = fab.at(1), Zab = fab.at(2), Mab = fab.at(3);
   this->vM.at(1) = -Mab;
  this->vN.at(1) = -Xab;
  FloatArray dM(3), dN(3), dkappa(3), dphi_mid(3), deps_mid(3),
  u_prev(3);
  dM = {0., 0., -1.};
  dN = {-1., 0., 0.};
  double delta_kappa = computeDeltaCurvatureFromInternalForces(vM.at(1),vN.at(1));
  double dkappa_dM = computeDerCurvatureMoment(vM.at(1),vN.at(1));
  double dkappa_dN = computeDerCurvatureNormalForce(vM.at(1),vN.at(1));
  dkappa = dM;
  dkappa.times(dkappa_dM);
  dkappa.add(dkappa_dN, dN);                                                      
  this->x.at(1) = 0;
  ub = {0., 0., 0.};
  jacobi.resize(3,3);
  double dx = curvedbeamLength/NIP;
  // basic loop over spatial steps
  for (int i=2; i<=NIP+1; i++){
  this->x.at(i) = this->x.at(i-1) + dx;
  u_prev = ub;
  // rotation at midstep and its derivatives with respect to the left-end forces
  double phi0_mid = eval_phi0(x.at(i)-dx/2.);
  double delta_phi_mid = u_prev.at(3) + delta_kappa * dx/2.;
  double phi_mid = phi0_mid + delta_phi_mid;
  FloatMatrix jacobiprime;
  jacobiprime.beTranspositionOf(jacobi);
  dphi_mid.beColumnOf(jacobiprime, 3);
  dphi_mid.add(dx/2, dkappa);  
  // normal force at midstep and its derivatives with respect to the left-end forces
  double N_mid = -Xab * cos(phi_mid) + Zab *  sin(phi_mid);
  vN.at(i) = N_mid;
  dN = ( Xab * sin(phi_mid) + Zab * cos(phi_mid)) * dphi_mid;
  dN.at(1) -= cos(phi_mid);
  dN.at(2) += sin(phi_mid);
  // centerline strain at midstep
  double eps_mid = computeCenterlineStrainFromInternalForces(vM.at(i-1),N_mid);
  double deps_dM = computeDerStrainMoment(vM.at(i-1),N_mid);
  double deps_dN = computeDerStrainNormalForce(vM.at(i-1),N_mid);
  deps_mid = dM;
  deps_mid.times(deps_dM);
  deps_mid.add(deps_dN, dN);   							
  // horizontal displacement at the end of the step
  ub.at(1) = u_prev.at(1)+dx*((1.+eps_mid)*cos(phi_mid)-cos(phi0_mid));
  this->u.at(i) = ub.at(1);
  // vertical displacement at the end of the step
  ub.at(2) = u_prev.at(2)+dx*(sin(phi0_mid)-(1.+eps_mid)*sin(phi_mid));
  this->w.at(i) = ub.at(2);
  // first two rows jacobi
  FloatArray j_row1, j_row2;
  double s1, s2;
  s1 = (1.+eps_mid)*(-sin(phi_mid));
  j_row1.beScaled(s1, dphi_mid);
  j_row1.add(cos(phi_mid), deps_mid);
  j_row1.times(dx);
  s2 = (1.+eps_mid)*(-cos(phi_mid));
  j_row2.beScaled(s2, dphi_mid);
  j_row2.add(-sin(phi_mid), deps_mid);
  j_row2.times(dx);
  jacobi.addSubVectorRow(j_row1, 1,1);
  jacobi.addSubVectorRow(j_row2, 2,1);
  // bending moment and curvature at the end of the step and their derivatives with respect to the left-end forces
  double M = -Mab+Xab*(eval_w0(x.at(i))+ub.at(2))-Zab*(x.at(i)+eval_u0(x.at(i))+ub.at(1));
  vM.at(i) = M;
  jacobiprime.beTranspositionOf(jacobi);
  dM.beColumnOf(jacobiprime, 2);
  dM.times(Xab);
  FloatArray aux;
  aux.beColumnOf(jacobiprime, 1);
  dM.add(-Zab, aux);
  dM.at(1) += eval_w0(x.at(i))+ub.at(2);
  dM.at(2) += -(x.at(i)+eval_u0(x.at(i))+ub.at(1));
  dM.at(3) += -1.;
  delta_kappa = computeDeltaCurvatureFromInternalForces(M,N_mid);
  dkappa_dM = computeDerCurvatureMoment(M,N_mid);
  dkappa_dN = computeDerCurvatureNormalForce(M,N_mid);
  dkappa = dM;
  dkappa.times(dkappa_dM);
  dkappa.add(dkappa_dN, dN);                                                      
  // rotation at the end of the step
  ub.at(3) = delta_phi_mid+delta_kappa*dx/2.;
  this->phi.at(i) = ub.at(3);
  // update Jacobi matrix
    FloatArray j_row3;
    j_row3 = dphi_mid;
    j_row3.add(dx/2., dkappa);
    jacobi.copySubVectorRow(j_row3, 3,1);
  } // end of loop over spatial steps
}

void
NlBeam_SM2 :: integrateAlongStraightBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi)
{
  ub.resize(3);
  FloatArray du(3), dw(3), dphi(3), dM(3), dkappa(3), dphi_mid(3), dN_mid(3);
  double Xab = fab.at(1), Zab = fab.at(2), Mab = fab.at(3);
  this->vM.at(1) = -Mab + Xab*w.at(1) - Zab * u.at(1);
  this->vN.at(1) = -Xab * cos(phi.at(1)) + Zab * sin(phi.at(1));
  double dx = beamLength/NIP;
  for (int i=2; i <= NIP+1; i++) {
    this->x.at(i) = this->x.at(i-1) + dx;
    double M = -Mab + Xab*w.at(i-1) - Zab * (x.at(i-1) + u.at(i-1));
    dM = {w.at(i-1)+Xab*dw.at(1)-Zab*du.at(1), -x.at(i-1) - u.at(i-1)+Xab*dw.at(2)-Zab*du.at(2), -1.+Xab*dw.at(3) - Zab*du.at(3)};
    double kappa = computeCurvatureFromMoment(M);
    double dMdkappa = computeDerMomentFromCurvature(kappa);
    dkappa = dM;
    dkappa.times(1./dMdkappa);
    double phi_mid = phi.at(i-1) + kappa * dx/2.;
    dphi_mid = dphi;
    dphi_mid.add(dx/2 ,dkappa);
    double N_mid = -Xab * cos(phi_mid) + Zab * sin(phi_mid);
    vN.at(i) = N_mid;
    dN_mid = ( Xab * sin(phi_mid) + Zab * cos(phi_mid)) * dphi_mid;
    dN_mid.at(1) -= cos(phi_mid);
    dN_mid.at(2) += sin(phi_mid);
    u.at(i) = u.at(i-1) + dx * ( ( 1. + N_mid / EA ) * cos(phi_mid) - 1. );
    du.add((dx/EA) * cos(phi_mid),dN_mid);
    du.add(- dx * (1. + N_mid/EA) * sin(phi_mid), dphi_mid);
    w.at(i) = w.at(i-1) -  dx * ( 1. + N_mid / EA ) * sin(phi_mid);
    dw.add( -dx / EA * sin(phi_mid),dN_mid);
    dw.add( -dx * ( 1. + N_mid / EA) * cos(phi_mid), dphi_mid);
    M = -Mab + Xab * w.at(i) - Zab * ( x.at(i) + u.at(i) );
    vM.at(i) = M;
    dM = {w.at(i) + Xab*dw.at(1) - Zab * du.at(1), -x.at(i)-u.at(i) + Xab * dw.at(2) - Zab * du.at(2), -1. + Xab * dw.at(3) - Zab * du.at(3)};
    kappa = computeCurvatureFromMoment(M);
    dMdkappa = computeDerMomentFromCurvature(kappa);
    dkappa = dM;
    dkappa.times(1./dMdkappa);
    phi.at(i) = phi_mid + kappa * dx / 2.;
    dphi = dphi_mid;
    dphi.add(dx/2 ,dkappa);
  }
  ub.at(1) = u.at(NIP+1); ub.at(2) = w.at(NIP+1); ub.at(3) = phi.at(NIP+1);

  jacobi.copySubVectorRow(du, 1,1);
  jacobi.copySubVectorRow(dw, 2,1);
  jacobi.copySubVectorRow(dphi, 3,1);

}


/*
Find forces and moment at the left end that lead to given displacements and rotations
at the right end, assuming that the displacement and rotation at the left end are zero.
The solution is computed iteratively, using the Newton-Raphson method. 
Everything is done here in the local coordinate system aligned with the undeformed beam.
*/
void
NlBeam_SM2 :: findLeftEndForcesLocal(FloatArray &ub_target, FloatArray &fab_loc)
{
  FloatArray fab_init(fab_loc);
  FloatArray res, dforces, ub_loc;
  FloatMatrix jacobi(3,3);
  int iter = 0;
  double tolerance = beam_tol* ub_target.computeNorm();
  this->integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, jacobi);
  res = ub_target-ub_loc; 
  double error = res.computeNorm();
  
  while ((iter==0) || (error>tolerance && iter < beam_maxit && error ==error) ){
    iter++;
    jacobi.solveForRhs(res, dforces);
    fab_loc.add(dforces);
    this->integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, jacobi);
    res = ub_target - ub_loc; 
    error = res.computeNorm();
  }

  if (iter >= beam_maxit || error != error) {
    //@todo: cut the step
    
    domain->giveEngngModel()->setAnalysisCrash(true);
    fab_loc = fab_init;
    //OOFEM_ERROR("No convergence in findLeftEndForcesLocal\n");
  }
}

/*
Auxiliary matrix for transformation of displacements (or in transposed form of forces).
It corresponds to rotation from global axes to the auxiliary system aligned with the deformed beam.
*/

void
NlBeam_SM2 :: construct_T(FloatMatrix &T, const double phia)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = (cosAlpha*cosBeta + sinAlpha*sinBeta)*cos(phia) + (sinAlpha*cosBeta-cosAlpha*sinBeta)*sin(phia);
  T.at(1,2) =  (sinAlpha*cosBeta-cosAlpha*sinBeta) *cos(phia) - (cosAlpha*cosBeta+sinAlpha*sinBeta) * sin(phia);


  T.at(2,1) = - T.at(1,2);
  T.at(3,3) = 1.;
}

/*
Auxiliary matrix for transformation of stiffness.
It is obtained by differentiating T with respect to phia.
*/


void
  NlBeam_SM2 :: construct_Tprime(FloatMatrix &T, const double phia)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = (sinAlpha*cosBeta - cosAlpha * sinBeta) *cos(phia) - (cosAlpha*cosBeta+sinAlpha*sinBeta)*sin(phia);
  T.at(1,2) = -(cosAlpha*cosBeta + sinAlpha*sinBeta)*cos(phia) - (sinAlpha*cosBeta - cosAlpha*sinBeta)*sin(phia);

  T.at(2,1) = -T.at(1,2);
}


void
NlBeam_SM2 :: construct_l(FloatArray &l, double phia)
{
  l.resize(3);
  l.at(1) = beamLength * (cos(phia)*cosBeta - sin(phia) * sinBeta - cosBeta);
  l.at(2) = beamLength * (sinBeta * cos(phia) + cosBeta * sin(phia) - sinBeta);
}


void
NlBeam_SM2 :: construct_l(FloatArray &l, double phia, double L)
{
  l.resize(3);
  l.at(1) = L * (cos(phia)*cosBeta - sin(phia) * sinBeta - cosBeta);
  l.at(2) = L * (sinBeta * cos(phia) + cosBeta * sin(phia) - sinBeta);
 
}



void
NlBeam_SM2 :: construct_lprime(FloatArray &l, const double phia)
{
  l.resize(3);
  l.at(1) = -beamLength * (sin(phia)*cosBeta + cos(phia) * sinBeta);
  l.at(2) =  beamLength * (cosBeta * cos(phia) -  sinBeta * sin(phia));
}



/*
Find forces and moment at the left end that lead to given displacements and rotations
at both ends.
The input displacements and output forces are taken with respect to the global coordinate system.
Note that the transformation matrix T is affected by angle alpha that specifies the initial beam geometry. 
*/
void
NlBeam_SM2 :: findLeftEndForces(const FloatArray &u, FloatArray &fab)
{
  FloatArray ub_loc, l, fab_loc; 
  FloatMatrix T; 
  // compute displacements of the right end with respect to the auxiliary coordinate system
  construct_l(ub_loc, u.at(3));
  construct_T(T, u.at(3));

  FloatArray u_b, u_a, u_ba, temp;
  u_a.beSubArrayOf(u, {1,2,3});
  u_b.beSubArrayOf(u, {4,5,6});
  u_ba.beDifferenceOf(u_b, u_a);
  temp.beProductOf(T,u_ba);
  ub_loc.add(temp);
  // transform initial guess to local coordinates
  fab_loc.beProductOf(T, fab);
  // find end forces in the local coordinate syste
  findLeftEndForcesLocal(ub_loc, fab_loc);
  // transform local end forces to the global coordinate syste
  fab.beTProductOf(T, fab_loc);
}

/*
Find end forces and moments that lead to given displacements and rotations
at both ends.
The input displacements and output forces are taken with respect to the global coordinate system.
*/

 void
NlBeam_SM2 :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  // solution vector
  answer.resize(6);
  FloatArray u, f_a(3);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
  // only the first three entries of f are computed
  this->findLeftEndForces(u, this->internalForces);
  answer.at(1) = this->internalForces.at(1);
  answer.at(2) = this->internalForces.at(2);
  answer.at(3) = this->internalForces.at(3);  
  answer.at(4) = -answer.at(1);
  answer.at(5) = -answer.at(2);
  double c1 =  beamLength*sinAlpha + u.at(5) - u.at(2);
  double c2 = -beamLength*cosAlpha - u.at(4) + u.at(1);
  answer.at(6) = c1 * answer.at(1) + c2 * answer.at(2) - answer.at(3);
}


 void
 NlBeam_SM2 :: giveInternalForcesVector_from_u(FloatArray &answer, TimeStep *tStep, const FloatArray &u)
{
  // solution vector
  answer.resize(6);
  FloatArray f_a(3);
  // only the first three entries of f are computed
  this->findLeftEndForces(u, this->internalForces);
  answer.at(1) = this->internalForces.at(1);
  answer.at(2) = this->internalForces.at(2);
  answer.at(3) = this->internalForces.at(3);  
  answer.at(4) = -answer.at(1);
  answer.at(5) = -answer.at(2);
  double c1 =  beamLength*sinAlpha + u.at(5) - u.at(2);
  double c2 = -beamLength*cosAlpha - u.at(4) + u.at(1);
  answer.at(6) = c1 * answer.at(1) + c2 * answer.at(2) - answer.at(3);
}

void
NlBeam_SM2 :: computeStiffnessMatrix_num(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
  // solution vector
  FloatArray u, iF, iFn;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
  this->giveInternalForcesVector_from_u(iF, tStep, u);
  double eps = 1.e-8;
  answer.resize(6,6);
    
  FloatArray du;
  /////////////
  for(int i = 1; i <= 6; i++) {
    du = u;
    du.at(i) += eps;
    this->giveInternalForcesVector_from_u(iFn, tStep, du);
    ////////////
    answer.at(1,i) = iFn.at(1) - iF.at(1);
    answer.at(2,i) = iFn.at(2) - iF.at(2);
    answer.at(3,i) = iFn.at(3) - iF.at(3);
    answer.at(4,i) = iFn.at(4) - iF.at(4);
    answer.at(5,i) = iFn.at(5) - iF.at(5);
    answer.at(6,i) = iFn.at(6) - iF.at(6);
  }

  answer.times(1./eps);
}




/*
Evaluate the tangent stiffness matrix based on given displacements at both ends
and given end forces at the left end (must be provided, which means that typically findLeftEndForces is run first).
*/
void
NlBeam_SM2 :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

  
  
  FloatArray lprime, fab_loc, ub_loc, Tu;
  FloatMatrix T, Tprime, G, Ginv, TtGinv;
  answer.resize(6,6);
  // solution vector
  FloatArray u;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u); 
  // compute auxiliary matrices
  construct_T(T, u.at(3)); 
  construct_Tprime(Tprime, u.at(3)); 
  construct_lprime(lprime, u.at(3));
  // transform left-end forces to local coordinates
  fab_loc.beProductOf(T, this->internalForces);
  // get Jacobi matrix in local coordinates (ub_loc is dummy, will not be used)
  integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, G);
  Ginv.beInverseOf(G);
  // compute product Ttransposed*Ginverse
  TtGinv.beTProductOf(T,Ginv);
  // compute product Ttransposed*Ginverse*T and store it in the upper stiffness block 
  FloatMatrix aux;
  aux.beProductOf(TtGinv, T);
  answer.setSubMatrix(aux, 1,1);
  answer.times(-1);
  answer.setSubMatrix(aux, 1,4);
  
  // compute Tprime*(ub-ua)
  FloatArray u_b, u_a, u_ba;
  u_a.beSubArrayOf(u, {1,2,3});
  u_b.beSubArrayOf(u, {4,5,6});
  u_ba.beDifferenceOf(u_b, u_a);
  Tu.beProductOf(Tprime, u_ba);
  // compute additional stiffness terms in the third column (associated with phia)
  FloatArray col3_1,col3_2;
  Tu.add(lprime);
  col3_1.beProductOf(TtGinv,Tu);
  col3_2.beTProductOf(Tprime, fab_loc);
  answer.addSubVectorCol(col3_1, 1,3);
  answer.addSubVectorCol(col3_2, 1,3);
  // construct the sixth row (complete formula)
  double c1 = beamLength*sinAlpha + u.at(5) - u.at(2);
  double c2 = -beamLength*cosAlpha - u.at(4) + u.at(1);
  answer.at(6,1) =  this->internalForces.at(2);
  answer.at(6,2) = -this->internalForces.at(1);
  answer.at(6,4) = -this->internalForces.at(2);
  answer.at(6,5) =  this->internalForces.at(1);
 
  for(int j = 1; j <= 6; j++) {
    answer.at(4,j) = -answer.at(1,j);
    answer.at(5,j) = -answer.at(2,j);
    answer.at(6,j) += c1 * answer.at(1,j) + c2 * answer.at(2,j) - answer.at(3,j);
  }

  /*
  FloatMatrix K;
  this->computeStiffnessMatrix_num(K, rMode, tStep);
  int huhu = 1;
  */

}
  



void
NlBeam_SM2 :: printOutputAt(FILE *file, TimeStep *tStep)
{
  if(tangentPoint.giveSize()) {
    this-> printOutputAt_CurvedBeam(file, tStep);
   } else {
    this-> printOutputAt_StraightBeam(file, tStep);
   }
}


void
NlBeam_SM2 :: printOutputAt_StraightBeam(FILE *file, TimeStep *tStep)
{

  FILE *FID;
  

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


  
  fileName += ".m";
  if ( ( FID = fopen(fileName.c_str(), "w") ) == NULL ) {
    OOFEM_ERROR("failed to open file %s", fileName.c_str() );
  }



  //transform displacements to global coordinate system
  FloatArray uab, ug(NIP+1), wg(NIP+1), phig(NIP+1), u_l, u_g;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, uab);
  ug.at(1) = uab.at(1);
  wg.at(1) = uab.at(2);
  phig.at(1) = uab.at(3);
  double L = 0;
  double dx = beamLength/NIP;
  for (int i=2; i <= NIP+1; i++) {
    L = L + dx;
    FloatArray l;
    FloatMatrix T;
    this->construct_T(T, uab.at(3));
    this->construct_l(l, uab.at(3), L);
    u_l = {this->u.at(i), this->w.at(i), 0};
    u_l.subtract(l);
    u_g.beTProductOf(T, u_l);
    ug.at(i) = u_g.at(1) + ug.at(1);
    wg.at(i) = u_g.at(2) + wg.at(1);
    phig.at(i) = this->phi.at(i) + phig.at(1);
    
  }

  
  

  fprintf( FID, " function [x z u w phi N M]= %s \n", functionname.c_str() );

   fprintf(FID, "x=[");
   for ( double val: x ) {
     fprintf( FID, "%f,", val * cosAlpha );
   }   
   fprintf(FID, "];\n");
   
   fprintf(FID, "z=[");
   for ( double val: x ) {
     fprintf( FID, "%f,", val * sinAlpha );
   }   
   fprintf(FID, "];\n");

   
   
   fprintf(FID, "u=[");
   for ( double val: ug ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");
   
   
   fprintf(FID, "w=[");
   for ( double val: wg ) {
     fprintf( FID, "%f,", val );
   }     
   fprintf(FID, "];\n");
   
   fprintf(FID, "phi=[");
   for ( double val: phig ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");


   fprintf(FID, "N=[");
   for ( double val: vN ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");


   fprintf(FID, "M=[");
   for ( double val: vM ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");

   
   fprintf(FID, "end\n");

   fclose(FID);

}


void
NlBeam_SM2 :: printOutputAt_CurvedBeam(FILE *file, TimeStep *tStep)
{

  FILE *FID;
  

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


  
  fileName += ".m";
  if ( ( FID = fopen(fileName.c_str(), "w") ) == NULL ) {
    OOFEM_ERROR("failed to open file %s", fileName.c_str() );
  }



  //transform displacements to global coordinate system
  FloatArray uab, ug(NIP+1), wg(NIP+1), phig(NIP+1), u_l, u_g;
  FloatArray xg(NIP+1), zg(NIP+1);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, uab);
  ug.at(1) = uab.at(1);
  wg.at(1) = uab.at(2);
  phig.at(1) = uab.at(3);
  double L = 0;
  double dx = curvedbeamLength/NIP;
  for (int i=2; i <= NIP+1; i++) {
    L = L + dx;
    FloatArray l;
    FloatMatrix T;
    double chLength = sqrt((L + eval_u0(L))*(L + eval_u0(L)) + eval_w0(L) * eval_w0(L));
    this->construct_T(T, uab.at(3));
    this->construct_l(l, uab.at(3),  chLength);
    
    // u_l = {this->u.at(i)-eval_u0(L), this->w.at(i)-eval_w0(L), 0};


    xg.at(i) = L + this->eval_u0(L);
    zg.at(i) = this->eval_w0(L);

    u_l = {this->u.at(i), this->w.at(i), 0};
    u_l.subtract(l);
    u_g.beTProductOf(T, u_l);    
    ug.at(i) = u_g.at(1) + ug.at(1);
    wg.at(i) = u_g.at(2) + wg.at(1);
    //phig.at(i) = this->phi.at(i) -eval_phi0(L)  + phig.at(1);
    phig.at(i) = this->phi.at(i)  + phig.at(1);
  }

 
  

  fprintf( FID, " function [x z u w phi N M]= %s \n", functionname.c_str() );

   fprintf(FID, "x=[");
   for ( double val: xg ) {
     fprintf( FID, "%f,", val);
   }   
   fprintf(FID, "];\n");
   
   fprintf(FID, "z=[");
   for ( double val: zg ) {
     fprintf( FID, "%f,", val );
   }   
   fprintf(FID, "];\n");

   
   
   fprintf(FID, "u=[");
   for ( double val: ug ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");
   
   
   fprintf(FID, "w=[");
   for ( double val: wg ) {
     fprintf( FID, "%f,", val );
   }     
   fprintf(FID, "];\n");
   
   fprintf(FID, "phi=[");
   for ( double val: phig ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");


   fprintf(FID, "N=[");
   for ( double val: vN ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");


   fprintf(FID, "M=[");
   for ( double val: vM ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");

   
   fprintf(FID, "end\n");

   fclose(FID);

}


FILE *
NlBeam_SM2 :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;

    char fext[100];
    sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );
    std :: string fileName;
    fileName = this->giveDomain()->giveEngngModel()->giveOutputBaseFileName();
    fileName += fext;
    fileName += ".m";
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }
    return answer;
}
  

} // end namespace oofem
