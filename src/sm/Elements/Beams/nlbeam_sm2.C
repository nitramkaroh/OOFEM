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
	//pitch  = atan2(yB + w_0.at(NIP+1) - yA - w_0.at(1) , xB + u_0.at(NIP+1) - xA - u_0.at(1) );
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
    tangentVector.clear();  
    IR_GIVE_OPTIONAL_FIELD(ir, tangentVector, _IFT_NlBeam_SM2_tangentVector);
      
    if(tangentVector.giveSize()) {
      IR_GIVE_FIELD(ir, u_0, _IFT_NlBeam_SM2_u0);
      IR_GIVE_FIELD(ir, w_0, _IFT_NlBeam_SM2_w0);
      IR_GIVE_FIELD(ir, phi_0, _IFT_NlBeam_SM2_phi0);
      IR_GIVE_FIELD(ir, kappa_0, _IFT_NlBeam_SM2_kappa0);
      IR_GIVE_FIELD(ir, curvedbeamLength, _IFT_NlBeam_SM2_curvedbeamLength);
      int n_cf = 0;
      IR_GIVE_OPTIONAL_FIELD(ir, n_cf,  _IFT_NlBeam_SM2_coordinateFlag);

      this->s.resize(NIP+1);
      this->u0.resize(NIP+1);
      this->w0.resize(NIP+1);
      this->phi0.resize(NIP+1);
      this->kappa0.resize(NIP+1);
      this->phi0mid.resize(NIP);

      
      if(n_cf == 0) {
	this->cf = CF_s;
	for(int i = 1; i <= NIP+1; i++) {
	  s.at(i) = double(i-1.)/NIP * curvedbeamLength;
	  kappa0.at(i) = this->eval_kappa0(s.at(i));
	  phi0.at(i) = this->eval_phi0(s.at(i));
	  u0.at(i) = this->eval_u0(s.at(i));
	  w0.at(i) = this->eval_w0(s.at(i));
	  double ds = curvedbeamLength/NIP;
	  if(i > 1) {
	    phi0mid.at(i-1) = this->eval_phi0(s.at(i)-ds/2.);
	  }
	}
      } else {
	this->cf = CF_x;
	IR_GIVE_FIELD(ir, sx, _IFT_NlBeam_SM2_s);

	for(int i = 1; i <= NIP+1; i++) {
	  //////////////////////
	  //	  double sd = double(i-1.)/NIP * curvedbeamLength;
	  //double x = 100 * sin(sd/100);
	  double x = double(i-1.)/NIP * fabs(nodeB->giveCoordinate(1)-nodeA->giveCoordinate(1));
	  //  double test = x - asinh(2.*0.02*x)/4./0.02 - x*sqrt(4.*0.02*0.02*x*x)/2.;
	  
	  s.at(i) = this->eval_s(x) ;
	  kappa0.at(i) = this->eval_kappa0(x);
	  phi0.at(i) = this->eval_phi0(x);
	  u0.at(i) = this->eval_u0(x);
	  w0.at(i) = this->eval_w0(x);
	  double dx = 1./NIP * nodeB->giveCoordinate(1);
	  if(i > 1) {
	    phi0mid.at(i-1) = this->eval_phi0(x-dx);
	  }

	  
	}
      
      }

      IR_GIVE_FIELD(ir, EI, _IFT_NlBeam_SM2_EI);
      cosBeta = (curvedbeamLength+ u0.at(NIP+1))/beamLength;
      sinBeta = w0.at(NIP+1)/beamLength;      
    }
   
    cosAlpha = (pointB.at(1) - pointA.at(1))/beamLength;
    sinAlpha = (pointB.at(2) - pointA.at(2))/beamLength;

    
    // relative tolerance for iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_tol, _IFT_NlBeam_SM2_Beam_Tolerance);
    // maximum number of iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_maxit, _IFT_NlBeam_SM2_Beam_MaxIteration);

    
    this->s.resize(NIP+1);
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


double
NlBeam_SM2 :: eval_s(double x)
{  
  if(sx.isDefined()) {
    return (& sx)->eval( { { "x", x } }, this->giveDomain() );
  } else {
    return 0;
  }
  
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
NlBeam_SM2 :: computeDeltaCurvatureFromInternalForces(double M, double N, double curvature)
{
  double eps = (N+M*curvature)/EA;
  double delta_kappa = eps*curvature + M/(EI*(1.+0.15*DEPTH*DEPTH*curvature*curvature));
  return delta_kappa;
}
double 
NlBeam_SM2 :: computeDerCurvatureMoment(double M, double N, double curvature)
{
  return curvature * curvature/(EA) + 1./(EI*(1.+0.15*DEPTH*DEPTH*(curvature*curvature)));
}
double 
NlBeam_SM2 :: computeDerCurvatureNormalForce(double M, double N, double curvature)
{
  return curvature/(EA);
}
double 
NlBeam_SM2 :: computeCenterlineStrainFromInternalForces(double M, double N, double curvature)
{
  return (N+M*curvature)/EA;
}
double 
NlBeam_SM2 :: computeDerStrainMoment(double M, double N, double curvature)
{
  return curvature/(EA);
}
double 
NlBeam_SM2 :: computeDerStrainNormalForce(double M, double N)
{
  return 1./EA;
}


void
NlBeam_SM2 :: integrateAlongBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi)
{
  if(tangentVector.giveSize()) {
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
  double delta_kappa = computeDeltaCurvatureFromInternalForces(vM.at(1),vN.at(1), kappa0.at(1));
  double dkappa_dM = computeDerCurvatureMoment(vM.at(1),vN.at(1), kappa0.at(1));
  double dkappa_dN = computeDerCurvatureNormalForce(vM.at(1),vN.at(1), kappa0.at(1));
  dkappa = dM;
  dkappa.times(dkappa_dM);
  dkappa.add(dkappa_dN, dN);                                                      
  this->s.at(1) = 0;
  ub = {0., 0., 0.};
  jacobi.resize(3,3);
  // basic loop over spatial steps
  for (int i=2; i<=NIP+1; i++){
    double ds = this->s.at(i) - this->s.at(i-1);
      //    this->s.at(i) = this->s.at(i-1) + ds.at(i-1);
    u_prev = ub;
    // rotation at midstep and its derivatives with respect to the left-end forces
    //double phi0_mid = eval_phi0(x.at(i)-ds/2.);
    double phi0_mid = phi0mid.at(i-1);
    double delta_phi_mid = u_prev.at(3) + delta_kappa * ds/2.;
    double phi_mid = phi0_mid + delta_phi_mid;
    FloatMatrix jacobiprime;
    jacobiprime.beTranspositionOf(jacobi);
    dphi_mid.beColumnOf(jacobiprime, 3);
    dphi_mid.add(ds/2., dkappa);  
    // normal force at midstep and its derivatives with respect to the left-end forces
    double N_mid = -Xab * cos(phi_mid) + Zab *  sin(phi_mid);
    vN.at(i) = N_mid;
    dN = ( Xab * sin(phi_mid) + Zab * cos(phi_mid)) * dphi_mid;
    dN.at(1) -= cos(phi_mid);
    dN.at(2) += sin(phi_mid);
    // centerline strain at midstep
    double eps_mid = computeCenterlineStrainFromInternalForces(vM.at(i-1),N_mid, kappa0.at(i));
    double deps_dM = computeDerStrainMoment(vM.at(i-1),N_mid, kappa0.at(i));
    double deps_dN = computeDerStrainNormalForce(vM.at(i-1),N_mid);
    deps_mid = dM;
    deps_mid.times(deps_dM);
    deps_mid.add(deps_dN, dN);   							
    // horizontal displacement at the end of the step
    ub.at(1) = u_prev.at(1)+ds*((1.+eps_mid)*cos(phi_mid)-cos(phi0_mid));
    this->u.at(i) = ub.at(1);
    // vertical displacement at the end of the step
    ub.at(2) = u_prev.at(2)+ds*(sin(phi0_mid)-(1.+eps_mid)*sin(phi_mid));
    this->w.at(i) = ub.at(2);
    // first two rows jacobi
    FloatArray j_row1, j_row2;
    double s1, s2;
    s1 = (1.+eps_mid)*(-sin(phi_mid));
    j_row1.beScaled(s1, dphi_mid);
    j_row1.add(cos(phi_mid), deps_mid);
    j_row1.times(ds);
    s2 = (1.+eps_mid)*(-cos(phi_mid));
    j_row2.beScaled(s2, dphi_mid);
    j_row2.add(-sin(phi_mid), deps_mid);
    j_row2.times(ds);
    jacobi.addSubVectorRow(j_row1, 1,1);
    jacobi.addSubVectorRow(j_row2, 2,1);
    // bending moment and curvature at the end of the step and their derivatives with respect to the left-end forces
    //double M = -Mab+Xab*(eval_w0(x.at(i))+ub.at(2))-Zab*(x.at(i)+eval_u0(x.at(i))+ub.at(1));
    double M = -Mab+Xab*(w0.at(i)+ub.at(2))-Zab*(s.at(i)+ u0.at(i)+ub.at(1));
    vM.at(i) = M;
    jacobiprime.beTranspositionOf(jacobi);
    dM.beColumnOf(jacobiprime, 2);
    dM.times(Xab);
    FloatArray aux;
    aux.beColumnOf(jacobiprime, 1);
    dM.add(-Zab, aux);
    //dM.at(1) = eval_w0(x.at(i))+ub.at(2);
    //dM.at(2) = -(x.at(i)+eval_u0(x.at(i))+ub.at(1));
    dM.at(1) += w0.at(i)+ub.at(2);
    dM.at(2) += -(s.at(i)+u0.at(i)+ub.at(1)); 
    dM.at(3) += -1.;
    delta_kappa = computeDeltaCurvatureFromInternalForces(M,N_mid, kappa0.at(i));
    dkappa_dM = computeDerCurvatureMoment(M,N_mid, kappa0.at(i));
    dkappa_dN = computeDerCurvatureNormalForce(M,N_mid, kappa0.at(i));
    dkappa = dM;
    dkappa.times(dkappa_dM);
    dkappa.add(dkappa_dN, dN);                                                      
    // rotation at the end of the step
    ub.at(3) = delta_phi_mid+delta_kappa*ds/2.;
    this->phi.at(i) = ub.at(3);
    // update Jacobi matrix
    FloatArray j_row3;
    j_row3 = dphi_mid;
    j_row3.add(ds/2., dkappa);
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
    this->s.at(i) = this->s.at(i-1) + dx;
    double M = -Mab + Xab*w.at(i-1) - Zab * (s.at(i-1) + u.at(i-1));
    dM = {w.at(i-1)+Xab*dw.at(1)-Zab*du.at(1), -s.at(i-1) - u.at(i-1)+Xab*dw.at(2)-Zab*du.at(2), -1.+Xab*dw.at(3) - Zab*du.at(3)};
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
    M = -Mab + Xab * w.at(i) - Zab * ( s.at(i) + u.at(i) );
    vM.at(i) = M;
    dM = {w.at(i) + Xab*dw.at(1) - Zab * du.at(1), -s.at(i)-u.at(i) + Xab * dw.at(2) - Zab * du.at(2), -1. + Xab * dw.at(3) - Zab * du.at(3)};
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
bool
NlBeam_SM2 :: findLeftEndForcesLocal(FloatArray &ub_target, FloatArray &fab_loc)
{
  FloatArray fab_init(fab_loc);
  FloatArray res, dforces, ub_loc;
  FloatMatrix jacobi(3,3);
  int iter = 0;
  double tolerance = beam_tol;
  if(ub_target.computeNorm() != 0) {
    tolerance *= ub_target.computeNorm();
  }
  
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
    fab_loc = fab_init;
    return false;
  }
  return true;
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
  FloatArray ub_loc, l, fab_loc, fab_init(fab); 
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
  bool success = findLeftEndForcesLocal(ub_loc, fab_loc);
  if (!success) { // multi-level substepping (note that fab_loc remained unchanged if the local iteration failed)
    FloatArray ub_loc_substep, ub_loc_intermediate(3);
    FloatMatrix jacobi(3,3);
    this->integrateAlongBeamAndGetJacobi(fab_loc, ub_loc_intermediate, jacobi);
    ub_loc_substep.beDifferenceOf(ub_loc, ub_loc_intermediate);
    int refinement_factor = 4;
    ub_loc_substep.times(1./refinement_factor);
    int isubstep = 0, maxnsubsteps = 1024;
    int nsubsteps = this->nsubsteps_init;
    while (isubstep < nsubsteps && nsubsteps <= maxnsubsteps){
      ub_loc_intermediate.add(ub_loc_substep);
      success = findLeftEndForcesLocal(ub_loc_intermediate, fab_loc);
      if (success){
	isubstep++;
      } else { // further refinement needed
	ub_loc_intermediate.subtract(ub_loc_substep);
	ub_loc_substep.times(1./refinement_factor);
	nsubsteps = isubstep + 4*(nsubsteps-isubstep);
      }
    }
    if (isubstep < nsubsteps){
      fab = fab_init;
      domain->giveEngngModel()->setAnalysisCrash(true);
      return;
    }
  }

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
  double eps = 1.e-10;
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
  FloatArray lprime, fab_loc, ub_loc, Tu, fab;
  FloatMatrix T, Tprime, G, Ginv, TtGinv;
  answer.resize(6,6);
  // solution vector
  FloatArray u;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
  //  this->findLeftEndForces(u, fab);

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
  


Interface *NlBeam_SM2 :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case VTKXMLExportModuleElementInterfaceType:
      return static_cast< VTKXMLExportModuleElementInterface * >( this );
    default:
      return StructuralElement :: giveInterface(it);
    }
}

void
NlBeam_SM2 :: giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
  if(tangentVector.giveSize()) {
    this->giveCompositeExportData_curved(vtkPieces, primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep );
  } else {
    this->giveCompositeExportData_straight(vtkPieces, primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep );
   }

}

void
NlBeam_SM2 :: giveCompositeExportData_curved(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{

  vtkPieces.resize(1);
 
    int numCells = this->NIP;
    const int numCellNodes  = 2; // linear line
    int nNodes = numCells * numCellNodes;

    vtkPieces.at(0).setNumberOfCells(numCells);
    vtkPieces.at(0).setNumberOfNodes(nNodes);

    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);
 
    Node *nodeA = this->giveNode(1); 
    FloatArray uab, ug(3);
    this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, uab);
    
    FloatMatrix T,T0;
    this->construct_T(T0, 0);
    this->construct_T(T, uab.at(3));

    int nodeNum = 1;
    FloatArray nodeCoords(3);
    IntArray connectivity(2);
    for ( int iElement = 1; iElement <= numCells; iElement++ ) {      
      for (int iNode = 1; iNode <= numCellNodes; iNode++) {
	int index = iElement  + iNode - 1;
	double Ls = this->s.at(index);
	FloatArray xs(3), xg;
	xs.at(1) = Ls + this->u0.at(index);
	xs.at(2) = this->w0.at(index);
	xs.at(3) = 0;
	xg.beTProductOf(T0, xs);
	nodeCoords.at(1) = nodeA->giveCoordinate(1) + xg.at(1);
	nodeCoords.at(3) = nodeA->giveCoordinate(3) + xg.at(2);
	vtkPieces.at(0).setNodeCoords(nodeNum, nodeCoords);
	nodeNum++;
	connectivity.at(iNode) = val++;
      }
      vtkPieces.at(0).setConnectivity(iElement, connectivity);
      offset += 2;
      vtkPieces.at(0).setOffset(iElement, offset);
      vtkPieces.at(0).setCellType(iElement, 3);
    }


    int n = primaryVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfPrimaryVarsToExport(n, nNodes);
    for ( int i = 1; i <= n; i++ ) {
        UnknownType utype = ( UnknownType ) primaryVarsToExport.at(i);
        if ( utype == DisplacementVector ) {
	  FloatArray l;
	  for ( int nN = 1; nN <= nNodes; nN++ ) {
	    int lN = nN % 2;
	    int iNode;
	    if(lN == 0) {
	      iNode = nN/2 + 1;
	    } else {
	      iNode = (nN + 1) / 2;
	    }
	  
	    double Ls = s.at(iNode);
	    double L = sqrt((Ls + this->u0.at(iNode)) * (Ls + this->u0.at(iNode)) + this->w0.at(iNode) * this->w0.at(iNode));
	    this->construct_l(l, uab.at(3), L);
	    FloatArray u_l, u_g;
	    u_l = {this->u.at(iNode), this->w.at(iNode), 0};
	    u_l.subtract(l);
	    u_g.beTProductOf(T, u_l);
	    ug.at(1) = u_g.at(1) + uab.at(1);
	    ug.at(3) = u_g.at(2) + uab.at(2);
	    ug.at(2) = 0;
	    vtkPieces.at(0).setPrimaryVarInNode(i, nN, ug);
	  }
        }
    }


}


void
NlBeam_SM2 :: giveCompositeExportData_straight(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
    vtkPieces.resize(1);
 
    int numCells = this->NIP;
    const int numCellNodes  = 2; // linear line
    int nNodes = 2 * numCellNodes;

    vtkPieces.at(0).setNumberOfCells(numCells);
    vtkPieces.at(0).setNumberOfNodes(nNodes);

    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);
 

    FloatArray uab, ug(3);
    this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, uab);
    int nodeNum = 1;
    double ds = beamLength/NIP;
    FloatMatrix T;
    this->construct_T(T, uab.at(3));
    FloatArray nodeCoords;
    IntArray connectivity(2);
    for ( int iElement = 1; iElement <= numCells; iElement++ ) {
      for (int iNode = 1; iNode <= numCellNodes; iNode++) {
	double L = (iElement-1) * ds + (iNode -1) * ds;
	nodeCoords.at(1) = L * cosAlpha;
	nodeCoords.at(2) = L * sinAlpha;      
	vtkPieces.at(0).setNodeCoords(iNode, nodeCoords);
	nodeNum++;
	connectivity.at(iNode) = val++;
      }
      vtkPieces.at(0).setConnectivity(iElement, connectivity);
      offset += 2;
      vtkPieces.at(0).setOffset(iElement, offset);
      vtkPieces.at(0).setCellType(iElement, 3);
    }


    int n = primaryVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfPrimaryVarsToExport(n, nNodes);
    for ( int i = 1; i <= n; i++ ) {
        UnknownType utype = ( UnknownType ) primaryVarsToExport.at(i);
        if ( utype == DisplacementVector ) {
	  FloatArray l;
	  FloatMatrix T;
	  for ( int nN = 1; nN <= nNodes; nN++ ) {
	    int lN = nN % 2;
	    int iNode;
	    if(lN == 0) {
	      iNode = nN/2 + 1;
	    } else {
	      iNode = (nN + 1) / 2;
	    }
	    double L = (iNode-1) * ds;
	    this->construct_l(l, uab.at(3), L);
	    this->construct_T(T, uab.at(3));
	    FloatArray u_l, u_g;
	    u_l = {this->u.at(nN), this->w.at(nN), 0};
	    u_l.subtract(l);
	    u_g.beTProductOf(T, u_l);
	    ug.at(1) = u_g.at(1) + uab.at(1);
	    ug.at(2) = u_g.at(2) + uab.at(2);
	    ug.at(3) = this->phi.at(i) + uab.at(3);
	    vtkPieces.at(0).setPrimaryVarInNode(i, nN, ug);
	  }
        }
    }

}
     


void
NlBeam_SM2 :: printOutputAt(FILE *file, TimeStep *tStep)
{
  if(tangentVector.giveSize()) {
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
   for ( double val: s ) {
     fprintf( FID, "%f,", val * cosAlpha );
   }   
   fprintf(FID, "];\n");
   
   fprintf(FID, "z=[");
   for ( double val: s ) {
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
    //double chLength = sqrt((L + eval_u0(L))*(L + eval_u0(L)) + eval_w0(L) * eval_w0(L));
    double chLength = sqrt((L + u0.at(NIP+1))*(L + u0.at(NIP+1)) + w0.at(NIP+1) * w0.at(NIP+1));
    this->construct_T(T, uab.at(3));
    this->construct_l(l, uab.at(3),  chLength);
    
    // u_l = {this->u.at(i)-eval_u0(L), this->w.at(i)-eval_w0(L), 0};

    //@todo:fix this
    //xg.at(i) = L + this->eval_u0(L);
    //zg.at(i) = this->eval_w0(L);

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
