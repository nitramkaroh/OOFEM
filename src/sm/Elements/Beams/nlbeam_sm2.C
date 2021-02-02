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
	//pitch  = atan2(yB + this->eval_w0(curvedbeamLength) - yA - this->eval_w0(0) , xB + this->eval_u0(curvedbeamLength)  - xA - this->eval_u0(0) );
	pitch  = atan2(yB + w_0.at(NIP+1) - yA - w_0.at(1) , xB + u_0.at(NIP+1) - xA - u_0.at(1) );
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
    /// input by array of points
    IR_GIVE_OPTIONAL_FIELD(ir, arcLengthCoordinate, _IFT_NlBeam_SM2_alCoord);

    u_0.resize(NIP+1);
    w_0.resize(NIP+1);
    phi_0.resize(NIP+1);
    phi_12.resize(NIP+1);
    kappa_0.resize(NIP+1);
    
    if(arcLengthCoordinate.giveSize()) {
      //here I assume the parabolic arch for now
      pu_0 = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
      pw_0 = {0.000000, 0.089849, 0.356869, 0.793846, 1.389885, 2.131722, 3.005029, 3.995474, 5.089454, 6.274517, 7.539542, 8.874760, 10.271677, 11.722956, 13.222276, 14.764197, 16.344035, 17.957748, 19.601839, 21.273274, 22.969411, 24.687944, 26.426848, 28.184344, 29.958859, 31.748999, 33.553524, 35.371330, 37.201425, 39.042921, 40.895018, 42.756993, 44.628191, 46.508017, 48.395930, 50.291436, 52.194084, 54.103459, 56.019182, 57.940901, 59.868295, 61.801065, 63.738935, 65.681649, 67.628968, 69.580673, 71.536556, 73.496426, 75.460102, 77.427417, 79.398213, 81.372342, 83.349667, 85.330055, 87.313385, 89.299540, 91.288411, 93.279896, 95.273897, 97.270321, 99.269082, 101.270096, 103.273284, 105.278573, 107.285891, 109.295171, 111.306348, 113.319362, 115.334154, 117.350668, 119.368852, 121.388655, 123.410027, 125.432924, 127.457300, 129.483113, 131.510323, 133.538890, 135.568776, 137.599946, 139.632366, 141.666002, 143.700822, 145.736796, 147.773894, 149.812088, 151.851350, 153.891654, 155.932975, 157.975288, 160.018570, 162.062797, 164.107948, 166.154002, 168.200937, 170.248734, 172.297374, 174.346839, 176.397109, 178.448169, 180.500000};
      pw_0.times(-1);
      pphi_0 = {0.000000, 0.042403, 0.084580, 0.126308, 0.167385, 0.207625, 0.246867, 0.284980, 0.321858, 0.357425, 0.391631, 0.424446, 0.455864, 0.485896, 0.514563, 0.541902, 0.567956, 0.592772, 0.616403, 0.638904, 0.660331, 0.680740, 0.700184, 0.718718, 0.736392, 0.753256, 0.769356, 0.784737, 0.799440, 0.813505, 0.826969, 0.839867, 0.852232, 0.864093, 0.875480, 0.886419, 0.896935, 0.907052, 0.916790, 0.926171, 0.935213, 0.943935, 0.952353, 0.960482, 0.968337, 0.975932, 0.983280, 0.990392, 0.997280, 1.003955, 1.010426, 1.016703, 1.022794, 1.028708, 1.034454, 1.040037, 1.045465, 1.050745, 1.055883, 1.060885, 1.065755, 1.070501, 1.075126, 1.079635, 1.084033, 1.088325, 1.092513, 1.096602, 1.100595, 1.104497, 1.108310, 1.112038, 1.115682, 1.119248, 1.122736, 1.126150, 1.129492, 1.132765, 1.135970, 1.139111, 1.142188, 1.145205, 1.148163, 1.151063, 1.153908, 1.156699, 1.159438, 1.162126, 1.164765, 1.167356, 1.169901, 1.172401, 1.174858, 1.177271, 1.179644, 1.181976, 1.184269, 1.186523, 1.188741, 1.190923, 1.193069, 1.195181, 1.197260, 1.199306, 1.201320, 1.203304, 1.205257, 1.207181, 1.209076, 1.210943, 1.212783, 1.214596, 1.216382, 1.218144, 1.219880, 1.221592, 1.223280, 1.224945, 1.226587, 1.228207, 1.229805, 1.231382, 1.232938, 1.234473, 1.235989, 1.237485, 1.238962, 1.240420, 1.241860, 1.243281, 1.244686, 1.246072, 1.247442, 1.248796, 1.250133, 1.251454, 1.252759, 1.254049, 1.255324, 1.256585, 1.257831, 1.259062, 1.260280, 1.261484, 1.262674, 1.263851, 1.265016, 1.266167, 1.267306, 1.268433, 1.269547, 1.270650, 1.271741, 1.272821, 1.273889, 1.274947, 1.275993, 1.277029, 1.278054, 1.279069, 1.280073, 1.281068, 1.282053, 1.283028, 1.283993, 1.284949, 1.285896, 1.286834, 1.287763, 1.288683, 1.289594, 1.290497, 1.291392, 1.292278, 1.293156, 1.294025, 1.294887, 1.295742, 1.296588, 1.297427, 1.298258, 1.299083, 1.299899, 1.300709, 1.301512, 1.302308, 1.303096, 1.303879, 1.304654, 1.305423, 1.306185, 1.306941, 1.307691, 1.308435, 1.309172, 1.309904, 1.310629, 1.311349, 1.312062, 1.312770, 1.313473};
      pkappa_0 =  {0.040000, 0.039573, 0.038346, 0.036471, 0.034149, 0.031585, 0.028955, 0.026386, 0.023963, 0.021731, 0.019706, 0.017888, 0.016268, 0.014828, 0.013551, 0.012418, 0.011411, 0.010517, 0.009719, 0.009007, 0.008368, 0.007795, 0.007279, 0.006812, 0.006390, 0.006006, 0.005656, 0.005337, 0.005045, 0.004777, 0.004531, 0.004304, 0.004094, 0.003900, 0.003720, 0.003553, 0.003397, 0.003252, 0.003117, 0.002990, 0.002871, 0.002760, 0.002656, 0.002557, 0.002465, 0.002377, 0.002295, 0.002217, 0.002143, 0.002073, 0.002007, 0.001944, 0.001884, 0.001827, 0.001773, 0.001721, 0.001672, 0.001625, 0.001580, 0.001537, 0.001496, 0.001457, 0.001419, 0.001383, 0.001348, 0.001315, 0.001283, 0.001253, 0.001223, 0.001195, 0.001167, 0.001141, 0.001116, 0.001091, 0.001068, 0.001045, 0.001023, 0.001002, 0.000981, 0.000961, 0.000942, 0.000924, 0.000906, 0.000888, 0.000871, 0.000855, 0.000839, 0.000824, 0.000809, 0.000794, 0.000780, 0.000766, 0.000753, 0.000740, 0.000728, 0.000716, 0.000704, 0.000692, 0.000681, 0.000670, 0.000659};
      
      //IR_GIVE_FIELD(ir, beta, _IFT_NlBeam_SM2_beta);
      IR_GIVE_FIELD(ir, curvedbeamLength, _IFT_NlBeam_SM2_curvedbeamLength);
      
      /*
      IR_GIVE_FIELD(ir, pu_0, _IFT_NlBeam_SM2_pu0);
      IR_GIVE_FIELD(ir, pw_0, _IFT_NlBeam_SM2_pw0);
      IR_GIVE_FIELD(ir, pphi_0, _IFT_NlBeam_SM2_pphi0);
      IR_GIVE_FIELD(ir, pkappa_0, _IFT_NlBeam_SM2_pkappa0);
     
      if(pu_0.giveSize() != NIP+1 || pw_0.giveSize() != NIP+1 || pphi_0.giveSize () != 2.*NIP+1 || pkappa_0.giveSize() != NIP+1) {
	OOFEM_ERROR("Inconsistent size of the input fields....");
      }
      */
      cosBeta = (curvedbeamLength+ pu_0.at(NIP+1))/beamLength;
      sinBeta = pw_0.at(NIP+1)/beamLength;

      
      
      u_0 = pu_0;
      w_0 = pw_0; 
      kappa_0 = pkappa_0;
      for(int i = 1; i <= NIP+1; i++) {
	phi_0.at(i) = pphi_0.at(2*i-1);
	if(i < NIP+1) {
	  phi_12.at(i) = pphi_0.at(2*i);
	}
      }
      
    } else if(tangentPoint.giveSize()) {
      IR_GIVE_FIELD(ir, fu_0, _IFT_NlBeam_SM2_fu0);
      IR_GIVE_FIELD(ir, fw_0, _IFT_NlBeam_SM2_fw0);
      IR_GIVE_FIELD(ir, fphi_0, _IFT_NlBeam_SM2_fphi0);
      IR_GIVE_FIELD(ir, fkappa_0, _IFT_NlBeam_SM2_fkappa0);
      
      FloatArray tangentLine(tangentPoint);
      tangentLine.subtract(pointA);
      curvedbeamLength = tangentLine.computeNorm();
      cosBeta = (curvedbeamLength+ this->eval_u0(curvedbeamLength))/beamLength;
      sinBeta = this->eval_w0(curvedbeamLength)/beamLength;
  
      double dx = curvedbeamLength/NIP;
      double x = 0;
      for(int i = 1; i <= NIP+1; i++) {
	u_0.at(i) = this->eval_u0(x);
	w_0.at(i) = this->eval_w0(x);
	phi_0.at(i) = this->eval_phi0(x);
	kappa_0.at(i) = this->eval_kappa0(x);	
	x +=dx;
	phi_12.at(i) = this->eval_phi0(x - dx/2.);
      }
      
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
  if(fkappa_0.isDefined()) {
    return (& fkappa_0)->eval( { { "x", x } }, this->giveDomain() );
  } else {
    return 0;
  }
  
}


  
double
NlBeam_SM2 :: eval_phi0(double x)
{  
  if(fphi_0.isDefined()) {
    return (& fphi_0)->eval( { { "x", x } }, this->giveDomain() );
  } else {
    return 0;
  }
  
}
double
NlBeam_SM2 :: eval_u0(double x)
{
  if(fu_0.isDefined()) {
    return (&fu_0)->eval( { { "x", x } }, this->giveDomain() );
  } else {
    return 0;
  }
}
double
NlBeam_SM2 :: eval_w0(double x)
{
  if(fw_0.isDefined()) {
    return (&fw_0)->eval( { { "x", x } }, this->giveDomain() );
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
  double delta_kappa = computeDeltaCurvatureFromInternalForces(vM.at(1),vN.at(1), kappa_0.at(1));
  double dkappa_dM = computeDerCurvatureMoment(vM.at(1),vN.at(1), kappa_0.at(1));
  double dkappa_dN = computeDerCurvatureNormalForce(vM.at(1),vN.at(1), kappa_0.at(1));
  dkappa = dM;
  dkappa.times(dkappa_dM);
  dkappa.add(dkappa_dN, dN);                                                      
  this->x.at(1) = 0;
  ub = {0., 0., 0.};
  jacobi.resize(3,3);
  // basic loop over spatial steps
  for (int i=2; i<=NIP+1; i++){
    double ds;
    if(arcLengthCoordinate.giveSize()) {
      ds = arcLengthCoordinate.at(i) - arcLengthCoordinate.at(i-1);
    } else {
      ds = curvedbeamLength/NIP;
    }
  this->x.at(i) = this->x.at(i-1) + ds;
  u_prev = ub;
  // rotation at midstep and its derivatives with respect to the left-end forces
  //double phi0_mid = eval_phi0(x.at(i)-ds/2.);
  double phi0_mid = phi_12.at(i-1);
  double delta_phi_mid = u_prev.at(3) + delta_kappa * ds/2.;
  double phi_mid = phi0_mid + delta_phi_mid;
  FloatMatrix jacobiprime;
  jacobiprime.beTranspositionOf(jacobi);
  dphi_mid.beColumnOf(jacobiprime, 3);
  dphi_mid.add(ds/2, dkappa);  
  // normal force at midstep and its derivatives with respect to the left-end forces
  double N_mid = -Xab * cos(phi_mid) + Zab *  sin(phi_mid);
  vN.at(i) = N_mid;
  dN = ( Xab * sin(phi_mid) + Zab * cos(phi_mid)) * dphi_mid;
  dN.at(1) -= cos(phi_mid);
  dN.at(2) += sin(phi_mid);
  // centerline strain at midstep
  double eps_mid = computeCenterlineStrainFromInternalForces(vM.at(i-1),N_mid, kappa_0.at(i));
  double deps_dM = computeDerStrainMoment(vM.at(i-1),N_mid, kappa_0.at(i));
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
  double M = -Mab+Xab*(w_0.at(i)+ub.at(2))-Zab*(x.at(i)+ u_0.at(i)+ub.at(1));
  vM.at(i) = M;
  jacobiprime.beTranspositionOf(jacobi);
  dM.beColumnOf(jacobiprime, 2);
  dM.times(Xab);
  FloatArray aux;
  aux.beColumnOf(jacobiprime, 1);
  dM.add(-Zab, aux);
  //double tdM1 = eval_w0(x.at(i))+ub.at(2);
  //double tdM2 = -(x.at(i)+eval_u0(x.at(i))+ub.at(1));
  dM.at(1) += w_0.at(i)+ub.at(2);
  dM.at(2) += -(x.at(i)+u_0.at(i)+ub.at(1));
 
  dM.at(3) += -1.;
  delta_kappa = computeDeltaCurvatureFromInternalForces(M,N_mid, kappa_0.at(i));
  dkappa_dM = computeDerCurvatureMoment(M,N_mid, kappa_0.at(i));
  dkappa_dN = computeDerCurvatureNormalForce(M,N_mid, kappa_0.at(i));
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
    //double chLength = sqrt((L + eval_u0(L))*(L + eval_u0(L)) + eval_w0(L) * eval_w0(L));
    double chLength = sqrt((L + u_0.at(NIP+1))*(L + u_0.at(NIP+1)) + w_0.at(NIP+1) * w_0.at(NIP+1));
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
