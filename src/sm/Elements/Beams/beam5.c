#include "../sm/Elements/Beams/nlbeam_sm.h"
#include "material.h"
#include "crosssection.h"

#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "engngm.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(NlBeam_SM);


NlBeam_SM :: NlBeam_SM(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    numberOfDofMans = 2;

    length = 0.;
    pitch = 10.;  // a dummy value
    numberOfIntegrationPoints = 4;

}


void
NlBeam_SM :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_w, R_v
    };
}


 double
NlBeam_SM :: givePitch()
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
        pitch  = atan2(yB - yA, xB - xA);
    }

    return pitch;
}

 IRResultType
NlBeam_SM :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    // first call parent
    StructuralElement :: initializeFrom(ir);

    // Numerical parameters
    // 1. number of segments for numerical integration along the beam, default value 100
    IR_GIVE_OPTIONAL_FIELD(ir, nIntegrationPoints, _IFT_FbarElementExtensionInterface_fbarflag);
    // relative tolerance for iterations at the section level
    IR_GIVE_OPTIONAL_FIELD(ir, section_tol, _IFT_FbarElementExtensionInterface_fbarflag);
    // maximum number of iterations at the section level
    IR_GIVE_OPTIONAL_FIELD(ir, section_maxit, _IFT_FbarElementExtensionInterface_fbarflag);

    // relative tolerance for iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_tol, _IFT_FbarElementExtensionInterface_fbarflag);
    // maximum number of iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_maxit, _IFT_FbarElementExtensionInterface_fbarflag);

    return IRRT_OK;
}


 

 

double
NlBeam_SM :: computeMomentFromCurvature(double kappa)
{
   return EI*kappa;
}

double
NlBeam_SM :: computeDerMomentFromCurvature(double kappa)
{
  return EI;
}

double
NlBeam_SM :: computeCurvatureFromMoment(double M)
{
  return M/EI;
}


void
NlBeam_SM :: integrateAlongBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi)
{
  FloatArray x(NN+1), u(NN+1), w(NN+1), phi(NN+1);
  double du(3), dw[3], dphi[3], dM[3], dkappa[3], dphi_mid[3], dN_mid[3];
  FloatArray du(3), dw(3), dphi(3), dM(3), dkappa(3), dphi_mid(3), dN_mid(3);
  double Xab = fab.at(1), Zab = fab.at(2), Mab = fab.at(3);


  
  double dx = beamLength/NN;
  for (int i=2; i <= NN+1; i++) {
    x.at(i) = x.at(i-1) + dx;
    double M = -Mab + Xab*w.at(i) - Zab * (x.at(i) + u.at(i));
    dM = {w.at(i-1)+Xab*dw.at(1)-Zab*du.at(1), -x.at(i) - u.at(i)+Xab*dw.at(2)-Zab*du.at(2), -1.+Xab*dw.at(3) - Zab*du.at(3)};
    double kappa = computeCurvatureFromMoment(M);
    double dMdkappa = computeDerMomentFromCurvature(kappa);
    dkappa = dM/dMdkappa;
    double phi_mid = phi.at(i) + kappa * dx/2.;
    dphi_mid = dphi+ dkappa * dx/2.;
    double N_mid = -Xab*cos(phi_mid)+Zab*sin(phi_mid);
    dN_mid = ( Xab * sin(phi_mid) + Zab * cos(phi_mid)) * dphi_mid;
    dN_mid.at(1) -= cos(phi_mid);
    dN_mid.at(2) += sin(phi_mid);
    u.at(i) = u.at(i-1) + dx * ( ( 1. + N_mid / EA ) * cos(phi_mid) - 1. );
    du += dx*(dN_mid/EA) * cos(phi_mid) - dx * (1. + N_mid/EA) * sin(phi_mid) * dphi_mid;
    w.at(i) = w.at(i-1) -  dx * ( 1. + N_mid / EA ) * sin(phi_mid);
    dw.at -= dx * (dN_mid / EA) * sin(phi_mid) + dx * ( 1. + N_mid / EA) * cos(phi_mid) * dphi_mid;
    M = -Mab+Xab * w.at(i) - Zab* ( x.at(i) + u.at(i) );
    dM = w.at(i) + Xab*dw.at(1) - Zab * du.at(1), -x.at(i)-u.at(i) + Xab * dw.at(2) - Zab * du.at(2), -1. + Xab * dw.at(3) - Zab * du.at(3);
    kappa = computeCurvatureFromMoment(M);
    dMdkappa = computeDerMomentFromCurvature(kappa);
    dkappa = dM/dMdkappa;
    phi.at(i) = phi_mid + kappa * dx / 2.;
    dphi = dphi_mid +dkappa * dx/2.;
  }
  ub.at(1) = u.at(NN+1); ub.at(2) = w.at(NN+1); ub.at(3) = phi.at(NN+1);
  for (j=0; j<=2; j++) {jacobi[0][j]=du[j]; jacobi[1][j]=dw[j]; jacobi[2][j]=dphi[j];}
}


/*
Find forces and moment at the left end that lead to given displacements and rotations
at the right end, assuming that the displacement and rotation at the left end are zero.
The solution is computed iteratively, using the Newton-Raphson method. 
Everything is done here in the local coordinate system aligned with the undeformed beam.
*/
bool
NlBeam_SM :: findLeftEndForcesLocal(FloatArray &ub_target, FloatArray &fab_loc)
{
  double ub, wb, phib;
  FloatArray res, dforces, ub_loc;
  FloatMatrix jacobi(3,3);
  int iter = 0;
  double tolerance = beam_tol* ub_target.computeNorm();
  this->integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, jacobi);
  res = ub_target-ub_loc; 
  double error = res.computeNorm();
  
  while ((iter==0) || (error>tolerance && iter < beam_maxit)){
    iter++;
    jacobi.solveForRhs(res, dforces, );
    fab_loc.add(dforces);
    this->integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, jacobi);
    res = ub_target - ub_loc; 
    error = res.computeNorm();
  }

  if (iter > beam_maxit) {
    //@todo: cut the step
    printf("No convergence in findLeftEndForcesLocal\n");
    return false;
  }
  return true;
}

/*
Auxiliary matrix for transformation of displacements (or in transposed form of forces).
It corresponds to rotation from global axes to the auxiliary system aligned with the deformed beam.
*/
void
  NlBeam_SM :: construct_T(FloatMatrix &T, const double phia)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = cos(alpha-phia);
  T.at(1,2) = sin(alpha-phia);
  T.at(2,1) = - T.at(1,2);
  T.at(3,3) = 1.;
}

/*
Auxiliary matrix for transformation of stiffness.
It is obtained by differentiating T with respect to phia.
*/
void
  NlBeam_SM :: construct_Tprime(FloatMatrix &T, const double phia)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = sin(alpha-phia);
  T.at(1,2) = -cos(alpha-phia);
  T.at(2,1) = -T.at(1,2);
}

void
NlBeam_SM :: construct_l(FloatArray &l, double phia)
{
  l.resize(3);
  l.at(1) = beamLength * (cos(phia)-1.);
  l.at(2) = beamLength * sin(phia);
}

void
NlBeam_SM :: construct_lprime(FloatArray &l, const double phia)
{
  l.resize(3);
  l.at(1) = -beamLength * sin(phia);
  l.at(2) = beamLength * cos(phia);
}


/*
Find forces and moment at the left end that lead to given displacements and rotations
at both ends.
The input displacements and output forces are taken with respect to the global coordinate system.
Note that the transformation matrix T is affected by angle alpha that specifies the initial beam geometry. 
*/
bool
NlBeam_SM :: findLeftEndForces(const FloatArray &u, FloatArray &fab)
{
  FloatArray ub_loc(3), l, fab_loc(3); 
  FloatMatrix T;
  double phia = u.at(3);
  
  // compute displacements of the right end with respect to the auxiliary coordinate system
  construct_l(ub_loc, phia);
  construct_T(T, phia); 
  for (i=0; i<=2; i++)
    for (j=0; j<=2; j++)
      ub_loc[i] += T[i][j]*(u[3+j]-u[j]);

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
NlBeam_SM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  // solution vector
  this->computeVectorOf({D_u, D_v, D_w}, modeType, tStep, u); 
  // only the first three entries of f are computed
  this->findLeftEndForces(u, f); 
  f.at(4) = -f.at(1);
  f.at(5) = -f.at(2);
  double c1 = beamLength*sin(alpha) + u.at(5) - u.at(2);
  double c2 = -beamLength*cos(alpha) - u.at(4) + u.at(1);
  f.at(6) = c1 * f.at(1) + c2 * f.at(2) - f.at(3);
}


/*
Evaluate the tangent stiffness matrix based on given displacements at both ends
and given end forces at the left end (must be provided, which means that typically findLeftEndForces is run first).
*/
void
NlBeam_SM :: computeStiffnessMatrix(FloatMatrix &answer, double u[6], double fab[3])
{

  FloatArray lprime, fab_loc, ub_loc, Tu;
  FloatMatrix T, Tprime, G, Ginv, TtGinv;
  // compute auxiliary matrices
  construct_T(T, u.at(3)); 
  construct_Tprime(Tprime, u.at(3)); 
  construct_lprime(lprime, u.at(3));
  // transform left-end forces to local coordinates
  fab_loc.beProductOf(T, fab);
  // get Jacobi matrix in local coordinates (ub_loc is dummy, will not be used)
  integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, G);
  Ginv.beInverseOf(G)
  // compute product Ttransposed*Ginverse
  TtGinv.beTProductOf(T,Ginv);
  // compute product Ttransposed*Ginverse*T and store it in the upper stiffness block 
  FloatMatrix aux;
  aux.beProductTOf(TtGinv, T);
  answer.setSubMatrix(aux, 1,1);
  answer.setSubMatrix(aux, 1,4);
  
  // compute Tprime*(ub-ua)
  FloatArray u_b, u_a;
  u_a.beSubArrayOf(u, {1,2,3});
  u_b.beSubArrayOf(u, {4,5,6});
  u_ba.beDiffecrenceOf(u_b, u_a);
  Tu.beProductOf(Tprime, u_ba);
  // compute additional stiffness terms in the third column (associated with phia)
  FloatArray col3_1,col3_2;
  col3_1.beProductOf(TtGinv,(Tu[j]+lprime[j]));
  col3_2.beProductOf(TtGinv,Tprime[j][i]*fab_loc[j]);

  // copy minus first row into fourth row and minus second row into fifth row
  for(j = 1; j <= 6; j++){
  }
  
  // construct the sixth row (complete formula)
  double c1 = beamLength*sin(alpha) + u.at(5) - u.at(2);
  double c2 = -beamLength*cos(alpha) - u.at(4) + u.at(1);
  K.at(6,1) = fab[1];
  K.at(6,2) = -fab[0];
  K.at(6,4) = -fab[1];
  K.at(6,5) = fab[0];
 
  for(j=1; j<=6; j++) {
    K.at(4,j) = -K.at(1,j);
    K.at(5,j) = -K.at(2,j) + c1*K.at(1,k) + c2*K.at(2,j) - K.at(3,j);
  }

}
  
     

} // end namespace oofem
