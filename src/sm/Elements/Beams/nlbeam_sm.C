#include "../sm/Elements/Beams/nlbeam_sm.h"
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
REGISTER_Element(NlBeam_SM);


NlBeam_SM :: NlBeam_SM(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    this->internalForces.resize(3);

}


void
NlBeam_SM :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_w, R_v
    };
}



double
NlBeam_SM :: computeLength()
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
    IR_GIVE_OPTIONAL_FIELD(ir, NIP, _IFT_NlBeam_SM_NIP);
    /*    // relative tolerance for iterations at the section level
    IR_GIVE_OPTIONAL_FIELD(ir, section_tol, _IFT_FbarElementExtensionInterface_fbarflag);
    // maximum number of iterations at the section level
    IR_GIVE_OPTIONAL_FIELD(ir, section_maxit, _IFT_FbarElementExtensionInterface_fbarflag);
    */
    
    // relative tolerance for iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_tol, _IFT_NlBeam_SM_Beam_Tolerance);
    // maximum number of iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_maxit, _IFT_NlBeam_SM_Beam_MaxIteration);

    
    this->givePitch();
    this->computeLength();

    x.resize(NIP+1);
    u.resize(NIP+1);
    w.resize(NIP+1);
    phi.resize(NIP+1);


    vN.resize(NIP+1);
    vV.resize(NIP+1);
    vM.resize(NIP+1);

    
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
  ub.resize(3);
  FloatArray du(3), dw(3), dphi(3), dM(3), dkappa(3), dphi_mid(3), dN_mid(3);
  double Xab = fab.at(1), Zab = fab.at(2), Mab = fab.at(3);
  vM.at(1) = -Mab + Xab*w.at(1) - Zab * u.at(1);
  vN.at(1) = -Xab * cos(phi.at(1)) + Zab * sin(phi.at(1));
  double dx = beamLength/NIP;
  for (int i=2; i <= NIP+1; i++) {
    x.at(i) = x.at(i-1) + dx;
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
NlBeam_SM :: findLeftEndForcesLocal(FloatArray &ub_target, FloatArray &fab_loc)
{
  FloatArray res, dforces, ub_loc;
  FloatMatrix jacobi(3,3);
  int iter = 0;
  double tolerance = beam_tol* ub_target.computeNorm();
  this->integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, jacobi);
  res = ub_target-ub_loc; 
  double error = res.computeNorm();
  
  while ((iter==0) || (error>tolerance && iter < beam_maxit)){
    iter++;
    jacobi.solveForRhs(res, dforces);
    fab_loc.add(dforces);
    this->integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, jacobi);
    res = ub_target - ub_loc; 
    error = res.computeNorm();
  }

  if (iter > beam_maxit) {
    //@todo: cut the step
    OOFEM_ERROR("No convergence in findLeftEndForcesLocal\n");
  }
}

/*
Auxiliary matrix for transformation of displacements (or in transposed form of forces).
It corresponds to rotation from global axes to the auxiliary system aligned with the deformed beam.
*/
void
NlBeam_SM :: construct_T(FloatMatrix &T, const double phia)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = cos(pitch-phia);
  T.at(1,2) = sin(pitch-phia);
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
  T.at(1,1) = T.at(2,2) = sin(pitch-phia);
  T.at(1,2) = -cos(pitch-phia);
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
void
NlBeam_SM :: findLeftEndForces(const FloatArray &u, FloatArray &fab)
{
  FloatArray ub_loc, l, fab_loc; 
  FloatMatrix T; 
  // compute displacements of the right end with respect to the auxiliary coordinate system
  construct_l(ub_loc, u.at(3));
  construct_T(T, u.at(3));

  FloatArray u_b, u_a, u_ba;
  u_a.beSubArrayOf(u, {1,2,3});
  u_b.beSubArrayOf(u, {4,5,6});
  u_ba.beDifferenceOf(u_b, u_a);
  ub_loc.beProductOf(T,u_ba);
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
  answer.resize(6);
  FloatArray u, f_a(3);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
  //u.at(5) = 0.1;
  // only the first three entries of f are computed
  this->findLeftEndForces(u, this->internalForces);
  answer.at(1) = this->internalForces.at(1);
  answer.at(2) = this->internalForces.at(2);
  answer.at(3) = this->internalForces.at(3);  
  answer.at(4) = -answer.at(1);
  answer.at(5) = -answer.at(2);
  double c1 = beamLength*sin(pitch) + u.at(5) - u.at(2);
  double c2 = -beamLength*cos(pitch) - u.at(4) + u.at(1);
  answer.at(6) = c1 * answer.at(1) + c2 * answer.at(2) - answer.at(3);
}


/*
Evaluate the tangent stiffness matrix based on given displacements at both ends
and given end forces at the left end (must be provided, which means that typically findLeftEndForces is run first).
*/
void
NlBeam_SM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
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
  aux.beProductTOf(TtGinv, T);
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
  double c1 = beamLength*sin(pitch) + u.at(5) - u.at(2);
  double c2 = -beamLength*cos(pitch) - u.at(4) + u.at(1);
  answer.at(6,1) =  this->internalForces.at(2);
  answer.at(6,2) = -this->internalForces.at(1);
  answer.at(6,4) = -this->internalForces.at(2);
  answer.at(6,5) =  this->internalForces.at(1);
 
  for(int j = 1; j <= 6; j++) {
    answer.at(4,j) = -answer.at(1,j);
    answer.at(5,j) = -answer.at(2,j);
    answer.at(6,j) += c1 * answer.at(1,j) + c2 * answer.at(2,j) - answer.at(3,j);
  }

}
  

void
NlBeam_SM :: printOutputAt(FILE *file, TimeStep *tStep)
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


  fprintf( FID, " function [x u w phi N M]= %s \n", functionname.c_str() );
  

   fprintf(FID, "x=[");
   for ( double val: x ) {
     fprintf( FID, "%f,", val );
   }   
   fprintf(FID, "];\n");
   
   
   fprintf(FID, "u=[");
   for ( double val: u ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");
   
   
   fprintf(FID, "w=[");
   for ( double val: w ) {
     fprintf( FID, "%f,", val );
   }     
   fprintf(FID, "];\n");
   
   fprintf(FID, "phi=[");
   for ( double val: phi ) {
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
   

}
  

FILE *
NlBeam_SM :: giveOutputStream(TimeStep *tStep)
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
