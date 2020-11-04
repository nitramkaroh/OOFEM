/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2020   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "../sm/Elements/Beams/beam2dnl.h"
#include "../sm/Materials/structuralms.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "engngm.h"
#include "mathfem.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(Beam2dNl);


Beam2dNlnl :: Beam2dNl(int n, Domain *aDomain) : BeamBaseElement(n, aDomain)
{
    numberOfDofMans = 2;
    length = 0.;
    pitch = 10.;  // a dummy value
}


Beam2dNl :: ~Beam2dNl()
{
}




void
Beam2dNl :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    // compute clamped stifness
    this->computeLocalStiffnessMatrix(answer, rMode, tStep);
}



void
Beam2dNl :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->give2dBeamStiffMtrx(answer, rMode, gp, tStep);
}


void
Beam2dNl :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveGeneralizedStress_Beam2dNl(answer, gp, strain, tStep);
}


bool
Beam2dNl :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver.
{
    double sine, cosine;

    int ndofs = computeNumberOfGlobalDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    sine = sin( this->givePitch() );
    cosine  = cos(pitch);
    answer.at(1, 1) =  cosine;
    answer.at(1, 2) =  sine;
    answer.at(2, 1) = -sine;
    answer.at(2, 2) =  cosine;
    answer.at(3, 3) =  1.;
    answer.at(4, 4) =  cosine;
    answer.at(4, 5) =  sine;
    answer.at(5, 4) = -sine;
    answer.at(5, 5) =  cosine;
    answer.at(6, 6) =  1.;

    for ( int i = 7; i <= ndofs; i++ ) {
        answer.at(i, i) = 1.0;
    }

    return true;
}




void
Beam2dNl :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_w, R_v
    };
}


double
Beam2dNl :: computeLength()
// Returns the length of the receiver.
{
    double dx, dy;
    Node *nodeA, *nodeB;

    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        length  = sqrt(dx * dx + dy * dy);
    }

    return length;
}


double
Beam2dNl :: givePitch()
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



int
Beam2dNl :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    double sine, cosine;

    answer.resize(3, 3);
    answer.zero();

    sine = sin( this->givePitch() );
    cosine = cos(pitch);

    answer.at(1, 1) = cosine;
    answer.at(1, 2) = sine;
    answer.at(2, 1) = -sine;
    answer.at(2, 2) = cosine;
    answer.at(3, 3) = 1.0;

    return 1;
}


IRResultType
Beam2dNl :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    // first call parent
    BeamBaseElement :: initializeFrom(ir);
    return IRRT_OK;
}

double
Beam2dNl :: computeMomentFromCurvature(double kappa)
{
  if (material==1) { // linear moment-curvature
    return EI*kappa;
  } else {
    // moment-curvature based on arctan stress-strain law
    double M = 0.;
    if (kappa!=0.){
      double aux=kappa/kappam;
      M = 1.5*EI*kappam*((1./aux+aux)*atan(aux)-1.)/aux;
    }
    return M;
  }
}

double
Beam2dNl :: computeDerMomentFromCurvature(double kappa)
{
  if (material==1) { // linear moment-curvature
   return EI;
  } else {
    // moment-curvature based on arctan stress-strain law
    double dM = EI;
    if (kappa!=0.){
      double aux=kappa/kappam;
      dM = (3.*EI/(aux*aux))*(1.-atan(aux)/aux);
    }
    return dM;
  }
}

double computeCurvatureFromMoment(double M)
{
 double kappa_lin = M/EI;
 if (material==1) { // linear moment-curvature
   return kappa_lin;
 } else {
   // moment-curvature based on arctan stress-strain law
   int iter = 1;
   double kappa = kappa_lin;
   double res = computeMomentFromCurvature(kappa)-M;
   while (fabs(res)>TOL_SECTION*fabs(M) && iter++<=MAXIT_SECTION){
     double der = computeDerMomentFromCurvature(kappa);
     kappa -= res/der;
     res = computeMomentFromCurvature(kappa)-M;
   }
   if (iter>MAXIT_SECTION){
     printf("No convergence in computeCurvatureFromMoment for M = %g\n",M);
     exit(0);
   }	
   return kappa;
 }
}

void
Beam2dNl :: integrateAlongBeamAndGetJacobi(FloatArray fab, FloatArray ub, FloatMatrix jacobi)
{

  FloatArray x(NN+1), u(NN+1), w(NN+1), phi(NN+1);
  FloatMatrix du(3), dw(3), dphi(3), dM(3), dkappa(3), dphi_mid(3), dN_mid(3);
  double Xab = fab.at(1), Zab = fab.at(2), Mab = fab.at(3);
 
  double dx = beamLength/NN;
  
  for (int i = 2; i<= NN+1; i++){
    x.at(i) = x.at(i-1) + dx;
    double M = -Mab + Xab * w.at(i) - Zab * ( x.at(i) + u.at(i) );
    dM.at(1) = w.at(i-1);
    dM.at(2) = - ( x.at(i-1) + u.at(i-1) );
    dM.at(3) = -1.;
    double kappa = computeCurvatureFromMoment(M);
    double dMdkappa = computeDerMomentFromCurvature(kappa);
    dkappa = dM / dMdkappa;
    double phi_mid = phi.at(i) + kappa * dx/2.;
    dphi_mid = dphi + dkappa * dx / 2.;
    double N_mid = - Xab * cos(phi_mid) + Zab * sin(phi_mid);
    dN_mid = (Xab * sin(phi_mid) + Zab * cos(phi_mid)) * dphi_mid.at(j);
    dN_mid.at(1) -= cos(phi_mid);
    dN_mid.at(2) += sin(phi_mid);
    
    u.at(i) = u.at(i-1) + dx * ( (1. + N_mid / EA) * cos(phi_mid) -1. );
    du.at(j) += dx * (dN_mid[j] / EA) * cos(phi_mid) - dx * (1. + N_mid / EA) * sin(phi_mid) * dphi_mid.at(j);
    w.at(i) = w.at(i-1) - dx * ( 1. + N_mid/EA ) * sin(phi_mid);
    dw.at(j) -= dx * ( dN_mid[j] / EA ) * sin(phi_mid) + dx * ( 1. + N_mid/EA ) * cos(phi_mid) * dphi_mid.at(j);
    M = -Mab+Xab*w[i]-Zab*(x[i]+u[i]);
    dM[0] = w[i]; dM[1] = -(x[i]+u[i]); dM[2] = -1.;
    kappa = computeCurvatureFromMoment(M);
    dMdkappa = computeDerMomentFromCurvature(kappa);
    dkappa[j] = dM[j]/dMdkappa;
    phi[i] = phi_mid+kappa*dx/2.;
    dphi[j] = dphi_mid[j]+dkappa[j]*dx/2.;
  }
  ub[0] = u[NN];
  ub[1] = w[NN];
  ub[2] = phi[NN];
 
  jacobi[0][j]=du[j];
  jacobi[1][j]=dw[j];
  jacobi[2][j]=dphi[j];
}

/*
Find forces and moment at the left end that lead to given displacements and rotations
at the right end, assuming that the displacement and rotation at the left end are zero.
The solution is computed iteratively, using the Newton-Raphson method. 
Everything is done here in the local coordinate system aligned with the undeformed beam.
*/
bool findLeftEndForcesLocal(double ub_target[3], double fab_loc[3])
{
  double ub, wb, phib;
  FloatArray res(3), dforces(3), ub_loc(3);
  double jacobi[3][3];
  int iter = 0, i, j;
  double tolerance = TOL_BEAM*sqrt(ub_target.computeNorm());
  integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, jacobi);
  //
  residuum(ub_target);
  residiuum.subtract(ub_loc);
  double error = residuum.computeNorm();  
  while (error>tolerance && iter++ < MAXIT_BEAM) {    
    jacobi.solveForRhs(res, dforces);
    fab_loc.add(dforces);
    integrateAlongBeam(fab_loc, G, b_loc, );
    residuum = ub_target - ub_loc; 
    error = residuum.computeNorm();
    if (iter>MAXIT_BEAM) {
      this->giveDomain()->giveEngngModel()->setAnalysisCrash(true);
    }
  }
  
}


this->giveLocalInternalForcesVector(f_loc, u_loc, false);

  
/*
Auxiliary matrix for transformation of displacements (or in trasposed form of forces).
It corresponds to rotation from global axes to the auxiliary system aligned with the deformed beam.
*/
void construct_T(FloatMatrix &T, double phia)
{
  
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = cos(alpha-phia);
  T.at(3,3) = 1.;
  
  T.at(1,2) = sin(alpha-phia);
  T.at(2,1) = -T.at(1,2);
}

void construct_l(FloatArray &l, double phia)
{
  l.resize(3);
  l.at(1) = beamLength * (cos(phia)-1.);
  l.at(2) = sin(phia);
}

void
Beam2dNl :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  this->computeLocalDisplacementVector(u_loc, T, u_glob);
  this->giveLocalInternalForcesVector(f_loc, u_loc, false);
  answer.beProductOf(T, f_loc);

}



void
Beam2dNl :: printOutputAt(FILE *File, TimeStep *tStep)
{
    FloatArray rl, Fl;

    fprintf(File, "beam element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    // ask for global element displacement vector
    this->computeVectorOf(VM_Total, tStep, rl);

    // ask for global element end forces vector
    this->giveEndForcesVector(Fl, tStep);

    fprintf(File, "  local displacements ");
    for ( auto &val : rl ) {
        fprintf(File, " %.4e", val);
    }

    fprintf(File, "\n  local end forces    ");
    for ( auto &val : Fl ) {
        fprintf(File, " %.4e", val);
    }

}



} // end namespace oofem
