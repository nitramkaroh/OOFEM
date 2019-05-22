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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#include "../sm/Elements/PlaneStress/quadmembraneSE.h"
#include "../sm/Materials/structuralms.h"
#include "fei2dquadlin.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "strainvector.h"
#include "classfactory.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Elements/nlstructuralelement.h"
#include "../sm/Materials/structuralms.h"

namespace oofem {
REGISTER_Element(QuadMembraneSE);

FEI2dQuadLin QuadMembraneSE :: interpolation(1, 2);

QuadMembraneSE :: QuadMembraneSE(int n, Domain *aDomain) :
    PlaneStress2d(n, aDomain), PressureFollowerLoadElementInterface(this)
    // Constructor.
{
    numberOfDofMans  = 4;
    numberOfGaussPoints = 4;
    nlGeometry = 1;
}

QuadMembraneSE :: ~QuadMembraneSE()
// Destructor
{ }

FEInterpolation *QuadMembraneSE :: giveInterpolation() const { return & interpolation; }

void
QuadMembraneSE :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, int li, int ui)
//
// Returns the [3x8] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );

    answer.resize(3, 12);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(i, 1);
        answer.at(2, 3 * i - 1) = dnx.at(i, 2);
	answer.at(3, 3 * i - 2) = dnx.at(i, 2);
	answer.at(3, 3 * i - 1) = dnx.at(i, 1);
    }

 
}


void
QuadMembraneSE :: computeNlBmatrixAt(GaussPoint *gp, FloatMatrix &Bnl, FloatMatrix &G,TimeStep *tStep, int li, int ui)
//
// Returns the [3x8] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
  FloatArray vU, dU;
  FloatMatrix dnx, A;
  

    this->interpolation.evaldNdx( dnx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
    this->computeVectorOf({D_u, D_v, D_w}, VM_Total, tStep, vU); // solution vector    


    G.resize(6, 12);
    G.zero();

    for ( int i = 1; i <= 4; i++ ) {
        G.at(1, 3 * i - 2) = dnx.at(i, 1);
        G.at(2, 3 * i - 1) = dnx.at(i, 2);
	G.at(3, 3 * i - 2) = dnx.at(i, 2);

	G.at(4, 3 * i - 0) = dnx.at(i, 2);
        G.at(5, 3 * i - 0) = dnx.at(i, 1);
	G.at(6, 3 * i - 1) = dnx.at(i, 1);
    }

    dU.beProductOf(G,vU);
    A = {{dU.at(1), 0, dU.at(3)},{0, dU.at(2), dU.at(6)},{0, dU.at(3),dU.at(1)},{0, dU.at(4),dU.at(5)},{dU.at(5), 0, dU.at(4)},{dU.at(6), 0, dU.at(2)}};
    Bnl.beProductOf(A,G);

}


void
QuadMembraneSE :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, double alpha)
//
// Returns the [4x8] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// @todo not checked if correct
{
    FloatMatrix dNdx;

    this->interpolation.evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper(tStep, alpha));
    answer.resize(9, 12);

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);     // dv/dy
        answer.at(6, 3 * i - 2) = dNdx.at(i, 2);     // du/dy
	answer.at(7, 3 * i - 0) = dNdx.at(i, 2);     // dw/dy
        answer.at(8, 3 * i - 0) = dNdx.at(i, 1);     // dw/dx
        answer.at(9, 3 * i - 1) = dNdx.at(i, 1);     // dv/dx
    }

}



void
QuadMembraneSE :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
      D_u, D_v, D_w
    };
}

  

IRResultType
QuadMembraneSE :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 4;
    IRResultType result = PlaneStressElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    
    if ( numberOfGaussPoints != 1 && numberOfGaussPoints != 4 && numberOfGaussPoints != 9 && numberOfGaussPoints != 16 && numberOfGaussPoints != 25 ) {
        numberOfGaussPoints = 4;
        OOFEM_WARNING("Number of Gauss points enforced to 4");
    }

    return IRRT_OK;
}


void
QuadMembraneSE :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType modeType, double alpha)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    this->computeVectorOf({D_u, D_v, D_w}, modeType, tStep, u); // solution vector
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt(gp, B, tStep, alpha);
    answer.beProductOf(B, u);

    answer.at(1) += 1.0;
    answer.at(2) += 1.0;
    answer.at(3) += 1.0;
    
}


void
QuadMembraneSE :: computeStiffnessMatrix(FloatMatrix &answer,MatResponseMode rMode, TimeStep *tStep)
{

  answer.clear();


  double dV;
  FloatArray u, strain, S;
  FloatMatrix B, BE, Bnl, d, dbj, G;
  //testing
  FloatMatrix ism(6,6);
  ism.zero();
  //bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
  for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
    this->computeBmatrixAt(gp, B, tStep);
    this->computeNlBmatrixAt(gp, Bnl, G, tStep);
    BE = Bnl;
    BE.times(0.5);
    BE.add(B);
    B.add(Bnl);
    this->giveStructuralCrossSection()->giveStiffnessMatrix_PlaneStress(d, rMode, gp, tStep);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStructuralCrossSection()->giveMaterial(gp)->giveStatus(gp) );
    FloatArray redvS = status->giveTempStressVector();
    FloatArray vS;
    StructuralMaterial :: giveFullSymVectorForm(vS, redvS, _PlaneStress);
    dV = this->computeVolumeAround(gp);
    dbj.beProductOf(d, B);
    answer.plusProductUnsym(B, dbj, dV);    

    if(vS.at(1) == 0 && vS.at(2) == 0) {
	vS.at(1) = 1000000;
	vS.at(2) = 1000000;
    }

    FloatMatrix mS;
    mS.resize(6,6);
   
    mS.zero();
    mS = {{vS.at(1), 0, vS.at(6), 0, 0, 0},{ 0, vS.at(2),0, 0, 0,vS.at(6)},{vS.at(6), 0, vS.at(2), 0, 0, 0},{0, 0, 0, vS.at(2), vS.at(6), 0},{ 0, 0, 0, vS.at(6), vS.at(1), 0},{0, vS.at(6), 0, 0, 0, vS.at(1)}};

    FloatMatrix SG;
    SG.beProductOf(mS,G);
    answer.plusProductUnsym(G, SG, dV);    
    ism.plusProductUnsym(G,SG,dV);
    

  }


}

void
QuadMembraneSE :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  FloatMatrix B, BE, Bnl, G;
  FloatArray u, stress, strain;
  //@todo:test to remove
  FloatMatrix k1(12,12), k2(12,12),k3(12,12), k4(12,12), D, DBl, DBnl;
  k1.zero();
  k2.zero();
  k3.zero();
  k4.zero();
  /////////////////////////////////////////
  // This function can be quite costly to do inside the loops when one has many slave dofs.
  this->computeVectorOf(VM_Total, tStep, u);
  // subtract initial displacements, if defined
  if ( initialDisplacements ) {
    u.subtract(* initialDisplacements);
  }

  // zero answer will resize accordingly when adding first contribution
  answer.clear();

  for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
    this->computeBmatrixAt(gp, B, tStep);
    this->computeNlBmatrixAt(gp, Bnl, G, tStep);
    //@todo: test- matrix to remove
    FloatMatrix Bl;
    Bl = B;
    //
    BE = Bnl;
    BE.times(0.5);
    BE.add(B);
    B.add(Bnl);
    if ( !this->isActivated(tStep) ) {
      strain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
      strain.zero();
    }
    strain.beProductOf(BE, u);
    this->giveStructuralCrossSection()->giveRealStress_PlaneStress(stress, gp, strain, tStep);    

    // updates gp stress and strain record  acording to current
    // increment of displacement
    if ( stress.giveSize() == 0 ) {
      break;
    }
    // compute nodal representation of internal forces using f = B^T*Sigma dV
    double dV = this->computeVolumeAround(gp);
    answer.plusProduct(B, stress, dV);

    //@todo: test to remove
    /*   this->giveStructuralCrossSection()->giveStiffnessMatrix_PlaneStress(D, TangentStiffness, gp, tStep);
    DBl.beProductOf(D,Bl);
    DBnl.beProductOf(D,Bnl);
    k1.plusProductUnsym(Bl, DBl, dV);
    k2.plusProductUnsym(Bnl, DBl, dV);
    k3.plusProductUnsym(Bl, DBnl, 0.5*dV);
    k4.plusProductUnsym(Bnl, DBnl, 0.5*dV);      */
  }
  
  // if inactive update state, but no contribution to global system
  if ( !this->isActivated(tStep) ) {
    answer.zero();
    return;
  }
  
}

 
// support for pressure follower load interface
Interface*
QuadMembraneSE :: giveInterface(InterfaceType interface)
{
    if ( interface == PressureFollowerLoadElementInterfaceType) {
        return static_cast< PressureFollowerLoadElementInterface* >(this);
    }
    return NULL;
}


void
QuadMembraneSE :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    int numNodes = this->giveNumberOfDofManagers();
    FloatArray N(numNodes);

    int dim = 3;

    answer.resize(dim, dim * numNodes);
    answer.zero();
    giveInterpolation()->evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );

    answer.beNMatrixOf(N, dim);
}
void
QuadMembraneSE::computeSurfaceNMatrix (FloatMatrix &answer, int boundaryID, const FloatArray& lcoords)
{
  FloatArray n_vec;
  this->giveInterpolation()->boundarySurfaceEvalN(n_vec, boundaryID, lcoords, FEIElementGeometryWrapper(this) );
  answer.beNMatrixOf(n_vec, 3);
}
  
void
QuadMembraneSE ::  surfaceEvaldNdxi(FloatMatrix &answer, int iSurf, GaussPoint *gp)
{

  FEI2dQuadLin *interp2d = static_cast<FEI2dQuadLin*> (this->giveInterpolation());
  FloatMatrix dn;
  interp2d->giveDerivatives(answer, gp->giveNaturalCoordinates());
  /*  answer.resize(2, dn.giveNumberOfRows()*3);
  for ( int i = 1; i <= dn.giveNumberOfRows(); i++ ) {
    answer.at(1, i * 3 - 2) = dn.at(i, 1);
    answer.at(2, i * 3 - 1) = dn.at(i, 2);
    }*/

}




  
void
QuadMembraneSE ::  surfaceEvalDeformedNormalAt(FloatArray &answer, FloatArray &dxdeta, FloatArray &dxdxi, int iSurf, GaussPoint *gp, TimeStep *tStep)
{
  IntArray bNodes;
  FloatArray lcoords, vU;
  FloatMatrix dNdxi, dx, x;

  lcoords = gp->giveNaturalCoordinates();
  this->giveBoundarySurfaceNodes (bNodes, iSurf);

  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);
  double nNodes = bNodes.giveSize();
  x.resize(nNodes,3);
  if(this->domain->giveEngngModel()->giveFormulation() != AL) {
    // compute actual node positions for Total Lagrangean formulation
    this->computeBoundaryVectorOf(bNodes, {D_u, D_v, D_w}, VM_Total, tStep, vU); // solution vector    
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(i,1) = node->giveCoordinate(1) + vU.at( (i-1) * 3 + 1);
      x.at(i,2) = node->giveCoordinate(2) + vU.at( (i-1) * 3 + 2);
      x.at(i,3) = node->giveCoordinate(3) + vU.at( (i-1) * 3 + 3);
    }
  } else {
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(i,1) = node->giveCoordinate(1);
      x.at(i,2) = node->giveCoordinate(2);
      x.at(i,3) = node->giveCoordinate(3);
    }
    
  }

  
  dx.beTProductOf(dNdxi,x);
  dx.copyRow(dxdeta, 1);
  dx.copyRow(dxdxi, 2);
  answer.beVectorProductOf(dxdeta, dxdxi);

  //test
  double J;
  FloatArray redvF, vF, N, n, answer2;
  FloatMatrix F, invF;
  N = {0,0,1};
  this->computeDeformationGradientVector(vF, gp, tStep,VM_Total);
  //  StructuralMaterial :: giveFullVectorFormF(vF, redvF, _Membrane2d);
  F.beMatrixForm(vF);
  J = F.giveDeterminant();
  invF.beInverseOf(F);
  answer2.beTProductOf(invF,N);
  answer2.times(J);
  //  answer = answer2;
  
  FloatArray gcoords;
  this->interpolation.local2global(gcoords, lcoords, FEIElementGeometryWrapper(this));
  double r = sqrt(gcoords.at(1)*gcoords.at(1) + gcoords.at(2)*gcoords.at(2));
  FloatArray dU, U, grad;
  FloatArray answer3(3);
  FloatArray answer4(3);
  FloatMatrix Bh, Nm;
  this->computeNmatrixAt(lcoords, Nm);
  computeBHmatrixAt(gp, Bh, tStep);
  dU.beProductOf(Bh, vU);
  U.beProductOf(Nm, vU);


  double ur = sqrt(U.at(1)*U.at(1) + U.at(2)*U.at(2));

  double cos = gcoords.at(1)/r;
  double sin = gcoords.at(2)/r;
  
  double du = dU.at(1) * cos + dU.at(6) * sin;
  double dw = dU.at(8) * cos + dU.at(7) * sin;
  


  
  
  answer3.at(1) = -(r+ur)/r*dw * cos;
  answer3.at(2) = -(r+ur)/r*dw * sin;
  answer3.at(3) = (r+ur)/r*(1+du);


  answer4.at(1) = -dw * cos;
  answer4.at(2) = -dw * sin;
  answer4.at(3) = (1+du);

  
  this->computeDeformationGradientVector(vF, gp, tStep, VM_Incremental, 1);
  F.beMatrixForm(vF);
  invF.beInverseOf(F);
  J = F.giveDeterminant();
  answer4.beTProductOf(invF,N);
  answer4.times(J);


  
  

  

  FloatArray dudksi, dudeta, a1, a2;
  FloatMatrix u, dudxi;
  u.resize(nNodes,3);
  for(int i = 1; i <= nNodes; i++) {
    u.at(i,1) = vU.at( (i-1) * 3 + 1);
    u.at(i,2) = vU.at( (i-1) * 3 + 2);
    u.at(i,3) = vU.at( (i-1) * 3 + 3);
  }
  

  dudxi.beTProductOf(dNdxi,u);
  dudxi.copyRow(dudksi, 2);
  dudxi.copyRow(dudeta, 1);

  a1.beVectorProductOf(dudksi, dxdeta);
  a2.beVectorProductOf(dudeta, dxdxi);

  /*  answer.add(a1);
  answer.subtract(a2);
  */
  /*

  FloatArray tu, tF, strain, nlstrain, linstrain;
  FloatMatrix tB, tBh, tBnl, tBE, tE, tG;
  this->computeVectorOf(VM_Total, tStep, tu);
  // subtract initial displacements, if defined
  
  for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
    this->computeBmatrixAt(gp, tB, tStep);
    this->computeNlBmatrixAt(gp, tBnl, tG, tStep);
    //@todo: test- matrix to remove
    FloatMatrix Bl, FF;
    this->computeDeformationGradientVector(tF, gp, tStep, VM_Total);
    this->computeBHmatrixAt(gp, tBh, tStep, 0);
    FF.beMatrixForm(tF);
   
    tE.beTProductOf(FF,FF);
    tE.at(1,1) -= 1;
    tE.at(2,2) -= 1;
    tE.at(3,3) -= 1;
    tE.times(0.5);

    FloatArray grad;
    FloatMatrix A, AA;
    grad.beProductOf(tG, tu);
    A = {{grad.at(1), 0, grad.at(3)},{0, grad.at(2), grad.at(6)},{0, grad.at(3),dU.at(1)},{0, grad.at(4),grad.at(5)},{grad.at(5), 0, grad.at(4)},{grad.at(6), 0, grad.at(2)}};
    AA.beProductOf(A,grad);
    nlstrain.beProductOf(tBnl, tu);
    nlstrain.times(0.5);
    //    tE.subtract(nlStrain);

    linstrain.beProductOf(tB, tu);
    
    //
    tBE = tBnl;
    tBE.times(0.5);
    tBE.add(tB);
    tB.add(tBnl);
    strain.beProductOf(tBE, tu);
  }
  */


  
  
  
}
  


void
QuadMembraneSE ::  surfaceEvalNumericalStiffMatrixAt(FloatMatrix &answer, FloatArray &dxdeta, FloatArray &dxdxi, int iSurf, GaussPoint *gp, TimeStep *tStep)
{
  IntArray bNodes;
  FloatArray lcoords, vU;
  FloatMatrix dNdxi, dx, x;

  lcoords = gp->giveNaturalCoordinates();
  this->giveBoundarySurfaceNodes (bNodes, iSurf);

  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);
  double nNodes = bNodes.giveSize();
  x.resize(nNodes,3);
  if(this->domain->giveEngngModel()->giveFormulation() != AL) {
    // compute actual node positions for Total Lagrangean formulation
    this->computeBoundaryVectorOf(bNodes, {D_u, D_v, D_w}, VM_Total, tStep, vU); // solution vector    
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(i,1) = node->giveCoordinate(1) + vU.at( (i-1) * 3 + 1);
      x.at(i,2) = node->giveCoordinate(2) + vU.at( (i-1) * 3 + 2);
      x.at(i,3) = node->giveCoordinate(3) + vU.at( (i-1) * 3 + 3);
    }
  } else {
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(i,1) = node->giveCoordinate(1);
      x.at(i,2) = node->giveCoordinate(2);
      x.at(i,3) = node->giveCoordinate(3);
    }
    
  }

  FloatArray n;
  dx.beTProductOf(dNdxi,x);
  dx.copyRow(dxdeta, 1);
  dx.copyRow(dxdxi, 2);
  n.beVectorProductOf(dxdeta, dxdxi);
  answer.resize(12,12);
  double pert = 1.e-9;
  FloatArray v;
  FloatMatrix xp(x), N, Ki;
  this->computeNmatrixAt(lcoords, N);
  int index = 0;
  for (int i = 1; i<=4; i++) {
    for (int j = 1; j <= 3; j++) {
      index++;
      xp.at(i,j) += pert;      
      dx.beTProductOf(dNdxi,xp);
      dx.copyRow(dxdeta, 1);
      dx.copyRow(dxdxi, 2);
      v.beVectorProductOf(dxdeta, dxdxi);
      v.subtract(n);
      v.times(1/pert);
      Ki.beTProductOf(N, v);
      for (int k = 1; k <= 12; k++) {
	answer.at(k,index) = Ki.at(k,1);
      }
      xp = x;
    }
  }
  
}
  

} // end namespace oofem
