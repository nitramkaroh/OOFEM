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

#include "Elements/3D/lspaceSE.h"
#include "fei3dhexalin.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "crosssection.h"
#include "classfactory.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Elements/nlstructuralelement.h"
#include "../sm/Materials/structuralms.h"

namespace oofem {
REGISTER_Element(LSpaceSE);

FEI3dHexaLin LSpaceSE :: interpolation;

  LSpaceSE :: LSpaceSE(int n, Domain *aDomain) : LSpace(n, aDomain)
    // Constructor.
{

}

void
LSpaceSE :: computeNlBmatrixAt(GaussPoint *gp, FloatMatrix &Bnl, FloatMatrix &G,TimeStep *tStep, int li, int ui)
//
// Returns the [3x8] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
  FloatArray vU, dU;
  FloatMatrix dNdx, A;
  

    this->interpolation.evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
    this->computeVectorOf({D_u, D_v, D_w}, VM_Total, tStep, vU); // solution vector    


    G.resize(9,24);
    G.zero();

    for ( int i = 1; i <= 8; i++ ) {
	G.at(1, 3 * i - 2) = dNdx.at(i, 1);     // du/dx
        G.at(2, 3 * i - 1) = dNdx.at(i, 2);     // dv/dy
        G.at(3, 3 * i - 0) = dNdx.at(i, 3);     // dw/dz
        G.at(4, 3 * i - 1) = dNdx.at(i, 3);     // dv/dz
        G.at(5, 3 * i - 2) = dNdx.at(i, 3);     // du/dz
	G.at(6, 3 * i - 2) = dNdx.at(i, 2);     // du/dy
        G.at(7, 3 * i - 0) = dNdx.at(i, 2);     // dw/dy
        G.at(8, 3 * i - 0) = dNdx.at(i, 1);     // dw/dx
        G.at(9, 3 * i - 1) = dNdx.at(i, 1);     // dv/dx	
    }

    dU.beProductOf(G,vU);
    A = {{dU.at(1), 0, 0, 0, dU.at(5), dU.at(6)}, {0, dU.at(2), 0, dU.at(4), 0, dU.at(9)}, {0, 0, dU.at(3), dU.at(7), dU.at(8), 0}, {0, 0, dU.at(4), dU.at(2), dU.at(9), 0},{0, 0, dU.at(5), dU.at(6), dU.at(1), 0},{0, dU.at(6), 0, dU.at(5), 0, dU.at(1)},{0, dU.at(7), 0, dU.at(2), 0, dU.at(8)},{dU.at(8), 0, 0, 0, dU.at(3), dU.at(7)},{dU.at(9), 0, 0, 0, dU.at(4), dU.at(2)}};
    Bnl.beProductOf(A,G);

}


void
LSpaceSE :: computeStiffnessMatrix(FloatMatrix &answer,MatResponseMode rMode, TimeStep *tStep)
{

  answer.clear();


  double dV;
  FloatArray u, strain, S;
  FloatMatrix B, BE, Bnl, d, dbj, G;
  //  bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
  for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
    this->computeBmatrixAt(gp, B, tStep);
    this->computeNlBmatrixAt(gp, Bnl, G, tStep);
    BE = Bnl;
    BE.times(0.5);
    BE.add(B);
    B.add(Bnl);
    this->giveStructuralCrossSection()->giveStiffnessMatrix_3d(d, rMode, gp, tStep);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStructuralCrossSection()->giveMaterial(gp)->giveStatus(gp) );
    FloatArray S = status->giveTempStressVector();
    dV = this->computeVolumeAround(gp);
    dbj.beProductOf(d, B);
    answer.plusProductUnsym(B, dbj, dV);    

    FloatMatrix mS;
    mS.resize(9,9);
    mS.zero();
    //    mS = {{S.at(1), 0, S.at(3), 0, 0, 0},{ 0, S.at(2),0, 0, 0,S.at(3)},{S.at(3), 0, S.at(2), 0, 0, 0},{0, 0, 0, S.at(2), S.at(3), 0},{ 0, 0, 0, S.at(3), S.at(1), 0},{0, S.at(3), 0, 0, 0, S.at(1)}};

    mS.at(1, 1) = S.at(1);
    mS.at(1, 5) = S.at(5);
    mS.at(1, 6) = S.at(6);
    
    mS.at(2, 2) = S.at(2);
    mS.at(2, 4) = S.at(4);
    mS.at(2, 9) = S.at(6);
    
    mS.at(3, 3) = S.at(3);
    mS.at(3, 7) = S.at(4);
    mS.at(3, 8) = S.at(5);
    
    mS.at(4, 2) = S.at(4);
    mS.at(4, 4) = S.at(3);
    mS.at(4, 9) = S.at(5);
    
    mS.at(5, 1) = S.at(5);
    mS.at(5, 5) = S.at(3);
    mS.at(5, 6) = S.at(4);
    
    mS.at(6, 1) = S.at(6);
    mS.at(6, 5) = S.at(4);
    mS.at(6, 6) = S.at(2);
    
    mS.at(7, 3) = S.at(4);
    mS.at(7, 7) = S.at(2);
    mS.at(7, 8) = S.at(6);
    
    mS.at(8, 3) = S.at(5);
    mS.at(8, 7) = S.at(6);
    mS.at(8, 8) = S.at(1);
    
    mS.at(9, 2) = S.at(6);
    mS.at(9, 4) = S.at(5);
    mS.at(9, 9) = S.at(1);
    
    
    
    
    FloatMatrix SG;
    SG.beProductOf(mS,G);
    answer.plusProductUnsym(G, SG, dV);    
        

  }


}

void
LSpaceSE :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
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

    this->giveStructuralCrossSection()->giveRealStress_3d(stress, gp, strain, tStep);    

    // updates gp stress and strain record  acording to current
    // increment of displacement
    if ( stress.giveSize() == 0 ) {
      break;
    }
    // compute nodal representation of internal forces using f = B^T*Sigma dV
    double dV = this->computeVolumeAround(gp);
    answer.plusProduct(B, stress, dV);

    //@todo: test to remove
    this->giveStructuralCrossSection()->giveStiffnessMatrix_3d(D, TangentStiffness, gp, tStep);
    DBl.beProductOf(D,Bl);
    DBnl.beProductOf(D,Bnl);
    k1.plusProductUnsym(Bl, DBl, dV);
    k2.plusProductUnsym(Bnl, DBl, dV);
    k3.plusProductUnsym(Bl, DBnl, 0.5*dV);
    k4.plusProductUnsym(Bnl, DBnl, 0.5*dV);      
  }
  
  // if inactive update state, but no contribution to global system
  if ( !this->isActivated(tStep) ) {
    answer.zero();
    return;
  }
  
}					       

} // end namespace oofem
