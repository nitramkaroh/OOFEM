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
 *
 *               Copyright (C) 1993 - 2015   Borek Patzak
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


#include "../sm/Elements/Micromorphic/basesecondgradientelement.h"
#include "../sm/Materials/Micromorphic/secondgradientms.h"
#include "../sm/Materials/Micromorphic/secondgradientmaterialextensioninterface.h"


#include "../sm/CrossSections/structuralcrosssection.h"
#include "material.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "domain.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "unknownnumberingscheme.h"


#include <cstdio>

namespace oofem {
BaseSecondGradientElement :: BaseSecondGradientElement()
{
}

    

void
BaseSecondGradientElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &B)
{
  this->giveElement()->computeBHmatrixAt(gp, B);
}

  
 
void
BaseSecondGradientElement :: computeGeneralizedStressVectors(FloatArray &vP, FloatArray &vM, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    SecondGradientMaterialExtensionInterface *secondGradientMat = static_cast< SecondGradientMaterialExtensionInterface * >(cs->giveMaterialInterface(SecondGradientMaterialExtensionInterfaceType, gp) );
    if ( !secondGradientMat ) {
        OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
    }

    FloatArray dudx, d2udx2;
    this->computeDisplacementGradients(dudx,d2udx2, gp, tStep);
    secondGradientMat->giveGeneralizedStressVectors(vP, vM, dudx, d2udx2, gp, tStep);
    
}


void
BaseSecondGradientElement :: computeDisplacementGradients(FloatArray &dudx, FloatArray &d2udx2, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    FloatMatrix B, G;
    NLStructuralElement *elem = this->giveElement();
    elem->computeVectorOf({D_u, D_v,D_w}, VM_Total, tStep, u);
    elem->computeDeformationGradientVector(dudx, gp, tStep, VM_Total);
    this->computeGmatrixAt(gp, G);
    d2udx2.beProductOf(G, u);
}

 

void
BaseSecondGradientElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    answer.clear();
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatArray BS, GM, vP, vM, vP_full;
    FloatMatrix B, G;
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      if ( useUpdatedGpRecord == 1 ) {
        vP_full = static_cast< SecondGradientMaterialStatus * >( gp->giveMaterialStatus())->giveTempPVector();
	vM = static_cast< SecondGradientMaterialStatus * >( gp->giveMaterialStatus())->giveTempMVector();
	StructuralMaterial :: giveReducedVectorForm( vP, vP_full, gp->giveMaterialMode()); 
      } else {
	this->computeGeneralizedStressVectors(vP, vM, gp, tStep);
      }

    // Compute nodal internal forces at nodes as f = B^T*Stress dV
    double dV  = elem->computeVolumeAround(gp);
    this->computeBmatrixAt(gp, B);
    this->computeGmatrixAt(gp, G);
    //
    BS.beTProductOf(B, vP);
    GM.beTProductOf(G, vM);
    //
    answer.add(dV, BS);
    answer.add(dV, GM);
    }
}


void
BaseSecondGradientElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    NLStructuralElement *elem = this->giveElement();
    int ndofs = elem->computeNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatMatrix B,G;
    FloatMatrix A_FF, A_GF, A_FG, A_GG, AFF_B, AGF_B, AFG_G, AGG_G;
    answer.clear();
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      this->computeBmatrixAt(gp,B);
      this->computeGmatrixAt(gp,G);
      SecondGradientMaterialExtensionInterface *secondGradientMat = static_cast< SecondGradientMaterialExtensionInterface * >(cs->giveMaterialInterface(SecondGradientMaterialExtensionInterfaceType, gp) );
      if ( !secondGradientMat ) {
        OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
      }

      secondGradientMat->giveSecondGradientMatrix_dPdF(A_FF, rMode, gp, tStep);
      secondGradientMat->giveSecondGradientMatrix_dMdF(A_GF, rMode, gp, tStep);
      secondGradientMat->giveSecondGradientMatrix_dPdG(A_FG, rMode, gp, tStep);
      secondGradientMat->giveSecondGradientMatrix_dMdG(A_GG, rMode, gp, tStep);
      double dV = elem->computeVolumeAround(gp);
      AFF_B.beProductOf(A_FF, B);
      if(A_GF.giveNumberOfRows()) {
	AGF_B.beProductOf(A_GF, B);
      }
      if(A_FG.giveNumberOfRows()) {
	AFG_G.beProductOf(A_FG, G);
      }
      AGG_G.beProductOf(A_GG, G);

      answer.plusProductUnsym(B, AFF_B, dV);
      if(A_GF.giveNumberOfRows()) {
	answer.plusProductUnsym(G, AGF_B, dV);
      }
      if(A_FG.giveNumberOfRows()) {
	answer.plusProductUnsym(B, AFG_G, dV);
      }
      answer.plusProductUnsym(G, AGG_G, dV);
    }
}







 

IRResultType
BaseSecondGradientElement :: initializeFrom(InputRecord *ir)
{
  // @todo Is this function necessary???

    return IRRT_OK;
}

void
BaseSecondGradientElement :: updateInternalState(TimeStep *tStep)
// Updates the receiver at end of step.
{
    FloatArray stress, strain;
    /*
    // force updating strains & stresses
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            this->computeStrainVector(strain, gp, tStep);
            this->computeStressVector(stress, strain, gp, tStep);
        }
    }
    */
}







} // end namespace oofem


