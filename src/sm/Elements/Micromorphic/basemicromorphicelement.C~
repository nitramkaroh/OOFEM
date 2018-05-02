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


#include "../sm/Elements/Micromorphic/basemicromorphicelement.h"
#include "../sm/Materials/Micromorphic/micromorphicms.h"
#include "../sm/Materials/Micromorphic/micromorphicmaterialextensioninterface.h"


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
BaseMicromorphicElement :: BaseMicromorphicElement()
{
}

    

void
BaseMicromorphicElement :: giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_m, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u,const IntArray &dofIdArray_m )
{
    // Routine to extract the location array of an element for given dofid array.
    locationArray_u.clear();
    locationArray_m.clear();
    NLStructuralElement *el = this->giveElement();
    int k = 0;
    IntArray nodalArray;
    for(int i = 1; i <= el->giveNumberOfDofManagers(); i++) {
      DofManager *dMan = el->giveDofManager( i );
      int itt = 1;
      for(int j = 1; j <= dofIdArray_u.giveSize( ); j++) {
	if(dMan->hasDofID( (DofIDItem) dofIdArray_u.at( j ) )) {
	  //  Dof *d = dMan->giveDofWithID( dofIdArray_u.at( j ) );
	  locationArray_u.followedBy( k + itt);
	}
	itt++;
      }
      for(int j = 1; j <= dofIdArray_m.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_m.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_m.followedBy( k + itt);
	}
	itt++;
      }
      k += dMan->giveNumberOfDofs( );
    }
}

void
BaseMicromorphicElement :: computeB_uMatrixAt(GaussPoint *gp, FloatMatrix &B, NLStructuralElement *element, bool isStressTensorSymmetric)
{
    if(isStressTensorSymmetric && element->giveGeometryMode() == 0) {
      element->computeBmatrixAt(gp, B);
    } else {
      element->computeBHmatrixAt(gp, B);
    }
}
 
void
BaseMicromorphicElement :: computeGeneralizedStressVectors(FloatArray &stress, FloatArray &s, FloatArray &M, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    MicromorphicMaterialExtensionInterface *micromorphMat = static_cast< MicromorphicMaterialExtensionInterface * >(cs->giveMaterialInterface(MicromorphicMaterialExtensionInterfaceType, gp) );
    if ( !micromorphMat ) {
        OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
    }

    IntArray IdMask_m;
    FloatArray micromorphicVar, micromorphicVarGrad;
    
    if(elem->giveGeometryMode() == 0) {
      FloatArray displacementGradient;
      this->computeDisplacementGradient(displacementGradient, gp, tStep, micromorphMat->isStressTensorSymmetric());
    this->giveDofManDofIDMask_m( IdMask_m );
    this->computeMicromorphicVars(micromorphicVar,micromorphicVarGrad,IdMask_m, gp, tStep);   
    micromorphMat->giveGeneralizedStressVectors(stress, s, M, gp, displacementGradient, micromorphicVar, micromorphicVarGrad, tStep);
    } else {
      FloatArray deformationGradient;
      elem->computeDeformationGradientVector(deformationGradient, gp, tStep);
      this->giveDofManDofIDMask_m( IdMask_m );
      this->computeMicromorphicVars(micromorphicVar,micromorphicVarGrad,IdMask_m, gp, tStep);   
      micromorphMat->giveFiniteStrainGeneralizedStressVectors(stress, s, M, gp, deformationGradient, micromorphicVar, micromorphicVarGrad, tStep);
    }
    
}


void
BaseMicromorphicElement :: computeDisplacementGradient(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, bool isStressTensorSymmetric)
{
    FloatArray u;
    FloatMatrix b;
    NLStructuralElement *elem = this->giveElement();
    elem->computeVectorOf({D_u, D_v,D_w}, VM_Total, tStep, u);
    this->computeB_uMatrixAt(gp, b, elem, isStressTensorSymmetric);
    answer.beProductOf(b, u);
}


void
BaseMicromorphicElement :: computeMicromorphicVars(FloatArray &micromorphicVar,FloatArray &micromorphicVarGrad, IntArray IdMask_m, GaussPoint *gp, TimeStep *tStep)
{
  FloatMatrix N_m, B_m;
  FloatArray u_m;
 
    this->computeMicromorphicNMatrixAt(gp, N_m);
    this->computeMicromorphicBMatrixAt(gp, B_m);
    /// @todo generalization for general micromorphic continua -- should be parameter of this function
    this->giveElement()->computeVectorOf(IdMask_m, VM_Total, tStep, u_m);
    micromorphicVar.beProductOf(N_m, u_m);
    micromorphicVarGrad.beProductOf(B_m, u_m);
}



void
BaseMicromorphicElement :: giveMicromorphicInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveElement();
    FloatMatrix N_m, B_m;
    answer.resize(this->giveNumberOfMicromorphicDofs());
    FloatArray vStress, vMicromorphicStress, vMicromorphicStressGrad;



    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      if ( useUpdatedGpRecord == 1 ) {
        vMicromorphicStress = static_cast< MicromorphicMaterialStatus * >( gp->giveMaterialStatus() )->giveTempMicromorphicStress();
        vMicromorphicStressGrad = static_cast< MicromorphicMaterialStatus* >( gp->giveMaterialStatus() )->giveTempMicromorphicStressGrad();
                               
      } else {
	this->computeGeneralizedStressVectors(vStress, vMicromorphicStress, vMicromorphicStressGrad, gp, tStep);      
      }

        // Compute nodal internal forces at nodes as f = \int (B^T*S-N^T*s) dV
      this->computeMicromorphicNMatrixAt(gp, N_m);
      this->computeMicromorphicBMatrixAt(gp, B_m);
      FloatArray NmSm, BmSmg;
      double dV  = elem->computeVolumeAround(gp);
      NmSm.beTProductOf(N_m, vMicromorphicStress);
      answer.add(dV, NmSm);
      BmSmg.beTProductOf(B_m, vMicromorphicStressGrad);
      answer.add(dV,BmSmg);
    }
}


void
BaseMicromorphicElement :: giveStandardInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatArray BS, vStress, s, M;
    FloatMatrix B;
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      if ( useUpdatedGpRecord == 1 ) {
	if(this->giveElement()->giveGeometryMode() == 0) {
	  vStress = static_cast< MicromorphicMaterialStatus * >( gp->giveMaterialStatus() )->giveTempStressVector();
	} else {
	  vStress = static_cast< MicromorphicMaterialStatus * >( gp->giveMaterialStatus() )->giveTempPVector();
	}	
      } else {
	this->computeGeneralizedStressVectors(vStress, s, M, gp, tStep);      
      }
    
      MicromorphicMaterialExtensionInterface *micromorphicMat = dynamic_cast< MicromorphicMaterialExtensionInterface * >(cs->giveMaterialInterface(MicromorphicMaterialExtensionInterfaceType, gp) );
      if ( !micromorphicMat ) {
	OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
      }
      // Compute nodal internal forces at nodes as f = B^T*Stress dV
      double dV  = elem->computeVolumeAround(gp);
      this->computeB_uMatrixAt(gp, B, elem, micromorphicMat->isStressTensorSymmetric());    
      BS.beTProductOf(B, vStress);
      answer.add(dV, BS);
    }
}

void
BaseMicromorphicElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    answer.resize(this->giveNumberOfDofs());
    answer.zero();
    FloatArray answerU(this->giveNumberOfDisplacementDofs());
    answer.zero();
    FloatArray answerM(this->giveNumberOfMicromorphicDofs());
    answerM.zero();
   
    this->giveMicromorphicInternalForcesVector(answerM, tStep, 0);
    this->giveStandardInternalForcesVector(answerU, tStep, 1);

    
   
    answer.assemble(answerU, locationArray_u);
    answer.assemble(answerM, locationArray_m);

    //answer.assemble(answerU, displacementDofsOrdering);
    //answer.assemble(answerM, micromorphicDofsOrdering);
}

void
BaseMicromorphicElement :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
 
    FloatArray localForces(this->giveNumberOfDisplacementDofs());
    FloatArray nlForces(this->giveNumberOfMicromorphicDofs());
    answer.resize(this->giveNumberOfDofs());

    this->computeLocForceLoadVector(localForces, tStep, mode);

    // answer.assemble(localForces, displacementDofsOrdering);
    //answer.assemble(nlForces, micromorphicDofsOrdering);
    answer.assemble(localForces, locationArray_u);
    answer.assemble(nlForces, locationArray_m);
}


/************************************************************************/
void
BaseMicromorphicElement :: computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    FloatMatrix T;
    NLStructuralElement *elem = this->giveElement();
    elem->computeLocalForceLoadVector(answer, tStep, mode);

    // transform result from global cs to nodal cs. if necessary
    if ( answer.isNotEmpty() ) {
        if ( elem->computeGtoLRotationMatrix(T) ) {
            // first back to global cs from element local
            answer.rotatedWith(T, 't');
        }
    } else {
        answer.resize(this->giveNumberOfDisplacementDofs());
        answer.zero();
    }
}


void
BaseMicromorphicElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
  answer.resize(this->giveNumberOfDofs(), this->giveNumberOfDofs());
  answer.zero();

    FloatMatrix Kuu, Kum, Kmu, Kmm;
    this->computeStiffnessMatrix_uu(Kuu, rMode, tStep);
    this->computeStiffnessMatrix_um(Kum, rMode, tStep);
    if(symmetricFormulation) {
      Kmu.beTranspositionOf(Kum);
    } else {
      this->computeStiffnessMatrix_mu(Kmu, rMode, tStep);
    }
    
    this->computeStiffnessMatrix_mm(Kmm, rMode, tStep);
 
    answer.assemble(Kuu, locationArray_u);
    answer.assemble(Kum, locationArray_u, locationArray_m);
    answer.assemble(Kmu, locationArray_m, locationArray_u);
    answer.assemble(Kmm, locationArray_m);


    //answer.assemble(answer1, displacementDofsOrdering);
    //answer.assemble(answer2, displacementDofsOrdering, micromorphicDofsOrdering);
    //answer.assemble(answer3, micromorphicDofsOrdering, displacementDofsOrdering);
    //answer.assemble(answer4, micromorphicDofsOrdering);
}


void
BaseMicromorphicElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatMatrix B, D, DB;
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.clear();
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {

        MicromorphicMaterialExtensionInterface *micromorphicMat = dynamic_cast< MicromorphicMaterialExtensionInterface * >(
            cs->giveMaterialInterface(MicromorphicMaterialExtensionInterfaceType, gp) );
        if ( !micromorphicMat ) {
            OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }
	bool isStressTensorSymmetric = micromorphicMat->isStressTensorSymmetric();
	this->computeB_uMatrixAt(gp,B, elem, isStressTensorSymmetric);

        micromorphicMat->giveMicromorphicMatrix_dSigdUgrad(D, rMode, gp, tStep);
        double dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D, B);

        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(B, DB, dV);
        } else {
            answer.plusProductUnsym(B, DB, dV);
        }


    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}



void
BaseMicromorphicElement :: computeStiffnessMatrix_um(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    double dV;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatMatrix B, N_m, D, DN_m;

    answer.clear();

    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
        MicromorphicMaterialExtensionInterface *micromorphicMat = dynamic_cast< MicromorphicMaterialExtensionInterface * >(cs->giveMaterialInterface(MicromorphicMaterialExtensionInterfaceType, gp) );
        if ( !micromorphicMat ) {
            OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }
	bool isStressTensorSymmetric = micromorphicMat->isStressTensorSymmetric();

	this->computeB_uMatrixAt(gp,B, elem, isStressTensorSymmetric);


        micromorphicMat->giveMicromorphicMatrix_dSigdPhi(D, rMode, gp, tStep);
        this->computeMicromorphicNMatrixAt(gp, N_m);

	dV = elem->computeVolumeAround(gp);
        DN_m.beProductOf(D, N_m);
        answer.plusProductUnsym(B, DN_m, dV);
    }

}




void
BaseMicromorphicElement :: computeStiffnessMatrix_mu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    double dV;
    NLStructuralElement *elem = this->giveElement();
    FloatMatrix B, DB, D, N_m;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();

    answer.clear();


    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
        MicromorphicMaterialExtensionInterface *micromorphicMat = dynamic_cast< MicromorphicMaterialExtensionInterface * >(cs->giveMaterialInterface(MicromorphicMaterialExtensionInterfaceType, gp) );
        if ( !micromorphicMat ) {
	    OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }
	bool isStressTensorSymmetric = micromorphicMat->isStressTensorSymmetric();
	this->computeB_uMatrixAt(gp, B, elem, isStressTensorSymmetric);

        micromorphicMat->giveMicromorphicMatrix_dSdUgrad(D, rMode, gp, tStep);
        this->computeMicromorphicNMatrixAt(gp, N_m);

	dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D,B);

        answer.plusProductUnsym(N_m,DB, dV);
    }

}


void
BaseMicromorphicElement :: computeStiffnessMatrix_mm(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    double dV;
    FloatMatrix lStiff;
    FloatMatrix B_m, N_m, dSdPhi, dMdPhiGrad, dSdPhi_N, dMdPhiGrad_B;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    answer.clear();

    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
        MicromorphicMaterialExtensionInterface *micromorphicMat = dynamic_cast< MicromorphicMaterialExtensionInterface * >(cs->giveMaterialInterface(MicromorphicMaterialExtensionInterfaceType, gp) );
        if ( !micromorphicMat ) {
            OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }

        micromorphicMat->giveMicromorphicMatrix_dSdPhi(dSdPhi, rMode, gp, tStep);
        micromorphicMat->giveMicromorphicMatrix_dMdPhiGrad(dMdPhiGrad, rMode, gp, tStep);
        this->computeMicromorphicNMatrixAt(gp, N_m);
        this->computeMicromorphicBMatrixAt(gp, B_m);
        dV = elem->computeVolumeAround(gp);

	dSdPhi_N.beProductOf(dSdPhi, N_m);
	answer.plusProductUnsym(N_m, dSdPhi_N, dV);	

	dMdPhiGrad_B.beTProductOf(dMdPhiGrad, B_m);
	answer.plusProductUnsym(B_m,dMdPhiGrad_B,dV);

    }
}

void
BaseMicromorphicElement :: computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    //set displacement and nonlocal location array
  answer.resize(this->giveNumberOfDofs(), this->giveNumberOfDofs());
  answer.zero();

    FloatMatrix Muu, Mmm;
    this->computeMassMatrix_uu(Muu, tStep);
    this->computeMassMatrix_mm(Mmm, tStep);
 
    answer.assemble(Muu, locationArray_u);
    answer.assemble(Mmm, locationArray_m);

}


void
BaseMicromorphicElement :: computeMassMatrix_uu(FloatMatrix &answer, TimeStep *tStep)
{
 
  double dV;
  FloatMatrix n;
  NLStructuralElement *elem = this->giveElement();
  
  answer.clear();

  for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
    elem->computeNmatrixAt(gp->giveSubPatchCoordinates(), n);
    dV = elem->computeVolumeAround(gp);
    answer.plusProductSymmUpper(n, n, dV);
  }
  
  answer.symmetrized();
}


void
BaseMicromorphicElement :: computeMassMatrix_mm(FloatMatrix &answer, TimeStep *tStep)
{
  
  double dV;
  FloatMatrix n;
  NLStructuralElement *elem = this->giveElement();
 
  answer.clear();

  for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
    this->computeMicromorphicNMatrixAt(gp, n);
    dV = elem->computeVolumeAround(gp);
    answer.plusProductSymmUpper(n, n, dV);
  }
  
  answer.symmetrized();
}

 

IRResultType
BaseMicromorphicElement :: initializeFrom(InputRecord *ir)
{
  // @todo Is this function necessary???

    return IRRT_OK;
}

void
BaseMicromorphicElement :: updateInternalState(TimeStep *tStep)
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


void
BaseMicromorphicElement :: postInitialize()
{
  IntArray IdMask_u, IdMask_m;
  this->giveDofManDofIDMask_u( IdMask_u );
  this->giveDofManDofIDMask_m( IdMask_m );
  this->giveLocationArrayOfDofIDs(locationArray_u,locationArray_m, EModelDefaultEquationNumbering(), IdMask_u, IdMask_m);
  //  this->giveLocationArrayOfDofIDs(locationArray_m,EModelDefaultEquationNumbering(), IdMask_m);
}




} // end namespace oofem


