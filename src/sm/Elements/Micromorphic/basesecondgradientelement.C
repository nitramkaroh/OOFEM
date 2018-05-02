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
BaseSecondGradientElement :: giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_m, IntArray &locationArray_l, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u,const IntArray &dofIdArray_m, const IntArray &dofIdArray_l )
{
    // Routine to extract the location array of an element for given dofid array.
    locationArray_u.clear();
    locationArray_m.clear();
    locationArray_l.clear();
    NLStructuralElement *el = this->giveElement();
    int k = 0;
    IntArray nodalArray;
    for(int i = 1; i <= el->giveNumberOfDofManagers(); i++) {
      DofManager *dMan = el->giveDofManager( i );
      int itt = 1;
      for(int j = 1; j <= dofIdArray_u.giveSize( ); j++) {
	if(dMan->hasDofID( (DofIDItem) dofIdArray_u.at( j ) )) {
	  locationArray_u.followedBy( k + itt);
	}
	itt++;
      }
      for(int j = 1; j <= dofIdArray_m.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_m.at( j ) )) {
	  locationArray_m.followedBy( k + itt);
	}
	itt++;
      }
      for(int j = 1; j <= dofIdArray_l.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_l.at( j ) )) {
	  locationArray_l.followedBy( k + itt);
	}
	itt++;
      }
      k += dMan->giveNumberOfDofs( );
    }
}

void
BaseSecondGradientElement :: computeB_uMatrixAt(GaussPoint *gp, FloatMatrix &B, NLStructuralElement *element, bool isStressTensorSymmetric)
{
    if(isStressTensorSymmetric)
      element->computeBmatrixAt(gp, B);
    else
      element->computeBHmatrixAt(gp, B);
}
 
void
BaseSecondGradientElement :: computeGeneralizedStressVectors(FloatArray &sigma, FloatArray &s, FloatArray &M, FloatArray &relativeStrain, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    SecondGradientMaterialExtensionInterface *secondGradientMat = static_cast< SecondGradientMaterialExtensionInterface * >(cs->giveMaterialInterface(SecondGradientMaterialExtensionInterfaceType, gp) );
    if ( !secondGradientMat ) {
        OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
    }

    IntArray IdMask_m, IdMask_l;
    FloatArray displacementGradient, micromorphicVar, micromorphicVarGrad, lagrangianMultipliers;
    
    this->computeDisplacementGradient(displacementGradient, gp, tStep, secondGradientMat->isStressTensorSymmetric());

    this->giveDofManDofIDMask_m( IdMask_m );
    this->computeMicromorphicVars(micromorphicVar,micromorphicVarGrad,IdMask_m, gp, tStep);   

    this->giveDofManDofIDMask_m( IdMask_m );
    this->computeLagrangianMultipliers(lagrangianMultipliers,IdMask_l, gp, tStep);   

    secondGradientMat->giveGeneralizedStressVectors(sigma, s, M, relativeStrain, gp, displacementGradient, micromorphicVar, micromorphicVarGrad, tStep);
    
}


void
BaseSecondGradientElement :: computeDisplacementGradient(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, bool isStressTensorSymmetric)
{
    FloatArray u;
    FloatMatrix b;
    NLStructuralElement *elem = this->giveElement();
    elem->computeVectorOf({D_u, D_v,D_w}, VM_Total, tStep, u);
    this->computeB_uMatrixAt(gp, b, elem, isStressTensorSymmetric);
    answer.beProductOf(b, u);
}

 

void
BaseSecondGradientElement :: computeMicromorphicVars(FloatArray &micromorphicVar,FloatArray &micromorphicVarGrad, IntArray IdMask_m, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix N_m, B_m;
    FloatArray d_m;
 
    this->computeMicromorphicNMatrixAt(gp, N_m);
    this->computeMicromorphicBMatrixAt(gp, B_m);
    /// @todo generalization for general micromorphic continua -- should be parameter of this function?
    this->giveElement()->computeVectorOf(IdMask_m, VM_Total, tStep, d_m);
    micromorphicVar.beProductOf(N_m, d_m);
    micromorphicVarGrad.beProductOf(B_m, d_m);
}


void
BaseSecondGradientElement :: computeLagrangianMultipliers(FloatArray &lagrangianMultipliers, IntArray IdMask_l, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix N_l;
    FloatArray d_l;
 
    this->computeLagrangianMultiplierNMatrixAt(gp, N_l);
    this->giveElement()->computeVectorOf(IdMask_l, VM_Total, tStep, d_l);
    lagrangianMultipliers.beProductOf(N_l, d_l);
}


void
BaseSecondGradientElement :: giveMicromorphicInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveElement();
    FloatMatrix N_m, B_m;
    answer.resize(this->giveNumberOfMicromorphicDofs());
    FloatArray vStress, vMicromorphicStress, vMicromorphicStressGrad, dE;

    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      if ( useUpdatedGpRecord == 1 ) {
        vMicromorphicStressGrad = static_cast< SecondGradientMaterialStatus* >( gp->giveMaterialStatus() )->giveTempMicromorphicStressGrad();
                               
      } else {
	this->computeGeneralizedStressVectors(vStress, vMicromorphicStress, vMicromorphicStressGrad, dE, gp, tStep);      
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
BaseSecondGradientElement :: giveStandardInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatArray BS, vStress, s, S, dE;
    FloatMatrix B;
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      if ( useUpdatedGpRecord == 1 ) {
        vStress = static_cast< SecondGradientMaterialStatus * >( gp->giveMaterialStatus() )->giveTempStressVector();
      } else {
	computeGeneralizedStressVectors(vStress, s, S, dE, gp, tStep);      
      }
      SecondGradientMaterialExtensionInterface *secondGradientMat = dynamic_cast< SecondGradientMaterialExtensionInterface * >(cs->giveMaterialInterface(SecondGradientMaterialExtensionInterfaceType, gp) );
        if ( !secondGradientMat ) {
            OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }
        // Compute nodal internal forces at nodes as f = B^T*Stress dV
      double dV  = elem->computeVolumeAround(gp);
      bool isStressTensorSymmetric = secondGradientMat->isStressTensorSymmetric();
      this->computeB_uMatrixAt(gp, B, elem, isStressTensorSymmetric);
      
      BS.beTProductOf(B, vStress);
      answer.add(dV, BS);
    }
}



void
BaseSecondGradientElement :: giveLagrangianMultipliersInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatArray BS, vStress, s, S, dE;
    FloatMatrix N_l, B_u, N_m;
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      SecondGradientMaterialExtensionInterface *secondGradientMat = dynamic_cast< SecondGradientMaterialExtensionInterface * >(cs->giveMaterialInterface(SecondGradientMaterialExtensionInterfaceType, gp) );
      if ( !secondGradientMat ) {
	OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
      }

      /* if ( useUpdatedGpRecord == 1 ) {
        vStress = static_cast< MicromorphicMaterialStatus * >( gp->giveMaterialStatus() )->giveTempStressVector();
	} else {*/
	computeGeneralizedStressVectors(vStress, s, S, dE, gp, tStep);      
	//      }
      // Compute nodal internal forces at nodes as f = \int (Nl^T(dE) dV
      this->computeLagrangianMultiplierNMatrixAt(gp, N_l);
      double dV  = elem->computeVolumeAround(gp);
      FloatArray NldE;
      NldE.beTProductOf(N_l, dE);
      answer.add(dV,NldE);
    }
}


void
BaseSecondGradientElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    answer.resize(this->giveNumberOfDofs());
    answer.zero();
    FloatArray answerU(this->giveNumberOfDisplacementDofs());
    answer.zero();
    FloatArray answerM(this->giveNumberOfMicromorphicDofs());
    answerM.zero();
    FloatArray answerL(this->giveNumberOfLagrangianMultipliersDofs());
    answerL.zero();


    this->giveLagrangianMultipliersInternalForcesVector(answerL, tStep, 0);   
    this->giveMicromorphicInternalForcesVector(answerM, tStep, 0);
    this->giveStandardInternalForcesVector(answerU, tStep, 1);

    
   
    answer.assemble(answerU, locationArray_u);
    answer.assemble(answerM, locationArray_m);
    answer.assemble(answerL, locationArray_l);

}

void
BaseSecondGradientElement :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
 
    FloatArray localForces(this->giveNumberOfDisplacementDofs());
    FloatArray nlForces(this->giveNumberOfMicromorphicDofs());
    answer.resize(this->giveNumberOfDofs());

    this->computeLocForceLoadVector(localForces, tStep, mode);

    answer.assemble(localForces, locationArray_u);
    answer.assemble(nlForces, locationArray_m);
}


/************************************************************************/
void
BaseSecondGradientElement :: computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    FloatMatrix T;
    NLStructuralElement *elem = this->giveElement();
    //@todo check this
    //elem->computeLocalForceLoadVector(answer, tStep, mode);

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
BaseSecondGradientElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    answer.resize(this->giveNumberOfDofs(), this->giveNumberOfDofs());
    answer.zero();

    FloatMatrix Kuu, Kul, Kmm, Kml,  Klu, Klm;

    this->computeStiffnessMatrix_uu(Kuu, rMode, tStep);
    this->computeStiffnessMatrix_ul(Kul, rMode, tStep);

    this->computeStiffnessMatrix_mm(Kmm, rMode, tStep);
    this->computeStiffnessMatrix_ml(Kml, rMode, tStep);
 

    this->computeStiffnessMatrix_lu(Klu, rMode, tStep);
    this->computeStiffnessMatrix_lm(Klm, rMode, tStep);



    answer.assemble(Kuu, locationArray_u);
    answer.assemble(Kul, locationArray_u, locationArray_l);

    answer.assemble(Kmm, locationArray_m);
    answer.assemble(Kml, locationArray_m, locationArray_l);

    answer.assemble(Klu, locationArray_l, locationArray_u);
    answer.assemble(Klm, locationArray_l, locationArray_m);


}


void
BaseSecondGradientElement :: computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatMatrix B, D, DB;
    bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.clear();
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {

        SecondGradientMaterialExtensionInterface *secondGradientMat = dynamic_cast< SecondGradientMaterialExtensionInterface * >(
            cs->giveMaterialInterface(SecondGradientMaterialExtensionInterfaceType, gp) );
        if ( !secondGradientMat ) {
            OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }
	bool isStressTensorSymmetric = secondGradientMat->isStressTensorSymmetric();
	this->computeB_uMatrixAt(gp,B, elem, isStressTensorSymmetric);

        secondGradientMat->giveSecondGradientMatrix_dSigdUgrad(D, rMode, gp, tStep);
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
BaseSecondGradientElement :: computeStiffnessMatrix_ul(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    double dV;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatMatrix B, N_l, D, DN_l;

    answer.clear();

    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
        SecondGradientMaterialExtensionInterface *secondGradientMat = dynamic_cast< SecondGradientMaterialExtensionInterface * >(cs->giveMaterialInterface(SecondGradientMaterialExtensionInterfaceType, gp) );
        if ( !secondGradientMat ) {
            OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }

        secondGradientMat->giveSecondGradientMatrix_dEdPhi(D, rMode, gp, tStep);

	bool isStressTensorSymmetric = secondGradientMat->isStressTensorSymmetric();
	this->computeB_uMatrixAt(gp,B, elem, isStressTensorSymmetric);
        this->computeLagrangianMultiplierNMatrixAt(gp, N_l);

	dV = elem->computeVolumeAround(gp);
        DN_l.beProductOf(D, N_l);
        answer.plusProductUnsym(B, DN_l, dV);
    }

}




void
BaseSecondGradientElement :: computeStiffnessMatrix_mm(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    double dV;
    FloatMatrix lStiff;
    FloatMatrix B_m, dMdPhiGrad,dMdPhiGrad_B;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    answer.clear();

    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
        SecondGradientMaterialExtensionInterface *secondGradientMat = dynamic_cast< SecondGradientMaterialExtensionInterface * >(cs->giveMaterialInterface(SecondGradientMaterialExtensionInterfaceType, gp) );
        if ( !secondGradientMat ) {
            OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }

	secondGradientMat->giveSecondGradientMatrix_dMdPhiGrad(dMdPhiGrad, rMode, gp, tStep);
        this->computeMicromorphicBMatrixAt(gp, B_m);
        dV = elem->computeVolumeAround(gp);

	dMdPhiGrad_B.beTProductOf(dMdPhiGrad, B_m);
	answer.plusProductUnsym(B_m,dMdPhiGrad_B,dV);

    }
}



void
BaseSecondGradientElement :: computeStiffnessMatrix_ml(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    double dV;
    NLStructuralElement *elem = this->giveElement();
    FloatMatrix N_l, N_m;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();

    answer.clear();


    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
        SecondGradientMaterialExtensionInterface *secondGradientMat = dynamic_cast< SecondGradientMaterialExtensionInterface * >(cs->giveMaterialInterface(SecondGradientMaterialExtensionInterfaceType, gp) );
        if ( !secondGradientMat ) {
	    OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }

        this->computeMicromorphicNMatrixAt(gp, N_m);
	this->computeLagrangianMultiplierNMatrixAt(gp, N_m);

	dV = elem->computeVolumeAround(gp);
        answer.plusProductUnsym(N_m,N_l, dV);
    }

}






void
BaseSecondGradientElement :: computeStiffnessMatrix_lu(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{


NLStructuralElement *elem = this->giveElement();
    double dV;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    FloatMatrix B, N_l, D, DB;

    answer.clear();

    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
        SecondGradientMaterialExtensionInterface *secondGradientMat = dynamic_cast< SecondGradientMaterialExtensionInterface * >(cs->giveMaterialInterface(SecondGradientMaterialExtensionInterfaceType, gp) );
        if ( !secondGradientMat ) {
            OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
        }

        secondGradientMat->giveSecondGradientMatrix_dEdPhi(D, rMode, gp, tStep);

	bool isStressTensorSymmetric = secondGradientMat->isStressTensorSymmetric();
	this->computeB_uMatrixAt(gp,B, elem, isStressTensorSymmetric);
        this->computeLagrangianMultiplierNMatrixAt(gp, N_l);

	dV = elem->computeVolumeAround(gp);
        DB.beProductOf(D, B);
        answer.plusProductUnsym(N_l, DB, dV);
    }


   

}


void
BaseSecondGradientElement :: computeStiffnessMatrix_lm(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    double dV;
    FloatMatrix N_m, N_l;
    answer.clear();

    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
        this->computeMicromorphicNMatrixAt(gp, N_m);
        this->computeLagrangianMultiplierNMatrixAt(gp, N_l);
        dV = elem->computeVolumeAround(gp);
	answer.plusProductUnsym(N_l,N_m,-dV);
    }
}







void
BaseSecondGradientElement :: computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
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
BaseSecondGradientElement :: computeMassMatrix_uu(FloatMatrix &answer, TimeStep *tStep)
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
BaseSecondGradientElement :: computeMassMatrix_mm(FloatMatrix &answer, TimeStep *tStep)
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


void
BaseSecondGradientElement :: postInitialize()
{
  IntArray IdMask_u, IdMask_m, IdMask_l;
  this->giveDofManDofIDMask_u( IdMask_u );
  this->giveDofManDofIDMask_m( IdMask_m );
  this->giveDofManDofIDMask_m( IdMask_l );
  this->giveLocationArrayOfDofIDs(locationArray_u,locationArray_m, locationArray_l, EModelDefaultEquationNumbering(), IdMask_u, IdMask_m, IdMask_l);
}




} // end namespace oofem


