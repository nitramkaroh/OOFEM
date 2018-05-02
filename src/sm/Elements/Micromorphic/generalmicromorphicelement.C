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


#include <cstdio>

namespace oofem {
BaseMicromorphicElement :: BaseMicromorphicElement()
{
}

    

void
BaseMicromorphicElement :: computeLocationArrayOfDofIDs( const IntArray &dofIdArray, IntArray &answer )
{
    // Routine to extract compute the location array an element given an dofid array.
    answer.clear();
    NLStructuralElement *el = this->giveElement();
    int k = 0;
    for(int i = 1; i <= el->giveNumberOfDofManagers(); i++) {
        DofManager *dMan = el->giveDofManager( i );
        for(int j = 1; j <= dofIdArray.giveSize( ); j++) {

            if(dMan->hasDofID( (DofIDItem) dofIdArray.at( j ) )) {
                Dof *d = dMan->giveDofWithID( dofIdArray.at( j ) );
                answer.followedBy( k + d->giveNumber( ) );
            }
        }
        k += dMan->giveNumberOfDofs( );
    }
}
  
void
BaseMicromorphicElement :: computeGeneralizedStressVectors(FloatArray &sigma, FloatArray &s, FloatArray &S, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    MicromorphicMaterialExtensionInterface *micromorphMat = static_cast< MicromorphicMaterialExtensionInterface * >(cs->giveMaterialInterface(MicromorphicMaterialExtensionInterfaceType, gp) );
    if ( !micromorphMat ) {
        OOFEM_ERROR("Material doesn't implement the required Micromorphic interface!");
    }

    FloatArray displacementGradient, micromorphicVar, micromorphicVarGrad;
    
    this->computeDisplacementGradient(displacementGradient, gp, tStep);
    this->computeMicromorphicVars(micromorphicVar,micromorphicVarGrad,gp,tStep);   

    micromorphMat->giveGeneralizedStressVectors(sigma, s, S, gp, displacementGradient, micromorphicVar, micromorphicVarGrad, tStep);
    
}


void
BaseMicromorphicElement :: computeDisplacementGradient(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    FloatMatrix b;
    NLStructuralElement *elem = this->giveElement();
    elem->computeVectorOf({D_u, D_v,D_w}, VM_Total, tStep, u);
    if(isStressTensorSymmetric)
      elem->computeBmatrixAt(gp, b);
    else
      elem->computeBHmatrixAt(gp, b);

    answer.beProductOf(b, u);
}

 

void
BaseMicromorphicElement :: computeMicromorphicVars(FloatArray &micromorphicVar,FloatArray &micromorphicVarGrad, GaussPoint *gp, TimeStep *tStep)
{
  FloatMatrix N_m, B_m;
  FloatArray u_m;
 
    this->computeMicromorphicNMatrixAt(gp, N_m);
    this->computeMicromorphicBMatrixAt(gp, B_m);
    /// @todo generalization for general micromorphic continua -- should be parameter of this function
    this->giveElement()->computeVectorOf({G_0}, VM_Total, tStep, u_m);
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
        vMicromorphicStressGrad = static_cast< MicromorphicMaterialStatus * >( gp->giveMaterialStatus() )->giveTempMicromorphicStressGrad();
                               
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
    FloatArray BS, vStress, s, S;
    FloatMatrix B;
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      if ( useUpdatedGpRecord == 1 ) {
        vStress = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveTempStressVector();
      } else {
	computeGeneralizedStressVectors(vStress, s, S, gp, tStep);      
      }

        // Compute nodal internal forces at nodes as f = B^T*Stress dV
      double dV  = elem->computeVolumeAround(gp);

      if(isStressTensorSymmetric)
	elem->computeBmatrixAt(gp, B);
      else
      elem->computeBHmatrixAt(gp, B);
      
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
   
    this->giveMicromorphicInternalForcesVector(answerM, tStep, useUpdatedGpRecord);
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
	
	if(isStressTensorSymmetric)
	  elem->computeBmatrixAt(gp, B);
	else
	  elem->computeBHmatrixAt(gp, B);
        
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

	if(isStressTensorSymmetric)
	  elem->computeBmatrixAt(gp, B);
	else
	  elem->computeBHmatrixAt(gp, B);


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

	if(isStressTensorSymmetric)
	  elem->computeBmatrixAt(gp, B);
	else
	  elem->computeBHmatrixAt(gp, B);


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



IRResultType
BaseMicromorphicElement :: initializeFrom(InputRecord *ir)
{
  // @todo Is this function necessary???

    return IRRT_OK;
}




void
BaseMicromorphicElement :: postInitialize()
{
  IntArray IdMask_u, IdMask_m;
  this->giveDofManDofIDMask_u( IdMask_u );
  this->giveDofManDofIDMask_m( IdMask_m );
  this->giveElement()->computeLocationArrayOfDofIDs(locationArray_u,EModelDefaultEquationNumbering(),&IdMask_u);
  this->giveElement()->computeLocationArrayOfDofIDs(locationArray_m,EModelDefaultEquationNumbering(),&IdMask_m);
}




} // end namespace oofem


