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


#include "../sm/Elements/ElectroMechanics/baseelectromechanicalelement.h"
#include "../sm/Materials/ElectroMechanics/electromechanicalmaterialextensioninterface.h"

#include "../sm/Materials/structuralms.h"

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
BaseElectroMechanicalElement :: BaseElectroMechanicalElement(int n, Domain *domain)
{
}

    

void
BaseElectroMechanicalElement :: giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_e, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u,const IntArray &dofIdArray_e )
{
    // Routine to extract the location array of an element for given dofid array.
    locationArray_u.clear();
    locationArray_e.clear();
    NLStructuralElement *el = this->giveStructuralElement();
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
      for(int j = 1; j <= dofIdArray_e.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_e.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_e.followedBy( k + itt);
	}
	itt++;
      }
      k += dMan->giveNumberOfDofs( );
    }

    for ( int i = 1; i <= el->giveNumberOfInternalDofManagers(); i++ ) {
      DofManager *dMan = el->giveInternalDofManager( i );
      int itt = 1;
      for(int j = 1; j <= dofIdArray_u.giveSize( ); j++) {
	if(dMan->hasDofID( (DofIDItem) dofIdArray_u.at( j ) )) {
	  //  Dof *d = dMan->giveDofWithID( dofIdArray_u.at( j ) );
	  locationArray_u.followedBy( k + itt);
	  itt++;
	}

      }
      for(int j = 1; j <= dofIdArray_e.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_e.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_e.followedBy( k + itt);
	  itt++;
	}

      }
      k += dMan->giveNumberOfDofs( );
    }
    
    

    
}

 

void
BaseElectroMechanicalElement :: compute_FirstPKStressVector_ElectricDisplacementVector(FloatArray &P, FloatArray &D, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveStructuralElement();
    SimpleElectroMechanicalCrossSection *cs = this->giveCrossSection();

    FloatArray electricField;
    if( elem->giveGeometryMode() == 0) {
      /*FloatArray strain;
      this->computeStrainVector(strain, gp, tStep);
      this->computePressure(pressure,gp, tStep);   
      mixedPressureMat->giveRealStressVector(stress, gp, strain,pressure, tStep);
      */
    } else {
      FloatArray F;
      elem->computeDeformationGradientVector(F, gp, tStep);
      this->computeElectricField(electricField, gp, tStep);   
      cs->give_FirstPKStressVector_ElectricDisplacementVector(P, D, gp, F, electricField, tStep);
    }
    
}


void
BaseElectroMechanicalElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray d_u;
    FloatMatrix B;
    IntArray IdMask_u;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->giveStructuralElement()->computeVectorOf(IdMask_u, VM_Total, tStep, d_u);
    this->computeDisplacementFieldBmatrixAt(gp,B);
    answer.beProductOf(B, d_u);
}


void
BaseElectroMechanicalElement :: computeElectricField(FloatArray &answer,GaussPoint *gp, TimeStep *tStep)
{
    IntArray IdMask_e;
    FloatArray d_e;
    FloatMatrix B_e;
    this->giveDofManDofIDMask_e( IdMask_e );
    this->giveStructuralElement()->computeVectorOf(IdMask_e, VM_Total, tStep, d_e);
    this->computeElectricFieldBmatrixAt(gp, B_e);  
    answer.beProductOf(B_e,d_e);
}



void
BaseElectroMechanicalElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveStructuralElement();
    FloatArray BS, vStress, electricDisplacement, BD;
    FloatMatrix B_u, B_e;
   
    answer.resize(this->giveNumberOfDofs());
    answer.zero();
    
    FloatArray answer_u(this->giveNumberOfDisplacementDofs());
    answer_u.zero();
    FloatArray answer_e(this->giveNumberOfElectricDofs());
    answer_e.zero();
   
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      this->compute_FirstPKStressVector_ElectricDisplacementVector(vStress, electricDisplacement, gp, tStep);
      // Compute nodal internal forces at nodes as f = B^T*Stress dV
      double dV  = elem->computeVolumeAround(gp);
      elem->computeBHmatrixAt(gp, B_u, tStep, 0);
      this->computeElectricFieldBmatrixAt(gp, B_e);      
      BS.beTProductOf(B_u, vStress);
      answer_u.add(dV, BS);
      BD.beTProductOf(B_e, electricDisplacement);
      answer_e.add(dV, BD);
    }

    answer.assemble(answer_u, locationArray_u);
    answer.assemble(answer_e, locationArray_e);

}




void
BaseElectroMechanicalElement :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
 
    FloatArray localForces(this->giveNumberOfDisplacementDofs());
    answer.resize(this->giveNumberOfDofs());
    this->computeLocForceLoadVector(localForces, tStep, mode);
    answer.assemble(localForces, locationArray_u);

}


/************************************************************************/
void
BaseElectroMechanicalElement :: computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    FloatMatrix T;
    NLStructuralElement *elem = this->giveStructuralElement();
    //@todo check this
    //    elem->computeLocalForceLoadVector(answer, tStep, mode);

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
BaseElectroMechanicalElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    answer.resize(this->giveNumberOfDofs(), this->giveNumberOfDofs());
    answer.zero();


    NLStructuralElement *elem = this->giveStructuralElement();
    SimpleElectroMechanicalCrossSection *cs = this->giveCrossSection();

    FloatMatrix B_u, B_e, Duu, Due, Dee, DuuBu, DueBe, DeeBe;
    FloatMatrix Kuu, Kue, Keu, Kee;
    
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      this->computeElectricFieldBmatrixAt(gp, B_e);  
      if( elem->giveGeometryMode() == 0) {
	/*FloatArray strain;
	  this->computeStrainVector(strain, gp, tStep);
	  this->computePressure(pressure,gp, tStep);   
	  mixedPressureMat->giveRealStressVector(stress, gp, strain,pressure, tStep);
	*/
      } else {
	elem->computeBHmatrixAt(gp, B_u);
	cs->giveElectroMechanicalConstitutiveMatrix_dPdF(Duu, rMode, gp, tStep);
	cs->giveElectroMechanicalConstitutiveMatrix_dPdE(Due, rMode, gp, tStep);
	cs->giveElectroMechanicalConstitutiveMatrix_dDdE(Dee, rMode, gp, tStep);

      }
      double dV  = elem->computeVolumeAround(gp);
      DuuBu.beProductOf(Duu, B_u);
      DueBe.beProductOf(Due, B_e);
      DeeBe.beProductOf(Dee, B_e);      
      Kuu.plusProductUnsym(B_u, DuuBu, dV);
      Kue.plusProductUnsym(B_u, DueBe, dV);
      Kee.plusProductUnsym(B_e, DeeBe, dV);

    }
      

    Keu.beTranspositionOf(Kue);    
    
    answer.assemble(Kuu, locationArray_u);
    answer.assemble(Kue, locationArray_u, locationArray_e);
    answer.assemble(Keu, locationArray_e, locationArray_u);
    answer.assemble(Kee, locationArray_e);

}


SimpleElectroMechanicalCrossSection*
BaseElectroMechanicalElement :: giveCrossSection()
// Returns the crossSection of the receiver.
{
  //NLStructuralElement *elem = this->giveElement();
  return static_cast< SimpleElectroMechanicalCrossSection* >( this->giveStructuralElement()->giveCrossSection() );
}



IRResultType
BaseElectroMechanicalElement :: initializeFrom(InputRecord *ir)
{
  // @todo Is this function necessary???

    return IRRT_OK;
}

void
BaseElectroMechanicalElement :: updateInternalState(TimeStep *tStep)
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
BaseElectroMechanicalElement :: postInitialize()
{
  IntArray IdMask_u, IdMask_e;
  this->giveDofManDofIDMask_u( IdMask_u );
  this->giveDofManDofIDMask_e( IdMask_e );
  this->giveLocationArrayOfDofIDs(locationArray_u,locationArray_e, EModelDefaultEquationNumbering(), IdMask_u, IdMask_e);
  
}




} // end namespace oofem


