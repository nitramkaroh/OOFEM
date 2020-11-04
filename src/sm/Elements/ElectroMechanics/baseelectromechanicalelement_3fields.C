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


#include "../sm/Elements/ElectroMechanics/baseelectromechanicalelement_3fields.h"
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
BaseElectroMechanicalElement_3Fields :: BaseElectroMechanicalElement_3Fields(int n, Domain *domain)
{
}

    

void
BaseElectroMechanicalElement_3Fields :: giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_phi, IntArray &locationArray_d, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u,const IntArray &dofIdArray_phi, const IntArray &dofIdArray_d )
{
    // Routine to extract the location array of an element for given dofid array.
    locationArray_u.clear();
    locationArray_phi.clear();
    locationArray_d.clear();

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
      for(int j = 1; j <= dofIdArray_phi.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_phi.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_phi.followedBy( k + itt);
	}
	itt++;
      }
      for(int j = 1; j <= dofIdArray_d.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_d.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_d.followedBy( k + itt);
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
      for(int j = 1; j <= dofIdArray_phi.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_phi.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_phi.followedBy( k + itt);
	  itt++;
	}

      }
      for(int j = 1; j <= dofIdArray_d.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_d.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_d.followedBy( k + itt);
	  itt++;
	}

      }

      k += dMan->giveNumberOfDofs( );
    }
    
    

    
}

 

void
BaseElectroMechanicalElement_3Fields :: compute_FirstPKStressVector_ElectricFieldVector(FloatArray &P, FloatArray &E, GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveStructuralElement();
    SimpleElectroMechanicalCrossSection_3Fields *cs = this->giveCrossSection();

    FloatArray electricDisplacemenet;
    if( elem->giveGeometryMode() == 0) {
      /*FloatArray strain;
      this->computeStrainVector(strain, gp, tStep);
      this->computePressure(pressure,gp, tStep);   
      mixedPressureMat->giveRealStressVector(stress, gp, strain,pressure, tStep);
      */
    } else {
      FloatArray F, electricDisplacement;
      elem->computeDeformationGradientVector(F, gp, tStep);
      this->computeElectricDisplacementVector(electricDisplacement, gp, tStep);   
      cs->give_FirstPKStressVector_ElectricFieldVector(P, E, gp, F, electricDisplacement, tStep);
    }
    
}


void
BaseElectroMechanicalElement_3Fields :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
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
BaseElectroMechanicalElement_3Fields :: computeElectricPotentialGradientVector(FloatArray &answer,GaussPoint *gp, TimeStep *tStep)
{
    IntArray IdMask_phi;
    FloatArray d_phi;
    FloatMatrix B_phi;
    this->giveDofManDofIDMask_phi( IdMask_phi );
    this->giveStructuralElement()->computeVectorOf(IdMask_phi, VM_Total, tStep, d_phi);
    this->computeElectricPotentialBmatrixAt(gp, B_phi);  
    answer.beProductOf(B_phi,d_phi);
}



void
BaseElectroMechanicalElement_3Fields :: computeElectricDisplacementVector(FloatArray &answer,GaussPoint *gp, TimeStep *tStep)
{
    IntArray IdMask_d;
    FloatArray d_d;
    FloatMatrix N_d;
    this->giveDofManDofIDMask_d( IdMask_d );
    this->giveStructuralElement()->computeVectorOf(IdMask_d, VM_Total, tStep, d_d);
    this->computeElectricDisplacementNmatrixAt(gp, N_d);  
    answer.beProductOf(N_d,d_d);
}
  



void
BaseElectroMechanicalElement_3Fields :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveStructuralElement();
    FloatArray BP, BD, vP, vE, NE, gradPhi, electricDisplacement;
    FloatMatrix B_u, B_phi, N_d;
   
    answer.resize(this->giveNumberOfDofs());
    answer.zero();
    
    FloatArray answer_u(this->giveNumberOfDisplacementDofs());
    answer_u.zero();
    FloatArray answer_phi(this->giveNumberOfElectricPotentialDofs());
    answer_phi.zero();
    FloatArray answer_d(this->giveNumberOfElectricDisplacementDofs());
    answer_d.zero();
    
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      this->compute_FirstPKStressVector_ElectricFieldVector(vP, vE, gp, tStep);
      this->computeElectricPotentialGradientVector(gradPhi, gp, tStep);
      this->computeElectricDisplacementVector(electricDisplacement, gp, tStep);
      //
      double dV  = elem->computeVolumeAround(gp);
      // Compute nodal internal forces at nodes as f_u = \int_V B^T*vP dV
      elem->computeBHmatrixAt(gp, B_u, tStep, 0);
      BP.beTProductOf(B_u, vP);
      answer_u.add(dV, BP);
      // Compute nodal internal forces at nodes as f_\phi = \int B^T* vD dV     
      this->computeElectricPotentialBmatrixAt(gp, B_phi);      
      BD.beTProductOf(B_phi, electricDisplacement);
      answer_phi.add(dV, BD);
      // Compute nodal internal forces at nodes as f_d = \int N^T(E + gradPhi) dB
      this->computeElectricDisplacementNmatrixAt(gp, N_d);
      vE.add(gradPhi);
      NE.beTProductOf(N_d, vE);
      answer_d.add(dV, NE);
      
    }

    answer.assemble(answer_u, locationArray_u);
    answer.assemble(answer_phi, locationArray_phi);
    answer.assemble(answer_d, locationArray_d);

}




void
BaseElectroMechanicalElement_3Fields :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
 
    FloatArray localForces(this->giveNumberOfDisplacementDofs());
    answer.resize(this->giveNumberOfDofs());
    this->computeLocForceLoadVector(localForces, tStep, mode);
    answer.assemble(localForces, locationArray_u);

}


/************************************************************************/
void
BaseElectroMechanicalElement_3Fields :: computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
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
BaseElectroMechanicalElement_3Fields :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
  //set displacement and nonlocal location array
    answer.resize(this->giveNumberOfDofs(), this->giveNumberOfDofs());
    answer.zero();


    NLStructuralElement *elem = this->giveStructuralElement();
    SimpleElectroMechanicalCrossSection_3Fields *cs = this->giveCrossSection();

    FloatMatrix B_u, B_phi, N_d, Duu, Dud, Ddd, DuuBu, DueNd, DudNd, DddNd, DudBd;
    FloatMatrix Kuu, Kud, Kdu, Kdd, Kde, Ked;
    
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      this->computeElectricPotentialBmatrixAt(gp, B_phi);
      this->computeElectricDisplacementNmatrixAt(gp, N_d);  
      if( elem->giveGeometryMode() == 0) {
	/*FloatArray strain;
	  this->computeStrainVector(strain, gp, tStep);
	  this->computePressure(pressure,gp, tStep);   
	  mixedPressureMat->giveRealStressVector(stress, gp, strain,pressure, tStep);
	*/
      } else {
	elem->computeBHmatrixAt(gp, B_u);
	cs->giveElectroMechanicalConstitutiveMatrix_dPdF(Duu, rMode, gp, tStep);
	cs->giveElectroMechanicalConstitutiveMatrix_dPdD(Dud, rMode, gp, tStep);
	cs->giveElectroMechanicalConstitutiveMatrix_dEdD(Ddd, rMode, gp, tStep);
      }
      double dV  = elem->computeVolumeAround(gp);
      DuuBu.beProductOf(Duu, B_u);
      DudNd.beProductOf(Dud, N_d);
      DddNd.beProductOf(Ddd, N_d);      
      Kuu.plusProductUnsym(B_u, DuuBu, dV);
      Kud.plusProductUnsym(B_u, DudNd, dV);
      Kdd.plusProductUnsym(N_d, DddNd, dV);
      Ked.plusProductUnsym(B_phi, N_d, dV);

    }
      

    Kdu.beTranspositionOf(Kud);
    Kde.beTranspositionOf(Ked);    
    
    answer.assemble(Kuu, locationArray_u);
    answer.assemble(Kud, locationArray_u, locationArray_d);
    answer.assemble(Kdu, locationArray_d, locationArray_u);
    answer.assemble(Ked, locationArray_phi, locationArray_d);
    answer.assemble(Kde, locationArray_d, locationArray_phi);    
    answer.assemble(Kdd, locationArray_d);

}


SimpleElectroMechanicalCrossSection_3Fields*
BaseElectroMechanicalElement_3Fields :: giveCrossSection()
// Returns the crossSection of the receiver.
{
  //NLStructuralElement *elem = this->giveElement();
  return static_cast< SimpleElectroMechanicalCrossSection_3Fields* >( this->giveStructuralElement()->giveCrossSection() );
}



IRResultType
BaseElectroMechanicalElement_3Fields :: initializeFrom(InputRecord *ir)
{
  // @todo Is this function necessary???

    return IRRT_OK;
}

void
BaseElectroMechanicalElement_3Fields :: updateInternalState(TimeStep *tStep)
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
BaseElectroMechanicalElement_3Fields :: postInitialize()
{
  IntArray IdMask_u, IdMask_phi, IdMask_d;
  this->giveDofManDofIDMask_u( IdMask_u );
  this->giveDofManDofIDMask_phi( IdMask_phi );
  this->giveDofManDofIDMask_d( IdMask_d );
  this->giveLocationArrayOfDofIDs(locationArray_u,locationArray_phi,locationArray_d, EModelDefaultEquationNumbering(), IdMask_u, IdMask_phi, IdMask_d);
  
}




} // end namespace oofem


