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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include "meandilatationalmethodelementinterface.h"
#include "domain.h"
#include "structuralcrosssection.h"
#include "structuralms.h"

		  


namespace oofem {


// constructor
MeanDilatationalMethodElementExtensionInterface :: MeanDilatationalMethodElementExtensionInterface(Domain *d)  : Interface()
{
	
}


void 
MeanDilatationalMethodElementExtensionInterface :: computeMeanDilatationalBHmatrix(FloatMatrix &answer, TimeStep *tStep, NLStructuralElement *elem)
{

  


}

void
MeanDilatationalMethodElementExtensionInterface :: computeStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseMode rMode, TimeStep *tStep, NLStructuralElement *elem)
{
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);

    answer.clear();

    if ( !this->isActivated(tStep) ) {
        return;
    }

    // add mean dilatational method contribution
    // first, compute mean Bmatrix and mean J
    FloatMatrix meanBh;
    doubel Jbar = this->computeMeanDilatationalBHmatrixAt(meanBh, tStep, elem)



    // Compute matrix from material stiffness (total stiffness for small def.) - B^T * dS/dE * B
    if ( elem->numberOfIntegrationRules == 1 ) {
        FloatMatrix B, D, Dp, Dvol, DB;
        for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {
          if ( nlGeometry == 1 ) {
	    if ( elem->giveDomain()->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
	      elem->computeBmatrixAt(gp, B, tStep);
	      IncompressibleMaterialExtensionInterface *imat = dynamic_cast< IncompressibleMaterialExtensionInterface* >(cs->giveMaterialInterface(IncompressibleMaterialExtensionInterfaceType, gp) );
	      imat->giveDeviatoric3dMaterialStiffnessMatrix_dCde(D, rMode, gp, tStep);
	      imat-> givePressure3dMaterialStiffnessMatrix_dCde(Dp, rMode, gp, tStep);
	      D.add(Dp)     
	    } else { // Material stiffness dP/dF
	      OOFEM_ERROR("Only actualized lagrangian formulation is currently supported.")
	    }
	  } else { // small strains
	    OOFEM_ERROR("Only large strain formulation is currently supported.")
	  }
	    
	  
	  double dV = elem->computeVolumeAround(gp);	  
	  // add initial stress stiffness 
	  FloatMatrix Bh, sI_Bh, sI;
	  cs->giveInitialStiffnessMatrix_Cauchy(sI, rMode, gp, tStep);
	  elem->computeBHmatrixAt(gp, Bh,tStep,1);
	  sI_Bh.beProductOf(sI, Bh);
	  if ( matStiffSymmFlag ) {
	    answer.plusProductSymmUpper(Bh, sI_Bh, dV);
	  } else {
	    answer.plusProductUnsym(Bh, sI_Bh, dV);
	  }
	  
	  DB.beProductOf(D, B);
	  if ( matStiffSymmFlag ) {
	    answer.plusProductSymmUpper(B, DB, dV);
	  } else {
	    answer.plusProductUnsym(B, DB, dV);
	  }




        }
    
        

    } else {
      OOFEM_ERROR("More than one integration rule is not currently supported.")
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}











} // end namespace oofem
