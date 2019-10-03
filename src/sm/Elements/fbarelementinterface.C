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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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


#include "fbarelementinterface.h"
#include "domain.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "gausspoint.h"		  
#include "engngm.h"
#include "../sm/Materials/structuralms.h"
#include "mathfem.h"

namespace oofem {


// constructor
FbarElementExtensionInterface :: FbarElementExtensionInterface(Domain *d)  : Interface()
{
  FbarFlag = 0;
  initialStepFlag = 1;
}

IRResultType
FbarElementExtensionInterface :: initializeFrom(InputRecord *ir)
{  
  IRResultType result;              // Required by IR_GIVE_FIELD macro
  // read the characteristic length
  IR_GIVE_OPTIONAL_FIELD(ir, FbarFlag, _IFT_FbarElementExtensionInterface_fbarflag);
  return IRRT_OK;
}



void
FbarElementExtensionInterface :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, NLStructuralElement *elem, ValueModeType modeType)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    elem->computeVectorOf({D_u, D_v, D_w}, modeType, tStep, u); // solution vector
    //@todo missing initial displacement
    // Displacement gradient H = du/dX
    FloatMatrix B;
    elem->computeBHmatrixAt(gp, B, tStep, 0);
    answer.beProductOf(B, u);

    // Deformation gradient F = H + I
    MaterialMode matMode = gp->giveMaterialMode();
    if ( matMode == _3dMat || matMode == _PlaneStrain ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
        answer.at(3) += 1.0;
    } else if ( matMode == _PlaneStress ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
    } else if ( matMode == _1dMat ) {
        answer.at(1) += 1.0;
    } else {
        OOFEM_ERROR("MaterialMode is not supported yet (%s)", __MaterialModeToString(matMode) );
    }
}


double
FbarElementExtensionInterface :: computeFbarDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, NLStructuralElement *nlElem, ValueModeType valueModeType)
{
  // create a fictitious integration point in the element centroid
  GaussPoint *gp0;
  FloatArray centroidCoords;
  MaterialMode matMode = nlElem->giveMaterialMode();
  nlElem->giveElementParametricCentroid(centroidCoords);
  gp0 = new GaussPoint( NULL, 1, centroidCoords, 1.0, matMode);
  // compute compatible deformation gradient
  FloatArray reducedvF,vF,vFbar;;
  FloatMatrix F, Fbar;
  this->computeDeformationGradientVector(reducedvF, gp, tStep, nlElem, valueModeType);
  StructuralMaterial :: giveFullVectorFormF(vF, reducedvF, matMode); // 9 components
  F.beMatrixForm(vF);
  // compute deformation gradient in the element centroid
  FloatArray reducedvF0,vF0;
  FloatMatrix F0;
  this->computeDeformationGradientVector(reducedvF0, gp0, tStep, nlElem, valueModeType);
  StructuralMaterial :: giveFullVectorFormF(vF0, reducedvF0, matMode); // 9 components
  F0.beMatrixForm(vF0);

  double J0, J;
  
  if(matMode == _PlaneStrain) {
    FloatMatrix Fp;
    Fp.beSubMatrixOf(F,1,2,1,2);
    double J = Fp.giveDeterminant();
    FloatMatrix Fp0;
    Fp0.beSubMatrixOf(F0,1,2,1,2);
    double J0 = Fp0.giveDeterminant();
  // compute Fbar
    Fbar = Fp;
    Fbar.times(sqrt(J0/J));
    Fbar.resizeWithData(3,3);
    Fbar.at(3,3) = 1;  
  } else {
    J = F.giveDeterminant();
    Fbar = F;
    
    if(J < 0){
      nlElem->giveDomain()->giveEngngModel()->setAnalysisCrash(true);
    }
    J0 = F0.giveDeterminant();
    if (J0 < 0) {
      nlElem->giveDomain()->giveEngngModel()->setAnalysisCrash(true);
    }
    // compute Fbar
    if(J > 0 && J0 > 0) {
      Fbar.times(pow(J0/J,1./3.));
    }

    
  }
  // convert matrix to vector
  vFbar.beVectorForm(Fbar);
  // convert full vector form to reduced vector form
  StructuralMaterial :: giveReducedVectorForm(answer, vFbar, matMode);

  return J0/J;
  
}


void
FbarElementExtensionInterface :: computeFbarStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep, NLStructuralElement *elem)
{
    FloatMatrix B, B0, D1, D2, AB, DB, DB0, dPdF;
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();
    answer.clear();

    // create a fictitious integration point in the element centroid
    GaussPoint *gp0;
    FloatArray centroidCoords;
    MaterialMode matMode = elem->giveMaterialMode();
    elem->giveElementParametricCentroid(centroidCoords);
    gp0 = new GaussPoint( NULL, 1, centroidCoords, 1.0, matMode);
    // loop over gauss points
   for ( auto &gp : *elem->giveDefaultIntegrationRulePtr() ) {
     if ( elem->nlGeometry == 1 ) {
       elem->computeBHmatrixAt(gp, B, tStep,1);
       elem->computeBHmatrixAt(gp0,B0, tStep,1);
       if ( elem->domain->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
	   FloatMatrix D;
	   B0.subtract(B);
	   /// @todo We probably need overloaded function (like above) here as well.
	   cs->giveSpatialStiffnessMatrix(D, rMode, gp, tStep);
	   this->computeSpatialFbarStiffnessCorrections(D1,D,gp,tStep, elem);
	   double dV = elem->computeVolumeAround(gp);
	   DB.beProductOf(D,B);
	   DB0.beProductOf(D1, B0);
	   /*FloatMatrix DB1;
	   DB1.beProductOf(D1, B);
	   */
	   // Stiffness matrix for F-bar formulation is unsymmetric !!!
	   answer.plusProductUnsym(B, DB,  dV);
	   //	   answer.plusProductUnsym(B, DB1, -dV);
	   answer.plusProductUnsym(B, DB0, dV);

	   FloatMatrix answer1, answer2;
	   answer1.beTProductOf(B,B);
	   /*
	   answer1.clear();
	   answer2.clear();
	   answer1.plusProductUnsym(B, DB1, -dV);
	   answer2.plusProductUnsym(B, DB0, dV);
	   */

       } else { 
	   // Material stiffness dP/dF
	   cs->giveStiffnessMatrix_dPdF(dPdF, rMode, gp, tStep);
	   FloatMatrix fulldPdF;
	   if(gp->giveMaterialMode() == _3dMat) {
	     fulldPdF = dPdF;
	   } else if(gp->giveMaterialMode() == _PlaneStrain) {
	     fulldPdF.resize(9,9);
	     fulldPdF.at(1,1) = dPdF.at(1,1);
	     fulldPdF.at(1,2) = dPdF.at(1,2);
	     fulldPdF.at(1,3) = dPdF.at(1,3);
	     fulldPdF.at(1,6) = dPdF.at(1,4);
	     fulldPdF.at(1,9) = dPdF.at(1,5);
	     
	     fulldPdF.at(2,1) = dPdF.at(2,1);
	     fulldPdF.at(2,2) = dPdF.at(2,2);
	     fulldPdF.at(2,3) = dPdF.at(2,3);
	     fulldPdF.at(2,6) = dPdF.at(2,4);
	     fulldPdF.at(2,9) = dPdF.at(2,5);
	     
	     fulldPdF.at(3,1) = dPdF.at(3,1);
	     fulldPdF.at(3,2) = dPdF.at(3,2);
	     fulldPdF.at(3,3) = dPdF.at(3,3);
	     fulldPdF.at(3,6) = dPdF.at(3,4);
	     fulldPdF.at(3,9) = dPdF.at(3,5);
	     
	     fulldPdF.at(6,1) = dPdF.at(4,1);
	     fulldPdF.at(6,2) = dPdF.at(4,2);
	     fulldPdF.at(6,3) = dPdF.at(4,3);
	     fulldPdF.at(6,6) = dPdF.at(4,4);
	     fulldPdF.at(6,9) = dPdF.at(4,5);
	     
	     fulldPdF.at(9,1) = dPdF.at(5,1);
	     fulldPdF.at(9,2) = dPdF.at(5,2);
	     fulldPdF.at(9,3) = dPdF.at(5,3);
	     fulldPdF.at(9,6) = dPdF.at(5,4);
	     fulldPdF.at(9,9) = dPdF.at(5,5);
	   } else {
	     OOFEM_ERROR("Fbar interface can be used with large strain formulation only");
	   }

	   FloatArray invF;
	   double J0_J = this->computeFbarStiffnessCorrections(D1,D2,fulldPdF,gp,tStep, elem);
	   dPdF.times(pow(J0_J,-1./3.));
	   //	   dPdF.times(J0_J);
	   double dV = elem->computeVolumeAround(gp);
	   AB.beProductOf(dPdF,B);
	   DB.beProductOf(D1, B);
	   DB0.beProductOf(D2, B0);

	   // Stiffness matrix for F-bar formulation is unsymmetric !!!
	   answer.plusProductUnsym(B, AB, dV);
	   answer.plusProductUnsym(B, DB, dV);
	   answer.plusProductUnsym(B, DB0, dV);
   
	   
	 }
     }
   }
}


void 
FbarElementExtensionInterface :: computeSpatialFbarStiffnessCorrections(FloatMatrix &answer, FloatMatrix &a, GaussPoint *gp, TimeStep *tStep, NLStructuralElement *elem)
{
  
  StructuralMaterialStatus *status = static_cast<   StructuralMaterialStatus * >( elem->giveMaterial()->giveStatus(gp) );

  FloatArray vC(status->giveTempCVector());
  if(!vC.giveSize())
    vC.resize(6);

  answer.resize(9,9);
  /// axismetric case need to be extended to full 3d case
  answer.at(1,1) = (a.at(1,1) + a.at(1,2) + a.at(1,3))/3. - 2./3. *  vC.at(1);  
  answer.at(1,2) = (a.at(1,1) + a.at(1,2) + a.at(1,3))/3. - 2./3. *  vC.at(1);  
  answer.at(1,3) = (a.at(1,1) + a.at(1,2) + a.at(1,3))/3. - 2./3. *  vC.at(1);  

  answer.at(2,1) = (a.at(2,1) + a.at(2,2) + a.at(2,3))/3. - 2./3. *  vC.at(2);  
  answer.at(2,2) = (a.at(2,1) + a.at(2,2) + a.at(2,3))/3. - 2./3. *  vC.at(2);  
  answer.at(2,3) = (a.at(2,1) + a.at(2,2) + a.at(2,3))/3. - 2./3. *  vC.at(2);  

  answer.at(3,1) = (a.at(3,1) + a.at(3,2) + a.at(3,3))/3. - 2./3. *  vC.at(3);  
  answer.at(3,2) = (a.at(3,1) + a.at(3,2) + a.at(3,3))/3. - 2./3. *  vC.at(3);  
  answer.at(3,3) = (a.at(3,1) + a.at(3,2) + a.at(3,3))/3. - 2./3. *  vC.at(3);  

  answer.at(6,1) = (a.at(6,1) + a.at(6,2) + a.at(6,3))/3. - 2./3. *  vC.at(6);  
  answer.at(6,2) = (a.at(6,1) + a.at(6,2) + a.at(6,3))/3. - 2./3. *  vC.at(6);  
  answer.at(6,3) = (a.at(6,1) + a.at(6,2) + a.at(6,3))/3. - 2./3. *  vC.at(6);  

  answer.at(9,1) = (a.at(9,1) + a.at(9,2) + a.at(9,3))/3. - 2./3. *  vC.at(6);  
  answer.at(9,2) = (a.at(9,1) + a.at(9,2) + a.at(9,3))/3. - 2./3. *  vC.at(6);  
  answer.at(9,3) = (a.at(9,1) + a.at(9,2) + a.at(9,3))/3. - 2./3. *  vC.at(6);
  
}
  

double
FbarElementExtensionInterface :: computeFbarStiffnessCorrections(FloatMatrix &K1, FloatMatrix &K2, const FloatMatrix &dPdF, GaussPoint *gp, TimeStep *tStep, NLStructuralElement *elem)
{
  
  // create a fictitious integration point in the element centroid
  StructuralMaterialStatus *status = static_cast<   StructuralMaterialStatus * >( elem->giveMaterial()->giveStatus(gp) );

  FloatArray vP(status->giveTempPVector());
  if(!vP.giveSize())
    vP.resize(9);



  GaussPoint *gp0;
  FloatArray centroidCoords;
  MaterialMode matMode = gp->giveMaterialMode();

  elem->giveElementParametricCentroid(centroidCoords);
  gp0 = new GaussPoint( NULL, 1, centroidCoords, 1.0, matMode);
  // compute compatible deformation gradient
  FloatArray reducedvF, vF, u;
  FloatMatrix F, invF;
  elem->computeVectorOf({D_u, D_v, D_w}, VM_Total, tStep, u); // solution vector
  this->computeDeformationGradientVector(reducedvF,gp,tStep, elem);
  StructuralMaterial :: giveFullVectorFormF(vF, reducedvF, matMode); // 9 components
  F.beMatrixForm(vF);  
  invF.beInverseOf(F);
  // compute deformation gradient in the element centroid
  FloatArray reducedvF0,vF0;
  FloatMatrix invF0, F0;
  this->computeDeformationGradientVector(reducedvF0, gp0, tStep, elem);
  StructuralMaterial :: giveFullVectorFormF(vF0, reducedvF0, matMode); // 9 components
  F0.beMatrixForm(vF0);
  invF0.beInverseOf(F0);

  double J = F.giveDeterminant();
  double J0 = F0.giveDeterminant();
  double J0_J = J0/J;
  double factor = pow(J0_J,-1./3.);



  FloatArray AF;
  AF.beProductOf(dPdF,vF);
  

  K1.resize(9,9);
  K2.resize(9, 9);

  
  for ( int i = 1; i <= 3; i++ ) {
    for ( int j = 1; j <= 3; j++ ) {
      for ( int k = 1; k <= 3; k++ ) {
	for ( int l = 1; l <= 3; l++ ) {
	  int IJ = StructuralMaterial :: giveVI(i,j);
	  int KL = StructuralMaterial :: giveVI(k,l);
	   K1.at( IJ, KL ) = 2./3. * factor * factor * vP.at(IJ) * invF.at(l,k)-1./3.*factor*AF.at(IJ) * invF.at(l,k);;
	   K2.at( IJ, KL ) = - 2./3. * factor * factor * vP.at(IJ) * invF0.at(l,k)+ 1./3.*factor*AF.at(IJ) * invF0.at(l,k);;
	}
      }
    }
  }

  if( matMode == _3dMat) {
  } else if (matMode == _PlaneStrain) {
    FloatMatrix K13, K23;

    K13 = K1;
    K23 = K2;

    K1.resize(5,5);
    K2.resize(5,5);

    K1.at(1, 1) = K13.at(1, 1);
    K1.at(1, 2) = K13.at(1, 2);    
    K1.at(1, 4) = K13.at(1, 6);
    K1.at(1, 5) = K13.at(1, 9);

    K1.at(2, 1) = K13.at(2, 1);
    K1.at(2, 2) = K13.at(2, 2);
    K1.at(2, 4) = K13.at(2, 6);
    K1.at(2, 5) = K13.at(2, 9);

    K1.at(3, 1) = K13.at(3, 1);
    K1.at(3, 2) = K13.at(3, 2);
    K1.at(3, 4) = K13.at(3, 6);
    K1.at(3, 5) = K13.at(3, 9);

    K1.at(4, 1) = K13.at(6, 1);
    K1.at(4, 2) = K13.at(6, 2);
    K1.at(4, 4) = K13.at(6, 6);
    K1.at(4, 5) = K13.at(6, 9);

    K1.at(5, 1) = K13.at(9, 1);
    K1.at(5, 2) = K13.at(9, 2);
    K1.at(5, 5) = K13.at(9, 9);
    K1.at(5, 5) = K13.at(9, 9);


    K2.at(1, 1) = K23.at(1, 1);
    K2.at(1, 2) = K23.at(1, 2);    
    K2.at(1, 4) = K23.at(1, 6);
    K2.at(1, 5) = K23.at(1, 9);

    K2.at(2, 1) = K23.at(2, 1);
    K2.at(2, 2) = K23.at(2, 2);
    K2.at(2, 4) = K23.at(2, 6);
    K2.at(2, 5) = K23.at(2, 9);

    K2.at(3, 1) = K23.at(3, 1);
    K2.at(3, 2) = K23.at(3, 2);
    K2.at(3, 4) = K23.at(3, 6);
    K2.at(3, 5) = K23.at(3, 9);

    K2.at(4, 1) = K23.at(6, 1);
    K2.at(4, 2) = K23.at(6, 2);
    K2.at(4, 4) = K23.at(6, 6);
    K2.at(4, 5) = K23.at(6, 9);

    K2.at(5, 1) = K23.at(9, 1);
    K2.at(5, 2) = K23.at(9, 2);
    K2.at(5, 5) = K23.at(9, 9);
    K2.at(5, 5) = K23.at(9, 9);


  } else {
    OOFEM_ERROR("F-bar interface supports only _3dMat and _PlaneStrain modes");
  }

  

return J0_J;

}


} // end namespace oofem
