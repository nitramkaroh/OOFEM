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

#include "../sm/Materials/Micromorphic/ContactGradientMaterial/airmaterialcontactrot.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
  REGISTER_Material(AirMaterialContactRot);

AirMaterialContactRot :: AirMaterialContactRot(int n, Domain *d) : AirMaterial(n, d), SecondGradientMaterialExtensionInterface(d)
{
  Ak = 0.;
}

AirMaterialContactRot :: ~AirMaterialContactRot()
{ }



void
AirMaterialContactRot :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}



void
AirMaterialContactRot :: giveSecondGradientMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    AirMaterialContactRotStatus *status = static_cast< AirMaterialContactRotStatus * >( this->giveStatus(gp) );
    FloatMatrix m3d;
    AirMaterial :: givePlaneStrainStiffMtrx_dPdF(answer, mode,gp, tStep);
}


void
AirMaterialContactRot :: giveSecondGradientMatrix_dPdG(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  answer.clear();
}



void
AirMaterialContactRot :: giveSecondGradientMatrix_dMdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    AirMaterialContactRotStatus *status = static_cast< AirMaterialContactRotStatus* >( this->giveStatus(gp) );
    answer.resize(8,5);
    FloatArray vF, vG, vCofF;
    FloatMatrix F, cofF, invF;
    vG = status->giveTempGVector();
    F.beMatrixForm(vF = status->giveTempFVector());
    invF.beInverseOf(F);
    double J = F.giveDeterminant();
    cofF.beTranspositionOf(invF);
    cofF.times(J);
    vCofF.beVectorForm(cofF);
    FloatArray redCofF;
    StructuralMaterial :: giveReducedVectorForm(redCofF, vCofF, _PlaneStrain);
    //
    answer.beProductTOf(vG, redCofF);
    answer.times(-a * Ak * exp(-a * J));        

}

void
AirMaterialContactRot :: giveSecondGradientMatrix_dMdG(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  AirMaterialContactRotStatus *status = static_cast< AirMaterialContactRotStatus* >( this->giveStatus(gp) );
  FloatMatrix F;
  F.beMatrixForm(status->giveTempFVector());
  double J = F.giveDeterminant(); 
  answer.resize(8,8);
  answer.beUnitMatrix();
  answer.times( Ak * exp(-a * J));        
}  




void
AirMaterialContactRot :: giveFiniteStrainGeneralizedStressVectors(FloatArray &vP, FloatArray &vM, const FloatArray &vF, const FloatArray &vG, GaussPoint *gp, TimeStep *tStep)
{
    ///@todo Move this to StructuralCrossSection ?
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
      this->giveFiniteStrainGeneralizedStressVectors_3d(vP, vM, vF, vG, gp, tStep);
    } else if ( mode == _PlaneStrain ) {
      this->giveFiniteStrainGeneralizedStressVectors_PlaneStrain(vP, vM, vF, vG, gp, tStep);
    } else {
      OOFEM_ERROR("Unknown material mode for the gradient polyconvex formulation");
    }
}
  




void
AirMaterialContactRot :: giveFiniteStrainGeneralizedStressVectors_3d(FloatArray &vP, FloatArray &vM, const FloatArray &vF, const FloatArray &vG, GaussPoint *gp, TimeStep *tStep)
{
  
  OOFEM_ERROR("3d not implemented");
      
}

  


void
AirMaterialContactRot :: giveFiniteStrainGeneralizedStressVectors_PlaneStrain(FloatArray &vP, FloatArray &vM, const FloatArray &reducedvF, const FloatArray &reducedvG, GaussPoint *gp, TimeStep *tStep)
{ 

  AirMaterialContactRotStatus *status = static_cast< AirMaterialContactRotStatus* >( this->giveStatus(gp) );
  FloatArray gradRn, gradFn, Fn;
  /*  status->give_lnRn(gradRn);
  status->give_gradF(gradFn);
  status->giveF(Fn);

  //
  FloatArray f(reducedvF);
  f.subtract(F);
  //
  FloatArray gradF(reducedvG);
  //
  FloatArray grad_f(reducedvG);
  f.subtract(gradFn);
  ///
  FloatArray grad_w(reducedvG);
  grad_w.add(gradR);
  //
  this->compute_d_gradw_d_F(dwdF, f, grad_f, F, gradF);
  */
  //deformation gradient, its inverse, cofactor, and determinant
  FloatArray vCofF, vF;
  FloatMatrix F, invF, cofF;
  StructuralMaterial :: giveFullVectorForm( vF, reducedvF, gp->giveMaterialMode() );
  F.beMatrixForm(vF);
  invF.beInverseOf(F);
  double J = F.giveDeterminant();
  cofF.beTranspositionOf(invF);
  cofF.times(J);
  vCofF.beVectorForm(cofF);
  FloatArray redCofF;
  StructuralMaterial :: giveReducedVectorForm(redCofF, vCofF, _PlaneStrain);
  // higher order stress
  vM.resize(8);
  vM.zero();
  vM.add(Ak*exp(-a*J), reducedvG);
  // first PK stress
  AirMaterial :: giveFirstPKStressVector_PlaneStrain(vP, gp, reducedvF, tStep);
  FloatArray fullvP;
  StructuralMaterial :: giveFullVectorForm( fullvP, vP, gp->giveMaterialMode() );
    
    
  status->letTempPVectorBe(fullvP);
  status->letTempFVectorBe(vF);
  //
  status->letTempMVectorBe(vM);
  status->letTempGVectorBe(reducedvG);


  
      
}
  

  
  

IRResultType
AirMaterialContactRot :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    AirMaterial :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, Ak, _IFT_SecondGradientMaterialExtensionInterface_Ak);
    IR_GIVE_FIELD(ir, a, _IFT_AirMaterialContactRot_a);
   
    
    return IRRT_OK;
}

int
AirMaterialContactRot :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    AirMaterialContactRotStatus *status = static_cast< AirMaterialContactRotStatus * >( this->giveStatus(gp) );
    

    if( type == IST_MicromorphicStress) {
      //StructuralMaterial :: giveFullVectorForm( answer, status->giveStressVector(), gp->giveMaterialMode() );

    } else if( type == IST_MicromorphicStrain ) {
      //      StructuralMaterial :: giveFullVectorForm( answer, status->giveStrainVector(), gp->giveMaterialMode() );
    }  else if ( type == IST_MicromorphicHigherOrderStress ) {
      /*FloatArray M = status->giveMicromorphicStressGrad(); 
      answer.resize(9);
      answer.zero();
      answer.at(6) = M.at(1);
      answer.at(9) = M.at(2);
      */
    } else if ( type == IST_MicromorphicHigherOrderStrain ) {
      /*
      FloatArray kappa = status->giveMicromorphicVarGrad();
      answer.resize(9);
      answer.zero();
      answer.at(6) = kappa.at(1);
      answer.at(9) = kappa.at(2);
      */
    } else if( type == IST_DeformationGradientTensor) {
      answer = status->giveFVector();
    } else {
      AirMaterial::giveIPValue(answer, gp, type, tStep);
    }
    return 1;
}
    

  AirMaterialContactRotStatus :: AirMaterialContactRotStatus(int n, Domain *d, GaussPoint *g) : SecondGradientMaterialStatus(n, d, g)
{
  // 2d plane strain
  /*    micromorphicVar.resize(5);
    micromorphicVar.at(1) = micromorphicVar.at(2) = micromorphicVar.at(3) = 1.;
    micromorphicStress.resize(5);
  */
    //2d plane strain
  /* micromorphicVarGrad.resize(8);
    micromorphicStressGrad.resize(8);
  */
  /*
    tempMicromorphicVar = micromorphicVar;
    tempMicromorphicVarGrad = micromorphicVarGrad;
    tempMicromorphicStress = micromorphicStress;
    tempMicromorphicStressGrad = micromorphicStressGrad;
  */
}


} // end namespace oofem
