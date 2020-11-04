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

#include "../sm/Materials/HyperelasticMaterials/ElectroMechanics/mooneyrivlin_idealdielectric_transverselyisotropic.h"
#include "../sm/Materials/HyperelasticMaterials/mooneyrivlin.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"




namespace oofem {
  REGISTER_Material(MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial);

  MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial(int n, Domain *d):Material(n, d), ElectroMechanicalMaterialExtensionInterface_3Field(d)
{
    this->hyperelasticMaterial = new MooneyRivlinMaterial(n, d);
}


IRResultType
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro


    result = this->hyperelasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) {
      return result;
    }

    IR_GIVE_FIELD(ir, this->epsilon_m, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial_epsilon);
    IR_GIVE_FIELD(ir, this->epsilon_f, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial_epsilonFibre);
    IR_GIVE_FIELD(ir, this->N, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial_director);
      
    return IRRT_OK;
}



void
MooneyRivlin_IdealDielectricMaterial :: give_FirstPKStressVector_ElectricalField_3d(FloatArray &vP, FloatArray &E, GaussPoint *gp, const FloatArray &vF, const FloatArray &D, TimeStep *tStep)
{
    // First Piol-Kirchhoff stress
    double Sigma_J;
    FloatArray vSigma_H;
    FloatMatrix F, Sigma_H;
    hyperelasticMaterial->giveFirstPKStressVector_3d(vP, gp, vF, tStep);    
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );   
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    //
    FloatMatrix H;
    H.beMatrixForm(vH);
    //
    this->compute_I7(I7, F);
    this->compute_I10(I10, F);
    //derivatives wrt F
    this->compute_dI4_dF(dI4_dF, F);
    this->compute_dI7_dF(dI7_dF, F);
    this->compute_dI10_dF(dI10_dF, F);
    //derivatives wrt D_0
    this->compute_dI7_dD(dI7_dD, D);
    this->compute_dI10_dD(dI10_dD, D);
    //
    vP.add(dI4, Cm);
    vP.add(dI7, 1./this->eps_f / J);
    vP.add(cofF, - I7 / J / J / this->eps_m);
    vP.add(dI10, 2. * I10 / /this->eps_f / J);
    vP.add(cofF, - I10 / J / J / this->eps_f);

    E.add(dI7_dD, 1. / J / this->eps_m);
    E.add(dI10_dD, 2. * I10 / J / this->eps_f);

    
    
    status->letTempPVectorBe(vP);
    status->letTempFVectorBe(vF);
    status->letTempEVectorBe(E);
    status->letTempDVectorBe(D);
}


  
  

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
 
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    hyperelasticMaterial->give3dMaterialStiffnessMatrix_dPdF(answer, mode, gp, tStep); 
    FloatArray vF = status->giveTempFVector();
    FloatArray E = status->giveTempEVector();

    // First Piola-Kirchhoff stress
    FloatMatrix F;
    // Kronecker delta
    FloatMatrix delta(3,3);
    delta.beUnitMatrix();
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    //
    FloatMatrix H;
    H.beMatrixForm(vH);
    ////////////////////////////////////    
    FloatMatrix Fx;
    F.beMatrixForm(vF);
    hyperelasticMaterial->compute_tensor_cross_product_tensor( Fx, F );
    ///////////////////////////////////
    double I7, I10;
    this->compute_I7(I7, F);
    this->compute_I7(I10, F);
    FloatMatrix d2I4_dF2, dI7_dF, d2I7_dF2, dI10_dF;
    this->compute_d2I4_dF2(d2I4_dF2, F);
    this->compute_dI7_dF(dI7_dF, F);
    this->compute_d2I7_dF2(d2I7_dF2, F);
    this->compute_dI10_dF(dI10_dF, F);
    FloatMatrix cFcF, cFdI7, dI7cF, Fx, dI10dI10, dI10cF, cFdI10;
    cFcF.beDyadicProductOf(cofF, cofF);
    
    answer.add(d2I4_dF2, this->Cm);
    answer.add(d2I7_dF2, 1. / this->epsilon_m / J);
    answer.add(cFcF,( I7/this->epsilon_m + I10 * I10 /this->epsilon_f) / J / J / J);
    answer.add(Fx, - ( I7/this->epsilon_m + I10 * I10 / this->epsilon_f) / J / J );
    answer.add(dI10dI10, 2 / this->epsilon_f / J);
    answer.add(cFI7s, - 1 / this->eps_m / J / J);
    answer.add(cFI10s, - 2 * I10 / this->epsilon_f / J / J);  

    
    
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: give3dMaterialStiffnessMatrix_dPdD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

      ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    hyperelasticMaterial->give3dMaterialStiffnessMatrix_dPdF(answer, mode, gp, tStep); 
    FloatArray vF = status->giveTempFVector();
    FloatArray E = status->giveTempEVector();

    // First Piola-Kirchhoff stress
    FloatMatrix F;
    // Kronecker delta
    FloatMatrix delta(3,3);
    delta.beUnitMatrix();
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    //
    FloatMatrix H;
    H.beMatrixForm(vH);
    ////////////////////////////////////    
    FloatMatrix Fx;
    F.beMatrixForm(vF);
    hyperelasticMaterial->compute_tensor_cross_product_tensor( Fx, F );
    ///////////////////////////////////
    double I10;
    this->compute_I10(I10, F);
    FloatMatrix d2I4_dF2, dI7_dF, d2I7_dF2, dI10_dF;
    this->compute_d2I7_dFdD(d2I4_dF2, F, D);
    this->compute_d2I10_dFdD(dI7_dF, F, D);
    this->compute_dI7_dD(dI10_dF, F);
    this->compute_dI10_dF(dI10_dF, F);
    this->compute_dI10_dD(dI10_dD, D);
    FloatMatrix cFcF, cFdI7, dI7cF, dI10cF, cFdI10;

    
    answer.add(d2I7_dFdD, 1. / this->epsilon_m / J);
    answer.add(d2I10_dFdD, 2. * I10 / this->epsilon_f / J);
    answer.add(dI10dF_dI10dDF, 2 / J / this->epsilon_f);
    answer.add(cF_dI7dD, - 1. / J / J / this->epsilon_m);
    answer.add(cF_dI10dD, - 2. * I10 / J / J / this->epsilon_f);
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: give3dMaterialStiffnessMatrix_dEdD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    FloatArray vF = status->giveTempFVector();
    FloatArray D = status->giveTempDVector();
    FloatMatrix F;
    F.beMatrixForm(vF);

    FloatMatrix d2I7_dD2, dI10_dD, d2I10_dD2;
    this->compute_d2I7_dD2(d2I7_dD2, D);
    this->compute_dI10_dD(dI10_dD, D);
    dI10dI10.beDyadicProductOf(dI10, dI10);
    double J = F.giveDeterminant(F);   
    answer.add(d2I7_dD2, 1./J / this->epsilon_m);
    answer.add(dI10dI10, 2 * I10 / J / this->epsilon_f);
}



double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: computeI4(const FloatMatrix &F)
{
  
  FloatArray CN;
  FloatMatrix C;
  C.beTProductOf(F,F);
  CN.beProductOf(C,this->N);
  return this->N.dotProduct(CN);  
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: computeI7(const FloatMatrix &F, const FloatArray &D)
{
  
  FloatArray CD;
  FloatMatrix C;
  C.beTProductOf(F,F);
  CD.beProductOf(C,D);
  return D.dotProduct(CD);  
}


double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: computeI10(const FloatMatrix &F, const FloatArray &D_0)
{
  
  FloatArray CD;
  FloatMatrix C;
  C.beTProductOf(F,F);
  CN.beProductOf(C,D);
  return this->N.dotProduct(CD);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI4dF(FloatMatrix &answer, const FloatMatrix &F)
{
  
  FloatArray FN;
  FN.beProductOf(F,N)

  return this->N.dotProduct(CN);  
}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI7dF(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  
  FloatArray FN;
  FN.beProductOf(F,N)

  return this->N.dotProduct(CN);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI7dD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  
  FloatArray FN;
  FN.beProductOf(F,N)

  return this->N.dotProduct(CN);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI10dF(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  
  FloatArray FN;
  FN.beProductOf(F,N)

  return this->N.dotProduct(CN);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI10dD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  
  FloatArray FN;
  FN.beProductOf(F,N)

  return this->N.dotProduct(CN);  
}

///
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I4dF2(FloatMatrix &answer, const FloatMatrix &F)
{
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {
	  
	  answer.at(giveVI(i,j), giveVI(k,l)) = delta.at(i,k) * N.at(j) * N.at(l) +  N.at(i) * N.at(l) * delta.at(j,k);
	}
      }
    }
  }   

}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I7dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {  
	  answer.at(giveVI(i,j), giveVI(k,l)) = delta.at(i,k) * D.at(j) * D.at(l) +  D.at(i) * D.at(l) * delta.at(j,k);
	}
      }
    }
  }
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I10dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{


  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {  
	  answer.at(giveVI(i,j), giveVI(k,l)) = delta.at(i,k) * D.at(j) * N.at(l) +  N.at(i) * D.at(l) * delta.at(j,k);
	}
      }
    }
  }

}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I7dD2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  answer.beTProductOf(F,F);
  answer.times(2);
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I7dFdF(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	answer.at(giveVI(i,j), k) = 2. * delta.at(i,k) * D.at(j) +  F.at(j,i) * D.at(k);
      }
    }
  }  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I7dFdF(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	answer.at(giveVI(i,j), k) = 2. * delta.at(i,k) * N.at(j) +  F.at(j,i) * N.at(k);
      }
    }
  }  
}





int
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

  //  ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    
   
    return 1;
    
}
    
  


} // end namespace oofem
