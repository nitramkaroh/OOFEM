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

#include "../sm/Materials/Micromorphic/Microdil/straindivergencematerial_elastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "strainvector.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"

#include "dynamicinputrecord.h"



namespace oofem {
  REGISTER_Material(StrainDivergenceMaterial_Elastic);

StrainDivergenceMaterial_Elastic :: StrainDivergenceMaterial_Elastic(int n, Domain *d) :IsotropicLinearElasticMaterial(n, d), SwcondGradientMaterialExtensionInterface(d)
{
  Hk = Ak = 0.;
  equivStrainType = EST_MicroDil;

}

StrainDivergenceMaterial_Elastic :: ~StrainDivergenceMaterial_Elastic()
{ }



void
StrainDivergenceMaterial_Elastic :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}



void
StrainDivergenceMaterial_Elastic :: giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();



  if(matMode == _3dMat){
    IsotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
  } else if(matMode == _PlaneStrain) {
    IsotropicLinearElasticMaterial :: givePlaneStrainStiffMtrx(answer, mode, gp, tStep);   
  } else if(matMode == _PlaneStress) {
    IsotropicLinearElasticMaterial :: givePlaneStressStiffMtrx(answer, mode, gp, tStep);   
  }


}

void
StrainDivergenceMaterial_Elastic :: giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
 MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );    
  FloatArray tempStrain, dEeqdE;
  FloatMatrix Iv, d2EqdE2;
  tempStrain = status->giveTempStrainVector();
   
  MaterialMode matMode = gp->giveMaterialMode();
  this->compute_dEeqdE(dEeqdE, tempStrain, matMode);

  if(matMode == _3dMat) {
    answer = {{dEeqdE.at(1), dEeqdE.at(2), dEeqdE.at(3), dEeqdE.at(4), dEeqdE.at(5), dEeqdE.at(6)}};
    //answer = {{1, 1, 1, 0, 0, 0}};
  } else if(matMode == _PlaneStrain) {
    //answer = {{1, 1, 1, 0}};
    answer = {{dEeqdE.at(1), dEeqdE.at(2), dEeqdE.at(3), dEeqdE.at(4)}};
  }

  answer.times(-Hk);

}



{
  MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );    
  FloatArray tempStrain, tempMicromorphicVar, dEeqdE;
  FloatMatrix Iv, d2EeqdE2;

  tempStrain = status->giveTempStrainVector();
  tempMicromorphicVar = status->giveTempMicromorphicVar();
 // @todo resizovat ve statusus
  if( tempMicromorphicVar.giveSize() != 1)
    tempMicromorphicVar.resize(1);

 double eEq = this->computeEquivalentStrain(tempStrain, matMode);
  this->compute_dEeqdE(dEeqdE, tempStrain, matMode);
  this->compute_d2EeqdE2(d2EeqdE2, tempStrain, matMode);


  /*
  if(matMode == _3dMat){
      I = {1, 1, 1, 0, 0, 0};
  } else if(matMode == _PlaneStrain) {
      I = {1, 1, 1, 0};
  } else if(matMode == _PlaneStress) {
      I = {1, 1, 0};
  }
  */

  // Iv.beDyadicProductOf(dEeq_dE,dEeq_dE);
  Iv.beDyadicProductOf(dEeqdE,dEeqdE);
  //Iv.beDyadicProductOf(I,I);
  Iv.times(Hk);

   double e = eEq - tempMicromorphicVar.at(1);
  d2EeqdE2.times(Hk * e);



  
  answer.add(Iv);
  answer.add(d2EeqdE2);


}


void
StrainDivergenceMaterial_Elastic :: giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
 FloatArray tempStrain, dEeqdE;

  MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );    

  tempStrain = status->giveTempStrainVector();
   
  MaterialMode matMode = gp->giveMaterialMode();
  this->compute_dEeqdE(dEeqdE, tempStrain, matMode);

  if(matMode == _3dMat) {
    //answer = {{1},{1},{1},{0},{0},{0}};
    answer = {{dEeqdE.at(1)},{dEeqdE.at(2)},{dEeqdE.at(3)},{dEeqdE.at(4)},{dEeqdE.at(5)},{dEeqdE.at(6)}};
  } else if(matMode == _PlaneStrain) {
    answer = {{dEeqdE.at(1)},{dEeqdE.at(2)},{dEeqdE.at(3)},{dEeqdE.at(4)}};
    //answer = {{1},{1},{1},{0}};
  }  
  answer.times(-Hk);

}



void
StrainDivergenceMaterial_Elastic :: giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _3dMat){
    answer.resize(3,3);
    answer.at(1,1) = Hk;
    answer.at(2,2) = Hk;
    answer.at(3,3) = Hk;
  } else if(matMode == _PlaneStrain) {
    answer.resize(1,1);
    answer.at(1,1) = Hk;
    //    answer.at(2,2) = Hk;
  }

}


void
StrainDivergenceMaterial_Elastic :: giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _3dMat){
    answer.resize(3,3);
    answer.at(1,1) = Ak;
    answer.at(2,2) = Ak;
    answer.at(3,3) = Ak;
  } else if(matMode == _PlaneStrain) {
    answer.resize(2,2);
    answer.at(1,1) = Ak;
    answer.at(2,2) = Ak;
  }

}




void
StrainDivergenceMaterial_Elastic :: giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &S, GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{
  
    MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );
    double eV;
    FloatMatrix De;
    MaterialMode matMode = gp->giveMaterialMode();
    FloatArray I;
    if(matMode == _3dMat){
      IsotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix(De, TangentStiffness, gp, tStep);
    
    } else if(matMode == _PlaneStrain) {
      IsotropicLinearElasticMaterial :: givePlaneStrainStiffMtrx(De, TangentStiffness, gp, tStep);   
    } 

    eV = this->computeEquivalentStrain(totalStrain,matMode);
    this->compute_dEeqdE(I, totalStrain, matMode);
    
    sigma.beProductOf(De,totalStrain);
    s = {micromorphicVar.at(1)-eV};
    s.times(Hk);
    I.times(-s.at(1));
    sigma.add(I);

    S = micromorphicVarGrad;
    S.times(Ak);

    status->letTempStrainVectorBe(totalStrain);
    status->letTempMicromorphicVarBe(micromorphicVar);
    status->letTempMicromorphicVarGradBe(micromorphicVarGrad);

    status->letTempStressVectorBe(sigma);
    status->letTempMicromorphicStressBe(s);
    status->letTempMicromorphicStressGradBe(S); 
      
}



double
StrainDivergenceMaterial_Elastic :: computeEquivalentStrain(const FloatArray &strain, MaterialMode matMode)
{
  
  if(matMode == _3dMat || matMode == _PlaneStrain){
    if(this->equivStrainType == EST_MicroDil){
      return (strain.at(1)+strain.at(2)+strain.at(3));
    } else if(this->equivStrainType == EST_MicroStrain) {
	StrainVector strainVec(strain,matMode);
	return strainVec.computeStrainNorm();
    } else
      return 0.;
  } else { 
  return 0.;
  }

}


void
StrainDivergenceMaterial_Elastic :: compute_dEeqdE(FloatArray &answer, const FloatArray &strain, MaterialMode matMode)
{
   if(this->equivStrainType == EST_MicroDil){
      if(matMode == _3dMat){
	answer = {1, 1, 1, 0, 0, 0};
      } else if(matMode == _PlaneStrain) {
	answer = {1, 1, 1, 0};
      }
    } else if(this->equivStrainType == EST_MicroStrain) {
	double norm = computeEquivalentStrain(strain, matMode);
	answer = strain;
	if(matMode == _3dMat){
	  answer.at(4) = answer.at(4)/2.; 
	  answer.at(5) = answer.at(5)/2.; 
	  answer.at(6) = answer.at(6)/2.;
      } else if(matMode == _PlaneStrain) {
	  answer.at(4) = answer.at(4)/2.; 
      }
	if(norm == 0)
	  answer.zero();
	else
	  answer.times(1./norm);
      }
    
  
}

void
StrainDivergenceMaterial_Elastic :: compute_d2EeqdE2(FloatMatrix &answer, const FloatArray &strain, MaterialMode matMode)
{
   if(this->equivStrainType == EST_MicroDil){
      answer.resize(0,0);
    } else if(this->equivStrainType == EST_MicroStrain) {
	double norm = this ->computeEquivalentStrain(strain,matMode);
	FloatArray Pe;
	FloatMatrix P, PePe;
	if(matMode == _3dMat){
	  P.resize(6,6);
	  P.at(1,1) = P.at(2,2) = P.at(3,3) = 1.;
	  P.at(4,4) = P.at(5,5) = P.at(6,6) = 0.5;
      } else if(matMode == _PlaneStrain) {
	  P.resize(4,4);	  
	  P.at(1,1) = P.at(2,2) = P.at(3,3) = 1.;
	  P.at(4,4) = 0.5;
      }
	
	answer = P;
	
	Pe.beProductOf(P,strain);
	PePe.beDyadicProductOf(Pe,Pe);
	PePe.times(-1./norm/norm);

	answer.add(PePe);
	if(norm == 0)
	  answer.zero();
	else
	  answer.times(1./norm);
      }
  
  
}



void
StrainDivergenceMaterial_Elastic :: giveInputRecord(DynamicInputRecord &input)
{
    IsotropicLinearElasticMaterial::giveInputRecord(input);

    input.setField(this->Hk, _IFT_MicromorphicMaterialExtensionInterface_Hk);
    input.setField(this->Ak, _IFT_MicromorphicMaterialExtensionInterface_Ak);
    //    input.setField(this->propertyDictionary.at(tAlpha), _IFT_IsotropicLinearElasticMaterial_talpha);
}


IRResultType
StrainDivergenceMaterial_Elastic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    IsotropicLinearElasticMaterial :: initializeFrom(ir);
    int epsEq = 0;

    IR_GIVE_OPTIONAL_FIELD(ir, epsEq, _IFT_StrainDivergenceMaterial_Elastic_equivstraintype);
    IR_GIVE_FIELD(ir, Hk, _IFT_MicromorphicMaterialExtensionInterface_Hk);
    IR_GIVE_FIELD(ir, Ak, _IFT_MicromorphicMaterialExtensionInterface_Ak);

    //    propertyDictionary.add(Hk, Hk);
    //propertyDictionary.add(Ak, Ak);


    return IRRT_OK;
}

int
StrainDivergenceMaterial_Elastic :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    SecondGradientcMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );
    
    if ( type == IST_MaxEquivalentStrainLevel ) {
        FloatArray mV = status->giveTempMicromorphicVar();
        answer.resize(1);
	answer.at(1) = mV.at(1);
        return 1;
    } else if( type == IST_MicromorphicRelativeStrain) {
      FloatArray strain = status->giveStrainVector();
      FloatArray strainDev;
      double eV = computeDeviatoricVolumetricSplit(strainDev, strain);
      FloatArray mV = status->giveTempMicromorphicVar();
      
      answer.resize(9);
      answer.at(1) = eV - mV.at(1);
      answer.zero();
      return 1;

    } else if (type == IST_MicromorphicRelativeStress) {
      FloatArray s = status->giveMicromorphicStress();
      answer.resize(9);
      answer.zero();
      answer.at(1) = s.at(1);
      return 1;
    } else if (type == IST_MicromorphicHigherOrderStress) {
      answer.resize(9);
      answer.zero();
      FloatArray stress = status->giveMicromorphicStressGrad();
      answer.at(1) = stress.at(1);
      answer.at(4) = stress.at(2);
      return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }    
}
    



} // end namespace oofem
