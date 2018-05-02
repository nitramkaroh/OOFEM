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

#include "../sm/Materials/Micromorphic/Micropolar/micropolarmaterial_elastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
  REGISTER_Material(MicropolarMaterial_Elastic);

MicropolarMaterial_Elastic :: MicropolarMaterial_Elastic(int n, Domain *d) :IsotropicLinearElasticMaterial(n, d), MicromorphicMaterialExtensionInterface(d)
{
  Hk = Ak = 0.;
}

MicropolarMaterial_Elastic :: ~MicropolarMaterial_Elastic()
{ }



void
MicropolarMaterial_Elastic :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}



void
MicropolarMaterial_Elastic :: giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  /*  
MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );  


  double perturbation = 1.e-8;
  
  FloatArray uGrad, mV, mVG;
  uGrad.resize(5);
  mV.resize(1);
  mVG.resize(1);  
  uGrad.zero();
  mV.zero();
  mVG.zero();

  mV = status->giveMicromorphicVar();
  mVG = status->giveMicromorphicVarGrad();
  uGrad = status->giveStrainVector();
  uGrad.resize(5);
  
  FloatArray uGp(5);
  uGp = uGrad;
  FloatArray sigma, s, S;
  FloatMatrix stiff(5,5);
  stiff.zero();

  uGp.at(1) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGp, mV, mVG, tStep);

  for(int i = 1; i <= 5; i++) {
    stiff.at(i,1) = sigma.at(i);    
  }

  uGp = uGrad;
  uGp.at(2) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGp, mV, mVG, tStep);

  for(int i = 1; i <= 5; i++) {
    stiff.at(i,2) = sigma.at(i);    
  }


  uGp = uGrad;
  uGp.at(3) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGp, mV, mVG, tStep);

  for(int i = 1; i <= 5; i++) {
    stiff.at(i,3) = sigma.at(i);    
  }

  uGp = uGrad;
  uGp.at(4) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGp, mV, mVG, tStep);

  for(int i = 1; i <= 5; i++) {
    stiff.at(i,4) = sigma.at(i);    
  }


  uGp = uGrad;
  uGp.at(5) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGp, mV, mVG, tStep);

  for(int i = 1; i <= 5; i++) {
    stiff.at(i,5) = sigma.at(i);    
  }


  stiff.times(1./perturbation);

  */






  MaterialMode matMode = gp->giveMaterialMode();
  
  FloatMatrix Iskew, Isym, I, De;
  
  
  Iskew.beSkewProjectionMatrix(); 
  Isym.beSymProjectionMatrix();
  if (matMode == _PlaneStrain) {

    IsotropicLinearElasticMaterial :: givePlaneStrainStiffMtrx(De, mode, gp, tStep); 

    FloatMatrix IskewFull, IsymFull;
 

    IsymFull = Isym;
    StructuralMaterial :: giveReducedSymMatrixForm(Isym, IsymFull, matMode);
    answer.beProductOf(Isym,De);
    answer = De;
 
    answer.resizeWithData(5,5);
    answer.at(4,5) = answer.at(5,4) = answer.at(5,5) = answer.at(4,4);

    
        
    IskewFull = Iskew;
   
    StructuralMaterial :: giveReducedMatrixForm(Iskew, IskewFull, matMode);
  
  } else { //@todo check that all other modes are threated as 3d modes
    IsotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix(De, mode, gp, tStep);
    answer.beProductOf(De,Isym);
    answer.resizeWithData(9,9);
    answer.at(7,7) = answer.at(4,4);
    answer.at(8,8) = answer.at(5,5);
    answer.at(9,9) = answer.at(6,6);
  }
  
  // I.beTProductOf(Iskew,Iskew);
  Iskew.times(Hk);
  answer.add(Iskew);
}

void
MicropolarMaterial_Elastic :: giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  /*
  MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );  
  double perturbation = 1.e-8;
  
  FloatArray uGrad, mV, mVG;
  uGrad.resize(5);
  mV.resize(1);
  mVG.resize(1);  
  uGrad.zero();
  mV.zero();
  mVG.zero();

  mV = status->giveMicromorphicVar();
  mVG = status->giveMicromorphicVarGrad();
  uGrad = status->giveStrainVector();
  uGrad.resize(5);
  FloatArray mVp(1);
  mVp = mV;
  FloatArray sigma, s, S;
  FloatMatrix stiff(5,1);
  stiff.zero();

  mVp.at(1) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGrad, mVp, mVG, tStep);

  for(int i = 1; i <= 5; i++) {
    stiff.at(i,1) = sigma.at(i);    
  }

  stiff.times(1./perturbation);
  

  */







  MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _PlaneStrain) {
    answer = {{0, 0, 0, 1, -1}};
    //answer = {{0, 0, 0, 0.5, -0.5}};
    answer.times(Hk);
  } else { //@todo check that all other modes are treated as 3d modes
    answer = {{0, 0, 0, -1, 0,  0, 1,  0,  0},
              {0, 0, 0,  0, 1,  0, 0, -1,  0},
              {0, 0, 0,  0, 0, -1, 0,  0,  1}};
    answer.times(-Hk);
  }
  
}


void
MicropolarMaterial_Elastic :: giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  /*
    
  double perturbation = 1.e-8;
  MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );  
  FloatArray uGrad, mV, mVG;
  uGrad.resize(5);
  mV.resize(1);
  mVG.resize(1);  
  uGrad.zero();
  mV.zero();
  mVG.zero();
  mV = status->giveMicromorphicVar();
  mVG = status->giveMicromorphicVarGrad();
  uGrad = status->giveStrainVector();
  uGrad.resize(5);

  FloatArray uGp(5);
  uGp = uGrad;
  FloatArray sigma, s, S;
  FloatMatrix stiff(1,5);
  stiff.zero();

  uGp.at(1) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGp, mV, mVG, tStep);
  stiff.at(1,1) = s.at(1);    
  

  uGp = uGrad;
  uGp.at(2) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGp, mV, mVG, tStep);
  stiff.at(1,2) = s.at(1);    
 


  uGp = uGrad;
  uGp.at(3) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGp, mV, mVG, tStep);
  stiff.at(1,3) = s.at(1);    
  

  uGp = uGrad;
  uGp.at(4) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGp, mV, mVG, tStep);

  stiff.at(1,4) = s.at(1);    
  


  uGp = uGrad;
  uGp.at(5) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGp, mV, mVG, tStep);
  stiff.at(1,5) = s.at(1);    
  

  stiff.times(1./perturbation);
  */
  





  MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _PlaneStrain) {
    //answer = {{0}, {0}, {0}, {0.5}, {-0.5}};
    answer = {{0}, {0}, {0}, {1.}, {-1.0}};
    answer.times(Hk);
  }  else {//@todo check that all other modes are treated as 3d modes
    FloatMatrix Iskew, redIskew, Iw;
    Iskew.beSkewProjectionMatrix(); 
    /*Iw = {{0, 0, 0, 1, 0, 0, -1, 0, 0},
      {0, 0, 0, 0, -1, 0, 0, 1, 0},
      {0, 0, 0, 0, 0, 1, 0, 0, -1}};
      answer.beTProductOf(Iw,Iskew);*/
  }


  //answer.zero()
}



void
MicropolarMaterial_Elastic :: giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _PlaneStrain) {
    answer.resize(1,1);
    answer.at(1,1) = 2.*Hk;
  } else {
    answer.resize(3,3);
    answer.at(1,1) = Hk;
    answer.at(2,2) = Hk;
    answer.at(3,3) = Hk;
  } 

}


void
MicropolarMaterial_Elastic :: giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  if (matMode == _PlaneStrain) {
    answer.resize(2,2);
  } else {
    answer.resize(9,9);
  } 
  answer.beUnitMatrix();
  answer.times(Ak);

}




void
MicropolarMaterial_Elastic :: giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &S, GaussPoint *gp, const FloatArray &displacementGradient, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{
  
    MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );
    MaterialMode matMode = gp->giveMaterialMode();
    //symmetric and skew part of the strain
    FloatArray strain, vRotation;
    // stiffnes matrix and stress in matrix form
    FloatMatrix De;

    strain = displacementGradient;    
    if(matMode == _PlaneStrain) {
      IsotropicLinearElasticMaterial :: givePlaneStrainStiffMtrx(De, TangentStiffness, gp, tStep);   
      
      strain.resize(4);   
      strain.at(4) =  (displacementGradient.at(4) + displacementGradient.at(5));
      
      vRotation.resize(1);
      vRotation.at(1) = 0.5*(displacementGradient.at(4) - displacementGradient.at(5));	

    } else {
      IsotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix(De, TangentStiffness, gp, tStep);
      strain.resize(6);
      strain.at(4) = (displacementGradient.at(4) + displacementGradient.at(7));
      strain.at(5) = (displacementGradient.at(5) + displacementGradient.at(8));
      strain.at(6) = (displacementGradient.at(6) + displacementGradient.at(9));
      vRotation.resize(3);
      vRotation.at(1) = (displacementGradient.at(4) - displacementGradient.at(7))/2;
      vRotation.at(2) = (displacementGradient.at(5) - displacementGradient.at(8))/2;
      vRotation.at(3) = (displacementGradient.at(6) - displacementGradient.at(9))/2;
    } 
    FloatArray reducedSigma, fullSigma;
    FloatMatrix mSigma;
    reducedSigma.beProductOf(De,strain);     
    StructuralMaterial :: giveFullSymVectorForm(fullSigma, reducedSigma, matMode);
    mSigma.beMatrixForm(fullSigma);
    fullSigma.beVectorForm(mSigma);
    
   

    FloatArray vCoupledStress;
    FloatMatrix mRotation, microRotation, coupledStress, higherStress; 
    
   
    S = micromorphicVarGrad;
    S.times(Ak);



    microRotation.giveMatrixOfAxialVector(micromorphicVar);
    mRotation.beSkewMatrixForm(vRotation);
    coupledStress = microRotation;
    coupledStress.subtract(mRotation);
    coupledStress.times(Hk);
    vCoupledStress.beVectorForm(coupledStress);

    if(matMode == _PlaneStrain) {
      s.resize(1);
      s.at(1) = (vCoupledStress.at(9) - vCoupledStress.at(6));
      //s.at(1) = 0.5*(vCoupledStress.at(9) - vCoupledStress.at(6));
    }





    FloatMatrix Iskew;
    FloatArray IskewS;
    //Iskew.beSkewProjectionMatrix(); 
    //IskewS.beProductOf(Iskew, vCoupledStress);
    fullSigma.subtract(vCoupledStress);
    StructuralMaterial :: giveReducedVectorForm(sigma, fullSigma, matMode);
    


    status->letTempStrainVectorBe(displacementGradient);
    status->letTempMicromorphicVarBe(micromorphicVar);
    status->letTempMicromorphicVarGradBe(micromorphicVarGrad);

    status->letTempStressVectorBe(sigma);
    status->letTempMicromorphicStressBe(s);
    status->letTempMicromorphicStressGradBe(S); 
      
}



IRResultType
MicropolarMaterial_Elastic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    IsotropicLinearElasticMaterial :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, Hk, _IFT_MicromorphicMaterialExtensionInterface_Hk);
    IR_GIVE_FIELD(ir, Ak, _IFT_MicromorphicMaterialExtensionInterface_Ak);

    return IRRT_OK;
}

int
MicropolarMaterial_Elastic :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );
    

    if( type == IST_MicromorphicStress) {
      StructuralMaterial :: giveFullVectorForm( answer, status->giveStressVector(), gp->giveMaterialMode() );

    } else if( type == IST_MicromorphicStrain ) {
      StructuralMaterial :: giveFullVectorForm( answer, status->giveStrainVector(), gp->giveMaterialMode() );
    } else if( type == IST_MicromorphicRelativeStress ) {
      FloatArray s = status->giveMicromorphicStress();
      // @todo this is correct only for 2d problems, something like giveFullVectorForm for MicromorphicMaterial class would be necessary
      answer.resize(9);
      answer.zero();
      answer.at(6) = s.at(1);
      answer.at(9) = -s.at(1);
    } else if ( type == IST_MicromorphicRelativeStrain ) {
      answer.resize(9);
      answer.zero();
    } else if ( type == IST_MicromorphicHigherOrderStress ) {
      FloatArray M = status->giveMicromorphicStressGrad(); 
      answer.resize(9);
      answer.zero();
      answer.at(6) = M.at(1);
      answer.at(9) = M.at(2);
    } else if ( type == IST_MicromorphicHigherOrderStrain ) {
      FloatArray kappa = status->giveMicromorphicVarGrad();
      answer.resize(9);
      answer.zero();
      answer.at(6) = kappa.at(1);
      answer.at(9) = kappa.at(2);
    } else {
      OOFEM_ERROR("Unknown InternalStateType");
    }
    return 1;
}
    



} // end namespace oofem
