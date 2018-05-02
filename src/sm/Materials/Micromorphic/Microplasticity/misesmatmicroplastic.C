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

#include "misesmatmicroplastic.h"
#include "stressvector.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"
#include "Materials/isolinearelasticmaterial.h"


namespace oofem {
  REGISTER_Material(MisesMatMicroplastic);

  MisesMatMicroplastic :: MisesMatMicroplastic(int n, Domain *d) : MisesMat(n, d),MicromorphicMaterialExtensionInterface(d)//,IsotropicLinearElasticMaterial(n, d)
{
  linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);

  Hk = Ak = 0.;
}

MisesMatMicroplastic :: ~MisesMatMicroplastic()
{ }

// creates a new material status  corresponding to this class
MaterialStatus *
MisesMatMicroplastic :: CreateStatus(GaussPoint *gp) const
{
    return new MisesMatMicroplasticStatus(1, this->giveDomain(), gp);
}

void
MisesMatMicroplastic :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}




void
MisesMatMicroplastic :: giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  
  if(gp->giveMaterialMode() == _1dMat) {
    this->giveLinearElasticMaterial()->give1dStressStiffMtrx(answer, mode, gp, tStep);
    MisesMatMicroplasticStatus *status = static_cast< MisesMatMicroplasticStatus * >( this->giveStatus(gp) );
    double kappa = status->giveCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double E = answer.at(1, 1);
    if ( mode != TangentStiffness ) {
        return;
    }
    if ( tempKappa <= kappa || tempKappa < 1.e-15 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
      return;
    }
    // === plastic loading ===
    answer.resize(1, 1);
    answer.at(1, 1) =  E * (H+Hk) / ( E + H + Hk );    

    // start from the elastic stiffness
  }  else {
    this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
  }
  if ( mode != TangentStiffness ) {
    if(gp->giveMaterialMode() == _PlaneStrain) {
      FloatMatrix plS;
      plS.resize(4, 4);
      plS.zero();
      //answer.beSubMatrixOf(m3d, indx, indx);
      
      plS.at(1, 1) = answer.at(1, 1);
      plS.at(1, 2) = answer.at(1, 2);
      plS.at(1, 4) = answer.at(1, 6);
      
      plS.at(2, 1) = answer.at(2, 1);
      plS.at(2, 2) = answer.at(2, 2);
      plS.at(2, 4) = answer.at(2, 6);
      
      plS.at(3, 1) = answer.at(3, 1);
      plS.at(3, 2) = answer.at(3, 2);
      plS.at(3, 4) = answer.at(3, 6);
      
      plS.at(4, 1) = answer.at(6, 1);
      plS.at(4, 2) = answer.at(6, 2);
      plS.at(4, 4) = answer.at(6, 6);
      
      answer = plS;
    }
    return;
  }
   
  MisesMatMicroplasticStatus *status = static_cast< MisesMatMicroplasticStatus * >( this->giveStatus(gp) );
    double kappa = status->giveCumulativePlasticStrain();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = tempKappa - kappa;
    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
      if(gp->giveMaterialMode() == _PlaneStrain) {
	FloatMatrix plS;
	plS.resize(4, 4);
	plS.zero();
	//answer.beSubMatrixOf(m3d, indx, indx);
	
	plS.at(1, 1) = answer.at(1, 1);
	plS.at(1, 2) = answer.at(1, 2);
	plS.at(1, 3) = answer.at(1, 3);
	plS.at(1, 4) = answer.at(1, 6);
	
	plS.at(2, 1) = answer.at(2, 1);
	plS.at(2, 2) = answer.at(2, 2);
	plS.at(2, 3) = answer.at(2, 3);
	plS.at(2, 4) = answer.at(2, 6);
	
	plS.at(3, 1) = answer.at(3, 1);
	plS.at(3, 2) = answer.at(3, 2);
	plS.at(3, 3) = answer.at(3, 3);
	plS.at(3, 4) = answer.at(3, 6);
	
	plS.at(4, 1) = answer.at(6, 1);
	plS.at(4, 2) = answer.at(6, 2);
	plS.at(4, 3) = answer.at(6, 3);
	plS.at(4, 4) = answer.at(6, 6);
	
	answer = plS;
      }
      
      return;
    }
    // === plastic loading ===
    // yield stress at the beginning of the step
    FloatArray micromorphicVar;
    micromorphicVar = status->giveTempMicromorphicVar();
    
    double sigmaY = sig0 + (H+Hk) * kappa - Hk*micromorphicVar.at(1);
    // trial deviatoric stress and its norm
    const FloatArray &trialStressDev = status->giveTrialStressDev();
    double trialS = computeStressNorm(trialStressDev);

    /////////////////////////////////////////////
    double perturbation = 1.e-8;
    FloatArray uGrad, mV, mVG, oldS;
    mV = status->giveTempMicromorphicVar();
    mVG = status->giveTempMicromorphicVarGrad();
    uGrad = status->giveTempStrainVector();
    oldS = status->giveTempMicromorphicStress();
    FloatArray sigma = status->giveTempStressVector();    
    FloatArray uGradp;
    uGradp = uGrad;
    FloatArray sigmaP, s, Sp;
    FloatMatrix stiff(4,4);
    stiff.zero();
    uGradp.at(1) += (perturbation);
    this-> giveGeneralizedStressVectors(sigmaP, s, Sp, gp, uGradp, mV, mVG, tStep);
    for(int i =1; i <=4; i++) {
      stiff.at(i,1) = sigmaP.at(i) - sigma.at(1);    
    }
    uGradp = uGrad;
    uGradp.at(2) += (perturbation);
    this-> giveGeneralizedStressVectors(sigmaP, s, Sp, gp, uGradp, mV, mVG, tStep);
    for(int i =1; i <=4; i++) {
      stiff.at(i,2) = sigmaP.at(i) - sigma.at(1);    
    }
    
    uGradp = uGrad;
    uGradp.at(3) += (perturbation);
    this-> giveGeneralizedStressVectors(sigmaP, s, Sp, gp, uGradp, mV, mVG, tStep);
    for(int i =1; i <=4; i++) {
      stiff.at(i,3) = sigmaP.at(i) - sigma.at(1);    
    }

    uGradp = uGrad;
    uGradp.at(4) += (perturbation);
    this-> giveGeneralizedStressVectors(sigmaP, s, Sp, gp, uGradp, mV, mVG, tStep);
    for(int i =1; i <=4; i++) {
      stiff.at(i,4) = sigmaP.at(i) - sigma.at(1);    
    }
   
    stiff.times(1./perturbation);
    /////////////////////////////////////////////



    // one correction term
    FloatMatrix stiffnessCorrection;
    stiffnessCorrection.beDyadicProductOf(trialStressDev, trialStressDev);
    double factor = -2. * sqrt(6.) * G * G / trialS;
    double factor1 = factor * sigmaY / ( ( H + Hk + 3. * G ) * trialS * trialS );
    answer.add(factor1, stiffnessCorrection);
    // another correction term
    stiffnessCorrection.bePinvID();
    double factor2 = factor * dKappa;
    answer.add(factor2, stiffnessCorrection);

    
    if(gp->giveMaterialMode() == _PlaneStrain) {
	FloatMatrix plS;
	plS.resize(4, 4);
	plS.zero();
	//answer.beSubMatrixOf(m3d, indx, indx);
	
	plS.at(1, 1) = answer.at(1, 1);
	plS.at(1, 2) = answer.at(1, 2);
	plS.at(1, 3) = answer.at(1, 3);
	plS.at(1, 4) = answer.at(1, 6);
	
	plS.at(2, 1) = answer.at(2, 1);
	plS.at(2, 2) = answer.at(2, 2);
	plS.at(2, 3) = answer.at(2, 3);
	plS.at(2, 4) = answer.at(2, 6);
	
	plS.at(3, 1) = answer.at(3, 1);
	plS.at(3, 2) = answer.at(3, 2);
	plS.at(3, 3) = answer.at(3, 3);
	plS.at(3, 4) = answer.at(3, 6);
	
	plS.at(4, 1) = answer.at(6, 1);
	plS.at(4, 2) = answer.at(6, 2);
	plS.at(4, 3) = answer.at(6, 3);
	plS.at(4, 4) = answer.at(6, 6);
	
	answer = plS;
      }

}

void
MisesMatMicroplastic :: giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  MisesMatMicroplasticStatus *status = static_cast< MisesMatMicroplasticStatus * >( this->giveStatus(gp) );
  MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _1dMat) {
    answer.resize(1,1);
  } else if(matMode == _PlaneStrain) {
    answer.resize(4,1);
  } else {
    answer.resize(6,1);
  }
  answer.zero();

  double kappa = status->giveCumulativePlasticStrain();
  // increment of cumulative plastic strain as an indicator of plastic loading
  double tempKappa = status->giveTempCumulativePlasticStrain();
  if ( tempKappa <= kappa || tempKappa < 1.e-15 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
    return;
  }

  
  LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
  double E = lmat->give('E', gp);
  FloatArray sigma = status->giveStressVector();
    // === plastic loading ===
  if(matMode == _1dMat) {
    answer.resize(1,1);
    answer.at(1,1) = -Hk*E/(E+H+Hk)*signum(sigma.at(1));
  } else {
    
    double perturbation = 1.e-10;
    FloatArray uGrad, mV, mVG, oldSigma;    
    mV = status->giveTempMicromorphicVar();
    mVG = status->giveTempMicromorphicVarGrad();
    uGrad = status->giveTempStrainVector();
    oldSigma = status->giveTempStressVector(); 
    FloatArray mVp(1);
    mVp = mV;
    FloatArray sigma, s, S;
    FloatMatrix stiff(4,1);
    stiff.zero();    
    mVp.at(1) += (perturbation);
    this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGrad, mVp, mVG, tStep);
    stiff = sigma;
    stiff.subtract(oldSigma);        
    stiff.times(1./perturbation);
  

    FloatArray reducedStress;
    FloatArray fullStress = status->giveTrialStressDev();    
    double stressNorm = computeStressNorm(fullStress);
    StructuralMaterial :: giveReducedSymVectorForm(reducedStress, fullStress, _PlaneStrain);
    answer = reducedStress;
    answer.times(1./stressNorm);
    answer.times(-sqrt(6)*G*Hk/(3.*G+Hk+H));
  }
  
}


void
MisesMatMicroplastic :: giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  MisesMatMicroplasticStatus *status = static_cast< MisesMatMicroplasticStatus * >( this->giveStatus(gp) );
  MaterialMode matMode = gp->giveMaterialMode();  
  if(matMode == _1dMat) {
    answer.resize(1,1);
  } else if(matMode == _PlaneStrain) {
    answer.resize(1,4);
  } else {
    answer.resize(6,1);
  }
  answer.zero();


  double kappa = status->giveCumulativePlasticStrain();
  // increment of cumulative plastic strain as an indicator of plastic loading
  double tempKappa = status->giveTempCumulativePlasticStrain();
  if ( tempKappa <= kappa || tempKappa < 1.e-15 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
    return;
  }
  // === plastic loading ===
  LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
  double E = lmat->give('E', gp);


  if(matMode == _1dMat) {
    answer.at(1,1) = -Hk*E/(E+H+Hk);
  } else {

    FloatArray reducedStress;
    FloatArray fullStress = status->giveTrialStressDev();    
    double stressNorm = computeStressNorm(fullStress);
    StructuralMaterial :: giveReducedSymVectorForm(reducedStress, fullStress, _PlaneStrain);
    for(int i = 1; i<= reducedStress.giveSize();i++) {
      if(i <= 3) {
	answer.at(1,i) = reducedStress.at(i);
      } else {
	answer.at(1,i) = reducedStress.at(i);
      }
    }
    answer.times(1./stressNorm);
    answer.times(-sqrt(6)*G*Hk/(3.*G+Hk+H));
  }


  
    double perturbation = -1.e-8;
    FloatArray uGrad, mV, mVG, oldS;
    mV = status->giveTempMicromorphicVar();
    mVG = status->giveTempMicromorphicVarGrad();
    uGrad = status->giveTempStrainVector();
    oldS = status->giveTempMicromorphicStress();
    
    FloatArray uGradp(1);
    uGradp = uGrad;
    FloatArray sigma, s, S;
    FloatMatrix stiff(1,4);
    stiff.zero();
    //  this-> giveGeneralizedStressVectors(sigma, oldS, S, gp, uGrad, mV, mVG, tStep);
    uGradp.at(1) += (perturbation);
    this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGradp, mV, mVG, tStep);
    stiff.at(1,1) = s.at(1) - oldS.at(1);    
    
    uGradp = uGrad;
    uGradp.at(2) += (perturbation);
    this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGradp, mV, mVG, tStep);
    stiff.at(1,2) = s.at(1) - oldS.at(1);    
    
    uGradp = uGrad;
    uGradp.at(3) += (perturbation);
    this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGradp, mV, mVG, tStep);
    stiff.at(1,3) = s.at(1) - oldS.at(1);    

    uGradp = uGrad;
    uGradp.at(4) += (perturbation);
    this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGradp, mV, mVG, tStep);
    stiff.at(1,4) = s.at(1) - oldS.at(1);    
    stiff.times(1./perturbation);
    




   
  
}



void
MisesMatMicroplastic :: giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  
  MisesMatMicroplasticStatus *status = static_cast< MisesMatMicroplasticStatus * >( this->giveStatus(gp) );
  MaterialMode matMode = gp->giveMaterialMode();
  
  double kappa = status->giveCumulativePlasticStrain();
  // increment of cumulative plastic strain as an indicator of plastic loading
  double tempKappa = status->giveTempCumulativePlasticStrain();
 
  answer.resize(1,1);
  answer.at(1,1) = Hk;
  
  if ( tempKappa <= kappa || tempKappa < 1.e-15 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
    return;
  }

  LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
  double E = lmat->give('E', gp);
  
  double perturbation = 1.e-8;
  FloatArray uGrad, mV, mVG, oldS;
  
  mV = status->giveTempMicromorphicVar();
  mVG = status->giveTempMicromorphicVarGrad();
  uGrad = status->giveTempStrainVector();
  oldS = status->giveTempMicromorphicStress();
  
  FloatArray mVp(1);
  mVp = mV;
  FloatArray sigma, s, S;
  FloatMatrix stiff(1,1);
  stiff.zero();
  mVp.at(1) += (perturbation);
  this-> giveGeneralizedStressVectors(sigma, s, S, gp, uGrad, mVp, mVG, tStep);
  stiff.at(1,1) = s.at(1) - oldS.at(1);    
  
  stiff.times(1./perturbation);
  
  if(matMode == _1dMat){
    answer.at(1,1) -= Hk*Hk/(E+H+Hk);
  } else {
    answer.at(1,1) -= Hk*Hk/(3.*G + H + Hk);
  }
}


void
MisesMatMicroplastic :: giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _1dMat) {
    answer.resize(1,1);
  } else if (matMode == _PlaneStrain) {
    answer.resize(2,2);
  } else {
    answer.resize(9,9);
  } 
  answer.beUnitMatrix();
  answer.times(Ak);

}




void
MisesMatMicroplastic :: giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &S, FloatArray &M, GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{
  
    MisesMatMicroplasticStatus *status = static_cast< MisesMatMicroplasticStatus * >( this->giveStatus(gp) );
    //this->initGpForNewStep(gp);
    this->initTempStatus(gp);
    this->performPlasticityReturn(gp, totalStrain, micromorphicVar, micromorphicVarGrad);
    
    sigma = status->giveTempStressVector();
    S = status->giveTempMicromorphicStress();
    M = status->giveTempMicromorphicStressGrad();
    

    status->letTempStrainVectorBe(totalStrain);
    status->letTempMicromorphicVarBe(micromorphicVar);
    status->letTempMicromorphicVarGradBe(micromorphicVarGrad);
  
     
      
}





void
MisesMatMicroplastic :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &micromorphicVar, const FloatArray &micromorphicVarGrad)
{
    MisesMatMicroplasticStatus *status = static_cast< MisesMatMicroplasticStatus * >( this->giveStatus(gp) );
    double kappa;
    FloatArray plStrain;
    FloatArray fullStress, answer;
    // get the initial plastic strain and initial kappa from the status
    plStrain = status->givePlasticStrain();
    kappa = status->giveCumulativePlasticStrain();

    /// micromorphic stresses
    double Str = Hk * (micromorphicVar.at(1) - kappa);
    FloatArray S(1), M(1);

    // === radial return algorithm ===
    if ( totalStrain.giveSize() == 1 ) {
        LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
        double E = lmat->give('E', gp);
        /*trial stress*/
        fullStress.resize(1);
	answer.resize(1);
        fullStress.at(1) = E * ( totalStrain.at(1) - plStrain.at(1) );
        double trialS = fabs(fullStress.at(1));
        /*yield function*/
        double yieldValue = trialS - (sig0 + H * kappa) + Str;
        // === radial return algorithm ===
        if ( yieldValue > 0 ) {
            double dKappa = yieldValue / ( H + Hk +  E );
            kappa += dKappa;
            plStrain.at(1) += dKappa * signum( fullStress.at(1) );
            plStrain.at(2) -= 0.5 * dKappa * signum( fullStress.at(1) );
            plStrain.at(3) -= 0.5 * dKappa * signum( fullStress.at(1) );
            answer.at(1) -= dKappa * E * signum( fullStress.at(1) );
	    //check
	    /*
	    Str = Hk * (micromorphicVar.at(1) - kappa);
	    yieldValue = answer.at(1) - (sig0 + H * kappa) + Str;  
	    int ahoj = 1;
	    */
	}
    } else {
        // elastic predictor
      FloatArray elStrain = totalStrain;
      if(gp->giveMaterialMode() == _PlaneStrain) {
	StructuralMaterial :: giveFullSymVectorForm(elStrain, totalStrain, _PlaneStrain);                                              
      }
        elStrain.subtract(plStrain);
        FloatArray elStrainDev;
        double elStrainVol;
        elStrainVol = computeDeviatoricVolumetricSplit(elStrainDev, elStrain);
        FloatArray trialStressDev;
        applyDeviatoricElasticStiffness(trialStressDev, elStrainDev, G);
        /**************************************************************/
        double trialStressVol = 3 * K * elStrainVol;
        /**************************************************************/

        // store the deviatoric and trial stress (reused by algorithmic stiffness)
        status->letTrialStressDevBe(trialStressDev);
        // check the yield condition at the trial state
        double trialS = computeStressNorm(trialStressDev);
        double yieldValue = sqrt(3./2.) * trialS - (sig0 + H * kappa) + Str;
        if ( yieldValue > 0. ) {
            // increment of cumulative plastic strain
            double dKappa = yieldValue / ( H + Hk + 3. * G );
            kappa += dKappa;
            FloatArray dPlStrain;
            // the following line is equivalent to multiplication by scaling matrix P
            applyDeviatoricElasticCompliance(dPlStrain, trialStressDev, 0.5);
            // increment of plastic strain
            plStrain.add(sqrt(3. / 2.) * dKappa / trialS, dPlStrain);
            // scaling of deviatoric trial stress
            trialStressDev.times(1. - sqrt(6.) * G * dKappa / trialS);
        }

        // assemble the stress from the elastically computed volumetric part
        // and scaled deviatoric part

        computeDeviatoricVolumetricSum(fullStress, trialStressDev, trialStressVol);
	
	if(gp->giveMaterialMode() == _PlaneStrain) {
	  StructuralMaterial :: giveReducedSymVectorForm(answer, fullStress, _PlaneStrain);                                         
	} else {

	  answer = fullStress;
	}
    }
    //S.resize(1);
    //M.resize(1);
    //S.at(1) = Hk * (micromorphicVar.at(1) - kappa);
    //M.at(1) = Ak * micromorphicVarGrad.at(1);
    S.at(1) = Hk * (micromorphicVar.at(1) - kappa);
    M = Ak * micromorphicVarGrad;

    // store the plastic strain and cumulative plastic strain
    status->letTempPlasticStrainBe(plStrain);
    status->setTempCumulativePlasticStrain(kappa);
    // store stresses into the status
    status->letTempStressVectorBe(answer);
    status->letTempMicromorphicStressBe(S);
    status->letTempMicromorphicStressGradBe(M);
    
}



IRResultType
MisesMatMicroplastic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    MisesMat :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, Hk, _IFT_MicromorphicMaterialExtensionInterface_Hk);
    IR_GIVE_FIELD(ir, Ak, _IFT_MicromorphicMaterialExtensionInterface_Ak);

    return IRRT_OK;
}

  /*
int
MisesMatMicroplastic :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
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
  */    

//=============================================================================

MisesMatMicroplasticStatus :: MisesMatMicroplasticStatus(int n, Domain *d, GaussPoint *g) :
    MicromorphicMaterialStatus(n, d, g), plasticStrain(6), tempPlasticStrain(), trialStressDev()
{
    stressVector.resize(6);
    strainVector.resize(6);
    kappa = tempKappa = 0.;
   
}

MisesMatMicroplasticStatus :: ~MisesMatMicroplasticStatus()
{ }

  /*void
MisesMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    
    }*/


// initializes temporary variables based on their values at the previous equlibrium state
void MisesMatMicroplasticStatus :: initTempStatus()
{
    MicromorphicMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
      plasticStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
      plasticStrain.zero();
    }
    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
    trialStressDev.clear(); // to indicate that it is not defined yet
}


// updates internal variables when equilibrium is reached
void
MisesMatMicroplasticStatus :: updateYourself(TimeStep *tStep)
{
  MicromorphicMaterialStatus :: updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    kappa = tempKappa;
    trialStressDev.clear(); // to indicate that it is not defined any more
}



} // end namespace oofem
