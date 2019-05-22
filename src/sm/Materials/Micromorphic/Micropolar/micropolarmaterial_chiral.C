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

#include "../sm/Materials/Micromorphic/Micropolar/micropolarmaterial_chiral.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
  REGISTER_Material(MicropolarMaterial_Chiral);

MicropolarMaterial_Chiral :: MicropolarMaterial_Chiral(int n, Domain *d) :IsotropicLinearElasticMaterial(n, d), MicromorphicMaterialExtensionInterface(d)
{
  
}

MicropolarMaterial_Chiral :: ~MicropolarMaterial_Chiral()
{ }



void
MicropolarMaterial_Chiral :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}



void
MicropolarMaterial_Chiral :: giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
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
 
    answer.beIsotropicMatrix(lambda, mu+mu_c, mu);

}

void
MicropolarMaterial_Chiral :: giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
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
    answer = {{0, 0, 0, -1, 0,  0, 1,  0,  0},
              {0, 0, 0,  0, 1,  0, 0, -1,  0},
              {0, 0, 0,  0, 0, -1, 0,  0,  1}};
    answer.times(Hk);

  
}



void
MicropolarMaterial_Chiral :: giveMicromorphicMatrix_dSigdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
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
  answer.beIsotropicMatrix(C1, C2, C3);
    
  
}


  

void
MicropolarMaterial_Chiral :: giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
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
  


    answer = {{0, 0, 0, -1, 0,  0, 1,  0,  0},
              {0, 0, 0,  0, 1,  0, 0, -1,  0},
              {0, 0, 0,  0, 0, -1, 0,  0,  1}};
    answer.times(Hk);

}



void
MicropolarMaterial_Chiral :: giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    answer.resize(3,3);
    answer.beUnitMatrix();
    answer.times(Hk);
}


void
MicropolarMaterial_Chiral :: giveMicromorphicMatrix_dSdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    answer = {{0, 0, 0, -1, 0,  0, 1,  0,  0},
              {0, 0, 0,  0, 1,  0, 0, -1,  0},
	    {0, 0, 0,  0, 0, -1, 0,  0,  1}};
    answer.times(C2-C3);

    
}



void
MicropolarMaterial_Chiral :: giveMicromorphicMatrix_dMdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    answer.beIsotropicMatrix(C1, C2, C3);
}


void
MicropolarMaterial_Chiral :: giveMicromorphicMatrix_dMdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    answer = {{0, 0, 0, -1, 0,  0, 1,  0,  0},
              {0, 0, 0,  0, 1,  0, 0, -1,  0},
	    {0, 0, 0,  0, 0, -1, 0,  0,  1}};
    answer.times(C2-C3);

    
}


  
  
void
MicropolarMaterial_Chiral :: giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  answer.beIsotropicMatrix(alpha, beta, gamma);

}




void
MicropolarMaterial_Chiral :: giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &m, GaussPoint *gp, const FloatArray &displacementGradient, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{
  
    MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );
    MaterialMode matMode = gp->giveMaterialMode();
    //symmetric and skew part of the strain
    FloatArray strain, vRotation, vMicroRotation;
    FloatArray micropolarStrain(displacementGradient), aux;
    FloatMatrix microRotation, LeviCivita;


    LeviCivita = {{0, 0, 0, -1, 0,  0, 1,  0,  0},
              {0, 0, 0,  0, 1,  0, 0, -1,  0},
	    {0, 0, 0,  0, 0, -1, 0,  0,  1}};

    
    microRotation.giveMatrixOfAxialVector(micromorphicVar);
    vMicroRotation.beVectorForm(microRotation);
    micropolarStrain.add(vMicroRotation);
    FloatMatrix E, C, A;
    E.beIsotropicMatrix(lambda, mu+mu_c, mu);
    A.beIsotropicMatrix(alpha, beta, gamma);
    C.beIsotropicMatrix(C1, C2, C3);
   
    sigma.beProductOf(E,micropolarStrain);
    aux.beProductOf(C,micromorphicVarGrad);
    sigma.add(aux);

    

    

    s.beProductOf(LeviCivita, sigma);
    
    m.beProductOf(A, micromorphicVarGrad);
    aux.beTProductOf(C, micropolarStrain);
    m.add(aux);
      
    status->letTempStrainVectorBe(displacementGradient);
    status->letTempMicromorphicVarBe(micromorphicVar);
    status->letTempMicromorphicVarGradBe(micromorphicVarGrad);

    status->letTempStressVectorBe(sigma);
    status->letTempMicromorphicStressBe(s);
    status->letTempMicromorphicStressGradBe(m); 
      
}



IRResultType
MicropolarMaterial_Chiral :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    //IsotropicLinearElasticMaterial :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, lambda, _IFT_MicromorphicMaterialExtensionInterface_Hk);
    IR_GIVE_FIELD(ir, mu, _IFT_MicromorphicMaterialExtensionInterface_Ak);
    IR_GIVE_FIELD(ir, mu_c, _IFT_MicromorphicMaterialExtensionInterface_Ak);


    IR_GIVE_FIELD(ir, alpha, _IFT_MicropolarMaterial_Chiral_alpha);
    IR_GIVE_FIELD(ir, beta, _IFT_MicropolarMaterial_Chiral_beta);
    IR_GIVE_FIELD(ir, gamma, _IFT_MicropolarMaterial_Chiral_gamma);

    IR_GIVE_FIELD(ir, C1, _IFT_MicropolarMaterial_Chiral_C1);
    IR_GIVE_FIELD(ir, C2, _IFT_MicropolarMaterial_Chiral_C2);
    IR_GIVE_FIELD(ir, C3, _IFT_MicropolarMaterial_Chiral_C3);

    return IRRT_OK;
}

int
MicropolarMaterial_Chiral :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
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
