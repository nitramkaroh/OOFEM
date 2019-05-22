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

#include "../sm/Materials/Micromorphic/GradientPolyconvex/gradientpolyconvexmaterial_defgrad.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
  REGISTER_Material(GradientF_PolyconvexMaterial);

GradientF_PolyconvexMaterial :: GradientF_PolyconvexMaterial(int n, Domain *d) :IsotropicLinearElasticMaterial(n, d), MicromorphicMaterialExtensionInterface(d)
{
  Hk = Ak = 0.;
  alpha = 1.e6;
}

GradientF_PolyconvexMaterial :: ~GradientF_PolyconvexMaterial()
{ }



void
GradientF_PolyconvexMaterial :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}



void
GradientF_PolyconvexMaterial :: giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    GradientF_PolyconvexMaterialStatus *status = static_cast< GradientF_PolyconvexMaterialStatus * >( this->giveStatus(gp) );

    answer.resize(9,9);
    answer.zero();
    //deformation gradient, its inverse, cofactor, and determinant
    double normC_tC1, normC_tC2;
    FloatArray vF, vB, delta, vC_tC1, vC_tC2;
    FloatMatrix F, invF, cofF, B, C, C_tC1, C_tC2, dCdF;
    delta = {1, 1, 1, 0, 0, 0, 0, 0, 0};
    F.beMatrixForm(status->giveTempFVector());
    B.beProductTOf(F,F);
    vB.beVectorForm(B);
    invF.beInverseOf(F);
    double J = F.giveDeterminant();
    cofF.beTranspositionOf(invF);
    cofF.times(J);
    vF.beVectorForm(F);


    FloatArray vInvF;
    vInvF.beVectorForm(invF);
    

    ///////////
    C.beTProductOf(F,F);

    C_tC1 = C;
    C_tC1.subtract(tC1);

    C_tC2 = C;
    C_tC2.subtract(tC2);
    

    normC_tC1 = C_tC1.computeFrobeniusNorm();
    normC_tC2 = C_tC2.computeFrobeniusNorm();

    vC_tC1.beVectorForm(C_tC1);
    vC_tC2.beVectorForm(C_tC2);
    this->compute_dC_dF(dCdF,vF);
    
   

    FloatArray vC_tC1_dCdF, vC_tC2_dCdF;
    vC_tC1_dCdF.beTProductOf(dCdF, vC_tC1);
    vC_tC2_dCdF.beTProductOf(dCdF, vC_tC2);

    
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
	for (int p = 1; p <= 3; p++) {
	  for (int q = 1; q <= 3; q++) {
	    // double well material stiffness
	    answer.at(giveVI(i,j),giveVI(p,q)) += 4. * ( vB.at(giveVI(i,p)) * delta.at(giveVI(j,q)) + vF.at(giveVI(i,q)) *  vF.at(giveVI(p,j))) * (normC_tC1*normC_tC1 + normC_tC2 * normC_tC2) + 4. * vC_tC1.at(giveVI(j,q)) * delta.at(giveVI(i,p)) * (normC_tC2 * normC_tC2) + 4. * vC_tC2.at(giveVI(j,q)) * delta.at(giveVI(i,p)) * (normC_tC1 * normC_tC1) + 4.*vC_tC1_dCdF.at(giveVI(i,j))*vC_tC2_dCdF.at(giveVI(p,q)) + 4.*vC_tC2_dCdF.at(giveVI(i,j))*vC_tC1_dCdF.at(giveVI(p,q));
    	    /// part with J^-1
	    //answer.at(giveVI(i,j),giveVI(p,q)) += 1./J * (vInvF.at(giveVI(j,i)) *  vInvF.at(giveVI(q,p)) + vInvF.at(giveVI(j,p)) *  vInvF.at(giveVI(q,i)));
	  }
	}
      }
    }

    answer.times(alpha);


    FloatMatrix unitMatrix(9,9);
    unitMatrix.beUnitMatrix();
    unitMatrix.times(Hk);
    answer.add(unitMatrix);
    

    

    if(gp->giveMaterialMode() == _PlaneStrain) {
      FloatMatrix m3d = answer;
      answer.resize(5,5);
      answer.zero();
      answer.at(1, 1) = m3d.at(1, 1);
      answer.at(1, 2) = m3d.at(1, 2);
      //answer.at(1, 3) = m3d.at(1, 3);
      answer.at(1, 4) = m3d.at(1, 6);
      answer.at(1, 5) = m3d.at(1, 9);
      
      answer.at(2, 1) = m3d.at(2, 1);
      answer.at(2, 2) = m3d.at(2, 2);
      //answer.at(2, 3) = m3d.at(2, 3);
      answer.at(2, 4) = m3d.at(2, 6);
      answer.at(2, 5) = m3d.at(2, 9);
      
      answer.at(3, 1) = m3d.at(3, 1);
      answer.at(3, 2) = m3d.at(3, 2);
      //answer.at(3, 3) = m3d.at(3, 3);
      answer.at(3, 4) = m3d.at(3, 6);
      answer.at(3, 5) = m3d.at(3, 9);
      
      answer.at(4, 1) = m3d.at(6, 1);
      answer.at(4, 2) = m3d.at(6, 2);
      //answer.at(4, 3) = m3d.at(6, 3);
      answer.at(4, 4) = m3d.at(6, 6);
      answer.at(4, 5) = m3d.at(6, 9);
      
      answer.at(5, 1) = m3d.at(9, 1);
      answer.at(5, 2) = m3d.at(9, 2);
      //answer.at(5, 3) = m3d.at(9, 3);
      answer.at(5, 4) = m3d.at(9, 6);
      answer.at(5, 5) = m3d.at(9, 9);
    }



    /*
    FloatMatrix testA(9,9);
    // relative stress for cofactor
    FloatArray s, micromorphicVar, mvF;
    mvF = status->giveTempMicromorphicVar();
    s = vF;
    s.subtract(mvF);
    s.times(-Hk);
    FloatArray vP, pvP, M, pertF, ps, reducedvF, reducedMicromorphicVar(9), micromorphicVarGrad(1);
    FloatArray col;
    double e = 1.e-6;
    reducedvF.beVectorForm(F);
    reducedMicromorphicVar = status->giveTempMicromorphicVar();
    this->giveFiniteStrainGeneralizedStressVectors_3d(vP, s, M, gp, reducedvF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
    pertF = reducedvF;
    pertF.at(1) += e;
    this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, s, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
    col = pvP;
    col.subtract(vP);
    col.times(1./e);
    testA.setColumn(col, 1);
    pertF = reducedvF;
    pertF.at(2) += e;
    this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, s, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
    col = pvP;
    col.subtract(vP);
    col.times(1./e);
    testA.setColumn(col, 2);
    pertF = reducedvF;
    pertF.at(3) += e;
    this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, s, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
    col = pvP;
    col.subtract(vP);
    col.times(1./e);
    testA.setColumn(col, 3);
    pertF = reducedvF;
    pertF.at(4) += e;
    this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, s, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
    col = pvP;
    col.subtract(vP);
    col.times(1./e);
    testA.setColumn(col, 4);
    pertF = reducedvF;
    pertF.at(5) += e;
    this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, s, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
    col = pvP;
    col.subtract(vP);
    col.times(1./e);
    testA.setColumn(col, 5);
    pertF = reducedvF;
    pertF.at(6) += e;
    this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, s, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
    col = pvP;
    col.subtract(vP);
    col.times(1./e);
    testA.setColumn(col, 6);
    pertF = reducedvF;
    pertF.at(7) += e;
    this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, s, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
    col = pvP;
    col.subtract(vP);
    col.times(1./e);
    testA.setColumn(col, 7);
    pertF = reducedvF;
    pertF.at(8) += e;
    this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, s, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
    col = pvP;
    col.subtract(vP);
    col.times(1./e);
    testA.setColumn(col, 8);
    pertF = reducedvF;
    pertF.at(9) += e;
    this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, s, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
    col = pvP;
    col.subtract(vP);
    col.times(1./e);
    testA.setColumn(col, 9);
    pertF = reducedvF;    
    this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, s, M, gp, reducedvF, reducedMicromorphicVar, micromorphicVarGrad, tStep);


    FloatMatrix m3d = testA;
    testA.resize(5,5);
    
    testA.zero();
    testA.at(1, 1) = m3d.at(1, 1);
    testA.at(1, 2) = m3d.at(1, 2);
    testA.at(1, 4) = m3d.at(1, 6);
    testA.at(1, 5) = m3d.at(1, 9);
    
    testA.at(2, 1) = m3d.at(2, 1);
    testA.at(2, 2) = m3d.at(2, 2);
    testA.at(2, 4) = m3d.at(2, 6);
    testA.at(2, 5) = m3d.at(2, 9);
    
    testA.at(3, 1) = m3d.at(3, 1);
    testA.at(3, 2) = m3d.at(3, 2);
    testA.at(3, 4) = m3d.at(3, 6);
    testA.at(3, 5) = m3d.at(3, 9);
    
    testA.at(4, 1) = m3d.at(6, 1);
    testA.at(4, 2) = m3d.at(6, 2);
    testA.at(4, 4) = m3d.at(6, 6);
    testA.at(4, 5) = m3d.at(6, 9);
    
    testA.at(5, 1) = m3d.at(9, 1);
    testA.at(5, 2) = m3d.at(9, 2);
    testA.at(5, 4) = m3d.at(9, 6);
    testA.at(5, 5) = m3d.at(9, 9);
    */
    
    
  
}

void
GradientF_PolyconvexMaterial :: giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{


    GradientF_PolyconvexMaterialStatus *status = static_cast< GradientF_PolyconvexMaterialStatus* >( this->giveStatus(gp) );

    
   if(gp->giveMaterialMode() == _PlaneStrain) {
     answer.resize(5,5);
     answer.beUnitMatrix();
     answer.times(-Hk);
     
   } else {
     answer.resize(9,9);
     answer.beUnitMatrix();
     answer.times(-Hk);
   }
   

   /*
   FloatMatrix testA(9,9);
   FloatMatrix F;
   F.beMatrixForm(status->giveTempFVector());
   FloatArray vP, pvP, M, pertF, ps, s, reducedvF, reducedMicromorphicVar(9), micromorphicVarGrad(1);
   FloatArray col;
   double e = 1.e-6;
   reducedvF.beVectorForm(F);
   reducedMicromorphicVar = status->giveTempMicromorphicVar();
   this->giveFiniteStrainGeneralizedStressVectors_3d(vP, s, M, gp, reducedvF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
   pertF = reducedvF;
   pertF.at(1) += e;
   this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, ps, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
   col = ps;
   col.subtract(s);
   col.times(1./e);
   testA.setColumn(col, 1);
   pertF = reducedvF;
   pertF.at(2) += e;
   this->giveFiniteStrainGeneralizedStressVectors_3d(pvP,ps, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
   col = ps;
   col.subtract(s);
   col.times(1./e);
   testA.setColumn(col, 2);
   pertF = reducedvF;
   pertF.at(3) += e;
   this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, ps, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
   col = ps;
   col.subtract(s);
   col.times(1./e);
   testA.setColumn(col, 3);
   pertF = reducedvF;
   pertF.at(4) += e;
   this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, ps, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
   col = ps;
   col.subtract(s);
   col.times(1./e);
   testA.setColumn(col, 4);
   pertF = reducedvF;
   pertF.at(5) += e;
   this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, ps, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
   col = ps;
   col.subtract(s);
   col.times(1./e);
   testA.setColumn(col, 5);
   pertF = reducedvF;
   pertF.at(6) += e;
   this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, ps, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
   col = ps;
   col.subtract(s);
   col.times(1./e);
   testA.setColumn(col, 6);
   pertF = reducedvF;
   pertF.at(7) += e;
   this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, ps, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
   col = ps;
   col.subtract(s);
   col.times(1./e);
   testA.setColumn(col, 7);
   pertF = reducedvF;
   pertF.at(8) += e;
   this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, ps, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
   col = ps;
   col.subtract(s);
   col.times(1./e);
   testA.setColumn(col, 8);
   pertF = reducedvF;
   pertF.at(9) += e;
   this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, ps, M, gp, pertF, reducedMicromorphicVar, micromorphicVarGrad, tStep);
   col = ps;
   col.subtract(s);
   col.times(1./e);
   testA.setColumn(col, 9);
   pertF = reducedvF;    
   this->giveFiniteStrainGeneralizedStressVectors_3d(pvP, ps, M, gp, reducedvF, reducedMicromorphicVar, micromorphicVarGrad, tStep);


   FloatMatrix m3d = testA;
   testA.resize(5,5);
    
   testA.zero();
   testA.at(1, 1) = m3d.at(1, 1);
   testA.at(1, 2) = m3d.at(1, 2);
   testA.at(1, 4) = m3d.at(1, 6);
   testA.at(1, 5) = m3d.at(1, 9);
    
   testA.at(2, 1) = m3d.at(2, 1);
   testA.at(2, 2) = m3d.at(2, 2);
   testA.at(2, 4) = m3d.at(2, 6);
   testA.at(2, 5) = m3d.at(2, 9);
    
   testA.at(3, 1) = m3d.at(3, 1);
   testA.at(3, 2) = m3d.at(3, 2);
   testA.at(3, 4) = m3d.at(3, 6);
   testA.at(3, 5) = m3d.at(3, 9);
    
   testA.at(4, 1) = m3d.at(6, 1);
   testA.at(4, 2) = m3d.at(6, 2);
   testA.at(4, 4) = m3d.at(6, 6);
   testA.at(4, 5) = m3d.at(6, 9);
    
   testA.at(5, 1) = m3d.at(9, 1);
   testA.at(5, 2) = m3d.at(9, 2);
   testA.at(5, 4) = m3d.at(9, 6);
   testA.at(5, 5) = m3d.at(9, 9);
   

   */



   
}


void
GradientF_PolyconvexMaterial :: giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

      GradientF_PolyconvexMaterialStatus *status = static_cast< GradientF_PolyconvexMaterialStatus* >( this->giveStatus(gp) );

   if(gp->giveMaterialMode() == _PlaneStrain) {
     answer.resize(5,5);
     answer.beUnitMatrix();
     answer.times(-Hk);
     
   } else {
     answer.resize(9,9);
     answer.beUnitMatrix();
     answer.times(-Hk);
   }


}
 
  

void
GradientF_PolyconvexMaterial :: giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    if(gp->giveMaterialMode() == _PlaneStrain) {
     answer.resize(5,5);
     answer.beUnitMatrix();
     answer.times(Hk);
     
   } else {
     answer.resize(9,9);
     answer.beUnitMatrix();
     answer.times(Hk);
   }

}


void
GradientF_PolyconvexMaterial :: giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  if (matMode == _PlaneStrain) {
    answer.resize(10,10);
  } else {
    answer.resize(27,27);
  } 
  answer.beUnitMatrix();
  answer.times(Ak);

}



void
GradientF_PolyconvexMaterial :: giveFiniteStrainGeneralizedStressVectors(FloatArray &vP, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &vF, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{
    ///@todo Move this to StructuralCrossSection ?
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
      this->giveFiniteStrainGeneralizedStressVectors_3d(vP, s, M, gp, vF, micromorphicVar, micromorphicVarGrad, tStep);
    } else if ( mode == _PlaneStrain ) {
      this->giveFiniteStrainGeneralizedStressVectors_PlaneStrain(vP, s, M, gp, vF, micromorphicVar, micromorphicVarGrad, tStep);
    } else {
      OOFEM_ERROR("Unknown material mode for the gradient polyconvex formulation");
    }
}
  


void
GradientF_PolyconvexMaterial :: giveFiniteStrainGeneralizedStressVectors_3d(FloatArray &vP, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &vF, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{
  
    GradientF_PolyconvexMaterialStatus *status = static_cast< GradientF_PolyconvexMaterialStatus* >( this->giveStatus(gp) );
    //deformation gradient, its inverse, cofactor, and determinant
    FloatArray vCofF, vInvF;
    FloatMatrix F, invF, cofF;
    F.beMatrixForm(vF);
    double J = F.giveDeterminant();
    // relative stress for cofactor
    s = vF;
    s.subtract(micromorphicVar);
    s.times(-Hk);
    // higher order stress
    M = micromorphicVarGrad;
    M.times(Ak);   
    // first PK stress
    FloatArray vPm(9);
    // double well material model    
    double normC_tC1, normC_tC2;
    FloatArray arb1, arb2, vC_tC1, vC_tC2;
    FloatMatrix C, dCdF, C_tC1, C_tC2;
    C.beTProductOf(F,F);
    this->compute_dC_dF(dCdF,vF);

    C_tC1 = C;
    C_tC1.subtract(tC1);

    C_tC2 = C;
    C_tC2.subtract(tC2);
    

    normC_tC1 = C_tC1.computeFrobeniusNorm();
    normC_tC2 = C_tC2.computeFrobeniusNorm();

    vC_tC1.beVectorForm(C_tC1);
    vC_tC2.beVectorForm(C_tC2);   


    arb1.beTProductOf(dCdF, vC_tC1);
    arb1.times(2. * normC_tC2 * normC_tC2);
      
    arb2.beTProductOf(dCdF, vC_tC2);
    arb2.times(2. * normC_tC1 * normC_tC1);

    
    vP = arb1;
    vP.add(arb2);
    vP.times(alpha);
    
    //subtract the relative stress
    vP.subtract(s);

    /// part with J^-1
    //vInvF.times(1./J);
    //vP.subtract(vInvF);

    
    status->letTempMicromorphicVarBe(micromorphicVar);
    status->letTempMicromorphicVarGradBe(micromorphicVarGrad);

    status->letTempPVectorBe(vP);
    status->letTempFVectorBe(vF);
    status->letTempMicromorphicStressBe(s);
    status->letTempMicromorphicStressGradBe(M); 
      
}

  


void
GradientF_PolyconvexMaterial :: giveFiniteStrainGeneralizedStressVectors_PlaneStrain(FloatArray &vP, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &reducedvF, const FloatArray &reducedMicromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{

  FloatArray vF, vMV, fullvP, fullS;
  StructuralMaterial :: giveFullVectorFormF(vF, reducedvF, _PlaneStrain);
  StructuralMaterial :: giveFullVectorFormF(vMV, reducedMicromorphicVar, _PlaneStrain);
  this->giveFiniteStrainGeneralizedStressVectors_3d(fullvP, fullS, M, gp, vF, vMV, micromorphicVarGrad, tStep) ;
  StructuralMaterial :: giveReducedVectorForm(vP, fullvP, _PlaneStrain);
  StructuralMaterial :: giveReducedVectorForm(s, fullS, _PlaneStrain);

      
}
  

void
GradientF_PolyconvexMaterial :: compute_dC_dF(FloatMatrix &dCdF,const FloatArray &vF)
{

  dCdF.resize(9,9);
  FloatArray delta = {1, 1, 1, 0, 0, 0, 0, 0, 0};
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int m = 1; m <=3; m++) {
	for (int n = 1; n<=3; n++) {
	  dCdF.at(giveVI(i,j),giveVI(m,n)) = vF.at(giveVI(m,j)) * delta.at(giveVI(i,n)) + vF.at(giveVI(m,i)) * delta.at(giveVI(j,n));
	}
      }
    }
  }
  
}



  
  

IRResultType
GradientF_PolyconvexMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    IsotropicLinearElasticMaterial :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, Hk, _IFT_MicromorphicMaterialExtensionInterface_Hk);
    IR_GIVE_FIELD(ir, Ak, _IFT_MicromorphicMaterialExtensionInterface_Ak);

    IR_GIVE_FIELD(ir, eps, _IFT_GradientF_PolyconvexMaterial_eps);
    IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_GradientF_PolyconvexMaterial_alpha );

    
    /*    tC1 = {{(1. + eps) * (1. + eps), 0, 0},{0, 1. / (1. + eps) / (1. + eps), 0}, {0, 0, 1}};
    tC2 = {{1. / (1. + eps) / (1. + eps), 0, 0},{0, (1. + eps)*(1. + eps), 0}, {0, 0, 1}};
    */

    
    tC1 = {{ 1.0, eps, 0},{eps, 1.+eps*eps, 0}, {0, 0, 1}};
    tC2 = {{1.0, -eps, 0},{-eps, 1.+eps*eps, 0}, {0, 0, 1}};
    

    /*    tC1 = {{ 1.0, eps, 0},{eps, 1., 0}, {0, 0, 1.}};
    tC2 = {{1.0, -eps, 0},{-eps, 1., 0}, {0, 0, 1.}};
    */
    /*
    tC1 = {{eps*eps, 0, 0},{0, 1., 0}, {0, 0, 1/eps/eps}};
    tC2 = {{1/eps/eps, 0, 0},{0, 1., 0}, {0, 0, eps*eps}};
    */
    
    
    return IRRT_OK;
}

int
GradientF_PolyconvexMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    GradientF_PolyconvexMaterialStatus *status = static_cast< GradientF_PolyconvexMaterialStatus * >( this->giveStatus(gp) );
    

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
    } else if(type == IST_MaxEquivalentStrainLevel) {
      FloatMatrix F, C, C_tC1, C_tC2;
      F.beMatrixForm(status->giveTempFVector());
      C.beTProductOf(F,F);
      C_tC1 = C;
      C_tC2 = C;
      C_tC1.subtract(tC1);
      C_tC2.subtract(tC2);
      double normC_tC1, normC_tC2;
      normC_tC1 = C_tC1.computeFrobeniusNorm();
      normC_tC2 = C_tC2.computeFrobeniusNorm();
      answer.resize(1);
      answer.at(1) = normC_tC1 * normC_tC1/(normC_tC1 * normC_tC1 + normC_tC2 * normC_tC2);
    } else {
      OOFEM_ERROR("Unknown InternalStateType");
    }
    return 1;
}
    

  GradientF_PolyconvexMaterialStatus :: GradientF_PolyconvexMaterialStatus(int n, Domain *d, GaussPoint *g, bool sym) : MicromorphicMaterialStatus(n, d, g, sym)
{
    micromorphicVar.resize(9);
    micromorphicVar.at(1) = micromorphicVar.at(2) = micromorphicVar.at(3) = 1.;
    micromorphicStress.resize(9);
    micromorphicVarGrad.resize(18);
    micromorphicStressGrad.resize(18);


    tempMicromorphicVar = micromorphicVar;
    tempMicromorphicVarGrad = micromorphicVarGrad;
    tempMicromorphicStress = micromorphicStress;
    tempMicromorphicStressGrad = micromorphicStressGrad;
}


} // end namespace oofem
