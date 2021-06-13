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

#include "martensitemat.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {
REGISTER_Material(MartensiteMaterial);

MartensiteMaterial :: MartensiteMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{ }


void
MartensiteMaterial :: giveSecondPKStressVector_3d(FloatMatrix &answer, GaussPoint *gp, const FloatMatrix &C, TimeStep *tStep)
// returns 6 components of the second Piola Kirchhoff stress corresponding to the given deformation gradient
{

    MartensiteMaterialStatus *status = static_cast< MartensiteMaterialStatus * >( this->giveStatus(gp) );

    if(!status->isDissipationVectorInitialized()) {
      FloatArray array_wellDistance,z;
      array_wellDistance.resize(M+1);
      z.resize(M+1);
      FloatArray vC,Ci, dC;
      vC = {1,1,1,0,0,0};
      double sum_norm = 0;
      for (int i = 0; i<= M; i++) {
	dC = vC;
	Ci = vector_Ci.at(i);
	dC.subtract(Ci);
	double norm = dC.computeNorm();
	array_wellDistance.at(i+1) = norm  -  arrayEpsilon.at(i+1);
	sum_norm += norm;
      }

	for (int i = 0; i<= M; i++) {
	  if(sum_norm == 0) {
	    z.at(i+1) = 1.  / M;
	  } else {
	    z.at(i+1) = (1. - array_wellDistance.at(i+1) / sum_norm) / M;
	  }
	}
    
      status->setDissipationVector(z);
      status->setDissipationVectorInitialized();
    }
  

    int index = 1;
    double Energy = 0;

    FloatArray vEi, vS;
    FloatMatrix iFi, Ei, I(3,3), D;
    I.beUnitMatrix();
    // Elastic part
    for(int i = 0; i <= M; i++) {
      FloatMatrix iFiF;
      iFi = vector_iF.at(i);
      this->give3dMaterialStiffnessMatrix(D, TangentStiffness, gp, tStep);
      this->compute_sym_dyadic_product_reduced(iFiF, iFi, iFi);
       
      FloatMatrix junk;
      junk.beTProductOf(iFi,C);
      Ei.beProductOf(junk, iFi);
      Ei.subtract(I);
      Ei.times(0.5);
      vEi.beSymVectorForm(Ei);
      // energy
      FloatArray E;
      E.beProductOf(D,vEi);
      double Ei = E.dotProduct(vEi);
      FloatArray vSi;
      vSi.beProductOf(iFiF, E);    

      if(i == 0) {
	Energy = Ei;
	index = 0;
	vS = vSi;
      } else {
	if(Ei < Energy) {
	  Energy = Ei;
	  index = i;
	  vS = vSi;
	}
      }    
    }
    status->setActiveWellIndex(index);
    //Dissipation
    double dissipation = 0;
    double sum_distance = 0;

    FloatArray array_wellDistance;
    array_wellDistance.resize(M+1);
    FloatArray z, Ci, vC, dC;
    vC.beSymVectorForm(C);
    z = status->giveDissipationVector();
    for(int i = 0; i <= M; i++) {
      dC = vC;
      Ci = vector_Ci.at(i);
      double epsilon_i = arrayEpsilon.at(i+1);
      dC.subtract(Ci);
      double wellDistance = dC.computeNorm() - epsilon_i;
      if(wellDistance < 0) {
	wellDistance = 0;
      }
      array_wellDistance.at(i+1) = wellDistance; 
      sum_distance += wellDistance;
    }
    if(sum_distance != 0) {
      array_wellDistance.times(1./sum_distance);
    }
    status->setArrayWellDistance(array_wellDistance);
    FloatArray lambda;
    lambda.resize(4);
    for(int i = 0; i <= M; i++) {
      lambda.at(i+1) = (1 - array_wellDistance.at(i+1)) / M;
      dissipation += (lambda.at(i+1) - z.at(i+1)) * (lambda.at(i+1) - z.at(i+1)) ;
    }
    dissipation = sqrt(dissipation);
    double diss = status->giveDissipation();
    status->setTempDissipation(diss + dissipation);
    status->setTempDissipationVector(lambda);
    FloatArray vSd;
    if(dissipation > 0) {
      for(int i = 0; i <= M; i++) { 
	if(array_wellDistance.at(i+1) > 0) {
	  FloatArray sum_dDistance_dC, dDistance_dC, vSdi;
	  for(int j = 0; j <= M; j++) {
	    if(array_wellDistance.at(j+1) > 0) {
	      FloatArray dC, Cj;
	      dC = vC;
	      Cj = vector_Ci.at(j);
	      dC.subtract(Cj);
	      sum_dDistance_dC.add(dC);
	    }
	  }
	  sum_dDistance_dC.times(2. * array_wellDistance.at(i+1) / sum_distance / sum_distance / M);

	  FloatArray Ci;
	  Ci = vector_Ci.at(i);
	  dDistance_dC = vC;
	  dDistance_dC.subtract(Ci);
	  dDistance_dC.times( - 2./sum_distance/M);
	  vSdi = sum_dDistance_dC;
	  vSdi.add(dDistance_dC);
	  vSdi.times(lambda.at(i+1) - z.at(i+1));
	  vSd.add(vSdi);
	}
      }
    } else {
      vSd.resize(6);
    }
    vS.add(vSd);
    answer.beMatrixForm(vS);

    
}


  
void
MartensiteMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{

    FloatArray vS;
    FloatMatrix F, C, S, P;
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    //
    C.beTProductOf(F, F);    
    this->giveSecondPKStressVector_3d(S, gp, C, tStep);
    P.beProductOf(F,S);
    vS.beSymVectorForm(S);
    answer.beVectorForm(P);
    // update gp
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    status->letTempFVectorBe(vF);
    status->letTempStressVectorBe(vS);
    status->letTempPVectorBe(answer);
}



void
MartensiteMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  answer.resize(6,6);
  answer.at(1,1) = answer.at(2,2) = answer.at(3,3) = D11;
  answer.at(1,2) = answer.at(1,3) = D12;
  answer.at(2,1) = answer.at(2,3) = D12;
  answer.at(3,1) = answer.at(3,2) = D12;
  answer.at(4,4) = answer.at(5,5) = answer.at(6,6) = D44;  

}

  
void
MartensiteMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    MartensiteMaterialStatus *status = static_cast< MartensiteMaterialStatus * >( this->giveStatus(gp) );

    FloatMatrix iFi, Di, dSdE(6,6);
    double i = status->giveActiveWellIndex();
    iFi = vector_iF.at(i);
    this->give3dMaterialStiffnessMatrix(Di, TangentStiffness, gp, tStep);
    FloatMatrix Qi, QD, QDQ;
    this->compute_sym_dyadic_product_reduced(Qi, iFi, iFi);
    QD.beProductOf(Qi,Di);
    QDQ.beProductTOf(QD,Qi);

    dSdE.add(QDQ);

    FloatArray vF, vC, vCi;
    FloatMatrix F, C;
    vF = status->giveTempFVector();
    F.beMatrixForm(vF);
    //
    C.beTProductOf(F, F);
    vC.beSymVectorForm(C);

    
    double diss = status->giveDissipation();
    double temp_diss =     status->giveTempDissipation();
    if(temp_diss - diss > 0) {
      double sum_distance = 0;
      FloatArray array_wellDistance;
      array_wellDistance = status->giveArrayWellDistance();
      for(int i = 0; i <= M; i++) { 
	if(array_wellDistance.at(i+1) > 0) {
	  FloatArray sum_dDistance_dC, dDistance_dC, vSdi;
	  for(int j = 0; j <= M; j++) {
	    FloatArray dC, Cj;
	    if(array_wellDistance.at(j+1) > 0) {
	      dC = vC;
	      Cj = vector_Ci.at(j);
	      dC.subtract(Cj);
	      sum_dDistance_dC.add(dC);
	      sum_distance += array_wellDistance.at(j+1);
	    } 
	  }
	  FloatArray dC, Cj;
	  dC = vC;
	  vCi = vector_Ci.at(i);
	  dC.subtract(vCi);
	  /// first derivative of lambda	     
	  FloatArray dLambda(6);
	  dLambda.add(2. * array_wellDistance.at(i+1) / sum_distance/sum_distance / M, sum_dDistance_dC);
	  dLambda.add(-2./ sum_distance/M, dC);
	  /// second derivative of lambda
	  FloatMatrix d2Lambda(6,6), Is, dC_dSumD,dSumD_dC,dSumD_dSumD  ;
	  FloatArray delta;
	  delta = {1,1,1,0,0,0};
	  this->compute_dyadic_product_reduced( Is, delta, delta);
	  this->compute_dyadic_product_reduced( dC_dSumD, dC, sum_dDistance_dC);
	  this->compute_dyadic_product_reduced( dSumD_dC, sum_dDistance_dC, dC);
	  this->compute_dyadic_product_reduced( dSumD_dSumD, sum_dDistance_dC, sum_dDistance_dC);
	  d2Lambda.add(-2./sum_distance/M, Is);
	  d2Lambda.add(2. / sum_distance / sum_distance / M, dC_dSumD);
	  d2Lambda.add(2. / sum_distance / sum_distance / M, dSumD_dC);
	  d2Lambda.add(2. / sum_distance / sum_distance / M, dSumD_dC);
	  d2Lambda.add(2. * array_wellDistance.at(i+1)/sum_distance / sum_distance, Is);
	  d2Lambda.add(-2. * array_wellDistance.at(i+1)/sum_distance / sum_distance / sum_distance, dSumD_dSumD);  
	}
      }
    } 
    
    // transform the second elasticity tensor to the first elasticity tensor    
    this->give_dPdF_from(dSdE, answer, gp, _3dMat);     

}




MaterialStatus *
MartensiteMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new MartensiteMaterialStatus(1, this->giveDomain(), gp);
}


IRResultType
MartensiteMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro


    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;


    double eta1, eta2;
    IR_GIVE_FIELD(ir, eta1, _IFT_MartensiteMaterial_eta1);
    IR_GIVE_FIELD(ir, eta2, _IFT_MartensiteMaterial_eta2);

    IR_GIVE_FIELD(ir, M, _IFT_MartensiteMaterial_M);


    FloatMatrix iF0, iF1, iF2, iF3;
      
    iF0 = {{1.,0,0},{0,1.,0},{0,0,1.}};
    iF1 = {{1./eta2,0,0},{0,1./eta1,0},{0,0,1./eta1}};
    iF2 = {{1./eta1,0,0},{0,1./eta2,0},{0,0,1./eta1}};
    iF3 = {{1./eta1,0,0},{0,1./eta1,0},{0,0,1./eta2}};

    vector_iF.push_back(iF0);
    vector_iF.push_back(iF1);
    vector_iF.push_back(iF2);
    vector_iF.push_back(iF3);


    FloatMatrix F0, F1, F2, F3, C0, C1, C2, C3;
      
    F0 = {{1.,0,0},{0,1.,0},{0,0,1.}};
    F1 = {{eta2,0,0},{0,eta1,0},{0,0,eta1}};
    F2 = {{eta1,0,0},{0,eta2,0},{0,0,eta1}};
    F3 = {{eta1,0,0},{0,eta1,0},{0,0,eta2}};
    C0.beTProductOf(F0,F0);
    C1.beTProductOf(F1,F1);
    C2.beTProductOf(F2,F2);
    C3.beTProductOf(F2,F2);

    FloatArray vC0, vC1, vC2, vC3;
    vC0.beSymVectorForm(C0);
    vC1.beSymVectorForm(C1);
    vC2.beSymVectorForm(C2);
    vC3.beSymVectorForm(C3);
    
    vector_Ci.push_back(vC0);
    vector_Ci.push_back(vC1);
    vector_Ci.push_back(vC2);
    vector_Ci.push_back(vC3);

    

    
    
    
    IR_GIVE_FIELD(ir, D11, _IFT_MartensiteMaterial_D11);
    IR_GIVE_FIELD(ir, D12, _IFT_MartensiteMaterial_D12);
    IR_GIVE_FIELD(ir, D44, _IFT_MartensiteMaterial_D44);

    arrayEpsilon.resize(4);
    arrayEpsilon = {10.07, 10.07, 10.07, 10.07};
    

    
    return IRRT_OK;
}




int
MartensiteMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    MartensiteMaterialStatus *status = static_cast< MartensiteMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
	FloatArray z;
	z = status->giveDissipationVector();
        answer.at(1) = z.at(1);
	answer.at(2) = z.at(2);
	answer.at(3) = z.at(3);
	answer.at(4) = z.at(4);
        return 1;
    } else if ( type == IST_DissWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveDissipation();
        return 1;
    } else if (type == IST_MaxEquivalentStrainLevel) {
        answer.resize(1);
	answer.at(1) = status->giveActiveWellIndex();
	return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }

    return 1; // to make the compiler happy
}









MartensiteMaterialStatus :: MartensiteMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    dissipation = 0;
    tempDissipation = 0;
    dissipationVector.clear();
    tempDissipationVector.clear();
    array_wellDistance.clear();

}


void
MartensiteMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempDissipation = this->dissipation;
    this->tempDissipationVector = this->dissipationVector;
}


void
MartensiteMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->dissipation = this->tempDissipation;
    this->dissipationVector = this->tempDissipationVector;

}



} // end namespace oofem
