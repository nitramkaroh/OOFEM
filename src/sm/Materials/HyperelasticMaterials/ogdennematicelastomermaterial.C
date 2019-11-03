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

#include "ogdennematicelastomermaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {
REGISTER_Material(OgdenNematicElastomerMaterial);

OgdenNematicElastomerMaterial :: OgdenNematicElastomerMaterial(int n, Domain *d) : OgdenMaterial(n, d)
{ }


 

void
OgdenNematicElastomerMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{

    FloatArray eVals, vS;
    FloatMatrix P, F, S, C, eVecs;
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    //
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    // phase of nematic elastomer
    int nep = 0;
    if(!qcEnvelop) {    
      this->giveSecondPKStressVector_3d(S, C, eVals, eVecs, tStep);
      // add the terms according to nematic elastomer
    } else {
      int minIndex = eVals.giveIndexMinElem();
      int maxIndex = eVals.giveIndexMaxElem();    

      if(1./sqrt(eVals.at(minIndex)) <= pow(a,1./6.)) {
	this->giveSecondPKStressVector_3d(S, C, eVals, eVecs, tStep);
      } else if(1./sqrt(eVals.at(minIndex)) >= pow(a, -0.5) * eVals.at(maxIndex)) {
	// solid phase
	this->giveSecondPKStressVector_3d(S, C, eVals, eVecs, tStep);
	nep = 1;
      } else {
	// smectic phase
	this->giveSecondPKStressVectorSmectic_3d(S, C, eVals, eVecs, tStep);
	nep = 2;
      }
    }
    
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
OgdenNematicElastomerMaterial :: giveSecondPKStressVector_3d(FloatMatrix &answer, const FloatMatrix &C, const FloatArray &eVals, const FloatMatrix &eVecs,TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    double J, lnJ;
    FloatMatrix invC, I1aInvC;
    //
    J = sqrt(C.giveDeterminant());
    //compute inverse of the right Cauchy-Green tensor(C)
    invC.beInverseOf(C);
    int maxIndex = eVals.giveIndexMaxElem();
    FloatArray a_coeff(3);
   
    
    for (int i = 1; i <= N; i++) {
      double l1a, l2a, l3a, Ja;
      FloatMatrix Si;
      // compute a coefficinets
      for (int j = 1; j <= 3; j++) {
	if(j == maxIndex) {
	  a_coeff.at(j) = pow(a, -alpha.at(i)/3.);
	} else {
	  a_coeff.at(j) = pow(a, alpha.at(i)/6.);
	}
      }

      l1a = a_coeff.at(1) * pow(eVals.at(1),alpha.at(i)/2.);
      l2a = a_coeff.at(2) * pow(eVals.at(2),alpha.at(i)/2.);
      l3a = a_coeff.at(3) * pow(eVals.at(3),alpha.at(i)/2.);
      Ja = pow(J, -alpha.at(i) / 3.);
      double I1a = (l1a + l2a + l3a) / 3.;
      I1aInvC = invC;
      I1aInvC.times(I1a);

      this->computeMatrixPower(Si, eVals, eVecs, ( alpha.at(i)-2. ) / 2., a_coeff );
      
      Si.subtract(I1aInvC);
      Si.times(mu.at(i) * Ja );
      answer.add(Si);
    }
    
    // compute jacobian and its logarithm
    lnJ = log(J);
    answer.add(K * lnJ, invC);

}


void
OgdenNematicElastomerMaterial :: giveSecondPKStressVectorSmectic_3d(FloatMatrix &answer, const FloatMatrix &C, const FloatArray &eVals, const FloatMatrix &eVecs,TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    
    double J, lnJ;
    FloatMatrix invC, I1aInvC, I1bInvC;
    //
    J = sqrt(C.giveDeterminant());
    //compute inverse of the right Cauchy-Green tensor(C)
    invC.beInverseOf(C);
    int minIndex = eVals.giveIndexMinElem();
    FloatArray a_coeff(3), b_coeff(3);
    FloatArray beta(alpha);
    beta.times(-0.5);
    
    for (int i = 1; i <= N; i++) {
      double l1a, l2a, l3a, Ja;
      double l1b, l2b, l3b, Jb;
      
      FloatMatrix Sia, Sib;
      // compute a coefficinets
      for (int j = 1; j <= 3; j++) {
	if(j == minIndex) {
	  a_coeff.at(j) = pow(a, alpha.at(i) / 6.);
	  b_coeff.at(j) = pow(a, beta.at(i) / 6.);
	} else {
	  a_coeff.at(j) = 0;
	  b_coeff.at(j) = 0;
	}
      }
      
      //
      l1a = a_coeff.at(1) * pow(eVals.at(1),alpha.at(i)/2.);
      l2a = a_coeff.at(2) * pow(eVals.at(2),alpha.at(i)/2.);
      l3a = a_coeff.at(3) * pow(eVals.at(3),alpha.at(i)/2.);
      Ja = pow(J, -alpha.at(i) / 3.);
      double I1a = (l1a + l2a + l3a) / 3.;
      I1aInvC = invC;
      I1aInvC.times(I1a);
      //
      this->computeMatrixPower(Sia, eVals, eVecs, ( alpha.at(i)-2. ) / 2., a_coeff );
      //
      Sia.subtract(I1aInvC);
      Sia.times(mu.at(i) * Ja );
      answer.add(Sia);
      /////////////////////////////////
      l1b = b_coeff.at(1) * pow(eVals.at(1),beta.at(i)/2.);
      l2b = b_coeff.at(2) * pow(eVals.at(2),beta.at(i)/2.);
      l3b = b_coeff.at(3) * pow(eVals.at(3),beta.at(i)/2.);
      Jb = pow(J, -beta.at(i) / 3.);
      double I1b = (l1b + l2b + l3b) / 3.;
      I1bInvC = invC;
      I1bInvC.times(I1b);

      this->computeMatrixPower(Sib, eVals, eVecs, ( beta.at(i)-2. ) / 2., b_coeff );
      
      Sib.subtract(I1bInvC);
      Sib.times(beta.at(i) / alpha.at(i) * mu.at(i) * Jb );
      answer.add(Sib);      
    }
    
    // compute jacobian and its logarithm
    lnJ = log(J);
    answer.add(K * lnJ, invC);



    
}

  


void
OgdenNematicElastomerMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    FloatArray vF, eVals;
    FloatMatrix P, F, S, C, invC, iCiC, dSdE;
    FloatMatrix eVecs;
    
    //deformation gradient from the status
    vF = status->giveTempFVector();
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    //
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);

    if(!qcEnvelop) {    
      this->give3dMaterialStiffnessMatrix_dSdE(dSdE, mode, gp, tStep);
      // add the terms according to nematic elastomer
    } else {
      int minIndex = eVals.giveIndexMinElem();
      int maxIndex = eVals.giveIndexMaxElem();    
      if(1./sqrt(eVals.at(minIndex)) <= pow(a,1./6.)) {
	this->give3dMaterialStiffnessMatrix_dSdE(dSdE, mode, gp, tStep);
      } else if(1./sqrt(eVals.at(minIndex)) >= pow(a, -0.5) * eVals.at(maxIndex)) {
	this->give3dMaterialStiffnessMatrix_dSdE(dSdE, mode, gp, tStep);
      } else {
	this->give3dMaterialStiffnessMatrixSmectic_dSdE(dSdE, mode, gp, tStep);
      }
    }
   
    this->give_dPdF_from(dSdE, answer, gp, _3dMat);
    
}




void
OgdenNematicElastomerMaterial :: give3dMaterialStiffnessMatrix_dSdE(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double J, lnJ;
    FloatArray vF, eVals;
    FloatMatrix P, F, S, C, invC, iCiC;
    FloatMatrix eVecs;
    FloatArray a_coeff(3);
    //deformation gradient from the status
    vF = status->giveTempFVector();
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    J = F.giveDeterminant();
    lnJ = log(J);
    //
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    int maxIndex = eVals.giveIndexMaxElem();   
    //compute inverse of the right Cauchy-Green tensor(C)
    invC.beInverseOf(C);
    //compute symetric dyadic product 
    this->compute_sym_dyadic_product_reduced(iCiC, invC, invC);
    //    
    for (int i = 1; i <= N; i++) {
      double l1a, l2a, l3a, Ja;
      //
      for (int j = 1; j <= 3; j++) {
	if(j == maxIndex) {
	  a_coeff.at(j) = pow(a, -alpha.at(i)/3.);
	} else {
	  a_coeff.at(j) = pow(a, alpha.at(i)/6.);
	}
      }
      
      Ja = pow(J, -alpha.at(i) / 3.);
      l1a = a_coeff.at(1) * pow(eVals.at(1),alpha.at(i)/2.);
      l2a = a_coeff.at(2) * pow(eVals.at(2),alpha.at(i)/2.);
      l3a = a_coeff.at(3) * pow(eVals.at(3),alpha.at(i)/2.);
      //
      double I1a = (l1a + l2a + l3a) / 3.;
      //
      FloatMatrix I1aInvC;
      I1aInvC = invC;
      I1aInvC.times(I1a);
      //
      FloatMatrix Sdev, powC, SiC;
      this->computeMatrixPower(powC, eVals, eVecs, ( alpha.at(i)-2. ) / 2., a_coeff );
      Sdev = powC;
      Sdev.subtract(I1aInvC);
      this->compute_dyadic_product_reduced(SiC, Sdev, invC);
      //
      FloatMatrix iCpowC, dCm_dC;
      this->compute_dyadic_product_reduced(iCpowC, invC, powC);
      this->compute_dCm_dC(dCm_dC, (alpha.at(i) - 2.) / 2., eVals, eVecs, a_coeff);
      //
      answer.add( 2. * mu.at(i) * Ja, dCm_dC);
      answer.add( -mu.at(i) * alpha.at(i) * Ja / 3., iCpowC);
      answer.add( 2. * mu.at(i)  * Ja * I1a, iCiC);
      answer.add( -alpha.at(i) * mu.at(i) * Ja / 3., SiC);
    }    
    // volumetric part
    answer.add(- 2. * K * lnJ, iCiC);
    FloatMatrix iCxiC;
    this->compute_dyadic_product_reduced(iCxiC, invC, invC);
    answer.add(K, iCxiC);   
}



void
OgdenNematicElastomerMaterial :: give3dMaterialStiffnessMatrixSmectic_dSdE(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double J, lnJ;
    FloatArray vF, eVals;
    FloatMatrix P, F, S, C, invC, iCiC;
    FloatMatrix eVecs;
    FloatArray a_coeff(3), b_coeff(3);
    //deformation gradient from the status
    vF = status->giveTempFVector();
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    J = F.giveDeterminant();
    lnJ = log(J);
    //
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    //compute inverse of the right Cauchy-Green tensor(C)
    invC.beInverseOf(C);
    //compute symetric dyadic product 
    this->compute_sym_dyadic_product_reduced(iCiC, invC, invC);
    //
    int minIndex = eVals.giveIndexMinElem();
    FloatArray beta(alpha);
    beta.times(-0.5);
    //
    for (int i = 1; i <= N; i++) {
      double l1a, l2a, l3a, Ja;
      double l1b, l2b, l3b, Jb;
      //
      // compute a and b coefficinets
      for (int j = 1; j <= 3; j++) {
	if(j == minIndex) {
	  a_coeff.at(j) = pow(a, alpha.at(i) / 6.);
	  b_coeff.at(j) = pow(a, beta.at(i) / 6.);
	} else {
	  a_coeff.at(j) = 0;
	  b_coeff.at(j) = 0;
	}
      }
      
      Ja = pow(J, -alpha.at(i) / 3.);
      l1a = a_coeff.at(1) * pow(eVals.at(1),alpha.at(i)/2.);
      l2a = a_coeff.at(2) * pow(eVals.at(2),alpha.at(i)/2.);
      l3a = a_coeff.at(3) * pow(eVals.at(3),alpha.at(i)/2.);
      //
      double I1a = (l1a + l2a + l3a) / 3.;
      //
      FloatMatrix I1aInvC;
      I1aInvC = invC;
      I1aInvC.times(I1a);
      //
      FloatMatrix Sdev, powC, SiC;
      this->computeMatrixPower(powC, eVals, eVecs, ( alpha.at(i)-2. ) / 2., a_coeff );
      Sdev = powC;
      Sdev.subtract(I1aInvC);
      this->compute_dyadic_product_reduced(SiC, Sdev, invC);
      //
      FloatMatrix iCpowC, dCm_dC;
      this->compute_dyadic_product_reduced(iCpowC, invC, powC);
      this->compute_dCm_dC(dCm_dC, (alpha.at(i) - 2.) / 2., eVals, eVecs, a_coeff);
      //
      answer.add( 2. * mu.at(i) * Ja, dCm_dC);
      answer.add( -mu.at(i) * alpha.at(i) * Ja / 3., iCpowC);
      answer.add( 2. * mu.at(i)  * Ja * I1a, iCiC);
      answer.add( -alpha.at(i) * mu.at(i) * Ja / 3., SiC);
      /////////////////////////////////////////////////////
      Jb = pow(J, -beta.at(i) / 3.);
      l1b = b_coeff.at(1) * pow(eVals.at(1),beta.at(i)/2.);
      l2b = b_coeff.at(2) * pow(eVals.at(2),beta.at(i)/2.);
      l3b = b_coeff.at(3) * pow(eVals.at(3),beta.at(i)/2.);
      //
      double I1b = (l1b + l2b + l3b) / 3.;
      //
      FloatMatrix I1bInvC;
      I1bInvC = invC;
      I1bInvC.times(I1b);
      //
      this->computeMatrixPower(powC, eVals, eVecs, ( beta.at(i)-2. ) / 2., b_coeff );
      Sdev = powC;
      Sdev.subtract(I1bInvC);
      this->compute_dyadic_product_reduced(SiC, Sdev, invC);
      //
      this->compute_dyadic_product_reduced(iCpowC, invC, powC);
      this->compute_dCm_dC(dCm_dC, (beta.at(i) - 2.) / 2., eVals, eVecs, b_coeff);
      //
      answer.add( 2. * beta.at(i)/ alpha.at(i) * mu.at(i) * Jb, dCm_dC);
      answer.add( -mu.at(i) * beta.at(i)/ alpha.at(i) * beta.at(i) * Jb / 3., iCpowC);
      answer.add( 2. * beta.at(i)/ alpha.at(i) * mu.at(i)  * Jb * I1b, iCiC);
      answer.add( -beta.at(i)/ alpha.at(i) * beta.at(i) * mu.at(i) * Jb / 3., SiC);


    }    
    // volumetric part
    answer.add(- 2. * K * lnJ, iCiC);
    FloatMatrix iCxiC;
    this->compute_dyadic_product_reduced(iCxiC, invC, invC);
    answer.add(K, iCxiC);   

}
  



  


MaterialStatus *
OgdenNematicElastomerMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}


IRResultType
OgdenNematicElastomerMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro


    result = OgdenMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, a, _IFT_OgdenNematicElastomerMaterial_a);
    this->qcEnvelop = ir->hasField(_IFT_OgdenNematicElastomerMaterial_qce);
    
    return IRRT_OK;
}



int
OgdenNematicElastomerMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_NematicElastomerPhase ) {
      StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
      FloatArray vF, eVals;
      FloatMatrix F, C, eVecs;
      //deformation gradient from the status
      vF = status->giveTempFVector();
      //store deformation gradient into matrix
      F.beMatrixForm(vF);
      //
      C.beTProductOf(F, F);
      // compute eigen values and eigen vectors of C
      C.jaco_(eVals, eVecs, 15);
      int minIndex = eVals.giveIndexMinElem();
      int maxIndex = eVals.giveIndexMaxElem();
      answer.resize(1);
      if(1./sqrt(eVals.at(minIndex)) <= pow(a,1./6.)) {
	answer.at(1) = 0;
      } else if(1./sqrt(eVals.at(minIndex)) >= pow(a, -0.5) * eVals.at(maxIndex)) {
	answer.at(1) = 1;
      } else {
	answer.at(1) = 2;
      }
      answer.at(1) = 1./sqrt(eVals.at(minIndex)) - sqrt(a) * eVals.at(maxIndex);
      return 1;
    } else {
      return OgdenMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}




  

} // end namespace oofem
