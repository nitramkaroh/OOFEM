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

#include "../sm/Materials/HyperelasticMaterials/doublewellmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
  REGISTER_Material(DoubleWellMaterial);

DoubleWellMaterial :: DoubleWellMaterial(int n, Domain *d) :StructuralMaterial(n, d)
{
  alpha = 1.e3;
}


void
DoubleWellMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
{
  
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus* >( this->giveStatus(gp) );
    //deformation gradient, its inverse, cofactor, and determinant
    FloatMatrix F;
    
    F.beMatrixForm(vF);

    // double well material model    
    double normC_tC1, normC_tC2;
    FloatArray arb1, arb2, vC_tC1, vC_tC2;
    FloatMatrix C, dCdF, C_tC1, C_tC2;
    C.beTProductOf(F,F);
    this->compute_dC_dF(dCdF,vF);

    FloatMatrix tC1, tC2;
    tC1 = this->givetC1(tStep);
    tC2 = this->givetC2(tStep);
    
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
    
    answer = arb1;
    answer.add(arb2);
    answer.times(alpha);

    status->letTempPVectorBe(answer);
    status->letTempFVectorBe(vF);
      
}


  

// returns the 9x9 tangent stiffness matrix - dP/dF
void
DoubleWellMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus* >( this->giveStatus(gp) );

    answer.resize(9,9);
    answer.zero();
    //deformation gradient, its inverse, cofactor, and determinant
    double normC_tC1, normC_tC2;
    FloatArray vF, vB, delta, vC_tC1, vC_tC2;
    FloatMatrix F, B, C, C_tC1, C_tC2, dCdF;
    delta = {1, 1, 1, 0, 0, 0, 0, 0, 0};
    F.beMatrixForm(status->giveTempFVector());
    B.beProductTOf(F,F);
    vB.beVectorForm(B);
    vF.beVectorForm(F);
    C.beTProductOf(F,F);

    FloatMatrix tC1, tC2;
    tC1 = this->givetC1(tStep);
    tC2 = this->givetC2(tStep);
    
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
	    answer.at(giveVI(i,j),giveVI(p,q)) += alpha * (4. * ( vB.at(giveVI(i,p)) * delta.at(giveVI(j,q)) + vF.at(giveVI(i,q)) *  vF.at(giveVI(p,j))) * (normC_tC1*normC_tC1 + normC_tC2 * normC_tC2) + 4. * vC_tC1.at(giveVI(j,q)) * delta.at(giveVI(i,p)) * (normC_tC2 * normC_tC2) + 4. * vC_tC2.at(giveVI(j,q)) * delta.at(giveVI(i,p)) * (normC_tC1 * normC_tC1) + 4.*vC_tC1_dCdF.at(giveVI(i,j))*vC_tC2_dCdF.at(giveVI(p,q)) + 4.*vC_tC2_dCdF.at(giveVI(i,j))*vC_tC1_dCdF.at(giveVI(p,q)) );
	  }
	}
      }
    }


}

  


  


 

void
DoubleWellMaterial :: compute_dC_dF(FloatMatrix &dCdF,const FloatArray &vF)
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

MaterialStatus *
DoubleWellMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}


  

FloatMatrix&
DoubleWellMaterial :: givetC1(TimeStep *tStep)
 {
   double t = tStep->giveIntrinsicTime();
   double tEps = (t+1) * eps;
     
   tC1_0 = {{ 1.0, tEps, 0},{tEps, 1.+tEps*tEps, 0}, {0, 0, 1}};
   // tC1_0 = {{ 1.0, tEps, 0},{tEps, 1., 0}, {0, 0, 1.}};
   //   tC1_0 = {{(1. + tEps) * (1. + tEps), 0, 0},{0, 1. / (1. + tEps) / (1. + tEps), 0}, {0, 0, 1}};
   //tC1_0 = {{(1. + tEps) * (1. + tEps), tEps, 0},{tEps, (1. + tEps) * (1. + tEps), 0}, {0, 0, 1}};
   
      
   return tC1_0;
 }
  

FloatMatrix &
DoubleWellMaterial :: givetC2(TimeStep *tStep)
 {
   double t = tStep->giveIntrinsicTime();
   double tEps = (t+1) * eps;
     
   tC2_0 = {{ 1.0, -tEps, 0},{-tEps, 1.+tEps*tEps, 0}, {0, 0, 1}};
   //tC2_0 = {{1.0, -tEps, 0},{-tEps, 1., 0}, {0, 0, 1.}};
   //tC2_0 = {{1. / (1. + tEps) / (1. + tEps), 0, 0},{0, (1. + tEps)*(1. + tEps), 0}, {0, 0, 1}};
   //tC1_0 = {{(1. + tEps) * (1. + tEps), -tEps, 0},{-tEps, (1. + tEps) * (1. + tEps), 0}, {0, 0, 1}};

   
   return tC2_0;
 }  
  

IRResultType
DoubleWellMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    //IsotropicLinearElasticMaterial :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, eps, _IFT_DoubleWellMaterial_eps);
    IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_DoubleWellMaterial_alpha );
    
    /*    tC1 = {{(1. + eps) * (1. + eps), 0, 0},{0, 1. / (1. + eps) / (1. + eps), 0}, {0, 0, 1}};
    tC2 = {{1. / (1. + eps) / (1. + eps), 0, 0},{0, (1. + eps)*(1. + eps), 0}, {0, 0, 1}};
    */

    tC1_0 = {{ 1.0, eps, 0},{eps, 1.+eps*eps, 0}, {0, 0, 1}};
    tC2_0 = {{1.0, -eps, 0},{-eps, 1.+eps*eps, 0}, {0, 0, 1}};
    

    //tC1 = {{ 1.0, eps, 0},{eps, 1., 0}, {0, 0, 1.}};
    //tC2 = {{1.0, -eps, 0},{-eps, 1., 0}, {0, 0, 1.}};
    
    
    
    return IRRT_OK;
}




int
DoubleWellMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    

    if(type == IST_MaxEquivalentStrainLevel) {
      FloatMatrix F, C, C_tC1, C_tC2;
      F.beMatrixForm(status->giveTempFVector());
      C.beTProductOf(F,F);
      C_tC1 = C;
      C_tC2 = C;
      FloatMatrix tC1, tC2;
      tC1 = this->givetC1(tStep);
      tC2 = this->givetC2(tStep);
      C_tC1.subtract(tC1);
      C_tC2.subtract(tC2);
      double normC_tC1, normC_tC2;
      normC_tC1 = C_tC1.computeFrobeniusNorm();
      normC_tC2 = C_tC2.computeFrobeniusNorm();
      answer.resize(1);
      if(normC_tC1 == 0 && normC_tC2 == 0) {
	answer = 0;
      } else {
	answer.at(1) = normC_tC1 * normC_tC1/(normC_tC1 * normC_tC1 + normC_tC2 * normC_tC2);
      }
    } else {
      return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
    return 1;
}




  


} // end namespace oofem
