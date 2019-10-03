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

#include "../sm/Materials/Micromorphic/GradientPolyconvex/gradientpolyconvexmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"
#include "../sm/Materials/HyperelasticMaterials/doublewellmaterial.h"


namespace oofem {
  REGISTER_Material(GradientPolyconvexMaterial);

GradientPolyconvexMaterial :: GradientPolyconvexMaterial(int n, Domain *d) :StructuralMaterial(n, d), MicromorphicMaterialExtensionInterface(d)
{
  Hk = Ak = 0.;
  gamma = 1.e-5;
  
}

GradientPolyconvexMaterial :: ~GradientPolyconvexMaterial()
{
  delete hyperelasticMaterial;
}


IRResultType
GradientPolyconvexMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    
    IR_GIVE_FIELD(ir, Hk, _IFT_MicromorphicMaterialExtensionInterface_Hk);
    IR_GIVE_FIELD(ir, Ak, _IFT_MicromorphicMaterialExtensionInterface_Ak);
    
    gamma = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, gamma, _IFT_GradientPolyconvexMaterial_gamma );
    
    hyperElasticMaterialType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, hyperElasticMaterialType, _IFT_GradientPolyconvexMaterial_hmt);
    
    if ( hyperElasticMaterialType == 0 ) {
      hyperelasticMaterial = new DoubleWellMaterial(this->giveNumber(), this->giveDomain());
    } else {
      OOFEM_WARNING("Unknown hyperelasticmaterial type %d", hyperElasticMaterialType);
      return IRRT_BAD_FORMAT;
    }


    hyperelasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) {
      return result;
    }
      
    return IRRT_OK;
}



  

void
GradientPolyconvexMaterial :: giveFiniteStrainGeneralizedStressVectors(FloatArray &vP, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &vF, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
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
GradientPolyconvexMaterial :: giveFiniteStrainGeneralizedStressVectors_3d(FloatArray &vP, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &vF, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{

    hyperelasticMaterial->giveFirstPKStressVector_3d(vP, gp, vF, tStep);
    
    GradientPolyconvexMaterialStatus *status = static_cast< GradientPolyconvexMaterialStatus* >( this->giveStatus(gp) );
    
    FloatArray vCofF, vPm;
    // compute cofactor using tensor cross product
    this->compute_2order_tensor_cross_product(vCofF, vF, vF);
    vCofF.times(0.5);  

    // relative stress for cofactor
    s = vCofF;
    s.subtract(micromorphicVar);
    s.times(-Hk);
    // higher order stress
    M = micromorphicVarGrad;
    M.times(Ak);
    
    this->compute_2order_tensor_cross_product(vPm, s, vF);
    vP.subtract(vPm);
    
   
    status->letTempMicromorphicVarBe(micromorphicVar);
    status->letTempMicromorphicVarGradBe(micromorphicVarGrad);

    status->letTempPVectorBe(vP);
    status->letTempFVectorBe(vF);
    status->letTempMicromorphicStressBe(s);
    status->letTempMicromorphicStressGradBe(M); 
      
}

  


void
GradientPolyconvexMaterial :: giveFiniteStrainGeneralizedStressVectors_PlaneStrain(FloatArray &vP, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &reducedvF, const FloatArray &reducedMicromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{

  FloatArray vF, vMV, fullvP, fullS;
  StructuralMaterial :: giveFullVectorFormF(vF, reducedvF, _PlaneStrain);
  StructuralMaterial :: giveFullVectorFormF(vMV, reducedMicromorphicVar, _PlaneStrain);
  this->giveFiniteStrainGeneralizedStressVectors_3d(fullvP, fullS, M, gp, vF, vMV, micromorphicVarGrad, tStep) ;
  StructuralMaterial :: giveReducedVectorForm(vP, fullvP, _PlaneStrain);
  StructuralMaterial :: giveReducedVectorForm(s, fullS, _PlaneStrain);     
}








void
GradientPolyconvexMaterial :: giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    GradientPolyconvexMaterialStatus *status = static_cast< GradientPolyconvexMaterialStatus * >( this->giveStatus(gp) );
    hyperelasticMaterial->give3dMaterialStiffnessMatrix_dPdF(answer, mode, gp, tStep);
    //deformation gradient, its inverse, cofactor, and determinant
    FloatArray vS, vF;
    FloatMatrix xF,xS, xFxF;

    vF = status->giveTempFVector();
    vS = status->giveTempMicromorphicStress();

    this->compute_tensor_cross_product_tensor(xF, vF);
    this->compute_tensor_cross_product_tensor(xS, vS);
    this->compute_4order_tensor_cross_product(xFxF, vF, xF);
    xFxF.times(Hk);

    answer.add(xFxF);
    answer.subtract(xS);


    


    if(gp->giveMaterialMode() == _PlaneStrain) {
      FloatMatrix m3d = answer;
      answer.resize(5,5);
      answer.zero();
      answer.at(1, 1) = m3d.at(1, 1);
      answer.at(1, 2) = m3d.at(1, 2);
      answer.at(1, 3) = m3d.at(1, 3);
      answer.at(1, 4) = m3d.at(1, 6);
      answer.at(1, 5) = m3d.at(1, 9);
      
      answer.at(2, 1) = m3d.at(2, 1);
      answer.at(2, 2) = m3d.at(2, 2);
      answer.at(2, 3) = m3d.at(2, 3);
      answer.at(2, 4) = m3d.at(2, 6);
      answer.at(2, 5) = m3d.at(2, 9);
      
      answer.at(3, 1) = m3d.at(3, 1);
      answer.at(3, 2) = m3d.at(3, 2);
      answer.at(3, 3) = m3d.at(3, 3);
      answer.at(3, 4) = m3d.at(3, 6);
      answer.at(3, 5) = m3d.at(3, 9);
      
      answer.at(4, 1) = m3d.at(6, 1);
      answer.at(4, 2) = m3d.at(6, 2);
      answer.at(4, 3) = m3d.at(6, 3);
      answer.at(4, 4) = m3d.at(6, 6);
      answer.at(4, 5) = m3d.at(6, 9);
      
      answer.at(5, 1) = m3d.at(9, 1);
      answer.at(5, 2) = m3d.at(9, 2);
      answer.at(5, 3) = m3d.at(9, 3);
      answer.at(5, 4) = m3d.at(9, 6);
      answer.at(5, 5) = m3d.at(9, 9);
    }

  
}

void
GradientPolyconvexMaterial :: giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{


    GradientPolyconvexMaterialStatus *status = static_cast< GradientPolyconvexMaterialStatus* >( this->giveStatus(gp) );
    this->compute_tensor_cross_product_tensor(answer, status->giveTempFVector());
    answer.times(-Hk);
    

    if(gp->giveMaterialMode() == _PlaneStrain) {
      FloatMatrix m3d = answer;
      answer.resize(5,5);
      answer.zero();
      answer.at(1, 1) = m3d.at(1, 1);
      answer.at(1, 2) = m3d.at(1, 2);
      answer.at(1, 3) = m3d.at(1, 3);
      answer.at(1, 4) = m3d.at(1, 6);
      answer.at(1, 5) = m3d.at(1, 9);
      
      answer.at(2, 1) = m3d.at(2, 1);
      answer.at(2, 2) = m3d.at(2, 2);
      answer.at(2, 3) = m3d.at(2, 3);
      answer.at(2, 4) = m3d.at(2, 6);
      answer.at(2, 5) = m3d.at(2, 9);
      
      answer.at(3, 1) = m3d.at(3, 1);
      answer.at(3, 2) = m3d.at(3, 2);
      answer.at(3, 3) = m3d.at(3, 3);
      answer.at(3, 4) = m3d.at(3, 6);
      answer.at(3, 5) = m3d.at(3, 9);
      
      answer.at(4, 1) = m3d.at(6, 1);
      answer.at(4, 2) = m3d.at(6, 2);
      answer.at(4, 3) = m3d.at(6, 3);
      answer.at(4, 4) = m3d.at(6, 6);
      answer.at(4, 5) = m3d.at(6, 9);
      
      answer.at(5, 1) = m3d.at(9, 1);
      answer.at(5, 2) = m3d.at(9, 2);
      answer.at(5, 3) = m3d.at(9, 3);
      answer.at(5, 4) = m3d.at(9, 6);
      answer.at(5, 5) = m3d.at(9, 9);
      
    }
     
}


void
GradientPolyconvexMaterial :: giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  this->giveMicromorphicMatrix_dSigdPhi(answer, mode, gp, tStep);   
}



void
GradientPolyconvexMaterial :: giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  MaterialMode matMode = gp->giveMaterialMode();
  if (matMode == _PlaneStrain) {
    answer.resize(5,5);
  } else {
    answer.resize(9,9);
  } 
  answer.beUnitMatrix();
  answer.times(Hk);

}


void
GradientPolyconvexMaterial :: giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
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




int
GradientPolyconvexMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    GradientPolyconvexMaterialStatus *status = static_cast< GradientPolyconvexMaterialStatus * >( this->giveStatus(gp) );
    

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
    } else if( type ==  IST_DeformationGradientTensor ) {
      GradientPolyconvexMaterialStatus *status = static_cast< GradientPolyconvexMaterialStatus* >( this->giveStatus(gp) );
      //deformation gradient, its inverse, cofactor, and determinant
      FloatMatrix F, invF, cofF;
      F.beMatrixForm(status->giveTempFVector());
      invF.beInverseOf(F);
      double J = F.giveDeterminant();
      cofF.beTranspositionOf(invF);
      cofF.times(J);
      answer.beVectorForm(cofF);      
    } else {
      return hyperelasticMaterial->giveIPValue(answer, gp, type, tStep);
    }
    return 1;
}
    

  GradientPolyconvexMaterialStatus :: GradientPolyconvexMaterialStatus(int n, Domain *d, GaussPoint *g, bool sym) : MicromorphicMaterialStatus(n, d, g, sym)
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
