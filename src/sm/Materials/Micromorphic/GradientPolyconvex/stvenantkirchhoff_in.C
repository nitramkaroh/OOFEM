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

#include "../sm/Materials/Micromorphic/GradientPolyconvex/stvenantkirchhoff_gradpolyconvex.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
  REGISTER_Material(StVenantKirchhoffGradientPolyconvexMaterial);

StVenantKirchhoffGradientPolyconvexMaterial :: StVenantKirchhoffGradientPolyconvexMaterial(int n, Domain *d) :IsotropicLinearElasticMaterial(n, d)
{
}

StVenantKirchhoffGradientPolyconvexMaterial :: ~StVenantKirchhoffGradientPolyconvexMaterial()
{ }



void
StVenantKirchhoffGradientPolyconvexMaterial :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}


giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{


    // first PK stress
    FloatArray vPm(9);
    // St. Venant-Kirchhoff material model
    FloatArray vE, vS;
    FloatMatrix E;
    E.beTProductOf(F, F);
    E.at(1, 1) -= 1.0;
    E.at(2, 2) -= 1.0;
    E.at(3, 3) -= 1.0;
    E.times(0.5);
    vE.beSymVectorFormOfStrain(E);      // 6


    LinearElasticMaterial::giveRealStressVector_3d(vS, gp, vE, tStep);
    // Compute first PK stress from second PK stress
    FloatMatrix P, S;
    S.beMatrixForm(vS);
    P.beProductOf(F, S); 
    vP.beVectorForm(P);
    
    // micromorphic contribution
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
	for (int m = 1; m <=3; m++) {
	  for (int n = 1; n<=3; n++) {
	    vPm.at(giveVI(i,j)) -= (vInvF.at(giveVI(j,i))*vInvF.at(giveVI(n,m))-vInvF.at(giveVI(n,i))*vInvF.at(giveVI(j,m))) * J * s.at(giveVI(m,n));
	  }
	}
      }
    }
    
    vP.add(vPm);


  
    double J, lnJ, I1, I2, barI1, barI2;
    FloatMatrix F, C, Cpow2, invFt, FC, invF;

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    // compute jacobian and its logarith
    J = F.giveDeterminant();
    lnJ = log(J);

    FloatArray dI1_Cdev_dF, dI2_Cdev_dF, dJ_dF;
    this->compute_dI1_Cdev_dF(dI1_Cdev_dF, F);
    //dI1_Cdev_dF.times(C1);
    this->compute_dI2_Cdev_dF(dI2_Cdev_dF, F);
    //dI2_Cdev_dF.times(C2);
    this->compute_dJ_dF(dJ_dF, F);
    
    answer.zero();
    answer.add(C1, dI1_Cdev_dF);
    answer.add(C2, dI2_Cdev_dF);
    answer.add(K * lnJ / J, dJ_dF);
    

    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);
}



void
MooneyRivlinMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double J, lnJ;
    FloatArray vF;
    FloatMatrix F, invF, invFt, d2I1dF2, d2I2dF2, dinvF_dF, iFtxiFt, dInvF_dF;

    //deformation gradient from the status
    vF = status->giveTempFVector();
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    invF.beInverseOf(F);
    invFt.beTranspositionOf(invF);
    J = F.giveDeterminant();
    lnJ = log(J);
    this->compute_d2I1_Cdev_dF2_and_d2I2_Cdev_dF2(d2I1dF2, d2I2dF2, F);
    FloatMatrix stiff1, stiff2;
    this->computeNumerical_d2I1_Cdev_dF2(stiff1, stiff2, F);


    this->compute_dInvFt_dF(dInvF_dF, invF);
    this->compute_cross_product(iFtxiFt, invFt, invFt);
    
    d2I1dF2.times(C1);
    d2I2dF2.times(C2);

    iFtxiFt.times(K);
    dInvF_dF.times(K * lnJ);

    answer1 = d2I1dF2;
    answer1.add(d2I2dF2);
    answer1.add(iFtxiFt);
    answer1.add(dInvF_dF);

}
  


  
  

IRResultType
StVenantKirchhoffGradientPolyconvexMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    IsotropicLinearElasticMaterial :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, K, _IFT_MooneyRivlinMaterial_k);
    
    
    return IRRT_OK;
}

int
StVenantKirchhoffGradientPolyconvexMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    StVenantKirchhoffGradientPolyconvexMaterialStatus *status = static_cast< StVenantKirchhoffGradientPolyconvexMaterialStatus * >( this->giveStatus(gp) );
    

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
    } else if( type == IST_DeformationGradientTensor) {
      answer = status->giveFVector();
    } else {
      OOFEM_ERROR("Unknown InternalStateType");
    }
    return 1;
}
    

  StVenantKirchhoffGradientPolyconvexMaterialStatus :: StVenantKirchhoffGradientPolyconvexMaterialStatus(int n, Domain *d, GaussPoint *g, bool sym) : MicromorphicMaterialStatus(n, d, g, sym)
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
