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

#include "../sm/Materials/Micromorphic/Microstretch/microstretch_elastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
  REGISTER_Material(Microstretch_Elastic);

Microstretch_Elastic :: Microstretch_Elastic(int n, Domain *d) :IsotropicLinearElasticMaterial(n, d), MicromorphicMaterialExtensionInterface(d)
{
  Hk = Ak = 0.;
}

Microstretch_Elastic :: ~Microstretch_Elastic()
{ }



void
Microstretch_Elastic :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}



void
Microstretch_Elastic :: giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  FloatMatrix I;
  I.beSymProjectionMatrix();
  if (matMode == _PlaneStrain) {
    IsotropicLinearElasticMaterial :: givePlaneStrainStiffMtrx(answer, mode, gp, tStep);  
    FloatMatrix IFull;
    IFull = I;
    StructuralMaterial :: giveReducedMatrixForm(I, IFull, matMode);
  } else { //@todo check that all other modes are threated as 3d modes
    IsotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
  }
  
  I.times(Hk);
  answer.add(I);
}

void
Microstretch_Elastic :: giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  FloatMatrix I;
  I.beSymProjectionMatrix();
  
  if(matMode == _PlaneStrain) {
    FloatMatrix IFull;
    IFull = I;
    StructuralMaterial :: giveReducedMatrixForm(I, IFull, matMode);
  }
  
  answer = I;
  answer.times(-Hk);

}


void
Microstretch_Elastic :: giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  this->giveMicromorphicMatrix_dSigdPhi(answer, mode, gp, tStep);
}



void
Microstretch_Elastic :: giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();

  if(matMode == _PlaneStrain) {
    answer.resize(4,4);
  } else {
    answer.resize(6,6);
  } 
  answer.beUnitMatrix();
  answer.times(Hk);

  
}


void
Microstretch_Elastic :: giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  if (matMode == _PlaneStrain) {
    answer.resize(8,8);
  } else {
    answer.resize(18,18);
  } 
  answer.beUnitMatrix();
  answer.times(Ak);

}




void
Microstretch_Elastic :: giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &S, GaussPoint *gp, const FloatArray &strain, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{
  
    MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );

    MaterialMode matMode = gp->giveMaterialMode();    
    FloatMatrix De;
    
    if(matMode == _PlaneStrain) {
      IsotropicLinearElasticMaterial :: givePlaneStrainStiffMtrx(De, TangentStiffness, gp, tStep);   
    } else {
      IsotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix(De, TangentStiffness, gp, tStep);
    } 
    
    sigma.beProductOf(De,strain);     
        
    s = micromorphicVar;
    s.subtract(strain);
    s.times(Hk);

    sigma.subtract(s);

    S = micromorphicVarGrad;
    S.times(Ak);



    status->letTempStrainVectorBe(strain);
    status->letTempMicromorphicVarBe(micromorphicVar);
    status->letTempMicromorphicVarGradBe(micromorphicVarGrad);

    status->letTempStressVectorBe(sigma);
    status->letTempMicromorphicStressBe(s);
    status->letTempMicromorphicStressGradBe(S); 
      
}



IRResultType
Microstretch_Elastic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    IsotropicLinearElasticMaterial :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, Hk, _IFT_MicromorphicMaterialExtensionInterface_Hk);
    IR_GIVE_FIELD(ir, Ak, _IFT_MicromorphicMaterialExtensionInterface_Ak);

    return IRRT_OK;
}

int
Microstretch_Elastic :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );
    
    if ( type == IST_MaxEquivalentStrainLevel ) {
        FloatArray mV = status->giveTempMicromorphicVar();
        answer.resize(1);
	answer.at(1) = mV.at(1);
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }    
}
    



} // end namespace oofem
