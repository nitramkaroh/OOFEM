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

#include "lsmastermat2.h"
#include "Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "stressvector.h"
#include "strainvector.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LargeStrainMasterMaterial2);

// constructor
LargeStrainMasterMaterial2 :: LargeStrainMasterMaterial2(int n, Domain *d) : StructuralMaterial(n, d)
{
    slaveMat = 0;
}

// destructor
LargeStrainMasterMaterial2 :: ~LargeStrainMasterMaterial2()
{ }

// reads the model parameters from the input file
IRResultType
LargeStrainMasterMaterial2 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    IR_GIVE_OPTIONAL_FIELD(ir, slaveMat, _IFT_LargeStrainMasterMaterial2_slaveMat); // number of slave material
    IR_GIVE_OPTIONAL_FIELD(ir, m, _IFT_LargeStrainMasterMaterial2_m); // type of Set-Hill strain tensor

    return IRRT_OK;
}

// creates a new material status  corresponding to this class
MaterialStatus *
LargeStrainMasterMaterial2 :: CreateStatus(GaussPoint *gp) const
{
  return new LargeStrainMasterMaterial2Status(1, this->giveDomain(), gp, slaveMat);
}

  
void
LargeStrainMasterMaterial2 :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
{
    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );
    this->initTempStatus(gp);

    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );

    FloatArray eVals, vSethHillStrain, vSethHillStress, lambda, eps(3);
    FloatMatrix F, C, eVecs, SethHillStrain, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    F.beMatrixForm(vF);
    //compute right Cauchy-Green tensor(C)
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    // compute Seth - Hill's strain measure, it depends on parameter m
    lambda = eVals;
    if ( m == 0 ) {
      eps.at(1) = 1. / 2. * log(lambda.at(1));
      eps.at(2) = 1. / 2. * log(lambda.at(2));
      eps.at(3) = 1. / 2. * log(lambda.at(3));
    } else {
      eps.at(1) = 1. /  m  * ( pow(lambda.at(1), m / 2. ) - 1. );
      eps.at(2) = 1. /  m  * ( pow(lambda.at(2), m / 2. ) - 1. );
      eps.at(3) = 1. /  m  * ( pow(lambda.at(3), m / 2. ) - 1. );
    }

    SethHillStrain.resize(3, 3);
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
	  SethHillStrain.at(i, j) = eps.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + eps.at(2) *eVecs.at(i, 2) * eVecs.at(j, 2) + eps.at(3) *eVecs.at(i, 3) * eVecs.at(j, 3);  
        }
    }

    vSethHillStrain.beSymVectorFormOfStrain(SethHillStrain);
    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->giveRealStressVector_3d(vSethHillStress, slaveGp, vSethHillStrain, tStep);
    
    SethHillStress.beMatrixForm(vSethHillStress);
    // transformation matrices
    this->giveTransformationMatrices(PP,TL, F, SethHillStress, eVals, eVecs);
    answer.beTProductOf(PP, vSethHillStress);

    status->letTempStressVectorBe(vSethHillStress);
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);

}



void
LargeStrainMasterMaterial2 :: giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &redvE, TimeStep *tStep)
{
    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );

    FloatArray vE, eVals, vSethHillStrain, redvSethHillStrain, vSethHillStress, redvSethHillStress, lambda, eps(3);
    FloatMatrix E, I, C, eVecs, SethHillStrain, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    StructuralMaterial :: giveFullSymVectorForm(vE, redvE, _PlaneStress);
    E.beMatrixFormOfStrain(vE);
    // identity matrix
    I.resize(3,3);
    I.beUnitMatrix();
    //compute right Cauchy-Green tensor(C)
    C = E;
    C.times(2.); 
    C.add(I);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    // compute Seth - Hill's strain measure, it depends on parameter m
    lambda = eVals;
    if ( m == 0 ) {
      eps.at(1) = 1. / 2. * log(lambda.at(1));
      eps.at(2) = 1. / 2. * log(lambda.at(2));
      eps.at(3) = 1. / 2. * log(lambda.at(3));
    } else {
      eps.at(1) = 1. /  m  * ( pow(lambda.at(1), m / 2. ) - 1. );
      eps.at(2) = 1. /  m  * ( pow(lambda.at(2), m / 2. ) - 1. );
      eps.at(3) = 1. /  m  * ( pow(lambda.at(3), m / 2. ) - 1. );
    }

    SethHillStrain.resize(3, 3);
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
	  SethHillStrain.at(i, j) = eps.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + eps.at(2) *eVecs.at(i, 2) * eVecs.at(j, 2) + eps.at(3) *eVecs.at(i, 3) * eVecs.at(j, 3);  
        }
    }

    vSethHillStrain.beSymVectorFormOfStrain(SethHillStrain);
    StructuralMaterial :: giveReducedSymVectorForm(redvSethHillStrain, vSethHillStrain, _PlaneStress);
    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->giveRealStressVector_PlaneStress(redvSethHillStress, slaveGp, redvSethHillStrain, tStep);
    StructuralMaterial :: giveFullSymVectorForm(vSethHillStress, redvSethHillStress, _PlaneStress);
    SethHillStress.beMatrixForm(vSethHillStress);
    // transformation matrices
    this->giveTransformationMatrices_PlaneStress_dSdE(PP,TL, SethHillStress, eVals, eVecs);
    // analytical calculation for m = -2
    FloatMatrix aPP, aTL;
    this->giveTransformationMatrices_PlaneStressAnalyticalM_2_dSdE(PP,TL, C, SethHillStress, eVals, eVecs);
    /*    FloatMatrix redPP;
    redPP.beSubMatrixOf(PP, {1, 2,3,4,5, 6}, {1,2,3,4,5, 6});
    */
    FloatArray fullvSethHillStress(9);
    //StructuralMaterial :: giveFullVectorForm(fullvSethHillStress, redvSethHillStress, _PlaneStress);
    fullvSethHillStress.at(1) = redvSethHillStress.at(1);
    fullvSethHillStress.at(2) = redvSethHillStress.at(2);
    fullvSethHillStress.at(6) = redvSethHillStress.at(3);
    fullvSethHillStress.at(9) = redvSethHillStress.at(3);

    
    FloatArray vS;
    vS.beProductOf(PP, fullvSethHillStress);
    StructuralMaterial :: giveReducedSymVectorForm(answer, vS, _PlaneStress);
    StructuralMaterial :: giveFullSymVectorForm(vS, answer, _PlaneStress);
    
    /*
    FloatMatrix DP, PDP;
    DP.beProductOf(fullStiffness, PP);
    PDP.beProductOf(PP, DP);

    FloatMatrix redTL;
    answer.beSubMatrixOf(PDP, {1, 2, 6}, {1, 2, 6});
    redTL.beSubMatrixOf(TL, {1, 2, 6}, {1, 2, 6});
    answer.add(redTL);     
    
    FloatArray vS;
    vS.beTProductOf(redPP, vSethHillStress);
    StructuralMaterial :: giveReducedSymVectorForm(answer, vS, _PlaneStress);

    */
    
    status->letTempStrainVectorBe(vE);
    status->letTempStressVectorBe(vS);

}
  

void
LargeStrainMasterMaterial2 :: giveFirstPKStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &redvF, TimeStep *tStep)
{
    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );
    this->initTempStatus(gp);

    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );

    FloatArray vF, eVals, vSethHillStrain, redvSethHillStrain, vSethHillStress, redvSethHillStress, lambda, eps(3);
    FloatMatrix F, C, eVecs, SethHillStrain, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    StructuralMaterial :: giveFullVectorFormF(vF, redvF, _PlaneStress);
    F.beMatrixForm(vF);
    //compute right Cauchy-Green tensor(C)
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    // compute Seth - Hill's strain measure, it depends on parameter m
    lambda = eVals;
    if ( m == 0 ) {
      eps.at(1) = 1. / 2. * log(lambda.at(1));
      eps.at(2) = 1. / 2. * log(lambda.at(2));
      eps.at(3) = 1. / 2. * log(lambda.at(3));
    } else {
      eps.at(1) = 1. /  m  * ( pow(lambda.at(1), m / 2. ) - 1. );
      eps.at(2) = 1. /  m  * ( pow(lambda.at(2), m / 2. ) - 1. );
      eps.at(3) = 1. /  m  * ( pow(lambda.at(3), m / 2. ) - 1. );
    }

    SethHillStrain.resize(3, 3);
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
	  SethHillStrain.at(i, j) = eps.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + eps.at(2) *eVecs.at(i, 2) * eVecs.at(j, 2) + eps.at(3) *eVecs.at(i, 3) * eVecs.at(j, 3);  
        }
    }

    vSethHillStrain.beSymVectorFormOfStrain(SethHillStrain);
    StructuralMaterial :: giveReducedSymVectorForm(redvSethHillStrain, vSethHillStrain, _PlaneStress);
    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->giveRealStressVector_PlaneStress(redvSethHillStress, slaveGp, redvSethHillStrain, tStep);
    StructuralMaterial :: giveFullSymVectorForm(vSethHillStress, redvSethHillStress, _PlaneStress);
    SethHillStress.beMatrixForm(vSethHillStress);
    // transformation matrices
    this->giveTransformationMatrices_PlaneStress(PP,TL, F, SethHillStress, eVals, eVecs);
    FloatMatrix PPtest, TLtest;
    //    this->giveTransformationMatrices_PlaneStressTest(PPtest, TLtest,  F, SethHillStress, eVals, eVecs);
    FloatMatrix PP3d, TL3d;
    this->giveTransformationMatrices(PP3d,TL3d, F, SethHillStress, eVals, eVecs);
    
    
    
    FloatArray vP;
    vP.beTProductOf(PP, vSethHillStress);
    StructuralMaterial :: giveReducedVectorForm(answer, vP, _PlaneStress);



    
    

    status->letTempStressVectorBe(redvSethHillStress);
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(vP);

}



void LargeStrainMasterMaterial2 :: giveFirstPKStressVector_Membrane2d(FloatArray &answer, GaussPoint *gp, const FloatArray &redvF, TimeStep *tStep)
{
    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );
    this->initTempStatus(gp);

    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );

    FloatArray vF, eVals, vSethHillStrain, redvSethHillStrain, vSethHillStress, redvSethHillStress, lambda, eps(3);
    FloatMatrix F, C, eVecs, SethHillStrain, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    StructuralMaterial :: giveFullVectorFormF(vF, redvF, _Membrane2d);
    F.beMatrixForm(vF);
    //compute right Cauchy-Green tensor(C)
    FloatMatrix Cred;
    C.beTProductOf(F, F);
    Cred.beSubMatrixOf(C,{1,2},{1,2});
    // compute eigen values and eigen vectors of C
    Cred.jaco_(eVals, eVecs, 25);
    eVals.resize(3);
    eVals.at(3) = 1.;
    eVecs.resizeWithData(3,3);
    eVecs.at(3,3) = 1.;
    // compute Seth - Hill's strain measure, it depends on parameter m
    lambda = eVals;
    if ( m == 0 ) {
      eps.at(1) = 1. / 2. * log(lambda.at(1));
      eps.at(2) = 1. / 2. * log(lambda.at(2));
      eps.at(3) = 1. / 2. * log(lambda.at(3));
    } else {
      eps.at(1) = 1. /  m  * ( pow(lambda.at(1), m / 2. ) - 1. );
      eps.at(2) = 1. /  m  * ( pow(lambda.at(2), m / 2. ) - 1. );
      eps.at(3) = 1. /  m  * ( pow(lambda.at(3), m / 2. ) - 1. );
    }

    SethHillStrain.resize(3, 3);
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
	  SethHillStrain.at(i, j) = eps.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + eps.at(2) *eVecs.at(i, 2) * eVecs.at(j, 2) + eps.at(3) *eVecs.at(i, 3) * eVecs.at(j, 3);  
        }
    }

    vSethHillStrain.beSymVectorFormOfStrain(SethHillStrain);
    StructuralMaterial :: giveReducedSymVectorForm(redvSethHillStrain, vSethHillStrain, _PlaneStress);
    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->giveRealStressVector_PlaneStress(redvSethHillStress, slaveGp, redvSethHillStrain, tStep);
    StructuralMaterial :: giveFullSymVectorForm(vSethHillStress, redvSethHillStress, _Membrane2d);
    SethHillStress.beMatrixForm(vSethHillStress);
    // get full strain including out-of-plane component
    StructuralMaterialStatus *slaveStatus = static_cast< StructuralMaterialStatus * > (sMat->giveStatus(slaveGp));
    FloatArray fullSethHillStrain;
    fullSethHillStrain  = slaveStatus->giveStrainVector();
    eVals.at(3) = fullSethHillStrain.at(3) + 1.;
    // transformation matrices
    this->giveTransformationMatrices_PlaneStress(PP,TL, F, SethHillStress, eVals, eVecs);
    FloatMatrix PPtest, TLtest, PPh, TLh;
    //    this->giveTransformationMatrices_PlaneStressTest(PPtest, TLtest,  F, SethHillStress, eVals, eVecs);
    FloatMatrix PP3d, TL3d;
    this->giveTransformationMatrices(PP3d,TL3d, F, SethHillStress, eVals, eVecs);
       
    
    FloatArray vP;
    FloatMatrix invF, S, P;
    invF.beInverseOf(F);
    //vP.beTProductOf(PP, vSethHillStress);
    vP.beTProductOf(PP3d, vSethHillStress);
    P.beMatrixForm(vP);
    S.beProductOf(invF, P);
    //vP.beTProductOf(PPtest, vSethHillStress);
    StructuralMaterial :: giveReducedVectorForm(answer, vP, _Membrane2d);    
    

    status->letTempStressVectorBe(redvSethHillStress);
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(vP);

}

  
  
 


void
LargeStrainMasterMaterial2 :: giveCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
{
    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );
    this->initTempStatus(gp);

    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );

    FloatArray eVals, vSethHillStrain, vSethHillStress, lambda, eps(3);
    FloatMatrix F, b, eVecs, SethHillStrain, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    F.beMatrixForm(vF);
    //compute right Cauchy-Green tensor(C)
    b.beProductTOf(F, F);
    // compute eigen values and eigen vectors of C
    b.jaco_(eVals, eVecs, 15);
    // compute Seth - Hill's strain measure, it depends on parameter m
    lambda = eVals;
    if ( m == 0 ) {
      eps.at(1) = 1. / 2. * log(lambda.at(1));
      eps.at(2) = 1. / 2. * log(lambda.at(2));
      eps.at(3) = 1. / 2. * log(lambda.at(3));
    } else {
      eps.at(1) = 1. /  m  * ( pow(lambda.at(1), m / 2. ) - 1. );
      eps.at(2) = 1. /  m  * ( pow(lambda.at(2), m / 2. ) - 1. );
      eps.at(3) = 1. /  m  * ( pow(lambda.at(3), m / 2. ) - 1. );
    }

    SethHillStrain.resize(3, 3);
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
	  SethHillStrain.at(i, j) = eps.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + eps.at(2) *eVecs.at(i, 2) * eVecs.at(j, 2) + eps.at(3) *eVecs.at(i, 3) * eVecs.at(j, 3);
        }
    }

    vSethHillStrain.beSymVectorFormOfStrain(SethHillStrain);
    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->giveRealStressVector_3d(vSethHillStress, slaveGp, vSethHillStrain, tStep);
    
    SethHillStress.beMatrixForm(vSethHillStress);
    // transformation matrices
    this->giveTransformationMatrices(PP,TL, F, SethHillStress, eVals, eVecs);
    answer.beTProductOf(PP, vSethHillStress);
    // store F, P, and Seth-Hill stress into status
    status->letTempStressVectorBe(vSethHillStress);
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);
}



void
LargeStrainMasterMaterial2 :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    FloatArray eVals;
    FloatMatrix stiffness, F, C, P, eVecs, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    F.beMatrixForm(status->giveTempFVector());
    //compute right Cauchy-Green tensor(C), its eigenvalues and eigenvectors
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    

    FloatArray vSethHillStress, vP(status->giveTempPVector());

    if(gp->giveMaterialMode() != _3dMat) {
      StructuralMaterial :: giveFullSymVectorForm( vSethHillStress, status->giveTempStressVector(), gp->giveMaterialMode() );
    } else {
      vSethHillStress =  status->giveTempStressVector();
    }


    P.beMatrixForm(vP);    
    SethHillStress.beMatrixForm(vSethHillStress);
    FloatMatrix invF, S, EP, delta_S;;
    // compute 2-PK stress
    invF.beInverseOf(F);
    S.beProductOf(invF, P);

    this->giveTransformationMatrices(PP,TL, F, SethHillStress, eVals, eVecs);
    this->giveDeltaS_Product(delta_S, S);


    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->give3dMaterialStiffnessMatrix(stiffness, mode, slaveGp, tStep);
    EP.beProductOf(stiffness, PP);

    answer.beTProductOf(PP, EP);
    answer.add(TL);   
    answer.add(delta_S);  

}


void
LargeStrainMasterMaterial2 :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    FloatArray eVals, vE;
    FloatMatrix stiffness, E, C, I, P, eVecs, SethHillStress, PP, TL;
    //store of Green-Lagrangian sstrain into into 3x3 matrix
    StructuralMaterial :: giveFullSymVectorForm(vE, status->giveTempStrainVector(), _PlaneStress);
    E.beMatrixFormOfStrain(vE);
    // identity matrix
    I.resize(3,3);
    I.beUnitMatrix();
    //compute right Cauchy-Green tensor(C)
    C = E;
    C.times(2.); 
    C.add(I);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    // slave status and slave gp
    GaussPoint *slaveGp = status->giveSlaveGp();
    StructuralMaterialStatus *slaveStatus = static_cast< StructuralMaterialStatus * > (sMat->giveStatus(slaveGp));    

    FloatArray vSethHillStress;
    StructuralMaterial :: giveFullSymVectorForm(vSethHillStress, slaveStatus->giveTempStressVector(), _PlaneStress);
    SethHillStress.beMatrixForm(vSethHillStress);
    // get full strain including out-of-plane component
    FloatArray fullSethHillStrain;
    fullSethHillStrain  = slaveStatus->giveStrainVector();
    eVals.at(3) = fullSethHillStrain.at(3) + 1.;
    
    FloatMatrix EP, PPt, TLt, PP2, TL2;
    this->giveTransformationMatrices_PlaneStress_dSdE(PP,TL, SethHillStress, eVals, eVecs);
    this->giveTransformationMatrices_PlaneStress2_dSdE(PP2,TL2, SethHillStress, eVals, eVecs);
    this->giveTransformationMatrices(PPt,TLt, I, SethHillStress, eVals, eVecs);
    // analytical calculation for m = -2
    FloatMatrix aPP, aTL;
    this->giveTransformationMatrices_PlaneStressAnalyticalM_2_dSdE(PP,TL, C, SethHillStress, eVals, eVecs);


    
    sMat->givePlaneStressStiffMtrx(stiffness, mode, slaveGp, tStep);

    FloatMatrix fStiff;
    StructuralMaterial :: giveFullSymMatrixForm(fStiff, stiffness, _PlaneStress);

    
    FloatMatrix fullStiffness(9,9);
    for (int k = 1; k <= 2; k++) {
      for (int l = 1; l <= 2; l++) {
	for (int m = 1; m <=2; m++) {
	  for (int n = 1; n<=2; n++) {
	    fullStiffness.at(giveVI(k,l),giveVI(m,n)) = fStiff.at(giveSymVI(k,l),giveSymVI(m,n));
	  }
	}
      }
    }

    FloatMatrix DP, PDP, dSdE;
    DP.beProductOf(fullStiffness, PP);
    PDP.beProductOf(PP, DP);
    dSdE = PDP;
    dSdE.add(TL);

    FloatMatrix tF, F;
    tF = {{1.008196805297511e+00,-9.515121653015688e-03,0.000000000000000e+00},{-7.972061560193748e-03,9.961966158206870e-01,0.000000000000000e+00},{6.139202182449857e-04,6.139202182449857e-04,1.000000000000000e+00}};
    tF = {{1.000000000000000e+00,0.000000000000000e+00,-5.355047384154995e-01},{0.000000000000000e+00,1.000000000000000e+00,-8.032169467759106e-01},{2.052402044974549e+00,3.078449145004574e+00,2.609161005889273e-01}};
    
    F.beTranspositionOf(tF);
    FloatArray vS, vF, vO(6);
    StructuralMaterial :: giveFullSymVectorForm(vS, status->giveTempStressVector(), _PlaneStress);
    vF.beVectorForm(F);
    FloatMatrix redTL;
    answer.beSubMatrixOf(PDP, {1, 2, 6}, {1, 2, 6});
    redTL.beSubMatrixOf(TLt, {1, 2, 6}, {1, 2, 6});
    answer.add(redTL);

    
    FloatMatrix dPdF, TLf1, TLf2;
    StructuralMaterial :: convert_dSdE_2_dPdF(dPdF, dSdE, vS, vF, _3dMat);
    StructuralMaterial :: convert_dSdE_2_dPdF(TLf1, TLt, vO, vF, _3dMat);

    
    
    /*
    FloatMatrix redPP, redTL;
    redPP.beSubMatrixOf(PP, {1, 2, 6}, {1, 2, 6});
    redTL.beSubMatrixOf(TL, {1, 2, 6}, {1, 2, 6});
    
    EP.beProductOf(stiffness, redPP);
    answer.beTProductOf(redPP, EP);
    answer.add(redTL);   
    */
}



void
LargeStrainMasterMaterial2 :: givePlaneStressStiffMtrx_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    FloatArray eVals;
    FloatMatrix stiffness, F, C, P, eVecs, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    F.beMatrixForm(status->giveTempFVector());
    //compute right Cauchy-Green tensor(C), its eigenvalues and eigenvectors
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    

    FloatArray vSethHillStress, vP(status->giveTempPVector());

    vSethHillStress =  status->giveTempStressVector();
    P.beMatrixForm(vP);    
    SethHillStress.beMatrixForm(vSethHillStress);
    FloatMatrix invF, S, EP, delta_S;;
    // compute 2-PK stress
    invF.beInverseOf(F);
    S.beProductOf(invF, P);

    this->giveTransformationMatrices_PlaneStress(PP,TL, F, SethHillStress, eVals, eVecs);
    this->giveDeltaS_Product_PlaneStress(delta_S, S);

   

    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->givePlaneStressStiffMtrx(stiffness, mode, slaveGp, tStep);
    EP.beProductOf(stiffness, PP);

    answer.beTProductOf(PP, EP);
    answer.add(TL);   
    answer.add(delta_S);  

}

  


void
LargeStrainMasterMaterial2 :: giveMembrane2dStiffMtrx_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep)
{

    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    FloatArray eVals;
    FloatMatrix stiffness, F, C, P, eVecs, SethHillStress, PP, TL;
    //store of deformation gradient into 3x3 matrix
    F.beMatrixForm(status->giveTempFVector());
    //compute right Cauchy-Green tensor(C)
    FloatMatrix Cred;
    C.beTProductOf(F, F);
    Cred.beSubMatrixOf(C,{1,2},{1,2});
    // compute eigen values and eigen vectors of C
    Cred.jaco_(eVals, eVecs, 25);
    eVals.resize(3);
    eVals.at(3) = 1.;
    eVecs.resizeWithData(3,3);
    eVecs.at(3,3) = 1.;
    
    FloatArray fullSethHillStrain;
    GaussPoint *slaveGp = status->giveSlaveGp();
    StructuralMaterialStatus *slaveStatus = static_cast< StructuralMaterialStatus * > (sMat->giveStatus(slaveGp));    
    fullSethHillStrain  = slaveStatus->giveStrainVector();
    eVals.at(3) = fullSethHillStrain.at(3) + 1.;

    double l3 = sqrt(eVals.at(3));
    FloatArray n, N;
    FloatMatrix invF;
    invF.beInverseOf(F);
    N = {0,0,1};
    n.beTProductOf(invF,N);
    n.normalize();
    F.at(1,3) = l3*n.at(1);
    F.at(2,3) = l3*n.at(2);
    F.at(3,3) = l3*n.at(3);
    
    FloatArray vSethHillStress, vP(status->giveTempPVector());    
    StructuralMaterial :: giveFullSymVectorForm( vSethHillStress, status->giveTempStressVector(), _PlaneStress );    
    P.beMatrixForm(vP);
    /////////////////
    //vSethHillStress = {1.407420628432634e+04,1.407420264388231e+04,0.000000000000000e+00,0.000000000000000e+00,0.000000000000000e+00,3.152716993387060e-03};
    //////////////////////////    
    SethHillStress.beMatrixForm(vSethHillStress);    
    FloatMatrix S, EP, delta_S;;
    // compute 2-PK stress
    invF.beInverseOf(F);
    S.beProductOf(invF, P);

    
    FloatMatrix PPt, TLt, PPps, TLps, Pc;

    this->giveTransformationMatrices_PlaneStress(PPt,TLt, F, SethHillStress, eVals, eVecs);
    Pc.beTProductOf(PPt, vSethHillStress);
    FloatMatrix PPtest, TLtest, PPh, TLh;
    //    this->giveTransformationMatrices_PlaneStressTest(PPtest, TLtest,  F, SethHillStress, eVals, eVecs);
    this->giveTransformationMatrices_huhu(PPh,TLh, F, SethHillStress, eVals, eVecs);
    this->giveTransformationMatrices(PP,TL, F, SethHillStress, eVals, eVecs);
    this->giveTransformationMatrices_PlaneStress(PPps,TLps,F, SethHillStress, eVals, eVecs);
    if(S.at(1,1) == 0 && S.at(2,2) == 0) {
      S.at(1,1) = 10;
      S.at(2,2) = 10;
    }
    ///////////////
    /*    vP.beTProductOf(PP, vSethHillStress);
    P.beMatrixForm(vP);	
    S.beProductOf(invF, P);
    */
    /////////////

    
    this->giveDeltaS_Product(delta_S, S);

    FloatMatrix h(TLps);
    h.add(delta_S);


    //    this->giveTransformationMatrices_PlaneStress(PP,TL, F, SethHillStress, eVals, eVecs);    
    FloatMatrix redPP(3,6), redTL(6,6), redDelta_S(6,6), redH(6,6), redTLtest(6,6);
    redPP.beSubMatrixOf(PP, {1, 2, 6},{1, 2, 6, 7, 8, 9});
    redTL.beSubMatrixOf(TLps, {1, 2, 6, 7, 8, 9}, {1, 2, 6, 7, 8, 9});
    redDelta_S.beSubMatrixOf(delta_S, {1, 2, 6, 7, 8, 9}, {1, 2, 6, 7, 8, 9});
    redH.beSubMatrixOf(h, {1, 2, 6, 7, 8, 9}, {1, 2, 6, 7, 8, 9});
    redTLtest.beSubMatrixOf(TLtest, {1, 2, 6, 7, 8, 9}, {1, 2, 6, 7, 8, 9});
    //redTL = redTLtest;
    sMat->givePlaneStressStiffMtrx(stiffness, mode, slaveGp, tStep);
    EP.beProductOf(stiffness, redPP);
   
    answer.beTProductOf(redPP, EP);
    answer.add(redTL);   
    answer.add(redDelta_S);



 /////////////////////////////////////////////
    /* double perturbation = 1.e-6;  
    FloatArray vF, vFp, Pp, vPp;
    vF = status->giveTempFVector();
    vFp.beSubArrayOf(vF, {1, 2, 6, 7, 8, 9});
    vPp.beSubArrayOf(vP, {1, 2, 6, 7, 8, 9});
    FloatMatrix stiff(6,6);
    stiff.zero();
    vFp.at(1) += (perturbation);
    this->giveFirstPKStressVector_Membrane2d(Pp, gp, vFp, tStep);   
    for(int i =1; i <=6; i++) {
      stiff.at(i,1) = Pp.at(i) - vPp.at(i);    
    }

    vFp.beSubArrayOf(vF, {1, 2, 6, 7, 8, 9});
    vFp.at(2) += (perturbation);
    this->giveFirstPKStressVector_Membrane2d(Pp, gp, vFp, tStep);  
    for(int i =1; i <=6; i++) {
      stiff.at(i,2) = Pp.at(i) - vPp.at(i);    
    }


    vFp.beSubArrayOf(vF, {1, 2, 6, 7, 8, 9});
    vFp.at(3) += (perturbation);
    this->giveFirstPKStressVector_Membrane2d(Pp, gp, vFp, tStep);   
    for(int i =1; i <=6; i++) {
      stiff.at(i,3) = Pp.at(i) - vPp.at(i);    
    }


    vFp.beSubArrayOf(vF, {1, 2, 6, 7, 8, 9});
    vFp.at(4) += (perturbation);
    this->giveFirstPKStressVector_Membrane2d(Pp, gp, vFp, tStep);  
    for(int i =1; i <=6; i++) {
      stiff.at(i,4) = Pp.at(i) - vPp.at(i);    
    }


    vFp.beSubArrayOf(vF, {1, 2, 6, 7, 8, 9});
    vFp.at(5) += (perturbation);
    this->giveFirstPKStressVector_Membrane2d(Pp, gp, vFp, tStep);   
    for(int i =1; i <=6; i++) {
      stiff.at(i,5) = Pp.at(i) - vPp.at(i);    
    }


    vFp.beSubArrayOf(vF, {1, 2, 6, 7, 8, 9});
    vFp.at(6) += (perturbation);
    this->giveFirstPKStressVector_Membrane2d(Pp, gp, vFp, tStep);   
    for(int i =1; i <=6; i++) {
      stiff.at(i,6) = Pp.at(i) - vPp.at(i);    
    }


    stiff.times(1./perturbation);*/
    /////////////////////////////////////////////

    

}



  
void
LargeStrainMasterMaterial2 :: giveSpatial3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    FloatArray eVals;
    FloatMatrix stiffness, F, b, C, eVecs, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    F.beMatrixForm(status->giveTempFVector());
    //compute right Cauchy-Green tensor(C), its eigenvalues and eigenvectors
    b.beProductTOf(F, F);
    // compute eigen values and eigen vectors of C
    b.jaco_(eVals, eVecs, 15);
    

    FloatArray vSethHillStress, vC(status->giveTempCVector());

    if(gp->giveMaterialMode() != _3dMat) {
      StructuralMaterial :: giveFullSymVectorForm( vSethHillStress, status->giveTempStressVector(), gp->giveMaterialMode() );
    } else {
      vSethHillStress =  status->giveTempStressVector();
    }

    C.beMatrixForm(vC);    
    SethHillStress.beMatrixForm(vSethHillStress);
    FloatMatrix invF, S, EP, delta_S;;

    this->giveTransformationMatrices(PP,TL, F, S, eVals, eVecs);
    this->giveDeltaS_Product(delta_S, S);


    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->give3dMaterialStiffnessMatrix(stiffness, mode, slaveGp, tStep);
    EP.beProductOf(stiffness, PP);

    answer.beTProductOf(PP, EP);
    answer.add(TL);   
    answer.add(delta_S);
    
}





void
LargeStrainMasterMaterial2 :: giveDeltaS_Product(FloatMatrix &answer, const FloatMatrix &S)
{   

    answer.resize(9,9);
  
    answer.at(1,1) = S.at(1,1);
    answer.at(1,5) = S.at(1,3);
    answer.at(1,6) = S.at(1,2);
    
    answer.at(2,2) = S.at(2,2);
    answer.at(2,4) = S.at(2,3);
    answer.at(2,9) = S.at(1,2);
    
    answer.at(3,3) = S.at(3,3);
    answer.at(3,7) = S.at(2,3);
    answer.at(3,8) = S.at(2,3);
    
    
    answer.at(4, 2) = S.at(3,2);
    answer.at(4, 4) = S.at(3,3);
    answer.at(4, 9) = S.at(3,1);
    
    answer.at(5, 1) = S.at(3,1);
    answer.at(5, 5) = S.at(3,3);
    answer.at(5, 6) = S.at(3,2);
    
    answer.at(6, 1) = S.at(2,1);
    answer.at(6, 5) = S.at(2,3);
    answer.at(6, 6) = S.at(2,2);
    
    answer.at(7, 3) = S.at(2,3);
    answer.at(7, 7) = S.at(2,2);
    answer.at(7, 8) = S.at(2,1);
    
    answer.at(8, 3) = S.at(1,3);
    answer.at(8, 7) = S.at(1,2);
    answer.at(8, 8) = S.at(1,1);
    
    answer.at(9, 2) = S.at(1,2);
    answer.at(9, 4) = S.at(1,3);
    answer.at(9, 9) = S.at(1,1);
}

void
LargeStrainMasterMaterial2 :: giveDeltaS_Product_PlaneStress(FloatMatrix &answer, const FloatMatrix &S)
{   

    answer.resize(9,9);
  
    answer.at(1,1) = S.at(1,1);
    answer.at(1,5) = S.at(1,3);
    answer.at(1,6) = S.at(1,2);
    
    answer.at(2,2) = S.at(2,2);
    answer.at(2,4) = S.at(2,3);
    answer.at(2,9) = S.at(1,2);
    
    answer.at(3,3) = S.at(3,3);
    answer.at(3,7) = S.at(2,3);
    answer.at(3,8) = S.at(2,3);
    
    
    answer.at(4, 2) = S.at(3,2);
    answer.at(4, 4) = S.at(3,3);
    answer.at(4, 9) = S.at(3,1);
    
    answer.at(5, 1) = S.at(3,1);
    answer.at(5, 5) = S.at(3,3);
    answer.at(5, 6) = S.at(3,2);
    
    answer.at(6, 1) = S.at(2,1);
    answer.at(6, 5) = S.at(2,3);
    answer.at(6, 6) = S.at(2,2);
    
    answer.at(7, 3) = S.at(2,3);
    answer.at(7, 7) = S.at(2,2);
    answer.at(7, 8) = S.at(2,1);
    
    answer.at(8, 3) = S.at(1,3);
    answer.at(8, 7) = S.at(1,2);
    answer.at(8, 8) = S.at(1,1);
    
    answer.at(9, 2) = S.at(1,2);
    answer.at(9, 4) = S.at(1,3);
    answer.at(9, 9) = S.at(1,1);
}

void
LargeStrainMasterMaterial2 :: giveTransformationMatrices(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &F, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N)
{   

    FloatMatrix n;
    n.beProductOf(F,N);
    
    FloatArray  d(3),f(3), eps(3);
    if(this->m == 0) { // log strain formulation (m = 0) has a special treatment
      d.at(1)= 1./lam.at(1);
      d.at(2)= 1./lam.at(2);
      d.at(3)= 1./lam.at(3);
      
      f.at(1) = -2./lam.at(1)/lam.at(1);
      f.at(2) = -2./lam.at(2)/lam.at(2);
      f.at(3) = -2./lam.at(3)/lam.at(3);

      eps.at(1) = log(lam.at(1))/2.;
      eps.at(2) = log(lam.at(2))/2.;
      eps.at(3) = log(lam.at(3))/2.;
    } else {  // all other strain measures (m != 0)
      d.at(1)= pow(lam.at(1), m/2.-1);
      d.at(2)= pow(lam.at(2), m/2.-1);
      d.at(3)= pow(lam.at(3), m/2.-1);
      
      f.at(1) = 2. * ( m/2. - 1 ) * pow( lam.at(1), m/2. - 2. );
      f.at(2) = 2. * ( m/2. - 1 ) * pow( lam.at(2), m/2. - 2. );
      f.at(3) = 2. * ( m/2. - 1 ) * pow( lam.at(3), m/2. - 2. );
      
      eps.at(1) = 1./m * ( pow(lam.at(1),m/2.) - 1. );
      eps.at(2) = 1./m * ( pow(lam.at(2),m/2.) - 1. );
      eps.at(3) = 1./m * ( pow(lam.at(3),m/2.) - 1. );

    }

    double eta = 0.;
    FloatMatrix TN, dzeta(3,3), ksi(3,3), theta(3,3);

    TN.beProductOf(SethHillStress,N);
    dzeta.beTProductOf(N, TN);
    




    // compute auxiliary variables 
    // the computation differes depends on if the eigenvalues of C are equal or not
    if(lam.at(1) != lam.at(2)) {
      if(lam.at(2) != lam.at(3)) {
	if(lam.at(1) != lam.at(3)) {
	  // all eigenvalues are different
	  for(int i = 1; i <= 3; i++) {
	    for(int j = 1; j <= 3; j++) {
	      if(i == j) {
		continue;
	      } else {
		theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
		ksi.at(i,j) = (theta.at(i,j) - 1./2.*d.at(j))/(lam.at(i) - lam.at(j));
		for(int k = 1; k <= 3; k++) {
		  if((k != i) && (k != j)) {
		    eta += eps.at(i)/2./(lam.at(i) - lam.at(j))/(lam.at(i)-lam.at(k));
		  }
		}
	      }
	    }
	  }
	} else { //l1 == l3 && l1 != l2
	  for(int i = 1; i <= 3; i++) {
	    for(int j = 1; j <= 3; j++) {
	      if( i == j ) {
		continue;
	      } else {
		if((i == 1 && j == 3) || (i == 3 && j == 1) ) {
		  theta.at(i,j) = 1./2.*d.at(i);
		  ksi.at(i,j) = 1./8.*f.at(i);
		} else {
		  theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
		  ksi.at(i,j) = (theta.at(i,j) - 1./2.*d.at(j))/(lam.at(i) - lam.at(j));
		}
	      }
	    }
	  }	  
	  eta = ksi.at(1,2);
	}
      } else { //l2 == l3 && l1 != l2
	for(int i = 1; i <= 3; i++) {
	  for(int j = 1; j <= 3; j++) {
	    if( i == j ) {
	      continue;
	    } else {
	      if((i == 2 && j == 3) || (i == 3 && j == 2) ) {
		theta.at(i,j) = 1./2.*d.at(i);
		ksi.at(i,j) = 1./8.*f.at(i);
	      } else {
		theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
		ksi.at(i,j) = (theta.at(i,j) - 1./2.*d.at(j))/(lam.at(i) - lam.at(j));
	      }
	    }
	  }
	} 
	eta = ksi.at(1,2);    
      }
    } else if(lam.at(1) != lam.at(3)) { // l1 == l2  && l1 != l3
      for(int i = 1; i <= 3; i++) {
	for(int j = 1; j <= 3; j++) {
	  if( i == j ) {
	    continue;
	  } else {
	    if((i == 1 && j == 2) || (i == 2 && j == 1) ) {
	      theta.at(i,j) = 1./2.*d.at(i);
	      ksi.at(i,j) = 1./8.*f.at(i);
	    } else {
	      theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
	      ksi.at(i,j) = (theta.at(i,j) - 1./.2*d.at(j))/(lam.at(i) - lam.at(j));
	    }
	  }
	}
      } 
      eta = ksi.at(1,3);
    } else {  // l1 == l2 == l3
      for(int i = 1; i <= 3; i++) {
	for(int j = 1; j <= 3; j++) {
	  theta.at(i,j) = 1./2. * d.at(i);
	  ksi.at(i,j) = 1./8. * f.at(i);
	}
      }
      eta = 1./8. * f.at(1);
    }


  FloatMatrix M(9,9);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
	for (int l = 1; l <=3; l++) {
	  M.at(giveVI(i,j),giveVI(k,l)) = n.at(k,i)*N.at(l,j) + n.at(k,j)*N.at(l,i);
	}
      }
    }
  }
    
  PP.resize(6,9), TL.resize(9,9);
  for (int k = 1; k <= 3; k++) {
    for (int l = 1; l <= 3; l++) {
      for (int m = 1; m <=3; m++) {
	for (int n = 1; n<=3; n++) {
	  for (int i = 1; i <= 3; i++) {
	    if(k == l ) {
	      PP.at(giveSymVI(k,l),giveVI(m,n)) += 0.5 * d.at(i) * N.at(k,i) * N.at(l,i) * M.at(giveVI(i,i),giveVI(m,n));
	    } else if(k < l) {
	      PP.at(giveSymVI(k,l),giveVI(m,n)) += d.at(i) * N.at(k,i) * N.at(l,i) * M.at(giveVI(i,i),giveVI(m,n));
	    }

	    TL.at(giveVI(k,l),giveVI(m,n)) += 1./4.*f.at(i)*dzeta.at(i,i)*M.at(giveVI(i,i),giveVI(k,l))*M.at(giveVI(i,i),giveVI(m,n));
	    for (int j  = 1; j <= 3; j++) {
	      if(j != i) {
		if( k == l) {
		  PP.at(giveSymVI(k,l),giveVI(m,n)) += theta.at(i,j)*N.at(k,i)*N.at(l,j)*M.at(giveVI(i,j),giveVI(m,n));
		} else if(k < l) {
		  PP.at(giveSymVI(k,l),giveVI(m,n)) += 2. * theta.at(i,j)*N.at(k,i)*N.at(l,j)*M.at(giveVI(i,j),giveVI(m,n));
		}
		TL.at(giveVI(k,l),giveVI(m,n)) += 2.*ksi.at(i,j) * ( dzeta.at(i,j)*(M.at(giveVI(i,j),giveVI(k,l))*M.at(giveVI(j,j),giveVI(m,n)) + M.at(giveVI(j,j),giveVI(k,l))*M.at(giveVI(i,j),giveVI(m,n))) + dzeta.at(j,j)*M.at(giveVI(i,j),giveVI(k,l))*M.at(giveVI(i,j),giveVI(m,n) ));
		for (int q = 1; q <= 3; q++) {
		  if(q == i) {
		    continue;
		  } else if(q == j) {
		    continue;
		  } else {
		    TL.at(giveVI(k,l),giveVI(m,n)) += 2.*eta*dzeta.at(i,j)*M.at(giveVI(i,q),giveVI(k,l))*M.at(giveVI(j,q),giveVI(m,n));
		  }
		}
	      }
	    }
	  }
	}	    
      }
    }
  }   



}



void
LargeStrainMasterMaterial2 :: giveTransformationMatrices_PlaneStress(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &F, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N)
{  

  FloatMatrix n;
  n.beProductOf(F,N);
    
  FloatArray  d(3),f(3), eps(3);
  if(this->m == 0) { // log strain formulation (m = 0) has a special treatment
    d.at(1)= 1./lam.at(1);
    d.at(2)= 1./lam.at(2);
    d.at(3)= 1./lam.at(3);
    
    f.at(1) = -2./lam.at(1)/lam.at(1);
    f.at(2) = -2./lam.at(2)/lam.at(2);
    f.at(3) = -2./lam.at(3)/lam.at(3);
    
    eps.at(1) = log(lam.at(1))/2.;
    eps.at(2) = log(lam.at(2))/2.;
    eps.at(3) = log(lam.at(3))/2.;
    
  } else {  // all other strain measures (m != 0)
    d.at(1)= pow(lam.at(1), m/2.-1);
    d.at(2)= pow(lam.at(2), m/2.-1);
    d.at(3)= pow(lam.at(3), m/2.-1);
      
    f.at(1) = 2. * ( m/2. - 1 ) * pow( lam.at(1), m/2. - 2. );
    f.at(2) = 2. * ( m/2. - 1 ) * pow( lam.at(2), m/2. - 2. );
    f.at(3) = 2. * ( m/2. - 1 ) * pow( lam.at(3), m/2. - 2. );
      
    eps.at(1) = 1./m * ( pow(lam.at(1),m/2.) - 1. );
    eps.at(2) = 1./m * ( pow(lam.at(2),m/2.) - 1. );
    eps.at(3) = 1./m * ( pow(lam.at(3),m/2.) - 1. );

  }

  double eta = 0.;
  FloatMatrix TN, dzeta(3,3), ksi(3,3), theta(3,3);

  TN.beProductOf(SethHillStress,N);
  dzeta.beTProductOf(N, TN);
    




  // compute auxiliary variables 
  // the computation differes depends on if the eigenvalues of C are equal or not
  if(lam.at(1) != lam.at(2)) {
    // all eigenvalues are different
    for(int i = 1; i <= 2; i++) {
      for(int j = 1; j <= 2; j++) {
	if(i == j) {
	  continue;
	} else {
	  theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
	  ksi.at(i,j) = (theta.at(i,j) - 1./2.*d.at(j))/(lam.at(i) - lam.at(j));
	  for(int k = 1; k <= 2; k++) {
	    if((k != i) && (k != j)) {
	      eta += eps.at(i)/2./(lam.at(i) - lam.at(j))/(lam.at(i)-lam.at(k));
	    }
	  }
	}
      }
    }
  } else {  // l1 == l2 == l3
    for(int i = 1; i <= 2; i++) {
      for(int j = 1; j <= 2; j++) {
	theta.at(i,j) = 1./2. * d.at(i);
	ksi.at(i,j) = 1./8. * f.at(i);
      }
    }
    eta = 1./8. * f.at(1);
  }


  FloatMatrix M(9,9);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
	for (int l = 1; l <=3; l++) {
	  M.at(giveVI(i,j),giveVI(k,l)) = n.at(k,i)*N.at(l,j) + n.at(k,j)*N.at(l,i);
	}
      }
    }
  }
    
  PP.resize(6,9), TL.resize(9,9);
  for (int k = 1; k <= 2; k++) {
    for (int l = 1; l <= 2; l++) {
      for (int m = 1; m <=2; m++) {
	for (int n = 1; n<=2; n++) {
	  for (int i = 1; i <= 2; i++) {
	    if(k == l ) {
	      PP.at(giveSymVI(k,l),giveVI(m,n)) += 0.5 * d.at(i) * N.at(k,i) * N.at(l,i) * M.at(giveVI(i,i),giveVI(m,n));
	    } else if(k < l) {
	      PP.at(giveSymVI(k,l),giveVI(m,n)) += d.at(i) * N.at(k,i) * N.at(l,i) * M.at(giveVI(i,i),giveVI(m,n));
	    }

	    TL.at(giveVI(k,l),giveVI(m,n)) += 1./4.*f.at(i)*dzeta.at(i,i)*M.at(giveVI(i,i),giveVI(k,l))*M.at(giveVI(i,i),giveVI(m,n));
	    for (int j  = 1; j <= 2; j++) {
	      if(j != i) {
		if( k == l) {
		  PP.at(giveSymVI(k,l),giveVI(m,n)) += theta.at(i,j)*N.at(k,i)*N.at(l,j)*M.at(giveVI(i,j),giveVI(m,n));
		} else if(k < l) {
		  PP.at(giveSymVI(k,l),giveVI(m,n)) += 2. * theta.at(i,j)*N.at(k,i)*N.at(l,j)*M.at(giveVI(i,j),giveVI(m,n));
		}
		TL.at(giveVI(k,l),giveVI(m,n)) += 2.*ksi.at(i,j) * ( dzeta.at(i,j)*(M.at(giveVI(i,j),giveVI(k,l))*M.at(giveVI(j,j),giveVI(m,n)) + M.at(giveVI(j,j),giveVI(k,l))*M.at(giveVI(i,j),giveVI(m,n))) + dzeta.at(j,j)*M.at(giveVI(i,j),giveVI(k,l))*M.at(giveVI(i,j),giveVI(m,n) ));
		for (int q = 1; q <= 2; q++) {
		  if(q == i) {
		    continue;
		  } else if(q == j) {
		    continue;
		  } else {
		    TL.at(giveVI(k,l),giveVI(m,n)) += 2.*eta*dzeta.at(i,j)*M.at(giveVI(i,q),giveVI(k,l))*M.at(giveVI(j,q),giveVI(m,n));
		  }
		}
	      }
	    }
	  }
	}	    
      }
    }
  }   
  
}









void
LargeStrainMasterMaterial2 :: giveTransformationMatrices_huhu(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &F, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N)
{  

  FloatMatrix n;
  n.beProductOf(F,N);
    
  FloatArray  d(3),f(3), eps(3);
  if(this->m == 0) { // log strain formulation (m = 0) has a special treatment
    d.at(1)= 1./lam.at(1);
    d.at(2)= 1./lam.at(2);
    d.at(3)= 1./lam.at(3);
    
    f.at(1) = -2./lam.at(1)/lam.at(1);
    f.at(2) = -2./lam.at(2)/lam.at(2);
    f.at(3) = -2./lam.at(3)/lam.at(3);
    
    eps.at(1) = log(lam.at(1))/2.;
    eps.at(2) = log(lam.at(2))/2.;
    eps.at(3) = log(lam.at(3))/2.;
    
  } else {  // all other strain measures (m != 0)
    d.at(1)= pow(lam.at(1), m/2.-1);
    d.at(2)= pow(lam.at(2), m/2.-1);
    d.at(3)= pow(lam.at(3), m/2.-1);
      
    f.at(1) = 2. * ( m/2. - 1 ) * pow( lam.at(1), m/2. - 2. );
    f.at(2) = 2. * ( m/2. - 1 ) * pow( lam.at(2), m/2. - 2. );
    f.at(3) = 2. * ( m/2. - 1 ) * pow( lam.at(3), m/2. - 2. );
      
    eps.at(1) = 1./m * ( pow(lam.at(1),m/2.) - 1. );
    eps.at(2) = 1./m * ( pow(lam.at(2),m/2.) - 1. );
    eps.at(3) = 1./m * ( pow(lam.at(3),m/2.) - 1. );

  }

  double eta = 0.;
  FloatMatrix TN, dzeta(3,3), ksi(3,3), theta(3,3);

  TN.beProductOf(SethHillStress,N);
  dzeta.beTProductOf(N, TN);
    




  // compute auxiliary variables 
  // the computation differes depends on if the eigenvalues of C are equal or not
  if(lam.at(1) != lam.at(2)) {
    // all eigenvalues are different
    for(int i = 1; i <= 2; i++) {
      for(int j = 1; j <= 2; j++) {
	if(i == j) {
	  continue;
	} else {
	  theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
	  ksi.at(i,j) = (theta.at(i,j) - 1./2.*d.at(j))/(lam.at(i) - lam.at(j));
	  for(int k = 1; k <= 2; k++) {
	    if((k != i) && (k != j)) {
	      eta += eps.at(i)/2./(lam.at(i) - lam.at(j))/(lam.at(i)-lam.at(k));
	    }
	  }
	}
      }
    }
  } else {  // l1 == l2 == l3
    for(int i = 1; i <= 2; i++) {
      for(int j = 1; j <= 2; j++) {
	theta.at(i,j) = 1./2. * d.at(i);
	ksi.at(i,j) = 1./8. * f.at(i);
      }
    }
    eta = 1./8. * f.at(1);
  }


  FloatMatrix M(9,9);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
	for (int l = 1; l <=3; l++) {
	  M.at(giveVI(i,j),giveVI(k,l)) = n.at(k,i)*N.at(l,j) + n.at(k,j)*N.at(l,i);
	}
      }
    }
  }
    
  PP.resize(6,9), TL.resize(9,9);
  for (int k = 1; k <= 3; k++) {
    for (int l = 1; l <= 3; l++) {
      for (int m = 1; m <=3; m++) {
	for (int n = 1; n<=3; n++) {
	  for (int i = 1; i <= 3; i++) {
	    if(k == l ) {
	      PP.at(giveSymVI(k,l),giveVI(m,n)) += 0.5 * d.at(i) * N.at(k,i) * N.at(l,i) * M.at(giveVI(i,i),giveVI(m,n));
	    } else if(k < l) {
	      PP.at(giveSymVI(k,l),giveVI(m,n)) += d.at(i) * N.at(k,i) * N.at(l,i) * M.at(giveVI(i,i),giveVI(m,n));
	    }

	    //  TL.at(giveVI(k,l),giveVI(m,n)) += 1./4.*f.at(i)*dzeta.at(i,i)*M.at(giveVI(i,i),giveVI(k,l))*M.at(giveVI(i,i),giveVI(m,n));
	    for (int j  = 1; j <= 3; j++) {
	      if(j != i) {
		if( k == l) {
		  PP.at(giveSymVI(k,l),giveVI(m,n)) += theta.at(i,j)*N.at(k,i)*N.at(l,j)*M.at(giveVI(i,j),giveVI(m,n));
		} else if(k < l) {
		  PP.at(giveSymVI(k,l),giveVI(m,n)) += 2. * theta.at(i,j)*N.at(k,i)*N.at(l,j)*M.at(giveVI(i,j),giveVI(m,n));
		}
		//	TL.at(giveVI(k,l),giveVI(m,n)) += 2.*ksi.at(i,j) * ( dzeta.at(i,j)*(M.at(giveVI(i,j),giveVI(k,l))*M.at(giveVI(j,j),giveVI(m,n)) + M.at(giveVI(j,j),giveVI(k,l))*M.at(giveVI(i,j),giveVI(m,n))) + dzeta.at(j,j)*M.at(giveVI(i,j),giveVI(k,l))*M.at(giveVI(i,j),giveVI(m,n) ));
		for (int q = 1; q <= 3; q++) {
		  if(q == i) {
		    continue;
		  } else if(q == j) {
		    continue;
		  } else {
		    TL.at(giveVI(k,l),giveVI(m,n)) += 2.*eta*dzeta.at(i,j)*M.at(giveVI(i,q),giveVI(k,l))*M.at(giveVI(j,q),giveVI(m,n));
		  }
		}
	      }
	    }
	  }
	}	    
      }
    }
  }   
  
}










  

void
LargeStrainMasterMaterial2 :: giveTransformationMatrices_PlaneStress_dSdE(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N)
{  

    
  FloatArray  d(3),f(3), eps(3);
  if(this->m == 0) { // log strain formulation (m = 0) has a special treatment
    d.at(1)= 1./lam.at(1);
    d.at(2)= 1./lam.at(2);
    d.at(3)= 1./lam.at(3);
    
    f.at(1) = -2./lam.at(1)/lam.at(1);
    f.at(2) = -2./lam.at(2)/lam.at(2);
    f.at(3) = -2./lam.at(3)/lam.at(3);
    
    eps.at(1) = log(lam.at(1))/2.;
    eps.at(2) = log(lam.at(2))/2.;
    eps.at(3) = log(lam.at(3))/2.;
    
  } else {  // all other strain measures (m != 0)
    d.at(1)= pow(lam.at(1), m/2.-1);
    d.at(2)= pow(lam.at(2), m/2.-1);
    d.at(3)= pow(lam.at(3), m/2.-1);
      
    f.at(1) = 2. * ( m/2. - 1 ) * pow( lam.at(1), m/2. - 2. );
    f.at(2) = 2. * ( m/2. - 1 ) * pow( lam.at(2), m/2. - 2. );
    f.at(3) = 2. * ( m/2. - 1 ) * pow( lam.at(3), m/2. - 2. );
      
    eps.at(1) = 1./m * ( pow(lam.at(1),m/2.) - 1. );
    eps.at(2) = 1./m * ( pow(lam.at(2),m/2.) - 1. );
    eps.at(3) = 1./m * ( pow(lam.at(3),m/2.) - 1. );

  }

  double eta = 0.;
  FloatMatrix TN, dzeta(3,3), ksi(3,3), theta(3,3);

  TN.beProductOf(SethHillStress,N);
  dzeta.beTProductOf(N, TN);
    




  // compute auxiliary variables 
  // the computation differes depends on if the eigenvalues of C are equal or not
  if(lam.at(1) != lam.at(2)) {
    // all eigenvalues are different
    for(int i = 1; i <= 2; i++) {
      for(int j = 1; j <= 2; j++) {
	if(i == j) {
	  continue;
	} else {
	  theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
	  ksi.at(i,j) = (theta.at(i,j) - 1./2.*d.at(j))/(lam.at(i) - lam.at(j));
	  for(int k = 1; k <= 2; k++) {
	    if((k != i) && (k != j)) {
	      eta += eps.at(i)/2./(lam.at(i) - lam.at(j))/(lam.at(i)-lam.at(k));
	    }
	  }
	}
      }
    }
  } else {  // l1 == l2 == l3
    for(int i = 1; i <= 2; i++) {
      for(int j = 1; j <= 2; j++) {
	theta.at(i,j) = 1./2. * d.at(i);
	ksi.at(i,j) = 1./8. * f.at(i);
      }
    }
    eta = 1./8. * f.at(1);
  }


  FloatMatrix M(9,9);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
	for (int l = 1; l <=3; l++) {
	  M.at(giveVI(i,j),giveVI(k,l)) = N.at(k,i)*N.at(l,j) + N.at(k,j)*N.at(l,i);
	}
      }
    }
  }
    
  PP.resize(9,9), TL.resize(9,9);
  for (int k = 1; k <= 3; k++) {
    for (int l = 1; l <= 3; l++) {
      for (int m = 1; m <=3; m++) {
	for (int n = 1; n<=3; n++) {
	  for (int i = 1; i <= 3; i++) {
	    PP.at(giveVI(k,l),giveVI(m,n)) += 0.5 * d.at(i) * N.at(k,i) * N.at(l,i) * M.at(giveVI(i,i),giveVI(m,n));
	    TL.at(giveVI(k,l),giveVI(m,n)) += 1./4.*f.at(i)*dzeta.at(i,i)*M.at(giveVI(i,i),giveVI(k,l))*M.at(giveVI(i,i),giveVI(m,n));
	    for (int j  = 1; j <= 3; j++) {
	      if(j != i) {
		PP.at(giveVI(k,l),giveVI(m,n)) += theta.at(i,j)*N.at(k,i)*N.at(l,j)*M.at(giveVI(i,j),giveVI(m,n));
		TL.at(giveVI(k,l),giveVI(m,n)) += 2.*ksi.at(i,j) * ( dzeta.at(i,j)*(M.at(giveVI(i,j),giveVI(k,l))*M.at(giveVI(j,j),giveVI(m,n)) + M.at(giveVI(j,j),giveVI(k,l))*M.at(giveVI(i,j),giveVI(m,n))) + dzeta.at(j,j)*M.at(giveVI(i,j),giveVI(k,l))*M.at(giveVI(i,j),giveVI(m,n) ));
		for (int q = 1; q <= 3; q++) {
		  if(q == i) {
		    continue;
		  } else if(q == j) {
		    continue;
		  } else {
		    TL.at(giveVI(k,l),giveVI(m,n)) += 2.*eta*dzeta.at(i,j)*M.at(giveVI(i,q),giveVI(k,l))*M.at(giveVI(j,q),giveVI(m,n));
		  }
		}
	      }
	    }
	  }
	}	    
      }
    }
  }   
  
}


  

void
LargeStrainMasterMaterial2 :: giveTransformationMatrices_PlaneStress2(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &F, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N)
{   

    FloatMatrix n;
    n.beProductOf(F,N);
    FloatArray  d(2),f(2), eps(2);
    if(this->m == 0) { // log strain formulation (m = 0) has a special treatment
      d.at(1)= 1. / lam.at(1);
      d.at(2)= 1. / lam.at(2);
      
      f.at(1) = -2. / lam.at(1) / lam.at(1);
      f.at(2) = -2. / lam.at(2) / lam.at(2);

      eps.at(1) = 0.5 * log( lam.at(1) );
      eps.at(2) = 0.5 * log( lam.at(2) );
    } else {  // all other strain measures (m != 0)
      d.at(1)= pow( lam.at(1), m / 2. - 1 );
      d.at(2)= pow( lam.at(2), m / 2. - 1 );
      
      f.at(1) = 2. * ( m / 2. - 1 ) * pow( lam.at(1), m / 2. - 2. );
      f.at(2) = 2. * ( m / 2. - 1 ) * pow( lam.at(2), m / 2. - 2. );
      
      eps.at(1) = 1. / m * ( pow( lam.at(1), m / 2. ) - 1. );
      eps.at(2) = 1. / m * ( pow( lam.at(2), m / 2. ) - 1. );
    }

    //    double eta = 0.;
    FloatArray vSethHillStress(6);
    vSethHillStress.beVectorForm(SethHillStress);
    FloatMatrix TN, dzeta(3,3), ksi(3,3), theta(3,3);
    // check dzeta
    TN.beProductOf(SethHillStress,N);
    dzeta.beTProductOf(N, TN);
    
    double beta;
    FloatArray gamma(2), kappa(2);
    // compute auxiliary variables 
    // the computation differes depends on if the eigenvalues of C are equal or not
    if(lam.at(1) != lam.at(2)) {
      // all eigenvalues are different
      beta = 2. * ( eps.at(1) - eps.at(2) ) / ( lam.at(1) - lam.at(2) );
      gamma.at(1) = d.at(1) - beta;
      gamma.at(2) = d.at(2) - beta;
      gamma.at(3) = d.at(3) - beta;
      kappa.at(1) = 2. * gamma.at(1) / ( lam.at(1) - lam.at(2) );
      kappa.at(2) = 2. * gamma.at(2) / ( lam.at(1) - lam.at(2) );
      
    } else { //l1 == l2
      beta = d.at(1);
      gamma.at(1) = d.at(1) - beta;
      gamma.at(2) = d.at(2) - beta;
      gamma.at(3) = d.at(3) - beta;
      kappa.at(1) = ( 1.5 - 1. ) * f.at(1);
      kappa.at(2) = ( 1.5 - 2. ) * f.at(1);


    }


    FloatMatrix I(9,9);
    FloatMatrix M(9,9);
    FloatArray delta(6);
    delta.zero();
    delta.at(1) = 1;
    delta.at(2) = 1;
    delta.at(3) = 1;
    for (int k = 1; k <= 3; k++) {
      for (int l = 1; l <= 3; l++) {
	for (int o = 1; o <= 3; o++) {
	  for (int p = 1; p <= 3; p++) {
	    I.at(giveVI(k,l),giveVI(o,p)) = 0.5 * ( delta.at(giveSymVI(k,p)) * F.at(o,l) + delta.at(giveSymVI(l,p)) * F.at(o,k) ) ;
	    M.at(giveVI(k,l),giveVI(o,p)) = n.at(o,k)*N.at(p,l) + n.at(o,l)*N.at(p,k);
	  }
	}
      }
    }
	    
    

    
    PP.resize(6,9), TL.resize(9,9);
    for (int k = 1; k <= 3; k++) {
      for (int l = 1; l <= 3; l++) {
	for (int m = 1; m <=3; m++) {
	  for (int n = 1; n<=3; n++) {
	    if(k <= l) {
	      PP.at(giveSymVI(k,l),giveVI(m,n)) += beta * I.at(giveVI(k,l),giveVI(m,n));
	    }
	    for (int i = 1; i <= 3; i++) {
	      // first  T:L term
	      TL.at(giveVI(k,l),giveVI(m,n)) += f.at(i)*dzeta.at(i,i)*M.at(giveVI(i,i),giveVI(k,l))*M.at(giveVI(i,i),giveVI(m,n));
	      // end of first  T:L term
	      // second  T:L term
	      if(i <= 2) {
		TL.at(giveVI(k,l),giveVI(m,n)) += pow( -1, i+1 ) * kappa.at(i) * ( dzeta.at(i,i) * I.at(giveVI(k,l),giveVI(m,n)) + vSethHillStress.at(giveVI(k,l)) * M.at(giveVI(i,i),giveVI(m,n)) +  M.at(giveVI(i,i),giveVI(k,l)) * vSethHillStress.at(giveVI(m,n) ) );
	      }
	      //end of second T:L term
	      if( k <= l ) {
		PP.at(giveSymVI(k,l),giveVI(m,n)) += gamma.at(i)*N.at(k,i)*N.at(l,i)*M.at(giveVI(i,i),giveVI(m,n));
	      }
	      for (int j = 1; j <= 2; j++) {      
		// third T:L term
		TL.at(giveVI(k,l),giveVI(m,n)) += pow( -1, j+1 ) * kappa.at(j) * ( dzeta.at(j,j) * M.at(giveVI(i,i),giveVI(k,l))*M.at(giveVI(i,i),giveVI(m,n))  + dzeta.at(i,i) *  M.at(giveVI(j,j),giveVI(k,l))*M.at( giveVI(i,i),giveVI(m,n) ) + dzeta.at(i,i) *  M.at(giveVI(i,i),giveVI(k,l) ) * M.at(giveVI(j,j),giveVI(m,n) ) );
		//end of second T:L term
	      }
	    }
	  }
	}	    
      }
    }

}


void
LargeStrainMasterMaterial2 :: giveTransformationMatrices_PlaneStress2_dSdE(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N)
{   
  
    FloatArray  d(3),f(3), eps(3);
    if(this->m == 0) { // log strain formulation (m = 0) has a special treatment
      d.at(1)= 1. / lam.at(1);
      d.at(2)= 1. / lam.at(2);
      d.at(3)= 1. / lam.at(3);
      
      f.at(1) = -2. / lam.at(1) / lam.at(1);
      f.at(2) = -2. / lam.at(2) / lam.at(2);
      f.at(3) = -2. / lam.at(3) / lam.at(3);

      eps.at(1) = 0.5 * log( lam.at(1) );
      eps.at(2) = 0.5 * log( lam.at(2) );
      eps.at(3) = 0.5 * log( lam.at(3) );
    } else {  // all other strain measures (m != 0)
      d.at(1)= pow( lam.at(1), m / 2. - 1 );
      d.at(2)= pow( lam.at(2), m / 2. - 1 );
      d.at(3)= pow( lam.at(3), m / 2. - 1 );
      
      f.at(1) = 2. * ( m / 2. - 1 ) * pow( lam.at(1), m / 2. - 2. );
      f.at(2) = 2. * ( m / 2. - 1 ) * pow( lam.at(2), m / 2. - 2. );
      f.at(3) = 2. * ( m / 2. - 1 ) * pow( lam.at(3), m / 2. - 2. );
      
      eps.at(1) = 1. / m * ( pow( lam.at(1), m / 2. ) - 1. );
      eps.at(2) = 1. / m * ( pow( lam.at(2), m / 2. ) - 1. );
      eps.at(3) = 1. / m * ( pow( lam.at(3), m / 2. ) - 1. );
    }

    //    double eta = 0.;
    FloatArray vSethHillStress(6);
    vSethHillStress.beVectorForm(SethHillStress);
    FloatMatrix TN, dzeta(3,3), ksi(3,3), theta(3,3);
    // check dzeta
    TN.beProductOf(SethHillStress,N);
    dzeta.beTProductOf(N, TN);
    
    double beta;
    FloatArray gamma(3), kappa(2);
    // compute auxiliary variables 
    // the computation differes depends on if the eigenvalues of C are equal or not
    if(lam.at(1) != lam.at(2)) {
      // all eigenvalues are different
      beta = 2. * ( eps.at(1) - eps.at(2) ) / ( lam.at(1) - lam.at(2) );
      gamma.at(1) = d.at(1) - beta;
      gamma.at(2) = d.at(2) - beta;
      gamma.at(3) = d.at(3) - beta;
      kappa.at(1) = 2. * gamma.at(1) / ( lam.at(1) - lam.at(2) );
      kappa.at(2) = 2. * gamma.at(2) / ( lam.at(1) - lam.at(2) );
      
    } else { //l1 == l2
      beta = d.at(1);
      gamma.at(1) = d.at(1) - beta;
      gamma.at(2) = d.at(2) - beta;
      gamma.at(3) = d.at(3) - beta;
      kappa.at(1) = ( 1.5 - 1. ) * f.at(1);
      kappa.at(2) = ( 1.5 - 2. ) * f.at(1);


    }


    FloatMatrix I(9,9);
    FloatMatrix M(9,9);
    FloatArray delta(6);
    delta.zero();
    delta.at(1) = 1;
    delta.at(2) = 1;
    delta.at(3) = 1;
    for (int k = 1; k <= 3; k++) {
      for (int l = 1; l <= 3; l++) {
	for (int o = 1; o <= 3; o++) {
	  for (int p = 1; p <= 3; p++) {
	    I.at(giveVI(k,l),giveVI(o,p)) = 0.5 * ( delta.at(giveSymVI(k,p)) * delta.at(giveSymVI(o,l)) + delta.at(giveSymVI(l,p)) * delta.at(giveSymVI(l,p)) );
	    M.at(giveVI(k,l),giveVI(o,p)) = N.at(o,k)*N.at(p,l) + N.at(o,l)*N.at(p,k);
	  }
	}
      }
    }
	    
    

    
    PP.resize(6,9), TL.resize(9,9);
    for (int k = 1; k <= 3; k++) {
      for (int l = 1; l <= 3; l++) {
	for (int m = 1; m <=3; m++) {
	  for (int n = 1; n<=3; n++) {
	    if(k <= l) {
	      PP.at(giveSymVI(k,l),giveVI(m,n)) += beta * I.at(giveVI(k,l),giveVI(m,n));
	    }
	    for (int i = 1; i <= 3; i++) {
	      // first  T:L term
	      TL.at(giveVI(k,l),giveVI(m,n)) += f.at(i)*dzeta.at(i,i)*M.at(giveVI(i,i),giveVI(k,l))*M.at(giveVI(i,i),giveVI(m,n));
	      // end of first  T:L term
	      // second  T:L term
	      if(i <= 2) {
		TL.at(giveVI(k,l),giveVI(m,n)) += pow( -1, i+1 ) * kappa.at(i) * ( dzeta.at(i,i) * I.at(giveVI(k,l),giveVI(m,n)) + vSethHillStress.at(giveVI(k,l)) * M.at(giveVI(i,i),giveVI(m,n)) +  M.at(giveVI(i,i),giveVI(k,l)) * vSethHillStress.at(giveVI(m,n) ) );
	      }
	      //end of second T:L term
	      if( k <= l ) {
		PP.at(giveSymVI(k,l),giveVI(m,n)) += gamma.at(i)*N.at(k,i)*N.at(l,i)*M.at(giveVI(i,i),giveVI(m,n));
	      }
	      for (int j = 1; j <= 2; j++) {      
		// third T:L term
		TL.at(giveVI(k,l),giveVI(m,n)) += pow( -1, j+1 ) * kappa.at(j) * ( dzeta.at(j,j) * M.at(giveVI(i,i),giveVI(k,l))*M.at(giveVI(i,i),giveVI(m,n))  + dzeta.at(i,i) *  M.at(giveVI(j,j),giveVI(k,l))*M.at( giveVI(i,i),giveVI(m,n) ) + dzeta.at(i,i) *  M.at(giveVI(i,i),giveVI(k,l) ) * M.at(giveVI(j,j),giveVI(m,n) ) );
		//end of second T:L term
	      }
	    }
	  }
	}	    
      }
    }

}

  


void
LargeStrainMasterMaterial2 :: giveTransformationMatrices_PlaneStressAnalyticalM_2(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &Fm, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N)
{

  FloatMatrix iC, invF,invFm, F;
  FloatArray n, vN;
  vN = {0,0,1};
  double l3 = Fm.at(3,3);
  invFm.beInverseOf(Fm);
  n.beProductOf(invFm, vN);
  n.normalize();
  F = Fm;
  F.at(1,3) = l3 * n.at(1);
  F.at(2,3) = l3 * n.at(2);
  F.at(3,3) = l3 * n.at(3);
  invF.beInverseOf(F);
  iC.beProductTOf(invF,invF);
   
  
  PP.resize(6,9), TL.resize(9,9);
  for (int k = 1; k <= 3; k++) {
    for (int l = 1; l <= 3; l++) {
      for (int p = 1; p <= 3; p++) {
	for (int q = 1; q <= 3; q++) {
	  for (int i = 1; i <=3; i++) {
	    for (int j = 1; j <=3; j++) {
	      if((giveSymVI(i,j) > 3) & (giveVI(k,l) > 3)) {
		PP.at(giveSymVI(i,j),giveVI(k,l)) = ( invF.at(i,k) * iC.at(j,l) + invF.at(j,k)*iC.at(i,l) );
	      } else {
		PP.at(giveSymVI(i,j),giveVI(k,l)) = 0.5 * ( invF.at(i,k) * iC.at(j,l) + invF.at(j,k)*iC.at(i,l) );
	      }
	      for (int n = 1; n <= 3; n++) {
		TL.at(giveVI(k,l),giveVI(p,q)) -= 0.5 * ( invF.at(i,p) * invF.at(q,k) * invF.at(l,n) * invF.at(j,n) + invF.at(i,k)  * invF.at(l,p) * invF.at(q,n) * invF.at(j,n) + invF.at(i,k) * invF.at(l,n) * invF.at(j,p) * invF.at(q,n) + invF.at(i,p)*invF.at(q,n)*invF.at(j,k)*invF.at(l,n) + invF.at(i,n)*invF.at(j,p)*invF.at(q,k)*invF.at(l,n) + invF.at(i,n)*invF.at(j,k)*invF.at(l,p)*invF.at(q,n) ) * SethHillStress.at(i,j);
	      }
	    }
	  }
	}
      }
    }
  }
  
      
      
}



void
LargeStrainMasterMaterial2 :: giveTransformationMatrices_PlaneStressAnalyticalM_2_dSdE(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &C, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N)
{

  FloatMatrix iC;
  iC.beInverseOf(C);
  
  
  PP.resize(9,9), TL.resize(9,9);
  for (int k = 1; k <= 3; k++) {
    for (int l = 1; l <= 3; l++) {
      for (int p = 1; p <= 3; p++) {
	for (int q = 1; q <= 3; q++) {
	  for (int i = 1; i <=3; i++) {
	    for (int j = 1; j <=3; j++) {
	      PP.at(giveVI(i,j),giveVI(k,l)) = 0.5 * ( iC.at(i,k) * iC.at(j,l) + iC.at(j,k)*iC.at(i,l) );
	      TL.at(giveVI(k,l),giveVI(p,q)) -= 0.5 * ( ( iC.at(i,p) * iC.at(q,k) + iC.at(i,q) * iC.at(p,k) ) * iC.at(l,j) +  ( iC.at(l,p) * iC.at(q,j) + iC.at(l,q) * iC.at(p,j) ) * iC.at(i,k) +  ( iC.at(i,p) * iC.at(q,l) + iC.at(i,q) * iC.at(p,l) ) * iC.at(k,j) + ( iC.at(k,p) * iC.at(q,j) + iC.at(k,q) * iC.at(p,j) ) * iC.at(i,l) ) * SethHillStress.at(i,j);
	    }
	  }
	}
      }
    }
  }    
      
}  




  
  

int
LargeStrainMasterMaterial2 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    LargeStrainMasterMaterial2Status *status = static_cast< LargeStrainMasterMaterial2Status * >( this->giveStatus(gp) );


    GaussPoint *slaveGp = status->giveSlaveGp();
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    return sMat->giveIPValue(answer, slaveGp, type, tStep);
}


//=============================================================================

  LargeStrainMasterMaterial2Status :: LargeStrainMasterMaterial2Status(int n, Domain *d, GaussPoint *g, int s) : StructuralMaterialStatus(n, d, g),
    Pmatrix(6, 6),
    TLmatrix(6, 6),
    transformationMatrix(6, 6),
    slaveMat(s)
{
    Pmatrix.beUnitMatrix();
    FloatArray coords = {0,0,0};
    slaveGp = new GaussPoint(NULL, 1, coords,1.0, g->giveMaterialMode() );
    //    slaveGp = new GaussPoint(NULL, 1, coords,1.0, _3dMat );
    
}


LargeStrainMasterMaterial2Status :: ~LargeStrainMasterMaterial2Status()
{ }


void
LargeStrainMasterMaterial2Status :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    MaterialStatus *mS = sMat->giveStatus(slaveGp);

    mS->printOutputAt(file, tStep);
    //  StructuralMaterialStatus :: printOutputAt(file, tStep);
}


// initializes temporary variables based on their values at the previous equlibrium state
void LargeStrainMasterMaterial2Status :: initTempStatus()
{
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    MaterialStatus *mS = sMat->giveStatus(slaveGp);
    mS->initTempStatus();

    if ( this->giveCVector().giveSize() == 0 ) {
      stressVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
    }

    StructuralMaterialStatus :: initTempStatus();
    tempCVector      = CVector;
      
    //
}


// updates internal variables when equilibrium is reached
void
LargeStrainMasterMaterial2Status :: updateYourself(TimeStep *tStep)
{
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    MaterialStatus *mS = sMat->giveStatus(slaveGp);
    mS->updateYourself(tStep);
    StructuralMaterialStatus :: updateYourself(tStep);
    CVector      = tempCVector;
	
    
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
LargeStrainMasterMaterial2Status :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    StructuralMaterial *sMat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    MaterialStatus *mS = sMat->giveStatus(slaveGp);
    // save parent class status
    if ( ( iores = mS->saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data

    return CIO_OK;
}


contextIOResultType
LargeStrainMasterMaterial2Status :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem




