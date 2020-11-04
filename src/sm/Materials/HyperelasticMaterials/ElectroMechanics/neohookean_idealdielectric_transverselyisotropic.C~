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

#include "../sm/Materials/HyperelasticMaterials/ElectroMechanics/mooneyrivlin_idealdielectric.h"
#include "../sm/Materials/HyperelasticMaterials/mooneyrivlin.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"




namespace oofem {
  REGISTER_Material(MooneyRivlin_IdealDielectricMaterial);

  MooneyRivlin_IdealDielectricMaterial :: MooneyRivlin_IdealDielectricMaterial(int n, Domain *d):Material(n, d), ElectroMechanicalMaterialExtensionInterface(d), ElectroMechanicalMaterialExtensionInterface_3Field(d)
{
    this->hyperelasticMaterial = new MooneyRivlinMaterial(n, d);
}


IRResultType
MooneyRivlin_IdealDielectricMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro


    result = this->hyperelasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) {
      return result;
    }

    IR_GIVE_FIELD(ir, this->epsilon, _IFT_MooneyRivlin_IdealDielectricMaterial_epsilon);
      
    return IRRT_OK;
}

    
void
MooneyRivlin_IdealDielectricMaterial :: give_FirstPKStressVector_ElectricalDisplacementVector_3d(FloatArray &vP, FloatArray &D, GaussPoint *gp, const FloatArray &vF, const FloatArray &E, TimeStep *tStep)
{
    // First Piol-Kirchhoff stress
    double Sigma_J;
    FloatArray vSigma_H;
    FloatMatrix F, Sigma_H;
    hyperelasticMaterial->giveFirstPKStressVector_3d(vP, gp, vF, tStep);    
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );

    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    //
    FloatMatrix H;
    H.beMatrixForm(vH);
    // Sigma_J
    FloatArray H_E, vP_Sigma_H;
    H_E.beProductOf(H, E);
    Sigma_J = -  0.5 * epsilon / J / J * H_E.dotProduct( H_E );
    // Sigma_H
    Sigma_H.beDyadicProductOf( H_E, E );
    Sigma_H.times( epsilon / J );
    vSigma_H.beVectorForm(Sigma_H);
    hyperelasticMaterial->compute_2order_tensor_cross_product( vP_Sigma_H, vSigma_H, vF);
    //
    vP.add( Sigma_J, vH );
    vP.add( vP_Sigma_H );
    

    // electric displacement
    FloatMatrix C, G;
    F.beMatrixForm(vF);
    C.beTProductOf(F,F);
    hyperelasticMaterial->compute_2order_tensor_cross_product( G, C, C );
    G.times(0.5);
    D.beProductOf(G, E);
    D.times( epsilon * J);
    
    
    status->letTempPVectorBe(vP);
    status->letTempFVectorBe(vF);
    status->letTempEVectorBe(E);
    status->letTempDVectorBe(D);
}

  

void
MooneyRivlin_IdealDielectricMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
 
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    hyperelasticMaterial->give3dMaterialStiffnessMatrix_dPdF(answer, mode, gp, tStep); 
    FloatArray vF = status->giveTempFVector();
    FloatArray E = status->giveTempEVector();

    // First Piola-Kirchhoff stress
    double Sigma_J;
    FloatArray vSigma_H;
    FloatMatrix F, Sigma_H;
    // Kronecker delta
    FloatMatrix delta(3,3);
    delta.beUnitMatrix();
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    //
    FloatMatrix H;
    H.beMatrixForm(vH);
    // Sigma_J
    FloatArray H_E, vP_Sigma_H;
    H_E.beProductOf(H, E);
    Sigma_J = - 0.5 * epsilon / J / J * H_E.dotProduct( H_E );
    // Sigma_H
    Sigma_H.beDyadicProductOf( H_E, E );
    Sigma_H.times( epsilon / J );
    // dSigmaH_dH
    FloatMatrix dSigmaH_dH, EE;
    EE.beDyadicProductOf(E,E);
    hyperelasticMaterial->compute_lower_dyadic_product(dSigmaH_dH, delta, EE);
    dSigmaH_dH.times( epsilon / J );
    // dSigmaH_dJ
    FloatMatrix dSigmaH_dJ, HEE;   
    HEE.beProductOf(H, EE);
    dSigmaH_dJ = HEE;
    dSigmaH_dJ.times( - epsilon / J / J );
    // dSigmaJ_dH
    FloatMatrix dSigmaJ_dH;
    //HEE.beProductOf(H, EE);
    //    hyperelasticMaterial->compute_dyadic_product(dSigmaJ_dH, HEE, H);
    dSigmaJ_dH = HEE;
    dSigmaJ_dH.times( - epsilon / J / J );
    // dSigmaJ_dJ
    double dSigmaJ_dJ;
    FloatArray HE;
    HE.beProductOf(H,E);
    dSigmaJ_dJ = HE.dotProduct(HE);
    dSigmaJ_dJ *= epsilon / J / J /J;
    ////////////////////////////////////    
    //
    FloatMatrix answer1;
    hyperelasticMaterial->compute_tensor_cross_product_tensor(answer1, Sigma_H);
    //
    FloatMatrix Fx;
    F.beMatrixForm(vF);
    hyperelasticMaterial->compute_tensor_cross_product_tensor( Fx, F );
    FloatMatrix FxdSigmaH_dH;
    FxdSigmaH_dH.beProductOf( Fx, dSigmaH_dH );
    FloatMatrix answer2;
    answer2.beProductOf( FxdSigmaH_dH, Fx );
    //
    FloatMatrix dSigmaJ_dHxF;
    hyperelasticMaterial->compute_2order_tensor_cross_product( dSigmaJ_dHxF,  dSigmaJ_dH, F );
    FloatMatrix answer3;
    hyperelasticMaterial->compute_dyadic_product( answer3,  H, dSigmaJ_dHxF );   
    //
    FloatMatrix answer4;
    hyperelasticMaterial->compute_tensor_cross_product_tensor( answer4, F );
    answer4.times(Sigma_J);
    //
    FloatMatrix FxdSigmaH_dJ;
    hyperelasticMaterial->compute_2order_tensor_cross_product( FxdSigmaH_dJ,  F, dSigmaH_dJ );
    FloatMatrix answer5;
    hyperelasticMaterial->compute_dyadic_product(answer5, FxdSigmaH_dJ, H);
    //
    FloatMatrix answer6;
    hyperelasticMaterial->compute_dyadic_product(answer6, H, H);
    answer6.times(dSigmaJ_dJ);

    answer.add(answer1);
    answer.add(answer2);
    answer.add(answer3);
    answer.add(answer4);
    answer.add(answer5);
    answer.add(answer6);
}

void
MooneyRivlin_IdealDielectricMaterial :: give3dMaterialStiffnessMatrix_dPdE(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    FloatArray vF = status->giveTempFVector();
    FloatArray E = status->giveTempEVector();
    FloatMatrix F;
    F.beMatrixForm(vF);
    /* FloatMatrix Ft = {{3.560443081703946e+01,3.125513097186706e+00,2.532976909651305e+00},{-3.293652638284074e+00,3.804298360058239e+01,1.678442129540754e-01},{-7.993619384767280e+00,1.105875866534277e+01,3.476737600181701e+01}};
    F.beTranspositionOf(Ft);
    vF.beVectorForm(F);
    E = {2645.0914179429201, 7738.6133364184934, -3645.0914179429196};
    */
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    // compute deformation gradient
    double J = 1./3. * vH.dotProduct( vF );
    // dD_dJ
    FloatMatrix H, G, Ht;
    H.beMatrixForm( vH );
    Ht.beTranspositionOf( H );
    G.beTProductOf( H, H );
    FloatArray dD_dJ;
    dD_dJ.beProductOf( G, E );
    dD_dJ.times( - epsilon / J / J );
    // dD_dH
    FloatMatrix delta(3,3);
    delta.beUnitMatrix();
    FloatMatrix dD_dH, dD_dH1, dD_dH2;
    hyperelasticMaterial->compute_lower_dyadic_product( dD_dH1, delta, H );
    hyperelasticMaterial->compute_dyadic_product( dD_dH2, Ht, delta );
    dD_dH1.add(dD_dH2);
    hyperelasticMaterial->compute_dot_product(dD_dH, dD_dH1, E, 4);
    dD_dH.times( epsilon / J );
    //
    FloatMatrix answer1, answer2;
    hyperelasticMaterial->compute_3order_dyadic_product( answer1, dD_dJ, H );
    FloatMatrix Fx;
    hyperelasticMaterial->compute_tensor_cross_product_tensor( Fx, F );
    FloatMatrix dD_dH_39;
    hyperelasticMaterial->give_3order_tensor_39(dD_dH_39, dD_dH);
    FloatMatrix answer2_39;
    answer2_39.beProductOf( dD_dH_39, Fx );
    hyperelasticMaterial->give_3order_tensor_93(answer2, answer2_39);
    
    //
    answer.add(answer1);
    answer.add(answer2);
    answer.times(0);
    /*
    ///
    double pert = 1.e-6;
    FloatArray strain;
    FloatArray Ep(E);
    FloatArray D;

    FloatMatrix A(9,3);
    /// numerical stiffness
    FloatArray vP, vPp;
    //
    Ep.at(1) += pert;
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vPp, D, gp, vF, Ep, tStep);
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vP, D, gp, vF, E, tStep);
    A.at(1,1) = vPp.at(1) - vP.at(1);
    A.at(2,1) = vPp.at(2) - vP.at(2);
    A.at(3,1) = vPp.at(3) - vP.at(3);
    A.at(4,1) = vPp.at(4) - vP.at(4);
    A.at(5,1) = vPp.at(5) - vP.at(5);
    A.at(6,1) = vPp.at(6) - vP.at(6);
    A.at(7,1) = vPp.at(7) - vP.at(7);
    A.at(8,1) = vPp.at(8) - vP.at(8);
    A.at(9,1) = vPp.at(9) - vP.at(9);
    //
    vPp.zero();
    vP.zero();
    Ep = E;
    Ep.at(2) += pert;
        this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vPp, D, gp, vF, Ep, tStep);
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vP, D, gp, vF, E, tStep);
    A.at(1,2) = vPp.at(1) - vP.at(1);
    A.at(2,2) = vPp.at(2) - vP.at(2);
    A.at(3,2) = vPp.at(3) - vP.at(3);
    A.at(4,2) = vPp.at(4) - vP.at(4);
    A.at(5,2) = vPp.at(5) - vP.at(5);
    A.at(6,2) = vPp.at(6) - vP.at(6);
    A.at(7,2) = vPp.at(7) - vP.at(7);
    A.at(8,2) = vPp.at(8) - vP.at(8);
    A.at(9,2) = vPp.at(9) - vP.at(9);
    //
    Ep = E;
    Ep.at(3) += pert;
    vPp.zero();
    vP.zero();
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vPp, D, gp, vF, Ep, tStep);
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vP, D, gp, vF, E, tStep);
    A.at(1,3) = vPp.at(1) - vP.at(1);
    A.at(2,3) = vPp.at(2) - vP.at(2);
    A.at(3,3) = vPp.at(3) - vP.at(3);
    A.at(4,3) = vPp.at(4) - vP.at(4);
    A.at(5,3) = vPp.at(5) - vP.at(5);
    A.at(6,3) = vPp.at(6) - vP.at(6);
    A.at(7,3) = vPp.at(7) - vP.at(7);
    A.at(8,3) = vPp.at(8) - vP.at(8);
    A.at(9,3) = vPp.at(9) - vP.at(9);

    A.times(1./pert);
    answer.times(0);
    */


}

void
MooneyRivlin_IdealDielectricMaterial :: give3dMaterialStiffnessMatrix_dDdE(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    
    FloatArray vF = status->giveTempFVector();
    FloatMatrix F, C;
    F.beMatrixForm(vF);
    double J = F.giveDeterminant();
    C.beTProductOf(F,F);
    hyperelasticMaterial->compute_2order_tensor_cross_product( answer, C, C );
    answer.times( 0.5 * epsilon / J );


    /*
     */
    double pert = 1.e-6;
    FloatArray strain;
    FloatArray E = status->giveTempEVector();
    FloatArray Ep(E);
    FloatArray D, Dp;
    /*
    FloatMatrix A(3,3);
    /// numerical stiffness
    FloatArray vP;
    //
    Ep.at(1) += pert;
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vP, Dp, gp, vF, Ep, tStep);
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vP, D, gp, vF, E, tStep);
    A.at(1,1) = Dp.at(1) - D.at(1);
    A.at(2,1) = Dp.at(2) - D.at(2);
    A.at(3,1) = Dp.at(3) - D.at(3);
    //
    Dp.zero();
    D.zero();
    Ep = E;
    Ep.at(2) += pert;
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vP, Dp, gp, vF, Ep, tStep);
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vP, D, gp, vF, E, tStep);
    A.at(1,2) = Dp.at(1) - D.at(1);
    A.at(2,2) = Dp.at(2) - D.at(2);
    A.at(3,2) = Dp.at(3) - D.at(3);
    //
    Ep = E;
    Ep.at(3) += pert;
    Dp.zero();
    D.zero();
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vP, Dp, gp, vF, Ep, tStep);
    this->give_FirstPKStressVector_ElectricalDisplacementVector_3d(vP, D, gp, vF, E, tStep);
    A.at(1,3) = Dp.at(1) - D.at(1);
    A.at(2,3) = Dp.at(2) - D.at(2);
    A.at(3,3) = Dp.at(3) - D.at(3);

    A.times(1./pert);

    */

    
}


  ////////////////////////////////////////////////// Functions for 3-Field formulation

void
MooneyRivlin_IdealDielectricMaterial :: give3dMaterialStiffnessMatrix_dPdD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  /* ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    FloatArray vF = status->giveTempFVector();
    FloatArray D = status->giveTempDVector();
    FloatMatrix F;
    F.beMatrixForm(vF);
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    // compute deformation gradient
    double J = 1./3. * vH.dotProduct( vF );
    // dD_dJ
    FloatMatrix H, G, Ht;
    H.beMatrixForm( vH );
    Ht.beTranspositionOf( H );
    G.beTProductOf( H, H );
    FloatArray dD_dJ;
    dD_dJ.beProductOf( G, E );
    dD_dJ.times( - epsilon / J / J );
    // dD_dH
    FloatMatrix delta(3,3);
    delta.beUnitMatrix();
    FloatMatrix dD_dH, dD_dH1, dD_dH2;
    hyperelasticMaterial->compute_lower_dyadic_product( dD_dH1, delta, H );
    hyperelasticMaterial->compute_dyadic_product( dD_dH2, Ht, delta );
    dD_dH1.add(dD_dH2);
    hyperelasticMaterial->compute_dot_product(dD_dH, dD_dH1, E, 4);
    dD_dH.times( epsilon / J );
    //
    FloatMatrix answer1, answer2;
    hyperelasticMaterial->compute_3order_dyadic_product( answer1, dD_dJ, H );
    FloatMatrix Fx;
    hyperelasticMaterial->compute_tensor_cross_product_tensor( Fx, F );
    FloatMatrix dD_dH_39;
    hyperelasticMaterial->give_3order_tensor_39(dD_dH_39, dD_dH);
    FloatMatrix answer2_39;
    answer2_39.beProductOf( dD_dH_39, Fx );
    hyperelasticMaterial->give_3order_tensor_93(answer2, answer2_39);
    
    answer.add(answer1);
    answer.add(answer2);
  */

}

void
MooneyRivlin_IdealDielectricMaterial :: give3dMaterialStiffnessMatrix_dEdD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    
    FloatArray vF = status->giveTempFVector();
    FloatMatrix F, C;
    F.beMatrixForm(vF);
    double J = F.giveDeterminant();
    answer.beTProductOf(F,F);
    answer.times( 1. / epsilon /  J );

}





  
int
MooneyRivlin_IdealDielectricMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

  //  ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    
   
    return 1;
    
}
    
  


} // end namespace oofem
