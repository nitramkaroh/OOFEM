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

#include "../sm/Materials/HyperelasticMaterials/ElectroMechanics/mooneyrivlin_idealdielectric_transverselyisotropic2.h"
#include "../sm/Materials/structuralmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"

using sm = oofem::StructuralMaterial;


namespace oofem {
  REGISTER_Material(MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2);

  MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2(int n, Domain *d) : Material(n, d), ElectroMechanicalMaterialExtensionInterface_3Fields(d)
{
    this->hyperelasticMaterial = new MooneyRivlinMaterial(n, d);
}


IRResultType
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro


    this->alpha = 2;
    this->beta = 2;
    result = this->hyperelasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) {
      return result;
    }

    double epsilon1 = 0, epsilon2 = 0, epsilon3 = 0, epsilon4 = 0, epsilon5 = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon1, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon1);
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon2, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon2);
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon3, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon3);
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon4, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon4);
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon5, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon5);

    if(epsilon1 == 0) {
      this->iEps1 = 0.;
    } else {
      this->iEps1 = 1. / epsilon1;

    }

    if(epsilon2 == 0) {
      this->iEps2 = 0.;
    } else {
      this->iEps2 = 1. / epsilon2;

    }

    if(epsilon3 == 0) {
      this->iEps3 = 0.;
    } else {
      this->iEps3 = 1. / epsilon3;

    }

    if(epsilon4 == 0) {
      this->iEps4 = 0.;
    } else {
      this->iEps4 = 1. / epsilon4;

    }

    if(epsilon5 == 0) {
      this->iEps5 = 0.;
    } else {
      this->iEps5 = 1. / epsilon5;

    }

    IR_GIVE_FIELD(ir, this->N, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_n);
    N.normalize();
    IR_GIVE_FIELD(ir, this->mu1, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_mu1);
    IR_GIVE_FIELD(ir, this->mu2, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_mu2);
    IR_GIVE_FIELD(ir, this->mu3, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_mu3);
    IR_GIVE_FIELD(ir, this->lambda, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_lambda);
      
    return IRRT_OK;
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: give_FirstPKStressVector_ElectricalFieldVector_3d(FloatArray &vP, FloatArray &E, GaussPoint *gp, const FloatArray &vF, const FloatArray &D, TimeStep *tStep)
{

    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );   
    //////////////////////
    vP.zero();
    E.zero();
    FloatMatrix F;
    F.beMatrixForm(vF);
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
    ////////////////    
    // First Piol-Kirchhoff stress
    FloatArray dI1dF, dI2dF;
    this->compute_dI1dF(dI1dF, F);
    this->compute_dI2dF(dI2dF, F);
    ///////////////////////////////
    double I4 = this->compute_I4(F);
    double I5 = this->compute_I5(F);
    double I9 = this->compute_I9(D);
    double I10 = this->compute_I10(F,D);
    ///////////////////////////////////
    FloatArray dI4dF, dI5dF, dI7dF, dI8dF, dI10dF;
    this->compute_dI4dF(dI4dF, F);
    this->compute_dI5dF(dI5dF, F);
    this->compute_dI7dF(dI7dF, F, D);
    this->compute_dI8dF(dI7dF, F, D);
    this->compute_dI10dF(dI7dF, F, D);
    ///////////////////////////////////
    FloatArray dI6dD, dI7dD,dI8dD,dI9dD, dI10dD;
    this->compute_dI6dD(dI7dD, D);
    this->compute_dI7dD(dI7dD, F, D);
    this->compute_dI8dD(dI7dD, F, D);
    this->compute_dI9dD(dI9dD, D);
    this->compute_dI10dD(dI10dD, F, D);
    ///////////////////////////////////
    // isotropic mechanical invariants
    vP.add(0.5 * mu1, dI1dF);
    vP.add(0.5 * mu2, dI2dF);
    vP.add(0.5 * lambda * (J - 1.), vH);
    ///////////////////////////////////
    // anisotropic mechanical invariants
    vP.add(0.5 * mu3 * alpha * pow(I4, alpha - 1.), dI4dF);
    vP.add(0.5 * mu3 * beta  * pow(I5, beta  - 1.), dI5dF);
    vP.add(0.5 * mu3 / J, vH);
    ///////////////////////////////////
    // electro-mechanical invariants
    vP.add(   0.5 * this->iEps1 / J, dI7dF);
    vP.add( - 0.5 * this->iEps1 / J / J, vH);
    vP.add(   this->iEps5 * I10 , dI10dF);
    ///////////////////////////////////
    // mix of invariants I1 * I6
    auto I6 = this->compute_I6(D);
    vP.add(iEps4 * I6, vF);


    ////
    double FF = vF.dotProduct(vF);
    E.add( 0.5 * this->iEps2, dI6dD);    
    E.add( 0.5 * this->iEps1 / J , dI7dD);
    E.add( I9  * this->iEps3, dI9dD);    
    E.add( I10 * this->iEps5, dI10dD);    
    E.add( 0.5 * FF * this->iEps4, dI6dD);    
    
    status->letTempPVectorBe(vP);
    status->letTempFVectorBe(vF);
    status->letTempEVectorBe(E);
    status->letTempDVectorBe(D);
}


  
  

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
 
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    // hyperelasticMaterial->give3dMaterialStiffnessMatrix_dPdF(answer, mode, gp, tStep); 
    FloatArray vF = status->giveTempFVector();
    //  FloatArray E = status->giveTempEVector();
    FloatArray D = status->giveTempDVector();

    // First Piola-Kirchhoff stress
    FloatMatrix F;
    F.beMatrixForm(vF);
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
    FloatMatrix H, HH;
    H.beMatrixForm(vH);
    ////////////////////////////////////    
    FloatMatrix Fx;
    F.beMatrixForm(vF);
    hyperelasticMaterial->compute_tensor_cross_product_tensor( Fx, F );
    HH.beDyadicProductOf(vH, vH);
    ///////////////////////////////////
    auto I4 = this->compute_I4(F);
    auto I5 = this->compute_I5(F);
    auto I6 = this->compute_I6(D);
    auto I7 = this->compute_I7(F, D);
    auto I9 = this->compute_I9(D); 
    auto I10 = this->compute_I10(F, D);
    FloatArray dI4dF, dI5dF, dI7dF, dI10dF;
    FloatMatrix d2I1dF2, d2I2dF2,d2I4dF2, d2I5dF2, d2I7dF2, d2I10dF2;
    this->compute_d2I1dF2(d2I1dF2, F);
    this->compute_d2I2dF2(d2I2dF2, F);
    //   
    this->compute_dI4dF(dI4dF, F);
    this->compute_d2I4dF2(d2I4dF2, F);
    //
    this->compute_dI5dF(dI5dF, F);
    this->compute_d2I5dF2(d2I5dF2, F);
    //
    this->compute_dI7dF(dI7dF, F, D);
    this->compute_d2I7dF2(d2I7dF2, F, D);
    //
    this->compute_dI10dF(dI10dF, F, D);
    this->compute_d2I10dF2(d2I10dF2, F, D);
    // Isotropic stiffness - mechanical only
    answer.add(0.5 * mu1, d2I1dF2);
    answer.add(0.5 * mu2, d2I2dF2);
    answer.add( lambda * ( J - 1 ) ,Fx);
    answer.add( lambda , HH);
    // Anisotropic stiffness - mechanical only
    FloatMatrix dI4dI4, dI5dI5;
    dI4dI4.beDyadicProductOf(dI4dF, dI4dF);
    dI5dI5.beDyadicProductOf(dI5dF, dI5dF);
    answer.add( 0.5 * mu3 * pow(I4, alpha-1), d2I4dF2);
    answer.add( 0.5 * mu3 * (alpha - 1.) * pow(I4, alpha-2), dI4dI4);
    answer.add( 0.5 * mu3 * pow(I5, beta-1.), d2I5dF2);
    answer.add( 0.5 * mu3 * (beta - 1) * pow(I4, beta-2), dI5dI5);
    answer.add( 0.5 * mu3 / J / J, HH);
    answer.add( - 0.5 * mu3  / J, Fx);
    // Isotropic stiffness - electromechanical
    FloatMatrix HdI7, dI7H;
    dI7H.beDyadicProductOf(dI7dF, vH);
    HdI7.beDyadicProductOf(vH, dI7dF);
    answer.add(0.5 * iEps1 / J, d2I7dF2);
    answer.add( iEps1 * I7 / J / J/ J,  HH);
    answer.add( - 0.5 * iEps1 / J / J, dI7H);
    answer.add( - 0.5 * iEps1 / J / J, HdI7);
    // Anisotropic stiffness - electromechanical
    // FloatMatrix dI9dI9;
    // dI9dI9.beDyadicProductOf(dI9dF, dI9dF);
    //    answer.add( iEps3 * I9, d2I9dF2);
    //answer.add( iEps3 ,  dI9dI9);
    FloatMatrix dI10dI10;
    dI10dI10.beDyadicProductOf(dI10dF, dI10dF);
    answer.add( iEps4 * I10, d2I10dF2);
    answer.add( iEps4 ,  dI10dI10);
    // Mix of invariants I1 * I6
    FloatMatrix dd;
    hyperelasticMaterial->compute_lower_dyadic_product(dd, delta, delta);
    answer.add(iEps4 * I6, dd);
    
    FloatMatrix num_d2I4dF2,num_d2I7dF2,num_d2I10dF2;
    FloatArray num_dI10dF;
    this->compute_d2I4dF2_num(num_d2I4dF2, F);
    this->compute_d2I7dF2_num(num_d2I7dF2, F, D);
    this->compute_d2I10dF2_num(num_d2I10dF2, F, D);
    this->compute_dI10dF_num(num_dI10dF, F, D);

    /*
    
    FloatMatrix dPdF, dEdF;
    this->compute_dPdF_dEdF(dPdF, dEdF, F, D, gp, tStep);

    //this->compute_dPdF_dEdF_fake(dPdF, dEdF, F, D, gp, tStep);
    */

    
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: give3dMaterialStiffnessMatrix_dPdD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(9,3);
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    FloatArray vF = status->giveTempFVector();
    FloatArray D = status->giveTempDVector();
    // First Piola-Kirchhoff stress
    FloatMatrix F;
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
    ////////////////////////////////////    
    FloatMatrix Fx;
    F.beMatrixForm(vF);
    hyperelasticMaterial->compute_tensor_cross_product_tensor( Fx, F );
    ///////////////////////////////////
    auto I10 = this->compute_I10(F,D);
    FloatArray dI7dD, dI10dF, dI10dD;
    FloatMatrix d2I7dFdD, d2I10dFdD;
    this->compute_d2I7dFdD(d2I7dFdD, F, D);
    this->compute_d2I10dFdD(d2I10dFdD, F, D);
    this->compute_dI7dD(dI7dD, F, D);
    this->compute_dI10dF(dI10dF, F, D);
    this->compute_dI10dD(dI10dD, F, D);
    FloatMatrix  H_dI7dD;
    H_dI7dD.beDyadicProductOf(vH,dI7dD);
    //
    answer.add(0.5 * iEps1 / J, d2I7dFdD);
    answer.add(- 0.5 * iEps1 / J / J, H_dI7dD);
    //
    FloatMatrix  dI10dI10;
    dI10dI10.beDyadicProductOf(dI10dF,dI10dD);
    answer.add(iEps5, dI10dI10);
    answer.add(iEps5, d2I10dFdD);
    //
    FloatMatrix FoD;
    FoD.beDyadicProductOf(vF,D);
    answer.add(2. * iEps4, FoD);
    if(mode == SecantStiffnessMatrix) {
      answer.times(0);
    }

    /*
    FloatMatrix num_d2I7dFdD,num_d2I10dFdD;
    this->compute_d2I7dFdD_num(num_d2I7dFdD, F, D);
    this->compute_d2I10dFdD_num(num_d2I10dFdD, F, D);

    FloatMatrix num_d2I4dF2,num_d2I7dF2,num_d2I10dF2;
    this->compute_d2I4dF2_num(num_d2I4dF2, F);
    this->compute_d2I7dF2_num(num_d2I7dF2, F, D);
    this->compute_d2I10dF2_num(num_d2I10dF2, F, D);
    
    FloatMatrix dPdF, dEdF;
    this->compute_dPdF_dEdF(dPdF, dEdF, F, D, gp, tStep);

    FloatMatrix dPdD, dEdD;
    this->compute_dPdD_dEdD(dPdD, dEdD, F, D, gp, tStep);

    FloatMatrix test(num_d2I7dFdD);
    test.times(1. / J / this->epsilon_m);
    //test.add(- 1. / J / J / this->epsilon_m, cF_dI7dD);
    */		     
    
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: give3dMaterialStiffnessMatrix_dEdD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  answer.resize(3,3);
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    FloatArray vF = status->giveTempFVector();
    FloatArray D = status->giveTempDVector();
    FloatMatrix F;
    F.beMatrixForm(vF);
    double I9 = this->compute_I9(D);
    double I10 = this->compute_I10(F,D);
    
    FloatArray dI9dD, dI10dD;
    FloatMatrix d2I6dD2, d2I7dD2, d2I9dD2,d2I10dD2, dI9dI9, dI10dI10;
    this->compute_d2I6dD2(d2I7dD2, D);
    this->compute_d2I7dD2(d2I7dD2, F, D);
    this->compute_dI9dD(dI9dD, D);
    // this->compute_d2I9dD2(d2I9dD2, D);
    this->compute_dI10dD(dI10dD, F, D);
    //this->compute_d2I10dD2(d2I10dD2, F, D);
    dI9dI9.beDyadicProductOf(dI9dD, dI9dD);
    dI10dI10.beDyadicProductOf(dI10dD, dI10dD);
    double J = F.giveDeterminant();
    //
    answer.add( 0.5 * iEps2, d2I6dD2);
    //
    answer.add( 0.5 * iEps1 / J, d2I7dD2);
    //
    answer.add( iEps3 * I9, d2I9dD2);
    answer.add( iEps3, dI9dI9);
    //
    answer.add( iEps5 * I10, d2I10dD2);
    answer.add( iEps5, dI10dI10);
    //
    ////mix of invariants
    double FF = vF.dotProduct(vF);
    answer.add(iEps4 * FF, d2I6dD2);

      /*
    FloatMatrix num_d2I7dD2;
    this->compute_d2I7dD2_num(num_d2I7dD2, F, D);

    FloatMatrix dPdD, dEdD;
    this->compute_dPdD_dEdD(dPdD, dEdD, F, D, gp, tStep);
    FloatArray num_dI10dD;
    this->compute_dI10dD_num(num_dI10dD, F, D);
    */

    

}



double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I4(const FloatMatrix &F)
{
  
  FloatArray CN;
  FloatMatrix C;
  C.beTProductOf(F,F);
  CN.beProductOf(C,this->N);
  return this->N.dotProduct(CN);  
}


double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I5(const FloatMatrix &F)
{

  FloatMatrix H;
  hyperelasticMaterial->compute_2order_tensor_cross_product( H, F, F );

  FloatArray CN;
  FloatMatrix C;
  C.beTProductOf(H,H);
  CN.beProductOf(C,this->N);
  return this->N.dotProduct(CN);  
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I6(const FloatArray &D)
{
  return D.dotProduct(D);  
}


double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I7(const FloatMatrix &F, const FloatArray &D)
{
  
  FloatArray CD;
  FloatMatrix C;
  C.beTProductOf(F,F);
  CD.beProductOf(C,D);
  return D.dotProduct(CD);  
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I9(const FloatArray &D)
{
  return this->N.dotProduct(D);  
}


  
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I10(const FloatMatrix &F, const FloatArray &D)
{
  
  FloatArray CD;
  FloatMatrix C;
  C.beTProductOf(F,F);
  CD.beProductOf(C,D);
  return this->N.dotProduct(CD);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI1dF(FloatArray &answer, const FloatMatrix &F)
{

  double I1 = hyperelasticMaterial->compute_I1_C_from_F(F);
  double J_23 = pow(F.giveDeterminant(), -2./3.);
  FloatMatrix invF, invFt;
  invF.beInverseOf(F);
  invFt.beTranspositionOf(invF);
  invFt.times(2./3.*I1);

  FloatMatrix mAnswer;
  mAnswer = F;
  mAnswer.times(2.);
  mAnswer.subtract(invFt);
  mAnswer.times(J_23);
  answer.beVectorForm(mAnswer);
  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI2dF(FloatArray &answer, const FloatMatrix &F)
{
  answer.resize(9);
  FloatArray vF, vH;
  vF.beVectorForm(F);
  hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
  vH.times(0.5);
  double J = 1./3. * vH.dotProduct( vF );

  double HH = vH.dotProduct( vH );

  FloatArray vFxH;
  hyperelasticMaterial->compute_2order_tensor_cross_product( vFxH, vF, vH );

  answer.add( 3. / J / J * sqrt(HH), vFxH );
  answer.add( - 2. / J / J / J * HH * sqrt(HH), vH ); 
  
}






void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI4dF(FloatArray &answer, const FloatMatrix &F)
{
  answer.resize(9);
  
  FloatArray FN;
  FN.beProductOf(F,N);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {  
      answer.at(sm::giveVI(i,j)) = 2. * FN.at(i) * N.at(j);
    }
  }
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI5dF(FloatArray &answer, const FloatMatrix &F)
{
  answer.resize(9);
  FloatMatrix H;
  hyperelasticMaterial->compute_2order_tensor_cross_product( H, F, F );
  H.times(0.5);
  FloatArray HN;
  FloatArray HNN(9);
  HN.beProductOf(H,N);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {  
      HNN.at(sm::giveVI(i,j)) = 2. * HN.at(i) * N.at(j);
    }
  }

  FloatArray vF;
  vF.beVectorForm(F);
  hyperelasticMaterial->compute_2order_tensor_cross_product( answer, HNN, vF );
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI7dF(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(9);
  
  FloatArray FD;
  FD.beProductOf(F,D);
  
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {  
      answer.at(sm::giveVI(i,j)) = 2. * FD.at(i) * D.at(j);
    }
  }
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI8dF(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(9);
  //missing
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI10dF(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(9);
  
  FloatArray FN, FD;
  FN.beProductOf(F,N); 
  FD.beProductOf(F,D);
  
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {  
      answer.at(sm::giveVI(i,j)) = FN.at(i) * D.at(j)  + FD.at(i) * N.at(j);
    }
  }
  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI6dD(FloatArray &answer, const FloatArray &D)
{

   answer.resize(3);
   answer.add(2., D);

}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI7dD(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(3);

  FloatMatrix C;
  C.beTProductOf(F,F);
  
  answer.beProductOf(C,D);
  answer.times(2.);
  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI8dD(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(3);
  FloatMatrix H;
  hyperelasticMaterial->compute_2order_tensor_cross_product( H, F, F );
  H.times(0.5);
  FloatArray HD;
  HD.beProductOf(H,D);
  
  answer.beTProductOf(H,HD);
  answer.times(2.);
  
}






void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI9dD(FloatArray &answer, const FloatArray &D)
{
  double ND = N.dotProduct(D);
  answer.add(2.* ND, N);
}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI10dD(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(9);

  FloatMatrix C;
  C.beTProductOf(F,F);  
  answer.beProductOf(C,N);
}

///////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I1dF2(FloatMatrix &answer, const FloatMatrix &F)
{


    answer.resize(9, 9);
    answer.zero();

    double I1;
    FloatMatrix C, invF, I(3, 3);
    double J_23 = pow(F.giveDeterminant(), -2./3.);
    C.beTProductOf(F,F);
    invF.beInverseOf(F);
    I.beUnitMatrix();
    I1 = C.at(1,1) + C.at(2,2) + C.at(3,3); 

    
    for ( int i = 1; i <= 3; i++ ) {
      for ( int j = 1; j <= 3; j++ ) {
	for ( int k = 1; k <= 3; k++ ) {
	  for ( int l = 1; l <= 3; l++ ) {
	    answer.at( sm::giveVI(i, j), sm::giveVI(k, l) ) += 3. * I.at( i, k ) * I.at( j, l ) + I1 * invF.at( j, k )* invF.at( l, i ) - 2./3. * I1 * invF.at( j, i )* F.at( l, k )- 2. * invF.at( l, k )* F.at( i, j )- 2. * invF.at( j, i )* F.at( k, l );
	  }
	}
      }
    }  

    answer.times(2./3.*J_23);
  

}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I2dF2(FloatMatrix &answer, const FloatMatrix &F)
{
  answer.resize(9,9);
  FloatArray vF, vH;
  vF.beVectorForm(F);
  hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
  vH.times(0.5);
  double J = 1./3. * vH.dotProduct( vF );
  FloatMatrix H;
  H.beMatrixForm(vH);

  FloatArray vFH;
  FloatMatrix Fx, FxFx, HoH, FxHoH, FxHoHxF;
  FloatMatrix HoFH, FHoH;
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, F);
  hyperelasticMaterial->compute_2order_tensor_cross_product(vFH, vF, vH);
  hyperelasticMaterial->compute_dyadic_product(FHoH, vFH, vH);
  hyperelasticMaterial->compute_dyadic_product(HoFH, vH, vFH);
  hyperelasticMaterial->compute_dyadic_product(HoH, vH, vH);
  FxFx.beProductOf(Fx,Fx);
  hyperelasticMaterial->compute_dyadic_product(HoH, H, H);
  FxHoH.beProductOf(Fx, HoH);
  FxHoHxF.beProductTOf(HoH, Fx);
  /////
  double HH = vH.dotProduct( vH );
  answer.add(3. / J / J / sqrt(HH), FxHoHxF);
  answer.add(3. / J / J * sqrt(HH), FxFx);
  answer.add( - 6. / J / J / J * sqrt(HH), FHoH);
  answer.add( - 6. / J / J / J * sqrt(HH), HoFH);
  answer.add(   6. / J / J / J / J * HH * sqrt(HH), HoH);

}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I4dF2(FloatMatrix &answer, const FloatMatrix &F)
{
  answer.resize(9,9);
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {	  
	  answer.at(sm::giveVI(i,j), sm::giveVI(k,l)) = 2. * delta.at(i,k) * N.at(j) * N.at(l);
	}
      }
    }
  }   

}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I5dF2(FloatMatrix &answer, const FloatMatrix &F)
{
  answer.resize(9,9);
  FloatMatrix H;
  hyperelasticMaterial->compute_2order_tensor_cross_product( H, F, F );
  H.times(0.5);
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  FloatArray HN, HNN(9);
  FloatMatrix delta_NN(9,9);
  HN.beProductOf(H,N);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {  
      HNN.at(sm::giveVI(i,j)) = 2. * HN.at(i) * N.at(j);
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {
	  delta_NN.at(sm::giveVI(i,j), sm::giveVI(k,l)) = delta.at(i,k) * this->N.at(j) * this->N.at(l);
	}
      }
    }
  }

  FloatMatrix Fx, HNNx;
  FloatMatrix Fx_delta_NN, Fx_delta_NN_xF;

  
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, F);
  hyperelasticMaterial->compute_tensor_cross_product_tensor(HNNx, HNN);


  Fx_delta_NN.beProductOf(Fx, delta_NN);
  Fx_delta_NN_xF.beProductTOf(delta_NN, Fx);

  answer.add(Fx_delta_NN_xF);
  answer.add(HNNx);

}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  answer.resize(9,9);
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {  
	  answer.at(sm::giveVI(i,j), sm::giveVI(k,l)) = 2. * delta.at(i,k) * D.at(j) * D.at(l);
	}
      }
    }
  }
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  answer.resize(9,9);
  FloatMatrix H;
  hyperelasticMaterial->compute_2order_tensor_cross_product( H, F, F );
  H.times(0.5);

  FloatMatrix delta(3,3);
  delta.beUnitMatrix();

  FloatArray HD,  HDD(9);
  FloatMatrix delta_DD(9,9);
  HD.beProductOf(H,D);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {  
      HDD.at(sm::giveVI(i,j)) = 2. * HD.at(i) * D.at(j);
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {
	  delta_DD.at(sm::giveVI(i,j), sm::giveVI(k,l)) = delta.at(i,k) * D.at(j) * D.at(l);
	}
      }
    }
  }

  FloatMatrix Fx, HDDx;
  FloatMatrix Fx_delta_DD, Fx_delta_DD_xF;

  
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, F);
  hyperelasticMaterial->compute_tensor_cross_product_tensor(HDDx, HDD);

  Fx_delta_DD.beProductOf(Fx, delta_DD);
  Fx_delta_DD_xF.beProductTOf(delta_DD, Fx);

  answer.add(Fx_delta_DD_xF);
  answer.add(HDDx);

}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I10dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{


  answer.resize(9,9);
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {  
	  answer.at(sm::giveVI(i,j), sm::giveVI(k,l)) = delta.at(i,k) * D.at(l) * N.at(j) +  N.at(l) * D.at(j) * delta.at(i,k);
	}
      }
    }
  }

}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I6dD2(FloatMatrix &answer, const FloatArray &D)
{
  answer.resize(3,3);
  answer.beUnitMatrix();
  answer.times(2.);
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7dD2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  answer.beTProductOf(F,F);
  answer.times(2);
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8dD2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  FloatMatrix H;
  hyperelasticMaterial->compute_2order_tensor_cross_product( H, F, F );
  H.times(0.5);
  answer.beTProductOf(H,H);
  answer.times(2);
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7dFdD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(9,3);
  FloatArray FD;
  FD.beProductOf(F,D);
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();

  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	answer.at(sm::giveVI(i,j), k) = 2. * (F.at(i,k)*D.at(j) + FD.at(i)*delta.at(j,k));
      }
    }
  }  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I10dFdD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(9,3);
  FloatArray FN;
  FN.beProductOf(F,N);
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();

  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	answer.at(sm::giveVI(i,j), k) = FN.at(i) * delta.at(j,k) + F.at(i,k) * N.at(j) ;
      }
    }
  }  
}


int
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: computeAcousticTensorMinEigenvalue(GaussPoint *gp, TimeStep *tStep)
{
   FloatMatrix Duu, Ddd, Dud;
   this->give3dMaterialStiffnessMatrix_dPdF(Duu, TangentStiffness, gp, tStep);
   this->give3dMaterialStiffnessMatrix_dEdD(Ddd, TangentStiffness, gp, tStep);
   this->give3dMaterialStiffnessMatrix_dPdD(Dud, TangentStiffness, gp, tStep);
   double minE;
   FloatArray n(3);
   FloatMatrix nn, I(3,3);
   int nStepsI, nStepsJ;
   nStepsI = nStepsJ = 20;
   int index = 0;
   for(int iStep = 1; iStep <= nStepsI; iStep ++) {
     double theta = (iStep - 1.)/(nStepsI-1.) * 2. * 3.14159265358979323; 
     for(int jStep = 1; jStep <= nStepsJ; jStep++) {
       double phi = (jStep - 1.)/(nStepsJ-1.) * 2. * 3.14159265358979323; 
       n.at(1) = sin(theta) * cos(phi);
       n.at(2) = sin(theta) * sin(phi);
       n.at(3) = cos(theta);
       index++;
       I.beUnitMatrix();
       nn.beProductTOf(n,n);
       I.subtract(nn);
       FloatMatrix a,b, iDdd;
       a.beProductOf(Ddd, I);
       b.beProductOf(I, a);
       this->computePseudoInverse(iDdd, b);
       FloatMatrix G(3,3), Q(3,3);
       for(int i = 1; i <= 3; i ++) {
	 for(int j = 1; j <= 3; j ++) {
	   for(int k = 1; k <= 3; k ++) {
	     G.at(j,k) += Dud.at(sm::giveVI(i,j),k) * n.at(i);
	     for(int l = 1; l <= 3; l ++) {
	       Q.at(j, l) += Duu.at(sm::giveVI(i,j), sm::giveVI(k,l)) * n.at(i) * n.at(k);
	     }
	   }
	 }
       }
       FloatMatrix iDdd_G, G_iDdd_G, At;
       iDdd_G.beProductTOf(iDdd, G);
       G_iDdd_G.beProductOf(G, iDdd_G);
       At = Q;
       At.subtract(G_iDdd_G);
       FloatArray eval;
       FloatMatrix v;
       At.jaco_(eval, v, 1.e-10);
       if(index == 1) {
	 minE = eval.at(eval.giveIndexMinElem());
       } else {
	 minE = min(minE, eval.at(eval.giveIndexMinElem()));
       }
       
     }
   }
   
   return minE;
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: computePseudoInverse(FloatMatrix &iDdd, FloatMatrix &b)
{
  FloatArray eval;
  FloatMatrix v;
  iDdd.resize(3,3);
  b.jaco_(eval, v, 1.e-8);
  int min = eval.giveIndexMinElem();
  for ( int i = 1; i <= 3; i++ ) {
    for ( int j = 1; j <= 3; j++ ) {
      iDdd.at(i, j) = 1./eval.at(1) * v.at(i, 1) * v.at(j, 1) + 1./eval.at(2) * v.at(i, 2) * v.at(j, 2);  
    }
  }
   
}

int
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

  ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
  
  if ( type == IST_ElectricDisplacementVector) {
      answer.resize(3);
      answer = status->giveDVector();
      return 1;
  } else if ( type == IST_ElectricFieldVector ) {
      answer.resize(3);
      answer = status->giveEVector();
      return 1;
  } else if ( type == IST_AcousticTensorMinEigenvalue ) {
    answer.resize(1);
    answer.at(1) = this->computeAcousticTensorMinEigenvalue(gp, tStep);
  } else {
      return  this->hyperelasticMaterial->giveIPValue(answer, gp, type, tStep);
  }
   
}
    



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI4dF_num(FloatArray &answer, const FloatMatrix &F)
{

  FloatMatrix Fp;
  auto I = this->compute_I4(F);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp = F;
      Fp.at(i,j) += pert;
      auto Ip = this->compute_I4(Fp);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 ::  compute_dI7dF_num(FloatArray &answer,const FloatMatrix &F, const FloatArray &D)
{

  FloatMatrix Fp;
  auto I = this->compute_I7(F,D);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp = F;
      Fp.at(i,j) += pert;
      auto Ip = this->compute_I7(Fp,D);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI10dF_num(FloatArray &answer,const FloatMatrix &F, const FloatArray &D)
{

  FloatMatrix Fp;
  auto I = this->compute_I10(F,D);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp = F;
      Fp.at(i,j) += pert;
      auto Ip = this->compute_I10(Fp,D);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI7dD_num(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_I7(F,D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_I7(F,Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI10dD_num(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_I10(F,D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_I10(F,Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I4dF2_num(FloatMatrix &answer,const FloatMatrix &F)
{


  FloatArray dI, dIp;
  this->compute_dI4dF(dI, F);

  FloatMatrix Fp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp = F;
      Fp.at(i,j) += pert;
      this->compute_dI4dF(dIp, Fp);
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {
	  answer.at(sm::giveVI(i,j), sm::giveVI(k,l)) = dIp.at(sm::giveVI(k,l)) -  dI.at(sm::giveVI(k,l));
	}
      }
    }
  }

  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7dF2_num(FloatMatrix &answer,const FloatMatrix &F, const FloatArray&D)
{


  FloatArray dI, dIp;
  this->compute_dI7dF(dI, F, D);

  FloatMatrix Fp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp = F;
      Fp.at(i,j) += pert;
      this->compute_dI7dF(dIp, Fp, D);
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {
	  answer.at(sm::giveVI(i,j), sm::giveVI(k,l)) = dIp.at(sm::giveVI(k,l)) -  dI.at(sm::giveVI(k,l));
	}
      }
    }
  }

  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I10dF2_num(FloatMatrix &answer,const FloatMatrix &F, const FloatArray&D)
{


  FloatArray dI, dIp;
  this->compute_dI10dF(dI, F, D);

  FloatMatrix Fp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp = F;
      Fp.at(i,j) += pert;
      this->compute_dI10dF(dIp, Fp, D);
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {
	  answer.at(sm::giveVI(k,l), sm::giveVI(i,j)) = dIp.at(sm::giveVI(k,l)) -  dI.at(sm::giveVI(k,l));
	}
      }
    }
  }

  answer.times(1./pert);  
}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 ::  compute_d2I7dFdD_num(FloatMatrix &answer, const FloatMatrix &F, const FloatArray&D)
{


  FloatArray dI, dIp;
  this->compute_dI7dF(dI, F,D);
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dI7dF(dIp, F,Dp);
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {
	  answer.at(sm::giveVI(k,l),i) = dIp.at(sm::giveVI(k,l)) -  dI.at(sm::giveVI(k,l));
	}
      }
  }
  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 ::  compute_d2I10dFdD_num(FloatMatrix &answer, const FloatMatrix &F, const FloatArray&D)
{


  FloatArray dI, dIp;
  this->compute_dI10dF(dI, F,D);
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dI10dF(dIp, F,Dp);
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {
	  answer.at(sm::giveVI(k,l),i) = dIp.at(sm::giveVI(k,l)) -  dI.at(sm::giveVI(k,l));
	}
      }
  }
  answer.times(1./pert);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7dD2_num(FloatMatrix &answer, const FloatMatrix &F, const FloatArray&D)
{


  FloatArray dI, dIp;
  this->compute_dI7dD(dI, F,D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(3,3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dI7dD(dIp, F,Dp);
      for(int k = 1; k <= 3; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);  
}









void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dPdF_dEdF(FloatMatrix &dPdF,FloatMatrix &dEdF, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep)
{


  FloatArray vP, vPp, vF, vFp, Ep, E;
  vF.beVectorForm(F);
  this->give_FirstPKStressVector_ElectricalFieldVector_3d(vP, E, gp, vF, D, tStep);

  FloatMatrix Fp;
  auto pert = 1.e-8;
  dPdF.resize(9,9);
  dEdF.resize(3,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp = F;
      Fp.at(i,j) += pert;
      vFp.beVectorForm(Fp);
      this->give_FirstPKStressVector_ElectricalFieldVector_3d(vPp, Ep, gp, vFp, D, tStep);

      for(int k = 1; k <= 3; k++) {
	dEdF.at(k, sm::giveVI(i,j)) = Ep.at(k) -  E.at(k);
	for(int l = 1; l <= 3; l++) {
	  dPdF.at(sm::giveVI(k,l), sm::giveVI(i,j)) = vPp.at(sm::giveVI(k,l)) -  vP.at(sm::giveVI(k,l));
	}
      }
    }
  }
  dPdF.times(1./pert);
  dEdF.times(1./pert);
  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dPdF_dEdF_fake(FloatMatrix &dPdF,FloatMatrix &dEdF, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep)
{

  //double em = this->epsilon_m;
  /*  this->epsilon_m = 0;
  this->epsilon_f = 5.e-12;
  */
  FloatArray vP, vPp, vF, vFp, Ep, E;
  vF.beVectorForm(F);
  this->give_FirstPKStressVector_ElectricalFieldVector_3d(vP, E, gp, vF, D, tStep);

  FloatMatrix Fp;
  auto pert = 1.e-8;
  dPdF.resize(9,9);
  dEdF.resize(3,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp = F;
      Fp.at(i,j) += pert;
      vFp.beVectorForm(Fp);
      this->give_FirstPKStressVector_ElectricalFieldVector_3d(vPp, Ep, gp, vFp, D, tStep);

      for(int k = 1; k <= 3; k++) {
	dEdF.at(k, sm::giveVI(i,j)) = Ep.at(k) -  E.at(k);
	for(int l = 1; l <= 3; l++) {
	  dPdF.at(sm::giveVI(k,l), sm::giveVI(i,j)) = vPp.at(sm::giveVI(k,l)) -  vP.at(sm::giveVI(k,l));
	}
      }
    }
  }
  dPdF.times(1./pert);
  dEdF.times(1./pert);

  /*  this->epsilon_m = em;
  this->epsilon_f = 0;
  */
  
}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dPdD_dEdD(FloatMatrix &dPdD,FloatMatrix &dEdD, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep)
{


  FloatArray vP, vPp, vF, Dp, Ep, E;
  vF.beVectorForm(F);
  this->give_FirstPKStressVector_ElectricalFieldVector_3d(vP, E, gp, vF, D, tStep);

  FloatMatrix Fp;
  auto pert = 1.e-15;
  dPdD.resize(9,3);
  dEdD.resize(3,3);
  for(int i = 1; i <= 3; i++) {
    Dp = D;
    Dp.at(i) += pert;
    this->give_FirstPKStressVector_ElectricalFieldVector_3d(vPp, Ep, gp, vF, Dp, tStep);
    for(int k = 1; k <= 3; k++) {
      dEdD.at(k, i) = Ep.at(k) -  E.at(k);
      for(int l = 1; l <= 3; l++) {
	dPdD.at(sm::giveVI(k,l), i) = vPp.at(sm::giveVI(k,l)) -  vP.at(sm::giveVI(k,l));
      }
    }
  }
  dPdD.times(1./pert);
  dEdD.times(1./pert);
  
}





} // end namespace oofem
