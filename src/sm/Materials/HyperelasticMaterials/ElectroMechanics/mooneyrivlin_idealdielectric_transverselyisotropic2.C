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
    IR_GIVE_OPTIONAL_FIELD(ir, this->alpha, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_alpha);
    IR_GIVE_OPTIONAL_FIELD(ir, this->beta, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_beta);

    
    double epsilon1 = 0, epsilon2 = 0, epsilon3 = 0, epsilon4 = 0, epsilon5 = 0, epsilon6 = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon1, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon1);
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon2, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon2);
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon3, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon3);
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon4, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon4);
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon5, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon5);
    IR_GIVE_OPTIONAL_FIELD(ir, epsilon6, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon6);
    IR_GIVE_OPTIONAL_FIELD(ir, a, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_a);    
    
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

    if(epsilon6 == 0) {
      this->iEps6 = 0.;
    } else {
      this->iEps6 = 1. / epsilon6;
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
    this->hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    //
    FloatMatrix H;
    H.beMatrixForm(vH);    
    ////////////////    
    // First Piol-Kirchhoff stress
    double I4 = this->compute_I4(J, vH, vF);
    double I5 = this->compute_I5(J, vH, vF);
    double I7 = this->compute_I7(J, vH, vF, D);
    ////////////////////////////////////////////
    FloatArray dI1devdF, dI2devpoldF;
    this->compute_dI1dev_dF(dI1devdF, J, vH, vF);
    this->compute_dI2devpol_dF(dI2devpoldF, J, vH, vF);
    /////////////////////////////////////////////////
    FloatArray dI4dF, dI5dF;
    this->compute_dI4_dF(dI4dF, J, vH, vF);
    this->compute_dI5_dF(dI5dF, J, vH, vF);
    /////////////////////////////////////////////////
    FloatArray dI7dF,dI7dD;
    this->compute_dI7_dF_dD(dI7dF, dI7dD, J, vH, vF, D);
    FloatArray dK2DdF,dK2DdD, dK2DpoldF, dK2DpoldD;
    this->compute_dK2Dinf_dF_dD(dK2DdF,dK2DdD, J, vH, vF, D);
    this->compute_dK2Dinfpol_dF_dD(dK2DpoldF,dK2DpoldD,J, vH, vF, D);
    ///////////////////////////////////
    FloatArray dI6dD, dK1DdD;
    this->compute_dI6_dD(dI6dD, D);
    this->compute_dK1Dinf_dD(dK1DdD, D);
    ///////////////////////////////////
    FloatArray dK2CdF, dK2CdD;
    double K2C = this->compute_K2_Cinf(J, vH, vF, D);
    this->compute_dK2Cinf_dF_dD(dK2CdF, dK2CdD, J, vH, vF, D);
    FloatArray dK1CdD;
    this->compute_dK1Cinf_dD(dK1CdD, D);
    //////////////////////////////////////
    // isotropic mechanical invariants
    vP.add(0.5 * mu1, dI1devdF);
    vP.add(0.5 * mu2, dI2devpoldF);
    vP.add( lambda * (J - 1.), vH);
    ///////////////////////////////////
    // anisotropic mechanical invariants
    vP.add(0.5 * mu3 * pow(I4, alpha - 1.), dI4dF);
    vP.add(0.5 * mu3 * pow(I5, beta  - 1.), dI5dF);
    vP.add( - mu3 / J, vH);
    ///////////////////////////////////
    // electro-mechanical invariants
    // I7/J
    vP.add(   0.5 * this->iEps1 / J, dI7dF);
    vP.add( - 0.5 * this->iEps1 / J / J * I7, vH);
    //
    vP.add(   0.5 * this->iEps4 , dK2DpoldF);
    vP.add(   0.5 * this->iEps5 , dK2DdF);
    ////////////////////////////////////
    //// contribution of K2C
    ///////////////////////////////////
    vP.add(this->a, dK2CdF);
    ///////////////////////////////////
     // Electric field
    E.add( 0.5 * this->iEps2, dI6dD);    
    E.add( 0.5 * this->iEps1 / J , dI7dD);
    E.add( 0.5 * this->iEps3, dK1DdD);
    E.add( 0.5 * this->iEps4, dK2DpoldD);   
    E.add( 0.5 * this->iEps5, dK2DdD);
    ////////////////////////////////////
    //// contribution of K2C^2
    E.add(this->a, dK2CdD);
    E.add(-this->a, dK1CdD);
    ////////////////////////////////////
    
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
    answer.resize(9,9);
    //
    double I4 = this->compute_I4(J, vH, vF);
    double I5 = this->compute_I5(J, vH, vF);
    double I7 = this->compute_I7(J, vH, vF, D);
    //
    FloatArray dI4dF, dI5dF, dI7dF, dI7dD;
    FloatMatrix d2I1devdF2, d2I2devpoldF2,d2I3dF2, d2I4dF2, d2I5dF2, d2I7dF2, d2K2DdF2, d2K2DpoldF2;
    HH.beDyadicProductOf(vH,vH);
    this->compute_d2I1dev_dF2(d2I1devdF2, J, vH, vF);
    this->compute_d2I2devpol_dF2(d2I2devpoldF2, J, vH, vF);
    this->compute_d2I3_dF2(d2I3dF2, J, vH, vF);
    //   
    this->compute_dI4_dF(dI4dF, J, vH, vF);
    this->compute_d2I4_dF2(d2I4dF2, J, vH, vF);
    //
    this->compute_dI5_dF(dI5dF, J, vH, vF);
    this->compute_d2I5_dF2(d2I5dF2, J, vH, vF);
    //
    this->compute_dI7_dF_dD(dI7dF, dI7dD, J, vH, vF, D);
    this->compute_d2I7_dF2(d2I7dF2, J, vH, vF, D);
    //
    this->compute_d2K2Dinfpol_dF2(d2K2DpoldF2, J, vH, vF, D);
    this->compute_d2K2Dinf_dF2(d2K2DdF2, J, vH, vF, D);
    // Isotropic stiffness - mechanical only
    answer.add(0.5 * mu1, d2I1devdF2);
    answer.add(0.5 * mu2, d2I2devpoldF2);
    answer.add( lambda * ( J - 1 ) , d2I3dF2);
    answer.add( lambda , HH);
    // Anisotropic stiffness - mechanical only
    FloatMatrix dI4dI4, dI5dI5;
    dI4dI4.beDyadicProductOf(dI4dF, dI4dF);
    dI5dI5.beDyadicProductOf(dI5dF, dI5dF);
    answer.add( 0.5 * mu3 * pow(I4, alpha-1), d2I4dF2);
    answer.add( 0.5 * mu3 * (alpha - 1.) * pow(I4, alpha-2), dI4dI4);
    answer.add( 0.5 * mu3 * pow(I5, beta-1.), d2I5dF2);
    answer.add( 0.5 * mu3 * (beta - 1) * pow(I5, beta-2), dI5dI5);
    answer.add(       mu3 / J / J, HH);
    answer.add( -     mu3  / J, d2I3dF2);
    // Isotropic stiffness - electromechanical
    //2nd derivative of I7/J
    FloatMatrix HdI7, dI7H;
    dI7H.beDyadicProductOf(dI7dF, vH);
    HdI7.beDyadicProductOf(vH, dI7dF);
    answer.add(   0.5 * iEps1 / J, d2I7dF2);
    answer.add(         iEps1 * I7 / J / J/ J,  HH);
    answer.add( - 0.5 * iEps1 / J / J, dI7H);
    answer.add( - 0.5 * iEps1 / J / J, HdI7);
    answer.add( - 0.5 * iEps1 * I7 / J / J, d2I3dF2);
    ///////////////////////////////////////////////////
    FloatArray dI7dF_num;
    FloatMatrix d2I7dF2_num;
    this->compute_dI7dF_num(dI7dF_num, J, vH, vF, D);
    this->compute_d2I7dF2_num(d2I7dF2_num,J, vH, vF, D);
    ///////////////////////////////////////////////////
    // Anisotropic stiffness - electromechanical
    answer.add( 0.5 * iEps4, d2K2DpoldF2);
    // nonpolyconvex invariant
    answer.add( 0.5 * iEps5, d2K2DdF2);
    ////////////////////////////////////
    //// contribution of K2C
    ///////////////////////////////////
    FloatMatrix d2K2CdF2;
    this->compute_d2K2Cinf_dF2(d2K2CdF2, J, vH, vF, D);
    //
    answer.add( this->a, d2K2CdF2);
    /*    
    FloatMatrix d2K2DpoldF2num;
    this->compute_d2K2DpoldF2_num(d2K2DpoldF2num,J, vH, vF, D);
    FloatMatrix di(d2K2DpoldF2num);
    di.subtract(d2K2DpoldF2);


    FloatMatrix dPdF, dEdF;
    this->compute_dPdF_dEdF(dPdF, dEdF, F, D, gp, tStep);

    FloatMatrix diff(answer);
    diff.subtract(dPdF);
    double n = diff.computeFrobeniusNorm();
    if(n > 10000.0) {
      int test = 1;
    }


    */

    
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: give3dMaterialStiffnessMatrix_dPdD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(9,3);
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    FloatArray vF = status->giveTempFVector();
    FloatArray D = status->giveTempDVector();
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    ////////////////////////////////////
    FloatArray dI7dF, dI7dD;
    this->compute_dI7_dF_dD(dI7dF, dI7dD, J, vH, vF, D);
    //
    FloatMatrix d2I7dFdD, H_dI7;
    H_dI7.beDyadicProductOf(vH, dI7dD);
    //
    this->compute_d2I7_dFdD(d2I7dFdD, J, vH, vF, D);
    //
    answer.add(  0.5 * this->iEps1 / J, d2I7dFdD);
    answer.add(- 0.5 * this->iEps1 / J / J, H_dI7);
    /////////////////////////////////////////////////
    FloatArray dI7dD_num;
    FloatMatrix d2I7dFdD_num;
    this->compute_dI7dD_num(dI7dD_num, J, vH, vF, D);
    this->compute_d2I7dFdD_num(d2I7dFdD_num, J, vH, vF, D);
    ///////////////////////////////////////////////// 
    FloatMatrix d2K2DdFdD, d2K2DpoldFdD;
    this->compute_d2K2Dinf_dFdD(d2K2DdFdD, J, vH, vF, D);
    this->compute_d2K2Dinfpol_dFdD(d2K2DpoldFdD,J, vH, vF, D);
    //
    answer.add(  0.5 * this->iEps4, d2K2DpoldFdD);
    answer.add(  0.5 * this->iEps5, d2K2DdFdD);
    ////////////////////////////////////
    //// contribution of K2C
    ///////////////////////////////////
    FloatMatrix d2K2CdFdD;
    this->compute_d2K2Cinf_dFdD(d2K2CdFdD, J, vH, vF, D);
    //
    answer.add( this->a, d2K2CdFdD);   
    if(mode == SecantStiffness) {
      answer.times(0);
    }

    /*
    FloatMatrix dPdD, dEdD, dPdF, dEdF;
    FloatMatrix F;
    F.beMatrixForm(vF);
    this->compute_dPdD_dEdD(dPdD, dEdD, F, D, gp, tStep);
    this->compute_dPdF_dEdF(dPdF, dEdF, F, D, gp, tStep);
    FloatMatrix diff(answer);
    diff.subtract(dPdD);
    double n = diff.computeFrobeniusNorm();
    FloatMatrix d2K2DpoldFdD_num;
    this->compute_d2K2DpoldFdD_num(d2K2DpoldFdD_num,J, vH, vF, D);
    if(n > 10000.0) {
      FloatArray dI7dF_num;
      FloatMatrix d2I7dFdD_num, d2K2DpoldFdD_num;
      this->compute_dI7dF_num(dI7dF_num, J, vH, vF, D);
      this->compute_d2I7dFdD_num(d2I7dFdD_num, J, vH, vF, D);
      this->compute_d2K2Dinfpol_dFdD(d2K2DpoldFdD_num,J, vH, vF, D);

      int test = 1;
    }
    */
    
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: give3dMaterialStiffnessMatrix_dEdD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    answer.resize(3,3);
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    FloatArray vF = status->giveTempFVector();
    FloatArray D = status->giveTempDVector();
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    //
    FloatMatrix d2I6dD2, d2I7dD2, d2K1DdD2, d2K2DpoldD2, d2K2DdD2;
    //
    this->compute_d2I6_dD2(d2I6dD2, D);
    this->compute_d2I7_dD2(d2I7dD2, J, vH, vF, D);
    ////////////////////////////////////////////////
    FloatMatrix d2I7dD2_num;
    this->compute_d2I7dD2_num(d2I7dD2_num,J, vH, vF, D);
    ///////////////////////////////////////////////////
    this->compute_d2K1Dinf_dD2(d2K1DdD2, D);
    this->compute_d2K2Dinfpol_dD2(d2K2DpoldD2, J, vH, vF, D);
    this->compute_d2K2Dinf_dD2(d2K2DdD2, J, vH, vF, D);
    //
    answer.add( 0.5 * iEps2, d2I6dD2);
    //
    answer.add( 0.5 * iEps1 / J, d2I7dD2);
    //
    answer.add( 0.5 * iEps3, d2K1DdD2);
    //
    answer.add( 0.5 * iEps4, d2K2DpoldD2);
    //
    answer.add( 0.5 * iEps5, d2K2DdD2);

    /*
    FloatMatrix dPdD, dEdD;
    FloatMatrix F;
    F.beMatrixForm(vF);
    this->compute_dPdD_dEdD(dPdD, dEdD, F, D, gp, tStep);
    
    FloatMatrix d2K2DpoldD2_num;
    this->compute_d2K2DpoldD2_num(d2K2DpoldD2_num,J, vH, vF, D);
    
    FloatMatrix diff(answer);
    diff.subtract(dEdD);
    double n = diff.computeFrobeniusNorm();
    if(n > 10000) {
      FloatMatrix d2I6dD2_num, d2I7dD2_num, d2K1DdD2_num, d2K2DpoldD2_num;
      this->compute_d2I6dD2_num(d2I6dD2_num, D);
      this->compute_d2I7dD2_num(d2I7dD2_num, J, vH, vF, D);
      this->compute_d2K1DdD2_num(d2K1DdD2_num, D);
      this->compute_d2K2DpoldD2_num(d2K2DpoldD2_num,J, vH, vF, D);
      int test = 1;
    }
    */
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: give3dMaterialStiffnessMatrix_dPdF_from(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep, const FloatArray &vF, const FloatArray &D)
{
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    // First Piola-Kirchhoff stress
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
    FloatMatrix H, HH;
    H.beMatrixForm(vH);
    ////////////////////////////////////    
    answer.resize(9,9);
    //
    double I4 = this->compute_I4(J, vH, vF);
    double I5 = this->compute_I5(J, vH, vF);
    double I7 = this->compute_I7(J, vH, vF, D);
    //
    FloatArray dI4dF, dI5dF, dI7dF, dI7dD;
    FloatMatrix d2I1devdF2, d2I2devpoldF2,d2I3dF2, d2I4dF2, d2I5dF2, d2I7dF2, d2K2DdF2, d2K2DpoldF2;
    HH.beDyadicProductOf(vH,vH);
    this->compute_d2I1dev_dF2(d2I1devdF2, J, vH, vF);
    this->compute_d2I2devpol_dF2(d2I2devpoldF2, J, vH, vF);
    this->compute_d2I3_dF2(d2I3dF2, J, vH, vF);
    //   
    this->compute_dI4_dF(dI4dF, J, vH, vF);
    this->compute_d2I4_dF2(d2I4dF2, J, vH, vF);
    //
    this->compute_dI5_dF(dI5dF, J, vH, vF);
    this->compute_d2I5_dF2(d2I5dF2, J, vH, vF);
    //
    this->compute_dI7_dF_dD(dI7dF, dI7dD, J, vH, vF, D);
    this->compute_d2I7_dF2(d2I7dF2, J, vH, vF, D);
    //
    this->compute_d2K2Dinfpol_dF2(d2K2DpoldF2, J, vH, vF, D);
    this->compute_d2K2Dinf_dF2(d2K2DdF2, J, vH, vF, D);
    // Isotropic stiffness - mechanical only
    answer.add(0.5 * mu1, d2I1devdF2);
    answer.add(0.5 * mu2, d2I2devpoldF2);
    answer.add( lambda * ( J - 1 ) , d2I3dF2);
    answer.add( lambda , HH);
    // Anisotropic stiffness - mechanical only
    FloatMatrix dI4dI4, dI5dI5;
    dI4dI4.beDyadicProductOf(dI4dF, dI4dF);
    dI5dI5.beDyadicProductOf(dI5dF, dI5dF);
    answer.add( 0.5 * mu3 * pow(I4, alpha-1), d2I4dF2);
    answer.add( 0.5 * mu3 * (alpha - 1.) * pow(I4, alpha-2), dI4dI4);
    answer.add( 0.5 * mu3 * pow(I5, beta-1.), d2I5dF2);
    answer.add( 0.5 * mu3 * (beta - 1) * pow(I5, beta-2), dI5dI5);
    answer.add(       mu3 / J / J, HH);
    answer.add( -     mu3  / J, d2I3dF2);
    // Isotropic stiffness - electromechanical
    //2nd derivative of I7/J
    FloatMatrix HdI7, dI7H;
    dI7H.beDyadicProductOf(dI7dF, vH);
    HdI7.beDyadicProductOf(vH, dI7dF);
    answer.add(   0.5 * iEps1 / J, d2I7dF2);
    answer.add(         iEps1 * I7 / J / J/ J,  HH);
    answer.add( - 0.5 * iEps1 / J / J, dI7H);
    answer.add( - 0.5 * iEps1 / J / J, HdI7);
    answer.add( - 0.5 * iEps1 * I7 / J / J, d2I3dF2);
    // Anisotropic stiffness - electromechanical
    answer.add( 0.5 * iEps4, d2K2DpoldF2);
    // nonpolyconvex invariant
    answer.add( 0.5 * iEps5, d2K2DdF2);
    ////////////////////////////////////
    //// contribution of K2C^2/J
    ///////////////////////////////////
    double K2C = this->compute_K2_Cinf(J, vH, vF, D);
    FloatArray dK2CdF, dK2CdD;
    //
    this->compute_dK2Cinf_dF_dD(dK2CdF, dK2CdD, J, vH, vF, D);
    FloatMatrix d2K2CdF2, dK2CdK2C, dK2C_H, H_dK2C;
    //
    dK2CdK2C.beDyadicProductOf(dK2CdF,dK2CdF);
    dK2C_H.beDyadicProductOf(dK2CdF, vH);
    H_dK2C.beDyadicProductOf(vH,dK2CdF);
    //
    this->compute_d2K2Cinf_dF2(d2K2CdF2, J, vH, vF, D);
    //
    answer.add( iEps6 * K2C / J, d2K2CdF2);
    //
    answer.add( iEps6 *K2C / J / J, dK2C_H);
    //
    answer.add( iEps6 *K2C / J / J, H_dK2C);
    //
    answer.add( 0.5 * iEps6 * K2C * K2C / J / J / J, HH);
    //
    FloatMatrix Fx;
    hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, vF);
    answer.add( 0.5 * iEps6 * K2C * K2C / J / J , Fx);

    
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: give3dMaterialStiffnessMatrix_dPdD_from(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep, const FloatArray &vF, const FloatArray &D)
{
    answer.resize(9,3);
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    ////////////////////////////////////
    FloatArray dI7dF, dI7dD;
    this->compute_dI7_dF_dD(dI7dF, dI7dD, J, vH, vF, D);
    //
    FloatMatrix d2I7dFdD, H_dI7;
    H_dI7.beDyadicProductOf(vH, dI7dD);
    //
    this->compute_d2I7_dFdD(d2I7dFdD, J, vH, vF, D);
    //
    answer.add(  0.5 * this->iEps1 / J, d2I7dFdD);
    answer.add(- 0.5 * this->iEps1 / J / J, H_dI7);
    ///////////////////////////////////////////////// 
    FloatMatrix d2K2DdFdD, d2K2DpoldFdD;
    this->compute_d2K2Dinf_dFdD(d2K2DdFdD, J, vH, vF, D);
    this->compute_d2K2Dinfpol_dFdD(d2K2DpoldFdD,J, vH, vF, D);
    //
    answer.add(  0.5 * this->iEps4, d2K2DpoldFdD);
    answer.add(  0.5 * this->iEps5, d2K2DdFdD);

    ////////////////////////////////////
    //// contribution of K2C^2/J
    ///////////////////////////////////
    double K2C = this->compute_K2_Cinf(J, vH, vF, D);
    //
    FloatArray dK2CdF, dK2CdD;
    this->compute_dK2Cinf_dF_dD(dK2CdF, dK2CdD, J, vH, vF, D);
    //
    FloatMatrix dK2CdK2dD, dK2CdFdD,H_dK2dD;
    dK2CdK2dD.beDyadicProductOf(dK2CdF,dK2CdD);
    H_dK2dD.beDyadicProductOf(vH,dK2CdD);
    //
    this->compute_d2K2Cinf_dFdD(dK2CdFdD, J, vH, vF, D);
    //
    answer.add( this->iEps6 * K2C / J, dK2CdFdD);
    answer.add( this->iEps6 / J, dK2CdK2dD);
    answer.add( - this->iEps6 * K2C / J / J, H_dK2dD);

    
    if(mode == SecantStiffness) {
      answer.times(0);
    }
    
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: give3dMaterialStiffnessMatrix_dEdD_from(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep, const FloatArray &vF, const FloatArray &D)
{

    answer.resize(3,3);
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    //
    FloatMatrix d2I6dD2, d2I7dD2, d2K1DdD2, d2K2DpoldD2, d2K2DdD2;
    //
    this->compute_d2I6_dD2(d2I6dD2, D);
    this->compute_d2I7_dD2(d2I7dD2, J, vH, vF, D);
    this->compute_d2K1Dinf_dD2(d2K1DdD2, D);
    this->compute_d2K2Dinfpol_dD2(d2K2DpoldD2, J, vH, vF, D);
    this->compute_d2K2Dinf_dD2(d2K2DdD2, J, vH, vF, D);
    //
    answer.add( 0.5 * iEps2, d2I6dD2);
    //
    answer.add( 0.5 * iEps1 / J, d2I7dD2);
    //
    answer.add( 0.5 * iEps3, d2K1DdD2);
    //
    answer.add( 0.5 * iEps4, d2K2DpoldD2);
    //
    answer.add( 0.5 * iEps5, d2K2DdD2);


}




int
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: computeAcousticTensorMinEigenvalue(GaussPoint *gp, TimeStep *tStep)
{

  FloatMatrix Duu, Ddd, Dud, minQ;
  /*if(tStep->giveNumber() >= 1 && gp->giveElement()->giveNumber()==6601) {

    FloatArray vF, D;
    vF={1.258406, 1.020216, 0.934341, -0.003643, 0.064165, -0.011682, -0.032953, -0.735274, 0.482873};
    D={0.000933, 0.001340, -0.007123};  
    
    this->give3dMaterialStiffnessMatrix_dPdF_from(Duu, TangentStiffness, gp, tStep, vF, D);
    this->give3dMaterialStiffnessMatrix_dEdD_from(Ddd, TangentStiffness, gp, tStep, vF, D);
    this->give3dMaterialStiffnessMatrix_dPdD_from(Dud, TangentStiffness, gp, tStep, vF, D);
  } else {
  */
    this->give3dMaterialStiffnessMatrix_dPdF(Duu, TangentStiffness, gp, tStep);
    this->give3dMaterialStiffnessMatrix_dEdD(Ddd, TangentStiffness, gp, tStep);
    this->give3dMaterialStiffnessMatrix_dPdD(Dud, TangentStiffness, gp, tStep);
    
    //}
   FloatMatrix iDdd;
   iDdd.beInverseOf(Ddd);
   double minE = 0;
   double minEig = 0;
   FloatArray n(3);
   FloatArray minN;
   FloatMatrix nn, I(3,3);
   int nStepsI, nStepsJ;
   nStepsI = 80;
   nStepsJ = 40;
   int index = 0;
   for(int iStep = 1; iStep <= nStepsI; iStep ++) {
     double phi = (iStep - 1.)/(nStepsI-1.) * 2. * 3.14159265358979323; 
     for(int jStep = 1; jStep <= nStepsJ; jStep++) {
       double theta = (jStep - 1.)/(nStepsJ-1.) * 3.14159265358979323; 
       n.at(1) = cos(phi) * sin(theta);
       n.at(2) = sin(phi) * sin(theta);
       n.at(3) = cos(theta);
       index++;
       I.beUnitMatrix();
       FloatArray iDddn;
       iDddn.beProductOf(iDdd, n);
       FloatMatrix Omega;
       Omega.beDyadicProductOf(n, iDddn);
       double niDddn = n.dotProduct(iDddn);
       Omega.times(1./ niDddn);
       Omega.subtract(I);
       FloatMatrix iDddOmega;
       iDddOmega.beProductOf(iDdd,Omega);

       FloatMatrix C, Dud_Omega_DudT, Dud_Omega;
       Dud_Omega.beProductOf(Dud,iDddOmega);
       Dud_Omega_DudT.beProductTOf(Dud_Omega,Dud);
       C = Duu;
       C.add(Dud_Omega_DudT);
       

       
       FloatMatrix Q(3,3);
       Q.zero();
       for(int i = 1; i <= 3; i ++) {
	 for(int j = 1; j <= 3; j ++) {
	   for(int k = 1; k <= 3; k ++) {
	     for(int l = 1; l <= 3; l ++) {
	       Q.at(i,k) += C.at(sm::giveVI(i,j), sm::giveVI(k,l)) * n.at(j) * n.at(l);

	     }
	   }
	 }
       }
       double a = Q.at(1,1)/this->mu1;
       double b = (Q.at(1,1)*Q.at(2,2) - Q.at(1,2) * Q.at(2,1))/this->mu1/this->mu1;
       double c = Q.giveDeterminant()/this->mu1/this->mu1/this->mu1;
       double m = min(min(a,b),c);
       FloatArray eval;
       FloatMatrix v;
       this->computeEigs(eval, Q);
       if(index == 1) {
	 minE = m;
	 minEig = min(eval.at(1),min(eval.at(2),eval.at(3)));
	 minQ = Q;
       } else {
	 if( m < minE) {
	   minN = n;
	   minE = m;
	 } 
	 double minEigQ = min(eval.at(1),min(eval.at(2),eval.at(3)));
	 if(minEigQ < minEig) {
	   minEig = minEigQ;
	   minQ = Q;
	 }
       }
       
     }
   }

   //   if(minE < 0) {
     
   if(tStep->giveNumber() >= 1 && gp->giveElement()->giveNumber()==6601) {
       ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
       FloatArray vF = status->giveTempFVector();
       FloatArray D = status->giveTempDVector();
       FloatArray vP = status->giveTempPVector();
       FloatArray E = status->giveTempEVector();
       char s_gp[10];
       sprintf( s_gp, "_gp_%d", gp->giveNumber() );
       //       sprintf(gp, "_gp_%d", gp->giveNumber());
       FILE *FID;  
       std :: string fileName, functionname, temp;
       if(tStep->giveNumber()== 50) {
	 fileName = "Stiffness50";
	 fileName += s_gp;
	 functionname = "Stiffness50";
	 functionname += s_gp;
       } else if(tStep->giveNumber()== 51) {
	 fileName = "Stiffness51";
	 fileName += s_gp;
	 functionname = "Stiffness51";
	 functionname += s_gp;
       } else if(tStep->giveNumber()== 52) {
	 fileName = "Stiffness52";
	 fileName += s_gp;
	 functionname = "Stiffness52";
	 functionname += s_gp;
       } else if(tStep->giveNumber()== 1) {
	 fileName = "Stiffness1";
	 fileName += s_gp;
	 functionname = "Stiffness1";
	 functionname += s_gp;
       }
       fileName += ".m";
       if ( ( FID = fopen(fileName.c_str(), "w") ) == NULL ) {
	 OOFEM_ERROR("failed to open file %s", fileName.c_str() );
       }
       
       fprintf( FID, " function [P, E, F, D, Kuu, Kud, Kdd, n, minEig]= %s \n", functionname.c_str() );
       fprintf(FID, "P=[");
       for(int i = 1; i <= 9; i++) {
	 if(i == 9) {
	   fprintf( FID, "%f;", vP.at(i) );
	 } else {
	   fprintf( FID, "%f, ", vP.at(i) );
	 }
       }
       fprintf(FID, "];\n");
       
       fprintf(FID, "E=[");
       for(int i = 1; i <= 3; i++) {
	 if(i == 3) {
	   fprintf( FID, "%f;", E.at(i) );
	 } else {
	   fprintf( FID, "%f, ", E.at(i) );
	 }
       }
       fprintf(FID, "];\n");
       
       fprintf(FID, "F=[");
       for(int i = 1; i <= 9; i++) {
	 if(i == 9) {
	   fprintf( FID, "%f;", vF.at(i) );
	 } else {
	   fprintf( FID, "%f, ", vF.at(i) );
	 }
       }
       fprintf(FID, "];\n");     
       
       fprintf(FID, "D=[");
       for(int i = 1; i <= 3; i++) {
	 if(i == 3) {
	   fprintf( FID, "%f;", D.at(i) );
	 } else {
	   fprintf( FID, "%f, ", D.at(i) );
	 }
       }
       fprintf(FID, "];\n");
       
       
       fprintf(FID, "Kuu=[");
       for(int i = 1; i <= Duu.giveNumberOfRows(); i++) {
	 for(int j = 1; j <= Duu.giveNumberOfColumns(); j++) {
	   if(j == Duu.giveNumberOfRows()) {
	     fprintf( FID, "%f;", Duu.at(i,j) );
	   } else {
	     fprintf( FID, "%f, ", Duu.at(i,j) );
	   }
	 }
       } 
       fprintf(FID, "];\n");
       
       
       fprintf(FID, "Kud=[");
       for(int i = 1; i <= Dud.giveNumberOfRows(); i++) {
	 for(int j = 1; j <= Dud.giveNumberOfColumns(); j++) {
	   if(j == Dud.giveNumberOfRows()) {
	     fprintf( FID, "%f;", Dud.at(i,j) );
	   } else {
	     fprintf( FID, "%f, ", Dud.at(i,j) );
	   }
	 }
       }
       fprintf(FID, "];\n");
       
       fprintf(FID, "Kdd=[");
       for(int i = 1; i <= Ddd.giveNumberOfRows(); i++) {
	 for(int j = 1; j <= Ddd.giveNumberOfColumns(); j++) {
	   if(j == Duu.giveNumberOfRows()) {
	     fprintf( FID, "%f;", Ddd.at(i,j) );
	   } else {
	     fprintf( FID, "%f, ", Ddd.at(i,j) );
	   }
	 }
       }
       fprintf(FID, "];\n");

       fprintf(FID, "Q=[");
       for(int i = 1; i <= 3; i++) {
	 for(int j = 1; j <= 3; j++) {
	   if(j == 3) {
	     fprintf( FID, "%f;", minQ.at(i,j) );
	   } else {
	     fprintf( FID, "%f, ", minQ.at(i,j) );
	   }
	 }
       }

       
       fprintf(FID, "];\n");
       
       fprintf(FID, "n=[");
       for ( double val: minN ) {
	 fprintf( FID, "%f,", val );
       }
       fprintf(FID, "];\n");

       
       fprintf(FID, "minE=[%f]\n;", minE);
       fprintf(FID, "minEig=[%f]\n;", minEig);
       
       
       fprintf(FID, "end\n");
       
       fclose(FID);   
   }
   if(minE < 0 && minEig < 0) {
     //minE = minEig;
   }
   return minE;

}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: computeAcousticTensor(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{

  FloatMatrix Duu, Ddd, Dud;
  this->give3dMaterialStiffnessMatrix_dPdF(Duu, TangentStiffness, gp, tStep);
  this->give3dMaterialStiffnessMatrix_dEdD(Ddd, TangentStiffness, gp, tStep);
  this->give3dMaterialStiffnessMatrix_dPdD(Dud, TangentStiffness, gp, tStep);
  
  FloatMatrix iDdd;
  iDdd.beInverseOf(Ddd);
  double minE = 0;
  double minEig = 0;
  FloatArray n(3);
  FloatArray minN;
  FloatMatrix nn, I(3,3);
  int nStepsI, nStepsJ;
  nStepsI = 50;
  nStepsJ = 25;
  int index = 0;
  for(int iStep = 1; iStep <= nStepsI; iStep ++) {
    double phi = (iStep - 1.)/(nStepsI-1.) * 2. * 3.14159265358979323; 
    for(int jStep = 1; jStep <= nStepsJ; jStep++) {
      double theta = (jStep - 1.)/(nStepsJ-1.) * 3.14159265358979323; 
      n.at(1) = cos(phi) * sin(theta);
      n.at(2) = sin(phi) * sin(theta);
      n.at(3) = cos(theta);
      index++;
      I.beUnitMatrix();
      FloatArray iDddn;
      iDddn.beProductOf(iDdd, n);
      FloatMatrix Omega;
      Omega.beDyadicProductOf(n, iDddn);
      double niDddn = n.dotProduct(iDddn);
      Omega.times(1./ niDddn);
      Omega.subtract(I);
      FloatMatrix iDddOmega;
      iDddOmega.beProductOf(iDdd,Omega);
      
      FloatMatrix C, Dud_Omega_DudT, Dud_Omega;
      Dud_Omega.beProductOf(Dud,iDddOmega);
      Dud_Omega_DudT.beProductTOf(Dud_Omega,Dud);
      C = Duu;
      C.add(Dud_Omega_DudT);
      
      
      
      FloatMatrix Q(3,3);
      Q.zero();
      for(int i = 1; i <= 3; i ++) {
	for(int j = 1; j <= 3; j ++) {
	  for(int k = 1; k <= 3; k ++) {
	    for(int l = 1; l <= 3; l ++) {
	      Q.at(i,k) += C.at(sm::giveVI(i,j), sm::giveVI(k,l)) * n.at(j) * n.at(l);
	      
	    }
	  }
	}
      }
       double a = Q.at(1,1);
       double b = Q.at(1,1)*Q.at(2,2) - Q.at(1,2) * Q.at(2,1);
       double c = Q.giveDeterminant();
       double m = min(min(a,b),c);
       if(index == 1) {
	 minE = m;
       } else {
	 if(m < minE) {
	   minN = n;
	   minE = m;
	   answer.at(1) = Q.at(1,1);
	   answer.at(2) = Q.at(1,2);
	   answer.at(3) = Q.at(1,3);
	   answer.at(4) = Q.at(1,2);
	   answer.at(5) = Q.at(2,2);
	   answer.at(6) = Q.at(2,3);
	   answer.at(7) = Q.at(1,3);
	   answer.at(8) = Q.at(2,3);
	   answer.at(9) = Q.at(3,3);
	   
	 }
       }
    }
  }
}





void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: computeEigs(FloatArray &eval, const FloatMatrix &Q)
{
  eval.resize(3);
  eval.zero();
  double pi = 3.14159265358979323;
  double p1 = Q.at(1,2)*Q.at(1,2) + Q.at(1,3)*Q.at(1,3) + Q.at(2,3)*Q.at(2,3);
  if(p1==0){
    eval.at(1) = Q.at(1,1);
    eval.at(2) = Q.at(2,2);
    eval.at(3) = Q.at(3,3);
  } else {
    double q = (Q.at(1,1) + Q.at(2,2) + Q.at(3,3))/3.;
    double p2 = (Q.at(1,1)-q)*(Q.at(1,1)-q) + (Q.at(2,2)-q)*(Q.at(2,2)-q) + (Q.at(3,3)-q)*(Q.at(3,3)-q) + 2.*p1;
    double p = sqrt(p2/6);
    FloatMatrix B(Q), I(3,3);
    I.beUnitMatrix();
    I.times(q);
    B.subtract(I);
    B.times(1./p);
    double r = B.giveDeterminant()/2.;
    double phi;
    if(r <= -1) {
      phi = pi/3.;
    }else if(r>=1) {
      phi = 0;
    } else {
      phi = acos(r) / 3.;
    }
    eval.at(1) = q + 2.*p *cos(phi);
    eval.at(2) = q + 2.*p *cos(phi+2./3.*pi);
    eval.at(3) = 3.*q -eval.at(1) - eval.at(2);
    
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
  } else if (type == IST_DamageTensor) {
    answer.resize(9);
    this->computeAcousticTensor(answer, gp, tStep);
  }else {
      return  this->hyperelasticMaterial->giveIPValue(answer, gp, type, tStep);
  }
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dPdD_dEdD(FloatMatrix &dPdD,FloatMatrix &dEdD, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep)
{


  FloatArray vP, vPp, vF, Dp, Ep, E;
  vF.beVectorForm(F);
  this->give_FirstPKStressVector_ElectricalFieldVector_3d(vP, E, gp, vF, D, tStep);

  FloatMatrix Fp;
  auto pert = 1.e-8;
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Invariant F:F
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I1(double J, const FloatArray &vH, const FloatArray &vF)
{
  return vF.dotProduct(vF);    
}



  
/// Invariant H:H
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 ::compute_I2(double J, const FloatArray &vH, const FloatArray &vF)
{
  return vH.dotProduct(vH);    

}

/// Invariant J^{-2/3}F:F
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I1_dev(double J, const FloatArray &vH, const FloatArray &vF)
{
  double I1_dev = this->compute_I1(J,vH,vF);
  I1_dev *= pow(J, -2./3.);
  return I1_dev;    
}

/// Invariant J^{-4/3}F:F
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I2_dev(double J, const FloatArray &vH, const FloatArray &vF)
{
  double I2_dev = this->compute_I2(J,vH,vF);
  I2_dev *= pow(J, -4./3.);
  return I2_dev;
}

/// Invariant (I2_dev)^{3/2}
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I2_dev_pol(double J, const FloatArray &vH, const FloatArray &vF)
{
  double I2_dev_pol = compute_I2_dev(J,vH,vF);
  return pow(I2_dev_pol, 1.5);

}
/// Invariant det F = 1/3 F:H
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I3(double J, const FloatArray &vH, const FloatArray &vF)
{
  return J;
}
    /// Invariant FN \cdot FN
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I4(double J, const FloatArray &vH, const FloatArray &vF)
{
  FloatMatrix F;
  F.beMatrixForm(vF);
  FloatArray FN;
  FN.beProductOf(F, this->N);
  return FN.dotProduct(FN);
}

/// Invariant HN \cdot HN
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I5(double J, const FloatArray &vH, const FloatArray &vF)
{
  FloatMatrix H;
  H.beMatrixForm(vH);
  FloatArray HN;
  HN.beProductOf(H, this->N);
  return HN.dotProduct(HN);
}

/// Invariant D_0 \cdot  D_0
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I6(const FloatArray &D)
{
  return D.dotProduct(D);
}

/// Invariant FD \cdot FD
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I7(double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix F;
  F.beMatrixForm(vF);
  FloatArray d;
  d.beProductOf(F,D);
  return d.dotProduct(d);
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I8(double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix H;
  H.beMatrixForm(vH);
  FloatArray HD;
  HD.beProductOf(H,D);
  return HD.dotProduct(HD);
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_I8_pol(double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  double I2 = this->compute_I2(J, vH, vF);
  double I6 = this->compute_I6(D);
  double I8 = this->compute_I8(J, vH, vF, D);
  return (this->a_I8pol * I2 * I2 + this->b_I8pol * I6 * I6 + this->c_I8pol * I8);
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_K1_Cinf(const FloatArray &D)
{
  return D.dotProduct(this->N);  
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_K2_Cinf(double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix F;
  F.beMatrixForm(vF);
  FloatArray FN;
  FN.beProductOf(F,N);
  //
  FloatArray d;
  d.beProductOf(F,D);
  //
  return d.dotProduct(FN);  
  
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_K2_Cinf_pol(double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  double I4 = this->compute_I4(J, vH, vF);
  double I5 = this->compute_I5(J, vH, vF);
  double I7 = this->compute_I7(J, vH, vF, D);
  double K2_Cinf = this->compute_K2_Cinf(J, vH, vF, D);
  //
  return (this->a_K2_Cinf_pol * this->a_K2_Cinf_pol  * (I4 + I5 - 2.* log(J))  + this->b_K2_Cinf_pol * this->b_K2_Cinf_pol * I7  + 2. * this->a_K2_Cinf_pol * this->b_K2_Cinf_pol * K2_Cinf);
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_K1_Dinf(const FloatArray &D)
{
  double DN = D.dotProduct(this->N);
  return (DN * DN);
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_K2_Dinf(double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix F;
  F.beMatrixForm(vF);
  double DN = D.dotProduct(N);
  FloatArray FD, FN;
  FN.beProductOf(F,N);
  //
  FloatArray d;
  d.beProductOf(F,D);
  //
  double dFN = d.dotProduct(FN);
  return (DN * dFN);

}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_K2_Dinf_pol(double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  double K1_Cinf = compute_K1_Dinf(D);
  double K2_CinfPol = compute_K2_Cinf_pol(J, vH, vF, D);
  return ((this->a_K2_Dinf_pol * K1_Cinf + this->b_K2_Dinf_pol * K2_CinfPol)*(this->a_K2_Dinf_pol * K1_Cinf + this->b_K2_Dinf_pol * K2_CinfPol));
  
}


//////////////////////////////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI1_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF)
{
  dIdF.add(2., vF);
}
  


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI2_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF)
{
  FloatArray FxH;
  hyperelasticMaterial->compute_2order_tensor_cross_product(FxH, vF, vH);
  dIdF.add(2.,FxH);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI3_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF)
{
  dIdF.add(vH);
}
////////////////////////////////////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI1dev_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF)
{
  double FF = vF.dotProduct(vF);
  dIdF.add(2. * pow( J, -2. / 3. ), vF );
  dIdF.add(- 2. / 3. * FF * pow( J , - 5. / 3. ), vH);
}
  
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI2dev_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF)
{
  double HH = vH.dotProduct(vH);
  FloatArray FxH;
  hyperelasticMaterial->compute_2order_tensor_cross_product(FxH, vF, vH);
  ///
  dIdF.add(2. * pow(J , - 4. / 3. ), FxH);
  dIdF.add(-4. / 3. * pow( J , - 7. / 3. ) * HH, vH);  
}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI2devpol_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF)
{
  double I2dev = compute_I2_dev(J, vH, vF);
  //
  FloatArray dI2;
  this->compute_dI2dev_dF(dI2, J, vH, vF);
  //
  dIdF.add(1.5 * sqrt(I2dev), dI2);
}

  
  ////////////////////////////////////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI4_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF)
{
  FloatMatrix F;
  F.beMatrixForm(vF);
  FloatArray FN;
  FN.beProductOf(F,N);
  FloatMatrix FNoN;
  FNoN.beDyadicProductOf(FN, N);
  dIdF.beVectorForm(FNoN);
  dIdF.times(2.);
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI5_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF)
{
  FloatMatrix H;
  H.beMatrixForm(vH); 
  FloatArray HN;
  HN.beProductOf(H,N);
  FloatMatrix HNoN;
  HNoN.beDyadicProductOf(HN, N);
  FloatArray vHNoN;
  vHNoN.beVectorForm(HNoN);
  ///
  hyperelasticMaterial->compute_2order_tensor_cross_product(dIdF, vF, vHNoN);
  dIdF.times(2.);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI6_dD(FloatArray &dIdF, const FloatArray &D)
{
  dIdF.add(2., D);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI7_dF_dD(FloatArray &dIdF, FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix F;
  F.beMatrixForm(vF);
  //
  FloatArray d;
  d.beProductOf(F,D);
  //
  FloatMatrix dD, Ftd;
  dD.beDyadicProductOf(d,D);
  //
  dIdF.beVectorForm(dD);
  dIdF.times(2.);
  //
  dIdD.beTProductOf(F,d);
  dIdD.times(2.);

}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI8_dF_dD(FloatArray &dIdF, FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix H;
  H.beMatrixForm(vH);
  FloatArray HD;
  HD.beProductOf(H,D);
  FloatMatrix HDoD;
  HDoD.beDyadicProductOf(HD, D);
  ///
  FloatArray SigH;
  SigH.beVectorForm(HDoD);
  SigH.times(2.);
  hyperelasticMaterial->compute_2order_tensor_cross_product(dIdF, vF, SigH);
  ///
  FloatMatrix G;
  G.beTProductOf(H,H);
  dIdD.beProductOf(G,D);
  dIdD.times(2.);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI8pol_dF_dD(FloatArray &dIdF, FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  double I2 = this->compute_I2(J, vH, vF);
  double I6 = this->compute_I6(D);
  FloatArray dI2dF, dI6dD, dI8dF, dI8dD;
  this->compute_dI2_dF(dI2dF, J, vH, vF);
  this->compute_dI6_dD(dI6dD, D);
  this->compute_dI8_dF_dD(dI8dF,dI8dD, J, vH, vF, D);
  //
  dIdF.add(2. * this->a_I8pol * I2, dI2dF);
  dIdF.add(this->c_I8pol, dI8dF);
  //
  dIdD.add(2. * this->b_I8pol * I6, dI6dD);
  dIdD.add(this->c_I8pol, dI8dD);

  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK1Cinf_dD(FloatArray &answer, const FloatArray &D)
{
  answer = this->N;
}


//
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2Cinf_dF_dD(FloatArray &dIdF,FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix F;
  F.beMatrixForm(vF);
  //
  FloatArray d;
  d.beProductOf(F,D);
  //
  FloatMatrix dN;
  dN.beDyadicProductOf(d, this->N);
  //
  dIdF.beVectorForm(dN);
  ///
  FloatArray Sigd;
  Sigd.beProductOf(F, this->N);
  //
  FloatMatrix SigdD;
  SigdD.beDyadicProductOf(Sigd,D);
  FloatArray vSigdD;
  vSigdD.beVectorForm(SigdD);
  //
  dIdF.add(vSigdD);
  //
  dIdD.beTProductOf(F,Sigd);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2Cinfpol_dF_dD(FloatArray &dIdF,FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI4dF, dI5dF, dI7dF, dK2CdF;
  FloatArray dI7dD, dK2CdD;
  this->compute_dI4_dF(dI4dF, J, vH, vF);
  //new, is it correct?
  this->compute_dI5_dF(dI5dF, J, vH, vF);
  FloatArray cofF(vH);
  cofF.times(-2./J);
  ///
  this->compute_dI7_dF_dD(dI7dF, dI7dD, J, vH, vF, D);
  this->compute_dK2Cinf_dF_dD(dK2CdF, dK2CdD, J, vH, vF, D);
  //
  dIdF.add(this->a_K2_Cinf_pol * this->a_K2_Cinf_pol, dI4dF);
  /////////////////new/////////////////////////////////////
  dIdF.add(this->a_K2_Cinf_pol * this->a_K2_Cinf_pol, dI5dF);
  dIdF.add(this->a_K2_Cinf_pol * this->a_K2_Cinf_pol, cofF);
  ////////////////////
  //
  dIdF.add(this->b_K2_Cinf_pol * this->b_K2_Cinf_pol, dI7dF);
  dIdF.add(2. * this->a_K2_Cinf_pol * this->b_K2_Cinf_pol, dK2CdF);
  //
  dIdD.add(this->b_K2_Cinf_pol * this->b_K2_Cinf_pol, dI7dD);
  dIdD.add(2. * this->a_K2_Cinf_pol * this->b_K2_Cinf_pol, dK2CdD);

}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK1Dinf_dD(FloatArray &answer, const FloatArray &D)
{
  double DN = D.dotProduct(this->N);
  answer = this->N;
  answer.times(2.*DN); 
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2Dinf_dF_dD(FloatArray &dIdF,FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix F;
  F.beMatrixForm(vF);
  //
  double DN =  D.dotProduct(this->N);
  FloatArray FN;
  FN.beProductOf(F, this->N);
  FloatMatrix FNoD;
  FNoD.beDyadicProductOf(FN, D);
  //
  FloatArray vFNoD;
  vFNoD.beVectorForm(FNoD);
  //
  FloatArray d;
  d.beProductOf(F,D);
  //
  FloatMatrix doN;
  doN.beDyadicProductOf(d, N);
  //
  dIdF.beVectorForm(doN);
  dIdF.add(vFNoD);
  //
  dIdF.times(DN);
  //  
  double dFN = d.dotProduct(FN);
  ///
  dIdD.beTProductOf(F, FN);
  dIdD.times(DN);
  dIdD.add(dFN, this->N);
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2Dinfpol_dF_dD(FloatArray &dIdF,FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  double K1_Cinf = compute_K1_Cinf(D);
  double K2_CinfPol = compute_K2_Cinf_pol(J, vH, vF, D);
  double f = this->a_K2_Dinf_pol * K1_Cinf + this->b_K2_Dinf_pol * K2_CinfPol;
  //
  FloatArray dK2dF, dK2dD, dK1dD;
  this->compute_dK2Cinfpol_dF_dD(dK2dF, dK2dD, J, vH, vF, D);
  //
  this->compute_dK1Cinf_dD(dK1dD, D);
  //
  dIdF.add(2. * f * this->b_K2_Dinf_pol, dK2dF);
  //
  dIdD.add(2. * f * this->a_K2_Dinf_pol, dK1dD);
  dIdD.add(2. * f * this->b_K2_Dinf_pol, dK2dD);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I1_dF2(FloatMatrix &d2IdF2,const double J, const FloatArray &vH, const FloatArray &vF)
{
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  this->hyperelasticMaterial->compute_lower_dyadic_product(d2IdF2, delta, delta);
  d2IdF2.times(2.);
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I2_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF)
{

  FloatMatrix Fx,FxFx,Hx;
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, vF);
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Hx, vH);
  FxFx.beProductOf(Fx,Fx);
  //
  d2IdF2.add(FxFx);
  d2IdF2.add(Hx);
  d2IdF2.times(2.); 
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I3_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF)
{
  hyperelasticMaterial->compute_tensor_cross_product_tensor(d2IdF2, vF);  
}
////////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I1dev_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF)
{
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  //
  double FF = vF.dotProduct(vF);
  //
  this->hyperelasticMaterial->compute_lower_dyadic_product(d2IdF2, delta, delta);
  d2IdF2.times(2. * pow( J, -2. / 3. ));
  //
  FloatMatrix vHvH, vHvF, vFvH, Fx;
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, vF);
  vHvH.beDyadicProductOf(vH,vH);
  vHvF.beDyadicProductOf(vH,vF);
  vFvH.beDyadicProductOf(vF,vH);
  //
  d2IdF2.add(- 4. / 3. * pow( J, -5. / 3. ), vFvH);
  d2IdF2.add(- 4. / 3. * pow( J, -5. / 3. ), vHvF);
  d2IdF2.add( 10. / 9. * FF * pow( J, -8. / 3. ), vHvH);
  d2IdF2.add(- 2. / 3. * FF * pow( J, -5. / 3. ), Fx);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I2dev_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF)
{
  double HH = vH.dotProduct(vH);
  //
  FloatMatrix Hx, FxFx, FxHoH, HoFxH, HoH, Fx;
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, vF);
  FxFx.beProductOf(Fx,Fx);
  //
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Hx, vH);
  //
  FloatArray vFxH;
  hyperelasticMaterial->compute_2order_tensor_cross_product(vFxH, vF, vH);
  //
  FxHoH.beDyadicProductOf(vFxH, vH);
  HoFxH.beDyadicProductOf(vH, vFxH);
  //
  HoH.beDyadicProductOf(vH,vH);
  //  
  d2IdF2.add(2. * pow( J, -4. / 3. ), Hx);
  d2IdF2.add(2. * pow( J, -4. / 3. ), FxFx);
  //
  d2IdF2.add(-8. / 3. * pow( J, -7. / 3. ), FxHoH);
  //
  d2IdF2.add(-8. / 3. * pow( J, -7. / 3. ), HoFxH);
  //
  d2IdF2.add( 28. / 9. * HH * pow( J, -10. / 3. ), HoH);
  //
  d2IdF2.add(- 4. / 3. * HH * pow( J, -7. / 3. ), Fx);

}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I2devpol_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF)
{

  double I2dev = compute_I2_dev(J, vH, vF);
  //
  FloatArray dI2;
  this->compute_dI2dev_dF(dI2, J, vH, vF);
  //
  FloatMatrix dI2dI2;
  dI2dI2.beDyadicProductOf(dI2, dI2);
  //
  FloatMatrix d2I2dF2;
  this->compute_d2I2dev_dF2(d2I2dF2, J, vH, vF);
  //
  d2IdF2.add(1.5 * sqrt(I2dev), d2I2dF2);
  d2IdF2.add(0.75 / sqrt(I2dev), dI2dI2);
}


//////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I4_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF)
{
  //
  FloatMatrix NN;
  NN.beDyadicProductOf(this->N, this->N);
  //
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  //
  this->hyperelasticMaterial->compute_lower_dyadic_product(d2IdF2, delta, NN);
  d2IdF2.times(2.);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I5_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF)
{
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  //
  FloatMatrix NN;
  NN.beDyadicProductOf(this->N, this->N);
  //
  FloatMatrix d2IdH2;
  this->hyperelasticMaterial->compute_lower_dyadic_product(d2IdH2, delta, NN);
  //
  FloatMatrix Fx,FxI;
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, vF);
  FxI.beProductOf(Fx,d2IdH2);
  //
  d2IdF2.beProductOf(FxI,Fx);
  d2IdF2.times(2.);
  //
  FloatMatrix H;
  H.beMatrixForm(vH); 
  FloatMatrix HNoN;
  HNoN.beProductOf(H, NN);
  FloatMatrix HNoNx;
  hyperelasticMaterial->compute_tensor_cross_product_tensor(HNoNx, HNoN);
  //
  d2IdF2.add(2., HNoNx);

}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I6_dD2(FloatMatrix &answer, const FloatArray &D)
{
  answer.resize(3,3);
  answer.beUnitMatrix();
  answer.times(2);
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  //
  FloatMatrix DD;
  DD.beDyadicProductOf(D, D);
  //
  this->hyperelasticMaterial->compute_lower_dyadic_product(d2IdF2, delta, DD);
  d2IdF2.times(2.);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix F;
  F.beMatrixForm(vF);
  d2IdD2.beTProductOf(F,F);
  d2IdD2.times(2.);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  //
  FloatMatrix F;
  F.beMatrixForm(vF);
  FloatArray FD;
  FD.beProductOf(F,D);
  //
  FloatMatrix FD_delta;
  hyperelasticMaterial->compute_3order_dyadic_product(FD_delta, FD, delta);
  //
  FloatMatrix FoD;
  hyperelasticMaterial->compute_3order_lower_dyadic_product(FoD, F, D);
  //
  d2IdFdD.add(2., FD_delta);
  d2IdFdD.add(2., FoD);

  /*  FloatMatrix test, test1, test2;
  test.resize(9,3);
  test1.resize(9,3);
  test2.resize(9,3);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      for(int k = 1; k <= 3; k++) {
	test1.at(sm::giveVI(i,j), k) = (F.at(i,k)*D.at(j));
	test2.at(sm::giveVI(i,j), k) = FD.at(i)*delta.at(j,k);
      }
    }
  }
  test = test1;
  test.add(test2);
  */

}


 void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  FloatMatrix DD;
  DD.beDyadicProductOf(D,D);
  FloatMatrix d2IdH2;
  this->hyperelasticMaterial->compute_lower_dyadic_product(d2IdH2, delta, DD);
  //
  FloatMatrix Fx,FxI;
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, vF);
  FxI.beProductOf(Fx,d2IdH2);
  //
  d2IdF2.beProductOf(FxI,Fx);
  d2IdF2.times(2.);
  //
  FloatMatrix H, HDD;
  H.beMatrixForm(vH);
  HDD.beProductOf(H,DD);
  HDD.times(2.);
  FloatMatrix dIx;
  hyperelasticMaterial->compute_tensor_cross_product_tensor(dIx, HDD);
  //
  d2IdF2.add(dIx);
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix H;
  H.beMatrixForm(vH);
  d2IdD2.beTProductOf(H,H);
  d2IdD2.times(2);
    
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix Fx;
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, vF);
  //
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  //
  FloatMatrix H;
  H.beMatrixForm(vH);
  FloatArray HD;
  HD.beProductOf(H,D);
  //
  FloatMatrix HD_delta;
  hyperelasticMaterial->compute_3order_dyadic_product(HD_delta, HD, delta);
  //
  FloatMatrix HoD;
  hyperelasticMaterial->compute_3order_lower_dyadic_product(HoD, H, D);
  //
  FloatMatrix d2IdHdD;
  d2IdHdD.add(2., HD_delta);
  d2IdHdD.add(2., HoD);
  //
  d2IdFdD.beProductOf(Fx, d2IdHdD);
}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8pol_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
    double I2 = compute_I2(J, vH, vF);
    FloatArray dI2dF;
    this->compute_dI2_dF(dI2dF, J, vH, vF);
    FloatMatrix dI2_dI2;
    dI2_dI2.beDyadicProductOf(dI2dF, dI2dF);
    //
    d2IdF2.add(2. * this->a_I8pol, dI2_dI2);
    //
    FloatMatrix d2I2dF2;
    this->compute_d2I2_dF2(d2I2dF2, J, vH, vF);
    //
    d2IdF2.add(2. * this->a_I8pol * I2, d2I2dF2);
    //
    FloatMatrix d2I8dF2;
    this->compute_d2I8_dF2(d2I8dF2, J, vH, vF, D);
    //
    d2IdF2.add(this->c_I8pol, d2I8dF2);
    
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8pol_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  double I6 = this->compute_I6(D);
  FloatArray dI6dD;
  this->compute_dI6_dD(dI6dD, D);
  FloatMatrix dI6_dI6;
  dI6_dI6.beDyadicProductOf(dI6dD, dI6dD);
  //
  d2IdD2.add(2. * this->b_I8pol, dI6_dI6);
  //
  FloatMatrix d2I6dD2;
  this->compute_d2I6_dD2(d2I6dD2, D);
   //
  d2IdD2.add(2. * this->b_I8pol * I6, d2I6dD2);
  //
  FloatMatrix d2I8dD2;
  this->compute_d2I8_dD2(d2I8dD2, J, vH, vF, D);
  //
  d2IdD2.add(this->c_I8pol, d2I8dD2);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8pol_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  this->compute_d2I8_dFdD(d2IdFdD, J, vH, vF, D);
}

////////////////////////////////////////////////////////////////

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Cinf_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  //
  FloatMatrix ND, DN;
  ND.beDyadicProductOf(N,D);
  DN.beDyadicProductOf(D,N);
  //
  FloatMatrix d2IdF2_1, d2IdF2_2;
  this->hyperelasticMaterial->compute_lower_dyadic_product(d2IdF2_1, delta, ND);
  this->hyperelasticMaterial->compute_lower_dyadic_product(d2IdF2_2, delta, DN);
  d2IdF2 = d2IdF2_1;
  d2IdF2.add(d2IdF2_2);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Cinf_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  d2IdD2.resize(3,3);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Cinf_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix delta(3,3);
  delta.beUnitMatrix();
  //
  FloatMatrix F;
  F.beMatrixForm(vF);
  //
  FloatArray FN;
  FN.beProductOf(F,N);
  //
  FloatMatrix FN_delta;
  hyperelasticMaterial->compute_3order_dyadic_product(FN_delta, FN, delta);
  //
  FloatMatrix FoN;
  hyperelasticMaterial->compute_3order_lower_dyadic_product(FoN, F, N);
  //
  d2IdFdD = FN_delta;
  d2IdFdD.add(FoN);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Cinfpol_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix d2I4dF2, d2I5dF2, Fx, HH, d2I7dF2, d2K2CdF2;
  this->compute_d2I4_dF2(d2I4dF2, J, vH, vF);
  this->compute_d2I5_dF2(d2I5dF2, J, vH, vF);
  hyperelasticMaterial->compute_tensor_cross_product_tensor(Fx, vF);
  HH.beDyadicProductOf(vH,vH);
  this->compute_d2I7_dF2(d2I7dF2, J, vH, vF, D);
  this->compute_d2K2Cinf_dF2(d2K2CdF2, J, vH, vF, D);
  //
  d2IdF2.add(this->a_K2_Cinf_pol * this->a_K2_Cinf_pol, d2I4dF2);
  d2IdF2.add(this->a_K2_Cinf_pol * this->a_K2_Cinf_pol, d2I5dF2);
  d2IdF2.add(-2.0 * this->a_K2_Cinf_pol * this->a_K2_Cinf_pol / J, Fx);
  d2IdF2.add(2.0 * this->a_K2_Cinf_pol * this->a_K2_Cinf_pol / J / J, HH);

  d2IdF2.add(this->b_K2_Cinf_pol * this->b_K2_Cinf_pol, d2I7dF2);
  d2IdF2.add(2. * this->a_K2_Cinf_pol * this->b_K2_Cinf_pol, d2K2CdF2);    
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Cinfpol_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix d2I7dD2, d2K2CdD2;
  this->compute_d2I7_dD2(d2I7dD2, J, vH, vF, D);
  this->compute_d2K2Cinf_dD2(d2K2CdD2, J, vH, vF, D);
  //
  d2IdD2.add(this->b_K2_Cinf_pol * this->b_K2_Cinf_pol, d2I7dD2);
  d2IdD2.add(2. * this->a_K2_Cinf_pol * this->b_K2_Cinf_pol, d2K2CdD2);
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Cinfpol_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix d2I7dFdD, d2K2CdFdD;
  this->compute_d2I7_dFdD(d2I7dFdD, J, vH, vF, D);
  this->compute_d2K2Cinf_dFdD(d2K2CdFdD, J, vH, vF, D);
  //
  d2IdFdD.add(this->b_K2_Cinf_pol * this->b_K2_Cinf_pol, d2I7dFdD);
  d2IdFdD.add(2. * this->a_K2_Cinf_pol * this->b_K2_Cinf_pol, d2K2CdFdD);    
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K1Dinf_dD2(FloatMatrix &d2IdD2,  const FloatArray &D)
{
  d2IdD2.beDyadicProductOf(N,N);
  d2IdD2.times(2.);
}




//
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Dinf_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  double DN = D.dotProduct(N);
  this->compute_d2K2Cinf_dF2(d2IdF2, J, vH, vF, D);
  d2IdF2.times(DN);
}
//
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Dinf_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix F, C;
  F.beMatrixForm(vF);
  C.beTProductOf(F,F);
  FloatArray CN;
  CN.beProductOf(C,N);
  FloatMatrix NoCN, CNoN;
  NoCN.beDyadicProductOf(N,CN);
  CNoN.beDyadicProductOf(CN,N);
  //
  d2IdD2.add(NoCN);
  d2IdD2.add(CNoN);
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Dinf_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatMatrix d2K2dFdD;
  this->compute_d2K2Cinf_dFdD(d2K2dFdD, J, vH, vF, D);
  double DN = D.dotProduct(N);
  //
  d2IdFdD.add(DN, d2K2dFdD);
  //
  FloatArray dKdF, dKdD;
  this->compute_dK2Cinf_dF_dD(dKdF,dKdD, J, vH, vF, D);
  //
  FloatMatrix dKdFoN;
  dKdFoN.beDyadicProductOf(dKdF, N);
  //
  d2IdFdD.add(dKdFoN);
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Dinfpol_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatArray dK2dF,dK2dD;
  FloatMatrix d2K2dF2, dI_dI;
  double K1_Cinf = compute_K1_Dinf(D);
  double K2_CinfPol = compute_K2_Cinf_pol(J, vH, vF, D);
  double f = this->a_K2_Dinf_pol * K1_Cinf + this->b_K2_Dinf_pol * K2_CinfPol;
  //
  this->compute_dK2Cinfpol_dF_dD(dK2dF,dK2dD,J,vH,vF, D);
  dI_dI.beDyadicProductOf(dK2dF,dK2dF);
  //
  d2IdF2.add(2. * this->b_K2_Dinf_pol, dI_dI);
  //
  this->compute_d2K2Cinfpol_dF2(d2K2dF2,J,vH,vF, D);
  d2IdF2.add(2. * f * this->b_K2_Dinf_pol, d2K2dF2);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Dinfpol_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix d2K2dD2, dI_dI;
  double K1_Cinf = compute_K1_Dinf(D);
  double K2_CinfPol = compute_K2_Cinf_pol(J, vH, vF, D);
  double f = this->a_K2_Dinf_pol * K1_Cinf + this->b_K2_Dinf_pol * K2_CinfPol;
  //
  this->compute_d2K2Cinfpol_dD2(d2K2dD2,J,vH,vF,D);
  d2IdD2.add(2. * f * this->b_K2_Dinf_pol, d2K2dD2);
  //
  FloatArray dK1CdD;
  this->compute_dK1Cinf_dD(dK1CdD, D);
  //
  FloatArray dK2dF, dK2CdD;
  this->compute_dK2Cinfpol_dF_dD(dK2dF,dK2CdD,J,vH,vF,D);
  //
  FloatArray dKdD;
  dKdD.add(this->a_K2_Dinf_pol, dK1CdD);
  dKdD.add(this->b_K2_Dinf_pol, dK2CdD);
  //
  FloatMatrix dKdK;
  dKdK.beDyadicProductOf(dKdD,dKdD);
  d2IdD2.add(2., dKdK); 
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2Dinfpol_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix d2K2dFdD;
  double K1_Cinf = compute_K1_Cinf(D);
  double K2_CinfPol = compute_K2_Cinf_pol(J, vH, vF, D);
  double f = this->a_K2_Dinf_pol * K1_Cinf + this->b_K2_Dinf_pol * K2_CinfPol;
  //
  this->compute_d2K2Cinfpol_dFdD(d2K2dFdD,J,vH,vF,D);
  d2IdFdD.add(2. * f * this->b_K2_Dinf_pol, d2K2dFdD);
  //
  FloatArray dK1dD;
  this->compute_dK1Cinf_dD(dK1dD, D);
  //
  FloatArray dK2dF, dK2dD;
  this->compute_dK2Cinfpol_dF_dD(dK2dF,dK2dD,J,vH,vF,D);
  //
  FloatArray dKdD;
  dKdD.add(this->a_K2_Dinf_pol, dK1dD);
  dKdD.add(this->b_K2_Dinf_pol, dK2dD);
  //
  FloatMatrix dKdK;
  dKdK.beDyadicProductOf(dK2dF, dKdD);
  d2IdFdD.add(2. * this->b_K2_Dinf_pol, dKdK); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI1dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_I1(J,vH,vF);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_I1(Jp, vHp, vFp);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI2dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_I2(J,vH,vF);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_I2(Jp, vHp, vFp);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI1devdF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_I1_dev(J,vH,vF);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_I1_dev(Jp, vHp, vFp);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI2devdF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_I2_dev(J,vH,vF);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_I2_dev(Jp, vHp, vFp);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI2devpoldF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_I2_dev_pol(J,vH,vF);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_I2_dev_pol(Jp, vHp, vFp);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI4dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_I4(J,vH,vF);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_I4(Jp, vHp, vFp);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI5dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_I5(J,vH,vF);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_I5(Jp, vHp, vFp);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI7dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_I7(J,vH,vF, D);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_I7(Jp, vHp, vFp, D);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI8dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_I8(J,vH,vF, D);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_I8(Jp, vHp, vFp, D);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI8poldF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_I8_pol(J,vH,vF, D);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_I8_pol(Jp, vHp, vFp, D);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2CdF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_K2_Cinf(J,vH,vF, D);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_K2_Cinf(Jp, vHp, vFp, D);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2CpoldF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_K2_Cinf_pol(J,vH,vF, D);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_K2_Cinf_pol(Jp, vHp, vFp, D);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2DdF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_K2_Dinf(J,vH,vF, D);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_K2_Dinf(Jp, vHp, vFp, D);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2DpoldF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatMatrix Fp, Hp;
  auto I = this->compute_K2_Dinf_pol(J,vH,vF, D);
  auto pert = 1.e-8;
  answer.resize(9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      auto Ip = this->compute_K2_Dinf_pol(Jp, vHp, vFp, D);
      answer.at(sm::giveVI(i,j)) = Ip - I;
    }
  }

  answer.times(1./pert);  
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI6dD_num(FloatArray &answer, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_I6(D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_I6(Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI7dD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_I7(J, vH, vF,D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_I7(J, vH, vF, Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI8dD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_I8(J, vH, vF,D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_I8(J, vH, vF, Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dI8poldD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_I8_pol(J, vH, vF,D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_I8_pol(J, vH, vF, Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK1CdD_num(FloatArray &answer, const FloatArray &D)
{
  FloatArray Dp;
  auto I = this->compute_K1_Cinf(D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_K1_Cinf(Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2CdD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_K2_Cinf(J, vH, vF, D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_K2_Cinf(J, vH, vF, Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2CpoldD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_K2_Cinf_pol(J, vH, vF,D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_K2_Cinf_pol(J, vH, vF, Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK1DdD_num(FloatArray &answer, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_K1_Dinf(D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_K1_Dinf(Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2DdD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_K2_Dinf(J, vH, vF,D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_K2_Dinf(J, vH, vF, Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_dK2DpoldD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray Dp;
  auto I = this->compute_K2_Dinf_pol(J, vH, vF,D);
  auto pert = 1.e-8;
  answer.resize(3);
  for(int i = 1; i <= 3; i++) {
      Dp = D;
      Dp.at(i) += pert;
      auto Ip = this->compute_K2_Dinf_pol(J, vH, vF, Dp);
      answer.at(i) = Ip - I;
  }
  answer.times(1./pert);  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I1dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{


  FloatArray dI, dIp;
  this->compute_dI1_dF(dI, J, vH, vF);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dI1_dF(dIp,Jp, vHp, vFp);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I1devdF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{


  FloatArray dI, dIp;
  this->compute_dI1dev_dF(dI, J, vH, vF);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero(); 
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dI1dev_dF(dIp,Jp, vHp, vFp);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I2dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{


  FloatArray dI, dIp;
  this->compute_dI2_dF(dI, J, vH, vF);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dI2_dF(dIp,Jp, vHp, vFp);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I2devdF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{


  FloatArray dI, dIp;
  this->compute_dI2dev_dF(dI, J, vH, vF);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dI2dev_dF(dIp,Jp, vHp, vFp);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I2devpoldF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{


  FloatArray dI, dIp;
  this->compute_dI2devpol_dF(dI, J, vH, vF);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dI2devpol_dF(dIp,Jp, vHp, vFp);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I4dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{


  FloatArray dI, dIp;
  this->compute_dI4_dF(dI, J, vH, vF);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dI4_dF(dIp,Jp, vHp, vFp);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I5dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF)
{


  FloatArray dI, dIp;
  this->compute_dI5_dF(dI, J, vH, vF);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dI5_dF(dIp,Jp, vHp, vFp);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{


  FloatArray dI, dIp, dIdD;
  this->compute_dI7_dF_dD(dI, dIdD, J, vH, vF, D);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dI7_dF_dD(dIp,dIdD, Jp, vHp, vFp, D);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI, dIp, dIdD;
  this->compute_dI8_dF_dD(dI, dIdD, J, vH, vF, D);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dI8_dF_dD(dIp,dIdD, Jp, vHp, vFp, D);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8poldF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI, dIp, dIdD;
  this->compute_dI8pol_dF_dD(dI, dIdD, J, vH, vF, D);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dI8pol_dF_dD(dIp,dIdD, Jp, vHp, vFp, D);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2CdF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI, dIp, dIdD;
  this->compute_dK2Cinf_dF_dD(dI, dIdD, J, vH, vF, D);
  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dK2Cinf_dF_dD(dIp,dIdD, Jp, vHp, vFp, D);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2CpoldF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI, dIp, dIdD;
  this->compute_dK2Cinfpol_dF_dD(dI, dIdD, J, vH, vF, D);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dK2Cinfpol_dF_dD(dIp,dIdD, Jp, vHp, vFp, D);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2DdF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI, dIp, dIdD;
  this->compute_dK2Dinf_dF_dD(dI, dIdD, J, vH, vF, D);
  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dK2Dinf_dF_dD(dIp,dIdD,Jp, vHp, vFp, D);
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2DpoldF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI, dIp, dIdD;
  this->compute_dK2Dinfpol_dF_dD(dI, dIdD, J, vH, vF, D);

  FloatMatrix Fp, Hp;
  auto pert = 1.e-8;
  answer.resize(9,9);
  for(int i = 1; i <= 3; i++) {
    for(int j = 1; j <= 3; j++) {
      dIp.zero();
      Fp.beMatrixForm(vF);
      Fp.at(i,j) += pert;
      hyperelasticMaterial->compute_2order_tensor_cross_product(Hp,Fp,Fp);
      Hp.times(0.5);
      FloatArray vFp,vHp;
      double Jp;
      vFp.beVectorForm(Fp);
      vHp.beVectorForm(Hp);
      Jp = 1./3. * vFp.dotProduct(vHp);
      this->compute_dK2Dinfpol_dF_dD(dIp,dIdD, Jp, vHp, vFp, D);
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {
	  answer.at(sm::giveVI(i,j), sm::giveVI(k,l)) = dIp.at(sm::giveVI(k,l)) -  dI.at(sm::giveVI(k,l));
	}
      }
    }
  }

  answer.times(1./pert);  
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I6dD2_num(FloatMatrix &answer, const FloatArray &D)
{

  FloatArray dI, dIp;
  this->compute_dI6_dD(dI, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(3,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dI6_dD(dIp, Dp);
      for(int k = 1; k <= 3; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert); 
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7dD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatArray dI, dIp, dIdF;
  this->compute_dI7_dF_dD(dIdF, dI, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(3,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dI7_dF_dD(dIdF, dIp, J, vH, vF, Dp);
      for(int k = 1; k <= 3; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8dD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI, dIp, dIdF;
  this->compute_dI8_dF_dD(dIdF, dI, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-5;
  answer.resize(3,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dI8_dF_dD(dIdF, dIp, J, vH, vF, Dp);
      for(int k = 1; k <= 3; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8poldD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI, dIp, dIdF;
  this->compute_dI8pol_dF_dD(dIdF, dI, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(3,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dI8pol_dF_dD(dIdF, dIp, J, vH, vF, Dp);
      for(int k = 1; k <= 3; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2CdD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI, dIp, dIdF;
  this->compute_dK2Cinf_dF_dD(dIdF, dI, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(3,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dK2Cinf_dF_dD(dIdF, dIp, J, vH, vF, Dp);
      for(int k = 1; k <= 3; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2CpoldD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatArray dI, dIp, dIdF;
  this->compute_dK2Cinfpol_dF_dD(dIdF, dI, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(3,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dK2Cinfpol_dF_dD(dIdF, dIp, J, vH, vF, Dp);
      for(int k = 1; k <= 3; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);

}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K1DdD2_num(FloatMatrix &answer, const FloatArray &D)
{
  FloatArray dI, dIp;
  this->compute_dK1Dinf_dD(dI, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(3,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dK1Dinf_dD(dIp, Dp);
      for(int k = 1; k <= 3; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2DdD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatArray dI, dIp, dIdF;
  this->compute_dK2Dinf_dF_dD(dIdF, dI, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(3,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dK2Dinf_dF_dD(dIdF, dIp, J, vH, vF, Dp);
      for(int k = 1; k <= 3; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2DpoldD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatArray dI, dIp, dIdF;
  this->compute_dK2Dinfpol_dF_dD(dIdF, dI, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(3,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dK2Dinfpol_dF_dD(dIdF, dIp, J, vH, vF, Dp);
      for(int k = 1; k <= 3; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);  
}
///////////////////////////////////////////////////////////////////////////////////////////////////

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I7dFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatArray dI, dIp, dIdD;
  this->compute_dI7_dF_dD(dI, dIdD, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dI7_dF_dD(dIp, dIdD, J, vH, vF, Dp);
      for(int k = 1; k <= 3; k++) {
	for(int l = 1; l <= 3; l++) {
	  answer.at(sm::giveVI(k,l),i) = dIp.at(sm::giveVI(k,l)) -  dI.at(sm::giveVI(k,l));
        }
      }
  }
  
  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8dFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatArray dI, dIp, dIdD;
  this->compute_dI8_dF_dD(dI, dIdD, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dI8_dF_dD(dIp, dIdD, J, vH, vF, Dp);
      for(int k = 1; k <= 9; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2I8poldFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatArray dI, dIp, dIdD;
  this->compute_dI8pol_dF_dD(dI, dIdD, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dI8pol_dF_dD(dIp, dIdD, J, vH, vF, Dp);
      for(int k = 1; k <= 9; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);  
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2CdFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  
  FloatArray dI, dIp, dIdD;
  this->compute_dK2Cinf_dF_dD(dI, dIdD, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dK2Cinf_dF_dD(dIp, dIdD, J, vH, vF, Dp);
      for(int k = 1; k <= 9; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert); 
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2CpoldFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  
  FloatArray dI, dIp, dIdD;
  this->compute_dK2Cinfpol_dF_dD(dI, dIdD, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dK2Cinfpol_dF_dD(dIp, dIdD, J, vH, vF, Dp);
      for(int k = 1; k <= 9; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert); 
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2DdFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{

  FloatArray dI, dIp, dIdD;
  this->compute_dK2Dinf_dF_dD(dI, dIdD, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dK2Dinf_dF_dD(dIp, dIdD, J, vH, vF, Dp);
      for(int k = 1; k <= 9; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);  
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :: compute_d2K2DpoldFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D)
{
  FloatArray dI, dIp, dIdD;
  this->compute_dK2Dinfpol_dF_dD(dI, dIdD, J, vH, vF, D);
  
  FloatArray Dp;
  auto pert = 1.e-8;
  answer.resize(9,3);
  for(int i = 1; i <= 3; i++) {
      dIp.zero();
      Dp = D;
      Dp.at(i) += pert;
      this->compute_dK2Dinfpol_dF_dD(dIp, dIdD, J, vH, vF, Dp);
      for(int k = 1; k <= 9; k++) {
	answer.at(k,i) = dIp.at(k) -  dI.at(k);
      }
  }
  
  answer.times(1./pert);  
}





/*
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );   
    //////////////////////
    vP.zero();
    E.zero();
    FloatMatrix F;
    F.beMatrixForm(vF2);
    //@todo: testing
    F = {{1.2, 0.3, -0.1},{-0.5, 0.2, 0.456},{-0.65,0.2, 1.4}};
    FloatArray vF;
    vF.beVectorForm(F);
    FloatArray D;
    D = {0.1, 4, 51};
    // H stands for cofactor of F
    FloatArray vH;
    // compute cofactor using tensor cross product
    this->hyperelasticMaterial->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double J = 1./3. * vH.dotProduct( vF );
    //
    FloatMatrix H;
    H.beMatrixForm(vH);
    
    ////////////////    
    // First Piol-Kirchhoff stress
    FloatArray dI1dF, dI2dF;
    this->compute_dI1_dF(dI1dF, J, vH, vF);
    this->compute_dI2_dF(dI2dF, J, vH, vF);
    ///
    FloatArray dI1dFnum, dI2dFnum;
    this->compute_dI1dF_num(dI1dFnum, J, vH, vF);
    this->compute_dI2dF_num(dI2dFnum, J, vH, vF);
    ////////////////////////////////////////////
    FloatArray dI1devdF, dI2devdF, dI2devpoldF;
    this->compute_dI1dev_dF(dI1devdF, J, vH, vF);
    this->compute_dI2dev_dF(dI2devdF, J, vH, vF);
    this->compute_dI2devpol_dF(dI2devpoldF, J, vH, vF);
    ///
    FloatArray dI1devdFnum, dI2devdFnum, dI2devpoldFnum;
    this->compute_dI1devdF_num(dI1devdFnum, J, vH, vF);
    this->compute_dI2devdF_num(dI2devdFnum, J, vH, vF);
    this->compute_dI2devpoldF_num(dI2devpoldFnum, J, vH, vF);
    /////////////////////////////////////////////////
    FloatArray dI4dF, dI5dF;
    this->compute_dI4_dF(dI4dF, J, vH, vF);
    this->compute_dI5_dF(dI5dF, J, vH, vF);
    //
    FloatArray dI4dFnum, dI5dFnum;
    this->compute_dI4dF_num(dI4dFnum, J, vH, vF);
    this->compute_dI5dF_num(dI5dFnum, J, vH, vF);
    /////////////////////////////////////////////////
    FloatArray dI7dF,dI7dD, dI8dF,dI8dD, dI8poldF,dI8poldD;
    this->compute_dI7_dF_dD(dI7dF, dI7dD, J, vH, vF, D);
    this->compute_dI8_dF_dD(dI8dF, dI8dD, J, vH, vF, D);
    this->compute_dI8pol_dF_dD(dI8poldF, dI8poldD, J, vH, vF, D);
    //
    FloatArray dI7dFnum,dI7dDnum, dI8dFnum,dI8dDnum, dI8poldFnum,dI8poldDnum;
    this->compute_dI7dF_num(dI7dFnum, J, vH, vF, D);
    this->compute_dI8dF_num(dI8dFnum, J, vH, vF, D);
    this->compute_dI8poldF_num(dI8poldFnum, J, vH, vF, D);
    this->compute_dI7dD_num(dI7dDnum, J, vH, vF, D);
    this->compute_dI8dD_num(dI8dDnum, J, vH, vF, D);
    this->compute_dI8poldD_num(dI8poldDnum, J, vH, vF, D);
    ///////////////////////////////////////////////// 
    FloatArray dK2CdF,dK2CdD, dK2CpoldF, dK2CpoldD, dK2DdF,dK2DdD, dK2DpoldF, dK2DpoldD;
    this->compute_dK2Cinf_dF_dD(dK2CdF,dK2CdD, J, vH, vF, D);
    this->compute_dK2Cinfpol_dF_dD(dK2CpoldF,dK2CpoldD, J, vH, vF, D);
    this->compute_dK2Dinf_dF_dD(dK2DdF,dK2DdD, J, vH, vF, D);
    this->compute_dK2Cinfpol_dF_dD(dK2DpoldF,dK2DpoldD,J, vH, vF, D);
    //
    FloatArray dK2CdFnum,dK2CdDnum, dK2CpoldFnum, dK2CpoldDnum, dK2DdFnum,dK2DdDnum, dK2DpoldFnum, dK2DpoldDnum;
    this->compute_dK2CdF_num(dK2CdFnum, J, vH, vF, D);
    this->compute_dK2CpoldF_num(dK2CpoldFnum, J, vH, vF, D);
    this->compute_dK2DdF_num(dK2DdFnum, J, vH, vF, D);
    this->compute_dK2CpoldF_num(dK2DpoldFnum,J, vH, vF, D);
    //
    this->compute_dK2CdD_num(dK2CdDnum, J, vH, vF, D);
    this->compute_dK2CpoldD_num(dK2CpoldDnum, J, vH, vF, D);
    this->compute_dK2DdD_num(dK2DdDnum, J, vH, vF, D);
    this->compute_dK2CpoldD_num(dK2DpoldDnum,J, vH, vF, D);
    ///////////////////////////////////
    FloatArray dI6dD, dK1CdD, dK1DdD;
    this->compute_dI6_dD(dI6dD, D);
    this->compute_dK1Cinf_dD(dK1CdD, D);
    this->compute_dK1Dinf_dD(dK1DdD, D);
    //
    FloatArray dI6dDnum, dK1CdDnum, dK1DdDnum;
    this->compute_dI6dD_num(dI6dDnum, D);
    this->compute_dK1CdD_num(dK1CdDnum, D);
    this->compute_dK1DdD_num(dK1DdDnum, D);
    ////////////////////////////////////////////////////////////////////////
    //second derivatives
    FloatMatrix d2I1dF2, d2I2dF2;
    this->compute_d2I1_dF2(d2I1dF2, J, vH, vF);
    this->compute_d2I2_dF2(d2I2dF2, J, vH, vF);
    ///
    FloatMatrix d2I1dF2num, d2I2dF2num;
    this->compute_d2I1dF2_num(d2I1dF2num, J, vH, vF);
    this->compute_d2I2dF2_num(d2I2dF2num, J, vH, vF);
    ////////////////////////////////////////////
    FloatMatrix d2I1devdF2, d2I2devdF2, d2I2devpoldF2;
    this->compute_d2I1dev_dF2(d2I1devdF2, J, vH, vF);
    this->compute_d2I2dev_dF2(d2I2devdF2, J, vH, vF);
    this->compute_d2I2devpol_dF2(d2I2devpoldF2, J, vH, vF);
    ///
    FloatMatrix d2I1devdF2num, d2I2devdF2num, d2I2devpoldF2num;
    this->compute_d2I1devdF2_num(d2I1devdF2num, J, vH, vF);
    this->compute_d2I2devdF2_num(d2I2devdF2num, J, vH, vF);
    this->compute_d2I2devpoldF2_num(d2I2devpoldF2num, J, vH, vF);
    /////////////////////////////////////////////////
    FloatMatrix d2I4dF2, d2I5dF2;
    this->compute_d2I4_dF2(d2I4dF2, J, vH, vF);
    this->compute_d2I5_dF2(d2I5dF2, J, vH, vF);
    //
    FloatMatrix d2I4dF2num, d2I5dF2num;
    this->compute_d2I4dF2_num(d2I4dF2num, J, vH, vF);
    this->compute_d2I5dF2_num(d2I5dF2num, J, vH, vF);
    /////////////////////////////////////////////////
    FloatMatrix d2I7dF2,d2I8dF2,d2I8poldF2;
    this->compute_d2I7_dF2(d2I7dF2, J, vH, vF, D);
    this->compute_d2I8_dF2(d2I8dF2, J, vH, vF, D);
    this->compute_d2I8pol_dF2(d2I8poldF2, J, vH, vF, D);
    //
    FloatMatrix d2I7dF2num, d2I8dF2num, d2I8poldF2num;
    this->compute_d2I7dF2_num(d2I7dF2num, J, vH, vF, D);
    this->compute_d2I8dF2_num(d2I8dF2num, J, vH, vF, D);
    this->compute_d2I8poldF2_num(d2I8poldF2num, J, vH, vF, D);
    ///////////////////////////////////////////////// 
    FloatMatrix d2K2CdF2, d2K2CpoldF2, d2K2DdF2, d2K2DpoldF2;
    this->compute_d2K2Cinf_dF2(d2K2CdF2, J, vH, vF, D);
    this->compute_d2K2Cinfpol_dF2(d2K2CpoldF2, J, vH, vF, D);
    this->compute_d2K2Dinf_dF2(d2K2DdF2, J, vH, vF, D);
    this->compute_d2K2Dinfpol_dF2(d2K2DpoldF2,J, vH, vF, D);
    //
    FloatMatrix d2K2CdF2num, d2K2CpoldF2num, d2K2DdF2num, d2K2DpoldF2num;
    this->compute_d2K2CdF2_num(d2K2CdF2num, J, vH, vF, D);
    this->compute_d2K2CpoldF2_num(d2K2CpoldF2num, J, vH, vF, D);
    this->compute_d2K2DdF2_num(d2K2DdF2num, J, vH, vF, D);
    this->compute_d2K2DpoldF2_num(d2K2DpoldF2num,J, vH, vF, D);
    /////////////////////////////////////////////////
    ///second derivatie wrt D
    FloatMatrix d2I7dD2,d2I8dD2,d2I8poldD2;
    this->compute_d2I7_dD2(d2I7dD2, J, vH, vF, D);
    this->compute_d2I8_dD2(d2I8dD2, J, vH, vF, D);
    this->compute_d2I8pol_dD2(d2I8poldD2, J, vH, vF, D);
    //
    FloatMatrix d2I7dD2num, d2I8dD2num, d2I8poldD2num;
    this->compute_d2I7dD2_num(d2I7dD2num, J, vH, vF, D);
    this->compute_d2I8dD2_num(d2I8dD2num, J, vH, vF, D);
    this->compute_d2I8poldD2_num(d2I8poldD2num, J, vH, vF, D);
    ///////////////////////////////////////////////// 
    FloatMatrix d2K2CdD2, d2K2CpoldD2, d2K2DdD2, d2K2DpoldD2;
    this->compute_d2K2Cinf_dD2(d2K2CdD2, J, vH, vF, D);
    this->compute_d2K2Cinfpol_dD2(d2K2CpoldD2, J, vH, vF, D);
    this->compute_d2K2Dinf_dD2(d2K2DdD2, J, vH, vF, D);
    this->compute_d2K2Cinfpol_dD2(d2K2DpoldD2,J, vH, vF, D);
    //
    FloatMatrix d2K2CdD2num, d2K2CpoldD2num, d2K2DdD2num, d2K2DpoldD2num;
    this->compute_d2K2CdD2_num(d2K2CdD2num, J, vH, vF, D);
    this->compute_d2K2CpoldD2_num(d2K2CpoldD2num, J, vH, vF, D);
    this->compute_d2K2DdD2_num(d2K2DdD2num, J, vH, vF, D);
    this->compute_d2K2CpoldD2_num(d2K2DpoldD2num,J, vH, vF, D);
    ///////////////////////////////////
    FloatMatrix d2I6dD2, d2K1CdD2, d2K1DdD2;
    this->compute_d2I6_dD2(d2I6dD2, D);
    // this->compute_d2K1Cinf_dD2(d2K1CdD2, D);
    this->compute_d2K1Dinf_dD2(d2K1DdD2, D);
    //
    FloatMatrix d2I6dD2num, d2K1CdD2num, d2K1DdD2num;
    this->compute_d2I6dD2_num(d2I6dD2num, D);
    //this->compute_d2K1CdD2_num(d2K1CdD2num, D);
    this->compute_d2K1DdD2_num(d2K1DdD2num, D);
    /////////////////////////////////////////////////
    ///second derivatie wrt FD
    FloatMatrix d2I7dFdD,d2I8dFdD,d2I8poldFdD;
    this->compute_d2I7_dFdD(d2I7dFdD, J, vH, vF, D);
    this->compute_d2I8_dFdD(d2I8dFdD, J, vH, vF, D);
    this->compute_d2I8pol_dFdD(d2I8poldFdD, J, vH, vF, D);
    //
    FloatMatrix d2I7dFdDnum, d2I8dFdDnum, d2I8poldFdDnum;
    this->compute_d2I7dFdD_num(d2I7dFdDnum, J, vH, vF, D);
    this->compute_d2I8dFdD_num(d2I8dFdDnum, J, vH, vF, D);
    this->compute_d2I8poldFdD_num(d2I8poldFdDnum, J, vH, vF, D);
    ///////////////////////////////////////////////// 
    FloatMatrix d2K2CdFdD, d2K2CpoldFdD, d2K2DdFdD, d2K2DpoldFdD;
    this->compute_d2K2Cinf_dFdD(d2K2CdFdD, J, vH, vF, D);
    this->compute_d2K2Cinfpol_dFdD(d2K2CpoldFdD, J, vH, vF, D);
    this->compute_d2K2Dinf_dFdD(d2K2DdFdD, J, vH, vF, D);
    this->compute_d2K2Dinfpol_dFdD(d2K2DpoldFdD,J, vH, vF, D);
    //
    FloatMatrix d2K2CdFdDnum, d2K2CpoldFdDnum, d2K2DdFdDnum, d2K2DpoldFdDnum;
    this->compute_d2K2CdFdD_num(d2K2CdFdDnum, J, vH, vF, D);
    this->compute_d2K2CpoldFdD_num(d2K2CpoldFdDnum, J, vH, vF, D);
    this->compute_d2K2DdFdD_num(d2K2DdFdDnum, J, vH, vF, D);
    this->compute_d2K2DpoldFdD_num(d2K2DpoldFdDnum,J, vH, vF, D);
*/


} // end namespace oofem
