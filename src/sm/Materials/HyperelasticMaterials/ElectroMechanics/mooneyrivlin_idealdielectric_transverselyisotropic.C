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

#include "../sm/Materials/HyperelasticMaterials/ElectroMechanics/mooneyrivlin_idealdielectric_transverselyisotropic.h"
#include "../sm/Materials/structuralmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"

using sm = oofem::StructuralMaterial;


namespace oofem {
  REGISTER_Material(MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial);

  MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial(int n, Domain *d) : Material(n, d), ElectroMechanicalMaterialExtensionInterface_3Fields(d)
{
    this->hyperelasticMaterial = new MooneyRivlinMaterial(n, d);
}


IRResultType
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro


    result = this->hyperelasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) {
      return result;
    }

    IR_GIVE_FIELD(ir, this->epsilon_m, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial_epsilon_m);
    IR_GIVE_FIELD(ir, this->epsilon_f, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial_epsilon_f);
    IR_GIVE_FIELD(ir, this->N, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial_n);
    N.normalize();
    IR_GIVE_FIELD(ir, this->C4, _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial_C4);
      
    return IRRT_OK;
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: give_FirstPKStressVector_ElectricalFieldVector_3d(FloatArray &vP, FloatArray &E, GaussPoint *gp, const FloatArray &vF, const FloatArray &D, TimeStep *tStep)
{

    vP.zero();
    E.zero();
    FloatMatrix F;
    F.beMatrixForm(vF);
    // First Piol-Kirchhoff stress
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
    //
    double I4 = this->compute_I4(F);
    double I7 = this->compute_I7(F, D);
    double I9 = this->compute_I9(D);
    //double I10 = this->compute_I10(F, D);
    //derivatives wrt F
    FloatArray dI4dF,dI7dF,dI7dD,dI9dD;
    this->compute_dI4dF(dI4dF, F);
    this->compute_dI7dF(dI7dF, F, D);
    //this->compute_dI10dF(dI10dF, F, D);
    //derivatives wrt D_0
    this->compute_dI7dD(dI7dD, F, D);
    this->compute_dI9dD(dI9dD, D);
    //this->compute_dI10dD(dI10dD, F, D);
    //
    vP.add( C4 * (I4 - 1.), dI4dF );
    if(this->epsilon_m != 0) {
      vP.add( 1. / this->epsilon_m / J, dI7dF );
      vP.add( - I7 / J / J / this->epsilon_m, vH);
      E.add( 1. / J / this->epsilon_m , dI7dD);
    }

    if(this->epsilon_f != 0) {
      E.add( 2. * I9 / this->epsilon_f, dI9dD);    
    }
    

    
    status->letTempPVectorBe(vP);
    status->letTempFVectorBe(vF);
    status->letTempEVectorBe(E);
    status->letTempDVectorBe(D);
}


  
  

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
 
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    hyperelasticMaterial->give3dMaterialStiffnessMatrix_dPdF(answer, mode, gp, tStep); 
    FloatArray vF = status->giveTempFVector();
    FloatArray E = status->giveTempEVector();
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
    auto I4 = this->compute_I4(F);
    auto I7 = this->compute_I7(F, D);
    //auto I10 = this->compute_I10(F, D);
    FloatArray dI4dF, dI7dF;
    FloatMatrix d2I4dF2, d2I7dF2;
    this->compute_dI4dF(dI4dF, F);
    this->compute_d2I4dF2(d2I4dF2, F);
    //
    this->compute_dI7dF(dI7dF, F, D);
    this->compute_d2I7dF2(d2I7dF2, F, D);
    //this->compute_dI10dF(dI10dF, F, D);
    //this->compute_d2I10dF2(d2I10dF2, F, D);
    FloatMatrix cFcF, cFdI7, dI7cF, dI4dI4;
    dI4dI4.beDyadicProductOf(dI4dF, dI4dF);
    cFcF.beDyadicProductOf(vH, vH);

    dI7cF.beDyadicProductOf(dI7dF, vH);
    cFdI7.beDyadicProductOf(vH, dI7dF);
    
    answer.add(this->C4, dI4dI4);
    answer.add(this->C4 * (I4 - 1.), d2I4dF2);

        
    if(this->epsilon_m != 0) {
      answer.add(1. / this->epsilon_m / J, d2I7dF2);
      answer.add(2. * I7 / this->epsilon_m / J / J / J, cFcF);
      answer.add(- I7 / this->epsilon_m / J / J , Fx);
      answer.add(- 1. / this->epsilon_m / J / J , cFdI7);
      answer.add(- 1. / this->epsilon_m / J / J , dI7cF);
    }

    /*if(this->epsilon_f != 0) {
      answer.add(2. / this->epsilon_f / J, dI10dI10);
      answer.add( 2. * I10 / this->epsilon_f / J, d2I10dF2);
      answer.add(2. * I10 * I10 / this->epsilon_f / J / J / J, cFcF);
      answer.add(- I10 * I10 / this->epsilon_f / J / J , Fx);
      answer.add(- 2. * I10 / this->epsilon_f / J / J , cFdI10);
      answer.add(- 2. * I10 / this->epsilon_f / J / J , dI10cF);

      }*/

    /*
    
    FloatMatrix num_d2I4dF2,num_d2I7dF2,num_d2I10dF2;
    FloatArray num_dI10dF;
    this->compute_d2I4dF2_num(num_d2I4dF2, F);
    this->compute_d2I7dF2_num(num_d2I7dF2, F, D);
    this->compute_d2I10dF2_num(num_d2I10dF2, F, D);
    this->compute_dI10dF_num(num_dI10dF, F, D);
    
    FloatMatrix dPdF, dEdF;
    this->compute_dPdF_dEdF(dPdF, dEdF, F, D, gp, tStep);

    //this->compute_dPdF_dEdF_fake(dPdF, dEdF, F, D, gp, tStep);
    */

    
}



void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: give3dMaterialStiffnessMatrix_dPdD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
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
    //auto I10 = this->compute_I10(F,D);
    FloatArray dI7dD;
    FloatMatrix d2I7dFdD;
    this->compute_d2I7dFdD(d2I7dFdD, F, D);
    //this->compute_d2I10dFdD(d2I10dFdD, F, D);
    this->compute_dI7dD(dI7dD, F, D);
    //this->compute_dI10dF(dI10dF, F, D);
    //this->compute_dI10dD(dI10dD, F, D);
    FloatMatrix  cF_dI7dD;
    cF_dI7dD.beDyadicProductOf(vH,dI7dD);
    
    if(this->epsilon_m != 0) {
      answer.add(1. / this->epsilon_m / J, d2I7dFdD);
      answer.add(- 1. / J / J / this->epsilon_m, cF_dI7dD);
    }
    

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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: give3dMaterialStiffnessMatrix_dEdD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  answer.resize(3,3);
    ElectroMechanicalMaterialStatus *status = static_cast< ElectroMechanicalMaterialStatus* >( this->giveStatus(gp) );
    FloatArray vF = status->giveTempFVector();
    FloatArray D = status->giveTempDVector();
    FloatMatrix F;
    F.beMatrixForm(vF);
    double I9 = this->compute_I9(D);
    FloatArray dI9dD;
    FloatMatrix d2I7dD2, d2I9dD2, dI9dI9;
    this->compute_d2I7dD2(d2I7dD2, F, D);
    this->compute_dI9dD(dI9dD, D);
    dI9dI9.beDyadicProductOf(dI9dD, dI9dD);
    double J = F.giveDeterminant();
    if(this->epsilon_m != 0) {
      answer.add(1./J / this->epsilon_m, d2I7dD2);
    }
    if(this->epsilon_f != 0) {
      answer.add( 2. / this->epsilon_f, dI9dI9);
    }

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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_I4(const FloatMatrix &F)
{
  
  FloatArray CN;
  FloatMatrix C;
  C.beTProductOf(F,F);
  CN.beProductOf(C,this->N);
  return this->N.dotProduct(CN);  
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_I7(const FloatMatrix &F, const FloatArray &D)
{
  
  FloatArray CD;
  FloatMatrix C;
  C.beTProductOf(F,F);
  CD.beProductOf(C,D);
  return D.dotProduct(CD);  
}

double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_I9(const FloatArray &D)
{
  return this->N.dotProduct(D);  
}


  
double
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_I10(const FloatMatrix &F, const FloatArray &D)
{
  
  FloatArray CD;
  FloatMatrix C;
  C.beTProductOf(F,F);
  CD.beProductOf(C,D);
  return this->N.dotProduct(CD);  
}





void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI4dF(FloatArray &answer, const FloatMatrix &F)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI7dF(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI7dD(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(9);

  FloatMatrix C;
  C.beTProductOf(F,F);
  
  answer.beProductOf(C,D);
  answer.times(2.);
  
}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI9dD(FloatArray &answer, const FloatArray &D)
{
  answer = N;
}


void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI10dF(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI10dD(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
{

  answer.resize(9);

  FloatMatrix C;
  C.beTProductOf(F,F);  
  answer.beProductOf(C,N);
}

///
void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I4dF2(FloatMatrix &answer, const FloatMatrix &F)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I7dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I10dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I7dD2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
{
  answer.beTProductOf(F,F);
  answer.times(2);
}

void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I7dFdD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I10dFdD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: computeAcousticTensorMinEigenvalue(GaussPoint *gp, TimeStep *tStep)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: computePseudoInverse(FloatMatrix &iDdd, FloatMatrix &b)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI4dF_num(FloatArray &answer, const FloatMatrix &F)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial ::  compute_dI7dF_num(FloatArray &answer,const FloatMatrix &F, const FloatArray &D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI10dF_num(FloatArray &answer,const FloatMatrix &F, const FloatArray &D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI7dD_num(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dI10dD_num(FloatArray &answer, const FloatMatrix &F, const FloatArray &D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I4dF2_num(FloatMatrix &answer,const FloatMatrix &F)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I7dF2_num(FloatMatrix &answer,const FloatMatrix &F, const FloatArray&D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I10dF2_num(FloatMatrix &answer,const FloatMatrix &F, const FloatArray&D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial ::  compute_d2I7dFdD_num(FloatMatrix &answer, const FloatMatrix &F, const FloatArray&D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial ::  compute_d2I10dFdD_num(FloatMatrix &answer, const FloatMatrix &F, const FloatArray&D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_d2I7dD2_num(FloatMatrix &answer, const FloatMatrix &F, const FloatArray&D)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dPdF_dEdF(FloatMatrix &dPdF,FloatMatrix &dEdF, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep)
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
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dPdF_dEdF_fake(FloatMatrix &dPdF,FloatMatrix &dEdF, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep)
{

  double em = this->epsilon_m;
  this->epsilon_m = 0;
  this->epsilon_f = 5.e-12;
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

  this->epsilon_m = em;
  this->epsilon_f = 0;
  
}




void
MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial :: compute_dPdD_dEdD(FloatMatrix &dPdD,FloatMatrix &dEdD, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep)
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
