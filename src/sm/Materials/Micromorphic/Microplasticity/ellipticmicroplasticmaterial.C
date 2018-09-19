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

#include "ellipticmicroplasticmaterial.h"
#include "stressvector.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"
#include "Materials/isolinearelasticmaterial.h"


namespace oofem {
  REGISTER_Material(EllipticMicroplasticMaterial);

  EllipticMicroplasticMaterial :: EllipticMicroplasticMaterial(int n, Domain *d) : EllipticPlasticMaterial(n, d),MicromorphicMaterialExtensionInterface(d)//,IsotropicLinearElasticMaterial(n, d)
{
  //linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);

  Hk = Ak = 0.;
}

EllipticMicroplasticMaterial :: ~EllipticMicroplasticMaterial()
{
  //delete linearElasticMaterial;
  
}

// creates a new material status  corresponding to this class
MaterialStatus *
EllipticMicroplasticMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new EllipticMicroplasticMaterialStatus(1, this->giveDomain(), gp);
}

void
EllipticMicroplasticMaterial :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}




void
EllipticMicroplasticMaterial :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &micromorphicVar, const FloatArray &micromorphicVarGrad, TimeStep *tStep)
{
    double tempKappa, toSolveScalar;
    FloatArray tempPlasticStrain, tempStress, trialStress, tempTensor2, toSolveTensor, PS, N, incStress;
    FloatMatrix De, Ce, PP, DalgoTensor;

    EllipticMicroplasticMaterialStatus *status = static_cast< EllipticMicroplasticMaterialStatus * >( this->giveStatus(gp) );
    FloatMatrix invP(6,6);
    invP.at(1,1) = invP.at(2,2) = invP.at(3,3) = 1;
    invP.at(4,4) = invP.at(5,5) = invP.at(6,6) = 0.5;
    // initialize the plastic strain and cumulative plastic strain
    // by values after the previous step
    tempPlasticStrain = status->givePlasticStrain();
    tempKappa = status->giveCumulativePlasticStrain();
    // evaluate the trial stress
    this->linearElasticMaterial->give3dMaterialStiffnessMatrix(De, TangentStiffness, gp, tStep);
    Ce.beInverseOf(De);
    FloatArray elStrain;
    if(gp->giveMaterialMode() == _PlaneStrain) {
      StructuralMaterial :: giveFullSymVectorForm(elStrain, totalStrain, _PlaneStrain);                                              
    }
    elStrain.subtract(tempPlasticStrain);
    trialStress.beProductOf(De, elStrain);
    this->constructYieldFunctionTensor(PP);
    // apply the iterative procedure that solves the system of nonlinear equations
    // consisting of the yield condition and discretized flow rule
    // and evaluates:
    // tempKappa   ... cumulative plastic strain at the end of the substep
    // tempStress  ... stress at the end of the substep
    // tempPlasticStrain ... plastic strain at the end of the substep
    tempTensor2.beProductOf(PP, trialStress);
    double SPS = sqrt( trialStress.dotProduct(tempTensor2) );
    double yieldValue = SPS - this->evaluateCurrentYieldStress(tempKappa, micromorphicVar.at(1));
    tempStress = trialStress;
    double rel_yield_tol = 1.e-6;
    if ( yieldValue < rel_yield_tol ) {
      // trial stress in elastic domain
    } else {
        // return to the yield surface needed
        // Initial valuesr
        toSolveTensor.resize(6);
        toSolveTensor.zero();
        toSolveScalar = yieldValue;
        double errorF = yieldValue;
        double errorR = 0;
        DalgoTensor = De;
	PS.beProductOf(PP, tempStress);
	SPS = sqrt(PS.dotProduct(tempStress));
	FloatArray dFdS;
	dFdS = PS;
	dFdS.times(1./sqrt(tempStress.dotProduct(PS)));
        double deltaKappa = 0.;
	bool convergence = false;
	double strain_tol = 1.e-6;
	int flagLoop = 0;
	int max_num_iter = 100;
        do {
  	    double plasModulus = evaluateCurrentPlasticModulus(tempKappa + deltaKappa);
            //*************************************
            //Evaluation of the Recursive Equations
            //*************************************
            tempTensor2.beProductOf(DalgoTensor, dFdS);
            double beta = dFdS.dotProduct(tempTensor2);
            beta +=  plasModulus;
            // Construction of the equation of Delta Kappa
            tempTensor2.beProductOf(DalgoTensor, toSolveTensor);
            double tempScalar = dFdS.dotProduct(tempTensor2);
            double incKappa =  ( toSolveScalar - tempScalar ) / beta;
            tempTensor2 = dFdS;
            tempTensor2.times(incKappa);
            tempTensor2 += toSolveTensor;
            incStress.beProductOf(DalgoTensor, tempTensor2);
	    ////////////////////////////////////////////////////////
	    deltaKappa += incKappa;
	    tempKappa += incKappa;
	    tempStress -= incStress;
	    //*************************************
	    // Evaluation of the f and R
	    //*************************************
	    // Construction of the derivative of the plastic flow
	    PS.beProductOf(PP, tempStress);
	    dFdS = PS;
	    dFdS.times(1./sqrt(tempStress.dotProduct(PS)));
	    //////////////////////////////
	    FloatMatrix d2FdS2, temp, tempTensor4;
	     d2FdS2 = PP;
	     temp.beDyadicProductOf(PS, PS);
	     temp.times(1./(tempStress.dotProduct(PS)));
	     d2FdS2.subtract(temp);
	     d2FdS2.times(1./sqrt(tempStress.dotProduct(PS)));
	    // Construction of the gradient Nabla_S of R and SSa tensor
	    tempTensor4 = d2FdS2;
	    tempTensor4.times(deltaKappa);
	    tempTensor4.add(Ce);
	    DalgoTensor.beInverseOf(tempTensor4);
	    // Evaluation of R
	    tempTensor2.beProductOf( Ce, ( tempStress - trialStress ) );
	    toSolveTensor = tempTensor2 + deltaKappa * dFdS;
	    // Evaluation of f
	    tempTensor2.beProductOf(PP, tempStress);
	    SPS = sqrt( tempStress.dotProduct(tempTensor2) );
	    toSolveScalar = SPS - this->evaluateCurrentYieldStress(tempKappa, micromorphicVar.at(1));
	    //*************************************
	    // Evaluation of the error
	    //*************************************
	    FloatArray errLoop;
	    errLoop = toSolveTensor;
	    errorR = sqrt( errLoop.dotProduct(errLoop) );
	    errorF = fabs(toSolveScalar);
	    
	    flagLoop++;
            convergence = ( fabs(errorF) < rel_yield_tol && errorR < strain_tol );
        } while ( flagLoop <= max_num_iter && !convergence );
	tempPlasticStrain += deltaKappa * dFdS;
    }
    FloatArray redStress;
    if(gp->giveMaterialMode() == _PlaneStrain) {
      StructuralMaterial :: giveReducedSymVectorForm(redStress, tempStress, _PlaneStrain);                                         
    }
    status->letTempPlasticStrainBe(tempPlasticStrain);
    status->setTempCumulativePlasticStrain(tempKappa);
    status->letTempStressVectorBe(redStress);
 
}
  


double
EllipticMicroplasticMaterial :: evaluateCurrentYieldStress(double kappa, double microKappa)
{
  double Str = Hk * (microKappa - kappa);
  return sig0 + H * kappa - Str;
}

double
EllipticMicroplasticMaterial :: evaluateCurrentPlasticModulus(const double kappa)
{
  return H + Hk;
}

  
void
EllipticMicroplasticMaterial :: giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    double tempKappa, kappa;
    FloatArray tempTensor2, N;
    EllipticMicroplasticMaterialStatus *status = static_cast< EllipticMicroplasticMaterialStatus * >( this->giveStatus(gp) );
    
    FloatArray PS, tempStress, DaN;
    FloatMatrix dN, PP, tempTensor4, Dalgo, cor, De, Ce;
    this->linearElasticMaterial->give3dMaterialStiffnessMatrix(De, TangentStiffness, gp, tStep);
    if ( mode == ElasticStiffness ) {
      answer = De;
    } else if ( mode == TangentStiffness ) {
        kappa = status->giveCumulativePlasticStrain();
        tempKappa = status->giveTempCumulativePlasticStrain();
	double dKappa = tempKappa - kappa;
        if ( dKappa > 1.e-13 ) {
	  //dN_dSig
	  FloatArray redTempStress;
	  redTempStress = status->giveTempStressVector();
	  this->constructYieldFunctionTensor(PP);
	  // Construction of the derivative of the plastic flow
	  StructuralMaterial :: giveFullSymVectorForm(tempStress, redTempStress, _PlaneStrain); 
 
	  PS.beProductOf(PP, tempStress);
	  FloatArray dFdS;
	  dFdS = PS;
	  dFdS.times(1./sqrt(tempStress.dotProduct(PS)));
	  //////////////////////////////
	  FloatMatrix d2FdS2, temp, tempTensor4;
	  d2FdS2 = PP;
	  temp.beDyadicProductOf(PS, PS);
	  temp.times(1./(tempStress.dotProduct(PS)));
	  d2FdS2.subtract(temp);
	  d2FdS2.times(1./sqrt(tempStress.dotProduct(PS)));
	  // Construction of the gradient Nabla_S of R and SSa tensor
	  Ce.beInverseOf(De);
	  tempTensor4 = d2FdS2;
	  tempTensor4.times(dKappa);
	  tempTensor4.add(Ce);
	  Dalgo.beInverseOf(tempTensor4);
	  double norm = sqrt(tempStress.dotProduct(PS));
	  /// N
	  N = PS;
	  N.times(1./norm);
	  DaN.beProductOf(Dalgo,N);
	  cor.beDyadicProductOf(DaN, DaN);
	  double den = N.dotProduct(DaN);
	  double Hmod = evaluateCurrentPlasticModulus(tempKappa);
	  den += Hmod;
	  cor.times( 1./ den);
	  answer = Dalgo;
	  answer.subtract(cor);	  
        } else {
	  answer = De;
	}
    }

    if(gp->giveMaterialMode() == _PlaneStrain) {
      FloatMatrix plS;
      plS.resize(4, 4);
      plS.zero();
      //answer.beSubMatrixOf(m3d, indx, indx);
      
      plS.at(1, 1) = answer.at(1, 1);
      plS.at(1, 2) = answer.at(1, 2);
      plS.at(1, 3) = answer.at(1, 3);
      plS.at(1, 4) = answer.at(1, 6);
	
      plS.at(2, 1) = answer.at(2, 1);
      plS.at(2, 2) = answer.at(2, 2);
      plS.at(2, 3) = answer.at(2, 3);
      plS.at(2, 4) = answer.at(2, 6);
	
      plS.at(3, 1) = answer.at(3, 1);
      plS.at(3, 2) = answer.at(3, 2);
      plS.at(3, 3) = answer.at(3, 3);
      plS.at(3, 4) = answer.at(3, 6);
	
      plS.at(4, 1) = answer.at(6, 1);
      plS.at(4, 2) = answer.at(6, 2);
      plS.at(4, 3) = answer.at(6, 3);
      plS.at(4, 4) = answer.at(6, 6);
	
      answer = plS;
      }



    
}

void
EllipticMicroplasticMaterial :: giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  EllipticMicroplasticMaterialStatus *status = static_cast< EllipticMicroplasticMaterialStatus * >( this->giveStatus(gp) );
  MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _1dMat) {
    answer.resize(1,1);
  } else if(matMode == _PlaneStrain) {
    answer.resize(4,1);
  } else {
    answer.resize(6,1);
  }
  answer.zero();

  double kappa = status->giveCumulativePlasticStrain();
  // increment of cumulative plastic strain as an indicator of plastic loading
  double tempKappa = status->giveTempCumulativePlasticStrain();
  if ( tempKappa <= kappa || tempKappa < 1.e-15 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
    return;
  }
  double dKappa = tempKappa - kappa;

  FloatArray redTempStress = status->giveTempStressVector();
  //////////////////////////////
  FloatArray PS, N, DaN, tempStress;
  FloatMatrix d2FdS2, temp, tempTensor4, PP, De, Ce, Dalgo;

  this->linearElasticMaterial->give3dMaterialStiffnessMatrix(De, TangentStiffness, gp, tStep);
  Ce.beInverseOf(De);
  
  this->constructYieldFunctionTensor(PP);
  StructuralMaterial :: giveFullSymVectorForm(tempStress, redTempStress, _PlaneStrain); 
  PS.beProductOf(PP, tempStress);

  d2FdS2 = PP;
  temp.beDyadicProductOf(PS, PS);
  temp.times(1./(tempStress.dotProduct(PS)));
  d2FdS2.subtract(temp);
  d2FdS2.times(1./sqrt(tempStress.dotProduct(PS)));
  // Construction of the gradient Nabla_S of R and SSa tensor
  Ce.beInverseOf(De);
  tempTensor4 = d2FdS2;
  tempTensor4.times(dKappa);
  tempTensor4.add(Ce);
  Dalgo.beInverseOf(tempTensor4);
  double norm = sqrt(tempStress.dotProduct(PS));
  /// N
  N = PS;
  N.times(1./norm);
  DaN.beProductOf(Dalgo,N);
  double den = N.dotProduct(DaN);
  double Hmod = evaluateCurrentPlasticModulus(tempKappa);
  den += Hmod;
  answer.initFromVector(DaN, false);
  answer.times( -Hk/ den);

  if(gp->giveMaterialMode() == _PlaneStrain) {
    FloatMatrix plS;
    plS.resize(4, 1);
    plS.zero();
    //answer.beSubMatrixOf(m3d, indx, indx);    
    plS.at(1, 1) = answer.at(1, 1);
    plS.at(2, 1) = answer.at(2, 1);
    plS.at(3, 1) = answer.at(3, 1);
    plS.at(4, 1) = answer.at(6, 1);
        
    answer = plS;
  }



  
}


void
EllipticMicroplasticMaterial :: giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  FloatMatrix m;
  this->giveMicromorphicMatrix_dSigdPhi(m, mode, gp, tStep);
  answer.beTranspositionOf(m);
  
}



void
EllipticMicroplasticMaterial :: giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  
  EllipticMicroplasticMaterialStatus *status = static_cast< EllipticMicroplasticMaterialStatus * >( this->giveStatus(gp) );
  
  double kappa = status->giveCumulativePlasticStrain();
  // increment of cumulative plastic strain as an indicator of plastic loading
  double tempKappa = status->giveTempCumulativePlasticStrain();
 
  answer.resize(1,1);
  answer.at(1,1) = Hk;
  
  if ( tempKappa <= kappa || tempKappa < 1.e-15 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
    return;
  }
  double dKappa = tempKappa - kappa;
  // Construction of the derivative of the plastic flow
  FloatArray tempPlasticStrain, tempStress, trialStress, tempTensor2, toSolveTensor, dFdS, PS, N, incStress, redTempStress;
  FloatMatrix De, Ce, PP, DalgoTensor;
  this->linearElasticMaterial->give3dMaterialStiffnessMatrix(De, TangentStiffness, gp, tStep);
  Ce.beInverseOf(De);
  
  redTempStress = status->giveTempStressVector();
  StructuralMaterial :: giveFullSymVectorForm(tempStress, redTempStress, _PlaneStrain);
  this->constructYieldFunctionTensor(PP);
  PS.beProductOf(PP, tempStress);
  dFdS = PS;
  dFdS.times(1./sqrt(tempStress.dotProduct(PS)));
  //////////////////////////////
  FloatMatrix d2FdS2, temp, tempTensor4;
  d2FdS2 = PP;
  temp.beDyadicProductOf(PS, PS);
  temp.times(1./(tempStress.dotProduct(PS)));
  d2FdS2.subtract(temp);
  d2FdS2.times(1./sqrt(tempStress.dotProduct(PS)));
  // Construction of the gradient Nabla_S of R and SSa tensor
  tempTensor4 = d2FdS2;
  tempTensor4.times(dKappa);
  tempTensor4.add(Ce);
  DalgoTensor.beInverseOf(tempTensor4);
  tempTensor2.beProductOf(DalgoTensor, dFdS);
  
  double plasModulus = evaluateCurrentPlasticModulus(tempKappa);  
  double beta = dFdS.dotProduct(tempTensor2);
  beta +=  plasModulus;

  answer.at(1,1) -= Hk*Hk/beta;

  

}


void
EllipticMicroplasticMaterial :: giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _1dMat) {
    answer.resize(1,1);
  } else if (matMode == _PlaneStrain) {
    answer.resize(2,2);
  } else {
    answer.resize(9,9);
  } 
  answer.beUnitMatrix();
  answer.times(Ak);

}




void
EllipticMicroplasticMaterial :: giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &S, FloatArray &M, GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{
  
    EllipticMicroplasticMaterialStatus *status = static_cast< EllipticMicroplasticMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    this->performPlasticityReturn(gp, totalStrain, micromorphicVar, micromorphicVarGrad, tStep);
    double kappa = status->giveTempCumulativePlasticStrain();
    S.resize(1);
    S.at(1) = Hk * (micromorphicVar.at(1) - kappa);
    M = Ak * micromorphicVarGrad;
        // store the plastic strain and cumulative plastic strain
    //    status->letTempPlasticStrainBe(plStrain);
    //status->setTempCumulativePlasticStrain(kappa);
    //status->setDKappa(dKappa);
    //store stresses into the status
    //status->letTempStressVectorBe(answer);
    status->letTempMicromorphicStressBe(S);
    status->letTempMicromorphicStressGradBe(M);    
    
    status->letTempStrainVectorBe(totalStrain);
    status->letTempMicromorphicVarBe(micromorphicVar);
    status->letTempMicromorphicVarGradBe(micromorphicVarGrad);
  
    
      
}




int
EllipticMicroplasticMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    EllipticMicroplasticMaterialStatus *status = static_cast< EllipticMicroplasticMaterialStatus * >( this->giveStatus(gp) );
    if  ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveCumulativePlasticStrain();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

IRResultType
EllipticMicroplasticMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    EllipticPlasticMaterial :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, Hk, _IFT_MicromorphicMaterialExtensionInterface_Hk);
    IR_GIVE_FIELD(ir, Ak, _IFT_MicromorphicMaterialExtensionInterface_Ak);

    return IRRT_OK;
}


//=============================================================================

EllipticMicroplasticMaterialStatus :: EllipticMicroplasticMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MicromorphicMaterialStatus(n, d, g), plasticStrain(6), tempPlasticStrain()
{
    stressVector.resize(6);
    strainVector.resize(6);
    kappa = tempKappa = 0.;
    dKappa = 0.;
   
}

EllipticMicroplasticMaterialStatus :: ~EllipticMicroplasticMaterialStatus()
{ }



// initializes temporary variables based on their values at the previous equlibrium state
void EllipticMicroplasticMaterialStatus :: initTempStatus()
{
    MicromorphicMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
      plasticStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
      plasticStrain.zero();
    }
    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
 }


// updates internal variables when equilibrium is reached
void
EllipticMicroplasticMaterialStatus :: updateYourself(TimeStep *tStep)
{
    MicromorphicMaterialStatus :: updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    kappa = tempKappa;
 
}



} // end namespace oofem
