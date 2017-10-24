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

#include "misesmatfinitestrain.h"
#include "Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(MisesMatFiniteStrain);

// constructor
  MisesMatFiniteStrain :: MisesMatFiniteStrain(int n, Domain *d) : StructuralMaterial(n, d), IncompressibleMaterialExtensionInterface()
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    H = 0.;
    sig0 = 0.;
    G = 0.;
    K = 0.;
}

// destructor
MisesMatFiniteStrain :: ~MisesMatFiniteStrain()
{
    delete linearElasticMaterial;
}

// reads the model parameters from the input file
IRResultType
MisesMatFiniteStrain :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    /*      result = linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    if ( result != IRRT_OK ) return result;
    */


    G = static_cast< IsotropicLinearElasticMaterial * >(linearElasticMaterial)->giveShearModulus();
    K = static_cast< IsotropicLinearElasticMaterial * >(linearElasticMaterial)->giveBulkModulus();
    
    IR_GIVE_FIELD(ir, K, _IFT_MisesMatFiniteStrain_K); // uniaxial yield stress
    IR_GIVE_FIELD(ir, mu, _IFT_MisesMatFiniteStrain_mu); // uniaxial yield stress

    IR_GIVE_FIELD(ir, sig0, _IFT_MisesMatFiniteStrain_sig0); // uniaxial yield stress

    H = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, H, _IFT_MisesMatFiniteStrain_h); // hardening modulus


    return IRRT_OK;
}

// creates a new material status  corresponding to this class
MaterialStatus *
MisesMatFiniteStrain :: CreateStatus(GaussPoint *gp) const
{
    return new MisesMatFiniteStrainStatus(1, this->giveDomain(), gp);
}



void
MisesMatFiniteStrain :: giveCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
{

 MisesMatFiniteStrainStatus *status = static_cast< MisesMatFiniteStrainStatus * >( this->giveStatus(gp) );
 
  FloatArray devC;

  this->giveDeviatoricCauchyStressVector_3d(devC, gp, vF, tStep);

  FloatArray tempvF = status->giveTempFVector();
  FloatMatrix F;
  F.beMatrixForm(tempvF);
  double J  = F.giveDeterminant();
  this->giveVolumetricCauchyStressVector_3d(answer, gp, J);
  
  answer.add(devC);

  status->letTempCVectorBe(answer);
  

}




void
MisesMatFiniteStrain :: giveDeviatoricCauchyStressVector_3d(FloatArray &answer,
                                       GaussPoint *gp,
                                       const FloatArray &redvF,
                                       TimeStep *tStep)
{
    MisesMatFiniteStrainStatus *status = static_cast< MisesMatFiniteStrainStatus * >( this->giveStatus(gp) );

    double J, kappa;
    FloatArray vF, vInvCp, vTempInvCp, vFn, vtempFn;
    FloatMatrix f, Fn, F, invCp, be_trial, junk;
    // initilizing data from status
    kappa = status->giveCumulativePlasticStrain();
    vInvCp = status->giveInvLeftCauchyGreenPl();
    vFn = status->giveFVector();
    Fn.beMatrixForm(vFn);
    StructuralMaterial :: giveFullVectorFormF( vF, redvF, gp->giveMaterialMode() );
    f.beMatrixForm(vF);
    F.beProductOf(f,Fn);
    // start computation
    invCp.beMatrixForm(vInvCp);
    junk.beProductTOf(invCp,F);
    J = F.giveDeterminant();
    junk.beProductTOf(invCp, F);
    // trial elastic left Cauchy-Green
    be_trial.beProductOf(F,junk);
    // principal value decomposition
    FloatArray eVals, lambdae_tr(3);
    FloatMatrix eVecs;
    be_trial.jaco_(eVals, eVecs, 8);
    lambdae_tr.at(1) = sqrt(eVals.at(1));
    lambdae_tr.at(2) = sqrt(eVals.at(2));
    lambdae_tr.at(3) = sqrt(eVals.at(3));

    // eigenvalues of Kirchhoff stress
    FloatArray evalsTau_tr(3);
    evalsTau_tr.at(1) = 2*mu*log(lambdae_tr.at(1)) - (2*mu/3)*log(J);
    evalsTau_tr.at(2) = 2*mu*log(lambdae_tr.at(2)) - (2*mu/3)*log(J);
    evalsTau_tr.at(3) = 2*mu*log(lambdae_tr.at(3)) - (2*mu/3)*log(J);

    FloatArray vTau_tr;
    FloatMatrix tau_tr(3, 3);
    for ( int i = 1; i < 4; i++ ) {
      for ( int j = 1; j < 4; j++ ) {
	tau_tr.at(i, j) = evalsTau_tr.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + evalsTau_tr.at(2)*eVecs.at(i, 2) * eVecs.at(j, 2) + evalsTau_tr.at(3)*eVecs.at(i, 3) * eVecs.at(j, 3);
      }
    }

    vTau_tr.beSymVectorForm(tau_tr);

    double dKappa, yieldValue;
    double trialS = computeStressNorm(vTau_tr);
    double sigmaY = sig0 + H * kappa;
    yieldValue = sqrt(3./2.)*trialS-sigmaY;
    FloatArray N(3);
    N.zero();
    //store deviatoric trial stress(reused by algorithmic stiffness)
    status->letTrialStressDevBe(vTau_tr);
    status->letTrialLambdaBe(lambdae_tr);
    //the radial return-mapping algorithm
    if ( yieldValue > 0 ) {
        dKappa = yieldValue / ( 3.*mu + H );
        kappa = kappa + dKappa;
	N =  evalsTau_tr;
	double factor = sqrt(2./3.) * trialS;
        N.times(1/factor);
	evalsTau_tr.times(1-2*mu*dKappa/(sqrt(2./3.)*trialS)); 
    }



    FloatArray lambda(3);
    lambda.at(1) = exp(log(lambdae_tr.at(1)) - dKappa*N.at(1));
    lambda.at(2) = exp(log(lambdae_tr.at(2)) - dKappa*N.at(2));
    lambda.at(3) = exp(log(lambdae_tr.at(3)) - dKappa*N.at(3));


    FloatMatrix tau(3,3), be(3,3);

  for ( int i = 1; i < 4; i++ ) {
      for ( int j = 1; j < 4; j++ ) {
	tau.at(i, j) = evalsTau_tr.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + evalsTau_tr.at(2)*eVecs.at(i, 2) * eVecs.at(j, 2) + evalsTau_tr.at(3)*eVecs.at(i, 3) * eVecs.at(j, 3);
	be.at(i, j) = lambda.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + lambda.at(2)*eVecs.at(i, 2) * eVecs.at(j, 2) + lambda.at(3)*eVecs.at(i, 3) * eVecs.at(j, 3);

      }
    }

  FloatMatrix invF, tempInvCp;
  invF.beInverseOf(F);
  junk.beProductTOf(be,F);
  tempInvCp.beProductOf(F,junk);
  vTempInvCp.beVectorForm(tempInvCp);

  answer.beSymVectorForm(tau);
  answer.times(1./J);
    
  // store trial values used for tangent stiffness
  

  status->letTempInvLeftCauchyGreenPlBe(vTempInvCp);
  status->setTempCumulativePlasticStrain(kappa);    

  status->letTempDeviatoricCauchyStressVectorBe(answer);
  status->letTempFVectorBe(vF);
    
}



void
MisesMatFiniteStrain :: giveVolumetricCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, double J)
{
  MisesMatFiniteStrainStatus *status = static_cast< MisesMatFiniteStrainStatus * >( this->giveStatus(gp) );

    double pressure = this->computePressure(J);
    answer.resize(6);
    answer.at(1) = answer.at(2) = answer.at(3) = 1.;
    answer.times(pressure);

    status->letTempVolumetricCauchyStressVectorBe(answer);
  

}


double
MisesMatFiniteStrain :: computePressure(double J)
{
    return K*log(J)/J;
}



void 
MisesMatFiniteStrain :: giveInitialStiffnessMatrix_Cauchy(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MisesMatFiniteStrainStatus *status = static_cast< MisesMatFiniteStrainStatus * >( this->giveStatus(gp) );
  FloatArray volC = status->giveTempVolumetricCauchyStressVector();
  FloatArray devC = status->giveTempDeviatoricCauchyStressVector();
  FloatArray redvCauchy, vCauchy;

  redvCauchy = volC;
  redvCauchy.add(devC);

  FloatMatrix stress_ident(9,9);

  if ( redvCauchy.giveSize() ) {
    StructuralMaterial :: giveFullSymVectorForm( vCauchy, redvCauchy, gp->giveMaterialMode() );
    // product Cauchy_il delta_jk
     /*
       [ sig11,     0,     0,     0, sig13, sig12,     0,     0,     0]
       [     0, sig22,     0, sig23,     0,     0,     0,     0, sig21]
       [     0,     0, sig33,     0,     0,     0, sig32, sig31,     0]
       [     0, sig32,     0, sig33,     0,     0,     0,     0, sig31]
       [ sig31,     0,     0,     0, sig33, sig32,     0,     0,     0]
       [ sig21,     0,     0,     0, sig23, sig22,     0,     0,     0]
       [     0,     0, sig23,     0,     0,     0, sig22, sig21,     0]
       [     0,     0, sig13,     0,     0,     0, sig12, sig11,     0]
       [     0, sig12,     0, sig13,     0,     0,     0,     0, sig11]
     */
    stress_ident.at(1, 1) = vCauchy.at(1);
    stress_ident.at(1, 5) = vCauchy.at(5);
    stress_ident.at(1, 6) = vCauchy.at(6);
    
    stress_ident.at(2, 2) = vCauchy.at(2);
    stress_ident.at(2, 4) = vCauchy.at(4);
    stress_ident.at(2, 9) = vCauchy.at(6);
    
    stress_ident.at(3, 3) = vCauchy.at(3);
    stress_ident.at(3, 7) = vCauchy.at(4);
    stress_ident.at(3, 8) = vCauchy.at(5);
    
    stress_ident.at(4, 2) = vCauchy.at(4);
    stress_ident.at(4, 4) = vCauchy.at(3);
    stress_ident.at(4, 9) = vCauchy.at(5);
    
    stress_ident.at(5, 1) = vCauchy.at(5);
    stress_ident.at(5, 5) = vCauchy.at(3);
    stress_ident.at(5, 6) = vCauchy.at(4);
    
    stress_ident.at(6, 1) = vCauchy.at(6);
    stress_ident.at(6, 5) = vCauchy.at(4);
    stress_ident.at(6, 6) = vCauchy.at(2);
    
    stress_ident.at(7, 3) = vCauchy.at(4);
    stress_ident.at(7, 7) = vCauchy.at(2);
    stress_ident.at(7, 8) = vCauchy.at(6);
    
    stress_ident.at(8, 3) = vCauchy.at(5);
    stress_ident.at(8, 7) = vCauchy.at(6);
    stress_ident.at(8, 8) = vCauchy.at(1);
    
    stress_ident.at(9, 2) = vCauchy.at(6);
    stress_ident.at(9, 4) = vCauchy.at(5);
    stress_ident.at(9, 9) = vCauchy.at(1);
  } else {
    stress_ident.zero();
  }
  IntArray indx;
  StructuralMaterial :: giveVoigtVectorMask(indx,gp->giveMaterialMode());
  answer.beSubMatrixOf(stress_ident, indx, indx);
  answer.symmetrized();
  
}
  



void MisesMatFiniteStrain :: computeCumPlastStrain(double &tempKappa, GaussPoint *gp, TimeStep *tStep)
{
    MisesMatFiniteStrainStatus *status = static_cast< MisesMatFiniteStrainStatus * >( this->giveStatus(gp) );
    tempKappa = status->giveTempCumulativePlasticStrain();
}


void 
MisesMatFiniteStrain :: give3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  FloatMatrix devC, pC;
  MisesMatFiniteStrainStatus *status = static_cast< MisesMatFiniteStrainStatus * >( this->giveStatus(gp) );
  FloatArray vF = status->giveTempFVector();
  FloatMatrix F;
  F.beMatrixForm(vF);
  double J  = F.giveDeterminant();



  this->giveDeviatoric3dMaterialStiffnessMatrix_dCde(devC, mode, gp, tStep);
  this->givePressure3dMaterialStiffnessMatrix_dCde(pC, J);
  this->giveVolumetric3dMaterialStiffnessMatrix_dCde(answer, J);

  answer.add(devC);
  answer.add(pC);
}

void 
MisesMatFiniteStrain :: giveDeviatoric3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  
  answer.resize(6,6);
  answer.zero();  

  MisesMatFiniteStrainStatus *status = static_cast< MisesMatFiniteStrainStatus * >( this->giveStatus(gp) );
  FloatArray vF = status->giveTempFVector();
  FloatMatrix F;
  F.beMatrixForm(vF);
  double J  = F.giveDeterminant();
  /*
   * Trial values from status
   */
  double trialS = 0, factor, d1, d2, d3, dKappa;
  FloatArray trialStressDev, vCauchyDev, eVals_tr, eVals, n, lambda_tr;
  FloatMatrix tau, eVecs_tr, eVecs, I(3,3), Ones(3,3), N, T, CauchyDev;
  Ones.beOnesMatrix();
  
  lambda_tr = status->giveTrialLambda();
  trialStressDev = status->giveTrialStressDev();
  vCauchyDev = status->giveTempDeviatoricCauchyStressVector();
  CauchyDev.beMatrixForm(vCauchyDev);
  CauchyDev.jaco_(eVals, eVecs, 8);
  

  if(trialStressDev.giveSize()) {
    tau.beMatrixForm(trialStressDev);
    trialS = computeStressNorm(trialStressDev);
  }  else {
    tau.resize(3,3);
    tau.zero();
  }


  tau.jaco_(eVals_tr, eVecs_tr, 8);
  dKappa = status->giveTempCumulativePlasticStrain()-status->giveCumulativePlasticStrain();
  
  n = eVals_tr;
  factor = sqrt(2./3.) * trialS;
  n.times(1./factor);
  
  I.beUnitMatrix();
  Ones.beOnesMatrix();
  
  FloatMatrix C1(3,3);
  C1.zero();

  if( dKappa > 0 ) {
    d1 = (1.-2.*mu*dKappa/factor)*2.*mu;
    d2 = (1.-2.*mu*dKappa/factor)*(-2./3.*mu);
    d3 =  4.*mu*mu*(sqrt(2./3.)*dKappa/trialS -(1./(3.*mu+H)));
    I.times(d1);
    Ones.times(d2);
    N.beDyadicProductOf(n,n);    
    N.times(d3);  
    C1.add(I);
    C1.add(Ones);
    C1.add(N);
    C1.times(1./J);
  } else {
    d1 = 2*mu;
    d2 = - 2./3.*mu;
    d3 = 0;
    I.times(d1);
    Ones.times(d2);
    C1.add(I);
    C1.add(Ones);
    C1.times(1./J);
  }


  answer.resize(6,6);
  answer.zero();


  for( int l = 1; l <= 3; l++ ) {
    for(int k = 1; k <= 3; k++ ) {
      for ( int j = 1; j <= 3; j++) {
	for ( int i = 1; i <= 3; i++) {
	  if(i >= j && k >= l) {
	    double sum = 0;
	    for ( int alpha = 1; alpha <= 3; alpha++) {
	      sum = sum - 2*eVals.at(alpha)*eVecs.at(i,alpha)*eVecs.at(j,alpha)* eVecs.at(k,alpha)*eVecs.at(l,alpha);
	      for ( int beta = 1; beta <= 3; beta ++) { 
		sum = sum + (C1.at(alpha,beta))*(eVecs.at(i,alpha)*eVecs.at(j,alpha)*eVecs.at(k,beta)*eVecs.at(l,beta));
		double lambda_a = lambda_tr.at(alpha);
		double lambda_b = lambda_tr.at(beta);
		double sigma_a  = eVals.at(alpha);
		double sigma_b  = eVals.at(beta);
		if  (alpha != beta) {
		  sum = sum + this->giveSl(lambda_a,lambda_b,sigma_a,sigma_b,J)*(eVecs.at(i,alpha)*eVecs.at(j,beta)*(eVecs.at(k,alpha)*eVecs.at(l,beta)+eVecs.at(k,beta)*eVecs.at(l,alpha)));                        
		}
	      }
	    }
	    answer.at(giveSymVI(i, j), giveSymVI(k, l)) += sum;
	  }
	}
      }
    }    
  }

}




double 
MisesMatFiniteStrain :: giveSl(double lambda_alpha,double lambda_beta, double sigma_alpha, double sigma_beta,double J)
{ 
  if (fabs(lambda_alpha-lambda_beta)<1.e-5) {
    return (mu/J - sigma_alpha);     
  } else
    return (sigma_alpha*lambda_beta*lambda_beta-sigma_beta*lambda_alpha*lambda_alpha)/(lambda_alpha*lambda_alpha-lambda_beta*lambda_beta); 
}


void 
MisesMatFiniteStrain :: givePressure3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, double J)
{
    
    double pressure = this->computePressure(J);
    
    
    answer.resize(6,6);
    answer.zero();
    
    answer.at(1,1) = answer.at(2,2) = answer.at(3,3) = -1;
    answer.at(4,4) = answer.at(5,5) = answer.at(6,6) = -1;
    
    answer.at(1,2) = answer.at(1,3) = answer.at(2,3) = 1;
    answer.at(2,1) = answer.at(3,1) = answer.at(3,2) = 1;
        
    answer.times(pressure);
}

void 
MisesMatFiniteStrain :: giveVolumetric3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, double J)
{

    answer.resize(6,6);
    answer.at(1,1) = answer.at(2,2) = answer.at(3,3) = 1;
    answer.at(1,2) = answer.at(1,3) = answer.at(2,3) = 1;
    answer.at(2,1) = answer.at(3,1) = answer.at(3,2) = 1;
        
    double factor = K/(J) -  K*log(J)/J;
    answer.times(factor);
}




int
MisesMatFiniteStrain :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    MisesMatFiniteStrainStatus *status = static_cast< MisesMatFiniteStrainStatus * >( this->giveStatus(gp) );
    /*    if ( type == IST_PlasticStrainTensor ) {
        answer = status->givePlasDef();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveCumulativePlasticStrain();
        return 1;
    } else if ( ( type == IST_DamageScalar ) || ( type == IST_DamageTensor ) ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
	}*/
}




MisesMatFiniteStrainStatus :: MisesMatFiniteStrainStatus(int n, Domain *d, GaussPoint *g) :
  StructuralMaterialStatus(n, d, g), trialStressDev(), trialLambda(), devCauchyStressVector(6), volCauchyStressVector(6)
{
   
    kappa = tempKappa = 0.;
    invCp.resize(6);
    invCp.at(1) = invCp.at(2) = invCp.at(3) = 1;
    trialLambda.resize(3);
    trialLambda.at(1) = trialLambda.at(2) = trialLambda.at(3) = 1;

    tempDevCauchyStressVector = devCauchyStressVector;
    tempVolCauchyStressVector = volCauchyStressVector;

    
}

MisesMatFiniteStrainStatus :: ~MisesMatFiniteStrainStatus()
{ }

void
MisesMatFiniteStrainStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);

}


// initializes temporary variables based on their values at the previous equlibrium state
void MisesMatFiniteStrainStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    tempDevCauchyStressVector = devCauchyStressVector;
    tempVolCauchyStressVector = volCauchyStressVector;
    tempInvCp = invCp;
    tempKappa = kappa;
    //    trialStressDev.clear(); // to indicate that it is not defined yet
}


// updates internal variables when equilibrium is reached
void
MisesMatFiniteStrainStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);

    kappa = tempKappa;
    trialStressDev.clear(); // to indicate that it is not defined any more
    invCp = tempInvCp;
    devCauchyStressVector = tempDevCauchyStressVector;
    volCauchyStressVector = tempVolCauchyStressVector;

    
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
MisesMatFiniteStrainStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data


    // write cumulative plastic strain (scalar)
    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    return CIO_OK;
}



contextIOResultType
MisesMatFiniteStrainStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    // read cumulative plastic strain (scalar)
    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
