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

#include "camclaymat.h"
#include "isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "function.h"
#include "intarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "fieldmanager.h"
#include "Elements/structuralelement.h"
#include "engngm.h"

namespace oofem {

REGISTER_Material(CamClayMat);


CamClayMat :: CamClayMat(int n, Domain *d) : StructuralMaterial(n, d),
    linearElasticMaterial(n, d)
{
}


IRResultType
CamClayMat :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    result = linearElasticMaterial.initializeFrom(ir); // takes care of elastic constants
    if ( result != IRRT_OK ) return result;

    G = linearElasticMaterial.giveShearModulus();
    K = linearElasticMaterial.giveBulkModulus();

    // parameter describing the ration of the minor to major axis of the yield function
    IR_GIVE_FIELD(ir, M, _IFT_CamClayMat_M);   
	// initial preconsolidation pressure
    IR_GIVE_FIELD(ir, pc0, _IFT_CamClayMat_pc); 
	

	e = 0.2;
	IR_GIVE_OPTIONAL_FIELD(ir, e, _IFT_CamClayMat_e);

	consolidationIndex = 0.066;
	IR_GIVE_OPTIONAL_FIELD(ir, consolidationIndex, _IFT_CamClayMat_consolidationIndex);

	swellingIndex = 0.0077;
	IR_GIVE_OPTIONAL_FIELD(ir, swellingIndex, _IFT_CamClayMat_swellingIndex);

	//the minus included here in the definition of theta is consistent with the idea
	//that a negative change of volumetric deformation should lead to a positive change
	//in preconsolidation pressure and vice versa; incorporating this negativity here eases
	//the handling of hardening law with respect to its various derivatives used in the algorithm
	theta =  -(1 + e) / (consolidationIndex - swellingIndex);
	
    yieldTolerance = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTolerance, _IFT_CamClayMat_yieldTol);
    plasticStrainTolerance = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, plasticStrainTolerance, _IFT_CamClayMat_plStrainTol);
	pcTolerance = 1.e-6;
	IR_GIVE_OPTIONAL_FIELD(ir, pcTolerance, _IFT_CamClayMat_pcTol);
    maxIter = 20;
    IR_GIVE_OPTIONAL_FIELD(ir, maxIter, _IFT_CamClayMat_maxIter);   
    return IRRT_OK;
}

// creates a new material status  corresponding to this class
MaterialStatus *
CamClayMat :: CreateStatus(GaussPoint *gp) const
{
  return new CamClayMatStatus(1, this->giveDomain(), gp, pc0);
}

double CamClayMat::giveYieldValueAtStress(FloatArray &stress, double &pc)
{
	if (stress.giveSize() != 6) {
		OOFEM_ERROR("wrong size of stress vector")
	}

	double sigmaM = (1. / 3.)*(stress.at(1) + stress.at(2) + stress.at(3));
	double p = -sigmaM;

	FloatArray devStress = stress;
	devStress.at(1) -= sigmaM;
	devStress.at(2) -= sigmaM;
	devStress.at(3) -= sigmaM;

	double q = sqrt(3. / 2.) * computeStressNorm(devStress);
	
	return (q*q - M*M*p*(pc - p));
}

void
CamClayMat :: giveRealStressVector_3d(FloatArray &answer,
                                 GaussPoint *gp,
                                 const FloatArray &totalStrain,
                                 TimeStep *tStep)
{
    CamClayMatStatus *status = static_cast< CamClayMatStatus * >( this->giveStatus(gp) );
    // subtract stress independent part
    FloatArray strainR(6);
    this->giveStressDependentPartOfStrainVector(strainR, gp, totalStrain,
                                                tStep, VM_Total);

    this->performPlasticityReturn(gp, strainR, tStep);
    answer = status->giveTempStressVector();
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}




void
CamClayMat :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep)
{
  
    CamClayMatStatus *status = static_cast< CamClayMatStatus * >( this->giveStatus(gp) );
    double pc = status->givePreconsolidationPressure();
    FloatArray plasticStrain;
    FloatArray fullStress;
    // get the initial plastic strain and initial kappa from the status
    plasticStrain = status->givePlasticStrain();
    // elastic predictor
    FloatArray elStrain = totalStrain;
    elStrain.subtract(plasticStrain);
    FloatArray elStrainDev;
    double elStrainVol;
    elStrainVol = computeDeviatoricVolumetricSplit(elStrainDev, elStrain);
    FloatArray trialStressDev;
    applyDeviatoricElasticStiffness(trialStressDev, elStrainDev, G);
    double trialStressVol = 3 * K * elStrainVol;   
    // check the yield condition at the trial state
    double trialQ = sqrt(3. / 2.) * computeStressNorm(trialStressDev); //because q = sqrt(3/2)*||s||
    double trialP = -trialStressVol; //because p=-sigma_m
    double yieldValue = trialQ*trialQ - M*M*trialP*(pc-trialP);
    double lambda = status->givePlasticMultiplier();
	if (yieldValue > 0.) {
		// closest-point projection

		//retrieving elasticity matrix
		FloatMatrix dE;
		give3dMaterialStiffnessMatrix(dE, MatResponseMode::ElasticStiffness, gp, tStep); //for linear elasticity
		FloatMatrix inverseDe;
		inverseDe.beInverseOf(dE);

		//initialization of parametres for the projection algorithm
		int iterationNumber = 0; //step counter
		double deltaLambda = 0.0;
		FloatArray nextPlasticStrain; //eps_p,n+1
		nextPlasticStrain.copySubVector(plasticStrain, 1);//initially there is no increment in plastic strain (equals trial state)
		FloatArray lastStress = status->giveStressVector();//sigma_n
		FloatArray nextStress; //sigma_n+1
		nextStress.beProductOf(dE, elStrain);//initially equals trial stress
		double nextPc = pc;//(p_c)_(n+1), pc will represent (p_c)_n from now on
		bool convergence;

		//computing strain difference needed to achieve
		FloatArray lastStrain = status->giveStrainVector();
		FloatArray strainDifference;
		strainDifference.add(totalStrain);
		strainDifference.times(-1.);
		strainDifference.add(lastStrain);//note - it is already as -delta_eps

		//initial computation of residuals
		//Ri_sigma
		FloatArray rISigma;
		FloatArray sigmaDifference;
		sigmaDifference.add(lastStress);
		sigmaDifference.times(-1.0);
		sigmaDifference.add(nextStress);
		givePlasticGradientAtStress(rISigma, nextStress, nextPc);
		rISigma.times(deltaLambda);
		rISigma.add(strainDifference);
		FloatArray elasticDeformationDifference;
		elasticDeformationDifference.beProductOf(inverseDe, sigmaDifference);
		elasticDeformationDifference.times(-1.0);
		rISigma.add(elasticDeformationDifference);

		//Ri_f
		double rIF = giveYieldValueAtStress(nextStress, nextPc);

		//Ri_p
		double nextP = (-1. / 3.)*(nextStress.at(1) + nextStress.at(2) + nextStress.at(3));
		double rIP = nextPc - theta*deltaLambda*(2 * nextP - nextPc)*nextPc - pc;
      
		do {
			//1) compute all elements of the equation matrix

			//R^sigma_,sigma
			FloatMatrix rSigmaBySigma;
			rSigmaBySigma.add(inverseDe);
			FloatMatrix yieldFunctionDoubleDerivative;
			giveYieldFunctionDoubleDerivative(yieldFunctionDoubleDerivative);
			yieldFunctionDoubleDerivative.times(deltaLambda);
			rSigmaBySigma.add(yieldFunctionDoubleDerivative);
			sigmaDifference.times(0);//cleaning sigmaDifference from previous iteration
			sigmaDifference.add(lastStress);
			sigmaDifference.times(-1.0);
			sigmaDifference.add(nextStress);
			FloatMatrix dEderivative;
			dEderivative.resize(6,6); //this also fills the array with zeros, which is sufficient currently, as stiffness is not stress-dependent
			FloatArray lastPart;
			lastPart.beProductOf(dEderivative, sigmaDifference);
			rSigmaBySigma.add(lastPart);

			//R^sigma_,f = R^f_,sigma
			FloatArray rSigmaByF;
			givePlasticGradientAtStress(rSigmaByF, nextStress, nextPc);

			//R^sigma_,p
			FloatArray rSigmaByP;
			FloatArray doubleDerivative;
			doubleDerivative.resize(6);
			doubleDerivative.at(1) = M*M / 3.;
			doubleDerivative.at(2) = M*M / 3.;
			doubleDerivative.at(3) = M*M / 3.;
			rSigmaByP.beProductOf(dE, doubleDerivative);

			//R^f_,sigma = R^sigma_,f

			//R^f_,f = 0

			//R^f_,p
			double rFByP = 2 * nextPc;

			//R^p_,sigma
			FloatArray rPBySigma;
			rPBySigma.resize(6);
			rPBySigma.at(1) = (2. / 3.)*theta*deltaLambda*nextPc;
			rPBySigma.at(2) = (2. / 3.)*theta*deltaLambda*nextPc;
			rPBySigma.at(3) = (2. / 3.)*theta*deltaLambda*nextPc;

			//R^p_,f
			double rPByF;
			double nextP = (-1. / 3.)*(nextStress.at(1) + nextStress.at(2) + nextStress.at(3));
			rPByF = -theta*(2 * nextP - nextPc)*nextPc;

			//R^p_,p
			double rPByP = 2 * theta*deltaLambda*nextPc + 1.0;

			//2) compute the elements of the right side vector

			//already done at the end of last session/start of cpp

			//3) solve for the three unknowns

			//assembling C
			FloatMatrix inverseCMatrix;
			inverseCMatrix.beDyadicProductOf(rSigmaByP, rPBySigma);
			inverseCMatrix.times((-1. / rPByP));
			inverseCMatrix.add(rSigmaBySigma);
			FloatMatrix cMatrix;
			cMatrix.beInverseOf(inverseCMatrix);

			//solving for deltaDeltaLambda
			FloatArray upperLeftBracket;
			upperLeftBracket.add(rPBySigma);
			upperLeftBracket.times(rFByP / rPByP);
			upperLeftBracket.add(rSigmaByF);
			FloatArray upperRightBracket;
			upperRightBracket.add(rSigmaByP);
			upperRightBracket.times(-rIP / rPByP);
			upperRightBracket.add(rISigma);
			FloatArray upperRightTimesC;
			upperRightTimesC.beProductOf(cMatrix, upperRightBracket);
			double upper = upperLeftBracket.dotProduct(upperRightTimesC);
			upper += -rIF - (rFByP / rPByP)*rIP;

			FloatArray lowerLeftBracket;
			lowerLeftBracket.add(upperLeftBracket); //they are identical
			FloatArray lowerRightBracket;
			lowerRightBracket.add(rSigmaByP);
			lowerRightBracket.times(-rPByF / rPByP);
			lowerRightBracket.add(rSigmaByF);
			FloatArray lowerRightTimesC;
			lowerRightTimesC.beProductOf(cMatrix, lowerRightBracket);
			double lower = lowerLeftBracket.dotProduct(lowerRightTimesC);
			lower += -(rFByP / rPByP)*rPByF;

			double deltaDeltaLambda = upper / lower;

			//solving for deltaSigma
			FloatArray sigmaBracket;
			sigmaBracket.add(lowerRightBracket); //it is identical
			sigmaBracket.times(deltaDeltaLambda);
			sigmaBracket.add(upperRightBracket); //again - the cluster has already been computed
			FloatArray deltaSigma;
			deltaSigma.beProductOf(cMatrix, sigmaBracket);
			deltaSigma.times(-1.);

			//solving for deltaPc
			double deltaPc = -(1. / rPByP)*(rIP + rPBySigma.dotProduct(deltaSigma) + rPByF*deltaDeltaLambda);

			//4) update stress, pc, deltaLambda,
			nextStress.add(deltaSigma);
			nextPc += deltaPc;
			deltaLambda += deltaDeltaLambda;

			//?? what about updating dE??

			//5) update residuals
			//Ri_sigma
			rISigma.times(0);//cleaning
			sigmaDifference.times(0);//cleaning SigmaDifference from the last iteration
			sigmaDifference.add(lastStress);
			sigmaDifference.times(-1.0);
			sigmaDifference.add(nextStress);
			givePlasticGradientAtStress(rISigma, nextStress, nextPc);
			rISigma.times(deltaLambda);
			rISigma.add(strainDifference);
			elasticDeformationDifference.times(0);
			elasticDeformationDifference.beProductOf(inverseDe, sigmaDifference);
			elasticDeformationDifference.times(-1.0);
			rISigma.add(elasticDeformationDifference);

			//Ri_f
			rIF = giveYieldValueAtStress(nextStress, nextPc);

			//Ri_p
			nextP = (-1. / 3.)*(nextStress.at(1) + nextStress.at(2) + nextStress.at(3));
			rIP = nextPc - theta*deltaLambda*(2 * nextP - nextPc)*nextPc - pc;

			//6) check for convergence
			iterationNumber++;
			double plasticStrainError = sqrt(rISigma.dotProduct(rISigma));
			convergence = (fabs(rIF) < yieldTolerance && plasticStrainError < plasticStrainTolerance && fabs(rIP) < pcTolerance);
	
		} while (iterationNumber <= maxIter && !convergence);

		//convergence reached
		//updating all values to values reached by algorithm
		pc = nextPc;
		fullStress = nextStress;

		FloatArray plasticStrainIncrement;
		givePlasticGradientAtStress(plasticStrainIncrement, fullStress, pc);
		plasticStrainIncrement.times(deltaLambda);
		plasticStrain.add(plasticStrainIncrement);

		lambda += deltaLambda;
		if(!convergence) {
			OOFEM_WARNING("Local equlibrium of CPP algorithm not reached in %d iterations, Element number %d, gp %d, continuing", iterationNumber, gp->giveElement()->giveNumber(), gp->giveNumber() );
		}
	
    } else {
		// assemble the stress from the elastically computed volumetric part
		// and scaled deviatoric part    
		computeDeviatoricVolumetricSum(fullStress, trialStressDev, trialStressVol);
    }
    
    
    // store the stress in status
    status->letTempStressVectorBe(fullStress);
    // store the plastic strain, preconsolidation stress, and plastic multiplier
    status->letTempPlasticStrainBe(plasticStrain);
    status->setTempPreconsolidationPressure(pc);
    status->setTempPlasticMultiplier(lambda);
}




// returns the consistent (algorithmic) tangent stiffness matrix
void
CamClayMat :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{
    // start from the elastic stiffness
	FloatMatrix dE;
    this->linearElasticMaterial.give3dMaterialStiffnessMatrix(dE, mode, gp, tStep);
    if ( mode != TangentStiffness ) {
		answer = dE;
        return;
    }
	
    CamClayMatStatus *status = static_cast< CamClayMatStatus * >( this->giveStatus(gp) );
	double tempLambda = status->giveTempPlasticMultiplier();
	double lastLambda = status->givePlasticMultiplier();
	if (tempLambda - lastLambda > 0) {
		//need to assemble the tangent matrix

		//retrieving important parametres

		FloatMatrix inverseDe;
		inverseDe.beInverseOf(dE);

		double deltaLambda = tempLambda - lastLambda;

		FloatArray lastStress;
		lastStress = status->giveStressVector();

		FloatArray nextStress;
		nextStress = status->giveTempStressVector();

		double nextPc = status->giveTempPreconsolidationPressure();

		double pc = status->givePreconsolidationPressure();

		//elements of matrix are the same as in plasticity return
		//therefore the code is too

		//R^sigma_,sigma
		FloatMatrix rSigmaBySigma;
		rSigmaBySigma.add(inverseDe);
		FloatMatrix yieldFunctionDoubleDerivative;
		giveYieldFunctionDoubleDerivative(yieldFunctionDoubleDerivative);
		yieldFunctionDoubleDerivative.times(deltaLambda);
		rSigmaBySigma.add(yieldFunctionDoubleDerivative);
		FloatArray sigmaDifference;
		sigmaDifference.add(lastStress);
		sigmaDifference.times(-1.0);
		sigmaDifference.add(nextStress);
		FloatMatrix dEderivative;
		dEderivative.resize(6,6); //this fills the array with zeros, which is sufficient currently, as stiffness is not stress-dependent
		FloatArray lastPart;
		lastPart.beProductOf(dEderivative, sigmaDifference);
		rSigmaBySigma.add(lastPart);

		//R^sigma_,f = R^f_,sigma
		FloatArray rSigmaByF;
		givePlasticGradientAtStress(rSigmaByF, nextStress, nextPc);

		//R^sigma_,p
		FloatArray rSigmaByP;
		FloatArray doubleDerivative;
		doubleDerivative.resize(6);
		doubleDerivative.at(1) = M*M / 3.;
		doubleDerivative.at(2) = M*M / 3.;
		doubleDerivative.at(3) = M*M / 3.;
		rSigmaByP.beProductOf(dE, doubleDerivative);

		//R^f_,sigma = R^sigma_,f

		//R^f_,f = 0

		//R^f_,p
		double rFByP = 2 * nextPc;

		//R^p_,sigma
		FloatArray rPBySigma;
		rPBySigma.resize(6);
		rPBySigma.at(1) = (2. / 3.)*theta*deltaLambda*nextPc;
		rPBySigma.at(2) = (2. / 3.)*theta*deltaLambda*nextPc;
		rPBySigma.at(3) = (2. / 3.)*theta*deltaLambda*nextPc;

		//R^p_,f
		double rPByF;
		double nextP = (-1. / 3.)*(nextStress.at(1) + nextStress.at(2) + nextStress.at(3));
		rPByF = -theta*(2 * nextP - nextPc)*nextPc;

		//R^p_,p
		double rPByP = 2*theta*deltaLambda*nextPc + 1.0;

		//now solving for D_ep

		//assembling C - copied from plasticityReturn
		FloatMatrix inverseCMatrix;
		inverseCMatrix.beDyadicProductOf(rSigmaByP, rPBySigma);
		inverseCMatrix.times((-1. / rPByP));
		inverseCMatrix.add(rSigmaBySigma);
		FloatMatrix cMatrix;
		cMatrix.beInverseOf(inverseCMatrix);

		//assembling upper part of fraction
		FloatArray upperLeftBracket;
		upperLeftBracket.add(rSigmaByP);
		upperLeftBracket.times(-(rPByF / rPByP));
		upperLeftBracket.add(rSigmaByF);
		FloatArray upperRightBracket;
		upperRightBracket.add(rPBySigma);
		upperRightBracket.times(-(rFByP / rPByP));
		upperRightBracket.add(rSigmaByF);
		FloatArray upperLeftBracketTimesC; //the matrix multiplication is done from left to right in sequence here
		upperLeftBracketTimesC.beProductOf(cMatrix, upperLeftBracket);
		FloatMatrix upperBrackets;
		upperBrackets.beDyadicProductOf(upperLeftBracketTimesC, upperRightBracket);
		FloatMatrix upper;
		upper.beProductOf(upperBrackets, cMatrix);

		//assembling lower part of fraction
		FloatArray lowerLeftBracket;
		lowerLeftBracket.add(upperRightBracket); //identical expressions
		FloatArray lowerRightBracket;
		lowerRightBracket.add(upperLeftBracket); //the same
		FloatArray lowerRightBracketTimesC;
		lowerRightBracketTimesC.beProductOf(cMatrix, upperRightBracket);
		double lower = lowerLeftBracket.dotProduct(lowerRightBracketTimesC);
		lower -= ((rFByP / rPByP) * rPByF); //controversial element

		//finalising fraction
		FloatMatrix fractionResult;
		fractionResult.add(upper);
		fractionResult.times((1 / lower));

		//finalising Dep
		FloatMatrix dEp;
		dEp.add(cMatrix);
		dEp.add(fractionResult);

		answer = dEp;

	}
	else {
		//no change in plastic multiplier -> no plastic step occured -> no need to compute tangent matrix
		answer = dE;
	}
}

void CamClayMat::givePlasticGradientAtStress(FloatArray &answer, const FloatArray &stress, double pc)
{
	if (stress.giveSize() != 6) {
		OOFEM_ERROR("wrong size of stress vector")
	}
	double sigmaM = (1. / 3.)*(stress.at(1) + stress.at(2) + stress.at(3));
	FloatArray sVector = stress;

	sVector.at(1) -= sigmaM;
	sVector.at(2) -= sigmaM;
	sVector.at(3) -= sigmaM;
	sVector.at(4) *= 2;
	sVector.at(5) *= 2;
	sVector.at(6) *= 2;

	answer.resize(6);
	answer.at(1) = 3 * sVector.at(1) + (M*M / 3.) * (2 * sigmaM + pc);
	answer.at(2) = 3 * sVector.at(2) + (M*M / 3.) * (2 * sigmaM + pc);
	answer.at(3) = 3 * sVector.at(3) + (M*M / 3.) * (2 * sigmaM + pc);
	answer.at(4) = 3 * sVector.at(4);
	answer.at(5) = 3 * sVector.at(5);
	answer.at(6) = 3 * sVector.at(6);
}

void CamClayMat :: giveYieldFunctionDoubleDerivative(FloatMatrix &answer)
{

	//creating derivative of svector by sigma (=double derivative of J2 by sigma)
	FloatMatrix sMatrix;
	sMatrix.resize(6, 6);
	sMatrix.zero();

	sMatrix.at(1, 1) = (2. / 3.);
	sMatrix.at(2, 2) = (2. / 3.);
	sMatrix.at(3, 3) = (2. / 3.);
	sMatrix.at(1, 2) = (-1. / 3.);
	sMatrix.at(1, 3) = (-1. / 3.);
	sMatrix.at(2, 1) = (-1. / 3.);
	sMatrix.at(2, 3) = (-1. / 3.);
	sMatrix.at(3, 1) = (-1. / 3.);
	sMatrix.at(3, 2) = (-1. / 3.);

	sMatrix.at(4, 4) = 2.;
	sMatrix.at(5, 5) = 2.;
	sMatrix.at(6, 6) = 2.;

	//creating derivative of sigma_m*delta = {sigma_m, sigma_m, sigma_m, 0, 0, 0} by sigma
	FloatMatrix deltaMatrix;
	deltaMatrix.resize(6, 6);
	deltaMatrix.zero();

	deltaMatrix.at(1, 1) = (1. / 3.);
	deltaMatrix.at(1, 2) = (1. / 3.);
	deltaMatrix.at(1, 3) = (1. / 3.);
	deltaMatrix.at(2, 1) = (1. / 3.);
	deltaMatrix.at(2, 2) = (1. / 3.);
	deltaMatrix.at(2, 3) = (1. / 3.);
	deltaMatrix.at(3, 1) = (1. / 3.);
	deltaMatrix.at(3, 2) = (1. / 3.);
	deltaMatrix.at(3, 3) = (1. / 3.);

	//derivative = 3*S + (2/3)*M^2*DELTA
	answer = sMatrix;
	answer.times(3.);
	deltaMatrix.times((2. / 3.)*M*M);
	answer.add(deltaMatrix);
}


int
CamClayMat :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    CamClayMatStatus *status = static_cast< CamClayMatStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        answer = status->givePlasticStrain();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

//=============================================================================

  CamClayMatStatus :: CamClayMatStatus(int n, Domain *d, GaussPoint *g, double pc0) :
    StructuralMaterialStatus(n, d, g), plasticStrain(6), tempPlasticStrain()
{
    stressVector.resize(6);
    strainVector.resize(6);
    pc = tempPc = pc0;
    plasticMultiplier = tempPlasticMultiplier = 0.;
}

CamClayMatStatus :: ~CamClayMatStatus()
{ }

void
CamClayMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "              plastic  ");
    for ( auto &val : this->plasticStrain ) {
        fprintf(file, "%.4e ", val);
    }
    fprintf(file, "}\n");

}


// initializes temporary variables based on their values at the previous equlibrium state
void CamClayMatStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    
    tempPlasticStrain = plasticStrain;
    tempPc = pc;
    tempPlasticMultiplier = plasticMultiplier;
    
   
}


// updates internal variables when equilibrium is reached
void
CamClayMatStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    pc = tempPc;
    plasticMultiplier = tempPlasticMultiplier;
   
}





// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
CamClayMatStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data

    // write plastic strain (vector)
    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write pre-consolidation pressure
    if ( !stream.write(pc) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    return CIO_OK;
}



contextIOResultType
CamClayMatStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read plastic strain (vector)
    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read pre-consolidation pressure
    if ( !stream.read(pc) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem