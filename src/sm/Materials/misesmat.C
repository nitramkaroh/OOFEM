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

#include "misesmat.h"
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
REGISTER_Material(MisesMat);

// constructor
MisesMat :: MisesMat(int n, Domain *d) : StructuralMaterial(n, d)
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    H = 0.;
    sig0 = 0.;
    G = 0.;
    K = 0.;
}

// destructor
MisesMat :: ~MisesMat()
{
    delete linearElasticMaterial;
}

// reads the model parameters from the input file
IRResultType
MisesMat :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    result = linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    if ( result != IRRT_OK ) return result;

    G = static_cast< IsotropicLinearElasticMaterial * >(linearElasticMaterial)->giveShearModulus();
    K = static_cast< IsotropicLinearElasticMaterial * >(linearElasticMaterial)->giveBulkModulus();

    H = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, H, _IFT_MisesMat_h); // hardening modulus


    // isotropic hardening parameters
    sigInf = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, sigInf, _IFT_MisesMat_sigInf); // hardening modulus
    sig0 = 0.;
    IR_GIVE_FIELD(ir, sig0, _IFT_MisesMat_sig0); // uniaxial yield stress
    expD = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, expD, _IFT_MisesMat_expD); // hardening modulus


    
    /*********************************************************************************************************/




    
    omega_crit = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, omega_crit, _IFT_MisesMat_omega_crit); // critical damage

    a = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, a, _IFT_MisesMat_a); // exponent in damage law
    /********************************************************************************************************/

    return IRRT_OK;
}

// creates a new material status  corresponding to this class
MaterialStatus *
MisesMat :: CreateStatus(GaussPoint *gp) const
{
    return new MisesMatStatus(1, this->giveDomain(), gp);
}

void
MisesMat :: giveRealStressVector_1d(FloatArray &answer,
                                 GaussPoint *gp,
                                 const FloatArray &totalStrain,
                                 TimeStep *tStep)
{
    /// @note: One should obtain the same answer using the iterations in the default implementation (this is verified for this model).
#if 1
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );

    this->performPlasticityReturn(gp, totalStrain);
    double omega = computeDamage(gp, tStep);
    answer = status->giveTempEffectiveStress();
    answer.times(1 - omega);

    // Compute the other components of the strain:
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    double E = lmat->give('E', gp), nu = lmat->give('n', gp);

    FloatArray strain = status->getTempPlasticStrain();
    strain[0] = totalStrain[0];
    strain[1] -= nu / E * status->giveTempEffectiveStress()[0];
    strain[2] -= nu / E * status->giveTempEffectiveStress()[0];

    status->letTempStrainVectorBe(strain);
    status->setTempDamage(omega);
    status->letTempStressVectorBe(answer);
#else
    StructuralMaterial :: giveRealStressVector_1d(answer, gp, totalStrain, tStep);
#endif
}

void
MisesMat :: giveRealStressVector_3d(FloatArray &answer,
                                 GaussPoint *gp,
                                 const FloatArray &totalStrain,
                                 TimeStep *tStep)
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
   // initialization
    this->initTempStatus(gp);
    this->performPlasticityReturn(gp, totalStrain);
    double omega = computeDamage(gp, tStep);
    answer = status->giveTempEffectiveStress();
    answer.times(1 - omega);
    status->setTempDamage(omega);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}





void
MisesMat :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain)
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    double kappa;
    FloatArray plStrain;
    FloatArray fullStress;
    // get the initial plastic strain and initial kappa from the status
    plStrain = status->givePlasticStrain();
    kappa = status->giveCumulativePlasticStrain();

    // === radial return algorithm ===
    if ( totalStrain.giveSize() == 1 ) {
        LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
        double E = lmat->give('E', gp);
        /*trial stress*/
        fullStress.resize(6);
        fullStress.at(1) = E * ( totalStrain.at(1) - plStrain.at(1) );
        double trialS = fabs(fullStress.at(1));
        /*yield function*/
        double yieldValue = trialS - (sig0 + H * kappa);
        // === radial return algorithm ===
        if ( yieldValue > 0 ) {
            double dKappa = yieldValue / ( H + E );
            kappa += dKappa;
            plStrain.at(1) += dKappa * signum( fullStress.at(1) );
            plStrain.at(2) -= 0.5 * dKappa * signum( fullStress.at(1) );
            plStrain.at(3) -= 0.5 * dKappa * signum( fullStress.at(1) );
            fullStress.at(1) -= dKappa * E * signum( fullStress.at(1) );
        }
    } else {
        // elastic predictor
        FloatArray elStrain = totalStrain;
        elStrain.subtract(plStrain);
        FloatArray elStrainDev;
        double elStrainVol;
        elStrainVol = computeDeviatoricVolumetricSplit(elStrainDev, elStrain);
        FloatArray trialStressDev;
        applyDeviatoricElasticStiffness(trialStressDev, elStrainDev, G);
        /**************************************************************/
        double trialStressVol = 3 * K * elStrainVol;
        /**************************************************************/

        // store the deviatoric and trial stress (reused by algorithmic stiffness)
        status->letTrialStressDevBe(trialStressDev);
        status->setTrialStressVol(trialStressVol);
        // check the yield condition at the trial state
        double trialS = computeStressNorm(trialStressDev);
        double yieldValue = sqrt(3./2.) * trialS - (sig0 + H * kappa);
	double sigmaY = this->computeYieldStress(kappa);
	double qtrial = sqrt(3./2.) * trialS;
	yieldValue = sqrt(3./2.) * trialS - sigmaY;
	if ( yieldValue > 0. ) {
	  double dKappa = 0;
	  while(true) {
	    double HiP = this->computeYieldStressPrime(kappa + dKappa);
	    double Hi = this->computeYieldStress(kappa + dKappa);
	    double g =  sqrt(3./2.) * trialS - 3. * G * dKappa - Hi;
	    double Dg = - 3. * G - HiP;
	    // increment of cumulative plastic strain
	    dKappa -= g/Dg;
	    if(fabs(g) < 1.e-10) {
	      break;
	    }
	  }
	  kappa += dKappa;
	  FloatArray dPlStrain;
	  // the following line is equivalent to multiplication by scaling matrix P
	  applyDeviatoricElasticCompliance(dPlStrain, trialStressDev, 0.5);
	  // increment of plastic strain
	  plStrain.add(sqrt(3. / 2.) * dKappa / trialS, dPlStrain);
	  // scaling of deviatoric trial stress
	  trialStressDev.times(1. - sqrt(6.) * G * dKappa / trialS);
	}
	
    
	/*        if ( yieldValue > 0. ) {
            // increment of cumulative plastic strain
            double dKappa = yieldValue / ( H + 3. * G );
            kappa += dKappa;
            FloatArray dPlStrain;
            // the following line is equivalent to multiplication by scaling matrix P
            applyDeviatoricElasticCompliance(dPlStrain, trialStressDev, 0.5);
            // increment of plastic strain
            plStrain.add(sqrt(3. / 2.) * dKappa / trialS, dPlStrain);
            // scaling of deviatoric trial stress
            trialStressDev.times(1. - sqrt(6.) * G * dKappa / trialS);

	    trialS = computeStressNorm(trialStressDev);
	    yieldValue = sqrt(3./2.) * trialS - (sig0 + H * kappa);

	    }*/

        // assemble the stress from the elastically computed volumetric part
        // and scaled deviatoric part

        computeDeviatoricVolumetricSum(fullStress, trialStressDev, trialStressVol);
    }

    // store the effective stress in status
    status->letTempEffectiveStressBe(fullStress);
    // store the plastic strain and cumulative plastic strain
    status->letTempPlasticStrainBe(plStrain);
    status->setTempCumulativePlasticStrain(kappa);
}

double
MisesMat :: computeDamageParam(double tempKappa)
{
    if ( tempKappa > 0. ) {
        return omega_crit * ( 1.0 - exp(-a * tempKappa) );
    } else {
        return 0.;
    }
}

double
MisesMat :: computeDamageParamPrime(double tempKappa)
{
    if ( tempKappa >= 0. ) {
        return omega_crit * a * exp(-a * tempKappa);
    } else {
        return 0.;
    }
}



double
MisesMat :: computeDamage(GaussPoint *gp,  TimeStep *tStep)
{
    double tempKappa, dam;
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    dam = status->giveDamage();
    computeCumPlastStrain(tempKappa, gp, tStep);
    double tempDam = computeDamageParam(tempKappa);
    if ( dam > tempDam ) {
        tempDam = dam;
    }

    return tempDam;
}


void MisesMat :: computeCumPlastStrain(double &tempKappa, GaussPoint *gp, TimeStep *tStep)
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    tempKappa = status->giveTempCumulativePlasticStrain();
}



// returns the consistent (algorithmic) tangent stiffness matrix
void
MisesMat :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{

    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );

    // start from the elastic stiffness
    this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
    if ( mode != TangentStiffness )  {
      double omega = status->giveTempDamage();
      answer.times(1. - omega);
      return;
    }

    double kappa = status->giveCumulativePlasticStrain();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = tempKappa - kappa;

    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return;
    }

    // === plastic loading ===    
    // yield stress at the beginning of the step
    double sigmaY = this->computeYieldStress(kappa);
    double HiP = this->computeYieldStressPrime(kappa + dKappa);
      
      //sig0 + H * kappa;

    // trial deviatoric stress and its norm
    const FloatArray &trialStressDev = status->giveTrialStressDev();
    //double trialStressVol = status->giveTrialStressVol();
    double trialS = computeStressNorm(trialStressDev);

    double qtrial = sqrt(3./2.)*trialS;

    double afact = 2.*G * (1-3.*G*dKappa/qtrial);
    double bfact = 6.*G*G*(dKappa/qtrial-1/(3.*G+HiP))/trialS/trialS;

    FloatArray delta;
    FloatMatrix sc, stiffnessCorrection1, stiffnessCorrection2, stiffnessCorrection3;
    stiffnessCorrection1.bePinvID();
    stiffnessCorrection1.times(afact);
    stiffnessCorrection2.beDyadicProductOf(trialStressDev, trialStressDev);
    stiffnessCorrection2.times(bfact);
    delta = {1,1,1,0,0,0};
    stiffnessCorrection3.beDyadicProductOf(delta, delta);
    stiffnessCorrection3.times(K);
    sc = stiffnessCorrection1;
    sc.add(stiffnessCorrection2);
    sc.add(stiffnessCorrection3);
    
    // one correction term
    FloatMatrix stiffnessCorrection;
    stiffnessCorrection.beDyadicProductOf(trialStressDev, trialStressDev);
    double factor = -2. * sqrt(6.) * G * G / trialS;
    //    double factor1 = factor * sigmaY / ( ( H + 3. * G ) * trialS * trialS );
    double factor1 = factor * sigmaY / ( ( 3.*G + HiP ) * trialS * trialS );
    answer.add(factor1, stiffnessCorrection);

    // another correction term
    stiffnessCorrection.bePinvID();
    double factor2 = factor * dKappa;
    answer.add(factor2, stiffnessCorrection);

    //influence of damage
    //    double omega = computeDamageParam(tempKappa);
    double omega = status->giveTempDamage();
    answer.times(1. - omega);
    const FloatArray &effStress = status->giveTempEffectiveStress();
    double omegaPrime = computeDamageParamPrime(tempKappa);
    //    double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + H ) / trialS;
    double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + HiP ) / trialS;
    stiffnessCorrection.beDyadicProductOf(effStress, trialStressDev);
    stiffnessCorrection.times(scalar);
    answer.add(stiffnessCorrection);

    //answer = sc;

    
}

void
MisesMat :: give1dStressStiffMtrx(FloatMatrix &answer,
                                  MatResponseMode mode,
                                  GaussPoint *gp,
                                  TimeStep *tStep)
{
    this->giveLinearElasticMaterial()->give1dStressStiffMtrx(answer, mode, gp, tStep);
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    double kappa = status->giveCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double omega = status->giveTempDamage();
    double E = answer.at(1, 1);
    if ( mode != TangentStiffness ) {
        return;
    }


    if ( tempKappa <= kappa ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        answer.times(1 - omega);
        return;
    }


    // === plastic loading ===
    const FloatArray &stressVector = status->giveTempEffectiveStress();
    double stress = stressVector.at(1);
    answer.resize(1, 1);
    answer.at(1, 1) = ( 1 - omega ) * E * H / ( E + H ) - computeDamageParamPrime(tempKappa) * E / ( E + H ) * stress * signum(stress);
}





double
MisesMat :: computeYieldStress(double kappa)
{

  /*  FloatArray kA = {0.000,0.001,0.002,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.12,0.14,0.16,0.20,0.25,0.3,0.5,1.0,11.0};
  FloatArray hA = {0.45,0.454577,0.459081,0.472155,0.492564,0.528701,0.559411,0.585539,0.607799,0.626794,0.643032,0.656942,0.668887,0.679172,0.695760,0.708326,0.718025,0.731879,0.743463,0.752122,0.779564,0.84424,2.13664};
  int i = 0;
  while(true) {
    i++;
    if(kappa >= kA.at(i) && kappa <= kA.at(i+1)) {
      break;
    }
  }
  double sY = hA.at(i) + (hA.at(i+1)-hA.at(i))/(kA.at(i+1)-kA.at(i))*(kappa-kA.at(i));
  
  return sY;
  */
  
  
  return this->sig0 + this-> H * kappa;// + ( this->sigInf - this->sig0 ) * (1. - exp(-expD*kappa));
}

double 
MisesMat :: computeYieldStressPrime(double kappa)
{
  /*
  FloatArray kA = {0.000,0.001,0.002,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.12,0.14,0.16,0.20,0.25,0.3,0.5,1.0,11.0};
  FloatArray hA = {0.45,0.454577,0.459081,0.472155,0.492564,0.528701,0.559411,0.585539,0.607799,0.626794,0.643032,0.656942,0.668887,0.679172,0.695760,0.708326,0.718025,0.731879,0.743463,0.752122,0.779564,0.84424,2.13664};
  
  int i = 0;
  while(true) {
    i++;
    if(kappa >= kA.at(i) && kappa <= kA.at(i+1)) {
      break;
    }
   
  }
  double sY = (hA.at(i+1)-hA.at(i))/(kA.at(i+1)-kA.at(i));
  return sY;
  */
  
  return this->H;// + ( this->sigInf - this->sig0 ) * expD * exp(-expD*kappa); 
}





int
MisesMat :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    MisesMatStatus *status = static_cast< MisesMatStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
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
    }
}

//=============================================================================

MisesMatStatus :: MisesMatStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), plasticStrain(6), tempPlasticStrain(), trialStressD(6)
{
    stressVector.resize(6);
    strainVector.resize(6);
    
    damage = tempDamage = 0.;
    kappa = tempKappa = 0.;
    effStress.resize(6);
    tempEffStress.resize(6);

}

MisesMatStatus :: ~MisesMatStatus()
{ }

void
MisesMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "              plastic  ");
    for ( auto &val : this->plasticStrain ) {
        fprintf(file, "%.4e ", val);
    }
    fprintf(file, "\n");

    fprintf(file, "status { ");
    /*
     * // this would not be correct, since the status is already updated and kappa is now the "final" value
     * if ( tempKappa > kappa ) {
     * fprintf(file, " Yielding, ");
     * } else {
     * fprintf(file, " Unloading, ");
     * }
     */

    /*********************************************************************************/

    fprintf(file, "damage %.4e", tempDamage);

    /********************************************************************************/
    /*
     * // print the plastic strain
     *  n = plasticStrain.giveSize();
     *  fprintf(file, " plastic_strains ");
     *  for ( i = 1; i <= n; i++ ) {
     *      fprintf( file, " %.4e", plasticStrain.at(i) );
     *  }
     */
    // print the cumulative plastic strain
    fprintf(file, ", kappa ");
    fprintf(file, " %.4e", kappa);

    fprintf(file, "}\n");
    /*
     * //print Left Cauchy - Green deformation tensor
     * fprintf(file," left Cauchy Green");
     * fprintf( file, " %.4e",tempLeftCauchyGreen.at(1,1) );
     * fprintf( file, " %.4e",tempLeftCauchyGreen.at(2,2) );
     * fprintf( file, " %.4e",tempLeftCauchyGreen.at(3,3) );
     * fprintf( file, " %.4e",tempLeftCauchyGreen.at(2,3) );
     * fprintf( file, " %.4e",tempLeftCauchyGreen.at(1,3) );
     * fprintf( file, " %.4e",tempLeftCauchyGreen.at(1,2) );
     *
     * //print deformation gradient
     * fprintf(file," Deformation Gradient");
     * fprintf( file, " %.4e",tempDefGrad.at(1,1) );
     * fprintf( file, " %.4e",tempDefGrad.at(1,2) );
     * fprintf( file, " %.4e",tempDefGrad.at(1,3) );
     * fprintf( file, " %.4e",tempDefGrad.at(2,1) );
     * fprintf( file, " %.4e",tempDefGrad.at(2,2) );
     * fprintf( file, " %.4e",tempDefGrad.at(2,3) );
     * fprintf( file, " %.4e",tempDefGrad.at(3,1) );
     * fprintf( file, " %.4e",tempDefGrad.at(3,2) );
     * fprintf( file, " %.4e",tempDefGrad.at(3,3) );
     */

    fprintf(file, "}\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void MisesMatStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    tempDamage = damage;
    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
    trialStressD.clear(); // to indicate that it is not defined yet
}


// updates internal variables when equilibrium is reached
void
MisesMatStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    kappa = tempKappa;
    damage = tempDamage;
    trialStressD.clear(); // to indicate that it is not defined any more
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
MisesMatStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
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

    // write cumulative plastic strain (scalar)
    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write damage (scalar)
    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}



contextIOResultType
MisesMatStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
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

    // read cumulative plastic strain (scalar)
    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read damage (scalar)
    if ( !stream.read(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
