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

#include "variationalbaseddamage.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "Materials/isolinearelasticmaterial.h"
#include "error.h"
#include "nonlocalmaterialext.h"
#include "datastream.h"
#include "contextioerr.h"
#include "stressvector.h"
#include "strainvector.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(VarBasedDamageMaterial);

//q Syntax?

VarBasedDamageMaterial :: VarBasedDamageMaterial(int n, Domain *d) : IsotropicDamageMaterial1(n, d), GradientDamageMaterialExtensionInterface(d)
    //
    // constructor
    //
{
}

VarBasedDamageMaterial :: ~VarBasedDamageMaterial()
//
// destructor
//
{ }


//q IRResultType --- syntax?  
 
IRResultType
VarBasedDamageMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    result = linearElasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    result = GradientDamageMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    IR_GIVE_FIELD(ir, gf, _IFT_IsotropicDamageMaterial1_gf);
    this->phaseFieldModelType = phaseFieldModel_JZ;

    int phaseFieldModelTypeRecord = 0; // default
    //non zero value corresponds to the Miehe phase-field model or Wu model
    IR_GIVE_OPTIONAL_FIELD(ir, phaseFieldModelTypeRecord, _IFT_VarBasedDamageMaterial_phaseFieldModelType);
    if ( phaseFieldModelTypeRecord == 0 ) {
        this->phaseFieldModelType = phaseFieldModel_JZ;
    } else if ( phaseFieldModelTypeRecord == 1 ) {
        this->phaseFieldModelType = phaseFieldModel_Miehe;
    } else if ( phaseFieldModelTypeRecord == 2 ) {
        this->phaseFieldModelType = phaseFieldModel_Wu;
    }  else {
        OOFEM_ERROR("Unknown phase-field model typed");
    }
   
    
    if (this->phaseFieldModelType == phaseFieldModel_JZ){
      this->p = 0.5;
      IR_GIVE_OPTIONAL_FIELD(ir, this->p, _IFT_VarBasedDamageMaterial_p);
      this->damageLaw = 0; // linear softening - default
      IR_GIVE_OPTIONAL_FIELD(ir, this->damageLaw, _IFT_VarBasedDamageMaterial_damageLaw);
       beta = 1.; // elastic-brittle model - default
      IR_GIVE_OPTIONAL_FIELD(ir, this->beta, _IFT_VarBasedDamageMaterial_beta);
    }
    else if (this->phaseFieldModelType == phaseFieldModel_Miehe) {
      this->damageLaw = 1; // damage law used by Miehe 
    }
    else if (this->phaseFieldModelType == phaseFieldModel_Wu) {
      this->p = 2.; // linear softening case for Wu
      IR_GIVE_OPTIONAL_FIELD(ir, this->p, _IFT_VarBasedDamageMaterial_p);
      this->a1 = 1.;
      IR_GIVE_OPTIONAL_FIELD(ir, this->a1, _IFT_VarBasedDamageMaterial_a1);
      this->a2 = -2.;
      IR_GIVE_OPTIONAL_FIELD(ir, this->a2, _IFT_VarBasedDamageMaterial_a2); 
      this->a3 = 0.;
      IR_GIVE_OPTIONAL_FIELD(ir, this->a3, _IFT_VarBasedDamageMaterial_a3);
      this->damageLaw = 5; // damage law used by Wu       
    }

    int equivStrainTypeRecord = 0; // default
    IR_GIVE_OPTIONAL_FIELD(ir, equivStrainTypeRecord, _IFT_VarBasedDamageMaterial_equivstraintype);
    if ( equivStrainTypeRecord == 0 ) {
        this->equivStrainType = EST_Mazars;
    } else if ( equivStrainTypeRecord == 1 ) {
        this->equivStrainType = EST_Rankine_Smooth;
    } else if ( equivStrainTypeRecord == 2 ) {
        this->equivStrainType = EST_ElasticEnergy;
    } else if ( equivStrainTypeRecord == 3 ) {
        this->equivStrainType = EST_Mises;
        IR_GIVE_FIELD(ir, k, _IFT_IsotropicDamageMaterial1_k);
    } else if ( equivStrainTypeRecord == 4 ) {
        this->equivStrainType = EST_Rankine_Standard;
    } else if ( equivStrainTypeRecord == 5 ) {
        this->equivStrainType = EST_ElasticEnergyPositiveStress;
    } else if ( equivStrainTypeRecord == 6 ) {
        this->equivStrainType = EST_ElasticEnergyPositiveStrain;
    } else if ( equivStrainTypeRecord == 7 ) {
        this->equivStrainType = EST_Griffith;
        IR_GIVE_OPTIONAL_FIELD(ir, griff_n, _IFT_IsotropicDamageMaterial1_n);
    } else {
        OOFEM_WARNING("Unknown equivStrainType %d", equivStrainType);
        return IRRT_BAD_FORMAT;
    }

    return this->mapper.initializeFrom(ir);
}

/////////////////////////////////////////////////////////////////////////////



int
VarBasedDamageMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _1dMat || mode == _PlaneStress || mode == _PlaneStrain || mode == _3dMat;
}

void
VarBasedDamageMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
                                   MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}


void
VarBasedDamageMaterial :: giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
  double tempDamage;
  if ( mode == ElasticStiffness ) {
    tempDamage = 0.0;
  } else {
    tempDamage = status->giveTempDamage();
    if ( tempDamage > 0.0 ) {
      tempDamage = min(tempDamage, maxOmega);
    }
  }
  
  this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, mode, gp, tStep);
  answer.times(1.0 - tempDamage);  

}

  

void
VarBasedDamageMaterial :: giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );

  double dDamage, damageDrivingVariable;
  FloatArray stress = status->giveTempEffectiveStressVector();
  damageDrivingVariable = status->giveTempNonlocalDamageDrivingVariable();
  this->computeDamagePrime(dDamage, damageDrivingVariable, gp);
  answer = stress;

  double deltaDamage = status->giveTempDamage() - status->giveDamage();
  if(deltaDamage <= 0) {
    answer.times(0);
  } else {
    answer.times(-dDamage);
  }
  /// zero block for now
  answer.times(0);
  
}


  // this function is not necessary for variational, i.e., a symmetric formulation
  // it is implemented anyway for possibility of further extensions
void
VarBasedDamageMaterial :: giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );

  FloatArray stress = status->giveTempEffectiveStressVector();
  double dDamage, damageDrivingVariable;
  damageDrivingVariable = status->giveTempNonlocalDamageDrivingVariable();
  this->computeDamagePrime(dDamage, damageDrivingVariable, gp);
  answer = stress;

  double deltaDamage = status->giveTempDamage() - status->giveDamage();
  if(deltaDamage <= 0) {
    answer.times(0);
  } else {
    answer.times(-dDamage);
  }
  /// zero block for now
  answer.times(0);

}

void
VarBasedDamageMaterial :: giveGradientDamageStiffnessMatrix_dd_NN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{  
  double localDamageDrivingVariable, damageDrivingVariable, dDamage,ddDamage, dDiss, ddDiss, damage;

    VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
    damageDrivingVariable = status->giveTempNonlocalDamageDrivingVariable();
    damage = status->giveTempDamage();

  // stored elastic energy
  this-> computeLocalDamageDrivingVariable(localDamageDrivingVariable, gp, tStep);
  
  this->computeDamagePrime(dDamage, damageDrivingVariable, gp);
  this->computeDamagePrime2(ddDamage, damageDrivingVariable, gp);
  this->computeDissipationFunctionPrime(dDiss, damage, damageDrivingVariable, gp);
  this->computeDissipationFunctionPrime2(ddDiss, damage, damageDrivingVariable, gp);

  answer.resize(1,1);
  answer.at(1,1) = (ddDiss*dDamage*dDamage) + (dDiss-localDamageDrivingVariable)*ddDamage;


}

void
VarBasedDamageMaterial :: giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _1dMat) {
    answer = {{1}};
  } else if(matMode == _PlaneStress || matMode == _PlaneStrain) {
    answer = {{1, 0},{0, 1}};
  } else if(matMode == _3dMat) {
    answer = {{1, 0, 0},{0, 1, 0}, {0, 0, 1}};
  }
  answer.times(gf*internalLength*internalLength);  
}

  



void
VarBasedDamageMaterial :: giveNonlocalInternalForces_N_factor(double &answer, double nonlocalDamageDrivingVariable, GaussPoint *gp, TimeStep *tStep)

  
{
  double localDamageDrivingVariable, dDamage, dDiss, damage;

  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
    //  damageDrivingVariable = status->giveTempNonlocalDamageDrivingVariable();
  damage = status->giveTempDamage();
  
  this-> computeLocalDamageDrivingVariable(localDamageDrivingVariable, gp, tStep);
  this->computeDamagePrime(dDamage, nonlocalDamageDrivingVariable, gp);
  this->computeDissipationFunctionPrime(dDiss, damage, nonlocalDamageDrivingVariable,  gp);  
  answer = ( dDiss - localDamageDrivingVariable ) * dDamage;
 
}

void
VarBasedDamageMaterial :: giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nonlocalDamageDrivingVariableGrad, GaussPoint *gp, TimeStep *tStep)
{
  answer = nonlocalDamageDrivingVariableGrad;
  answer.times(gf * internalLength * internalLength);
#ifdef keep_track_of_dissipated_energy
  this->computeRegulirizingWork(gp, nonlocalDamageDrivingVariableGrad);
#endif
}

double
VarBasedDamageMaterial :: computeQ_Wu(double damageDrivingVariable)
{
  return this->a1*damageDrivingVariable + this->a2*pow(damageDrivingVariable,2) + this->a3*pow(damageDrivingVariable,3);
}

double
VarBasedDamageMaterial :: computeQ_Wu_prime1(double damageDrivingVariable)
{
  return this->a1 + this->a2*2*damageDrivingVariable + this->a3*3*pow(damageDrivingVariable,2);
}

double
VarBasedDamageMaterial :: computeQ_Wu_prime2(double damageDrivingVariable)
{
  return this->a2*2 + this->a3*6*damageDrivingVariable;
}


void
VarBasedDamageMaterial :: computeDamage(double &answer, double damageDrivingVariable, GaussPoint *gp)
{
  /// for now, class of models using linear softening with exponent p is used

  if (this->phaseFieldModelType == phaseFieldModel_JZ) { 
    if(this->p == 1) {
      answer = 1-exp( -damageDrivingVariable );
    } else {   
      answer = 1. - pow( ( 1. + ( this->p - 1. ) * damageDrivingVariable ), 1. / ( 1. - this->p ) );
    }
  }
  else if(this->phaseFieldModelType == phaseFieldModel_Miehe) {
    //corresponds to the Miehe phase-field model
    answer = 2*damageDrivingVariable-damageDrivingVariable*damageDrivingVariable;
  }
  else if(this->phaseFieldModelType == phaseFieldModel_Wu) {
    //corresponds to the Wu phase-field model
    double Q = this->computeQ_Wu(damageDrivingVariable);  
    answer = 1. - pow(1 - damageDrivingVariable, this->p)/(Q + pow(1 - damageDrivingVariable,this->p)); 
  }
  else {
    OOFEM_ERROR("Unknown phase-field model type");
  }
  
  if(answer > 1.) {
    answer = 1.;
  }
  if(answer < 0.) {
    answer = 0.;
  }
}

void
VarBasedDamageMaterial :: computeDamagePrime(double &answer, double damageDrivingVariable, GaussPoint *gp)
{
  /// for now, only class of models using exponent p is used
  if (this->phaseFieldModelType == phaseFieldModel_JZ) {
    if(this->p == 1) {
      answer = exp( -damageDrivingVariable );
    } else {   
      answer = pow( ( 1. + ( this->p - 1. ) * damageDrivingVariable ), this->p / ( 1. - this->p ) );
    }
  }
  else if(this->phaseFieldModelType == phaseFieldModel_Miehe) {
    //corresponds to the Miehe phase-field model
    answer = 2 * ( 1. - damageDrivingVariable );
  }
  else if(this->phaseFieldModelType == phaseFieldModel_Wu) {
    //corresponds to the Wu phase-field model
    double Q = this->computeQ_Wu(damageDrivingVariable);
    double prime1_Q = this->computeQ_Wu_prime1(damageDrivingVariable);
    answer = ((1.-damageDrivingVariable)*prime1_Q+p*Q)*pow(1-damageDrivingVariable,p-1.)/pow(Q+pow(1.-damageDrivingVariable,p),2);
  }
  else {
    OOFEM_ERROR("Unknown phase-field model type");
  }
}
  
void
VarBasedDamageMaterial :: computeDamagePrime2(double &answer, double damageDrivingVariable, GaussPoint *gp)
{
  /// for now, only class of models using exponent p is used
  if (this->phaseFieldModelType == phaseFieldModel_JZ) {
    if(this->p == 1) {
      answer = - exp( - damageDrivingVariable );
    } else {   
      answer = - this->p * pow( ( 1. + ( this->p - 1. ) * damageDrivingVariable ), (2. * this->p - 1) / ( 1. - this->p ) );
    }
  }
  else if(this->phaseFieldModelType == phaseFieldModel_Miehe) {
    //corresponds to the Miehe phase-field model
    answer = -2;
  }

  else if(this-> phaseFieldModelType == phaseFieldModel_Wu) {
    //corresponds to the Wu phase-field model
    double Q = this->computeQ_Wu(damageDrivingVariable);
    double prime1_Q = this->computeQ_Wu_prime1(damageDrivingVariable);
    double prime2_Q = this->computeQ_Wu_prime2(damageDrivingVariable);
    answer = pow(1. - damageDrivingVariable,p-2.)*(-(this->p-1.)*this->p*pow(Q,2)+Q*(pow(1.-damageDrivingVariable,this->p)*this->p*(1.+this->p) + 2*(damageDrivingVariable - 1.)*this->p*prime1_Q + pow(damageDrivingVariable -1.,2)*prime2_Q)+(damageDrivingVariable -1.)*(-2*(damageDrivingVariable -1.)*pow(prime1_Q,2) + pow(1. - damageDrivingVariable,this->p)*(-2*this->p*prime1_Q + (-1. + damageDrivingVariable)*prime2_Q)))/pow(pow(1. - damageDrivingVariable, this->p) + Q,3);
  }
  else {
    OOFEM_ERROR("Unknown phase-field model type");
  }
}


double
VarBasedDamageMaterial :: solveExpLaw(double dam, double c)
{
  double answer = 0.;
  double f, df;
  for (int i=1; i<1000; i++){
    f = (1.-dam)*(answer+1.)-exp(-c*answer);
    df = 1.-dam+c*exp(-c*answer);
    answer -= f/df;
    if (fabs(f)<1.e-8)
      return answer;
  }
  printf("No convergence in VarBasedDamageMaterial :: solveExpLaw(%g,%g)\n",dam,c);
  exit(0);
}

double
VarBasedDamageMaterial :: compute_dissipation_Wu_prime1_in_gamma(double damageDrivingVariable)
{
  return gf*(2 - 2*damageDrivingVariable);
}

double
VarBasedDamageMaterial :: compute_dissipation_Wu_prime2_in_gamma(double damageDrivingVariable)
{
  return -2*gf;
}

void
VarBasedDamageMaterial :: computeDissipationFunctionPrime(double &answer, double damage, double damageDrivingVariable,  GaussPoint *gp)
{
  /// softening - currently using 0=linear (default), 1=Miehe, 2=Fremond, 3=pseudo-exponential, 4=exponential

  double aux, c;
  switch (this->damageLaw){
  case 1: { // damage law used by Miehe phase-field model
    answer = gf/2. * (1./sqrt(1.-damage) - 1.);
  }break;
  case 2: {// damage law used by Fremond-Nedjar
    answer = gf*(1.+2.*(1.-beta)*damage/(1.-damage));
  }  break;
  case 3: {// pseudo-exponential law, not successful
    aux = 1.-(1.-beta)*log(1.-damage);
    answer = gf*aux*aux;
  }  break;
  case 4: {// exponential law, implicit
    c = beta/(1.-beta);
    aux = 1.+solveExpLaw(damage,c);
    answer = gf*aux*aux;
  }  break;
  case 5: {//Wu
    double dissPrime1inGamma = compute_dissipation_Wu_prime1_in_gamma(damageDrivingVariable);
    double damagePrime1inGamma;
    this->computeDamagePrime(damagePrime1inGamma, damageDrivingVariable, gp);     
    answer = dissPrime1inGamma/damagePrime1inGamma;
  }  break;    
  default: {// linear softening
    aux = 1. + ( beta - 1. )*damage;
    answer = gf/(aux*aux);
  }
  }
}

void
VarBasedDamageMaterial :: computeDissipationFunctionPrime2(double &answer, double damage, double damageDrivingVariable, GaussPoint *gp)
{
  /// softening - currently using 0=linear (default), 1=Miehe, 2=Fremond, 3=pseudo-exponential, 4=exponential
  
  double aux, c;
  switch (this->damageLaw){
  case 1: {// damage law used by Miehe phase-field model
    answer = gf/4. * 1./pow( 1- damage, 3./2. );
  }  break;
  case 2: {// damage law used by Fremond-Nedjar
    answer = 2.*gf*(1.-beta)/(1.-damage)/(1.-damage);
  }  break;
  case 3: {// pseudo-exponential law, not successful
    aux = 1.-(1.-beta)*log(1.-damage);
    answer = 2.*gf*aux*(1.-beta)/(1.-damage);
  }  break;
  case 4: {// exponential law, implicit
    c = beta/(1.-beta);
    aux = 1.+solveExpLaw(damage,c);
    answer = 2.*gf*aux*aux/(1.-damage+c*exp(-c*(aux-1.)));
  }  break;
  case 5: {//Wu
    double dissPrime1inGamma = compute_dissipation_Wu_prime1_in_gamma(damageDrivingVariable);
    double dissPrime2inGamma = compute_dissipation_Wu_prime2_in_gamma(damageDrivingVariable);
    double damagePrime1inGamma;
    double damagePrime2inGamma;
    this->computeDamagePrime(damagePrime1inGamma, damageDrivingVariable, gp);
    this->computeDamagePrime2(damagePrime2inGamma, damageDrivingVariable, gp);
    answer = dissPrime2inGamma/pow(damagePrime1inGamma,2) - (dissPrime1inGamma*damagePrime2inGamma)/pow(damagePrime1inGamma,3);
  }  break;    
  default: {// linear softening
    answer = 2. * (1. - beta) * gf/(pow( 1. + ( beta - 1. ) * damage, 3. ));
  }
  }
}

void
VarBasedDamageMaterial :: giveRealStressVectorGradientDamage(FloatArray &stress, double &localDamageDrivingVariable, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalDamageDrivingVariable, TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();

    FloatArray reducedTotalStrainVector;

    FloatMatrix de;
    double damage = 0.0;
    this->initTempStatus(gp);
    
    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, tStep, VM_Total);
    status->letTempStrainVectorBe(reducedTotalStrainVector);
    lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
    stress.beProductOf(de, reducedTotalStrainVector);
    // store the effective stress into the status
    status->letTempEffectiveStressVectorBe(stress);
    // compute equivalent strain
    this->computeLocalDamageDrivingVariable(localDamageDrivingVariable, gp, tStep);
    // compute damage from the nonlocal field
    this->computeDamage(damage, nonlocalDamageDrivingVariable, gp);
    stress.times(1-damage);
    // update gp
    status->letTempStressVectorBe(stress);
    status->setTempNonlocalDamageDrivingVariable(nonlocalDamageDrivingVariable);
    status->setTempLocalDamageDrivingVariable(localDamageDrivingVariable);
    status->setTempDamage(damage);

#ifdef keep_track_of_dissipated_energy
    status->computeWork(gp);
#endif

    
}




  
void
VarBasedDamageMaterial :: computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep)
{

  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus* >( gp->giveMaterialStatus() );
  FloatArray strain = status->giveTempStrainVector();
  FloatArray redStrain;
  StructuralMaterial :: giveReducedSymVectorForm(redStrain, strain, gp->giveMaterialMode());
  double eps, E;
  E = linearElasticMaterial->give('E', gp);
  this->computeEquivalentStrain(eps, redStrain, gp, tStep);
  answer = 0.5*E*eps*eps;    
}


MaterialStatus *
VarBasedDamageMaterial :: CreateStatus(GaussPoint *gp) const
{
  return new VarBasedDamageMaterialStatus(1, VarBasedDamageMaterial :: domain, gp, initialDamage);
}




#ifdef keep_track_of_dissipated_energy
void
VarBasedDamageMaterial :: computeRegulirizingWork(GaussPoint *gp, const FloatArray &nonlocalDamageDrivingVariableGrad)
{
    double tempReg = 0.5*gf*internalLength*internalLength*nonlocalDamageDrivingVariableGrad.computeSquaredNorm();
    VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus* >( gp->giveMaterialStatus() );
    status->setTempRegularizingEnergy(tempReg);
}
#endif

int
VarBasedDamageMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

  
    VarBasedDamageMaterialStatus* status = static_cast< VarBasedDamageMaterialStatus* >( this->giveStatus(gp) );
    if ( type == IST_RegularizingEnergy ) {
        answer.resize(1);
        answer.at(1) = status->giveRegularizingEnergy();
    }  else {
        return IsotropicDamageMaterial1 :: giveIPValue(answer, gp, type, tStep);
    }

    return 1; // to make the compiler happy
}

  
  VarBasedDamageMaterialStatus :: VarBasedDamageMaterialStatus(int n, Domain *d, GaussPoint *g, double initialDamage) : IsotropicDamageMaterial1Status(n, d, g), GradientDamageMaterialStatusExtensionInterface()
{
  damage = initialDamage;
  int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
  effectiveStressVector.resize(rsize);
  tempEffectiveStressVector.resize(rsize);
  strainEnergy = 0;
  tempRegularizingEnergy = regularizingEnergy = 0;
}


VarBasedDamageMaterialStatus :: ~VarBasedDamageMaterialStatus()
{

}



void
VarBasedDamageMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    IsotropicDamageMaterial1Status :: initTempStatus();
    GradientDamageMaterialStatusExtensionInterface :: initTempStatus();
    tempEffectiveStressVector = effectiveStressVector;
    tempStrainEnergy = strainEnergy;
    this->tempDamage = this->damage;
    tempRegularizingEnergy = regularizingEnergy;
}





  

void
VarBasedDamageMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    IsotropicDamageMaterial1Status :: updateYourself(tStep);
    GradientDamageMaterialStatusExtensionInterface :: updateYourself(tStep);
    effectiveStressVector = tempEffectiveStressVector;
    strainEnergy =  tempStrainEnergy;
    regularizingEnergy = tempRegularizingEnergy;
      
}



contextIOResultType
VarBasedDamageMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = IsotropicDamageMaterial1Status :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream.write(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
VarBasedDamageMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = IsotropicDamageMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream.read(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}



}     // end namespace oofem

