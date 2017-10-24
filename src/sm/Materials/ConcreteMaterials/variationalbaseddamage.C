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

VarBasedDamageMaterial :: VarBasedDamageMaterial(int n, Domain *d) : IsotropicDamageMaterial1(n, d), GradientDamageMaterialExtensionInterface(d)
    //
    // constructor
    //
{
    internalLength = 0;
    initialDamage = 0;
}


VarBasedDamageMaterial :: ~VarBasedDamageMaterial()
//
// destructor
//
{ }

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
    IR_GIVE_OPTIONAL_FIELD(ir, initialDamage, _IFT_VarBasedDamageMaterial_initDamage);


    
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
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMat:
        give1dStressStiffMtrx(answer, mode, gp, tStep);
        break;
    case _PlaneStress:
        givePlaneStressStiffMtrx(answer, mode, gp, tStep);
        break;
    case _PlaneStrain:
        givePlaneStrainStiffMtrx(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("mMode = %d not supported\n", mMode);
    }
}

  

void
VarBasedDamageMaterial :: giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );

  FloatArray stress = status->giveTempEffectiveStressVector();
  double damage = status->giveTempDamage();
  double dDamage = damage - status->giveDamage();

    answer = stress;
    answer.times(-2.*(1-damage));
    //@todo: experimental staff
    answer.times(0.);

 
}

void
VarBasedDamageMaterial :: giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );

  FloatArray stress = status->giveTempEffectiveStressVector();
  double damage = status->giveTempDamage();
  double dDamage = damage - status->giveDamage();

    answer = stress;
    answer.times(-2.*(1-damage));
    //@todo: experimental staff
    answer.times(0.);

  /*
  FloatArray s, os, r;
  double eps = 1.e-8;
  double ld;
  FloatArray totalStrain = status->giveTempStrainVector();
  os = status->giveTempStressVector();
  double nonlocalDamageDrivingVariable = status->giveNonlocalDamageDrivingVariable();
  double dnddv = nonlocalDamageDrivingVariable + eps;  
  this->giveRealStressVectorGradientDamage(s, ld, gp, totalStrain, dnddv, tStep);
  r = s;
  r.subtract(os);
  r.times(1./eps);

  this->giveRealStressVectorGradientDamage(s, ld, gp, totalStrain, nonlocalDamageDrivingVariable, tStep);

  */



  
}

void
VarBasedDamageMaterial :: giveGradientDamageStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{  
  double localDamageDrivingVariable;
  this-> computeLocalDamageDrivingVariable(localDamageDrivingVariable, gp, tStep);
  
  answer.resize(1,1);
  double factor = gf/internalLength + localDamageDrivingVariable;
  // regularization
  double reg;
  this->computeStiffnessRegularizationTerm(reg,gp,tStep);
  factor +=  reg;  
  answer.at(1,1) = (factor);  
}

void
VarBasedDamageMaterial :: giveGradientDamageStiffnessMatrix_dd_l(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  //answer = {{1, 0, 0,0 },{0, 1,0, 0},{0, 0, 1,0}, {0, 0, 0,1}};
  answer = {{1, 0},{0, 1}};
  answer.times(gf*internalLength);
  
}

  

void
VarBasedDamageMaterial :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp,  TimeStep *tStep)
{
    IsotropicDamageMaterialStatus *status = static_cast< IsotropicDamageMaterialStatus * >( this->giveStatus(gp) );
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    double om;
    om = status->giveTempDamage();
    om = min(om, maxOmega);
    answer.resize(1, 1);
    answer.at(1, 1) = lmat->give('E', gp);
    answer.times(1.0 - om);
}




void
VarBasedDamageMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
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
    //@todo change this
    answer.times(1.0 - tempDamage);
}



void
VarBasedDamageMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
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
    //@todo change this
    answer.times(1.0 - tempDamage);
}


void
VarBasedDamageMaterial :: giveNonlocalInternalForces_N_factor(double &answer,GaussPoint *gp, TimeStep *tStep)
{
  
   double localDamageDrivingVariable;
   this-> computeLocalDamageDrivingVariable(localDamageDrivingVariable, gp, tStep);
   answer = gf/internalLength + localDamageDrivingVariable;
 
}

  void
VarBasedDamageMaterial :: giveNonlocalInternalForces_B_factor(double &answer,GaussPoint *gp, TimeStep *tStep)
{
  answer = gf * internalLength;
}

void  
VarBasedDamageMaterial :: computeInternalForcesRegularizationTerm(double &answer,GaussPoint *gp, TimeStep *tStep)
{

    VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
    double dTime = tStep->giveTimeIncrement();
    double damage = status->giveDamage();
    answer = status->giveTempDamage() - damage;
    if(answer < 0) {
      answer *= answer;
      answer *= -penalty/dTime/2;
    } else {
      answer = 0;
    }   
}

  void  
VarBasedDamageMaterial :: computeStiffnessRegularizationTerm(double &answer,GaussPoint *gp, TimeStep *tStep)
{

    VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
    double dTime = tStep->giveTimeIncrement();
    double damage = status->giveDamage();
    answer = status->giveTempDamage() - damage;
    if(answer < 0) {
      answer *= -penalty/dTime;
    } else {
      answer = 0;
    }   
}
  

void
VarBasedDamageMaterial :: computeDamage(double &answer, double damageDrivingVariable, GaussPoint *gp)
{

  double damage;
  damage = damageDrivingVariable;
  if(damage > maxOmega) {
    damage = maxOmega;
  } else if(damage < 0) {
    damage = 0;
    /*    VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
	  damage = status->giveDamage();*/
  } else {
    answer = damage;
  }
  /*
  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
  if(damage < status->giveDamage()) {

    damage = status->giveDamage();
  }

  if(damage != damage)
    int ahoj = 1;
  */
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
    if(initialDamage != 0)
      int ahoj = 0;
    stress.times(1-damage);
    // @todo: remove this
    stress.times(1-damage);
    if(stress.at(1) != stress.at(1))
      int ahoj = 1;


    // update gp
    status->letTempStressVectorBe(stress);
    status->setNonlocalDamageDrivingVariable(nonlocalDamageDrivingVariable);
    status->setTempLocalDamageDrivingVariable(localDamageDrivingVariable);
    status->setTempDamage(damage);
}
void
VarBasedDamageMaterial :: computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep)
{
  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus* >( gp->giveMaterialStatus() );
    FloatArray strain, stress;
    stress = status->giveTempEffectiveStressVector();
    strain = status->giveTempStrainVector();
    answer = stress.dotProduct( strain );
    double strainEnergy = status->giveStrainEnergy();
    if(answer < strainEnergy) {
      answer = strainEnergy;
    }
    status->setTempStrainEnergy(answer);
}


MaterialStatus *
VarBasedDamageMaterial :: CreateStatus(GaussPoint *gp) const
{
  return new VarBasedDamageMaterialStatus(1, VarBasedDamageMaterial :: domain, gp, initialDamage);
}


  VarBasedDamageMaterialStatus :: VarBasedDamageMaterialStatus(int n, Domain *d, GaussPoint *g, double initialDamage) : IsotropicDamageMaterial1Status(n, d, g), GradientDamageMaterialStatusExtensionInterface()
{
  damage = initialDamage;
  int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
  effectiveStressVector.resize(rsize);
  tempEffectiveStressVector.resize(rsize);
  strainEnergy = 0;
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
