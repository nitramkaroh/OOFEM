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

#include "idmmicromorphic.h"
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
REGISTER_Material(IsotropicDamageMaterialMicromorphic);

IsotropicDamageMaterialMicromorphic :: IsotropicDamageMaterialMicromorphic(int n, Domain *d) : VarBasedDamageMaterial(n, d)
//
// constructor
//
{

}


IsotropicDamageMaterialMicromorphic :: ~IsotropicDamageMaterialMicromorphic()
//
// destructor
//
{ }

IRResultType
IsotropicDamageMaterialMicromorphic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    result = VarBasedDamageMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    //    IR_GIVE_OPTIONAL_FIELD(ir, di_eta, _IFT_IsotropicDamageMaterialMicromorphic_di_eta);
    IR_GIVE_FIELD(ir, k1, _IFT_IsotropicDamageMaterialMicromorphic_k);
    k2 = 2.*gf*internalLength*internalLength;
    return IRRT_OK;
}









void
IsotropicDamageMaterialMicromorphic :: giveGradientDamageStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  answer.resize(1,1);
  answer.at(1,1) = k1;  
}


  
void
IsotropicDamageMaterialMicromorphic :: giveGradientDamageStiffnessMatrix_dd_l(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode matMode = gp->giveMaterialMode();
  if(matMode == _1dMat) {
    answer = {{1}};
  } else if(matMode == _PlaneStress || matMode == _PlaneStrain) {
    answer = {{1, 0},{0, 1}};
  } else if(matMode == _3dMat) {
    answer = {{1, 0, 0},{0, 1, 0}, {0, 0, 1}};
  }
  answer.times(k2);  
}



void
IsotropicDamageMaterialMicromorphic :: giveNonlocalInternalForces_N_factor(double &answer, double micromorphicDamageDrivingVariable, GaussPoint *gp, TimeStep *tStep)
{

  IsotropicDamageMaterialMicromorphicStatus *status = static_cast< IsotropicDamageMaterialMicromorphicStatus * >( this->giveStatus(gp) );

  double damageDrivingVariable = status->giveTempLocalDamageDrivingVariable();

  answer = k1 * ( micromorphicDamageDrivingVariable - damageDrivingVariable);
 
}

void
IsotropicDamageMaterialMicromorphic :: giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &micromorphicDamageDrivingVariableGrad, GaussPoint *gp, TimeStep *tStep)
{
  answer = micromorphicDamageDrivingVariableGrad;
  answer.times(k2);
}


  
  
void
IsotropicDamageMaterialMicromorphic :: giveRealStressVectorGradientDamage(FloatArray &stress, double &tempDamageDrivingVariable, GaussPoint *gp, const FloatArray &totalStrain, double micromorphicDamageDrivingVariable, TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    IsotropicDamageMaterialMicromorphicStatus *status = static_cast< IsotropicDamageMaterialMicromorphicStatus * >( this->giveStatus(gp) );
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    FloatArray strain;
    
    FloatMatrix de;
    double damage, dDamage, f, dDiss, damageDrivingVariable,tempDamage;
    //    this->initTempStatus(gp);
    
    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(strain, gp, totalStrain, tStep, VM_Total);
    status->letTempStrainVectorBe(strain);
    lmat->giveStiffnessMatrix(de, SecantStiffness, gp, tStep);
    stress.beProductOf(de, strain);
    // compute equivalent strain
    damageDrivingVariable = status->giveLocalDamageDrivingVariable();
    // compute damage and micromorphic damage
    this->computeDamage(damage, damageDrivingVariable, gp);    
    // compute equivalent strain
    // compute equivalent strain
    double epsEq, E, storedEnergy;
    E = linearElasticMaterial->give('E', gp);
    this->computeEquivalentStrain(epsEq, strain, gp, tStep);
    storedEnergy = 0.5*E*epsEq*epsEq;    
    //damage = status->giveDamage();
    this->computeDamagePrime(dDamage, damageDrivingVariable, gp);
    this->computeDissipationFunctionPrime(dDiss, damage, gp);
    /*damage based version*/
    //    f = storedEnergy + k1 * ( micromorphicDamage - damage) - dDiss;
    this->computeEquivalentStrain(epsEq, strain, gp, tStep);
    storedEnergy = 0.5*E*epsEq*epsEq;    
    f = storedEnergy*dDamage + k1 * ( micromorphicDamageDrivingVariable - damageDrivingVariable) - dDiss*dDamage;
    if(f <= 0) {
      tempDamageDrivingVariable = damageDrivingVariable;
      tempDamage = damage;
    } else {
      // compute damage from the nonlocal field
      //      this->computeDamage(tempDamage, damage, micromorphicDamage, storedEnergy, gp);
      this->computeDamageDrivingVariable(tempDamageDrivingVariable, damageDrivingVariable, micromorphicDamageDrivingVariable, storedEnergy, gp);
      this->computeDamage(tempDamage, tempDamageDrivingVariable, gp);    

    }
    
    stress.times(1-tempDamage);
    // update gp
    status->letTempStressVectorBe(stress);
    status->setTempLocalDamageDrivingVariable(tempDamageDrivingVariable);
    status->setTempNonlocalDamageDrivingVariable(micromorphicDamageDrivingVariable);
    status->setTempDamage(tempDamage);
    //    status->setTempMicromorphicDamage(micromorphicDamage);

#ifdef keep_track_of_dissipated_energy
    status->computeWork(gp);
#endif

}


  /*
void
IsotropicDamageMaterialMicromorphic :: computeDamage(double &answer, double damage_n, double micromorphicDamage, double storedEnergy, GaussPoint *gp)
{

  int index  = 0;
  double f, df, dDiss, ddDiss;
  answer = damage_n;
  while(true) {
    index++;
    this->computeDissipationFunctionPrime(dDiss,answer, gp);
    f = storedEnergy + k1 * ( micromorphicDamage - answer) - dDiss;
    this->computeDissipationFunctionPrime2(ddDiss, answer, gp);
    df = -k1 - ddDiss;
    answer -= f/df;
    if(answer > 1.) {
      answer = 1.;
    } else if(answer < 0.) {
      answer = 0.;
    }
    
    if(fabs(f) < 1.e-8*k1) {     
      break;
    }
    if(index > 1000) {
      OOFEM_ERROR("Error in compute damage function");
    }

  }
  
}
  */

void
IsotropicDamageMaterialMicromorphic :: computeDamageDrivingVariable(double &answer, double localDamageDrivingVariable_n, double micromorphicDamageDrivingVariable, double storedEnergy, GaussPoint *gp)
{

  int index  = 0;
  double f, df, dDiss, ddDiss, dDamage, ddDamage, damage;
  answer = localDamageDrivingVariable_n;
  this->computeDamage(damage, answer, gp);    
  
  while(true) {
    index++;
    this->computeDissipationFunctionPrime(dDiss,damage, gp);
    this->computeDissipationFunctionPrime2(ddDiss,damage, gp);
    this->computeDamagePrime(dDamage, answer, gp);
    this->computeDamagePrime2(ddDamage, answer, gp);
    f = dDamage * storedEnergy + k1 * ( micromorphicDamageDrivingVariable - answer) - dDamage * dDiss;  
    df = storedEnergy * ddDamage - ddDiss * dDamage * dDamage - dDiss * ddDamage -k1;
    answer -= f/df;
    this->computeDamage(damage, answer, gp);    

    
    if(fabs(f) < 1.e-8*k1) {     
      break;
    }
    if(index > 1000) {
      OOFEM_ERROR("Error in compute damage function");
    }

  }
  
}
  

  


MaterialStatus *
IsotropicDamageMaterialMicromorphic :: CreateStatus(GaussPoint *gp) const
{
    return new IsotropicDamageMaterialMicromorphicStatus(1, IsotropicDamageMaterialMicromorphic :: domain, gp);
}


  IsotropicDamageMaterialMicromorphicStatus :: IsotropicDamageMaterialMicromorphicStatus(int n, Domain *d, GaussPoint *g) : VarBasedDamageMaterialStatus(n, d, g, 0)
{

}


IsotropicDamageMaterialMicromorphicStatus :: ~IsotropicDamageMaterialMicromorphicStatus()
{ }



  
void
IsotropicDamageMaterialMicromorphicStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    VarBasedDamageMaterialStatus :: initTempStatus();
}

  
void
IsotropicDamageMaterialMicromorphicStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    VarBasedDamageMaterialStatus :: updateYourself(tStep);
}
 

  /*
contextIOResultType
IsotropicDamageMaterialMicromorphicStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
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
  */
  /*
contextIOResultType
IsotropicDamageMaterialMicromorphicStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
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
  */

}     // end namespace oofem
