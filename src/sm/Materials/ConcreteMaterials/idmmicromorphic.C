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
#include "engngm.h"

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
    IR_GIVE_FIELD(ir, k1, _IFT_IsotropicDamageMaterialMicromorphic_k1);
    k2 = 2.*gf*internalLength*internalLength;
    return IRRT_OK;
}




void
IsotropicDamageMaterialMicromorphic :: giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
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
  
  if ( mode == TangentStiffness ) {
    if(tempDamage > status->giveDamage() ) {
      double localDamageDrivingVariable, dDamage, ddDamage, dDiss, ddDiss, equivStrain, storedEnergy, E;
      FloatArray reducedStrain, eta, effectiveStress;
      FloatMatrix correctionTerm;


      FloatArray totalStrain = status->giveTempStrainVector();
      StructuralMaterial :: giveReducedSymVectorForm( reducedStrain, totalStrain, gp->giveMaterialMode() );

      localDamageDrivingVariable = status->giveTempLocalDamageDrivingVariable();
      this->computeDissipationFunctionPrime(dDiss, tempDamage, localDamageDrivingVariable, gp);
      this->computeDissipationFunctionPrime2(ddDiss, tempDamage, localDamageDrivingVariable, gp);
      this->computeEquivalentStrain(equivStrain, reducedStrain, gp, tStep);      
      this->computeDamagePrime(dDamage, localDamageDrivingVariable, gp);
      this->computeDamagePrime2(ddDamage, localDamageDrivingVariable, gp);
      this->computeEta(eta, reducedStrain, gp, tStep);
      effectiveStress = status->giveTempStressVector();
      effectiveStress.times(1. / ( 1. - tempDamage ));
      E = linearElasticMaterial->give('E', gp);
      storedEnergy = 0.5 * E * equivStrain * equivStrain;
      // dyadic product of eff stress and eta
      correctionTerm.beDyadicProductOf(effectiveStress, eta);

      // times minus derivative of damage function
      correctionTerm.times( -  E * equivStrain * dDamage * dDamage / ( k1 + dDamage * dDamage * ddDiss + ddDamage * ( dDiss - storedEnergy ) ) );
      // add to secant stiffness
      answer.add(correctionTerm);


      
      /*
      FloatArray strainP, stressP, oldStrain, oldStress;
      oldStress = status->giveTempStressVector();
      oldStrain = status->giveTempStrainVector();
      strainP = oldStrain;
      double pert = 1.e-6 * strainP.at(1);
      strainP.at(1) += pert;
      double lddv;
      double mddv = status->giveTempNonlocalDamageDrivingVariable();
      this->giveRealStressVectorGradientDamage(stressP, lddv, gp, strainP, mddv, tStep);
      FloatMatrix stiff;
      stiff.resize(1,1);
      stiff.at(1,1) = (stressP.at(1) - oldStress.at(1))/pert;
      this->giveRealStressVectorGradientDamage(stressP, lddv, gp, oldStrain, mddv, tStep);
      */
    }
  }
  

}

  



void
IsotropicDamageMaterialMicromorphic :: giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
  answer.initFromVector(status->giveTempStressVector(), false);

  if ( mode == TangentStiffness ) {

    double tempDamage = status->giveTempDamage();
    if(tempDamage > status->giveDamage() ) {
      answer.times( 1. / (1. - tempDamage ) );
      double dDiss, ddDiss, dDamage, ddDamage, localDamageDrivingVariable, equivStrain;
      FloatArray totalStrain, reducedStrain;
      
      totalStrain = status->giveTempStrainVector();
      StructuralMaterial :: giveReducedSymVectorForm( reducedStrain, totalStrain, gp->giveMaterialMode() );
      this->computeEquivalentStrain(equivStrain, reducedStrain, gp, tStep);
      double E = linearElasticMaterial->give('E', gp);
      double storedEnergy = 0.5 * E * equivStrain * equivStrain;

      localDamageDrivingVariable = status->giveTempLocalDamageDrivingVariable();
      this->computeDissipationFunctionPrime(dDiss, tempDamage, localDamageDrivingVariable, gp);
      this->computeDissipationFunctionPrime2(ddDiss, tempDamage, localDamageDrivingVariable, gp);
      

      this->computeDamagePrime(dDamage, localDamageDrivingVariable, gp);
      this->computeDamagePrime2(ddDamage, localDamageDrivingVariable, gp);
      answer.times( - k1 * dDamage / ( k1 + dDamage * dDamage * ddDiss + ddDamage * ( dDiss - storedEnergy ) ) );


      /*
      double lddv;
      FloatArray strainP, stressP, oldStrain, oldStress;
      oldStress = status->giveTempStressVector();
      oldStrain = status->giveTempStrainVector();
      double mddv = status->giveTempNonlocalDamageDrivingVariable();
      double pert = 1.e-6 * mddv;
      strainP = oldStrain;
      double mddvp = mddv + pert;
      this->giveRealStressVectorGradientDamage(stressP, lddv, gp, strainP, mddvp, tStep);
      FloatMatrix stiff;
      stiff.resize(1,1);
      stiff.at(1,1) = (stressP.at(1) -oldStress.at(1))/pert;
      this->giveRealStressVectorGradientDamage(stressP, lddv, gp, oldStrain, mddv, tStep);
      */
    } else {
      answer.times(0);
    }
  } else {
    answer.times(0);
  }

}  





void
IsotropicDamageMaterialMicromorphic :: giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
  double equivStrain;
  FloatArray reducedStrain, eta;
  FloatArray totalStrain = status->giveTempStrainVector();
  StructuralMaterial :: giveReducedSymVectorForm( reducedStrain, totalStrain, gp->giveMaterialMode() );
  double tempDamage = status->giveTempDamage();
  this->computeEquivalentStrain(equivStrain, reducedStrain, gp, tStep);
  this->computeEta(eta, reducedStrain, gp, tStep);
  answer.initFromVector(eta, false);
  if ( mode == TangentStiffness ) {
    if(tempDamage > status->giveDamage() ) {
      double dDamage, ddDamage, dDiss, ddDiss;
      double E = linearElasticMaterial->give('E', gp);
      double localDamageDrivingVariable = status->giveTempLocalDamageDrivingVariable();
      double storedEnergy = 0.5 * E * equivStrain * equivStrain;
   
      this->computeDissipationFunctionPrime(dDiss, tempDamage, localDamageDrivingVariable, gp);
      this->computeDissipationFunctionPrime2(ddDiss, tempDamage, localDamageDrivingVariable, gp);
      this->computeDamagePrime(dDamage, localDamageDrivingVariable, gp);
      this->computeDamagePrime2(ddDamage, localDamageDrivingVariable, gp);
      // times minus derivative of damage function
      answer.times( - k1 * E * equivStrain * dDamage / ( k1 + dDamage * dDamage * ddDiss + ddDamage * ( dDiss - storedEnergy ) ) );

      /*
      FloatArray strainP, stressP, oldStrain, oldStress;
      oldStress = status->giveTempStressVector();
      oldStrain = status->giveTempStrainVector();
      double mddv = status->giveTempNonlocalDamageDrivingVariable();
      double lddv = status->giveTempLocalDamageDrivingVariable();
      double lddvp;
      strainP = oldStrain;
      double pert = 1.e-6 * strainP.at(1);
      strainP.at(1) += pert;
      this->giveRealStressVectorGradientDamage(stressP, lddvp, gp, strainP, mddv, tStep);
      FloatMatrix stiff;
      stiff.resize(1,1);
      stiff.at(1,1) = -k1 * (lddvp - lddv) / pert;
      this->giveRealStressVectorGradientDamage(stressP, lddv, gp, oldStrain, mddv, tStep);
      */
      
    } else {
      answer.times(0);
    }
  } else {
    answer.times(0);
  }  
}

  


void
IsotropicDamageMaterialMicromorphic :: giveGradientDamageStiffnessMatrix_dd_NN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

  answer.resize(1,1);
  answer.at(1,1) = k1;  

  if ( mode == TangentStiffness ) {
    VarBasedDamageMaterialStatus *status = static_cast< VarBasedDamageMaterialStatus * >( this->giveStatus(gp) );
    double tempDamage = status->giveTempDamage();
    if(tempDamage > status->giveDamage() ) {
      double dDiss, ddDiss, dDamage, ddDamage, localDamageDrivingVariable, equivStrain;
      FloatArray totalStrain, reducedStrain;

      localDamageDrivingVariable = status->giveTempLocalDamageDrivingVariable();      
      totalStrain = status->giveTempStrainVector();
      StructuralMaterial :: giveReducedSymVectorForm( reducedStrain, totalStrain, gp->giveMaterialMode() );
      this->computeEquivalentStrain(equivStrain, reducedStrain, gp, tStep);
      double E = linearElasticMaterial->give('E', gp);
      double storedEnergy = 0.5 * E * equivStrain * equivStrain;

      
      this->computeDissipationFunctionPrime(dDiss, tempDamage, localDamageDrivingVariable, gp);
      this->computeDissipationFunctionPrime2(ddDiss, tempDamage, localDamageDrivingVariable, gp);
      this->computeDamagePrime(dDamage, localDamageDrivingVariable, gp);
      this->computeDamagePrime2(ddDamage, localDamageDrivingVariable, gp);      
      answer.at(1,1) *= ( dDamage * dDamage * ddDiss + ddDamage * ( dDiss - storedEnergy ) ) / ( k1 + dDamage * dDamage * ddDiss + ddDamage * ( dDiss - storedEnergy ) );

      /*
      FloatArray strainP, stressP, oldStrain, oldStress;
      oldStress = status->giveTempStressVector();
      oldStrain = status->giveTempStrainVector();
      double mddv = status->giveTempNonlocalDamageDrivingVariable();
      double lddv = status->giveTempLocalDamageDrivingVariable();
      double lddvp;
      strainP = oldStrain;
      double pert = 1.e-6 * mddv;
      double mddvp = mddv + pert;
      this->giveRealStressVectorGradientDamage(stressP, lddvp, gp, strainP, mddvp, tStep);
      FloatMatrix stiff;
      stiff.resize(1,1);
      stiff.at(1,1) = k1 * (1.- (lddvp - lddv) / pert );
      this->giveRealStressVectorGradientDamage(stressP, lddv, gp, oldStrain, mddv, tStep);
      */
      
    }
  } else {
    answer.times(0);
  }  

}


  
void
IsotropicDamageMaterialMicromorphic :: giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
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
#ifdef keep_track_of_dissipated_energy
  this->computeRegulirizingWork(gp,  micromorphicDamageDrivingVariableGrad);
#endif
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
    double epsEq, E, storedEnergy;
    E = linearElasticMaterial->give('E', gp);
    this->computeEquivalentStrain(epsEq, strain, gp, tStep);
    storedEnergy = 0.5*E*epsEq*epsEq;    
    //damage = status->giveDamage();
    this->computeDamagePrime(dDamage, damageDrivingVariable, gp);
    this->computeDissipationFunctionPrime(dDiss, damage, damageDrivingVariable, gp);
    /*damage based version*/
    //    f = storedEnergy + k1 * ( micromorphicDamage - damage) - dDiss;
    this->computeEquivalentStrain(epsEq, strain, gp, tStep);
    storedEnergy = 0.5 * E * epsEq * epsEq;    
    f = storedEnergy * dDamage + k1 * ( micromorphicDamageDrivingVariable - damageDrivingVariable) - dDiss * dDamage;
    if(f <= 0) {
      tempDamageDrivingVariable = damageDrivingVariable;
      tempDamage = damage;
    } else {
      // compute damage from the nonlocal field
      //      this->computeDamage(tempDamage, damage, micromorphicDamage, storedEnergy, gp);
      this->computeDamageDrivingVariable(tempDamageDrivingVariable, damageDrivingVariable, micromorphicDamageDrivingVariable, storedEnergy, gp);
      this->computeDamage(tempDamage, tempDamageDrivingVariable, gp);    

    }
    
    stress.times( 1. - tempDamage );
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
    this->computeDissipationFunctionPrime(dDiss,damage, localDamageDrivingVariable_n, gp);
    this->computeDissipationFunctionPrime2(ddDiss,damage, localDamageDrivingVariable_n, gp);
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
      //OOFEM_WARNING("Error in compute damage function");
      this->giveDomain()->giveEngngModel()->setAnalysisCrash(true);
      break;
    }

  }
  
}
  

  


MaterialStatus *
IsotropicDamageMaterialMicromorphic :: CreateStatus(GaussPoint *gp) const
{
    return new IsotropicDamageMaterialMicromorphicStatus(1, IsotropicDamageMaterialMicromorphic :: domain, gp);
}


#ifdef keep_track_of_dissipated_energy
void
IsotropicDamageMaterialMicromorphic :: computeRegulirizingWork(GaussPoint *gp, const FloatArray &nonlocalDamageDrivingVariableGrad)
{
    IsotropicDamageMaterialMicromorphicStatus *status = static_cast< IsotropicDamageMaterialMicromorphicStatus* >( gp->giveMaterialStatus() );
    double nlddv = status->giveTempNonlocalDamageDrivingVariable();
    double lddv = status->giveTempLocalDamageDrivingVariable();
      
    double tempReg = 0.5 * k2 * nonlocalDamageDrivingVariableGrad.computeSquaredNorm() + 0.5 * k1 * (nlddv - lddv) * (nlddv - lddv);
    
    status->setTempRegularizingEnergy(tempReg);
} 
#endif



  
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
