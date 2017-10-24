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

#include "mplasticdamagematerial2.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"



namespace oofem {
#define YIELD_TOL 1.e-6
#define YIELD_TOL_SECONDARY 1.e-4
#define RES_TOL   1.e-8
#define PLASTIC_MATERIAL_MAX_ITERATIONS 120

MPlasticDamageMaterial2 :: MPlasticDamageMaterial2(int n, Domain *d)  : MPlasticMaterial2(n, d)
    //
    // constructor
    //
{
}



MPlasticDamageMaterial2 :: ~MPlasticDamageMaterial2()
//
// destructor
//
{
}



MaterialStatus *
MPlasticDamageMaterial2 :: CreateStatus(GaussPoint *gp) const
{
    return new MPlasticDamageMaterial2Status(1, this->giveDomain(), gp, this->giveSizeOfReducedHardeningVarsVector(gp));
}


void
MPlasticDamageMaterial2 :: giveRealStressVector(FloatArray &answer,
                                          GaussPoint *gp,
                                          const FloatArray &totalStrain,
                                          TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
// completely formulated in Reduced stress-strain space
{
    FloatArray strainSpaceHardeningVariables;
    FloatArray fullStressVector;
    FloatArray strainVectorR, plasticStrainVectorR;
    FloatArray helpVec;
    IntArray activeConditionMap(this->nsurf);
    FloatArray gamma;

    MPlasticDamageMaterial2Status *status = static_cast< MPlasticDamageMaterial2Status * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    MPlasticDamageMaterial2 :: giveRealStressVector(answer, gp, totalStrain, tStep);
    //perform isotropic damage on effective stress
    double omega = computeDamage(gp, strainSpaceHardeningVariables, tStep);
    status->letTempDamageBe(omega);
    answer.times(1. - omega);
    //end of damage part
}


void
MPlasticDamageMaterial2 :: give_dLambda_dEps_Matrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    //
    // returns receiver material matrix for given reached state
    // described by temp variables.
    //

    FloatMatrix consistentModuli, elasticModuli, umat, vmat;
    FloatMatrix elasticModuliInverse, kl, ks;
    FloatMatrix gradientMatrix, gmat, gmatInv, gradMat, helpMtrx, helpMtrx2, answerR;
    FloatArray gradientVector, stressVector, fullStressVector;
    FloatArray strainSpaceHardeningVariables, helpVector;
    std :: vector< FloatArray > yieldGradSigVec(this->nsurf), loadGradSigVec(this->nsurf), * loadGradSigVecPtr;
    std :: vector< FloatArray > yieldGradKVec(this->nsurf), loadGradKVec(this->nsurf);
    FloatArray helpVector2;

    IntArray activeConditionMap, mask;
    FloatArray gamma;
    int strSize, indx, actSurf = 0;
    bool hasHardening = this->hasHardening() > 0;

    MPlasticDamageMaterial2Status *status = static_cast< MPlasticDamageMaterial2Status * >( this->giveStatus(gp) );

    if ( this->plType == associatedPT ) {
        loadGradSigVecPtr = & yieldGradSigVec;
    } else {
        loadGradSigVecPtr  = & loadGradSigVec;
    }

    // ask for plastic consistency parameter
    gamma = status->giveTempGamma();
    activeConditionMap = status->giveTempActiveConditionMap();
    //
    // check for elastic cases
    //
    if ( ( status->giveTempStateFlag() == MPlasticMaterial2Status :: PM_Elastic ) ||
        ( status->giveTempStateFlag() == MPlasticMaterial2Status :: PM_Unloading ) ) {
        this->giveStiffnessMatrix(answer, ElasticStiffness, gp, tStep);
        return;
    }

    //
    // plastic case
    //
    // determine number of active surfaces
    for ( int i = 1; i <= nsurf; i++ ) {
        if ( activeConditionMap.at(i) ) {
            actSurf++;
        }
    }

    // compute elastic moduli and its inverse
    this->computeReducedElasticModuli(elasticModuli, gp, tStep);
    elasticModuliInverse.beInverseOf(elasticModuli);
    strSize = elasticModuliInverse.giveNumberOfRows();

    stressVector = status->giveTempStressVector();
    StructuralMaterial :: giveFullSymVectorForm( fullStressVector, stressVector, gp->giveMaterialMode() );
    strainSpaceHardeningVariables = status->giveTempStrainSpaceHardeningVarsVector();

    //
    // compute consistent moduli
    //
    this->computeAlgorithmicModuli(consistentModuli, gp, elasticModuliInverse, gamma,
                                   activeConditionMap, fullStressVector, strainSpaceHardeningVariables);

    //computee gmatInv
    for ( int i = 1; i <= nsurf; i++ ) {
        this->computeReducedStressGradientVector(yieldGradSigVec [ i - 1 ], yieldFunction, i, gp, fullStressVector,
                                                 strainSpaceHardeningVariables);
        if ( hasHardening ) {
            this->computeKGradientVector(yieldGradKVec [ i - 1 ], yieldFunction, i, gp, fullStressVector,
                                         strainSpaceHardeningVariables);
        }
    }

    if ( this->plType == nonassociatedPT ) {
        for ( int i = 1; i <= nsurf; i++ ) {
            this->computeReducedStressGradientVector(loadGradSigVec [ i - 1 ], loadFunction, i, gp, fullStressVector,
                                                     strainSpaceHardeningVariables);
            if ( hasHardening ) {
                this->computeKGradientVector(loadGradKVec [ i - 1 ], loadFunction, i, gp, fullStressVector,
                                             strainSpaceHardeningVariables);
            }
        }
    }


    umat.resize(strSize, actSurf);
    umat.zero();
    vmat.resize(actSurf, strSize);
    vmat.zero();
    gmat.resize(actSurf, actSurf);
    gmat.zero();

    if ( hasHardening ) {
        this->computeReducedHardeningVarsSigmaGradient(ks, gp, activeConditionMap, fullStressVector,
                                                       strainSpaceHardeningVariables, gamma);
        this->computeReducedHardeningVarsLamGradient(kl, gp, actSurf, activeConditionMap,
                                                     fullStressVector, strainSpaceHardeningVariables, gamma);
    }

    for ( int i = 1; i <= nsurf; i++ ) {
        if ( ( indx = activeConditionMap.at(i) ) ) {
            if ( hasHardening ) {
                helpVector.beTProductOf(ks, yieldGradKVec [ i - 1 ]);
                helpVector.add(yieldGradSigVec [ i - 1 ]);
                for ( int j = 1; j <= strSize; j++ ) {
                    vmat.at(indx, j) = helpVector.at(j);
                }
            } else {
                for ( int j = 1; j <= strSize; j++ ) {
                    vmat.at(indx, j) = yieldGradSigVec [ i - 1 ].at(j);
                }
            }

            if ( hasHardening ) {
                this->computeReducedSKGradientMatrix(gradientMatrix,  i, gp, fullStressVector,
                                                     strainSpaceHardeningVariables);
                helpMtrx.beProductOf(gradientMatrix, kl);
                helpMtrx.times( gamma.at(i) );
                umat.add(helpMtrx);
            }

            for ( int j = 1; j <= strSize; j++ ) {
                umat.at(j, indx) += ( ( * loadGradSigVecPtr ) [ i - 1 ] ).at(j);
            }


            if ( hasHardening ) {
                helpVector2.beTProductOf(kl, yieldGradKVec [ i - 1 ]);
                for ( int j = 1; j <= actSurf; j++ ) {
                    gmat.at(indx, j) = ( -1.0 ) * helpVector2.at(j);
                }
            }
        }
    }

    helpMtrx.beProductOf(vmat, consistentModuli); // V S
    helpMtrx2.beProductOf(helpMtrx, umat);
    gmat.add(helpMtrx2);
    /////////////////////////////
    gmatInv.beInverseOf(gmat);


    answer.beProductOf(gmatInv, helpMtrx);
    answer.negated();
}




// overloaded from structural material

void
MPlasticDamageMaterial2 :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                   MatResponseMode mode,
                                                   GaussPoint *gp,
                                                   TimeStep *tStep)
//
//
//
// computes full 3d constitutive matrix for case of 3d stress-strain state.
// it returns elasto-plastic stiffness material matrix.
// if strainIncrement == NULL a loading is assumed
// for detailed description see (W.F.Chen: Plasticity in Reinforced Concrete, McGraw-Hill, 1982,
// chapter 6.)
//
// if derived material would like to implement failure behaviour
// it must redefine basic Give3dMaterialStiffnessMatrix function
// in order to take possible failure (tension cracking) into account
//
//
//
{ 
    MPlasticDamageMaterial2Status *status = static_cast< MPlasticDamageMaterial2Status * >( this->giveStatus(gp) );
    double tempDamage = status->giveTempDamage();
    
    MaterialMode originalMode = gp->giveMaterialMode();
    if ( originalMode != _3dMat ) {
        OOFEM_ERROR("Different stressStrain mode encountered");
    }

    // we can force 3d response, and we obtain correct 3d tangent matrix,
    // but in fact, stress integration algorithm will not work
    // because in stress integration algorithm we are unable to recognize
    // which reduction from 3d case should be performed to obtain correct result.
    // so for new stressStrain state, instead of programming 3d reduction,
    // you should enhance imposeConstraints functions for ne state, and
    // then programming simple inteface function for you stressstrain state
    // calling GiveMaterailStiffenssMatrix, which imposes constrains correctly.
    if ( mode == TangentStiffness ) {
        if ( rmType == mpm_ClosestPoint ) {
            this->giveConsistentStiffnessMatrix(answer, mode, gp, tStep);
        } else {
            this->giveElastoPlasticStiffnessMatrix(answer, mode, gp, tStep);
        }
	/*	answer.times(1.0 - tempDamage);

	/// add contribution of damage to the stiffness matrix
	FloatArray effectiveStress =  status->giveTempStressVector();
	effectiveStress.times(1./(1.0- tempDamage));
	// ask for plastic consistency parameter
	FloatArray gamma = status->giveTempGamma();
	// ask for hardening variables
	FloatArray strainSpaceHardeningVariables = status->giveTempStrainSpaceHardeningVarsVector();
	FloatMatrix dLam_dEps, dKappa_dEps, dKappa_dLam;	
	// compute derivative of damage driving variable wrt strain
	this->give_dLambda_dEps_Matrix(dLambda_dEps, _3dMat, gp, tStep);
	this->computeReducedHardeningVarsLamGradient(dKappa_dLam, gp, actSurf, activeConditionMap, effectiveStress, strainSpaceHardeningVars, gamma);
	dKappa_dEps.beProductOf(dKappa_dLam, dLam_dEps);
	// stiffness correction  = - omega_Prime * sigma  \otimes dKappa_dEps
	int size = this->giveSizeOfReducedHardeningVarsVector(gp);
	double damagePrime = this->give_dDamageParam_dKappa(strainSpaceHardeningVariables.at(size));
	FloatMatrix stiffnessCorrection;
	stiffnessCorrection.beDyadicProductOf(effectiveStress, dKappa_dEps);
	stiffnessCorrection.times(damagePrime);
	answer.subtract(stiffnessCorrection);
	*/
    } else if ( mode == SecantStiffness ) {
        this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
	answer.times(1.0 - tempDamage);
    } else {
        this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
    }
}


void
MPlasticDamageMaterial2 :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *tStep)

//
// return receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialStiffnessMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
//
{

    MPlasticDamageMaterial2Status *status = static_cast< MPlasticDamageMaterial2Status * >( this->giveStatus(gp) );
    double tempDamage = status->giveTempDamage(); 
  
    if ( mode == TangentStiffness ) {
        if ( rmType == mpm_ClosestPoint ) {
            this->giveConsistentStiffnessMatrix(answer, mode, gp, tStep);
        } else {
            this->giveElastoPlasticStiffnessMatrix(answer, mode, gp, tStep);
        }
	
	/*	answer.times(1.0 - tempDamage);
	/// add contribution of damage to the stiffness matrix
	FloatArray effectiveStress =  status->giveTempStressVector();
	effectiveStress.times(1./(1.0- tempDamage));
	// ask for plastic consistency parameter
	FloatArray gamma = status->giveTempGamma();
	// ask for hardening variables
	FloatArray strainSpaceHardeningVariables = status->giveTempStrainSpaceHardeningVarsVector();
	FloatMatrix dLam_dEps, dKappa_dEps, dKappa_dLam;	
	// compute derivative of damage driving variable wrt strain
	this->give_dLambda_dEps_Matrix(dLambda_dEps, _3dMat, gp, tStep);
	this->computeReducedHardeningVarsLamGradient(dKappa_dLam, gp, actSurf, activeConditionMap, effectiveStress, strainSpaceHardeningVars, gamma);
	dKappa_dEps.beProductOf(dKappa_dLam, dLam_dEps);
	// stiffness correction  = - omega_Prime * sigma  \otimes dKappa_dEps
	int size = this->giveSizeOfReducedHardeningVarsVector(gp);
	double damagePrime = this->give_dDamageParam_dKappa(strainSpaceHardeningVariables.at(size));
	FloatMatrix stiffnessCorrection;
	stiffnessCorrection.beDyadicProductOf(effectiveStress, dKappa_dEps);
	stiffnessCorrection.times(damagePrime);
	answer.subtract(stiffnessCorrection);*/

    } else if ( mode == SecantStiffness ) {
      this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
      answer.times(1.0 - tempDamage);
    } else {
        this->giveLinearElasticMaterial()->givePlaneStrainStiffMtrx(answer, mode, gp, tStep);
    }
}



  MPlasticDamageMaterial2Status :: MPlasticDamageMaterial2Status(int n, Domain *d, GaussPoint *g, int statusSize) : MPlasticMaterial2Status(n, d, g, statusSize), damage(0.), tempDamage(0.)
{ }

MPlasticDamageMaterial2Status :: ~MPlasticDamageMaterial2Status()
{ }

void
MPlasticDamageMaterial2Status :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( ( state_flag == MPlasticMaterial2Status :: PM_Yielding ) || ( state_flag == MPlasticMaterial2Status :: PM_Unloading ) ) {
        if ( state_flag == MPlasticMaterial2Status :: PM_Yielding ) {
            fprintf(file, " Yielding,");
        } else {
            fprintf(file, " Unloading,");
        }

        fprintf(file, " Plastic strains");
        for ( auto &val : plasticStrainVector ) {
            fprintf( file, " %.4e", val );
        }

	fprintf(file, " Damage %.4e", damage);



        if ( strainSpaceHardeningVarsVector.giveSize() ) {
            fprintf(file, " Strain space hardening vars");
            for ( auto &val : strainSpaceHardeningVarsVector ) {
                fprintf( file, " %.4e", val );
            }
        }

        fprintf(file, " ActiveConditionMap");
        for ( auto &val : activeConditionMap ) {
            fprintf( file, " %d", val );
        }
    }

    fprintf(file, "}\n");
}

int
MPlasticDamageMaterial2 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    MPlasticDamageMaterial2Status *status = static_cast< MPlasticDamageMaterial2Status * >( this->giveStatus(gp) );
    if ( type == IST_DamageScalar ) {
      answer.resize(1);
        answer.at(1) = status->giveTempDamage();
        return 1;
    } else {
        return MPlasticMaterial2 :: giveIPValue(answer, gp, type, tStep);
    }
}


  
void MPlasticDamageMaterial2Status :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
//
{
    MPlasticMaterial2Status :: initTempStatus();
    this->tempDamage = this->damage;
}


void
MPlasticDamageMaterial2Status :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    MPlasticMaterial2Status :: updateYourself(tStep);
    this->damage = this->tempDamage;
}




contextIOResultType
MPlasticDamageMaterial2Status :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = MPlasticMaterial2Status :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    
    return CIO_OK;
}



contextIOResultType
MPlasticDamageMaterial2Status :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = MPlasticMaterial2Status :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( !stream.read(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    return CIO_OK; // return succes
}
} // end namespace oofem
