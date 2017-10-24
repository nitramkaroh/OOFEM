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

#ifndef variationalbaseddamage_h
#define variationalbaseddamage_h

#include "Materials/ConcreteMaterials/idm1.h"
#include "../sm/Materials/graddamagematerialextensioninterface.h"

#define _IFT_VarBasedDamageMaterial_Name "varbaseddamagematerial"
#define _IFT_VarBasedDamageMaterial_initDamage "initdamage"


namespace oofem {
/**
 * Gradient-enhanced Isotropic Damage model for concrete in tension,
 */
class VarBasedDamageMaterial : public IsotropicDamageMaterial1, GradientDamageMaterialExtensionInterface
{
protected:
  double initialDamage;

public:
    /// Constructor
    VarBasedDamageMaterial(int n, Domain * d);
    /// Destructor
    virtual ~VarBasedDamageMaterial();

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "VarBasedDamageMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_VarBasedDamageMaterial_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == GradientDamageMaterialExtensionInterfaceType ) {
            return static_cast< GradientDamageMaterialExtensionInterface * >(this);
        } else {
            return NULL;
        }
    }
    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual void giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_dd_l(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep);

    void giveStiffnessMatrix(FloatMatrix &answer,  MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    
    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode, GaussPoint *gp,  TimeStep *tStep);
    void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep);

    virtual void giveNonlocalInternalForces_N_factor(double &answer,GaussPoint *gp, TimeStep *tStep);
    virtual void giveNonlocalInternalForces_B_factor(double &answer,GaussPoint *gp, TimeStep *tStep);
    virtual void computeInternalForcesRegularizationTerm(double &answer,GaussPoint *gp, TimeStep *tStep);
    virtual void computeStiffnessRegularizationTerm(double &answer,GaussPoint *gp, TimeStep *tStep);
    virtual void computeDamage(double &answer, double damageDrivingVariable, GaussPoint *gp);
};


/**
 * Material status for gradient-enhanced isotropic damage model for concrete in tension.
 */
class VarBasedDamageMaterialStatus : public IsotropicDamageMaterial1Status, public GradientDamageMaterialStatusExtensionInterface
{
 protected:
    FloatArray effectiveStressVector;
    FloatArray tempEffectiveStressVector;
    double strainEnergy, tempStrainEnergy;

 public:
    VarBasedDamageMaterialStatus(int n, Domain * d, GaussPoint * g, double initialDamage);
    virtual ~VarBasedDamageMaterialStatus();

    virtual const char *giveClassName() const { return "VarBasedDamageMaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

      virtual void letTempEffectiveStressVectorBe(const FloatArray &effectiveStress){this->tempEffectiveStressVector = effectiveStress;}
    const FloatArray &giveTempEffectiveStressVector() const { return tempEffectiveStressVector; }
	
    double giveStrainEnergy(){return strainEnergy;}
    void setTempStrainEnergy(double sE){strainEnergy = sE;}

};
} // end namespace oofem
#endif // variationalbaseddamage_h
