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

#ifndef idmmicromorphic_h
#define idmmicromorphic_h


#include "variationalbaseddamage.h"


#define _IFT_IsotropicDamageMaterialMicromorphic_Name "idmmicromorphic"
#define _IFT_IsotropicDamageMaterialMicromorphic_k "k"


namespace oofem {
/**
 * Micromorphic formulation of Variationally-based Gradient Isotropic Damage models,

 */
class IsotropicDamageMaterialMicromorphic : public VarBasedDamageMaterial
{
protected:
  double k1;
  double k2;

public:
    /// Constructor
    IsotropicDamageMaterialMicromorphic(int n, Domain * d);
    /// Destructor
    virtual ~IsotropicDamageMaterialMicromorphic();

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "IsotropicDamageMaterialMicromorphic"; }
    virtual const char *giveInputRecordName() const { return _IFT_IsotropicDamageMaterialMicromorphic_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    //    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual void giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);     
    virtual void giveGradientDamageStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_dd_l(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep);

    //    virtual void computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep);

    virtual void giveNonlocalInternalForces_N_factor(double &answer, double nlddv, GaussPoint *gp, TimeStep *tStep);
    virtual void giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlddv, GaussPoint *gp, TimeStep *tStep);
    
    
    //    virtual void computeDamage(double &answer, double micromorphicDamage, double damage, double storedEnergy, GaussPoint *gp);
    void computeDamageDrivingVariable(double &answer, double localDamageDrivingVariable_n, double micromorphicDamageDrivingVariable, double storedEnergy, GaussPoint *gp);

    //    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);


    
};


/**
 *  * Variatinally-based Gradient Isotropic Damage models,
 */
class IsotropicDamageMaterialMicromorphicStatus : public VarBasedDamageMaterialStatus
{
 protected:


 public:
    IsotropicDamageMaterialMicromorphicStatus(int n, Domain * d, GaussPoint * g);
    virtual ~IsotropicDamageMaterialMicromorphicStatus();

    virtual const char *giveClassName() const { return "IsotropicDamageMaterialMicromorphicStatus"; }
     
    

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *);

    
};
} // end namespace oofem
#endif // idmmicromorphic_h
