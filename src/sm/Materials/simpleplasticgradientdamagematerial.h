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

#ifndef simpleplasticgradientmaterial_h
#define simpleplasticgradientmaterial_h


#include "Materials/simpleplasticmaterial.h"
#include "Materials/graddamagematerialextensioninterface.h"

///@name Input fields for SimplePlasticGradientDamageMaterial
//@{
#define _IFT_SimplePlasticGradientDamageMaterial_m "m"
//@}





namespace oofem {
class Domain;
/**
 * This class implements a general class of simple plastic materials with damage
 * it assumes decomposition of energy in form psi(epsilon) + psi(kappa).
 * only isotropic hardening and isotropic damage are supported now 
 */
 class SimplePlasticGradientDamageMaterial : public SimplePlasticMaterial, public GradientDamageMaterialExtensionInterface
 {
protected:
   double mParam;
   double l;

public:
    SimplePlasticGradientDamageMaterial(int n, Domain * d);
    virtual ~SimplePlasticGradientDamageMaterial();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "SimplePlasticGradientDamageMaterial"; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
protected:
    virtual double computeDamage(GaussPoint *gp, const FloatArray &strainSpaceHardeningVariables, TimeStep *tStep);
    virtual double compute_dDamage_dKappa(GaussPoint *gp, const FloatArray &strainSpaceHardeningVariables, TimeStep *tStep);
    

  /// Left upper block
    virtual void giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    /// Left lower block
    virtual void giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    /// Right upper block
    virtual void giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    /// Right lower block
    virtual void giveGradientDamageStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    /// Stress-based averaging
    virtual void giveGradientDamageStiffnessMatrix_dd_l(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    
    /// gradient - based giveRealStressVector
    virtual void giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalDamageDrivningVariable, TimeStep *tStep);
    virtual void giveFirstPKStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalDamageDrivningVariable, TimeStep *tStep) { OOFEM_ERROR("not implemented") }


    virtual void giveInternalLengthMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual Interface *giveInterface(InterfaceType t); 

    virtual void giveNonlocalInternalForces_N_factor(double &answer,double nlddv, GaussPoint *gp, TimeStep *tStep);
  virtual void giveNonlocalInternalForces_B_factor(FloatArray &answer,const FloatArray &nlddv, GaussPoint *gp, TimeStep *tStep);
    void computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep);

    
};


    ////dodelat status
class SimplePlasticGradientDamageMaterialStatus : public MPlasticMaterial2Status, public GradientDamageMaterialStatusExtensionInterface
{
protected:

public:
  SimplePlasticGradientDamageMaterialStatus(int n, Domain * d, GaussPoint * g,int statusSize);
  virtual ~SimplePlasticGradientDamageMaterialStatus(){;}
    // definition
    virtual const char *giveClassName() const { return "SimplePlasticGradientDamageMaterialStatus"; }


    virtual double giveNonlocalCumulatedStrain() { return nonlocalDamageDrivingVariable; }
    virtual void setNonlocalCumulatedStrain(double nonlocalCumulatedStrain) { this->nonlocalDamageDrivingVariable = nonlocalCumulatedStrain; }
};
 
} // end namespace oofem
#endif // simpleplasticgradientmaterial_h
