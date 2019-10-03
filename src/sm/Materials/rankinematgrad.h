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

#ifndef RankineMatGrad_h

#include "rankinemat.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "graddamagematerialextensioninterface.h"
#include "cltypes.h"

///@name Input fields for RankineMatGrad
//@{
#define _IFT_RankineMatGrad_Name "rankinematgrad"
#define _IFT_RankineMatGrad_m "m"
#define _IFT_RankineMatGrad_negligibleDamage "negligible_damage"
#define _IFT_RankineMatGrad_formulationType "formtype"
//@}

namespace oofem {
/**
 * Gradient rankine material status.
 */
class RankineMatGradStatus : public RankineMatStatus
{
protected:

    double kappa_hat;
    double tempNonlocalCumulativePlasticStrain;
    double nonlocalCumulativePlasticStrain;
    
    FloatArray nonlocalCumulativePlasticStrainGrad;
    FloatArray tempNonlocalCumulativePlasticStrainGrad;


public:

    RankineMatGradStatus(int n, Domain *d, GaussPoint *g);
    virtual ~RankineMatGradStatus() { }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "RankineMatGradStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    void setKappa_hat(double kap) { kappa_hat = kap; }
    double giveKappa_hat() { return kappa_hat; }
    
    virtual double giveTempNonlocalCumulativePlasticStrain() { return tempNonlocalCumulativePlasticStrain; }
    virtual void setTempNonlocalCumulativePlasticStrain(double nonlocalCumulativeStrain) { this->tempNonlocalCumulativePlasticStrain = nonlocalCumulativeStrain; }


    const FloatArray &giveTempNonlocalCumulativePlasticStrainGrad() const { return tempNonlocalCumulativePlasticStrainGrad; }
    virtual void letTempNonlocalCumulativePlasticStrainGradBe(const FloatArray &ncpsg) { this->tempNonlocalCumulativePlasticStrainGrad = ncpsg; }

    
    
};


/**
 * Gradient Rankine material.
 */
class RankineMatGrad : public RankineMat, GradientDamageMaterialExtensionInterface
{
protected:
    double mParam;
    double negligible_damage;

    enum GradientDamageFormulationType {
        GDFT_Standard = 0,
        GDFT_Eikonal = 2
    };

    GradientDamageFormulationType gradientDamageFormulationType;


public:
    RankineMatGrad(int n, Domain *d);
    virtual ~RankineMatGrad() { }

    virtual const char *giveClassName() const { return "RankineMatGrad"; }
    virtual const char *giveInputRecordName() const { return _IFT_RankineMatGrad_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == GradientDamageMaterialExtensionInterfaceType ) {
            return static_cast< GradientDamageMaterialExtensionInterface * >( this );
        } else {
            return NULL;
        }
    }

    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual void giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_du_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_dd_NN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveGradientDamageStiffnessMatrix_dd_BN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual void computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep);

    virtual void giveNonlocalInternalForces_N_factor(double &answer, double nlddv, GaussPoint *gp, TimeStep *tStep);
    virtual void giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlddv, GaussPoint *gp, TimeStep *tStep);


    virtual void giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep);

    virtual void givePlaneStressStiffMtrx(FloatMatrix & answer, MatResponseMode, GaussPoint * gp,  TimeStep * tStep);
    void givePlaneStressGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

protected:
    double computeInternalLength(GaussPoint *gp);
    int giveDimension(GaussPoint *gp);

    double computeEikonalInternalLength_a(GaussPoint *gp);
    double computeEikonalInternalLength_b(GaussPoint *gp);
    double computeEikonalInternalLength_aPrime(GaussPoint *gp);
    double computeEikonalInternalLength_bPrime(GaussPoint *gp);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new RankineMatGradStatus(1, this->giveDomain(), gp); }
};
} // end namespace oofem
#define RankineMatGrad_h
#endif
