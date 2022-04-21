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

#ifndef misesmat_h
#define misesmat_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "Materials/linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for MisesMat
//@{
#define _IFT_MisesMat_Name "misesmat"
#define _IFT_MisesMat_sig0 "sig0"
#define _IFT_MisesMat_sigInf "siginf"
#define _IFT_MisesMat_h "h"
#define _IFT_MisesMat_expD "expd"

#define _IFT_MisesMat_omega_crit "omega_crit"
#define _IFT_MisesMat_a "a"

// MJ 3.2.2022
#define _IFT_MisesMat_kappa1 "kappa1"
#define _IFT_MisesMat_dsig1 "dsig1"
#define _IFT_MisesMat_p1 "p1"
#define _IFT_MisesMat_kappa2 "kappa2"
#define _IFT_MisesMat_dsig2 "dsig2"
#define _IFT_MisesMat_p2 "p2"


//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * This class implements an isotropic elastoplastic material
 * with Mises yield condition, associated flow rule
 * and linear isotropic hardening.
 *
 * It differs from other similar materials (such as J2Mat)
 * by implementation - here we use the radial return, which
 * is the most efficient algorithm for this specific model.
 * Also, an extension to large strain will be available.
 * 
 * The model also exemplifies how to implement non-3d material modes, in this case 1D, 
 * by overloading the default implementations that iterates over the 3D method.
 */
class MisesMat : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    LinearElasticMaterial *linearElasticMaterial;

    /// Elastic shear modulus.
    double G;

    /// Elastic bulk modulus.
    double K;

    /// Hardening modulus.
    double sigInf;
    double sig0;

    
    double H;
    double expD;


    double omega_crit;
    double a;

    double yieldTol;

  // MJ 3.2.2022 - added parameters
  double kap1;
  double dsig1;
  double p1;
  double kap2;
  double dsig2;
  double p2;
    
public:
    MisesMat(int n, Domain * d);
    virtual ~MisesMat();

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);
    void performPlasticityReturn_PlaneStress(GaussPoint *gp, const FloatArray &totalStrain);
    double computeDamage(GaussPoint *gp, TimeStep *tStep);
    double computeDamageParam(double tempKappa);
    double computeDamageParamPrime(double tempKappa);
    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual const char *giveInputRecordName() const { return _IFT_MisesMat_Name; }
    virtual const char *giveClassName() const { return "MisesMat"; }

    /// Returns a reference to the basic elastic material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);
    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer,
	  MatResponseMode mmode, GaussPoint *gp,
                                          TimeStep *tStep);


    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep); 

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);

    double computeYieldStress(double kappa);
    double computeYieldStressPrime(double kappa);


    
protected:

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};

//=============================================================================


class MisesMatStatus : public StructuralMaterialStatus
{
protected:
    /// Plastic strain (initial).
    FloatArray plasticStrain;

    /// Plastic strain (final).
    FloatArray tempPlasticStrain;

    /// Deviatoric trial stress - needed for tangent stiffness.
    FloatArray trialStressD;

    /**************************************************/
    double trialStressV;
    /**************************************************/

    FloatArray effStress;
    FloatArray tempEffStress;

    /// Cumulative plastic strain (initial).
    double kappa;

    /// Cumulative plastic strain (final).
    double tempKappa;

    /************************/
    double tempDamage, damage;
    /******************************/

    double dGamma;


public:
    MisesMatStatus(int n, Domain * d, GaussPoint * g);
    virtual ~MisesMatStatus();

    const FloatArray &givePlasticStrain() { return plasticStrain; }

    const FloatArray &giveTrialStressDev() { return trialStressD; }

    /*******************************************/
    double giveTrialStressVol() { return trialStressV; }
    /*******************************************/
    double giveDamage() { return damage; }
    double giveTempDamage() { return tempDamage; }

    double giveCumulativePlasticStrain() { return kappa; }
    double giveTempCumulativePlasticStrain() { return tempKappa; }


    const FloatArray & giveTempEffectiveStress() { return tempEffStress; }
    const FloatArray & giveEffectiveStress() { return effStress; }

    void letTempPlasticStrainBe(const FloatArray &values) { tempPlasticStrain = values; }
    const FloatArray &getTempPlasticStrain() const { return tempPlasticStrain; }

    void letTrialStressDevBe(const FloatArray &values) { trialStressD = values; }

    void letEffectiveStressBe(const FloatArray &values) { effStress = values; }

    void letTempEffectiveStressBe(const FloatArray &values) { tempEffStress = values; }


    void setTrialStressVol(double value) { trialStressV = value; }

    void setDGamma(double dg) {dGamma = dg;}
    double giveDGamma() { return dGamma;}
    
    void setTempCumulativePlasticStrain(double value) { tempKappa = value; }
    /****************************************/
    void setTempDamage(double value) { tempDamage = value; }
    /************************************************/


    const FloatArray &givePlasDef() { return plasticStrain; }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "MisesMatStatus"; }
};
} // end namespace oofem
#endif // misesmat_h
