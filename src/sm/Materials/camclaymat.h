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

#ifndef camclay_h
#define camclay_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "isolinearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "scalarfunction.h"

///@name Input fields for CamClay
//@{

#define _IFT_CamClayMat_Name "camclay"
#define _IFT_CamClayMat_M "m"
#define _IFT_CamClayMat_pc "pc"
#define _IFT_CamClayMat_yieldTol "yieldtol"
#define _IFT_CamClayMat_plStrainTol "plstraintol"
#define _IFT_CamClayMat_pcTol "pctol"
#define _IFT_CamClayMat_maxIter "maxiter"

#define _IFT_CamClayMat_e "voidratio"
#define _IFT_CamClayMat_consolidationIndex "lambda"
#define _IFT_CamClayMat_swellingIndex "kappa"

//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * This class implements an isotropic elastoplastic material
 * with Cam-Clay yield condition and associated flow rule
 * 
 *
 * It has two nostandard features
 * 1) The elastic stress-strain law is nonlinear with two possibilities of stress-strain law
 * 2) The yield function does not have standard isotropic-like hardening form
 */
class CamClayMat : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    IsotropicLinearElasticMaterial linearElasticMaterial;

    /// initial pre-consolidation pressure
    double pc0;

    /// parameter describing the ration of the minor to major axis of the yield function
    double M;

    /// elastic constants
    double G, K;

   
    // relative tolerance of the yield condition in the CPP algorithm
    double yieldTolerance;
    
    // relative tolerance of the plastic strain in the CPP algorithm
    double plasticStrainTolerance;

	// relative tolerance of the preconsolidation pressure in the CPP algorithm
	double pcTolerance;

    // maximum number of iterations of the CPP algorithm
    int maxIter;

	// void ratio
	double e;

	// virgin compression index lambda
	double consolidationIndex;

	// recompression/swelling index kappa
	double swellingIndex;

	// hardening parameter
	double theta;

public:
    CamClayMat(int n, Domain * d);
    virtual ~CamClayMat() {}

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep);
    IRResultType initializeFrom(InputRecord *ir) override;
    int hasNonLinearBehaviour() override { return 1; }
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) override { return true; }

    const char *giveInputRecordName() const override { return _IFT_CamClayMat_Name; }
    const char *giveClassName() const override { return "CamClayMat"; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep) override;

	void givePlasticGradientAtStress(FloatArray &answer, const FloatArray &stress, double pc);

	void giveYieldFunctionDoubleDerivative(FloatMatrix &answer);

	double giveYieldValueAtStress(FloatArray &stress, double &pc);

    void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep) override;


protected:
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};


class CamClayMatStatus : public StructuralMaterialStatus
{
protected:
    /// Plastic strain (initial).
    FloatArray plasticStrain;

    /// Plastic strain (final).
    FloatArray tempPlasticStrain;

    /// pre-consolidation pressure(initial).
    double pc;

    /// pre-consolidation pressure(initial).
    double tempPc;

    /// plastic multiplier (initial)
    double plasticMultiplier;
    
    /// plastic multiplier (final)
    double tempPlasticMultiplier;
    
public:
    CamClayMatStatus(int n, Domain * d, GaussPoint * g, double p0);
    virtual ~CamClayMatStatus();


    void letTempPlasticStrainBe(const FloatArray &values) { tempPlasticStrain = values; }
    void setTempPreconsolidationPressure(double pp) { tempPc = pp; }
    void setTempPlasticMultiplier(double pm) {tempPlasticMultiplier = pm;}


    double givePlasticMultiplier() {return plasticMultiplier;}
    double giveTempPlasticMultiplier() {return tempPlasticMultiplier;}
    

    
    const FloatArray &givePlasticStrain() { return plasticStrain; }
    const FloatArray &giveTempPlasticStrain() const { return tempPlasticStrain; }

    
    double givePreconsolidationPressure() { return pc; }
    double giveTempPreconsolidationPressure() { return tempPc; }
    
    
    



    void printOutputAt(FILE *file, TimeStep *tStep) override;

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;

    contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL) override;
    contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL) override;
    const char *giveClassName() const override { return "CamClayMatStatus"; }
};
} // end namespace oofem
#endif // camclay_h