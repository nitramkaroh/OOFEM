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

#ifndef lourencomasonrymat_h
#define lourencomasonrymat_h

#include "Materials/ConcreteMaterials/mplasticmaterial2.h"

///@name Input fields for LouencoMasonryMat
//@{
#define _IFT_LourencoMasonryMat_Name "lourencomasonrymat"
#define _IFT_LourencoMasonryMat_fcx "fcx" ///< Compressive strength in direction x
#define _IFT_LourencoMasonryMat_fcy "fcy" ///< Compressive strength in direction y
#define _IFT_LourencoMasonryMat_ftx "ftx" ///< Tensile strength in direction x
#define _IFT_LourencoMasonryMat_fty "fty" ///< Tensile strength in direction y
#define _IFT_LourencoMasonryMat_tauU "tauu" ///< Shear strength
#define _IFT_LourencoMasonryMat_Gcx "gcx" ///< Fracture energy in compression in direction x
#define _IFT_LourencoMasonryMat_Gcy "gcy" ///< Fracture energy in compression in direction y
#define _IFT_LourencoMasonryMat_Gtx "gtx" ///< Fracture energy in tension in direction x
#define _IFT_LourencoMasonryMat_Gty "gty" ///< Fracture energy in tension in direction y
#define _IFT_LourencoMasonryMat_beta "beta" ///< Coupling parameter between compressive stresses in two directions
#define _IFT_LourencoMasonryMat_kappap "kappap" ///< Value of equivalent plastic strain corresponding to peak of compressive yield stress
//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * This class implements an orthotropic elasto-plastic material in the framework of multi-surface plasticity
 * with a Rankine-type yield condition in tension, non-associated flow rule, kinematic orthotropic hardening, and
 * with a Hill-type yield condition in compression, associated flow rule, mixed orthotropic hardening.
 */
class LourencoMasonryMat : public MPlasticMaterial2
{
protected:
    // Reference to the basic elastic material.
    //LinearElasticMaterial *linearElasticMaterial;

    /// Compressive strength in x direction
    double fcx;

    /// Compressive strength in y direction
    double fcy;

    /// Tensile strength in x direction
    double ftx;

    /// Tensile strength in y direction
    double fty;

    /// Shear strength
    double tauU;

    ///  Fracture energy in compression in x direction
    double Gcx;

    ///  Fracture energy in compression in y direction
    double Gcy;

    ///  Fracture energy in tension in x direction
    double Gtx;

    ///  Fracture energy in tension in y direction
    double Gty;

    /// Coupling parameter between compressive stresses in two directions
    double beta;

    /// Value of equivalent plastic strain corresponding to peak of compressive yield stress
    double kappap;
 



public:
    LourencoMasonryMat(int n, Domain * d);
    virtual ~LourencoMasonryMat();

    int hasMaterialModeCapability(MaterialMode mode) override;

    IRResultType initializeFrom(InputRecord *ir) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) override { return false; }

    const char *giveClassName() const override { return "LourencoMasonryMat"; }
    const char *giveInputRecordName() const override { return _IFT_LourencoMasonryMat_Name; }

    /// Returns a reference to the basic elastic material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    int giveSizeOfFullHardeningVarsVector() override { return 2; }
    int giveSizeOfReducedHardeningVarsVector(GaussPoint *) const override { return 2; } //cummulative strain = one per each surface

protected:
    int giveMaxNumberOfActiveYieldConds(GaussPoint *gp) override { return 2; } //normally one less than number of all conditions //but in this case no, both yields function could be active in the same time

    double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector, const FloatArray &strainSpaceHardeningVariables) override;

   
    void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector, const FloatArray &stressSpaceHardeningVars) override;

    /// Computes second derivative of yield/loading function with respect to stress
    void computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables) override;

    void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) override;

    /// Functions related to hardening
    int hasHardening() override { return 1; }

    /// Compute dot(kappa_1), dot(kappa_2) etc.
    void computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp, const FloatArray &stress, const FloatArray &dlambda, const FloatArray &dplasticStrain, const IntArray &activeConditionMap) override;

    /// Computes the derivative of yield/loading function with respect to kappa_1, kappa_2 etc.
    void computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables) override;

    /// computes mixed derivative of load function with respect to stress and hardening variables
    void computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix, int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables) override;

    /// computes dk(i)/dsig(j) gradient matrix
    void computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &dlambda) override;

    /// computes dKappa_i/dLambda_j
    void computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &dlambda) override;



    double computeTensileYieldStress( double E, double f, double G, double kappa, double le );
    double computeTensileYieldStressKGradient( double E, double f, double G, double kappa, double le );
    double computeCompressiveYieldStress( double E, double f, double G, double kappa, double le );
    double computeCompressiveYieldStressKGradient( double E, double f, double G, double kappa, double le );

    void computePiVector(FloatArray &pi);
    void computeTensileBackStressVector(FloatArray& answer, double kappa, double le, GaussPoint *gp  );
    void computeCompressiveBackStressVector(FloatArray& answer, double kappa, double le, GaussPoint *gp  );
    void computeTensileBackStressKGradientVector(FloatArray& answer,  double kappa, double le, GaussPoint *gp  );
    void computeCompressiveBackStressKGradientVector(FloatArray& answer, double kappa, double le, GaussPoint *gp  );
    void computePcMatrix(FloatMatrix&, double kappa, double le, GaussPoint *gp );
    void computePtMatrix(FloatMatrix&);
    void computePgMatrix(FloatMatrix&);


      
      



    
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // lourencomasonrymat_h
