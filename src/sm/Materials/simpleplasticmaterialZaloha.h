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

#ifndef simpleplasticmaterial_h
#define simpleplasticmaterial_h


#include "Materials/ConcreteMaterials/mplasticmaterial2.h"
#include "Materials/isohardeningmaterialextensioninterface.h"


namespace oofem {
class Domain;
/**
 * This class implements a general class of simple plastic materials
 * it assumes decomposition of energy in form psi(epsilon) + psi(kappa).
 * only isotropic hardening is supported now.
 */
 class SimplePlasticMaterial : public MPlasticMaterial2
{
protected:
    int kinematicHardeningFlag, isotropicHardeningFlag;
    double kinematicModuli, isotropicModuli;
    double k;

public:
    SimplePlasticMaterial(int n, Domain * d);
    virtual ~SimplePlasticMaterial();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveInputRecordName() const { return "huhu"; }
    virtual const char *giveClassName() const { return "huhu"; }

    virtual int giveSizeOfFullHardeningVarsVector();
    virtual int giveSizeOfReducedHardeningVarsVector(GaussPoint *gp) const;
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    virtual int giveMaxNumberOfActiveYieldConds(GaussPoint *gp) { return 2; }

    virtual double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                                       const FloatArray &strainSpaceHardeningVars);

    virtual void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp,
                                             const FloatArray &stressVector, const FloatArray &strainSpaceHardeningVars);

    virtual void computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp,
                                                     const FloatArray &stress, const FloatArray &dlambda,
                                                     const FloatArray &dplasticStrain, const IntArray &activeConditionMap);
    virtual void computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector,
                                        const FloatArray &strainSpaceHardeningVariables);

    virtual void computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap,
                                                          const FloatArray &fullStressVector,
                                                          const FloatArray &strainSpaceHardeningVars,
                                                          const FloatArray &gamma);
    virtual void computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf,
                                                        const IntArray &activeConditionMap,
                                                        const FloatArray &fullStressVector,
                                                        const FloatArray &strainSpaceHardeningVars,
                                                        const FloatArray &gamma);
    virtual int hasHardening();
    /* virtual void  computeReducedGradientMatrix (FloatMatrix& answer, int isurf,
     *                                          GaussPoint *gp,
     *                                          const FloatArray& stressVector,
     *                                          const FloatArray& stressSpaceHardeningVars) = 0;*/
    virtual void  computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                                 const FloatArray &strainSpaceHardeningVariables);
    virtual void  computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                                 const FloatArray &strainSpaceHardeningVariables);

    // auxiliary function
    static double computeJ2InvariantAt(const FloatArray &stressVector);
    double giveIsotropicHardeningVar(GaussPoint *gp, const FloatArray &strainSpaceHardeningVars);
    void giveStressBackVector(FloatArray &answer, GaussPoint *gp, const FloatArray &strainSpaceHardeningVars);
};
} // end namespace oofem
#endif // mplasticmaterial2_h
