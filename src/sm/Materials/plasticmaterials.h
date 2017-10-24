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

#ifndef plasticmaterials_h
#define plasticmaterials_h


#include "Materials/simpleplasticmaterial.h"
///@name Input fields for VonMisesMaterial
//@{
#define _IFT_VonMisesMaterial_Name "vonmisesmat"
#define _IFT_VonMisesMaterial_sig0 "sig0"
//@}

namespace oofem {
class Domain;
/**
 * This class implements 
 */
 class VonMisesMaterial : public SimplePlasticMaterial
{
protected:
  double sig0;
public:
    VonMisesMaterial(int n, Domain * d);
    virtual ~VonMisesMaterial(){;}

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveInputRecordName() const { return _IFT_VonMisesMaterial_Name; }
    virtual const char *giveClassName() const { return "VonMisesMaterial"; }
    
protected:  
    virtual double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                                       const FloatArray &strainSpaceHardeningVars);

    virtual void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector, const FloatArray &strainSpaceHardeningVars);

    virtual void  computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables);


    virtual void computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &gamma);

    virtual void computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp,
                                             const FloatArray &stress, const FloatArray &dlambda,
						     const FloatArray &dplasticStrain, const IntArray &activeConditionMaps);
    
};
 
} // end namespace oofem
#endif // plasticmaterials_h
