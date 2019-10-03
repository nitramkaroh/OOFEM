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

#include "plasticmaterials.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "isolinearelasticmaterial.h"
#include "mathfem.h"



namespace oofem {
  REGISTER_Material(VonMisesMaterial);
                 
VonMisesMaterial :: VonMisesMaterial(int n, Domain *d)  : SimplePlasticMaterial(n, d)
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);

}

  IRResultType
VonMisesMaterial :: initializeFrom(InputRecord *ir)
{

    // Required by IR_GIVE_FIELD macro
    IRResultType result;
    // call the corresponding service of simple plastic material
    result = SimplePlasticMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    result = linearElasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    
    IR_GIVE_FIELD(ir, sig0, _IFT_VonMisesMaterial_sig0); // uniaxial yield stress
    
    return IRRT_OK;
    
}
  
double
VonMisesMaterial :: computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                             const FloatArray &strainSpaceHardeningVars)
{
    double f = sqrt( 3. * this->computeJ2StressInvariantAt(stressVector) );
    return f - this->sig0  - this->giveIsotropicHardeningStress(strainSpaceHardeningVars.at(strainSpaceHardeningVars.giveSize()));
}

  

void
VonMisesMaterial :: computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector, const FloatArray &strainSpaceHardeningVars)
{
    //stress gradient of yield function in full stress - strain space

    answer.resize(6);
    answer.zero();
    this->compute_dSqrtJ2StressInvariant_dStressAt(answer, stressVector);
    answer.times(sqrt(3.));
    
}


void
VonMisesMaterial :: computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int isurf, GaussPoint *gp, const FloatArray &fullStressVector,
                                        const FloatArray &strainSpaceHardeningVars)
{
  IntArray mask;
  FloatMatrix fullGradientMatrix;
  //full form of gradient matirx
  this->compute_d2SqrtJ2StressInvariant_dStress2At(fullGradientMatrix, fullStressVector);
  fullGradientMatrix.times(sqrt(3.));
  // transforming to reduced form
  // put it maybe in the simple plastic material instead ??
  StructuralMaterial :: giveVoigtSymVectorMask( mask, gp->giveMaterialMode());
  gradientMatrix.beSubMatrixOf(fullGradientMatrix, mask, mask);

}



void
VonMisesMaterial :: computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf,
                                                const IntArray &activeConditionMap,
                                                const FloatArray &fullStressVector,
                                                const FloatArray &strainSpaceHardeningVars,
                                                const FloatArray &gamma)
{

 int size = this->giveSizeOfReducedHardeningVarsVector(gp);
    answer.resize(size, 1);

    if ( this->kinematicHardeningFlag ) {
    }

    if ( isotropicHardeningFlag ) {
      answer.at(size,1) = 1;
    }
  
}

void
VonMisesMaterial :: computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp,
                                             const FloatArray &stress, const FloatArray &dlambda,
                                             const FloatArray &dplasticStrain, const IntArray &activeConditionMaps)
{

  int size = this->giveSizeOfReducedHardeningVarsVector(gp);
  answer.resize(size);
  FloatArray strainSpaceHardeningVariables(size);
  strainSpaceHardeningVariables.zero();
  
    if ( this->kinematicHardeningFlag ) {
    }

    if ( isotropicHardeningFlag ) {
      answer.at(size) = dlambda.at(1);
    }
}





  


  

} // end namespace oofem
