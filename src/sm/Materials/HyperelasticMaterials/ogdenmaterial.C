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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#include "ogdenmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {
REGISTER_Material(OgdenMaterial);

OgdenMaterial :: OgdenMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{ }


void
OgdenMaterial :: giveSecondPKStressVector_3d(FloatMatrix &answer, const FloatMatrix &C, const FloatArray &eVals, const FloatMatrix &eVecs,TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    double J, lnJ;
    FloatMatrix invC, I1aInvC;
    //
    J = sqrt(C.giveDeterminant());
    //compute inverse of the right Cauchy-Green tensor(C)
    invC.beInverseOf(C);
    for (int i = 1; i <= N; i++) {
      double l1a, l2a, l3a, Ja;
      FloatMatrix Si;
      //
      l1a = pow(eVals.at(1),alpha.at(i)/2.);
      l2a = pow(eVals.at(2),alpha.at(i)/2.);
      l3a = pow(eVals.at(3),alpha.at(i)/2.);
      Ja = pow(J, -alpha.at(i) / 3.);
      double I1a = (l1a + l2a + l3a) / 3.;
      I1aInvC = invC;
      I1aInvC.times(I1a);
      this->computeMatrixPower(Si, eVals, eVecs, ( alpha.at(i)-2. ) / 2. );
      
      Si.subtract(I1aInvC);
      Si.times(mu.at(i) * Ja );
      answer.add(Si);
    }
    
    // compute jacobian and its logarithm
    lnJ = log(J);
    answer.add(K * lnJ, invC);

}
  

void
OgdenMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{
  FloatArray eVals, vS;
    FloatMatrix P, F, S, C, eVecs;
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    //
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    
    this->giveSecondPKStressVector_3d(S, C, eVals, eVecs, tStep);
    P.beProductOf(F,S);
    vS.beSymVectorForm(S);
    answer.beVectorForm(P);
    // update gp
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);
    status->letTempStressVectorBe(vS);
}



void
OgdenMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double J, lnJ;
    FloatArray vF, eVals;
    FloatMatrix P, F, S, C, invC, iCiC, dSdE;
    FloatMatrix eVecs, invF, invFt, d2I1dF2, d2I2dF2, dinvF_dF, iFtxiFt, dInvF_dF;
    
    //deformation gradient from the status
    vF = status->giveTempFVector();
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    J = F.giveDeterminant();
    lnJ = log(J);
    //
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    //compute inverse of the right Cauchy-Green tensor(C)
    invC.beInverseOf(C);
    //compute symetric dyadic product 
    this->compute_sym_dyadic_product_reduced(iCiC, invC, invC);
    //
    for (int i = 1; i <= N; i++) {
      double l1a, l2a, l3a, Ja;
      //
      Ja = pow(J, -alpha.at(i) / 3.);
      l1a = pow(eVals.at(1),alpha.at(i)/2.);
      l2a = pow(eVals.at(2),alpha.at(i)/2.);
      l3a = pow(eVals.at(3),alpha.at(i)/2.);
      //
      double I1a = (l1a + l2a + l3a) / 3.;
      //
      FloatMatrix I1aInvC;
      I1aInvC = invC;
      I1aInvC.times(I1a);
      //
      FloatMatrix Sdev, powC, SiC;
      this->computeMatrixPower(powC, eVals, eVecs, ( alpha.at(i)-2. ) / 2. );
      Sdev = powC;
      Sdev.subtract(I1aInvC);
      this->compute_dyadic_product_reduced(SiC, Sdev, invC);

      FloatMatrix test;
      this->compute_dyadic_product(test, Sdev, invC);
      //
      FloatMatrix iCpowC, dCm_dC;
      this->compute_dyadic_product_reduced(iCpowC, invC, powC);
      this->compute_dCm_dC(dCm_dC, (alpha.at(i) - 2.) / 2., eVals, eVecs);

      FloatMatrix dC1_dC;
      this->compute_dCm_dC(dC1_dC, -1, eVals, eVecs);

      
      dSdE.add( 2. * mu.at(i) * Ja, dCm_dC);
      dSdE.add( -mu.at(i) * alpha.at(i) * Ja / 3., iCpowC);
      dSdE.add( 2. * mu.at(i)  * Ja * I1a, iCiC);
      dSdE.add( -alpha.at(i) * mu.at(i) * Ja / 3., SiC);
    }

    
    this->give_dPdF_from(dSdE, answer, gp, _3dMat);
    
    
    invF.beInverseOf(F);
    invFt.beTranspositionOf(invF);
    this->compute_dInvFt_dF(dInvF_dF, invF);
    answer.add(K * lnJ, dInvF_dF);
  
    this->compute_dyadic_product(iFtxiFt, invFt, invFt);
    answer.add(K, iFtxiFt);   

}



  


MaterialStatus *
OgdenMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}


IRResultType
OgdenMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro


    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, K, _IFT_OgdenMaterial_k);
    IR_GIVE_FIELD(ir, alpha, _IFT_OgdenMaterial_alpha);
    IR_GIVE_FIELD(ir, mu, _IFT_OgdenMaterial_mu);

    N = alpha.giveSize();
    int M = mu.giveSize();
    if(N != M) {
      OOFEM_ERROR("Inconsistent size of alpha and mu");
    }
    
    return IRRT_OK;
}

} // end namespace oofem
