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

#include "mooneyrivlin.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {
REGISTER_Material(MooneyRivlinMaterial);

MooneyRivlinMaterial :: MooneyRivlinMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{ }


void
MooneyRivlinMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{
    double J, lnJ;
    FloatMatrix F, C, Cpow2, invFt, FC, invF;

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    // compute jacobian and its logarith
    J = F.giveDeterminant();
    lnJ = log(J);

    FloatArray dI1_Cdev_dF, dI2_Cdev_dF, dJ_dF;
    this->compute_dI1_Cdev_dF(dI1_Cdev_dF, F);
    //dI1_Cdev_dF.times(C1);
    this->compute_dI2_Cdev_dF(dI2_Cdev_dF, F);
    //dI2_Cdev_dF.times(C2);
    this->compute_dJ_dF(dJ_dF, F);
    //    dJ_dF.times(K * lnJ / J);
   
    /*    invF.beInverseOf(F);
    invFt.beTranspositionOf(invF);

    // compute right Cauchy-Green tensor
    C.beTProductOf(F, F);
    // compute C*C
    Cpow2.beProductOf(C, C);
    // compute F*C
    FC.beProductOf(F, C);
    //compute first invariant of deviatoric part of C;
    I1 = ( C.at(1, 1) + C.at(2, 2) + C.at(3, 3) );
    I2 = 1. / 2. * ( I1 * I1 - Cpow2.at(1, 1) - Cpow2.at(2, 2) - Cpow2.at(3, 3) );
    barI1 = I1 * pow(J, -2. / 3.);
    barI2 = I2 * pow(J, -4. / 3.);


    FloatMatrix P;
    //first part of stress tensor : C1 * \frac{\partial \bar{I}_1}{\partial F_ij}
    P.add(2 * C1 / pow(J, 2 / 3.), F);
    P.add(-2. / 3. * C1 * barI1, invFt);
    // second part of stress tensor : C2 * \frac{\partial \bar{I}_2}{\partial F_ij}
    P.add( 2. * C2 * barI1 / pow(J, 2. / 3.), F );
    P.add(-4. / 3. * C2 * barI2, invFt);
    P.add(-2. * C2 / pow(J, 4. / 3.), FC );
    // third part of stress tensor : K * \frac{\partial ln J }{F_ij}
    P.add(K * lnJ, invFt);
    FloatMatrix test(invFt);
    test.times(K* lnJ);

    FloatArray answer1;
    answer1.beVectorForm(P);
    */
    
    answer.zero();
    answer.add(C1, dI1_Cdev_dF);
    answer.add(C2, dI2_Cdev_dF);
    answer.add(K * lnJ / J, dJ_dF);
    
    /*answer = dI1_Cdev_dF;
    answer.add(dI2_Cdev_dF);
    answer.add(dJ_dF);
    */
    
    //answer.beVectorForm(P);

    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);
}



void
MooneyRivlinMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double J, lnJ;
    FloatArray vF;
    FloatMatrix F, invF, invFt, d2I1dF2, d2I2dF2, dinvF_dF, iFtxiFt, dInvF_dF;

    //deformation gradient from the status
    vF = status->giveTempFVector();
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    invF.beInverseOf(F);
    invFt.beTranspositionOf(invF);
    J = F.giveDeterminant();
    lnJ = log(J);
    this->compute_d2I1_Cdev_dF2_and_d2I2_Cdev_dF2(d2I1dF2, d2I2dF2, F);
    FloatMatrix stiff1, stiff2;
    this->computeNumerical_d2I1_Cdev_dF2(stiff1, stiff2, F);


    this->compute_dInvFt_dF(dInvF_dF, invF);
    this->compute_dyadic_product(iFtxiFt, invFt, invFt);
    
    d2I1dF2.times(C1);
    d2I2dF2.times(C2);

    iFtxiFt.times(K);
    dInvF_dF.times(K * lnJ);

    answer = d2I1dF2;
    answer.add(d2I2dF2);
    answer.add(iFtxiFt);
    answer.add(dInvF_dF);


    FloatMatrix KK;
    KK.add(iFtxiFt);
    KK.add(dInvF_dF);
    
    

    

}



  


MaterialStatus *
MooneyRivlinMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}


IRResultType
MooneyRivlinMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro


    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, K, _IFT_MooneyRivlinMaterial_k);
    IR_GIVE_FIELD(ir, C1, _IFT_MooneyRivlinMaterial_c1);
    IR_GIVE_FIELD(ir, C2, _IFT_MooneyRivlinMaterial_c2);

    return IRRT_OK;
}

} // end namespace oofem
