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

namespace oofem {
REGISTER_Material(MooneyRivlinMaterial);

MooneyRivlinMaterial :: MooneyRivlinMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{ }


void
MooneyRivlinMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{
    double J, lnJ, I1, I2, barI1, barI2;
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
    this->compute_cross_product(iFtxiFt, invFt, invFt);
    
    d2I1dF2.times(C1);
    d2I2dF2.times(C2);

    iFtxiFt.times(K);
    dInvF_dF.times(K * lnJ);

    FloatMatrix answer1;
    answer1 = d2I1dF2;
    answer1.add(d2I2dF2);
    answer1.add(iFtxiFt);
    answer1.add(dInvF_dF);
    

    

    double I1;
    FloatMatrix C, Cpow2, FC;
    // compute right Cauchy-Green tensor
    C.beTProductOf(F, F);
    //compute first invariant of deviatoric part of C;
    I1 = ( C.at(1, 1) + C.at(2, 2) + C.at(3, 3) );
    answer.resize(9,9);
    answer.zero();

    answer.at(1, 1) = ( 2. * C1 * ( 5. * I1 * invF.at(1, 1) * invF.at(1, 1) - 12. * F.at(1, 1) * invF.at(1, 1) + 9. ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(1, 2) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(2, 2) + 6. * F.at(2, 2) * invF.at(1, 1) - 2. * I1 * invF.at(1, 1) * invF.at(2, 2) - 3. * I1 * invF.at(1, 2) * invF.at(2, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(1, 3) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(1, 1) - 2. * I1 * invF.at(1, 1) * invF.at(3, 3) - 3. * I1 * invF.at(1, 3) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(1, 4) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(3, 2) + 6. * F.at(2, 3) * invF.at(1, 1) - 2. * I1 * invF.at(1, 1) * invF.at(3, 2) - 3. * I1 * invF.at(1, 2) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(1, 5) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(3, 1) + 6. * F.at(1, 3) * invF.at(1, 1) - 5. * I1 * invF.at(1, 1) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(1, 6) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(2, 1) + 6. * F.at(1, 2) * invF.at(1, 1) - 5. * I1 * invF.at(1, 1) * invF.at(2, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(1, 7) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(1, 1) - 2. * I1 * invF.at(1, 1) * invF.at(2, 3) - 3. * I1 * invF.at(2, 1) * invF.at(1, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(1, 8) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(1, 1) - 5. * I1 * invF.at(1, 1) * invF.at(1, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(1, 9) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(1, 2) + 6. * F.at(2, 1) * invF.at(1, 1) - 5. * I1 * invF.at(1, 1) * invF.at(1, 2) ) ) / ( 9. * pow(J, 2. / 3.) );




    answer.at(2, 1) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(2, 2) + 6. * F.at(2, 2) * invF.at(1, 1) - 2. * I1 * invF.at(1, 1) * invF.at(2, 2) - 3. * I1 * invF.at(1, 2) * invF.at(2, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(2, 2) = ( 2. * C1 * ( 5. * I1 * invF.at(2, 2) * invF.at(2, 2) - 12. * F.at(2, 2) * invF.at(2, 2) + 9 ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(2, 3) = -( 2. * C1 * ( 6. * F.at(2, 2) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(2, 2) - 2. * I1 * invF.at(2, 2) * invF.at(3, 3) - 3. * I1 * invF.at(2, 3) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(2, 4) = -( 2. * C1 * ( 6. * F.at(2, 2) * invF.at(3, 2) + 6. * F.at(2, 3) * invF.at(2, 2) - 5. * I1 * invF.at(2, 2) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(2, 5) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(2, 2) + 6. * F.at(2, 2) * invF.at(3, 1) - 3. * I1 * invF.at(2, 1) * invF.at(3, 2) - 2. * I1 * invF.at(2, 2) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(2, 6) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(2, 2) + 6. * F.at(2, 2) * invF.at(2, 1) - 5. * I1 * invF.at(2, 1) * invF.at(2, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(2, 7) = -( 2. * C1 * ( 6. * F.at(2, 2) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(2, 2) - 5. * I1 * invF.at(2, 2) * invF.at(2, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(2, 8) = -( 2. * C1 * ( 6. * F.at(2, 2) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(2, 2) - 3. * I1 * invF.at(1, 2) * invF.at(2, 3) - 2. * I1 * invF.at(1, 3) * invF.at(2, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(2, 9) = -( 2. * C1 * ( 6. * F.at(2, 1) * invF.at(2, 2) + 6. * F.at(2, 2) * invF.at(1, 2) - 5. * I1 * invF.at(1, 2) * invF.at(2, 2) ) ) / ( 9. * pow(J, 2. / 3.) );


    answer.at(3, 1) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(1, 1) - 2. * I1 * invF.at(1, 1) * invF.at(3, 3) - 3. * I1 * invF.at(1, 3) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(3, 2) = -( 2. * C1 * ( 6. * F.at(2, 2) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(2, 2) - 2. * I1 * invF.at(2, 2) * invF.at(3, 3) - 3. * I1 * invF.at(2, 3) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(3, 3) = ( 2. * C1 * ( 5. * I1 * invF.at(3, 3) * invF.at(3, 3) - 12. * F.at(3, 3) * invF.at(3, 3) + 9 ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(3, 4) = -( 2. * C1 * ( 6. * F.at(2, 3) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(3, 2) - 5. * I1 * invF.at(3, 2) * invF.at(3, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(3, 5) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(3, 1) - 5. * I1 * invF.at(3, 1) * invF.at(3, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(3, 6) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(2, 1) - 2. * I1 * invF.at(2, 1) * invF.at(3, 3) - 3. * I1 * invF.at(3, 1) * invF.at(2, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(3, 7) = -( 2. * C1 * ( 6. * F.at(3, 2) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(2, 3) - 5. * I1 * invF.at(2, 3) * invF.at(3, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(3, 8) = -( 2. * C1 * ( 6. * F.at(3, 1) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(1, 3) - 5. * I1 * invF.at(1, 3) * invF.at(3, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(3, 9) = -( 2. * C1 * ( 6. * F.at(2, 1) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(1, 2) - 2. * I1 * invF.at(1, 2) * invF.at(3, 3) - 3. * I1 * invF.at(1, 3) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );




    answer.at(4, 1) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(3, 2) + 6. * F.at(2, 3) * invF.at(1, 1) - 2. * I1 * invF.at(1, 1) * invF.at(3, 2) - 3. * I1 * invF.at(1, 2) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(4, 2) = -( 2. * C1 * ( 6. * F.at(2, 2) * invF.at(3, 2) + 6. * F.at(2, 3) * invF.at(2, 2) - 5. * I1 * invF.at(2, 2) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(4, 3) = -( 2. * C1 * ( 6. * F.at(2, 3) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(3, 2) - 5. * I1 * invF.at(3, 2) * invF.at(3, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(4, 4) = ( 2. * C1 * ( 5. * I1 * invF.at(3, 2) * invF.at(3, 2) - 12. * F.at(2, 3) * invF.at(3, 2) + 9 ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(4, 5) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(3, 2) + 6. * F.at(2, 3) * invF.at(3, 1) - 5. * I1 * invF.at(3, 1) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(4, 6) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(3, 2) + 6. * F.at(2, 3) * invF.at(2, 1) - 2. * I1 * invF.at(2, 1) * invF.at(3, 2) - 3. * I1 * invF.at(2, 2) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(4, 7) = -( 2. * C1 * ( 6. * F.at(2, 3) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(3, 2) - 3. * I1 * invF.at(2, 2) * invF.at(3, 3) - 2. * I1 * invF.at(2, 3) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(4, 8) = -( 2. * C1 * ( 6. * F.at(2, 3) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(3, 2) - 3. * I1 * invF.at(1, 2) * invF.at(3, 3) - 2. * I1 * invF.at(1, 3) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(4, 9) = -( 2. * C1 * ( 6. * F.at(2, 1) * invF.at(3, 2) + 6. * F.at(2, 3) * invF.at(1, 2) - 5. * I1 * invF.at(1, 2) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );



    answer.at(5, 1) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(3, 1) + 6. * F.at(1, 3) * invF.at(1, 1) - 5. * I1 * invF.at(1, 1) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(5, 2) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(2, 2) + 6. * F.at(2, 2) * invF.at(3, 1) - 3. * I1 * invF.at(2, 1) * invF.at(3, 2) - 2. * I1 * invF.at(2, 2) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(5, 3) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(3, 1) - 5. * I1 * invF.at(3, 1) * invF.at(3, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(5, 4) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(3, 2) + 6. * F.at(2, 3) * invF.at(3, 1) - 5. * I1 * invF.at(3, 1) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(5, 5) = ( 2. * C1 * ( 5. * I1 * invF.at(3, 1) * invF.at(3, 1) - 12. * F.at(1, 3) * invF.at(3, 1) + 9 ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(5, 6) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(3, 1) + 6. * F.at(1, 3) * invF.at(2, 1) - 5. * I1 * invF.at(2, 1) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(5, 7) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(3, 1) - 3. * I1 * invF.at(2, 1) * invF.at(3, 3) - 2. * I1 * invF.at(3, 1) * invF.at(2, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(5, 8) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(3, 1) - 3. * I1 * invF.at(1, 1) * invF.at(3, 3) - 2. * I1 * invF.at(1, 3) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(5, 9) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(1, 2) + 6. * F.at(2, 1) * invF.at(3, 1) - 3. * I1 * invF.at(1, 1) * invF.at(3, 2) - 2. * I1 * invF.at(1, 2) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );



    answer.at(6, 1) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(2, 1) + 6. * F.at(1, 2) * invF.at(1, 1) - 5. * I1 * invF.at(1, 1) * invF.at(2, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(6, 2) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(2, 2) + 6. * F.at(2, 2) * invF.at(2, 1) - 5. * I1 * invF.at(2, 1) * invF.at(2, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(6, 3) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(2, 1) - 2. * I1 * invF.at(2, 1) * invF.at(3, 3) - 3. * I1 * invF.at(3, 1) * invF.at(2, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(6, 4) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(3, 2) + 6. * F.at(2, 3) * invF.at(2, 1) - 2. * I1 * invF.at(2, 1) * invF.at(3, 2) - 3. * I1 * invF.at(2, 2) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(6, 5) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(3, 1) + 6. * F.at(1, 3) * invF.at(2, 1) - 5. * I1 * invF.at(2, 1) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(6, 6) = ( 2. * C1 * ( 5. * I1 * invF.at(2, 1) * invF.at(2, 1) - 12. * F.at(1, 2) * invF.at(2, 1) + 9 ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(6, 7) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(2, 1) - 5. * I1 * invF.at(2, 1) * invF.at(2, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(6, 8) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(2, 1) - 3. * I1 * invF.at(1, 1) * invF.at(2, 3) - 2. * I1 * invF.at(2, 1) * invF.at(1, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(6, 9) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(1, 2) + 6. * F.at(2, 1) * invF.at(2, 1) - 3. * I1 * invF.at(1, 1) * invF.at(2, 2) - 2. * I1 * invF.at(1, 2) * invF.at(2, 1) ) ) / ( 9. * pow(J, 2. / 3.) );




    answer.at(7, 1) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(1, 1) - 2. * I1 * invF.at(1, 1) * invF.at(2, 3) - 3. * I1 * invF.at(2, 1) * invF.at(1, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(7, 2) = -( 2. * C1 * ( 6. * F.at(2, 2) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(2, 2) - 5. * I1 * invF.at(2, 2) * invF.at(2, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(7, 3) = -( 2. * C1 * ( 6. * F.at(3, 2) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(2, 3) - 5. * I1 * invF.at(2, 3) * invF.at(3, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(7, 4) = -( 2. * C1 * ( 6. * F.at(2, 3) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(3, 2) - 3. * I1 * invF.at(2, 2) * invF.at(3, 3) - 2. * I1 * invF.at(2, 3) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(7, 5) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(3, 1) - 3. * I1 * invF.at(2, 1) * invF.at(3, 3) - 2. * I1 * invF.at(3, 1) * invF.at(2, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(7, 6) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(2, 1) - 5. * I1 * invF.at(2, 1) * invF.at(2, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(7, 7) = ( 2. * C1 * ( 5. * I1 * invF.at(2, 3) * invF.at(2, 3) - 12. * F.at(3, 2) * invF.at(2, 3) + 9 ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(7, 8) = -( 2. * C1 * ( 6. * F.at(3, 1) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(1, 3) - 5. * I1 * invF.at(1, 3) * invF.at(2, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(7, 9) = -( 2. * C1 * ( 6. * F.at(2, 1) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(1, 2) - 2. * I1 * invF.at(1, 2) * invF.at(2, 3) - 3. * I1 * invF.at(1, 3) * invF.at(2, 2) ) ) / ( 9. * pow(J, 2. / 3.) );





    answer.at(8, 1) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(1, 1) - 5. * I1 * invF.at(1, 1) * invF.at(1, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(8, 2) = -( 2. * C1 * ( 6. * F.at(2, 2) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(2, 2) - 3. * I1 * invF.at(1, 2) * invF.at(2, 3) - 2. * I1 * invF.at(1, 3) * invF.at(2, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(8, 3) = -( 2. * C1 * ( 6. * F.at(3, 1) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(1, 3) - 5. * I1 * invF.at(1, 3) * invF.at(3, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(8, 4) = -( 2. * C1 * ( 6. * F.at(2, 3) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(3, 2) - 3. * I1 * invF.at(1, 2) * invF.at(3, 3) - 2. * I1 * invF.at(1, 3) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(8, 5) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(3, 1) - 3. * I1 * invF.at(1, 1) * invF.at(3, 3) - 2. * I1 * invF.at(1, 3) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(8, 6) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(2, 1) - 3. * I1 * invF.at(1, 1) * invF.at(2, 3) - 2. * I1 * invF.at(2, 1) * invF.at(1, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(8, 7) = -( 2. * C1 * ( 6. * F.at(3, 1) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(1, 3) - 5. * I1 * invF.at(1, 3) * invF.at(2, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(8, 8) = ( 2. * C1 * ( 5. * I1 * invF.at(1, 3) * invF.at(1, 3) - 12. * F.at(3, 1) * invF.at(1, 3) + 9 ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(8, 9) = -( 2. * C1 * ( 6. * F.at(2, 1) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(1, 2) - 5. * I1 * invF.at(1, 2) * invF.at(1, 3) ) ) / ( 9. * pow(J, 2. / 3.) );





    answer.at(9, 1) = -( 2. * C1 * ( 6. * F.at(1, 1) * invF.at(1, 2) + 6. * F.at(2, 1) * invF.at(1, 1) - 5. * I1 * invF.at(1, 1) * invF.at(1, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(9, 2) = -( 2. * C1 * ( 6. * F.at(2, 1) * invF.at(2, 2) + 6. * F.at(2, 2) * invF.at(1, 2) - 5. * I1 * invF.at(1, 2) * invF.at(2, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(9, 3) = -( 2. * C1 * ( 6. * F.at(2, 1) * invF.at(3, 3) + 6. * F.at(3, 3) * invF.at(1, 2) - 2. * I1 * invF.at(1, 2) * invF.at(3, 3) - 3. * I1 * invF.at(1, 3) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(9, 4) = -( 2. * C1 * ( 6. * F.at(2, 1) * invF.at(3, 2) + 6. * F.at(2, 3) * invF.at(1, 2) - 5. * I1 * invF.at(1, 2) * invF.at(3, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(9, 5) = -( 2. * C1 * ( 6. * F.at(1, 3) * invF.at(1, 2) + 6. * F.at(2, 1) * invF.at(3, 1) - 3. * I1 * invF.at(1, 1) * invF.at(3, 2) - 2. * I1 * invF.at(1, 2) * invF.at(3, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(9, 6) = -( 2. * C1 * ( 6. * F.at(1, 2) * invF.at(1, 2) + 6. * F.at(2, 1) * invF.at(2, 1) - 3. * I1 * invF.at(1, 1) * invF.at(2, 2) - 2. * I1 * invF.at(1, 2) * invF.at(2, 1) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(9, 7) = -( 2. * C1 * ( 6. * F.at(2, 1) * invF.at(2, 3) + 6. * F.at(3, 2) * invF.at(1, 2) - 2. * I1 * invF.at(1, 2) * invF.at(2, 3) - 3. * I1 * invF.at(1, 3) * invF.at(2, 2) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(9, 8) = -( 2. * C1 * ( 6. * F.at(2, 1) * invF.at(1, 3) + 6. * F.at(3, 1) * invF.at(1, 2) - 5. * I1 * invF.at(1, 2) * invF.at(1, 3) ) ) / ( 9. * pow(J, 2. / 3.) );
    answer.at(9, 9) = ( 2. * C1 * ( 5. * I1 * invF.at(1, 2) * invF.at(1, 2) - 12. * F.at(2, 1) * invF.at(1, 2) + 9 ) ) / ( 9. * pow(J, 2. / 3.) );

    

    
    answer.at(1, 1) = answer.at(1, 1) - K *invF.at(1, 1) * invF.at(1, 1) * ( lnJ - 1. );
    answer.at(1, 2) = answer.at(1, 2) + K * ( invF.at(1, 1) * invF.at(2, 2) - invF.at(1, 2) * invF.at(2, 1) * lnJ );
    answer.at(1, 3) = answer.at(1, 3) + K * ( invF.at(1, 1) * invF.at(3, 3) - invF.at(1, 3) * invF.at(3, 1) * lnJ );
    answer.at(1, 4) = answer.at(1, 4) + K * ( invF.at(1, 1) * invF.at(3, 2) - invF.at(1, 2) * invF.at(3, 1) * lnJ );
    answer.at(1, 5) = answer.at(1, 5) - K *invF.at(1, 1) * invF.at(3, 1) * ( lnJ - 1. );
    answer.at(1, 6) = answer.at(1, 6) - K *invF.at(1, 1) * invF.at(2, 1) * ( lnJ - 1. );
    answer.at(1, 7) = answer.at(1, 7) + K * ( invF.at(1, 1) * invF.at(2, 3) - invF.at(2, 1) * invF.at(1, 3) * lnJ );
    answer.at(1, 8) = answer.at(1, 8) - K *invF.at(1, 1) * invF.at(1, 3) * ( lnJ - 1. );
    answer.at(1, 9) = answer.at(1, 9) - K *invF.at(1, 1) * invF.at(1, 2) * ( lnJ - 1. );




    answer.at(2, 1) = answer.at(2, 1) + K * ( invF.at(1, 1) * invF.at(2, 2) - invF.at(1, 2) * invF.at(2, 1) * lnJ );
    answer.at(2, 2) = answer.at(2, 2) - K *invF.at(2, 2) * invF.at(2, 2) * ( lnJ - 1. );
    answer.at(2, 3) = answer.at(2, 3) + K * ( invF.at(2, 2) * invF.at(3, 3) - invF.at(2, 3) * invF.at(3, 2) * lnJ );
    answer.at(2, 4) = answer.at(2, 4) - K *invF.at(2, 2) * invF.at(3, 2) * ( lnJ - 1. );
    answer.at(2, 5) = answer.at(2, 5) + K * ( invF.at(2, 2) * invF.at(3, 1) - invF.at(2, 1) * invF.at(3, 2) * lnJ );
    answer.at(2, 6) = answer.at(2, 6) - K *invF.at(2, 1) * invF.at(2, 2) * ( lnJ - 1. );
    answer.at(2, 7) = answer.at(2, 7) - K *invF.at(2, 2) * invF.at(2, 3) * ( lnJ - 1. );
    answer.at(2, 8) = answer.at(2, 8) + K * ( invF.at(1, 3) * invF.at(2, 2) - invF.at(1, 2) * invF.at(2, 3) * lnJ );
    answer.at(2, 9) = answer.at(2, 9) - K *invF.at(1, 2) * invF.at(2, 2) * ( lnJ - 1. );




    answer.at(3, 1) = answer.at(3, 1) + K * ( invF.at(1, 1) * invF.at(3, 3) - invF.at(1, 3) * invF.at(3, 1) * lnJ );
    answer.at(3, 2) = answer.at(3, 2) + K * ( invF.at(2, 2) * invF.at(3, 3) - invF.at(2, 3) * invF.at(3, 2) * lnJ );
    answer.at(3, 3) = answer.at(3, 3) - K *invF.at(3, 3) * invF.at(3, 3) * ( lnJ - 1. );
    answer.at(3, 4) = answer.at(3, 4) - K *invF.at(3, 2) * invF.at(3, 3) * ( lnJ - 1. );
    answer.at(3, 5) = answer.at(3, 5) - K *invF.at(3, 1) * invF.at(3, 3) * ( lnJ - 1. );
    answer.at(3, 6) = answer.at(3, 6) + K * ( invF.at(2, 1) * invF.at(3, 3) - invF.at(3, 1) * invF.at(2, 3) * lnJ );
    answer.at(3, 7) = answer.at(3, 7) - K *invF.at(2, 3) * invF.at(3, 3) * ( lnJ - 1. );
    answer.at(3, 8) = answer.at(3, 8) - K *invF.at(1, 3) * invF.at(3, 3) * ( lnJ - 1. );
    answer.at(3, 9) = answer.at(3, 9) + K * ( invF.at(1, 2) * invF.at(3, 3) - invF.at(1, 3) * invF.at(3, 2) * lnJ );




    answer.at(4, 1) = answer.at(4, 1) + K * ( invF.at(1, 1) * invF.at(3, 2) - invF.at(1, 2) * invF.at(3, 1) * lnJ );
    answer.at(4, 2) = answer.at(4, 2) - K *invF.at(2, 2) * invF.at(3, 2) * ( lnJ - 1. );
    answer.at(4, 3) = answer.at(4, 3) - K *invF.at(3, 2) * invF.at(3, 3) * ( lnJ - 1. );
    answer.at(4, 4) = answer.at(4, 4) - K *invF.at(3, 2) * invF.at(3, 2) * ( lnJ - 1. );
    answer.at(4, 5) = answer.at(4, 5) - K *invF.at(3, 1) * invF.at(3, 2) * ( lnJ - 1. );
    answer.at(4, 6) = answer.at(4, 6) + K * ( invF.at(2, 1) * invF.at(3, 2) - invF.at(2, 2) * invF.at(3, 1) * lnJ );
    answer.at(4, 7) = answer.at(4, 7) + K * ( invF.at(2, 3) * invF.at(3, 2) - invF.at(2, 2) * invF.at(3, 3) * lnJ );
    answer.at(4, 8) = answer.at(4, 8) + K * ( invF.at(1, 3) * invF.at(3, 2) - invF.at(1, 2) * invF.at(3, 3) * lnJ );
    answer.at(4, 9) = answer.at(4, 9) - K *invF.at(1, 2) * invF.at(3, 2) * ( lnJ - 1. );



    answer.at(5, 1) = answer.at(5, 1) - K *invF.at(1, 1) * invF.at(3, 1) * ( lnJ - 1. );
    answer.at(5, 2) = answer.at(5, 2) + K * ( invF.at(2, 2) * invF.at(3, 1) - invF.at(2, 1) * invF.at(3, 2) * lnJ );
    answer.at(5, 3) = answer.at(5, 3) - K *invF.at(3, 1) * invF.at(3, 3) * ( lnJ - 1. );
    answer.at(5, 4) = answer.at(5, 4) - K *invF.at(3, 1) * invF.at(3, 2) * ( lnJ - 1. );
    answer.at(5, 5) = answer.at(5, 5) - K *invF.at(3, 1) * invF.at(3, 1) * ( lnJ - 1. );
    answer.at(5, 6) = answer.at(5, 6) - K *invF.at(2, 1) * invF.at(3, 1) * ( lnJ - 1. );
    answer.at(5, 7) = answer.at(5, 7) + K * ( invF.at(3, 1) * invF.at(2, 3) - invF.at(2, 1) * invF.at(3, 3) * lnJ );
    answer.at(5, 8) = answer.at(5, 8) + K * ( invF.at(1, 3) * invF.at(3, 1) - invF.at(1, 1) * invF.at(3, 3) * lnJ );
    answer.at(5, 9) = answer.at(5, 9) + K * ( invF.at(1, 2) * invF.at(3, 1) - invF.at(1, 1) * invF.at(3, 2) * lnJ );




    answer.at(6, 1) = answer.at(6, 1) - K *invF.at(1, 1) * invF.at(2, 1) * ( lnJ - 1. );
    answer.at(6, 2) = answer.at(6, 2) - K *invF.at(2, 1) * invF.at(2, 2) * ( lnJ - 1. );
    answer.at(6, 3) = answer.at(6, 3) + K * ( invF.at(2, 1) * invF.at(3, 3) - invF.at(3, 1) * invF.at(2, 3) * lnJ );
    answer.at(6, 4) = answer.at(6, 4) + K * ( invF.at(2, 1) * invF.at(3, 2) - invF.at(2, 2) * invF.at(3, 1) * lnJ );
    answer.at(6, 5) = answer.at(6, 5) - K *invF.at(2, 1) * invF.at(3, 1) * ( lnJ - 1. );
    answer.at(6, 6) = answer.at(6, 6) - K *invF.at(2, 1) * invF.at(2, 1) * ( lnJ - 1. );
    answer.at(6, 7) = answer.at(6, 7) - K *invF.at(2, 1) * invF.at(2, 3) * ( lnJ - 1. );
    answer.at(6, 8) = answer.at(6, 8) + K * ( invF.at(2, 1) * invF.at(1, 3) - invF.at(1, 1) * invF.at(2, 3) * lnJ );
    answer.at(6, 9) = answer.at(6, 9) + K * ( invF.at(1, 2) * invF.at(2, 1) - invF.at(1, 1) * invF.at(2, 2) * lnJ );





    answer.at(7, 1) = answer.at(7, 1) + K * ( invF.at(1, 1) * invF.at(2, 3) - invF.at(2, 1) * invF.at(1, 3) * lnJ );
    answer.at(7, 2) = answer.at(7, 2) - K *invF.at(2, 2) * invF.at(2, 3) * ( lnJ - 1. );
    answer.at(7, 3) = answer.at(7, 3) - K *invF.at(2, 3) * invF.at(3, 3) * ( lnJ - 1. );
    answer.at(7, 4) = answer.at(7, 4) + K * ( invF.at(2, 3) * invF.at(3, 2) - invF.at(2, 2) * invF.at(3, 3) * lnJ );
    answer.at(7, 5) = answer.at(7, 5) + K * ( invF.at(3, 1) * invF.at(2, 3) - invF.at(2, 1) * invF.at(3, 3) * lnJ );
    answer.at(7, 6) = answer.at(7, 6) - K *invF.at(2, 1) * invF.at(2, 3) * ( lnJ - 1. );
    answer.at(7, 7) = answer.at(7, 7) - K *invF.at(2, 3) * invF.at(2, 3) * ( lnJ - 1. );
    answer.at(7, 8) = answer.at(7, 8) - K *invF.at(1, 3) * invF.at(2, 3) * ( lnJ - 1. );
    answer.at(7, 9) = answer.at(7, 9) + K * ( invF.at(1, 2) * invF.at(2, 3) - invF.at(1, 3) * invF.at(2, 2) * lnJ );




    answer.at(8, 1) = answer.at(8, 1) - K *invF.at(1, 1) * invF.at(1, 3) * ( lnJ - 1. );
    answer.at(8, 2) = answer.at(8, 2) + K * ( invF.at(1, 3) * invF.at(2, 2) - invF.at(1, 2) * invF.at(2, 3) * lnJ );
    answer.at(8, 3) = answer.at(8, 3) - K *invF.at(1, 3) * invF.at(3, 3) * ( lnJ - 1. );
    answer.at(8, 4) = answer.at(8, 4) + K * ( invF.at(1, 3) * invF.at(3, 2) - invF.at(1, 2) * invF.at(3, 3) * lnJ );
    answer.at(8, 5) = answer.at(8, 5) + K * ( invF.at(1, 3) * invF.at(3, 1) - invF.at(1, 1) * invF.at(3, 3) * lnJ );
    answer.at(8, 6) = answer.at(8, 6) + K * ( invF.at(2, 1) * invF.at(1, 3) - invF.at(1, 1) * invF.at(2, 3) * lnJ );
    answer.at(8, 7) = answer.at(8, 7) - K *invF.at(1, 3) * invF.at(2, 3) * ( lnJ - 1. );
    answer.at(8, 8) = answer.at(8, 8) - K *invF.at(1, 3) * invF.at(1, 3) * ( lnJ - 1. );
    answer.at(8, 9) = answer.at(8, 9) - K *invF.at(1, 2) * invF.at(1, 3) * ( lnJ - 1. );





    answer.at(9, 1) = answer.at(9, 1) - K *invF.at(1, 1) * invF.at(1, 2) * ( lnJ - 1. );
    answer.at(9, 2) = answer.at(9, 2) - K *invF.at(1, 2) * invF.at(2, 2) * ( lnJ - 1. );
    answer.at(9, 3) = answer.at(9, 3) + K * ( invF.at(1, 2) * invF.at(3, 3) - invF.at(1, 3) * invF.at(3, 2) * lnJ );
    answer.at(9, 4) = answer.at(9, 4) - K *invF.at(1, 2) * invF.at(3, 2) * ( lnJ - 1. );
    answer.at(9, 5) = answer.at(9, 5) + K * ( invF.at(1, 2) * invF.at(3, 1) - invF.at(1, 1) * invF.at(3, 2) * lnJ );
    answer.at(9, 6) = answer.at(9, 6) + K * ( invF.at(1, 2) * invF.at(2, 1) - invF.at(1, 1) * invF.at(2, 2) * lnJ );
    answer.at(9, 7) = answer.at(9, 7) + K * ( invF.at(1, 2) * invF.at(2, 3) - invF.at(1, 3) * invF.at(2, 2) * lnJ );
    answer.at(9, 8) = answer.at(9, 8) - K *invF.at(1, 2) * invF.at(1, 3) * ( lnJ - 1. );
    answer.at(9, 9) = answer.at(9, 9) - K *invF.at(1, 2) * invF.at(1, 2) * ( lnJ - 1. );
    



    FloatMatrix test(answer);
    test.subtract(answer1);
    answer = answer1;
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
